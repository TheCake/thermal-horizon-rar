"""FITPOOL -- the Tier-1 CPU worker pool for the 10-series fit engine.

Nothing rides on this pool until `py calcs/fitpool.py gate` has
certified, ON THIS MACHINE, that pooled results are bit-identical to
the in-process serial reference at two different worker counts.

Determinism contract:
  1. Workers are separate PROCESSES launched via subprocess, NOT the
     multiprocessing module -- Windows spawn re-imports the parent
     __main__, and our stage scripts are unguarded module-level flows
     (a spawn child would re-execute the whole stage).
  2. A task is a PURE function of its content (kind, sub, member,
     th0). All randomness is drawn in the PARENT and baked into the
     sub before submission; workers never touch an rng.
  3. Results are returned KEYED and never read by order; assignment
     is round-robin (worker i takes tasks[i::n]), so the worker count
     cannot change any result -- and the gate PROVES it by running
     one battery serial + two pool widths and requiring bit-equality.
  4. Worker BLAS/OMP threading is pinned to 1 via the child env (and
     in-process before numpy loads), matching the serial reference.

Worker engine = exec-inheritance of calcs/stage10v_strat.py truncated
at the baselines marker (the flowboot route): SPARC/CF4 parse, FAMS,
build_sub, fit_deep, contest_fit -- BIT-VERBATIM, read-only at load.
A worker that dies OR silently STOPs (the engine's census STOP exits
0) is detected by its missing output pickle; the parent raises with
that worker's log tail.

API (parent side):
    sys.path.insert(0, 'calcs'); from fitpool import run_tasks
    res = run_tasks(tasks, subs, nworkers=None, tag='x')
      tasks: [{'key': str, 'kind': 'fit_deep'|'contest_fit',
               'sub': <name in subs>, 'member': m, 'th0': [4 floats]}]
      subs:  {name: build_sub(...) dict}  (deduplicated shipping)
      res:   {key: {'fun': float, 'x': [floats], 'gap': float}}
"""
import os, sys

THREAD_ENV = {k: '1' for k in (
    'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
    'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')}
os.environ.update(THREAD_ENV)          # before any numpy import

import math, pickle, shutil, subprocess, tempfile, time  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MARKER = "# ---------------- baselines (both modes)"


def load_engine():
    """Exec-inherit stage10v_strat.py to the baselines marker (the
    flowboot route); returns the truncated namespace dict."""
    src = open(os.path.join(HERE, 'stage10v_strat.py'),
               encoding='utf-8').read()
    ns = {'__name__': 'stage10v_trunc',
          '__file__': os.path.join(HERE, 'stage10v_strat.py')}
    argv = sys.argv
    sys.argv = ['stage10v_strat.py', 'gates']
    exec(compile(src[:src.index(MARKER)], 'stage10v_trunc', 'exec'), ns)
    sys.argv = argv
    return ns


def exec_task(ns, task, subs):
    """The ONE task-execution path, shared by worker and serial
    reference -- same source, same bits."""
    sub = subs[task['sub']]
    m = task['member']
    if task['kind'] == 'fit_deep':
        bf = ns['fit_deep'](sub, ns['FAMS'][m], th0=list(task['th0']))
        gap = 0.0
    elif task['kind'] == 'contest_fit':
        bf, gap = ns['contest_fit'](sub, m, list(task['th0']))
    else:
        raise ValueError(f"unknown kind {task['kind']!r}")
    return {'fun': float(bf.fun), 'x': [float(v) for v in bf.x],
            'gap': float(gap)}


# ---------------- worker mode --------------------------------------
def worker_main(subs_path, tasks_path, out_path):
    t0 = time.time()
    subs = pickle.load(open(subs_path, 'rb'))
    tasks = pickle.load(open(tasks_path, 'rb'))
    ns = load_engine()
    print(f"[worker] engine loaded in {time.time()-t0:.1f}s; "
          f"{len(tasks)} tasks", flush=True)
    out = {}
    for i, t in enumerate(tasks):
        out[t['key']] = exec_task(ns, t, subs)
        with open(out_path + '.prog', 'w') as f:
            f.write(f"{i+1}\n")
        print(f"[worker] {t['key']} done  [{time.time()-t0:.0f}s]",
              flush=True)
    with open(out_path + '.tmp', 'wb') as f:
        pickle.dump(out, f)
    os.replace(out_path + '.tmp', out_path)
    print(f"[worker] all {len(tasks)} tasks done "
          f"[{time.time()-t0:.0f}s]", flush=True)


# ---------------- parent side --------------------------------------
def default_workers():
    n = os.cpu_count() or 4
    return max(2, min(n - 2, 16))


def run_tasks(tasks, subs, nworkers=None, tag='pool', log=print):
    keys = [t['key'] for t in tasks]
    assert len(keys) == len(set(keys)), "duplicate task keys"
    for t in tasks:
        assert t['sub'] in subs, f"task {t['key']}: unknown sub"
    nw = min(nworkers or default_workers(), len(tasks))
    tmp = tempfile.mkdtemp(prefix=f'fitpool_{tag}_')
    t0 = time.time()
    try:
        subs_path = os.path.join(tmp, 'subs.pkl')
        with open(subs_path, 'wb') as f:
            pickle.dump(subs, f)
        env = dict(os.environ, **THREAD_ENV)
        procs, outs, logs = [], [], []
        for i in range(nw):
            part = tasks[i::nw]
            tp = os.path.join(tmp, f'tasks_{i}.pkl')
            op = os.path.join(tmp, f'out_{i}.pkl')
            lp = os.path.join(tmp, f'log_{i}.txt')
            with open(tp, 'wb') as f:
                pickle.dump(part, f)
            lf = open(lp, 'w', encoding='utf-8')
            procs.append(subprocess.Popen(
                [sys.executable, os.path.abspath(__file__),
                 'worker', subs_path, tp, op],
                cwd=ROOT, env=env, stdout=lf, stderr=subprocess.STDOUT))
            outs.append(op); logs.append((lp, lf))
        log(f"  [{tag}] {len(tasks)} tasks -> {nw} workers")
        last = 0.0
        while any(p.poll() is None for p in procs):
            time.sleep(5)
            if time.time() - last >= 30:
                done = 0
                for op in outs:
                    try:
                        done += int(open(op + '.prog').read())
                    except Exception:
                        pass
                log(f"  [{tag}] {done}/{len(tasks)} done "
                    f"[{(time.time()-t0)/60:.1f} min]")
                last = time.time()
        for _, lf in logs:
            lf.close()
        bad = []
        for i, (p, op) in enumerate(zip(procs, outs)):
            if p.returncode != 0 or not os.path.exists(op):
                tail = open(logs[i][0], encoding='utf-8',
                            errors='replace').read()[-2000:]
                bad.append(f"worker {i} rc={p.returncode} "
                           f"out={'ok' if os.path.exists(op) else 'MISSING'}"
                           f"\n--- log tail ---\n{tail}")
        if bad:
            raise RuntimeError(
                f"[{tag}] {len(bad)} worker(s) failed (tmp kept: {tmp})\n"
                + "\n".join(bad))
        res = {}
        for op in outs:
            res.update(pickle.load(open(op, 'rb')))
        assert set(res) == set(keys), "result keys != task keys"
        log(f"  [{tag}] complete in {(time.time()-t0)/60:.1f} min")
    except Exception:
        raise
    else:
        shutil.rmtree(tmp, ignore_errors=True)
    return res


# ---------------- gate mode ----------------------------------------
def gate_main():
    L = []

    def P(s=""):
        print(s, flush=True)
        L.append(s)

    P("FITPOOL DETERMINISM GATE (pool == serial, two widths)")
    P(f"cpu_count = {os.cpu_count()}; default workers = "
      f"{default_workers()}")
    P("")
    t0 = time.time()
    ns = load_engine()
    P(f"parent engine loaded [{time.time()-t0:.0f}s]")
    import numpy as np
    W78, LEGA, FLOW = ns['W78'], ns['LEGA'], ns['FLOW']
    ARC, MEMBERS = ns['ARC'], ns['MEMBERS']
    build_sub = ns['build_sub']

    subs = {'anch': build_sub(W78, LEGA), 'flow': build_sub(W78, FLOW)}
    # parent-drawn jittered replicate (flowboot construction, seed 424)
    rng = np.random.default_rng(424)
    pick = rng.choice(FLOW, size=len(FLOW), replace=True)
    docc = []
    for g in pick:
        s_g = 1.0 + rng.normal(0, W78['rel'][int(g)])
        docc.append(-math.log10(min(max(s_g, 0.5), 1.5)))
    subs['jit'] = build_sub(W78, [int(g) for g in pick],
                            dlg_per_occ=docc)

    tasks = []
    for leg in ('anch', 'flow'):
        for m in MEMBERS:
            tasks.append(dict(key=f'deep:{leg}:{m}', kind='fit_deep',
                              sub=leg, member=m,
                              th0=[math.log10(ARC[m]), 1.0, 0.08, 0.0]))
    for m in MEMBERS:
        tasks.append(dict(key=f'contest:flow:{m}', kind='contest_fit',
                          sub='flow', member=m,
                          th0=[math.log10(ARC[m]), 1.0, 0.08, 0.0]))
    for m in MEMBERS:
        tasks.append(dict(key=f'deep:jit:{m}', kind='fit_deep',
                          sub='jit', member=m,
                          th0=[math.log10(ARC[m]), 1.0, 0.08, 0.0]))
    P(f"battery: {len(tasks)} tasks "
      f"(8 deep baselines + 4 flow contests + 4 jittered deeps)")
    P("")

    t1 = time.time()
    ser = {t['key']: exec_task(ns, t, subs) for t in tasks}
    t_ser = time.time() - t1
    P(f"serial reference: {t_ser/60:.1f} min "
      f"({t_ser/len(tasks):.1f} s/task)")

    verdict = True
    times = {}
    for nw in (4, default_workers()):
        t1 = time.time()
        res = run_tasks(tasks, subs, nworkers=nw, tag=f'gate{nw}',
                        log=P)
        times[nw] = time.time() - t1
        nfun = nx = ngap = 0
        for k in ser:
            if res[k]['fun'] == ser[k]['fun']: nfun += 1
            if res[k]['x'] == ser[k]['x']: nx += 1
            if res[k]['gap'] == ser[k]['gap']: ngap += 1
        ok = (nfun == nx == ngap == len(tasks))
        verdict &= ok
        P(f"pool@{nw}: fun exact {nfun}/{len(tasks)}, x exact "
          f"{nx}/{len(tasks)}, gap exact {ngap}/{len(tasks)}, wall "
          f"{times[nw]/60:.1f} min (speedup x{t_ser/times[nw]:.1f}) "
          f"-> {'BIT-IDENTICAL' if ok else 'MISMATCH'}")
        if not ok:
            for k in sorted(ser):
                if res[k] != ser[k]:
                    P(f"  DIFFER {k}: serial {ser[k]['fun']!r} pool "
                      f"{res[k]['fun']!r}")
    P("")
    msg = ("bit-identical to serial at both widths -- the pool may "
           "carry science" if verdict else
           "NOT identical -- the pool ships nothing")
    P(f"GATE {'PASS' if verdict else 'FAIL'}: pooled fits are {msg}")
    P(f"wall-clock total: {(time.time()-t0)/60:.1f} min")
    with open(os.path.join(ROOT, 'data', 'fitpool_gate.txt'), 'w',
              encoding='utf-8') as f:
        f.write("\n".join(L) + "\n")
    print("\nsaved: data/fitpool_gate.txt")
    sys.exit(0 if verdict else 1)


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'gate'
    if mode == 'worker':
        worker_main(sys.argv[2], sys.argv[3], sys.argv[4])
    elif mode == 'gate':
        gate_main()
    else:
        raise SystemExit(f"unknown mode {mode!r}")
