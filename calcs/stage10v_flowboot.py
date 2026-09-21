"""10V VERIFICATION ANNEX -- the flow-leg paired bootstrap (R48 SS11
successor item 2, executed before the flow leg is quoted anywhere
again, incl. Proposal F).

Computes the refit-grade uncertainty of the flow-71 contest deltas:
200 paired galaxy-bootstrap replicates (resample + per-galaxy
distance jitter at each galaxy's own rel, the 10U G10U-5 stream
pattern, seed 202; no UMa shared draw -- the flow leg has no UMa),
all four members refit warm per replicate. OUTPUT ONLY: P(A beats B)
matrix + the d(BE,boot) distribution vs the frozen SD_block 68.24 and
the fixed-global surrogate (z -0.71/-0.78). No letter, no bar, no
credence content -- this measures an uncertainty, it tests nothing.
"""
import math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)

SRCV = open(os.path.join(HERE, 'stage10v_strat.py'),
            encoding='utf-8').read()
cutv = SRCV.index("# ---------------- baselines (both modes)")
NSV = {'__name__': 'stage10v_trunc',
       '__file__': os.path.join(HERE, 'stage10v_strat.py')}
_argv = sys.argv
sys.argv = ['stage10v_strat.py', 'gates']
exec(compile(SRCV[:cutv], 'stage10v_trunc', 'exec'), NSV)
sys.argv = _argv

W78, FLOW, FAMS = NSV['W78'], NSV['FLOW'], NSV['FAMS']
MEMBERS, ARC = NSV['MEMBERS'], NSV['ARC']
build_sub, fit_deep = NSV['build_sub'], NSV['fit_deep']

L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

t0 = time.time()
P("10V FLOW-LEG PAIRED BOOTSTRAP (verification annex; R48 SS11-2)")
P("")
subF = build_sub(W78, FLOW)
base = {}
for m in MEMBERS:
    base[m] = fit_deep(subF, FAMS[m],
                       th0=[math.log10(ARC[m]), 1.0, 0.08, 0.0])
fmin = min(base[m].fun for m in MEMBERS)
P("baseline (warm, GB R-i convention): " +
  ", ".join(f"{m} {base[m].fun-fmin:+.1f}" for m in MEMBERS))
TH0 = {m: list(base[m].x) for m in MEMBERS}

rng = np.random.default_rng(202)
funs = {m: [] for m in MEMBERS}
NREP = 200
for rep in range(NREP):
    pick = rng.choice(FLOW, size=len(FLOW), replace=True)
    docc = []
    for g in pick:
        s_g = 1.0 + rng.normal(0, W78['rel'][int(g)])
        s_g = min(max(s_g, 0.5), 1.5)
        docc.append(-math.log10(s_g))
    subr = build_sub(W78, [int(g) for g in pick], dlg_per_occ=docc)
    for m in MEMBERS:
        funs[m].append(fit_deep(subr, FAMS[m], th0=TH0[m]).fun)
    if (rep + 1) % 10 == 0:
        P(f"  rep {rep+1}/{NREP}  [{(time.time()-t0)/60:.1f} min]")
for m in MEMBERS:
    funs[m] = np.array(funs[m])

P("")
P("P(A beats B) over 200 paired replicates:")
for A in MEMBERS:
    row = []
    for B in MEMBERS:
        row.append("  .  " if A == B else
                   f"{float(np.mean(funs[B] - funs[A] > 0)):.3f}")
    P(f"  {A:5s}: " + "  ".join(f"{c:>6s}" for c in row))

d = funs['boot'] - funs['BE']
P("")
P(f"d(BE,boot) replicates: mean {d.mean():+.1f}, SD {d.std(ddof=1):.1f}, "
  f"pct16/50/84 = {np.percentile(d,16):+.1f}/{np.percentile(d,50):+.1f}/"
  f"{np.percentile(d,84):+.1f}")
d_sky = base['boot'].fun - base['BE'].fun
z = d_sky/d.std(ddof=1)
P(f"sky d = {d_sky:+.1f} -> paired-boot z = {z:+.2f} "
  f"(frozen SD_block 68.24 -> 0.71; fixed-global surrogate -0.71/-0.78)")
P(f"P(BE better than boot) = {float(np.mean(d > 0)):.3f}")
P("")
P("NOTE: uncertainty measurement only -- no letter, no bar, no "
  "credence content; quotes of the flow leg now carry THIS grade.")
P(f"wall-clock: {(time.time()-t0)/60:.1f} min")
with open('data/stage10v_flowboot.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage10v_flowboot.txt")
