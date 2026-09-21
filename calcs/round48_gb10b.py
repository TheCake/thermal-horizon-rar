"""GB-10b -- corrected frozen-bulge probe (adjudication of the GB-10
DIFFER): the first GB attempt used gd+gb (bulge overweighted by the
UB/UD ratio 1.4) and warm-only fits; the reviewer's recipe converts
the bulge to the disk Upsilon scale (gb/0.7*0.5) and uses the
contest-grade estimator. This run matches his recipe exactly, in the
stage's own inherited engine."""
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
build_sub, contest_fit = NSV['build_sub'], NSV['contest_fit']

t0 = time.time()
subF = build_sub(W78, FLOW)
WARM = {'BE': [math.log10(1.1132e-10), 1.426, 0.0509, 0.0],
        'boot': [math.log10(1.1221e-10), 1.747, 0.0498, 0.0]}
base = {m: contest_fit(subF, m, WARM[m])[0] for m in ('BE', 'boot')}
d_base = base['BE'].fun - base['boot'].fun
print(f"baseline flow d = {d_base:+.2f}  [{time.time()-t0:.0f}s]",
      flush=True)
subB = dict(subF)
subB['gd'] = subF['gd'] + subF['gb']/0.7*0.5
subB['gb'] = np.zeros_like(subF['gb'])
fb = {m: contest_fit(subB, m, WARM[m])[0] for m in ('BE', 'boot')}
d_free = fb['BE'].fun - fb['boot'].fun
print(f"bulge-M/L-freed flow d = {d_free:+.2f} (rev +37.93; "
      f"a0 BE {10**fb['BE'].x[0]:.3e}, f {fb['BE'].x[1]:.3f})  "
      f"[{time.time()-t0:.0f}s]", flush=True)
ok = abs(d_free - 37.93) <= 3.0
line = (f"GB 10b (corrected recipe gb/0.7*0.5, contest grade): "
        f"baseline {d_base:+.2f}, freed {d_free:+.2f} "
        f"(rev 48.73 -> 37.93, -22%) -> "
        f"{'AGREE' if ok else 'DIFFER'}; the first GB-10 attempt is "
        f"RETIRED (wrong bulge scale gd+gb = x1.4 overweight + "
        f"warm-only starts)")
print(line, flush=True)
with open('data/round48_addendum.txt', 'a', encoding='utf-8') as f:
    f.write("\n" + line + "\n")
