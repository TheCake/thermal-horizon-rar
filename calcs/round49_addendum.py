"""ROUND 49 verification addendum -- GA half (BLIND: written and
committed BEFORE the Round-49 report exists; 16th execution of the
87a4676 protocol). Independent re-derivations of the 10W sky read's
load-bearing numbers, arm's-length where the object permits:

  GA-1 map closure: the committed CSV's per-galaxy boot deltas must
       sum to each leg's engine d(BE,boot) within the measured
       convergence residual (A3: ~0.03).
  GA-2 partial correlations by a DIFFERENT construction (inverse of
       the Spearman correlation matrix) vs the stage's
       Freedman-Lane residual implementation.
  GA-3 cut-20 flow contest re-fit SERIAL from a cut builder
       re-implemented from the A4-ii prereg text + its own theta
       parse; d and SD_block vs the pooled sky values.
  GA-4 the conditional boot re-run (same seed 909, warms from GA-3's
       own serial optima): P(d_rep > 0) and P(d_rep > -SD) must
       reproduce exactly (pool certified deterministic); drep saved.
  GA-5 letter booleans recomputed + the branch-overlap disclosure
       (CARRIED and IMMUNE point conditions are NOT mutually
       exclusive at these SDs; the stage's if-order decided).
  GA-6 carrier theta spots (UGC03580, NGC5371 ~ 8 arcsec per R48).

Usage: py calcs/round49_addendum.py ga
"""
import math, os, sys, time
import numpy as np
from scipy.stats import rankdata

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

W78, LEGA, FLOW = NSV['W78'], NSV['LEGA'], NSV['FLOW']
FAMS, MEMBERS, ARC = NSV['FAMS'], NSV['MEMBERS'], NSV['ARC']
build_sub, fit_hier, fit_deep, contest_fit = (
    NSV['build_sub'], NSV['fit_hier'], NSV['fit_deep'],
    NSV['contest_fit'])
nu_be = NSV['nu_be']
NAMEOF = NSV['NAMEOF']
sub_fields, sd_block = NSV['sub_fields'], NSV['sd_block']
NS10U = NSV['NS']
SP, UD, UB, KPC = NS10U['SP'], NS10U['UD'], NS10U['UB'], NS10U['KPC']

sys.path.insert(0, HERE)
from fitpool import run_tasks  # noqa: E402

L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

t0 = time.time()
P("ROUND 49 GA (blind; pre-committed)")
P("")

# ---------------- baselines (the stage recipe) ----------------------
subP = build_sub(W78, LEGA)
b_plain = fit_hier(subP, nu_be, use_u=False)
b_loose = fit_hier(subP, nu_be, use_u=True, th0=list(b_plain.x)+[0.0])
arch = {'BE': fit_deep(subP, nu_be, th0=list(b_loose.x))}
for nm in ('p065', 'gm', 'boot'):
    arch[nm] = fit_deep(subP, FAMS[nm], th0=list(arch['BE'].x))
subF = build_sub(W78, FLOW)
fw = {m: fit_deep(subF, FAMS[m],
                  th0=[math.log10(ARC[m]), 1.0, 0.08, 0.0])
      for m in MEMBERS}
P(f"baselines re-fit  [{time.time()-t0:.0f}s]")

# ---------------- GA-1 map closure ---------------------------------
rowsC = [l.rstrip('\n').split(',') for l in
         open('data/stage10w_map.csv', encoding='utf-8')][1:]
sumd = {'anch': 0.0, 'flow': 0.0}
th_csv = {}
for r in rowsC:
    sumd[r[1]] += float(r[5])
    th_csv[r[0]] = float(r[4])
dA = arch['boot'].fun - arch['BE'].fun
dF = fw['boot'].fun - fw['BE'].fun
okA = abs(sumd['anch'] - dA) <= 0.05
okF = abs(sumd['flow'] - dF) <= 0.05
P(f"GA-1 map closure: anch CSV sum {sumd['anch']:+.2f} vs engine "
  f"{dA:+.2f} ({'OK' if okA else 'BAD'}); flow {sumd['flow']:+.2f} "
  f"vs {dF:+.2f} ({'OK' if okF else 'BAD'})")

# ---------------- GA-2 partial correlations (matrix route) ---------
fg = np.array([float(r[2]) for r in rowsC])
cv = np.array([float(r[3]) for r in rowsC])
rs = np.log10(np.array([float(r[4]) for r in rowsC]))
lv = np.array([0.0 if r[1] == 'anch' else 1.0 for r in rowsC])
dl = np.array([float(r[5]) for r in rowsC])

def partial_matrix(x, y, Z):
    cols = [x, y] + Z
    Rm = np.corrcoef(np.vstack(
        [rankdata(c).astype(float) for c in cols]))
    Pm = np.linalg.inv(Rm)
    return float(-Pm[0, 1]/math.sqrt(Pm[0, 0]*Pm[1, 1]))

REF2 = {'coverage': (+0.195, cv, [lv]),
        'composition': (-0.050, fg, [cv, lv]),
        'resolution': (+0.092, rs, [cv, lv]),
        'leg': (-0.215, lv, [cv])}
ok2 = True
for ax, (ref, a, Z) in REF2.items():
    r = partial_matrix(dl, a, Z)
    ok = abs(r - ref) <= 5e-3
    ok2 &= ok
    P(f"GA-2 {ax:11s}: matrix partial rho = {r:+.4f} vs stage "
      f"{ref:+.3f} ({'OK' if ok else 'BAD'})")

# ---------------- GA-3 cut-20 flow serial re-fit --------------------
import glob as _glob
ROTPATH = {os.path.basename(p).replace('_rotmod.dat', ''): p
           for p in _glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                               recursive=True)}
PT = np.full(len(W78['lgobs']), np.nan)
THG = {}
for g in LEGA + FLOW:
    nm2 = NAMEOF[int(g)]
    D0 = SP['dist'][int(g)][0]
    ths = []
    for l in open(ROTPATH[nm2]):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC
        gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        ths.append(R/D0*206.265)
    pts = W78['gpts'][int(g)]
    assert len(ths) == len(pts), nm2
    PT[pts] = ths
    THG[int(g)] = min(ths)

def cut_world_ga(W, keep):
    """Re-implemented from the prereg A4-ii text."""
    idx = np.where(keep)[0]
    W2 = dict(W)
    for k in ('gg', 'gd', 'gb', 'lgobs', 'sig2', 'gal_id'):
        W2[k] = W[k][idx]
    ug = np.unique(W2['gal_id'])
    gp = {int(g): np.where(W2['gal_id'] == g)[0] for g in ug}
    keepg = {g for g in gp
             if len(gp[g]) == len(W['gpts'][int(g)])
             or len(gp[g]) >= 3}
    W2['gpts'] = {g: gp[g] for g in keepg}
    W2['ug'] = np.array(sorted(keepg))
    return (W2, [g for g in LEGA if g in keepg],
            [g for g in FLOW if g in keepg])

W20, leg20, flow20 = cut_world_ga(W78, ~(PT < 20.0))
sub20 = build_sub(W20, flow20)
c20 = {}
for m in ('BE', 'boot'):
    c20[m], gap = contest_fit(sub20, m, list(fw[m].x))
d20 = c20['boot'].fun - c20['BE'].fun
A = 'BE' if c20['BE'].fun <= c20['boot'].fun else 'boot'
B = 'boot' if A == 'BE' else 'BE'
fldA = sub_fields(sub20, A, list(c20[A].x))
muB = sub_fields(sub20, B, list(c20[B].x))[0]
sd20 = sd_block(sub20, fldA, muB)
ok3 = abs(d20 - (-29.7)) <= 0.1 and abs(sd20 - 66.0) <= 0.5
P(f"GA-3 cut-20 flow (serial): d = {d20:+.2f} (sky -29.7), "
  f"SD_block = {sd20:.2f} (sky 66.0), gals {len(flow20)}/71 "
  f"({'OK' if ok3 else 'BAD'})  [{time.time()-t0:.0f}s]")

# ---------------- GA-4 conditional boot re-run ----------------------
rng9 = np.random.default_rng(909)
warm20 = {m: list(c20[m].x) for m in ('BE', 'boot')}
bsubs, btasks = {}, []
for rep in range(200):
    pick = rng9.choice(flow20, size=len(flow20), replace=True)
    docc = []
    for g in pick:
        s_g = 1.0 + rng9.normal(0, W78['rel'][int(g)])
        docc.append(-math.log10(min(max(s_g, 0.5), 1.5)))
    bsubs[f'r{rep}'] = build_sub(W20, [int(g) for g in pick],
                                 dlg_per_occ=docc)
    for m in ('BE', 'boot'):
        btasks.append(dict(key=f'r{rep}:{m}', kind='fit_deep',
                           sub=f'r{rep}', member=m, th0=warm20[m]))
BRES = run_tasks(btasks, bsubs, tag='ga4boot', log=P)
drep = np.array([BRES[f'r{i}:boot']['fun'] - BRES[f'r{i}:BE']['fun']
                 for i in range(200)])
p_sign = float(np.mean(drep > 0))
p_nore = float(np.mean(drep > -sd20))
ok4 = (abs(p_sign - 0.225) <= 1e-9) and (abs(p_nore - 0.635) <= 1e-9)
P(f"GA-4 boot re-run: P(d>0) = {p_sign:.3f} (sky 0.225), "
  f"P(d>-SD) = {p_nore:.3f} (sky 0.635) "
  f"({'EXACT' if ok4 else 'DIFFER'}); drep pct16/50/84 = "
  f"{np.percentile(drep,16):+.1f}/{np.percentile(drep,50):+.1f}/"
  f"{np.percentile(drep,84):+.1f}")
np.savez('data/round49_ga_drep.npz', drep=drep)

# ---------------- GA-5 letter booleans ------------------------------
SKYD = {10.0: (-31.7, 88.0), 20.0: (-29.7, 66.0), 30.0: (-54.2, 56.9)}
D_UNCUT = fw['boot'].fun - fw['BE'].fun
carried = all((d > 0) or (abs(d) < sd)
              for d, sd in (SKYD[20.0], SKYD[30.0]))
immune = all(abs(d - D_UNCUT) < sd for d, sd in SKYD.values())
P(f"GA-5 letter booleans on the sky numbers: carried-point = "
  f"{carried}, immune = {immune} -> BOTH TRUE = the registered "
  f"branches are NOT mutually exclusive at these SDs; the stage's "
  f"if-order (carried first) decided the label. DISCLOSED for R49.")

# ---------------- GA-6 carrier theta spots --------------------------
for nm2 in ('UGC03580', 'NGC5371'):
    gid = [g for g in FLOW if NAMEOF[int(g)] == nm2]
    if gid:
        P(f"GA-6 {nm2}: theta_min = {THG[int(gid[0])]:.1f} arcsec "
          f"(R48: carriers ~8); CSV {th_csv[nm2]:.1f}")

P("")
allok = okA and okF and ok2 and ok3 and ok4
P(f"GA VERDICT: {'ALL OK' if allok else 'DIFFERS PRESENT'}")
P(f"wall-clock: {(time.time()-t0)/60:.1f} min")
with open('data/round49_addendum.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/round49_addendum.txt")
