"""ROUND 49 verification -- GB half (post-report): re-compute every
load-bearing REVIEWER number in own implementations before adoption
(the standing verify-reviewer-math rule). Appends to
data/round49_addendum.txt.

Legs:
  GB-1 the engine-stall gaps on cut-flow worlds (his error #2):
       indep-profiled totals at the stage optima vs engine funs;
       profiled-at-stage-x d for 20/30 vs his -33.99/-56.21.
  GB-2 SD_block(flow,20) with deep-profiled fields vs his 64.30;
       the NGC5985 dml-wall diagnosis.
  GB-3 the population-fixed control (20-arcsec points on the
       30-arcsec 64-galaxy population) vs his -30.28.
  GB-4 symmetrised partials: leg|cov,res and res|cov (his -0.133
       p 0.107 / +0.192 p 0.021), FL p at own seed.
  GB-5 the joint single-normalisation re-run (149-galaxy fit,
       per-galaxy deltas, coverage/leg axes; his coverage p 0.0072,
       leg p 0.070, map rank-corr 0.926).
  GB-6 the PAIRED co-requirement (his C13 construction, own seed
       911): P(Delta > 0) vs his 0.620.
  GB-7 dml box-pinning census at the cut-20 flow optima (his 6/69
       pinned carrying 84% of the Delta-sum; NGC5985 -0.432/+0.700;
       IC2574 pinned both members anchored).
  GB-8 estimator-box inspection (done by source read: trigger
       unreachable CONFIRMED) + H0(anch20, BE) vs his 64.54.
  GB-9 the inner-minus-outer residual contrast at 20 arcsec
       (his anch -0.086 +/- 0.030, flow -0.004 +/- 0.018).

Usage: py calcs/round49_gb.py
"""
import math, os, sys, time
import numpy as np
from scipy.optimize import minimize_scalar as msc
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
S_ML, U_PRIOR, LN10 = NSV['S_ML'], NSV['U_PRIOR'], NSV['LN10']
NAMEOF = NSV['NAMEOF']
sub_fields, sd_block = NSV['sub_fields'], NSV['sd_block']
h0_of_a0 = NSV['h0_of_a0']
NS10U = NSV['NS']
SP, UD, UB, KPC = NS10U['SP'], NS10U['UD'], NS10U['UB'], NS10U['KPC']

sys.path.insert(0, HERE)
from fitpool import run_tasks  # noqa: E402

L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

t0 = time.time()
P("")
P("=" * 68)
P("ROUND 49 GB (post-report; reviewer numbers re-computed)")
P("")

# ---------------- shared: baselines + theta + cut worlds -----------
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
WARM_A = {m: list(arch[m].x) for m in MEMBERS}
WARM_F = {m: list(fw[m].x) for m in MEMBERS}
P(f"baselines  [{time.time()-t0:.0f}s]")

import glob as _glob
ROTPATH = {os.path.basename(p).replace('_rotmod.dat', ''): p
           for p in _glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                               recursive=True)}
PT = np.full(len(W78['lgobs']), np.nan)
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
    PT[W78['gpts'][int(g)]] = ths

def cut_world(W, keep):
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

W20, leg20, flow20 = cut_world(W78, ~(PT < 20.0))
W30, leg30, flow30 = cut_world(W78, ~(PT < 30.0))
sub20 = build_sub(W20, flow20)
sub30 = build_sub(W30, flow30)

# the stage's cut fits (deterministic re-derivation, pooled)
tasks, subs = [], {'f20': sub20, 'f30': sub30,
                   'a20': build_sub(W20, leg20)}
for m in MEMBERS:
    tasks.append(dict(key=f'f20:{m}', kind='contest_fit', sub='f20',
                      member=m, th0=WARM_F[m]))
    tasks.append(dict(key=f'f30:{m}', kind='contest_fit', sub='f30',
                      member=m, th0=WARM_F[m]))
tasks.append(dict(key='a20:BE', kind='contest_fit', sub='a20',
                  member='BE', th0=WARM_A['BE']))
CR = run_tasks(tasks, subs, tag='gbcuts', log=P)

# ---------------- GB-1 engine-stall gaps ---------------------------
def indep_pergal(sub, nm, th, cap=None):
    la0, f, s_int, u = th
    a0 = 10**la0
    nu = FAMS[nm]
    per, dmls = [], []
    for k in range(sub['n']):
        pts = sub['gpts'][k]
        lg = sub['lg'][pts]
        gg = sub['gg'][pts]; gd = sub['gd'][pts]; gb = sub['gb'][pts]
        s2 = sub['s2'][pts]
        sv = sub['sv'][k]
        isu = sub['isuma_g'][k]
        se2 = s2 + s_int*s_int
        dml, dv = 0.0, 0.0
        prev = None
        for _ in range(cap or 400):
            gN = gg + f*math.exp(dml)*gd + gb
            r0 = lg - np.log10(gN*nu(gN/a0)) - u*isu
            w = 1.0/se2
            dv = float(np.sum(w*r0)/(np.sum(w) + 1.0/sv**2))

            def od(dl):
                gN2 = gg + f*math.exp(dl)*gd + gb
                rr = lg - np.log10(gN2*nu(gN2/a0)) - dv - u*isu
                return float(np.sum(rr*rr/se2) + dl*dl/(S_ML*S_ML))
            dml = msc(od, bounds=(-0.7, 0.7), method='bounded').x
            gN = gg + f*math.exp(dml)*gd + gb
            r = lg - np.log10(gN*nu(gN/a0)) - dv - u*isu
            cur = (float(np.sum(r*r/se2 + np.log(se2)))
                   + dml*dml/(S_ML*S_ML) + dv*dv/(sv*sv))
            if prev is not None and abs(prev - cur) < 1e-10:
                prev = cur
                break
            prev = cur
        per.append(prev)
        dmls.append(float(dml))
    return per, sum(per) + (u/U_PRIOR)**2, dmls

P("")
P("-- GB-1 engine-stall gaps (indep-profiled vs engine, cut flow) --")
IND = {}
for w, sub in (('f20', sub20), ('f30', sub30)):
    for m in MEMBERS:
        per, tot, dmls = indep_pergal(sub, m, CR[f'{w}:{m}']['x'])
        IND[f'{w}:{m}'] = (per, tot, dmls)
        gap = CR[f'{w}:{m}']['fun'] - tot
        P(f"  {w}/{m}: engine - profiled = {gap:+.2f} "
          f"(rev claims 5.5-11.8 member-dependent)")
d20p = IND['f20:boot'][1] - IND['f20:BE'][1]
d30p = IND['f30:boot'][1] - IND['f30:BE'][1]
P(f"  profiled-at-stage-x d: 20 = {d20p:+.2f} (rev -33.99), "
  f"30 = {d30p:+.2f} (rev -56.21)")
ord20 = {m: IND[f'f20:{m}'][1] for m in MEMBERS}
o20 = min(ord20.values())
P("  profiled 20-arcsec ordering: " +
  ", ".join(f"{m} {ord20[m]-o20:+.1f}" for m in MEMBERS) +
  "  (rev: BE +34.0, p065 +0.9, gm/boot tie)")

# ---------------- GB-2 SD_block(flow,20) deep -----------------------
P("")
P("-- GB-2 SD_block(flow,20) at deep-profiled fields --")
def deep_fields(sub, nm, th, iters=3000):
    la0, f, s_int, u = th
    a0 = 10**la0
    nu = FAMS[nm]
    npt = len(sub['lg'])
    MU = np.empty(npt); RW = np.empty(npt); RES = np.empty(npt)
    nonconv = []
    for k in range(sub['n']):
        pts = sub['gpts'][k]
        lg = sub['lg'][pts]
        gg = sub['gg'][pts]; gd = sub['gd'][pts]; gb = sub['gb'][pts]
        s2 = sub['s2'][pts]
        sv = sub['sv'][k]
        isu = sub['isuma_g'][k]
        se2 = s2 + s_int*s_int
        dml, dv = 0.0, 0.0
        conv = False
        for _ in range(iters):
            gN = gg + f*math.exp(dml)*gd + gb
            r0 = lg - np.log10(gN*nu(gN/a0)) - u*isu
            w = 1.0/se2
            dvn = float(np.sum(w*r0)/(np.sum(w) + 1.0/sv**2))

            def od(dl):
                gN2 = gg + f*math.exp(dl)*gd + gb
                rr = lg - np.log10(gN2*nu(gN2/a0)) - dvn - u*isu
                return float(np.sum(rr*rr/se2) + dl*dl/(S_ML*S_ML))
            dmln = msc(od, bounds=(-0.7, 0.7), method='bounded').x
            if abs(dmln - dml) < 1e-12 and abs(dvn - dv) < 1e-12:
                dml, dv = dmln, dvn
                conv = True
                break
            dml, dv = dmln, dvn
        if not conv:
            nonconv.append(k)
        gN = gg + f*math.exp(dml)*gd + gb
        mu = np.log10(gN*nu(gN/a0)) + dv + u*isu
        MU[pts] = mu; RW[pts] = 1.0/se2; RES[pts] = lg - mu
    return (MU, RW, RES), nonconv

r20 = {m: CR[f'f20:{m}'] for m in MEMBERS}
A = 'BE' if r20['BE']['fun'] <= r20['boot']['fun'] else 'boot'
B = 'boot' if A == 'BE' else 'BE'
fldA_e = sub_fields(sub20, A, r20[A]['x'])
muB_e = sub_fields(sub20, B, r20[B]['x'])[0]
sd_e = sd_block(sub20, fldA_e, muB_e)
fldA_d, nc1 = deep_fields(sub20, A, r20[A]['x'])
muB_d, nc2 = deep_fields(sub20, B, r20[B]['x'])
sd_d = sd_block(sub20, fldA_d, muB_d[0])
P(f"  engine-fields SD = {sd_e:.2f} (stage 66.0); deep-profiled SD "
  f"= {sd_d:.2f} (rev 64.30); non-converged galaxies at 200 iters: "
  f"engine-side n/a, deep {len(nc1)}/{len(nc2)}")
k5985 = flow20.index([g for g in flow20
                      if NAMEOF[int(g)] == 'NGC5985'][0])
for m in ('BE', 'boot'):
    dml_g = IND[f'f20:{m}'][2][k5985]
    P(f"  NGC5985 dml({m}) at 20-arcsec optimum = {dml_g:+.3f} "
      f"(rev: BE -0.432, boot +0.700 = the box wall)")

# ---------------- GB-3 population-fixed control ---------------------
P("")
P("-- GB-3 the 20-arcsec cut on the 30-arcsec population --")
subPF = build_sub(W20, flow30)          # 64 galaxies, 20-arcsec points
cpf = {m: contest_fit(subPF, m, WARM_F[m])[0] for m in ('BE', 'boot')}
dpf = cpf['boot'].fun - cpf['BE'].fun
P(f"  d(64 gal, 20-arcsec points) = {dpf:+.2f} (rev -30.28; stage "
  f"20-arcsec 69-gal -29.70, 30-arcsec 64-gal -54.16)")
P(f"  -> population share of the 20->30 move = "
  f"{abs(dpf-(-29.70)):.2f}, resolution share = "
  f"{abs(-54.16-dpf):.2f} (rev 0.58 / 23.88)")

# ---------------- GB-4 symmetrised partials -------------------------
P("")
P("-- GB-4 symmetrised partial correlations --")
rowsC = [l.rstrip('\n').split(',') for l in
         open('data/stage10w_map.csv', encoding='utf-8')][1:]
fg = np.array([float(r[2]) for r in rowsC])
cv = np.array([float(r[3]) for r in rowsC])
rs = np.log10(np.array([float(r[4]) for r in rowsC]))
lv = np.array([0.0 if r[1] == 'anch' else 1.0 for r in rowsC])
dl = np.array([float(r[5]) for r in rowsC])

def _resid(v, Z):
    A_ = (np.column_stack([np.ones(len(v)), Z]) if Z.shape[1]
          else np.ones((len(v), 1)))
    beta, *_ = np.linalg.lstsq(A_, v, rcond=None)
    return v - A_ @ beta

def partial_p(delta, axis, Zcols, rng, nperm=10000):
    rd = rankdata(delta).astype(float)
    ra = rankdata(axis).astype(float)
    rZ = (np.column_stack([rankdata(z).astype(float) for z in Zcols])
          if Zcols else np.empty((len(delta), 0)))
    ed, ea = _resid(rd, rZ), _resid(ra, rZ)
    den_a = math.sqrt(float(np.sum(ea*ea)))
    obs = float(np.sum(ed*ea)) / (math.sqrt(float(np.sum(ed*ed)))
                                  * den_a)
    hits = 0
    for _ in range(nperm):
        e = ed[rng.permutation(len(ed))]
        r = float(np.sum(e*ea)) / (math.sqrt(float(np.sum(e*e)))
                                   * den_a)
        if abs(r) >= abs(obs) - 1e-15:
            hits += 1
    return obs, (1 + hits)/(1 + nperm)

rngb = np.random.default_rng(4951)
r_leg, p_leg = partial_p(dl, lv, [cv, rs], rngb)
r_res, p_res = partial_p(dl, rs, [cv], rngb)
P(f"  leg | cov,res: rho = {r_leg:+.3f}, p = {p_leg:.3f} "
  f"(rev -0.133 / 0.107)")
P(f"  res | cov    : rho = {r_res:+.3f}, p = {p_res:.3f} "
  f"(rev +0.192 / 0.021)")

# ---------------- GB-5 joint single-normalisation -------------------
P("")
P("-- GB-5 the joint 149-galaxy re-run --")
subJ = build_sub(W78, LEGA + FLOW)
cj = {}
for m in ('BE', 'boot'):
    cj[m], _ = contest_fit(subJ, m,
                           [math.log10(ARC[m]), 1.0, 0.08, 0.0])
perJ_BE, _, _ = indep_pergal(subJ, 'BE', list(cj['BE'].x))
perJ_bt, _, _ = indep_pergal(subJ, 'boot', list(cj['boot'].x))
dJ = np.array(perJ_bt) - np.array(perJ_BE)
rc = float(np.corrcoef(rankdata(dJ), rankdata(dl))[0, 1])
r_cvJ, p_cvJ = partial_p(dJ, cv, [lv], rngb)
r_lgJ, p_lgJ = partial_p(dJ, lv, [cv], rngb)
P(f"  joint-fit deltas vs pooled map: rank-corr = {rc:.3f} "
  f"(rev 0.926)")
P(f"  coverage | leg (joint): rho = {r_cvJ:+.3f}, p = {p_cvJ:.4f} "
  f"(rev +0.222 / 0.0072)")
P(f"  leg | cov (joint): rho = {r_lgJ:+.3f}, p = {p_lgJ:.3f} "
  f"(rev p 0.070)")
P(f"  [{time.time()-t0:.0f}s]")

# ---------------- GB-6 the PAIRED co-requirement --------------------
P("")
P("-- GB-6 the paired co-requirement (own construction, seed 911) --")
common = [g for g in FLOW if g in set(flow20)]
rng11 = np.random.default_rng(911)
bsubs, btasks = {}, []
for rep in range(200):
    pick = rng11.choice(common, size=len(common), replace=True)
    docc = []
    for g in pick:
        s_g = 1.0 + rng11.normal(0, W78['rel'][int(g)])
        docc.append(-math.log10(min(max(s_g, 0.5), 1.5)))
    pick = [int(g) for g in pick]
    bsubs[f'u{rep}'] = build_sub(W78, pick, dlg_per_occ=docc)
    bsubs[f'c{rep}'] = build_sub(W20, pick, dlg_per_occ=docc)
    for m in ('BE', 'boot'):
        btasks.append(dict(key=f'u{rep}:{m}', kind='fit_deep',
                           sub=f'u{rep}', member=m, th0=WARM_F[m]))
        btasks.append(dict(key=f'c{rep}:{m}', kind='fit_deep',
                           sub=f'c{rep}', member=m,
                           th0=list(r20[m]['x'])))
BR = run_tasks(btasks, bsubs, tag='gb6paired', log=P)
dU = np.array([BR[f'u{i}:boot']['fun'] - BR[f'u{i}:BE']['fun']
               for i in range(200)])
dC = np.array([BR[f'c{i}:boot']['fun'] - BR[f'c{i}:BE']['fun']
               for i in range(200)])
DD = dC - dU
P(f"  Delta = d_cut - d_uncut: sky {(-29.70) - (-48.73):+.2f}; "
  f"replicates mean {DD.mean():+.2f}, SD {DD.std(ddof=1):.2f} "
  f"(rev +11.35 / 37.62)")
P(f"  P(Delta > 0) = {float(np.mean(DD > 0)):.3f} (rev 0.620; "
  f"bar 0.95)")

# ---------------- GB-7 dml box-pinning census -----------------------
P("")
P("-- GB-7 dml box-pinning census --")
for m in ('BE', 'boot'):
    per, tot, dmls = IND[f'f20:{m}']
    dml_arr = np.array(dmls)
    pinned = np.abs(dml_arr) >= 0.6999
    perBE = np.array(IND['f20:BE'][0])
    perBT = np.array(IND['f20:boot'][0])
    dsum = perBT - perBE
    share = float(np.sum(np.abs(dsum[pinned]))
                  / np.sum(np.abs(dsum))) if pinned.any() else 0.0
    P(f"  f20/{m}: pinned {int(pinned.sum())}/{len(dmls)}; "
      f"|Delta|-share of pinned rows = {share:.0%} "
      f"(rev: 6/69 carrying 84%)")
ic = [g for g in LEGA if NAMEOF[int(g)] == 'IC2574']
kic = LEGA.index(ic[0])
for m in ('BE', 'boot'):
    perA, totA, dmlsA = indep_pergal(subP, m, WARM_A[m])
    P(f"  IC2574 dml({m}, anchored baseline) = {dmlsA[kic]:+.3f} "
      f"(rev: pinned -0.700 both members)")

# ---------------- GB-8 H0 lever -------------------------------------
P("")
h0a20 = h0_of_a0(10**CR['a20:BE']['x'][0])
P(f"-- GB-8: C15 trigger outside the estimator box CONFIRMED by "
  f"source read (f in (0.3,2.5), s_int in [1e-3,0.4), dml +-0.7); "
  f"H0(anch 20-arcsec, BE) = {h0a20:.2f} (rev 64.54)")

# ---------------- GB-9 inner-minus-outer contrast -------------------
P("")
P("-- GB-9 inner-minus-outer residual contrast at 20 arcsec --")
for legnm, sub, seq, opt in (('anch', subP, LEGA, arch),
                             ('flow', subF, FLOW, fw)):
    A2 = 'BE' if opt['BE'].fun <= opt['boot'].fun else 'boot'
    fld = sub_fields(sub, A2, list(opt[A2].x))
    res = fld[2]
    diffs = []
    for k, g in enumerate(seq):
        pts = sub['gpts'][k]
        th_g = PT[W78['gpts'][int(g)]]
        inner = res[pts][th_g < 20.0]
        outer = res[pts][th_g >= 20.0]
        if len(inner) and len(outer):
            diffs.append(float(inner.mean() - outer.mean()))
    d_ = np.array(diffs)
    P(f"  {legnm}: mean {d_.mean():+.4f} +/- "
      f"{d_.std(ddof=1)/math.sqrt(len(d_)):.4f} dex over {len(d_)} "
      f"galaxies (rev anch -0.086 +/- 0.030, flow -0.004 +/- 0.018)")

P("")
P(f"GB complete; wall-clock {(time.time()-t0)/60:.1f} min")
with open('data/round49_addendum.txt', 'a', encoding='utf-8') as f:
    f.write("\n".join(L) + "\n")
print("\nappended: data/round49_addendum.txt")
