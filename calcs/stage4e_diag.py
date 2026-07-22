"""
STAGE 4E diagnostics v2 (autopsy before any claim is logged):
 D1: G2 with the correct 4B window (y_fid < 30) -- confirm the +54 offset was
     the 7 dropped y>=30 points, not a data/code drift.
 D2: per-regime decomposition of Delta(-2lnL) at the fiducial joint best fits,
     with COMMON bin membership (binned by fiducial y = gN(f=1)/1.2e-10 so the
     same points sit in the same bins for every family; v1 binned by each
     family's own y and points migrated between columns). Columns now sum
     exactly to the SPARC-only part of the fit deltas.
 D3: deep-window (y<0.5) fits: parameters + common-bin decomposition +
     SPARC-only totals (reconciliation).
 D4: per-galaxy influence on (a) the fiducial BE-simple delta, (b) the
     deep-regime (y_fid<0.1) BE-simple delta -- few-galaxy fragility check
     on BOTH the BE joint lead and the simple deep-regime lead.
Writes data/stage4e_diag.txt.
"""
import glob, math, os, re
import numpy as np
from scipy.optimize import minimize

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25

meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        name, inc, q = t[0], float(t[5]), int(t[17])
        meta[name] = (inc, q)
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id, gal_names = [], [], [], [], [], [], []
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q = meta.get(name, (0, 3))
    if inc < 30 or q > 2: continue
    kept += 1
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC; gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        g_gas.append(gg); g_dsk.append(gd); g_bul.append(gb)
        gobs.append(Vo*Vo/R*KPC)
        sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi); gal_names.append(name)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs = np.log10(gobs)
id2name = {}
for i, n in zip(gal_id, gal_names): id2name[int(i)] = n

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask_fid = l_gbar >= LENS_CUT_FID

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
FAMS = {'BE': nu_be, 'simple': nu_simple, 'standard': nu_standard}

def m2ll(th, nu, sel, w_gal, lmask, lens_obs, use_lens=True):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    gN = (g_gas + f*g_dsk + g_bul)[sel]
    gm = gN*nu(gN/a0)
    se2 = sig2[sel] + s_int*s_int
    r = lgobs[sel] - np.log10(gm)
    out = np.sum(w_gal[gal_id][sel]*(r*r/se2 + np.log(se2)))
    if use_lens:
        lg = l_gbar[lmask] + dlt
        yl = 10**lg/a0
        lgm = lg + np.log10(nu(yl))
        rl = lens_obs[lmask] - lgm
        out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
        out += (dlt/DELTA_PRIOR)**2
    return out

def fit(nu, sel, w_gal, lmask, lens_obs):
    best = None
    for th0 in ([math.log10(A0_FID), 1.0, 0.08, 0.0],
                [math.log10(A0_FID)+0.1, 0.8, 0.12, -0.1]):
        b = minimize(lambda t: m2ll(t, nu, sel, w_gal, lmask, lens_obs),
                     th0, method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
        if best is None or b.fun < best.fun: best = b
    return best

def perpoint(nu, th):
    la0, f, s_int, dlt = th
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gN*nu(gN/a0))
    return r*r/se2 + np.log(se2)

L = ["STAGE 4E DIAGNOSTICS v2 (common-bin membership by y_fid)"]
allsel = np.ones(len(gobs), bool)
ones = np.ones(gal_id.max()+1)
y_fid = (g_gas + g_dsk + g_bul)/A0_FID

# ---- D1
def chi2_4b(nu, sel):
    def c2(th):
        la0, f = th
        gN = (g_gas + f*g_dsk + g_bul)[sel]
        r = (lgobs[sel] - np.log10(gN*nu(gN/10**la0)))/sig[sel]
        return np.sum(r*r)
    return minimize(c2, [math.log10(A0_FID), 1.0], method='Nelder-Mead',
                    options=dict(maxiter=2000, xatol=1e-7, fatol=1e-7)).fun
stored = {}
for m in re.finditer(r'\[all y\]\s+(\w+): chi2=\s*([\d.]+)',
                     open('data/stage4b_branch.txt').read()):
    stored[m.group(1)] = float(m.group(2))
L += ["", f"D1: 4B regression, sel = y_fid<30 (n={int((y_fid<30).sum())}):"]
for name, nu in FAMS.items():
    c = chi2_4b(nu, y_fid < 30.0)
    L.append(f"  {name:>8}: {c:.1f} vs stored {stored[name]:.1f} -> "
             f"{'OK' if abs(c-stored[name]) < 1.0 else 'FAIL'}")

# ---- D2
L += ["", "D2: fiducial joint fits; Delta(-2lnL) by COMMON y_fid bin:"]
best = {n: fit(nu, allsel, ones, lmask_fid, l_gobs) for n, nu in FAMS.items()}
pp = {n: perpoint(FAMS[n], best[n].x) for n in FAMS}
sp_tot = {n: m2ll(best[n].x, FAMS[n], allsel, ones, lmask_fid, l_gobs,
                  use_lens=False) for n in FAMS}
edges = [0, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 1e9]
L.append("  y_fid bin        n     BE-simple   BE-standard")
sums = np.zeros(2)
for lo_, hi_ in zip(edges[:-1], edges[1:]):
    m = (y_fid >= lo_) & (y_fid < hi_)
    dbs = pp['BE'][m].sum() - pp['simple'][m].sum()
    dbt = pp['BE'][m].sum() - pp['standard'][m].sum()
    sums += (dbs, dbt)
    L.append(f"  [{lo_:7.2f},{hi_:7.2f}) {m.sum():5d}  {dbs:+10.2f}  {dbt:+11.2f}")
L.append(f"  column sums      {int(len(y_fid)):5d}  {sums[0]:+10.2f}  "
         f"{sums[1]:+11.2f}")
L.append(f"  SPARC-only fit deltas:      {sp_tot['BE']-sp_tot['simple']:+10.2f}"
         f"  {sp_tot['BE']-sp_tot['standard']:+11.2f}  (must match sums)")
L.append(f"  full joint deltas:          "
         f"{best['BE'].fun-best['simple'].fun:+10.2f}  "
         f"{best['BE'].fun-best['standard'].fun:+11.2f}")

# ---- D3
L += ["", "D3: y_fid<0.5 window joint fits:"]
sel05 = y_fid < 0.5
best05 = {n: fit(nu, sel05, ones, lmask_fid, l_gobs) for n, nu in FAMS.items()}
for name, b in best05.items():
    la0, f, s_int, dlt = b.x
    flag = " <-- f_ML AT BOUND" if f > 2.45 or f < 0.35 else ""
    L.append(f"  {name:>8}: -2lnL={b.fun:9.2f} a0={10**la0:.3e} f_ML={f:.2f} "
             f"s_int={s_int:.3f} dlt={dlt:+.3f}{flag}")
pp05 = {n: perpoint(FAMS[n], best05[n].x) for n in FAMS}
sp05 = {n: m2ll(best05[n].x, FAMS[n], sel05, ones, lmask_fid, l_gobs,
                use_lens=False) for n in FAMS}
edges2 = [0, 0.01, 0.03, 0.1, 0.2, 0.35, 0.5]
L.append("  y_fid bin        n     BE-simple   BE-standard")
sums = np.zeros(2)
for lo_, hi_ in zip(edges2[:-1], edges2[1:]):
    m = sel05 & (y_fid >= lo_) & (y_fid < hi_)
    dbs = pp05['BE'][m].sum() - pp05['simple'][m].sum()
    dbt = pp05['BE'][m].sum() - pp05['standard'][m].sum()
    sums += (dbs, dbt)
    L.append(f"  [{lo_:7.2f},{hi_:7.2f}) {m.sum():5d}  {dbs:+10.2f}  {dbt:+11.2f}")
L.append(f"  column sums      {int(sel05.sum()):5d}  {sums[0]:+10.2f}  "
         f"{sums[1]:+11.2f}")
L.append(f"  SPARC-only fit deltas:      {sp05['BE']-sp05['simple']:+10.2f}"
         f"  {sp05['BE']-sp05['standard']:+11.2f}  (must match sums)")
L.append(f"  full joint deltas:          "
         f"{best05['BE'].fun-best05['simple'].fun:+10.2f}  "
         f"{best05['BE'].fun-best05['standard'].fun:+11.2f}")

# ---- D4
L += ["", "D4a: per-galaxy contribution to fiducial BE-simple (all points):"]
dfull = pp['BE'] - pp['simple']
per_gal = {}
for gi in np.unique(gal_id):
    per_gal[int(gi)] = dfull[gal_id == gi].sum()
diffs = sorted((v, k) for k, v in per_gal.items())
L.append("  5 most BE-favoring:  " + "  ".join(
    f"{id2name[gi]} {d:+.2f}" for d, gi in diffs[:5]))
L.append("  5 most simple-favoring:  " + "  ".join(
    f"{id2name[gi]} {d:+.2f}" for d, gi in diffs[-5:]))
for drop_n in (1, 3):
    w_drop = ones.copy()
    for d, gi in diffs[:drop_n]: w_drop[gi] = 0.0
    rB = fit(nu_be, allsel, w_drop, lmask_fid, l_gobs).fun
    rS = fit(nu_simple, allsel, w_drop, lmask_fid, l_gobs).fun
    L.append(f"  refit dropping top-{drop_n} BE-favoring: Delta(BE-simple) = "
             f"{rB-rS:+.2f}  (was {best['BE'].fun-best['simple'].fun:+.2f})")

L += ["", "D4b: per-galaxy contribution to the DEEP-regime (y_fid<0.1) "
      "BE-simple delta (at fiducial best fits):"]
mdeep = y_fid < 0.1
ddeep = (pp['BE'] - pp['simple'])*mdeep
per_gal_d = {}
for gi in np.unique(gal_id[mdeep]):
    per_gal_d[int(gi)] = ddeep[gal_id == gi].sum()
diffs_d = sorted((v, k) for k, v in per_gal_d.items())
tot_d = sum(v for v in per_gal_d.values())
L.append(f"  deep-regime total = {tot_d:+.2f} over {int(mdeep.sum())} points, "
         f"{len(per_gal_d)} galaxies")
L.append("  5 most BE-favoring:  " + "  ".join(
    f"{id2name[gi]} {d:+.2f}" for d, gi in diffs_d[:5]))
L.append("  5 most simple-favoring:  " + "  ".join(
    f"{id2name[gi]} {d:+.2f}" for d, gi in diffs_d[-5:]))
top_share = sum(d for d, _ in diffs_d[-5:])/tot_d if tot_d != 0 else 0
L.append(f"  top-5 simple-favoring share of deep total: {100*top_share:.0f}%")

out = "\n".join(L)
print(out)
with open('data/stage4e_diag.txt', 'w') as f_:
    f_.write(out+"\n")
