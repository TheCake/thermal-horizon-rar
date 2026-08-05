"""STAGE 8S — THE GAS-DOMINATED c1 (referee-queue T5, the pincer's M/L
defense).  Pre-registered BEFORE any run.

The measured galaxy dial c1 = 0.26-0.52 (4S flat / 4Z hier) carries a
standing objection (the Hees-catalog pincer's W10): the deep-window
coefficient could be manufactured by disk M/L freedom.  Gas-dominated
galaxies are the standard immunizer: where g_gas carries the baryonic
budget, the fitted f_ML is nearly inert and c1 is pinned by physics.

Instrument (4S machinery lifted VERBATIM: the nu_lam family c1 = lam/2,
the SPARC data block, the marginalized objective): galaxies are split
by the fraction of their kept points with g_gas > g_dsk + g_bul at
f = 1 (GDFRAC >= 0.5 => gas-dominated).  Fits, all SPARC-ONLY
marginalized (lam, a0, f_ML, s_int profiled; NO lensing leg — its
0.2-dex stellar-mass systematic would contaminate the M/L-immunity
claim) + the raw y<1 co-read:
  (1) full-sample SPARC-only lam profile   (the internal reference)
  (2) gas-dominated subsample profile + 200-rep galaxy bootstrap
  (3) disk-dominated complement profile    (the contrast control)
  (4) THE IMMUNITY DEMONSTRATION: lam-hat with f_ML FORCED to 0.5 and
      2.0 on both subsamples — gas-dominated lam-hat should barely
      move; the complement shows the contrast (pre-stated band:
      |d lam-hat| <= 0.2 for GD; informational for the complement).

Gates (any FAIL => STOP; amendment pre-quote): G8S-0 selector sanity
(GD galaxy count >= 15, else POWER-STOP; counts + GDFRAC distribution
printed); G8S-1 endpoint regression — the 4S JOINT objective (lifted
verbatim, lensing included ONLY here) reproduces the 4F/4S references
at lam=1 (-8397.72) and lam=0 (-8341.95) to +-0.5; G8S-2 additivity —
the SPARC-only objective at one theta splits exactly (1e-6) into
GD + complement; G8S-3 series c1(lam) = lam/2 (the 4S G1, verbatim).

Bars (locked): T5-DEFENDED iff GD bootstrap P(lam > 0) >= 0.95 AND the
GD Delta-1 interval overlaps the full-sample interval.  T5-POWER-
LIMITED iff P(lam > 0) < 0.95 with the GD interval covering both 0 and
the full-sample lam-hat (quote the interval; no verdict on the dial).
T5-TENSION iff the GD Delta-1 interval EXCLUDES the full-sample
lam-hat.  The immunity band is co-reported (informational).
NO credence movement (measurement round; pre-stated); the galaxy-side
credence map is untouched.
Output: data/stage8s_gasc1.txt
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize
import time

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

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
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
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs = np.log10(gobs)

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_lam(y, lam):
    return (1.0-lam)*nu_standard(y) + lam*nu_be(y)

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

t00 = time.time()
P("8S THE GAS-DOMINATED c1 (T5; pre-reg committed BEFORE any run; "
  "measurement round; NO credence movement)")

# ---------------- G8S-0: the selector ----------------
allg = np.unique(gal_id)
NGAL = gal_id.max()+1
gdfrac = np.zeros(NGAL)
for g_ in allg:
    m = gal_id == g_
    gdfrac[g_] = float(np.mean(g_gas[m] > g_dsk[m] + g_bul[m]))
gd_set = np.array([g_ for g_ in allg if gdfrac[g_] >= 0.5])
dd_set = np.array([g_ for g_ in allg if gdfrac[g_] < 0.5])
n_gd, n_dd = len(gd_set), len(dd_set)
qs = np.percentile(gdfrac[allg], [10, 25, 50, 75, 90])
g0_ok = n_gd >= 15
P(f"G8S-0 selector: {kept} galaxies, {len(gobs)} points; GDFRAC "
  f"(share of points with g_gas > g_dsk+g_bul at f=1) percentiles "
  f"10/25/50/75/90 = " + "/".join(f"{v:.2f}" for v in qs)
  + f"; GAS-DOMINATED (>=0.5): {n_gd} galaxies "
  f"({int(np.isin(gal_id, gd_set).sum())} points); complement: "
  f"{n_dd} galaxies -> {'PASS' if g0_ok else 'POWER-STOP'}")

# ---------------- G8S-3: series (4S G1 verbatim) ----------------
xg = np.linspace(0.005, 0.05, 40)
yg = xg*xg
g3_ok = True
for lam in (0.0, 0.5, 1.0):
    h = nu_lam(yg, lam)*xg - 1.0
    c2n, c1n = np.polyfit(xg, h/xg, 1)
    oks = abs(c1n - lam/2) < 2e-3
    g3_ok &= oks
    P(f"G8S-3 series lam={lam:.1f}: c1 = {c1n:+.4f} (pred {lam/2:+.4f})"
      f" -> {'OK' if oks else 'FAIL'}")

# ---------------- objectives ----------------
LN10 = math.log(10)
def m2ll_joint(th, lam, w_gal):
    # 4S objective VERBATIM (SPARC + lensing) — used ONLY for G8S-1
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    gm = gN*nu_lam(gN/a0, lam)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm)
    out = np.sum(w_gal[gal_id]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    yl = 10**lg/a0
    lgm = lg + np.log10(nu_lam(yl, lam))
    rl = l_gobs[lmask] - lgm
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    return out

def m2ll_sparc(th, lam, w_gal, f_forced=None):
    la0, f, s_int = th
    if f_forced is not None:
        f = f_forced
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    gm = gN*nu_lam(gN/a0, lam)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm)
    return np.sum(w_gal[gal_id]*(r*r/se2 + np.log(se2)))

def fit_sparc_at(lam, w_gal, th_warm=None, f_forced=None):
    starts = [[math.log10(A0_FID), 1.0, 0.08],
              [math.log10(A0_FID)+0.1, 0.8, 0.12]]
    if th_warm is not None: starts.insert(0, list(th_warm))
    best = None
    for th0 in starts:
        b = minimize(lambda t: m2ll_sparc(t, lam, w_gal, f_forced), th0,
                     method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
        if best is None or b.fun < best.fun: best = b
    return best

LGRID = np.round(np.arange(-0.30, 1.501, 0.05), 3)
def profile_sparc(w_gal, f_forced=None):
    prof, th = [], None
    for lam in LGRID:
        b = fit_sparc_at(lam, w_gal, th, f_forced)
        prof.append(b.fun); th = b.x
    prof = np.array(prof)
    i = int(np.argmin(prof))
    lam_hat = LGRID[i]
    if 0 < i < len(LGRID)-1:
        x3, y3 = LGRID[i-1:i+2], prof[i-1:i+2]
        c2_, c1_, _ = np.polyfit(x3, y3, 2)
        if c2_ > 0: lam_hat = -c1_/(2*c2_)
    lo = hi = None
    for j in range(i, -1, -1):
        if prof[j] > prof[i]+1.0:
            lo = np.interp(prof[i]+1.0, [prof[j+1], prof[j]][::1],
                           [LGRID[j+1], LGRID[j]][::1]) \
                if prof[j] != prof[j+1] else LGRID[j]
            break
    for j in range(i, len(LGRID)):
        if prof[j] > prof[i]+1.0:
            hi = np.interp(prof[i]+1.0, [prof[j-1], prof[j]],
                           [LGRID[j-1], LGRID[j]])
            break
    return prof, lam_hat, lo, hi, i

def wvec(gset):
    w = np.zeros(NGAL)
    w[gset] = 1.0
    return w

# ---------------- G8S-1: joint endpoint regression ----------------
ones = np.ones(NGAL)
def fit_joint_at(lam, th_warm=None):
    starts = [[math.log10(A0_FID), 1.0, 0.08, 0.0],
              [math.log10(A0_FID)+0.1, 0.8, 0.12, -0.1]]
    if th_warm is not None: starts.insert(0, list(th_warm))
    best = None
    for th0 in starts:
        b = minimize(lambda t: m2ll_joint(t, lam, ones), th0,
                     method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
        if best is None or b.fun < best.fun: best = b
    return best
be_end = fit_joint_at(1.0).fun
std_end = fit_joint_at(0.0, None).fun
g1_ok = abs(be_end - (-8397.72)) < 0.5 and abs(std_end - (-8341.95)) < 0.5
P(f"G8S-1 joint endpoints (4S/4F refs): lam=1 {be_end:.2f} "
  f"(ref -8397.72), lam=0 {std_end:.2f} (ref -8341.95) "
  f"-> {'PASS' if g1_ok else 'FAIL'}")

# ---------------- G8S-2: additivity ----------------
th_t = [math.log10(A0_FID), 1.0, 0.08]
va = m2ll_sparc(th_t, 0.7, ones)
vg = m2ll_sparc(th_t, 0.7, wvec(gd_set))
vd = m2ll_sparc(th_t, 0.7, wvec(dd_set))
g2_ok = abs(va - (vg+vd)) < 1e-6
P(f"G8S-2 additivity at probe theta: full {va:.6f} = GD {vg:.6f} + "
  f"DD {vd:.6f} (d = {va-(vg+vd):.2e}) -> {'PASS' if g2_ok else 'FAIL'}")

if not (g0_ok and g1_ok and g2_ok and g3_ok):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    out = "\n".join(L)
    with open('data/stage8s_gasc1.txt', 'w') as f:
        f.write(out+"\n")
    raise SystemExit(0)
P("GATES: G8S-0/1/2/3 ALL PASS")
P("")

# ---------------- the three profiles ----------------
res = {}
for tag, gset in (('FULL', allg), ('GD', gd_set), ('DD', dd_set)):
    w = wvec(gset)
    prof, lam_hat, lo, hi, i = profile_sparc(w)
    dz = prof - prof.min()
    z0 = dz[np.argmin(np.abs(LGRID-0.0))]
    zh = dz[np.argmin(np.abs(LGRID-1.0))]
    edge = 'INTERIOR' if 0 < i < len(LGRID)-1 else 'EDGE'
    res[tag] = dict(lam=lam_hat, lo=lo, hi=hi, z0=z0)
    P(f"[{tag:4}] SPARC-only marginalized: lam_hat = {lam_hat:.3f} "
      f"(D1 {None if lo is None else round(lo,3)}.."
      f"{None if hi is None else round(hi,3)}) -> c1_hat = "
      f"{lam_hat/2:.3f}; D(-2lnL) at c1=0: {z0:+.1f}, at c1=1/2: "
      f"{zh:+.1f}; {edge}")

# raw y<1 co-read on GD
y_fid = (g_gas + g_dsk + g_bul)/A0_FID
sel1 = y_fid < 1.0
def fit_raw_at(lam, w_pt, th_warm=None):
    def chi2(th):
        la0, f = th
        if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
        gN = (g_gas + f*g_dsk + g_bul)[sel1]
        r = (lgobs[sel1] - np.log10(gN*nu_lam(gN/10**la0, lam)))/sig[sel1]
        return np.sum(w_pt[sel1]*r*r)
    starts = [[math.log10(A0_FID), 1.0]]
    if th_warm is not None: starts.insert(0, list(th_warm))
    best = None
    for th0 in starts:
        b = minimize(chi2, th0, method='Nelder-Mead',
                     options=dict(maxiter=2000, xatol=1e-7, fatol=1e-7))
        if best is None or b.fun < best.fun: best = b
    return best
wpt_gd = wvec(gd_set)[gal_id]
prof, th = [], None
for lam in LGRID:
    b = fit_raw_at(lam, wpt_gd, th)
    prof.append(b.fun); th = b.x
prof = np.array(prof)
i = int(np.argmin(prof))
lam_raw = LGRID[i]
if 0 < i < len(LGRID)-1:
    x3, y3 = LGRID[i-1:i+2], prof[i-1:i+2]
    c2_, c1_, _ = np.polyfit(x3, y3, 2)
    if c2_ > 0: lam_raw = -c1_/(2*c2_)
ngd1 = int((wpt_gd[sel1] > 0).sum())
P(f"[GD  ] raw y<1 co-read ({ngd1} points): lam_hat = {lam_raw:.3f} "
  f"-> c1_hat = {lam_raw/2:.3f}; Dchi2 at c1=0 = "
  f"{(prof-prof.min())[np.argmin(np.abs(LGRID-0.0))]:+.1f}")

# ---------------- the immunity demonstration ----------------
P("")
for tag, gset in (('GD', gd_set), ('DD', dd_set)):
    w = wvec(gset)
    lams = {}
    for fv in (0.5, 2.0):
        _, lam_f, _, _, _ = profile_sparc(w, f_forced=fv)
        lams[fv] = lam_f
    dl = abs(lams[2.0]-lams[0.5])
    P(f"[{tag:4}] IMMUNITY: lam_hat(f=0.5) = {lams[0.5]:.3f}, "
      f"lam_hat(f=2.0) = {lams[2.0]:.3f}, |d lam| = {dl:.3f}"
      + ("  (band <= 0.2)" if tag == 'GD' else "  (contrast control)"))
    if tag == 'GD':
        gd_immune = dl <= 0.2

# ---------------- GD bootstrap ----------------
rng = np.random.default_rng(41)
def joint4_sparc(th, w_gal):
    lam = th[0]
    if not (-0.4 < lam < 1.6): return 1e12
    return m2ll_sparc(th[1:], lam, w_gal)
lam_b = []
for k in range(200):
    pick = rng.choice(gd_set, len(gd_set), replace=True)
    wg = np.zeros(NGAL)
    for g_ in pick: wg[g_] += 1
    best = None
    for th0 in ([1.0, math.log10(A0_FID), 1.0, 0.08],
                [0.4, math.log10(A0_FID), 0.9, 0.10]):
        b = minimize(joint4_sparc, th0, args=(wg,), method='Nelder-Mead',
                     options=dict(maxiter=6000, xatol=1e-6, fatol=1e-6))
        if best is None or b.fun < best.fun: best = b
    lam_b.append(best.x[0])
lam_b = np.array(lam_b)
p16, p50, p84 = np.percentile(lam_b, [16, 50, 84])
pgt0 = float((lam_b > 0).mean())
P("")
P(f"[GD  ] bootstrap (200 reps over {n_gd} galaxies): lam 16/50/84 = "
  f"{p16:.3f}/{p50:.3f}/{p84:.3f}; P(lam>0) = {pgt0:.3f}; "
  f"P(lam>0.5) = {float((lam_b>0.5).mean()):.3f}")
P(f"[GD  ] c1 (bootstrap): {p50/2:.3f} +{(p84-p50)/2:.3f} "
  f"-{(p50-p16)/2:.3f}")

# ---------------- verdict ----------------
P("")
full_lo, full_hi = res['FULL']['lo'], res['FULL']['hi']
gd_lo, gd_hi = res['GD']['lo'], res['GD']['hi']
overlap = not (gd_hi is not None and full_lo is not None
               and gd_hi < full_lo) \
      and not (gd_lo is not None and full_hi is not None
               and gd_lo > full_hi)
excl_full = (gd_lo is not None and gd_hi is not None
             and (res['FULL']['lam'] < gd_lo or res['FULL']['lam'] > gd_hi))
if pgt0 >= 0.95 and overlap:
    P(f"==> 8S VERDICT (locked grammar): T5-DEFENDED — the "
      f"gas-dominated subsample alone excludes c1 = 0 at "
      f"P(lam>0) = {pgt0:.3f} with its interval overlapping the "
      f"full-sample dial; the deep-window coefficient is NOT a disk-M/L "
      f"artifact (immunity |d lam| "
      + ("<= 0.2 shown" if gd_immune else "band MISSED, disclosed")
      + ").")
elif excl_full:
    P("==> 8S VERDICT (locked grammar): T5-TENSION — the gas-dominated "
      "interval excludes the full-sample dial; the subsample disagrees "
      "with the population (investigate before any dial quote).")
else:
    P(f"==> 8S VERDICT (locked grammar): T5-POWER-LIMITED — "
      f"P(lam>0) = {pgt0:.3f} < 0.95; the gas-dominated subsample "
      f"cannot alone exclude c1 = 0; interval quoted, dial unchanged.")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t00)/60:.1f} min)")

out = "\n".join(L)
with open('data/stage8s_gasc1.txt', 'w') as f:
    f.write(out+"\n")
print("\nsaved: data/stage8s_gasc1.txt")
