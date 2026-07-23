"""
STAGE 4S: the zero-point coefficient MEASURED -- c1 promoted from a branch
contest (4B/4F discrete families) to a continuous parameter with an error bar.

Family: nu_lam(y) = (1-lam)*nu_standard(y) + lam*nu_BE(y) -- a one-parameter
slice through coefficient space with c1(lam) = lam/2 exactly (G1 verifies this
numerically; c2(lam) rides along and is reported at the fit). lam=0 is the
c1=0 branch (standard-mu), lam=1 the RAR/occupation law. The 1/4-branch's c1
sits at lam=1/2 -- same c1, different c2 than 4F's bootstrap function: the
slice pins c1, not the full function, and is labeled as such.

Both treatments, per the paper's dual-likelihood discipline:
 (a) raw chi2, (a0,f_ML) profiled per lam; y<1 window (4B's strongest) + y<30;
 (b) scatter-marginalized joint SPARC+lensing (4E fiducial objective;
     (a0,f_ML,s_int,dlt) profiled per lam).
Errors: profile-likelihood Delta=1 intervals + 200-rep galaxy bootstrap
(full joint refit with lam free per rep).

Gates: G1 series c1(lam)=lam/2 at lam=0,1/2,1; G2 endpoint regression vs the
4F fiducial -2lnL (BE -8397.72, standard -8341.95, +-0.5); G3 lam-hat interior.
Writes data/stage4s_c1fit.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25

# ---------------- SPARC + lensing (identical to 4E/4F) ----------------
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

# ---------------- the lam family ----------------
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_lam(y, lam):
    return (1.0-lam)*nu_standard(y) + lam*nu_be(y)

L = [f"STAGE 4S continuous-c1 fit: {kept} galaxies, {len(gobs)} SPARC points "
     f"+ {int(lmask.sum())} lensing points (4E fiducial config)",
     "family: nu_lam = (1-lam)*standard + lam*BE;  c1 = lam/2"]

# ---------------- G1: numerical c1(lam) ----------------
xg = np.linspace(0.005, 0.05, 40)
yg = xg*xg
g1ok = True
L.append("")
for lam in (0.0, 0.5, 1.0):
    h = nu_lam(yg, lam)*xg - 1.0
    c2n, c1n = np.polyfit(xg, h/xg, 1)
    ok = abs(c1n - lam/2) < 2e-3
    g1ok &= ok
    L.append(f"G1 series lam={lam:.1f}: c1 = {c1n:+.4f} (pred {lam/2:+.4f}), "
             f"c2 = {c2n:+.4f} -> {'OK' if ok else 'FAIL'}")

# ---------------- objectives (verbatim 4E/4F) ----------------
LN10 = math.log(10)
def m2ll(th, lam, w_gal):
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

def fit_joint_at(lam, w_gal, th_warm=None):
    starts = [[math.log10(A0_FID), 1.0, 0.08, 0.0],
              [math.log10(A0_FID)+0.1, 0.8, 0.12, -0.1]]
    if th_warm is not None: starts.insert(0, list(th_warm))
    best = None
    for th0 in starts:
        b = minimize(lambda t: m2ll(t, lam, w_gal), th0, method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
        if best is None or b.fun < best.fun: best = b
    return best

def fit_raw_at(lam, sel, w_pt, th_warm=None):
    def chi2(th):
        la0, f = th
        if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
        gN = (g_gas + f*g_dsk + g_bul)[sel]
        r = (lgobs[sel] - np.log10(gN*nu_lam(gN/10**la0, lam)))/sig[sel]
        return np.sum(w_pt[sel]*r*r)
    starts = [[math.log10(A0_FID), 1.0]]
    if th_warm is not None: starts.insert(0, list(th_warm))
    best = None
    for th0 in starts:
        b = minimize(chi2, th0, method='Nelder-Mead',
                     options=dict(maxiter=2000, xatol=1e-7, fatol=1e-7))
        if best is None or b.fun < best.fun: best = b
    return best

y_fid = (g_gas + g_dsk + g_bul)/A0_FID
NG = gal_id.max()+1
allg = np.unique(gal_id)
ones = np.ones(NG)
LGRID = np.round(np.arange(-0.30, 1.501, 0.05), 3)

def profile(fit_fn):
    prof, th = [], None
    for lam in LGRID:
        b = fit_fn(lam, th)
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
                           [LGRID[j+1], LGRID[j]][::1]) if prof[j] != prof[j+1] else LGRID[j]
            break
    for j in range(i, len(LGRID)):
        if prof[j] > prof[i]+1.0:
            hi = np.interp(prof[i]+1.0, [prof[j-1], prof[j]],
                           [LGRID[j-1], LGRID[j]])
            break
    return prof, lam_hat, lo, hi, i

# ---------------- (a) raw treatment ----------------
L.append("")
for ymax, tag in ((1.0, 'y<1'), (30.0, 'y<30')):
    sel = y_fid < ymax
    prof, lam_hat, lo, hi, i = profile(
        lambda lam, th: fit_raw_at(lam, sel, np.ones(len(gobs)), th))
    dz = prof - prof.min()
    z0  = dz[np.argmin(np.abs(LGRID-0.0))]
    zq  = dz[np.argmin(np.abs(LGRID-0.5))]
    zh  = dz[np.argmin(np.abs(LGRID-1.0))]
    L.append(f"(a) RAW [{tag}]: lam_hat = {lam_hat:.3f} "
             f"(D1 interval {lo if lo is None else round(lo,3)}"
             f"..{hi if hi is None else round(hi,3)}) -> "
             f"c1_hat = {lam_hat/2:.3f}")
    L.append(f"    Dchi2 at c1=0: {z0:+.1f} | c1=1/4: {zq:+.1f} | "
             f"c1=1/2: {zh:+.1f}   (edge: {'INTERIOR' if 0 < i < len(LGRID)-1 else 'EDGE'})")

# raw bootstrap (y<1): 3-param (lam, la0, f) per rep
rng = np.random.default_rng(41)
sel1 = y_fid < 1.0
def raw3(th, w_pt):
    lam, la0, f = th
    if not (-0.4 < lam < 1.6): return 1e12
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    gN = (g_gas + f*g_dsk + g_bul)[sel1]
    r = (lgobs[sel1] - np.log10(gN*nu_lam(gN/10**la0, lam)))/sig[sel1]
    return np.sum(w_pt[sel1]*r*r)
lam_raw_b = []
for k in range(200):
    pick = rng.choice(allg, len(allg), replace=True)
    wg = np.zeros(NG)
    for g_ in pick: wg[g_] += 1
    w_pt = wg[gal_id]
    best = None
    for th0 in ([1.0, math.log10(A0_FID), 1.0], [0.4, math.log10(A0_FID), 0.9]):
        b = minimize(raw3, th0, args=(w_pt,), method='Nelder-Mead',
                     options=dict(maxiter=3000, xatol=1e-6, fatol=1e-6))
        if best is None or b.fun < best.fun: best = b
    lam_raw_b.append(best.x[0])
lam_raw_b = np.array(lam_raw_b)
L.append(f"    raw bootstrap (y<1, 200 reps): lam 16/50/84 = "
         f"{np.percentile(lam_raw_b,16):.3f}/{np.percentile(lam_raw_b,50):.3f}/"
         f"{np.percentile(lam_raw_b,84):.3f}; P(lam>0.5) = "
         f"{(lam_raw_b>0.5).mean():.3f}; P(lam>0) = {(lam_raw_b>0).mean():.3f}")

# ---------------- (b) marginalized treatment ----------------
prof, lam_hat, lo, hi, i = profile(lambda lam, th: fit_joint_at(lam, ones, th))
dz = prof - prof.min()
z0 = dz[np.argmin(np.abs(LGRID-0.0))]
zq = dz[np.argmin(np.abs(LGRID-0.5))]
zh = dz[np.argmin(np.abs(LGRID-1.0))]
be_end  = prof[np.argmin(np.abs(LGRID-1.0))]
std_end = prof[np.argmin(np.abs(LGRID-0.0))]
g2 = abs(be_end - (-8397.72)) < 0.5 and abs(std_end - (-8341.95)) < 0.5
L += ["",
      f"(b) MARGINALIZED joint: lam_hat = {lam_hat:.3f} "
      f"(D1 interval {lo if lo is None else round(lo,3)}"
      f"..{hi if hi is None else round(hi,3)}) -> c1_hat = {lam_hat/2:.3f}",
      f"    D(-2lnL) at c1=0: {z0:+.1f} | c1=1/4: {zq:+.1f} | c1=1/2: {zh:+.1f}"
      f"   (edge: {'INTERIOR' if 0 < i < len(LGRID)-1 else 'EDGE'})",
      f"    G2 endpoints: lam=1 {be_end:.2f} (ref -8397.72), "
      f"lam=0 {std_end:.2f} (ref -8341.95) -> {'OK' if g2 else 'FAIL'}"]
# c2 at the fit
h = nu_lam(yg, min(max(lam_hat,0),1.5))*xg - 1.0
c2n, c1n = np.polyfit(xg, h/xg, 1)
L.append(f"    at lam_hat: c1 = {c1n:+.4f}, c2 = {c2n:+.4f} "
         f"(BE c2 = +0.0833, simple +0.1250, boot +0.0729)")

# joint bootstrap: 5-param (lam, la0, f, s_int, dlt) per rep
def joint5(th, w_gal):
    lam = th[0]
    if not (-0.4 < lam < 1.6): return 1e12
    return m2ll(th[1:], lam, w_gal)
lam_j_b = []
for k in range(200):
    pick = rng.choice(allg, len(allg), replace=True)
    wg = np.zeros(NG)
    for g_ in pick: wg[g_] += 1
    best = None
    for th0 in ([1.0, math.log10(A0_FID), 1.0, 0.08, 0.0],
                [0.4, math.log10(A0_FID), 0.9, 0.10, -0.05]):
        b = minimize(joint5, th0, args=(wg,), method='Nelder-Mead',
                     options=dict(maxiter=6000, xatol=1e-6, fatol=1e-6))
        if best is None or b.fun < best.fun: best = b
    lam_j_b.append(best.x[0])
lam_j_b = np.array(lam_j_b)
L.append(f"    joint bootstrap (200 reps): lam 16/50/84 = "
         f"{np.percentile(lam_j_b,16):.3f}/{np.percentile(lam_j_b,50):.3f}/"
         f"{np.percentile(lam_j_b,84):.3f}; P(lam>0.5) = {(lam_j_b>0.5).mean():.3f}; "
         f"P(lam>0) = {(lam_j_b>0).mean():.3f}")
L.append(f"    c1_hat (joint, bootstrap): "
         f"{np.percentile(lam_j_b,50)/2:.3f} "
         f"+{(np.percentile(lam_j_b,84)-np.percentile(lam_j_b,50))/2:.3f} "
         f"-{(np.percentile(lam_j_b,50)-np.percentile(lam_j_b,16))/2:.3f}")

L.append("")
L.append("profile tables (lam, raw chi2 y<1, marginalized -2lnL):")
sel = y_fid < 1.0
pr_raw, _, _, _, _ = profile(lambda lam, th: fit_raw_at(lam, sel, np.ones(len(gobs)), th))
for lam, a_, b_ in zip(LGRID, pr_raw, prof):
    L.append(f"  {lam:+.2f}  {a_:9.1f}  {b_:10.2f}")

out = "\n".join(L)
print(out)
with open('data/stage4s_c1fit.txt', 'w') as f:
    f.write(out+"\n")
