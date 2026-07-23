"""
STAGE 4Z (O3): the hierarchical c1 -- the Stage 4S continuous-coefficient fit
with per-galaxy disk-M/L offsets (0.1-dex lognormal prior) profiled at every
lambda node. Goal: cut the galaxy-population variance that made the 4S
bootstrap interval (+0.29/-0.25) three times wider than the profile interval.

Grid lam in [-0.2, 1.4] step 0.1, warm-started; joint SPARC+lensing
marginalized objective (4S/4E fiducial); 100-rep galaxy bootstrap with
(lam, a0, f, s_int, dlt) + offsets refit per rep.

Gates: G1 delta-prior->0 endpoint regression vs the stored 4S values
(lam=0: -8341.95, lam=1: -8397.72, +-1.0); G2 injection at lam-truth = 1
with drawn offsets recovers lam-hat in [0.8, 1.2].
Writes data/stage4z_hierc1.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25
LN10 = math.log(10)
S_ML = 0.1*LN10

meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta[t[0]] = (float(t[5]), int(t[17]))
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
ug = np.unique(gal_id)
NGal = len(ug)
gmap = {g: i for i, g in enumerate(ug)}
gidx = np.array([gmap[g] for g in gal_id])
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]

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

def m2h(th, lam, dml, w_g, pd):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm = gN*nu_lam(gN/a0, lam)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm)
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    yl = 10**lg/a0
    rl = l_gobs[lmask] - (lg + np.log10(nu_lam(yl, lam)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(pd*pd)
    return out

def fit_node(lam, w_g, pd=S_ML, th0=None, dml0=None, rounds=3,
             lg_override=None):
    global lgobs
    lg_save = lgobs
    if lg_override is not None: lgobs = lg_override
    dml = np.zeros(NGal) if dml0 is None else dml0.copy()
    best = None
    try:
        for rd in range(rounds):
            starts = ([list(best.x)] if best is not None else []) + \
                     ([list(th0)] if th0 is not None else []) + \
                     [[math.log10(A0_FID), 1.0, 0.08, 0.0]]
            bb = None
            for t0 in starts:
                b = minimize(lambda t: m2h(t, lam, dml, w_g, pd), t0,
                             method='Nelder-Mead',
                             options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
                if bb is None or b.fun < bb.fun: bb = b
            best = bb
            la0, f, s_int, dlt = best.x
            if pd < 1e-3: continue
            se2c = s_int*s_int
            for gi2 in range(NGal):
                if w_g[gi2] == 0: dml[gi2] = 0.0; continue
                mm = GIDXS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = lgobs[mm] - np.log10(gN2*nu_lam(gN2/10**la0, lam))
                    s2 = sig2[mm] + se2c
                    return (w_g[gi2]*np.sum(rr*rr/s2)
                            + w_g[gi2]*dl*dl/(pd*pd))
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
    finally:
        lgobs = lg_save
    return best, dml

ones = np.ones(NGal)
LG = np.round(np.arange(-0.2, 1.401, 0.1), 3)
L = [f"STAGE 4Z hierarchical c1: {kept} galaxies, {len(gobs)} points + "
     f"{int(lmask.sum())} lensing; delta prior 0.1 dex; grid {LG[0]}..{LG[-1]}",
     ""]

prof, th, dm = [], None, None
for lam in LG:
    b, dm = fit_node(lam, ones, th0=th, dml0=dm)
    th = b.x
    prof.append(b.fun)
prof = np.array(prof)
i = int(np.argmin(prof))
lam_hat = LG[i]
if 0 < i < len(LG)-1:
    x3, y3 = LG[i-1:i+2], prof[i-1:i+2]
    c2_, c1_, _ = np.polyfit(x3, y3, 2)
    if c2_ > 0: lam_hat = -c1_/(2*c2_)
dz = prof - prof.min()
lo = hi = None
for j in range(i, -1, -1):
    if dz[j] > 1.0:
        lo = np.interp(1.0, [dz[j+1], dz[j]], [LG[j+1], LG[j]]); break
for j in range(i, len(LG)):
    if dz[j] > 1.0:
        hi = np.interp(1.0, [dz[j-1], dz[j]], [LG[j-1], LG[j]]); break
z0 = dz[np.argmin(np.abs(LG-0.0))]
zh = dz[np.argmin(np.abs(LG-1.0))]
L.append(f"profile: lam_hat = {lam_hat:.3f} (D1 "
         f"{None if lo is None else round(lo,3)}..{None if hi is None else round(hi,3)})"
         f" -> c1_hat = {lam_hat/2:.3f}")
L.append(f"D(-2lnL): c1=0 at {z0:+.1f} | c1=1/2 at {zh:+.1f} | "
         f"edge: {'INTERIOR' if 0 < i < len(LG)-1 else 'EDGE'}")
L.append("profile table (lam, obj):")
for lam, p in zip(LG, prof):
    L.append(f"  {lam:+.2f}  {p:10.2f}")
L.append("")

# bootstrap
rng = np.random.default_rng(43)
allg = np.arange(NGal)
lam_b = []
for k in range(100):
    pick = rng.choice(allg, NGal, replace=True)
    w = np.zeros(NGal)
    for g_ in pick: w[g_] += 1
    def j5(t):
        lam = t[0]
        if not (-0.4 < lam < 1.6): return 1e12
        return m2h(t[1:], lam, dml_b, w, S_ML)
    dml_b = np.zeros(NGal)
    best = None
    for rd in range(2):
        bb = None
        for t0 in ([best.x.tolist()] if best is not None else []) + \
                  [[lam_hat, math.log10(A0_FID), 1.0, 0.08, 0.0],
                   [0.5, math.log10(A0_FID), 1.0, 0.10, 0.0]]:
            b = minimize(j5, t0, method='Nelder-Mead',
                         options=dict(maxiter=5000, xatol=1e-6, fatol=1e-6))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        lamr, la0, f, s_int, dlt = best.x
        for gi2 in range(NGal):
            if w[gi2] == 0: dml_b[gi2] = 0.0; continue
            mm = GIDXS[gi2]
            def od(dl):
                fc = f*math.exp(dl)
                gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                rr = lgobs[mm] - np.log10(gN2*nu_lam(gN2/10**la0, lamr))
                s2 = sig2[mm] + s_int*s_int
                return w[gi2]*(np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML))
            dml_b[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                         method='bounded').x
    lam_b.append(best.x[0])
lam_b = np.array(lam_b)
p16, p50, p84 = np.percentile(lam_b, [16, 50, 84])
L.append(f"bootstrap (100 reps): lam 16/50/84 = {p16:.3f}/{p50:.3f}/{p84:.3f}")
L.append(f"c1_hat (hier, bootstrap): {p50/2:.3f} +{(p84-p50)/2:.3f} "
         f"-{(p50-p16)/2:.3f}   [4S comparator: 0.427 +0.290/-0.246]")
L.append(f"P(lam>0) = {(lam_b>0).mean():.3f}; P(lam>0.5) = {(lam_b>0.5).mean():.3f}")
L.append("")

# gates
b0, _ = fit_node(0.0, ones, pd=1e-4, rounds=1)
b1, _ = fit_node(1.0, ones, pd=1e-4, rounds=1)
ok1 = abs(b0.fun-(-8341.95)) < 1.0 and abs(b1.fun-(-8397.72)) < 1.0
L.append(f"G1 prior->0 endpoints: lam=0 {b0.fun:.2f} (4S -8341.95), "
         f"lam=1 {b1.fun:.2f} (4S -8397.72) -> {'PASS' if ok1 else 'FAIL'}")
d_t = rng.normal(0, S_ML, NGal)
gN_t = g_gas + 1.0*np.exp(d_t[gidx])*g_dsk + g_bul
lg_inj = (np.log10(gN_t*nu_be(gN_t/A0_FID))
          + rng.normal(0, np.sqrt(sig2 + 0.08**2)))
def j5i(t):
    lam = t[0]
    if not (-0.4 < lam < 1.6): return 1e12
    return m2h(t[1:], lam, np.zeros(NGal), ones, S_ML)
lg_save = lgobs; lgobs = lg_inj
bi = None
for t0 in ([1.0, math.log10(A0_FID), 1.0, 0.10, 0.0],
           [0.5, math.log10(A0_FID), 1.0, 0.10, 0.0]):
    b = minimize(j5i, t0, method='Nelder-Mead',
                 options=dict(maxiter=6000, xatol=1e-6, fatol=1e-6))
    if bi is None or b.fun < bi.fun: bi = b
lgobs = lg_save
ok2 = 0.8 <= bi.x[0] <= 1.2
L.append(f"G2 injection (truth lam=1 + offsets): lam_rec = {bi.x[0]:.3f} -> "
         f"{'PASS' if ok2 else 'FAIL'}")

out = "\n".join(L)
print(out)
with open('data/stage4z_hierc1.txt', 'w') as f:
    f.write(out+"\n")
