"""
STAGE 5N (O11b): galaxy bootstrap of the surviving tail vote.

5M left Delta(gm - BE) = -40.6 with the measured vertical channel active --
a point estimate. The 5D-grade galaxy-population sigma was ~65, so no
significance may be quoted without a bootstrap. This runs 40 paired
galaxy-resample reps of the 5M machinery (delta_d NM-profiled + delta_v
closed-form, joint lensing with noise redraws), each function warm-started
from its full-fit optimum (globals + offsets), adaptive-lite (tol 0.5,
max 5 rounds, 2 sweeps/round).

Output: mean +- sd of (gm - BE), sign count, percentiles.
Writes data/stage5n_dvboot.txt.
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
        meta[t[0]] = (float(t[5]), int(t[17]), float(t[2]), float(t[3]),
                      float(t[6]))
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
sigv_g_map = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 10.0, 1.0, 3.0))
    if inc < 30 or q > 2: continue
    kept += 1
    irad = math.radians(inc)
    sv = math.hypot((eD/max(D, 1e-3))/LN10,
                    2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10)
    sigv_g_map[gi] = max(sv, 0.01)
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
SIGV = np.array([sigv_g_map[g] for g in ug])
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def nu_gm(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    a = y**0.75
    w = np.sqrt(nu_simple(y))
    for _ in range(30):
        u = np.minimum(a*w, 60.0)
        eu = np.exp(u)
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        H = w*w - 1.0 - n
        dH = 2.0*w + a*eu/(em1*em1)
        w = np.maximum(w - H/dH, 1e-8)
    return w*w

def m2hv(th, nu, dml, dv, w_g, lens_obs):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*nu(gN/a0)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm_) - dv[gidx]
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    rl = lens_obs[lmask] - (lg + np.log10(nu(10**lg/a0)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    out += np.sum(w_g*dv*dv/(SIGV*SIGV))
    return out

def fit_v(nu, w_g, lens_obs, th0=None, dml0=None, dv0=None,
          tol=0.05, max_rounds=15, sweeps=3):
    dml = np.zeros(NGal) if dml0 is None else dml0.copy()
    dv = np.zeros(NGal) if dv0 is None else dv0.copy()
    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08, 0.0]]
                  if rd == 0 and th0 is None else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2hv(t, nu, dml, dv, w_g, lens_obs), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int, dlt = best.x
        se2c = s_int*s_int
        for _ in range(sweeps):
            fac = f*np.exp(dml[gidx])
            gN = g_gas + fac*g_dsk + g_bul
            r0 = lgobs - np.log10(gN*nu(gN/10**la0))
            for gi2 in range(NGal):
                if w_g[gi2] == 0: dv[gi2] = 0.0; continue
                mm = GIDXS[gi2]
                w = 1.0/(sig2[mm] + se2c)
                dv[gi2] = np.sum(w*r0[mm])/(np.sum(w) + 1.0/SIGV[gi2]**2)
            for gi2 in range(NGal):
                if w_g[gi2] == 0: dml[gi2] = 0.0; continue
                mm = GIDXS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = (lgobs[mm] - np.log10(gN2*nu(gN2/10**la0))
                          - dv[gi2])
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2hv(best.x, nu, dml, dv, w_g, lens_obs)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2hv(t, nu, dml, dv, w_g, lens_obs),
                 list(best.x), method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml, dv

ones = np.ones(NGal)
L = [f"STAGE 5N dv-ON bootstrap: {kept} galaxies, {len(gobs)} points + "
     f"{int(lmask.sum())} lensing; 40 paired reps, warm-started lite"]

# full-fit warm anchors (regression vs 5M)
fb, dmlb, dvb = fit_v(nu_be, ones, l_gobs)
fg, dmlg, dvg = fit_v(nu_gm, ones, l_gobs)
d_full = fg.fun - fb.fun
ok = abs(fb.fun-(-12152.49)) < 1.5 and abs(fg.fun-(-12193.04)) < 1.5
L.append(f"full-fit regression: BE {fb.fun:.2f} (5M -12152.49), gm "
         f"{fg.fun:.2f} (5M -12193.04), d = {d_full:+.2f} -> "
         f"{'PASS' if ok else 'FAIL'}")

rng = np.random.default_rng(53)
allg = np.arange(NGal)
ds = []
for k in range(40):
    pick = rng.choice(allg, NGal, replace=True)
    w = np.zeros(NGal)
    for g_ in pick: w[g_] += 1
    lo = l_gobs + rng.normal(0, np.sqrt(l_sig2))
    vb, _, _ = fit_v(nu_be, w, lo, th0=fb.x, dml0=dmlb, dv0=dvb,
                     tol=0.5, max_rounds=5, sweeps=2)
    vg, _, _ = fit_v(nu_gm, w, lo, th0=fg.x, dml0=dmlg, dv0=dvg,
                     tol=0.5, max_rounds=5, sweeps=2)
    ds.append(vg.fun - vb.fun)
    if (k+1) % 10 == 0:
        print(f"  rep {k+1}/40: running mean {np.mean(ds):+.1f}", flush=True)
ds = np.array(ds)
L.append(f"bootstrap (40 reps): gm - BE = {ds.mean():+.2f} +/- "
         f"{ds.std(ddof=1):.2f}  (gm better in {int((ds<0).sum())}/40)")
L.append(f"percentiles 16/50/84: "
         f"{np.percentile(ds,[16,50,84]).round(1).tolist()}")

out = "\n".join(L)
print(out)
with open('data/stage5n_dvboot.txt', 'w') as f:
    f.write(out+"\n")
print("\nsaved: data/stage5n_dvboot.txt")
