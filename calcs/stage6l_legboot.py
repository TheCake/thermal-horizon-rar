"""
STAGE 6L (the round's missing error bar): bootstrap grade of the
leg-count measurement.

6I measured L = 2 against both flanking integers (vertical, fiducial
gate s = 0.8681): L1 -52.76 / L2 -59.05 / L3 -54.52 vs BE, i.e.
d12 = +6.29 and d32 = +4.53 in L2's favor. Point estimates only.
This stage grades them with the 6C/6J machinery (40 paired
galaxy-resample reps, dv-ON vertical, lensing noise redraws, rng 53):
per rep, all three L-laws are fit warm-lite and the paired deltas
recorded. Report: mean +/- sd of d12, d32; count of reps where L2 is
strictly best. Full-fit regression gates vs the 6I absolutes (<3).
Writes data/stage6l_legboot.txt.
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

def s_of(e):
    x = math.sqrt(e)
    n = 1.0/(math.exp(x) - 1.0)
    return n/(1.0 + n)

S_GLOB = s_of(0.02)

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

def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))

def nu_L(y, Lc):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    sL = S_GLOB**Lc
    ly = np.log(y)
    nu = nu_simple(y)
    for _ in range(80):
        d1 = 2.0*nu - 1.0
        b = 0.5*sL/d1**Lc
        db = -Lc*sL/d1**(Lc + 1)
        u = np.exp(np.minimum(0.5*(1.0+b)*ly + b*np.log(nu), 60.0))
        eu = np.exp(np.minimum(u, 60.0))
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        F = nu - 1.0 - n
        dudnu = u*(db*(0.5*ly + np.log(nu)) + b/nu)
        dF = 1.0 + (eu/(em1*em1))*dudnu
        step = F/np.where(np.abs(dF) > 1e-12, dF, 1e-12)
        nu = np.maximum(nu - np.clip(step, -0.5*nu, 0.5*nu), 1.0 + 1e-15)
    return nu

def m2hv(th, Lc, dml, dv, w_g, lens_obs):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*nu_L(gN/a0, Lc)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm_) - dv[gidx]
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    rl = lens_obs[lmask] - (lg + np.log10(nu_L(10**lg/a0, Lc)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    out += np.sum(w_g*dv*dv/(SIGV*SIGV))
    return out

def fit_v(Lc, w_g, lens_obs, th0=None, dml0=None, dv0=None,
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
            b = minimize(lambda t: m2hv(t, Lc, dml, dv, w_g, lens_obs), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int, dlt = best.x
        se2c = s_int*s_int
        for _ in range(sweeps):
            fac = f*np.exp(dml[gidx])
            gN = g_gas + fac*g_dsk + g_bul
            r0 = lgobs - np.log10(gN*nu_L(gN/10**la0, Lc))
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
                    rr = (lgobs[mm] - np.log10(gN2*nu_L(gN2/10**la0, Lc))
                          - dv[gi2])
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2hv(best.x, Lc, dml, dv, w_g, lens_obs)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2hv(t, Lc, dml, dv, w_g, lens_obs),
                 list(best.x), method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml, dv

REF6I = {1: -12205.26, 2: -12211.54, 3: -12207.01}
ones = np.ones(NGal)
L = [f"STAGE 6L leg-count bootstrap: {kept} galaxies, {len(gobs)} points"
     f" + {int(lmask.sum())} lensing; 40 paired reps; gate s = "
     f"{S_GLOB:.4f}; 6I point estimates d12 = +6.29, d32 = +4.53"]

FULL = {}
for Lc in (2, 1, 3):
    th0 = FULL[2][0].x if Lc != 2 and 2 in FULL else None
    dm0 = FULL[2][1] if Lc != 2 and 2 in FULL else None
    dv0 = FULL[2][2] if Lc != 2 and 2 in FULL else None
    bb, dm, dvv = fit_v(Lc, ones, l_gobs, th0=th0, dml0=dm0, dv0=dv0)
    FULL[Lc] = (bb, dm, dvv)
    d = bb.fun - REF6I[Lc]
    L.append(f"full-fit regression L{Lc}: {bb.fun:.2f} (6I {REF6I[Lc]:.2f},"
             f" d={d:+.2f}) -> {'PASS' if abs(d) < 3.0 else 'FAIL'}")
    print(L[-1], flush=True)
    assert abs(d) < 3.0

rng = np.random.default_rng(53)
allg = np.arange(NGal)
d12s, d32s, best2 = [], [], 0
for k in range(40):
    pick = rng.choice(allg, NGal, replace=True)
    w = np.zeros(NGal)
    for g_ in pick: w[g_] += 1
    lo = l_gobs + rng.normal(0, np.sqrt(l_sig2))
    vals = {}
    for Lc in (1, 2, 3):
        bb, _, _ = fit_v(Lc, w, lo, th0=FULL[Lc][0].x, dml0=FULL[Lc][1],
                         dv0=FULL[Lc][2], tol=0.5, max_rounds=5, sweeps=2)
        vals[Lc] = bb.fun
    d12s.append(vals[1] - vals[2])
    d32s.append(vals[3] - vals[2])
    if vals[2] < vals[1] and vals[2] < vals[3]: best2 += 1
    if (k+1) % 10 == 0:
        print(f"  rep {k+1}/40: d12 {np.mean(d12s):+.1f}, d32 "
              f"{np.mean(d32s):+.1f}, L2-best {best2}/{k+1}", flush=True)
d12s, d32s = np.array(d12s), np.array(d32s)
L.append(f"bootstrap (40 reps): d12 = L1 - L2 = {d12s.mean():+.2f} +/- "
         f"{d12s.std(ddof=1):.2f} (L2 better in {int((d12s>0).sum())}/40)")
L.append(f"                     d32 = L3 - L2 = {d32s.mean():+.2f} +/- "
         f"{d32s.std(ddof=1):.2f} (L2 better in {int((d32s>0).sum())}/40)")
L.append(f"L2 strictly best of three: {best2}/40")
L.append(f"percentiles d12 16/50/84: "
         f"{np.percentile(d12s,[16,50,84]).round(1).tolist()}; d32: "
         f"{np.percentile(d32s,[16,50,84]).round(1).tolist()}")

out = "\n".join(L)
print(out)
with open('data/stage6l_legboot.txt', 'w') as f:
    f.write(out+"\n")
print("\nsaved: data/stage6l_legboot.txt")
