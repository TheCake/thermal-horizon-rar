"""
STAGE 4W (O1): the full-hierarchy second moment. Per-galaxy nuisances:
  delta_d  disk-M/L multiplier (lognormal prior, 0.1 dex)
  delta_b  bulge-M/L multiplier (0.1 dex; only galaxies with bulge light)
  delta_v  vertical offset in dex (distance + inclination; per-galaxy prior
           sigma_v measured from the SPARC table: (e_D/D)/ln10 (+) 2*cot(i)*e_i/ln10
           -- both errors move g_obs vertically at fixed g_bar under our loader)
Scatter models recompete on the residual:
  M0w const | M1bw oscillator+floor | M2w free 6-bin | M3w area-scaled+floor.
Questions: does the thermal trend survive the FULL hierarchy? does the 4U
V-shape flatten (bulge arm) once bulge M/L and vertical offsets are absorbed?
Gates:
  G1 prior->0 reproduces Stage 4T M0/M1b within 1.0.
  G2 CALIBRATED injection: inject delta_d (0.1 dex) + delta_v (per-galaxy
     sigma_v); expected recovery slope estimated from a pure-noise run
     (shrinkage-aware); PASS if |slope_obs - slope_exp| < 0.25.
Writes data/stage4w_fullhier.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)
S_ML = 0.1*LN10                      # ln-prior on delta_d, delta_b

meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        name = t[0]
        D, eD = float(t[2]), float(t[3])
        inc, einc = float(t[5]), float(t[6])
        q = int(t[17])
        meta[name] = (inc, q, D, eD, einc)
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc = [], [], [], [], [], [], []
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
        gal_id.append(gi); R_kpc.append(R)
g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc))
sig2 = sig*sig
lgobs = np.log10(gobs)

ug = np.unique(gal_id)
NGal = len(ug)
gmap = {g: i for i, g in enumerate(ug)}
gidx = np.array([gmap[g] for g in gal_id])
SIGV = np.array([sigv_g_map[g] for g in ug])
HASB = np.array([bool(np.any(g_bul[gidx == i] > 0)) for i in range(NGal)])
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]

dR = np.zeros(len(R_kpc))
for gi in range(NGal):
    m = GIDXS[gi]; Rg = R_kpc[m]
    d = np.zeros(len(Rg))
    if len(Rg) == 1: d[0] = Rg[0]
    else:
        d[1:-1] = 0.5*(Rg[2:]-Rg[:-2]); d[0] = Rg[1]-Rg[0]; d[-1] = Rg[-1]-Rg[-2]
    dR[m] = np.maximum(d, 1e-3)
AREA = 2*np.pi*R_kpc*dR

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def osc_shape(x):
    xc = np.minimum(np.clip(x, 1e-9, None), 80.0)
    n = 1.0/np.expm1(xc)
    return np.sqrt(n/(1.0+n))

x_fid = np.sqrt((g_gas + g_dsk + g_bul)/A0_FID)
QE = np.quantile(x_fid, np.linspace(0, 1, 7))
QE[0], QE[-1] = 0.0, np.inf
BIN_FID = np.clip(np.searchsorted(QE, x_fid, side='right')-1, 0, 5)

def s_of(model, th, x, idx=None):
    if model == 'M0w':
        return np.full(len(x), th[0])
    if model == 'M1bw':
        lnN, sf = th
        return np.sqrt((osc_shape(x)/(math.sqrt(math.exp(lnN))*LN10))**2 + sf*sf)
    if model == 'M2w':
        b = BIN_FID if idx is None else BIN_FID[idx]
        return np.asarray(th)[b]
    lnsA, sf = th
    A = AREA if idx is None else AREA[idx]
    stt = osc_shape(x)/(np.sqrt(math.exp(lnsA)*A)*LN10)
    return np.sqrt(stt*stt + sf*sf)

S_BOUNDS = {'M0w': [(1e-4, 0.5)], 'M1bw': [(-2, 12), (0.0, 0.5)],
            'M2w': [(1e-4, 0.5)]*6, 'M3w': [(-8, 10), (0.0, 0.5)]}
S_START = {'M0w': [[0.07]], 'M1bw': [[math.log(60.0), 0.05]],
           'M2w': [[0.07]*6], 'M3w': [[math.log(2.0), 0.05]]}

def resid_and_se(la0, f, dd, db, model, sth):
    fac_d = f*np.exp(dd[gidx])
    fac_b = np.exp(db[gidx])
    gN = g_gas + fac_d*g_dsk + fac_b*g_bul
    x = np.sqrt(gN/10**la0)
    r = lgobs - np.log10(gN*nu_be(gN/10**la0))
    s = s_of(model, sth, x)
    return r, sig2 + s*s

def total_obj(la0, f, dd, db, dv, model, sth, pd, pv):
    r, se2 = resid_and_se(la0, f, dd, db, model, sth)
    rr = r - dv[gidx]
    out = np.sum(rr*rr/se2 + np.log(se2))
    out += np.sum(dd*dd)/(pd*pd) + np.sum(dv*dv/(pv*pv))
    out += np.sum(db[HASB]**2)/(pd*pd)
    return out

def fit_w(model, pd=S_ML, pv=None, rounds=4, inject_lg=None):
    global lgobs
    pv = SIGV if pv is None else pv
    lg_data = lgobs if inject_lg is None else inject_lg
    dd = np.zeros(NGal); db = np.zeros(NGal); dv = np.zeros(NGal)
    gb = None
    lg_save = lgobs; lgobs = lg_data
    try:
        for rd in range(rounds):
            def gobj(th):
                la0, f = th[0], th[1]
                if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
                sth = th[2:]
                for v, (lo, hi) in zip(sth, S_BOUNDS[model]):
                    if not (lo <= v <= hi): return 1e12
                return total_obj(la0, f, dd, db, dv, model, sth, pd, pv)
            starts = ([list(gb.x)] if gb is not None else []) + \
                     [[math.log10(A0_FID), 1.0] + s0 for s0 in S_START[model]]
            gbest = None
            for t0 in starts:
                b = minimize(gobj, t0, method='Nelder-Mead',
                             options=dict(maxiter=6000, xatol=1e-6, fatol=1e-6))
                if gbest is None or b.fun < gbest.fun: gbest = b
            gb = gbest
            la0, f, sth = gb.x[0], gb.x[1], gb.x[2:]
            if pd < 1e-3:
                continue
            for sweep in range(2):
                r, se2 = resid_and_se(la0, f, dd, db, model, sth)
                for gi in range(NGal):
                    m = GIDXS[gi]
                    w = 1.0/se2[m]
                    dv[gi] = np.sum(w*r[m])/(np.sum(w) + 1.0/pv[gi]**2)
                for gi in range(NGal):
                    m = GIDXS[gi]
                    def od(d):
                        fac = f*math.exp(d)
                        gN = g_gas[m] + fac*g_dsk[m] + np.exp(db[gi])*g_bul[m]
                        xx = np.sqrt(gN/10**la0)
                        rr = (lgobs[m] - np.log10(gN*nu_be(gN/10**la0))
                              - dv[gi])
                        ss = s_of(model, sth, xx, idx=m)
                        s2 = sig2[m] + ss*ss
                        return np.sum(rr*rr/s2 + np.log(s2)) + d*d/(pd*pd)
                    dd[gi] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                             method='bounded').x
                    if HASB[gi]:
                        def ob(d):
                            fac = f*math.exp(dd[gi])
                            gN = (g_gas[m] + fac*g_dsk[m]
                                  + math.exp(d)*g_bul[m])
                            xx = np.sqrt(gN/10**la0)
                            rr = (lgobs[m] - np.log10(gN*nu_be(gN/10**la0))
                                  - dv[gi])
                            ss = s_of(model, sth, xx, idx=m)
                            s2 = sig2[m] + ss*ss
                            return (np.sum(rr*rr/s2 + np.log(s2))
                                    + d*d/(pd*pd))
                        db[gi] = minimize_scalar(ob, bounds=(-0.7, 0.7),
                                                 method='bounded').x
    finally:
        lgobs = lg_save
    return gb, dd, db, dv

L = [f"STAGE 4W full hierarchy: {kept} galaxies ({int(HASB.sum())} with bulges), "
     f"{len(gobs)} points; priors: M/L 0.1 dex, vertical = measured "
     f"(sigma_v 16/50/84 = {np.percentile(SIGV,[16,50,84]).round(3).tolist()} dex)",
     ""]

res = {}
for model in ('M0w', 'M1bw', 'M2w', 'M3w'):
    res[model] = fit_w(model)
g0, g1, g2, g3 = (res[m][0] for m in ('M0w', 'M1bw', 'M2w', 'M3w'))
L.append(f"M0w  const:        obj = {g0.fun:10.2f}  s0 = {g0.x[2]:.4f}")
L.append(f"M1bw osc+floor:    obj = {g1.fun:10.2f}  N-hat = "
         f"{math.exp(g1.x[2]):8.1f}  floor = {g1.x[3]:.4f}")
L.append(f"M2w  free 6-bin:   obj = {g2.fun:10.2f}  s_b = "
         f"{np.round(g2.x[2:8],4).tolist()}")
L.append(f"M3w  local (area): obj = {g3.fun:10.2f}  sigA = "
         f"{math.exp(g3.x[2]):.3f}/kpc^2  floor = {g3.x[3]:.4f}")
L.append(f"D: M1bw-M0w = {g1.fun-g0.fun:+.2f} | M2w-M0w = {g2.fun-g0.fun:+.2f} "
         f"| M3w-M1bw = {g3.fun-g1.fun:+.2f}")
L.append("4U comparators (disk-M/L only): M1bh-M0h = -22.77; "
         "s_b = [0.1045,0.0778,0.0545,0.0272,0.0646,0.0837] (V-shaped)")
L.append("")
dd1, db1, dv1 = res['M1bw'][1], res['M1bw'][2], res['M1bw'][3]
L.append(f"offsets (M1bw): std(dd) = {np.std(dd1)/LN10:.4f} dex; "
         f"std(db|bulge) = {np.std(db1[HASB])/LN10:.4f}; "
         f"std(dv) = {np.std(dv1):.4f} dex (prior median {np.median(SIGV):.3f})")
medx_g = np.array([np.median(x_fid[GIDXS[i]]) for i in range(NGal)])
L.append(f"corr(|dd|, med x_g) = {np.corrcoef(np.abs(dd1), medx_g)[0,1]:+.3f}; "
         f"corr(|dv|, med x_g) = {np.corrcoef(np.abs(dv1), medx_g)[0,1]:+.3f}")
L.append("")

# G1: priors -> 0
g1a, *_ = fit_w('M0w', pd=1e-4, pv=np.full(NGal, 1e-4), rounds=1)
g1b, *_ = fit_w('M1bw', pd=1e-4, pv=np.full(NGal, 1e-4), rounds=1)
ok1 = abs(g1a.fun-(-8338.12)) < 1.0 and abs(g1b.fun-(-8363.17)) < 1.0
L.append(f"G1 priors->0: M0w {g1a.fun:.2f} (4T -8338.12), M1bw {g1b.fun:.2f} "
         f"(4T -8363.17) -> {'PASS' if ok1 else 'FAIL'}")

# G2 calibrated injection
rng = np.random.default_rng(13)
dd_t = rng.normal(0, S_ML, NGal)
dv_t = rng.normal(0, SIGV)
gN_t = g_gas + 1.0*np.exp(dd_t[gidx])*g_dsk + g_bul
lg_true = np.log10(gN_t*nu_be(gN_t/A0_FID)) + dv_t[gidx]
noise = rng.normal(0, np.sqrt(sig2 + 0.06**2))
_, dd_r, _, _ = fit_w('M0w', rounds=3, inject_lg=lg_true + noise)
lg_pure = np.log10((g_gas+g_dsk+g_bul)
                   * nu_be((g_gas+g_dsk+g_bul)/A0_FID))
_, dd_n, _, _ = fit_w('M0w', rounds=3, inject_lg=lg_pure
                      + rng.normal(0, np.sqrt(sig2 + 0.06**2)))
var_noise = np.var(dd_n)
slope = np.sum(dd_r*dd_t)/np.sum(dd_t*dd_t)
slope_exp = S_ML**2/(S_ML**2 + var_noise)
ok2 = abs(slope - slope_exp) < 0.25
L.append(f"G2 calibrated injection: slope_obs = {slope:.3f}, slope_exp "
         f"(shrinkage) = {slope_exp:.3f}, corr = "
         f"{np.corrcoef(dd_r, dd_t)[0,1]:.3f} -> {'PASS' if ok2 else 'FAIL'}")

out = "\n".join(L)
print(out)
with open('data/stage4w_fullhier.txt', 'w') as f:
    f.write(out+"\n")
