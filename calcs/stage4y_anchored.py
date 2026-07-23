"""
STAGE 4Y (O1b): the anchored-distance second moment. The 4W identifiability
boundary broke on per-galaxy vertical freedom; here the contest reruns on the
subsample whose distances are externally anchored (SPARC e_D/D <= 0.10 --
TRGB/Cepheid/SN-grade), where the vertical priors pin small and the
thermal-vs-vertical degeneracy should break.

Machinery = 4W (disk M/L + bulge M/L + measured-prior vertical offsets),
restricted to the anchored subsample; x-bins FROZEN at the full-sample
sextile edges so bin-by-bin comparisons carry over.

Gates: G1 priors->0 on the subsample is self-consistent (M2 free >= others);
G2 calibrated injection ON the subsample -- the claim that anchoring restores
identifiability is itself tested (slope_exp should rise vs 4W's 0.93-vs-0.23
failure).
Writes data/stage4y_anchored.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)
S_ML = 0.1*LN10
EDD_MAX = 0.10

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

rows = []
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept_all = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 10.0, 1.0, 3.0))
    if inc < 30 or q > 2: continue
    kept_all += 1
    anchored = (eD/max(D, 1e-3)) <= EDD_MAX
    irad = math.radians(inc)
    sv = math.hypot((eD/max(D, 1e-3))/LN10,
                    2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10)
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC; gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        rows.append((gi, anchored, max(sv, 0.01), gg, gd, gb,
                     Vo*Vo/R*KPC, 2*eV/Vo/math.log(10)))
rows = np.array(rows, dtype=object)
gid_all = np.array([r[0] for r in rows])
anch_pt = np.array([r[1] for r in rows], bool)
g_gas_a = np.array([r[3] for r in rows]); g_dsk_a = np.array([r[4] for r in rows])
g_bul_a = np.array([r[5] for r in rows]); gobs_a = np.array([r[6] for r in rows])
sig_a = np.array([r[7] for r in rows]); sv_pt = np.array([r[2] for r in rows])

# full-sample frozen bins
x_all = np.sqrt((g_gas_a + g_dsk_a + g_bul_a)/A0_FID)
QE = np.quantile(x_all, np.linspace(0, 1, 7)); QE[0], QE[-1] = 0.0, np.inf

# restrict to the anchored subsample
m = anch_pt
g_gas, g_dsk, g_bul = g_gas_a[m], g_dsk_a[m], g_bul_a[m]
gobs, sig = gobs_a[m], sig_a[m]
sig2 = sig*sig
lgobs = np.log10(gobs)
gid_s = gid_all[m]
ug = np.unique(gid_s)
NGal = len(ug)
gmap = {g: i for i, g in enumerate(ug)}
gidx = np.array([gmap[g] for g in gid_s])
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]
SIGV = np.array([sv_pt[m][GIDXS[i][0]] for i in range(NGal)])
HASB = np.array([bool(np.any(g_bul[GIDXS[i]] > 0)) for i in range(NGal)])
BIN_FID = np.clip(np.searchsorted(
    QE, np.sqrt((g_gas+g_dsk+g_bul)/A0_FID), side='right')-1, 0, 5)

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def osc_shape(x):
    xc = np.minimum(np.clip(x, 1e-9, None), 80.0)
    n = 1.0/np.expm1(xc)
    return np.sqrt(n/(1.0+n))

def s_of(model, th, x, idx=None):
    if model == 'M0':
        return np.full(len(x), th[0])
    if model == 'M1b':
        lnN, sf = th
        return np.sqrt((osc_shape(x)/(math.sqrt(math.exp(lnN))*LN10))**2 + sf*sf)
    b = BIN_FID if idx is None else BIN_FID[idx]
    return np.asarray(th)[b]

S_BOUNDS = {'M0': [(1e-4, 0.5)], 'M1b': [(-2, 12), (0.0, 0.5)],
            'M2': [(1e-4, 0.5)]*6}
S_START = {'M0': [[0.08]], 'M1b': [[math.log(40.0), 0.05]], 'M2': [[0.08]*6]}

def fit_y(model, pd=S_ML, pv=None, rounds=4, inject_lg=None):
    global lgobs
    pv = SIGV if pv is None else pv
    lg_save = lgobs
    if inject_lg is not None: lgobs = inject_lg
    dd = np.zeros(NGal); db = np.zeros(NGal); dv = np.zeros(NGal)
    gb = None
    try:
        for rd in range(rounds):
            def gobj(th):
                la0, f = th[0], th[1]
                if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
                sth = th[2:]
                for v, (lo, hi) in zip(sth, S_BOUNDS[model]):
                    if not (lo <= v <= hi): return 1e12
                fac = f*np.exp(dd[gidx]); fb = np.exp(db[gidx])
                gN = g_gas + fac*g_dsk + fb*g_bul
                x = np.sqrt(gN/10**la0)
                r = lgobs - np.log10(gN*nu_be(gN/10**la0)) - dv[gidx]
                s = s_of(model, sth, x)
                se2 = sig2 + s*s
                return (np.sum(r*r/se2 + np.log(se2))
                        + np.sum(dd*dd)/(pd*pd) + np.sum(dv*dv/(pv*pv))
                        + np.sum(db[HASB]**2)/(pd*pd))
            starts = ([list(gb.x)] if gb is not None else []) + \
                     [[math.log10(A0_FID), 1.0] + s0 for s0 in S_START[model]]
            gbest = None
            for t0 in starts:
                b = minimize(gobj, t0, method='Nelder-Mead',
                             options=dict(maxiter=6000, xatol=1e-6, fatol=1e-6))
                if gbest is None or b.fun < gbest.fun: gbest = b
            gb = gbest
            la0, f, sth = gb.x[0], gb.x[1], gb.x[2:]
            if pd < 1e-3: continue
            for sweep in range(2):
                fac = f*np.exp(dd[gidx]); fb = np.exp(db[gidx])
                gN = g_gas + fac*g_dsk + fb*g_bul
                x = np.sqrt(gN/10**la0)
                r0 = lgobs - np.log10(gN*nu_be(gN/10**la0))
                s = s_of(model, sth, x)
                se2 = sig2 + s*s
                for gi2 in range(NGal):
                    mm = GIDXS[gi2]
                    w = 1.0/se2[mm]
                    dv[gi2] = np.sum(w*r0[mm])/(np.sum(w) + 1.0/pv[gi2]**2)
                for gi2 in range(NGal):
                    mm = GIDXS[gi2]
                    def od(dl):
                        fc = f*math.exp(dl)
                        gN2 = (g_gas[mm] + fc*g_dsk[mm]
                               + math.exp(db[gi2])*g_bul[mm])
                        xx = np.sqrt(gN2/10**la0)
                        rr = (lgobs[mm] - np.log10(gN2*nu_be(gN2/10**la0))
                              - dv[gi2])
                        ss = s_of(model, sth, xx, idx=mm)
                        s2 = sig2[mm] + ss*ss
                        return np.sum(rr*rr/s2 + np.log(s2)) + dl*dl/(pd*pd)
                    dd[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                              method='bounded').x
                    if HASB[gi2]:
                        def ob(dl):
                            fc = f*math.exp(dd[gi2])
                            gN2 = (g_gas[mm] + fc*g_dsk[mm]
                                   + math.exp(dl)*g_bul[mm])
                            xx = np.sqrt(gN2/10**la0)
                            rr = (lgobs[mm]
                                  - np.log10(gN2*nu_be(gN2/10**la0)) - dv[gi2])
                            ss = s_of(model, sth, xx, idx=mm)
                            s2 = sig2[mm] + ss*ss
                            return (np.sum(rr*rr/s2 + np.log(s2))
                                    + dl*dl/(pd*pd))
                        db[gi2] = minimize_scalar(ob, bounds=(-0.7, 0.7),
                                                  method='bounded').x
    finally:
        lgobs = lg_save
    return gb, dd, db, dv

L = [f"STAGE 4Y anchored second moment: e_D/D <= {EDD_MAX} subsample = "
     f"{NGal}/{kept_all} galaxies, {len(gobs)} points "
     f"({int(HASB.sum())} with bulges)",
     f"sigma_v 16/50/84 = {np.percentile(SIGV,[16,50,84]).round(3).tolist()} dex "
     f"(4W full sample: [0.046, 0.097, 0.139])", ""]

res = {}
for model in ('M0', 'M1b', 'M2'):
    res[model] = fit_y(model)
g0, g1, g2 = (res[m2][0] for m2 in ('M0', 'M1b', 'M2'))
L.append(f"M0  const:      obj = {g0.fun:9.2f}  s0 = {g0.x[2]:.4f}")
L.append(f"M1b osc+floor:  obj = {g1.fun:9.2f}  N-hat = {math.exp(g1.x[2]):8.1f}"
         f"  floor = {g1.x[3]:.4f}")
L.append(f"M2  free 6-bin: obj = {g2.fun:9.2f}  s_b = "
         f"{np.round(g2.x[2:8],4).tolist()}")
L.append(f"D: M1b-M0 = {g1.fun-g0.fun:+.2f} | M2-M0 = {g2.fun-g0.fun:+.2f}")
L.append("(4W full-sample comparators: M1bw-M0w = +24.36; bump bin4 0.054)")
L.append("")

# G2 calibrated injection on the subsample
rng = np.random.default_rng(17)
dd_t = rng.normal(0, S_ML, NGal)
dv_t = rng.normal(0, SIGV)
gN_t = g_gas + 1.0*np.exp(dd_t[gidx])*g_dsk + g_bul
lg_true = np.log10(gN_t*nu_be(gN_t/A0_FID)) + dv_t[gidx]
noise = rng.normal(0, np.sqrt(sig2 + 0.06**2))
_, dd_r, _, _ = fit_y('M0', rounds=3, inject_lg=lg_true + noise)
lg_pure = np.log10((g_gas+g_dsk+g_bul)*nu_be((g_gas+g_dsk+g_bul)/A0_FID))
_, dd_n, _, _ = fit_y('M0', rounds=3,
                      inject_lg=lg_pure + rng.normal(0, np.sqrt(sig2+0.06**2)))
var_noise = np.var(dd_n)
slope = np.sum(dd_r*dd_t)/np.sum(dd_t*dd_t)
slope_exp = S_ML**2/(S_ML**2 + var_noise)
ok2 = abs(slope - slope_exp) < 0.25
L.append(f"G2 calibrated injection (anchored subsample): slope_obs = {slope:.3f}"
         f", slope_exp = {slope_exp:.3f}, corr = "
         f"{np.corrcoef(dd_r, dd_t)[0,1]:.3f} -> {'PASS' if ok2 else 'FAIL'}")
L.append("(4W full-sample gate for contrast: slope 0.225 vs expected 0.931 FAIL)")
ok_nest = g2.fun <= min(g0.fun, g1.fun) + 1e-6
L.append(f"G1 nesting: {'PASS' if ok_nest else 'FAIL'}")

out = "\n".join(L)
print(out)
with open('data/stage4y_anchored.txt', 'w') as f:
    f.write(out+"\n")
