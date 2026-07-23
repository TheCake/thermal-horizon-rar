"""
STAGE 4U: the second moment, hierarchically M/L-marginalized (the decisive
version of 4T) -- plus the mode-locality (annulus-area) axis.

4T found the RAR intrinsic scatter x-dependent, with the oscillator+floor
capturing the monotone decline (N ~ 21) and an x~1 bump flagged as probable
per-galaxy M/L scatter. Here every galaxy gets a free disk-M/L offset
delta_g (lognormal prior, 0.1 dex), profiled jointly with the globals; the
scatter models then compete on what REMAINS:

  M0h  constant s0
  M1bh oscillator + floor:  s^2 = [shape(x)/(sqrt(N) ln10)]^2 + sf^2,
       shape = sqrt(n/(n+1)), n = 1/(e^x - 1)
  M2h  free per x-bin (the 4T sextiles, frozen)
  M3h  local-mode variant: thermal term ~ 1/sqrt(sigA * 2*pi*R*dR) per point
       (modes with surface density sigA per kpc^2 averaged over the annulus).
       Beam smearing scales the OPPOSITE way with annulus size, so the sign
       of the preference is a discriminator.

Questions: does the monotone decline survive M/L marginalization? does the
x~1 bump vanish (M/L identified) or persist? what happens to N-hat? does
the data prefer area-scaled (local) or x-only (global) thermal noise?
Also reported: Var(delta_hat_g) vs prior, and corr(|delta_hat_g|, median x_g)
-- a galaxy-level thermal draw would inflate low-x galaxies' offsets.

Gates:
  G1 prior->0 regression: with sigma_ML = 1e-4 the M0h/M1bh fits reproduce
     Stage 4T's M0/M1b -2lnL within 1.0.
  G2 injection: synthetic per-galaxy offsets (0.1 dex) recovered with
     corr(delta_hat, delta_true) > 0.8.
Writes data/stage4u_mlmarg.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)
SIG_ML = 0.1*LN10          # 0.1 dex lognormal prior on per-galaxy disk M/L

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

g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc = [], [], [], [], [], [], []
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
        gal_id.append(gi); R_kpc.append(R)
g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc))
sig2 = sig*sig
lgobs = np.log10(gobs)

# annulus widths dR per point (within each galaxy, midpoint spacing)
ug = np.unique(gal_id)
NGal = len(ug)
gmap = {g: i for i, g in enumerate(ug)}
gidx = np.array([gmap[g] for g in gal_id])
dR = np.zeros(len(R_kpc))
for gi in range(NGal):
    m = np.where(gidx == gi)[0]
    Rg = R_kpc[m]
    d = np.zeros(len(Rg))
    if len(Rg) == 1:
        d[0] = Rg[0]
    else:
        d[1:-1] = 0.5*(Rg[2:]-Rg[:-2])
        d[0] = Rg[1]-Rg[0]; d[-1] = Rg[-1]-Rg[-2]
    dR[m] = np.maximum(d, 1e-3)
AREA = 2*np.pi*R_kpc*dR          # kpc^2 per annulus

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

def s_model(model, th, x):
    if model == 'M0h':
        return np.full(len(x), th[0])
    if model == 'M1bh':
        lnN, sf = th
        return np.sqrt((osc_shape(x)/(math.sqrt(math.exp(lnN))*LN10))**2 + sf*sf)
    if model == 'M2h':
        return np.asarray(th)[BIN_FID]
    if model == 'M3h':
        lnsA, sf = th
        sth = osc_shape(x)/(np.sqrt(math.exp(lnsA)*AREA)*LN10)
        return np.sqrt(sth*sth + sf*sf)
    raise ValueError(model)

S_BOUNDS = {'M0h': [(1e-4, 0.5)],
            'M1bh': [(-2, 12), (0.0, 0.5)],
            'M2h': [(1e-4, 0.5)]*6,
            'M3h': [(-8, 10), (0.0, 0.5)]}
S_START = {'M0h': [[0.10], [0.06]],
           'M1bh': [[math.log(30.0), 0.05], [math.log(8.0), 0.10]],
           'M2h': [[0.10]*6, [0.06]*6],
           'M3h': [[math.log(3.0), 0.05], [math.log(0.1), 0.10]]}

def m2ll_all(la0, f, dml, model, sth, prior_sig):
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    x = np.sqrt(gN/10**la0)
    gm = gN*nu_be(gN/10**la0)
    s = s_model(model, sth, x)
    se2 = sig2 + s*s
    r = lgobs - np.log10(gm)
    return (np.sum(r*r/se2 + np.log(se2))
            + np.sum(dml*dml)/(prior_sig*prior_sig))

def fit_h(model, prior_sig=SIG_ML, rounds=4, lg0=None):
    dml = np.zeros(NGal)
    best_glob = None
    for rd in range(rounds):
        def gobj(th):
            la0, f = th[0], th[1]
            if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
            sth = th[2:]
            for v, (lo, hi) in zip(sth, S_BOUNDS[model]):
                if not (lo <= v <= hi): return 1e12
            return m2ll_all(la0, f, dml, model, sth, prior_sig)
        starts = ([list(best_glob)] if best_glob is not None else []) + \
                 [[math.log10(A0_FID), 1.0] + s0 for s0 in S_START[model]]
        gb = None
        for t0 in starts:
            b = minimize(gobj, t0, method='Nelder-Mead',
                         options=dict(maxiter=6000, xatol=1e-6, fatol=1e-6))
            if gb is None or b.fun < gb.fun: gb = b
        best_glob = gb.x
        la0, f, sth = gb.x[0], gb.x[1], gb.x[2:]
        if prior_sig < 1e-3:
            continue                      # G1 mode: keep dml pinned at zero
        for gi in range(NGal):
            m = gidx == gi
            def lobj(d):
                fac = f*math.exp(d)
                gN = g_gas[m] + fac*g_dsk[m] + g_bul[m]
                x = np.sqrt(gN/10**la0)
                gm = gN*nu_be(gN/10**la0)
                if model == 'M0h':
                    s = np.full(int(m.sum()), sth[0])
                elif model == 'M1bh':
                    lnN, sf = sth
                    s = np.sqrt((osc_shape(x)/(math.sqrt(math.exp(lnN))*LN10))**2
                                + sf*sf)
                elif model == 'M2h':
                    s = np.asarray(sth)[BIN_FID[m]]
                else:
                    lnsA, sf = sth
                    stt = osc_shape(x)/(np.sqrt(math.exp(lnsA)*AREA[m])*LN10)
                    s = np.sqrt(stt*stt + sf*sf)
                se2 = sig2[m] + s*s
                r = lgobs[m] - np.log10(gm)
                return np.sum(r*r/se2 + np.log(se2)) + d*d/(prior_sig*prior_sig)
            dml[gi] = minimize_scalar(lobj, bounds=(-0.7, 0.7),
                                      method='bounded').x
    return gb, dml

L = [f"STAGE 4U hierarchical M/L second moment: {kept} galaxies, "
     f"{len(gobs)} points; prior sigma_ML = 0.1 dex; family = BE",
     f"annulus area percentiles (kpc^2): "
     f"{np.round(np.percentile(AREA,[16,50,84]),2).tolist()}", ""]

res = {}
for model in ('M0h', 'M1bh', 'M2h', 'M3h'):
    gb, dml = fit_h(model)
    res[model] = (gb, dml)
gb0 = res['M0h'][0]; gb1 = res['M1bh'][0]
gb2 = res['M2h'][0]; gb3 = res['M3h'][0]
L.append(f"M0h  const:        -2lnL+prior = {gb0.fun:10.2f}  s0 = {gb0.x[2]:.4f}")
L.append(f"M1bh osc+floor:    -2lnL+prior = {gb1.fun:10.2f}  N-hat = "
         f"{math.exp(gb1.x[2]):8.1f}  floor = {gb1.x[3]:.4f}")
L.append(f"M2h  free 6-bin:   -2lnL+prior = {gb2.fun:10.2f}  s_b = "
         f"{np.round(gb2.x[2:8],4).tolist()}")
L.append(f"M3h  local (area): -2lnL+prior = {gb3.fun:10.2f}  sigA = "
         f"{math.exp(gb3.x[2]):.3f} modes/kpc^2  floor = {gb3.x[3]:.4f}")
L.append(f"D: M1bh-M0h = {gb1.fun-gb0.fun:+.2f} | M2h-M0h = {gb2.fun-gb0.fun:+.2f} "
         f"| M3h-M1bh = {gb3.fun-gb1.fun:+.2f} (local-vs-global thermal)")
L.append("")
L.append("4T (unmarginalized) comparators: M0 -8338.12, M1b -8363.17 (N=21.5, "
         "floor 0.1014), M2 [0.1438,0.1268,0.1146,0.1112,0.1262,0.1073]")
L.append("")

# delta_g diagnostics (from the M1bh fit)
dml1 = res['M1bh'][1]
medx_g = np.array([np.median(x_fid[gidx == gi]) for gi in range(NGal)])
sd_d = np.std(dml1)/LN10
cc = np.corrcoef(np.abs(dml1), medx_g)[0, 1]
L.append(f"delta_g diagnostics (M1bh): std(delta_hat) = {sd_d:.4f} dex "
         f"(prior 0.1); corr(|delta_hat|, median x_g) = {cc:+.3f}")
L.append("")

# G1: prior -> 0 regression vs 4T
g1a, _ = fit_h('M0h', prior_sig=1e-4, rounds=1)
g1b, _ = fit_h('M1bh', prior_sig=1e-4, rounds=1)
ok1 = abs(g1a.fun - (-8338.12)) < 1.0 and abs(g1b.fun - (-8363.17)) < 1.0
L.append(f"G1 prior->0: M0h {g1a.fun:.2f} (4T -8338.12), M1bh {g1b.fun:.2f} "
         f"(4T -8363.17) -> {'PASS' if ok1 else 'FAIL'}")

# G2: injection of known per-galaxy offsets
rng = np.random.default_rng(7)
d_true = rng.normal(0, SIG_ML, NGal)
fac_t = 1.0*np.exp(d_true[gidx])
gN_t = g_gas + fac_t*g_dsk + g_bul
lg_true = np.log10(gN_t*nu_be(gN_t/A0_FID))
lgobs_save = lgobs.copy()
lgobs = lg_true + rng.normal(0, np.sqrt(sig2 + 0.08**2))
_, d_rec = fit_h('M0h', rounds=3)
cinj = np.corrcoef(d_rec, d_true)[0, 1]
lgobs = lgobs_save
L.append(f"G2 injection: corr(delta_hat, delta_true) = {cinj:.3f} -> "
         f"{'PASS' if cinj > 0.8 else 'FAIL'}")

out = "\n".join(L)
print(out)
with open('data/stage4u_mlmarg.txt', 'w') as f:
    f.write(out+"\n")
