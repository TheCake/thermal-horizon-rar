"""
STAGE 5D (O6, convergence-hardened): the bath matrix under hierarchical M/L,
adaptive coordinate descent.

5C found the ordering FLIP (boot < BE < standard < simple hierarchically)
but its fixed 3-round descent under-converges by ~10 lnL vs 4Z's
warm-chained grid (4Z lam=0: -10423.72, lam=1: -10435.06 vs 5C standard
-10413.70, BE -10423.69). The boot-BE gap (-73) exceeds that slop, but the
claim needs converged fits and truth-calibrated margins.

Upgrades over 5C:
  - adaptive rounds: alternate (global NM refit, 3x delta_d sweeps) until
    the objective moves < 0.05, max 15 rounds; round trajectory printed;
  - CONV gate: converged hier BE must reach <= 4Z lam=1 (-10435.06) + 1.0,
    hier standard <= 4Z lam=0 (-10423.72) + 1.0 (the warm-chained 4Z values
    are the convergence benchmark);
  - truth-calibrated injections: synthetic data at BE truth, boot truth,
    simple truth (same offsets prior + noise) -> converged {boot, BE,
    simple} margins, so the observed gap has a ruler;
  - 50-rep galaxy bootstrap for the boot-BE sign, warm-started from the
    converged offsets, adaptive-lite (tol 0.5, max 6 rounds).
Writes data/stage5d_hierbath_conv.txt.
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
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_boot(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    u = 0.5*(y + np.sqrt(y*y + 4.0*y))
    for _ in range(14):
        uc = np.minimum(u, 45.0)
        eu = np.exp(uc)
        em1 = eu - 1.0
        F = u - y - y/em1
        dF = 1.0 + y*eu/(em1*em1)
        u = np.maximum(u - F/dF, 1e-13)
    nu = u/y
    big = y > 45.0
    if np.any(big):
        yb = np.minimum(y, 700.0)
        nu = np.where(big, 1.0 + 1.0/np.expm1(yb), nu)
    return nu

FAMS = {'BE': nu_be, 'simple': nu_simple, 'boot': nu_boot,
        'standard': nu_standard}

def m2h(th, nu, dml, w_g, lens_obs):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm = gN*nu(gN/a0)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm)
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    rl = lens_obs[lmask] - (lg + np.log10(nu(10**lg/a0)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    return out

def sweep_dml(nu, la0, f, s_int, dml, w_g, n_sweep=3):
    se2c = s_int*s_int
    for _ in range(n_sweep):
        for gi2 in range(NGal):
            if w_g[gi2] == 0: dml[gi2] = 0.0; continue
            mm = GIDXS[gi2]
            def od(dl):
                fc = f*math.exp(dl)
                gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                rr = lgobs[mm] - np.log10(gN2*nu(gN2/10**la0))
                s2 = sig2[mm] + se2c
                return w_g[gi2]*(np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML))
            dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                       method='bounded').x
    return dml

def fit_conv(nu, w_g=None, lens_obs=None, th0=None, dml0=None,
             tol=0.05, max_rounds=15, lg_override=None, trace=None):
    global lgobs
    lg_save = lgobs
    if lg_override is not None: lgobs = lg_override
    if w_g is None: w_g = np.ones(NGal)
    if lens_obs is None: lens_obs = l_gobs
    dml = np.zeros(NGal) if dml0 is None else dml0.copy()
    best = None
    prev = None
    try:
        for rd in range(max_rounds):
            starts = ([list(best.x)] if best is not None else []) + \
                     ([list(th0)] if th0 is not None and best is None else []) + \
                     ([[math.log10(A0_FID), 1.0, 0.08, 0.0]] if rd == 0 else [])
            bb = None
            for t0 in starts:
                b = minimize(lambda t: m2h(t, nu, dml, w_g, lens_obs), t0,
                             method='Nelder-Mead',
                             options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
                if bb is None or b.fun < bb.fun: bb = b
            best = bb
            la0, f, s_int, dlt = best.x
            dml = sweep_dml(nu, la0, f, s_int, dml, w_g)
            cur = m2h(best.x, nu, dml, w_g, lens_obs)
            if trace is not None: trace.append(cur)
            if prev is not None and abs(prev - cur) < tol:
                prev = cur
                break
            prev = cur
    finally:
        lgobs = lg_save
    # final global refit at converged offsets
    b = minimize(lambda t: m2h(t, nu, dml, w_g,
                               lens_obs if lg_override is None else lens_obs),
                 list(best.x), method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if lg_override is not None:
        lgobs = lg_override
        b = minimize(lambda t: m2h(t, nu, dml, w_g, lens_obs), list(best.x),
                     method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
        lgobs = lg_save
    return b, dml

L = [f"STAGE 5D bath matrix, converged: {kept} galaxies, {len(gobs)} points "
     f"+ {int(lmask.sum())} lensing; adaptive descent tol 0.05, max 15 rounds",
     ""]

# ---------------- converged real-data ladder ----------------
res, dms = {}, {}
for name in ('boot', 'BE', 'simple', 'standard'):
    tr = []
    b, dm = fit_conv(FAMS[name], trace=tr)
    res[name] = b; dms[name] = dm
    la0, f, s_int, dlt = b.x
    L.append(f"  {name:>8}: {b.fun:10.2f}  rounds={len(tr)}  a0={10**la0:.3e}"
             f"  f_ML={f:.2f}  s_int={s_int:.3f}  dlt={dlt:+.3f}")
    L.append(f"           trajectory: "
             + " ".join(f"{v:.1f}" for v in tr[:12]))
hb, he, hs, hstd = (res[k].fun for k in ('boot', 'BE', 'simple', 'standard'))
L.append(f"  Delta: boot-BE = {hb-he:+.2f} | boot-simple = {hb-hs:+.2f} | "
         f"BE-simple = {he-hs:+.2f} | boot-standard = {hb-hstd:+.2f}")
L.append(f"  [5C 3-round comparators: -73.28 / -170.56 / -97.28 / -83.26]")
ok_conv = he <= -10435.06 + 1.0 and hstd <= -10423.72 + 1.0
L.append(f"CONV gate vs 4Z warm-chained: BE {he:.2f} (<= -10434.06?), "
         f"standard {hstd:.2f} (<= -10422.72?) -> "
         f"{'PASS' if ok_conv else 'FAIL'}")
L.append("")

# ---------------- truth-calibrated injections ----------------
rng = np.random.default_rng(11)
L.append("truth-calibrated injections (offsets 0.1 dex + noise 0.08, "
         "converged fits, margins = fitted-family minus truth-family):")
for tname in ('BE', 'boot', 'simple'):
    d_t = rng.normal(0, S_ML, NGal)
    gN_t = g_gas + 1.0*np.exp(d_t[gidx])*g_dsk + g_bul
    lg_inj = (np.log10(gN_t*FAMS[tname](gN_t/A0_FID))
              + rng.normal(0, np.sqrt(sig2 + 0.08**2)))
    vals = {}
    for fname in ('boot', 'BE', 'simple'):
        b, _ = fit_conv(FAMS[fname], lg_override=lg_inj, tol=0.1,
                        max_rounds=10, th0=res[fname].x,
                        dml0=dms[fname])
        vals[fname] = b.fun
    margins = {k: vals[k]-vals[tname] for k in vals}
    L.append(f"  truth={tname:>6}: " + " | ".join(
        f"{k} {vals[k]:9.2f} ({margins[k]:+7.2f})" for k in
        ('boot', 'BE', 'simple'))
        + f"  -> truth wins: {'YES' if min(vals, key=vals.get)==tname else 'NO'}")
L.append("  [observed real-data boot-BE for scale: "
         f"{hb-he:+.2f}]")
L.append("")

# ---------------- bootstrap for the flip sign ----------------
rng2 = np.random.default_rng(47)
allg = np.arange(NGal)
d_bb = []
for k in range(50):
    pick = rng2.choice(allg, NGal, replace=True)
    w = np.zeros(NGal)
    for g_ in pick: w[g_] += 1
    lo = l_gobs + rng2.normal(0, np.sqrt(l_sig2))
    vb, _ = fit_conv(nu_boot, w_g=w, lens_obs=lo, th0=res['boot'].x,
                     dml0=dms['boot'], tol=0.5, max_rounds=6)
    ve, _ = fit_conv(nu_be, w_g=w, lens_obs=lo, th0=res['BE'].x,
                     dml0=dms['BE'], tol=0.5, max_rounds=6)
    d_bb.append(vb.fun - ve.fun)
d_bb = np.array(d_bb)
L.append(f"bootstrap (50 reps, converged-warm): boot - BE = {d_bb.mean():+.2f}"
         f" +/- {d_bb.std(ddof=1):.2f} (boot better in "
         f"{int((d_bb<0).sum())}/50)")

out = "\n".join(L)
print(out)
with open('data/stage5d_hierbath_conv.txt', 'w') as f:
    f.write(out+"\n")
print("\nsaved: data/stage5d_hierbath_conv.txt")
