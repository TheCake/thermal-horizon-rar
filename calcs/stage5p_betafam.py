"""
STAGE 5P (O5 theory instrument): the mixing family -- "how self-consistent
is the bath frequency?" as a measurable dial.

The three bath cells tested so far are members of a one-parameter family:
    omega = omega_source^(1-beta) * omega_total^beta
    =>  nu = 1 + n_BE(u),  u = y^((1+beta)/2) * nu^beta
beta = 0 is the occupation law (source-driven), beta = 1 the quantum
bootstrap (self-consistent), beta = 1/2 the geometric mean. Derived exact
coefficients (verified numerically below):
    c1(beta) = 1/(2(1+beta))
    c2(beta) = 1/(12(1+beta)) + beta/(8(1+beta)^2)
    tail exponent p_tail = (1+beta)/2
so the family carries the parameter-free RELATION c1 * p_tail = 1/4: the
zero-point coefficient and the screening sharpness are locked together by
the mixing exponent. (beta=0: 1/2,1/12; beta=1/2: 1/3,1/12 -- the Bernoulli
1/12 recurs; beta=1: 1/4,7/96.) Under deep-MOND scale invariance both
frequencies scale identically, so every beta is scale-covariant; beta = 1/2
is the unique symmetric (source <-> response exchange) member.

This script: G0 series gates (c1, c2 vs the formulas on a beta grid) +
member regressions (beta = 0/0.5/1 reproduce nu_be/nu_gm/nu_boot to 1e-9);
then the hierarchical galaxy profile over beta -- plain (5D machinery) and
with the measured vertical channel (5M machinery) -- giving beta-hat under
both treatments with D1 intervals.
Writes data/stage5p_betafam.txt.
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
def nu_gm_ref(y):
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
def nu_boot_ref(y):
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

def nu_beta(y, beta):
    """nu = 1 + n_BE(y^((1+beta)/2) * nu^beta), Newton on nu (F monotone)."""
    y = np.clip(np.asarray(y, float), 1e-14, None)
    A = y**(0.5*(1.0+beta))
    nu = nu_simple(y)
    for _ in range(40):
        u = np.minimum(A*nu**beta, 60.0)
        eu = np.exp(u)
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        F = nu - 1.0 - n
        dF = 1.0 + (eu/(em1*em1))*u*beta/nu
        nu = np.maximum(nu - F/dF, 1.0 + 1e-15)
    return nu

def c1c2(beta):
    return (1.0/(2*(1+beta)),
            1.0/(12*(1+beta)) + beta/(8*(1+beta)**2))

L = [f"STAGE 5P mixing family: {kept} galaxies, {len(gobs)} points + "
     f"{int(lmask.sum())} lensing", ""]

# G0: series + member regressions
xg = np.linspace(0.005, 0.05, 40)
yg = xg*xg
ok0 = True
L.append("G0 series gates (numeric vs derived c1 = 1/(2(1+b)), "
         "c2 = 1/(12(1+b)) + b/(8(1+b)^2)):")
for b in (0.0, 0.25, 0.5, 0.75, 1.0):
    h = nu_beta(yg, b)*xg - 1.0
    cf = np.polyfit(xg, h/xg, 1)
    c1p, c2p = c1c2(b)
    ok = abs(cf[1]-c1p) < 2e-3 and abs(cf[0]-c2p) < 4e-3
    ok0 &= ok
    L.append(f"  b={b:4.2f}: c1 = {cf[1]:.4f} ({c1p:.4f}), c2 = {cf[0]:.4f}"
             f" ({c2p:.4f}) {'OK' if ok else 'FAIL'}")
yv = np.logspace(-6, 1.5, 200)
r0 = np.max(np.abs(nu_beta(yv, 0.0)/nu_be(yv) - 1))
r5 = np.max(np.abs(nu_beta(yv, 0.5)/nu_gm_ref(yv) - 1))
r1 = np.max(np.abs(nu_beta(yv, 1.0)/nu_boot_ref(yv) - 1))
okm = max(r0, r5, r1) < 1e-9
ok0 &= okm
L.append(f"member regressions: |b=0 vs BE| {r0:.1e}, |b=.5 vs gm| {r5:.1e},"
         f" |b=1 vs boot| {r1:.1e} -> {'OK' if okm else 'FAIL'}")
L.append(f"family relation: c1*p_tail = 1/4 exactly (p_tail = (1+b)/2)")
L.append(f"G0 -> {'PASS' if ok0 else 'FAIL'}")
L.append("")

# ---------------- hierarchical profiles over beta ----------------
def m2hv(th, nu, dml, dv, use_v):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*nu(gN/a0)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm_) - (dv[gidx] if use_v else 0.0)
    out = np.sum(r*r/se2 + np.log(se2))
    lg = l_gbar[lmask] + dlt
    rl = l_gobs[lmask] - (lg + np.log10(nu(10**lg/a0)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(dml*dml)/(S_ML*S_ML)
    if use_v: out += np.sum(dv*dv/(SIGV*SIGV))
    return out

def fit_conv(nu, use_v, th0=None, dml0=None, dv0=None, tol=0.05,
             max_rounds=15):
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
            b = minimize(lambda t: m2hv(t, nu, dml, dv, use_v), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int, dlt = best.x
        se2c = s_int*s_int
        for _ in range(3):
            if use_v:
                fac = f*np.exp(dml[gidx])
                gN = g_gas + fac*g_dsk + g_bul
                r0_ = lgobs - np.log10(gN*nu(gN/10**la0))
                for gi2 in range(NGal):
                    mm = GIDXS[gi2]
                    w = 1.0/(sig2[mm] + se2c)
                    dv[gi2] = np.sum(w*r0_[mm])/(np.sum(w) + 1.0/SIGV[gi2]**2)
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = (lgobs[mm] - np.log10(gN2*nu(gN2/10**la0))
                          - (dv[gi2] if use_v else 0.0))
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2hv(best.x, nu, dml, dv, use_v)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2hv(t, nu, dml, dv, use_v), list(best.x),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml, dv

BGRID = [0.0, 0.25, 0.5, 0.75, 1.0]
REF = {(False, 0.0): -10435.00, (False, 0.5): -10519.78,
       (False, 1.0): -10510.60,
       (True, 0.0): -12152.49, (True, 0.5): -12193.04,
       (True, 1.0): -12161.04}
for use_v, tag in ((False, "plain hier (5D machinery)"),
                   (True, "with vertical channel (5M machinery)")):
    L.append(f"beta profile, {tag}:")
    th, dm, dvv = None, None, None
    prof = []
    for b in BGRID:
        nu = lambda y, b=b: nu_beta(y, b)
        bb, dm, dvv = fit_conv(nu, use_v, th0=th, dml0=dm, dv0=dvv)
        th = bb.x
        prof.append(bb.fun)
        ref = REF.get((use_v, b))
        L.append(f"  b={b:4.2f}: {bb.fun:10.2f}"
                 + (f"  [ref {ref:.2f} {'OK' if abs(bb.fun-ref) < 2.0 else 'DIFF'}]"
                    if ref else ""))
    prof = np.array(prof)
    ib = int(np.argmin(prof))
    bhat = BGRID[ib]
    if 0 < ib < len(BGRID)-1:
        x3, y3 = np.array(BGRID[ib-1:ib+2]), prof[ib-1:ib+2]
        c2_, c1_, _ = np.polyfit(x3, y3, 2)
        if c2_ > 0: bhat = -c1_/(2*c2_)
    dz = prof - prof.min()
    L.append(f"  beta-hat = {bhat:.2f} (edge={'no' if 0 < ib < len(BGRID)-1 else 'YES'}); "
             f"D vs b=0: {dz[0]:+.1f}; implied c1 = {c1c2(bhat)[0]:.3f}, "
             f"p_tail = {(1+bhat)/2:.3f}")
    L.append("")

out = "\n".join(L)
print(out)
with open('data/stage5p_betafam.txt', 'w') as f:
    f.write(out+"\n")
print("\nsaved: data/stage5p_betafam.txt")
