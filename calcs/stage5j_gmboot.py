"""
STAGE 5J (O5 galaxy side): the geometric-mean bootstrap on SPARC+lensing.

The O5 candidate nu_gm: nu = 1 + n_BE(y^(3/4) sqrt(nu)) -- mode frequency
omega = sqrt(omega_source * omega_total), the geometric mean of the
source-driven and self-consistent prescriptions; equivalently the occupation
argument is x^2 = [T_U(g_N)/T_dS]*[T_U(g_obs)/T_dS]. Parameter-free.
Derived series: nu*sqrt(y) = 1 + x/3 + x^2/12 (c1 = 1/3 -- inside both
measured bands; c2 = 1/12 -- the occupation law's Bernoulli coefficient
survives); transition nu(1) = 1.433 (binary-acceptable grade per 5H);
Newtonian tail e^{-y^(3/4)} (p-equivalent 3/4 -- at the top of the 5G hier
band 0.65-0.75).

This script: (a) gates the implementation (series c1/c2, solver residual,
deep/Newtonian limits); (b) runs the converged hierarchical fit (5D
machinery verbatim; comparators boot -10510.60, p065 -10491.38, BE
-10435.00); (c) runs the FLAT-M/L joint fits (prior->0, 4E objective;
BE regression -8397.72) for nu_gm and nu_p(0.65) to supply the a0 ladder's
galaxy leg under the new functions (O10b).
Writes data/stage5j_gmboot.txt.
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
def nu_p065(y):
    yc = np.clip(np.asarray(y, float), 1e-14, None)
    ex = np.exp(-np.minimum(yc**0.65, 60.0))
    return (1.0-ex)**(-1.0/1.3)
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

# ---------------- gates ----------------
xg = np.linspace(0.005, 0.05, 40)
yg = xg*xg
h = nu_gm(yg)*xg - 1.0
cfit = np.polyfit(xg, h/xg, 1)
c1n, c2n = cfit[1], cfit[0]
yv = np.logspace(-8, 2, 400)
w_ = np.sqrt(nu_gm(yv))
u_ = np.minimum(yv**0.75*w_, 60.0)
n_ = np.where(u_ < 60, 1.0/np.expm1(u_), 0.0)
res = np.max(np.abs(w_*w_ - 1.0 - n_)/(1.0 + n_))   # relative: deep n ~ 1e8
lo_lim = float(nu_gm(np.array([1e-10]))[0])*1e-5
hi_lim = float(nu_gm(np.array([100.0]))[0]) - 1.0
n1 = float(nu_gm(np.array([1.0]))[0])
g0ok = (abs(c1n-1/3) < 2e-3 and abs(c2n-1/12) < 4e-3 and res < 1e-12
        and abs(lo_lim-1) < 1e-3 and hi_lim < 1e-12)
L = [f"STAGE 5J geometric-mean bootstrap: {kept} galaxies, {len(gobs)} "
     f"points + {int(lmask.sum())} lensing",
     f"G0 nu_gm: c1 = {c1n:.4f} (1/3), c2 = {c2n:.4f} (1/12 = 0.0833); "
     f"max|H| = {res:.1e}; nu*sqrt(y)@1e-10 = {lo_lim:.5f}; nu(100)-1 = "
     f"{hi_lim:.1e}; nu(1) = {n1:.3f} -> {'PASS' if g0ok else 'FAIL'}", ""]

def m2h(th, nu, dml, w_g, pd):
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
    rl = l_gobs[lmask] - (lg + np.log10(nu(10**lg/a0)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(pd*pd)
    return out

def fit_conv(nu, pd=S_ML, tol=0.05, max_rounds=15, flat=False):
    w_g = np.ones(NGal)
    dml = np.zeros(NGal)
    best = None
    prev = None
    rounds = 1 if flat else max_rounds
    pdu = 1e-4 if flat else pd
    for rd in range(rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08, 0.0],
                   [math.log10(A0_FID)+0.1, 0.8, 0.12, -0.1]] if rd == 0
                  else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2h(t, nu, dml, w_g, pdu), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        if flat: break
        la0, f, s_int, dlt = best.x
        se2c = s_int*s_int
        for _ in range(3):
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = lgobs[mm] - np.log10(gN2*nu(gN2/10**la0))
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(pdu*pdu)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2h(best.x, nu, dml, w_g, pdu)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    if not flat:
        b = minimize(lambda t: m2h(t, nu, dml, w_g, pdu), list(best.x),
                     method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
        if b.fun < best.fun: best = b
    return best

# ---------------- flat fits (a0 ladder legs + regression) ----------------
L.append("FLAT-M/L joint fits (4E objective; a0 ladder legs):")
fb = fit_conv(nu_be, flat=True)
okf = abs(fb.fun - (-8397.72)) < 1.0
L.append(f"  BE regression: {fb.fun:.2f} (4E -8397.72) -> "
         f"{'PASS' if okf else 'FAIL'}")
for name, nu in (("p065", nu_p065), ("gm", nu_gm)):
    b = fit_conv(nu, flat=True)
    la0, f, s_int, dlt = b.x
    L.append(f"  {name:>5}: -2lnL = {b.fun:10.2f}  a0 = {10**la0:.3e}  "
             f"f_ML = {f:.2f}  s_int = {s_int:.3f}  dlt = {dlt:+.3f}")
L.append("")

# ---------------- converged hierarchical fit ----------------
L.append("HIERARCHICAL converged (5D machinery; comparators: boot -10510.60,"
         " p065 -10491.38, BE -10435.00, standard -10424.49):")
bh = fit_conv(nu_gm)
la0, f, s_int, dlt = bh.x
L.append(f"  gm hier: {bh.fun:10.2f}  a0 = {10**la0:.3e}  f_ML = {f:.2f}  "
         f"s_int = {s_int:.3f}  dlt = {dlt:+.3f}")
L.append(f"  gm - boot = {bh.fun-(-10510.60):+.2f} | gm - p065 = "
         f"{bh.fun-(-10491.38):+.2f} | gm - BE = {bh.fun-(-10435.00):+.2f}")

out = "\n".join(L)
print(out)
with open('data/stage5j_gmboot.txt', 'w') as f_:
    f_.write(out+"\n")
print("\nsaved: data/stage5j_gmboot.txt")
