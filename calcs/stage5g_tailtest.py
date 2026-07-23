"""
STAGE 5G (O8): tail or coefficient? Disambiguating the hierarchical flip.

5D found the hierarchical bath ladder boot < BE < standard < simple --
exactly monotone in Newtonian-tail sharpness (e^-y, e^-sqrt(y), ~y^-2,
~y^-1). 5F found the binaries REJECT boot's transition (nu(1)=1.35 too
weak; alpha edge-rides) while holding the 1/2-branch. Hypothesis: the
galaxies' boot preference is a SCREENING-TAIL vote, not a deep-c1=1/4
vote. Probe: the sec.-3 screening family nu_p = (1-e^{-y^p})^{-1/(2p)},
which IS the occupation law at p=1/2 (same deep c1=1/2 there) and
sharpens the Newtonian tail as p grows with only a mild transition change
(nu_p(1): 1.582 at p=.5 -> 1.487 at p=.578 -> 1.376 at p=.75; boot 1.350;
NLO exponent drifts off pure x as 2p, disclosed). If nu_p at p>1/2 climbs
toward boot's hierarchical -10510.6, the flip is a tail story and the
binaries' 1/2 transition survives everywhere; if it stays at BE's
-10435.0, the flip is genuinely the deep coefficient.

Converged hier fits (5D machinery verbatim: adaptive descent, tol 0.05,
max 15 rounds, per-galaxy disk-M/L 0.1-dex prior, joint SPARC+lensing).
Grid p in {0.5, 0.578, 0.65, 0.75, 0.9}, warm-chained.
Gate G1: p=0.5 reproduces converged hier BE (-10435.00) within 0.5.
Writes data/stage5g_tailtest.txt.
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

def nu_p(y, p):
    yc = np.clip(y, 1e-14, None)
    ex = np.exp(-np.minimum(yc**p, 60.0))
    return (1.0-ex)**(-1.0/(2.0*p))

def m2h(th, nu, dml, w_g):
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
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    return out

def fit_conv(nu, th0=None, dml0=None, tol=0.05, max_rounds=15):
    w_g = np.ones(NGal)
    dml = np.zeros(NGal) if dml0 is None else dml0.copy()
    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08, 0.0]] if rd == 0 else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2h(t, nu, dml, w_g), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
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
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2h(best.x, nu, dml, w_g)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2h(t, nu, dml, w_g), list(best.x),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    return (b if b.fun < best.fun else best), dml

L = [f"STAGE 5G tail-vs-coefficient: {kept} galaxies, {len(gobs)} points + "
     f"{int(lmask.sum())} lensing; nu_p family, hierarchical converged",
     f"comparators (5D converged): boot -10510.60, BE -10435.00, "
     f"standard -10424.49, simple -10336.28", ""]

PGRID = [0.5, 0.578, 0.65, 0.75, 0.9]
th, dm = None, None
vals = []
for p in PGRID:
    nu = lambda y, p=p: nu_p(y, p)
    b, dm = fit_conv(nu, th0=th, dml0=dm)
    th = b.x
    vals.append(b.fun)
    la0, f, s_int, dlt = b.x
    n1 = float(nu_p(np.array([1.0]), p)[0])
    L.append(f"  p={p:5.3f}: {b.fun:10.2f}  nu(1)={n1:.3f}  a0={10**la0:.3e}"
             f"  f_ML={f:.2f}  s_int={s_int:.3f}  dlt={dlt:+.3f}")
vals = np.array(vals)
ok1 = abs(vals[0] - (-10435.00)) < 0.5
L.append("")
L.append(f"G1 p=0.5 == hier BE: {vals[0]:.2f} (5D -10435.00) -> "
         f"{'PASS' if ok1 else 'FAIL'}")
ib = int(np.argmin(vals))
L.append(f"minimum at p={PGRID[ib]}; gain vs BE = {vals[ib]-vals[0]:+.2f}; "
         f"boot's gain was -75.59")
frac = (vals[0]-vals[ib])/75.59
L.append(f"fraction of the boot flip reproduced by the tail dial alone: "
         f"{100*frac:.0f}%")

out = "\n".join(L)
print(out)
with open('data/stage5g_tailtest.txt', 'w') as f:
    f.write(out+"\n")
print("\nsaved: data/stage5g_tailtest.txt")
