"""
STAGE 5T (O13b): configuration-dependent beta -- where does the galaxies'
beta = 1/2 vote live, and what do the two acceleration arms read alone?

The two-system puzzle (5P): galaxies read beta-hat = 0.45-0.64, binaries
hold beta ~ 0. The binaries sit AT the transition (y ~ 1) inside a
dominant external field; the galaxy discrimination could live anywhere in
x. The reconciliation hypothesis says beta depends on configuration --
if so, the galaxy vote should be structured in depth. Two instruments:

A. The plain-hier beta profile re-run (5P machinery verbatim, params
   CAPTURED per beta) -- regression-gated against the 5P grid.
B. The DECOMPOSITION: at each beta's full-sample best fit, the Delta(-2lnL)
   vs beta=0 is split into fixed y-bins (assigned once, at the beta=0 fit)
   + lensing + priors. Pure evaluation -- zero convergence risk. Shows
   which accelerations carry the beta vote.
C. FREE subset fits: the low arm (y_fid < 1, with the deep lensing points)
   and the high arm (y_fid >= 1) each get their own full profile over
   beta with all nuisances refit. If the arms disagree on beta-hat, the
   dial is configuration-dependent within a single dataset.

Gates: G1 part-A grid vs 5P (|d| < 2); G2 decomposition rows sum to the
profile totals (< 1e-6); subset la0 edge-riding flagged. Plain-hier only
(the vertical channel on subsets would absorb freely; scope stated).
Writes data/stage5t_betasplit.txt.
"""
import glob, math, os, time
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
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 10.0, 1.0, 3.0))
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
GIDXS_ALL = [np.where(gidx == i)[0] for i in range(NGal)]

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))

def nu_beta(y, beta):
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

def m2h(th, nu, dml, pts, lens_on):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx[pts]])
    gN = g_gas[pts] + fac*g_dsk[pts] + g_bul[pts]
    gm_ = gN*nu(gN/a0)
    se2 = sig2[pts] + s_int*s_int
    r = lgobs[pts] - np.log10(gm_)
    out = np.sum(r*r/se2 + np.log(se2))
    if lens_on:
        lg = l_gbar[lmask] + dlt
        rl = l_gobs[lmask] - (lg + np.log10(nu(10**lg/a0)))
        out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(dml*dml)/(S_ML*S_ML)
    return out

def fit_conv(nu, pts, lens_on, th0=None, dml0=None, tol=0.05, max_rounds=15):
    GS = [np.intersect1d(GIDXS_ALL[i], pts) for i in range(NGal)]
    act = [i for i in range(NGal) if len(GS[i])]
    dml = np.zeros(NGal) if dml0 is None else dml0.copy()
    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08, 0.0]]
                  if rd == 0 and th0 is None else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2h(t, nu, dml, pts, lens_on), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int, dlt = best.x
        se2c = s_int*s_int
        for _ in range(3):
            for gi2 in act:
                mm = GS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = lgobs[mm] - np.log10(gN2*nu(gN2/10**la0))
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2h(best.x, nu, dml, pts, lens_on)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2h(t, nu, dml, pts, lens_on), list(best.x),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml

BGRID = [0.0, 0.25, 0.5, 0.75, 1.0]
REF5P = {0.0: -10435.00, 0.25: -10497.63, 0.5: -10519.79,
         0.75: -10520.59, 1.0: -10510.62}
ALL = np.arange(len(gobs))
L = [f"STAGE 5T: configuration-dependent beta -- {kept} galaxies, "
     f"{len(gobs)} points + {int(lmask.sum())} lensing", ""]

# ---------- A: full-sample profile with captured params ----------
L.append("A. full-sample plain-hier profile (5P regression gate):")
FITS = {}
th, dm = None, None
g1ok = True
for b in BGRID:
    nu = lambda y, b=b: nu_beta(y, b)
    t0 = time.time()
    bb, dm = fit_conv(nu, ALL, True, th0=th, dml0=dm)
    th = bb.x
    FITS[b] = (bb.x.copy(), dm.copy(), bb.fun)
    d = bb.fun - REF5P[b]
    ok = abs(d) < 2.0
    g1ok &= ok
    L.append(f"  b={b:4.2f}: {bb.fun:10.2f}  [5P {REF5P[b]:.2f}, "
             f"d={d:+.2f} {'OK' if ok else 'FAIL'}]  "
             f"({(time.time()-t0)/60:.1f} min)")
    print(L[-1], flush=True)
L.append(f"G1 -> {'PASS' if g1ok else 'FAIL'}")
L.append("")

# ---------- B: decomposition into fixed y-bins ----------
th0_, dml0_, _ = FITS[0.0]
a00 = 10**th0_[0]
gN0 = g_gas + th0_[1]*np.exp(dml0_[gidx])*g_dsk + g_bul
y0 = gN0/a00
EDGES = [0, 0.03, 0.1, 0.3, 1.0, 3.0, np.inf]
BINL = [f"y in [{EDGES[i]:g},{EDGES[i+1]:g})" for i in range(len(EDGES)-1)]
BIN_OF = np.digitize(y0, EDGES[1:-1])

def rows(b):
    thb, dmlb, tot = FITS[b]
    la0, f, s_int, dlt = thb
    a0 = 10**la0
    nu = lambda y: nu_beta(y, b)
    gN = g_gas + f*np.exp(dmlb[gidx])*g_dsk + g_bul
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gN*nu(gN/a0))
    c = r*r/se2 + np.log(se2)
    out = [float(np.sum(c[BIN_OF == k])) for k in range(len(BINL))]
    lg = l_gbar[lmask] + dlt
    rl = l_gobs[lmask] - (lg + np.log10(nu(10**lg/a0)))
    out.append(float(np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
                     + (dlt/DELTA_PRIOR)**2))
    out.append(float(np.sum(dmlb*dmlb)/(S_ML*S_ML)))
    return np.array(out), tot

R0, T0 = rows(0.0)
L.append("B. decomposition of Delta(-2lnL) vs beta=0 (fixed beta=0 bins; "
         "negative = that bin prefers this beta):")
hdr = "   bin (N)                " + "".join(f"  b={b:4.2f}" for b in BGRID[1:])
L.append(hdr)
g2ok = True
DROWS = {}
for b in BGRID[1:]:
    Rb, Tb = rows(b)
    DROWS[b] = Rb - R0
    g2 = abs((Rb - R0).sum() - (Tb - T0))
    g2ok &= g2 < 1e-6
labels = [f"{BINL[k]} ({int((BIN_OF == k).sum())})" for k in range(len(BINL))]
labels += [f"lensing+dlt prior ({int(lmask.sum())})", "M/L prior"]
for k, lab in enumerate(labels):
    L.append(f"  {lab:<24}" + "".join(f" {DROWS[b][k]:+7.2f}"
                                      for b in BGRID[1:]))
L.append(f"  {'TOTAL':<24}" + "".join(
    f" {FITS[b][2]-FITS[0.0][2]:+7.2f}" for b in BGRID[1:]))
L.append(f"G2 decomposition identity -> {'PASS' if g2ok else 'FAIL'}")
L.append("")

# ---------- C: free subset fits ----------
gN_fid = g_gas + 1.0*g_dsk + g_bul
y_fid = gN_fid/A0_FID
LOW = np.where(y_fid < 1.0)[0]
HIGH = np.where(y_fid >= 1.0)[0]
L.append(f"C. free subset fits: LOW arm y_fid<1 ({len(LOW)} pts + lensing), "
         f"HIGH arm y_fid>=1 ({len(HIGH)} pts, no lensing)")
for name, pts, lens_on in (("LOW", LOW, True), ("HIGH", HIGH, False)):
    th, dm = None, None
    prof = []
    for b in BGRID:
        nu = lambda y, b=b: nu_beta(y, b)
        t0 = time.time()
        bb, dm = fit_conv(nu, pts, lens_on, th0=th, dml0=dm)
        th = bb.x
        prof.append(bb.fun)
        la0 = bb.x[0]
        edge = " LA0-EDGE" if (la0 < -10.55 or la0 > -9.45) else ""
        L.append(f"  {name} b={b:4.2f}: {bb.fun:10.2f}  la0={la0:+.3f} "
                 f"f={bb.x[1]:.3f} s_int={bb.x[2]:.3f}{edge}  "
                 f"({(time.time()-t0)/60:.1f} min)")
        print(L[-1], flush=True)
    prof = np.array(prof)
    ib = int(np.argmin(prof))
    bhat = BGRID[ib]
    if 0 < ib < len(BGRID)-1:
        x3, y3 = np.array(BGRID[ib-1:ib+2]), prof[ib-1:ib+2]
        cc2, cc1, _ = np.polyfit(x3, y3, 2)
        if cc2 > 0: bhat = -cc1/(2*cc2)
    dz = prof - prof.min()
    L.append(f"  {name}: beta-hat = {bhat:.2f} "
             f"(edge={'no' if 0 < ib < len(BGRID)-1 else 'YES'}); "
             f"Delta vs b=0: {prof[0]-prof.min():+.1f}; grid Delta: "
             + " ".join(f"{v - prof[0]:+.1f}" for v in prof))
    L.append("")

out = "\n".join(L)
print("\n" + out)
with open('data/stage5t_betasplit.txt', 'w') as f:
    f.write(out + "\n")
print("\nSTAGE 5T done -> data/stage5t_betasplit.txt")
