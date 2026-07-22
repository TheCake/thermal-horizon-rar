"""
STAGE 4A: the NLO kill test of the Bose-Einstein identity (TODO #9).
Deep-MOND expansion of nu(y) = 1/(1-exp(-sqrt(y))):
    g = sqrt(g_N a0) * (1 + c1*x + c2*x^2 + ...),  x = sqrt(g_N/a0)
Parameter-free predictions:
    BE:        c1 = 1/2,  c2 = 1/12 = 0.083
    simple-nu: c1 = 1/2,  c2 = 1/8  = 0.125   (same branch at NLO)
    standard:  c1 = 0,    c2 = 1/4              (different branch)
A measured c1 far from 1/2 kills the BE identity; c1 ~ 0 favors the standard
branch. Fit on SPARC low-acceleration points (window scan y < ymax for
truncation stability), free (a0, c1, c2, f_ML) with f_ML a global disk-M/L
scale; errors via galaxy-level bootstrap. Writes data/stage4a_nlo.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10

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
print(f"galaxies {kept}, points {len(gobs)}")

L = [f"STAGE 4A NLO test v2 (two-step asymptotic matching): {kept} galaxies, "
     f"{len(gobs)} points",
     "step A: a0 from ultra-deep points (y<Y_DEEP), corrections negligible",
     "step B: a0 FIXED, linear WLS for (c1,c2) on Y_DEEP<y<ymax",
     "predictions: BE c1=0.5 c2=0.083 | simple c1=0.5 c2=0.125 | "
     "standard c1=0 c2=0.25"]

Y_DEEP = 0.05
F_GRID = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])

def two_step(ymax, gal_weights, fml):
    """returns (a0, c1, c2, chi2) for given galaxy weights and M/L factor"""
    gN = g_gas + fml*g_dsk + g_bul
    y_fid = gN/A0_FID
    w_pt = gal_weights[gal_id]
    # step A: deep points -> a0 (slope-1/2 intercept), correct for the small
    # residual NLO term iteratively using c1=0 first pass
    deep = (y_fid < Y_DEEP) & (w_pt > 0)
    if deep.sum() < 30: return None
    lw = w_pt[deep]/sig[deep]**2
    la0 = np.sum(lw*(2*np.log10(gobs[deep]) - np.log10(gN[deep]))) / np.sum(lw)
    a0 = 10**la0
    # step B: intermediate window, linear WLS for (c1, c2)
    y = gN/a0
    win = (y >= Y_DEEP) & (y < ymax) & (w_pt > 0)
    if win.sum() < 50: return None
    x = np.sqrt(y[win])
    u = gobs[win]/np.sqrt(gN[win]*a0) - 1.0
    # error on u ~ ln(10)*sig*(1+u)
    su = math.log(10)*sig[win]*np.maximum(1+u, 0.1)
    ww = w_pt[win]/su**2
    X = np.stack([x, x*x], axis=1)
    XtW = X.T*ww
    beta = np.linalg.solve(XtW@X, XtW@u)
    resid = u - X@beta
    chi2 = np.sum(ww*resid**2)
    return a0, beta[0], beta[1], chi2, int(deep.sum()), int(win.sum())

allg = np.unique(gal_id)
NG = gal_id.max()+1
ones = np.zeros(NG); ones[allg] = 1.0
rng = np.random.default_rng(5)
for ymax in (0.3, 0.5, 1.0):
    # profile f_ML on the full sample
    fits = [(two_step(ymax, ones, f), f) for f in F_GRID]
    fits = [(r, f) for r, f in fits if r is not None]
    (a0, c1, c2, chi2, nd, nw), fbest = min(fits, key=lambda t: t[0][3])
    # bootstrap galaxies (f_ML re-profiled per replicate)
    c1b, c2b, a0b = [], [], []
    for k in range(400):
        pick = rng.choice(allg, len(allg), replace=True)
        wg = np.zeros(NG)
        for g_ in pick: wg[g_] += 1
        cand = [(two_step(ymax, wg, f), f) for f in F_GRID]
        cand = [(r, f) for r, f in cand if r is not None]
        if not cand: continue
        (a0_, c1_, c2_, _, _, _), _ = min(cand, key=lambda t: t[0][3])
        a0b.append(a0_); c1b.append(c1_); c2b.append(c2_)
    c1b, c2b, a0b = map(np.array, (c1b, c2b, a0b))
    e1, e2 = c1b.std(ddof=1), c2b.std(ddof=1)
    z_be = (c1-0.5)/e1; z_std = (c1-0.0)/e1
    L.append(f"ymax={ymax}: deep n={nd}, window n={nw}, f_ML={fbest:.1f}, "
             f"a0={a0:.3e} +/- {a0b.std(ddof=1):.2e}")
    L.append(f"  c1 = {c1:.3f} +/- {e1:.3f}   c2 = {c2:.3f} +/- {e2:.3f}")
    L.append(f"  c1 vs BE(0.5): {z_be:+.1f} sigma | vs standard(0): "
             f"{z_std:+.1f} sigma")
for l in L: print(l)
with open('data/stage4a_nlo.txt', 'w') as f_:
    f_.write("\n".join(L)+"\n")
print("\nsaved: data/stage4a_nlo.txt")
