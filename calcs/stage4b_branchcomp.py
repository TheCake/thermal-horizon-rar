"""
STAGE 4B: branch comparison without expansion truncation. The 4A expansion
estimator is power-limited (c1 to +/-0.4-0.6; windows see-saw). Direct test:
fit the FULL nu functions on the same low/intermediate-acceleration points
with only (a0, f_ML) free, compare chi2:
   BE:       nu = 1/(1-exp(-sqrt(y)))          [c1=1/2 branch]
   simple:   nu = 1/2+sqrt(1/4+1/y)            [c1=1/2 branch]
   standard: nu = sqrt((1+sqrt(1+4/y^2))/2)    [c1=0 branch]
Galaxy bootstrap gives the sign-stability of Delta-chi2 (BE - standard).
Writes data/stage4b_branch.txt.
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

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
FAMS = {'BE': nu_be, 'simple': nu_simple, 'standard': nu_standard}

L = [f"STAGE 4B branch comparison: {kept} galaxies, {len(gobs)} points",
     "full-nu fits, (a0, f_ML) free; windows in y (fiducial a0)"]

def fit_fam(nu, sel, w_pt):
    def chi2(th):
        la0, f = th
        a0 = 10**la0
        gN = (g_gas + f*g_dsk + g_bul)[sel]
        gm = gN*nu(gN/a0)
        r = (np.log10(gobs[sel]) - np.log10(gm))/sig[sel]
        return np.sum(w_pt[sel]*r*r)
    b = minimize(chi2, [math.log10(A0_FID), 1.0], method='Nelder-Mead',
                 options=dict(maxiter=2000, xatol=1e-7, fatol=1e-7))
    return b.fun, 10**b.x[0], b.x[1]

allg = np.unique(gal_id)
NG = gal_id.max()+1
rng = np.random.default_rng(11)
y_fid = (g_gas + g_dsk + g_bul)/A0_FID
for ymax in (0.5, 1.0, 30.0):
    sel = y_fid < ymax
    tag = f"y<{ymax}" if ymax < 10 else "all y"
    ones = np.ones(len(gobs))
    res = {}
    for name, nu in FAMS.items():
        c2, a0, f = fit_fam(nu, sel, ones)
        res[name] = (c2, a0, f)
        L.append(f"  [{tag}] {name:>8}: chi2={c2:9.1f}  a0={a0:.3e}  "
                 f"f_ML={f:.2f}  (n={sel.sum()})")
    dbs = res['BE'][0]-res['standard'][0]
    dsim = res['BE'][0]-res['simple'][0]
    # bootstrap sign stability of BE-vs-standard
    wins = 0; tot = 0
    for k in range(200):
        pick = rng.choice(allg, len(allg), replace=True)
        wg = np.zeros(NG)
        for g_ in pick: wg[g_] += 1
        w_pt = wg[gal_id]
        cb,_,_ = fit_fam(nu_be, sel, w_pt)
        cs,_,_ = fit_fam(nu_standard, sel, w_pt)
        wins += int(cb < cs); tot += 1
    L.append(f"  [{tag}] Delta-chi2 (BE - standard) = {dbs:+.1f} "
             f"(negative = BE better); BE better in {wins}/{tot} bootstraps")
    L.append(f"  [{tag}] Delta-chi2 (BE - simple)   = {dsim:+.1f}")
for l in L: print(l)
with open('data/stage4b_branch.txt', 'w') as f_:
    f_.write("\n".join(L)+"\n")
print("\nsaved: data/stage4b_branch.txt")
