"""
STAGE 4T: the bath's second moment -- does the RAR's intrinsic scatter carry
the Planck-oscillator fluctuation signature?

If nu - 1 = n_BE(x) is a thermal occupation, the occupation FLUCTUATES:
Var(n) = n(n+1) for a Bose mode. A per-point draw from N effective modes
predicts relative scatter on g_obs of
    sigma_rel(x) = sqrt(n(n+1))/(1+n)/sqrt(N) = sqrt(n/(n+1))/sqrt(N),
i.e. a SHAPE: scatter rising monotonically toward deep MOND (-> 1/sqrt(N))
and dying as e^(-x/2)/sqrt(N) in the Newtonian regime. The RAR's observed
tightness then BOUNDS the mode count N from below; the scatter's x-shape
tests the reading.

Models for the intrinsic scatter s(x) (dex), fitted jointly with (a0, f_ML)
on SPARC (2,700 points; per-point reported errors always included):
  M0  constant s0                                   (3 params)
  M1  pure oscillator: s = sqrt(n/(n+1))/(sqrt(N)*ln10)      (3 params)
  M1b oscillator + constant floor (quadrature)      (4 params)
  M2  free per x-bin (6 fixed quantile bins)        (8 params)
Family = occupation law (BE); simple-nu repeated as robustness.

Gates:
  G1 nesting: -2lnL(M2) <= -2lnL(M1b) <= min(M0, M1) (optimizer sanity).
  G2 injection: synthetic residuals at N=30 recover N-hat within a factor 1.5.
Caveat stated: per-point independence assumed; a galaxy-level draw (modes
shared along a rotation curve) is the natural next refinement.
Writes data/stage4t_bathnoise.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)

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
sig2 = sig*sig
lgobs = np.log10(gobs)

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))

# fixed x-bins from the fiducial acceleration map (sextiles)
x_fid = np.sqrt((g_gas + g_dsk + g_bul)/A0_FID)
QE = np.quantile(x_fid, np.linspace(0, 1, 7))
QE[0], QE[-1] = 0.0, np.inf
BIN_FID = np.clip(np.searchsorted(QE, x_fid, side='right')-1, 0, 5)

def osc_shape(x):
    """sqrt(n/(n+1)) for n = 1/(e^x - 1): -> 1 deep, -> e^(-x/2) Newtonian."""
    xc = np.minimum(np.clip(x, 1e-9, None), 80.0)
    n = 1.0/np.expm1(xc)
    return np.sqrt(n/(1.0+n))

def m2ll(la0, f, s_dex, nu):
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    gm = gN*nu(gN/a0)
    se2 = sig2 + s_dex*s_dex
    r = lgobs - np.log10(gm)
    return np.sum(r*r/se2 + np.log(se2))

def fit(model, nu, x0):
    def obj(th):
        la0, f = th[0], th[1]
        if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
        x = np.sqrt((g_gas + f*g_dsk + g_bul)/10**la0)
        if model == 'M0':
            s0 = th[2]
            if not (1e-4 <= s0 < 0.5): return 1e12
            s = np.full(len(gobs), s0)
        elif model == 'M1':
            lnN = th[2]
            if not (-2 < lnN < 12): return 1e12
            s = osc_shape(x)/(math.sqrt(math.exp(lnN))*LN10)
        elif model == 'M1b':
            lnN, sf = th[2], th[3]
            if not (-2 < lnN < 12) or not (0 <= sf < 0.5): return 1e12
            s = np.sqrt((osc_shape(x)/(math.sqrt(math.exp(lnN))*LN10))**2 + sf*sf)
        else:  # M2
            sb = np.asarray(th[2:8])
            if np.any(sb < 1e-4) or np.any(sb > 0.5): return 1e12
            s = sb[BIN_FID]
        return m2ll(la0, f, s, nu)
    best = None
    for t0 in x0:
        b = minimize(obj, t0, method='Nelder-Mead',
                     options=dict(maxiter=8000, xatol=1e-6, fatol=1e-6))
        if best is None or b.fun < best.fun: best = b
    return best

L = [f"STAGE 4T bath second moment: {kept} galaxies, {len(gobs)} SPARC points",
     f"x-bin edges (fiducial): {np.round(QE[1:-1],3).tolist()}",
     ""]
LA, FM = math.log10(A0_FID), 1.0
for famname, nu in (('BE', nu_be), ('simple', nu_simple)):
    r = {}
    r['M0']  = fit('M0',  nu, [[LA, FM, 0.10], [LA, FM, 0.06]])
    r['M1']  = fit('M1',  nu, [[LA, FM, math.log(30.0)], [LA, FM, math.log(8.0)]])
    r['M1b'] = fit('M1b', nu, [[LA, FM, math.log(30.0), 0.05],
                               [LA, FM, math.log(8.0), 0.10]])
    r['M2']  = fit('M2',  nu, [[LA, FM]+[0.10]*6, [LA, FM]+[0.06]*6])
    L.append(f"[{famname}]")
    m0, m1, m1b, m2 = (r[k].fun for k in ('M0', 'M1', 'M1b', 'M2'))
    L.append(f"  M0  const:          -2lnL = {m0:10.2f}  s0 = {r['M0'].x[2]:.4f}")
    L.append(f"  M1  oscillator:     -2lnL = {m1:10.2f}  N-hat = "
             f"{math.exp(r['M1'].x[2]):8.1f}")
    L.append(f"  M1b osc + floor:    -2lnL = {m1b:10.2f}  N-hat = "
             f"{math.exp(r['M1b'].x[2]):8.1f}  floor = {r['M1b'].x[3]:.4f}")
    sb = r['M2'].x[2:8]
    L.append(f"  M2  free 6-bin:     -2lnL = {m2:10.2f}  s_b = "
             f"{np.round(sb,4).tolist()}")
    L.append(f"  D(-2lnL): M1-M0 = {m1-m0:+.2f} | M1b-M0 = {m1b-m0:+.2f} "
             f"(1 extra param) | M2-M0 = {m2-m0:+.2f} (5 extra)")
    g1 = (m2 <= m1b + 1e-6) and (m1b <= min(m0, m1) + 1e-6)
    L.append(f"  G1 nesting: {'PASS' if g1 else 'FAIL'}")
    # oscillator curve at M1b fit vs the free bins
    la0b, fb = r['M1b'].x[0], r['M1b'].x[1]
    xb = np.sqrt((g_gas + fb*g_dsk + g_bul)/10**la0b)
    Nb, sfb = math.exp(r['M1b'].x[2]), r['M1b'].x[3]
    L.append(f"  per-bin: median x | osc+floor pred | free-fit s_b")
    for b in range(6):
        m = BIN_FID == b
        xm = np.median(xb[m])
        pred = math.sqrt((osc_shape(np.array([xm]))[0]/(math.sqrt(Nb)*LN10))**2
                         + sfb*sfb)
        L.append(f"    bin{b}: {xm:6.3f} | {pred:.4f} | {sb[b]:.4f}")
    L.append("")

# G2 injection: synthesize residuals with oscillator scatter at N=30 (BE truth)
rng = np.random.default_rng(4)
la0t, ft, Ntrue = math.log10(A0_FID), 1.0, 30.0
gNt = g_gas + ft*g_dsk + g_bul
xt = np.sqrt(gNt/10**la0t)
s_true = osc_shape(xt)/(math.sqrt(Ntrue)*LN10)
lgobs_true = np.log10(gNt*nu_be(gNt/10**la0t))
lgobs_save = lgobs.copy()
lgobs = lgobs_true + rng.normal(0, np.sqrt(sig2 + s_true**2))
binj = fit('M1', nu_be, [[la0t, ft, math.log(10.0)], [la0t, ft, math.log(60.0)]])
Nrec = math.exp(binj.x[2])
lgobs = lgobs_save
g2 = 1/1.5 < Nrec/Ntrue < 1.5
L.append(f"G2 injection (truth N=30): N-hat = {Nrec:.1f} -> "
         f"{'PASS' if g2 else 'FAIL'}")
L.append("")
L.append("caveat: per-point independence assumed; galaxy-level correlated draw")
L.append("(shared modes along a curve) is the next refinement and would loosen N.")

out = "\n".join(L)
print(out)
with open('data/stage4t_bathnoise.txt', 'w') as f:
    f.write(out+"\n")
