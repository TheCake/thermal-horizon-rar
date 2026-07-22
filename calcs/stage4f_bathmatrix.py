"""
STAGE 4F: the bath matrix -- testing the self-consistent quantum bootstrap
(the 1/4-branch) that the Stage 4D/4E thermal reading generates.

The Planck-oscillator reading of the RAR admits a 2x2 matrix of bath laws
(mode frequency prescription x occupation statistics):

              | source-driven, w ~ sqrt(gN*a0) | self-consistent, w ~ g_tot
  quantum     | BE/RAR: nu = 1+n_BE(sqrt(y))   | BOOT: nu = 1+n_BE(nu*y)  <- NEW
  (Planck n)  |   c1=1/2, c2=1/12              |   c1=1/4, c2=7/96
  classical   | CSD: nu = 1+1/sqrt(y)          | nu = 1/2+sqrt(1/4+1/y)
  (n=kT/hbw)  |   c1=1, c2=0 (Cassini-hostile) |   == EXACTLY simple-nu

The classical/self-consistent cell IS simple-nu (two lines: nu = 1+1/(nu*y)
=> y*nu^2 - y*nu - 1 = 0 => nu = 1/2+sqrt(1/4+1/y)) -- so the BE-vs-simple
stalemate is secretly "quantum vs classical bath". The untested cell is the
quantum bootstrap, defined implicitly per point: nu = 1 + 1/(exp(nu*y)-1).
High-acceleration behavior: BOOT approaches Newton as e^-y (fastest
screening of all families); CSD approaches as 1/sqrt(y) (slowest -- dead by
Cassini/solar-system independently of any fit; included for completeness).

This script ranks all four cells + standard-mu (dead-branch control) under
BOTH likelihood treatments:
  (a) raw chi2, 4B objective (windows in y, (a0,f_ML) free), and
  (b) scatter-marginalized joint SPARC+lensing, 4E objective
      ((a0,f_ML,s_int,dlt) free, Mistele+24 block, dlt ~ N(0,0.2)).
Bootstrap sign-stability for the new pairings in both treatments.

Gates:
  G1 series: implemented solver must reproduce c1=1/4, c2=7/96 numerically.
  G2 solver: Newton residual < 1e-9 everywhere; nu(y=100)-1 < 1e-40;
     nu*sqrt(y) -> 1 at y=1e-10.
  G3 regression: BE/simple/standard fiducial -2lnL must reproduce
     stage4e_lensing.txt (same code path, same numbers).
Writes data/stage4f_bathmatrix.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25

# ---------------- SPARC (identical to 4A/4B/4E) ----------------
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

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask_fid = l_gbar >= LENS_CUT_FID

# ---------------- families ----------------
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_csd(y):
    return 1.0+1.0/np.sqrt(np.clip(y, 1e-14, None))
def nu_boot(y):
    """self-consistent quantum bath: nu = 1 + n_BE(nu*y), solved per point.
    Root of F(u) = u - y - y/(e^u - 1), u = nu*y. Seed u0 = simple-nu's u
    (exact root of the classical bootstrap) satisfies F(u0) > 0 with F
    monotone increasing => Newton descends cleanly onto the unique root."""
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
        'standard': nu_standard, 'CSD': nu_csd}
PAIR_NOTE = {'BE': 'quantum/source (c1=1/2, c2=1/12)',
             'simple': 'classical/self-consistent (c1=1/2, c2=1/8)',
             'boot': 'quantum/self-consistent (c1=1/4, c2=7/96)',
             'standard': 'dead-branch control (c1=0)',
             'CSD': 'classical/source (c1=1; Cassini-dead a priori)'}

L = [f"STAGE 4F bath matrix: {kept} galaxies, {len(gobs)} SPARC points + "
     f"{int(lmask_fid.sum())} lensing points (4E fiducial config)"]
for n_, note in PAIR_NOTE.items(): L.append(f"  {n_:>8}: {note}")

# ---------------- G1/G2 solver gates ----------------
xg = np.linspace(0.005, 0.05, 40)
yg = xg*xg
h = nu_boot(yg)*xg - 1.0          # = c1*x + c2*x^2 + ...
cfit = np.polyfit(xg, h/xg, 1)    # [c2, c1]
c1n, c2n = cfit[1], cfit[0]
g1ok = abs(c1n-0.25) < 2e-3 and abs(c2n-7/96) < 4e-3
u_test = nu_boot(np.logspace(-8, 1.6, 300))*np.logspace(-8, 1.6, 300)
yv = np.logspace(-8, 1.6, 300)
res = np.max(np.abs(u_test - yv - yv/np.expm1(np.minimum(u_test, 45))))
lo_lim = nu_boot(np.array([1e-10]))[0]*1e-5
hi_lim = nu_boot(np.array([100.0]))[0] - 1.0
L += ["", f"G1 series: c1 = {c1n:.4f} (pred 0.2500), c2 = {c2n:.4f} "
      f"(pred {7/96:.4f}) -> {'OK' if g1ok else 'FAIL'}",
      f"G2 solver: max|F| = {res:.2e}; nu*sqrt(y)@1e-10 = {lo_lim:.6f} "
      f"(->1); nu(100)-1 = {hi_lim:.1e} -> "
      f"{'OK' if res < 1e-9 and abs(lo_lim-1) < 1e-3 and hi_lim < 1e-40 else 'FAIL'}"]

# ---------------- objectives (verbatim 4E) ----------------
LN10 = math.log(10)
def m2ll(th, nu, sel, w_gal, lmask, lens_obs, use_lens=True):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    gN = (g_gas + f*g_dsk + g_bul)[sel]
    gm = gN*nu(gN/a0)
    se2 = sig2[sel] + s_int*s_int
    r = lgobs[sel] - np.log10(gm)
    out = np.sum(w_gal[gal_id][sel]*(r*r/se2 + np.log(se2)))
    if use_lens:
        lg = l_gbar[lmask] + dlt
        yl = 10**lg/a0
        lgm = lg + np.log10(nu(yl))
        rl = lens_obs[lmask] - lgm
        out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
        out += (dlt/DELTA_PRIOR)**2
    return out

def fit_joint(nu, sel, w_gal, lmask, lens_obs):
    best = None
    for th0 in ([math.log10(A0_FID), 1.0, 0.08, 0.0],
                [math.log10(A0_FID)+0.1, 0.8, 0.12, -0.1]):
        b = minimize(lambda t: m2ll(t, nu, sel, w_gal, lmask, lens_obs),
                     th0, method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
        if best is None or b.fun < best.fun: best = b
    return best

def fit_raw(nu, sel, w_pt):
    def chi2(th):
        la0, f = th
        if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
        gN = (g_gas + f*g_dsk + g_bul)[sel]
        r = (lgobs[sel] - np.log10(gN*nu(gN/10**la0)))/sig[sel]
        return np.sum(w_pt[sel]*r*r)
    return minimize(chi2, [math.log10(A0_FID), 1.0], method='Nelder-Mead',
                    options=dict(maxiter=2000, xatol=1e-7, fatol=1e-7))

allsel = np.ones(len(gobs), bool)
ones = np.ones(gal_id.max()+1)
y_fid = (g_gas + g_dsk + g_bul)/A0_FID
NG = gal_id.max()+1
allg = np.unique(gal_id)

# ---------------- (a) raw 4B treatment ----------------
L += ["", "(a) RAW chi2 (4B objective, (a0,f_ML) free):"]
raw = {}
for ymax, tag in ((0.5, 'y<0.5'), (1.0, 'y<1'), (30.0, 'y<30')):
    sel = y_fid < ymax
    raw[tag] = {}
    for name, nu in FAMS.items():
        b = fit_raw(nu, sel, np.ones(len(gobs)))
        raw[tag][name] = b
        L.append(f"  [{tag:5}] {name:>8}: chi2={b.fun:9.1f}  a0={10**b.x[0]:.3e}"
                 f"  f_ML={b.x[1]:.2f}"
                 + ("  <-- f_ML AT BOUND" if b.x[1] > 2.45 else ""))
    L.append(f"  [{tag:5}] boot - BE = {raw[tag]['boot'].fun-raw[tag]['BE'].fun:+9.1f} | "
             f"boot - simple = {raw[tag]['boot'].fun-raw[tag]['simple'].fun:+9.1f} | "
             f"boot - standard = {raw[tag]['boot'].fun-raw[tag]['standard'].fun:+9.1f}")
# raw bootstrap at y<1 (4B's strongest window)
rng = np.random.default_rng(29)
sel1 = y_fid < 1.0
w_bb, w_bs = 0, 0
for k in range(200):
    pick = rng.choice(allg, len(allg), replace=True)
    wg = np.zeros(NG)
    for g_ in pick: wg[g_] += 1
    w_pt = wg[gal_id]
    cb = fit_raw(nu_boot, sel1, w_pt).fun
    w_bb += int(cb < fit_raw(nu_be, sel1, w_pt).fun)
    w_bs += int(cb < fit_raw(nu_simple, sel1, w_pt).fun)
L.append(f"  raw bootstrap (y<1, 200 reps): boot beats BE in {w_bb}/200; "
         f"boot beats simple in {w_bs}/200")

# ---------------- (b) honest 4E treatment ----------------
L += ["", "(b) SCATTER-MARGINALIZED joint SPARC+lensing (4E fiducial):"]
fid = {}
for name, nu in FAMS.items():
    b = fit_joint(nu, allsel, ones, lmask_fid, l_gobs)
    fid[name] = b
    la0, f, s_int, dlt = b.x
    L.append(f"  {name:>8}: -2lnL = {b.fun:10.2f}  a0={10**la0:.3e}  "
             f"f_ML={f:.2f}  s_int={s_int:.3f}  dlt={dlt:+.3f}"
             + ("  <-- f_ML AT BOUND" if f > 2.45 else ""))
# G3 regression vs stage4e_lensing.txt values
REF = {'BE': -8397.72, 'simple': -8379.01, 'standard': -8341.95}
g3 = all(abs(fid[k].fun - v) < 0.5 for k, v in REF.items())
L.append(f"  G3 regression vs 4E fiducial: "
         + ", ".join(f"{k} {fid[k].fun:+.2f} (ref {v:+.2f})"
                     for k, v in REF.items())
         + f" -> {'OK' if g3 else 'FAIL'}")
L.append(f"  Delta(-2lnL): boot - BE = {fid['boot'].fun-fid['BE'].fun:+.2f} | "
         f"boot - simple = {fid['boot'].fun-fid['simple'].fun:+.2f} | "
         f"boot - standard = {fid['boot'].fun-fid['standard'].fun:+.2f} | "
         f"CSD - BE = {fid['CSD'].fun-fid['BE'].fun:+.2f}")

# honest bootstrap for {boot, BE, simple}
rng2 = np.random.default_rng(31)
d_bb, d_bs = [], []
for k in range(200):
    pick = rng2.choice(allg, len(allg), replace=True)
    wg = np.zeros(NG)
    for g_ in pick: wg[g_] += 1
    lo = l_gobs + rng2.normal(0, np.sqrt(l_sig2))
    rB = fit_joint(nu_boot, allsel, wg, lmask_fid, lo).fun
    rE = fit_joint(nu_be, allsel, wg, lmask_fid, lo).fun
    rS = fit_joint(nu_simple, allsel, wg, lmask_fid, lo).fun
    d_bb.append(rB - rE); d_bs.append(rB - rS)
d_bb, d_bs = np.array(d_bb), np.array(d_bs)
L += [f"  bootstrap (200 reps): boot - BE = {d_bb.mean():+.2f} +/- "
      f"{d_bb.std(ddof=1):.2f} (boot better in {int((d_bb<0).sum())}/200)",
      f"                        boot - simple = {d_bs.mean():+.2f} +/- "
      f"{d_bs.std(ddof=1):.2f} (boot better in {int((d_bs<0).sum())}/200)"]

out = "\n".join(L)
print(out)
with open('data/stage4f_bathmatrix.txt', 'w') as f_:
    f_.write(out+"\n")
print("\nsaved: data/stage4f_bathmatrix.txt")
