"""
STAGE 10D -- O5-KAPPA: THE JOINT SINGLE-KAPPA FIT (ROUND-24 condition 8;
author call 2026-08-09 "Go for the next arc close... Hunt continues").

QUESTION. 10C's transmission map is nu = 1 + kappa*n_BE(x) with kappa =
4g^2/(Om om) the mechanism's ONE free constant. The sky reads kappa
two-pole: the a0-lock (deep normalization vs the horizon temperature)
says 1.00 +/- 0.05; the 4Z hier transition-constant meter says 1.484
(1.382-1.584) at profile grade (bootstrap band 0.6-1.8 contains 1).
Both meters are PROJECTIONS of fits made in OTHER families (BE fixed-
form; the 4S/4Z lambda-mixture family) translated through the series
map kappa = 2(1-c1). This stage runs the FIRST DIRECT kappa-family
fit: deep + transition, one kappa, one likelihood -- the two channels
(deep amplitude kappa^2*a0; transition constant c1 = 1 - kappa/2)
confront each other inside one chi^2, on the SAME hierarchical
apparatus that produced the tension (4Z plain hier + lensing; 5M
vertical-hardened co-read).

FAMILY. nu_kap(y; kappa) = 1 + kappa*n_BE(sqrt(y)), n_BE(x) =
1/(e^x - 1), y = g_N/a0. At kappa = 1 this is EXACTLY the 4Z BE
endpoint (identity point; G0). Deep: nu -> kappa/x, so the deep arm
measures kappa^2*a0; transition (renormalized to unit deep
coefficient, xt = x/kappa): c1 = 1 - kappa/2, c2 = kappa^2/12 --
the series map used by the 10C meters, validated here at G1.

FITS (plain-hier apparatus = 4Z bit-verbatim model expressions;
per-galaxy disk-M/L offsets, 0.1-dex lognormal prior, joint SPARC +
lensing marginalized objective; fitting WRAPPER passes lgobs
explicitly instead of the 4Z global swap -- wrapper-level deviation,
model expressions untouched):
  F3  kappa = 1 fixed, a0 free ......... the C&T baseline (= archived
      4Z lam=1.00 hier node -10435.06; regression G0).
  F1  kappa free, a0 free .............. THE PRIMARY MEASUREMENT
      (kappa identified by shape alone; deep arm kappa-blind at free
      a0 by the kappa^2*a0 degeneracy -- disclosed).
  F2  kappa free, a0 horizon-locked .... the PHYSICS fit (temperature
      prior: la0 ~ N(log10 1.043e-10 [cH0/2pi Planck, the 4V/5L
      convention], 0.018 dex = half the Planck-SH0ES spread); SH0ES-
      centered alternate row printed).
  F4  split-kappa (kappa_d, kappa_t, kappa_l) over fiducial regime
      bins x_fid = sqrt(g_N,fid/1.2e-10) with g_N,fid at f=1, dml=0
      (DATA-DEFINED, fit-independent): deep x<0.7 / transition
      0.7-2.5 / tail >=2.5. Lensing points (all deep) join the deep
      bin. kappa_l is a PROXY row (the tail rides the tail-EXPONENT
      physics, 5G/9U business, NOT part of the R24-cond-8 contest --
      disclosed iv). Nested: F4 at equal kappas == F1 (identity G3).

INSTRUMENTS AND GATES (trap #12: every letter clause gated):
  G0  IDENTITY (STOP-class): (a) max rel |nu_kap(y,1) - nu_be(y)| over
      the data y-range <= 5e-12; (b) co-evaluated objective at one
      theta: |m2k(kappa=1) - m2h_BE| <= 0.01 (lnL grade, the 8P/8Q
      rule); (c) F3 refit reproduces -10435.06 within 1.0 AND the
      pd->1e-4 endpoint reproduces the 4S lam=1 value -8397.72 within
      1.0 (4Z G1 bars).
  G1  SERIES (sympy, exact): renormalized ladder of nu_kap: c1 =
      1 - kappa/2 and c2 = kappa^2/12 exact => the 10C translation
      map is validated at series level.
  G2  INJECTIONS (4Z G2 recipe verbatim-adapted; real lensing points
      retained -- disclosed v): (a) truth kappa=1 (BE + drawn offsets
      + noise): F1 recovers kappa in [0.85, 1.15]; (b) truth
      kappa=1.5 (nu_kap, a0 = 1.2e-10): F1 recovers in [1.35, 1.65];
      (c) SPLIT-NULL: on the kappa=1 truth, D(F1-F4) < 11.8 (the
      false-positive control for the K-SPLIT clause). (a) or (c)
      fail => STOP-class (machinery bug). (b) fail => kappa-family
      recovery flagged, letter blocked to TENSE.
  G3  NESTING/OPT: F4(equal kappas) == F1 objective at the same theta
      (<= 1e-6); F1 multi-start (kappa0 in {0.7, 1.0, 1.5}) best-two
      spread <= 1.0; F4 best <= F1 best + 0.01 (nesting). Fail =>
      letters blocked (report-only).
  G4  BOOTSTRAP (PRIMARY UNCERTAINTY -- trap #7: no grid-point
      bounds): 40 paired galaxy-resample reps (rng seed 47, fresh,
      disclosed), each rep fits F1 AND F4 (2 rounds, warm-started);
      report kappa-hat 16/50/84, (kappa_t - kappa_d) 16/50/84,
      D(F1-F4) percentiles + counts.
  G5  DECOMP: per-bin -2lnL decomposition of F1 vs F4 at the optima
      (where does the one-kappa form pay?) + boundary variants
      {0.5, 2.0} and {1.0, 3.0} (round-10 risk-axis rule).
  G6  VERTICAL CO-READ (point estimates only, no bootstrap --
      disclosed vi): 5M apparatus (measured sigma_v, closed-form dv)
      with the kappa family; regression: kappa=1 dv-ON reproduces the
      archived -12152.49 within 1.5; then F1v, F4v. Lean-grade.

LETTER GRAMMAR (pre-registered; bars locked):
  STOP-class: G0 or G2(a)/(c) fail => abort, no verdict (bug).
  K-SPLIT: point D(F1-F4) >= 11.8 AND boot median D >= 11.8 AND
      >= 32/40 reps D >= 6.2 => the one-kappa form is REJECTED
      between deep and transition -- kappa RUNS with regime (the
      running profile is the measurement; kill for a constant-kappa
      closure; the kappa=1 derivation target re-scopes to the deep
      limit).
  K-TENSE: point D >= 6.2 but K-SPLIT bars not met => strain
      reported, one-kappa carried (lean language only).
  K-ONE: K-SPLIT/K-TENSE not fired AND F1 bootstrap 16-84 of
      kappa-hat CONTAINS 1.0 => one kappa suffices and the closure
      value is inside the measurement.
  K-ONE-OFF: split not fired AND the 16-84 band EXCLUDES 1.0 =>
      one kappa suffices but NOT the closure value (the kappa = 1
      derivation target dies at this grade; report the value).
  POWER flag: boot 16-84 width of (kappa_t - kappa_d) > 1.0 =>
      POWER-LIMITED printed alongside any letter (disclosed).
  The F2 lock row and the G6 vertical row print as co-reads in every
  cell (mandatory tension rows, 10C G5 precedent).

PRE-SIGNED CREDENCE MAP (arc-level, shared with 10E -- mechanical):
  10D K-ONE (contains 1) + ROUND-25 no unpatched hole => bath-
      mechanism conditional 15 -> 17. 10E V-RABI + clean round =>
      15 -> 18. BOTH => 15 -> 20.
  10D K-ONE-OFF or K-SPLIT => HOLD 15 (structure informative; the
      kappa=1 target re-scopes; NO strike from 10D alone).
  10E V-WEAK or V-DIAGONAL (reviewer-affirmed) => 15 -> 8 (the
      R22-cond-4 pre-registered strike) REGARDLESS of 10D.
  Any unpatched hole blocks that stage's cell. anomaly-real 53
  UNTOUCHED by both stages (no binary-anomaly evidence in play).

DISCLOSURES: (i) first direct kappa-family fit; archived "hier kappa
= 1.484" was a cross-family series translation. (ii) deep-arm
kappa^2*a0 degeneracy at free a0; F1's kappa is shape-identified;
the lock enters only via F2. (iii) bins data-defined at fiducial;
boundary variants shipped. (iv) kappa_l = proxy row. (v) injection
retains real lensing (4Z-verbatim recipe). (vi) vertical co-read
point-grade. (vii) bootstrap primary everywhere (4S/9U population-
variance lessons).

AMENDMENT A1 (2026-08-09, pre-quote, bug-class; run-1 archived as
data/stage10d_kappa_run1.txt; committed BEFORE the re-run): run-1
fired STOP on G0(c) (F3 cold refit -10432.80 vs archived -10435.06,
2.26 past the 1.0 bar) and G2(a) (kappa=1 truth recovered 1.167 vs
[0.85, 1.15]). Diagnosis: fixed-round coordinate descent
under-converges on the known a0-f_ML ridge (the archived node was
warm-started along the 4Z sweep; the 5M cold BE fit needed its
tol-looped 12+-round descent to land -10434.99). Fixes (machinery
only, NO bar moved): (1) fit_plain gains the 5M-style tol-looped
convergence (tol 0.05; max_rounds 15 point / 10 injection / 4
bootstrap) + a final NM polish; (2) theta multi-starts gain the
second ridge corner (a0 ~ 0.8e-10, f ~ 1.4) alongside the fiducial;
(3) the kappa guard widens 0.2 -> 0.02 because run-1's tail proxy
kappa_l = 0.200 SAT ON the guard (grid-edge riding, 7I discipline;
the letter bars never referenced the guard). Injection draws are
rng-stream-identical (no rng calls added or removed). Bars, bins,
letters, credence map: UNCHANGED.

Writes data/stage10d_kappa.txt.
"""
import glob, math, os, time
import numpy as np
import sympy as sp
from scipy.optimize import minimize, minimize_scalar

T00 = time.time()
KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25
LN10 = math.log(10)
S_ML = 0.1*LN10
A0_HOR_PL = 1.043e-10
A0_HOR_SH = 1.130e-10
S_HOR_DEX = 0.018
XB1, XB2 = 0.7, 2.5

OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

# ---------------- loader (5M form = 4Z + sigma_v extension) ----------------
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
LGOBS = np.log10(gobs)
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

# fiducial regime bins (data-defined; fit-independent)
gN_fid = g_gas + 1.0*g_dsk + g_bul
x_fid = np.sqrt(gN_fid/A0_FID)
def make_bins(b1, b2):
    bsp = np.where(x_fid < b1, 0, np.where(x_fid < b2, 1, 2))
    xl = np.sqrt(10**l_gbar[lmask]/A0_FID)
    ble = np.where(xl < b1, 0, np.where(xl < b2, 1, 2))
    return bsp, ble
BSP, BLE = make_bins(XB1, XB2)

# ---------------- families ----------------
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def occ_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 0.0, 1.0/np.expm1(np.minimum(x, 40)))
def nu_kap(y, kap):
    return 1.0 + kap*occ_be(y)
def nu_kap3(y, b, kd, kt, kl):
    ka = np.where(b == 0, kd, np.where(b == 1, kt, kl))
    return 1.0 + ka*occ_be(y)

# ---------------- objectives (4Z m2h bit-verbatim expressions) -------------
def m2core(kap_sp, kap_le, th, dml, w_g, lgv, prior_hor=None):
    """kap_sp: per-SPARC-point kappa array (or scalar); kap_le: lensing."""
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm = gN*(1.0 + kap_sp*occ_be(gN/a0))
    se2 = sig2 + s_int*s_int
    r = lgv - np.log10(gm)
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    yl = 10**lg/a0
    rl = l_gobs[lmask] - (lg + np.log10(1.0 + kap_le*occ_be(yl)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    if prior_hor is not None:
        out += ((la0 - math.log10(prior_hor))/S_HOR_DEX)**2
    return out

def kaps_from(mode, tk, bsp, ble):
    if mode == 'one':
        return tk[0], tk[0]
    if mode == 'fix1':
        return 1.0, 1.0
    kd, kt, kl = tk
    ksp = np.where(bsp == 0, kd, np.where(bsp == 1, kt, kl))
    kle = np.where(ble == 0, kd, np.where(ble == 1, kt, kl))
    return ksp, kle

NK = {'one': 1, 'fix1': 0, 'three': 3}

def m2k(t, mode, dml, w_g, lgv, bsp, ble, prior_hor=None):
    nk = NK[mode]
    tk, th = t[:nk], t[nk:]
    if mode != 'fix1' and not all(0.02 < k < 3.0 for k in tk): return 1e12
    ksp, kle = kaps_from(mode, tk, bsp, ble)
    return m2core(ksp, kle, th, dml, w_g, lgv, prior_hor)

def fit_plain(mode, w_g, pd=S_ML, starts=None, rounds=3, lgv=None,
              bsp=None, ble=None, prior_hor=None, tol=0.05,
              max_rounds=15):
    if lgv is None: lgv = LGOBS
    if bsp is None: bsp, ble = BSP, BLE
    nk = NK[mode]
    dml = np.zeros(NGal)
    best, prev = None, None
    for rd in range(max_rounds):
        ss = ([list(best.x)] if best is not None else []) + (starts or [])
        bb = None
        for t0 in ss:
            b = minimize(lambda t: m2k(t, mode, dml, w_g, lgv, bsp, ble,
                                       prior_hor), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000+1500*nk, xatol=1e-6,
                                      fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        tk, (la0, f, s_int, dlt) = best.x[:nk], best.x[nk:]
        if pd < 1e-3:
            if rd >= rounds - 1: break
            continue
        ksp, _ = kaps_from(mode, tk, bsp, ble)
        se2c = s_int*s_int
        for gi2 in range(NGal):
            if w_g[gi2] == 0: dml[gi2] = 0.0; continue
            mm = GIDXS[gi2]
            kk = ksp if np.isscalar(ksp) else ksp[mm]
            def od(dl):
                fc = f*math.exp(dl)
                gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                rr = lgv[mm] - np.log10(gN2*(1.0 + kk*occ_be(gN2/10**la0)))
                s2 = sig2[mm] + se2c
                return (w_g[gi2]*np.sum(rr*rr/s2)
                        + w_g[gi2]*dl*dl/(pd*pd))
            dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                       method='bounded').x
        cur = m2k(best.x, mode, dml, w_g, lgv, bsp, ble, prior_hor)
        if rd >= rounds - 1 and prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2k(t, mode, dml, w_g, lgv, bsp, ble,
                               prior_hor), list(best.x),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000+1500*nk, xatol=1e-6,
                              fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml

ONES = np.ones(NGal)
say("="*78)
say(f"STAGE 10D -- O5-KAPPA: joint single-kappa fit; {kept} galaxies, "
    f"{len(gobs)} points + {int(lmask.sum())} lensing")
say(f"bins at x_fid = sqrt(g_N,fid/1.2e-10): deep <{XB1} "
    f"(N={int((BSP==0).sum())}+{int((BLE==0).sum())}L), trans {XB1}-{XB2} "
    f"(N={int((BSP==1).sum())}+{int((BLE==1).sum())}L), tail >={XB2} "
    f"(N={int((BSP==2).sum())}+{int((BLE==2).sum())}L)")
say("="*78)

# ---------------- G0 identity ----------------
say("")
say("G0 IDENTITY (STOP-class):")
yy = np.concatenate([gN_fid/A0_FID, 10**l_gbar[lmask]/A0_FID])
rel = np.max(np.abs(nu_kap(yy, 1.0) - nu_be(yy))/nu_be(yy))
ok0a = rel <= 5e-12
say(f"  (a) family identity at kappa=1: max rel diff = {rel:.2e} "
    f"(bar 5e-12) -> {'PASS' if ok0a else 'FAIL'}")
th_test = np.array([math.log10(A0_FID), 1.0, 0.08, 0.0])
dml0 = np.zeros(NGal)
v_kap = m2k(np.array([1.0] + list(th_test)), 'one', dml0, ONES, LGOBS,
            BSP, BLE)
# 4Z m2h co-evaluation (verbatim expressions, nu_be path)
def m2h_be(th, dml, w_g):
    la0, f, s_int, dlt = th
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm = gN*nu_be(gN/a0)
    se2 = sig2 + s_int*s_int
    r = LGOBS - np.log10(gm)
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    yl = 10**lg/a0
    rl = l_gobs[lmask] - (lg + np.log10(nu_be(yl)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    return out
v_be = m2h_be(th_test, dml0, ONES)
ok0b = abs(v_kap - v_be) <= 0.01
say(f"  (b) objective identity (lnL grade): |{v_kap:.4f} - {v_be:.4f}| = "
    f"{abs(v_kap-v_be):.2e} (bar 0.01) -> {'PASS' if ok0b else 'FAIL'}")
b3, dml3 = fit_plain('fix1', ONES, rounds=6,
                     starts=[[math.log10(A0_FID), 1.0, 0.08, 0.0],
                             [math.log10(8e-11), 1.4, 0.07, 0.0]])
ok0c1 = abs(b3.fun - (-10435.06)) < 1.0
say(f"  (c) F3 (kappa=1) hier refit: {b3.fun:.2f} (4Z lam=1 node "
    f"-10435.06) -> {'PASS' if ok0c1 else 'FAIL'}  [a0 = "
    f"{10**b3.x[0]:.3e}, f_ML = {b3.x[1]:.2f}, s_int = {b3.x[2]:.3f}]")
b3e, _ = fit_plain('fix1', ONES, pd=1e-4, rounds=1,
                   starts=[[math.log10(A0_FID), 1.0, 0.08, 0.0]])
ok0c2 = abs(b3e.fun - (-8397.72)) < 1.0
say(f"      pd->1e-4 endpoint: {b3e.fun:.2f} (4S lam=1 -8397.72) -> "
    f"{'PASS' if ok0c2 else 'FAIL'}")
G0 = ok0a and ok0b and ok0c1 and ok0c2
say(f"  G0 -> {'PASS' if G0 else 'FAIL (STOP)'}")

# ---------------- G1 series ----------------
say("")
say("G1 SERIES (sympy):")
xs, ks = sp.symbols('x kappa', positive=True)
nu_s = 1 + ks/(sp.exp(xs) - 1)
xt = sp.symbols('xt', positive=True)
ser = sp.series(nu_s.subs(xs, ks*xt), xt, 0, 3).removeO().expand()
c1_s = sp.simplify(ser.coeff(xt, 0) - (1 - ks/2))
c2_s = sp.simplify(ser.coeff(xt, 1) - ks**2/12)
deep = sp.simplify(sp.limit(nu_s*xs, xs, 0) - ks)
ok1 = (c1_s == 0) and (c2_s == 0) and (deep == 0)
say(f"  deep coeff = kappa: {'exact' if deep == 0 else 'FAIL'}; "
    f"renormalized c1 - (1-kappa/2) = {c1_s}; c2 - kappa^2/12 = {c2_s}")
say(f"  G1 -> {'PASS' if ok1 else 'FAIL'}  [the kappa = 2(1-c1) meter "
    f"map is series-exact]")

# ---------------- point fits ----------------
say("")
say("POINT FITS (plain hier):")
S1 = [[k, math.log10(A0_FID), 1.0, 0.08, 0.0] for k in (0.7, 1.0, 1.5)]
S1R = [[k, math.log10(8e-11), 1.4, 0.07, 0.0] for k in (0.7, 1.0, 1.5)]
f1_runs = []
for s0, s0r in zip(S1, S1R):
    b, dm = fit_plain('one', ONES, starts=[s0, s0r], rounds=4)
    f1_runs.append((b, dm))
f1_runs.sort(key=lambda r: r[0].fun)
b1, dml1 = f1_runs[0]
spread_f1 = f1_runs[1][0].fun - f1_runs[0][0].fun
kap1 = b1.x[0]
say(f"  F1 (kappa, a0 free): kappa = {kap1:.3f}, a0 = {10**b1.x[1]:.3e}, "
    f"f_ML = {b1.x[2]:.2f}, s_int = {b1.x[3]:.3f}, obj = {b1.fun:.2f}")
say(f"     multi-start best-two spread = {spread_f1:.2f} "
    f"(bar 1.0: {'OK' if spread_f1 <= 1.0 else 'AMBIG'})")
say(f"     implied c1 = 1 - kappa/2 = {1-kap1/2:.3f}  "
    f"[4Z hier c1 = 0.258 (0.208-0.309); 4S flat c1 = 0.450]")
say(f"     implied deep norm kappa^2*a0 = {kap1**2*10**b1.x[1]:.3e}")

b2, dml2 = fit_plain('one', ONES, starts=[[kap1] + list(b1.x[1:]),
                                          [1.0, math.log10(A0_HOR_PL),
                                           1.0, 0.08, 0.0]],
                     prior_hor=A0_HOR_PL)
say(f"  F2 (horizon-locked, Planck): kappa = {b2.x[0]:.3f}, a0 = "
    f"{10**b2.x[1]:.3e}, obj = {b2.fun:.2f} (D vs F1 {b2.fun-b1.fun:+.2f})")
b2s, _ = fit_plain('one', ONES, starts=[[b2.x[0]] + list(b2.x[1:]),
                                        [1.0, math.log10(A0_HOR_SH),
                                         1.0, 0.08, 0.0]],
                   prior_hor=A0_HOR_SH)
say(f"      (SH0ES-centered alternate: kappa = {b2s.x[0]:.3f}, a0 = "
    f"{10**b2s.x[1]:.3e})")
say(f"  F3 (kappa=1): a0 = {10**b3.x[0]:.3e}, obj = {b3.fun:.2f} "
    f"(D vs F1 {b3.fun-b1.fun:+.2f})")

S4 = [[kap1, kap1, kap1] + list(b1.x[1:]),
      [1.0, 1.5, 1.0, math.log10(A0_FID), 1.0, 0.08, 0.0],
      [1.5, 1.5, 1.5, math.log10(A0_FID), 1.0, 0.08, 0.0],
      [1.3, 1.0, 0.5, math.log10(8e-11), 1.4, 0.07, 0.0]]
b4, dml4 = fit_plain('three', ONES, starts=S4, rounds=4)
kd, kt, kl = b4.x[:3]
say(f"  F4 (split): kappa_d = {kd:.3f}, kappa_t = {kt:.3f}, "
    f"kappa_l = {kl:.3f} (proxy), a0 = {10**b4.x[3]:.3e}, obj = {b4.fun:.2f}")
D_split = b1.fun - b4.fun
ok3n = b4.fun <= b1.fun + 0.01
v4eq = m2k(np.array([kap1, kap1, kap1] + list(b1.x[1:])), 'three', dml1,
           ONES, LGOBS, BSP, BLE)
ok3i = abs(v4eq - m2k(np.array([kap1] + list(b1.x[1:])), 'one', dml1,
                      ONES, LGOBS, BSP, BLE)) <= 1e-6
say(f"  D(F1 - F4) = {D_split:+.2f}   [K-SPLIT point bar 11.8; "
    f"K-TENSE bar 6.2]")
say(f"  G3: F4-equal==F1 identity {'PASS' if ok3i else 'FAIL'}; nesting "
    f"{'PASS' if ok3n else 'FAIL'}; F1 opt {'OK' if spread_f1 <= 1.0 else 'AMBIG'}")
G3 = ok3i and ok3n and (spread_f1 <= 1.0)

# per-bin decomposition at the two optima
say("")
say("G5 DECOMP (per-bin SPARC -2lnL contribution, F1 minus F4):")
def bin_terms(bx, dmlx, mode):
    nk = NK[mode]
    tk, (la0, f, s_int, dlt) = bx.x[:nk], bx.x[nk:]
    ksp, _ = kaps_from(mode, tk, BSP, BLE)
    fac = f*np.exp(dmlx[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm = gN*(1.0 + ksp*occ_be(gN/10**la0))
    se2 = sig2 + s_int*s_int
    r = LGOBS - np.log10(gm)
    tt = r*r/se2 + np.log(se2)
    return np.array([tt[BSP == k].sum() for k in range(3)])
t1b, t4b = bin_terms(b1, dml1, 'one'), bin_terms(b4, dml4, 'three')
for k, nm in enumerate(('deep', 'trans', 'tail')):
    say(f"    {nm:>5}: F1 {t1b[k]:9.2f}  F4 {t4b[k]:9.2f}  "
        f"D {t1b[k]-t4b[k]:+7.2f}")
say("  boundary variants (F4 point refits):")
for bb1, bb2 in ((0.5, 2.0), (1.0, 3.0)):
    bspv, blev = make_bins(bb1, bb2)
    bv, _ = fit_plain('three', ONES, starts=[list(b4.x), S4[1]],
                      bsp=bspv, ble=blev, max_rounds=8)
    say(f"    bins ({bb1},{bb2}): kappa_d/t/l = {bv.x[0]:.3f}/"
        f"{bv.x[1]:.3f}/{bv.x[2]:.3f}, D vs F1 = {b1.fun-bv.fun:+.2f}")

# ---------------- G2 injections ----------------
say("")
say("G2 INJECTIONS (4Z recipe; real lensing retained -- disclosed):")
rng = np.random.default_rng(47)
def inject(nu_truth, a0_truth):
    d_t = rng.normal(0, S_ML, NGal)
    gN_t = g_gas + 1.0*np.exp(d_t[gidx])*g_dsk + g_bul
    return (np.log10(gN_t*nu_truth(gN_t/a0_truth))
            + rng.normal(0, np.sqrt(sig2 + 0.08**2)))
lg_i1 = inject(nu_be, A0_FID)
bi1, _ = fit_plain('one', ONES, starts=S1[:2] + S1R[:2], rounds=2,
                   max_rounds=10, lgv=lg_i1)
ok2a = 0.85 <= bi1.x[0] <= 1.15
say(f"  (a) truth kappa=1: recovered {bi1.x[0]:.3f} "
    f"(bar [0.85, 1.15]) -> {'PASS' if ok2a else 'FAIL'}")
lg_i2 = inject(lambda y: nu_kap(y, 1.5), A0_FID)
bi2, _ = fit_plain('one', ONES, starts=[[1.5, math.log10(A0_FID), 1.0,
                                         0.08, 0.0], S1[1]],
                   rounds=2, max_rounds=10, lgv=lg_i2)
ok2b = 1.35 <= bi2.x[0] <= 1.65
say(f"  (b) truth kappa=1.5: recovered {bi2.x[0]:.3f}, a0 "
    f"{10**bi2.x[1]:.3e} (bar [1.35, 1.65]) -> {'PASS' if ok2b else 'FAIL'}")
bi4, _ = fit_plain('three', ONES, starts=[[bi1.x[0]]*3 + list(bi1.x[1:]),
                                          S4[1]], rounds=2, max_rounds=10,
                   lgv=lg_i1)
D_null = bi1.fun - bi4.fun
ok2c = D_null < 11.8
say(f"  (c) split-null on kappa=1 truth: D(F1-F4) = {D_null:+.2f} "
    f"(bar < 11.8) -> {'PASS' if ok2c else 'FAIL'}")
G2 = ok2a and ok2c
say(f"  G2 -> {'PASS' if G2 else 'FAIL (STOP)'}"
    f"{'' if ok2b else '  [(b) fail -> letters blocked to TENSE]'}")

# ---------------- G4 bootstrap ----------------
say("")
say("G4 BOOTSTRAP (40 paired reps, seed 47):")
allg = np.arange(NGal)
K1B, KDB, KTB, KLB, DB = [], [], [], [], []
for rep in range(40):
    pick = rng.choice(allg, NGal, replace=True)
    w = np.zeros(NGal)
    for g_ in pick: w[g_] += 1
    br1, _ = fit_plain('one', w, starts=[[kap1] + list(b1.x[1:])],
                       rounds=2, max_rounds=4)
    br4, _ = fit_plain('three', w,
                       starts=[[br1.x[0]]*3 + list(br1.x[1:]),
                               list(b4.x)], rounds=2, max_rounds=4)
    K1B.append(br1.x[0]); KDB.append(br4.x[0]); KTB.append(br4.x[1])
    KLB.append(br4.x[2]); DB.append(br1.fun - br4.fun)
    if (rep+1) % 10 == 0:
        say(f"    rep {rep+1}/40 done ({(time.time()-T00)/60:.1f} min)")
K1B, KDB, KTB, KLB, DB = map(np.array, (K1B, KDB, KTB, KLB, DB))
q1 = np.percentile(K1B, [16, 50, 84])
qd = np.percentile(KTB-KDB, [16, 50, 84])
qD = np.percentile(DB, [16, 50, 84])
n62 = int((DB >= 6.2).sum()); n118 = int((DB >= 11.8).sum())
say(f"  F1 kappa-hat 16/50/84 = {q1[0]:.3f}/{q1[1]:.3f}/{q1[2]:.3f}")
say(f"  (kappa_t - kappa_d) 16/50/84 = {qd[0]:.3f}/{qd[1]:.3f}/{qd[2]:.3f}"
    f"   width = {qd[2]-qd[0]:.3f}")
say(f"  D(F1-F4) 16/50/84 = {qD[0]:.2f}/{qD[1]:.2f}/{qD[2]:.2f}; "
    f"reps >= 6.2: {n62}/40; >= 11.8: {n118}/40")
POWER = (qd[2]-qd[0]) > 1.0

# ---------------- G6 vertical co-read ----------------
say("")
say("G6 VERTICAL CO-READ (5M apparatus, kappa family; point grade):")
def m2kv(t, mode, dml, dv, sv, w_g, bsp, ble):
    nk = NK[mode]
    tk, th = t[:nk], t[nk:]
    if mode != 'fix1' and not all(0.02 < k < 3.0 for k in tk): return 1e12
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    ksp, kle = kaps_from(mode, tk, bsp, ble)
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*(1.0 + ksp*occ_be(gN/a0))
    se2 = sig2 + s_int*s_int
    r = LGOBS - np.log10(gm_) - dv[gidx]
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    rl = l_gobs[lmask] - (lg + np.log10(1.0 + kle*occ_be(10**lg/a0)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    out += np.sum(w_g*dv*dv/(sv*sv))
    return out

def fit_vert(mode, sv, starts, tol=0.05, max_rounds=12):
    w_g = np.ones(NGal)
    nk = NK[mode]
    dml = np.zeros(NGal)
    dv = np.zeros(NGal)
    best, prev = None, None
    for rd in range(max_rounds):
        ss = ([list(best.x)] if best is not None else []) + \
             (starts if rd == 0 else [])
        bb = None
        for t0 in ss:
            b = minimize(lambda t: m2kv(t, mode, dml, dv, sv, w_g, BSP, BLE),
                         t0, method='Nelder-Mead',
                         options=dict(maxiter=4000+1500*nk, xatol=1e-6,
                                      fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        tk, (la0, f, s_int, dlt) = best.x[:nk], best.x[nk:]
        ksp, _ = kaps_from(mode, tk, BSP, BLE)
        se2c = s_int*s_int
        for _ in range(3):
            fac = f*np.exp(dml[gidx])
            gN = g_gas + fac*g_dsk + g_bul
            r0 = LGOBS - np.log10(gN*(1.0 + ksp*occ_be(gN/10**la0)))
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                w = 1.0/(sig2[mm] + se2c)
                dv[gi2] = np.sum(w*r0[mm])/(np.sum(w) + 1.0/sv[gi2]**2)
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                kk = ksp if np.isscalar(ksp) else ksp[mm]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = (LGOBS[mm]
                          - np.log10(gN2*(1.0 + kk*occ_be(gN2/10**la0)))
                          - dv[gi2])
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2kv(best.x, mode, dml, dv, sv, w_g, BSP, BLE)
        if prev is not None and abs(prev - cur) < tol: break
        prev = cur
    b = minimize(lambda t: m2kv(t, mode, dml, dv, sv, w_g, BSP, BLE),
                 list(best.x), method='Nelder-Mead',
                 options=dict(maxiter=4000+1500*nk, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best

bv3 = fit_vert('fix1', SIGV, [[math.log10(A0_FID), 1.0, 0.08, 0.0]],
               max_rounds=15)
okv = abs(bv3.fun - (-12152.49)) < 1.5
say(f"  F3v (kappa=1) regression: {bv3.fun:.2f} (5M dv-ON BE -12152.49) "
    f"-> {'PASS' if okv else 'FAIL'}")
bv1 = fit_vert('one', SIGV, [[kap1] + list(b1.x[1:]),
                             [1.0, math.log10(A0_FID), 1.0, 0.08, 0.0]])
say(f"  F1v: kappa = {bv1.x[0]:.3f}, a0 = {10**bv1.x[1]:.3e}, obj = "
    f"{bv1.fun:.2f} (D vs F3v {bv1.fun-bv3.fun:+.2f})")
bv4 = fit_vert('three', SIGV, [[bv1.x[0]]*3 + list(bv1.x[1:]),
                               [1.0, 1.5, 1.0, math.log10(A0_FID), 1.0,
                                0.08, 0.0]])
say(f"  F4v: kappa_d/t/l = {bv4.x[0]:.3f}/{bv4.x[1]:.3f}/{bv4.x[2]:.3f}, "
    f"D(F1v-F4v) = {bv1.fun-bv4.fun:+.2f}  [lean-grade co-read]")
G6 = okv

# ---------------- letter ----------------
say("")
say("="*78)
gates_ok = G0 and G2 and G3 and G6 and ok1
blocked = not ok2b
if not (G0 and G2):
    letter = 'STOP'
elif not G3:
    letter = 'K-REPORT-ONLY (G3 fail)'
else:
    fired_split = (D_split >= 11.8) and (qD[1] >= 11.8) and (n62 >= 32)
    if fired_split:
        letter = 'K-SPLIT'
    elif D_split >= 6.2:
        letter = 'K-TENSE'
    elif q1[0] <= 1.0 <= q1[2]:
        letter = 'K-ONE'
    else:
        letter = 'K-ONE-OFF'
    if blocked and letter in ('K-ONE', 'K-ONE-OFF', 'K-SPLIT'):
        letter = 'K-TENSE (G2b blocked)'
if POWER:
    letter += ' +POWER-LIMITED'
say(f"TOKENS: G0:{'P' if G0 else 'F'} G1:{'P' if ok1 else 'F'} "
    f"G2:{'P' if G2 else 'F'}{'' if ok2b else '(b-)'} "
    f"G3:{'P' if G3 else 'F'} G6:{'P' if G6 else 'F'}  "
    f"kappa1={kap1:.3f} boot[{q1[0]:.2f},{q1[2]:.2f}] "
    f"Dsplit={D_split:+.1f} bootD[{qD[0]:.1f},{qD[2]:.1f}]")
say(f"LETTER: {letter}")
say("")
say("MANDATORY CO-READS (print in every cell): F2 lock row kappa = "
    f"{b2.x[0]:.3f} (Planck) / {b2s.x[0]:.3f} (SH0ES); vertical F1v "
    f"kappa = {bv1.x[0]:.3f}, D(F1v-F4v) = {bv1.fun-bv4.fun:+.2f}"
    f"{'' if G6 else '  [VERTICAL ROW UNLICENSED -- regression fail]'}")
say("")
say("CREDENCE (pre-signed arc map): K-ONE + clean round -> mech 15->17; "
    "K-ONE-OFF/K-SPLIT -> HOLD 15; strike cells live in 10E only; "
    "anomaly-real 53 UNTOUCHED.")
say("="*78)
say(f"done ({(time.time()-T00)/60:.1f} min)")

with open('data/stage10d_kappa.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage10d_kappa.txt")
