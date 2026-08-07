"""
STAGE 9U (O5-AVERAGING): THE GAMMA METER -- where is the galaxy tail
exponent p, CONTINUOUSLY, and with what error?

ROUND 17 left one load-bearing question: all couplings are O(H), so the
universe has completed ~one Rabi swing (gamma = Omega*T ~ 1) -- the 1/2
in p = 1/2 + r*g/2 is a time-average the universe may not have taken.
The three pre-registered landing bands (fiducial gate g = 0.7536):
  ONE-SWING   p_one   = 1/2 + sin^2(1/2)*g/2  = 0.5866   (gamma = 1, r = 0.230)
  FLOOR EDGE  p_floor = 1/2 + 0.454*g/2       = 0.6711   (9T scan floor r = 0.454)
  FULL-AVG    p_full  = 1/2 + g/4             = 0.6884   (secular r = 1/2)
The 5G anchor p-hat = 0.65 sits BETWEEN one-swing and the floor, on a
0.05-step grid with no sigma.  The 9T pre-registered kill clause is armed:
sigma_p <= 0.02 demanding p < 0.67 kills the full-averaging resonance
reading.  This stage is the instrument that can fire or clear it.

DESIGN (no new function forms -- nu_p is world-table family; freeze-
compliant refinement, the successor 9S/9T both queued by name):
  A1 plain-hier fine profile: 5G machinery VERBATIM (fit_conv/m2h),
     p in {0.52..0.58 step .02} + {0.60..0.72 step .01} + {0.74..0.84
     step .02} (23 pts), warm-chained ascending.  Continuity anchor.
  A2 vertical-hardened fine profile (PRIMARY per the 7I freeze): 5M/6J
     machinery VERBATIM (m2hv generic-nu + fit_v with w_g/lens_obs
     threading = the 6J-validated merge), same 23-pt grid.
  B  bootstrap sigma_p on A2: 40 paired galaxy-resample reps (rng 53,
     warm-lite tol 0.5 / max_rounds 5 / sweeps 2 = the 6J spec), per-rep
     5-pt local grid {p_c +/- 0.06 step 0.03}, parabolic minimum.
     sigma_p = SD(ddof=1) of the 40 per-rep minima (bootstrap PRIMARY
     per the 4S lesson; profile-curvature sigma co-quoted).
  C  inversions (arithmetic, mpmath/sympy-gated): r = 2(p-1/2)/g per
     gate treatment g in {0.7536 fid, 0.8681 maxclust, 0.9518 noclust};
     gamma-dial under BOTH readings: instantaneous r = sin^2(gamma/2)
     -> gamma = 2 asin(sqrt r); window-average r = (1 - sin(gamma)/
     gamma)/2 inverted on (0, pi].  Error propagated from sigma_p.

GATES (bars locked at this commit, BEFORE any run):
  G9U-0a  plain regression: re-run the exact 5G chain (PGRID .5/.578/
          .65/.75/.9, same warm order); each within 0.5 of the values
          parsed from data/stage5g_tailtest.txt.  5/5 required.
  G9U-0b  vertical regression: fresh-start p=0.5 within 1.0 of the 5M
          dv-ON BE value parsed from data/stage5m_hierv.txt (-12152.49);
          fresh-start p=0.65 within 1.0 of the 5M dv-ON p065 value.
          (nu_p(y,0.65) is BIT-IDENTICAL to 5M's nu_p065: 2.0*0.65
          rounds to the same double as the literal 1.3.)
  G9U-1   sky profile minima INTERIOR on both treatments + locally
          convex (both neighbors above the minimum).
  G9U-2   injection power gate (go/no-go for location-grade language):
          two mock skies generated from the p=0.65 vertical fit's
          nuisances at p_true = 0.62 and 0.70 (SPARC noise sqrt(sig2 +
          s_int^2), lensing noise sqrt(l_sig2), rng 202), refit with a
          7-pt profile (p_true +/- 0.06 step 0.02, medium convergence
          tol 0.1 / rounds 8 / sweeps 3).  Bar: |p_rec - p_true| <=
          max(0.02, 2*sigma_curv_inj) for BOTH.  FAIL -> the verdict
          letter is locked to U-POWER regardless of bootstrap sigma.
          (Unlike 9R the sky is NOT unread -- 5G published the coarse
          read; phases A/B run either way, language degrades.)
  G9U-3   bootstrap sanity: per-rep edge-hits <= 8/40; else amendment
          A-widen fires ONCE (7-pt grid +/- 0.09 step 0.03, re-run
          edged reps), logged in output.  Percentiles reported.
  G9U-4   arithmetic regressions (sympy/mpmath): p_full/p_floor/p_one
          at fid = {0.6884, 0.6711, 0.5866} within 5e-4; r(p=0.65) =
          {0.3980, 0.3456, 0.3152} per gate within 5e-4 (the 9S rows);
          window-average r(gamma) -> 1/2 as gamma -> inf (series check)
          and equals sin^2 average identity at gamma = 2 pi.

VERDICT LETTERS (primary = A2 + bootstrap sigma_p, fiducial gate for
cell edges, all three gates disclosed in the table):
  U-POWER  G9U-2 FAIL or sigma_p > 0.02.  Location claims lean-grade.
  U-FULL   powered AND [p-hat - 2 sigma_p, p-hat + 2 sigma_p]
           intersects [0.6711, 0.6884 + 2 sigma_p]: floor satisfied,
           the 9T honest tension DISSOLVES as a grid artifact.
  U-ONE    powered AND |p-hat - 0.5866| <= 2 sigma_p: the one-swing
           prediction LANDS quantitatively.
  U-GAMMA  powered AND p-hat + 2 sigma_p < 0.6711 AND p-hat - 2 sigma_p
           > 0.5866: the 9T kill FIRES on the full-averaging reading;
           gamma-hat measured interior (both readings quoted) = the
           O5-AVERAGING answer "partially averaged".
  U-LOW    powered AND p-hat + 2 sigma_p < 0.5866: below one-swing --
           the transfer picture itself is strained.

CREDENCE MAP (pre-signed; bath-mechanism conditional currently 12;
anomaly-real 53 UNTOUCHED in every cell -- this is a mechanism
instrument):
  U-FULL  -> 12 -> 16  (tension dissolved; resonance + averaging cohere)
  U-ONE   -> 12 -> 16  (the O(H) un-averaged prediction lands)
  U-GAMMA -> HOLD 12   (full-averaging sub-reading retires per the 9T
                        clause; gamma-hat booked; successor question
                        named O5-GAMMA-WHY: derive the measured gamma)
  U-LOW   -> 12 -> 8   (strike)
  U-POWER -> HOLD 12   (sigma_p itself booked as the instrument spec
                        successors must beat)
Treatment-split rule: if |p-hat_plain - p-hat_vert| > 2 sigma_p, flag
TREATMENT-SPLIT, primary stays vertical, split booked as a finding.

Writes data/stage9u_gammameter.txt (progressive) + data/stage9u_profiles.npz.
Compute: ~4-6 h wall-clock (23+23 full hier fits, 14 injection fits,
200 warm-lite bootstrap fits).
"""
import csv, glob, math, os, re, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar

T0 = time.time()
KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25
LN10 = math.log(10)
S_ML = 0.1*LN10

OUT_TXT = 'data/stage9u_gammameter.txt'
OUT_NPZ = 'data/stage9u_profiles.npz'
LINES = []
def emit(s=""):
    LINES.append(s)
    print(s, flush=True)
    with open(OUT_TXT, 'w') as f:
        f.write("\n".join(LINES) + "\n")

# ---------- data (5M loader verbatim; superset of 5G's, same cuts) ----------
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

# ---------- the family (5G verbatim) ----------
def nu_p(y, p):
    yc = np.clip(y, 1e-14, None)
    ex = np.exp(-np.minimum(yc**p, 60.0))
    return (1.0-ex)**(-1.0/(2.0*p))

# ---------- plain-hier objective + fitter (5G verbatim) ----------
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

# ---------- vertical objective + fitter (5M body, 6J threading) ----------
def m2hv(th, nu, dml, dv, w_g, sky_obs, lens_obs):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*nu(gN/a0)
    se2 = sig2 + s_int*s_int
    r = sky_obs - np.log10(gm_) - dv[gidx]
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    rl = lens_obs[lmask] - (lg + np.log10(nu(10**lg/a0)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    out += np.sum(w_g*dv*dv/(SIGV*SIGV))
    return out

def fit_v(nu, w_g, sky_obs, lens_obs, th0=None, dml0=None, dv0=None,
          tol=0.05, max_rounds=15, sweeps=3):
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
            b = minimize(lambda t: m2hv(t, nu, dml, dv, w_g, sky_obs,
                                        lens_obs), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int, dlt = best.x
        se2c = s_int*s_int
        for _ in range(sweeps):
            fac = f*np.exp(dml[gidx])
            gN = g_gas + fac*g_dsk + g_bul
            r0 = sky_obs - np.log10(gN*nu(gN/10**la0))
            for gi2 in range(NGal):
                if w_g[gi2] == 0: dv[gi2] = 0.0; continue
                mm = GIDXS[gi2]
                w = 1.0/(sig2[mm] + se2c)
                dv[gi2] = np.sum(w*r0[mm])/(np.sum(w) + 1.0/SIGV[gi2]**2)
            for gi2 in range(NGal):
                if w_g[gi2] == 0: dml[gi2] = 0.0; continue
                mm = GIDXS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = sky_obs[mm] - np.log10(gN2*nu(gN2/10**la0)) - dv[gi2]
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2hv(best.x, nu, dml, dv, w_g, sky_obs, lens_obs)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2hv(t, nu, dml, dv, w_g, sky_obs, lens_obs),
                 list(best.x), method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml, dv

ONES = np.ones(NGal)

def parab_min(ps, vs):
    """vertex of the quadratic through the argmin and neighbors.
    returns (p_hat, sigma_curv, edge_flag)"""
    ps = np.asarray(ps, float); vs = np.asarray(vs, float)
    i = int(np.argmin(vs))
    if i == 0 or i == len(ps)-1:
        return float(ps[i]), None, True
    x0, x1, x2 = ps[i-1], ps[i], ps[i+1]
    y0, y1, y2 = vs[i-1], vs[i], vs[i+1]
    d1 = (y2-y0)/(x2-x0)
    d2 = ((y2-y1)/(x2-x1) - (y1-y0)/(x1-x0))/(0.5*(x2-x0))
    if d2 <= 0:
        return float(x1), None, False
    ph = x1 - d1/d2
    ph = min(max(ph, x0), x2)
    sig = math.sqrt(2.0/d2)
    return float(ph), float(sig), False

# ================= PHASE 0: bands + arithmetic gates =================
emit("STAGE 9U THE GAMMA METER: %d galaxies, %d points + %d lensing" %
     (kept, len(gobs), int(lmask.sum())))
emit("pre-reg bands (fiducial g = 0.7536): one-swing 0.5866 | floor "
     "0.6711 | full-avg 0.6884; 5G anchor p-hat = 0.65 (grid 0.05)")
emit("")

GATES = {'fid': 0.7536, 'maxclust': 0.8681, 'noclust': 0.9518}
R_FLOOR = 0.454
S_HALF = math.sin(0.5)**2

def p_of_r(r, g): return 0.5 + r*g/2.0
def r_of_p(p, g): return 2.0*(p-0.5)/g

def gamma_inst(r):
    if not (0.0 < r < 1.0): return float('nan')
    return 2.0*math.asin(math.sqrt(r))

def gamma_wind(r):
    # (1 - sin(g)/g)/2 = r on (0, pi]
    if not (0.0 < r < 0.5): return float('nan')
    lo, hi = 1e-9, math.pi
    for _ in range(200):
        mid = 0.5*(lo+hi)
        val = 0.5*(1.0 - math.sin(mid)/mid)
        if val < r: lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

g4_rows = []
g4_rows.append(abs(p_of_r(0.5, GATES['fid']) - 0.6884) < 5e-4)
g4_rows.append(abs(p_of_r(R_FLOOR, GATES['fid']) - 0.6711) < 5e-4)
g4_rows.append(abs(p_of_r(S_HALF, GATES['fid']) - 0.5866) < 5e-4)
g4_rows.append(abs(r_of_p(0.65, GATES['fid']) - 0.3980) < 5e-4)
g4_rows.append(abs(r_of_p(0.65, GATES['maxclust']) - 0.3456) < 5e-4)
g4_rows.append(abs(r_of_p(0.65, GATES['noclust']) - 0.3152) < 5e-4)
# window-average identities: gamma -> inf limit 1/2; at 2 pi the window
# average over exactly one full period equals sin^2 average = 1/2 - 0 (sinc(2pi)=0)
g4_rows.append(abs(0.5*(1.0 - math.sin(50.0)/50.0) - 0.5) < 0.005)
g4_rows.append(abs(0.5*(1.0 - math.sin(2*math.pi)/(2*math.pi)) - 0.5) < 1e-12)
ok4 = all(g4_rows)
emit("G9U-4 arithmetic regressions (9S/9T rows): %s  [%s]" %
     ('PASS' if ok4 else 'FAIL', ''.join('1' if x else '0' for x in g4_rows)))
emit("")

# ================= PHASE 1: G9U-0a plain 5G chain regression ==========
tgt5g = {}
with open('data/stage5g_tailtest.txt') as f:
    for m in re.finditer(r'p=(\d\.\d+):\s+(-\d+\.\d+)', f.read()):
        tgt5g[float(m.group(1))] = float(m.group(2))
emit("PHASE 1: 5G chain regression (targets: %s)" %
     {k: v for k, v in sorted(tgt5g.items())})
PGRID_5G = [0.5, 0.578, 0.65, 0.75, 0.9]
th, dm = None, None
reg_ok = []
plain_chain = {}
for p in PGRID_5G:
    nu = lambda y, p=p: nu_p(y, p)
    b, dm = fit_conv(nu, th0=th, dml0=dm)
    th = b.x
    plain_chain[p] = (b, dm.copy())
    d = b.fun - tgt5g.get(round(p, 3), np.nan)
    ok = abs(d) < 0.5
    reg_ok.append(ok)
    emit("  p=%5.3f: %10.2f  (5G %10.2f, d=%+.3f) %s  [%.1f min]" %
         (p, b.fun, tgt5g.get(round(p, 3), np.nan), d,
          'OK' if ok else 'FAIL', (time.time()-T0)/60))
ok0a = all(reg_ok)
emit("G9U-0a: %s" % ('PASS' if ok0a else 'FAIL'))
emit("")

# ================= PHASE 2: G9U-0b vertical regressions ==============
tgt5m = {}
with open('data/stage5m_hierv.txt') as f:
    txt5m = f.read()
for nm in ('BE', 'p065'):
    m = re.search(nm + r':\s+(-\d+\.\d+)', txt5m)
    if m: tgt5m[nm] = float(m.group(1))
emit("PHASE 2: 5M dv-ON regressions (targets: %s)" % tgt5m)
bv50, dml50, dv50 = fit_v(lambda y: nu_p(y, 0.5), ONES, lgobs, l_gobs)
d50 = bv50.fun - tgt5m['BE']
bv65, dml65, dv65 = fit_v(lambda y: nu_p(y, 0.65), ONES, lgobs, l_gobs)
d65 = bv65.fun - tgt5m['p065']
ok0b = abs(d50) < 1.0 and abs(d65) < 1.0
emit("  p=0.50: %10.2f (5M BE   %10.2f, d=%+.3f)" %
     (bv50.fun, tgt5m['BE'], d50))
emit("  p=0.65: %10.2f (5M p065 %10.2f, d=%+.3f)" %
     (bv65.fun, tgt5m['p065'], d65))
emit("G9U-0b: %s  [%.1f min]" %
     ('PASS' if ok0b else 'FAIL', (time.time()-T0)/60))
emit("")

# ================= PHASE 3: G9U-2 injection power gate ===============
emit("PHASE 3: injection power gate (p_true 0.62 / 0.70; rng 202)")
rng_inj = np.random.default_rng(202)
la0_65, f_65, sint_65, dlt_65 = bv65.x
a0_65 = 10**la0_65
inj_results = []
for p_true in (0.62, 0.70):
    fac = f_65*np.exp(dml65[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    model_sky = np.log10(gN*nu_p(gN/a0_65, p_true)) + dv65[gidx]
    mock_sky = model_sky + rng_inj.normal(
        0, np.sqrt(sig2 + sint_65*sint_65))
    lg = l_gbar + dlt_65
    model_lens = lg + np.log10(nu_p(10**lg/a0_65, p_true))
    mock_lens = model_lens + rng_inj.normal(0, np.sqrt(l_sig2))
    ps = [round(p_true + k*0.02, 3) for k in range(-3, 4)]
    thi, dmi, dvi = None, None, None
    vs = []
    for p in ps:
        nu = lambda y, p=p: nu_p(y, p)
        b, dmi, dvi = fit_v(nu, ONES, mock_sky, mock_lens, th0=thi,
                            dml0=dmi, dv0=dvi, tol=0.1, max_rounds=8,
                            sweeps=3)
        thi = b.x
        vs.append(b.fun)
    ph, sc, edge = parab_min(ps, vs)
    err = ph - p_true
    bar = max(0.02, 2*sc if sc else 0.02)
    ok = (not edge) and abs(err) <= bar
    inj_results.append((p_true, ph, sc, err, ok))
    emit("  p_true=%.2f: rec %.4f (sigma_curv %s, err %+.4f, bar %.3f)"
         " %s  [%.1f min]" %
         (p_true, ph, ('%.4f' % sc) if sc else 'n/a', err, bar,
          'OK' if ok else 'FAIL', (time.time()-T0)/60))
ok2 = all(r[4] for r in inj_results)
emit("G9U-2: %s%s" % ('PASS' if ok2 else 'FAIL',
     '' if ok2 else '  -> verdict locked to U-POWER'))
emit("")

# ================= PHASE 4: A1 plain fine profile ====================
FINE = ([0.52, 0.54, 0.56, 0.58] +
        [round(0.60 + 0.01*k, 2) for k in range(13)] +
        [0.74, 0.76, 0.78, 0.80, 0.82, 0.84])
emit("PHASE 4: plain-hier fine profile (%d pts)" % len(FINE))
th, dm = None, None
b58, dm58 = plain_chain[0.578]
th, dm = b58.x, dm58
plain_vals = []
for p in FINE:
    nu = lambda y, p=p: nu_p(y, p)
    b, dm = fit_conv(nu, th0=th, dml0=dm)
    th = b.x
    plain_vals.append(b.fun)
    emit("  p=%.2f: %10.2f  [%.1f min]" % (p, b.fun, (time.time()-T0)/60))
ph_p, sc_p, edge_p = parab_min(FINE, plain_vals)
emit("plain minimum: p-hat = %.4f (sigma_curv %s)%s" %
     (ph_p, ('%.4f' % sc_p) if sc_p else 'n/a',
      ' EDGE' if edge_p else ''))
emit("")

# ================= PHASE 5: A2 vertical fine profile (PRIMARY) =======
emit("PHASE 5: vertical-hardened fine profile (%d pts, PRIMARY)" %
     len(FINE))
thv, dmv, dvv = bv50.x, dml50, dv50
vert_vals = []
vert_fits = {}
for p in FINE:
    nu = lambda y, p=p: nu_p(y, p)
    b, dmv, dvv = fit_v(nu, ONES, lgobs, l_gobs, th0=thv, dml0=dmv,
                        dv0=dvv)
    thv = b.x
    vert_vals.append(b.fun)
    vert_fits[p] = (b.x.copy(), dmv.copy(), dvv.copy())
    emit("  p=%.2f: %10.2f  [%.1f min]" % (p, b.fun, (time.time()-T0)/60))
ph_v, sc_v, edge_v = parab_min(FINE, vert_vals)
ok1g = (not edge_p) and (not edge_v)
emit("vertical minimum: p-hat = %.4f (sigma_curv %s)%s" %
     (ph_v, ('%.4f' % sc_v) if sc_v else 'n/a',
      ' EDGE' if edge_v else ''))
emit("G9U-1 interior minima: %s" % ('PASS' if ok1g else 'FAIL'))
emit("")

# ================= PHASE 6: bootstrap sigma_p on A2 ==================
emit("PHASE 6: bootstrap (40 paired reps, rng 53, warm-lite)")
p_c = round(ph_v, 2)
BOOT_PS = [round(p_c + k*0.03, 2) for k in range(-2, 3)]
BOOT_PS = [min(max(p, 0.51), 0.90) for p in BOOT_PS]
emit("  per-rep grid: %s" % BOOT_PS)
warm = {}
for p in BOOT_PS:
    if p in vert_fits:
        warm[p] = vert_fits[p]
    else:
        nu = lambda y, p=p: nu_p(y, p)
        b, dmw, dvw = fit_v(nu, ONES, lgobs, l_gobs, th0=thv, dml0=dmv,
                            dv0=dvv, tol=0.1, max_rounds=8)
        warm[p] = (b.x.copy(), dmw.copy(), dvw.copy())
rng = np.random.default_rng(53)
allg = np.arange(NGal)
boot_ph = []
edge_hits = 0
rep_rows = []
for k in range(40):
    pick = rng.choice(allg, NGal, replace=True)
    w = np.zeros(NGal)
    for g_ in pick: w[g_] += 1
    lo = l_gobs + rng.normal(0, np.sqrt(l_sig2))
    vs = []
    for p in BOOT_PS:
        tw, dw, vw = warm[p]
        nu = lambda y, p=p: nu_p(y, p)
        b, _, _ = fit_v(nu, w, lgobs, lo, th0=tw, dml0=dw, dv0=vw,
                        tol=0.5, max_rounds=5, sweeps=2)
        vs.append(b.fun)
    ph, _, edge = parab_min(BOOT_PS, vs)
    boot_ph.append(ph)
    rep_rows.append((k, ph, edge))
    if edge: edge_hits += 1
    if (k+1) % 5 == 0:
        emit("  rep %d/40: running mean %.4f sd %.4f edges %d  [%.1f min]"
             % (k+1, np.mean(boot_ph),
                np.std(boot_ph, ddof=1) if k else 0.0, edge_hits,
                (time.time()-T0)/60))
boot_ph = np.array(boot_ph)
ok3 = edge_hits <= 8
if not ok3:
    emit("  AMENDMENT A-widen FIRED (%d edge hits > 8): widening to "
         "+/- 0.09 step 0.03 for edged reps" % edge_hits)
    WIDE_PS = [round(p_c + k*0.03, 2) for k in range(-3, 4)]
    WIDE_PS = [min(max(p, 0.51), 0.90) for p in WIDE_PS]
    for p in WIDE_PS:
        if p not in warm:
            nu = lambda y, p=p: nu_p(y, p)
            b, dmw, dvw = fit_v(nu, ONES, lgobs, l_gobs, th0=thv,
                                dml0=dmv, dv0=dvv, tol=0.1, max_rounds=8)
            warm[p] = (b.x.copy(), dmw.copy(), dvw.copy())
    rng = np.random.default_rng(53)
    boot2 = []
    for k in range(40):
        pick = rng.choice(allg, NGal, replace=True)
        w = np.zeros(NGal)
        for g_ in pick: w[g_] += 1
        lo = l_gobs + rng.normal(0, np.sqrt(l_sig2))
        if not rep_rows[k][2]:
            boot2.append(rep_rows[k][1])
            continue
        vs = []
        for p in WIDE_PS:
            tw, dw, vw = warm[p]
            nu = lambda y, p=p: nu_p(y, p)
            b, _, _ = fit_v(nu, w, lgobs, lo, th0=tw, dml0=dw, dv0=vw,
                            tol=0.5, max_rounds=5, sweeps=2)
            vs.append(b.fun)
        ph, _, _ = parab_min(WIDE_PS, vs)
        boot2.append(ph)
        emit("  widened rep %d: %.4f  [%.1f min]" %
             (k, ph, (time.time()-T0)/60))
    boot_ph = np.array(boot2)
sig_boot = float(np.std(boot_ph, ddof=1))
pcts = np.percentile(boot_ph, [16, 50, 84])
emit("bootstrap: p-hat SD = %.4f; percentiles 16/50/84 = %s; "
     "edge hits %d/40" % (sig_boot, pcts.round(4).tolist(), edge_hits))
emit("G9U-3: %s" % ('PASS' if ok3 else
     'AMENDED (A-widen applied, logged)'))
emit("")

# ================= PHASE 7: inversions + verdict =====================
POWERED = ok2 and sig_boot <= 0.02
p_hat = ph_v
s_p = sig_boot
emit("=" * 68)
emit("PRIMARY (vertical-hardened + bootstrap): p-hat = %.4f +/- %.4f"
     % (p_hat, s_p))
emit("  (curvature sigma %.4f; plain-hier p-hat = %.4f; split %+.4f)"
     % (sc_v if sc_v else float('nan'), ph_p, ph_p - ph_v))
split_flag = abs(ph_p - ph_v) > 2*s_p
if split_flag:
    emit("  TREATMENT-SPLIT flag: |plain - vert| > 2 sigma_p "
         "(primary stays vertical)")
emit("")
emit("gate-resolved inversion table (r = 2(p-1/2)/g; gamma dial):")
emit("  %-9s %-8s %-8s %-8s %-10s %-10s" %
     ("gate", "g", "r-hat", "sigma_r", "gam_inst", "gam_wind"))
inv_rows = {}
for gname, g in GATES.items():
    r_hat = r_of_p(p_hat, g)
    s_r = 2*s_p/g
    gi_ = gamma_inst(r_hat)
    gw_ = gamma_wind(r_hat)
    dgi = (gamma_inst(min(r_hat + 0.5*s_r, 0.999)) -
           gamma_inst(max(r_hat - 0.5*s_r, 1e-4)))/1.0 if s_r > 0 else 0
    inv_rows[gname] = (r_hat, s_r, gi_, gw_)
    emit("  %-9s %-8.4f %-8.4f %-8.4f %-10.4f %-10.4f" %
         (gname, g, r_hat, s_r, gi_, gw_))
emit("")
emit("band distances (fiducial, in sigma_p): one-swing %+.2f | floor "
     "%+.2f | full-avg %+.2f" %
     ((p_hat - 0.5866)/s_p, (p_hat - 0.6711)/s_p, (p_hat - 0.6884)/s_p))
emit("")

if not POWERED:
    letter = 'U-POWER'
elif (p_hat + 2*s_p >= 0.6711):
    letter = 'U-FULL'
elif abs(p_hat - 0.5866) <= 2*s_p:
    letter = 'U-ONE'
elif (p_hat + 2*s_p < 0.6711) and (p_hat - 2*s_p > 0.5866):
    letter = 'U-GAMMA'
elif (p_hat + 2*s_p < 0.5866):
    letter = 'U-LOW'
else:
    letter = 'U-POWER'
emit("VERDICT LETTER: %s" % letter)
cred = {'U-FULL': 'bath-mechanism conditional 12 -> 16',
        'U-ONE': 'bath-mechanism conditional 12 -> 16',
        'U-GAMMA': 'HOLD 12; full-averaging sub-reading retires per the '
                   '9T clause; gamma-hat booked; O5-GAMMA-WHY named',
        'U-LOW': 'bath-mechanism conditional 12 -> 8',
        'U-POWER': 'HOLD 12; sigma_p booked as the instrument spec'}
emit("pre-signed credence cell: %s" % cred[letter])
emit("anomaly-real 53 UNTOUCHED (pre-signed, all cells)")
emit("")
emit("gates: G9U-0a %s | G9U-0b %s | G9U-1 %s | G9U-2 %s | G9U-3 %s | "
     "G9U-4 %s" %
     tuple('PASS' if x else ('AMENDED' if i == 4 and not ok3 else 'FAIL')
           for i, x in enumerate([ok0a, ok0b, ok1g, ok2, ok3, ok4])))
emit("total wall-clock %.1f min" % ((time.time()-T0)/60))

np.savez(OUT_NPZ,
         fine_p=np.array(FINE), plain=np.array(plain_vals),
         vert=np.array(vert_vals), boot=boot_ph,
         boot_ps=np.array(BOOT_PS),
         inj=np.array([(a, b, c if c else -1, d) for a, b, c, d, _ in
                       inj_results]),
         p_hat=p_hat, sig_boot=sig_boot,
         sc_v=sc_v if sc_v else -1, ph_plain=ph_p)
print("\nsaved:", OUT_TXT, "+", OUT_NPZ)
