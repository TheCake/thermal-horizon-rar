"""
STAGE 9V (O5-AVERAGING successor): THE R-LADDER -- the direct exchange-
weight fit, with the per-galaxy measured gates INSIDE the likelihood.

9U measured the global tail exponent (p-hat = 0.6471 +/- >=0.075,
U-POWER) and demoted the 9S point-bound.  The queued successor (named
by 9S, freeze-labeled CONSISTENCY) is the two-sided r MEASUREMENT:
    p_i = 1/2 + r * g_i / 2       (one parameter r, per-galaxy gate g_i)
with g_i = s_i^2 the measured Chae-ambient gates (6I machinery), fitted
through the nu_p family (world-table member; no new function form).
This differs from 9U in information structure: each galaxy's tail is
predicted at its OWN gate (gate heterogeneity modeled, not averaged),
and the gate-treatment systematic becomes an explicit treatment axis.

FEASIBILITY VERDICT BOOKED IN-STAGE (primary read 2026-08-07, trap #6
executed): the void-TAIL channel is DATA-BLOCKED at public grade --
Pustilnik+2020 (MNRAS 491,4993; Table 1 read from the PDF today) void
dwarfs have V_max = 31.5-80.3 km/s = deep-MOND objects with NO
Newtonian arm (the tail lives at y >~ 1); VGS masses are the same
class.  The FEAS block below quantifies with SPARC dwarf analogs.
Reopen conditions: a void-crossmatched sample WITH y >~ 1 coverage
(WALLABY DR2 + environment tags; DR4-era), or the unified LV database
gaining environment columns.

R BANDS (locked): r = 0 no-exchange (== BE exactly); r = 0.2298
one-swing sin^2(1/2); r = 0.454 the 9T theorem floor; r = 0.500 full
secular averaging.  The 9U conversion expectation: sigma_r ~ 2
sigma_p / g ~ 0.17 -- the stage's own bootstrap decides the letter.

DESIGN:
  Treatments: PRIMARY = vertical-hardened + per-galaxy maxclust gates
  (the 6J 'pgmax' prescription); co-reads = global-fiducial gate
  (g = 0.7536 all galaxies) and per-galaxy noclust.  Unmatched
  galaxies take the matched median (6J pattern); lensing leg runs at
  the matched-median gate of its treatment.
  Fine r-grid 0.00..1.10 step 0.05 (23 pts, warm-chained ascending;
  headroom above the r = 1 instantaneous ceiling so edge-riding is
  detectable rather than manufactured).
  Bootstrap sigma_r on the PRIMARY: 40 paired galaxy reps (rng 53,
  6J warm-lite), per-rep 5-pt grid r_c +/- 0.15 step 0.075, parabolic
  minimum; SD(ddof=1).

GATES (bars locked at this commit, BEFORE any run):
  G9V-0  code-path identity: at the global-fid treatment, r = 0.5
         evaluated through the per-point-p path must equal the 9U
         scalar-p path at identical (th, dml, dv): |d| <= 1e-9.
  G9V-1  r = 0 fresh-start vertical fit reproduces the 5M dv-ON BE
         value (-12152.49, parsed) within 1.0 (r=0 => p_i = 1/2 = BE
         exactly, gate-independent).
  G9V-2  injection power gate (rng 202): mock skies at r_true = 0.25
         and 0.50 from the nearest fitted nuisances (primary
         treatment), 7-pt refit profile r_true +/- 0.15 step 0.05,
         medium convergence.  Bar: |r_rec - r_true| <= max(0.05,
         2*sigma_curv_inj) BOTH.  FAIL -> letter locked to V-POWER.
  G9V-3  sky minima interior + locally convex (all three treatments).
  G9V-4  arithmetic regressions: p(0.5, 0.7536) = 0.6884; p(0.2298,
         0.7536) = 0.5866; r_of_p(0.6471, fid) = 0.3904 (the 9U tie);
         within 5e-4 each.
  G9V-5  bootstrap edge rule: <= 8/40 edge hits, else pre-authorized
         A-widen fires ONCE (+/- 0.225 step 0.075, re-run edged reps);
         persisting edge pileup => sigma_r quoted as CLIPPED LOWER
         BOUND (the 9U language).

VERDICT LETTERS (primary treatment + bootstrap sigma_r):
  V-LOW       powered (sigma_r <= 0.15, G9V-2 PASS) AND r-hat +
              2 sigma_r < 0.2298: below one-swing -- the exchange
              picture itself strained.
  V-SHARP     sigma_r <= 0.07 AND the r-hat +/- 2 sigma_r interval
              excludes at least one of {0.2298, 0.500} (name which).
  V-MEASURED  sigma_r <= 0.15: the two-sided r measurement replaces
              the demoted 9S bound; r = 0 exclusion quoted at its
              achieved sigma.
  V-POWER     sigma_r > 0.15 or G9V-2 FAIL: r stays lean-grade;
              sigma_r booked as spec.
  OVERSHOOT disclosure (any letter): r-hat - 0.5 > sigma_r -> note
  the instantaneous-reading region (sin^2 can exceed 1/2); no verdict
  weight.

CREDENCE MAP (pre-signed): ONLY V-LOW moves -- bath-mechanism
conditional 12 -> 8.  Every other cell HOLD 12.  anomaly-real 53
UNTOUCHED in all cells.  r-hat +/- sigma_r becomes the operative
averaging-flank statement (superseding the 9S numeric bound) in
EVERY cell including V-POWER (lean-labeled there).

Writes data/stage9v_rladder.txt (progressive) + data/stage9v_profiles.npz.
Compute: ~5 min wall-clock at 9U-measured fit speed.
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

OUT_TXT = 'data/stage9v_rladder.txt'
OUT_NPZ = 'data/stage9v_profiles.npz'
LINES = []
def emit(s=""):
    LINES.append(s)
    print(s, flush=True)
    with open(OUT_TXT, 'w') as f:
        f.write("\n".join(LINES) + "\n")

# ---------- data (9U loader verbatim + per-galaxy V_max tracking) ----------
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
vmax_map = {}
gal_name = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 10.0, 1.0, 3.0))
    if inc < 30 or q > 2: continue
    kept += 1
    gal_name[gi] = name
    irad = math.radians(inc)
    sv = math.hypot((eD/max(D, 1e-3))/LN10,
                    2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10)
    sigv_g_map[gi] = max(sv, 0.01)
    vm = 0.0
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC; gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        vm = max(vm, Vo)
        g_gas.append(gg); g_dsk.append(gd); g_bul.append(gb)
        gobs.append(Vo*Vo/R*KPC)
        sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi)
    vmax_map[gi] = vm
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs = np.log10(gobs)
ug = np.unique(gal_id)
NGal = len(ug)
gmap = {g: i for i, g in enumerate(ug)}
gidx = np.array([gmap[g] for g in gal_id])
SIGV = np.array([sigv_g_map[g] for g in ug])
VMAX = np.array([vmax_map[g] for g in ug])
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

# ---------- gates (6J pattern, both Chae treatments) ----------
def n_amb_of(e):
    x = math.sqrt(e)
    return 1.0/(math.exp(x) - 1.0)
def s_of(e):
    n = n_amb_of(e)
    return n/(1.0 + n)

G_FID = s_of(0.02)**2

chae_mx, chae_nc = {}, {}
with open('data/chae2021_table3.csv') as f:
    for row in csv.DictReader(l for l in f if not l.startswith('#')):
        chae_mx[row['galaxy']] = 10.0**float(row['log_eN_maxclust'])
        chae_nc[row['galaxy']] = 10.0**float(row['log_eN_noclust'])

def gate_arrays(cat):
    raw = np.full(NGal, np.nan)
    for i in range(NGal):
        nm = gal_name[ug[i]]
        if nm in cat:
            raw[i] = s_of(cat[nm])**2
    matched = np.isfinite(raw)
    med = float(np.median(raw[matched]))
    return np.where(matched, raw, med), med, int(matched.sum())

G_MX, G_MX_MED, n_mx = gate_arrays(chae_mx)
G_NC, G_NC_MED, n_nc = gate_arrays(chae_nc)
TREATMENTS = {
    'pgmax': (G_MX, G_MX_MED),
    'fid': (np.full(NGal, G_FID), G_FID),
    'noclust': (G_NC, G_NC_MED),
}

# ---------- the family (9U verbatim) ----------
def nu_p(y, p):
    yc = np.clip(y, 1e-14, None)
    ex = np.exp(-np.minimum(yc**p, 60.0))
    return (1.0-ex)**(-1.0/(2.0*p))

# ---------- objective + fitter (9U fit_v with per-point p) ----------
def m2hvr(th, P_PTS, p_lens, dml, dv, w_g, sky_obs, lens_obs):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*nu_p(gN/a0, P_PTS)
    se2 = sig2 + s_int*s_int
    r = sky_obs - np.log10(gm_) - dv[gidx]
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    rl = lens_obs[lmask] - (lg + np.log10(nu_p(10**lg/a0, p_lens)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    out += np.sum(w_g*dv*dv/(SIGV*SIGV))
    return out

def fit_vr(P_PTS, p_lens, w_g, sky_obs, lens_obs, th0=None, dml0=None,
           dv0=None, tol=0.05, max_rounds=15, sweeps=3):
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
            b = minimize(lambda t: m2hvr(t, P_PTS, p_lens, dml, dv, w_g,
                                         sky_obs, lens_obs), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int, dlt = best.x
        se2c = s_int*s_int
        for _ in range(sweeps):
            fac = f*np.exp(dml[gidx])
            gN = g_gas + fac*g_dsk + g_bul
            r0 = sky_obs - np.log10(gN*nu_p(gN/10**la0, P_PTS))
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
                    rr = (sky_obs[mm]
                          - np.log10(gN2*nu_p(gN2/10**la0, P_PTS[mm]))
                          - dv[gi2])
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2hvr(best.x, P_PTS, p_lens, dml, dv, w_g, sky_obs, lens_obs)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2hvr(t, P_PTS, p_lens, dml, dv, w_g, sky_obs,
                                 lens_obs), list(best.x),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml, dv

# scalar-p objective (9U code path) for the G9V-0 identity
def m2hv_scalar(th, nu, dml, dv, w_g, sky_obs, lens_obs):
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

ONES = np.ones(NGal)

def p_pts_of(r, gates):
    gv, gl = gates
    return 0.5 + r*gv[gidx]/2.0, 0.5 + r*gl/2.0

def parab_min(ps, vs):
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
    rh = x1 - d1/d2
    rh = min(max(rh, x0), x2)
    return float(rh), float(math.sqrt(2.0/d2)), False

# ================= PHASE 0: FEAS block + arithmetic gates ============
emit("STAGE 9V THE R-LADDER: %d galaxies, %d points + %d lensing" %
     (kept, len(gobs), int(lmask.sum())))
emit("gates: fid g = %.4f; maxclust matched %d (median g %.4f); "
     "noclust matched %d (median g %.4f)" %
     (G_FID, n_mx, G_MX_MED, n_nc, G_NC_MED))
emit("r bands: 0 = BE | 0.2298 one-swing | 0.454 floor | 0.500 full")
emit("")
emit("FEAS (void-tail channel; Pustilnik+2020 Table 1 primary-read "
     "2026-08-07):")
LC = [("KK246", 6.85, 42.0), ("UGC4115", 7.73, 56.5),
      ("J0926+3343", 10.6, 31.5), ("UGC5288", 11.4, 72.4),
      ("UGC4148", 13.5, 63.9), ("J0630+23", 22.9, 80.2),
      ("J0626+24", 23.2, 80.3), ("J0929+1155", 24.3, 62.3)]
emit("  Lynx-Cancer 8: V_max = %s km/s" %
     sorted(v for _, _, v in LC))
gN0 = g_gas + g_dsk + g_bul
ymax_gal = np.array([np.max(gN0[GIDXS[i]])/A0_FID for i in range(NGal)])
dw = VMAX <= 85.0
emit("  SPARC analogs (V_max <= 85 km/s, N = %d): per-galaxy max "
     "y_bar percentiles 16/50/84 = %s -> NO tail coverage (tail "
     "needs y >~ 1)" %
     (int(dw.sum()), np.percentile(ymax_gal[dw], [16, 50, 84]).round(3)
      .tolist()))
emit("  VERDICT (feasibility): the void-TAIL r-meter is DATA-BLOCKED "
     "at public grade; reopen = void-crossmatched y>~1 sample "
     "(WALLABY DR2 + env tags / DR4-era) or LV-database env columns")
emit("")
S_HALF = math.sin(0.5)**2
g4 = []
g4.append(abs((0.5 + 0.5*0.7536/2) - 0.6884) < 5e-4)
g4.append(abs((0.5 + S_HALF*0.7536/2) - 0.5866) < 5e-4)
g4.append(abs(2*(0.6471-0.5)/0.7536 - 0.3904) < 5e-4)
ok4 = all(g4)
emit("G9V-4 arithmetic regressions: %s  [%s]" %
     ('PASS' if ok4 else 'FAIL', ''.join('1' if x else '0' for x in g4)))
emit("")

# ================= PHASE 1: G9V-1 r=0 == BE regression ==============
tgt_be = None
with open('data/stage5m_hierv.txt') as f:
    m = re.search(r'BE:\s+(-\d+\.\d+)', f.read())
    if m: tgt_be = float(m.group(1))
P0, pl0 = p_pts_of(0.0, TREATMENTS['pgmax'])
b0, dml0, dv0 = fit_vr(P0, pl0, ONES, lgobs, l_gobs)
d0 = b0.fun - tgt_be
ok1 = abs(d0) < 1.0
emit("PHASE 1  G9V-1 r=0 == BE: %10.2f (5M %10.2f, d=%+.3f) -> %s "
     "[%.1f min]" % (b0.fun, tgt_be, d0, 'PASS' if ok1 else 'FAIL',
                     (time.time()-T0)/60))
emit("")

# ================= PHASE 2: sky profiles, three treatments ===========
RGRID = [round(0.05*k, 2) for k in range(23)]  # 0.00 .. 1.10
prof = {}
fits = {}
for tname in ('pgmax', 'fid', 'noclust'):
    gates = TREATMENTS[tname]
    emit("PHASE 2  %s profile (%d pts):" % (tname, len(RGRID)))
    th, dm, dvv = b0.x, dml0, dv0
    vals = []
    for rr in RGRID:
        PP, pl = p_pts_of(rr, gates)
        b, dm, dvv = fit_vr(PP, pl, ONES, lgobs, l_gobs, th0=th,
                            dml0=dm, dv0=dvv)
        th = b.x
        vals.append(b.fun)
        fits[(tname, rr)] = (b.x.copy(), dm.copy(), dvv.copy())
    prof[tname] = np.array(vals)
    rh, sc, edge = parab_min(RGRID, vals)
    emit("  min: r-hat = %.4f (sigma_curv %s)%s; d(-2lnL) vs r=0: "
         "%+.2f  [%.1f min]" %
         (rh, ('%.4f' % sc) if sc else 'n/a', ' EDGE' if edge else '',
          vals[int(np.argmin(vals))] - vals[0], (time.time()-T0)/60))
    prof[tname + '_min'] = (rh, sc, edge)
emit("")
ok3 = all(not prof[t + '_min'][2] for t in ('pgmax', 'fid', 'noclust'))
emit("G9V-3 interior minima (3 treatments): %s" %
     ('PASS' if ok3 else 'FAIL'))

# G9V-0 identity at fid r=0.5 (evaluation-level, lnL grade)
thI, dmI, dvI = fits[('fid', 0.5)]
PPI, plI = p_pts_of(0.5, TREATMENTS['fid'])
vA = m2hvr(thI, PPI, plI, dmI, dvI, ONES, lgobs, l_gobs)
p_scalar = 0.5 + 0.5*G_FID/2.0
vB = m2hv_scalar(thI, lambda y: nu_p(y, p_scalar), dmI, dvI, ONES,
                 lgobs, l_gobs)
ok0 = abs(vA - vB) <= 1e-9
emit("G9V-0 code-path identity (fid, r=0.5): |d| = %.2e -> %s" %
     (abs(vA - vB), 'PASS' if ok0 else 'FAIL'))
emit("")

# ================= PHASE 3: G9V-2 injection power gate ===============
emit("PHASE 3  injections (r_true 0.25 / 0.50; rng 202; primary "
     "treatment):")
rng_inj = np.random.default_rng(202)
gates_p = TREATMENTS['pgmax']
inj_results = []
for r_true in (0.25, 0.50):
    thT, dmT, dvT = fits[('pgmax', r_true)]
    la0T, fT, sintT, dltT = thT
    a0T = 10**la0T
    PPt, plt_ = p_pts_of(r_true, gates_p)
    fac = fT*np.exp(dmT[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    mock_sky = (np.log10(gN*nu_p(gN/a0T, PPt)) + dvT[gidx]
                + rng_inj.normal(0, np.sqrt(sig2 + sintT*sintT)))
    lg = l_gbar + dltT
    mock_lens = (lg + np.log10(nu_p(10**lg/a0T, plt_))
                 + rng_inj.normal(0, np.sqrt(l_sig2)))
    rs = [round(r_true + 0.05*k, 2) for k in range(-3, 4)]
    thi, dmi, dvi = None, None, None
    vs = []
    for rr in rs:
        PP, pl = p_pts_of(rr, gates_p)
        b, dmi, dvi = fit_vr(PP, pl, ONES, mock_sky, mock_lens, th0=thi,
                             dml0=dmi, dv0=dvi, tol=0.1, max_rounds=8,
                             sweeps=3)
        thi = b.x
        vs.append(b.fun)
    rh, sc, edge = parab_min(rs, vs)
    err = rh - r_true
    bar = max(0.05, 2*sc if sc else 0.05)
    ok = (not edge) and abs(err) <= bar
    inj_results.append((r_true, rh, sc if sc else -1, err, ok))
    emit("  r_true=%.2f: rec %.4f (sigma_curv %s, err %+.4f, bar %.3f)"
         " %s  [%.1f min]" %
         (r_true, rh, ('%.4f' % sc) if sc else 'n/a', err, bar,
          'OK' if ok else 'FAIL', (time.time()-T0)/60))
ok2 = all(x[4] for x in inj_results)
emit("G9V-2: %s%s" % ('PASS' if ok2 else 'FAIL',
     '' if ok2 else '  -> letter locked to V-POWER'))
emit("")

# ================= PHASE 4: bootstrap sigma_r (primary) ==============
emit("PHASE 4  bootstrap (40 paired reps, rng 53, warm-lite, pgmax):")
rh_p, sc_p, _ = prof['pgmax_min']
r_c = round(rh_p*20)/20.0
BOOT_RS = [round(r_c + k*0.075, 3) for k in range(-2, 3)]
BOOT_RS = [min(max(x, 0.0), 1.10) for x in BOOT_RS]
emit("  per-rep grid: %s" % BOOT_RS)
warm = {}
for rr in BOOT_RS:
    key = round(rr, 2)
    if (('pgmax', key) in fits) and abs(key - rr) < 1e-9:
        warm[rr] = fits[('pgmax', key)]
    else:
        PP, pl = p_pts_of(rr, gates_p)
        thn, dmn, dvn = fits[('pgmax', r_c)]
        b, dmw, dvw = fit_vr(PP, pl, ONES, lgobs, l_gobs, th0=thn,
                             dml0=dmn, dv0=dvn, tol=0.1, max_rounds=8)
        warm[rr] = (b.x.copy(), dmw.copy(), dvw.copy())
rng = np.random.default_rng(53)
allg = np.arange(NGal)
boot_rh = []
rep_rows = []
edge_hits = 0
for k in range(40):
    pick = rng.choice(allg, NGal, replace=True)
    w = np.zeros(NGal)
    for g_ in pick: w[g_] += 1
    lo = l_gobs + rng.normal(0, np.sqrt(l_sig2))
    vs = []
    for rr in BOOT_RS:
        tw, dw_, vw = warm[rr]
        PP, pl = p_pts_of(rr, gates_p)
        b, _, _ = fit_vr(PP, pl, w, lgobs, lo, th0=tw, dml0=dw_, dv0=vw,
                         tol=0.5, max_rounds=5, sweeps=2)
        vs.append(b.fun)
    rh, _, edge = parab_min(BOOT_RS, vs)
    boot_rh.append(rh)
    rep_rows.append((k, rh, edge))
    if edge: edge_hits += 1
    if (k+1) % 10 == 0:
        emit("  rep %d/40: mean %.4f sd %.4f edges %d  [%.1f min]" %
             (k+1, np.mean(boot_rh),
              np.std(boot_rh, ddof=1) if k else 0.0, edge_hits,
              (time.time()-T0)/60))
ok5 = edge_hits <= 8
clipped = False
if not ok5:
    emit("  AMENDMENT A-widen FIRED (%d > 8): +/- 0.225 step 0.075" %
         edge_hits)
    WIDE = [round(r_c + k*0.075, 3) for k in range(-3, 4)]
    WIDE = [min(max(x, 0.0), 1.10) for x in WIDE]
    for rr in WIDE:
        if rr not in warm:
            PP, pl = p_pts_of(rr, gates_p)
            thn, dmn, dvn = fits[('pgmax', r_c)]
            b, dmw, dvw = fit_vr(PP, pl, ONES, lgobs, l_gobs, th0=thn,
                                 dml0=dmn, dv0=dvn, tol=0.1,
                                 max_rounds=8)
            warm[rr] = (b.x.copy(), dmw.copy(), dvw.copy())
    rng = np.random.default_rng(53)
    boot2 = []
    still_edge = 0
    for k in range(40):
        pick = rng.choice(allg, NGal, replace=True)
        w = np.zeros(NGal)
        for g_ in pick: w[g_] += 1
        lo = l_gobs + rng.normal(0, np.sqrt(l_sig2))
        if not rep_rows[k][2]:
            boot2.append(rep_rows[k][1])
            continue
        vs = []
        for rr in WIDE:
            tw, dw_, vw = warm[rr]
            PP, pl = p_pts_of(rr, gates_p)
            b, _, _ = fit_vr(PP, pl, w, lgobs, lo, th0=tw, dml0=dw_,
                             dv0=vw, tol=0.5, max_rounds=5, sweeps=2)
            vs.append(b.fun)
        rh, _, edge2 = parab_min(WIDE, vs)
        boot2.append(rh)
        if edge2: still_edge += 1
        emit("  widened rep %d: %.4f%s  [%.1f min]" %
             (k, rh, ' EDGE' if edge2 else '', (time.time()-T0)/60))
    boot_rh = boot2
    clipped = still_edge > 0
    emit("  post-widen persistent edges: %d -> sigma_r %s" %
         (still_edge, 'CLIPPED LOWER BOUND' if clipped else 'clean'))
boot_rh = np.array(boot_rh)
sig_r = float(np.std(boot_rh, ddof=1))
pcts = np.percentile(boot_rh, [16, 50, 84])
emit("bootstrap: sigma_r = %.4f%s; percentiles 16/50/84 = %s; "
     "edge hits %d/40" %
     (sig_r, ' (clipped lower bound)' if clipped else '',
      pcts.round(4).tolist(), edge_hits))
emit("G9V-5: %s" % ('PASS' if ok5 else 'AMENDED (A-widen applied)'))
emit("")

# ================= PHASE 5: verdict ==================================
r_hat = rh_p
emit("=" * 68)
emit("PRIMARY (vertical + pgmax + bootstrap): r-hat = %.4f +/- %.4f%s"
     % (r_hat, sig_r, ' (clipped)' if clipped else ''))
emit("  curvature sigma %.4f; co-reads: fid r-hat = %.4f, noclust "
     "r-hat = %.4f" % (sc_p if sc_p else float('nan'),
                       prof['fid_min'][0], prof['noclust_min'][0]))
emit("  d(-2lnL) r=0 -> min: pgmax %+.2f | fid %+.2f | noclust %+.2f" %
     (prof['pgmax'][int(np.argmin(prof['pgmax']))] - prof['pgmax'][0],
      prof['fid'][int(np.argmin(prof['fid']))] - prof['fid'][0],
      prof['noclust'][int(np.argmin(prof['noclust']))]
      - prof['noclust'][0]))
emit("  band distances (sigma_r): one-swing %+.2f | floor %+.2f | "
     "full %+.2f | r=0 %+.2f" %
     ((r_hat-0.2298)/sig_r, (r_hat-0.454)/sig_r, (r_hat-0.500)/sig_r,
      r_hat/sig_r))
emit("")
POWERED = ok2 and sig_r <= 0.15
if not POWERED:
    letter = 'V-POWER'
elif r_hat + 2*sig_r < 0.2298:
    letter = 'V-LOW'
elif sig_r <= 0.07 and ((r_hat - 2*sig_r > 0.2298 or
                         r_hat + 2*sig_r < 0.2298) or
                        (r_hat - 2*sig_r > 0.500 or
                         r_hat + 2*sig_r < 0.500)):
    letter = 'V-SHARP'
else:
    letter = 'V-MEASURED'
overshoot = (r_hat - 0.5) > sig_r
emit("VERDICT LETTER: %s%s" %
     (letter, '  [OVERSHOOT disclosure: r-hat > 1/2 by > 1 sigma = '
      'instantaneous-reading region]' if overshoot else ''))
cred = {'V-LOW': 'bath-mechanism conditional 12 -> 8',
        'V-SHARP': 'HOLD 12 (excluded band named in NOTES)',
        'V-MEASURED': 'HOLD 12; r-hat +/- sigma_r supersedes the 9S '
                      'numeric bound as the operative averaging-flank '
                      'statement',
        'V-POWER': 'HOLD 12; r lean-grade; sigma_r booked as spec'}
emit("pre-signed credence cell: %s" % cred[letter])
emit("anomaly-real 53 UNTOUCHED (pre-signed, all cells)")
emit("")
emit("gates: G9V-0 %s | G9V-1 %s | G9V-2 %s | G9V-3 %s | G9V-4 %s | "
     "G9V-5 %s" %
     ('PASS' if ok0 else 'FAIL', 'PASS' if ok1 else 'FAIL',
      'PASS' if ok2 else 'FAIL', 'PASS' if ok3 else 'FAIL',
      'PASS' if ok4 else 'FAIL', 'PASS' if ok5 else 'AMENDED'))
emit("total wall-clock %.1f min" % ((time.time()-T0)/60))

np.savez(OUT_NPZ, rgrid=np.array(RGRID),
         pgmax=prof['pgmax'], fid=prof['fid'], noclust=prof['noclust'],
         boot=boot_rh, boot_rs=np.array(BOOT_RS),
         inj=np.array([x[:4] for x in inj_results]),
         r_hat=r_hat, sig_r=sig_r)
print("\nsaved:", OUT_TXT, "+", OUT_NPZ)
