"""
STAGE 10I (2026-08-09): THE ISOLATION LADDER -- the 9V void-channel
reopen at Local-Volume grade (O5-AVERAGING / dispersive r = 1/2 axis).

QUESTION.  The registered r-meter (PREDICTIONS.md P1 annotations 9S/9U/
10B) reads the tail exponent at the ambient gate: p_i = 1/2 + r*g_i/2.
9V measured r-hat = 0.3365 +/- 0.1869 (V-POWER, clipped) with Chae
clustering gates (range s^2 ~ 0.6-0.8).  The 9V pre-registered REOPEN
clause: "the unified LV database gaining environment columns."  The
Karachentsev Updated Nearby Galaxy Catalog (UNGC, Karachentsev, Makarov
& Kaisina 2013, AJ 145, 101; VizieR J/AJ/145/101, fetched 2026-08-09 by
calcs/fetch_flynn_corpus.py -> data/karachentsev_ungc.csv, 869 rows)
carries those columns for the Local Volume: tidal indices Ti1/Ti5, main
disturber, K luminosities, distances, coordinates.  Cross-matching
SPARC x UNGC extends the gate axis toward the open-gate (isolated) end
AND the closed-gate (group) end -- the first isolation-based ambient
gates in the program.  Two readouts:
  R1 (the r-ladder, 9V machinery verbatim): r-hat with the extended
      gate range; the registered dispersive kill (r < 1/2 at > 2 sigma)
      evaluates ONLY if the instrument is powered.
  R2 (the ordering, 6Y-P2 out-of-sample axis): Dp = p_iso - p_nonisoLV
      stratum contrast (scalar-p fits per stratum), with its honest
      power spec.  Under r = 1/2 the predicted contrast is small
      (+0.02-ish at LV isolation grade -- printed by the stage); the
      power bar below decides whether ANY ordering claim is quotable.

DATA.  SPARC (data/sparc/, our standing files, 9V loader verbatim) +
Mistele lensing leg (9V verbatim) + Chae Table 3 gates (9V verbatim) +
UNGC (data/karachentsev_ungc.csv).  NO new sky data enters the
likelihood; the UNGC axis modifies per-galaxy GATES only.

CONSTRUCTION (all constants locked here, BEFORE any census number is
computed -- the census bars below were set blind):
  Name match: canon(name) = uppercase, strip space/dash/underscore,
    strip leading zeros in the numeric tail of NGC/UGC/IC/DDO/ESO/PGC
    prefixes; try UNGC Name, SimbadName, NEDname.  Distance gate:
    |D_SPARC - D_UNGC| <= max(2.0 Mpc, 0.35*D_UNGC).
  Isolation strata: PRIMARY iso = Ti1 < 0 (Karachentsev's isolated
    class); co-reads Ti1 < -1.0 (deep-iso) and Ti1 < +0.5.
  Computed ambient (the gate dial): e_N,i = [ sum_{j != i, D_ij <= 3
    Mpc} G * (10^KLum_j) M_sun / D_ij^2 ] / A0_FID + FLOOR, with
    M/L_K = 1.0 (the UNGC Theta convention), 3D positions from
    (RAJ2000, DEJ2000, Dist), FLOOR in {0.005, 0.01, 0.02} (PRIMARY
    0.01; the large-scale-structure term the LV sum cannot see;
    additive, disclosed).  Gate g_i = s_of(e_N,i)^2, s_of = the 6E/9V
    machinery verbatim.  ISO treatment arrays: UNGC-matched galaxies
    take UNGC gates; unmatched keep the 9V pgmax arrays.
  Transition-crosser: >= 3 points with y_bar = gN0/A0_FID >= 0.8
    (gN0 at the loader's fiducial Y_disk, 9V verbatim).

GATES (bars locked at this commit; STOP = later phases do not run):
  G10I-0  CENSUS (STOP): N_match >= 20 AND N_iso(Ti1<0) >= 6 AND
          N_iso,transition >= 3.  Else letter I-FEASIBILITY-LIMITED.
  G10I-1  READER IDENTITY (lnL grade, the 8P/8Q rule): evaluate the
          10I code path on the 9V pgmax treatment over the archived
          23-point r-grid; max |d(-2lnL)| vs data/stage9v_profiles.npz
          pgmax <= 0.5 per point (warm-chained refits; the archive was
          produced by the same fitter class).  ALSO the 9V G9V-1
          anchor: r = 0 fresh vertical fit reproduces -12152.49
          within 1.0.  FAIL -> STOP (wiring, not physics).
  G10I-2  GATE CONSTRUCTION: (a) Theta re-derivation: regress our
          max_j log10(M_j/D_ij^3) against published Ti1 on the matched
          set; |slope - 1| <= 0.25 (offset printed, recorded not
          barred).  (b) UNGC-computed e_N vs Chae e_N (maxclust) on
          the overlap: Spearman rho > 0.  (c) iso-stratum gate median
          > non-iso LV gate median (the axis premise).  (a) or (b)
          FAIL -> STOP.  (c) FAIL -> ordering readout disabled
          (disclosed); ladder may still run.
  G10I-3  INJECTION POWER, SKY-BLIND ON THE NEW AXIS (the firewall
          improvement over 9V's own order): mock skies at r_true =
          0.34 and 0.50 from the G10I-1 regenerated pgmax nuisances;
          REFIT with the ISO treatment (7-pt profile r_true +/- 0.15
          step 0.05).  BAR (9V G9V-2 verbatim): |r_rec - r_true| <=
          max(0.05, 2*sigma_curv) BOTH.  FAIL -> letter I-POWER; the
          iso-treatment SKY profile is NEVER computed.  Dp power leg
          (same mocks): stratum-contrast recovery sigma via 12-rep
          paired galaxy bootstrap on the r_true = 0.50 mock; the
          ordering instrument is POWERED only if sigma(Dp) <= 0.011
          (2 sigma for the r = 1/2 LV-grade prediction).  Dp leg
          failing does NOT block the ladder (disclosed co-read).
  G10I-4  ARITHMETIC (9V G9V-4 verbatim): p(0.5,0.7536) = 0.6884;
          p(0.2298,0.7536) = 0.5866; r_of_p(0.6471,fid) = 0.3904;
          within 5e-4.
  G10I-5  BOOTSTRAP EDGE RULE (9V G9V-5 verbatim): 40 paired reps,
          rng 53; <= 8/40 edge hits, else pre-authorized A-widen fires
          ONCE (+/- 0.225 step 0.075); persisting pileup -> sigma_r
          quoted CLIPPED LOWER BOUND.
  G10I-6  sky minima interior + locally convex (iso treatment).

VERDICT LETTERS:
  I-FEASIBILITY-LIMITED  G10I-0 STOP or G10I-2(a/b) STOP.
  I-POWER      G10I-3 recovery FAIL (either r_true) or bootstrap
               sigma_r > 0.15: census + gate calibration + power spec
               ship; no sky iso read (recovery fail) or lean-grade
               only (sigma fail).
  I-MEASURED   powered, sigma_r <= 0.15: r-hat_iso +/- sigma_r
               CO-QUOTED with 9V (supersedes only if sigma improves).
  I-SHARP      sigma_r <= 0.07 AND the +/- 2 sigma interval excludes
               0.2298 or 0.500 (name which).
  I-KILL-FIRED powered per G10I-3 AND r-hat + 2*sigma_r < 0.500: the
               REGISTERED dispersive kill (PREDICTIONS 10B annotation)
               fires -> PREDICTIONS row flips at the verdict commit.
  I-ORDERED / I-ANTI  Dp/sigma(Dp) >= 2 with the Dp leg POWERED
               (else Dp is a disclosed co-read, no letter weight).

CREDENCE MAP (pre-signed; the ONLY movable cells):
  I-FEASIBILITY-LIMITED / I-POWER / I-MEASURED unresolved: HOLD
    mech-conditional 8, HOLD anomaly-real 53.
  I-KILL-FIRED (round-affirmed): mech 8 -> 4 (the dispersive
    realization struck by its own registered test).
  I-SHARP excluding 0.500 from below = the kill cell (8 -> 4).
  I-SHARP excluding 0.2298 from above: mech 8 -> 10.
  I-ORDERED (powered, 2 sigma): mech 8 -> 10; AND p_iso inside
    [0.727, 0.750] with sigma_p <= 0.02: mech 8 -> 12 (cap).
  I-ANTI (powered, 2 sigma): mech 8 -> 5 (the 6Y-P2 ordering failing
    its first powered out-of-sample galaxy test).
  anomaly-real 53 UNTOUCHED in ALL cells (family-internal instrument;
  no Newton contest).
Bankable regardless of letter: the isolated-transition-crosser census
(the number DR4-era/void-survey designers need), the UNGC-vs-Chae
ambient cross-validation (first two-meter calibration of the gate
program), the 9V reopen clause formally executed.

FIREWALL.  The iso-gate likelihood is evaluated ONLY inside G10I-3
mock refits until G10I-3 passes; the iso SKY profile is computed only
after.  G10I-1 re-evaluates the ARCHIVED 9V treatment (reproduction of
a published fit, not a new read).  Amendments, if any, are logged in
this docstring pre-quote with failed runs preserved.

COMPUTE ~25-45 min (9V fit speed; fits are warm-chained).
Writes data/stage10i_isoladder.txt (progressive) +
data/stage10i_profiles.npz.
"""
import csv, glob, math, os, re, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import spearmanr

T0 = time.time()
KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25
LN10 = math.log(10)
S_ML = 0.1*LN10

G_SI = 6.674e-11
MSUN = 1.989e30
MPC_M = 3.0857e22
ML_K = 1.0
NEIGH_MPC = 3.0
FLOORS = (0.005, 0.01, 0.02)
FLOOR_PRIMARY = 0.01
TI_PRIMARY = 0.0
TI_CO = (-1.0, 0.5)
YTRANS, NTRANS = 0.8, 3

OUT_TXT = 'data/stage10i_isoladder.txt'
OUT_NPZ = 'data/stage10i_profiles.npz'
LINES = []
def emit(s=""):
    LINES.append(s)
    print(s, flush=True)
    with open(OUT_TXT, 'w') as f:
        f.write("\n".join(LINES) + "\n")

def mins():
    return f"[{(time.time()-T0)/60:.1f} min]"

# ---------- data (9V loader verbatim) ----------
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
GNAMES = [gal_name[ug[i]] for i in range(NGal)]
DSPARC = np.array([meta[gal_name[ug[i]]][2] for i in range(NGal)])

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

# ---------- gates (6J/9V machinery verbatim) ----------
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
        nm = GNAMES[i]
        if nm in cat:
            raw[i] = s_of(cat[nm])**2
    matched = np.isfinite(raw)
    med = float(np.median(raw[matched]))
    return np.where(matched, raw, med), med, int(matched.sum())

G_MX, G_MX_MED, n_mx = gate_arrays(chae_mx)

# ---------- family + objective + fitters (9V verbatim) ----------
def nu_p(y, p):
    yc = np.clip(y, 1e-14, None)
    ex = np.exp(-np.minimum(yc**p, 60.0))
    return (1.0-ex)**(-1.0/(2.0*p))

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

# ==================== PHASE 0: UNGC parse + census (G10I-0) ================
emit("STAGE 10I THE ISOLATION LADDER: %d galaxies, %d pts + %d lensing" %
     (kept, len(gobs), int(lmask.sum())))
emit("pre-reg: this commit. 9V archive: data/stage9v_profiles.npz")
emit("")

PFX = ('NGC', 'UGC', 'IC', 'DDO', 'ESO', 'PGC', 'UGCA')
def canon(nm):
    s = re.sub(r'[\s\-_]', '', str(nm)).upper()
    for p in PFX:
        if s.startswith(p):
            rest = s[len(p):]
            m = re.match(r'0*(\d.*)', rest)
            if m: return p + m.group(1)
    return s

ungc = []
with open('data/karachentsev_ungc.csv', encoding='utf-8',
          errors='replace') as f:
    for row in csv.DictReader(f):
        try:
            ra = float(row['RAJ2000']); de = float(row['DEJ2000'])
            dist = float(row['Dist'])
        except (ValueError, KeyError):
            continue
        def ff(k):
            try: return float(row.get(k, ''))
            except ValueError: return np.nan
        ungc.append(dict(name=row['Name'].strip(),
                         simbad=row.get('SimbadName', '').strip(),
                         ned=row.get('NEDname', '').strip(),
                         ra=ra, de=de, dist=dist,
                         klum=ff('KLum'), ti1=ff('Ti1'), ti5=ff('Ti5'),
                         md=row.get('MD', '').strip()))
emit("UNGC parsed: %d rows with (RA, Dec, Dist)" % len(ungc))

# 3D positions (Mpc)
U_POS = np.array([[u['dist']*math.cos(math.radians(u['de']))
                   * math.cos(math.radians(u['ra'])),
                   u['dist']*math.cos(math.radians(u['de']))
                   * math.sin(math.radians(u['ra'])),
                   u['dist']*math.sin(math.radians(u['de']))]
                  for u in ungc])
U_M = np.array([10.0**u['klum']*ML_K if np.isfinite(u['klum']) else 0.0
                for u in ungc])

ucanon = {}
for k, u in enumerate(ungc):
    for nm in (u['name'], u['simbad'], u['ned']):
        if nm:
            ucanon.setdefault(canon(nm), k)

match_idx = np.full(NGal, -1)
for i in range(NGal):
    k = ucanon.get(canon(GNAMES[i]), -1)
    if k >= 0:
        du = ungc[k]['dist']
        if abs(DSPARC[i] - du) <= max(2.0, 0.35*du):
            match_idx[i] = k
MATCHED = match_idx >= 0
N_MATCH = int(MATCHED.sum())

TI1 = np.full(NGal, np.nan)
for i in range(NGal):
    if MATCHED[i]:
        TI1[i] = ungc[match_idx[i]]['ti1']
ISO = MATCHED & (TI1 < TI_PRIMARY)
NONISO_LV = MATCHED & ~(TI1 < TI_PRIMARY)

gN0 = g_gas + g_dsk + g_bul
ntrans_pts = np.array([int(np.sum(gN0[GIDXS[i]]/A0_FID >= YTRANS))
                       for i in range(NGal)])
TRANS = ntrans_pts >= NTRANS
N_ISO = int(ISO.sum())
N_ISO_TR = int((ISO & TRANS).sum())

emit("CENSUS: matched %d / %d SPARC (name+distance); iso (Ti1<%.1f) %d; "
     "iso transition-crossers (>=%d pts y>=%.1f) %d" %
     (N_MATCH, NGal, TI_PRIMARY, N_ISO, NTRANS, YTRANS, N_ISO_TR))
emit("  co-strata: Ti1<-1.0: %d | Ti1<+0.5: %d | non-iso LV: %d "
     "(transition-crossers among them: %d)" %
     (int((MATCHED & (TI1 < -1.0)).sum()),
      int((MATCHED & (TI1 < 0.5)).sum()), int(NONISO_LV.sum()),
      int((NONISO_LV & TRANS).sum())))
iso_names = sorted((GNAMES[i], round(float(TI1[i]), 2))
                   for i in range(NGal) if ISO[i])
emit("  iso list: %s" % iso_names)
unmatched_lv = sorted(GNAMES[i] for i in range(NGal)
                      if not MATCHED[i] and DSPARC[i] <= 12.0)
emit("  unmatched SPARC with D<=12 Mpc (for the record): %s" %
     unmatched_lv)

g0 = (N_MATCH >= 20) and (N_ISO >= 6) and (N_ISO_TR >= 3)
emit("G10I-0 census: %s (bars 20/6/3)" % ("PASS" if g0 else
     "STOP -> I-FEASIBILITY-LIMITED"))
emit("")

if not g0:
    emit("LETTER: I-FEASIBILITY-LIMITED (pre-signed HOLD 8 / HOLD 53)")
    raise SystemExit(0)

# ==================== PHASE 1: gate construction (G10I-2) ==================
D2ALL = np.sum((U_POS[None, :, :] - U_POS[:, None, :])**2, axis=2)
np.fill_diagonal(D2ALL, np.inf)

def e_ngc_of(k, floor):
    d2 = D2ALL[k]
    sel = d2 <= NEIGH_MPC**2
    e = np.sum(G_SI*U_M[sel]*MSUN/(d2[sel]*MPC_M**2))/A0_FID
    return e + floor

EN_UNGC = {}
for fl in FLOORS:
    EN_UNGC[fl] = np.array([e_ngc_of(match_idx[i], fl) if MATCHED[i]
                            else np.nan for i in range(NGal)])

# G10I-2a: Theta re-derivation slope
tt_x, tt_y = [], []
for i in range(NGal):
    if not MATCHED[i] or not np.isfinite(TI1[i]): continue
    k = match_idx[i]
    d2 = D2ALL[k]
    sel = (d2 <= 5.0**2) & (U_M > 0)
    if not sel.any(): continue
    val = np.max(np.log10(U_M[sel]) - 1.5*np.log10(d2[sel]))
    tt_x.append(val); tt_y.append(TI1[i])
tt_x, tt_y = np.array(tt_x), np.array(tt_y)
A = np.vstack([tt_x, np.ones_like(tt_x)]).T
slope, off = np.linalg.lstsq(A, tt_y, rcond=None)[0]
g2a = abs(slope - 1.0) <= 0.25
emit("G10I-2a Theta re-derivation: slope %.3f (bar |s-1|<=0.25) offset "
     "%.2f (recorded) N=%d -> %s" %
     (slope, off, len(tt_x), "PASS" if g2a else "FAIL/STOP"))

# G10I-2b: UNGC e_N vs Chae e_N on the overlap
ov = [i for i in range(NGal) if MATCHED[i] and GNAMES[i] in chae_mx]
en_u = np.array([EN_UNGC[FLOOR_PRIMARY][i] for i in ov])
en_c = np.array([chae_mx[GNAMES[i]] for i in ov])
if len(ov) >= 8:
    rho, pval = spearmanr(en_u, en_c)
else:
    rho, pval = np.nan, np.nan
g2b = (len(ov) >= 8) and (rho > 0)
emit("G10I-2b UNGC-vs-Chae ambient cross-validation: overlap N=%d, "
     "Spearman rho=%.3f (p=%.3g), median ratio UNGC/Chae=%.3f -> %s" %
     (len(ov), rho, pval,
      float(np.median(en_u/en_c)) if len(ov) else np.nan,
      "PASS" if g2b else "FAIL/STOP"))
if len(ov):
    emit("  (calibration row: e_N percentiles 16/50/84  UNGC %s | "
         "Chae %s, a0 units)" %
         (np.percentile(en_u, [16, 50, 84]).round(4).tolist(),
          np.percentile(en_c, [16, 50, 84]).round(4).tolist()))

# ISO treatment gate arrays
G_ISO = {}
for fl in FLOORS:
    arr = G_MX.copy()
    for i in range(NGal):
        if MATCHED[i]:
            arr[i] = s_of(EN_UNGC[fl][i])**2
    G_ISO[fl] = arr
GATES_ISO = G_ISO[FLOOR_PRIMARY]
med_iso = float(np.median(GATES_ISO[ISO]))
med_non = float(np.median(GATES_ISO[NONISO_LV])) if NONISO_LV.any() else np.nan
g2c = med_iso > med_non
emit("G10I-2c axis premise: gate median iso %.4f vs non-iso LV %.4f -> %s"
     % (med_iso, med_non,
        "PASS" if g2c else "FAIL (ordering readout disabled)"))
pred_dp = 0.5*(med_iso - med_non)/2.0
emit("  r=1/2 predicted stratum contrast Dp = %.4f (LV grade)" % pred_dp)
emit("  gate range: iso 16/50/84 %s | Chae pgmax matched median %.4f" %
     (np.percentile(GATES_ISO[ISO], [16, 50, 84]).round(4).tolist(),
      G_MX_MED))
emit("")
if not (g2a and g2b):
    emit("LETTER: I-FEASIBILITY-LIMITED (gate construction failed; "
         "pre-signed HOLD 8 / HOLD 53)")
    raise SystemExit(0)

# ==================== PHASE 2: reader identity (G10I-1) ====================
emit("PHASE 2  reader identity vs the 9V archive (lnL grade)")
Z9V = np.load('data/stage9v_profiles.npz')
RGRID = Z9V['rgrid']
arch = Z9V['pgmax']

def p_pts_of_gates(r, gv, gl):
    return 0.5 + r*gv[gidx]/2.0, 0.5 + r*gl/2.0

prof_re = []
th_w = dml_w = dv_w = None
best_pack = {}
for r in RGRID:
    P_PTS, p_lens = p_pts_of_gates(float(r), G_MX, G_MX_MED)
    b, dml_w, dv_w = fit_vr(P_PTS, p_lens, ONES, lgobs, l_gobs,
                            th0=th_w, dml0=dml_w, dv0=dv_w)
    th_w = b.x
    prof_re.append(b.fun)
    best_pack[float(r)] = (b.x.copy(), dml_w.copy(), dv_w.copy(), b.fun)
prof_re = np.array(prof_re)
d_arch = prof_re - arch
g1a = float(np.max(np.abs(d_arch))) <= 0.5
emit("G10I-1a archive reproduction: max|d| = %.4f over %d pts (bar 0.5) "
     "-> %s  %s" % (float(np.max(np.abs(d_arch))), len(RGRID),
                    "PASS" if g1a else "FAIL/STOP", mins()))

b0, dml0_, dv0_ = fit_vr(np.full(len(gobs), 0.5), 0.5, ONES, lgobs, l_gobs)
g1b = abs(b0.fun - (-12152.49)) <= 1.0
emit("G10I-1b r=0 anchor: %.2f vs -12152.49 (bar 1.0) -> %s  %s" %
     (b0.fun, "PASS" if g1b else "FAIL/STOP", mins()))
emit("")
if not (g1a and g1b):
    emit("LETTER: STOP (wiring). No sky read.")
    raise SystemExit(0)

# G10I-4 arithmetic
S_HALF = math.sin(0.5)**2
g4 = []
g4.append(abs((0.5 + 0.5*0.7536/2) - 0.6884) < 5e-4)
g4.append(abs((0.5 + S_HALF*0.7536/2) - 0.5866) < 5e-4)
g4.append(abs(2*(0.6471-0.5)/0.7536 - 0.3904) < 5e-4)
emit("G10I-4 arithmetic: %s" % ("PASS" if all(g4) else "FAIL"))
emit("")

# ==================== PHASE 3: injections (G10I-3, sky-blind) ==============
emit("PHASE 3  injections, SKY-BLIND on the iso axis (rng 202)")
rng = np.random.default_rng(202)
i_arch = int(np.argmin(prof_re))
r_arch = float(RGRID[i_arch])
th_a, dml_a, dv_a, _ = best_pack[r_arch]

def mock_sky(r_true, rng_):
    la0, f, s_int, dlt = th_a
    P_PTS, p_lens = p_pts_of_gates(r_true, G_MX, G_MX_MED)
    fac = f*np.exp(dml_a[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    mu = np.log10(gN*nu_p(gN/10**la0, P_PTS)) + dv_a[gidx]
    sky = mu + rng_.normal(0, np.sqrt(sig2 + s_int*s_int))
    lg = l_gbar + dlt
    mul = lg + np.log10(nu_p(10**lg/10**la0, p_lens))
    lens = mul + rng_.normal(0, np.sqrt(l_sig2))
    return sky, lens

def profile_fit(gv, gl, sky, lens, r_grid, th0=None):
    vals, packs = [], {}
    thw = th0; dmlw = dvw = None
    for r in r_grid:
        P_PTS, p_lens = p_pts_of_gates(float(r), gv, gl)
        b, dmlw, dvw = fit_vr(P_PTS, p_lens, ONES, sky, lens,
                              th0=thw, dml0=dmlw, dv0=dvw)
        thw = b.x
        vals.append(b.fun)
        packs[float(r)] = (b.x.copy(), dmlw.copy(), dvw.copy(), b.fun)
    return np.array(vals), packs

GL_ISO = float(np.median(GATES_ISO[np.isfinite(GATES_ISO)]))
inj_res = []
for r_true in (0.34, 0.50):
    grid = np.arange(r_true-0.15, r_true+0.151, 0.05)
    sky_m, lens_m = mock_sky(r_true, rng)
    vals, _ = profile_fit(GATES_ISO, GL_ISO, sky_m, lens_m, grid,
                          th0=th_a)
    rrec, scurv, edge = parab_min(grid, vals)
    bar = max(0.05, 2*(scurv if scurv else 0.05))
    ok = (abs(rrec - r_true) <= bar) and not edge
    inj_res.append(ok)
    emit("  r_true=%.2f: rec %.4f (sigma_curv %s, err %+.4f, bar %.3f)%s "
         "%s  %s" % (r_true, rrec,
                     f"{scurv:.4f}" if scurv else "edge", rrec-r_true,
                     bar, " EDGE" if edge else "",
                     "OK" if ok else "FAIL", mins()))
g3 = all(inj_res)
emit("G10I-3 recovery: %s" % ("PASS" if g3 else
     "FAIL -> letter I-POWER (no iso sky read)"))

# Dp power leg on the r_true = 0.50 mock (12 paired reps)
emit("  Dp power leg (r_true = 0.50 mock, 12 paired galaxy reps):")
sky_m, lens_m = mock_sky(0.50, rng)
iso_ids = np.where(ISO)[0]
non_ids = np.where(~ISO)[0]
def strat_p(w_g, sky, lens, pgrid, th0=None):
    vals = []
    thw = th0; dmlw = dvw = None
    for p in pgrid:
        b, dmlw, dvw = fit_vr(np.full(len(gobs), p), p, w_g, sky, lens,
                              th0=thw, dml0=dmlw, dv0=dvw,
                              tol=0.2, max_rounds=6, sweeps=2)
        thw = b.x
        vals.append(b.fun)
    ph, sc, ed = parab_min(pgrid, vals)
    return ph, ed
PGRID_I = np.arange(0.55, 0.851, 0.075)
PGRID_N = np.arange(0.50, 0.801, 0.075)
dps = []
rngb = np.random.default_rng(53)
for rep in range(12):
    wi = np.zeros(NGal); wn = np.zeros(NGal)
    for i in rngb.choice(iso_ids, len(iso_ids), replace=True): wi[i] += 1
    for i in rngb.choice(non_ids, len(non_ids), replace=True): wn[i] += 1
    pi, e1 = strat_p(wi, sky_m, lens_m, PGRID_I, th0=th_a)
    pn, e2 = strat_p(wn, sky_m, lens_m, PGRID_N, th0=th_a)
    dps.append(pi - pn)
    if rep % 4 == 3:
        emit("    rep %d/12: Dp mean %.4f sd %.4f  %s" %
             (rep+1, float(np.mean(dps)), float(np.std(dps, ddof=1)),
              mins()))
sd_dp = float(np.std(dps, ddof=1))
dp_powered = sd_dp <= 0.011
emit("  sigma(Dp) = %.4f (power bar 0.011; r=1/2 predicted Dp %.4f) -> "
     "ordering %s" % (sd_dp, pred_dp,
                      "POWERED" if dp_powered else
                      "UNPOWERED (co-read only)"))
emit("")

# ==================== PHASE 4: sky (only if G10I-3 passed) =================
np.savez(OUT_NPZ, rgrid=RGRID, pgmax_re=prof_re, gates_iso=GATES_ISO,
         iso_mask=ISO.astype(int), ti1=TI1, en_ungc=EN_UNGC[FLOOR_PRIMARY],
         dp_mock=np.array(dps))

if not g3:
    emit("=" * 68)
    emit("LETTER: I-POWER (G10I-3 recovery FAIL; iso sky profile never")
    emit("computed -- the pre-registered firewall). SHIPPED: the census")
    emit("(%d iso / %d iso transition-crossers at public SPARC-LV grade)," %
         (N_ISO, N_ISO_TR))
    emit("the UNGC-vs-Chae calibration row, the Dp power spec "
         "(sigma %.4f vs signal %.4f), the 9V reopen clause executed." %
         (sd_dp, pred_dp))
    emit("pre-signed credence: HOLD mech 8; HOLD anomaly-real 53")
    emit("gates: G0 PASS | G1 PASS | G2 %s/%s/%s | G3 FAIL | G4 %s" %
         ("PASS" if g2a else "FAIL", "PASS" if g2b else "FAIL",
          "PASS" if g2c else "FAIL", "PASS" if all(g4) else "FAIL"))
    emit("total wall-clock %s" % mins())
    raise SystemExit(0)

emit("PHASE 4  iso-treatment sky profile (23 pts)")
vals_iso, packs_iso = profile_fit(GATES_ISO, GL_ISO, lgobs, l_gobs,
                                  RGRID, th0=th_a)
rh, sc, edge = parab_min(RGRID, vals_iso)
d0 = float(vals_iso[int(np.argmin(vals_iso))] - vals_iso[0])
emit("  min: r-hat = %.4f (sigma_curv %s)%s; d(-2lnL) vs r=0: %.2f  %s" %
     (rh, f"{sc:.4f}" if sc else "n/a", " EDGE" if edge else "", d0,
      mins()))
g6 = (not edge) and (sc is not None)
emit("G10I-6 interior/convex: %s" % ("PASS" if g6 else "FAIL"))

for fl in [f for f in FLOORS if f != FLOOR_PRIMARY]:
    gv = G_ISO[fl]
    gl = float(np.median(gv[np.isfinite(gv)]))
    vals_f, _ = profile_fit(gv, gl, lgobs, l_gobs,
                            np.arange(max(0.0, rh-0.15), rh+0.151, 0.075),
                            th0=th_a)
    rf, sf, ef = parab_min(np.arange(max(0.0, rh-0.15), rh+0.151, 0.075),
                           vals_f)
    emit("  floor %.3f co-read: r-hat = %.4f%s  %s" %
         (fl, rf, " EDGE" if ef else "", mins()))

# strata Dp on the sky (co-read; letter weight only if dp_powered)
pi_sky, _ = strat_p(ISO.astype(float), lgobs, l_gobs, PGRID_I, th0=th_a)
pn_sky, _ = strat_p((~ISO).astype(float), lgobs, l_gobs, PGRID_N,
                    th0=th_a)
emit("  sky strata: p_iso = %.4f | p_field = %.4f | Dp = %+.4f "
     "(sigma %.4f, %s)  %s" %
     (pi_sky, pn_sky, pi_sky-pn_sky, sd_dp,
      "POWERED" if dp_powered else "co-read", mins()))

# bootstrap (G10I-5, 9V verbatim rules)
emit("PHASE 5  bootstrap (40 paired reps, rng 53, warm-lite, iso)")
rngb2 = np.random.default_rng(53)
grid_b = np.array([rh-0.15, rh-0.075, rh, rh+0.075, rh+0.15])
boot, edges = [], 0
for rep in range(40):
    w = np.zeros(NGal)
    for i in rngb2.choice(np.arange(NGal), NGal, replace=True): w[i] += 1
    vals_b, _ = profile_fit(GATES_ISO, GL_ISO, lgobs, l_gobs, grid_b,
                            th0=packs_iso[min(packs_iso,
                                              key=lambda r:
                                              packs_iso[r][3])][0])
    rb, _, eb = parab_min(grid_b, vals_b)
    if eb: edges += 1
    boot.append(rb)
    if rep % 10 == 9:
        emit("  rep %d/40: mean %.4f sd %.4f edges %d  %s" %
             (rep+1, float(np.mean(boot)), float(np.std(boot, ddof=1)),
              edges, mins()))
clipped = False
if edges > 8:
    emit("  AMENDMENT A-widen FIRED (%d > 8): +/- 0.225 step 0.075" %
         edges)
    grid_w = np.array([max(0.0, rh-0.225), rh-0.15, rh-0.075, rh,
                       rh+0.075, rh+0.15, rh+0.225])
    boot2, edges2 = [], 0
    rngb3 = np.random.default_rng(53)
    for rep in range(40):
        w = np.zeros(NGal)
        for i in rngb3.choice(np.arange(NGal), NGal, replace=True):
            w[i] += 1
        keep = boot[rep]
        if abs(keep - grid_b[0]) > 1e-9 and abs(keep - grid_b[-1]) > 1e-9:
            boot2.append(keep); continue
        vals_b, _ = profile_fit(GATES_ISO, GL_ISO, lgobs, l_gobs, grid_w,
                                th0=th_a)
        rb, _, eb = parab_min(grid_w, vals_b)
        if eb: edges2 += 1
        boot2.append(rb)
        emit("    widened rep %d: %.4f%s  %s" %
             (rep, rb, " EDGE" if eb else "", mins()))
    boot = boot2
    clipped = edges2 > 8
    if clipped:
        emit("  post-widen persistent edges: %d -> sigma_r CLIPPED "
             "LOWER BOUND" % edges2)
sig_r = float(np.std(boot, ddof=1))
emit("bootstrap: sigma_r = %.4f%s; percentiles 16/50/84 = %s" %
     (sig_r, " (clipped lower bound)" if clipped else "",
      np.percentile(boot, [16, 50, 84]).round(4).tolist()))
emit("")

# ==================== verdict ====================
emit("=" * 68)
emit("PRIMARY (vertical + iso gates + bootstrap): r-hat = %.4f +/- %.4f%s"
     % (rh, sig_r, " (clipped)" if clipped else ""))
emit("  9V co-quote: 0.3365 +/- 0.1869 (clipped); Chae-gate treatments "
     "fid 0.3895 / noclust 0.3085")
kill = g3 and (not clipped) and (sig_r <= 0.15) and (rh + 2*sig_r < 0.5)
sharp = (sig_r <= 0.07) and ((rh+2*sig_r < 0.2298) or (rh-2*sig_r >
        0.2298) or (rh+2*sig_r < 0.5) or (rh-2*sig_r > 0.5))
measured = (sig_r <= 0.15) and not clipped
ordered = dp_powered and abs(pi_sky-pn_sky)/sd_dp >= 2
if kill:
    letter = "I-KILL-FIRED"
elif sharp:
    letter = "I-SHARP"
elif ordered:
    letter = "I-ORDERED" if pi_sky > pn_sky else "I-ANTI"
elif measured:
    letter = "I-MEASURED"
else:
    letter = "I-POWER"
emit("VERDICT LETTER: %s" % letter)
emit("pre-signed credence cell: see docstring map")
emit("gates: G0 PASS | G1 PASS | G2 %s/%s/%s | G3 PASS | G4 %s | "
     "G5 %s | G6 %s" %
     ("PASS" if g2a else "FAIL", "PASS" if g2b else "FAIL",
      "PASS" if g2c else "FAIL", "PASS" if all(g4) else "FAIL",
      "CLIPPED" if clipped else "PASS", "PASS" if g6 else "FAIL"))
emit("total wall-clock %s" % mins())
