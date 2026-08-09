"""
STAGE 10K (2026-08-09): THE GROUP-END POWER SPEC -- pricing the R30
inverted dividend before any contest is built.

QUESTION.  ROUND 30 (10I adoption): the point-source neighbor meter is
ASYMMETRIC -- null at the isolated end (<= 1.7% of any LSS floor) but
with real dynamic range at the group end (NGC2976's neighbor term =
110% of the 0.01 floor).  The named dividend: "a 'does the gate CLOSE
(p -> 1/2) toward dense environments' test is powered here ... a real
out-of-sample test of p = 1/2 + r*g/2 at the low-g end."  Powered for
the METER is not powered for the LIKELIHOOD: this stage prices the
group-end contest against the program's two measured walls -- the 9U
realization wall sigma_p = 0.075 and the 9V bootstrap wall sigma_r =
0.1869 -- WITHOUT ever reading the sky's own preference.  It also
executes R30 condition 7 (the alias-table census fix, counts only).

DATA.  SPARC + Mistele lensing + Chae Table 3 + UNGC (all 10I
verbatim) + data/sparc_simbad_aliases.csv (calcs/fetch_sparc_aliases.py,
SIMBAD TAP cross-IDs for the 28 unmatched-at-D<=12 SPARC names, fetched
at this pre-reg).  data/stage9v_profiles.npz = the archived 9V primary
profile (reader-identity target).

CONSTRUCTION (locked here, before any number is computed):
  Matcher modes: m0 = 10I name-only (regression target 31/13/3);
    m1 = m0 + alias table (canon() applied to every SIMBAD alias; the
    SAME 10I distance gate |D_SPARC - D_UNGC| <= max(2, 0.35 D_UNGC)).
  Strata (matched set, m1): ISO Ti1 < 0 (10I verbatim); GROUP Ti1 > 0
    primary, Ti1 >= 1 co-read.  Transition-crosser: >= 3 pts y >= 0.8
    (10I verbatim).
  Gate worlds: A = the 9V pgmax arrays (Chae maxclust, median-filled)
    VERBATIM; B = for m1-matched galaxies with a Chae row, e_B =
    e_Chae(maxclust) + neighbor sum (floorless, 10I e_ngc_of verbatim);
    all other galaxies identical to A.  B is an UPPER BOUND on the
    available group-end signal: Chae's 2M++ field may already contain
    part of the neighbor term (double-count direction disclosed; a
    dead upper bound is decisive, a live one is not).
  Power statistic (Fisher with explicit dv-profiling): at r in
    {0.3365 (9V operative), 0.50 (co-read)}, Dmu_i = log10 nu_p(y_i,
    p_B) - log10 nu_p(y_i, p_A) at the FITTED (la0, f, dml) of the
    licensed reproduction fit; per galaxy the constant part of Dmu is
    profiled against the dv channel exactly (absorbed fraction
    W/(W+P), W = sum 1/sig_eff^2, P = 1/SIGV^2), leaving S_raw =
    sum_g [ sum (Dmu - m_g)^2/sig_eff^2 + m_g^2 W P/(W+P) ].  Global
    nuisance absorption is bounded by the MEASURED deflation of the
    same statistic on the 9V r-axis: D_f = C_arch/C_fish with C_arch
    = 2/0.073^2 (the archived 9V curvature) and C_fish the identical
    Fisher construction for a small r-step on gate world A.  The
    quoted signal is S_def = S_raw * min(1, D_f).
  Stratum leg: Dp_group = r * mean_{i in GROUP}(g_B - g_A)/2 quoted
    against sigma_p = 0.075 (the 9U wall), both r values.

GATES (bars locked at this commit; STOP = wiring, no power quote):
  G10K-1  alias integrity: CSV present, >= 1 alias row, and the R30
          exemplar resolves (canon of some UGC05721 alias == NGC3274).
          FAIL -> alias legs DEGRADE to name-only (disclosed); power
          legs unaffected; census rows keep the 10I UNDERCOUNT label.
  G10K-2  matcher regression: m1 with an EMPTY alias table reproduces
          the 10I census EXACTLY (31 matched / 13 iso / 3 iso-TR and
          the identical iso list).  FAIL -> STOP.
  G10K-3  reader identity (lnL grade, the 8P/8Q rule): warm-chained
          fit_vr along the archived RGRID from 0 to the point nearest
          r = 0.3365 on gate world A; |fun - archived| <= 0.5 at that
          point.  FAIL -> STOP (no Fisher quotes).
  G10K-4  Fisher calibration: C_fish/C_arch in [0.8, 6.0].  Outside
          -> STOP (the approximation is not licensed).
  G10K-5  arithmetic (9V G9V-4 verbatim): p(0.5, 0.7536) = 0.6884;
          p(0.2298, 0.7536) = 0.5866; r_of_p(0.6471, fid) = 0.3904.
  G10K-6  the power bars (the verdict gate; see letters).

VERDICT LETTERS (mechanical on S_def at r = 0.3365 and the Dp ratio):
  K-POWERED      S_def >= 9 OR Dp_group/sigma_p >= 2: the group-end
                 contest is licensed as a NAMED SUCCESSOR (own pre-reg;
                 NOT run in this stage).
  K-MARGINAL     4 <= S_def < 9 OR 1 <= Dp/sigma_p < 2.
  K-POWER-DEAD   S_def < 4 AND Dp/sigma_p < 1: the group-end axis is
                 unpowered at SPARC x LV grade; parked with numbers.
  K-FEASIBILITY-LIMITED  any of G10K-2/3/4 STOP (wiring).

CREDENCE MAP (pre-signed): ALL cells HOLD mech-conditional 8, HOLD
anomaly-real 53.  This is a census + power-spec stage; the sky's
preference between gate worlds is NEVER read (firewall below), so no
cell can move a sky credence.  Bankable regardless of letter: the
alias-fixed census (R30 condition 7, counts only -- the 10I letter is
closed and does not reopen), the group-end power spec, the measured
nuisance-deflation factor D_f (reusable by any future gate-axis
design), the powered-galaxy list.

FIREWALL.  The likelihood is evaluated ONLY on gate world A (the
archived 9V treatment -- a licensed reproduction), warm-chained along
the archived grid.  No fit, profile, or lnL difference is ever
computed with B gates: the Fisher sums use model derivatives and the
fitted variance only.  The one sky-derived input is the archived 9V
curvature 0.073 (published).  Amendments, if any, are logged in this
docstring pre-quote with failed runs preserved.

COMPUTE ~15-30 min (8 warm-chained fits + arithmetic).
Writes data/stage10k_groupend.txt.
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
FLOOR_PRIMARY = 0.01

R_OP = 0.3365          # 9V operative r-hat
R_CO = 0.50            # dispersive co-read
SIG_P = 0.075          # 9U realization wall
C_ARCH = 2.0/0.073**2  # archived 9V curvature (-2lnL units)
TI_GROUP = 0.0         # group = Ti1 > 0 (primary); >= 1 co-read
POW_FRAC = 0.5         # "powered row" = neighbor term >= 50% of floor

OUT_TXT = 'data/stage10k_groupend.txt'
LINES = []
def emit(s=""):
    LINES.append(s)
    print(s, flush=True)
    with open(OUT_TXT, 'w') as f:
        f.write("\n".join(LINES) + "\n")

def mins():
    return f"[{(time.time()-T0)/60:.1f} min]"

# ---------- data (10I/9V loader verbatim) ----------
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
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]
GNAMES = [gal_name[ug[i]] for i in range(NGal)]
DSPARC = np.array([meta[gal_name[ug[i]]][2] for i in range(NGal)])

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

# ---------- gates (verbatim) ----------
def n_amb_of(e):
    x = math.sqrt(e)
    return 1.0/(math.exp(x) - 1.0)
def s_of(e):
    n = n_amb_of(e)
    return n/(1.0 + n)

chae_mx = {}
with open('data/chae2021_table3.csv') as f:
    for row in csv.DictReader(l for l in f if not l.startswith('#')):
        chae_mx[row['galaxy']] = 10.0**float(row['log_eN_maxclust'])

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

# ---------- family + objective + fitter (9V verbatim) ----------
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

ONES = np.ones(NGal)

# ==================== PHASE 0: matcher m0 regression + m1 ================
emit("STAGE 10K THE GROUP-END POWER SPEC: %d galaxies, %d pts + %d "
     "lensing" % (kept, len(gobs), int(lmask.sum())))
emit("pre-reg: this commit.  9V archive: data/stage9v_profiles.npz")
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

U_POS = np.array([[u['dist']*math.cos(math.radians(u['de']))
                   * math.cos(math.radians(u['ra'])),
                   u['dist']*math.cos(math.radians(u['de']))
                   * math.sin(math.radians(u['ra'])),
                   u['dist']*math.sin(math.radians(u['de']))]
                  for u in ungc])
U_M = np.array([10.0**u['klum']*ML_K if np.isfinite(u['klum']) else 0.0
                for u in ungc])
D2ALL = np.sum((U_POS[None, :, :] - U_POS[:, None, :])**2, axis=2)
np.fill_diagonal(D2ALL, np.inf)

def e_ngc_of(k, floor):
    d2 = D2ALL[k]
    sel = d2 <= NEIGH_MPC**2
    e = np.sum(G_SI*U_M[sel]*MSUN/(d2[sel]*MPC_M**2))/A0_FID
    return e + floor

ucanon = {}
for k, u in enumerate(ungc):
    for nm in (u['name'], u['simbad'], u['ned']):
        if nm:
            ucanon.setdefault(canon(nm), k)

gN0 = g_gas + g_dsk + g_bul
ntrans_pts = np.array([int(np.sum(gN0[GIDXS[i]]/A0_FID >= 0.8))
                       for i in range(NGal)])
TRANS = ntrans_pts >= 3

def run_match(alias_map):
    mi = np.full(NGal, -1)
    for i in range(NGal):
        cands = [canon(GNAMES[i])]
        cands += [canon(a) for a in alias_map.get(GNAMES[i], [])]
        for c in cands:
            k = ucanon.get(c, -1)
            if k >= 0:
                du = ungc[k]['dist']
                if abs(DSPARC[i] - du) <= max(2.0, 0.35*du):
                    mi[i] = k
                    break
        # name matched but distance failed on the first candidate:
        # keep scanning remaining candidates (loop above already does)
    return mi

# G10K-2: empty-alias regression vs the 10I census
mi0 = run_match({})
M0 = mi0 >= 0
TI1_0 = np.array([ungc[mi0[i]]['ti1'] if M0[i] else np.nan
                  for i in range(NGal)])
iso0 = M0 & (TI1_0 < 0.0)
n_m0, n_i0, n_it0 = int(M0.sum()), int(iso0.sum()), int((iso0 & TRANS).sum())
iso_list0 = sorted((GNAMES[i], round(float(TI1_0[i]), 2))
                   for i in range(NGal) if iso0[i])
ISO_LIST_10I = [('NGC0024', -1.0), ('NGC0891', -0.9), ('NGC2915', -1.1),
                ('NGC3741', -0.7), ('NGC4068', -0.5), ('NGC4559', -0.4),
                ('NGC5055', -0.1), ('NGC5585', -0.8), ('NGC6503', -1.1),
                ('NGC6789', -1.3), ('UGC01281', -1.1), ('UGC07690', -0.5),
                ('UGCA442', -0.1)]
g2 = (n_m0, n_i0, n_it0) == (31, 13, 3) and iso_list0 == ISO_LIST_10I
emit("G10K-2 matcher regression (empty alias table): %d/%d/%d vs 31/13/3,"
     " iso list %s -> %s" %
     (n_m0, n_i0, n_it0, "identical" if iso_list0 == ISO_LIST_10I else
      "DIFFERS", "PASS" if g2 else "FAIL/STOP"))
if not g2:
    emit("LETTER: K-FEASIBILITY-LIMITED (wiring)")
    raise SystemExit(0)

# G10K-1: alias table
ALIAS = {}
alias_ok = False
if os.path.exists('data/sparc_simbad_aliases.csv'):
    with open('data/sparc_simbad_aliases.csv') as f:
        rd = csv.reader(l for l in f if not l.startswith('#'))
        hdr = next(rd, None)
        for row in rd:
            if len(row) >= 2:
                ALIAS.setdefault(row[0], []).append(row[1])
    ex = any(canon(a) == 'NGC3274' for a in ALIAS.get('UGC05721', []))
    alias_ok = bool(ALIAS) and ex
emit("G10K-1 alias integrity: %d SPARC names with aliases, exemplar "
     "UGC05721->NGC3274 %s -> %s" %
     (len(ALIAS), "present" if alias_ok else "MISSING",
      "PASS" if alias_ok else "DEGRADE (name-only census, disclosed)"))

mi1 = run_match(ALIAS) if alias_ok else mi0
M1 = mi1 >= 0
TI1 = np.array([ungc[mi1[i]]['ti1'] if M1[i] else np.nan
                for i in range(NGal)])
ISO = M1 & (TI1 < 0.0)
GRP = M1 & (TI1 > TI_GROUP)
GRP1 = M1 & (TI1 >= 1.0)
new_matches = sorted((GNAMES[i], ungc[mi1[i]]['name'],
                      round(float(TI1[i]), 2))
                     for i in range(NGal) if M1[i] and not M0[i])
emit("")
emit("CENSUS (m1, alias-fixed; R30 condition 7 -- counts only):")
emit("  matched %d (was 31); iso %d (was 13); iso transition-crossers %d"
     " (was 3)" % (int(M1.sum()), int(ISO.sum()),
                   int((ISO & TRANS).sum())))
emit("  group Ti1>0: %d (transition-crossers %d); Ti1>=1: %d (TR %d)" %
     (int(GRP.sum()), int((GRP & TRANS).sum()), int(GRP1.sum()),
      int((GRP1 & TRANS).sum())))
emit("  new matches (SPARC -> UNGC, Ti1): %s" % new_matches)
iso_new = sorted((GNAMES[i], round(float(TI1[i]), 2))
                 for i in range(NGal) if ISO[i] and not M0[i])
emit("  newly matched isolated: %s" % iso_new)

# Spearman re-quote under m1 (counts-only row, recorded not barred)
EN_FL = {}
for i in range(NGal):
    if M1[i]:
        EN_FL[i] = e_ngc_of(mi1[i], FLOOR_PRIMARY)
ov = [i for i in range(NGal) if M1[i] and GNAMES[i] in chae_mx]
if len(ov) >= 8:
    rho, pval = spearmanr([EN_FL[i] for i in ov],
                          [chae_mx[GNAMES[i]] for i in ov])
else:
    rho, pval = np.nan, np.nan
emit("  UNGC-vs-Chae Spearman under m1: N=%d rho=%.3f (recorded; the 10I"
     " letter is closed)" % (len(ov), rho))

# neighbor terms (floorless) across the m1 matched set
neigh = {i: e_ngc_of(mi1[i], 0.0) for i in range(NGal) if M1[i]}
npow = sorted(((GNAMES[i], round(float(TI1[i]), 2),
                round(neigh[i]/FLOOR_PRIMARY*100, 1))
               for i in neigh if neigh[i] >= POW_FRAC*FLOOR_PRIMARY),
              key=lambda t: -t[2])
emit("  neighbor-powered rows (term >= 50%% of the 0.01 floor), "
     "(name, Ti1, %% of floor): %s" % npow)
nonzero = sorted(neigh.values())
emit("  neighbor-term percentiles 16/50/84 over matched: %s (%% of floor)"
     % [round(v/FLOOR_PRIMARY*100, 2)
        for v in np.percentile(nonzero, [16, 50, 84])])
emit("")

# ==================== PHASE 1: gate worlds ====================
G_A = G_MX.copy()
G_B = G_MX.copy()
aff = []
for i in range(NGal):
    if M1[i] and GNAMES[i] in chae_mx and neigh[i] > 0:
        eB = chae_mx[GNAMES[i]] + neigh[i]
        G_B[i] = s_of(eB)**2
        if G_A[i] - G_B[i] > 1e-6:
            aff.append(i)
emit("PHASE 1  gate worlds: %d galaxies with any gate change; largest:"
     % len(aff))
for i in sorted(aff, key=lambda j: G_B[j]-G_A[j])[:8]:
    emit("    %-10s Ti1 %+0.1f  e_Chae %.4f + neigh %.4f  g: %.4f -> "
         "%.4f  (Dp at r=%.4f: %+.5f)" %
         (GNAMES[i], TI1[i], chae_mx[GNAMES[i]], neigh[i], G_A[i],
          G_B[i], R_OP, R_OP*(G_B[i]-G_A[i])/2))
emit("  DOUBLE-COUNT CAVEAT: e_B adds the neighbor sum to Chae's 2M++ "
     "value; overlap makes B an UPPER BOUND on the signal.")

for lab, mask in (("group Ti1>0", GRP), ("group Ti1>=1", GRP1)):
    if mask.any():
        dpb = R_OP*float(np.mean(G_B[mask]-G_A[mask]))/2
        emit("  stratum %s (N=%d): mean Dgate %+.5f -> Dp(r=%.4f) = "
             "%+.6f = %.3f sigma_p" %
             (lab, int(mask.sum()), float(np.mean(G_B[mask]-G_A[mask])),
              R_OP, dpb, abs(dpb)/SIG_P))
GRP_DP = R_OP*float(np.mean(G_B[GRP]-G_A[GRP]))/2 if GRP.any() else 0.0
GRP_DP_CO = R_CO*float(np.mean(G_B[GRP]-G_A[GRP]))/2 if GRP.any() else 0.0
emit("  co-read r=0.50: Dp(group Ti1>0) = %+.6f = %.3f sigma_p" %
     (GRP_DP_CO, abs(GRP_DP_CO)/SIG_P))
emit("")

# G10K-5 arithmetic
S_HALF = math.sin(0.5)**2
g5 = []
g5.append(abs((0.5 + 0.5*0.7536/2) - 0.6884) < 5e-4)
g5.append(abs((0.5 + S_HALF*0.7536/2) - 0.5866) < 5e-4)
g5.append(abs(2*(0.6471-0.5)/0.7536 - 0.3904) < 5e-4)
emit("G10K-5 arithmetic: %s" % ("PASS" if all(g5) else "FAIL"))
emit("")

# ==================== PHASE 2: reader identity (G10K-3) ====================
emit("PHASE 2  licensed reproduction fit (gate world A, archived grid)")
Z9V = np.load('data/stage9v_profiles.npz')
RGRID = Z9V['rgrid']
arch = Z9V['pgmax']
i_t = int(np.argmin(np.abs(RGRID - R_OP)))
r_t = float(RGRID[i_t])
emit("  target grid point: r = %.4f (archived -2lnL %.4f)" %
     (r_t, arch[i_t]))

th_w = dml_w = dv_w = None
fun_t = None
pack_t = None
for j in range(i_t+1):
    r = float(RGRID[j])
    P_PTS = 0.5 + r*G_A[gidx]/2.0
    p_lens = 0.5 + r*G_MX_MED/2.0
    b, dml_w, dv_w = fit_vr(P_PTS, p_lens, ONES, lgobs, l_gobs,
                            th0=th_w, dml0=dml_w, dv0=dv_w)
    th_w = b.x
    if j == i_t:
        fun_t = b.fun
        pack_t = (b.x.copy(), dml_w.copy(), dv_w.copy())
    if j % 3 == 0 or j == i_t:
        emit("    r=%.2f: %.4f  %s" % (r, b.fun, mins()))
g3 = abs(fun_t - arch[i_t]) <= 0.5
emit("G10K-3 reader identity at r=%.4f: %.4f vs archived %.4f "
     "(|d| = %.4f, bar 0.5) -> %s" %
     (r_t, fun_t, arch[i_t], abs(fun_t - arch[i_t]),
      "PASS" if g3 else "FAIL/STOP"))
emit("")
if not g3:
    emit("LETTER: K-FEASIBILITY-LIMITED (wiring; no Fisher quotes)")
    raise SystemExit(0)

th_a, dml_a, dv_a = pack_t
la0, fml, s_int, dlt = th_a
se2 = sig2 + s_int*s_int
fac = fml*np.exp(dml_a[gidx])
gN = g_gas + fac*g_dsk + g_bul
Y = gN/10**la0

def dv_profiled_S(dmu):
    """Fisher signal with the per-galaxy constant profiled vs dv."""
    S = 0.0
    for i in range(NGal):
        mm = GIDXS[i]
        if len(mm) == 0: continue
        w = 1.0/se2[mm]
        W = float(np.sum(w))
        P = 1.0/SIGV[i]**2
        m_g = float(np.sum(w*dmu[mm]))/W
        S += float(np.sum(w*(dmu[mm]-m_g)**2)) + m_g*m_g*W*P/(W+P)
    return S

# ==================== PHASE 3: Fisher calibration (G10K-4) =================
DR = 0.02
P_A = 0.5 + r_t*G_A[gidx]/2.0
P_A2 = 0.5 + (r_t+DR)*G_A[gidx]/2.0
dmu_r = np.log10(nu_p(Y, P_A2)) - np.log10(nu_p(Y, P_A))
C_fish = 2.0*dv_profiled_S(dmu_r)/DR**2
ratio = C_fish/C_ARCH
D_f = min(1.0, C_ARCH/C_fish)
g4 = 0.8 <= ratio <= 6.0
emit("PHASE 3  Fisher calibration on the 9V r-axis:")
emit("  C_fish = %.1f vs C_arch = %.1f (ratio %.2f; band [0.8, 6.0]) "
     "-> %s" % (C_fish, C_ARCH, ratio, "PASS" if g4 else "FAIL/STOP"))
emit("  deflation factor D_f = min(1, C_arch/C_fish) = %.3f" % D_f)
emit("")
if not g4:
    emit("LETTER: K-FEASIBILITY-LIMITED (calibration outside band)")
    raise SystemExit(0)

# ==================== PHASE 4: the group-end power spec ====================
emit("PHASE 4  the group-end contest signal (B vs A, never read on sky)")
res = {}
for r_use in (R_OP, R_CO):
    P_Av = 0.5 + r_use*G_A[gidx]/2.0
    P_Bv = 0.5 + r_use*G_B[gidx]/2.0
    dmu = np.log10(nu_p(Y, P_Bv)) - np.log10(nu_p(Y, P_Av))
    S_raw = dv_profiled_S(dmu)
    S_def = S_raw*D_f
    res[r_use] = (S_raw, S_def)
    emit("  r = %.4f: S_raw = %.3f -> S_def = %.3f  (2sig bar 4, 3sig "
         "bar 9)" % (r_use, S_raw, S_def))
    percont = []
    for i in aff:
        mm = GIDXS[i]
        w = 1.0/se2[mm]
        W = float(np.sum(w)); P = 1.0/SIGV[i]**2
        m_g = float(np.sum(w*dmu[mm]))/W
        c = float(np.sum(w*(dmu[mm]-m_g)**2)) + m_g*m_g*W*P/(W+P)
        percont.append((GNAMES[i], round(c*D_f, 4)))
    emit("    per-galaxy S_def: %s" %
         sorted(percont, key=lambda t: -t[1])[:6])
S_OP = res[R_OP][1]
emit("")

# ==================== verdict (G10K-6) ====================
emit("=" * 68)
dp_ratio = abs(GRP_DP)/SIG_P
if (S_OP >= 9) or (dp_ratio >= 2):
    letter = "K-POWERED"
elif (S_OP >= 4) or (dp_ratio >= 1):
    letter = "K-MARGINAL"
else:
    letter = "K-POWER-DEAD"
emit("VERDICT LETTER: %s" % letter)
emit("  contest signal S_def(r=%.4f) = %.3f (bars 4/9); stratum "
     "Dp = %+.6f = %.3f sigma_p (bars 1/2)" %
     (R_OP, S_OP, GRP_DP, dp_ratio))
emit("  the group-end axis at SPARC x LV grade: %s" %
     ("licensed as a successor contest" if letter == "K-POWERED" else
      "borderline; successor spec printed, not licensed"
      if letter == "K-MARGINAL" else
      "unpowered; parked with numbers (the R30 dividend priced)"))
emit("BANKED: the alias-fixed census (R30 cond 7), the power spec, "
     "D_f = %.3f, the powered-row list." % D_f)
emit("pre-signed credence: HOLD mech 8; HOLD anomaly-real 53 (all "
     "cells; the sky's A-vs-B preference was never read)")
emit("gates: G1 %s | G2 PASS | G3 PASS | G4 PASS | G5 %s | G6 %s" %
     ("PASS" if alias_ok else "DEGRADE", "PASS" if all(g5) else "FAIL",
      letter))
emit("total wall-clock %s" % mins())
