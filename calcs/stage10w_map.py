"""STAGE 10W -- THE PER-GALAXY FORM-PREFERENCE MAP + THE
ANGULAR-RESOLUTION CUT.

The R48-C18 successor instrument: one per-galaxy map at fixed parent
globals replacing stratum contests (Part 1 descriptive, Part 2 the
joint axis instrument with coverage-conditioned nulls), plus the R48
named live rival -- inner angular resolution -- given its registered
test as NEW SKY (Part 3 cut refits). Pre-registration:
PREREG-MAP-10W.md + amendment A1 (both committed before any gate).

Engine = exec-inheritance of calcs/stage10v_strat.py truncated at the
baselines marker (the flowboot/GB route). Margins ride the GB R-d
realized-galaxy-block SD, RECOMPUTED on every cut world (the 10V
bar-mismatch lesson). NO injection-calibrated bar anywhere (trap #28
moot). Fits may ride calcs/fitpool.py ONLY under a recorded GATE
PASS certificate; otherwise serial (G10W-6).

Modes: gates (sky-blind: engine identity, decomposition closure,
arcsec census, cut wiring at theta = 0, null machinery on synthetic
data, pool spot-check; NO Part-2/3 statistic) | sky (map -> axes ->
cuts -> conditional boot -> letter).
"""
import math, os, sys, time
import numpy as np
from scipy.optimize import minimize_scalar as msc
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
MODE = sys.argv[1] if len(sys.argv) > 1 else 'gates'
SKY = (MODE == 'sky')
OUTFILE = ('data/stage10w_skyread.txt' if SKY
           else 'data/stage10w_gates.txt')

# ---------------- engine inheritance (flowboot route) --------------
SRCV = open(os.path.join(HERE, 'stage10v_strat.py'),
            encoding='utf-8').read()
cutv = SRCV.index("# ---------------- baselines (both modes)")
NSV = {'__name__': 'stage10v_trunc',
       '__file__': os.path.join(HERE, 'stage10v_strat.py')}
_argv = sys.argv
sys.argv = ['stage10v_strat.py', 'gates']
exec(compile(SRCV[:cutv], 'stage10v_trunc', 'exec'), NSV)
sys.argv = _argv

W78, LEGA, FLOW = NSV['W78'], NSV['LEGA'], NSV['FLOW']
FAMS, MEMBERS, ARC = NSV['FAMS'], NSV['MEMBERS'], NSV['ARC']
build_sub, fit_hier, fit_deep, contest_fit = (
    NSV['build_sub'], NSV['fit_hier'], NSV['fit_deep'],
    NSV['contest_fit'])
nu_be = NSV['nu_be']
S_ML, U_PRIOR, LN10 = NSV['S_ML'], NSV['U_PRIOR'], NSV['LN10']
NAMEOF = NSV['NAMEOF']
sub_fields, sd_block, lag1 = (NSV['sub_fields'], NSV['sd_block'],
                              NSV['lag1'])
window_world = NSV['window_world']
NS10U = NSV['NS']
SP, UD, UB, KPC = NS10U['SP'], NS10U['UD'], NS10U['UB'], NS10U['KPC']

sys.path.insert(0, HERE)
from fitpool import run_tasks  # noqa: E402

Lst = []
def P(s=""):
    print(s, flush=True)
    Lst.append(s)

def save():
    with open(OUTFILE, 'w', encoding='utf-8') as f:
        f.write("\n".join(Lst) + "\n")

t00 = time.time()
P(f"STAGE 10W MAP + RESOLUTION CUT -- mode = {MODE} "
  f"({'sky read' if SKY else 'sky-blind gates'})")
P("")

# ---------------- pool certificate (G10W-6 half 1) -----------------
POOL_OK = False
try:
    _gt = open('data/fitpool_gate.txt', encoding='utf-8').read()
    POOL_OK = 'GATE PASS' in _gt
except OSError:
    pass
P(f"pool certificate: {'GATE PASS on file' if POOL_OK else 'ABSENT/FAIL'}"
  f" -> Part-3 fits {'pooled' if POOL_OK else 'SERIAL (disclosed)'}")

def pool_or_serial(tasks, subs, tag):
    if POOL_OK:
        return run_tasks(tasks, subs, tag=tag, log=P)
    res = {}
    for i, t in enumerate(tasks):
        sub = subs[t['sub']]
        if t['kind'] == 'contest_fit':
            bf, gap = contest_fit(sub, t['member'], list(t['th0']))
        else:
            bf = fit_deep(sub, FAMS[t['member']], th0=list(t['th0']))
            gap = 0.0
        res[t['key']] = {'fun': float(bf.fun),
                         'x': [float(v) for v in bf.x],
                         'gap': float(gap)}
        P(f"    [serial {tag}] {t['key']} ({i+1}/{len(tasks)}) "
          f"[{time.time()-t00:.0f}s]")
    return res

# ---------------- per-galaxy decomposition (R48 indep_obj) ---------
def indep_pergal(sub, nm, th):
    """R48 indep_obj with the per-galaxy contributions captured.
    Returns (list of per-galaxy obj_g, total incl. the global u-prior
    term). obj_g = residual block + ln(se2) + its own dml/dv priors;
    the single global (u/U_PRIOR)^2 is EXCLUDED from shares and
    included in the total (disclosed)."""
    la0, f, s_int, u = th
    a0 = 10**la0
    nu = FAMS[nm]
    per = []
    for k in range(sub['n']):
        pts = sub['gpts'][k]
        lg = sub['lg'][pts]
        gg = sub['gg'][pts]; gd = sub['gd'][pts]; gb = sub['gb'][pts]
        s2 = sub['s2'][pts]
        sv = sub['sv'][k]
        isu = sub['isuma_g'][k]
        se2 = s2 + s_int*s_int
        dml, dv = 0.0, 0.0
        prev = None
        for _ in range(400):
            gN = gg + f*math.exp(dml)*gd + gb
            r0 = lg - np.log10(gN*nu(gN/a0)) - u*isu
            w = 1.0/se2
            dv = float(np.sum(w*r0)/(np.sum(w) + 1.0/sv**2))

            def od(dl):
                gN2 = gg + f*math.exp(dl)*gd + gb
                rr = lg - np.log10(gN2*nu(gN2/a0)) - dv - u*isu
                return float(np.sum(rr*rr/se2) + dl*dl/(S_ML*S_ML))
            dml = msc(od, bounds=(-0.7, 0.7), method='bounded').x
            gN = gg + f*math.exp(dml)*gd + gb
            r = lg - np.log10(gN*nu(gN/a0)) - dv - u*isu
            cur = (float(np.sum(r*r/se2 + np.log(se2)))
                   + dml*dml/(S_ML*S_ML) + dv*dv/(sv*sv))
            if prev is not None and abs(prev - cur) < 1e-10:
                prev = cur
                break
            prev = cur
        per.append(prev)
    return per, sum(per) + (u/U_PRIOR)**2

# ---------------- angular radii (two conventions) ------------------
import glob as _glob
ROTPATH = {os.path.basename(p).replace('_rotmod.dat', ''): p
           for p in _glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                               recursive=True)}

def theta_min_gb17(gals):
    """GB-17 VERBATIM (eV-filter only; catalog-D) -- the R48 census
    convention, reproduced for G10W-3."""
    th = []
    for g in gals:
        nm2 = NAMEOF[int(g)]
        pth = [ROTPATH[nm2]] if nm2 in ROTPATH else []
        if not pth: continue
        D0 = SP['dist'][int(g)][0]
        rmin = None
        for l in open(pth[0]):
            if l.startswith('#'): continue
            t = l.split()
            if len(t) < 6: continue
            R, Vo, eV = float(t[0]), float(t[1]), float(t[2])
            if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
            rmin = R if rmin is None else min(rmin, R)
        if rmin is not None and D0 > 0:
            th.append(rmin/D0*206.265)
    return float(np.median(th))

def build_theta():
    """World-aligned per-point angular radii: the SP parse loop
    VERBATIM (both filters: eV/Vo <= 0.10 AND gg+gd+gb > 0) so the
    kept lines align one-to-one with W78's per-galaxy points; catalog
    D. Returns (PT_THETA aligned to W78 arrays, per-galaxy innermost
    dict, per-galaxy count-match report)."""
    PT = np.full(len(W78['lgobs']), np.nan)
    THG = {}
    mism = []
    for g in LEGA + FLOW:
        nm2 = NAMEOF[int(g)]
        assert nm2 in ROTPATH, f"no rotmod for {nm2}"
        D0 = SP['dist'][int(g)][0]
        ths = []
        for l in open(ROTPATH[nm2]):
            if l.startswith('#'): continue
            t = l.split()
            if len(t) < 6: continue
            R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
            if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
            gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC
            gb = UB*Vb*Vb/R*KPC
            if gg+gd+gb <= 0: continue
            ths.append(R/D0*206.265)
        pts = W78['gpts'][int(g)]
        if len(ths) != len(pts):
            mism.append((nm2, len(ths), len(pts)))
            continue
        PT[pts] = ths
        THG[int(g)] = min(ths)
    return PT, THG, mism

# ---------------- axes (10V conventions, verbatim) -----------------
def fgas_of(g):
    pts = W78['gpts'][int(g)]
    den = W78['gg'][pts] + W78['gd'][pts] + W78['gb'][pts]
    return float(np.median(W78['gg'][pts]/den))

def cov_of(g):
    pts = W78['gpts'][int(g)]
    gN1 = W78['gg'][pts] + 1.0*W78['gd'][pts] + W78['gb'][pts]
    return float(np.median(np.log10(gN1/ARC['BE'])))

# ---------------- rank statistic + stratified permutation ----------
def spear(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    return float(np.sum(ra*rb)
                 / math.sqrt(np.sum(ra*ra)*np.sum(rb*rb)))

def _resid(v, Z):
    A = (np.column_stack([np.ones(len(v)), Z]) if Z.shape[1]
         else np.ones((len(v), 1)))
    beta, *_ = np.linalg.lstsq(A, v, rcond=None)
    return v - A @ beta

def partial_p(delta, axis, Zcols, rng, nperm=10000):
    """A4-i machinery: partial rank correlation (Spearman on
    least-squares rank residuals given the conditioning set) with a
    Freedman-Lane null (unrestricted permutation of the
    delta-residuals). Two-sided."""
    rd = rankdata(delta).astype(float)
    ra = rankdata(axis).astype(float)
    rZ = (np.column_stack([rankdata(z).astype(float) for z in Zcols])
          if Zcols else np.empty((len(delta), 0)))
    ed, ea = _resid(rd, rZ), _resid(ra, rZ)
    den_a = math.sqrt(float(np.sum(ea*ea)))
    obs = float(np.sum(ed*ea)) / (math.sqrt(float(np.sum(ed*ed)))
                                  * den_a)
    hits = 0
    for _ in range(nperm):
        e = ed[rng.permutation(len(ed))]
        r = float(np.sum(e*ea)) / (math.sqrt(float(np.sum(e*e)))
                                   * den_a)
        if abs(r) >= abs(obs) - 1e-15:
            hits += 1
    return obs, (1 + hits)/(1 + nperm)

# ---------------- the A4-ii cut builder ----------------------------
def cut_world(W, keep):
    """Point cut with the touched-galaxy rule: a galaxy the cut left
    UNTOUCHED is kept whatever its point count (it is in the
    baseline contests); a touched galaxy needs >= 3 surviving
    points. theta = 0 is therefore the exact identity."""
    idx = np.where(keep)[0]
    W2 = dict(W)
    for k in ('gg', 'gd', 'gb', 'lgobs', 'sig2', 'gal_id'):
        W2[k] = W[k][idx]
    ug = np.unique(W2['gal_id'])
    gp = {int(g): np.where(W2['gal_id'] == g)[0] for g in ug}
    keepg = set()
    for g in gp:
        if (len(gp[g]) == len(W['gpts'][int(g)])
                or len(gp[g]) >= 3):
            keepg.add(g)
    W2['gpts'] = {g: gp[g] for g in keepg}
    W2['ug'] = np.array(sorted(keepg))
    legw = [g for g in LEGA if g in keepg]
    floww = [g for g in FLOW if g in keepg]
    return W2, legw, floww

# ---------------- baselines (both modes; the 10V recipe) -----------
P("")
P("== baselines (the 10V G10V-1/G10V-2a recipe, verbatim) ==")
subP = build_sub(W78, LEGA)
b_plain = fit_hier(subP, nu_be, use_u=False)
b_loose = fit_hier(subP, nu_be, use_u=True, th0=list(b_plain.x)+[0.0])
arch = {'BE': fit_deep(subP, nu_be, th0=list(b_loose.x))}
for nm in ('p065', 'gm', 'boot'):
    arch[nm] = fit_deep(subP, FAMS[nm], th0=list(arch['BE'].x))
P(f"  anchored refits done  [{time.time()-t00:.0f}s]")
subF = build_sub(W78, FLOW)
fw = {}
for m in MEMBERS:
    fw[m] = fit_deep(subF, FAMS[m],
                     th0=[math.log10(ARC[m]), 1.0, 0.08, 0.0])
P(f"  flow refits done  [{time.time()-t00:.0f}s]")
REF_A = {'gm': 24.0, 'p065': 25.9, 'boot': 63.3}
REF_F = {'BE': 48.7, 'p065': 7.6, 'gm': 2.5, 'boot': 0.0}
fminF = min(fw[m].fun for m in MEMBERS)
D_UNCUT_FLOW = fw['boot'].fun - fw['BE'].fun     # = -48.7 archived
D_UNCUT_ANCH = arch['boot'].fun - arch['BE'].fun  # = +63.3 archived
WARM_A = {m: list(arch[m].x) for m in MEMBERS}
WARM_F = {m: list(fw[m].x) for m in MEMBERS}

g1_ok = True
P("  G10W-1 engine identity (A2 bars: anch delta 0.5 = the 10V bar; "
  "flow delta 0.05 vs GB R-i; a0 0.1%):")
for nm in MEMBERS:
    da = abs(10**arch[nm].x[0] - ARC[nm])/ARC[nm]
    dd = 0.0 if nm == 'BE' else abs(
        (arch[nm].fun - arch['BE'].fun) - REF_A[nm])
    df = abs((fw[nm].fun - fminF) - REF_F[nm])
    ok = (da <= 1e-3 and dd <= 0.5 and df <= 0.05)
    g1_ok &= ok
    P(f"    {nm:5s}: anch a0 d = {100*da:.3f}%; |d_anch - arch| = "
      f"{dd:.3f}; |d_flow - arch| = {df:.3f} "
      f"{'ok' if ok else 'OFF'}")
P(f"  G10W-1: {'PASS' if g1_ok else 'FAIL'}")

# ==================================================================
if not SKY:
    # ---------------- G10W-2 decomposition closure -----------------
    P("")
    P("-- G10W-2 decomposition closure (A3 sign-aware form) --")
    g2_ok = g1_ok
    for leg, sub, opt in (('anch', subP, arch), ('flow', subF, fw)):
        for nm in ('BE', 'boot'):
            per, tot = indep_pergal(sub, nm, list(opt[nm].x))
            gap = tot - opt[nm].fun
            ok = -0.1 <= gap <= 1e-6
            g2_ok &= ok
            P(f"  {leg}/{nm}: indep_total - engine fun = {gap:+.3e} "
              f"(A3 band [-0.1, 1e-6]) {'ok' if ok else 'OFF'}")
    P(f"  G10W-2: {'PASS' if g2_ok else 'FAIL'}")

    # ---------------- G10W-3 arcsec census -------------------------
    P("")
    P("-- G10W-3 arcsec census --")
    thF17 = theta_min_gb17(FLOW)
    thA17 = theta_min_gb17(LEGA)
    g3_ok = abs(thF17 - 15.0) <= 0.05 and abs(thA17 - 34.8) <= 0.05
    P(f"  GB-17 convention medians: flow {thF17:.2f} (rev 15.0), "
      f"anchored {thA17:.2f} (rev 34.8) "
      f"{'ok' if g3_ok else 'OFF'}")
    PT_THETA, TH_G, mism = build_theta()
    if mism:
        g3_ok = False
        for nm2, a, b in mism:
            P(f"  COUNT MISMATCH {nm2}: parse {a} vs world {b}")
    thFw = float(np.median([TH_G[g] for g in FLOW]))
    thAw = float(np.median([TH_G[g] for g in LEGA]))
    P(f"  world-aligned (both filters, the CUT convention): flow "
      f"{thFw:.2f}, anchored {thAw:.2f} (operative; difference vs "
      f"GB-17 disclosed)")
    rng31 = np.random.default_rng(31)
    ANCH_NONUMA = [g for g in LEGA if g not in set(int(x)
                   for x in W78['uma'])]
    UMA = [int(x) for x in W78['uma']]
    P("  trap-#30 spot-checks (2 per class, seeded 31):")
    for cls, lst in (('anch-nonUMa', ANCH_NONUMA), ('UMa', UMA),
                     ('flow', FLOW)):
        for g in rng31.choice(lst, size=2, replace=False):
            g = int(g)
            P(f"    {cls:11s} {NAMEOF[g]:12s}: n_pts = "
              f"{len(W78['gpts'][g])}, theta_min = {TH_G[g]:.1f} "
              f"arcsec")
    P(f"  G10W-3: {'PASS' if g3_ok else 'FAIL'}")

    # ---------------- G10W-4 cut wiring ----------------------------
    P("")
    P("-- G10W-4 cut wiring (A4-ii builder) --")
    keep0 = ~(PT_THETA < 0.0)          # NaN -> True; theta=0 cut
    W0, leg0, flow0 = cut_world(W78, keep0)
    g4_ok = (list(leg0) == list(LEGA) and list(flow0) == list(FLOW))
    for k in ('gg', 'gd', 'gb', 'lgobs', 'sig2', 'gal_id'):
        g4_ok &= bool(np.array_equal(W0[k], W78[k]))
    s0 = build_sub(W0, leg0)
    for k in ('lg', 'gg', 's2'):
        g4_ok &= bool(np.array_equal(s0[k], subP[k]))
    P(f"  theta = 0 identity (world arrays + anchored sub): "
      f"{'exact' if g4_ok else 'BROKEN'}")
    keep20 = ~(PT_THETA < 20.0)
    W20, leg20, flow20 = cut_world(W78, keep20)
    nptsA = sum(len(W78['gpts'][g]) for g in LEGA)
    nptsF = sum(len(W78['gpts'][g]) for g in FLOW)
    npts20A = sum(len(W20['gpts'][g]) for g in leg20)
    npts20F = sum(len(W20['gpts'][g]) for g in flow20)
    P(f"  20-arcsec census (NO fits, sky-blind): anchored "
      f"{npts20A}/{nptsA} pts, {len(leg20)}/78 gal; flow "
      f"{npts20F}/{nptsF} pts, {len(flow20)}/71 gal")
    dropF20 = [NAMEOF[g] for g in FLOW if g not in set(flow20)]
    P(f"  flow galaxies dropped at 20: "
      f"{', '.join(dropF20) if dropF20 else 'none'}")
    P(f"  G10W-4: {'PASS' if g4_ok else 'FAIL'}")

    # ---------------- G10W-5 null machinery (A2 form) --------------
    P("")
    P("-- G10W-5 null machinery (synthetic; A2 median over 25 "
      "draws; A4-i partial-rank machinery) --")
    gal_all = LEGA + FLOW
    cov = np.array([cov_of(g) for g in gal_all])
    fgas = np.array([fgas_of(g) for g in gal_all])
    legv = np.array([0.0]*len(LEGA) + [1.0]*len(FLOW))
    rng5 = np.random.default_rng(1313)
    p_un, p_st = [], []
    for _ in range(25):
        syn = 2.0*cov + rng5.normal(0.0, float(np.std(cov)),
                                    len(cov))
        _, pu = partial_p(syn, cov, [legv], rng5, nperm=2000)
        _, ps = partial_p(syn, fgas, [cov, legv], rng5, nperm=2000)
        p_un.append(pu); p_st.append(ps)
    med_un = float(np.median(p_un)); med_st = float(np.median(p_st))
    g5_ok = (med_un <= 1e-3) and (med_st >= 0.2)
    P(f"  planted-coverage synthetic: median coverage-axis p = "
      f"{med_un:.4f} (bar <= 0.001); median conditioned "
      f"composition-axis p = {med_st:.3f} (bar >= 0.2)")
    P(f"  G10W-5: {'PASS' if g5_ok else 'FAIL'}")

    # ---------------- G10W-6 pool spot-check -----------------------
    P("")
    P("-- G10W-6 pool spot-check (2 tasks, pool == serial) --")
    g6_ok = True
    if POOL_OK:
        tasks = [dict(key=f'spot:{m}', kind='fit_deep', sub='anch',
                      member=m, th0=WARM_A[m]) for m in ('BE', 'boot')]
        resP = run_tasks(tasks, {'anch': subP}, nworkers=2,
                         tag='spot', log=P)
        for m in ('BE', 'boot'):
            sref = fit_deep(subP, FAMS[m], th0=WARM_A[m])
            ok = (resP[f'spot:{m}']['fun'] == float(sref.fun)
                  and resP[f'spot:{m}']['x'] ==
                  [float(v) for v in sref.x])
            g6_ok &= ok
            P(f"  {m}: {'bit-identical' if ok else 'MISMATCH'}")
    else:
        P("  pool certificate absent -> stage runs serial "
          "(science unaffected; disclosed)")
    P(f"  G10W-6: {'PASS' if g6_ok else 'FAIL (serial fallback)'}")

    P("")
    allok = g1_ok and g2_ok and g3_ok and g4_ok and g5_ok
    P(f"== GATES {'ALL PASS' if allok else 'FAIL -> STOP'} "
      f"(G10W-6 informational) ==")
    P(f"wall-clock: {(time.time()-t00)/60:.1f} min")
    save()
    print("\nsaved:", OUTFILE)
    sys.exit(0 if allok else 1)

# ==================================================================
# SKY MODE
# ==================================================================
if not g1_ok:
    P("STOP: engine identity failed in sky mode")
    save(); sys.exit(1)

PT_THETA, TH_G, mism = build_theta()
assert not mism, "count mismatches must be resolved at gate time"

# ---------------- trap-#29 baseline grading (FIRST) ----------------
P("")
P("== trap-#29 baseline grading (before any Part-2/3 statistic) ==")
BASE = {}
for leg, sub, opt, dunc in (('anch', subP, arch, D_UNCUT_ANCH),
                            ('flow', subF, fw, D_UNCUT_FLOW)):
    A = 'BE' if opt['BE'].fun <= opt['boot'].fun else 'boot'
    B = 'boot' if A == 'BE' else 'BE'
    fldA = sub_fields(sub, A, list(opt[A].x))
    muB = sub_fields(sub, B, list(opt[B].x))[0]
    sdb = sd_block(sub, fldA, muB)
    BASE[leg] = dict(d=dunc, sd=sdb)
    P(f"  {leg}: d(BE,boot) = {dunc:+.1f}; SD_block = {sdb:.1f}; "
      f"|d|/SD = {abs(dunc)/sdb:.2f} sigma "
      f"(rho_lag1 = {lag1(sub, fldA):+.2f})")
P("  no Part-2/3 grammar may claim what these baselines cannot "
  "support (trap #29)")

# ---------------- PART 1: the map ---------------------------------
P("")
P("== PART 1: the per-galaxy form-preference map (fixed parent "
  "globals) ==")
PG = {}          # leg -> member -> per-galaxy list
for leg, sub, opt, seq in (('anch', subP, arch, LEGA),
                           ('flow', subF, fw, FLOW)):
    PG[leg] = {}
    for m in MEMBERS:
        per, tot = indep_pergal(sub, m, list(opt[m].x))
        PG[leg][m] = per
    P(f"  {leg}: decomposed 4 members  [{time.time()-t00:.0f}s]")

# per-galaxy SD^2 shares for the boot pair (fldA convention)
SD2SH = {}
for leg, sub, opt, seq in (('anch', subP, arch, LEGA),
                           ('flow', subF, fw, FLOW)):
    A = 'BE' if opt['BE'].fun <= opt['boot'].fun else 'boot'
    B = 'boot' if A == 'BE' else 'BE'
    muA, wA, resA = sub_fields(sub, A, list(opt[A].x))
    muB = sub_fields(sub, B, list(opt[B].x))[0]
    D = muB - muA
    terms = []
    for k in range(sub['n']):
        pts = sub['gpts'][k]
        terms.append(4.0*float(np.sum(wA[pts]*resA[pts]*D[pts]))**2)
    tot = sum(terms)
    SD2SH[leg] = [t/tot for t in terms]

CENSUS_PENDING = {'D564-8', 'D631-7'}
rows = []
for leg, seq in (('anch', LEGA), ('flow', FLOW)):
    dts = {m: [PG[leg][m][k] - PG[leg]['BE'][k]
               for k in range(len(seq))] for m in MEMBERS}
    tot_abs = sum(abs(v) for v in dts['boot'])
    for k, g in enumerate(seq):
        rows.append(dict(
            name=NAMEOF[int(g)], leg=leg,
            fgas=fgas_of(g), cov=cov_of(g), theta=TH_G[int(g)],
            d_boot=dts['boot'][k], d_p065=dts['p065'][k],
            d_gm=dts['gm'][k],
            abs_share=abs(dts['boot'][k])/tot_abs,
            sd2_share=SD2SH[leg][k],
            flag=('census-pending' if NAMEOF[int(g)]
                  in CENSUS_PENDING else '')))
with open('data/stage10w_map.csv', 'w', encoding='utf-8') as f:
    cols = ['name', 'leg', 'fgas', 'cov', 'theta', 'd_boot',
            'd_p065', 'd_gm', 'abs_share', 'sd2_share', 'flag']
    f.write(','.join(cols) + '\n')
    for r in rows:
        f.write(','.join(
            f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c])
            for c in cols) + '\n')
P(f"  map written: data/stage10w_map.csv ({len(rows)} rows)")
for leg in ('anch', 'flow'):
    sub_rows = [r for r in rows if r['leg'] == leg]
    top = sorted(sub_rows, key=lambda r: -abs(r['d_boot']))[:5]
    P(f"  {leg} top-5 carriers (boot pair): " + "; ".join(
        f"{r['name']} d = {r['d_boot']:+.1f} (sd2 {r['sd2_share']:.0%})"
        + (f" [{r['flag']}]" if r['flag'] else '') for r in top))

# ---------------- PART 2: the joint axis instrument ----------------
P("")
P("== PART 2: axes (partial-rank primary, A4-i machinery; seed 707;"
  " Bonferroni x4) ==")
gal_all = LEGA + FLOW
delta = np.array([r['d_boot'] for r in rows])
fgas = np.array([r['fgas'] for r in rows])
cov = np.array([r['cov'] for r in rows])
res = np.log10(np.array([r['theta'] for r in rows]))
legv = np.array([0.0 if r['leg'] == 'anch' else 1.0 for r in rows])
AXES = {'coverage': cov, 'composition': fgas, 'resolution': res,
        'leg': legv}
COND = {'coverage': [legv], 'composition': [cov, legv],
        'resolution': [cov, legv], 'leg': [cov]}
P(f"  collinearity (disclosed): spear(fgas,cov) = "
  f"{spear(fgas, cov):+.2f}; spear(res,cov) = {spear(res, cov):+.2f};"
  f" spear(res,leg) = {spear(res, legv):+.2f}")

rng7 = np.random.default_rng(707)
P2 = {}
for ax in ('coverage', 'composition', 'resolution', 'leg'):
    rho, p = partial_p(delta, AXES[ax], COND[ax], rng7, nperm=10000)
    P2[ax] = (rho, p)
    P(f"  {ax:11s}: partial rho = {rho:+.3f}; p = {p:.4f} "
      f"(x4 Bonferroni = {min(4*p, 1.0):.4f})")

# OLS co-read (descriptive only)
Z = np.column_stack([(v - v.mean())/v.std() for v in
                     (cov, fgas, res, legv)])
beta, *_ = np.linalg.lstsq(
    np.column_stack([np.ones(len(delta)), Z]), delta, rcond=None)
P(f"  OLS co-read (DESCRIPTIVE; collinear): std slopes cov "
  f"{beta[1]:+.1f}, fgas {beta[2]:+.1f}, res {beta[3]:+.1f}, leg "
  f"{beta[4]:+.1f}")

# carrier-drop battery
DROPS = [('UGC03580',), ('NGC5371',), ('IC2574',), ('NGC0891',),
         ('UGC03580', 'NGC5371')]
names_arr = np.array([r['name'] for r in rows])
P("  carrier-drop battery (rho; p at 2000 perms):")
DROP_P = {ax: [] for ax in P2}
for dnames in DROPS:
    keep = ~np.isin(names_arr, dnames)
    parts = []
    for ax in ('coverage', 'composition', 'resolution', 'leg'):
        rho, p = partial_p(delta[keep], AXES[ax][keep],
                           [z[keep] for z in COND[ax]], rng7,
                           nperm=2000)
        DROP_P[ax].append((rho, p))
        parts.append(f"{ax[:4]} {rho:+.2f}/{p:.3f}")
    P(f"    -{'+'.join(dnames):22s}: " + "  ".join(parts))

P("  grammar per axis (registered SS4 + A1):")
for ax in P2:
    rho, p = P2[ax]
    pb = min(4*p, 1.0)
    drops_ok = all(min(4*pp, 1.0) <= 0.05
                   and np.sign(rr) == np.sign(rho)
                   for rr, pp in DROP_P[ax])
    if pb <= 0.01 and drops_ok:
        verdict = 'ESTABLISHED'
    elif pb <= 0.05:
        verdict = 'SUPPORTED-AT-POINT-GRADE'
    else:
        verdict = 'NOT ESTABLISHED'
    P(f"    {ax:11s}: {verdict} (corrected p = {pb:.4f}; "
      f"carrier-robust = {drops_ok})")

# ---------------- PART 3: the angular-resolution cut ---------------
P("")
P("== PART 3: the cut instrument (theta_min in 10/20/30 arcsec; "
  "contest grade; SD_block recomputed per cut) ==")
CUTS = (10.0, 20.0, 30.0)
cutW, tasks, subs = {}, [], {}
for th in CUTS:
    keep = ~(PT_THETA < th)
    W2, legc, flowc = cut_world(W78, keep)
    cutW[th] = (W2, legc, flowc)
    for legnm, seq in (('anch', legc), ('flow', flowc)):
        sname = f'{legnm}{int(th)}'
        subs[sname] = build_sub(W2, seq)
        warm = WARM_A if legnm == 'anch' else WARM_F
        for m in MEMBERS:
            tasks.append(dict(key=f'{sname}:{m}', kind='contest_fit',
                              sub=sname, member=m, th0=warm[m]))
    dropA = [NAMEOF[g] for g in LEGA if g not in set(legc)]
    dropF = [NAMEOF[g] for g in FLOW if g not in set(flowc)]
    P(f"  cut {int(th)}: anchored {len(legc)}/78 gal "
      f"({sum(len(W2['gpts'][g]) for g in legc)} pts), flow "
      f"{len(flowc)}/71 gal "
      f"({sum(len(W2['gpts'][g]) for g in flowc)} pts)")
    if dropA: P(f"    dropped anch: {', '.join(dropA)}")
    if dropF: P(f"    dropped flow: {', '.join(dropF)}")
P(f"  {len(tasks)} contest tasks -> "
  f"{'pool' if POOL_OK else 'serial'}")
CRES = pool_or_serial(tasks, subs, 'cuts')

CUTRES = {}
for th in CUTS:
    for legnm in ('anch', 'flow'):
        sname = f'{legnm}{int(th)}'
        sub = subs[sname]
        r = {m: CRES[f'{sname}:{m}'] for m in MEMBERS}
        fmin = min(r[m]['fun'] for m in MEMBERS)
        d = r['boot']['fun'] - r['BE']['fun']
        A = 'BE' if r['BE']['fun'] <= r['boot']['fun'] else 'boot'
        B = 'boot' if A == 'BE' else 'BE'
        fldA = sub_fields(sub, A, r[A]['x'])
        muB = sub_fields(sub, B, r[B]['x'])[0]
        sdb = sd_block(sub, fldA, muB)
        winner = min(MEMBERS, key=lambda m: r[m]['fun'])
        CUTRES[(th, legnm)] = dict(d=d, sd=sdb, winner=winner, r=r)
        P(f"  cut {int(th):2d} {legnm}: winner = {winner}; "
          f"d(BE,boot) = {d:+.1f}; SD_block = {sdb:.1f}; deltas " +
          ", ".join(f"{m} {r[m]['fun']-fmin:+.1f}" for m in MEMBERS))
        # C15: box-pin flags
        for m in MEMBERS:
            x = r[m]['x']
            if x[2] <= 1e-4 or not (0.2 <= x[1] <= 3.0):
                P(f"    C15 FLAG {m}: s_int = {x[2]:.5f}, f_ML = "
                  f"{x[1]:.3f} (profile check owed)")

# ---------------- letter -------------------------------------------
P("")
P("== letter (registered SS5 + A1-iii) ==")
sdF = {th: CUTRES[(th, 'flow')]['sd'] for th in CUTS}
dF = {th: CUTRES[(th, 'flow')]['d'] for th in CUTS}
dA = {th: CUTRES[(th, 'anch')]['d'] for th in CUTS}
sdA = {th: CUTRES[(th, 'anch')]['sd'] for th in CUTS}

carried_pt = all((dF[th] > 0) or (abs(dF[th]) < sdF[th])
                 for th in (20.0, 30.0))
anch_stable = all(abs(dA[th] - D_UNCUT_ANCH) <= sdA[th]
                  for th in (20.0, 30.0))
immune = all(abs(dF[th] - D_UNCUT_FLOW) < sdF[th] for th in CUTS)
P(f"  flow: d_uncut = {D_UNCUT_FLOW:+.1f}; d_cut = " +
  ", ".join(f"{int(th)}: {dF[th]:+.1f} (SD {sdF[th]:.1f})"
            for th in CUTS))
P(f"  anch: d_uncut = {D_UNCUT_ANCH:+.1f}; d_cut = " +
  ", ".join(f"{int(th)}: {dA[th]:+.1f} (SD {sdA[th]:.1f})"
            for th in CUTS))
P(f"  point conditions: carried = {carried_pt} (flow sign-flip or "
  f"sub-SD at 20 AND 30); anch-stable = {anch_stable}; "
  f"immune = {immune}")

LETTER = None
if carried_pt and anch_stable:
    LETTER = 'W-RESOLUTION-CARRIED (point grade)'
elif immune:
    LETTER = 'W-RESOLUTION-IMMUNE'
else:
    LETTER = 'W-MIXED'

# conditional boot (fires only on the carried point conditions)
if carried_pt and anch_stable:
    P("")
    P("  conditional co-requirement boot: 200 paired reps, seed 909,"
      " 20-arcsec cut flow leg (flowboot construction)")
    W2, legc, flowc = cutW[20.0]
    warm20 = {m: CRES[f'flow20:{m}']['x'] for m in MEMBERS}
    rng9 = np.random.default_rng(909)
    bsubs, btasks = {}, []
    for rep in range(200):
        pick = rng9.choice(flowc, size=len(flowc), replace=True)
        docc = []
        for g in pick:
            s_g = 1.0 + rng9.normal(0, W78['rel'][int(g)])
            docc.append(-math.log10(min(max(s_g, 0.5), 1.5)))
        bsubs[f'r{rep}'] = build_sub(W2, [int(g) for g in pick],
                                     dlg_per_occ=docc)
        for m in ('BE', 'boot'):
            btasks.append(dict(key=f'r{rep}:{m}', kind='fit_deep',
                               sub=f'r{rep}', member=m,
                               th0=warm20[m]))
    BRES = pool_or_serial(btasks, bsubs, 'cutboot')
    drep = np.array([BRES[f'r{i}:boot']['fun']
                     - BRES[f'r{i}:BE']['fun'] for i in range(200)])
    p_sign = float(np.mean(drep > 0))
    p_nore = float(np.mean(drep > -sdF[20.0]))
    P(f"  boot: P(d_rep > 0) = {p_sign:.3f}; "
      f"P(d_rep > -SD_block(20)) = {p_nore:.3f} (C5 bar 0.95)")
    p_used = p_sign if dF[20.0] > 0 else p_nore
    if p_used >= 0.95:
        LETTER = 'W-RESOLUTION-CARRIED (ATTRIBUTED; C5 co-req met)'
    else:
        LETTER = ('W-RESOLUTION-CARRIED (point grade only; C5 '
                  'co-requirement NOT met)')

P("")
P(f"LETTER: {LETTER}")
P("credences: every outcome pre-signed HOLD 53 / HOLD 8 (SS8); the "
  "meter (band 65.4-70.4, stat 5.0) untouched by construction")
P(f"fits pooled = {POOL_OK}; wall-clock: "
  f"{(time.time()-t00)/60:.1f} min")
save()
print("\nsaved:", OUTFILE)
