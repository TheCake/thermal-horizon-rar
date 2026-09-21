"""STAGE 10V -- THE STRATIFICATION INSTRUMENT.

Why does the identical four-member contest prefer BE on the anchored
78 and boot on the flow 71 (BE LAST, +48.7)? Three registered axes:
H1 distance treatment / H2 composition / H3 acceleration coverage.
Pre-registration: PREREG-STRAT-10V.md (bars, seeds, thresholds and
letters frozen before any counterfactual was computed).

Engine = exec-inheritance of calcs/stage10u_nuform.py truncated at
the G10U-2 marker (the R47-addendum route): SPARC/CF4 parse, census,
make_world, FAMS, build_sub, fit_hier, fit_deep, contest_fit --
BIT-VERBATIM. No injection-calibrated bar anywhere (trap #28 moot);
margins ride the GB R-d realized-galaxy-block SD construction.

Modes: gates (sky-blind: regressions + contest-grade flow baseline
[STOP-bearing] + wiring + census; NO counterfactual) | sky (the
counterfactual battery + letter).
"""
import math, os, sys, time
import numpy as np
from scipy.optimize import minimize_scalar as msc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
MODE = sys.argv[1] if len(sys.argv) > 1 else 'gates'
SKY = (MODE == 'sky')
OUTFILE = ('data/stage10v_skyread.txt' if SKY
           else 'data/stage10v_gates.txt')

# ---------------- engine inheritance (R47-addendum route) ----------
SRC = open(os.path.join(HERE, 'stage10u_nuform.py'),
           encoding='utf-8').read()
cut = SRC.index("# ================= G10U-2: member identity")
NS = {'__name__': 'stage10u_trunc',
      '__file__': os.path.join(HERE, 'stage10u_nuform.py')}
_argv = sys.argv
sys.argv = ['stage10u_nuform.py', 'gates']
exec(compile(SRC[:cut], 'stage10u_trunc', 'exec'), NS)
sys.argv = _argv
NS['OUTFILE'] = os.path.join(ROOT, 'data', 'stage10v_scratch.txt')

W78, LEGA, FAMS = NS['W78'], NS['LEGA'], NS['FAMS']
build_sub, fit_hier, fit_deep, contest_fit = (
    NS['build_sub'], NS['fit_hier'], NS['fit_deep'], NS['contest_fit'])
nu_be = NS['nu_be']
MEMBERS, MIDX, ARC = NS['MEMBERS'], NS['MIDX'], NS['ARC']
S_ML, U_PRIOR, LN10 = NS['S_ML'], NS['U_PRIOR'], NS['LN10']
h0_of_a0 = NS['h0_of_a0']
rc_names = NS['rc_names']
FLOW = [int(g) for g in W78['flow']]
LEGA = [int(g) for g in LEGA]
NAMEOF = {int(g): W78['names'][int(g)] for g in LEGA + FLOW}

Lst = []
def P(s=""):
    print(s, flush=True)
    Lst.append(s)

def save():
    with open(OUTFILE, 'w', encoding='utf-8') as f:
        f.write("\n".join(Lst) + "\n")

t00 = time.time()
NFIT = [0]
P(f"STAGE 10V STRATIFICATION -- mode = {MODE} "
  f"({'sky read' if SKY else 'sky-blind gates'})")
P("")

# ---------------- counterfactual world builders --------------------
def sv_swap(W, gset, rel_new):
    """Prereg SS6 formula: strip the stored rel from sigv, install
    rel_new, keep the inclination term."""
    W2 = dict(W)
    W2['sigv'] = dict(W['sigv'])
    for g in gset:
        g = int(g)
        sv, rel = W['sigv'][g], W['rel'][g]
        incl = math.sqrt(max(sv*sv - (rel/LN10)**2, 0.0))
        W2['sigv'][g] = max(math.hypot(rel_new/LN10, incl), 0.01)
    return W2

def rescale_world(W, gset, s):
    """Pin-3 coherent distance rescale: lgobs -= log10 s on the set."""
    W2 = dict(W)
    W2['lgobs'] = W['lgobs'].copy()
    for g in gset:
        W2['lgobs'][W['gpts'][int(g)]] -= math.log10(s)
    return W2

def window_world(W, keep):
    """Point-filtered world; galaxies keep >= 3 surviving points."""
    idx = np.where(keep)[0]
    W2 = dict(W)
    for k in ('gg', 'gd', 'gb', 'lgobs', 'sig2', 'gal_id'):
        W2[k] = W[k][idx]
    ug = np.unique(W2['gal_id'])
    gp = {int(g): np.where(W2['gal_id'] == g)[0] for g in ug}
    keepg = {g for g in gp if len(gp[g]) >= 3}
    W2['gpts'] = {g: gp[g] for g in keepg}
    W2['ug'] = np.array(sorted(keepg))
    legw = [g for g in LEGA if g in keepg]
    floww = [g for g in FLOW if g in keepg]
    return W2, legw, floww

# ---------------- fields / block SD / rho (GB R-d verbatim) --------
def sub_fields(sub, nm, th):
    la0, f, s_int, u = th
    a0 = 10**la0
    nu = FAMS[nm]
    npt = len(sub['lg'])
    MU = np.empty(npt); RW = np.empty(npt); RES = np.empty(npt)
    for k in range(sub['n']):
        pts = sub['gpts'][k]
        lg = sub['lg'][pts]
        gg = sub['gg'][pts]; gd = sub['gd'][pts]; gb = sub['gb'][pts]
        s2 = sub['s2'][pts]
        sv = sub['sv'][k]
        isu = sub['isuma_g'][k]
        se2 = s2 + s_int*s_int
        dml, dv = 0.0, 0.0
        for _ in range(200):
            gN = gg + f*math.exp(dml)*gd + gb
            r0 = lg - np.log10(gN*nu(gN/a0)) - u*isu
            w = 1.0/se2
            dvn = float(np.sum(w*r0)/(np.sum(w) + 1.0/sv**2))

            def od(dl):
                gN2 = gg + f*math.exp(dl)*gd + gb
                rr = lg - np.log10(gN2*nu(gN2/a0)) - dvn - u*isu
                return float(np.sum(rr*rr/se2) + dl*dl/(S_ML*S_ML))
            dmln = msc(od, bounds=(-0.7, 0.7), method='bounded').x
            if abs(dmln - dml) < 1e-12 and abs(dvn - dv) < 1e-12:
                dml, dv = dmln, dvn
                break
            dml, dv = dmln, dvn
        gN = gg + f*math.exp(dml)*gd + gb
        mu = np.log10(gN*nu(gN/a0)) + dv + u*isu
        MU[pts] = mu; RW[pts] = 1.0/se2; RES[pts] = lg - mu
    return MU, RW, RES

def sd_block(sub, fldA, muB):
    muA, wA, resA = fldA
    D = muB - muA
    v = 0.0
    for k in range(sub['n']):
        pts = sub['gpts'][k]
        v += 4.0*float(np.sum(wA[pts]*resA[pts]*D[pts]))**2
    return math.sqrt(v)

def lag1(sub, fld):
    res = fld[2]
    num = den = 0.0
    for k in range(sub['n']):
        r = res[sub['gpts'][k]]
        if len(r) < 3: continue
        r = r - r.mean()
        num += float(np.sum(r[:-1]*r[1:]))
        den += float(np.sum(r*r))
    return num/den if den > 0 else float('nan')

# ---------------- contest + evaluation -----------------------------
def run_contest(sub, warm_map, label, grade='contest'):
    fits, gaps = {}, {}
    for m in MEMBERS:
        if grade == 'contest':
            bf, gap = contest_fit(sub, m, warm_map[m])
        else:
            bf = fit_deep(sub, FAMS[m], th0=warm_map[m])
            gap = 0.0
        fits[m] = bf; gaps[m] = gap
        NFIT[0] += 2 if grade == 'contest' else 1
        P(f"    [{label}] {m} fit  [{time.time()-t00:.0f}s]")
    fmin = min(fits[m].fun for m in MEMBERS)
    P(f"  {label} ({grade}; n_gal = {sub['n']}):")
    P("    member  a0          f_ML   s_int   u        d      gap")
    for m in MEMBERS:
        x = fits[m].x
        P(f"    {m:5s}  {10**x[0]:.4e}  {x[1]:.3f}  {x[2]:.4f}  "
          f"{x[3]:+.4f}  {fits[m].fun-fmin:6.1f}  {gaps[m]:5.2f}")
    return fits, gaps

def evaluate(sub, fits, gaps, incumbent, cross, label,
             base_d=None, sd_base=None):
    """Prereg SS4: winner, decisive-pair margin at block grade,
    REVERSAL / MATERIAL. d_pair convention: fun_boot - fun_BE
    (positive = BE better)."""
    funs = {m: fits[m].fun for m in MEMBERS}
    winner = min(MEMBERS, key=lambda m: funs[m])
    gapmax = max(gaps.values())
    d_pair = funs['boot'] - funs['BE']
    A = 'BE' if funs['BE'] <= funs['boot'] else 'boot'
    B = 'boot' if A == 'BE' else 'BE'
    fldA = sub_fields(sub, A, list(fits[A].x))
    muB = sub_fields(sub, B, list(fits[B].x))[0]
    sdb = sd_block(sub, fldA, muB)
    lead = funs[incumbent] - funs[cross]
    bar = max(sdb, gapmax)
    reversal = (winner == cross) and (lead >= bar)
    material = False
    if base_d is not None and not reversal:
        material = abs(d_pair - base_d) >= 2.0*sd_base
    tag = ('REVERSAL' if reversal else
           'MATERIAL' if material else
           ('winner-neither' if winner not in (incumbent, cross)
            else 'no-fire'))
    P(f"    -> winner = {winner}; d(BE,boot) = {d_pair:+.1f}; "
      f"SD_block = {sdb:.1f}; cross-lead = {lead:+.1f} vs bar "
      f"{bar:.1f}; rho = {lag1(sub, fldA):+.2f} -> {tag}")
    return dict(winner=winner, reversal=reversal, material=material,
                d_pair=d_pair, sd=sdb, gap=gapmax, funs=funs)

# ---------------- baselines (both modes) ---------------------------
P("== baselines (deterministic archived conventions) ==")
subP = build_sub(W78, LEGA)
b_plain = fit_hier(subP, nu_be, use_u=False)
b_loose = fit_hier(subP, nu_be, use_u=True, th0=list(b_plain.x)+[0.0])
arch = {'BE': fit_deep(subP, nu_be, th0=list(b_loose.x))}
for nm in ('p065', 'gm', 'boot'):
    arch[nm] = fit_deep(subP, FAMS[nm], th0=list(arch['BE'].x))
NFIT[0] += 6
P(f"  anchored refits done  [{time.time()-t00:.0f}s]")

REF_A = {'gm': 24.0, 'p065': 25.9, 'boot': 63.3}
RB_REF = {'BE': (1.050, 0.0324, -0.0211), 'p065': (1.172, 0.0331, -0.0208),
          'gm': (1.187, 0.0331, -0.0211), 'boot': (1.266, 0.0341, -0.0215)}
g1_ok = True
P("  G10V-1 anchored regression:")
for nm in MEMBERS:
    a0 = 10**arch[nm].x[0]
    da = abs(a0 - ARC[nm])/ARC[nm]
    dd = 0.0 if nm == 'BE' else abs(
        (arch[nm].fun - arch['BE'].fun) - REF_A[nm])
    f_, s_, u_ = arch[nm].x[1], arch[nm].x[2], arch[nm].x[3]
    rf, rs, ru = RB_REF[nm]
    ok = (da <= 1e-3 and dd <= 0.5 and abs(f_-rf) <= 0.01
          and abs(s_-rs) <= 0.002 and abs(u_-ru) <= 0.002)
    g1_ok &= ok
    P(f"    {nm:5s}: a0 d = {100*da:.3f}% (bar 0.1%); |d - arch| = "
      f"{dd:.2f} (bar 0.5); nuis {'ok' if ok else 'OFF'}")
P(f"  G10V-1: {'PASS' if g1_ok else 'FAIL'}")
if not g1_ok:
    P("STOP: engine regression failed (vetoes everything)")
    save(); sys.exit(0)

subF = build_sub(W78, FLOW)
fw = {}
for m in MEMBERS:
    fw[m] = fit_deep(subF, FAMS[m],
                     th0=[math.log10(ARC[m]), 1.0, 0.08, 0.0])
    NFIT[0] += 1
P(f"  flow warm refits done  [{time.time()-t00:.0f}s]")
REF_F = {'BE': 48.7, 'p065': 7.6, 'gm': 2.5, 'boot': 0.0}
fmin_w = min(fw[m].fun for m in MEMBERS)
g2a_ok = True
P("  G10V-2a flow warm replication (vs GB R-i):")
for m in MEMBERS:
    d = fw[m].fun - fmin_w
    ok = abs(d - REF_F[m]) <= 1.0
    g2a_ok &= ok
    P(f"    {m:5s}: d = {d:+.1f} vs GB {REF_F[m]:+.1f} "
      f"-> {'ok' if ok else 'OFF'}")
P(f"  G10V-2a: {'PASS' if g2a_ok else 'FAIL'}")
if not g2a_ok:
    P("STOP: flow replication failed (vetoes all flow clauses)")
    save(); sys.exit(0)

flow_fit, flow_gap = {}, {}
for m in MEMBERS:
    bc = fit_deep(subF, FAMS[m], th0=None)
    NFIT[0] += 1
    flow_fit[m] = bc if bc.fun < fw[m].fun else fw[m]
    flow_gap[m] = abs(bc.fun - fw[m].fun)
    P(f"    flow cold {m} fit  [{time.time()-t00:.0f}s]")
fminF = min(flow_fit[m].fun for m in MEMBERS)
bestF = min(MEMBERS, key=lambda m: flow_fit[m].fun)
d_be_F = flow_fit['BE'].fun - fminF
P("  G10V-2b flow contest-grade baseline:")
P("    member  a0          f_ML   s_int   u        d      gap")
for m in MEMBERS:
    x = flow_fit[m].x
    P(f"    {m:5s}  {10**x[0]:.4e}  {x[1]:.3f}  {x[2]:.4f}  "
      f"{x[3]:+.4f}  {flow_fit[m].fun-fminF:6.1f}  {flow_gap[m]:5.2f}")
g2b_ok = (bestF == 'boot') and (d_be_F >= 20.0)
P(f"  G10V-2b: best = {bestF}, d(BE) = {d_be_F:+.1f} "
  f"-> {'PASS' if g2b_ok else 'FAIL'}")
if not g2b_ok:
    P("")
    P("== LETTER ==")
    P("V-BASELINE-DISSOLVED: the banked cross-leg inversion does not")
    P("  survive the contest-grade estimator (pre-registered STOP;")
    P("  this is itself the finding and goes to review). No")
    P("  counterfactual was computed.")
    save(); sys.exit(0)

D_PAIR_A = arch['boot'].fun - arch['BE'].fun
D_PAIR_F = flow_fit['boot'].fun - flow_fit['BE'].fun

# ---------------- G10V-3 wiring ------------------------------------
P("")
P("-- G10V-3 wiring --")
Wr1 = rescale_world(W78, FLOW, 1.0)
w3a = np.array_equal(Wr1['lgobs'], W78['lgobs'])
S_RS = 73.0/65.4
Wrs = rescale_world(W78, FLOW, S_RS)
fpts = np.concatenate([W78['gpts'][g] for g in FLOW])
apts = np.concatenate([W78['gpts'][g] for g in LEGA])
sh = Wrs['lgobs'][fpts] - W78['lgobs'][fpts]
w3b = (np.allclose(sh, -math.log10(S_RS), atol=1e-12)
       and np.array_equal(Wrs['lgobs'][apts], W78['lgobs'][apts]))
w3c = True
for g in LEGA[:20]:
    W1 = sv_swap(W78, [g], W78['rel'][g])
    w3c &= abs(W1['sigv'][g] - max(W78['sigv'][g], 0.01)) <= 1e-12
sub0 = build_sub(W78, LEGA, dlg_per_occ=[0.0]*len(LEGA))
w3d = np.array_equal(sub0['lg'], subP['lg'])
g3_ok = w3a and w3b and w3c and w3d
P(f"  s=1 identity {w3a}; s={S_RS:.4f} exact-shift {w3b}; "
  f"sv self-swap {w3c}; zero-injection identity {w3d}")
P(f"  G10V-3: {'PASS' if g3_ok else 'FAIL'}")
if not g3_ok:
    P("STOP: wiring failed (vetoes H1 clauses; deterministic, "
      "identical in both modes)")
    save(); sys.exit(0)

# ---------------- G10V-5 block-SD regression -----------------------
P("")
P("-- G10V-5 block-SD regression --")
fldA_a = sub_fields(subP, 'BE', list(arch['BE'].x))
muB_a = sub_fields(subP, 'boot', list(arch['boot'].x))[0]
SD_A = sd_block(subP, fldA_a, muB_a)
g5_ok = abs(SD_A - 22.96)/22.96 <= 0.02
P(f"  anchored SD_block(BE->boot) = {SD_A:.2f} vs GB 22.96 "
  f"(bar 2%) -> {'PASS' if g5_ok else 'FAIL'}")
fldA_f = sub_fields(subF, 'boot', list(flow_fit['boot'].x))
muB_f = sub_fields(subF, 'BE', list(flow_fit['BE'].x))[0]
SD_F = sd_block(subF, fldA_f, muB_f)
P(f"  flow SD_block(boot->BE) = {SD_F:.2f} (FROZEN at this print; "
  f"no bar -- new number)")
P(f"  flow lag-1 rho (boot fields) = {lag1(subF, fldA_f):+.3f}")
if not g5_ok:
    P("G10V-5 FAIL: margin rules vetoed; letter capped at "
      "V-UNRESOLVED with point-delta disclosures")

# ---------------- G10V-6 census ------------------------------------
P("")
P("-- G10V-6 census (values FROZEN at this print) --")
REL_MED_FLOW = float(np.median([W78['rel'][g] for g in FLOW]))
REL_MED_ANCH = float(np.median([W78['rel'][g] for g in LEGA]))
SIG_INJ = REL_MED_FLOW/LN10
P(f"  rel_med_flow = {REL_MED_FLOW:.4f}; rel_med_anch = "
  f"{REL_MED_ANCH:.4f}; sigma_inj = {SIG_INJ:.4f} dex")

def fgas_of(g):
    pts = W78['gpts'][int(g)]
    den = W78['gg'][pts] + W78['gd'][pts] + W78['gb'][pts]
    return float(np.median(W78['gg'][pts]/den))

FGAS = {int(g): fgas_of(g) for g in LEGA + FLOW}
STRATA = dict(
    flow_gasrich=[g for g in FLOW if FGAS[g] >= 0.5],
    flow_gaspoor=[g for g in FLOW if FGAS[g] < 0.5],
    anch_gasrich=[g for g in LEGA if FGAS[g] >= 0.5],
    anch_gaspoor=[g for g in LEGA if FGAS[g] < 0.5])
for k, v in STRATA.items():
    P(f"  {k}: {len(v)} galaxies "
      f"{'(UNPOWERED < 15)' if len(v) < 15 else ''}")
H2_POWERED = (len(STRATA['flow_gasrich']) >= 15
              and len(STRATA['flow_gaspoor']) >= 15)

gN1_all = W78['gg'] + 1.0*W78['gd'] + W78['gb']
ly = np.log10(gN1_all/ARC['BE'])
lyA, lyF = ly[apts], ly[fpts]
WLO = max(np.percentile(lyA, 10), np.percentile(lyF, 10))
WHI = min(np.percentile(lyA, 90), np.percentile(lyF, 90))
WMID = 0.5*(WLO + WHI)
P(f"  log10 y percentiles: legA P10/P90 = {np.percentile(lyA,10):.2f}"
  f"/{np.percentile(lyA,90):.2f}; flow = {np.percentile(lyF,10):.2f}"
  f"/{np.percentile(lyF,90):.2f}")
P(f"  window = [{WLO:.2f}, {WHI:.2f}] in log10 y (y = gN1/ARC_BE)")
keepmask = (ly >= WLO) & (ly <= WHI)
WW, LEGW, FLOWW = window_world(W78, keepmask)
P(f"  windowed legs: anchored {len(LEGW)}/78, flow {len(FLOWW)}/71 "
  f"galaxies survive (>= 3 in-window points)")
dropA = [NAMEOF[g] for g in LEGA if g not in set(LEGW)]
P(f"  dropped anchored: {', '.join(dropA) if dropA else 'none'}")
H3_POWERED = (len(LEGW) >= 15 and len(FLOWW) >= 15)
if not H2_POWERED:
    P("  H2 clause: UNPOWERED (min-15 rule)")
if not H3_POWERED:
    P("  H3 clause: UNPOWERED (min-15 rule)")

agree = 0
for g in FLOW:
    gaspoor = FGAS[g] < 0.5
    hy = float(np.median(ly[W78['gpts'][g]])) > WMID
    agree += int(gaspoor == hy)
A_CONF = agree/len(FLOW)
P(f"  a_conf (flow: gas-poor vs high-y agreement) = {A_CONF:.2f} "
  f"{'-> CONFOUNDED suffix armed' if A_CONF >= 2/3 else ''}")

t_fit = (time.time() - t00)/max(NFIT[0], 1)
P(f"  measured cost: {NFIT[0]} fits so far, {t_fit:.0f} s/fit "
  f"(alternation fields included) -> sky ETA ~ "
  f"{158*t_fit/60:.0f} min for ~158 fits")

# ---------------- gate record --------------------------------------
P("")
P(f"GATES: G10V-1 PASS  G10V-2a PASS  G10V-2b PASS (best = {bestF}, "
  f"BE {d_be_F:+.1f})  G10V-3 {'PASS' if g3_ok else 'FAIL'}  "
  f"G10V-5 {'PASS' if g5_ok else 'FAIL'} (SD_A {SD_A:.1f}, SD_F "
  f"{SD_F:.1f})  G10V-6 censused (H2 "
  f"{'powered' if H2_POWERED else 'UNPOWERED'}, H3 "
  f"{'powered' if H3_POWERED else 'UNPOWERED'}, a_conf {A_CONF:.2f})")
if not SKY:
    P("")
    P("gates mode complete -- NO counterfactual contest was computed;")
    P("all census values above are FROZEN by the pre-sky commit.")
    P(f"wall-clock: {(time.time()-t00)/60:.1f} min")
    save()
    print("\nsaved:", OUTFILE)
    sys.exit(0)

# ================= SKY: the counterfactual battery ==================
WARM_A = {m: list(arch[m].x) for m in MEMBERS}
WARM_F = {m: list(flow_fit[m].x) for m in MEMBERS}
RES = {}

P("")
P("== H1 DISTANCE ==")
P("")
P("-- H1c anchored-loosened (PRIMARY; rel -> rel_med_flow, "
  "data unchanged) --")
W_loose = sv_swap(W78, LEGA, REL_MED_FLOW)
sub_h1c = build_sub(W_loose, LEGA)
f_h1c, g_h1c = run_contest(sub_h1c, WARM_A, 'H1c')
RES['H1c'] = evaluate(sub_h1c, f_h1c, g_h1c, 'BE', 'boot', 'H1c',
                      base_d=D_PAIR_A, sd_base=SD_A)

P("")
P("-- H1b(ii) anchored + flow-grade scatter + loosened priors "
  "(PRIMARY; 4 seeds) --")
h1b2 = []
for seed in (7001, 7002, 7003, 7004):
    rng = np.random.default_rng(seed)
    dlg = [float(rng.normal(0, SIG_INJ)) for _ in LEGA]
    subI = build_sub(W_loose, LEGA, dlg_per_occ=dlg)
    fi, gi = run_contest(subI, WARM_A, f'H1b2-s{seed}')
    r = evaluate(subI, fi, gi, 'BE', 'boot', f'H1b2-s{seed}',
                 base_d=D_PAIR_A, sd_base=SD_A)
    P(f"    [trap-#28 vetting] refit s_int: " +
      ", ".join(f"{m} {fi[m].x[2]:.4f}" for m in MEMBERS))
    h1b2.append(r)
RES['H1b2'] = h1b2

P("")
P("-- H1b(i) anchored + scatter, priors UNCHANGED (disclosure; "
  "2 seeds) --")
h1b1 = []
for seed in (7001, 7002):
    rng = np.random.default_rng(seed)
    dlg = [float(rng.normal(0, SIG_INJ)) for _ in LEGA]
    subI = build_sub(W78, LEGA, dlg_per_occ=dlg)
    fi, gi = run_contest(subI, WARM_A, f'H1b1-s{seed}')
    r = evaluate(subI, fi, gi, 'BE', 'boot', f'H1b1-s{seed}',
                 base_d=D_PAIR_A, sd_base=SD_A)
    P(f"    [trap-#28 vetting] refit s_int: " +
      ", ".join(f"{m} {fi[m].x[2]:.4f}" for m in MEMBERS))
    h1b1.append(r)
RES['H1b1'] = h1b1

P("")
P("-- H1a flow-rescale (CO-READ; coherent D x 73/H0') --")
h1a = []
for s in (73.0/65.4, 73.0/70.4):
    Wr = rescale_world(W78, FLOW, s)
    subR = build_sub(Wr, FLOW)
    fr, gr = run_contest(subR, WARM_F, f'H1a-s{s:.4f}')
    h1a.append(evaluate(subR, fr, gr, 'boot', 'BE', f'H1a-s{s:.4f}',
                        base_d=D_PAIR_F, sd_base=SD_F))
RES['H1a'] = h1a

P("")
P("-- H1d flow-tightened (CO-READ; rel -> rel_med_anch, "
  "data unchanged) --")
W_tight = sv_swap(W78, FLOW, REL_MED_ANCH)
sub_h1d = build_sub(W_tight, FLOW)
f_h1d, g_h1d = run_contest(sub_h1d, WARM_F, 'H1d')
RES['H1d'] = evaluate(sub_h1d, f_h1d, g_h1d, 'boot', 'BE', 'H1d',
                      base_d=D_PAIR_F, sd_base=SD_F)

h1_fire = RES['H1c']['reversal'] or \
    (sum(r['reversal'] for r in RES['H1b2']) >= 3)
P(f"  H1 FIRES: {h1_fire} (H1c reversal {RES['H1c']['reversal']}; "
  f"H1b(ii) reversals {sum(r['reversal'] for r in RES['H1b2'])}/4)")

P("")
P("== H2 COMPOSITION (gas dominance; fgas >= 0.5) ==")
h2 = {}
for k in ('flow_gasrich', 'flow_gaspoor', 'anch_gasrich',
          'anch_gaspoor'):
    gl = STRATA[k]
    if len(gl) < 3:
        P(f"  {k}: only {len(gl)} galaxies -- skipped")
        continue
    warm = WARM_F if k.startswith('flow') else WARM_A
    inc, cro = (('boot', 'BE') if k.startswith('flow')
                else ('BE', 'boot'))
    subS = build_sub(W78, gl)
    fs, gs = run_contest(subS, warm, k)
    h2[k] = evaluate(subS, fs, gs, inc, cro, k)
RES['H2'] = h2
h2_fire = False
if H2_POWERED and 'flow_gasrich' in h2 and 'flow_gaspoor' in h2:
    r1, r2 = h2['flow_gasrich'], h2['flow_gaspoor']
    h2_fire = (r1['winner'] == 'BE'
               and (r1['funs']['boot'] - r1['funs']['BE'])
               >= max(r1['sd'], r1['gap'])
               and r2['winner'] == 'boot'
               and (r2['funs']['BE'] - r2['funs']['boot'])
               >= max(r2['sd'], r2['gap']))
P(f"  H2 FIRES: {h2_fire} "
  f"{'' if H2_POWERED else '(UNPOWERED by census)'}")

P("")
P("== H3 COVERAGE (window-matched legs) ==")
h3 = {}
if H3_POWERED:
    subWA = build_sub(WW, LEGW)
    fa, ga = run_contest(subWA, WARM_A, 'H3-anchored-windowed')
    h3['anch'] = evaluate(subWA, fa, ga, 'BE', 'boot',
                          'H3-anchored-windowed')
    subWF = build_sub(WW, FLOWW)
    ff, gf = run_contest(subWF, WARM_F, 'H3-flow-windowed')
    h3['flow'] = evaluate(subWF, ff, gf, 'boot', 'BE',
                          'H3-flow-windowed')
    agree_w = h3['anch']['winner'] == h3['flow']['winner']
    mA = h3['anch']; mF = h3['flow']
    other = {'BE': 'boot', 'boot': 'BE'}
    okA = mA['winner'] in other and \
        (mA['funs'][other[mA['winner']]] - mA['funs'][mA['winner']]
         ) >= max(mA['sd'], mA['gap'])
    okF = mF['winner'] in other and \
        (mF['funs'][other[mF['winner']]] - mF['funs'][mF['winner']]
         ) >= max(mF['sd'], mF['gap'])
    h3_fire = agree_w and okA and okF
else:
    h3_fire = False
    P("  UNPOWERED by census -- not run")
RES['H3'] = h3
P(f"  H3 FIRES: {h3_fire}")

P("")
P("== SCREENING (warm-only; no letter weight) ==")
RC11 = [g for g in LEGA if NAMEOF[g] in rc_names
        and NAMEOF[g] != NS['ZERO_PT']]
UMA26 = [g for g in LEGA if W78['fd'][g] == 4]
AN41 = [g for g in LEGA if g not in set(RC11) | set(UMA26)]
for k, gl in (('anchored-direct-41', AN41), ('UMa-26', UMA26),
              ('reclassified-11', RC11)):
    subS = build_sub(W78, gl)
    fs, gs = run_contest(subS, WARM_A, k, grade='screening')
    fm = min(fs[m].fun for m in MEMBERS)
    P(f"    -> winner = {min(MEMBERS, key=lambda m: fs[m].fun)} "
      f"(screening grade)")

P("")
P("-- flow jackknife (standing; top-2 fixed-parameter carriers) --")
contribF = {}
muA_f, wA_f, resA_f = fldA_f
for g in FLOW:
    pts_g = subF['gpts'][FLOW.index(g)]
    DD = muB_f[pts_g] - muA_f[pts_g]
    contribF[NAMEOF[g]] = float(
        np.sum(wA_f[pts_g]*(resA_f[pts_g] - DD)**2)
        - np.sum(wA_f[pts_g]*resA_f[pts_g]**2))
top2 = sorted(contribF.items(), key=lambda kv: -kv[1])[:2]
P("  fixed-parameter top carriers of d(boot,BE) on flow "
  "(CONSTRUCTION-SENSITIVE): " +
  ", ".join(f"{n} {v:+.1f}" for n, v in top2))
for nmgal, _v in top2:
    gi = [g for g in FLOW if NAMEOF[g] == nmgal][0]
    gl = [g for g in FLOW if g != gi]
    sub2 = build_sub(W78, gl)
    f2 = {}
    for m in ('BE', 'boot'):
        f2[m] = fit_deep(sub2, FAMS[m], th0=WARM_F[m]).fun
        NFIT[0] += 1
    dj = (flow_fit['BE'].fun - flow_fit['boot'].fun) - \
        (f2['BE'] - f2['boot'])
    P(f"  drop {nmgal}: d(boot beats BE) changes by {dj:+.1f} "
      f"of {-D_PAIR_F:+.1f}")

# ---------------- letter -------------------------------------------
P("")
P("== LETTER ==")
fired = []
if h1_fire: fired.append('DISTANCE')
if h2_fire: fired.append('COMPOSITION')
if h3_fire: fired.append('COVERAGE')
conf_sfx = ('-CONFOUNDED' if A_CONF >= 2/3 and
            any(x in fired for x in ('COMPOSITION', 'COVERAGE'))
            else '')
if not g5_ok:
    letter = 'V-UNRESOLVED (margin rules vetoed by G10V-5)'
elif len(fired) == 0:
    letter = 'V-UNRESOLVED'
elif len(fired) == 1:
    letter = f'V-{fired[0]}-CARRIED{conf_sfx}'
else:
    letter = f"V-MIXED({'+'.join(fired)}){conf_sfx}"
P(f"{letter}")
P("")
P("MANDATORY DISCLOSURES:")
P("  multiplicity: 4 fire opportunities (H1c, H1b(ii), H2, H3);")
P("  attribution-grade only -- NO meter quantity moves on any")
P("  branch (band 65.4-70.4, stat 5.0, all banked clauses")
P("  untouched); fixed-parameter decompositions are")
P("  construction-sensitive (refit jackknife primary); screening")
P("  rows are warm-only grade.")
mat_rows = [k for k in ('H1c', 'H1d') if RES[k].get('material')] + \
    [f'H1a-{i}' for i, r in enumerate(RES['H1a']) if r['material']] + \
    [f'H1b2-{i}' for i, r in enumerate(RES['H1b2'])
     if r['material']] + \
    [f'H1b1-{i}' for i, r in enumerate(RES['H1b1']) if r['material']]
P(f"  MATERIAL rows: {mat_rows if mat_rows else 'none'}")
P(f"  a_conf = {A_CONF:.2f}; H2 "
  f"{'powered' if H2_POWERED else 'UNPOWERED'}; H3 "
  f"{'powered' if H3_POWERED else 'UNPOWERED'}")
P("")
P("credences: NO cell moves (pre-signed SS8); 53 / 8 untouched")
P(f"total fits: {NFIT[0]}; wall-clock: {(time.time()-t00)/60:.1f} min")
save()
print("\nsaved:", OUTFILE)
