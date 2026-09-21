"""ROUND-48 ADDENDUM -- verification for stage 10V (stratification).

GA half (BLIND, committed before the review round exists; 15th
execution of the 87a4676 protocol): six independent legs against the
10V sky read, none reusing the stage's optimiser path where a fresh
implementation is possible.
GB half (post-report): re-computes every load-bearing reviewer number
before adoption (standing memory rule); filled in after the round.

Engine + builders inherited by exec of stage10v_strat.py truncated at
the baselines marker (so NO stage fit runs on import); stage10v itself
exec-inherits stage10u_nuform.py -- one chain, zero drift.
"""
import math, os, re, sys, time
import numpy as np
from scipy.optimize import minimize_scalar as msc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
MODE = sys.argv[1] if len(sys.argv) > 1 else 'ga'

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
FAMS, MEMBERS = NSV['FAMS'], NSV['MEMBERS']
build_sub, fit_deep = NSV['build_sub'], NSV['fit_deep']
sv_swap, rescale_world, window_world = (
    NSV['sv_swap'], NSV['rescale_world'], NSV['window_world'])
S_ML, U_PRIOR, LN10 = NSV['S_ML'], NSV['U_PRIOR'], NSV['LN10']
ARC = NSV['ARC']
NAMEOF = NSV['NAMEOF']

OUT = []
def Q(s=""):
    print(s, flush=True)
    OUT.append(s)

def save():
    with open('data/round48_addendum.txt', 'a', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

GATES_TXT = open('data/stage10v_gates.txt', encoding='utf-8').read()
SKY_TXT = open('data/stage10v_skyread.txt', encoding='utf-8').read()

Q("")
Q("=" * 72)
Q(f"ROUND-48 ADDENDUM -- mode = {MODE}")
Q("=" * 72)

# ---------------- shared parsers -----------------------------------
def contest_table(txt, label):
    """Parse a printed contest table -> {member: (a0,f,s,u,d,gap)}."""
    pat = (re.escape(f"  {label} (") +
           r"[^\n]*\):\n[^\n]*\n((?:\s+\w+\s+[\d.e+-]+\s+[\d.]+\s+"
           r"[\d.]+\s+[+-][\d.]+\s+[\d.-]+\s+[\d.]+\n){4})")
    m = re.search(pat, txt)
    if not m: return None
    rows = {}
    for line in m.group(1).strip().split("\n"):
        t = line.split()
        rows[t[0]] = tuple(float(v) for v in t[1:7])
    return rows

def eval_line(txt, label):
    """Parse the evaluate print following a contest block."""
    pat = (re.escape(f"  {label} (") +
           r".*?-> winner = (\w+); d\(BE,boot\) = ([+-][\d.]+); "
           r"SD_block = ([\d.]+); cross-lead = ([+-][\d.]+) vs bar "
           r"([\d.]+); rho = ([+-][\d.]+) -> ([\w-]+)")
    m = re.search(pat, txt, re.S)
    if not m: return None
    return dict(winner=m.group(1), d=float(m.group(2)),
                sd=float(m.group(3)), lead=float(m.group(4)),
                bar=float(m.group(5)), rho=float(m.group(6)),
                tag=m.group(7))

# ---------------- independent profiled objective -------------------
def indep_obj(sub, nm, th):
    """Fresh implementation (R47-GA-1 pattern generalized to a sub):
    profiled -2lnL at fixed (la0, f, s_int, u), own alternation over
    per-galaxy dml/dv to 1e-10."""
    la0, f, s_int, u = th
    a0 = 10**la0
    nu = FAMS[nm]
    tot = 0.0
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
        tot += prev
    return tot + (u/U_PRIOR)**2

def indep_fields(sub, nm, th):
    """Own-alternation per-point fields for the block-SD leg."""
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
        for _ in range(400):
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

def blk_sd(sub, fldA, muB):
    muA, wA, resA = fldA
    D = muB - muA
    v = 0.0
    for k in range(sub['n']):
        pts = sub['gpts'][k]
        v += 4.0*float(np.sum(wA[pts]*resA[pts]*D[pts]))**2
    return math.sqrt(v)

def th_of(row):
    a0, f, s, u = row[0], row[1], row[2], row[3]
    return [math.log10(a0), f, s, u]

# ---------------- census re-derivation (shared) --------------------
REL_MED_FLOW = float(np.median([W78['rel'][g] for g in FLOW]))
REL_MED_ANCH = float(np.median([W78['rel'][g] for g in LEGA]))
SIG_INJ = REL_MED_FLOW/LN10

def fgas_of(g):
    pts = W78['gpts'][int(g)]
    den = W78['gg'][pts] + W78['gd'][pts] + W78['gb'][pts]
    return float(np.median(W78['gg'][pts]/den))

FGAS = {int(g): fgas_of(g) for g in LEGA + FLOW}
gN1_all = W78['gg'] + 1.0*W78['gd'] + W78['gb']
ly = np.log10(gN1_all/ARC['BE'])
apts = np.concatenate([W78['gpts'][g] for g in LEGA])
fpts = np.concatenate([W78['gpts'][g] for g in FLOW])
WLO = max(np.percentile(ly[apts], 10), np.percentile(ly[fpts], 10))
WHI = min(np.percentile(ly[apts], 90), np.percentile(ly[fpts], 90))
WW, LEGW, FLOWW = window_world(W78, (ly >= WLO) & (ly <= WHI))

if MODE == 'ga':
    ok_all = True
    t0 = time.time()

    # ---- GA-1: independent objective at printed optima ----
    Q("")
    Q("-- GA-1 independent profiled objective at PRINTED optima --")
    ga1_ok = True
    W_loose = sv_swap(W78, LEGA, REL_MED_FLOW)
    sub_h1c = build_sub(W_loose, LEGA)
    rng = np.random.default_rng(7001)
    dlg = [float(rng.normal(0, SIG_INJ)) for _ in LEGA]
    sub_b2 = build_sub(W_loose, LEGA, dlg_per_occ=dlg)
    sub_wf = build_sub(WW, FLOWW)
    for label, sub in (('H1c', sub_h1c), ('H1b2-s7001', sub_b2),
                       ('H3-flow-windowed', sub_wf)):
        tab = contest_table(SKY_TXT, label)
        if tab is None:
            Q(f"  {label}: TABLE NOT PARSED -> MISMATCH")
            ga1_ok = False
            continue
        iv = {m: indep_obj(sub, m, th_of(tab[m])) for m in MEMBERS}
        ivmin = min(iv.values())
        ds = []
        for m in MEMBERS:
            d_ind = iv[m] - ivmin
            ok = abs(d_ind - tab[m][4]) <= 0.5
            ga1_ok &= ok
            ds.append(f"{m} {d_ind:+.1f}/{tab[m][4]:+.1f}"
                      f"{'' if ok else ' MISMATCH'}")
        Q(f"  {label}: indep/printed d = " + ", ".join(ds) +
          f"  [{time.time()-t0:.0f}s]")
    Q(f"GA-1: {'OK' if ga1_ok else 'MISMATCH'}")
    ok_all &= ga1_ok

    # ---- GA-2: independent block-SD ----
    Q("")
    Q("-- GA-2 independent block-SD at printed optima --")
    ga2_ok = True
    for label, sub in (('H1c', sub_h1c),
                       ('H3-flow-windowed', sub_wf)):
        ev = eval_line(SKY_TXT, label)
        tab = contest_table(SKY_TXT, label)
        if ev is None or tab is None:
            Q(f"  {label}: PARSE FAIL"); ga2_ok = False; continue
        funs = {m: tab[m][4] for m in MEMBERS}
        A = 'BE' if funs['BE'] <= funs['boot'] else 'boot'
        B = 'boot' if A == 'BE' else 'BE'
        fldA = indep_fields(sub, A, th_of(tab[A]))
        muB = indep_fields(sub, B, th_of(tab[B]))[0]
        sd = blk_sd(sub, fldA, muB)
        ok = abs(sd - ev['sd'])/ev['sd'] <= 0.10
        ga2_ok &= ok
        Q(f"  {label}: SD_block indep {sd:.1f} vs printed "
          f"{ev['sd']:.1f} -> {'OK' if ok else 'MISMATCH'}")
    Q(f"GA-2: {'OK' if ga2_ok else 'MISMATCH'}")
    ok_all &= ga2_ok

    # ---- GA-3: letter logic from prints ----
    Q("")
    Q("-- GA-3 letter logic re-derived from the printed record --")
    fires = dict(re.findall(r"(H[123]) FIRES: (\w+)", SKY_TXT))
    fired = []
    if fires.get('H1') == 'True': fired.append('DISTANCE')
    if fires.get('H2') == 'True': fired.append('COMPOSITION')
    if fires.get('H3') == 'True': fired.append('COVERAGE')
    aconf = float(re.search(r"a_conf .* = ([\d.]+)", GATES_TXT)
                  .group(1))
    sfx = ('-CONFOUNDED' if aconf >= 2/3 and
           any(x in fired for x in ('COMPOSITION', 'COVERAGE'))
           else '')
    if len(fired) == 0:
        lett = 'V-UNRESOLVED'
    elif len(fired) == 1:
        lett = f'V-{fired[0]}-CARRIED{sfx}'
    else:
        lett = f"V-MIXED({'+'.join(fired)}){sfx}"
    printed = re.search(r"== LETTER ==\n(\S[^\n]*)", SKY_TXT).group(1)
    ga3_ok = printed.startswith(lett)
    Q(f"  re-derived: fired = {fired}, letter = {lett}; printed = "
      f"{printed} -> {'OK' if ga3_ok else 'MISMATCH'}")
    ok_all &= ga3_ok

    # ---- GA-4: census re-derivation ----
    Q("")
    Q("-- GA-4 census re-derivation (own code) vs frozen gates --")
    gtxt = GATES_TXT
    checks = [
        ('rel_med_flow', REL_MED_FLOW,
         float(re.search(r"rel_med_flow = ([\d.]+)", gtxt).group(1)),
         1e-3),
        ('rel_med_anch', REL_MED_ANCH,
         float(re.search(r"rel_med_anch = ([\d.]+)", gtxt).group(1)),
         1e-3),
        ('sigma_inj', SIG_INJ,
         float(re.search(r"sigma_inj = ([\d.]+)", gtxt).group(1)),
         1e-3),
        ('window_lo', WLO,
         float(re.search(r"window = \[([-\d.]+), ([-\d.]+)\]",
                         gtxt).group(1)), 0.01),
        ('window_hi', WHI,
         float(re.search(r"window = \[([-\d.]+), ([-\d.]+)\]",
                         gtxt).group(2)), 0.01),
    ]
    ga4_ok = True
    for nm_, mine, prt, tol in checks:
        ok = abs(mine - prt) <= tol
        ga4_ok &= ok
        Q(f"  {nm_}: mine {mine:.4f} vs printed {prt:.4f} -> "
          f"{'OK' if ok else 'MISMATCH'}")
    cnt = dict(flow_gasrich=len([g for g in FLOW if FGAS[g] >= 0.5]),
               flow_gaspoor=len([g for g in FLOW if FGAS[g] < 0.5]),
               anch_gasrich=len([g for g in LEGA if FGAS[g] >= 0.5]),
               anch_gaspoor=len([g for g in LEGA if FGAS[g] < 0.5]))
    for k, v in cnt.items():
        prt = int(re.search(k + r": (\d+) galaxies", gtxt).group(1))
        ok = v == prt
        ga4_ok &= ok
        Q(f"  {k}: mine {v} vs printed {prt} -> "
          f"{'OK' if ok else 'MISMATCH'}")
    wcnt = re.search(r"windowed legs: anchored (\d+)/78, flow "
                     r"(\d+)/71", gtxt)
    ok = (len(LEGW) == int(wcnt.group(1))
          and len(FLOWW) == int(wcnt.group(2)))
    ga4_ok &= ok
    Q(f"  windowed counts: mine {len(LEGW)}/{len(FLOWW)} vs printed "
      f"{wcnt.group(1)}/{wcnt.group(2)} -> "
      f"{'OK' if ok else 'MISMATCH'}")
    agree = 0
    WMID = 0.5*(WLO + WHI)
    for g in FLOW:
        hy = float(np.median(ly[W78['gpts'][g]])) > WMID
        agree += int((FGAS[g] < 0.5) == hy)
    ac = agree/len(FLOW)
    ok = abs(ac - aconf) <= 0.01
    ga4_ok &= ok
    Q(f"  a_conf: mine {ac:.2f} vs printed {aconf:.2f} -> "
      f"{'OK' if ok else 'MISMATCH'}")
    Q(f"GA-4: {'OK' if ga4_ok else 'MISMATCH'}")
    ok_all &= ga4_ok

    # ---- GA-5: injection stream fingerprints ----
    Q("")
    Q("-- GA-5 H1b rng fingerprints + injection identity --")
    ga5_ok = True
    for seed in (7001, 7002, 7003, 7004):
        r = np.random.default_rng(seed)
        d3 = [float(r.normal(0, SIG_INJ)) for _ in range(3)]
        Q(f"  seed {seed}: first-3 deltas = "
          f"[{d3[0]:+.4f}, {d3[1]:+.4f}, {d3[2]:+.4f}]")
    base = build_sub(W_loose, LEGA)
    dsub = sub_b2['lg'] - base['lg']
    exp = np.concatenate([np.full(len(base['gpts'][k]), dlg[k])
                          for k in range(base['n'])])
    ok = np.allclose(dsub, exp, atol=1e-12)
    ga5_ok &= ok
    Q(f"  seed-7001 sub lg == base + per-galaxy delta exactly -> "
      f"{'OK' if ok else 'MISMATCH'}")
    Q(f"GA-5: {'OK' if ga5_ok else 'MISMATCH'}")
    ok_all &= ga5_ok

    # ---- GA-6: arithmetic pack over every evaluate line ----
    Q("")
    Q("-- GA-6 evaluate-line arithmetic (all printed contests) --")
    ga6_ok = True
    labels = re.findall(r"\n  ([\w().=-]+(?:-s[\d.]+)?) \(contest; ",
                        SKY_TXT)
    for label in labels:
        ev = eval_line(SKY_TXT, label)
        tab = contest_table(SKY_TXT, label)
        if ev is None or tab is None: continue
        funs = {m: tab[m][4] for m in MEMBERS}
        d_pair = funs['boot'] - funs['BE']
        winner = min(MEMBERS, key=lambda m: funs[m])
        gapmax = max(tab[m][5] for m in MEMBERS)
        bar = max(ev['sd'], gapmax)
        ok = (abs(d_pair - ev['d']) <= 0.15
              and winner == ev['winner']
              and abs(bar - ev['bar']) <= 0.15)
        ga6_ok &= ok
        if not ok:
            Q(f"  {label}: MISMATCH (d {d_pair:+.1f}/{ev['d']:+.1f}, "
              f"winner {winner}/{ev['winner']}, bar "
              f"{bar:.1f}/{ev['bar']:.1f})")
    Q(f"  {len(labels)} contest blocks checked")
    Q(f"GA-6: {'OK' if ga6_ok else 'MISMATCH'}")
    ok_all &= ga6_ok

    Q("")
    Q(f"GA BLIND HALF: {'ALL OK' if ok_all else 'MISMATCHES PRESENT'}")
    save()
    sys.exit(0)

# ================= GB half (post-report verification) ================
# Re-computes every load-bearing Round-48 reviewer number in own code
# (standing memory rule). Construction-sensitive quantities are graded
# AGREE-RANGE per the R47 GB convention.
import itertools
from scipy.optimize import minimize
from scipy.stats import spearmanr

contest_fit = NSV['contest_fit']
sub_fields = NSV['sub_fields']
sd_block_f = NSV['sd_block']
NS10U = NSV['NS']
GI = {NAMEOF[g]: g for g in LEGA + FLOW}
ok_all = True
t0 = time.time()

def gb(tag, ok, txt):
    global ok_all
    ok_all &= bool(ok)
    Q(f"GB {tag}: {txt} -> {'AGREE' if ok else 'DIFFER'}")

# ---- flow baseline (warm, deterministic archived convention) ----
subF = build_sub(W78, FLOW)
fF = {}
for m in ('BE', 'boot'):
    fF[m] = fit_deep(subF, FAMS[m],
                     th0=[math.log10(ARC[m]), 1.0, 0.08, 0.0])
d_flow = fF['boot'].fun - fF['BE'].fun
Q(f"  flow baseline d(BE,boot) = {d_flow:+.2f}  [{time.time()-t0:.0f}s]")

# ---- GB-1: post-drop refits ----
Q("")
Q("-- GB-1 post-drop flow values (refits) --")
ref1 = {('UGC03580',): +7.2, ('NGC5371',): -31.2,
        ('UGC03580', 'NGC5371'): +27.5}
gb1 = True
for drops, want in ref1.items():
    gl = [g for g in FLOW if NAMEOF[g] not in drops]
    s2 = build_sub(W78, gl)
    fb = {m: fit_deep(s2, FAMS[m], th0=list(fF[m].x)).fun
          for m in ('BE', 'boot')}
    d = fb['boot'] - fb['BE']
    win = 'BE' if d > 0 else 'boot'
    ok = abs(d - want) <= 1.5
    gb1 &= ok
    Q(f"  drop {'+'.join(drops)}: d(BE,boot) = {d:+.1f} (rev {want:+.1f}), "
      f"winner {win} -> {'ok' if ok else 'OFF'}")
gb('1 post-drop flips', gb1, "UGC03580 alone flips the flow leg"
   if gb1 else "post-drop values differ")

# ---- GB-2: baseline margin + SD_F carrier decomposition ----
Q("")
Q("-- GB-2 flow baseline vs margin rule + SD_F decomposition --")
fldB = sub_fields(subF, 'boot', list(fF['boot'].x))
muBE = sub_fields(subF, 'BE', list(fF['BE'].x))[0]
D = muBE - fldB[0]
sg = {}
for k in range(subF['n']):
    pts = subF['gpts'][k]
    sg[NAMEOF[FLOW[k]]] = 2.0*float(
        np.sum(fldB[1][pts]*fldB[2][pts]*D[pts]))
v = sum(s*s for s in sg.values())
SDF = math.sqrt(v)
sh1 = sg['UGC03580']**2/v
sh2 = (sg['UGC03580']**2 + sg['NGC5371']**2)/v
SDF_wo = math.sqrt(v - sg['UGC03580']**2 - sg['NGC5371']**2)
z_base = abs(d_flow)/SDF
gb('2 SD_F + shares',
   abs(SDF - 68.24) <= 1.5 and abs(sh1 - 0.52) <= 0.06
   and abs(sh2 - 0.64) <= 0.06 and abs(SDF_wo - 40.85) <= 3.0
   and abs(z_base - 0.71) <= 0.05,
   f"SD_F {SDF:.2f} (rev 68.24), UGC share {100*sh1:.0f}% (52), "
   f"top2 {100*sh2:.0f}% (64), SD w/o {SDF_wo:.1f} (40.85), "
   f"baseline z {z_base:.2f} (0.71)")

# ---- GB-3: per-galaxy fixed-global map + robust scores ----
Q("")
Q("-- GB-3 per-galaxy delta map at parent optima (both legs) --")
def pergal_map(sub, gals, thA, thB):
    fA = sub_fields(sub, 'BE', thA)
    fB = sub_fields(sub, 'boot', thB)
    out = {}
    for k in range(sub['n']):
        pts = sub['gpts'][k]
        cA = float(np.sum(fA[2][pts]**2*fA[1][pts]
                          - np.log(fA[1][pts])))
        cB = float(np.sum(fB[2][pts]**2*fB[1][pts]
                          - np.log(fB[1][pts])))
        out[NAMEOF[gals[k]]] = cB - cA
    return out
subP = build_sub(W78, LEGA)
bA = {m: fit_deep(subP, FAMS[m],
                  th0=[math.log10(ARC[m]), 1.0, 0.08, 0.0])
      for m in ('BE', 'boot')}
MAP_F = pergal_map(subF, FLOW, list(fF['BE'].x), list(fF['boot'].x))
MAP_A = pergal_map(subP, LEGA, list(bA['BE'].x), list(bA['boot'].x))
vals = np.array(list(MAP_F.values()))
med, mad = float(np.median(vals)), float(
    np.median(np.abs(vals - np.median(vals))))
srob = 1.4826*mad
zU = (MAP_F['UGC03580'] - med)/srob
zN = (MAP_F['NGC5371'] - med)/srob
opp = max(MAP_F.items(), key=lambda kv: kv[1])
n69 = sum(v for k2, v in MAP_F.items()
          if k2 not in ('UGC03580', 'NGC5371'))
gb('3 robust scores',
   zU < -30 and zN < -10 and n69 > 5,
   f"UGC {MAP_F['UGC03580']:+.1f} = {zU:+.1f} rob-sig (rev -57.4), "
   f"N5371 {zN:+.1f} (rev -19.9), max-opp {opp[0]} {opp[1]:+.1f} "
   f"(rev NGC5985 +19.7), 69-gal sum {n69:+.1f} BE-ward (rev +15.8) "
   f"[construction-sensitive: AGREE-RANGE]")

# ---- GB-2b: fixed-global galaxy bootstrap both legs ----
rngb = np.random.default_rng(9901)
def fgboot(mp):
    v2 = np.array(list(mp.values()))
    obs = v2.sum()
    n = len(v2)
    sums = np.array([v2[rngb.integers(0, n, n)].sum()
                     for _ in range(200000)])
    return obs, float(sums.std()), obs/float(sums.std())
oF, sF, zF = fgboot(MAP_F)
oA, sA, zA = fgboot(MAP_A)
gb('2b fixed-global boots',
   abs(zF) <= 1.3 and 1.6 <= abs(zA) <= 3.0,
   f"flow obs {oF:+.1f} SD {sF:.1f} z {zF:+.2f} (rev -0.78); "
   f"anch obs {oA:+.1f} SD {sA:.1f} z {zA:+.2f} (rev +2.22) "
   f"[AGREE-RANGE]")

# ---- GB-4: permutation battery on the anchored map ----
Q("")
Q("-- GB-4 label permutations (anchored map) --")
def fgas_med(g):
    pts = W78['gpts'][int(g)]
    den = W78['gg'][pts] + W78['gd'][pts] + W78['gb'][pts]
    return float(np.median(W78['gg'][pts]/den))
GAS_A = np.array([FGAS[g] >= 0.5 for g in LEGA])
LYm = np.array([float(np.median(ly[W78['gpts'][g]])) for g in LEGA])
tA = np.array([MAP_A[NAMEOF[g]] for g in LEGA])
rngp = np.random.default_rng(9902)
def perm_p(mask, t, cond=None, nrep=200000):
    obs = t[mask].sum()
    k = int(mask.sum())
    n = len(t)
    cnt = 0
    sims = np.empty(nrep)
    if cond is None:
        for i in range(nrep):
            idx = rngp.permutation(n)[:k]
            sims[i] = t[idx].sum()
    else:
        qs = np.quantile(cond, [0.25, 0.5, 0.75])
        binid = np.digitize(cond, qs)
        for i in range(nrep):
            pm = np.empty(n, dtype=bool)
            for b in range(4):
                sel = np.where(binid == b)[0]
                kk = int(mask[sel].sum())
                pick = rngp.permutation(len(sel))[:kk]
                pm[sel] = False
                pm[sel[pick]] = True
            sims[i] = t[pm].sum()
    p = float(np.mean(np.abs(sims - sims.mean())
                      >= abs(obs - sims.mean())))
    return obs, float(sims.mean()), float(sims.std()), p
o1, m1, s1, p1 = perm_p(GAS_A, tA)
o2, m2, s2v, p2 = perm_p(GAS_A, tA, cond=LYm)
deep = LYm < np.median(LYm)
o3, m3, s3, p3 = perm_p(deep, tA)
GAS_F = np.array([FGAS[g] >= 0.5 for g in FLOW])
tF = np.array([MAP_F[NAMEOF[g]] for g in FLOW])
o4, m4, s4, p4 = perm_p(GAS_F, tF)
gb('4 permutations',
   p1 <= 0.06 and 0.04 <= p2 <= 0.25 and p3 <= 0.02 and p4 >= 0.15,
   f"anch gasrich {o1:+.1f} vs {m1:+.1f}+-{s1:.1f} p={p1:.4f} "
   f"(rev 0.0144); conditioned p={p2:.3f} (rev 0.092); coverage "
   f"deep-half {o3:+.1f} vs {m3:+.1f}+-{s3:.1f} p={p3:.4f} "
   f"(rev 0.0039); flow gasrich p={p4:.3f} (rev 0.373) "
   f"[AGREE-RANGE]")

# ---- GB-5: anch_gasrich re-gradings ----
Q("")
Q(f"-- GB-5 anch_gasrich re-grades  [{time.time()-t0:.0f}s] --")
AGR = [g for g in LEGA if FGAS[g] >= 0.5]
subG = build_sub(W78, AGR)
thw = {'BE': [math.log10(9.9741e-11), 0.904, 0.0532, 0.0228],
       'boot': [math.log10(1.0618e-10), 0.990, 0.0518, 0.0229]}
fG = {m: fit_deep(subG, FAMS[m], th0=thw[m]) for m in ('BE', 'boot')}
dG = fG['boot'].fun - fG['BE'].fun
fldGB = sub_fields(subG, 'boot', list(fG['boot'].x))
muGBE = sub_fields(subG, 'BE', list(fG['BE'].x))[0]
DG = muGBE - fldGB[0]
sgG = np.array([2.0*float(np.sum(
    fldGB[1][subG['gpts'][k]]*fldGB[2][subG['gpts'][k]]
    * DG[subG['gpts'][k]])) for k in range(16)])
T0 = sgG.sum()
cnt = 0
for signs in itertools.product((1, -1), repeat=16):
    if abs(sum(s*x for s, x in zip(signs, sgG))) >= abs(T0) - 1e-12:
        cnt += 1
p_sf = cnt/65536.0
loo = []
for k in range(16):
    gl = [g for i, g in enumerate(AGR) if i != k]
    s2 = build_sub(W78, gl)
    f2 = {m: fit_deep(s2, FAMS[m], th0=list(fG[m].x)).fun
          for m in ('BE', 'boot')}
    loo.append(f2['boot'] - f2['BE'])
loo = np.array(loo)
se_j = math.sqrt((16 - 1)/16.0*np.sum((loo - loo.mean())**2))
zj = abs(dG)/se_j
rng5 = np.random.default_rng(202)
rel_sig = {g: (2.3/18.0 if W78['fd'][g] == 4 else W78['rel'][g])
           for g in AGR}
uma_f = {g: W78['fd'][g] == 4 for g in AGR}
wins = 0
for _ in range(200):
    pick = rng5.choice(AGR, size=16, replace=True)
    scs = 1.0 + rng5.normal(0, 0.9/18.0)
    docc = []
    for g in pick:
        s_g2 = 1.0 + rng5.normal(0, rel_sig[int(g)])
        s_g2 = min(max(s_g2, 0.5), 1.5)
        if uma_f[int(g)]: s_g2 *= scs
        docc.append(-math.log10(s_g2))
    sr = build_sub(W78, [int(g) for g in pick], dlg_per_occ=docc)
    fb2 = {m: fit_deep(sr, FAMS[m], th0=list(fG[m].x)).fun
           for m in ('BE', 'boot')}
    wins += int(fb2['BE'] - fb2['boot'] > 0)
p_boot = wins/200.0
gb('5 reversal re-grades',
   0.10 <= p_sf <= 0.35 and 1.0 <= zj <= 2.0 and p_boot < 0.97,
   f"d {dG:+.1f}; sign-flip p={p_sf:.4f} (rev 0.2163); LOO z={zj:.2f} "
   f"(rev 1.46); P_boot(200,seed202) = {p_boot:.3f} (rev 0.930 "
   f"seed4801); co-requirement 0.95 {'FAILED' if p_boot < 0.95 else 'passed'}"
   f"  [{time.time()-t0:.0f}s]")

# ---- GB-6: confound census ----
Q("")
rc_s, _p = spearmanr([FGAS[g] for g in LEGA], LYm)
rc_f, _p2 = spearmanr([FGAS[g] for g in FLOW],
                      [float(np.median(ly[W78['gpts'][g]]))
                       for g in FLOW])
WMID = 0.5*(WLO + WHI)
ac_a = np.mean([(FGAS[g] < 0.5) ==
                (float(np.median(ly[W78['gpts'][g]])) > WMID)
                for g in LEGA])
dropped = [g for g in LEGA if g not in set(LEGW)]
dr_rich = sum(1 for g in dropped if FGAS[g] >= 0.5)
dr_poor = len(dropped) - dr_rich
gb('6 confound census',
   rc_s <= -0.8 and rc_f <= -0.7 and abs(ac_a - 0.59) <= 0.03
   and dr_rich == 8 and dr_poor == 6,
   f"spearman anch {rc_s:.2f} (rev -0.89) flow {rc_f:.2f} (rev "
   f"-0.81); a_conf(anch) {ac_a:.2f} (rev 0.59); window drops "
   f"{dr_rich}/16 gasrich vs {dr_poor}/62 gaspoor (rev 8, 6)")

# ---- GB-7: flow_gasrich census + carrier + ridge ----
Q("")
Q("-- GB-7 flow_gasrich: coverage, carrier, ridge --")
FGR = [g for g in FLOW if FGAS[g] >= 0.5]
fgr_pts = np.concatenate([W78['gpts'][g] for g in FGR])
y_own = gN1_all[fgr_pts]/6.2684e-11
npts = len(fgr_pts)
fgas_1003 = FGAS[GI['NGC1003']]
subFG = build_sub(W78, FGR)
thg = {'BE': [math.log10(6.2684e-11), 1.9813, 0.0348, 0.0],
       'boot': [math.log10(6.2804e-11), 2.3487, 0.0324, 0.0]}
fFG = {m: fit_deep(subFG, FAMS[m], th0=thg[m])
       for m in ('BE', 'boot')}
dFG = fFG['boot'].fun - fFG['BE'].fun
gl = [g for g in FGR if NAMEOF[g] != 'NGC1003']
s2 = build_sub(W78, gl)
f2 = {m: fit_deep(s2, FAMS[m], th0=list(fFG[m].x)).fun
      for m in ('BE', 'boot')}
cost1003 = dFG - (f2['boot'] - f2['BE'])
def prof_at(la0_fix, th0):
    """FULL independent profiled objective (indep_obj incl. priors)
    minimised over (f, s_int, u) at fixed la0 -- the ridge probe."""
    def objx(t3):
        f_, s_, u_ = t3
        if not (0.3 < f_ < 5.0) or not (1e-3 <= s_ < 0.4) \
           or abs(u_) > 0.2:
            return 1e12
        return indep_obj(subFG, 'BE', [la0_fix, f_, s_, u_])
    r = minimize(objx, th0, method='Nelder-Mead',
                 options=dict(maxiter=120, xatol=1e-4, fatol=1e-4))
    return r.fun
pr = {}
for a0v in (6.2684e-11, 8.0e-11, 1.0e-10):
    pr[a0v] = prof_at(math.log10(a0v), [2.0, 0.035, 0.0])
    Q(f"    ridge probe a0 {a0v:.3e}: profiled {pr[a0v]:.2f} "
      f"[{time.time()-t0:.0f}s]")
flatness = pr[6.2684e-11] - pr[8.0e-11]
rise = pr[1.0e-10] - pr[8.0e-11]
gb('7 flow_gasrich',
   float(y_own.max()) < 1.0 and npts == 216
   and abs(fgas_1003 - 0.52) <= 0.01
   and abs(cost1003 - 15.46) <= 2.0
   and abs(flatness) <= 3.0 and 0.0 <= rise <= 12.0,
   f"max y(own a0) {float(y_own.max()):.3f} (rev 0.833), pts {npts} "
   f"(rev 216, 0 above y=1: {int((y_own > 1).sum())}); NGC1003 fgas "
   f"{fgas_1003:.2f} (rev 0.52), drop costs {cost1003:+.1f} (rev "
   f"15.46); ridge surrogate d(6.3-8.0) {flatness:+.1f} ~flat, "
   f"d(10-8) {rise:+.1f} rising (rev +0.9/+3.7 shape) [SHAPE-ONLY]")

# ---- GB-8: pin-3 SPARC-convention arithmetic ----
Q("")
import glob as _glob
path_u = [p for p in _glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                                recursive=True) if 'UGC03580' in p][0]
gbar0, gobs0, gbar1, gobs1 = [], [], [], []
s_rs = 1.25
for l in open(path_u):
    if l.startswith('#'): continue
    t = l.split()
    if len(t) < 6: continue
    R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
    if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
    for (Rx, sc, gb_, go_) in ((R, 1.0, gbar0, gobs0),
                               (R*s_rs, s_rs, gbar1, gobs1)):
        vg = Vg*math.sqrt(sc); vd = Vd*math.sqrt(sc)
        vb = Vb*math.sqrt(sc)
        gb_.append((vg*abs(vg) + 0.5*vd*abs(vd) + 0.7*vb*vb)/Rx)
        go_.append(Vo*Vo/Rx)
dbar = max(abs(math.log10(a/b)) for a, b in zip(gbar1, gbar0))
dobs = [math.log10(a/b) for a, b in zip(gobs1, gobs0)]
gb('8 pin-3 exactness',
   dbar <= 1e-12 and abs(dobs[0] + math.log10(s_rs)) <= 1e-12
   and (max(dobs) - min(dobs)) <= 1e-12,
   f"g_bar invariance {dbar:.2e} (rev 1.9e-16); g_obs shift "
   f"{dobs[0]:+.6f} = -log10(1.25) exactly, spread "
   f"{max(dobs)-min(dobs):.2e}")

# ---- GB-9: 2x unmodelled-scatter erosion ----
Q("")
Q(f"-- GB-9 2x scatter erosion (2 seeds, BE/boot)  "
  f"[{time.time()-t0:.0f}s] --")
SIG2 = 2.0*SIG_INJ
ds9 = []
for seed in (7001, 7002):
    rng9 = np.random.default_rng(seed)
    # first draw the 1x stream shape then scale: use fresh draws at 2x
    dlg9 = [float(rng9.normal(0, SIG2)) for _ in LEGA]
    s9 = build_sub(W78, LEGA, dlg_per_occ=dlg9)
    f9 = {m: fit_deep(s9, FAMS[m], th0=list(bA[m].x)).fun
          for m in ('BE', 'boot')}
    ds9.append(f9['boot'] - f9['BE'])
m9 = float(np.mean(ds9))
gb('9 erosion at 2x', 15.0 <= m9 <= 55.0,
   f"mean d(BE better) = {m9:+.1f} at 2x0.1086 dex (rev ~32.0, "
   f"49% erosion; seeds differ from reviewer's -> RANGE)")

# ---- GB-10: frozen-bulge probe ----
Q("")
subFb = build_sub(W78, FLOW)
subFb = dict(subFb)
subFb['gd'] = subFb['gd'] + subFb['gb']
subFb['gb'] = np.zeros_like(subFb['gb'])
fb10 = {m: fit_deep(subFb, FAMS[m],
                    th0=[math.log10(ARC[m]), 1.0, 0.08, 0.0]).fun
        for m in ('BE', 'boot')}
d10 = fb10['BE'] - fb10['boot']
gb('10 bulge-in-free-slot', 30.0 <= d10 <= 46.0,
   f"flow lead with bulge in the free-M/L slot = {d10:+.1f} "
   f"(rev 37.93 from 48.73, -22%)")

# ---- GB-11: sv_swap UMa audit + gate coverage ----
Q("")
meta10 = NS10U['meta']
recov, true_i = [], []
for g in LEGA:
    if W78['fd'][g] != 4: continue
    sv, rel = W78['sigv'][g], W78['rel'][g]
    recov.append(math.sqrt(max(sv*sv - (rel/LN10)**2, 0.0)))
    inc, q, D0, eD0, einc, fdv = meta10[NAMEOF[g]]
    irad = math.radians(inc)
    true_i.append(2.0*(math.radians(max(einc, 1.0))
                       / math.tan(irad))/LN10)
med_r, med_t = float(np.median(recov)), float(np.median(true_i))
sw_ok = abs(med_r - 0.0252) <= 0.004 and abs(med_t - 0.0087) <= 0.004
uma_in20 = sum(1 for g in LEGA[:20] if W78['fd'][g] == 4)
gb('11 sv_swap audit', sw_ok and uma_in20 == 0,
   f"UMa recovered-incl median {med_r:.4f} (rev 0.0252) vs true "
   f"{med_t:.4f} (rev 0.0087); UMa in LEGA[:20] = {uma_in20} (rev 0)")

# ---- GB-12: UMa-26 u=0 refits ----
Q("")
fit_hier = NSV['fit_hier']
UMA = [g for g in LEGA if W78['fd'][g] == 4]
subU = build_sub(W78, UMA)
r12 = {}
for m, (wf, wu0, wa0) in (('BE', (-1417.505, -1401.457, 1.0038e-10)),
                          ('boot', (-1429.708, -1408.837,
                                    1.0675e-10))):
    bu = fit_deep(subU, FAMS[m],
                  th0=[math.log10(ARC[m]), 1.3, 0.03, -0.08])
    b0 = fit_hier(subU, FAMS[m], use_u=False,
                  th0=[math.log10(ARC[m]), 1.3, 0.03],
                  tol=5e-4, max_rounds=60)
    r12[m] = (bu.fun, b0.fun, 10**b0.x[0], bu.fun - b0.fun)
gb('12 UMa u-degeneracy',
   abs(r12['BE'][3] + 16.05) <= 2.0
   and abs(r12['boot'][3] + 20.87) <= 2.0
   and abs(r12['BE'][2] - 1.0038e-10)/1.0038e-10 <= 0.02,
   f"freeing u buys BE {-r12['BE'][3]:.1f} (rev 16.05) boot "
   f"{-r12['boot'][3]:.1f} (rev 20.87); u=0 a0(BE) = "
   f"{r12['BE'][2]:.4e} (rev 1.0038e-10)")

# ---- GB-13: P15/P85 window variants ----
Q("")
Q(f"-- GB-13 P15/P85 windows  [{time.time()-t0:.0f}s] --")
WLO2 = max(np.percentile(ly[apts], 15), np.percentile(ly[fpts], 15))
WHI2 = min(np.percentile(ly[apts], 85), np.percentile(ly[fpts], 85))
WW2, LEGW2, FLOWW2 = window_world(W78, (ly >= WLO2) & (ly <= WHI2))
res13 = {}
for tag, gl2 in (('anch', LEGW2), ('flow', FLOWW2)):
    s13 = build_sub(WW2, gl2)
    f13 = {}
    for m in MEMBERS:
        f13[m] = fit_deep(s13, FAMS[m],
                          th0=[math.log10(ARC[m]), 1.3, 0.04, 0.0])
    fmin13 = min(f13[m].fun for m in MEMBERS)
    win13 = min(MEMBERS, key=lambda m: f13[m].fun)
    d13 = f13['boot'].fun - f13['BE'].fun
    fldw = sub_fields(s13, 'boot' if d13 < 0 else 'BE',
                      list(f13['boot' if d13 < 0 else 'BE'].x))
    muw = sub_fields(s13, 'BE' if d13 < 0 else 'boot',
                     list(f13['BE' if d13 < 0 else 'boot'].x))[0]
    sd13 = sd_block_f(s13, fldw, muw)
    res13[tag] = (win13, d13, sd13, len(gl2))
gb('13 P15/P85',
   res13['flow'][0] == 'boot' and abs(res13['flow'][1] + 30.8) <= 3.0
   and res13['flow'][2] <= 22.0
   and res13['anch'][0] in ('p065', 'boot', 'gm')
   and abs(res13['anch'][1] + 8.6) <= 3.0,
   f"flow: winner {res13['flow'][0]} d {res13['flow'][1]:+.1f} "
   f"(rev -30.8) SD {res13['flow'][2]:.1f} (rev 16.0 -> ratio "
   f"{abs(res13['flow'][1])/res13['flow'][2]:.2f} vs rev 1.93); "
   f"anch: winner {res13['anch'][0]} (rev p065) d "
   f"{res13['anch'][1]:+.1f} (rev -8.6)")

# ---- GB-14: fgas threshold spots ----
Q("")
r14 = {}
for thr, leg, gals, want_d in ((0.4, 'anch', LEGA, -10.5),
                               (0.6, 'flow', FLOW, -5.1)):
    gl3 = [g for g in gals if FGAS[g] >= thr]
    s14 = build_sub(W78, gl3)
    f14 = {m: fit_deep(s14, FAMS[m],
                       th0=[math.log10(ARC[m]), 1.5, 0.04, 0.0])
           for m in ('BE', 'boot')}
    d14 = f14['boot'].fun - f14['BE'].fun
    r14[(leg, thr)] = (len(gl3), d14,
                       f14['BE'].x[1], f14['boot'].x[1])
gb('14 threshold spots',
   r14[('anch', 0.4)][0] == 21
   and abs(r14[('anch', 0.4)][1] - (-10.5)) <= 2.5
   and r14[('flow', 0.6)][0] == 15
   and r14[('flow', 0.6)][2] >= 2.45
   and r14[('flow', 0.6)][3] >= 2.45,
   f"anch@0.4: n {r14[('anch',0.4)][0]} (21) d "
   f"{r14[('anch',0.4)][1]:+.1f} (rev -10.5); flow@0.6: n "
   f"{r14[('flow',0.6)][0]} (15) f_ML {r14[('flow',0.6)][2]:.3f}/"
   f"{r14[('flow',0.6)][3]:.3f} (rev pinned 2.500)")

# ---- GB-15: arithmetic pack ----
Q("")
Q("-- GB-15 arithmetic pack --")
from scipy.stats import norm
packs = [
    ('flow z 48.7/68.24 = 0.71', abs(48.7/68.24 - 0.71) <= 0.005),
    ('anch z 63.3/22.96 = 2.76', abs(63.3/22.96 - 2.757) <= 0.01),
    ('p 0.2163 two-sided -> 1.24 sigma',
     abs(norm.ppf(1 - 0.2163/2) - 1.237) <= 0.01),
    ('p 0.0039 -> 2.66 sigma; rev quotes 2.66',
     abs(norm.ppf(1 - 0.0039) - 2.66) <= 0.02),
    ('7.77/5.33 = 1.46', abs(7.77/5.33 - 1.458) <= 0.01),
    ('multiplicity 1-(1-p)^N: p=.08 N=4 -> 0.28',
     abs(1 - (1 - 0.08)**4 - 0.284) <= 0.005),
    ('p=.16 N=8 -> 0.75', abs(1 - (1 - 0.16)**8 - 0.752) <= 0.005),
    ('P15/85 flow ratio 30.8/16.0 = 1.93',
     abs(30.8/16.0 - 1.925) <= 0.01),
    ('erosion 49.8/63.3 -> 21%', abs(1 - 49.8/63.3 - 0.213) <= 0.01),
    ('bulge -22%: 37.93/48.73', abs(1 - 37.93/48.73 - 0.2216) <= 0.005),
    ('Ups flow 1.426*0.5/1.747*0.5 = 0.71/0.87',
     abs(1.426*0.5 - 0.713) <= 0.002 and abs(1.747*0.5 - 0.8735)
     <= 0.002),
    ('sv 2.4%: 0.1115/0.1089', abs(0.1115/0.1089 - 1.0239) <= 0.003),
]
p15 = True
for txt, ok in packs:
    p15 &= ok
    Q(f"  {'OK ' if ok else 'BAD'} {txt}")
gb('15 arithmetic', p15, "12-item pack")

# ---- GB-17: innermost-point angular census (disclosure-grade) ----
Q("")
SPd = NS10U['SP']
ROTPATH = {os.path.basename(p).replace('_rotmod.dat', ''): p
           for p in _glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                               recursive=True)}
def theta_min(gals):
    th = []
    for g in gals:
        nm2 = NAMEOF[int(g)]
        pth = [ROTPATH[nm2]] if nm2 in ROTPATH else []
        if not pth: continue
        D0 = SPd['dist'][int(g)][0]
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
thF, thA2 = theta_min(FLOW), theta_min(LEGA)
gb('17 inner-arcsec census',
   10.0 <= thF <= 20.0 and 25.0 <= thA2 <= 45.0,
   f"innermost kept point median: flow {thF:.1f} arcsec (rev 15.0), "
   f"anchored {thA2:.1f} (rev 34.8; catalog-D convention, "
   f"disclosure-grade RANGE)")

# ---- GB-18: family-curve separation profile ----
Q("")
ys = np.array([0.01, 1.0, 2.0, 30.0])
sep18 = np.log10(FAMS['BE'](ys)/FAMS['boot'](ys))
gb('18 curve separation',
   abs(sep18[0] - 0.011) <= 0.004
   and 0.06 <= max(sep18[1], sep18[2]) <= 0.08
   and abs(sep18[3]) <= 0.005,
   f"log10(nuBE/nuboot) at y=0.01/1/2/30 = "
   + "/".join(f"{v:+.3f}" for v in sep18)
   + " (rev +0.011, peak +0.069-0.072, +0.002)")

# ---- GB-19: carrier metadata + bulge-dominated counts ----
Q("")
inc_u, q_u, D_u, eD_u, einc_u, fdv_u = meta10['UGC03580']
npt_u = len(W78['gpts'][GI['UGC03580']])
pts_u = W78['gpts'][GI['UGC03580']]
bd_u = int(np.sum(W78['gb'][pts_u] >
                  W78['gg'][pts_u] + W78['gd'][pts_u]))
bd_f = int(np.sum(W78['gb'][fpts] > W78['gg'][fpts] + W78['gd'][fpts]))
gb('19 carrier metadata',
   abs(D_u - 20.7) <= 0.2 and abs(inc_u - 63) <= 1 and q_u == 2
   and npt_u == 39 and abs(bd_u - 13) <= 3 and abs(bd_f - 322) <= 40,
   f"UGC03580: D {D_u} eD {eD_u} inc {inc_u} Q {q_u} pts {npt_u} "
   f"(rev 20.7/5.2/63/2/39); bulge-dominated {bd_u}/39 (rev 13), "
   f"flow {bd_f}/1437 (rev 322) [dominance def approximate: RANGE]")

Q("")
Q(f"GB POST-REPORT HALF: "
  f"{'ALL REVIEWER NUMBERS CONFIRMED' if ok_all else 'DIFFERENCES PRESENT (see DIFFER rows)'}")
Q(f"wall-clock: {(time.time()-t0)/60:.1f} min")
save()
