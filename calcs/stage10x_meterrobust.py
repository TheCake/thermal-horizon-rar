"""
STAGE 10X -- THE METER ROBUSTNESS BATTERY (pre-reg PREREG-METER-10X.md,
committed 8f11fe7; amendments A1/A2 5ae0d13, both before any run).

Referee-grade robustness axes on the anchored central, never run before:
A1 leave-one-out influence, A2 inclination ladder (40/45/50), A3 Q=1
only, A4 bulge-free, A5 nuisance-prior width x1.5 -- each against a
null-scatter-scaled bar max(2 sigma_exp, 1.0), sigma_exp = sigma_stat *
sqrt(N_full/N_sub - 1). DIAGNOSTIC grade: no variant replaces the
primary; band 65.4-70.4 and stat 5.0 untouched by construction; every
credence cell pre-signed HOLD. Census leg (no letter clause): D564-8 /
D631-7 CF4 cross-check + KK98-251 resolution via session-documented
positions (data/stage10x_pgc.json).

Estimator inheritance: exec of calcs/stage10t_legregrow.py to the
G10T-1 marker (data load, worlds, nu families, fit_hier/fit_deep).
Modes:
  py calcs/stage10x_meterrobust.py gates -> data/stage10x_gates.txt
     (reproductions of ARCHIVED numbers + axis counts; NO variant
     central computed)
  py calcs/stage10x_meterrobust.py run   -> data/stage10x_battery.txt
     (re-verifies gates in-process, then fires the battery)
"""
import json, math, os, re, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
MODE = sys.argv[1] if len(sys.argv) > 1 else 'gates'
RUN = (MODE == 'run')
OUTFILE = 'data/stage10x_battery.txt' if RUN else 'data/stage10x_gates.txt'

LX = []
def PX(s=""):
    print(s, flush=True)
    LX.append(s)

def savex():
    with open(OUTFILE, 'w', encoding='utf-8') as f:
        f.write("\n".join(LX) + "\n")

t00 = time.time()
PX(f"STAGE 10X METER ROBUSTNESS BATTERY -- mode = {MODE}")
PX("")

# ---------------- inherit the 10T engine (verbatim, to the marker) ----
src = open('calcs/stage10t_legregrow.py', encoding='utf-8').read()
MARKER = "# ---------------- G10T-1"
assert MARKER in src, "inheritance marker missing"
header = src[:src.index(MARKER)]
ns = {'__name__': 'stage10t_inherited',
      '__file__': os.path.abspath('calcs/stage10t_legregrow.py')}
_argv = sys.argv
sys.argv = ['stage10t_legregrow.py', 'gates']   # inherited MODE constants
exec(compile(header, 'stage10t_legregrow.py[header]', 'exec'), ns)
sys.argv = _argv
PX("[inherited 10T header executed: SPARC/CF4/LT load, worlds, engine]")

make_world = ns['make_world']; build_sub = ns['build_sub']
fit_hier = ns['fit_hier']; fit_deep = ns['fit_deep']
nu_be = ns['nu_be']; h0_of_a0 = ns['h0_of_a0']
meta = ns['meta']; LN10 = ns['LN10']

SIG_STAT = 4.99          # archived primary-engine bootstrap (10T, R46-C1)
N_FULL = 78

# archived reference values (parsed, not retyped)
t10t = open('data/stage10t_skyread.txt', encoding='utf-8').read()
m = re.search(r"PRIMARY hier BE \(DEEP, R46-C2\): a0 = ([\d.e+-]+), .*?"
              r"-> H0_A = ([\d.]+)", t10t)
A0_REF, H0_REF = float(m.group(1)), float(m.group(2))
m = re.search(r"A-var2 .*?: a0 = ([\d.e+-]+) \(H0 ([\d.]+)\)", t10t)
A0_V2REF = float(m.group(1))

# ---------------- primary world + reproduction chain -----------------
W = make_world(lt_slot=None)
legA = [int(g) for g in W['legA']]
assert len(legA) == N_FULL
subT = build_sub(W, legA)

PX("")
PX("-- G10X-1 primary reproduction (archived chain: plain -> loose -> deep) --")
b_plain = fit_hier(subT, nu_be, use_u=False)
b_loose = fit_hier(subT, nu_be, use_u=True, th0=list(b_plain.x) + [0.0])
bT = fit_deep(subT, nu_be, th0=list(b_loose.x))
a0_T = 10**bT.x[0]
H0_T = h0_of_a0(a0_T)
d1 = abs(a0_T / A0_REF - 1)
g1 = d1 <= 1e-3 and abs(H0_T - H0_REF) <= 0.02
PX(f"  a0 = {a0_T:.4e} vs archived {A0_REF:.4e} (rel {d1:.2e} vs 1e-3); "
   f"H0 = {H0_T:.2f} vs {H0_REF:.2f}")
PX(f"G10X-1: {'PASS' if g1 else 'STOP'}  [{time.time()-t00:.0f}s]")
if not g1:
    savex(); sys.exit("G10X-1 STOP")

PX("")
PX("-- G10X-2 variant-harness identity (A-var2 world) --")
Wv2 = make_world(uma_moved=True, lt_slot=None)
bv2 = fit_deep(build_sub(Wv2, [int(g) for g in Wv2['legA']]), nu_be,
               th0=list(bT.x))
d2 = abs(10**bv2.x[0] / A0_V2REF - 1)
g2 = d2 <= 1e-3
PX(f"  a0 = {10**bv2.x[0]:.4e} vs archived {A0_V2REF:.4e} (rel {d2:.2e} "
   f"vs 1e-3)")
PX(f"G10X-2: {'PASS' if g2 else 'STOP'}")
if not g2:
    savex(); sys.exit("G10X-2 STOP")

# ---------------- axis construction (counts; deterministic) ----------
def inc_of(g):
    return meta[W['names'][g]][0]

def q_of(g):
    return meta[W['names'][g]][1]

def bulge_free(g):
    return bool(np.all(W['gb'][W['gpts'][g]] == 0.0))

AXES = {}
for cut in (40, 45, 50):
    AXES[f"A2 incl>={cut}"] = [g for g in legA if inc_of(g) >= cut]
AXES["A3 Q=1"] = [g for g in legA if q_of(g) == 1]
AXES["A4 bulge-free"] = [g for g in legA if bulge_free(g)]
AXES["A5 svx1.5"] = list(legA)

def bar_of(nsub):
    if nsub >= N_FULL:
        return 1.0
    return max(2.0 * SIG_STAT * math.sqrt(N_FULL / nsub - 1.0), 1.0)

PX("")
PX("-- G10X-3 axis census (bars from counts alone) --")
ELIG, THIN = [], []
BAR_LOO = max(2.0 * SIG_STAT * math.sqrt(N_FULL / (N_FULL - 1) - 1.0), 1.0)
PX(f"  A1 LOO: N_sub = {N_FULL-1} per fit, 78 fits; bar = "
   f"{BAR_LOO:.2f} km/s/Mpc (on the max)")
ELIG.append("A1 LOO")
for nm, members in AXES.items():
    nsub = len(members)
    thin = nsub < 20
    (THIN if thin else ELIG).append(nm)
    PX(f"  {nm}: N_sub = {nsub}; bar = {bar_of(nsub):.2f} km/s/Mpc"
       + ("  [N-THIN: co-read only, excluded from letter]" if thin else ""))
    assert set(members) <= set(legA)
PX(f"  letter-eligible set E = {ELIG}")
PX("G10X-3: PASS")

PX("")
PX("-- G10X-4 warm-start determinism (A1 amendment bar: rel 1e-4) --")
b_id = fit_deep(build_sub(W, legA), nu_be, th0=list(bT.x))
d4 = abs(10**b_id.x[0] / a0_T - 1)
g4 = d4 <= 1e-4
PX(f"  identity refit rel delta = {d4:.2e} vs 1e-4")
PX(f"G10X-4: {'PASS' if g4 else 'STOP'}")
if not g4:
    savex(); sys.exit("G10X-4 STOP")

if not RUN:
    PX("")
    PX("GATES: ALL PASS (no variant central computed in gates mode)")
    PX(f"total wall-clock: {(time.time()-t00)/60:.1f} min")
    savex()
    sys.exit(0)

# ===================== THE BATTERY (run mode) ========================
PX("")
PX("== 10X BATTERY ==")
fired = []      # (axis, delta, bar)
rows = []

# ---- A1 leave-one-out ----
PX("")
PX(f"-- A1 leave-one-out (78 deep refits; bar {BAR_LOO:.2f} on the max) --")
t1 = time.time()
infl = []
for g in legA:
    sub_g = build_sub(W, [x for x in legA if x != g])
    bg = fit_deep(sub_g, nu_be, th0=list(bT.x))
    infl.append((W['names'][g], h0_of_a0(10**bg.x[0]) - H0_T))
infl.sort(key=lambda t: -abs(t[1]))
mx_nm, mx_d = infl[0]
PX(f"  max |dH0| = {abs(mx_d):.3f} ({mx_nm}); top influences: "
   + ", ".join(f"{nm} {d:+.3f}" for nm, d in infl[:6]))
PX(f"  distribution |dH0|: median {np.median([abs(d) for _, d in infl]):.3f}, "
   f"90th pct {np.percentile([abs(d) for _, d in infl], 90):.3f}")
PX(f"  [{(time.time()-t1)/60:.1f} min]")
rows.append(("A1 LOO max", abs(mx_d), BAR_LOO))
if abs(mx_d) > BAR_LOO:
    fired.append(("A1 LOO", abs(mx_d), BAR_LOO))

# ---- A2/A3/A4 subset refits ----
for nm in list(AXES.keys()):
    if nm == "A5 svx1.5":
        continue
    members = AXES[nm]
    nsub = len(members)
    bar = bar_of(nsub)
    if nsub < 20:
        PX(f"-- {nm}: N-THIN (N = {nsub}) -- co-read only")
    b = fit_deep(build_sub(W, members), nu_be, th0=list(bT.x))
    h0 = h0_of_a0(10**b.x[0])
    d = h0 - H0_T
    tag = "N-THIN co-read" if nsub < 20 else \
          ("FIRES" if abs(d) > bar else "pass")
    PX(f"-- {nm}: N = {nsub}, H0 = {h0:.2f}, dH0 = {d:+.2f} "
       f"(bar {bar:.2f}) -> {tag}")
    rows.append((nm, abs(d), bar))
    if nsub >= 20 and abs(d) > bar:
        fired.append((nm, abs(d), bar))

# ---- A5 nuisance-prior width x1.5 ----
sub5 = build_sub(W, legA)
sub5['sv'] = sub5['sv'] * 1.5
b5 = fit_deep(sub5, nu_be, th0=list(bT.x))
h05 = h0_of_a0(10**b5.x[0])
d5 = h05 - H0_T
bar5 = 1.0
tag5 = "FIRES" if abs(d5) > bar5 else "pass"
PX(f"-- A5 svx1.5: H0 = {h05:.2f}, dH0 = {d5:+.2f} (bar {bar5:.2f}) "
   f"-> {tag5}")
rows.append(("A5 svx1.5", abs(d5), bar5))
if abs(d5) > bar5:
    fired.append(("A5 svx1.5", abs(d5), bar5))

# ---- the letter (pre-registered grammar; code-cap in the same block) ----
PX("")
if fired:
    PX("LETTER: X-METER-SENSITIVE(" +
       "; ".join(f"{nm} |d| {d:.2f} > bar {b:.2f}" for nm, d, b in fired)
       + ")")
    PX("  [DIAGNOSTIC GRADE -- a SENSITIVE letter triggers a review round "
       "before anything is quoted from it]")
else:
    PX("LETTER: X-METER-ROBUST -- no letter-eligible axis exceeds its bar")
    PX("  [DIAGNOSTIC GRADE -- quotable as robustness co-reads only]")
PX("  No variant reading replaces the primary; the band 65.4-70.4 and "
   "stat 5.0 are uncut-primary properties and do not move (pre-registered "
   "cap, this block).")

# ---- census leg (informational; NO letter clause) ----
PX("")
PX("-- CENSUS LEG (clearing rule frozen in prereg section 6) --")
pos = json.load(open('data/stage10x_pgc.json'))
# CF4 coordinates re-parse (RAdeg 138-145, DEdeg 147-154, 1-indexed)
cf4pos = {}
for l in open('data/cf4/table2.dat'):
    try:
        pgc = int(l[0:7])
    except ValueError:
        continue
    s1, s2 = l[137:145].strip(), l[146:154].strip()
    if s1 and s2:
        cf4pos[pgc] = (float(s1), float(s2))
PX(f"  CF4 coordinate rows parsed: {len(cf4pos)}")

def match_cf4(ra, dec, rad_deg=0.05):
    best, bd = None, rad_deg
    cd = math.cos(math.radians(dec))
    for pgc, (r, d) in cf4pos.items():
        dd = math.hypot((r - ra) * cd, d - dec)
        if dd < bd:
            best, bd = pgc, dd
    return best, bd

adopt = ns['adopt']
for nm in ('D564-8', 'D631-7'):
    ra, dec = pos[nm]['ra'], pos[nm]['dec']
    pgc, dist_deg = match_cf4(ra, dec)
    gi = ns['SPNAME2GI'][nm]
    D_sp, eD_sp = ns['SP']['dist'][gi]
    if pgc is None:
        PX(f"  {nm}: NO CF4 row within 3 arcmin -> flag STANDS "
           f"(no per-method anchor to check against)")
        continue
    a = adopt(pgc)
    if a is None:
        PX(f"  {nm}: CF4 PGC {pgc} ({dist_deg*60:.2f} arcmin) has NO "
           f"anchor-grade per-method modulus -> flag STANDS")
        continue
    mth, D_cf, eD_cf = a
    gap = D_cf / D_sp - 1
    tol = max(2 * eD_sp / D_sp, 0.15)
    ok = abs(gap) <= tol
    PX(f"  {nm}: CF4 PGC {pgc} ({dist_deg*60:.2f} arcmin), method {mth}: "
       f"D_CF4 = {D_cf:.2f} +- {eD_cf:.2f} vs SPARC D = {D_sp:.2f} +- "
       f"{eD_sp:.2f}; gap {gap:+.1%} vs tol {tol:.0%} -> flag "
       f"{'CLEARED' if ok else 'STANDS (gap reported)'}")
    PX(f"    [no membership or distance change either way; membership "
       f"moves only through a registered census stage]")
ra, dec = pos['KK98-251']['ra'], pos['KK98-251']['dec']
pgc, dist_deg = match_cf4(ra, dec)
if pgc is None:
    PX(f"  KK98-251: no CF4 row within 3 arcmin of the NED position; NED "
       f"itself gives only an anonymous field designation (PGC1 0065001 "
       f"NED007, the NGC 6946 field) -> RESOLVED AS UNLISTABLE: no proper "
       f"PGC, no CF4 per-method entry; it cannot be CF4-reclassified and "
       f"stays a flow member at its SPARC distance. Census question CLOSED.")
else:
    PX(f"  KK98-251: nearest CF4 row PGC {pgc} at {dist_deg*60:.2f} arcmin "
       f"-> candidate identification; report for the next census stage "
       f"(no change here).")

PX("")
PX("CREDENCES: every cell pre-signed HOLD (53 / 8); diagnostic stage.")
PX(f"total wall-clock: {(time.time()-t00)/60:.1f} min")
savex()
