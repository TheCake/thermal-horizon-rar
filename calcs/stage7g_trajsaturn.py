"""
STAGE 7G: the trajectory formulation vs Saturn -- making the 6Y corollary
QUANTITATIVE.

The 4K/5S/5I record: every MG (field-formulation) member of the measured
function family sources a solar quadrupole 4.0-5.8x the Cassini bound via
the Sun's transition shell at r_M = sqrt(GM/a0) ~ 7000 AU. The 6Y
derivation reached the one escape BY structure: coupling to the system's
collective/barycentric mode = a trajectory functional (the
modified-inertia class). 7G asks the two quantitative questions that
convert "MI might do it" into numbers:

  Q1 (binaries): does the TRAJECTORY formulation of the AMB function fit
     the wide binaries as well as the field formulation? Machinery = the
     4L mi_t prescription verbatim (per-orbit boost at the Kepler
     time-average <1/r^2>, EFE-respecting via the same table, exact
     Newton engine with G_eff) on the corrected-velocity v7 pipeline,
     with the 6G MG-AMB rows as same-seed comparators.
  Q2 (Saturn): what does the SAME function predict on Saturn's worldline
     in the trajectory formulation? The anomaly attaches to the
     trajectory's own occupation: y_Sat = g_Sun(r_Sat)/a0 ~ 5.4e5, so
     the occupation argument u = y^((1+beta)/2) nu^beta is ~1e3 and
     nu - 1 ~ exp(-u): the transition shell that sources the 4K
     quadrupole NEVER EXISTS on Saturn's trajectory. The residual
     quadrupole is the ordinary galactic tide (~7.5e-31 s^-2, the 6W
     moot block). CONTRAST (order-of-magnitude, labeled): a POWER-LAW
     tail (simple-nu) in the trajectory formulation leaves an
     r-dependent effective G at the (nu-1) ~ a0*r^2/GM ~ 2e-6 level at
     Saturn -- ephemeris-dead by orders; the trajectory door is open
     ONLY for Boltzmann-screened functions, i.e. exactly the screening
     class the data selected.

Gates: G0 the alpha=0 (Newton) rows of the MI runs must equal the 6G
MG runs' Newton lnL per seed (< 0.5; same data, same seed, same
nuisance machinery -- the 4L bit-identity property); table = the
existing gated efe_boost_amb_g1p2.npy (6G).

PRE-REGISTERED BARS:
  B1 binary tie: mean D(lnL, mi_t-AMB minus MG-AMB) over seeds 31/101
     within +-5/seed with alpha-hat interior both seeds -> the
     trajectory formulation fits AS WELL. Worse than -10/seed -> the
     door NARROWS (reported as such). Better than +10 -> flag.
  B2 Saturn margin: the trajectory-AMB solar observable must sit >= 2
     orders under Cassini (expected ~450 orders); the galactic-tide
     residual is quoted alongside.
  VERDICT "SATURN RESOLVED WITHIN THE DERIVED CLASS (formulation
  grade)" iff B1 tie AND B2 -- with the carried costs stated: mi_t is
  a PROXY (per-orbit functional, not a full nonlocal MI theory);
  time-nonlocality/conservation caveats; Petersen & Lelli 2020
  rotation-curve constraints on MI classes; and the MG-formulation
  tension STANDS in the world table (this opens a door, it does not
  delete the row).
Writes data/stage7g_trajsaturn.txt (+ data/stage7g_summary.txt runs).
"""
import math, os, re, sys, time, gc
import numpy as np

L = []
def say(s=''):
    L.append(s); print(s, flush=True)
def save():
    with open('data/stage7g_trajsaturn.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

say("STAGE 7G: the trajectory formulation vs Saturn")
say("=" * 72)

# ---------------- Q2: the Saturn worldline (analytic) -----------------------
GMSUN = 1.32712440018e20
A0 = 1.2e-10
R_SAT = 9.5826*1.495978707e11
G_SAT = GMSUN/R_SAT**2
Y_SAT = G_SAT/A0
G_AMB = 0.1118                      # the solar galactic gate (e_N = 1.2, 6G)
# AMB at deep-Newtonian y: nu -> 1, q -> 1, beta -> g/2
BETA_T = G_AMB/2.0
LNU = 0.5*(1.0 + BETA_T)*math.log(Y_SAT)          # ln u, nu^beta -> 1
U_SAT = math.exp(LNU)
LOG10_ANOM = -U_SAT/math.log(10.0)                # log10(nu - 1) ~ -u/ln10
# Cassini-equivalent fractional acceleration
Q2_CASSINI = 3e-27
DAA_CASSINI = Q2_CASSINI*R_SAT/G_SAT
MARGIN_ORDERS = math.log10(DAA_CASSINI) - LOG10_ANOM
Q_TIDE = 1.9e-10/2.55e20                          # 6W moot block
R_M = math.sqrt(GMSUN/A0)
say("Q2 the Saturn worldline (trajectory formulation of AMB):")
say(f"  y_Sat = g_Sun(r_Sat)/a0 = {Y_SAT:.3e}")
say(f"  occupation argument u = y^((1+beta)/2) with tail beta = g/2 = "
    f"{BETA_T:.4f}: u = {U_SAT:.1f}")
say(f"  trajectory anomaly nu - 1 ~ 10^({LOG10_ANOM:.1f})")
say(f"  Cassini-equivalent da/a = Q2*r/g = {DAA_CASSINI:.2e} "
    f"(Q2 = 3e-27 s^-2, Hees+14)")
say(f"  MARGIN: the trajectory observable sits {MARGIN_ORDERS:.0f} ORDERS "
    f"OF MAGNITUDE under Cassini")
say(f"  residual quadrupole = the galactic tide ~ {Q_TIDE:.1e} s^-2 "
    f"= {Q2_CASSINI/Q_TIDE:.0f}x under the bound (6W)")
say(f"  transition radius r_M = sqrt(GM/a0) = {R_M/1.495978707e11:.0f} AU:")
say(f"  trajectories crossing ~7 kAU (Oort bodies, extreme-aphelion")
say(f"  comets) DO sample the transition -- the residual solar probe.")
nu1_simple = 1.0/Y_SAT                            # simple-nu tail: 1/y
say(f"  CONTRAST (order-of-magnitude, labeled): a power-law tail leaves")
say(f"  nu-1 ~ 1/y = {nu1_simple:.2e} at Saturn as an r-DEPENDENT G_eff --")
say(f"  ephemeris-dead by ~3-4 orders; the trajectory door is open ONLY")
say(f"  for Boltzmann-screened functions (the measured screening class).")
b2 = MARGIN_ORDERS >= 2.0
say(f"B2: {'PASS' if b2 else 'FAIL'} (>= 2 orders; got {MARGIN_ORDERS:.0f})")
say('')
save()

# ---------------- Q1: mi_t-AMB on the binaries ------------------------------
PRESC, ISO = 'mi_t', False
src = open('calcs/stage3p_v7budget.py', encoding='utf-8-sig').read()

OLD_DATA = """vt_d = (4.74047/plx[ok]*np.hypot(d['pmra1'][ok]-d['pmra2'][ok],
        d['pmdec1'][ok]-d['pmdec2'][ok]))/(0.9417*np.sqrt(Mtot_d/s_d))
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
vx_ = d['pmra2']-d['pmra1']; vy_ = d['pmdec2']-d['pmdec1']"""
NEW_DATA = """dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
sn_ = np.maximum(np.hypot(sx_, sy_), 1e-12)
sux_, suy_ = sx_/sn_, sy_/sn_
_cn = list(d.columns.names)
def _pick(*names):
    for n in names:
        if n in _cn: return np.asarray(d[n], dtype=np.float64)
    raise KeyError(names)
r1_ = _pick('radial_velocity1', 'dr2_radial_velocity1', 'rv1')
r2_ = _pick('radial_velocity2', 'dr2_radial_velocity2', 'rv2')
try:
    er1_ = _pick('radial_velocity_error1', 'dr2_radial_velocity_error1')
    er2_ = _pick('radial_velocity_error2', 'dr2_radial_velocity_error2')
except KeyError:
    er1_ = np.full(len(r1_), 2.0); er2_ = np.full(len(r2_), 2.0)
h1_, h2_ = np.isfinite(r1_), np.isfinite(r2_)
w1_ = np.where(h1_, 1.0/np.maximum(er1_, 0.5)**2, 0.0)
w2_ = np.where(h2_, 1.0/np.maximum(er2_, 0.5)**2, 0.0)
rvs_ = (np.where(h1_, r1_, 0.0)*w1_ + np.where(h2_, r2_, 0.0)*w2_) \\
       / np.maximum(w1_+w2_, 1e-12)
th_ = sep/(2.06265e8/plx)
pmcor_ = rvs_*th_*plx/4.74047
vx_ = d['pmra2']-d['pmra1'] + pmcor_*sux_
vy_ = d['pmdec2']-d['pmdec1'] + pmcor_*suy_
vt_d = (4.74047/plx[ok]*np.hypot(vx_[ok], vy_[ok]))/(0.9417*np.sqrt(Mtot_d/s_d))"""
assert src.count(OLD_DATA) == 1
src = src.replace(OLD_DATA, NEW_DATA)

OLD_TABS = """if GTAG == '1p9':
    PS, PB = 'data/efe_boost_simple.npy', 'data/efe_boost_be.npy'
else:
    PS = f'data/efe_boost_simple_g{GTAG}.npy'
    PB = f'data/efe_boost_be_g{GTAG}.npy'
LNY_U, TAB_S = load_tab(PS)
_,     TAB_B = load_tab(PB)"""
NEW_TABS = """LNY_U, TAB_S = load_tab(LAMPATH)
TAB_B = TAB_S"""
assert src.count(OLD_TABS) == 1
src = src.replace(OLD_TABS, NEW_TABS)

OLD_LOOP = 'for law, TAB in (("simple", TAB_S), ("BE", TAB_B)):'
NEW_LOOP = 'for law, TAB in ((LAMLAW, TAB_S),):'
assert src.count(OLD_LOOP) == 1
src = src.replace(OLD_LOOP, NEW_LOOP)

OLD_SUM = """    P(f"  seed {seed}: BE-minus-simple best lnL = "
      f"{best_lnl['BE']-best_lnl['simple']:+.1f}  "
      f"({(time.time()-t0)/60:.1f} min)")"""
NEW_SUM = """    P(f"  seed {seed}: best lnL = {best_lnl[LAMLAW]:+.3f}  "
      f"({(time.time()-t0)/60:.1f} min)")"""
assert src.count(OLD_SUM) == 1
src = src.replace(OLD_SUM, NEW_SUM)

OLD_OUT = "with open('data/stage3u_summary.txt', 'a') as f:"
NEW_OUT = "with open('data/stage7g_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1
src = src.replace(OLD_OUT, NEW_OUT)

# the 4L engine block, verbatim
OLD_ENG = """                    e_s = e_of(p, eta, wr)
                    vp = vp_c(p, e_s, tab_a) if al > 0 else None
                    mode = 5 if al > 0 else 1
                    kw = dict(a0=A0_CAN, tab=tab_a, lny0=LNY0, dlny=DLNY,
                              vp=vp) if al > 0 else {}
                    o = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'],
                            p['uph'], 8, 2500, mode, **kw)"""
NEW_ENG = """                    e_s = e_of(p, eta, wr)
                    gch = GM*p['M_s']/p['a_s']**2
                    if PRESC.startswith('mi_t'):
                        gch = gch/np.sqrt(np.maximum(1.0-e_s**2, 1e-4))
                    if ISO:
                        yv = np.maximum(gch/A0_CAN, 1e-12)
                        if law == "BE":
                            nub = 1.0/(1.0-np.exp(-np.minimum(np.sqrt(yv),
                                                              40.0)))
                        else:
                            nub = 0.5+np.sqrt(0.25+1.0/yv)
                        Bp = 1.0 + al*(nub-1.0)
                    else:
                        Bp = np.interp(np.log(gch/A0_CAN), LNY_U, tab_a,
                                       right=1.0)
                    o = run(p['a_s'], e_s, p['psi0'], p['f_ip'],
                            p['M_s']*Bp, p['uph'], 8, 2500, 1)"""
assert src.count(OLD_ENG) == 1
src = src.replace(OLD_ENG, NEW_ENG)

SEEDS = ['31', '101']
have = os.path.exists('data/stage7g_summary.txt') and \
    open('data/stage7g_summary.txt').read().count(' ambmi: ') >= len(SEEDS)
if not have:
    ns2 = {'__name__': '__main__', 'PRESC': PRESC, 'ISO': ISO,
           'LAMPATH': 'data/efe_boost_amb_g1p2.npy', 'LAMLAW': 'ambmi'}
    sys.argv = ['stage7g', '1p2'] + SEEDS
    print("\n===== mi_t-AMB (trajectory formulation) on binaries =====",
          flush=True)
    exec(compile(src, 'stage3p_patched_7g_ambmi', 'exec'), ns2)
    del ns2; gc.collect()
say("binary runs done")
save()

# ---------------- parse + verdict -------------------------------------------
def parse(path):
    txt = open(path).read()
    out = {}
    for m in re.finditer(
            r'seed (\d+) (\w+): a_hat=([0-9.]+) \(grid ([0-9.]+), '
            r'interior=(\w+)\), dlnL\(Newton\)=\+([0-9.]+), wr=([0-9.]+),'
            r'.*?\n\s*seed \1: best lnL = (-?[0-9.]+)', txt):
        s, law, ah, agrid, inter, dn, wr, lnl = m.groups()
        out[(law, int(s))] = dict(a=float(ah), interior=(inter == 'True'),
                                  lnl=float(lnl), dn=float(dn))
    return out

REC = {}
for p in ('data/stage6g_summary.txt', 'data/stage7g_summary.txt'):
    if os.path.exists(p): REC.update(parse(p))

say('')
say("Q1 binaries (corrected velocities; MG-AMB comparators = 6G same-seed):")
ds, ahs, ints, g0ok = [], [], 0, True
for s in (31, 101):
    if ('ambmi', s) not in REC or ('amb', s) not in REC:
        say(f"  seed {s}: MISSING ROW"); continue
    mi, mg = REC[('ambmi', s)], REC[('amb', s)]
    newton_mi = mi['lnl'] - mi['dn']
    newton_mg = mg['lnl'] - mg['dn']
    g0 = abs(newton_mi - newton_mg) < 0.5
    g0ok &= g0
    d = mi['lnl'] - mg['lnl']
    ds.append(d); ahs.append(mi['a']); ints += mi['interior']
    say(f"  seed {s}: mi_t-AMB lnL {mi['lnl']:+.2f} vs MG-AMB "
        f"{mg['lnl']:+.2f} -> D = {d:+.2f}; a_hat = {mi['a']:.2f} "
        f"(interior={mi['interior']}); dlnL(Newton) = +{mi['dn']:.1f}; "
        f"G0 Newton-row match d = {newton_mi-newton_mg:+.3f} "
        f"{'OK' if g0 else 'FAIL'}")
say(f"G0 (bit-identity of the alpha=0 rows): {'PASS' if g0ok else 'FAIL'}")
ds = np.array(ds)
b1_tie = len(ds) == 2 and np.all(np.abs(ds) <= 5.0) and ints == 2
b1_narrow = len(ds) == 2 and np.mean(ds) <= -10.0
say(f"  mean D = {ds.mean():+.2f}; interior {ints}/2")
say(f"B1: {'TIE (fits as well)' if b1_tie else ('NARROWS' if b1_narrow else 'see rows')}")
say('')
if b1_tie and b2 and g0ok:
    verd = ("SATURN RESOLVED WITHIN THE DERIVED CLASS (formulation grade): "
            "the trajectory formulation of the SAME measured function fits "
            "the binaries equally well and passes Cassini by "
            f"~{MARGIN_ORDERS:.0f} orders of magnitude")
elif b1_narrow:
    verd = "THE TRAJECTORY DOOR NARROWS (binary penalty >= 10/seed)"
else:
    verd = "PARTIAL -- see components"
say(f"VERDICT (pre-registered): {verd}")
say("carried costs: mi_t is a per-orbit PROXY, not a full nonlocal MI")
say("theory (time-nonlocality + conservation need the full formalism);")
say("Petersen & Lelli 2020 constrain specific MI models on rotation")
say("curves; the MG-formulation quadrupole tension STANDS in the world")
say("table -- this stage opens a door by numbers, it does not delete the")
say("row. The r-dependent-G contrast is order-of-magnitude, labeled.")
save()
print("\nSTAGE 7G done -> data/stage7g_trajsaturn.txt")
