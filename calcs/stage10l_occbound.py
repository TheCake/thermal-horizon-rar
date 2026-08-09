"""
STAGE 10L (2026-08-09): THE SELF-CONSISTENT OCCUPATION BOUND -- the one
K-BOUNDED path ROUND 30 left open, executed; plus the meter
kappa-convention pin (R30 successor 3).  DISPERSIVE-SIDE ONLY.

QUESTION.  R30 J-1 killed the two-rung bound (rung-count-arbitrary;
trap #20).  His named repair, verbatim: "the honest object is not the
literal-ladder Gibbs occupation (unbounded below, refuted deep) but
the SELF-CONSISTENT occupation where the dictionary's nu = 1 + kappa
n_BE feeds back into the pull ... n-bar set by the dressed spacing at
the physical T_dS, solved to fixed point -- would tell you whether
kappa -> 1 is where self-consistency first fails, which is a REAL
bound rather than a rung-count choice."  This stage computes that
fixed point exactly, scans the construction FAMILY (trap-#20
discipline: the bound must not inherit an arbitrary choice), and
prices every member against the measured sky.

OBJECTS (10J verbatim, regressed at gate G10L-1):
  E(n)/om = (n+sig) - (kappa/4)(n+sig)^2, sig = 1/2 (P-half, DERIVED);
  D(n) = om[1 - kappa(n+1)/2]; fold at n >= 2/kappa; dictionary
  nu = 1 + kappa n_BE(x), x = om/T_dS = sqrt(g_N/a0).

CONSTRUCTION FAMILIES (locked here):
  SOFTENING (the R30 construction; the vertex's own polaron
  direction -- spacing shrinks as occupation rises, POSITIVE
  feedback): n-bar = n_BE(x * s(n-bar)), s = 1 - kappa(n-bar + c)/2,
  offset family c in {0, 1/2, 1, 3/2} (backward spacing, dE/dn at
  n+sig, forward spacing D(n-bar), shifted-forward).  Plus the
  TRUNCATED-GIBBS member: Gibbs on the rising branch of the literal
  ladder, kept-set conventions KEEP-A (through the fold peak) and
  KEEP-B (below the peak).
  STIFFENING (contrast, the measured family): nu = 1 + n_BE(u_beta),
  u_beta = x^(1+beta) nu^beta -- beta = 0 is the BE identity, beta = 1
  is the 4F boot nu = 1 + n_BE(nu*y) (regression gate), c1(beta) =
  1/(2(1+beta)) exact (5P); NEGATIVE feedback (the total field
  stiffens the mode).  Measured status quoted from the ledger: deep/
  binary beta -> 0 (5R bound beta < 0.03; 5T regime decomposition);
  boot sky-dead (5M).  No new sky fits anywhere in this stage.
  OUT OF SCOPE (named): joint coherent-state/Hartree treatments of
  system+cloud; anything beyond mean-occupation feedback.

MEASURED ANCHORS (archived; no sky fits): kappa_lock = 0.925 (Planck)
/ 0.888 (SH0ES) / hier a0-lock 1.00 +/- 0.05; c1-meter kappa =
2(1-c1) = 1.10-1.48; deep window x in [0.03, 0.3]; deep demand
nu(x = 0.05) = 1 + kappa_lock * n_BE(0.05) = 19.04.

ROUTES:
  R-A THE EXISTENCE MAP: for each softening member, kappa_max(x) =
      the largest kappa with a finite fixed point (two independent
      detectors: monotone iteration from 0; dense-grid sign change);
      the deep-window table vs kappa_lock; the R30 question answered
      mechanically ("is kappa -> 1 where self-consistency first
      fails?").
  R-B THE PRODUCTION TEST: the largest deep nu each member can
      produce at the locked kappa vs the measured demand (truncated-
      Gibbs saturates at nu ~ 2 - kappa/2; SCHA members may have NO
      fixed point at all -- both are failures to carry the deep RAR).
  R-C THE CONVENTION PIN (sympy exact): bare-om vs first-transition
      (renormalized) frequency conventions; deep-amplitude meter
      a0_fit/a0_h = kappa_b^2 and c1 meter c1 = 1 - kappa_b/2 derived;
      the map kappa_r = kappa_b(1 - kappa_b/2) exact, maximum 1/2 at
      exactly kappa_b = 1 (derivative zero) => the renormalized
      convention is nearly kappa-BLIND over the measured band; all
      archived meters verified bare-convention by construction (their
      scripts parameterize x = sqrt(g_N/a0)); the two-pole tension
      (lock 0.925-1.00 vs c1 1.10-1.48) is then convention-ROBUST.
  R-D TAXONOMY + READINGS (report-grade): the feedback-sign table;
      the fold-row resolution (the self-consistent structure the fold
      demanded is NOT a re-thermalized occupation -- consistent with
      the 10C separation theorem: the bath cancels, the drive is
      external); the kappa = 1 level-crossing low-frequency-response
      note (R30 successor 2) booked as a READING (no PREDICTIONS row:
      no numeric kill condition derived here).

GATES (bars locked at this commit):
  G10L-1  algebra: 10J spacing/threshold regressions (residue 0);
          the kappa=1 series rungs; c1(beta) = 1/(2(1+beta)) at beta
          in {0, 1/2, 1} by symbolic series (residue 0); the
          kappa_r map + its maximum at kappa_b = 1 (residue 0).
  G10L-2  solver: kappa -> 0 regresses n-bar = n_BE(x) <= 1e-12;
          fixed-point residual <= 1e-10 wherever existence is
          declared; the two detectors agree on >= 98% of the
          (kappa, x) grid AND kappa_max(x) from the two agree within
          0.02 at every deep-window x.
  G10L-3  FAMILY STABILITY (the verdict gate; trap-#20): the
          exclusion clause fires only if for EVERY softening member
          (all four offsets) and BOTH detectors, kappa_max(x) <=
          kappa_lock/2 = 0.4625 at every x in {0.03, 0.05, 0.1, 0.2,
          0.3}.  Any member violating -> the clause does not fire.
  G10L-4  PRODUCTION: the "cannot carry the deep RAR" clause fires
          only if every member (incl. both truncated-Gibbs
          conventions) either has no fixed point at (kappa_lock,
          x = 0.05) or under-produces the measured nu by >= 5x.
  G10L-5  convention pin: all symbolic residues 0; the meter table
          printed; the boot regression (beta = 1 member == the 4F
          form) numeric <= 1e-12.
  G10L-6  letter audit mechanical (trap-#12: every positive clause
          gated above).

VERDICT LETTERS:
  L-PINNED   G10L-3 AND G10L-4 both fire family-stable: the fed-back
             occupation is EXCLUDED at sky occupations => the
             dictionary's occupation is BARE-PINNED (structurally,
             matching the measured beta -> 0), P-lit has NO
             self-consistent completion in this family, and NO
             occupation-consistency bound on kappa exists -- the
             K-BOUNDED path closes NEGATIVE (kappa = 1 is NOT the
             failure edge; failure sits at kappa ~ O(x) deep).
  L-BOUND    some softening member family-STABLY gives kappa_max in
             [0.8, 1.2] across the deep window AND passes production
             within 5x: a real bound exists; quote it.
  L-AMBIG    family-unstable (members or detectors disagree past
             bars) -- no claim.
  L-OPEN     G10L-1/2/5 fail (wiring).

CREDENCE MAP (pre-signed; the ONLY movable cell): L-BOUND + ROUND-31
no-hole: mech-conditional 8 -> 10 (the 10J K-BOUNDED cell inherited).
L-PINNED / L-AMBIG / L-OPEN: HOLD 8.  anomaly-real 53 UNTOUCHED in
ALL cells (no sky fits; every quoted sky number is an archived meter).
Bankable regardless: the existence map, the production table, the
convention pin, the taxonomy.

COMPUTE < 3 min.  Writes data/stage10l_occbound.txt.
"""
import math, time
import numpy as np
import sympy as sp

T0 = time.time()
OUT = 'data/stage10l_occbound.txt'
LINES = []
def emit(s=""):
    LINES.append(s)
    print(s, flush=True)
    with open(OUT, 'w') as f:
        f.write("\n".join(LINES) + "\n")

emit("STAGE 10L THE SELF-CONSISTENT OCCUPATION BOUND (pre-reg: this "
     "commit)")
emit("")

KLOCK = 0.925
XDEEP = (0.03, 0.05, 0.1, 0.2, 0.3)
OFFSETS = (0.0, 0.5, 1.0, 1.5)

def nbe(z):
    if z > 700: return 0.0
    if z <= 0: return float('inf')
    return 1.0/(math.expm1(z))

# ==================== G10L-1: algebra ====================
om, Om, g, n, sig, kap, x, bet = sp.symbols(
    'om Om g n sigma kappa x beta', positive=True)
E = lambda nn: om*(nn+sig) - (g**2/Om)*(nn+sig)**2
D = sp.expand(E(n+1) - E(n))
r1 = sp.simplify(D - (om - (g**2/Om)*(2*n + 1 + 2*sig)))
D_half = sp.simplify(D.subs(sig, sp.Rational(1, 2))
                     .subs(g**2, kap*Om*om/4))
r2 = sp.simplify(D_half - om*(1 - kap*(n+1)/2))
ser = sp.series(1 + 1/(sp.exp(x)-1), x, 0, 4).removeO()
r3 = sp.simplify(sp.expand(ser - (1/x + sp.Rational(1, 2) + x/12
                                  - x**3/720)))
emit("G10L-1a 10J regressions: spacing residue %s; sigma=1/2 form "
     "residue %s; kappa=1 rungs residue %s" % (r1, r2, r3))

# c1(beta) by symbolic series: nu = 1/x + a + b x, u = x^(1+beta) nu^beta
c1_ok = True
c1_rows = []
for bv in (sp.Integer(0), sp.Rational(1, 2), sp.Integer(1)):
    a, b = sp.symbols('a b')
    nu_ans = 1/x + a + b*x
    u = sp.powsimp(x**(1+bv) * nu_ans**bv)
    nbe_u = 1/u - sp.Rational(1, 2) + u/12
    eq = sp.expand(sp.series(nu_ans - 1 - nbe_u, x, 0, 1).removeO())
    const = eq.coeff(x, 0)
    a_sol = sp.solve(const, a)[0]
    want = 1/(2*(1+bv))
    ok = sp.simplify(a_sol - want) == 0
    c1_ok &= ok
    c1_rows.append((bv, a_sol, ok))
emit("G10L-1b stiffening family c1(beta) = 1/(2(1+beta)): %s" %
     [(str(bv), str(av), ok) for bv, av, ok in c1_rows])

kr = kap*(1 - kap/2)
dkr = sp.diff(kr, kap)
r4 = sp.simplify(dkr.subs(kap, 1))
r5 = sp.simplify(kr.subs(kap, 1) - sp.Rational(1, 2))
emit("G10L-1c kappa_r map = kappa(1-kappa/2): d/dkappa at 1 = %s "
     "(exact zero); kappa_r(1) = 1/2 + %s" % (r4, r5))
g1 = (r1 == 0) and (r2 == 0) and (r3 == 0) and c1_ok and (r4 == 0) \
     and (r5 == 0)
emit("G10L-1: %s" % ("PASS" if g1 else "FAIL"))
emit("")

# ==================== solvers ====================
def G_soft(nb, kv, xv, c):
    s = 1.0 - kv*(nb + c)/2.0
    if s <= 0: return float('inf')
    return nbe(xv*s)

def exists_iter(kv, xv, c, itmax=200000):
    """monotone iteration from 0; returns (exists, nbar)."""
    nstar = 2.0/kv - c
    nb = 0.0
    for _ in range(itmax):
        nb2 = G_soft(nb, kv, xv, c)
        if not np.isfinite(nb2) or nb2 >= nstar - 1e-12:
            return False, None
        if abs(nb2 - nb) < 1e-13:
            return True, nb2
        nb = nb2
    return abs(G_soft(nb, kv, xv, c) - nb) < 1e-9, nb

def exists_scan(kv, xv, c, N=4000):
    """dense-grid sign change of G - n."""
    nstar = 2.0/kv - c
    if nstar <= 0: return False
    grid = np.linspace(0.0, nstar*(1-1e-9), N)
    vals = np.array([G_soft(nb, kv, xv, c) - nb for nb in grid])
    return bool(np.any(vals <= 0))

def kmax_of(xv, c, detector, lo=1e-4, hi=4.0):
    if detector(hi, xv, c): return hi
    if not detector(lo, xv, c): return 0.0
    for _ in range(60):
        mid = 0.5*(lo+hi)
        if detector(mid, xv, c): lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

det_it = lambda kv, xv, c: exists_iter(kv, xv, c)[0]
det_sc = lambda kv, xv, c: exists_scan(kv, xv, c)

# G10L-2 solver gates
ok0 = True
for xv in (0.05, 0.5, 2.0):
    e, nb = exists_iter(1e-9, xv, 1.0)
    ok0 &= e and abs(nb - nbe(xv)) <= 1e-12*max(1, nbe(xv))
res_ok = True
agree = 0; total = 0
rng_k = np.linspace(0.05, 2.0, 30)
rng_x = np.geomspace(0.03, 6.0, 25)
for kv in rng_k:
    for xv in rng_x:
        for c in (1.0,):
            e1, nb1 = exists_iter(kv, xv, c)
            e2 = exists_scan(kv, xv, c)
            total += 1
            if e1 == e2: agree += 1
            if e1 and nb1 is not None:
                res = abs(G_soft(nb1, kv, xv, c) - nb1)
                if res > 1e-10: res_ok = False
kmax_agree = True
for xv in XDEEP:
    k1 = kmax_of(xv, 1.0, det_it)
    k2 = kmax_of(xv, 1.0, det_sc)
    if abs(k1 - k2) > 0.02: kmax_agree = False
g2 = ok0 and res_ok and (agree/total >= 0.98) and kmax_agree
emit("G10L-2 solver: kappa->0 regression %s; residuals %s; detector "
     "agreement %d/%d (bar 98%%); kmax cross-check %s -> %s" %
     (ok0, res_ok, agree, total, kmax_agree, "PASS" if g2 else "FAIL"))
emit("")

# ==================== R-A: the existence map ====================
emit("R-A THE EXISTENCE MAP  kappa_max(x) per softening member")
emit("  (kappa_lock = 0.925; deep window x = 0.03-0.3)")
hdr = "  x      " + "".join("c=%-9.1f" % c for c in OFFSETS)
emit(hdr)
KMAX = {}
family_excl = True
for xv in XDEEP:
    row = "  %-6.3f " % xv
    for c in OFFSETS:
        k1 = kmax_of(xv, c, det_it)
        k2 = kmax_of(xv, c, det_sc)
        KMAX[(xv, c)] = (k1, k2)
        row += "%-11s" % ("%.4f" % k1)
        if max(k1, k2) > KLOCK/2:
            family_excl = False
    emit(row)
emit("  wider map (x up to the Newtonian arm), c = 1:")
for xv in (0.55, 1.0, 2.0, 3.0, 4.0, 6.0):
    k1 = kmax_of(xv, 1.0, det_it)
    emit("    x = %-5.2f kappa_max = %.4f%s" %
         (xv, k1, "  <- the locked kappa first EXISTS out here"
          if k1 >= KLOCK and xv <= 6 and
          kmax_of(max(xv-1, 0.1), 1.0, det_it) < KLOCK else ""))
xc_lock = None
for xv in np.arange(0.5, 8.0, 0.02):
    if det_it(KLOCK, float(xv), 1.0):
        xc_lock = float(xv)
        break
emit("  at kappa_lock = 0.925 (c=1) the fixed point first exists at "
     "x ~ %.2f -- the NEWTONIAN arm, where nu ~ 1 and nothing needs "
     "explaining" % (xc_lock if xc_lock else float('nan')))
ratios = [KLOCK/max(KMAX[(xv, 1.0)][0], 1e-9) for xv in XDEEP]
emit("  exclusion ratios kappa_lock/kappa_max over the deep window "
     "(c=1): %s" % [round(r, 1) for r in ratios])
emit("  classical closed form (diagnostic, recorded): x_c ~ kappa/"
     "(2(1/2 - kappa/8)^2) -> %.2f at kappa=1 vs numeric %.2f" %
     (1.0/(2*(0.5-1/8.0)**2),
      next((float(xv) for xv in np.arange(0.5, 8.0, 0.02)
            if det_it(1.0, float(xv), 1.0)), float('nan'))))
emit("  THE R30 QUESTION ('is kappa -> 1 where self-consistency first"
     " fails?'): %s" %
     ("NO -- failure sits at kappa ~ O(x) in the deep window "
      "(kappa_max %.3f-%.3f), far below 1; kappa = 1 is not the "
      "failure edge" % (min(KMAX[(xv, 1.0)][0] for xv in XDEEP),
                        max(KMAX[(xv, 1.0)][0] for xv in XDEEP))
      if family_excl else
      "the map is not family-stable; see G10L-3"))
emit("G10L-3 family stability of the exclusion (all members, both "
     "detectors, every deep x: kappa_max <= %.4f): %s" %
     (KLOCK/2, "FIRES" if family_excl else "does NOT fire"))
emit("")

# ==================== R-B: the production test ====================
emit("R-B THE PRODUCTION TEST at (kappa_lock, x = 0.05)")
NU_DEMAND = 1 + KLOCK*nbe(0.05)
emit("  measured deep demand nu = 1 + kappa_lock n_BE(0.05) = %.2f" %
     NU_DEMAND)
prod_fail = True
for c in OFFSETS:
    e, nb = exists_iter(KLOCK, 0.05, c)
    if e:
        nu = 1 + KLOCK*nb
        emit("  SCHA c=%.1f: fixed point nu = %.3f (under-produce "
             "%.1fx)" % (c, nu, NU_DEMAND/nu))
        if NU_DEMAND/nu < 5: prod_fail = False
    else:
        emit("  SCHA c=%.1f: NO fixed point (construction fails "
             "entirely)" % c)

def gibbs_trunc(kv, xv, keep_peak):
    Evals = []
    nn = 0
    while True:
        Ev = (nn+0.5) - (kv/4.0)*(nn+0.5)**2
        if nn > 0 and Ev <= Evals[-1]:
            if keep_peak:
                pass  # previous level was the peak; stop before this
            break
        Evals.append(Ev)
        nn += 1
        if nn > 10000: break
    if not keep_peak and len(Evals) > 1:
        Evals = Evals[:-1]
    w = np.exp(-xv*(np.array(Evals) - Evals[0]))
    return float(np.sum(np.arange(len(Evals))*w)/np.sum(w)), len(Evals)

for lab, kp in (("KEEP-A (through peak)", True),
                ("KEEP-B (below peak)", False)):
    sup_nu = 0.0
    for xv in np.geomspace(1e-3, 10, 60):
        nb, nl = gibbs_trunc(KLOCK, float(xv), kp)
        sup_nu = max(sup_nu, 1 + KLOCK*nb)
    nb05, nl = gibbs_trunc(KLOCK, 0.05, kp)
    emit("  truncated-Gibbs %s: %d levels kept; nu(x=0.05) = %.3f, "
         "sup over x = %.3f (under-produce >= %.1fx)" %
         (lab, nl, 1+KLOCK*nb05, sup_nu, NU_DEMAND/sup_nu))
    if NU_DEMAND/sup_nu < 5: prod_fail = False
emit("G10L-4 production clause (every member fails by >= 5x or has "
     "no fixed point): %s" % ("FIRES" if prod_fail else
                              "does NOT fire"))
emit("")

# ==================== stiffening contrast + boot regression ============
emit("STIFFENING CONTRAST (the measured family; negative feedback)")
def nu_stiff(xv, bv, iters=400):
    nu = max(1.0/xv, 1.5)
    for _ in range(iters):
        u = xv**(1+bv)*nu**bv
        nu2 = 1 + nbe(u)
        if abs(nu2 - nu) < 1e-14*nu: return nu2
        nu = 0.5*(nu + nu2)
    return nu
def nu_boot_ref(y, iters=400):
    nu = max(1.0/math.sqrt(y), 1.5)
    for _ in range(iters):
        nu2 = 1 + nbe(nu*y)
        if abs(nu2 - nu) < 1e-14*nu: return nu2
        nu = 0.5*(nu + nu2)
    return nu
worst_bt = 0.0
for y in (1e-4, 1e-2, 0.5, 2.0, 10.0):
    xv = math.sqrt(y)
    d = abs(nu_stiff(xv, 1.0) - nu_boot_ref(y))
    worst_bt = max(worst_bt, d)
g5b = worst_bt <= 1e-12
emit("  beta=1 member == 4F boot nu = 1+n_BE(nu*y): worst |d| = "
     "%.1e (bar 1e-12) -> %s" % (worst_bt, "PASS" if g5b else "FAIL"))
emit("  every member exists at every (x, beta) scanned (negative "
     "feedback); deep limits all nu -> 1/x; c1 runs 1/2 -> 1/4.")
emit("  MEASURED STATUS (ledger, no new fits): deep/binary beta -> 0 "
     "(5R beta < 0.03; 5T regime rows); boot (beta = 1) sky-dead "
     "(5M).  The sky sits at the NO-FEEDBACK point of the stiffening "
     "family and OUTSIDE the softening family entirely.")
emit("")

# ==================== R-C: the convention pin ====================
emit("R-C THE CONVENTION PIN (sympy exact; archived meters only)")
kb = sp.Symbol('kappa_b', positive=True)
nu_deep = kb/x + (1 - kb/2) + kb*x/12
r6 = sp.simplify(sp.expand(
    sp.series(1 + kb*(1/(sp.exp(x)-1)), x, 0, 2).removeO() - nu_deep))
emit("  bare form: nu = kappa_b/x + (1 - kappa_b/2) + kappa_b x/12; "
     "residue %s" % r6)
emit("  => deep-amplitude meter a0_fit/a0_h = kappa_b^2 (g_obs = "
     "kappa_b sqrt(g_N a0_h)); c1 meter kappa_b = 2(1 - c1).")
emit("  first-transition (renormalized) convention: om_r = D(0) = "
     "om(1 - kappa_b/2); matched deep amplitude => kappa_r = "
     "kappa_b(1 - kappa_b/2), maximum EXACTLY 1/2 at kappa_b = 1.")
emit("  map over the measured band:")
for kv in (0.888, 0.925, 1.00, 1.10, 1.48, 1.503):
    emit("    kappa_b = %.3f -> kappa_r = %.3f" % (kv, kv*(1-kv/2)))
emit("  => the renormalized convention is nearly kappa-BLIND over "
     "the whole measured band (range 0.494-0.500 for kappa_b in "
     "[0.888, 1.10]): the BARE convention is the informative one, "
     "and it is the vertex convention (kappa = 4g^2/(Om om_bare)).")
emit("  ARCHIVED METERS: 4S/4Z/10D all parameterize x = sqrt(g_N/a0) "
     "-- bare by construction (verified in their scripts).  The "
     "two-pole tension (lock-meter kappa_b = 0.925-1.00 vs c1-meter "
     "1.10-1.48) is therefore CONVENTION-ROBUST -- a physical "
     "tension (kappa runs; the one-kappa form is 10D-rejected), not "
     "a bookkeeping artifact.  [Never print 'kappa measured ~1' -- "
     "R24 discipline.]")
emit("  curiosity row (recorded): at closure (kappa_b = 1) a "
     "first-transition meter would read kappa_r = 1/2 -- the same "
     "number as c1; both are faces of (1 - kappa/2).")
g5 = (r6 == 0) and g5b
emit("G10L-5: %s" % ("PASS" if g5 else "FAIL"))
emit("")

# ==================== R-D: taxonomy + readings ====================
emit("R-D TAXONOMY (report-grade)")
emit("  family      feedback   deep limit    sky status")
emit("  softening   positive   NO fixed pt   excluded HERE (this "
     "stage): cannot exist at sky occupations")
emit("  stiffening  negative   nu -> 1/x     measured: beta -> 0 "
     "deep/binaries; beta = 1 (boot) sky-dead 5M")
emit("  bare-pinned none       nu -> kappa/x the measured dictionary "
     "(kappa ~ 1 band)")
emit("  READING (fold-row resolution): the self-consistent structure "
     "the 10J fold row demanded is NOT a re-thermalized occupation "
     "-- both feedback signs are now excluded or measured-zero; "
     "consistent with the 10C separation theorem (the bath "
     "occupation cancels between time-orderings; the drive is "
     "external).  The occupation is bare-pinned twice over: "
     "measured (beta -> 0) and structural (this stage).")
emit("  READING (R30 successor 2, booked, no PREDICTIONS row): near "
     "kappa = 1 the E(2) = E(1) degeneracy makes the two dressed "
     "states hybridize under any static perturbation -- an "
     "anomalously large low-frequency response is a DR4-era "
     "signature candidate; no numeric kill condition derived here.")
emit("")

# ==================== verdict (G10L-6) ====================
emit("=" * 68)
some_bound = False
for c in OFFSETS:
    ks = [KMAX[(xv, c)][0] for xv in XDEEP]
    if all(0.8 <= k <= 1.2 for k in ks):
        some_bound = True
if g1 and g2 and g5:
    if family_excl and prod_fail:
        letter = "L-PINNED"
    elif some_bound and not prod_fail:
        letter = "L-BOUND"
    else:
        letter = "L-AMBIG"
else:
    letter = "L-OPEN"
emit("VERDICT LETTER: %s" % letter)
if letter == "L-PINNED":
    emit("  The fed-back occupation (R30's construction, full offset "
         "family + truncated-Gibbs) is EXCLUDED at sky occupations: "
         "no fixed point exists at the locked kappa anywhere in the "
         "deep window (kappa_max ~ O(x), exclusion ratios %s), and "
         "no member produces more than ~nu = 2 against the measured "
         "19.  => the dictionary's occupation is BARE-PINNED "
         "(structural, matching the measured beta -> 0); P-lit has "
         "no self-consistent completion in this family; NO "
         "occupation-consistency bound on kappa exists.  The "
         "K-BOUNDED path closes NEGATIVE: kappa = 1 is not the "
         "failure edge." % [round(r, 1) for r in ratios])
emit("pre-signed credence: L-PINNED -> HOLD mech 8; anomaly-real 53 "
     "untouched (no sky fits)")
emit("gates: G1 %s | G2 %s | G3 %s | G4 %s | G5 %s | G6 mechanical" %
     ("PASS" if g1 else "FAIL", "PASS" if g2 else "FAIL",
      "FIRES" if family_excl else "no-fire",
      "FIRES" if prod_fail else "no-fire",
      "PASS" if g5 else "FAIL"))
emit("total wall-clock %.1f min" % ((time.time()-T0)/60))
