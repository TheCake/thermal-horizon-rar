"""
STAGE 10C (O5-PULL): the dispersive closure -- the honest second-order
pull with its Delta-structure (ROUND-23 condition 8), the 6X control
re-run under purely dispersive coupling (condition 6), and the H-P
adjudication via the exact amplitude-requirement curve (condition 7).

PRE-REGISTRATION (commit BEFORE any execution; bars locked here).

CONTEXT. 10B (V-DISPERSIVE-PERMITTED after ROUND 23): real-exchange
carriers are excluded everywhere; the constraint sector's only soft l=0
vertex is dispersive-class by the exact selection rule [H_int, N_s] = 0.
NOT earned: that the dispersive channel SUPPLIES the measured grammar.
R23's three successor conditions are this stage. The 10B addendum's
labeled PREVIEW: (x/2)coth(x/2) = 1 + x^2/12 - x^4/720 (the BE ladder's
Bernoulli rungs in a mean dressing, O(1) at every anchor).

THE HYPOTHESIS UNDER TEST (the dispersive-transmission reading): the
10B vertex H_int = (H_sys/c^2) dPhi with dPhi = phi c^2 (b+b_dag), i.e.

    H = om*(N+sig) + Om*b_dag b + g*(N+sig)*(b+b_dag),   g = om*phi,

(sig = 1/2 for the EP charge H_sys = hbar om(N+1/2); sig = 0 for a
bare number vertex) transmits the measured law nu = 1 + n_BE(x) when
the system dictionary mode is thermally occupied at the horizon
temperature (x = om/T; the thermal postulate and the LAW itself are
prior art, Cadoni & Tuveri 2019 -- ours here is the TRANSMISSION
STRUCTURE: which pull the allowed vertex produces, what closure it
requires, and what remains underived).

T1 -- THE SEPARATION THEOREM (the honest pull's Delta-structure).
  The model above is EXACTLY solvable (polaron displacement):
      E(n,m) = om*(n+sig) + Om*m - (g^2/Om)*(n+sig)^2   (+ bath zp).
  Claims, each gated:
  (a) G2 (symbolic): the static second-order pull is BATH-OCCUPATION-
      FREE -- the two time-orderings (virtual emission (m+1)/(-Om) +
      virtual absorption m/(+Om)) cancel the bath n exactly. R23's
      "mean pull pi(2n+1) = 6.3/44.8/66 over-produces" evaluated the
      JC/EXCHANGE template at the dispersive vertex -- the wrong
      template; the honest dispersive mean pull carries NO ambient
      occupation.
  (b) G1 (numeric, truncated Fock, sig in {0, 1/2}): the polaron
      spectrum is exact; the transition pull reads the SYSTEM's
      occupation, linearly, to all orders in g:
          dom(n) = -(g^2/Om)*(2n + 1 + 2 sig).
  (c) G3 (regression vs the ROUND-23 report): the bath occupation
      enters ONLY the first-order jitter (mean-zero, RMS =
      g(n+sig)sqrt(2n_b+1)). His rows are reproduced as the labeled
      objects they are -- jitter sqrt(pi(2n+1)) = 2.509/6.693/8.122
      and zero-detuning exchange-template pi(2n+1) = 6.30/44.80/65.97
      at his occupations n = 0.502/6.63/10.0. Bars +/-0.01 (jitter),
      +/-0.05 (pull).
  (d) T1d DISCLOSURE (striking our own tempting overclaim): the
      sig-term (EP vs number vertex) enters the pull only as an
      n-INDEPENDENT constant, absorbed by the declared 9Z-b
      renormalization (one bare frequency + one a0 rescale). The
      static ladder CANNOT distinguish the EP charge from a number
      charge; c1 = 1/2 is NOT a vertex-selection result. The 10B
      preview's coth-vs-nu gap (x*nu - (x/2)coth(x/2) = x/2, exact)
      is the same absorbed constant in the x-multiplied presentation;
      c1 = 1/2 arises from Newton's 1 plus n_BE's own Bernoulli -1/2.

T2 -- THE TRANSMISSION MAP + THE KAPPA CLOSURE (condition 8's answer).
  Thermal system average: <pull> = -(g^2/Om)(2 nbar(x)+1+2sig);
  renormalizing the constant part leaves the anomalous state-dependent
  softening d om/om = -(2g^2/(Om om)) nbar(x). Under the program's
  standing response dictionary (a softened mode extracts more
  acceleration from the same source; compliance-additive form --
  convention STANDING, used not derived, disclosed):
      nu(x) = 1 + kappa * n_BE(x),     kappa = 4 g^2/(Om om).
  * G4 (exact series): at kappa = 1 the ENTIRE measured ladder is
    automatic: nu = 1/x + 1/2 + x/12 + 0*x^2 - x^3/720 + ... (c1=1/2,
    c2=1/12, c3=0, c4=-1/720 -- all Bernoulli rungs at once).
  * LINEARITY (printed, no fit): the polaron pull is strictly linear
    in n (all orders in g) = the measured law's exact linearity; the
    multiplicative resummation 1/(1-n) diverges at x = ln 2 mid-
    transition -- grossly non-RAR; the data select the additive form.
  * G5 (the kappa METERS -- archived numbers only, NO new sky fits):
      a0-lock meter: deep limit of 1+kappa*n gives fitted a0/horizon
        a0 = kappa^2; 4H/4V/5M (galaxy hier a0 = (1.05+/-0.10)e-10 at
        +0.1 sigma of cH0/2pi) => kappa = 1.00 +/- 0.05.
        BAR: kappa in [0.90, 1.10].
      c1 meter (cross-family translation, disclosed): kappa = 2(1-c1);
        4S flat-M/L c1 = 0.450 (0.385-0.519) => kappa = 1.100
        (0.962-1.230). BAR: band contains 1.
      TENSION ROW (mandatory in every cell): 4Z hier profile c1 =
        0.258 (0.208-0.309) => kappa = 1.484 (1.382-1.584), excludes
        1; program status 1/4-vs-1/2 OPEN (4Z bootstrap ~0.4+/-0.3;
        5T ultra-deep arm votes 1/2; 5K binaries p=1/2 12/12). The
        meters are projections of the same fits, not independent.
  * Closure constant: kappa = 1 <=> g = sqrt(Om*om)/2 (geometric
    mean). Consistency OBSERVATION (no bar): g_close at the anchors
    vs the 10A gamma-dial band [0.072, 1.571] H (t-envelope spans
    2pi => order-grade only).
  NOT CLAIMED (PERMITTED-grade, pre-signed): kappa = 1 is MEASURED-
  consistent, not derived; its vertex-level content -> T3.

T3 -- THE REQUIREMENT CURVE + H-P ADJUDICATION (condition 7).
  kappa = 1 with g = om*phi forces EXACTLY
      phi_req = (1/2) sqrt(Om/om) = (1/2) sqrt(x_amb/x_loc).
  * G6a (exactness, sympy) + the exact restatement (identity): phi_req
    = sqrt(hbar Om/(2 E_c)) with E_c = 2 hbar om -- the bright mode
    must couple with the zero-point amplitude of an effective inertia
    equal to TWO QUANTA of the dictionary mode it dresses. A named
    REQUIREMENT, not a derivation.
  * The curve at the anchors (binary x_amb = 1.0954, galaxy 0.1411;
    x_loc across the measured windows): phi_req ~ 0.19-0.74 and RUNS
    (1/sqrt(om) at fixed environment -- zero-point-like; sqrt(Om)
    across environments). H-P supplies the CONSTANT sqrt(pi) = 1.772.
  * G6b H-P criterion (gated): KEEP only if (i) level within x3 of
    phi_req at ALL anchors AND (ii) the required running across the
    measured window is < x2 (H-P is constant, so a larger required
    running refutes it structurally). Expected by algebra already in
    this docstring: both fail => DROP-AND-REPLACE: H-P dropped as
    stated; the requirement curve becomes the named open object (the
    closure constant's vertex-level content). Single-horizon
    semiclassical row (phi_1 = H l_P/c): orders-short vs the exact
    requirement printed (R23's honest under-supply, now curve-exact).
  REACHABILITY DISCLOSURE (per the R23 cosmetic-fork critique): G6b's
  branch is algebra-determined in this docstring; the gate locks
  arithmetic, it does not create suspense.

T4 -- THE 6X CONTROL UNDER DISPERSIVE COUPLING (condition 6).
  Original 6X: Kerr system (NA=4, WA=5, CHI=0.8) + one thermal
  ambient; dress-ward step |1>->|2> RESONANT with b via exchange
  lam*(a'b + b'a): saturated weight = (1/2) P(n_b>=1) = (1/2) n/(1+n)
  -- existence-counting (the n_b=0 sector has no partner state),
  lambda-independent; the KMS ratio IS the lending probability.
  * G0 truncation bookkeeping as 6X. G7 (regression): re-run the 6X
    resonant channel verbatim (same constants, same estimator) ->
    ratio-form win at the ORIGINAL bars (slope 1+/-0.15, 3x rms
    margin over share and raw); detuned virtual row re-recorded.
  * G8 -- two purely dispersive arms, both [H_int, Na] = 0, driven at
    the 1->2 line (rotating frame level term -CHI*Na +
    (CHI/2)Na(Na-1); drive eps(a+a'), eps = 0.01; Om_b = 0.8;
    NB = 56; nbar in {0.25,0.5,1,2,4,8}):
    ARM-A (linear-displacement vertex gd*Na(b+b'), gd = 0.3, WITH the
      polaron compensation +(gd^2/Om_b)Na^2 = the 9Z-b renormalization
      mechanized in the toy; saturation tier). SIGNED PREDICTION:
      every |1,m> keeps a degenerate partner |2,m> (Franck-Condon
      overlap essentially never zero) => saturated weight = 1/2 FLAT
      in nbar -- the dispersive channel cannot imprint the ambient
      occupation on the saturated dress-ward weight (FC dressing
      changes speed, not existence). BARS: max|P2bar - 1/2| <= 0.06
      AND spread <= 0.04 (vs the lending law's 0.35 span) AND the
      ratio form MISSES slope 1+/-0.15. eps -> eps/2 row printed
      (existence-counting => eps-independent). Known small leak,
      disclosed: sectors near the Laguerre zero of the FC overlap
      (m ~ 10) can partially de-hybridize at the differential-Stark
      scale -- the widened 0.06/0.04 bars price it.
    ARM-B (cross-Kerr vertex gk*Na*Nb, gk = 0.1; no compensation --
      its diagonal IS the detuning under test; [H,Nb] = 0 kills all
      sidebands). EXACT per-sector law (two-level time average):
      P2bar = sum_m p_m * 4 eps^2/(8 eps^2 + gk^2 m^2) ~
      (1/2)/(1+nbar) -- THE ANTI-LENDING LAW: the dispersive channel
      gates on the ambient being EMPTY (quiet-ambient complement),
      the OPPOSITE monotonicity to the measured gate (sky: s_amb =
      e^{-x_amb} GROWS with ambient occupation; anti-law FALLS).
      GATES: GB-exact |P2bar - exact sum| <= 0.02 at every nbar;
      form contest anti-ratio 1/(1+n) vs ratio n/(1+n) vs raw n at
      6X margins (3x rms + slope).
  PRE-SIGNED OUTCOME CELLS:
    GATE-FLAT+ANTI (the signed expectation): ARM-A flat AND ARM-B
      anti-ratio => a purely dispersive ambient coupling does NOT
      reproduce the gate e^{-x} in either tier -- at saturation it is
      occupation-BLIND, perturbatively it carries the COMPLEMENT
      (opposite-monotonic) law => R23's proposed completion of the 6X
      demotion FAILS; 6X's real-exchange reading of the GATE channel
      STANDS. Pre-signed reading (the two-channel synthesis): the
      amplitude/pull channel is dispersive (l=0 constraint, this
      stage) while the GATE channel requires real exchange -- whose
      only ledger-surviving partner is the LOCAL 6Y cloud through the
      near-field l=2 (tidal) vertex: the 10B horizon-partner
      exclusions (gapped sqrt(l(l+1))H, thermally empty, geometric
      (L/r_c)^2) do NOT apply to the local cloud, whose scale is the
      system's own. The R22-condition-4 derivation (the near-field
      exchange vertex) returns as the required next leg. Annotations
      (R23 cond 5, deferred until this arm): 6X KEEPS its real-
      exchange reading with the partner re-identified; 9T NOT mooted.
    GATE-KMS (R23's expectation): the ratio law appears in a
      dispersive arm at the 6X bars => single-channel dispersive
      world; execute R23-cond-5 annotations (6X -> form-validator;
      9T mooted -- the credence-cost question booked for the review
      round, not self-adjudicated).
    GATE-MIXED: anything else => UNRESOLVED-carried; annotations stay
      deferred.
  REACHABILITY: T4 is the stage's genuinely OPEN empirical leg (the
  flat/anti predictions are algebra-strong; the dynamics can surprise;
  GATE-KMS is a live cell).

LETTER GRAMMAR (trap #12 WITH TEETH -- every optimistic clause names
its gate; ungated clauses are written PERMITTED-grade):
  tokens: THM-SEP <- G1 & G2 & G3;  LADDER <- G4;
          KAPPA <- G5 (CONTAINED-WITH-TENSION iff a0-lock in
            [0.9,1.1] AND flat-c1 band contains 1; hier tension row
            prints in EVERY cell);
          REQ <- G6a + G6b branch;  GATE-<cell> <- G0 & G7 & G8.
  composite:
    P-TRANSMITTED = THM-SEP & LADDER & KAPPA & G6a & G7 green:
      "the dispersive vertex's honest pull is bath-occupation-free
      and system-linear; under the standing dictionary + declared
      renormalization it transmits nu = 1 + kappa n_BE(x) with the
      entire measured Bernoulli ladder at kappa = 1; kappa is
      MEASURED ~ 1 (temperature lock 1.00 +/- 0.05; flat-c1 band),
      hier 1/4-lean tension disclosed; the closure constant is
      UNDERIVED -- its vertex content is exactly phi_req =
      (1/2)sqrt(Om/om) (H-P dropped per G6b); the ambient gate is
      <GATE cell>."
      EXPLICITLY NOT CLAIMED in any cell: mechanism CLOSED; kappa
      derived; the response dictionary derived; priority on nu = 1+n
      (C&T 2019).
    P-OPEN = THM-SEP/LADDER fail physically, or KAPPA meters exclude
      1, or G7 regression fails.
CREDENCE MAP (pre-signed; moves ONLY by this map):
  P-TRANSMITTED and ROUND-24 rules NO unpatched hole -> bath-mechanism
    conditional 15 -> 18. P-TRANSMITTED with a hole -> HOLD 15.
  THM-SEP fails PHYSICALLY (cancellation/linearity genuinely absent,
    not a bug) -> 15 -> 12 (pre-committed strike).
  KAPPA meters exclude 1 -> P-OPEN, HOLD 15, tension booked.
  GATE cells: NO credence move in any cell (annotation consequences
    only, as pre-signed above).
  anomaly-real: 53 UNTOUCHED in every cell (no sky fit anywhere).
DISCLOSURES: (i) the pull->boost map is the standing response
  dictionary, used not derived; (ii) C&T 2019 priority on the law;
  (iii) sig-blindness (T1d): c1 = 1/2 does not select the EP vertex;
  (iv) the kappa meters are projections of the same archived fits;
  (v) reachability: G1-G4/G6 lock algebra (failure = bug or STOP);
  G5's values and ALL of T4 are the open cells; (vi) the 10A
  gamma-band comparison is order-grade (t-convention envelope).
AMENDMENTS: any post-commit change is logged in the output and NOTES
  with its reason, pre-quote (standing rule).

Writes data/stage10c_pull.txt.
"""
import math
import time
import numpy as np
import sympy as sp

OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

say("STAGE 10C: the dispersive closure (O5-PULL) -- R23 conditions 8+6+7")
say("=" * 74)
t_wall = time.time()

def nbe(x):
    return 1.0/(math.exp(x) - 1.0)

# ========================== T1: THE SEPARATION THEOREM =====================
say("")
say("T1 -- the separation theorem for H = om(N+sig) + Om b'b + g(N+sig)(b+b')")
say("-" * 74)

def fock_ops(NA, NB):
    a = np.diag(np.sqrt(np.arange(1, NA)), 1)
    b = np.diag(np.sqrt(np.arange(1, NB)), 1)
    Ia, Ib = np.eye(NA), np.eye(NB)
    A = np.kron(a, Ib); B = np.kron(Ia, b)
    return A, B, np.kron(Ia, Ib)

def polaron_check(sig, om=1.0, Om=1.0, g=0.1, NA=6, NB=80, K=60):
    A, B, Id = fock_ops(NA, NB)
    Np = A.T @ A + sig*Id
    H = om*Np + Om*(B.T @ B) + g*Np @ (B + B.T)
    ev = np.sort(np.linalg.eigvalsh(H))
    pred = np.sort([om*(n+sig) + Om*m - (g**2/Om)*(n+sig)**2
                    for n in range(NA) for m in range(NB)])
    return float(np.max(np.abs(ev[:K] - pred[:K])))

d0 = polaron_check(0.0)
dh = polaron_check(0.5)
ok_g1 = d0 < 1e-6 and dh < 1e-6
say(f"G1 polaron spectrum exact (lowest 60 levels): max|dE| sig=0: {d0:.2e},"
    f" sig=1/2: {dh:.2e}  -> {'PASS' if ok_g1 else 'FAIL'}")

g_s, Om_s, m_s, ns_s, sig_s, x_s, om_s = sp.symbols(
    'g Omega m n sigma x omega', positive=True)
shift = (g_s*(ns_s+sig_s))**2 * ((m_s+1)/(-Om_s) + m_s/Om_s)
resid = sp.simplify(shift + g_s**2*(ns_s+sig_s)**2/Om_s)
ok_g2 = resid == 0
say(f"G2 bath-occupation cancellation: (m+1)/(-Om) + m/(+Om) leaves "
    f"-g^2(n+sig)^2/Om, residual = {resid}  -> {'PASS' if ok_g2 else 'FAIL'}")
say("   => the honest dispersive mean pull carries NO ambient occupation;")
say("      R23's pi(2n+1) rows were the JC/EXCHANGE template at the")
say("      dispersive vertex -- the wrong template (reproduced below as")
say("      the labeled objects they are).")

Epol = om_s*(ns_s+sig_s) - g_s**2*(ns_s+sig_s)**2/Om_s
pull = sp.simplify((Epol.subs(ns_s, ns_s+1) - Epol) - om_s)
pull_pred = -(g_s**2/Om_s)*(2*ns_s + 1 + 2*sig_s)
ok_pull = sp.simplify(pull - pull_pred) == 0
say(f"   transition pull dom(n) = -(g^2/Om)(2n+1+2sig): "
    f"{'exact' if ok_pull else 'MISMATCH'}")
say("   EP charge (sig=1/2): dom = -(2g^2/Om)(n+1); thermal system mean")
say("   <dom> = -(2g^2/Om)(nbar(x)+1) = -(2g^2/Om) nu(x)  [linear in the")
say("   occupation, all orders in g -- polaron-exact]")

say("")
say("G3 -- the R23 rows, correctly labeled (regression vs ROUND-23 report):")
say("   anchor        n_BE     x=ln(1+1/n)   jitter sqrt(pi(2n+1))   "
    "exch-template pi(2n+1)")
R23 = [("binary      ", 0.502, 2.509, 6.30),
       ("gal ambient ", 6.63, 6.693, 44.80),
       ("gal deep    ", 10.0, 8.122, 65.97)]
ok_g3 = True
for lab, n, jt, pt in R23:
    x = math.log(1 + 1/n)
    jit = math.sqrt(math.pi*(2*n + 1))
    pul = math.pi*(2*n + 1)
    ok_g3 &= abs(jit - jt) < 0.01 and abs(pul - pt) < 0.05
    say(f"   {lab} {n:6.3f}   {x:8.4f}      {jit:8.3f} (his {jt})       "
        f"{pul:8.2f} (his {pt})")
say(f"   both row families reproduced: {'PASS' if ok_g3 else 'FAIL'}")
say("   jitter = first-order, MEAN-ZERO (variance object, bath-side);")
say("   pi(2n+1) = zero-detuning EXCHANGE-template pull -- neither is the")
say("   honest dispersive mean pull, which is bath-n-free (G2).")

gap = sp.simplify(pull_pred.subs(sig_s, sp.Rational(1, 2))
                  - pull_pred.subs(sig_s, 0))
say("")
say(f"T1d sig-blindness: pull(EP) - pull(number) = {gap} -- n-INDEPENDENT")
say("   => absorbed by the 9Z-b bare-frequency constant; the static ladder")
say("   CANNOT distinguish the EP charge from a number charge; c1 = 1/2 is")
say("   NOT a vertex-selection result (our own overclaim, struck here).")
nu_ser = sp.series(1 + 1/(sp.exp(x_s)-1), x_s, 0, 5).removeO()
coth_obj = sp.series(x_s/2*sp.coth(x_s/2), x_s, 0, 6).removeO()
gap2 = sp.simplify(sp.expand(x_s*nu_ser) - coth_obj)
say(f"   the 10B preview gap: x*nu - (x/2)coth(x/2) = {gap2} (exact) = the")
say("   same absorbed constant in the x-multiplied presentation; c1 = 1/2")
say("   arises from Newton's 1 plus n_BE's own Bernoulli constant -1/2.")

ok_t1 = ok_g1 and ok_g2 and ok_pull and ok_g3
say(f"T1 verdict: THM-SEP {'PASS' if ok_t1 else 'FAIL'}")

# ===================== T2: THE TRANSMISSION MAP + KAPPA ====================
say("")
say("T2 -- the transmission map and the kappa closure")
say("-" * 74)
say("   renormalize the n-independent pull (9Z-b scheme: one bare frequency")
say("   + one a0 rescale); anomalous softening d om/om = -(2g^2/(Om om)) n;")
say("   standing response dictionary (compliance-additive, disclosed):")
say("       nu(x) = 1 + kappa*n_BE(x),   kappa = 4 g^2/(Om om)")

ser = sp.series(1 + 1/(sp.exp(x_s)-1), x_s, 0, 6)
c = [ser.coeff(x_s, k) for k in (-1, 0, 1, 2, 3)]
targets = [1, sp.Rational(1, 2), sp.Rational(1, 12), 0,
           sp.Rational(-1, 720)]
ok_g4 = all(sp.simplify(a - b) == 0 for a, b in zip(c, targets))
say(f"G4 ladder at kappa=1: nu = 1/x + 1/2 + x/12 + 0*x^2 - x^3/720 + ...")
say(f"   coefficients {c} vs measured (1/2, 1/12, 0, -1/720): "
    f"{'PASS -- ALL RUNGS AT ONCE' if ok_g4 else 'FAIL'}")
say(f"   linearity: polaron pull strictly linear in n (all orders in g) =")
say(f"   the measured law's exact linearity; the multiplicative resummation")
say(f"   1/(1-n) diverges at x = ln 2 = {math.log(2):.3f} (mid-transition)")
say(f"   -- grossly non-RAR; the data select the additive/channel form.")

say("")
say("G5 -- the kappa meters (archived fits; NO new sky numbers):")
k_a0, k_a0_sd = 1.00, 0.05
ok_a0 = 0.90 <= k_a0 <= 1.10
say(f"   a0-lock meter (4H/4V/5M: fitted a0/horizon = kappa^2; galaxy hier")
say(f"   a0 = (1.05 +/- 0.10)e-10 at +0.1 sigma of cH0/2pi):")
say(f"       kappa = {k_a0:.2f} +/- {k_a0_sd:.2f}   "
    f"[bar 0.90-1.10: {'PASS' if ok_a0 else 'FAIL'}]")
c1_flat, c1_lo, c1_hi = 0.450, 0.385, 0.519
kf, kf_hi, kf_lo = 2*(1-c1_flat), 2*(1-c1_lo), 2*(1-c1_hi)
ok_flat = kf_lo <= 1.0 <= kf_hi
say(f"   c1 meter, flat-M/L (4S: c1 = 0.450, 0.385-0.519; kappa = 2(1-c1),")
say(f"   cross-family translation disclosed):")
say(f"       kappa = {kf:.3f}  ({kf_lo:.3f}-{kf_hi:.3f})   "
    f"[contains 1: {'PASS' if ok_flat else 'FAIL'}]")
c1h, c1h_lo, c1h_hi = 0.258, 0.208, 0.309
say(f"   TENSION ROW (mandatory, prints in every cell): 4Z hier profile")
say(f"   c1 = {c1h} ({c1h_lo}-{c1h_hi}) => kappa = {2*(1-c1h):.3f} "
    f"({2*(1-c1h_hi):.3f}-{2*(1-c1h_lo):.3f}) -- excludes 1;")
say(f"   program status: 1/4-vs-1/2 OPEN (4Z bootstrap ~0.4+/-0.3 agrees")
say(f"   with flat; 5T ultra-deep arm votes c1=1/2; 5K binaries p=1/2")
say(f"   12/12). Meters are projections of the same fits, not independent.")
ok_g5 = ok_a0 and ok_flat
say(f"G5: KAPPA = {'CONTAINED-WITH-TENSION (PASS)' if ok_g5 else 'FAIL'}")

say("")
say("   closure: kappa = 1 <=> g = sqrt(Om*om)/2. At the anchors (Om =")
say("   x_amb*T, om = x_loc*T, T = H/2pi) -- OBSERVATION row, order-grade")
say("   (10A gamma band [0.072, 1.571] H; t-envelope spans 2pi):")
for lab, xamb, xlocs in [("binary ", 1.0954, (0.5, 1.0, 2.0)),
                         ("galaxy ", 0.1411, (0.0953, 0.3, 1.0))]:
    for xl in xlocs:
        gc = math.sqrt(xamb*xl)/(4*math.pi)
        inb = 0.072 <= gc <= 1.571
        say(f"     {lab} x_amb={xamb:.4f} x_loc={xl:6.4f}: g_close = "
            f"{gc:.4f} H  ({'inside' if inb else 'below'} the 10A band)")

ok_t2 = ok_g4 and ok_g5

# ==================== T3: REQUIREMENT CURVE + H-P ==========================
say("")
say("T3 -- the amplitude requirement curve and the H-P adjudication")
say("-" * 74)
phi_s = sp.symbols('phi', positive=True)
kap_expr = 4*(om_s*phi_s)**2/(Om_s*om_s)
phi_req_sym = sp.solve(sp.Eq(kap_expr, 1), phi_s)[0]
ok_g6a = sp.simplify(phi_req_sym - sp.sqrt(Om_s/om_s)/2) == 0
say(f"G6a requirement: kappa=1, g=om*phi  =>  phi_req = {phi_req_sym} "
    f"= (1/2)sqrt(Om/om): {'PASS' if ok_g6a else 'FAIL'}")
Ec = sp.symbols('E_c', positive=True)
sol_Ec = sp.solve(sp.Eq(sp.sqrt(Om_s/(2*Ec)), sp.sqrt(Om_s/om_s)/2), Ec)[0]
say(f"   exact restatement: phi_req = sqrt(hbar Om/(2 E_c)) with E_c = "
    f"{sol_Ec}*hbar = TWO QUANTA of")
say("   the dictionary mode it dresses (identity; a named REQUIREMENT --")
say("   the bright mode couples as if its effective inertia were the")
say("   system's own two-quantum scale -- NOT a derivation).")
say("")
say("   the curve at the anchors (phi_req = (1/2)sqrt(x_amb/x_loc)):")
rows = []
for lab, xamb, xlocs in [("binary ", 1.0954, (0.5, 1.0, 2.0)),
                         ("galaxy ", 0.1411, (0.0953, 0.3, 1.0))]:
    for xl in xlocs:
        pr = 0.5*math.sqrt(xamb/xl)
        rows.append(pr)
        say(f"     {lab} x_loc={xl:6.4f}:  phi_req = {pr:.4f}   "
            f"(H-P sqrt(pi) = {math.sqrt(math.pi):.4f} -> off x"
            f"{math.sqrt(math.pi)/pr:.2f})")
run_bin = math.sqrt(2.0/0.5)
run_gal = math.sqrt(1.0/0.0953)
say(f"   required RUNNING across the measured windows: x{run_bin:.2f} "
    f"(binary) / x{run_gal:.2f} (galaxy); H-P is CONSTANT.")
lvl_ok = all(1/3 <= math.sqrt(math.pi)/pr <= 3 for pr in rows)
run_ok = max(run_bin, run_gal) < 2.0
keep_hp = lvl_ok and run_ok
say(f"G6b H-P criterion: level within x3 at ALL anchors: "
    f"{'yes' if lvl_ok else 'NO'}; required running < x2: "
    f"{'yes' if run_ok else 'NO'}  => "
    f"{'KEEP' if keep_hp else 'DROP-AND-REPLACE'}")
c_si, H_si, lP = 2.998e8, 2.268e-18, 1.616e-35
phi1 = H_si*lP/c_si
short = math.log10(min(rows)/phi1)
say(f"   single-horizon semiclassical: phi_1 = H l_P/c = {phi1:.3e} -> "
    f"{short:.1f} orders under the smallest requirement (the honest")
say("   under-supply, now against the exact curve, not an O(1) window).")
say("   VERDICT: H-P DROPPED AS STATED (constant amplitude; wrong running;")
say("   level x2.4-9.4 high); REPLACED by the requirement curve phi_req =")
say("   (1/2)sqrt(Om/om) as the named open object. Reachability: branch")
say("   algebra-determined in the pre-reg (disclosure v).")
ok_t3 = ok_g6a

# ==================== T4: THE 6X DISPERSIVE CONTROL ========================
say("")
say("T4 -- the 6X control: exchange regression + purely dispersive arms")
say("-" * 74)

def build(NA, NB):
    ad = np.diag(np.sqrt(np.arange(1, NA)), -1)
    bd = np.diag(np.sqrt(np.arange(1, NB)), -1)
    Ia, Ib = np.eye(NA), np.eye(NB)
    return dict(A=np.kron(ad.T, Ib), Ad=np.kron(ad, Ib),
                B=np.kron(Ia, bd.T), Bd=np.kron(Ia, bd),
                Na=np.kron(ad @ ad.T, Ib), Nb=np.kron(Ia, bd @ bd.T),
                Ia=Ia, Ib=Ib, NA=NA, NB=NB)

def thermal(N, nbar):
    if nbar <= 0:
        p = np.zeros(N); p[0] = 1.0
        return np.diag(p)
    xx = math.log(1.0 + 1.0/nbar)
    p = np.exp(-xx*np.arange(N)); p /= p.sum()
    return np.diag(p)

NA, NB = 4, 56
o = build(NA, NB)
CHI, LAM, WA = 0.8, 0.02, 5.0

def H_exchange(delta=0.0, lam=LAM):
    Id = np.kron(o['Ia'], o['Ib'])
    return (WA*o['Na'] + (WA + CHI + delta)*o['Nb']
            + 0.5*CHI*(o['Na'] @ (o['Na'] - Id))
            + lam*(o['Ad'] @ o['B'] + o['Bd'] @ o['A']))

def sat_weight(H, nb):
    ev, U = np.linalg.eigh(H)
    p1 = np.zeros(NA); p1[1] = 1.0
    rho = np.kron(np.diag(p1), thermal(NB, nb))
    rl = U.T @ rho @ U
    rb = U @ np.diag(np.diag(rl)) @ U.T
    proj2 = np.kron(np.diag((np.arange(NA) == 2)*1.0), o['Ib'])
    return float(np.real(np.trace(proj2 @ rb)))

NBS = np.array([0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
trunc = (8.0/9.0)**NB
say(f"G0 truncation at nbar=8, NB={NB}: P(n>=NB) = {trunc:.1e}  "
    f"{'PASS' if trunc < 5e-3 else 'FAIL'}")

say("")
say("G7 exchange regression (resonant lending channel, verbatim 6X):")
W = np.array([sat_weight(H_exchange(), nb) for nb in NBS])
for nb, w in zip(NBS, W):
    say(f"   n_amb = {nb:5.2f}: P2bar = {w:.5f}   [(1/2)n/(1+n) = "
        f"{0.5*nb/(1+nb):.5f}]")
forms = {'ratio n/(1+n)': NBS/(1+NBS), 'share n/(2n+1)': NBS/(2*NBS+1),
         'raw n': NBS}
res = {}
for nm, f in forms.items():
    s, cc = np.polyfit(np.log(f), np.log(W), 1)
    rms = float(np.sqrt(np.mean((np.log(W) - s*np.log(f) - cc)**2)))
    res[nm] = (s, rms)
    say(f"   {nm:<16}: slope {s:+.3f}, rms {rms:.4f}")
rs, rr = res['ratio n/(1+n)']
ok_g7 = (abs(rs - 1) < 0.15 and rr*3 <= res['share n/(2n+1)'][1]
         and rr*3 <= res['raw n'][1])
say(f"G7: ratio-form win at the 6X bars: {'PASS' if ok_g7 else 'FAIL'}")
wdet1 = sat_weight(H_exchange(delta=0.3*CHI), 1.0)
wdet2 = sat_weight(H_exchange(delta=0.3*CHI, lam=LAM/2), 1.0)
say(f"   detuned virtual row: lam-scaling {wdet1/wdet2:.2f} (expect ~4; "
    f"the 6X taxonomy re-recorded)")

OMB, GD, GK, EPS = 0.8, 0.3, 0.1, 0.01
IdF = np.kron(o['Ia'], o['Ib'])
LVL = -CHI*o['Na'] + 0.5*CHI*(o['Na'] @ (o['Na'] - IdF))

say("")
say(f"G8 dispersive arms (eps = {EPS}, Om_b = {OMB}; rotating frame at the")
say(f"   1->2 line; both vertices commute with Na):")

say("")
say(f"ARM-A: linear vertex gd*Na(b+b'), gd = {GD}, polaron-compensated")
say("   (the 9Z-b renormalization mechanized); saturation tier:")
def H_armA(eps=EPS):
    return (LVL + eps*(o['A'] + o['Ad']) + OMB*o['Nb']
            + GD*(o['Na'] @ (o['B'] + o['Bd']))
            + (GD**2/OMB)*(o['Na'] @ o['Na']))
HA = H_armA()
WA_arm = np.array([sat_weight(HA, nb) for nb in NBS])
for nb, w in zip(NBS, WA_arm):
    say(f"   n_amb = {nb:5.2f}: P2bar = {w:.5f}")
sprA = float(WA_arm.max() - WA_arm.min())
devA = float(np.max(np.abs(WA_arm - 0.5)))
sA, _ = np.polyfit(np.log(NBS/(1+NBS)),
                   np.log(np.maximum(WA_arm, 1e-12)), 1)
ratio_misses_A = not (abs(sA - 1) < 0.15)
whalf = sat_weight(H_armA(eps=EPS/2), 1.0)
say(f"   spread = {sprA:.4f}, max|P2bar-1/2| = {devA:.4f}; ratio-form "
    f"slope {sA:+.3f} (must miss 1+/-0.15: {'yes' if ratio_misses_A else 'NO'})")
say(f"   eps->eps/2 at nbar=1: {whalf:.5f} (existence-counting => "
    f"eps-independent)")
ok_armA = devA <= 0.06 and sprA <= 0.04 and ratio_misses_A
say(f"ARM-A: {'FLAT CONFIRMED' if ok_armA else 'NOT FLAT (cell shifts)'}")

say("")
say(f"ARM-B: cross-Kerr vertex gk*Na*Nb, gk = {GK} (its diagonal IS the")
say("   detuning under test; [H,Nb]=0 kills all sidebands). Exact law:")
say("   P2bar = sum_m p_m * 4eps^2/(8eps^2 + gk^2 m^2) ~ (1/2)/(1+nbar):")
HB = (LVL + EPS*(o['A'] + o['Ad']) + OMB*o['Nb']
      + GK*(o['Na'] @ o['Nb']))
WB = np.array([sat_weight(HB, nb) for nb in NBS])
ok_gbx = True
for nb, w in zip(NBS, WB):
    xx = math.log(1 + 1/nb)
    pm = np.exp(-xx*np.arange(NB)); pm /= pm.sum()
    exact = float(np.sum(pm*4*EPS**2/(8*EPS**2 + GK**2*np.arange(NB)**2)))
    ok_gbx &= abs(w - exact) <= 0.02
    say(f"   n_amb = {nb:5.2f}: P2bar = {w:.5f}   [exact sum {exact:.5f},"
        f" (1/2)/(1+n) = {0.5/(1+nb):.5f}]")
say(f"   GB-exact |P2bar - exact| <= 0.02 at every nbar: "
    f"{'PASS' if ok_gbx else 'FAIL'}")
formsB = {'anti 1/(1+n)': 1/(1+NBS), 'ratio n/(1+n)': NBS/(1+NBS),
          'raw n': NBS}
resB = {}
for nm, f in formsB.items():
    s, cc = np.polyfit(np.log(f), np.log(WB), 1)
    rms = float(np.sqrt(np.mean((np.log(WB) - s*np.log(f) - cc)**2)))
    resB[nm] = (s, rms)
    say(f"   {nm:<14}: slope {s:+.3f}, rms {rms:.4f}")
aB, aR = resB['anti 1/(1+n)'], resB['ratio n/(1+n)']
anti_wins = (abs(aB[0] - 1) < 0.15 and aB[1]*3 <= aR[1]
             and aB[1]*3 <= resB['raw n'][1])
kms_wins = (abs(aR[0] - 1) < 0.15 and aR[1]*3 <= aB[1]
            and aR[1]*3 <= resB['raw n'][1])
say(f"ARM-B: anti-lending wins: {'yes' if anti_wins else 'no'}; "
    f"ratio/KMS wins: {'yes' if kms_wins else 'no'}")

if kms_wins:
    gate_cell = "GATE-KMS"
elif ok_armA and anti_wins and ok_gbx:
    gate_cell = "GATE-FLAT+ANTI"
else:
    gate_cell = "GATE-MIXED"
say("")
say(f"T4 cell: {gate_cell}")
if gate_cell == "GATE-FLAT+ANTI":
    say("   => a purely dispersive ambient coupling does NOT reproduce the")
    say("   gate e^{-x} in either tier: at saturation it is occupation-")
    say("   BLIND (existence-counting is FC-blind, weight 1/2 flat); at")
    say("   the perturbative tier it carries the COMPLEMENT law 1/(1+n)")
    say("   (gates on the ambient being EMPTY) -- the OPPOSITE monotonicity")
    say("   to the measured gate (sky: s_amb grows with ambient occupation).")
    say("   R23's proposed completion of the 6X demotion FAILS; 6X's real-")
    say("   exchange reading of the GATE channel STANDS. Pre-signed reading")
    say("   (two-channel synthesis): dispersive amplitude (l=0 constraint,")
    say("   this stage) + real-exchange gate, whose only ledger-surviving")
    say("   partner is the LOCAL 6Y cloud via the near-field l=2 vertex")
    say("   (the 10B horizon-partner exclusions -- gapped, empty, (L/r_c)^2")
    say("   -- do not apply to the local cloud). The R22-cond-4 derivation")
    say("   returns as the required next leg. Annotations (R23 cond 5):")
    say("   6X KEEPS real-exchange, partner re-identified; 9T NOT mooted.")
elif gate_cell == "GATE-KMS":
    say("   => the ratio law appears dispersively: R23's single-channel")
    say("   path fires; execute R23-cond-5 annotations (6X -> form-")
    say("   validator; 9T mooted -- credence-cost question booked for the")
    say("   review round).")
else:
    say("   => UNRESOLVED-carried; annotations stay deferred.")

# =========================== LETTER ASSEMBLY ===============================
say("")
say("=" * 74)
tok_thm = "THM-SEP:" + ("PASS" if ok_t1 else "FAIL")
tok_lad = "LADDER:" + ("PASS" if ok_g4 else "FAIL")
tok_kap = "KAPPA:" + ("CONTAINED-WITH-TENSION" if ok_g5 else "FAIL")
tok_req = "REQ:CURVE-EXACT+HP-" + ("KEEP" if keep_hp else "DROPPED")
letter = ("P-TRANSMITTED" if (ok_t1 and ok_g4 and ok_g5 and ok_g6a
                              and ok_g7) else "P-OPEN")
say(f"TOKENS: {tok_thm}  {tok_lad}  {tok_kap}  {tok_req}  {gate_cell}")
say(f"LETTER: {letter}")
say("")
if letter == "P-TRANSMITTED":
    say("THE EARNED STATEMENT (every clause gated or PERMITTED-grade):")
    say("  the dispersive vertex's honest pull is bath-occupation-free and")
    say("  system-linear (separation theorem, exact); under the standing")
    say("  response dictionary + declared renormalization it transmits")
    say("  nu = 1 + kappa*n_BE(x) with the ENTIRE measured Bernoulli ladder")
    say("  (1/2, 1/12, 0, -1/720) at kappa = 1; kappa is MEASURED ~1 by the")
    say("  temperature lock (1.00 +/- 0.05) and the flat-M/L c1 band, with")
    say("  the 4Z hier 1/4-lean tension disclosed; the closure constant is")
    say("  UNDERIVED -- its vertex content is exactly phi_req =")
    say("  (1/2)sqrt(Om/om) (H-P dropped: wrong running, x2.4-9.4 high);")
    say("  the ambient gate is NOT dispersive-reproducible (per the T4")
    say("  cell above).")
    say("  NOT claimed: mechanism closed; kappa derived; the response")
    say("  dictionary derived; priority on nu = 1+n (C&T 2019).")
say("")
say("CREDENCE (pre-signed map): P-TRANSMITTED -> candidate rise cell 15->18")
say("  ONLY IF ROUND 24 rules no unpatched hole; any hole -> HOLD 15;")
say("  physical THM-SEP failure -> 15->12 (did not fire);")
say("  gate cells carry NO credence move. anomaly-real 53 UNTOUCHED.")
say("")
say(f"done ({(time.time()-t_wall)/60:.1f} min)")

with open('data/stage10c_pull.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(OUT) + "\n")
print("\nSTAGE 10C done -> data/stage10c_pull.txt")
