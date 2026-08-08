"""
STAGE 10B -- O5-CARRIER: THE CARRIER LEDGER AND THE CONSTRAINT VERTEX
(round-22 conditions 4 + 5-residual; author call 2026-08-08 "Let's
chase it. All you got.").

QUESTION: 10A excluded the radiative graviton sector as the carrier of
the measured O(H) coupling (U = 7e-45/1e-34, confinement-closed) and
round 22 demanded the positive half: DERIVE the constraint-sector
coupling as a Rabi-capable real-exchange channel compatible with 6X
(real exchange, virtual = null) and 9U (Rabi phase) -- or identify the
actual channel and compute its normalization. This stage builds the
COMPLETE carrier ledger: every candidate soft partner a bound system
could exchange with, each computed to its ceiling, plus the one vertex
gravity's constraint sector actually provides.

THE CANDIDATES (the ledger rows):
  C1 radiative free gravitons ..... DEAD (10A: 34-44 orders; banked)
  C2 confined/cavity radiative .... DEAD (10A addendum: 5-11 orders)
  C3 matter ELASTIC modes ......... this stage T1 (frequency-forbidden)
  C4 matter SECULAR modes ......... this stage T2 (dictionary-mismatch)
  C5 horizon soft mode, monopole .. this stage T3 (THE vertex: exists,
                                    O(1) amplitude -- but dispersive
                                    CLASS, exchange-forbidden by an
                                    exact selection rule)
  C6 local-expansion mode delta-H . constraint-slaved to local matter
                                    (Friedmann) -> reduces to C4
  C7 other systems' soft modes .... monopole-monopole G hbar w1 w2 /
                                    (c^4 R) ~ 1e-80 H -- one line, dead

INSTRUMENTS:
  T1 THE EPICYCLIC DICHOTOMY + THE CLOUD-CONTINUUM CLOSURE. sympy: for
     a Kepler orbit in a weak isotropic external tide T_ext,
     Omega_tot^2 = GM/r^3 + T and kappa^2 = (1/r^3) d(r^4
     Omega_tot^2)/dr = Omega_K^2 + 4T EXACTLY, so the matter-side
     frequency ladder is {m kappa + k Omega} with minimum nonzero
     splitting |Omega - kappa| = (3/2) T/Omega_K + O(T^2) (the secular
     rate) and everything else >= kappa ~ Omega_orb. THE DICHOTOMY:
     matter modes are either ELASTIC (>= the measured orbital floors:
     SPARC min Omega/H = 30.2, binaries ~5e3 -- ledger mech-9q-frozen)
     or SECULAR (the near-degenerate combinations = C4). The 6Y
     cloud's internal continuum is therefore DOUBLY dead as a leak
     (R22 condition 5-residual): detuning >= 29 H AND thermal
     occupation n_BE(2 pi * 30.2) = e^-189.7 -- compute the admixture
     bound exactly (log-form).
  T2 THE SECULAR-CARRIER EXCLUSION. The C4 frequencies CAN reach ~H
     (computed at the anchors: the Galactic-tide apsidal rate of a
     10 kAU binary is ~2.5 H -- the one matter family in the soft
     band). But the DICTIONARY kills it: omega_sec = (3/2)
     Omega_g^2/Omega_orb runs as s^{+3/2} while omega_dict =
     sqrt(g_N a0)/c runs as s^{-1}: the ratio sweeps s^{5/2} ~ 10^6.0
     across the measured 0.2-50 kAU sample, and the implied effective
     a0 (if x tracked the secular frequency, a0_eff ~ ratio^2) sweeps
     ~12 orders where the measured law holds with ONE a0 to ~10%
     (5M temperature lock; v7/3U fits interior at the physical
     field). Galaxy-side co-read: an axisymmetric-disk star has no
     comparable slow tide at all. C4 is excluded AS THE CARRIER
     (sub-leading admixture not excluded -- labeled).
  T3 THE CONSTRAINT VERTEX (the round-22 demand, answered with a
     different sign than asked):
     (a) UNIQUENESS+FORM: an l=0 potential fluctuation sourced outside
         a region is CONSTANT inside it (interior Laplace, sympy:
         the regular l=0 solution of Laplace's equation is const) --
         so the horizon's soft l=0 coordinate couples to a bound
         system ONLY through its total energy: H_int = (H_sys/c^2)
         delta-Phi. The charge is equivalence-principle-forced
         (energy gravitates; the system's soft quantum carries
         hbar omega_s/c^2 -- model-independent). Horizon l >= 1
         modes: gapped at sqrt(l(l+1)) H (R16) AND (L/r_c)^l
         geometry-suppressed -- dead (computed).
     (b) CLASS (the selection rule, exact): H_int proportional to
         H_sys commutes with the system number operator, [H_int,
         N_s] = 0 (sympy, operator algebra on the oscillator) =>
         the vertex is DISPERSIVE/longitudinal (JC-dispersive
         a-dagger-a form), NOT exchange. Gravity's constraint
         sector provides NO off-diagonal soft vertex. Combined with
         C1-C4/C7 dead: THE REAL-EXCHANGE REALIZATION OF THE
         GRAMMAR HAS NO CARRIER anywhere in physics as enumerated.
     (c) SINGLE-PORT => M=1 MECHANIZED: the energy couples to ONE
         operator (the sum of microstate displacements) => the 9W
         bright-mode theorem applies EXACTLY (not analogically):
         the S microstates present exactly one collective mode to
         any monopole-coupled system -- the R16 patch ("one soft
         collective coordinate; the microstates are the reservoir")
         acquires its mechanism.
     (d) AMPLITUDE (hypothesis H-P, stated as the one load-bearing
         assumption): each horizon microstate contributes an
         independent Planck-grade potential quantum delta-Phi_1 =
         c H l_P (one Planck length of horizon displacement). Then
         the bright-mode amplitude is EXACT BY IDENTITY: delta-Phi_B
         = sqrt(S) c H l_P = sqrt(pi) c^2 (since S = pi r_c^2/l_P^2
         and r_c = c/H; sympy-exact). The linearity cap (delta-Phi
         <= c^2) is saturated at O(1) -- the collective zero-point
         sits AT the geometric-linearity edge, not orders past it.
     (e) THE AMPLITUDE CONDITION: the fractional dressing of the
         system's soft frequency is then delta-omega/omega =
         sqrt(pi * (2 n_h + 1))-grade = O(1), occupation-structured
         exactly as the 6H JC-dispersive anchor (lam^2(2n+1)/Delta
         selects the zero-point share 1/(2 nu - 1)) -- the O(1)
         grammar amplitude the whole O5-NORM arc demanded, carried
         by the constraint sector WITHOUT real exchange. Semiclassical
         per-microstate amplitude (no sqrt-S collectivization) would
         under-supply by sqrt(S) ~ 5e60 -- the bracket is printed.
     (f) CONSISTENCY SWEEP (each row computed/cited): the dispersive
         realization must carry every measured grammar element
         without real exchange -- the gate s^L = e^{-Lx} is a
         Boltzmann/KMS state-statistic (6U derived it as detailed
         balance, no borrowing needed); the local factor 1/(2nu-1)
         is the dispersive zero-point share (6H); the p-family's r
         re-reads as a structural weight with the dispersive-fixed
         point r = 1/2 sitting -0.55 sigma from the measured
         r-hat = 0.39 +/- 0.20 (9U/9V, archived); the 9T detuning
         unresolvability and 9W reduction are state-level and
         survive; the 6X toy remains the FORM-validator (its real
         exchange was the toy's realization, not the sky's); ROUND
         19's own retraction ("tail VIRTUAL" -> the locality data
         re-read) already pointed this direction. ANNOTATION LIST
         (if the letter fires): the resonance-arc INTERPRETIVE layer
         (9T "resonance", 9U "swing count", 9S detuning-bound
         readings) demotes to effective-parameter language; rows
         annotated, never deleted; the L3 ceiling-1/2 becomes the
         zero-point share + exchange-symmetric fixed point (5P
         beta = 1/2 reading), not a resonant transfer ceiling.

PRE-REGISTERED BARS (locked before any computation):
  B-T1: the cloud-continuum admixture bound <= 1e-30 x g^2 (log-exact
        arithmetic) AND the epicyclic identity kappa^2 = Omega_K^2 +
        4T exact => the clause "cloud continuum closed / elastic
        carrier frequency-forbidden" is licensed.
  B-T2: the secular/dictionary ratio sweeps >= 3 orders across the
        measured binary s-range (expect ~6) AND the implied a0 sweep
        >= 6 orders (expect ~12) vs the measured ~10% constancy =>
        "secular carrier excluded" licensed. The anchor secular rates
        are REPORTED (one lands ~H -- the honest near-miss stated).
  B-T3: (i) interior-constancy exact; (ii) [H_int, N_s] = 0 exact;
        (iii) the sqrt(S) l_P = sqrt(pi) r_c identity exact; (iv) the
        amplitude condition lands in [0.1, 10] (O(1) window) under
        H-P. Each clause of the letter carries its own gate (trap
        #12).
LETTER GRAMMAR (trap #11/#12 -- every gate AND every positive clause
wired):
  STOP-class (bug): GB-1 epicyclic identity, GB-3a interior
    constancy, GB-3b selection rule, GB-3c sqrt-pi identity fail =>
    abort, no verdict.
  V-DISPERSIVE-FORCED: B-T1 AND B-T2 AND B-T3(i-iv) all fire => the
    carrier ledger is EMPTY for real exchange (C1-C7 all dead) AND
    the constraint sector supplies the dispersive vertex at O(1)
    amplitude under H-P => the grammar's realization is
    DISPERSIVE-FORCED: the mechanism's exchange layer is excluded,
    its dispersive layer acquires a derived vertex + amplitude with
    ONE named hypothesis (H-P) and one exact identity. The
    consistency sweep (f) rows print; the annotation list executes at
    booking (rows annotated, interpretive layers demoted).
  V-EMPTY: B-T1/B-T2 fire but B-T3(iv) lands outside [0.1, 10] =>
    every carrier AND the dispersive amplitude fail => the literal
    mechanism layer (both realizations) loses its last candidate =>
    pre-signed 15 -> 12 (fourth realization-class strike: rates
    6K/6M/6N, oscillator-exchange + dispersive-amplitude here); the
    equilibrium/phenomenological family survives as the carrier of
    the anomaly program (unaffected: the measured grammar, the
    E-ladder, T1' first law, the dissolved FATAL corner).
  V-QUOTE: any rival alive (B-T1 or B-T2 misfire) or a mixed cell =>
    quote everything, no reframe executed.
CREDENCE MAP (pre-signed): V-DISPERSIVE-FORCED -> HOLD 15 pending
  ROUND 23 (rise cell: R23 no-unpatched-hole AND V-DISPERSIVE-FORCED
  -> 15 -> 18 -- the constructive cell the hunt was for; any hole ->
  hold + conditions). V-EMPTY -> 15 -> 12 immediately (pre-signed
  strike; the map does NOT wait for the round -- the 6K precedent).
  V-QUOTE -> HOLD 15. anomaly-real 53 UNTOUCHED in every cell (no
  sky fits anywhere; T2 uses ledger-archived facts).

DISCLOSURE (the full in-session design trail, none verdict-grade):
  this stage's shape CHANGED during design -- the original plan
  (derive a Rabi-capable near-field exchange coupling, per R22
  condition 4 as literally worded) hit an exact wall: the
  energy-monopole vertex commutes with N_s (dispersive class), and
  the in-session scratch triage killed every exchange candidate
  (mechanical oscillators at soft frequencies do not exist except
  secular ones; the secular dictionary mismatches; per-quantum
  horizon-monopole exchange at semiclassical amplitude is
  Planck-dead ~60 orders; even coherent-enhanced variants fall ~19
  orders short). The sqrt(S)-collectivization insight (9W applied to
  the microstate ensemble; sqrt(S) l_P = sqrt(pi) r_c identity
  noticed in-session) converts the same vertex to O(1) amplitude in
  the DISPERSIVE channel. Scratch magnitudes disclosed: binary
  Galactic-tide apsidal rate ~2.5 H at 10 kAU (the C4 near-miss);
  secular/dictionary ratio ~s^{5/2} (~6 orders across the sample);
  the amplitude condition lands
  at sqrt(pi (2n+1)) ~ 2.5-8 for n = 0.5-10. Expected letter:
  V-DISPERSIVE-FORCED (stated so the firing is not read as a
  surprise); the genuinely open cells are the exact coefficients
  and whether any gate breaks. NOT computed pre-commit: the exact
  epicyclic series coefficient (expected -3/2), the exact admixture
  log-bound, the exact a0-sweep, the T3 numeric rows.
Output: data/stage10b_carrier.txt. Wall-clock: ~1-3 min (sympy).
"""
import math
import os
import time

import sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

OUT = []
def P(s=""):
    print(s, flush=True)
    OUT.append(s)

def save():
    with open('data/stage10b_carrier.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

t00 = time.time()
P("=" * 78)
P("STAGE 10B -- O5-CARRIER: THE CARRIER LEDGER AND THE CONSTRAINT "
  "VERTEX")
P("=" * 78)
gates = {}

TWO_PI = 2*math.pi
H_SI = 2.27e-18; C = 2.998e8; GSI = 6.674e-11
HBAR = 1.0546e-34
MSUN = 1.989e30; AU = 1.496e11; KPC = 3.0857e19
A0 = 1.05e-10

# =====================================================================
# PART 1 -- T1: the epicyclic dichotomy + the cloud-continuum closure
# =====================================================================
P("")
P("PART 1 -- T1 THE EPICYCLIC DICHOTOMY")

r, GM, T = sp.symbols('r GM T', positive=True)
Om2 = GM/r**3 + T                      # isotropic external tide
kap2 = sp.simplify(sp.diff(r**4*Om2, r)/r**3)
ok_kap = sp.simplify(kap2 - (GM/r**3 + 4*T)) == 0
# secular splitting: Omega - kappa to leading order in T
eps = sp.Symbol('epsilon', positive=True)
OmK = sp.sqrt(GM/r**3)
split = sp.series(sp.sqrt(GM/r**3 + eps) - sp.sqrt(GM/r**3 + 4*eps),
                  eps, 0, 2).removeO()
coef = sp.simplify(split/eps)
ok_split = sp.simplify(coef + sp.Rational(3, 2)/OmK) == 0
gates['GB-1'] = ok_kap and ok_split
P("GB-1 epicyclic identity: kappa^2 = Omega_K^2 + 4T %s; secular "
  "splitting Omega - kappa = -(3/2) T/Omega_K + O(T^2) %s => %s"
  % (ok_kap, ok_split, "PASS" if gates['GB-1'] else "FAIL"))
if not gates['GB-1']:
    P("STOP: epicyclic identity failed (bug-class).")
    save(); raise SystemExit(0)
P("  THE DICHOTOMY: matter-side frequencies are combinations "
  "{m kappa + k Omega}; every combination is either ELASTIC "
  "(>= kappa ~ Omega_orb: measured floors Omega/H = 30.2 SPARC / "
  "~5e3 binaries, ledger mech-9q-frozen) or SECULAR (the "
  "near-degenerate |Omega - kappa| family = candidate C4). The band "
  "(secular, elastic-floor) contains NO matter modes.")

# the cloud-continuum leak bound (R22 condition 5-residual)
x_floor = TWO_PI*30.2
ln_n = -x_floor                       # ln n_BE ~ -x for x >> 1
det = 30.2 - 0.035                    # detuning in H units (>= 29 H)
ln_adm = math.log(4.0) + ln_n - 2*math.log(det)
gates['B-T1'] = (ln_adm < math.log(1e-30))
P("  cloud-continuum closure: admixture <= 4 g^2 n(2 pi 30.2)/"
  "(29 H)^2; ln n = %.1f => admixture <= 10^%.1f x g^2 "
  "(bar 1e-30) -> %s  [the cloud's internal modes are Wien-dead AND "
  "far-detuned; the leak cannot re-open inside the hosting sector]"
  % (ln_n, ln_adm/math.log(10),
     "PASS" if gates['B-T1'] else "FAIL"))

# =====================================================================
# PART 2 -- T2: the secular-carrier exclusion
# =====================================================================
P("")
P("PART 2 -- T2 THE SECULAR CARRIER (C4)")

OM_GAL_ROT = 2*math.pi/(230e6*3.156e7)     # Galactic rotation rate
def om_orb(s_m, M_kg): return math.sqrt(GSI*M_kg/s_m**3)
def om_sec(s_m, M_kg): return 1.5*OM_GAL_ROT**2/om_orb(s_m, M_kg)
def om_dict(s_m, M_kg):
    gN = GSI*M_kg/s_m**2
    return math.sqrt(gN*A0)/C

P("  anchor rates (M_tot = 1 M_sun):")
for s_kau in (0.2, 2.0, 10.0, 50.0):
    s_m = s_kau*1e3*AU
    osec = om_sec(s_m, 2*MSUN*0.5)
    odic = om_dict(s_m, 2*MSUN*0.5)
    P("    s = %5.1f kAU : omega_sec = %.2e /s (%.2e H)   omega_dict "
      "= %.2e /s (%.2e H)   ratio %.2e"
      % (s_kau, osec, osec/H_SI, odic, odic/H_SI, osec/odic))
s_lo, s_hi = 0.2e3*AU, 50e3*AU
ratio_lo = om_sec(s_lo, MSUN)/om_dict(s_lo, MSUN)
ratio_hi = om_sec(s_hi, MSUN)/om_dict(s_hi, MSUN)
sweep = abs(math.log10(ratio_hi/ratio_lo))
a0_sweep = 2*sweep
gates['B-T2'] = (sweep >= 3.0 and a0_sweep >= 6.0)
P("  the dictionary mismatch: omega_sec/omega_dict runs s^{5/2} -> "
  "sweeps %.1f orders across the measured 0.2-50 kAU sample; the "
  "implied a0_eff (if x tracked the secular frequency) sweeps %.1f "
  "orders vs the measured ~10%% constancy (5M lock; v7/3U interior "
  "fits) -> %s"
  % (sweep, a0_sweep, "PASS" if gates['B-T2'] else "FAIL"))
P("  honest near-miss stated: the 10-kAU Galactic-tide apsidal rate "
  "IS ~H-scale (%.1f H) -- the secular family is the ONE matter "
  "family in the soft band; it dies by dictionary, not by frequency. "
  "Galaxy-side co-read: an axisymmetric-disk star carries no "
  "comparable slow tide (the family barely exists there). Sub-leading "
  "secular ADMIXTURE is not excluded -- labeled; the carrier is."
  % (om_sec(10e3*AU, MSUN)/H_SI))
P("  C6 (local delta-H mode): constraint-slaved to local matter by "
  "the Friedmann constraint -> reduces to C4 -> excluded with it.")
c7 = GSI*HBAR*(0.174*H_SI)*(0.035*H_SI)/(C**4*2.5e20*H_SI)
P("  C7 (system-system monopole-monopole, R = 8 kpc): g_c ~ "
  "G hbar w1 w2/(c^4 R) = %.1e H -- dead." % c7)

# =====================================================================
# PART 3 -- T3: the constraint vertex
# =====================================================================
P("")
P("PART 3 -- T3 THE CONSTRAINT VERTEX")

# (a) interior constancy of the l=0 exterior-sourced mode (sympy):
# regular l=0 solution of Laplace inside a source-free ball:
# Phi'' + (2/r) Phi' = 0 -> Phi = C1 + C2/r; regularity kills C2.
Phi = sp.Function('Phi')
sol = sp.dsolve(sp.diff(Phi(r), r, 2) + 2*sp.diff(Phi(r), r)/r,
                Phi(r))
basis = sol.rhs.free_symbols - {r}
const_part = sol.rhs.subs(dict.fromkeys(
    [s_ for s_ in basis if sp.limit(sp.diff(sol.rhs, s_), r, 0)
     in (sp.oo, -sp.oo, sp.zoo)], 0))
ok_const = sp.diff(const_part, r) == 0
gates['GB-3a'] = ok_const
P("GB-3a interior constancy: regular l=0 Laplace solution inside a "
  "source-free region = %s (constant; the 1/r branch removed by "
  "regularity) -> %s  => the horizon's soft l=0 coordinate couples "
  "to a bound system ONLY through its total energy: H_int = "
  "(H_sys/c^2) deltaPhi  [charge = energy: equivalence-principle-"
  "forced; the soft quantum's charge hbar omega_s/c^2 is "
  "model-independent]"
  % (const_part, "PASS" if gates['GB-3a'] else "FAIL"))

# (b) the selection rule: [H_int, N_s] = 0 (oscillator algebra)
a_op = sp.Symbol('a', commutative=False)
ad_op = sp.Symbol('ad', commutative=False)
b_op = sp.Symbol('b', commutative=False)
bd_op = sp.Symbol('bd', commutative=False)
def normal(expr):
    # reduce using [a, ad] = 1, [b, bd] = 1, a-sector commutes with b
    e = sp.expand(expr)
    changed = True
    while changed:
        changed = False
        e2 = e.subs(a_op*ad_op, ad_op*a_op + 1)
        e2 = sp.expand(e2.subs(b_op*bd_op, bd_op*b_op + 1))
        # move b-sector left of a-sector (they commute)
        e2 = sp.expand(e2.subs(a_op*bd_op, bd_op*a_op)
                       .subs(a_op*b_op, b_op*a_op)
                       .subs(ad_op*bd_op, bd_op*ad_op)
                       .subs(ad_op*b_op, b_op*ad_op))
        if e2 != e:
            e = e2; changed = True
    return e
N_s = ad_op*a_op
H_int_op = ad_op*a_op*(b_op + bd_op)     # (H_s/c^2) deltaPhi form
comm = normal(H_int_op*N_s - N_s*H_int_op)
gates['GB-3b'] = (comm == 0)
P("GB-3b the selection rule: [(a+a)(b+bd)-form energy-monopole "
  "vertex, N_s] = %s -> %s  => the vertex is DISPERSIVE-class "
  "(JC-longitudinal a†a(b+b†)); it CANNOT transfer system quanta -- "
  "gravity's constraint sector provides NO off-diagonal soft vertex"
  % (comm, "PASS" if gates['GB-3b'] else "FAIL"))
# contrast row: an exchange vertex does NOT commute (displayed)
H_ex = ad_op*b_op + bd_op*a_op
comm_ex = normal(H_ex*N_s - N_s*H_ex)
P("  contrast (exchange vertex a†b + b†a): [H_ex, N_s] = %s != 0 -- "
  "the class distinction is exact" % comm_ex)

# (c) single-port -> 9W: stated (the energy couples to ONE operator =
# the sum of microstate displacements = the bright-mode projection)
P("  single-port: the system's energy couples to ONE operator (the "
  "collective sum of microstate displacements) => the 9W bright-mode "
  "theorem applies EXACTLY: the S microstates present ONE collective "
  "mode to any monopole-coupled system (the R16 'one soft "
  "coordinate; microstates = reservoir' patch, mechanized).")

# (d) the sqrt(pi) identity (exact)
S_sym, rc, lP, c_s, H_sym = sp.symbols('S r_c l_P c H', positive=True)
ident = sp.sqrt(sp.pi*rc**2/lP**2)*c_s*H_sym*lP
ident = sp.simplify(ident.subs(rc, c_s/H_sym))
ok_id = sp.simplify(ident - sp.sqrt(sp.pi)*c_s**2) == 0
gates['GB-3c'] = ok_id
P("GB-3c the amplitude identity: sqrt(S) x (c H l_P) with S = pi "
  "r_c^2/l_P^2 and r_c = c/H  =  %s  == sqrt(pi) c^2 -> %s"
  % (ident, "PASS" if ok_id else "FAIL"))
if not (gates['GB-3a'] and gates['GB-3b'] and gates['GB-3c']):
    P("STOP: a T3 exactness gate failed (bug-class).")
    save(); raise SystemExit(0)
P("  HYPOTHESIS H-P (the one load-bearing assumption, named): each "
  "of the S = %.2e horizon microstates contributes an independent "
  "Planck-grade potential quantum deltaPhi_1 = c H l_P = %.2e m^2/s^2 "
  "(one Planck length of horizon displacement). Under 9W the bright "
  "mode then carries deltaPhi_B = sqrt(pi) c^2 EXACTLY -- the "
  "collective zero-point saturates the geometric-linearity edge "
  "(deltaPhi ~ c^2), an O(1) statement, not orders past it."
  % (math.pi*(C/H_SI)**2/(1.616e-35)**2, C*H_SI*1.616e-35))

# (e) the amplitude condition
P("")
P("  THE AMPLITUDE CONDITION (B-T3 iv): fractional dressing of the "
  "system's soft frequency, sqrt(pi (2 n + 1))-grade:")
amps = {}
for nm, n_amb in (('binary (n_amb = 0.502)', 0.502),
                  ('galaxy ambient (n_amb = 6.63)', 6.63),
                  ('galaxy deep local (n ~ 10)', 10.0)):
    amp = math.sqrt(math.pi*(2*n_amb + 1))
    amps[nm] = amp
    P("    %-28s delta-omega/omega ~ %.2f" % (nm, amp))
ok_amp = all(0.1 <= v <= 10.0 for v in amps.values())
gates['B-T3iv'] = ok_amp
P("  O(1) window [0.1, 10]: %s  -- the constraint vertex at the "
  "collective amplitude delivers ORDER-UNITY, occupation-structured "
  "((2n+1)) dressing = the amplitude the O5-NORM arc demanded, in "
  "the DISPERSIVE channel (exactly the 6H JC anchor's lam^2(2n+1)/"
  "Delta structure whose zero-point share 1/(2nu-1) the grammar "
  "measures)" % ("PASS" if ok_amp else "FAIL"))
P("  semiclassical bracket (no collectivization): per-microstate "
  "amplitude alone under-supplies by sqrt(S) = %.1e -- the same "
  "60-order gap 10A's C5 scratch found; the sqrt(S) bright-mode "
  "factor IS the bridge, and it is the program's own banked 9W "
  "theorem doing the bridging." % math.sqrt(math.pi*(C/H_SI)**2
                                            / (1.616e-35)**2))
# horizon l >= 1 co-read
P("  horizon l >= 1 co-read: gapped at sqrt(l(l+1)) H (R16) AND "
  "(L_sys/r_c)^l geometry-suppressed: galaxy (10 kpc/horizon) = "
  "%.1e, binary %.1e -- dead both ways."
  % (10*KPC/(C/H_SI), 1e4*AU/(C/H_SI)))

# (f) the consistency sweep (computed/cited rows)
P("")
P("  CONSISTENCY SWEEP -- the dispersive realization vs every "
  "measured grammar element:")
P("    gate s^L = e^-Lx: Boltzmann/KMS state-statistic (6U derived "
  "it AS detailed balance; no real borrowing required) ... CARRIES")
P("    local factor 1/(2nu-1): the dispersive zero-point share (6H "
  "JC anchor -- dispersive BY CONSTRUCTION) ................ CARRIES")
P("    p-family weight r: dispersive-fixed r = 1/2 (the zero-point "
  "share at resonance-language face value) sits %.2f sigma from "
  "the measured r-hat = 0.390 +/- 0.198 (9U) ............... CARRIES"
  % ((0.5 - 0.390)/0.198))
P("    9T detuning unresolvability + 9W reduction: state-level "
  "theorems, realization-independent ....................... CARRY")
P("    6X toy: FORM-validator (its real exchange was the toy's "
  "realization; R19 already retracted 'tail VIRTUAL' the other "
  "direction) ............................................. RE-READS")
P("    M=1 (6Y measured): mechanized by (c) single-port ....... GAINS")

# =====================================================================
# PART 4 -- the verdict (locked grammar)
# =====================================================================
P("")
led_dead = gates['B-T1'] and gates['B-T2']
if led_dead and gates['B-T3iv']:
    P("==> 10B VERDICT (locked grammar): V-DISPERSIVE-FORCED -- the "
      "carrier ledger is EMPTY for real exchange (C1/C2 dead by 10A; "
      "C3 frequency-forbidden; C4 dictionary-excluded; C6 reduces to "
      "C4; C7 dead by 80 orders; and the one vertex gravity's "
      "constraint sector provides is exchange-FORBIDDEN by the exact "
      "selection rule [H_int, N_s] = 0). The grammar's realization "
      "is therefore DISPERSIVE: the constraint vertex at the "
      "9W-collective amplitude (sqrt(S) l_P = sqrt(pi) r_c, exact) "
      "delivers O(1), (2n+1)-structured dressing under the single "
      "named hypothesis H-P -- the amplitude condition the O5-NORM "
      "arc demanded, met WITHOUT real exchange. The round-22 "
      "condition-4 demand is answered with the opposite sign: a "
      "Rabi-capable constraint coupling does not exist, and the "
      "measured grammar never needed one -- every measured element "
      "(the KMS gate, the zero-point share, r = 1/2 at -0.55 sigma) "
      "is state-statistics or dispersive-class. The exchange "
      "INTERPRETIVE layer (9T resonance language, 9U swing-count, "
      "9S detuning-bound readings, 6X's realization) demotes to "
      "effective/form-validator status -- rows annotated at booking, "
      "never deleted.")
elif led_dead:
    P("==> 10B VERDICT (locked grammar): V-EMPTY -- every exchange "
      "carrier AND the dispersive amplitude fail. The literal "
      "mechanism layer loses its last candidate. Pre-signed strike: "
      "bath-mechanism conditional 15 -> 12 (fourth realization-class "
      "strike). Surviving: the equilibrium/phenomenological family, "
      "the measured grammar, the E-ladder, the first law, the "
      "dissolved FATAL corner.")
else:
    P("==> 10B VERDICT (locked grammar): V-QUOTE -- a rival carrier "
      "cell misfired; quote everything, no reframe executed. Cells: "
      "B-T1=%s B-T2=%s B-T3iv=%s"
      % (gates['B-T1'], gates['B-T2'], gates['B-T3iv']))
P("    CREDENCE: per the pre-signed map (V-DISPERSIVE-FORCED -> HOLD "
  "15 pending ROUND 23, rise cell 15 -> 18 on no-hole; V-EMPTY -> "
  "15 -> 12 immediately; V-QUOTE -> HOLD 15). anomaly-real 53 "
  "untouched in every cell (no sky fits).")

# ledger legs
import csv
need = {'mech-10a-dprov': False, 'mech-9w-multimode': False,
        'mech-9q-frozen': False, 'gal-9u-ptail': False}
for row in csv.reader(open('LEDGER.csv', encoding='utf-8')):
    if row and row[0] in need and row[1] in ('CURRENT', 'CO-QUOTED'):
        need[row[0]] = True
gates['GB-6'] = all(need.values())
P("")
P("GB-6 ledger legs: %s -> %s"
  % (need, "PASS" if gates['GB-6'] else "FAIL"))
P("")
P("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
  for k, v in sorted(gates.items())))
P("")
P("done (%.1f min)" % ((time.time() - t00)/60))
save()
print("\nsaved: data/stage10b_carrier.txt")
