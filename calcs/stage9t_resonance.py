"""
STAGE 9T -- THE RESONANCE THEOREM (O5-DETUNING answered at toy grade).

THE QUESTION (ROUND-16 A9, the unpatched hole; 9S bounded it at
r >= 0.315 from data): why would the system's dressing mode and the
ambient collective mode be IN TUNE? Off resonance the exchange weight
is r = (1/2) * g_c^2/(g_c^2 + delta^2) < 1/2, generically
field-dependent -- yet the sky measures the 1/2.

THE ANSWER PROPOSED HERE (three legs; the first two exact, the third
a witnessed premise):

T1 (THE PHASE BUDGET, exact): by Lemma 1 (9Q) every mode in the
anomaly sector has omega = x*H/(2*pi) with x in [0, 1] -- the WHOLE
sector spans H/(2*pi). Any two soft modes are detuned by
delta = dx*H/(2*pi) with dx <= 1. No system has existed longer than
t_max = eta/H (eta ~ 1; age*H0 = 0.99 in LCDM). Therefore the
accumulated detuning phase obeys

    phi = delta * t  <=  eta * dx / (2*pi)  <=  eta/(2*pi) ~ 0.16 rad.

EPOCH-INVARIANT: at redshift z the sector width is H(z)/(2*pi) and
the available time is ~1/H(z) -- H cancels; the budget 1/(2*pi) is
the same at every epoch. The soft sector is NEVER phase-resolved.

T2 (THE TRANSFER IDENTITY, exact): the two-mode secular weight is
<P> = (1/2) * g_c^2/(g_c^2 + delta^2) = (1/2) / (1 + (phi/gamma)^2)
with phi = delta*t, gamma = g_c*t (same t). Under the budget
phi <= eta/(2*pi) and any coupling history achieving gamma >= gamma_min,

    r >= (1/2) / (1 + eta^2/(4*pi^2*gamma_min^2)).

At (eta, gamma_min) = (1, 1): r >= 0.4875. The finite-time
(non-averaged) transfer stays within O(phi^2) of resonant as well
(numeric envelope gate G9T-3).

T3 (THE COUPLING WITNESS, argument-grade, disclosed): gamma >= ~1 is
witnessed by the anomaly's EXISTENCE: if g_c*t << 1 the dressing
transfer is O((g_c t)^2) << 1, nu -> 1, and there is no anomaly at
all. Any measured boost of order the thermal occupation certifies at
least ~one completed exchange cycle. (Not circular: it uses only
boost > 0, not the value 1/2.) This is the one non-derived premise;
gamma_min is scanned UPWARD only (stronger witnesses tighten r).

CONSEQUENCES (printed, consistency-only this stage):
  - O5-DETUNING answered-toy: resonance is not assumed; detuning is
    UNRESOLVABLE (the universe is too young to resolve the soft
    sector). The ROUND-16 A9 scenario cannot arise inside the sector.
  - Parameter-free window r in [floor, 1/2] -> the 9S data bound
    (r >= 0.315) is consistent and strictly weaker.
  - Void kill-band: p_void = 1/2 + r/2 in [~0.72, 0.75] across the
    scan (sharp cell [0.744, 0.75]) -- the P1 ceiling is predicted
    SATURATED in voids (P1 annotation upgrade at booking).
  - Fiducial-gate galaxy band: p(fid) in [1/2 + floor*G/2, 0.6884];
    the 5G point p-hat = 0.65 sits ~0.02-0.04 below the band edge =
    CONSISTENT AT CURRENT RESOLUTION (grid step 0.05, no sigma
    printed) -- named future kill: a hier tail measurement with
    sigma_p <= 0.02 demanding p < 0.67 at the fiducial gate kills
    the theorem.

PRE-REGISTERED BARS (locked before execution):
  T-THEOREM    : T1 and T2 sympy chains EXACT (zero residuals; the
                 epoch-invariance identity exact) AND the scan floor
                 min over (eta, gamma_min) in {0.5, 1, 2} x {1, pi,
                 2*pi} of r_floor is >= 0.45 AND G9T-3 envelope
                 holds (<= 0.03).
  T-CONSISTENT : chains exact but the scan floor drops below 0.45 in
                 some cell (quote the map; the theorem holds with a
                 softer window).
  T-BROKEN     : any exact chain fails or the envelope breaks.
  The letter names the claim precisely: THEOREM-WITH-WITNESSED-
  PREMISE (T3 is argument-grade and says so everywhere).

CREDENCE MAP (pre-signed): anomaly-real 53 NO MOVE. bath-mechanism
conditional (~8 since 6N; held at 9Q by the exactly-one-derived
cell): T-THEOREM AND a ROUND-17 red-team (fresh Opus agent;
authorized by the author 2026-08-07 "use more Opus agents...to
bounce your ideas") leaves NO unpatched hole in T1/T2/T3 -> 8 -> 12
(the bump 9Q withheld becomes earnable: L1+L2+T close the chain the
L3 cap was waiting for). T-THEOREM with an unpatched hole -> HOLD 8,
theorem booked with the hole named. T-BROKEN -> 8 -> 6 (a failed
rescue leaves A9 standing sharper). L3's ledger grade does NOT flip
this stage (grade flips only through a review round per the 9Q
precedent); a note is added.

DISCLOSURE (honesty block): the ENTIRE theorem sketch -- the 1/(2*pi)
budget, the r >= 0.4875 sharp cell, the epoch-invariance, the
p-hat = 0.65 vs band-edge 0.68 resolution note -- was derived
in-session BEFORE this commit. The bars bind on what was NOT
computed: the exact sympy residuals, the full (eta, gamma_min) scan
corners, the numeric envelope, and the mpmath recompute. Grades move
only down post-run; the single upgrade path is the pre-stated
ROUND-17 rule in the credence map (not a grade change of this
stage's letter).

GATES:
  G9T-0  input regression: the 6I fiducial gate and the 6E gal
         postdiction re-parsed and re-verified (p(1/2, G_fid) =
         0.6884), tying 9T to the 9S provenance chain.
  G9T-1  sympy exact: the secular-weight identity; the ratio
         monotonicities; the phase-budget identity with symbolic
         H(z) cancellation; zero residuals.
  G9T-2  mpmath 50-digit recompute of every scan cell, rel <= 1e-10.
  G9T-3  numeric envelope: max over gamma in {1, 3, 10, 30, 100},
         phi in {0, 0.08, 0.159, 0.318} of |<P>(T) - <P>_res(T)|
         (cumulative averages) <= 0.03.
  G9T-4  ledger legs live: mech-9q-onemode, mech-9q-frozen,
         mech-9q-ceiling, mech-9s-detuning all CURRENT.

Output: data/stage9t_resonance.txt. Wall-clock: seconds.
"""
import math
import os
import re

import sympy as sp
import mpmath as mp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

OUT = []
def say(s=""):
    print(s)
    OUT.append(s)

say("=" * 78)
say("STAGE 9T -- THE RESONANCE THEOREM (unresolvable detuning)")
say("=" * 78)

gates = {}

# ---------------- G9T-0: input regression (9S provenance chain) -----
t6e = open('data/stage6e_ambgate.txt', encoding='utf-8').read()
t6i = open('data/stage6i_chaegate.txt', encoding='utf-8').read()
pred_gal = float(re.search(r'G4 gal: tail p_hat = [\d.]+ \(pred ([\d.]+)\)',
                           t6e).group(1))
g_fid = float(re.search(r'fiducial s = [\d.]+ \(g = ([\d.]+)\)',
                        t6i).group(1))
ok_0 = abs((0.5 + 0.5*g_fid/2) - pred_gal) <= 5e-4
gates['G9T-0'] = ok_0
say("G9T-0 input regression: p(1/2, G_fid=%.4f) = %.4f vs 6E pred "
    "%.4f -> %s" % (g_fid, 0.5 + 0.5*g_fid/2, pred_gal,
    "PASS" if ok_0 else "FAIL"))
say("")

# ---------------- G9T-1: the sympy chains ---------------------------
gc, dl, tt, phi, gam, eta, dx = sp.symbols(
    'g_c delta t phi gamma eta dx', positive=True)
Hz = sp.Function('H', positive=True)
z = sp.symbols('z', nonnegative=True)

# (a) secular weight identity: <P> = 1/2 * gc^2/(gc^2+dl^2), and in
#     (phi, gamma) variables = 1/2 / (1 + (phi/gamma)^2)
w = sp.sqrt(gc**2 + dl**2)
Ptr = (gc**2/(gc**2 + dl**2))*sp.sin(w*tt)**2
avg = sp.simplify(sp.integrate(Ptr, (tt, 0, 2*sp.pi/w))/(2*sp.pi/w))
residA = sp.simplify(avg - gc**2/(2*(gc**2 + dl**2)))
subbed = sp.simplify(avg.subs([(dl, phi/tt), (gc, gam/tt)]))
residB = sp.simplify(subbed - sp.Rational(1, 2)/(1 + phi**2/gam**2))
# (b) monotonicities: decreasing in phi, increasing in gamma
r_expr = sp.Rational(1, 2)/(1 + phi**2/gam**2)
mono_phi = sp.simplify(sp.diff(r_expr, phi))
mono_gam = sp.simplify(sp.diff(r_expr, gam))
ok_mono = (mono_phi.is_nonpositive or
           bool(sp.simplify(-mono_phi).is_positive)) and \
          bool(sp.simplify(mono_gam).is_positive)
# (c) the phase budget with symbolic H(z): delta = dx*H/(2 pi),
#     t_max = eta/H  ->  phi = eta*dx/(2 pi), H cancels at every z
phi_budget = sp.simplify((dx*Hz(z)/(2*sp.pi))*(eta/Hz(z)))
residC = sp.simplify(phi_budget - eta*dx/(2*sp.pi))
ok_1 = (residA == 0) and (residB == 0) and ok_mono and (residC == 0)
gates['G9T-1'] = ok_1
say("G9T-1 sympy: secular-weight residual = %s; (phi,gamma)-form "
    "residual = %s; monotone (d/dphi < 0, d/dgamma > 0): %s; "
    "phase-budget H(z)-cancellation residual = %s -> %s"
    % (residA, residB, ok_mono, residC, "PASS" if ok_1 else "FAIL"))
say("  T1: phi = delta*t <= eta*dx/(2*pi); dx <= 1 in the anomaly "
    "sector (Lemma 1); EPOCH-INVARIANT (H cancels).")
say("  T2: r = (1/2)/(1 + (phi/gamma)^2) with gamma = g_c*t.")
say("  T3 (witnessed premise, argument-grade): gamma >= ~1 because "
    "the anomaly EXISTS (boost O(thermal occupation) requires at "
    "least ~one completed exchange cycle; uses boost > 0 only).")
say("")

# ---------------- the scan ------------------------------------------
say("THE FLOOR SCAN  r_floor = (1/2)/(1 + eta^2/(4 pi^2 gamma_min^2))")
say("-" * 78)
cells = []
for eta_v in (0.5, 1.0, 2.0):
    for gam_v in (1.0, math.pi, 2*math.pi):
        rf = 0.5/(1.0 + eta_v**2/(4*math.pi**2*gam_v**2))
        cells.append((eta_v, gam_v, rf))
        say("  eta = %.1f  gamma_min = %5.3f  ->  r_floor = %.4f  "
            "p_void in [%.4f, 0.7500]  phi_max = %.3f rad"
            % (eta_v, gam_v, rf, 0.5 + rf/2,
               eta_v/(2*math.pi)))
r_floor_min = min(c[2] for c in cells)
r_floor_sharp = 0.5/(1.0 + 1.0/(4*math.pi**2))
say("  scan floor min = %.4f (worst cell eta=2, gamma_min=1); sharp "
    "cell (1, 1) = %.4f" % (r_floor_min, r_floor_sharp))
say("")

# consistency rows (printed, not verdict-bearing)
say("CONSISTENCY (not verdict-bearing):")
say("  vs 9S data bound: r >= 0.315 (measured) is IMPLIED-and-weaker "
    "(theorem floor %.3f)" % r_floor_min)
p_lo = 0.5 + r_floor_min*g_fid/2
say("  fiducial-gate galaxy band: p in [%.4f, %.4f]; 5G point "
    "p-hat = 0.65 (grid 0.05, no sigma) sits %.3f below the band "
    "edge = CONSISTENT AT CURRENT RESOLUTION; named kill: a hier "
    "tail measurement with sigma_p <= 0.02 demanding p < 0.67 at "
    "the fiducial gate kills the theorem"
    % (p_lo, 0.5 + 0.5*g_fid/2, max(0.0, p_lo - 0.65)))
say("  void kill-band: p_void in [%.3f, 0.750] across the scan "
    "(sharp cell [%.4f, 0.7500]) -- the P1 ceiling is predicted "
    "SATURATED in voids" % (0.5 + r_floor_min/2,
                            0.5 + r_floor_sharp/2))
say("")

# ---------------- G9T-2: mpmath recompute ---------------------------
mp.mp.dps = 50
ok_2 = True
for (eta_v, gam_v, rf) in cells:
    rm = mp.mpf('0.5')/(1 + mp.mpf(repr(eta_v))**2 /
                        (4*mp.pi**2*mp.mpf(repr(gam_v))**2))
    ok_2 &= abs(rf/float(rm) - 1.0) <= 1e-10
gates['G9T-2'] = ok_2
say("G9T-2 mpmath 50-digit recompute of all 9 scan cells: %s"
    % ("PASS" if ok_2 else "FAIL"))

# ---------------- G9T-3: numeric envelope (finite-time) -------------
def cumavg(g_, d_, T):
    Om = math.sqrt(g_*g_ + d_*d_)
    A = g_*g_/(g_*g_ + d_*d_)
    return A*(0.5 - math.sin(2*Om*T)/(4*Om*T))

env_max = 0.0
for gam_v in (1.0, 3.0, 10.0, 30.0, 100.0):
    for phi_v in (0.0, 0.08, 1.0/(2*math.pi), 0.318):
        T = 1.0
        dev = abs(cumavg(gam_v, phi_v, T) - cumavg(gam_v, 0.0, T))
        env_max = max(env_max, dev)
ok_3 = env_max <= 0.03
gates['G9T-3'] = ok_3
say("G9T-3 finite-time envelope: max |<P> - <P>_res| over the "
    "(gamma, phi) grid = %.4f (bar 0.03) -> %s"
    % (env_max, "PASS" if ok_3 else "FAIL"))

# ---------------- G9T-4: ledger legs --------------------------------
need = ['mech-9q-onemode', 'mech-9q-frozen', 'mech-9q-ceiling',
        'mech-9s-detuning']
led = {}
for line in open('LEDGER.csv', encoding='utf-8'):
    for rid in need:
        if line.startswith(rid + ','):
            led[rid] = line.split(',')[1]
ok_4 = all(led.get(rid) == 'CURRENT' for rid in need)
gates['G9T-4'] = ok_4
say("G9T-4 ledger legs: " + "; ".join("%s=%s" % (rid,
    led.get(rid, 'MISSING')) for rid in need) +
    " -> %s" % ("PASS" if ok_4 else "FAIL"))
say("")

# ---------------- verdict -------------------------------------------
allok = all(gates.values())
if allok and r_floor_min >= 0.45:
    verdict = ("T-THEOREM (with witnessed premise): detuning inside "
               "the soft sector is UNRESOLVABLE -- phi <= eta/(2*pi) "
               "at every epoch -- so r in [%.3f, 0.500] "
               "parameter-free; O5-DETUNING answered at toy grade "
               "pending ROUND-17" % r_floor_min)
elif allok:
    verdict = ("T-CONSISTENT: chains exact; scan floor %.3f < 0.45 "
               "in some cell" % r_floor_min)
else:
    verdict = "T-BROKEN: a chain or gate failed; no claim"
say("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
    for k, v in sorted(gates.items())) +
    "  -> ALL %s" % ("PASS" if allok else "FAIL"))
say("")
say("VERDICT (locked grammar): " + verdict)
say("")
say("CREDENCE MAP (pre-signed): anomaly-real 53 NO MOVE. "
    "bath-mechanism conditional: T-THEOREM + ROUND-17 "
    "no-unpatched-hole -> 8 -> 12; T-THEOREM with a hole -> HOLD 8 "
    "(hole named); T-BROKEN -> 6. Final number booked in NOTES "
    "after ROUND-17.")

with open('data/stage9t_resonance.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(OUT) + "\n")
say("")
say("written: data/stage9t_resonance.txt")
