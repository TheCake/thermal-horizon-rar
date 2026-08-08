"""
9X/9Y ROUND-19 ADOPTION ADDENDUM (disclosed post-hoc; executes the
reviewer's condition 3 (the decisive soft-regime linearity test) and
condition 2 (the tail-VIRTUAL retraction arithmetic); records every
relabel.  The as-fired letters stand in their own outputs; the
OPERATIVE letters after ROUND-19 are printed here.)

Reviewer conditions (adopted verbatim, REVIEW-ROUND19-OPUS.md):
 1  X-CLOSED downgraded (soft leak not closed by 9X)
 2  "tail VIRTUAL" RETRACTED (pure-lending Lorentzian null = 3.8-3.9
    at 16 lam >= my bar 3; the lambda-ratio diagnostic cannot separate)
 3  soft-regime response: P2(n2) at large detuning, NB >= 5-6 n2 --
    pre-registered interpretation: LINEAR => the leak lives
 4  the single missing scale g_c^sky/H + explicit rho(omega) =
    the successor computation (O5-LEAK, named)
 5  cumulative dispersive shift not bounded by 9T -- separate budget
 6  locality = data-supported (6Y), not mechanism-derived
 7  9Y renamed NON-THERMALITY METER (classical super-thermal mixtures
    gate-degenerate with squeezing on the sky scalar; parity lab-only)
 8  gate functional off-thermal ambiguous (P(n>=2) vs [P(n>=1)]^2
    split ~50% under squeezing; <a^2> pair channel unmodeled) ->
    Y-THEOREM relabeled Y-DIAL (conditional); DR4 lever conditional
 9  frozen-window consistency condition stated for 9Y

Writes data/stage9xy_addendum.txt.
"""
import math
import numpy as np

OUT = 'data/stage9xy_addendum.txt'
L = []
def say(s=''):
    L.append(s); print(s, flush=True)
    with open(OUT, 'w') as f:
        f.write("\n".join(L) + "\n")

say("9X/9Y ROUND-19 ADOPTION ADDENDUM")
say("=" * 72)

def low(d):
    return np.diag(np.sqrt(np.arange(1, d)), 1)

def thermal_vec(N, nbar):
    x = math.log(1.0 + 1.0/nbar)
    p = np.exp(-x*np.arange(N)); p /= p.sum()
    return p

CHI, LAM, WA = 0.8, 0.02, 5.0

def dressward1(NB, n, delta):
    a4 = low(4)
    b = low(NB)
    A = np.kron(a4, np.eye(NB))
    B = np.kron(np.eye(4), b)
    Na = A.T @ A
    H = WA*Na + 0.5*CHI*(Na @ (Na - np.eye(4*NB))) \
        + (WA + CHI + delta)*(B.T @ B) + LAM*(A.T @ B + B.T @ A)
    ev, U = np.linalg.eigh(H)
    p_a = np.zeros(4); p_a[1] = 1.0
    rho = np.diag(np.kron(p_a, thermal_vec(NB, n)))
    rl = U.T @ rho @ U
    rb = U @ np.diag(np.diag(rl)) @ U.T
    proj = np.kron(np.diag((np.arange(4) == 2)*1.0), np.eye(NB))
    return float(np.real(np.trace(proj @ rb)))

say("CONDITION 3 -- the soft-regime linearity test (pre-registered "
    "interpretation: LINEAR => the leak lives):")
for delta in (2.5, 5.0):
    say("  delta = %.1f (= %.0f lam):" % (delta, delta/LAM))
    rows = []
    for n2 in (2, 8, 20, 50, 100):
        NB = max(24, 6*n2)
        trunc = (n2/(1.0 + n2))**NB
        P2 = dressward1(NB, float(n2), delta)
        pred = 4*LAM*LAM*n2/(delta*delta)
        rows.append((n2, P2))
        say("    n2 = %4d (NB %4d, trunc %.0e): P2 = %.6f | "
            "P2/n2 = %.2e | 4 lam^2 n/delta^2 = %.6f" %
            (n2, NB, trunc, P2, P2/n2, pred))
    slope = (math.log(rows[-1][1]) - math.log(rows[0][1])) / \
            (math.log(rows[-1][0]) - math.log(rows[0][0]))
    say("    log-log slope n2 = 2 -> 100: %.3f (1.0 = linear)" % slope)
say("  VERDICT (condition 3): the response is LINEAR in occupation at "
    "soft detunings across n2 = 2-100 -- the reviewer's numbers "
    "reproduce; per the pre-registered interpretation THE LEAK LIVES "
    "as an open question: the 9X S-ratio was measured at the one "
    "detuning (10 lam) where saturation onsets early (k* = 12.5).")
say('')

say("CONDITION 2 -- the tail-VIRTUAL retraction (arithmetic):")
for n, gc in ((0.5, 0.0656), (2.0, 0.0822)):
    d16 = 16*LAM
    R_null = (1 + (2*d16/gc)**2)/(1 + (d16/gc)**2)
    say("  n=%.1f: pure-lending Lorentzian null at 16 lam (width "
        "halves with lambda): R = %.2f" % (n, R_null))
say("  measured R = 3.37 < null 3.8-3.9: the lambda-ratio diagnostic "
    "CANNOT separate lending from virtual in a Lorentzian tail; "
    "'tail VIRTUAL' is RETRACTED.  (The 9X core fit, the S-ratio "
    "values, and the kernel rows themselves all stand -- the "
    "retraction is of one interpretive claim.)")
say('')

say("OPERATIVE LETTERS AFTER ROUND-19:")
say("  9X RELABELED: X-OPEN (was X-CLOSED as fired).  The soft-sector "
    "leak is NOT closed: the decisive object is the continuum "
    "integral int rho(omega) n(omega) (g_c/delta)^2 domega plus the "
    "single missing scale g_c^sky/H -- both outside a discrete "
    "NB<=56 toy.  Locality stays DATA-SUPPORTED (6Y), mechanism-"
    "consistent-in-tested-regime, NOT mechanism-derived.  The l-gap "
    "single-channel corollary is PENDING the same scale.  The "
    "cumulative dispersive shift (sum rho n lam^2/delta) needs its "
    "own budget -- 9T does not bound it.")
say("  9Y RELABELED: Y-DIAL (conditional) (was Y-THEOREM as fired). "
    "The instrument is a NON-THERMALITY METER: classical super-"
    "thermal mixtures (e.g. mixture-of-thermals, reviewer demo "
    "+0.014 at the binary anchor) move the gate the same direction "
    "as squeezing and are degenerate on the sky's scalar; parity "
    "would separate them but is lab-only.  The dial magnitude is "
    "conditional on the g = P(n>=2) reading ([P(n>=1)]^2 splits "
    "~50% under squeezing; the <a^2> coherent-pair channel is "
    "unmodeled).  The DR4 lever inherits the conditional.  The "
    "frozen-window consistency condition (squeezing persistence vs "
    "averaging convergence vs 6N dephasing) is REQUIRED and "
    "unestablished.")
say("  SURVIVES REVIEW UNCONDITIONALLY: the 9Y Gaussian-classical "
    "absorbable-vs-squeezing-anisotropy covariance mathematics "
    "(exact); the P1 void state-independence line (reviewer "
    "reproduced 0.9678 vs 0.9675); the 9X kernel rows, core fits, "
    "and the thermal-mixture-of-Rabi-Lorentzians width mechanism "
    "(his model reproduces our rows to <1% in the core).")
say('')
say("THE NAMED SUCCESSOR (O5-LEAK): (i) an explicit horizon-side "
    "mode density rho(omega) for the soft sector; (ii) the scale "
    "g_c^sky/H; (iii) the continuum lending-leak integral and the "
    "dispersive-shift budget against P2_sys.  The l-gap corollary "
    "and the soft-leak closure BOTH collapse to (i)+(ii).")
say('')
say("CREDENCE (pre-signed package map, executed on the ruling): "
    "any-hole -> HOLD bath-mechanism conditional 15.  "
    "anomaly-real 53 UNTOUCHED.")
print("\nsaved:", OUT)
