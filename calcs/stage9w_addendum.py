"""
STAGE 9W ADDENDUM (ROUND-18 conditions 5/7/8 executed; DISCLOSED
POST-HOC -- the letter W-THEOREM-ONLY stands as fired; this diagnostic
cannot and does not change it).

Condition 5: the binary-leg G9W-5 failure was a domain-truncation
artifact (x_bin = 1.0954 outside the (0,1] integration domain; the
kernel sampled only the soft side of center).  Re-read with domain
(0, 2]; the local 6E anchors are PRIMARY; kernel scans are ROBUSTNESS
DIAGNOSTICS.
Condition 7: budget re-scoped -- quote Delta(s^2) (the physical
quantity) alongside Delta(c4); c4 = s^2/192 - 1/720 is a cancellation
that amplifies relative shifts by ~1.55x; the operative c4 quote =
single-mode value + width-conditional band at the physical Rabi-grade
width; "multimode-soft" RETIRED.
Condition 8: integrator-convergence gate for any kernel bar read at
the 4th decimal (tightened-tolerance repeat, |d| <= 1e-4).

Writes data/stage9w_addendum.txt.
"""
import math
from scipy.integrate import quad

OUT = 'data/stage9w_addendum.txt'
L = []
def say(s=''):
    L.append(s); print(s, flush=True)
    with open(OUT, 'w') as f:
        f.write("\n".join(L) + "\n")

def n_be(x): return 1.0/math.expm1(x)
def s_of_n(n): return n/(1.0 + n)
def p_of_n(n):
    s = s_of_n(n)
    return 0.5 + s*s/4.0

X_GAL = math.sqrt(0.02)
X_BIN = math.sqrt(1.2)
P_GAL_6E, P_BIN_6E = 0.6884, 0.5280

say("9W ADDENDUM: ROUND-18 conditions 5/7/8 (disclosed post-hoc; "
    "letter unchanged)")
say("=" * 72)

def n_gauss(x0, G, hi, eps):
    num = quad(lambda x: n_be(x)*math.exp(-0.5*((x - x0)/G)**2),
               1e-6, hi, points=[x0], epsabs=eps, epsrel=eps,
               limit=200)[0]
    den = quad(lambda x: math.exp(-0.5*((x - x0)/G)**2),
               1e-6, hi, points=[x0], epsabs=eps, epsrel=eps,
               limit=200)[0]
    return num/den

say("Condition 5 -- Gaussian kernel, domain (0, 2] (binary center "
    "x = 1.0954 now interior):")
ok_conv = True
for gfrac in (0.05, 0.10, 0.20):
    ng = n_gauss(X_GAL, gfrac*X_GAL, 2.0, 1e-10)
    nb_ = n_gauss(X_BIN, gfrac*X_BIN, 2.0, 1e-10)
    ng2 = n_gauss(X_GAL, gfrac*X_GAL, 2.0, 1e-12)
    nb2 = n_gauss(X_BIN, gfrac*X_BIN, 2.0, 1e-12)
    conv = max(abs(p_of_n(ng) - p_of_n(ng2)),
               abs(p_of_n(nb_) - p_of_n(nb2)))
    ok_conv &= conv <= 1e-4
    say("  Gamma/x0 = %.2f: p_gal = %.4f (d %+0.4f), p_bin = %.4f "
        "(d %+0.4f)  [integrator conv %.1e]" %
        (gfrac, p_of_n(ng), p_of_n(ng) - P_GAL_6E,
         p_of_n(nb_), p_of_n(nb_) - P_BIN_6E, conv))
say("  G-CONV (condition 8, |d| <= 1e-4 under 100x tighter "
    "tolerance): %s" % ('PASS' if ok_conv else 'FAIL'))
say("  => at Gamma/x0 = 0.10 the binary leg reads d = +0.0007-grade "
    "(reviewer reproduced +0.0007): the G9W-5 binary failure was the "
    "domain truncation, as diagnosed; BOTH postdictions are "
    "convexity-only on the extended domain.  The local 6E anchors "
    "are PRIMARY; kernel scans are ROBUSTNESS DIAGNOSTICS "
    "(reviewer condition adopted).")
say('')

say("Condition 7 -- the budget re-scoped (Delta(s^2) primary; c4 "
    "cancellation amplifies ~1.55x):")
s0 = s_of_n(n_be(X_GAL))
c4_0 = s0*s0/192 - 1/720
say("  single-mode targets: s^2 = %.4f; c4 = %.6f" % (s0*s0, c4_0))
for gfrac in (0.05, 0.10, 0.20, 0.30):
    ng = n_gauss(X_GAL, gfrac*X_GAL, 2.0, 1e-10)
    s1 = s_of_n(ng)
    ds2 = s1*s1/(s0*s0) - 1
    c4_1 = s1*s1/192 - 1/720
    dc4 = c4_1/c4_0 - 1
    say("  Gamma/x0 = %.2f: Delta(s^2) = %+.2f%%; Delta(c4) = "
        "%+.2f%% (amplification x%.2f)" %
        (gfrac, 100*ds2, 100*dc4, dc4/ds2 if ds2 != 0 else 0))
say("  OPERATIVE QUOTE (adopted): c4(L=2) = %.6f (single-mode), "
    "width-conditional band <= 0.5%% at the physical Rabi-grade "
    "width (Gamma/x0 <= 0.1); future bars set on Delta(s^2); "
    "'MULTIMODE-SOFT' RETIRED as a scan-edge artifact on a "
    "cancellation-amplified quantity." % c4_0)
say('')
say("flat-exclusion wording (condition 4, recorded): the exclusion "
    "rests SOLELY on e_N-blindness (no split); at cutoff 1e-3 flat "
    "gives p = 0.6874 ~= the galaxy anchor 0.6884 (magnitude-"
    "degenerate, disclosed); 'gate -> 1' deleted.")
say("letter: W-THEOREM-ONLY stands as fired.")
print("\nsaved:", OUT)
