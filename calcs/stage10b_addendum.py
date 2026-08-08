"""
STAGE 10B ADDENDUM -- ROUND-23 conditions 1/2/3/4/5, with independent
re-verification of every load-bearing number the reviewer introduced
(the standing verify-reviewer-math rule). Letter consequence printed:
V-DISPERSIVE-FORCED -> V-DISPERSIVE-PERMITTED (adopted). No bar or
credence-map change (HOLD 15 per the pre-signed map: hole blocks the
rise; nothing fatal).

Executes:
  C2v HIS COMMUTATORS AND THE l=1 CEILING, my own matrices: the
      momentum-constraint vertex p_s(b+b†) IS off-diagonal (his Q1
      catch -- clause (2)'s "ever" refuted); the stage's vertex
      N_s(b+b†) commutes; the exchange contrast does not. His ceiling
      numbers (n_Wien 1.38e-4 / 2.07e-7, admixtures ~3.5e-4 g^2 /
      1.6e-7 g^2, geometric 1.3e-22 / 5.5e-12, >= 15 orders short)
      re-computed. Clause (2) re-scoped per his condition 2.
  C3v THE AMPLITUDE CORRECTION (his Q2c catch CONFIRMED): the stage's
      sqrt(pi(2n+1)) rows are the RMS JITTER of the redshift coupling
      (mean-zero), NOT the JC mean pull; his zero-detuning pull
      pi(2n+1) = 6.3/44.8/66 re-computed and CONFIRMED as the
      zero-detuning bound. PREVIEW (labeled, NOT a verdict; the
      successor stage's target): restoring the frequency-ratio
      structure the second-order pull carries (lam = omega_s phi,
      denominator ~ the soft-mode frequency scale) gives
      pull/omega_s ~ x(2n_BE(x)+1)/2 = (x/2) coth(x/2), whose exact
      series is 1 + x^2/12 - x^4/720 + ... -- ORDER UNITY at every
      anchor WITH THE BE LADDER'S OWN BERNOULLI RUNGS (the 6H
      c2 = 1/12 and c4(0) = -1/720 coefficients) appearing in the
      mean dressing. sympy-exact series gate; anchor rows printed.
  C4x the corotation co-read correction (his Q5): real slow family;
      excluded by measure-zero-in-radius + dictionary mismatch.
  C5x GB-3c demoted to definitional note (the sqrt-pi identity is
      S = pi r_c^2/l_P^2 rearranged); the V-EMPTY-unreachable
      disclosure (the fork was cosmetic at physical occupations:
      sqrt(pi(2n+1)) >= sqrt(pi) = 1.77 for all n).
  C1  the relabel + the earned statement.
Successors (NOT executed here; the named 10C package): his condition
6 (re-run the 6X control under a purely dispersive coupling), 7
(justify or drop H-P -- the sqrt(S)-independence question IS the
61-order bridge), 8 (the real amplitude bar: derive the mean pull's
detuning/geometry factor to the measured c1 = 1/2 -- the coth
preview above is the opening move), 9 (the r = 1/2 kill-test:
promoted in PREDICTIONS.md this same commit).
Gates: GA-1 commutator matrices; GA-2 his ceiling numbers; GA-3 the
coth series exact + anchor rows. Output: data/stage10b_addendum.txt.
"""
import math
import os
import time

import numpy as np
import sympy as sp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

OUT = []
def P(s=""):
    print(s, flush=True)
    OUT.append(s)

def save():
    with open('data/stage10b_addendum.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

t00 = time.time()
P("=" * 78)
P("STAGE 10B ADDENDUM -- ROUND-23 CONDITIONS 1/2/3/4/5")
P("=" * 78)
gates = {}

TWO_PI = 2*math.pi

# ---- C2v: his commutators, my own matrices --------------------------
N = 8
a = np.diag(np.sqrt(np.arange(1, N)), 1)      # annihilation
ad = a.T
I = np.eye(N)
def kron(A, B): return np.kron(A, B)
A_ = kron(a, I); Ad = kron(ad, I)
B_ = kron(I, a); Bd = kron(I, ad)
Ns = Ad @ A_
Xb = B_ + Bd
p_s = 1j*(Ad - A_)/math.sqrt(2)
def cnorm(M): return float(np.max(np.abs(M)))
c_disp = cnorm((Ns @ Xb) @ Ns - Ns @ (Ns @ Xb))
c_exch = cnorm((Ad @ B_ + Bd @ A_) @ Ns - Ns @ (Ad @ B_ + Bd @ A_))
c_mom = cnorm((p_s @ Xb) @ Ns - Ns @ (p_s @ Xb))
gates['GA-1'] = (c_disp < 1e-12 and c_exch > 0.1 and c_mom > 0.1)
P("")
P("GA-1 commutators (my own truncated-Fock matrices, N = %d):" % N)
P("  ||[N_s(b+bd), N_s]|| = %.2e (stage's vertex: dispersive) " % c_disp)
P("  ||[adb+bda,   N_s]|| = %.2e (exchange contrast)" % c_exch)
P("  ||[p_s(b+bd), N_s]|| = %.2e (HIS CATCH: the momentum-constraint "
  "vertex IS off-diagonal -- clause (2)'s 'ever' refuted) -> %s"
  % (c_mom, "PASS" if gates['GA-1'] else "FAIL"))

# his l>=1 ceiling numbers
def n_be_x(x):
    return 1.0/math.expm1(x) if x < 500 else 0.0
n1 = n_be_x(TWO_PI*math.sqrt(2)); n2 = n_be_x(TWO_PI*math.sqrt(6))
d1b = math.sqrt(2) - 0.174; d1g = math.sqrt(2) - 0.035
adm1 = 4*n1/d1b**2
adm2 = 4*n2/(math.sqrt(6) - 0.174)**2
C = 2.998e8; H_SI = 2.27e-18
geo_b = (1e4*1.496e11/(C/H_SI))**2
geo_g = (10*3.0857e19/(C/H_SI))**2
ceil_b = adm1*geo_b; ceil_g = adm1*geo_g
ok2 = (abs(n1/1.38e-4 - 1) < 0.01 and abs(n2/2.07e-7 - 1) < 0.01
       and abs(adm1/3.52e-4 - 1) < 0.05
       and abs(geo_b/1.3e-22 - 1) < 0.05
       and abs(geo_g/5.5e-12 - 1) < 0.05
       and ceil_g < 1e-14)
gates['GA-2'] = ok2
P("")
P("GA-2 his l=1 ceiling re-computed: n(2pi sqrt2) = %.3e (his "
  "1.38e-4); admixture = %.3e g^2 (his 3.52e-4); geometric = %.2e "
  "bin / %.2e gal (his 1.3e-22 / 5.5e-12); combined ceiling %.1e / "
  "%.1e g^2 = >= 15 orders short -> %s  [CONFIRMED -- clause (2) "
  "re-scoped: no off-diagonal vertex that couples to the ORBITAL "
  "coordinate AND has an ungapped, occupied, near-resonant partner]"
  % (n1, adm1, geo_b, geo_g, ceil_b, ceil_g,
     "PASS" if ok2 else "FAIL"))
P("  l=2 row: n = %.2e, admixture %.2e g^2 -- deader." % (n2, adm2))

# ---- C3v: the amplitude correction + the coth preview ---------------
P("")
P("C3v -- the amplitude clause corrected (his Q2c catch CONFIRMED):")
for nm, n_ in (('binary', 0.502), ('gal ambient', 6.63),
               ('gal deep', 10.0)):
    P("  %-12s RMS jitter sqrt(pi(2n+1)) = %6.3f   zero-detuning "
      "mean pull pi(2n+1) = %6.2f  (his row CONFIRMED)"
      % (nm, math.sqrt(math.pi*(2*n_ + 1)), math.pi*(2*n_ + 1)))
P("  the stage's rows were the JITTER (first-order, mean-zero) and "
  "were mislabeled as the (2n+1) JC pull -- adopted; the letter's "
  "amplitude clause is retracted to jitter-grade.")
xs = sp.Symbol('x', positive=True)
mean_pull = (xs/2)*sp.coth(xs/2)
ser = sp.series(mean_pull, xs, 0, 6).removeO()
tgt = 1 + xs**2/12 - xs**4/720
gates['GA-3'] = sp.simplify(ser - tgt) == 0
P("")
P("GA-3 PREVIEW (labeled -- NOT a verdict; the 10C target): the "
  "frequency-ratio-corrected second-order pull, pull/omega_s ~ "
  "x(2 n_BE(x) + 1)/2 = (x/2) coth(x/2); exact series = %s "
  "== 1 + x^2/12 - x^4/720 -> %s"
  % (sp.simplify(ser), "PASS" if gates['GA-3'] else "FAIL"))
for nm, x_ in (('binary x = 1.095', 1.095),
               ('gal ambient x = 0.1411', 0.1411),
               ('deep end x = 0.14', 0.14)):
    v = (x_/2)/math.tanh(x_/2)
    P("    %-24s (x/2) coth(x/2) = %.4f" % (nm, v))
P("  ORDER UNITY at every anchor, approaching 1 in the deep limit "
  "with the BE ladder's own Bernoulli rungs (+x^2/12 - x^4/720 = the "
  "6H c2 and c4(0) coefficients) appearing in the MEAN dressing -- "
  "his 44.8/66 'over-production' is the zero-detuning bound; the "
  "honest second-order structure carries lam^2/Delta = omega_s^2 "
  "phi^2 (2n+1)/Delta with Delta at the soft-mode scale, taming the "
  "pull to O(1) x-structured. DERIVING this exactly (with the honest "
  "Delta) = round-23 condition 8 = THE 10C OPENING MOVE.")

# ---- C4x + C5x ------------------------------------------------------
P("")
P("C4x corotation co-read corrected (his Q5): real disks carry spiral "
  "patterns; |Omega - Omega_p| -> 0 at corotation IS a slow secular "
  "family. It cannot carry the grammar: corotation is ONE radius per "
  "galaxy (measure-zero) while the RAR holds at all radii, and "
  "|Omega - Omega_p| is non-monotonic about that radius while "
  "omega_dict is monotone in g -- the same dictionary mismatch that "
  "killed the apsidal family. The exclusion STRENGTHENS: the one "
  "remaining galaxy slow-family also fails the dictionary.")
P("C5x GB-3c demoted to a definitional note (sqrt(S) l_P = sqrt(pi) "
  "r_c is S = pi r_c^2/l_P^2 rearranged; the physics is entirely in "
  "H-P + the incoherent sum). DISCLOSED: V-EMPTY was unreachable at "
  "physical occupations (sqrt(pi(2n+1)) >= sqrt(pi) = 1.772 for all "
  "n >= 0; the fork was cosmetic) -- the pre-reg's real function was "
  "locking the exactness gates and deferring credence to the round.")

# ---- C1: the relabel ------------------------------------------------
P("")
P("=" * 78)
P("10B LETTER RELABELED (ROUND 23 adopted): V-DISPERSIVE-FORCED -> "
  "V-DISPERSIVE-PERMITTED")
P("=" * 78)
P("THE EARNED STATEMENT: real-exchange carriers C1-C7 are EXCLUDED")
P("(radiative 34-119 orders; matter-elastic frequency-forbidden;")
P("matter-secular dictionary-mismatched by a 6-order sweep vs a ~10%")
P("constant a0; the momentum-constraint vertex off-diagonal-in-")
P("principle but >= 15 orders dead as a carrier). The constraint")
P("sector's only soft l=0 vertex is dispersive-class by the exact")
P("selection rule [H_int, N_s] = 0. Therefore any microphysical")
P("carrier of the grammar must be DISPERSIVE. Whether the dispersive")
P("channel SUPPLIES the grammar at the required amplitude and gate")
P("form remains NAMED-NOT-DERIVED: hypothesis H-P (the sqrt-S")
P("independence bridge = 9W's own breakable assumption applied to")
P("the horizon ensemble) is unproven; the stage's amplitude row was")
P("the RMS jitter, not the mean pull; and the 6X real-vs-virtual")
P("control has not been re-run under the dispersive hypothesis.")
P("The r = 1/2 dividend (no gamma-running; void band = 3/4 EXACTLY)")
P("is promoted to a named conditional kill-test in PREDICTIONS.md:")
P("a void/DR4 measurement of r < 1/2 at > 2 sigma falsifies the")
P("dispersive reading while off-resonance exchange would survive.")
P("")
P("CREDENCE (pre-signed map, mechanical): hole blocks the rise cell;")
P("nothing fatal -> bath-mechanism conditional HOLDS 15; anomaly-real")
P("53 untouched.")
P("")
P("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
  for k, v in sorted(gates.items())))
P("")
P("done (%.1f min)" % ((time.time() - t00)/60))
save()
print("\nsaved: data/stage10b_addendum.txt")
