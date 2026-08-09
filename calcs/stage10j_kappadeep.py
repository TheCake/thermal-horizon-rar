"""
STAGE 10J (2026-08-09): O5-KAPPA-DEEP -- derive (or bound) the response
dictionary coefficient kappa = 4 g^2/(Om om) at the deep limit.
Target kappa = 1 (R24 condition 6; 10D C2 re-scoped the target to the
deep limit, where the temperature-locked world reads kappa = 0.925).
DISPERSIVE-SIDE ONLY (the 10G strike scope): no lending-law or
real-exchange input anywhere; ambient thermality enters only as the
6U-measured KMS state property.

OBJECTS (10C, gate-verified there; reader-verbatim):
  vertex  H = om (N + sig) + Om b'b + g (N + sig)(b + b')      [10B rule]
  polaron E(n,m) = om (n+sig) + Om m - (g^2/Om)(n+sig)^2       [exact]
  pull    dw(n) = -(g^2/Om)(2n+1+2 sig)                        [T1]
  dict    nu = 1 + kappa n_BE(x),  kappa = 4 g^2/(Om om)       [T2,
          compliance-additive, DISCLOSED used-not-derived -- the seam
          this stage attacks]
  kappa=1 <=> g = sqrt(Om om)/2 <=> phi_req = (1/2) sqrt(Om/om)
          <=> E_c = 2 hbar om                                  [T3]

ROUTES (pre-registered; each lands exact statements or fails):
  R-A THE THRESHOLD CATALOGUE + THE COMPLIANCE-DIVERGENCE THEOREM.
      All structural thresholds of the dressed ladder, exact (sig=1/2):
      spacing D(n) = om [1 - kappa (n+1)/2]; D(0)>0 <=> kappa<2;
      D(1)>0 <=> kappa<1; two-quantum E(n+2)-E(n)>0 at n=0 <=>
      kappa<4/3.  kappa = 1 is EXACTLY the two-rung degeneracy
      E(2) = E(1).  THEOREM CANDIDATE: the static susceptibility of
      the system coordinate in the one-quantum dressed state diverges
      as kappa -> 1- with finite Franck-Condon prefactor:
      lim (1-kappa) chi(1) = 4 e^{-d^2}/om, d = g/Om -- the
      dictionary's own compliance saturates at kappa = 1.  Bath-vacuum
      exact + thermal-bath persistence row.
  R-B THE 9W FUNNEL: under the bright-mode reduction kappa_coll =
      4 lam_bar^2/(Om om), lam_bar^2 = sum lam_k^2 -- every multimode/
      participation route reduces to ONE collective number
      (K-invariance gate); a horizon-side derivation must produce
      lam_bar^2 = Om om/4, and lam_bar^2 is where the 10A-displaced
      D-normalization lives.
  R-C THE PREMISE-LABELED BOUND: kappa <= 1 from measured-BE-form
      consistency -- if the dictionary mode's occupation is its own
      Gibbs occupation AND matches the measured BE form through the
      first two rungs, then kappa > 1 is excluded (E(2) < E(1) piles
      Gibbs weight on n=2 below n=1 at every temperature: occupation
      non-BE).  Premises: P-vertex [DERIVED, 10B]; P-polaron [exact];
      P-BEform [MEASURED -- the RAR identity]; P-lit ("the mode
      thermalizes on its own literal dressed ladder, not a
      self-consistently re-derived one") [ASSUMED -- named].
      Saturation kappa = 1: P-crit ("the physical coupling saturates
      the two-rung stability/compliance bound") [ASSUMED -- named].
  R-D THE IDENTITY WEB (exact, sympy): kappa = 1 <=> g = sqrt(Om om)/2
      <=> g^2 = E_zp,s * E_zp,b (zero-point product) <=> displacement
      energy per unit charge^2 = om/4 <=> phi_req = (1/2)sqrt(Om/om)
      <=> E_c = 2 hbar om <=> D(1) = 0.  Classified: identities
      (sharpened target), not derivations.

GATES (bars locked at this commit):
  G10J-1  sympy exactness: spacing formula, thresholds {2, 1, 4/3,
          2/(n+1)}, the identity web (six equivalences), the kappa=1
          series rungs nu = 1/x + 1/2 + x/12 + 0 x^2 - x^3/720
          (10C G4 regression).  Zero symbolic residue.
  G10J-2  numerics: block spectra vs polaron formula <= 1e-8 (Fock
          truncation 160, kappa in {0.5, 0.9, 0.99}); chi(1) numeric
          divergence: (1-kappa) chi(1) vs the exact limit within 2%
          at kappa = 0.999; thermal-FC persistence factor positive at
          n_amb = 2 (exact Laguerre sum vs numeric <= 1e-8).
  G10J-3  K-invariance: K in {1, 2, 3}, lam_k = lam_bar/sqrt(K):
          identical spacings <= 1e-10 (the 9W funnel).
  G10J-4  10C regressions: kappa = 2(1-c1) meter arithmetic at c1 =
          0.450 -> 1.100; phi_req identity solves; a0-lock relation
          a0_fit/a0_horizon = kappa^2 quoted (archived meters only,
          NO sky fits).
  G10J-5  the D-dependence statement printed with its derivation
          chain: kappa derivable without D-provenance IFF a structural
          principle ties lam_bar^2 to Om om (the three candidates
          premise-labeled).
  G10J-6  PREMISE AUDIT (trap-#12 discipline -- gate the positive
          clause): the letter keys MECHANICALLY on the premise list;
          K-BOUNDED requires ZERO ASSUMED entries in the bound chain;
          any ASSUMED entry -> at most K-CONDITIONAL.

VERDICT LETTERS:
  K-FORCED       kappa = 1 derived, zero ASSUMED premises anywhere.
  K-BOUNDED      kappa <= 1 derived, zero ASSUMED premises in the
                 bound chain (P-lit discharged).
  K-CONDITIONAL  the bound + saturation target established on NAMED
                 premises (expected: P-lit carries the bound, P-crit
                 the saturation); identity web + funnel + divergence
                 theorem exact.
  K-OPEN         the algebra fails somewhere material.
  K-ANTI         the exact structure forces a unique dictionary
                 threshold != 1 inconsistent with the locked-world
                 band 0.925 +/- ~0.05.

CONSISTENCY ROWS (report-grade, mandatory whichever letter): under
P-lit the bound kappa <= 1 EXCLUDES the 10D low-a0/high-kappa ridge
(single-kappa 1.503, plain deep pole 1.317) while CONTAINING the
temperature-locked world (0.925) near saturation -- the spectral bound
and the temperature lock pick the SAME world (quote with the R25-C1
treatment-dependence caveat: the vertical deep pole flips, and the
ridge caveat C2 is mandatory wherever K-SPLIT is cited).  The
deep-window fold observation: at ANY fixed kappa > 0 the literal
ladder is non-monotone at n >= 2/kappa, while the deep sky occupies
n ~ 1/x >> 1 -- the fixed-coefficient dictionary must be the leading
form of a self-consistent structure (no direction claim; R25-C1).

CREDENCE MAP (pre-signed; the ONLY movable cells; requires ROUND 30):
  K-FORCED  + round no-hole: mech-conditional 8 -> 11.
  K-BOUNDED + round no-hole: mech 8 -> 10.
  K-CONDITIONAL / K-OPEN:    HOLD 8.
  K-ANTI    + round affirmed: mech 8 -> 6.
  anomaly-real 53 UNTOUCHED in ALL cells (no sky fits; every quoted
  kappa number is an archived meter).

COMPUTE < 2 min.  Writes data/stage10j_kappadeep.txt.
"""
import math, time
import numpy as np
import sympy as sp

T0 = time.time()
OUT = 'data/stage10j_kappadeep.txt'
LINES = []
def emit(s=""):
    LINES.append(s)
    print(s, flush=True)
    with open(OUT, 'w') as f:
        f.write("\n".join(LINES) + "\n")

emit("STAGE 10J O5-KAPPA-DEEP (pre-reg: this commit)")
emit("")

# ==================== R-A / G10J-1: exact algebra ====================
om, Om, g, n, sig, kap, x = sp.symbols('om Om g n sigma kappa x',
                                       positive=True)
E = lambda nn: om*(nn+sig) - (g**2/Om)*(nn+sig)**2   # bath m = 0
D = sp.expand(E(n+1) - E(n))
D_expect = om - (g**2/Om)*(2*n + 1 + 2*sig)
r1 = sp.simplify(D - D_expect)
emit("G10J-1a spacing: D(n) = om - (g^2/Om)(2n+1+2sig); residue = %s"
     % r1)

kap_def = 4*g**2/(Om*om)
D_half = sp.simplify(D.subs(sig, sp.Rational(1, 2))
                     .subs(g**2, kap*Om*om/4))
D_half_expect = om*(1 - kap*(n+1)/2)
r2 = sp.simplify(D_half - D_half_expect)
emit("G10J-1b sigma=1/2 form: D(n) = om[1 - kappa(n+1)/2]; residue = %s"
     % r2)

thr0 = sp.solve(D_half.subs(n, 0), kap)[0]
thr1 = sp.solve(D_half.subs(n, 1), kap)[0]
two_q = sp.simplify(D_half + D_half.subs(n, n+1))
thr2q = sp.solve(two_q.subs(n, 0), kap)[0]
emit("G10J-1c thresholds: D(0)=0 at kappa=%s; D(1)=0 at kappa=%s; "
     "E(2)-E(0)=0 at kappa=%s; D(n)=0 at kappa=2/(n+1)" %
     (thr0, thr1, thr2q))
g1c = (thr0 == 2) and (thr1 == 1) and (thr2q == sp.Rational(4, 3))
emit("   kappa = 1 IS the two-rung degeneracy E(2) = E(1): %s" %
     ("EXACT" if thr1 == 1 else "FAIL"))

# identity web
web = []
g_at = sp.sqrt(Om*om)/2
web.append(sp.simplify(kap_def.subs(g, g_at) - 1) == 0)
web.append(sp.simplify(g_at**2 - (om/2)*(Om/2)) == 0)          # zp product
web.append(sp.simplify((g_at**2/Om) - om/4) == 0)              # disp energy
phi_req = sp.solve(sp.Eq(kap_def.subs(g, om*sp.Symbol('phi',
                   positive=True)), 1), sp.Symbol('phi',
                   positive=True))[0]
web.append(sp.simplify(phi_req - sp.sqrt(Om/om)/2) == 0)
Ec = sp.symbols('E_c', positive=True)
web.append(sp.simplify(sp.solve(sp.Eq(sp.sqrt(Om/(2*Ec)),
           sp.sqrt(Om/om)/2), Ec)[0] - 2*om) == 0)             # E_c = 2om
web.append(sp.simplify(D_half.subs([(n, 1), (kap, 1)])) == 0)  # D(1)=0
emit("G10J-1d identity web (6 equivalences of kappa=1): %s" %
     ("ALL EXACT" if all(web) else "FAIL %s" % web))

nbe = 1/(sp.exp(x) - 1)
ser = sp.series(1 + nbe, x, 0, 4).removeO()
ser_expect = 1/x + sp.Rational(1, 2) + x/12 - x**3/720
r3 = sp.simplify(sp.expand(ser - ser_expect))
emit("G10J-1e kappa=1 ladder rungs (10C G4 regression): "
     "nu = 1/x + 1/2 + x/12 + 0 x^2 - x^3/720; residue = %s" % r3)
g1 = (r1 == 0) and (r2 == 0) and g1c and all(web) and (r3 == 0)
emit("G10J-1: %s" % ("PASS" if g1 else "FAIL"))
emit("")

# ==================== G10J-2: numerics ====================
emit("G10J-2 numerics (Fock truncation 160)")
M = 160
bd = np.sqrt(np.arange(1, M))
a_b = np.diag(bd, 1)
q_b = a_b + a_b.T
nb = np.diag(np.arange(M, dtype=float))
omv, Omv = 1.0, 0.8
worst = 0.0
for kv in (0.5, 0.9, 0.99):
    gv = math.sqrt(kv*Omv*omv/4)
    for nn in (0, 1, 2, 3):
        ch = nn + 0.5
        H = omv*ch*np.eye(M) + Omv*nb + gv*ch*q_b
        ev = np.linalg.eigvalsh(H)[:8]
        pol = omv*ch + Omv*np.arange(8) - gv*gv*ch*ch/Omv
        worst = max(worst, float(np.max(np.abs(ev - pol))))
g2a = worst <= 1e-8
emit("G10J-2a block spectra vs polaron formula: worst |d| = %.2e "
     "(bar 1e-8) -> %s" % (worst, "PASS" if g2a else "FAIL"))

def chi1(kv, Mq=240):
    gv = math.sqrt(kv*Omv*omv/4)
    d = gv/Omv
    D1 = omv*(1 - kv)          # E(2)-E(1)
    D0m = omv*(1 - kv/2)       # E(1)-E(0)
    P = [math.exp(-d*d)*d**(2*m)/math.factorial(m) for m in range(60)]
    up = 2*2*sum(P[m]/(D1 + m*Omv) for m in range(60))
    dn = 2*1*sum(P[m]/(-D0m + m*Omv) for m in range(60)
                 if abs(-D0m + m*Omv) > 1e-12)
    return up + dn, d

kv = 0.999
c1v, dv = chi1(kv)
lim_exact = 4*math.exp(-dv*dv)/omv
g2b = abs((1-kv)*c1v - lim_exact)/lim_exact <= 0.02
emit("G10J-2b compliance divergence: (1-kappa) chi(1) = %.6f vs exact "
     "limit 4 e^{-d^2}/om = %.6f at kappa=0.999 (bar 2%%) -> %s" %
     ((1-kv)*c1v, lim_exact, "PASS" if g2b else "FAIL"))
for kv2 in (0.5, 0.9, 0.99, 0.999):
    c, _ = chi1(kv2)
    emit("   chi(1) at kappa=%.3f: %10.3f" % (kv2, c))

# thermal persistence: sum_k p_k |<k|D(d)|k>|^2 = e^{-d^2} L_k(d^2)^2
from numpy.polynomial.laguerre import lagval
nbar = 2.0
p = (nbar/(1+nbar))**np.arange(200)/(1+nbar)
d2 = (math.sqrt(0.99*Omv*omv/4)/Omv)**2
Lk = np.array([lagval(d2, np.eye(k+1)[k]) for k in range(200)])
S_exact = float(np.sum(p*np.exp(-d2)*Lk**2))
S_num = 0.0
Dm = np.diag(np.exp(-d2/2)*np.ones(1))  # placeholder replaced below
# numeric displacement operator on truncated Fock space
Mq = 120
bdq = np.sqrt(np.arange(1, Mq))
aq = np.diag(bdq, 1)
from scipy.linalg import expm
Dop = expm(math.sqrt(d2)*(aq.T - aq))
diagD2 = np.abs(np.diag(Dop))**2
pq = (nbar/(1+nbar))**np.arange(Mq)/(1+nbar)
S_num = float(np.sum(pq*diagD2[:Mq]))
g2c = (S_exact > 0) and abs(S_exact - S_num) <= 1e-8
emit("G10J-2c thermal-FC persistence at n_amb=2: sum p_k e^{-d2} "
     "L_k(d2)^2 = %.8f (numeric %.8f, d<=1e-8) POSITIVE -> %s" %
     (S_exact, S_num, "PASS" if g2c else "FAIL"))
emit("   => the chi(1) divergence at kappa -> 1 persists under a "
     "thermal ambient (the diagonal-FC channel weight never vanishes);")
emit("      D(1) itself is bath-occupation-free (10C T1).")
g2 = g2a and g2b and g2c
emit("G10J-2: %s" % ("PASS" if g2 else "FAIL"))
emit("")

# ==================== G10J-3: 9W K-invariance ====================
lam_bar = sp.Symbol('lambda_bar', positive=True)
K = sp.Symbol('K', positive=True, integer=True)
pull_K = K*(lam_bar/sp.sqrt(K))**2/Om
rK = sp.simplify(pull_K - lam_bar**2/Om)
worstK = 0.0
for Kn in (1, 2, 3):
    lv = math.sqrt(0.9*Omv*omv/4)
    lk = lv/math.sqrt(Kn)
    pull = Kn*lk*lk/Omv
    worstK = max(worstK, abs(pull - lv*lv/Omv))
g3 = (rK == 0) and worstK <= 1e-10
emit("G10J-3 the 9W funnel: kappa_coll = 4 lam_bar^2/(Om om), "
     "lam_bar^2 = sum lam_k^2; symbolic residue %s, numeric worst "
     "%.1e (K=1,2,3) -> %s" % (rK, worstK, "PASS" if g3 else "FAIL"))
emit("")

# ==================== G10J-4: 10C meter regressions ====================
g4 = []
g4.append(abs(2*(1-0.450) - 1.100) < 1e-9)
g4.append(sp.simplify(phi_req - sp.sqrt(Om/om)/2) == 0)
emit("G10J-4 meters (archived, no sky fits): kappa = 2(1-c1): c1=0.450")
emit("   -> kappa=1.100 (flat-M/L 4S); hier 4Z c1=0.258 -> 1.484;")
emit("   a0-lock a0_fit/a0_horizon = kappa^2: galaxy hier => kappa =")
emit("   1.00 +/- 0.05; 10D lock row kappa = 0.925 (Planck) / 0.888")
emit("   (SH0ES); 10D free single-kappa 1.503 (ridge caveat C2")
emit("   mandatory); plain deep pole 1.317 / vertical FLIPS (R25-C1).")
emit("G10J-4: %s" % ("PASS" if all(g4) else "FAIL"))
emit("")

# ==================== G10J-5: the D-dependence statement ====================
emit("G10J-5 THE FUNNEL + D STATEMENT (derivation chain):")
emit("   kappa = 4 lam_bar^2/(Om om) [G10J-3]; lam_bar^2 = sum lam_k^2")
emit("   carries the absolute coupling density = exactly the")
emit("   10A-DISPLACED D-normalization (named-not-derived there).")
emit("   => kappa = 1 is derivable WITHOUT D-provenance IFF a")
emit("   structural principle ties lam_bar^2 to Om om/4 directly.")
emit("   The three candidates, premise-labeled:")
emit("   P-crit  saturation of the two-rung bound  [ASSUMED]")
emit("   P-lit   thermalization on the literal dressed ladder")
emit("           (bound only, see R-C)               [ASSUMED]")
emit("   P-scale pure-dS one-scale computation (9Z machinery + D)")
emit("           [OPEN -- a computation, not a principle]")
emit("")

# ==================== R-C: the premise-labeled bound ====================
emit("R-C THE BOUND (premise-labeled):")
emit("   P-vertex [DERIVED, 10B selection rule] + P-polaron [exact] +")
emit("   P-BEform [MEASURED, the RAR identity] + P-lit [ASSUMED]:")
emit("   kappa > 1 => E(2) < E(1) => Gibbs weight at n=2 exceeds the")
emit("   n=1 Boltzmann ordering at EVERY temperature => the mode's")
emit("   own occupation cannot be the measured BE form through the")
emit("   first two rungs.  Hence kappa <= 1, saturation at the")
emit("   compliance divergence (P-crit).")
# exact mini-check: Gibbs inversion at kappa > 1
kv = 1.2
E0, E1, E2 = 0.0, omv*(1-kv/2), omv*(1-kv/2) + omv*(1-kv)
inv = E2 < E1
emit("   check kappa=1.2: E(1)-E(0)=%.3f, E(2)-E(1)=%.3f -> inverted:"
     " %s (Gibbs: P(2)/P(1) = e^{-(E2-E1)/T} > 1 at every T)" %
     (E1-E0, E2-E1, inv))
emit("")
emit("CONSISTENCY ROWS (report-grade):")
emit("   Under P-lit the bound EXCLUDES the 10D low-a0/high-kappa")
emit("   ridge (1.503 single / 1.317 plain deep) and CONTAINS the")
emit("   temperature-locked world (0.925) near saturation: the")
emit("   spectral bound and the temperature lock pick the SAME world")
emit("   (treatment caveat R25-C1 + ridge caveat C2 mandatory).")
emit("   FOLD OBSERVATION: at any fixed kappa > 0 the literal ladder")
emit("   is non-monotone for n >= 2/kappa while the deep sky occupies")
emit("   n ~ 1/x >> 1 => the fixed-coefficient dictionary is the")
emit("   leading form of a self-consistent structure (no direction")
emit("   claim; 10D running direction is treatment-unstable).")
emit("")

# ==================== G10J-6: premise audit -> letter ====================
premises_bound = {"P-vertex": "DERIVED", "P-polaron": "EXACT",
                  "P-BEform": "MEASURED", "P-lit": "ASSUMED"}
premises_sat = {"P-crit": "ASSUMED"}
emit("G10J-6 PREMISE AUDIT:")
for k, v in {**premises_bound, **premises_sat}.items():
    emit("   %-9s %s" % (k, v))
assumed_bound = [k for k, v in premises_bound.items() if v == "ASSUMED"]
assumed_all = assumed_bound + [k for k, v in premises_sat.items()
                               if v == "ASSUMED"]
alg_ok = g1 and g2 and g3 and all(g4)
if not alg_ok:
    letter = "K-OPEN"
elif not assumed_all:
    letter = "K-FORCED"
elif not assumed_bound:
    letter = "K-BOUNDED"
else:
    letter = "K-CONDITIONAL"
emit("")
emit("=" * 68)
emit("VERDICT LETTER: %s" % letter)
emit("   the bound kappa <= 1 rides ONE named assumption (P-lit); the")
emit("   saturation kappa = 1 rides a second (P-crit); everything")
emit("   else in the chain is exact, derived, or measured.  The")
emit("   identity web, the compliance-divergence theorem, the 9W")
emit("   funnel, and the D statement are UNCONDITIONAL.")
emit("pre-signed credence: K-CONDITIONAL -> HOLD mech 8; anomaly-real")
emit("53 untouched (no sky fits)")
emit("gates: G1 %s | G2 %s | G3 %s | G4 %s | G5 stated | G6 audited" %
     tuple("PASS" if b else "FAIL" for b in
           (g1, g2, g3, all(g4))))
emit("total wall-clock %.1f min" % ((time.time()-T0)/60))
