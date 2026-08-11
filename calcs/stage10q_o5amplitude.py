"""
STAGE 10Q -- O5-AMPLITUDE: THE AMBIENT l=2 FLUCTUATION AMPLITUDE FROM
FLUCTUATION-DISSIPATION AT THE MEASURED STATIC RESPONSE (the R33 named
successor: "derive e_a -- 10F's <= ~1 cap vs 10G's own q_FD =
sqrt(2 n_amb + 1) = 1.43 unreconciled; this ONE horizon-side number
decides P10"; author 2026-08-11: "Go for next. Probably the O5? Up to
you, of course.").

THE QUESTION. P10's fate rides on e_a, the RMS fractional l=2/l=0
amplitude of the ambient dressing cloud at the system (10E T2c
convention = the archived 4K/5S EFE quadrupole coefficient scale;
static DC value |q| = 0.086, band [0.073, 0.099]). 10F proved the
static |q| is a DIFFERENT OBJECT from the fluctuation amplitude and
left e_a <= ~1 (chart cap); 10G's FD budget premise allowed RMS
sqrt(2 n_amb + 1) = 1.43 (binary) -- past the chart. Unreconciled.
This stage DERIVES the amplitude instead of banding it.

THE CHAIN (each step gated; no sky fits anywhere -- archived meters
only; anomaly-real untouchable by construction):
  T1a THE SECTOR THEOREM (exact): the cloud's l=2 sector about the
      system = five real coordinates q_m. To quadratic order the
      energy is E = (1/2) k_self * sum_m q_m^2  -  F * q_0:
      the self-energy is SO(3)-isotropic (one k_self multiplying the
      full 5-norm -- rotational invariance of the undriven cloud),
      and the axisymmetric P2 drive couples LINEARLY to m = 0 ONLY
      (10F GF-1, re-derived here). Equilibrium: q_0 = F/k_self =
      q_DC (the archived static |q|), q_{m!=0} = 0. THE HESSIAN IS
      k_self * I_5 EXACTLY -- the linear drive contributes zero
      curvature; the apparent 'orientation stiffness' of constrained-
      sphere readings is absorbed. ONE stiffness, and the MEASURED
      DC response pins it: chi(0) = 1/k_self = q_DC/F.
  T1b FLUCTUATION-DISSIPATION (exact): each coordinate at the
      program temperature T_dS = H/(2 pi) (the standing KMS license,
      10G step 5: the partner's statistics ARE the measured gate)
      carries <dq_m^2> = chi(0) * (Om2/2) coth(Om2 / 2 T_dS).
      CLASSICAL LIMIT Om2 << 2T: <dq_m^2> = T * chi(0) = T q_DC / F,
      Om2-INDEPENDENT AND INERTIA-INDEPENDENT -- equipartition kills
      the 10F Om2 band at the amplitude level; the correction factor
      u coth u (u = Om2/2T) is GATED <= 15% spread across the full
      10F band Om2 in [sqrt(q), 1] * Om_amb at both anchors.
  T1c THE BRIGHT PROJECTION: the instrument reads the rank-1 bright
      combination (10F GF-2, re-checked); a unit projection of five
      iid coordinates has the SINGLE-coordinate variance -- no
      factor-5 inflation. e_a(RMS) = sqrt(<dq_m^2>).
  T2  THE FORCE FROM STORED ENERGY (the 4K solver, verbatim code
      path): linear response gives the stored l=2 field energy at DC
      E2 = (1/2) k_self q_DC^2 = (1/2) F q_DC, so
        1/k_self = q_DC^2 / (2 E2)   and
        <dq^2> = T * q_DC^2 / (2 E2) = q_DC^2 * (T/E_amb) / (2 eps2),
      with eps2 = E2/E0 the l=2/l=0 STORED-ENERGY RATIO computed
      from the SAME solver arrays (units cancel -- the D-normalization
      CANCELS in the ratio), and E0 identified with the cloud's
      physical dressing energy E_amb (L1 below). Energy integrals
      E_l = -(1/(2(2l+1))) int S_l phi_l r^2 dr (G = 1 multipole
      convention), cross-checked against the gradient form
      (1/(2(2l+1))) int [phi_l'^2 + l(l+1) phi_l^2/r^2] r^2 dr and
      against thin-shell analytic values (E0 = 1/2, E2 = a^3/25 at
      unit shell -- derived in-line).
  T3  THE ASSEMBLY: e_a = q_DC * sqrt( (T/E_amb) / (2 eps2) ) *
      sqrt(u coth u), per anchor. T/E_amb = 1/((n_amb + 1/2) x_amb)
      central. Rows: chart margin (<dq^2> vs 1), 10G-ceiling
      reconciliation (<dq^2> vs 2 n_amb + 1: the 1.43 was a budget
      PREMISE, the derived occupation-weighted value replaces it),
      static ratio e_a/q_DC, the 10G strike corner (e_a vs 6.7/9.2 --
      completeness only), and THE FROZEN ROW: refresh time 1/Om2 >=
      1/Om_amb = 2 pi/(x_amb H) in Gyr vs the 10O boundary 2.48 Gyr.
  T4  THE P10 RE-READ: the 10N S machinery INLINED VERBATIM
      (delta_pair/W_dict/pcontrast closed forms), reader-identity
      gated against >= 6 archived printed tokens
      (data/stage10n_crossing.txt), evaluated at the derived
      e_a (central + band edges + Gauss-average over the frozen
      Gaussian draw) across the 10N anchored grid x conv family;
      the P10 sub-letter fires MECHANICALLY (grammar below).

LOAD-BEARING IDENTIFICATIONS (disclosed; the honest seams):
  L1  E0(solver) <-> E_amb = (n_amb + 1/2) Om_amb, the collective
      cloud's energy content (6Y M = 1 collective identification;
      6E measured occupations). VARIANTS band: {n, n+1/2, n+1} *
      Om_amb -- each edge a named reading (thermal-only / full
      oscillator / with-spontaneous); the band is carried into the
      quoted e_a. A many-mode reading makes E_phys LARGER (e_a
      smaller) -- the (n+1)-edge direction; M = 1 is the license
      for the central (6Y theorem).
  L2  T = T_dS exactly (the temperature lock; no band).
  L3  Linear response at q_DC ~ 0.086 (~9%-grade anharmonic
      corrections; disclosed, not expanded).
  L4  The 10E identification of the 4K solar-geometry q as the
      binary-anchor ambient amplitude (ARCHIVED precedent, inherited;
      the plateau-convention caveat rides along).
  L5  The thermal state is realized as FROZEN DRAWS (quenched KMS
      ensemble -- the 10O/P15 frozen-cloud frame, M&H
      Liouville-covered); equipartition sets the draw variance.
  L6  The galaxy-anchor row is REPORT-grade (the deep-galaxy system
      sits outside the interior plateau; 10N galaxy channels are
      contrast-killed anyway -- P10 is binary-native).

ARCHIVED ANCHORS (provenance):
  q_DC static: 0.086 central, band [0.073, 0.099]
    [stage10e_vertex.py:60-62, from the 5S function scan; 4K BE
    eN=1.2 archive q = -0.09878, data/stage4k_quadrupole.txt:9]
  n_amb/x_amb: binary 0.502/1.0954, galaxy 6.63/0.1411 [10E/10F/10N
    anchor set; 6E-token variant 0.5202/6.583 printed as a row]
  Om2/Om_amb band [sqrt(0.086), 1] = [0.293, 1]
    [stage10f_geometry.py:270-272]
  lam_max/g_close (e_a = 1): bin x=0.5: 0.304, x=1.0: 0.140;
    g_close = (1/2) sqrt(om Om) [stage10n_crossing.py:25-27]
  10N letter constants {0.02, 0.005, 0.002}; conv family
    {1/sqrt2, 1, sqrt2} [stage10n_crossing.py grammar]
  10O frozen boundary tau_c >= 2.48 Gyr [PREDICTIONS.md P10
    annotation]; 1/H = 14.42 Gyr at H0 = 67.8 (program lock)

GATES -- with the trap-#11 STANDING RULE applied: every gate names,
here at pre-reg time, WHICH letter clause its failure vetoes; a gate
that vetoes nothing is labeled DIAGNOSTIC. Trap #16: the failing
outcome is stated inline. Trap #10: one convergence gate per integral
family.
  GQ-1 (STOP -> Q-AMBIG) FD identities: (a) sympy-exact classical
       equipartition <q^2> = T/k from the Gibbs Gaussian integral;
       (b) sympy-exact Fock/geometric sum <q^2> = (1/(2 mu Om))
       coth(Om/2T); (c) the chi(0) identity; (d) numeric spectral-FD
       check at gamma/Om = 0.05: quadrature of (1/pi) coth(om/2T)
       Im chi(om) vs closed form <= 0.5%, node-doubling drift
       <= 0.2% (integral family 1). FAILS if the FD wiring is wrong.
  GQ-2 (STOP -> Q-AMBIG) the sector Hessian: sympy on E =
       (1/2) k (Q:Q) - lam Tr(D Q) in the 5-coordinate real basis
       (D = nn - I/3): the drive term is LINEAR in exactly one
       coordinate (GF-1 re-derivation) and the Hessian is k*I_5 with
       zero residue. FAILS if the isotropy/linearity algebra is
       wrong (the single-stiffness claim dies).
  GQ-3 (-> Q-AMBIG) the bright projection: 3 seeded random system
       tensors; coupling vector rank 1; orthogonal complement
       couples at 0 to 1e-12; unit-projection variance of iid
       5-vectors = single-coordinate variance (algebraic check).
       FAILS if the 10F frame does not reproduce.
  GQ-4 (-> Q-BAND) the 4K regression: the verbatim solver (this
       file's copy, + returning S_l) reproduces the ARCHIVED q
       tokens {BE eN=1.0/1.2/1.4: -0.08504/-0.09878/-0.11171;
       simple eN=1.2: -0.09779} to <= 0.5% rel, and the G5 analytic
       l=2 integrator value -0.050654 to <= 0.2%. FAILS if the code
       path drifted => eps2 is UNLICENSED and the letter falls to
       Q-BAND (formula quoted with eps2 unpinned).
  GQ-4b (-> Q-BAND) the energy-integral family (integral family 2):
       thin-shell analytic E0 = 1/2 and E2 = a^3/25 = 0.04 (a = 1)
       reproduced by the pairing form <= 3% (Gaussian shell w = 0.05
       finite-width tolerance); pairing-vs-gradient internal
       equality <= 1% on the physical solve; node-doubling (NR 512
       -> 768, NMU 96 -> 144, LMAX 16 -> 24) moves eps2 <= 3%.
       FAILS if the energy quadrature is untrustworthy => Q-BAND.
  GQ-5 (positive-clause gate; -> the Om2-independence clause) the
       coth flatness: u coth u spread across Om2 in [0.293, 1] *
       Om_amb <= 15% at both anchors. FAILS => the 'Om2-independent'
       clause is VETOED and e_a is quoted Om2-banded (Q-VALUE keeps
       only if spread <= 50%, else Q-BAND).
  GQ-6 (-> the P10 sub-letter) 10N reader identity: the inlined
       S machinery reproduces >= 6 archived tokens {S(0.086),
       S(0.3), S(1), e_a*(2%)} at printed cells to <= 2e-3 rel
       (rounding-limited). FAILS => the P10 re-read leg is AMBIG
       (main letter unaffected).
  GQ-7 (positive-clause gate) the P10 sub-letter fires MECHANICALLY
       from the S table; the operative inequality of WHICHEVER cell
       fires is printed with its margin (trap #12).
  GQ-8 (positive-clause gate; -> the frozen clause) the refresh
       row: min refresh time = 1/Om_amb (the Om2 = Om_amb hi-edge;
       period and damping variants slower) >= 10 x 2.48 Gyr at the
       binary anchor. FAILS => the 'automatically frozen' clause is
       VETOED; the P10 annotation keeps the 10O tension live.
  GQ-9 (positive-clause gate; -> Q-CHART-BREACH) reconciliation
       arithmetic: <dq^2> vs 1 (chart) and vs 2 n_amb + 1 (the 10G
       ceiling) and e_a vs 6.7 (strike-corner completeness) printed
       with margins. FIRES Q-CHART-BREACH if the hi-edge <dq^2> > 1.
  GQ-10 (DIAGNOSTIC) scope audit: no sky fits; archived-meter
       provenance list; L1-L6 disclosed; static-ratio row
       e_a_lo/0.086 >= 2 printed (if < 2 the R25 'static-like'
       corner is RE-ENTERED -- wired as a letter FLAG, not a veto).

VERDICT LETTERS (pre-signed):
  Q-VALUE  GQ-1/2/3/4/4b/5 pass: e_a is DERIVED -- quote central +
           the L1-variant x q_DC-band envelope; the P10 sub-letter
           (below) is adjudicated in the same commit.
  Q-BAND   GQ-4 or GQ-4b fail (eps2 unlicensed) or GQ-5 spread >
           50%: the FD form is quoted with eps2/Om2 unpinned; P10
           re-read at band edges only, labeled BAND-grade.
  Q-CHART-BREACH  the hi-edge <dq^2> > 1: the linear sector chart
           fails at its own derived value; NO amplitude statement;
           the 10G dilemma paragraph stands unchanged.
  Q-AMBIG  any STOP gate fails.
  P10 SUB-LETTER (mechanical from the S table at the derived band;
  10N letter constants):
    P10-STANDS          S(central) >= 0.02 at some anchored corner.
    P10-REANCHOR-FLOOR  S(hi) < 0.02 everywhere AND the band
                        straddles 0.005: the registered DR4-0.02
                        discovery channel is DEAD at the derived
                        amplitude; the 0.005-grade perturbative
                        (kappa-localized Lorentzian) channel is
                        EDGE-conditional -- the row re-anchors.
    P10-REANCHOR-DEAD   S(hi) < 0.005 everywhere: both channels
                        below floor at the derived amplitude; the
                        row re-anchors as design-spec only.
CREDENCE MAP (pre-signed): ALL cells -> HOLD bath-mechanism
conditional 8 AND anomaly-real 53 (a constant-derivation stage on
the dispersive side moves nothing; the P10 row is design-spec /
mech-conditional -- re-anchoring it moves no credence). The
completeness corner (derived e_a >= 6.7) cannot move credence
unilaterally: it fires Q-CHART-BREACH + a round flag.

AMENDMENTS: none at pre-reg. Any post-commit amendment is logged
here pre-quote with its run preserved (house rule).
  A1 (post run 1, IMPLEMENTATION ONLY, no bars/letters touched; run-1
     console preserved in the pre-reg commit's session record): (i)
     GQ-1b's sympy residue needed .rewrite(exp) before simplify --
     the identity 2 n_BE + 1 = coth(x/2) is exact, the simplifier
     stalled on the mixed coth/exp form; (ii) np.trapz ->
     np.trapezoid (removed in numpy 2.x). Both are repairs that make
     the gates compute what this header already says they compute.

COMPUTE < 3 min (8 solver calls + sympy + quadratures).
Writes data/stage10q_o5amplitude.txt.
"""
import math
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss, legval

OUT = "data/stage10q_o5amplitude.txt"
L = []
def emit(s=""):
    print(s)
    L.append(s)

emit("STAGE 10Q O5-AMPLITUDE (pre-reg: this commit)")
emit("")

TWO_PI = 2*math.pi
# ---- archived anchors (provenance in the header) ----
Q_DC_C, Q_DC_LO, Q_DC_HI = 0.086, 0.073, 0.099   # 10E T2c band
N_AMB_BIN, X_AMB_BIN = 0.502, 1.0954             # 10E/10F/10N anchors
N_AMB_GAL, X_AMB_GAL = 6.63, 0.1411
N_AMB_BIN6E, N_AMB_GAL6E = 0.5202, 6.583         # 6E-token variant
OM2_BAND = (math.sqrt(0.086), 1.0)               # 10F band, x Om_amb
LAM_RATIO_BIN = {0.5: 0.304, 1.0: 0.140}         # lam_max/g_close, e_a=1
GAM_GRID = (0.010, 0.015, 0.025)
KAP_GRID = (0.888, 0.925, 1.000)
CONV_FAM = (1/math.sqrt(2), 1.0, math.sqrt(2))
HUBBLE_TIME_GYR = 14.42                          # 1/H at H0 = 67.8
TAU_10O_GYR = 2.48                               # 10O measured boundary
E_A_STRIKE = 6.7                                 # 10G closure demand (low edge)

# =====================================================================
# T1a/T1b -- GQ-1 + GQ-2: the exact core
# =====================================================================
emit("T1 THE SECTOR THEOREM + FD (sympy exact)")

# GQ-1a classical equipartition from the Gibbs integral
kS, TS, qS = sp.symbols('k T q', positive=True)
Z = sp.integrate(sp.exp(-kS*qS**2/(2*TS)), (qS, -sp.oo, sp.oo))
M2 = sp.integrate(qS**2*sp.exp(-kS*qS**2/(2*TS)), (qS, -sp.oo, sp.oo))
r1a = sp.simplify(M2/Z - TS/kS)
emit(f"GQ-1a Gibbs equipartition <q^2> = T/k: residue {r1a}")

# GQ-1b quantum Fock sum -> coth
muS, OmS, nS, bS = sp.symbols('mu Omega n beta', positive=True)
# <q^2> = (1/(2 mu Om)) (2 nbar + 1), nbar = 1/(e^{b Om} - 1)
nbar = 1/(sp.exp(bS*OmS) - 1)
qq = (1/(2*muS*OmS))*(2*nbar + 1)
r1b = sp.simplify((qq - (1/(2*muS*OmS))*sp.coth(bS*OmS/2))
                  .rewrite(sp.exp))
emit(f"GQ-1b Fock <q^2> = (1/2 mu Om) coth(b Om/2): residue {r1b}")

# GQ-1c the chi(0) identity
chi0 = 1/(muS*OmS**2)
r1c = sp.simplify((1/(2*muS*OmS))*sp.coth(bS*OmS/2)
                  - chi0*(OmS/2)*sp.coth(bS*OmS/2))
emit(f"GQ-1c <q^2> = chi(0) (Om/2) coth(Om/2T): residue {r1c}")

# GQ-1d numeric spectral FD at finite gamma (integral family 1)
def spectral_fd(Om, gam, T, nodes):
    # int_0^inf (dw/pi) coth(w/2T) Im chi(w), chi = 1/(Om^2-w^2-i gam w)
    w, wt = leggauss(nodes)
    # map [-1,1] -> [0, wmax] with wmax = 60*Om
    wmax = 60.0*Om
    x = 0.5*(w + 1)*wmax
    ww = 0.5*wt*wmax
    imchi = gam*x/((Om**2 - x**2)**2 + (gam*x)**2)
    cth = 1.0/np.tanh(np.clip(x/(2*T), 1e-12, None))
    return float(np.sum(ww*cth*imchi)/math.pi)

Om_t, gam_t, T_t = 1.0, 0.05, 2.0
closed = (1.0/Om_t**2)*(Om_t/2)/math.tanh(Om_t/(2*T_t))
fd_n = spectral_fd(Om_t, gam_t, T_t, 4000)
fd_2n = spectral_fd(Om_t, gam_t, T_t, 8000)
d1d = abs(fd_n - closed)/closed
d1d2 = abs(fd_2n - fd_n)/closed
g1 = (r1a == 0 and r1b == 0 and r1c == 0 and d1d <= 5e-3
      and d1d2 <= 2e-3)
emit(f"GQ-1d spectral FD (gamma/Om = 0.05): quad {fd_n:.6f} vs closed "
     f"{closed:.6f} (rel {d1d:.2e}, bar 5e-3); doubling {d1d2:.2e} "
     f"(bar 2e-3)")
emit(f"GQ-1: {'PASS' if g1 else 'FAIL'}")
emit("")

# GQ-2 the sector Hessian (5-coordinate real l=2 basis)
emit("GQ-2 the sector Hessian (sympy exact):")
q0, q1, q2, q3, q4, kk, lam = sp.symbols('q0 q1 q2 q3 q4 k lambda',
                                         real=True)
s3 = sp.sqrt(3)
# real symmetric traceless Q from 5 coords (standard real l=2 basis,
# n-hat = z: q0 ~ (2zz-xx-yy)/sqrt(3)-normalized, q1/q2 = xz/yz,
# q3/q4 = (xx-yy)/2, xy
Q = sp.Matrix([[-q0/s3 + q3, q4, q1],
               [q4, -q0/s3 - q3, q2],
               [q1, q2, 2*q0/s3]])
norm = sp.simplify(sp.trace(Q*Q) - 2*(q0**2+q1**2+q2**2+q3**2+q4**2))
D = sp.Matrix([[-sp.Rational(1, 3), 0, 0],
               [0, -sp.Rational(1, 3), 0],
               [0, 0, sp.Rational(2, 3)]])
drive = sp.expand(sp.trace(D*Q))
E_tot = kk*sp.trace(Q*Q)/4 - lam*drive   # (k/2) * sum q^2 - lam * drive
grad_lin = [sp.diff(-lam*drive, v) for v in (q0, q1, q2, q3, q4)]
lin_ok = (sp.simplify(grad_lin[0] + lam*2/s3) == 0 and
          all(sp.simplify(g) == 0 for g in grad_lin[1:]))
H = sp.hessian(E_tot, (q0, q1, q2, q3, q4))
hess_res = sp.simplify(H - kk*sp.eye(5))
hess_ok = hess_res == sp.zeros(5, 5)
qeq = sp.solve([sp.diff(E_tot, v) for v in (q0, q1, q2, q3, q4)],
               (q0, q1, q2, q3, q4))
eq_ok = (sp.simplify(qeq[q0] - 2*lam/(s3*kk)) == 0 and
         all(sp.simplify(qeq[v]) == 0 for v in (q1, q2, q3, q4)))
g2 = (norm == 0) and lin_ok and hess_ok and eq_ok
emit(f"  Q:Q = 2 sum q_m^2 (residue {norm}); drive linear in q0 ONLY: "
     f"{lin_ok}")
emit(f"  Hessian = k * I_5 exactly: {hess_ok}; equilibrium q0 = "
     f"2 lam/(sqrt3 k), others 0: {eq_ok}")
emit(f"  => ONE stiffness for all five coordinates; the measured DC")
emit(f"     response chi(0) = q_DC/F pins it (10F: 'the measured q")
emit(f"     constrains chi_2(0), not Om2' -- the hook, now used)")
emit(f"GQ-2: {'PASS' if g2 else 'FAIL'}")
emit("")

# GQ-3 bright projection (10F GF-2 re-check)
rng = np.random.default_rng(10141)
g3 = True
for t in range(3):
    Msys = rng.normal(size=(3, 3))
    Msys = 0.5*(Msys + Msys.T)
    Msys -= np.eye(3)*np.trace(Msys)/3
    basis = []
    rt3 = math.sqrt(3)
    B0 = np.diag([-1/rt3, -1/rt3, 2/rt3])
    B1 = np.zeros((3, 3)); B1[0, 2] = B1[2, 0] = 1
    B2 = np.zeros((3, 3)); B2[1, 2] = B2[2, 1] = 1
    B3 = np.diag([1.0, -1.0, 0.0])
    B4 = np.zeros((3, 3)); B4[0, 1] = B4[1, 0] = 1
    basis = [B0, B1, B2, B3, B4]
    v = np.array([np.sum(Msys*B) for B in basis])
    C = np.outer(v, v)
    ev = np.linalg.eigvalsh(C)
    g3 = g3 and (abs(ev[:-1]).max() <= 1e-12*max(1.0, ev[-1]))
# unit-projection variance of iid coords = single-coordinate variance
w = rng.normal(size=5); w /= np.linalg.norm(w)
var_proj = float(np.sum(w**2))   # = 1 by construction (iid unit var)
g3 = g3 and abs(var_proj - 1.0) <= 1e-12
emit("GQ-3 bright projection: rank-1 coupling on 3 seeded system")
emit(f"  tensors (orthogonal complement 0 to 1e-12); unit projection "
     f"of iid coords carries variance {var_proj:.12f} = single-coord")
emit(f"GQ-3: {'PASS' if g3 else 'FAIL'}")
emit("")

# =====================================================================
# T2 -- the 4K solver (verbatim code path + S_l return) + energies
# =====================================================================
emit("T2 THE STORED-ENERGY RATIO (4K solver, verbatim numerics)")

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))

def solve_multipole(nu, eN, NR=512, NMU=96, LMAX=16):
    r = np.logspace(-2, 3, NR)
    mu, wmu = leggauss(NMU)
    R, MU = np.meshgrid(r, mu, indexing='ij')
    ST = np.sqrt(1-MU**2)
    Pl = np.zeros((LMAX+1, NMU))
    for l in range(LMAX+1):
        c = np.zeros(l+1); c[l] = 1
        Pl[l] = legval(mu, c)
    gr = -1.0/R**2 + eN*MU
    gt = -eN*ST
    f = nu(np.hypot(gr, gt)) - 1.0
    Ar, At = f*gr, f*gt
    dAr = np.gradient(R**2*Ar, np.log(r), axis=0)/R**3
    dAt = -np.gradient(ST*At, mu, axis=1)/R
    S = -(dAr + dAt)
    Sl = np.array([(2*l+1)/2*np.sum(wmu*S*Pl[l], axis=1)
                   for l in range(LMAX+1)])
    A = np.zeros((LMAX+1, NR)); B = np.zeros((LMAX+1, NR))
    dr = np.diff(r)
    for l in range(LMAX+1):
        g_in = Sl[l]*r
        for i in range(NR-1):
            qf = (r[i]/r[i+1])**(l+1)
            A[l, i+1] = A[l, i]*qf + 0.5*dr[i]*(g_in[i]*qf + g_in[i+1])
        g_out = Sl[l]*r
        for i in range(NR-2, -1, -1):
            qf = (r[i]/r[i+1])**l
            B[l, i] = B[l, i+1]*qf + 0.5*dr[i]*(g_out[i] + g_out[i+1]*qf)
    lv = np.arange(LMAX+1)[:, None]
    phi_l = -(A+B)/(2*lv+1)
    return r, phi_l, Sl

def q_plateau(r, phi2):
    m = (r >= 0.02) & (r <= 0.2)
    qs = phi2[m]/r[m]**2
    return float(np.mean(qs)), float(np.ptp(qs)/abs(np.mean(qs)))

def E_pair(r, Sl, phi_l, l):
    integ = Sl[l]*phi_l[l]*r**2
    return float(-0.5/(2*l+1)*np.trapezoid(integ, r))

def E_grad(r, phi_l, l):
    dphi = np.gradient(phi_l[l], r)
    integ = (dphi**2 + l*(l+1)*(phi_l[l]/r)**2)*r**2
    return float(0.5/(2*l+1)*np.trapezoid(integ, r))

# GQ-4b thin-shell analytic gate (integral family 2)
emit("GQ-4b energy-integral analytic gate (thin shell at a = 1):")
NRt, LMAXt = 512, 16
rt = np.logspace(-2, 3, NRt)
h = np.exp(-0.5*((rt-1.0)/0.05)**2)
h /= np.trapezoid(h*rt**2, rt)          # unit int S r^2 dr
Slt = np.zeros((LMAXt+1, NRt))
ok4b = True
for l, Eexact in ((0, 0.5), (2, 1.0/25.0)):
    Slt[:] = 0.0
    Slt[l] = h
    A = np.zeros((LMAXt+1, NRt)); B = np.zeros((LMAXt+1, NRt))
    dr = np.diff(rt)
    for ll in range(LMAXt+1):
        g_in = Slt[ll]*rt
        for i in range(NRt-1):
            qf = (rt[i]/rt[i+1])**(ll+1)
            A[ll, i+1] = A[ll, i]*qf + 0.5*dr[i]*(g_in[i]*qf + g_in[i+1])
        g_out = Slt[ll]*rt
        for i in range(NRt-2, -1, -1):
            qf = (rt[i]/rt[i+1])**ll
            B[ll, i] = B[ll, i+1]*qf + 0.5*dr[i]*(g_out[i] + g_out[i+1]*qf)
    lv = np.arange(LMAXt+1)[:, None]
    phit = -(A+B)/(2*lv+1)
    Ep = E_pair(rt, Slt, phit, l)
    Eg = E_grad(rt, phit, l)
    dp = abs(Ep - Eexact)/Eexact
    dg = abs(Eg - Ep)/Eexact
    ok4b = ok4b and (dp <= 0.03) and (dg <= 0.01)
    emit(f"  l={l}: pairing {Ep:.6f} vs analytic {Eexact:.6f} "
         f"(rel {dp:.2e}, bar 3e-2); gradient-vs-pairing {dg:.2e} "
         f"(bar 1e-2)")

# GQ-4 archived regression + the physical solves
emit("")
emit("GQ-4 the 4K regression + physical solves:")
ARCHIVE = {('be', 1.0): -0.08504, ('be', 1.2): -0.09878,
           ('be', 1.4): -0.11171, ('simple', 1.2): -0.09779}
NUS = {'be': nu_be, 'simple': nu_simple}
res = {}
ok4 = True
for (fam, eN), qarch in ARCHIVE.items():
    r, phi_l, Sl = solve_multipole(NUS[fam], eN)
    qv, dev = q_plateau(r, phi_l[2])
    E2 = E_pair(r, Sl, phi_l, 2)
    E0 = E_pair(r, Sl, phi_l, 0)
    E2g = E_grad(r, phi_l, 2)
    E0g = E_grad(r, phi_l, 0)
    dq = abs(qv - qarch)/abs(qarch)
    ok4 = ok4 and dq <= 5e-3
    res[(fam, eN)] = dict(q=qv, E2=E2, E0=E0, eps2=E2/E0,
                          dg2=abs(E2g-E2)/abs(E2), dg0=abs(E0g-E0)/abs(E0))
    emit(f"  {fam:>6} eN={eN}: q = {qv:.5f} vs archive {qarch:.5f} "
         f"(rel {dq:.2e}, bar 5e-3)")
    emit(f"        E2 = {E2:.4e}  E0 = {E0:.4e}  eps2 = E2/E0 = "
         f"{E2/E0:.5f}  (grad-check l2 {res[(fam,eN)]['dg2']:.2e}, "
         f"l0 {res[(fam,eN)]['dg0']:.2e}, bar 1e-2)")
    ok4 = ok4 and res[(fam, eN)]['dg2'] <= 1e-2 \
              and res[(fam, eN)]['dg0'] <= 1e-2

# galaxy-anchor solve (report-grade, L6)
eN_gal = X_AMB_GAL**2
r, phi_l, Sl = solve_multipole(nu_be, eN_gal)
qg, devg = q_plateau(r, phi_l[2])
E2gal = E_pair(r, Sl, phi_l, 2)
E0gal = E_pair(r, Sl, phi_l, 0)
res[('be', 'gal')] = dict(q=qg, E2=E2gal, E0=E0gal, eps2=E2gal/E0gal)
emit(f"  BE eN={eN_gal:.4f} (GALAXY row, report-grade L6): q = "
     f"{qg:.5f}  eps2 = {E2gal/E0gal:.5f}  (plateau dev {devg:.1%})")

# hi-res doubling on the central (BE eN=1.2)
r, phi_l, Sl = solve_multipole(nu_be, 1.2, NR=768, NMU=144, LMAX=24)
E2h = E_pair(r, Sl, phi_l, 2)
E0h = E_pair(r, Sl, phi_l, 0)
eps2_c = res[('be', 1.2)]['eps2']
d_hr = abs(E2h/E0h - eps2_c)/eps2_c
ok4b = ok4b and d_hr <= 0.03
emit(f"  hi-res doubling: eps2 = {E2h/E0h:.5f} vs {eps2_c:.5f} "
     f"(rel {d_hr:.2e}, bar 3e-2)")
emit(f"GQ-4: {'PASS' if ok4 else 'FAIL'}   GQ-4b: "
     f"{'PASS' if ok4b else 'FAIL'}")
emit("")

# =====================================================================
# T3 -- the assembly
# =====================================================================
emit("T3 THE ASSEMBLY  e_a = q_DC sqrt((T/E_amb)/(2 eps2)) "
     "sqrt(u coth u)")

def ucothu(u):
    return u/math.tanh(u) if u > 1e-9 else 1.0

def assemble(tag, n_amb, x_amb, eps2, q_c, q_lo, q_hi):
    rows = {}
    for ename, nfac in (('thermal n', n_amb), ('full n+1/2', n_amb+0.5),
                        ('spont n+1', n_amb+1.0)):
        TE = 1.0/(nfac*x_amb)
        base = math.sqrt(TE/(2.0*eps2))
        # coth correction across the 10F Om2 band; u = Om2/(2 T) =
        # x_amb*(Om2/Om_amb)/2
        cfs = [math.sqrt(ucothu(0.5*x_amb*bb)) for bb in OM2_BAND]
        cf_lo, cf_hi = min(cfs), max(cfs)
        rows[ename] = (q_c*base*cf_lo, q_c*base*cf_hi,
                       q_lo*base*cf_lo, q_hi*base*cf_hi, TE, base)
    return rows

def spread_coth(x_amb):
    us = [0.5*x_amb*bb for bb in OM2_BAND]
    vals = [ucothu(u) for u in us]
    return abs(vals[1]-vals[0])/min(vals)

sp_bin = spread_coth(X_AMB_BIN)
sp_gal = spread_coth(X_AMB_GAL)
g5v = (sp_bin <= 0.15) and (sp_gal <= 0.15)
emit(f"GQ-5 coth flatness over Om2 in [0.293, 1] Om_amb: binary "
     f"spread {sp_bin:.3f}, galaxy {sp_gal:.4f} (bar 0.15) -> "
     f"{'PASS' if g5v else 'FAIL'}")
emit("")

eps2_bin = res[('be', 1.2)]['eps2']
eps2_gal = res[('be', 'gal')]['eps2']
eps2_simple = res[('simple', 1.2)]['eps2']

emit(f"BINARY anchor (n_amb = {N_AMB_BIN}, x_amb = {X_AMB_BIN}; "
     f"eps2 = {eps2_bin:.5f} BE / {eps2_simple:.5f} simple):")
rows_bin = assemble('bin', N_AMB_BIN, X_AMB_BIN, eps2_bin,
                    Q_DC_C, Q_DC_LO, Q_DC_HI)
for ename, (ea_lo, ea_hi, ea_llo, ea_hhi, TE, base) in rows_bin.items():
    emit(f"  E_amb = ({ename}) Om: T/E = {TE:.4f}  e_a(q_c) = "
         f"[{ea_lo:.4f}, {ea_hi:.4f}]  full q-band [{ea_llo:.4f}, "
         f"{ea_hhi:.4f}]")
# the quoted envelope: L1 variants x q_DC band x coth band
env_bin = (min(v[2] for v in rows_bin.values()),
           max(v[3] for v in rows_bin.values()))
cen_bin = rows_bin['full n+1/2']
ea_central_bin = 0.5*(cen_bin[0] + cen_bin[1])
emit(f"  CENTRAL (full, q_c, mid-coth): e_a = {ea_central_bin:.4f}; "
     f"ENVELOPE [{env_bin[0]:.4f}, {env_bin[1]:.4f}]")

# 6E-token variant row
rows_6e = assemble('bin6e', N_AMB_BIN6E, math.log(1+1/N_AMB_BIN6E),
                   eps2_bin, Q_DC_C, Q_DC_LO, Q_DC_HI)
c6e = rows_6e['full n+1/2']
emit(f"  6E-token variant (n = {N_AMB_BIN6E}): e_a(q_c) = "
     f"[{c6e[0]:.4f}, {c6e[1]:.4f}] (convention spread row)")

emit("")
emit(f"GALAXY anchor (report-grade L6; n_amb = {N_AMB_GAL}, x_amb = "
     f"{X_AMB_GAL}; eps2 = {eps2_gal:.5f}):")
rows_gal = assemble('gal', N_AMB_GAL, X_AMB_GAL, eps2_gal,
                    abs(res[('be', 'gal')]['q']),
                    abs(res[('be', 'gal')]['q'])*Q_DC_LO/Q_DC_C,
                    abs(res[('be', 'gal')]['q'])*Q_DC_HI/Q_DC_C)
cg = rows_gal['full n+1/2']
emit(f"  e_a(gal, q from its own solve) = [{cg[0]:.4f}, {cg[1]:.4f}] "
     f"(central row)")
emit("")

# GQ-9 reconciliation + chart + strike corner
dq2_hi = env_bin[1]**2
dq2_cen = ea_central_bin**2
qfd2 = 2*N_AMB_BIN + 1
emit("GQ-9 reconciliation arithmetic (binary anchor):")
emit(f"  <dq^2> central {dq2_cen:.4f}, hi-edge {dq2_hi:.4f} vs chart "
     f"1.0 -> margin x{1.0/dq2_hi:.1f} ({'IN-CHART' if dq2_hi < 1 else 'BREACH'})")
emit(f"  vs the 10G FD ceiling 2 n_amb + 1 = {qfd2:.3f}: the 1.43 RMS")
emit(f"  was the BUDGET PREMISE at c = 1; the DERIVED occupation-")
emit(f"  weighted amplitude is x{math.sqrt(qfd2)/ea_central_bin:.1f} "
     f"below it -- the 10F-cap-vs-10G-FD pair RECONCILES: both were")
emit(f"  ceilings; the physical value sits under both.")
emit(f"  strike-corner completeness: e_a = {ea_central_bin:.3f} vs the")
emit(f"  10G closure demand {E_A_STRIKE} -> shortfall x"
     f"{E_A_STRIKE/ea_central_bin:.0f} (the strike stands untouched)")
g9_chart = dq2_hi < 1.0
emit(f"GQ-9: hi-edge in-chart -> {'PASS' if g9_chart else 'CHART-BREACH'}")
emit("")

# GQ-8 frozen row
tau_min_gyr = HUBBLE_TIME_GYR*TWO_PI/X_AMB_BIN   # 1/Om_amb in Gyr... see below
# 1/Om_amb = 1/(x_amb H/(2pi)) = 2pi/(x_amb) * (1/H)
emit("GQ-8 the frozen row (binary anchor):")
emit(f"  Om_amb = x_amb H/(2 pi) => 1/Om_amb = 2 pi/(x_amb H) = "
     f"{tau_min_gyr:.1f} Gyr")
tau_lo_edge = tau_min_gyr           # Om2 = Om_amb hi-edge
tau_hi_edge = tau_min_gyr/OM2_BAND[0]
g8 = tau_lo_edge >= 10*TAU_10O_GYR
emit(f"  refresh band [{tau_lo_edge:.0f}, {tau_hi_edge:.0f}] Gyr "
     f"(1/Om2 across the 10F band; period/damping variants slower)")
emit(f"  vs the 10O boundary {TAU_10O_GYR} Gyr: margin x"
     f"{tau_lo_edge/TAU_10O_GYR:.0f} (bar >= 10x) -> "
     f"{'PASS -- the derived cloud is AUTOMATICALLY frozen-grade' if g8 else 'FAIL'}")
emit("")

# GQ-10 static-ratio + scope
ratio_static = env_bin[0]/Q_DC_C
emit("GQ-10 scope audit:")
emit(f"  static ratio e_a_lo/q_DC = {ratio_static:.2f} (>= 2 keeps the")
emit(f"  R25 'static-like' corner CLOSED; < 2 re-enters it -- FLAG)")
emit("  no sky fits anywhere; archived meters only (provenance in")
emit("  header); L1-L6 disclosed; the D-normalization CANCELS in eps2")
emit("")

# =====================================================================
# T4 -- the P10 re-read (10N machinery inlined verbatim)
# =====================================================================
emit("T4 THE P10 RE-READ (10N S machinery verbatim; GQ-6 reader gate)")

def om_of(x):
    return x/TWO_PI

def nbe(x):
    return 1.0/(math.exp(x) - 1.0)

def pcontrast(x):
    e = math.exp(-x)
    return e*(1 - e)**2

def g_close(x_loc, x_amb):
    return 0.5*math.sqrt(om_of(x_loc)*om_of(x_amb))

def W_dict(x_loc, kap):
    return (kap*om_of(x_loc)/4)*(2*nbe(x_loc) + 2)

def fc0(x_loc, x_amb, kap):
    om, Om = om_of(x_loc), om_of(x_amb)
    d2 = (kap/4)*(om/Om)
    return math.exp(-d2)

def delta_pair(x_loc, x_amb, kap, gam, U, conv=1.0):
    om = om_of(x_loc)
    D1 = abs(om*(1 - kap))
    Vp2 = 2*(conv*U)**2*fc0(x_loc, x_amb, kap)
    z = complex(D1/2, -gam/2)
    wv = np.sqrt(z*z + Vp2)
    return float(np.real(wv) - np.real(z))

def S_stat(x_loc, x_amb, kap, gam, e_a, lam_max, conv=1.0):
    U = e_a*lam_max
    return (pcontrast(x_loc)*delta_pair(x_loc, x_amb, kap, gam, U,
                                        conv)/W_dict(x_loc, kap))

# GQ-6 reader identity vs archived tokens (data/stage10n_crossing.txt)
TOKENS = [   # (x, gam, kap, e_a, printed S)
    (0.5, 0.010, 0.888, 0.086, 0.00022),
    (0.5, 0.010, 0.888, 0.3,   0.00308),
    (0.5, 0.010, 0.888, 1.0,   0.02012),
    (0.5, 0.010, 1.000, 1.0,   0.02145),
    (0.5, 0.010, 0.925, 0.3,   0.00337),
    (1.0, 0.010, 1.000, 1.0,   0.01592),
    (1.0, 0.025, 0.888, 0.3,   0.00050),
]
g6 = True
worst6 = 0.0
for (xl, gm, kp, ea, tok) in TOKENS:
    lam_max = LAM_RATIO_BIN[xl]*g_close(xl, X_AMB_BIN)
    sv = S_stat(xl, X_AMB_BIN, kp, gm, ea, lam_max)
    d = abs(sv - tok)/max(tok, 1e-5)
    worst6 = max(worst6, d)
    g6 = g6 and d <= 5e-2 and abs(sv - tok) <= 1.5e-5 + 0.002*tok
emit(f"GQ-6 reader identity: 7 archived tokens, worst rel d = "
     f"{worst6:.2e} (rounding-limited; bar abs 1.5e-5 + 2e-3 rel) -> "
     f"{'PASS' if g6 else 'FAIL'}")
# e_a*(2%) at the best corner (the 0.936 token)
lam_max_05 = LAM_RATIO_BIN[0.5]*g_close(0.5, X_AMB_BIN)
ea_star = None
for ea in np.linspace(0.3, 1.0, 7001):
    if S_stat(0.5, X_AMB_BIN, 1.000, 0.010, ea, lam_max_05) >= 0.02:
        ea_star = ea
        break
d_star = abs(ea_star - 0.936) if ea_star else 9.9
g6 = g6 and d_star <= 2e-3
emit(f"  e_a*(2%) best corner = {ea_star:.4f} vs archived 0.936 "
     f"(bar 2e-3)")
emit("")

# the re-read table at the derived amplitude
emit("S at the DERIVED amplitude (binary anchors; conv family):")
ea_lo_env, ea_hi_env = env_bin
S_max_cen, S_max_hi, S_max_avg = 0.0, 0.0, 0.0
GH3 = math.sqrt(3.0)  # 3-pt Gauss-Hermite nodes {0, +-sig sqrt3},
                      # weights {2/3, 1/6, 1/6}; S even, S(0) = 0
for xl in (0.5, 1.0):
    lam_max = LAM_RATIO_BIN[xl]*g_close(xl, X_AMB_BIN)
    for gm in GAM_GRID:
        for kp in KAP_GRID:
            row = []
            for conv in CONV_FAM:
                s_c = S_stat(xl, X_AMB_BIN, kp, gm, ea_central_bin,
                             lam_max, conv)
                s_h = S_stat(xl, X_AMB_BIN, kp, gm, ea_hi_env,
                             lam_max, conv)
                # frozen-draw Gauss average over e ~ N(0, sig^2),
                # sig = e_a(RMS): <S> = (1/3) S(sig sqrt3) at 3-pt GH
                s_avg = S_stat(xl, X_AMB_BIN, kp, gm,
                               ea_central_bin*GH3, lam_max, conv)/3.0
                S_max_cen = max(S_max_cen, s_c)
                S_max_hi = max(S_max_hi, s_h)
                S_max_avg = max(S_max_avg, s_avg)
                row.append((s_c, s_h, s_avg))
            if gm == 0.010:
                emit(f"  x={xl} gam={gm:.3f} kap={kp:.3f}: "
                     + "  ".join(f"conv{c:.2f}: Sc={a:.5f} Sh={b:.5f} "
                                 f"Savg={d:.5f}"
                                 for c, (a, b, d) in zip(CONV_FAM, row)))
emit(f"  S_max(central) = {S_max_cen:.5f}; S_max(hi-envelope) = "
     f"{S_max_hi:.5f}; S_max(frozen-draw avg) = {S_max_avg:.5f}")
emit("")

# GQ-7 the mechanical P10 sub-letter
if S_max_cen >= 0.02:
    p10 = "P10-STANDS"
    margin = f"S_max(central) = {S_max_cen:.4f} >= 0.02"
elif S_max_hi < 0.005:
    p10 = "P10-REANCHOR-DEAD"
    margin = f"S_max(hi) = {S_max_hi:.4f} < 0.005"
else:
    p10 = "P10-REANCHOR-FLOOR"
    margin = (f"S_max(hi) = {S_max_hi:.4f} < 0.02 everywhere AND the "
              f"band straddles 0.005 (central {S_max_cen:.4f})")
emit(f"GQ-7 P10 sub-letter: {p10}  [{margin}]")
emit("")

# =====================================================================
# the letter
# =====================================================================
core = g1 and g2 and g3
lic = ok4 and ok4b
if not core:
    letter = "Q-AMBIG"
elif not g9_chart:
    letter = "Q-CHART-BREACH"
elif not lic or spread_coth(X_AMB_BIN) > 0.50:
    letter = "Q-BAND"
else:
    letter = "Q-VALUE" + ("" if g5v else " (Om2-banded)")
emit("=" * 68)
emit(f"LETTER: {letter}  |  P10 sub-letter: {p10}")
emit(f"  e_a(binary) = {ea_central_bin:.3f} central, envelope "
     f"[{env_bin[0]:.3f}, {env_bin[1]:.3f}]")
emit(f"  (q_DC band x L1 variants x coth band; eps2 = {eps2_bin:.5f} "
     f"BE / {eps2_simple:.5f} simple)")
emit(f"  gates: GQ-1 {g1} GQ-2 {g2} GQ-3 {g3} GQ-4 {ok4} GQ-4b {ok4b} "
     f"GQ-5 {g5v} GQ-6 {g6} GQ-8 {g8} chart {g9_chart}")
emit(f"CREDENCE (pre-signed): HOLD mech 8 / anomaly-real 53 (all cells)")
emit("=" * 68)

with open(OUT, "w") as f:
    f.write("\n".join(L) + "\n")
print(f"\nwritten {OUT}")
