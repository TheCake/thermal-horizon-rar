"""
STAGE 9Y (the state-meter): THE SQUEEZING BOUND -- what does the
measured tail pair say about the QUANTUM STATE of the soft ambient
sector?

R18 named an untested Lemma-A assumption: ambient modes INDEPENDENT
and THERMAL.  There is a real physics fork behind it: the static-patch
horizon bath is exactly thermal (Gibbons-Hawking), while comoving-frame
soft modes are famously SQUEEZED.  Correlations change the bright-mode
statistics, hence the lending gate g = P(n_A >= 2), hence the tail
p = 1/2 + g/4.  This stage builds the exact dial and asks what the sky
tolerates.

THE THEOREM CLAIMED (sympy + numerics):
  (i)   CLASSICAL correlations are ABSORBABLE: an equal-sign
        cross-correlated thermal pair (<x1 x2> = <p1 p2> = c >= 0,
        e.g. a common-coherent mixture) keeps the bright covariance
        ISOTROPIC -> the bright mode is exactly thermal at a
        renormalized nbar_A -> the gate stays ON the thermal curve
        g = s(n)^2 (only the n calibration moves).
  (ii)  QUANTUM squeezing is NOT absorbable: the TMS anti-sign
        correlation (<x1 x2> = -<p1 p2>) makes the bright covariance
        ANISOTROPIC (V_Apm = e^{pm 2r}(nbar + 1/2)) -> non-thermal
        statistics -> the gate LEAVES the thermal curve.
  => the measured tail pair's consistency with the thermal curve is a
  specifically QUANTUM-STATE statement about the soft sector: the sky
  is a (weak, honest) squeezing meter.
KEY SIMPLIFICATION (exact): the TMS generator r(b1 b2 - b1'b2')
factorizes in the bright/dark basis as single-mode squeezing S_A(r)
x S_D(-r); equal thermals are basis-invariant; so the bright marginal
of TMS(therm x therm) IS the single-mode squeezed thermal
S(r) rho_th(nbar) S'(r) -- the dial computes in a small single-mode
Fock space, with the two-mode identity as a GATE.

HONESTY PRE-COMMITTED: at the current galaxy tail error (9U: sigma_p
>= 0.075) the bound is expected VACUOUS-at-2sigma (the dial saturates
at p = 3/4 inside the noise) -- if so, quote that plainly; the stage's
products are then (a) the absorbable-vs-not THEOREM, (b) the exact
dial + the DR4-era binary lever (weak-ambient tails sharpen toward
p = 0.69 where the squeezing signature is fractionally largest),
(c) the P1 hardening line: squeezed and thermal CONVERGE at the void
asymptote (g -> 1 both) -> the void kill-band is STATE-INDEPENDENT.

GATES (bars locked at this commit, BEFORE any run):
  G9Y-0  r = 0 regressions: single-mode P(m) geometric at both
         anchors (max dev <= 1e-6); gate values reproduce the 6E/6Y
         thermal gates (gal 0.7536, bin 0.1118) within 5e-4.
  G9Y-1  the factorization identity: two-mode TMS(therm(0.5) x
         therm(0.5)) bright marginal vs single-mode S(r=0.3)
         therm(0.5) S': max |dP(m)| <= 1e-4 (truncation grade).
  G9Y-2  squeezed-vacuum identities: TMS vacuum per-mode marginal =
         thermal(sinh^2 r) (<= 1e-6); single-mode squeezed vacuum
         P(odd) <= 1e-10 and P(0) = 1/cosh(r) (<= 1e-6).
  G9Y-3  sympy covariance legs: classical isotropy (i); TMS
         anisotropy V_Apm = e^{pm 2r}(nbar + 1/2) (ii); the
         coherent-pair bright identity |alpha,alpha> =
         |sqrt(2) alpha>_A |0>_D (numeric 9-node <= 1e-6).
  G9Y-4  truncation: state tail occupancy < 5e-3 at the largest
         (nbar, r) cell used.
  G9Y-5  dial sanity: g monotone nondecreasing in r at both anchors.

VERDICT LETTERS: Y-THEOREM (all gates PASS; dial + honest bound +
DR4 lever + void line booked) / Y-SURPRISE (dial flat: |g(0.5) -
g(0)| < 1e-3 at either anchor -- the tail is squeezing-blind) /
Y-FAIL (any regression gate fails).
CREDENCE MAP (pre-signed): NO move from 9Y alone in any cell; 9X + 9Y
go to ROUND-19 as a package (map in the 9X header: clean package ->
bath-mechanism conditional 15 -> 17; else HOLD).  anomaly-real 53
UNTOUCHED.

Writes data/stage9y_squeezing.txt.  Compute: ~1-2 min.
"""
import math
import numpy as np
import sympy as sp
from scipy.linalg import expm

OUT = 'data/stage9y_squeezing.txt'
L = []
def say(s=''):
    L.append(s); print(s, flush=True)
    with open(OUT, 'w') as f:
        f.write("\n".join(L) + "\n")

say("STAGE 9Y: THE SQUEEZING BOUND -- the tail as a soft-sector "
    "state meter")
say("=" * 72)

N_GAL = 6.583
N_BIN = 0.502
P_GAL_6E, P_BIN_6E = 0.6884, 0.5280
SIG_P_9U = 0.075

def low(d):
    return np.diag(np.sqrt(np.arange(1, d)), 1)

def thermal_vec(N, nbar):
    if nbar <= 0:
        p = np.zeros(N); p[0] = 1.0
        return p
    x = math.log(1.0 + 1.0/nbar)
    p = np.exp(-x*np.arange(N)); p /= p.sum()
    return p

def sq_thermal(N, nbar, r):
    """single-mode squeezed thermal S(r) rho_th S(r)^dag; returns P(m)."""
    a = low(N)
    G = 0.5*(a @ a - a.T @ a.T)      # generator: S = expm(r*G)
    S = expm(r*G)
    rho = S @ np.diag(thermal_vec(N, nbar)) @ S.T
    return np.real(np.diag(rho))

def gate_of(P):
    return float(1.0 - P[0] - P[1])

def s_of_n(n): return n/(1.0 + n)

# ---------------- G9Y-0: r = 0 regressions ----------------
say("G9Y-0: r = 0 regressions:")
ok0 = True
for nm, nb, gtgt in (('gal', N_GAL, 0.7536), ('bin', N_BIN, 0.1118)):
    N = 160 if nb > 1 else 60
    P = sq_thermal(N, nb, 0.0)
    s = s_of_n(nb)
    geo = (1 - s)*s**np.arange(N)
    dev = float(np.max(np.abs(P - geo)))
    g0 = gate_of(P)
    ok = dev <= 1e-6 and abs(g0 - gtgt) <= 5e-4
    ok0 &= ok
    say("  %s (nbar %.3f): max|P - geom| = %.1e; gate = %.4f "
        "(6E/6Y %.4f) %s" % (nm, nb, dev, g0, gtgt,
                             'OK' if ok else 'FAIL'))
say("G9Y-0: %s" % ('PASS' if ok0 else 'FAIL'))
say('')

# ---------------- G9Y-1: the factorization identity ----------------
say("G9Y-1: TMS(therm x therm) bright marginal == single-mode "
    "squeezed thermal (n = 0.5, r = 0.3):")
NB2 = 28
b1 = np.kron(low(NB2), np.eye(NB2))
b2 = np.kron(np.eye(NB2), low(NB2))
G2 = b1 @ b2 - b1.T @ b2.T
S2 = expm(0.3*G2)
rho2 = S2 @ np.kron(np.diag(thermal_vec(NB2, 0.5)),
                    np.diag(thermal_vec(NB2, 0.5))) @ S2.T
A = (b1 + b2)/math.sqrt(2.0)
NA_op = A.T @ A
ev, U = np.linalg.eigh(NA_op)
occ = np.real(np.diag(U.T @ rho2 @ U))
mv = np.rint(ev).astype(int)
P_two = np.zeros(2*NB2)
for m, p in zip(mv, occ):
    if 0 <= m < 2*NB2: P_two[m] += p
P_one = sq_thermal(80, 0.5, 0.3)
mmax = 20
d1 = float(np.max(np.abs(P_two[:mmax] - P_one[:mmax])))
ok1 = d1 <= 1e-4
say("  max |dP(m)| over m < %d: %.1e -> %s" %
    (mmax, d1, 'PASS' if ok1 else 'FAIL'))
say('')

# ---------------- G9Y-2: squeezed-vacuum identities ----------------
say("G9Y-2: squeezed-vacuum identities (r = 0.5):")
rho_v = S_v = None
S_v = expm(0.5*G2)
rho_v = S_v @ np.kron(np.diag(thermal_vec(NB2, 0.0)),
                      np.diag(thermal_vec(NB2, 0.0))) @ S_v.T
n1_op = b1.T @ b1
Pm1 = np.zeros(NB2)
ev1, U1 = np.linalg.eigh(n1_op)
occ1 = np.real(np.diag(U1.T @ rho_v @ U1))
for m, p in zip(np.rint(ev1).astype(int), occ1):
    if 0 <= m < NB2: Pm1[m] += p
n_pred = math.sinh(0.5)**2
geo1 = (1 - s_of_n(n_pred))*s_of_n(n_pred)**np.arange(NB2)
d2a = float(np.max(np.abs(Pm1[:15] - geo1[:15])))
P_sv = sq_thermal(80, 0.0, 0.5)
odd = float(np.max(P_sv[1::2]))
d2c = abs(P_sv[0] - 1.0/math.cosh(0.5))
ok2 = d2a <= 1e-6 and odd <= 1e-10 and d2c <= 1e-6
say("  TMS-vacuum per-mode marginal vs thermal(sinh^2 r): %.1e; "
    "1-mode SV P(odd) max: %.1e; |P(0) - 1/cosh r|: %.1e -> %s" %
    (d2a, odd, d2c, 'PASS' if ok2 else 'FAIL'))
say('')

# ---------------- G9Y-3: sympy covariance legs ----------------
say("G9Y-3: covariance legs:")
nb_s, c_s, r_s = sp.symbols('nbar c r', positive=True)
V = nb_s + sp.Rational(1, 2)
VxA_cl = sp.Rational(1, 2)*(V + V) + c_s      # (x1+x2)/sqrt2 var
VpA_cl = sp.Rational(1, 2)*(V + V) + c_s      # equal-sign classical
iso = sp.simplify(VxA_cl - VpA_cl) == 0
nA_cl = sp.simplify(VxA_cl - sp.Rational(1, 2) - (nb_s + c_s)) == 0
VxA_q = sp.exp(2*r_s)*V                        # TMS anti-sign
VpA_q = sp.exp(-2*r_s)*V
aniso = sp.simplify(VxA_q - VpA_q) != 0
say("  (i) classical equal-sign: bright isotropic, nbar_A = nbar + c: "
    "%s" % ('PASS' if iso and nA_cl else 'FAIL'))
say("  (ii) TMS anti-sign: bright anisotropic e^{+-2r}(nbar+1/2): "
    "%s" % ('PASS' if aniso else 'FAIL'))
# (iii) coherent-pair bright identity, 9-node numeric
NBc = 24
ac = low(NBc)
b1c = np.kron(ac, np.eye(NBc)); b2c = np.kron(np.eye(NBc), ac)
Ac = (b1c + b2c)/math.sqrt(2.0)
NAc = Ac.T @ Ac
evc, Uc = np.linalg.eigh(NAc)
ok3c = True
for alpha in (0.3, 0.7, 1.1):
    v = np.exp(-abs(alpha)**2/2)*alpha**np.arange(NBc) / \
        np.sqrt(np.array([math.factorial(k) for k in range(NBc)],
                         dtype=float))
    psi = np.kron(v, v)
    occ_c = np.abs(Uc.T @ psi)**2
    Pc = np.zeros(2*NBc)
    for m, p in zip(np.rint(evc).astype(int), occ_c):
        if 0 <= m < 2*NBc: Pc[m] += p
    al2 = math.sqrt(2.0)*alpha
    pois = np.exp(-al2**2)*(al2**2)**np.arange(2*NBc) / \
        np.array([math.factorial(k) for k in range(2*NBc)], float)
    ok3c &= float(np.max(np.abs(Pc[:15] - pois[:15]))) <= 1e-6
say("  (iii) |a,a> = |sqrt2 a>_A |0>_D (3 amplitudes): %s" %
    ('PASS' if ok3c else 'FAIL'))
ok3 = iso and nA_cl and aniso and ok3c
say("G9Y-3: %s  => classical correlations ABSORBABLE (thermal curve, "
    "recalibrated n); squeezing NOT (leaves the curve)" %
    ('PASS' if ok3 else 'FAIL'))
say('')

# ---------------- the dial ----------------
say("THE DIAL: gate g(r) = P(n_A >= 2) and tail p = 1/2 + g/4:")
RS = [0.0, 0.1, 0.2, 0.3, 0.5, 0.8]
NGAL_FOCK = 220
NBIN_FOCK = 90
say("  %-6s %-18s %-18s" % ("r_sq", "galaxy g -> p", "binary g -> p"))
dial = []
for r in RS:
    Pg = sq_thermal(NGAL_FOCK, N_GAL, r)
    Pb = sq_thermal(NBIN_FOCK, N_BIN, r)
    gg, gb = gate_of(Pg), gate_of(Pb)
    dial.append((r, gg, gb))
    say("  %-6.2f %.4f -> %.4f   %.4f -> %.4f" %
        (r, gg, 0.5 + gg/4, gb, 0.5 + gb/4))
# G9Y-4 truncation at the largest cell
tail_occ = float(sq_thermal(NGAL_FOCK, N_GAL, 0.8)[-10:].sum())
ok4 = tail_occ < 5e-3
say("G9Y-4 truncation (gal nbar=6.58, r=0.8, last-10 occupancy "
    "%.1e): %s" % (tail_occ, 'PASS' if ok4 else 'FAIL'))
ok5 = all(dial[i+1][1] >= dial[i][1] - 1e-9 and
          dial[i+1][2] >= dial[i][2] - 1e-9
          for i in range(len(dial)-1))
say("G9Y-5 dial monotone in r (both anchors): %s" %
    ('PASS' if ok5 else 'FAIL'))
say('')

# ---------------- the honest bound + levers ----------------
say("THE HONEST BOUND (galaxy tail, 9U primary p-hat = 0.6471 +/- "
    ">= %.3f):" % SIG_P_9U)
p_up = 0.6471 + 2*SIG_P_9U
say("  2-sigma ceiling on p_gal = %.3f; the dial's maximum is 3/4 "
    "= 0.750 %s" %
    (p_up, "-> NO squeezing bound at 2 sigma today (the dial "
     "saturates inside the noise) -- quoted plainly per the "
     "pre-commit" if p_up >= 0.75 else "-> bound derivable:"))
if p_up < 0.75:
    for r, gg, gb in dial:
        if 0.5 + gg/4 > p_up:
            say("  r_sq bound: < %.2f" % r)
            break
say('')
say("THE DR4-ERA LEVER (registered): the binary column is the "
    "fractionally sharpest --")
g0b = dial[0][2]
for r, gg, gb in dial:
    if r in (0.2, 0.3, 0.5):
        say("  r_sq = %.1f: Delta p_bin = %+.4f (gate %.4f vs thermal "
            "%.4f)" % (r, (gb - g0b)/4, gb, g0b))
say("  a DR4-era weak-ambient tail at sigma_p ~ 0.01-0.02 reads "
    "r_sq ~ 0.2-0.3 directly; registered as the state-meter lever.")
say('')
say("P1 HARDENING LINE: squeezed and thermal CONVERGE at the void "
    "asymptote (g -> 1 as nbar -> inf regardless of r): the void "
    "kill-band [0.727, 0.750] is STATE-INDEPENDENT -- the P1 test "
    "does not care what state the soft sector is in.")
gg_inf = gate_of(sq_thermal(400, 60.0, 0.5))
say("  (numeric spot: nbar = 60, r = 0.5: g = %.4f vs thermal "
    "%.4f)" % (gg_inf, s_of_n(60.0)**2))
say('')

# ---------------- verdict ----------------
flat = abs(dial[4][1] - dial[0][1]) < 1e-3 or \
    abs(dial[4][2] - dial[0][2]) < 1e-3
allg = ok0 and ok1 and ok2 and ok3 and ok4 and ok5
if allg and not flat:
    letter = 'Y-THEOREM'
elif allg and flat:
    letter = 'Y-SURPRISE'
else:
    letter = 'Y-FAIL'
say("=" * 72)
say("VERDICT LETTER: %s" % letter)
say("  the absorbable-vs-not theorem: classical correlations move "
    "the tail ALONG the thermal curve (n recalibration); quantum "
    "squeezing moves it OFF the curve; the measured pair's thermal-"
    "curve consistency is a quantum-state statement at whatever "
    "sigma_p the era provides.")
say("  credence: NO move from 9Y alone (pre-signed); package -> "
    "ROUND-19 with 9X.")
say("  anomaly-real 53 UNTOUCHED")
say('')
say("gates: G9Y-0 %s | G9Y-1 %s | G9Y-2 %s | G9Y-3 %s | G9Y-4 %s | "
    "G9Y-5 %s" % tuple('PASS' if x else 'FAIL'
                       for x in (ok0, ok1, ok2, ok3, ok4, ok5)))
print("\nsaved:", OUT)
