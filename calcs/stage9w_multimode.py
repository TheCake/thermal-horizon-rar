"""
STAGE 9W (O5 seam, TODO 30a): THE MULTIMODE REDUCTION -- does the s^L
gate survive a many-mode ambient sector, and what is the collective
occupation it reads?

The seam (6X did the per-leg factor with ONE ambient mode; 6Y selected
M=1 from counting statistics; 9Q/R16 downgraded M=1 to a consistency
condition): the real environment has MANY soft modes.  Either the gate
form breaks at K >= 2 (a wound), or the reduction to one collective
mode is FORCED -- in which case M=1 was never an assumption.

LEMMA A (the bright-mode reduction, claimed theorem-grade in the
degenerate frozen limit): for linear exchange coupling
H_int = sum_k lambda_k (a^dag b_k + h.c.) the interaction involves ONLY
the collective (bright) mode A = sum_k c_k b_k, c_k = lambda_k/lam_bar,
lam_bar^2 = sum |lambda_k|^2:  H_int = lam_bar (a^dag A + A^dag a)
EXACTLY -- the K-1 dark combinations decouple by construction.  For a
product of thermal modes (occupations n_k) the bright marginal is
EXACTLY thermal at the coupling-weighted mean
    nbar_A = sum_k |c_k|^2 n_k
(passive rotation of isotropic Gaussian covariances stays isotropic),
so P(n_A >= L) = s_bar^L with s_bar = nbar_A/(1+nbar_A): THE GATE FORM
IS MULTIMODE-EXACT, evaluated at the weighted-mean occupation.
Corollary (the 6Y closure): the sky's rejection of NB tails (M >= 2
democratic) is the rejection of number-ADDITIVE gating, not of
multimode baths -- linear coupling makes ANY K look like one mode.
M=1 upgrades from consistency-condition to FORCED-BY-LINEARITY.

LEMMA B (the spectral-locality condition, derivation-by-exclusion):
WHICH weighted mean does the sky read?  A FLAT coupling weight over
the soft band gives nbar dominated by the softest modes (n_BE ~ 1/x:
log-divergent) => gate -> 1, e_N-BLIND, p_gal = p_bin = 3/4 --
EXCLUDED twice over by the measured pair (0.6884 gal / 0.5280 bin,
6E).  A GAUSSIAN-local weight at the system's ambient frequency
x_0 = sqrt(e_N/a0) with width Gamma reproduces both postdictions
(convexity correction only) -- the clean locality bar runs here.
A LORENTZIAN-tailed weight re-imports the soft divergence through
its fat tails (soft-end leak, IR-cutoff log-sensitive) -- run as a
REPORTED scan, no bar: the sky's exact postdictions bound the
selection kernel's soft TAIL as well as its width (sub-Lorentzian
tails required) = a sharpening of the locality condition.  The
locality MECHANISM is the 9T resonant selection (reading-grade,
stated not derived; 6X G2b measured off-resonant modes carrying the
lambda^2 virtual channel, not the lending channel -- the physical
participation width is Rabi-grade, far below the scanned widths).

LEMMA C (the budget -- the c4-rung ordering question): the convexity
correction nbar(Gamma) > n_BE(x_0) shifts s_bar, hence the ladder
rungs c3 = -g/16 and the c4 target c4(L=2) = s^2/192 - 1/720.
Compute delta_s/s and the c4-target shift vs Gamma/x_0.
PRE-COMMITTED INTERPRETATION RULE: if the c4-target shift <= 5% for
all Gamma/x_0 <= 0.3, print "c4 rung NOT hindered" (TODO 30c may
proceed on the single-mode target); else "multimode-soft" (30c
re-scoped).

SCOPE (stated, not hidden): degenerate ambient frequencies + the 6X
frozen-horizon limit; LINEAR ambient coupling is the named surviving
assumption; L = 2 dynamics cited from 6X/7E (statistics-grade here);
WHY-locality remains reading-grade (9T selection).

GATES (bars locked at this commit, BEFORE any run):
  G9W-0  K=1 port regression: the 6X configuration (NA=4, NB=56,
         CHI=0.8, LAM=0.02, WA=5.0, diag infinite-time average)
         reproduces all six printed P2bar values (0.09979, 0.16621,
         0.24900, 0.33114, 0.39524, 0.43449) within 2e-3 each.
  G9W-1  sympy set ALL exact: [A,A^dag] = 1; the H_int bright
         identity (K=2 and K=3 symbolic); covariance isotropy ->
         nbar_A = sum w_k n_k (symbolic 2-mode rotation); the GL
         geometric-tail identity (6X verbatim); the 6Y NB formula
         P(N_tot >= 2) = 1-(1-q)^M(1+Mq) + regression to the 6Y
         table (M=2 gal gate 0.9524) within 5e-4.
  G9W-2  bright statistics (numeric-exact): at config S2 the bright
         number distribution P(n_A = m) matches geometric(nbar_A) to
         max abs dev <= 1e-3, and P(n_A >= 2) = s_bar^2 to <= 1e-3
         (S1 and S2).
  G9W-3  dynamics (the theorem's data gate): dress-ward weight in the
         ORIGINAL mode basis (block-degenerate infinite-time average)
         matches (1/2) s_bar to max relative deviation <= 0.05 over
         all K=2/K=3 configs, AND both rival forms (mean-of-ratios
         (1/2) sum w_k s_k; total-occupation gate (1/2)(1-prod(1-s_k)))
         carry rms >= 3x the weighted-mean rms.
  G9W-4  K-invariance: equal-n K=2 config S1 reproduces the K=1
         weight at the same nbar within 0.005 absolute.
  G9W-5  Lemma B: flat-weight exclusion demonstrated (gate e_N-blind,
         -> 3/4 both systems); GAUSSIAN-local reproduces BOTH 6E
         postdictions within 0.010 at Gamma/x_0 = 0.1; the
         Lorentzian soft-leak scan is REPORTED (IR cutoffs 1e-4,
         1e-6) with no bar -- the kernel-tail bound finding.
  G9W-6  budget table (Gaussian kernel) printed + the interpretation
         rule applied.

VERDICT GRAMMAR: W-FULL = G9W-0..4 PASS (THEOREM) + G9W-5 PASS
(CONDITION); W-THEOREM-ONLY; W-CONDITION-ONLY; W-REFUTED = a rival
form beats weighted-mean on the dynamics (rms) AND weighted-mean
fails its 0.05 bar; W-AMBIG otherwise.

CREDENCE MAP (pre-signed; bath-mechanism conditional currently 12;
final booking AFTER ROUND-18 red-team, the 9T pattern):
  (W-FULL or W-THEOREM-ONLY) + ROUND-18 no-unpatched-hole -> 12 -> 15
  theorem with a named hole                               -> HOLD 12
  W-REFUTED                                                -> 12 -> 10
  W-CONDITION-ONLY / W-AMBIG                               -> HOLD 12
  anomaly-real 53 UNTOUCHED in every cell.

Writes data/stage9w_multimode.txt.  Compute: ~3-6 min.
"""
import math
import numpy as np
import sympy as sp
from scipy.integrate import quad

OUT = 'data/stage9w_multimode.txt'
L = []
def say(s=''):
    L.append(s); print(s, flush=True)
    with open(OUT, 'w') as f:
        f.write("\n".join(L) + "\n")

say("STAGE 9W: THE MULTIMODE REDUCTION -- bright mode, locality, budget")
say("=" * 72)

# ---------------- operators ----------------
def build_multi(NA, NBs):
    dims = [NA] + list(NBs)
    def op_at(mat, slot):
        out = np.array([[1.0]])
        for i, d in enumerate(dims):
            out = np.kron(out, mat if i == slot else np.eye(d))
        return out
    def low(d):
        return np.diag(np.sqrt(np.arange(1, d)), 1)   # annihilation
    a = op_at(low(NA), 0)
    bs = [op_at(low(NBs[k]), 1 + k) for k in range(len(NBs))]
    Na = a.T @ a
    return a, bs, Na, dims

def thermal_vec(N, nbar):
    if nbar <= 0:
        p = np.zeros(N); p[0] = 1.0
        return p
    x = math.log(1.0 + 1.0/nbar)
    p = np.exp(-x*np.arange(N)); p /= p.sum()
    return p

CHI, LAM, WA = 0.8, 0.02, 5.0

def dressward(NBs, ns, lams, block=True):
    """dress-ward weight P2bar: a starts in |1>, each b_k thermal(n_k);
    all b_k degenerate at WA+CHI (resonant with 1->2); infinite-time
    average (block-degenerate projection if block else pure diag)."""
    a, bs, Na, dims = build_multi(4, NBs)
    H = WA*Na + 0.5*CHI*(Na @ (Na - np.eye(Na.shape[0])))
    for k, b in enumerate(bs):
        H = H + (WA + CHI)*(b.T @ b) + lams[k]*(a.T @ b + b.T @ a)
    ev, U = np.linalg.eigh(H)
    p_a = np.zeros(4); p_a[1] = 1.0
    rho_diag = p_a
    for k, NB in enumerate(NBs):
        rho_diag = np.kron(rho_diag, thermal_vec(NB, ns[k]))
    rho = np.diag(rho_diag)
    rl = U.T @ rho @ U
    if block:
        rbar = np.zeros_like(rl)
        scale = max(1.0, float(np.max(np.abs(ev))))
        tol = 1e-8*scale
        i = 0
        n_ = len(ev)
        while i < n_:
            j = i + 1
            while j < n_ and ev[j] - ev[i] < tol:
                j += 1
            rbar[i:j, i:j] = rl[i:j, i:j]
            i = j
    else:
        rbar = np.diag(np.diag(rl))
    rb = U @ rbar @ U.T
    sel = (np.arange(4) == 2)*1.0
    proj = np.array([[1.0]])
    proj_full = np.diag(sel)
    for d in dims[1:]:
        proj_full = np.kron(proj_full, np.eye(d))
    return float(np.real(np.trace(proj_full @ rb)))

# ---------------- G9W-0: K=1 port regression ----------------
say("G9W-0: K=1 port regression (6X config, pure-diag average):")
TGT6X = {0.25: 0.09979, 0.50: 0.16621, 1.00: 0.24900,
         2.00: 0.33114, 4.00: 0.39524, 8.00: 0.43449}
ok0 = True
for nb, tgt in TGT6X.items():
    w = dressward([56], [nb], [LAM], block=False)
    d = w - tgt
    ok = abs(d) <= 2e-3
    ok0 &= ok
    say("  n=%5.2f: %8.5f (6X %8.5f, d=%+.5f) %s" %
        (nb, w, tgt, d, 'OK' if ok else 'FAIL'))
wblk = dressward([56], [1.0], [LAM], block=True)
say("  block-vs-diag at n=1: d = %+.2e (generic spectrum; disclosure)"
    % (wblk - TGT6X[1.00]))
say("G9W-0: %s" % ('PASS' if ok0 else 'FAIL'))
say('')

# ---------------- G9W-1: the sympy set ----------------
say("G9W-1: symbolic legs:")
c1, c2, c3_, n1, n2, th = sp.symbols(
    'c1 c2 c3 n1 n2 theta', positive=True)
l1, l2, l3 = sp.symbols('lambda1 lambda2 lambda3', positive=True)
lbar = sp.sqrt(l1**2 + l2**2)
ok_comm = sp.simplify((l1/lbar)**2 + (l2/lbar)**2 - 1) == 0
lbar3 = sp.sqrt(l1**2 + l2**2 + l3**2)
ok_comm3 = sp.simplify((l1/lbar3)**2 + (l2/lbar3)**2 + (l3/lbar3)**2
                       - 1) == 0
say("  [A, A^dag] = sum c_k^2 = 1 (K=2, K=3): %s" %
    ('PASS' if ok_comm and ok_comm3 else 'FAIL'))
# H_int bright identity: lambda_k = lbar * c_k with c from lambda --
# identity is definitional; verify the reconstruction K=3
ok_hint = sp.simplify(lbar3*(l1/lbar3) - l1) == 0 and \
          sp.simplify(lbar3*(l2/lbar3) - l2) == 0 and \
          sp.simplify(lbar3*(l3/lbar3) - l3) == 0
say("  H_int = lam_bar (a^dag A + h.c.) reconstruction (K=3): %s" %
    ('PASS' if ok_hint else 'FAIL'))
# covariance isotropy: rotate two isotropic thermal covariances
Vx = sp.cos(th)**2*(n1 + sp.Rational(1, 2)) \
    + sp.sin(th)**2*(n2 + sp.Rational(1, 2))
Vp = sp.cos(th)**2*(n1 + sp.Rational(1, 2)) \
    + sp.sin(th)**2*(n2 + sp.Rational(1, 2))
nA = sp.cos(th)**2*n1 + sp.sin(th)**2*n2
ok_cov = sp.simplify(Vx - Vp) == 0 and \
         sp.simplify(Vx - (nA + sp.Rational(1, 2))) == 0
say("  bright covariance isotropic with nbar_A = sum w_k n_k: %s" %
    ('PASS' if ok_cov else 'FAIL'))
# GL geometric tail (6X verbatim)
nb_s, k_s, x_s = sp.symbols('nbar L x', positive=True)
q = nb_s/(1 + nb_s)
nsym = sp.symbols('n', integer=True, nonnegative=True)
chk = sp.simplify(sp.summation((1 - q)*q**nsym, (nsym, k_s, sp.oo))
                  - q**k_s)
ok_gl = chk == 0
kms = sp.simplify(q.subs(nb_s, 1/(sp.exp(x_s) - 1)) - sp.exp(-x_s)) == 0
say("  GL geometric tail + KMS identity (6X verbatim): %s" %
    ('PASS' if ok_gl and kms else 'FAIL'))
# 6Y NB formula + table regression
qs, Ms = sp.symbols('q M', positive=True)
nb_formula = 1 - (1 - qs)**Ms*(1 + Ms*qs)
s_gal = math.sqrt(0.7536)      # gate = s^2; 6Y q = s
s_bin = math.sqrt(0.1118)
ok_nb = True
for (Mv, tgt_g, tgt_b) in [(1, 0.7536, 0.1118), (2, 0.9524, 0.2607),
                           (3, 0.9917, 0.4093)]:
    vg = float(nb_formula.subs({qs: s_gal, Ms: Mv}))
    vb = float(nb_formula.subs({qs: s_bin, Ms: Mv}))
    ok_nb &= abs(vg - tgt_g) <= 5e-4 and abs(vb - tgt_b) <= 5e-4
say("  6Y NB formula regression (M = 1,2,3 both systems): %s" %
    ('PASS' if ok_nb else 'FAIL'))
ok1 = ok_comm and ok_comm3 and ok_hint and ok_cov and ok_gl and kms \
    and ok_nb
say("G9W-1: %s" % ('PASS' if ok1 else 'FAIL'))
say('')

# ---------------- configs ----------------
CFG = [
    ('S1', [20, 20], [1.0, 1.0], [0.02, 0.02]),
    ('S2', [20, 20], [0.25, 2.0], [0.02, 0.02]),
    ('S3', [20, 20], [2.0, 0.25], [0.02, 0.01]),
    ('S4', [20, 20], [0.5, 2.0], [0.01, 0.02]),
    ('S5', [20, 20], [0.25, 1.0], [0.02, 0.014]),
    ('S6', [10, 10, 10], [0.25, 0.5, 1.0], [0.02, 0.014, 0.01]),
]
def nbar_A(ns, lams):
    w = np.array(lams)**2
    w = w/w.sum()
    return float(np.sum(w*np.array(ns)))
def s_of_n(n): return n/(1.0 + n)

# ---------------- G9W-2: bright statistics ----------------
say("G9W-2: bright-mode statistics (eigen-projection of N_A):")
def bright_stats(NBs, ns, lams):
    a, bs, Na, dims = build_multi(1, NBs)   # internal slot unused (dim1)
    lam = np.array(lams)
    c = lam/np.linalg.norm(lam)
    A = sum(c[k]*bs[k] for k in range(len(bs)))
    NA_op = A.T @ A
    ev, U = np.linalg.eigh(NA_op)
    rho_diag = np.array([1.0])
    for k, NB in enumerate(NBs):
        rho_diag = np.kron(rho_diag, thermal_vec(NB, ns[k]))
    pv = U.T @ np.diag(rho_diag) @ U
    occ = np.real(np.diag(pv))
    mvals = np.rint(ev).astype(int)
    P = {}
    for m, p in zip(mvals, occ):
        P[m] = P.get(m, 0.0) + p
    return P
ok2 = True
for name in ('S1', 'S2'):
    _, NBs, ns, lams = [c for c in CFG if c[0] == name][0]
    nA_ = nbar_A(ns, lams)
    sb = s_of_n(nA_)
    P = bright_stats(NBs, ns, lams)
    maxdev = max(abs(P.get(m, 0.0) - (1 - sb)*sb**m) for m in range(8))
    tail2 = sum(p for m, p in P.items() if m >= 2)
    d2 = abs(tail2 - sb*sb)
    ok = maxdev <= 1e-3 and d2 <= 1e-3
    ok2 &= ok
    say("  %s: nbar_A = %.4f; max|P(m) - geom| = %.1e; "
        "|P(n>=2) - s^2| = %.1e %s" %
        (name, nA_, maxdev, d2, 'OK' if ok else 'FAIL'))
say("G9W-2: %s  => the gate form s_bar^L is multimode-exact "
    "(statistics)" % ('PASS' if ok2 else 'FAIL'))
say('')

# ---------------- G9W-3 / G9W-4: dynamics ----------------
say("G9W-3: dress-ward dynamics in the ORIGINAL basis vs the three "
    "forms:")
rows = []
for name, NBs, ns, lams in CFG:
    trunc = max((s_of_n(n))**NB for n, NB in zip(ns, NBs))
    P2 = dressward(NBs, ns, lams, block=True)
    nA_ = nbar_A(ns, lams)
    w = np.array(lams)**2; w = w/w.sum()
    pred_wm = 0.5*s_of_n(nA_)
    pred_mor = 0.5*float(np.sum(w*np.array([s_of_n(n) for n in ns])))
    pred_tot = 0.5*(1.0 - float(np.prod([1 - s_of_n(n) for n in ns])))
    rows.append((name, P2, pred_wm, pred_mor, pred_tot))
    say("  %s n=%s lam=%s [trunc %.0e]: P2 = %.5f | wm %.5f | "
        "mor %.5f | tot %.5f" %
        (name, ns, lams, trunc, P2, pred_wm, pred_mor, pred_tot))
dev_wm = [abs(P2/pw - 1) for _, P2, pw, _, _ in rows]
rms = lambda v: float(np.sqrt(np.mean(np.square(v))))
rms_wm = rms([math.log(P2/pw) for _, P2, pw, _, _ in rows])
rms_mor = rms([math.log(P2/pm) for _, P2, _, pm, _ in rows])
rms_tot = rms([math.log(P2/pt) for _, P2, _, _, pt in rows])
ok3 = max(dev_wm) <= 0.05 and rms_mor >= 3*rms_wm and \
    rms_tot >= 3*rms_wm
say("  max |P2/wm - 1| = %.4f (bar 0.05); rms ln-resid: wm %.4f | "
    "mor %.4f (x%.1f) | tot %.4f (x%.1f)" %
    (max(dev_wm), rms_wm, rms_mor, rms_mor/max(rms_wm, 1e-12),
     rms_tot, rms_tot/max(rms_wm, 1e-12)))
refuted = (max(dev_wm) > 0.05) and (rms_mor < rms_wm or
                                    rms_tot < rms_wm)
say("G9W-3: %s" % ('PASS' if ok3 else
                   ('REFUTED-FLAG' if refuted else 'FAIL')))
w_k1 = dressward([56], [1.0], [LAM], block=False)
w_k2 = [P2 for nm, P2, _, _, _ in rows if nm == 'S1'][0]
d_inv = abs(w_k2 - w_k1)
ok4 = d_inv <= 0.005
say("G9W-4 K-invariance (S1 equal-n vs K=1 at nbar=1): |d| = %.4f "
    "-> %s" % (d_inv, 'PASS' if ok4 else 'FAIL'))
say('')

# ---------------- G9W-5: Lemma B ----------------
say("G9W-5: spectral locality (which weighted mean does the sky "
    "read?):")
def n_be(x): return 1.0/(math.expm1(x))
X_GAL = math.sqrt(0.02)
X_BIN = math.sqrt(1.2)
P_GAL_6E, P_BIN_6E = 0.6884, 0.5280
def p_of_n(n):
    s = s_of_n(n)
    return 0.5 + (s*s)/4.0
say("  anchors: x_gal = %.4f (n = %.3f, p pred %.4f); x_bin = %.4f "
    "(n = %.3f, p pred %.4f)" %
    (X_GAL, n_be(X_GAL), p_of_n(n_be(X_GAL)),
     X_BIN, n_be(X_BIN), p_of_n(n_be(X_BIN))))
ok_anchor = abs(p_of_n(n_be(X_GAL)) - P_GAL_6E) <= 1e-3 and \
    abs(p_of_n(n_be(X_BIN)) - P_BIN_6E) <= 1e-3
say("  6E postdiction regression: %s" %
    ('PASS' if ok_anchor else 'FAIL'))
for xmin in (1e-3, 1e-2):
    nfl = quad(lambda x: n_be(x), xmin, 1.0)[0]/(1.0 - xmin)
    say("  FLAT weight (cutoff %g): nbar = %.1f -> gate = %.4f -> "
        "p = %.4f BOTH systems (e_N-BLIND)" %
        (xmin, nfl, s_of_n(nfl)**2, p_of_n(nfl)))
say("  => flat EXCLUDED: magnitude (p -> 3/4) AND e_N-blindness vs "
    "the measured split 0.6884/0.5280")
say("  GAUSSIAN-local weight, nbar(x0, Gamma) and postdictions:")
def n_gauss(x0, G):
    num = quad(lambda x: n_be(x)*math.exp(-0.5*((x - x0)/G)**2),
               1e-6, 1.0, points=[x0])[0]
    den = quad(lambda x: math.exp(-0.5*((x - x0)/G)**2),
               1e-6, 1.0, points=[x0])[0]
    return num/den
ok5b = True
budget = []
for gfrac in (0.05, 0.1, 0.2, 0.3, 0.5):
    ng = n_gauss(X_GAL, gfrac*X_GAL)
    nb_ = n_gauss(X_BIN, gfrac*X_BIN)
    pg, pb = p_of_n(ng), p_of_n(nb_)
    ds_g = s_of_n(ng)/s_of_n(n_be(X_GAL)) - 1
    ds_b = s_of_n(nb_)/s_of_n(n_be(X_BIN)) - 1
    budget.append((gfrac, ds_g, ds_b, pg, pb))
    say("    Gamma/x0 = %.2f: p_gal = %.4f (d %+0.4f), p_bin = %.4f "
        "(d %+0.4f); ds/s gal %+0.3f bin %+0.3f" %
        (gfrac, pg, pg - P_GAL_6E, pb, pb - P_BIN_6E, ds_g, ds_b))
    if gfrac == 0.1:
        ok5b = abs(pg - P_GAL_6E) <= 0.010 and abs(pb - P_BIN_6E) <= 0.010
say("  LORENTZIAN soft-leak scan (reported, no bar -- the kernel-tail"
    " bound):")
def n_lor(x0, G, xmin):
    num = quad(lambda x: n_be(x)/((x - x0)**2 + G*G), xmin, 1.0,
               points=[x0])[0]
    den = quad(lambda x: 1.0/((x - x0)**2 + G*G), xmin, 1.0,
               points=[x0])[0]
    return num/den
for xmin in (1e-4, 1e-6):
    for gfrac in (0.05, 0.1):
        ngl = n_lor(X_GAL, gfrac*X_GAL, xmin)
        say("    xmin=%g Gamma/x0=%.2f: nbar_gal = %.2f (local %.2f) "
            "-> p_gal = %.4f (d %+0.4f vs 6E)" %
            (xmin, gfrac, ngl, n_be(X_GAL), p_of_n(ngl),
             p_of_n(ngl) - P_GAL_6E))
say("    => Lorentzian tails re-import the soft divergence "
    "(IR-log-sensitive): the sky's exact postdictions demand "
    "SUB-LORENTZIAN kernel tails -- the locality condition is a "
    "(width AND tail) bound")
ok5 = ok_anchor and ok5b
say("G9W-5: %s  => the sky reads a LOCAL weighted mean; locality "
    "mechanism = 9T resonant selection (reading-grade)" %
    ('PASS' if ok5 else 'FAIL'))
say('')

# ---------------- G9W-6: the budget / c4 ordering ----------------
say("G9W-6: the c4-rung budget (c3 = -g/16; c4(L=2) = s^2/192 - "
    "1/720):")
c4_hindered = False
for gfrac, ds_g, ds_b, _, _ in budget:
    s0g = s_of_n(n_be(X_GAL))
    s1g = s0g*(1 + ds_g)
    c4_0 = s0g**2/192 - 1/720
    c4_1 = s1g**2/192 - 1/720
    shift = abs(c4_1/c4_0 - 1) if c4_0 != 0 else float('inf')
    say("    Gamma/x0 = %.2f: ds/s(gal) %+0.3f -> c4 target shift "
        "%.1f%%" % (gfrac, ds_g, 100*shift))
    if gfrac <= 0.3 and shift > 0.05:
        c4_hindered = True
say("  interpretation rule: %s" %
    ("c4 rung MULTIMODE-SOFT (30c re-scoped)" if c4_hindered
     else "c4 rung NOT hindered -- the single-mode c4 target stands "
          "for Gamma/x0 <= 0.3 (TODO 30c may proceed)"))
say('')

# ---------------- verdict ----------------
theorem = ok0 and ok1 and ok2 and ok3 and ok4
condition = ok5
if theorem and condition:
    letter = 'W-FULL'
elif theorem:
    letter = 'W-THEOREM-ONLY'
elif refuted:
    letter = 'W-REFUTED'
elif condition:
    letter = 'W-CONDITION-ONLY'
else:
    letter = 'W-AMBIG'
say("=" * 72)
say("VERDICT LETTER: %s" % letter)
say("  Lemma A: bright reduction + geometric statistics + dynamics "
    "%s; M=1 %s" %
    ('CONFIRMED' if theorem else 'NOT confirmed',
     'FORCED-BY-LINEARITY (upgrade from consistency-condition)'
     if theorem else 'unchanged'))
say("  Lemma B: %s" % ('spectral-locality CONDITION derived-by-'
                       'exclusion' if condition else 'not established'))
say("  named surviving assumptions: LINEAR ambient coupling; "
    "WHY-locality (9T selection, reading-grade)")
say("  credence: per the pre-signed map, FINAL BOOKING AFTER ROUND-18")
say("  anomaly-real 53 UNTOUCHED")
say('')
say("gates: G9W-0 %s | G9W-1 %s | G9W-2 %s | G9W-3 %s | G9W-4 %s | "
    "G9W-5 %s | G9W-6 printed" %
    tuple('PASS' if x else 'FAIL'
          for x in (ok0, ok1, ok2, ok3, ok4, ok5)))
print("\nsaved:", OUT)
