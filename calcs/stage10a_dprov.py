"""
STAGE 10A -- O5-NORM: THE D-PROVENANCE DISCHARGE (the sector ledger +
the normalization bridge). Round-20 condition 8 / round-21 condition 5
("pin down D and the closure becomes real"); author call 2026-08-08:
"let's hunt down the derivation, give it all you got."

QUESTION: 9Z fixed the soft-sector coupling density by CONVENTION
(D_conf = omega/4pi^2, charge g = H, window = g) and 9Z-b showed the
whole dispersive closure scales linearly in that normalization, with a
joint pessimal corner (g = 0.3H x D3) past FATAL. What is the PHYSICAL
normalization -- and which channel carries both the O(1) grammar and
the suppressed leak (the round-21 one-seam objection)?

THE PHYSICS, STATED PRE-COMPUTE (the sector-ledger hypothesis):
  Gravity splits into two sectors with different coupling laws.
  (i) CONSTRAINT sector: near-field/first-law couplings. The enclosed
      mass is thermodynamically conjugate to the cosmological-horizon
      entropy (SdS first law dM = -T_c dS_c, EXACT along the family) at
      O(1) conjugacy; Birkhoff forbids any l=0 radiative channel, and
      the horizon's own multipole tower is gapped at sqrt(l(l+1))H
      (R16). This sector can carry O(H) exchange couplings and has NO
      soft radiative continuum attached.
  (ii) RADIATIVE sector: bulk graviton continuum, l >= 2 only, coupling
      to any system coordinate through its PER-QUANTUM multipole charge
      -- suppressed both by eta^(2l) geometry AND by the macroscopic
      effective mass in q_t^2 = hbar/(2 M_eff omega).
  If the measured exchange coupling g_c (the 9U gamma dial) CANNOT be
  supplied by sector (ii), the grammar necessarily rides sector (i),
  the 9Z-b budget's physical normalization is the (tiny) radiative one,
  the closure becomes real, and the round-21 objection dissolves: the
  O(1) grammar and the suppressed leak are carried by DIFFERENT sectors
  of the same field (the cavity-QED geometry: strong coupling to the
  confined collective mode, weak coupling to the free continuum).

INSTRUMENTS:
  T1 CONSTRAINT CREDENTIAL (exact): sympy Birkhoff-with-Lambda (the
     R_tr => static step; the static SdS family solves Einstein+Lambda
     exactly) and the SdS first law dM = -T_c dS_c verified EXACTLY
     (all orders along the family, parameterized by r_c). Corollary:
     the constraint-tower admixture budget (horizon multipoles at
     sqrt(l(l+1))H, KMS occupations, detuning >= sqrt2 H - Omega):
     reported coefficient x g_c^2.
  T2 THE l=2 RADIAL PROBLEM: the l-general minimal potential V_l =
     H^2[l(l+1)/sinh^2(Hx) - 2 sech^2(Hx)] derived symbolically from
     the metric reduction; l=0 must regress to the 9Z well AND the
     numeric amplitude machinery must reproduce the EXACT l=0 ratio
     1 + H^2/omega^2 (the decisive pipeline regression); then l=2 is
     a dial turn: E_l2(omega) = D_l2^dS/D_l2^flat from float64 RK4
     with a 2x-resolution gate (the amplitude invariant u^2 +
     (u'/omega)^2 needs no phase coverage; V decays as e^{-2x}, so
     x_max = 25 suffices at any omega). Soft exponent reported.
  T3 THE AB-INITIO RADIATIVE NUMBER: per-quantum l=2 coupling density
     of the soft coordinate. Charge^2 = |dQbar/ds|^2 q_t^2 with
     |dQbar/ds|^2 = (2/3) mu^2 s^2 exact (sympy, static-averaged
     circular-orbit quadrupole; the full-modulation OVERESTIMATE
     stated), q_t^2 = hbar/(2 M_eff omega), physical branch M_eff >=
     mu (the coordinate moves the orbit's mass; U scales as 1/M_eff
     -- the scaling row is printed, with the crossover M_eff* where
     the radiative channel WOULD suffice). J_rad = Gamma_spont/(2 pi)
     with Gamma_spont fixed by the harmonic-cascade correspondence
     (coherent-state radiated power P_cl = Gamma_spont hbar omega
     n_bar EXACTLY for an oscillator => Gamma = G omega^4 Lam_c^2 /
     (5 c^5 M_eff); hbar cancels -- the audit is printed), gated by
     the sympy circular-orbit regression P = (32/5) G mu^2 s^4
     Omega_orb^6 / c^5 (exact tensor sum). E_l2 from T2 multiplies
     (conservative: enhancements included). U := J_rad(Omega_sys) /
     J_req(Omega_sys) per system, where J_req = A_min_band x
     Dhat_l0(Omega_sys) -- the WEAKEST bridged normalization on the
     measured g_c band, so the split-forcing conclusion is
     band-robust.
  T4 THE NORMALIZATION BRIDGE (the cancellation bound): under 9W
     (bright-mode reduction lam_bar^2 = Sum lam_k^2) + the 9X one-scale
     window, the resonant coupling is the window integral of the SAME
     density: g_c^2 = INT_window J. The absolute normalization then
     CANCELS out of every budget observable: dc1 = g_c^2 x F(geometry,
     shape). Fixed point A(g_c) at W = g_c (self-consistent Rabi
     width; W in {g_c/2, g_c, 2 g_c} robustness rows); the 9Z-b D2
     machinery re-run bit-verbatim with (A, W) freed; the measured
     g_c band from the 9U dial VERBATIM (r = 2(p-1/2)/gate; inst
     r = sin^2(gamma/2), wind r = (1-sin(gamma)/gamma)/2 clipped to
     (0, pi]; gamma = 2 g_c t_hat; t_hat in {1/H, 2pi/H} = an
     irreducible convention ENVELOPE, quoted as such, never
     averaged). Outputs: the dc1_T4(g_c) curve (mean AND deep-end
     per R21 cond 1), the leak eps_T4(g_c) (own gate + smallest
     gate), the kill-window {g_c: |dc1_mean| >= 0.15}, and the four
     convention-central readings. The l=2-quintic shape leg is
     CUTOFF-DOMINATED (secular approximation fails at omega ~
     Omega_orb ~ 1e5 H) -- noted, EXCLUDED from verdict cells by
     pre-registration.

PRE-REGISTERED BARS:
  B-T1: first-law residual == 0 symbolically (else STOP -- textbook
        identity, failure = bug).
  B-T2: l=0 exact regression of the amplitude machinery AND l=2
        2x-resolution convergence.
  B-T3: U_phys(M_eff = mu) <= 1e-6 in BOTH systems => the radiative
        sector CANNOT supply the measured exchange => the sector split
        is FORCED. U in (1e-6, 0.1] => split-supported (GRAY). U > 0.1
        => one-sector reading LIVE (arms the T4 verdict cell).
  B-T4: verdict-relevant ONLY if U > 0.1 (else reported as the
        all-radiative conditional + kill-window). If armed: |dc1_T4|
        mean >= 0.15 at ALL FOUR convention centrals => FATAL-grade.
  B-PHYS: dc1(A_rad, proxy shape) and eps(A_rad) inside the 9Z-b/9Z
        PASS bars by >= 10x (expected by ~30 orders if B-T3 fires).

LETTER GRAMMAR (trap #11 -- every gate wired):
  STOP-class (bug => abort, fix-or-abort disclosed, NO verdict):
    GN-0a Birkhoff steps, GN-0b first law, GN-1 potential + l=0
    regression, GN-3a quadrupole-formula regression, GN-3b charge
    tensor, GN-4a/b 9Z-b/9Z reader identities.
  Wiring:
    GN-2a/b FAIL -> E_l2 := 1 carried into T3 (labeled), letter cap
        N-GRAY.
    GN-5 dial regression FAIL -> g_c band := envelope gamma in
        [1.0, 2.55] x t_hat in {1, 2pi}; labeled ENVELOPE-grade;
        letters otherwise unaffected.
    GT5 fixed-point pathology -> affected sub-band excluded +
        reported; whole band pathological -> T4 unavailable, cap
        N-GRAY.
  LETTERS:
    N-SPLIT-CLOSED: STOPs green, B-T3 fires (U <= 1e-6 both), B-PHYS
        holds, T2 gates green. The D-provenance DISCHARGES: D is a
        radiative-sector normalization, computed ab initio; the
        budget closes at the physical channel; the grammar's O(1)
        coupling lives in the constraint sector; the T4 bound is
        quoted as the all-radiative conditional with its gamma
        kill-window.
    N-GRAY: any mid cell (U in (1e-6, 0.1]; T2 dead; T4 unavailable;
        B-PHYS short of its 10x margin without breaching FATAL).
    N-FATAL: (U > 0.1 AND |dc1_T4| >= 0.15 at all four centrals) OR
        B-PHYS breached at FATAL grade (|dc1_phys| >= 0.15 or
        eps_phys >= 3 x smallest gate).
CREDENCE MAP (pre-signed): N-FATAL -> bath-mechanism conditional
  15 -> 8 immediately. N-SPLIT-CLOSED or N-GRAY -> HOLD 15 pending
  ROUND 22 (rise cell: R22 no-unpatched-hole AND N-SPLIT-CLOSED ->
  15 -> 18, the round-21 "closure becomes real" cell; any hole ->
  hold and execute conditions). anomaly-real UNTOUCHED (no sky number
  moves; the gamma band re-reads archived 9U outputs, no new fits).

DISCLOSURE (in-session sketches, stated in full; none verdict-grade):
  the sector-split hypothesis and the cavity-QED analogy were formed
  in-session BEFORE this pre-reg; the 9Z Part-D eta^4 weights
  (1.5e-47 / 4.5e-29) make B-T3's firing EXPECTED (stated so the
  "forced" cell is not read as a surprise); scratch T3 magnitudes
  U ~ 1e-50 (bin) / 1e-41 (gal) at M_eff = mu were estimated by hand,
  including catching a wrong first pass that used the CLASSICAL
  quadrupole amplitude as the charge (off by the macroscopic/quantum
  amplitude ratio -- the per-quantum lesson is wired into T3's
  cascade-correspondence audit; a second in-session catch fixed a
  factor 2 in Gamma_spont via the same correspondence); the T4
  direction sketch says small-g_c is the dangerous end (window
  shrinks AND curvature grows; scratch vac-scale at g_c ~ 0.107 was
  ~2x fiducial with much curvier shape) and large-g_c is safe
  (residual ~ 1/W^3), with rough scale multipliers 3-11x vs fiducial
  BEFORE affine refit -- the refit changes these materially, no
  verdict number exists pre-run. The fiducial-implied bridged
  coupling was sketched at g_c ~ 0.26-0.33 H. NOT computed
  pre-commit: any E_l2 value, any U value, any dc1_T4 curve point,
  the kill-window, eps values, the tower sum.

AMENDMENTS (post-first-firing, pre-quote; as-fired run archived in
data/stage10a_dprov_run1.txt -- letter as fired: N-GRAY with
GN-2b/GT5 red, U_forced/B-PHYS/dial all green; no bar, letter cell,
or credence-map change):
  A1 (GN-2b instrument): the fixed uniform RK4 step is unconverged
     for l >= 1 (the l(l+1)/x^2 barrier makes the near-origin region
     stiff; first firing d = 8.7e-2 uniformly). Fix: exact Frobenius
     series start at x0 = 0.05 (b_1..b_5 from the recurrence
     b_k 2k(2l+2k+1) = Sum b_j w_{k-1-j}, potential series
     coefficients from sympy) + two-phase integration (dense to
     x = 1, then coarse to 25). The l=0 exact regression re-gates
     the machinery end-to-end.
  A2 (GT5 metric): first firing conflated fixed-point pathology with
     the STATISTIC saturating its domain (at small g_c the residual
     exceeds x in the deep arm, nu(x + r) leaves its domain =
     past-FATAL-by-construction, printed inf). Redefined: GT5 gates
     A(g_c) finiteness; saturated cells print SAT and count as
     FATAL-grade in kill-window and centrals (they already did
     numerically). Also: band-min selection gains a 1e-9 rounding
     guard (first firing dropped the lowest band point, quoting U
     against A = 0.196 instead of 0.131 -- anti-conservative by
     1.5x, irrelevant at 41 orders, fixed); B-PHYS's dc1 is quoted
     as the analytic bound max|nu'| x max|resid| (the direct
     difference underflows double precision at A_rad ~ 1e-42).
  A3 (bonus display, no gate): l=1 E-ratio rows + IR-exponent
     probes, and a closed-form candidate scan E_l =? Prod_k (1 +
     k^2/omega^2) over small integer sets (the run-1 numbers sat a
     uniform ~1.5% from (1+1/om^2)(1+9/om^2) -- if a candidate
     locks at 1e-8 across four frequencies it is printed as a
     numerically-verified conjecture, else the measured values
     stand alone).
Output: data/stage10a_dprov.txt. Wall-clock: ~5-15 min (sympy Ricci
is the slow part).
"""
import csv
import math
import os
import time

import numpy as np
import sympy as sp
from scipy.integrate import quad

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

OUT = []
def P(s=""):
    print(s, flush=True)
    OUT.append(s)

def save():
    with open('data/stage10a_dprov.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

t00 = time.time()
P("=" * 78)
P("STAGE 10A -- O5-NORM: THE D-PROVENANCE DISCHARGE")
P("=" * 78)
gates = {}
labels = []

TWO_PI = 2*math.pi

# =====================================================================
# PART 0 -- T1: the constraint-sector credential (exact)
# =====================================================================
P("")
P("PART 0 -- T1 CONSTRAINT CREDENTIAL")

t_, r_, th_, ph_ = sp.symbols('t r theta phi')
M_, H_ = sp.symbols('M H', positive=True)
X = [t_, r_, th_, ph_]

def christoffel(g, ginv):
    n = 4
    Gam = [[[0]*n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                s_ = 0
                for d in range(n):
                    if ginv[a, d] == 0:
                        continue
                    s_ += ginv[a, d]*(sp.diff(g[d, b], X[c])
                                      + sp.diff(g[d, c], X[b])
                                      - sp.diff(g[b, c], X[d]))
                v = sp.simplify(s_/2)
                Gam[a][b][c] = v
                Gam[a][c][b] = v
    return Gam

def ricci(Gam):
    n = 4
    R = sp.zeros(n, n)
    for b in range(n):
        for c in range(b, n):
            s_ = 0
            for a in range(n):
                s_ += (sp.diff(Gam[a][b][c], X[a])
                       - sp.diff(Gam[a][b][a], X[c]))
                for d in range(n):
                    s_ += (Gam[a][a][d]*Gam[d][b][c]
                           - Gam[a][c][d]*Gam[d][b][a])
            v = sp.simplify(s_)
            R[b, c] = v
            R[c, b] = v
    return R

# ---- GN-0a step 1: general spherical ansatz -> R_tr forces static ---
A_ = sp.Function('A', positive=True)(t_, r_)
B_ = sp.Function('B', positive=True)(t_, r_)
g1 = sp.diag(-A_, B_, r_**2, r_**2*sp.sin(th_)**2)
Gam1 = christoffel(g1, g1.inv())
Ric1 = ricci(Gam1)
Bdot = sp.diff(B_, t_)
ratio_tr = sp.simplify(Ric1[0, 1]/Bdot)
ok_tr = (ratio_tr != 0 and not ratio_tr.has(Bdot)
         and not ratio_tr.has(sp.Derivative))
# ---- GN-0a step 2: the static SdS family solves Einstein+Lambda -----
fS = 1 - 2*M_/r_ - H_**2*r_**2
g2 = sp.diag(-fS, 1/fS, r_**2, r_**2*sp.sin(th_)**2)
Gam2 = christoffel(g2, g2.inv())
Ric2 = ricci(Gam2)
Rs2 = sp.simplify(sum(g2.inv()[i, i]*Ric2[i, i] for i in range(4)))
ok_static = True
for i in range(4):
    for j in range(4):
        e = sp.simplify(Ric2[i, j] - sp.Rational(1, 2)*g2[i, j]*Rs2
                        + 3*H_**2*g2[i, j])
        if e != 0:
            ok_static = False
gates['GN-0a'] = bool(ok_tr and ok_static)
P("GN-0a Birkhoff-Lambda: R_tr/(dB/dt) = %s (nonzero, derivative-free "
  "=> vacuum+Lambda forces static B); SdS solves Einstein+Lambda "
  "exactly: %s  => %s"
  % (ratio_tr, ok_static, "PASS" if gates['GN-0a'] else "FAIL"))

# ---- GN-0b: the SdS first law, exact along the family ---------------
rc = sp.Symbol('r_c', positive=True)
m_of_rc = rc*(1 - H_**2*rc**2)/2               # from f(r_c) = 0
f_r = 1 - 2*m_of_rc/r_ - H_**2*r_**2
kappa = sp.simplify(-sp.diff(f_r, r_).subs(r_, rc)/2)   # > 0 at r_c
T_c = kappa/(2*sp.pi)
S_c = sp.pi*rc**2                               # A/4, G = hbar = c = 1
first_law = sp.simplify(sp.diff(m_of_rc, rc) + T_c*sp.diff(S_c, rc))
gates['GN-0b'] = (first_law == 0)
P("GN-0b SdS first law dM = -T_c dS_c: kappa = %s; residual dM/dr_c + "
  "T_c dS_c/dr_c = %s  => %s"
  % (kappa, first_law, "PASS" if gates['GN-0b'] else "FAIL"))
if not (gates['GN-0a'] and gates['GN-0b']):
    P("STOP: constraint-credential gate failed (bug-class; textbook "
      "identity). Nothing downstream quoted.")
    save(); raise SystemExit(0)

HBAR = 1.0546e-34; C = 2.998e8; GSI = 6.674e-11
H_SI = 2.27e-18
MSUN = 1.989e30
dS_dM = 2*math.pi*C**2/(HBAR*H_SI)
P("  conjugacy magnitude: |dS_c/dM| = 2 pi c^2/(hbar H) = %.3e /kg "
  "-> %.3e nats for 1 M_sun (O(1) thermodynamic conjugacy; NOT "
  "multipole-suppressed)" % (dS_dM, dS_dM*MSUN))
labels.append("T1 BANKED (first law + Birkhoff steps exact)")

# ---- T1 corollary: the constraint-tower admixture coefficient -------
X_BIN = 1.095; X_GAL_LOC = 0.22; X_GAL_AMB = 0.1411
OM_BIN = X_BIN/TWO_PI; OM_GAL = X_GAL_LOC/TWO_PI

def n_be_x(xv):
    if xv > 500: return 0.0
    return 1.0/math.expm1(xv)

tower = 0.0
for l in range(1, 60):
    wl = math.sqrt(l*(l + 1))
    tower += (2*l + 1)*4*n_be_x(TWO_PI*wl)/(wl - OM_BIN)**2
P("  constraint-tower budget: Sum_(l>=1) (2l+1) 4 n(sqrt(l(l+1))) / "
  "(sqrt(l(l+1)) - Omega_bin)^2 = %.3e x g_c^2  (vs smallest gate "
  "0.1117: negligible for any g_c <= 2H)" % tower)

# =====================================================================
# PART 1 -- T2: the l-general radial problem
# =====================================================================
P("")
P("PART 1 -- T2 THE l=2 RADIAL PROBLEM (H = 1 units)")

x, w, Hs = sp.symbols('x w H', positive=True)
l_ = sp.Symbol('l', positive=True, integer=True)
r_of_x = sp.tanh(Hs*x)/Hs
f_of_x = 1 - Hs**2*r_of_x**2
V_l_general = sp.simplify(f_of_x*(l_*(l_ + 1)/r_of_x**2 - 2*Hs**2))
V_target = Hs**2*(l_*(l_ + 1)/sp.sinh(Hs*x)**2 - 2/sp.cosh(Hs*x)**2)
ok_pot = sp.simplify(V_l_general - V_target) == 0
V0 = sp.simplify(V_target.subs(l_, 0))
ok_l0well = sp.simplify(V0 + 2*Hs**2/sp.cosh(Hs*x)**2) == 0
v_min = sp.tanh(Hs*x)*sp.cos(w*x) + (w/Hs)*sp.sin(w*x)
res0 = sp.simplify(sp.diff(v_min, x, 2)
                   + (w**2 + 2*Hs**2/sp.cosh(Hs*x)**2)*v_min)
gates['GN-1'] = ok_pot and ok_l0well and (res0 == 0)
P("GN-1 potential + l=0 regression: V_l = H^2[l(l+1)/sinh^2 - "
  "2 sech^2] %s; l=0 well %s; 9Z solution residual %s => %s"
  % (ok_pot, ok_l0well, res0, "PASS" if gates['GN-1'] else "FAIL"))
if not gates['GN-1']:
    P("STOP: generalized pipeline fails its l=0 regression.")
    save(); raise SystemExit(0)

_WREG_CACHE = {}
def _wreg_coeffs(l, nmax=5):
    """even-series coefficients [c0, c2, c4, ...] of the regular part
    W(x) = l(l+1)(1/sinh^2 x - 1/x^2) - 2 sech^2 x (exact, sympy)."""
    if l in _WREG_CACHE:
        return _WREG_CACHE[l]
    xs = sp.Symbol('x', positive=True)
    expr = l*(l + 1)*(1/sp.sinh(xs)**2 - 1/xs**2) - 2/sp.cosh(xs)**2
    ser = sp.series(expr, xs, 0, 2*nmax).removeO()
    cs = [float(ser.coeff(xs, 2*m)) for m in range(nmax)]
    _WREG_CACHE[l] = cs
    return cs

def dS_amp_ratio(l, om, N=40000, x_max=25.0):
    """1/amplitude^2 for the dS minimal-potential solution normalized
    to near-origin coefficient 1 on x^{l+1}. Amendment A1: exact
    Frobenius series start at x0 = 0.05 (recurrence b_k 2k(2l+2k+1) =
    Sum_j b_j w_{k-1-j}) + two-phase float64 RK4 (dense to x = 1,
    then coarse to x_max; the amplitude invariant u^2 + (u'/om)^2 is
    exact once V ~ 0, and V decays as e^{-2x})."""
    x0 = 0.05
    cs = _wreg_coeffs(l)
    wv = [cs[0] - om*om] + cs[1:]
    b = [1.0]
    for k in range(1, 6):
        acc_ = 0.0
        for j in range(k):
            if k - 1 - j < len(wv):
                acc_ += b[j]*wv[k - 1 - j]
        b.append(acc_/(2.0*k*(2*l + 2*k + 1)))
    u = sum(bk*x0**(l + 1 + 2*k) for k, bk in enumerate(b))
    up = sum((l + 1 + 2*k)*bk*x0**(l + 2*k) for k, bk in enumerate(b))
    def acc(xx, uu):
        sh = math.sinh(xx); ch = math.cosh(xx)
        V = l*(l + 1)/(sh*sh) - 2.0/(ch*ch)
        return (V - om*om)*uu
    for (xa, xb, Nph) in ((x0, 1.0, N//2), (1.0, x_max, N//2)):
        h = (xb - xa)/Nph
        xx = xa
        for _ in range(Nph):
            k1u = up;            k1v = acc(xx, u)
            k2u = up + h/2*k1v;  k2v = acc(xx + h/2, u + h/2*k1u)
            k3u = up + h/2*k2v;  k3v = acc(xx + h/2, u + h/2*k2u)
            k4u = up + h*k3v;    k4v = acc(xx + h, u + h*k3u)
            u, up = (u + h/6*(k1u + 2*k2u + 2*k3u + k4u),
                     up + h/6*(k1v + 2*k2v + 2*k3v + k4v))
            xx += h
    amp2 = u*u + (up/om)**2
    return 1.0/amp2

def flat_amp_ratio(l, om):
    dfact = 1.0
    for k in range(2*l + 1, 0, -2):
        dfact *= k
    return (om**(l + 1)/dfact)**2

P("")
P("  E_l(omega) = dS/flat near-origin density ratio (float64 RK4, "
  "2x-resolution gated):")
ok_regl0 = True
for omv in (0.0553, OM_GAL, OM_BIN, 0.5, 2.0):
    e0 = dS_amp_ratio(0, omv)/flat_amp_ratio(0, omv)
    tgt = 1 + 1/omv**2
    d = abs(e0/tgt - 1)
    ok_regl0 &= (d < 1e-6)
    P("    l=0 om=%.4f : E = %.8f  vs exact 1+1/om^2 = %.8f  (d=%.1e)"
      % (omv, e0, tgt, d))
gates['GN-2a'] = bool(ok_regl0)
P("  GN-2a l=0 exact regression of the amplitude machinery: %s"
  % ("PASS" if ok_regl0 else "FAIL"))

E_l2 = {}
ok_2x = True
for name, omv in (('binary', OM_BIN), ('galaxy', OM_GAL),
                  ('mid', 0.5), ('uv-check', 3.0)):
    e1 = dS_amp_ratio(2, omv, N=40000)
    e2 = dS_amp_ratio(2, omv, N=80000)
    d2x = abs(e1/e2 - 1)
    ok_2x &= (d2x < 1e-6)
    E_l2[name] = e2/flat_amp_ratio(2, omv)
    P("    l=2 om=%.4f : E_l2 = %.6f   (2x-resolution d = %.1e)  [%s]"
      % (omv, E_l2[name], d2x, name))
gates['GN-2b'] = bool(ok_2x)
P("  GN-2b l=2 2x-resolution convergence: %s"
  % ("PASS" if ok_2x else "FAIL"))
eA = dS_amp_ratio(2, 0.02)/flat_amp_ratio(2, 0.02)
eB = dS_amp_ratio(2, 0.01)/flat_amp_ratio(2, 0.01)
soft_slope = math.log(eA/eB)/math.log(2.0)
P("  l=2 soft exponent: d ln E / d ln om between om = 0.01 and 0.02 "
  "= %+.3f  (compare l=0's -2)" % soft_slope)
# A3 bonus: l=1 rows + exponent (display, no gate)
E_l1 = {}
for omv in (OM_BIN, OM_GAL, 0.5, 3.0):
    E_l1[omv] = dS_amp_ratio(1, omv)/flat_amp_ratio(1, omv)
    P("    l=1 om=%.4f : E_l1 = %.6f  (bonus row)" % (omv, E_l1[omv]))
e1A = dS_amp_ratio(1, 0.02)/flat_amp_ratio(1, 0.02)
e1B = dS_amp_ratio(1, 0.01)/flat_amp_ratio(1, 0.01)
P("  l=1 soft exponent: %+.3f"
  % (math.log(e1A/e1B)/math.log(2.0)))
# A3 candidate scan: E_l =? Prod_k (1 + k^2/om^2)
def cand(om, ks):
    return float(np.prod([1 + k*k/(om*om) for k in ks]))
probe2 = [(OM_BIN, E_l2['binary']), (OM_GAL, E_l2['galaxy']),
          (0.5, E_l2['mid']), (3.0, E_l2['uv-check'])]
best2 = None
for ks in ((1, 3), (1, 2), (2, 3), (1, 1), (2, 2), (3, 3)):
    dev = max(abs(cand(om, ks)/ev - 1) for om, ev in probe2)
    if best2 is None or dev < best2[1]:
        best2 = (ks, dev)
if best2[1] < 1e-8:
    P("  A3 CONJECTURE (numerically verified, 4 frequencies, "
      "maxdev %.1e): E_l2 = Prod_k(1 + k^2/om^2), k in %s"
      % (best2[1], str(best2[0])))
else:
    P("  A3 candidate scan: best l=2 candidate k in %s at maxdev "
      "%.1e -- NOT a lock; measured values stand alone"
      % (str(best2[0]), best2[1]))
probe1 = [(om, ev) for om, ev in E_l1.items()]
best1 = None
for ks in ((1,), (2,), (3,), (1, 2)):
    dev = max(abs(cand(om, ks)/ev - 1) for om, ev in probe1)
    if best1 is None or dev < best1[1]:
        best1 = (ks, dev)
P("  A3 l=1 candidate: k in %s at maxdev %.1e%s"
  % (str(best1[0]), best1[1],
     " (LOCK)" if best1[1] < 1e-8 else " (no lock)"))
t2_ok = gates['GN-2a'] and gates['GN-2b']
if not t2_ok:
    E_l2 = {k: 1.0 for k in E_l2}
    P("  (T2 degraded: E_l2 := 1 carried into T3, letter cap N-GRAY)")
labels.append("T2 %s" % ("BANKED numeric-grade (l=0-exact-regressed, "
                         "2x-gated)" if t2_ok else "DEGRADED (E := 1)"))

# =====================================================================
# PART 2 -- T3: the ab-initio radiative number
# =====================================================================
P("")
P("PART 2 -- T3 THE AB-INITIO RADIATIVE CHANNEL")

mu_s, s_s, Om_s, ts = sp.symbols('mu s Omega t', positive=True)
x1 = s_s*sp.cos(Om_s*ts); y1 = s_s*sp.sin(Om_s*ts)
XY = [x1, y1, 0]
Q = sp.zeros(3, 3)
for i in range(3):
    for j in range(3):
        Q[i, j] = mu_s*(XY[i]*XY[j] - (s_s**2 if i == j else 0)/3)
Qddd = Q.applyfunc(lambda e: sp.diff(e, ts, 3))
ssum = sp.simplify(sum(Qddd[i, j]**2 for i in range(3)
                       for j in range(3)))
avg = sp.simplify(sp.integrate(ssum, (ts, 0, 2*sp.pi/Om_s))
                  / (2*sp.pi/Om_s))
P_gw = sp.simplify(avg/5)          # G = c = 1
ok_325 = sp.simplify(P_gw - sp.Rational(32, 5)*mu_s**2*s_s**4*Om_s**6) == 0
gates['GN-3a'] = ok_325
P("GN-3a circular-orbit tensor sum: <Qddd_ij Qddd_ij>/5 = %s == "
  "(32/5) mu^2 s^4 Omega^6 -> %s"
  % (sp.simplify(P_gw), "PASS" if ok_325 else "FAIL"))

Qbar = sp.zeros(3, 3)
for i in range(3):
    for j in range(3):
        Qbar[i, j] = sp.simplify(
            sp.integrate(Q[i, j], (ts, 0, 2*sp.pi/Om_s))
            / (2*sp.pi/Om_s))
dQ = Qbar.applyfunc(lambda e: sp.diff(e, s_s))
dQ2 = sp.simplify(sum(dQ[i, j]**2 for i in range(3) for j in range(3)))
ok_dq = sp.simplify(dQ2 - sp.Rational(2, 3)*mu_s**2*s_s**2) == 0
gates['GN-3b'] = ok_dq
P("GN-3b static-averaged charge: |dQbar/ds|^2 = %s == (2/3) mu^2 s^2 "
  "-> %s" % (sp.simplify(dQ2), "PASS" if ok_dq else "FAIL"))
if not (ok_325 and ok_dq):
    P("STOP: quadrupole/charge regression failed (bug-class).")
    save(); raise SystemExit(0)

# hbar audit (the cascade correspondence): an oscillator level n decays
# at n Gamma_spont; a coherent state radiates P_cl = Gamma hbar omega
# n_bar EXACTLY. Classical: P_cl = (G/5c^5) <Qddd^2> = (G/10c^5)
# omega^6 Lam_c^2 q0^2 with q0^2 = 2 hbar n_bar/(M_eff omega)
# => Gamma_spont = G omega^4 Lam_c^2/(5 c^5 M_eff)  [hbar CANCELS]
# => J_rad = Gamma/(2 pi) = G omega^4 Lam_c^2/(10 pi c^5 M_eff)
def J_rad_overH(mu_kg, s_m, om_tilde, Meff_kg, El2):
    om_si = om_tilde*H_SI
    Lam2 = (2.0/3.0)*mu_kg**2*s_m**2
    J_si = GSI*om_si**4*Lam2/(10*math.pi*C**5*Meff_kg)
    return (J_si/H_SI)*El2

AU = 1.496e11; KPC = 3.0857e19
ANCH = {
 'binary': dict(mu=0.5*MSUN, s=1e4*AU, Om=OM_BIN,
                El2=E_l2['binary']),
 'galaxy': dict(mu=1.0*MSUN, s=10*KPC, Om=OM_GAL,
                El2=E_l2['galaxy']),
}
P("")
P("  per-quantum radiative density at the window center (hbar "
  "cancels via the cascade correspondence; E_l2 included):")
Jrad = {}
for nm, a in ANCH.items():
    Jrad[nm] = J_rad_overH(a['mu'], a['s'], a['Om'], a['mu'],
                           a['El2'])
    P("    %-6s J_rad/H = %.3e  at M_eff = mu  (E_l2 = %.3f)"
      % (nm, Jrad[nm], a['El2']))

# =====================================================================
# PART 3 -- T4: the normalization bridge
# =====================================================================
P("")
P("PART 3 -- T4 THE NORMALIZATION BRIDGE")

GAP_GAL = X_GAL_AMB/TWO_PI
GAP_BIN = X_BIN/TWO_PI
G_GAL_GATE = 0.7536
G_BIN_GATE = 0.1117

def n_be(wv):
    z = TWO_PI*wv
    if z > 500: return 0.0
    return 1.0/math.expm1(z)

def Dhat(wv):                      # l=0 minimal shape, 9Z convention
    return (wv**2 + 1.0)/(4*math.pi**2*wv)

def window_integral(Om_, gap_, W_):
    lo = max(gap_, Om_ - W_)
    hi = Om_ + W_
    v, _ = quad(Dhat, lo, hi, limit=400)
    return v

def A_of_gc(gc, Om_, gap_, Wfac=1.0):
    W_ = Wfac*gc
    return gc*gc/window_integral(Om_, gap_, W_), W_

def delta_vac_up(Om_, W_, A_):
    return -(A_/(4*math.pi**2))*math.log(1 + Om_/W_)/Om_

def delta_vac_band(Om_, W_, gap_, A_):
    if Om_ - W_ <= gap_:
        return 0.0
    fb = lambda wv: math.log(wv/(Om_ - wv))
    return (A_/(4*math.pi**2))*(fb(Om_ - W_) - fb(gap_))/Om_

def delta_thermal(Om_, W_, gap_, A_, UV=30.0):
    def fth(wv):
        return Dhat(wv)*2*n_be(wv)/(Om_ - wv)
    tot = 0.0
    if Om_ - W_ > gap_:
        v1, _ = quad(fth, gap_, Om_ - W_, limit=400)
        tot += v1
    v2, _ = quad(fth, Om_ + W_, UV, limit=400)
    return (tot + v2)*A_

def Delta_gen(x_, A_, W_, gap_):
    Om_ = x_/TWO_PI
    return (delta_vac_up(Om_, W_, A_)
            + delta_vac_band(Om_, W_, gap_, A_)
            + delta_thermal(Om_, W_, gap_, A_))

def leak_gen(Om_, gap_, W_, A_):
    lo1, hi1 = gap_, Om_ - W_
    lo2 = Om_ + W_
    def fl(wv):
        d = Om_ - wv
        return Dhat(wv)*n_be(wv)*4*A_/(d*d)
    tot = 0.0
    if hi1 > lo1:
        v1, _ = quad(fl, lo1, hi1, limit=400)
        tot += v1
    UV = 6.0
    v2, _ = quad(fl, lo2, UV, limit=400)
    v2b, _ = quad(fl, UV, 2*UV, limit=400)
    return tot + v2 + v2b

def nu_f(x_):
    if x_ <= 0: return float('inf')
    if x_ > 500: return 1.0
    return 1.0 + 1.0/math.expm1(x_)

XLO, XHI = 0.14, 1.73
DEEP_LO, DEEP_HI = 0.14, 0.45
xg = np.linspace(XLO, XHI, 160)
Amat = np.vstack([np.ones_like(xg), xg]).T
mdeep = (xg >= DEEP_LO) & (xg <= DEEP_HI)

def budget(A_, W_):
    """returns (mean_deep, deep_end, split, coef, sat, max_resid).
    sat = True when nu(x + resid) leaves its domain anywhere in the
    deep window (residual exceeds x: past-FATAL-by-construction,
    amendment A2)."""
    dxv = np.array([TWO_PI*Delta_gen(float(x_), A_, W_, GAP_GAL)
                    for x_ in xg])
    coef, *_ = np.linalg.lstsq(Amat, dxv, rcond=None)
    resid = dxv - Amat @ coef
    sat = bool(np.any((xg + resid) <= 0))
    dnu = np.array([nu_f(float(x_) + float(r2)) - nu_f(float(x_))
                    for x_, r2 in zip(xg, resid)])
    with np.errstate(invalid='ignore'):
        mean_deep = float(np.mean(dnu[mdeep]))
        deep_end = float(dnu[0])
    dxb = TWO_PI*Delta_gen(X_BIN, A_, W_, GAP_BIN)
    rb = dxb - (coef[0] + coef[1]*X_BIN)
    return (mean_deep, deep_end, 2*abs(rb), coef, sat,
            float(np.max(np.abs(resid))))

m_r, de_r, split_r, _, _, _ = budget(1.0, 1.0)
ok4a = (abs(m_r - 0.004925) <= 2e-5 and abs(split_r - 0.0004) <= 2e-4)
gates['GN-4a'] = ok4a
P("GN-4a 9Z-b reader identity at (A=1, W=1): dc1_eff = %+.6f "
  "(archive +0.004925), split = %.5f (archive 0.0004) -> %s"
  % (m_r, split_r, "PASS" if ok4a else "FAIL"))
e_bin_r = leak_gen(OM_BIN, GAP_BIN, 1.0, 1.0)
e_gal_r = leak_gen(OM_GAL, GAP_GAL, 1.0, 1.0)
ok4b = (abs(e_bin_r/1.627e-5 - 1) <= 0.005
        and abs(e_gal_r/3.822e-5 - 1) <= 0.005)
gates['GN-4b'] = ok4b
P("GN-4b 9Z leak reader identity: eps_bin = %.3e (archive 1.627e-5), "
  "eps_gal = %.3e (archive 3.822e-5) -> %s"
  % (e_bin_r, e_gal_r, "PASS" if ok4b else "FAIL"))
if not (ok4a and ok4b):
    P("STOP: reader identity failed; nothing downstream quoted.")
    save(); raise SystemExit(0)

# ---- GN-5: the 9U dial, verbatim reimplementation -------------------
GATE_FID = 0.7536
def r_of_p(p, g_): return 2.0*(p - 0.5)/g_
def gamma_inst(r2): return 2*math.asin(math.sqrt(max(min(r2, 1.0),
                                                     0.0)))
def gamma_wind(r2):
    lo, hi = 1e-6, math.pi
    for _ in range(200):
        mid = 0.5*(lo + hi)
        if (1 - math.sin(mid)/mid)/2 < r2: lo = mid
        else: hi = mid
    return 0.5*(lo + hi)

r_hat = r_of_p(0.6471, GATE_FID)
gi = gamma_inst(r_hat); gw = gamma_wind(r_hat)
ok5 = (abs(r_of_p(0.65, GATE_FID) - 0.3980) < 5e-4
       and abs(gi - 1.35) <= 0.02 and abs(gw - 2.55) <= 0.03)
gates['GN-5'] = ok5
P("")
P("GN-5 dial regression: r(0.65) = %.4f (9U 0.3980); r_hat(0.6471) "
  "= %.4f -> gamma_inst = %.3f (9U 1.35), gamma_wind = %.3f "
  "(9U 2.55) -> %s"
  % (r_of_p(0.65, GATE_FID), r_hat, gi, gw,
     "PASS" if ok5 else "FAIL"))

SIG_R = 0.198
if ok5:
    r_lo = max(r_hat - SIG_R, 0.02); r_hi = min(r_hat + SIG_R, 0.98)
    gam_set = {'inst': (gamma_inst(r_lo), gi, gamma_inst(r_hi)),
               'wind': (gamma_wind(r_lo), gw, gamma_wind(r_hi))}
    band_label = "measured (9U r_hat +/- 1 sigma)"
else:
    gam_set = {'inst': (1.0, 1.35, 2.0), 'wind': (1.6, 2.55, 3.14)}
    band_label = "ENVELOPE-grade fallback"
T_CONV = {'t=1/H': 1.0, 't=2pi/H': TWO_PI}
centrals = {}
band_lo, band_hi = 1e9, 0.0
for rd, (glo, gc0, ghi) in gam_set.items():
    for tn, tv in T_CONV.items():
        centrals['%s %s' % (rd, tn)] = gc0/(2*tv)
        band_lo = min(band_lo, glo/(2*tv))
        band_hi = max(band_hi, ghi/(2*tv))
P("  g_c = gamma/(2 t_hat) [%s]: centrals %s" % (band_label,
  "  ".join("%s: %.3f" % (k, v) for k, v in centrals.items())))
P("  honest g_c band (union over readings x t-conventions x r_hat "
  "+/- 1 sigma): [%.3f, %.3f] H" % (band_lo, band_hi))

P("")
P("  dc1_T4(g_c) with A from the fixed point (W = g_c, galaxy "
  "geometry; the D-normalization cancels):")
P("  %8s %10s %12s %12s %10s %12s %12s" % ("g_c/H", "A(g_c)",
  "dc1_mean", "dc1_deepend", "split", "eps_g/own", "eps_g/minG"))
gc_grid = sorted(set([round(v, 6) for v in
                      list(np.geomspace(max(band_lo, 0.03),
                                        band_hi, 13))
                      + list(centrals.values())
                      + [0.05, 0.333, 2.0]]))
curve = {}
for gc in gc_grid:
    Ag, Wg = A_of_gc(gc, OM_GAL, GAP_GAL)
    m_, de_, sp_, _, sat_, _ = budget(Ag, Wg)
    ep_g = leak_gen(OM_GAL, GAP_GAL, Wg, Ag)
    curve[gc] = (Ag, m_, de_, sp_, ep_g, sat_)
    if sat_:
        P("  %8.3f %10.4f %12s %12s %10.5f %12.3e %12.3e"
          % (gc, Ag, "SAT(>FATAL)", "SAT(>FATAL)", sp_,
             ep_g/G_GAL_GATE, ep_g/G_BIN_GATE))
    else:
        P("  %8.3f %10.4f %+12.5f %+12.5f %10.5f %12.3e %12.3e"
          % (gc, Ag, m_, de_, sp_, ep_g/G_GAL_GATE, ep_g/G_BIN_GATE))
def _fatal_mean(gc):
    return curve[gc][5] or (math.isfinite(curve[gc][1])
                            and abs(curve[gc][1]) >= 0.15)
def _fatal_de(gc):
    return curve[gc][5] or (math.isfinite(curve[gc][2])
                            and abs(curve[gc][2]) >= 0.15)
kill = [gc for gc in gc_grid if _fatal_mean(gc)]
kill_de = [gc for gc in gc_grid if _fatal_de(gc)]
P("  kill-window (|mean| >= 0.15): %s"
  % ("{%.3f .. %.3f}" % (min(kill), max(kill)) if kill
     else "EMPTY on grid"))
P("  deep-end >= 0.15 window:     %s"
  % ("{%.3f .. %.3f}" % (min(kill_de), max(kill_de)) if kill_de
     else "EMPTY on grid"))
P("  four convention-central readings:")
cent_fatal = 0
for k, gc in centrals.items():
    Ag, Wg = A_of_gc(gc, OM_GAL, GAP_GAL)
    m_, de_, sp_, _, sat_, _ = budget(Ag, Wg)
    fat = sat_ or (math.isfinite(m_) and abs(m_) >= 0.15)
    cent_fatal += int(fat)
    P("    %-14s g_c = %.3f : dc1_mean = %s  deep-end = %s  %s"
      % (k, gc,
         "SAT" if sat_ else "%+.5f" % m_,
         "SAT" if sat_ else "%+.5f" % de_,
         "[FATAL-grade]" if fat else ""))
gc_gm = float(np.exp(np.mean(np.log(list(centrals.values())))))
P("  W-convention rows at g_c = %.3f (gm of centrals):" % gc_gm)
for wf in (0.5, 1.0, 2.0):
    Ag, Wg = A_of_gc(gc_gm, OM_GAL, GAP_GAL, Wfac=wf)
    m_, de_, sp_, _, sat_, _ = budget(Ag, Wg)
    P("    W = %.1f g_c : dc1_mean = %s  deep-end = %s"
      % (wf, "SAT" if sat_ else "%+.5f" % m_,
         "SAT" if sat_ else "%+.5f" % de_))
gates['GT5'] = bool(all(np.isfinite(curve[gc][0]) for gc in gc_grid))

# ---- U: the under-supply factor ------------------------------------
P("")
A_band = [curve[gc][0] for gc in gc_grid
          if band_lo - 1e-9 <= gc <= band_hi + 1e-9]
A_min_band = min(A_band) if A_band else min(c[0] for c in
                                            curve.values())
U = {}
for nm, a in ANCH.items():
    J_req_min = A_min_band*Dhat(a['Om'])
    U[nm] = Jrad[nm]/J_req_min
    P("  U(%s) = J_rad / J_req(weakest bridged) = %.3e / %.3e = %.3e"
      % (nm, Jrad[nm], J_req_min, U[nm]))
Mstar = {nm: a['mu']*U[nm] for nm, a in ANCH.items()}
P("  M_eff scaling: U ~ 1/M_eff; radiative supply would suffice only "
  "below M_eff* = %.2e kg (bin) / %.2e kg (gal) vs mu = %.2e / "
  "%.2e kg  -- the physical branch M_eff >= mu sits %.0f / %.0f "
  "orders past the crossover"
  % (Mstar['binary'], Mstar['galaxy'], ANCH['binary']['mu'],
     ANCH['galaxy']['mu'], -math.log10(U['binary']),
     -math.log10(U['galaxy'])))
P("  display row: U at M_eff/mu = 1e-10 would be %.2e (bin) / "
  "%.2e (gal)" % (U['binary']*1e10, U['galaxy']*1e10))

# ---- B-PHYS: the physical-normalization budget ----------------------
A_rad_gal = Jrad['galaxy']/Dhat(OM_GAL)
m_p, de_p, sp_p, _, sat_p, maxres_p = budget(A_rad_gal, 1.0)
ep_p_gal = leak_gen(OM_GAL, GAP_GAL, 1.0, A_rad_gal)
A_rad_bin = Jrad['binary']/Dhat(OM_BIN)
ep_p_bin = leak_gen(OM_BIN, GAP_BIN, 1.0, A_rad_bin)
# amendment A2: the direct nu-difference underflows double precision
# at A_rad ~ 1e-42..1e-35; quote the analytic bound |dc1| <=
# max|nu'| x max|resid| instead (nu' ~ -1/x^2 dominates at x = 0.14)
nup_max = max(abs(-math.exp(x_)/math.expm1(x_)**2) for x_ in
              (0.14, 0.20, 0.45))
dc1_bound = nup_max*maxres_p
P("")
P("  B-PHYS (proxy shape, physical normalization, W = 1 display "
  "geometry): |dc1_phys| <= %.2e (analytic bound max|nu'| x "
  "max|resid|; direct difference underflows double precision; PASS "
  "bar 0.05); eps_phys gal/bin = %.2e / %.2e (PASS bar 0.3 x 0.1117)"
  % (dc1_bound, ep_p_gal, ep_p_bin))
gates['B-PHYS'] = bool((not sat_p) and dc1_bound <= 0.005
                       and ep_p_gal <= 0.03*G_BIN_GATE
                       and ep_p_bin <= 0.03*G_BIN_GATE)
P("  l=2-quintic shape leg: EXCLUDED from verdict by pre-registration "
  "(secular approximation fails at omega ~ Omega_orb ~ 1e5 H; an "
  "all-radiative bridged budget under that shape is UV-cutoff-"
  "dominated -- itself evidence the all-radiative reading is not "
  "well-posed).")

# =====================================================================
# PART 4 -- the verdict (locked grammar)
# =====================================================================
P("")
u_forced = (U['binary'] <= 1e-6 and U['galaxy'] <= 1e-6)
u_live = (U['binary'] > 0.1 or U['galaxy'] > 0.1)
phys_fatal = (sat_p or dc1_bound >= 0.15 or ep_p_gal >= 3*G_BIN_GATE
              or ep_p_bin >= 3*G_BIN_GATE)
if u_live and cent_fatal == 4:
    P("==> 10A VERDICT (locked grammar): N-FATAL -- the radiative "
      "sector could carry the exchange AND the bridged budget is "
      "FATAL at every convention central. Pre-signed strike: "
      "bath-mechanism conditional 15 -> 8.")
elif phys_fatal:
    P("==> 10A VERDICT (locked grammar): N-FATAL -- the physical "
      "budget itself breaches. Pre-signed strike: 15 -> 8.")
elif u_forced and gates['B-PHYS'] and t2_ok:
    P("==> 10A VERDICT (locked grammar): N-SPLIT-CLOSED -- the "
      "D-provenance DISCHARGES.%s"
      % ("" if ok5 else " [g_c band at ENVELOPE grade]"))
    P("    (1) The radiative graviton sector CANNOT supply the "
      "measured exchange: per-quantum under-supply U = %.1e (bin) / "
      "%.1e (gal) on the physical branch M_eff >= mu; the crossover "
      "masses M_eff* = %.1e / %.1e kg sit %.0f / %.0f orders below "
      "any coordinate that moves the orbit's mass."
      % (U['binary'], U['galaxy'], Mstar['binary'], Mstar['galaxy'],
         -math.log10(U['binary']), -math.log10(U['galaxy'])))
    P("    (2) The O(1) grammar therefore rides the CONSTRAINT "
      "sector, whose credential is exact: SdS first law dM = -T_c "
      "dS_c (all orders), Birkhoff steps, and a gapped multipole "
      "tower costing only %.1e x g_c^2 in admixture." % tower)
    P("    (3) The PHYSICAL dispersive budget and leak are the "
      "radiative-sector ones: |dc1_phys| <= %.1e (analytic bound), "
      "eps_phys <= %.1e -- inside every band by tens of orders. The "
      "9Z-b joint pessimal corner evaporates: it priced a radiative "
      "normalization the physics does not supply."
      % (dc1_bound, max(ep_p_gal, ep_p_bin)))
    P("    (4) The all-radiative CONDITIONAL is quoted above with "
      "its own gamma kill-window (D cancels under the bridge; that "
      "budget rides the measured gamma alone). The round-21 one-seam "
      "objection resolves by sector separation: one field carries an "
      "O(1) constraint coupling and an eta^4-per-quantum radiative "
      "coupling at once -- the cavity-QED geometry.")
else:
    P("==> 10A VERDICT (locked grammar): N-GRAY -- quote everything; "
      "cells: U_forced=%s U_live=%s B-PHYS=%s T2=%s dial=%s "
      "cent_fatal=%d/4." % (u_forced, u_live, gates['B-PHYS'],
                            t2_ok, ok5, cent_fatal))
P("    CREDENCE: per the pre-signed map (N-FATAL -> 15 to 8; "
  "otherwise HOLD 15 pending ROUND 22; rise cell 15 -> 18 requires "
  "R22 no-unpatched-hole AND N-SPLIT-CLOSED).")

need = {'mech-9z-leak': False, 'mech-9zb-diff': False,
        'mech-9w-multimode': False, 'mech-9t-resonance': False,
        'gal-9u-ptail': False}
for row in csv.reader(open('LEDGER.csv', encoding='utf-8')):
    if row and row[0] in need and row[1] in ('CURRENT', 'CO-QUOTED'):
        need[row[0]] = True
gates['GN-6'] = all(need.values())
P("")
P("GN-6 ledger legs: %s -> %s"
  % (need, "PASS" if gates['GN-6'] else "FAIL"))
P("")
P("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
  for k, v in sorted(gates.items())))
P("LABELS: " + " | ".join(labels))
P("")
P("done (%.1f min)" % ((time.time() - t00)/60))
save()
print("\nsaved: data/stage10a_dprov.txt")
