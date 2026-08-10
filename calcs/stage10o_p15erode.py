"""
STAGE 10O -- P15-ERODE: the eccentricity-erosion ledger
========================================================
Pre-registered 2026-08-11 (this file committed BEFORE any run; bars, gates,
letters, credence map locked below). The P15-stage queued by TODO 32(a) and
the DIARY 2026-08-11 same-day update ("pre-reg bars blind; gates: reproduce
M&H eq. 26 as G0, Hwang alpha(s) table primary re-read, statistics bar set
before the integral runs").

QUESTION
--------
P15 (PREDICTIONS.md section E, registered 61767d4) leg (ii): the frozen
ambient-cloud l=2 fluctuation adds a random tide that ERODES wide-binary
superthermality; what does the SURVIVAL of the measured alpha(s)
(Hwang-Ting-Zakamska 2022, Table 1, primary-read 2026-08-11) bound?

DESIGN-LEVEL SHARPENING (disclosed pre-run, before any computation)
-------------------------------------------------------------------
1. ERRATUM to the registered parenthetical: P15 leg (ii) printed
   "delta g / g_bin ~ s^2". The interior-regular l=2 solution (10E: r^2 Y_2)
   gives a UNIFORM tidal curvature Gamma_c inside the cloud: delta g =
   Gamma_c * s, so delta g / g_bin = Gamma_c s^3 / (G m_b) ~ s^3. The
   registered exponent was imprecise; direction and instrument unchanged.
   Disclosed here before any number is computed.
2. THE ONCE-ONLY LEMMA: phase mixing in a FIXED field is idempotent (the
   phase-mixed DF is invariant under further evolution in the same field;
   M&H 2023 sec 2). A strictly frozen cloud (tau_c ~ 1/H ~ age) therefore
   applies the erosion map AT MOST ONCE, and a once-only halving is largely
   REABSORBABLE into the unknown primordial alpha_0(s). The registered
   survival-bound hope is therefore sharpened to THREE instruments:
   (a) THE WEDGE: a joint exclusion in (Gamma_c, tau_c) -- amplitude x
       refresh time. A refreshing cloud applies k_c ~ age/tau_c mixings;
       the measured alpha(s) bins cap the total reachable interval after
       k mixings (superthermal side: primordial ceiling alpha_0 <= 2.5;
       subthermal side: the integrability floor alpha_0 > -1 makes the cap
       PRIOR-FREE).
   (b) THE STEP: a cloud with in-window mixing boundary (t_mix,c = age
       crossing 100 AU - 31.6 kAU) predicts a localized DOWNWARD step in
       alpha(s); the measured smooth run tests it.
   (c) ALPHA_FLOOR: one full mixing maps ANY initial DF into
       [alpha_floor, map(alpha_0max)]; if alpha_floor exceeds a measured
       bin's alpha + 2 sigma, even a FROZEN cloud that mixes that bin is
       excluded -> a genuine frozen-corner bound. alpha_floor is MEASURED
       in Part A (unknown at pre-reg time; the bars below do not depend
       on its value).
3. FIREWALL: this stage touches NO program binary data and NO program
   likelihood machinery. Inputs = Hwang+22 Table 1 (transcribed below from
   the 2026-08-11 primary read) + literature constants (H22/M&H23). The
   anomaly-real credence CANNOT move here (no sky statistic of ours is
   read). Our w_rad = 0.20 and the 4G cross-validation are cited as
   context only, never used in a gate.

INPUTS (primary-read 2026-08-11, ar5iv)
---------------------------------------
- M&H 2023 (2303.15531): eq. 26 f_inf = contour average of f_0 weighted by
  residence time; eq. 30 delta H = (A a^2 / 8) H_Gamma with
  H_Gamma = (1/j^2)[(j^2 - 3 Gamma jz^2)(5 - 3 j^2)
                    - 15 Gamma (j^2 - jz^2)(1 - j^2) cos 2w];
  Gamma = 1/3 for the Galactic disc tide (H22 sec 4); eq. 33
  alpha_f ~ (1 + alpha_i)/2 near unity; eq. 18 P(e) = (1+alpha) e^alpha;
  t0 ~ 4 Gyr (rho0/0.2 Msun pc^-3)^-1 (m_b/Msun)^(1/2) (a/10^4 AU)^(-3/2).
- H22 (2202.01307): t_sec = T_Z^2/T_b with (2pi/T_Z)^2 = 4 pi G rho0,
  rho0 = 0.2 Msun/pc^3 fiducial (Widmark 2019); phase-mixed after ~5 t_sec;
  outer orbit fixed (no re-randomization modeled).
- Hwang+22 (2111.01789) Table 1 (verbatim transcription, symmetric sigma =
  max of the +/- errors):
    log10 s bin      N       alpha    sigma
    [1.50, 2.00]     686     0.08     0.15
    [2.00, 2.25]     2814    0.59     0.10
    [2.25, 2.50]     12357   0.82     0.05
    [2.50, 2.75]     18379   0.94     0.04
    [2.75, 3.00]     14884   1.20     0.05
    [3.00, 3.50]     18496   1.30     0.06
    [3.50, 4.00]     7570    1.32     0.09
    [4.00, 4.50]     2301    1.17     0.15
  Convention identical to M&H eq. 18. Projected s taken as a (disclosed;
  bins are 0.25-0.5 dex wide, the s-vs-a scatter is subdominant).

PARTS
-----
A. THE MAP MACHINERY (the DIARY's G0): sympy-derived secular EOM from
   H_Gamma; batch RK4 ensemble mixing at Gamma in {1/3, 1}; measured map
   table alpha_f(alpha_i) on alpha_i in {-0.99, 0.1, 0.5, 0.8, 1.0, 1.2,
   1.6, 2.0, 2.5}; alpha_floor := alpha_f(-0.99); refresh chains (random
   re-orientation between mixings) at DF level for k = 1, 2, 3.
B. TIMESCALES: t0(a) implemented and gated against H22/M&H normalizations;
   Galactic mixing boundary; cloud t_mix,c(a) = t0_Gal(a) * Gamma_Gal /
   Gamma_c; conversion Gamma_c = (e_a kappa_c) * a0 / r_M (r_M =
   sqrt(G m_b / a0); kappa_c = the unpinned O5-NORM conversion constant --
   all bounds quoted on the PRODUCT e_a * kappa_c, P10 comparisons at
   kappa_c = 1 labeled conditional).
C. THE WEDGE: grid tau_c in [0.05, 10] Gyr x Gamma_c/Gamma_Gal in
   [1e-2, 1e6]; per-bin reachable intervals after k_tot mixings (DF-level
   chains for k <= 3, interpolated measured map beyond); excluded cell iff
   any bin's [alpha_j +/- 2 sigma_j] misses the reachable interval under
   the rules below; frozen-corner readout at tau_c = age.
D. THE STEP: quadratic-in-log-s null vs null + downward step, location
   scanned; verdict per the pre-fixed bars.
E. CONSISTENCY READ (no novelty claim): the measured Gamma=1/3 map applied
   across the Galactic mixing boundary vs the observed bin-7/bin-8 values,
   at rho0 = 0.2 and 0.1.

WEDGE RULES (locked)
--------------------
- k_c(tau_c, Gamma_c; a_j) = floor(age/tau_c) if [t_mix,c(a_j) <= tau_c
  AND Gamma_c >= 2 Gamma_Gal] else 0 (a refresh only re-randomizes the
  combined field if the cloud is the dominant term; factor 2 = disclosed
  threshold, robustness row at 1 and 4).
- k_gal(a_j) = 1 if t0_Gal(a_j) <= age else 0 (partial Galactic mixing of
  narrow bins ignored = CONSERVATIVE for exclusion: any partial Galactic
  erosion only shrinks the cloud's room).
- k_tot = k_c + k_gal; reachable interval [A_k, B_k] with A_k =
  map^(k-1)(alpha_floor) for k >= 1 (A_0 = -1), B_k = map^k(2.5)
  (B_0 = 2.5). alpha_0 ceiling 2.5 = disclosed prior axis (rows also at
  2.0 and 1.6); the subthermal side needs NO prior (alpha_0 > -1 is
  integrability).
- Exclusion at 2 sigma per bin; PRIMARY quote uses bins 2-8 (bin 1
  [32-100 AU] = the most selection-exposed bin; its strengthened numbers
  reported as a labeled SECONDARY row).
- age = 10 Gyr (rows at 8 and 12); m_b = 1 Msun (row at 2); rho0 = 0.2
  (row at 0.1).

BARS (all fixed before any computation; the step bars were fixed in the
session record BEFORE the Hwang table was fetched)
--------------------------------------------------
G10O-0a  measured map vs eq. 33 near unity: |alpha_f - (1+alpha_i)/2|
         <= 0.05 at alpha_i in {0.8, 1.2}, both Gamma.
G10O-0b  thermal fixed point: |alpha_f(1.0) - 1| <= 0.02, both Gamma.
G10O-0c  idempotence: second-window pooled e-histogram vs first-window,
         L1 <= 0.02 (mixedness).
G10O-0d  H_Gamma conservation: 99th pct |Delta H|/(|H|+1) <= 2e-3
         (convergence protocol: halve dtau and re-run once if violated;
         a second violation = G-FAIL).
G10O-0f  vZLK exact turning point (Gamma = 1): the e_max(i_0) relation
         j_min^2 = (5/3) jz^2 for circular-start librators reproduced on
         contours to |Delta H| <= 1e-10 (symbolic check via sympy).
G10O-1   transcription: the table above reproduces the paper's headline
         ("superthermal at > 10^3 AU", alpha = 1.32 +/- 0.09 at
         3.16-10 kAU) -- checked by inspection against the primary read;
         any mismatch = G-FAIL.
G10O-2   timescale: implemented t0(10^4 AU, 1 Msun, rho0 = 0.2) within
         x2.5 of BOTH the M&H 4 Gyr and the H22 5 x 1 Gyr mixing figures.
G10O-3   resolution: map values at alpha_i in {-0.99, 1.0, 2.5} stable to
         <= 0.03 under (2x ensemble, 2x tau_max, seed 101).
G10O-4   wedge positive clause (trap #12): at the 3 exclusion-boundary
         cells adjacent to tau_c = {0.3, 1, 3} Gyr, a forward simulation
         (alpha_0(s) = observed alpha(s), cloud erosion applied by the
         rules) must violate >= 1 bin at 2 sigma -- i.e. the claimed
         exclusion actually excludes. Fail = wedge quoted CANDIDATE-only.
G10O-5a  step power: a -0.15 step injected at x0 = 2.75 into the smooth
         null best-fit with the measured sigmas is recovered at
         Delta chi2 >= 9 in >= 80% of 200 draws (seed 31). Fail = step
         verdict POWER-LIMITED regardless of data Delta chi2.
G10O-5b  step null calibration: stepless draws give Delta chi2 >= 9 in
         <= 10% of 200 draws. Fail = the SUPPORTED letter unavailable.
STEP VERDICT BARS: SUPPORTED iff data Delta chi2 >= 9 (and 5b holds);
ABSENT (band exclusion quotable) iff data Delta chi2 <= 4 (and 5a holds);
else UNRESOLVED.

LETTERS (grammar wired to gates, trap #11)
------------------------------------------
O-WEDGE            core gates (0a-0f, 1, 2, 3) PASS and the wedge is
                   non-empty in the physical domain: quote the wedge, the
                   frozen-corner row, the step verdict, the consistency
                   read. (G4 fail -> wedge CANDIDATE-only; step sub-letter
                   per its own bars.)
O-DEGENERATE       core gates PASS but the wedge is EMPTY everywhere:
                   the instrument is power-dead at public-catalog grade;
                   quote the once-only degeneracy statement + alpha_floor.
O-UNCONVERGED      G10O-3 fails: no number quoted; re-run queued.
O-FEASIBILITY-LIMITED  G10O-0a fails after the convergence protocol, or
                   G10O-1 fails: stop at the named blocker (9R/10H
                   discipline).
NOTE: leg (i) of P15 (the no-pump prediction) is NOT tested here and
cannot fire in either direction: this data cannot evidence a pump
(no formation ceiling exists that the measured values exceed -- checked
mechanically in Part C output); kill (a) remains DR4-era as registered.

PRE-SIGNED CREDENCE MAP (all cells)
-----------------------------------
EVERY letter above -> HOLD anomaly-real 53 / mech-conditional 8. This
stage is a constraint instrument + theorem formalization; no sky statistic
of the program is read; the frozen-corner outcome (if it fires) is a
CONSISTENCY dividend for the frozen reading, booked as annotation, not a
credence move. If the wedge bites the P10 e_a window, P10/PREDICTIONS get
annotations (no status flips from this stage).

Environment: py (Windows), numpy + sympy + scipy. Seed 31 primary, 101
resolution row. Wall-clock estimate ~5-15 min CPU.
"""

import json
import time

import numpy as np
import sympy as sp

T_START = time.time()
RNG = np.random.default_rng(31)

# ----------------------------------------------------------------------
# Constants (SI)
# ----------------------------------------------------------------------
G_SI = 6.674e-11
MSUN = 1.989e30
AU = 1.496e11
PC = 3.086e16
YR = 3.156e7
GYR = 1e9 * YR
A0 = 1.2e-10

RHO0_FID = 0.2 * MSUN / PC**3      # H22/M&H fiducial
RHO0_ALT = 0.1 * MSUN / PC**3
AGE_FID = 10.0                      # Gyr
AGES = [8.0, 10.0, 12.0]
MB_FID = 1.0                        # Msun
ALPHA0_CEILS = [1.6, 2.0, 2.5]
ALPHA0_CEIL_PRIMARY = 2.5

GAMMA_GAL = 4.0 * np.pi * G_SI * RHO0_FID     # s^-2, dominant zz component

# Hwang+22 Table 1 (see header)
BIN_EDGES = np.array([1.50, 2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50])
ALPHA_OBS = np.array([0.08, 0.59, 0.82, 0.94, 1.20, 1.30, 1.32, 1.17])
SIG_OBS = np.array([0.15, 0.10, 0.05, 0.04, 0.05, 0.06, 0.09, 0.15])
N_OBS = np.array([686, 2814, 12357, 18379, 14884, 18496, 7570, 2301])
BIN_MID = 0.5 * (BIN_EDGES[:-1] + BIN_EDGES[1:])          # log10 AU
A_J = 10.0**BIN_MID * AU                                   # separation ~ a (m)
NBINS = len(ALPHA_OBS)

GATES = {}
RESULTS = {}


def gate(name, passed, detail):
    GATES[name] = (bool(passed), detail)
    print(f"  GATE {name}: {'PASS' if passed else 'FAIL'}  ({detail})")


# ----------------------------------------------------------------------
# Part A: the map machinery
# ----------------------------------------------------------------------
print("=" * 72)
print("PART A: the phase-mixing map (M&H eq. 26/30/33 machinery)")
print("=" * 72)

j_s, jz_s, w_s, Gam_s = sp.symbols('j j_z omega Gamma', positive=False)
H_expr = (1 / j_s**2) * ((j_s**2 - 3 * Gam_s * jz_s**2) * (5 - 3 * j_s**2)
                         - 15 * Gam_s * (j_s**2 - jz_s**2) * (1 - j_s**2)
                         * sp.cos(2 * w_s))
dH_dj = sp.diff(H_expr, j_s)
dH_dw = sp.diff(H_expr, w_s)
H_fn = sp.lambdify((j_s, w_s, jz_s, Gam_s), H_expr, 'numpy')
dHdj_fn = sp.lambdify((j_s, w_s, jz_s, Gam_s), dH_dj, 'numpy')
dHdw_fn = sp.lambdify((j_s, w_s, jz_s, Gam_s), dH_dw, 'numpy')

# G10O-0f: exact vZLK turning point at Gamma = 1 (symbolic)
i0 = sp.symbols('i0', positive=True)
jz_c = sp.cos(i0)
H_circ = H_expr.subs([(j_s, 1), (jz_s, jz_c), (Gam_s, 1)])
j_turn = sp.sqrt(sp.Rational(5, 3)) * jz_c
H_turn = H_expr.subs([(j_s, j_turn), (w_s, sp.pi / 2), (jz_s, jz_c),
                      (Gam_s, 1)])
resid_0f = sp.simplify(H_circ - H_turn)
gate('G10O-0f', resid_0f == 0, f"symbolic vZLK turning-point residual = {resid_0f}")


def sample_initial(alpha_i, n, rng):
    """Isotropic ensemble with P(e) = (1+alpha) e^alpha."""
    u = rng.random(n)
    e = u**(1.0 / (1.0 + alpha_i))
    e = np.clip(e, 1e-4, 1 - 1e-6)
    cosi = rng.uniform(-1, 1, n)
    w = rng.uniform(0, 2 * np.pi, n)
    j = np.sqrt(1 - e**2)
    jz = j * cosi
    return j, w, jz


def mix_ensemble(j0, w0, jz, Gam, tau_max=80.0, dtau=1e-3, t_frac=0.30,
                 sample_every=100):
    """Batch RK4 of dj/dtau = -dH/dw, dw/dtau = +dH/dj.

    Returns pooled e-samples over [t_frac*tau_max, tau_max] in two half
    windows (for the idempotence gate), the final state, and the H drift.
    """
    j = j0.copy()
    w = w0.copy()
    H_init = H_fn(j, w, jz, Gam)
    nstep = int(round(tau_max / dtau))
    istart = int(nstep * t_frac)
    imid = istart + (nstep - istart) // 2
    jlo = np.abs(jz) + 1e-9
    pool1, pool2 = [], []

    def rhs(jj, ww):
        jj = np.clip(jj, jlo, 1.0 - 1e-12)
        return -dHdw_fn(jj, ww, jz, Gam), dHdj_fn(jj, ww, jz, Gam)

    for istep in range(nstep):
        k1j, k1w = rhs(j, w)
        k2j, k2w = rhs(j + 0.5 * dtau * k1j, w + 0.5 * dtau * k1w)
        k3j, k3w = rhs(j + 0.5 * dtau * k2j, w + 0.5 * dtau * k2w)
        k4j, k4w = rhs(j + dtau * k3j, w + dtau * k3w)
        j = j + (dtau / 6) * (k1j + 2 * k2j + 2 * k3j + k4j)
        w = w + (dtau / 6) * (k1w + 2 * k2w + 2 * k3w + k4w)
        j = np.clip(j, jlo, 1.0 - 1e-12)
        if istep >= istart and (istep - istart) % sample_every == 0:
            e_now = np.sqrt(np.clip(1 - j**2, 0, 1))
            (pool1 if istep < imid else pool2).append(e_now.copy())
    H_fin = H_fn(j, w, jz, Gam)
    drift = np.abs(H_fin - H_init) / (np.abs(H_init) + 1.0)
    return (np.concatenate(pool1), np.concatenate(pool2),
            (j, w, jz), float(np.percentile(drift, 99)))


def fit_alpha(e_samples):
    """MLE for P(e) = (1+alpha) e^alpha."""
    e = e_samples[(e_samples > 1e-6) & (e_samples < 1 - 1e-9)]
    return -len(e) / np.sum(np.log(e)) - 1.0


ALPHA_GRID = [-0.99, 0.1, 0.5, 0.8, 1.0, 1.2, 1.6, 2.0, 2.5]
N_ENS = 1500
MAP_TABLE = {}
drift_worst = 0.0
idem_worst = 0.0

for Gam in (1.0 / 3.0, 1.0):
    for a_i in ALPHA_GRID:
        j0, w0, jz = sample_initial(a_i, N_ENS, RNG)
        p1, p2, _, drift = mix_ensemble(j0, w0, jz, Gam)
        a_f1, a_f2 = fit_alpha(p1), fit_alpha(p2)
        h1, _ = np.histogram(p1, bins=25, range=(0, 1), density=True)
        h2, _ = np.histogram(p2, bins=25, range=(0, 1), density=True)
        l1 = np.mean(np.abs(h1 - h2)) / 2.0
        MAP_TABLE[(Gam, a_i)] = 0.5 * (a_f1 + a_f2)
        drift_worst = max(drift_worst, drift)
        idem_worst = max(idem_worst, l1)
        print(f"  Gamma={Gam:.3f} alpha_i={a_i:+.2f} -> alpha_f={MAP_TABLE[(Gam, a_i)]:+.4f}"
              f"  (halves w1/w2 {a_f1:+.3f}/{a_f2:+.3f}, drift99 {drift:.1e}, L1 {l1:.3f})")

gate('G10O-0d', drift_worst <= 2e-3, f"worst 99th-pct H drift = {drift_worst:.2e}")
gate('G10O-0c', idem_worst <= 0.02, f"worst window L1 = {idem_worst:.3f}")

err_0a = max(abs(MAP_TABLE[(g, a)] - 0.5 * (1 + a))
             for g in (1.0 / 3.0, 1.0) for a in (0.8, 1.2))
gate('G10O-0a', err_0a <= 0.05, f"max |map - (1+a)/2| near unity = {err_0a:.4f}")
err_0b = max(abs(MAP_TABLE[(g, 1.0)] - 1.0) for g in (1.0 / 3.0, 1.0))
gate('G10O-0b', err_0b <= 0.02, f"max |map(1) - 1| = {err_0b:.4f}")

ALPHA_FLOOR = {g: MAP_TABLE[(g, -0.99)] for g in (1.0 / 3.0, 1.0)}
print(f"  ALPHA_FLOOR: Gamma=1/3 -> {ALPHA_FLOOR[1/3]:.4f}, Gamma=1 -> {ALPHA_FLOOR[1.0]:.4f}")
RESULTS['map_table'] = {f"{g:.3f},{a:+.2f}": v for (g, a), v in MAP_TABLE.items()}
RESULTS['alpha_floor'] = {f"{g:.3f}": v for g, v in ALPHA_FLOOR.items()}

# Refresh chains at DF level (k = 1, 2, 3), Gamma = 1 (cloud geometry)
print("Refresh chains (random re-orientation between mixings), Gamma = 1:")


def elements_to_vectors(j, w, jz, rng):
    n = len(j)
    cosi = np.clip(jz / np.maximum(j, 1e-12), -1, 1)
    sini = np.sqrt(np.clip(1 - cosi**2, 0, 1))
    Om = rng.uniform(0, 2 * np.pi, n)   # node uniform (axisymmetric mixed state)
    e = np.sqrt(np.clip(1 - j**2, 0, 1))
    # unit angular momentum vector
    jhat = np.stack([sini * np.sin(Om), -sini * np.cos(Om), cosi], axis=1)
    # eccentricity vector: direction of periapsis in the orbit plane
    xhat = np.stack([np.cos(Om), np.sin(Om), np.zeros(n)], axis=1)
    yhat = np.cross(jhat, xhat)
    ehat = xhat * np.cos(w)[:, None] + yhat * np.sin(w)[:, None]
    return e * ehat[:, 0], e * ehat[:, 1], e * ehat[:, 2], j, jhat


def vectors_to_elements(evec_x, evec_y, evec_z, j, jhat, axis):
    """Re-express elements w.r.t. a new field axis."""
    e = np.sqrt(evec_x**2 + evec_y**2 + evec_z**2)
    cosi = jhat @ axis
    jz = j * cosi
    # omega: angle from the ascending node to periapsis, in the orbit plane
    node = np.cross(np.tile(axis, (len(j), 1)), jhat)
    nn = np.linalg.norm(node, axis=1)
    polar = nn < 1e-8
    node[polar] = np.array([1.0, 0, 0])
    nn[polar] = 1.0
    node /= nn[:, None]
    evec = np.stack([evec_x, evec_y, evec_z], axis=1)
    ehat = np.where(e[:, None] > 1e-12, evec / np.maximum(e, 1e-12)[:, None], node)
    cosw = np.clip(np.sum(node * ehat, axis=1), -1, 1)
    sign = np.sum(np.cross(node, ehat) * jhat, axis=1)
    w = np.where(sign >= 0, np.arccos(cosw), 2 * np.pi - np.arccos(cosw))
    return j, w, jz


CHAIN_ALPHAS = {}
for a_i in (-0.99, 1.32, 2.5):
    j, w, jz = sample_initial(a_i, N_ENS, RNG)
    chain = []
    for k in range(1, 4):
        p1, p2, (j, w, jz), _ = mix_ensemble(j, w, jz, 1.0)
        chain.append(fit_alpha(np.concatenate([p1, p2])))
        # refresh: random new axis
        ex, ey, ez, jj, jhat = elements_to_vectors(j, w, jz, RNG)
        v = RNG.normal(size=3)
        axis = v / np.linalg.norm(v)
        j, w, jz = vectors_to_elements(ex, ey, ez, jj, jhat, axis)
    CHAIN_ALPHAS[a_i] = chain
    print(f"  chain from alpha_0={a_i:+.2f}: k=1,2,3 -> "
          + ", ".join(f"{c:+.4f}" for c in chain))
RESULTS['chains'] = {f"{k:+.2f}": v for k, v in CHAIN_ALPHAS.items()}

# G10O-1: transcription check (mechanical restatement of the primary read)
ok_1 = (abs(ALPHA_OBS[6] - 1.32) < 1e-9 and abs(SIG_OBS[6] - 0.09) < 1e-9
        and ALPHA_OBS[5] > 1.0 and ALPHA_OBS[4] > 1.0)
gate('G10O-1', ok_1, "Table 1 headline (superthermal > 10^3 AU; 1.32+/-0.09 at 3.16-10 kAU)")

# ----------------------------------------------------------------------
# Part B: timescales
# ----------------------------------------------------------------------
print("=" * 72)
print("PART B: timescales and conversion")
print("=" * 72)


def t0_gyr(a_m, gamma=GAMMA_GAL, mb=MB_FID):
    """Mixing time t0 = 2 pi n / Gamma (M&H normalization), in Gyr."""
    n_mm = np.sqrt(G_SI * mb * MSUN / a_m**3)
    return 2 * np.pi * n_mm / gamma / GYR


t0_ref = t0_gyr(1e4 * AU)
ok_2 = (t0_ref / 4.0 < 2.5 and 4.0 / t0_ref < 2.5
        and t0_ref / 5.0 < 2.5 and 5.0 / t0_ref < 2.5)
gate('G10O-2', ok_2, f"t0(10^4 AU) = {t0_ref:.2f} Gyr vs M&H 4 / H22 5x1")

R_M = np.sqrt(G_SI * MB_FID * MSUN / A0)
CONV = A0 / R_M / GAMMA_GAL     # Gamma_c/Gamma_Gal per unit (e_a kappa_c)
print(f"  r_M(1 Msun) = {R_M/AU:.0f} AU; Gamma_Gal = {GAMMA_GAL:.3e} s^-2")
print(f"  conversion: Gamma_c/Gamma_Gal = {CONV:.3e} x (e_a kappa_c)")
print(f"  Galactic t0 per bin (Gyr): "
      + ", ".join(f"{t0_gyr(a):.1f}" for a in A_J))
RESULTS['t0_ref_gyr'] = t0_ref
RESULTS['conv_ratio_per_product'] = CONV

# ----------------------------------------------------------------------
# Part C: the wedge
# ----------------------------------------------------------------------
print("=" * 72)
print("PART C: the exclusion wedge")
print("=" * 72)

alpha_i_arr = np.array(ALPHA_GRID)


def map_interp(gam, a):
    vals = np.array([MAP_TABLE[(gam, x)] for x in ALPHA_GRID])
    return float(np.interp(a, alpha_i_arr, vals))


def reachable(k, gam, ceil):
    """[A_k, B_k] after k full mixings, primordial alpha_0 in (-1, ceil]."""
    if k == 0:
        return -1.0, ceil
    lo = ALPHA_FLOOR[gam]
    hi = map_interp(gam, ceil)
    for _ in range(k - 1):
        lo = map_interp(gam, lo)
        hi = map_interp(gam, hi)
    return lo, hi


def wedge(tau_grid, ratio_grid, ceil=ALPHA0_CEIL_PRIMARY, age=AGE_FID,
          use_bins=range(1, NBINS), nsig=2.0, dom_factor=2.0, gam=1.0):
    """Excluded[i, j] over (tau_c, Gamma_c/Gamma_Gal)."""
    exc = np.zeros((len(tau_grid), len(ratio_grid)), dtype=bool)
    kgal = np.array([1 if t0_gyr(a) <= age else 0 for a in A_J])
    for i, tau in enumerate(tau_grid):
        kc_base = int(np.floor(age / tau))
        for jx, ratio in enumerate(ratio_grid):
            bad = False
            for b in use_bins:
                tmix_c = t0_gyr(A_J[b]) / ratio
                kc = kc_base if (tmix_c <= tau and ratio >= dom_factor) else 0
                ktot = kc + kgal[b]
                lo, hi = reachable(ktot, gam, ceil)
                if (ALPHA_OBS[b] + nsig * SIG_OBS[b] < lo
                        or ALPHA_OBS[b] - nsig * SIG_OBS[b] > hi):
                    bad = True
                    break
            exc[i, jx] = bad
    return exc, kgal


TAU_GRID = np.geomspace(0.05, 10.0, 25)
RATIO_GRID = np.geomspace(1e-2, 1e6, 33)
EXC, KGAL = wedge(TAU_GRID, RATIO_GRID)
print(f"  k_gal per bin: {KGAL.tolist()}")

wedge_summary = {}
for tau_pick in (0.1, 0.3, 1.0, 3.0, AGE_FID):
    i = int(np.argmin(np.abs(TAU_GRID - tau_pick)))
    row = EXC[i]
    if row.any():
        lo_r = RATIO_GRID[row][0]
        prod = lo_r / CONV
        wedge_summary[tau_pick] = (float(lo_r), float(prod))
        print(f"  tau_c ~ {TAU_GRID[i]:.2f} Gyr: excluded ratio >= {lo_r:.3g}"
              f"  -> e_a*kappa_c >= {prod:.3g}")
    else:
        wedge_summary[tau_pick] = None
        print(f"  tau_c ~ {TAU_GRID[i]:.2f} Gyr: NO exclusion")
RESULTS['wedge_primary'] = {str(k): v for k, v in wedge_summary.items()}

# frozen corner: tau_c = age (one configuration)
i_frozen = len(TAU_GRID) - 1
frozen_row = EXC[i_frozen]
if frozen_row.any():
    lo_r = RATIO_GRID[frozen_row][0]
    print(f"  FROZEN CORNER (tau_c = age): excluded ratio >= {lo_r:.3g}"
          f" -> e_a*kappa_c >= {lo_r/CONV:.3g}")
    RESULTS['frozen_corner'] = float(lo_r / CONV)
else:
    print("  FROZEN CORNER (tau_c = age): NO exclusion (alpha_0-degenerate)")
    RESULTS['frozen_corner'] = None

# secondary row: bin 1 included
EXC_B1, _ = wedge(TAU_GRID, RATIO_GRID, use_bins=range(0, NBINS))
sec = {}
for tau_pick in (1.0, AGE_FID):
    i = int(np.argmin(np.abs(TAU_GRID - tau_pick)))
    row = EXC_B1[i]
    sec[tau_pick] = float(RATIO_GRID[row][0] / CONV) if row.any() else None
print(f"  SECONDARY (bin 1 included): tau=1 -> {sec[1.0]}, frozen -> {sec[AGE_FID]}")
RESULTS['wedge_bin1'] = {str(k): v for k, v in sec.items()}

# robustness rows
rob = {}
for label, kw in [
        ('ceil2.0', dict(ceil=2.0)), ('ceil1.6', dict(ceil=1.6)),
        ('age8', dict(age=8.0)), ('age12', dict(age=12.0)),
        ('dom1', dict(dom_factor=1.0)), ('dom4', dict(dom_factor=4.0)),
        ('gamthird', dict(gam=1.0 / 3.0))]:
    e2, _ = wedge(TAU_GRID, RATIO_GRID, **kw)
    i1 = int(np.argmin(np.abs(TAU_GRID - 1.0)))
    r = e2[i1]
    rob[label] = float(RATIO_GRID[r][0] / CONV) if r.any() else None
print("  robustness (tau_c = 1 Gyr bound on e_a*kappa_c): "
      + ", ".join(f"{k}={v if v is None else format(v, '.3g')}" for k, v in rob.items()))
RESULTS['wedge_robustness_tau1'] = rob

# no-pump non-test check (mechanical): can the data evidence a pump here?
# A pump would need alpha_obs above the reachable ceiling for k=0 (= the
# prior ceiling itself). Max observed + 2 sigma vs ceiling:
max_obs = float(np.max(ALPHA_OBS + 2 * SIG_OBS))
print(f"  no-pump non-test: max(alpha + 2 sigma) = {max_obs:.2f} < ceiling"
      f" {ALPHA0_CEIL_PRIMARY} -> pump direction CANNOT fire here (as pre-stated)")
RESULTS['nopump_nontest'] = max_obs < ALPHA0_CEIL_PRIMARY

# G10O-4: injection at three boundary cells
print("G10O-4 injection at the exclusion boundary:")
inj_ok = []
for tau_pick in (0.3, 1.0, 3.0):
    i = int(np.argmin(np.abs(TAU_GRID - tau_pick)))
    row = EXC[i]
    if not row.any():
        continue
    jx = int(np.argmax(row))            # first excluded ratio
    ratio = RATIO_GRID[jx]
    tau = TAU_GRID[i]
    kc_base = int(np.floor(AGE_FID / tau))
    violated = False
    for b in range(1, NBINS):
        tmix_c = t0_gyr(A_J[b]) / ratio
        kc = kc_base if (tmix_c <= tau and ratio >= 2.0) else 0
        ktot = kc + KGAL[b]
        if ktot == 0:
            continue
        a_now = ALPHA_OBS[b]
        for _ in range(ktot):
            a_now = map_interp(1.0, max(a_now, -0.99))
        if abs(a_now - ALPHA_OBS[b]) > 2 * SIG_OBS[b]:
            violated = True
    inj_ok.append(violated)
    print(f"  cell (tau={tau:.2f}, ratio={ratio:.3g}): forward-sim violates a bin: {violated}")
gate('G10O-4', len(inj_ok) > 0 and all(inj_ok),
     f"boundary cells verified: {inj_ok}")

# G10O-3: resolution row (2x ensemble, 2x tau_max, seed 101)
print("G10O-3 resolution row:")
rng2 = np.random.default_rng(101)
res_err = 0.0
for a_i in (-0.99, 1.0, 2.5):
    j0, w0, jz = sample_initial(a_i, 2 * N_ENS, rng2)
    p1, p2, _, _ = mix_ensemble(j0, w0, jz, 1.0, tau_max=160.0)
    a_hi = fit_alpha(np.concatenate([p1, p2]))
    res_err = max(res_err, abs(a_hi - MAP_TABLE[(1.0, a_i)]))
    print(f"  alpha_i={a_i:+.2f}: hi-res map {a_hi:+.4f} vs {MAP_TABLE[(1.0, a_i)]:+.4f}")
gate('G10O-3', res_err <= 0.03, f"max map shift = {res_err:.4f}")

# ----------------------------------------------------------------------
# Part D: the step statistic
# ----------------------------------------------------------------------
print("=" * 72)
print("PART D: the step statistic")
print("=" * 72)

X = BIN_MID
W_INV = 1.0 / SIG_OBS**2


def fit_null(y, x=X, wi=W_INV):
    A = np.vstack([np.ones_like(x), x, x**2]).T
    Aw = A * wi[:, None]
    coef, *_ = np.linalg.lstsq(Aw.T @ A, Aw.T @ y, rcond=None)
    resid = y - A @ coef
    return float(np.sum(wi * resid**2)), coef


def fit_step(y, x=X, wi=W_INV, x0_scan=None):
    if x0_scan is None:
        x0_scan = np.arange(2.375, 3.751, 0.125)
    best = (np.inf, None, None)
    for x0 in x0_scan:
        step_col = (x >= x0).astype(float)
        if step_col.sum() in (0, len(x)):
            continue
        A = np.vstack([np.ones_like(x), x, x**2, step_col]).T
        Aw = A * wi[:, None]
        coef, *_ = np.linalg.lstsq(Aw.T @ A, Aw.T @ y, rcond=None)
        if coef[3] > 0:      # downward steps only (pre-registered)
            continue
        resid = y - A @ coef
        c2 = float(np.sum(wi * resid**2))
        if c2 < best[0]:
            best = (c2, x0, coef[3])
    return best


chi2_null, null_coef = fit_null(ALPHA_OBS)
chi2_step, x0_best, step_best = fit_step(ALPHA_OBS)
DCHI2 = chi2_null - chi2_step if np.isfinite(chi2_step) else 0.0
print(f"  data: chi2_null = {chi2_null:.2f} (5 dof), best step at x0 = {x0_best},"
      f" size = {step_best if step_best is None else format(step_best, '.3f')},"
      f" Delta chi2 = {DCHI2:.2f}")
RESULTS['step_data'] = dict(chi2_null=chi2_null, dchi2=float(DCHI2),
                            x0=x0_best, size=step_best)

# G10O-5a/5b calibrations
A_null = np.vstack([np.ones_like(X), X, X**2]).T
y_smooth = A_null @ null_coef
n_draw = 200
rec, fp = 0, 0
rng5 = np.random.default_rng(31)
for _ in range(n_draw):
    noise = rng5.normal(0, SIG_OBS)
    y_inj = y_smooth + noise - 0.15 * (X >= 2.75)
    d1 = fit_null(y_inj)[0] - fit_step(y_inj)[0]
    if d1 >= 9:
        rec += 1
    y_nul = y_smooth + rng5.normal(0, SIG_OBS)
    d0 = fit_null(y_nul)[0] - fit_step(y_nul)[0]
    if d0 >= 9:
        fp += 1
gate('G10O-5a', rec / n_draw >= 0.80, f"step recovery {rec}/{n_draw}")
gate('G10O-5b', fp / n_draw <= 0.10, f"false positive {fp}/{n_draw}")

if not GATES['G10O-5a'][0]:
    step_verdict = 'POWER-LIMITED'
elif DCHI2 >= 9 and GATES['G10O-5b'][0]:
    step_verdict = 'SUPPORTED'
elif DCHI2 <= 4:
    step_verdict = 'ABSENT'
else:
    step_verdict = 'UNRESOLVED'
print(f"  STEP VERDICT: {step_verdict}")
RESULTS['step_verdict'] = step_verdict

# If ABSENT: the in-window band exclusion (theorem-size step would be seen)
if step_verdict == 'ABSENT':
    # band: clouds whose mixing boundary falls in-window at tau_c >= age
    # (static): t_mix,c(a) = age at a in [10^2.375, 10^3.75] AU
    lo_a, hi_a = 10**2.375 * AU, 10**3.75 * AU
    r_hi = t0_gyr(lo_a) / AGE_FID
    r_lo = t0_gyr(hi_a) / AGE_FID
    print(f"  in-window static band excluded: Gamma_c/Gamma_Gal in"
          f" [{r_lo:.3g}, {r_hi:.3g}] -> e_a*kappa_c in"
          f" [{r_lo/CONV:.3g}, {r_hi/CONV:.3g}]")
    RESULTS['step_band'] = [float(r_lo / CONV), float(r_hi / CONV)]

# ----------------------------------------------------------------------
# Part E: consistency read (no novelty claim)
# ----------------------------------------------------------------------
print("=" * 72)
print("PART E: Galactic-halving consistency read")
print("=" * 72)
for rho_lab, rho in (('0.2', RHO0_FID), ('0.1', RHO0_ALT)):
    gam_g = 4 * np.pi * G_SI * rho
    mixed = [t0_gyr(a, gamma=gam_g) <= AGE_FID for a in A_J]
    pred = map_interp(1.0 / 3.0, 1.32)
    d8 = (ALPHA_OBS[7] - pred) / SIG_OBS[7]
    print(f"  rho0={rho_lab}: bins mixed within age: "
          f"{[i+1 for i, m in enumerate(mixed) if m]};"
          f" map(1.32) = {pred:.3f} vs bin-8 1.17+/-0.15 -> d = {d8:+.2f} sigma")
RESULTS['consistency_map132'] = map_interp(1.0 / 3.0, 1.32)

# ----------------------------------------------------------------------
# Letter
# ----------------------------------------------------------------------
print("=" * 72)
# AMENDMENT A1 (2026-08-11, logged pre-quote, before any run output was
# read): the original letter block routed a 0b/0c/0d/0f core-gate failure
# to O-DEGENERATE, whose grammar means "gates PASS, wedge empty" -- a
# trap-#11 wiring gap. A1: ANY post-protocol core-gate failure other than
# the 0a/1 feasibility pair and the 3 convergence gate -> the instrument
# fails at O-FEASIBILITY-LIMITED grade with the failing gate named. The
# stage is seeded/deterministic; when all core gates pass, A1 and the
# original wiring produce the identical letter.
core = ['G10O-0a', 'G10O-0b', 'G10O-0c', 'G10O-0d', 'G10O-0f',
        'G10O-1', 'G10O-2', 'G10O-3']
core_ok = all(GATES[g][0] for g in core)
wedge_nonempty = EXC.any()
if not GATES['G10O-0a'][0] or not GATES['G10O-1'][0]:
    LETTER = 'O-FEASIBILITY-LIMITED'
elif not GATES['G10O-3'][0]:
    LETTER = 'O-UNCONVERGED'
elif not core_ok:
    failed = [g for g in core if not GATES[g][0]]
    LETTER = f'O-FEASIBILITY-LIMITED (gate {",".join(failed)})'
elif wedge_nonempty:
    LETTER = 'O-WEDGE' + ('' if GATES['G10O-4'][0] else ' (CANDIDATE-only)')
else:
    LETTER = 'O-DEGENERATE'
n_pass = sum(1 for v in GATES.values() if v[0])
print(f"LETTER: {LETTER}   gates {n_pass}/{len(GATES)}"
      f"   step: {step_verdict}   wall {time.time()-T_START:.0f} s")
RESULTS['letter'] = LETTER
RESULTS['gates'] = {k: [v[0], v[1]] for k, v in GATES.items()}

with open('data/stage10o_results.json', 'w') as f:
    json.dump(RESULTS, f, indent=1)
print("results -> data/stage10o_results.json")
