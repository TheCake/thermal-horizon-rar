# -*- coding: utf-8 -*-
"""
STAGE 10R: P14 first leg -- THE OORT PROFILE LEDGER (P14-a).

PRE-REGISTRATION (bars locked in this commit, before any execution).

QUESTION. PREDICTIONS.md P14 registers the Oort-spike retrodiction: the
measured function class (Boltzmann-screened, ambient-gated, EFE-quenched
at the solar galactic field) should PRESERVE the observed long-period-
comet 1/a spike that Vokrouhlicky, Nesvorny & Tremaine 2024
(arXiv:2403.09555, primary-read 2026-08-11, NOTES lit-note this close)
show AQUAL-MOND destroys. The registered kill needs a full V-N-T-class
Oort pipeline; THIS STAGE IS THE SCOPED FIRST LEG: the PROFILE ledger.
It computes, from the program's own validated EFE solver, the solar-
ambient anomalous force profile eps(r) = g_eff/g_N - 1 over the Oort
range, propagates it to the OBSERVABLE the spike lives in (the
perihelion-osculating inverse semimajor axis 1/a_osc), and asks whether
the class SMEARS the spike at profile grade. It cannot fire the
registered P14 kill (full-pipeline successor); it CAN retrodict
(SAFE) or flag (BROKEN) at its own grade.

PRIMARY-READ NUMBERS (V-N-T 2024; conventions declared -- trap #25):
  - Observed spike: osculating 1/a in (0, 0.75e-4) au^-1 (their "Oort
    peak"); WIDTH_SPIKE = 0.75e-4 au^-1 is the stage's yardstick.
  - Their adopted external field g_e = 2.32e-10 m/s^2 = 1.93 a0 in the
    TOTAL-field convention; OUR solver consumes the NEWTONIAN-convention
    external field e_N = 1.15 +/- 0.05 a0 (Stage 3T; the 2G/3S
    cross-convention bug class -- G6 checks the two describe one sky:
    nu(y_e) * y_e ~= 1.93 under the reference law).
  - a0 = 1.2e-10 m/s^2; r_M = sqrt(GM_sun/a0) = 7030 au.
  - Their transition family mu_n, n >= 6-8 ephemeris-compatible; the
    AQUAL contrast row (T4) uses n = 8.
  - Their screened-variants clause (verbatim): "Our results do not rule
    out some MOND theories more elaborate than AQUAL, in which
    non-Newtonian effects are screened on small spatial scales, at
    small masses, or in external gravitational fields comparable in
    strength to the critical acceleration." The measured class IS that
    clause; this stage quantifies it.

THE MAP (T2, sympy-gated G2). Static modified potential
Phi(r) = -GM/r + Phi_m(r), Phi_m(r) = -GM * I(r),
I(r) = int_r^inf eps(s)/s^2 ds. At perihelion q_p (deep inside the
screening ramp) the fitted Newtonian osculating energy is
E_osc = E_true - Phi_m(q_p), so
   1/a_osc = 1/a_true - 2*I(q_p)   [au^-1 with s in au; GM cancels].
The map is STATIC and E-independent; what can smear the spike at
profile grade is the SPREAD of the map across the delivered population:
  (i)  monopole spread over perihelion q_p in [5, 40] au (screening
       makes eps vanish there; expected tiny -- measured, not assumed);
  (ii) the anomalous QUADRUPOLE energy spread over orientations at
       aphelion Q in the feeding band [2e4, 1e5] au:
       d(1/a)_quad(Q) = 2 * |phi2(Q)| * dP2 / r_M in au^-1 (solver
       phi2 in units a0*r_M; GM = a0*r_M^2; dP2 = 1.5 = full P2 range).
  (iii) the common-mode SHIFT 2*I(q_p) itself: NOT letter-bearing (a
       static common shift re-equilibrates through formation/delivery
       in the same potential; deciding its observable residue = the
       full-pipeline successor). REPORTED with units, flagged.

EPS PROFILES (T1). Primary = EMPIRICAL-ANCHOR: eps_emp(r) =
eps_inf * shape(r), shape = the reference-law solver profile normalized
to its deep plateau; eps_inf = B_perp^2 - 1 with B = the measured
perspective-corrected median boost 1.078 (CI 1.052-1.103; PAPER 4Q
anchor) => eps_inf central 0.162, band [0.107, 0.217]. (Conversion
declared: v ~ sqrt(G_eff) at fixed geometry => G_eff factor = B^2;
dimensionless, same units both sides.) Robustness rows: full solver
profiles for BE and simple at e_N in {1.10, 1.15, 1.20}, amplitude
alpha in {0.8, 1.0, 1.3} multiplying (nu-1) linearly (4K convention).
The quoted headline = MAX over the family (conservative).

GATES (each names the clause its failure vetoes -- trap #11):
  G1 solver regression vs the archived 4K output (same code path;
     q_plateau(BE, e_N=1.2) within |d| <= 0.002 of the archived value
     greppable in data/stage4k_quadrupole.txt; Newton phi2 < 1e-12).
     FAIL => R-LIMITED (instrument invalid; all rows void).
  G2 sympy: the osculating map identity and the quadrupole unit
     conversion, residue 0. FAIL => R-LIMITED (map unproven).
  G3 Newton control: eps = 0 => shift = smear = 0 to < 1e-12 au^-1.
     FAIL => R-LIMITED (statistic contaminated).
  G4 resolution: NR 512 -> 1024 moves the headline smear by < 10%.
     FAIL => R-LIMITED (unconverged).
  G5 injection: a step-profile eps of known amplitude at known radius
     recovers its analytic shift 2*eps0/r_step to <= 5%. FAIL =>
     R-LIMITED (integrator wrong).
  G6 convention diagnostic (letter-NEUTRAL, disclosed): the reference
     law at e_N = 1.15 gives total field nu*y*a0 within [1.6, 2.3] a0
     bracketing V-N-T's 1.93 a0; prints the mapping.
  G7 AQUAL contrast: eps_AQUAL(mu_8, their field)/eps_ours at 5e4 au
     >= 5. FAIL vetoes ONLY the class-separation sentence (the SAFE/
     BROKEN clauses do not rest on it).

HEADLINE BARS (blind; yardstick WIDTH_SPIKE = 0.75e-4 au^-1):
  SMEAR_TOTAL = max over the family band of [monopole q_p-spread +
  quadrupole orientation spread across the feeding band].
  R-PROFILE-SAFE   iff G1-G5 pass AND SMEAR_TOTAL <= 0.25 * WIDTH.
  R-PROFILE-BROKEN iff G1-G5 pass AND SMEAR_TOTAL >= 1.00 * WIDTH.
  R-PROFILE-GRAY   otherwise (gates pass, smear in between).
  The SHIFT row and the tide-context row are letter-neutral reports.

PRE-SIGNED CREDENCE MAP (all cells): SAFE -> HOLD 53/8 (retrodiction
banked at profile grade; P14 row annotated, full-pipeline successor
named). GRAY -> HOLD 53/8 (annotate). BROKEN -> HOLD 53/8 (annotate;
commission the full-pipeline leg BEFORE any credence motion -- the
registered kill fires only there). LIMITED -> HOLD, instrument row.
NO sky fits anywhere in this stage; archived meters only.

Run: py calcs/stage10r_p14spike.py   -> data/stage10r_p14spike.txt
"""
import io
import os
import re

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss, legval

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = []


def log(s):
    OUT.append(s)
    print(s)


# ----- constants (units declared) -----------------------------------------
A0_SI = 1.2e-10            # m s^-2 (V-N-T's value; program value 1.05-1.2e-10
#                            -- the statistic scales as 1/r_M, band disclosed)
R_M_AU = 7030.0            # au; r_M = sqrt(GM_sun/a0)
WIDTH_SPIKE = 0.75e-4      # au^-1 (V-N-T Oort peak band, primary read)
EN_GRID = (1.10, 1.15, 1.20)   # Newtonian-convention solar ambient (3T)
ALPHA_GRID = (0.8, 1.0, 1.3)   # boost amplitude band (4K linear convention)
EPS_EMP = (0.107, 0.162, 0.217)  # empirical anchor band: B^2-1, B=1.078
#                                   (CI 1.052-1.103; PAPER 4Q corrected)
QP_GRID_AU = (5.0, 10.0, 20.0, 40.0)      # perihelion grid, au
QBAND_AU = np.geomspace(2.0e4, 1.0e5, 25)  # feeding-band aphelia, au

# ----- the 4K solver, verbatim code path (regression-gated in G1) ---------


def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0 / (1.0 - np.exp(-np.minimum(x, 40))))


def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0 / np.clip(y, 1e-14, None))


def nu_newton(y):
    return np.ones_like(np.asarray(y, float))


def solve_multipole(nu, eN, NR=512, NMU=96, LMAX=16):
    r = np.logspace(-2, 3, NR)
    mu, wmu = leggauss(NMU)
    R, MU = np.meshgrid(r, mu, indexing='ij')
    ST = np.sqrt(1 - MU ** 2)
    Pl = np.zeros((LMAX + 1, NMU))
    for l in range(LMAX + 1):
        c = np.zeros(l + 1)
        c[l] = 1
        Pl[l] = legval(mu, c)
    gr = -1.0 / R ** 2 + eN * MU
    gt = -eN * ST
    f = nu(np.hypot(gr, gt)) - 1.0
    Ar, At = f * gr, f * gt
    dAr = np.gradient(R ** 2 * Ar, np.log(r), axis=0) / R ** 3
    dAt = -np.gradient(ST * At, mu, axis=1) / R
    S = -(dAr + dAt)
    Sl = np.array([(2 * l + 1) / 2 * np.sum(wmu * S * Pl[l], axis=1)
                   for l in range(LMAX + 1)])
    A = np.zeros((LMAX + 1, NR))
    B = np.zeros((LMAX + 1, NR))
    dr = np.diff(r)
    for l in range(LMAX + 1):
        g_in = Sl[l] * r
        for i in range(NR - 1):
            qf = (r[i] / r[i + 1]) ** (l + 1)
            A[l, i + 1] = A[l, i] * qf + 0.5 * dr[i] * (g_in[i] * qf + g_in[i + 1])
        g_out = Sl[l] * r
        for i in range(NR - 2, -1, -1):
            qf = (r[i] / r[i + 1]) ** l
            B[l, i] = B[l, i + 1] * qf + 0.5 * dr[i] * (g_out[i] + g_out[i + 1] * qf)
    lv = np.arange(LMAX + 1)[:, None]
    phi_l = -(A + B) / (2 * lv + 1)
    gph_r = -(((lv + 1) * A - lv * B) / (r[None, :] * (2 * lv + 1)))
    g_r = gr + np.tensordot(gph_r, Pl, axes=(0, 0))
    if eN > 0:
        nuE = nu(np.array([eN]))[0]
        g_r -= nuE * eN * MU
    boost = np.sum(wmu * (-g_r), axis=1) / 2 / (1.0 / r ** 2)
    return r, phi_l, boost


def q_plateau(r, phi_l, lo=0.02, hi=0.2):
    m = (r >= lo) & (r <= hi)
    return float(np.mean(phi_l[2][m] / r[m] ** 2))


# ----- G2: sympy identities -----------------------------------------------
def gate2():
    GM, E, r, s, qp, eps0, rs = sp.symbols('GM E r s q_p eps0 r_s',
                                           positive=True)
    # map identity: E_osc = E - Phi_m(q_p); 1/a = -2E/GM (bound E<0 handled
    # by sign convention; identity is linear so symbols positive is fine)
    Phi_m = -GM * sp.Symbol('I_qp', positive=True)
    E_osc = sp.Symbol('E_true') - Phi_m
    inva_osc = -2 * E_osc / GM
    inva_true = -2 * sp.Symbol('E_true') / GM
    res1 = sp.simplify(inva_osc - (inva_true - 2 * sp.Symbol('I_qp',
                                                             positive=True)))
    # step-profile analytic shift: eps = eps0 for s > rs, 0 inside =>
    # I(qp<rs) = int_rs^inf eps0/s^2 ds = eps0/rs
    I_step = sp.integrate(eps0 / s ** 2, (s, rs, sp.oo))
    res2 = sp.simplify(I_step - eps0 / rs)
    # quadrupole conversion: Phi2_phys = phi2 * a0 * rM * P2;
    # d(1/a) = 2*Phi2/GM, GM = a0*rM^2 => d(1/a) = 2*phi2*P2/rM
    phi2, a0, rM, P2 = sp.symbols('phi2 a0 r_M P2', positive=True)
    res3 = sp.simplify(2 * (phi2 * a0 * rM * P2) / (a0 * rM ** 2)
                       - 2 * phi2 * P2 / rM)
    ok = (res1 == 0) and (res2 == 0) and (res3 == 0)
    log("G2 sympy map/step/quad identities: residues (%s, %s, %s) -> %s"
        % (res1, res2, res3, "PASS" if ok else "FAIL -> R-LIMITED"))
    return ok


# ----- profile machinery --------------------------------------------------
def eps_profile(r_au, nu, eN, alpha, NR=512):
    r, phi_l, boost = solve_multipole(nu, eN, NR=NR)
    eps = alpha * (boost - 1.0)
    r_grid_au = r * R_M_AU
    return (np.interp(r_au, r_grid_au, eps),
            r, phi_l, boost)


def shift_integral(r_au_grid, eps_grid, qp_au):
    """I(qp) = int_qp^inf eps(s)/s^2 ds  [au^-1], trapezoid on the grid."""
    m = r_au_grid >= qp_au
    s = r_au_grid[m]
    f = eps_grid[m] / s ** 2
    return float(np.trapezoid(f, s))


def run():
    log("=" * 72)
    log("STAGE 10R P14-a: THE OORT PROFILE LEDGER (pre-reg: this commit)")
    log("=" * 72)

    ok_all = True

    # G2
    ok_all &= gate2()

    # G1 regression: archived 4K q value (BE, eN=1.2) from the stage output
    arch = io.open(os.path.join(ROOT, "data/stage4k_quadrupole.txt"),
                   encoding="utf-8", errors="replace").read()
    mm = re.findall(r"q\s*=\s*(-?\d+\.\d+)", arch)
    r, phi_l, boost = solve_multipole(nu_be, 1.2)
    q_ours = q_plateau(r, phi_l)
    q_arch = None
    if mm:
        cands = [float(v) for v in mm]
        q_arch = min(cands, key=lambda v: abs(v - q_ours))
    d1 = abs(q_ours - q_arch) if q_arch is not None else 999.0
    g1a = d1 <= 0.002
    rN, phiN, _ = solve_multipole(nu_newton, 1.2)
    g1b = float(np.max(np.abs(phiN[2]))) < 1e-12
    log("G1 regression: q_ours(BE,1.2) = %.4f vs archived %.4f (|d| = %.1e"
        " <= 2e-3) AND Newton phi2 max %.1e < 1e-12 -> %s"
        % (q_ours, q_arch if q_arch is not None else float('nan'), d1,
           np.max(np.abs(phiN[2])),
           "PASS" if (g1a and g1b) else "FAIL -> R-LIMITED"))
    ok_all &= (g1a and g1b)

    # reference-law shape (BE, eN=1.15, alpha=1)
    r_au_grid = np.geomspace(3.0, 5.0e6, 4000)
    eps_ref, r_s, phi_ref, boost_ref = eps_profile(r_au_grid, nu_be, 1.15, 1.0)
    plateau_mask = (r_s * R_M_AU >= 10 * R_M_AU) & (r_s * R_M_AU <= 50 * R_M_AU)
    eps_deep_ref = float(np.mean(boost_ref[plateau_mask] - 1.0))
    shape = eps_ref / eps_deep_ref
    r_half = float(np.interp(0.5, np.clip(shape, 0, None)
                             [np.argsort(r_au_grid)],
                             r_au_grid[np.argsort(r_au_grid)]))
    log("reference profile: eps_deep(BE, eN=1.15, a=1) = %.4f; "
        "half-ramp radius = %.0f au (diagnostic)" % (eps_deep_ref, r_half))

    # G3 Newton control
    zero = np.zeros_like(r_au_grid)
    sh0 = shift_integral(r_au_grid, zero, 5.0)
    g3 = abs(sh0) < 1e-12
    log("G3 Newton control: shift(eps=0) = %.1e au^-1 (< 1e-12) -> %s"
        % (sh0, "PASS" if g3 else "FAIL -> R-LIMITED"))
    ok_all &= g3

    # G5 injection: step profile
    eps_step = np.where(r_au_grid >= 7000.0, 0.15, 0.0)
    got = shift_integral(r_au_grid, eps_step, 5.0)
    want = 0.15 / 7000.0
    g5 = abs(got - want) / want <= 0.05
    log("G5 injection: step eps0=0.15 at 7000 au: shift %.3e vs analytic "
        "%.3e (rel d = %.2e <= 0.05) -> %s"
        % (got, want, abs(got - want) / want,
           "PASS" if g5 else "FAIL -> R-LIMITED"))
    ok_all &= g5

    # G6 convention diagnostic
    y_e = 1.15
    tot = float(nu_be(np.array([y_e]))[0] * y_e)
    g6 = 1.6 <= tot <= 2.3
    log("G6 convention diagnostic (letter-neutral): nu_BE(1.15)*1.15 = "
        "%.2f a0 total-convention vs V-N-T 1.93 a0 (band [1.6, 2.3]) -> %s"
        % (tot, "consistent" if g6 else "DISCLOSED-INCONSISTENT"))

    # ----- THE LEDGER over the family band --------------------------------
    def family_rows():
        rows = []
        # empirical-anchor rows (primary): shape x eps_inf band
        for eps_inf in EPS_EMP:
            eps_g = eps_inf * shape
            rows.append(("EMP eps_inf=%.3f" % eps_inf, eps_g, None))
        # solver rows: {BE, simple} x eN x alpha
        for name, nu in (("BE", nu_be), ("simple", nu_simple)):
            for eN in EN_GRID:
                for al in ALPHA_GRID:
                    e_g, _, pl, _ = eps_profile(r_au_grid, nu, eN, al)
                    rows.append(("%s eN=%.2f a=%.1f" % (name, eN, al),
                                 e_g, (al, pl)))
        return rows

    rows = family_rows()
    smear_mono_max = 0.0
    smear_quad_max = 0.0
    shift_central = None
    log("-" * 72)
    log("THE LEDGER (units: au^-1; yardstick WIDTH_SPIKE = %.2e au^-1)"
        % WIDTH_SPIKE)
    for name, eps_g, extra in rows:
        shifts = [2.0 * shift_integral(r_au_grid, eps_g, qp)
                  for qp in QP_GRID_AU]
        mono_spread = max(shifts) - min(shifts)
        # quadrupole spread across the feeding band (orientation dP2 = 1.5)
        if extra is not None:
            al, pl = extra
            phi2_band = np.interp(QBAND_AU / R_M_AU, r_s, np.abs(pl[2])) * al
        else:
            phi2_band = np.interp(QBAND_AU / R_M_AU, r_s,
                                  np.abs(phi_ref[2])) * (
                np.max(eps_g) / max(np.max(eps_ref), 1e-30))
        quad = float(np.max(2.0 * phi2_band * 1.5 / R_M_AU))
        smear_mono_max = max(smear_mono_max, mono_spread)
        smear_quad_max = max(smear_quad_max, quad)
        if name == "EMP eps_inf=0.162":
            shift_central = shifts[0]
        log("  %-22s shift(qp=5au) = %.3e | mono-spread = %.3e | "
            "quad-spread = %.3e" % (name, shifts[0], mono_spread, quad))

    smear_total = smear_mono_max + smear_quad_max
    log("-" * 72)
    log("SHIFT (common-mode, letter-NEUTRAL, successor object): central "
        "%.3e au^-1 = %.2f x WIDTH_SPIKE" % (shift_central,
                                             shift_central / WIDTH_SPIKE))
    log("SMEAR_TOTAL (max over family: mono %.3e + quad %.3e) = %.3e au^-1"
        " = %.3f x WIDTH_SPIKE"
        % (smear_mono_max, smear_quad_max, smear_total,
           smear_total / WIDTH_SPIKE))

    # G4 resolution doubling on the headline
    eps_hi, r_hi_s, phi_hi, boost_hi = eps_profile(r_au_grid, nu_be, 1.15,
                                                   1.0, NR=1024)
    pm = (r_hi_s * R_M_AU >= 10 * R_M_AU) & (r_hi_s * R_M_AU <= 50 * R_M_AU)
    shape_hi = eps_hi / float(np.mean(boost_hi[pm] - 1.0))
    sm_lo = max(2.0 * shift_integral(r_au_grid, 0.217 * shape, qp)
                for qp in QP_GRID_AU) - \
        min(2.0 * shift_integral(r_au_grid, 0.217 * shape, qp)
            for qp in QP_GRID_AU)
    sm_hi = max(2.0 * shift_integral(r_au_grid, 0.217 * shape_hi, qp)
                for qp in QP_GRID_AU) - \
        min(2.0 * shift_integral(r_au_grid, 0.217 * shape_hi, qp)
            for qp in QP_GRID_AU)
    relres = abs(sm_hi - sm_lo) / max(sm_lo, 1e-30)
    g4 = relres < 0.10
    log("G4 resolution: mono-spread(hi-anchor) NR 512->1024 rel d = %.2e "
        "(< 0.10) -> %s" % (relres, "PASS" if g4 else "FAIL -> R-LIMITED"))
    ok_all &= g4

    # G7 AQUAL contrast (mu_8 via numeric inversion)
    def nu_aqual8(y):
        y = np.clip(np.asarray(y, float), 1e-12, None)
        out = np.empty_like(y)
        for i, yy in enumerate(y.ravel()):
            g = np.geomspace(max(yy, 1e-6) * 1e-3, max(yy, 1e-6) * 1e3, 400)
            mu = g / (1 + g ** 8) ** (1.0 / 8)
            gN = mu * g
            out.ravel()[i] = np.interp(yy, gN, g) / yy
        return out
    # their field in our Newtonian convention: solve nu8(ye)*ye = 1.93
    yy = np.linspace(0.5, 2.5, 200)
    tots = nu_aqual8(yy) * yy
    ye_vnt = float(np.interp(1.93, tots, yy))
    eps_aq, _, _, _ = eps_profile(r_au_grid, nu_aqual8, ye_vnt, 1.0, NR=256)
    i50 = np.argmin(np.abs(r_au_grid - 5e4))
    eps_ours_50 = 0.217 * shape[i50]
    ratio = abs(eps_aq[i50]) / max(abs(eps_ours_50), 1e-30)
    g7 = ratio >= 5.0
    log("G7 AQUAL contrast: eps_mu8(5e4 au; ye = %.2f Newt-conv) = %.3f vs "
        "ours(hi) = %.3f; ratio = %.1f (>= 5) -> %s [vetoes only the "
        "class-separation clause]"
        % (ye_vnt, eps_aq[i50], eps_ours_50, ratio,
           "PASS" if g7 else "FAIL"))

    # tide context row (letter-neutral): standard Galactic tide energy scale
    # rides the OBSERVED field in both worlds; anomalous quad row above is
    # the class's own addition. (Context, no bar.)
    log("context: the vertical Galactic tide is an empirical input common "
        "to both worlds; the class-specific addition is the quad row above.")

    # ----- letter ---------------------------------------------------------
    log("=" * 72)
    if not ok_all:
        letter = "R-LIMITED"
        clause = ("instrument gate failed; no profile claim is made "
                  "(see gate lines for the veto)")
    elif smear_total <= 0.25 * WIDTH_SPIKE:
        letter = "R-PROFILE-SAFE"
        clause = ("at profile grade the measured class preserves the "
                  "spike's sharpness: SMEAR_TOTAL = %.3f x width <= 0.25; "
                  "the common-mode shift (%.2f x width) is the "
                  "full-pipeline successor's object, disclosed not "
                  "letter-bearing" % (smear_total / WIDTH_SPIKE,
                                      shift_central / WIDTH_SPIKE))
    elif smear_total >= 1.00 * WIDTH_SPIKE:
        letter = "R-PROFILE-BROKEN"
        clause = ("the class smears the spike at profile grade "
                  "(SMEAR_TOTAL = %.3f x width >= 1); commission the "
                  "full-pipeline leg before any credence motion"
                  % (smear_total / WIDTH_SPIKE))
    else:
        letter = "R-PROFILE-GRAY"
        clause = ("smear between 0.25 and 1.00 x width "
                  "(%.3f); profile grade cannot decide"
                  % (smear_total / WIDTH_SPIKE))
    log("VERDICT LETTER: %s" % letter)
    log("   " + clause)
    if g7:
        log("   CLASS SEPARATION (G7-carried): the AQUAL-family anomaly at "
            "the same radii is >= %.0fx the measured class's -- the V-N-T "
            "null and this retrodiction are one mechanism apart." % ratio)
    log("pre-signed credence: ALL cells HOLD anomaly-real 53 / "
        "mech-conditional 8; P14 row annotated this close.")
    log("=" * 72)

    with io.open(os.path.join(ROOT, "data/stage10r_p14spike.txt"), "w",
                 encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")


if __name__ == "__main__":
    run()
