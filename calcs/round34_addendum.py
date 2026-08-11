"""
ROUND 34 verification addendum (stage 10O P15-ERODE).

Protocol (the 87a4676 two-halves standard, 5th execution):
  GA half -- BLIND: independent re-computation of the stage's own
  load-bearing numbers, run and committed BEFORE the ROUND 34 report
  exists. Methods are implementation-independent from the stage:
    GA-1  the phase-mixing map via the MICROCANONICAL CONTOUR QUADRATURE
          (no time integration; the M&H eq. 26 residence measure
          dw dj delta(H - H0) / |grad H| restricted to the connected
          component of the initial condition), with the Hamiltonian
          TRANSCRIBED FRESH from the primary read (catches a stage-side
          transcription typo). Uses the same A4 control-variate form
          (alpha_i + alpha_hat_weighted - alpha_hat_0 on the same ICs)
          so the method comparison is not drowned by IC-draw noise.
          Bars: |Delta| <= 0.06 near unity, 0.10 at the wing cells
          (-0.99, 2.5) -- two finite-ensemble methods, stage cells are
          R = 6 means (A6). Because the quadrature carries NO time
          integration and NO window operator, agreement also closes the
          A5 operator-asymmetry question independently.
    GA-2  the wedge boundary arithmetic re-derived in fresh code from
          the JSON map table + Table 1 (bar: exact to 1e-9 on boundary
          ratio values at tau_c in {0.3, 1, 3} Gyr + the frozen corner).
    GA-3  t0 normalization and the e_a*kappa_c conversion from scratch
          (bar: 1e-6 relative).
    GA-4  the step statistic re-fit via explicit normal equations
          (bar: chi2_null, Delta chi2, x0 to 1e-6).
    GA-5  the consistency read re-derived (bar: 1e-6 on map(1.32) reuse
          + d-sigma arithmetic).
  GB half -- POST-REPORT: re-computation of every load-bearing REVIEWER
  number (memory rule feedback-verify-reviewer-math). Filled in only
  after REVIEW-ROUND34-OPUS.md exists. Do not run GB before the report.

Usage: py calcs/round34_addendum.py GA   (or GB after the report)
"""

import json
import sys

import numpy as np
import sympy as sp
from scipy import ndimage

HALF = sys.argv[1] if len(sys.argv) > 1 else 'GA'
RES = json.load(open('data/stage10o_results.json'))

# constants, re-derived from scratch (GA-3 checks against the stage)
G_SI = 6.674e-11
MSUN = 1.989e30
AU = 1.496e11
PC = 3.086e16
GYR = 1e9 * 3.156e7
A0 = 1.2e-10
RHO0 = 0.2 * MSUN / PC**3
GAMMA_GAL = 4 * np.pi * G_SI * RHO0
AGE = 10.0

BIN_EDGES = np.array([1.50, 2.00, 2.25, 2.50, 2.75, 3.00, 3.50, 4.00, 4.50])
ALPHA_OBS = np.array([0.08, 0.59, 0.82, 0.94, 1.20, 1.30, 1.32, 1.17])
SIG_OBS = np.array([0.15, 0.10, 0.05, 0.04, 0.05, 0.06, 0.09, 0.15])
BIN_MID = 0.5 * (BIN_EDGES[:-1] + BIN_EDGES[1:])
A_J = 10.0**BIN_MID * AU

checks = []


def check(name, ok, detail):
    checks.append(ok)
    print(f"  {name}: {'PASS' if ok else 'FAIL'}  ({detail})")


def t0_gyr(a_m, gamma=GAMMA_GAL, mb=1.0):
    return 2 * np.pi * np.sqrt(G_SI * mb * MSUN / a_m**3) / gamma / GYR


MAP = {}
for k, v in RES['map_table'].items():
    g, a = k.split(',')
    MAP[(round(float(g), 3), float(a))] = v
GRID_A = sorted({a for (_, a) in MAP})


def map_interp(gam, a):
    key = round(gam, 3)
    vals = np.array([MAP[(key, x)] for x in GRID_A])
    return float(np.interp(a, GRID_A, vals))


if HALF == 'GA':
    print("=" * 70)
    print("GA (blind half, pre-report)")
    print("=" * 70)

    # ---------------- GA-1: contour-quadrature map ----------------
    # Hamiltonian transcribed FRESH from M&H 2023 eq. 30 (primary read
    # 2026-08-11): H_G = (1/j^2)[(j^2 - 3 G jz^2)(5 - 3 j^2)
    #                            - 15 G (j^2 - jz^2)(1 - j^2) cos 2w]
    j_s, jz_s, w_s, G_s = sp.symbols('j jz w G')
    H_e = (1 / j_s**2) * ((j_s**2 - 3 * G_s * jz_s**2) * (5 - 3 * j_s**2)
                          - 15 * G_s * (j_s**2 - jz_s**2) * (1 - j_s**2)
                          * sp.cos(2 * w_s))
    Hf = sp.lambdify((j_s, w_s, jz_s, G_s), H_e, 'numpy')
    dHj = sp.lambdify((j_s, w_s, jz_s, G_s), sp.diff(H_e, j_s), 'numpy')
    dHw = sp.lambdify((j_s, w_s, jz_s, G_s), sp.diff(H_e, w_s), 'numpy')

    def contour_map(alpha_i, gam, n_ic=600, ngrid=520, seed=202):
        rng = np.random.default_rng(seed)
        u = rng.random(n_ic)
        e0 = np.clip(u**(1 / (1 + alpha_i)), 1e-4, 1 - 1e-6)
        cosi = rng.uniform(-1, 1, n_ic)
        w0 = rng.uniform(0, 2 * np.pi, n_ic)
        j0 = np.sqrt(1 - e0**2)
        jz0 = j0 * cosi
        ee = e0[e0 > 1e-6]
        a_hat0 = -len(ee) / np.sum(np.log(ee)) - 1.0
        num = 0.0
        den = 0.0
        wgrid = np.linspace(0, 2 * np.pi, ngrid, endpoint=False)
        for m in range(n_ic):
            jlo = abs(jz0[m]) + 1e-6
            jgrid = np.linspace(jlo, 1 - 1e-9, ngrid)
            WW, JJ = np.meshgrid(wgrid, jgrid)
            Hg = Hf(JJ, WW, jz0[m], gam)
            H0 = Hf(j0[m], w0[m], jz0[m], gam)
            eps = 0.015 * max(np.std(Hg), 1e-6)
            band = np.abs(Hg - H0) < eps
            if not band.any():
                continue
            # connected component containing the IC (with omega wrapping)
            lab, nl = ndimage.label(band)
            # merge labels that touch across the omega seam
            seam_pairs = set()
            for r in range(ngrid):
                a_l, b_l = lab[r, 0], lab[r, -1]
                if a_l > 0 and b_l > 0 and a_l != b_l:
                    seam_pairs.add((a_l, b_l))
            parent = list(range(nl + 1))

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x

            for a_l, b_l in seam_pairs:
                ra, rb = find(a_l), find(b_l)
                if ra != rb:
                    parent[ra] = rb
            iw = int(np.argmin(np.abs((wgrid - w0[m]) % (2 * np.pi))))
            ij = int(np.clip(np.searchsorted(jgrid, j0[m]), 0, ngrid - 1))
            lab0 = lab[ij, iw]
            if lab0 == 0:
                # IC pixel missed the band (rare): nearest band pixel
                idx = np.argwhere(band)
                d2 = (idx[:, 0] - ij)**2 + (idx[:, 1] - iw)**2
                lab0 = lab[tuple(idx[int(np.argmin(d2))])]
            root0 = find(lab0)
            comp = np.zeros_like(band)
            for l_ in range(1, nl + 1):
                if find(l_) == root0:
                    comp |= (lab == l_)
            # GA-1 v2 FIX (caught by the blind cross-check itself, v1 log
            # archived as data/round34_ga_v1.log): the delta-band CELL
            # CENSUS already carries one 1/|grad H| factor (band thickness
            # ~ eps/|grad H|), so equal weight per band cell IS the
            # microcanonical/residence measure; v1's extra 1/|grad H|
            # double-counted it and over-weighted slow (librating-center)
            # segments -- biasing every cell HIGH, worst at high e
            # (contour 2.05 vs stage 1.71 at alpha_i = 2.5).
            wgt = np.where(comp, 1.0, 0.0)
            tot = np.sum(wgt)
            if tot <= 0:
                continue
            e_g = np.sqrt(np.clip(1 - JJ**2, 1e-12, 1))
            num += np.sum(wgt * np.log(e_g)) / tot
            den += 1.0
        a_raw = -den / num - 1.0
        return alpha_i + a_raw - a_hat0     # A4-form control variate

    targets = [(1 / 3, 1.0, 0.06), (1 / 3, 1.2, 0.06),
               (1.0, -0.99, 0.10), (1.0, 1.0, 0.06), (1.0, 2.5, 0.10)]
    ga1_ok = True
    for gam, a_i, bar in targets:
        a_ind = contour_map(a_i, gam)
        a_stage = MAP[(round(gam, 3), a_i)]
        d = abs(a_ind - a_stage)
        ga1_ok = ga1_ok and (d <= bar)
        print(f"  GA-1 cell Gamma={gam:.3f} a_i={a_i:+.2f}: contour "
              f"{a_ind:+.4f} vs stage {a_stage:+.4f} (d={d:.4f}, bar {bar})")
    check('GA-1', ga1_ok, "per-cell bars 0.06/0.10")

    # independent map(1.32) for the consistency read
    a132 = contour_map(1.32, 1 / 3)
    d132 = abs(a132 - RES['consistency_map132'])
    check('GA-5', d132 <= 0.06,
          f"map(1.32): contour {a132:.4f} vs stage {RES['consistency_map132']:.4f}")

    # ---------------- GA-2: wedge boundary re-derivation ----------------
    ALPHA_FLOOR = {round(float(k), 3): v for k, v in RES['alpha_floor'].items()}

    def reachable(k, gam, ceil):
        if k == 0:
            return -1.0, ceil
        lo = ALPHA_FLOOR[round(gam, 3)]
        hi = map_interp(gam, ceil)
        for _ in range(k - 1):
            lo = map_interp(gam, lo)
            hi = map_interp(gam, hi)
        return lo, hi

    TAU_GRID = np.geomspace(0.05, 10.0, 25)
    RATIO_GRID = np.geomspace(1e-2, 1e6, 33)
    CONV = A0 / np.sqrt(G_SI * MSUN / A0) / GAMMA_GAL
    kgal = np.array([1 if t0_gyr(a) <= AGE else 0 for a in A_J])
    worst2 = 0.0
    for tau_pick in (0.3, 1.0, 3.0, 10.0):
        i = int(np.argmin(np.abs(TAU_GRID - tau_pick)))
        tau = TAU_GRID[i]
        kc_base = int(np.floor(AGE / tau))
        bound = None
        for ratio in RATIO_GRID:
            bad = False
            for b in range(1, 8):
                tmix = t0_gyr(A_J[b]) / ratio
                kc = kc_base if (tmix <= tau and ratio >= 2.0) else 0
                lo, hi = reachable(kc + kgal[b], 1.0, 2.5)
                if (ALPHA_OBS[b] + 2 * SIG_OBS[b] < lo
                        or ALPHA_OBS[b] - 2 * SIG_OBS[b] > hi):
                    bad = True
                    break
            if bad:
                bound = ratio / CONV
                break
        key = str(float(tau_pick)) if tau_pick != 10.0 else str(AGE)
        stage_val = RES['wedge_primary'].get(key)
        stage_b = None if stage_val is None else stage_val[1]
        if bound is None and stage_b is None:
            d = 0.0
        elif bound is None or stage_b is None:
            d = np.inf
        else:
            d = abs(bound - stage_b) / stage_b
        worst2 = max(worst2, d)
        print(f"  GA-2 tau={tau_pick}: bound {bound} vs stage {stage_b}")
    check('GA-2', worst2 <= 1e-9, f"worst rel d = {worst2:.2e}")

    # GA-2b (informational, quantifies the x1.78 grid conservatism):
    # the same boundary on a 12-pts/decade ratio grid
    RATIO_FINE = np.geomspace(1e-2, 1e6, 97)
    for tau_pick in (0.3, 1.0, 3.0):
        i = int(np.argmin(np.abs(TAU_GRID - tau_pick)))
        tau = TAU_GRID[i]
        kc_base = int(np.floor(AGE / tau))
        bound = None
        for ratio in RATIO_FINE:
            bad = False
            for b in range(1, 8):
                tmix = t0_gyr(A_J[b]) / ratio
                kc = kc_base if (tmix <= tau and ratio >= 2.0) else 0
                lo, hi = reachable(kc + kgal[b], 1.0, 2.5)
                if (ALPHA_OBS[b] + 2 * SIG_OBS[b] < lo
                        or ALPHA_OBS[b] - 2 * SIG_OBS[b] > hi):
                    bad = True
                    break
            if bad:
                bound = ratio / CONV
                break
        print(f"  GA-2b tau={tau_pick}: fine-grid boundary product = {bound}")

    # ---------------- GA-3: normalization ----------------
    t0_ref = t0_gyr(1e4 * AU)
    d3a = abs(t0_ref - RES['t0_ref_gyr']) / RES['t0_ref_gyr']
    d3b = abs(1 / CONV - 1 / RES['conv_ratio_per_product']) * RES['conv_ratio_per_product']
    check('GA-3', d3a <= 1e-6 and d3b <= 1e-6,
          f"t0 rel d = {d3a:.2e}, conv rel d = {d3b:.2e}; t0 = {t0_ref:.3f} Gyr")

    # ---------------- GA-4: step statistic ----------------
    X = BIN_MID
    WI = 1 / SIG_OBS**2

    def wls(cols, y):
        A = np.vstack(cols).T
        M = (A * WI[:, None]).T @ A
        v = (A * WI[:, None]).T @ y
        c = np.linalg.solve(M, v)
        r = y - A @ c
        return float(np.sum(WI * r**2)), c

    chi2_null, _ = wls([np.ones_like(X), X, X**2], ALPHA_OBS)
    best = (np.inf, None)
    for x0 in np.arange(2.375, 3.751, 0.125):
        s_col = (X >= x0).astype(float)
        if s_col.sum() in (0, len(X)):
            continue
        c2, c = wls([np.ones_like(X), X, X**2, s_col], ALPHA_OBS)
        if c[3] > 0:
            continue
        if c2 < best[0]:
            best = (c2, x0)
    dchi2 = chi2_null - best[0] if np.isfinite(best[0]) else 0.0
    d4 = max(abs(chi2_null - RES['step_data']['chi2_null']),
             abs(dchi2 - RES['step_data']['dchi2']))
    x0_ok = (best[1] == RES['step_data']['x0'])
    check('GA-4', d4 <= 1e-6 and x0_ok,
          f"chi2_null {chi2_null:.3f}, dchi2 {dchi2:.3f}, x0 {best[1]}")

    n_ok = sum(checks)
    print(f"GA SUMMARY: {n_ok}/{len(checks)} PASS")

elif HALF == 'GB':
    print("GB half: filled after REVIEW-ROUND34-OPUS.md exists.")
    print("(re-compute every load-bearing reviewer number here)")
