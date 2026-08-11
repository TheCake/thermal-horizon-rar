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
            # ROUND-34 ANNOTATION (adoption condition 3): v2's REMAINING
            # wing-cell disagreement (2.022 at alpha_i = 2.5; -0.774 at
            # -0.99) is ALSO this method's artifact -- a residual global-
            # eps band-census bias at the 1/j^2 corner. The round's two
            # independent re-derivations (1-D Gauss-Chebyshev residence
            # quadrature 1.725 and independent ODE 1.741 / -0.903) plus
            # the session's aliasing probe (1.731) all confirm the STAGE
            # values (1.709 / -0.905). Do NOT treat GA-1 v2 wing values
            # as a competing measurement; near-unity cells (where this
            # method is clean) remain the valid cross-validation.
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
    print("=" * 70)
    print("GB (post-report): re-computing ROUND 34 load-bearing numbers")
    print("=" * 70)

    def load_map(path):
        rr = json.load(open(path))
        mm = {}
        for k, v in rr['map_table'].items():
            g, a = k.split(',')
            mm[(round(float(g), 3), float(a))] = v
        return mm

    def mi(mm, gam, a):
        key = round(gam, 3)
        gr = sorted({aa for (g2, aa) in mm if g2 == key})
        vals = np.array([mm[(key, x)] for x in gr])
        return float(np.interp(a, gr, vals))

    M6 = MAP

    # GB-1: B10 pinning from starts 1.71/2.02/2.5/3.0 (his 0.978/0.979/
    # 0.983/0.983; interp clamps above the 2.5 table edge as his did)
    his_b10 = {1.71: 0.978, 2.02: 0.979, 2.5: 0.983, 3.0: 0.983}
    ok1 = True
    for start, hv in his_b10.items():
        b = start
        for _ in range(10):
            b = mi(M6, 1.0, min(b, 2.5))
        ok1 = ok1 and abs(b - hv) <= 0.02
        print(f"  GB-1 B10 from {start}: {b:.4f} vs his {hv}")
    check('GB-1', ok1, "B10 pinning")

    # GB-2: the run-6 map fixed point (his ~0.973)
    fp = 1.5
    for _ in range(60):
        fp = mi(M6, 1.0, fp)
    check('GB-2', abs(fp - 0.973) <= 0.01, f"fixed point {fp:.4f} vs his 0.973")

    # GB-3: the tau=3.3 flip arithmetic (his: run-5 B4 = 1.076 < 1.14;
    # run-6 B4 = 1.147 > 1.14; and B4 from 1.71/2.02 = 1.094/1.112)
    M5 = load_map('data/stage10o_results_run5.json')
    b45 = 2.5
    for _ in range(4):
        b45 = mi(M5, 1.0, min(b45, 2.5))
    b46 = 2.5
    for _ in range(4):
        b46 = mi(M6, 1.0, min(b46, 2.5))
    ok3 = abs(b45 - 1.076) <= 0.01 and abs(b46 - 1.147) <= 0.01
    for start, hv in ((1.71, 1.094), (2.02, 1.112)):
        b = start
        for _ in range(4):
            b = mi(M6, 1.0, min(b, 2.5))
        ok3 = ok3 and abs(b - hv) <= 0.01
        print(f"  GB-3 B4 from {start}: {b:.4f} vs his {hv}")
    check('GB-3', ok3, f"run5 B4 {b45:.4f} (<1.14), run6 B4 {b46:.4f} (>1.14)")

    # GB-4: the P10 exclusion boundary in tau (his 2.48 Gyr, amplitude-
    # saturated over 0.086/0.5/1.43)
    ALPHA_FLOOR = {round(float(k), 3): v for k, v in RES['alpha_floor'].items()}

    def reachable6(k, ceil=2.5):
        if k == 0:
            return -1.0, ceil
        lo = ALPHA_FLOOR[1.0]
        hi = mi(M6, 1.0, ceil)
        for _ in range(k - 1):
            lo = mi(M6, 1.0, min(lo, 2.5))
            hi = mi(M6, 1.0, min(hi, 2.5))
        return lo, hi

    kgal = np.array([1 if t0_gyr(a) <= AGE else 0 for a in A_J])
    CONV = A0 / np.sqrt(G_SI * MSUN / A0) / GAMMA_GAL
    ok4 = True
    for prod in (0.086, 0.5, 1.43):
        ratio = prod * CONV
        tau_excl = None
        for tau in np.geomspace(1.0, 5.0, 400):
            kc_base = int(np.floor(AGE / tau))
            bad = False
            for b in range(1, 8):
                tmix = t0_gyr(A_J[b]) / ratio
                kc = kc_base if (tmix <= tau and ratio >= 2.0) else 0
                lo, hi = reachable6(kc + kgal[b])
                if (ALPHA_OBS[b] + 2 * SIG_OBS[b] < lo
                        or ALPHA_OBS[b] - 2 * SIG_OBS[b] > hi):
                    bad = True
                    break
            if bad:
                tau_excl = tau
        ok4 = ok4 and tau_excl is not None and abs(tau_excl - 2.48) <= 0.05
        print(f"  GB-4 product {prod}: excluded up to tau = "
              f"{tau_excl if tau_excl is None else round(tau_excl, 3)} Gyr vs his 2.48")
    check('GB-4', ok4, "P10 boundary 2.48 Gyr amplitude-saturated")

    # GB-5: frozen-corner reachable interval contains every bin at 2 sigma
    lo1, hi1 = reachable6(1)
    inside = all(ALPHA_OBS[b] + 2 * SIG_OBS[b] >= lo1
                 and ALPHA_OBS[b] - 2 * SIG_OBS[b] <= hi1 for b in range(8))
    check('GB-5', inside, f"one-mixing interval [{lo1:.3f}, {hi1:.3f}] contains all bins")

    # GB-6: rho0 = 0.1 timescales (his 16.5 / 3.1 Gyr for bins 7/8)
    gam01 = 4 * np.pi * G_SI * (0.1 * MSUN / PC**3)
    t7 = t0_gyr(A_J[6], gamma=gam01)
    t8 = t0_gyr(A_J[7], gamma=gam01)
    check('GB-6', abs(t7 - 16.5) <= 0.2 and abs(t8 - 3.1) <= 0.1,
          f"bin7 {t7:.2f} / bin8 {t8:.2f} Gyr")

    # GB-7 (informational): the main session's aliasing probe on
    # (Gamma=1, 2.5) gave CV 1.7317/1.7312 (neval 350/1400, fresh seed 9)
    # -- inside his 1.71-1.74 cluster; a fourth independent confirmation
    # of the stage wing value.
    print("  GB-7 (info): alias-probe CV 1.7312-1.7317 in the 1.71-1.74 cluster")

    n_ok = sum(checks)
    print(f"GB SUMMARY: {n_ok}/{len(checks)} PASS")
