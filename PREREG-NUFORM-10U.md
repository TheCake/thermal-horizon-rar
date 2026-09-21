# PRE-REGISTRATION — STAGE 10U: PIN THE NU-FORM (the anchored-leg function-family contest)

Date: 2026-09-21 (early morning, the session after the 10T/R46 close).
Author go: "You are compacted now. Lets go" on the standing recommendation
(pin-the-nu-form first, dipole after — the 10T verdict's named successor).
Committed BEFORE any contest statistic, power calibration, or paired
bootstrap exists. Bars in this file were written blind to all of them.

## 0. Question and blindness inventory

The 10T operative verdict (data/stage10t_verdict.txt) left the anchored
leg SYSTEMATICS-LIMITED: the nu-form family band 65.4–70.4 km/s/Mpc is
the dominant term over the stat 5.0. This stage asks the one in-catalog
question that could shrink it: **can the 78-galaxy anchored leg DEMOTE
any member of the registered four-member family at calibrated grade —
or is the band the honest floor at this N?**

PRE-KNOWN (public in data/stage10t_skyread.txt, commit 87e689e/d4e1652):
the four members' deep-fit centrals on the 78 leg —
BE a0 = 1.0106e-10 (H0 65.4), p065 1.0892e-10 (70.4),
gm 1.0250e-10 (66.3), boot 1.0446e-10 (67.6); BE params
f_ML = 1.050, s_int = 0.032, u = −0.0211; primary boot sigma 4.99.

BLIND (never computed anywhere): the four members' profiled −2lnL
VALUES and their pairwise deltas on the 78 leg; the injection power
calibration of those deltas; the paired-bootstrap distribution of the
deltas; any demotion ruling. The gates mode computes sky member fits
only inside the G10U-1 regression and MASKS their objective values
(a0 regressions printed, funs never printed — the R46-C12 masking
pattern).

EXPECTATION DISCLOSURE (so the record shows what we knew): on the FULL
153-galaxy catalog the vertical ladder gave lean-grade preferences for
sharper-tail members (5M/6C/6L record). Those leans are NOT imported
here as priors, bars, or tie-breakers; the anchored-leg contest is
standalone. The honest prior stated in advance: at N = 78 the contest
may well be power-limited — and a measured power limit is itself the
deliverable (it prices the meter's bottleneck for Proposal F and sets
the BIG-SPARC refire requirement).

## 1. Pins (inherited, frozen)

1. **World**: `make_world(reclass=True, lt_slot=None)` from
   calcs/stage10t_legregrow.py — the R46 primary 78-galaxy leg A
   (41 SPARC-anchored + 26 UMa + 11 CF4-reclassified; UGC09992 zero
   points; LT excluded by the standing G10T-3b ruling). No census
   changes; KK98-251 / D564-8 / D631-7 ride at the NEXT census as
   booked. Machinery copied BIT-VERBATIM into calcs/stage10u_nuform.py
   (the 10T module exits in gates mode on import, so verbatim copy +
   identity regression G10U-1 is the inheritance route; LT loader
   stubbed inert — lt_slot=None always, the LT block of make_world is
   never evaluated).
2. **Family**: exactly the four registered meter members — BE
   (primary), p065, gm, boot — implementations verbatim from 10T
   (FAMS). NO new functions (PREDICTIONS §0 freeze). AMB is NOT in the
   meter family (requires per-galaxy ambient gates not uniformly
   available on the anchored leg) — stated, not contested here.
3. **Engine**: `fit_deep` (tol 5e-4, max_rounds 60, use_u=True) = the
   R46 operative estimator. The **contest estimator** per member =
   best-of-two-starts fit_deep (cold th0=None AND warm; warm = the
   archived convention: BE warm from the plain→loose chain, rivals
   warm from the BE deep optimum in sky mode / from truth params in
   injection worlds), keeping the smaller objective. Rationale
   pre-stated: R46-F3 showed start-dependence is the known hazard and
   the deeper profiler arbitrates. Objective values are profiled −2lnL
   up to one common additive constant (identical data, identical
   parameter count per member: la0, f, s_int, u + per-galaxy dml/dv),
   so cross-member deltas are meaningful.
4. **a0 → H0**: H0 = 2*pi*a0/c / KMS_MPC (the lock; constants pinned in
   the script).
5. **Precision rule** (R46 standing): −2lnL deltas printed at 0.1
   precision; never quote a cross-member delta below the measured
   convergence-repeat spread; never quote 0.01-grade deltas.
6. **Credences**: NO cell moves on ANY letter (pre-signed; 53 / 8
   stand). A pinned or narrowed band changes the three-world SCORING
   (PREDICTIONS SF annotation, append-only) but no credence cell.
7. **Verdict file rule**: data/stage10u_verdict.txt (post-round) may
   supersede ONLY the family-band clause of stage10t_verdict.txt; the
   10T file is never edited.
8. **Outputs**: data/stage10u_gates.txt (sky-blind), then
   data/stage10u_skyread.txt; every run preserved on amendment.

## 2. Statistic

For members A, B on a world w: d_AB(w) = fun_B(w) − fun_A(w), fun = the
contest-estimator profiled −2lnL. Positive d_AB = A beats B.
Sky observables: the six unordered deltas (printed 0.1 grade) and the
sky-best member B_sky = argmin fun.

## 3. Gates (each names what its failure vetoes)

- **G10U-1 REGRESSION** (STOP-grade; vetoes everything): the copied
  world + engine reproduce the ARCHIVED 10T deep fits on the 78 leg
  under the ARCHIVED start convention: |a0 − archive|/archive ≤ 0.1%
  for all four members (archive values in §0); BE params match at
  printed precision. Gates mode prints a0 regressions ONLY (funs
  masked). World counts asserted: legA = 78, flow = 71, UMa = 26.
- **G10U-2 MEMBER IDENTITY** (STOP-grade; vetoes everything): each
  FAMS member vs an INDEPENDENT scalar implementation (mpmath
  root-solve of the defining relation) at y ∈ {0.03, 0.1, 0.3, 1, 3,
  10}: |d nu|/nu ≤ 1e-6, and monotonicity nu' < 0 checked on the grid.
- **G10U-3 POWER** (the letter arm; an unpowered pair can never
  produce a demotion): 4 truths × 8 seeds (11, 23, 42, 101, 202, 303,
  404, 505; per-truth independent rng streams, the 10T G10T-4
  pattern). Generative model VERBATIM G10T-4 noisy arm with nu_T and
  a0_T = the member's own public 10T central: base =
  log10(gN1 * nu_T(gN1/a0_T)), gN1 = gg + 1.0*gd + gb on the 78 world;
  per-point noise N(0, sqrt(sig2 + 0.08^2)); per-galaxy distance
  offsets N(0, sigv_g). All four members fit per world (contest
  estimator).
  - **Self-recovery clause** (STOP-grade, trap #26): for every truth T
    and rival R, mean_seeds d_TR ≥ −2·SE — the truth member never
    SIGNIFICANTLY loses to a rival on its own worlds. Failure = the
    contest instrument is broken; no sky read.
  - **Separation** (R43 standing rule — NP separation measured before
    any bar is applied): sep(A,B) = (mean d_AB|truth A − mean
    d_AB|truth B) / sqrt((SD_A^2 + SD_B^2)/2). Ordered pair POWERED
    iff sep ≥ 2.0. With 8+8 seeds sep itself carries ~0.35·sep
    sampling noise — disclosed; powered/unpowered is therefore a
    MEASURED, lean-grade boundary and is printed with its inputs.
- **G10U-4 CONVERGENCE REPEAT** (STOP-grade above 2.0; also feeds the
  demotion bar): on the truth-BE seed-11 world, every member fit from
  both starts; c_conv = max over members |fun_cold − fun_warm|. If
  c_conv > 2.0: cross-member deltas at few-unit grade are unquotable —
  STOP, fix the instrument. Otherwise c_conv is ADDED to every
  demotion bar.
- **G10U-5 PAIRED BOOTSTRAP + HEALTH** (sky mode; its failure vetoes
  every demotion clause, capping the letter at U-FORM-LEAN): 200 reps,
  seed 202, the EXACT sigma_boot_hier draw stream (same resamples as
  the R46 primary boot rep-for-rep); per rep ALL FOUR members fit
  (warm from each member's own sky optimum — engine-matched). 
  - **G10U-5a stream identity**: the BE marginal of the paired boot
    must reproduce the archived primary boot: sigma within 0.5% of
    4.99 and percentiles [60.8, 65.4, 70.5] within 0.2 each.
  - **G10U-5b health (EVALUATED, printed PASS/FAIL — R46 rule)**: for
    the sky-best pair, |median_boot d − d_sky| ≤ 1.0 · SD_boot(d).
- **G10U-6 N-FORECAST** (diagnostic ONLY; vetoes nothing; labeled
  FORECAST): the band-decisive pair (BE, p065) at N-multipliers 1, 2,
  4 (world duplicated with independent noise), 4 seeds per truth side,
  cold fits: sep_hat(N) + a sqrt(N) extrapolation → the N at which
  sep = 2 and sep = 3. The BIG-SPARC requirement number for Proposal F.

## 4. Demotion rule and letters

Member B is **DEMOTED** iff some member A satisfies ALL of:
  (i) the ordered pair (A, B) is POWERED (G10U-3);
 (ii) sky d_AB > m_B + 2·s_B + c_conv, where (m_B, s_B) = mean/SD of
      d_AB over the truth-B injection seeds — i.e., the observed lead
      is inconsistent at the 2-SD grade with B being true, with the
      convergence spread charged against the demotion;
(iii) paired bootstrap: frac(reps with d_AB > 0) ≥ 0.95, health PASS.

Multiplicity disclosed: three rival tests per member at 2-SD grade;
the bootstrap co-requirement is the family-wise control; any demotion
is quoted as REPORT-grade, sharpenable by seed extension at review.

Letters (decision table, exhaustive):
- **U-FORM-PINNED**: exactly one survivor. The quoted meter band
  collapses to that member's H0 ± stat (verdict file supersedes the
  10T band clause).
- **U-FORM-NARROWED**: ≥ 1 demotion, ≥ 2 survivors. Band = survivor
  span.
- **U-FORM-LEAN**: ≥ 1 powered pair, no demotion. Band 65.4–70.4
  stands; leans quoted at bootstrap grade (9V rule: always quote the
  bootstrap grade, never the raw delta alone).
- **U-POWER-LIMITED**: no powered pair. Band stands as the measured
  floor; the power table + N-forecast are the deliverable. Sky deltas
  print WITH the power table adjacent, labeled "point deltas, no
  calibrated grade".
- Instrument STOPs (G10U-1/2/3-self/4): NO sky read; fix first.

In every branch: the flat treatment variant, membership band, peg
clause, and validity domain of the 10T verdict are UNTOUCHED.

## 5. Successor mandates (whichever letter fires)

- Proposal F absorbs the letter as the meter's bottleneck statement
  (a pinned/narrowed band, or the measured N-requirement).
- BIG-SPARC refire inherits THIS harness (contest + power arm) as-is.
- The external form-pinning levers stay registered where they are:
  a0(z) sign (P3 lens), DR4 binaries, void asymptote (P1).

## 6. Review

Round 47, fresh Opus reviewer, full package. GA blind half committed
BEFORE the report is read (the 87a4676 protocol, 14th execution); GB
half re-computes every load-bearing reviewer number. Amendments only
pre-quote, runs preserved.

---
## AMENDMENT A1 (2026-09-21, PRE-SKY, pre-quote; run 1 preserved in
data/stage10u_gates_run1.txt)

The run-1 power matrix put two unordered pairs INSIDE the sep
estimator's own sampling noise of the 2.0 bar: BE–p065 at 1.93 and
gm–boot at 2.21 (8+8 seeds carry ~0.35·sep noise; the diagnostic
forecast arm scored BE–p065 at 3.22 on independent seeds — the
powered/unpowered call for these pairs is currently a coin toss on the
SD estimate). Rule, fixed before any extension world is fit:

1. Any unordered pair with |sep8 − 2.0| ≤ 0.5 is re-measured with 16
   ADDITIONAL seeds per truth side (seeds 606, 707, 808, 909, 1010,
   1111, 1212, 1313, 1414, 1515, 1616, 1717, 1818, 1919, 2020, 2121;
   same generative model and rng scheme). Applies to BE–p065 and
   gm–boot and to no other pair.
2. For those pairs the POOLED 24-seed separation and truth-side
   calibration stats (m_B, s_B) REPLACE the 8-seed values everywhere
   (powered call, demotion bar). The bar stays 2.0. Direction-neutral:
   the extension can promote OR demote a pair.
3. Diagnostic re-scope, disclosed: the run-1 ×4 forecast point (sep
   1.97 against the ×1→×2 trend 3.22→5.48, with 4-seed SD blow-up and
   rounds-cap risk on 312-galaxy cold fits) is flagged
   convergence-suspect; the sqrt(N) forecast coefficient is computed
   from the ×1–×2 points only, the ×4 row printed as-is with the flag.
   (FORECAST grade throughout; vetoes nothing.)

---
GATE RECORD (appended after the sky-blind gates run, before the sky
read): full record = data/stage10u_gates.txt (A1 form; run 1 =
data/stage10u_gates_run1.txt). Summary: G10U-1 PASS (all four archived
a0 reproduced <= 0.004%, BE params exact), G10U-2 PASS (independent
solves <= 2e-15), G10U-3 self-recovery PASS (truth never loses),
G10U-4 PASS (c_conv = 0.00 — no convergence penalty in the bars),
A1 pooled separations: BE–p065 1.93 -> 1.80 UNPOWERED, gm–boot
2.21 -> 1.61 UNPOWERED (the extension DEMOTED a run-1 powered pair —
direction-neutrality exercised). Final powered ordered pairs (6/12):
(BE,gm), (gm,BE), (BE,boot), (boot,BE), (p065,boot), (boot,p065).
Structural pre-sky consequence, stated before the contest is read:
the band EDGES (BE, p065) cannot separate each other at N = 78; the
only demotions that could move the quoted band are boot->p065 (top
edge to gm/boot) or gm/boot->BE (low edge up). N-forecast (x1-x2):
sep ~ 3.55 sqrt(N/78) FORECAST-grade, in tension with the direct
24-seed sep(78) = 1.80 — the direct measurement is operative; the
tension is flagged for Round 47.
