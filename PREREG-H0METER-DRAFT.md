# PRE-REGISTRATION — stage 10S "H0METER" (2026-09-20)
# Status: COMMITTED on the author's word ("let's run it", 2026-09-20
# evening), AFTER the full instrument build and BEFORE any sky fit.
# History: drafted 2026-09-20 morning; instrument-side start authorized
# same day ("let's start the Hubble-meter anyway... we don't need to
# commit anything"); gates G1-G6 + UMa desk check executed sky-blind
# (record below); final pins 6-10 written pre-commit, pre-sky.

## Question
If a₀ = cH₀/2π (the lock), the measured acceleration scale is an H₀ meter.
SPARC's distance provenance contaminates the naive reading: 97/175 galaxies
(55%) carry Hubble-flow distances that assume H₀ = 73 (SPARC table note f_D=1
verbatim); 50 (29%) are H₀-independent anchors (45 TRGB, 3 Cepheid, 2 SN);
28 (16%) share one Ursa Major cluster distance of unverified provenance.
Two legs, one agreement test:
- LEG A (anchored): refit a₀ on the anchor subsample only → H₀_A = 2π a₀_A/c.
- LEG B (flow, self-consistent): for flow galaxies, distances scale as
  D ∝ 73/H₀'; solve a₀_fit(H₀') = cH₀'/2π for the unique fixed point → H₀_B.
- THE AGREEMENT TEST: under the lock, H₀_A = H₀_B is forced; a numerical
  coincidence carries no such constraint. Agreement is the paper's spine.

## Convention pins (fixed before any fit)
1. Function: a₀ is quoted under the BE form (ν = 1 + n_BE) as PRIMARY, with
   the ledger-function family spread reported as a systematic band alongside
   (P2 Table 3 row-3 machinery). Never quote a single-function H₀ without the
   family band: the function spread (1.044–1.125) is itself of order the
   Planck–SH0ES gap.
2. Treatment: hierarchical with measured distance/inclination priors (the
   stage5m_hierv.py vertical treatment) = PRIMARY; flat-M/L sparc_rar_fit.py
   = co-read. Distance priors for flow galaxies re-CENTER with H₀' in leg B;
   prior widths keep the SPARC fractional errors.
3. Scaling laws for distance rescaling s = D'/D: R' = sR, V_bar'² = sV_bar²,
   V_obs unchanged ⇒ g_bar invariant, g_obs ∝ 1/s (derived; VERIFY by
   injection, gate G3 below — do not trust the hand-wave).
4. UMa subsample (28 gal): classification rule = read the SPARC/Verheijen
   provenance of the adopted UMa distance BEFORE any fit; if its calibration
   is anchor-based (Cepheid-calibrated TF zero point), UMa joins leg A with a
   shared-distance nuisance; if flow-informed, joins leg B with its own
   scaling; if ambiguous, EXCLUDED from both legs (report all three counts).
5. H₀ conversions: H₀ = 2π a₀/c, 1 km/s/Mpc = 3.240779e-20 s⁻¹; all H₀ in
   km/s/Mpc; machine-to-machine constants (never retype).

## Gates (instrument-side; all must pass before sky legs)
- G1 (subsample regression): leg-A machinery on the FULL sample must
  reproduce the archived a₀ = 1.05 ± 0.10 (flat) and the hier row point
  estimates to 0.5% — the meter inherits, not reinvents.
- G2 (anchor census): the 50/97/28 split reproduced from the .mrt by two
  independent parsers (token + anchored regex; NOTE 2026-09-20: the .mrt's
  true layout deviates from its byte-by-byte header spec, so fixed-width
  slicing is invalid — the census finding, reproduced at harness time);
  any disagreement resolved by eye before proceeding.
- G3 (scaling injection): rescale ALL flow galaxies by s ∈ {0.85, 1.0, 1.15},
  refit; measured d(ln a₀)/d(ln s) must match the analytic −2 within the
  fit's own regime mixing (report the measured exponent; bar: |exponent|
  ∈ [1.5, 2.5], else the leg-B solve uses the MEASURED response curve — the
  bar gates the ANALYTIC SHORTCUT, not the stage).
- G4 (leg-B solver) — AMENDED 2026-09-20 pre-commit, first firing: the
  original single-noisy-mock 1% bar conflated solver wiring with mock
  realization noise (the fixed-point geometry amplifies any fractional a₀
  offset by 1/(γ−1); at the G3-measured γ = 1.585 that is ×1.71, so one
  noisy realization sits at ~3–5% by construction — the bar was mis-posed,
  trap-#22 family: measure between-draw scatter before signing a bar).
  Split:
  - G4a (wiring): NOISELESS synthetic recovers H₀_true to 0.2% at three
    values spanning [60, 80].
  - G4b (calibration): noisy replications (≥5 seeds at the outer truths;
    measured per-point + 0.08 dex intrinsic noise): recovery UNBIASED
    (|mean error| ≤ max(1.5%, 2·SE of the mean)); the replication SD is
    REPORTED as the leg-B single-realization noise floor.
  Leg-B uncertainty at sky time: galaxy bootstrap of the FULL fixed-point
  solve (≥100 reps) — this carries the ×1/(γ−1) amplification
  automatically; the G4b floor is its cross-check.
- G5 (power): leg-A anchored subsample (N≈50): bootstrap σ(a₀_A) measured
  BEFORE the sky read; if σ(H₀_A) > 15 km/s/Mpc the leg is POWER-LIMITED and
  the letter says so (no silent quoting of a useless leg).
- G6 (literature bar): one targeted search round for any published
  a₀→H₀-as-measurement execution (beyond the 2026-09-20 scout NOT-FOUND);
  logged in the stage docstring with date.
- G7 (composition disclosure — the dissident-dwarf interaction): the anchored
  subsample is nearby-dominated and therefore dwarf/gas-dominated-heavy
  (DDO154-class members confirmed at census). Leg A quotes a₀ BOTH with and
  without the gas-dominated class, both printed; if the two differ by more
  than the leg's own 1σ, the letter carries a named GD-composition caveat and
  no headline H₀_A is quoted without it. (G2 census result: 97/45/3/28/2,
  two parsers exact match, 2026-09-20.)

## Letters (grammar fixed in advance; every clause gated — trap #11/#12)
- M-AGREE: |H₀_A − H₀_B| within the joint 1σ AND both legs pass their gates
  → "the two legs agree; the lock passes its first self-consistency test as
  an instrument; H₀(meter) = <joint>." (Positive clause gated by G1–G5.)
- M-SPLIT: legs disagree > 2σ joint → "the lock fails the agreement test at
  current grade OR a distance-provenance systematic is unmodeled; the meter
  does not return an H₀." (Named-suspect list mandatory.)
- M-GRAY: between 1σ and 2σ, or exactly one leg power-limited → report both
  numbers, no headline H₀, successor named.
- M-POWER-DEAD: G5 fails on A and the B-solve is degenerate → the meter
  awaits better anchors (report the roadmap only).
- Prohibitions: never print a single-leg H₀ as "the" measurement; never
  print any H₀ without the function-family band; never call the result
  "ladder-independent" (leg A USES ladder rungs — the honest term is
  "H₀-assumption-independent"); the two-grade rule applies to any
  "tension-relevant" phrasing.

## Credence
No cell of the anomaly-real or mechanism maps moves on ANY outcome of this
stage (it measures a consequence of the lock, not the anomaly). If M-AGREE
fires, the lock's standing improves at REPORT grade only; any credence
implication waits for a dedicated pre-signed map (author + review round).

## Desk-check record (pin 4 executed 2026-09-20, primary read, pre-fit)
The SPARC paper (Lelli+2016, AJ 152, 157; PDF read directly) Sec. 2: "28
objects lie in the Ursa Major cluster (Verheijen & Sancisi 2001), which has
an average distance of 18 ± 0.9 Mpc (Sorce et al. 2013b)... estimated depth
of ~2.3 Mpc (Verheijen 2001), hence we adopt an error of sqrt(2.3^2+0.9^2)
= 2.5 Mpc for individual galaxies." Sorce et al. 2013b = ApJ 765, 94
(arXiv:1301.4833), the mid-IR Tully-Fisher CALIBRATION paper; its zero point
= "26 galaxies with Cepheid or tip of the red giant branch distances that
define the zero point" (abstract, verbatim). The adopted UMa distance is a
TFR distance on a Cepheid/TRGB-calibrated zero point — ANCHOR-BASED, not
flow-informed. RULING (per pin 4): **UMa JOINS LEG A** with a
shared-distance nuisance (shared component 0.9 Mpc on 18 Mpc; per-galaxy
depth 2.3 Mpc independent). Leg-A raw membership = 50 + 28 = 78; leg B = 97.
Cross-checks: SPARC's own reassignments NGC3992 -> f_D=5 (SN, 23.7 Mpc) and
UGC06446 -> f_D=1 (flow, 12 Mpc) verified in the .mrt (the SN count of 2
includes NGC3992). Provenance note for the paper: SPARC chose H0 = 73
explicitly to put flow distances "on a similar zero-point scale" as the
anchored groups (Tully et al. 2013 EDD scale) — the flow leg's H0
assumption is exactly the dial leg B inverts.

## Pre-commit final pins (2026-09-20, written BEFORE commit and BEFORE any
## sky fit; all direction-neutral, set by principle)
6. LENSING EXCLUDED from both legs: the 5M joint-lensing block inherits its
   own cosmology's H0 through the lensing catalog's distances; including it
   would break the H0-assumption-independence claim. Leg fits = the 5M
   machinery MINUS the lensing term (dlt fixed at 0). The G1b regression
   (full sample WITH lensing) stands as the machinery-identity gate; the
   no-lensing delta is one deleted term, and the leg-A runner carries an
   identity check: the UMa shared-nuisance parameter at prior width 1e-6
   must reproduce the plain no-nuisance fit.
7. UMa shared-distance nuisance (leg A, hier): one shared log-offset
   applied to the 26 UMa galaxies' vertical channel, Gaussian prior width
   (0.9/18)/ln10 dex; per-galaxy sigma_v keeps the .mrt e_D (depth) as is.
8. Leg-B solve: PRIMARY central = hier (no-lensing) a0 fitted on the flow
   subsample at trial H0' grid {60, 65, 70, 75, 80} (distances re-centered,
   fractional prior widths kept), crossing of a0_hier(H0') vs cH0'/2pi by
   linear interpolation between bracketing grid nodes; if unbracketed,
   extend by 5 to at most [50, 95], else the no-fixed-point outcome is
   reported per the letter grammar. ENGINE + uncertainty = the gated flat
   BE solve (G4 machinery, bracket [55, 90]): sigma_B = galaxy bootstrap
   of the full flat solve (150 reps), cross-checked against the G4b floor.
   sigma_A = the G5 bootstrap machinery (flat, 200 reps, distance jitter +
   UMa shared draw). Joint sigma for the agreement test = quadrature.
9. Function-family band (pin 1 executed as): hier family {p065, gm, boot}
   re-fit alongside BE on leg A; flat-solve family variants on leg B; both
   bands quoted next to every H0.
10. G7 selector inherited verbatim from stage 8S: GDFRAC = per-galaxy share
   of points with g_gas > g_dsk + g_bul at f = 1; gas-dominated iff
   GDFRAC >= 0.5.

## Gate record (instrument build executed 2026-09-20, sky-blind)
Harness: calcs/stage10s_h0meter.py -> data/stage10s_h0meter_gates.txt.
ALL GATES GREEN; no sky central computed or printed anywhere.
- G2 PASS: two parsers agree, 97/45/3/28/2; NGC3992->SN and UGC06446->flow
  reassignments verified in the .mrt.
- G1a PASS: flat 4H machinery reproduces archived a0 to 0.01% (p, f_d
  exact). G1b PASS: hier 5M machinery reproduces -12152.49 to 0.00 and
  a0 = 1.044e-10 to 0.03%.
- G3 PASS: measured d(ln a0)/d(ln s) = -1.585 on the flow subsample
  (in-bar; regime mixing + f_d freedom soften the analytic -2 as
  anticipated). Only ratios printed.
- G4a PASS (wiring): noiseless recovery exact to 0.001% at 62/70/78.
  G4b PASS (calibration): unbiased (mean -0.07%/-0.05% at the outer
  truths); single-realization floor SD ~2.5% of H0 (~1.7 km/s/Mpc);
  fixed-point amplification 1/(gamma-1) = 1.71.
- G5 PASS (powered): leg A after quality cuts = 67 galaxies (41 anchored
  + 26 UMa), sigma(H0_A) = 12.2 km/s/Mpc < 15 (central MASKED). The
  agreement test will be LEG-A-LIMITED; expect a wide joint sigma.
- G6 PASS: second targeted round found no a0->H0 lock-inversion
  execution. NEAREST NEIGHBOR (new, must-cite): Schombert, McGaugh &
  Lelli 2020 (AJ 160, 71) — bTFR-as-ladder H0 = 75.1 +- 2.3 +- 1.5 on
  the same catalog and same anchor/flow split, no lock in the method;
  their 50 anchors cross-validate the census.
- G7 armed for sky time (GD split of leg A).
Amendment log: G4 split into G4a/G4b pre-commit, first firing (the
original 1% single-noisy-mock bar was mis-posed; dated note at the gate).
G2 parser wording corrected to token+regex (the .mrt byte-spec deviation,
census finding reproduced).

## Order of work (firewall)
1. Instrument build + G1–G5 on synthetic/archived inputs (sky-blind). DONE
   2026-09-20 (all green; record above).
2. UMa provenance desk check (pin 4). DONE 2026-09-20 (UMa -> leg A;
   record above).
3. COMMIT this pre-reg (author's word) →
4. Sky legs fire in one run; letters per grammar; review round; then the
   P2 journal-pass paragraph and/or Proposal F drafting.
   Leg-B uncertainty = galaxy bootstrap of the full fixed-point solve
   (>=100 reps); leg-A = the G5 bootstrap machinery unmasked, hier
   co-read per pin 2.
