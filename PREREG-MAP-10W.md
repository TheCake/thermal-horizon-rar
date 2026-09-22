# PREREG — STAGE 10W: THE PER-GALAXY FORM-PREFERENCE MAP + THE ANGULAR-RESOLUTION CUT

Registered 2026-09-22, sky-blind: every bar, seed, grid, letter and
veto below is frozen before any Part-2 or Part-3 statistic exists.
This is the R48-C18 successor instrument (one map replacing stratum
contests) plus the R48 named live rival (inner angular resolution)
given its registered test. Round 49 reviews the fired stage.

## 1. Objects and prior exposure (owned up front)

The 10U cross-leg inversion (anchored 78 prefer BE; flow 71 prefer
boot, BE last +48.7) was dissolved by 10V/R48 to: a 0.71-sigma flow
baseline against its own SD_block = 68.24, carried by UGC03580 (52%
of the null variance; drop -> flow flips to BE +7.2; drop + NGC5371
-> 69 galaxies prefer BE +27.5), with coverage the best-supported
axis at point grade and composition inseparable from coverage
(Spearman -0.89/-0.81).

PRIOR EXPOSURE: R48's GB verification already computed per-galaxy
contribution maps and permutation batteries ON THE FLOW LEG. 10W is
therefore a CONSOLIDATION instrument on that leg and a FIRST LOOK
only where registered below as such (the anchored-leg map, the joint
axis regression with conditioning, the resolution-cut refits = new
sky). Known carriers are named as known: UGC03580, NGC5371 (flow);
IC2574, NGC0891 (the 10U anchored-lead carriers). Discovery grammar
is restricted to what R48 did not already inspect.

NOTHING HERE CAN MOVE THE METER: band 65.4-70.4, stat 5.0, and the
credences 53/8 are out of scope by construction (SS8).

## 2. Engine

Exec-inheritance of calcs/stage10v_strat.py truncated at the
baselines marker (the flowboot/GB route): SPARC/CF4 parse, census,
FAMS, build_sub, fit_deep, contest_fit, sub_fields, sd_block —
BIT-VERBATIM. Margins ride the GB R-d realized-galaxy-block SD; NO
injection-calibrated bar anywhere (trap #28 moot by design). Contest
grade = contest_fit best-of-two, warm from parent optima. The CPU
pool (calcs/fitpool.py) may carry fits ONLY if data/fitpool_gate.txt
records GATE PASS; otherwise the stage runs serial and discloses.
Frozen conventions inherited and disclosed, not re-opened: bulge
Upsilon_bul = 0.7; global f_ML without prior; sv_swap UMa 2.4% note;
census-pending D564-8/D631-7 flagged wherever they appear.

## 3. Part 1 — the map (descriptive; fixed parent globals)

For each leg (anchored-78, flow-71): at the leg's own contest optima
per member (warm-replicated 10V baselines), the per-galaxy
contribution obj_g (residual block + its profiled per-galaxy
nuisance priors, the R48 indep_obj decomposition) and
Delta_g(m) = obj_g(m) - obj_g(BE) for m in {boot, p065, gm}.

Products (no bars, no letters): data/stage10w_map.csv — one row per
galaxy: name, leg, fgas, median log10 y, innermost-point angular
radius (arcsec), Delta_g for the three pairs; carrier shares
(fraction of sum|Delta| and of SD_block^2); top-5 carriers per leg.
This file is Proposal F's systematics-figure source.

TRAP-#29 COMPLIANCE (printed BEFORE any Part-2/3 statistic): each
leg's baseline d(BE,boot) graded against its own SD_block from this
run — the known flow 0.71 sigma restated, the anchored equivalent
printed. No Part-2/3 grammar may claim what these baselines cannot
support.

## 4. Part 2 — the joint axis instrument (the attribution question)

Response: Delta_g(boot vs BE) per galaxy, both legs (leg indicator
carried). Registered axes (4): composition (fgas), coverage (median
log10 y per galaxy), resolution (log10 innermost-point arcsec), leg
(anchored/flow).

PRIMARY per-axis statistic: Spearman rank correlation between
Delta_g and the axis (carrier-robust by construction). CO-READ:
the OLS slope (carrier-sensitive; collinearity disclosed — fgas vs
coverage is -0.89/-0.81 here, so joint-fit partials are DESCRIPTIVE
ONLY, never a finding).

Nulls (10,000 permutations, seed 707):
- coverage axis: unrestricted permutation of Delta_g across galaxies
  (within leg).
- composition, resolution, leg axes: permutation WITHIN coverage
  terciles (within leg; terciles on median log10 y) — the tested
  axis must carry information BEYOND coverage (the R48 conditioning
  lesson made structural).

Multiplicity: 4 registered tests; Bonferroni x4 stated with every p.

Grammar (C5-compliant):
- ESTABLISHED: Bonferroni-corrected stratified p <= 0.01 AND the
  Spearman sign+p (<= 0.05 corrected) survive every registered
  carrier drop (UGC03580; NGC5371; IC2574; NGC0891; the flow pair
  jointly).
- SUPPORTED-AT-POINT-GRADE: corrected p <= 0.05, not carrier-robust
  or not reaching 0.01.
- NOT ESTABLISHED: otherwise.
No credence content in any outcome (SS8).

## 5. Part 3 — the angular-resolution cut instrument (NEW SKY)

The R48 named live rival: flow innermost points median 15.0 arcsec
vs anchored 34.8; both flow carriers at ~8. Registered grid (no
tuning): theta_min in {10, 20, 30} arcsec. For each cut: drop every
point with angular radius < theta_min (angular radius = the R48
arcsec construction, gate-pinned below); galaxies keep >= 3 points
(window_world rule); dropped points/galaxies censused per leg.

Contests: 4 members x 2 legs x 3 cuts, contest grade, warm from
parent optima. At each cut: d(BE,boot), SD_block RECOMPUTED ON THE
CUT WORLD (the 10V bar-mismatch lesson: the frozen 68.24 is the
uncut flow object and is NOT a cut-world bar), winner, and the
anchored-lead co-read (does the anchored BE lead survive the same
cut — both-ways honesty).

Letters (registered; d = fun_boot - fun_BE, positive = BE better):
- W-RESOLUTION-CARRIED: on the flow leg, at BOTH theta_min = 20 and
  30, d changes sign (flow prefers BE) OR |d| < 1 x SD_block(cut);
  AND the anchored d stays within 1 x SD_block(cut) of its uncut
  value at those cuts. Grammar: "the flow boot-residue is
  RESOLUTION-CONSISTENT at point grade." The word ATTRIBUTED
  additionally requires the conditional co-requirement: a 200-rep
  paired galaxy bootstrap (seed 909, flowboot construction) on the
  20-arcsec cut flow leg with P(sign agrees with the point claim)
  >= 0.95 (the C5 bar). The boot fires ONLY if the point conditions
  are met.
- W-RESOLUTION-IMMUNE: flow d within 1 x SD_block(cut) of the uncut
  +(-48.7) at ALL three cuts -> the resolution rival is RETIRED as
  the carrier of the flow residue.
- W-MIXED: anything else — descriptive; feeds F as a measured
  sensitivity, no attribution word.
C15 harness rule: any box-pinned or one-sided-coverage stratum gets
a profile check before its parameters are printed.

## 6. Gates (sky-blind; each names what its failure vetoes — trap #11)

- G10W-1 engine identity: warm replication of the 10V baselines
  (anchored BE-best ordering and flow boot-best +48.7) to <= 0.05 in
  the deltas. FAIL -> vetoes everything; STOP.
- G10W-2 decomposition closure: sum_g obj_g = total objective to
  <= 1e-6, both legs, at BE and boot optima. FAIL -> vetoes Parts
  1+2 (and the stage: STOP — Part 3 shares the engine's trust).
- G10W-3 arcsec census: reproduce R48's innermost medians 15.0
  (flow) / 34.8 (anchored) exactly; spot-checks SPAN both legs
  including UMa members (trap #30: never a sorted prefix). FAIL ->
  vetoes Part 3 and the resolution axis; STOP.
- G10W-4 cut wiring: at theta_min = 0 the cut pipeline returns the
  uncut world bit-exactly (array equality, no fits); at 20 arcsec
  print the dropped census only (sky-blind: no cut FITS in gates).
  FAIL -> vetoes Part 3; STOP.
- G10W-5 null machinery: on a synthetic Delta with a planted
  coverage gradient (seed 1313), the unrestricted null detects it
  (p <= 0.001) and the coverage-stratified null does not (p >= 0.2).
  FAIL -> vetoes every Part-2 p-value; STOP.
- G10W-6 pool certificate: data/fitpool_gate.txt exists and records
  GATE PASS, plus an in-stage 2-task pool==serial spot re-check.
  FAIL -> science proceeds SERIAL (pool is convenience, never a
  premise); disclosed in the read.

## 7. Seeds, sizes, order of execution

Permutations seed 707 (10,000); synthetic-null seed 1313;
conditional cut-boot seed 909 (200 reps). Execution order: gates
(commit record) -> sky Part 1 (baseline grading first) -> Part 2 ->
Part 3 -> conditional boot if triggered -> letter -> Round 49.

## 8. Credence map (pre-signed) and never-quote

EVERY outcome of every letter above = HOLD 53 / HOLD 8. This stage
measures attribution structure and a rival; it cannot move the
meter (band 65.4-70.4, stat 5.0 untouched by construction) and it
carries no anomaly-level or mechanism-level content. Inherited
NEVER-QUOTE stands: flow_gasrich a0 6.27e-11 / "H0 ~ 40"; UMa-26
screening a0; windowed-leg a0s as meter readings; injection-bar
significances; "N ~ 25".

## 9. Outputs

calcs/stage10w_map.py (modes: gates | sky);
data/stage10w_gates.txt; data/stage10w_skyread.txt;
data/stage10w_map.csv (committed — F's figure source);
data/stage10w_verdict.txt after Round 49 adoption.
