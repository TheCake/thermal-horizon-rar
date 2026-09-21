# PRE-REGISTRATION — STAGE 10V: THE STRATIFICATION INSTRUMENT
2026-09-21. Written and committed BEFORE any counterfactual contest is
computed. Successor to 10U per the banked R47 fork (the cross-leg
inversion "is itself the sharpest description of the meter's form
systematic"); author's blanket go 2026-09-21 ("Do whatever you think
is needed").

## §0 Expectation disclosure (on record, feeds NO bar)
The anchored leg is dwarf/gas/deep-heavy; the flow leg is
spiral/transition-heavy; the two composition and coverage axes are
expected CORRELATED and possibly confounded, with the
distance-marginalization axis a live rival (flow galaxies carry
e_D 15–30% vs the anchors' 4–9%, and the model's per-galaxy dv priors
ride on exactly that width). No expectation enters any bar or letter.

## §1 The registered question
The identical four-member contest prefers BE on the anchored 78, boot
on the flow 71 (BE LAST, +48.7), gm on the union 149 (banked, GB-exact,
R47-C9). WHICH AXIS CARRIES THE INVERSION? Three registered suspects:
- H1 DISTANCE TREATMENT (error magnitude + its dv-marginalization),
- H2 COMPOSITION (gas dominance),
- H3 ACCELERATION COVERAGE (the y-range the leg samples).
This stage ATTRIBUTES; it moves NO meter quantity: the 10T/10U band
65.4–70.4, stat 5.0, and every banked clause are untouched on every
branch. Union contests are OUT OF SCOPE (no clause uses them).
KK98-251 / D564-8 / D631-7 stay out of scope (next census).

## §2 World, engine, inheritance (frozen)
- World = the 10U frozen census: W78 via make_world(reclass=True,
  lt_slot=None); legA = 78 (41 anchored-direct + 26 UMa + 11
  CF4-reclassified), flow = 71. NO membership changes.
- Engine = exec-inheritance of calcs/stage10u_nuform.py source
  truncated at the "G10U-2: member identity" marker (the R47-addendum
  route): SPARC/CF4 parse, census, make_world, FAMS, build_sub,
  fit_hier, fit_deep, contest_fit — BIT-VERBATIM, zero drift.
- Estimators: contest grade = contest_fit (best-of-two fit_deep,
  cold + warm); screening grade = warm-only fit_deep (labeled, no
  letter weight). Deltas printed at 0.1; any delta at or below its
  contest's max start gap is sub-convergence, not quotable (R46).
- NO injection-calibrated bars anywhere in this stage (trap #28 is
  moot by design): the only synthetic element is H1b's scatter ADDED
  TO THE REAL SKY (a data transform, not a generative mock); its
  refit s_int and lag-1 rho shifts are printed as disclosures.

## §3 Blindness inventory
PUBLIC PRE-KNOWN (inputs, not blind): the anchored contest table
(BE best; gm +24.0, p065 +25.9, boot +63.3) + per-member nuisances;
the GB warm-only flow deltas (boot best; gm +2.5, p065 +7.6, BE
+48.7); the anchored block SDs (BE→boot 22.96 own-construction);
rho = 0.62; carriers IC2574/NGC0891.
COMPUTED AT GATES, FROZEN PRE-SKY: the contest-grade flow baseline
(STOP-bearing, §5 G10V-2b); census values (rel medians, strata
counts, window bounds, a_conf, flow-side SD_block, flow baseline
nuisances); measured per-fit cost.
BLIND UNTIL SKY: every counterfactual contest, all margins and
reversal/material evaluations, strata and windowed winners, the flow
jackknife, the letter.

## §4 Statistics (registered)
- Winner of a contest = argmin profiled −2lnL over the four members.
- DECISIVE PAIR = (BE, boot) everywhere. Incumbent = that leg's sky
  winner (anchored: BE; flow: boot). Cross-winner = the other leg's
  sky winner.
- SD_block = the GB R-d realized-galaxy-block construction, verbatim:
  v_blk = Σ_g 4·(Σ_pts w·r·D)², D = mu_B − mu_A per point, w and r
  from member A's profiled fields; SD_block = sqrt(v_blk). Convention:
  A = that contest's better-fitting member of the decisive pair,
  fields at that contest's own optima.
- MARGIN RULE: a lead is margin-passing iff
  lead ≥ max(1 × SD_block, that contest's max start gap).
- REVERSAL (fires a clause): the counterfactual's overall winner
  equals the CROSS-winner AND its lead over the incumbent is
  margin-passing. A counterfactual won by gm or p065 is NOT a
  reversal (disclosed as winner-neither).
- MATERIAL (disclosure tier only, never fires): the decisive-pair
  delta moves from its leg's baseline by ≥ 2 × SD_block(baseline leg)
  without reversal.
- If the decisive pair ties within the start gap, margins are
  unquotable: clause no-fire, disclosed.
- Refit jackknife = standing harness part (R47): drop-one refits on
  the flow baseline's top-2 fixed-parameter carriers, (BE, boot) pair.
  Fixed-parameter decompositions are printed but labeled
  CONSTRUCTION-SENSITIVE (R47 R-f); the refit jackknife is primary.

## §5 Gates (sky-blind; each names what its failure vetoes)
- G10V-1 anchored regression: inherited engine re-produces the
  archived anchored contest — a0 within 0.1% of ARC per member;
  deltas within ±0.5 of (24.0/25.9/63.3); nuisances within R-b
  tolerances (f ±0.01, s_int ±0.002, u ±0.002).
  FAIL VETOES: everything (no letter; STOP).
- G10V-2a flow warm replication: warm-only fit_deep per member
  (warm = [log10 ARC, 1.0, 0.08, 0.0]) reproduces the GB flow deltas
  within ±1.0 each. FAIL VETOES: all flow-referencing clauses (STOP).
- G10V-2b flow contest-grade baseline: add cold fit_deep per member;
  best-of-two. PROCEED requires: boot best AND d(BE) − best ≥ 20.
  Otherwise STOP with letter V-BASELINE-DISSOLVED (the banked
  inversion was estimator-grade; that is itself the finding and goes
  to review). The 10V baseline deltas for all comparisons = these
  contest-grade values, superseding warm-only spots.
- G10V-3 wiring: (a) rescale at s = 1 leaves lgobs bit-identical;
  (b) at s = 1.1162 every rescaled point shifts by exactly
  −log10(1.1162); (c) sigv swap with rel_new = own rel reproduces
  sigv to 1e-12 (non-degenerate galaxies); (d) H1b at sigma = 0
  leaves the sub's lg bit-identical.
  FAIL VETOES: H1a/H1b/H1c/H1d clauses.
- G10V-5 block-SD regression: the generalized SD_block code on the
  anchored (BE→boot) pair reproduces 22.96 within 2%. The flow-side
  SD_block(decisive pair) is printed and FROZEN (new number, no bar).
  FAIL VETOES: all margin rules (letter capped at V-UNRESOLVED with
  point-delta disclosures only).
- G10V-6 census (values frozen at the gates print):
  rel_med_flow = median of W['rel'] over flow galaxies;
  rel_med_anch = median over legA; sigma_inj = rel_med_flow/LN10;
  fgas_g = median over the galaxy's points of gg/(gg+gd+gb) at
  f_ML = 1; gas-dominated iff fgas_g ≥ 0.5; strata counts for
  flow-gasrich/flow-gaspoor/anch-gasrich/anch-gaspoor; the y-window
  (§6 H3) bounds and surviving-galaxy counts; a_conf (§6). MIN-N
  RULE: any letter-feeding stratum or windowed leg with < 15
  galaxies renders that clause UNPOWERED (auto-no-fire, disclosed).
  FAIL of the census print itself VETOES H2/H3 clauses.

## §6 The counterfactual battery (all constants frozen here)
H1 DISTANCE — four probes; sigv swap formula: incl_g =
sqrt(max(sigv_g² − (rel_g/LN10)², 0)); sigv_new = max(hypot(
rel_new/LN10, incl_g), 0.01).
- H1c anchored-loosened (PRIMARY): all 78 legA galaxies get
  rel_new = rel_med_flow; data unchanged; contest.
- H1b anchored + flow-grade scatter (PRIMARY): per-galaxy
  delta_g ~ N(0, sigma_inj) via rng = default_rng(seed), drawn in
  LEGA order, applied via dlg_per_occ. Variant (ii) [primary]:
  priors ALSO loosened as H1c; seeds (7001, 7002, 7003, 7004).
  Variant (i) [disclosure]: priors unchanged (understated-error
  world); seeds (7001, 7002).
- H1a flow-rescale (CO-READ): lgobs[flow points] −= log10 s,
  s = 73/65.4 = 1.1162 (co-read s = 73/70.4 = 1.0369); sigv
  unchanged; contest on flow.
- H1d flow-tightened (CO-READ): flow galaxies get rel_new =
  rel_med_anch; data unchanged; contest.
H1 FIRES iff H1c reverses, OR H1b(ii) reverses in ≥ 3 of 4 seeds
(margin rule each seed). H1a/H1d/H1b(i) can reach MATERIAL only.

H2 COMPOSITION — contests on the four strata (galaxies keep all
their points). H2 FIRES iff, on the FLOW leg: the gas-rich stratum's
winner = BE with a margin-passing lead over boot AND the gas-poor
stratum's winner = boot with a margin-passing lead over BE (the
inversion reproduces INSIDE the flow leg along composition).
Anchored strata = mirrored co-read (disclosure only).

H3 COVERAGE — y_i = gN1_i/ARC_BE with gN1 = gg + 1.0·gd + gb;
per-leg point distributions of log10 y; window = [max(P10_A, P10_F),
min(P90_A, P90_F)] (percentiles over points); keep in-window points;
drop galaxies with < 3 surviving points; contest each windowed leg.
H3 FIRES iff the two windowed legs AGREE on the overall winner and
each windowed leg's winner has a margin-passing lead over the other
decisive-pair member (the inversion DISSOLVES under matched
coverage).

CONFOUND METRIC: on the flow leg, a_conf = fraction of galaxies
where the labels (gas-poor) and (median point log10 y > window
midpoint) agree; a_conf ≥ 2/3 adds the suffix -CONFOUNDED to any
letter in which H2 or H3 participates.

SCREENING (warm-only, no letter weight): sub-block contests on
anchored-direct 41 / UMa 26 / reclassified 11; per-member f_ML and
s_int tables for every counterfactual; lag-1 rho per counterfactual;
H1b refit s_int (the trap-#28 vetting print).

## §7 Letters (grammar; positive words gated)
- V-DISTANCE-CARRIED: H1 fires, H2 and H3 do not.
- V-COMPOSITION-CARRIED: H2 fires alone.
- V-COVERAGE-CARRIED: H3 fires alone.
- V-MIXED(list): two or more fire (suffix -CONFOUNDED per §6).
- V-UNRESOLVED: none fires — the inversion stands measured with no
  attributed carrier; MATERIAL rows and winner-neither rows disclosed.
- V-BASELINE-DISSOLVED: the G10V-2b STOP branch (gates-stage letter).
MANDATORY DISCLOSURES on every branch: multiplicity (4 fire
opportunities: H1c, H1b(ii), H2, H3 — attribution-grade, not
demotion-grade); construction-sensitivity labels; the
no-meter-quantity-moves pin; nuisance tables; rho co-reads.

## §8 Credence map (pre-signed)
EVERY branch of §7 → HOLD anomaly-real 53 / mech 8. No exceptions;
attribution of an instrument systematic moves no world-credence cell.

## §9 Verification & review
GA blind half committed BEFORE the review round (15th protocol
execution); ROUND 48 = fresh reviewer session; GB re-computes every
load-bearing reviewer number in calcs/round48_addendum.py before
adoption (standing memory rule). Amendments: instrument-side only,
pre-quote, runs preserved.

## §10 Successor mandates
Whichever axis fires (or UNRESOLVED): (a) the BIG-SPARC refire
inherits the finding as its stratification spec alongside the
R47-corrected arm; (b) Proposal F's systematics section quotes the
letter verbatim; (c) if H2 fires, the dissident-dwarfs stage (
Proposal D) inherits the strata tables as its starting census.

## §11 Runtime (honest, from measured cost)
Per-fit cost is MEASURED at gates and the sky ETA quoted from it
(the 10U lesson: never estimate from a cheaper arm). Registered fit
budget: gates ≈ 14 fits + census; sky ≈ 144 fits (H1c 8, H1d 8,
H1a 16, H1b(ii) 32, H1b(i) 16, H2 32, H3 16, screening 12 warm-only,
jackknife 4). Heartbeat prints in every loop over 10+ fits.
