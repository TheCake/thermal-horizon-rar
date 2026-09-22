# PRE-REGISTRATION — STAGE 10Y: THE RESIDUAL-CONTRAST ORGANIZER

Status: committed BEFORE any new sky statistic is computed. The only
residual-derived numbers consulted in writing this document are the
published Round-49 values (data/stage10w_verdict.txt C17;
data/round49_addendum.txt GB-9). Bars, branch grammar, constants,
seeds and gate designs are locked here. Sky legs run only after this
file and the gates output are committed.

## 1. The question and the object

Round 49 banked a direction-finder (C17, GB-confirmed): at each
arm's own parent optimum, the galaxy-paired inner-minus-outer
residual contrast at a 20-arcsec split is

    anchored:  -0.0857 +/- 0.0301 dex over 20 galaxies  (2.9 sigma)
    flow:      -0.0055 +/- 0.0176 dex over 45 galaxies  (0.2 sigma)

(the four-decimal values are the GB own-stream prints, the operative
regression targets). The registered successor (C17: "no fits
needed") must break the direction-finder's three named caveats:
y-composition of inner points, within-curve autocorrelation
(rho_lag1 0.61-0.72), and the misfit-vs-smearing degeneracy.

THE QUESTION: what ORGANIZES the anchored inner residual depression?
Three registered candidates:

  (a) ACCELERATION y — inner points sample different y than outer;
      the depression is the parent form's y-shape misfit wearing a
      radius costume (the form-band story; no new axis).
  (b) ANGULAR RADIUS (arcsec) — beam smearing / finite resolution
      suppresses inner velocities (the 10V named live rival;
      instrument axis). Smearing prediction: the effect organizes by
      ANGLE, not by physical or disk-scaled radius, and worsens for
      galaxies whose inner points sit at smaller angles (and, at
      fixed physical radius, for more distant galaxies).
  (c) DISK-SCALED RADIUS R/R_d — inner-disk astrophysics
      (non-circular motions, streaming, bars) organized by the
      disk's own scale. Frozen-record precedent (P2's 7B/7C arc,
      PAPER.md): SPARC point scatter at fixed x is inner-disk
      organized, ~2.4x above outer inside R < 1.5 R_d, and the x~1
      excess vanishes on outer points (b_clean = 0.0000). Boundary
      1.5 R_d is INHERITED from that record, not tuned here.

The candidates are collinear within a curve (y falls with R; theta
and R/R_d both grow with R). The instrument therefore reads the
DISAGREEMENT cells, where the organizers separate.

## 2. Data, parents, and the no-new-fit rule

World: the 10V frozen census (W78; LEGA 78 anchored, FLOW 71 flow),
via exec-inherit of calcs/stage10v_strat.py truncated at the
"baselines (both modes)" marker — the round49_gb.py pattern, copied
verbatim. Parent optima: the GB construction copied verbatim
(anchored: plain -> loose hier -> fit_deep BE -> warm chain to
p065/gm/boot; flow: fit_deep from the standard th0 per member).
These re-derive published archives (regression grade). The stage
fits NOTHING else; every new statistic is a function of per-point
residuals RES (= lg - mu with per-galaxy dml, dv profiled at the
parent globals — the sub_fields object, extended verbatim to also
return per-point y = gN/a0 at the same optimum), per-point angles
theta (rotmod R/D x 206.265, the GB-9 census cuts verbatim),
per-point physical radii R_kpc (same parse), and per-galaxy R_d
(SPARC_Lelli2016c.mrt bytes 62-66, matched by name).

Per leg, the residual field is evaluated at the leg's own best form
(BE vs boot by engine fun at the parent optima — the GB-9 A2 rule,
copied verbatim). ANCHORED is the primary arm (the depression lives
there); flow is a co-read at its banked 0.71-sigma-baseline grade
(R48 C3) — no flow letter clause.

Gate/sky separation rule: the gates may use (i) published numbers as
regression targets, (ii) the DESIGN at the published parent optima
(point census, y, theta, R_kpc, R_d, galaxy memberships), and (iii)
SYNTHETIC residuals; they may not compute or print any new statistic
of the real residual values. Sky legs run post-commit only.

## 3. Statistics (all constants locked here)

Constants: theta* = 20 arcsec (inherited from C17; not re-tuned).
Disk boundary R* = 1.5 R_d (inherited, section 1c). Seed rng =
numpy default_rng(20260923) everywhere. Permutations: 4000.
Circular-shift null draws: 2000 (size), 1000 (branch-rate tables).
Synthetic noise: per-galaxy AR(1), rho = 0.66 (the C17 mid),
sigma_pt = 0.10 dex (the frozen-record point-floor scale). Holm
family for S2: {S2a, S2b}. Contrast sign convention: inner minus
outer; the published depression is NEGATIVE.

S0 IDENTITY (gate G10Y-1, pre-sky): reproduce GB-9 on both legs to
|d(mean)| <= 1.5e-4, |d(SE)| <= 1.5e-4, galaxy counts exactly 20/45.

S1 THE Y-MATCHED CONTRAST (the caveat-1 breaker; primary
construction = overlap matching, assumption-free): per galaxy with
points on both sides of theta*, the y-overlap window is
[max(min y_in, min y_out), min(max y_in, max y_out)]; keep points
inside the window on both sides; require >= 1 point each side;
contrast c_g = mean RES(in) - mean RES(out) over kept points.
Report mean +/- SE over qualifying galaxies (the GB-9 convention:
one number per galaxy; between-galaxy SE) + the overlap census.
N_min(S1) = 8 qualifying galaxies; below it S1 is UNPOPULATED and
branch B1 cannot fire (grammar falls through).
Declared-bias co-read (no branch weight): pooled detrend — fit
RES vs log10 y (linear / quadratic / cubic, three co-reads) over
the leg's points, subtract, recompute the theta* contrast. The
detrend is biased TOWARD collapse (y and radius are collinear
within curves; the trend can absorb a genuine radius effect) — it
is quoted with that label and cannot fire a branch. If the
quadratic detrend and the overlap-match disagree in verdict class,
the stage says so in the letter parenthesis.

S2 THE ACROSS-GALAXY ANGULAR ORGANIZER (smearing's between-galaxy
prediction): with c_g = the UNMATCHED per-galaxy contrasts (the
GB-9 diffs) on the anchored leg,
  S2a: Spearman rank correlation of c_g with theta_min,g (the
       galaxy's innermost kept angle). Smearing sign: POSITIVE
       (more negative contrast at smaller angle).
  S2b: Spearman of c_g with distance D_g. Smearing sign: NEGATIVE
       (more negative contrast farther away).
p by galaxy-label permutation (4000, two-sided), Holm over the
pair; the smearing SIGN is required for a fire. N_min(S2) = 15
galaxies with defined c_g.

S3 THE DISAGREEMENT CELLS (the caveat-3 breaker; within-galaxy):
classify each point of a qualifying galaxy by the two organizers:
I/O by theta (< / >= theta*) and i/o by disk radius (< / >= R*).
Cells: Ii (inner-both), Io (inner-angle only), Oi (inner-disk
only), Oo (outer-both). Per galaxy, cell deltas vs the Oo
baseline: d_Io = mean RES(Io) - mean RES(Oo), d_Oi likewise.
Leg-level: mean +/- SE over galaxies contributing the cell.
N_min(cell) = 5 galaxies; an under-populated cell reads
UNPOPULATED and its signature cannot fire.
  ANGULAR signature: d_Io <= -2 SE(d_Io) AND |d_Oi| < 1 SE(d_Oi).
  DISK signature:    d_Oi <= -2 SE(d_Oi) AND |d_Io| < 1 SE(d_Io).
(Mutually exclusive by construction.) Alignment p for each fired
signature: within-galaxy CIRCULAR-SHIFT null — per draw, each
galaxy's radius-ordered residual sequence is circularly shifted by
an independent uniform offset (preserves the curve's
autocorrelation, breaks alignment with the radius classes; the
caveat-2 answer), 2000 draws; report P(null as extreme).

S4 DESCRIPTIVE PRODUCTS (no branch weight, banked for Proposal F):
the banded residual profile (theta bands 0-10/10-20/20-40/40+ x
R/R_d bands <1/1-1.5/1.5-3/3+ x y terciles; mean +/- SE + counts);
rho_lag1 restated per leg; the flow-leg replication of S1-S3; the
carrier note (UGC03580, NGC5371 innermost angles ~8 arcsec).

## 4. Letter grammar (trap #31 compliant)

Evaluated in this fixed order (if / elif; mutually exclusive by
construction). "Survives" = S1 overlap-match mean <= -2 SE (sign
required); "collapses" = |S1 mean| < 1 SE.

  B1 Y-CARRIED: S1 populated AND collapses.
     Reading: the depression is y-composition riding the parent
     form; no radius-specific axis; the resolution rival gains
     nothing from C17. FALSE in-grid if |S1| >= 1 SE.
  B2 RADIUS-SPECIFIC-ANGULAR: S1 populated AND survives AND S3
     fires the ANGULAR signature (populated cells) AND S2 Holm
     fires (p <= 0.05, smearing sign, either member). Reading: the
     resolution rival gains a measured footprint on the anchored
     arm (still bounded for the meter — section 6). FALSE in-grid
     if any conjunct fails.
  B3 RADIUS-SPECIFIC-DISKSCALE: S1 populated AND survives AND S3
     fires the DISK signature AND S2 does not fire in smearing
     sign. Reading: the depression is inner-disk astrophysics on
     the disk's own scale (the 7B/7C object on the anchored arm);
     the resolution rival loses its remaining anchored foothold.
     FALSE in-grid if any conjunct fails.
  B4 UNRESOLVED-DEGENERATE: everything else — S1 in [1,2) SE, S1
     unpopulated, signatures unpopulated or not fired, S3 fired
     without its required S2 state, or conflicting axes. Named
     sub-labels printed for the record (which conjunct failed),
     letter weight none.

Tie-break: none needed — B2/B3's S3 signatures are exclusive by
construction and their S2 conditions are complementary; any
configuration outside the three defined conjunctions IS B4.

Branch priors (trap #33, recorded before gates): B1 0.25, B2 0.15,
B3 0.20, B4 0.40. All >= 10% — every branch stays a letter clause.

P(fire | synthetic world) per branch is COMPUTED AT GATES (G10Y-3)
on both registered worlds and printed before the sky runs (trap
#32). If a branch's gate power fails (section 5), that branch is
demoted to descriptive pre-sky and the letter says so.

## 5. Gates (all pre-sky; failures stop the stage)

  G10Y-1 IDENTITY: S0 as section 3 (the GB-9 regression).
  G10Y-2 R_d CENSUS: every LEGA + FLOW member matched by name to
     SPARC_Lelli2016c.mrt with R_d > 0; print three spot values
     (IC2574, DDO154, NGC2403) for the record; print the theta*
     <-> 1.5 R_d crossing census (which cells CAN populate).
  G10Y-3 NULL MACHINERY + BRANCH-RATE TABLES: two synthetic worlds
     on the real design (real point census, y, theta, R_kpc, R_d;
     synthetic residuals; AR(1) rho 0.66, sigma_pt 0.10):
       WORLD-Y (pure y-carried): residuals = b x log10 y + noise,
       with b solved so the composition-induced theta* contrast on
       the anchored design equals -0.0857 (the published value).
       WORLD-R (pure radius-carried): residuals = noise + a flat
       -0.0857 step on inner-by-theta* points.
     (i) SIZE: on noise-only draws (2000), each firing statistic's
     P(fire) printed; bars: S3 signatures and S2 Holm each within
     [0.01, 0.10] nominal (the C15-ii standard).
     (ii) BRANCH RATES (1000 draws/world): the full grammar
     evaluated per draw; the trap-#32 table printed. POWER BARS:
     P(B1 | WORLD-Y) >= 0.60 and P(B2 or B3 or B4-with-S1-surviving
     | WORLD-R) >= 0.60, i.e. the instrument must usually collapse
     a pure-y world at B1 and usually keep a pure-radius world's S1
     alive. A failed power bar demotes the dependent branch(es) to
     descriptive, stated in the gates output before the sky.
  G10Y-4 Y-MATCH WIRING: on WORLD-Y, the S1 overlap-match mean must
     read within 1 SE of zero in >= 80% of 500 draws (the matching
     actually removes a pure y-trend); on WORLD-R, the S1 mean must
     stay <= -1 SE in >= 60% of 500 draws.

## 6. Meter impact and guardrails

No branch moves the meter. H0 = 65.4 +/- 5.0 (stat) +/- 0.5, the
form band 65.4-70.4, and the validity domain 62-75 are untouched by
construction (no fit is run off the parent optima). The stage's
maximum conceivable meter exposure under B2 is ALREADY BANKED as the
10W C14 cut lever (anchored H0 moves at most -0.8 = 0.16 sigma_stat
under inner cuts); B2 would change the LABEL on part of that
sensitivity, not its size. No new H0 number in this stage is a
reading. Inherited never-quote strings remain (cut-flow a0 values;
"W-RESOLUTION-CARRIED"; "winner = gm at 20"; the Part-2 leg-axis p
bare; "paired t = +4.27"). Flow-leg sentences carry the R48 C3
baseline grade (0.71 sigma) whenever quoted.

Credence map: EVERY cell of this stage pre-signs HOLD (anomaly-real
53 / mechanism 8). An attribution instrument on a residual feature
moves no world-level credence regardless of branch.

Amendments: instrument-side only, each committed pre-quote with the
prior run preserved under a suffixed filename (house standard).

## A1 — AMENDMENT (pre-quote; gates run 1 preserved as
## data/stage10y_gates_r1.txt; committed BEFORE the amended gates run)

Run 1 (G10Y-1 PASS, everything else FAIL) found one bug and two
DESIGN-DEAD constructions. No sky statistic was touched. The fixes:

A1-i PARSER: the .mrt is not at the header's byte layout on disk;
  Rdisk is parsed by whitespace token (index 11), name = token 0.
  Spot truths: IC2574 2.78, DDO154 0.37, NGC2403 1.39 kpc.

A1-ii S1 REDESIGN (the per-galaxy y-overlap match qualifies ZERO
  galaxies — y is near-monotone in radius, so inner and outer
  y-ranges barely overlap within one galaxy; run-1 print):
  S1 primary = the Y-BANDED matched contrast. Band edges = quintiles
  of log10 y over the pooled points of both-sides galaxies (DESIGN
  quantities, fixed once). Per band with >= 5 inner and >= 5 outer
  points: c_b = mean RES(in, b) - mean RES(out, b); S1 = the
  count-weighted mean (w_b = n_in n_out/(n_in + n_out)); SE = the
  leave-one-galaxy-out JACKKNIFE over the both-sides galaxies
  (block-honest, deterministic). Populated: >= 3 usable bands AND
  >= 8 both-sides galaxies. Survive/collapse bars unchanged
  (-2 SE / |.| < 1 SE). The per-galaxy overlap match demotes to a
  census-reported co-read; the pooled detrends stay labeled
  co-reads. Rationale: within a y band, inner and outer points
  carry (nearly) the same y, so a pure y-trend cancels; a
  radius-specific depression does not.

A1-iii S3 REDESIGN (the Io cell is UNPOPULATABLE on the anchored
  arm — 20 arcsec at anchored distances is a small physical radius,
  so no point is inner-by-angle but outer-by-disk; run-1 census
  Io = 0, Oi = 20): the populated grammar reads THREE cells:
  Ii (th < 20 AND R < 1.5 Rd), Oi (th >= 20 AND R < 1.5 Rd),
  Oo (th >= 20 AND R >= 1.5 Rd). Per galaxy with Oo populated:
  d_Ii and d_Oi vs the Oo baseline. The organizers separate on
  d_Oi — the points BEYOND the angular split but INSIDE the disk
  scale:
    ANGULAR signature: d_Ii <= -2 SE(d_Ii) AND |d_Oi| < 1 SE(d_Oi)
      (the depression stops at the angle).
    DISK signature: d_Ii <= -2 SE(d_Ii) AND d_Oi <= -2 SE(d_Oi)
      (the depression fills the disk scale).
  Mutually exclusive via the d_Oi clause; the gap (d_Oi between
  -2 SE and -1 SE, or d_Ii shallow) is B4. Io is reported as a
  census line only. Named confound: Oi sits at lower y than Ii; the
  y-banded d_Oi co-read is printed (descriptive). N_min(cell) = 5
  unchanged; circular-shift p unchanged, now for d_Ii/d_Oi.

A1-iv SIZE BARS: the [0.01, 0.10] band applies to each CONSTITUENT
  one-sided 2-SE test (d_Ii, d_Oi, S2 Holm+sign); the compound
  signature rates are PRINTED without bars (a conjunction is
  legitimately rarer than 1% under noise; the run-1 bar was
  mis-posed for conjunctions).

A1-v POWER/WIRING: same bars (P(B1|WORLD-Y) >= 0.60;
  P(radius kept|WORLD-R) >= 0.60; wiring >= 0.80 / >= 0.60), now
  evaluated on the banded S1 and the A1-iii grammar.

Branch priors unchanged (B1 0.25 / B2 0.15 / B3 0.20 / B4 0.40).
Nothing in this amendment reads a real residual number beyond the
published targets and the run-1 gate prints.

## A2 — AMENDMENT (pre-quote; gates run 2 preserved as
## data/stage10y_gates_r2.txt; committed BEFORE the amended gates run)

Run 2: G10Y-1/2 PASS (parser fixed; censuses healthy: Ii 19, Oi 55,
Oo 77 galaxies; constituent sizes 0.031/0.028/0.018 all in-bar), but
the A1-ii banded S1 is DESIGN-DEAD too: the 20 both-sides anchored
galaxies hold only 27 inner points total (band counts 3/4/4/4/12
against 75/73/73/73/66 outer), so at the >= 5-per-band rule a single
band qualifies, deterministically. A y-banding fine enough to cancel
the trend cannot be populated at 27 points. Design fact banked: the
C17 depression is carried by 27 sub-20-arcsec points.

A2-i S1 REDESIGN (the caveat-1 breaker, third form): point-level
  NEAREST-NEIGHBOR y-matching. Pool the qual (both-sides) galaxies'
  points; for each INNER point a, its window = all OUTER points
  (any qual galaxy) with |log10 y - log10 y_a| <= 0.15 dex; a
  qualifies if the window holds >= 5 points. S1 = mean over
  qualifying inner points of (RES_a - mean RES(window)); SE = the
  leave-one-galaxy-out jackknife over the union of contributing
  galaxies (a left-out galaxy's inner points leave the average AND
  its outer points leave every window; qualification is
  re-evaluated per leave-out). Populated: >= 15 qualifying inner
  points from >= 8 galaxies. Windows are design objects (y is
  design), precomputed once. Bars unchanged (survive <= -2 SE;
  collapse |.| < 1 SE). Residual leak under a pure trend is
  bounded by |b| x 0.15 ~ 0.025 with near-cancelling sign across
  points — G10Y-4 measures it. Cross-galaxy matching is licensed
  by the per-galaxy dv profiling already removing galaxy offsets
  (stated limitation: any residual galaxy-level structure enters
  both sides of the match).
  The A1-ii banded form demotes to a printed co-read (no branch
  weight); the per-galaxy overlap census stays a census line.

A2-ii The same matcher, applied to Oi points vs Oo windows, is the
  y-matched d_Oi co-read named in A1-iii (descriptive).

A2-iii Constants: dy = 0.15 dex, window minimum 5 points, floors
  15 points / 8 galaxies. Everything else unchanged, including the
  A1-iii cell grammar, bars, priors, and worlds.

## A3 — AMENDMENT (pre-quote; gates run 3 preserved as
## data/stage10y_gates_r3.txt; committed BEFORE the amended gates run)

Run 3: everything passes except the G10Y-4 WORLD-Y clause, and that
bar was ARITHMETICALLY IMPOSSIBLE as posed: it demanded
P(|v| < 1 SE) >= 0.80, but that event has ~0.68 probability for an
exactly unbiased estimator (the run-3 noise-only context prints
0.669; WORLD-Y prints 0.624 — the matcher removes the trend to
within Monte-Carlo resolution). Bar-design lesson for the record
(the A1-iv class, second occurrence in this stage): a gate bar is a
fire-event too — compute its rate under the IDEAL instrument before
signing it.

A3-i G10Y-4 WORLD-Y clause re-posed as a BIAS gate: over 500
  WORLD-Y draws, |mean of the matched value| <= 0.010 dex
  (one-tenth of the target depression), with the Monte-Carlo SE
  printed. The WORLD-R keep clause is unchanged (run-3 value 0.970,
  PASS). The branch-rate power bars (G10Y-3) already carry the
  "usually collapses" requirement and passed (0.625 / 0.968).

## 7. Outputs

calcs/stage10y_organizer.py (modes: gates / sky) ->
data/stage10y_gates.txt, data/stage10y_skyread.txt; verdict file
data/stage10y_verdict.txt after the review round; ledger row
cos-10y-organizer. Review protocol: GA blind pre-committed before
the round report is read; fresh-session review round; GB
re-computation of every load-bearing reviewer number (standing
rules).

## POST-ROUND ANNOTATION (2026-09-23, Round 52 conditions 9/10;
## append-only — the registered text above is not edited)

(i) A1-ii's "qualifies ZERO galaxies ... run-1 print" is corrected:
the per-galaxy y-overlap match qualifies ONE anchored galaxy
(ESO563-G021, contrast -0.1379, GB-13 exact), and no overlap count
appears in data/stage10y_gates_r1.txt — the actual run-1 evidence is
that the run-1 construction returned NaN below 2 qualifying galaxies,
driving every run-1 wiring and branch rate to 0.000. The conclusion
(1 << N_min = 8, construction design-dead) is unaffected.

(ii) A1-iii's Io = 0 census was decided on the RUN-1 (broken-parse)
Rdisk values; the bug inflated Rdisk and therefore biased TOWARD
Io = 0. It was re-confirmed on the corrected run-2 census before the
sky ran, and Round 52 confirmed it structural (every sub-20-arcsec
anchored point sits at R/Rd 0.078-1.089, median 0.485). Conclusion
sound; the decision order was not, and is recorded here.

(iii) Section 5's WORLD-R power clause says "B4-with-S1-surviving"
(text: survive = -2 SE) while the code counts S1 <= -1 SE. Measured
both ways (R52 + GB-12): 0.815-0.863 at -2 SE, 0.968-0.970 at -1 SE;
the >= 0.60 bar passes under either reading. The coded -1 SE form is
the as-fired one.

(iv) The registered letter string for branch B4 is superseded by the
Round-52 replacement letter in data/stage10y_verdict.txt
(Y-EXCLUDED-ANGLE-FAVOURED ...); "B4" remains the grammar's branch
output (metadata), never the stage summary.
