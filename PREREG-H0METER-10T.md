# PRE-REGISTRATION — stage 10T "LEG-A REGROW" (2026-09-20 night)
# Status: to be COMMITTED after the sky-blind instrument build passes its
# gates and BEFORE any extended-set sky central is computed (the author's
# blanket go for the hammering queue, 2026-09-20: 10T is item 1).
# Parent: PREREG-H0METER-DRAFT.md (stage 10S, commits 2a461d5/fb74226) —
# all 10S pins 1–10 are INHERITED verbatim unless amended below. The 10S
# verdict (data/stage10s_verdict.txt, Round 45) stands: leg A is the
# durable product; leg B is POWER-LIMITED at SPARC grade and IS NOT
# REFIRED here.

## Question
The 10S meter's leg A read H0_A = 65.87 ± 12.52 (hier BE primary; family
65.9–70.9) on 67 galaxies. CosmicFlows-4 (Tully+23, J/ApJ/944/94) now
supplies anchor-grade distances for a set of SPARC flow galaxies, and
LITTLE THINGS (Oh+15, 9R-validated port) supplies independent anchored
dwarfs. Does the regrown anchored leg (expected 82 galaxies, the tightest
e_D members replacing 15–30% flow errors with 4–9% TRGB errors) sharpen
H0_A, and do the additions integrate (composition splits inside 1σ) or
split (a named systematic)?

## Distance-source pins
- CF4 source = table2 (per-galaxy) PER-METHOD moduli ONLY. The CF4
  combined DM and the TF/FP columns are EXCLUDED: TF/FP distances presume
  a luminosity–kinematics relation — dynamically circular for a RAR
  instrument (and the combined DM mixes them in). Admissible methods:
  DMmas, DMceph, DMtrgb, DMsbf, DMsnIa, DMsnII.
- Provenance (CF4 ReadMe Note 4, verbatim): per-method moduli are given
  "after registration to a common scale with the MCMC analysis"; the
  common scale is set by "Cepheid period-luminosity relation and tip of
  the red giant branch observations founded on local stellar parallax
  measurements along with the geometric maser distance to NGC 4258."
  These moduli are H0-FREE but PEG-SHARED with the ladder (the
  three-world gearing registered in PREDICTIONS.md §F this commit). CF4's
  own H0 = 74.6 is an output of their flow analysis and is NOT imported.
- Method priority when several exist: mas > ceph > trgb > sbf > snia >
  snII (geometric, then primary stellar rungs, then anchored secondary
  rungs). Adopted D = 10^((DM−25)/5) Mpc; e_D = D·(ln10/5)·e_DM.
- PGC keying: data/scout10t_pgc_cache.json (SIMBAD, canonical alias
  "LEDA n"), frozen this commit after the alias round (CamB resolved
  PGC 166084; WLM/Haro29/Haro36 resolved; F564-V3 unresolved but
  3.6µm-flagged = excluded regardless). Still-unresolved SPARC names
  (20: the F5xx/D5xx LSB set + KK98-251) are DISCLOSED census
  incompleteness; all are flow LSBs at ≳45 Mpc, beyond TRGB reach —
  expected CF4 anchor yield nil. The stage runs OFFLINE from the frozen
  cache (deterministic; no network at stage time).

## Membership rules (R1–R8, pinned before any extended-set fit)
- R1 (SPARC reclassification): a SPARC f_D=1 galaxy in the 10S harness
  kept-set moves flow → leg A iff its cached PGC has ≥1 admissible CF4
  method. Census (frozen 2026-09-20, harness-derived): 12 reclassify on
  paper — ESO563-G021 (snII, 60.8→68.9), NGC2903 (trgb, 6.6→8.95),
  NGC4559 (trgb, 9.0→8.63), NGC5585 (trgb, 7.06→6.82), UGC02487 (sbf,
  69.1→52.5), UGC05918 (trgb, 7.66→8.32), UGC05986 (trgb, 8.63→8.91),
  UGC07323 (trgb, 8.0→7.45), UGC07603 (trgb, 4.7→6.46), UGC08550 (trgb,
  6.7→8.75), UGC09992 (trgb, 10.7→11.02), UGC12632 (trgb, 9.77→8.83).
  UGC09992 has ZERO points surviving the verbatim eV/V ≤ 0.10 cut and
  therefore enters no fit (11 effective). UGC02455 (scout hit) fails the
  10S quality cuts and stays OUT (R6).
- R2 (rescale): reclassified members' points rebuilt by the 10S pin-3
  scaling laws at s = D_CF4/D_SPARC: lgobs → lgobs − log10(s), g_bar
  invariant; vertical prior sv rebuilt with CF4 e_D.
- R3 (leg B): the reclassified members LEAVE the flow set (82 → 71
  with-points). Leg B is not refired (R45 standing); the flow census
  change is disclosure-only.
- R4 (UMa): the 26 UMa members KEEP the 10S treatment (shared-u nuisance,
  18 ± 0.9 Mpc block, per-galaxy 2.3 Mpc depth) in the PRIMARY. The two
  CF4-Cepheid members are a pre-registered CONSISTENCY CO-READ, frozen
  here from the census: NGC3972 DMceph 31.635 → 21.23 Mpc (s = 1.180 vs
  the block), NGC4051 DMceph 31.103 → 16.62 Mpc (s = 0.923). The co-read
  compares the fitted shared u against these two (disclosure; N=2 inside
  a 2.3 Mpc depth — annotation grade). Family variant A-var2 moves both
  to direct anchors at CF4 Cepheid D (u-group 26 → 24).
- R5 (LITTLE THINGS): an LT galaxy enters leg A iff ALL of: not
  SPARC-overlap (PGC-keyed), not 3.6µm-flagged (t2 f_alphamin/f_alpha3.6
  — the released DM curve absorbs the stars, port invalid; the 9R rule),
  admissible CF4 method exists, i ≥ 30, and ≥3 points survive the
  verbatim eV/V ≤ 0.10 cut after the 9R port (V_bar² = V_tot² − V_DM²,
  ring tolerance max(0.005 kpc, 1% R), V_bar² ≤ 0 dropped) at Hunter+12
  D rescaled to CF4 D by pin-3. Census (frozen): 4 enter — CVnIdwA
  (trgb 3.8 Mpc, 5 pts), DDO_52 (trgb 9.6, 12), DDO_133 (trgb 4.8, 17),
  DDO_210 (trgb 1.0, 5); all four genuinely gas-dominated (Oh+15 t2
  Mgas > Mstar). Excluded with reasons printed: no-CF4-anchor (DDO_70,
  DDO_101, NGC_1569, UGC_8508, HARO_36), flagged (DDO_43/46/47,
  F564-V3, HARO_29), i < 30 (DDO_53, NGC_3738, DDO_46), zero/thin points
  at the cut (IC_10, IC_1613, DDO_216), overlap (7 galaxies).
- R6 (no quality bypass): SPARC galaxies failing the 10S cuts (inc < 30
  or Q > 2) never enter, whatever their CF4 status.
- R7 (LT likelihood wiring): LT members enter with the lumped ported
  baryons in the GAS SLOT (g_gas = V_bar²/R, g_dsk = g_bul = 0): the
  M/L nuisance (f, dml) is inert on them BY CONSTRUCTION — Oh+15's own
  stellar mass stands, physically motivated for gas-dominated dwarfs.
  Family variant A-var3 puts the lumped term in the DISK slot (f + dml
  active on the total; the 9R G9R-W-licensed wiring).
- R8 (overlap galaxies): SPARC∩LT galaxies contribute ONLY their SPARC
  points (native decomposition beats the port); the LT port for them is
  the G10T-3 cross-validation instrument, nothing else.

## Primary and family
PRIMARY: hier BE (verbatim 10S machinery: 5M minus lensing, UMa u
active, U_PRIOR (0.9/18)/ln10), extended membership (82 = 41 anchored +
26 UMa + 11 reclass + 4 LT), CF4 distances on new members only,
LT-in-gas-slot. H0 = 2π a0/c, machine constants verbatim.
FAMILY/VARIANTS (every quoted H0 carries bands; no single-function H0):
  (i)   function family {p065, gm, boot} hier on the extended set;
  (ii)  flat co-read (fit_flat_pts) on the extended set;
  (iii) A-var2: NGC3972/4051 → direct Cepheid anchors (u-group 24);
  (iv)  A-var3: LT disk-slot wiring;
  (v)   A-var4: CF4-uniform — every leg-A member whose cached PGC has an
        admissible CF4 method swaps to the CF4 priority distance
        (original anchors included; unresolved/absent keep SPARC D;
        count printed);
  (vi)  A-var5: LT point cut relaxed to 0.20 (LT members only; the 9R
        ARM-R precedent; heterogeneous by design, co-read only);
  (vii) the 10S-67 subset re-read (= G10T-1 regression; the old
        baseline, always printed alongside).
σ_A: the 10S G5 bootstrap machinery verbatim (flat, 200 reps, seed 202,
galaxy resample + per-galaxy distance jitter + UMa shared draw);
rel_sig per member from its operative e_D source (CF4 for new members,
SPARC for old, UMa 2.3/18 + shared 0.9/18). Percentiles 16/50/84
printed (bootstrap-health disclosure; |median − central| < 0.5σ
expected — R45 hole-2 lesson).

## Gates (each names the letter clause its failure vetoes — trap #11)
- G10T-1 (regression, STOP): the 10T harness in 10S mode (67 galaxies,
  SPARC distances, no LT) reproduces the archived sky read
  (data/stage10s_skyread.txt): hier BE a0_A = 1.0186e-10, f_ML = 1.07,
  s_int = 0.031, u = −0.0258; flat 1.1522e-10; family {p065 1.0969e-10,
  gm 1.0334e-10, boot 1.0537e-10}; σ(H0_A) = 12.52 (seed-pinned); GD
  split hier 1.0056e-10/1.0287e-10 — point estimates to ≤0.1%, σ and
  splits to ≤0.5% (same code path, deterministic). FAIL ⇒ STOP, no sky
  read.
- G10T-2 (census, STOP): the R1–R8 tables derived at run time must
  match the frozen lists above exactly (names, methods, counts); CF4
  parse re-verified on known rows (M31 PGC 2557: ceph 24.397 ± 0.066,
  trgb 24.53 ± 0.10, sbf 24.214 ± 0.08; NGC4258 PGC 39600 DMmas ≈ 29.4
  = the maser anchor); full reclassification table printed (name, PGC,
  method, D_old → D_new, s, e_D both sources). FAIL ⇒ STOP.
- G10T-3 (port, STOP for the LT clause only): (a) port-code regression:
  the 9R ARM-V census reproduced on the 9R LT-only configuration (65
  points / 7 galaxies at the 0.10 cut); (b) overlap cross-validation:
  for SPARC∩LT galaxies present in both reductions, LT-ported points
  rescaled to the SPARC distance vs SPARC-native points at matched
  rings — median |Δ lgobs| and median |Δ lgbar| printed per galaxy and
  pooled. Bar (measured-then-pinned at gates time, disclosed): pooled
  medians ≤ 0.06 dex for lgobs and ≤ 0.12 dex for lgbar (two independent
  reductions of the same sky: velocities shared, baryon models
  independent — lgbar is the loose axis by construction). FAIL ⇒ the 4
  LT members are EXCLUDED from the primary (SPARC-reclass-only primary;
  LT quoted as co-read only); the T-letter carries the port caveat.
- G10T-4 (injection, STOP; truths beyond every edge — the R45 lesson):
  on the EXTENDED design, (a) noiseless lock-mocks at H0_true ∈ {50, 62,
  73, 85, 100} (both flanks outside any plausible quote) through the
  full hier BE fit: recovery ≤ 0.5% each; (b) noisy: 5 seeds × truths
  {62, 85}, noise = per-point ⊕ 0.08 dex intrinsic ⊕ per-galaxy vertical
  draws: |mean bias| ≤ max(1%, 2·SE); the replication SD is REPORTED as
  the leg's single-realization floor. FAIL ⇒ STOP.
- G10T-5 (power, letter-selecting): σ(H0_A) on the extended set measured
  with the CENTRAL MASKED (gates mode prints σ only). σ ≤ 12.0 ⇒ the
  T-SHARPENED path is open; σ > 12.0 ⇒ T-NULL-GROWTH (the additions
  bought no power; report why — point counts, e_D distribution).
  Forecast for the record: ~9.5–11.
- G10T-6 (provenance/lit): the Note-4 quote printed in the stage output;
  TF/FP exclusion stated; one delta literature round for any published
  a0→H0 inversion since 2026-09-20 logged in the stage docstring
  (three NOT-FOUNDs on file; SML20 = nearest neighbor, must-cite; CF4's
  74.6 named-not-imported).
- G10T-7 (composition, letter-selecting; sky time): (a) GD/non-GD split
  on the extended leg A (10S G7 machinery; the additions are GD-heavy
  BY DESIGN — both directions disclosed); (b) old-vs-new split: hier BE
  on the 10S-67 vs on the 15 additions alone, with a 200-rep bootstrap
  σ on the additions subset; splits compared against their own joint σ.
  Both inside 1σ ⇒ additions integrate. Either in 1–2σ ⇒ the combined
  headline carries the named caveat sentence. Either > 2σ ⇒
  T-COMPOSITION-SPLIT: no combined headline; the 10S-67 read stays
  operative; suspects named (GD composition = the dissident-dwarf
  interaction; port systematics; CF4-vs-SPARC zero-point seam).
- G10T-8 (UMa Cepheid co-read, annotation): the fitted shared u vs the
  frozen NGC3972/4051 Cepheid shifts (+0.072 dex and −0.035 dex vs the
  block); printed with the depth context (2.3 Mpc); annotation grade,
  never verdict-bearing.

## Letters (grammar fixed in advance; clause vetoes named above)
- T-SHARPENED: all STOP gates green, G10T-5 σ ≤ 12.0, both G10T-7 splits
  ≤ 1σ ⇒ quote the extended-set H0_A ± σ_A as the meter's updated leg-A
  reading, with the function-family band, the variant spread
  (ii)–(vi), the 10S-67 baseline alongside, and the mandatory
  disclosures. This is a LEG-A-ONLY update of the 10S M-GRAY state: the
  meter still returns no headline H0 (leg B remains power-limited); the
  honest sentence stays "H0(leg A) = X ± Y (stat) with a Z-wide function
  family band, consistent with both camps."
- T-COMPOSITION-SPLIT: any G10T-7 split > 2σ ⇒ both subset numbers
  quoted, no combined number; 10S-67 stays operative; successor named
  (port audit / CF4 zero-point seam stage).
- T-NULL-GROWTH: σ > 12.0 ⇒ census + σ reported; the 10S read stands
  unchanged; the growth path re-routes to BIG-SPARC/WALLABY (registered
  watch triggers).
- T-CENSUS-LIMITED: additions-with-points < 8 galaxies at gates time ⇒
  census-only report, no extended sky read (moot at the frozen census =
  15, kept for grammar completeness).
- Prohibitions (inherited 10S + R45): never a single-function H0; never
  "ladder-independent" ("H0-assumption-independent" only — leg A USES
  ladder rungs and now shares the ladder's pegs explicitly); never quote
  leg B at SPARC grade; the retracted 10S strings are never quoted; any
  quoted H0 carries the family band in the same sentence.

## Credence
NO cell of any map moves on ANY outcome (instrument-growth stage; same
clause as 10S). If T-SHARPENED fires, the improvement is REPORT-grade.

## Registered same commit
PREDICTIONS.md §F: the three-world discriminant table (P/N/L + lock-false
control; the ×2 peg-gearing statement with the measured 1.3–1.6 response
band) and world V (KBC-void transduction: radial a0 gradient + dipole)
— signed BEFORE any meter read sharper than 10S exists. The 10T result
is scored against §F at report grade only (σ ~ 10 separates nothing;
the table's discrimination threshold is σ(meter) ≈ 2–3).

## Order of work (firewall)
1. Instrument build; gates G10T-1..6 sky-blind (no extended-set central
   printed anywhere; asserts guard it). Output: data/stage10t_gates.txt.
2. Gate record appended below; commit prereg + stage script + frozen
   cache + census + PREDICTIONS §F (pre-sky commit).
3. Sky mode fires in one run (re-verifies gates in-process): extended
   leg A + family + variants + σ + G10T-7/8 + letter per grammar.
   Output: data/stage10t_skyread.txt.
4. Review round (fresh reviewer session, round 46); every reviewer
   number re-verified in an addendum script before adoption (standing
   memory rule); ledger row cos-10t-legregrow; NOTES/LOG/DIARY/CLAUDE.md.

## Gate record (instrument build executed 2026-09-20 night, sky-blind)
Harness: calcs/stage10t_legregrow.py gates → data/stage10t_gates.txt.
No extended-set central computed or printed anywhere.
- G10T-2 PASS: census matches the frozen lists exactly (12 on-paper
  reclassifications, UGC09992 zero-point confirmed; LT4 = CVnIdwA/
  DDO_52/DDO_133/DDO_210, all Mgas > Mstar); CF4 parse verified on M31
  (ceph 24.397±0.066, trgb 24.53±0.10, sbf 24.214±0.08) and the NGC4258
  maser anchor (DMmas 29.40±0.03).
- G10T-1 PASS: the 10S archived sky read reproduced on the 10S world —
  hier a0 d=0.001%, flat 0.003%, family ≤0.004%, σ 12.515 vs 12.52
  (0.038%), GD/nonGD ≤0.003%, u-identity exact, params (1.07/0.031/
  −0.0258) exact.
- G10T-3a PASS: 9R ARM-V census identity 65 pts / 7 galaxies exact.
- G10T-3b FAIL → THE PRE-REGISTERED CLAUSE FIRED: pooled overlap
  cross-validation (46 matched rings, 5 galaxy pairs) gives
  |med Δlgobs| = 0.008 dex (bar 0.06 — the two surveys' velocity
  fields agree) but |med Δlgbar| = 0.159 dex (bar 0.12): Oh+15's
  baryonic models sit systematically below SPARC's photometric ones
  (port V_bar²/V_tot² ~ 0.09–0.20 = DM-dominated subtraction regime;
  per-galaxy Δlgbar −0.41…−0.08). THE 4 LT DWARFS ARE EXCLUDED FROM
  THE PRIMARY (co-reads only). Primary leg A = 78 (41 anchored + 26
  UMa + 11 reclassified).
- G10T-4 PASS (on the operative 78-galaxy design): noiseless truths
  {50, 62, 73, 85, 100} recovered to ≤0.25% (both flanks beyond any
  plausible quote — the R45 lesson executed); noisy 5-seed sets at
  62/85: mean bias −1.62% within bars (2·SE 4.4/4.6%), realization
  floor SD 4.9–5.1% reported.
- G10T-5: σ(H0_A, extended, central masked) = 10.71 km/s/Mpc < 12.0 →
  the T-SHARPENED path is OPEN (10S baseline 12.52).
- G10T-6 PASS: Note-4 provenance quote printed; TF/FP exclusion stated;
  three NOT-FOUNDs logged (2026-09-20).
- G10T-7/8 armed for sky time.
Membership note for the record: the prereg's "expected 82" resolves to
PRIMARY 78 + 4 LT co-reads via the G10T-3b clause — the exclusion is
the gate working as designed, not an amendment (no rule changed).

## ROUND 46 ADOPTION RECORD (2026-09-20 night; report in
## REVIEW-ROUND46-OPUS.md, uncommitted; verification =
## calcs/round46_addendum.py GA (blind, pre-committed 508d2de) + GB;
## operative corrected output = data/stage10t_skyread.txt re-fired
## from the patched harness, with the AS-FIRED original preserved in
## data/stage10t_skyread_asfired.txt and at commit 87e689e; operative
## verdict = data/stage10t_verdict.txt)
Ruling adopted IN FULL: UNPATCHED HOLE YES, letter-scoped, not
verdict-reversing; T-SHARPENED → T-SHARPENED (ENGINE-QUALIFIED); all
14 conditions executed (C1–C13 in the patched harness + addendum;
C14 = no credence move, as pre-registered).
The two severe holes, in the record's own words:
- F1: the as-fired "± 10.7 (stat)" was the FLAT engine's bootstrap
  dispersion attached to the HIER central; the pre-registered
  bootstrap-health check (|median − central| < 0.5σ), which this
  prereg wrote down and the stage printed inputs for but never
  EVALUATED, fails against the primary central. The corrected quote
  carries the primary-engine bootstrap σ, with the flat dispersion
  retained only as an envelope spanning the treatment split. The leg
  is thereby SYSTEMATICS-LIMITED (ν-form band > stat), and the
  successor priority flips from more-anchors to PIN-THE-ν-FORM.
- F2: the G10T-3b cross-validation pool was silently restricted to
  SPARC-kept overlap galaxies (dropping the two Q=3 pairs, NGC_2366
  with 44 matched rings and DDO_50); under the literal registered
  pool the ring-median statistic PASSES. The LT exclusion STANDS on
  the robust statistics (all seven galaxy medians negative, ring mean
  and galaxy median both beyond the 0.12 bar), now the pinned decision
  rule in the harness. Materiality either way ≤ 0.2 km/s/Mpc.
Also adopted: F3 (the loose alternation is a coordinate-descent fixed
point ~1 km/s/Mpc above the profile optimum; the operative estimator
is the deep-converged fit, the loose fit is kept only as the 10S
regression tie; never quote 0.01-precision deltas across estimators);
F4 (variant band repaired: primary-eligible variants run LT-free,
A-var4 swaps all four CF4-bearing UMa members, the flat treatment
variant is named at its own value in the band sentence); F5 (the
composition disclosure was INVERTED — the additions are GD-LIGHT,
1/11; corrected to the measured fractions, with the three zero-point
anchored dwarfs named); F6 (old-vs-new relabeled a diagnostic — the
gate was powerless, trap #16); F7 (injection arm: per-truth
independent rng streams, 8 seeds, interior truth added; the
self-widening bar named; the realization floor carried into the
letter); F8/C9 (every "H0-assumption-independent" now carries the
peg-sharing clause; gearing quoted as the measured 1.3–1.6 with ×2 as
the deep-limit asymptote); F9/C10 (the absolute-scale provenance
sentence is from the CF4 ReadMe ABSTRACT, not Note 4 — CORRECTION to
this prereg's "Distance-source pins" section, which mis-attributed
it); F10/C11 (CORRECTION to the pin's unresolved-name sentence: NOT
all unresolved names are flow LSBs ≳45 Mpc — D564-8 and D631-7 are
f_D=2 at ~8–9 Mpc, KK98-251 is flow at 6.8 Mpc and D512-2 at 15.2;
KK98-251's PGC must be resolved before the next refire); F11 (var4
count 29-with-points; the pin-7 UMa 0.9-double-count seam disclosed —
inherited from 10S, immaterial; gates_all wired into the letter;
"computed anywhere" → "REPORTED").
R46 AMENDMENT A1 (2026-09-20, logged after the corrected harness's
run 1 STOPPED at the strengthened G10T-4 and BEFORE any corrected
verdict; the stopped run is preserved verbatim in
data/stage10t_skyread_r46run1.txt): executing C7 (8 independent seeds
per truth, interior truth added) exposed a REAL one-sided estimator
bias at the EDGE truth only — truth 85 recovers −2.96 ± 0.73% (all 8
seeds negative) while truths 62 (+0.87 ± 1.01%) and 73 (−0.66 ±
1.31%) are unbiased and the noiseless arm is exact everywhere: a
noise-induced calibration property of the deep estimator ~20 km/s/Mpc
above the actual reading, invisible to the as-fired 5-shared-seed arm
(which passed its pre-registered bar). Amended semantics,
direction-neutral for the reading: quote-region truths {62, 73} (plus
the addendum's 66) are LETTER-VETOING; the edge truth 85 is a
CALIBRATION-MAP leg whose failure adds a VALIDITY-DOMAIN clause to
the letter (the quote is valid where it sits; any future reading
above ~75 requires recalibration first) instead of stopping a reading
deep inside the validated region. The finding itself is banked as
instrument knowledge for the BIG-SPARC-era refire.
Standing additions to the program's rules (from this round):
- A pre-registered health check must be EVALUATED by the harness with
  a printed PASS/FAIL, never merely have its inputs printed.
- A dispersion is quoted only for the estimator that produced the
  central ("engine-matched σ"); a different engine's dispersion may be
  co-quoted as a labeled envelope only.
- A cross-validation pool's implicit joins (name→index maps, kept-set
  filters) are part of the gate definition and must be enumerated in
  the prereg.
