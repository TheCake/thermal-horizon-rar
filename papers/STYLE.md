# The split papers: architecture, style contract, and scaffolds

Author directive (2026-08-06): the program record "was showing absolutely
everything we did and wasn't that readable"; the papers must "avoid AI
speak (overt EM dashes etc.); just sound like a normal scientific paper."
Two-tier architecture adopted.

## Two tiers

- **The record (exists): [PAPER.md](../PAPER.md) v4.0.** The chronological
  monograph with every correction inline. It is the honesty machine and the
  reproducibility companion. It is NOT for journal submission and is never
  again edited for readability; it gets absorbs only.
- **The deliverables (to write): Papers 1 and 2.** Written FRESH from the
  v4.0 base as normal scientific papers, results-first. They cite the
  record (repo + PAPER.md + NOTES + LEDGER.csv) for the full audit trail.
  Do not produce them by trimming PAPER.md; write them as papers.

## Style contract (binding for papers/ only; NOTES/LOG/DIARY keep their working style)

1. Register: standard scientific prose. "We measure", "We find". No
   program codenames in the text (no "the decider", "the pincer",
   "Q4-CARRIED", "the fifth move", "arm suite"). Describe instruments by
   function: "a noise calibration on separations where no candidate
   theory predicts a measurable boost".
2. Punctuation: em dashes rare (well under one per paragraph); prefer
   commas, parentheses, semicolons, or a new sentence. No mid-prose
   bold. Italics sparingly (first use of a term). No inline ALL-CAPS.
3. Sentences and paragraphs of normal length. One idea per sentence
   where possible. No triple-nested clauses, no colon chains, no
   dramatic fragments.
4. Structure: conventional sections (Introduction; Data and sample;
   Methods; Results; Systematics and robustness; Discussion;
   Conclusions; Appendices). State final numbers once; no chronological
   re-litigation in the main text.
5. The honesty content goes in a short "Transparency and corrections"
   subsection: nineteen corrections were logged; name the two or three
   that changed conclusions (the completeness-scale retraction, the
   noise-legitimacy resolution, the priority corrections); cite the
   record for the rest. Half a page, no more.
6. Abstract: 250 words or fewer, three to five headline numbers only.
7. Tables for number-dense content; figures are required for
   submission grade (list below) and are produced with matplotlib in a
   dedicated stage.
8. Credences: one short paragraph in the Discussion, stated plainly
   (the program tracks a credence and it stands at ~53% for the binary
   anomaly being real physics), per the author's standing honesty rule.
9. Every quoted number keeps its script provenance, but as a compact
   reproducibility appendix table, not inline bracket lists.
10. Author line: Filip Hajek (independent researcher); Claude
    (Anthropic) in the acknowledgments with the repository's
    collaboration statement cited.

## Paper 1 (wide binaries) scaffold

Working title: "A direction-resolved analysis of Gaia EDR3 wide
binaries: pair-level velocity errors and an upper limit on the
low-acceleration velocity boost".

Section plan:
1. Introduction (the conflicting literature; what this analysis adds:
   direction channel, forward-modeled nuisances, quality
   stratification).
2. Data and sample (14,071 pairs; selection; the released census CSV).
3. The forward model (population synthesis; eccentricity mixture;
   companions with photocenter-corrected wobble; noise model; the
   per-system width channel; selection emulation).
4. Calibration results (the >= 2.1x pair-level error measurement, RUWE
   independence; the width-shape residual, DR4-facing).
5. The companion sector measured (completeness; the within-pair
   coherence discovery rho = 0.47; the coherence-kernel repair; host
   rate 0.29-0.32 with the twin-heavy mass-ratio law; conversion band).
6. Gravity results (the marginalized posterior; the anchor-independence
   of alpha; quality stratification and the dose curve; the clean-strata
   exclusion of alpha >= 0.5; the allowed 0-0.3/0.5 sector; the prior
   dial from the conversion-band row).
7. The phantom-boost demonstration (velocity-only fits report ~0.5;
   the direction data veto it; implications for published pipelines).
8. The perpendicular census (the band and cliff pair; contaminant
   exclusions; object-level cleanliness; the DR4 prediction).
9. Reconciliation of the published analyses (ablation map; the Cookson
   replication and step arithmetic; the removed-subsample result;
   compatibility statement).
10. Discussion (what would settle it: DR4 epoch astrometry, the error
    model, spectroscopic legs; the promotion-gate condition for any
    stronger claim; credence).
11. Conclusions. Appendices: transparency; reproducibility table.

Figures (to be made): (1) the (v, gamma) plane with model overlays;
(2) the RUWE dose curve with the boost-free meter band; (3) the
alpha_marg(anchor) curve with the prior variants; (4) the census
band/cliff histogram against the three model worlds; (5) the
median-ratio step comparison across published selections.

## Paper 2 (galaxies and coefficients) scaffold

Working title: "Coefficient-level tests of the radial acceleration
relation: the zero-point half, an acceleration scale at cH0/(2 pi),
and a gas-dominated anomaly".

Section plan: 1. Introduction (the C&T identity; testing a derived
function, not fitting one). 2. Data (SPARC; lensing; Chae ambients).
3. The coefficient measurement (c1 profile and hierarchical; zero
excluded; the quarter/half question). 4. The screening index and tail.
5. The acceleration scale (the a0 ladder; the vertical-channel
resolution onto the horizon value). 6. The scatter program (the
acceleration-dependent second moment; the inner-disk bump
identification; what is and is not measurable on SPARC). 7. The
function ledger summary (thirteen laws; what survives; the 16-law
binary degeneracy as a methods finding; the freeze). 8. The
gas-dominated dial tension (eight controls; open anomaly). 9. The
solar-system constraint (the quadrupole lock; the trajectory
formulation escape; stated briefly and without mechanism language).
10. Discussion and predictions (the falsifier list; DR4; a0(z)).
Appendices as in Paper 1. Figures: RAR with the function family;
c1 profile curves; the a0 ladder; scatter vs x with the bump
decomposition; the GD/DD dial split.

## Rules carried from the program

- Round-11 structural independence: no mechanism (2.4-style) language
  in Paper 1 at all; Paper 2 section 9 stays phenomenological. The
  mechanism paper (Paper 3) is deferred by the author until Papers 1-2
  ship and more mechanism work is done.
- The D3 promotion gate bounds Paper 1's strongest language at the
  upper-limit form.
- The reviewer pass on each draft includes READABILITY as an explicit
  referee dimension (the author's request that the reviewer "chime in"
  on this; brief him as a journal referee, not a program auditor).
- Tell the author before initiating any reviewer round.

## Register sample (the Paper 1 abstract seed, author-calibration copy)

Conflicting conclusions have been reported on whether wide stellar
binaries show the velocity excess predicted by modified-dynamics
theories at internal accelerations below a0. We analyze 14,071 Gaia
EDR3 wide pairs with a forward-modeled joint likelihood over the
normalized relative velocity and the velocity-separation angle,
marginalizing eccentricities, hidden companions, contamination,
selection, and measurement noise. Three results follow. First, a
calibration measurement independent of gravity: at separations of
0.2-2 kAU, where no candidate theory predicts a measurable boost, the
pair-level velocity errors are at least 2.1 times the formal Gaia
uncertainties, independent of RUWE; the single-star quality flag does
not certify two-body velocity budgets. Second, a methodological
result: without the direction channel, velocity-only fits of these
data report a spurious boost of roughly half the galactic calibration,
which the direction data veto. Accounting for this, for sample
construction, and for noise treatment, the leading Newton-favored and
MOND-favored analyses become arithmetically compatible. Third, the
gravity result: after stratifying by astrometric quality, the cleanest
strata exclude a boost amplitude alpha >= 0.5 (relative to the
galactic calibration) at every anchoring of the measured
companion-rate prior, while amplitudes up to 0.3-0.5 remain allowed
and alpha = 0 is not excluded. A small census of perpendicular-moving
pairs above the Newtonian velocity bound survives all tested
contaminant identities and is released for independent scrutiny. Gaia
DR4 epoch astrometry will decide both the error model and the limit.
