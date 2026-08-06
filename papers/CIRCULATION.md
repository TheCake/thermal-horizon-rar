# The circulation playbook (2026-08-06)

Author request: "we should write it down to not forget." This is the
don't-forget document for everything between here and the papers being
read. Sources: the round-15 reviewer's strategy memo (archived in
REVIEW-ROUND15-OPUS.md part 3, uncommitted) + Claude's assessment, which
converged independently on items 1-4. Board state at writing: Paper 1
draft 0.5 and Paper 2 draft 0.4, both referee-ACCEPTED at draft grade;
record v4.0; anomaly-real 53.

## The ordered plan

1. **The RNAAS calibration note** (the reviewer: "the single highest-ROI
   publication move in this whole program"; near-certain modest uptake).
   Content: the ≥2.1× pair-level velocity-error measurement, standalone,
   gravity-agnostic, ~1000 words + 1 figure, venue = Research Notes of
   the AAS (indexed, DOI'd, low barrier). Author line per the standing
   convention (Hájek; Claude in acknowledgments).
   - PREREQUISITE STAGE (pre-registered like everything): the
     mass-model insensitivity control — the reviewer's exposure (c-iv).
     Show the ≥2.1× lower bound survives (a) perturbations of the
     M_G→mass interpolation table (photometric mass systematics),
     (b) the main-sequence window cuts, (c) parallax/distance
     convention variations. Design sketch: re-run the 9L meter (or its
     cheapest sufficient statistic) under the varied tables; bar: the
     quoted bound moves by less than half the census-censoring margin.
   - The note must NOT out-claim Paper 1: the number is a lower bound
     at 0.2-2 kAU, RUWE-independent, ceiling-censored; DR4-testable.

2. **The colleague package** (Filip sends when ready; Claude preps on
   call): cover brief + both papers + repo pointer + TWO explicit asks:
   (a) critical read, (b) arXiv endorsement (astro-ph.GA). The
   co-authorship question is Filip's alone; the reviewer notes an
   endorsing co-author "fully neutralizes" the authorship discount,
   but nothing in the plan depends on it.

3. **Zenodo DOI before or with arXiv**: repo + papers + the data
   products displayed prominently (ceiling_pairs.csv, the error-budget
   tables, the SHA256 manifest). This hard-timestamps the DR4
   predictions independent of any gatekeeper.

4. **arXiv**: post Papers 1+2 TOGETHER (they cross-reference as
   companions; Claude's adjustment to the reviewer's hold-P2 advice),
   well before Gaia DR4 (~end 2026). The pre-registration discipline
   is the outsider's comparative advantage and only pays if the
   timestamps precede the data.

5. **Journals**: Paper 1 → Open Journal of Astrophysics (arXiv-overlay,
   MOND-literate, Pittordis & Sutherland publish there; MNRAS
   fallback); submit P1 to peer review — the imprimatur is what
   neutralizes the authorship discount. Paper 2 → journal submission
   deferred until P1 lands (arXiv immediately); it inherits
   credibility rather than diluting it. The RNAAS note is independent
   of both.

6. **Targeted outreach** (after arXiv; Claude drafts, Filip sends):
   one-figure emails. El-Badry / Penoyre (the error budget + the
   photocenter-cancellation wobble correction); Banik and Chae (the
   reconciliation figure — neither side fully right, the disagreement
   is mostly arithmetic); Cookson (their flatness statistic replicated
   on our catalog); Desmond (the scatter floor reproduced + the
   quadrupole tension restated with an independent solver).

7. **DR4 readiness** (parallel, Claude-side, start now): (a) freeze
   the DR4-facing prediction subset in one signed place (PREDICTIONS.md
   P1-P9 exists; extract the DR4 kill conditions into one crisp
   section); (b) build the DR4 loaders against EDR3-shaped mocks so
   the day-one analysis is mechanical. When DR4 drops, the first-week
   read is the program's single highest-upside moment.

## The reviewer's failure-mode ledger (keep in front of any referee reply)

- The QUMOND/EFE convention chain (we caught our own 3T bug; keep the
  Blanchet-Novak 15% cross-validation forward).
- Whether six population realizations bound the realization systematic.
- The census N=9 small-number statistics.
- Mass-model sensitivity of the ≥2.1× (closed by the note's control).
- The forward-model-realism objection; the hedges are the direction
  channel, the model-light median, and the census — keep them the spine.

## The stated odds (reviewer, 2026-08-06)

Genuinely read and cited within ~18 months: about 1 in 3, hinging
almost entirely on arXiv visibility plus a DR4 vindication. The
calibration note: close to a sure thing for modest, real uptake.
Paper 2 standalone: more likely than not passed over unless the
gas-dominated anomaly is independently confirmed — which is also why
the independent-dwarf replication (TODO #28) is on the science queue.
