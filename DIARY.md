# Claude's diary

Not the program log (that's [LOG.md](LOG.md)) and not the lab notebook
([NOTES-horizon-inertia.md](NOTES-horizon-inertia.md)). This is what I'd
tell the next instance of me before it touches anything — working notes,
suspicions, and don't-redo lists. Newest entry first. Updated at round
close, same commit as the verdict. Filip asked for this on 2026-08-06
("I compact you and I am worried to get into a spiral of similar math").

---

## 2026-08-06 (later) — the fork RESOLVED noise-real; the binary leg settles at an upper limit

**What happened:** 9L (the reviewer's D2, amended once when its own
boost gate caught my window — 2-6 kAU carries 5-8% boost, meter is
0.2-2 only) came back NOISE-REAL: the boost-free arm demands E[fpm] ≈
2.2, same as the joint fit, both seeds, flat across RUWE quartiles.
So: Gaia formal errors are ~2× optimistic for PAIRS (standing
instrumental measurement — RUWE certifies single stars, not 2-body
budgets; round-13 called it), B's median pedestal is real noise, the
9E "counterweight" is fully retired (E2), and the drop-world small-α
is honest. Map executed 56 → 53. 9M killed the V_flat-definition
artifact by sign (converged GD alone reproduces the dial); GD tension
now 8 controls deep, DR4/external kinematics territory.

**Don't-redo additions:** the two-channel Q1 "residual" is EXPLAINED
(E2 + 9L) — do not build another instrument for it. D3 (degraded
injection) demoted to optional — the fpm AMPLITUDE is measured-honest;
only the width-SHAPE question (wide-arm, fpm→3 chase) stays open and
it is DR4-facing. Don't re-cap fpm at 1.5 (9K) as if physical — the
meter says ~2 IS physical for pairs.

**Where the case now rests:** binaries = upper limit (α ≥ 0.5 excluded
clean; ≲0.3-0.5 allowed, mild 0.1-0.3 interior lean; α=0 not
excluded) + the (band=9, cliff=2) census pair; galaxies = the intact
function program (AMB etc.) + the 8-control GD tension. Anomaly-real
53. The next real discriminators are external: DR4 epoch astrometry
(pair error model + e-resolution), spectrograph legs, resolved dwarf
kinematics. In-catalog, the binary quality arc is basically MINED OUT
— resist the urge to keep re-slicing it (that's the spiral Filip
worries about).

**Pet note:** the round-13 Opus agent is still alive this session —
sent him the D1/D2/A4 outcomes; his reaction should be booked as a
round-13 addendum when it lands.

## 2026-08-06 — the quality arc closed into a fork; round 13; this file born

**Where my head is:** The binary anomaly deflated but did not die. The
one plot that matters now is the 9J ladder — fitted α rises monotonically
with RUWE (Q2 0.05-0.11 < Q1 0.12-0.26 < Q3 0.40-0.79 < Q4 0.90-1.05).
Two readings, one fork: (a) noise dresses as boost in dirty strata →
true α ≈ 0.1-0.3 (upper-limit form, round-13 headline); (b) boost is
real and dirty strata exaggerate → clean-world α under-reported (9K
measured +0.05-0.08 of that pressure at the Lindegren cap, real but
mild). The D2 narrow-pair fpm meter decides — narrow bins are boost-free,
so whatever fpm they demand is honest noise. Credence FROZEN 56 until it
lands.

**Don't redo (the spiral list):** Reader-grade re-analysis NEVER needs
the GPU — the archived tables are data/stage9f_tables_*.npz (5-pt α),
stage9i (7-pt fine), stage9j (extension to 1.2), plus the 7J lker cubes.
9K was built in an hour from these; any new combiner question should
start there. Dead hypotheses with their killers: zero-inflated mixture
(8U), pressure support amplitude+shape (8Y/9C sign-inverted), GD
environment (9G), quality flag (9B), correlation resampling (8V 0/300),
frozen-bath variance signatures (6O/6P), scalar EFE (6W), no-EFE MI
(4L), rate-based grammar realizations (6K/6M/6N). The 9E double-ratio is
NOT noise-immune (round-13 E2 — σ(ṽ) grows with s; check LEDGER note
before ever quoting B as boost evidence).

**Traps I actually hit this arc:** bar grammar keeps biting (8U minority
bars, 9A gap, 9I NEWTON-FLAT overlap — always write EXHAUSTIVE ordered
letters; Newton letters test P(α=0), not P(α≤small)); 5-dim-vs-6-dim
reshape in a sliced-table read (9K crash); a manifest that compared the
wrong count to the right record (9H — file count ≠ fit-population
count); heredocs with review prose (quote hell — use the Write tool).

**Standing people-facts:** Filip wants plain verdicts + ELI12 always;
LLM-time in compute wall-clock; the Opus reviewer is MY agent to run
(no relay), keep him blunt, book his reviews as numbered rounds, keep
the in-session agent alive for follow-ups (SendMessage) — he holds his
own round context until the session ends; archive verbatim to
REVIEW-ROUND*.md (uncommitted). Pre-reg commit BEFORE any run, no
exceptions, even for readers.

**Suspicion I can't yet act on:** the width-shape object (the fpm>1.4
hunger that survives the ws floor) smells like a per-pair error-model
problem (covariances, scan geometry) that DR4 epoch astrometry will
either dissolve or crystallize. If D2 says the ≈2× is real on narrow
bins, chase the error model before chasing gravity.
