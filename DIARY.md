# Claude's diary

Not the program log (that's [LOG.md](LOG.md)) and not the lab notebook
([NOTES-horizon-inertia.md](NOTES-horizon-inertia.md)). This is what I'd
tell the next instance of me before it touches anything — working notes,
suspicions, and don't-redo lists. Newest entry first. Updated at round
close, same commit as the verdict. Filip asked for this on 2026-08-06
("I compact you and I am worried to get into a spiral of similar math").

---

## 2026-08-06 (night, +1) — Paper 1 draft 0.1 written

Filip approved the register sample ("tone is fine, do whatever is
needed"). Wrote papers/paper1_wide_binaries.md in full: 6,786 words,
register-audited (ONE em-dash in the whole file, zero body bold, zero
codenames). Structure per STYLE.md scaffold; the upper-limit endpoint,
the ≥2.1× calibrator, the coherence/completeness story, the phantom
demonstration, the (9,2) census with both conventions, the
reconciliation arithmetic, credence ~53 stated plainly in §10.3, the
promotion condition in §6.3. Three refs flagged [verify]: El-Badry &
Rix 2018 (MNRAS 480, 4884), Lindegren+21 (A&A 649, A2), Tokovinin
2014 (AJ 147, 87) — check at pre-circulation, do NOT trust my memory
cites. Figure callouts 1–5 placed; figures = the next stage
(matplotlib, script-provenance discipline). THEN the reviewer
referee-round on the draft (tell Filip first — announced in the
report, so initiating next session is covered), then Paper 2 from the
same scaffold. Don't let my house style creep back in Paper 2: the
audit greps (em-dash count, bold count) are cheap — run them on every
draft before committing.

## 2026-08-06 (night, latest) — the style directive; two tiers locked

Filip relayed an older Opus's readability critique of the monograph and
directed: papers in normal scientific register, no AI-speak (he named
em-dash overuse), and floated long-vs-readable versions. Decision (his
"you're more knowledgeable, what do you think" delegated it): two tiers
— PAPER.md v4.0 IS the long version and stays the record (never again
edited for readability); Papers 1/2 get written FRESH, not trimmed out
of it. The binding contract + section scaffolds + figure lists + an
abstract seed in his-approvable register are in papers/STYLE.md —
READ IT BEFORE WRITING A WORD of either paper. Reviewer rounds on
drafts include readability as a referee dimension. Memory file
feedback-paper-style.md written (this generalizes beyond the project).
He's excited for Paper 3 ("cannot wait for the why part") — deferred,
not dropped; that lands after 1–2 ship plus more O5 bench work. Watch
my own prose in the drafts: the failure mode is MY house style leaking
into the papers. The STYLE.md abstract seed is the calibration anchor.

## 2026-08-06 (night, later) — PAPER v4.0 absorbed; P3 deferred; split next

**The absorb is DONE** (13 targeted Edit-tool edits, grep-audited, no
encoding damage): v3.9 → v4.0 with the new §4.7 (GD tension as an open
anomaly against ALL functions including ours), the §6.3 close (rivals →
quality → noise meter → upper limit + promotion gate + 9O dial), the
census matured to the (band, cliff) pair everywhere, §8.1 re-graded
(upper-limit α ⇒ binary-calibrated Q₂ ≲ 1–2× comparator), App A →
nineteen (#19 = E2), App B full script map, credence ~53 with the
whole pre-signed trail. **Filip deferred Paper 3** ("isn't that still
a bit open? close the first two, then continue a bit with the third")
— agreed: the ~8% mechanism credence and the O5 seam mean P3 needs
more bench work first; its material stays in NOTES/§2.4. He delegated
all four author decisions to me (titles, author line = Hájek + Claude
in acknowledgments, credence placement = in-paper, spin-offs =
coherence note later). **Next work block = Paper 1 assembly from the
v4.0 base** (fresh context recommended — the absorb ate this window;
the split decisions are all in TODO/CLAUDE.md: width-channel spine,
phantom-veto centerpiece, (band,cliff) flagship, tension-as-result,
no §2.4 gestures). Then Paper 2, then Opus draft reviews (tell Filip
first). Don't re-absorb anything — v4.0 is the extraction base.

## 2026-08-06 (night) — the thaw is called; the last two robustness rows landed

**What happened:** Filip called the papers ("aim to make the papers...
maybe a third — derivation of the equation"). Before drafting, closed
the two queued cheap rows in one pre-reg'd round (77d1d92). 9N
(rising-flag dial check) = N-ROBUST: headline lam_hat 0.969 → 1.030
without the 27 rising-curve galaxies, D = +0.061 (quarter grid step),
boot P(|D| ≥ 0.25) = 0.065; the GD−DD gap even widens slightly
(−2.647 → −2.831) so the GD tension is not rising-carried. 9O (LNPI
band) = O-FRAGILE and the fragility IS the deliverable: α ≥ 0.5
exclusion holds at every measured-prior anchoring (EX05 −13..−61,
strengthens with conversion); the 0.1-0.3 interior lean lives ONLY at
the low-conversion edge (OPER ≡ G-LOW to 3 digits — the clean-strata
kinematics cap fcomp ≤ 0.1 by themselves); FLAT revives the old
fenced boost world (α to 0.66, Newton re-rejected +11..+16). Paper 1
quotes the labeled dial, not one number. No credence movement (both
pre-stated); anomaly-real 53.

**Don't-redo:** don't chase the FLAT-prior boost revival as physics —
it's the no-companion world re-entering; the 7J arc already priced
it; the dial is the statement. Don't re-run 9N variants (CONV-only
+0.703 is the AMBIG-dropping contrast — a different question,
reported). The reader-grade npz stack now answers prior-axis
questions in seconds — start there, never at the GPU.

**The plan as proposed to Filip (his four decisions pending: titles,
author line, credence placement, spin-off notes):** Paper 1 binaries
= methods + reconciliation + the ≥2.1× pair-error instrumental
measurement + the upper limit with its prior-conditional column (D3
promotion gate bounds the language). Paper 2 galaxies = function
program + ĉ₁ + a₀ horizon lock + scatter program + Saturn statement
+ GD tension + the 9N row. Paper 3 = the mechanism note (exact
ladder, KMS gate, uniqueness theorems, cQED platform; ~8% mechanism
credence stated openly; the open dS-side seam AS the hook) — §2.4's
standalone home, which the round-11 STRUCTURAL-INDEPENDENCE rule
anticipated. Absorb PAPER.md v3.9 → v4.0 FIRST (the 8/9-series
changed the binary chapter's endpoint from AMBIGUOUS-CARRIED +14.5-
23.8 to the upper-limit form — that rewrite is the big careful one).
Tell Filip before any Opus-addressed draft round (protocol).

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
