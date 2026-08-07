# Claude's diary

Not the program log (that's [LOG.md](LOG.md)) and not the lab notebook
([NOTES-horizon-inertia.md](NOTES-horizon-inertia.md)). This is what I'd
tell the next instance of me before it touches anything — working notes,
suspicions, and don't-redo lists. Newest entry first. Updated at round
close, same commit as the verdict. Filip asked for this on 2026-08-06
("I compact you and I am worried to get into a spiral of similar math").

---

## 2026-08-07 (close, +8) — 9U the gamma meter: built it, gated it, and it says the sky can't answer yet

Filip: "continue, we gotta get to the bottom of this." The bottom of
O5-AVERAGING is one number — where p sits continuously and with what σ.
So I built the instrument 9S/9T queued by name (pre-reg 21376aa BEFORE
run): fine ν_p profile on both hier treatments + injection power gate +
40-rep paired galaxy bootstrap. It came back U-POWER in 2.3 minutes of
compute (my 4–6 h estimate was ~100× conservative — warm-chained
Nelder-Mead hier fits converge in seconds; remember that for planning).

What the next me must not redo or misread:
- p̂ = 0.6471 ± 0.0746 vertical-primary (plain 0.6927; split 0.046 =
  sub-σ, no flag). σ_p is a CLIPPED LOWER BOUND — 31/40 reps edged the
  ±0.06 grid, and after the pre-authorized A-widen they piled at the
  ±0.09 edges too. Percentiles 16/50/84 = 0.56/0.615/0.74. The per-rep
  tail wells are ~7–10 lnL over ±0.09; galaxy resampling tilts more
  than that. The tail LOCATION is realization-dominated — the 3A lesson
  wearing a galaxy coat. DO NOT quote any sharp p location from SPARC
  hier, ever, without this row.
- All three bands inside 1σ (one-swing +0.81, floor −0.32, full-avg
  −0.55). NOBODY won. The 9T kill did NOT fire (needed σ_p ≤ 0.02, got
  ≥ 0.075). Credence: pre-signed U-POWER cell = HOLD 12, mechanically.
  anomaly-real untouched. No Opus round — nothing verdict-bearing.
- The honest-updating pair that DID land: (1) 9T's booked tension
  (r = 0.398 < floor 0.454) DISSOLVES at measured error — it was the 5G
  grid point pretending to be a measurement; row annotated. (2) 9S's
  "r ≥ 0.315 bound" demotes to point-conditional — r > 0 holds at ~2σ
  lean; mech-9s-detuning CO-QUOTED → gal-9u-ptail. If a future round
  cites "the 9S bound" as sharp, correct it.
- Trap fingerprint (new, general): a bound inherited from a GRID POINT
  with no σ is not a bound — 9S inverted 5G's p̂ = 0.65 as if exact;
  the first real error bar moved its floor by half its value. Same
  family as round-10's "flat in center, collapse on width".
- The warm-lite caveat is disclosed in the row: per-rep location noise
  may inflate the scatter, but even 2× overestimate leaves U-POWER.
  Don't spend a round re-running the bootstrap at full tolerance
  hoping for 0.02 — the full-fit curvature alone is 0.031.
- THE LIVE SUCCESSOR is the VOID ASYMPTOTE (P1 annotation): p_void =
  ½ + r/2 is gate-free — the signal isn't divided by g ≈ 0.75, so the
  same absolute σ_p buys ~2× the r-resolution, and the 9T kill-band
  [0.727, 0.750] sits AT the ceiling where detection is cleanest.
  Next O5-AVERAGING move on Filip's go: scope what void-RAR data
  exist (Desmond? DR4-era?) — a literature scout job first, not a fit.
- Board: ledger 190 rows (gal-9u-ptail), worldtable 232 tokens six
  gates PASS, everything through the 9U verdict commit.

Where the fork stands after today: the mechanism's §2 skeleton for
Paper 3 is unchanged (L1 + L2 + 9T + 9S-as-calibrated + void band),
but the averaging pillar now has a MEASURED data-grade wall on the
SPARC side. The γ question is real, open, and correctly parked with
its instrument spec — exactly the state Paper 3 can print honestly.

Filip's directive: "Do whatever you think should be next. Can be even
concurrently if possible." I ran the detuning bound (9S, cheap reader)
and the GD replication attempt (9R, the TODO-28 falsifier) in one
block. Both pre-registered before execution (4b46cd5 / 7f32d19).

What the next me should know:
- 9S landed S-BOUNDED, gates 4/4 first run: the measured galaxy tail
  (5G hier p-hat = 0.65) forces exchange weight r >= 0.315 in EVERY
  gate treatment => |delta|/g_c <= 0.77. The round-16 A9 objection now
  carries a measured answer: far-off-resonance exchange is excluded by
  data we already had. HONESTY SHAPE: the bound is hier-anchored (the
  4H flat treatment is too soft at its -1 sigma edge — printed, not
  hidden); binaries are gate-blind. Successors: r-ladder consistency
  re-fit (label it consistency — the 7I freeze forbids new function
  SEARCH), and the void r-meter (P1 annotation: p_void = 1/2 + r/2,
  gate-free).
- 9R refused to answer, correctly: R-POWER-LIMITED. The construction
  is GOOD and reusable (V_bar^2 = V_tot^2 - V_DM^2 matched-radius
  subtraction from Oh+15's rotdmbar/rotdm; 95% ring match; port
  regression hit -1.542 exactly; wiring license passed) — but both
  arms failed their own injection power gates (GD world -1.31 recovers
  at -0.5; the relaxed arm rides the -2.00 edge). THE SKY WAS NEVER
  FIT — the script's control flow reads the sky only after power
  passes. Keep that design. Next instrument: 9R-b = two-world CONTEST
  (LR of +0.97 vs -1.31, own injection calibration) — the ~1.1 world
  separations in ARM-V say a coarser statistic might be powered where
  the profile is not.
- Trap fingerprint (the 6th scout overstatement): Scout B claimed
  Oh+15 has "full baryonic decomposition" — the VizieR byte tables
  say scaled TOTAL + DM-only curves. Primary-source byte descriptions
  decided it; the subtraction construction rescued the stage. Never
  design a stage off a scout's data claim without reading the ReadMe.
- Day totals: anomaly-real 53 untouched all day; mech conditional ~8
  (held by the 9Q map); ledger 187 rows; worldtable 226 tokens.
- LATE ADDENDUM (9R-b, pre-reg 025b135): the contest form ALSO landed
  B-POWER-LIMITED — by ONE injection (ARM-R headline side 15/20 vs
  the locked 16/20; GD side 20/20; 5/40 total). The bars held; no
  post-hoc softening; the sky was never read by ANY of the three
  instruments. TODO 28 is now CLOSED at the Oh+15 public-data grade
  per the pre-registered clause. Next me: do NOT re-run these
  instruments on the same input hoping for different luck — the
  reopen triggers are Iorio+2017 curves (likely tips the contest) or
  DR4-era data. Ledger 188 rows / 228 tokens after this booking.
- NIGHT ADDENDUM (9T + ROUND 17, the resonance theorem): the
  round-16 A9 detuning flank is CLOSED and survived its own
  red-team SCOPED (r >= 0.454; delta small vs EVERY rate — keep his
  separability framing); mech conditional 8 -> 12 by the pre-signed
  map. THREE things the next me must never do: (a) quote 9T as
  "resonance established" — the 1/2 is a time-average the universe
  may not have had time to take; (b) forget that the measured
  p = 0.65 implies r = 0.398 < the 0.454 floor (tension-at-
  resolution, honestly booked; the kill is sigma_p <= 0.02 with
  p < 0.67); (c) apply the tail claim to binaries (their boost
  witnesses only gamma ~ 0.29 — excluded by scope). O5-AVERAGING is
  now THE question: gamma ~ 1 predicts p ~ 0.587, full averaging
  predicts ~0.69, the sky says 0.65 (grid-limited) — designing the
  gamma-measurement (a sigma_p <= 0.02 hier tail instrument, or the
  finite-gamma p(gamma) curve as a new fit dial) is the single
  highest-value P3 move. Keep his gems: Landau-Zener adiabaticity
  39.5 REINFORCES equal sharing; the phase-budget framing is partly
  decorative (real content: delta <= H/2pi AND g_c >~ H). Ledger
  189 / 230 tokens after this booking.

## 2026-08-07 (+5) — 9Q the lemma round + ROUND 16: the red-team earned its keep

Filip opened the theory gate ("go for it even via your intuition —
though of course measure everything") and explicitly authorized Opus
idea-bouncing agents. Shape of the round: pre-reg 471bef4 with
bars/map/disclosure BEFORE execution, all six gates first-run PASS,
then the adversarial review, then mechanical map application.

What the next me should know:
- The one-mode arithmetic is real and easy to re-derive (lambda_mean =
  2pi*R_H EXACT; ~4e-3 thermal quanta per horizon volume; a lab cavity
  holds 5.5e14). But do NOT quote L1 as deriving M=1 — ROUND 16's
  relabel is binding: regime = derived, M=1 = consistency condition
  anchored by 6Y's counting statistics. Patches to keep using: l >= 1
  horizon multipoles are gapped at sqrt(l(l+1))*H above omega_dress;
  the ~1e122 microstates are the reservoir setting ONE collective
  coordinate's occupation, not 1e122 partners.
- L2 is the round's keeper: kappa^2 = Omega^2 - Lambda*c^2 (sympy
  exact), threshold 1.449, SPARC measured floor 30.2, and the
  reviewer's OWN strengthening tau_B = 2pi/H (KMS time; our 1/H was
  conservative). Quote it as "the theorem gives the O(1) floor; the
  data give the 30-5000x freezing." The 6K/6M/6N analog failures are
  now RETRODICTED — never build a rate-based analog again.
- Trap fingerprint this round: a bars-mechanical verdict can overreach
  the lemma SENTENCE. My pre-reg said DERIVED-toy iff bars pass — the
  bars were about counts, the sentence claimed M=1. Next derivation
  pre-reg: pin the verdict letter to the CLAIM, not to the sub-checks.
- O5-detuning is the new sharpest mechanism question (his A9): the 1/2
  prefactor is only guaranteed at resonance; generic detuning
  (x_int - x_amb)*H/2pi would make it field-dependent; the sky
  measures 1/2 anyway (c1 at ~7 sigma, 7E 0.480). A detuning-resolved
  lending law beta(delta) is the successor instrument. Do not re-open
  rate-based analogs to answer it.
- Credence: bath-mechanism conditional ~8 HELD by the map's
  exactly-one-derived cell (L1 downgraded, L2 derived). No wiggle
  taken; the map's ambiguity clause was NOT needed. anomaly-real 53
  untouched all day.
- GD replication (TODO 28) is data-READY: Oh+2015 VizieR J/AJ/149/180
  first (26 dwarfs, ~8-10 slow, full baryonic decomposition),
  Iorio+2017 supplement. Next session can build the loader straight
  onto the 8S-c dial instrument.
- Round-16 file: REVIEW-ROUND16-OPUS.md, uncommitted like all of them.

## 2026-08-06 (night, +4) — P2 figures + ROUND 15 + the reference pass; both papers now referee-hardened

Figures: calcs/paper2_figures.py, seven gates PASS first run. The
pattern that worked again: parse the committed stage outputs, re-derive
their own printed headline numbers from their own tables as the gate,
assert every paper-quoted value, THEN plot. Nothing freehand; the 4T
osc+floor curve is the per-bin PRINTED values joined, not a re-fitted
model. Loader counts resolved a would-be confusion: 153 galaxies pass
cuts but only 149 contribute points (38 GD + 111 DD) — both numbers
are right, they count different things; the draft's 149 stands.

ROUND 15 (fresh Opus, same journal-referee brief): MAJOR-REVISION
light, ~40 checks, 38 PASS, figures 5/5, register certified. His two
real catches were both MY DIRECTION/GRADE SLIPS, the exact class the
round-14 lesson warned about (check orderings against the stage's own
verdict lines — I did it for orderings and still slipped on a drift
DIRECTION and an exclusion GRADE):
- "drifted high" in §5: the hier a0 drift is LOW (5D/5L: 0.79–0.99,
  f_ML-traded). I conflated the binary-translation row (which sits
  high) with the galaxy-hier rows. Trap fingerprint: when two ladder
  rows strain in opposite directions, name the row per clause.
- "quarter variant excluded once vertical is active": 5M has boot
  still 8.55 AHEAD dv-ON; the channel collapses its lead 76→9. The
  outright veto was the FENCED binary result, and 7J-d dissolved it.
  Trap fingerprint: "X died" claims must carry WHICH instrument
  killed X and whether that instrument survived the landed posterior.
Also adopted: Methods §2.2 (his "not self-contained" item — correct;
P1 had methods, P2 didn't), profile-vs-ledger distinction sentence,
discovery→anomaly softening (his argument from our own NO-credence-
movement stage banners — use that argument myself next time),
alpha defined, Newton +7.9 qualified full-sample. Draft → 0.3,
mean 21.6 words. Report archived REVIEW-ROUND15-OPUS.md (uncommitted).

Reference pass: Haiku scout + primary checks. P1's three flags
cleared. ONE catch: Desmond 2023 coordinates were Stiskalek &
Desmond's (525, 6130); the 0.034-dex paper is Desmond solo MNRAS
526, 3342 — abstract has sigma_int verbatim. Fixed in P2 + PAPER.md.
Same-author-same-year-wrong-volume = a phantom-citation subclass the
INSPIRE pass can miss; check the NUMBER against the CLAIM, not just
the name. Chae 921, 104 confirmed against our own archived PDF —
always prefer the on-disk primary over a scout's search.

State: P1 draft 0.4 (round-14-hardened, refs verified), P2 draft 0.3
(round-15-hardened, refs verified), both pushed. Next: (a) give the
referee a quick re-look at P2 0.3 or proceed straight to the joint
pre-circulation polish pass (both papers, one sweep: residual long
sentences, cross-references, abstract word counts); (b) then the
author decides circulation (CERN colleague first per the standing
path). Paper 3 stays deferred. Don't reopen: the figure scripts
(gated), the reference lists (verified 2026-08-06), the GD control
table (all eight numbers now referee-checked on disk).

CLOSE (+4, same night): the verification pass came back ACCEPT —
ten/ten resolved, cross-paper coherence clean, one one-word residual
adopted (symbolically→numerically; lesson: describe what the GATE
did, not what the math is). Sweep done in the same arc: abstracts
capped (P2 249w), final splits (mean 21.2), P1 given its reciprocal
companion cross-refs (it never pointed at P2 before — check
BIDIRECTIONAL cross-refs on any paper pair), worldtable six gates
green, CLAUDE.md handoff current. Both papers referee-accepted at
draft grade: P1 0.5, P2 0.4, all pushed. The paper block is CLOSED;
nothing paper-side left to self-start. Waiting on Filip: the
circulation call.

## 2026-08-06 (night, +3) — Paper 2 draft 0.1 written clean

Written fresh from v4.0 with the round-14 lessons from sentence one:
tables built in (data / c1 treatments / a0 ladder / GD controls),
zero prose em-dashes, first-draft mean 22.6 words (P1 started 28.8 —
the lesson took). 4,500 words. Structural independence HELD: no
mechanism language anywhere; the C&T identity framed as "fixes the
coefficients in advance"; the ledger summary names functions
descriptively (occupation / simple / quarter-coefficient variant);
the quadrupole section is pure phenomenology + the formulation fork.
The binary no-discrimination result (7J-d) framed as the field-level
methods finding per round 11. Next: P2 figures (5, gated — RAR+family
needs the SPARC loader + nu curves; c1 profiles + a0 ladder parse
from outputs; scatter figure needs care — parse 4T free-bin values,
do NOT freehand a scatter estimator; GD split re-uses the 9N engine),
then the referee round on P2 (tell Filip first), then the joint
pre-circulation reference pass on both papers.

## 2026-08-06 (night, +2) — figures built; ROUND 14 landed; draft at 0.3

Figures: calcs/paper1_figures.py, three provenance gates all PASS
first run (anchor 1.0784 vs 1.0779; census 9/2 from loader AND CSV;
9O regression). Fig 1 overlays = real GPU runs at the 7K-a printed
cells. ROUND 14 (fresh Opus, journal-referee brief): MAJOR-REVISION
light; 30/38 spot-checks PASS; endpoint discipline certified. HIS
TWO REAL CATCHES, both my drafting errors, both verified before
adopting: (1) I attached the eleven-pair 3.8e-9 to the nine-pair
census — P(≥9|0.9) = 4.8e-7; (2) "monotone dose curve in RUWE" is
false — Q1/Q2 inverted, and 8Z's own verdict line said
NON-MONOTONE all along. Lesson for Paper 2: when I compress a
ladder into prose, CHECK THE ORDERING against the stage verdict
line, and never move a probability between counts. All his
correctness/calibration items applied (0.3, commit 52784b4);
PAPER.md's two ambiguous dose clauses clarified. QUEUED for 0.4:
Table 1 cut ladder (needs a per-cut counting script), the
sentence-length pass (mean 27.9 → ~20 — my semicolon chains are
the em-dash workaround the contract also forbids; write SHORTER,
not differently-punctuated), the phantom-veto figure, full script
filenames in App B. Then Paper 2. The referee agent is still
alive this session for follow-ups (SendMessage).

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
