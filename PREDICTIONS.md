# PREDICTIONS.md — the signed prediction ledger

Registered predictions of the thermal-horizon reading (the ambient-gated bath and its
parents), separated from postdictions, with kill conditions stated in numbers.

**Rules.** (1) Every live prediction is registered by git commit BEFORE its test
executes; the commit hash is its timestamp. (Terminology: this is *bar-locking* —
committed thresholds with knowledge of all previous stages. It guards against post-hoc
bar-moving and outcome-shopping, and several locked bars have fired against us; it is
not prospective pre-registration and does not control the sequential-search error rate —
see §0.) (2) Entries are never edited after
registration — status flips only (LIVE → CONFIRMED / KILLED / MOVED, with a pointer),
mirroring [LEDGER.md](LEDGER.md) discipline. (3) "Parameter-free" means no dial of the
rule can absorb a failure; every kill condition is quantitative. (4) Postdictions are
labeled as such: they are consistency records, not predictions. (5) Uniqueness class is
stated per entry: **DM-immune** (no dark-matter effective theory produces it),
**classical-immune** (no classical-statistics bath produces it), or **MOND-generic**
(shared with classical MOND — listed for the record, claimed for discrimination only
against Newton/DM).

---

## 0. The frozen package (registered 2026-07-25, this commit)

Adopted from an external review solicited by the author (2026-07-24; the point-by-point
response is logged in NOTES Stage 7I): the successive-refinement construction of the
function family — each step motivated by the previous round's results on the same data —
is legitimate exploration but cannot certify its own survivor. **The function search is
therefore CLOSED on the in-sample data as of this commit.**

**Frozen function list** (exact forms, all previously registered; y = g_N/a₀):
- **Newton** (control): ν = 1.
- **simple-ν** (control; the classical self-consistent bath, 4F): ν = ½ + √(¼ + 1/y).
- **BE** — the occupation law (Cadoni & Tuveri 2019): ν = 1 + n_BE(√y).
- **gm** — the β = ½ exchange-symmetric fixed point (5J): ν = 1 + n_BE(y^¾ √ν).
- **AMB** — the derived ambient-gated bath (6E/6U, post-hoc flag carried):
  β = ½·[q_loc·s_amb]², q_loc = 1/(2ν−1), s_amb = n_amb/(1+n_amb),
  n_amb = n_BE(√(e_N/a₀)); solver form u = y^((1+β)/2)·ν^β; trajectory formulation
  per 7G/7H.

**Frozen conventions:** a₀ = 1.2e−10 m s⁻² fiducial (fitted where stated); binary
external field = Newtonian g_N,ext = 1.2 a₀ (3S/3T RAR inversion); QUMOND EFE tables
from [calcs/qumond_efe_solver.py](calcs/qumond_efe_solver.py); galaxy primary
likelihood = the vertical-hardened hierarchical ladder
([calcs/stage5m_hierv.py](calcs/stage5m_hierv.py) machinery); binary primary
likelihood = the corrected-velocity six-seed v7 budget
([calcs/stage3p_v7budget.py](calcs/stage3p_v7budget.py) + the 4R data patch, seeds
31/101/202/303/404/505, g1p2 tables).

The in-sample function comparisons carry an uncorrected sequential-search multiplicity
(~30 forms tested along a data-guided path); no post-hoc trials correction can repair
that, and none is attempted — closure plus out-of-sample routing is the repair.

**The rule.** (1) No further functional-form iteration is scored on the in-sample data
(SPARC 153 + its hierarchical treatments, EDR3 14,071 + the ceiling census, the Chae
fields, the lensing set) *as evidence*. Theory work continues; a newly constructed form
may be registered here with its predictions, but its in-sample fit joins the record as a
consistency row only. (2) Function-level claims from here on route through out-of-sample
tests (Section C; the anchored subsample; held-out selections — Cookson-catalog overlap
check pending; DR4). (3) The world table (Stage 6Q) is the closed in-sample record of
the thirteen laws tested.

**Primary analyses (multiplicity control, same review).** Each headline claim has ONE
pre-declared primary treatment; every other treatment is a sensitivity row. Galaxies:
the vertical-hardened hierarchical ladder. Binaries: the six-seed corrected-velocity v7
budget. Scatter channel: outer-disk / anchored subsets only (7C); the mode-count N is
not identifiable at SPARC depth (4W) and is retired from the claims list (the measured
x-shape and its environmental control, 4T/5B, stand). Within-program grading language
("largest lead", "top mark") is notebook vocabulary, not evidence — capped in PAPER.md.

---

## A. The record (postdictions and survived kill tests)

| item | value | uniqueness class | stage |
|---|---|---|---|
| Wide-binary velocity boost | 1.086 (CI 1.064–1.110); perspective-corrected 1.078 (1.052–1.103) | DM-immune (nothing dark clusters at 10 kAU; required local density ~10⁶× the DM background) | 2C/2D, 4Q |
| EFE amplitude demanded by binaries | no-EFE MI dead 12/12 (−20…−28 lnL); data demand the external-field suppression | DM-immune (strong-equivalence-principle violation; no halo mimic) | 4L |
| Perpendicular-velocity ceiling census | 11 pairs in the Newton-forbidden band [√2, 1.67); cliff at the boosted edge √(2·1.36) ≈ 1.65; leakage null P = 3.8e−9 | DM-immune | 4J/4M |
| Exponential (Boltzmann) screening + half-quantum NLO | ĉ₁: galaxies 0.26–0.45 (hier–flat range), binaries 0.37–0.50; c₁ = 0 excluded everywhere; ½-branch over standard-μ 198–200/200 raw, Δ−2lnL −56 honest | **classical-immune** (the classical self-consistent bath — simple-ν, Stage 4F — has a power-law tail: no exponential rungs, no ½ℏω share) | 4A/4B, 4S, 4X, 4Z |
| Two-system tail postdiction | p_gal = 0.6884 (measured 0.65–0.75); p_bin = 0.5280 (measured band [0.45, 0.60]) | postdiction; gate-specific | 6E/6G/6U |
| Temperature identity a₀ = cH₀/2π | galaxies +0.1σ (4V); all hier fits return a₀ = 1.04–1.13e−10 (5M); best binary row (AMB) 1.28 ± 0.15 = +1.6σ | MOND-generic as a number; the *identity* is horizon-specific | 4V, 5M, 6G |
| Radial-orbit excess w_rad = 0.20 | matches Hwang+22 superthermal law (20–22% at e > 0.9) — external cross-validation | MOND-agnostic consistency | 3N/4G |
| Untied-exponent falsifier | tied gate (2,2) survives the untied contest (margins within noise) | — | 6V |
| Cassini screening bound | p > 0.234 required; measured 16th pct 0.462 — passed | — | 3 (RAR fit) |
| Vacuum share in the gate's local factor | classical-pull replacement (no "+1") rejected at +556.5 on the galaxy ladder (cap-robust lower bound; worse than no admixture at all); binary direct leg instrument-N/A, indirect support via 6G acceptance + β < 0.03 | **classical-immune** (the "+1" is spontaneous emission — absent from any classical field bath; conditional on the admixture grammar) | 7D |

**7J annotation (2026-07-25; rows above not edited, per Rule 2):** Stage 7J measured
the companion completeness (flag C = 0.410; true host multiplicity f_host 0.42–0.57,
Raghavan-consistent) and ran the multiplicity-marginalized fit with the
photocenter-corrected amplitude law. The pre-registered verdict fired **COMPANION-WIN**:
at the measured host rate the forward likelihood prefers Newton + companions
(α_marg = 0.00 both laws, both seeds; NOT an f_pm-grid artifact — unchanged at
f_pm = 1.5; the winning cell itself misfits by −60…−116, so no cell fits both the
kinematics and the multiplicity). The quadrant map registered at d2dc7eb resolved to
**quadrant C**; the pre-committed credence move executed: anomaly-real → ~35%.
Consequences for this table: the three binary rows (boost, ceiling census, EFE
amplitude) are model-light *measurements* whose gravity attribution is SUSPENDED
pending 7K-a (forward median at the 7J winning cell) and 7K-b (census leakage null
under measured multiplicity); correction #15 retracts the 7I cut-survival immunity
inference. Galaxy rows untouched. Full record: NOTES Stage 7J; ledger rows
bin-7j-completeness / bin-7j-marginal / bin-7j-rawcond / ret-7i-median-immune.

Corrections record: eighteen logged in [PAPER.md](PAPER.md) Appendix A. In-sample
environmental ordering WITHDRAWN (correction #14, Stage 6Z) — moved to C/P2.
7I cut-survival immunity inference RETRACTED (correction #15, Stage 7J; #18
moves it to undetermined-pending-repair). §7.3 "manufactured significance"
claim RETRACTED (correction #16, Stage 7J-x; one clause re-graded under #18).
Prior-axis exemption from our own AMBIGUOUS criterion FIXED (correction #17,
Stage 7J-y): the α_marg(prior-anchor) curve is the primary fit-level statement,
labels are annotations. **Part A's absolute scale RETRACTED (correction #18,
Stage 7J-y): within-pair overluminosity correlation ρ ≈ +0.47 = common-mode
over-attribution; C biased low, f_host biased high; literature subsystem rates
(~0.2 per pair, Tokovinin-class, primary-source verification pending) now
primary; the 7J verdict label at that anchoring is AMBIGUOUS, decided by the
pre-registered six-realization budget; the anchoring-invariant durable finding:
the fenced +99–110 Newtonian rejection is dead at every defensible companion
prior.** 7J standing additions: the per-arm injection standard (amendment 6,
pre-registered at 0045eea) with its committed downgrade rule; the 7J-y
discriminant pre-registered at 3f7212e and fired against our own Part A.

---

## B. Live — in-sample, executable now

### P6 — the Einstein particle term in the RAR scatter (γ = 1)
*Registered 2026-07-24, this commit, BEFORE execution of Stage 7A.*

**Statement.** If ν − 1 = n̄ is a thermal occupation with **quantum** statistics, its
fluctuation law is fixed by Einstein's 1909 formula, Var(n) = n̄ + n̄² (particle +
wave), with zero new shape freedom. Through δlog g = δn/((1+n̄) ln10) this forces the
intrinsic RAR variance to be a pure Boltzmann factor in x = √(g_N/a₀):

- **quantum (Bose): s²(x) ∝ e^(−x)** — i.e. γ = 1 in s² = S₀e^(−γx),
- classical wave bath (Var = n̄², the Rayleigh–Jeans/continuous-field limit; also the
  fluctuation law of the 4F classical self-consistent bath): **γ = 2**,
- pure corpuscular/shot bath (Var = n̄): s² ∝ e^(−x)(1 − e^(−x)) — dies in the deep
  limit (distinct shape, not a γ).

The a₀ inside x is LOCKED to the mean-law fit (fluctuation–dissipation: one a₀, both
moments). Same parameter count for all rivals (amplitude + floor).

**Uniqueness.** Classical-immune by construction: a classical effective ν(y) can mimic
any *mean* curve, but its self-consistent bath fluctuates with the wave term only —
γ = 1 requires the particle term, i.e. quanta. DM predicts no universal x-organized
intrinsic scatter at all (scatter organizes by formation history), and no mean↔variance
lock. This is the same argument by which Einstein 1909 first inferred light quanta from
fluctuation statistics — pointed at the radial-acceleration relation.

**Prior state.** Stage 4T fit ONLY the quantum shape (beats constant by −25.1 at 1
extra param; N̂ ≈ 21 over a 0.101-dex floor); the rivals were never contested; γ was
never measured. 4U/4W caveat carried: the thermal amplitude deflates under M/L and
vertical marginalization — this test is 4T-grade (per-point), calibration-gated, and
any SUPPORT stays one-channel until hier-hardened.

**Test.** [calcs/stage7a_einstein.py](calcs/stage7a_einstein.py) (same commit):
free-γ profile + three fixed laws at matched complexity; instrument calibrated by
paired injections (quantum truth vs classical truth — verdicts quoted only if the
recovered γ̂'s separate); x ≈ 1 bump handled by mask and by explicit bump component
(the 4W point-level bump is a known non-monotone excess).

**Kill condition (pre-committed).** If γ = 1 is excluded at Δ−2lnL ≥ 9 with the
calibration gate passing and the exclusion surviving the bump treatments → the
quantum-statistical reading of the occupation takes a REAL strike (the mean-law
function record is untouched; what dies is quantum statistics in the scatter channel,
and the "quantum bath" language gets demoted program-wide). Symmetrically: γ = 2
excluded at ≥ 9 under the same gates = the particle term is present = the first
in-sample uniquely-quantum signature. Both directions pre-stated.

**STATUS (flipped 2026-07-24, Stage 7A executed same day): LIVE — UNRESOLVED at
per-point (4T) grade.** The instrument is validated (calibration gate PASSED: paired
quantum/classical injections recovered at γ̂ = 1.03 / 2.05), the corpuscular rival is
EXCLUDED (+25.05, collapses onto the floor), but quantum-vs-classical returned no vote
(EQ−EC = −0.43; γ profile bimodal, whole axis ≤ 4.4 units). Structural reason
identified: the two laws differ maximally at x ≈ 0.4–1.2, and the point-level x ≈ 1
bump (4W: survives every marginalization) occupies the discriminant window. Unlock:
bump-source identification/subtraction, a hier-hardened γ instrument on
distance-anchored data, or DR4-era samples. Both kill directions remain pre-committed.
[data/stage7a_einstein.txt]

**STATUS 2 (flipped 2026-07-24, Stages 7B/7C — the bump hunt): LIVE, MOVED fully
out-of-sample.** The bump was identified in locality: it is an inner-disk (R < 1.5 R_d)
phenomenon — its fitted amplitude on outer-disk points is exactly zero, and the
radius-mix shift explains two-thirds of the window excess (PARTIAL vs the 0.75 bar;
candidate drivers scout-anchored: non-circular streaming 10–40 km/s, beam smearing,
decomposition — inner-disk astrophysics, not the law). With the inner term modeled, the
γ profile de-bimodalizes (interior γ̂ = 0.52, γ = 2 nominally +12) BUT the calibrated
injection gate FAILS on the informed design, the clean-subset variant leans the
opposite way (+4.5), and 7B measured within-curve residual coherence ρ ≈ 0.87 — the
quantum-classical difference sits below SPARC's correlated-systematics floor.
**Verdict: γ is NOT MEASURABLE on SPARC at any grade tried; neither pre-committed kill
direction fired. P6 now awaits anchored-distance / DR4-era / IFU samples with explicit
non-circular modeling.** [data/stage7b_bumphunt.txt, data/stage7c_gammaclean.txt]

---

## C. Live — out-of-sample (the DR4-era ledger)

### P1 — the tail ceiling p ≤ ¾ (parameter-free)
The gate exponent p = ½ + g/4 with g = s_amb² ≤ 1 caps every system's tail index at
¾; void galaxies asymptote ≈ 0.72, never beyond. **Kill:** one clean void galaxy with
tail index beyond ¾ at 3σ. (Stage 6Y; ledger P1.)

*Annotation 2026-07-31 (Stage 8C, pre-reg a5b4816 — first in-sample ceiling
contest; no status change): verdict AMBIG by the locked grammar with the BREAK
side clean — pooled p = 0.617 ± 0.133 (bootstrap; profile 0.647, Δ1 0.628–0.667,
p = 0.75 at +9.7), ZERO calibrated exceedances (40 quotable galaxies;
false-positive calibration clean on two p_true = 0.72 null skies); the
per-galaxy census is POWER-LIMITED at SPARC grade (0/34 exceedances fired even
at injected p_true = 0.90 — the literal one-galaxy kill clause needs
anchored/DR4-era data); the eN-HIGH arm sits nominally above (0.904,
0.79–1.03 = 1.3σ arm-local, inside the locked bars; axis 6Z-unreadable, no
ordering claim either direction). P1 stays LIVE; kill conditions unchanged.*

*Annotation 2026-08-07 (Stage 9S, pre-reg 4b46cd5; no status change): the
ceiling test doubles as a DIRECT meter of the exchange prefactor r (the
ROUND-16 detuning object): in the void asymptote g → 1 the tail reads
p_void = ½ + r/2 with no gate uncertainty — r = ½ ⇔ ¾. The in-catalog
bound r ≥ 0.315 (9S, hier-anchored) implies p_void ≥ 0.66 within the AMB
form; a clean void measurement below that would strike the resonance
reading before it ever touches the ¾ ceiling itself.*

*Annotation 2026-08-07 later (Stage 9T + ROUND 17; no status change):
conditional on the averaging flank (O5-AVERAGING), the resonance theorem
predicts the void asymptote SATURATED: p_void ∈ [0.727, 0.750] (edge
0.718 for z_form = 10). The fiducial-gate kill stands: a hier tail
measurement with σ_p ≤ 0.02 demanding p < 0.67 kills the resonance
reading. NEW falsification lever (ROUND-17): if the exchange has
completed only ~one Rabi radian (γ ~ 1, the un-averaged regime), the
tail should read p ≈ 0.587 — a γ-measurement decides between ~0.59 and
the ~0.69 band; the measured 0.65 (grid-limited) sits between them.*

*Annotation 2026-08-07 close (Stage 9U, the γ-meter; no status change):
the in-catalog γ-measurement is POWER-LIMITED — p̂ = 0.6471 ± 0.0746
(vertical-hardened primary, 40-rep galaxy bootstrap, SD clipped at the
instrumented window: honest σ_p ≥ 0.075 vs the 0.02 requirement). All
three bands sit inside 1σ (one-swing +0.81σ, floor −0.32σ, full-avg
−0.55σ): the ROUND-17 kill did NOT fire, and neither reading is
vindicated — SPARC-hier cannot resolve the swing count (per-rep tail
wells ~7–10 lnL; location realization-dominated). Calibration
propagated: the 9S in-catalog bound r ≥ 0.315 is point-conditional
(error-calibrated r > 0 at ~2σ); p_void ≥ 0.66 above softens
accordingly. CONSEQUENCE: this void channel is now the LIVE r-meter —
gate-free (the signal is not divided by g), with the kill-band at the
ceiling — and carries the O5-AVERAGING question alone until anchored/
DR4-era tails or a void-RAR sample exist.*

*Annotation 2026-08-08 (Stage 10B + ROUND 23; REGISTERED conditional
kill-test, reviewer-promoted): under the DISPERSIVE reading of the
grammar (10B: real-exchange carriers excluded C1–C7; the constraint
vertex is dispersive-class by the exact selection rule; reading
PERMITTED-not-forced), the exchange weight is structurally FIXED at
r = ½ with no γ-running — the void asymptote is p_void = ¾ EXACTLY
(not a band), and the tail carries no swing-count dependence. KILL
(dispersive reading only): a void/DR4-era measurement of r < ½ at
> 2σ falsifies the dispersive realization while off-resonance
exchange would survive; conversely p_void = 0.750 ± small with no
γ-running signature is the dispersive reading's signature. The
resonance-reading band [0.727, 0.750] and the parameter-free p ≤ ¾
ceiling are unchanged.*

*Annotation 2026-08-08 (Stage 9Y + ROUND 19; HARDENING, no status
change): the void kill-band [0.727, 0.750] is STATE-INDEPENDENT —
squeezed-thermal and thermal ambient states converge at the void
asymptote (both gates → 1; reviewer-reproduced 0.9678 vs 0.9675 at
n̄ = 60, r = 0.5). The P1 test does not depend on the quantum state of
the soft sector. Separately REGISTERED (conditional): the ambient
NON-THERMALITY dial — super-thermal ambient statistics (quantum
squeezing OR classical thermal mixtures; degenerate on this scalar)
raise the binary-side tail by Δp_bin ≈ +0.006/+0.013/+0.032 at
r_sq-equivalent 0.2/0.3/0.5 under the g = P(n≥2) gate reading (the
reading is a named conditional; [P(n≥1)]² splits ~50% from it under
squeezing). A DR4-era weak-ambient tail at σ_p ~ 0.01–0.02 reads this
dial; today's bound is vacuous at 2σ (9U σ_p ≥ 0.075).*

*Annotation 2026-08-07 second close (Stage 9V; no status change): the
direct in-catalog r fit (per-galaxy measured gates) confirms and closes
the arc: r̂ = 0.3365 ± 0.1869 (clipped bootstrap; fid treatment 0.3895
= the 9U conversion 0.3904 at d = 0.001 — cross-instrument coherence);
r = 0 disfavored ~31 nominal −2lnL units in every gate treatment but
only +1.8σ at honest grade; the r_true = ½ injection reads 0.407
(biased low, disclosed — full averaging is NOT more disfavored than
face value). VOID-SAMPLE FEASIBILITY: the only public void RC sample
(Pustilnik+2020, 8 dwarfs, V_max 31.5–80.3 km/s) has no Newtonian arm
— this prediction's test requires void galaxies WITH y ≳ 1 coverage
(WALLABY DR2 + environment tags, or DR4-era). The prediction stands
armed and untested.*

*Annotation 2026-08-09 (Stage 10I, pre-reg 630516e + ROUND-30 rescope;
no status change): the LV-ISOLATION route to this test is CLOSED at
public grade — instrument-null, channel untested. The pre-registered
null-axis gate fired first-run: the point-source neighbor-sum meter
has NO dynamic range at the isolated end (iso neighbor terms ≤ 1.7%
of any LSS floor; the meter is asymmetric — at the group end it has
real signal, NGC2976 = 110% of floor), and Chae's LSS-aware e_N puts
LV galaxies, isolated and not, already near-open (0.005–0.009 a₀,
gates s² = 0.825–0.869) — LV tidal isolation does not reach the void
regime (g → 1 needs e_N ≪ 0.005). Census (an UNDERCOUNT — the
matcher misses UGC-designated aliases, exemplar UGC05721 = NGC3274,
isolated, letter alias-robust at ρ → −0.126): ≥3 isolated
transition-crossers in SPARC × UNGC. The test's operative data
requirement sharpens to: DENSITY-FIELD-selected (true void) galaxies
with y ≳ 1 coverage — neighbor-isolation catalogs cannot substitute.
Inverted dividend: the meter's powered end is the GROUP end — a
gate-CLOSING (p → ½) test toward dense environments is the axis this
data supports. No sky read was performed; r̂ = 0.3365 ± 0.1869 (9V)
stands operative.*

### P2 — the environmental ordering p(e_N)
Weaker ambient → sharper tail (p rises toward 0.72–0.75 in voids; galaxy-population
median 0.689). In-sample version WITHDRAWN (correction #14: gate heterogeneity, not
ordering — 6Z). **Out-of-sample failure on environment-resolved DR4-era samples is a
REAL strike — pre-stated.** (Stages 6I/6Z.)

### P3 — the z-locked pair: a₀(z) AND p(z) run together
The temperature identity forces a₀(z) = cH(z)/2π (parameter-free leg); the gate then
forces the galaxy tail index to run with it (environment-conditional leg, fiducial
galaxy ambient n_amb = 6.6 held comoving):

| z | H(z)/H₀ (flat ΛCDM 0.3/0.7) | a₀(z)/a₀ | p_gal(z) |
|---|---|---|---|
| 0 | 1.00 | 1.00 | 0.689 |
| 0.5 | 1.31 | 1.31 | 0.695 |
| 1 | 1.76 | 1.76 | 0.702 |
| 2 | 2.97 | 2.97 | 0.712 |

**The signature is the LOCK** — one function H(z) moves two observables in fixed
ratio. Classical MOND holds a₀ constant; DM holds the RAR normalization roughly
constant; nothing else ties the two drifts. **Kill:** measured a₀(z) flat at z ≳ 1
at 3σ, or the two drifts measured with incompatible ratio. (Stages 4V/6U; data era:
high-z rotation curves with pressure-support control.)

*Annotation 2026-07-29 (8A scout round; post-registration bookkeeping, no status
change): Magneticum (Mayer+22, arXiv:2206.04333, scout-level) reports
ΛCDM-simulated a₀ RISING ≈×3 by z = 2.3 — approximately the H(z) trajectory
(×2.97 at z = 2). The a₀(z) leg ALONE is therefore not ΛCDM-immune, and the "DM
holds the normalization roughly constant" clause above is too strong at
simulation grade. P3's discriminating content is the LOCK — the paired p_gal(z)
drift in fixed ratio — exactly as registered; the kill conditions are unchanged.
Primary read of Mayer+22 owed before any use.*

*Annotation 2026-07-31 (8E primary reads; no status change): the a₀(z) data era
has OPENED. MUSE-DARK III (arXiv:2604.22613, primary-read): a₀ RISES — binned
1.99→2.71e−10 over 0.33 < z < 1.44; linear a₀(0) = 1.00±0.04e−10 (essentially ON
cH₀/2π = 1.08), a₁ = 1.59e−10/z = ~2× the lock's 0→1 secant (+0.82e−10) at face
value, with their own M/L systematic budget (+0.2–0.45 dex) spanning the
difference and their linear form flagged phenomenological. Limbach+09's
"cH(z)-coupling excluded" DEFLATES on primary read (formal-errors-only,
self-disavowed at systematics level in their own conclusion). Shachar+23 BTFR
stays flat to z≈2.5 — the field is internally contradictory. The registered KILL
("a₀(z) flat at z ≳ 1 at 3σ") is NOT triggered; the LOCK (p_gal(z) paired drift)
remains unmeasured by anyone and is still the unique discriminant. NOTES
lit-note 2026-07-31 carries the full numbers.*

*Annotation 2026-08-09 (re-read with digit verification; no status change): the
lit-note's "their MOND-framework refit agrees within errors" clause CORRECTED —
the slope is fitting-framework-conditional: DM decomposition a₁ = 1.59 ± 0.10
(95% CI, verbatim) vs their own MOND-framework fit a₁ = 1.20 ± 0.10 (their Eq.
E.5) vs ΛCDM-halo variants 1.63. Lock linearized over their window = 1.04–1.09;
vs the MOND-framework row that is 1.6σ (if 1σ) / 3.1σ (if 95%) — the
apples-to-apples comparison is NOT far from the lock, and the framework split
itself evidences form-coupling of the recovered slope (the lens the LOCK
predicts: p_gal drifts too, their form is held fixed). Successor sharpened to a
two-leg contest (lock-vs-linear on released per-galaxy products + a lens-bias
leg on locked-pair truth); their RCs are downloadable (their footnote 3). NOTES
addendum 2026-08-09 carries the digits.*

*Annotation 2026-08-09 (stage 10H executed; no status change): the two-leg
contest was BUILT and pre-registered (546641f + three amendments) and closed
H-FEASIBILITY-LIMITED — the public release underdetermines their slope
(per-point errors, M_HI/gas extent, and the 'regular' flags are unreleased;
their gas text is internally ×2-and-curvature ambiguous), so neither leg was
read; the firewall held. What did land: the z∼1 LEVEL reproduces from the
release (report-grade, a₀ ≈ 2.2–2.6), and at their own Planck-2015 cosmology
the lock intercept cH₀/2π = 1.047e-10 lands on their MOND-framework and
per-galaxy rows (0.7σ / 0.1σ) and 2.3σ — just outside 95% — from their DM
headline row (ROUND-29 phrasing). The flip criteria are untouched —
model-conditioned tracks are not clean-grade external data. Reopeners: the
author-contact route (Desmond is a co-author; outreach list), galpak source
read, MIGHTEE z<0.08, DR4-era. Ledger ext-10h-musedark-recon.*

### P4 — the DR4 weak-ambient pair
Weak-ambient wide binaries sharpen toward p ≈ 0.69, and the source-vs-dressed
convention split becomes resolvable (Δp ≈ 0.025 at e_N = 0.4). (Stage 6U.)

### P5 — the nested-ambient rule (M = 1)
The gate reads ONE collective ambient mode — the largest relevant environment's; ambients
do not compose pointwise. Cluster-member vs field galaxies at matched g_N and local e_N
differ per the parent ambient. (Stage 6Y exclusion theorem.)

### P7 — the MI eccentricity signature (the Saturn kill test)
The reservoir coupling is trajectory-state (collective-mode = EFE-respecting MI class,
6Y corollary): eccentricity-resolved DR4 binaries separate MI from MG where 4L ties.
MG-with-quadrupole predicts the Cassini violation stands (4.0–5.8× for every member,
5S); MI predicts Saturn clean + e-dependent binary residuals. **Kill (of the MI door):**
DR4 e-resolved fits reproduce the MG shape with no e-dependence beyond MG's. (4L/6Y.)

### P8 — the c₄ rung sign flip
At L = 2 the fourth ladder coefficient is c₄ = s²/192 − 1/720: sign flips with ambient
strength (positive for s² > 192/720 ≈ 0.267, i.e. all but the deepest-ambient systems).
Beyond current population grade (6L: the deep arm of 153 galaxies cannot read c₃⁺
rungs); listed for the DR4+/LSST era. (Stage 6H.)

### P9 — the DR4 width-object discriminator (instrument-level; registered before DR4 exists)
The v7/7J sky demands two width objects: the noise-scale chase (P(fpm = 3.0) =
0.54/0.97 — 3× formal PM errors, past the ~1.4× validation ceiling) and the per-system
scatter sq = 0.2. The 8F-c fc10 world (registered 604a759, executed f97e712) shows a
~10% KT = 4 error-tail population reproduces the chase with α uncorrupted — while the
7J-z6 sky-side rejection of an explicit tail axis (+0.0) argues against the literal-tail
reading. THE DISCRIMINATOR: if the width object is an ERROR TAIL, DR4's improved
astrometry quenches the fpm chase AND the sq demand TOGETHER (one mis-modeled noise
object); if it is ASTROPHYSICAL (unmodeled orbital/population broadening), sq persists
at ≈ 0.2 while the fpm posterior normalizes toward ≤ 1.4×. Registered NOW, before DR4.
Uniqueness class: instrument-level (not a gravity claim — either outcome leaves the α
measurement standing per 8F/8F-c: manufacture excluded for symmetric variance-inflation
tails; the BE coupling bounded below ~15% severity). (Stages 8F/8F-c; this commit.)

---

## D. What we do NOT claim as unique

- Any single mean-law shape ν(y) — mimicable by construction (a phenomenologist may
  simply postulate the curve).
- The boost amplitude alone — MOND-generic.
- Thermal *bunching* alone — the wave term n̄² is fully classical (thermal light
  bunches classically); the quantum claim rides ONLY on the particle term (P6).
- Laboratory signatures — our own cancellation theorem (6N) says buildable flat-bath
  configurations measure zero; we do not predict a tabletop detection.
