# PREDICTIONS.md — the signed prediction ledger

Registered predictions of the thermal-horizon reading (the ambient-gated bath and its
parents), separated from postdictions, with kill conditions stated in numbers.

**Rules.** (1) Every live prediction is registered by git commit BEFORE its test
executes; the commit hash is its timestamp. (2) Entries are never edited after
registration — status flips only (LIVE → CONFIRMED / KILLED / MOVED, with a pointer),
mirroring [LEDGER.md](LEDGER.md) discipline. (3) "Parameter-free" means no dial of the
rule can absorb a failure; every kill condition is quantitative. (4) Postdictions are
labeled as such: they are consistency records, not predictions. (5) Uniqueness class is
stated per entry: **DM-immune** (no dark-matter effective theory produces it),
**classical-immune** (no classical-statistics bath produces it), or **MOND-generic**
(shared with classical MOND — listed for the record, claimed for discrimination only
against Newton/DM).

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

Corrections record: fourteen logged in [PAPER.md](PAPER.md) Appendix A. In-sample
environmental ordering WITHDRAWN (correction #14, Stage 6Z) — moved to C/P2.

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

---

## C. Live — out-of-sample (the DR4-era ledger)

### P1 — the tail ceiling p ≤ ¾ (parameter-free)
The gate exponent p = ½ + g/4 with g = s_amb² ≤ 1 caps every system's tail index at
¾; void galaxies asymptote ≈ 0.72, never beyond. **Kill:** one clean void galaxy with
tail index beyond ¾ at 3σ. (Stage 6Y; ledger P1.)

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

---

## D. What we do NOT claim as unique

- Any single mean-law shape ν(y) — mimicable by construction (a phenomenologist may
  simply postulate the curve).
- The boost amplitude alone — MOND-generic.
- Thermal *bunching* alone — the wave term n̄² is fully classical (thermal light
  bunches classically); the quantum claim rides ONLY on the particle term (P6).
- Laboratory signatures — our own cancellation theorem (6N) says buildable flat-bath
  configurations measure zero; we do not predict a tabletop detection.
