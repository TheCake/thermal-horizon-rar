# DR4.md — the Gaia DR4 day-one plan (v2, 2026-08-12)

**STATUS: ADOPTED 2026-08-12** (the author's free-hand session). Supersedes the
2026-08-06 v1 in place (v1 = git history). Produced by an adversarial red-team
round (findings F1–F18, report REVIEW-DR4-OPUS.md, uncommitted per convention)
plus a build session that dry-ran the day-one ADQL against the live DR3 archive
(11/11 probes) and verified the DR4 facts against ESA at source; reviewed in
full and merged by the main session.

**Nothing in this file is pre-registered.** Every threshold below is marked
**PROPOSED** and requires its own pre-registration commit, with bar-locking, in
advance of 2026-12-02. The red team's finding F4/F5 was precisely that DR4.md
claimed bars existed where they did not; this draft does not repeat that.

**No credence numbers are invented here.** The credence map (§9) is a skeleton
with empty outcome cells. Signing it is a DEDICATED pre-DR4 session task (it
gates every pre-registration commit and must land well before 2026-12-02); it
was deliberately NOT signed in the adopting session, which was at the end of a
long multi-arc day — pre-signing deserves fresh judgment, not leftover budget.

---

## 1. What this document is

DR4.md's purpose was to freeze WHAT gets tested, in WHAT order, so the first
week is mechanical. The red team's verdict was that the frozen list is
day-one executable for **zero of eleven** rows as written, and that its largest
single omission is the one thing DR4 uniquely delivers: a **direct measurement
of the hidden-companion population** that the whole binary chapter currently
infers photometrically.

This revision therefore does four things.

1. Puts the companion census first, as **Test 0** (F8), with its
   false-positive calibration requirement (F9).
2. Separates **Gaia-native** tests from **DR4-era, non-Gaia** tests (F2), and
   marks one row **NOT-GAIA-FEASIBLE** with the arithmetic attached (F3).
3. Gives every surviving Gaia-native row an exact table/column list, an ADQL
   sketch dry-run against DR3, a named statistic, and a PROPOSED bar (F4, F10).
4. Adds the guards the program's own standing rules require and the frozen plan
   omitted: the sky-blind firewall (F5), the solution-heterogeneity
   stratification (F7), the trap-#26 known-outcome scoring (F11), the loader
   identity regression (F16), and the credence-map skeleton (F14).

Companion artifacts produced with this draft:

- `calcs/dr4_dryrun_adql.py` — the day-one ADQL, dry-run against the public
  Gaia DR3 TAP service (11 probes, all returned).
- `data/dr4_dryrun_adql.txt` — the banked query text with DR3 row counts as the
  day-one regression target.

---

## 2. External facts, verified against ESA (2026-08-12)

Where ESA and the red team differ, **ESA is followed** and the difference is
noted.

| Fact | Value | Source |
|---|---|---|
| Release date | **2 December 2026** | [ESA DR4 page](https://www.cosmos.esa.int/web/gaia/data-release-4) |
| Data span | 25 Jul 2014 10:30 UTC → 20 Jan 2020 22:00 UTC, **66 months (5.5 yr)** | [ESA DR4 content](https://www.cosmos.esa.int/web/gaia/dr4) |
| Total volume | ~400 TB | ESA DR4 content |
| Sources processed | ~2.8 × 10⁹ | ESA DR4 content |
| `gaia_source` | **2,000,033,276 rows, 586 GB** — the *high-quality subset*, "high-quality astrometry and high-quality photometry" | ESA DR4 content |
| Model consolidation | `gaia_source` parameters are "consolidated from several different processing modules such that the best-suited model has been applied for a given source (e.g. a binary-star model instead of the default single-star model as was the case in Gaia DR3)" | ESA DR4 content (verbatim) |
| All-source tables | `all_source_astrometry` (610 GB), `all_source_photometry` (755 GB), `all_source_rvs` (28 GB), `all_source_flags` | ESA DR4 content |
| NSS | **10 tables**, named ones include `nss_two_body_orbit`, `nss_acceleration_astro`, `nss_resolved_pair`, `nss_multiple_orbits`, `nss_masses`, `nss_non_linear_spectro`, `nss_multiplicity`, `nss_vim_fl`, `optical_pair` | ESA DR4 content |
| Epoch data | `epoch_astrometry` **62 TB**, `epoch_photometry` 7.2 TB, `epoch_photometry_ccd` 23 TB | ESA DR4 content |
| RVS | `rvs_mean_spectrum` 1.8 TB, `rvs_epoch_spectrum` 49 TB, `rvs_epoch_parameters_single/_double` | ESA DR4 content |
| RUWE | **not listed on the DR4 content page** | ESA DR4 content |
| DR3 NSS false positives | A software error found *during DR4 epoch-astrometry validation* invalidated four DR3 `nss_two_body_orbit` astrometric-orbit solutions: `4698424845771339520` (WD 0141-675), `5765846127180770432` (HIP 64690), `522135261462534528` (54 Cas), `1712614124767394816` (HIP 66074); three corresponding `binary_masses` rows are void. `gaia_source` parameters for all four remain valid. | [DR3 known issues](https://www.cosmos.esa.int/web/gaia/dr3-known-issues) |
| Residual images | new DR4 product, ~7 × 10⁶ sources | [ESA IoW 2025-12-18](https://www.cosmos.esa.int/web/gaia/iow_20251218) |
| Content detail | ESA states further specification (source counts per product) was expected end of June 2026; the content page is still partly a placeholder | ESA DR4 pages |

### 2.1 ESA-vs-red-team discrepancies found

- **Agreements (verified, not merely repeated):** the release date; the ~2.8 ×
  10⁹ / ~2 × 10⁹ / ~400 TB figures; the model-consolidation sentence *verbatim*;
  the absence of RUWE from the content page; the ten NSS tables; epoch
  astrometry; and the four-source DR3 NSS false-positive class. F7, F8, F9 and
  F15's date premise all survive external check.
- **Addition the red team missed, and it matters.** ESA publishes the four
  **`all_source_*` tables** alongside `gaia_source`. `gaia_source` in DR4 is
  **not** the catalogue — it is a **quality-selected 2.0 × 10⁹ subset of 2.8 ×
  10⁹**. Which of the two is our parent sample is therefore a *day-one
  pre-commitment*, not a detail, because the quality selection is correlated
  with binarity — the single nuisance the whole Paper-1 result is exposed to
  (7J). See §6.1. `all_source_flags` is also the most likely home of the
  applied-model flag that F7 needs.
- **Numbers in circulation that do not match the content page.** Secondary ESA
  text states "all ~2.5 billion sources will be published" with no quality
  filtering. The content page's own figures (2.8 × 10⁹ processed, 2.0 × 10⁹ in
  `gaia_source`) are taken as authoritative here; the 2.5 × 10⁹ figure is
  flagged as unreconciled and is not used in any yield estimate.
- **Measured, not quoted.** DR3 `nss_two_body_orbit` has **443,205** rows and
  `nss_acceleration_astro` **338,215** rows (measured this session, Q7/Q8 —
  see §4.2). Figures near 800k that circulate for "DR3 NSS" count all NSS
  solution products together, not this table.
- **Not verifiable today.** Whether DR4 publishes RUWE or a named analogue,
  whether the applied-model flag is exposed in `gaia_source` or only in
  `all_source_flags`, and the DR4 proper-motion improvement factor are all
  **unknown until the DR4 data model is published**. Every use of them below is
  conditional and labelled.

---

## 3. Timeline, and the F15 collision

The circulation priority is standing and is the author's call, not this
document's: **P1/P2 to arXiv before 2026-12-02 OUTRANKS the DR4 queue.**

| Date | Item | Owner | Blocking? |
|---|---|---|---|
| — | GitHub repo visibility flip (`thermal-horizon-rar` still answers 404 anonymously; all three papers cite its URL) | author (a click) | **YES** — blocks circulation, which outranks DR4 |
| when published | Read the DR4 **data model**, **known issues**, astrometric-solution description, NSS validation papers; pre-commit the error convention and the strata definition *from the documentation, before any kinematics* | main session | YES for every meter |
| 2026-11-01 | **DATE GATE (PROPOSED).** If `papers/note_pair_errors.md` and P1/P2 are not public by this date, the DR4 day-one queue reorders behind them. | author | — |
| 2026-12-02 | DR4 release | ESA | — |
| release + 0 | Documentation read, manifest audit, loader identity regression (§8.4). **No sky likelihood.** | main session | YES |
| release + days | Test 0 (§4), then the Gaia-native queue (§5) | main session | — |

**Rollback rule (PROPOSED, absent from DR4.md).** If DR4 slips, or ships
partially, or ESA amends the release post-launch: the manifest six-gate audit
(stage9h) re-runs, every affected stage re-runs, and no verdict letter written
against the pre-amendment release is quoted without the re-run. A partial
release is treated as a *different dataset*, not a delayed one.

**Compute reality (F18).** The 7J-z cube batches ran ~4 h GPU each; the arm
suite ~91 min. A sample larger by a factor N scales those roughly linearly in
the likelihood evaluation. A first-week claim on Test 3 (§5.3) is a **compute
claim** and must carry its own wall-clock forecast before it is scheduled.

---

## 4. TEST 0 (NEW, FIRST) — the NSS companion census

*Red-team F8: this is the collapse of the measurement's largest degeneracy, and
it is available day one from published tables with no modelling. It is not in
DR4.md at all.*

### 4.1 The question it answers

The entire binary chapter's dominant systematic is the **undetected companion
population**. Everything the program knows about it is a *photometric*
inference from overluminosity: the 7J→7J-z→7J-z2c chain, the retracted Part-A
scale (correction #18), the within-pair common-mode ρ ≈ +0.47, the v2c host
rate **0.29–0.32 per component**, the twin-heavy q-table (twin-convention host
≈ 0.23), and the *kinematic* preference of **~0.10–0.15**. The fork this feeds
is the operative one: **α_marg = 0.74/0.70 at +23.8/+23.2** versus **α_marg =
0.00** in the forced-multiplicity world, with the fifth-move exposure alive
(7J-z7: MATERIAL 4/4) and the exposure contained only by the census pair
(band = 9, cliff = 2; 7J-z8 rejects the twin-forced world at 5.5 × 10⁻¹⁴).

DR4 replaces that inference with a measurement over the period range that
dominates the wobble budget.

### 4.2 Tables and columns

**Primary (per component of every pair):**

| Table | Columns |
|---|---|
| `gaiadr4.nss_two_body_orbit` | `source_id`, `nss_solution_type`, `period`, `period_error`, `eccentricity`, `eccentricity_error`, `significance`, `mass_ratio`, `mass_ratio_error`, `a_thiele_innes`(+`b,f,g,c,h` and errors), `goodness_of_fit`, `astrometric_n_good_obs_al`, `flags`, `bit_index` — **all 77 DR3 columns enumerated in `data/dr4_dryrun_adql.txt` (probe Q9)**; the DR4 column set must be re-verified against the DR4 data model |
| `gaiadr4.nss_acceleration_astro` | `source_id`, acceleration terms and significances — the long-period companions that never close an orbit inside the baseline |
| `gaiadr4.nss_resolved_pair`, `nss_multiplicity`, `nss_masses`, `optical_pair` | disambiguation: which "companions" are resolved, which are catalogue artifacts |
| `gaiadr4.gaia_source` | `astrometric_excess_noise`, `astrometric_excess_noise_sig`, `ipd_frac_multi_peak`, `ipd_gof_harmonic_amplitude`, `phot_bp_rp_excess_factor`, `phot_g_mean_mag` |
| `gaiadr4.all_source_flags` | the applied-astrometric-model flag (see §6.1) and any RUWE analogue |
| `gaiadr4.epoch_astrometry` | per-epoch along-scan residual summaries, **source_id-list-driven server-side only** (62 TB; never bulk) |

**DR3 measurements banked as the regression target** (this session, `data/dr4_dryrun_adql.txt`):
`nss_two_body_orbit` = **443,205** rows; `nss_acceleration_astro` = **338,215**
rows; in a 28.3 deg² test patch at ϖ > 5 mas, **5** of 1,064 sources carry an
NSS orbit (≈ 0.5%). The DR4 rate must exceed this by a large factor for the
census to bite; **if it does not, that is itself the day-one result** and must
be reported as such rather than absorbed.

### 4.3 ADQL

Probe **Q10** in `calcs/dr4_dryrun_adql.py` is the cross-match grammar,
dry-run and returning. The day-one form joins the *pair table* (built by Q4)
to the NSS tables on `source_id`, once per component, producing a per-pair
4-state flag (neither / primary only / secondary only / both).

### 4.4 Statistic

1. The **measured companion rate per component**, `f_NSS`, in the wobble-relevant
   (q, P) window, with its **detection efficiency surface** ε(q, P, G, ϖ) — the
   rate alone is not the object; the program's round-10 lesson is that scalar
   conversions between differently-weighted moments are a trap, so the
   **(q, P) resolution function is passed, never a scalar**.
2. The **7J-z marginal re-run conditioned on the flag**, not marginalized over a
   prior: `f_comp` becomes a conditioning variable on the flagged stratum and a
   (much smaller) marginalized nuisance on the unflagged one.
3. A **companion-free stratum**: pairs with clean epoch residuals and no NSS
   solution in either component — the sample the α measurement has never had.
4. The standing cross-check: the census pair **(band = 9, cliff = 2)** recomputed
   in the companion-free stratum.

### 4.5 PROPOSED bars — *require pre-registration commit before DR4*

> **PROPOSED (not registered).** Report the completeness-corrected per-component
> host rate as an **interval, not a verdict word** (F8's explicit instruction).
> Pre-registered letter clauses to be locked in advance, e.g.:
> - the 90% interval lies **below 0.20** ⇒ the kinematic anchoring (~0.10–0.15)
>   is corroborated and the photometric 0.29–0.32 is the over-attributed one;
> - the 90% interval lies **above 0.26** ⇒ the photometric anchoring is
>   corroborated and the fifth move fires on its own terms;
> - otherwise **UNRESOLVED-CARRIED**, and the α fork stays open.
>
> **PROPOSED power gate:** the census may not be read unless the detection
> efficiency ε over the wobble-relevant (q, P) window exceeds a pre-stated floor
> on injected mocks; below it the stage STOPs with the firewall intact (§8.2).
>
> **PROPOSED trap-#11 clause:** each gate names, at pre-reg time, which clause
> of which letter its failure vetoes. A gate that vetoes nothing is labelled a
> diagnostic and cannot carry the letter (trap #24).

### 4.6 Mandatory false-positive calibration (F9)

The companion flag becomes a **conditioning variable**, so an uncalibrated
false-positive mode biases α in an unknown direction. Required before any
conditioning:

- read the DR4 **known-issues** page and the NSS validation papers on day one,
  before any cross-match;
- calibrate the flag's **false-positive rate and completeness** against the
  **redundant validator set** — the resolved-companion population, the
  photometric δ-channel already built at 7J-z2c, and the acceleration solutions
  — never a single hand-built cross-check (**trap #24**);
- carry the four named DR3 false positives (§2) as the known template of this
  failure mode;
- **standing clause:** any ESA post-release amendment triggers a stage9h
  manifest re-audit and a re-run.

---

## 5. The Gaia-native queue

Order rationale unchanged in spirit: instrument-level first, then paper-level,
then theory-facing. But **Test 0 precedes all of them**, and nothing runs before
the §8 guards pass.

### 5.1 Test 1 — P9, the width-object discriminator (instrument-level)

**Registered form:** PREDICTIONS.md §C P9 (immutable). Quoted, not altered.

**Staleness restatement required (F12).** P9 was registered against a
*monolithic* sq = 0.2 and a chase at P(f_pm = 3.0) = 0.54/0.97. Since
registration: 8Z/9A/9D/9E/9J dissolved the monolith into the **RUWE dose curve**
(Q1 ≈ 0 width, Q2/Q3 ≈ 0.1), and **9L landed NOISE-REAL** — pair velocity errors
≈ 2× formal at 0.2–2 kAU, **quality-independent**, which partly pre-decides the
error branch in-sample. The two branches must therefore be restated in
**per-stratum dose-curve language**, with the DR4 strata **re-derived** (§8.1),
and the decision tree phrased against the **two-sided upper limit** headline
(α ≥ 0.5 excluded on clean strata; α ≲ 0.3–0.5 allowed; α = 0 not excluded), not
the retired 0.68–0.74 operative band.

**Tables/columns.** `gaia_source` (or `all_source_astrometry`, per §6.1):
`source_id, ra, dec, parallax, parallax_error, parallax_over_error, pmra, pmdec,
pmra_error, pmdec_error, pmra_pmdec_corr, phot_g_mean_mag`, the RUWE analogue,
`astrometric_params_solved` **and the applied-model flag**;
`radial_velocity`, `radial_velocity_error` for the perspective correction (4Q).

**ADQL.** Probes Q1/Q4 of `calcs/dr4_dryrun_adql.py` (dry-run, returning);
the `_pair_body` selection is the pair construction, the census probe is Q1.

**Statistic.** Per-stratum posterior of `f_pm` and per-stratum width, read in the
boost-free 0.2–2 kAU bin, **within a fixed solution type** (§8.1).

> **PROPOSED bars — require pre-registration commit before DR4.**
> - Error-tail branch: the f_pm posterior expectation falls to ≤ 1.4 × formal
>   **within every solution-type stratum** *and* the width term's posterior mass
>   below 0.05 exceeds a pre-stated fraction.
> - Astrophysical branch: the width term persists at ≥ 0.15 in at least the
>   cleanest stratum while f_pm normalizes to ≤ 1.4 ×.
> - **Artifact branch (new, F7):** a quench that appears only in the *mixed*
>   sample and not within a fixed solution type is an **artifact**, not the
>   error-tail branch, and vetoes the letter.
> - Every number quoted at **bootstrap grade** (9U/9V realization wall: ~31
>   nominal −2lnL units = +1.8σ), with between-draw estimator scatter measured
>   by control variates first (**trap #22**).
>
> **PROPOSED trap-#26 requirement (F11) — blocking, before any letter:**
> score **known-outcome worlds** through the P9 discriminator under
> **DR4-scaled errors**: the 8F-c `fc10` error-tail world and a matched
> `sq`-true astrophysical world, formal errors scaled by the forecast DR4
> factor. **If the two branches do not separate at bootstrap grade, P9 needs a
> different statistic before the data land, not after.** The specific worry is
> concrete: an f_pm sitting at "≈ 2× formal" against a formal budget that has
> itself shrunk by ~3× is not obviously the same object, and both branches may
> drift toward 1. Nobody has checked. (ROUND 39 / 10R: a statistic that scored
> AQUAL μ10/μ20 as SAFE when they fail at 3 × 10⁻¹⁶.)

### 5.2 Test 2 — the pair-error bound re-measurement

Same sample and columns as Test 1, restricted to 0.2–2 kAU. The prior
measurement is the RNAAS note's **2.1–2.3 × formal**, envelope floor **1.77**,
RUWE-quartile-flat. That is a *prior measurement*, **not a DR4 bar**.

> **PROPOSED bar:** the DR4 re-measurement is compared to the envelope
> [1.77, 2.52] as a pre-stated consistency interval per stratum; a value below
> 1.4 in every stratum is the error-tail branch of Test 1 and must agree with it
> — a disagreement between Test 1 and Test 2 vetoes both letters.

**F15 dependency:** the note should be public before DR4 lands, or the program
re-measures its own unpublished result.

### 5.3 Test 3 — the operative upper limit

The 7J-z marginal machinery on the DR4 catalogue, with the companion prior
**replaced by Test 0's conditioning variable**. This is the row Test 0 exists to
serve.

> **PROPOSED bar:** the α upper limit quoted as the **anchor curve**, primary,
> per correction #17 — labels are annotations. Bootstrap grade mandatory.
> **PROPOSED STOP:** not scheduled until its wall-clock forecast is written
> (F18); it is a 4 h-GPU-class batch per cube at EDR3 size.

### 5.4 Test 4 — the perpendicular ceiling census

**F13 correction adopted: "~10× effective clean sample" is unsourced and is not
carried into this draft**, and "sharpens if real, dissolves if noise" is **not a
statistic** — a noise-driven band also grows with N.

The operative object is the self-defending **pair (band = 9, cliff = 2)**; the
raw count was retired ABSORBER-CONDITIONAL at 7K-b.

> **PROPOSED bar:** the **band/cliff ratio**, with the leakage null recomputed
> under DR4 errors and the pre-committed counts under each hypothesis stated
> before the read. The forecast yield must be computed from the EDR3 σ_v
> distribution under a stated DR4 improvement factor and printed **with its risk
> axis**, not asserted.

**Measured yield input available now.** The dry-run's full-sky pair-volume
estimate (§7) is the parent-population anchor for that forecast.

### 5.5 Test 5 — P7, the MI eccentricity signature (**rewritten, F6**)

**DR4.md's "eccentricity-resolved orbits" is wrong and is deleted.** Wide
binaries at 0.2–50 kAU have periods of 10⁴–10⁷ yr; DR4 NSS orbital solutions
cover periods within roughly the mission baseline. The program's eccentricity
information has always been **statistical**, through the v–r angle γ (3L/3N/4J)
and the radial-excess weight w_rad. DR4 improves the **precision of that 2D
statistic**, not the resolution of individual orbits.

Correct form: *improved γ and ṽ precision on a larger clean sample → re-run the
2D joint law × e machinery.*

> **PROPOSED bar:** the smallest α(s) step distinguishable at **bootstrap
> grade**, forecast from DR4-scaled errors against the **10O baseline of 35/200
> recovery of a −0.15 step** at Hwang-table grade. If the forecast does not beat
> that baseline by a pre-stated margin, the test is POWER-LIMITED and STOPs
> before the sky is read.

**One run, two registered rows:** the same instrument fires **P15 leg (i)**
(frozen-at-formation vs actively-pumped). Its registered kill (a) — e-distributions
shown actively pumped in a depth-tracking way — is a DR4-era clause.

### 5.6 Test 6 — P12, stream coherence (**added, F1**)

Registered 2026-08-11, Gaia-native, **absent from DR4.md**, and it is a
registered kill of the MI door's galaxy-scale leg. Data: Gaia stream members
(Pal 5, GD-1, Orphan) — astrometry plus RVs; the MG contrast class is Thomas+18
(arXiv:1709.01934, scout-grade, primary read owed).

**No instrument exists.** This row needs one built before DR4, or it will not be
a day-one item.

> **PROPOSED:** build the coherent-vs-local-field discriminator pre-release and
> score **known-outcome worlds** through it (trap #26) before any letter. No
> numeric bar is proposed here; the instrument does not exist yet.

### 5.7 Test 7 — P8 / P16, the ladder rungs

Listed for the DR4+/LSST era at population grade. 6L already measured that the
deep arm of 153 galaxies cannot read the c₃⁺ rungs; the binary rungs are P8's
DR4+ note. **Carried as a listed row with no day-one work**, so it is not
silently dropped again (F1).

### 5.8 P10 — the ambient-l2 amplitude (**added as an explicitly NON-DISCOVERY row, F17**)

P10 is the only registered row carrying an explicit numeric DR4 threshold, and
it was missing. It is added **with its bad news attached**, per F17:

- **σ_S = 0.02, the registered coarse discovery channel, is DEAD** at the
  mechanism's own derived amplitude. 10Q/ROUND 36 derived **e_a = 0.234 central,
  envelope [0.159, 0.388]** (matched-solve co-quote 0.268); S_max at the
  hi-envelope is **0.011 < 0.02** everywhere. **Revival requires e_a ≥ 0.662 =
  ×1.7 the hi-envelope.**
- The surviving channel is **σ_S ≈ 0.005**, and it is **edge**: central 0.0056
  — which ROUND 37 disclosed is the **convention-family maximum**, nominal
  convention ~0.002 — with a frozen-draw average of 0.0039.
- **This is not a κ = 1 crossing test.** 10N/ROUND 33: at detectable amplitude
  the response is a κ-**plateau** (S(1)/S(0.6) = 1.09); the detectable and
  diagnostic regimes are **mutually exclusive**.
- 10O/ROUND 34: any P10-amplitude cloud must be **frozen, τ_c ≳ 2.5 Gyr** — the
  derived cloud satisfies this automatically (×33).

**Coverage check owed before any effort is budgeted.** P10's instrument is an
**ambient-tide-correlated (R_gal / Z-height) split** in the transition window
x ∈ [0.5, 1.2]. Within the parallax-limited wide-binary volume there is very
little ambient contrast to split on — the same geometry that kills P4 (§6.2).
Scoping arithmetic (§6.2) gives ≲ ±10% in the radial term over the reachable
volume, with the vertical term adding, never subtracting.

> **PROPOSED:** run the P10 split-axis coverage check as a short arithmetic
> stage **now**, exactly as F17 asks. If the reachable ambient contrast cannot
> produce a split at σ_S ≈ 0.005, say so **in the plan**, not in December.
> No bar is proposed here until that check runs.

---

## 6. Reclassified rows

### 6.1 The parent-sample pre-commitment (blocking, all tests)

DR4's `gaia_source` is a **quality-selected 2.0 × 10⁹ subset** of 2.8 × 10⁹
processed sources; the complete astrometry is in `all_source_astrometry`. The
quality selection is **correlated with binarity**. This must be pre-committed
before any kinematics, with the choice and its selection bias stated:

- **(a) `gaia_source` only** — cleaner, but binarity-selected; the bias must be
  named in every downstream claim; **or**
- **(b) `all_source_astrometry`, stratified** — complete, but heterogeneous.

*(This item is not in the red-team review; it follows from the ESA content page.)*

### 6.2 P4 — the weak-ambient pair: **NOT-GAIA-FEASIBLE** (F3)

P4 is registered (PREDICTIONS.md §C) and stays registered, unaltered. What
changes here is its **data requirement**, written down with numbers, as 10I/10K
did for P2.

The scoping arithmetic (inputs stated; this is **not** a stage and is not
quotable until one runs):

| Quantity | Value | Source of the number |
|---|---|---|
| P4's requirement | binaries at ambient **e_N ≈ 0.4 a₀** (Δp ≈ 0.025) | PREDICTIONS.md P4 |
| Pipeline precision cut | σ_v < 0.03 km s⁻¹ | Paper 1 Table 1 |
| ⇒ EDR3 parallax floor | ϖ > 5 mas (**200 pc**) | Paper 1 Table 1; σ_v = 4.74 σ_μ/ϖ |
| DR4 baseline gain | 66 months vs EDR3's 34; σ_μ ∝ T^(−3/2) ⇒ **×2.7**, more with calibration gains — call it **×2.7–4** | ESA DR4 span; **conditional, unverifiable until the data model** |
| ⇒ DR4 reach at fixed σ_v | **540–800 pc**, ϖ > ~1.25–1.85 mas | linear in the pm gain |
| Solar-radius ambient | g ≈ V_c²/R₀ = (233 km/s)²/8.12 kpc = **2.17 × 10⁻¹⁰ m s⁻²** total; Newtonian g_ext = **1.2 a₀** after RAR inversion | 3S/3T |
| Radial contrast over ±0.8 kpc | flat curve ⇒ g ∝ 1/R ⇒ **∓9.8%**, i.e. e_N ∈ ~[1.09, 1.33] a₀ | arithmetic |
| Vertical term at \|z\| = 0.5 kpc | K_z = 2πGΣ ≈ **4.4 × 10⁻¹¹ m s⁻² = 0.36 a₀**, and it **adds** | Σ(\|z\|<0.5) ≈ 50 M⊙ pc⁻² |
| ⇒ off-plane pairs | sit at **higher**, never lower, ambient | arithmetic |
| To reach e_N = 0.4 a₀ | requires **R_gal ≳ 24 kpc** (factor-3 drop on a flat curve) or the halo | arithmetic |

**Conclusion (F3, confirmed by this arithmetic).** The ambient axis P4 needs
**does not exist in the volume where P4's own precision cut can be met**, in DR4
or DR5. The achievable lever arm is ~±10%; P4 needs a factor ~3.

This is exactly the **9V standing design check** — *verify y-coverage overlaps
the discriminating window before designing any tail/p/r instrument* — which
DR4.md skipped.

> **PROPOSED disposition:** either (a) **retire P4 from the day-one list with
> the number attached**, or (b) redesign it as a **low-precision, large-volume**
> statistic with a different likelihood, a different bar, and a new
> pre-registration — admitting the precision trade explicitly. Annotate P4 in
> PREDICTIONS.md with its measured data requirement (annotation only; the
> registered form is immutable). **This draft recommends (a).**

### 6.3 Moved out: DR4-era, non-Gaia (F2)

PREDICTIONS.md uses "DR4-era" as a **chronological** marker meaning "future,
better data of whatever kind". DR4.md read it as *the Gaia DR4 dataset*. These
rows are not Gaia tests and belong in a separate `DATA-ERA.md` (or a non-Gaia
section):

| Row | Why it is not a Gaia test | The actual data requirement |
|---|---|---|
| **P2** — environmental ordering (**strike-bearing**) | a *galaxy* tail-index ordering p(e_N); Gaia DR4 publishes no external-galaxy rotation curves | open arm: **density-field-selected (true void) galaxies with y ≳ 1 kinematic coverage** (10I); closed arm: **low-g_bar kinematics in genuinely dense environments** (10K) |
| **P6** — the γ discriminator | the 4Y/7C anchored subsample is Cepheid/TRGB-anchored **SPARC galaxies at Mpc distances**; Gaia measures Milky Way parallaxes | anchored-distance external-galaxy RCs with explicit non-circular modelling; IFU |
| **P1** void annotations | void asymptote / r-meter needs void-galaxy kinematics | density-field-selected void RCs with y ≳ 1 (9V feasibility close) |

**Why this matters more than tidiness (F2's own point).** P2 is the
**strike-bearing** row: DR4.md books "failure = REAL strike, pre-stated" against
a dataset that can produce neither outcome. The realistic failure mode is worse
than doing nothing — an underpowered environment split on whatever DR4-adjacent
sample is at hand returns a null, and a real strike is collected from an
instrument with no power. That is the 10I/10K pattern the program spent two
stages learning to avoid.

### 6.4 What does NOT get run (F18)

Stated explicitly, because a plan that lists only things to do is half a plan:

- **No** eccentricity-resolved wide-binary orbit fitting (does not exist, §5.5).
- **No** environment-ordering test on a Gaia sample (§6.3).
- **No** γ-discriminator work on Gaia data (§6.3).
- **No** P4 weak-ambient tail read at σ_v < 0.03 km s⁻¹ (§6.2).
- **No** σ_S = 0.02 P10 discovery claim (§5.8 — the channel is dead).
- **No** bulk `epoch_astrometry` pull (62 TB; source_id-list-driven only).
- **No** transfer of the 9J RUWE dose curve to DR4 — it is re-derived or not used.

---

## 7. Data access (F10)

**No community DR4 wide-binary catalogue will exist on day one.** EDR3 was
released 2020-12-03; the El-Badry/Rix/Heintz catalogue the entire pipeline
consumes appeared months later. The day-one state is 2.8 × 10⁹ sources, ~400 TB
of products, a loaded archive and per-query row limits.

**Dry-run executed, 2026-08-12** (`calcs/dr4_dryrun_adql.py`,
`data/dr4_dryrun_adql.txt`): **11/11 probes returned** against
`https://gea.esac.esa.int/tap-server/tap` in 9.5 s wall clock. The ADQL grammar
works: a self-join pair construction with parallax consistency, a
separation window mapped through the parallax, and an El-Badry-style proper-motion
ceiling.

**Full-sky pair-volume estimate.**

| Patch | Position | Pairs | Density | Full-sky |
|---|---|---|---|---|
| A | (150.0, +10.0), r = 3° (28.3 deg²), l ≈ 228, b ≈ +47 | 54 | 1.91 deg⁻² | 78,787 |
| B | (100.0, +10.0), r = 3° (28.3 deg²), l ≈ 204, b ≈ +12 | 57 | 2.02 deg⁻² | 83,164 |

**Bracket: 7.9 × 10⁴ – 8.3 × 10⁴ pairs** at ϖ > 5 mas, s ∈ [0.2, 50] kAU, on DR3
astrometry, **before** the main-sequence and chance-alignment cuts.

**Cross-check, and it lands.** Paper 1 Table 1 reaches **74,502 pairs** at the
same parallax and separation window *after* the main-sequence cut, from the
El-Badry catalogue. Our independent ADQL, which does **not** apply that cut,
returns 79k–83k — larger by the right amount in the right direction. The
grammar reproduces the parent-population scale.

**Risk axis shipped with the number** (standing rule):

- two cones only; the true distribution is latitude- and crowding-dependent, so
  this is a bracket of two samples, **not a confidence interval**;
- no chance-alignment subtraction — the 4I `R_chance` ladder is not applied;
- edge truncation: pairs straddling the cone boundary are lost (perimeter
  annulus ≈ 4.6% of area, partial loss within it ⇒ deficit ≲ 2%);
- the DR4 yield scales with the reachable volume (~f³ in the pm-improvement
  factor f) **only if** the magnitude and main-sequence cuts do not bind. **F13
  stands: compute it, do not assume "~10×".**

**Fallback order (PROPOSED, pre-committed here):**
(a) our own ADQL against the DR4 archive; (b) a community DR4 wide-binary
catalogue **if** one exists **and** it reproduces our cut ladder within a
pre-stated tolerance; (c) bulk mirror download (AIP `gaia.aip.de`,
ARI/Heidelberg, VizieR/CDS) only if (a) is rate-limited past week one.

**Storage/quota (PROPOSED):** the pair table and its ~20 pipeline columns are
the only bulk product; everything else is server-side. Products go under
`data/` with a `.gitignore` entry and a SHA256-pinned manifest audited by
stage9h, exactly as the 9H chain does today.

---

## 8. The guards

### 8.1 Solution-heterogeneity stratification — MANDATORY, day one (F7)

The 9L/9F pipeline reads `pmra_error1/2`, `pmdec_error1/2` as a homogeneous
formal budget on a **single-star fit**, and the entire quality-stratification arc
(8W/8Z/9A/9D/9E, the 9J dose curve) strata on **RUWE, a single-star fit
statistic**. In DR4, `gaia_source` parameters are consolidated across processing
modules and **a binary-star model may have been applied** (§2). RUWE is not on
the content page.

Blindly re-running the f_pm meter on DR4 measures the **mixture of solution
types**, not the noise — and the mixture is **correlated with binarity**, the one
nuisance the whole Paper-1 result is exposed to.

**Requirements (blocking, before any f_pm or RUWE-based meter):**

1. **Solution type / applied astrometric model is a first-class stratum axis**,
   alongside any RUWE analogue.
2. **Pre-commit** whether the primary sample is single-star-solution-only
   (clean but binarity-selected — flag the bias explicitly) or
   all-solutions-stratified. Ties to §6.1.
3. The **9J dose curve is re-derived on DR4**, never transferred.
4. **Stated in advance:** an f_pm quench that appears only in the *mixed* sample
   and not within a fixed solution type is an **artifact**, not the P9 error-tail
   branch (§5.1).

**DR3 template, measured this session** (probes Q2/Q3, 28.3 deg², ϖ > 5 mas,
n = 1,290):

| `astrometric_params_solved` | n | RUWE mean | RUWE min | RUWE max |
|---|---|---|---|---|
| 31 (5-parameter) | 970 (75.2%) | 1.397 | 0.750 | 26.06 |
| 95 (6-parameter) | 320 (24.8%) | 1.365 | 0.831 | 14.74 |

RUWE-null count at ϖ > 5 mas: **0** (2-parameter solutions carry no parallax, so
they never enter this sample — in DR4 the analogous exclusion is **not**
guaranteed and must be re-measured). **A quarter of the nearby sample already
sits in a different solution class in DR3**; DR4 adds the binary-model class on
top.

### 8.2 The sky-blind firewall — standing clause (F5)

Every recent stage that lacked power **never read the sky, by design**: 9R,
9R-b, 9R-c, 10I, 10K all halted at a pre-registered STOP with the firewall
intact. DR4.md contains no firewall clause and an explicit instruction to run
the meters first.

> **PROPOSED standing clause:** *No sky likelihood is evaluated on DR4 until
> (a) the loader identity regression passes (§8.4), (b) the manifest six-gate
> audit passes, and (c) each test's injection power gate passes on
> DR4-error-scaled mocks.* Each STOP letter is named in the test's own
> pre-registration, in advance.

Day one is the highest-temptation moment the program will ever face, and the
moment with the most unknown systematics. A plan that does not name the firewall
will not have one at 3 a.m. on 2026-12-02.

### 8.3 Trap #26 — score the known-outcome worlds first

Applies to **every discriminator statistic** in this plan, and **blocking for
P9** (§5.1) and P12 (§5.6). ROUND 39's lesson: the 10R statistic scored AQUAL
μ10/μ20 as SAFE when they fail V-N-T at 3 × 10⁻¹⁶. A statistic that has not been
run against constructed truths under a DR4-grade error model cannot license a
letter.

### 8.4 Loader identity regression — blocking build item (F16)

Standing rule: *every new pipeline MODE gets a regression against the old mode
at its identity point* (the GB0w precedent, which caught a silent amplitude-law
fallthrough on first firing). A DR4 loader **plus** an error-model adapter
**plus** a new strata definition is **three new modes**, and DR4.md regresses
none of them.

> **PROPOSED blocking item:** feed the EDR3 catalogue through the DR4 code path
> and reproduce the cut ladder (**14,071 pairs**) and the anchor (**1.078,
> CI 1.052–1.103**) bit-exactly, plus an **lnL-grade** reader identity against a
> stored cube. No DR4 sky read until it passes.

### 8.5 PREDICTIONS.md as the single source of truth (F1)

DR4.md was a hand copy that went stale in five days: P10 registered 2026-08-10
and P11–P16 on 2026-08-11 were silently absent, including the only row with a
numeric DR4 threshold.

> **PROPOSED rule:** the day-one plan becomes a **generated view over
> PREDICTIONS.md** — one row per registered prediction with an explicit
> `data source = Gaia DR4 / other / not applicable` column — and **any new §C/§E
> registration adds its row in the same commit**.

---

## 9. Credence map SKELETON (F14)

**The program's standing rule is: no credence movement outside pre-signed maps.**
DR4.md restates the rule and supplies no map. Pre-signing is only credible
*before* the data; writing this in December 2026 is post-hoc by construction.

**Current values (context, from the record — not a proposal):** anomaly-real
**53**; bath-mechanism conditional **8**.

**This build session must not invent credence numbers.** The table below is the
structure only.

| Test | Outcome cell | anomaly-real (from 53) | mech-conditional (from 8) |
|---|---|---|---|
| **Test 0** companion census | rate interval **below 0.20** (kinematic anchoring corroborated) | *TO BE PRE-SIGNED BY THE MAIN SESSION* | *TO BE PRE-SIGNED* |
| | rate interval **above 0.26** (photometric anchoring corroborated; fifth move fires) | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | **UNRESOLVED-CARRIED** | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | census **POWER-FAIL** (efficiency floor not met, firewall holds) | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| **Test 1** P9 | **error-tail** branch | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | **astrophysical-width** branch | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | **artifact** branch (quench only in the mixed sample) | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | branches **do not separate** (trap-#26 pre-check fails) | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| **Test 2** pair-error bound | confirmed in [1.77, 2.52] | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | revised **below 1.4** | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | Test 1 / Test 2 **disagree** | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| **Test 3** upper limit | α ≥ 0.5 **still excluded**, α = 0 still not excluded | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | α **interior and resolved away from 0** | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | α = 0 **excluded from above** (limit collapses to null) | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| **Test 4** census | band/cliff **ratio preserved** at larger N | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | ratio **dissolves** (band scales, cliff fills) | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| **Test 5** P7 / P15(i) | e-statistics **frozen-at-formation** | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | **actively pumped**, depth-tracking (P15 kill (a)) | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | **POWER-LIMITED**, STOP before sky read | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| **Test 6** P12 streams | coherent component **present** | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | MG-EFE local-field shape, **no coherent component** (MI galaxy leg killed) | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | instrument **not built in time** | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| **P10** (§5.8) | σ_S ≈ 0.005 channel **runs and is null** | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | σ_S ≈ 0.005 channel **finds the excess** | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| | **coverage check fails**, channel never runs | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| **Any test** | firewall STOP fires; sky never read | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |
| **DR4 slips / partial** | rollback rule invoked (§3) | *TO BE PRE-SIGNED* | *TO BE PRE-SIGNED* |

**Note for the signer.** The cells where **nothing moves** are the ones that
protect the program, and they must be filled in explicitly rather than left
implicit.

---

## 10. Build list (revised)

Blocking items are marked ⛔ — no DR4 sky read until they pass.

- [x] **ADQL written and dry-run against DR3** (F10) — `calcs/dr4_dryrun_adql.py`,
      11/11 probes returning; DR3 row counts banked in `data/dr4_dryrun_adql.txt`
      as the day-one regression target.
- [ ] ⛔ **Loader identity regression** (F16, §8.4): EDR3 through the DR4 code
      path → 14,071 pairs and 1.078 (1.052–1.103) bit-exact, plus lnL-grade
      reader identity against a stored cube.
- [ ] ⛔ **Solution-type stratification** built and pre-committed (F7, §8.1).
- [ ] ⛔ **Firewall clause** adopted with named STOP letters (F5, §8.2).
- [ ] ⛔ **Trap-#26 known-outcome scoring of the P9 discriminator** under
      DR4-scaled errors (F11, §5.1). If the branches do not separate, P9 needs a
      new statistic **before** December.
- [ ] **Test 0 instrument**: NSS cross-match + efficiency surface + the
      conditioned 7J-z marginal (F8, §4), with the false-positive calibration
      (F9, §4.6).
- [ ] **Coverage checks**: P4 (F3, §6.2 — arithmetic drafted here, needs a stage
      to be quotable) and P10's ambient-split axis (F17, §5.8).
- [ ] **Pre-registrations** for each surviving test: letter grammar,
      gate → clause mapping (trap #11), injection power gate with a pre-stated
      STOP, between-draw estimator scatter measured by control variates
      (trap #22), bootstrap-grade quoting mandated in the letter (F4).
- [ ] **Credence map signed** (F14, §9) — main session only.
- [ ] **Documentation read** with an owner and a date: DR4 data model, known
      issues, astrometric solution description, NSS validation papers, release
      scenario page (F18, §3).
- [ ] **Error-model adapter**: convention pre-committed from the published data
      model, *before* any kinematics.
- [ ] **Census re-count script**: 4J/7K-b machinery, (band, cliff) convention
      frozen, leakage null recomputed under DR4 errors (F13, §5.4).
- [ ] **P12 instrument** built, or the row is honestly not a day-one item (§5.6).
- [ ] **DATA-ERA.md** split for the non-Gaia rows (F2, §6.3), with each row's
      named data requirement.
- [ ] **PREDICTIONS.md annotations** (annotations only — registered forms are
      immutable): P4's measured data requirement; P6's; P10's coverage note.
- [ ] **Repo visibility flip** (author's click) — blocks circulation, which
      outranks DR4 (§3).

---

## 11. Standing rules that survive into the DR4 era

Unchanged from DR4.md, restated: pre-reg before any run; bars locked before
data; the six-gate manifest audit after any refetch (stage9h); **no credence
movement outside pre-signed maps**; NOTES/LOG/LEDGER discipline unchanged. The
in-sample function freeze (PREDICTIONS §0) remains in force — DR4 is
out-of-sample territory, where function-level claims are allowed to move again.

Added by this revision: the firewall (§8.2), the trap-#26 pre-check (§8.3), the
identity regression (§8.4), the generated-view rule (§8.5), the rollback rule
(§3), and the explicit not-run list (§6.4).

---

## 12. Finding traceability

| Finding | Where addressed |
|---|---|
| F1 stale vs PREDICTIONS.md | §5.6 (P12), §5.7 (P8/P16), §5.8 (P10), §5.5 (P15 i), §8.5 (generated-view rule) |
| F2 non-Gaia rows | §6.3, §6.4 |
| F3 P4 infeasible | §6.2 (arithmetic attached) |
| F4 no numbers | every PROPOSED bar; §10 pre-registration item |
| F5 no firewall | §8.2 |
| F6 P7 misdescribed | §5.5 (rewritten) |
| F7 solution heterogeneity | §8.1 (with measured DR3 template), §6.1 |
| F8 NSS companion census | §4 (Test 0, first) |
| F9 NSS false positives | §4.6, §2 (four named DR3 sources) |
| F10 no data-access plan | §7, `calcs/dr4_dryrun_adql.py` |
| F11 trap #26 on P9 | §5.1 (blocking), §8.3 |
| F12 staleness of the width object | §5.1 (dose-curve restatement, upper-limit vocabulary) |
| F13 "~10×" unsourced | §5.4 (dropped; ratio statistic proposed), §7 (measured parent volume + risk axis) |
| F14 no credence map | §9 (skeleton, unsigned) |
| F15 timeline collision | §3 (date gate, repo dependency) |
| F16 no identity regression | §8.4 |
| F17 P10 missing + its bad news | §5.8 |
| F18 operational scaffolding | §3 (compute, rollback, documentation), §6.4 (not-run list), §7 (storage/quota) |

---

*Plan v2 ends. Nothing above is pre-registered or credence-signed; every bar
marked PROPOSED awaits its own bar-locking commit.*
