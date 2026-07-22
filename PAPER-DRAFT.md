# Testing the Thermal-Horizon Reading of the Radial Acceleration Relation, from Galaxies to Wide Binaries
## (working draft v2 — rewritten 2026-07-23 after Stages 3B–4D; do not circulate)

**Candidate titles:**
1. *The radial acceleration relation as a Planck-oscillator law: first coefficient
   tests and an independent wide-binary measurement*
2. *Gravity's low-acceleration expansion: the Bernoulli ladder, its zero-point
   term, and a strength-one wide-binary anomaly*

---

### Abstract (draft v2 — numbers final, prose to polish)

The interpolating function preferred by the radial acceleration relation (RAR),
ν = (1−e^(−√y))⁻¹, was shown by Cadoni & Tuveri (2019) to equal 1 + n_BE(x),
x = √(g_bar/a₀) — the Bose–Einstein occupation of thermal excitations of the
de Sitter horizon, with a₀ = H/2π derived. That framework has remained untested
beyond rotation-curve fitting. We present its first empirical stress-tests.
(1) We note two exact structural consequences, apparently unremarked: ν equals
½ + ½·coth(x/2) — the Planck oscillator mean-energy law, whose frozen (x → ∞)
limit ν → 1 identifies Newtonian gravity with the pure zero-point response —
and the low-acceleration expansion coefficients are the Bernoulli numbers,
ν = Σ B_n⁺ x^(n−1)/n!, a ladder of parameter-free predictions. (2) We test the
first nontrivial rung on 2,700 SPARC points: the NLO-coefficient-½ branch
(shared by the occupation law and the simple-ν family) is preferred over the
c₁ = 0 branch (standard-μ family) in 198–200 of 200 galaxy-level bootstrap
resamplings under a raw-χ² likelihood, deflating to a sign-robust strong lean
(Δ(−2lnL) = 56, 166/200) when the relation's intrinsic scatter is
marginalized — both treatments are reported; the truncated-expansion
estimator alone is power-limited (σ(c₁) ≈ 0.5), making the full-function
comparison the honest instrument. Extending to the weak-lensing RAR
(KiDS-1000, two public reductions, ~2 dex deeper) anchors a₀ and rejects the
Newtonian control on 15 points alone, but cannot reach rung 2: the resolving
power for 1/12-vs-1/8 is 0.1σ at the surveys' 0.2-dex stellar-mass
calibration — the ladder's next rung is measurable in principle, not yet in
practice. (3) We independently measure the low-acceleration anomaly in
14,071 Gaia EDR3 wide binaries with a joint 2D (velocity, direction-angle)
likelihood, a physically-identified nuisance set (a ~20% near-parabolic orbit
sub-population; companions bounded jointly by photometry and kinematics; two
contamination channels with distinct 2D signatures; the catalog's velocity-
consistency selection emulated — the fitted near-parabolic fraction
w_rad = 0.20 independently matching the e>0.9 mass fraction, 20–22%, of the
superthermal eccentricity law of Hwang, Ting & Zakamska 2022), and full
seed + bootstrap error budgets:
the boost strength is α = 1.18 ± 0.11 (simple-ν) / 1.13 ± 0.13 (occupation
law) relative to the galactic calibration at the RAR-inverted Newtonian
external field g_N,ext = 1.15 ± 0.05 a₀, with Newtonian dynamics excluded in
all 2000 bootstrap contests (min ΔlnL +53). Two further observational results
stand alone: the wide-binary velocity-direction distributions are U-shaped —
inconsistent with any single power-law eccentricity family — and an
orbit-population "realization systematic" larger than the Newton-vs-MOND gap,
which plausibly contributes to the published disagreement between wide-binary
analyses. The velocity-direction (v-r angle) distributions confirm and extend
— to 50 kAU and inside a joint gravity-law fit — the superthermal trend of
Hwang, Ting & Zakamska (2022), executing the gravity application those
authors proposed; a folded γ ≈ 82° perpendicular excess remains open and is
disclosed. The framework's remaining falsifiers are stated: the rung-2
coefficient (1/12 vs the simple family's 1/8 — shown here to require lensing
mass cross-calibration at the ~0.02-dex level, or Gaia DR4), and a₀ ∝ H(z).

### Claims map (post-retraction, Stage 4C/4D — the paper must respect this)

- **Cadoni & Tuveri 2019 (cite prominently; one-line independent-arrival note
  permitted)**: the identity ν = 1 + n_BE, its derivation from dS-horizon
  thermal excitations, a₀ = H/2π.
- **Hwang, Ting & Zakamska 2022 (cite prominently — Stage 4G, correction #9)**:
  the v-r angle (γ) method on this catalog; the superthermal eccentricity
  distribution at >10³ AU (α ≈ 1.2–1.3 to 31.6 kAU; e>0.9 enhanced); the
  PROPOSAL of v-r angles as a gravity test (with Banik & Zhao 2018, 2021).
  Our w_rad = 0.20 is their confirmation (20–22% implied e>0.9 fraction),
  not our discovery of the phenomenon.
- **Ours (scout-clean, final INSPIRE pass owed at submission)**: the coth /
  zero-point form and its reading (Newton = frozen vacuum response); the
  Bernoulli-ladder statement; the first coefficient test (4B) and its
  power-analysis (4A); the dual-likelihood branch re-test and the rung-2
  lensing power analysis with its honest null (4E); the bath matrix — the
  simple-ν-as-classical-bath identity and the ¼-branch test (4F, scout
  pending); the entire wide-binary program (Stages 2–3V); the EXECUTED joint
  2D (ṽ×γ) gravity-law × e-mixture fit (Hwang+22 proposed it; we ran it,
  fenced, to 50 kAU) and the w_rad ↔ superthermal-α cross-validation; the
  γ≈82° residual; the realization systematic; the screening-index p; the
  g_N,ext convention analysis (AQUAL-total vs QUMOND-Newtonian inputs).

### Section skeleton → assets (assemble from NOTES next session)

1. **Introduction** — anomaly framing; the Chae/Banik/P&S three-way split as
   motivation for a systematics-first re-measurement.
2. **The structure of the RAR** — C&T attribution; coth/zero-point; Bernoulli
   ladder (sympy-verified); the bath matrix as a table (simple-ν = exact
   classical self-consistent bath ⇒ BE-vs-simple = quantum-vs-classical;
   the ½'s two readings honestly stated); what is and is not claimed.
   [Stage 4C/4D/4F notes; INSPIRE pass owed on the simple-ν derivation]
3. **Screening index** — ν_p family; PRIMARY: M/L-marginalized
   p = 0.578 +0.121/−0.115 (4H; fixed-M/L 0.443 +0.063/−0.050 as the
   comparison row; p = ½ at 0.7σ; Cassini p > 0.234 passed at the 16th pct).
   [calcs/sparc_rar_fit.py, calcs/stage4h_p_ml.py]
4. **Coefficient tests, SPARC + lensing RAR** — 4A estimator + honest
   power-null; 4B branch verdict 198–200/200 (raw χ²); 4E scatter-marginalized
   re-test (−56, 166/200; deep window agnostic — quote both) + the
   lensing-anchored rung-2 null (resolving power 0.1σ; the 0.2-dex mass wall);
   a₀ = 1.206e−10 (raw) / 1.00e−10 (marginalized), f_ML ≈ 1.1 recovered free.
   [calcs/stage4a_nlo_test.py, calcs/stage4b_branchcomp.py,
   calcs/stage4e_lensing_rar.py + stage4e_diag.py]
5. **Wide binaries: pipeline and model** — data cuts; joint 2D (ṽ, γ)
   likelihood; nuisance identities (w_rad = 0.20 stable across 8 model
   variations; companion sector photometry+kinematics-bounded; contaminants;
   catalog selection; sigv/√2 noise convention). [Stages 2B–3O]
6. **Wide binaries: results** — α = 1.18 ± 0.11 / 1.13 ± 0.13 at
   g_N,ext = 1.15 ± 0.05a₀; Newton excluded in 2000/2000; α(g_ext) degeneracy
   mapped; the external-field convention analysis (our 2G bug owned; the
   AQUAL-total vs QUMOND-Newtonian distinction made explicit for the field).
  [Stages 3P–3V, 3S/3T]
7. **The eccentricity sector and the realization systematic** — the radial
   excess as CONFIRMATION + extension of Hwang, Ting & Zakamska 2022
   (w_rad = 0.20 vs their implied 20–22% e>0.9 fraction — lead with the
   cross-validation); circular vetoed; the γ≈82° residual disclosed as open;
   the realization systematic (still apparently ours). [Stages 3D, 3L–3N, 4G]
8. **Transparency appendix** — the nine logged corrections, verbatim from
   NOTES (axial 2/π; non-conservative EFE; a₀ conventions; grid-edge bullseye
   retraction; mass-error hypothesis; circular-mixture reading; identity
   priority; raw-χ² branch-verdict inflation + the retracted SPARC simple-lean,
   4E; U-shape priority → Hwang+22, 4G). This appendix is the paper's
   credibility spine — do not trim it.
9. **Predictions & program** — rung 2: tested and shown to need ~0.02-dex
   lensing mass cross-calibration or Gaia DR4 (4E — a kill test the authors
   aim at their own preferred reading, with the current power stated);
   a₀ ∝ H(z); DR4 forecasts.

### To do before any submission
- [ ] Assemble sections from NOTES (fresh session; NOTES Stages 1→4D are the source)
- [ ] Final INSPIRE pass on all priority claims (coth, Bernoulli, coefficient tests,
      U-shape, realization systematic)
- [ ] Verify Chae 2023's exact external-field sentence against full text
- [x] M/L marginalization on the p measurement (Stage 4H — p = 0.58 ± 0.12 primary)
- [x] Chance-alignment stress test (Stage 4I — stable, direction-check clean)
- [x] Rung-2 lensing-RAR test (Stage 4E) — landed: honest null with power
      quantified; folded into §2/§4/§9 and the abstract
- [ ] Colleague review (COLLEAGUE-BRIEF.md — update numbers to 3V finals first)
- [ ] Zenodo DOI snapshot; then arXiv (endorsement via colleague)
