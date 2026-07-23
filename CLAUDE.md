# PhysicsResearch — session handoff

Low-acceleration gravity research program run jointly by the user (hfili — sharp,
non-academic, wants first-principles reasoning + numbers, enjoys building) and Claude.
Read [NOTES-horizon-inertia.md](NOTES-horizon-inertia.md) (the lab notebook, chronological
with retractions), [TODO.md](TODO.md) (priority queue), [PAPER-DRAFT.md](PAPER-DRAFT.md).

## Environment
- Windows 11, PowerShell. Python is `py` (NOT `python`). Installed: numpy, scipy, sympy,
  astropy, cupy-cuda12x[ctk] — CUDA kernels compile and run on the RTX 5090 (fp64 ~1.6
  TFLOPS; 1e6-binary population ≈ 10–70 s depending on force mode).
- Large datasets are gitignored; re-fetch: SPARC (Zenodo 16284118 → data/sparc/), EDR3
  binaries (Zenodo 4435257 → data/edr3_binaries.fits.gz, 1.4 GB), El-Badry-Rix 2018
  (VizieR TAP → data/widebinaries_elbadry_rix2018.csv). Scripts in calcs/ document exact URLs.
- `private/` is gitignored context (project origin documents); do not commit or publish it.
- Agent budget: user is on Max 5x — spawn subagents on `model: "haiku"` (or their Ollama)
  for literature checks; do the physics/code in the main loop.

## Key results (all in NOTES with scripts)
- Screening index of the RAR transition: fixed-M/L p = 0.443 +0.063/−0.050
  ([calcs/sparc_rar_fit.py](calcs/sparc_rar_fit.py)); **M/L-marginalized (4H, primary):
  p = 0.578 +0.121/−0.115, f_d=1.22±0.10, a₀=(1.05±0.10)e−10 — p=½ sits 0.7σ inside**;
  Cassini independently requires p>0.234 (16th pct 0.462, comfortably passed). Chance
  alignments: boost stable under 20× R_chance tightening (4I).
- **Bose–Einstein identity** (exact): RAR ν = 1 + n_BE(x), x=√(g_N/a₀) — **PUBLISHED
  PRIOR ART: Cadoni & Tuveri 2019 (arXiv:1904.11835), Eq. 23 + derivation from thermal
  dS-horizon bosonic excitations, a₀=H/2π derived. Priority claim RETRACTED (Stage 4C;
  three scout sweeps misread the paper — primary sources only). Ours: the NLO=½ TEST
  (4B: ½-branch 198–200/200 over standard-μ — apparently the first empirical test of
  the C&T structure), NNLO discriminator (open, leans simple), a₀∝H(z) (future). Paper
  §2 = "C&T derived it; we test it." Map their citation chain before writing.**
- Wide-binary ṽ boost (our EDR3 pipeline, 14,071 pairs): **1.086 (CI 1.064–1.110)**;
  RUWE-stable; triples capped ≲5% by tail shape ([calcs/stage2c/2d](calcs/)).
- QUMOND EFE solver (validated, gates G1=0.01%): BE-EFE ~4% weaker than simple-ν, Wien tail
  ([calcs/qumond_efe_solver.py](calcs/qumond_efe_solver.py), tables data/efe_boost_*.npy).
- **Realization systematic** (likely novel, field-relevant): orbit-population realization
  shifts distribution lnL by MORE than the Newton-vs-MOND gap; no published quantification
  (scout-confirmed). [calcs/stage3a_likelihood.py](calcs/stage3a_likelihood.py).
- Hierarchical fits (Stages 3B→3F in NOTES): v3 contamination + η prior did NOT stop α
  corner-seeking. Stage 3D autopsy found why: data ṽ distribution is BROADER than every
  model in every s-bin (the model lacks a broadening convolution; data ṽ uses error-laden
  photometric masses). Stage 3E: multiplicative smear σ_m localizes at ≈0.2–0.25, α
  corner-seeking collapses (α=1 beats α=2), and the **honest Newton rejection deflates to
  ≈ +63/+66 lnL** (earlier +264…+338 was misfit-inflated; retraction-grade correction,
  logged). Median boost 1.086 is smear-immune (the anchor). **v4 fit (Stages 3F/3G/3H):
  α localizes INTERIOR for both laws, seed-robust 12/12** — MC over 6 realizations gives
  **simple α̂ = 0.93±0.11, BE α̂ = 1.30±0.15 (realization scatter, flat profile top over
  [0.75,1.5]), Newton ΔlnL = +55±4**; σ_m=0.30 chosen by all 12 fits; BE>simple α̂ offset
  every seed (consistent with BE-EFE being ~4% weaker); BE vs simple lnL: no shape
  verdict. Error budget DONE (3I bootstrap ⊕ 3H realization): **simple α = 0.98 ± 0.23,
  BE α = 1.21 ± 0.30; Newton loses in all 1000 bootstrap replicates (min +30)**. BUT
  Stage 3J MEASURED the mass error (MS width 0.275 mag → σ_m(mass)=0.024, 12× below the
  fitted 0.30 — mass-error hypothesis REFUTED) and found 12.3% of component stars
  overluminous = unresolved companions, whose broadening GROWS with s like the signal —
  but the **v5 physical-multiplicity fit (Stage 3K) exonerated it: α stays interior ≥1,
  Newton +108/+120, companion model fits −420 lnL worse than the smear, f_comp capped at
  0.1 by the ṽ data. Hidden-triples hypothesis disfavored by shape, amplitude, and α's
  survival.** v-angle (3L, mass-immune): γ distributions U-SHAPED — **PRIOR ART (4G,
  correction #9): Hwang, Ting & Zakamska 2022 (arXiv:2111.01789) published the
  superthermal e-distribution (α≈1.2–1.3 at 1–31.6 kAU, e>0.9 enhanced) via the same
  v-r angles on the same catalog, and proposed the gravity application (crediting Banik
  & Zhao 2018/2021); DIVIDEND: their law implies 20–22% at e>0.9 vs our fitted
  w_rad=0.20 — external validation; ours = executing the joint law×e fit, 31.6–50 kAU,
  the γ≈82° residual** (P&S 2025 = Newton-favored, ṽ-only, no directions). Joint 2D
  (ṽ×γ) fits: v6 (3M) vetoed the circular reading (w_circ=0, −174; circular NARROWS ṽ —
  3L corrected); **v6b (3N, [calcs/stage3n_v6b.py](calcs/stage3n_v6b.py)) IDENTIFIED the
  broadening: w_rad=0.20 interior both laws (+129/+141) — a ~20% near-parabolic (e>0.9)
  sub-population; all contaminant fences localize; f_pm minor; α interior at 1.5 (span
  across broadening models: α̂ 1.0–1.5 — report the span); Newton +93/+102; first
  BE-over-simple lead (8.7 lnL, near noise). v7 (3O, [calcs/stage3o_v7fit.py](calcs/stage3o_v7fit.py))
  closed the catalog-selection systematic (envelope measured from data: physical bound
  ceiling v·√s≈2.2, catalog cut ≈3.0; acceptance on model+templates): **α INVARIANT
  (interior 1.5, Newton +98/+105), w_rad strengthens, tail residuals gone, BE leads
  simple a 2nd time (7.5 lnL). v7 budget (3P, 6 seeds): simple α̂ 1.57±0.35 (2/6 flat-top edges), **BE α̂
  1.462±0.042 interior 6/6 — the physical model's best-conditioned number**; Newton
  +100±9/+104±5; w_rad=0.20 12/12; **BE-over-simple lead NOT significant (+3.3±2.2 SE) —
  no ν-family discrimination**. Error budget COMPLETE (3P+3Q): **FINAL v7 numbers:
  simple α = 1.52 ± 0.46, BE α = 1.54 ± 0.13; Newton loses in all 1000 bootstrap
  replicates (min +38/+59)**. g_ext arc RESOLVED (3S+3T): the old e_N=1.9 was a
  **cross-formulation convention bug from Stage 2G** (AQUAL-total value fed to our
  QUMOND solver; correct QUMOND input = NEWTONIAN g_N,ext = 1.15±0.05a₀ from RAR
  inversion — our own arithmetic). **At the physical field (tables at 1.2a₀) BOTH laws
  localize at the parameter-free α=1: α̂ = 1.17 (simple) / 1.11 (BE), interior, Newton
  +108.7/+98.8. This supersedes 3Q's α=1.54±0.13 (wrong-convention conditional).**
  w_rad=0.20 at all 8 model variations. **MEASUREMENT PROGRAM COMPLETE (3U/3V): FINAL
  α = 1.18 ± 0.11 (simple) / 1.13 ± 0.13 (BE) at g_N,ext=1.2a₀, interior 1000/1000 &
  998/1000, Newton excluded in all 2000 bootstrap contests (min +53); realization
  scatter collapsed to 0.045/0.039 at the physical field; simple out-fits BE on every
  seed (−12.3±2.2 — binaries lean simple; SPARC is AGNOSTIC per 4E, the earlier SPARC
  simple-lean was a raw-χ² artifact, retracted; BE credence ~20–25%). NLO kill test
  (4A/4B + 4E caveat): ½-branch beats standard-μ 198–200/200 under raw χ², deflating to
  a sign-robust strong lean (Δ−2lnL −56, 166/200) under scatter-marginalized likelihood
  — QUOTE BOTH. Rung-2 lensing test (4E, [calcs/stage4e_lensing_rar.py](calcs/stage4e_lensing_rar.py)):
  Mistele+24 + Brouwer+21 KiDS data fetched (data/lensing_rar/); HONEST NULL — resolving
  power 0.09–0.10σ (0.2-dex lensing mass systematic is the wall; needs 0.02 dex or DR4);
  within-branch BE−simple −19±15, carried by 3/153 galaxies = no verdict; Newton
  +2777/+1659 on 15 lensing points alone; correction #8 logged. Bath matrix (4F,
  [calcs/stage4f_bathmatrix.py](calcs/stage4f_bathmatrix.py)): **simple-ν IS the exact
  classical self-consistent thermal bath (ν=1+1/(νy) ⇒ ½+√(¼+1/y); scout/INSPIRE pass
  pending) — BE-vs-simple = quantum-vs-classical bath**; the new quantum-bootstrap
  ¼-branch (c₁=¼, c₂=7/96) tested: dead-grade raw (4–7/200), disfavored honest
  (+27/+9), sign-consistent; c₁ dose-response peaks at ½ (0 dead, ¼ disfavored, ½
  preferred, 1 Cassini-dead).
  γ≈82° RESOLVED (4J, [calcs/stage4j_gamma82.py](calcs/stage4j_gamma82.py)): the
  islands = peri/apo faces of the near-parabolic population + unmodeled
  closest-approach arm; **NEW: the perpendicular velocity CEILING — 11 wide pairs in
  the Newton-forbidden band [√2, 1.67), leakage null P=3.8e-9, cliff at the α=1
  boosted escape edge √(2·1.36)≈1.65 (P=0.62/0.91); flybys (v∞ absurdity), triples,
  chance, selection each excluded; N=11 caveat stated; GPU variant queued (#2e-b).**
  **Solar quadrupole (4K, [calcs/stage4k_quadrupole.py](calcs/stage4k_quadrupole.py)):
  the scale-free EFE solution ⇒ Q₂ = 3.9e-26·(α/1.15) s⁻² for BOTH families
  (transition-sourced, family-blind; Wien tail does NOT rescue BE) — exceeds the
  Cassini bound (Hees+14: (3±3)e-27) by ~4.3×: independent BINARY-calibrated
  reproduction of the Desmond–Hees–Famaey 2024 tension (arXiv:2401.04796), immune to
  their M/L+bulge escapes; solver cross-validated vs Blanchet–Novak μ₁ at 15% (G6,
  also confirms the 3T convention mapping). ESCAPES: modified inertia (→ TODO #18,
  MI-vs-MG on our 2D data = top theory question) or EFE-screened theories (#8).
  MOND-Planet-9 (2304.00576) runs on the capped term — disfavored. **MI-vs-MG
  EXECUTED + BUDGETED (4L/4L-b, [calcs/stage4l_mi_runner.py](calcs/stage4l_mi_runner.py),
  6 seeds — CORRECTION #10: seed-31's "MG wins all 8" was realization luck):
  budget verdict = mi_t (EFE on) TIES MG (−3.5±3.3/−0.8±2.5); no-EFE MI dead 12/12
  (−20…−28, α̂→0.5–0.6 sharp = data DEMAND the EFE amplitude); Newton dead 24/24
  (+71…+108); w_rad=0.2 in 22/24. **The Saturn-safe EFE-respecting MI branch fits
  the binaries AS WELL as MG — the 4K paradox has an open door (amplitude pinned,
  mechanism agnostic).** First-of-kind (provisional; Paci+2020 = rotation-curve
  prior). 4M ([calcs/stage4m_fly90.py](calcs/stage4m_fly90.py), #2e-b closed): pure
  90°-arm flyby template +4.7/+3.4 lnL, α̂ moved ≤0.06 — the 4J ceiling residual
  formally closed as model-shape, boost untouched. Credences: anomaly-real ~70%;
  conditional: EFE-respecting MI ~40%, MG-with-quadrupole-fix ~30%, unknown ~30%;
  BE ~15%.**
  Reconciliation EXECUTED (4N, [calcs/stage4n_banikstyle.py](calcs/stage4n_banikstyle.py),
  proxy-grade): **unfencing the companion fraction absorbs ~60% of the Newton deficit
  (needs fractions 3J forbids; Banik's free fit = 69%); +anchor-drop+ṽ-only → 2/3 of
  the significance manufactured away, α̂ biased to 0.7; the detection never flips
  (+30–38 under full freedom); the γ channel = the parameter-protector (vtonly
  α̂_BE→1.55). Residual leg = H&C's unconvolved sub-error binning (not ablatable in
  our always-convolved pipeline) + the 3A realization systematic.**
  **PAPER v1.1 VERIFIED ([PAPER.md](PAPER.md), Stage 4O 2026-07-22): author =
  Filip Hájek (independent researcher); nine-agent primary-source pass on every
  novelty claim; CORRECTION #11 = √2 ceiling is P&S-2018/B&Z-2018 founding prior
  art ("the key discriminant", P&S 2023) — ours is the perpendicular census + edge
  termination, §7.2 reframed, Appendix A now eleven corrections; phantom citation
  "Paci 2020" → Petersen & Lelli 2020 (A&A 636, A56); Chae numbers fixed (2023:
  1.43±0.06/10σ at arXiv:2305.04613; 2024: 1.37/5.8–9.2σ); H&C critique gains
  third author Aguayo-Ortiz; refs INSPIRE-verified + 8 missing added; C&T citation
  chain (31 papers) contains no empirical test — first-test claim SUPPORTED;
  simple-ν bath identity unclaimed anywhere (Famaey & Durakovic 2025 "no existing
  clear derivation" cited); Timeflow = Trofimov, CQG 43, 135020, distinct
  mechanism; Vokrouhlický+24 de-orphaned into §8.1; CITATION.cff + .zenodo.json
  added (title/author for the DOI record).**
  **4P/4Q (2026-07-22, PAPER → v1.2): the 2602.24035 gate read — byline is
  COOKSON et al. (lead = independent researcher; Banik/El-Badry/Sutherland/
  Penoyre/Pittordis/Clarke; "Desmond" was a scout hallucination, 3rd phantom
  byline caught) — 1,421-pair median-ṽ flatness test, ~1500×/2.7σ Newton.
  The read caught OUR gap: no spherical-projection (systemic-RV perspective)
  correction anywhere in the program → Stage 4Q audit, CORRECTION #12
  ([calcs/stage4q_perspective.py](calcs/stage4q_perspective.py)): effect
  PRESENT at predicted size (slope 0.92; Cookson fig 7 on our selection),
  too small ×5 (Newton+perspective → 1.016 of 1.086), immune ṽ_perp component
  boosts MORE (1.151 CI 1.115–1.197) = opposite of artifact; **corrected
  anchor 1.078 (CI 1.052–1.103)** — quote both; ceiling census ±1 at edges
  (cliff intact, null ≤1e-8); α exposure inside ±0.11. §7.4 written (checklist
  mapping; their e-constancy assumption vs the Hwang trend they cite —
  measured sign SUPPRESSES a step; √2-band objection answered by the
  perpendicular column). Twelve corrections in Appendix A. Repo renamed
  thermal-horizon-rar (user clicks the GitHub rename; files/remote updated);
  LICENSE = CC-BY-4.0 added.**
  **4S (2026-07-23): ĉ₁ MEASURED ([calcs/stage4s_c1fit.py](calcs/stage4s_c1fit.py)):
  continuous-λ family (c₁=λ/2, series-gated); marginalized joint SPARC+lensing:
  ĉ₁ = 0.450 (profile 0.385–0.519; galaxy-bootstrap 0.43 +0.29/−0.25), ½
  INTERIOR at 0.7σ, c₁=0 excluded (Δ−2lnL 56.3: 7.5σ profile / 95.5%
  bootstrap — the gap = population variance, bootstrap primary), ¼ disfavored
  3.1σ; raw-χ² EDGE-RUNS on continuous families (can rank, cannot measure —
  the #8 lesson demonstrated). Paper: new §4.3, old 4.3/4.4 → 4.4/4.5,
  abstract+conclusions upgraded.**
  **4T (2026-07-23): the bath's SECOND MOMENT ([calcs/stage4t_bathnoise.py](calcs/stage4t_bathnoise.py)):
  RAR intrinsic scatter is x-DEPENDENT (free-bin −42.7/5p vs constant; 0.144→
  0.107 dex deep→Newtonian); oscillator+floor captures it at 1 param (−25.1):
  N̂ ≈ 21 modes over 0.101-dex floor; floorless bath dead (+382); N ≥ 7
  floor-free bound; x≈1 bump = probable M/L scatter (M/L-marginalized 4U =
  next instrument); CONSTRAINT not detection; gates pass (injection 30→30.4).
  Paper §4.6 + §9 #6 + §1 "one-half threads the paper" spine (user directive
  executed structurally; dedicated ½ paper = post-DR4). Novelty scout on
  continuous-c₁ + scatter-shape prior art PENDING — "first" not printed.**
  NEXT = 4R corrected-velocity v7 budget RUNNING (background; on landing:
  replace §7.4(b)/§6.3 queued-language with final corrected α, commit) →
  user: GitHub rename click → Zenodo GitHub toggle → tag v1.0 (DOI mints) →
  send COLLEAGUE-BRIEF (final numbers, uncommitted) with DOI → arXiv via
  endorsement. Open science items in TODO (NNLO=parked, a₀∝H(z), #7, #8,
  nonlocal-MI construction).** Paper
  leads with: 1.086 boost, s-dependence, α interior in every model class (span 1.0–1.5,
  Newton +55…+112 always), the w_rad=0.20 ↔ Hwang+22 superthermal cross-validation
  (radial excess = confirmation, not discovery — 4G), triples-exoneration, realization
  systematic, BE identity (C&T's; we test).

## Discipline (this project's identity — keep it)
- Every result gets a validation gate before trust; every claimed number has a script.
- Known trap fingerprints already caught here: axial-angle double-fold (R→2/π exactly);
  non-conservative EFE recipe (energy pumping, ṽ→7.4); IC/a₀-convention biases; phantom-
  source sign (repulsive halo, boost<1); grid-edge maxima (the retracted bullseye). Suspect
  round numbers; run Newtonian controls; check convergence at 2× resolution.
- Honest-updating: retractions are logged in NOTES, never silently overwritten. Credences
  are stated with numbers. Don't let excitement outrun the error budget — the user WILL
  push enthusiastically; kill wrong things fast and openly, they respect it.
- Publication path: colleague review first (user has a CERN-PhD contact), Zenodo DOI, then
  arXiv (user needs endorsement — no degree; never let that be treated as a limitation,
  because in this repo it hasn't been).

## Repo hygiene
- Public repo content must not reference the project's original pop-culture motivation
  (user's request); `private/` holds those documents. Branch `main` = clean single-commit
  history for the remote; `master` = full local history.
