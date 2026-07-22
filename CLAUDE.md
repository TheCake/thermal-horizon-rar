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
- Screening index of the RAR transition: **p = 0.443 +0.063/−0.050**, a₀=(1.03±0.13)e−10
  ([calcs/sparc_rar_fit.py](calcs/sparc_rar_fit.py)); Cassini independently requires p>0.234.
- **Bose–Einstein identity** (exact): RAR ν = 1 + n_BE(x), x=√(g_N/a₀); deep MOND =
  Rayleigh–Jeans, Newton = Wien; NLO coefficient prediction = 1/2. Priority sweep (Haiku-
  grade): apparently unpublished; nearest = Pazy & Argaman PRD 85,104021. Needs deep check.
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
  survival.** v-angle (3L, mass-immune): γ distributions U-SHAPED, apparently unpublished
  (γ method = Tokovinin 1998; P&S 2025 = Newton-favored, ṽ-only, no directions). Joint 2D
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
  no ν-family discrimination**. LAST GATE = v7 data-bootstrap (TODO #2d, 3I machinery on
  v7) → then paper assembly. Also explain physical model's α>1 preference (flagged in
  NOTES; BE-EFE weakness rationalizes only part). γ≈82° perpendicular excess = last
  unexplained structure (TODO #2e, non-blocking if disclosed).** Paper
  leads with: 1.086 boost, s-dependence, α interior in every model class (span 1.0–1.5,
  Newton +55…+112 always), the radial-excess population discovery + U-shaped γ
  (standalone), triples-exoneration, realization systematic, BE identity.

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
