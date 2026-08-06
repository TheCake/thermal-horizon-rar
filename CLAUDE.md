# PhysicsResearch — session handoff

Low-acceleration gravity research program run jointly by the user (hfili — sharp,
non-academic, wants first-principles reasoning + numbers, enjoys building) and Claude.
Read [NOTES-horizon-inertia.md](NOTES-horizon-inertia.md) (the lab notebook, chronological
with retractions), [TODO.md](TODO.md) (priority queue), [PAPER-DRAFT.md](PAPER-DRAFT.md).

## CURRENT STATE (2026-08-06 — read this first, then LOG.md, then DIARY.md tail)
The fastest catch-up is [LOG.md](LOG.md) (one line per stage, current);
[DIARY.md](DIARY.md) (newest-first) is Claude's own working diary — the
don't-redo lists, trap fingerprints, and the current fork in plain terms;
update it at every round close (same commit as the verdict). The giant
"Key results" narrative below is accurate through the 7-series/PAPER v3.9 but STOPS
there; the 8/9-series (the rivals arc, the control arcs, the QUALITY arc) lives in
LOG.md + NOTES + [LEDGER.csv](LEDGER.csv) (175 rows, six-gate audited, 204 tokens).
Where things stand: **anomaly-real 53** (trail 50→58→53→56→53, every move by
pre-signed map; ROUND 13 = the Opus agent's audit: combiner math CLEAN; his D2
meter then EXECUTED as 9L and resolved the fork NOISE-REAL — EDR3 pair velocity
errors are ≈2× formal at 0.2-2 kAU, quality-independent, a standing
measurement; the binary forward-fit leg = an UPPER LIMIT: α ≥ 0.5 excluded on
clean strata, ≲0.3-0.5 allowed with a mild 0.1-0.3 interior lean, α = 0 not
excluded; B-median statistics carry a real noise pedestal (E2); the in-catalog
binary quality arc is substantially mined out — next discriminators are
external: DR4 epoch astrometry, T2 spectrograph legs). Binary side: the
QUALITY-CONCENTRATION OBJECT is now a measured RUWE DOSE CURVE (9J Q-alone
ladder: Q2 0.05-0.11 < Q1 0.12-0.26 < Q3 0.40-0.79 < Q4 0.90-1.05; M-STD
full-grid 0.39-0.79 interior, seed-split disclosed). OPERATIVE HEADLINE
(round-13 form): on clean strata α ≥ 0.5 is EXCLUDED (−13..−25); α is bounded
above at a deflated value (~0.3-0.5 under fpm-ceiling/fcomp-conversion
treatments); small-real-boost vs Newton-plus-width UNRESOLVED. ROUND-13 E2
CORRECTION: the 9E double-ratio does NOT control differential noise (σ(ṽ)
rises with s; wide arm inflates ×1.155 vs narrow ×1.126 Q1→Q4) — B(Q1)=1.122
= un-decomposed sum, NOT independent boost corroboration (values stand,
counterweight reading withdrawn). 9K (D1): capping fpm at 1.5 raises drop-α
+0.045..+0.081 in 4/4 (real, half the fragile bar); under the cap ALL
law-seeds prefer small α over Newton (BE-101 P(α=0) 0.404→0.026). Galaxy
side: the GD dial-tension survived SEVEN controls — reads genuine slow-dwarf
physics or an untested-class systematic; strongly V-ordered. Data chain: 9H
manifest = SHA256-pinned + 15/15 invariants (re-run stage9h after any
refetch). The "Opus reviewer" = an Opus-model AGENT I spawn (user corrected
2026-08-06: they do NOT relay); brief him fresh each round (no memory), book
reviews as numbered rounds; round-13 text in REVIEW-ROUND13-OPUS.md
(uncommitted). QUEUE (round-13-ranked): D2 the narrow-pair fpm meter (top
binary — the fpm-legitimacy decider), D3 degraded injection (conditional),
LNPI conversion-band row, A4 the GD rotation-curve convergence split (top
galaxy, candidate eighth control), DR4 external, paper thaw at author's call.

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
  **4R DONE (2026-07-23): corrected-velocity budget, six paired seeds —
  α̂ shift −0.012±0.012 (simple) / −0.022±0.015 (BE) = zero within errors;
  Newton cedes ~5 lnL (+103±9/+90±7, min +84); w_rad 12/12; simple-over-BE
  lean unchanged (−12.6±2.4). §6.3 numbers STAND, correction EXECUTED, all
  queued-language deleted; PAPER → v1.3 (Stages 1–4T). PAPER IS TAG-READY.**
  **4U/4V/§2.4 (2026-07-23, PAPER v1.4): 4U hierarchical disk-M/L second
  moment — thermal trend SURVIVES (−22.8/1p; floor 0.101→0.059; N̂ 21.5→63;
  constraint now N ~ 20–60), no galaxy-level draw (r=−0.02); V-shaped free
  bins (Newtonian arm = probable un-marginalized BULGE M/L); G2 injection
  gate FAILED r=0.74 → area/local-mode contest REPORTED UNRESOLVED (honest
  fail, logged). 4V scorecard (§9.1 table): galaxy a₀ vs cH₀/2π +0.1σ/−0.5σ
  (the temperature check PASSES); binary a₀ translation (κ=0.78/0.92) =
  1.48±0.18/1.37±0.17e−10 → +2.5σ/+1.9σ = the sharpest internal tension,
  carried openly; five dials ≲1.6σ. §2.4 = the six-point mechanism
  specification; Unruh scout: Milgrom-1999 + Deser–Levin-1997 VERIFIED and
  cited; the occupation-ratio reading (x = T_U(g_obs)/T_dS) and ω = √(g_N
  a₀)/c = g_obs/c NOT FOUND anywhere = ours, stated as observations.**
  **4W (O1 DONE): full hierarchy = IDENTIFIABILITY BOUNDARY — vertical
  (distance/inclination, measured priors) freedom absorbs the thermal term
  (N̂→189, +24 vs const) but calibrated G2 FAILS (slope 0.23 vs 0.93) =
  channels degenerate on SPARC; floor converged to Desmond's 0.034; **the
  x≈1 transition bump SURVIVES all marginalizations** (0.054 vs 0.026–0.035
  dex) → environmental-g_ext-scatter hypothesis feeds O4; unlock = O1b
  (Cepheid/TRGB-anchored subsample). 4X (O2 DONE): **binary ĉ₁ = 0.37–0.50
  ([calcs/stage4x_binlam.py](calcs/stage4x_binlam.py); λ-tables gated;
  seeds 31/101 peak λ=0.75/1.00; c₁=0 rejected ΔlnL≈20/seed, α edge-riding
  at low λ = shape rejection) — TWO DISCONNECTED SYSTEMS READ THE SAME
  DIAL (galaxies 0.43±0.27) = the program's strongest oscillator statement.
  Paper v1.5: §6.4 new (old 6.4→6.5), scorecard 10 rows.**
  **4Y/4Z/5A (2026-07-23, PAPER v1.6): 4Z hierarchical ĉ₁ RESHAPED the
  claim — hier profile ĉ₁ = 0.258 (0.208–0.309) vs flat 0.450, gates PASS
  (injection unbiased → the ±0.2 spread is real structure); bootstraps
  agree (≈0.4±0.3); **galaxy dial = c₁ ∈ 0.26–0.45, 0 excluded everywhere,
  ¼-vs-½ OPEN (hier profile peaks AT ¼; §4.4 ¼-function contest now
  flat-M/L-conditional → O6)**; binaries stay 0.37–0.50 independent.
  4Y anchored subsample (43 gal): thermal not required there (+0.02),
  identifiability improves not restored (slope 0.42), Newtonian-end spike
  flagged. 5A template contest self-disqualified (injection produces no
  bump — per-galaxy templates degenerate with vertical), yielding the
  sharp fact: **the x≈1 bump is POINT-LEVEL** (survives every per-galaxy
  channel); EFE template absorbs deep bins (0.063→0.045) = the
  thermal-vs-environmental rivalry for the deep trend awaits real g_ext
  data (agent hunting Chae tables).**
  **5B–5H (2026-07-23, the O4→O9 cascade, PAPER → v1.7): O4 (5B/5E,
  convergence-hardened after a self-caught 80-lnL coordinate-descent trap —
  ALL verdict fits now adaptive+nesting-gated): Chae+21 Table 3 per-galaxy
  g_N,ext extracted (94 matched, data/chae2021_table3.csv committed) —
  **thermal SURVIVES the environmental control: β̂=0.044, β=1 (collinear
  max-clustering) excluded +168, credit scramble-generic, correlation null;
  deep-end EFE visibility ≲2–5%** (Desmond-23's "weak EFE evidence"
  convergent; scout: monotone deep-scatter growth still unpublished). O6
  (5C/5D): bath matrix under hier M/L FLIPS — **boot (¼, ν=1+n_BE(νy))
  beats BE by 75.6, truth-calibrated injections (obs gap ≈ boot-truth
  −98, not BE-truth +38), 43/50 bootstrap; simple falls LAST (−99 behind
  BE)**; scout: implicit boot function NOT FOUND in literature (nor
  simple's implicit identity — both stay ours, "apparently unpublished").
  O7 (5F): **binaries VETO boot** (+17–24 behind ½-branch, α̂ edge-rides
  = shape rejection of ν(1)=1.35 transition; Newton dead +75/+81). O8
  (5G): **the decomposition — 5D's ladder is monotone in Newtonian-tail
  sharpness; ν_p (≡BE at p=½) hier-converged gives p̂≈0.65, −56.4 = 75%
  of boot's flip from the TAIL alone** at ½-grade transition (ν(1)=1.42);
  convergent with §3's p=0.578±0.12 → the hier galaxies vote for SHARPER
  SCREENING (Cassini-friendly direction); ¼-vs-½ deep digit stays OPEN.
  O9 (5H, unification): binaries ACCEPT ν_p(0.65) — α̂ interior, Newton
  +93/+91; firmed by O10.
  **O10/O5/O11 (5I–5M, 2026-07-23 night, PAPER → v1.8): 5I quadrupole =
  AMPLITUDE-LOCKED (raw q drops 0.81×/0.85× but binary α̂ rises in exact
  proportion; Q₂·α̂ ≈ 4.2–4.4× Cassini for EVERY function — no
  interpolating-function escape; MI door unchanged). 5J THE CONSTRUCTION:
  geometric-mean bath ν_gm = 1+n_BE(y^¾√ν) (ω=√(ω_source·ω_total), ZERO
  params, derived c₁=⅓ c₂=1/12 gate-verified, ν(1)=1.433, tail p≡¾) —
  LEADS hier galaxy ladder −10519.8 (−85 vs BE), ties flat (Δ2.0);
  priority scout NOT FOUND + Pazy–Argaman closed by direct read (FD
  stats, O(e^{−1/x}) deep corrections — cannot anticipate; gitignored
  data/pazy_argaman_text.txt). 5K binary budget (4 seeds × 4 functions,
  regression exact): binaries prefer p=½ SIGN-CONSISTENT 12/12 (+5.8…
  +10.4 per seed); sharp functions viable (α̂ interior: p065 1.40±0.04,
  gm 1.44±0.08); Newton dead +77…+98. 5L a₀ ladder (κ gates pass): sharp
  binary a₀ = 1.58±0.11 = +4.9σ off cH₀/2π (vs +1.9σ at ½); galaxy flat
  legs ON horizon (gm −0.2σ). 5M vertical-channel control (measured σ_v,
  G1 exact, nesting OK): **the a₀/f_ML anomaly RESOLVES — all hier fits
  return to a₀ = 1.04–1.13e-10 = horizon (the program's strongest
  temperature-lock statement); tail vote HALVES but survives (gm −41,
  p065 −32); boot COLLAPSES (−9) = the ¼ cell dead everywhere.** END
  STATE: the ladder digit is three-cornered — ½ (binary-anchored, mildest
  a₀), ⅓ (gm: best cross-system profile, leads every galaxy treatment),
  ¼ (excluded).**
  **5N/5O/5P (2026-07-23 close, PAPER → v1.9): THE MIXING FAMILY (5P, the
  "why" instrument) — ω = ω_src^(1−β)·ω_tot^β ⇒ EXACT c₁ = 1/(2(1+β)),
  c₂ = 1/(12(1+β)) + β/(8(1+β)²), p_tail = (1+β)/2, relation c₁·p_tail =
  ¼; β=0/½/1 = BE/gm/boot (member regressions 1e-13); Bernoulli 1/12
  recurs at β=½; β=½ = unique exchange-symmetric fixed point (impedance-
  matching / Curzon–Ahlborn / source-stiffness+response-inertia readings —
  stated as readings, not derivations). MEASURED: galaxy β̂ = 0.64 (plain
  hier, Δ+86) / 0.45 (vertical-hardened, Δ+43, implied c₁ = 0.35 ≈ ⅓ —
  ON the symmetric point); binaries hold β ≈ 0. Final two-system form:
  β ≈ ½ (galaxies) vs β ≈ 0 (binaries) — one falsifiable parameter;
  configuration-dependent β (e_N-dependent mixing) = reconciliation
  hypothesis → O13. 5N bootstrap: vertical-robust gm lead = −29 ± 53,
  29/40 (72%) = a LEAN, quoted as such. 5O six-seed: occupation preferred
  18/18 (p065 +5.2±0.9, gm +8.5±2.1); α̂ 1.078±0.023 / 1.41 / 1.36 all
  interior; sharp a₀ 1.59 (+5.1σ) / 1.51 (+4.3σ) vs ≈ +1.8σ at ½.
  O5-remaining = derive WHICH β from horizon microphysics.** Near-miss logged: 5A's tanh-EFE template had an
  ambiguous curvature sign (5A unaffected — free-sign coefficients);
  5B/5E use Chae+21's exact Eq. (2), gated.
  **5Q–5T (2026-07-23, O13 COMPLETE, PAPER → v2.0): 5Q c₃(β) =
  β(3β+1)/(24(1+β)³) — β=0 is the family's UNIQUE Bernoulli zero
  (c₄(0)=−1/720 exact); rescaled ladder polynomial (c₁p=¼, c₂p²=
  1/48+5β/96, c₃p³=β(3β+1)/192); third log-cumulant c₃−c₁c₂+c₁³/3 =
  β(2β−1)/(16(1+β)³) VANISHES at β=½. 5R binary β-profile (b025/b075
  tables gated 1.3e-15/0.01%, boot completed; 6 seeds × 5 β): **β=0
  beats every β>0 24/24; mean Δ −4.1/−8.5/−7.9/−15.4; bounds β<0.030
  (1σ), β<0.121 (2σ)** — sharp; b075 edge-rides 5/6 (shape rejection).
  5T galaxy decomposition (regression d=0.00 all nodes; identity 1e-6):
  **the galaxy β̂=0.5–0.64 is a COMPROMISE — ultra-deep votes AGAINST
  β>0 (+9.3 at ½, wants c₁=½), tail carries it (−61 at y>1), transition
  flips against at β=1; HIGH-arm free β̂=0.76 interior, LOW-arm free
  fit ridge-flagged (f→1.94, not interpreted). β RUNS WITH REGIME:
  deep+transition→0 EVERYWHERE (binaries, SPARC-deep, SPARC-transition),
  tail→½–¾ = the lone β>0 voter. The c₁·p_tail=¼ LOCK is what the data
  strain against (want c₁≈0.4–0.5 AND p≈0.65–0.75 → ~0.3). MI-mimicry
  reading of the split DEMOTED (intra-SPARC split is all circular
  orbits = MI≡MG there).** 5S quadrupole scan: q monotone −0.099→−0.073
  but **every member 4.0–5.8× Cassini (edge members = lower bounds) —
  Saturn's veto is β-blind**. NEW TODO O14: the lock-breaking function
  (running β(n): 0 at n≫1 → ½ at n≪1, tail p=¾) — testable in both
  pipelines as they stand. O5 sharpened: derive the classical→quantum
  β-switch.
  **5U–5X (2026-07-23, O14 EXECUTED = the derived function, PAPER →
  v2.1): 5U the spontaneous-fraction bath — response admixture weighted
  by the spontaneous channel (Einstein coefficients): β = ½/(1+n) =
  1/(2ν) (F1) / energy-share β = 1/(2(2ν−1)) (F2); ZERO parameters;
  c₁=½ EXACT, p=¾ EXACT, lock→⅜, c₂ = −1/6/−1/24, ν(1)=1.470/1.494;
  all gates (sympy exact, mpmath 50-digit, residual 9e-13, uniqueness);
  C&T's β=0 = the CLASSICAL limit; exchange-symmetric ½ = the quantum
  endpoint (only spontaneous processes remain). Scout: running
  exponent / Einstein-coefficient derivation / (½,¾) combo / C&T
  follow-up all NOT FOUND (scout-level). **5V: F2/F1 LEAD THE ENTIRE
  GALAXY LADDER both treatments (−107.2/−101.3 plain, −51.9/−44.3
  vertical vs BE; gm −84.8/−42.7), BE regressions d=0.00, a₀
  horizon-adjacent. 5W: binaries −5.5±2.3 (F2) / −5.8±1.9 (F1) vs gm
  −8.5±2.1 — a third of the gap closed, interior 12/12, α̂ 1.487±0.017
  / 1.395±0.039, pre-stated bands straddled = improved NOT accepted
  (residual ~2.4 SE, localized at the transition = the refinement
  target). 5X: Q₂ 4.8×/5.3× Cassini — lock holds.** The
  spontaneous-fraction bath = the program's best cross-system function
  (indicative joint −41 vs gm −26), its two predictions pre-registered.
  OPEN: O15 (κ/a₀ ladder for F1/F2 — temperature-strain adjudication;
  faster-die-off weighting refinement; 5N-grade bootstrap of the F2
  lead).
  **5Y/5Z/6A–6D (2026-07-23 night, O15+O16 EXECUTED, PAPER → v2.2):
  THE SPLIT IS SYSTEM-LEVEL. 5Z two-leg functions (β=1/(2ν²),
  1/(2(2ν−1)²)): c₁=½ AND c₂=1/12 EXACT, ν(1)=1.503/1.537, all gates.
  6A: galaxies reward again — F3 −111.4/−50.9, **F4 −108.7/−64.2 =
  biggest controlled lead ever**, a₀ horizon-adjacent; 6C bootstrap
  **F4−BE = −57.4±38.3, 37/40 (92.5%)** = strongest function grade in
  program (still a lean). 6B: binaries REJECT the two-leg too (−6.84
  interior / −6.96 with α̂-EDGE 6/6) — across EIGHT sharpened functions
  the binary penalty −5…−8.5 is ~independent of ν(1): **they reject the
  screening behavior under the dominant external field, not the
  transition.** 5Y κ ladder: κ=1.26–1.35 all running fns → a₀ +5.2σ/
  +6.3σ/+6.2σ (F4 edge-n/a); ONLY pure BE passes (+1.9σ). 6D the
  decisive control: pointwise drive-weighted rule (β·g_int/(g_int+g_ext),
  solver-level, gates exact: isolated→F4 0.01%, wide→BE 0.35%)
  **EXCLUDED — −9.64±1.49, 0/6, punished exactly at its predicted
  mid-separation sag — while PASSING temperature (α̂=1.098±0.021,
  κ=1.033, a₀=1.31±0.13 = +2.1σ): the two binary vetoes are
  independent.** END STATE: two systems carry DIFFERENT effective baths
  (BE under dominant ambient field; two-leg-sharpened when isolated),
  difference NOT a function of local field configuration. Surviving
  (post-hoc-flagged): ambient-GATED bath (test: Chae high-e_N galaxies,
  power ≲2–5%; real test DR4 weak-ambient binaries → should sharpen) vs
  MI trajectory functional (split follows orbit class; DR4
  eccentricity-resolved). TODO O15 closed, O16 opened, O5 sharpened to
  "why would the ambient gate at system level — or is the trajectory
  the carrier".
  **6E/6F/6G (2026-07-23 night, O5→O16 RESOLVED FIRST-PASS, PAPER →
  v2.3): THE AMBIENT-GATED BATH — the theory question answered with a
  SIGN LESSON (weak fields = deep = classical occupied ambients;
  n_amb = n_BE(√(e_N/a₀)): gal 6.6, bin 0.52) and the rule with the
  right sign: self-consistent dressing needs STIMULATED reservoir
  assistance → **β = [n_amb/(1+n_amb)]²·½/(2ν−1)² — admixture = local
  quantumness × ambient classicality, ZERO parameters, system-level
  (explains the 6D no-sag pattern)**. Exact: c₁=½, c₂=1/12
  g-INDEPENDENT, c₃=−g/16; **p = ½+g/4 POST-DICTS both measured tails
  (0.689 gal vs 0.65–0.75; 0.529 bin vs ½)**; ν(1)=1.548/1.577; g=0→BE
  1e-12. Scout: env-dependent SHAPE / decoherence gating / Einstein-
  coeff modulation / two-function MOND all NOT FOUND. Pre-registered
  tests: **6G binaries ACCEPT (−0.88±2.66, 2/6 prefer, interior 6/6,
  α̂=1.060±0.024, κ=0.924, a₀=1.28±0.15 = +1.6σ — BEST temperature row
  in program, better than BE's +1.9σ); 6F galaxies −59.05 vertical
  (PASS ≤−50 bar, 2nd-best ever, dilution cost only 5.2 vs F4) /
  −92.4 plain (MISS ≤−100 bar by 8, disclosed; band −82.6…−97.4).
  JOINT −57.3 vs BE = THE FIRST SINGLE FUNCTION TO PASS BOTH SYSTEMS**
  (F4 −50.3 binary-vetoed, gm −25.7, BE 0). Post-hoc flag carried
  everywhere; non-guaranteed content that landed: vertical survival,
  κ/a₀ +1.6σ, both p-postdictions. Out-of-sample: Chae g_amb(e)
  softening dial, DR4 weak-ambient sharpening toward p=0.69, MI
  eccentricity discrimination. O16 remaining: Chae leg, AMB-lead
  bootstrap, AMB quadrupole (expect BE-grade), DR4. O5 remaining:
  horizon-side derivation of the share-squared grammar. PUBLICATION
  still DEFERRED; O1b anchors parked; Zenodo/colleague (#11/#12)
  parked until called.**
  **6H/6I/6J (2026-07-23, the grammar round, PAPER → v2.4): 6H
  ([calcs/stage6h_grammar.py](calcs/stage6h_grammar.py), all gates
  first-run): grammar β = ½·[q_loc·s_amb]^L; JC dispersive pull =
  λ²(2n+1)/Δ EXACT ⇒ zero-point share 1/(2ν−1) is derivation-grade
  (selects F2/F4 over F1/F3; "local quantumness" = vacuum fraction of
  the frequency pull, deterministic, user's randomness worry answered);
  **THE BERNOULLI-BREAK RUNG IS THE LEG COUNT** (L=1 kills c₂=1/12−s/8;
  L=2 keeps c₂, c₃=−s²/16; L=3 keeps c₂+c₃; c₄(0)=−1/720 all L); NEW
  rung c₄(L=2) = s²/192 − 1/720 (sign flip at galaxy gates, future
  falsifier); tails p=½+s^L/4; Chae 109-galaxy gates MEASURED sharper
  than fiducial (median g 0.868 maxclust/0.952 noclust vs 0.754);
  6I bars pre-registered on disk. 6I
  ([calcs/stage6i_chaegate.py](calcs/stage6i_chaegate.py), regressions
  d=−0.00/dd=+0.00 exact): **L = 2 point-preferred (vertical: L1
  −52.76, L2 −59.05, L3 −54.52) — but see 6L/CORRECTION #13 below; measured ambients IMPROVE both
  treatments (vertical PGmax −61.68; plain −100.51 maxclust/−105.94
  noclust = the 6F ≤−100 bar RESOLVED — the miss was the fiducial-e
  artifact); joint ledger ≈ −59.9.** Honest residuals: galaxy leg
  alone still monotone toward F4 (+2.5 ahead vertical) — the gate
  amplitude is pinned by the BINARIES (two-system content);
  noclust-vs-maxclust ≈ 5 lnL input ambiguity carried. AMB quadrupole
  RECORDED: Q₂ = 3.60e−26 = 4.0× Cassini (amb q within 0.4% of BE; MI
  door standing). Scout (Haiku): Bernoulli-break↔leg-count and (2n+1)
  selection both NOT FOUND (scout-level). 6J bootstrap of the PGmax
  lead ([calcs/stage6j_ambboot.py](calcs/stage6j_ambboot.py), 40
  paired reps, rng 53, comparators 6C F4 −57.4±38.3 37/40 / 5N gm
  −29.3±53.0 29/40): **DONE — AMB(pgmax) − BE = −56.71 ± 35.65, 37/40
  (92.5%), percentiles [−88.7, −59.9, −15.7] = the program's top
  bootstrap grade, tied with F4, binary-compatible where F4 is not; a
  strong lean, not a detection (3/40 flip).**
  **6K/6L (2026-07-23 close, PAPER → v2.5): 6K desktop analog
  ([calcs/stage6k_analog.py](calcs/stage6k_analog.py), pre-registered
  at ddc83bc BEFORE execution; lab-native identity β =
  ½·tanh²(x_loc/2)·e^(−2x_amb) gated): **rate-based realizations of
  the grammar EXCLUDED** — vanilla two-channel = constant β (sky-dead
  per 5T), mediated jump class = wrong-sign running + G4 estimator
  pathology (exponents unquotable); κ-dependence 0.3–0.5/decade both
  configs = the grammar's rate-freedom does NOT emerge from jump
  competition; **pre-committed strike LOGGED: bath-microphysics
  conditional ~20–25% → ~15%, surviving mass on the coherent-pull
  reading; NEXT = 6K-v2 coherent susceptibility calculation (driven
  Kerr Liouvillian + quantum regression — the object the JC anchor
  actually describes), now the sharpest O5 item.** 6L leg-count
  bootstrap ([calcs/stage6l_legboot.py](calcs/stage6l_legboot.py)):
  **CORRECTION #13 — the "L=2 measured" claim DEFLATES: d(L1−L2) =
  +9.9±21.4 (29/40 = lean), d(L3−L2) = +1.6±13.1 (21/40 =
  UNRESOLVED), L2 strictly best 10/40; reframed everywhere as
  point-preference (the deep arm of 153 galaxies too thin to read the
  c_{L+1} rung at population grade). Unaffected: 6H correspondence
  (exact math), Part B measured-ambients + 6J 37/40, quadrupole
  record, 6K strike. Paper → v2.5, Appendix A = thirteen
  corrections.**
  **6M (= 6K-v2, 2026-07-24, pre-reg 2b3bf78 before execution): the
  structured-bath calculation — ONE thermal bath, finite spectral
  resolution b (Ohmic-Gaussian golden-rule kernel), κ cancels EXACTLY
  (fixes v1's rate-ratio objection); endpoints gated (Davies 0.0e+0,
  source 0.67%). VERDICT: AMBIG per bands — share test FAILED 205×
  where clean, resolution collapse failed 5.3× because the estimator
  drowns (endpoint gap K-order vs kernel systematic b²-order;
  overshoot columns flagged). Beyond-band observations (no committed
  credence move on AMBIG; ~15% HELD, lean noted): (a) clean-region
  monotonicity OPPOSITE to grammar — occupied modes sit fully
  dressed, sparse modes source-pin first (resolution physics:
  spread K·√(n(n+1))); (b) resolution-β in gravity = ρ = x(√ν−1),
  transition-peaked (0.26 at x≈1–1.5), tail-zero = wrong sky shape
  in both non-deep regimes. Pull lemma exact: thermal comb centroid
  = 2K·n̄ (statistics, not shares). 6N (=6K-v3, pre-reg 6724597):
  **THE LAB LEG CLOSED — CLOSE-OPP. The cancellation theorem: flat
  bath ⇒ exact source-locking (full two-mode Lindblad confirms:
  max|λ| = 0.0065, g-independent — the buildable fridge config
  measures ZERO); which occupation a self-shifting mode takes =
  whether the bath's KMS structure is resolved across the self-shift
  (C&T-vs-boot in one line). KMS-contrast carrier works but is
  bath-geographic (34× off tanh²; λ rises with x₀ — recorded, not a
  rescue); contrast dial's sky shape transition-peaked/tail-vanishing
  ⇒ the sky's tail-only β is ANTI-STANDARD. Three classes (6K/6M/6N)
  all fail ⇒ bath-microphysics conditional 15% → ~8% per pre-commit;
  NO v4; fridge moot. Survivors: horizon-specific non-Markovian, MI/
  trajectory, or pure sky-phenomenology (~92% complement — function
  record untouched). Paper §2.4 closed-out in place (still v2.5).**
  **6O/6P/6P-b (2026-07-24, the frozen-bath tests, both pre-reg'd):
  the frozen sub-reading (bath correlation ~1/H ≫ orbital times ⇒
  quenched draws; mode wavelength ≫ galaxy ⇒ coherent) given its two
  accessible tests. 6O galaxy coherence: AMBIG per gates — estimator
  corrected post-commit (exact z-marginalization Occam determinant;
  the MAP version was monotone-degenerate, its own G0 caught it),
  then the data-side failure: profile edge-runs −835, injection
  over-recovers, SHUFFLED template still gains −82 ⇒ per-galaxy
  radial systematics dwarf any draw (x-shape beats shuffled 10× but
  non-specific). 6P/6P-b binary s-shaped scatter (v7 patch-runner,
  paired SQN grid, seeds 31+101): edge-PREFERRED +11.1/+13.5,
  monotone to N=2 (+23.5), BUT the s-FLAT matched control WINS by
  25.9 lnL (flat-8 = +37.0) ⇒ GENERIC broadening absorption, NOT the
  frozen fingerprint. v7 error-model finding logged: an s-flat
  per-system scatter worth ~+37 lnL exists (3E/3J echo); α̂ exposure
  ≤ +0.08 = inside the ±0.11 systematic — the measurement stands.
  Frozen sub-branch: no support, no formal strike, both near-term
  signatures spent; ~8% holds. O5 remaining: horizon-side
  non-Markovian derivation; MI-trajectory (DR4); phenomenology.**
  **6Q (2026-07-24): THE MEASUREMENT LEDGER + WORLD TABLE —
  [LEDGER.csv](LEDGER.csv) (75 hand-curated rows: value in NOTES
  convention, stage, script, output, data-deps, status CURRENT/
  CO-QUOTED/SUPERSEDED/RETRACTED + pointer; spec [LEDGER.md](LEDGER.md);
  hand-curation BY RULE — scraping would resurrect α=1.54/+264/e_N=1.9/
  "SPARC leans simple"/"L=2 measured", all kept as marked rows) +
  [calcs/stage6q_worldtable.py](calcs/stage6q_worldtable.py) (six gates
  ALL PASS: provenance-on-disk, supersession resolution, 15 verbatim
  value greps vs stage outputs; gates caught 3 transcription slips
  first run). World table (13 laws × 7 tests): formal vetoes = Newton,
  F4, boot, drive-weighted, no-EFE MI; no-veto = simple (hier-rejected
  −99 strong-lean), BE, p065, gm, F1/F2, F3 (all binary-lean +
  a₀-strained), AMB, mi_t; **mechanical sentence: AMB uniquely (a)
  binary-unvetoed (tie), (b) top galaxy bootstrap 37/40, (c) best a₀
  row +1.6σ, (d) postdicts both tails — post-hoc flag + shared 4.0×
  Cassini carried; independence explicit (2 data votes + Cassini);
  ladder digit ½/⅓/sharper stays OPEN**. STANDING RULE: every stage
  with a headline number adds its ledger row in the same commit;
  supersede = flip status + pointer, never delete.**
  **6R/6S/6T (2026-07-24, O5 horizon-side round, pre-reg 85dcc72): THE
  RESOLUTION BATH — 6N corollary + one-scale dS bath (KMS width =
  temperature = the a₀ scale) FORCE R = νy−√y; β = ½·R²/(1+R²) in the
  5P mixing (Lorentzian = flagged representative choice). Exact+gated
  ([calcs/stage6r_resolution.py](calcs/stage6r_resolution.py)): ladder
  preserved through c₄, **break c₅ = −1/16** (deepest preservation;
  AMB breaks at c₃); tail = gm argument (p=¾); ν_R(1)=1.539; crossover
  y* = 1.88 = the 5T β turn-on arm; Deser–Levin matrix CLOSED (source
  branch: no finite fixed point; total: boot + invisible n_BE(2π)
  floor; **geodesic detector T_U = 0 exactly = the free-fall
  derivation of binary β<0.03**). SKY: 6S galaxies STRONG — **−58.59
  vertical (AMB grade) / −113.72 plain = largest plain lead on
  record**, BE regressions d=−0.00; 6T binaries **SHAPE REJECTION**
  (α̂ edge 2.00 in 5/6, 0/6 prefer, −7.38±2.08 = the eight-function
  band, Newton +82.7±1.8, κ=1.351 a₀-row edge-invalidated). **THE
  TRIANGULATION COMPLETE: local-resolution-only dead (6T) +
  ambient-pointwise-only dead (6D) + two-factor accepted (6G) ⇒ the
  ambient gate is REQUIRED — the two-system split is system-level,
  not local-in-y. O5 sharpened to: derive the s_amb² gate
  horizon-side (the local factor's deep half is now derivation-grade:
  free-fall + one-scale resolution).** Ledger rows galfn-resn/
  binfn-resn/mech-resolution added; world table has the RESN row +
  triangulation line (17 spot-checks, gates green); PAPER §2.4
  paragraph + App B (ledger cited, "thirteen" fixed).**
  **6U (2026-07-24, same day): THE GATE DERIVED
  ([calcs/stage6u_gatederiv.py](calcs/stage6u_gatederiv.py), gates
  first-run) — path-resolving the 6H JC loop splits the (2n+1) pull
  into absorption (n) / emission (n+1) time-orderings; per-leg ratio
  n/(1+n) = e^(−x) EXACT = the KMS detailed-balance weight; **s^L =
  the Boltzmann cost of L borrowed ambient quanta, L=2 = the loop's
  own order — ONE loop, two vertices, both grammar factors per vertex
  (local zero-point share × ambient KMS ratio)**. Uniqueness on
  measured facts, no new fits: rate-balance s²/(1+s²) (p_gal 0.607)
  + absorption-SHARE² (0.554: **data choose detailed balance over the
  share from the same decomposition**) + inverse² (both tails, wrong
  sign) + n² (absurd) excluded by tail bands; pointwise by 6D;
  local-only by 6T; **s² = unique survivor**; postdictions regress to
  6E (0.6884/0.5280). Dividends: DR4 source-vs-dressed x_amb rung
  (Δp≈0.025 at e_N=0.4); **the tail exponent RUNS with H(z)** (p_gal
  0.688→0.702 at z=1; gate now inside kill test #14). Scout (Haiku):
  KMS gate / Boltzmann EFE / reservoir assistance / path-resolved JC
  in gravity / C&T env follow-ups — ALL NOT FOUND (scout-level;
  near-misses logged). Honest seam carried: the borrowing NARRATIVE
  stays reading-grade (non-Markovian horizon microphysics); 6O/6P
  variance nulls do not strike it (mean-level). Ledger row mech-gate;
  world-table 18 tokens green; PAPER §2.4 gate paragraph + App B
  script map. AMB = a DERIVED structure with one reading-grade seam,
  not a post-hoc rule.**
  **6V (2026-07-24, the 6U falsifier, pre-reg 6ebf545,
  [calcs/stage6v_untied.py](calcs/stage6v_untied.py)): untied contest
  β = ½q^L1·s^L2 on the vertical ladder (BE d=−0.00, tied-cell
  dd=+0.00) — margins vs tied (2,2): U21 +3.33, U32 +0.74, U23 −3.72,
  U12 −7.11. **B1 TOLERATED (tied survives — no strike, no proof; grid
  flat inside 6L's ±13–21 noise). B2 FAILED as pre-registered and
  LOGGED (my instrument model wrong: the likelihood reads L1 through
  the TRANSITION-suppression channel, not only the deep rung —
  dividend: one-leg disfavored −7.1 via a second, rung-independent
  channel, convergent with 6L's 29/40 lean; L1=3-vs-2 blind).**
  Ledger row gal-untied (19 tokens green); PAPER → v2.6 (header names
  the full 2026-07-24 arc: lab closures, frozen tests, LEDGER,
  triangulation, KMS gate + rungs, untied contest; abstract §7 §2.4
  conclusions all carry the derived-gate status with the one seam).**
  **6X (2026-07-24 afternoon, pre-reg 67626f7, user-prompted
  recalibration "decompose the months-scale item and try",
  [calcs/stage6x_borrow.py](calcs/stage6x_borrow.py)): THE BORROWING
  DYNAMICS — SUPPORT, the first CONSTRUCTIVE mechanism result (after
  three exclusion rounds). Frozen-horizon limit = closed dynamics
  (no bath during evolution); Kerr mode + one thermal ambient mode;
  dress-ward transition resonant (must absorb a real quantum).
  **Dress-ward weight = ½·P(n≥1) = ½·n/(1+n) to 0.2–2.2% over n_amb
  0.25–8; ratio slope +0.989 rms 0.0046 (6.2×/29× over share/raw) —
  SUPPORT bar cleared. GL exact: P(n≥L) = s^L = e^(−Lx) — the KMS
  ratio IS the lending probability; the gate = P(ambient can supply
  the L quanta). Resonant channel λ-independent (0.997) = REAL
  exchange; detuned virtual control λ² + raw-n — the sky's ratio²
  selects the real channel (frozen bath cannot absorb virtual
  imbalances). Dephasing closes the channel (6N direction).**
  Amendments post-commit PRE-RESULTS logged (linear modes carry no
  pull — measured; per-Fock anchor now 0.998–1.000; NB budget).
  Per-leg factor DERIVED at toy grade; seam = dS-side reservoir
  identification + multimode; ceiling ½ candidate reading = resonant
  equal-time-sharing (flagged). Credence HELD ~8% (only the STRIKE
  branch was pre-committed); lean now constructive. Ledger row
  mech-borrow (20 tokens green); NOTES entry; PAPER §2.4 paragraph.**
  **6W (2026-07-24, pre-reg 2193672, the Saturn round,
  [calcs/stage6w_scalarefe.py](calcs/stage6w_scalarefe.py)): SCALAR
  vs VECTOR EFE on the binaries — my scalar-EFE hypothesis (direction-
  blind thermodynamic EFE ⇒ spherical phantom shell ⇒ shell theorem ⇒
  Q₂ = galactic tide ~8e-31 = 12,000× under Cassini) **EXCLUDED by
  the data: scalar-BE trails vector-BE −10.35±2.03, 0/6, interior 6/6
  = valid shape rejection (the program's strongest single-model
  exclusion); scalar-AMB −11.24±1.90.** The moot twist recorded: the
  scalar temperature row would have been the best ever (α̂=1.002±0.049,
  a₀=1.20±0.16=+1.0σ). Scout: scalar-EFE class = apparently novel,
  never contested anywhere; Chae+21 PRIMARY READ: their EFE detection
  is magnitude-level, direction "has a minor effect" (1/6 modulation)
  — the scout's "directional detection" claim was a mischaracterization
  (4th scout misread caught by primary source). **THREE DIVIDENDS:
  (1) FIRST measurement of the EFE composition character (strictly:
  the vector composition's angle-averaged radial profile — our
  likelihood is direction-blind); (2) the Cassini tension is now
  COMPOSITION-locked as well as amplitude-locked (data-forced, not
  formulation-chosen); (3) #8's job description sharpened: reproduce
  the angle-averaged vector profile at e≈1.2a₀ while suppressing P₂ —
  ad hoc unless derived; MI remains the standing escape.** Caveat:
  ν(y+e) class tested, not every magnitude-only rule; cluster-
  lopsidedness literature = paper-level engagement item. PAPER §8.1
  composition paragraph + App B; ledger row binfn-scalarefe (21
  tokens green); NOTES verdict: my hypothesis EXCLUDED in ninety
  minutes — Saturn SHARPENED, not solved.**
  **6Y/6Z (2026-07-24, pre-reg 0e0f4fc): THE RESERVOIR IDENTIFICATION
  ([calcs/stage6y_reservoir.py](calcs/stage6y_reservoir.py), gates
  first-run): EXCLUSION THEOREM — one collective thermal mode gives
  the measured s^L; M≥2 democratic modes give the negative-binomial
  tail 1−(1−q)^M(1+Mq) which saturates the galaxy gate at M=2 (0.952
  vs measured 0.754), washes out the 6I-measured e_N-dependence, and
  pushes the binary postdiction to the rejected edge (0.565) ⇒ **the
  gate data select M=1: one collective ambient mode per system**
  (counting-statistics route converging with the dynamical 6D/6T/6G
  route). Identification (reading): the BARYCENTRIC mode = the
  environment's dressing cloud, n_amb = ν(e_N)−1; e→0 rejoins the
  horizon soft sector (gate→1); ambient GAPS the soft sector.
  **PREDICTIONS: P1 tail ceiling p ≤ ¾ EXACT (parameter-free; void
  asymptote ~0.72; one galaxy beyond ¾ kills the gate); P2 p(e_N)
  ordering (→6Z); P3 z-rung; P4 DR4 pair; P5 nested-ambient rule.
  SATURN COROLLARY (reading, 3 measured legs): collective-mode
  coupling = trajectory-state = the EFE-respecting MI class — ties
  MG on binaries (4L), meets the 6W composition demand, no capped
  quadrupole ⇒ the thermal rule lands in the one door Saturn left
  open BY DERIVATION; kill test = DR4 eccentricity resolution.**
  Ledger row mech-reservoir (22 tokens green); PAPER §2.4
  identification paragraph. 6Z = the ordering shuffle test (TRUE
  gates vs 8 permutation nulls, distribution-preserving): **ANTI
  FIRED — true −61.62 loses to ALL 8 shuffles (null −62.91±0.80;
  empirical p≈0.11) ⇒ CORRECTION #14: the 6I measured-ambient gains
  are gate HETEROGENEITY (generic under permutation), NOT ordering —
  the shuffle control does to 6I what the s-flat control did to 6P.
  In-sample environmental-ordering leg WITHDRAWN; paper claim sites
  annotated (App A = fourteen corrections); amplitude-level legs
  (two-system split, 6G acceptance, postdictions, p≤¾ ceiling, M=1
  theorem, 6X mechanics) UNTOUCHED; ordering prediction now fully
  out-of-sample (DR4) where failure = REAL strike, pre-stated.
  Ledger rows gal-ordering + gal-chae-ambients note (23 tokens
  green).**
  **7A (2026-07-24, pre-reg 21e7d28): THE EINSTEIN FLUCTUATION TEST +
  [PREDICTIONS.md](PREDICTIONS.md) (the signed prediction ledger:
  registered-before-test, status-flips only; P1–P8 with uniqueness
  classes + numeric kill conditions; A = the record incl. banked
  anti-classical votes — Boltzmann tail + c₁=½ exclude the power-law
  classical bath; D = non-claims, bunching alone is classical). P6 =
  Einstein 1909 on the RAR scatter: Var(n) = n̄+n̄² (quantum γ=1) vs
  n̄² (classical wave γ=2) vs n̄ (shot), a₀ locked to the mean, same
  param count ([calcs/stage7a_einstein.py](calcs/stage7a_einstein.py)).
  EXECUTED same day: **shot bath EXCLUDED (−25.05, collapses onto the
  floor); quantum-vs-classical NO VOTE (EQ−EC = −0.43; γ profile
  BIMODAL, both edges beat the middle, axis ≤4.4 units); G2
  calibration PASS (paired injections recovered γ̂ 1.03/2.05) ⇒
  misspecification not power — THE x≈1 BUMP OCCUPIES THE DISCRIMINANT
  WINDOW (laws differ maximally x≈0.4–1.2; bump 0.8–1.4). Scatter
  x-shape = deep-concentrated excess + bump + floor (bump-modeled
  γ̂→bound, b=0.083, beats 6-bin free by 10).** P6 LIVE both
  directions; unlocks: bump source (the 4W mystery now ALSO the γ
  blocker), 7B hier-hardened γ (4W caveat: may stay unidentifiable at
  SPARC depth), DR4. Flagship out-of-sample unique prediction = P3
  z-LOCKED PAIR: a₀(z) = cH(z)/2π AND p_gal(z) run together (0.689→
  0.702 at z=1; table in PREDICTIONS.md) — one function, two locked
  drifts, no classical-MOND/DM mimic. Ledger rows galfn-shotbath +
  gal-scatter-gamma (25 tokens green); NOTES 7A; PAPER §4.6 closing
  paragraph; plain verdict NEEDS REFINEMENT (instrument+ledger =
  SUCCESS-grade; headline question blocked by the bump).**
  **7B/7C (2026-07-24, THE BUMP HUNT, pre-reg d82cc4b/59c582a):
  the bump is CAUGHT — INNER-DISK (R < 1.5 R_d) astrophysics, not
  the law. 7B matrix (composition/coherence/geometry at fixed x,
  gates all PASS): verdict MIXED per bars but two discoveries —
  (1) within-curve lag-1 ρ = 0.876/0.844 MEASURED (perm null ±0.10):
  smooth radial misfit dominates raw scatter, every point-level
  −2lnL in the program nominally calibrated only, offsub point floor
  0.04–0.07 dex ≈ Desmond 0.034; (2) scatter INNER-DISK-organized at
  fixed x (R<1.5R_d ~2.4× outer, both slices; window = 49% inner vs
  16% deep-control → radius mix); composition structurally powerless
  (window f_* terciles 0.986/0.998). 7C confirmation: **b_clean =
  0.0000 on outer points (locality DECISIVE; 7A full 0.083);
  EXPLAINED = 0.67 PARTIAL (miss vs 0.75 bar disclosed; residual =
  x-dependence of the inner excess); γ contest with inner term
  (ĉ_in=0.105 ≈ predicted 0.114, worth ~83 lnL/1p): profile
  de-bimodalizes, interior γ̂=0.52, γ=2 nominally +12 — but G2 FAILS
  on the informed design, clean-subset variant leans OPPOSITE
  (+4.48), thinning +0.14 ⇒ NO VERDICT (systematics, not physics).
  THE CLOSE: γ NOT MEASURABLE on SPARC at any grade tried; P6 fully
  out-of-sample (anchored/DR4/IFU + non-circular modeling); neither
  kill direction fired — the quantum reading UNTESTED in the scatter
  channel, not struck.** Scout: the acceleration-binned scatter
  measurement (4T/4W) apparently unpublished; streaming 10–40 km/s
  documented. Ledger: gal-curve-coherence + gal-bump-inner +
  gal-gamma-final, gal-scatter-gamma→SUPERSEDED (27 tokens green);
  PAPER §4.6 hunt paragraph; plain verdicts: hunt = SUCCESS (mystery
  resolved in locality+driver class), in-sample Einstein = CLOSED
  HONESTLY / NEEDS DIFFERENT DATA.**
  **7D/7E (2026-07-24, the vacuum rung + the platform, pre-reg
  c53f2e2/ff1de37): 7D THE VACUUM-SHARE CONTEST — the classical
  vacuum-free local share q_cl = 1/(2(ν−1)) (pull 2n, no "+1")
  ANNIHILATED on the galaxy ladder: **+556.51 vertical vs AMB (bar
  ≥+25 kill ⇒ ×22 past), cap16 +580.06 = cap-robust LOWER bound,
  plain +274.06, worse than pure BE by +497 — the vacuum +1 is
  LOAD-BEARING (caps q ≤ 1, keeps β finite where the tail data
  demand it)**; binary direct leg instrument-N/A (capped-cliff
  function fails spherical identity 10.61% at NR=2048 → pre-reg'd
  abort; indirect: binaries accept quantum-share AMB, bound β<0.03);
  RJA classical-ambient rung null gal (+0.24, BE→RJ converge
  deep-ambient) / UNRESOLVED bin (+1.68 mean, scatter −2.7…+10.3,
  4 seeds); formal verdict PARTIAL (two-system headline not
  claimable — pre-registered honesty); three gate amendments
  pre-results + one stale-table hazard caught (run-1 table written
  before its gate fired, silently reused — purged). Quantum-ladder
  status: occupation banked / vacuum-share gal-direct+bin-indirect /
  ambient-stats open / fluctuation out-of-sample. 7E THE PLATFORM
  TRANSLATION: the lending gate as a buildable 3D-cQED experiment
  (transmon χ/2π=250 MHz, cavity κ/2π=5 kHz, beam-splitter λ/2π=2
  MHz, thermal n̄ 0.25–8 calibrated): **lending law ½n/(1+n) to
  0.25% at platform params, dissipative shift +1.5%/+1.1%,
  saturation ×7.2, THE L=2 GEOMETRIC RUNG [n/(1+n)]² slope +0.987
  rms 0.0077 (raw-n 34× worse) — the sky's s^L as a bench curve;
  L=2 prefactor 0.480 ≈ ½ = 2nd data point for equal-time-sharing
  ceiling (reading-grade)**; honest frame: mechanism-class
  validation, NOT a gravity test (6N theorem); clean lab negative
  would strike 6X; cQED is THE platform (ions marginal, optomech
  unsuitable); 6N resolution crossover switchable on-chip. Ledger
  rows galfn-qcl + binfn-rja + mech-platform (30 tokens green);
  PREDICTIONS.md §A vacuum-share row; PAPER §2.4 two closing
  passages; NOTES 7D+7E; plain verdicts: 7D galaxy leg SUCCESS /
  two-system NEEDS REFINEMENT / RJA UNRESOLVED; 7E SUCCESS.**
  **README.md ADDED (f5beebd, user-requested): the public front door
  — ELI12 physics, the honesty machine, the open human+AI
  collaboration statement (user wants the repo "open for all" incl.
  the Claude-developer audience); no private/-origin references.
  7F (2026-07-24, pre-reg a4696c1, seam-iv contest): THE MIXING-MEAN
  UNIQUENESS — power-mean family u_p vs the geometric mixing.
  **Asymptotic uniqueness PROVEN (sympy exact, symbolic β: p>0 →
  tail 1 = boot-dead; p<0 → ½ = the measured two-system split
  impossible; ONLY p=0 continuous → ½+g/4 = the measured pair);
  deep ladder p-BLIND through c₃ (symbolic a₄ = −gp/64 + g/192 −
  1/720, da₄/dp = −g/64 — and the p=0 value = 6H's c₄(L=2) EXACT
  cross-validation, two derivations one rung); window bound (y∈
  [3,30] vs measured bands): positive flank excluded at p ≥ +0.1
  (boot direction dies), negative flank bounded at −0.5, both-in-
  band [−0.50, 0.00] ⇒ p_mix ∈ [−0.5, +0.1), p=0 interior both.**
  Formal verdict PARTIAL per the (too-symmetric) pre-reg bar —
  negative flank window-degenerate (snap beyond y~30; coincidence-
  of-scales objection stated soft); seam (iv) upgraded chosen→
  measured; full closure = deeper tails or the c₄ rung. One
  pre-results estimator amendment (far-tail underflow → index from
  u, logged). GF4 regression: window estimator reproduces the known
  postdictions 0.6885/0.5280 at 4 digits. Ledger mech-mixmean (31
  tokens green); PAPER §2.4 passage; plain verdict: theorem+cross-
  validation+positive-flank SUCCESS / negative flank NEEDS
  REFINEMENT (window-limited). **7G EXECUTED (pre-reg 8ce8466): THE
  SATURN CANCELLATION QUANTIFIED — trajectory formulation of AMB:
  Saturn worldline y = 5.38e5 ⇒ u = 1060.8 ⇒ anomaly 10^(−460.7) =
  451 ORDERS under Cassini (B2 PASS; residual = galactic tide
  4026× under; r_M = 7030 AU ⇒ Oort/extreme-aphelion trajectories
  = the residual solar probe; CONTRAST: power-law tails leave
  r-dependent G_eff ~1.9e-6 = ephemeris-dead ⇒ the trajectory door
  is open ONLY for Boltzmann-screened functions = the measured
  screening class — the two program halves point at each other).
  Binary equivalence PARTIAL: mi_t-AMB proxy trails MG-AMB −8.43
  mean (seeds −11.33/−5.54; α̂ 1.55 interior / 2.00 EDGE; G0 α=0
  bit-identity PASS to 0.07) — neither tie (|D|≤5 + interior both)
  nor formal narrows (mean ≤ −10); proxy-crudeness vs genuine
  class-cost UNRESOLVED (α̂ inflation points at discretization
  loss, unproven); the 4L standard-law tie does NOT auto-transfer
  to AMB; refinement = 6 seeds + beyond-proxy implementation
  (unpromised); no post-hoc seed extension (pre-reg fixed 31/101).
  MG-formulation tension STANDS in the world table — a door opened
  by numbers (10^451 one side, −8 lnL the other), not a row
  deleted. Ledger sol-trajmargin + binfn-ambmi (33 tokens green);
  PAPER §8.2 closing passage; NOTES 7G; plain verdicts: Saturn leg
  SUCCESS / binary-equivalence NEEDS REFINEMENT.**
  **7H (2026-07-24, pre-reg e9db7c8): TIE-RESTORED — the 7G gap was
  the proxy's JENSEN DISCRETIZATION. Signed hypothesis committed
  pre-run (mi_t uses B(⟨g⟩); apoapsis-lingering orbits make the
  honest adiabatic functional ⟨B(g)⟩ strictly boostier — +0.078 at
  y=1, e=0.9): all three predictions landed. Six-seed budget (2-seed
  hit the pre-registered ambiguity band → auto-extended):
  **D(mi_avg-AMB − MG-AMB) = −2.92 ± 3.46 SE, sign-mixed 2/6
  positive, interior 6/6, α̂ = 1.272 ± 0.038 (deflated from
  1.55/2.00)**; G0 bit-identity PASS; one pre-results gate amendment
  (G2 quadrature bar 1e-6→1e-5 at e=0.999, 3 orders below likelihood
  resolution). THE COMBINED SATURN STATEMENT (7G+7H, strongest ever):
  the trajectory formulation of the derived function fits the
  binaries AT FIELD GRADE and passes Cassini by 451 orders — the
  quadrupole tension is a FORMULATION property, not a property of
  the measured function; door tail-selective (Boltzmann-screened
  only = the measured class). Honest residuals: α̂ recalibration
  1.27 vs 1.06 between formulations (prescription-level); still an
  MI bracket not a full nonlocal theory; MG-formulation tension
  stands for the field reading; DR4 e-resolution = the formulation
  kill test. Ledger binfn-ambmiavg (34 tokens green); PAPER §8.2
  completion passage; NOTES 7H; plain verdict SUCCESS — maximal-
  demand prong (b) COMPLETE at achievable grade.**
  **7I (2026-07-25, THE EXTERNAL-REVIEW ROUND, pre-reg 2aa6ef9/a17d3a5/
  3ef57c7, PAPER → v2.7): two solicited reviews (LLM + Opus 5.0) adopted
  point-by-point. THE FREEZE (PREDICTIONS.md §0): in-sample function search
  CLOSED at the world-table 13; new forms = consistency rows only;
  primary analyses declared (vertical-hardened hier / six-seed v7corr);
  scatter channel outer/anchored-only, mode-count N retired; ~30-form
  search multiplicity stated uncorrectable. Language: "pre-registered" =
  BAR-LOCKING defined at first use; AMB two-system pass = DESIGN PROPERTY
  not measurement (4 sites); binary ĉ₁ 0.37–0.50 → TWO-REALIZATION
  INDICATION (6 sites); superlatives ledger-scoped; correction counts
  synced (fourteen). Census RELEASED: data/ceiling_pairs.csv (23 rows,
  11/11 reproduced, corrected 9, both conventions per pair; gitignore
  exception). Robustness rows ([calcs/stage7i_ablations.py](calcs/stage7i_ablations.py)):
  **W (w_rad frozen to external Hwang 0.21): CLOSED — Δα̂ ≤ 0.01, the
  eccentricity-nuisance objection dead. S (strict multiplicity, RUWE<1.2
  both + no overluminous, 8047 pairs = 57%): MATERIAL FIRED — fit α̂
  COLLAPSES (Δ −0.740±0.053/−0.457±0.048, Newton +9.4/+12.5 = 16–24% of
  N-scaled, w_rad grid-edge 6/6). Three bar-locked instruments: (1)
  model-light median boost SURVIVES the cut (1.086→1.103 raw, 1.151→1.185
  perp; [calcs/stage7i_material.py](calcs/stage7i_material.py)) =
  companions-as-median-carrier REFUTED; (2) FPM-mismatch hypothesis
  REFUTED on its own bar (strict prefers 1.5 4/4;
  [calcs/stage7i_fpm.py](calcs/stage7i_fpm.py)); (3) α↔w_rad degeneracy
  REFUTED (strict+wr=0.21 → α̂ 0.55/0.65;
  [calcs/stage7i_sw.py](calcs/stage7i_sw.py)) ⇒ COMPANION-DIRECTION per
  the pre-committed tree. Band survives cut at base rate (5/11). THE
  TENSION IS OPEN: model-free statistics see the boost on the cleaned 57%
  where the forward likelihood does not — both cannot be right; credence
  anomaly-real 70% → ~60–65% pending 7J. Pipeline rules: FPM/nuisance
  grids are sample properties (profile them; watch grid-edge riding).
  Ledger rows bin-ablation-wrad/-strict/bin-ceilingpairs (38 tokens
  green); PAPER §6.3 full disclosure + §7.2 CSV citation + §4.6
  end-state-first + freeze/primary declarations (abstract+§1).
  QUEUE (Opus-designed): 7J completeness + MARGINALIZED (f_comp, w_rad,
  f_pm, α) posterior = THE DECISIVE INSTRUMENT; 7K ceiling tail control
  (empirical error tail from γ≥75 anchor pairs); 7L v7 on the Cookson
  selection w/ ablation ladder (after 7J); ĉ₁ six-seed budget; THE SPLIT
  (Paper 1 binaries+ablation+completeness, Paper 2 coefficients+
  quadrupole; §2.4 stays in the log; STRUCTURAL-INDEPENDENCE drafting
  rule — no §2.4 gestures in Paper 1).**
  **7J EXECUTED (2026-07-25, pre-reg d2dc7eb + amendments c893827/a2ecf2e/
  a9e8b16, PAPER → v2.8): THE DECISIVE INSTRUMENT FIRED AGAINST US.
  Part A ([calcs/stage7j_completeness.py](calcs/stage7j_completeness.py)):
  flag completeness C=0.410 (false-flag 0.023), f_photo 0.22–0.30 →
  f_host = 0.42–0.57 (peak 0.51, Raghavan-consistent); fence f≤0.1
  priced at ln π = −768/−608; Banik 0.69 at −36. As-published wobble law
  (q/(1+q), no photocenter cancellation — twins wobble ZERO; correct =
  |q/(1+q) − ℓ/(1+ℓ)|) INTERNALLY INCONSISTENT with measured
  multiplicity: forcing it costs −486…−2010 lnL (why 3K capped at 0.1).
  Part B ([calcs/stage7j_marginal.py](calcs/stage7j_marginal.py), photo
  law + κ_w + fpm→2.4 + host prior; GB0/GB0p exact 8/8): **full sample
  COMPANION-WIN — α_marg = 0.00, dN = +0.0, both laws 2/2 seeds; strict
  COMPANION-CONFIRMED but formally weightless (power gate POWER-FAIL:
  BE arm rode 2.00 edge; simple arm RECOVERED injected α=1.18 at
  1.27/+48.8 = sensitivity proven).** Diagnostics
  ([calcs/stage7j_diag.py](calcs/stage7j_diag.py)): NOT fpm-edge-carried
  (unchanged at fpm=1.5, 8/8); posterior 100% at fcomp=0.35/0.20,
  kw=0.7; α=1 disfavored 38–87 in the marginal; **winning cell itself
  misfits −60…−116 (NO cell fits kinematics AND multiplicity); fcomp=0
  rejected −660…−811 (data demand companions, just fewer than the sky
  hosts)**. CORRECTION #15: 7I instrument-1 "companions-as-median-carrier
  refuted" RETRACTED — under-powered at C=0.41 (cut barely moves the
  companion FRACTION 0.51→0.47); median VALUES stand, immunity inference
  does not; ceiling census leakage null inherits the same audit.
  Quadrant C fired per pre-commit: **credence anomaly-real ~60-65% →
  ~35%** (low end: #15 removed the floor argument; not lower: the −100
  residual misfit + galaxy legs untouched). Binary DATA VOTE SUSPENDED
  in the world table (banner added); all binfn-* contests fence-
  conditional pending re-run under the 7J posterior (7J-d). Ledger:
  bin-7j-completeness/-marginal/-rawcond + ret-7i-median-immune; α/
  Newton/boost rows → CO-QUOTED; audit 44 tokens green. PAPER v2.8:
  §6.3 completion, §7.3 ablation-conditional relabel, abstract/§1/§10
  reframed, App A = FIFTEEN corrections. Forward-modelability caveat
  carried BOTH ways (injection validates against own companion model
  only). NEW DECISIVE QUEUE: **7K-a forward 2C median at the winning
  cell (α=0, fcomp=0.35, kw=0.7) vs observed 1.078 (CI 1.052–1.103) —
  if Newton+companions reproduces it, the median falls; if not, it's
  the surviving anomaly. 7K-b census leakage null at the same cell**
  (folds in the tail-calibration design). Then 7J-d, 7L, ĉ₁ budget,
  THE SPLIT (Paper 1 binary chapter now leads with completeness +
  marginalized null).**
  **7J-x/7J-y (2026-07-25, the Opus review rounds 2-3, PAPER → v2.9,
  App A = EIGHTEEN): round 2 = amendment 6 (per-arm injection standard,
  pre-reg 0045eea: fullpow/fullpowbe/strictpowbe + seeds 202-505;
  standing rules: both-full-arms-fail → PROVISIONAL + ~50%; seed gap
  > −5 or interior α → AMBIGUOUS) + low-end prior stress (verdict
  survives, zero un-pins) + CORRECTION #16 (the "manufactured/
  photometry-forbids" claim retracted hard) + §8.1/§2.4/abstract
  propagations. Round 3 = **CORRECTION #17** (our boundary-free
  criterion applied to the prior axis too: the α_marg(prior-anchor)
  CURVE is the primary statement, labels = annotations) and the big
  one: **7J-y (pre-reg 3f7212e, [calcs/stage7j_paircorr.py](calcs/stage7j_paircorr.py))
  — within-pair overluminosity correlation ρ_core = +0.465 (bar 0.30)
  ⇒ C-FAIL ⇒ CORRECTION #18: Part A's scale RETRACTED (common-mode
  metallicity/age/extinction over-attributed as companions; C = 0.41
  biased LOW, f_host = 0.42-0.57 biased HIGH = the verdict-
  manufacturing direction; the scout's 2-3× literature tension —
  Tokovinin ~0.1/component — resolved in the literature's favor).**
  END STATE: literature-anchored prior (per-pair ~0.22) PRIMARY;
  verdict label = COMPANION-WIN at the retracted measured prior /
  **AMBIGUOUS at the primary anchoring** (seed 31 interior 0.37/0.50,
  seed 101 zero; the in-flight 6-seed budget decides — cubes are
  prior-independent, so the literature-prior marginal comes free on
  completion); DURABLE at every anchoring: **the fenced +99-110
  Newton rejection is DEAD (α_marg ≤ ~0.5, Newton within single
  digits)**. Credence anomaly-real ~35% → **~45%** (the collapse was
  partly manufactured by our own over-scaled prior). #16 one clause
  re-graded (at literature rates our 0.1 fence is NEARER truth than
  Banik's 0.69); #15 → undetermined-pending-repair. fpm note: Opus's
  round-3 premise wrong (grid WAS extended; marginal rides the NEW
  edge 82-99% at 2.4) but diagnosis right = missing variance channel
  (the 6P s-flat object) → add per-system scatter to the marginal
  model (queued). QUEUE: **7J-z = the repaired Part A (joint 2D
  (δ₁,δ₂) mixture: common-mode diagonal vs companion axis-aligned —
  THE gating item)**; fullpowlit (literature-anchored power test,
  AFTER batch — never edit stage7j_marginal.py mid-batch); 7K-b now
  carries the per-pair (q,ℓ) wobble-bound instrument (Opus design:
  max companion-induced ṽ per band pair from per-component δ); the
  paircorr ρ = +0.47 is ALSO a standing astrophysical measurement
  (wide-pair photometric coherence). Ledger: completeness→RETRACTED,
  marginal→CO-QUOTED, +bin-7j-paircorr; audit 49 tokens/106 rows
  green.**
  **ROUNDS 4-5 (2026-07-25, same day): blending audit CLEARS #18
  (ρ(θ) flat 2-250″; blending-safe bar θ≥8″ = +0.448 → C-FAIL STANDS;
  d-gradient astrophysical, with the mass-sampling second-reading
  caveat logged); lit center 0.166 (insensitive); prior = UNMEASURED
  at 0.2-50 kAU (7J-z = first measurement, a contribution); **THE
  CADENCE RULE: credence FROZEN ~45% until 7J-z + 7J-g land;
  measurements not verdicts in interim rounds**; the four-absorber
  diagnosis adopted (boost/near-parabolic/companions/noise share ONE
  width budget; verdicts track the tightest current external
  constraint; only γ separates them) → **7J-g = γ-collapsed four-way
  contest = the Paper-1-deciding test**; arms: fullpow ARM-VALIDATED
  but AT THE RETRACTED CELL (no arm validated at the operative
  literature anchor — cube reuse is seed-axis-only, arms need new
  injected data); fullpowbe ARM-FAIL (BE own-truth 0.73, "both laws"
  DROPPED from verdict phrasing — simple quoted, BE pending);
  7J-w armdiag ([calcs/stage7j_armdiag.py](calcs/stage7j_armdiag.py)):
  own-truth offsets +0.30/−0.40 are PRIOR-FREE (PROF=MARG both arms;
  lit-prior re-marg moves nothing — mismatch does not manufacture α
  at informative likelihood) = likelihood/realization layer at
  single-injection grade; my multiplicative-map clause retracted (an
  additive −0.40 fits the same point; BE α=0.4 injection = the
  discriminator); at the lit anchor the interior fit is NAMED A
  RESULT: α ≈ 0.4-0.5 = HALF the galactic calibration = the
  split-budget signature 7J-g decodes; 7J-z slopes must be PER-PAIR
  (mass-dependent metallicity response from the components' colors;
  amplitude+scatter free — global slope would recreate the
  over-attribution). **ORDER FIXED: 7J-z → 7J-g → arm suite at the
  landed anchor.** PAPER v2.9 updated in place; ledger
  +bin-7j-armdiag (107 rows, 51 tokens green).**
  **ROUND 6 (Opus's last — exchange CLOSED): (1) THE KNEE MEASURED
  OUT (7J-e3): my "no pending instrument un-kills it" RETRACTED —
  but the fine curve (anchors 0.06–0.30, both σ) is nearly FLAT:
  α_marg 0.18–0.31, Newton +2…+4 everywhere incl. below the old
  fence, because the likelihood alone pays +12…+28 to sit at
  fcomp=0.2; **the fenced +99–110 returns ONLY behind a hard
  exclusion of f≥0.2 — no rate measurement revives it; the revival
  mechanism is the MISSING WIDTH CHANNEL** (= the −60…−116 misfit =
  the fpm edge = 6P's s-flat scatter, all one object) → the
  variance repair folds into 7J-z's modeling round; grid caveat: no
  cube cells in (0.1, 0.2), denser cells queued. (2) INJECTION
  OPTIMISM QUANTIFIED (7J-w2): injections 4–30× more informative
  about fcomp than the sky (real: tens across 0.1–0.35; injected:
  cliffs of hundreds) — post-7J-z arm suite must degrade injections
  with the measured common-mode. Ledger +bin-7j-anchorcurve (108
  rows, 54 tokens green). 7J-s6 closed the batch (c06ecb2): NO seed
  break 16/16; lit-anchor 0.13±0.07/0.30±0.09, Newton +1/+4.**
  **THE TWO DECIDERS EXECUTED (2026-07-26, 7J-z/7J-z2/7J-z2b/7J-g,
  PAPER → v3.0, review rounds 7-8 absorbed in-flight): pre-reg chain
  f991372→96a72ad→ffc0195→e6bf62e→08d2ec5 (every bar before every
  run). 7J-z v1 ([calcs/stage7jz_mixture.py](calcs/stage7jz_mixture.py)):
  estimator gates ALL ~exact (GZ1a/b = 0.000 — the #18 mechanism dead
  even under mis-spec; GZ2 0.103/0.251), f̂ = 0.166 blended/component
  — GZ6 postdiction FAILED informatively → THE VECTOR COMMON MODE
  (ρ(dcol) decays 0.86→0.36 smoothly; twins project identically;
  0.084 mag excess coherence at dcol<0.1) — no anchor shipped, npz
  quarantined. v2 (kernel, GH quadrature) KILLED by its own GZQ
  (d=18.7 — twin-spike integrand; + start-gap 43) mid-run; v2b
  ([calcs/stage7jz2b_exact.py](calcs/stage7jz2b_exact.py), EXACT BVN
  likelihood): +476 over v1, GZ6′ FULL PASS incl. twin slice,
  **GZ1c-b ACQUITTAL (kernel-blind fit on kernel-truth sky = 0.000)**,
  f̂ chain 0.166→0.162→0.159 (3 model classes, 2 integrators, 5
  basins), **host 0.29-0.32/component = ~3× field (prices lit 0.09 at
  −1634; first measurement at 0.2-50 kAU, scout-confirmed unmeasured
  range)** — UNSHIPPED on two numerical-protocol bars (GZS11 2.08,
  start-gap 2.8; no-third-iteration rule held) → **v2c queued =
  certificate job (profiled-λ + resolution-gated continuous
  q-integral)**. THE WIDTH CHANNEL (amendment 7 photow mode,
  [calcs/stage7j_marginal.py](calcs/stage7j_marginal.py); GB0w caught
  a wobble-law wiring bug first firing, fixed, 4/4 retests 0.00e+00;
  8 cubes + γ-collapsed twins committed): sq = 0.2 INTERIOR demanded
  (P(sq>0)=1.00, zero edge mass, 4/4) = the 3E/6P object
  parameterized. **DECIDER 1 ([calcs/stage7jz_read.py](calcs/stage7jz_read.py)):
  α_marg = 0.74/0.70 at ΔlnL +23.8/+23.2, ANCHOR CURVE FLAT to
  0.01-0.03 across 0.06-0.34 (α companion-prior-INDEPENDENT; at the
  unshipped measured anchor: 0.67-0.74 at +12..+16); the old
  COMPANION-WIN lives in the sq=0 slice; D2 INVERTED (multiplicity
  cost ×1.4 — with sq the kinematics want fcomp=0.1; tension vs
  measured ~0.3 host = named object); D3 BE rides fpm=2.4 at 0.97
  (extension to 3.0 queued); VERDICT AMBIGUOUS-CARRIED
  (LIT-CONDITIONAL) — missed the pre-registered α≥0.5 AND dN≥+25 bar
  by 1.2-1.8 lnL. DECIDER 2 ([calcs/stage7jg_read.py](calcs/stage7jg_read.py),
  dual-config amendment 8): PRIMARY (sq free) = ABSORBER-LEVEL (γ
  pins SD(w_rad) 0.015-0.035→0.000; α channel-robust 0.70-0.74 vs
  vt 0.62-0.66); CO-READ (sq=0, four-absorber) = SEPARATION-CONFIRMED
  α grade — **ṽ-only reports a PHANTOM α≈0.50/0.52 at +12.8/+12.2
  that the γ data veto to 0.20/0.25 at +2.4/+3.7** (the width-less
  configuration = every published ṽ-only pipeline; reading). Tier
  difference = the finding (γ and sq are partially interchangeable
  absorber-police). PAPER-1 RULE: PROCEEDS (2D claim operates;
  phantom-veto = the methods centerpiece). CREDENCE MAP APPLIED:
  AMBIGUOUS → **anomaly-real HELD ~45% by the map** (freeze released,
  value re-affirmed; miss direction noted). Ledger +bin-7jz-v1/-v2b/
  -width/+bin-7jg-gamma, anchorcurve→CO-QUOTED; audit 113 rows/64
  tokens green; worldtable banner = width-channel form. NEXT:
  **v2c-plus DONE + anchored re-read DONE (see the certificate-round
  block below) → the ARM SUITE (pre-reg f6916e0; implementation =
  the immediate next work block; E-arm fpm→3.0 first, then the A–D
  twin-mismatch injection arms)**;
  then 7K-a/b, 7J-d, 7L, ĉ₁ budget, THE SPLIT
  (Paper 1 leads: completeness+kernel discovery+width channel+
  phantom-veto+census). Review rounds 7-8 (Opus): axis chain written
  (P/H/PAIR; the 1−(1−0.24)²=0.42 trap named), his pair-map +
  fence/anchor claims withdrawn by him, dual-config + anchor-curve
  re-read adopted from him.**
  **ROUND 9 (2026-07-26, same day, PAPER → v3.1): the adequacy
  decomposition (7J-z3, [calcs/stage7j_sqclose.py](calcs/stage7j_sqclose.py),
  read-only, frame+bands pre-committed 5024ecb, G0/G1 regressions vs
  the shipped read exact 4/4): **S1 anchor strain = 0.0 in 4/4 —
  α_marg 0.74/0.70 IS the width-complete model's own free optimum
  (the LIT16 anchor does no work); S3 FORCED fcomp≥0.35 → α̂ = 0.00,
  dN +0.0 in 4/4 = FIFTH-MOVE-LIVE — α is FULLY exposed to the
  multiplicity tension** (if the ~0.3 host rate certificates AND the
  per-companion signature is as modeled, α collapses — at a cell that
  still misfits by 135–153; escapes: v2c certificate fails, or the
  kernel direction — common-mode masquerading as companion light —
  which the DEGRADED-injection arm suite decides; the D2 NOT-CLOSED
  token had already shipped — Opus's "unreported" premise wrong, his
  question right). Also booked: sq IDENTITY = open flag in §6.3
  (existence measured 4× — 3E/6P/misfit/fpm-edge — but no traced
  channel makes 20%; size concession: sq=0.2 > the ~12% signal;
  median 1.078 sq-IMMUNE); structural argument made explicit
  (γ-invisible absorber soaks γ-invisible systematics, α identified
  by direction shape; w_rad externally anchored, 7I-W Δα̂≤0.01);
  grid-step phrasing corrections (WR/SQ SD "0.000" → single-node,
  steps 0.10/0.1; cosmetic SD(sq)=0.927 = sqrt(0.86) pinned-axis
  broadcasting artifact, explained); Paper-1 reframe ADOPTED (spine =
  width channel + phantom-veto + absorber accounting as the field's
  mechanistic reconciliation, stated as READING; α conditional on its
  flat curve; Cookson median-flatness NOT explained by symmetric
  width — stays independent, 7L unchanged); "3× field from triple
  dynamics" NOT printed pending scout. No bar moved, credence HELD
  ~45% (cadence: v2c + arm suite are the deciders). Ledger
  +bin-7j-sqclose (114 rows/66 tokens green).**
  **ROUND 10 (2026-07-26, Opus's exit — he owns the D2 misclaim, 3rd
  false execution claim in 4 rounds, ends round-by-round audit; door
  open for 7L/split reads; PAPER → v3.2): 7J-z4
  ([calcs/stage7j_qmoments.py](calcs/stage7j_qmoments.py), pre-framed
  3cbedbe, amendments A1/A2 logged post-fail pre-quote at 2a3821d,
  amended run ALL GATES PASS): **(A) anchor-STRENGTH curve — his knee
  mechanics wrong (collapse axis is width, not center-range: σ=0.03
  anchors hold α even centered 0.35+), substance right and MEASURED:
  σ* = 0.02 both laws (α holds 0.67–0.75 at σ≥0.03; collapse routes
  via the fcomp=0.2 cell; 0.35 only at σ=0.005); unshipped-rate
  precision ~0.015 < σ* ⇒ a face-precision certificate FIRES the
  fifth move on the flat-q axis. (C) kw attribution: the 135–153
  forced-multiplicity rejection is WOBBLE-BINDING (Dwob +314/+306 —
  kw 0.7→1.4 triples the cost). (B) his two-moments objection
  VINDICATED-and-completed: detected companions carry 0.13–0.17 of
  undetected wobble variance (6–8×; hard-flag q_min 0.70–0.88; wfac
  peaks q≈0.4–0.6); his gap = the hidden-mass channel √(1+q/2),
  twin-maximal, π-stable (0.25–0.30), non-binding; JOINT conversion
  band fce_joint(wobble) = [0.10, 0.39] — the detection-shaped
  bracket makes host 0.30 ≡ 0.10–0.12 kinematically (his "~0.12"
  reproduced); formal letter = MIXED (A2: first band missed the
  completeness rescaling — logged, joint reported alongside).
  CONSEQUENCE: the exposure lives on (rate, precision, π_q) JOINTLY —
  the ×3 multiplicity tension is NOT YET KINEMATIC; **v2c UPGRADED
  to v2c-plus BY REQUIREMENT: ship a (q,P)-resolved rate + extract
  the v2b posterior's q-information, never a scalar**; his
  "indifferent at +0.0" corrected (that's α-indifference WITHIN the
  forced world; the likelihood pays 135–153 to avoid it);
  "tension-as-result" adopted for Paper 1's spine. Ledger
  +bin-7j-qmoments (115 rows/68 tokens); credence HELD ~45%
  (cadence).**
  **THE CERTIFICATE ROUND (7J-z2c + anchored re-read, 2026-07-26
  night, PAPER → v3.3): v2c-plus SHIPPED, ALL GATES
  ([calcs/stage7jz2c_cert.py](calcs/stage7jz2c_cert.py), 135.8 min;
  pre-reg 4c3fc89; amendments A3/A4/A5 = gate-metric designs each
  caught by its own STOP, logged pre-quote): both v2b protocol
  objects closed STRUCTURALLY — segment-erf exact q-integration
  (certified 1.7e-8 vs 30001-node reference; f=0 identity exact) +
  profiled-λ (start gap 0.03 vs v2b's failing 2.8; λ*=1.64 interior;
  GV1 outer-doubling 0.010 vs GZS11's 2.08). **f̂ chain 0.166 →
  0.159 → 0.159 six-way stable; host [0.29, 0.32] flat/qalt; the
  acquittal repeats (kernel-truth f=0 → 0.000); npz v2c-cert carries
  the GV7 q-table + conv_band + face-invalidity note. THE GV7
  DISCOVERY: the photometry measures the subsystems TWIN-HEAVY —
  twin t=5 beats flat +162.0 (q^-0.5 −267; t=2 −19.9 behind t=5),
  f̂ = 0.124 → host ≈ 0.23 ≈ field-adjacent; twins = the
  wobble-quiet population ⇒ twin-convention host ≡ flat-q kinematic
  ~0.10–0.15 = the kinematic preference — THE ×3 TENSION DISSOLVES
  in round-10's predicted direction** (reading-grade; joint fit =
  future). Scout (pre-committed with the round): Opus's "3×
  expected from triple dynamics" UNSUPPORTED — Tokovinin 2014 finds
  field-like wide-binary subsystem rates ("does not affect
  inner-subsystem formation"); the 3× literature = contact-binary/
  2+2-specific; El-Badry-Rix ~20%/component flagged for primary
  read. ANCHORED RE-READ (amendment 11 = LANDED-CONV operative,
  d437921 pre-committed; 11b = OPER wiring fix c272bef — first
  anchored run's verdict block read the face anchor, caught in
  output, fixed to spec, re-run): **OPERATIVE α_marg = 0.74/0.70 @
  +23.8/+23.2 = AMBIGUOUS-CARRIED at the LANDED anchor** (no longer
  LIT-conditional; extension rule silent, BE gap 0.10; = LIT16 to
  the last digit — flat curve); **FACE (flat-q face-precision
  diagnostic): THE FIFTH MOVE FIRED LIVE (α → 0.00, fcomp → 0.2,
  exactly per σ* = 0.02)** — quoted as the conditional whose premise
  the q-table rejects at −162; 7J-g tiers unchanged at LANDED-CONV;
  D3: BE rides fpm 2.4 at 0.95–0.98 (E-arm). Ledger +bin-7jz-v2c/
  -qshape/-anchored, width→CO-QUOTED (118 rows/73 tokens); credence
  HELD ~45% (cadence: the ARM SUITE = the remaining decider).
  ARM SUITE pre-reg committed f6916e0 (bars locked pre-anchor):
  A null ≤0.3 / B simple 0.74±0.25 / C BE 0.70±0.25 / D BE 0.40
  discriminator / E sky fpm→3.0 + identity regression; injections
  q-MISMATCHED (draw twin t=5 = the GV7 winner, fit flat-q) + sq_true
  0.2 + fpm_true at posterior mode. **E-ARM EXECUTED (7J-z5-E,
  [calcs/stage7jz5_eread.py](calcs/stage7jz5_eread.py), FPME wiring
  fca062e, GB0w+GB0e 0.00e+00 ×4, GE0 exact 4/4): THE NOISE EDGE
  CHASES — P(fpm=3.0) = 0.54/0.97, mode at the new edge both laws,
  α 0.74→0.68 / 0.70→0.67, dN +23.8→+14.5 / +23.2→+16.1 (~8 lnL
  ceded); Dwob stable +305/+310. fpm=3.0 = 3× Gaia formal errors =
  past the Lindegren+21 ceiling (~1.4) ⇒ the missing-width-SHAPE
  signature (the 3E→6P→sq ladder's next rung); pre-registered
  decision clause: NO further extension (unphysical direction);
  operative α CO-QUOTED with the band 0.68–0.74 / +14.5–23.8
  (ledger bin-7jz5-earm, 119 rows/75 tokens).
  **THE ARM SUITE (7J-z5 A–D, 2026-07-27, PAPER → v3.4): 4/4
  VALIDATED ([calcs/stage7jz5_armread.py](calcs/stage7jz5_armread.py),
  implementation bad6770, bars f6916e0, chase bands pre-stated in
  the reader commit before the batch; ~91 min GPU): **ARM A (null,
  fcomp 0.20 twin-mismatch): α_marg = 0.00 BOTH laws — the
  machinery manufactures NOTHING (7J-w2 answered at the operative
  anchor); B 0.74→0.65/+20.8, C 0.70→0.64/+38.7, D 0.40→0.48
  PROF=MARG — all within bars, no extension; nuisances sharp 8/8
  (sq at injected 0.2 everywhere). CHASE-UNEXPLAINED: P(fpm=3.0) ≤
  0.02 in 0/4 arms vs sky 0.54/0.97 ⇒ the sky's noise hunger =
  REAL width-SHAPE incompleteness (not companions/twins/sq) —
  width-shape refinement = the named next binary instrument;
  matched-truth arms score +21…+39 vs the sky's +23 (residual
  width eats sky significance, coherent). Own-truth biases mildly
  conservative (−0.06…−0.09, single-injection grade, not applied).
  CREDENCE: both cadence deciders landed → anomaly-real ~45% →
  ~50% (reasoned move, disclosed: no mechanical map existed; up on
  the companion-alternative collapse + machinery validation,
  capped by AMBIGUOUS-CARRIED + the open width-shape systematic).
  Ledger +bin-7jz5-arms (120 rows/77 tokens); worldtable banner =
  validated form.
  **7J-z6 (2026-07-27, THE WIDTH-SHAPE CONTEST, pre-reg 66a8045 +
  amendments 2/3 pre-stated ec7514c/ab0c399, PAPER → v3.5): verdict
  UNRESOLVED-CARRIED (locked grammar) — floor B1 one-law
  (+7.12/+10.61, edge-rides 0.045 km/s WITH fpm 3.0 = 4.3× formal),
  tail DEAD (+0.0, axis unused); band 0.68–0.74/+14.5–23.8 STANDS.
  THE ROUND'S PRODUCT = the NOISE-CEILING PROFILE ("ship the risk
  axis" executed): PHYS envelope (fpm≤1.8 Lindegren, ws≤0.015) **α
  = 0.80/0.80 at dN +35.2/+32.3 — the operative quote is
  noise-DILUTED** → operative band → unphysical corner α = 0.00
  (the floor-leg B3 MATERIAL-LOW). EXCHANGEABILITY ARM INFORMATIVE
  (pre-stated bar): injected-0.74 sky under the floor fitter KEEPS
  its boost (α 0.60, floor quiet, P(fpm=3.0)=0.00) ⇒ the sky's
  corner preference is DATA-DRIVEN (~7–11 lnL genuinely unmodeled
  width). VOLUME MECHANISM named: the tail axis unused at peak
  still drains dN to +10.1/+6.2 (sub-B1 axes zero α_marg by
  posterior volume — transferable caution). Shape EXCLUSION SET:
  floor/tail/fpm/sq/perspective (perspective closed by composition:
  100% RV coverage); Part-A fingerprint (G0-diag exact): mid-
  shoulder ṽ 0.7–1.7 + radial γ=8° column + inner bins ⇒ successor
  hypothesis BOOKED not opened: the inner-bin eccentricity/radial-
  population sector (eta grid = 2 values, e-sector nearly frozen).
  PHYS = reported conditional, non-operative (promotion = a
  decider-grade question, flagged for the next external review).
  All gates exact (GW0 0.00e+00 ×4; ws==0 branch verbatim-legacy).
  NO credence move (pre-stated). Ledger +bin-7jz6-widthshape/-xchg
  (122 rows/81 tokens). E&R 2018 PRIMARY READ done (18bb2ca):
  §2.4.1 verified — 10.5% measured (q≳0.5, twin-adjacent) / ~20%
  per component extrapolated / ~36% per pair ⇒ brackets the v2c
  host rate; the ×3-vs-literature tension = wrong-comparator,
  externally confirmed.
  **THE ACCOUNTING ROUND (7K-a/7K-b/7J-d, 2026-07-27, PAPER →
  v3.6): 7K-a forward median (pre-reg 3c931d5 + two pre-quote gate
  amendments — G0-i/iii FAILs = MEASURED model facts: the fitted
  e-run makes model-Newton's 2C ratio ≈0.98 ("Newton predicts
  1.000" is scale-free-conditional) + a_s-prior edge truncation;
  G0-iv interior control PASS 1.005/1.002): **GRAY both seeds —
  the landed Newton-best cell forward-produces R = 1.033/1.043 of
  the 1.078 anchor (CI 1.052–1.103): absorbers cover 52–59% of the
  model-referenced gap, unabsorbed excess +0.035–0.045; boost
  cells OVER-produce (1.114–1.127); the median = model-light
  EXCESS, no longer a standalone Newton-killer** (anchor rows
  note-updated, no flips). 7K-b census null (pre-reg eb6084d;
  post-hoc legs LABELED): **NULL-BROKEN by the letter (landed cell
  floods the band 26–38 vs 9) and SELF-REFUTED by the cliff (17–32
  predicted vs 2 observed above 1.67, P ≤ 7e-6; no config incl.
  sq=0/PHYS/boost reproduces 9-with-cliff-2, min P ≈ 9e-5) ⇒ THE
  RELOCATION: the (band=9, cliff=2) PAIR is the operative census
  statistic, self-defending (cliff bounds tail noise to ≈formal,
  where the band is unleakable at 4J's 3.8e-9); count retired
  ABSORBER-CONDITIONAL; third mid-shoulder-localization witness;
  wobble-lite bound logged VACUOUS (ignores catalog survival);
  7K-b2 folded into the population successor.** 7J-d THE
  UNSUSPENSION (pre-reg fdb4bfc; G0a0 exact ×16, GD1 lam100≡BE
  cube BIT-EXACT both seeds; 16 fn × seed 31 + 12 × 101, ~4h GPU
  each): **NO FUNCTION DISCRIMINATION at the landed anchor — all
  16 within ±8 lnL of BE, seed scatter comparable; four fenced
  vetoes (rb4/boot/resn/dwf) FLIP per the letter = instrument
  resolution loss, not rehabilitation (fenced vetoes stand as
  fenced-model results); simple's fenced −12.6 lean over BE →
  no-verdict (3P generalized to all 16); λ profile FLAT (peak c₁
  0.375 seed 31 / 0.000 seed 101) ⇒ the 4X binary ĉ₁ = 0.37–0.50
  ERASED under the landed posterior — bin-c1 CO-QUOTED, the §6.4
  "two systems read the dial" claim now GALAXY-LED (supersedes the
  six-seed fenced budget by construction); NEWTON FUNCTION-ROBUSTLY
  DEAD (dN ≥ +7.9 on all 32 fn×seed rows; seed-mean +14.5/+16.1 =
  the co-quoted band); P(fpm=3.0) ≈ 1.0 function-blind = the width
  object is upstream of the ν-choice; SPARC = the program's only
  function discriminator (4E generalized); sbe/qcl excluded-
  disclosed, trigger not fired.** Ledger +bin-7ka-median/
  +bin-7kb-census/+bin-7jd-unsuspension, bin-c1→CO-QUOTED (126
  rows/88 tokens green); worldtable banner = unsuspension form;
  PAPER v3.6 (header, abstract-(3) rewrite, §6.3 7K-a passage,
  §6.4 erasure passage, §7.2 7K-b passage, spine/§2.3/§4.3/
  scorecard c₁ touches).
  **7L (2026-07-27, THE COOKSON QUANTIFICATION, pre-reg ad1e79e,
  PAPER → v3.7): proxy mask N = 1194 vs their 1421 (16%, inside the
  2× bar; |ΔRV|<10 halves both-RV; s∈1–30+130pc removes the deep
  anchor). **L3 (their statistic, G0 scale-free PASS 0.995/1.002,
  e-run baseline ≈1.00 on THIS statistic): their flatness
  REPLICATES on our catalog (data step 0.985, CI 0.909–1.060);
  the landed boost cells predict 1.092–1.101 = HALF their tested
  20% step ⇒ their 2.7σ ≈ 1.3–1.4σ vs the landed model;
  boost-vs-Newton ≈ 0.9σ_boot at their N; verdict GRAY by the
  letter. L1: CONSISTENT 4/4 (d_op 1.2–3.3 ≤ 4 — the cook sample
  cannot reject the operative α; α=2 excluded −25…−39); marginal
  dN +0.0–0.1 = the N-scaling expectation. THE FORK IS ARITHMETIC:
  their cleaning discards 92% of the joint sample incl. the whole
  deep anchor; the remainder carries ~1σ of α-information — their
  null and our +14.5–23.8 are COMPATIBLE; the disagreement is
  entirely whether the removed 92% is data or contamination.**
  Ledger +bin-7l-cookson (127 rows); PAPER §7.4(d) new; worldtable
  banner + tokens.
  **ROUNDS 11–12 (2026-07-28, the Opus structural read + the
  round-12 arc, PAPER → v3.9): round 11 adopted 1–6 with three
  corrections (A: neither noise-curve endpoint is a candidate truth
  unidentified — brackets; B: "0.0–0.8 each defensible" corrected;
  C: the +0.0-indifference misread repeated + tension staleness);
  curve-as-primary, conditional load-path abstract, named
  multiplicity-interface section, Paper 2 spine REPLACED (two-dial
  device retired; 16-law ±8 degeneracy = field-level methods
  finding), 7K-b flagship + counts/search rules, 7L = well-posed
  question (TODO 17/18). Round 12: he owns C; B-refinement adopted
  (rejection basis IN the abstract); his solver question → **7J-z7
  (pre-reg a03f3d0, G0-q EXACT ×4 stream-preserving): FIFTH-MOVE-
  ALIVE (MATERIAL 4/4) — the twin world trades wobble for HIDDEN
  MASS (cost halves 135–153→69–87); his premise corrected (S3 ran
  on the corrected law), his instinct vindicated; round-11
  rejected-premise defense self-corrected (q-table PRICES, does
  not exclude); DIVIDEND: free twin-q marginal 0.61–0.84 =
  q-law-ROBUST. 7J-z8 (pre-reg 285045a): EXPOSURE-CONTAINED — the
  cliff rejects the twin-forced world at 5.5e-14 (overshoot 37 vs
  2) in both laws; the MEDIAN REPRODUCES it (1.085) and RETIRES
  from this defense (free-world 7K-a role unchanged); both
  expectations missed, both logged. FINAL EXPOSURE FORM: the
  (band=9, cliff=2) pair holds BOTH collapse flanks — the load-path
  abstract clause is final.** Ledger +bin-7jz7-twinforce/
  +bin-7jz8-contained (129 rows); split plan final abstract form +
  decision-4 recommendations (coherence vector → standalone note;
  rate stays abstract-as-result). NEXT: THE SPLIT prose assembly
  awaits the author (titles, credence placement, author line,
  spin-off notes) — tell the user before any Opus-addressed
  draft.** Then the DM
  mimicry-ledger paper
  section.** Paper
  leads with: 1.086 boost, s-dependence, α interior in every model class (span 1.0–1.5,
  Newton +55…+112 always), the w_rad=0.20 ↔ Hwang+22 superthermal cross-validation
  (radial excess = confirmation, not discovery — 4G), triples-exoneration, realization
  systematic, BE identity (C&T's; we test).

## Discipline (this project's identity — keep it)
- Every result gets a validation gate before trust; every claimed number has a script.
- Record-keeping (user directive 2026-08-05): [LOG.md](LOG.md) is the plain program
  log — ONE short line per stage ("did -> found"), added in the same commit as the
  stage's verdict. Keep NOTES entries lean; no verbosity. PAPER.md stays FROZEN
  between author-called rewrites — stages are absorbed into the paper only when the
  author calls a rewrite, never per stage.
- Every new pipeline MODE gets a regression against the old mode at its identity
  point (the GB0w precedent — it caught a silent amplitude-law fallthrough on first
  firing; review-adopted corollary: gates only catch what someone thought to write).
- Never edit PAPER.md (or any Unicode-heavy file) via PowerShell Get/Set-Content —
  PS 5.1 double-encodes; one full-file corruption caught by git diff --stat
  (2026-07-26). Use the Edit tool.
- A flatness/robustness claim is only as strong as the reason the scanned range
  or axis ended where it did — ship the risk axis with the curve (round-10
  standing note: the anchor curve was flat in CENTER while the collapse lived on
  WIDTH at σ* = 0.02, one measurement away; scalar conversions between
  differently-weighted moments of a distribution are the same trap — pass the
  distribution, per the q-moments audit).
- Reader/lift scripts copy model expressions BIT-VERBATIM including the
  arithmetic FORM of constants (0.995−erf ≠ literal 0.045 by 4e-17; the
  near-parabolic integrator amplifies ulp-level input differences to ~10 lnL
  at floor-amplified likelihood cells while census stays half-count-blind —
  8P/8Q amendment 1, caught by the program's first lnL-grade identity gate).
  Wherever a stored cube exists, gate reader identity at lnL grade, not
  census grade.
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
