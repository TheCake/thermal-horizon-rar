# TODO / Next Steps

*Priority-ordered. Each item names its blocking question and the script it extends.*

## Hardening the hierarchical fit (the α number)
1. ~~**v3 fit: contamination + priors**~~ — DONE ([calcs/stage3c_v3fit.py](calcs/stage3c_v3fit.py),
   Stage 3C in NOTES). Outcome: hypothesis FALSIFIED — ridge did NOT collapse; α still
   corner-seeks even with η pinned at 1.3; contamination worth only +8/+12 lnL. Newton
   rejection robust (+296/+264). No α headline exists; #2 is now the blocking item.
2. ~~**Control-bin autopsy**~~ — DONE (Stages 3D+3E in NOTES,
   [calcs/stage3d_autopsy.py](calcs/stage3d_autopsy.py), [calcs/stage3e_smear.py](calcs/stage3e_smear.py)).
   Verdict: data ṽ distribution is broader than every model in every bin; a multiplicative
   smear σ_m ≈ 0.2–0.25 localizes interior, collapses the α corner-seeking (α=1 beats α=2),
   and deflates the Newton rejection to ≈ +63/+66 lnL. v4 fit with the σ_m axis =
   [calcs/stage3f_v4fit.py](calcs/stage3f_v4fit.py).
2b. ~~**Decompose σ_m physically**~~ — RESOLVED (Stage 3N/v6b): the broadening is a
   **~20% near-parabolic (e>0.9) sub-population** (w_rad=0.20 interior for both laws,
   +129/+141 lnL) plus small localized contaminants (all fences now interior); f_pm minor.
   Full refutation chain: masses (3J) → companions (3K) → circular (3M) → radial excess
   CONFIRMED (3N). Remaining caveat: α̂ is model-dependent across broadening identities
   (1.0 with σ_m → 1.5 with w_rad); report the span.
2c. ~~**v7: emulate the catalog's velocity-consistency selection**~~ — DONE (Stage 3O,
   [calcs/stage3o_v7fit.py](calcs/stage3o_v7fit.py)): envelope measured from data first
   (physical bound ceiling v·√s≈2.2 well inside the 5 M⊙ catalog cut ≈3.0); acceptance
   applied to model + templates; tail residuals vanish (+376 lnL); **α INVARIANT
   (interior 1.5, Newton +98/+105); w_rad strengthens; BE leads simple 2nd time (7.5)**.
2d. ~~**v7 error budget**~~ — DONE (3P seeds + 3Q bootstrap). **FINAL v7 numbers:
   simple α = 1.52 ± 0.46; BE α = 1.54 ± 0.13; Newton loses in all 1000 replicates
   (min +38/+59).** Fit program v1→v7 complete.
2f. ~~**g_ext scan**~~ — DONE (Stage 3S). **Simple law: α=1 recovered sharply at
   g_ext = 1.4–1.6a₀ with the scan's best likelihoods — the α>1 tension is an
   EFE-calibration degeneracy, not a data demand. BE law: NOT rescued at any scanned
   g_ext (α̂ 1.5–2.0 everywhere) — first α-structure separation between the ν-families.**
2g. ~~**Pin down the MW's Newtonian external field**~~ — DONE (Stage 3T). RAR inversion:
   g_N,ext = 1.15±0.05a₀. **At the physical field BOTH laws localize at the
   parameter-free α=1 (α̂ 1.17/1.11, interior, Newton +109/+99).** Root cause of the
   old 1.9: Stage 2G imported the AQUAL-total convention into our QUMOND solver —
   OUR cross-formulation bug, owned in NOTES; published Chae (AQUAL+total) and Banik
   are each internally consistent.
2h. **Re-run the 3P/3Q error budgets at g_N,ext=1.2** ← final mechanical gate — 6 seeds +
   1000-replicate bootstrap on the corrected tables; produces the paper's headline
   α ± (expected ≈1.1±0.15). Also verify Chae 2023's exact external-field sentence
   (full-text read, not scout) before writing the comparison paragraph.
2e. **The γ≈82° perpendicular excess** (ṽ≈0.07 and ≈1.66 cells, wide bins, z≈+5, present
   in every model) — the last unexplained structure. Candidates: genuinely circular
   sub-population entering only at wide s (formation?), resolved-triple outer pairs,
   or a projection/selection artifact. Standalone puzzle; not blocking the paper if
   disclosed.
3. ~~**Error budget on lnL**~~ — DONE (Stages 3H+3I). Realization (6 seeds): simple
   0.93±0.11, BE 1.30±0.15, interior 12/12, Newton +55±4. Bootstrap (1000 replicates):
   simple 0.98±0.20, BE 1.21±0.26, Newton +60±11 (min +30). **Combined: simple
   α = 0.98 ± 0.23; BE α = 1.21 ± 0.30** — conditional on the σ_m nuisance being benign
   (see #2b/v5).
4. **Mass-model systematic** — marginalize M/L in the Stage-1 screening-index fit and
   propagate photometric-mass errors into ṽ (affects both p and α).
5. **Chance-alignment stress test** — scan R_chance thresholds (0.01 → 0.001) in the
   20–50 kAU bin (N=214); check boost stability.

## Reconciliation (the paper's credibility keystone)
6. **Reproduce the Banik-style statistic** (binomial pixel likelihood, their cuts) on OUR
   sample, and our statistic on their cuts — localize exactly where the 10σ-vs-16σ
   disagreement is created. The realization-systematic measurement (Stage 3A) is the frame.

## Theory
7. **Literature hardening** — citation-graph walk of McGaugh+16 (~1000 citing papers) for
   any prior statement of the Bose–Einstein occupation identity; deepen the Pazy & Argaman
   comparison (their statistics-on-screens vs our occupation reading).
8. **Two-field BE theory** — a proper Lagrangian whose EFE structure follows from the
   occupation postulate (the sphericalized QUMOND solve is a stand-in); predict the
   anisotropic (curl) component; feed back through the orbit engine.
9. ~~**NLO test**~~ — DONE (Stages 4A/4B). Truncated-expansion estimator is power-limited
   on SPARC (c₁ ± 0.4–0.6, honest null); the truncation-free full-function branch
   comparison is decisive: **c₁=½ branch (BE/simple) beats c₁=0 branch (standard-μ) in
   198–200/200 galaxy bootstraps — the BE identity's parameter-free NLO prediction
   PASSED its kill test**; within-branch (NNLO 1/12 vs 1/8) unresolved, slight simple
   lean. a₀=1.206e-10, f_ML=1.10 recovered free. Remaining BE kill test: a₀∝H(z) (#14).

## Publication path
10. **Paper assembly** from [PAPER-DRAFT.md](PAPER-DRAFT.md) once items 1–5 land.
11. **External review** — hand the repo to a professional (ask them to break Stages 1, 2C,
    3B specifically); they are also the arXiv endorsement path.
12. **Zenodo DOI snapshot** for a citable timestamp; arXiv preprint after review.

## Longer horizon (data that decides things)
13. Gaia DR4: per-mille orientation statistic; α at distribution level with 10× pairs.
14. High-z rotation curves (JWST/ALMA): a₀ ∝ H(z) vs constant — the BE reading's kill test.
15. Wide-binary eccentricity distribution at 10–50 kAU (the η=edge hint): measurable with
    Hwang's v-angle method on the EDR3 catalog — possible standalone paper.
