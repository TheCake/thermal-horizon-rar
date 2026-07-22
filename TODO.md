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
2b. **Decompose σ_m physically** ← **NEW BLOCKING** — 0.2–0.25 implies ~50% effective mass
   error; the accounted photometric budget is far smaller. Candidates: M_G→mass table
   scatter (age/metallicity/extinction), unresolved binarity beyond f_t, e-model mismatch,
   intra-bin selection. Measure it independently (main-sequence width of the sample;
   FLAME/spectroscopic masses for a subsample) and turn σ_m into a PRIOR, not a free float —
   the α interval inherits the (α,σ_m) degeneracy until then.
3. **Error budget on lnL** ← **gates the α interval** — MC repeats of the v4 fit
   ([calcs/stage3f_v4fit.py](calcs/stage3f_v4fit.py)) over seeds + data bootstrap; the v4
   near-peak profile differences (1–3 lnL) are within single-realization noise, so the
   [0.87, 1.05]-style intervals are NOT quotable until this lands. Finer α grid near 1.0.
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
9. **NLO test** — the identity predicts the next-to-leading RAR coefficient = 1/2 exactly;
   test on SPARC with M/L marginalized (extends [calcs/sparc_rar_fit.py](calcs/sparc_rar_fit.py)).

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
