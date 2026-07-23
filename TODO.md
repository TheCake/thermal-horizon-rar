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
2h. ~~**Re-run the error budgets at g_N,ext=1.2**~~ — DONE (Stages 3U/3V). **FINAL:
   α = 1.18 ± 0.11 (simple) / 1.13 ± 0.13 (BE), interior 1000/1000 & 998/1000, Newton
   excluded in all 2000 bootstrap contests (min +53). THE MEASUREMENT PROGRAM IS
   COMPLETE — remaining work is writing (#10) + disclosed open items.** (Chae-text
   verification still owed before the comparison paragraph.)
2e. ~~**The γ≈82° perpendicular excess**~~ — RESOLVED IN IDENTITY (Stage 4J,
   [calcs/stage4j_gamma82.py](calcs/stage4j_gamma82.py)): the islands are the
   pericenter/apocenter faces of the near-parabolic population + an unmodeled
   closest-approach arm (v7's flyby template is γ∈[0,30°]-only). **NEW RESULT: the
   perpendicular velocity CEILING — 11 pairs in the Newton-forbidden band
   [√2, 1.67) (leakage null P=3.8e-9), cliff exactly at the α=1 boosted escape edge
   √(2·1.36)≈1.65 (P=0.62/0.91 consistent); flybys/triples/chance/selection each
   excluded by their own numbers. Small-N caveat (11 pairs) stated.**
2e-b. ~~**GPU model variant**~~ — DONE (Stage 4M, [calcs/stage4m_fly90.py](calcs/stage4m_fly90.py)).
   **Pure 90°-arm flyby template: +4.7/+3.4 lnL, α̂ moved ≤0.06 (1.03–1.16 interior
   8/8, w_rad=0.2 8/8); 50/50 mix and e-ceiling lift rejected/null. The
   γ≈82°/ceiling residual formally closed as a model-shape gap; boost untouched.**
3. ~~**Error budget on lnL**~~ — DONE (Stages 3H+3I). Realization (6 seeds): simple
   0.93±0.11, BE 1.30±0.15, interior 12/12, Newton +55±4. Bootstrap (1000 replicates):
   simple 0.98±0.20, BE 1.21±0.26, Newton +60±11 (min +30). **Combined: simple
   α = 0.98 ± 0.23; BE α = 1.21 ± 0.30** — conditional on the σ_m nuisance being benign
   (see #2b/v5).
4. ~~**Mass-model systematic**~~ — DONE (Stage 4H, [calcs/stage4h_p_ml.py](calcs/stage4h_p_ml.py)).
   **p = 0.578 +0.121/−0.115 with M/L marginalized** (f_d = 1.22±0.10; error doubles;
   p=½ now 0.7σ inside; Cassini floor comfortably passed). The ṽ mass-error leg was
   closed by 3J (σ_m(mass) = 0.024, 12× below the broadening). Quote the marginalized
   p as primary.
5. ~~**Chance-alignment stress test**~~ — DONE (Stage 4I, [calcs/stage4i_rchance.py](calcs/stage4i_rchance.py)).
   R_chance scanned 0.1 → 0.0005: boost ratio stable (max drift 0.033, all CIs
   overlap); wide-bin median drift runs OPPOSITE to a contamination bias. Closed.

## Reconciliation (the paper's credibility keystone)
6. ~~**Reproduce the Banik-style statistic**~~ — EXECUTED at proxy grade (Stage 4N,
   [calcs/stage4n_banikstyle.py](calcs/stage4n_banikstyle.py)): **two modeling choices
   manufacture ~2/3 of the "16σ Newton" result on our own data — unfencing the
   companion fraction absorbs ~60 lnL (needing fractions 3J photometry forbids;
   Banik's free fit lands at 69%), dropping the deep anchor + ṽ-only adds ~5–14 and
   biases α̂ to 0.7; the detection never flips (+30–38 survives full freedom); the γ
   channel is what protects the α measurement (vtonly: α̂_BE wanders to 1.55).**
   Residual future work: line-by-line reproduction incl. their sub-error unconvolved
   binning (H&C's third defect — not honestly ablatable in our always-convolved
   pipeline); multi-seed budget of the ablations if referees ask.

## Theory
7. **Literature hardening** — citation-graph walk of McGaugh+16 (~1000 citing papers) for
   any prior statement of the Bose–Einstein occupation identity; deepen the Pazy & Argaman
   comparison (their statistics-on-screens vs our occupation reading).
8. **Two-field BE theory** — a proper Lagrangian whose EFE structure follows from the
   occupation postulate (the sphericalized QUMOND solve is a stand-in); predict the
   anisotropic (curl) component; feed back through the orbit engine. **Post-4K job
   description: produce the binary boost (α≈1.15 at e_N=1.2) while keeping the solar
   quadrupole under the Cassini cap Q₂ ≤ 9e-27 s⁻² — a 4.3× suppression that
   AQUAL/QUMOND cannot deliver with RAR-compatible ν.**
9. ~~**NLO test**~~ — DONE (Stages 4A/4B; caveat added by 4E). Truncated-expansion
   estimator is power-limited on SPARC (c₁ ± 0.4–0.6, honest null); the truncation-free
   full-function branch comparison: **c₁=½ branch (BE/simple) beats c₁=0 branch
   (standard-μ) in 198–200/200 galaxy bootstraps under raw χ²; Stage 4E's
   scatter-marginalized likelihood deflates this to a strong sign-robust lean (Δ−2lnL
   −56, 166/200, deep-window alone agnostic) — quote BOTH treatments.** Within-branch
   (NNLO 1/12 vs 1/8): no verdict in either treatment (the earlier "slight simple lean"
   was a raw-χ² artifact — retracted in 4E). a₀=1.206e-10 (raw) / 1.00e-10 (honest
   likelihood), f_ML ≈ 1.1 recovered free. Remaining BE kill test: a₀∝H(z) (#14).

18. **MI-vs-MG on our own data — THE post-4K question.** Stage 4K showed the
    modified-gravity (AQUAL/QUMOND) reading of our boost predicts a solar quadrupole
    4× over the Cassini cap (binary-calibrated DHF-2024 tension; escapes: modified
    inertia or EFE-screened theories). Milgrom's modified-inertia predicts
    TRAJECTORY-DEPENDENT boosts (circular vs radial orbits differ) and no EFE
    quadrupole of the capped type. Our 2D (ṽ×γ) likelihood + fitted e-mixture +
    the 4J perpendicular ceiling are exactly the discriminating instrument: implement
    an e-dependent/phase-dependent boost variant in the orbit engine, refit v7,
    compare lnL vs the MG (field-side) boost. Scout verdicts (2026-07-23): the
    binary-eccentricity MI-vs-MG test = NOT FOUND in the literature (nearest: Paci+
    2020 arXiv:2001.03348, rotation-curve MG-vs-MI, MG favored 6.9σ — cite as the
    general-question prior); Milgrom 2011 (arXiv:1111.1611) gives the eccentric-orbit
    MI principle QUALITATIVELY only (no closed form → we bracket with mi_a/mi_t
    prescriptions, labeled as representative); MI-evades-Cassini NOT claimed in print
    (DHF frame the quadrupole as the diagnostic). **EXECUTED (Stage 4L,
    [calcs/stage4l_mi_runner.py](calcs/stage4l_mi_runner.py), seed 31, 4 MI brackets
    × 2 laws): MG WINS ALL 8 CONTESTS — MI-EFE −9…−16 lnL, MI-no-EFE −26…−31 with
    α̂→0.5 sharply (the data DEMAND the EFE suppression); Newton loses under every
    model class (+71…+109); w_rad=0.20 in all 8; BE−simple ≈ −6…−8 everywhere. The
    4K paradox fully quantified: MG fits best but breaks Saturn; MI brackets survive
    Saturn but fit worse; naive no-EFE MI dead. Surviving theory space: EFE-screened
    MG (#8), time-nonlocal MI (unconstructed), or systematics. First-of-kind
    (provisional). **4L-b BUDGET DONE — CORRECTION #10: "MG wins all 8" was seed-31
    luck. Budget verdict: mi_t (EFE on) TIES MG (−3.5±3.3 simple / −0.8±2.5 BE);
    mi_a mildly MG-leaning (simple only, −9.3±2.3); no-EFE MI dead 12/12 (−20…−28,
    α̂→0.5–0.6 sharp = the data demand the EFE amplitude); Newton dead 24/24
    (+71…+108). The Saturn-safe EFE-respecting MI branch fits the binaries AS WELL
    as MG — the 4K paradox has an open door (amplitude pinned, mechanism agnostic).**
19. ~~**Solar quadrupole from our solver**~~ — DONE (Stage 4K,
    [calcs/stage4k_quadrupole.py](calcs/stage4k_quadrupole.py)): Q₂ = 3.9e-26·(α/1.15)
    s⁻² for BOTH families (transition-sourced, family-blind); exceeds Hees+ 2014
    Cassini (3±3)e-27 by ~4.3× — independent binary-calibrated reproduction of the
    Desmond–Hees–Famaey 2024 tension, immune to their M/L + bulge escapes; solver
    cross-validated against Blanchet–Novak μ₁ at 15% (G6). MOND-Planet-9 mechanism
    (2304.00576) runs on the same capped term — disfavored. Full credence
    bookkeeping in NOTES 4K.

## Oscillator program (the current front, opened 2026-07-23)
O1. ~~**Full-hierarchy second moment**~~ — DONE (Stage 4W): **identifiability
    boundary demonstrated** — vertical (distance/inclination) freedom absorbs the
    thermal term (N̂→189, osc+floor +24 vs const) but the calibrated injection
    gate FAILS (slope 0.23 vs 0.93 expected) = the channels are degenerate on
    SPARC; decision needs external distance anchors (Cepheid/TRGB). **The x≈1
    transition bump SURVIVES everything** (0.054 vs 0.026–0.035 dex; not bulge,
    not vertical) → the environmental-g_ext-scatter hypothesis feeds O4. Floor
    converged to 0.035 dex ≈ Desmond-2023's 0.034 (cross-check ✓).
O1b. **External-anchor second moment** (the unlock): import Cepheid/TRGB distance
    subsample (SPARC f_D flags) → refit 4W with those galaxies' σ_v pinned small;
    the thermal-vs-vertical degeneracy breaks where distances are known.
O2. ~~**Binary-side c₁**~~ — DONE (Stage 4X): **binary ĉ₁ = 0.37–0.50** (peaks
    λ=0.75/1.00 across seeds 31/101; σ_λ≈0.08 curvature); **c₁=0 rejected ΔlnL≈20
    per seed, α grid-edge-riding at low λ = shape rejection**; two disconnected
    systems now read the same dial (galaxies 0.43±0.27). Optional deepening:
    the full 6-seed budget on the λ grid; a finer grid near 0.75–1.0.
O3. ~~**Hierarchical-M/L ĉ₁**~~ — DONE (Stage 4Z), and it RESHAPED the claim
    instead of tightening it: hierarchical profile λ̂=0.516 → **ĉ₁ = 0.258
    (0.208–0.309)** vs flat-M/L 0.450; gates PASS (injection λ=1 → 1.03, so the
    ±0.2 profile spread is real structure, not bias); bootstraps agree across
    treatments (≈0.4±0.3). **Galaxy dial: c₁ ∈ 0.26–0.45, 0 excluded everywhere,
    ¼-vs-½ OPEN** (hier profile peaks AT ¼). Paper updated throughout.
O6. ~~**Bath-matrix contest under hierarchical M/L**~~ — DONE (Stages 5C/5D,
    convergence-hardened): **the matrix FLIPS — boot (¼, quantum
    self-consistent, ν=1+n_BE(νy)) beats BE by 75.6 and simple by 174;
    truth-calibrated injections all recover (observed gap ≈ the boot-truth
    calibration −98, not the BE-truth +38); bootstrap 43/50 (86%) boot.**
    Simple falls to LAST hierarchically. Scout: the implicit boot function is
    NOT FOUND in the literature (nor the simple-ν implicit identity) —
    apparently unpublished, "first" withheld pending deeper pass. → opened O7.
O4. ~~**EFE-curve test**~~ — DONE (Stages 5B/5E, convergence-hardened): Chae+21
    Table 3 extracted (94/153 matched); exact Eq.-(2) collinear template; all
    gates pass (nesting, β-injections 3→2.97 / 0→0.00). **Thermal SURVIVES the
    environmental control: β̂=0.044, upper D1 ≈0.02; β=1 EXCLUDED at +168; the
    −7 credit is scramble-generic; correlation channel null (depth-partialed
    ρ=−0.12, p=0.34); the x≈1 bump untouched.** Deep-end EFE visibility limit
    ≲2–5% of collinear max-clustering. Caveats logged (collinear=maximal;
    fixed SPARC distances; S-sample thermal credit −4.6 vs −22.8 full).
O7. ~~**Boot on the binaries**~~ — DONE (Stage 5F): **the binaries VETO boot**
    (+17–24 lnL behind the ½-branch, α̂ grid-edge-riding both seeds = shape
    rejection of the weak transition ν(1)=1.35; Newton still dead +75/+81).
O8. ~~**Tail-vs-coefficient decomposition**~~ — DONE (Stage 5G): the 5D
    hier ladder is monotone in tail sharpness; ν_p (≡ occupation at p=½)
    through the converged hier machinery gives **p̂_hier ≈ 0.65, −56.4 over
    p=½ = 75% of boot's flip from the tail dial alone** at ½-branch-grade
    transition; convergent with §3's p=0.578±0.12. ¼-vs-½ deep digit stays
    OPEN; the hierarchical galaxies vote for SHARPER SCREENING.
O9. ~~**Unification check**~~ — DONE (Stage 5H): **the binaries ACCEPT
    ν_p(0.65)** — α̂ 1.53/1.38 INTERIOR (no edge-ride), Newton +93/+91,
    concession to BE only 4.6/7.6 (realization grade) vs the −56.4
    hier-galaxy gain: **the only function tested viable on BOTH systems**.
    The §2.4 construction spec is now: ½-grade transition + p≈0.6–0.65
    tail + thermal deep structure.
O10. **p-family follow-ups** (opened by 5G/5H): (a) joint-treatment p
    measurement adding the binary channel (finer p grid on binaries, 6
    seeds); (b) self-consistent p-family a₀ scorecard (binary α̂≈1.4–1.5
    under p065 moves the naive translation worse, but the galaxy-side a₀
    moves too — recompute the ladder jointly before any tension claim);
    (c) ν_p(0.65) solar quadrupole through the 4K machinery (sharper tail
    should relieve; verify); (d) 5F/5H six-seed budgets.
O5. Theory question (parked with #8/#18): a construction meeting the §2.4
    six-point specification; the Unruh-ratio reading (x = T_U/T_dS) is the hint.

## Publication path (DEFERRED by user 2026-07-23 — oscillator program first; #11/#12 parked until called)
10. ~~**Paper assembly + verification pass**~~ — DONE. Draft v1 ([PAPER.md](PAPER.md)):
    full manuscript from NOTES Stages 1→4N. **Verification pass DONE (Stage 4O,
    2026-07-22) → v1.1: author = Filip Hájek (independent researcher); nine-agent
    primary-source sweep of every novelty claim; all flagged full-text checks closed
    (Chae sentence verified + quoted, Timeflow = Trofimov CQG 43, 135020 — distinct
    mechanism, 2304.00576 = Jones-Smith & Mathur); references INSPIRE-verified;
    CORRECTION #11 (√2-ceiling crediting) and the phantom "Paci 2020" → Petersen &
    Lelli 2020 fix applied; CITATION.cff + .zenodo.json added for the DOI record.**
10b. ~~**Cookson et al. 2026 (arXiv:2602.24035) full methodological read**~~ — DONE
    (Stage 4P; the "Desmond" byline was a scout hallucination — lead author is
    Cookson, an independent researcher, with Banik/El-Badry/Sutherland/Penoyre/
    Pittordis/Clarke). Median-ṽ flatness test, 1,421 RV-clean pairs ≤130 pc, ~1500×
    Newton preference (2.7σ). Engaged in PAPER §7.4: their e-constancy assumption
    contradicts the Hwang+22 trend they cite (measured sign SUPPRESSES a step);
    their √2-band-impractical remark is answered by our perpendicular column; their
    g_N,e=1.184a₀ confirms our 3T conversion to 3%. **The read caught a real gap in
    OUR pipeline → 4Q/correction #12 (below).**
2i. ~~**Corrected-velocity v7 re-fit**~~ — DONE (Stage 4R,
    [calcs/stage4r_corrected_refit.py](calcs/stage4r_corrected_refit.py), six paired
    seeds, both laws): **α̂ shift −0.012±0.012 (simple) / −0.022±0.015 (BE) —
    consistent with zero; Newton cedes ~5 lnL, stands at +103±9/+90±7; w_rad=0.2
    12/12; simple-over-BE lean unchanged (−12.6±2.4). §6.3 numbers stand,
    correction executed.** Residual sub-item: reconcile the 4J-vs-4Q ceiling census
    boundary convention (11 vs 10 raw; ±1 at edges) — cosmetic, pre-arXiv.
11. **External review** — hand the repo to a professional (ask them to break the
    ceiling null, the quadrupole conversions, and the ablation fairness — the three
    targeted asks now in COLLEAGUE-BRIEF.md, final numbers, held uncommitted);
    they are also the arXiv endorsement path.
12. **Zenodo DOI snapshot** — metadata files in place; next action is the user's:
    enable the repo in Zenodo's GitHub integration, tag release v1.0, DOI mints
    automatically; then send the brief with the DOI. arXiv preprint after review.

## Longer horizon (data that decides things)
13. Gaia DR4: per-mille orientation statistic; α at distribution level with 10× pairs.
14. High-z rotation curves (JWST/ALMA): a₀ ∝ H(z) vs constant — the BE reading's kill test.
15. Wide-binary eccentricity distribution at 10–50 kAU (the η=edge hint): measurable with
    Hwang's v-angle method on the EDR3 catalog — possible standalone paper.
17. ~~**The bath matrix / ¼-branch test**~~ — DONE (Stage 4F,
    [calcs/stage4f_bathmatrix.py](calcs/stage4f_bathmatrix.py)). The thermal reading's
    2×2 matrix (statistics × frequency prescription): **simple-ν IS the classical
    self-consistent bath (exact two-line identity — priority scout pending)**; the new
    quantum-bootstrap cell (c₁=¼, c₂=7/96) is **dead-grade under raw χ² (4–7/200),
    disfavored under the honest likelihood (+27/+9, 62–89/200), sign-consistent**;
    c₁ dose-response peaks at ½ (0 dead, ¼ disfavored, ½ preferred, 1 screened out).
    BE-vs-simple = quantum-vs-classical bath = rung 2, still parked. Scout: no exact
    prior for either result; residue = write-time full-text checks on Timeflow-2024
    (IOP, claims simple-μ as thermodynamic equilibrium) + Zhao astro-ph/0512425, and
    the INSPIRE pass.
16. ~~**The Bernoulli-ladder test on the lensing RAR**~~ — DONE (Stage 4E,
    [calcs/stage4e_lensing_rar.py](calcs/stage4e_lensing_rar.py)). Both public datasets
    fetched (Mistele+24 Table 1 exact-deprojection = primary; Brouwer+21 ESD+covariance
    = cross-check). **Verdict: honest null — rung 2 is NOT reachable: resolving power
    0.09–0.10σ (the 0.2-dex lensing stellar-mass systematic is the wall; needs ~0.02
    dex or Gaia DR4). Full-function within-branch: BE −19 ± 15, carried by 3/153
    galaxies = no discrimination. Bonus corrections (#8a/b): the SPARC "simple lean"
    was a raw-χ² artifact (retracted — binaries lean simple, SPARC agnostic), and the
    4B branch kill deflates to a strong sign-robust lean (166/200 scatter-marginalized
    vs 198–200/200 raw — paper must quote both).** Newton +2777/+1659 on 15 lensing
    points alone; joint a₀ = (1.00±0.09)e−10 under the honest likelihood.
