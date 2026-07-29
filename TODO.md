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
   AQUAL/QUMOND cannot deliver with RAR-compatible ν.** **SHARPENED by Stage 6W
   (2026-07-24): the natural scalar-sum EFE ν(y+e) is EXCLUDED by the binaries
   (−10.35±2.03, 0/6, interior 6/6 — the composition is now MEASURED; tension
   composition-locked). The surviving MG target: reproduce the ANGLE-AVERAGED
   vector-EFE radial profile at e≈1.2a₀ while suppressing its P₂ component —
   ad hoc unless derived (can the thermal reading derive the isotropized
   kernel? O5-side). MI (#18) meets it naturally and remains the standing
   escape.**
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
O10. ~~**p-family follow-ups**~~ — DONE (Stages 5I/5K/5L, 2026-07-23):
    (a) binary p-profile at FOUR seeds: **the binaries prefer p=½
    sign-consistently 12/12** (p0578 +7.7±0.7, p065 +5.8±1.3, gm +10.4±2.8);
    sharp functions stay viable (α̂ interior; p065 1.40±0.04, gm 1.44±0.08).
    (b) a₀ ladder (κ gates pass; κ=1.22/1.34): **sharp-function binary a₀ =
    1.58/1.57±0.11 = +4.9σ off cH₀/2π (vs +1.9σ at ½); galaxy FLAT legs ON
    the horizon (gm −0.2σ)**. (c) quadrupole: **AMPLITUDE-LOCKED** — raw q
    drops (0.81×/0.85×) but α̂ rises in proportion; Q₂·α̂ ≈ 4.2–4.4× Cassini
    for every function; no interpolating-function escape; MI door stands.
O5. **Construction** — MAJOR PROGRESS (Stage 5J): the geometric-mean bath
    **ν_gm = 1 + n_BE(y^¾√ν)** (ω = √(ω_source·ω_total); zero parameters;
    derived c₁=⅓, c₂=1/12 both gate-verified; ν(1)=1.433; tail p≡¾) —
    **LEADS every hierarchical galaxy treatment** (−85 vs BE plain-hier,
    −41 with the vertical channel), ties flat, binary-viable. Priority:
    scout NOT FOUND (exact form, geometric-mean argument, coefficient ⅓);
    Pazy–Argaman closed by direct read (FD statistics, O(e^{−1/x}) deep
    corrections — cannot anticipate). REMAINING: the microphysics of WHY
    the frequency is the geometric mean (the sharpened construction
    problem); ν_gm NNLO/Bernoulli structure vs the C&T ladder.
O11. ~~**Vertical-channel disambiguation**~~ — DONE (Stage 5M): G1/nesting
    pass; **the a₀/f_ML anomaly RESOLVES (all hier fits return to a₀ =
    1.04–1.13e-10 = horizon); the tail vote HALVES but survives (gm −41,
    p065 −32); boot COLLAPSES (−9, coherent with its binary veto)**.
    Verdict: partial absorption; residual sharp-tail preference carried
    with the 4W identifiability caveat.
O11b. ~~**dv-ON galaxy bootstrap**~~ — DONE (Stage 5N): **gm − BE = −29 ±
    53, gm ahead 29/40 (72%)** — the vertical-robust lead is a LEAN, not a
    detection; quoted as such everywhere. Regression vs 5M exact.
O12. ~~**Binary 6-seed completion**~~ — DONE (Stage 5O): **occupation law
    preferred 18/18 (p065 +5.2±0.9, gm +8.5±2.1)**; α̂ six-seed 1.078±0.023
    / 1.41±0.04 / 1.36±0.07, all interior; sharp-function a₀ translations
    1.59 (+5.1σ) / 1.51 (+4.3σ) vs occupation ≈ +1.8σ.
O13. ~~**The β-family follow-ons**~~ — DONE (Stages 5Q/5R/5S/5T, 2026-07-23):
    (a) binary β-profile (5R): **SHARP — β=0 beats every β>0 on every seed
    (24/24); bounds β<0.030 (1σ), β<0.121 (2σ)**; b075/boot edge-ride α̂→2
    (shape rejection); b025 interior α̂=1.57±0.11. (b) configuration test
    (5T): **the galaxy β-vote DECOMPOSES — ultra-deep votes AGAINST β>0
    (+9.3 at ½: wants c₁=½), tail carries it (−61 at y>1), transition flips
    against at β=1; HIGH-arm free fit β̂=0.76 interior; LOW-arm free fit
    ridge-flagged (f→1.94). β RUNS WITH REGIME: deep+transition→0, tail→½–¾
    — same pattern as the two-system split, now intra-SPARC. The family's
    c₁·p_tail=¼ lock is what the data strain against (they ask ~0.3).**
    (c) c₃(β)=β(3β+1)/(24(1+β)³) (5Q): β=0 the unique Bernoulli zero;
    c₄(0)=−1/720 exact; rescaled ladder polynomial (c₂p²=1/48+5β/96,
    c₃p³=β(3β+1)/192); third log-cumulant vanishes at β=½. (d) quadrupole
    scan (5S): **every member 4.0–5.8× Cassini (edge members = lower
    bounds) — Saturn's veto is β-blind, no family escape.**
O14. ~~**The lock-breaking function**~~ — DONE (Stages 5U/5V/5W/5X,
    2026-07-23): DERIVED, not engineered — the spontaneous-fraction
    weighting β = 1/(2ν) (F1; energy-share sibling F2 = 1/(2(2ν−1))),
    c₁=½ and p=¾ EXACT and pre-registered (5U, all gates). **5V: F2/F1
    lead the ENTIRE galaxy ladder both treatments (−107.2/−101.3 plain,
    −51.9/−44.3 vertical vs BE; gm was −84.8/−42.7) at horizon-adjacent
    a₀. 5W: binaries close a third of the gap (−5.5±2.3 / −5.8±1.9 vs
    gm's −8.5), interior 12/12, one seed flips positive — pre-stated
    bands straddled: improved, NOT accepted (residual lean ~2.4 SE).
    5X: Q₂ = 4.8×/5.3× Cassini — amplitude lock holds.** The
    spontaneous-fraction bath = the program's best cross-system function
    (indicative joint −41 vs gm's −26). Scout: running exponent /
    Einstein-coefficient derivation / (½,¾) combo all NOT FOUND.
O15. ~~**The running-β completion**~~ — DONE (Stages 5Y/5Z/6A/6B/6C,
    2026-07-23 late), verdict = THE SPLIT IS SYSTEM-LEVEL: (a) κ ladder:
    κ = 1.26–1.35 for ALL running functions → **a₀ translations +5.2σ/
    +6.3σ/+6.2σ (F4 edge-invalidated); only pure BE passes (+1.9σ)** —
    the temperature wall. (b) two-leg refinement (5Z: β=1/(2ν²),
    1/(2(2ν−1)²); c₁=½ AND c₂=1/12 EXACT; ν(1)=1.503/1.537): galaxies
    reward it AGAIN (6A: F3 −111.4/−50.9, **F4 −108.7/−64.2 = biggest
    controlled lead ever**) but binaries REJECT it like every sharpened
    function (6B: −6.8/−7.0, rb4 α̂-edge 6/6) — the penalty is ~constant
    −5…−8.5 across eight functions INDEPENDENT of ν(1): **the binaries
    reject the screening-region behavior under the dominant external
    field, not the transition value.** (c) 6C bootstrap: **F4 − BE =
    −57.4 ± 38.3, 37/40 (92.5%)** — strongest function-lead grade in the
    program, still a lean. (d) bath-matrix re-run: skipped (6A subsumes
    the ladder question for the new functions).
O16. **The system-level split — FIRST PASS RESOLVED** (6D exclusion +
    6E/6F/6G, 2026-07-23 night): the pointwise drive-weighted rule
    EXCLUDED (6D: −9.64±1.49, 0/6; temperature row passes +2.1σ). Then
    the AMBIENT-GATED bath DERIVED (6E, sign lesson logged: weak fields
    = classical occupied ambients; the gate is stimulated reservoir
    assistance): **β = [n_amb/(1+n_amb)]²·½/(2ν−1)², zero params;
    c₁=½, c₂=1/12 g-independent, c₃=−g/16 exact; p = ½+g/4 post-dicts
    BOTH measured tails (0.689 gal / 0.529 bin).** Tested vs
    pre-registered bars: **6G binaries ACCEPT (−0.88±2.66, interior
    6/6, κ=0.924, a₀ = 1.28±0.15 = +1.6σ — best temperature row in
    program); 6F galaxies hold −59.05 vertical (PASS ≤−50 bar;
    second-best ever) / −92.4 plain (MISS ≤−100 bar by 8, disclosed;
    ambient band −82.6…−97.4). Joint −57.3 = FIRST single function to
    pass both systems.** Post-hoc flag carried. UPDATE (6H/6I,
    2026-07-23): **(a) Chae per-galaxy leg DONE — measured ambients
    IMPROVE both treatments (vertical −61.68 vs −59.05; plain −100.51
    maxclust / −105.94 noclust = the 6F ≤−100 bar RESOLVED; zero
    params; noclust-vs-maxclust ≈ 5-lnL ambiguity carried); (d) AMB
    quadrupole DONE — Q₂ = 3.60e−26 = 4.0× Cassini, no rescue,
    recorded.** **(b) 6J bootstrap DONE — AMB(pgmax) − BE = −56.7 ± 35.7,
    37/40 (92.5%), percentiles [−88.7, −59.9, −15.7]: same top grade
    as 6C's F4 (37/40) while binary-compatible; a strong lean, not a
    detection.** REMAINING: (c) DR4 weak-ambient binaries (sharpen
    toward p=0.69) and eccentricity-resolved MI discrimination — the
    out-of-sample deciders.
O5. Theory, status after 6H/6I/6K: the grammar β = ½·[q_loc·s_amb]^L is
    now PART-MEASURED — the JC (2n+1) anchor makes the energy-share
    branch derivation-grade (6H G1), and the leg count is READ but not
    measured: **L = 2 point-preferred (6I: −52.8/−59.1/−54.5), 6L
    bootstrap DEFLATED the measurement claim — CORRECTION #13 (L1
    rejection 29/40 = lean; L2-vs-L3 21/40 = unresolved; strictly-best
    10/40; the deep arm of 153 galaxies is too thin to read the
    c_{L+1} rung at population grade — deeper samples read the dial
    directly)**. New
    exact rung on file: c₄(L=2) = s²/192 − 1/720 (sign flip at galaxy
    gates; future-data falsifier). Lab-native form (6K G0):
    **β = ½·tanh²(x_loc/2)·e^(−2x_amb)** — the squared gate = the
    Boltzmann cost of two borrowed ambient quanta. **6K v1 (desktop
    analog, pre-registered at ddc83bc BEFORE execution): rate-based
    realizations EXCLUDED — vanilla two-channel competition gives
    CONSTANT β (sky-dead per 5T) and the mediated jump class gives
    wrong-sign running (λ rises with occupation) + a flagged estimator
    pathology; κ-dependence 0.3–0.5/decade everywhere = the grammar's
    rate-freedom does not emerge from jump competition. Pre-committed
    strike LOGGED: bath-microphysics conditional ~20–25% → ~15%,
    surviving mass on the coherent-pull reading.** 6K-v2 EXECUTED as
    **6M (structured-bath, pre-reg 2b3bf78): AMBIG per bands** — the
    single-bath finite-resolution interpolation is κ-free by
    construction (fixes v1's objection) but the estimator drowned
    (endpoint gap K-order vs kernel systematic b²-order; overshoot
    columns); share test failed 205× where clean. Beyond-band
    observations logged: (a) clean-region monotonicity OPPOSITE to
    the grammar (occupied → fully dressed, sparse → source-pinned);
    (b) resolution-β's gravity shape = transition-peaked, tail-zero =
    wrong sky shape both non-deep regimes. Credence HELD ~15% per
    discipline (AMBIG carries no committed move), lean noted. 6K-v3 EXECUTED as
    **6N (pseudomode, pre-reg 6724597): CLOSE-OPP — THE LAB LEG IS
    CLOSED.** The cancellation theorem (flat bath ⇒ exact
    source-locking; confirmed to |λ| ≤ 0.0065 in the full two-mode
    Lindblad = the buildable fridge configuration measures ZERO);
    the KMS-contrast carrier works but is bath-geographic, not
    share-organized (34× off tanh²); the contrast dial's sky shape is
    transition-peaked/tail-vanishing = the sky's tail-only β is
    ANTI-STANDARD. Three realization classes (6K jump / 6M structured
    / 6N exact) all fail to produce the grammar. **Pre-committed:
    bath-microphysics conditional 15% → ~8%; NO v4; fridge moot.**
    O5 REMAINING (the ~8% + alternatives): (a) horizon-specific
    non-Markovian derivation (dS bath = causal structure, not a
    filter bank; bath-width ~T loophole) — the only surviving bath
    home; (b) the MI/trajectory reading (never bath-dependent; DR4
    eccentricity test stands); (c) accept the function as sky
    phenomenology (the ~92% complement — costs nothing observational).
    **(a) FIRST-PASS EXECUTED (Stages 6R/6S/6T, 2026-07-24, pre-reg
    85dcc72): the RESOLUTION BATH — the 6N corollary + the one-scale
    dS bath (KMS width = temperature) force R = νy−√y; β = ½R²/(1+R²);
    exact: ladder preserved through c₄, break c₅ = −1/16, tail = gm
    (p=¾), crossover y* = 1.88 = the 5T turn-on arm; Deser–Levin
    matrix closed analytically (source = no fixed point; total = boot
    + invisible floor; geodesic T_U = 0 = the free-fall derivation of
    binary β<0.03). SKY: galaxies STRONG (−58.6 vertical / −113.7
    plain = plain record) but binaries SHAPE-REJECT (α̂ edge 5/6, 0/6,
    −7.4±2.1 = the eight-function band) ⇒ THE TRIANGULATION:
    local-only dead (6T) + ambient-pointwise dead (6D) + two-factor
    accepted (6G) — the ambient gate is REQUIRED, not chosen.
    **THE GATE DERIVED same day (Stage 6U,
    [calcs/stage6u_gatederiv.py](calcs/stage6u_gatederiv.py)):
    path-resolving the 6H loop splits the (2n+1) pull into absorption
    (n) and emission (n+1) paths; per-leg ratio n/(1+n) = e^(−x) EXACT
    (KMS detailed balance); s^L = Boltzmann cost of L borrowed ambient
    quanta, L=2 = the loop's order — one loop, two vertices, both
    grammar factors per vertex. Uniqueness on measured facts (no new
    fits): rate-balance, absorption-share, inverse-ratio, n², pointwise
    (6D), local-only (6T) ALL excluded; s² unique survivor. Scout: five
    structures NOT FOUND (scout-level). New falsifiable rungs: DR4
    source-vs-dressed x_amb split (Δp≈0.03 at weak ambients); the tail
    exponent RUNS with H(z) (p_gal 0.688→0.702 at z=1) → folded into
    #14. O5 remaining (reading-grade seams): ~~the borrowing
    microphysics~~ — **FIRST CONSTRUCTIVE RESULT (Stage 6X, pre-reg
    67626f7): exact closed-system dynamics in the frozen-horizon limit
    produces the per-leg KMS ratio as the LENDING PROBABILITY
    (ratio slope +0.989, 6.2×/29× margins; P(n≥L) = s^L = e^(−Lx)
    exact; the sky's ratio² = the real-exchange channel; detuned
    virtual control carries raw-n = the discriminating physics).
    ~~Remaining: the dS-side reservoir identification + multimode~~ —
    **DONE (Stage 6Y, pre-reg 0e0f4fc): the EXCLUSION THEOREM — the
    measured gate selects ONE collective ambient mode (M≥2 democratic
    lending saturates the galaxy gate 0.95 and pushes the binary
    postdiction to the rejected edge; negative-binomial tail exact);
    identification (reading): the barycentric mode = the environment's
    dressing cloud (n_amb = ν(e_N)−1); e→0 rejoins the horizon soft
    sector. NEW FALSIFIERS: tail ceiling p ≤ ¾ EXACT (one void galaxy
    beyond ¾ kills the gate); p(e_N) ordering — **6Z EXECUTED: ANTI at
    weak grade (true loses to all 8 shuffles, p≈0.11) ⇒ CORRECTION
    #14: the 6I measured-ambient gains are gate HETEROGENEITY, not
    ordering; the prediction is now fully out-of-sample (DR4-era),
    where failure would be a REAL strike — pre-stated;**
    SATURN COROLLARY: the coupling is trajectory-state = the
    EFE-respecting MI class — the 4K door reached BY DERIVATION (kill
    test: DR4 eccentricity resolution). Remaining reading-grade: the
    barycentric identification itself; ω = g/c (standing §2.4
    postulate); the ½ ceiling candidate (resonant equal-time-sharing).**
    ~~untied-L contest~~ — EXECUTED same day (Stage 6V, pre-reg
    6ebf545): **tied form SURVIVES (B1 TOLERATED — worst untied
    margin +3.33, inside 6L population noise; no strike, no proof);
    B2 instrument-split FAILED as pre-registered and logged (the
    likelihood reads L1 via the TRANSITION channel too; dividend:
    one-leg disfavored −7.1 through a rung-independent channel).**
    Optional refinement: AMB-vs-[resolution×ambient] local-factor
    contest (w(R) vs 1/(2ν−1)² at fixed gate).**
    FROZEN-BATH SUB-READING (the Hubble-correlation-time argument):
    both accessible tests EXECUTED 2026-07-24 and closed — 6O galaxy
    coherence AMBIG (unmeasurable vs per-galaxy systematics even with
    shape leverage; shuffled-null −82); 6P/6P-b binary s-shaped
    scatter: edge-preferred +11/+13 (2 seeds) but the s-FLAT control
    WINS by 26 lnL = generic broadening, NOT the frozen fingerprint.
    No support, no formal strike; sub-branch weakened, ~8% holds.
    NEW low-priority item: v7 error-model refinement — an s-flat
    per-system scatter worth ~+37 lnL exists (3E/3J echo); α̂
    exposure ≤ +0.08 (inside the ±0.11 systematic; measurement
    stands). ~~The curated MEASUREMENT LEDGER + world-table audit~~ —
    DONE (Stage 6Q, 2026-07-24): [LEDGER.csv](LEDGER.csv) (75 rows,
    supersession discipline, spec in [LEDGER.md](LEDGER.md)) +
    [calcs/stage6q_worldtable.py](calcs/stage6q_worldtable.py) (six
    gates all PASS; 13 laws × 7 tests). **Standing rule: every future
    stage that produces a headline number adds its ledger row in the
    same commit; superseding = flipping status + pointer, never
    deleting.**
    The corollary worth keeping: WHICH occupation a self-shifting
    mode takes = whether the bath's KMS structure is resolved across
    the self-shift (the C&T-vs-boot dichotomy in one line). Scout-level: Bernoulli-break↔leg-count
    and (2n+1) selection both NOT FOUND (post-6H sweep).

O18. **The vacuum rung + the platform translation** (Stages 7D/7E,
    2026-07-24, pre-reg c53f2e2/ff1de37). 7D: the classical vacuum-free
    local share (q_cl = 1/(2(ν−1)), no "+1") REJECTED on the galaxy
    ladder at kill-grade ×22 (+556.5 vertical, cap-robust lower bound,
    worse than no admixture at all); binary direct leg instrument-N/A
    (the cliff defeats the solver's identity gate at 10.6% — refinement
    item: a cliff-capable 1D spherical integrator, NOT promised); RJA
    classical-ambient rung UNRESOLVED (weak as pre-stated). Quantum-
    ladder status: occupation rung banked, vacuum-share rung galaxy-
    direct + binary-indirect, ambient-statistics open, fluctuation
    (P6) out-of-sample. 7E: the lending gate as a buildable cQED
    experiment (χ/2π = 250 MHz transmon + thermal 3D cavity + 2 MHz
    beam-splitter): lending law to 0.25%, dissipative shift < 2%,
    saturation ×7.2, the L=2 geometric rung slope 0.987 — mechanism-
    class validation, not a gravity test (6N frame); L=2 prefactor
    0.480 ≈ ½ = 2nd data point for the equal-time-sharing ceiling
    reading.

O17. **The Einstein fluctuation test + the signed prediction ledger**
    (Stage 7A, 2026-07-24, pre-reg 21e7d28). [PREDICTIONS.md](PREDICTIONS.md)
    created: registered-before-test rule, status-flips only; P1–P8 with
    uniqueness classes (DM-immune / classical-immune / generic) + numeric
    kill conditions. P6 (the particle term: scatter γ = 1 quantum vs 2
    classical, a₀ locked to the mean — Einstein 1909 pointed at the RAR)
    EXECUTED same day: **shot/corpuscular bath EXCLUDED (−25.05)**;
    quantum-vs-classical **NO VOTE** (EQ−EC = −0.43; γ profile bimodal,
    axis ≤ 4.4 units); instrument CALIBRATED (G2: paired injections
    recovered 1.03/2.05) ⇒ data-side misspecification, not power — **the
    x≈1 bump occupies the discriminant window (x ≈ 0.4–1.2)**. P6 LIVE,
    both kill directions pre-committed. Flagship
    out-of-sample unique prediction = P3, the z-LOCKED PAIR (a₀(z) =
    cH(z)/2π AND p_gal(z) rise together — one function moves two
    observables; numbers tabled in PREDICTIONS.md).
    **THE HUNT EXECUTED same day (7B/7C, pre-reg d82cc4b/59c582a):
    the bump is CAUGHT — an INNER-DISK (R < 1.5 R_d) phenomenon
    (b_clean = 0.0000 on outer points; radius-mix EXPLAINED = 0.67
    PARTIAL vs 0.75 bar; driver class = streaming/beam-smearing/
    decomposition, scout-anchored). Standing methodological numbers:
    within-curve lag-1 ρ ≈ 0.876/0.844 (every point-level −2lnL in
    the program nominally calibrated only; offsub point floor
    0.04–0.07 dex) and the inner-disk term c_in ≈ 0.105 worth ~83
    lnL/1p. THE CLOSE: γ NOT MEASURABLE on SPARC (informed-design G2
    fails; variants sign-disagree ±5; thinning +0.14) — P6 fully
    out-of-sample (anchored/DR4/IFU with non-circular modeling);
    neither kill direction fired. Scout: the 4T/4W acceleration-binned
    scatter measurement itself apparently unpublished (Stone &
    Courteau/Desmond quote global only — scout-level).**

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

17. **The 7L settling instrument (round 11, named not opened)** — run the
    completeness + width machinery ON the subsample Cookson-grade cuts
    REMOVE (the ΔRV>10 / wide / distant 92%): either its kinematics
    decompose into the measured companion + noise sectors, or they carry
    the boost. Settles "data or contamination" — the now-well-posed form
    of the field's fork (PAPER §7.4(d), ledger bin-7l-cookson).

18. **The width-object identification (the 7J-z6 successor, the program's
    top binary instrument)** — the inner-bin eccentricity/radial-population
    sector (fingerprint: mid-shoulder ṽ, radial γ=8° column, inner bins;
    floor/tail/fpm/sq/perspective excluded as shapes). Identifying it
    collapses the noise-ceiling curve to a point (round-11 form) and
    restores function-discrimination power lost at 7J-d.

## The rivals arc (opened 2026-07-29 — user directive: "shoot down more
## theories and strengthen or shoot down ours"; PAPER.md FROZEN at v3.9
## for the arc's duration, absorption at arc close)
19. ~~**8A the rivals' ladder**~~ — DONE ([calcs/stage8a_ladder.py](calcs/stage8a_ladder.py),
    all gates first-run; NOTES 8A; ledger galfn-rivals-ladder /
    mech-expmu-boot / galfn-hees-pincer): standard/n≥2 family c₁=0
    excluded at the 56.3 grade; additive-analytic class theorem (c₁=1
    exactly: Verlinde exact member, superfluid flagged); **exp-μ ≡ boot
    identity** (¼, 7/96, ν(1)=1.350 — apparently the function's first
    data contest, dead per 5M); Hees+16 catalog ladders (ν̃ c₁=a, ν̄
    1−1/(2a), ν̂₁≡BE gated); **THE PINCER: empty intersection between
    the measured window and their Cassini verdicts across the whole
    catalog.** Residue: primary read of the Hees table before paper
    use; superfluid composition equation primary check.
20. ~~**8B the simulation ladder — the shoot-down-OURS instrument**~~ —
    EXECUTED (2026-07-29, pre-reg c03604f; EAGLE SQL access granted by
    the team; [calcs/stage8b_simladder.py](calcs/stage8b_simladder.py);
    NOTES 8B; ledger galfn-eagle-ladder): **G3 GATE-FAIL / OFF-FAMILY —
    the ΛCDM attractor does NOT print the digit.** The EAGLE aperture
    RAR (7,239 centrals + the Recal box) rides the family boundary at
    every c₁; where localizable it prints c₁ ≈ 1.1; the D1 fingerprint
    contest has the ADDITIVE class beating BE by −6,005 = the 8A
    theorem observed in silico; deep slope runs 0.51–0.62; a₀_sim
    unlocked (1.2–2.5e-10). The sky's ½ stays unique against the
    attractor at aperture grade. Successor (commissioned-only): the
    FIRE-2 snapshot pipeline for beyond-aperture grade; the 4T
    scatter-shape leg on sims likewise waits for disk-plane profiles.
21. **8C ours-sharpening**: (a) the P1 tail-ceiling test (p ≤ ¾,
    parameter-free) on the deepest SPARC galaxies; (b) the c₄-rung
    Fisher forecast (which survey reads the next digit).
22. **8E lead (logged, primary reads owed)**: observational a₀(z) fits
    exist (a₁ ≈ 1.6e−10/z, scout-reported MUSE-DARK era) AND Magneticum
    ΛCDM produces a₀ rising ×3 by z=2.3 ≈ the H(z) direction — P3's
    a₀ leg alone is NOT ΛCDM-immune; the LOCK (a₀ + p paired drift) is
    the discriminator (P3 annotated 2026-07-29).
