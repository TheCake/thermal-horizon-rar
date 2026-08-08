# A horizon-thermal mechanism for the radial acceleration relation: exact results and open constants

**Filip Hájek** (independent researcher) — hfilip11@gmail.com

*Draft 0.1 (2026-08-08), written from the program record (PAPER.md v4.0 and the post-v4.0 ledger through the soft-sector stages); companion papers: Hájek 2026a (wide binaries) and 2026b (coefficients). Figures 1–4 are produced by calcs/paper3_figures.py under provenance gates (data/paper3_figs.txt). References inherit the program's verified lists; the three entries new to this paper (Gibbons & Hawking 1977; Jaynes & Cummings 1963; Pöschl & Teller 1933) are standard classics. Not for circulation.*

*Acknowledgments: the computational analysis, the adversarial review rounds, and the manuscript drafting were performed in collaboration with Claude (Anthropic). Two external adversarial reviews of the soft-sector calculations, run as independent referee sessions, reproduced every load-bearing number and downgraded two verdicts; both downgrades are adopted in the text. The full chronological record, including all logged corrections and the pre-registration commit chain, is public in the repository (Appendix B).*

## Abstract

Two companion papers measure the radial acceleration relation at coefficient level and bound the corresponding velocity boost in wide binaries. This note asks a narrower question: can a thermal-horizon exchange mechanism carry those measurements, and what would kill it? We work in the Cadoni–Tuveri (2019) framework, where the observed interpolating function is a Bose–Einstein occupation, ν = 1 + n_BE(√(g_bar/a₀)), at the de Sitter temperature. We report four results. First, reduction theorems: for linear single-port coupling, a many-mode ambient sector reduces exactly to one collective mode, and every bound orbit in Schwarzschild–de Sitter satisfies Ω/H > 1.449, licensing a frozen-bath treatment. Second, the soft sector is solvable: the l = 0 minimally coupled radial problem in the de Sitter static patch is the Pöschl–Teller well, giving an exact infrared enhancement 1 + H²/ω² of the coupling density. The continuum admixture diverges without an infrared gap, so the data-preferred ambient gap becomes a structural necessity; with it, the leak at physical coupling is four orders of magnitude below the smallest measured exchange weight. Third, the continuum's dispersive pull is a common-mode de Sitter Lamb shift; after a declared two-parameter renormalization the residual ladder distortion is within the measured coefficient band at the fiducial point (deep-arm 0.023 against a band of 0.05), with an order-unity margin controlled by one unpinned normalization constant. Fourth, the exchange law (1/2)·n/(1+n) is verified in closed models and is testable in circuit QED. The conditional credence we assign the mechanism is about 15%.

## 1. Introduction

The companion papers of this series establish two empirical anchors. On galaxies, the radial acceleration relation is measured at coefficient level: the deep-side expansion of the interpolating function carries a nonzero next-to-leading coefficient in a band around one half, the acceleration scale locks onto cH₀/2π under hierarchical treatment, and the Newtonian-side screening is exponential (Hájek 2026b). On wide binaries, the corresponding velocity boost is bounded from above at half the galactic calibration on the cleanest strata, with a small allowed sector remaining and Newton not excluded (Hájek 2026a).

Those papers are deliberately phenomenological. This note is the mechanism companion, and its scope is stated at the outset. We do not claim a derivation of low-acceleration dynamics. We ask whether a specific microphysical reading, gravitational dressing by exchange with the thermal soft sector of the de Sitter horizon, is internally consistent with the measured numbers, which parts of it can be established at theorem grade, which parts fail, and which single quantities now carry the burden. The reading follows Cadoni & Tuveri (2019), who derived the observed interpolating function as a Bose–Einstein occupation law of thermally excited horizon degrees of freedom, with the acceleration scale a₀ = cH₀/2π emerging rather than being fitted. Famaey & Durakovic (2025) note that no clear derivation of the interpolating function otherwise exists in the literature; the Cadoni–Tuveri identity is, in our assessment, the strongest available candidate, and the program of this series treats it as a hypothesis to be tested rather than adopted.

Our method is the same discipline as in the companion papers: verdict thresholds are committed to the public repository before each deciding computation runs, and failed conclusions are retracted in the open. Adversarial review rounds were run against the load-bearing calculations by an independent referee session, whose numerical reproductions and two verdict downgrades are adopted here. Every number in this paper is produced by a named script (Appendix B).

## 2. Measured targets and the dictionary

The identity at the center of the framework is algebraic. With y = g_bar/a₀ and x = √y, the interpolating function preferred by the rotation-curve data satisfies

  ν(y) = 1/(1 − e^(−x)) = 1 + n_BE(x),  n_BE(x) = 1/(e^x − 1).

Reading x as a ratio of frequency to temperature fixes the dictionary: a system whose Newtonian field is g_N carries a characteristic mode frequency ω = √(g_N a₀)/c, the horizon temperature in the same units is T = a₀/c = H/2π, and x = ω/T. Deep systems are soft (ω below T), Newtonian-regime systems are hard, and the occupations implied are large exactly where the anomaly lives. This dictionary is not decoration; every quantitative statement below is made in it.

The mechanism must then carry, at minimum, the following measured facts from the companion papers and the program record. The deep-side coefficient band c₁ ≈ 0.26–0.45 with the one-half branch preferred in direct contests. The tail exponents of the screening, p ≈ 0.65 ± 0.075 on galaxies and p ≈ 1/2 on binaries. The two-system split itself: galaxies and binaries prefer measurably different effective functions, and the difference survives every local-field control, pointing at the ambient environment as the discriminating variable. The exchange-weight reading compresses this split into one curve. Suppose the ambient sector at occupation n couples through its Boltzmann ratio s = n/(1+n) = e^(−x), and the observable weight is the probability of supplying two quanta, s². Then the measured galaxy and binary values sit on one exponential at their measured ambient fields (Figure 1). The same weight postdicts both tail exponents through p = 1/2 + s²/4: 0.688 and 0.528, against the measured 0.65 ± 0.075 and ≈ 1/2. These postdictions were obtained before the tail measurements existed at their current precision and are the strongest single reason to take the exchange reading seriously.

![Figure 1](figs/p3_fig1_gate.png)

*Figure 1. The exchange weight s² = e^(−2x) against the ambient occupation variable. The two points are the measured galaxy and binary values from the program record; the annotations give the tail-exponent postdictions p = 1/2 + s²/4 against the measured values.*

## 3. Reduction theorems

Three structural results hold at theorem grade within the model class and set the terms for everything else.

### 3.1 One collective mode

For linear coupling of the system mode a to K ambient modes b_k through a single port, H_int = Σ_k λ_k(a†b_k + b_k†a), the interaction is exactly H_int = λ̄(a†A + A†a) with A = Σ_k c_k b_k, c_k = λ_k/λ̄, λ̄² = Σ_k λ_k². The orthogonal (dark) combinations decouple from the dynamics, including in the presence of the system's Kerr nonlinearity, and the bright-mode marginal of independent thermal modes is exactly thermal at the weighted occupation. Consequently the two-quantum exchange weight keeps the form s̄^L for the collective mode: the counting statistics of a many-mode thermal ambient are indistinguishable from one effective mode. The theorem is defensive, in the reviewer's apt phrase: it does not produce the mechanism, but it closes the objection that a realistic multimode environment would replace the measured geometric statistics with a many-mode alternative. The data independently reject the leading such alternative (a which-path-distinguishable, number-additive gating produces a negative-binomial tail that saturates at 0.95 where 0.754 is measured).

### 3.2 Bound orbits see a frozen bath

In Schwarzschild–de Sitter spacetime, every bound orbit satisfies Ω/H > 1.449, where Ω is the orbital frequency: the stability region for bound motion terminates before the orbital and cosmological rates become comparable. The measured systems sit far inside: the slowest SPARC orbits have Ω/H ≈ 30, and the wide-binary sample sits near 5 × 10³. On the mechanism side this licenses the quenched treatment used throughout, in which the horizon sector is frozen on orbital timescales and exchanges quanta rather than acting as a Markovian damping bath. The corollary for free fall is exact: a geodesic detector registers the pure horizon temperature with no acceleration supplement, which is the derivation-grade root of the binaries' sharp screening in this framework.

### 3.3 Detuning is unresolvable

The exchange reading requires near-resonant lending between the system mode and the collective ambient mode. Within the model class this cannot be spoiled by mistuning: all couplings are of order H, and no system in the sample has existed for longer than a Hubble time, so the accumulated phase from any detuning in the soft sector is bounded by roughly 0.16 radians. Under every extraction rule tested (coherent, decoherent, finite-time), a soft-sector detuning perturbs the exchange weight at or below the ten percent level. This closes the mistuning objection. It does not establish resonance as a positive mechanism; that distinction is maintained throughout.

## 4. The soft sector is solvable

The objection that survived longest is the leak: a continuum of soft horizon modes, each contributing a small off-resonant admixture that grows linearly with occupation, could swamp the resonant exchange weight. Closing or confirming this required the horizon-side mode density, not another toy model. It turns out the required problem is exactly solvable.

In the de Sitter static patch, with tortoise coordinate x defined by Hr = tanh(Hx), the l = 0 radial equation for a massless scalar has potential V = −2H² sech²(Hx) for minimal coupling and V = 0 for conformal coupling. The minimal case is the λ = 1 Pöschl–Teller well, whose regular solution is elementary, v(x) = tanh(Hx) cos(ωx) + (ω/H) sin(ωx). Flux normalization then gives the exact ratio of near-origin coupling densities,

  D_min(ω)/D_conf(ω) = 1 + H²/ω²,

which localizes the de Sitter infrared enhancement in closed form (Figure 2). Transverse-traceless gravitational perturbations obey the minimally coupled equation, so the enhanced case is the graviton-like one. An independent referee session re-derived the potential from the metric (including the curvature coupling, R = 12H²) and reproduced the ratio exactly.

![Figure 2](figs/p3_fig2_well.png)

*Figure 2. Left: the l = 0 static-patch potentials; the minimally coupled (graviton-like) case is the reflectionless Pöschl–Teller well. Right: the exact infrared enhancement of the coupling density.*

Three consequences follow.

First, a necessity theorem. The thermally weighted continuum admixture for the graviton-like sector diverges linearly in the infrared (the integrand approaches [1/(8π³Ω²)]/ω² exactly), while the conformal sector converges. An infrared gap is therefore required, not optional, for the enhanced sector. The gap the data had already selected, in which each system's ambient environment lifts the deepest modes into a collective mode at the ambient dressing frequency, is thereby promoted from an interpretive reading to a structural requirement. We emphasize the asymmetry of this result: the geometry does not merely tolerate the gap, it demands one, and the ambient gap is the candidate the data independently prefer. Gap sufficiency is a separate question and is not established; at the physical coupling the safety margin is carried by the resonant window rather than the gap itself.

Second, with the gap at the measured ambient frequencies, the leak is bounded. At the physical coupling scale g ≈ H, the continuum admixture is 1.6 × 10⁻⁵ (binaries) and 3.8 × 10⁻⁵ (galaxies), four orders of magnitude below the smallest measured exchange weight (0.112), and it preserves the measured two-system ratio to a part in 10⁴. The small-coupling corner of the scan, where the admixture formally approaches the exchange weight, lies where the linear-response formula itself fails (per-mode admixtures of order ten) and is excluded from any claim in either direction.

Third, locality. Near the origin the l-th partial wave couples as (ωR/c)^l, which evaluates to weights of 10⁻²⁴ (binaries, l = 1) and 10⁻¹⁵ (galaxies, l = 1) at the measured scales; the thermal population of the gapped angular channels contributes at the 4 × 10⁻⁴ level. The single-channel structure of the mechanism is therefore derived rather than assumed, for field-value coupling. Since spin-2 perturbations begin at l = 2, the computed l = 0 channel is a conservative proxy, and careful multipole accounting (the standard long-wavelength law: the 2^l-pole moment against the ω^l field gradient) confirms the physical graviton channel is more suppressed than the proxy, not less. The caveat this leaves is the normalization question of Section 5.

## 5. The dispersive budget and the two open axes

The same continuum exerts a dispersive pull on every system's frequency. Computed exactly at fiducial coupling, the pull has the closed form Δ = −(g²/4π²)·ln(1 + Ω/g)/Ω and is, to within one part in 10³, a pure vacuum effect: a common-mode de Sitter Lamb shift of about −0.02 H, nearly identical for galaxies and binaries. Taken literally it is excluded by the data themselves (it would shift the deep relation catastrophically), so it must renormalize; the honest question is what remains after renormalization.

We declare the scheme rather than assume it. Exactly two universal transformations are absorbable: one additive constant in the bare frequency (the grammar is formulated in dressed frequencies; a constant pull with no system dependence redefines the unobservable zero point) and one multiplicative rescale of x, which is algebraically an a₀ recalibration. Together these are the affine transformations of x, and nothing beyond affine is absorbable. After removing the fitted affine part, the residual is propagated exactly through the occupation function and compared to the measured coefficient band.

The result is two-sided and we state both sides. At the fiducial point the residual is inside every measured band. The deep-arm distortion is 0.005 as a window mean and 0.023 at the most amplified deep end, against a measured coefficient band of width 0.05. The implied binary exchange-weight distortion is 4 × 10⁻⁴ against a 25% measurement precision. The absorbed rescale spends 2.1% of the acceleration-scale calibration, or 0.22σ of the measured a₀ band, a price the temperature lock can afford. The referee session verified, additionally, that the result is independent of the ambient anchor at fiducial coupling (identical for every tested gap frequency).

Against this, the margin is order unity, not comfortable, and it is controlled by two quantities the framework has not pinned (Figure 3). The first is the coupling scale: the budget is comfortable for g ≳ 0.3H and degrades below, though the degraded region is also where linear response fails. The second, and the sharper of the two, is the absolute normalization of the horizon coupling density. The residual scales linearly with it, and a factor-three normalization at g = 0.3H crosses the measured band; the deep-end value crosses the exclusion grade. The same constant governs whether the suppressed physical graviton channel (l = 2) can simultaneously carry the order-unity exchange grammar. These are two faces of one unresolved constant, the absolute strength with which a bound system couples to the horizon's soft spectrum. Pinning it is the single computation that would convert this section's accounting from bounded to closed, or falsify it.

![Figure 3](figs/p3_fig3_budget.png)

*Figure 3. The deep-arm ladder distortion after the declared renormalization, as a function of coupling. Solid: fiducial normalization; dashed: normalization times three (the unpinned axis). The horizontal lines mark the measured coefficient band width and the exclusion grade; the shaded strip marks where linear response becomes marginal.*

## 6. The exchange mechanism at model grade

The positive content of the mechanism, as opposed to its consistency, currently rests at the level of exactly solved closed models plus one derivation-grade decomposition.

In the frozen-horizon limit, a Kerr-shifted system mode exchanging with one thermal ambient mode acquires its dressed configuration with probability weight (1/2)·n/(1+n): one half from equal-time sharing at resonance, and the Boltzmann ratio n/(1+n) = e^(−x) as the probability that the ambient can supply the quantum. The closed model reproduces this law to a slope of 0.989 over the relevant occupation range, the generalization P(supply L quanta) = s^L is exact, and the resonant channel is coupling-independent while a detuned virtual channel is not: the exchange is real, not virtual (Figure 4). A path-resolved decomposition of the standard dispersive pull (Jaynes & Cummings 1963) supplies the derivation-grade reading of the two factors: the per-vertex ratio of absorption to emission orderings is exactly the Kubo–Martin–Schwinger weight e^(−x), and the loop order fixes the exponent at two. The laboratory translation is direct: with a transmon's dispersive shift standing in for the Kerr term and a calibrated thermal input, the same law is measurable at the fraction-of-a-percent level in circuit QED. The two-quantum prefactor observed in the translated model, 0.480, is the closest thing the one-half ceiling currently has to an independent check. A clean laboratory negative would strike the exchange dynamics; we consider this a feature.

![Figure 4](figs/p3_fig4_lending.png)

*Figure 4. The exchange-ward weight against ambient occupation, with the closed-model and circuit-QED-translation checks and the two sky operating points.*

Equally important is what failed. Three families of rate-based realizations were built and excluded in pre-registered rounds. Plain two-channel jump competition produces a constant weight where the sky demands a running one. Mediated-jump constructions run with the wrong sign. A flat-bath configuration obeys an exact source-locking theorem: the buildable analog measures zero, with the full two-mode Lindblad confirming decoupling at the 10⁻² level. A pointwise field-weighted variant is excluded by the binary data directly, and the many-mode counting alternative by the tail shape. These exclusions each carried a pre-committed credence penalty and are the main reason the mechanism's conditional credence is low despite the theorem-grade consistency results above.

## 7. What is not established

We list the open items with the same prominence as the results.

The averaging question. The measured galaxy tail sits between the un-averaged and fully averaged limits of the exchange weight, and the natural discriminator, the tail-exponent precision, hits a realization wall: the honest uncertainty on the galaxy tail exponent is σ_p ≈ 0.075 where 0.02 would decide. A direct fit of the exchange weight gives r̂ = 0.34 ± 0.19: nonzero at two-sigma lean, cross-validated by two estimators that agree to one part in 10³, but not a detection. The one-half ceiling of the weight is consistent with, and no better than, the two model-grade checks above.

The reservoir identification. The collective ambient mode is identified with the environment's dressing cloud by counting statistics and by exclusion, not by construction from horizon microphysics. The state prescription for the enhanced sector (the thermal weighting of a field with no de Sitter-invariant vacuum) is the standard static-patch choice and is named as a choice.

The normalization constant. As Section 5 states, the absolute coupling-density normalization is unpinned, it now carries the closure of the dispersive budget, and it is entangled with the spin-2 channel question. This is, in our judgment, the sharpest open problem the mechanism has.

The non-thermality channel. A dial connecting ambient-state non-thermality to the exchange weight exists but is conditional on an untested functional choice; its DR4-era lever is labeled accordingly in the prediction list and claims nothing today.

## 8. Predictions and falsifiers

The framework's registered predictions, maintained with numeric kill conditions in the repository, include the following.

A ceiling on the galaxy tail exponent, p ≤ 3/4 exactly, approached in low-ambient environments; the void asymptote is a gate-free meter with a kill band of [0.727, 0.750] at the ceiling. The band is state-independent: thermal and strongly squeezed ambient states move the band edge by less than 0.1%, so the prediction survives ignorance of the horizon state.

A locked pair of drifts with redshift: a₀(z) = cH(z)/2π and the galaxy tail exponent rising together (0.689 to 0.702 at z = 1). Either drift alone can be mimicked; the locked pair is the framework's signature. Current intermediate-redshift measurements trend in the predicted direction at low significance.

Gaia DR4 items: the source-versus-dressed distinction in the ambient variable (a tail-exponent difference of about 0.025 at intermediate external fields), the environment-ordering of the exchange weight (fully out-of-sample after an in-sample version was withdrawn), and the conditional non-thermality lever.

A laboratory test: the exchange law and its two-quantum rung in circuit QED at stated parameters. This tests the mechanism class, not gravity; a negative would remove the dynamical basis of the exchange reading.

The kill levers are equally explicit. A galaxy-tail measurement at σ_p ≤ 0.02 finding p < 0.67 falsifies the averaged exchange weight. A normalization computation finding the coupling density a factor of a few above fiducial breaks the dispersive budget of Section 5. A void-environment tail beyond 3/4 kills the gate structure outright.

## 9. Discussion

The mechanism's standing after this investigation can be stated in one sentence: the exchange reading is now consistent at theorem grade where it used to be conjectural, it remains underived where it always was, and it is falsifiable on named axes rather than in principle. The theorem-grade half is the collective mode, the frozen bath, the unresolvable detuning, and the solvable soft sector with its required gap and bounded leak. The underived half is the averaging step and the absolute normalization.

Relation to the literature. The identity and the temperature are Cadoni & Tuveri's (2019); this series contributes the coefficient-level tests (Hájek 2026b) and the mechanism consistency program reported here. The trajectory-formulation escape from the solar-system quadrupole, which the companion papers require on data grounds, is congenial to modified-inertia proposals in the line of Milgrom (1999). The frozen-bath and free-fall results above give that direction a concrete footing in horizon thermodynamics (Unruh 1976; Gibbons & Hawking 1977; Deser & Levin 1997). We make no claim about relativistic completion; the constraint structure there is set by Skordis & Złośnik (2021). Entropic and emergent-gravity approaches (Verlinde 2017) share the horizon-thermodynamic vocabulary but make different quantitative commitments; the coefficient measurements of the companion paper are the discriminating data.

Credences. Following the program's standing rule, we state ours plainly. The probability that the wide-binary anomaly of Paper I is real physics stands at about 53%. Conditional on the low-acceleration anomaly being real, the probability we assign to the thermal-exchange mechanism in roughly its present form is about 15%. The number is low by construction: three excluded realization families each carried a pre-committed penalty, and the two adversarial reviews of the soft-sector calculations each ended in an adopted downgrade rather than a confirmation. It has also been rising slowly, from 8% after the exclusion rounds, on the strength of the reduction theorems and the solvable soft sector. We consider a mechanism program that publishes its own kill record more informative than one that reports only survivals.

## 10. Conclusions

The soft sector of the de Sitter horizon, treated exactly, neither destroys the thermal-exchange reading of the radial acceleration relation nor yet completes it. The geometry requires the infrared gap the data had already chosen, and the leak the gap controls is negligible at physical coupling. The continuum's dispersive pull renormalizes to a residual inside the measured bands at the fiducial point, and the exchange law itself is exact in closed models and testable on a bench. What remains is not a fog but a list: one normalization constant, one averaging step, one laboratory curve, and three sky measurements with dates attached. The companion papers carry the measurements; this note has tried to make the mechanism carry its own weight, and to make its failure modes as legible as its successes.

## Appendix A. Transparency and corrections

The mechanism line of the program logged its failures with the same discipline as the measurement line; the record holds the full list. Four items changed conclusions and belong here. First, three families of rate-based bath realizations were excluded in pre-registered rounds, each carrying a pre-committed credence penalty; the surviving reading (frozen-horizon exchange) is what this paper reports. Second, an early claim that the soft-sector tail was virtual was retracted when review showed the diagnostic had been calibrated in a saturating regime; the corrected treatment is Section 4. Third, the two soft-sector verdicts were both downgraded on adversarial review (the leak stage to ambiguous pending the dispersive question; the budget stage from closed to gray on definition-sensitivity and the unpinned normalization), and the downgrades are adopted verbatim in Sections 4 and 5. Fourth, an in-sample environment-ordering claim was withdrawn after a permutation control showed the gain was generic; the prediction is now fully out-of-sample (Section 8). The pre-registration commit chain for every verdict, including the bars that were set before each computation ran, is public in the repository.

## Appendix B. Reproducibility

Every number quoted in this paper is produced by a named script in the public repository (github.com/TheCake/thermal-horizon-rar), with outputs archived alongside. The program record (PAPER.md, NOTES, LEDGER.csv) carries the full audit trail.

| Quantity | Script | Output |
|---|---|---|
| Exchange weights, tail postdictions | calcs/stage6e_ambgate.py; calcs/stage6i_chaegate.py | program record |
| Collective-mode reduction theorem | calcs/stage9w_multimode.py | data/stage9w_multimode.txt |
| Bound-orbit stability bound | calcs/stage9q_lemmas.py | program record |
| Detuning phase budget | calcs/stage9t_resonance.py | program record |
| Pöschl–Teller solution, D-ratio, gap necessity, leak | calcs/stage9z_leak.py | data/stage9z_leak.txt |
| Leak scan validity, ambient co-read, UV completion | calcs/stage9z_addendum.py | data/stage9z_addendum.txt |
| Renormalization scheme, dispersive residual, envelope | calcs/stage9zb_diff.py; calcs/stage9zb_addendum.py | data/stage9zb_diff.txt; data/stage9zb_addendum.txt |
| Exchange law, two-quantum rung | calcs/stage6x_borrow.py | program record |
| Circuit-QED translation | calcs/stage7e_platform.py | program record |
| Rate-based exclusions | calcs/stage6k_analog.py; calcs/stage6m_analog2.py; calcs/stage6n_pseudomode.py | program record |
| Void ceiling and state-independence | calcs/stage6y_reservoir.py; calcs/stage9y_squeezing.py | data/stage9y_squeezing.txt |
| Redshift pair prediction | calcs/stage6u_gatederiv.py | program record |
| Figures 1–4 | calcs/paper3_figures.py | data/paper3_figs.txt |

## References

- Cadoni, M., & Tuveri, M. 2019, PRD 100, 024029
- Chae, K.-H., et al. 2021, ApJ 921, 104
- Deser, S., & Levin, O. 1997, Class. Quantum Grav. 14, L163
- Desmond, H., Hees, A., & Famaey, B. 2024, MNRAS (arXiv:2401.04796)
- Famaey, B., & Durakovic, A. 2025 (arXiv:2501.17006)
- Gibbons, G. W., & Hawking, S. W. 1977, PRD 15, 2738
- Hájek, F. 2026a, companion paper (wide binaries), this repository
- Hájek, F. 2026b, companion paper (coefficients), this repository
- Jaynes, E. T., & Cummings, F. W. 1963, Proc. IEEE 51, 89
- McGaugh, S. S., Lelli, F., & Schombert, J. M. 2016, PRL 117, 201101
- Milgrom, M. 1983, ApJ 270, 365
- Milgrom, M. 1999, Phys. Lett. A 253, 273
- Pöschl, G., & Teller, E. 1933, Z. Phys. 83, 143
- Skordis, C., & Złośnik, T. 2021, PRL 127, 161302
- Unruh, W. G. 1976, PRD 14, 870
- Verlinde, E. 2017, SciPost Phys. 2, 016
