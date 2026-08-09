# Horizon-Coupled Inertia: Research Notes

*Second arc of the PhysicsResearch sessions (2026-07-20/21). Status: parked, awaiting data.*
*Calculations: [calcs/one_dial.py](calcs/one_dial.py), [calcs/toy_inertia.py](calcs/toy_inertia.py),
[calcs/predictors.py](calcs/predictors.py), [calcs/run_to_limits.py](calcs/run_to_limits.py),
[calcs/coincidences.py](calcs/coincidences.py), [calcs/hypothesis_audit.py](calcs/hypothesis_audit.py).*

## The claim under investigation

The standing anomalies of gravity (rotation curves, Λ, why-now, Hubble tension) share one
suspect assumption: **local dynamics is blind to cosmic time.** The empirical anchors:

- a₀ ≈ 1.2×10⁻¹⁰ m/s² (rotation curves) = cH₀/2π to 10%. Equivalently a₀ ≈ c/t_universe.
- Baryonic Tully–Fisher: v⁴ = G·a₀·M_b, ~zero scatter over 5 decades — Mercury-perihelion-grade.
- Freeman's law: disc surface densities saturate at a₀/(2πG) = 0.29 kg/m² (exact match).
- ρ_Λ ≈ 0.7 × 3H₀²c²/8πG — the "constant" sits at today's value of a time-dependent scale.
- The Sun's galactocentric acceleration is 1.8 a₀ — we live in the transition zone.

## The toy and its verdicts

Postulate: inertia = reaction to Unruh-temperature excess over the de Sitter floor
(Deser–Levin): F = m[√(a² + a_c²) − a_c].

1. Functional form: reproduces flat curves; BTFR becomes an identity; tracks the empirical
   RAR to ≤12%. Coefficient: naive a_c = cH₀ overshoots 11×; data demand ~cH₀/4π
   (Verlinde's entropy route lands closest, cH₀/6).
2. **Falsified at Saturn**: the 1/x tail leaves a constant ~a_c residual; Cassini-class
   ephemerides exclude it by 2–3 orders. Real mechanism must recover Newton exponentially
   (phase transition/screening, not a thermal tail). Caveat: linear-acceleration Unruh
   formula applied to circular orbits is itself suspect (circular Unruh differs).
3. LISA Pathfinder (proper accel ~10⁻¹⁴ below a₀, behaved conventionally) kills
   proper-acceleration triggers; survivors reference the cosmic/galactic frame.

## Predictor battery (discriminates H-coupled vs Λ-coupled vs superfluid-DM vs ΛCDM)

| Test | Status |
|---|---|
| High-z Tully–Fisher zero point (a₀ ∝ H(z) → +31% v_flat at z=2; Λ-coupled → 0%) | open; needs outer-disc cold-gas curves (ALMA/SKA), current inner-disc data inconclusive (probe 3–30 a₀) |
| Wide binaries (~10–15% velocity excess with EFE) | Gaia dispute ongoing |
| External field effect: identical galaxies rotate slower in crowded environments (SEP violation; no DM theory mimics it) | Chae+ claim detection in SPARC |
| Tidal dwarf galaxies (DM-free in ΛCDM: 29 vs 64 km/s for 10⁹ M_sun at 5 kpc) | few known TDGs favor boost; JWST can grow sample |
| Lensing-vs-kinematics RAR split (kills or crowns superfluid DM — phonons don't pull light) | lensing RAR tracks kinematic so far |
| Ignatiev SHLEM: ~cm spot, ~0.5 ms, 2/yr, ~80° lat; displacement ~2×10⁻¹⁷ m (PRL 98, 101101 (2007)) | never attempted |
| Saddle-point flythrough (Earth–Sun saddle 259,000 km out; bubble 2.6 m, halo ~100s km; Bekenstein–Magueijo 2006) | LPF flyby proposed, never funded; cubesat + cold-atom gradiometer would do |
| Planet-9 zone: EFE perturbations 0.5–2% at 500–1000 AU shepherd TNOs; test = Rubin finds no planet but clustering aligns with galactic field | open |
| Big-G metrology: sweep test-body acceleration within one apparatus, look for trend in extracted G | never done as protocol |
| Deep-space navigation: unmodeled anomaly at 550 AU ⇒ ~27,000 km drift over 20 yr for a solar-lens mission | future |

## Underlying-physics ladder

Entropic/emergent (Jacobson→Verlinde; couples to Λ, no z-drift) · H-coupled horizon inertia
(our toy; full z-drift; needs exponential screening) · superfluid dark matter
(Berezhiani–Khoury; hybrid, eats the Bullet Cluster; killed by lensing/kinematics identity)
· vacuum-drag inertia (Rueda–Haisch/QI; cornered by MICROSCOPE + binding-energy universality).

## Nearby-space model results ([calcs/space_model.py](calcs/space_model.py), corrected EFE recipe)

- Solar system acquires an **edge at ~5,240 AU** where the galactic field overtakes the Sun's;
  the EFE floor (1.8 a₀) means nearby space is never deep-MOND — only anisotropic transition zone.
- **Sedna-like TNOs** (a=500 AU, e=0.7): apsidal drift ~−0.013 to −0.019°/orbit with a sin(2ω)-shaped
  orientation dependence (textbook quadrupole/Stark signature, max at 45° to the galactic field).
  Differential drift rate ⇒ steady-state apsidal clustering ∝ 1/|dω/dt| with NO planet — clustering
  axis tied to the galactic field, the discriminator vs Planet Nine.
- **Wide binaries** (10 kAU): strongly non-Keplerian clocks and precessing rosettes (~−3.8°/orbit
  in-sim; magnitude IC-dependent, standard statistic is ~15–20% velocity excess). NEW statistical
  predictor: wide-binary apsidal orientations should be non-uniform in galactic coordinates — Gaia
  has the data; test apparently never run.
- **Saddle flyby, honest version**: for the Cassini-surviving exponential family, the anomaly is
  confined to ~100–150 m around the saddle (~7×10⁻¹¹ m/s², LPF-detectable) — a hard-navigation,
  small-corridor mission; Cassini already excluded the families with big saddle halos.
- Methodological note: the first model version had a spurious constant EFE force (bad frame
  correction) that ejected a bound binary and would have re-violated Cassini — caught and fixed
  by exactly those two symptoms. Algebraic EFE recipes are treacherous; full solvers needed for
  publication-grade numbers.

## Gaia dig, first pass (2026-07-21, [calcs/gaia_widebinary_orientation.py](calcs/gaia_widebinary_orientation.py))

Ran the apsidal-orientation test on the El-Badry & Rix 2018 catalog (55,128 pairs, d<200 pc,
downloaded via VizieR TAP to [data/](data/)). Axial Rayleigh statistic of separation-vector
position angles vs the projected galactic-field direction, weighted by in-plane field fraction,
in four separation bins with 100–3,000 AU as selection-systematics controls.
**Result: clean null in all bins** (R = 0.005–0.015, all p > 0.15; v–r angle flat at ~50°).
Sensitivity: alignment amplitude >~2% (2σ) excluded at 8–50 kAU. Constrains strongly-librating
modified-dynamics scenarios; converting to a force-law bound needs population synthesis
(precession-rate modulation + 230 Myr field rotation vs ~66 Myr precession period).
Next step: same pipeline on the El-Badry 2021 EDR3 catalog (~1M pairs) → per-mille sensitivity.

## Interactive sim

[sim/orbit-lab.html](sim/orbit-lab.html) — 3D orbit lab (published as a private artifact):
modified vs Newtonian ghost orbits, EFE vector with adjustable strength/tilt/rotation period,
RAR vs toy vs Newton force families, live precession/period readouts and apsis strip-chart.

## Stage 1: measured screening index ([calcs/sparc_rar_fit.py](calcs/sparc_rar_fit.py), 2026-07-21)

Rebuilt the radial-acceleration relation from raw SPARC data (Zenodo mirror; 153 galaxies
after i≥30°, Q≤2, 10% velocity cuts → 2,700 points; standard M/L 0.5/0.7). Scatter 0.133 dex —
matches the published RAR. Introduced a one-parameter **screening family**
ν_p(y) = (1−e^(−y^p))^(−1/2p) (p=1/2 ≡ the RAR function; smaller p = softer tail):

- **p = 0.443 +0.063/−0.050**, a₀ = (1.03 ± 0.13)×10⁻¹⁰ m/s² (galaxy-level bootstrap).
- Cassini (anomalous Saturn accel < 2×10⁻¹⁴ m/s²) independently requires **p > 0.234**.
- **Coherence finding:** galaxies (10⁻¹²–10⁻⁹ m/s²) and Saturn radio tracking (6.5×10⁻⁵ m/s²),
  eight decades apart, select overlapping regions of the same parameter — no tuning.
- The instanton form p = 1/2 sits 1.1σ from best fit — consistent.
- SPARC alone cannot separate RAR from "simple ν" (ΔAIC 0.4) — but simple ν has a 1/y tail
  that Cassini kills; the two datasets are complementary within one family.
- Soft/algebraic families rejected by galaxies alone: standard ν ΔAIC +245, toy +80.
- Caveats: fixed M/L, no distance/inclination error propagation (matches original RAR
  methodology); marginalizing M/L would widen the p interval somewhat.

Next: Stage 2 — GPU forward-model of wide-binary populations (velocity excess + orientation
amplitude prediction vs our 2% null); optional PySR free-form symbolic regression as a check
that no non-nested family beats ν_p.

## Stage 2A: orientation-verdict forward model (2026-07-21)

Scripts: [calcs/stage2a_orientation_verdict.py](calcs/stage2a_orientation_verdict.py),
[calcs/stage2a_validate.py](calcs/stage2a_validate.py). Agent cost: one Haiku literature scout
(ṽ conventions; Chae 2023: boost 1.43±0.06 at 10^−10.15 m/s², 10σ claim; Banik 2024: 19σ
Newtonian claim; Hwang+22 superthermal p(e)∝e^1.25 at >1 kAU; g_ext ≈ 1.9 a₀ incl. vertical).

- **Methodology bug caught by fingerprint:** first run predicted R = 0.6364 ≈ 2/π exactly —
  the signature of double-folding an axial angle (arccos|cos| then ×2). Fixed with signed
  sky-plane angles. Post-fix prediction under fast-circulation washout: **R ≈ 0.4%**, below
  our 2.2% (2σ) Gaia bound.
- **Real finding (validated at 2× resolution, Newtonian control = 0 exactly):** under the
  standard EFE recipe with our measured law (p=0.443, a₀=1.03×10⁻¹⁰), wide-binary orbits
  beyond ~8 kAU are NOT quasi-Keplerian rosettes. Apsidal drift depends violently on
  orientation to the galactic field: −0.7°/orbit aligned, −0.2 at 45°, **−40°/orbit
  perpendicular, with aphelion pumped 15→29 kAU** (EFE-driven Kozai-like cycles). Secular
  "A + B·sin2ψ" theory is inadequate; several grid cells even circularize (apsis undefined).
- Consequences: (1) the orientation-anisotropy prediction is uncertain between ~0.4% (washout)
  and a few % (rate-contrast clustering) — **the EDR3 million-pair test (~0.3–0.5% sensitivity)
  is decisive either way**; (2) at the widest separations the ṽ statistics in the Chae–Banik
  fight sit on top of these EFE cycles — full orbit-population modeling (GPU) is not optional.
- Next (Stage 2B, GPU): evolve ~10⁶ binaries over Gyr under the recipe (torch/CUDA on the 5090),
  extract both the stationary orientation statistic and the ṽ distributions; download EDR3
  catalog (~1M pairs) and run both measurements at full sensitivity.

## Stage 2C: our independent EDR3 ṽ measurement (2026-07-21, [calcs/stage2c_vtilde_data.py](calcs/stage2c_vtilde_data.py))

El-Badry+21 catalog (1.82M pairs, Zenodo) → 14,071 after strict cuts (R_chance<0.01, d<200 pc,
plx SNR>20, MS band, σ_v<0.03 km/s). Median ṽ vs separation: 0.551 / 0.557 / 0.580 / 0.595 /
0.667 (0.2→50 kAU) — monotonic rise. **Boost ratio (6–30 kAU vs 0.2–2 kAU): 1.086, 68% CI
1.064–1.110** (~3.7σ above Newton). Measured noise inflation: only ~1% in the 6–20 kAU bin
(3% at 20–50) → noise-corrected boost ≈ 1.075 ± 0.023.

**Interpretation: the raw statistic sits in no-man's land** — 3σ above Newton (contradicting a
naive "nothing there"), but ~half the AQUAL-EFE prediction of ~1.20 (Chae–Milgrom 2022 boost
1.37–1.43 in acceleration; our Haiku scout retrieved their fitting function:
g/g_N = ν(y_β)[1+tanh((1.1e_N/y)^1.2)·ν̂(y_β)/3], ν=simple, y_β=√(y²+(1.1e_N)²)).
The tie-breaker is hierarchical-triple contamination (inflates ṽ via unseen mass + PM wobble):
exactly where Chae (boost survives deconvolution) and Banik (triples absorb it) diverge.
Our raw number explains how both camps read the same sky differently.

Next (Stage 2D): triple-star forward model in the GPU engine (triple fraction, inner-orbit
PM wobble, photometric mass bias) → does 1.075 deconvolve toward 1.0 or 1.2? Plus the
dr2_radial_velocity subset enables a 3D-velocity cross-check. C&M22 AQUAL formula ready for
the kernel; orientation channel awaits its anisotropic (curl) part.

## Stage 2D: the tie-breaker runs (2026-07-21)

**RUWE diagnostic ([calcs/stage2d_ruwe_variant.py](calcs/stage2d_ruwe_variant.py)):** the data
boost is RUWE-stable — 1.086 / 1.088 / 1.087 / 1.069 as cuts tighten to RUWE<1.1 (halving the
sample). Short-period wobble triples are NOT driving the signal. (RUWE is blind to inner
companions at ~2–100 AU, so long-period triples remain viable suspects.)

**Triple forward model ([calcs/stage2d_forward.py](calcs/stage2d_forward.py)):** GPU populations
under Newton and the Chae–Milgrom AQUAL formula, with RUWE-invisible triples injected
(a_in ∈ 2–100 AU, photocenter wobble + hidden-mass bias):

| f_triple | Newton boost | AQUAL boost |
|---|---|---|
| 0.00 | 0.965 | 1.316 |
| 0.15 | 0.987 | 1.307 |
| 0.30 | 1.039 | 1.316 |

DATA: 1.086 (68% 1.065–1.105). **Both pure hypotheses fail at face value:** Newton+triples
needs f_t ≳ 0.45 (implausible after quality cuts; multiplicity surveys give ~15–25% total);
full AQUAL-EFE overshoots the data by ~10σ in-model (medians are triple-robust — the wobble
makes outliers, not median shifts; only the hidden-mass bias moves medians, weakly).

**Reading:** the measured anomaly is real but ~half-strength relative to AQUAL — consistent
with (a) a partially-screened/intermediate law, or (b) unmodeled astrophysics (empirical
e(s) trend, IC self-consistency, selection edges — our Newton baseline is 0.965 not 1.000
from projection/truncation effects, showing these matter at the few-% level).
Refinements queued: empirical α(s) eccentricity, self-consistent ICs in the modified
potential, distribution-shape (not median) fits — tails separate triples from boosts —
and the dr2 radial-velocity 3D subset.

## Stage 2E: hardened comparison — the verdict moved (2026-07-21, [calcs/stage2e_refined.py](calcs/stage2e_refined.py))

Upgrades: empirical α(s) eccentricity law (Hwang+22); **self-consistent ICs in the modified
potential** (energy+angular-momentum quadrature; apoapsis fidelity check = 1.000); per-bin
distribution shapes. Findings:

1. **The IC bias was the big one:** with consistent ICs the AQUAL prediction in the 6–20 kAU
   bin drops from 1.32 to ~1.16. Stage 2D's "half-strength anomaly" was partly our own artifact.
2. **Median profile (normalized to 0.2–2 kAU):** DATA 1.045 / 1.073 / 1.202 across 2–6 / 6–20 /
   20–50 kAU. AQUAL(f_t=0–0.1): 1.03–1.05 / 1.16–1.18 / 1.15–1.21 — matches at 2–6 and 20–50,
   data ~8% (≈4σ) below AQUAL at 6–20. Newton: 0.95–0.97 declining — excluded at ~5σ in the
   wide bins (26% gap at 20–50 kAU).
3. **Tail statistic is decisive against triples:** DATA P90/P50 stays flat at 1.83→1.95 across
   all bins. Injected wobble triples explode the tails (f_t=0.2 → P90/P50 of 5–16). The data's
   thin tails cap effective RUWE-invisible triple contamination at ≲5%, killing the Newton
   rescue (which needs ~45%) independently of the RUWE test.
4. Remaining wiggle: photometric-mass systematics, chance alignments in the N=214 widest bin,
   crude wobble amplitudes, the 6–20 kAU dip (noise wiggle or transition-shape information?).

**Updated reading: our independent pipeline now favors a real, AQUAL-consistent velocity
excess in wide binaries, with triple contamination bounded small by distribution shape.**
The 6–20 kAU deficit vs AQUAL is the one surviving tension — worth chasing: it is where a
sharper screening function (our Stage-1 p≈0.44 vs the soft simple-ν used in C&M) would
naturally bite hardest. Next: swap the C&M interpolating function for our measured ν_p family
and see if the 6–20 kAU dip is *predicted*.

## Stage 2F: the consilience test (2026-07-21, [calcs/stage2f_nup_test.py](calcs/stage2f_nup_test.py))

Three self-consistent GPU populations vs data, per-bin normalized medians:

| bin [kAU] | DATA (68%) | Newton | AQUAL-CM (canonical) | ν_p-EFE (our law) |
|---|---|---|---|---|
| 2–6 | 1.045 (1.028–1.063) | 0.967 | 1.030 | 1.026 |
| 6–20 | 1.073 (1.036–1.091) | 0.961 | 1.097 | 1.113 |
| 20–50 | 1.202 (1.161–1.277) | 0.951 | 1.158 | 1.172 |

χ² (3 wide bins): **Newton 55 (≈7σ excluded) · AQUAL-CM 2.0 · ν_p-EFE 3.6.**

- Honest correction: Stage 2E's "6–20 kAU dip below AQUAL" was mostly an a₀-convention
  artifact (mode 2 had run with our a₀ instead of canonical). With conventions clean, both
  modified laws fit the profile well; the residual 6–20 offset is ~1σ.
- The consilience question answered: at current precision the wide-binary channel CANNOT
  distinguish our sharper screening family from simple-ν (Δχ²≈1.6) — both are galaxy-
  calibrated and both fit. Discriminating them needs ~±1.5% precision in the 6–20 kAU bin
  (Gaia DR4-scale, or RV-enhanced samples).
- What stands after the full chain (SPARC→Cassini→GPU→EDR3→tails→consilience): a ~7σ
  preference for modified dynamics over Newton in our own end-to-end pipeline, with triples
  capped ≲5% by distribution shape and noise ~1%. Remaining systematics: photometric masses,
  e(s) law, selection function, chance alignments at 20–50 kAU; and the unresolved conflict
  with Banik+24's Newtonian preference (different sample/statistic) — reconciling the two
  methodologies on one sample is the outstanding community question our pipeline could tackle.

## Stage 2G: the Bose–Einstein identity and its first test (2026-07-22, [calcs/stage2g_be_efe.py](calcs/stage2g_be_efe.py))

**The identity (exact, verified symbolically):** the p=1/2 interpolating function that
galaxies select is ν(y) = 1/(1−e^(−√y)) ≡ **1 + n_BE(x)** with x = √(g_N/a₀) = r_M/r,
where n_BE = 1/(e^x −1) is a Bose–Einstein occupation. Consequences:
- Deep MOND = the **Rayleigh–Jeans limit** (n → 1/x → g = √(g_N a₀), scale invariance =
  classical equipartition); Newtonian regime = the **Wien tail** (mode frozen out);
  the RAR is a Planck curve in the variable r_M/r.
- Series prediction at NLO: g = √(g_N a₀) + g_N/2 + ... (coefficient 1/2, testable in
  principle in RAR curvature).
- Mode energy ε ∝ √(g_N a₀) = ħ/2πc × (deep-MOND acceleration's Unruh quantum); bath
  temperature = T(a₀), and a₀ = cH₀/2π ties the bath to the horizon — the 2π bookkeeping
  from the first toy session reappears exactly here. **If the bath is the current horizon,
  a₀ ∝ H(z): the BE reading selects the z-drift branch — high-z Tully–Fisher kills or
  crowns it.**
- Our Stage-1 measurement p = 0.443 ± 0.06 is 1.1σ from the BE-required p = 1/2.

**First confrontation:** naive BE-EFE embedding (external field inside the mode energy,
quadrature): wide-binary profile 1.030/1.140/1.152 → χ² = 7.2 vs AQUAL-CM 2.0, Newton 53.5.
The data now discriminate EFE *embeddings* (Δχ² ≈ 5) even where ν-families were degenerate.
Naive BE embedding mildly disfavored (~2σ) → the BE reading needs its own two-field theory,
not a quadrature shortcut.

**Priority sweep verdict (2026-07-22, two Haiku agents):** the occupation identity is NOT
explicitly published. Milgrom: phenomenological; McGaugh+16: empirical, no thermal reading;
Desmond+23 exhaustive symbolic regression: no BE remark; nearest neighbors to cite: Pazy &
Argaman (PRD 85, 104021 — quantum statistics + freeze-out on holographic screens, different
function), q-deformed heat-capacity MOND (arXiv:2010.03530), Smolin, Ho–Minic–Ng. Caveat:
sweep is Haiku-grade, not exhaustive — do a citation-graph walk of McGaugh+16 before
submission; frame as "we note that…" regardless. Section-2 draft text now in
[PAPER-DRAFT.md](PAPER-DRAFT.md).

**Theory-paper skeleton now on the table:** (1) the exact identity + Planck-curve reading;
(2) the measured screening index as its test (p = 1/2 within 1.1σ); (3) the SPARC+Cassini
joint constraint; (4) the wide-binary EFE-embedding discrimination; (5) predictions: NLO
coefficient 1/2, a₀ ∝ H(z), occupation-fluctuation noise (uncomputed).

## Stage 2H: the true BE-EFE, and a methodological discovery (2026-07-22)

**QUMOND solver built and validated** ([calcs/qumond_efe_solver.py](calcs/qumond_efe_solver.py)):
point mass in uniform external field, stable multipole solve. Gates: spherical identity
reproduced to 0.01%; simple-ν tracks the C&M formula to 1–2% at y≥3 (10% gap at y≲1 =
AQUAL-vs-QUMOND formulation difference, expected). One sign bug caught (phantom source built
from +g_N instead of ∇φ_N → repulsive phantom, boost<1 — the fingerprint).

**Theory result:** the Bose–Einstein ν's true (QUMOND-solved, sphericalized) EFE is ~4%
weaker than simple-ν's in the wide-binary window and has a Wien-fast tail (boost dies by
y~30 where simple-ν keeps a soft 1/y floor). Tables saved: [data/efe_boost_*.npy](data/).

**Population verdicts** ([calcs/stage2h_true_be.py](calcs/stage2h_true_be.py) + controls):
Newton 53.5 / AQUAL-CM(formula) 2.0 / BE-true(table) 10.1 / simple-ν-QUMOND(table) 17.9.
**BUT: interpretation suspended** — the isolation experiment
([calcs/stage2h_isolate.py](calcs/stage2h_isolate.py)) showed the median-ṽ statistic
responds sub-linearly to force boosts under (rp,ra)-preserving self-consistent ICs: a
designed 12.8% velocity boost (verified at fixed radius) appears as only ~6% in median
v/v_K, because boosted orbits redistribute time-weighting across radii. Additionally the
normalization bin is contaminated differently by soft-tail vs Wien-tail functions. So the
χ² ranking mixes: EFE saturation level, tail softness (control-bin normalization), IC
conventions, and statistic response — each at the several-% level, i.e. the same order as
the model differences. **The BE reading is NOT yet fairly disfavored; nor is AQUAL fairly
crowned.** This sensitivity applies equally to published forward models in the
Chae–Banik debate — a paper-worthy methodological point on its own.

**Correct next step:** distribution-level (not median) likelihood with IC-convention
systematics marginalized, on one sample, for all models — the reconciliation framework.

## Stage 3A: distribution-level likelihoods + the field's undocumented systematic (2026-07-22)

Engine: [calcs/stage3a_likelihood.py](calcs/stage3a_likelihood.py) — full ṽ-distribution
log-likelihood per bin, model populations noise-convolved with the data's own per-bin error
distribution; each modified law run under two IC conventions (A: self-consistent (rp,ra);
B: Keplerian-then-relax). Results (ΔlnL vs Newton; higher better):

| model | ΔlnL |  | model | ΔlnL |
|---|---|---|---|---|
| aqual-A | **+730** | | aqual-B | −27 |
| BE-A | **+455** | | BE-B | +75 |

- **Modified gravity with self-consistent populations beats Newton decisively (+455 to +730).**
- **But the IC/realization spread (A−B: 757 and 380) EXCEEDS the model-to-model differences
  (aqual-A vs BE-A: 275) and even the Newton-vs-modified gap.** Caveat: convention B is
  partly a straw man — it fails to preserve the intended e-distribution (orbits circularize),
  so A−B overstates the legitimate band; the fair statement is that the *realized orbit
  population* must be controlled/marginalized, and its influence is first-order.
- Suspicious detail to chase: models differ by +143 lnL in the CONTROL bin (y≫1), where law
  differences should vanish — soft-tail control-bin contamination and/or data spread
  (triples) mimicry. Needs a dedicated look before any lnL number is quoted externally.

**Methodology scout (Haiku) on the published camps:** Chae = binned medians of orthogonal
residuals, mocks sampled in the NEWTONIAN potential even when testing MOND; Banik = binomial
likelihood in (r_sky, ṽ) pixels with a MOND-integrated orbit library + nuisance
marginalization (eccentricity, triples). **Confirmed critical gap: no published work
quantifies the sensitivity of these statistics to orbit-population/IC realization within a
fixed force law — the systematic our engine just measured at ΔlnL ~ hundreds.** The
Chae–Banik disagreement plausibly lives inside it (arXiv:2312.03162 reviews the dispute).

**Next milestone (the definitive machine):** hierarchical inference — fit the gravity-law
parameter JOINTLY with the orbit-population priors (e-distribution, a-distribution, triple
fraction) at distribution level, realization-consistent per law. Banik marginalizes some of
this; nobody marginalizes the realization procedure itself. That is the reconciliation paper.

## Stage 3B: hierarchical inference v1 — first results (2026-07-22, [calcs/stage3b_hierarchical.py](calcs/stage3b_hierarchical.py))

Joint fit (α_grav × η_ecc × f_triple), populations re-realized self-consistently per grid
point, full-distribution multinomial likelihood, 14,071 pairs. Cubes: [data/hier_cube_*.npy](data/).

| law | α̂ | η̂ | f̂_t | ΔlnL(α̂ vs Newton) |
|---|---|---|---|---|
| simple (QUMOND) | 1.25 (GRID EDGE) | 2.0 (edge) | 0.05 | **+243.8** |
| Bose–Einstein | **1.00 (interior!)** | 2.0 (edge) | 0.05 | **+220.3** |

Findings:
1. **Newton remains rejected (ΔlnL ≈ +220–244) even with population and triple nuisances
   free** — the anomaly is NOT absorbed by the realization systematic when that systematic
   is marginalized. (Naive Wilks ~20σ; do not quote σ literally — histogram-floor and
   control-bin caveats pending.)
2. **The BE law peaks exactly at its parameter-free prediction α = 1.0** (over-strength
   α=1.25 disfavored by 14.7); the simple law runs off the grid edge wanting >1.25× its own
   strength. However simple-at-edge still edges BE-at-1.0 by ~24 lnL total — no clean shape
   verdict until the α grid is extended.
3. Both laws push η to the 2.0 edge — the fit wants eccentricities more superthermal than
   Hwang's 1.3 at the widest separations (either real, untested regime >10 kAU — or misfit
   absorption). Extend grid; check against Hwang's actual validity range.
4. f_t = 0.05 preferred — consistent with our independent tail bound.

Refinement queue before quoting externally: extend α (→2.0) and η (→3.0) grids; resolve the
Stage-3A control-bin anomaly; bootstrap lnL noise (N=500k MC ±few units); chance-alignment
stress test in the widest bin; mass-model systematic.

### Stage 3B v2: extended grids (2026-07-22, later) — RETRACTION + diagnosis

With α ≤ 2.0 and η ≤ 3.0, BOTH laws run to the new corner (α=2.0, η=3.0, f_t=0.05):
simple ΔlnL(Newton) = +338.2; BE = +306.6; BE α-profile non-monotonic (local max at α=1.0,
dip at 1.25, rise to the edge).

- **RETRACTED: the v1 "BE bullseye at α=1.0" was a grid-boundary artifact.** With wider
  fences the maximum migrates along an (α,η) degeneracy ridge — higher e-populations lower
  typical ṽ (apoapsis dwell), higher α raises it; the two trade off and the fit climbs to
  implausible corners (p(e)∝e³, double-strength gravity).
- **ROBUST across everything so far: Newton's rejection (+306 to +338, growing with grid
  size), and f_t = 0.05.**
- Diagnosis: corner-seeking = missing model ingredient absorbing into (α,η). Prime suspects:
  (1) chance-alignment/unbound contamination at wide s (high-ṽ tail) — needs its own
  nuisance component (Banik has one; we don't yet); (2) flat η prior is wrong — Hwang's
  measurement warrants η ~ N(1.3, 0.3) as a Bayesian prior; (3) the 3A control-bin shape
  issue. Expect the ridge to collapse and α to localize once (1)+(2) are in.
- Meta-note: this happened AFTER the git commit and BEFORE any public claim — the discipline
  working as intended.

## Stage 3C: v3 fit — contamination + η prior did NOT collapse the ridge (2026-07-22, [calcs/stage3c_v3fit.py](calcs/stage3c_v3fit.py))

TODO #1 implemented as designed: chance-alignment component flat in ṽ with weight
f_c0·(s_b/31.6 kAU)², Hwang-anchored prior η ~ N(1.3, 0.3), denser η grid (0.8–1.9),
α ≤ 2.0. 119 GPU population realizations, ~13 min total. Cubes (prior included):
[data/hier3_cube_*.npy](data/), summary: [data/stage3c_summary.txt](data/stage3c_summary.txt).

**Result: the collapse hypothesis of Stage 3B v2 is FALSIFIED.** Both laws again maximize
at a corner — (α=2.0, η=1.9, f_t=0.05, f_c0=0.1), interior max = False for both, α "1σ
range" degenerate at the fence. Newton rejection: simple +296.4, BE +263.8 (fourth
consecutive configuration in the +220…+340 band).

Cube diagnostics (why the fix failed):
1. **The η prior is powerless against the likelihood slope.** lnL still gains ≈ +17–21 per
   grid step from η=1.6→1.9 while the prior charges only −2.0 at η=1.9. A N(1.3,0.3) prior
   cannot hold a ridge this steep; the fit's preference for superthermal e at wide s is a
   likelihood feature, not a flat-prior artifact.
2. **Contamination is real but small: +8.3 (simple) / +12.0 (BE) lnL total** — it soaks a
   little tail but was never the ~300-lnL missing ingredient. (f_c0 sits at its 0.1 fence;
   extending it is cosmetic at this gain scale.)
3. **The killer: α runs to the 2.0 corner even with η pinned at 1.3** (conditional profile
   monotonic, gaining ≈ +89 (simple) / +51 (BE) from α=1→2). So this is NOT just an (α,η)
   trade — the fit wants ever-more boost (or is absorbing a shape misfit through it)
   independent of the eccentricity population. Double-strength MOND out-fitting MOND is a
   misspecification fingerprint, not physics.
4. f_t = 0.05 again sharply preferred (hundreds of lnL penalty either side) — the most
   stable nuisance in the problem.

Interpretation: the missing ingredient is upstream of the nuisance layer — prime suspect is
the Stage-3A control-bin anomaly (~140 lnL of model discrimination generated in the 0.2–2
kAU bin where the force effect is per-mille; the likelihood is evidently not shape-faithful
in the Newtonian regime, and α gains cannot be read as gravity until it is).
**Consequence for the paper: there is no α ± interval headline. Do NOT extend the α grid
further (v2's lesson: chasing the ridge outward just moves the corner).** The paper leads
with what is robust: the 1.086 boost measurement, Newton rejection with nuisances free, the
realization systematic, and the BE identity as theory motivation.

Next (blocking): TODO #2 control-bin autopsy, with one instrumentation change — modify the
fit to save per-bin lnL contributions so we can localize which s-bins generate the α
preference before touching any more nuisances.

## Stage 3D: control-bin autopsy (2026-07-22, same night, [calcs/stage3d_autopsy.py](calcs/stage3d_autopsy.py))

Five models on the identical seed-31 realization (Newton, simple/BE at α=1, simple/BE at
α=2; η=1.3, f_t=0.05, f_c0=0.1), decomposed per s-bin and per ṽ-region
(R0: ṽ<0.3, R1: 0.3–0.8, R2: 0.8–1.5, R3: >1.5). Full tables:
[data/stage3d_autopsy.txt](data/stage3d_autopsy.txt). Two findings, one good, one bad:

1. **The Newton rejection is generated where physics says it should be**: the 2–6 and
   6–20 kAU bins, R2 region. At α=1 both laws reproduce the R2 data fraction almost
   exactly where Newton fails badly (6–20 kAU: data 0.251, Newton 0.158, simple 0.238,
   BE 0.260). The rejection is NOT a far-tail artifact (R3 contributes little at α=1).
2. **The α>1 preference is NOT generated by the signature — it is misfit compensation.**
   In EVERY bin (including the Newtonian control bin) the data ṽ distribution is BROADER
   than every model: excess at R0 and R2, deficit at R1 and R3. All models sit at
   χ²/dof ≈ 15–19 in the control bin — no α cures it; raising α just shifts the model peak
   rightward (R1→R2), buying lnL against the baseline deficit while worsening R0. The
   Stage-3A "~140 lnL control-bin discrimination" largely evaporates at η=1.3 (→ +18–31
   at α=1): it was mostly the η-ridge in action.

Diagnosis: the model ṽ distribution is too narrow. Prime physical suspect: the data's ṽ is
normalized by error-laden photometric masses (multiplicative smear ≈ ½·δM/M), the model's
by its own true masses — the model is missing a broadening convolution. This is the engine
under the (α,η) ridge: superthermal η was a broadening proxy.

## Stage 3E: the smear test — ridge broken, rejection deflated (2026-07-22, same night, [calcs/stage3e_smear.py](calcs/stage3e_smear.py))

Multiplicative model-ṽ smear exp(σ_m·g) added at evaluation level, σ_m ∈ [0…0.40], same
five models/realization. Results ([data/stage3e_smear.txt](data/stage3e_smear.txt)):

- **σ_m localizes interior: ≈ 0.20–0.25 for every model** (lnL turnover confirmed; Newton
  gains up to +394, α=1 laws +277 at their optima). The broadening is real and enormous.
- **The α corner-seeking collapses.** (α=2 − α=1) lnL: simple +89 → −23, BE +51 → −22 at
  each law's preferred σ_m. With broadening in the model, α≈1 finally beats both Newton
  AND over-strength gravity — the fit behaves like a measurement for the first time.
- **The honest Newton rejection deflates from +296/+264 to ≈ +63/+66** (simple/BE at α=1
  vs Newton, both at their preferred σ_m ≈ 0.25, one realization). Newton benefits MORE
  from broadening than MOND does — a symmetric s-independent smear partially mimics the
  excess — but cannot close the gap because the real signature is s-DEPENDENT. The
  distribution-level rejections quoted in Stages 3A–3C were inflated by this misfit.
- The median-based boost 1.086 (CI 1.064–1.110) is immune by construction (multiplicative
  lognormal smear is median-preserving) — it stands as the model-independent anchor.
- **Honest-updating**: credence in "excess is real physics, not systematics" moves DOWN,
  not up, on this result (~65% → ~60%): a single symmetric nuisance ate 3/4 of the lnL gap,
  and σ_m ≈ 0.25 implies ~50% effective mass error — far beyond the accounted photometric
  budget, so the smear is a proxy for something not yet decomposed (mass errors + e-model
  mismatch + selection effects?). Until σ_m is independently measured (main-sequence width,
  FLAME/spectroscopic masses for a subsample) and PRIORED rather than floated, the α
  measurement inherits the (α,σ_m) partial degeneracy.

Next: v4 fit ([calcs/stage3f_v4fit.py](calcs/stage3f_v4fit.py)) with σ_m as a full nuisance
axis (α × η × σ_m × f_t × f_c0, η prior kept), reporting α both with σ_m free and capped
at 0.15 — running.

## Stage 3F: v4 fit — α localizes interior at ≈1.0 for BOTH laws (2026-07-22, same night, [calcs/stage3f_v4fit.py](calcs/stage3f_v4fit.py))

Full grid α(9) × η(5) × σ_m(6) × f_t(3) × f_c0(3), η prior kept, seed 31. Cubes:
[data/hier4_cube_*.npy](data/), summary: [data/stage3f_summary.txt](data/stage3f_summary.txt).

With σ_m free:
| law | α̂ (parabolic) | grid 1σ | interior? | Newton ΔlnL |
|---|---|---|---|---|
| simple | 0.94 | [0.87, 1.05] | **YES** | +59.6 |
| BE | 1.05 | [0.98, 1.06] | **YES** | +58.1 |

Supporting structure: η peaks interior at 1.05 (prior-consistent — the superthermal-η
pathology is gone, σ_m replaced it, confirming the 3D diagnosis); BE vs simple at each's
best: −1.5 lnL — statistically indistinguishable, no shape verdict; with σ_m capped ≤0.15
α runs back to the 2.0 fence (the (α,σ_m) trade made explicit — the bracket is the honest
statement until σ_m is independently measured).

**Caveats before anyone gets excited (v1-bullseye silhouette — treat accordingly):**
1. **σ_m best = 0.30 = its grid edge** in this cube (the richer nuisance freedom pushed it
   past 3E's 0.25 turnover). Edge-extension run (Stage 3G) required before the result
   stands. [RESOLVED in Stage 3G below: σ_m=0.30 is a clean interior optimum.]
2. **The near-peak profile differences (−1.0, −2.7 lnL) are within single-realization
   noise** — the tight 1σ intervals are NOT yet quotable; MC repeats (TODO #3) must set the
   real interval width.
3. f_t flipped from the previously stable 0.05 to 0.00 (−311 lnL for 0.05!): σ_m and
   triple-wobble broadening are strongly degenerate; the broadening sector's internal
   decomposition is unconstrained by this data. α appears robust to it (marginalized), but
   quote nothing about triples from this fit.
4. Same-realization caveat as always: one seed (31), one noise-pick convention.

Reading: after mechanically breaking the ridge with the broadening nuisance, the fit's
preferred boost strength is the parameter-free prediction α=1 — for both interpolating-
function families, from a distribution-level fit with five nuisances marginalized. This is
what the retracted v1 "bullseye" would have looked like if it had been real. It is now
allowed to be interesting — after 3G and the MC error budget it may be allowed to be true.

## Stage 3G: σ_m edge extension — the localization stands (2026-07-22, same night, [calcs/stage3g_smedge.py](calcs/stage3g_smedge.py))

Extended σ_m grid [0.20…0.45] (α and η grids coarsened for runtime; same seed). Summary:
[data/stage3g_summary.txt](data/stage3g_summary.txt), cubes data/hier4b_cube_*.npy.

- **σ_m = 0.30 is a clean INTERIOR optimum**: profile [−110, −15, 0, −42, −153, −310]
  (simple), [−127, −22, 0, −35, −146, −308] (BE). Falls off steeply on both sides. Caveat 1
  of Stage 3F is closed — the v4 best point is a real maximum, not a fence artifact.
- **The α result is unchanged with the extra σ_m room**: grid best α=1.0 for both laws,
  interior, Newton ΔlnL +59.6/+58.1 — identical profile values to v4 (consistency check on
  the shared realization passes). (The parabolic α̂ shifts to 0.88/1.11 purely because the
  coarser α grid changes the interpolation neighbors — grid artifact, ignore; the statement
  is "α consistent with 1", not a third decimal.)
- Every axis is now interior or at a physical (not artificial) boundary: α=1.0, η=1.05,
  σ_m=0.30, f_t=0.00 (physical floor), f_c0=0.02.

**Standing result of the night (2026-07-22):** the distribution-level hierarchical fit,
with broadening + eccentricity + triples + contamination marginalized and the η prior
active, prefers boost strength α ≈ 1 (the parameter-free MOND-scale prediction) for both
ν-families, and disfavors Newton by ΔlnL ≈ +58–60. Remaining gates before the α interval
is quotable: MC error budget (TODO #3 — near-peak lnL differences are within
single-realization noise), and the physical decomposition/independent measurement of
σ_m ≈ 0.30 (TODO #2b — it implies ~60% effective mass error, which photometry alone cannot
supply; whatever it really is, the fit needs its prior).

Credence update (2026-07-22, post-3G): excess is real physics ~65% (back up from 60: the
α=1 localization with all nuisances free is a genuine positive; the undecomposed σ_m keeps
it from going higher). Anomaly pattern is one law: ~55% (unchanged pending MC + σ_m
physics). BE-vs-simple shape discrimination: none (−1.5 lnL).

## Stage 3H: MC realization budget — localization is seed-robust; the interval is realization-dominated (2026-07-22, same night, [calcs/stage3h_mcbudget.py](calcs/stage3h_mcbudget.py))

Six independent population realizations (seeds 31/101/202/303/404/505), reduced v4 grid
around the optimum. Full lines: [data/stage3h_summary.txt](data/stage3h_summary.txt).

| seed | simple α̂ | BE α̂ | ΔlnL(Newton) s/BE |
|---|---|---|---|
| 31 | 0.94 | 1.05 | +59.6 / +58.1 |
| 101 | 1.13 | 1.27 | +50.5 / +52.4 |
| 202 | 0.83 | 1.47 | +58.6 / +55.1 |
| 303 | 0.84 | 1.45 | +52.2 / +52.9 |
| 404 | 0.94 | 1.24 | +51.2 / +54.8 |
| 505 | 0.87 | 1.30 | +54.7 / +55.6 |

**Aggregates: simple α̂ = 0.93 ± 0.11, BE α̂ = 1.30 ± 0.15 (realization scatter);
interior 12/12; Newton ΔlnL = +54.5 ± 3.9 (simple), +54.8 ± 2.0 (BE).**

Findings:
1. **The interior localization is realization-robust: 12/12 fits, never a corner.** The α
   profile has a genuinely flat top over ≈[0.75, 1.5] with per-seed maxima decided by
   sub-lnL wiggles — the per-seed parabolic "1σ" intervals of Stage 3F (±0.05) are
   confirmed meaningless; the realization systematic dominates the α error budget by ×2–3.
   This closes the loop on the Stage-3A discovery: the realization systematic is now
   measured ON the headline parameter itself.
2. **Newton rejection is the stable number: ΔlnL ≈ +55 ± 4 across all seeds and laws.**
3. **BE α̂ > simple α̂ on every seed** (offsets +0.11…+0.64): systematic, direction
   consistent with the Stage-2G finding that the BE-EFE boost is ~4% weaker than simple-ν
   (a weaker per-α law needs higher α̂ to match the same excess; the boost enters
   nonlinearly so the offset amplifies). Both laws bracket α=1 within ~2σ_realization.
4. σ_m = 0.30 selected by every one of the 12 fits — eerily stable; whatever it is, it is
   a property of the data, not of the realization.
5. Caveat: the reduced η grid's floor (0.8) was active in several fits — mildly subthermal
   η preferred with σ_m free; harmless for α given the flat top, but the full-grid fit
   should extend η below 0.8 once σ_m has a physical prior.

**End-of-night statement (2026-07-22, ~4am):** α = 1 is consistent with both ν-families at
realization-dominated precision ~±0.15–0.2; Newton is disfavored by ΔlnL ≈ +55 ± 4 with
five nuisances marginalized; the two remaining gates on a quotable α ± interval are the
physical decomposition of σ_m ≈ 0.30 (TODO #2b) and the data-bootstrap half of TODO #3.

## Stage 3I: data bootstrap — the α error budget completes (2026-07-22 morning, [calcs/stage3i_bootstrap.py](calcs/stage3i_bootstrap.py))

1000 bootstrap replicates of the 14,071 pairs, re-scored against a saved seed-31 model
grid (resampling at the evaluation layer — no orbit re-runs; model noise-draw convention
held fixed, a stated second-order approximation). [data/stage3i_summary.txt](data/stage3i_summary.txt):

- simple: α̂ = 0.98 ± 0.20 (16–84%: [0.84, 1.15]); interior 984/1000
- BE: α̂ = 1.21 ± 0.26 (16–84%: [1.00, 1.51]); interior 978/1000
- Newton ΔlnL = +60.7 ± 10.6 (simple) / +58.9 ± 10.3 (BE); **minimum over 1000
  replicates: +30** — Newton loses in every resampling of the data.

**Combined error budget (bootstrap ⊕ realization, quadrature): simple α = 0.98 ± 0.23;
BE α = 1.21 ± 0.30.** Both consistent with the parameter-free α = 1 within 1σ. TODO #3 is
done. This is the α ± interval — CONDITIONAL on the σ_m nuisance being benign, which
Stage 3J immediately puts under fire:

## Stage 3J: σ_m mass-budget — mass errors REFUTED as the source; multiplicity promoted (2026-07-22 morning, [calcs/stage3j_massbudget.py](calcs/stage3j_massbudget.py))

Independent measurement of the photometric mass error from the sample's own photometry
(main-sequence ridge in (BP−RP, M_G) from the 2×14,071 component stars; the binary trick:
corr(δ₁,δ₂)=0.47 across pairs splits shared vs per-star scatter). Full output:
[data/stage3j_summary.txt](data/stage3j_summary.txt). (One bug caught: the catalog encodes
missing photometry as 1e20 sentinels, not NaN — first pass gave a 2.86-mag "MS width";
suspicious round-trip numbers remain confessions.)

- MS robust width: **0.275 mag** at fixed color → per-star σ_lnM ≈ **0.056** (5.6% mass
  error — photometric masses are GOOD) → **σ_m(mass) = 0.024 (16–84%: 0.020–0.031)**.
- **The fitted σ_m = 0.30 is 12× larger. Mass errors supply 8% of the variance; the
  residual (0.299) is the entire effect. The Stage-3E "error-laden photometric masses"
  hypothesis is REFUTED** — logged as such; that was this model's best guess at 4am and
  the measurement killed it in 40 lines. Honest-updating cuts both ways.
- The same measurement hands us the replacement suspect: **12.3% of component stars are
  overluminous by >0.4 mag at fixed color = unresolved companions, counted in OUR sample.**

Why this is the critical suspect and not a detail: a hidden companion's wobble is constant
in km/s while vc falls with s, so multiplicity broadening GROWS with separation — the same
direction as the MOND signal (the classic hidden-triple concern, now with a measured
anchor). Our earlier f_t ≲ 5% tail bound used the crude two-population wobble model, which
the σ_m fit rejects in favor of a smooth broadening — a realistic companion CONTINUUM
(lognormal period distribution à la Raghavan, anchored to f_comp = 0.123) may reproduce
the σ_m = 0.30 smear. **Until the v5 fit replaces σ_m with that physical model, α ≈ 1 is
NOT safe — the s-dependence that survived the s-INDEPENDENT smear could be partly eaten by
an s-DEPENDENT one.** This is now the sharpest known threat to the result and the next
mandatory computation (supersedes the generic TODO #2b phrasing).

## Stage 3K: v5 physical-multiplicity fit — α SURVIVES; hidden triples disfavored (2026-07-22 morning, [calcs/stage3k_v5fit.py](calcs/stage3k_v5fit.py))

The decisive test from 3J: replace σ_m entirely with a physical companion continuum
(per-star fraction f_comp gridded 0→0.6; q~U(0.1,1); Raghavan lognormal periods; Kepler
a_in; unresolved <130 AU + stability <a/5 cuts; PM-averaging suppression S=min(1,P/17.8yr);
signed reflex wobble from both components; hidden-mass boost). This channel's broadening
GROWS with s — if multiplicity were faking the signal, f_comp would run high and α → 0.
Summary: [data/stage3k_summary.txt](data/stage3k_summary.txt), cubes data/hier5_cube_*.npy.

| law | α̂ | interior? | Newton ΔlnL | f̂_comp | total lnL vs v4-σ_m |
|---|---|---|---|---|---|
| simple | 1.25 | YES | +119.8 | 0.10 | **−417** |
| BE | 1.75 (flat 1.0–2.0) | YES | +108.1 | 0.10 | **−427** |

Findings:
1. **α does NOT collapse under the s-dependent broadening channel.** It stays interior at
   ≥1 (and drifts UP, because companions under-supply the broadening the data need, mildly
   re-opening the old absorption direction). Newton rejection grows to +108/+120.
2. **The companion model fits ~420 lnL WORSE than the phenomenological σ_m smear** on the
   same data/seed — the wobble+hidden-mass shape (heavy tail, s-growing) is the wrong
   shape; the data want smooth, s-independent, log-symmetric broadening.
3. **The ṽ data cap f_comp at 0.10** (0.20 costs −380, 0.30 costs −920) while the
   photometric census (12.3% overluminous) would need f_c ≈ 1.1 under this model's
   q-uniform assumption — i.e., real companions are mostly in configurations this wobble
   channel suppresses (short-P, PM-averaged), and the velocity data independently forbid
   the rest from being large. **The Clarke-style "hidden triples fake the wide-binary
   anomaly" hypothesis is disfavored by shape, by amplitude, and by α's survival.**
4. Model-crudeness caveats: S(P) kernel is rough, M_h=M_s/2, q-uniform, resolved-companion
   selection approximated. A determined skeptic can tune the population; the reply is that
   the fit was free to use the channel at any strength and declined at −420 lnL.

**σ_m's physical identity remains open.** Refuted: mass errors (3J), companion wobble (3K).
Leading remaining candidate: eccentricity-distribution SHAPE beyond our power-law family
(a smooth ṽ-widener at all s) — measurable directly on the data via Hwang's v-angle method
(TODO #15, now promoted: it closes the last systematic AND is a standalone result).

Credence update (2026-07-22, post-3K): wide-binary excess is real physics: **~70%** (the
canonical killer hypothesis was given a physical channel and the data rejected it).
Low-acceleration anomaly is one law (galaxies+binaries): ~60%. α ≈ 1 (given excess is
real): ~75% within ±0.3. Horizon/BE microphysics specifically: ~15–20% (unchanged — needs
NLO=½ and a₀∝H(z) tests, not more binary fits).

## Stage 3L: v-angle eccentricity measurement — the e-distribution is a MIXTURE; σ_m's identity found? (2026-07-22 morning, [calcs/stage3l_vangle.py](calcs/stage3l_vangle.py))

γ = sky angle between separation and relative-PM vectors, folded [0°,90°] — uses only
DIRECTIONS, hence immune to the mass normalization by construction: the clean channel for
p(e|s). Forward-modeled with our orbit engine (Newton + BE α=1), matched per-pair noise,
S/N>3 on both sides. 13,929 usable pairs. [data/stage3l_summary.txt](data/stage3l_summary.txt).

Gates: G2 mock-recovery PASSES (η=1.3 recovered sharply). G1 required recalibration — my
"circular ⇒ γ≈90°" expectation was WRONG (projection does not preserve perpendicularity;
an edge-on circular orbit gives γ=0). Pure-geometry MC: circular ⇒ 47% in [70°,90°], mean
γ=60° — the engine reproduces it (50%). Gate passes with the corrected target. (Logged
as a lesson: calibrate gates from first principles, not intuition.)

Naive per-bin power-law fits: η̂ = 2.4(EDGE), 2.4/2.1, 1.7–1.8, 0.6–1.0 by s-bin — the
close bin runs to the fence, in tension with Hwang's ~1.3. Histograms show why the
power-law family cannot fit AT ALL:

**The data γ distributions are U-SHAPED** (excess at BOTH γ<20° and γ>80°, deficit in the
middle; strongest in the wide bins — widest: 0.169 at [0,10] and 0.159 at [80,90] vs
model ~0.09/0.11), while EVERY single-power-law model is a smooth gentle tilt. A U-shape
requires a MIXTURE: (a) a near-circular subpopulation (projection piles circular orbits
at 90°: 31% in [80,90], pure geometry); (b) a radial/aligned component — e→1 orbits give
γ≈0 at ANY viewing angle (parallel 3D vectors project parallel), and unbound flybys/
chance-alignments do the same, growing toward wide s exactly like the f_c0 ∝ s²
contamination scaling. The close-bin η̂=2.4 edge was the power-law family chasing the
low-γ arm with its only knob.

**Interpretation: the missing σ_m broadening now has a face.** A circular+radial e-mixture
broadens ṽ smoothly at all s (circular ⇒ ṽ concentrated; radial ⇒ ṽ spread wide; mixture ⇒
broader than any single power law ≈ a multiplicative smear) — and the direction channel
demands the mixture INDEPENDENTLY of everything the previous stages used. All three σ_m
candidates now: mass errors REFUTED (3J), companion wobble REFUTED (3K), e-mixture
SUPPORTED (3L, direction-only evidence).

Also logged: the main pipeline applies sigv per sky component; the data's per-component
error is sigv/√2 — the model has been mildly OVER-noised throughout (conservative: works
AGAINST the too-narrow-model finding, cannot have caused σ_m; fix in v6).

**Next (the definitive fit): v6 = JOINT (ṽ, γ) likelihood** with p(e) = w_circ·(low-e
component) + (1−w)·power-law(α_e), plus the aligned/unbound γ-0 component tied to the
existing f_c0 contamination. The γ data pin the e-mixture; the ṽ data then measure α_grav
with the broadening PHYSICAL instead of phenomenological. If α ≈ 1 survives v6, every
nuisance in the fit has a measured identity.

Priority scouts (two Haiku-grade sweeps, 2026-07-22, 53 searches total): the SPECIFIC
claim — measured U-shaped γ distributions at >6 kAU decomposed as a circular+radial
e-mixture, resolving the ṽ-broadening systematic — appears UNPUBLISHED. But the novelty
is narrower than "e-bimodality is unknown": (1) **El-Badry 2024 review (arXiv:2403.12146)
reportedly uses "U-shaped eccentricity distribution" language for wide binaries** —
VERIFY the exact statement before citing/claiming; if p(e) U-shape is review-level known,
our contribution is the γ-space measurement + the σ_m connection, not the bimodality
itself; (2) the wide-TWIN population is published as extremely eccentric with v-r angle
enhancement (Hwang/Andrews 2022 — scout citations may be confabulated, verify DOIs);
(3) **Pittordis & Sutherland 2023/2025 (arXiv:2504.07569)** fit binary+triple+flyby+
chance mixtures on similar data — the framework prior art; flyby-as-contaminant is
published, flyby-at-γ≈0 as a measured signature is not; (4) all Hwang-method follow-ups
reportedly assume single power-law families. Action unchanged: READ El-Badry 2024 and
P&S 2025 properly before building v6 or claiming anything.

## Stage 3M: v6 joint (ṽ,γ) fit — α robust again; the circular-mixture reading of 3L REFUTED; the broadening's identity narrows further (2026-07-22 night, [calcs/stage3m_v6fit.py](calcs/stage3m_v6fit.py))

First joint 2D likelihood (per s-bin 20×6 histograms in ṽ×γ — no independence
approximation), NO σ_m anywhere, e-mixture w_circ·U(0,0.2)+(1−w)·power(η), 3K companion
sector, two contaminants with distinct 2D signatures (chance: flat×flat; flyby:
ṽ~U(0.5,3) × γ triangular [0°,30°]), corrected sigv/√2 noise. Grids: α∈[0,2] (5), η
∈{1.05,1.3}, w_circ∈[0,0.3], f_comp∈{0,0.1}, f_c0∈[0,0.1], f_fly∈[0,0.05]. Summary:
[data/stage3m_summary.txt](data/stage3m_summary.txt), cubes data/hier6_cube_*.npy.

| law | α profile [0,0.5,1,1.5,2] | Newton ΔlnL | w_circ | fences |
|---|---|---|---|---|
| simple | [−111.5, −42.1, −3.1, 0, −6.2] | +111.5 | **0.0** | f_comp, f_c0, f_fly at max |
| BE | [−107.3, −50.5, −5.8, 0, −17.3] | +107.3 | **0.0** | same |

1. **α remains interior with a flat top over [1.0, 1.5]** (1.0-vs-1.5 differences of 3–6
   lnL are within single-realization noise) under a likelihood that sees directions —
   Newton loses by +107/+112 even with both contaminants and companions at their fences.
   BE−simple = −4.2: still no shape verdict.
2. **RETRACTION-GRADE CORRECTION of the 3L interpretation: the joint fit VETOES the
   circular component** (w_circ=0.1 costs −174/−179; 0.3 costs ~−890). The physics I got
   backwards in 3L: circular orbits NARROW the ṽ distribution — they cannot be the σ_m
   broadening, and the ṽ channel punishes them harder than the γ 90°-arm rewards them.
   The U-shape's 0°-arm is being absorbed by flybys+chance (both at fence values); whether
   the 90°-arm is genuinely unfit at the best point needs a residual dump (next session).
3. Without σ_m the broadening deficit persists (α's flat top drifts up, all contaminant
   fences engaged). **Identity scoreboard: NOT mass errors (3J), NOT companion wobble at
   capped levels (3K), NOT a circular e-component (3M).** Still standing: larger
   contaminant populations than gridded (extend fences), radial-side e-structure beyond
   the power family (e.g., excess weight at e>0.9), PM-error underestimation (a factor
   f_pm on the noise), intra-bin a/s-distribution.
4. The quotable α stays the 3F–3I number — **α = 0.98 ± 0.23 (simple) / 1.21 ± 0.30 (BE)**
   with the σ_m caveat now sharpened to: "a smooth broadening of undetermined origin,
   shown NOT to be masses, companions, or circular orbits, and which does not remove the
   anomaly under any tested physical identity."

Next-session queue (v6b): extend f_c0/f_fly/f_comp fences; dump best-fit 2D residuals
(is the γ 90°-arm unfit? which cells drive Newton's loss?); add f_pm noise-scale nuisance;
try a radial-excess e-component instead of the circular one.

## Stage 3N: v6b — the broadening IDENTIFIED: a ~20% near-parabolic sub-population; all fences localize; α model-dependence made explicit (2026-07-22 night, [calcs/stage3n_v6b.py](calcs/stage3n_v6b.py))

Four queue items in one joint-2D fit: radial-excess e-component w_rad·[e~0.9–0.995]
(replacing the vetoed circular one), extended contaminant/companion fences, f_pm noise
scale, best-fit residual dump. [data/stage3n_summary.txt](data/stage3n_summary.txt),
cubes data/hier6b_cube_*.npy.

| law | α̂ | α profile [0,.5,1,1.5,2] | Newton | w_rad | best lnL |
|---|---|---|---|---|---|
| simple | 1.5 | [−92.9,−29.2,−6.5,0,−4.3] | +92.9 | **0.20 (interior)** | −56715.7 |
| BE | 1.5 | [−101.6,−46.8,−26.3,0,−24.3] | +101.6 | **0.20 (interior)** | **−56707.0** |

1. **σ_m's identity found (primary component): w_rad = 0.20 is a clean interior optimum
   for BOTH laws** (+129/+141 lnL over none; 0.3 costs −51/−58). A ~20% near-parabolic
   (e>0.9) sub-population at wide separations, ON TOP of the e^1.3 family — beyond
   Hwang's published 22%-of-e>0.9 (which our power family already contains). This is the
   3L U-shape's low-γ arm made quantitative, and it's what the σ_m smear was mimicking.
   The scoreboard closes: not masses (3J), not companions (3K), not circular (3M),
   **radial excess (3N) — plus small localized contaminants**.
2. **Every v6 fence localizes with room to spare**: f_comp=0.1 (of {0,.1,.2}), f_c0=0.1
   (of {0,.1,.2}), f_fly=0.05/0.10 (of {0,.05,.1,.15}). f_pm nearly flat (+5–7 for 1.5×):
   noise scale is a minor player.
3. **α localizes interior at 1.5 for both laws** — but note the model-dependence honestly:
   α̂ = 0.94–1.05 under the σ_m model (3F), 1.25–1.75 under v5, 1.0–1.5 flat under v6,
   1.5 under v6b. **Across broadening identities the anomaly's preferred strength spans
   α ≈ 1.0–1.5; the invariants are: interior localization every time, α > 0 at enormous
   significance, Newton ΔlnL ≈ +55…+112 in every model class.** The α̂-vs-nuisance-model
   covariance is itself a reportable systematic (realization/bootstrap budget for the
   v6b model still owed).
4. **First time ever: BE out-fits simple (−56707.0 vs −56715.7, ΔlnL 8.7)** — near
   realization noise, but the sign is new; track it.
5. Residual map (the v7 agenda): (a) the γ≈82° column stays underfit (+2.7/+2.8 marginal
   z) — a genuinely perpendicular-excess population the ṽ channel won't let us buy with
   circular orbits at current weights; (b) the model OVERpredicts the extreme tail
   (ṽ≈5.2, z≈−6): almost certainly the catalog's own R_chance/boundedness selection
   culling high-Δμ pairs — **a data-side selection we don't emulate, sitting exactly
   where boost signatures live; modeling it could move α̂ either way and is the single
   most important remaining systematic (v7)**. (c) Newton's failure cells are at
   ṽ≈1.25–1.66 with ORBITAL-looking angles (γ 38°–82°), not flyby-like ones — the
   anomaly's directions are bound-orbit directions.

## Stage 3O: v7 — catalog selection closed; α invariant; the BE lead persists (2026-07-22 night, [calcs/stage3o_v7fit.py](calcs/stage3o_v7fit.py))

First, the selection was MEASURED before modeled: the data's velocity envelope sits at
v·√s ≤ 2.2 km/s·kAU^½ — the PHYSICAL escape ceiling of the heaviest pairs (bound orbits
obey ṽ ≤ √2 identically; a virial identity worth writing down), well inside El-Badry's
nominal Δμ_orbit cut (5 M⊙ escape + 2σ ≈ 3.0; every data pair below 80% of it). So the
v6b ṽ≈5 overprediction was companion-wobble pairs the real catalog would have CULLED.
v7 = v6b + that acceptance applied to model pairs and both contamination templates.
[data/stage3o_summary.txt](data/stage3o_summary.txt), cubes data/hier7_cube_*.npy.

| law | α profile [0,.5,1,1.5,2] | Newton | w_rad | best lnL |
|---|---|---|---|---|
| simple | [−97.7,−29.8,−5.0,0,−4.9] | +97.7 | 0.20 | −56339.7 |
| BE | [−105.1,−48.3,−22.3,0,−19.8] | +105.1 | 0.20 | **−56332.2** |

1. **The selection systematic is closed and α is INVARIANT under it.** Total lnL improves
   +376/+375 (the phantom tail mass removed); the ṽ≈5 residual cells vanish; the α
   profiles are statistically unchanged from v6b. The feared bias (catalog cut punishing
   high-boost models) is empirically ~nil at this grid resolution.
2. **w_rad = 0.20 strengthens** (−196/−207 without it): the near-parabolic sub-population
   is more demanded, not less, once the tail is handled honestly.
3. **BE out-fits simple for the second consecutive model class (7.5 lnL after 8.7 in
   v6b).** Same-realization, so not independent evidence — but the sign is stable under
   two structural model changes. If it survives the v7 error budget over seeds, it is the
   first shape-level discrimination between ν-families in the program.
4. Newton's failure cells remain at ṽ≈1.25–1.66 with orbital angles — the anomaly is
   untouched by selection modeling. Sixth→seventh model class, same verdict:
   **Newton ΔlnL ≈ +93…+112, α interior every time.**
5. Last unexplained structure: the γ≈82° perpendicular excess (cells at ṽ≈0.07 AND
   ṽ≈1.66, wide bins, z≈+5) — present in every model generation. Everything else in the
   2D space now fits.

**v7 is the program's most complete model** — physical e-mixture, physical companions,
two contaminants, catalog selection, noise scale, all localized interior — and the
anomaly stands. Remaining before paper assembly: the v7-model error budget (seeds ×
bootstrap, the 3H/3I machinery re-run on this model — the LAST gate) and the γ≈82° chase.

## Stage 3P: v7 realization budget — BE's α̂ is startlingly stable; the BE-vs-simple lead is NOT significant (2026-07-22 night, [calcs/stage3p_v7budget.py](calcs/stage3p_v7budget.py))

Six seeds on the v7 model (radial mixture + companions + contaminants + catalog
acceptance; f_pm fixed 1.5). [data/stage3p_summary.txt](data/stage3p_summary.txt).

| seed | simple α̂ | BE α̂ | Newton s/BE | BE−simple |
|---|---|---|---|---|
| 31 | 1.50 | 1.52 | +98.4/+105.8 | +7.4 |
| 101 | 1.44 | 1.45 | +114.2/+111.8 | −2.4 |
| 202 | 1.30 | 1.46 | +93.4/+104.0 | +10.6 |
| 303 | 1.19 | 1.50 | +94.5/+98.1 | +3.6 |
| 404 | 2.00(edge) | 1.43 | +109.0/+106.1 | −2.9 |
| 505 | 2.00(edge) | 1.41 | +93.8/+97.4 | +3.6 |

**Aggregates: simple α̂ = 1.57 ± 0.35 (interior 4/6; the two edge hits are flat tops with
α=1.0 only 7–9 lnL down); BE α̂ = 1.462 ± 0.042 (interior 6/6); Newton = +100.5 ± 8.9 /
+103.9 ± 5.4; w_rad = 0.20 in 12/12 fits.**

1. **The BE-vs-simple lead is NOT significant: mean +3.3 ± 2.2 (SE), sign flips 2/6.**
   The v6b/v7 same-seed leads were realization noise. No ν-family discrimination —
   honest verdict logged before anyone got attached.
2. **Under the fully-physical v7 model, the BE law yields an extraordinarily
   well-conditioned α measurement: α̂ = 1.46 ± 0.04 realization scatter** (vs the σ_m
   model's ±0.15). The simple law's profile is flatter and wanders to fences on 2/6
   seeds. (Realization scatter only; the data-bootstrap half for v7 — 3I machinery on
   this model — still owed before quoting a total ±.)
3. Newton rejection: +94…+114 across every seed and law. w_rad = 0.20 is as stable as
   σ_m = 0.30 was — same underlying structure, now physically named.
4. **Model-dependence, final statement for the paper: α̂ = 0.98±0.23 / 1.21±0.30 under
   the phenomenological σ_m model; α̂ ≈ 1.5 / 1.46±0.04(real.) under the physical v7
   model. Across every model class: α interior, α > 0 overwhelming, Newton ΔlnL
   +55…+114.** The α>1 preference of the physical model (both laws) is now the
   interesting open question: either the boost really is stronger than the galactic
   calibration (EFE mismatch? a₀ environment dependence?) or a residual model ingredient
   still leans on α. For the BE reading specifically: α=1.46 with the BE-EFE being ~4% weaker
   than simple-ν partially rationalizes BE α̂ > 1, but not to 1.46 — flag, don't spin.

## Stage 3Q: v7 bootstrap — the error budget completes; α = 1.54 ± 0.13 (BE), and the α>1 tension becomes the headline question (2026-07-23 morning, [calcs/stage3q_v7boot.py](calcs/stage3q_v7boot.py))

1000 bootstrap replicates against the stored seed-31 v7 grid ([data/stage3q_summary.txt](data/stage3q_summary.txt)):

- simple: α̂ = 1.52 ± 0.30 (boot) ⊕ 0.35 (real.) → **α = 1.52 ± 0.46**; interior 792/1000;
  Newton +101.0 ± 19.2 (min +38.5).
- BE: α̂ = 1.54 ± 0.12 (boot) ⊕ 0.042 (real.) → **α = 1.54 ± 0.13**; interior 953/1000;
  Newton +106.3 ± 14.6 (min +59.3).

**Newton loses in all 1000 replicates for both laws.** TODO #2d is closed; the fit
program (v1→v7, Stages 3B→3Q) is complete.

The interpretive centerpiece the paper must own: under the most physical model the BE-law
α is measured PRECISELY, and it is ~4σ above the parameter-free α=1 — while the
phenomenological σ_m model gave 1.21 ± 0.30 (consistent with 1). Candidate explanations,
in order of checkability: (1) **g_ext = 1.9a₀ has never been varied** — the EFE table's
external-field strength was fixed at Stage 2G; a weaker true local field ⇒ stronger boost
⇒ α̂ falls toward 1. A g_ext scan (1.4–2.4a₀) is the single highest-leverage next
computation. (2) Residual model ingredient leaning on α (the γ≈82° excess is still
unfit). (3) Real physics: the binary-regime boost genuinely exceeds the
galactic-calibrated law. Do not pick a favorite in print; scan g_ext first.

## Stage 3S: the g_ext scan — the α>1 tension dissolves for the simple law, persists for BE (2026-07-23, [calcs/stage3r_gext_tables.py](calcs/stage3r_gext_tables.py), [calcs/stage3s_gext_fit.py](calcs/stage3s_gext_fit.py))

EFE tables re-solved (validated stage-2G machinery) at g_ext/a₀ ∈ {1.4, 1.6, 2.2, 2.4};
v7 fit (seed 31) at each; Newton is table-independent so ΔlnL(Newton) is comparable
across the scan. [data/stage3s_summary.txt](data/stage3s_summary.txt).

| g_ext/a₀ | simple α̂ | simple ΔlnL | α=1 penalty | BE α̂ | BE ΔlnL |
|---|---|---|---|---|---|
| 1.4 | **1.02 (sharp)** | **+109.5** | 0 | 2.0 (edge) | +96.5 |
| 1.6 | **1.03 (sharp)** | +108.4 | 0 | 1.55 | +99.6 |
| 1.9 | 1.50 | +98.4 | −5.0 | 1.52 | +105.8 |
| 2.2 | 1.31 | +109.3 | −6.1 | 2.0 (edge) | +101.9 |
| 2.4 | 1.41 | +108.3 | −9.7 | 2.0 (edge) | +108.1 |

1. **For the simple-ν law, the α>1 tension is NOT demanded by the data.** At
   g_ext = 1.4–1.6a₀ the fit prefers α = 1.0 SHARPLY (α=1.5 penalized 25–27 lnL — far
   above the ~5-lnL jitter) and achieves the best overall likelihoods of the whole scan.
   (α≈1, g_ext_eff≈1.5a₀) is the global sweet spot: **the parameter-free galactic
   strength is recovered if the effective EFE suppression is ~25% weaker than the
   C&M-convention 1.9a₀.** Elsewhere the (α, g_ext) surface is degenerate at ±1–2 lnL.
2. **The BE law is NOT rescued by g_ext anywhere in the scan**: its α̂ stays 1.5–2.0
   (edges at three of five points). The BE-EFE (Wien-suppressed) is intrinsically weaker,
   and no scanned external field lets it reach the data at α=1. Under the v7 model, the
   parameter-free BE reading is in genuine tension — flag honestly; this is the first
   quantitative separation between the two ν-families the program has produced (the lnL
   totals still don't discriminate; the α-structure does).
3. **The decisive next fact is not a fit — it is the true Newtonian external field of
   the solar neighborhood in our convention.** If MW baryonic models give
   g_N,ext ≈ 1.4–1.6a₀, the simple-law story closes clean at α=1; if 1.9a₀ is robust,
   the tension returns. Calculable/literature number — queued.
4. Caveats: single seed (the α̂ lobe-hopping between 1.0 and 1.5 at g≥1.9 is decided by
   ~5 lnL ≈ jitter; the SHARP α=1 preferences at 1.4–1.6 are not); scan grid is coarse;
   w_rad=0.20 selected at every point (7th consecutive model variation).

## Stage 3T: the external-field convention resolved — α = 1 for BOTH laws at the physical field; the 2G convention bug found and owned (2026-07-23)

Three converging pieces ([data/stage3s_summary.txt](data/stage3s_summary.txt) g=1p2 rows;
scout + Chae 2023 abstract):

1. **Our own arithmetic**: inverting the RAR at the solar circle (v_c = 233±4 km/s,
   R₀ = 8.178 kpc) gives the QUMOND-relevant NEWTONIAN external field
   **g_N,ext = 1.15 ± 0.05 a₀** (simple-ν: 1.10–1.21; BE: 1.13–1.25); direct baryonic
   estimates (v_bary 165–185 km/s) give 0.9–1.13a₀. The 1.9 value cannot be a baryonic
   Newtonian field (it would need v_bary ≈ 240 km/s > total).
2. **The fit at the physical field (g_N,ext = 1.2a₀, tables re-solved)**:

   | law | α̂ | grid best | Newton ΔlnL | α=1.5 penalty | α=2 penalty |
   |---|---|---|---|---|---|
   | simple | 1.17 | **1.0 (interior)** | +108.7 | −4.2 | −76.6 |
   | BE | 1.11 | **1.0 (interior)** | +98.8 | −10.5 | −31.6 |

   **Both ν-families localize at the parameter-free α = 1.** The Stage-3S "BE not
   rescued" verdict was an artifact of scanning 1.4–2.4 — the physical value lies BELOW
   the old floor. Stage-3S conclusions superseded accordingly.
3. **Where the 1.9 came from — the bug is OURS, not the field's**: Chae 2023 is an AQUAL
   analysis (abstract explicitly), where the self-consistent external input is the TOTAL
   field ≈ 1.8–1.9a₀; Banik uses the kinematic total in an AQUAL-family library —
   internally consistent too. Stage 2G imported the AQUAL-appropriate 1.9 into our
   QUMOND solver, whose correct input is the NEWTONIAN 1.2 — a cross-formulation
   convention error (the a₀/IC-convention trap family; caught in our own kitchen, again).
   Sanity check: QUMOND boost @g_N=1.2 (1.40 at y=1) ≈ what AQUAL@total-1.8 should give,
   consistent with formulation near-equivalence. (Scout's claim that Chae labels 1.9
   "purely baryonic" remains to be verified against the full text — Haiku citation risk;
   non-blocking for our pipeline, which is now self-consistent.)

**REVISED HEADLINE (supersedes 3Q's α=1.54±0.13, which was conditional on the wrong
convention): at the physical external field, the wide-binary data prefer the
parameter-free galactic-calibrated boost, α ≈ 1.1 ± ~0.15 (error budget to be re-run at
g=1.2), for BOTH interpolating families, with Newton rejected by ΔlnL ≈ +99–109.** The
w_rad = 0.20 radial population persists (8th consecutive model variation). Remaining
before the paper: re-run 3P/3Q budgets at g=1.2 (mechanical), verify Chae's exact text,
γ≈82° disclosure. Credence updates: anomaly-is-real unchanged ~70–75%; one-law-spanning-
scales UP to ~65% (the α=1 consilience at the physical field is exactly its prediction);
BE microphysics recovers to ~15–20% (the 3S "BE tension" was the convention artifact).

## Stage 4A/4B: the NLO kill test — the ½-branch survives; the 0-branch dies (2026-07-23, [calcs/stage4a_nlo_test.py](calcs/stage4a_nlo_test.py), [calcs/stage4b_branchcomp.py](calcs/stage4b_branchcomp.py))

The BE identity's first kill test (TODO #9): deep-MOND expansion g = √(g_N a₀)(1 + c₁x +
c₂x² + …), x=√(g_N/a₀); BE predicts c₁ = ½ exactly (shared with simple-ν; the standard-μ
branch predicts c₁ = 0; they differ at c₂: BE 1/12, simple 1/8, standard 1/4).

**4A (truncated-expansion estimator): INCONCLUSIVE by design flaw then by power.** First
attempt was ill-conditioned (a₀–c₁–c₂ collinearity; fitted a₀ 3e-11 nonsense — caught
immediately). The two-step asymptotic-matching v2 is honest but power-limited: c₁ to
±0.4–0.6 with window see-saw (x vs x² are ~99% collinear over one window) and f_ML
edge-hopping. SPARC cannot measure c₁ through a truncated series. [data/stage4a_nlo.txt](data/stage4a_nlo.txt).

**4B (truncation-free full-function comparison): DECISIVE.** Fit complete ν-functions
with only (a₀, f_ML) free; galaxy-bootstrap sign-stability:

| window | ΔlnL-proxy Δχ²(BE−standard) | BE better in | Δχ²(BE−simple) |
|---|---|---|---|
| y<0.5 | −7439 | 198/200 | +374 |
| y<1.0 | −13899 | **200/200** | +797 |
| all | −20059 | 199/200 | +129 |

- **The c₁=½ branch defeats the c₁=0 branch at 198–200/200 bootstrap stability** (raw Δχ²
  magnitudes inflated by underestimated per-point errors — quote the bootstrap fraction,
  not σ). Standard-μ also drags f_ML to 1.7–2.0 (unphysical) trying to fake the shape.
  **The BE identity's parameter-free NLO prediction PASSES its kill test; the rival
  branch is excluded on the low-acceleration side** (localizing what the Wien-tail/
  high-y arguments said before to the deep-MOND expansion structure itself).
- Within-branch (BE vs simple, NNLO 1/12 vs 1/8): small sign-consistent lean toward
  simple (+129…+797) — NOT resolved; stays open (consistent with Stage-1's
  p = 0.443 +0.063/−0.050 sitting 1.1σ below the BE p=½).
- Bonus consilience: full-sample fit gives a₀ = 1.206e-10, f_ML = 1.10 — canonical
  values recovered with nothing pinned.

Kill-test scoreboard for the BE identity: **NLO=½: PASSED (branch-level).**
a₀ ∝ H(z): awaiting high-z data (TODO #14). NNLO (1/12 vs 1/8): open, needs better
low-y data or a dedicated estimator. Credence (BE microphysics specifically): 15–20% →
**~20–25%** (it faced a falsifier and the alternative branch died instead; capped
because the test doesn't separate BE from simple within the branch).

## Stage 3U/3V: final error budgets at the physical field — THE MEASUREMENT PROGRAM CLOSES (2026-07-23, [calcs/stage3p_v7budget.py](calcs/stage3p_v7budget.py) g=1p2, [calcs/stage3q_v7boot.py](calcs/stage3q_v7boot.py))

Six seeds + 1000-replicate bootstrap on the v7 model with the g_N,ext = 1.2a₀ tables.
[data/stage3u_summary.txt](data/stage3u_summary.txt), [data/stage3v_boot.txt](data/stage3v_boot.txt).

Seeds (α̂ simple / BE): 1.17/1.11, 1.28/1.17, 1.23/1.09, 1.28/1.10, 1.29/1.07, 1.26/1.06
→ realization scatter 0.045 / 0.039 (vs 0.35/0.04 at the wrong field — the physical model
is beautifully conditioned for BOTH laws), interior 12/12, w_rad = 0.20 in 12/12.

**FINAL HEADLINE NUMBERS:**
| law | α (boot ⊕ realization) | interior | Newton ΔlnL (min over 1000) |
|---|---|---|---|
| simple | **1.18 ± 0.11** | 1000/1000 | +110.3 ± 18.3 (+53.0) |
| BE | **1.13 ± 0.13** | 998/1000 | +100.0 ± 14.7 (+57.8) |

Both consistent with the parameter-free α = 1 (1.6σ / 1.0σ); the mild common upward lean
(~10–18%) is partially covered by the g_N,ext uncertainty (±0.05a₀ ⇒ ~±0.04 on α̂ from
the scan slope) and is disclosed, not interpreted. **Newton is excluded in all 2000
bootstrap contests across both laws.**

Also logged: at the physical field the simple law out-fits BE on every seed
(mean −12.3 ± 2.2 SE), echoing the galaxy-side 4B lean — two independent datasets now
mildly prefer simple-ν over BE within the ½-branch. BE-specific credence trimmed to
~15–20% (branch passed, within-branch leaning away); does not affect the α result.

**The measurement program (Stages 1 → 3V, 22 stages) is COMPLETE. Every gate is closed:**
screening index p = 0.443+0.063/−0.050; boost 1.086 (1.064–1.110); α = 1.18 ± 0.11 /
1.13 ± 0.13 at g_N,ext = 1.15 ± 0.05a₀; Newton rejected universally; the w_rad = 0.20
near-parabolic population and the U-shaped γ measurement as standalone discoveries; the
NLO ½-branch confirmed on SPARC; six retraction-grade corrections caught and logged.
Remaining work is WRITING (paper assembly, TODO #10) plus the disclosed open items
(γ≈82° excess, Chae-text verification, NNLO, a₀∝H(z), Banik reproduction #6, deep
priority check #7, two-field theory #8).

## Stage 4C: PRIORITY RETRACTION — the BE identity IS published: Cadoni & Tuveri 2019 (2026-07-23)

Deep priority sweep (25+ searches) surfaced arXiv:1904.11835 (Cadoni & Tuveri, "Galactic
dynamics and long-range quantum gravity", May 2019) with the verdict "close but does not
state the identity." **Reading the primary source shows the verdict was wrong — the
scouts misread it, twice.** Direct quotes from the paper:

- Their Eq. (23): "The McGaugh form for F(x) leads to the additional acceleration term
  a_DF = a_B/(e^√(a_B/a₀) − 1)" — **this IS a_B · n_BE(√(a_B/a₀)), i.e. exactly our
  identity ν = 1 + n_BE(x).**
- Their Eq. (22), DERIVED from physics (soft spin-2 bosonic excitations of the dark
  energy in thermal equilibrium with the dS horizon at T_dS = ℏ/2πL, occupation N(ε)):
  a_DF = 2π·a_B/(e^{2π√a_B/H} − 1), matching Eq. (23) with **a₀ = H/2π derived** via a
  generalized thermal equivalence principle (their §4).
- Even the limits-reading is theirs in substance: "the Newtonian regime corresponds to
  hard DE bosonic excitations… N goes to zero exponentially [Wien]; the MOND regime
  corresponds to a huge number N≫1 of extremely soft excitations [Rayleigh–Jeans]."

**RETRACTED: every "apparently unpublished" claim about the identity (Stages from
2026-07-22 onward, both memory files, PAPER-DRAFT §2 framing). The identity, the
horizon-bath reading, and a₀ = cH₀/2π are Cadoni & Tuveri's (and companion papers —
citation chain to be mapped: their refs [22,23,39,73]).** Seventh correction of the
project; caught by our own primary-source read after three scout sweeps failed to.

What remains OURS: (1) the NLO=½ expansion coefficient as an explicit falsifiable
prediction AND its SPARC branch test (4B: ½-branch 198–200/200 over standard-μ) — C&T
do not extract or test expansion coefficients; (2) the entire wide-binary program
(α = 1.18±0.11/1.13±0.13, independent of them); (3) the NNLO 1/12-vs-1/8 discriminator;
(4) the screening-index p measurement. **The paper's §2 becomes: "Cadoni & Tuveri derived
the RAR as a Bose–Einstein occupation effect; we subject that framework to two new
falsification tests and an independent low-acceleration measurement." Honestly, this
STRENGTHENS the paper (the theory has refereed pedigree; we are its testers, not its
parents) while ending our novelty claim on the identity itself.**

Credence bookkeeping: horizon/BE microphysics as the true story: ~15–20% → **~25%**
(independent professional derivation existed all along — the idea has more legs than we
knew); OUR contribution being the identity: 0% (it never was). The project's empirical
results are untouched.

Addendum — citation-tree walk (scout-grade, 24 searches): C&T's companion papers
(corpuscular-gravity lineage: Cadoni/Casadio/Giusti/Tuveri 2017–2020, Symmetry 2020,
PRD 102) and ~65 citing works stay on rotation curves; **no expansion-coefficient
extraction, no coefficient tests, no wide-binary application anywhere in the tree.**
Claims map for the paper (provisional until a final INSPIRE pass at write-time):
identity + bath + a₀=H/2π = C&T 2019 (cite; one-line independent-arrival note);
NLO=½ extraction + SPARC branch test (4B), NNLO discriminator, wide-binary program
(α=1.18±0.11/1.13±0.13), radial population, U-shaped γ, realization systematic = ours.

## Stage 4D: the Bernoulli ladder and the zero-point half — the expansion program becomes a structure (2026-07-23)

Two exact statements, both verified symbolically/numerically (sympy series; coth identity
to 5e-13; sympy's simplify fails on it — algebra: ½+½coth(x/2) = e^x/(e^x−1) = ν ✓):

1. **ν_RAR(x) = ½ + ½·coth(x/2)** — the RAR interpolating function IS the Planck
   quantum-oscillator mean-energy law E/ħω = ½coth(ħω/2kT), zero-point term included.
   Readings: deep MOND = classical equipartition; **Newton = the FROZEN oscillator:
   ν(∞) = ½ + ½ = 1 — Newtonian gravity as the pure zero-point response, the thermal
   part switched off. The NLO constant ½ that Stage 4B tested and confirmed (198–200/200
   over the 0-branch) is the ZERO-POINT OCCUPATION of the bath modes.** "We empirically
   detected the vacuum half-quantum of the gravitational bath" is the sharp form of what
   4B measured.
2. **The low-acceleration expansion coefficients are exactly the Bernoulli numbers**:
   ν = Σ_{n≥0} B_n⁺ x^{n−1}/n! = 1/x + 1/2 + x/12 − x³/720 + x⁵/30240 − … (all even
   powers of x vanish with the odd Bernoullis). The ladder vs simple-ν (1/x + 1/2 + x/8
   − x³/128 + x⁵/1024): rung 0 (BTFR) shared; rung 1 (½) shared — TESTED ✓; rung 2:
   **1/12 vs 1/8 (ratio 2:3)**; rung 3: **−1/720 vs −1/128 (5.6×)**; rung 4: 1/30240 vs
   1/1024 (~30×). **Discrimination grows exponentially with depth** — the within-branch
   BE-vs-simple question is decidable in principle by a deep-enough RAR.
3. Program: SPARC lacks depth for rung 2 (4A: σ(c₁)≈0.5). The dataset with the reach is
   the **weak-lensing RAR (Brouwer et al. 2021, KiDS: ~2 dex deeper in acceleration)** —
   rung-2/3 measurement there = the within-branch decider = new TODO #16. (Current lean:
   both SPARC full-function fits and the binary fits mildly favor simple's 1/8 — the
   ladder test could KILL the BE reading; that is what makes it worth running.)
4. Priority: dedicated sweep (15 searches) — **both the coth/zero-point statement and
   the Bernoulli-ladder statement: NOT FOUND anywhere** (nearest: Pazy & Argaman, no coth,
   no expansion; C&T themselves never expand — verified by our own read). Honest
   qualifier for the paper: since C&T published ν = 1+n_BE, the coth form is an
   ELEMENTARY COROLLARY of their identity that nobody stated — our claim is the
   statement, the zero-point READING (Newton = frozen vacuum response; the tested ½ IS
   the zero-point occupation), and the LADDER as a falsification program with its first
   two rungs run. Frame as sharpening C&T, not as independent theory. (Negative scout
   claims stay provisional until the write-time INSPIRE pass — thrice-learned lesson.)

## Stage 4E: the rung-2 lensing test — an honest null with structure; correction #8 (2026-07-23, [calcs/stage4e_lensing_rar.py](calcs/stage4e_lensing_rar.py), [calcs/stage4e_diag.py](calcs/stage4e_diag.py))

TODO #16 executed: can the weak-lensing RAR (2 dex deeper than SPARC) decide rung 2 of
the Bernoulli ladder (BE 1/12 vs simple 1/8)? Answer: **no — and the attempt taught us
two things about our own earlier claims.** [data/stage4e_lensing.txt](data/stage4e_lensing.txt),
[data/stage4e_diag.txt](data/stage4e_diag.txt).

**Data acquired** (data/lensing_rar/, URLs + provenance in script header):
- PRIMARY: Mistele–McGaugh–Lelli–Schombert–Li 2024 (JCAP 04, 020) Table 1 — 15 stacked
  exact-deprojection points, isolated KiDS-1000 lenses, log g_bar ∈ [−14.86, −11.41],
  extracted verbatim from the arXiv LaTeX source (primary-source rule). Their global
  0.2-dex stellar-mass systematic modeled as a log-g_bar offset nuisance δ ~ N(0, 0.2).
- CROSS-CHECK: Brouwer+21 KiDS-1000 release (ESD + full covariance; their Eq.-7
  SIS-approx conversion — the one Mistele+24 showed distorts the shape at both ends).

**Method upgrade over 4B (and the caveat it forces).** 4B compared families with raw χ²
(χ²/dof ≈ 57 — no intrinsic scatter; the most precise points carry everything). 4E uses
−2lnL with a per-family profiled intrinsic scatter s_int (≈ 0.12 dex, matching the known
RAR scatter) + the lensing block + the δ prior; free (a₀, f_ML, s_int, δ) per family.
Gates: G1 a₀ sane ✓; G2 exact regression to stored 4B numbers (the first run's "FAIL"
was 4B's 'all y' actually being y<30 — config, not code) ✓; G3 Newton catastrophic on
lensing alone (+2777 Mistele / +1659 Brouwer GLS over 15 points) ✓; G4/G5 below.

**Results:**
1. **Rung 2 is UNRESOLVED — an honest null with the power quantified.** The direct c₂
   estimator, even with the lensing-anchored a₀: c₂ = +0.14 ± 0.23 (stat) ± 0.37 (mass
   syst) in the best window — **resolving power for 1/12-vs-1/8: 0.09–0.10σ.** The
   0.2-dex lensing stellar-mass calibration is the wall (rung 2 needs ~0.01–0.02 dex).
   The c₁-free variant reproduces 4A's collinearity pathology (c₁ ≈ −0.5 ± nonsense) —
   consistent with 4A's power-null, now WITH the anchor: the truncated estimator is dead
   as an instrument at current calibration, full stop.
2. **Full-function within-branch: Δ(−2lnL) BE − simple = −18.7 (BE better), bootstrap
   −19.5 ± 15.4, BE better in 183/200.** Sign stable across all 14 variants (SPARC
   window × lensing depth cut × EFE-nuisance grid e_N ∈ {0.01–0.05} × face-value
   masses), magnitude −3…−23. BUT the diagnostics kill any verdict reading: the entire
   lead is carried by **3 of 153 galaxies** (UGC03580 −7.1, NGC4217 −5.9, UGC02916 −5.4;
   refit without them: +0.8), and the deep regime where c₂ actually lives (y<0.1,
   960 points) contributes only −7.8. **No within-branch discrimination — again.**
3. **CORRECTION #8a (to 4B/3V's "two datasets lean simple"):** the SPARC-side "slight
   sign-consistent simple lean" (+129…+797 raw Δχ²) does NOT survive scatter
   marginalization — it flips to the small, 3-galaxy-fragile BE lean above. The SPARC
   within-branch comparison is likelihood-model-dependent noise, not a lean. RETRACTED
   to: "the wide binaries lean simple (3V: −12.3 ± 2.2 SE, stands); SPARC is agnostic."
4. **CORRECTION #8b (the 4B branch kill gets its likelihood caveat):** under the honest
   likelihood the ½-branch vs 0-branch verdict deflates from 198–200/200 to
   **Δ(−2lnL) = −56, BE better in 166/200** joint — and the deep window (y<0.5) alone
   flips to +8…+10 FOR standard-μ (it buys the shape back with f_ML=1.69 + scatter
   freedom once precise points lose their veto). The ½-branch preference is sign-robust
   in every joint variant and both likelihood treatments — **the conclusion stands as a
   strong lean, not a kill** — and the paper must quote both treatments (raw-χ²
   198–200/200 AND scatter-marginalized 166/200) or it is overclaiming.
5. a₀ (joint, BE) = (1.000 ± 0.094)e-10 under the honest likelihood vs 1.206e-10 raw —
   a ~20% likelihood-model sensitivity (f_ML-correlated), bracketing Stage-1's
   (1.03 ± 0.13)e-10 and C&T's H₀/2π ≈ 1.08e-10. Not a claim; a disclosed sensitivity.
6. Lensing alone is family-agnostic (Mistele: −60.6/−60.6/−61.4 across BE/simple/
   standard) — it anchors, it does not discriminate. The Brouwer-release GLS cross-check
   mildly prefers standard (38.3 vs 46.6, a₀ ≈ 1.9e-10) — attributed (not proven) to
   their SIS-approx conversion + face-value masses; noted as a caution, outweighed by
   the primary deprojection dataset.
7. **In-session artifact caught (logged as method discipline):** the first diagnostic's
   per-decade table binned points by each family's OWN fitted y — bin migration between
   columns manufactured a spectacular fake story ("simple wins the deep side +136, BE
   wins the Wien tail −113"). Common-membership binning (v2) shows every regime delta is
   single-digit-to-±11 and columns reconcile exactly with the fit totals. Decomposition
   tables must bin on common membership; the pretty version was the wrong version, again.

**Kill-test scoreboard after 4E:** NLO = ½: PASSED as a robust-sign strong lean (both
likelihoods, every window, plus lensing consistency) — no longer quoted as 198–200/200
alone. NNLO (1/12 vs 1/8): **open and NOT reachable with current data** — needs either
lensing mass cross-calibration at the 0.02-dex level (a survey problem, not ours) or the
Gaia DR4 distribution-level binary measurement. a₀ ∝ H(z): future. Credences: BE
microphysics stays ~20–25% (the feared SPARC simple-lean evaporated — relaxing pressure —
but no positive discrimination appeared either; net wash). All wide-binary numbers
untouched (this stage is galaxy-side).

## Stage 4F: the bath matrix — simple-ν IS the classical bath; the ¼-branch tested and disfavored (2026-07-23, [calcs/stage4f_bathmatrix.py](calcs/stage4f_bathmatrix.py))

Theory-driven stage: taking the Planck-oscillator reading seriously generates a 2×2
matrix of bath laws (occupation statistics × mode-frequency prescription), and it turns
out we had been testing three of its four cells without knowing:

| | source-driven ω ~ √(g_N a₀) | self-consistent ω ~ g_tot |
|---|---|---|
| **quantum (Planck n)** | BE/RAR: ν = 1+n_BE(√y) — c₁=½, c₂=1/12 | **boot: ν = 1+n_BE(νy) — c₁=¼, c₂=7/96 (NEW)** |
| **classical (n = kT/ℏω)** | ν = 1+1/√y — c₁=1 (Cassini-dead a priori) | **ν = ½+√(¼+1/y) = EXACTLY simple-ν** |

**The identity in the lower-right cell is exact and (pending scout) possibly unremarked:
solving classical equipartition self-consistently, ν = 1 + 1/(νy), gives νy² − νy… →
ν = ½ + √(¼ + 1/y) — the Famaey–Binney "simple" function, the empirical workhorse of
the MOND literature, derived in two lines as a classical thermal bath.** Consequence:
the BE-vs-simple stalemate (4B/4E, wide binaries 3V) is secretly the question "is the
bath quantum or classical at the transition?" — and the ½ that 4B confirmed has TWO
readings (zero-point occupation, or self-consistency algebra), which tempers 4D's
"detected the vacuum half-quantum" rhetoric into "confirmed the ½-branch, whose two
members are the two thermal-bath readings."

The untested cell — the self-consistent QUANTUM bath ("boot", solved per point by
seeded Newton iteration; gates: series c₁ = 0.2500 exact / c₂ = 0.0741 vs 7/96 = 0.0729
✓, solver residual 5e-15 ✓, exact regression to 4E fiducials ✓) — was fit under BOTH
likelihood treatments ([data/stage4f_bathmatrix.txt](data/stage4f_bathmatrix.txt)):

- **Raw χ² (4B objective): boot is dead-grade.** +4693…+12313 vs BE/simple across all
  windows; galaxy bootstrap at y<1: boot beats BE in 4/200, simple in 7/200. It does
  beat standard-μ everywhere (−2746…−7874) — ¼ > 0.
- **Scatter-marginalized joint (4E objective): boot is disfavored, not dead.** Fiducial
  Δ(−2lnL): +27.3 vs BE, +8.6 vs simple, −28.5 vs standard; bootstrap +25.8 ± 61.6
  (boot better in 62/200) and +5.6 ± 74.1 (89/200). Spreads are 4× the BE−simple
  spread — partly genuine (boot's residual pattern is more galaxy-lumpy), possibly
  optimizer-jitter-inflated (caveat logged; point estimates match bootstrap means, so
  no systematic bias). Boot also drags f_ML up (1.40 fiducial; 1.58–1.84 raw windows) —
  the same buying-shape-with-M/L tell that marked standard-μ.
- Both treatments AGREE IN SIGN: the ¼-branch loses to both ½-branch members.
- CSD (classical/source-driven, c₁=1): best raw fit in the deep windows (the known
  deep-end upturn likes extra boost) but worst overall (−8169 honest, +228 vs BE) and
  solar-system-dead independently — included only to complete the matrix.

**The c₁ dose-response now has four tested points: c₁ = 0 dead (4B), ¼ disfavored (4F),
½ preferred (4B/4E, both members), 1 screened out (4F + Cassini). The likelihood along
the branch axis peaks at ½.** (Caveat: these functions differ beyond c₁ — high-y
behavior varies from e^−y to 1/√y — so this is a function-family ranking, not a pure
coefficient scan.)

Fantasy scoreboard (the reading stays a reading): within the thermal picture, the
SOURCE-DRIVEN frequency prescription — C&T's original — survives its first structured
alternative; the quantum-vs-classical bath question is exactly the rung-2 discriminator
and stays parked with it (4E: 0.1σ reach). Priority scout (Haiku, 12 searches): **no
exact match found for either result** — not the ν = 1+1/(νy) ⇒ simple-ν identity, not
the BE bootstrap/¼-branch. Nearest art: the equipartition-on-holographic-screen family
(Pazy & Argaman arXiv:1106.4108; Debye-entropic arXiv:1206.1030, 1302.4411 — modified
equipartition on a screen, NOT self-consistent mode frequency) and a 2024 "Timeflow
Gravity" paper (IOP) claiming simple-μ as a thermodynamic equilibrium form — **full-text
checks owed at write time on Timeflow-2024 and Zhao astro-ph/0512425** (whose √(1+4x)
is probably just simple-family algebra, but verify). C&T's own derivation is
source-driven (our 4C primary read), so the bootstrap cell is not theirs. All novelty
claims stay PROVISIONAL until the write-time INSPIRE pass (thrice-learned rule).
Credences: BE microphysics unchanged ~20–25% (its prescription survived a rival; its
statistics remain untested vs classical).

## Stage 4G: priority audit of the U-shape — Hwang, Ting & Zakamska 2022 have the superthermal trend; our w_rad becomes their confirmation (2026-07-23, correction #9)

Primary-source read of arXiv:2111.01789 (MNRAS 512, 3383) — the "verify" item Stage 3L
flagged and never closed (we even quoted "Hwang's ~1.3" secondhand while claiming the
U-shape unpublished — the C&T failure mode, fourth instance; the rule is now: primary
reads at CLAIM time, not write time). Verbatim findings from the PDF:

- **Method**: "v-r angles" = our γ (their convention unfolded 0–180°), on the SAME
  El-Badry EDR3 catalog, s = 10^1.5–10^4.5 AU. Method lineage: Tokovinin 1998/2020,
  Shatsky 2001.
- **Their headline**: e-distribution "close to uniform at ~100 AU… thermal at
  ~10^2.7 AU… superthermal (α > 1) at > 10³ AU". Table 1: α = 1.30±0.05 (1–3.16 kAU),
  1.32 +0.09/−0.08 (3.16–10 kAU), 1.17 +0.14/−0.15 (10–31.6 kAU). Multi-step variant:
  "e < 0.3 suppressed, e > 0.9 enhanced."
- They explicitly attribute the enhanced 0°/180° arms of the v-r distributions to the
  superthermal e (their Fig. 6) — **the radial arm of our U-shape, published, with the
  same explanation.**
- They explicitly propose the gravity application (their §4, crediting Banik & Zhao
  2018, 2021): "For non-Newtonian gravity, wide binaries at ≳ 7000 AU would deviate
  from Keplerian orbits and thus the v-r angle distribution can be an independent test
  on gravity theory."

**RETRACTED/AMENDED (correction #9): "U-shaped γ apparently unpublished" and "U-shaped
γ = standalone discovery" (3L, 3V summary, CLAUDE.md, paper-leads).** The superthermal/
radial-heavy e-distribution at wide separations and its v-angle signature are Hwang,
Ting & Zakamska's. The 53-search scout sweeps missed content sitting in a paper whose
headline number 3L itself quoted.

**What survives as ours:** (1) the joint 2D (ṽ×γ) gravity-law × e-mixture likelihood
EXECUTED — they proposed the test, we ran it, with contaminant fences, at 2–50 kAU
(including 31.6–50 kAU, beyond their range); (2) w_rad quantified INSIDE a gravity fit;
(3) the σ_m-identity chain (what the ṽ-broadening IS); (4) the folded γ≈82°
perpendicular excess (no 90° structure reported by them; still open, TODO #2e).

**THE DIVIDEND — the collision is a validation.** The e>0.9 mass fraction of their
superthermal law is 1 − 0.9^(α+1) = **20.4–21.7%** across their three wide bins —
against our independently-fitted **w_rad = 0.20** (12/12 fits, 8 model variations,
2–50 kAU). Two pipelines with disjoint assumptions — theirs directions-only, Bayesian,
Newton-assumed, ≤31.6 kAU; ours joint 2D, gravity-law free — agree on the
near-parabolic fraction to ≲1.5 percentage points. **The paper's largest nuisance
parameter is now externally validated.** (Also imported: their <3% contamination
estimate at s > 10⁴ AU — consistent with our fenced contaminant fractions.)

Paper effects: abstract and §7 recast from "standalone discovery" to "confirms and
extends Hwang+22 (to 50 kAU, inside a joint gravity fit) with the w_rad ↔ α
cross-validation as the new headline of that section"; cite Banik & Zhao 2018/2021 for
the γ-gravity proposal. Credences: w_rad-is-real-physics UP (hostile-assumptions
replication); the α measurement untouched (γ pinned the nuisances; ṽ carries α).

## Stage 4H: the M/L-marginalized screening index — p = 0.58 ± 0.12 (2026-07-23, [calcs/stage4h_p_ml.py](calcs/stage4h_p_ml.py), TODO #4 closed)

Stage 1 held M/L fixed at (0.5, 0.7). Three-step closure
([data/stage4h_p_ml.txt](data/stage4h_p_ml.txt)):
(a) regression: the Stage-1 fit reproduces EXACTLY (p = 0.443 +0.063/−0.050, seed 42);
(b) global disk-M/L scale free: **p = 0.578 +0.103/−0.097, f_d = 1.22 ± 0.10
(M/L_disk = 0.61 — consistent with the 4B/4E f_ML ≈ 1.1–1.2 and stellar-population
values), a₀ = (1.05 ± 0.10)e−10** — the fit is M/L-degenerate at nearly constant rms
(0.1316 vs 0.1324 dex), so freeing M/L moves p by +0.14 while barely improving the fit;
(c) per-galaxy 0.1-dex M/L jitter: +0.041 shift, 0.047 scatter.
**FINAL: p = 0.578 +0.121/−0.115 (stat ⊕ jitter syst ⊕ shift).** Consequences:
p = ½ (the RAR/BE form) now sits 0.7σ INSIDE the band (the fixed-M/L 0.443 sat 1.1σ
below — the marginalization moves the data TOWARD the ½-branch); Cassini stays
comfortably passed (16th pct 0.462 vs required 0.234). Honest reading: p is a soft
discriminator (M/L-sensitive at the ±0.1 level; quote the marginalized number as
primary); the Cassini screening floor is the robust part of the Stage-1 result.

## Stage 4I: chance-alignment stress test — STABLE under a 20× tightening (2026-07-23, [calcs/stage4i_rchance.py](calcs/stage4i_rchance.py), TODO #5 closed)

R_chance threshold scanned 0.1 → 0.0005 (baseline 0.01), full re-derivation of ṽ per
cut ([data/stage4i_rchance.txt](data/stage4i_rchance.txt)):
- Boost ratio median ṽ(6–30)/median ṽ(0.2–2): 1.091, 1.090, **1.086 (baseline)**,
  1.094, 1.077, 1.088, 1.119 — max drift 0.033; every tightened subsample's 68% CI
  overlaps the baseline CI. VERDICT: STABLE.
- Direction check: the 20–50 kAU median RISES as the cut tightens (0.667 → 0.810 at
  N = 45) — the OPPOSITE of a chance-contamination bias (flat-ṽ chance pairs inflate
  medians; removing them harder would lower it). The (noisy, N=45) drift direction if
  anything strengthens the wide-bin signal.
- Loosening to 0.1 (N_20–50: 214 → 358) moves the ratio by < 0.005 — the baseline cut
  already sits past the dose-response knee.
Chance alignments cannot be driving the boost; consistent with v7's fenced f_chance.

## Stage 4J: the γ≈82° excess resolved in identity — and a perpendicular velocity CEILING at the boosted escape speed (2026-07-23, [calcs/stage4j_gamma82.py](calcs/stage4j_gamma82.py), TODO #2e)

Data-side autopsy of the last unexplained structure (top-γ cells at ṽ≈0.07 and ≈1.66,
wide bins, z≈+5 vs every model). Framing identity (vis-viva at γ_3D = 90°): apocenter
ṽ² = 1−e (ṽ=0.07 ⟺ e≈0.995); pericenter ṽ² = 1+e ≤ 2 — **a bound Newtonian pair can
NEVER exceed ṽ = √2 = 1.414 (any phase, any e, any geometry); the boosted α=1 law
raises the ceiling to √(2·1.36) ≈ 1.65 at wide s.** The observed island at 1.66 sits
exactly there. [data/stage4j_gamma82.txt](data/stage4j_gamma82.txt).

Findings (wide pairs s ≥ 6 kAU):
1. **T1 identity:** both islands are γ-STRUCTURED, not flat — apo island 35% in
   [75,90°] (z=+2.9 vs uniform); peri window bimodal (40% at [0,30°] = asymptote arm,
   24% at [75,90°] = closest-approach arm, [60,75°] DEPLETED z=−2.5) — the geometry of
   extreme-phase orbits, not of comoving noise. T3: slow top-γ pairs imply e median
   0.957 (16–84%: 0.90–0.99), at 4.9× the noise floor — the APO FACE of the w_rad
   population, real velocities. (The S/N>3 cut used in the γ channel removes 92% of
   the apo island — the fits saw it only through the noise-convolved ṽ channel.)
2. **T2 THE CEILING (the new result):** top-γ column (γ≥75, N=256): 11 pairs in
   [1.2,1.414), **11 pairs in the Newton-forbidden band [1.414,1.67) at nearly equal
   density (43 vs 51 per unit), then a CLIFF: 1 pair in [1.67,2.2), 0 beyond.** The
   γ<30 comparison column (the flyby/unbound continuum) shows NO edge: 17, 14, 9, 4, 1
   out to ṽ=6 — and the catalog's own selection admits pairs to ṽ≈3, so the cliff is
   not the catalog cut.
3. **T2b leakage null (the decisive number):** with the measured per-pair σ_ṽ = 0.044,
   a TRUE Newtonian edge at √2 predicts 0.9 pairs in the forbidden band —
   **11 observed, P(≥11) = 3.8e-9.** A true edge at the boosted 1.65 predicts 11.6 in
   the band (11 observed, P=0.62) and 0.5 beyond (1 observed, P=0.91). **The
   perpendicular speed distribution terminates at the α=1 boosted escape ceiling, not
   the Newtonian one — with no parameter tuned to this test** (1.36 is the α=1 table
   boost at wide s from the physical-field fits).
4. **Rival identities, each killed by its own numbers:** field-star flybys at closest
   approach arrive at ṽ ~ v∞/v_c ~ 100, not 1.5, and have no cliff mechanism; the
   fast island's astrometry is PRISTINE (RUWE_max median 1.06 vs 1.28 for the same-ṽ
   low-γ continuum; S/N=30; Rch ~ 8e-4 — physical pairs at 99.9%); 3J's measured
   2.4% mass error cannot produce +15% ṽ shifts; hidden companions inflate by ≤ q-cubed
   amounts on 12% of pairs — order 1–2 pairs, not 11, and would smear past 1.67.
5. **Model-gap diagnosis (why the residual existed):** v7's flyby template puts ALL
   unbound weight at γ∈[0,30°] (asymptote regime) — the 90° closest-approach arm is
   unmodeled by construction; and the bound radial component under-produces the
   extreme-phase cells (e-band capped at 0.995; peri dwell under-sampled). The excess
   was a MODEL-SHAPE gap, and its content is PRO-boost, not a threat.

Honest caveats, stated with the claim: N=11 in the forbidden band (small); the ceiling
location inherits the ±few-% spread of boost(s, g_ext) per pair (consistent with the
1 pair at 1.67–2.2); the islands were FLAGGED by earlier residual maps (the
ceiling-EDGE test itself is new and was not used to find them — partial independence);
per-pair line-of-sight projection can only LOWER ṽ below the 3D value, so projection
cannot fake the band (it thins it). **Verdict: TODO #2e resolved in identity; the
"excess" is the near-parabolic population's pericenter/apocenter faces plus an
unmodeled closest-approach arm — and it contains the program's most direct single
signature of the boost: eleven perpendicular-moving pairs that Newton forbids and the
α=1 law places exactly at its own escape edge.** Queued (#2e-b): GPU model variant
(closest-approach arm + e-ceiling lift) to confirm the cell closes with α invariant;
DR4 multiplies the band's N by ~10. Paper: new §7 subsection "the perpendicular
ceiling" with the caveats verbatim.

## Stage 4K: the solar quadrupole — our own boost, fed to our own solver, breaks Saturn (2026-07-23, [calcs/stage4k_quadrupole.py](calcs/stage4k_quadrupole.py))

The Planet-9/Oort question ("can our numbers predict outer-solar-system dynamics?")
turned into the program's most consequential NEGATIVE result. The EFE solve is
scale-free: the same (GM=1, a₀=1) solution that produced the wide-binary boost tables
IS the Sun-in-the-galactic-field solution, with r_M = 7,032 AU and the anomalous
interior quadrupole δφ = q·r²·P₂(cosθ) along the galactic axis. Extracting the ℓ=2
moment from the validated Stage-2G solver ([data/stage4k_quadrupole.txt](data/stage4k_quadrupole.txt)):

**Gates (all pass):** G1 exact regression of the ℓ-averaged boost against the stored
g1p2 tables (0.00%); G2 Newton control φ₂ ≡ 0; G3 interior plateau flat to 1–2% over
r ∈ [0.02, 0.2] r_M; G4 resolution doubling 0.0%; G5 analytic ℓ=2 integrator test
exact to machine precision; **G6 external cross-validation: Blanchet & Novak 2011
(arXiv:1010.1349, AQUAL at g_ext = 1.9e-10 TOTAL) report Q₂ up to 4.1e-26 s⁻² for
μ₁ = simple — our QUMOND solve at the matched physical config (e_N = 1.2 NEWTONIAN,
the Stage-3T mapping) gives 3.35e-26 raw: ~15% agreement across formulations AND an
independent validation of our AQUAL-total ↔ QUMOND-Newtonian convention resolution.**

**Numbers.** q̂ = −0.0978 (simple) / −0.0988 (BE) at e_N = 1.2, in units a₀/r_M =
1.141e-25 s⁻²; B&N convention Q₂ = 3|q̂|·(a₀/r_M)·α:

| law | Q₂ (α-scaled) | Cassini (Hees+ 2014: (3±3)e-27, 2σ cap 9e-27) |
|---|---|---|
| simple | (3.95 ± 0.39)e-26 s⁻² | **exceeds the cap ~4.4× (~12σ from the measured Q₂)** |
| BE | (3.82 ± 0.46)e-26 s⁻² | **exceeds the cap ~4.2× (~12σ)** |

1. **BE = simple here (1% apart).** The quadrupole is sourced in the TRANSITION region
   (r ~ r_M), where the two RAR-compatible ν's are near-identical — the same fact that
   made rung 2 hard makes the quadrupole family-blind. **The Wien tail does NOT rescue
   BE**: exponential screening kills the local (ℓ=0) anomaly at Saturn (the p > 0.234
   test we pass) but not the long-range ℓ=2 moment of the transition region.
2. **Priority: the tension is Desmond, Hees & Famaey 2024's (arXiv:2401.04796 — "8.7σ
   under fiducial assumptions; requires a sharper transition than the RAR allows";
   verified via abstract).** Ours is an INDEPENDENT, differently-calibrated
   reproduction: their amplitude came from RAR fits (escapable via M/L freedom →
   1.9σ, and bulge removal); **ours comes from the wide binaries — which have no
   bulges and whose mass errors we MEASURED at 2.4% (3J). The two DHF mitigation
   routes do not apply to the binary-calibrated version.** What survives as escapes:
   (a) **modified inertia** (NEED NOT produce this quadrupole — MI's time-nonlocal
   EFE depends on frequency ratios, Milgrom 2011 arXiv:1111.1611; 4L-scout
   correction: evasion is NOT demonstrated in print, and DHF 2024 explicitly frame
   the quadrupole as the MI-vs-MG *diagnostic* — and our trajectory-sensitive data
   (w_rad e-mixture, the 2D ṽ×γ likelihood, the 4J perpendicular ceiling) are
   exactly the discriminating instrument — this becomes the program's top theory
   question, TODO #18); (b) EFE-screened /
   two-field formulations (TODO #8 now has a sharp job description: produce the
   binary boost with Q₂ ≤ 9e-27); (c) the boost is not gravity (the tension is
   honest evidence FOR the systematics reading of our own measurement — logged).
3. **Consequence for the MOND-Planet-9 story** (arXiv:2304.00576 — byline verify at
   write time, Jones-Smith vs Brown & Mathur; critique arXiv:2403.09555): their ETNO
   clustering mechanism RUNS on this same capped term — the P9-alternative inherits
   the Cassini tension (as Vokrouhlický+ noted). Our α-calibrated amplitude
   quantifies it.
4. Systematics on our side, stated: the α-scaling of the anomaly is linear-order
   (±10%-grade); QUMOND-vs-AQUAL formulation difference ~15% (G6); g_ext ± 0.05a₀ →
   ±4%. None approach the 4.3× gap.

**Credence bookkeeping (the honest hit):** anomaly-is-real: unchanged (~70–75% — the
binary data are untouched). One-law-spanning-scales as PHENOMENOLOGY: 65% → ~55%
(coherence hit). **The law realized as modified GRAVITY (AQUAL/QUMOND-type): major
haircut — the same field configuration that boosts the binaries misses Saturn by 4× —
conditional credence ≤ 30%; modified-INERTIA realization correspondingly up (~50%
conditional); BE microphysics (a force-side framework) inherits the hit: ~20–25% →
~15%.** The next decisive computation is ours to run: MI-vs-MG on our own 2D data.

## Stage 4L: MI-vs-MG on our own data — the binaries vote for modified gravity; the paradox is now fully ours (2026-07-23, [calcs/stage4l_mi_runner.py](calcs/stage4l_mi_runner.py), TODO #18 executed)

The post-4K question, run same-day. Design: modified inertia = the boost is a
per-ORBIT functional (Milgrom 2011 gives the principle, no closed form — we BRACKET),
so an MI orbit is exact Kepler dynamics with G_eff = B·G: implemented as the engine's
Newton mode with M_eff = M·B(y_char), ṽ normalized downstream by the TRUE mass —
**the entire v7 nuisance machinery is bit-identical between MI and MG; the lnL
difference is pure theory.** Four brackets: y_char at the a-scale (mi_a) or the exact
Kepler time-average ⟨1/r²⟩ (mi_t), each with the EFE tables (conservative) or bare
isolated ν (Milgrom's frequency-decoupling guess). Same seed-31 population and grids
as the stored MG run (3U) ⇒ identical Newton baseline; differences subtract exactly.
[data/stage4l_summary.txt](data/stage4l_summary.txt).

| model | α̂ (simple/BE) | ΔlnL(Newton) s/BE | vs MG s/BE |
|---|---|---|---|
| **MG (3U, stored)** | 1.17 / 1.11 | **+108.7 / +98.8** | — |
| MI-t (EFE) | 1.40 / 1.55 | +97.6 / +89.4 | −11.1 / −9.4 |
| MI-a (EFE) | 1.18 / 1.36 | +92.6 / +85.0 | −16.1 / −13.8 |
| MI-t (no-EFE) | 0.59 / 0.59 | +78.8 / +72.6 | −29.9 / −26.2 |
| MI-a (no-EFE) | 0.51 / 0.51 | +77.5 / +71.2 | −31.2 / −27.6 |

1. **MG wins all 8 contests.** Ordering: MG > MI-EFE (−9…−16) > MI-no-EFE (−26…−31)
   > Newton (−71…−109). The 2D (ṽ×γ) shape wants the LOCAL, phase-dependent boost —
   the apocenter-boosted structure of the radial population (the same physics as the
   4J islands) is where per-orbit-global MI underperforms even with α̂ free (it
   inflates to 1.4–1.55 chasing the shape and still loses).
2. **The data demand the EFE.** Within MI, EFE-on beats EFE-off by 14–19 lnL, and the
   no-EFE variants localize α̂ ≈ 0.5 sharply — the data force the isolated boost DOWN
   to the EFE-suppressed amplitude. Milgrom's frequency-decoupling guess in its naive
   (fully isolated) form is rejected by the binaries themselves.
3. Newton loses under every model class yet tested (MG + 4 MI brackets): the anomaly's
   model-class robustness extends again. w_rad = 0.20 selected in all 8 fits (11th and
   12th consecutive model variations). BE−simple ≈ −6…−8 everywhere: the within-branch
   stalemate is formulation-independent.
4. **The paradox, fully quantified on one dataset + one ephemeris:** MG fits the
   binaries best (+9…+16 over the MI brackets) but predicts the Cassini-forbidden
   quadrupole (4K, ×4.3); the MI brackets need not produce the quadrupole but fit the
   binaries worse; naive no-EFE MI is dead. Surviving corners: (i) EFE-screened MG
   (TODO #8's sharpened target: boost with Q₂ ≤ 9e-27), (ii) time-NONLOCAL MI that
   mimics MG's phase structure (not yet constructed by anyone; our per-orbit-average
   brackets are its adiabatic limit), (iii) an unmodeled systematic (Newton still
   loses by ≥71 under every class — but 4K keeps this live).
5. Caveats: single seed (the 3U cross-seed scatter on ΔlnL was ±4–9; differences share
   the population so are stabler — the 6-seed budget is queued as 4L-b); the two
   prescriptions are representative members of an unbounded MI family, labeled as
   such; iso-variant α-grid coarse near 0.5 (interior, sharp profiles).
6. Priority (scout, provisional): the binary-eccentricity MI-vs-MG test appears to be
   FIRST-OF-KIND (nearest art: Paci+ 2020 rotation-curve MG-vs-MI, MG favored 6.9σ —
   convergent verdict from an independent observable; cite).

**Credences after 4K+4L:** anomaly-is-real ~70% (unchanged). Conditional on real:
EFE-screened-MG-like ~40%, nonlocal-MI ~30%, neither/unknown ~30%. BE microphysics
~15%. The theory space is now pinched from BOTH sides by our own results — which is
exactly what a measurement program is for.

## Stage 4L-b: the seed budget — CORRECTION #10: "MG wins all 8" does not survive; EFE-respecting MI ties (2026-07-23, [calcs/stage4l_mi_runner.py](calcs/stage4l_mi_runner.py) seeds 101–505)

All four MI brackets re-fit on the five remaining seeds; stored MG twins subtracted
per seed ([data/stage4l_summary.txt](data/stage4l_summary.txt)). Mean ΔlnL(MI − MG)
± SE over 6 seeds:

| bracket | simple | BE | verdict |
|---|---|---|---|
| mi_t (EFE on) | −3.5 ± 3.3 | −0.8 ± 2.5 | **STATISTICAL TIE with MG** |
| mi_a (EFE on) | −9.3 ± 2.3 | −3.1 ± 2.6 | mild MG lean (simple only) |
| mi_t no-EFE | −25.2 ± 3.4 | −19.5 ± 2.6 | dead |
| mi_a no-EFE | −27.9 ± 2.5 | −21.5 ± 2.2 | dead |

**CORRECTION #10: the seed-31 "MG wins all 8 contests" (4L) was realization luck —
seed 31 sat at the pessimistic edge (mi_t beats MG outright on seeds 303/404/505).**
The surviving statements, budget-grade: (1) the binaries DEMAND the external-field
suppression (no-EFE MI: −20…−28 mean, 12/12 seeds, α̂ collapsing to 0.51–0.63 sharply
interior — the isolated boost forced to the field-suppressed amplitude); (2) Newton
loses all 24 fits (+71…+108); (3) local-vs-global boost character: NOT distinguished
once the EFE amplitude is present (mi_t is a tie; only the a-scale prescription leans
MG, simple law only). α̂(mi_t) runs 1.4–1.7 (BE edge-pinned at seed 101 — flagged);
w_rad = 0.2 in 22/24 (two 0.3 at seed 101). **Consequence for the 4K paradox: the
Saturn-safe branch — modified inertia WITH proper external-field suppression — fits
the binaries as well as modified gravity. The paradox has an open door.** Caveat: our
"MI + EFE" imports the MG-derived suppression amplitude phenomenologically; a true MI
theory must generate that suppression from its own (frequency-based) EFE — what the
data pin is the AMPLITUDE, agnostic on mechanism. Credences (conditional on
anomaly-real): EFE-respecting MI (incl. nonlocal) ~40%, screened/plain-MG-with-
unknown-quadrupole-fix ~30%, neither/unknown ~30%.

## Stage 4M: the 4J residual closes — the 90° closest-approach arm is real model content; α untouched (2026-07-23, [calcs/stage4m_fly90.py](calcs/stage4m_fly90.py), TODO #2e-b closed)

Four v7 variants patching what the 4J autopsy prescribed, seed 31, vs the stored
baseline (+108.7/+98.8) ([data/stage4m_summary.txt](data/stage4m_summary.txt)):

| variant | ΔlnL vs baseline (s/BE) | α̂ (s/BE) |
|---|---|---|
| **fly90 (pure 90°-arm)** | **+4.7 / +3.4** | 1.11 / 1.08 |
| flymix (50/50) | −4.8 / −4.9 | 1.10 / 1.06 |
| erad (ceiling 0.9995) | −0.8 / +0.9 | 1.16 / 1.08 |
| both | −6.0 / −4.7 | 1.10 / 1.03 |

The pure closest-approach template WINS (+4.7/+3.4) — the data want the unbound
weight concentrated at γ≈90° in the ceiling band, exactly as 4J diagnosed (and NOT
split with the old low-γ asymptote arm; the e-ceiling was never the issue). **α̂
moves by ≤0.06 across all four variants (1.03–1.16, interior 8/8; w_rad = 0.2 in
8/8): the γ≈82°/ceiling residual is formally closed as a MODEL-SHAPE gap with zero
impact on the boost measurement.** TODO #2e-b done; the paper's §7 ceiling subsection
gains its model-side confirmation.

## Stage 4N: the reconciliation — two modeling choices manufacture two-thirds of the "16σ Newton" result (2026-07-23, [calcs/stage4n_banikstyle.py](calcs/stage4n_banikstyle.py), TODO #6)

Scout-verified context (primary IDs in the scout log): Banik+ 2024 (arXiv:2311.03436)
fit 540 (r_sky, ṽ) count cells on 8,611 pairs (2–30 kAU, <250 pc) with the
hidden-triple fraction FREE — it lands at **f_HT ≈ 69% under Newton** — and quote
16–19σ for Newton. Hernandez & Chae (arXiv:2312.03162) name three defects: the
fitted f_HT exceeds every independent calibration (≤50%; OUR 3J photometry: 12.3%
overluminous, kinematic fence 0.1); ṽ bins narrower than the measurement errors with
no noise convolution in the comparison; and no deep Newtonian anchor (their window
STARTS at the 2 kAU transition). Chae's own analyses get boost 1.37–1.49, Newton
rejected 5.8–9.2σ. No side-by-side exists in print. Stage 4N runs the ablations on
OUR pipeline (Banik-STYLE, honestly labeled — not line-by-line), seed 31, stored
baselines +108.7/+98.8 ([data/stage4n_summary.txt](data/stage4n_summary.txt)):

| ablation | ΔlnL(Newton) s/BE | α̂ s/BE | reading |
|---|---|---|---|
| vtonly (γ channel removed) | +105.3 / +98.4 | 1.18 / 1.55 | detection INTACT; the measurement degrades (BE α̂ 1.11→1.55, w_rad→0.3 both laws) — **the direction channel protects the parameters, not the detection** |
| freecomp (fence 0.1→0.8) | +42.8 / +43.8 | 0.71 / 0.83 | **THE MECHANISM: unfencing the companion fraction absorbs ~60% of the Newton deficit** — while requiring fractions the photometry forbids (3J: 12%; Banik's fit: 69%) |
| banikproxy (both + 2–30 kAU, anchor bin dropped) | +37.5 / +29.8 | 0.68 / 0.73 | **~2/3 of the significance manufactured away; α̂ biased low (~0.7)** |

Conclusions:
1. **The detection never flips.** Even under full Banik-style freedom our data retain
   ΔlnL = +30–38 for the boost — Newton does not win here under any ablation.
2. **The disagreement's anatomy, measured:** the unfenced companion fraction is the
   dominant manufactured component (~60 lnL); dropping the 0.2–2 kAU anchor + going
   ṽ-only adds ~5–14 more and biases α̂ to ~0.7. The residual distance to an actual
   "Newton wins" verdict plausibly lives in the H&C-documented noise-binning defect
   (sub-error bins, models compared unconvolved — our pipeline ALWAYS convolves
   per-pair noise, so that leg cannot be honestly ablated here) plus the Stage-3A
   realization systematic.
3. **The fix is measurement, not argument:** companion fractions are photometrically
   boundable (3J) and deep anchor bins exist. An analysis that uses both cannot land
   on Newton with this catalog.
4. Bonus: vtonly shows the γ dimension is what pins α near the parameter-free 1
   (without it BE's α̂ wanders to 1.55) — the direction channel is the field's
   missing instrument, quantified.
Caveats: single seed; α-grid coarse near the low α̂ of the ablated fits; the
freecomp Newton fits presumably sit at fence-forbidden fcomp (cube stores it;
summary prints wr only — noted). TODO #6 EXECUTED at proxy grade; a line-by-line
Banik reproduction remains future work if referees demand it.

## Stage 4O — the verification day (pre-circulation pass, 2026-07-22)

Purpose: every novelty claim and flagged quote in PAPER.md checked against primary
sources before the Zenodo tag (nine Haiku scout agents + one direct fetch; arXiv /
ar5iv / INSPIRE; verbatim retrieval only, NOT-FOUND reported as such). Author line
set: **Filip Hájek (independent researcher)**; CITATION.cff + .zenodo.json added so
the DOI record carries the paper title and author regardless of repo metadata.

**CORRECTION #11 (crediting, logged):** the √2 ceiling is the FIELD'S founding
bound — P&S 2018 (1711.10867: "the well-known result that u_3D < √2 for any bound
orbit"), Banik & Zhao 2018 (1805.12273: ṽ "can't exceed √2 in Newtonian gravity. In
MOND, we expect the upper limit to be somewhat higher"), and P&S 2023 (2205.02846)
name the band above the Newtonian decline "the key discriminant between models";
Banik 2019 (1902.01857) uses the above-√2 population as a contamination diagnostic.
Our 4J language ("the test itself is new") was overbroad. What IS new: the
perpendicular-column census and the edge-termination statistic. §7.2 now leads with
the crediting; Appendix A gains item 11; all "ten corrections" counts → eleven.

**Phantom citation caught:** "Paci et al. 2020" never existed — arXiv:2001.03348 is
**Petersen & Lelli 2020, A&A 636, A56** ("A first attempt to differentiate between
modified gravity and modified inertia with galaxy rotation curves"; MG mildly
favored within 1.5σ, 15 galaxies). Also found and now cited: Chae 2022 (2207.11069;
rotation curves disfavor MI at 6.9σ — complementary regime to our binary MI tie),
McCulloch & Lucio 2019 (1908.01434; quantised inertia on binaries — nearest miss to
the first-of-kind claim, different framework, no MI-vs-MG contrast), Milgrom 2023
(2310.14334: MI predicts "possibly a stronger external-field effect" on very wide
binaries — the qualitative anticipation of our no-EFE-dead result).

**Verified verbatim:** Chae 2023 (2305.04613) AQUAL-external-field sentence quoted
in §6.2; his numbers are γ = 1.43 ± 0.06 at 10σ — the 5.8–9.2σ our draft carried is
Chae 2024 (2402.05720, γ_g = 1.37), intro now quotes both correctly (our ref list
had the wrong 2023 arXiv ID, fixed). P&S 2025 triple-population caveat exact;
byline gains a third author (Pittordis, Sutherland & Shepherd, OJAp). 2304.00576 =
Jones-Smith & Mathur (AJ) — flag closed. DHF 2024: "an 8.7σ tension under fiducial
model assumptions" + both escapes verified (M/L → 7.2σ; bulge removal → 2.7σ) +
their own remark that Cassini already beats the wide-binary constraint — nobody has
run the binary-CALIBRATED direction; ours stands.

**Priority sweeps:** C&T citation chain (INSPIRE recid 1731589, 31 citations, 24
retrieved): no empirical test of their form anywhere in it, no author follow-up
fitting data — "apparently the first empirical test" SUPPORTED (7 unfetchable
entries disclosed in §2.2). Simple-ν self-consistency identity: no prior derivation
found; Famaey & Durakovic 2025 (2501.17006) review states "There is no existing
clear derivation of such a transition from first principles" — now cited as the
anchor; the "Timeflow" flag resolved = Trofimov 2026, CQG 43, 135020, a vacuum
phase-interference mechanism, distinct from our identity; its paywalled full text
is the one unexecuted read (disclosed in Appendix A). Realization systematic:
re-scouted, still no published quantification. Reference pass: 10 arXiv IDs
verified, 8 missing citations resolved (QUMOND = Milgrom 2010, MNRAS 403, 886;
El-Badry 2101.05282; Raghavan 1007.0414; Pecaut & Mamajek 1307.2657; Chae & Milgrom
2022 = 2201.02109; Bekenstein & Milgrom 1984, ApJ 286, 7; P&S/B&Z volumes); byline
fixes: Hernandez, Chae & **Aguayo-Ortiz** 2024 (third author we had dropped), Pazy
& Argaman is **2012**, PRD 85, 104021.

**De-orphaned:** Vokrouhlický, Nesvorný & Tremaine 2024 (2403.09555) — AQUAL fails
Oort-cloud comet binding energies and the detached disk; integrated into §8.1 as
the population-dynamics bracket of our ephemerides veto (their escape clause =
screened theories, same as our §8.2). **New arrival flagged:** Desmond et al. 2026
(2602.24035), quality-framework WB analysis, "No Evidence for MOND" — cited in §1
with an explicit pending-full-read label; the full methodological read is a new
TODO item that must close before arXiv.

PAPER.md → **v1.1** (verification pass stamped in the header). COLLEAGUE-BRIEF.md
rewritten to the final numbers (still uncommitted by design). Verification scripts:
none (literature pass); agent reports summarized here are the record.

## Stage 4P — the Cookson et al. 2026 full read (the pre-arXiv gate, 2026-07-22)

**Byline correction first (notebook discipline):** the 4O entry called arXiv:2602.24035
"Desmond et al. 2026" — that byline was a scout-agent hallucination (the third phantom
byline caught this cycle). The full read (PDF, all 22 pages) gives: **Cookson, Banik,
El-Badry, Sutherland, Penoyre, Pittordis & Clarke 2026, MNRAS** — lead author Stephen
A. Cookson, *Independent Researcher, Crawley, UK* (the co-author list spans the
Newton-favored wing of the field). PAPER/TODO/CLAUDE fixed; the error stands here in
4O as logged.

Content: a quality checklist for the WBT + a median-ṽ flatness test on 1,421 RV-clean
DR3 pairs within 130 pc (both stars need RVs; ΔRV<10; RUWE<1.25; HR parallelogram;
ipd_frac_multi_peak=0; ṽ<2.5; degrouping; 1–30 kAU). Result: flat medians, Newton
χ²=2.48 vs MOND 15.10 (5 dof) → ~500–1500× likelihood ratio (2.6–2.7σ equivalent).
Their meta-analysis: MOND signals shrink as checklist scores rise. Useful gifts:
g_N,e = 1.184a₀ (their eq 4) — independently confirms our 3T conversion (1.15±0.05)
to 3%; exact volumes for several refs; and their §7.2 demonstration that omitting
the **spherical projection correction** inflates median ṽ by ~0.15 beyond r_M.
Their vulnerabilities (for §7.4): the flatness test assumes e-distribution constant
with s while citing Hwang+22 — who MEASURED it to steepen (uniform→superthermal);
the measured trend's sign SUPPRESSES a step (apocenter lingering) rather than faking
one; no forward model by design. Their "counting ṽ>√2 is not a practical MOND test
(contaminants)" = precisely the objection our perpendicular column answers.

**The gate finding: our pipeline never applied the spherical projection correction
(their item 4.1.6) — no mention anywhere in calcs/ or NOTES.** Triggered 4Q immediately.

## Stage 4Q — the perspective audit (correction #12, 2026-07-22)

The omitted term: Δv_spur = −RV_sys·θ·ŝ (receding pairs appear to shrink) — purely
RADIAL in the (ṽ,γ) plane, growing ~s^1.5 in ṽ units: shaped like the signal, aimed
at the w_rad arm. [calcs/stage4q_perspective.py](calcs/stage4q_perspective.py),
output data/stage4q_perspective.txt; gates G1/G2/G3 all PASS.

- **Exposure (G1):** median κ=25 km/s·θ/v_c by bin: 0.0006 (anchor — immune), 0.006,
  0.029 (6–20), 0.119 (20–50 kAU; 84th pct 0.22). The anchor bins cannot carry it;
  the widest bin genuinely could.
- **Q1 (Newton+perspective-only):** predicted anchor ratio 1.016 vs observed 1.086 —
  the effect can supply ≤1.6 of the 8.6 points.
- **Q2 (component split — the kill test):** ṽ ratio 1.086 (reproduces 2C exactly);
  **ṽ_perp (immune) 1.151 (CI 1.115–1.197)** — the immune component boosts MORE, the
  OPPOSITE of the artifact signature — monotonic per-bin rise 0.298→0.301→0.332→0.384
  against the radialization headwind; ṽ_rad 1.052.
- **Q3 (direct correction, catalog RVs ~100% on this bright subsample):** slope of
  observed widening vs −RV·θ = **0.923 (CI 0.795–1.082; predict 1)** — Cookson's
  fig 7 reproduced on our own selection, the systematic is REAL and now measured;
  exact per-pair correction: anchor **1.086 → 1.078 (CI 1.052–1.103)** — the honest
  haircut. G2 injection round-trip passes to <0.5%.
- **Ceiling addendum** ([calcs/stage4q_ceilingcheck.py](calcs/stage4q_ceilingcheck.py)):
  independent reimplementation of the γ machinery counts 10 raw in-band vs 4J's 11
  (boundary-convention fuzz, ±1); correction moves ONE pair 1.662→1.675±0.044 (edge
  pile-up, cliff-consistent) and tightens another to γ=89.8°; corrected census 9–10,
  leakage null degrades to no worse than ~1e-8. Conclusion intact.
- α exposure: bounded by the 1.6%-of-ratio contribution — inside ±0.11 and the g_ext
  systematic. **Corrected-velocity v7 re-fit queued as TODO #2i (pre-arXiv hardening).**

Verdict: the systematic is present, measured, and too small by 5×; the immune-component
boost is stronger than the full-statistic boost. The anomaly stands at **1.078
(CI 1.052–1.103)** on the corrected anchor. Correction #12 logged; PAPER → v1.2 with
new §7.4 (checklist mapping, audit, and the e-trend/forward-model reply to their null).
Twelve corrections now in Appendix A.

## Stage 4S — the zero-point coefficient measured (2026-07-23)

User directive: stop treating the ½ as a contest, make it a measurement ("if this
thing is genuinely new, that must be our opus magnum"). Instrument: promote c₁ to a
continuous parameter via ν_λ = (1−λ)·ν_standard + λ·ν_RAR, c₁ = λ/2 EXACTLY on the
slice (G1: numerical series extraction reproduces λ/2 to 2e-3 at λ=0/½/1; c₂ rides
along: 0.25 → 0.167 → 0.0833). Grid λ ∈ [−0.30, 1.50]; (a₀, f_ML[, s_int, δ_lens])
profiled per node, warm-started; 200-galaxy-bootstrap × full 5-param refit.
[calcs/stage4s_c1fit.py](calcs/stage4s_c1fit.py) → data/stage4s_c1fit.txt.

**Marginalized joint SPARC+lensing (the honest instrument): λ̂ = 0.900, INTERIOR,
parabolic profile; ĉ₁ = 0.450, profile interval 0.385–0.519; galaxy-bootstrap
ĉ₁ = 0.427 +0.290/−0.246.** Placements: c₁ = ½ at Δ(−2lnL) = 0.5 (0.7σ — dead
center); c₁ = 0 at +56.3 (7.5σ profile; 95.5% of bootstraps ≈ 1.7σ one-sided);
c₁ = ¼ at +9.7 (3.1σ profile; 73.5% of bootstraps above ¼). c₂(λ̂) = 0.100 —
between BE 1/12 and simple 1/8 (rung 2 still unresolvable, §4.5). G2 endpoint
regression: λ=1 → −8397.72, λ=0 → −8341.95, both EXACT vs the stored 4E/4F fits.

Two disclosures logged with the number:
1. **Profile-vs-bootstrap gap (7.5σ vs 1.7σ)** = galaxy-population variance — the
   SPARC sibling of the 3A realization systematic. Bootstrap quoted as primary.
2. **Raw-χ² FAILS STRUCTURALLY on continuous families**: λ̂ edge-runs to the grid
   boundary in both windows (Δχ² thousands; bootstrap pinned at the parameter
   bound) — at χ²/dof ≈ 57 the unmodeled scatter rewards unbounded shape-bending.
   The 4B/4E deflation was the warning; this is the demonstration. Raw χ² can rank
   discrete families; it cannot measure a coefficient. (No correction number — the
   raw treatment was never quoted as a coefficient measurement — but the §4.2
   caveat is now a theorem-by-example.)

Paper: new §4.3 (old 4.3/4.4 → 4.4/4.5), abstract (2) upgraded, conclusions
bullet upgraded. The ½ now enters the paper as a measured quantity with the
prediction sitting 0.7σ from center — the coefficient-level spine the user asked
for. Rung 2 remains the wall (0.1σ); no enthusiasm spent pretending otherwise.

## Stage 4T — the bath's second moment (2026-07-23)

The Planck-oscillator reading's first NEW falsifiable channel beyond the mean:
Var(n) = n(n+1) per Bose mode ⇒ the RAR's intrinsic scatter should carry the
shape σ_rel(x) = √(n/(n+1))/√N — rising into deep MOND, dying ~e^(−x/2) on the
Newtonian side, amplitude 1/√N = the effective mode count. Four scatter models
fitted jointly with (a₀, f_ML) on 2,700 SPARC points (per-point errors always
in): constant / pure oscillator / oscillator+floor / six free x-sextiles.
[calcs/stage4t_bathnoise.py](calcs/stage4t_bathnoise.py) → data/stage4t_bathnoise.txt.

Results (BE; simple identical):
- Scatter is NOT constant: free 6-bin beats constant by Δ(−2lnL)=42.7/5 params;
  monotone decline 0.144 dex (deepest sextile) → 0.107 (most Newtonian).
- Floorless oscillator DEAD (+382): the data demand a floor everywhere.
- **Oscillator+floor: −25.1 for ONE extra param over constant; N̂ = 21.5 modes
  over a 0.101-dex floor.** Floor-free bound from the deepest bin alone: N ≥ 7.
- Residual: a bump at x≈1 keeps free-shape ahead of osc+floor by 17.6/4 params —
  exactly where per-galaxy M/L scatter projects maximally onto g_obs. Flagged as
  the probable mundane term; the M/L-marginalized second-moment fit = the
  required next instrument (queued as the 4U candidate).
- Gates: G1 nesting PASS; G2 injection (truth N=30) recovers N̂=30.4 PASS.
- Framing enforced: CONSTRAINT, not detection. Caveat: per-point independence;
  a galaxy-level correlated draw loosens N.

Paper: new §4.6, §9 falsifier #6, abstract clause, conclusions bullet; §1 gains
the "one-half threads the paper" spine paragraph (user directive: the ½ is the
through-line — executed structurally, not by rebuilding the measurement paper
around a 1.7σ-bootstrap-grade exclusion; the dedicated ½ paper is the post-DR4
cycle). Scout launched on prior continuous-c₁ fits + scatter-vs-acceleration
literature (the §4.3/§4.6 novelty wording is scout-conditional; "first" is not
printed until it returns).

**4T addendum — scout returned (2026-07-23): BOTH channels open.** (1) c₁: no
published fit of the deep-expansion coefficient; nearest art = DHF24's δ/γ/n
families (continuous, but transition-SHARPNESS — same axis as our §3 p; their
family definitions fetched and checked) — §4.3 now credits and distinguishes.
(2) Scatter: Lelli+17 = global scatter, residuals property-uncorrelated;
Desmond 2023 (MNRAS 525, 6130) = hierarchically-marginalized intrinsic scatter
0.034±0.002 dex — as a ceiling on the thermal term this sharpens the mode bound
to N ≳ 10² (added to §4.6); NO acceleration-binned intrinsic-scatter
measurement and NO fluctuation-statistics reading found anywhere. **Phantom #4
caught:** the scout attributed a covariant-error quote to arXiv:1907.04501; a
direct fetch shows that ID is Tian & Ko 2019 "Halo Acceleration Relation" —
no such content. The covariant-error point (distance/inclination errors enter
both RAR axes and project harder where the relation flattens) is retained as
our own unattributed methodological observation in §4.6; the phantom ID is not
cited. Four phantom attributions this cycle, all caught by direct fetch —
primary-sources-only discipline validated a fourth time.

## Stage 4R — the corrected-velocity budget (correction #12 executed, 2026-07-23)

TODO #2i run to completion instead of footnoted (user directive: "they are gonna
stop reading the damn thing if there's a footnote that we will fix something" —
correct, and now moot). Patch-runner [calcs/stage4r_corrected_refit.py](calcs/stage4r_corrected_refit.py)
applies the 4Q spherical-projection correction at the stage3p data build (both ṽ
and γ inherit it; RV column names resolved dynamically — the FITS carries
dr2_radial_velocity*), same six seeds as the stored 3U baseline, physical-field
tables, both laws. Output data/stage4r_summary.txt; baseline untouched.

**Paired per-seed verdict (corrected − baseline):**
- simple α̂: Δ = −0.012 ± 0.012 (per-seed: +0.01, −0.02, +0.04, −0.03, −0.03, −0.04)
- BE α̂:     Δ = −0.022 ± 0.015 (−0.01, −0.04, 0.00, +0.03, −0.08, −0.03)
- Newton ΔlnL cedes ~5 (≈5%): corrected +102.8 ± 9.2 (simple) / +90.1 ± 6.5 (BE),
  minimum across 12 fits +83.6 — exclusion untouched.
- w_rad = 0.2 in 12/12; α̂ interior 12/12.
- **BE-minus-simple lean UNCHANGED: −12.6 ± 2.4 (was −12.3 ± 2.2)** — the binaries'
  classical-bath lean survives the perspective correction.

Verdict: exactly what the 4Q audit bounded — the correction is real, small, and
changes no conclusion. §6.3 final numbers (α = 1.18 ± 0.11 / 1.13 ± 0.13) STAND,
now correction-executed; §7.4(b) queued-language deleted; Appendix A item 12
closed with the shift numbers; PAPER → v1.3. (A corrected-data 3V-style 1000-
bootstrap would move the center by ≤0.02 against ±0.11 quoted — completeness
theater, not run; stated here for the record.)

## Stages 4U + 4V — the hierarchical second moment, the scorecard, and §2.4 (2026-07-23)

User redirected: publication DEFERRED, oscillator = the program's front. Three
instruments executed same-day.

**4U ([calcs/stage4u_mlmarg.py](calcs/stage4u_mlmarg.py)):** per-galaxy disk-M/L
offsets (0.1-dex lognormal prior) profiled jointly with globals; scatter models
recompete on the residual. **The thermal-direction trend SURVIVES marginalization**
(osc+floor −22.8/1 param vs constant; floor drops 0.101 → 0.059 dex; N̂ 21.5 → 62.9
— 4T was floor-limited; constraint now quoted N ~ 20–60, nuisance-depth-dependent).
No galaxy-level draw (corr(|δ̂_g|, ⟨x⟩_g) = −0.02). std(δ̂_g) = 0.168 dex > the 0.1
prior (δ absorbs distance/inclination too — expected). DISCLOSURES: (1) the
marginalized free shape is V-SHAPED — s_b = [0.105, 0.078, 0.055, 0.027, 0.065,
0.084] — the Newtonian-side rise sits where BULGE M/L (not marginalized; only disk
was) surfaces; (2) **G2 injection gate FAILED at r = 0.74 vs the 0.8 target** —
partial recovery; therefore cross-model comparisons that re-optimize δ_g (the M2h
−349 and the M3h area-scaling −229, σ_A = 2.4 modes/kpc²) are REPORTED AS
UNRESOLVED, not measured. G1 prior→0 reproduces 4T exactly (PASS). Decisive next:
the full hierarchy (bulge+distance+inclination) = TODO O1.

**4V ([calcs/stage4v_scorecard.py](calcs/stage4v_scorecard.py)):** the scorecard +
the binary a₀ translation. κ = dln(B−1)/dln a₀ from the g-scan tables: 0.784
(simple) / 0.916 (BE), gates pass (α=1 → 1.2e−10 exact). **a₀(binaries) =
(1.48 ± 0.18) / (1.37 ± 0.17)e−10** — the disclosed α>1 lean in temperature
language: +2.5σ/+1.9σ vs Planck cH₀/2π, +2.0σ/+1.4σ vs SH0ES — the reading's
sharpest internal tension, carried openly (either g_ext calibration owns it or
the deep amplitude is not fully H₀-locked → a₀∝H(z) adjudicates). Galaxy-side
rows: +0.1σ / −0.5σ vs Planck — the temperature check passes where it is
cleanest. Five dials within ≲1.6σ. Paper: §9 restructured (9.1 scorecard table,
9.2 falsifiers).

**§2.4 written (the mechanism specification):** six measured requirements any
construction must clear. Unruh scout (primary sources): Milgrom 1999 (astro-ph/
9805346; T(a)−T(0) "depends on a in the same way that MOND inertia does";
mechanism "still a far cry off") and Deser & Levin 1997 (gr-qc/9706018; 2πT =
√(a²+Λ/3)) VERIFIED and now cited; the OCCUPATION-RATIO reading (x = T_U(g_obs)/
T_dS, exact in the deep limit) and the geometric-mean frequency ω = √(g_N a₀)/c
= g_obs/c: NOT FOUND anywhere — stated in §2.4 as observations, not derivations.
(Scout's sloppy "nobody wrote ν=1+n_BE" side-claim discarded — it contradicts
C&T Eq. 23, which correction #7 established; only the Unruh-ratio FRAMING is
unclaimed.) Nearest adjacent: arXiv:0908.4239 (Unruh-like MOND), cited by ID.

## Stage 4W — the full hierarchy: an identifiability boundary (2026-07-23)

O1 executed ([calcs/stage4w_fullhier.py](calcs/stage4w_fullhier.py)): per-galaxy
disk M/L + bulge M/L (31 bulge galaxies) + a VERTICAL offset with MEASURED priors
(SPARC's own e_D/D and e_inc: both move g_obs at fixed g_bar under our loader;
median σ_v = 0.097 dex). Scatter contest on the residual.

**Verdict: NOT a survival, NOT a refutation — a demonstrated identifiability
boundary.** With vertical freedom on, osc+floor LOSES to constant (+24.4), N̂
inflates to 189 (amplitude ≈ 0), floor → 0.035 dex (≈ Desmond-2023's 0.034 —
our hierarchy converges to his intrinsic number, good cross-check). BUT the
CALIBRATED injection gate FAILS decisively: slope_obs 0.225 vs shrinkage-expected
0.931 — the machinery itself proves per-galaxy M/L and vertical offsets are
degenerate on SPARC (a disk-M/L shift IS mostly a vertical shift for
disk-dominated galaxies at low x), and any between-galaxy part of a real
x-shaped signal is eaten with them. G1 (priors→0 → 4T repro) PASS.

Honest chain across the ladder: 4T raw trend (real) → 4U survives disk-M/L
(−22.8) → 4W not identifiable vs distance/inclination freedom. The channel's
decision needs EXTERNAL distance anchors (Cepheid/TRGB — Desmond's inputs) or an
independent dataset. Paper §4.6 rewritten to say exactly this; scorecard N row
= "20–60 if thermal; identity unresolved"; §2.4 item 4 conditioned.

**The surviving fact — the transition bump:** free bins STILL beat constant
(−84) with all offsets free, and the x≈1 excess persists at 0.054 dex vs
0.026–0.035 neighbors — it survived disk-M/L (4U), bulge-M/L AND vertical
freedom (4W). NOT bulge, NOT distance/inclination, NOT the monotone thermal
term. New hypothesis logged: environmental EXTERNAL-FIELD scatter projects onto
g_obs maximally at the transition — exactly where the bump sits. O4 (the
environment-split RAR) inherits a sharp, falsifiable target: the excess must
CORRELATE with per-galaxy g_ext estimates. M3w area-scaling verdict remains
unresolved (same G2 caveat as 4U, now demonstrated rather than suspected).

## Stage 4X — the binaries read the same dial (2026-07-23)

O2 executed ([calcs/stage4x_binlam.py](calcs/stage4x_binlam.py)): six λ-family
tables through the QUMOND solver at e_N=1.2 (gates: λ=1 reproduces the stored
BE table <1%; λ=½ isolated-sphere identity <2%), then the v7 likelihood
(perspective-corrected data, all nuisances, α refit per node) profiled over λ,
seeds 31+101. Absolute best-lnL per (λ, seed) → data/stage4x_summary.txt.

Combined profile (Δ from max): λ=0: −41.1; 0.25: −13.0; 0.50: −4.6; 0.75: 0;
1.00: −4.7; 1.25: −9.5. Per-seed peaks: seed 31 at λ=1.00, seed 101 at 0.75 —
realization scatter ≈ the grid step. **Binary ĉ₁ = 0.37–0.50** (curvature σ_λ ≈
0.08; honest span quoted). **c₁ = 0 rejected at ΔlnL ≈ 20 PER SEED (−23.0/−20.5)
— sign-robust**; and at λ ≤ 0.25 the fits ride the α grid edge (α̂→2,
interior=False) and STILL lose = shape rejection, not amplitude trade. At λ ≥
0.5, α̂ interior (1.1–1.4, consistent with the standing lean). w_rad = 0.2 in
12/12.

**The program's strongest single oscillator statement: two disconnected
systems — rotation curves (ĉ₁ = 0.43 +0.29/−0.25, 4S) and wide binaries at the
external-field transition (0.37–0.50, this stage) — independently return the
same leading coefficient, each consistent with ½, each excluding 0 on its own.**
Caveats logged: 2 seeds (full budget queued if wanted), Δλ=0.25 grid, g_N,ext
systematic shared with 6.3, family-conditional as 4S. Paper: new §6.4 (old
6.4→6.5, cross-refs fixed), abstract (3) + conclusions + scorecard row added.

## Stages 4Y + 5A — the anchored subsample and the bump's identity (2026-07-23)

**4Y (O1b, [calcs/stage4y_anchored.py](calcs/stage4y_anchored.py)):** e_D/D ≤ 0.10
subsample = 43/153 galaxies, 695 points; σ_v halves (median 0.044 vs 0.097) —
anchoring works. Verdicts: thermal term NOT required there (M1b−M0 = +0.02,
N̂→∞); identifiability IMPROVES but is not restored (calibrated injection slope
0.42 vs expected 0.93 — cf. 0.225 full-sample; corr 0.80) — 43 galaxies is
underpowered, gate honestly FAILED again; deep bins mildly elevated
(0.049/0.041 vs 0.020 mid); NEW oddity: Newtonian-end spike s_b[5] = 0.116 on
the anchored subsample (few-galaxy inner/bulge points, 7 bulge galaxies —
flagged, not interpreted). Conclusion: the unlock needs MORE anchors, not
better machinery — parked until a distance-ladder expansion.

**5A (O4a, [calcs/stage5a_bumpid.py](calcs/stage5a_bumpid.py)):** the template
contest DISQUALIFIED ITSELF — G1 injection FAILED in the informative way:
injected per-galaxy p-scatter produces NO transition bump (within one galaxy
t_p is nearly constant in x → the vertical channel eats it; dp̂ recovery 0.66).
Per-galaxy templates cannot answer the bump question, PROVEN. The residual
value: (1) **the bump survives EVERY per-galaxy channel (vertical, disk M/L,
bulge M/L, sharpness, environment) → it is POINT-LEVEL structure at the
transition** — no galaxy-wide nuisance can make it; candidates now: beam
smearing at the RAR knee, intrinsic transition structure, or (thermal reading)
the transition shell. (2) The deep-weighted EFE template DOES absorb the deep
bins (0.063→0.045) — environmental scatter stays alive as the thermal term's
rival for the monotone trend; only real per-galaxy environment estimates
resolve it (agent hunting the Chae 2020/2021 tables). (3) Template overlap
corr(t_p,t_e) = −0.61 (contest would have been valid had the channels been
identifiable). Paper §4.6 + §9.2 #6 updated; hypothesis-churn logged openly:
EFE-at-transition (wrong, corrected same day) → sharpness-scatter (instrument
impossible, proven) → point-level (current, characterization not yet identity).

## Stage 4Z — the hierarchical c₁ RELOCATES the profile (2026-07-23)

O3 executed ([calcs/stage4z_hierc1.py](calcs/stage4z_hierc1.py)): the 4S λ-profile
with per-galaxy disk-M/L offsets (0.1-dex prior) profiled at every node; 100-rep
galaxy bootstrap with full joint refits. **BOTH GATES PASS** (endpoints reproduce
4S exactly; injection at λ-truth=1 WITH offsets recovers 1.032 — the machinery
does not bias λ̂ down).

**The result reshapes the ĉ₁ claim: profile λ̂ = 0.516 → ĉ₁ = 0.258 (interval
0.208–0.309)** — vs 4S's flat-M/L 0.450 (0.385–0.519). The two PROFILE intervals
don't overlap: the profile estimator is nuisance-model-dependent at ±0.2, and
since the injection gate passed, that spread is REAL per-galaxy structure doing
work in the global fit, not machinery. Bootstraps agree across treatments
(0.427 +0.290/−0.246 flat; 0.377 +0.234/−0.311 hier) — population variance
dominates both. c₁=0: +28.9 profile (5.4σ), P(λ>0)=0.89. c₁=½: +17.6 from the
hierarchical peak — **the hierarchical galaxy profile peaks AT the ¼ value.**

Honest synthesis (paper updated everywhere — abstract, §1 spine, §4.3, §2.4,
§6.4, scorecard, conclusions): **galaxy ĉ₁ ∈ 0.26–0.45 (profile, treatment-
spanning), ≈ 0.4 ± 0.3 (bootstrap); c₁ = 0 excluded in EVERY treatment (Δ 29–56,
≥89% bootstraps); the second digit — ¼ vs ½ — is OPEN on the galaxy side;
binaries independently 0.37–0.50.** The 4S headline "ĉ₁ = 0.450, ½ at 0.7σ" was
flat-M/L-conditional; corrected the same day it was minted — the dial stays
measured, its second digit does not. The §4.4 ¼-FUNCTION contest (c₂=7/96
function, not the λ-slice) is now labeled flat-M/L-conditional; re-running the
bath-matrix contest under hierarchical M/L = new TODO O6.

## Stage 5B/5E — the environmental control: THERMAL SURVIVES (2026-07-23)

O4 executed. Chae et al. 2021 (ApJ 921, 104) Table 3 — per-galaxy environmental
Newtonian fields log₁₀ e_N,env for 109 SPARC galaxies (SDSS footprint), "max
clustering" and "no clustering" columns — extracted from the PDF
(data/chae2021.pdf → data/chae2021_text.txt; the two columns are a near-constant
0.9-dex offset apart, so one global amplitude β × the max-clustering pattern
spans both). Matched 94/153 of our kept galaxies (1347 points, 52% of the deep
bins — proportionate). EFE model: **the exact 1D QUMOND/AQUAL collinear ratio,
their Eq. (2)** (ν_e = ½ + [D−C]/2y), applied as a suppression template on the
BE mean; formula gated (e→0 identity 9e-16, y→0 cap analytic, monotone).
NOTE: 5A's tanh-form EFE template (Chae–Milgrom fitting-function style) had an
ambiguous sign in its curvature term as we transcribed it; 5A's conclusions are
unaffected (its per-galaxy template coefficients were free-sign) but 5B/5E use
the exact Eq. (2). Near-miss logged, no published number affected.

First pass (5B, [calcs/stage5b_envtest.py](calcs/stage5b_envtest.py)) exposed a
convergence trap in our own machinery: fixed 3-round coordinate descent leaves
up to ~80-lnL slop (the β=0 profile node sat 84 above the identical E0 model;
"scrambled" nested fits landed BELOW the base model — impossible at
convergence). Same trap quantified on the 5C side by the 4Z-vs-5C benchmark
(~10 lnL at 3 rounds, full sample). **All verdict-bearing numbers therefore
come from the convergence-hardened re-runs** (adaptive descent, tol 0.05, max
15 rounds, nesting-inequality gates): 5E
([calcs/stage5e_envtest_conv.py](calcs/stage5e_envtest_conv.py)), and 5D for
the bath matrix. Instrument gates all PASS: nesting inequalities hold; the
warm-chained β=0 node reproduces E0 (+0.32); injections recover β_true=3 →
2.97 and β_true=0 → 0.00; full-sample machinery reproduces 4U exactly (5B G1).

**The verdict: the real environmental pattern takes NO significant credit on
the hierarchical-M/L RAR.** β̂ = 0.044 with upper D1 ≈ 0.02–0.05; **β=1 (max
clustering at face value, collinear) is excluded at Δ(−2lnL) = +168**; even
β=0.25 costs +11. The tiny E2−E0 = −7.2 credit sits INSIDE the scramble
distribution (8 perms: 0.03…10.13, one scramble beats it) — generic freedom,
not pattern. Correlation channel (model-light): raw Spearman(deep residual,
log e_N) = −0.165 (p=0.20); depth-partialed −0.123 (p=0.34) — null. The x≈1
bump: untouched by the EFE term (expected — EFE dead at x≈1), still
point-level. Thermal credit on S: E1−E0 = −4.6 (vs −22.8 full-sample; S holds
52% of deep points, so the credit is sample-dependent beyond point counting —
the S galaxies' deep ends are less scatter-elevated; logged as a caveat, the
4T/4U claim was always CONSTRAINT not detection).

Read: **the deep-end scatter trend is NOT environmental at published
amplitudes — the thermal (oscillator shot-noise) candidate survives its
sharpest control to date.** Secondary finding: our pipeline sets a tight upper
limit on EFE visibility in SPARC deep ends (≲2–5% of the collinear
max-clustering amplitude, at fixed SPARC distances/inclinations with
hierarchical disk-M/L). Stated caveats: collinear = maximal geometry (random
orientation ~halves the effective e, softening the β=1 exclusion to roughly
β_eff=0.5, still +55); e_N log-errors 0.29 dex attenuate β̂ but cannot
manufacture a +168 rejection; our channel (point-level RAR with measured
distances) ≠ Chae's channel (within-galaxy outer-decline shapes with free
per-galaxy distance/inclination/e) — we do NOT adjudicate their detection,
we close OUR deep trend's environmental escape.

## Stage 5C/5D — the bath matrix FLIPS under hierarchical M/L (2026-07-23)

O6 executed. The 4F four-cell contest (+ dead-branch control) re-ranked with
per-galaxy disk-M/L profiled (4Z treatment, joint SPARC+lensing objective):
5C ([calcs/stage5c_hierbath.py](calcs/stage5c_hierbath.py)) found the flip at
3 rounds; 5D ([calcs/stage5d_hierbath_conv.py](calcs/stage5d_hierbath_conv.py))
confirmed it CONVERGED (adaptive descent; CONV gate: BE reaches −10435.00 vs
4Z's warm-chained λ=1 benchmark −10435.06 ✓).

**Converged ladder: boot −10510.6 < BE −10435.0 < standard −10424.5 < simple
−10336.3.** Under flat M/L (4F) the order was BE > simple > boot > standard;
hierarchically the quantum SELF-CONSISTENT cell — ν = 1 + n_BE(ν·y), c₁ = ¼,
c₂ = 7/96, e^−y screening, frequency set by the TOTAL acceleration (the §2.4
Unruh-ratio reading ω = g_obs/c) — **wins by 75.6 over BE and 174 over
simple.** Coherent with 4Z's continuous profile peaking at ĉ₁ = 0.258 ≈ ¼:
two independent instruments (mixture-family slice; exact function contest)
relocate to the same cell once per-galaxy M/L is real.

Truth-calibrated injections (all three truths recovered by the converged
machinery): BE-truth → BE wins (boot +38 behind); boot-truth → boot wins (BE
+98 behind); simple-truth → simple wins. **The observed real-data gap (−75.6)
sits near the boot-truth calibration (−98), far from the BE-truth one (+38).**
Galaxy bootstrap (50 reps, converged-warm): boot−BE = −69 ± 65, **boot better
in 43/50 (86%)** — a strong lean, not yet decisive (population variance
dominates, as everywhere in this program). Also notable: simple falls to LAST
(even standard beats it by 88) — the flat-M/L "simple-lean" of the galaxy side
does not survive hierarchy either.

Scout verdict (Haiku, arXiv/INSPIRE, direct-quote discipline): the implicit
self-consistent function ν = 1 + n_BE(ν·y) — **NOT FOUND anywhere**; no
self-consistent-argument variant of Cadoni–Tuveri exists; even the simple-ν
implicit identity ν = 1 + 1/(ν·y) is unpublished as a stated claim (consistent
with Famaey & Durakovic 2025's "no clear derivation" line we already cite).
Both stay "apparently unpublished" pending a deeper pass — no "first" printed.

Hier caveats, stated: the hier optimum trades a₀ down (boot 8.7e−11) against
f_ML up (1.64) — the a₀-vs-cH₀/2π scorecard row is defined on the fiducial
treatments and a hier-boot scorecard is future work; the binary side read
c₁ = 0.37–0.50 through the λ-MIXTURE family whose low-c₁ members carry
standard-μ's shape — boot is a different function at the same c₁, so the
binary cross-check needs boot's own EFE tables → executed as Stage 5F (O7).

## Stage 5F/5G — the binaries veto boot; the flip decomposes as a TAIL story (2026-07-23)

**5F (O7, [calcs/stage5f_bootbin.py](calcs/stage5f_bootbin.py)): the binaries
REJECT the boot function.** QUMOND-EFE boot table generated at g_N,ext=1.2a₀
(spherical-identity gate 0.01%; ν_boot(1)=1.350 vs BE 1.582); v7
perspective-corrected 2D fit, seeds 31/101, on the 4X grid footing. Result:
boot lnL = −56360.6 / −56369.1 — **+22.3 / +17.1 behind BE, +19.9 / +24.1
behind λ0.75, with α̂ riding the grid edge (2.00, non-interior) in both
seeds** — the same edge-riding-and-still-losing signature that rejected the
low-λ mixtures in 4X: a SHAPE rejection, not an amplitude trade. Newton stays
dead under boot (+75/+81). The x≈1 sample refuses boot's weak transition.

**5G (O8, [calcs/stage5g_tailtest.py](calcs/stage5g_tailtest.py)): what the
hierarchical galaxies actually vote for is the TAIL.** Observation: the 5D
ladder (boot < BE < standard < simple) is monotone in Newtonian-tail
sharpness (e^−y, e^−√y, ~y⁻², ~y⁻¹). Probe: the §3 screening family ν_p —
which IS the occupation law at p=½ — run through the converged hierarchical
machinery (G1: p=½ reproduces hier BE −10435.00 exactly): **p̂_hier ≈ 0.65,
gain −56.4 over p=½ — 75% of boot's −75.6 flip from the tail dial alone**,
with the transition kept at ½-branch grade (ν_p(0.65)(1)=1.423) and the deep
NLO structure off-boot (exponent 2p, disclosed). Convergent with the §3
flat-treatment measurement p = 0.578 ± 0.12: two treatments bracket
**p ≈ 0.58–0.65, sharper than the pure occupation tail** — and Cassini wants
sharp. The ¼-vs-½ deep digit therefore stays OPEN (boot's residual non-tail
gain is only ~−19, and its ¼-carrying transition is binary-vetoed); what
moved tonight is the SCREENING: three independent instruments (p-fit, hier
function ladder, Cassini floor) now point the same direction. Paper: §3
sharpened, §4.4 carries the decomposition, §6.4 the binary veto. Unification
check (does ν_p(0.65) fit the binaries too?) = Stage 5H.

## Stage 5H — the unification check: ν_p(0.65) is VIABLE on both systems (2026-07-23)

O9 executed ([calcs/stage5h_punified.py](calcs/stage5h_punified.py)): ν_p(0.65)
QUMOND-EFE table at g_N,ext=1.2a₀ (spherical identity 0.01%, ν(1)=1.423), v7
perspective-corrected 2D fit, seeds 31/101. **The binaries ACCEPT the
sharpened-screening occupation law: α̂ = 1.53 / 1.38 INTERIOR (no edge-riding
— the shape-rejection signature that killed boot is absent), Newton dead at
+92.8/+90.8, at a modest lnL concession to BE of 4.6/7.6** (comparable to the
per-seed realization scatter on the 4X footing; boot was +22.3/+17.1 WITH
edge-riding). Bookkeeping across systems, honestly stated side-by-side (the
two likelihoods are not summable): hierarchical SPARC+lensing prefers p065
over BE by **−56.4**; the binaries counter-lean toward BE by **+4.6/+7.6**.
ν_p(0.65) is thus the only function tested that is viable on BOTH systems —
BE loses the hierarchical galaxies by 56, boot loses the binaries by ~20 with
a shape rejection.

**Standing synthesis after the 5B–5H cascade: the data triangulate a
screening-sharpened occupation law — ½-branch-grade transition (binary-held),
p ≈ 0.58–0.65 Newtonian tail (hier-galaxy-held, §3-convergent,
Cassini-friendly direction), thermal deep structure with the ¼-vs-½ digit
open.** That is a sharper §2.4 construction spec than we had this morning.

Open flags carried forward (queued in TODO): (1) the binary amplitude under
p065 runs α̂ ≈ 1.4–1.5 interior — the a₀ translation under the sharpened tail
moves the binary-side high-tension direction naively worse, but the
galaxy-side a₀ moves too (hier p065 a₀ = 8.9e−11); the scorecard a₀ ladder
needs a self-consistent p-family recomputation before any tension claim.
(2) The exact p: hier profile 0.65, flat 0.578±0.12 — a joint-treatment p
measurement with the binary channel added is the natural next instrument.
(3) ν_p(0.65) quadrupole vs the 4K Cassini bound (sharper tail should
help; verify with the solver). (4) 6-seed budgets for 5F/5H (2 seeds
tonight, 4X footing).

## Stage 5I — the quadrupole is AMPLITUDE-LOCKED (2026-07-23)

O10a executed ([calcs/stage5i_quadrupole2.py](calcs/stage5i_quadrupole2.py),
4K multipole machinery verbatim; G1 regression on simple 0.1%, Newton
control exact, plateaus flat, resolution-doubling 0.1%). The hope: sharper
screening → smaller solar quadrupole → Cassini relief. The result kills the
hope cleanly: **ν_p(0.65) has q(1.2) = −0.0793 = 0.81× the ½-branch raw
quadrupole — but the binaries demand α̂ ≈ 1.45 under it, and 0.81 × 1.28 ≈
1.04: Q₂·α̂ = 3.9e-26 = 4.4× the Cassini cap, unchanged.** Same for the
geometric-mean bootstrap (raw 0.85×, α̂ pending 5K). **The Cassini tension
is amplitude-locked: the binaries pin the transition boost, the quadrupole
is transition-sourced, and the product is invariant under the tail-sharpness
dial.** No ν-function escape exists in this family direction; the 4K/4L
verdict stands — modified-gravity readings carry the 4× excess, the
EFE-respecting MI branch remains the Saturn-safe door. (In plain terms:
a null that strengthens the MI door rather than softening the paradox.)

## Stage 5J — the geometric-mean bootstrap LEADS the galaxies (2026-07-23)

O5 construction attempt, galaxy leg
([calcs/stage5j_gmboot.py](calcs/stage5j_gmboot.py)). The candidate from the
5B–5H spec: **ν_gm = 1 + n_BE(y^(3/4)√ν)** — mode frequency ω =
√(ω_source·ω_total), the geometric mean of the source-driven and
self-consistent prescriptions; equivalently the occupation argument is
x² = [T_U(g_N)/T_dS]·[T_U(g_obs)/T_dS]. Zero free parameters. Derived and
gate-verified series: **ν√y = 1 + x/3 + x²/12 — c₁ = ⅓ exactly (inside both
measured bands), c₂ = 1/12 (the occupation law's Bernoulli coefficient
SURVIVES the geometric mean)**; transition ν(1) = 1.433 (the 5H-accepted
grade); Newtonian tail e^(−y^(3/4)) (p-equivalent ¾, the top of the 5G hier
band). Numerics: Newton solver on w=√ν, monotone H, relative residual 1e-12
(first pass used an absolute-residual gate that false-FAILed on the deep
end's n~1e8 — machine precision misread; gate corrected to relative, logged).

Results: **flat M/L — gm concedes 2.0 to BE (−8395.72 vs −8397.72,
statistical silence; a₀ = 1.023e-10, f_ML = 1.31). Hierarchical converged —
gm = −10519.78, the NEW LEADER: −9.2 vs boot, −28.4 vs the fitted-p optimum,
−84.8 vs BE.** Honest sizing: gm-vs-boot (−9.2) is inside population
variance (the boot-BE bootstrap σ was ~65) — the robust statement is that
the sharp-tail class leads and gm is its parameter-free member sitting on
top; gm-vs-BE (−85) is the same grade as the 86%-bootstrap 5D flip. A
derived function with no dials matching the fitted-p family's optimum (and
nudging past it) is exactly what a correct construction should do. Binary
leg = 5K (running). Scout EXECUTED same day (Haiku sweep + primary read):
(a) the exact implicit form NOT FOUND; (c) NO published interpolating
function with deep NLO coefficient ⅓ (Milgrom's reviews + families searched);
(b) the one flagged candidate — Pazy & Argaman 1106.4108 ("a general
expression for the MOND interpolating function") — CLOSED BY DIRECT READ
(data/pazy_argaman_text.txt): their construction is FERMI–DIRAC (a₀ ↔ Fermi
energy), μ-side with total-acceleration argument, evaluated as a
dilog/numerical curve, and their own text states its deep-MOND corrections
are "exponentially small in 1/x, i.e., O(exp[−1/x])" — no polynomial NLO
at all, no implicit closure, no geometric mean. It cannot anticipate ν_gm.
Status: apparently unpublished (scout-grade negative, nearest candidate
primary-checked).

## Stage 5K/5L — the binary p-profile and the a₀ ladder: the tension gets two more axes (2026-07-23)

**5K (O10c + O5 binary, [calcs/stage5k_binprofile.py](calcs/stage5k_binprofile.py)):
four seeds × four functions**, all through the identical v7 machinery
(regression: p050 seeds 31/101 reproduce the 4X λ=1.00 values to the third
decimal; all six new EFE tables spherical-identity-gated at 0.01%). Per-seed
ΔlnL vs the occupation law (positive = worse): **p0578 +7.7±0.7, p065
+5.8±1.3, gm +10.4±2.8 (seed-31 outlier +18.6, others ~+7.6) — sign-consistent
12/12: the binaries measurably prefer the pure p=½ occupation transition.**
The sharper functions stay VIABLE (α̂ interior: p065 1.40±0.04 4/4, gm
1.44±0.08 4/4; p0578 mixed with one flat-top edge; no shape rejections) and
Newton is dead under every function (+77…+98). The 5H two-seed "realization
grade" phrasing is hereby FIRMED: at four seeds the binary counter-lean is a
real, modest, consistent preference — an order of magnitude smaller than the
hier-galaxy gains in the other direction (56–85).

**5L (O10b, [calcs/stage5l_ladder.py](calcs/stage5l_ladder.py)): the a₀
ladder under the sharp functions** (κ from the new e_N = 1.0/1.2/1.4 tables;
gates: BE-κ regression vs 4V exact at +0.916, α=1 identity, κ>0; κ(p065) =
+1.22, κ(gm) = +1.34). **Galaxy FLAT legs sit ON the horizon temperature:
gm a₀ = 1.02±0.09 (−0.2σ vs cH₀/2π), p065 1.08±0.09 (+0.5σ). Binary legs
under the sharp functions blow up: a₀ = 1.58/1.57 ± 0.11 = +4.9σ/+4.8σ**
(vs +1.9σ under the occupation law, 4V). The larger κ softens the exponent
but the α̂ ≈ 1.4 dominates.

**The synthesis after 5I+5K+5L: the two-system tension is now measured on
three axes and they cohere.** The hierarchical galaxies vote for sharper
tails (+56…+85); the binaries vote against (−6…−10), their amplitude under
sharp functions breaks the a₀–horizon lock (+4.9σ), and the quadrupole is
amplitude-locked regardless (5I). Meanwhile the FLAT galaxy treatment has no
tail preference (5J: BE best by 2) and sits exactly on cH₀/2π. The suspicious
alignment: everything anomalous lives in the hierarchical treatment, exactly
where a₀ runs low (0.84–0.89e-10) and f_ML runs high (1.5+) — and the 4W
lesson says the vertical (distance/inclination) channel is the great absorber
at SPARC depth. **Sober hypothesis: the hier tail vote is per-galaxy vertical
structure in disguise; if true, the coherent picture is the occupation law +
horizon a₀ (flat galaxies −0.2σ, binaries +1.9σ, binary function-preference
p=½).** Disambiguator = Stage 5M (function ladder WITH the measured-prior
vertical channel): collapse → sober reading wins; survive → the tension is
physical.

## Stage 5M — the vertical channel: HALF the tail vote absorbs, the a₀ anomaly RESOLVES (2026-07-23)

O11 executed ([calcs/stage5m_hierv.py](calcs/stage5m_hierv.py)): the function
ladder {BE, p065, gm, boot} re-run with the 4W vertical channel added to the
5D machinery — per-galaxy δv with MEASURED priors (σ_v from the SPARC table's
own distance & inclination errors, median 0.097 dex), closed-form profiled
alongside δd, converged. Gates: G1 σ_v→0 reproduces 5D's BE to 0.01;
nesting holds for all four; std(δv) = 0.106 ≈ the prior median (the channel
spends its budget, doesn't blow through it).

**Result 1 — the a₀/f_ML anomaly was vertical structure: with δv active,
every function's hierarchical a₀ returns to the horizon value** (BE 1.044,
gm 1.052, boot 1.073, p065 1.125 ×10⁻¹⁰ vs cH₀/2π = 1.042; f_ML back to
1.16–1.41; s_int → 0.044). The low-a₀/high-f_ML coupling that shadowed every
hier fit tonight is gone — the full hierarchy now AGREES with the flat
treatment and the binaries on the temperature. (Scorecard-grade: the
strongest a₀ = cH₀/2π statement the program has.)

**Result 2 — the tail vote HALVES but survives: Δ vs BE goes gm −84.8 →
−40.6, p065 −56.4 → −31.7. Result 3 — boot COLLAPSES: −75.6 → −8.6.** The
¼-cell's galaxy case was mostly vertical-degenerate; combined with its
binary veto (5F), the quantum self-consistent bootstrap is now out
everywhere. Only the ½-branch-transition sharp-tail functions (gm, p065)
keep real gains — exactly the two the binaries accept.

**End-of-night synthesis.** The geometric-mean bath ν = 1 + n_BE(y^¾√ν)
(c₁=⅓, c₂=1/12, parameter-free, apparently unpublished) leads the galaxy
ladder in EVERY treatment tried (hier −85, hier+vertical −41, flat tie) and
is binary-viable (α̂ interior 4/4). The occupation law (½) remains the
binaries' preferred function (sign-consistent +6…+10 over every sharper one)
and carries the mildest binary a₀ translation (+1.9σ vs +4.9σ). The residual
two-system tension: galaxies-with-full-hierarchy prefer gm by ~40 (point
estimate; a dv-ON galaxy bootstrap = O11b before quoting significance — the
5D-grade population σ was ~65); binaries prefer ½ by ~8 with the sharp-α̂
a₀ break as the harder objection. The ladder digit is now three-cornered —
½ (binary-anchored), ⅓ (geometric mean, best cross-system profile), ¼
(dead) — and the identifiability caveat of 4W applies to the surviving
galaxy-side preference (anchors = O1b remain the unlock). Quadrupole:
amplitude-locked at ~4× Cassini for ALL of them (5I); the MI door is the
program's standing escape, unchanged.

## Stage 5N/5P — THE MIXING DIAL: β measured; the why sharpened to one exponent (2026-07-23)

**5P theory ([calcs/stage5p_betafam.py](calcs/stage5p_betafam.py)): the three
bath cells are one family.** Mix the two frequency prescriptions
geometrically — ω = ω_source^(1−β)·ω_total^β, i.e. occupation argument
u = y^((1+β)/2)·ν^β — and derive (hand algebra, gate-verified to 2e-3 on the
coefficients and 1e-13 on the member functions):

  **c₁(β) = 1/(2(1+β));  c₂(β) = 1/(12(1+β)) + β/(8(1+β)²);
  p_tail = (1+β)/2;  hence c₁·p_tail = ¼ EXACTLY across the family.**

β=0 = occupation law (½, 1/12); β=½ = geometric mean (⅓, 1/12 — the
Bernoulli 1/12 RECURS at the symmetric point); β=1 = quantum bootstrap
(¼, 7/96). Structural facts: (i) deep scale invariance cannot choose β
(both frequencies coincide deep); (ii) **β=½ is the unique fixed point of
the source↔response exchange symmetry** (β→1−β) — symmetry forces c₁=⅓;
(iii) the geometric mean is where maximal coupling lands in canonical
settings: impedance matching (√(Z₁Z₂)), Curzon–Ahlborn maximum power
(√(T_hT_c)), and the normal mode with source-stiffness/response-inertia
(ω²=ω₁ω₂ — an MI-flavored reading). Stated as readings; the sharpened
construction problem = derive WHICH mix the horizon modes carry (C&T's
implicit choice was β=0). The "why ⅓" is now "why symmetric mixing" — one
exponent, measurable.

**And measured: galaxy β̂ = 0.64 (plain hier, interior, Δ+85.6 over β=0) /
0.45 (with the measured vertical channel — the most-controlled treatment —
interior, Δ+42.7; implied c₁ = 0.346 ≈ ⅓, p_tail = 0.72).** The
vertical-hardened galaxy measurement lands ON the exchange-symmetric point.
The binaries hold β ≈ 0 (5K/5O: the occupation law beats the β=½ and β=1
members sign-consistently). Endpoint regressions exact; the warm-chained
profile converged 2 lnL deeper than the 5M one-shot fits at b=0.5/1
(logged; direction benign). The two-system statement in final form: **the
degree of self-consistency of the bath frequency reads β ≈ ½ in galaxies
and β ≈ 0 in binaries** — a one-parameter, falsifiable disagreement.
Candidate physics for the split (noted, untested): the binaries sit in a
DOMINANT external field at the transition while the deep galaxy points are
quasi-isolated — if the response fraction β is configuration-dependent
(the bath sees the total field only where internal response dominates),
both readings can be right; alternatively the MI reading (trajectory
functional) naturally mimics β=0 for the binaries' statistic. Theory
target for O5. NNLO c₃(β) ladder: open.

**5O (O12, [calcs/stage5o_seeds.py](calcs/stage5o_seeds.py)): binary
six-seed completion.** Seeds 404/505 added for {p050, p065, gm}: the
counter-lean FIRMS — **occupation law preferred 18/18 across six seeds:
p065 +5.2 ± 0.9, gm +8.5 ± 2.1 (SE over seeds; gm's seed-31 +18.6 remains
the outlier, median ~+7)**; α̂ six-seed: p050 = 1.078 ± 0.023 (softens the
4V binary-a₀ pull to ≈ +1.8σ), p065 = 1.41 ± 0.04, gm = 1.36 ± 0.07 — all
interior 18/18; Newton dead 18/18 (+77…+98). Six-seed a₀ translations:
p065 → 1.59 ± 0.11 (+5.1σ), gm → 1.51 ± 0.11 (+4.3σ) — the sharp-function
temperature break stands at six seeds.

**5N (O11b, [calcs/stage5n_dvboot.py](calcs/stage5n_dvboot.py)): the
surviving lead bootstrapped.** Full-fit regression reproduces 5M exactly
(BE −12152.49, gm −12193.04). 40 paired galaxy reps (warm-started lite):
**gm − BE = −29.3 ± 53.0, gm better in 29/40 (72.5%), percentiles
[−72.8, −21.0, +21.4].** The vertical-robust gm lead is a LEAN, not a
detection — quote it as such everywhere (the 5D pre-vertical flip was 86%;
the vertical channel took half the effect and half the significance).

## Stage 5Q (2026-07-23): the c₃(β) NNLO ladder — the family's algebra one rung deeper (O13c)

[calcs/stage5q_c3ladder.py](calcs/stage5q_c3ladder.py) → data/stage5q_c3ladder.txt.
Writing x = √y and S = xν, the family's implicit equation is exactly
S = x/2 + S^(−β) + x²S^β/12 − x⁴S^(3β)/720 + O(x⁶). Symbolic series to x⁴
(sympy; G1 reproduces the 5P c₁, c₂ formulas EXACTLY; G2 mpmath 50-digit
roots, residual scaling exponent 5.00/4.99):

- **c₃(β) = β(3β+1)/(24(1+β)³)** — its ONLY zero is β = 0: the Bernoulli
  zero of the pure occupation law (the x² term of 1/(1−e^(−x)) vanishes
  because odd Bernoulli numbers beyond B₁ vanish) is UNIQUE in the family.
  Any admixture of self-consistency (β > 0) switches the NNLO rung on.
- c₄(β) = (502β³+51β²−54β−8)/(5760(1+β)⁴); c₄(0) = −1/720 = −B₄/4! ✓ (the
  exact BE check). Member values: β=¼ → c₃ = 7/750; β=½ → 5/324; β=¾ →
  13/686; β=1 → 1/48.
- **In the p_tail-rescaled variable the ladder goes polynomial**: c₁p = ¼
  (5P), c₂p² = 1/48 + 5β/96 (linear!), c₃p³ = β(3β+1)/192. The (1+β)
  denominators are an artifact of the x-variable; the natural-variable
  ladder is polynomial in β.
- **Bonus structure: the third log-cumulant d₃ = c₃ − c₁c₂ + c₁³/3 =
  β(2β−1)/(16(1+β)³) vanishes at β = ½** (and β=0) — the geometric-mean
  member is the unique β > 0 point where the deep expansion of ln(xν) has
  no cubic term. Another independent algebraic distinction of the
  symmetric point, logged as an observation (no derivation claimed).

Observational size: Δc₃(0→½) = 0.0154 ⇒ Δν/ν ~ 1.5e−5 at x = 0.1 —
structural algebra, not an instrument. The measurable discriminators stay
c₁, p_tail, ν(1).

**Verdict: SUCCESS (algebra complete, gates exact; two new structural
facts — BE's Bernoulli zero is unique in the family, and β=½ kills the
third log-cumulant).**

## Stage 5T (2026-07-23): configuration-dependent β — the galaxy vote decomposed by depth (O13b)

[calcs/stage5t_betasplit.py](calcs/stage5t_betasplit.py) → data/stage5t_betasplit.txt.
Three instruments on the plain-hier machinery. G1: the full-sample β
profile reproduces 5P to d = ±0.00 at every node (warm-chained, exact).
G2: the decomposition rows sum to the profile totals to 1e−6.

**B. The decomposition** (Δ(−2lnL) vs β=0 by fixed y-bin, evaluated at
each β's joint best fit — pure evaluation, zero convergence risk):
**the galaxies' β̂ ≈ 0.5–0.64 is a COMPROMISE between opposing arms.**
At β=½: ultra-deep y<0.03 (56 pts) votes AGAINST +9.3 (the deep zero
point wants the occupation c₁=½; the penalty grows monotonically to
+13.2 at β=1); 0.03–0.1 mild FOR (−4.3); 0.1–0.3 strong FOR (−20.0);
transition 0.3–1 mild FOR (−5.5) but flips AGAINST at β=1 (+6.9 — the
transition defends its 1.58 boost); the Newtonian side carries the vote
(y 1–3: −24.4; y>3: −36.8 — the screening tail again, convergent with
the 5G tail decomposition); lensing −1.0; M/L prior mild FOR (−2.2:
sharper functions need less per-galaxy M/L bending).

**C. Free arm fits**: HIGH arm alone (y_fid≥1, 627 pts): **β̂ = 0.76
INTERIOR, Δ +5.8** (la0 edge-rides at β=0 only — BE's soft tail is
unfittable there without pushing a₀ into the bound; interior at all
β≥0.25). LOW arm (y_fid<1, 2073 pts + lensing): monotone to the β=1
edge (−65.8) but **RIDGE-COMPROMISED — f inflates 1.57→1.94 across the
grid** (the deep-arm a₀·f degeneracy with no Newtonian anchor); flagged
and NOT interpreted. The decomposition is the trustworthy localizer.

**The finding: β is not one number for the galaxies — the effective
mixing exponent READS β→0 in the deep limit and β ≈ ½–¾ in the
screening tail, within the single SPARC dataset.** Cross-system, the
pattern is now uniform: every regime that measures the TRANSITION or
the DEEP ZERO POINT votes β ≈ 0 (binaries at the transition, 5K/5O;
galaxy ultra-deep; galaxy transition bin at high β), and the single
vote for β > 0 is the Newtonian screening tail. Structurally: the
family locks c₁·p_tail = ¼, so no member can deliver what the joint
data ask — c₁ ≈ 0.4–0.5 (deep + binaries, 4S/4X) AND p ≈ 0.65–0.75
(tail, 5G/5T) — a product of ~0.30–0.35. **The data strain against the
family's own lock; the construction target for O5 is now "break the ¼
lock upward": BE-grade deep zero point + sharpened screening,
decoupled.** (Caution registered: the ultra-deep objection is 56 points
and ~2σ-grade in bin terms; nuisance-locked evaluation.)

**Verdict: SUCCESS as an instrument — and the answer points at
DIFFERENT PHYSICS: β runs with regime rather than being one number.
The reconciliation hypothesis graduates from speculation to a measured
pattern, and the microphysics question sharpens from "which β" to "why
does the response admixture vanish deep and grow through the
screening tail".**

## Stage 5R (2026-07-23): the binary β-profile — β pinned at zero, 24/24 (O13a)

[calcs/stage5r_binbeta.py](calcs/stage5r_binbeta.py) → data/stage5r_summary.txt,
data/stage5r_profile.txt. Dedicated ν_β tables at β = 0.25/0.75
(e_N = 1.2a₀); G1 end-to-end member regression — the β=0 table through
the solver reproduces the stored BE table to 1.3e-15; G2 spherical
identities 0.01%; G3 the EFE boost is strictly monotone DECREASING in β
at every y (B(y=1) = 1.354/1.276/1.220/1.179/1.147 across the dial).
Boot completed on seeds 202–505. Full 6-seed × 5-β lnL matrix on one
machinery footing (4X/5K patch set, exact-count asserts).

- **β = 0 preferred over every β > 0 on every seed: 24/24
  sign-consistent.** Mean ΔlnL behind β=0: −4.14 ± 1.11 (β=¼), −8.50 ±
  2.13 (½), −7.90 ± 1.85 (¾), −15.35 ± 1.99 (1) (SE over 6 seeds).
- **Mean-profile crossings: β < 0.030 (1σ, ΔlnL=0.5), β < 0.121 (2σ,
  ΔlnL=2.0).** The binaries' β = 0 is SHARP, not shallow — O13a's
  question answered.
- Shape-rejection texture: b075 edge-rides α̂→2.0 in 5/6 seeds (boot
  3/6) — high-β members cannot buy their weakened transition back with
  amplitude. b025 stays interior 5/6 at α̂ = 1.57 ± 0.11 (interior-only
  mean 1.48; per-seed spread 1.16–2.0 — the α̂ flat-top scatter).
- Two-system table + an indicative additive joint (vertical-hardened
  galaxies + binaries in −2lnL): joint minimum at β ≈ 0.25 (−29.8) —
  carried with the 5N caveat (galaxy bootstrap yardstick ±~50 on this
  contrast) and SUPERSEDED by 5T: the galaxy Δ is itself a tail/deep
  compromise, so the additive joint double-counts a structural tension.
  Do not lean on it; the regime decomposition is the honest frame.

**Verdict: SUCCESS — the binary mixing dial reads β = 0 sharply (1σ
upper bound 0.03). Combined with 5T: the two-system split is a REGIME
split. Every probe of the transition and the deep zero point — binaries,
galaxy ultra-deep, galaxy transition bin — holds β ≈ 0; the galaxy
screening tail is the lone β > 0 voter in the program.**

## Stage 5S (2026-07-23): the β-family quadrupole scan — Saturn's veto is β-blind (O13d)

[calcs/stage5s_betaquad.py](calcs/stage5s_betaquad.py) →
data/stage5s_betaquad.txt (solver cache data/stage5s_q.npy). 4K/5I
multipole machinery verbatim; G2 Newton control 0; G1 simple regression
0.1% vs 4K; member gate β=½ reproduces the 5I gm q to 0.00%; G4 hi-res
≤ 0.1% all members. q(e_N=1.2) runs −0.0988 → −0.0901 → −0.0833 →
−0.0777 → −0.0731 across β = 0→1 (monotone; transition-sourced as
always). Lock join with the 5R/5O six-seed amplitudes:

    β      Q₂(α̂)      × Cassini   status
    0.00   3.64e-26    4.0×        interior 6/6 (measurement)
    0.25   4.83e-26    5.4×        interior 5/6
    0.50   3.88e-26    4.3×        interior 6/6 (measurement)
    0.75   5.19e-26    5.8×        interior 1/6 → LOWER BOUND
    1.00   4.61e-26    5.1×        interior 3/6 → lower-bound-grade

**Every member sits ≥ 4× the Cassini cap; the edge-ridden members are
lower bounds (α capped at 2.0). The mixing dial provides no
interpolating-function escape — the 4K/5I amplitude-lock is now
family-complete.** The ±20% wobble in Q₂·α̂ across members is α̂
realization scatter (flat profile tops), not structure. The MI and
EFE-screened doors (§8.2/#8) are unchanged — and strengthened, since
the last function-space escape hatch inside this family is now shut.

**Verdict: SUCCESS (the expected no-escape verified at family level;
Saturn's veto is β-blind).**

## Stage 5U (2026-07-23): the spontaneous-fraction bath — the running β DERIVED (O5 construction)

[calcs/stage5u_runbeta.py](calcs/stage5u_runbeta.py) → data/stage5u_runbeta.txt.
The 5T pattern (β = 0 where occupation is high, β → ½ where it is low)
gets a mechanism. The mixing weight β is the RESPONSE share of the mode
frequency; in a driven thermal mode the channel split is
Einstein-coefficient physics — stimulated processes (rate ∝ n) follow
the DRIVE (the source-prepared field), spontaneous emission (rate ∝ 1)
probes the mode's OWN structure. Weighting the response by the channel
share gives two zero-parameter candidates:

- **F1 (rate share): β = ½·1/(1+n) = 1/(2ν)**
- F2 (zero-point energy share): β = ½·(½)/(n+½) = 1/(2(2ν−1))

Both: classical limit n≫1 → β→0 (source-driven — Cadoni–Tuveri's
implicit choice becomes the CLASSICAL limit of the bath); quantum limit
n≪1 → β→½ (the exchange-symmetric point, now WITH a reason: it is where
only spontaneous processes remain). The ½ asymptote is inherited from
5P's exchange symmetry; the running is the new derived content. The
weighting CHOICE is a stated reading, not a proof — but the form's
consequences are exact and were written down before any fit:

- **c₁ = ½ EXACTLY for both** (the running dies fast enough deep that
  the occupation zero point survives untouched) — the 5T ultra-deep vote.
- **tail p = ¾ EXACTLY** (the β∞ = ½ asymptote) — the sharpened
  screening the galaxy tail votes for, at the symmetric-point value.
- **lock product c₁·p = ⅜ = 0.375** — the constant-β family's ¼ lock
  broken UPWARD into the band the data ask for (~0.3–0.375).
- c₂ goes negative: −1/6 (F1), −1/24 (F2) vs BE's +1/12 — a rung-2
  signature, unmeasurable at the current 0.1σ reach (disclosed);
  c₃ = 13/48 (F1), 1/24 (F2). ν(1) = 1.4702 (F1) / 1.4943 (F2) — between
  gm (1.433) and BE (1.582): the binaries decide (5W).

Gates: G1 sympy series exact (c₁ = ½ symbolic identity both); G2 solver
residual 8.4e-13/9.5e-13 relative over y ∈ [1e-8, 1e4] + root-uniqueness
scans PASS; G3 mpmath 50-digit series residual scaling 3.98/3.99
(expect 4); G4 numeric tail exponent 0.7500/0.7501. Priority scout
(Haiku, 4 targeted sweeps): running/occupation-dependent interpolating
exponent NOT FOUND; Einstein-coefficient (spontaneous/stimulated)
derivation of an interpolating function NOT FOUND; the (c₁=½, p=¾)
combination NOT FOUND; no Cadoni–Tuveri follow-up modifying the
occupation argument — all scout-level negatives, quoted as "apparently
unpublished" pending the pre-circulation primary-source pass.

**Verdict: SUCCESS (construction stage) — the O5 question "why would β
run?" has a derived, zero-parameter, exactly-gated candidate answer:
because the response enters through the spontaneous channel. Whether it
is TRUE is 5V/5W's job — c₁ = ½ and p = ¾ were fixed before the fits.**

## Stages 5V/5W/5X (2026-07-23): the derived function tested — galaxies crown it, the transition still tithes, Saturn unmoved (O14)

[calcs/stage5v_rungal.py](calcs/stage5v_rungal.py),
[calcs/stage5w_runbin.py](calcs/stage5w_runbin.py),
[calcs/stage5x_runquad.py](calcs/stage5x_runquad.py) →
data/stage5v_rungal.txt, stage5w_verdict.txt, stage5x_runquad.txt.

**5V galaxies (BE regression gates d = −0.00 both treatments):**

    treatment    F1        F2        [gm comparator]
    plain hier   −101.29   −107.22   [−84.79]   (Δ vs BE)
    vertical     −44.31    −51.89    [−42.69]

**F2 is the new outright leader of the galaxy ladder under BOTH
treatments** (beats gm by 22.4 plain / 9.2 vertical); F1 second. The
5T-decomposition prediction — keep the tail's gain, recover the
ultra-deep bins — landed. Vertical-hardened a₀ = 0.98–1.0e-10
(horizon-adjacent), f_ML = 1.32–1.34, s_int unchanged.

**5W binaries (six seeds, same machinery as 5K/5O/5R; all gates pass):**

    law   mean ΔlnL vs p050   worse   interior   α̂
    gm    −8.50 ± 2.13        6/6     6/6        1.358 ± 0.068
    rb1   −5.77 ± 1.93        5/6     6/6        1.395 ± 0.039
    rb2   −5.47 ± 2.27        5/6     6/6        1.487 ± 0.017

The pre-stated bands (−3 = accepted, −6 = gm-grade rejection) are
STRADDLED — the honest middle verdict: **the running functions close
about a third of the transition gap (no edge-riding anywhere, seed 303
flips positive) but the binaries still lean to the pure occupation law
by ~2.4 SE.** The residual is localized exactly where 5T pointed: the
transition wants even less admixture than β(ν(1)) ≈ 0.29–0.34 — the
weighting's refinement target is a faster-than-1/ν die-off through the
transition.

**Open ledger item (flagged, not computed): the a₀ translation.** κ
tables (e_N = 1.0/1.4) were not built for F1/F2. At gm-grade κ ≈ 1.34,
α̂ = 1.49 would translate to a₀ ≈ 1.6e-10 (+4–5σ — the sharp-function
temperature strain repeats); at BE-grade κ it would not. Queued (O15a);
until then the temperature lock's binary row still votes pure-BE.

**5X Saturn:** q(1.2) = −0.0903 (F1) / −0.0943 (F2), hi-res 0.0%;
with the 5W amplitudes: **Q₂ = 4.31e-26 (4.8×) / 4.80e-26 (5.3×
Cassini)** — inside the family band (4.0–5.8×), the amplitude lock
holds. No escape; the MI door unchanged.

Indicative additive two-system ledger (−2lnL vs BE, vertical-hardened
galaxies + binaries; 5N caveat applies): **F2 nets −41** (gm: −25.7) —
the best cross-system function the program has, and it was DERIVED:
c₁ = ½ and p = ¾ were on paper before any fit touched data.

**Verdict: SUCCESS with one honest reservation. The derivation held its
pre-registered predictions and took the galaxy ladder outright; the
binaries upgraded it from gm-grade to two-thirds-accepted; Saturn is
unmoved (as for every MG reading). NEEDS REFINEMENT on exactly one
localized number — the transition admixture — plus the pending κ/a₀
leg. The spontaneous-fraction bath is now the program's leading
candidate function.**

## Stages 5Y/5Z/6A/6B (2026-07-23 late): the two-leg refinement — galaxies converge, the binaries veto through the EFE, the a₀ wall (O15a/b)

[calcs/stage5z_twoleg.py](calcs/stage5z_twoleg.py),
[calcs/stage6a_twogal.py](calcs/stage6a_twogal.py),
[calcs/stage6b_twobin.py](calcs/stage6b_twobin.py),
[calcs/stage5y_kladder.py](calcs/stage5y_kladder.py).

**5Z (algebra, all gates exact):** two-leg spontaneity — the
self-frequency requires the round trip (emission AND reabsorption
spontaneous), squaring the share: F3: β = 1/(2ν²), F4: β = 1/(2(2ν−1)²).
**Both Bernoulli rungs now survive exactly (c₁ = ½ AND c₂ = 1/12)**;
c₃ = −¼ / −1/16; tail p = ¾ preserved; ν(1) = 1.5032/1.5373 (raised
toward BE's 1.582). Verdict: SUCCESS.

**6A (galaxies):** the ladder improves AGAIN — plain: F3 −111.4,
F4 −108.7 (F2 was −107.2, gm −84.8); vertical-hardened: F3 −50.9,
**F4 −64.2 — the largest controlled-treatment lead any function has
held** (F2 −51.9, gm −42.7); a₀ = 0.97e-10 horizon-adjacent. The galaxy
ladder is now MONOTONE toward "occupation-deep + sharp-tail" in both
treatments. Verdict: SUCCESS.

**6B (binaries): THE INFORMATIVE VETO.** rb3 −6.84 ± 1.58 (interior,
α̂ = 1.500); **rb4 −6.96 ± 2.08 with α̂ EDGE-RIDING 2.0 on ALL SIX
seeds** — despite rb4 having the closest-to-BE transition of any
sharpened function. The two-leg refinement did NOT close the 5W
residual; it worsened it. Across all eight sharpened functions now
tested (gm, p0578, p065, F1–F4), the binaries pay −5…−8.5 nearly
INDEPENDENT of ν(1): **what they reject is not the transition value but
the function's behavior in the screening region UNDER THE DOMINANT
EXTERNAL FIELD** (the EFE evaluates the function at argument
~e_N + y ≈ 1.2–1.5; the fits at e_N = 1.2a₀ read the early-screening
shape, and they want BE's SOFT p = ½ screening there — while the galaxy
Newtonian arm at the same y-range wants p ≈ ¾). Verdict: DIFFERENT
PHYSICS — no local ν(y) reshaping can satisfy both systems, because
they probe the same region under different ambient fields.

**5Y (the a₀ wall):** κ = 1.35/1.33/1.34/1.26 for F1–F4 (all gm-grade;
G3 BE regression +0.916 exact). Binary a₀ translations: **F1 1.54±0.10
(+5.2σ), F2 1.62±0.09 (+6.3σ), F3 1.63±0.09 (+6.2σ), F4 formally
+11.4σ but α̂ = 2.0 edge → n/a as measurement.** Only pure BE passes
(+1.9σ). The temperature lock is the binaries' second, independent veto
of every sharpened function under the external field. Verdict: SUCCESS
as instrument; the wall is real.

**Synthesis — the tension's final form and its one open door:** the
hierarchical galaxies reward sharp screening monotonically (best-ever
−64 with a derived function); the binaries, probing the SAME functions
through the dominant external field, reject every one (−5…−8.5 lnL AND
+4–6σ in temperature). The residual physics is IN THE EFE REGION.
Reading with a mechanism (O13's configuration-dependence, now
concrete): **the drive is whichever field dominates** — a mode prepared
by the ambient Galactic field is externally driven (β → 0, BE-like)
regardless of local occupation; a self-sourced galaxy point keeps the
spontaneous-share admixture. Simplest realization: β = w_self ·
(two-leg share), w_self = g_int/(g_int+g_ext). Isolated limit → F4
exactly (all 6A results carry over); external-dominated limit → pure
BE exactly. This abandons "ν is a local function of |g_N| alone" —
NON-LOCALITY, converging with the 4K MI door and the 4L MI-tie from
independent directions. Executed as Stage 6D.

## Stage 6C (2026-07-23 late): the F4 galaxy lead bootstrapped — the strongest grade yet (O15c)

[calcs/stage6c_f4boot.py](calcs/stage6c_f4boot.py) → data/stage6c_f4boot.txt.
40 paired galaxy-resample reps (5N machinery verbatim; full-fit
regression: BE −12152.49 exact, F4 −12215.07 vs 6A's −12216.70, d PASS):
**F4 − BE = −57.4 ± 38.3, F4 better in 37/40 (92.5%), percentiles
[−89.4, −60.5, −16.9].** Comparator (5N): gm − BE = −29.3 ± 53.0,
29/40 (72.5%). The derived two-leg function's vertical-hardened lead
survives the galaxy-population bootstrap at ~1.5σ with 92.5%
sign-consistency — the strongest bootstrap grade any function lead has
achieved in this program; still a LEAN by the 5N yardstick (not a
detection), and quoted as such. Verdict: SUCCESS.

## Stage 6D (2026-07-23 late): the drive-weighted bath — EXCLUDED by shape, passed by temperature (O16 first pass)

[calcs/stage6d_drivebath.py](calcs/stage6d_drivebath.py) →
data/stage6d_verdict.txt. The simplest configuration rule — β = w_self ·
(two-leg share), w_self = g_int/(g_int+g_ext) pointwise — implemented at
solver level (one changed line; barycenter response = BE(e_N) exactly at
w_self = 0). Gates: G1 isolated identity → F4 at 0.01%; G2 BE/F4
sandwich PASS; G3 wide-pair limit → BE at 0.35% (the deep pairs are
essentially pure occupation law, as designed).

**Binary six-seed: mean −9.64 ± 1.49 vs p050, better in 0/6 — WORSE
than F4 (−7.0) and gm (−8.5), despite BE-grade interior amplitudes
(α̂ = 1.098 ± 0.021, 6/6).** The likelihood punishes exactly the rule's
predicted signature: the mid-separation sag (B(1) = 1.324 vs BE 1.354;
B(3) = 1.143 vs 1.201) where w_self ≈ 0.5–0.7 pulls toward sharp
screening. **The binaries want the soft p = ½ screening at ALL their
separations, including where their internal field dominates the
ambient — the pointwise-field-ratio configuration rule is REFUTED.**
The two binary vetoes now split cleanly: κ(dwf) = 1.033 → a₀ = 1.31 ±
0.13 (+2.1σ) — dwf PASSES the temperature row while failing the shape
row; amplitude and shape are independent verdicts.

Surviving readings, both flagged as formed AFTER this exclusion:
(i) ambient-GATED mixing — the system's environment as a whole, not the
local ratio, sets the bath character (binary modes sit in the Galaxy's
e_N ≈ 1.2a₀ ambient at every separation → n_amb ≈ 0.5, quantum-ish →
source-driven; galaxy disk points sit in e ~ 0.02a₀ → n_amb ~ 7,
classical → admixture allowed). At the two anchors this reduces to
"binaries = BE, galaxies = F4" BY CONSTRUCTION — no in-sample
discrimination beyond them except the Chae high-e_N galaxy leg (power
≲ 2–5% per 5B/5E); parked as O16-remaining with that status stated.
(ii) The modified-inertia reading: a trajectory functional can carry
different effective screening for eccentric bound orbits in a strong
ambient field than for circular rotation at the same local field — the
cross-system function split is exactly where MI would show, and the 4L
tie keeps it open (while 5T's INTRA-galaxy tail-vs-deep split still
needs the function shape, MI cannot produce that part).

**Verdict: DIFFERENT PHYSICS, confirmed by exclusion — the two-system
function split is real, sharp (shape −5…−10 lnL AND temperature +4…+6σ
against every sharpened function on the binary side; −64 for the
sharpest derived function on the galaxy side), and NOT reconcilable by
any local ν(y) or by the pointwise drive-weighted rule. What remains:
system-level configuration dependence (ambient-gated, weakly testable
in-sample) or trajectory dependence (MI, the standing door).**

## Stage 6E (2026-07-23 night): the ambient-gated bath DERIVED — admixture = local quantumness × ambient classicality (O5)

[calcs/stage6e_ambgate.py](calcs/stage6e_ambgate.py) → data/stage6e_ambgate.txt.
The theory question attacked head-on, and the first guess FAILED ON
SIGN — instructively: "ambient traffic dilutes the spontaneous share"
predicts suppression for WEAK ambients, but weak fields mean deep,
highly occupied (classical) ambient baths (n_amb = n_BE(√(e_N/a₀)):
galaxies e~0.02 → n_amb ≈ 6.6; binaries e≈1.15 → n_amb ≈ 0.52). The
measured pattern is the REVERSE: sparse ambient → no admixture. The
reading that carries the right sign, same Einstein-coefficient grammar
as 5U/5Z: **a mode can only dress itself to its self-consistent
frequency if the ambient reservoir can assist the exchange — and
reservoir assistance is a STIMULATED process (∝ n_amb/(1+n_amb)). A
classical ambient is a floppy medium (mode re-equilibrates); a quantum
ambient is stiff (mode pinned at the drive frequency).** The unified
zero-parameter rule:

    β = ½ · [local zero-point share]² · [ambient stimulated share]²
      = g_amb · ½/(2ν−1)²,   g_amb = [n_amb/(1+n_amb)]²

Admixture = (local quantumness) × (ambient classicality). System-level
by construction (one n_amb per system — explains the 6D exclusion
pattern with no mid-separation sag). POST-HOC STATUS FLAGGED
throughout: built after 6B/6D; pre-registered content = the exact
algebra + two p-postdictions + the 6F/6G bars. Exact consequences
(all gates pass; sympy at SYMBOLIC g): c₁ = ½ and c₂ = 1/12
g-INDEPENDENT; **c₃ = −g/16 exactly — the entire ambient dependence
enters the deep series at one rung, linearly in the gate**; tail
p = ½ + g/4: **galaxies g = 0.754 → p = 0.688 (measured band
0.65–0.75 ✓), binaries g = 0.117 → p = 0.529 (their held ½ ✓)** — both
measured screening exponents post-dicted from zero parameters; ν(1) =
1.548/1.577; g = 0 → BE at 1e-12, g = 1 ≡ F4. Priority scout (4
sweeps): environment-dependent interpolating-function SHAPE, decoherence
gating, Einstein-coefficient environmental modulation, and two-function
MOND all NOT FOUND (nearest art = Pazy 1106.4108's holographic
quantum-classical crossover, closed by direct read in 5J).

**Verdict: SUCCESS (construction) — the O5 question "why would the
ambient gate at system level" has a derived answer with the right sign:
because self-consistent dressing needs stimulated reservoir assistance,
and only classical (deep, occupied) ambients provide it.**

## Stages 6F/6G (2026-07-23 night): the ambient gate tested — ONE FUNCTION PASSES BOTH SYSTEMS (O16)

[calcs/stage6f_ambgal.py](calcs/stage6f_ambgal.py),
[calcs/stage6g_ambbin.py](calcs/stage6g_ambbin.py) →
data/stage6f_ambgal.txt, stage6g_verdict.txt.

**6F galaxies (the falsifiable leg; BE regressions d = −0.00 both):**
vertical-hardened AMB(g=0.754) = **−59.05 vs BE — PASSES the
pre-registered ≤−50 bar**; second-best ever (F4's undiluted −64.2 above
it, F3/F2 below); dilution cost only 5.2 vs F4. Plain: −92.40 —
**MISSES the pre-registered ≤−100 bar by ~8** (the ambient-band ends:
e=0.01 → −97.4 grazes it; e=0.05 → −82.6); disclosed as a partial. The
sensitivity ordering (more isolated → sharper → better) is the
in-sample Chae dial.

**6G binaries (pre-registered: within ~−2 of p050): mean −0.88 ± 2.66,
better in 2/6, interior 6/6, α̂ = 1.060 ± 0.024 — ACCEPTED** (per-seed
−0.3/−0.8/−12.5/+4.7/+5.8/−2.2; seed 202 the lone objector; the mean
is noise-grade). κ = 0.924 (BE-grade, as predicted). **a₀ = 1.28 ±
0.15 → +1.6σ — the best binary temperature row in the program** (pure
BE itself: +1.9σ). Both binary vetoes — shape AND temperature — passed
by the same function that holds −59 on the galaxies.

**The joint ledger (vertical galaxies + binaries, additive −2lnL vs
BE): AMB nets −57.3** — best cross-system function ever measured here
(F4: −50.3 and binary-vetoed; BE: 0 and galaxy-inferior; gm: −25.7).
**For the first time in the program, ONE zero-parameter function
passes both systems.** Carried caveats: post-hoc construction (the
binary pass was near-guaranteed by design; the NON-guaranteed content
that landed: the vertical galaxy survival, the κ/a₀ = +1.6σ, the two
p-postdictions, the plain-bar shortfall honestly against); the
out-of-sample tests are the Chae per-galaxy softening dial and DR4
weak-ambient binaries (predicted to sharpen toward p = 0.69).

**Verdict: SUCCESS — with the post-hoc flag carried openly. The
two-system tension that 5Y/6B/6D sharpened into "different baths" is
reconciled by the derived ambient gate at every in-sample number, and
the reconciliation makes falsifiable out-of-sample predictions with a
quantitative dial. The program's function, as of tonight: ν = 1 +
n_BE(y^((1+β)/2)ν^β), β = [n_amb/(1+n_amb)]²·½/(2ν−1)².**

## Stage 6H — the grammar formalized: β = ½·[q_loc·s_amb]^L, the leg count as a measurable (2026-07-23, O5 second pass)

User asked two things: keep deriving ("that's the only way we might find
some other possible predictors... I really wanna see something truly
original"), and define "local quantumness" precisely ("to me that would
almost seem random and might not fit the simulations properly"). The
answer to the second: it is NOT stochastic — it is the zero-point share
of the mode's energy, ½/(n_loc+½) = 1/(2ν−1), a deterministic function
of the local field through the self-consistent occupation n_loc = ν−1;
it is already inside every amb table via the implicit Newton solve
(6E gates, residual ~1e-12). The genuinely prescriptive joint is the
AMBIENT split (one n_amb per system) — which this stage + 6I target.

[calcs/stage6h_grammar.py](calcs/stage6h_grammar.py) →
data/stage6h_grammar.txt. All gates PASS first run.

**G1 — branch selection is now derivation-grade.** The dispersive
frequency pull of a system coupled to a mode of occupation n is
λ²(2n+1)/Δ EXACTLY (JC manifold block, sympy) — vacuum 1 + thermal 2n,
an energy-type (n+½) weight, not a rate-type (n+1) weight. A
frequency-dressing process therefore carries the zero-point share
1/(2n+1) = 1/(2ν−1): this SELECTS the F2/F4 lineage over F1/F3
(a-posteriori explanation of the 6A/6B pattern) and fixes the local
factor of the 6E rule. Still reading-grade: that horizon dressing is
dispersive.

**G2 — THE BERNOULLI-BREAK RUNG IS THE LEG COUNT** (the stage's
discovery; symbolic-s ladders, c₄(s→0) = −1/720 = 5Q's BE value all L):

    L=1: c₂ = 1/12 − s/8 (Bernoulli c₂ KILLED), c₃ = s(3s+1)/96
    L=2: c₂ = 1/12 survives, c₃ = −s²/16 (= 6E), c₄ = s²/192 − 1/720  ← NEW rung
    L=3: c₂ AND c₃ = 0 survive, c₄ = −s³/32 − 1/720

The first deep coefficient the admixture touches is c_{L+1} — the
grammar's structural exponent is readable off the deep ladder, so the
5T instrument (the ultra-deep arm's Bernoulli vote) can MEASURE L.
Note the new L=2 rung flips c₄'s sign at galaxy gates (+0.0025 at
g=0.754). Tails: p = ½ + s^L/4 (numerics 0.717/0.689/0.665 at the
fiducial gate, all solver gates 1e-12, uniqueness scanned); ν(1) =
1.505/1.548/1.568. Binary-side tails 0.586/0.529/0.510 — the binary leg
is NOT the L discriminator (L=1's 0.586 sits between the 5K-tested ½
and 0.65); the deep rung is.

**G5 — the measured ambients (Chae+21 Table 3, 109 galaxies):**
maxclust median e_N = 0.0050 → g = 0.868, p = 0.717 (span 0.691–0.729);
noclust median e_N = 0.0006 → g = 0.952, p = 0.738. Both SHARPER than
the 6E/6F fiducial e = 0.02 (g = 0.754, p = 0.688) — the fiducial was
conservative. Honest note: the measured-ambient tail postdiction moves
to p ≈ 0.71–0.74, still inside the 5G/5T band 0.65–0.75 but nearer its
upper edge. (Value coincidence flagged to prevent confusion: median
per-galaxy g = 0.868 numerically equals the fiducial s = 0.8681.)

**Pre-registered 6I bars (on disk before any 6I fit):** A — grammar
survives iff L=2 ranks best of {L1, L2, L3} vertical (an L=1 win =
tail out-votes deep = grammar strike). B — per-galaxy measured gates
must not hurt (≥ global both treatments); plain reaching ≤−100
RESOLVES the 6F disclosed partial. C — quadrupole record expected ~4×
Cassini, no rescue.

**Verdict: SUCCESS (construction + two new exact structures: the
Bernoulli-break↔leg-count correspondence and the c₄ rung; the JC
(2n+1) anchor upgrades branch selection from reading to derivation).
The falsifiable content is deferred to 6I by construction.**

## Stage 6I — the leg count measured (L = 2) + the measured ambients resolve the plain bar (2026-07-23, O16a/O16c + the 6H contest)

[calcs/stage6i_chaegate.py](calcs/stage6i_chaegate.py) →
data/stage6i_chaegate.txt (+ cache data/stage6i_q.npy). 153 galaxies,
2700 points + 12 lensing; Chae matched 91/149 (unmatched + lensing
carry the matched median — disclosed). Regressions EXACT: BE both
treatments d = −0.00 vs 5P; AMBg deltas dd = +0.00 vs 6F (−59.05
vertical / −92.40 plain) — the unified per-point-gate solver path
reproduces 6F to the hundredth.

**Part A — the leg-count contest (vertical-hardened, fiducial gate
s = 0.8681, pre-registered in 6H): L = 2 WINS.**

    L1 −52.76 | L2 −59.05 | L3 −54.52   (Δ vs BE; lower = better)

Both flanking integers lose — a DISCRETE INTERIOR measurement of the
grammar's structural exponent. The deep arm resolved the c₂ rung: L1
(Bernoulli c₂ broken to −0.0252, sharpest tail p = 0.717) fell 6.3
lnL behind L2 — the tail did NOT out-vote the deep. L3 (c₂ and c₃
both Bernoulli, first break at c₄) fell 4.5 behind. **The round trip
is now an empirical statement: the dressing exchange has two legs.**

**Part B — the measured-ambient leg (O16a): both bars PASS, the 6F
partial RESOLVES.** Vertical PGmax −61.68 vs global −59.05 (measured
environments improve the rule; gap to undiluted F4 narrows 5.2 → 2.5).
Plain: **PGmax −100.51, PGno −105.94 — BOTH cross the original ≤ −100
bar** that 6F's fiducial-e run missed by 8. The miss is now explained
and closed by measurement: the 0.02 fiducial was too soft (measured
medians s = 0.9317 maxclust / 0.9756 noclust vs 0.8681). Zero new
parameters anywhere. Honest residuals carried: (i) the galaxy leg
alone still rewards sharper gates monotonically (F4 = s→1 stays its
formal ceiling, +2.5 ahead vertical) — the gate AMPLITUDE is pinned by
the binaries, not by the galaxy leg; the content is the two-system
joint. (ii) The noclust column (more isolated) fits better than
maxclust — same direction, sharper; the column choice is a ~5-lnL
ambiguity we carry, not hide. Updated joint two-system ledger vs BE:
**≈ −59.9** (gal −61.68 ⊕ binary +1.8 from 6G's −0.88 lnL).

**Part C — the AMB quadrupole record (O16c):** per-eN solar gates
g(1.0/1.2/1.4) = 0.135/0.112/0.094; q = −0.0855/−0.0992/−0.1120; amb
within +0.38% of BE at 1.2 (gates: Newton 0, simple-vs-4K 0.1%,
BE-vs-5S-cache 0.00%, hi-res 0.0%). Lock join with 6G's α̂ = 1.060 ±
0.024: **Q₂ = 3.60e−26 s⁻² = 4.0× Cassini** — BE-grade, no rescue, as
pre-stated; the MI door remains the standing escape (4L).

**Priority scout (Haiku, post-6H, scout-level):** the Bernoulli-break↔
leg-count correspondence and the (2n+1)-dispersive selection argument
both NOT FOUND (also re-confirmed: no environment-dependent exponent,
no mixing-family analogue; Famaey & Durakovic's "no existing clear
derivation" line re-surfaced). Quoted as scout-level negatives.

**Verdict: SUCCESS — the derivation round paid out in full. The
grammar's discrete exponent was measured (L = 2, both flanks rejected),
the measured environments strengthened every galaxy number at zero
parameters and retro-resolved the one bar the rule had missed, and the
quadrupole record is on file at 4.0× with no escape claimed. 6J
(bootstrap grade of the PGmax lead, 40 paired reps) launched — running
at close of entry.**

## Stage 6J — bootstrap grade of the measured-ambient lead (2026-07-23, O16b)

[calcs/stage6j_ambboot.py](calcs/stage6j_ambboot.py) →
data/stage6j_ambboot.txt. 40 paired galaxy-resample reps (dv-ON
vertical, lensing noise redraws, rng 53 — 6C/5N machinery verbatim),
referent = the rule as prescribed (per-galaxy Chae maxclust gates,
median s = 0.9317; lensing carries the matched median). Full-fit
regression PASS (BE exact −12152.49; AMB −61.62 vs 6I's −61.68).

**AMB(pgmax) − BE = −56.71 ± 35.65, AMB better in 37/40 (92.5%);
percentiles 16/50/84 = [−88.7, −59.9, −15.7].** Comparators: 6C F4
−57.4 ± 38.3 (37/40); 5N gm −29.3 ± 53.0 (29/40). The ambient-gated
rule with measured environments carries the SAME bootstrap grade as
the best undiluted function ever measured here (F4) — while, unlike
F4, also passing the binaries. Still a strong LEAN, not a detection
(3/40 reps flip; spread ±36), and quoted as such.

**Verdict: SUCCESS — O16b closed; the two-system rule's galaxy lead is
bootstrap-graded at the program's top mark (37/40), tied with F4 and
binary-compatible where F4 is not.**

## Stage 6K — the desktop analog, v1: rate-based realizations EXCLUDED, the pre-committed strike logged (2026-07-23, O5 lab leg)

User asked the right question before the run: "would open-quantum-systems
math actually work on this? What if the math won't predict this?" —
answered by pre-registration: mapping, estimator, fingerprint bands
(p_s ∈ [1.6,2.4], p_n ∈ [−2.4,−1.6], κ-slope ∈ [−0.2,0.2]), outcome
tree, and credence commitments were committed (ddc83bc) BEFORE
execution. [calcs/stage6k_analog.py](calcs/stage6k_analog.py) →
data/stage6k_analog.txt. One post-commit patch: a float-range guard
(e^700 overflow) — numerical, no design change, disclosed.

**The lab-native identity (G0, new, exact):** the rule rewrites as
**β = ½·tanh²(x_loc/2)·e^(−2x_amb)** — since 1/(2n+1) = tanh(x/2) and
n/(1+n) = e^(−x), the squared ambient gate IS the Boltzmann cost of
borrowing two ambient quanta, one per leg (gal 0.7536 / bin 0.1171
reproduced to 1e-16). Both factors are standard thermal functions; only
their per-vertex assignment is nonstandard.

**The model:** exactly solvable birth–death NESS of a Kerr mode
(E_n = n·w0 + (K/2)n(n−1)) between a frequency-blind source channel
(pins n₀ = n_S; the β=0 endpoint) and an ambient channel in two
configurations: CFG-V (vanilla Davies bath) and CFG-M (mediated
two-vertex jump: source quantum + ambient mismatch quantum; its
detailed balance alone → dressed Gibbs = the β=1 endpoint, verified
2.6e-13). All gates pass (endpoints exact, K=0 closed form 0.0e+0,
truncation 0.0e+0).

**RESULT — pre-registered verdict: ALT (with a G4 pathology in CFG-M):**
- CFG-V: λ ≈ 0.50–0.54 FLAT in both shares (p_s = +0.04, p_n = −0.00),
  κ-slope +0.50 — pure rate-weighted mixing. As a gravity law this is
  CONSTANT β — **already sky-excluded** (5T: constant-β family dead,
  binaries 24/24 + the regime decomposition).
- CFG-M: λ RISES with occupation (p_n = +0.57 — opposite sign to the
  grammar's −2 AND to the sky's deep→0 demand) and overshoots λ > 1 in
  the hot-ambient corner → **G4 estimator-pathology FAIL flagged** —
  its exponents are parameterization-contaminated and NOT quotable as
  measurements (the honest AMBIG sub-branch).
- Robust cross-config negative: **κ-dependence 0.3–0.5/decade in both —
  the grammar's rate-freedom does not emerge from jump-rate
  competition anywhere in the scan.**

**The insight that reshapes the lab program:** 6H's q_loc came from the
dispersive PULL — a COHERENT second-order Hamiltonian shift — while
v1's model class contains only INCOHERENT jump competition. The strike
therefore lands on the **rate-based reading** of the grammar
(excluded: it yields constant or wrong-sign running, both sky-dead —
the pre-committed −10-point strike applies to that reading), while the
**coherent-pull reading remains untested**: v2 must compute the
driven-Kerr susceptibility-peak admixture (full Liouvillian + quantum
regression, n_max ~ 100 — feasible on this machine) — the object the
JC (2n+1) anchor actually describes. Bath-microphysics conditional
credence: ~20–25% → **~15%**, surviving mass now concentrated on the
coherent-pull/non-Markovian corner. The user's "what if the math won't
predict this?" has its live answer: the first standard class doesn't,
the strike is logged per pre-commitment, and the question sharpened
into a specific next calculation instead of dying vaguely.

**Verdict: SUCCESS as a test (pre-registration honored end-to-end,
gates clean, every branch informative) — STRIKE LOGGED for the
rate-based mechanism reading; the analog program continues at v2
(coherent susceptibility), now the sharpest open item under O5.**

## Stage 6L — CORRECTION #13: the leg-count "measurement" deflates under its own bootstrap (2026-07-23)

The error bar the round owed came back and took the headline with it.
[calcs/stage6l_legboot.py](calcs/stage6l_legboot.py) →
data/stage6l_legboot.txt. 40 paired galaxy-resample reps (6C/6J
machinery, rng 53), all three L-laws per rep; full-fit regressions
PASS (d = +1.56/+0.36/+0.40 vs 6I).

**Result: d(L1−L2) = +9.94 ± 21.39 — L2 better in 29/40 (72.5%) = a
LEAN; d(L3−L2) = +1.56 ± 13.13 — L2 better in 21/40 (52.5%) =
UNRESOLVED; L2 strictly best of three in 10/40 (25%).** Percentiles:
d12 [−10.7, +11.1, +28.4]; d32 [−8.9, +0.4, +16.9].

**CORRECTION #13 (logged in paper Appendix A; v2.4 → v2.5):** the
claim "L = 2 MEASURED against both flanking integers" (6I, v2.4) was
too strong. Under galaxy-population resampling the three-way ordering
shuffles: the one-leg rejection survives at lean grade (carried by the
c₂ Bernoulli break), the two-vs-three-leg contest carries essentially
no population-grade information, and the strict three-way win is
realization-lucky. Diagnosis: the discrimination lives in the
ultra-deep arm, which a 153-galaxy sample populates thinly — the 3A
realization-systematic lesson recurring at the function-structure
level. **Reframed everywhere: L = 2 is the point-preferred process
order with a lean-grade one-leg rejection; the round trip is the
PREDICTED order with a first empirical vote, not a measured one.**

What is UNAFFECTED: the 6H Bernoulli-break↔leg-count correspondence
(exact mathematics — the instrument is sound; the sample can't read
the dial at grade); the measured-ambient results (6I Part B) and
their 37/40 bootstrap (6J — a different, stronger statement); the
quadrupole record; the 6K strike. Deeper rotation-curve samples read
the rung directly — added to the future-data ledger beside c₄.

**Verdict: the stage — SUCCESS (it did exactly what it was built to
do); the claim — DEFLATED, correction #13 logged. The program's
yardstick applied to its own newest headline, same day.**

## Stage 6M (= 6K-v2) — the structured-bath calculation: AMBIG per bands, with two clean beyond-band observations (2026-07-24)

The physical interpolation v1 lacked: ONE thermal bath of finite
spectral resolution b (Ohmic-weighted Gaussian golden-rule kernel) on
the Kerr ladder — b ≪ K resolves the comb (Davies → dressed Gibbs =
the β=1 endpoint), K·n ≪ b ≪ ω₀ cannot (singular → populates at the
source frequency = the β=0 endpoint) — and the overall rate κ cancels
EXACTLY in the NESS, removing v1's κ objection by construction.
Pre-registered at 2b3bf78 BEFORE execution (share band, resolution-
collapse band, outcome tree, credence commitments).
[calcs/stage6m_analog2.py](calcs/stage6m_analog2.py) →
data/stage6m_analog2.txt. Gates: Davies endpoint 0.0e+0; kernel-vs-
delta 0.00%; source endpoint 0.67%; truncation 0.0e+0; pull lemma
exact (thermal comb centroid = 2K·n̄ — the coherent pull weighs the
state's STATISTICS; tanh(x/2) lives in the (n+½) bookkeeping, not the
dynamics). Two numerics patches post-commit (np.trapezoid rename;
kernel/grid coverage with analytic branch beyond) — disclosed, no
design change.

**Pre-registered verdict: AMBIG.** The share test FAILED decisively
(λ/tanh² spread 205× at the most transitional column — no share
organization anywhere), but the resolution-collapse test ALSO failed
(5.3× vs the <2 band) — because at practical K the endpoint gap
(K-order) is smaller than the kernel's frequency-sampling systematic
(b²-order) through the mid-scan: λ overshoots (−1…−122) once the wide
kernel feeds the mode off-resonant thermal content. The estimator
drowned in exactly the way the AMBIG branch anticipated: kernel-detail
sensitivity → **a pseudomode v3 (mode + damped auxiliary mode, exact
Lindblad, true Lorentzian bath — no kernel heuristics) is required
before any lab-grade claim.**

**Beyond-band observations (logged as observations, not verdicts;
no pre-committed credence move fires on AMBIG):**
1. **In the clean region the monotonicity is OPPOSITE to the grammar:**
   λ(x₀) at fixed r runs 1.000 → 0.435 as occupation falls (K=0.005,
   r=2.8) — occupied/classical modes sit FULLY self-consistent (their
   larger level spread K·√(n(n+1)) is easier to resolve), sparse/
   quantum modes get source-pinned first. The grammar runs the other
   way (occupied → source-pinned, sparse → dressed). Sign-level, robust
   to the estimator issues.
2. **The gravity translation of the resolution dial (analytic, part
   iii): ρ = x(√ν−1) peaks ≈ 0.26 at x ≈ 1–1.5 and vanishes BOTH deep
   (√x) and in the tail (x·e^(−x)/2)** — a resolution-β would be
   maximal at the transition and zero in the tail; the sky (5T) wants
   β = 0 deep AND transition with the tail alone at ½–¾. Wrong shape
   in both non-deep regimes for every O(1) bath-width ζ.
Both lean against standard-bath realizations of the grammar; per
discipline the ~15% bath-microphysics conditional HOLDS (the AMBIG
branch carries no committed move) with the lean noted.

**Verdict: the stage — SUCCESS (pre-registration honored; the AMBIG
branch fired exactly as designed); the mechanism question — still
open, now with a fully specified closing instrument (pseudomode v3)
and two sign-level leans against, on the record.**

## Stage 6N (= 6K-v3) — the pseudomode calculation CLOSES the lab leg: CLOSE-OPP, credence 15% → ~8% (2026-07-24)

Pre-registered at 6724597 BEFORE execution.
[calcs/stage6n_pseudomode.py](calcs/stage6n_pseudomode.py) →
data/stage6n_pseudomode.txt. All gates pass (linear composite 2.7e-7;
theorem gate exact; truncation 4.7e-7; solve residual 2e-16;
positivity −3e-17; GB1/GB2 machine-zero; GB3 Davies endpoint 1.05).

**THE CANCELLATION THEOREM (stated in design, confirmed exactly):**
for any FLAT reservoir — occupation constant across its line — the
up/down ratio of every mediated transition is the same Boltzmann
factor, so a self-shifting mode stays source-locked EXACTLY: the
Lorentzian lineshape cancels in the weak-coupling ratio, and every
composite eigen-transition pays the same constant ratio in the
strong-coupling secular limit. **Part A (exact two-mode Lindblad, the
buildable circuit-QED configuration: Kerr mode + damped thermal
filter): max |λ| = 0.0065 over the FULL scan (three x₀ × three κ/K ×
three g/κ; g-independent to 4 digits)** — the non-secular residual is
sub-1% and unorganized; the single-filter fridge experiment measures
zero, by theorem. (Secular source-locking is presumably quantum-optics
folklore — stated as such, no novelty claimed without a scout.)
Corollary worth one line: **which occupation a self-shifting mode
takes is decided ENTIRELY by whether the bath's thermal (KMS)
structure is resolved across the self-shift** — the C&T-vs-boot
dichotomy in one sentence.

**Part B (the KMS-contrast carrier, two filters with local thermal
occupations — the minimal resolved-structure bath):** λ runs 0.22–3.06
across the scan — the carrier WORKS (GB1/GB2/GB3 establish the
organization: zero without contrast, zero with equal occupations,
Davies with a fine bank) — and it is **NOT share-organized: λ/tanh²
spread 34× in the best window (band was 1.69)**. Honest note: λ RISES
with x₀ in all four rows — superficially the grammar's direction —
but the organization is by BATH GEOGRAPHY (which lines exist, where,
at what occupations: contrast 1−e^(−δx₀) at fixed δ), not by the
mode's own occupation share; the grammar's β depends on the mode
alone. Recorded, not spun; it does not rescue the share form
(quantitatively 34× off) and 6M's resolution dial ran the OTHER way —
both are standard physics, neither is the grammar.

**Part C (analytic):** the physical contrast dial in gravity,
C = [n(x)−n(x+Δx)]/n(x) with Δx = x(√ν−1): 0.80 deep → 0.32 at x=1 →
0.02 at x=5 — again transition-weighted, VANISHING in the tail. The
sky's tail-only β ≈ ½–¾ (5T) is unreachable by any monotone map of
this dial: **within standard bath physics, the tail is the regime
LEAST able to sustain admixture — the sky's pattern is
anti-standard.**

**PRE-REGISTERED VERDICT: CLOSE-OPP. The lab leg CLOSES.** Three
realization classes — jump-rate competition (6K), golden-rule
structured bath (6M), exact-Lindblad filter dynamics (6N) — and none
produces the grammar's share-squared, rate-free, mode-local gating;
the one exact theorem says the flat-bath endpoint is source-locking,
and the one physical carrier (KMS contrast) is bath-geographic with
the wrong sky shape. Per the pre-commitment: **bath-microphysics
conditional 15% → ~8%. NO v4; a fridge experiment is moot — it would
re-measure the standard physics computed here (Part A IS the
buildable configuration, and its answer is zero).** The mechanism's
remaining homes: (a) horizon-specific non-Markovian physics (the dS
"bath" is causal structure, not a filter bank — explicitly not
fridge-buildable; the ~T bath-width loophole noted), (b) the
MI/trajectory reading (which never lived in a bath), or (c) no
microphysics at all — the ambient-gated function as pure sky
phenomenology, which is the ~92% complement and takes nothing from
its function-level record (joint two-system pass, 37/40 bootstrap,
measured environments, tail postdictions).

**Verdict: SUCCESS — the question the user asked two nights ago
("what if the math won't predict this?") is now fully answered by
three pre-registered calculations and one theorem: it doesn't, we
know exactly why, and the program knows precisely where the mechanism
must live if it lives anywhere.**

## Stage 6O — frozen-bath Test 1 (galaxy coherence): AMBIG — the identifiability boundary survives shape leverage (2026-07-24)

The frozen-bath reading (bath correlation time ~ 1/H ≫ orbital times;
mode wavelength ≫ galaxy ⇒ one coherent draw per galaxy) predicts a
per-galaxy MEAN channel with the thermal shape t(x) = e^(−x/2) —
against the x-FLAT astrometric channels (distance: g_bar-invariant
uniform log g_obs shift; inclination: uniform — both derived exactly).
Pre-registered at 53a7395; ESTIMATOR CORRECTED post-commit,
pre-results (disclosed in-script): the committed MAP z-channel was
monotone-degenerate (a new per-galaxy parameter always lowers a
penalized objective — the first launch crashed on its own G0 gate);
fixed by EXACT z-marginalization (per-galaxy 2×2 Occam determinant,
identically zero at A = 0 ⇒ the historic baseline objective and G0
regression preserved). [calcs/stage6o_cohtest.py](calcs/stage6o_cohtest.py)
→ data/stage6o_cohtest.txt.

**Result: G0 EXACT (−12152.49, d = −0.00) — and both injection gates
FAIL, so the pre-registered verdict is AMBIG.** The profile runs away
to the grid edge (Δ = −835 at N = 8.4, no interior minimum); the
recovery injection over-recovers (Â = 0.16 for A_inj = 0.10); and the
decisive null: the within-galaxy SHUFFLED template still gains −82.
The z-channel is soaking up real per-galaxy radial residual structure
(M/L gradients, warps, non-circular motions) that dwarfs any
thermal-draw amplitude. One informative number for the record: the
true x-shape beats the shuffled one ~10× (−835 vs −82) — the
residuals ARE strongly x-monotone across galaxies — but that is
exactly what ordinary radial systematics produce too, and the gates
correctly refused to let it be read as a detection.

**Verdict: AMBIG (the honest 4W echo, now established at shape
level): the galaxy-coherent channel is UNMEASURABLE on SPARC against
real per-galaxy systematics — radial-shape leverage does not unlock
the boundary. No credence move (pre-committed). The frozen reading's
galaxy-side test is closed as boundary-limited; its binary-side test
(6P, orthogonal noise model) is the live one.**

## Stage 6P/6P-b — frozen-bath Test 2: the edge preference is NOT shape-specific — the s-flat control WINS; no frozen evidence; a v7 error-model finding logged (2026-07-24)

Pre-registered at ac89f1b (6P) and 15c1e81 (6P-b), both before
execution. [calcs/stage6p_frozen.py](calcs/stage6p_frozen.py) +
[calcs/stage6pb_extend.py](calcs/stage6pb_extend.py) →
data/stage6p_summary.txt, stage6p_frozen.txt, stage6pb_extend.txt.
Implementation: per-system draw on the boost occupation (proxy-grade,
disclosed), patched into the full v7 machinery (exact-count patches;
fz = 1 at SQN→∞ recovers the unpatched model; same seed ⇒ identical
population/orbits/draws across the grid — perfectly paired).
G1 baseline: α̂ = 1.110 interior, dlnL(Newton) = +97.3, wr = 0.20
(3S-era regression PASS).

**6P: the likelihood WANTS the component — monotonically.** Δ vs off:
+1.2 (N=60), +6.3 (N=20), +11.1 (N=8, GRID EDGE); α̂ interior
throughout; Newton rejection strengthens (+97→+108). The PREFERRED
branch fired AT THE EDGE — flagged immediately per the correction-#4
edge discipline, and 6P-b supplied the three missing legs.

**6P-b: the discriminator lands hard.** (i) Extension: STILL no
turnover — +16.2 (N=4), +23.5 (N=2): the appetite for scatter runs to
70%-grade draws, deep into generic-broadening territory. (ii) THE
s-FLAT CONTROL (matched generously at the wide-bin amplitude): **+37.0
at matched-8 — BEATS the shaped +11.1 by 25.9 lnL** (flat-20: +19.0
vs shaped +6.3). Given equal strength, the data prefer the scatter
WITHOUT the frozen s-shape. (iii) Seed 101: shaped Δ = +13.5 —
sign-consistent: the appetite is real across realizations; its SHAPE
is not the frozen one. Formal verdict per the written bands: AMBIG
(spec = −25.9 sits outside both bands — the pre-registration did not
anticipate flat WINNING); substantive verdict: **GENERIC, decisively —
the 6P gain is s-blind broadening absorption, NOT frozen evidence.**

**The v7 error-model finding (real, useful, non-mechanism):** the
pipeline underfits an s-FLAT per-system velocity scatter worth ~+37
lnL (echoing the 3E/3J smear residual — the fitted broadening that
measured mass errors could not supply). α̂ exposure across all these
scatter variants: 1.110 → 1.190 max = **≤ +0.08, inside the quoted
±0.11 systematic** — the anchor and the α measurement stand as
published. Queued as a low-priority error-model refinement, not a
re-opening.

**The frozen-bath reading after both tests:** Test 1 (6O)
boundary-limited (unmeasurable vs galaxy systematics); Test 2
control-rejected (its fingerprint subdominant — the data prefer
non-frozen scatter shapes). No formal strike (the pre-registered
STRIKE band did not fire), no support anywhere; the ~8%
bath-microphysics conditional HOLDS with its frozen sub-branch
weakened — both of its accessible near-term signatures are now spent.
Remaining mechanism homes unchanged: the horizon-side non-Markovian
derivation (theory), MI/trajectory (DR4), or sky-phenomenology.

**Verdict: the stages — SUCCESS (the control did exactly its job: a
lesser protocol would have announced a two-seed, +11/+13 "detection"
tonight; ours converted it into an error-model fact in four hours);
the frozen reading — NO SUPPORT from either test, honestly closed at
current data.**

## Stage 6Q — the measurement ledger + the world-table consistency audit (2026-07-24)

**The user's question from last night ("are we keeping the values? could we
retroactively fit an answer equation?") executed as infrastructure.** No new
fits — this is the bookkeeping stage that turns four layers of record (stage
outputs, NOTES, paper scorecard, git) into ONE machine-checkable table with
supersession discipline, and then scores every candidate law against every
ledger row mechanically.

**Deliverables:**
- **[LEDGER.csv](../LEDGER.csv)** (repo root; spec in [LEDGER.md](../LEDGER.md)):
  75 hand-curated rows — every headline number with value (NOTES convention,
  direction words), stage, script, output file, data dependence, status
  (CURRENT / CO-QUOTED / SUPERSEDED / RETRACTED), and supersession pointer.
  Hand-curated BY RULE: auto-scraping stage outputs is exactly the mechanism
  that would resurrect a superseded number (the wrong-convention α=1.54, the
  inflated +264, the e_N=1.9 field, the "SPARC leans simple" artifact, the
  "L=2 measured" claim — all present as visibly-marked SUPERSEDED/RETRACTED
  rows so they can never quietly come back).
- **[calcs/stage6q_worldtable.py](../calcs/stage6q_worldtable.py)** →
  [data/stage6q_worldtable.txt](../data/stage6q_worldtable.txt): six
  mechanical gates, then the world table.

**Gates (all PASS on 75 rows):** G1 ids/status/fields; G2 every cited script
and output file exists on disk; G3 every supersession pointer resolves to a
live row; G4 no two CURRENT rows claim the same quantity; G5 fifteen ledger
values grepped verbatim out of their stage output files; G6 every world-table
cell cites a live (non-retracted) ledger row. The gates caught three of my
own transcription slips on the first run — a CSV quoting bug that shifted the
α-final row, a rounding that didn't match 5N's file (−29.27±52.98, now quoted
at file precision), and a wrong output-file pointer for 6G — which is the
audit doing for the ledger what the fit gates do for the physics.

**The world table** (13 laws × 7 test columns, per-cell ledger provenance):
- **Formal sky vetoes:** Newton (everywhere, 24/24 + 2000/2000), F4
  (binary α-edge 6/6), boot (binary edge + vertical collapse), the pointwise
  drive-weighted rule (0/6), no-EFE MI (12/12). These five are out.
- **No formal veto:** simple (but hier-galaxy-rejected at strong-lean grade
  −99 — the two-system tension in one row), BE, ν_p(0.65), gm, F1/F2, F3
  (all binary-lean-against and a₀-strained at +4.3…+6.3σ), **AMB**, and
  EFE-respecting MI (the Saturn door).
- **The mechanical sentence the table proves:** AMB is the only candidate
  simultaneously (a) un-vetoed on the binaries (−0.88±2.66 tie, interior
  6/6), (b) at the top galaxy bootstrap grade (37/40 — same as F4, which IS
  binary-vetoed), (c) best on the binary temperature row (+1.6σ), and (d)
  postdicting both measured tail exponents (0.689/0.529). It carries the
  post-hoc flag and the same 4.0× Cassini quadrupole as every MG member.
- **Independence made explicit:** all galaxy columns share SPARC-153
  (+Chae21) = one data vote; all binary columns share EDR3-14071 = one data
  vote; Cassini is the third, independent, and 4.0–5.8× against every MG
  member (escape: mi_t, which ties the binary contest). "Passes both
  systems" = exactly two independent votes, not five.
- The ladder digit (½ vs ⅓ vs sharper) is confirmed OPEN at population
  grade by the assembled table — no row closes it.

On the "answer equation" half of the user's question: the audit shows the
retro-fit is already done in the legitimate direction — the rigid candidate
(AMB + a₀ = cH₀/2π + α = 1) postdicts the ledger's sky rows; what it does
NOT absorb is exactly the residue the table isolates (Saturn ×4, the binary
a₀ +1.6σ, N ~ 20–60, the mechanism). Any future function proposal can now
be scored against the full table in one script run instead of a NOTES
archaeology session.

**Discipline note:** the ledger is append-only in spirit — replacing a number
means flipping the old row to SUPERSEDED with a pointer, never deleting.
Every future stage that produces a headline number should add its row in the
same commit.

**Verdict: SUCCESS — infrastructure stage; no physics changed, every number
now has one canonical home, and the "which function survives what" question
is answered by a script instead of memory.**

## Stage 6R/6S/6T — the resolution bath: the horizon-side round (2026-07-24)

**O5's top item executed: the 6N cancellation-theorem corollary applied to
the one thing that makes the dS bath unlike every lab bath — it has ONE
scale.** The corollary said: which occupation a self-shifting mode takes =
whether the bath's KMS structure is resolved across the self-shift. For lab
baths the resolution scale is arbitrary geography (why 6K/6M/6N all failed
with transition-peaked, tail-vanishing sky shapes). The dS bath's occupation
structure varies on the scale of its own temperature, and T_dS = H/2π IS the
a₀ scale — so the resolution criterion is a forced function of the field:

    R = (ω_tot − ω_src)/T_dS = νy − √y

**Stage 6R ([calcs/stage6r_resolution.py](../calcs/stage6r_resolution.py),
pre-registered at 85dcc72 BEFORE the fits; all gates PASS):**
- Deep: R → x²/2 (quadratic vanishing) ⇒ EXACT source-locking — the measured
  β = 0 of every deep+transition arm (5T, 5R) is the deep half of this
  reading. Tail: R → y ⇒ dressing on. **Crossover R = 1 at y* = 1.88 — the
  5T turn-on arm (y > 1).** No lab bath puts it there; this bath can't put
  it anywhere else.
- The zero-parameter function: β(y) = ½·R²/(1+R²) in the 5P frequency mixing
  (Lorentzian profile = the one representative choice, flagged). Exact:
  **Bernoulli ladder preserved through c₄; break at c₅ = −1/16** (sympy
  exact + mpmath 50-digit; deepest ladder preservation of any candidate —
  AMB breaks at c₃). **Tail = the gm argument exactly (p = ¾).** ν_R(1) =
  1.539 (F4-grade sharp). Solver gates: residual 8e-11, mpmath agreement
  1e-12, wfac=0 → BE at 9e-14.
- **The Deser–Levin exact-temperature matrix CLOSED analytically:**
  source-frequency branch has NO finite fixed point (runaway; anti-Newtonian
  catastrophe); total-frequency branch = boot deep (dead 5F/5M) + invisible
  constant floor n_BE(2π) = 0.00187 (pure G-renormalization); **a geodesic
  (free-fall) detector has T_U = 0 exactly ⇒ pure C&T — the binary β < 0.03
  (5R) is the free-fall answer, now derivation-grade.**
- Pre-registered bars + the honest prior risk (no ambient gate, F4-grade
  transition ⇒ leans binary-REJECT) committed before execution.

**Stage 6S (galaxies, [calcs/stage6s_resngal.py](../calcs/stage6s_resngal.py),
6I machinery, BE regression gates d = −0.00 both treatments): STRONG.**
Vertical-hardened **−58.59** vs BE (bar ≤ −55; AMB fiducial −59.05,
measured-ambients −61.68; F4 −64.2); plain hier **−113.72 = the largest
plain-treatment lead ever recorded** (F3 −111.4, F4 −108.7, AMB-measured
−100.5/−105.9). A zero-parameter pre-hoc function leads the ladder.

**Stage 6T (binaries, [calcs/stage6t_resnbin.py](../calcs/stage6t_resnbin.py),
6G machinery verbatim: tables gated 0.01% spherical identity, corrected
velocities, 6 seeds): SHAPE REJECTION — the pre-registered discrimination
fired.** α̂ edge-rides at 2.00 in 5/6 seeds (interior 1/6 at 1.67); 0/6
seeds prefer it; mean **−7.38 ± 2.08 vs p050** = squarely the
eight-function band (gm −8.50, rb4 −6.96, dwf −9.64; AMB −0.88); Newton
still dead (+82.7 ± 1.8); w_rad = 0.2 in 6/6. κ = 1.351 → the a₀ row is
edge-invalidated (formally +10σ at the edge value — quoted as invalid,
direction noted).

**THE TRIANGULATION IS NOW COMPLETE AND MECHANICAL:**
- local-resolution running alone: galaxy-excellent, binary-DEAD (6T);
- pointwise ambient weighting alone: binary-DEAD (6D);
- the two-factor rule (local quantumness × ambient classicality): the only
  structure that passes both systems (6G/6F/6I).
Three independent exclusions now *require* the grammar's two-factor form —
the ambient gate is not a modeling choice, it is what the binaries demand.
What survives of the resolution reading is its deep half (why deep+transition
source-lock everywhere, parameter-free); what died is "resolution alone,
system-blind." The mechanism ledger: the ~8% bath-microphysics conditional
now lives specifically in "one-scale non-Markovian dS bath supplies the
local factor; the ambient gate needs its own horizon-side derivation" — the
sharpest remaining O5 formulation.

Numbers appended to LEDGER.csv (per the 6Q standing rule): galfn-resn,
binfn-resn, mech-resolution; world table updated with the RESN row and the
triangulation line ([calcs/stage6q_worldtable.py](../calcs/stage6q_worldtable.py)
rerun, gates green).

**Verdicts: 6R (derivation) — SUCCESS (five exact structures, all gated;
two analytic closures; the free-fall reading of β_bin = 0). 6S (galaxies) —
SUCCESS at STRONG grade (−58.6/−113.7; plain record). 6T (binaries) —
SUCCESS as a test, REJECTION as a candidate (the pre-registered
discrimination did its job: the ambient gate is REQUIRED). The round —
DIFFERENT PHYSICS excluded cleanly: the two-system split is NOT local-in-y;
it is system-level, as AMB encodes.**

## Stage 6U — the gate derivation: the KMS cost of borrowed quanta (2026-07-24)

**The O5 question left by 6R/6T — derive the s_amb² gate — executed the
same day, from the same object that gave the local factor.**
([calcs/stage6u_gatederiv.py](../calcs/stage6u_gatederiv.py), all gates
first-run; no new fits — the empirical half is scored entirely against
LEDGER rows.)

**The algebra (exact, sympy):** path-resolving the 6H dispersive loop
separates the (2n+1) pull into its two time-ordered exchanges — the
absorption-side path (weight n: the partner mode must SUPPLY a quantum)
and the emission-side path (weight n+1: the vacuum rides along). Their
ratio is **n/(1+n) = e^(−x) EXACTLY — the KMS/detailed-balance weight.**
A per-leg net-borrowing channel is Boltzmann-taxed at the partner's
occupation; L legs cost **s^L = e^(−L·x_amb) = the Boltzmann cost of L
borrowed ambient quanta**, with L = 2 the loop's own (Lamb-shift) order.
Each vertex of ONE loop carries both measured factors: the local
zero-point share 1/(2ν−1) on one side (6H), the ambient KMS ratio on the
other (this stage). The user's original two questions ("why two?", "why
local quantumness?") now have one answer: one second-order loop, two
vertices, two faces per vertex.

**The uniqueness selection (the derivation's empirical half):** of the
loop's algebraically natural weights, measured facts exclude every rival:
- rate-balance s²/(1+s²): p_gal = 0.607, below the measured tail band
  [0.65, 0.75] — EXCLUDED;
- absorption SHARE [n/(2n+1)]²: p_gal = 0.554 — EXCLUDED (**the data
  choose detailed balance over the share from the same decomposition** —
  the sharpest new fact);
- inverse ratio [1/(1+n)]²: fails BOTH tails, wrong environmental sign
  (the 6E sign lesson formalized);
- raw amplitude n²: p_gal ≈ 11, absurd + ceiling violation — EXCLUDED;
- pointwise s(y_loc)²: EXCLUDED by measurement (6D);
- ambient-free local running: EXCLUDED by measurement (6T);
- s¹/s³: at their 6L bootstrap grades (disfavored-lean / open — the tail
  instrument is blind to the L digit; the Bernoulli rung is not).
**s² is the unique survivor.** Postdiction regression: p_gal = 0.6884 /
p_bin = 0.5280 reproduce 6E's 0.689/0.529.

**Dividends (both quantified):** (1) the source-vs-dressed ambient
frequency convention is DEGENERATE today (Δp = 0.004 at galaxy ambients)
but separates at weak-ambient binaries (Δp ≈ 0.025 at e_N = 0.4) — a DR4
discriminator; (2) the gate is a Boltzmann cost at T_dS = H(z)/2π ⇒ **the
tail exponent runs with redshift at fixed physical environment** (p_gal
0.688 → 0.702 at z = 1) — the gate is now part of kill test #14.

**Priority scout (Haiku, five structures, 15+ query variations):** the
KMS-ratio/detailed-balance gate, the Boltzmann-suppressed EFE, the
reservoir-assistance mechanism for gravity dressing, the path-resolved
JC decomposition in any gravitational context, and any C&T environmental
follow-up — **all NOT FOUND (scout-level)**. Near-misses disclosed:
2024–25 stimulated-graviton papers (2502.10221, 2407.11929 — graviton
scattering, not dressing/EFE); entropic-MOND EFE derivations (1201.4160
— entropic force, no occupation gating); modified-inertia environmental
dependence via time-averaging (2208.07073 — no Boltzmann structure).

**Honest labels, carried:** the KMS algebra and the uniqueness selection
are exact/measured; the borrowing NARRATIVE (frozen one-scale horizon
bath cannot supply dressing quanta ⇒ the occupied ambient reservoir
does) stays reading-grade — its microphysics lives with the
non-Markovian horizon formalism (6K/6M/6N showed standard OQS cannot
produce it). The 6O/6P frozen nulls do NOT strike it: they tested draw
VARIANCE; borrowing enters at mean level. Remaining open: the formal
per-leg pairing Hamiltonian; the untied-L contest (β = ½q^L1·s^L2,
L1 ≠ L2) as a future falsifier; the ceiling ½ stays at its 5P
exchange-symmetric reading.

**Verdict: SUCCESS — the grammar's last underived factor now has the
same grade as the local one: exact algebra from the shared loop, plus
unique survival against every natural rival under already-measured
facts; two new falsifiable rungs (DR4 convention split, z-running
tail); priority scout-clean. The ambient-gated bath is now a DERIVED
structure with one reading-grade seam (the borrowing microphysics),
not a post-hoc rule.**

## Stage 6V — the untied-exponent contest: the 6U falsifier (2026-07-24)

**The derivation was given its own way to fail the same day.** The loop
derivation REQUIRES tied exponents (one L: each vertex carries the local
share AND the ambient KMS ratio together). 6V contests the off-diagonal
cells β = ½·q^L1·s^L2 on the vertical-hardened ladder at the fiducial
gate ([calcs/stage6v_untied.py](../calcs/stage6v_untied.py), pre-reg
6ebf545 BEFORE execution; BE regression d = −0.00; tied-(2,2) regression
vs 6I dd = +0.00).

**Results (Δ vs BE, point grade):** tied (2,2) −59.05; U12 (q¹s²)
−51.94; U32 (q³s²) −59.79; U21 (q²s¹) −62.38; U23 (q²s³) −55.33.
Margins vs tied: U21 +3.33, U32 +0.74, U23 −3.72, U12 −7.11.

- **B1 (tied survival): TOLERATED — no strike, no proof.** Worst untied
  margin +3.33, far below the +5 strike bar and far inside the ±13–21
  population noise 6L measured for these very cells. The tied form
  survives; the grid is statistically flat at this sample (expected
  under correction #13).
- **B2 (instrument split): FAILED as pre-registered — logged.** I
  predicted L1-varied cells move less than L2-varied (deep rung reads
  L1, tail reads L2; vertical fit tail-dominated). Measured: 3.92 vs
  3.52 — the prediction is wrong. Autopsy: the deep-series mapping is
  exact algebra and stands, but the likelihood ALSO reads L1 through
  the TRANSITION-suppression profile (q^L1 sets where in ν the
  admixture turns on), and the transition votes β = 0 (5T). U12's
  −7.11 — the only large signal on the grid — is exactly that channel
  punishing one-leg. The failed bar impeaches my model of the
  instrument, not the loop pairing; it goes in the ledger as a failed
  pre-registered bar regardless.
- **Dividend:** L1 = 1 is now disfavored through a SECOND channel
  (transition shape, −7.1 point grade), independent of the deep
  Bernoulli rung the 6L bootstrap said this sample cannot read. The
  one-leg rejection (6L: 29/40 lean) gains a convergent line of
  evidence; L1 = 3 vs 2 stays blind (+0.74).

**Verdict: the contest — SUCCESS (both bars answered at their
pre-registered grades); the tied-loop form — SURVIVES (tolerated, not
proven; the untied family buys at most 3.3 lnL of point-grade noise);
B2 — my prediction FAILED and is logged as such (the instrument reads
L1 via the transition, a better instrument model than the one I
pre-registered).**

## Stage 6X — the borrowing dynamics: SUPPORT (2026-07-24, the first
## constructive mechanism result)

**The user's challenge ("your time estimates are human-calibrated —
decompose it and try") executed same-afternoon: the borrowing
narrative's per-leg factor tested in exact closed-system quantum
mechanics — and it PASSED.**
([calcs/stage6x_borrow.py](../calcs/stage6x_borrow.py), pre-reg
67626f7; amendments post-commit PRE-RESULTS logged in-script per the
6O precedent: (a) two linearly-coupled harmonic modes carry NO
occupation pull — the first crash MEASURED that, exactly λ²/Δ flat —
so the (2n+1) anchor was re-seated on the qubit sub-case where it is
now exact (ratios 1.000/0.999/0.998); (b) channel-physics split:
resonant = Rabi/λ-independent, detuned = virtual/λ²; (c) NB budget.)

**The setting the 6N survivors permit:** frozen horizon bath
(correlation ~1/H) ⇒ on orbital timescales the only dynamical
reservoir is the ambient — the extreme non-Markovian limit is CLOSED
dynamics. Kerr-anharmonic local mode, one thermal ambient mode,
dress-ward transition resonant (dressing must ABSORB a real ambient
quantum), source transition detuned.

**Result: the dress-ward weight = ½·P(n_amb ≥ 1) = ½·n/(1+n) to
0.2–2.2% over n_amb = 0.25–8** (0.09979 vs 0.10000 … 0.43449 vs
0.44444). Form regression: **ratio slope +0.989, rms 0.0046 — 6.2×
better than share, 29× better than raw n — the pre-registered SUPPORT
bar cleared decisively.** G2a: λ-independence 0.997 (real-exchange
Rabi channel). G2b: the detuned virtual control scales λ² (3.51) and
carries RAW-n weighting — **the two channels carry different forms,
and the sky's ratio² (6U) identifies the borrowing as the
REAL-EXCHANGE channel** — consistent with a frozen bath that cannot
absorb virtual imbalances. G3: dephasing the exchange suppresses the
channel ×1.7 (the 6N locking direction). GL (sympy exact): **P(n ≥ L)
= [n̄/(1+n̄)]^L = e^(−Lx) — the geometric-tail identity: the KMS ratio
IS the lending probability, and the gate s^L = P(the ambient can
supply the L quanta the L-leg dressing needs).** The Boltzmann-cost
(6U) and lending-availability readings are the same number; the
availability reading derives the exponent-to-quanta tie directly (the
structure 6V could only tolerate).

**The mechanism chain now:** sky selects ratio² (6U uniqueness) →
exact closed dynamics produces ratio = the per-leg lending
probability (6X) → s^L = P(n≥L) exact (GL) → L = 2 = the loop's own
order (6H exact) → the frozen one-scale horizon is what opens the
real channel and locks the deep regime (6N theorem + 6R). **The
per-leg factor of the borrowing reading is DERIVED at toy grade.**
Seam narrowed to: the dS-side identification of the lending reservoir
(why the ambient field's modes at x_amb are the resonant partners) +
the multimode generalization. Observation, reading-grade and flagged:
the toy's ½ prefactor = resonant equal-time-sharing between source
and dressed configurations — a candidate dynamical origin of the
grammar's ½ ceiling (the 5P exchange-symmetric point, given
dynamics).

**Credence discipline:** the pre-commit specified only the STRIKE
branch (→ ~5%); SUPPORT carries no pre-committed number. Bath-
microphysics conditional HELD at ~8% with the lean flipped from
negative to constructive; the next pre-registered test (dS-side) owns
any upward move. A toy result shows the structure is natural in
closed QM — not that the horizon does it.

**Verdict: SUCCESS — SUPPORT at the pre-registered margins; the
borrowing reading's central factor is now derived-at-toy-grade
mechanics, not narrative.**

## Stage 6W — scalar vs vector EFE: the binaries demand the vector
## composition (2026-07-24; the Saturn round's honest verdict)

**The non-MI door tested and CLOSED by the data.** The idea (pre-reg
2193672 with the consequence map stated both ways; scout: scalar/
direction-blind EFE = apparently novel, never contested on any data;
Chae+21 primary read: their EFE detection is MAGNITUDE-level, direction
"has a minor effect" [1/6 modulation, their words], orientation-resolved
tests = their future work): if the EFE is thermodynamic (scalar bath
state) rather than vector field composition, the solar transition-shell
phantom is spherical ⇒ shell theorem ⇒ Q₂ = the true galactic tide
≈ 7.5e-31 s⁻² = **12,000× below the Cassini cap** — the 4K tension
removed at theory-class level. The binaries are the only data probing
the EFE's shape at the transition: scalar tables B(y) = ν(y+e) vs our
sphericalized vector-QUMOND tables (profile difference: ln(B−1) spans
−0.33…+0.28 over the data window — a large discriminating signal).

**Result ([calcs/stage6w_scalarefe.py](../calcs/stage6w_scalarefe.py),
6 seeds × 2 laws, corrected velocities):**
- **sbe vs vector-BE: EXCLUDED — mean −10.35 ± 2.03, 0/6, interior 6/6
  (a VALID shape rejection, not edge pathology)** — beyond the
  eight-function band; the strongest single-model rejection of the
  program's function contests.
- samb vs vector-AMB: −11.24 ± 1.90, 0/6, interior 4/6 (shape
  rejection + edge).
- The twist, recorded: the scalar reading's temperature row would have
  been the program's best — **α̂ = 1.002 ± 0.049 (parameter-free α = 1
  dead center), κ = 0.934, a₀ = 1.20 ± 0.16 = +1.0σ** — moot under the
  shape exclusion, but a marker of how close the scalar reading came
  everywhere EXCEPT the s-profile. Newton dead +70…+93; w_rad 0.2 6/6.

**What the round actually bought (three things, all new):**
1. **The first measurement of the EFE's composition character.** The
   wide binaries measure the VECTOR (directional-composition) radial
   profile at ~10 lnL per realization over the natural scalar
   alternative — a novel result on the same data everyone uses for the
   boost amplitude (scout-level: no prior contest exists). Precision:
   our tables are angle-AVERAGED, so strictly the binaries demand the
   vector composition's radial profile, not its anisotropy directly.
2. **The Saturn tension is now COMPOSITION-LOCKED, not just
   amplitude-locked** (5I): the same data that calibrate the boost
   demand the vector structure in exactly the regime (ambient ~1.2a₀,
   transition-shell y) that sources the solar quadrupole. The tension
   is data-forced at a deeper level than DHF-2024's formulation-side
   statement.
3. **The sharpened job description for any Saturn-safe MG (#8):** it
   must reproduce the ANGLE-AVERAGED vector-EFE radial profile B(y) at
   e ≈ 1.2a₀ (the binaries demand it) while suppressing the P₂
   component (Cassini forbids it). The logically-open cell — an
   "isotropized-vector-kernel" theory that posits exactly the
   sphericalized profile as an isotropic response — is flagged as
   currently AD HOC (indistinguishable from vector on binaries by
   construction; distinguishable only at Saturn; interesting only if
   the thermal reading can DERIVE the isotropic kernel — an O5-side
   question, noted). Caveat: one scalar prescription was tested (the
   natural ν(y+e)); the exclusion covers the natural scalar-sum class,
   not every conceivable magnitude-only rule.

**Standing after the round: Saturn's doors are (a) EFE-respecting MI
(ties MG on binaries, 4L — unchanged), (b) the ad-hoc isotropized
kernel (needs a derivation to be physics), (c) unknown. The scalar-sum
door is closed. Engagement item logged: the cluster-lopsidedness
literature (1706.07825, 1808.05962) as claimed EFE-direction
signatures — to be engaged at paper level (statistical asymmetries,
contested interpretations; not constraint-grade).**

**Verdict: the contest — SUCCESS (decisively answered, both ways
pre-stated); my scalar-EFE hypothesis — EXCLUDED (killed in ninety
minutes, as it should be); the Saturn problem — SHARPENED, not solved:
the tension is now demonstrably data-forced by our own binaries, and
the honest remaining escape is still modified inertia.**

## Stage 6Y — the reservoir identification: the exclusion theorem, the
## predictions, and the Saturn corollary (2026-07-24)

**The dS-side seam attacked and substantially closed
([calcs/stage6y_reservoir.py](../calcs/stage6y_reservoir.py), pre-reg
0e0f4fc with 6Z, all gates first-run).**

**The exclusion theorem (exact, GY1/GY2):** who lends the L quanta? A
single collective thermal mode gives the measured P(n≥L) = s^L (the 6X
geometric tail). M democratic thermal modes give the negative-binomial
tail 1−(1−q)^M(1+Mq) — and already at **M = 2 the galaxy gate saturates
(0.952 vs the measured 0.754)**, the e_N-dependence that 6I measured
as improving every treatment washes out, and the binary p-postdiction
moves to 0.565 = the band edge the binaries reject; M = 10 is fully
open everywhere. **The measured gate SELECTS M = 1: one collective
ambient mode per system** — converging with the dynamically measured
system-level scalar structure (6D pointwise dead / 6T local dead / 6G
system-level accepted). Two independent routes, one structure.

**The identification (reading, stated plainly):** the collective mode
= the system's barycentric coordinate in the ambient field — one
degree of freedom per system, occupation n_amb = n_BE(x_amb) =
ν(e_N) − 1 = **the environment's own dressing cloud. Systems borrow
from their environment's boost.** Continuity: e_N → 0 sends the
collective mode into the horizon's soft (IR) sector (n → ∞, gate → 1)
— the isolated limit rejoins the one-bath picture; a strong ambient
GAPS the soft sector and closes the gate. The ω = g/c frequency
assignment remains the program's standing §2.4 postulate — inherited,
not re-derived.

**THE SATURN COROLLARY (the user's "different theory" question
answered; reading-grade, three measured legs):** a coupling to the
barycentric collective mode is a TRAJECTORY-STATE coupling — the
modified-inertia class. That class (i) ties vector-MG on the binaries
(4L: mi_t −3.5±3.3), (ii) thereby meets the 6W composition demand at
current resolving power, (iii) carries no capped-type solar quadrupole
(4K's standing escape). **The borrowing mechanics lands the thermal
rule in the one door Saturn left open — by derivation, not
preference.** This upgrades "MI might do it" to "the mechanism's
formalization IS MI-class." Kill test: DR4 eccentricity-resolved
boosts (trajectory-dependence vs field-dependence).

**The prediction ledger extracted (the falsifiers):**
- **P1 (exact): the tail-exponent CEILING p ≤ ¾** — gate = P(n≥2) ≤ 1;
  void asymptote p ≈ 0.72–0.73. One population-grade galaxy tail
  beyond ¾ kills the gate. New, quotable, parameter-free.
- **P2: the ambient ordering** p(e_N) monotone decreasing (curve:
  0.735 at e=0.001 → 0.528 at e=1.2) — being tested NOW at full
  sample power (6Z, pre-registered in the same commit).
- P3: the z-rung (6U) — p runs with H(z).
- P4: DR4 — weak-ambient sharpening toward 0.69 + the
  source-vs-dressed convention split (Δp ≈ 0.03).
- P5: nested ambients gate by the LOCAL total field (operationally
  what 6G/6I already used — consistency, not freedom).

**Verdict: SUCCESS — the reservoir is pinned by an exact exclusion
theorem plus a continuity statement; the mechanism chain (6U→6X→6Y)
now runs: sky selects ratio² → mechanics produces ratio = lending
probability → gate = P(n≥L) of ONE collective ambient mode → the
coupling class is MI = the Saturn door. Remaining reading-grade: the
barycentric identification itself and the ω = g/c assignment; the ½
ceiling candidate stands.**

## Stage 6Z — the ordering shuffle test: ANTI, and CORRECTION #14
## (2026-07-24; the shuffle control does to 6I what the s-flat control
## did to 6P)

**The 6Y ordering prediction (P2) given its full-power in-sample test
([calcs/stage6z_ordering.py](../calcs/stage6z_ordering.py), pre-reg
0e0f4fc; regressions exact: BE d=−0.00, TRUE delta −61.62 vs 6I's
−61.68 dd=+0.06): the pre-registered ANTI bar FIRED.** The true Chae
gate assignment loses to ALL EIGHT distribution-preserving shuffles:
true −61.62 vs null −62.91 ± 0.80, range [−64.59, −62.09], 0/8.

**CORRECTION #14 (logged here and in PAPER Appendix A):** the 6I
claim "measured ambients IMPROVE both treatments" stands as a NUMBER
but its ATTRIBUTION was wrong: the gain over the global gate (−59.05 →
−61.68 controlled; the plain-bar crossing) is **gate HETEROGENEITY —
generic under permutation — not the physical ordering.** Any
assignment of the same gate values does as well or better (shuffles
average 1.3 points BETTER). The in-sample environmental-ordering leg
of the ambient-gated rule's evidence is withdrawn; the plain-bar
"resolution by measurement" is reattributed (real, ordering-blind).
All paper claim sites annotated (abstract, §2.4 grammar passage, the
6Y paragraph, conclusions, Appendix A #14, App B count → fourteen).

**Grade honesty both ways:** ANTI at K=8 carries empirical p ≈ 0.11
(one-in-nine under exchangeability) — a lean against the ordering,
not a demonstration of anti-ordering; and Chae e_N noise dilutes a
true ordering toward null (though noise alone cannot produce ANTI).
Post-hoc shuffles were NOT added to soften the verdict (that would be
p-hacking; K=8 was the pre-registered bar).

**What is and is not struck:** struck — the intra-SPARC cross-galaxy
ordering channel at current grade, and 6I's ordering attribution.
Untouched — everything amplitude-level: the two-system split (the
gate's strongest evidence: galaxy-vs-binary gate difference), the 6G
binary acceptance, the tail postdictions (median-gate), the p ≤ ¾
ceiling (P1), the M=1 exclusion theorem (6Y), the 6X lending
mechanics, and the 6J function-level bootstrap (37/40, with the
heterogeneity asterisk on its measured-ambient component). The
ordering prediction P2 moves FULLY out-of-sample: environment-resolved
DR4-era samples remain its signed falsifier — and after today a
DR4-grade ordering failure would be a REAL strike, pre-stated.

**Verdict: the test — SUCCESS (the control did exactly its job,
second time today); the ordering prediction — ANTI at weak grade,
in-sample leg withdrawn (correction #14); the ambient-gated rule —
evidential basis NARROWED to its amplitude-level legs, which are the
strong ones.**

## Stage 7A (2026-07-24): the Einstein fluctuation test — UNRESOLVED at 4T grade; shot bath EXCLUDED; the discriminant window is occupied by the bump

Pre-reg 21e7d28 (PREDICTIONS.md P6 + [calcs/stage7a_einstein.py](calcs/stage7a_einstein.py)
committed BEFORE execution). The question: does the RAR's intrinsic scatter carry
Einstein's 1909 PARTICLE term? Var(n) = n̄ + n̄² (quantum) vs n̄² (classical wave — the
fluctuation law of ANY classical continuous-field bath, incl. 4F's simple-ν) vs n̄
(pure shot). Through δlogg = δn/((1+n̄)ln10): s² ∝ e^(−x) (γ=1) vs e^(−2x) (γ=2) vs
e^(−x)(1−e^(−x)). Same parameter count, a₀ locked to the mean fit. 4T had fit ONLY
the quantum shape; γ had never been measured.

Results (data/stage7a_einstein.txt; G0 regression: EQ reproduces 4T M1b to 0.001):
- **SHOT EXCLUDED**: its best fit collapses onto the constant floor (N → 1.2e6 = zero
  amplitude), −25.05 behind both rivals. A purely corpuscular bath cannot produce
  deep-end scatter that persists to x → 0. Real exclusion.
- **Quantum vs classical: NO VOTE** — EQ−EC = −0.43 (bar was ≥ 9).
- **Free γ does not localize**: the profile is BIMODAL — minima at BOTH edges (γ ≈
  0.25–0.5 and ≥ 3.5), the middle worst (+4.4 at γ = 1.5); the whole axis spans ≤ 4.4
  units. With the bump masked OR explicitly modeled, γ̂ runs to the 6.0 bound: the
  preferred structure is a deep-concentrated excess + the x≈1 bump (b = 0.083) + floor,
  which beats even the 6-bin free fit by 10.
- **G2 CALIBRATED PASS**: injected quantum truth recovered at γ̂ = 0.98/1.13/0.98
  (mean 1.03), classical truth at 2.18/2.19/1.77 (mean 2.05), min-gap 0.64. The
  instrument WORKS on clean truth — the non-verdict is data-side misspecification,
  not power. G4: simple-ν mean gives the same non-localization.

THE STRUCTURAL FINDING: the quantum–classical discriminant |e^(−x) − e^(−2x)| peaks at
x = ln 2 ≈ 0.69 and carries its weight over x ≈ 0.4–1.2 — and the 4W point-level bump
(the one scatter feature that survives every marginalization) sits at 0.8–1.4. **The
bump occupies the discriminant window.** SPARC at per-point grade cannot vote on the
particle term until the bump's source is identified or subtracted. The scatter x-shape
is not a single exponential: deep-concentrated excess + bump + floor (convergent with
the 4W/5A/5B autopsy: deep arm vertical-contested, bump point-level, environment
bounded ≲ 2–5%).

P6 status: LIVE, unkilled and unconfirmed — status-flipped in PREDICTIONS.md with this
stage as its first test; both kill directions stay pre-committed for the hardened
instrument. Unlocks named: (a) bump-source identification/subtraction; (b) hier-hardened
γ (vertical + M/L + measured ambients on distance-anchored data — with the honest 4W
expectation that SPARC depth may not resolve it); (c) DR4-era samples.

Deliverable alongside: **PREDICTIONS.md, the signed prediction ledger** (registered-
before-test rule, status-flips only, never edit; P1–P8 with uniqueness classes —
DM-immune / classical-immune / MOND-generic — and numeric kill conditions; section A =
the record, incl. the anti-classical votes already banked: the Boltzmann tail + c₁ = ½
exclude the power-law classical bath; section D = what we do NOT claim, incl. that
thermal bunching alone is classical — the quantum claim rides only on the particle
term).

PLAIN VERDICT: NEEDS REFINEMENT — the instrument is built and calibrated (that part is
SUCCESS-grade and new), the shot bath is dead (a real exclusion), but the headline
question is UNRESOLVED at this grade for a now-understood structural reason. No
credence move (no bar fired).

## Stages 7B/7C (2026-07-24): THE BUMP HUNT — caught in the inner disk; the Einstein γ closes at SPARC grade

Pre-reg d82cc4b (7B matrix) and 59c582a (7C confirmation; the 7B inner-disk finding
was post-hoc → its own pre-registered stage per house rule). Scout (Haiku, logged in
commit): transition-localized RAR scatter excess NOT FOUND at scout level — the 4T/4W
bump measurement itself appears unpublished; non-circular streaming documented at
10–40 km/s at 1–3 R_d (the mundane driver class is well-motivated).

**7B ([calcs/stage7b_bumphunt.py](calcs/stage7b_bumphunt.py)): pre-registered verdict
MIXED/UNRESOLVED — correctly, because the data broke the matrix's assumptions
informatively.** All gates PASS (G0 regression exact; G2a/b/c calibrated splits).
Discoveries in the descriptive cells:
- **Within-curve coherence MEASURED: lag-1 ρ = +0.876 (window) / +0.844 (controls),
  perm p = 0.0000, null ±0.10.** Smooth per-galaxy radial misfit dominates raw
  scatter; the 4T independence caveat now has a number; every point-level −2lnL gap
  in the program is nominally calibrated only. Offset-subtracted deciles put the
  genuinely point-level scatter at 0.04–0.07 dex (vs 0.10–0.15 raw) — convergent
  with Desmond's 0.034 hier floor.
- **The scatter is INNER-DISK-organized at fixed x**: R < 1.5 R_d cells carry ~2.4×
  the variance of outer cells in BOTH x-slices; within R/R_d strata the window-vs-C1
  excess vanishes or goes negative. The window samples 49% inner points vs 16% in
  the deep control → radius-mix candidate.
- Composition test structurally power-limited: in-window f_* terciles are
  0.986/0.998 — everything at x ≈ 1 is star-dominated; D bounds nothing.

**7C ([calcs/stage7c_gammaclean.py](calcs/stage7c_gammaclean.py)): the confirmation +
the unblocked contest.**
- **C2 (locality): DECISIVE — b_clean = 0.0000**: on outer-disk points (n = 1932) the
  explicit bump component fits to exactly zero (7A full-data: 0.083). The bump lives
  entirely in inner-disk points.
- **C1 (accounting): EXPLAINED = 0.67 — PARTIAL** (bar 0.75, miss disclosed): the mix
  explains two-thirds; the residual third is x-dependence of the inner-disk excess
  itself (inner points at x ≈ 1 are noisier than inner points elsewhere).
- **C3 (the γ contest, unblocked): UNRESOLVED by the gates, and rightly.** Variant B
  (full data + shared inner term, ĉ_in = 0.105 ≈ the 7B-predicted 0.114; the term is
  worth ~83 lnL for ONE param — standing error-model fact): profile DE-BIMODALIZES to
  interior γ̂ = 0.52, γ=2 nominally +12.1, γ=1 +7.3 — but **G2 FAILS on the informed
  design** (recoveries 0.7–6.0, no separation; the cin–floor–γ degeneracy eats the
  estimator). Variant A (clean subset): OPPOSITE lean (EQ−EC = +4.48, edge-running).
  G5 thinning: +0.14 = nothing. Two variants disagreeing in sign at few-lnL under
  measured ρ ≈ 0.87 = systematics, not physics. NO VERDICT QUOTED.

**THE CLOSE: the bump is caught — inner-disk astrophysics (streaming/beam-smearing/
decomposition class), not the law; and the Einstein particle term is NOT MEASURABLE
on SPARC** — the discriminant sits below the correlated-systematics floor at every
grade tried (raw, informed, clean-subset, thinned). P6 status-flipped: fully
out-of-sample (anchored distances / DR4-era / IFU with explicit non-circular
modeling). Neither pre-committed kill direction fired; the quantum-statistical
reading is UNTESTED in the scatter channel, not supported and not struck.

Ledger: gal-curve-coherence + gal-bump-inner + gal-gamma-final added;
gal-scatter-gamma flipped SUPERSEDED → gal-gamma-final. The 4T/4U/4W "N ~ 20–60
constraint" language stands (already hedged as constraint-not-detection) with the ρ
annotation now attached program-wide.

PLAIN VERDICTS: the hunt — SUCCESS (the 4W point-level mystery is resolved in
locality and driver class, and the program gained two standing methodological
numbers: ρ ≈ 0.87 and the 83-lnL inner term). The in-sample Einstein test — CLOSED
HONESTLY, NEEDS DIFFERENT DATA (P6 out-of-sample, pre-registered both ways).

## Stage 7E (2026-07-24): the platform translation — the lending gate is a buildable cQED experiment

Pre-reg ff1de37 ([calcs/stage7e_platform.py](calcs/stage7e_platform.py)); the queued
lab-platform stage (user directive: concrete platforms, real numbers). HONEST FRAME
FIRST: nothing here tests gravity — our own 6N cancellation theorem says flat-bath
configs measure zero, and the sky's tail-only β is anti-standard. What IS testable:
the MECHANISM CLASS — the 6X lending gate and the 6N resolution crossover. A clean
lab negative WOULD strike the 6X mechanics.

Platform card (3D circuit QED, all standard): transmon χ/2π = 250 MHz (T1 = 200 μs),
3D cavity κ/2π = 5 kHz, parametric beam-splitter λ/2π = 2 MHz, injected-noise thermal
n̄ = 0.25–8 (photon-number-splitting calibration); L=1 config ω_b = ω_a + χ; L=2
config ω_b = ω_a + 1.5χ (single-quantum channels detuned χ/2 ≫ λ). Protocol:
thermalize (~160 μs) → reset + π-pulse to |1⟩ → pump beam-splitter 2 μs → dispersive
readout of P(|2⟩) (or P(|3⟩) in the L=2 config).

ALL GATES PASS FIRST RUN (data/stage7e_platform.txt):
- GP1: closed-system weight at platform parameters = ½·n̄/(1+n̄) to 0.02–0.25%
  across the n̄ grid (tighter than 6X's band — smaller λ/χ).
- GP2: dissipative Lindblad at the real rates (cavity up/down + Γ1 + Γφ, sparse
  expm_multiply): shift +1.5% (n̄=1) / +1.1% (n̄=4) — the experiment survives
  decoherence with an order of magnitude of margin.
- GP3a: saturation discriminator — linear response extrapolates to 3.20 at n̄=8,
  the lending law measures 0.443: ×7.2 separation, unmissable at %-level readout.
- GP3b: **the L=2 geometric rung is a measurable curve** — level-3 weight regresses
  on [n̄/(1+n̄)]² with slope +0.987, rms 0.0077 (raw-n rms 0.264 = 34× worse).
  **NEW OBSERVATION: the L=2 prefactor = 0.480 ≈ ½ — the equal-time-sharing ½
  (the 6X ceiling candidate reading) persists at second order.** Logged as a second
  data point for that reading, still reading-grade.
- GP4: cQED is THE platform (ions: χ_eff/g ~ 10 marginal; optomech: no Kerr). The 6N
  resolution crossover is switchable on the same chip (broadband 50Ω ↔ Purcell-filter
  port; χ/bandwidth spans ~0.1–250).

PLAIN VERDICT: SUCCESS — the mechanism's building blocks (the lending law ½·P(n≥1),
its saturation, the s² geometric rung) are concrete, currently-buildable measurements
with quantified error budgets; framed honestly as mechanism-class validation, not a
gravity test.

## Stage 7D (2026-07-24): the VACUUM-SHARE contest — the galaxies demand the "+1" at kill grade ×22; the binary direct leg is instrument-limited

Pre-reg c53f2e2 ([calcs/stage7d_vacuumshare.py](calcs/stage7d_vacuumshare.py)); three
gate amendments logged pre-results (non-uniform g→0 limit; its O(g/(ν−1))
amplification; the resolution cliff → NR=2048 + graded identity gate) + one caught
stale-table hazard (run-1's g1p2 table written before its own gate fired, silently
reused by the resume check — purged, exists-branch now prints). All gates green on
the final run: GQ1 2.78e-10 with linear-in-g ratio 10.0; GQ2 5e-15; BE regressions
d = −0.00 both treatments; AMB node reproduces 6F to the hundredth (−59.05).

THE QUESTION: the AMB local factor q = 1/(2ν−1) is the vacuum share of the quantum
(2n+1) pull. A classical bath's pull carries 2n only — no zero-point. QCL = the same
admixture grammar on the classical share q_cl = 1/(2(ν−1)); RJA = classical-
equipartition ambient occupation (n_RJ = 1/x). 6U never contested the local factor.

RESULTS ([data/stage7d_vacuumshare.txt](data/stage7d_vacuumshare.txt)):
- **QCL galaxies: ANNIHILATED — Δ(QCL−AMB) = +556.51 vertical (bar: reject ≥ +10,
  kill ≥ +25 → exceeded 22×); cap16 = +580.06 (looser cap = WORSE ⇒ the rejection is
  a LOWER bound, not a cap artifact); plain +274.06; worse than pure BE by +497/+521/
  +182 — adding the classical admixture is worse than no admixture at all** (the R0
  structure, now measured directly). Mechanism visible in the numbers: without the
  vacuum floor, q_cl = 1/(2n) grows as e^x in the tail → β_cl diverges → screening
  sharpens toward a step (p_eff = 4.5 at cap 8) — and the RAR tail, which measured
  p ≈ 0.65–0.75, destroys it. **The quantum +1 is load-bearing: it caps q ≤ 1 and
  keeps the admixture finite exactly where the data demand finiteness.** The fit's
  compensation attempts (f → 1.456/1.725) fail.
- **QCL binaries: INSTRUMENT-N/A** — the e_N = 1.2 table failed the spherical
  identity at 10.61% even at NR = 2048 (the capped function's cliff, ν−1 ~
  e^(−y^4.5), defeats the multipole solver); the pre-registered ≥5% abort fired and
  the leg was skipped rather than fit on an untrustworthy table. Not a contrary
  vote; a solver limitation, disclosed. (Indirect binary support stands: binaries
  ACCEPT the quantum-share AMB (6G) and bound β < 0.03 (5R).)
- **RJA: null/unresolved both, as pre-stated for the weak rung** — galaxy +0.24
  (deep-ambient BE→RJ convergence, near-null by construction); binaries mean +1.68
  vs AMB with per-seed scatter −2.7…+10.3 (4 seeds; the +10 is one realization —
  the 4L-b lesson) = UNRESOLVED; α̂ = 1.070 ± 0.016 interior 4/4. The
  ambient-statistics digit is not readable at current precision.

FORMAL VERDICT per the pre-registered bars: **PARTIAL — "vacuum rung measured in
both systems" is NOT claimable** (the binary direct leg never ran). The claimable
statement: **the galaxy ladder demands the vacuum "+1" in the local factor at
kill-grade ×22, cap-robust, both treatments — conditional on the admixture grammar
(like every 6U rival exclusion)** — sitting alongside the banked rungs R0 (β=0 pure
BE excluded by the gate credit −59/−100) and R1 (c₁ = 0 excluded in both systems,
4S/4X). Quantum-ladder status: occupation rung (banked), vacuum-share rung (galaxy
direct + binary indirect), ambient-statistics rung (open), fluctuation rung (P6,
out-of-sample).

PLAIN VERDICTS: the galaxy vacuum-share leg — SUCCESS (a kill-grade exclusion,
measured). The two-system headline — NEEDS REFINEMENT (binary leg needs either a
cliff-capable 1D solver or a softer classical variant; queued as a refinement, not
promised). RJA — UNRESOLVED (weak rung stayed weak, honestly).

## Stage 7F (2026-07-24): the mixing-mean uniqueness contest — asymptotic uniqueness proven; window bound p_mix ∈ [−0.5, +0.1); the 6H c₄ cross-validated

Pre-reg a4696c1 ([calcs/stage7f_mixmean.py](calcs/stage7f_mixmean.py)); one
pre-results estimator amendment (far-tail underflow → index read stably from the
solution's own argument u; logged in-script). Seam (iv) of the derivation — WHY the
geometric mixing between the source argument A = √y and the total/boot argument
B = νy — contested against the full power-mean family u_p = [(1−β)A^p + βB^p]^(1/p)
with the AMB running gate. All gates PASS (data/stage7f_mixmean.txt):

- **GF1, the exact endpoint theorem (sympy, symbolic β)**: p > 0 snaps the tail
  exponent to 1 (boot grade — binary-vetoed program-wide); p < 0 snaps to ½ (ungated
  BE — the measured two-system tail SPLIT becomes impossible); **only p = 0
  interpolates continuously**, p_tail = (1+β)/2 → ½ + g/4 = the measured pair.
  Asymptotic uniqueness of the geometric mean: PROVEN.
- **GF-deep (symbolic p, order-by-order)**: a₁ = ½, a₂ = 1/12, a₃ = −g/16 all
  p-INDEPENDENT; **a₄ = −gp/64 + g/192 − 1/720** — the mixing choice first enters at
  the c₄ rung (da₄/dp = −g/64). Two dividends: (i) WHY the Bernoulli ladder never saw
  the mixing choice (it is p-blind through c₃ under the running gate — the tail must
  do the selecting); (ii) **at p = 0 the a₄ value reproduces 6H's independently
  derived c₄(L=2) = s²/192 − 1/720 EXACTLY** — two derivations, different routes,
  same rung: strong cross-validation.
- GF2/GF3/GF4 regressions: p=0 vs the 6G Newton solver 2.5e-14; BE/gm/boot members
  0–2e-16; far-tail index regresses to ½+g/4 at 4 digits (0.6885/0.5280 = the known
  postdictions).
- **The windowed bound (y ∈ [3,30] vs the measured bands)**: positive flank SHARPLY
  excluded — p ≥ +0.1 exits the galaxy band (the boot direction dies immediately);
  negative flank bounded at −0.5 (galaxy exits below); **both-in-band range
  [−0.50, 0.00]; p_mix = 0 comfortably inside both systems.**

FORMAL VERDICT per the pre-registered (symmetric) bar: **PARTIAL** — my bar demanded
|p| ≥ 0.5 excluded on BOTH signs; the negative flank retains window slack (a p < 0
member with |p| ≲ 0.5 shows geometric-grade in-window behavior because its snap to ½
lies beyond the measured window). The honest content: **asymptotic uniqueness proven
+ an asymmetric measured bound p_mix ∈ [−0.5, +0.1)** — the boot-ward side is
data-excluded, the BE-ward side is scale-degenerate today (choosing p < 0 requires
the observation window to sit before the snap — a coincidence-of-scales objection,
stated but soft). Seam (iv) upgraded: from "chosen algebra" to "measured index with
an exact uniqueness theorem at asymptotic grade"; fully closing it = a tail
measurement beyond y ~ 30 (deep-window galaxies / DR4-era) or the c₄ rung.

PLAIN VERDICT: SUCCESS on the theorem + the cross-validation + the positive-flank
exclusion; NEEDS REFINEMENT on the negative flank (window-limited, stated).

## Stage 7G (2026-07-24): the trajectory formulation vs Saturn — the cancellation QUANTIFIED (451 orders); the binary equivalence PARTIAL at proxy grade

Pre-reg 8ce8466 ([calcs/stage7g_trajsaturn.py](calcs/stage7g_trajsaturn.py)); the 6Y
Saturn corollary converted to numbers. Machinery: the 4L mi_t prescription verbatim
(per-orbit boost at the Kepler time-average, EFE-respecting via the gated AMB table,
exact-Newton engine) on the corrected-velocity v7 pipeline; 6G MG-AMB rows as
same-seed comparators; G0 α=0 bit-identity PASS (Newton rows match to 0.07).

**Q2 — SATURN, the headline ([data/stage7g_trajsaturn.txt](data/stage7g_trajsaturn.txt)):
in the trajectory formulation the anomaly attaches to the worldline's own occupation.
Saturn: y = 5.38e5 → occupation argument u = 1060.8 → trajectory anomaly ν−1 ~
10^(−460.7) vs the Cassini-equivalent δa/a = 6.7e-11 ⇒ THE TRAJECTORY OBSERVABLE
SITS 451 ORDERS OF MAGNITUDE UNDER CASSINI (bar: ≥2). B2 PASS.** The transition
shell that sources the 4K quadrupole never exists on Saturn's trajectory; the
residual observable is the ordinary galactic tide (7.5e-31 s⁻² = 4026× under the
bound, the 6W moot block). Dividends: (i) trajectories crossing r_M = 7030 AU (Oort
bodies, extreme-aphelion comets) DO sample the transition = the residual solar
probe; (ii) the CONTRAST row (order-of-magnitude, labeled): a power-law tail leaves
an r-dependent G_eff at ν−1 ~ 1.9e-6 at Saturn — ephemeris-dead by orders — **the
trajectory door is open ONLY for Boltzmann-screened functions, i.e. exactly the
screening class the data selected. The two halves of the program point at each
other.**

**Q1 — the binary equivalence: PARTIAL, honestly.** mi_t-AMB vs MG-AMB same-seed:
−11.33 (seed 31, α̂ = 1.55 interior) and −5.54 (seed 101, α̂ = 2.00 EDGE); mean
−8.43; interior 1/2. Per the pre-registered bars this is neither TIE (needed |D| ≤
5/seed + interior both) nor formally NARROWS (needed mean ≤ −10): **the trajectory
PROXY trails the field formulation by ~8.4 ± ~3 at 2 seeds** — in the same band as
the eight-function sharpened penalty, and with α̂ inflated (1.55/2.00 vs MG's 1.06:
the per-orbit averaging washes out boost and buys it back with amplitude —
consistent with proxy-discretization loss, but not provably so). Cannot distinguish
"the proxy is too crude for the AMB shape" from "the class genuinely fits worse"
at this grade. The 4L STANDARD-law tie (mi_t ties MG, −3.5±3.3, 6 seeds) does NOT
automatically transfer to AMB. Refinement (unpromised): budget-grade 6 seeds + a
better-than-proxy trajectory implementation. No post-hoc seed extension — the
pre-registration fixed 31/101.

Carried costs (stated in-output): mi_t is a per-orbit proxy, not a full nonlocal MI
theory (time-nonlocality + conservation need the real formalism); Petersen & Lelli
2020 constrain specific MI models on rotation curves; **the MG-formulation
quadrupole tension STANDS in the world table — this stage opened a door by numbers,
it did not delete the row.**

PLAIN VERDICTS: the Saturn cancellation — SUCCESS (quantified at 451 orders; the
4K paradox now has a numerical open door within the derived class). The binary
formulation-equivalence — NEEDS REFINEMENT (proxy trails −8.4 ± ~3, one α̂ edge).
Overall: "modified inertia might do it" has become "the derived coupling class
does it on Saturn's side by 10^451, at a measured, unresolved −8-lnL binary cost
for the proxy."

## Stage 7H (2026-07-24): TIE-RESTORED — 7G's gap was the proxy's Jensen discretization; the trajectory formulation fits the binaries at field grade

Pre-reg e9db7c8 ([calcs/stage7h_miavg.py](calcs/stage7h_miavg.py)); one pre-results
gate amendment (G2 quadrature bar 1e-6 → 1e-5 at the e = 0.999 extreme — measured
7.3e-6, three orders below likelihood resolution; disclosed). The SIGNED HYPOTHESIS
was committed before the run: mi_t evaluates B(⟨g⟩); eccentric orbits linger at
apoapsis (low g, high boost), so the honest adiabatic functional is ⟨B(g)⟩ —
strictly boostier by Jensen (+0.078 at y=1, e=0.9, measured in the gate block) —
predicting boost recovery, α̂ deflation from 1.55/2.00 toward 1.06, and gap closure.

RESULT ([data/stage7h_miavg.txt](data/stage7h_miavg.txt)): the 2-seed result hit the
pre-registered ambiguity band → auto-extended to six. **Six-seed budget: D(mi_avg-AMB
− MG-AMB) = −14.41 / −5.32 / +7.39 / −7.06 / −4.94 / +6.81 → mean −2.92 ± 3.46 SE,
SIGN-MIXED (2/6 positive), interior 6/6, α̂ = 1.272 ± 0.038** (edge-riding gone;
deflation confirmed), G0 α=0 bit-identity PASS. The pre-registered TIE bar (mean ∈
[−5,+5], all interior) FIRED. Same grade as the 4L standard-law tie (−3.5 ± 3.3).

**All three signed predictions landed: boost recovered, α̂ deflated, gap closed.
7G's −8.4 was the proxy's discretization, not a class cost.**

THE COMBINED SATURN STATEMENT (7G + 7H, the strongest the program has ever made):
the trajectory formulation of the derived ambient-gated function **fits the wide
binaries exactly as well as the field formulation (−2.9 ± 3.5, six realizations,
all amplitudes interior) AND passes Cassini by 451 orders of magnitude** — the
quadrupole tension is a property of the FIELD formulation, not of the measured
function; and the door is tail-selective (Boltzmann-screened functions only = the
measured screening class). Honest residuals carried: α̂_traj = 1.27 vs the field
1.06 (~20% amplitude recalibration between formulations — a nuisance-level
prescription dependence, noted, not hidden); the adiabatic average is still an MI
BRACKET, not a full nonlocal theory (time-nonlocality/conservation live in the
real formalism); the MG-formulation tension STANDS in the world table for the
field reading; DR4 eccentricity resolution remains the formulation kill test.

PLAIN VERDICT: SUCCESS — prong (b) of the maximal demand is complete at its
achievable grade: equal binary fit + 10^451 Saturn margin, reached by
pre-registered stages with a signed mechanistic hypothesis confirmed.

## Stage 7I (2026-07-25): the external-review round — the freeze, the robustness rows, the census release, and a MATERIAL branch worked to ground

CONTEXT. Two external reviews solicited by the user (an LLM review 2026-07-24;
an Opus 5.0 review + follow-up exchange 2026-07-25); point-by-point responses
in the session record; the adopted items became this stage and the language
pass (commits 2aa6ef9, a17d3a5, 3ef57c7 — each BEFORE its execution).
Conceded in full from review 2: (i) the ambient gate's two-system pass is a
DESIGN PROPERTY, not a measurement (the gate places AMB within 0.3% of the
occupation law at the binary ambient by construction — our own 6G pre-reg
called that leg "near-guaranteed"); reframed at all four claim sites. (ii)
Our "pre-registration" is BAR-LOCKING (same-session committed thresholds;
guards bar-moving and outcome-shopping, NOT the sequential-search error
rate) — defined at first use in PAPER §1 and PREDICTIONS.md rules. (iii) The
~30-form search multiplicity is uncorrected and uncorrectable post hoc — the
FREEZE + out-of-sample routing is the repair (PREDICTIONS.md §0: search
closed at the world-table thirteen; new forms = consistency rows only;
primary analyses declared: vertical-hardened hier / six-seed
corrected-velocity v7; scatter channel restricted to outer/anchored;
mode-count N retired). (iv) Binary c1 = 0.37-0.50 regraded to
TWO-REALIZATION INDICATION at all six claim sites (the thinnest headline
number in the paper; six-seed budget queued). Also adopted: section 2.4
stays out of anything circulated (research log until a Lagrangian); the
split (Paper 1 binaries+ablation, Paper 2 coefficients+quadrupole) with
STRUCTURAL INDEPENDENCE from 2.4 as a standing drafting rule. Pushbacks
that landed: 3K's kinematic cap is threshold-free (conceded by the
reviewer, with the right consequence — 7J must MARGINALIZE f_comp under a
completeness prior, not re-fence); the ceiling band pairs are
astrometrically CLEANER than the same-vt low-gamma continuum (RUWE 1.06 vs
1.28 — wrong direction for an error-tail identity; granted in full); F&B
2005 has no binary leg (granted; the galaxy-side c1 lean = the field's old
simple-over-standard preference promoted to coefficient language; F&B
cited). Correction in our own mouth, adopted from the exchange: the
"termination at the boosted edge" is non-contradiction (P = 0.62), not
positive evidence — the band count carries the census, which is where the
tail objection bites (-> 7K).

PART C — THE CEILING-PAIR RELEASE (GATE GC PASS). data/ceiling_pairs.csv
committed (gitignore exception): 23 rows — the 11-pair census + below-edge
context + above-cliff; per pair: source IDs, both velocity conventions, gamma
both ways, systemic RV + perspective exposure, direction S/N, sigma_vt, RUWE
both, R_chance, |dplx|/sigma. 4J reproduced 11/11 (both S/N conventions);
corrected recount 9 (4Q's sqrt2-edge construction gave 10->9; the +-1 edge
sensitivity was already logged). Re-countable under any convention without
our pipeline.

PART W — w_rad FROZEN TO THE EXTERNAL HWANG 0.21 (BW: CLOSED). d_alpha =
-0.010+-0.000 (simple) / +0.005+-0.005 (BE); Newton within 1.4 of baseline;
interior 2/2. The "eccentricity mixture fitted jointly with the boost"
objection is closed with the reviewer's own instrument (free fit had chosen
0.20 in 12/12; Hwang-implied 20-22%).

PART S — STRICT MULTIPLICITY CUT (BS: MATERIAL, extension fired, 6 seeds).
RUWE < 1.2 both + no overluminous component (3J delta < -0.4 ridge; star
fraction reproduced 0.123 = GATE GS PASS) retains 8,047/14,071 = 57.2%. The
v7 fits COLLAPSE: d_alpha = -0.740+-0.053 (simple) / -0.457+-0.048 (BE),
Newton +9.4/+12.5 = 16/24% of the N-scaled expectation, interior 6/6, w_rad
at the GRID EDGE 0.3 in 6/6 (vs 0.2 interior 12/12 full-sample).

THE INVESTIGATION (three bar-locked instruments; both rescue hypotheses
died on their own bars):
1. Model-light discriminator (stage7i_material.py): the 2C anchor-relative
   median boost — no forward model — SURVIVES and slightly strengthens
   under every cut component: raw 1.086 -> 1.103 [boot 1.071-1.142];
   overluminous-free alone 1.099; corrected 1.078 -> 1.099; perpendicular
   (mass-immune) 1.151 -> 1.185. Companions-as-the-median-carrier REFUTED:
   removing the flagged 43% RAISES the ratio.
2. FPM profile (stage7i_fpm.py, bars at a17d3a5): the named suspect —
   fixed noise inflation 1.5 over-broadening the clean subsample — REFUTED
   under its own locked bar: the strict sample PREFERS 1.5 in 4/4 rows
   (+2..+7 lnL over 1.2) and a_hat stays collapsed at every FPM. My miss,
   caught by the bar.
3. SW (stage7i_sw.py, bars at 3ef57c7): the surviving flag — the
   alpha<->w_rad degeneracy on thinned gamma statistics (wide bins
   776/131 pairs) — REFUTED: strict + w_rad frozen at the external 0.21
   gives a_hat 0.55 (0.46/0.64) / 0.65 (0.60/0.70), Newton +9.9/+13.4,
   interior 4/4 ==> COMPANION-DIRECTION per the pre-committed tree.

FACT CHECK ON THE CENSUS: 5/11 forbidden-band pairs survive the strict cut
vs 6.3 expected at 57% retention — base-rate, no multiplicity enrichment of
the band; the census stays a full-sample statement (a 5-pair strict band is
below census grade on its own).

END STATE (stated, not resolved): at full-likelihood level the cleanest 57%
of the catalog genuinely prefers ~half the boost amplitude and a Newton
margin of only +10-14, and no instrument correction we locked in advance
explains it away — while on the SAME cleaned pairs the model-light median
boost is 1.10, the perpendicular boost 1.19, and the band survives at base
rate. The model-free statistics see the boost where the forward model's
likelihood does not; both cannot be right. This is the reviews' companion
concern given empirical teeth at fit level — disclosed in PAPER 6.3 in
exactly these terms. DECISIVE INSTRUMENT = 7J: photometric completeness
forward model -> prior on the TRUE multiplicity fraction; v7 refit with
(f_comp, w_rad, f_pm, alpha) freed jointly; report the MARGINALIZED
(f_comp, alpha) posterior (the Opus design, adopted before this result made
it urgent). CREDENCE MOVED (honest pricing): anomaly-real 70% -> ~60-65%
pending 7J — the fit-level wound is real; the floor under it is that the
median cannot be companion-carried (instrument 1) and the galaxy legs are
untouched.

PIPELINE RULES ADOPTED: noise inflation and nuisance grids are sample
properties — any subsample row must profile FPM and watch grid-edge riding
(w_rad); every future robustness row inherits the three-instrument pattern
(model-light discriminator first).

PLAIN VERDICTS: freeze + census release + language pass: SUCCESS. Part W:
SUCCESS (the w_rad objection closed at d_alpha <= 0.01). Part S +
investigation: NEEDS REFINEMENT, honestly disclosed — three locked
instruments, three hypotheses dead (two of them mine), the tension stated
and handed to 7J. The measurement-grade statement that stands today: the
median boost is stable at 1.09-1.10 under the strictest multiplicity
cleaning this catalog supports.

QUEUE: 7J completeness + marginalized (f_comp, alpha) posterior [DECISIVE];
7K ceiling tail control (empirical error tail from gamma>=75 anchor pairs,
per-slice); 7L v7 on the Cookson selection carrying the ablation ladder
(after 7J); c1 six-seed budget; THE SPLIT (Paper 1 / Paper 2; 2.4 stays in
the log; structural-independence rule).

## Stage 7J (2026-07-25): the decisive instrument fired — COMPANION-WIN at the measured multiplicity; CORRECTION #15; the binary data vote suspended pending 7K

Pre-registered d2dc7eb; amendments c893827 (Part A restructure + envelope),
a2ecf2e (host-axis remap + photocenter law — Opus's units point, code-
verified), a9e8b16 (provenance disclosure, raw conditional record, ell(q)
validation leg, fpm extension to 2.4, single-seed rule). Scripts:
[calcs/stage7j_completeness.py], [calcs/stage7j_marginal.py],
[calcs/stage7j_diag.py]; outputs data/stage7j_completeness.txt,
data/stage7j_{full,strict}_photo.txt, data/stage7j_verdict.txt,
data/stage7j_diag.txt; cubes data/stage7j_cube_*.npy (raw + photo).

PART A — THE COMPLETENESS MEASUREMENT (the program's first): binned
Poisson mixture fit of the overluminosity delta distribution on the
sample's own MS locus, flux-weighted blend photometry, q^-0.5 systematic
folded into the envelope. Results: flag completeness C = 0.410
(false-flag 0.023); f_photo in [0.22, 0.30]; host-axis remap through
P_blend = 0.53-0.59 gives f_host in [0.42, 0.57] (peak 0.51 — consistent
with Raghavan 2010 ~46% solar-type multiplicity = external validation);
strict residual r_host in [0.37, 0.53]. Reference prices: the old fence
f_comp <= 0.1 sits at ln pi = -768 (full) / -608 (strict); Banik's free-fit
0.69 at -36. Gates GA1-GA3 PASS.

THE RAW CONDITIONAL (as-published amplitude law, q/(1+q) wobble): at the
measured host fraction the full-sample fit returns alpha = 0 with Newton
tying, at a likelihood cost of -486..-525 (fcomp 0.35) / -1034..-1078
(0.5) / -1967..-2010 (0.7) relative to its own (0.1, alpha~1.1) optimum
— the as-published law is INTERNALLY INCONSISTENT with the measured
multiplicity (it cannot host the companions the sky demands). This is
why 3K's physical-multiplicity fit capped f_comp at 0.1: the law it
tested overstates high-q wobble (no photocenter cancellation — twins
wobble zero; the astrometric amplitude is |q/(1+q) - ell/(1+ell)|).
Provenance carried honestly: the correction was registered before the
conditional was READ but not provably before it existed on disk.

PART B — THE PHOTO-MODE VERDICTS (photocenter-corrected amplitudes,
kappa_w in {0.7,1.0,1.4}, fpm extended to 2.4 per correction-#4, GB0p
cross-mode gate 0.00e+00 exact 8/8, host prior on axis fcomp):

  full   (2 seeds): a_marg = 0.00, dN_marg = +0.0, both laws
                    ==> COMPANION-WIN (bar: <=0.7 or <=+15 either)
  strict (2 seeds): a_marg = 0.00 ==> COMPANION-CONFIRMED — but
                    FORMALLY WEIGHTLESS: the power gate fired POWER-FAIL
  strictpow: injected truth alpha=1.18 (simple) + r_host companions;
                    simple arm RECOVERED a_marg=1.27, dN=+48.8 (in-band);
                    BE arm rode to the 2.00 grid edge -> the both-laws
                    bar failed. Honest reading: the instrument HAS power
                    (simple-arm recovery through the same prior), the BE
                    edge is a cross-family recalibration artifact at 1
                    seed; formally the strict quadrant carries no weight
                    and the verdict rests on the FULL sample alone.

7J-c DIAGNOSTICS (pre-registered pipeline rule "watch grid-edge riding";
[calcs/stage7j_diag.py], no refitting): (R1) the COMPANION-WIN is NOT
carried by the fpm edge — restricting to fpm <= 1.8 or fpm = 1.5 leaves
a_marg = 0.00, dN_marg = +0.0 and the posterior cell unchanged in 8/8
reads (the PROF rows do ride fpm to 2.4, but the marginal doesn't need
it). (R2) mechanism: the likelihood alone still prefers (alpha ~ 1,
fcomp ~ 0.1-0.2) — the prior prices that at -768/-268 — and at the
prior-allowed fcomp (posterior 100% at 0.35 full / 0.20 strict, kw=0.7)
the fitted alpha is 0. The posterior chooses "Newton + measured
companions, kinematic misfit -60..-116" over "alpha ~ 1 + forbidden
multiplicity". alpha = 0.5 / 1.0 are disfavored by 17-61 / 38-87 lnL in
the marginal. (R3) the winning cell is itself a POOR ABSOLUTE FIT: under
this forward model NO cell fits both the kinematics and the measured
multiplicity; the photocenter correction softened the multiplicity cost
~5x without closing it. Which channel carries the residual -60..-116
(the median offset? the gamma structure?) is not derivable from the
cubes — that is 7K-a. (R4) fcomp = 0 is REJECTED kinematically
(-660..-811 full, -111..-146 strict): the data demand companions; they
demand FEWER than the sky hosts; at the measured rate the boost dies.
(R5) kw pins to 0.7 everywhere — soft edge, second-order.

CORRECTION #15 (retraction of an inference, not a value): 7I
instrument-1 concluded "companions-as-median-carrier refuted" from the
median's survival of the strict cut. At measured C = 0.410 that argument
is UNDER-POWERED: the cut removes only flagged companions AND shrinks
the denominator, so the companion FRACTION barely moves (0.51 -> ~0.47)
— a fully companion-carried median would survive the cut within CI. The
median VALUES stand (1.078/1.086/1.151, model-light facts); the immunity
inference does not. Ledger row ret-7i-median-immune; the ceiling census
inherits the same audit (its leakage null predates the completeness
measurement) -> 7K-b.

SYNTHESIS (quadrant C fired, pre-committed): credence anomaly-real
~60-65% -> ~35% (low end of the pre-committed 35-45 band, because
correction #15 removes the "floor" argument the band was written with;
not lower, because the winning cell's own -60..-116 misfit means the
forward model still cannot actually explain these data, and the galaxy
legs are untouched). The binary DATA VOTE in the world table is
SUSPENDED (banner added): every BIN/A0 row was measured inside the
fence; all binfn-* function contests are fence-conditional pending
re-run under the 7J posterior (7J-d). Standing caveat BOTH ways: the
forward-modelability limit — the injection recovery validates the
instrument against ITS OWN companion model; a real companion population
differing in q-distribution / inner periods could drag alpha to 0
through mismatch the self-consistent gate cannot expose. And the same
limit cuts the other way: it also cannot manufacture a boost.

WHAT 7J CHANGES STRUCTURALLY: the anomaly's binary evidence now rests
entirely on the MODEL-LIGHT channels, and both need re-derivation under
measured multiplicity: 7K-a = forward 2C median at the 7J winning cell
(alpha=0, fcomp=0.35, kw=0.7, photo wobble) vs the observed 1.078 (CI
1.052-1.103) — if Newton+measured-companions reproduces the median, the
median falls as evidence; if not, it stands as the unexplained residual
R3 points at. 7K-b = ceiling-census leakage null under the same cell —
11 forbidden-band pairs with the cliff structure, or not. These are now
THE decisive binary instruments (the tail-calibration design from the
review folds into 7K-b). 7L (Cookson selection) unchanged. The galaxy
program (coefficients, screening index, scatter, ambient gate) does not
touch any of this and proceeds on its own grades.

LEDGER: rows bin-7j-completeness, bin-7j-marginal, bin-7j-rawcond,
ret-7i-median-immune added; bin-alpha-final, bin-newton-final,
bin-newton-v7corr -> CO-QUOTED (fence-conditional, pointer
bin-7j-marginal); bin-boost-corr, bin-boost-perp -> CO-QUOTED (pointer
bin-7j-completeness); bin-ceiling note gains the 7K-b audit; worldtable
+6 tokens (44 total) + the conditionality banner; all six gates PASS.

PLAIN VERDICTS: Part A (completeness measurement): SUCCESS — the
program's first measured companion function, externally validated.
Raw conditional: SUCCESS as an instrument finding (the as-published law
is internally inconsistent; the field's wide-binary companion modeling
inherits this). Part B (the decisive question "does the boost survive
the measured multiplicity at fit level?"): DIFFERENT PHYSICS — the
pre-registered answer is NO: at the measured companion rate the forward
likelihood prefers Newton + companions, and this is not an fpm-edge
artifact. The anomaly is not dead — the winning model misfits its own
data by -60..-116 and the model-light channels are unadjudicated — but
the burden has moved to 7K, and the credence moved with it (~35%).

QUEUE: 7K-a forward median at the winning cell + 7K-b census null under
measured multiplicity [BOTH DECISIVE]; 7J-d function-contest re-run
under the 7J posterior; 7L Cookson selection; c1 six-seed budget; THE
SPLIT after 7K (Paper 1's binary chapter now leads with the marginalized
null + the completeness measurement).

## Stage 7J-x (2026-07-25): the second-review hardening round — per-arm
## injections, the low-end prior stress, CORRECTION #16, and the wrong-
## population citation withdrawn

The second solicited review (Opus 5.0) audited the 7J close. Its points,
each resolved by computation or booked as a correction — not argued:

(1) POWER-GATE SCOPE (load-bearing, conceded): the strictpow injection
ran on the STRICT selection under the SIMPLE truth only — the full
sample carried the verdict without its own injection, and the BE arm
was never asked to recover its own truth anywhere (the 2.00-edge
POWER-FAIL was a cross-family read, not a BE-arm validation). Standard
adopted verbatim: any arm carrying a verdict needs its own passing
injection through the same prior on the same sample; a grid edge in a
power gate disqualifies exactly as in a fit. AMENDMENT 6 pre-registered
(commit 0045eea, BEFORE execution): modes fullpow (simple truth 1.18,
full selection, host-peak companions; bar a_marg in [0.9,1.5] AND
dN >= +25 on the simple row), fullpowbe + strictpowbe (BE truth 1.13;
bar BE a_marg in [0.85,1.45] AND dN >= +20). STANDING RULES committed:
COMPANION-WIN is an either-law bar -> stands if >= 1 full arm
validates; if BOTH full arms fail, verdict downgrades to
PROVISIONAL-INSTRUMENT and the credence move partially reverts
(~35% -> ~50%) pending repair.

(2) TWO SEEDS AT A BOUNDARY (conceded; correction-#10 is our own
precedent): agreement at alpha = 0 compresses the luck-detecting
scatter. Seeds 202/303/404/505 launched on full+strict photo; the
boundary-free check statistic is the alpha-marginal gap lm(0.5)-lm(0);
the verdict stands unless any new seed shows gap > -5 lnL or an
interior a_marg > 0 at the peak prior (either -> AMBIGUOUS per the
original tree). Pre-registered in amendment 6, same commit.

(3) LOW-END PRIOR (the reviewer's sharpest scenario — "if the prior is
20% too high the bound is where alpha goes regardless"): EXECUTED
immediately from the cached cubes ([calcs/stage7j_lowend.py],
data/stage7j_lowend.txt). Result: recentring the measured envelope at
its 1-sigma-low end (full peak 0.51 -> 0.42; strict 0.47 -> 0.37)
leaves COMPANION-WIN standing (full seed means a_marg 0.19 simple /
0.25 BE, bar <= 0.7 fires; strict unchanged) but UN-PINS the zero: the
full posterior moves to fcomp = 0.2 with a residual-boost scrap
(seed 31: a_marg 0.37/0.50, Newton within +4.9/+7.3; seed 101: 0.00) —
the exact zero is peak-prior-specific and is now co-quoted with the
low-end row everywhere. The 16-lnL seed scatter at the low end
independently supports the seed budget of (2).

(4) MISFIT SAME ORDER AS THE GAP (conceded): the winning cell's
-60..-116 misfit is the same order as the +99 Newton-rejection gap the
fenced fits had won by — the marginalization selects the less-bad of
two misspecified models, and that ordering is provisional against the
missing component. Promoted from a residuals paragraph into the same
breath as the verdict (PAPER abstract (3) + 6.3).

(5) PHOTOCENTER PROVENANCE (the reviewer asked twice; answered
plainly): the correction was registered (a2ecf2e) AFTER part A's host
inversion raised the amplitude question and AFTER an interim verdict
under the misspecified blended-axis prior had already shown the
companion direction (quarantined in
data/stage7j_verdict_stale_blendprior.txt), but BEFORE the as-published
conditional was explicitly read (a9e8b16). The physics is
sequence-independent (exact twin cancellation, validated by the ell(q)
leg); the adoption environment was NOT neutral; the as-published
conditional travels alongside everywhere as the guard. Sentence added
to PAPER App A item 15 and 6.3.

(6) CORRECTION #16 (the hard one, stated generously as demanded):
7.3's claim that the 16-sigma-Newton result was "manufactured" by a
companion fraction "the photometry forbids" is RETRACTED — the
completeness measurement shows the photometry never forbade ~0.5 (the
flagged 12% at C = 0.41 implies 0.42-0.57), our own f <= 0.1 fence was
the excluded assumption (ln pi = -768), and our marginalized verdict
landed on the Newton-favored side of the fit-level question. **Banik
et al.'s companion-dominated reading was substantially closer to the
truth than our fenced fits.** Retracted in the normative sense for ALL
THREE ablation legs (the arithmetic stands as an as-published-law map).
Still standing, each specific: their fitted 0.69 excluded (ln pi = -36;
~-2000 in our 2D likelihood); the shared as-published amplitude law
internally inconsistent without photocenter cancellation; H&C's
sub-error-binning critique (theirs); the direction channel's
information advantage. 7.3 retitled + rewritten; abstract (4), intro
Third, and the conclusions bullet rewritten hard; App A item 16.

(7) RAGHAVAN = WRONG POPULATION (conceded) — AND THE SCOUT RETURNED A
REAL EXTERNAL TENSION: 46% overall solar-type multiplicity counts the
wide companion itself; our f_host is HIGHER-ORDER multiplicity — the
subsystem fraction among components of already-selected wide pairs
(per-component 0.24-0.34). The scout (Haiku, scout-grade, primary-
source verification NOW REQUIRED per the 4-misread rule) reports the
published rates SIT A FACTOR ~2-3 BELOW US: Tokovinin 2014 (AJ 147,
86-87) 10.0%/7.3% per primary/secondary component; Tokovinin 2010
12+-4% (5-100 AU subsystems); Hwang 2022 field wide-tertiary baseline
5.35%; Moe & Di Stefano triples ~7-8% of systems. If the literature is
right, part A's completeness C = 0.41 is the suspect quantity (it
divides the hard 12.3% flagged rate; C ~ 0.8 would reconcile) and the
host prior is over-scaled. RESOLVED AT VERDICT LEVEL IMMEDIATELY
(7J-e2, same script, cached cubes): **the COMPANION-WIN does not hinge
on the prior scale — at the literature-anchored cell (fcomp = 0.2) the
full-sample conditional is alpha = 0.36/0.49 (seed 31) and 0.00 (seed
101) with Newton within +4.7/+7.3, and a literature-centred prior
(peak 0.22, sigma 0.08) gives seed-mean a_marg 0.20/0.25, dN +2.5/
+3.7 — COMPANION-WIN fires at every plausible prior anchoring
(measured peak 0.51 / measured low 0.42 / literature 0.22); what
changes is only the null's flavor (exact zero -> a one-seed scrap of
~0.4). Even in the most boost-friendly plausible cell, >= 93% of the
fenced +99 Newton rejection is gone.** Part A's ABSOLUTE scale carries
the external-tension flag until the primary-source requote (candidate
correction #17 if C is confirmed biased; the scout also reports NO
published subsystem-vs-separation curve in our 0.2-50 kAU regime =
part A may be measuring something genuinely new at these separations
— both branches stated, neither claimed). "Raghavan-consistent"
withdrawn from PAPER/LEDGER.

(8) PROPAGATIONS the reviewer demanded and got: 8.1's binary-calibrated
quadrupole novelty SUSPENDED with the binary amplitude (status headnote;
the tension reverts to the DHF24 galaxy-calibrated form; the 6W
composition contest and 8.2 MI/MG brackets carry the same flag);
the AMB "passes both systems" claim now carries the second flag (the
binary leg fence-conditional; the two-system split it reconciles is
suspended pending 7J-d) in the abstract and conclusions; PREDICTIONS.md
corrections record extended.

(9) The reviewer's Part-A-as-standalone-paper reading and the 7K
twin-check design (are the 11 band pairs' components in the
wobble-capable regime? — with the stated hole: low-q companions are
wobble-capable but photometrically dim, so the per-pair delta test
exonerates only the twin-flagged end) are ADOPTED into the 7K design
and THE SPLIT plan.

In flight at close of entry: fullpow/fullpowbe/strictpowbe + 4 seed
extension (pre-reg 0045eea). Credence: ~35% HELD pending the arms/seed
verdicts under the committed standing rules (both-full-arms-fail ->
~50%; seed break -> AMBIGUOUS).

PLAIN VERDICTS: low-end stress: SUCCESS (verdict robust, zero
un-pinned, honesty improved). Power-gate scope: NEEDS REFINEMENT
(conceded; runs in flight under pre-registered bars). Correction #16:
executed as a hard retraction, the kind that hurts and holds.

## Stage 7J-y (2026-07-25, review round 3): C-FAIL — the within-pair
## correlation retracts Part A's scale; CORRECTIONS #17 + #18; the
## verdict label moves to AMBIGUOUS at the now-primary anchoring

Round 3 of the Opus exchange. Two reviewer catches booked as
corrections, one reviewer premise corrected, and the reviewer's new
instrument executed — firing against our Part A.

CORRECTION #17 (label grammar; reviewer-caught inconsistency between
two same-day criteria of ours): the amendment-6 boundary-free
criterion (interior a_marg or gap > -5 -> AMBIGUOUS) was applied to
the seed axis but not the prior-anchor axis. Applied consistently:
"COMPANION-WIN fires at every plausible anchoring" (b5205e2, PAPER
6.3) was WRONG AS LABELED — at the literature anchoring seed 31 is
interior (0.37/0.50) with positive gaps. Correct statement: **the
fenced +99-110 Newton rejection is dead at every anchoring (the
durable result); the LABEL is prior-dependent — COMPANION-WIN at the
measured prior, AMBIGUOUS at the literature prior.** Adopted with the
reviewer's presentation fix: a_marg AS A FUNCTION OF THE PRIOR ANCHOR
is the primary result; verdict labels are annotations on that curve.

7J-y (pre-reg 3f7212e BEFORE execution;
[calcs/stage7j_paircorr.py](calcs/stage7j_paircorr.py),
data/stage7j_paircorr.txt): the reviewer's over-attribution
discriminant — wide-pair components are coeval/co-chemical, so
metallicity/age/extinction displace BOTH deltas (common-mode) while a
companion displaces ONE. Guards designed in: dcolor slicing (ridge-
error confound), attenuation correction (per-component parallax noise
DILUTES here — reviewer's partial-out premise inverted for our
construction, MG uses per-component parallaxes), tail-incidence
channel separated (Tokovinin correlated-subsystem confound). RESULT:
**rho_core_att = +0.465 on the bar slice (N=5423; bars 0.10/0.30);
+0.394 at |dcol| >= 0.4; tail incidence ratio 2.74; flagged star's
partner sits 0.114 mag bright of baseline ==> C-FAIL, decisively.**
~Half the core variance is common-mode astrophysics the 1D mixture
read as companions.

CORRECTION #18 (measurement retraction): **Part A's absolute scale is
RETRACTED — C = 0.41 biased LOW, f_photo = 0.22-0.30 and f_host =
0.42-0.57 biased HIGH, by common-mode over-attribution** (the round-2
scout tension — literature 2-3x below us — is hereby EXPLAINED in the
literature's favor; the round-2 "both branches open" closes onto the
scale-bias branch). The completeness FORMALISM survives; the
single-star model (independent per-star delta) does not. Repair path
= 7J-z: joint 2D (delta1, delta2) mixture — common-mode pair-level
component (diagonal) + per-star noise + companion shoulders
(axis-aligned) — the geometry that separates the channels. QUEUED,
CPU-grade. Until then the literature-anchored prior (peak ~0.22,
sigma ~0.08, scout-grade pending primary sources) is PRIMARY per the
pre-registered C-FAIL branch.

CASCADE, honestly propagated:
  (a) 7J verdict at the primary anchoring = **AMBIGUOUS** (a_marg
      seed-mean ~0.2, seed 31 interior 0.37-0.50, Newton within
      +5-7). Per the original 7J tree, AMBIGUOUS -> seed extension
      decides — and the 4-seed extension is ALREADY RUNNING; its
      cubes are prior-independent, so the literature-prior marginal
      over 6 seeds comes free on completion. THE DURABLE RESULT
      UNCHANGED: the fenced Newton rejection does not survive
      companion marginalization at ANY defensible multiplicity
      (alpha_marg <= ~0.5 everywhere; >= 93% of the +99 gap gone
      even in the most boost-friendly cell).
  (b) The "-768 fence pricing" DEFLATES: under the literature prior
      the old f <= 0.1 fence is ~1.5 sigma low (ln pi ~ -1.1), not
      absurd. 6.3/abstract lines rewritten.
  (c) Correction #16 SOFTENED in one clause (correction-on-
      correction): "Banik substantially closer to the truth than our
      fenced fits" was conditioned on Part A's scale; under
      literature rates (~0.2) our 0.1 fence is CLOSER to the truth
      than their 0.69. What stands of #16: the "photometry forbids"
      LOGIC was invalid regardless (no completeness accounting), the
      companion freedom was directionally legitimate, their 0.69
      stays excluded under every reading, the manufactured-framing
      retraction stands.
  (d) The raw-law "internal inconsistency" claim RESCALES: at
      literature f ~ 0.2 the as-published law's cost is ~ -100..-170
      (not -505..-1078); still a real defect (photocenter
      cancellation is exact physics), milder grade. bin-7j-rawcond
      annotated.
  (e) Correction #15 moves to UNDETERMINED-pending-7J-z: its
      retraction of the 7I cut-survival argument assumed C = 0.41;
      if true C is high (0.6-0.8), the cut-survival argument regains
      power and the median's companion-immunity partially returns.
      No flip-flop: status = pending the repaired completeness.
  (f) The injections in flight validate power at the (now-demoted)
      measured-prior operating point; the literature-anchored power
      test (fullpowlit: truth companions at 0.22 through the
      literature prior) is QUEUED for after batch completion (the
      running batch's script must not be edited mid-flight).
  (g) fpm: the reviewer's round-3 premise was factually wrong (the
      grid WAS extended to 2.4 in amendment 5) but its diagnosis is
      right and sharpened by the facts: the marginal puts 82-99% at
      the NEW edge — the error model lacks a variance channel (the
      6P s-flat scatter object, worth ~+37 lnL there). The proper
      instrument is the per-system scatter component in the marginal
      model, queued with 7J-z; the fpm=1.5-pinned invariance of the
      verdict is the standing guard.
  (h) CREDENCE: anomaly-real ~35% -> **~45%** (stated reasoning: the
      durable negative — fenced rejection dead at every anchoring —
      keeps it well below the pre-7J 60-65; the exact-zero collapse
      was partly manufactured by our own over-scaled prior, which
      returns some mass; pending seeds, arms, 7J-z, 7K).
  (i) Reviewer refinements adopted: the per-pair (q, l) wobble-bound
      instrument for the 11 band pairs (upper bound on companion-
      induced vtilde per pair — decides all eleven on measurement)
      folded into 7K-b; the abstract-(2) dial sentence and the
      one-half-thread/conclusions "measured twice in disconnected
      systems" framings get the fence-conditional flag explicitly
      (round-3 checklist item).

PLAIN VERDICTS: 7J-y instrument: SUCCESS (cheap, decisive, and it cut
against us — the discipline working). Part A: DIFFERENT PHYSICS in
the small — the overluminosity distribution is common-mode-dominated,
which is an astrophysical measurement of this sample in its own right
(metallicity/age coherence of wide pairs) and the reason its
companion reading failed. 7J verdict: NEEDS REFINEMENT — AMBIGUOUS at
the primary anchoring, seeds in flight, 7J-z is the repair, and the
durable finding (no fenced Newton rejection survives marginalization)
stands at every anchoring.

## Stage 7J-y2 + round 4 (2026-07-25): the blending audit CLEARS #18;
## the cadence rule; the four-absorber diagnosis; 7J-g queued

Round 4 (Opus). Executed and adopted:

(1) BLENDING DISCRIMINANT (pre-reg fa29f33; the reviewer's confound-
in-the-same-slot: BP/RP windows are arcsec-wide, so close pairs
cross-contaminate — both go bright — mimicking shared metallicity):
**rho(theta) is FLAT 2"-250"** (core: +0.62/+0.60/+0.52/+0.51/+0.54;
close-end elevation ~15% relative, far under the pre-stated 50%
blending signature); the extinction-distance guard shows the
d-dependence an astrophysical story predicts (d>=100pc rho ~0.6-0.69
vs d<100pc ~0.39-0.52 at ALL theta — blending would track theta, not
d); **THE BLENDING-SAFE BAR (|dcol|>=0.15 core, theta>=8"): rho =
+0.448 -> C-FAIL STANDS independent of blending.** Correction #18
survives its audit; the coherence measurement UPGRADES (genuinely
astrophysical, separation-independent). 7J-z design consequences: no
theta-dependent common-mode component needed; the reviewer's free-
slope note ADOPTED (unit-slope diagonal would push mass-dependent
metallicity/age response off-diagonal onto the companion axes =
recreating the over-attribution subtly; slope + off-diagonal scatter
free, checked against the mass-ratio distribution).

(2) LITERATURE CENTER RECENTERED 0.22 -> 0.166 (Tokovinin 10.0%/7.3%
per component combined per-pair; reviewer arithmetic accepted): all
conclusions center-insensitive across 0.16-0.22 (lit16 rows: full
posterior still fcomp=0.2, a_marg 0.40/0.50 seed 31, 0.00 seed 101;
strict at 0.1). Adopted framing: **the companion prior at 0.2-50 kAU
is UNMEASURED — not "measured high, retracted, defaulting to
literature" — and 7J-z is the FIRST measurement at these
separations: a contribution, not a repair.** Our old fence 0.1 vs
literature ~0.16: low but not dramatic; "nearer the truth than 0.69"
holds without being vindication.

(3) ARM RESULTS (interim, batch mid-flight): fullpow ARM-VALIDATED
(simple recovers own truth through the measured prior: a_marg=1.48,
dN=+94.5; note the ~+0.3 upward marginalization bias — the real-data
simple zero is not arm compression); fullpowbe ARM-FAIL (BE recovers
own truth at 0.73 vs band [0.85,1.45] — the BE-full null is formally
UNVALIDATED; the ~0.65 compression cannot arithmetically explain a
0.00, but the arm carries no formal weight until repaired). Per the
amendment-6 standing rule (>= 1 full arm validated), no PROVISIONAL
downgrade. strictpowbe + 4-seed extension still running.

(4) THE CADENCE RULE (adopted verbatim, effective immediately):
**credence FROZEN at ~45%. No further credence re-booking until BOTH
decisive instruments land (7J-z repaired completeness; 7J-g direction-
channel decomposition). Interim rounds record MEASUREMENTS, not
verdicts.** The reviewer's diagnosis accepted: eighteen corrections
and a 65->35->45 credence in two days means verdicts were being
booked faster than the quantities settled; correction-on-correction
is the signal. The durable statement (the fenced +99-110 Newton
rejection is dead at every anchoring) carries no credence and needs
none.

(5) THE FOUR-ABSORBER DIAGNOSIS (adopted as the structural frame —
the real finding of these three days): excess vtilde width has four
absorbers — boost, near-parabolic orbits, hidden companions, noise
inflation — and one width budget. The verdict has tracked whichever
absorber carries the tightest CURRENT external constraint (tighten
the companion prior -> companions win; retract it -> boost partly
returns; free fpm -> 82-99% to the new edge). Better priors on any
one absorber RELOCATE the degeneracy; only a shape observable breaks
it, and gamma is that observable (companions: direction uncorrelated
with separation vector -> flattens the angular distribution; near-
parabolic: concentrates it; noise: ~isotropic broadening; boost:
scales vtilde). **NEW STAGE 7J-g QUEUED (priority ABOVE further
injections): re-run the four-way contest with gamma COLLAPSED
(vtilde-only likelihood) and quantify the share of (alpha, w_rad,
f_comp, f_pm) discrimination carried by the direction channel. This
is the test that decides whether Paper 1 exists: if gamma carries
the separation, the 2D methodological claim is answered on
measurement; if it does not, no companion-prior work settles the
boost on this catalog.** Implementation after the running batch
(script frozen mid-flight).

PLAIN VERDICTS: blending audit: SUCCESS (#18 confirmed astrophysical,
the safe bar decisive). Recentering: SUCCESS (conclusions center-
insensitive). Cadence rule + four-absorber frame: adopted — the
program's verdict-booking discipline just got its own correction.

## Stage 7J-w + round 5 (2026-07-25): the arm diagnosis — prior-free
## offsets, no arm at the operative anchor, per-pair 7J-z slopes

Round 5 (Opus). All points executed or adopted; measurements only
(cadence rule in force — no credence motion, no new verdict labels).

(1) NO VALIDATED ARM AT THE OPERATIVE PRIOR (conceded, made explicit
everywhere): fullpow validated the simple arm THROUGH THE RETRACTED
MEASURED PRIOR — a cell that no longer exists. Scope facts stated in
PAPER 6.3 + here: cube re-use is free on the SEED axis only (the
prior enters at marginalization); power ARMS at a new anchor need
NEWLY INJECTED companion populations (different truth data = new
cubes). fullpowlit is the only injection bearing on the reportable
number — and per the reviewer's own ordering it runs at whatever
anchor 7J-z lands (GPU at soon-to-be-vacated cells is waste):
**ORDER FIXED: 7J-z -> 7J-g -> the arm suite at the landed anchor
(incl. the BE alpha=0.4 additive-vs-multiplicative discriminator).**

(2) THE +0.3 DIAGNOSED ([calcs/stage7j_armdiag.py], zero GPU, from
cached pow cubes): the answer is the dichotomy's third way — the
offsets are PRIOR-FREE. PROF (no prior) = 1.47 vs matched-prior MARG
1.48 (simple, truth 1.18); PROF = MARG = 0.73 (BE, truth 1.13). And
re-marginalizing the SAME injected data under the literature prior
(0.16/0.22) moves NOTHING (a-hat identical; the injected likelihood
is informative enough about fcomp to override a mismatched prior —
posterior stays at 0.50 against the penalty): prior-misspecification
does not mechanically manufacture alpha when the likelihood is
informative. Remaining candidates for the +0.30/-0.40 pair: the
Stage-3A realization layer (one truth pop, one fitting pop) or a
genuine estimator bias — SINGLE-INJECTION GRADE, multi-truth
recovery map queued with the post-7J-z arm suite. Cross-family
consistency note: BE-truth 1.13 read by the simple arm at 0.94 =
the standard BE/simple ratio — the arms disagree only about their
OWN truths, opposite signs.

(3) BE ARITHMETIC CLAUSE RETRACTED (reviewer-caught, unnumbered
interim-note fix): my "a 0.65 compression cannot produce 0.00 from
1.13" assumed a multiplicative recovery map; an additive -0.40
offset fits the same single point and sends any truth below 0.4 to
zero. One injection cannot distinguish them; the disputed region IS
alpha 0.2-0.5. The BE alpha=0.4 injection resolves it (queued at the
landed anchor). Meanwhile the phrasing rule adopted: "both laws" is
DROPPED from the verdict sentences — the simple arm is quoted, the
BE arm is flagged pending (PAPER abstract + 6.3 edited).

(4) 7J-z DESIGN SHARPENED (adopted verbatim): the common-mode
diagonal slope must be PER-PAIR, set from the two components' colors
via a mass-dependent metallicity response (0.1 dex moves a mid-M
dwarf several times further than an F star; M_G span 2.6-14.2 means
F+M diagonals far from unit slope, twins near it); only amplitude
and scatter free. A single global slope would land the unequal-mass
pairs' common-mode on the companion axes = the retracted over-
attribution returning through the door being closed. Implementation
note: the response ratio can be self-calibrated from the data
(slope field vs (col1, col2) on the core, iterated against the
mixture's posterior singles) with an isochrone cross-check.

(5) THE INTERIOR READING NAMED AS A RESULT (adopted): at the
literature anchor the fit sits at a SPECIFIC value — alpha ~ 0.4-0.5,
roughly HALF the galactic calibration, interior on one realization —
the split-width-budget signature, precisely the object 7J-g decodes.
Stated so in PAPER 6.3 (not as a label failure).

(6) Blending-audit caution clause logged (reviewer): the rho
d-gradient has a second reading — distance shifts the sampled mass
range through the magnitude limit, and metallicity sensitivity is
mass-dependent — which does not touch the blending conclusion but
belongs next to the coherence measurement's interpretation.

PLAIN VERDICTS: arm diagnosis: SUCCESS as a measurement (prior-free,
mismatch-immune at informative likelihood; single-injection grade
stated). Scope corrections: executed. The stage's standing summary
is unchanged by round 5 — AMBIGUOUS at the primary anchoring, the
durable Newton-rejection collapse at every anchoring, everything
else queued behind 7J-z -> 7J-g.

## Stage 7J-e3/7J-w2 + round 6 (2026-07-25): the review exchange
## closes — the knee measured OUT, the revival mechanism relocated,
## the injection optimism quantified

Round 6 (Opus's last: "the curve-sampling and the curvature check are
the last two things I'd do before 7J-z"). Both executed from cached
cubes; both returned answers sharper than the questions.

(1) THE KNEE (7J-e3, [calcs/stage7j_lowend.py]): Opus's catch — my
"no pending instrument un-kills it" was carried over a curve sampled
only at anchors >= 0.166, while the fenced fit (hard wall f <= 0.1)
sits at alpha = 1.06 / +99.5, so the whole detection appeared to live
between 0.10 and 0.166, with 7J-z able to revive it at 0.12. The
CLAUSE IS RETRACTED — and the fine curve replaces the interpolation
with a measurement that changes the geometry: **the curve is nearly
FLAT — alpha_marg = 0.18-0.31, Newton within +2..+4, for EVERY smooth
anchor 0.06-0.30 (sigma 0.05 AND 0.03, both laws), including anchors
BELOW the old fence — because the likelihood alone pays +12..+28 to
sit at the fcomp = 0.2 cell. The fenced detection returns ONLY behind
a hard exclusion of f >= 0.2, which no rate measurement can justify.
THERE IS NO KNEE IN SMOOTH-ANCHOR SPACE.** Consequences: (a) 7J-z at
0.12 revives nothing — its output reads position on a nearly flat
curve, so the precision-requirement question dissolves (any honest
sigma_f suffices for the lookup); (b) THE REVIVAL MECHANISM
RELOCATES: what could materially move alpha is the forward model's
MISSING WIDTH CHANNEL — the 0.2 cell buys its +12..+28 by supplying
variance the model otherwise lacks (the same object as the winning
cell's -60..-116 misfit, the fpm edge, and 6P's s-flat scatter) — so
the model-side variance repair, not the rate, is the pending
instrument that can un-kill the rejection; (c) grid caveat stated:
no cube cells between 0.1 and 0.2, so the curve is the exact marginal
output on this grid but cannot resolve conditional structure inside
the gap (denser-fcomp cells queued with the post-batch runner edit).
Every "dead at every defensible anchoring" claim in PAPER/banner
replaced by the measured form.

(2) INJECTION OVER-INFORMATIVENESS (7J-w2,
[calcs/stage7j_armdiag.py]): confirmed and quantified — real
full-sample fcomp profiles vary by TENS across 0.1-0.35 (0.1 vs 0.2:
~12; 0.2 vs 0.35: ~100) where the injected profiles CLIFF by
hundreds-to-thousands around their truth cell (0.5 vs 0.35: -60/-81;
vs 0.2: -368/-427; vs 0.1: -1014/-1129) — a 4-30x informativeness
excess. The model's own companions carry a cleaner fcomp signature
than the sky's width budget permits (the same missing-channel object
again, seen from a third side). Consequence, stated for the
post-7J-z arm suite: **arm validations at any anchor overstate power
on data; the suite must degrade its injections with the measured
common-mode or state the excess alongside every recovery.**

(3) THE EXCHANGE CLOSES (Opus's own framing accepted): five rounds
of checklist->execution->checklist is the 2.4 dynamic and it stops
here. The two runs that determine everything are 7J-z (per-pair-
slope 2D completeness = the first subsystem-rate measurement at
these separations, WITH the width-channel repair folded into the
same modeling round) and 7J-g (does gamma separate the four
absorbers = whether Paper 1's methodological claim operates in the
regime that matters). Everything else — arm suite at the landed
anchor with degraded injections, denser-fcomp cubes, Tokovinin
primary-source requote — queues behind them. Next conversation =
their results.

PLAIN VERDICTS: 7J-e3: SUCCESS — the knee question answered by
measurement (flat curve; revival mechanism relocated model-side; my
"no pending instrument" clause retracted and replaced by a sharper
true statement). 7J-w2: SUCCESS — the optimism quantified before it
could inflate an arm verdict. Review exchange: closed with the
program better than it entered — one Part A scale retracted, one
coherence measurement gained, one discipline rule adopted, and the
two decisive instruments named and ordered.

## Stage 7J-s6 (2026-07-25): the six-seed budget closes the batch —
## no seed break; the operative-anchor measurement gets error bars

Batch complete (pre-reg 0045eea; [calcs/stage7j_seed6.py],
data/stage7j_seed6.txt). Measurements only (cadence rule).

STRICTPOWBE (final arm leg): ARM-FAIL (BE-strict null unvalidated) —
but in the OTHER failure mode: position IN band (a_marg = 0.90 vs
truth 1.13, band [0.85, 1.45]) with margin below bar (dN = +14.7 vs
>= +20). Note logged: the +20 bar was not N-scaled to the 57% strict
sample, so this failure is partly bar calibration; the bar is the
bar. Arm scorecard final: simple-full VALIDATED (at the retracted
anchor), BE-full FAIL (position), BE-strict FAIL (margin); no arm at
the operative anchor (unchanged); suite re-runs post-7J-z with
degraded injections.

THE SIX-SEED READ (both anchorings, 6 seeds x 2 samples x 2 laws):
- Measured prior (continuity row; its prior retracted): a_marg = 0.00
  in 24/24 reads, gaps lm(0.5)-lm(0) = -17..-61. **AMENDMENT-6 SEED
  RULE: NO BREAK in 16 new-seed reads** — the measured-prior verdict
  stands on the seed axis.
- Literature anchor 0.16 (operative): **full a_marg = 0.13 +- 0.07 SE
  (simple) / 0.30 +- 0.09 (BE), Newton within +1.0/+4.1; interior
  scraps persist on 2/6 (simple) and 4/6 (BE) seeds (0.31-0.55, gaps
  to +12.9); strict = 0.00 in 12/12.** The flat-curve statement now
  carries six-seed error bars: the operative-anchor residual is a
  scrap at 0.1-0.3 — an order below the fenced 1.06-1.18 and below
  any detection grade — consistent with the split-width-budget
  reading 7J-g decodes. The BE-over-simple ordering of the scrap
  (0.30 vs 0.13) is noted, not interpreted (the BE arm is the
  unvalidated one).

No label motion, no credence motion (frozen ~45%): these numbers ARE
the interim record. NEXT = 7J-z -> 7J-g as fixed.

PLAIN VERDICT: batch SUCCESS — the seed axis is closed (no luck
break), the operative-anchor measurement is now error-barred, and
the program's full attention moves to the two deciders.

## Stage 7J-z / 7J-g PRE-REGISTRATION (2026-07-26, committed before
## execution; this commit's hash is the timestamp)

THE TWO DECIDERS, instrumented in one round:

7J-z part 1 — [calcs/stage7jz_mixture.py](calcs/stage7jz_mixture.py):
the repaired completeness measurement. Joint 2D (δ₁, δ₂) mixture with
PER-PAIR common-mode response slopes r(c₂)/r(c₁) (self-calibrated from
the cross-component covariance field, which is companion-immune under
centering + independence), companion sector = part A's combined-δ
amplitude law unchanged, noise = per-component σ_n(c) + parallax term.
Gates: GZ0 regression (ρ_core_att = +0.465, N = 13784); GZ1a/b the
over-attribution control (companion-free rank-1 AND rank-2 truths →
f̂ ≤ 0.03/0.05; GZ1b fail → rank-2 promoted, pre-registered); GZ2
recovery (0.10/0.25 within max(0.03, 25%)); GZ3 nesting/convergence;
GZ6 postdiction (ρ bar slice ±0.05, the dcol-slice TREND ±0.07 — the
rank-1 geometry's own test — flag fraction ±0.02; incidence ratio
reported, no bar). Output = the LANDED ANCHOR: f-profile envelope
(flat-q × q^-0.5) through part A's host remap verbatim →
data/stage7jz_prior.npz. Axes per component throughout (closes a
convention ambiguity in the lit-0.16 anchor: Tokovinin's 10.0%/7.3%
are per-decision-level population fractions; per-pair 0.166 ≈
per-component ~0.09 — the flat curve made 7J-s6's conclusion immune,
noted here, no numbered correction). Fail branch: anchor NOT shipped,
part 2 reads lit-only.

7J-z part 2 — amendment 7 in
[calcs/stage7j_marginal.py](calcs/stage7j_marginal.py): 'photow' mode
= photo + the per-system WIDTH CHANNEL (vtn → vtn·exp(sq·g_i), g_i
drawn once per system, applied after the selection cut =
normalization-error semantics, the 3E σ_m / 6P s-flat object;
γ untouched by construction). SQ_GRID = [0, 0.1, 0.2, 0.3]. GB0w:
sq = 0 slice must equal the cached photo cube ≤ 1e-3 (GB0p precedent
bar). Batch: FULL sample, seeds 31/101, both laws. Reader bars
([calcs/stage7jz_read.py](calcs/stage7jz_read.py)): D1 channel usage
(P(sq>0) ≥ 0.7; edge flag at P(sq=0.3) ≥ 0.5), D2 multiplicity-cost
closure (the −60…−116 object; ratio ≤ 0.5 CLOSED / ≥ 0.8 NOT-CLOSED),
D3 fpm-edge release; VERDICT at the landed anchor: BOOST-REVIVES
(a_marg ≥ 0.5 AND dN ≥ +25, either law) / NO-DETECTION (a_marg ≤ 0.3
AND dN ≤ +10, both) / else AMBIGUOUS-CARRIED; extension rule
|Δa_marg| > 0.25 between seeds → +202/303 first.

7J-g — [calcs/stage7jg_read.py](calcs/stage7jg_read.py): the same
batch emits γ-COLLAPSED TWIN cubes (same model histograms, γ summed
out of data and model — the collapse is of the same model, disclosed,
so the metrics isolate the γ CHANNEL). Metrics at the landed anchor,
seed means: M1 α-width shrink ≥ 30%; M2 dN(2D) − dN(vt) ≥ +15; M3
|a_marg(2D) − a_marg(vt)| ≥ 0.30; M4 absorber grade (SD of wr or
fcomp shrinks ≥ 30%). THREE-TIER VERDICT: SEPARATION-CONFIRMED (α
grade, M1-M3) / ABSORBER-LEVEL (M4 only) / SEPARATION-ABSENT (all
quiet) / else AMBIGUOUS. PAPER-1 RULE pre-committed: ABSENT → the 2D
methodological claim does not operate at the landed configuration →
Paper 1 reframed around the model-light channels; CONFIRMED → the
claim operates. CREDENCE MAP pre-committed (releases the ~45% freeze
when BOTH deciders have reported): BOOST-REVIVES + separation ~55-60%
(ABSENT: ~50-55%); NO-DETECTION + separation ~35-40% (ABSENT: ~30-35%
+ reframe); any AMBIGUOUS → hold ~45%.

Scout note (Haiku, primary-source pass, this morning): Tokovinin's
10.0%/7.3% CONFIRMED verbatim at AJ 147, 87 (products 0.464×0.214 and
0.464×0.157; full 67-pc F/G sample, all separations); nearest direct
wide-binary subsystem survey = Tokovinin 2011 (AJ 140, 510), 12±4% at
outer 150-1000 AU only; "UNMEASURED at 0.2-50 kAU" SURVIVES
(scout-grade), with the caveat that subsystem incidence is stated to
GROW with outer separation — the lit anchor may sit low for our range.
Part 1 is then the first measurement in its range if gates pass
(claim held at scout grade until the paper's own novelty pass).

AMENDMENTS (2026-07-26, in-flight, all logged BEFORE any decider
result was read; review round 7 = Opus's post-batch note):

7e — GB0w FIRST FIRING (gate catch, pre-results): the photow batch
aborted seed-31-simple at max|photow(sq=0)−photo| = 8.19e+02. The
diagnostic localizes it exactly: diff = 0.000 at fcomp = 0, growing
monotonically along the fcomp axis (−29/−117/−244/−377/−574),
kw-dependent — the companion-sector signature. Cause: build_pop's
photocenter amplitude branch still tested AMP == 'photo', so photow
fell through to the AS-PUBLISHED raw wobble law. One-line fix
(AMP in ('photo','photow')); the invalid cube pair DELETED before
relaunch (the exists-check would have resurrected it — the 7D
stale-table lesson applied); the diagnostic's peek at the invalid
cube's sq-axis means is disclosed and superseded. The gate class
earns its keep: this was invisible in the PROF row and would have
poisoned both deciders.

8 — 7J-g DUAL-CONFIGURATION READ (review round 7, adopted): the
width channel sq is γ-blind by PHYSICS (normalization errors do not
rotate velocities) and by construction, so its inclusion dilutes the
γ-carried share of total information and pulls the M-metrics toward
"no separation" partly by design. stage7jg_read now reads BOTH
configurations from the SAME cubes: sq FREE (pre-registered primary)
and sq PINNED to 0 (the original four-absorber diagnosis). Primary
tier from sq-free as pre-registered; the sq-pinned co-read is the
interpretation guard; if the tiers differ, the difference is itself
the finding. (Cheaper than the reviewer's request — no second batch:
the sq=0 slice of the same cubes IS the pinned configuration.)

THE AXIS CHAIN (review round 7's demand, written before GZ1 landed;
every companion number now carries its axis):
  (P) = blended-photometric fraction PER COMPONENT — what any
      δ-distribution fit measures (only companions inside Gaia's
      resolution blend and brighten);
  (H) = host fraction PER COMPONENT = (P) ÷ P_blend, with P_blend ≈
      0.52-0.65 from the model's own (q, logP) laws at the sample's
      distances (part A's four combos). THIS is v7's fcomp axis (the
      code draws each component's companion independently at fcomp).
  (PAIR) = 1 − (1−H)² — comparison-only; never enters v7.
  Literature (Tokovinin, scout-confirmed): subsystem rates are REAL
      companion rates → (H) ≈ 0.09-0.10; (PAIR) ≈ 0.166. The 7J-s6
      "lit-0.16" anchor therefore read a (PAIR)-flavored number onto
      the (H) axis — ~1.7× high. Conclusion immune (the 7J-e3 curve
      is flat 0.06-0.30 at smooth anchors), but the row gains this
      annotation; the operative anchor is about to be the measured
      one regardless.
  Retracted part A: f_photo 0.22-0.30 = (P); f_host 0.42-0.57 (peak
      0.51) = (H) PER COMPONENT via ÷P_blend — NOT via the pair map;
      the numerical coincidence 1−(1−0.24)² = 0.42 is a trap, named
      here so nobody falls into it again.
  Part 1 interim MAP (envelope + gates pending): f̂ = 0.166 = (P) →
      (H) ≈ 0.26-0.32 → v7-axis position between the 0.2 cell and the
      0.35 posterior cell; (PAIR) ≈ 0.45-0.54; vs literature on (H):
      ×2.6-3.2. CORRECTED interim comparison (supersedes the mixed-
      axis sentence "between the literature and the retracted scale"):
      on (P) the repair pulls the blended fraction DOWN from part A's
      0.22-0.30, as the #18 mechanism predicts; on (H) the landed
      position is in the companion-favored region of the grid and
      ~3× the literature host rate (the scout's incidence-grows-with-
      separation caveat and the evolved-interloper upper-bound caveat
      both push that same direction).
  The knee, stated precisely: 7J-e3 measured NO knee for smooth
      anchors (σ 0.03/0.05) anywhere in 0.06-0.30 — the +99.5 row
      lives behind a HARD WALL only. But a measured prior tighter
      than σ ≈ 0.02 landing LOW would begin to reconstruct the wall
      (penalty at 0.2 ≥ +16 vs the likelihood's +12..+28 payment), so
      the shipped envelope's WIDTH is quoted next to its position and
      the readers interpolate the actual envelope, never a label.

## Stage 7J-z part 1, v1 VERDICT (2026-07-26): honest GZ6 fail — the
## branch fires, and the failure is the discovery

v1 ([calcs/stage7jz_mixture.py](calcs/stage7jz_mixture.py)) ran to its
gates: GZ0 exact PASS; convergence perfect (two starts Δ0.0);
**rank-1 MAP f̂ = 0.1657 (blended, per component)**, r(c) 0.11-0.29
rising redward, σ_n collapses to 0.04-0.06 mid-MS (the pair
correlation carries the width); companions demanded at +4361 over
f=0; gray rank-2 preferred +334.5 (σ_g=0.058) with **f̂ STABLE at
0.168**; GZ3 PASS. ALL FOUR INJECTION GATES PASS AT ~EXACT (GZ1a
0.000, GZ1b 0.000 — the #18 over-attribution mechanism is dead in
this estimator even under rank-2 mis-specification; GZ2 0.103/0.251).
**GZ6 postdiction FAILED (ii)+(iii)** → per the pre-registered
branch: NO PRIOR SHIPS from v1; deciders read LIT16 unless a repaired
instrument passes its own gates.

The diagnosis (read-only, logged): (iii) was OUR BAR'S construction
mismatch — model simulates on the 12,084 windowed pairs, bar compared
to the full-sample 0.121; like-for-like the windowed data give 0.107
vs model 0.096 = INSIDE ±0.02. (ii) is REAL: the measured ρ(dcol)
profile decays smoothly +0.86/+0.77/+0.71/+0.71/+0.60/+0.48/+0.36
across the full range (ridge-error confound 80× too small to matter)
— **the pair-common displacement is a VECTOR (abundance pattern +
age) whose projections decorrelate as the component masses separate;
twins project identically.** A scalar-z model structurally cannot
make ρ(dcol→0) ≈ 0.86 (v1 postdicts 0.47). Unmodeled twin-end
correlation is an over-attribution channel the GZ1 injections could
not see (their truths lacked the structure) — the exact #18
direction, at the q≈1 corner. Excess correlated variance at
dcol<0.1: ~0.084 mag (20.7% of core pairs).

## Stage 7J-z2 PRE-REGISTRATION (v2, the coherence-kernel instrument;
## committed before execution)

[calcs/stage7jz2_twin.py](calcs/stage7jz2_twin.py): v1 + k(dcol) =
exp(−dcol/λ) coherence kernel — common displacement split into
pair-shared (√k·z) and component-private (√(1−k), merges EXACTLY into
per-star variance; quadrature stays 2D; 13 params). Gates: GZ0; GZQ
quadrature (13×5 vs 15×7 ≤ 1.0); GZ1c-a over-attribution on
kernel-truth (≤0.03); GZ1c-b DIAGNOSTIC = v1-style (λ pinned) fit on
kernel-truth → measures v1's kernel-blind over-attribution; GZ2
0.10/0.25; GZ3 nesting incl. L(v2) ≥ L(v1)=2497.9 and λ interior;
GZ6′ like-for-like (windowed flag 0.107 ±0.02; both ρ slices ±0.07 —
the failed one now reachable). DECISION RULE: all pass → v2 envelope
= THE LANDED ANCHOR; any fail → LIT16 stands, no third instrument
today. THE NUMBER: f̂(v2) vs v1's 0.166 = the data-level
over-attribution measurement either way. v1's npz quarantined as
_v1_unshipped on its completion (existence ≠ shipped — reader treats
existence as landed, so the quarantine is load-bearing).

10 — v2 KILLED BY ITS OWN GZQ MID-RUN; v2b = EXACT INTEGRATION
(logged before any v2 gate verdict printed; ordering disclosed): v2's
GZQ fired at d = 18.7 (bar 1.0) — with the kernel fitted, twin pairs'
effective variance collapses and the z-integrand becomes a spike
narrower than the GH node spacing; the quadrature carries O(10-20)
lnL error, and its five-start dispersion also tripped the convergence
bar (gap 43.1, short-λ basin). Run ABORTED (the GB0w precedent). The
fix is EXACT, not finer: every companion state's (z, gray) integral
is a closed-form bivariate-normal marginal — no nodes, ~10× faster
([calcs/stage7jz2b_exact.py](calcs/stage7jz2b_exact.py); GZQ replaced
by GZS11 state-11-thinning ≤ 0.5; five dispersed starts, top-two ≤
2.0; SAME model, SAME bars, SAME decision rule — not a third
instrument: the model never changed, its integration was corrected by
its own gate). RETRO-FLAG: v1's GH quadrature is the same class
(injection-validated and self-consistent, so its f̂ calibration
stands; its lnL SCALE carries undisclosed-magnitude quadrature error
— noted on the v1 ledger row at close-out). v2 interim record (approx,
nothing shipped): lnL 2847.9, f = 0.162, λ = 2.95 — the kernel worth
~+350 and f stable across all three model classes (0.166/0.168/0.162)
at approximate grade; v2b's exact numbers decide.

9 — THE ANCHOR CURVE RE-READ (review round 8, adopted pre-results):
the 7J-e3 flatness was measured under the OLD absorber configuration;
the width channel competes for the same budget, so stage7jz_read now
re-measures the smooth-anchor curve on the photow cubes (same centers,
extended to 0.34; σ 0.05/0.03; reading rules pre-registered:
KNEE-REAPPEARS if any 0.06-0.30 anchor reaches a_marg ≥ 0.5 with
dN ≥ +25; FLAT-PRESERVED if σ=0.03 span ≤ 0.25 with all dN ≤ +10;
else INTERMEDIATE). Round 8 also: the reviewer's pair-map arithmetic
withdrawn by the reviewer (the 1−(1−0.24)² = 0.42 match confirmed as
the named trap — his own warning's failure mode, executed one
paragraph after warning about it, conceded openly); the fence/anchor
conflation withdrawn; his standing caution adopted verbatim into
discipline: "the gate that catches a bug is always one someone
thought to write" — every new mode gets a regression against the old
mode at its identity point (goes to CLAUDE.md discipline at
close-out).

## Stage 7J-z2b + THE TWO DECIDERS (2026-07-26): the day closes

7J-z2b ([calcs/stage7jz2b_exact.py](calcs/stage7jz2b_exact.py), exact
BVN likelihood): the resolution cascade ended honestly — GZS11 failed
at d = 2.08 (with twins modeled correctly their kernels are ~0.04 mag
wide, comparable to the q-grid spacing; every discreteness now shows
at multi-lnL scale) and GZ3 failed on the pre-polish top-two start
gap 2.8 vs 2.0 (a λ-ridge multimodality; every OTHER GZ3 condition
passed: +475.9 over v1, companions +4918, λ interior). Everything
else PASSED: GZ0; GZ1c-a = 0.000; **GZ1c-b THE ACQUITTAL: a
kernel-BLIND (v1-style) fit on a full-kernel companion-free sky also
returns f̂ = 0.000 — the twin-coherence over-attribution channel is
measured ABSENT; v1's 0.166 was never inflated by it**; GZ2 exact
(0.099/0.254); **GZ6′ FULL PASS incl. the twin slice (+0.743 vs
+0.796) and like-for-like flag (0.093 vs 0.107) — the kernel model
postdicts every measured structure that killed v1**. Per the
pre-committed rule: NOTHING SHIPS today (two numerical-protocol
bars), no third same-day iteration; **v2c queued** (profiled-λ +
resolution-gated continuous q-integral — a certificate job: the
physics is already five-basin, two-integrator, three-model-class
stable). UNSHIPPED RECORD: **f̂ chain 0.166 → 0.162 → 0.159; blended
envelope [0.17, 0.17]; host axis [0.29, 0.32] peak 0.32** — pricing
the literature host center 0.09 at −1634: the wide-binary subsystem
rate at 0.2–50 kAU is ~3× the field rate (the direction Tokovinin's
separation trend predicted; a publishable measurement once
certificated).

DECIDER 1 — 7J-z part 2 ([calcs/stage7jz_read.py](calcs/stage7jz_read.py),
LIT-CONDITIONAL per the fallback rule): **AMBIGUOUS-CARRIED** — the
seed-mean missed the pre-registered detection bar by 1.2–1.8 lnL
(α_marg = 0.74/0.70, dN = +23.8/+23.2 vs the ≥ +25 bar). The
boundary-free statement (the #17 standard, and the day's headline):
**given the width channel, α_marg ≈ 0.70–0.74 at ΔlnL ≈ +23 with the
ANCHOR CURVE FLAT TO 0.01–0.03 ACROSS EVERY ANCHOR 0.06–0.34** —
the α answer is companion-prior-INDEPENDENT; at v2b's would-be
measured anchor (0.29–0.32) it reads 0.67–0.74 at +12…+16. The
width channel is demanded (D1: P(sq>0) = 1.00, sq = 0.2 interior,
zero edge mass — the 3E/6P object finally has its parameter and its
value) and it does NOT close the multiplicity cost — it INVERTS it
(D2 ratio 1.4: with sq available the kinematics want fcomp = 0.1 at
every anchor; the tension between measured multiplicity ~0.3 and
kinematic preference 0.1 is now a named object, absorbed as dN
attrition at high anchors, not α motion). D3: BE still rides
fpm = 2.4 at the marginal (P = 0.97; simple 0.43) — correction-#4
flag, FPM extension to 3.0 queued with the arm suite. Continuity:
the retracted-prior read still gives α = 0 (the old COMPANION-WIN
lives in the sq = 0 slice); the ledger keeps both. Anchor-curve
label INTERMEDIATE by the letter of the pre-registered rule — the
rule conflated "flat" with "null"; the measured shape is FLAT-α at
NON-null dN, the flattest curve in the program (disclosed as a
bar-design miss, no renumbering).

DECIDER 2 — 7J-g ([calcs/stage7jg_read.py](calcs/stage7jg_read.py),
dual-configuration per amendment 8): **PRIMARY (sq free) =
ABSORBER-LEVEL SEPARATION** (M4: SD(w_rad) 0.015–0.035 → 0.000 both
laws; α channel-robust: 2D 0.74/0.70 vs vt 0.62/0.66, dN diffs
+3.6/+4.6) — under the width model the boost is in the SPEEDS and
the direction channel's job is pinning the absorbers. **CO-READ
(sq = 0, the four-absorber configuration) = SEPARATION-CONFIRMED at
α grade (M3)**: ṽ-only fitting reports a PHANTOM α ≈ 0.50/0.52 at
+12.8/+12.2 which the direction data veto to 0.20/0.25 at +2.4/+3.7
— **in the width-less configuration every published ṽ-only pipeline
occupies, the γ channel is the difference between a fake detection
and honesty** (stated as a reading; published models differ in their
error sectors). THE TIER DIFFERENCE IS THE FINDING: γ and the
per-system width channel are partially interchangeable
absorber-police; γ's α-grade power exists exactly when the width
channel is absent. GVT: α = 0 rows law-blind at 0.00e+00 both
channels.

PAPER-1 RULE APPLIED (pre-committed): ABSORBER-LEVEL primary + α-grade
co-read ⇒ **Paper 1 PROCEEDS** — the 2D methodological claim operates:
at nuisance grade in the five-absorber model (and the α bottom line is
quoted as channel-independent, a STRONGER reality statement than
γ-dependence would have been), at α grade in the four-absorber
configuration (the phantom-veto demonstration = the methods chapter's
centerpiece). CREDENCE MAP APPLIED (pre-committed): the 7J-z verdict
is AMBIGUOUS ⇒ **anomaly-real HELD AT ~45% by the map's own rule**
(the cadence-rule freeze is RELEASED — both deciders reported — and
the value is re-affirmed, not inherited; the miss against the
detection bar was 1.2–1.8 lnL, direction noted, no motion booked).

PLAIN VERDICTS: 7J-z part 1 = NEEDS REFINEMENT (v2c certificate
queued; the measurement itself is stable and the acquittal is in);
7J-z part 2 = SUCCESS as an instrument (the width channel measured,
the anchor dependence dissolved), AMBIGUOUS-CARRIED as a verdict;
7J-g = SUCCESS (both tiers informative; Paper 1 exists). The day's
scoreboard: two phantom results prevented (the GB0w wobble-law bug;
the [sq0 vt] phantom boost), one channel measured (sq = 0.2), one
acquittal (v1's rate), one unshipped first-of-kind measurement
(subsystem rate at 0.2–50 kAU ≈ 3× field), five corrections of our
own instruments logged pre-results, zero bars moved.

## Round 9 (2026-07-26): the adequacy decomposition — the reviewer's completeness question lands (7J-z3)

Opus round 9, same day as the deciders' close. Adjudication point by
point, then the executed item.

**Adopted-and-executed: the absorber-set completeness question.** His
premise was wrong in the letter — the statistic he asked for ("does the
−60…−116 close?") was pre-registered and shipped as D2 (bars ratio ≤0.5
CLOSED / ≥0.8 NOT-CLOSED; fired NOT-CLOSED at 1.37–1.42, printed in the
verdict line) — but the question behind it had not been decomposed, and
the decomposition is the decision-relevant object. Stage 7J-z3
([calcs/stage7j_sqclose.py](calcs/stage7j_sqclose.py), read-only cube
arithmetic, frame + reading bands committed 5024ecb BEFORE execution;
G0/G1 regressions vs the shipped read exact, 4/4 PASS):

- S1 anchor strain = 0.0 in all four seed-law reads — α_marg = 0.74/0.70
  IS the width-complete model's own unconstrained optimum; the LIT16
  anchor does no work at all. The kinematic face of adequacy is CLOSED.
- S2 free-optimum cell: α interior (refined 0.65–0.75), wr = 0.20,
  fcomp = 0.10, sq = 0.2, fpm 2.1 (simple) / 2.4-edge (BE — extension
  to 3.0 stays queued).
- S3 THE FIFTH-MOVE EXPOSURE: forcing fcomp ≥ 0.35 (the measured-host
  grid cell) collapses α̂ to 0.00 at dN = +0.0 in ALL FOUR reads —
  **FIFTH-MOVE-LIVE by the pre-stated band.** α is fully exposed to the
  multiplicity tension. If the ~0.3 host rate certificates AND the
  per-companion kinematic signature is as modeled, α goes to zero — at
  a cell that still misfits the kinematics by 135–153, so nothing fits
  there and the model-light channels remain. Three resolutions, two
  already first in the queue: (a) the rate fails its v2c certificate;
  (b) rate + signature both right → α → 0 (the fifth move, openly
  priced now); (c) rate right but the per-companion ṽ signature
  over-predicted — the kernel discovery's own direction (common-mode
  vector displacement masquerading as companion light inflates the
  modeled wobble/ℓ per companion) — decided by the arm suite with
  kernel-degraded injections. Queue order unchanged (v2c → anchored
  re-read → arm suite), now carrying the exposure explicitly. NO
  credence move (cadence rule: v2c + the arm suite are the deciders;
  AMBIGUOUS-CARRIED already withheld detection status with D2 open).

**Adopted: sq's identity is open, and that is a real gap vs our own
§5.3 standard.** Concessions made explicit in PAPER §6.3: (i) his
arithmetic is right — sq = 0.2 (20% per-system velocity scale) exceeds
the α = 0.74 signal (~12% in ṽ at wide separations); (ii) identity
candidates with current bounds: photometric-mass channel measured 12×
too small (3J: 0.024), the v2b coherent common mode ≲3% in ṽ by the
same MS translation even taken fully coherent, distance errors 1–2%,
companion-shaped broadening rejected on shape (3K: −420), measurement-
error inflation separately carried (fpm) and STILL edge-riding — no
named channel produces 20%; identity = open item. What keeps it from
being "a free parameter that resurrected the boost": the channel's
EXISTENCE was measured four independent ways BEFORE it was
parameterized (3E σ_m ≈ 0.2–0.25; 6P s-flat +37; the −60…−116 misfit;
the fpm edge), the parameterization was queued by the reviewer himself
in round 3, pre-registered, and the data then demanded it interior with
zero edge mass. And the median boost 1.078/1.086 is sq-IMMUNE (a
symmetric per-system smear moves widths, not medians) — the same
immunity he grants Cookson's flatness null protects our model-light
anchor.

**Adopted-and-made-explicit: the structural-robustness argument.** With
a γ-invisible width absorber free, further γ-invisible systematics load
onto sq, not α — α is identified by direction-resolved shape. This does
make w_rad load-bearing, and w_rad has an EXTERNAL anchor: freezing it
to Hwang's 0.21 moved α̂ by ≤0.01 (7I ablation W, v7 pipeline; stated
with that provenance). Booked in §6.3.

**Adopted: grid-quantization phrasing.** "SD(w_rad) → 0.000 exactly" is
wrong — WR_GRID step is 0.10, SQ step 0.1; correct statement: ALL
POSTERIOR MASS ON ONE GRID NODE (SD below one step). Fixed in PAPER
§6.3, the worldtable banner, and the bin-7jg-gamma ledger note. Same
round: the constant SD(sq) = 0.927 in the sq0 rows of
stage7jg_read.txt explained = sqrt(0.86), a broadcasting artifact of
the pinned length-1 sq axis against the full length-4 grid (μ
spuriously 0.6) — cosmetic; no verdict metric (M1–M4) reads it; the
sqfree SDs are length-matched and sound.

**Adopted (framing): Paper 1 leads with the instrument, α as a
conditional result.** His reframe = our 7I SPLIT plan sharpened: spine
= the width channel + the phantom-veto + the absorber accounting (the
mechanistic reading of the field's decade of disagreement: one
unmodeled width budget of order the signal, and what each group
concludes depends on which absorber they let eat it — Banik's 69%
companions, the fenced detections, Chae — stated as a READING);
deliverables = completeness, the kernel discovery, the subsystem rate
(post-certificate), the wobble-law inconsistency, the census; α quoted
on its flat curve with the sq-identity flag and the S3 exposure. His
Cookson caveat adopted verbatim (a symmetric smear inflates width
without moving a median → their flatness null is NOT explained by the
width channel; stays the sharpest external tension; 7L unchanged in
queue). His "3× field is expected from triple-dynamics formation" =
plausible, NOT printed pending a scout (queued with v2c).

**Refuted/refined:** (1) "you have a direct test and you haven't
reported it" — D2 was pre-registered, computed, and printed with its
NOT-CLOSED token in the shipped read; what was missing was the
decomposition, now done. (2) "each move was caused by changing the
absorber set rather than by new data" — the 0.00 → 0.13–0.30 move was
caused by DATA (7J-y's ρ = +0.47 paircorr measurement → correction #18
→ literature anchoring); the trajectory is absorber accounting
converging under measurement pressure, every move logged. (3)
"everything rides on w_rad" — softened by the external anchor
(Δα̂ ≤ 0.01 frozen at Hwang's value).

PLAIN VERDICT: round 9 = SUCCESS as a review round (one new read
executed against pre-stated bands, gates 4/4, three phrasing/identity
gaps closed, zero bars moved) — and the read itself went AGAINST us in
the sharpest honest way: FIFTH-MOVE-LIVE. The α result is now
explicitly conditional on the multiplicity tension resolving toward
the kinematic preference, and the two instruments that decide it are
already first in the queue.

ELI12: The referee asked: "your new answer of 0.74 — is it just what
the model wants, or did your assumptions push it there?" We checked:
zero push — the model picks 0.74 entirely on its own. But his second
worry landed hard: if our own (not-yet-stamped) count of hidden
companion stars is right, and we force the model to use it, the gravity
boost drops to zero. So everything hangs on one question — is that
companion count right, and does each companion really shake the
velocities as hard as our model assumes? The next two instruments in
the queue answer exactly that. Our confidence number stays at 45%
because the rule says only those instruments get to change it.

## Round 10 (2026-07-26): the reviewer's close — the conversion audit (7J-z4)

Opus round 10 = his exit (he owns the D2 misclaim — his third false
execution claim in four rounds — and ends round-by-round auditing; door
left open for 7L and the split). Adjudication, then the executed
instrument.

**His knee mechanics: wrong in the letter, right in substance —
so the missing axis was MEASURED.** His claim: "the knee was one grid
step past the top of the sampled range." Not the mechanism: the anchor
curve varied the CENTER at fixed width (σ = 0.05/0.03), and a smooth
σ = 0.03 anchor centered even at 0.35+ would NOT collapse α (the
likelihood pays the ~30 lnL prior penalty to hold fcomp = 0.1; the
kinematic gain is 135–153). The collapse axis is anchor STRENGTH, and
Part A of 7J-z4 measured it: **α_marg holds 0.67–0.75 at every σ ≥
0.03 and collapses at σ* = 0.02 (both laws), the posterior jumping
0.10 → 0.20 (not 0.35 — the collapse routes through the cell the grid
caveat already flags) and reaching 0.35 only at σ = 0.005.** The
unshipped rate's apparent precision is ~0.015 < σ*: a certificate at
face precision FIRES the fifth move on the flat-q axis. His substance
("the flat-curve claim shipped without the axis that carries the
risk") is adopted as a standing discipline line.

**His two-moments objection: VINDICATED as mechanism, INCOMPLETE as
stated — both halves now measured** ([calcs/stage7j_qmoments.py](calcs/stage7j_qmoments.py),
pre-framed 3cbedbe, amendments A1/A2 logged post-fail pre-quote,
all gates PASS on the amended run; GA0/GA1/GC0 regressions exact):
- His anticorrelation is real and large: hard-flag threshold q_min =
  0.70–0.88 (M_h 0.35–0.75); wfac peaks at q ≈ 0.4–0.6; **detected
  companions carry 0.13–0.17 of the undetected wobble variance**
  (6–8× suppression) — the companions the photometry counts are the
  ones that wobble least, quantitatively.
- His gap: the model's second companion channel — the hidden-mass
  velocity inflation √(1+q/2), TWIN-MAXIMAL (+22% per system at
  q = 1), one-sided, γ-invisible, s-flat (the sq-shaped channel).
  Twins are wobble-invisible but mass-channel-loud.
- The channel attribution (the cubes' own kw axis, which scales only
  the wobble kicks): **the forced-multiplicity rejection is
  WOBBLE-BINDING — Dwob = +314/+306 (doubling the wobble amplitude
  triples the cost 153 → 470)** — so the flexible (wobble) moment is
  the operative one at the margin. Residual attribution AT kw = 0.7
  (the fitted edge) stays open for the arm suite.
- THE JOINT CONVERSION BAND (amendment A2's corrected object):
  fce_joint(wobble) = **[0.10, 0.39]** across π(q) brackets — the
  detection-shaped bracket makes measured host 0.30 kinematically
  equivalent to **0.10–0.12 = the kinematic preference exactly (his
  "~0.12" reproduced)**; the q^-0.5 bracket makes the tension WORSE
  (0.35–0.39). The mass-channel joint is π-STABLE (0.25–0.30;
  completeness and mass-weight co-vary and cancel) but non-binding.
- FORMAL reading by the pre-stated letter: **MIXED** (the original
  band read the det-shape moment ratio, 0.69 > 0.5, without the
  completeness rescaling any π re-attribution implies — a band-design
  miss, logged as A2, the joint reported alongside; no silent
  re-labeling). A1: GB0's dark-companion clause tested a limit the
  shared implementation intentionally clips away (MS-table floor
  0.102 M_sun; ℓ(q=0.1) = 0.001–0.024 printed as a model property);
  replaced by the law's true identities (twin-zero exact, ℓ(1) = 1
  exact, interior maximum).

**CONSEQUENCE — the fifth-move exposure is RE-SHAPED, not resolved:**
S3's forced fcomp ≥ 0.35 forced the FLAT-Q axis; identifying
"measured host 0.30 ≡ the model's 0.35 cell" is exactly the
scalar-passing the review flagged. The exposure's true coordinates
are (rate, precision, π_q) jointly: the fifth move fires iff the
certificate lands tight (σ ≤ 0.02) on the flat-q-EQUIVALENT axis,
which cannot be evaluated without the q-resolved conversion.
**v2c is UPGRADED to v2c-plus by requirement: the certificate must
ship a (q, P)-resolved rate (and the v2b posterior's q-information
extracted), not a scalar.** His "the tension is itself the result"
framing: adopted for Paper 1's spine (first group to measure its own
sample's companion rate AND check it against its own kinematics'
preference), with one correction — his "the likelihood is indifferent
at ΔlnL +0.0 either way" misreads dN(forced) = +0.0 (that is
α-indifference WITHIN the forced world; the likelihood pays 135–153
to avoid that world — the indifference is between laws there, not
between worlds). The tension statement that survives: the photometry
demands a world the kinematics price at 135–153, and whether that
price is real physics or q-weighting is precisely the open
conversion.

His exit note is answered in kind: three false execution claims in
four rounds is his honest self-diagnosis of prose-granularity review,
and the audit function inverting (our gates refuting his assertions)
is real — but rounds 7–10 also contributed the dual-config read, the
anchor-curve re-read, the width-channel push, the q-moments
question, and the Paper-1 spine. Net positive, gratefully logged.

PLAIN VERDICT: round 10 = SUCCESS (two instruments in one stage, all
gates pass on the amended run, two bar-design misses logged openly,
zero bars moved, no credence move — cadence holds at ~45%). The
program's binary headline is now fully coordinatized: α = 0.74/0.70
prior-independent for σ ≥ 0.03, σ* = 0.02, conversion band
[0.10, 0.39], and one required instrument (the q-resolved
certificate) that closes it.

ELI12: The referee's parting shot: "the companions your camera can
see are exactly the ones that barely jiggle the star — so your star
count and your jiggle model are counting different things." We
checked: he's right, 6–8× right. And the model's complaint about
"too many companions" is indeed about jiggling (turning the jiggle
knob up makes it much angrier). So the "your count says 0.30 but the
motion wants 0.10" fight might be no fight at all — if hidden
companions cluster near look-alike twins, 0.30 of them jiggle exactly
like 0.10 of the model's average ones. He missed one thing: twins DO
shake the system another way (extra hidden weight speeds everything
up 22%) — but that's not the channel the model is complaining
through. Verdict: our next instrument must report not just HOW MANY
companions, but WHICH KIND — and until it does, nobody gets to say
whether the boost survives. Confidence stays 45%.

## PRE-REGISTRATION: the v2c round + the degraded-arm suite (2026-07-26, committed before execution)

Queue execution ordered by the user ("do all the next parts"). Order:
v2c-plus certificate → anchored re-read (amendment-11 operative
anchor) → the arm suite. Bars locked here BEFORE any v2c number is
seen.

**7J-z2c (the certificate)**: pre-reg in the script docstring
(commit 4c3fc89; amendments A3/A4 = instrument-level, logged
pre-quote: the GV0 print-rounding reference fix + np.trapezoid API +
the GV0c three-clause metric after the far-tail float-cancellation
false alarm — the synthetic sweep proving seg_int exact is in the
scratchpad record). Ship rule: GV0–GV6 all pass → v2c-cert npz with
the GV7 q-resolved block; any fail → quarantine, LIT16 stands, STOP.

**Anchored re-read (amendment 11, committed BEFORE the anchored run,
commit d437921)**: operative anchor = LANDED-CONV = the
conversion-widened profile lnpi_conv(fc) = max_g lnpi_host(fc/g),
g ∈ conv_band/0.30 = [0.33, 1.30] (bin-7j-qmoments); the face
'LANDED' read is printed as the flat-q CONDITIONAL diagnostic.
Verdict bars UNCHANGED (BOOST-REVIVES α ≥ 0.5 AND dN ≥ +25 /
NO-DETECTION α ≤ 0.3 AND dN ≤ +10 / else AMBIGUOUS-CARRIED),
evaluated at OPER = LANDED-CONV; seed-extension rule unchanged.

**7J-z5 (the degraded-arm suite; implementation AFTER the re-read)**:
per-arm injection validation at the operative anchor per the 7I
standard and the round-5 rule (no arm is validated at the operative
anchor until injected data exist there). THE DEGRADATION (7J-w2
requirement, design locked now): injected skies carry the sky's
confusion channels — (i) sq_true = 0.2 (the demanded width channel);
(ii) fpm_true at the per-law posterior mode (2.1 simple / 2.4 BE);
(iii) companion q-MISMATCH: injected companions drawn from the
v2c-preferred q-law (GV7 winner; fallback twin t=2 if v2c does not
ship) while the fitting model stays flat-q as run — the
model-misspecification injection that prices 7J-w2's 4–30×
informativeness excess honestly.
ARMS (each = one injected dataset + full photow cube + read at OPER;
seed 31; extension to 101 pre-authorized if any bar lands within 0.1
of its edge):
  A (null/false-positive): α_true = 0, fcomp_true = 0.20 → bar:
    recovered α_marg ≤ 0.3, else ARM-FAIL.
  B (simple own-truth): α_true = 0.74, fcomp 0.10 → bar:
    |α̂ − 0.74| ≤ 0.25 AND dN ≥ +10.
  C (BE own-truth): α_true = 0.70 → same bars.
  D (BE discriminator): α_true = 0.40 (the 7J-w
    additive-vs-multiplicative probe) → bar: within 0.25; PROF vs
    MARG offsets reported.
  E (sky re-run, fpm → 3.0): FPM grid extended [1.2 … 2.4, 3.0] on
    the REAL sky, both laws both seeds; identity-point regression:
    the ≤2.4 slice of the new cube reproduces the old cube ≤ 1e-3
    (the GB0w rule); then D3 re-read (does the BE edge release?) +
    the kw-residual attribution at fpm = 3.0.
GATES: GB-regressions for every touched build path; injection
determinism; arms VALIDATE — no physics verdict moves from arms;
the sky verdict machinery stands as pre-registered.

Scout result (Haiku, primary-source, logged with this pre-reg): the
round-10 "3× field expected from triple dynamics" claim is
UNSUPPORTED for general wide binaries and actively contradicted —
Tokovinin 2014 (AJ 147, 87): wide-binary subsystem rates are
field-like ("the presence of a wide companion does not affect
inner-subsystem formation"); the literature's 3× enhancements are
contact-binary (Gao+21) and 2+2-quadruple (Fezenko+22) specific;
the 0.2–50 kAU subsystem-fraction range confirmed unmeasured. ⇒ if
the ~0.3 host rate certificates, it is a SURPRISE against the
Tokovinin trend, raising the stakes on the q-resolution; the
"expected from formation" framing is NOT printed. Flag for a
primary read before any paper use: El-Badry & Rix 2018's ~20%
per-component unresolved-companion characterization (scout-level).

## 7J-z2c + the anchored re-read (2026-07-26): THE CERTIFICATE SHIPS — and the photometry measures the subsystems TWIN-HEAVY

**The certificate (135.8 min, ALL GATES PASS, npz v2c-cert shipped).**
The two v2b blockers closed structurally, not by finer grids: the
twin-spiked q-direction integrated in closed form (segment-erf;
GV0c certifies it at 1.7e-8 relative against a 30001-node reference;
GV0b f=0 identity exact 0.0; GV0 identity-point regression 2973.83 vs
v2b's printed 2973.8), the only remaining numeric axis converged at
GV1 d = 0.010 (vs GZS11's failing 2.08), and profiling λ killed the
start-dispersion direction outright: **top-two gap 0.03** (vs v2b's
failing 2.8), λ* = 1.64 interior with a clean peak (neighbors −26,
extremes −103/−126). **f̂ chain: v1 0.166 → v2b 0.159 → v2c
(exact-cont) 0.159 — six-way stable across three model classes,
three integration schemes, and every basin tried.** Injections sharp
(f=0 truth → 0.000, the kernel-blind acquittal repeating; 0.10 →
0.099; 0.25 → 0.253); postdiction bars all pass. Host conversion
(flat/qalt convention, v2b machinery verbatim): **[0.29, 0.32] peak
0.32 per component.** Amendment trail on the way in (all
instrument-level, each caught by its own STOP or output, logged
pre-quote): A3 (GV0 print-rounding reference + np.trapezoid), A4/A5
(the GV0c metric — far-tail float cancellation is not a math error;
the reference's own discretization is not the erf's error; the
relevance floor is the full mixture, not P00 — third design stuck).

**THE Q-SHAPE DELIVERABLE (GV7) IS A DISCOVERY: the blended-δ
photometry itself measures the inner-subsystem mass-ratio
distribution TWIN-HEAVY — twin t=5 beats flat by 162.0 lnL (t=2:
−19.9 behind t=5; flat −162.0; q^-0.5 −267.4), with f̂ = 0.124
blended per component under the winning law** (host ≈ 0.124/0.53 ≈
0.23 — closer to the field than the flat-convention 0.29–0.32). The
population the photometry demands is exactly the wobble-quiet one
round 10 identified as conversion-flexible: twins are photometrically
LOUD (Δm = 0.75) and kinematically QUIET (photocenter cancellation
exact at q = 1). Reading (joint fit = future work, stated as such):
the twin-convention host ~0.23 carries a flat-q kinematic equivalent
of roughly 0.10–0.15 by the 7J-z4 moments — **consistent with the
kinematic preference fcomp ≈ 0.1: the ×3 multiplicity tension is
DISSOLVING under the measured q-shape, in the direction round 10
predicted.** The Tokovinin surprise softens the same way: the
field-like published subsystem rates were flat-q-convention
inferences; our twin-convention host ~0.23 sits nearer the field
than the 0.29–0.32 headline. FLAG: this q-shape is an in-window
blended-δ shape measurement (0.4 ≲ Δm window sensitivity) — the
low-q completion remains prior-bracketed; conv_band and the
face-invalidity note are stored in the shipped npz itself.

**The anchored re-read (amendment 11 operative; 11b = the OPER
wiring fix, c272bef — the amendment-11 edit added the LANDED-CONV
anchor but omitted the one-line OPER flip, so the first anchored run's
verdict block read the face anchor; caught in the output, fixed to
the d437921 pre-committed spec, re-run):**
- **OPERATIVE (LANDED-CONV): α_marg = 0.74 (simple) / 0.70 (BE) at
  dN = +23.8/+23.2 — VERDICT AMBIGUOUS-CARRIED at the LANDED
  anchor** (no longer LIT-conditional; the extension rule is silent
  at OPER, BE seed gap 0.10 ≤ 0.25; identical to the LIT16 read to
  the last digit — the widened prior and the literature anchor agree
  because the curve is flat). D1 CHANNEL-USED; D2 NOT-CLOSED
  (1.41/1.39); D3 BE still rides fpm = 2.4 (P = 0.95–0.98 → the
  E-arm extension queued).
- **FACE (flat-q at face precision, the diagnostic): THE FIFTH MOVE
  FIRED LIVE — α_marg = 0.00/0.00 (simple), 0.37/0.00 (BE), the
  posterior pinned at fcomp = 0.2** — exactly as σ* = 0.02
  predicted. Quoted as the flat-q conditional whose premise the
  certificate's own q-table rejects at −162.
- 7J-g at LANDED-CONV: tiers unchanged (ABSORBER-LEVEL primary,
  SEPARATION-CONFIRMED α-grade co-read) — now at the landed
  operative anchor.

CREDENCE: HELD ~45%. The cadence rule named v2c + the arm suite as
the deciders; one has landed (and it landed threat-side: the
sharpest Newton-rescue — many companions — lost its q-convention
teeth), the other is pending. No move until the arms report.

PLAIN VERDICTS: 7J-z2c certificate = SUCCESS (the resolution cascade
closed structurally; the rate certificated; the acquittal repeated).
The q-shape deliverable = SUCCESS and the round's discovery
(first mass-ratio-shape measurement for inner subsystems at 0.2–50
kAU; twin-heavy at +142 over the nearest alternative). The anchored
re-read = SUCCESS as machinery (both anchors read, the verdict
landed, one wiring amendment logged); the verdict itself remains
AMBIGUOUS-CARRIED — honest, at 1.2–1.8 below the bar, now with a
landed anchor under it.

ELI12: The star-counter got its certification stamp — every check
passed, and the count (about 1 in 6) has now come out the same six
different ways. But the certificate came with a twist we didn't
order: the flavor question. The hidden companions aren't a random
mix — they're mostly look-alike TWINS of their host stars (the data
prefer that reading by a landslide). Twins are the one kind of
companion that glows brightly in photos but barely wiggles anything
— so the old fight ("your photos say 0.30, the motions want 0.10")
was two people using different units: in wiggle-units, 0.23 worth of
twins IS about 0.10 worth of average companions. The fight is
dissolving. The gravity-boost answer at the properly-cautious anchor:
0.74 strength, 23 points, still two points shy of the official
detection bar — same as before, but now standing on a certified
foundation. Confidence stays 45% until the final validation round
(the "arm suite") reports.

## 7J-z5-E (2026-07-26 night): the fpm → 3.0 extension — THE NOISE EDGE CHASES

Four photow3 cubes built behind GB0w + GB0e identity gates (all
0.00e+00 — the ≤2.4 slice reproduces the operative decider cubes
bit-exactly; the E-arm wiring is clean); GE0 slice-identity
regressions vs the shipped anchored read exact 4/4. THE FINDING:
**the edge does not release — it chases. P(fpm = 3.0) = 0.54
(simple) / 0.97 (BE), the posterior mode AT the new edge in both
laws (D3-EXT STILL-RIDING), and this time it costs the signal:
α_marg 0.74 → 0.68 / 0.70 → 0.67, dN +23.8 → +14.5 / +23.2 →
+16.1 — the Newton gap cedes ~8 lnL to the extra noise freedom.**
Dwob unchanged (+310/+305; the wobble attribution is stable).

Reading: fpm = 3.0 claims Gaia's formal PM errors are underestimated
3× — beyond the published physical ceiling (Lindegren+21 ~1.1–1.4;
the 7J-c diagnostic flagged 2.4 as already outside anything
published). The fit is demanding noise the error model cannot
physically supply = the missing-width-SHAPE signature, the
3E → 6P → sq ladder's next rung (sq is a single global s-flat σ;
the data evidently want structured width). DECISION per the
pre-registered clause (a ride at 3.0 queues an extension DECISION,
not an automatic run): **NO further mechanical extension** — more
nodes chase an unphysical direction. The operative α stays the
pre-registered ≤2.4 anchored read, now CO-QUOTED with the
fpm-conditional band: **α 0.68–0.74, dN +14.5–23.8**. Attribution
→ the A–D mismatch arms (do twin-drawn injected skies reproduce the
edge-chasing?) + a width-shape refinement (queued behind the arms).

PLAIN VERDICT: SUCCESS as an instrument (gates exact, the question
answered decisively) — and the answer is adverse-leaning: the noise
sector is still telling us the width model is incomplete, and it is
willing to spend the boost's significance to say so. The honest
binary headline after tonight: α = 0.74/0.70 at +23.8/+23.2 on the
certificated anchor, with a −8 lnL / −0.06 α exposure to the
noise-model ceiling, and the arm suite as the named adjudicator.

ELI12: We gave the fit permission to blame even more of the spread
on "the telescope's error bars are too small." It took the
permission greedily — all the way past what the telescope's own
calibration papers allow — and paid for it out of the gravity
signal's pocket (23 points → 15). When a model keeps reaching for
an excuse physics forbids, that's not the excuse being right —
that's the model telling you its noise description is still missing
a shape. We stopped feeding it new excuses (rule was pre-agreed),
wrote down both readings, and the validation arms will now test
whether twin companions produce exactly this kind of greed.

## 7J-z5 A–D (2026-07-27): THE ARM SUITE — 4/4 VALIDATED, the chase is real structure

Four twin-mismatch injection arms (pre-reg f6916e0; implementation
bad6770; injections at sq_true = 0.2, per-law fpm_true, twin-t5
companions vs a flat-q fitter, fc/ff 0.10/0.05, truth kw = 1.0;
photow3 grid so the 3.0 node was available; seed 31; ~91 min GPU).
Read at the operative LANDED-CONV anchor, bars as pre-registered:

- **ARM A (null, α = 0, fcomp 0.20): ARM-PASS decisively** — α_marg
  = 0.00 BOTH laws, dN +0.0, PROF 0.00. Under the full mismatch +
  width + noise confusion the machinery manufactures NOTHING. The
  7J-w2 injection-optimism worry is answered at the operative
  anchor: honestly-degraded injections, honest machine.
- **ARM B (simple 0.74): ARM-VALIDATED** — recovered 0.65 (offset
  −0.09), dN +20.8. **ARM C (BE 0.70): ARM-VALIDATED** — 0.64
  (−0.06), dN +38.7. **ARM D (BE 0.40): ARM-VALIDATED** — 0.48
  (+0.08), PROF = MARG (offset 0.00): the 7J-w BE-arm pathology
  does NOT reproduce at operative-curve values. No bar within 0.1
  of its edge → no seed extension required.
- Nuisance recovery sharp in all eight reads: P(sq) = δ(0.2) =
  injected; P(fcomp) at injected; P(fpm) at injected mode.
- Calibration note (single-injection grade, NOT applied as a
  correction; the 3A realization layer covers part): own-truth
  recovery biases are mildly CONSERVATIVE (−0.06…−0.09) — if
  anything the sky's 0.74/0.70 slightly under-reads its truth.
  Also: matched-α injected skies show dN ≈ +21…+39 vs the sky's
  +23 — the sky sits at the low end, consistent with residual
  unmodeled width eating significance (the same object as the
  chase).
- **THE CHASE: UNEXPLAINED — truth-law P(fpm = 3.0) = 0.00/0.00/
  0.02/0.00, 0/4 arms ride the edge.** Twin-mismatch skies with the
  full known confusion set do NOT reproduce the real sky's noise
  hunger. The fpm chase (E-read: sky P(3.0) = 0.54/0.97, ~8 lnL)
  is REAL missing width-SHAPE structure — not companions, not
  twins, not the s-flat sq, not PM-error physics. Per the
  pre-stated band: the width-shape refinement is PROMOTED to the
  named next instrument.

CREDENCE (both named deciders now landed; the freeze releases; no
mechanical map was pre-committed for the arm outcomes — stated
openly, so this is a reasoned move, not a map application):
**anomaly-real ~45% → ~50%.** Up: the era's sharpest alternative
(companions absorb the boost) lost its quantitative case across the
certificate round (twin-heavy = kinematically quiet ≈ the kinematic
preference; the ×3 tension dissolved), AND the measurement machinery
survived adversarial validation 4/4 with a clean null. Capped: the
verdict itself remains AMBIGUOUS-CARRIED (+23.8/+23.2, below the
+25 bar; honest band down to +14.5 under the fpm question), and a
genuine model-incompleteness systematic (width shape) stands
unattributed. Not lower: that systematic EATS significance in the
validated direction (the arms show matched-truth skies scoring
higher than the sky does).

PLAIN VERDICT: SUCCESS — the strongest validation round the binary
program has produced (clean null under honest mismatch, three
own-truth recoveries, sharp nuisance identification), with one
honest adverse residue carried forward by name (the width-shape
object, now the sharpest open binary instrument).

ELI12: We built four fake skies — seeded with exactly the twin
companions, the width fuzz, and the noise the real sky has — and
handed them to our measuring machine without telling it the answers.
It aced all four: the sky with NO gravity boost read exactly zero
(the machine can't be tricked into inventing one), and the skies
with known boosts read them back correctly. One clue didn't
reproduce: the real sky keeps begging for extra noise headroom that
none of our fake skies need — so something about the real spread
still isn't in the model, and hunting its shape is now the top job.
Confidence: 45% → 50% — the "it's all hidden companions" story lost
its numbers this week, and the machine proved honest; but the
detection bar still hasn't been crossed, so only half-way it is.

## Stage 7J-z6 pre-registration (2026-07-27, committed BEFORE execution): THE WIDTH-SHAPE REFINEMENT

The arm suite's CHASE-UNEXPLAINED verdict promoted this instrument: the sky
rides fpm = 3.0 (3x Gaia formal errors, past the Lindegren+21 ceiling ~1.4)
with P = 0.54/0.97 while NO injected sky — twins, sq = 0.2, posterior-mode
noise — puts more than 0.02 there. The sky's noise hunger is real structure
the current width family (global error-proportional Gaussian fpm + s-flat
lognormal sq) cannot express. This stage asks WHAT SHAPE it wants.

**Composition probe (input, recorded before this pre-reg):** the
perspective-residual candidate is CLOSED — 100.0% of the 14,071 ok-pairs
carry at least one Gaia RV (nearby bright sample), and the RV-error-propagated
along-separation residual is <= 0.06 sg0 even at 20–50 kAU (median 0.0008
km/s vs sg0 0.015). The 4R correction is essentially complete; no likelihood
was run. Live candidates: ADDITIVE FLOOR (error-independent jitter), TAIL
(a fraction of pairs with badly underestimated errors), MAGNITUDE-dependence
(the Lindegren inflation axis), or none-of-these (carried).

**Part A — diagnostic (measurement, no verdict; conditional at the sky's
PROF mode cells wr=0.2, fcomp=0.1, sq=0.2, per-law alpha cell, seed 31;
the 5B coordinate-descent lesson applies: conditional previews steer the
contest, they do not decide anything):**
- A1 error-half split: within each s-bin, split pairs at the bin-median
  formal error sg0; data histograms and model noise pools split identically;
  read the per-half fpm profile. LOW-half-carried = the demand is
  error-INDEPENDENT (floor family); HIGH-half-carried = error-CORRELATED
  (tail/inflation family).
- A2 cell attribution: per-(v-row, gamma-column) decomposition of the
  data*log(pp) gain between fpm = 2.4 and 3.0 at the mode cell. Tail-rows
  (v > 1.5) carrying >= 60% = tail family; shoulder-carried = width proper.
- A3 magnitude split: halves by G_faint (the fainter component drives both
  the PM error and the Lindegren inflation). Faint-carried = magnitude-
  shaped inflation (tail-on-faint variant).
- OUTCOME MAP (locks Part B's first mode): A1-LOW -> FLOOR first; A1-HIGH
  or A3-faint-sharp -> TAIL first; ambiguous -> FLOOR first (declared
  default). The second shape runs only if the first fails B1.

**Part B — the shape contest (bars LOCKED NOW; grid VALUES for the shape
axes set from Part A scales and logged as an amendment before Part B runs;
kt = 4 fixed for tail, kt = 6 only as a disclosed sensitivity variant):**
- Modes (new WSHAPE env in stage7j_marginal.py, one at a time):
  FLOOR: per-pair sigma_eff = sqrt((sg0*fpm)^2 + floor^2), floor axis
  appended; TAIL: per-pair fpm_eff = fpm*(1+(kt-1)*[u_t < ft]), ft axis
  appended, u_t a NEW draw appended LAST in build_pop (the amendment-7
  precedent — every earlier draw unchanged). The fpm core grid KEEPS the
  3.0 node (the chase must dissolve voluntarily, not by amputation).
- Gates: GW0 identity — the shape-param = 0 slice must reproduce the
  photow3 cube <= 1e-3 (expect exact 0.0e+0 by stream construction; the
  GB0e precedent); GW1 — the B4 arm's own-truth shape parameter recovered
  within one grid step.
- B1 SHAPE-WIN: profile-max(full shape cube) - profile-max(shape=0 slice)
  >= +8.0 lnL in BOTH laws (one law only = PARTIAL). +8 = the scale the
  sky ceded to the fpm extension; >> the 1-param AIC cost.
- B2 CHASE-DISSOLVED: P(fpm = 3.0) <= 0.10 in the shape-mode marginal at
  the LANDED-CONV anchor, per law.
- B3 alpha-STABILITY: alpha_marg(shape) in [0.55, 0.90] AND dN >= +10 ->
  band-stable (the co-quoted band extends to include the new value).
  alpha_marg < 0.55 OR dN < +10 -> MATERIAL-LOW (the width-shape was
  carrying alpha; reported without band protection). alpha_marg > 0.90 ->
  MATERIAL-HIGH (the mis-shape was suppressing). Named, not protected.
- B4 REPRODUCTION: one arm per law at the winning shape's posterior mode
  (its own alpha_marg, shape params, fcomp, sq), FIT by the global-fpm
  photow3 machinery: P(fpm = 3.0) >= 0.4 = REPRODUCED (the sky showed
  0.54/0.97); 0.1–0.4 = PARTIAL; < 0.1 = FAILED (the shape is not the
  chase's cause even if it fits better).
- VERDICT GRAMMAR: RESOLVED = B1(both) AND B2 AND B4 >= 0.4, with B3
  reported; PARTIAL = exactly one of B1/B2/B4 misses; UNRESOLVED-CARRIED =
  neither shape clears B1 -> the co-quoted band 0.68–0.74 / +14.5–23.8
  stands, the item CLOSES (no grid-extension chase; pre-stated), queue
  proceeds to 7J-d.
- Seeds: contest at seed 31 both laws; the winner re-run at seed 101 for
  stability (reported, not a bar).
- CREDENCE: this is a systematics instrument, not a named decider — NO
  credence move this stage regardless of outcome. A MATERIAL B3 would
  reopen the band and route the question to the next decider (stated now).

**Part A EXECUTED (same day; G0-diag EXACT 0.000 both laws — the
recomputed mode cell reproduces the stored cube value bit-for-bit).
THE FINGERPRINT:**
- A1 (error split): AMBIGUOUS — simple's low-error half edge-rides
  (+0.39) while its high half peaks at 2.4 (−0.63); BE's BOTH halves
  edge-ride (+0.25/+1.74). Per the locked map: FLOOR RUNS FIRST.
- A2 (cell attribution): NOT tail-shaped — tail (v>1.5) share only 0.18
  both laws. The 2.4→3.0 gain lives in the MID-SHOULDER (v 0.7–1.7:
  +23/+21 lnL summed) and the RADIAL gamma column (8°: +9.4/+15.7),
  paid from the v 0.13–0.53 region (−26/−28); carried by the INNER/MID
  s-bins (simple: 2–6 kAU +6.8 of +6.75 total; BE: 0.2–2 kAU +4.8 of
  +2.87) — strongly s-GROWING shapes are dead (consistent with the
  perspective closure).
- A3 (magnitude): INCONSISTENT between laws (simple bright-carried
  +7.38, BE faint peaks 2.1) — no clean Lindegren signature.
- Honest reading logged BEFORE the contest: the fingerprint
  (mid-shoulder + radial-column + inner-bin, tail-poor) looks as much
  like a POPULATION shape (the 4J/4M radial-arm direction) as an
  instrument-noise shape; the contest stays locked to FLOOR/TAIL and
  the fingerprint is booked as a dividend regardless of verdict.

**AMENDMENT 2 (logged after the floor read, BEFORE the validation arm
runs; the tail leg proceeds per the locked flow):** the floor read
returned B1 ONE-LAW (+7.12 simple MISS / +10.61 BE PASS), B2 MISS both
(P(fpm=3.0) = 0.60/0.81 — the chase DEEPENS: floor edge-rides 0.045
km/s WITH fpm still at 3.0; total demanded noise ≈ 4.3× formal =
deeply unphysical; the correction-#4 grid-edge fingerprint on the new
axis), and B3 MATERIAL-LOW both laws: **α_marg → 0.00, dN → +0.0 at
LANDED-CONV — the +23.8/+23.2 detection is exchangeable against one
sub-bar isotropic width axis** (P(ws)=1.0 at the edge; fcomp pinned
0.1, sq 0.2). Interpretation REQUIRES the exchangeability arm (bars
pre-stated NOW, before it runs): the EXISTING 7J-z5 arm-B injection
(simple α=0.74 truth, twin-t5, sq=0.2, fpm=2.1, NO floor in truth)
read by the WSHAPE=floor fitter —
- INFORMATIVE: α_marg(simple arm) ≥ 0.5 AND P(ws=0.045) ≤ 0.5 (a real
  injected boost SURVIVES the floor axis → the sky's collapse would be
  data-driven, and the exposure is a genuine verdict-level finding);
- EATER: α_marg ≤ 0.3 AND P(ws=0.045) ≥ 0.5 (the axis eats its own
  injected boost → the sky collapse is CONSTRUCTIONAL — the floor is a
  boost-degenerate absorber, the collapse says nothing about the sky;
  the exposure is named as a degeneracy, the band stands);
- else AMBIG (reported as-is).
The arm interprets B3; it cannot flip the B1/B2 outcomes or the locked
verdict grammar. GW0 for the arm: its ws=0 slice vs the existing
fullarmb photow3 cube (exact expected).

**AMENDMENT 3 (logged BEFORE the tail/arm results are seen): the PHYS
conditional** — a reported (non-operative) read of every contest cube
under the PHYSICAL noise envelope: fpm ≤ 1.8 (Lindegren+21 inflation
ceiling ~1.4; grid cell above it retained for margin) and ws ≤ 0.015
km/s (the ~0.025–0.03 mas/yr angular-covariance systematic at ≤200
pc). Purpose: locate the α detection's dependence on ALLOWING
unphysical noise — if α rises under the physical envelope, the
operative 0.74/+23.8 was noise-diluted and the exposure runs the
OTHER way; if it collapses there too, the exchange is not
noise-mediated at all. Reported alongside the verdict; the operative
anchor/model does not change in an interim round (cadence rule).
Exchangeability-arm reader block added in the same commit.

**AMENDMENT 1 (grids from Part A scales, logged BEFORE Part B runs):**
floor axis [0, 0.015, 0.030, 0.045] km/s (0.045 ≈ sg0_med·√(3²−1.2²),
the value that mimics fpm=3.0 at the median error); tail axis
[0, 0.03, 0.08, 0.15] with kt=4 fixed. B4/GW1 arm design note: the
width arm injects with MODEL-MATCHED flat-q companions (channel
isolation — the twin-mismatch axis was validated separately by the
7J-z5 suite); truth passed via env WTRUTH, the SAME injection read
twice (global-fpm fitter for B4, shape fitter for GW1).

## Stage 7J-z6 EXECUTED (2026-07-27): THE WIDTH-SHAPE CONTEST — verdict UNRESOLVED-CARRIED, and the NOISE-CEILING PROFILE is the round's real product

All gates exact on every leg (GW0 = 0.00e+00 ×4 — the ws=0 slices
reproduce photow3 bit-for-bit; GB0w/GB0e 0.00e+00; the ws==0 branch
carries the legacy expression verbatim by construction).

**The locked-grammar verdict: UNRESOLVED-CARRIED.** B1: floor +7.12
(simple, MISS) / +10.61 (BE, PASS) = ONE-LAW; tail +0.00/+0.00 = DEAD
(the axis is unused at the profile peak — ws mode 0.0 both laws).
Neither shape clears B1 both laws; B2 never dissolved (P(fpm=3.0)
0.60–0.99 across legs); B4 never ran (no winner). The co-quoted band
0.68–0.74 / +14.5–23.8 STANDS; the item closes per the pre-registered
closure rule — no grid-extension chase into the (already unphysical)
width direction.

**Dividend 1 — THE NOISE-CEILING PROFILE (the round-10 "ship the risk
axis" rule, executed):** the α detection now has a measured
three-point profile along the admitted-noise axis:
- PHYSICAL envelope (fpm ≤ 1.8 Lindegren, ws ≤ 0.015; amendment-3
  conditional, non-operative): **α_marg = 0.80/0.80, dN = +35.2/+32.3**
  (photow3); floor-axis-open variant 0.71/0.72 at +33.5/+30.3; within
  PHYS the fpm posterior still rides its sub-ceiling edge (1.8) —
  the demand object persists at every cap; the cap sets who pays.
- OPERATIVE (flat to fpm 3.0, the co-quoted band): 0.68–0.74 /
  +14.5–23.8.
- UNPHYSICAL CORNER ADMITTED (floor 0.045 + fpm 3.0 ≈ 4.3× formal):
  **α_marg = 0.00, dN = +0.0 both laws** (the floor-leg B3
  MATERIAL-LOW) — the detection is fully exchangeable against
  isotropic noise ~3× beyond the published Gaia error budget.
The operative quote is therefore noise-DILUTED, not noise-protected:
every unphysical width cell admitted drains the Newton margin through
the α=0 row (the tail leg shows the volume mechanism cleanly: an
axis UNUSED at its peak still drains dN +23.8→+10.1 / +23.2→+6.2 by
posterior volume — sub-B1 axes can zero α_marg without ever winning;
B1-vs-B3 is peak-vs-volume). External physical anchors against the
corner: Lindegren+21 (inflation ~1.1–1.4, systematic floor
~0.025–0.03 mas/yr ≈ 0.012–0.024 km/s at ≤200 pc — the corner needs
0.045 km/s ON TOP of 3× inflation) and the census cliff (an isotropic
0.045 km/s floor would blur the 1.65 edge measured sharp in 7I/4J).

**Dividend 2 — THE EXCHANGEABILITY ARM: INFORMATIVE** (pre-stated bar,
amendment 2): the 7J-z5 arm-B sky (injected simple α = 0.74, twin-t5,
sq = 0.2, fpm = 2.1, NO floor in truth) read by the floor fitter
returns **α_marg = 0.60, floor QUIET (P(ws = 0) = 0.58), P(fpm = 3.0)
= 0.00** (BE cross-read 0.60/+11.7, floor quiet). The floor axis does
NOT eat real boosts — a matched-width sky keeps its α and rejects the
floor. The sky's floor-corner preference is therefore DATA-DRIVEN:
the real sky contains ~7–11 lnL of width the model lacks, and only
the real sky offers it. (Own-truth under-read 0.74→0.60 consistent
with the arm-suite conservative bias, single-injection grade.)

**Dividend 3 — THE SHAPE EXCLUSION SET + the fingerprint:** the
residual width object is now measured to be NOT floor-shaped (B1
one-law, edge-riding into unphysical amplitude), NOT tail-shaped
(+0.0 — kt=4 fraction axis wholly unused), NOT fpm-shaped (the chase
runs past the physical ceiling), NOT sq-shaped (sq = 0.2 already in,
demand persists), and NOT perspective (100% RV coverage, composition
probe). Its Part-A fingerprint: mid-shoulder ṽ 0.7–1.7, the RADIAL
γ = 8° column, inner/mid s-bins, spread across error and magnitude
halves. Five noise costumes excluded while the location is
population-like → **named successor hypothesis (booked, not opened):
the inner-bin eccentricity/radial-orbit sector — the fit's e-shape
freedom is nearly frozen (eta grid = 2 values; w_rad external), and
a radial-population surplus at mid-ṽ is exactly what the sky keeps
trying to buy as noise.** DR4-era or a dedicated stage; the freeze
(§0) is untouched (this is an error-model/population refinement, not
a new gravity function).

No credence move (pre-stated: not a named decider). The operative
anchor and band are unchanged; the PHYS conditional is REPORTED
alongside them (promoting it to operative would need a decider-grade
round — flagged as a candidate question for the next external
review, alongside the exposure it inverts).

**Plain verdict: SUCCESS as an exclusion-and-profile instrument**
(the promoted question — what shape does the sky's noise hunger
want? — is answered: NO physical noise shape; the risk axis is now a
measured curve; the machinery re-validated under the new axis), with
the object itself still unidentified — its identification is the
named successor hypothesis, not this stage's claim. NEEDS REFINEMENT
applies only to that successor.

ELI12: The sky keeps asking for more measurement-blur than Gaia's
own manual allows. We offered it two shaped blurs — a constant
"everyone shakes a little" blur and a "a few stars shake a lot"
blur. It ignored the second, grabbed the first only when we let the
blur grow into impossible territory, and — tellingly — once the blur
got impossible, the gravity signal vanished into it. So we ran the
honesty check: feed a fake sky with a KNOWN gravity boost to the
same machinery — it kept the boost and refused the blur. Real skies
only. And when we cap the blur at what Gaia's manual actually
permits, the gravity signal comes back STRONGER than our official
number. Bottom line: our official number is the cautious one; the
mystery extra spread is real, it isn't any kind of instrument noise
we can draw, and its footprint (medium speeds, radial directions,
closer pairs) smells like orbits we haven't drawn quite right — the
next thing to model, written down before we model it.

## Stage 7J-d pre-registration (2026-07-27, committed BEFORE execution): THE FUNCTION-CONTEST RE-RUN UNDER THE LANDED POSTERIOR (the DATA-VOTE unsuspension)

Since 7J the world table carries the suspension banner: every binfn-*
contest was run inside the multiplicity fence (f ≤ 0.1) without the
width channel. 7J-d re-scores the EXISTING world-table function set
under the landed width-complete model. The freeze (PREDICTIONS.md §0)
applies: NO new functions — a re-scoring of declared forms only.

**MODEL (resolved by the 7J-z6 verdict):** the photow3 9-dim grid
(fpm to 3.0, sq, kw, fcomp) at the LANDED-CONV anchor — no shape
axis (UNRESOLVED-CARRIED left the operative model unchanged); the
noise-ceiling exposure banner carries to every row.

**FUNCTIONS (16 new cubes; all tables grid-identical to the
reference, fingerprints printed at load):** p065, gm, rb1–rb4
(= F1–F4), boot, amb (AMB), resn (RESN), dwf (drive-weighted), and
the λ-family lam000/025/050/075/100/125 — **the λ leg IS the landed
ĉ₁ budget**, superseding the promised six-seed fenced budget
(supersession reasoning: two realizations under the landed model
outrank six under the retired fenced model; the 4X ledger row flips
to CO-QUOTED with a pointer on completion). simple/BE: the existing
photow3 cubes ARE their 7J-d cubes (reused, not rebuilt).
EXCLUDED-DISCLOSED: sbe (its solver grid lacks the 1e-4..1e-6 deep
end — re-interpolation would extrapolate; TRIGGER pre-stated: if any
veto MOVES by > 8 lnL between fenced and landed models, sbe gets a
dedicated deep-end-regenerated run before the unsuspension is
declared) and qcl (grid mismatch; its +556 annihilation is two
orders beyond fence-sensitivity). OUT OF SCOPE: b025/b075 (the 5R
β-bound stands as a fenced-model result; the β-direction is readable
through gm), mi_t (not a table function; the 4L/7H formulation
results are not fence results — DR4 item).

**GATES:** G0-alpha0 — every function cube's α = 0 row must equal
the simple sky cube's α = 0 row to ≤ 1e-9 (Newton is function-blind;
the newt_cache precedent), ABORT on fail; GD1 — the lam100 table is
bit-identical to BE (verified: array_equal True), so the lam100 cube
must reproduce the BE photow3 cube to ≤ 1e-9 = the free end-to-end
FUNC-path regression; per-function B(y=1) fingerprint printed at
load (manual cross-check vs the known ν(1) ladder); GB0w/GB0e print
SKIPPED-disclosed (no photo/photow references exist for function
laws — G0-alpha0 + GD1 substitute).

**METRICS per function (per seed):** α_marg + dN(Newton) at
LANDED-CONV; PROF α̂ + interiority (edge-riding = shape rejection,
the 6-series standard); the VOTE = ln-evidence(fn) − ln-evidence(BE)
(logsumexp over the full cube + eta prior + LANDED-CONV fcomp prior
+ flat α grid — identical priors for every function); nuisance
posteriors (does any function un-ride the fpm edge?).

**BARS (locked now):**
- VETO grammar: a previously vetoed member (rb4/F4, boot, resn, dwf)
  STAYS vetoed iff (PROF α̂ edge-riding) OR (vote ≤ −8); a flip is
  named and its world-table row updated with both model conditions
  quoted.
- LEAN grammar: |vote| < 5 = tie-grade (the 3P SE precedent); 5–15 =
  lean; > 15 = strong lean; no detection language.
- λ-family: the landed binary ĉ₁ = λ̂/2 from the evidence profile
  over the six λ nodes, quoted with the lam000 (c₁ = 0) rejection
  margin and the per-seed peak spread; two-realization grade stated.
- UNSUSPENSION: declared when the full set has seed-31 rows AND the
  no-veto set {p065, gm, rb1, rb2, rb3, amb} + λ-family has seed-101
  rows with sign-stable votes; the world-table banner then flips to
  the landed-vote form. Seed-31-only members carry "single-seed" in
  their rows.
- No credence move (not a named decider); the binary-vote statement
  feeds the world table, not the credence line.

## Stage 7K-a pre-registration (2026-07-27, committed BEFORE execution): THE FORWARD MEDIAN — does the width-complete Newton-best cell reproduce the model-free anchor?

The 2C statistic is model-light by construction: boost = median
ṽ(6–30 kAU) / median ṽ(0.2–2 kAU); Newton predicts 1.000 exactly
(ṽ's distribution is s-independent under Newton — that is the
normalization's purpose). The corrected anchor: 1.078 (CI
1.052–1.103; 4Q, correction #12). The absorbers are NOT ratio-blind:
companion wobble in ṽ units grows as w·√(s/M) and the noise floor as
sg0·√s — both rise toward the wide bin. 7K-a MEASURES how much of
the 1.078 the landed Newton-best cell produces.

**INSTRUMENT:** self-contained forward (the marginal machinery's
population, companions incl. kw, noise, sq; own copy — the batch
script stays frozen), seed 31 primary / 101 stability; the model
median ratio computed on the model pairs DIRECTLY (no histogram
coarsening), both naturally-weighted and s-REWEIGHTED to the data's
within-bin s-distribution (10 sub-bins per 2C bin; the reweighted
number is the bar-carrying one). Cells read from the photow3 cubes +
LANDED-CONV prior in-script:
- NEWTON-BEST: argmax of the α = 0 slice (nuisances free);
- SANITY legs: the α = 0.5 and α = 1.0 cells (the operative 0.74/0.70
  bracket) — expectation band R ∈ [1.04, 1.12]; outside → instrument
  flag, investigate before quoting.
**GATE G0-7ka (amended BEFORE execution — the baseline-noise clause):**
two variants: (i) NOISE-OFF (α = 0, fcomp = 0, sq = 0, fpm term
removed entirely) ⇒ R ∈ [0.995, 1.005] — the exact Newton-flat
identity; (ii) BASELINE (fpm = 1.2, others off) ⇒ R ∈ [0.995, 1.015]
— baseline measurement noise itself grows as sg0·√s in ṽ units and
is PREDICTED to lift R by ~0.005–0.010; the gap between (i) and (ii)
is the noise contribution to the ratio, reported (it also bounds how
much of the data's 1.078 baseline noise alone explains).
**BARS (locked now, on the s-reweighted Newton-best R_N, seed 31; 101
must not flip the category or the result is AMBIG-quoted):**
- R_N ≥ 1.052 (the data CI low edge) → MEDIAN-ABSORBED: the
  model-free anchor falls as independent evidence — the absorber
  budget covers it; the anchor's ledger rows flip to
  ABSORBER-CONDITIONAL.
- R_N ≤ 1.030 → MEDIAN-SURVIVES: the anchor is the surviving
  model-free anomaly; the gap (1.078 − R_N) is quoted as the
  unabsorbed excess with the data CI.
- 1.030 < R_N < 1.052 → GRAY (partial absorption, quoted as such;
  no status flips).
No credence move (not a named decider); the result feeds the anchor
rows and Paper 1's model-light chapter.

**AMENDMENT (G0 design error, logged pre-quote after the gate fired
FAIL on first run):** G0-i/ii read 0.969–0.992 — the identity
assumption was WRONG as designed: the FITTED population's
eccentricity distribution runs with separation (the v7 al ramp
0.6→1.0→eta over 100–1000 AU), and more-eccentric wide orbits have
lower median ṽ, so model-Newton's 2C ratio is genuinely ≈ 0.98, not
1.000 — "Newton predicts 1.000" holds only for scale-free
populations (the 2C narrative assumption, now measured). The wiring
identity is re-housed in G0-iii: a SCALE-FREE control (al ≡ 1.0,
noise off) must give R = 1.000 ± 0.005. The pre-registered VERDICT
bars are unaffected (they compare R_N directly to the data CI, no
1.000 reference). DIVIDEND recorded: relative to the model's own
Newton baseline (0.986/0.992 with noise), the data's 1.078 sits
+0.092/+0.086 high — the model-referenced excess is LARGER than the
naive 1.078/1.000 reading.

## Stage 7K-a EXECUTED (2026-07-27): THE FORWARD MEDIAN — verdict GRAY, seed-stable; the anchor is half-absorbed and survives as an excess

Instrument validated on the third gate design: **G0-iv interior
scale-free control PASSES both seeds (R = 1.005 / 1.002)** — the
wiring identity holds where it should; the G0-i/iii "failures" are
MEASURED MODEL PROPERTIES, not errors: (a) the fitted eccentricity
run (v7 al ramp) makes model-Newton's 2C ratio ≈ 0.98 — the "Newton
predicts 1.000" narrative is scale-free-conditional; (b) the a_s
prior's 0.15/60 kAU truncation distorts the edge-adjacent 2C bins.
Data-side regression: the instrument reproduces the 4Q anchor
exactly (1.078 from the corrected catalog).

**VERDICT (pre-registered bars, s-reweighted): GRAY both seeds** —
the landed Newton-best cell (fcomp = 0.1, fpm = 3.0, kw = 0.7,
sq = 0.2, its own eta/wr) forward-produces **R_N = 1.033 (seed 31) /
1.043 (seed 101)** vs the anchor 1.078 (CI 1.052–1.103): above the
SURVIVES bar (1.030), below the data CI low edge (1.052). No status
flips (per bars); the anchor rows gain the measured note.
- Model-referenced accounting: the Newton+noise baseline at the
  fitted e-run is 0.986/0.992, so the data sit +0.092/+0.086 above
  model-Newton and the absorbers lift +0.047/+0.051 — **the landed
  absorber budget covers ≈ 52–59% of the model-referenced gap;
  unabsorbed excess +0.045/+0.035.**
- The boost cells OVER-produce the ratio (prof cells: 1.127/1.119
  seed 31, 1.119/1.114 seed 101 — the seed-31 simple leg 0.007 above
  the [1.04, 1.12] sanity band = instrument flag per pre-reg,
  investigated: the fpm = 3.0 chase lifts R; the fit trades the 2C
  ratio against the 2D histograms it actually fits — the known
  residual-width story's median face). **The data's 1.078 sits
  almost equidistant between Newton-best (−0.035…−0.045) and the
  boost cells (+0.036…+0.049): the 2C ratio ALONE no longer
  discriminates at the landed absorber budget — the discrimination
  lives in the full 2D likelihood (where Newton loses +14.5–23.8).**
  The anchor's honest role: a model-light EXCESS over
  Newton+absorbers (half-covered), no longer a standalone
  Newton-killer.

Plain verdict: SUCCESS — the question the queue asked ("does
Newton+companions reproduce the median?") is answered with a
measured number and a category (GRAY: half of it, not all of it,
seed-stable), the instrument's three gate iterations are logged
pre-quote, and the anchor's evidential role is now stated precisely.

ELI12: Our cleanest single number says wide star-pairs move 7.8%
faster than gravity-as-usual predicts. We asked: could the boring
stuff we now know about — hidden third stars, measurement noise,
spread — fake that number with NO new gravity? Answer: it fakes
about half (3.3–4.3%), and can't reach the real 7.8%. Meanwhile the
new-gravity versions overshoot it slightly. So this one number
alone can't pick the winner anymore — the full detailed comparison
still can (and there gravity-as-usual keeps losing) — but the
number stays honest evidence that SOMETHING beyond the boring stuff
is there. Also, two of our test alarms rang and turned out to be
telling us true things about the model, not mistakes — logged, as
always, before we looked at the answer.

## Stage 7K-b pre-registration (2026-07-27, committed BEFORE execution): THE CENSUS LEAKAGE NULL AT THE LANDED CELL + the wobble-lite per-pair bounds

Correction #15 flagged the ceiling census's leakage null as
inheriting the companion audit. The landed-cell version is CLEANER
than 4J's analytic leak: the forward Newton-best sky at the landed
posterior mode ALREADY CONTAINS the companion wobble + hidden-mass
machinery, the noise at its demanded fpm, and the smear — so its
band occupancy IS the companion-marginalized null.

**INSTRUMENT (stage7kb_census.py, the 7K-a machinery):** per seed
(31/101), forward skies at (a) the landed NEWTON-BEST cell, (b)
sensitivity variants fcomp = 0.2 and kw = 1.4 (posterior-adjacent
cells), (c) the per-law PROF boost cells (descriptive consistency —
the boost model should populate the band and starve [1.67, 2.2)).
Band = ṽ ∈ [1.414, 1.67) ∧ γ ≥ 75° computed exactly as lnL_point
(γ pre-sq, ṽ post-sq, keep-cut applied); μ = Σ_bins N_data(bin) ×
f_band(bin) over the v7 s-bins (9950/2684/1223/214). Observed:
**n_obs = 9 (corrected convention, the operative census; raw-sn3 11
reported alongside).** The data-side census carries quality cuts the
model does not model — an UNDERCOUNT direction, making P(≥ n_obs)
conservative (stated).
**WOBBLE-LITE per-pair bounds (no photometry needed):** for each of
the 9 pairs, the UNCONDITIONAL maximum companion wobble shift:
Δṽ_max = kw_max(1.4) · max_q wfac(q) · v_orb(P = 17.8 yr) / 4.74047
/ vc_pair (the S·v_orb product peaks at the Gaia-baseline period;
a_in ≈ 8 AU ≪ every validity cap). SAFE iff Δṽ_max < (ṽ_corr −
1.414). The HIDDEN-MASS channel (boost √(1+m_h/M) — twin-maximal
and photometrically LOUD at δ = 0.75 mag) is explicitly NOT boundable
without per-pair photometry — disclosed; the photometric per-pair leg
(7K-b2) fires ONLY if the population null breaks.
**BARS (locked now, on the corrected n_obs = 9, seed 31; seed 101
and the sensitivity variants must not change the category, else the
weakest category is quoted):**
- NULL-INTACT: P(≥ 9 | μ_Newton-best) ≤ 1e-4 → the census stands
  companion-audited at the population level (the #15 inheritance
  resolved).
- NULL-BROKEN: P ≥ 0.01 → the landed absorbers explain the band;
  census rows flip to ABSORBER-CONDITIONAL; 7K-b2 fires.
- else GRAY (quoted; no flips).
- Wobble-lite: if ≥ 7/9 pairs SAFE, the sentence "the wobble channel
  cannot build the census pair-by-pair" is licensed (informational).
No credence move (not a named decider).

## Stage 7K-b EXECUTED (2026-07-27): verdict NULL-BROKEN by the letter — and the breaking cell SELF-REFUTES on the cliff; the census's power relocates to the (band, cliff) PAIR

**Bar verdict (locked): NULL-BROKEN, both seeds, all pre-registered
variants** — the landed Newton-best cell predicts μ_band = 26–38 vs
the observed 9 (P(≥9) ≈ 1): under the landed absorbers the raw band
COUNT carries no Newton-rejection power. But the same cell predicts
μ = 17–32 in the overshoot [1.67, 2.2) where the data show **2** —
the pre-registered descriptive statistic reads P(≤2) = 7.1e-6 …
2.9e-11: **the cell that absorbs the band is refuted by the cliff.**

**Attribution (post-hoc, labeled):** μ_band is fpm-INSENSITIVE at
sq = 0.2 (24–26 at physical fpm 1.5/1.8) — roughly half the flood is
the γ-blind lognormal smear lifting the perpendicular ṽ ≈ 1.1–1.3
population; the rest is landed noise scattering the near-parabolic
pericenter pile (w_rad = 0.2; pericenter velocity is PERPENDICULAR —
the pile sits at the √2 edge by geometry) plus companions. At sq = 0
+ physical noise, Newton still floods (μ_band 13–14, overshoot ≈ 13)
— and the BOOST cells at sq = 0 likewise (band 13–15, overshoot
13.5–14, P(≤2) ≈ 1e-4): **NO tried configuration reproduces the
observed 9-with-a-cliff-of-2; the cliff is sharper than any noise
≥ 1.5× formal allows.** The only regime that fits the cliff is
formal noise — exactly what 4J measured (cliff at the boosted edge
1.65, P = 0.62/0.91) — and at formal noise the band is unleakable
(the original 4J null, μ ≈ 1, P = 3.8e-9).

**THE RELOCATION (the stage's product):** the census is
SELF-DEFENDING as a pair: the cliff bounds the tail noise to
≈ formal; at formal noise the 9 band pairs cannot be leakage. You
cannot simultaneously hold "noise big enough to fake the band" and
"the observed cliff". The raw count retires as a standalone null
(ABSORBER-CONDITIONAL); **the (band = 9, cliff = 2) pair becomes
the operative census statistic**, and it rejects every landed-model
configuration tried (min P ≈ 9e-5) — a THIRD independent channel
(after the noise-ceiling profile and Lindegren) testifying that the
unmodeled width object is mid-shoulder-LOCALIZED and does not reach
the tail: the global lognormal smear and the fpm inflation both
over-extend into a region the data keep clean. Convergent with the
7J-z6 fingerprint and the population-sector successor hypothesis.

**Instrument lessons (logged):** the wobble-lite unconditional bound
is VACUOUS as designed (a P ≈ 18 yr companion's wobble is 8–30× vc —
such pairs leave the catalog; the bound ignores catalog survival;
0/9 SAFE is NOT quotable as vulnerability) — superseded by the
population null in the same stage. 7K-b2's firing condition was met
by the letter, but the measured mechanism is not companions —
smear + noise + pericenter geometry — so the per-pair photometric
leg is NOT the informative follow-up; it FOLDS into the
population-sector successor (no separate stage).

Plain verdict: SUCCESS — the census audit the fifth correction
demanded is done; the count falls, the pair stands stronger than
the count ever was, and the tail region now actively constrains the
width object's shape. No credence move (not a named decider).

ELI12: Nine star-pairs fly faster than plain gravity's speed limit,
in the "sideways" direction where faking is hardest. We asked our
newly-humbled model: could your fuzz have pushed ordinary pairs
past the limit? It said "sure — I'd push 26 past it… and 17 WAY
past it." But the sky shows only 9 a little past and just 2 way
past — a sharp cliff. Fuzz that could fake the nine would smear
away the cliff; the cliff is there, so the fuzz isn't. The nine
defend themselves by standing next to a cliff. (Also: one of our
back-up checks turned out to be a dud by design — noted so nobody
trusts it later.)

## Stage 7J-d EXECUTED (2026-07-27): THE UNSUSPENSION — and the honest vote is NO FUNCTION DISCRIMINATION; Newton stays function-robustly dead

Gates immaculate: G0-alpha0 = 0.00e+00 on every function cube
(Newton is function-blind — the wiring identity) and **GD1 = the
lam100 cube reproduces the BE cube BIT-FOR-BIT on both seeds** (the
free end-to-end FUNC-path regression the lam100 ≡ BE table identity
bought us). 16 functions × seed 31 + 12 × seed 101, ~4 h GPU
each batch.

**THE VERDICT: at the LANDED-CONV anchor, the binary function
contest is a WASH.** All 16 functions sit within ±8 lnL of BE
(seed 31: +6.15 lam075 … −5.35 resn; seed 101: +7.66 p065 … −6.08
simple) with realization scatter of the same size. Per the
pre-registered letter, all four fenced vetoes FLIP (rb4, boot,
resn, dwf now indistinguishable from BE) — read correctly, this is
**instrument resolution loss, not rehabilitation**: the same
width/companion absorbers that ate half the median (7K-a) eat
essentially all of the ν-shape differences. The fenced vetoes stand
as fenced-model results, quoted with both conditions. Sign
stability: the no-veto set is 6/6 stable (p065/amb positive,
gm/rb1-rb3 negative — all tie-grade); lam125 and simple SIGN-FLIP
across seeds (named). simple's fenced −12.6 lean over BE degrades
to no-verdict — **the 3P "no ν-family discrimination" is now
generalized to all 16 under the landed model. SPARC remains the
program's function discriminator** (the 4E lesson, completed).

**THE λ LEG: the landed binary ĉ₁ = NO CONSTRAINT.** The evidence
profile over λ is flat within realization scatter (peak c₁ = 0.375
at seed 31 with rejection margin +3.4; peak c₁ = 0.000 at seed 101
with margin +0.0 — both inside the ±5 tie band). The 4X fenced
result (c₁ = 0 rejected at ΔlnL ≈ 20/seed) was a FENCED-MODEL
result; ledger row bin-c1 → CO-QUOTED. **The §6.4 "two systems
read the same dial" claim now rests on the galaxy leg alone; the
binary leg is fence-conditional.** (This supersedes the promised
six-seed fenced budget by construction — more seeds of a retired
model answer nothing.)

**THE SURVIVING VOTE: Newton is function-robustly dead.** Every
function × seed row rejects Newton at dN ≥ +7.9 (32 rows; seed-mean
simple +14.5 / BE +16.1 on the extended grid = exactly the
co-quoted band). The unsuspended world-table binary column reads:
"votes against Newton at the band grade; agnostic among all 16
modified laws." P(fpm = 3.0) ≈ 1.0 on every function row (one
exception: simple seed 101 at 0.09) — **the noise chase is
function-blind: the width object sits upstream of the ν-choice**,
consistent with 7J-z6/7K-b.

sbe/qcl stayed excluded-disclosed; the pre-registered sbe trigger
did NOT fire (no veto moved by > 8 — the flips are compressions
toward zero, not reversals).

Plain verdict: SUCCESS — the unsuspension is executed and the
honest answer ("no discrimination") is itself the finding; the
binary program's evidential structure is now cleanly layered:
model-light excess (7K-a, half-unabsorbed) + the (band, cliff) pair
(7K-b) + the function-robust Newton rejection (7J-d) + the
α measurement conditional on its noise-ceiling curve (7J-z6) — and
NO binary claim about WHICH modified law. NEEDS REFINEMENT only in
the pre-named direction (the population-sector width object, whose
resolution would restore discrimination power).

ELI12: We re-ran the beauty contest between all sixteen gravity
formulas, this time letting the boring-stuff dials (hidden stars,
noise, spread) turn freely. Result: the judges can no longer tell
the formulas apart — the boring dials can imitate the differences
between them. Two things survive: plain gravity STILL loses to
every single formula, and the galaxies (where the dials don't
reach) remain the place where formulas get told apart. Also, the
earlier "the binaries measured the ½ knob" claim gets demoted
honestly: that reading needed the old fenced dials; with free
dials the binaries can't read that knob at all.

## Stage 7L pre-registration (2026-07-27, committed BEFORE execution): THE COOKSON SELECTION — our pipeline on their sample, their statistic on our model

Cookson et al. (2602.24035): 1,421 RV-clean pairs, median-ṽ flatness
vs r_sky/r_M, no 20% MOND step (Δχ² ≈ 12–16 ≈ 2.7σ, ~1500×). §7.4
engaged the checklist and argued (c) qualitatively — the measured
e-trend SUPPRESSES their step. 7L quantifies all of it on one
pipeline.

**THE PROXY MASK (their cuts on our EDR3 catalog; proxy-grade
DECLARED — their catalog construction differs):** d < 130 pc (plx >
7.6923), plx S/N > 40 both, RUWE < 1.25 both, both components with
finite RVs and |ΔRV| < 10 km/s, s ∈ 1–30 kAU, velocity uncertainty
< 0.1 in ṽ units, corrected ṽ < 2.5, both colors 0.5 < BP−RP < 3.5,
no overluminous component (the 3J ridge criterion, δ ≥ −0.4 both —
their HR filter's proxy). N reported against their 1,421; if N is
off by more than 2× the stage downgrades to qualitative
(PROXY-FAIL, pre-stated). ΔRV/RUWE cuts are data-side only (the
forward model carries no RV/RUWE channel — neutral by construction,
disclosed).

**THE LADDER:**
- L1 (2D likelihood on the cook sample, landed model, LANDED-CONV
  anchor, seeds 31+101): CONSISTENCY BAR — the cook α-marginal must
  not reject the operative band: d = lm_max − max(lm at the α = 0.5
  and 1.0 cells) ≤ 4 → CONSISTENT; > 8 → TENSION-NAMED (their
  cleaning genuinely removes our signal's carriers); else GRAY.
  dN(Newton) reported with the N-scaling expectation (~N/14071 ×
  the full-sample band ≈ +1.5–2.5) — explicitly NO detection claim
  on ~1.4k pairs, power stated first.
- L2 (descriptive): ṽ-only vs 2D on the cook sample (the 7J-g
  phantom structure on their information content).
- L3 (THEIR statistic on OUR model): step = median ṽ(r_sky/r_M ∈
  [1, 3.1]) / median ṽ([0.05, 0.5]), computed with their ṽ < 2.5
  ceiling on (i) our data under the cook mask (bootstrap CI), (ii)
  the landed operative boost cells' forwards, (iii) the landed
  Newton-best forwards. The 7K-a lessons are built in UP FRONT:
  G0-7L = the scale-free control (al ≡ 1) on the same statistic
  must give 1.000 ± 0.010 (wiring); the fitted-e-run model-Newton
  baseline is the honest zero point and is reported FIRST (the
  e-run is s-dependent and r_M mixes M with s — expect baseline ≠ 1).
  BARS: RECONCILED-BY-DESIGN if the boost cells' predicted step ≤
  1.11 (their ≈1.5σ sensitivity: 1σ ≈ 7.4% of the 20% step) AND the
  data step is consistent with the prediction within the bootstrap
  CI; GENUINE-TENSION if predicted ≥ 1.15 AND data step ≤ 1.05;
  else GRAY (quoted).
No credence move (not a named decider). The stage feeds §7.4 with
numbers in place of (c)'s qualitative argument, and Paper 1's
methods chapter (the modeling-philosophy fork, quantified).

## Stage 7L EXECUTED (2026-07-27): THE COOKSON QUANTIFICATION — their flatness replicates, their power against the landed model is ~1.4σ, and the fork is now arithmetic

Proxy mask: N = 1194 vs their 1,421 (16% — well inside the 2× bar;
cut ladder printed: the |ΔRV| < 10 cut halves the both-RV sample,
s ∈ 1–30 + d < 130 pc remove our deep anchor entirely).

**L3 (their statistic on our model; G0-7L scale-free control PASS
0.995/1.002; e-run baseline ≈ 1.00 on THIS statistic — their r/r_M
range excludes the bins where the e-run bites):**
- **Their flatness REPLICATES on our own catalog**: data step =
  0.985 (68% CI 0.909–1.060; N = 795/148 in the two r/r_M bins).
- **The landed boost cells predict step = 1.092–1.101, NOT 1.20**:
  the measured e-trend, the absorbers, their ṽ < 2.5 ceiling and
  their range cut the naive MOND step in half. Their 2.7σ against
  the 20% step is ≈ 1.3–1.4σ against the landed prediction.
- Newton-best predicts 1.029–1.031; boost-vs-Newton separation on
  this statistic at this N ≈ 0.9σ_boot — **their instrument cannot
  distinguish the landed boost from Newton at their sample size.**
- Verdict by the locked letter: GRAY (prediction ≤ 1.11 ✓ but just
  outside the data's 68% CI at ~1.5σ_boot; the ≥ 1.15 tension bar
  not reached). The data lean 0.6σ Newton-side of Newton-best —
  noise-grade.
**L1 (the 2D likelihood on their selection, landed model,
LANDED-CONV): CONSISTENT in ALL FOUR law × seed reads** (d_op =
1.2–3.3, bar ≤ 4): the cook sample does not reject the operative α
— it cannot distinguish 0 from 0.7 (α = 2 IS excluded at −25…−39,
so the sample has power, just not at the level that matters);
marginal dN = +0.0–0.1 (the N-scaling expectation was +1.2–2.0 —
consistent, sub-noise). No signal-carrier-removal signature beyond
the range restriction. **L2**: ṽ-only ≈ 2D at this N (nothing to
phantom on 1.2k pairs).

**THE FORK, QUANTIFIED (feeds §7.4):** their cleaning discards 92%
of the joint sample — including the entire deep anchor — and what
remains carries ~1σ of α-information on their statistic and ~3 lnL
on ours. "No MOND evidence on 1,421 RV-clean pairs" and "Newton
excluded at +14.5–23.8 on 14,071" are ARITHMETICALLY COMPATIBLE
statements; the disagreement is entirely about whether the removed
92% is data or contamination — the modeling-philosophy fork of
§7.4(c), now with numbers in place of the qualitative suppression
argument.

Plain verdict: SUCCESS — both their headline observation and its
compatibility with ours are now reproduced inside one pipeline;
no bar fired against either side. No credence move (not a named
decider).

ELI12: The other team cleaned the star-pair list down to the 1,400
most pristine pairs and saw no gravity boost. We rebuilt their
exact cleaned list from our data — and saw the same flat nothing
they did. Then we asked our model: what SHOULD the boost look like
on that cleaned list? Answer: about half as big as what they tested
for, and their list is too small to see even that. So both camps'
numbers are true at once: their pristine subset genuinely shows
nothing (it can't), and the full sample genuinely shows something.
The real argument is whether the 92% they threw away is treasure
or trash — and that's now a question with numbers on both sides.

## Review round 11 (2026-07-28): the Opus STRUCTURAL READ — adopted with four corrections

His verdict: architecture right, spine right, protect the
conditionals. Adjudication:

**ADOPTED (1) — the noise-ceiling curve as the primary form, PHYS
not promoted.** His structural argument is right and we sharpen it:
his sub-claim "if the object is astrophysical, the physical envelope
is correct and 0.80/+35 is the honest number" is TOO STRONG — with
the width object UNMODELED, BOTH endpoints are biased brackets: the
flat-to-3.0 branch lets noise absorb object AND boost (biased low);
the physical-cap branch forces the object into the remaining
channels, where masquerade-as-boost is unbounded a priori (the
radial-column fingerprint sits where boost signatures live). Under
the noise identification the low branch is right; under the
astrophysical identification the honest number awaits the MODELED
object — the curve brackets the noise axis and is silent on the
population axis. CORRECTION #A to his text; his conclusion (report
the curve, mark both envelopes, identification open, the successor
collapses it to a point) is STRENGTHENED by it. → split plan + §6.3
form updated.

**ADOPTED-CORRECTED (2) — the conditional abstract.** The form is
right (no bare "α = 0.74" in the abstract; AMBIGUOUS-CARRIED demands
the conditional). His range sentence is WRONG: "runs from 0.0 to 0.8
depending on assumptions each of which is separately defensible" —
the 0-reaching configurations are NOT separately defensible; they are
the configurations the data reject on ADJACENT statistics: the noise
corner (floor 0.045 + fpm 3.0) is excluded by the census cliff
(7K-b: tail noise ≈ formal) and by Lindegren; the forced-multiplicity
cell misfits by 135–153 and its flat-q premise is rejected by the
certificate's own q-table at −162. The defensible span is ≈ 0.6–0.8
(noise-cap choices within the physical envelope + grid extension),
with the 0-configurations quoted as REJECTED-PREMISE conditionals.
Conceding "separately defensible" would give away ground the program
measured. → abstract form fixed in the plan with this span.

**ADOPTED (3) — the orphans**, with one staleness correction. The
pair-common brightness vector and the 0.2–50 kAU subsystem rate get
standalone-note flags (author decision added to the plan). The
"tension needs a named home" point is structurally right —
Paper 1 gains a named section, "The multiplicity interface" — but
his content is STALE and repeats his own round-10 misread: (a) "the
likelihood is indifferent at +0.0" — corrected in round 10 and
re-corrected here: the likelihood PAYS 135–153 to avoid the forced
cell; the +0.0 is α-indifference WITHIN the forced world only
(CORRECTION #B, second occurrence of the same misread); (b) "the
tension has no home and quietly becomes a caveat" — the ×3 tension
DISSOLVED at measurement grade in the certificate round (twin-heavy
q-table +162; conversion band [0.10, 0.39]; twin-convention host
≈ 0.23 ≡ flat-q kinematic 0.10–0.15). What the named section houses
is the INTERFACE result — a scalar companion fraction is the wrong
interface between photometric and kinematic measurements; the
q-resolved conversion is the right one; the residual open item is
the joint population fit — stated as a resolved-then-residual
structure, not an open tension (CORRECTION #C, staleness).

**ADOPTED (4) — Paper 2's replacement spine.** The "two systems, one
dial" device is retired outright (no caveated thread); the new spine
is his cleaner statement, elevated to a finding: SPARC measures c₁
interior with zero excluded across treatments; the binaries, once
nuisances are honest, do not resolve it — an instrument-resolution
result. His second point adopted prominently: the sixteen-law ±8
degeneracy is a METHODS finding about every binary function contest
in the field (including other groups'), stated as such in Paper 2.

**ADOPTED (5) — 7K-b promoted to Paper 1's flagship empirical
object**, with both cautions folded into standing rules: (i) the
Poisson P is never quoted without the counts (9, 2) adjacent;
(ii) "no configuration reproduces it" is a SEARCH statement — the
searched space is enumerated wherever the claim appears (the landed
grid + variants + smear-off + PHYS + boost cells; the same
discipline as the flatness-axis rule).

**ADOPTED (6) — the 7L reframe**: "modeling-philosophy fork" →
"well-posed empirical question: is the removed 92% data or
contamination?" — with the settling measurement NAMED: run the
completeness + width machinery ON the removed subsample (do the
ΔRV > 10 / wide / distant pairs' kinematics decompose into the
measured companion + noise sectors, or carry the boost?). Logged as
TODO (new instrument, not opened; the instrument freeze is on
functions, and the queue is complete — it waits its turn).

No credence move (interim round; cadence rule). PAPER → v3.8
(three claim-hygiene edits in the monolith so the assembly source
carries the round: the curve-form sentence in the 7J-z6 passage,
the counts+search-space sentence in the 7K-b passage, the
well-posed-question sentence in §7.4(d)). Split plan updated in
place. Reply to Opus: NOT drafted (standing rule — awaiting the
user's go-ahead).

## Review round 12 (2026-07-28): Opus's close — B-refinement adopted, decision-#4 recommendation logged, and THE TWIN-Q QUESTION answered in two layers + pre-registered as 7J-z7

He owns C (second occurrence); takes A as stated (his restatement is
clean: "an unmodeled object doesn't vanish at the physical envelope,
it relocates, and if its fingerprint sits in the radial column the
relocation target is α itself"); accepts B's substance.

**ADOPTED — the B refinement (load-path abstract):** the rejection
basis travels WITH the span in the abstract, not in §6: "α ≈ 0.6–0.8,
with lower configurations excluded by the band–cliff consistency
requirement and the twin-weighted q-table." His reason is the right
one: α has been a confident number four times; what distinguishes
this one is that its floor is held by two statistics EXTERNAL to the
likelihood that produces it, and the abstract should show the load
path so a revision of either statistic visibly moves the number.
→ split-plan abstract form updated.

**Decision #4 — his recommendation logged (the DECISION stays the
author's):** photometric-coherence vector = standalone short note
(self-contained stellar astrophysics, no gravity context needed);
subsystem rate STAYS in Paper 1 (its credibility chain — completeness
instrument + certificate — lives there) but titled and abstracted AS
A RESULT, not a method. Both sensible; awaiting the author's yes.

**THE SUBSTANTIVE ITEM — answered in two layers:**
- LAYER 1 (his premise corrected): the S3/forced-fcomp collapse was
  NOT produced by the retired as-published amplitude law. The GB0w
  first-firing catch (amendment 7e) — the exact bug class he cites —
  happened BEFORE any photow cube shipped: the invalid pair was
  deleted, and every operative photow/photow3 cube since carries the
  photocenter-cancellation law |q/(1+q) − ℓ/(1+ℓ)| (GB0w = 0.00e+00
  vs the photo cubes on every run since, including all of
  yesterday's). 7J-z3's S3 ran on those cubes with exact G0/G1
  regressions. "The collapse was an artifact of a law you've since
  retired" is therefore false as stated.
- LAYER 2 (his underlying demand is right): the fitter's companion
  q-DISTRIBUTION is still FLAT (0.1–1). Twins are wobble-quiet
  through the law (wfac(q=1) = 0 exactly), but the FRACTION of
  wobble-quiet companions at a forced fcomp is set by the q-draw —
  and the twin-heavy conversion currently enters the operative
  pipeline only at PRIOR level (the LANDED-CONV widening) and at
  TRUTH level (the arm-suite injections). A solver-level forced scan
  with the measured twin-t5 q-law has never run. His "dead in the
  code, not just on paper" demand is well-posed and cheap.

**7J-z7 pre-registration (committed BEFORE execution): THE TWIN-Q
FORCED SCAN.**
- IMPLEMENTATION: env QLAW='t5' in stage7j_marginal.py — build_pop's
  companion q redrawn from the arm-suite twin-t5 marginal law (the
  GV7 winner convention: split 0.9/1.4, low branch uniform [0.1,0.9],
  spike uniform [0.9,1.0]) via SINGLE-UNIFORM inverse CDF from the
  SAME rng slot that fed the flat draw — every other draw's stream is
  untouched, so the fcomp = 0 slice must be BIT-IDENTICAL to the
  operative photow3 cube.
- GATE G0-q: max|twin-cube(fcomp=0) − photow3(fcomp=0)| ≤ 1e-9
  (expect exact 0.0), ABORT on fail. TAG '_qt5' (operative cubes
  never overwritten).
- READ (stage7jz7_read.py, pre-framed): the S3-style forced read —
  marginal at LANDED-CONV restricted to fcomp ≥ 0.35 — on the twin
  cubes, both laws, seeds 31 + 101 (matching S3's 4/4 grammar); plus
  the FREE marginal on the twin cubes (descriptive: does twin-q move
  the unforced answer?) and the forced cell's misfit (descriptive:
  S3's forced cell misfit 135–153 under flat-q).
- BARS (locked now): FIFTH-MOVE-DEAD-IN-CODE = forced α̂_marg ≥ 0.5
  AND dN ≥ +10 in ≥ 3/4 law × seed reads (the collapse was a
  flat-q-convention artifact; the conversion arithmetic verified at
  solver level). FIFTH-MOVE-ALIVE = forced α̂ ≤ 0.3 in ≥ 3/4 →
  MATERIAL (the exposure reopens, routed to the next decider). Else
  PARTIAL (quoted).
- EXPECTATION (stated, not a bar): host 0.35 under twin-t5 ≈
  kinematic-equivalent 0.12–0.18 (the conversion band) = inside the
  anchor-flat region → α ≈ 0.65–0.75 → DEAD-IN-CODE predicted.
- No credence move this round regardless (cadence; interim round).
His exit standing: ping on structural events (referee report; the 7L
settling measurement). Reply-to-Opus with this round's results:
drafted only on the user's go-ahead, per the rule.

## Stage 7J-z7 EXECUTED (2026-07-28): FIFTH-MOVE-ALIVE (MATERIAL, 4/4) — the expectation MISSED, the mechanism decomposes, and round-11's defense is self-corrected

G0-q EXACT 0.00e+00 ×4 (the stream-preserving redraw verified — the
fcomp = 0 slices are bit-identical to the operative cubes). **The
pre-registered ALIVE bar fired 4/4: forced fcomp ≥ 0.35 under the
measured twin-t5 q-law gives α = 0.00, dN = +0.0 in every law × seed
read.** The stated expectation (DEAD at α ≈ 0.65–0.75 via the
conversion band) was WRONG.

**THE MECHANISM (why the conversion arithmetic missed):** the twin
world trades channels. Twins are wobble-QUIET (photocenter
cancellation — the conversion band converted exactly this, the
wobble moment) but hidden-mass-MAXIMAL (boost = √(1+q·M_h/M_s), the
round-10 "non-binding" clause — whose scope was the FREE optimum
only). Mass inflation is multiplicative and γ-preserving — the most
boost-degenerate absorber in the model. Under twin-q forcing the
wobble damage the data reject shrinks (**cost-to-force HALVES:
135–153 flat-q → 69.2–86.6 twin-q**) while the mass channel absorbs
α entirely. Opus's round-12 instinct is VINDICATED in consequence
(the solver-level scan was needed) though his premise was wrong
(the S3 collapse never ran on the retired law).

**SELF-CORRECTION to round 11 (correction B's multiplicity clause
OVERSTATED):** "the forced-multiplicity cell is killed by the
q-table at −162" is true only of the flat-q FACE premise; the
twin-q forced conditional is NOT premise-rejected — it is
disfavored (69–87 lnL + the certificated rate 0.23–0.32 sitting
below the 0.35 threshold) but LIVE. The load-path abstract revises
accordingly (see the split plan): the α floor's supports are the
band–cliff consistency (noise side) and the measured-rate-vs-0.35
margin + the forcing cost (multiplicity side); the q-table clause
demotes from "excludes" to "prices".

**DIVIDEND (positive):** the FREE twin-q marginal reads α =
0.61/0.78/0.84/0.63 at dN +13.9–18.1 — the operative band under the
solver-level q-law swap: **the unforced answer is q-law-ROBUST.**

No credence move (pre-stated; the MATERIAL routes forward, below).

**7J-z8 pre-registration (committed BEFORE execution): THE ADJACENT
STATISTICS vs THE TWIN-FORCED WORLD.** The forced world's absorber
is now hidden mass — which is ratio-FLAT in separation (it cannot
produce an s-RISING median) and inflates the near-parabolic
pericenter pile (it must flood the census overshoot). Both adjacent
instruments exist; both are read AT the twin-forced cell (the
argmax of the forced slice of the qt5 cubes, per law, seed 31; the
7K machinery with the t5 q-draw):
- MEDIAN leg (7K-a statistic): R_forced ≥ 1.052 → the forced world
  reproduces the anchor; R ≤ 1.030 → REJECTED by the median; GRAY
  between.
- CENSUS leg (7K-b statistic): P(≤ 2 observed overshoot | μ_forced)
  < 1e-3 → REJECTED by the cliff; ≥ 0.01 → consistent; GRAY
  between.
- VERDICT: EXPOSURE-CONTAINED if ≥ 1 leg REJECTS in both laws (the
  fifth-move world is excluded by data external to the 2D
  likelihood — the 7K-b self-defense structure applied to the
  exposure itself); EXPOSURE-STANDS if both legs are consistent in
  either law (flagged for the next named decider; credence
  untouched either way this round — cadence).
- Gate: the z8 forward at the forced cell with fcomp set to 0 must
  reproduce the 7K-a G0 legs' behavior (the machinery identity —
  scale-free control PASS carried over; the t5 redraw is the
  z7-verified stream-preserving one).

## Stage 7J-z8 EXECUTED (2026-07-28): EXPOSURE-CONTAINED — the cliff holds both flanks, and the median honestly retires from this defense

**Verdict per the locked grammar: EXPOSURE-CONTAINED** — the CENSUS
leg REJECTS the twin-forced world in BOTH laws at **P(≤ 2) =
5.5e-14** (μ_overshoot = 37.1 vs the observed 2 — the worst cliff
violation of any configuration ever tried: the hidden-mass boost
pushes the near-parabolic pericenter pile past 1.67 wholesale;
μ_band = 36.1 vs 9). The forced cell (α = 0, fcomp = 0.35, wr = 0.3,
fpm = 3.0, kw = 0.7, sq = 0.0 — identical both laws, α = 0 is
law-blind ✓) is excluded by a statistic EXTERNAL to the 2D
likelihood, at fourteen orders.

**THE HONEST TWIST (expectation-miss #2 of the round, logged):** the
MEDIAN leg came out the other way — **R_forced = 1.085: the
twin-forced Newton world REPRODUCES the 1.078 anchor** (my "mass is
ratio-flat so the median rejects" reasoning ignored that at
fcomp = 0.35 the remaining ~57% non-twin companions DOUBLE the
√s-rising wobble channel relative to the free Newton-best cell's
fcomp = 0.1, and the noise at fpm = 3.0 adds its own rise). The
median therefore RETIRES from the fifth-move defense: it cannot
discriminate the forced-companion world from the boost world. Its
7K-a role (half-absorbed excess over the FREE Newton-best world)
stands unchanged — the retirement is specific to the forced
conditional.

**THE FINAL FORM OF THE EXPOSURE (the round-12 arc closed):** the α
claim's multiplicity floor is held by (i) the forcing cost (69–87
lnL) plus the certificated rate (0.23–0.32 by convention) sitting
below the 0.35 collapse threshold, and decisively (ii) **the census
cliff: ANY fcomp ≥ 0.35 world — flat-q via wobble texture, twin-q
via mass inflation — floods the overshoot the data keep at 2.** The
same statistic already held the noise flank (7K-b). ONE
model-light, internally-self-defending statistic now guards BOTH
collapse corners — the load-path abstract clause becomes exact:
"α ≈ 0.6–0.8, with the collapse configurations (noise beyond the
physical envelope; forced host multiplicity ≥ 0.35) each excluded
by the band–cliff consistency requirement." Opus's round-12 closing
("the census pair keeping its counts attached") is now true for a
measured reason: the pair is the paper's floor-holder, full stop.

Plain verdict (7J-z7 + 7J-z8 together): SUCCESS — the reviewer's
code-level question was answered by running it (his premise wrong,
his instinct right); the fifth move is ALIVE inside the 2D
likelihood (MATERIAL, honestly booked, with the mechanism
decomposed: wobble→mass channel swap, cost halved) and CONTAINED
outside it (the cliff at 5.5e-14); both of my stated expectations
missed and both misses are logged — the bars, locked before each
run, did exactly the work they exist for. No credence move
(cadence; the containment is content for the next named decider).

ELI12: The reviewer asked "does your hidden-companion escape hatch
stay shut when the companions are the twin-heavy kind you actually
measured?" We ran it. Surprise one: the hatch OPENS — twins fake
the boost through hidden weight instead of wobble, at half the
old price. Surprise two: we then walked that twin-world outside to
face the two independent witnesses, and the second witness
demolished it — a sky full of hidden twin-weight would hurl
thirty-seven pairs past the speed cliff where the real sky has
two. The first witness (the median) turned out to be fooled by
this world — noted honestly, it steps down from guard duty. End
state: one incorruptible witness, the cliff, now guards both doors
the signal could have escaped through — and both of my own wrong
guesses along the way are in the log, because the rules were
written down before each answer came back.

## Mechanism-consistency note (2026-07-28, ten-minute grade): VERLINDE EXCLUDED AT LADDER LEVEL

Prompted by the "other interpretations" discussion: the emergent-
gravity interpolation is a closed form, so its ladder is exact
algebra against measurements we already own (freeze-compliant
consistency row; [calcs/note_verlinde_c1.py](calcs/note_verlinde_c1.py),
gates G1/G2/G3 sympy-exact). The additive apparent-DM reading
(Brouwer+17-tested): ν_V = √(π/3)/x + 1 exactly — **c₁ = 1, c₂ = 0.
The measured c₁ (0.21–0.52 across treatments; 0 excluded at 7.5σ
grade at distance 0.45) puts c₁ = 1 farther outside every interval
than the already-excluded 0 — while the AMPLITUDE is viable: a_M =
cH₀/6 vs a₀ = cH₀/2π is a 4.7% shift, inside the a₀ error. The
1/6-vs-1/2π near-coincidence is exactly why amplitude-level tests
(lensing normalizations) could never see this; the digits see it.**
Caveats carried (spherical/quasi-static derivation; disk corrections
not computed; λ-grid ended at 0.625 → the Δ at c₁ = 1 is
curvature-extrapolation grade); scout on prior coefficient-level
Verlinde exclusions PENDING — no novelty claim printed. Plain
verdict: SUCCESS (one of the six rival readings converted from
philosophy to an excluded row, in ten minutes, on existing
measurements). ELI12: the "stretchy vacuum" theory writes its own
serial number — a 1 where our universe measures roughly one-half —
and its disguise was that its VOLUME knob (1/6) accidentally
matches ours (1/2π) to five percent, so everyone checking volume
saw nothing wrong; we checked the digits.

## Literature note (2026-07-27): El-Badry & Rix 2018 PRIMARY READ — the flagged ~20%/component claim VERIFIED

Flagged at scout level during the certificate round; primary read of
arXiv:1807.06011 (MNRAS 480, 4884), §2.4.1 "Higher-order multiples":
- MEASURED: "≈ 10.5% of primaries and secondaries with 1<(GBP−GRP)<2
  fall above this line and likely have an unresolved companion with
  q ≳ 0.5" (CMD cut MG = 2.8(GBP−GRP)+2.4 — the regime where the
  binary sequence separates cleanly = the TWIN-ADJACENT regime).
- EXTRAPOLATED: "∼ 20% of main-sequence primaries and secondaries have
  an unresolved main-sequence companion, such that at least one
  component has an unresolved main-sequence companion in ∼ 36% of the
  wide binaries" + a few percent unresolved WD companions.
CONSEQUENCE for Paper 1: the catalog authors' own photometric
characterization brackets our v2c host rate (0.29–0.32 blended,
0.23 twin-convention ≈ their 20% + WD few %) — the "×3 vs literature"
tension was a comparison against short-period spectroscopic subsystem
rates (Tokovinin ~0.1/component), not against the all-separation
photometric census; the certificate round's resolution direction is
externally confirmed. Their measured core (10.5% at q ≳ 0.5) is
exactly the twin-heavy regime the GV7 q-table found. Caveats to carry:
DR2 vs EDR3 catalog; their CMD method shares the photometric channel
with ours (common-mode caveat applies to both; our kernel measurement
7J-z quantified it); "unresolved" = all P below the ~1″ resolution
limit, wider than the wobble-relevant window. Consistency, not
circular confirmation — cite §2.4.1 in Paper 1's completeness
discussion. (3J's 12.3% overluminous ≈ their 10.5% likewise.)

## Stage 8A (2026-07-29): THE RIVALS' LADDER — the fingerprint table: two catalog kills, one identity, one class theorem, and the catalog-wide pincer

THE RIVALS ARC opened (user directive: "shoot down more theories and
strengthen or shoot down ours"). Standing arc rule per the same
directive: **PAPER.md stays FROZEN at v3.9 for the arc's duration** —
results book to NOTES/LEDGER as always, paper absorption happens once
at arc close (no per-stage version churn). 8A is freeze-compliant by
construction: consistency rows only, no new fit anywhere — the only
new content is exact series algebra (sympy-gated, all first-run PASS)
and three Haiku scout sweeps (sim-data / catalog+prior-art / analytic
ΛCDM), primary-source flags carried per the 4C discipline.

**Instrument ([calcs/stage8a_ladder.py](calcs/stage8a_ladder.py),
gates G1–G9 exact):** the deep ladder ν = A/x + c₁ + c₂x + … as a
CLASSIFIER, priced against the existing measurements (c₁ flat
0.385–0.519 [4S] / hier 0.208–0.309 [4Z], bootstrap ~0.4±0.3, c₁=0
excluded Δ(−2lnL)=56.3; binary c₁ ERASED [7J-d], not quoted).

1. **THE MASTER INVERSION (G4, exact):** for any AQUAL μ(z) = z +
   m₂z² + m₃z³, c₁ = −m₂/2 and c₂ = 5m₂²/8 − m₃/2 — **the measured
   c₁ is a measurement of μ''(0)**; the data demand m₂ ∈
   [−1.04, −0.42]. Four member cross-checks exact.
2. **THE STANDARD FAMILY DIES BY ITS ZERO-POINT (G5):** the whole
   n ≥ 2 family (incl. the field's historical "standard"
   μ = z/√(1+z²)) has **c₁ = 0 EXACTLY** — excluded at the c₁=0
   grade (56.3), with the direct contest already on file (4B:
   198–200/200 raw, −56 honest) — the ladder supplies the WHY.
3. **THE IDENTITY (G6, exact): exp-μ ≡ boot.** μ(z) = 1 − e^(−z) in
   AQUAL is EXACTLY the 4F quantum-bootstrap bath ν = 1 + n_BE(νy)
   (both solve e^u = ν/(ν−1)); fingerprint (¼, 7/96) = the 4F
   printed cell; ν(1) = 1.3500. Scout: the bare exponential μ is
   apparently NOWHERE named or tested (attribution NOT FOUND) ⇒
   **the program's boot adjudication (4F raw-dead → 5C hier-flip
   +75.6 → 5F binary veto → 5M vertical collapse) is apparently this
   function's FIRST data contest (scout-level), and it is dead on
   the primary treatment.** Companion identity re-gated (simple-μ ⇔
   classical bath). Reading: the μ-side's two natural members ARE
   the two self-consistent bath closures; the data kill both and
   keep the source-driven occupation forms — the same exponential in
   the temperature variable x (McGaugh RAR fit = BE, surviving) vs
   the acceleration variable νy (dead): **the data distinguish the
   exponential's ARGUMENT.**
4. **THE ADDITIVE CLASS THEOREM (G7, symbolic-F exact):** any
   apparent-DM model g_obs = g_N + √(a₀g_N)·F(g_N/a₀) with F
   analytic, F(0)=1, has **c₁ = 1 EXACTLY** (F-independent, odd
   rungs only). Members: Verlinde (exact — note-V), the
   superfluid-DM MOND-limit composition (membership FLAGGED —
   primary equation not retrieved by scout), any dark component
   tracking baryons analytically. c₁ ≠ 1 requires half-integer
   powers of g_N — the temperature variable. Note-V's pending scout
   RESOLVED: prior art = Lelli+17 (arXiv:1702.04355) fit-level EG
   exclusion (M/L amplitude + radius-residual); the
   coefficient-level exclusion apparently new (scout-level).
5. **THE HEES+16 CATALOG READ (G9, arXiv:1510.01369, forms
   scout-quoted, ladders exact):** ν̃_a: **c₁ = a — the family
   parameter IS the zero-point coefficient** (viable iff a ∈
   [0.21, 0.52]); ν̄_a: c₁ = 1 − 1/(2a) (viable a ∈ [0.63, 1.04]);
   ν̂₁ ≡ BE ≡ the RAR fit (gated identity — McGaugh's function is
   the a=1 hat member); ν̂_{a≥2}: c₁ = 0 (the standard fingerprint).
6. **THE PINCER (G10):** their Cassini verdicts (ν̃ dead at ALL a;
   ν_α/ν̂ need a ≥ 7–8 = c₁ = 0 members; ν̄ needs a ≥ 2 = c₁ ≥ ¾)
   × the measured window = **EMPTY INTERSECTION across the entire
   published catalog** — the amplitude-locked solar tension
   (4K/5S/5I) rederived in coefficient language through the field's
   own function zoo; escape doors structural as before (7G
   trajectory formulation, 451 orders; 6W scalar-EFE excluded).
   Caveat: their bounds scout-quoted, MG/AQUAL-conditional,
   fixed-a₀ — primary table read queued before any paper use.

**Scout residue (all scout-level):** (i) coefficient-level
expansion-and-constraint of the interpolating function NOT FOUND —
the ladder instrument apparently unpublished; (ii) NO closed-form
ΛCDM ν(y) exists (Navarro+17 sim-based; Paranjape–Sheth equations
inaccessible; Desmond semi-empirical, scatter overpredicted 3.5σ) ⇒
**8B must be data-side** — and the sim-data survey found nobody
publishes (g_bar, g_obs) points: FIRE-2 snapshots open (flathub, no
login; Mercado+23 = 20 galaxies with "hooks & bends"), TNG needs
registration, EAGLE/MUGS2/MassiveBlack-II unclear — pre-computed
-profile feasibility scout in flight, route decision after; (iii)
**Magneticum (Mayer+22, arXiv:2206.04333): ΛCDM-simulated a₀ RISES
≈×3 by z = 2.3 ≈ the H(z) trajectory (×2.97 at z=2)** — the a₀(z)
leg ALONE is not ΛCDM-immune; P3's discriminating content is the
LOCK (paired p_gal(z) drift), as registered — P3 annotated, no
status change; observational a₀(z) fits also exist (a₁ ≈
1.6e−10/z scout-reported) = an 8E lead, primary reads owed.

Plain verdict: SUCCESS — one afternoon of exact algebra plus the
measurements we already own prices the field's entire function
catalog: the standard family and the additive apparent-DM class are
excluded on the zero-point digit, an untested function is identified
as our already-dead boot cell, and the Saturn tension is shown
catalog-wide in coefficient language. The arc's shoot-down-OURS
instrument (8B, the sim ladder) is scoped: no algebra shortcut
exists, the data route is open at FIRE-2, decision pending the
feasibility scout.

ELI12: Every gravity formula carries a serial number — the first
digits of how it switches on. We measured digit #1 ≈ ½ months ago.
Today we machine-read the serial number of every formula in the
field's catalog plus the big rival theories: the classic "standard"
formula prints 0, and the whole invisible-matter-shadowing-the-stars
family (Verlinde's, the superfluid one) prints exactly 1 — no
tuning can change either, it's baked into their structure. Neither
matches ½. Bonus one: a formula nobody ever bothered to test turns
out to be secretly the SAME equation as the "self-heating bath" we
invented ourselves and already killed with data — so it arrives
pre-dead. Bonus two: crossing our digit with the field's own
Saturn-probe test, NO catalog formula passes both — same message
our Saturn work already carried: the fix must be structural, not a
different formula. Next up: check whether plain cold-dark-matter
simulations accidentally print a ½ — that's the test that could
shoot OUR reading down, which is exactly why we're running it.

## Stage 8B pre-registration (2026-07-29, committed BEFORE execution): THE SIMULATION LADDER — does the ΛCDM attractor print the zero-point digit?

THE SHOOT-DOWN-OURS INSTRUMENT. The measured c₁ ≈ ½-grade zero-point
is the program's sharpest galaxy digit. If ΛCDM galaxy formation
(halo response + baryon physics, no horizon input) produces a
simulated RAR carrying the same digit, the digit's evidential weight
for the thermal reading DEFLATES to galaxy-formation phenomenology;
if the sim RAR carries a different fingerprint, the digit becomes a
measured DISCRIMINANT against the ΛCDM attractor. Both directions
are live; the bars are locked now, before any fit.

**DATA (access granted by the EAGLE team 2026-07-29, John Helly;
credentials in private/, never committed):** the public EAGLE
database, RefL0100N1504 (the 100 Mpc flagship), SnapNum 28 (z≈0),
Aperture table (enclosed Mass_Star/Gas/DM/BH at r = 1,3,5,10,20,30,
40,50,70,100 pkpc — probe-verified) joined to Subhalo. SELECTION:
SubGroupNumber = 0 (centrals; probe showed satellites carry
tidally-truncated flat profiles), Spurious = 0, M*(30 pkpc) ∈
[1e9, 1e11] M☉ (SPARC-like window) — **N = 7,239 galaxies
(probe-counted)**. Radii used: r ∈ {5,10,20,30,40,50,70,100} pkpc
(1 and 3 pkpc dropped: softening ~0.7 pkpc). Points: (g_bar, g_obs)
= (G·(M_star+M_gas)(<r)/r², G·M_tot(<r)/r²) — the spherical-aperture
proxy, stated as such. Cross-check box: RecalL0025N0752 (the
high-resolution recalibrated 25 Mpc run), same selection. Server
politeness: ~6 chunked queries (mass-bin chunks + MorphoKinem +
Recal), seconds each, cached to data/eagle/ (gitignored), one-time.

**ESTIMATOR (the 4S machinery, faithfully adapted):** family ν_λ =
(1−λ)·standard + λ·BE, c₁ = λ/2 exactly; λ grid −0.30…1.50 step
0.05; per λ profile (log a₀, s_int) — no M/L nuisance (sim masses
exact; all scatter intrinsic, lognormal, s_int profiled); Δ=1
profile interval + 200-rep GALAXY bootstrap (resample the 7,239,
λ free per rep) = the primary error. PRIMARY fit: full sample,
stars+gas baryons, all 8 radii.

**GATES (ABORT on fail, no verdict quoted):**
- G0 units: hand-check V_c(10 pkpc) for one probed galaxy against
  the constant path (4.302e-6·(M/M☉)/(r/kpc) km²/s²); sane range.
- G1 the λ-family series regression (c₁ = λ/2 at λ = 0, ½, 1 —
  the verbatim 4S block).
- G2 descriptive sanity: the sim RAR must be TIGHT and MONOTONE
  (binned medians monotone in g_bar; raw scatter < 0.3 dex) — else
  the aperture proxy failed as an instrument.
- G3 λ̂ interior on the grid (edge-riding = no measurement).
- G4 injection: synthetic skies at the same point distribution with
  known λ_true = 0.0 and 0.75 (+ lognormal s = 0.1) must recover
  λ̂ within ±0.05 — estimator-bias control.

**BARS (locked now; on the PRIMARY fit's bootstrap 16–84 c₁
interval and median):**
- B1 REPRODUCES-THE-DIGIT: median c₁_sim ∈ [0.15, 0.60] AND
  P(c₁ > 0.05) ≥ 0.975 AND P(c₁ < 0.95) ≥ 0.975 → the ΛCDM
  attractor manufactures a half-grade zero-point; the sky digit's
  horizon-uniqueness DEFLATES (named honest strike direction;
  credence move deferred to the next named decider per cadence —
  this stage books the fact).
- B2 DISCRIMINATES: median c₁_sim outside [0.15, 0.60] AND the
  bootstrap interval excludes [0.21, 0.52] entirely → the measured
  digit separates the sky from the ΛCDM attractor (the rival-control
  the coefficient measurement has lacked).
- B3 AMBIG/UNSTABLE: anything else — OR any pre-named robustness
  leg moving median c₁ by > 0.25: (a) stars-only g_bar (hot-gas
  bracket), (b) r ≥ 10 pkpc (softening/inner-structure bracket),
  (c) disk subset KappaCoRot > 0.4 (MorphoKinem), (d) M* > 10^9.5
  (resolution bracket), (e) the Recal 25 Mpc box (volume/resolution
  cross-check). Then: instrument-limited at aperture grade, quoted;
  successor = the FIRE-2 snapshot pipeline (commissioned-only).
**Secondary (report-grade, no bars):** a₀_sim both boxes vs cH₀/2π
and vs Ludlow+17's published ballpark (~2.6e-10, scout-grade);
scatter-vs-x shape (does the sim produce the 4T deep-growing trend
and the x≈1 bump?); mass-tercile a₀ trend (the literature's
mass-dependence claim). PRE-STATED CONFOUND: the aperture c₁ folds
in sphericalization + hot-gas systematics — B1 is claimable only if
it survives the (a)–(e) stability clause; the confound direction is
NOT assumed benign.

No credence move this stage regardless (cadence: not a named
decider; a B1 outcome is flagged INTO the next named decider).

## Stage 8B EXECUTED (2026-07-29): THE SIMULATION LADDER — G3 GATE-FAIL / OFF-FAMILY: the ΛCDM attractor does not print the digit; it prints the ADDITIVE class, exactly as 8A's theorem predicts

**The amendment chain (all logged pre-quote; no run's numbers were
quoted as results before its gates were resolved):** AMENDMENT 1 —
the first run rode the inherited 4S grid's λ = 1.50 edge (G3-serving
grid extension to 3.0 + a family-validity guard; the guard's first
implementation truncated from the wrong end — negative-λ members are
tail-non-monotone — fixed to keep the contiguous valid run around
BE, λ ∈ [−0.10, 2.55]). AMENDMENT 2 — the extended run rode the
FAMILY-VALIDITY boundary itself (λ̂ → 2.55, still climbing):
restructured to a formal G3 GATE-FAIL report (no c₁_sim, B-bars
UNFIRED by the letter, bootstrap skipped — a boundary pile-up is not
an interval) + two LABELED DIAGNOSTICS. Both edge runs stand as
gate-fail records only.

**Gates:** G0 units PASS; G1/G1b series + validity PASS; G2 PASS —
the EAGLE aperture RAR is TIGHT (0.132 dex raw about the running
median) and monotone; **G4 injections PASS sharply (λ_true 0.00 →
−0.008; 0.75 → 0.784): the estimator recovers in-family truths —
the edge is the DATA, not the tool.** Support check: every aperture
point lies inside the sky measurement's joint support (deepest
2.8e-14 > the 4E lensing leg's 5.6e-15) — the off-family behavior
is physical, not a window artifact.

**THE RESULT:**
- PRIMARY (57,912 points): λ̂ = boundary-edge; along the slice the
  sky's digit region is rejected by Δ(−2lnL) c₁=½: +7,867 and
  c₁=0: +15,977 (NOMINAL point grade; 8 radially-correlated points
  per galaxy ⇒ effective deflation up to ~×8 ⇒ still ~+1,000-grade —
  decisive at any honest counting).
- Legs: (b) r ≥ 10, (c) disks-only (2,796 gal), (d) M* > 10^9.5,
  (e) the Recal 25 Mpc high-res box — ALL boundary-edge (the
  off-family pull is selection- and resolution-robust). TWO legs
  localize INTERIOR: **(a) stars-only g_bar → c₁ = 1.108 at a₀ =
  2.47e-10 (Ludlow+17's published ~2.6e-10 reproduced — external
  wiring consistency at scout grade); D2's kinematic window
  (g_bar ≥ 1e-12) → c₁ = 1.154 at a₀ = 1.85e-10.** Where the sim
  RAR can be localized on the slice at all, it prints ≈ 1.1 — the
  ADDITIVE digit, not the sky's ½.
- **D1, the fingerprint contest (the 8A classes head-to-head on the
  sim):** additive −6,004.8 vs BE (Ref100) / −217.7 (Recal25);
  simple −262.8; standard +8,110. **The ΛCDM attractor's RAR is
  best-described by the ADDITIVE class — the G7 theorem (a dark
  component tracking baryons ⇒ c₁ = 1 grammar) observed in silico.**
- D2: the sim's deep slope RUNS (0.51 → 0.58 → 0.62 → 0.61 per
  decade) — no fixed-a₀ MOND-form ν describes it globally (the
  Mercado "hooks and bends" statement reproduced at aperture grade).
- a₀_sim is UNLOCKED: 1.18–2.47e-10 across legs/boxes (2× spread)
  vs the sky's hier-locked 1.04–1.13e-10 (5M) — the attractor
  carries no temperature lock.

**CONSEQUENCE (the arc's question answered at this grade):** the
sky's zero-point digit c₁ ≈ ½ does NOT lose uniqueness to the ΛCDM
attractor — the attractor fails one level EARLIER (the functional
class): its RAR is additive-class, the sky's is occupation-class,
and 8A's algebra is exactly the classifier that separates them. The
sky's measured c₁ ∈ [0.21, 0.52] with c₁ = 0 and c₁ = 1 excluded is
now flanked by BOTH rivals on the c₁ = 1 side (Verlinde/superfluid
by theorem; the ΛCDM attractor by measurement) and the historical
catalog on the c₁ = 0 side. Formal verdict stays GATE-FAIL by the
pre-registered letter (B-bars unfired — no bootstrap interval
exists for an off-family optimum); the B3 successor stands
(FIRE-2 snapshot pipeline, commissioned-only, for beyond-aperture
grade). Caveats carried: spherical-aperture proxy (not disk-plane
kinematics — the class statement is aperture-grade); nominal
point Δs (×8 correlation deflation stated); EAGLE(+Recal)-specific;
sim scatter s_int ≈ 0.10–0.13 dex reported as-is (Desmond-direction,
aperture-inflated, not over-quoted). Server load for the whole
stage: ~15 queries, ~20 s total, one-time (cached).

No credence move (cadence: not a named decider; the result is
flagged INTO the next named decider alongside 8A).

Plain verdict: SUCCESS — the shoot-down-OURS instrument fired and
our digit survived it: the honest formal outcome is a gate-fail
(no c₁_sim exists because the sim never enters the measured family's
neighborhood), and the labeled diagnostics turn that into the
sharpest form of the answer — the attractor prints the OTHER
class's fingerprint, the one 8A proved belongs to
dark-stuff-tracking-baryons.

ELI12: We asked: if boring invisible-matter physics secretly runs
galaxies, would a universe simulated with exactly that physics
accidentally print our serial-number digit ½? We fed 7,239
simulated galaxies through the same measuring machine as the real
sky. Answer: the simulated universe doesn't print ½ — it can't
even be read on our dial without pegging the needle, and when we
check WHICH family of formulas it does match, it's the
"invisible-stuff-shadowing-the-stars" family with digit 1 — the
very family our theorem said dark matter must produce, and the
very digit the real sky rules out. So the real universe's ½ is
not an accident of galaxy formation — at least not in this
simulation, measured this way (the fine print: our quick method
reads each galaxy as a sphere, and it's one simulation code; a
heavier-duty check is queued if we ever want it). Our test alarms
all worked: two early runs hit the fence and were logged as fence
reports, and fake data with known answers came back exact.

## Literature notes (2026-07-30): the two 8A primary reads — Hees+16 VERIFIED (pincer premise holds + a NEW a₀ arm), BFK superfluid composition VERIFIED (additive member confirmed)

(1) **Hees+16 (arXiv:1510.01369, full-text read):** the Eq. (5a–5d)
family formulas are EXACTLY the 8A G9 forms (the scout was
verbatim-accurate), and the verdict sentence is verbatim: "the class
of transition functions ν̃α seems to be completely excluded by this
combined analysis. The functions να and ν̂α are excluded for low
values of α but begins to be marginally acceptable for large values
of α. The only class of functions that seem to be able to produce a
satisfactory fit to the galactic rotation curves without producing a
too large deviation in the Solar System is ν̄α for α ≥ 2." Their
Cassini bound: 0 ≤ Q₂ ≤ 6e-27 s⁻² (1σ, Hees+14); their external
field ge = 1.9/2.4e-10 m/s² (the AQUAL-side convention — the
3T-story echo, internally consistent). **THE PINCER UPGRADES
scout-quoted → PRIMARY-VERIFIED, and gains a SECOND ARM (their
Table 2, col 2): every Cassini-surviving member also fits rotation
curves only at an off-horizon a₀ — ν̄_{2…7}: a₀ = 0.72–0.82e-10
(25–35% LOW); ν8/ν̂6: 1.44–1.46e-10 (~35% HIGH) — while the sky's
occupation-class fits lock onto 1.04–1.13e-10 ≈ cH₀/2π (5M). The
catalog pincer is now two-armed: the digit AND the temperature.**
Subtlety caught (the 7F window-vs-asymptote lesson, logged): ν̄_{α=½}
≡ BE exactly (at α = ½ the first-term correction becomes the
constant ½; the asymptotic c₁ = 1 − 1/(2α) formula holds for α > ½,
so the viable band [0.63, 1.04] stands, but near the α → ½ edge the
window-effective c₁ exceeds the asymptotic value) — the pincer is
unaffected (its operative members α ≥ 2 are asymptotically clean);
the 8A row is annotated. Also verified: ν₁ = "simple", ν₂ =
"standard" (their words), ν̄_{0.5} "extensively used in Famaey &
McGaugh (2012)" = the RAR-fit/BE function appearing in the catalog
twice (≡ ν̂₁).

(2) **Berezhiani–Famaey–Khoury 2018 (arXiv:1711.05748, full-text
read):** Eq. (1) verbatim: **a = a_b + a_DM + a_phonon** — the
composition is ADDITIVE; and (their Sec. VII, verbatim): "we
nevertheless found that a_phonon ≃ √(a₀ab) to within a couple of
percent. In other words, the phonon force closely matches the
deep-MOND acceleration." ⇒ the phonon term alone is an
additive-class member with F = 1 + O(2%) ⇒ **c₁ = 1 by the G7
theorem — MEMBERSHIP CONFIRMED at primary grade** (the 8A flag
resolves), and the superfluid-core a_DM term adds FURTHER dark
response on top — the direction 8B measured for the attractor
class. Caveat kept: the phonon force operates only inside the
superfluid region (their EFE structure differs from MOND's); the
RAR-domain composition is Eq. (1).

Both PDFs + extracted texts cached in data/ (gitignored, chae
precedent). Plain verdict: SUCCESS — both pre-paper caveats closed
the day after they were opened, and the catalog pincer came back
stronger (two-armed). ELI12: instead of trusting summaries, we read
the two source papers ourselves. The rival-catalog paper says
exactly what we quoted — and its own table shows that every formula
surviving its Saturn test also needs a volume knob 25–35% off the
value the real sky picks; ours sits on it. The superfluid paper's
own first equation stacks its new force on top of ordinary gravity
— the digit-1 family, just as our theorem assumed. Two checkmarks
earned, one footnote fixed, nothing against us.

## Stage 8C pre-registration (2026-07-31, committed BEFORE execution): THE CEILING TEST — p ≤ ¾, parameter-free (P1's in-sample clause)

The gate's sharpest in-sample falsifier (6Y → PREDICTIONS.md P1):
the ambient-gated bath forces the screening-tail exponent p = ½ + g/4
with the gate g = [n_amb/(1+n_amb)]² ∈ [0, 1) ⇒ **p < ¾ EXACT for
every galaxy, parameter-free**; population median postdicted 0.689
(measured 0.65–0.75); void asymptote ~0.72. One galaxy measurably
beyond ¾ kills the gate. Correction #14 (6Z) withdrew the in-sample
ORDERING claim — this stage tests ONLY the ceiling (exceedance); the
e_N-split arm carries ceiling-only language by construction.

**INSTRUMENT (calcs/stage8c_ceiling.py):** the 5M vertical-hardened
hier machinery — the primary treatment: global (a₀, f_ML, s_int,
δ_lens) + per-galaxy dml (0.1-dex prior) + per-galaxy dv (measured
σ_v priors) + the Mistele lensing leg — with the tail dial
ν_p(y) = (1−e^(−y^p))^(−1/(2p)), the established 5G convention
(≡ BE at p = ½; the scale on which 0.65–0.75 was measured and 0.689
postdicted; convention-tied, stated).
- LAYER I (pooled + arms — the power-carrying layer): shared-p
  profile over p ∈ 0.40…1.10 step 0.05 (parabolic refine), Δ=1
  interval + 25-rep galaxy bootstrap (coarse-grid {0.5…0.9} +
  ±0.05/0.025 refine per rep, fit-lite rounds; design fixed now);
  ARMS: the Chae-matched subset split at MEDIAN log e_N(maxclust)
  (≈ 47/47; unmatched galaxies = a third reported arm; arms exclude
  the lensing leg — population-level, disclosed; the pooled fit
  keeps it, carrying the pooled p̂).
- LAYER II (per-galaxy census — the literal P1 clause): joint fit
  with per-galaxy p_g (scalar-profiled per round, bounds
  [0.30, 1.50]) alongside dml/dv; then per-galaxy p-scans (grid
  0.30–1.50 step 0.05) with LOCAL re-profiling of (dml_g, dv_g) per
  node → per-galaxy Δ(−2lnL) curves. QUOTABLE galaxy: Δ=1
  half-width ≤ 0.30 AND both interval edges interior. EXCEEDANCE:
  p̂_g > 0.75 with Δ(−2lnL at p = 0.75) ≥ 9 (≈3σ one-sided).

**GATES:** G1 identity + regression — ν_p(½) ≡ BE (grid assert
≤ 1e-10) and the all-p-pinned-½ converged fit reproduces 5M's dv-ON
BE value −12152.49 within 1.5; G2 nesting — every freed fit ≤ its
pinned comparator; G3 convergence traces; **G-INJ calibration on
three synthetic skies at the fitted globals/nuisances (dml, dv
drawn from their priors; noise from sig ⊕ s_int): two at p_true =
0.72 uniform (near-ceiling null) → the criterion may fire in ≤ 1
galaxy per sky, else the Δ bar is RAISED before the real census is
read (retune logged, real census untouched until then); one at
p_true = 0.90 uniform → POWER: the criterion must fire in ≥ 1/3 of
quotable galaxies, else LAYER II is declared POWER-LIMITED and
carries no bar** (the literal one-galaxy clause then waits for
better data — stated, never silently passed).

**BARS (locked now):**
- CEILING-BROKEN: any LAYER-I arm with p̂ − 0.75 > 2σ_boot; OR
  (LAYER II, only if powered) ≥ 2 calibrated exceedances, or one at
  Δ ≥ 16. → P1 flips to STRUCK in PREDICTIONS.md (in-sample kill);
  credence movement deferred to the next named decider, the strike
  booked immediately.
- CEILING-HOLDS: every arm ≤ 0.75 + 1σ_boot AND zero calibrated
  exceedances (if powered) AND ≥ 30 quotable census galaxies. P1
  stays LIVE-and-passing at SPARC grade (not "confirmed" — the void
  asymptote stays out-of-sample).
- Else AMBIG (quoted; the expected route if the census is
  power-limited and the arms are consistent).
No credence move regardless (cadence). Convention caveat carried:
p̂ is the ν_p-family dial (fractional-rung window effects near
p → ½, the 8A/ν̄ lesson); the AMB-structured free-gate variant is
the named successor refinement if AMBIG.

## Literature note (2026-07-31, 8E primary reads): THE a₀(z) LANDSCAPE — the P3 data era has opened, contradictory at systematics grade; the kill direction is NOT triggered

Primary reads (PDFs cached, gitignored): MUSE-DARK III (Ciocan,
Bouché, Fensch, Krajnović, Freundlich, Desmond, Famaey, Techi —
A&A 709, L16, arXiv:2604.22613) and Limbach, Psaltis & Özel 2009
(arXiv:0809.2790). Scout survey backing: MIGHTEE-HI (Vărășteanu+25,
arXiv:2504.20857), Mancera-Piña 2017 (arXiv:1703.06110), Shachar+23
BTFR (flat to z≈2.5), Li+18 (no galaxy-to-galaxy a₀ variation
locally).

**(1) MUSE-DARK III (verbatim numbers, 79 SFGs at 0.33 < z < 1.44,
pressure-support-corrected 3D forward modelling, MNR fits):**
a₀|z∼1 = 2.38 (+0.12/−0.10) e−10 (95% CI; ~19σ above canonical
SPARC 1.2); binned rise 1.99 → 2.71e−10 across their four z-bins;
linear fit **a₀(0) = 1.00 ± 0.04 e−10 (95%), a₁ = 1.59
(+0.10/−0.10) e−10 per z**; their MOND-framework refit (their
App. E) agrees within errors; scatter grows 0.13 → 0.19 dex with z.
Their own caveats (verbatim-grade): the linear parametrisation is
"a simple, phenomenological description... If the real z-dependence
of a₀ is non-linear or more complex, the resulting evolutionary
trend could be affected by systematic biases"; reconciling to the
z = 0 canonical value would need stellar masses +0.2–0.45 dex.
**READ AGAINST P3 (a₀(z) = cH(z)/2π):** the SIGN matches (a₀
rises); the INTERCEPT sits essentially ON the horizon value
(1.00 vs cH₀/2π = 1.08 — nearer the horizon than the canonical
1.2); the SLOPE at face value is ~2× the lock's secant (H-ratio
z 0→1 is ×1.76 ⇒ +0.82e−10 vs their +1.59e−10) — steeper than
the lock, with their systematic budget spanning the difference.
Notably Desmond AND Famaey (the DHF/Cassini authors) are on the
byline.
**(2) Limbach+09:** the scout-level "a₀ ∝ cH(z) excluded" claim
DEFLATES on primary read — exclusion is "within the formal
uncertainties" only, and their own conclusion self-disavows it:
systematics (IMF, stellar populations) "affect the data
significantly... Definitive conclusions can only be drawn after
the systematic errors are more precisely quantified." Face-value
TF-intercept direction is consistent with mild a₀ rise. NOT a
standing exclusion.
**(3) The landscape is internally CONTRADICTORY:** MUSE-DARK's
steep rise vs Shachar+23's flat BTFR (z to 2.5) vs MIGHTEE-HI's
a₁ = 4.47 ± 1.88 (2.4σ, consistent with MUSE-DARK within 1.5σ) —
method- and systematics-separated; no measurement cleanly brackets
the lock yet.

**P3 STATUS (annotation, no flip):** the kill condition ("measured
a₀(z) flat at z ≳ 1 at 3σ") is NOT triggered — the current
detections run the OTHER way; the a₀-leg's face-value slope
overshoots the lock ×2 at systematics grade; the pair's second
observable (p_gal(z)) is unmeasured everywhere — the LOCK remains
untested and is still the unique discriminant (Magneticum's ΛCDM
a₀(z) ≈ ×3 by z = 2.3 also rises — the a₀ leg alone separates
nothing). Successor instrument (8E, unopened): refit the
MUSE-DARK binned a₀(z) against cH(z)/2π with their M/L systematic
marginalized — a two-parameter contest (lock vs linear) their
published bins + errors already support.

Plain verdict: SUCCESS (reads done; P3's evidential environment
mapped; one scout overclaim deflated at source). ELI12: The first
telescope measurements of our "gravity's volume knob grows with
the universe's temperature" prediction are arriving. One big new
study says the knob DOES grow — even starting exactly at our
predicted today-value — but at face value it grows about twice as
fast as our formula says, while an older study says it doesn't
grow at all, and both admit their star-weighing could be off by
enough to flip their answers. So: nobody has killed the
prediction, nobody has confirmed it, and the tie-breaker we
registered (two dials moving together) hasn't been measured by
anyone yet. We wrote down exactly which follow-up computation
would settle our part of it.

### Addendum (2026-08-09, re-read with digit verification): one correction, four additions; the successor sharpens

A fresh scout surfaced Ciocan+ as "new"; the pre-action don't-redo
grep caught the existing 8E read only AFTER a duplicate fetch (the
duplicate PDF was deleted; process lesson logged in DIARY). The
re-read against the cached PDF produced one correction and four
additions to the note above; no status change, no credence move.

**CORRECTION (digit grade):** the clause "their MOND-framework
refit (their App. E) agrees within errors" is WRONG as written.
Their Eq. (4) DM-decomposition fit is a₀(0) = 1.0 ± 0.04, a₁ =
1.59 (+0.10/−0.10), "with the errors denoting the 95% CI"
(verbatim; so 1σ ≈ 0.05 on a₁). Their App. E MOND-framework
z-dependent fit (their Eq. E.5) is a₀(0) = 1.03 ± 0.05, a₁ = 1.20
(+0.10/−0.10) — 0.39 apart, NOT within errors; their App. D
ΛCDM-halo-profile variants give a₁ = 1.63 (+0.13/−0.12); their
individual-galaxy MOND regression gives a₁ = 1.42 (+0.94/−0.89)
(wide). THE MEASURED SLOPE IS FITTING-FRAMEWORK-CONDITIONAL at the
±25% level — direct evidence that the recovered a₁ is sensitive to
the assumed RAR form, which is exactly the axis the lock lives on.

**ADDITIONS:** (1) The lock compared correctly: linearizing
a₀(z) = cH(z)/2π over THEIR window (0.33 < z < 1.44) gives slope
1.04–1.09e−10 per z (H₀ 67.4–73). Against their DM-decomposition
a₁ = 1.59 (1σ ≈ 0.05): ~10σ face value. Against their OWN
MOND-framework a₁ = 1.20 (+0.10/−0.10): 1.6σ if the interval is
1σ, 3.1σ if 95% (convention for E.5 not stated in text; flagged).
The apples-to-apples row for this program is the MOND-framework
one, and it is NOT far from the lock. (2) Their Sect. 4 verbatim:
"our measured a₀(z) is faster than that of H(z) (Milgrom 1983a)" —
keyed to the 1.59 row; same paragraph notes Magneticum's ×3 rise
to z = 2 vs their inferred ×4. (3) Their Eq.-4 fit predicts a bTFR
zero-point shift ΔZP ≈ −0.2 dex by z ∼ 1; the bTFR literature they
cite is split (Übler+17 −0.44; Jeanneau+26 no evolution) — an
external consistency lever for the contest. (4) The data release
is per-galaxy: "All catalogues and data products from our
disk–halo decomposition, including the RCs" at
dark.univ-lyon1.fr/data-releases (their footnote 3) ⇒ the 8E
successor upgrades from published-bins grade to released-products
grade.

**THE SUCCESSOR, SHARPENED (still unopened):** the lock-vs-linear
contest should now be run IN-FORM with a lens leg — the locked
pair also drifts the tail exponent (p_gal 0.689 → 0.702 by z = 1),
their pipeline holds the RAR form FIXED at all z, and their own
1.59-vs-1.20 framework split shows form-coupling at the required
size. Leg A: contest lock vs linear on their released a₀(z) with
M/L marginalized (the 8E design). Leg B: generate locked-pair
truth, fit with their fixed-form + linear-a₀ pipeline, and measure
the recovered-a₁ bias — if the lens inflates 1.04 toward their
1.59, the tension dissolves and the lock EXPLAINS their excess
slope; if not, the tension is honest and stands. Pre-reg before
any run, power gates first (the 9R lesson).

Plain verdict: SUCCESS (one correction to our own lit-note, the
friendliest number in the paper recovered, the successor upgraded
to released-data grade). ELI12: We accidentally read the same
telescope paper twice — but the second read caught something the
first missed: the paper's "gravity knob grows twice as fast as
your formula" headline SHRINKS to "grows only a bit faster,
maybe compatible" when the growth is measured inside the same
kind of gravity formula we actually use. And how fast the knob
seems to grow depends on which formula they fit — which is
precisely the loophole our two-dials-locked prediction would
exploit. The follow-up computation is now fully specified and
their galaxy-by-galaxy data are downloadable.

## Stage 8C EXECUTED (2026-07-31): THE CEILING TEST — AMBIG by the locked grammar; the break side is clean, the per-galaxy clause is power-limited at SPARC grade

Runtime 2.5 min (the warm-started machinery; the exact G1b match
certifies convergence). **Gates: G1a identity ≤ 1e-10 PASS; G1b
pinned-½ reproduces 5M's dv-ON BE −12152.49 EXACTLY; G2 nesting OK;
false-exceedance calibration CLEAN (0/48 and 0/42 fired on the two
p_true = 0.72 null skies); POWER gate FAIL as a measured fact —
0/34 exceedances fired even at injected p_true = 0.90 ⇒ the census
carries no bar (pre-registered clause): the literal "one galaxy
beyond ¾" kill is UNTESTABLE at SPARC grade** (per-galaxy Δ(0.75)
≥ 9 needs constraining power individual SPARC galaxies do not
have).

**RESULTS:** LAYER I pooled p̂ = 0.647 (Δ1 0.628–0.667, INTERIOR;
p = 0.75 at +9.7 ≈ 3σ-profile above the minimum); galaxy bootstrap
16/50/84 = 0.534/0.617/0.799 (σ_boot = 0.133 — population variance
dominates, the 4S profile-vs-bootstrap precedent; P(p > 0.75) =
0.24). Arms (ceiling-only language BOTH directions per correction
#14 — the e_N axis is 6Z-unreadable in-sample, gate heterogeneity
+ type/mass confounds): eN-LOW 0.551 (0.495–0.614), eN-HIGH 0.904
(0.792–1.028), unmatched 0.628 (0.591–0.631). LAYER II census: 40
quotable galaxies, ZERO exceedances.

**VERDICT (locked bars): AMBIG** — CEILING-BROKEN fired nowhere
(pooled 3σ below; eN-HIGH's nominal excursion is 1.3σ arm-local,
inside the 2σ_boot bar; zero census exceedances); CEILING-HOLDS
not granted (the eN-HIGH arm exceeds 0.75 + 1σ_boot; census
power-limited). The pooled number doubles as a refined tail
measurement under the primary vertical-hardened treatment —
consistent with 5G's 0.65 and the 6E postdiction 0.689. P1 stays
LIVE, kill conditions unchanged (annotation added). Successor
(named, not opened): the AMB-structured free-gate variant (fit
the gate parameter itself, unclamped) + anchored/DR4-era data for
the per-galaxy clause. No credence move (cadence).

Plain verdict: NEEDS REFINEMENT — the instrument is honest and
calibrated, the ceiling took no hit, but SPARC cannot power the
sharp per-galaxy version of the test; the strong sub-results
(pooled 3σ below the ceiling with clean calibration) are banked.

ELI12: Our theory says every galaxy's "switch-off sharpness" dial
is capped at ¾ — find one clean galaxy past ¾ and the theory dies.
We built the test, calibrated it on fake universes (it never
cries wolf — and, honestly measured, it also can't hear a real
wolf at today's data quality: even a fake universe set beyond the
cap produced no detections, so single galaxies can't decide this
yet). The population as a whole reads ~0.62–0.65, comfortably
under the cap, right where the theory wants it. One subgroup
leans high but within the error bars we locked in advance — and
that particular sorting axis is one we already proved unreadable
with current data, so it gets noted, not spun. The cap survives;
the sharp version of the test waits for better telescope data.

## Stage 8D pre-registration (2026-07-31, committed BEFORE any cube): THE SETTLING INSTRUMENT — the completeness+width machinery ON the removed 92% (TODO-17; the round-11 named measurement)

The 7L fork, made operational: their cleaning discards 92% of the
joint sample (incl. the entire deep anchor); "is the removed 92%
data or contamination?" is the field's well-posed question and the
event Opus asked to be pinged on. 8D runs the landed machinery on
the COMPLEMENT of the 7L proxy mask (expected N ≈ 14071 − 1194 =
12,877 pairs) and asks whether its kinematics decompose into the
measured companion + noise sectors at externally-anchored levels,
or carry the boost.

**INSTRUMENT:** SAMPLE = 'anticook' in stage7j_marginal.py (the
complement mask; wiring committed WITH this pre-reg, before any
cube), operative photow3 config (photow amplitude law + FPME fpm →
3.0), seeds 31 + 101, both laws = 4 cubes. **MODE REGRESSION (the
GB0w rule): G17-ID — the partition identity.** The binned 2D
likelihood is pair-additive, so per cube cell lnL_full = lnL_cook
+ lnL_anti + (cell-independent const): gate = std over cells of
(full − cook − anti) ≤ 0.01 with max-dev from the median offset
≤ 0.05 (fp-summation tolerance; wiring errors are O(≫1)), per
(seed, law), on the existing full/cook photow3 cubes. ABORT on
fail.

**READ (stage8d_read.py, written after the batch; bars locked
HERE):** the anticook marginal at LANDED-CONV (the 7jz_read
construction) AND at the flat-fcomp prior (the anchor curve's two
ends, 7J-e3 discipline). Composition note pre-stated: the
complement is companion-ENRICHED by construction (the cook HR
filter removes overluminous pairs INTO it), so the landed anchor
is conservative there; both reads quoted, bars on the conjunction.

**BARS (locked now):**
- BOOST-CARRIED: α_marg(anti) ≥ 0.5 in ≥ 3/4 law × seed reads
  under BOTH anchors, with dN ≥ +10 in each → the removed 92%
  carries the boost at the landed absorber budget: the fork
  resolves toward DATA (their cleaning removed the signal's
  carriers). §7.4 gains the settling sentence; Opus-ping-worthy
  (after the user's go-ahead, per the standing rule).
- CONTAMINATION-DECOMPOSED: α_marg(anti) ≤ 0.2 in ≥ 3/4 reads
  under the LANDED anchor with the winning nuisance cell inside
  the external anchors on the COMPANION axis (fcomp ≤ 0.35; fpm
  EXCLUDED from the anchor clause — the width object is upstream
  and function-blind per 7J-d/7J-z6, disclosed) → MATERIAL against
  the α claim, routed to the next named decider.
- GRAY otherwise (quoted with the decomposition table).
**DIAGNOSTICS (no bars):** the price-of-contamination table (the
anti α = 0 slice's best cell vs the external anchors; its dN =
the cost-to-force); the census overlap statement (how many of the
(band = 9, cliff = 2) pairs the cook cuts remove — s-based if the
CSV lacks the distance/RV columns); the three-way partition
co-read (full 0.68–0.74 / cook ~0-powerless / anti = this stage).
**EXPECTATION (stated, not a bar):** BOOST-CARRIED — the cook
subsample carried ~3 lnL of the full +14.5–23.8, and the
complement holds the deep anchor and the census pairs.
No credence move this stage regardless (cadence; a
CONTAMINATION-DECOMPOSED outcome routes to a decider). PAPER stays
frozen at v3.9.

## Stage 8D AMENDMENT 1 (2026-07-31, logged PRE-QUOTE — the reader aborted at the gate; no marginal/verdict number was generated or seen): G17-ID's premise was WRONG for this likelihood; replaced by G17-ID2

The v1 gate (partition identity lnL_full = lnL_cook + lnL_anti +
const) FIRED FAIL on all 4 cubes (std ≈ 10–12 lnL, max-dev ≈ 50) —
and the diagnosis shows the PREMISE was mis-specified, not the
wiring: the residual is AXIS-UNIFORM (fpm slices: std 10.8–11.5
with no trend; sq slices 10.3–12.4; α slices 9.6–13.9) — the
signature of PER-SAMPLE TEMPLATE CONDITIONING. The v7/7J forward
model builds each s-bin's model distribution as a mixture over the
SAMPLE'S OWN pairs (masses, separations, errors enter the template
mix), so the full-sample model in a bin is the pair-weighted
mixture of the subsample models — Σ n ln(p_mix) ≠ Σ n ln(p_sub) +
Σ n ln(p_sub'): additivity is not a property of this likelihood BY
DESIGN. A partition identity therefore cannot certify a sample
mode here (it would fail for ANY correctly-wired complement).

**G17-ID2 (the valid certificate for a SAMPLE mode, pre-stated
before any read):** (a) the partition-COUNT identity — anticook
N = 12,877 = 14,071 − 1,194 exact (batch-printed); (b) the mask
certificate — data/stage7l_cookmask.npy dtype = bool (verified:
bool, sum 1,194, length 1,817,594 raw rows), so `ok & ~mask` is
the exact boolean complement (an int-typed mask would have failed
silently — checked); (c) the v1 residual-structure record kept as
the post-mortem (axis-uniform, seed-dependent const −3.5/−5.6 =
realization-mixed templates), explicitly NOT quoted as a pass;
(d) descriptive cross-read: the anticook PROF nuisance portrait
equals the landed one (wr = 0.2, fcomp = 0.1, fpm = 3.0, sq = 0.2
in 4/4) — no new-mode pathology. HONESTY: this certificate is
weaker than a GB0w-class bit-identity (none exists for
sample-conditioned models); the mode's entire new code is one
boolean line whose only observable effect (the count) is exact —
residual risk stated, not waved off. Verdict bars UNCHANGED.

## Stage 8D EXECUTED (2026-07-31): THE SETTLING INSTRUMENT — BOOST-CARRIED 4/4 under both anchors; the fork resolves toward DATA

Batch: 4 anticook photow3 cubes (12,877 pairs, both laws × seeds
31/101, 44 min GPU). G17-ID2 certificate PASS at its achievable
grade (amendment 1, logged pre-quote — see above; the v1
partition-identity premise was invalid for a per-sample
template-conditioned likelihood).

**THE VERDICT (locked bars): BOOST-CARRIED — the removed 92%
carries the boost at the landed absorber budget.**
- α_marg(anti) = 0.62 / 0.62 (simple, seeds 31/101) and 0.85 /
  0.81 (BE) at dN = +21.3…+25.9 — in EVERY law × seed read, and
  **anchor-INDEPENDENT: the LANDED-CONV and FLAT reads agree to
  0.2 lnL** (the likelihood alone drives it; the companion prior
  does no work).
- **The companion posterior is PINNED at fcomp = 0.10 in 8/8
  reads** — the removed sample does not want extra companions;
  contamination-reads 0/4. The pre-stated composition worry
  (complement companion-ENRICHED by the HR filter) turned out
  immaterial: even enriched, the kinematics refuse fractions
  above the measured rate.
- **Price-of-contamination: +21.5…+25.9 lnL** — and the α = 0
  world's own best cell sits at the SAME landed nuisance portrait
  (wr 0.2, fcomp 0.1, fpm 3.0, sq 0.2): the no-boost reading
  cannot even buy companions profitably; it just pays.
- Census overlap: ≥ 4/9 of the operative band pairs are removed
  by their s/d cuts ALONE (the ΔRV cut removes more — lower
  bound): **the (band = 9, cliff = 2) statistic lives in the
  complement.**
- The three-way partition closes the 7L arithmetic: full
  0.68–0.74 @ +14.5–23.8 | cook ~0-powerless (~3 lnL) | anti
  0.62–0.85 @ +21–26. "No signal in the clean 1,194" and "the
  signal lives in the removed 12,877 at the measured companion
  rate" are now BOTH measured inside one pipeline.

CONSEQUENCE: the round-11 settling measurement is EXECUTED — the
well-posed question "is the removed 92% data or contamination?"
now has a measured answer: at the landed absorber budget, with
companions at their externally-anchored rate, the removed sample
carries the boost. Their cleaning removed the signal's carriers.
This is the named Opus-ping event (his exit request); the note is
NOT drafted — awaiting the user's go-ahead per the standing rule.
No credence move (cadence — 8A/8B/8D all flag into the next named
decider). PAPER stays frozen at v3.9; §7.4's settling sentence
waits for the arc close.

Plain verdict: SUCCESS — the decider ran, the bars fired in the
pre-stated direction, and the one gate that failed was diagnosed,
amended pre-quote, and replaced with the certificate the
instrument can actually support.

ELI12: The other team kept only the 1,400 "cleanest" star pairs
and saw nothing — we already showed their clean subset is too
small to see anything. Today we tested the other side: we took
ONLY the 12,900 pairs they threw away and asked our full
machinery — which knows about hidden companions, noise, and every
smudge we've measured — does the extra speed live here, or is it
all just dirt? Answer, four times out of four: the boost is
there, big, and the machinery votes that the "dirt" level is
exactly the modest amount we measured photometrically — it
refuses to blame invisible companions even when allowed to.
Forcing "it's all dirt" costs the fit a fortune every time. And
the nine cliff-guarded pairs that anchor our cleanest evidence?
At least four of them are in the discard pile by their distance
rules alone. So the argument "their cleaning found nothing" now
reads: their cleaning threw away the signal — measured, not
asserted. One of our own test alarms turned out to be testing an
impossible thing; we wrote down why before looking at any answer,
and replaced it with the honest version.

## Review round (2026-07-31): the fresh-Opus blind two-phase review — grades, blind-wishlist validation, TWO ADOPTIONS (the credence split; T6 → Stage 8F)

The external Opus reviewer's context was lost. At the author's
request a FRESH Opus-class agent referee was commissioned under a
deliberately blind two-phase protocol: phase 1 = the v3.9-freeze
state (pre-rivals-arc) graded blind — no tools, no repo access,
briefing text only, all corrections/retractions/credences included,
explicitly invited to attack the briefing itself; phase 2 = the
rivals arc revealed, update requested. Full verbatim record:
REVIEW-FRESH-OPUS.md (untracked, never committed; the old
OPUS-NOTE.md draft is superseded — no recipient exists).

Grades blind → post-reveal: binary anomaly 5/10 → 6/10 ("soft +1":
8D a genuine win vs the strongest external challenge; W4 untouched
and inherited); galaxy program 7/10 → 8/10 ("clear winner of the
campaign" — it re-derived the 8A additive-class theorem BY HAND
before crediting it, and called 8B "a positive control for the c₁
instrument"); cross-system split + AMB 4/10 → 4/10 (hold: the
fence-conditional demotion of the binary veto leg and the 8C
ceiling survival offset); methodology 8/10 → 8.5/10 ("textbook").

BLIND-WISHLIST VALIDATION (the protocol's payoff): with zero
knowledge of the arc, its phase-1 test list contained T4 = "run
the full machinery on the Cookson-removed complement" = Stage 8D
verbatim; T1 = the landed-posterior function re-contest = 7J-d;
and its W7 ("only a₀ ∝ H(z) makes the lock discriminating") = the
8E axis. Independent confirmation the arc was aimed at the right
targets. Genuinely new items from it: T2 external-RV cross-check
(APOGEE/GALAH/LAMOST), T6 fat-tail injections ("cheapest decisive
move"), T3 empirical PM-tail characterization, T5 gas-dominated
c₁, T9 second sim family + mock rotation curves, T10 the joint
a₀(z)–p(z) pair-lock instrument, corner-plot release.

POST-CAMPAIGN ATTACK SURFACE (its ranking): W4 — the fpm = 3.0
width chase read as a fat Gaia relative-PM error tail (~40%) —
now "unambiguously #1: every other binary alternative was
stress-tested this round; this one was conceded untouched, and 8D
rides the same edge." W1 Cassini/MI (~50% MG-excluded) unchanged;
new W9 (the 8E read is favorable-lean on systematics-dominated
external data, ~40%) and W10 (the pincer is only as sharp as the
c₁ window edges; c₁ = 0 and c₁ = 1 exclusions robust, the ¼-vs-½
digit not, ~30%).

**ADOPTION 1 (author delegated the pick "do whatever you think is
best" — executed now): THE CREDENCE SPLIT.** Its §5 correction
adopted verbatim. Two ledgers henceforth: (i) BINARY
"anomaly-real" ~50% (HELD; 8D feeds here; the named decider class
= W4-grade tests — 8F below is the first, T2/T3 the external
legs); (ii) GALAXY-SIDE "the RAR function is BE-family/
thermal-class rather than additive/DM-like" INITIALIZED ~65%
(referee: 65–70; conservative end taken; 8A/8B/8C/8E feed here;
deciders = T5, T9, T10). RULE: no cross-ledger lifting — a galaxy
win never moves the binary number by association, nor vice versa.
The 8D-entry routing "8A/8B/8D flag into the next named decider
together" is SUPERSEDED by this split.

**ADOPTION 2: T6 commissioned as Stage 8F** (pre-reg below).
T2/T3/T5/T9/T10 + corner-plot queued in TODO; all other program
state unchanged; PAPER stays frozen at v3.9.

## Stage 8F PRE-REGISTRATION (2026-07-31, committed BEFORE any run): THE FAT-TAIL ARMS — the W4 manufacture test (referee T6)

QUESTION: if the sky's TRUE relative-velocity errors carry a fat
tail the operative fitter does not model (the photow3 model =
Gaussian single-scale fpm + per-system sq only), does the fitter
MANUFACTURE a phantom α on a Newton sky — and does that world
reproduce the sky's own symptom (the P(fpm = 3.0) chase)? This is
the referee's "cheapest decisive move" on W4 and is fully in-hand:
the 7J-z6 ARMW machinery already carries the fat-tail TRUTH model
(ftl = the fraction of pairs at KT_TAIL = 4× the core error,
per-pair via the ut stream, forward_pp line-level) and the WTRUTH
env interface; the operative fitter reads it with WSHAPE unset.
ZERO model-code changes.

CONTEXT HONESTY (pre-stated): 7J-z6 already offered the SKY a
fitter-side tail axis and the sky priced it at +0.0 (axis unused,
preferring flat fpm = 3.0) — the KT = 4 mixture is therefore
already disfavored AS the sky's own width explanation. 8F measures
the thing never measured: the manufacture RISK (instrument
calibration under tail mis-specification). Both outcomes
pre-stated: manufacture-fired = the α instrument is tail-breakable
(W4 viable); manufacture-clean = tail-robust, W4's α-arm narrows
to non-tail width shapes.

DESIGN (all injections SAMPLE = fullarmw, FPME = 1, photow;
companions MODEL-MATCHED flat-q at the pinned rate 0.10 — the tail
is the ONLY truth-model mismatch by construction; the companion
mismatch direction is already covered by the 7J-z5 A-arm):
- F0 CALIBRATION SWEEP (seed 31; each run fits BOTH laws): Newton
  truth, core fpm_true = 1.2 (the Lindegren-compliant "honest Gaia
  core" hypothesis), sq_true = 0 (the tail must manufacture the
  ENTIRE width demand), ftl ∈ {0.05, 0.10, 0.20, 0.35} → ARMTAGs
  t05/t10/t20/t35; plus one control member ftl = 0.10 with
  sq_true = 0.2 → ARMTAG tsq. WTRUTH =
  'simple,0.00,0.10,1.2,{sq},0,{ftl}'.
- SYMPTOM-MATCH RULE (locked): P₃ ≡ the LANDED-CONV posterior mass
  on the fpm = 3.0 node, averaged over the two law reads, seed 31.
  Calibrated member = the pure-tail member minimizing |P₃ − 0.75|
  (sky: 0.54/0.97). Symptom-matched requires P₃ ≥ 0.30; if no
  member reaches it → CHASE-UNREPRODUCIBLE-BY-TAILS (KT = 4 class)
  is logged and the verdict leg runs on the max-ftl member,
  labeled NON-MATCHED.
- F-N VERDICT (the manufacture test): the calibrated member's own
  seed-31 cubes ARE the Newton+tail verdict reads; plus a seed-101
  confirmation run of that member.
- F-B / F-C RECOVERY: boost truths at the calibrated ftl* —
  'simple,0.74,0.10,1.2,0.2,0,ftl*' and 'BE,0.70,0.10,1.2,0.2,0,
  ftl*' (seed 31; sq_true = 0.2 as in the 7J-z5 arms; ARMTAGs
  fb/fc) — does the boost survive tail mis-spec?
- All reads at LANDED-CONV (operative) + FLAT (cross-check), the
  stage8d_read.py construction verbatim; reader =
  calcs/stage8f_read.py.

BARS (locked now):
- B-MAN (MANUFACTURE-FIRED): calibrated F-N α_marg ≥ 0.5 AND
  dN ≥ +10 at LANDED-CONV in ≥ 1 law, either seed → the operative
  α instrument is tail-breakable at sky-symptom severity. CREDENCE
  MAP (binary ledger): anomaly-real ~50% → ~40% (not lower: the
  7J-z6 sky-side tail rejection, the census cliff, and the
  PHYS-envelope direction all still argue the sky's width is not
  this tail; α rows annotated tail-conditional).
- B-CLEAN (MANUFACTURE-EXCLUDED): F-N α_marg ≤ 0.3 in both laws ×
  both seeds AND F-B/F-C recover within ±0.25 of truth → the
  KT = 4 tail class cannot manufacture the boost. CREDENCE MAP: if
  symptom-matched → ~50% → ~55% (capped: KT = 4 is a
  representative, not exhaustive; T2/T3 remain the external legs);
  if NON-MATCHED → HOLD ~50% and log the W4 narrowing (the tail
  class cannot even produce the chase) as a finding, not a
  credence event.
- Otherwise UNRESOLVED-CARRIED; hold; report per the locked
  grammar.

GATES:
- G8F-WT (config certificate, reader-side): each config's OUT must
  contain the WIDTH-SHAPE injection line matching the registered
  WTRUTH exactly, plus the new injected-histogram checksum line.
  PRE-FLIGHT: the ARMTAG namespace is empty on disk (verified:
  zero matches) — the stale-cube resurrection hazard (the 7D
  lesson: the exists-check would silently reuse another config's
  cubes under the shared legacy name) is closed BY NAMING.
- WIRING (committed with this pre-reg; the 25b248e precedent):
  (1) ARMTAG env — a naming-only suffix on TAG/OUT/cube paths,
  asserted to ARMW with WSHAPE/FUNCS/QLAW off; no model path reads
  it; (2) an additive weighted-checksum print in the ARMW
  injection block (makes injection identity checkable across
  configs and reruns). Amendment protocol standard: any gate/bar
  change is logged here BEFORE results are quoted.
- QUEUED (not this stage): 8F-b = the census-side tail null (does
  a tail-truth sky populate the band WITH the cliff?) on the 7K-b
  forward machinery; the corner-plot release.

EXPECTATION (stated, non-binding): B-CLEAN with the chase only
partially reproduced (the 7J-z5 arms showed Gaussian truths do NOT
chase; a 10–20% tail at 4× should raise P₃ but plausibly not to
0.75) — i.e., the likely landing is "manufacture excluded, chase
partially tail-attributable, W4 narrows to width shapes beyond the
KT = 4 mixture." If B-MAN fires instead, that is the program
finding its own anomaly's source — logged as the better outcome
for truth, per the discipline.

## Stage 8F EXECUTED (2026-07-31): THE FAT-TAIL ARMS — manufacture EXCLUDED 12/12; recovery SPLIT (F-B pass / F-C missed HIGH) ⇒ UNRESOLVED-CARRIED by the letter; credence HELD ~50% per the map

Batch: 8 fullarmw ARMTAG runs (5 F0 configs seed 31 + t35 seed 101
+ fb/fc; ~3 h GPU); G8F-WT config certificates PASS 7/7 (every
injection line float-matched to the registry + checksums). All
numbers: data/stage8f_read.txt.

THE MANUFACTURE LEG (the referee's W4 kill scenario) — DECISIVELY
CLEAN:
- **α_marg = 0.00, dN = +0.0 in 12/12 Newton-truth reads** (ftl
  0.05/0.10/0.20/0.35 + the sq control; both laws, both anchors;
  the t35 verdict member 4/4 across seeds 31/101). The operative α
  instrument manufactures NOTHING from unmodeled error tails at
  any severity up to one pair in three carrying 4× errors.
- The strongest single line: **seed 101's t35 realization
  REPRODUCES the sky's chase (P₃ = 0.51 vs sky 0.54/0.97) and
  STILL reads α = 0.00** — a fake world wearing the sky's own
  noise symptom produces no boost.
- CHANNEL SEPARATION measured: pure-tail skies fit sq = 0 exactly
  (P(sq) = [1,0,0,0]); the sq_true = 0.2 control recovers sq = 0.2
  exactly; fcomp reads the injected 0.10 everywhere. Error tails
  do NOT masquerade as the width channel or as companions — the
  sky's demanded sq = 0.2 is not a disguised KT = 4 tail.

THE SYMPTOM LEG — REALIZATION-DEPENDENT (the honest correction to
my own F0 interim wording): the locked seed-31 rule printed
CHASE-UNREPRODUCIBLE (max P₃ = 0.072 at t20; non-monotone — at
ftl = 0.35 the fitter re-centers its global scale instead of
chasing) and by the letter the NON-MATCHED label stands; but seed
101's t35 chase (0.51) shows the narrowing claim must be SOFTENED:
heavy tails CAN produce a sky-grade chase in some realizations.
The claim that survives is the stronger one: chase or no chase, no
manufacture.

THE RECOVERY LEG — SPLIT, and the miss is a finding:
- F-B (simple truth 0.74 + tail): recovered 0.67 — PASS (±0.25).
- F-C (BE truth 0.70 + tail): own-law read **1.33 = MISSED HIGH by
  +0.63** (the same sky's simple-law read: 0.73; P(fpm = 3.0) =
  0.00 — the BE fit left the noise axis quiet and converted the
  tail width into amplitude). NAMED: **TAIL→α ANTI-CONSERVATIVE
  COUPLING (BE arm, heavy-tail regime)** — on boosted skies with a
  heavy unmodeled tail the BE amplitude can inflate. An
  amplitude-ACCURACY caveat, not an existence risk (Newton skies
  stay at zero), consistent with BE's softer α-conditioning seen
  before (8D complement BE 0.85/0.81 above the full band; 3P
  flat-top α̂_BE profiles).

VERDICT (locked bars): B-MAN not fired (no manufacture anywhere);
B-CLEAN not fired (F-C outside 0.70 ± 0.25) ⇒ **UNRESOLVED-CARRIED;
binary anomaly-real HELD ~50%** exactly per the pre-stated map.
EXPECTATION MISSED and logged: I pre-stated B-CLEAN; the F-C bar
caught a real coupling instead — the bar structure did its job.

CONSEQUENCES:
- The referee's W4 DETECTION arm (bad errors → fake discovery) is
  **EXCLUDED at KT = 4 grade** — the strongest instrument-
  validation statement the binary pipeline has (12/12 zeros, the
  symptom-wearing realization included).
- The W4 AMPLITUDE arm is now MEASURED and OPEN: worst-case +0.63
  (BE, 35% @ 4×; simple −0.07). If the sky carried a heavy real
  tail, α̂_BE would read high — noted alongside the opposing
  sky-side evidence (7J-z6: the sky declines a tail axis at +0.0;
  the PHYS envelope reads α HIGHER when noise is capped).
- Successors queued (post-hoc-flagged, unopened): 8F-c = the bias
  curve at realistic severities (fb/fc at ftl 0.05–0.20 — does
  the BE coupling vanish where tails are plausible?); 8F-b = the
  census-side tail null; T2/T3 = the external legs that would
  close W4 outright.

Plain verdict: NEEDS REFINEMENT by the locked grammar
(UNRESOLVED-CARRIED — the F-C recovery bar missed), with the core
question answered clean: fat error tails cannot fake the
detection; what they CAN do (inflate a real BE amplitude when
severe) is now a measured, named caveat with its follow-up queued.

ELI12: We ran the sting to the end. Fake universes with lying
measurements never produced a fake discovery — twelve out of
twelve said "no extra gravity," including one fake universe that
imitated the real sky's suspicious noise habit almost perfectly.
The software also cleanly told lying measurements apart from
genuinely jittery stars, which were the judge's two look-alike
suspects. One surprise: when we planted REAL extra gravity into a
heavily lying universe, one of the two gravity flavors over-read
its strength (1.33 where we planted 0.70). So bad measurements
cannot CREATE our discovery — but if they were severe, they could
EXAGGERATE its size in that one flavor. We wrote that down as a
new warning label, kept our confidence at exactly 50/50 because
that is what our pre-written rules said to do, and queued the
follow-up that checks whether the exaggeration disappears at
realistic error levels.

## Stage 8F-c PRE-REGISTRATION (2026-07-31, committed BEFORE any run): THE BIAS CURVE — does the TAIL→α coupling operate at realistic severities?

QUESTION: 8F measured a tail→α anti-conservative coupling on the
BE arm at heavy severity (truth 0.70 read 1.33 at ftl = 0.35,
KT = 4) while simple recovered (−0.07). The coupling matters for
the sky only if it operates at severities the sky could plausibly
harbor. 8F-c maps the own-law bias b(ftl) = α̂ − α_truth on
boosted-truth skies at ftl ∈ {0, 0.05, 0.10, 0.20} for both laws
(the 0.35 end already measured = 8F's fb/fc, joined not rerun;
the ftl = 0 members anchor the config's own baseline).

DESIGN: 8 new fullarmw ARMTAG runs, seed 31 (arm-suite precedent;
single-injection grade pre-stated), core fpm_true = 1.2,
sq_true = 0.2, model-matched flat-q companions at 0.10 —
identical to 8F's fb/fc except the ftl dial:
- fb00/fb05/fb10/fb20 = 'simple,0.74,0.10,1.2,0.2,0,{ftl}'
- fc00/fc05/fc10/fc20 = 'BE,0.70,0.10,1.2,0.2,0,{ftl}'
Reads at LANDED-CONV (operative) + FLAT; the OWN-LAW α̂ carries
the bars; cross-law reads, P₃, sq/fcomp portraits = diagnostics.
Reader: calcs/stage8fc_read.py (the 8F construction; G8F-WT
float-compare certificates; the 8F fb/fc rows joined as the 0.35
end).

BARS (locked now; verdict on the REALISTIC-zone members
ftl = 0.05 and 0.10, own-law, LANDED-CONV):
- C-VANISH: |b| ≤ 0.15 at both realistic members for both laws →
  the coupling is a heavy-tail-only phenomenon; the sky-relevant
  amplitude risk is bounded at the systematic scale.
- C-PRESENT: b_BE ≥ +0.25 at either realistic member → the
  coupling operates where tails are plausible. CONSEQUENCE =
  ANNOTATION not credence: BE α rows carry a tail-exposure caveat
  with the measured bound, and the simple-law amplitude is named
  the better-conditioned number (consistent with its clean
  recovery).
- Otherwise GRAY-CARRIED: the curve itself is the product, quoted
  with both flags.
- EDGE RULE: any deciding member within 0.05 of its threshold →
  a seed-101 confirmation of THAT member is REQUIRED before the
  verdict is quoted (pre-registered conditional extension; on
  confirm the member's b = the two-seed mean, bars re-applied; no
  other extension permitted).
- BASELINE RULE: if |b(0)| > 0.15 for a law (the config's own
  baseline bias), that law's bars operate on Δb = b(ftl) − b(0)
  (attribution mode), disclosed in the verdict line.
CREDENCE MAP: NO credence move in ANY branch — this is an
instrument-characterization stage; the existence question was
8F's and held ~50%; C-PRESENT changes annotations only.

GATES: G8F-WT reused (registry float-compare + checksum lines);
pre-flight = the fb00/fb05/fb10/fb20/fc00/fc05/fc10/fc20
namespace empty on disk; amendment protocol standard.

EXPECTATION (stated, non-binding): simple vanishes everywhere
(|b| ≤ 0.1); BE genuinely uncertain between C-VANISH and GRAY
(the 0 → +0.63 interpolation shape is unknown; point guess
b_BE(0.10) ≈ +0.1–0.25); C-PRESENT would be the surprising
outcome and would make the simple-law α the operative amplitude
quote going forward (annotation-grade).

## REVIEW ROUND 3 (2026-08-01 canonical; consult run 2026-07-31 night): W4 40% → 25%; the ranked queue; two derivables; three 8F-c patches demanded

The fresh-Opus referee, shown 8F + 8F-c in full (record:
REVIEW-FRESH-OPUS.md round 3): **W4 updated ~40% → ~25%** with a
floor ~20% ("injection is model-internal — it can show the
pathologies you wrote down don't fool the fitter, never that the
true error distribution is among them"). DEAD for the tested
class: manufacture-from-Newton; chase-as-manufacture-evidence
(the fc10 benign reading); tail-masquerade into sq/companions.
ALIVE: the BE amplitude coupling ("simple is the conservative
anchor — weight the higher BE reads down"); CENSUS-side leakage
(untested = 8F-b); un-tested error shapes (the floor); the
width-object identity (7J-z6 leans astrophysical). Ranking: #1
8F-b census tail null (pre-registered prediction: a tail world
reproduces the band count but FAILS the cliff); #2 register the
DR4 fc10 discriminator NOW (time-sensitive); #3 the graded
flag-correlation width hunt; #4 T2 (scout: marginal — archival
RV precision 0.2–2 km/s vs the 0.1–0.6 km/s signal; first leg =
overlap census); #5 corner-plot. Derivables: D1 analytic
band/cliff leakage under a KT-tail (generalizes 8F-b, removes
realization scatter); D2 the ∂α/∂tail-width BE-vs-simple
derivation (predict the onset analytically). Patches demanded
before leaning on 8F-c: H1 multi-seed the −0.21 simple baseline
("no α from this config quoted better than ±0.2 until
diagnosed"); H2 multi-seed the realistic-zone curve; H3
seed-average the onset (0.15/0.20). Scope phrasing adopted as
standing rule: manufacture-exclusion claims read "for symmetric
variance-inflation tails"; directional/scan-correlated shapes
and tail×companion interaction remain untested (closable only by
T2/T3-class external data). H1–H3 QUEUED (GPU ~3 h, next
session); its W4 number is ITS ledger — ours held per the maps.

## Stage 8F-b PRE-REGISTRATION (2026-08-01, committed BEFORE any run): THE CENSUS TAIL NULL (referee #1; the D1 derivation as the instrument) + P9 REGISTERED

QUESTION: can a fat-error-tail Newton world produce the operative
census pair (band = 9 in the Newton-forbidden [√2, 1.67), cliff =
2 beyond)? This is W4's last in-pipeline angle against the
program's model-light flagship ("the 9 pairs could be the 9
fattest error tails"). The instrument is ANALYTIC (the reviewer's
D1): the 4J T2b leakage null generalized to a per-pair mixture
kernel — no realization draws, no GPU.

DESIGN ([calcs/stage8fb_censustail.py](calcs/stage8fb_censustail.py)):
- Data + masks verbatim from stage4j_gamma82.py (γ ≥ 75 top
  column, WIDE s ≥ 6 kAU; primary convention S/N > 3 = the 3L/v7
  census convention; no-cut variant reported).
- PER-PAIR estimator (primary): every sub-edge pair i (0.2 ≤
  vt_i < √2) contributes leakage probability P_i(window) =
  (1−ftl)·Pwin(vt_i, σ_i) + ftl·Pwin(vt_i, KT·σ_i), σ_i = the 4J
  per-pair sig_vt; μ_band = Σ P_i([1.414, 1.67)), μ_cliff =
  Σ P_i([1.67, 2.2)); joint P(ftl) = P(N_b ≥ 9)·P(N_c ≤ 2) under
  Poisson. Scan ftl ∈ {0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30,
  0.40, 0.50} × KT ∈ {2, 4, 8} (KT = 4 carries the bars; 2/8 =
  D1 generality diagnostics).
- ESTIMATOR DIRECTION (pre-stated): observed-as-truth
  double-counts noise → OVER-estimates leakage → conservative
  toward EXPOSURE; using the corrected count 9 (not the raw 11)
  is likewise exposure-conservative. Both bias AGAINST the
  census defense.
- The 4J-structure median-σ variant retained with the mixture for
  continuity.

GATES: G0 IDENTITY — the ftl = 0 median-σ path must reproduce the
published 4J T2b numbers (parsed from data/stage4j_gamma82.txt if
on disk; both edges × both S/N conventions; else recomputed and
SKIPPED-disclosed). G1 MC — per-pair Bernoulli draws (2e5 reps at
three severities) must bracket the Poisson P's within MC error.
G2 — the per-pair ftl = 0 variant must agree with the median-σ
variant at order level (the per-pair refinement is D1's upgrade,
not a new physics choice).

BARS (locked now; KT = 4, S/N > 3 primary):
- B-EXPOSED: max over ftl of P_joint ≥ 1e-3 → the census is
  tail-fakeable; the census annotation flips EXPOSED; the
  reviewer is informed next round (his W4 climbs back).
- B-IMMUNE: max over ftl of P_joint ≤ 1e-5 AND (at the smallest
  ftl with μ_band ≥ 9, if any exists: μ_cliff ≥ 6 = 3× the
  observed 2 — the monotone-kernel mechanism line) → the census
  is tail-immune at the KT = 4 class; annotation CENSUS-IMMUNE
  (class-scoped per the standing phrasing rule).
- GRAY between. CREDENCE MAP: NO move in ANY branch (defense/
  annotation stage; pre-stated).
EXPECTATION (non-binding): B-IMMUNE via the cliff — the mixture
kernel is monotone, so any ftl large enough to push μ_band
toward 9 should overshoot μ_cliff far beyond 2 (the reviewer's
own pre-stated prediction; also the 7K-b Newton-flood precedent).

P9 REGISTERED in PREDICTIONS.md in this same commit (the DR4
width-object discriminator from the fc10 fingerprint;
time-sensitive — before DR4 exists; instrument-level, not a
gravity claim).

AMENDMENT 1 (2026-08-01, logged BEFORE any verdict is quoted;
the first run's Poisson-based verdict line is NOT quoted and is
preserved in the output as the failed-gate record): **G1 FIRED —
the Poisson approximation to the Poisson-binomial count is
ANTI-CONSERVATIVE at the deciding cell** (KT = 4, ftl = 0.50:
Poisson P(N_b ≥ 9) = 1.8e-3 vs MC 3.5e-4, a 5× overestimate far
outside MC error; the Poisson-binomial variance is below its
mean, so the Poisson tail is fat by construction) — and the MC
itself is resolution-limited below ~5e-6 at 2e5 reps. Per the
gate's own arbitration the verdict estimator is replaced by the
EXACT Poisson-binomial (dynamic-programming convolution over the
per-pair probabilities; no approximation), which G1's MC then
brackets at the resolvable cells. BARS UNCHANGED (1e-3 / 1e-5 on
P_joint). Two further pre-quote notes: (i) the B-IMMUNE mechanism
clause resolves VACUOUSLY IN THE DEFENSE'S FAVOR by its own "if
any exists" wording — no ftl in the class reaches μ_band ≥ 9
(max 2.7 at KT = 4, 5.1 at the KT = 8 diagnostic, vs observed 9):
the tails cannot fill the band in expectation at all, a stronger
defense than the pre-registered cliff mechanism; (ii) a
sensitivity DIAGNOSTIC row is added extending the source
population below the registered vt ≥ 0.2 floor (the apo-island
pairs CAN reach the band at KT ≥ 4 in rare fluctuations;
including them INCREASES leakage = the exposure-conservative
direction the pre-reg prefers) — reported alongside, the verdict
stays on the registered estimator, any bar-crossing by the
diagnostic is flagged per the ship-the-risk-axis rule.

## Stage 8F-b EXECUTED (2026-08-01): GRAY-CARRIED by the letter — and the substance: THE BAND IS UNPRODUCIBLE BY THE TAIL CLASS ("the 9 fattest error tails" dead at realistic severity)

Amendment 1 applied (above): the exact Poisson-binomial replaced
the failed Poisson approximation; G1 now BRACKETS (exact
3.397e-4 vs MC 3.340e-4 ± 1.8e-5 at the deciding cell; 1e6
reps). G0: the 4J T2b code path recomputed with the published
output's T2b lines printed verbatim alongside for the diff. All
numbers: data/stage8fb_censustail.txt (KT = 4 primary, S/N > 3,
210 sub-edge sources).

- **μ_band NEVER exceeds 2.7 (vs observed 9) at any severity to
  50%** — the band count is unreachable IN EXPECTATION by
  symmetric variance-inflation tails. The pre-registered cliff
  mechanism never engages (μ_cliff ≤ 0.33) because the band
  fails FIRST — the defense is STRONGER than the predicted
  mechanism (expectation missed in the informative direction,
  logged: the reviewer and I both predicted band-fills-cliff-
  overshoots; reality is band-never-fills).
- EXACT P_joint: **8.5e-8 at the fc10-fingerprint severity
  (ftl = 0.10)**; 2.9e-6 at 0.20; 3.4e-4 at the absurd class
  edge (0.50 = half the sample at 4× formal). The 4J leakage
  null (3.8e-9 at formal errors) degrades only to ~1e-7 at the
  fingerprint severity — the census remains a ≥ 5σ-grade
  statement even in the tail world that reproduces the sky's
  noise portrait.
- VERDICT by the letter: **GRAY-CARRIED** — the 1e-5 immunity
  bar is missed only by the Poisson-luck residual at the class
  edge (statistic edge-rising at ftl = 0.50; class-scoped,
  disclosed). B-EXPOSED did NOT fire under the valid estimator
  (the first run's Poisson-based B-EXPOSED line is the
  failed-gate record, preserved, not quoted).
- D1 generality: KT = 2 stays ≤ 1e-5-grade everywhere; KT = 8
  (diagnostic) reaches 3.8e-2 at ftl = 0.50 — the defense
  degrades only toward extreme-kurtosis + extreme-frequency
  worlds (8× formal errors on HALF the catalog — visible in any
  whole-catalog σ audit; stated qualitatively).
- Sensitivity (no lower source floor, exposure-conservative):
  identical to print precision — the registered floor does no
  work.
- Census annotation: NOT flipped. The census tail-defense
  statement, class-scoped per the standing phrasing rule: "the
  (band = 9, cliff = 2) pair is unproducible by symmetric
  variance-inflation tails at KT ≤ 4 to 50% severity;
  P ≤ 1e-7 at the fingerprint severity." NO credence move
  (pre-stated in every branch).

CONSEQUENCE: W4's last in-pipeline angle is closed at its
achievable grade — the reviewer's "the 9 pairs could be the 9
fattest error tails" now has a number: ≤ 8.5e-8 at the severity
his own fingerprint world implies. Remaining W4 legs are all
external (T2/T3) exactly as he said. H1–H3 (the 8F-c seed
patches) remain the queued GPU block.

Plain verdict: SUCCESS with the letter GRAY — the instrument
answered the substantive question decisively (the census cannot
be fat-tail leakage at any realistic severity) while the formal
immunity bar proved tighter than the class edge allows; both
stated, nothing reinterpreted.

## Stage 8F-d PRE-REGISTRATION (2026-08-02, committed BEFORE any run): THE PATCH BLOCK — the reviewer's H1–H3 hardening of 8F-c

LINEAGE (no bar-moving): 8F-c's C-VANISH verdict stands AS QUOTED
(single-seed grade, disclosed at quote time). The review round
demanded three patches before that verdict is leaned on (H1 the
−0.21 simple baseline, H2 the single-seed realistic-zone curve,
H3 the at-threshold onset point). 8F-d executes them at two-seed
grade; if the two-seed re-read lands differently, the 8F-c ledger
row flips to SUPERSEDED with a pointer per standing grammar —
never deleted.

RUNS (9 × ~23 min GPU, bar-deciding members first):
- H1 baselines: fb00 + fc00 at seed 101 (two-seed b(0) per law).
- H2 realistic zone: fb05, fb10, fc05, fc10 at seed 101.
- H3 onset: fc20 at seed 101 + a NEW member fc15 =
  'BE,0.70,0.10,1.2,0.2,0,0.15' at seeds 31 AND 101 (the
  reviewer's requested 0.15/0.20 seed-averaged onset bracket).
fb20/fb/fc stay single-seed (curve context, not bar-deciding).

READ (calcs/stage8fd_read.py — successor reader, the shipped
stage8fc_read.py untouched): the 8F-c construction verbatim;
member b = the MEAN over available seeds with the realization
half-spread reported; bars = the 8F-c bars re-applied at
two-seed grade (C-VANISH |b_eff| ≤ 0.15 at ftl 0.05 AND 0.10
both laws / C-PRESENT b_BE ≥ +0.25 at either / GRAY between;
BASELINE RULE unchanged: |b(0)| > 0.15 → that law's bars on
Δb = b − b(0), with the two-seed b(0)). EDGE handling: a
two-seed deciding member still within 0.05 of a threshold is
quoted with a SEED-LIMITED flag; NO third seed (pre-stated — the
reviewer's H2 asked for multi-seed = 2; no further extension).

H1 DIAGNOSIS GRAMMAR: |mean b(0)_simple| > 0.15 at two seeds →
the config baseline is REAL (not a realization fluke) — the
attribution-mode reading stands and the reviewer's quote-
precision rule (±0.2 for this config) is ADOPTED as a standing
annotation until the mechanism is diagnosed; ≤ 0.15 → fluke
component confirmed, raw-b bars apply, noted.

H3 GRAMMAR: the onset statement becomes "b_BE crosses +0.25
between the two-seed ftl = 0.15 and 0.20 means" (or wherever it
lands), with half-spreads; single-seed language retired.

GATES: G8F-WT per config (registry float-compare + checksums;
fc15 pre-flight namespace EMPTY — verified before launch); GD0
REGRESSION — the 8fd reader's seed-31 own-law â for every
existing member must reproduce the shipped stage8fc_read.txt
values to 0.01 (parse-compare; catches reader-transcription
drift). Amendment protocol standard.

CREDENCE MAP: NO move in ANY branch (hardening stage; even a
C-VANISH → C-PRESENT flip is annotation-only per the original
8F-c map: BE tail-exposure caveat + simple named the
better-conditioned amplitude).

EXPECTATION (non-binding): C-VANISH confirmed at two-seed grade;
b(0)_simple regresses toward ~−0.1 but may persist (genuinely
uncertain — this is H1's point); the onset lands between 0.12
and 0.20 with fc15 near or above +0.15.

## Stage 8F-d EXECUTED (2026-08-02): C-VANISH CONFIRMED-HARDENED at two-seed grade; both config baselines land AT the boundary (±0.16, opposite signs); the onset bracketed 0.15–0.20

Batch 9 runs (~3.5 h GPU); G8F-WT 11/11; **GD0 regression PASS
(all 10 shipped 8F-c seed-31 values reproduced)**. All numbers:
data/stage8fd_read.txt.

THE TWO-SEED CURVE (mean, half-spread):
- simple: −0.16(0.05) / −0.14(0.04) / −0.10(0.02) at ftl
  0/0.05/0.10 (0.20/0.35 single-seed −0.14/−0.07).
- BE: +0.16(0.02) / +0.02(0.08) / +0.06(0.10) / +0.23(0.03) /
  +0.27(0.01) at 0/0.05/0.10/0.15/0.20 (+0.63 at 0.35).

H1 — THE BASELINE DIAGNOSIS: the simple offset SHRANK (−0.21 →
−0.16 two-seed) but sits AT the 0.15 boundary → REAL by the
letter, and BE's baseline crossed it too (+0.14 → +0.16): **both
laws in attribution mode, the ±0.2 quote-precision annotation
ADOPTED — scoped to the 8F-family arm configs** (fpm_true = 1.2,
model-matched companions; the operative sky band 0.68–0.74 rests
on the 7J-z5 own-config arms at −0.06…−0.09, NOT on this
config). Notable structure, logged without interpretation: the
two baselines are SYMMETRIC (±0.16, simple under / BE over) —
a law-dependent own-truth conditioning at this cell, boundary-
grade, neither sharp nor noise (the reviewer's H1 worry lands
exactly on the fence).

H2 — THE VERDICT FOUNDATION: seed-101 moved individual
realistic-zone members by up to 0.20 (fc10: −0.04 → +0.06 mean,
half-spread 0.10) — the two-seed grade genuinely mattered — and
**C-VANISH HOLDS: |eff| ≤ 0.15 at both realistic members, both
laws, in BOTH framings** (raw BE +0.02/+0.06; attributed
−0.14/−0.10). One SEED-LIMITED flag (BE 0.05, attributed −0.14,
within 0.05 of the threshold; no third seed per pre-reg).

H3 — THE ONSET: **b_BE crosses +0.25 between the two-seed
ftl = 0.15 (+0.23 ± 0.03) and 0.20 (+0.27 ± 0.01)** — the
single-seed "at-threshold" language retired. ATTRIBUTION CO-READ
(disclosed, ship-the-risk-axis): with the +0.16 baseline
subtracted, the TAIL-ADDED bias at 0.15/0.20 is only +0.07/+0.11
and reaches +0.47 at 0.35 — the raw onset bracket partly rides
the config baseline; the heavy-tail blow-up is real in both
framings.

VERDICT (locked bars, two-seed): **C-VANISH CONFIRMED-HARDENED —
the coupling is heavy-tail-only; the 8F-c verdict stands,
upgraded from single-seed grade** (8F-c ledger row stays CURRENT;
this row records the upgrade). NO credence move (pre-stated).
EXPECTATION largely met (C-VANISH ✓, onset bracket ✓ incl. fc15
at +0.23; the baseline neither regressed fully nor persisted
hard — boundary, the honest H1 answer).

REVIEWER SCORECARD (for the standing calibration question): H2
vindicated (members moved at the half-spread scale; the verdict
needed the second seed to be leanable); H3 vindicated (clean
bracket); H1 landed ON the fence (the offset shrank but stays
boundary-real — neither its "~0.2 bias nobody understands" nor
my "realization fluke" reading; the measurement split the
difference, and the annotation is now scoped correctly).

Plain verdict: SUCCESS — all three patches measured, the 8F-c
verdict hardened rather than flipped, the annotation scoped, and
the onset now carries an error bar.

ELI12: We re-flipped the three coins the judge worried about.
(1) The weird baseline offset shrank but didn't vanish — it's a
small real quirk of the fake-universe setup, the same size in
both gravity flavors but opposite directions, and it does NOT
touch our real-sky numbers, which were validated on a different,
cleaner setup. We stamped those fake-universe results "good to
±0.2" as the judge asked. (2) The verdict that matters —
"realistic levels of lying measurements do not inflate the
signal" — survived the second coin-flip in both flavors. (3) The
switch-on point of the exaggeration is now pinned between 15%
and 20% lying-fraction with real error bars, far above anything
believable for Gaia. Net: the judge's demands made the numbers
sturdier and changed no conclusions — which is what good
refereeing is for.

ELI12: Our judge's last remaining in-house worry was: maybe the
nine "impossible-speed" star pairs are just the nine unluckiest
measurement errors. We did the math exactly — every slower pair
gets its own tiny probability of faking its way into the
forbidden speed zone, lies included, and we add them all up.
Verdict: even if one pair in ten had 4×-underestimated errors
(the level our best fake-sky needed), the chance of faking nine
is about one in ten MILLION. To fake it you'd need half of all
measurements to be wrong by 4× — cartoon territory. One formal
threshold we set in advance wasn't quite met at the cartoon end
of the dial, so the official label is "gray, carried" — but the
question the test was built to answer came back loud: those nine
pairs are not measurement luck. Also, our own approximation
failed its built-in alarm mid-run and was replaced by the exact
calculation before we read any verdict — the alarm did its job.

## Stage 8F-c EXECUTED (2026-07-31): THE BIAS CURVE — C-VANISH; the coupling's turn-on LOCATED (between 10% and 20% severity); the most sky-like fake world yet built

Batch: 8 configs seed 31 (~3 h, 22.7 min/config like clockwork);
G8F-WT certificates 10/10 (the joined 8F fb/fc ends included).
All numbers: data/stage8fc_read.txt.

THE CURVE (own-law bias b = α̂ − truth, LANDED-CONV, at ftl =
0 / 0.05 / 0.10 / 0.20 / 0.35):
- simple: −0.21 / −0.18 / −0.12 / −0.14 / −0.07. BASELINE RULE
  FIRED (|b(0)| = 0.21 > 0.15): bars ran on Δb = +0.03 / +0.09 in
  the realistic zone — the tail ADDS essentially nothing; the
  baseline under-recovery is the config's own conditioning at
  single-injection grade, and its direction is CONSERVATIVE
  (under-reads, never inflates). Disclosed per the rule.
- BE: +0.14 / −0.07 / −0.04 / +0.25 / +0.63 (raw bars; b(0)
  inside 0.15). REALISTIC ZONE CLEAN (−0.07 / −0.04): at
  plausible severities the BE amplitude recovers truth. The
  coupling TURNS ON between ftl = 0.10 and 0.20 (+0.25 at 0.20 —
  single-seed, numerically AT the C-PRESENT threshold but outside
  the pre-registered deciding zone; quoted as the located onset)
  and cliffs at 0.35 (+0.63).

VERDICT (locked bars; edge rule NOT fired — all deciding members
≥ 0.05 from their thresholds): **C-VANISH — the tail→α coupling
is a heavy-tail-only phenomenon; the sky-relevant amplitude risk
is bounded at the systematic scale. NO credence move (pre-stated
in every branch).** EXPECTATION MET (pre-run: simple vanishes,
BE between VANISH and GRAY — landed VANISH).

THE DIVIDEND — the chase window: P₃ on boosted skies is
NON-MONOTONE in ftl: 0.37 (0.05) → 0.88 (0.10) → 0.01 (0.20) →
0.00 (0.35) on the BE arm. **fc10 (BE truth 0.70 + a 10% tail) is
the most sky-like fake world ever built: it chases the noise edge
at 0.88 (the real sky's BE read: 0.97) while reading α̂ = 0.66 ≈
truth, sq = 0.2, fcomp = 0.10** — a boost plus a MODEST tail
reproduces the sky's entire Gaussian-fitter portrait with the
amplitude uncorrupted. CAVEAT carried both ways: this is a
FINGERPRINT for the width-object hunt (TODO-18) — "the sky's
chase is what a ~10% tail looks like to this fitter" — NOT a
resolution: 7J-z6 offered the sky the explicit tail axis (KT = 4,
ws up to 0.15) and the sky priced it +0.0; that standing
rejection argues against the literal-tail reading, and both
statements are quoted together.

CONSEQUENCES: the 8F amplitude caveat is BOUNDED — the BE
annotation downgrades from "open exposure" to "bounded: onset
requires ≥ ~15–20% severity at KT = 4"; the operative α band
0.68–0.74 and the 8D complement reads carry no live tail-bias
exposure at this class; W4's remaining live legs are the external
ones (T2 external velocities / T3 PM-tail measurement), plus the
8F-b census-side null, all queued.

Plain verdict: SUCCESS — the curve measured, the verdict by the
letter, the coupling's onset located, and the stage returned a
bonus fingerprint for the next hunt.

ELI12: The exaggeration is a cliff, not a slope. At believable
levels of measurement-lying (one pair in ten or twenty), both
gravity flavors read the planted signal correctly — the
exaggeration we caught yesterday only wakes up when a fifth of
all measurements lie by 4×, and only becomes huge at a third. So
yesterday's warning label gets a bound: "applies only in
cartoon-noise worlds." Bonus: the fake universe with a real
signal plus 10% lying measurements behaves almost exactly like
our real sky — same suspicious noise appetite, same healthy
signal reading. So the sky's odd noise habit might literally be
telling us "about one pair in ten has worse errors than
advertised" — and even if that's true, the signal measurement
doesn't care. One older test pushes back on that tidy story, so
it goes in as a lead for the hunt, not a conclusion — both sides
written down.

## Stage 8G PRE-REGISTRATION (2026-08-03, committed BEFORE any run): THE E-SECTOR CONTROL — the width-object program, leg 1 (mundane-first)

CONTEXT (the program, author-adopted this session): treat the
sq = 0.2 width object as a candidate hidden variable instead of
dirt (the author's Einstein-flip framing). Three legs, strictly
ordered: leg 1 = THIS (free the 7J-z6-booked e-sector — the
mundane candidate goes FIRST); leg 2 = the boundedness contest on
the (band=9, cliff=2) census pair (8H, designed only AFTER leg 1
lands — its premise and bars depend on this verdict); leg 3 = the
second-moment cross-system derivation (alongside, no pipeline).
EXPECTED OUTCOME, pre-stated per the published recommendation:
E-ABSORB (the mundane branch). An E-SURVIVE verdict fires AGAINST
the stated expectation and must be flagged as such.

WHAT IS FROZEN TODAY (the fingerprint's named sector,
stage7j_marginal e_of): the power-law arm's inner anchors are
HARD-CODED (al = 0.6 at 100 AU, 1.0 at 500 AU; only the >= 1 kAU
exponent eta is on a grid, 2-noded [1.05, 1.3],
Hwang-prior-anchored); the radial component's floor is HARD-CODED
(e_rad = 0.90 + 0.095 u). 7J-z6's Part-A fingerprint
(mid-shoulder ṽ 0.7–1.7 + radial γ=8° column + inner bins) points
at exactly this sector; sq's existence is measured 4× (3E / 6P /
the −60..−116 misfit / the fpm edge) but no traced channel makes
20%.

DESIGN (mode ESEC=1 on the operative landed grid, full photow
FPME=1): two new axes inserted after wr —
  EIN_GRID = [0.5, 1.0, 1.5, 2.0] — multiplier on the frozen
    inner anchors: al = interp(log10 a_s at [100, 500, 1000,
    50000] AU -> [0.6·ein, 1.0·ein, eta, eta]); ein = 1 = legacy,
    INTERIOR.
  ERF_GRID = [0.80, 0.90, 0.95] — the radial floor: e_rad = erf +
    (0.995 − erf)·u; erf = 0.90 = legacy, INTERIOR.
Wide-arm eta stays 2-noded Hwang-anchored — SCOPED freedom,
disclosed: this control frees the fingerprint's named sector
only. The identity cell (1.0, 0.90) is BRANCHED to the legacy
e_of verbatim (the ws=0 precedent). Cubes:
data/stage7j_cube_full_esec_{seed}_{law}.npy, 11-dim (A5 E2 WR5
EIN4 ERF3 FCOMP6 FC0·1 FFLY2 FPM6 KW3 SQ4); fresh names (the
ARMTAG/7D stale-cube rule). GB0w/GB0e SKIPPED-disclosed under
ESEC (G8G-0 substitutes and is strictly stronger). Seeds 31 +
101, bars on seed means (the amendment-5d standard); ~1080 GPU
cells per seed, estimate ~4–5 h/seed (non-binding).

GATES (abort-grade before any quote):
  G8G-0 (in-run, per law×seed): the (ein=1.0, erf=0.90) slice of
    the ESEC cube AND cubevt must equal the operative photow3
    cube EXACTLY — max|Δ| <= 1e-9, target 0.0e+0 (the GB0w
    new-mode standard; the identity cells are RECOMPUTED through
    the new loop, so this is an end-to-end regression of the
    restructured build).
  G8G-1 (reader): the identity-slice LANDED-CONV read must
    reproduce the shipped per-seed E-arm EXT rows (parsed from
    data/stage7jz5_eread.txt at runtime: simple 0.63/+16.2 and
    0.74/+12.7, BE 0.75/+12.8 and 0.59/+19.3) to |Δα| <= 0.01,
    |ΔdN| <= 0.1.
  G8G-2 (reader): finite-mass bookkeeping + EDGE report on the
    new axes (correction-#4 standard: freed-posterior mode AT a
    new-axis edge with >= 0.5 mass -> EDGE flag; extension is a
    decision, not an auto-run).

THE READ (calcs/stage8g_read.py; LANDED-CONV anchor + prior_eta,
the read6 convention extended to 11-dim): per law×seed — the
identity-slice read (α, dN, P(sq)) vs the freed read (α, dN,
P(sq), P(ein), P(erf)), the profile gain G_e = max lnL(freed) −
max lnL(identity), and the marginal ΔdN. P(sq) = the sq-axis
posterior marginalizing everything else including α.

BARS (locked; seed means, both laws; sq mode = argmax of the
seed-mean posterior):
  E-ABSORB:  P(sq > 0) <= 0.50 both laws -> the width demand is
    absorbable by the named e-sector; sq = e-sector-inadequacy
    CONFIRMED at control grade; the freed α/dN are reported as
    conditional candidates (promotion = a review-round question,
    NOT this stage's verdict).
  E-SURVIVE: P(sq > 0) >= 0.90 AND sq mode >= 0.1, both laws ->
    the width object survives the freed e-sector; leg 2 (8H
    boundedness) becomes the decider. AGAINST-EXPECTATION flag
    mandatory.
  else E-PARTIAL: the decomposition is the product; report
    per-law per-seed as-is, no interpretation beyond the printed
    numbers.
  SEED-SPLIT: seeds landing in different branches per law ->
    report as-is, verdict = the seed-mean branch with the split
    flagged; NO auto-extension (cost honesty, disclosed).
  α exposure: Δα = freed − identity (seed mean, per law);
    |Δα| > 0.11 (the standing systematic band) -> MATERIAL flag
    with the pre-scripted sentence: "α is exposed to the e-sector
    at the freed grid; the operative quote gains this as a
    systematic annotation" (no re-quote in this stage).
  NO credence move in ANY branch (pre-stated; the anomaly-real
  deciders remain the W4-external class + the census program).

Amendment rule: grid values and gate tolerances may be amended
pre-quote with logged reasons (the standing rule); bars may NOT
move.

## Stage 8G EXECUTED (2026-08-03): E-SURVIVE — the width object survives the freed e-sector, AGAINST the pre-stated expectation; the inner anchors are data-confirmed; the radial floor tightens to the 0.95 edge; α gains a −0.12 e-sector annotation

Gates first-run clean: G8G-0 bit-exact 8/8 (0.00e+00 including the
vt cubes, both seeds both laws — the restructured loop IS the legacy
code at the identity cell); G8G-1 4/4 (the identity-slice reads
reproduce the shipped E-arm rows to the digit: 0.63/+16.2,
0.74/+12.7, 0.75/+12.8, 0.59/+19.3); G8G-2 fires the EDGE flag
(below). Runtime 274 + 269 min (the 4–5 h/seed estimate borne out).

THE ANSWER: **P(sq > 0) = 1.00 in ALL FOUR law×seed freed reads,
sq-mode = 0.2 everywhere** (the single off-diagonal cell is BE-31's
P(0.1) = 0.01). The 7J-z6-named sector — freed on both its frozen
axes — does NOT absorb the width object. E-SURVIVE fired in both
laws with no seed split, AGAINST the pre-registered expectation
(E-ABSORB was the published bet). The sq exclusion table now reads:
not mass errors (3J, ×12 too small), not companions (7J posterior),
not KT=4-class error tails (8F family), NOT the e-sector (8G). Leg
2 — the 8H boundedness contest on the (band=9, cliff=2) census
pair — is the decider by pre-registration.

THE SECTOR ITSELF (secondary measurements): (a) **P(ein) = 1.00 at
the IDENTITY node in every read** — the published inner-anchor
profile [0.6, 1.0] is data-CONFIRMED, not merely assumed; the
inner-bin arm of the 7J-z6 fingerprint does not want a different
inner power law. (b) The radial floor runs to **erf = 0.95 at
P = 1.00 — EDGE flag** (correction-#4: extension is a decision,
not an auto-run; folded into the 8H design round, where the
tightened near-parabolic column also changes the census-forward
premise). The tighter column (e ∈ [0.95, 0.995]) buys G_e = +18.5
(simple) / +16.1 (BE) lnL at the best cell — a real model
improvement ORTHOGONAL to sq: it does not touch the width demand.
(c) Profile-level only (not a marginal read): the freed PROF rows
take fpm off the 3.0 edge in 3 of 4 rows (2.4/2.4/2.1; BE-101
stays 3.0) — the tightened column may eat part of the noise chase;
the freed P(fpm) marginal was not computed here (a 8H-round item).

THE α ANNOTATION (MATERIAL, the pre-scripted sentence applies):
Δα(freed − identity) = −0.131 (simple) / −0.117 (BE) — just past
the ±0.11 band. Freed α_marg per seed: 0.57/0.53 (simple 31/101),
0.54/0.57 (BE). "α is exposed to the e-sector at the freed grid;
the operative quote gains this as a systematic annotation" (no
re-quote in this stage). Newton under the freed sector: dN
+7.5/+18.4 (simple), +3.2/+20.8 (BE) — rejected in all four reads;
the seed-mean margin softens ~4 lnL from the operative band; the
BE-31 single read at +3.2 is the low-water mark, quoted openly.

NO credence move (pre-stated in every branch; the census program
and the external W4 legs remain the deciders).

Plain verdict: SUCCESS — the control answered decisively, and the
answer is E-SURVIVE against our stated bet: the width object is
real structure that this model family cannot explain away with
orbit shapes.

ELI12: We bet the boring answer would win — that our too-rigid
menu of orbit shapes was faking the 20% jitter. We unlocked the
menu and let the data order anything. The data said: "your inner
orbit shapes were RIGHT all along (they picked the exact published
values back), the plunging-orbit family should be slightly MORE
extreme (a real small improvement), and the 20% jitter? Still
there, at full certainty, in every check." The mystery just
survived its most likely boring explanation — losing our bet the
way that makes the mystery more interesting. One warning label
from the same run: with looser orbit menus the gravity-boost knob
reads ~0.55 instead of ~0.7, so that goes on the quote. Next up:
does the jitter respect the escape-speed ceiling? Physics obeys
speed limits; noise doesn't.

## Stage 8G-b (2026-08-03, read-only diagnostic, frame 983f9d9): the freed e-sector RELEASES the noise chase in 3 of 4 reads — MEASUREMENT ONLY

G8Gb-0 4/4 (the identity-slice P(fpm) reproduces the shipped E-arm
rows to <= 0.005). THE MEASUREMENT — marginal P(fpm = 3.0),
identity -> freed: simple-31 1.00 -> 0.16 (mode 2.4, interior);
simple-101 0.09 -> 0.02 (mode 2.1); BE-31 0.94 -> 0.14 (mode 2.4);
BE-101 1.00 -> 1.00 (the hold-out). Released modes sit at 2.1–2.4
= still above the ~1.4 Lindegren ceiling: the width-shape object
DEFLATES in three reads, does not vanish. Companion measurement:
P(wr) moves 0.20 -> 0.30 in ALL FOUR reads — more weight in the
tighter erf = 0.95 column (the radial channel restructures as a
unit). Numbers only tonight per the cadence rule and the committed
frame; the candidate reading — that the sky's noise-edge chase
(the arm suite's "chase-unexplained" token) was partly the
too-loose radial floor, while the sq width object survives
untouched — belongs to the 8H pre-reg, to be stated and TESTED
there, not assumed here. No verdict, no credence move.

ELI12: A quick bonus X-ray of today's big run. Our fit has always
had a weird habit: it kept cranking the "assume the telescope is
3x worse than advertised" knob to the max. With the orbit menu
unlocked, that habit mostly stopped — in 3 of 4 checks the knob
came down to normal-ish — but one check still cranks it, and
even "normal-ish" is twice the advertised errors. So part of the
weird habit was our orbit menu's fault, part still isn't
explained. What we did NOT do tonight: decide what that means.
That's tomorrow's carefully-worded test.

## Stage 8H PRE-REGISTRATION (2026-08-04, committed BEFORE any run): THE CENSUS SHAPE-AND-ATTRIBUTION CONTEST — the width-object program, leg 2 (the decider named by the 8G pre-reg)

PREMISE (set by 8G/8G-b): the width object survived the freed
e-sector (E-SURVIVE, P(sq>0)=1.00 all reads); the freed model
revision is real (erf=0.95, wr 0.2->0.3, fpm chase released in 3/4
reads). QUESTION: what SUPPORT SHAPE does the width object have in
the census channel, and what actually floods the (band=9, cliff=2)
pair? The shipped 7K-b rows show the overshoot flood is NOT purely
the smear: mu_hi ~ 13 at sq=0 (companion wobble spikes implicated)
vs 17 at sq=0.2, vs 2 observed. So this contest carries BOTH axes:
KERNEL SHAPE and CHANNEL ATTRIBUTION, at the NEW freed-e premise.

DERIVED NOTE (stated pre-run): the coherent per-system boost draw
is INDISTINGUISHABLE from the Gaussian smear in the (band, cliff)
statistic — both multiply vtilde by an unbounded lognormal factor;
a drawn boost moves each system's own ceiling with it, so the
fixed-window counts coincide. The census can therefore bound the
smear channel's SUPPORT (bounded vs unbounded) but cannot separate
draw-from-noise; K3 is retired into K1 by this equivalence.

DESIGN ([calcs/stage8h_censusshape.py](calcs/stage8h_censusshape.py),
the 7K-b machinery + the esec cubes; seeds 31 + 101, both laws):
per seed x law, the FREED MAP cell (argmax of the prior-augmented
11-dim esec cube; its alpha/eta/wr/ein/erf/fcomp/fpm/kw/sq printed)
gets the forward ladder:
  A  fcm=0,  sq=0      — pure orbit + noise base (attribution anchor)
  B  fcm=cell, sq=0    — + companions (K0)
  C  fcm=cell, sq=cell — + Gaussian smear (K1; = the operative)
  K2a = C with the smear BOUNDED at C(alpha_cell) = sqrt(2*(1 +
    alpha_cell*0.36)) (B-1 = 0.36 at the wide field, the 4J edge);
  K2b = C bounded at C(1) = 1.649; K2n = C bounded at C(0) = 1.414
    (the Newton-ceiling control).
  The bound acts on the SMEAR ONLY: vts_b = min(vts, max(vtn, C))
  — companion-carried and noise-carried values above C keep their
  legacy values (a photocenter wobble is real velocity, not
  escape-bounded; pre-registered as the physical choice).
Scoring: mu_band [1.414, 1.67) and mu_hi [1.67, 2.2) at gamma>=75
(the 7K-b windows verbatim, NDATA scaling verbatim); joint Poisson
LL at the observed pair: LL = ln Pois(9; mu_band) + ln Pois(2;
mu_hi). Population-Poisson convention stated (model-side forward
expectations; the 8F-b exact-PB estimator is the DATA-side tool).

GATES (abort-grade):
  G8H-0 (regression): at the LANDED photow3 cells the machinery
    must reproduce the shipped 7K-b rows — boost simple/BE mu_band
    and mu_hi (29.57/17.56, 30.06/18.85, 29.74/18.42, 29.00/17.86)
    and the POST-HOC boost sq=0 fpm=1.5 rows — to |d mu| <= 0.05.
  G8H-1 (arithmetic sanity): the bound can only remove — mu(K2x)
    <= mu(C) per window per config; K2n's band contribution from
    the smear must be zero by construction.
  Consistency note (not a gate): mu_hi at A under Newton vs the
    8F-b analytic leakage scale (~<= 2.7) — different constructions,
    printed for orientation.

BARS (locked; a config is ADMISSIBLE if joint P = Pois(9; mu_b) *
Pois(2; mu_c) >= 1e-3; evaluated per law x seed, 4 runs):
  RESOLVED(config): some ladder config admissible in >= 3 of 4 runs
    -> name it. If UNIQUELY a K2 member: SHAPE-SELECTED (the width
    object has bounded support; the selecting ceiling named — the
    ceiling position doubles as an alpha read, mechanical output).
    If B or A at the freed premise: FREED-MODEL-RESOLVED (the 8G
    revision itself repairs the census; the smear must then NOT
    live in this channel). If multiple admissible: RESOLVED-
    DEGENERATE (report the set, no selection sentence).
  ALL-FAIL: no config admissible in >= 3 of 4 -> the (9,2)
    reproduction problem STANDS at kernel grade; the ATTRIBUTION
    is then the product: share of mu_hi carried by the companion
    sector = (B-A)/C per window; if >= 50% in both laws, the
    census pair is named a COMPANION-SECTOR diagnostic (feeds the
    7J-z4 wobble-binding line), and the width-object reading of
    the census is DEMOTED (the cliff meters the wobble law, not
    the smear).
  MECHANICAL OUTPUT regardless of branch: the observed band spans
    vt 1.490-1.607 (six of nine above C(0.5) = 1.536) — where the
    bounded family's admissibility lands vs alpha is printed as a
    ceiling-position read, not interpreted here.
EXPECTED OUTCOME (pre-stated, honest): ALL-FAIL with companion-
  sector attribution (base rate: no config has ever reproduced the
  pair; the sq=0 rows already sat at P ~ 9e-5). The live upside
  branch — freed-premise admission (the 8G revision repairing the
  census) — is named but NOT expected.
NO credence move in ANY branch (pre-stated). 8H decides the
  width-object SHAPE question, not anomaly-real.

FOLD-INS (logged): the erf-grid extension is DEFERRED with reason
(the contest runs at the measured 0.95 mode; a ~6 h cube run the
band/cliff question does not need; revisit conditional on this
verdict). The 8G-b candidate reading (chase = partly the loose
floor) gets its census-side test embedded in the freed-premise
runs. Amendment rule: gate tolerances amendable pre-quote with
logged reasons; bars may NOT move.

AMENDMENT 1 (2026-08-04, logged PRE-QUOTE after G8H-1 fired on the
first run; verdict was WITHHELD by the gate as designed; first-run
output preserved as data/stage8h_censusshape_run1.txt = the
failed-gate record, never quoted as a result): the pre-registered
CLAMP bound (vts_b = min(vts, max(vtn, C))) creates a delta
pile-up AT the ceiling — it (a) manufactures band counts (the
clamped mass lands inside [1.414, 1.67), mu_band ROSE 28 -> 38,
exactly what G8H-1's bound-only-removes condition existed to
catch) and (b) makes every in-window ceiling equivalent for a
window-count statistic (K2a = K2b = K2n to the last digit = a
DEGENERATE kernel family; the design note had flagged the clamp's
pile-up and chosen it anyway — wrongly). REPLACEMENT (the 8F-b
amendment-1 precedent: a gate fired on instrument internals; the
component is replaced; the BARS do not move): the CONDITIONED-ON-
BOUND smear — per system the log-normal draw is truncated at
g_max = ln(max(vtn, C)/vtn)/sq and RENORMALIZED via the uniform
map g_t = PHI^-1(PHI(g) * PHI(g_max)) (no pile-up; the excess
mass redistributes proportionally below the cap; different
ceilings now differ). G8H-1 is restated to the truncation
kernel's true invariants: (i) cliff-only-removes, mu_hi(K2x) <=
mu_hi(C); (ii) non-degeneracy, the K2 band values must differ
across ceilings; (iii) mu_hi(K2n) <= mu_hi(B) (all smear-carried
cliff mass removed at cap 1.414). Admissibility bar, ladder,
scoring, and verdict grammar UNCHANGED.

## Stage 8H EXECUTED (2026-08-04, amended run, all gates PASS): ALL-FAIL with COMPANION ATTRIBUTION 0.81/0.78 — the census pair is a companion-sector diagnostic; the width-object reading of the census is DEMOTED

Gates: G8H-0 4/4 (landed-cell regression reproduces every shipped
7K-b row to the digit — machinery certified); G8H-1 4/4 on the
amendment-1 truncation kernel (the first-run clamp form was caught
by G8H-1's own bound-only-removes condition, logged pre-quote,
preserved in _run1.txt — the gate wrote the amendment, we did not
move a bar). Runtime 1.2 min.

THE LADDER (freed-e MAP cells, both laws both seeds — the cells
themselves: alpha=0.5, eta=1.05, wr=0.3, ein=1.0, erf=0.95,
fcomp=0.1, kw=0.7 = the GRID FLOOR of the amplitude nuisance,
fpm 2.1–3.0, sq=0.2; observed pair (band 9, cliff 2)):
  A  orbit+noise only:  mu_band 0.77–1.81, mu_hi 0.03–0.21 —
     ESSENTIALLY CLEAN (Newton-A cleaner still, 0.13–0.45 / 0.00,
     consistent with the 8F-b analytic leakage scale);
  B  +companions:       mu_band 13.4–15.4, mu_hi 13.4–14.4 —
     THE COMPANION WOBBLE CHANNEL FLOODS BOTH WINDOWS (×1.6 the
     observed band, ×7 the observed cliff) at the model's OWN
     preferred settings (fcomp = 0.1, kw at the grid floor);
  C  +smear (operative): mu_band 28.3–29.6, mu_hi 16.7–18.2;
  K2a/K2b/K2n (bounded smear, non-degenerate after amendment 1):
     cliff pinned at 8.6–9.2 in ALL variants = the companion-
     carried floor; band runs 24.6–26.8 / 32.1–34.2 / 8.5–9.6.
ADMISSIBILITY: 0/4 runs for EVERY config (best jointP 8.9e-4,
K2n, below the 1e-3 bar). VERDICT BY THE LOCKED BARS: ALL-FAIL —
the (9,2) reproduction problem STANDS at kernel grade; attribution
(companion share of the cliff at C) = 0.81 (simple) / 0.78 (BE) >=
0.5 both laws -> **the census pair is named a COMPANION-SECTOR
diagnostic (feeds the 7J-z4 wobble-binding line); the width-object
reading of the census is DEMOTED.** The pre-stated expected
outcome, landed.

WHAT THIS MEANS (scoped carefully): (1) THE FLOOD IS LOCALIZED —
the model's wobble-spike tail (the photocenter amplitude law ×
the S period-smearing × the log-normal period tail) puts ~13
phantom pairs into the census windows that the sky does not show,
even with the amplitude nuisance riding its grid FLOOR (kw = 0.7)
— the count-level sharpening of 7J-z4's Dwob +314 wobble-binding.
(2) THE CLIFF-METER IDEA (leg 3's "cliff bounds the coherent
draw") DIES HONESTLY: the cliff was never smear-dominated — at
sq = 0 it is already flooded ×7 by wobble spikes; the census
cannot bound the smear's support until the wobble tail is
repaired. (3) SCOPE OF THE DEMOTION: the DATA-SIDE census
defenses are UNTOUCHED (4J leakage null 3.8e-9, 8F-b tail null
8.5e-8 — error-model class, no companion physics in them); what
is suspended is the FORWARD-model gravity-evidence reading of the
band (the model's own account of the 9 pairs is companion-
dominated AND overshooting — broken either way). The 7K-b
"self-defending pair" sentence survives with its scope named:
the pair defends against error tails; it currently READS the
wobble model. (4) The pre-registered ceiling-position/alpha read
returns NO READ (confounded by the companion channel; K2n's
near-admissible (9, 9) at the Newton ceiling decomposes as
down-shuffled companion flood, stated with that caveat, not
selected — the bar held). (5) NEXT LEVER, named: the wobble-TAIL
instrument — the period/amplitude tail of the companion model is
now constrained by THREE independent strains (7J-z4 kinematic
binding, round-10 q-moments, the 8H census counts); repairing it
is the gate to re-opening the census as a width meter AND to the
fifth-move exposure chain. External legs (T2/T3, P9/DR4)
unchanged.

NO credence move (pre-stated in every branch). The in-pipeline
width-object program is now: leg 1 E-SURVIVE (8G), leg 2
ALL-FAIL-with-attribution (8H) — sq's identity remains open
in-pipeline; its census meter is confounded; the coherent-draw
reading stays alive by the K3=K1 equivalence (uncontested, not
supported).

Plain verdict: SUCCESS on the attribution question (decisive,
0.78–0.81, the flood localized to a named model component with a
named repair); NEEDS DIFFERENT DATA on the boundedness question
itself (the census cannot measure the width object's support at
the current companion-model grade).

ELI12: We asked our nine fastest star pairs to referee the jitter
mystery: "real physics respects the escape speed limit, noise
doesn't — which is the jitter?" The referee turned out to be
shouting about something else entirely. Our model of hidden THIRD
stars (little companions that make the light wobble) predicts
about fourteen pairs beyond the speed limit — the sky shows two.
That's not the jitter's fault: even with jitter switched off, the
companion model floods the forbidden zone all by itself. So the
verdict: the speed-limit test can't judge the jitter yet, because
the companion model is drowning out the courtroom — and that
model was ALREADY the prime suspect from two earlier
measurements. The silver lining is real: we now know exactly
which dial is broken (the companion model's loud tail), three
independent measurements point at it, and fixing it un-blocks
two other tests at once. The jitter keeps its secret a little
longer — but the list of places the secret can hide keeps
shrinking.

## Stage 8I-a PRE-REGISTRATION (2026-08-04, committed BEFORE any run): THE WOBBLE-SURVIVAL INSTRUMENT — the wobble-tail repair, measurement round

PREMISE (the three strains, all shipped): 7J-z4 kinematic
wobble-binding (Dwob +314/+306; forcing measured multiplicity
costs 135–153); round-10 q-moments (detected companions carry
0.13–0.17 of undetected wobble variance); 8H census counts (~13
phantom window pairs from the wobble channel at kw = the grid
floor). THE HYPOTHESIS (selection physics, not new dynamics):
EDR3 catalog quality cuts (RUWE-class; Penoyre–Belokurov-type
subsystem astrometric inflation) remove pairs with large
photocenter wobbles; the v7 model RETAINS those systems and
hands them clean PMs with km/s spikes. The repair: a survival
cap — systems whose summed active-companion wobble amplitude
kw·(w1+w2)·4.74047 km/s exceeds wcut LEAVE the model population
(as they leave the catalog).

DESIGN (mode WSRV=1 on the landed photow3 grid; e-sector PINNED
at the 8G-measured mode ein = 1.0, erf = 0.95; seeds 31 + 101,
both laws): WCUT_GRID = [0.05, 0.10, 0.20, 0.40, 1e9] km/s (the
1e9 node = survival OFF = the identity cell), one new axis
between kw and sq; cubes data/stage7j_cube_full_wsrv_{seed}_
{law}.npy, 10-dim; fresh names. Survival is applied to the
population BEFORE the catalog-cut/histogram chain (numerator and
normalization both — the pair is simply absent). The wcut = 1e9
branch passes the untouched arrays (no boolean indexing), so the
identity slice is bit-exact by construction. GB0w/GB0e
SKIPPED-disclosed under WSRV (G8I-0 substitutes).

GATES (abort-grade):
  G8I-0 (in-run, per law x seed): the wcut = 1e9 slice of the
    WSRV cube AND cubevt must equal the esec cube's (ein = 1.0,
    erf = 0.95) slice EXACTLY (max|d| <= 1e-9, target 0.0e+0) —
    the pinned-sector + survival-off model IS that slice.
  G8I-1 (reader): the identity-slice LANDED-CONV read computed
    from the wsrv cube must equal the same read computed from
    the esec slice directly (arithmetic identity, <= 1e-9).
  G8I-2 (reader): finite-mass bookkeeping + EDGE report on WCUT
    (mode at 0.05 with >= 0.5 mass -> EDGE flag, correction-#4).

THE READ (calcs/stage8i_read.py; LANDED-CONV anchor + prior_eta):
per law x seed — P(wcut) posterior; alpha_marg, dN (marginal over
everything incl. wcut); P(fcomp) (does the posterior move toward
the measured ~0.3 host rate once survival explains the kinematic
tolerance?); P(sq) (does the width demand deflate under the
repaired tail?); P(fpm); Dwob analogue at the extended grid; THE
CENSUS FORWARD at the repaired MAP cell (8H band_mu + the surv
cull at the cell's wcut/kw/fcomp; comparator = the same cell at
wcut = 1e9), scored at (band 9, cliff 2) by the 8H convention.

BARS (locked; seed means, both laws):
  W-DEMANDED: P(wcut < 1e9) >= 0.90 both laws -> the sky wants
    the survival cap; the wobble-tail repair is DATA-SUPPORTED.
  W-REFUSED:  P(wcut = 1e9) >= 0.50 both laws -> the sky refuses
    the repair; the strain stands unresolved; census stays
    confounded.
  else W-PARTIAL: decomposition reported as-is.
  CENSUS-REOPENED: repaired-MAP-cell jointP >= 1e-3 (the 8H bar,
    continuity) in >= 3 of 4 law x seed runs -> the census meter
    is CONDITIONALLY repaired (pending 8I-b power).
  Movement grammar (descriptive this round): alpha/fcomp/sq
    movements reported with the MATERIAL flag at |d alpha| >
    0.11; if alpha_marg <= 0.2 AND dN <= +5 both laws, print the
    pre-scripted sentence: "fifth-move-shaped movement observed —
    VERDICT DEFERRED to the powered 8I-b" and NOTHING stronger.
  CREDENCE FROZEN in 8I-a (pre-stated, every branch — the
    cadence rule): 8I-b (the powered round: own-truth injection
    arms through the WSRV fitter + the credence map) is the
    decider; its map is pre-registered there, not here.
  EXPECTED OUTCOME (pre-stated): W-DEMANDED with wcut interior
    (0.10–0.40) — the three strains all point at a finite cap;
    the census-reopen is the live uncertain branch (the band may
    stay smear-flooded at sq = 0.2 even with the tail culled).

Cost estimate (non-binding): ~60 min/seed (the WCUT axis is
CPU-side; GPU cells unchanged at 90/seed). Amendment rule: gate
tolerances amendable pre-quote with logged reasons; bars may NOT
move.

## Stage 8I-a EXECUTED (2026-08-04, all gates PASS): W-REFUSED 4/4 — the sky puts ZERO mass on the survival cap; the wobble mismatch is SHAPE-level, not survival-level

Gates: G8I-0 bit-exact 8/8 in-run (0.00e+00 incl. vt — the WSRV
loop restructure certified against the esec slice); G8I-1 4/4
(reader arithmetic identity). Runtime 55.6 + 56.2 min (estimate
borne out). EXPECTATION MISSED — the SECOND consecutive stage to
land against the pre-stated bet (8G expected absorb -> survived;
8I-a expected W-DEMANDED -> refused). Logged prominently: the
strain triangle made the cap look inevitable; the sky disagreed.

THE ANSWER: **P(wcut < inf) = 0.00 in ALL FOUR law x seed reads**
— not a mild preference: the posterior concentrates entirely at
survival-off at every severity offered (0.05–0.40 km/s; even the
0.40 node, which culls only the ~1–3% most extreme spike systems,
is refused). alpha/dN unchanged to the third digit (the marginal
never touches the new axis); fcomp holds at 0.1; sq holds at 0.2;
CENSUS-REOPENED 0/4 (the repaired MAP cell IS the off cell — the
8H demotion stands).

THE DIAGNOSIS (the informative part): the kinematic likelihood
WANTS the spike population — the retained wobblers carry real
fitting work in the mid-shoulder (the 7J-z6 fingerprint region
the smear also feeds), while the census counts reject the same
population's EXTREME tail (x7 cliff overshoot at sq=0, 8H). A
system-level cull removes both at once, so the likelihood vetoes
it. Together 8H + 8I-a localize the wobble channel's defect to
its DISTRIBUTION SHAPE: roughly-right mid-amplitude mass,
too-heavy extreme tail — the (P, q) corners that survive the
S-suppression with km/s photocenter velocities. MECHANISM NOTE
(profile-level, quotable as diagnostic): the Dwob analogue
computed WITH the escape axis available collapses from the
7J-z4 comparator +314/+306 to -6..-24 — forcing measured
multiplicity at kw=1.4 stops being catastrophic the moment the
fit may cull wobblers, even though the marginal never elects to
— the binding is confirmed wobble-amplitude-carried, and the
tension between the likelihood's mid-shoulder appetite and the
census's tail rejection is now the sharpest statement of the
companion-sector defect.

CONSEQUENCES: (1) the SURVIVAL reading of the strain triangle is
DEAD (a bluntly-culled catalog is not what the data describe);
(2) the named successor is the wobble TAIL-SHAPE instrument —
reshape the spike distribution (the S period-smearing factor
and/or the logP tail sector) so the mid-shoulder mass survives
while the extreme tail dies; design deferred (NOT pre-committed
tonight — the 8I-a lesson is that this sector needs a
distribution-level instrument, not another scalar dial); (3)
8I-b in its planned powered form is MOOT (zero movement =
nothing to power-validate; the powered round waits for an
instrument the sky accepts); (4) the census demotion (8H)
stands; external legs (T2/T3, P9/DR4) unchanged.

CREDENCE FROZEN (pre-stated in every branch of the pre-reg).

Plain verdict: SUCCESS as an instrument (a decisive, gate-clean
answer that killed our best repair hypothesis in two hours and
localized the defect to distribution shape); the repair itself:
NEEDS REFINEMENT (the tail-shape successor, design pending).

ELI12: We built the bouncer: "throw out every star pair that
wiggles harder than X" — and let the sky set X. The sky's
answer: NO BOUNCER, at any strictness. Surprise — but it
teaches us something sharp. The fit secretly NEEDS the wigglers:
they help explain the crowded middle of the speed histogram.
What the sky objects to is only their loudest few — the ones
screaming past the speed limit where almost nobody real is
found. So the companion model isn't too big — it's the wrong
SHAPE: right amount of medium wiggling, too much extreme
wiggling. You can't fix a shape with a bouncer; you need a
tailor. Designing the tailor is the next job — done carefully,
not tonight. Two bets in a row have now landed against us,
which is exactly why we write the bets down first.

## Stage 8J PRE-REGISTRATION (2026-08-04, committed BEFORE any run): THE WOBBLE SATURATION INSTRUMENT — the distribution-level tail-shape repair

PREMISE (set by 8H + 8I-a): the wobble channel's defect is
DISTRIBUTION-SHAPE-level — the likelihood needs the spike
population's mid-amplitude mass (survival cull refused 4/4 at
zero posterior mass), the census rejects its extreme tail (×7
cliff overshoot at sq = 0). The repair must keep the mid and
trim the tail WITHOUT removing systems.

THE LAW (one parameter, monotone, physical): per-companion
photocenter-wobble SATURATION —
  w_eff = (w0/4.74047) · tanh(w · 4.74047/w0),   w0 in km/s.
Small wobbles pass to first order (the mid-shoulder survives);
large wobbles plateau at w0 (the km/s spike corner compresses).
Physics: a large photocenter orbit partially ABSORBS into the
5-parameter astrometric solution and inflates the errors rather
than leaking fully into PM — a saturation, not a disappearance
(and not a cull: the pair stays in the catalog, consistent with
the 8I-a refusal). The saturation acts on the wobble VELOCITY
only; the hidden dynamical mass mh (and its boost) is untouched
(the amendment-4 principle: mass is real regardless of light).
Applied per companion slot, before the kw nuisance scale.

DESIGN (mode WSAT=1; the landed photow3 grid; e-sector pinned at
the 8G mode ein = 1.0, erf = 0.95; seeds 31 + 101, both laws):
W0SAT_GRID = [0.1, 0.2, 0.4, 0.8, 1e9] km/s sharing the wx axis
slot (WSRV and WSAT are exclusive modes); the 1e9 node is
BRANCHED to the verbatim legacy accumulation arrays (no
transform arithmetic), so the identity slice is bit-exact by
construction. Cubes data/stage7j_cube_full_wsat_{seed}_{law}.npy,
10-dim; fresh names; GB0w/GB0e SKIPPED-disclosed (G8J-0
substitutes).

GATES (abort-grade): G8J-0 (in-run): the w0 = 1e9 slice AND its
vt cube must equal the esec (1.0, 0.95) slice EXACTLY (<= 1e-9,
target 0.0e+0). G8J-1 (reader): the same identity at read level
(arithmetic). G8J-2: EDGE report on w0 (mode at 0.1 with >= 0.5
mass -> EDGE flag, correction-#4).

THE READ (calcs/stage8j_read.py, the 8I-a reader chassis): P(w0)
posterior; alpha_marg/dN; P(fcomp) (reconciliation watch), P(sq),
P(fpm); the Dwob analogue; THE CENSUS FORWARD at the repaired MAP
cell with the SATURATED wobble in band_mu (mh untouched),
comparator at the identity node, scored at (band 9, cliff 2) by
the 8H convention.

BARS (locked; seed means, both laws):
  T-DEMANDED: P(w0 < 1e9) >= 0.90 both laws -> the sky accepts
    the saturation; the tail-shape repair is DATA-SUPPORTED; the
    POWERED ROUND (own-truth injections through the WSAT fitter +
    the pre-registered credence map) REVIVES as the next decider.
  T-REFUSED:  P(w0 = 1e9) >= 0.50 both laws -> the third
    distribution-level refusal; the wobble sector goes to the
    REVIEWER ROUND with the 8H/8I-a/8J scorecard — NOT to a
    fourth dial (pre-committed stop rule).
  else T-PARTIAL: decomposition reported as-is.
  CENSUS-REOPENED: repaired-MAP-cell jointP >= 1e-3 in >= 3 of 4
    law x seed runs (the 8H bar) -> the census meter conditionally
    repaired pending the powered round.
  Movement grammar: alpha/fcomp/sq movements descriptive with the
    MATERIAL flag at |d alpha| > 0.11; the fifth-move-shaped
    sentence (alpha_marg <= 0.2 AND dN <= +5 both laws) is
    pre-scripted DEFER-only.
  CREDENCE FROZEN in every branch (pre-stated; only the revived
    powered round may move it, under its own pre-registered map).
  EXPECTED OUTCOME (pre-stated, with explicit LOW CONFIDENCE
    after two consecutive missed bets): T-DEMANDED — the design
    targets exactly the likelihood's measured preference (keep
    mid, trim tail); the refusal branch is fully live and would
    be the third distribution-level fact about this sector.

Cost: ~55 min/seed GPU (the wx axis is CPU-side); setup and read
are minutes (LLM time — calibration noted). Amendment rule: gate
tolerances amendable pre-quote with logged reasons; bars may NOT
move.

## Stage 8J EXECUTED (2026-08-04, all gates PASS): T-REFUSED BY THE LETTER — wrapping a VIOLENT SEED SPLIT; seed 31 produces the fifth move LIVE in the marginal; the stop rule fires

Gates: G8J-0 bit-exact 8/8 (0.00e+00 incl. vt), G8J-1 4/4 — the
saturation machinery is certified; nothing below is wiring.
Runtime 64.9 + 64.9 min. EXPECTATION MISSED — the THIRD
consecutive stage against the pre-stated bet.

THE LETTER: seed-mean P(w0 = off) = 0.52 (simple) / 0.500 (BE —
EXACTLY at the bar, boundary-riding disclosed: one realization
percent flips the letter, not the content) -> T-REFUSED both
laws; d_alpha = -0.334/-0.314 MATERIAL; CENSUS-REOPENED 0/4.
By the pre-committed STOP RULE the wobble sector goes to the
REVIEWER ROUND with the 8H/8I-a/8J scorecard — no fourth dial.

THE CONTENT (the split, per-seed, both laws concordant within
each seed):
  SEED 31 — THE SATURATED WORLD: P(w0) puts 0.97–1.00 on strong
    saturation (0.1–0.2 km/s); fcomp mass moves 0.1 -> 0.35/0.50
    (0.97 of posterior at the MEASURED host-rate territory);
    fpm -> 3.0 (0.97/0.99 — the noise chase RETURNS); Dwob -> +0.0
    (the wobble binding GONE, mechanism-consistent); sq holds 0.2;
    **alpha_marg = 0.00, dN = +0.0 — NEWTON TIES, both laws: the
    fifth-move configuration, voluntarily, in the marginal.** The
    7J-z3 S3 collapse (forced fcomp >= 0.35 at cost 135–153) now
    happens for free once the spikes saturate.
  SEED 101 — THE REFUSING WORLD: P(w0 = off) = 1.00 flat; fcomp
    stays 0.1; alpha_marg 0.44/0.50, dN +8.9/+11.3 — the 8I-a
    refusal repeated exactly.
The DATA are identical between these rows; only the model-side
Monte Carlo realization differs. The 3A realization systematic
(correction-#10 lineage) is now DECIDING a verdict-grade
question: the likelihood gap between the collapse world and the
boost world is smaller than the realization scatter at two
seeds. Nothing about this can be resolved without a seed
budget — which belongs to the review round per the stop rule,
not to a unilateral extension tonight.

THE CENSUS PUSHBACK (the one model-light leg standing against
the collapse world): 0/4 reopened; the seed-31 collapse cell
(fcomp = 0.35, w0 = 0.2, fpm = 3.0) still over-produces the band
x3 (mu = 29.27 vs 9; jointP 2.6e-9), and its saturation-OFF
comparator is catastrophic (67.49, 68.54) — the many-companion
world, saturated or not, cannot reproduce the (9, 2) pair that
the sky actually shows. The fifth-move world has an unanswered
census bill.

FIFTH-MOVE GRAMMAR: the pre-scripted trigger (both laws AND both
seeds) did not fire; the observation is booked as
fifth-move-shaped movement AT SEED 31 ONLY, deferred by
construction to a powered + seed-budgeted round.

CREDENCE FROZEN (pre-stated in every branch). The stop-rule
review round is where the map gets pre-registered; the inputs it
must weigh: (a) the seed-31 collapse world and its census bill,
(b) the seed-101 refusal, (c) the realization systematic as the
current decider-blocker, (d) the three-round distribution-shape
scorecard (8H attribution, 8I-a refusal, 8J split).

Plain verdict: SUCCESS as an instrument at the highest grade —
it exposed the sector's true state: the binary anomaly's
fifth-move alternative (saturation + measured multiplicity +
max noise) is ONE REALIZATION away from a Newton tie, two
realizations disagree violently, and the census still vetoes
the alternative. The question 8J asked (is the tailor demanded?)
returns UNRESOLVED-SPLIT; the pre-committed stop rule routes
everything to external review with a seed budget.

ELI12: We gave the fit the tailor. One copy of our simulated
universe LOVED it — tailored the wigglers, hired lots of
companions at exactly the rate the sky measures, cranked the
noise dial, and said "I don't need extra gravity at all" — a
perfect tie with plain Newton. The OTHER copy said "no tailor,
no extra companions, keep the gravity boost." Same real sky,
two different random casts of fake stars, opposite verdicts —
which means our random-cast wobble is currently louder than the
answer, and the honest fix is more casts, supervised by an
outside reviewer, exactly as our stop rule pre-ordered. One
referee still shouts from the corner: the nine-fastest-pairs
census says the companion-heavy story predicts three times too
many fast pairs — the boring explanation still owes the census
an answer even in its best world. Nobody's number moved
tonight; the process held. Third lost bet in a row — and this
is precisely the kind of night the betting ledger was built for.

## REVIEWER ROUND 4 RECEIVED + ADOPTION DECISIONS (2026-08-04): process 9/10; the collapse world graded an absorber-stack at the unphysical corner; the seed-budget resequenced behind a PHYSICAL astrometric-response model; the credence map PRE-REGISTERED

THE ROUND (fresh-Opus, context intact from rounds 1–3; full text
in the session record): (A) process 9/10, end state "honest
UNRESOLVED, correctly escalated"; the seed-31 collapse world =
an alternative-absorber stack pinned at the unphysical corner
(fpm 3.0 past the ~1.4 ceiling; fcomp 0.35–0.50 at 3–5× the
kinematic pin, reconciled only by the new saturation DOF; the
new-DOF-opens-an-expensive-corner shape flagged as
textbook-overfit-like) — "physically disfavored, not
likelihood-excluded, credibility gated on the saturation's
physicality." THE TIE-BREAKER: not more seeds — a DIRECT
inner-companion census of the boost-carrying and census pairs
(Gaia DR3 non-single-star solutions, RUWE/excess-noise
signatures, SB1 RV variability): the collapse world REQUIRES
~35–50% active-but-quiet companions in those systems, at the
object level. (B) seed budget: 12 interim / 24 cap, per-seed
CLASSIFICATION HISTOGRAM primary (pooled marginal-of-marginals
forbidden as headline — the bimodal-mean trap, our own round-10
rule), seed-mean threshold letters RETIRED, both-truth DEGRADED
injections (7J-w2 standard) through the same fitter, split = a
defined outcome. (C) the credence map (below, adopted). (D)
attacks: D1 the 8J letter mislabeled (UNRESOLVED-BIMODAL, not
T-REFUSED); D2 replace the tanh with the physical photocenter
law; D3 saturation freedom erodes the boost even where it
survives (seed-101 dN +8.9/+11.3 vs the operative +14.5–23.8);
D4 the 8I-a summed-amplitude cut ignores the period-vs-baseline
physics of RUWE (short-P averages out); D5 the 8H truncation
kernel is itself a choice — test a physical orbital pile-up.

OUR AUDIT (adopt designs, correct premises — the standing rule):
- D2 PREMISE SLIP, corrected: the model's wobble AMPLITUDE law
  ALREADY IS the physical photocenter law |q/(1+q) − l/(1+l)|
  (amendment 4, Part-A-validated, twin-quiet). What is
  phenomenological is the tanh SATURATION stacked on top. The
  reviewer's demand survives the correction in merged D2+D4
  form: replace BOTH proxies (survival cut, tanh) with the
  PERIOD-RESOLVED Gaia astrometric-response model — leakage vs
  absorption as a function of period against the baseline
  (short-P averages out, P~baseline partially absorbs,
  long-P leaks) — built from the published response curves
  (Penoyre–Belokurov class; scout queued). ADOPTED as the
  gating item before any seed budget.
- D1 ADOPTED AS A RELABEL: the locked 8J letter stands as
  printed (bars are never retro-moved); the OPERATIVE
  characterization becomes UNRESOLVED-BIMODAL; seed-mean-over-
  bimodal threshold bars are retired from future instruments.
  Ledger note appended.
- D3 ADOPTED as a conditional annotation: IF the response model
  proves a physical saturation channel, the honest Newton
  rejection under that freedom is ~+9–11 (seed-101 grade), down
  from the operative +14.5–23.8 — carried openly, gated on the
  response-model round.
- A/B ADOPTED with resequencing: (1) the DR3 companion census
  (the tie-breaker; converges with queued T2; cheap, external,
  object-level) + the response-curve scout; (2) the physical-
  response fitter replacing tanh/survival; (3) the 24-seed
  budget with both-truth degraded injections under the adopted
  map. D5 folds into (2).
- C ADOPTED WITH AMENDMENTS, and hereby PRE-REGISTERED as the
  map for the budgeted round (conditioning re-scoped per the D2
  correction: "with the period-resolved response model"):
    floor ~25% (the model-light census — 9 forbidden-band pairs
      no error model produces — survives every branch);
    BOOST WINS -> ~60%;
    COLLAPSE WINS -> ~30% (held above floor by the census bill
      falling mechanism-intrinsically on the collapse world; IF
      the enabling saturation is shown unphysical, this branch
      pre-emptively becomes ~40%);
    SPLIT PERSISTS -> ~48%, and thereafter credence is gated on
      EXTERNAL data only (T2/DR3-census/DR4) — no in-pipeline
      stage moves it.
  The structural-asymmetry argument (collapse-flood mechanism-
  intrinsic vs boost-flood parameterization-contingent) is
  priced at reviewer's +2–3% on the split branch only — "a
  hypothesis, not evidence" until the width object's correct
  parameterization is exhibited. Accepted.
- Calibration note: the round contains the program's 5th
  reviewer premise slip (D2's "you used tanh instead of the law
  you validated" — we use both, at different layers), alongside
  genuinely decisive contributions (the tie-breaker, the
  histogram-primary budget design, the map). The
  adopt-designs-correct-premises pattern holds.

CREDENCE: STILL FROZEN at ~50% (the map above governs the
FUTURE budgeted round; nothing moved today).

Plain verdict: the review round delivered — the sector's next
three moves are now externally designed, the fifth-move world
has a named object-level kill test, and the credence map is
locked before the instrument that will use it exists. That is
the discipline working exactly as built.

ELI12: The referee looked at our coin-flip and said three
things. One: your process is excellent — nine out of ten. Two:
the "no gravity needed" world is leaning on four cranked-up
dials at once, including a noise dial past its certified
ceiling, so treat it with suspicion — but the way to kill it
isn't more dice: go LOOK for the hidden companions it needs.
If half our star pairs secretly host them, Gaia's own wobble
flags and radial-velocity surveys will show them; if they
don't, that world dies at the object level. Three: before any
big dice campaign, replace our hand-drawn "how companions hide"
curve with the real physics curve — because the whole question
is whether real companions CAN hide. We also caught the referee
mis-remembering one thing (we already use the real wobble-size
law; it's the hiding-with-time piece that needs building) —
trust but verify runs in both directions here. The belief
number stayed frozen at 50-50; the rulebook for moving it is
now written and signed before the game.

## Stage 8K PRE-REGISTRATION (2026-08-04, committed BEFORE any run): THE WOBBLE CENSUS — the reviewer's object-level tie-breaker, run on columns already in the catalog

PREMISE: the collapse world (8J seed 31, review-graded an
absorber-stack) REQUIRES fcomp = 0.35–0.50 per component of
astrometrically ACTIVE-but-saturated companions among the pairs
that carry the boost. The light-channel census (7J-z2c: host
0.23 twin-heavy / 0.29–0.32 flat-q per component) cannot address
ACTIVITY (light counts all periods; the collapse needs the
P ~ 0.3–10 yr wobble window). The saturation story's own
signature: absorbed wobble = inflated residuals — elevated
RUWE / astrometric_excess_noise below the survival ceiling. The
catalog carries ruwe1/2, aen1/2 + significance, ipd flags, DR2
RVs + errors, and DR2-epoch RUWE — the census is fully local.

DESIGN (calcs/stage8k_wobblecensus.py, CPU-only; the 8F-b/4J
loader verbatim; differential where possible):
  Populations: WIDE = s >= 6 kAU (the boost-carrying bins);
  NARROW control = 0.2–2 kAU (same selection, boost-negligible);
  the 9 census pairs (ceiling_pairs.csv census_corr rows, matched
  back to catalog rows on (s, Mtot, vc) to float tolerance).
  HOT (per component, pre-registered constants): ruwe > 1.25 OR
  (aen_sig > 3). Pair-hot = either component hot.
  S1 the RATE: f_hot per component in WIDE, binomial CI.
    VALIDITY CONDITION (pre-stated): S1's absolute bar is valid
    only if the catalog shows NO hard RUWE ceiling — P99(ruwe) >
    1.6. If a ceiling is detected, S1 -> DISCLOSED-INVALID and
    the verdict rests on S3 + S2 (fallback pre-stated).
  S2 the CORRELATION: rank correlation of pair activity vs vt in
    WIDE (hot indicator and mag-binned RUWE percentile variants),
    permutation null 1e4. ASYMMETRIC interpretation pre-stated:
    positive = NECESSARY for collapse but degenerate with the
    error channel (not proof); null = collapse-killing for the
    vt-excess-as-wobble reading.
  S3 the NINE: hot-pair count among the 9 census pairs (their
    band membership IS wobble under the collapse+8H attribution).
  S4/S5 descriptive, no bars: DR2 RV-error inflation hot-vs-cold;
    DR2->EDR3 RUWE deltas (acceleration signature); census pairs
    individually printed.
GATES: G8K-0 loader identity (14,071 pairs; per-bin N = [9950,
  2684, 1223, 214]; 9/9 census rows matched). G8K-1 flag
  completeness (NaN fraction < 5% or disclosed). G8K-2
  mag-dependence report (if f_hot runs strongly with G, the
  mag-binned variant is primary — pre-stated).

BARS (locked):
  S1: f_hot(WIDE) < 0.20 -> RATE-SHORTFALL; >= 0.30 -> RATE-MET;
    else GRAY. (The collapse requirement is >= ~0.30 after
    activity-window inefficiency; the light census itself
    predicts ~0.10–0.15 hot under the boost world.)
  S3: <= 3/9 hot -> CENSUS-CLEAN (baseline-consistent; the
    collapse account of the band dies at object level and the
    4J/8F-b defense upgrades: the band pairs are neither error
    tails NOR wobblers); >= 6/9 -> CENSUS-WOBBLE-SUSPECT (a real
    hit AGAINST the census's gravity reading — priced below,
    honesty both directions); else MIXED.
  VERDICT: OBJECT-LEVEL-ABSENT if RATE-SHORTFALL AND
    CENSUS-CLEAN (S1 valid) or CENSUS-CLEAN + S2-null (S1
    invalid). WOBBLE-RICH if RATE-MET AND >= 6/9. Else
    MIXED-REPORTED.
CREDENCE MAP FOR THIS STAGE (pre-stated; the first authorized
movement since the freeze, both directions priced):
  OBJECT-LEVEL-ABSENT -> anomaly-real ~50% -> ~55%, AND the
    24-seed budget is DE-PRIORITIZED (the collapse basin loses
    its population at object level; reviewer to be notified);
  WOBBLE-RICH -> ~50% -> ~40%, AND the census's gravity-evidence
    reading is SUSPENDED pending the NSS leg; budget proceeds
    under the registered map;
  MIXED -> ~50% HELD; budget proceeds.
FOLLOW-UP (named, not run tonight): 8K-b = the Gaia DR3
non-single-star TAP crossmatch (orbital + acceleration
solutions) — the network leg.
EXPECTED OUTCOME (pre-stated, low confidence — three straight
misses): genuinely uncertain; mild lean OBJECT-LEVEL-ABSENT
(the light census already sits below the collapse requirement).
Amendment rule: gate tolerances amendable pre-quote; bars and
the map may NOT move.

AMENDMENT 1 (2026-08-04, logged PRE-QUOTE; the first-run verdict
line is a FAILED-VALIDITY record, preserved as
data/stage8k_wobblecensus_run1.txt, never quoted as a result):
the instrument's own G8K-2 gate measured the HOT flag
mis-calibrated — f_hot = 0.96–0.99 for ALL components brighter
than G ~ 12 (mag run 5.32, STRONG). Diagnosis: the aen_sig > 3
arm is a BRIGHTNESS flag, not a binarity flag — bright Gaia
sources have tiny formal errors, so calibration-level excess
noise is always "significant" (known DR3 property). The
first-run baselines it produced (f_hot 0.57–0.63 per component;
pair-hot 0.79) are physically impossible as companion rates,
and the S1/S3 bar constants were written for a binarity-grade
flag with a ~0.10–0.20 baseline. THE AMENDMENT (the 8H
amendment-1 precedent — a gate fired on instrument internals;
the component is replaced; bars and the map do NOT move):
HOT = ruwe > 1.25 ONLY (the literature-standard Belokurov-class
threshold; RUWE is brightness-normalized by construction —
verified in the amended run's own mag table). The aen_sig
channel is retained DESCRIPTIVELY (mag-matched percentile
comparison for the census pairs — not silently dropped). The
S2 mag-binned-percentile variant (already pre-declared primary
under the fired trigger) is implemented. S4 sentinel-value bug
fixed (RV errors filtered to finite < 1e3). Design-flaw
ownership: the pre-reg promoted the mag-robust variant for S2
only and left S1/S3 on the raw flag — that hole is ours, caught
by our own gate; both runs preserved.

## Stage 8K EXECUTED (2026-08-04, amended run, all gates PASS): OBJECT-LEVEL-ABSENT — the collapse world's companion population is not in the catalog; the census pairs are astrometrically clean; THE MAP EXECUTES: anomaly-real ~50% -> ~55%

Gates: G8K-0 loader identity (14,071; bins exact; census pairs
matched 9/9 uniquely), G8K-1 completeness 0.0000, S1 validity
VALID (P99 RUWE = 2.40, no ceiling). Amendment 1 (logged
pre-quote, above): run-1's aen_sig arm was a brightness flag
(0.96–0.99 hot at G < 12, caught by G8K-2's own mag table);
HOT = ruwe > 1.25 only; run-1 verdict line = the failed-validity
record, preserved, never quoted. The amended mag table is
binarity-sane (0.06–0.22, peaking in RUWE's known G 11–13
sensitivity window). EXPECTATION MET (the mild ABSENT lean) —
after three straight misses.

THE NUMBERS:
  S1 RATE-SHORTFALL: f_hot(WIDE, per component) = 0.090
    [0.080, 0.101] vs the collapse requirement >= ~0.30 — a 3x
    shortfall at tight CI. Even RUWE's most sensitive mag window
    peaks at 0.22. Bonus fact: WIDE (0.090) is CLEANER than the
    NARROW control (0.142) — the boost-carrying bins are the
    astrometrically quietest part of the catalog, the opposite
    of any wide-excess-as-astrometric-contamination story.
  S3 CENSUS-CLEAN: 2/9 census pairs hot vs the 0.166 pair-hot
    baseline (expected ~1.5) — exactly baseline; 16 of 18
    components have pristine RUWE (0.89–1.24). THE UPGRADE: the
    nine band pairs are now defended on THREE independent
    channels — not error tails (4J 3.8e-9, 8F-b 8.5e-8), not
    wobblers (8K), while the light channel (v2c) already sat
    below the collapse rate.
  S2 (asymmetric reading as pre-stated): a REAL positive
    correlation exists (mean vt hot-cold = +0.17, perm p < 1e-4;
    mag-binned rho = +0.12, p = 1e-5) — necessary-for-collapse
    but error-degenerate, and with S1 failing 3x it reads as the
    fcomp ~ 0.1 real companion sector doing its known
    mid-shoulder work, not as a 0.35–0.50 population. Carried,
    not interpreted further. S5 (desc): hot pairs show a
    POSITIVE DR2->EDR3 RUWE delta (+0.098 vs -0.030 cold) —
    epoch-evolving wobble, i.e., the flag is seeing real
    astrophysics where it fires.
  WINDOW-OVERLAP CAVEAT (carried): very-short-P companions hide
    from RUWE (average out) — but they also carry no wobble in
    the model (S = P/17.8), so they cannot power the collapse
    mechanism either; the collapse-relevant and RUWE-sensitive
    period windows overlap by construction. The 8K-b NSS TAP
    leg (orbital/acceleration solutions) remains named for the
    residual window.

VERDICT (locked bars): OBJECT-LEVEL-ABSENT — S1 shortfall AND
census-clean, S1 valid. THE PRE-REGISTERED MAP EXECUTES (the
first credence movement since the freeze, both directions were
priced): **binary anomaly-real ~50% -> ~55%**; the 24-seed
budget is DE-PRIORITIZED (the collapse basin loses its required
population at the object level — the fitter corner exists, its
stars do not); the reviewer is being notified with the two-run
story and the amendment open to attack. Galaxy-side credence
untouched (~65% provisional). The census's gravity-evidence
reading recovers only its DEFENSIVE role (three-channel-clean);
the forward-model demotion (8H) stands until the response-model
round.

Plain verdict: SUCCESS — the object-level kill test the
reviewer designed, run same-hour on columns we already owned,
with one self-caught flag amendment; the collapse world now
needs a companion population that three independent channels
say is not there.

ELI12: The referee said "stop arguing about dice — go count the
hidden companions the no-gravity story needs." We opened Gaia's
own quality flags and counted. First try, our companion
detector was accidentally wired to flag BRIGHT stars (96% of
them lit up — nonsense); our own sanity table caught it, we
fixed the wiring to the standard detector, kept both printouts.
The honest count: the no-gravity story needs almost a third of
these stars wobbling; the real number is nine percent. And our
nine fastest pairs — the crown jewels — are among the CLEANEST
stars in the whole catalog: not noisy, not wobbling, just fast.
So the "it's all hidden companions" world fails its own
head-count, at the level of actual objects, and for the first
time since we froze it, the belief number moves — from 50 to
55, exactly as the pre-signed rulebook said it must. Small
step, honestly earned, three lost bets and then a won one.

## REVIEWER ROUND 5 (targeted, 2026-08-04): the 8K amendment and S2 handling ENDORSED; one precise reservation — the collapse world's SELF-PREDICTED f_hot is the closing number; credence HOLDS at 55

THE RESPONSE: (1) the amendment SOUND — aen_sig's bright-star
pathology is documented (it is why RUWE became the standard),
the catching gate was pre-committed, the invalidation is
outcome-independent; two points in our favor we had not
pressed: ruwe > 1.25 is the MORE companion-inclusive threshold
(so 0.090 is conservative-high), and WIDE-cleaner-than-NARROW
is the wrong sign for a companion-driven anomaly. Fair nit
accepted: the bright-star failure was foreseeable — an ideal
pre-reg would have mag-gated aen_sig a priori. (2) S2 handling
DEFENSIBLE and stronger than we stated: the correlation is
measured among the same systems S1 counted, so it cannot imply
a larger RUWE-visible population; asked-to-show: the forward
number that +0.17 is what fcomp ~ 0.1 predicts through the
RUWE-hot selection. (3) THE RESERVATION (could rise to
correction grade): "windows overlap by construction" is airtight
only if the collapse world's velocity-faking wobble is nonzero
exclusively where RUWE is sensitive. THE DECISIVE CHECK: push
the seed-31 collapse world (fcomp 0.35–0.50, saturated) through
a RUWE forward model — does it SELF-PREDICT f_hot >= ~0.30?
Self-predicts >= 0.30 -> 0.090 falsifies it cleanly (reviewer
would support 55–58). Self-predicts ~ 0.10 -> the census is
blind exactly at the collapse's hiding place; the S1 leg gains
the qualifier "absent IF RUWE-visible" at correction grade.
S3 he grades IMMUNE and our strongest leg regardless (the nine
pairs are individually pristine — those specific velocities are
not inner-companion wobble, modulo the NSS sliver).

ADOPTIONS + ONE CORRECTION OF HIS: the S1 verdict carries the
conditional annotation AS OF NOW — "collapse companions absent
IF RUWE-visible; the self-consistency number closes it" — and
the check is ASSIGNED to the response-model round (D2+D4),
which now has its deliverable list: (a) the period-resolved
leakage/absorption curves R(P/T) (Penoyre–Belokurov-class;
scout queued); (b) THE self-consistency number f_hot(collapse
world) — the decisive check; (c) the S2 forward number at
fcomp = 0.1; (d) the honest Newton band under physical
saturation freedom (his D3); (e) the census-meter re-opening
test. HIS long-P amplitude claim corrected on our side (6th
premise slip, minor, logged): long-P wobble velocity is NOT
small — v_orb falls only as P^(-1/3) (a 30-yr companion at
~8.8 AU carries ~9 km/s orbital, wfac-scaled to ~0.5–2 km/s)
— the long-P window is velocity-faking AND RUWE-quiet, gated
in practice by the catalog velocity cut; exactly why the
forward number, not hand arguments, closes this. The substance
of his reservation is thereby STRENGTHENED, not weakened, by
the correction.

CREDENCE: HOLDS at ~55% (his recommendation and the cadence
rule agree): the 57–58 upgrade or the correction-grade
qualifier both wait on the self-consistency number. Budget
stays de-prioritized; 8K-b NSS leg queued for the long-P/NSS
sliver.

ELI12: The referee checked our star-count and signed off on
the fixed detector and the small belief bump — with one sharp
homework problem: "prove the hiding companions your rival
theory imagines would have shown up on the detector you used.
If they'd be invisible to it, your count didn't test them."
Fair. The physics needed to answer is the same curve we
already planned to build next. He also called our nine fast
pairs the strongest card in the deck — individually spotless,
no matter how the homework turns out. And we caught him in one
small slip (slow-orbit companions are NOT gentle — they carry
real speed), which actually makes his homework MORE important,
not less. Belief stays at 55 until the curve says otherwise.

## Stage 8L-a PRE-REGISTRATION (2026-08-04, committed BEFORE any run): THE RESPONSE MODEL, decisive half — the leakage/residual curves + the reviewer's self-consistency number

DELIVERABLES (a)–(c) of the response-model round (the fitter
half (d)–(e) = 8L-b, separate pre-reg): (a) the period-resolved
curves for a photocenter orbit against a 5-parameter fit over
baseline T — the PM LEAKAGE L(P/T) (fraction of instantaneous
photocenter velocity biasing the fitted PM) and the RESIDUAL
R(P/T) (RMS unfitted deviation as a fraction of the angular
amplitude, driving RUWE via the Penoyre–Belokurov-class
inflation sqrt(1 + (R a_phot / sigma_AL)^2)); derived
analytically (circular orbit, uniform sampling, phase/
orientation-averaged; eccentricity refinement deferred) with
exact limits as gates. (b) THE SELF-CONSISTENCY NUMBER: the
seed-31 collapse world (fcomp 0.35/0.50, its companion draws,
hosts sampled from the catalog's (plx, G) joint) pushed through
the RUWE forward -> predicted f_hot vs the measured 0.090.
(c) the S2 forward number: the boost world (fcomp = 0.1) ->
predicted Delta-vt(hot - cold) vs the measured +0.17.
Constants (sigma_AL(G), baseline T, N_obs, the RUWE u0 role)
are SCOUT-SOURCED (Haiku, primary: Penoyre & Belokurov 2020,
Lindegren+2021, the RUWE technical note) and inserted at
implementation with citations; the scout also hunts the
provenance of the legacy S = min(1, P/17.8 yr) factor.

GATES:
  GL0 (exact limits, sympy/quadrature): L -> 1 as P/T -> inf,
    L -> 0 as P/T -> 0; R -> 0 as P/T -> inf (curvature-
    suppressed), R -> orbit-RMS as P/T -> 0; numeric quadrature
    vs analytic <= 1e-6; continuity.
  GL1 (REPORTED deliverable, interpretation DEFERRED): the
    derived L(P/T) vs the operative model's legacy S = min(1,
    P/17.8 yr) — any material disagreement is a candidate
    model-side correction that goes through its OWN
    pre-registered stage (8L-b), not this one.
  GL2 (VALIDITY, the 8K pattern): the forward at the OPERATIVE
    world (fcomp = 0.1) must postdict the measured sky rate —
    f_hot(WIDE) prediction in [0.05, 0.15] (measured 0.090).
    FAIL -> the absolute scale is broken; the (b) bars are
    DISCLOSED-INVALID and the round returns instrument-repair,
    no verdict, no credence movement.

BARS + THE MAP for (b) (locked; executes only if GL2 passes;
this implements the reviewer's round-5 branch and the 8K
conditional):
  CLOSES-CLEAN: f_hot(collapse) >= 0.30 -> the collapse world's
    own companions would have tripped the 8K census; 0.090
    falsifies them cleanly; the S1 conditional annotation is
    DISCHARGED; anomaly-real ~55% -> ~57%.
  BLIND-SPOT: f_hot(collapse) < 0.15 -> the census was blind at
    the collapse's hiding place; the S1 leg rises to CORRECTION
    grade (the 8K map's ABSENT premise partially unwinds);
    anomaly-real ~55% -> ~50% (the 8K move retreats; S3 and the
    three-channel-clean census pairs stand on their own);
    the NSS leg (8K-b) becomes the decider.
  GRAY: 0.15-0.30 -> hold ~55%; the NSS leg decides.
BAR for (c) (descriptive consistency, no credence coupling):
  predicted Delta-vt within [0.5x, 2x] of +0.17 -> CONSISTENT
  (the round-5 "shown not asserted" item closed); outside ->
  reported as a forward-model residual, carried.
EXPECTED OUTCOME (pre-stated, low confidence): genuinely open;
  the physics cuts both ways (mid-P companions are RUWE-loud,
  but the collapse world's saturation-selected corner may
  concentrate mass at RUWE-quiet periods). No lean stated.
Amendment rule: gate tolerances amendable pre-quote with logged
reasons; bars and the map may NOT move.

## Stage 8L-a RUN 1: GL2 ABORT (the gate worked; no verdict, no movement) + AMENDMENT 1 + THE 8L-a2 PRE-REGISTRATION (2026-08-04, committed before the repaired run)

RUN 1 (preserved as data/stage8la_response_run1.txt): GL0 PASSED
after one gate-internal numeric-reference fix (the first GL0
firing caught the endpoint-inclusive grid's O(2/N) second-moment
offset; midpoint sampling -> 1.58e-09). THE CURVES ARE GOOD and
scout-shape-confirmed (Belokurov+20: excess peaks at P ~ T,
P^(2/3) below, strong suppression above). GL1 REPORTED (stands
regardless): the derived leakage L vs the legacy S = min(1,
P/17.8 yr) disagree strongly at mid-P (P = 6 yr: L = 0.80 vs
S = 0.34; P = 3 yr: L = 0.36 vs 0.17; and L is already 0.98 at
17.8 yr) — AND the scout finds the 17.8-yr constant NOT IN THE
LITERATURE (provenance likely internal). Interpretation stays
deferred to 8L-b, now with sharpened stakes. GL2 FAILED as
designed: the operative-world forward predicts f_hot = 0.019 vs
measured 0.090 — DIAGNOSIS: the absolute-rate observable was
mis-conceived; the sky's 0.090 is dominated by the SINGLES'
intrinsic RUWE spread (P90 = 1.294 puts ~10% of all components
near threshold from calibration scatter), which the base=unity
convention deliberately excludes. The S2 light forward
(+2.84 vs +0.174) exposed a second omission: no catalog velocity
cut (the sky removes the monster leaks the forward kept). Per
the pre-reg: bars DISCLOSED-INVALID, instrument-repair, NO
credence movement; the 8L-a map is VOID WITH ITS INSTRUMENT
(nothing ever moved on it).

AMENDMENT 1 (pre-quote): (i) numerical — R2 series branch
widened (u < 0.05) + core clipped at 0 (catastrophic-
cancellation NaNs at P > ~9000 yr, harmless but sloppy); (ii)
the S2 forward gains the catalog velocity ceiling (per-pair
2.978/sqrt(s) + 2.8284 sigma survival); (iii) the validity and
verdict OBSERVABLES are re-scoped to the calibration-immune
conditional form below.

## Stage 8L-a2 PRE-REGISTRATION (same commit): THE FAKER-CONDITIONAL — the reviewer's question in calibration-immune form

THE INSTRUMENT: for each of the 9 census pairs, over the model's
OWN companion law (q flat, logP N(5.03, 2.28), valid mask,
photocenter wfac), restricted to draws that could actually FAKE
the pair — projected leaked speed L(u)·wfac·v_orb·proj >=
0.5·vt_j·vc_j (the wobble supplies at least half the observed
velocity) AND total faked velocity below the pair's own catalog
ceiling — compute P(RUWE-hot | faker) via the response curves
(excess R(u)·a_phot vs sigma_AL(G_j), base = unity). Report
per-pair and the 9-pair mean; then the binomial read of the
OBSERVED 2/9 hot against the faker prediction. This conditions
on the collapse's own mechanism and cancels the absolute base:
a faker needs km/s-scale photocenter motion at these systems,
and the curves say whether such motion can hide from RUWE.
VALIDITY GL2' (the license): the CATALOG-CUT-REPAIRED S2 forward
at the operative world must land in [0.5x, 2x] of the measured
+0.174 — if the response model reproduces the sky's measured
activity-velocity correlation at the relevant amplitudes, its
(P, q) -> (leak, hot) mapping is sky-calibrated. GL2' FAIL ->
full stop, hold ~55%, the NSS leg decides (pre-stated).

BARS + MAP (locked; serves the round-5 branch the void 8L-a map
served, in conditional form):
  FAKER-LOUD: mean P(hot | faker) >= 0.6 AND binomial
    P(<= 2 of 9 hot) <= 0.05 -> the collapse account of the
    census pairs is object-level DEAD (they could not have
    faked it quietly); the S1 conditional is DISCHARGED at the
    census (S3-grade); anomaly-real ~55% -> ~57%.
  FAKER-QUIET: mean P(hot | faker) < 0.3 -> the fakers CAN hide
    from RUWE; the census blind spot is REAL at correction
    grade; anomaly-real ~55% -> ~50%; the NSS leg decides.
  GRAY: 0.3-0.6, or binomial P in (0.05, 0.5) under FAKER-LOUD
    means -> hold ~55%; NSS decides.
EXPECTED OUTCOME: none stated (the physics genuinely cuts both
ways: fakers need large a_phot at mid-P -> loud; but long-P
fakers ride L ~ 1 with R ~ 0 -> quiet; the (P, q) prior decides
the mix). Amendment rule: gate tolerances pre-quote; bars and
this map may NOT move.

## Stage 8L-a2 EXECUTED (2026-08-04): GL2' FAIL BY THE LETTER (2.6% over the band edge) — full stop, hold ~55%; the UNLICENSED faker table leans FAKER-QUIET and the NSS leg becomes THE decider

THE LETTER: the license gate required the catalog-cut-repaired
S2 forward to land in [0.5x, 2x] of the measured +0.174; it
predicted +0.357 = 2.05x — FAIL by 2.6% of the value (near-miss
disclosed; the amendment's catalog cut did massive work,
+2.84 -> +0.357: the mapping is roughly right at factor ~2).
Per the pre-stated branch: FULL STOP — no verdict from the
faker table, NO credence movement, hold ~55%, the NSS leg
decides.

THE UNLICENSED TABLE (reported as a lean, never a verdict):
P(hot | faker) = 0.11-0.23 across all nine census pairs (mean
0.142, min n_faker = 26,679); exact Poisson-binomial P(<= 2 of
9 hot | all nine are fakers) = 0.88 — the observed cleanliness
is CONSISTENT with all-fakers. THE PHYSICS (instrument-grade
understanding, quotable): the faking locus at these systems'
km/s-scale wobble requirements selects LONG periods (P >~ 5-10
yr), where the leakage curve L -> 1 (the wobble rides into the
fitted PM as a clean linear drift) and the residual curve
R -> 0 (nothing left for RUWE). The reviewer's round-5 blind
spot is CONFIRMED AT LEAN GRADE: RUWE cannot adjudicate the
census pairs' faker hypothesis; the channel that can is the
NSS/acceleration one (curvature-sensitive exactly at P ~ 3-15
yr) — 8K-b, exactly the pre-stated fallback.

THE HONEST ANNOTATION (carried, unlicensed): 8K's S3 leg
("census-clean = strongest, immune") is WEAKENED at lean grade
— its immunity assumed wobble-capable companions are
flag-visible; the faker table says the capable ones are quiet.
The 8K OBJECT-LEVEL-ABSENT verdict's S1 leg was already
conditionalized (round 5); S3 now shares a conditional of the
same shape. FORMAL handling: the 8K-b pre-registration MUST
price the branch where the quiet-faker reading is confirmed
(the +55 partially unwinds) and the branch where NSS finds the
census pairs orbit-clean at the faking periods (the collapse
account dies at the object level for real). No number moves
tonight.

STANDING VALUE regardless of the stop: the GL1 table — the
derived leakage L disagrees strongly with the legacy
S = min(1, P/17.8 yr) at mid-P (L = 0.80 vs S = 0.34 at 6 yr;
L = 0.98 at 17.8 yr where S first reaches 1) and the 17.8-yr
constant is NOT FOUND in the literature (scout; provenance
likely internal) — the 8L-b model-correction stakes are real
and now carry a candidate sign: the operative model
UNDER-weights mid-P wobble velocity.

Plain verdict: NEEDS REFINEMENT by the letter (the license
missed by 2.6% and the stop rule held) — with the round's real
product being instrument-grade physics: we now know WHICH
companions could fake the census (long-P), WHY the wobble flags
cannot see them (leak-as-PM), and WHICH instrument can (NSS
curvature). The chain 8K -> round 5 -> 8L-a -> 8L-a2 is the
honesty machine narrowing a question to a single named
measurement.

ELI12: We built the physics curve the referee ordered and asked
it: "could the companions that would fake our nine fast pairs
have hidden from the wobble flags?" The curve's first answer
failed its own driving test by a whisker (it predicted the
correlation twice as strong as measured, and our rule says
close only counts as a stop), so tonight's answer is formally
"no answer, belief stays 55." But the sketch it drew before
stopping is important and uncomfortable: the companions big and
slow enough to fake our pairs would wobble so smoothly that
Gaia reads it as ordinary drift — invisible to yesterday's
flags. So yesterday's "cleanest stars in the catalog" win is
weaker than it looked, and everything now rides on one last
channel: Gaia's orbit-curvature catalog, which sees exactly the
slow smooth wobbles the flags miss. That query is the next
thing we run — with the price of both outcomes written down
before we look.

## Stage 8K-b PRE-REGISTRATION (2026-08-05, committed BEFORE any query): THE NSS/ORBIT-CURVATURE CROSSMATCH — the decider named by 8L-a2, selection-aware

SCOUTED SELECTION FACTS (primary sources, locked into the design):
the DR3 astrometric-NSS pipeline processed ONLY sources with
RUWE > 1.4 and G < 19 (Halbwachs+23 sect. 1.2 verbatim);
acceleration solutions accepted at significance s > 20 (sect.
4.3); so for a RUWE-quiet source, NSS-absence is STRUCTURAL, not
evidence. The RV channel: radial velocities to G_RVS = 14 hard
limit; the documented variability criterion (Katz+23 sect. 3.7):
rv_nb_transits >= 10 AND rv_template_teff in [3900, 8000] K AND
rv_chisq_pvalue <= 0.01 AND rv_renormalised_gof > 4; RV trends
live in nss_non_linear_spectro (pulled as the 4th table).

DEFINITIONS (per component; pair = either component):
  astro-INFORMATIVE: ruwe > 1.4 (the pipeline's own entry gate —
    only there can NSS-absence inform);
  astro-ACTIVE: any nss_two_body_orbit or nss_acceleration_astro
    row, or non_single_star bit 1;
  rv-COVERED: the Katz applicability set (nb_transits >= 10, teff
    in window, finite gof/pvalue);
  rv-VARIABLE: covered AND pvalue <= 0.01 AND gof > 4;
  rv-TREND: any nss_non_linear_spectro row.
  Pair COVERED = rv-covered OR astro-informative on either
  component; pair ACTIVE = astro-active OR rv-variable OR
  rv-trend.

POPULATIONS: the 9 census pairs (decisive leg); WIDE (1437) and
a size-matched NARROW control (rng 41) as DESCRIPTIVE context
only (no bars — the selection functions are too opaque for rate
bars; stated).

BARS + MAP (locked; the census-pair leg):
  STRUCTURALLY-BLIND: covered < 6 of 9 -> the in-catalog
    channels are EXHAUSTED for the faker question; hold ~55%;
    T2 (external RVs of the nine) is named THE decider; the
    landing itself is the round's product.
  NINE-ACTIVE: covered >= 6 AND active >= 5 of 9 -> the
    quiet-faker account gains object-level support; anomaly-real
    ~55% -> ~50%; the seed budget re-prioritizes.
  NINE-CLEAN: covered >= 6 AND active <= 2 among the covered ->
    the faker account loses its last in-catalog hiding place at
    the covered periods; the S3 immunity partially RESTORES;
    anomaly-real ~55% -> ~57%. WINDOW DISCLOSURE carried: the
    RV channel covers the short-to-mid faking window plus slow
    trends at partial power; P >~ 15-20 yr fakers remain
    unprobed in-catalog (the all-nine-tuned-long-P conspiracy is
    priced at reading grade only), so 57 is this map's ceiling.
  MIXED: otherwise -> hold ~55%; T2 decides.
FETCH: TAP sync, chunked IN-lists, cached CSVs committed for
reproducibility; gates G8Kb-0 (population identity: 9/1437/1437,
unique ids), G8Kb-1 (pull completeness: gaia_source rows
returned for >= 99% of queried ids; else disclosed).
EXPECTED OUTCOME (pre-stated): no lean on active-vs-clean;
coverage expected partial (G_RVS = 14 vs our G ~ 8-15) — the
STRUCTURALLY-BLIND branch is live. Amendment rule: gate
tolerances pre-quote; bars and the map may NOT move.

## Stage 8K-b EXECUTED (2026-08-05, all gates PASS): NINE-CLEAN — the quiet-faker account loses its last in-catalog hiding place; THE MAP EXECUTES: anomaly-real ~55% -> ~57% (the map's ceiling)

Gates: G8Kb-0 populations exact (9 / 1437 / 1437-matched);
G8Kb-1 pull completeness 5748/5748 = 1.0000. The pulls (cached,
committed): 56 orbit rows, 8 acceleration rows, 2 RV-trend rows
across all three populations — the NSS products are SPARSE for
this catalog, exactly as the entry gates predict.

THE NINE: COVERED 7/9 (five via the Katz RV-applicability set;
two via astro-informativeness — the ruwe 1.48/1.53 pairs were
ENTRY-ELIGIBLE for the NSS pipeline and returned NO accepted
orbit/acceleration solution); ACTIVE 2 among the 7 covered =
exactly at the <= 2 bar -> NINE-CLEAN by the letter. The two
actives (s = 34.22 and 44.40 kAU, one RV-variable component
each) are DESCRIPTIVELY noted: at the WIDE baseline RV-variable
rate (0.074/pair) the expectation among 7 covered is ~0.5, so
2 is a mild excess (P ~ 0.09) — those two pairs may genuinely
host inner companions; THE CENSUS SURVIVES DROPPING THEM (a
7-pair band against the 4J/8F-b nulls remains astronomically
forbidden under Newton+errors). The two UNCOVERED pairs (8.82,
7.26 kAU) disclosed. Population context (descriptive, no bars):
WIDE astro-active 2.4% vs narrow-matched 2.0%; RV-variable
7.4% vs 5.8% — the fcomp ~ 0.1 sector's scale, nothing at the
collapse world's 35-50%.

THE VERDICT (locked bars + map): NINE-CLEAN — across every
in-catalog channel that CAN see them (RUWE flags where visible,
NSS orbits/accelerations where entry-eligible, RV variability
where covered), the census pairs show no companion activity at
the periods that could fake them. The S3 immunity partially
RESTORES; **anomaly-real ~55% -> ~57%**, the pre-stated ceiling,
with the window disclosure carried: P >~ 15-20 yr fakers and
the 2 uncovered pairs remain unprobed in-catalog (the
all-nine-tuned-long-P conspiracy stays priced at reading grade);
T2 external RVs remain the closure for that tail.

THE ARC LANDS: the collapse world (8J seed 31) needed a
0.35-0.50 companion population. Its light channel sits below
(v2c); its wobble-flag channel shows 9% where it needs 30 (8K,
conditionalized by 8L-a2); and its orbit/RV channels now return
the field rate, with the nine decisive pairs individually clean
at every covered period (8K-b). The fifth-move world is not
likelihood-excluded — it is POPULATION-STARVED at the object
level, in-catalog, on three independent channels, with the
residual windows named and priced.

Plain verdict: SUCCESS — the decider decided within its stated
window; the second earned credence movement of the program,
both under pre-registered maps.

ELI12: Final channel, final count. Gaia keeps two more lists we
had not checked: stars whose paths bend (hidden-companion
orbits) and stars whose speeds jitter (spectroscopic wobble).
Our nine fast pairs: seven are checkable, and five of those
seven are spotless; two show mild speed-jitter — and even if we
hand those two to the skeptics, the remaining seven still break
Newton's speed limit at astronomical odds. The two
highest-wobble pairs were exactly the ones Gaia's orbit-finder
was ALLOWED to inspect, and it found no orbit. So the
"companions faked everything" story now fails its head-count in
light, in wobble flags, and in orbits and speed-jitter — every
way the catalog can look. Belief moves 55 -> 57, the exact
ceiling we signed before looking, with the honest fine print:
very-slow companions (15+ year orbits) and two of the nine
remain uncheckable until we point real spectrographs at them.

## REVIEWER ROUND 6 (closing note, 2026-08-05): the closure ENDORSED — and the kernel debt named as the forward-model priority

THE NOTE: (1) the round-5 reservation was CONFIRMED and the
conditional bar worked as designed — holding at 55 then, and
weakening S3's immunity rather than banking the object test,
was the right call; the license-missed response model's lean is
CONSERVATIVE (a model that over-predicts the correlation and
still finds the fakers RUWE-quiet bounds true detectability
from above). (2) The fallback architecture is what makes the
closure robust: the verdict rests on channels that BYPASS the
RUWE blind spot (light; NSS/RV), not on the channel shown
blind. (3) SCOPE FORMALIZED (adopted into the record):
"NINE-CLEAN" is bar-defined, not literal — 2 of 9 carry
confirmed RV-variable companions, 2 are uncovered, and the
anomaly is robust to WORST-CASE REMOVAL OF ALL FOUR (>= 5 band
pairs against the 3.8e-9-class leakage null stays
astronomical). (4) THE LIVE THREAD WITH REAL STAKES = THE
KERNEL, not the census: the un-sourced 17.8-yr S-factor
conditions the ENTIRE forward-model machinery — the alpha band,
the collapse basin, the 8J saturation all ran on it — so the
8L-b correction is THE priority for anything forward-model,
and its headline deliverable is the SIGN: which direction the
corrected (boostier mid-P wobble) kernel pushes alpha. Our
note: the sign is not hand-derivable (two opposing arms — more
velocity per companion absorbs more excess, but the kinematic
binding per companion also strengthens, pinning fcomp lower);
that contest IS the 8L-b measurement. The standing hierarchy
reinforced: model-light/object-level results = the anchors;
forward-model alpha = the conditional layer with a NAMED
CALIBRATION DEBT (carried at every alpha quote until 8L-b).
(5) Credence: +2 within the committed map, "not something I'll
contest"; his view that a hold at 55 would have been equally
defensible (2-at-the-bar, license miss, kernel debt) is
CARRIED VERBATIM in the record — the 57 stands by the letter
of the pre-committed map, with that minority reading logged.

Consequence queue: 8L-b promoted to the top of the
forward-model queue (the derived L-curve replaces the legacy S;
headline = the alpha sign + the honest Newton band + the census
meter re-test); T2 = the external tail-closure; the 8-series
arc awaits paper absorption (PAPER still frozen at v3.9).

## Stage 8L-b PRE-REGISTRATION (2026-08-05, committed BEFORE any run): THE KERNEL ROUND — the derived leakage curve replaces the legacy S; the sign is the headline

PREMISE (round 6): the un-sourced S = min(1, P/17.8 yr)
conditions every forward-model number; the derived, gate-checked
leakage curve L(u) = 3(sin u − u cos u)/u^3, u = pi T/P,
T = 2.83 yr (GL0 1.6e-9; shape scout-confirmed) replaces it.
Parameter-free swap: w = wfac·v_orb·L/4.74047 in the companion
sector; rng streams untouched (S is deterministic in P).
Mode LKER=1 on the landed photow3 grid, e-sector pinned at the
8G mode (ein = 1.0, erf = 0.95); standard 9-dim cubes, TAG
'_lker', fresh names; seeds 31 + 101, both laws (~2 h GPU).

GATES (abort-grade):
  G8Lb-1 (identity sub-space): the fcomp = 0 column is
    kernel-independent — it must equal the esec (1.0, 0.95)
    slice's fcomp = 0 column EXACTLY (<= 1e-9; the G0-q
    precedent). GB0w/GB0e SKIPPED-disclosed (G8Lb-1
    substitutes).
  G8Lb-2 (kernel-magnitude report): population quantiles of
    L(P) vs S_legacy(P) over the actual drawn P distribution
    (the x2.4-at-6-yr claim measured on the draws).
  Reader G8Lb-3: the baseline read (the esec slice) must
    reproduce the shipped 8I-a OFF rows (0.57/+7.8, 0.54/+19.5,
    0.55/+3.4, 0.58/+21.9) to 0.01/0.1.

THE READS (calcs/stage8lb_read.py; LANDED-CONV anchor): per
law x seed — alpha_marg, dN, P(fcomp), P(kw), P(sq), P(fpm),
the Dwob analogue, all CORRECTED vs the esec-slice BASELINE;
the census forward at the corrected MAP cell (band_mu with the
corrected kernel; 8H convention and bar).

BARS (locked; seed means, both laws; NO credence movement in
any branch — the debt is being PAID, not adjudicated; the
in-pipeline credence path is closed per 8K-b/round 6):
  SIGN (the headline): d_alpha = corrected − baseline.
    |d_alpha| <= 0.11 -> WITHIN-SYSTEMATIC (the debt retires
    quietly; alpha quotes gain "kernel-corrected" status);
    d_alpha > +0.11 -> BOOSTIER-ALPHA (annotation, band shifts
    up); d_alpha < −0.11 -> QUIETER-ALPHA (the companion
    channel absorbs more; the honest band deflates —
    annotation; if alpha_marg <= 0.2 AND dN <= +5 both laws,
    the fifth-move sentence DEFERS to a powered round, as
    always).
  THE HONEST BAND: dN under the corrected kernel becomes the
    co-quotable Newton margin (supersedes the tanh-based
    "+9-11 under saturation freedom" conditional — that was
    the phenomenological form of this physical correction).
  CENSUS-RETEST: corrected-MAP-cell jointP >= 1e-3 in >= 3/4
    (the 8H bar) -> the meter REOPENS (context: the legacy
    C-row sat at ~1e-10).
  Dwob' reported (context: 7J-z4 +314; boostier mid-P wobble
    should strengthen per-companion binding — direction
    reported, not interpreted).
Context note (not a gate): the 8L-a light forward with this
same L over-predicted the S2 correlation 2.05x — the full
fitter carries the complete noise+selection model the light
forward lacked; its own gates are the identity/regression ones
above. EXPECTED OUTCOME (pre-stated, low confidence): no sign
lean — the two arms (more absorption per companion vs stronger
binding per companion) genuinely compete; that contest is the
measurement. Amendment rule: gate tolerances pre-quote; bars
may NOT move.

## Stage 8L-b EXECUTED (2026-08-05, all gates PASS): KERNEL-PAID — WITHIN-SYSTEMATIC both laws; the calibration debt retires; the honest Newton band FIRMS to ~+15

Gates: G8Lb-1 bit-exact 8/8 (the fcomp = 0 identity column,
0.00e+00 incl. vt); G8Lb-2 magnitude on the draws (mean L/S =
1.66 at P 3–18 yr; L q25 0.895 vs legacy 0.478 — the model had
been muting mid-P wobble by ~2x in the lower quartile); G8Lb-3
baselines 4/4 to the digit. Runtime ~2 h as estimated.

THE SIGN (the reviewer's round-6 demand, answered):
d_alpha = +0.006 (simple) / +0.044 (BE), seed means — WITHIN-
SYSTEMATIC by the locked bar. The two arms nearly cancel with a
mild positive residual (the binding arm edges the absorption
arm). CONSEQUENCE: every forward-model alpha quote GRADUATES to
kernel-corrected status; the round-6 calibration-debt annotation
RETIRES.

THE HONEST BAND: dN(corrected) seed means = +15.3 / +15.1 —
and the WEAK reads firmed (simple-31 +7.8 -> +11.0; BE-31 +3.4
-> +8.4 = the low-water mark more than doubled). This SUPERSEDES
the tanh-based "+9–11 under saturation freedom" conditional (the
physical correction replaces its phenomenological stand-in).
Mechanism notes: Dwob' eases 300-309 -> 249-259 (boostier wobble
per companion = fewer companions needed per unit kinematic work
— directionally sensible, still enormous); kw pins at the 0.7
floor at P = 1.0 (the amplitude nuisance absorbs part of the
kernel boost — grid-edge noted, the correction-#4 standard:
extension = a future decision); fcomp holds 0.1; sq holds 0.2.
CENSUS-RETEST: NOT reopened (0/4; the corrected kernel slightly
worsens the cliff, mu_hi 18.4–20.3 — as the boostier mid-P
spike flood predicts). The 8H companion-sector diagnosis stands.

NO credence movement (pre-stated; external-only per 8K-b/round
6). No fifth-move shape (alpha interior 0.54–0.63 everywhere).

THE 8-SERIES BINARY ARC IS COMPLETE: the width object measured
and defended (8G/8G-b), the census attributed (8H), the
survival and saturation repairs refused (8I-a/8J), the collapse
world population-starved at object level (8K/8K-b), the
response physics built and its blind spot mapped (8L-a/-a2),
and the kernel debt paid with the forward-model numbers intact
(8L-b). Remaining paths are external (T2/DR4) and galaxy-side;
the paper absorption of the arc awaits the author's call.

Plain verdict: SUCCESS — the correction with real stakes landed
within-systematic, which is the best kind of boring: the
program's numbers were never riding on the unfound constant.

ELI12: Our companion model had a hand-drawn dial nobody could
source — and the real physics curve says it muted mid-period
wobbles almost two-fold. Scary: every gravity number ran
through it. So we swapped in the real curve and re-ran
everything. Result: the gravity knob moved by less than a
hundredth (and the wobbliest-case readings actually got
STRONGER — Newton's worst deficit doubled). The scary dial
turned out to be load-bearing for nothing; our numbers now
stand on derived physics instead of an orphan constant. The
week's campaign is complete: every boring explanation tested,
each one either executed or priced, belief at 57 with two
earned moves, and the next battles belong to new telescopes
and the galaxy side.

## Stage 8M PRE-REGISTRATION (2026-08-05, committed BEFORE any run): THE JOINT-COHERENCE SCAN — the census-closure round (the reviewer's round-6 carried item)

THE QUESTION. Round 6 named one forward-model-internal
inconsistency blocking any "final" stamp on the α machinery: every
fitted world over-produces the perpendicular-ceiling census. At the
corrected-kernel MAP cells the forward pair is μ = (28.5–30.0,
18.4–20.3) against the observed (band = 9, cliff = 2), jointP
~1e-12. 8H probed six configs at the MAP cells (ALL-FAIL, companion
attribution 0.81/0.78); the CLASS was never swept. Three unmeasured
questions: (1) does ANY cell of the operative model class reproduce
(9, 2) at the 8H admissibility bar? (2) at what kinematic price?
(3) where does α sit when the kinematic likelihood and the census
face each other in one function?

THE INSTRUMENT. For every census-distinct cell of the operative
lker cube (α × η × w_rad × fcomp × fpm × kw × sq = 5·2·5·6·6·3·4 =
21,600 per law-seed; the fence axes fc0/ffly are data-side template
mixtures, census-blind by the shipped convention that G8M-0
anchors; wcut/ws slots singleton in photow3), forward the
corrected-kernel census pair (μ_band, μ_cliff) with the 8L-b
reader-verbatim band_mu machinery (orbit engine re-run per
(α, η, w_rad); one population per seed, lker=True; α=0 blocks
law-shared by construction and cached). Census term
CLL = ln Pois(9|μ_band) + ln Pois(2|μ_cliff).

ESTIMATOR CHOICE (locked): the PRIMARY floor for both admissibility
and CLL is the cell's own MC half-count resolution, μ_eff = max(μ,
0.5·Σ_bins NDATA/nk(cell)) — an MC zero is a resolution statement,
not a physical zero (every cell carries fpm ≥ 1.2 Gaussian noise;
its true tail rate is never exactly zero; the legacy 1e-12 floor
would convert MC resolution into −56/−261-grade lnL artifacts
exactly at the interesting boost-edge cells). The legacy 1e-12
convention (7K-b/8H verbatim) is CO-READ on every headline
quantity; divergence is itself reported. The 8H rows are floor-
insensitive (all probed μ measured nonzero) — no retro-tension.

PRODUCTS. (1) THE ADMISSIBLE SET: cells with jointP ≥ 1e-3 (the 8H
bar, primary floor): count, α composition, cheapest-kinematic-price
member (price = max cb − cb at the cell, fences profiled out).
Admissibility is evaluated WITHIN the operative prior support (Δkin
< 1e6 — LNPI-excluded fcomp cells are not the operative model); the
prior-blind census-only count is co-printed.
(2) THE TRADE FRONTIER: best CLL within kinematic ΔlnL ≤ {0, 2, 5,
10, 15, 20, 30, 40, 60, 100}, each with its α — the risk axis
shipped with the curve (the round-10 standing rule). (3) THE JOINT
DIAGNOSTIC: cbJ = cb + CLL (broadcast over fences); read9 →
α_joint, dN_joint, nuisance posteriors, joint-MAP cell + its
census pair; grid-edge flags mechanical (the correction-#4
standard; zero-floors of bounded axes are boundary-real per 8F-d).

THE DOUBLE-COUNT DISCLOSURE (locked label). The 9+2 census pairs
are members of the kinematic sample; their per-pair kinematic
weight is convolution-diluted but nonzero. The joint is therefore
an ACCOUNTING DIAGNOSTIC, not a new operative posterior: α_joint
and dN_joint are CO-QUOTED-DIAGNOSTIC, the operative band and the
8L-b corrected reads are untouched this round, and promotion of
any joint number to operative status is an external-review
decider-grade question (the PHYS-envelope precedent).

BARS (locked):
- B1 CLASS: N_adm ≥ 1 in ≥ 3/4 law-seeds → CLASS-CONTAINS; N_adm
  = 0 in ≥ 3/4 → CLASS-EMPTY (the wobble tail-SHAPE successor —
  the 8I-a successor, distribution-level — is then named MANDATORY
  before the α machinery is called final; the census flag stays
  permanent until repaired); 2/2 → CLASS-SPLIT, per-seed report,
  no closure claim.
- B2 PRICE (CONTAINS only; SPLIT reported ungraded): seed-mean
  cheapest Δkin to admissibility per law — ≤ 10 CHEAP / 10–40
  PRICED / > 40 SEVERE.
- B3 JOINT DIAGNOSTIC (always): α_joint, dN_joint per law-seed +
  seed means, CO-QUOTED-DIAGNOSTIC label only.
- B4 NEWTON FLANK: pre-stated expectation dN_joint ≥ dN_kin (the
  lker read) per law-seed; violated in ≥ 3/4 → AGAINST-EXPECTATION,
  quoted with the rescuing corner named.
- NO credence movement in ANY branch (pre-stated; the binary map
  is external-only per 8K-b/round 6 — internal-coherence
  accounting only).

EXPECTATION (pre-stated, two-sided): admissible cells exist in the
low-nuisance corner at PRICED grade; α_joint > α_kin (the census
demands the band be filled by something that dies at the cliff —
the boost's edge shape does that; the spike channel cannot: 8H's B
config made band 13–15 only with cliff 13–14); dN_joint ≥ dN_kin.
Openly uncertain: the α ≥ 1 bare cells were never forwarded — the
α=0.5 A-config starved at μ_band 0.86, and if the higher-α bare
cells under-fill too, the joint lands on companion-lite mid
corners and the α pull is weaker. Either branch is a measurement.

GATES: G8M-0 — recompute the 8L-b MAP cells from the cubes; the
direct band_mu call must match the four shipped CENSUS rows (±0.05
both μ; MAP params must match the shipped print). G8M-1 —
bare-Newton row (α=0, η=1.3, w_rad=0.3, fcomp=0, fpm=1.5, kw=0.7,
sq=0): μ_band < 3.0 (the 8F-b analytic scale; 8H notes 0.13–0.45),
all four law-seeds. G8M-2 — shape/finiteness/identity: cbJ.shape =
cb.shape; CLL finite; cbJ − cb at the joint MAP equals the
broadcast CLL exactly. G8M-3 (the GB0w-class wiring gate) — the
assembled scan block at the shipped MAP coordinates equals the
G8M-0 direct call to ≤ 1e-9 (identical expression sequence;
asymmetric axis sizes + the cell check detect transposition
scrambles).

Protocol: two seeds (31/101) per the standard; a 6-seed budget
only as a future decision if a verdict-grade object emerges (not
auto). Script [calcs/stage8m_jointcensus.py](calcs/stage8m_jointcensus.py);
output data/stage8m_jointcensus.txt; ledger row bin-8m-jointcensus
in the result commit; ~20–30 min compute.

## Stage 8M EXECUTED (2026-08-05, all gates PASS): CLASS-CONTAINS at SEVERE price — census-coherent worlds exist ONLY in the prior-forbidden no-companion corner; the operative α is census-ROBUST; B4 fires AGAINST-EXPECTATION with the mechanism named

Gates 16/16 first run: G8M-0 4/4 (direct calls reproduce the
shipped 8L-b census rows to the digit); G8M-1 4/4 (bare-Newton
μ_band = 0.00); G8M-2 4/4 (joint-cube identity); G8M-3 4/4
(scan-vs-direct bit-exact, e.g. 29.808979 ≡ 29.808979 — the
GB0w-class wiring gate did its job silently). Floor disclosure:
the pre-registered resolution-floor concern never bit — the
primary (half-count, median 0.103) and legacy (1e-12) conventions
return IDENTICAL results on every headline quantity (every
reachable cell's μ is measurably nonzero). 25.6 min compute.

B1 — CLASS-CONTAINS (4/4): N_adm = 1137–1209 of 21,600 (~5.3%)
per law-seed, all five α values represented. The composition is
the finding: the cheapest admissible cell sits at fcomp = 0 in
4/4 (α = 1 in three, α = 0.5 in one; sq = 0.2 in all; kw
irrelevant there), with μ ≈ (17.0, 4.4–5.0) and jointP
1.1–1.7e-3 — the best supported worlds still run ×1.9 over on
the band and clear the 8H bar only marginally (perfect-fit
jointP = 3.6e-2). NOTHING in the class actually fits (9, 2); the
admissible set is merely not-grossly-excluded.

B2 — SEVERE: admission price 5002.6 (simple) / 5010.3 (BE) lnL
seed means, dominated by the v2c host-rate prior + the kinematic
width demand against fcomp = 0. The census's ~24 lnL of leverage
cannot buy a ~5,000 lnL move: the joint MAP equals the kinematic
MAP in 4/4 (dkin = 0.0) and the model pays the census penalty in
full (CLL −25.3…−27.1 ≈ jointP ~1e-12).

THE FRONTIER (the risk axis shipped with the curve): flat where
it matters — within dkin ≤ 20 the census improves ≤ 3 lnL; the
−15-grade improvements need dkin ≥ 30–40; even dkin ≤ 100 reaches
only −13.5 (best possible −3.3). No cheap reconciliation exists
anywhere in the class.

B3 — JOINT-DIAGNOSTIC (CO-QUOTED-DIAGNOSTIC only, double-count
disclosed): α_J = 0.52 / 0.58 seed means (kinematic 0.56 / 0.61);
dN_J = +12.8 / +12.7 (kinematic +15.3 / +15.1). THE ROBUSTNESS
STATEMENT: folding the census into the likelihood moves α by
≤ 0.06 and costs Newton's rejection ~2.5 lnL seed-mean — the
operative numbers are census-robust; the inconsistency does not
destabilize the measurement. Nuisances unmoved (fcomp 0.1,
sq 0.2, kw 0.7 floor; BE-101 still rides fpm = 3.0 at P = 1.00,
edge flagged — the width-shape hunger is census-fold-invariant).

B4 — AGAINST-EXPECTATION (4/4): dN_J < dN_kin everywhere (−1.0 /
−2.1 / −3.9 / −2.7). Expectation-miss owned and resolved
mechanically: the pre-stated two-sided vice (band starves bare
Newton, cliff rejects spiked Newton) bites only in the
UNREACHABLE bare corner; within reachable mass every world
floods, and spiked-Newton floods marginally less than
spiked-boost (the 8L-b boostier-spike effect) — the frontier's
near-MAP best-census cells are α = 0 at several price points
(BE-31 at dkin 10–20, simple-101 at 20, BE-101 at 30). The
census credit at reachable distances goes to Newton; the
census's true preference (α ≈ 1, fcomp = 0) is priced out.

CLOSURE STATUS of the round-6 carried item: INSTRUMENTED AND
MEASURED, not repaired. The class reproduces (9, 2) only where
the measured host rate forbids it to stand; at its own optimum it
still over-produces the band ×3 with a ×10 cliff flood. The α
machinery GAINS the census-robustness annotation (≤ 0.06 / ~2.5
lnL exposure); the final-stamp REMAINS blocked pending the wobble
tail-SHAPE repair (the 8I-a successor, unchanged as the named
next binary instrument) or the external legs (T2/DR4). No 6-seed
budget triggered (4/4 uniform; price scatter 4986–5020; no
verdict-grade instability). NO credence movement (pre-stated;
external-only per 8K-b/round 6).

Worth keeping (mechanism note): the census's own preferred world
inside the class is the boosted no-companion one — α = 1 with
the width channel, the 4J reading reappearing in the forward
model — but at μ_band ≈ 17 it over-produces even there: the
model's tail machinery over-produces the band whenever it can
produce it at all. The (9, 2) pair stays a model-light data-side
statistic (defenses 4J 3.8e-9 / 8F-b 8.5e-8 untouched); the
forward model remains unable to make it, now at measured class
grade.

Plain verdict: SUCCESS as an instrument (the round-6
inconsistency is now a measured object with a price tag, a
mechanism, and an α exposure bound), NEEDS DIFFERENT PHYSICS as
a reconciliation (nothing in the class fits kinematics and
census together; the repair target is unchanged).

ELI12: Our star-pair model has a known problem: at its best-fit
settings it predicts ~30 pairs in the "faster than Newton
allows" speed band where the sky shows 9, and ~19 in the even-
faster zone where the sky shows 2. We had only ever checked a
handful of settings, so this time we checked ALL 21,600. Answer:
settings that get the sky's counts right DO exist — but only
with zero hidden companion stars, and our own measurements say
companions exist, so those settings carry a ~5,000-point
penalty. Folding the speed-band counts directly into the fit
therefore barely moves anything: the gravity knob shifts by less
than 0.06 and Newton claws back ~2.5 points of its deficit. Two
lessons: (1) our headline numbers don't secretly ride on this
flaw, and (2) the flaw is real and deep — the model's wobble
machinery overfills the fast band whenever it's allowed to fill
it at all, so the fix is reshaping HOW companion wobble is
distributed, not retuning knobs. And a nice touch: the settings
the speed band itself prefers have the gravity knob at 1 —
exactly the galaxy value — they're just unreachable while
companions exist.

## Stage 8N PRE-REGISTRATION (2026-08-05, committed BEFORE any run): THE FLOOD ANATOMY — anatomy-first replaces dial-guessing; the residual-survival candidate killed at design time

PART 0, THE DESIGN-TIME KILL (measured before this pre-reg, script
shipped in the stage for reproduction). The natural next repair
candidate after 8M was the RESIDUAL-SURVIVAL kernel: the catalog's
own quality cuts remove systems whose companion wobble sits in the
astrometric residual (R(u) from 8L-a), so the model should cull
them too — and the cut variable r = wfac·v_orb·R(u) has the
opposite P-shape to the leaked w = wfac·v_orb·L(u). The design-time
scale calc (population grade, rng 31) measured the premise DEAD:
**the band-flooding spike carriers (kw·w_leak > 0.3 km/s, 6.8% of
systems) have residual-to-noise ratio q50 = 0.5, q90 = 74.8 — the
flood is LONG-P (L→1, R→0), i.e. RUWE-SILENT by the same L/R
tradeoff that made 8L-a2's capable fakers quiet.** A survival cut
at any severity removes ≲10% of the flood carriers (the short-P
residual-loud ones, who are not the flood). SRVR is NOT BUILT; no
GPU is spent on a third refused dial. CONVERGENCE NAMED: the
model's flood population and the sky's RUWE-hiding population are
the same physics — so 8K's S3 cleanliness does NOT contradict the
model flood (the RUWE channel cannot arbitrate it), and the RV
channel (T2/8O — km/s offsets from exactly these long-P
companions) is sharpened as the arbiter.

THE ROUND. Before any further repair dial is proposed, MEASURE the
flood's anatomy at the operative cells (the four shipped 8L-b MAP
cells, corrected kernel):
- (a) THE P-LOCUS: among band and cliff members carrying an active
  companion, the dominant companion's period composition (bins
  < 3 / 3–10 / 10–30 / 30–100 / ≥ 100 yr + no-companion class).
  Pre-stated expectation (from the part-0 calc): ≥ 60% of band
  companion-carriers at P ≥ 10 yr.
- (b) THE CHANNELS: census pair recomputed with wobble-off (mass
  boost only), mass-off (wobble only), both-off — the 8H
  attribution refined to (band, cliff) × (wobble, mass, noise).
- (c) THE q-LAW DIAL: the paired flat-q vs TWIN-HEAVY t5 census
  (the v2c GV7-measured subsystem law; stream-preserving redraw,
  the 7J-z7 convention) at fcomp = 0.1 (MAP) and the freed 0.2 /
  0.35 — does the MEASURED wobble-quiet q-law reduce the flood, and
  does its mass-boost channel (twins boost ×1.2) re-flood what the
  wobble channel gives up (the 7J-z8 overshoot direction, there
  observed only in the FORCED fcomp ≥ 0.35 world)?

READING GRAMMAR (locked): R3 = REDUCES-MATERIALLY if at the MAP
fcomp the t5 pairing drops μ_cliff by ≥ 30% AND μ_band moves toward
9 without increasing — then the combined t5+lker operative mode is
NAMED as the next GPU instrument (own pre-reg). NEUTRAL/WORSENS →
the flood is population-structure in P (the successor = the
population-prior round anchored on the subsystem-period literature
scout now in flight, plus the external T2). Either way the
successor is evidence-named, not dial-guessed. NO credence movement
in any branch (measurement round; external-only per 8K-b/round 6).

GATES: G8N-0 — the full-channel census at the shipped MAP cells
reproduces the 8L-b rows (±0.05 both μ, 4/4; parse at runtime).
G8N-1 — t5 stream preservation: at fcomp = 0 the t5 and flat
censuses are EXACTLY equal (≤ 1e-9; companions inert, all other
draws share the stream). G8N-2 — channel identity: the both-off
channel equals the fcomp = 0 call exactly (≤ 1e-9). G8N-3 —
locus completeness: the P-class expected counts sum to the
full-channel band/cliff μ exactly (≤ 1e-9).

Protocol: seeds 31/101, both laws; ~10 min compute. Script
[calcs/stage8n_floodanatomy.py](calcs/stage8n_floodanatomy.py);
output data/stage8n_floodanatomy.txt; ledger row
bin-8n-floodanatomy in the result commit. Concurrent lane: 8O (the
T2 external-RV reconnaissance) pre-registers separately on the
scout-verified VizieR tables.

## Stage 8N EXECUTED (2026-08-05, all gates PASS): THE FLOOD HAS TWO CULPRITS — the band flank is the sq-TAIL, the cliff flank is the long-P companion spike; R1 CONFIRMED, R3 NEUTRAL; the successor is population-prior + width-shape

Gates 16/16 first run (G8N-0 shipped rows to the digit; G8N-1 t5
stream preservation exact; G8N-2 channel identity exact; G8N-3
locus completeness exact). 0.6 min compute.

PART 0 (in-record): spike carriers = 6.8% of systems, residual-
to-noise q50 = 0.46, q90 = 74.8; a wsr = 100 survival cut removes
7.7% of the flood carriers — SRVR dead as designed-killed.

R1 — CONFIRMED (mean 0.70 vs the ≥ 0.60 expectation): band
companion-carriers sit at P ≥ 10 yr in 0.67–0.72 of cases (cliff
similar) — the flood locus is the RUWE-silent long-P window in the
full forward, not just the population calc. The RV channel is the
arbiter (8O).

THE CHANNEL DECOMPOSITION (the round's discovery): the two census
flanks have DIFFERENT culprits.
- BAND: the no-companion channel alone makes μ_band ≈ 16.5–17.3
  (of the full ~29–30) at the MAP cell — vs 8H's A-config (sq = 0)
  at 0.86. **The sq = 0.2 width channel amplifies the bare band
  ×20 and carries ~57% of the band flood** — the lognormal smear's
  tail pushes sub-band pairs across √2 . The band flank indicts
  the WIDTH-CHANNEL SHAPE (the 3E→6P→sq→fpm-chase object's census
  face), not companions.
- CLIFF: companion share 0.77–0.81 (the 8H attribution reproduced,
  now flank-resolved): the long-P spike tail multiplies the bare
  cliff (~3.5–4.7) ×4–5 to 18.4–20.3.

R3 — NEUTRAL-OR-WORSENS by the locked grammar: the MEASURED
twin-heavy t5 law trims the cliff only −11…−17% (mean −14%, bar
≥ 30%) and the band −1…−5%; jointP improves ~20× and remains
~1e-10. The q-law is not the repair.

SUCCESSOR (named by the grammar + the scout, two-flank form):
(a) CLIFF flank — the population-prior round: the subsystem
P-distribution. Scout (Haiku, primary-flagged): Tokovinin 2014
(AJ 147, 87) measures inner subsystems in hierarchies SHORTER-P
than the field lognormal (excess < 30 d, nonmonotonic, tidal ×2
excess < 10 d); our model draws the RAGHAVAN FIELD lognormal for
subsystems — plausibly over-populating the exact P ≥ 10 yr locus
that carries the cliff flood; the inner-P distribution parameters
are NOT published at the needed grade (the MSC/Spectroscopic-
Orbits series holds the data) → the successor round starts with a
primary read/catalog fetch, not a dial. (b) BAND flank — the
width-channel SHAPE (bounded/truncated vs lognormal tail): the
census band is now a direct constraint on the sq distribution's
tail, connecting the width-object identification (TODO #18) to a
measured statistic. (c) The external arbiter — 8O/T2 (the RV
channel sees the long-P carriers both flanks share).

NO credence movement (pre-stated; measurement round). Ledger row
bin-8n-floodanatomy; worldtable tokens added.

Plain verdict: SUCCESS — the flood stopped being one mystery and
became two named, differently-shaped model defects with an
evidence-anchored repair queue.

ELI12: We dissected the model's overcrowded "too-fast" speed band.
Turns out two different things overfill it: the band itself is
mostly overfilled by our "measurement blur" knob (its bell curve
has too fat a tail — it smears ordinary pairs over the line),
while the even-faster cliff zone is overfilled by simulated hidden
companions on decade-long orbits. Neither knob we already tried
fixes this. But now we know exactly what to fix: use the REAL
measured orbit-size distribution for hidden companions (astronomy
books say they hug their stars closer than our recipe assumes),
and give the blur knob a physically-shaped curve instead of a
fat-tailed one. Also: those decade-orbit companions are invisible
to the wobble detector but LOUD to spectrographs — which is
exactly the instrument we point next.

## Stage 8O PRE-REGISTRATION (2026-08-05, committed BEFORE any query): THE T2 EXTERNAL-RV RECONNAISSANCE — the named external decider attempted on the public archives

THE QUESTION. 8K-b's map ended the in-catalog path at 57% with the
P ≳ 15–20 yr faker window + 2 uncovered pairs disclosed and "T2 =
the closure." The 8N part-0 kill SHARPENED the stakes: the model's
census-flood companions are RUWE-SILENT long-P systems — but they
are RV-LOUD (reflex amplitudes ~1–8 km/s at the faking locus). The
public RV archives (APOGEE multi-visit, GALAH, LAMOST, RAVE +
Gaia's own DR3 RV row) are the free first pass of T2.

TARGETS. All 23 ceiling-census pairs (46 components,
data/ceiling_pairs.csv); PRIMARY = the census region by the
corrected convention (band_corr ∈ {band: 9, above: 2} = 11 pairs,
22 components — includes 8K-b's 2 uncovered); the remaining pairs =
descriptive CONTROL rows. Positions from a Gaia DR3 TAP pull on the
46 source_ids; crossmatch via the CDS XMatch service (3″) against
the scout-verified VizieR tables: APOGEE DR17 III/286/catalog
(HRV, e_HRV, s_HRV scatter, Nvis, GaiaEDR3), GALAH DR3
J/MNRAS/506/150/stars (RVgalah, e_RVgalah), LAMOST DR7 V/156
(dr7slrs + dr7melrs; RV, e_RV, Nobs), RAVE DR6 III/283
(ravedr6/master fallback; HRV, e_HRV). Pulls cached and committed
(the 8K-b precedent; gitignore exceptions added).

CLASSIFICATION (locked): per covered component —
- RV-VARIABLE: within-survey multi-epoch scatter ≥ 1.0 km/s
  (APOGEE s_HRV with Nvis ≥ 2; LAMOST melrs if its epochs parse).
- OFFSET-FLAG: |RV_survey − RV_Gaia| > 3·√(e²+e²_G) + 1.0 km/s
  (cross-survey baselines span ~2003–2020 vs Gaia 2014–2017 —
  the only channel that can catch quasi-static long-P offsets).
- QUIET: covered, neither flag. Single-epoch rows with no Gaia RV
  = COVERED-SNAPSHOT (no read). Pair RV-ACTIVE if any component
  VARIABLE or OFFSET-FLAG.

BARS + MAP (locked; evaluated on the PRIMARY 11 only):
- C = # primary components with ≥ 1 external row. C ≤ 2 →
  T2-ATTEMPTED, NULL-COVERAGE: no move; T2 remains
  future-spectrograph work, now attempted-and-scoped.
- C ≥ 3: RV-CLEAN (0 active pairs) → +1 (57 → 58; small BY DESIGN:
  km/s thresholds + partial coverage + the within-survey-static
  blindness disclosed — a constant offset inside one survey is
  invisible; only cross-survey offsets catch the long end).
  RV-ACTIVE-MINOR (1 pair) → no move; the pair annotated.
  RV-ACTIVE-MAJOR (≥ 2 pairs) → −3 (57 → 54): object-level support
  for the faker channel inside the band; pairs conditionalized.
  Asymmetry argued in advance: weak-clean vs strong-hit at these
  thresholds and coverages. No other movement path.

GATES: G8O-0 target integrity (23 pairs; 9 band + 2 above by
band_corr; 46 distinct nonzero source_ids; Gaia position pull
46/46). G8O-1 pull integrity (every XMatch response parses;
per-survey row counts logged; caches committed). G8O-2 match
sanity (angDist ≤ 3″; nearest-hit dedup). G8O-3 column discipline
(classifications computed ONLY from the pre-named columns above;
an absent column → that survey SKIPPED-DISCLOSED; any renaming =
pre-quote amendment, never improvisation).

Protocol: network-only, minutes; no GPU. Script
[calcs/stage8o_extrv.py](calcs/stage8o_extrv.py); output
data/stage8o_extrv.txt; ledger row bin-8o-extrv in the result
commit. Known limits stated up front: APOGEE targets giants +
selected dwarfs (H-band), GALAH/RAVE southern, LAMOST northern —
expected coverage of 22 nearby field FGK/M components is LOW;
a null-coverage outcome is a booked scoping result, not a failure.

## Stage 8O EXECUTED (2026-08-05, amended run, all gates PASS): RV-CLEAN at thin coverage — THE MAP EXECUTES BY THE LETTER: anomaly-real 57 → 58; the control pair validates the instrument; T2's spectrograph leg remains the closure

TWO WIRING AMENDMENTS, both pre-quote, runs preserved: run 1 died
404 on a wrong XMatch endpoint path before any survey response
(_run1.txt; fixed to /api/v1/sync + per-candidate exception
handling). Run 2 printed NULL-COVERAGE and was VOIDED by its own
internal inconsistency (9 match rows returned, 0 registered): the
sid parse int(float(·)) TRUNCATED 19-digit Gaia source_ids
(float64 mantissa) so every hit landed under a corrupted id
(_run2.txt; exact-integer parse; skipped-survey rows now count as
coverage per the pre-reg's definition; G8O-2 hardened to parsed ==
returned — the vacuous-pass hole closed, the GB0w lesson again).
Bars/map/classification untouched throughout. Run 3: G8O-0 PASS
(23 pairs, 9+2; Gaia 46/46, RV for 45), G8O-1 all five tables
reached (APOGEE III/286 1 row, GALAH 1, LAMOST slrs 1, melrs 1,
RAVE III/283/ravedr6 5), G8O-2 9/9 parsed, G8O-3 three surveys
classification-skipped-disclosed on absent pre-named error columns
(GALAH e_RVgalah, LAMOST e_RV, melrs RV/e_RV — VizieR ships them
without; rows retained as coverage).

RESULT. Coverage C = 3/22 primary components (the bar's letter:
≥ 3), composition disclosed prominently: exactly ONE
classification-grade read — pair09 comp2 (band, 19.7 kAU):
**RAVE 2004 HRV = 54.6 ± 0.77 vs Gaia DR3 (2014–17) = QUIET under
the locked offset test — a ~decade-baseline external-instrument
stability of a census band component, the first non-Gaia datum in
the program** — plus two COVERED-SNAPSHOT rows (pair01 comp1
LAMOST; pair04 comp1 melrs with RV-pipeline sentinels — one of
8K-b's two uncovered pairs is now covered-but-unread). Control
context, right sign both ways: pair17 (BELOW-band) comp1 is
RV-ACTIVE (APOGEE s_HRV = 5.73 km/s over 4 visits) — the
machinery detects binarity when it is there, and found it outside
the census region; pair22 (below) both components RAVE-QUIET.

VERDICT (locked bars + map, executed by the letter): RV-CLEAN
(C = 3 ≥ 3, 0 active primary pairs) → **anomaly-real 57 → 58**,
the pre-priced small move (km/s thresholds; 1-of-3
classification-grade composition; within-survey-static blindness
— a constant offset inside one survey is invisible; 19/22 primary
components uncovered). The 8K-b window disclosure STANDS: the
P ≳ 15–20 yr locus and the two 8K-b pairs remain unread
(pair04's melrs row failed its RV pipeline). T2's real closure =
targeted spectrographs on the nine (two epochs, ~1 km/s, months
apart) — now attempted-and-scoped with the archive pass on
record.

Ledger row bin-8o-extrv; worldtable tokens; runs 1–2 preserved;
pulls cached and committed (targets + 5 xm CSVs).

Plain verdict: SUCCESS at reconnaissance grade — the archives
gave one real decade-baseline quiet read on a band pair, one
instrument-validating control detection, and a precise scoping of
what only new spectra can close; the map moved +1 by its letter.

ELI12: We asked every big public telescope archive: "have you
ever measured the speed-along-the-line-of-sight of our 22 special
stars?" Answer: mostly no (they're ordinary nearby stars big
surveys skip), but three had data. The best one: a speed
measured in 2004 agrees with the 2014–17 satellite value — that
star has NOT been yanked around by a hidden companion across a
decade. Bonus proof the method works: in a comparison pair
OUTSIDE our special band, the archive clearly caught a hidden
companion (its speed jumps by 6 km/s between visits). So: the
one deep look we got was clean, the detector demonstrably
detects, and belief ticks 57 → 58 exactly as the pre-signed
rulebook says. The full answer still needs us to point a
spectrograph at all nine pairs ourselves.

## Stage 8P PRE-REGISTRATION (2026-08-05, committed BEFORE any run): THE SQ-TAIL SHAPE CONTEST — the band-flank repair candidate

Premise (from 8N, measured): the census BAND flood is carried by the
WIDTH channel's lognormal tail (the no-companion channel alone floods
the band ~x20 at sq=0.2 vs the bare sq=0 forward), while the kinematic
likelihood demands sq>0 four independent ways (3E / 6P / the -60..-116
misfit / the fpm edge). The two-channel question this stage decides:
does the KINEMATIC demand extend to the lognormal's FAR TAIL, or only
to its bulk variance? If only the bulk, the band flank of the census
inconsistency is a shape-CONVENTION artifact and is repairable at zero
kinematic cost.

Instrument ([calcs/stage8p_sqshape.py](calcs/stage8p_sqshape.py),
committed with this entry): a single-config kinematic block evaluator —
the lnL_point ws=0/wcut-off legacy path lifted VERBATIM from
stage7j_marginal.py — plus the 8N census forward, both at the operative
lker MAP (alpha/eta/wr FROZEN at the cube argmax; no alpha
re-derivation this round, conditioning disclosed). Width draw
m = exp(sq*T(g)) with T a unit-variance DETERMINISTIC transform of the
same stream draw (stream-preserving; matched ln-m second moment):
logn (identity) / clip2 (clipped +-2sigma, rescaled) / ulog (bounded
uniform) / twopt (two-point) / lapl (HEAVIER tail = direction control).
sq axis extended to 0.5 for ALL shapes symmetrically (bounded shapes
may want more bulk variance; the extension is pre-registered here, not
post-hoc; cube identity gate uses the first four nodes).

Gates (any FAIL = STOP, no shape rows quoted, amendment pre-quote, run
preserved): G8P-0 in-evaluator identity (logn path = the verbatim
legacy arrays); G8P-1 cube identity (logn block sq<=0.3 vs the stored
lker cube at MAP, max|d| <= 0.05 lnL = keep-boundary grade, exact max
reported); G8P-2 census identity (vs stage8lb_read.txt mu, 0.05);
G8P-3 moment calibration (|var-1| <= 0.02, |mean| <= 0.01 per shape).

Bars (locked; law-seed majority >= 3/4):
- B1 SHAPE-ARTIFACT-CONFIRMED iff some bounded shape has profiled
  Dkin >= -2.0 vs logn AND Poisson pmf(9 | mu_band at its profiled
  cell) >= 1e-3 in >= 3/4 law-seeds, no law-seed at Dkin <= -5;
  winner = best mean Dkin among passing shapes (machine line
  'WINNER: <shape>' for the 8Q combined leg).
- B3 TAIL-DEMANDED iff EVERY bounded shape has Dkin <= -5.0 in >= 3/4
  law-seeds (then the kinematics and the census want OPPOSITE tails of
  the same object = a named model inconsistency, the sharpest possible
  successor statement). Else MIXED-CARRIED.
- B2 direction control: lapl floods >= logn at fixed nuisances in
  >= 3/4 (else DIRECTION-FLAG). Cliff cross-check: bounded shapes move
  mu_cliff <= 20% at fixed nuisances (8N attribution), breach = flag.
NO credence movement (measurement round, pre-stated). Successor if B1:
cube-grade re-run of the winner (own pre-reg) = operative-model
candidate. Expectation stated for the record: B1 leans plausible
(the profiled sq may climb to 0.4-0.5 for bounded shapes — flagged as
the compensation signature to watch); B3 would be the more important
result if it fires.

## Stage 8Q PRE-REGISTRATION (2026-08-05, committed BEFORE any run, same commit as 8P): THE SUBSYSTEM P-PRIOR BRACKET — the cliff-flank repair

Premise (from 8N, measured): the census CLIFF flood is long-period
companion wobble spikes (share 0.77-0.81; RUWE-silent locus 0.70 at
P >= 10 yr) drawn from the FIELD Raghavan lognormal logP[d] ~
N(5.03, 2.28). Source work (this session's scout, PRIMARY-SOURCE
grade at scout level; my own primary read owed before paper use):
Tokovinin 2014 (AJ 147:87) Table 2 MEASURES the inner-subsystem law
for components of wider pairs — L11: logP[d] ~ N(3.25+-0.12,
1.80+-0.09), ML detection-corrected (~80% completeness), median P
~4.9 yr vs the field's ~23 yr. His own caveat, carried verbatim: the
shortness is explainable by dynamical truncation and the sample is
NOT stratified by outer separation — for kAU outers (no stability
constraint) the true law is unconstrained between the two. Hence a
BRACKET, not a swap: raghavan (identity) / tokL11 (3.25, 1.80) /
mid (4.14, 2.04, labeled representative). The earlier scout claim
"params unpublished at needed grade" is hereby corrected (Table 2 has
them); the MSC fetch is SKIPPED by decision — Tokovinin 2018 states
the MSC "does not reflect the real statistics" (selection), so a
biased in-house histogram would be weaker than the published
ML-corrected law. MSC access route verified and logged for a future
stratified measurement (VizieR J/ApJS/235/6, systems.dat, Parent
field encodes hierarchy).

Instrument ([calcs/stage8q_pprior.py](calcs/stage8q_pprior.py),
committed with this entry): stream-preserving recast logP_new =
x0 + sg*z (same z-score draw; raghavan branch passes the untouched
legacy array bit-exactly); the same block evaluator + census forward
with the 8N P-locus decomposition, at the frozen lker MAP. Physics
pre-stated: the swap moves companion mass from the P >= 10 yr
spike/leak locus into the 1-10 yr mid-shoulder (P(P >= 10 yr | valid)
roughly halves) — the direction the kinematics may PREFER while the
cliff un-floods; it also shifts the fcomp<->host-rate mapping meaning
(LNPI kept at operative convention, pure-lnL co-printed, caveat
carried).

Gates (any FAIL = STOP + amendment pre-quote): G8Q-0 cube identity
(raghavan block vs stored cube <= 0.05); G8Q-1 census identity
(<= 0.05); G8Q-2 recast calibration (mean/std of logP within 0.02 of
target per slot); G8Q-4 shared non-companion streams (bit-identical
a_s/gs across priors); G8Q-5 locus completeness (1e-9).

Bars (locked; >= 3/4): Q1 CLIFF-REPAIRED iff tokL11 profiled
mu_cliff <= 0.5x raghavan's AND pmf(2|mu_cliff) >= 1e-3. Q2 tokL11
kinematics ACCEPTED iff Dkin >= -2.0 (REJECTED iff <= -5.0; else
CARRIED). Q3 secondary/informational (pre-stated fingerprint): the
kw=0.7 floor-riding and fpm=3.0 ceiling-riding should RELAX if the
companion-channel shape was their driver. Verdict = Q1 x Q2; mid rows
= the caveat envelope. NO credence movement (pre-stated).

COMBINED LEG (own invocation after the 8P record exists; rule
pre-stated, no discretion): width shape = 8P's WINNER line (none ->
logn), prior = tokL11; bar JOINT-COHERENT iff jointP >= 1e-3 AND
Dkin >= -5.0 vs the operative baseline in >= 3/4 law-seeds — this
re-asks 8M's admissibility question at the repaired model, POINT
grade. If coherent: 8M's inconsistency closes at point grade and the
cube-grade re-run (alpha under the repaired model) is the named
final-stamp decider — a successor round, NOT this one; no credence
move either way here.

### 8P/8Q AMENDMENT 1 (2026-08-05, logged PRE-QUOTE, all runs preserved): a 4e-17 constant — the program's first lnL-grade identity gate catches an ulp-divergent expression the census-grade gates were structurally blind to

Run-1 of both stages STOPPED by their own gates: G8P-1/G8Q-0 cube
identity missed at 2.6–10.1 lnL (law-seed-dependent, DETERMINISTIC —
bit-identical across the two concurrently-running scripts and across a
solo re-run; the concurrency hypothesis tested and excluded), while
census identity was simultaneously EXACT 4/4. Records preserved
(data/stage8p_sqshape_run1.txt, data/stage8q_pprior_run1.txt).

Diagnosis ([calcs/stage8pq_diag.py](calcs/stage8pq_diag.py) +
[calcs/stage8pq_diag2.py](calcs/stage8pq_diag2.py), records
committed): data-side objects, all population arrays, pick streams,
assembly, and orbit determinism are BIT-IDENTICAL between the marginal
and the reader; the marginal-driven lnL block reproduces the stored
lker cube at 0.000e+00 exactly; the SOLE remaining difference is one
expression — the reader lineage writes the radial-branch eccentricity
slope as the literal 0.045
([calcs/stage8lb_read.py:172](calcs/stage8lb_read.py#L172), inherited
verbatim by 8N and then 8P/8Q), while the marginal computes
erf+(0.995−erf) — and 0.995−0.95 differs from the literal 0.045 by
4.16e-17. The near-parabolic (e ≈ 0.95–0.995) orbit integration
amplifies that ulp-level eccentricity difference into macroscopic
phase shifts for a subset of systems → keep/bin-boundary flips →
whole-lnL-unit moves at probability-floor-amplified likelihood cells,
while every census statistic stays inside half-count tolerance. That
asymmetry is exactly WHY four consecutive census-grade gates (8L-b,
8M, 8N ×2) passed over it and the first lnL-grade identity gate caught
it on first firing — the GB0w corollary again: gates only catch what
someone thought to write.

FIX (this amendment): e_of_x in 8P/8Q rewritten to the bit-verbatim
marginal expression. NOT re-run: the shipped 8L-b/8N records STAND at
their stated tolerances — their conclusions (×3 flood factors,
~5,000-lnL admission prices, 0.70 locus fractions) sit orders above
boundary grade; this note is their provenance annotation, not a
supersession. NEW STANDING RULE (added to CLAUDE.md): reader scripts
copy model expressions bit-verbatim INCLUDING the arithmetic form of
constants; wherever a stored cube exists, identity gates are lnL-grade,
not census-grade. PREDICTION REGISTERED before the amended runs:
G8P-1/G8Q-0 → 0.00e+00 in all four law-seeds; census gates stay PASS
at ≤ 0.05 (the ship values carry the perturbed population, so a
nonzero last digit is allowed).

### 8P/8Q AMENDMENT 2 (2026-08-05, logged PRE-QUOTE, runs preserved as _run2): the census-identity bar re-baselined 0.05 → 0.10 — the amendment-1 prediction's census clause missed by one count-flip and is owned

The amendment-1 prediction landed BIT-EXACTLY on its primary clause —
amended G8P-1/G8Q-0 = 0.00e+00 in 4/4 law-seeds (the reader now
reproduces the operative lker cube exactly; the program's identity
chain reader↔cube is closed at machine precision for the first time) —
and HALF-MISSED its census clause: G8P-2/G8Q-1 failed on exactly ONE
of eight readings (BE-101 band: true-population 28.54 vs shipped
28.47 = 0.07 > the 0.05 bar; the other seven land 0.00–0.04). The
cause was already named in amendment 1: the SHIPPED census
(stage8lb_read.txt) carries the ulp-perturbed population, the amended
reader the true one; their difference is count-flip grade, and the
0.05 bar (inherited from G8N-0, where BOTH sides shared the perturbed
population) under-priced it by one flip — prediction miss owned.
Re-baseline: census-vs-ship bar → 0.10, explicitly a CROSS-POPULATION
consistency check; the bit-exact lnL gate is the primary identity
anchor. DIVIDEND BOOKED (the provenance annotation for 8L-b/8M/8N):
the TRUE-population MAP-cell census reads (29.81, 18.42) /
(29.11, 20.07) / (29.97, 20.29) / (28.54, 20.18) vs shipped
(29.81, 18.42) / (29.07, 20.07) / (29.96, 20.29) / (28.47, 20.22) —
max shift 0.07, three orders below every 8M/8N conclusion's margin.

## Stage 8P EXECUTED (2026-08-05, run 3 after amendments 1–2, ALL GATES PASS): SHAPE-INSENSITIVE — the census band flood is BULK-VARIANCE-carried, not tail-carried; the width-shape repair axis is DEAD and the width channel's distribution is measured out to ~2σ

Gates: G8P-0 2/2 (in-evaluator identity), G8P-1 cube identity =
0.00e+00 in 4/4 (the amendment-1 closure — reader ≡ operative lker
cube at machine precision), G8P-2 census 4/4 (0.10 cross-population
bar), G8P-3 moments 10/10. Verdict by the locked grammar:
MIXED-CARRIED, WINNER none — and the measurement inside the label is
sharp:

- THE BAND FLOOD IS THE SMEAR'S BULK, NOT ITS TAIL. Bounded shapes at
  matched Var(ln m) — clipped-2σ, hard-bounded uniform, two-point —
  leave μ_band at 27.3–30.7 vs lognormal's 28.5–30.0 (pmf9 never
  better than 3e-5; the B1 census leg fails 0/4 for every bounded
  shape) while clip2 is kinematically nearly FREE (Δkin −0.33 / −1.09
  / −1.99 / +0.26). At sq=0.2 the ±1σ factors (×1.22) walk the large
  sub-band population across ṽ=1.414 — no far tail needed. WORDING
  CORRECTION to 8N (mechanism word only, the measurement stands):
  8N's channel attribution (no-companion channel ≈57% of the band
  flood, ×20 amplification) is untouched; my interpretive word
  "sq-TAIL" is corrected to "sq BULK VARIANCE" — 8P was built to
  distinguish exactly this and did.
- THE LIKELIHOOD RESOLVES THE SMEAR'S SHAPE AT ±1σ AND IS BLIND
  BEYOND ~2σ: the two-point smear (continuum removed) is rejected
  −11.98/−15.03/−11.84/−16.57 in 4/4; clipping above 2σ costs ~0;
  the pre-registered sq-axis extension went unused (P(sq>0.3) = 0.00
  everywhere — no compensation hunger).
- THE LAPL CONTROL CONFIRMS THE MECHANISM BY ITS SIGN SPLIT: at fixed
  nuisances the heavier tail LOWERS the band (28.7 vs 29.8 — more
  central mass, less at ±1σ) and RAISES the cliff (21.0–22.3 vs
  18.4–20.3). B2's pre-registered direction (written under the
  tail-carrier premise) flags 0/4 — the flag marks the premise's
  failure, not the instrument's; the cliff cross-check breach (9/12,
  twopt −13…−21%) is the same physics (the cliff has a tail-fed
  margin the band lacks).
- CONSEQUENCE: at the variance the kinematics demand, NO reshaping of
  the width distribution reconciles the census band. The band-flank
  inconsistency is VARIANCE-level — the kinematic and census channels
  disagree about the same second moment. The TODO-18 width object
  gains a census SIDE-CONDITION: the physical width carrier must be
  LOCALIZED away from the band's feeder population (a universal
  multiplicative smear of the demanded size necessarily floods the
  band) — convergent with the 7J-z6 mid-shoulder/inner-bin
  fingerprint.

Plain verdict: SUCCESS as a measurement (bulk-carried, continuum
required at ±1σ, tail-blind beyond 2σ, repair axis closed) — NEEDS
DIFFERENT PHYSICS for the repair itself.

ELI12: Our simulated stars each get a random "blur" so they match
real measurement messiness. We suspected the blur's rare EXTREME
values were faking too many stars in our special speed band, so we
tried blurs with the extremes chopped off — same average blur, no
extremes. Surprise: the fake crowd didn't shrink at all. It's the
ORDINARY blur values doing it — tons of stars sit just under the
band, and a routine 20% bump pushes them in. Chopping extremes
changes nothing (and the fit barely notices); only shrinking the
whole blur would work, and the fit refuses that because it needs the
blur elsewhere. So the clash is real and deeper than a fat tail:
whatever blurs real stars must NOT touch the stars near our band — a
strong clue about where the blur physically lives.

## Stage 8Q EXECUTED (2026-08-05, run 3 + combined leg, ALL GATES PASS): NO-REPAIR × REJECTED — the kinematics REFUSE the measured subsystem period law, the cliff flood survives it, and the combined leg leaves the 8M price standing

Gates: G8Q-0 = 0.00e+00 4/4, G8Q-1 4/4 (0.10), G8Q-2 12/12 (recast
exact to <0.005), G8Q-4 2/2, G8Q-5 all; combined G8Qc-0 0.00e+00 4/4.

Results: Q1 = 0/4 → NO-REPAIR — the Tokovinin L11 law (3.25, 1.80)
does NOT halve the cliff: μ_cliff 18.4/18.4/20.9/26.3 vs raghavan
18.4/20.1/20.3/20.2; the P-locus shows the flood RELOCATING, not
dying (P≥100 yr mass drops, the 3–30 yr bins swell to ~5–6.6 each —
shorter-period companions spike the cliff just as hard through the
leak kernel's mid-P response). Q2 = 4/4 reject → REJECTED — Δkin =
−9.65/−11.34/−15.30/−15.99: the fitted model actively PREFERS the
field-like Raghavan long-P wobble mass. Q3 mixed (2/4): simple's
noise-edge hunger relaxes under tokL11 (P(fpm=3.0) 0.23→0.02,
0.07→0.07) while BE's does not (0.80–0.97); kw floor-pinned 12/12.
The mid bracket is rejected too (−9.29…−12.57) and floods MORE (band
30.7–33.6): the bracket is non-monotonic because P(valid) rises
toward shorter laws (0.563→0.728→0.879 — the a_in < 130 AU validity
cut was silently absorbing the field law's ultra-long tail). Under
tokL11 the profiled sq drops 0.2→0.1 in 3/4 — the added mid-P wobble
SUBSTITUTES for part of the width channel (a variance-conservation
echo of 8P's bulk finding: the likelihood defends a total width
budget whatever carries it).

Readings (labeled): (i) the rejection is P-SHAPE-carried, not
active-rate-carried — the doubled valid fraction at fixed fcomp could
cost only the known-flat fcomp-axis few-lnL, not −10…−16 (the
pre-reg's LNPI/fcomp-semantics caveat is carried; a semantics-matched
refit is the press-further instrument, unpromised). (ii) With 8N's
locus and 8I-a/8J, the wobble defect is bracketed as
AMPLITUDE-DISTRIBUTION-level within the long-P sector: survival cuts
(8I-a), saturation (8J), and now the period prior (8Q) all fail to
fix it.

COMBINED LEG (rule pre-stated; 8P WINNER = none → lognormal +
tokL11): JOINT-NOT-REACHED 0/4 — best jointP 4.6e-9 vs the 1e-3 bar
at Δkin −9.7…−16.0. THE ROUND'S SUMMARY SENTENCE: both evidence-named
census repairs are EXCLUDED at the operative posterior; the 8M
admission price STANDS; the (band=9, cliff=2) pair is
repair-resistant along both conventional axes and remains the
self-defending census statistic (7K-b form) — now hardened by two
more excluded alternatives.

Plain verdict: SUCCESS as a measurement (the evidence-named repair
candidate excluded on both its promised axes; the bracket's
non-monotonicity and the width-substitution effect are new model
facts) — the cliff flank stays OPEN.

ELI12: The textbook says hidden companions orbiting one star of a
wide pair tend to have SHORTER orbits than average field companions.
Our sim used the average-field rule, and we suspected that was why it
makes too many extreme fakes. So we plugged in the textbook's
measured shorter-orbit rule. Double surprise: (1) the fit got WORSE
by a lot — the data genuinely prefer the field-style wobble; (2) the
fakes didn't go away — shorter-orbit companions shake stars just as
hard, on a different rhythm. The suspect is cleared on both counts,
so the real problem stays where three earlier tests pointed: how
STRONG the simulated wobbles are, not how long their orbits take.

## Stage 8R PRE-REGISTRATION (2026-08-05, committed BEFORE any run): THE LOCALIZED-WIDTH CONTEST — the width-object instrument under the new census side-condition

Premise (both halves measured): 8P showed a UNIVERSAL smear of the
demanded variance necessarily floods the census band, whatever its
shape; the 7J-z6 Part-A fingerprint localizes the width DEMAND in the
mid-shoulder ṽ / radial (γ≈8°) column / inner bins — while the band
lives at γ ≥ 75°. Demand-sector and flood-sector are nearly ORTHOGONAL
in γ: a SECTOR-LOCALIZED width channel could satisfy both channels at
once. This is TODO-18's top instrument given its 8P side-condition.

Instrument ([calcs/stage8r_locwidth.py](calcs/stage8r_locwidth.py),
committed with this entry; the 8P machinery with a channel switch):
smear applied to (uni) all — the identity; (rad) the radial-orbit
subpopulation u_mix < wr — the PHYSICALLY NAMED carrier
(near-parabolic e>0.9 phase/eccentricity scatter); (gam) model γ<45° —
the fingerprint column at effective-description grade; (perp) γ≥75° —
the ANTI-control (maximal flood direction). sq grid extended to 0.9
symmetrically (a channel on fraction f needs sq ≈ 0.2/√f for
equivalent total width; pre-registered, not post-hoc).

Gates: G8R-0 structural uni-identity; G8R-1 cube identity = 0.00e+00
(the amendment-1 standard); G8R-2 census identity (0.10); G8R-3 rad
fraction = wr ± 0.01; G8R-4 mask-off wiring (rad forced-all ≡ uni,
bit). Bars (locked, ≥3/4): L1 LOCALIZED-VIABLE(rad/gam) iff Δkin ≥
−2.0 vs uni AND pmf(9|μ_band@profiled) ≥ 1e-3, no row ≤ −5 → THE
WIDTH OBJECT LOCALIZES (rad preferred as winner if both); successor =
cube-grade α re-run under the localized channel (the final-stamp path
REOPENS). L2 dose-response: perp floods ≥ uni (≥3/4). L3
GLOBAL-DEMAND iff rad AND gam both lose ≥5 in ≥3/4 (the sector
reading weakens; census inconsistency stays variance-level). Else
MIXED-CARRIED. Cliff co-printed (companion-carried; big response =
flag). NO credence movement. Expectation stated for the record: open;
rad is the motivated candidate; the honest risk is L3 (the width
demand may be bulk-global after all — that is what the instrument
decides).

## Stage 8S PRE-REGISTRATION (2026-08-05, committed BEFORE any run, same commit as 8R): THE GAS-DOMINATED c₁ — referee-queue T5, the pincer's M/L defense

Premise: the galaxy dial c₁ = 0.26–0.52 (4S flat / 4Z hier) carries
the standing W10-class objection that disk-M/L freedom manufactures
the deep-window coefficient. Gas-dominated galaxies are the standard
immunizer: where g_gas carries the budget, f_ML is inert and c₁ is
pinned by physics, not mass-to-light choices.

Instrument ([calcs/stage8s_gasc1.py](calcs/stage8s_gasc1.py); the 4S
machinery verbatim — ν_λ family with c₁ = λ/2, SPARC block,
marginalized objective): split at GDFRAC ≥ 0.5 (share of a galaxy's
kept points with g_gas > g_dsk+g_bul at f=1); fits SPARC-ONLY
marginalized (NO lensing — its 0.2-dex stellar-mass systematic would
contaminate the immunity claim): FULL reference / GD subsample +
200-rep galaxy bootstrap / DD complement contrast; the IMMUNITY
DEMONSTRATION (λ̂ with f_ML forced to 0.5 and 2.0; GD band |Δλ̂| ≤
0.2 pre-stated); raw y<1 GD co-read. Gates: G8S-0 selector sanity
(≥15 GD galaxies else POWER-STOP); G8S-1 the 4S joint endpoints
reproduced ±0.5 (−8397.72 / −8341.95); G8S-2 GD+DD additivity 1e-6;
G8S-3 c₁(λ) = λ/2 series. Bars (locked): T5-DEFENDED iff GD bootstrap
P(λ>0) ≥ 0.95 AND GD Δ1 interval overlaps the full-sample interval;
T5-POWER-LIMITED iff P(λ>0) < 0.95 with the interval covering both;
T5-TENSION iff the GD interval excludes the full-sample λ̂. NO
credence movement; the galaxy map untouched.

## Stage 8S-b PRE-REGISTRATION (2026-08-05, committed BEFORE any run; 8S executed first — its record ships as-is below): THE EDGE RESOLUTION — grid extension + the gas-budget knob

8S landed loud beneath its label: the gas-dominated subsample rides
the family's lower GRID EDGE (λ̂ = −0.30 = the first node; the
bootstrap slams its own −0.4 parameter bound; c₁=½ rejected +31
within the subsample; the raw y<1 co-read agrees at +38.7), while
the disk-dominated complement carries the entire positive dial
(λ̂ = 1.13) and swings bound-to-bound under M/L forcing — the
immunity contrast demonstrating exactly the alleged degeneracy. Two
process facts logged with the record: (i) by the correction-#4
standard an edge-riding optimum is UNQUOTABLE until the grid
extends — the 8S "dial" numbers are edge-invalid; (ii) the 8S
verdict grammar had a HOLE (a one-sided Δ1 interval cannot fire
T5-TENSION; it printed POWER-LIMITED) — the letter stands, the
label is annotated, the hole is closed in 8S-b's grammar.

8S-b ([calcs/stage8sb_gasedge.py](calcs/stage8sb_gasedge.py),
committed with this entry): (1) LGRID extended to −2.0 (below λ=0
the family is a DEEP-COEFFICIENT PROBE, c₁ = λ/2 < 0, not a
physical member set — labeled; ν-positivity guarded and reported);
(2) THE GAS KNOB — gN = fg·g_gas + f·g_dsk + g_bul with fg free in
[0.7, 1.4] (the helium/molecular budget envelope: the one coherent
systematic that shifts gas-dominated galaxies specifically), run on
FULL/GD/DD. Gates: G8Sb-0 shared-node regression vs the 8S profile;
G8Sb-1 fg-slice identity (1e-9); G8Sb-2 ν-positivity node report.
Bars (locked): E1 LOCALIZED iff the extended-grid GD minimum is
interior with a two-sided Δ1 inside (−2.0, 1.5). E2 GAS-EXPLAINED
iff with fg free the GD interval overlaps the FULL free-fg interval
(the 8S edge attributed to the fixed gas normalization; T5's
defense conditional on the gas systematic, both quoted). E3
DIAL-TENSION iff with fg free the GD interval still excludes the
full-sample λ̂ (a REAL subsample tension — named for the galaxy
program, not resolved here). POWER-CARRIED otherwise. NO credence
movement (pre-stated).

## Stage 8R-b PRE-REGISTRATION (2026-08-05, committed BEFORE any run; 8R executed — GLOBAL-DEMAND, its record ships below): THE COMPLEMENT CHANNELS — the γ-decomposition of the width demand

8R's only-sector channels all died (−31…−42, 4/4 each; dose-response
4/4; gates incl. cube identity 0.00e+00 all PASS) — the width demand
is distributed across sectors. But 8R's design tested "smear ONLY the
demand sector"; the 8P census side-condition asks for the COMPLEMENT:
smear everything EXCEPT the flood sector. 8R's gam channel is one
complement cut (sparing γ≥45 costs 36–42); the finer cuts are
unmeasured — in particular whether the band's own γ≥75 slice carries
kill-grade demand. 8R-b
([calcs/stage8rb_complement.py](calcs/stage8rb_complement.py),
committed with this entry) runs uni / x60 (spare γ≥60) / xperp
(spare γ≥75 = THE side-condition channel; the band's feeders stay
unsmeared so the no-companion band flood collapses by construction —
the question is purely the kinematic price). Gates: the 8R set
(cube identity 0.00e+00; census 0.10; fractions; mask-off wiring).
Bars (locked, ≥3/4): V1 VIABLE(xperp) iff Δkin ≥ −2 AND pmf9 ≥ 1e-3
(no ≤−5 row) → the reconciling channel EXISTS, successor =
cube-grade α under xperp (final-stamp path reopens). V2 NEAR-VIABLE
iff −5 < Δkin < −2 dominates with pmf9 ≥ 1e-3 (trade quantified,
carried). V3 SECTOR-DEMAND iff Δkin ≤ −5 in ≥3/4 → with 8P and 8R
this CLOSES the pincer on the multiplicative-smear class (no
γ-shaped per-system ṽ-broadening reconciles both channels; the
width object is NOT a multiplicative ṽ-smear — additive/noise-side
or data-side classes take the floor). x60 = the curve point.
NO credence movement (pre-stated).

## Stage 8R EXECUTED (2026-08-05, first run, ALL GATES PASS): GLOBAL-DEMAND — every only-sector width channel dies at kill grade; the kinematic width demand is distributed across γ

Gates: G8R-0 2/2, G8R-1 cube identity 0.00e+00 4/4 (the amendment-1
standard holds in fresh machinery on first firing), G8R-2 4/4,
G8R-3 fractions 4/4 (rad population fraction = wr exactly; kept-sector
fractions: rad ≈ 0.17, γ<45 ≈ 0.53, γ≥75 ≈ 0.15), G8R-4 mask-off
wiring 4/4.

Results: rad (the physically named near-parabolic carrier) loses
−33.6/−31.1/−37.9/−39.4; gam (γ<45, the fingerprint column — which is
simultaneously one COMPLEMENT cut: it spares γ≥45) loses
−37.7/−35.6/−39.3/−42.3; the perp anti-control loses −37.3/−35.5/
−40.3/−38.0 AND floods the band worst (L2 4/4 — the dose-response
confirms the mechanism). L1 0/4 for both localized channels; verdict
GLOBAL-DEMAND by the locked grammar. Notables: (i) simple-31's rad
cell (fpm=3.0, sq=0.3) reaches μ_band = 15.9, pmf9 = 2.2e-2 =
band-admissible — the localized channel CAN clear the band, at a
kinematic price the likelihood refuses (the ~30-lnL gap between
census-coherence and the kinematic optimum appears again, third
independent direction); (ii) crippled-width cells drift to
fcomp = 0.20 — the fit substitutes companions for width (the
width-budget defense's third appearance: sq ↔ wobble ↔ fcomp are
partially interchangeable absorbers). The 7J-z6 sector reading
WEAKENS: the fingerprint localized the residual misfit, not the
width channel's support.

Plain verdict: SUCCESS as a measurement (the γ-distribution of the
width demand is now half-mapped and the only-sector class is dead);
the reconciling question moves to the complement (8R-b).

ELI12: We tried letting the mystery blur touch ONLY the stars it
seemed aimed at (the ones moving along the line to their partner).
The fit hated it — it wants blur everywhere, not just there. We also
tried blur ONLY on the sideways-movers (where our too-crowded speed
band lives): hated too, and it crowds the band even worse — which at
least proves the dial turns the way we thought. So the blur the data
demand is everywhere-ish; the one cut we haven't priced is the exact
reverse: blur everyone EXCEPT the sideways band zone. That's the
next (and last) knob of this kind.

## Stage 8S EXECUTED (2026-08-05, first run, ALL GATES PASS) + 8S-b EXECUTED (the edge resolution): T5 NOT OBTAINED — the gas-dominated subsample lands a REAL DIAL-TENSION; the gas-budget suspect is cleared; the hier-GD leg is the named decider

8S gates: G8S-0 selector 38 GD galaxies / 422 points (PASS), G8S-1
the 4S joint endpoints reproduced to the printed digit (−8397.72 /
−8341.95), G8S-2 additivity 0.00e+00, G8S-3 series exact. Results:
FULL SPARC-only marginalized λ̂ = 0.919 (c₁ = 0.460 — reproduces the
4S joint dial internally); GD λ̂ = −0.300 = THE GRID EDGE (one-sided
D1 ..−0.258; c₁=½ rejected +31.0 and c₁=0 rejected +7.2 WITHIN the
subsample; raw y<1 co-read agrees, edge at +38.7); DD complement
λ̂ = 1.129 (c₁ = 0.565) carrying the entire positive dial. THE
IMMUNITY INSTRUMENT WORKED EXACTLY AS DESIGNED: GD |Δλ̂| = 0.000
under f_ML forced 0.5→2.0 (perfectly M/L-immune), DD swings 1.800 =
bound-to-bound (the alleged degeneracy demonstrated in the control).
Letter verdict T5-POWER-LIMITED — WRONG IN SPIRIT, two process flags
logged with the record: the edge-riding optimum is unquotable per
correction-#4, and the grammar had a one-sided-interval hole (could
not fire T5-TENSION); both annotated, the hole closed in 8S-b.

8S-b (pre-reg 95f752f; gates G8Sb-0/1/2 ALL PASS incl. the fg-slice
identity to 1e-9 and no ν-positivity exclusions): the extended grid
LOCALIZES the GD optimum INTERIOR at λ̂ = −1.542 (D1 −1.76..−1.31,
c₁ ≈ −0.77 — sharper than standard-μ, outside every physical family
member); THE GAS KNOB DOES NOT RESCUE — with fg free ∈ [0.7, 1.4]
the GD fit sits at fĝ = 1.171 interior and λ̂ = −1.732 (the
helium/molecular-budget suspect CLEARED), while FULL and DD ride the
fg = 0.7 edge (wanting less gas — noted). VERDICT: DIAL-TENSION by
the locked grammar — the M/L-immune subsample excludes the
full-sample dial under the FLAT treatment.

Readings (labeled, ranked): (i) the VERTICAL/DISTANCE channel is the
lead suspect — GD galaxies are precisely the worst-distance dwarfs,
the 4W identifiability boundary's home terrain, and the program's own
history is direction-consistent (flat 0.90 → hier 0.516 → GD-flat
deep-negative: more careful treatments pull λ̂ down); the 5T
ultra-deep POINT-level vote wanted c₁=½ under hier treatment —
point-level and galaxy-level splits disagree, which is itself the
vertical-channel signature. (ii) selection/inclination systematics of
dwarfs. (iii) genuine law shape — not adjudicable at flat grade.
NAMED DECIDER (queued, not run tonight — the hier machinery is the
program's heaviest and gets a fresh session): 8S-c = the λ profile
through the 4Z/5M vertical-hardened hierarchical machinery restricted
to the GD subsample. Until then the 4S/4Z dial quotes stand (they
carry the hier treatment), annotated with this tension.

Plain verdict: NEEDS REFINEMENT (the instrument worked, the immunity
demonstration is a keeper, but the T5 defense was NOT obtained and a
real tension is now on the books with its decider named).

ELI12: To prove our galaxy number isn't an artifact of guessing star
weights, we re-measured it using only gas-rich galaxies, where the
guess barely matters. The test itself worked beautifully — changing
the star-weight guess by 4× moves their answer by exactly zero. But
their answer came out very DIFFERENT from everyone else's — so
different it fell off our dial and we had to extend it. We checked
the obvious excuse (maybe we counted the gas itself slightly wrong):
freeing that changes nothing. What's left: these are the smallest,
farthest-blurriest galaxies, where distance errors — the thing our
fancier fitting machinery handles and this quick fit doesn't — hit
hardest. Next time: rerun them through the heavy machinery. Until
then, the number wears a warning sticker, exactly as the honesty
rules require.

## Stage 8R-b EXECUTED (2026-08-05, first run, ALL GATES PASS): SECTOR-DEMAND — THE PINCER CLOSES ON THE MULTIPLICATIVE-SMEAR CLASS

Gates: G8Rb-1 cube identity 0.00e+00 4/4, G8Rb-2 4/4, G8Rb-3
reported (spared fractions x60 ≈ 0.31, xperp ≈ 0.15), G8Rb-4 4/4.

Results: xperp (spare the band's γ≥75 sector) loses
−13.43/−12.54/−17.87/−18.58 = kill grade 4/4 — the band's own sector
demands width; its census CLEARS by construction (μ_band 15.5–16.9,
pmf9 = 1.4–2.6e-2, 4/4 band-admissible) — the geometry works exactly
as the side-condition said, and the kinematics refuse the price. x60
loses −27.8…−34.4. THE CURVE: kinematic cost ≈ 1 lnL per % of γ-mass
spared (15% → ~15, 31% → ~31, 47% [8R's gam] → ~38): the width
demand is ≈ UNIFORM IN γ — a new quantitative model fact.

THE THREE-STAGE PINCER (8P × 8R × 8R-b, every identity gate
bit-exact): a per-system multiplicative ṽ-broadening CANNOT satisfy
the kinematic likelihood and the (band=9, cliff=2) census
simultaneously — (i) any SHAPE at the demanded variance floods the
band (8P); (ii) any only-sector γ-localization is kinematically dead
(8R); (iii) any flood-sector-sparing complement is kinematically
dead too (8R-b). The sq = 0.2 object every fit demands (P = 1.00) is
therefore a MIS-SPECIFIED EFFECTIVE STAND-IN: whatever the real
broadening physics is, it is not multiplicative-with-a-γ-profile.
Surviving classes, named: noise-side shapes beyond the 7J-z6 grid
(its floor leg edge-rode one-law = the open crack), data-side ṽ
systematics, and the MIXTURE class (a second sub-population — the
w_rad precedent). The census pair's referee record now includes:
fitted forward worlds (7K-b), twin-forced worlds (7J-z8), analytic
tails (8F-b), shaped smears (8P), localized and complement smears
(8R/8R-b), and the measured subsystem P-law (8Q).

Plain verdict: SUCCESS — the width object's CLASS question is
answered in the negative at three-stage grade; the identification
hunt narrows to noise-side / mixture / data-side, with the census
pair as the standing referee.

### Stage 8T PRE-REGISTRATION (2026-08-05, committed BEFORE any run): THE ADDITIVE FLOOR UNDER THE CENSUS — the noise-side crack, census-tested

The surviving class #1 gets its instrument: 7J-z6's floor axis
(σ_eff = √((σ·fpm)² + (ws/4.74047)²), error-independent km/s jitter)
gained +7.1/+10.6 one-law and EDGE-RODE its 0.045 km/s grid top; its
census behavior was never evaluated. 8T
([calcs/stage8t_floorcensus.py](calcs/stage8t_floorcensus.py)):
ws grid extended ×9 to {0, 0.045, 0.09, 0.15, 0.25, 0.40} km/s;
joint block (fcomp, ffly, fpm, kw, sq, ws) at the frozen lker MAP;
census forwarded with the same floor. Gates: G8T-1 the ws=0 slice vs
the cube = 0.00e+00 (verbatim-branch, GW0-style); G8T-2 census
identity 0.10; G8T-3 shared-expression. Bars (≥3/4): F1
NOISE-SIDE-VIABLE iff wŝ>0 AND sq̂ ≤ 0.1 (the floor SUBSTITUTES the
smear) AND pmf9 ≥ 1e-3 → the width object identifies noise-side
(cube-grade successor; final-stamp path reopens). F2
ADDITIONAL-DEMAND iff gain ≥ +5 with sq̂ ≥ 0.2 and census still
inadmissible (the z6 reading — floor on top of smear). F3
CRACK-CLOSED iff P(ws>0) < 0.5. Edge flag at 0.40 (no further
extension — ~9× the z6 edge, ≫ the Gaia formal floor). NO credence
movement.

### Stage 8S-c PRE-REGISTRATION (2026-08-05, committed BEFORE any run, same commit): THE DISTANCE-CHANNEL CONTROL — the GD tension under measured vertical priors

The physics that makes this analytic: a distance error moves a SPARC
point VERTICALLY (g_bar is distance-invariant at fixed flux; g_obs ∝
1/D) and inclination errors likewise (V_bar components
photometry-side); so the vertical channel at measured-prior grade is
the conjugate random-intercept model — per-galaxy x_i ~ N(0, s_i²),
s_i² = (e_D/D)²/ln²10 + (2e_inc/tan i)²/ln²10, marginalized in closed
form; the s=0 branch is the verbatim flat expression (bit-identity).
This is the DISTANCE-GRADE leg of the named 8S-c decider — the full
hier M/L leg stays named if the tension survives.
([calcs/stage8sc_gddist.py](calcs/stage8sc_gddist.py)). Gates:
G8Sc-0 OFF-branch bit-identity + OFF GD profile reproduces 8S-b's
−1.542 (one grid step); G8Sc-1 analytic-vs-numeric marginal at 3
probes (1e-6); G8Sc-2 s_i distribution reported. Bars: C1
TENSION-DISSOLVES iff GD(ON) Δ1 overlaps FULL(ON) (attribution =
the measured vertical channel; T5 conditional); C2 TENSION-SURVIVES
iff GD(ON) still excludes FULL(ON) λ̂ (hardens; full-hier = the
remaining control); C3 POWER-CARRIED. Co-read: FULL's own λ̂ shift
under the channel. NO credence movement.

## Stage 8T EXECUTED (2026-08-05, first run, ALL GATES PASS): ADDITIONAL-DEMAND — the floor is REAL and INTERIOR at 0.045 km/s, does NOT substitute the smear, worsens the census — and IDENTIFIES THE fpm CHASE

Gates: G8T-1 ws=0 slice vs cube 0.00e+00 4/4; G8T-2 census 4/4;
G8T-3 shared-expression. Verdict F2 (3/4 rows; the fourth misses only
the gain≥5 clause at +3.96): the floor is DEMANDED — P(ws>0) =
0.99–1.00 in 4/4, gains +6.39/+13.56/+5.96/+3.96, and wŝ = 0.045
km/s INTERIOR on the ×9-extended grid (the 7J-z6 "edge-riding" was
its optimum; the old grid was too coarse to show the peak) — but
sq̂ = 0.2 in 4/4 (NO substitution: the floor rides on top of the
smear) and the census worsens monotonically in ws (band 30→110
across the axis; +1–3 even at 0.045).

THE DIVIDEND — THE fpm CHASE IDENTIFIED: P(fpm=3.0) collapses from
0.23/0.24/0.07/1.00 (operative) to 0.01/0.04/0.01/0.06 once the
floor axis exists. The sky's noise-shape hunger — the E-arm chase,
function-blind at 7J-d, the "real width-SHAPE incompleteness" of the
arm suite — is substantially a ~45 m/s error-INDEPENDENT additive
velocity floor, not a multiplicative scaling of the formal errors.
(The 7J-z6 PHYS envelope pre-dated this measurement and capped ws at
0.015 — pointer noted for the next noise-ceiling read, not
relitigated here.) Width-object standing after 8T: noise-side is
measured as ADDITIONAL (+4–14 lnL, census-worsening); the
sq-equivalent object remains unidentified; MIXTURE and DATA-SIDE are
the last classes standing.

Plain verdict: SUCCESS as a measurement (the floor pinned interior;
the fpm chase largely explained) — the width identity itself stays
open at two classes.

ELI12: We gave the simulation permission to add a tiny fixed amount
of extra measurement noise to every star — not proportional to the
official error bars, just a flat 45 m/s of "stuff happens." The fit
grabbed it immediately (it's been begging for something like this —
that's why it kept cranking the error-scaling knob to its ceiling;
that knob now relaxes). But it did NOT give back the mystery blur,
and extra noise makes the fake-crowding worse, not better. So: one
long-standing puzzle solved (the noise-knob chase), the main mystery
narrowed to two suspects — a hidden second population, or something
in the data itself.

## Stage 8S-c EXECUTED (2026-08-05, first run, ALL GATES PASS): TENSION-SURVIVES — the measured vertical channel barely moves the gas-dominated optimum; the dial-tension hardens

Gates: G8Sc-0 OFF-branch bit-identity (−1.14e-13) + OFF profile
reproduces 8S-b's −1.542 exactly; G8Sc-1 analytic = numeric marginal
to 8 decimals at 3 probes; G8Sc-2 s_vert distribution: GD median
0.132 dex vs DD 0.084 (the dwarfs carry ~2× the vertical freedom, as
expected — and it still doesn't matter). Results (vertical ON):
FULL λ̂ = 0.960 (D1 0.848..1.074), GD λ̂ = −1.315 (D1 −1.451..−1.178,
INTERIOR), DD λ̂ = 1.331 (now ABOVE BE). The GD interval still
excludes the full-sample dial by ~2 units of λ — TENSION-SURVIVES.

Suspect ledger: gas budget CLEARED (8S-b), distance/inclination at
measured-prior grade CLEARED (8S-c). Remaining, ranked: (i) the
WITHIN-CURVE CORRELATION object (7B: lag-1 ρ ≈ 0.87 — every
point-level −2lnL nominally calibrated only; the 8S galaxy-level
bootstrap partially covers this and agreed at P(λ>0) = 0.015, but it
was bound-limited at the old grid — an extended-grid galaxy
bootstrap is the cheap next control); (ii) the full-hier joint
treatment (4Z/5M machinery on GD — though for gas-dominated
galaxies the per-galaxy M/L channel is inert by construction, so
its marginal value over 8S-c is limited); (iii) dwarf
selection/inclination systematics; (iv) genuine law shape in the
deep regime — the GD↔DD split (−1.3 vs +1.3) is now a named,
hardened structure for the galaxy program.

Plain verdict: NEEDS DIFFERENT PHYSICS or the correlation control —
the tension is real at measured-vertical grade and is now the galaxy
program's sharpest open object.

ELI12: We gave every dwarf galaxy its honest wiggle room for "maybe
we measured your distance or tilt wrong" — for these little
galaxies that's a lot of room — and re-asked. Their answer moved a
hair and still disagrees hard with everyone else's. So it's not the
gas bookkeeping, not distances, not tilts. Either the points within
each galaxy are so correlated that we're overcounting the evidence
(checkable), or gas-rich dwarfs really do follow a different-shaped
deep-gravity curve than star-rich spirals — which would be a big
deal, and now has a name and a next test.

ELI12: Last blur experiment: blur every star EXCEPT the ones in our
special sideways band. The overcrowding disappears — geometry says
it must — but the fit says no AGAIN: even the band-zone stars need
their share of blur. Across the three experiments the trap has
closed: NO version of "multiply each star's speed by a random
factor" — any distribution, any zone map — can please both the fit
and the star count at once. So the extra spread in real data must be
a different KIND of thing: added noise, a second hidden population,
or something in the data itself. Our little 9-and-2 star count has
now refereed its sixth wrong idea — small statistic, big stick.

## Stage 8U + 8V PRE-REGISTRATION (2026-08-05, committed BEFORE any run): THE MIXTURE INSTRUMENT (the width object's last-classes contest) + THE GALAXY-LEVEL BOOTSTRAP OF THE GD TENSION (the correlation control)

Round 4 of the autonomy session; both stages are MEASUREMENT ROUNDS —
NO credence movement, pre-stated. Ledger rows on completion:
bin-8u-mixture, gal-8v-gdboot.

**8U THE MIXTURE INSTRUMENT** ([calcs/stage8u_mixture.py](calcs/stage8u_mixture.py),
output data/stage8u_mixture.txt). The pincer (8P/8R/8R-b) excluded
every GLOBAL multiplicative smear (any shape at matched variance, any
γ-localization); 8T excluded noise-side substitution. Named survivors:
MIXTURE and DATA-SIDE. The mixture hypothesis is the one shape class
8P never touched: a ZERO-INFLATED smear — a fraction f_mx of systems
carries per-system lognormal broadening of scale sm, the rest are
clean (the w_rad precedent: the last such contest found a real 20%
sub-population). Implementation: p['u_sq'] = rng.random(N) appended
AFTER p['gs'] in build_pop (stream-preserving prefix — every earlier
draw bit-identical); smear line vts = vtn*exp(sm*(gk*msk)) with
msk = (u_sq < f_mx) as float; at f_mx = 1.0, gk*1.0 = gk exactly, so
the (f_mx=1, sm∈SQ_GRID) slice IS the legacy global smear bit-for-bit
— the identity gate. Grids: F_MX = {0.05, 0.10, 0.20, 0.35, 0.65,
1.00}; SM = {0.0, 0.10, 0.20, 0.30, 0.45, 0.65}; ws ∈ {0, 0.045}
(the 8T-measured floor carried as a free binary axis — the mixture
must beat the global smear GIVEN the floor, not by mimicking it);
fcomp/ffly/fpm/kw grids and the LNPI host prior verbatim from 8T;
frozen lker MAP cells (α=0.5, η=1.05, wr=0.3); seeds 31/101 × both
laws. Var(ln m) of the mixture = f_mx·sm²; the grid covers the
equal-variance locus of the global sq=0.2 (var 0.040): (1.0,0.2)
0.040 / (0.35,0.3) 0.032 / (0.2,0.45) 0.041 / (0.1,0.65) 0.042 /
(0.05,0.65) 0.021 — census co-read along this locus at operative
nuisances, vars quoted.
GAIN := max(block+prior) − max(f_mx=1 slice + prior) — the mixture's
lnL preference over the global smear under identical freedom.
Gates (any FAIL ⇒ STOP, preserve run, amend pre-quote): G8U-1 cube
identity — (f_mx=1, ws=0, sm∈{0,.1,.2,.3}) slice vs the stored lker
cube = 0.00e+00 bit-grade, 4/4; G8U-2 census identity at shipped
nuisances |Δμ| ≤ 0.10 both components; G8U-3 8T regression at lnL
grade — gain_repro from the (f_mx=1, ws∈{0,0.045}, sm≤0.3) slice vs
the 8T printed gains (+6.39/+13.56/+5.96/+3.96), bar 0.011.
Bars (locked; law-seed majority ≥ 3/4):
- M1 MIXTURE-VIABLE iff profiled f_mx ≤ 0.35 AND GAIN ≥ +3 AND
  pmf(9|μ_band@prof) ≥ 1e-3 ⇒ the width object identifies as a real
  minority sub-population; successor = cube-grade re-run with the
  mixture axes (final-stamp path reopens) + the physical-ID round.
- M2 CENSUS-BLOCKED iff f_mx ≤ 0.35 AND GAIN ≥ +3 AND pmf9 < 1e-3 ⇒
  kinematics prefer a mixture but the census excludes it as
  reconciliation; DATA-SIDE is the last reconciliation class.
- M3 SHAPE-INDIFFERENT iff GAIN < +3 ⇒ the likelihood cannot
  separate a mixture from the global smear at this grade; sub-read
  (pre-stated): M3a if an equal-variance-locus cell with f_mx < 1
  sits within 2 lnL of the free max AND is census-admissible
  (pmf9 ≥ 1e-3) — a census-open escape EXISTS; M3b if no such cell —
  the census excludes the whole locus.
- Else MIXED-CARRIED. Edge flags per correction-#4: sm=0.65 top,
  f_mx=0.05 bottom, ws grid restricted by design (floor pinned at
  its 8T-measured interior value; pre-stated, not an edge).

**8V THE GALAXY-LEVEL BOOTSTRAP** ([calcs/stage8v_gdboot.py](calcs/stage8v_gdboot.py),
output data/stage8v_gdboot.txt). The 8S-c suspect ledger ranks the
within-curve correlation object (7B lag-1 ρ ≈ 0.87) first: point-level
D1 intervals overcount. The honest resampling unit is the GALAXY.
PAIRED design: NBOOT = 300 replicates, rng 71; each replicate draws
149 galaxies with replacement from the kept 149 (GD 38 + DD 111,
counts gated); per replicate fit λ̂ on the full draw AND on its
GD-member instances AND on its DD-member instances (duplicated
galaxies enter the marginal once per instance — independent-intercept
bootstrap semantics; m2ll loops over the instance list). Statistic:
Δ_rep = λ̂_GD_rep − λ̂_FULL_rep (paired). Estimator: coarse λ grid
−2.0..1.5 step 0.25 (15 pts), warm-chained Nelder-Mead (warm + one
cold start per λ), parabolic refine at interior minima; a FAST
bincount-vectorized m2ll (per-galaxy sums via np.bincount) is the
bootstrap engine, cross-gated against the VERBATIM 8S-c expression
(the GB0w new-MODE precedent). Skip rule: subset < 10 instances ⇒
NaN + count (expected never at Binomial(149, 38/149)).
Gates: G8V-0 verbatim-lift regression — the OFF-branch probe value
equals 8S-c's printed −647.874131278 at 1e-6, GD/DD counts 38/111
exact; G8V-1 fast-vs-verbatim cross-gate at 20 random θ probes incl.
a duplicated-instance list, max|d| ≤ 1e-6; G8V-2 identity-replicate
estimator check — coarse λ̂ within 0.10 of 8S-c's fine values
(FULL 0.960 / GD −1.315 / DD 1.331) for all three sets; G8V-3 rng
bookkeeping (first-replicate multiset fingerprint printed).
Bars (locked): B1 TENSION-ROBUST iff P(Δ_rep ≥ 0) ≤ 0.05 AND the
95th percentile of λ̂_GD_rep < 0.960 (ONE-SIDED, explicitly — the
8S-b grammar lesson). B2 CORRELATION-ABSORBED iff P(Δ ≥ 0) ≥ 0.20 ⇒
the point-level D1 overcounted; the dial-tension downgrades to gray
and TODO 26b re-ranks. B3 GRAY-CARRIED else. Edge rule (pre-stated
asymmetry): GD low-edge (−2.0) censoring RAISES λ̂_GD toward the
dial, i.e. biases AGAINST B1 — a B1 pass is conservative; report
edge fractions per set; if GD edge fraction > 30%, append
ESTIMATE-CENSORED to the verdict. Co-reads: percentiles (5/50/95)
of λ̂_FULL/λ̂_GD/λ̂_DD/Δ; P(λ̂_GD ≥ λ̂_DD); FULL high-edge fraction.

## Stage 8U EXECUTED (2026-08-05, first run, ALL GATES PASS): MIXED-CARRIED by the letter — the data pick a MAJORITY-mixture the minority-bars didn't cover (grammar miss owned); the pre-registered locus co-read measures the pincer running INSIDE the mixture family

Gates: G8U-1 cube bit-identity 0.00e+00 in 4/4 (the f_mx=1 slice IS
the legacy global smear); G8U-2 census identity ≤ 0.07; G8U-3 the
8T gains reproduced at lnL grade (+6.385/+13.556/+5.957/+3.961 vs
+6.39/+13.56/+5.96/+3.96 — the amendment-1 standard holding on a
fresh lift). 7.1 min.

Measurements (verdict letter = MIXED-CARRIED, M1 0/4 M2 0/4 M3 2/4;
rows stand): (1) the free optimum is the SAME cell in 4/4 —
**f_mx = 0.65, sm = 0.20, ws = 0.045** (var 0.026): a
majority-smeared / 35%-clean population over the 8T floor, gaining
+5.32/+5.28 (simple) and +1.14/+1.65 (BE) over the global smear;
P(f_mx<1) = 0.82–1.00; P(ws>0) = 0.99–1.00 (the floor demand
unconditional, 8T confirmed). My M1/M2 bars said MINORITY
(f_mx ≤ 0.35); the data went majority — the letter falls through to
MIXED-CARRIED and the miss is owned here (the 8S one-sided-interval
precedent: bar grammar must anticipate the whole outcome space).
(2) The profiled cell still FLOODS the census: band μ = 25.2–28.5
vs observed 9, pmf9 ≤ 1.3e-4, 4/4. (3) THE LOCUS TABLE (pre-stated
co-read): across the equal-variance locus the two demands are
monotone-OPPOSED — every kinematically-NEAR cell (d ≥ −2) floods
(pmf9 ≤ 9.5e-5); the only census-admissible cells (f_mx=0.05,
sm=0.65: pmf9 1.3e-3–3.8e-3) are kinematically dead (−39.2…−44.8);
**the pre-registered M3a flag (a NEAR + admissible + f_mx<1 cell
exists) fired 0/4.** No grid cell satisfies both demands in any
law-seed.

Standing: the mixture class as RECONCILIATION is measured-excluded
at every grid cell tried — stated as measurement, not bar-verdict
(the letter miss above); DATA-SIDE is the last reconciliation class
standing. The partial-mixture + floor cell (0.65, 0.20, 0.045)
replaces global sq=0.2 as the best kinematic-side effective
stand-in. The (9,2) census pair's referee count: SEVEN model
classes now vetoed.

Plain verdict: NEEDS REFINEMENT as a formal contest (the bars
mis-anticipated the optimum's location) / SUCCESS as a measurement
— the locus table answers the round's question: no zero-inflated
smear reconciles kinematics and census either.

ELI12: We tested the "hidden second family" idea: maybe only SOME
star pairs get the mystery blur. The fit likes a version of that —
blur about two-thirds, leave a third clean. But the trap door
closes: to stop the fake overcrowding of our 9-star corner, the
blurred family must be small and wildly blurred — and the fit HATES
that version (worst one tried). Every version the fit likes
overcrowds the corner; every version the corner allows, the fit
rejects. So "hidden family" joins the pile. What's left: the blur
isn't in the stars at all — it's something about the data itself.
(Honesty note: we wrote the rulebook expecting a SMALL hidden
family, the data picked a BIG one, so the referee sheet formally
reads "no ruling" — the measurements still say what they say, and
the rulebook bug is logged.)

## Stage 8V EXECUTED (2026-08-05, first run, ALL GATES PASS): TENSION-ROBUST — 0/300 galaxy-level replicates reach the dial; the correlation suspect is CLEARED and the GD↔DD split is bootstrap-hardened

Gates: G8V-0 OFF-probe −647.874131278 reproduced (d = −3.05e-10),
counts 38/111 exact; G8V-1 the fast bincount engine vs the verbatim
8S-c expression max|d| = 3.64e-12 over 20 probes incl.
duplicated-instance lists (the GB0w new-MODE rule); G8V-2 the
coarse estimator reproduces all three 8S-c fine optima
(0.969/−1.309/1.337 vs 0.960/−1.315/1.331); G8V-3 fingerprint
logged. 300 paired replicates, 0 skips, 6.3 min.

Results: **P(Δ ≥ 0) = 0.0000 — not one replicate lifts the
gas-dominated optimum to the full-sample dial**; Δ percentiles
5/50/95 = −3.00/−2.18/−1.22; λ̂_GD = −2.000/−1.254/−0.333 (median ≈
the 8S-c point −1.315; low-edge 8%, direction pre-stated
conservative — censoring can only weaken the pass); λ̂_FULL =
+0.31/+0.96/+1.50; 95th pct λ̂_GD = −0.333 vs the dial +0.960 ⇒ B1
passed with ~1.3 λ to spare. P(λ̂_GD ≥ λ̂_DD) = 0.0000 — the split
itself is galaxy-level robust. Co-read flag: the DD complement
rides its +1.5 grid top in 47% of replicates (λ̂_DD =
+0.53/+1.36/+1.50) — the star-rich side wants MORE than the grid
gives; the split is if anything wider than quoted (cosmetic for the
bars, which are GD-vs-FULL).

Suspect ledger after 8V: gas budget CLEARED (8S-b) → vertical
channel at measured grade CLEARED (8S-c) → within-curve correlation
CLEARED at galaxy-resampling grade (8V). REMAINING: dwarf selection
systematics; genuine deep-regime shape. (The full-hier joint leg
stays named, limited marginal value for GD.) Three controls
survived — the GD↔DD split is now a HARDENED structure.

Plain verdict: SUCCESS as a control — the cheap suspect died in one
afternoon; the tension itself now reads NEEDS DIFFERENT PHYSICS or
a selection story, and it is the galaxy program's sharpest open
object.

ELI12: Skeptic's worry: "each galaxy's points move together, so
maybe you counted one galaxy's opinion twenty times and the dwarf
disagreement is fake." Honest fix: treat whole GALAXIES as the
voting unit — reshuffle which galaxies are in the sample 300 times
and redo the whole fit each time. If the disagreement were a
vote-counting bug, some shuffles would erase it. Zero out of 300
did. Three suspects checked, three cleared — what's left is either
"how these dwarf galaxies got picked" or "deep gravity really is
different in gas-rich dwarfs."

## Stage 8W + 8X pre-registration (2026-08-05, committed before any run): the quality-strata tracking test (the data-side instrument) + the regime-vs-composition decomposition (the GD tension)

Both measurement rounds; no credence movement. Ledger rows on
completion: bin-8w-strata, gal-8x-regime. (First round under the
new lean-notes rule.)

**8W** ([calcs/stage8w_strata.py](calcs/stage8w_strata.py)): if the
width object is data-side, sq should track data-quality strata; if
flat everywhere, the last class loses its tracked handle. Four
median-split axes with pre-named physics shadows: A1 pair mean G
(shadow: mass model), A2 min parallax S/N (shadow: distance), A3
|ecliptic latitude| (scan-law proxy — physics-blind, the clean
discriminant), A4 max RUWE (shadow: multiplicity; skipped with a
note if the column is absent). Per stratum: rebuild data
histograms, noise pools, acceptance templates; model population
shared per law-seed at the frozen lker MAP; block (fcomp, ffly,
fpm, kw, sq, ws in {0, 0.045}); statistic = posterior-mean sq per
stratum; dsq = hi − lo per axis.
Gates: G8W-1 unsplit ws=0 slice vs the lker cube bit-identity;
G8W-2 stratum count conservation (exact); G8W-3 8T gain regression
at lnL grade (bar 0.011). Bars: D1 TRACKING iff some axis has
|dsq| ≥ 0.08 same-sign in ≥3/4 law-seeds (carrier named per the
shadow table; A3 = the clean data-side carrier). D2 FLAT iff all
axes |dsq| ≤ 0.05 in ≥3/4 (data-side loses its last tracked
handle; the width object graduates to "unattributed effective
absorber", DR4 the arbiter). D3 MIXED-CARRIED else. Co-read: the
observed band pairs' strata memberships (the in-script selector
must reproduce (9,2); else the co-read is skipped — pre-stated).

**8X** ([calcs/stage8x_regime.py](calcs/stage8x_regime.py)): is the
GD↔DD split composition (gas) or regime (deep)? Deep := gN(f=1) <
1.2e-10, fixed pre-fit. Six identity fits (GD, DD) × (all, deep,
nondeep), vertical-ON, the 8V fast engine with point masks; paired
galaxy bootstrap 200 reps (rng 71) over the three deciding subsets
(DD-deep, GD-deep, DD-nondeep — pre-stated cost cut). Decider =
DD-deep. Gates: G8X-0 unmasked GD/FULL reproduce 8V's printed
coarse values (bar 0.002); G8X-1 masked fast-vs-slow probes ≤
1e-6; G8X-2 mask accounting exact (deep + nondeep = all, per set).
Bars: X1 REGIME iff p95(DD-deep) < 0 AND P(DD-deep ≥ DD-nondeep)
≤ 0.05 — the split follows the deep regime (connects to 5T's
measured deep-arm vote; selection story weakens). X2 COMPOSITION
iff p5(DD-deep) > 0 AND P(DD-deep ≤ GD-deep) ≤ 0.05 — DD's deep
points vote with DD: the split is carried by galaxy type, not
regime. X3 GRAY-CARRIED else. Co-reads: GD-nondeep point fit,
counts, s_int per subset, edge fractions.

8X amendment 1 (wiring, logged pre-quote): run 1 died at the
GD-nondeep identity co-read fit — that subset has ZERO points (GD
galaxies are 100% deep; the confounding is total, which is the
measured fact that makes DD-deep the decider). Empty-subset guard
added; gates had all PASSED and no verdict quantity was generated;
bars untouched; run 1 preserved as the task log.

## Stage 8X EXECUTED (2026-08-05, run 2 after the wiring amendment, ALL GATES PASS): GRAY-CARRIED by the letter, COMPOSITION-leaning in substance — DD's deep points never vote like GD's (0/200); the regime hypothesis is dead

Gates: G8X-0 unmasked GD/FULL reproduce 8V exactly (−1.309/0.969);
G8X-1 masked fast-vs-slow 4.55e-12; G8X-2 counts exact — and the
accounting itself is a finding: **GD is 100% deep (422/422 points;
gas-domination and deep-regime are fully confounded on SPARC)**,
which is why DD-deep decides. DD: 1651 deep / 627 nondeep.

Identity fits: DD-all +1.337, DD-deep +0.772, DD-nondeep −0.507
(GD −1.309). Bootstrap (200 paired reps): DD-deep 5/50/95 =
−0.276/+0.697/+1.500; **P(DD-deep ≤ GD-deep) = 0.0000** — at
matched (deep) regime, galaxy type still separates completely;
P(DD-deep ≥ DD-nondeep) = 0.955 — within DD the deep arm pulls
POSITIVE, the opposite gradient from GD. X1 (regime) failed
outright; X2 (composition) missed only its p5>0 clause (−0.276) ⇒
GRAY-CARRIED by the letter, composition-leaning: the GD↔DD split
is carried by GALAXY TYPE, not by the acceleration regime.

Named next suspect (sharpened by this round): pressure-support /
asymmetric-drift corrections — gas-rich dwarfs are exactly the
class that needs them, SPARC rotmods don't include them, and the
bias direction matches (under-corrected V lowers g_obs at the deep
end ⇒ λ pulled negative). Cheap probe: split GD by V_flat — the
bias scales as (σ_gas/V)², so if GD-fast rejoins the dial it's
pressure support; if flat in V, genuine composition shape. That is
the 8Y candidate.

Plain verdict: NEEDS REFINEMENT by the letter (one clause short of
COMPOSITION) / the regime escape is closed as a measurement.

ELI12: We asked: maybe gas-rich dwarfs disagree just because we
only see their DEEP-gravity points, and everyone's deep points
would disagree? Test: the star-rich galaxies have deep points too
— do those disagree the same way? No. Never — 0 in 200 shuffles.
Same gravity depth, opposite votes: it really is about what KIND
of galaxy it is. Top remaining mundane suspect: little gas-rich
galaxies spin so slowly that gas pressure holds up part of their
weight, and the standard tables don't correct for that — testable
by checking whether the faster-spinning ones behave.

## Stage 8W EXECUTED (2026-08-05, first run, ALL GATES PASS): TRACKING — the width object is QUALITY-STRATIFIED: it concentrates in high-RUWE / low-parallax-S/N pairs; the physics-blind axis is flat

Gates: G8W-1 cube bit-identity 0.00e+00 4/4; G8W-2 count
conservation all axes; G8W-3 8T gains reproduced at lnL grade 4/4.
8.2 min. Census-pair co-read SKIPPED per pre-reg (crude selector
gives (39, 23) vs the operative (9, 2) — the operative census
carries cuts my mask lacks; honest skip).

dsq (hi − lo), 4 law-seeds: **ruwe +0.196/+0.116/+0.171/+0.152
(4/4)** — the clean-RUWE half (max RUWE ≤ 1.118) nearly ZEROES the
width object (sqbar 0.002–0.079; BE-101 even drops the floor,
P(ws>0) = 0.04); **plxsn −0.100 ×4 (4/4, one full grid step)** —
the high-S/N half wants sq = 0.1, the low half 0.2; eclat (the
scan-law, physics-blind axis) FLAT 4/4; gmag mostly flat (1/4
crosses). Verdict D1 TRACKING, carriers ruwe + plxsn. fpm was free
per stratum — the tracked width is not absorbable by error
scaling.

Reading (shadow table, stated carefully): the pure-geographic
data-side carrier (eclat) is CLEAN, so "generic data-side" is NOT
licensed; both tracking carriers have named shadows. What is
measured: the width object localizes in pairs with poor
astrometric solutions. Two surviving readings: (i) astrometric
noise beyond the formal errors in the poor half; (ii) unresolved
subsystems expressed through the ASTROMETRIC-SOLUTION channel
(photocenter motion corrupting the PM fit — exactly what RUWE
flags) rather than through the orbital-velocity wobble our
companion model implements — which would explain at once why the
velocity-space companion model cannot fit it (3K/7J: wrong
channel), the 8H companion-sector census attribution, and the
7I-S strict-cut behavior (RUWE < 1.2 ≈ this stage's median split:
the cut removes the width carrier, and the absorber field
reorganizes). Successors named: RUWE-continuous dose-response
(sq vs RUWE quantiles), and the cube-grade RUWE-stratified-sq
re-fit (α exposure under a localized width model); DR4 epoch
astrometry = the arbiter.

Plain verdict: SUCCESS — the width object has a measured ADDRESS
(poor-astrometry pairs); identity narrowed to astrometric noise vs
astrometric-channel multiplicity.

ELI12: The mystery blur isn't spread evenly after all — it lives
almost entirely in the star pairs whose position-tracking fits
were flagged as wobbly (and the far/faint ones), and NOT along the
satellite's scanning pattern. Two suspects left standing: those
pairs' measurements are just noisier than advertised, or hidden
companion stars are shaking the position fits themselves (the
flag is literally designed to catch that). Either way, we now
know WHERE the blur lives — that's a big step toward evicting it.

## Stage 8Y + 8Z pre-registration (2026-08-05, committed before any run): the pressure-support dose test (GD tension) + the RUWE dose-response (width object)

Measurement rounds; no credence movement. Ledger rows on
completion: gal-8y-pressure, bin-8z-dose.

**8Y** ([calcs/stage8y_pressure.py](calcs/stage8y_pressure.py)):
pressure-support bias scales as (σ/V)², so if it drives the GD
dial, the slow rotators must carry it. (a) V-split: GD halved at
the identity-sample median of V (Vflat from the SPARC table;
fallback max V_obs where Vflat = 0; counts printed); paired
bootstrap 300 reps over the 38 GD galaxies (threshold frozen;
halves fit if ≥6 galaxies; Δ_fs = λ̂_fast − λ̂_slow). (b) The
correction lever: refit λ̂_GD with g_obs → (V² + k·σ²)/R at σ = 10
km/s, k ∈ {0, 0.5, 1, 1.5, 2, 3} (flat-σ bound, crudeness
pre-stated: real AD corrections are radius-dependent; this brackets
the magnitude); DD at k ∈ {0, 2} = the bluntness control (DD is
fast — if DD moves much, the lever is too crude to read).
Gates: G8Y-0 the k=0 GD fit reproduces 8X's −1.309 (bar 0.002);
G8Y-1 fast-vs-slow engine probes ≤ 1e-6 (lg-argument variant);
G8Y-2 V-metric accounting. Bars: Y1 DIRECTION iff P(Δ_fs ≤ 0) ≤
0.05 (slow rotators carry the negative dial — pressure-support
direction confirmed); sub-clause Y1a REJOINS iff p95(λ̂_fast) ≥
0.960 (the fast half reaches the dial). Y2 FLAT iff P(Δ_fs ≤ 0) ≥
0.20 AND |median Δ_fs| ≤ 0.5 (no V-dependence; pressure support
disfavored; genuine composition shape strengthens). Y3
GRAY-CARRIED else. k-curve = co-read with pre-stated reading:
λ̂_GD crossing 0 by k ≤ 2 ⇒ plausible-magnitude corrections
neutralize the dial; no crossing by k = 3 ⇒ out of reach.

**8Z** ([calcs/stage8z_dose.py](calcs/stage8z_dose.py)): the 8W
carriers at quartile resolution — sqbar(Q1..Q4) for ruwe and
plxsn, blocks (fcomp, ffly, fpm, kw, sq, ws∈{0,0.045}) at the
frozen lker MAP. Gates: G8Z-1 unsplit cube bit-identity; G8Z-2
quartile count conservation; G8Z-3 8T gain regression (0.011);
G8Z-4 8W lineage regression — the (simple, 31) ruwe median halves
re-derived, sqbar vs 8W's printed 0.002/0.199 (bar 0.005). Bars
(per axis, ≥3/4 law-seeds): Z1 MONOTONE-DOSE iff sqbar(Q4) −
sqbar(Q1) ≥ 0.10 with no inversion > 0.02; sub-reads Z-STEP iff
Q4−Q3 ≥ 0.10 AND Q3−Q1 ≤ 0.05 (threshold behavior = a distinct
bad sub-population); Z-CLEAN iff sqbar(Q1) ≤ 0.05 (the clean
quartile carries no width — full attribution to the tracked
axes). Z2 NON-MONOTONE else ⇒ MIXED-CARRIED. The measured shape
parameterizes the cube-grade RUWE-stratified successor.

## Stage 8Y EXECUTED (2026-08-05, first run, ALL GATES PASS): DIRECTION + boundary-grade REJOINS — the GD dial is strongly V-ordered, but the direct correction lever is OUT OF REACH by k = 3

Gates: G8Y-0 k=0 GD fit = 8X's −1.309 exact; G8Y-0b lg_of_k(0)
bit-identity; G8Y-1 engine probes 1.82e-12; G8Y-2 V accounting
(Vflat 122 / fallback 27; GD threshold 66.2 km/s; 19/19). 3.9 min.

Two instruments, two answers: (a) V-SPLIT — Y1 DIRECTION fired,
P(Δ_fs ≤ 0) = 0.033 (conservative: GD-slow is edge-CENSORED at
−2.000, p50 = p5 = −2.000, which can only shrink Δ); identity
halves −2.000 (slow, LO-EDGE) vs −0.904 (fast). Y1a REJOINS fired
BY THE LETTER at the boundary: p95(GD-fast) = +0.971 vs the 0.960
bar — margin 0.011; the clause as designed encodes NON-EXCLUSION
of the dial in the bootstrap tail, not central rejoining (GD-fast
median −0.758; grammar-grading note, no retraction). (b) THE
k-CURVE — the flat-σ AD correction moves λ̂_GD only −1.309 →
−1.040 (k=2, σ_eff 14 km/s) → −0.925 (k=3, 17 km/s): **k* NOT
REACHED by k=3 for even λ = −0.5** (slope ≈ +0.13/k ⇒ λ=0 needs
k ≈ 10, σ_eff ≈ 32 km/s — the rotation speed itself; absurd); DD
control moves only −0.16 at k=2 (lever valid, not blunt).

Joint reading: whatever drives the GD dial RIDES THE V AXIS (the
slow half is maximally negative), but the pressure-support
correction at any physical σ CANNOT be the resolution — out of
reach by ~an order of magnitude in k. The V-ordering's carrier is
therefore not (only) the AD magnitude: candidates = the real
radius-dependent AD in the slowest dwarfs exceeding the flat-σ
bound locally (needs Σ_gas(R) — external data, bookable), a
V-correlated data-quality/inclination systematic, or genuine
shape concentrated in the slowest gas-rich dwarfs. The flat-σ
crudeness was pre-stated; the proper Σ(R)-weighted AD correction
is the definitive successor instrument.

Plain verdict: SUCCESS as a measurement — direction real,
resolution out of reach; the GD tension survives its FOURTH
control at correction-lever grade, now with a measured V-order.

ELI12: If gas pressure were faking the dwarfs' weird answer, the
slowest spinners should be weirdest — and they are! But when we
actually APPLY the pressure fix, even a generous dose moves the
needle only a fifth of the way home; you'd need gas pressure as
strong as the spinning itself, which is absurd. So: the weirdness
really does live in the slow spinners, but the textbook pressure
fix can't explain it. Either the fix is bigger than textbooks
allow in exactly those galaxies (checkable with better gas maps),
or slow gas-rich dwarfs genuinely behave differently.

## Stage 8Z EXECUTED (2026-08-05, first run, ALL GATES PASS): MIXED-CARRIED by the letter — and the quartile curves DISSOLVE the monolithic smear: bottom-step + plateau + an AGGREGATION component

Gates: G8Z-1 cube bit-identity 0.00e+00 4/4; G8Z-2 counts; G8Z-3
8T gains at lnL grade 4/4; G8Z-4 lineage — the (simple, 31) ruwe
median halves reproduce 8W's 0.002/0.199 EXACTLY. 8.5 min.

Measured curves (sqbar Q1→Q4, 4 law-seeds consistent):
ruwe 0.00 / 0.08 / 0.10–0.12 / 0.09–0.12 — a BOTTOM-STEP at
RUWE ≈ 1.05 then a plateau ≈ 0.10; **Z-CLEAN 4/4 (the clean
quartile carries NO width, ≤ 0.005)**. plxsn 0.12–0.16 / 0.20 /
0.12–0.19 / 0.02–0.09 — non-monotone (peak Q2, clean at best-S/N
Q4). Letter: ruwe Z1 2/4 (seed-31 rows miss Q4−Q1 ≥ 0.10 by
0.006–0.014), Z-STEP 0/4, plxsn 0/4 ⇒ MIXED-CARRIED. Two grammar
warts owned: Z-STEP was encoded at the TOP of the axis (the real
step is at the bottom), and Z1's sign convention cannot fire for
plxsn (quality rises with quartile there) — bars as designed,
findings quoted as measurements (8U precedent).

**THE AGGREGATION FINDING** (gate-anchored): the ruwe-hi HALF
demands sqbar = 0.199 (G8Z-4 exact vs 8W) while its constituent
quartiles want only 0.102 / 0.087 — the union demands ~2× the
width of any part; likewise the full catalog (0.2) vs the
quartile max (0.12). A single global error-scale (fpm) across
quality-heterogeneous strata MANUFACTURES per-system width. The
width object now decomposes: (a) zero in the clean quartile, (b)
a ≈ 0.1 plateau in the upper three, (c) an aggregation/
heterogeneity component from pooling, (d) the 8T floor (present
in every quartile incl. clean Q1 — floor and smear separate).
Successor sharpened: the cube-grade stratified re-fit should
carry fpm PER STRATUM (or an fpm(RUWE) law) + stratified sq —
expected from these curves: the global smear collapses toward
≤ 0.1; the α exposure is the number to measure.

Plain verdict: SUCCESS as a measurement (the object's anatomy:
clean-quartile zero + plateau + aggregation artifact) / the
formal dose grammar NEEDS REFINEMENT (both warts above).

ELI12: We split the star pairs into four quality grades and asked
each grade how much mystery blur it needs. Cleanest grade: NONE.
The other three: a modest amount, roughly equal. But here's the
trick finding: pour two grades into one bucket and fit them
together, and the bucket asks for DOUBLE — because one
noise-dial for a mixed crowd fits nobody, and the fit fakes the
difference as extra blur. So a chunk of our famous "0.2 blur"
was never real blur at all — it was mixing different-quality
data in one pot. The rest lives in the flagged pairs, and the
cleanest quarter of the sky needs none.

## Stage 9A + 9B pre-registration (2026-08-05, committed before any run): the stratified-noise α re-fit (the α-exposure decider) + the quality-flag control (GD)

Measurement rounds; no credence movement. Ledger rows on
completion: bin-9a-stratalpha, gal-9b-qflag.

**9A** ([calcs/stage9a_stratalpha.py](calcs/stage9a_stratalpha.py)):
the 8Z anatomy says the global sq = 0.2 is partly an aggregation
artifact of one error-scale over heterogeneous strata. This stage
measures what that repair does to α. RUWE quartiles (the measured
axis); per α ∈ {0, 0.5, 1, 1.5, 2} (η, wr frozen at MAP — BLOCK
grade, profile + LNPI; explicitly NOT comparable to the 7J
marginal numbers; direction-measurement, full marginal re-run =
successor if material): per-stratum blocks (fcomp, ffly, fpm, kw,
sq, ws∈{0,0.045}); TIED model = Σ_q blocks maximized with shared
(fpm, sq); FREE model = Σ_q max over per-stratum (fpm_q, sq_q) —
nested, same stratified data, no pooling artifact; ws and the
population axes shared (the floor is global per 8T/8Z). α̂ by
parabolic refine; report α̂_tied, α̂_free, Δα̂, ΔlnL(free−tied) at
α̂, per-stratum (fpm̂_q, sq̂_q), Newton contrast ΔlnL(α̂ vs 0)
under both models, edge flags.
Gates: G9A-1 unsplit block at MAP α, ws=0 vs the lker cube =
0.00e+00 (bit); G9A-2 stratum count conservation; G9A-3 8Z
regression — sqbar recomputed from the α=MAP stratum blocks vs
8Z's printed values, 16 checks, bar 0.0005.
Bars: A1 EXPOSURE-MATERIAL iff |Δα̂| ≥ 0.25 in ≥3/4 law-seeds
(the honest width model materially moves α; successor = the full
marginal re-run). A2 EXPOSURE-CONTAINED iff |Δα̂| ≤ 0.15 in ≥3/4
(the α measurement is robust to the width-model repair). A3
GRAY-CARRIED else.

**9B** ([calcs/stage9b_qflag.py](calcs/stage9b_qflag.py)): is the
GD dial carried by the SPARC quality flag? GD split Q=1 vs Q=2
(sizes printed; <8 galaxies either side ⇒ POWER-FLAG, pre-stated);
identity fits + paired bootstrap 300 reps (rng 71); Δ_Q = λ̂(Q1) −
λ̂(Q2). Co-read: the 8Y V-split repeated within each Q stratum
(identity fits only, power-limited). Gates: G9B-0 GD-all
reproduces 8X's −1.309 (bar 0.002); G9B-1 accounting incl. the 8Y
V-threshold regression (66.2). Bars: B1 QUALITY-CARRIED iff
P(Δ_Q ≤ 0) ≤ 0.05 AND median Δ_Q ≥ 1.0 (the high-quality subset
carries materially less dial). B2 QUALITY-BLIND iff P(Δ_Q ≤ 0) ≥
0.20 AND |median Δ_Q| ≤ 0.5 (the dial ignores the quality flag —
the data-quality suspect weakens). B3 GRAY-CARRIED else.

## Stage 9B EXECUTED (2026-08-05, first run, ALL GATES PASS): QUALITY-BLIND — the SPARC quality flag carries none of the GD dial; the V-ordering replicates inside both quality strata

Gates: G9B-0 GD-all = 8X's −1.309 exact; G9B-1 accounting (Q1/Q2
= 18/20, no power flag; V threshold 66.2 = 8Y exact). 3.8 min.

Identity: GD-Q1 λ̂ = −1.281, GD-Q2 = −1.222 — indistinguishable.
Bootstrap (300 paired): Δ_Q median = 0.000, P(Δ_Q ≤ 0) = 0.537 ⇒
B2 QUALITY-BLIND fired cleanly. Co-read (the finding): the 8Y
V-ordering REPLICATES within each flag stratum — Q1 slow/fast =
−2.000 (LO-EDGE)/−0.887; Q2 = −2.000/−0.638 — so the V-structure
is not a quality-flag artifact.

Suspect ledger after 9B: gas budget (8S-b) → vertical (8S-c) →
correlation (8V) → correction-lever (8Y) → quality flag (9B) all
CLEARED — FIFTH control survived. Remaining: radius-dependent AD
beyond the flat-σ bound (needs Σ_gas(R), external), or genuine
physics in the slow gas-rich dwarfs. The V-ordering is THE
structure to explain.

Plain verdict: SUCCESS as a control; the tension is now
five-controls hard, and the V-ordering is its measured shape.

ELI12: Maybe the astronomers' own "medium quality" sticker marks
the galaxies causing trouble? No — the high-quality and
medium-quality gas dwarfs give the same weird answer, and in BOTH
groups the slow spinners are the weird ones. Whatever this is, it
cares about how fast a galaxy spins, not about the sticker.

## Stage 9A EXECUTED (2026-08-05, first run, ALL GATES PASS): GRAY-CARRIED by the letter — and the direction is measured: the honest width model moves α UP or NOWHERE; the Newton contrast strengthens

Gates: G9A-1 cube bit-identity 0.00e+00 4/4 (at MAP α, fresh
orbit runs); G9A-2 counts; G9A-3 the 8Z sqbar regression 16/16 at
5e-4. 20.0 min.

Results (tied → free, per law-seed): α̂ 0.301→0.446 / 0.322→0.522
/ 0.272→0.472 / 0.422→0.420; **Δα̂ = +0.144/+0.199/+0.201/−0.002
— uniform-positive-or-zero**; letter GRAY (A1 0/4, A2 2/4 — three
rows landed in the 0.15–0.25 GAP between my bars; bar-design note
owned: the gap caught the actual values). **ΔlnL(free−tied) =
+88…+95 in 4/4** — the per-stratum error model is decisively
preferred (the 8Z aggregation artifact confirmed at full-block
grade). Per-stratum widths at α̂ replicate the anatomy EXACTLY:
sq = 0.0/0.1/0.1/0.3 in 4/4. Newton contrast (block grade, NOT
comparable to 7J marginal numbers): tied +0.6…+4.9 → free
+4.2…+9.1 — **the width repair does not erode the Newton
rejection; it strengthens it modestly.** Correction-#4 flag: Q4
(worst RUWE) rides BOTH noise axes at grid top (sq=0.3, fpm=3.0,
4/4) — the pathological stratum is grid-censored; successors:
extend Q4's axes in the cube-grade re-fit, and the Q4-DROPPED
robustness fit (α on Q1–Q3 only; named, not run).

The round's question answered at its pre-registered grade: **the
α measurement survives the width-model repair** — the feared
collapse direction did not occur; the honest model likes α
slightly MORE. The full marginal re-run stays available but its
trigger (A1) did not fire.

Plain verdict: SUCCESS as a measurement (direction + the ~90 lnL
heterogeneity confirmation) / the bar grammar NEEDS REFINEMENT
(the A1–A2 gap).

ELI12: Big worry: once we stop faking the blur (one noise dial
for everybody) and give each quality grade its own honest dial,
maybe the gravity signal was living in the fake blur and dies.
Tested: the honest version likes the data hugely more (as
expected), and the gravity knob moves a hair UP, not down — and
"no gravity boost" gets slightly MORE disfavored, not less. The
signal wasn't hiding in the blur. (One footnote: the very worst
data group wants even more noise than our dials allow — its
knobs are pinned at max; we flagged it for the big re-run.)

## Stage 9C + 9D pre-registration (2026-08-05, committed before any run): the radius-dependent AD correction (GD) + the Q4 robustness pair (binary)

Measurement rounds; no credence movement. Ledger rows on
completion: gal-9c-adcorr, bin-9d-q4robust.

**9C** ([calcs/stage9c_adcorr.py](calcs/stage9c_adcorr.py)): the
8Y flat-σ lever lacked the radius dependence; the real asymmetric
drift grows outward. Bounding model (crudeness pre-stated:
outer-disk exponential; HI central holes ignored): Σ_HI(R) =
Σ₀ e^(−R/R_g), R_g solved per galaxy from the catalog pair
(MHI, RHI) via Σ(RHI) = 1 M☉/pc² (brentq, fallback counted);
correction g_obs → (V² + σ²·R/R_g)/R (normalization-free slope;
constant σ). λ̂_GD(σ) at σ ∈ {0, 8, 10, 12, 15} km/s; DD control
at σ = 10; 100-rep GD bootstrap at σ = 10 (co-read interval).
Gates: G9C-0 σ=0 bit-identity to LGOBS0 + GD fit = 8X's −1.309
(0.002); G9C-1 engine probes ≤ 1e-6; G9C-2 R_g accounting (solved/
fallback counts, R_g percentiles, correction-fraction percentiles
at σ=10). Bars (at σ = 10, the physical value): C1 NEUTRALIZED
iff λ̂_GD ≥ 0.0 (the dial is pressure-support at catalog grade;
T5's defense collapses to conditional). C2 PARTIAL iff −1.0 <
λ̂_GD < 0.0 (material, not neutralizing). C3 OUT-OF-REACH iff
λ̂_GD ≤ −1.0 (the radius-dependent version also fails: the
pressure story is dead at catalog grade; remaining = external
Σ_gas(R) data or genuine physics). DD bluntness flag if
|Δλ̂_DD(σ=10)| > 0.3. Co-read: σ* crossings.

**9D** ([calcs/stage9d_q4robust.py](calcs/stage9d_q4robust.py)):
9A's Q4 rides both noise axes at grid top (correction-#4). Two
reads per law-seed on the free-model α profile: DROP (Q1–Q3 only)
and EXT (Q4 on extended grids sq ≤ 0.5, fpm ≤ 4.2; Q1–Q3
standard). Gates: G9D-0 cube identity at MAP; G9D-1 the
recomputed free_std α̂ equals 9A's printed values (bar 0.002,
4/4). Bars: D1 ROBUST iff |α̂_drop − α̂_std| ≤ 0.15 AND |α̂_ext −
α̂_std| ≤ 0.15 in ≥3/4 (the α measurement does not depend on the
pathological stratum). D2 Q4-CARRIED iff |α̂_drop − α̂_std| ≥
0.25 in ≥3/4. D3 GRAY-CARRIED else. Co-read: Q4-extended
(sq̂, fpm̂) interior-or-edge.

## Stage 9C EXECUTED (2026-08-05, first run, ALL GATES PASS): OUT-OF-REACH — and SIGN-INVERTED: the radius-weighted pressure correction makes the GD dial WORSE; the pressure story is dead in-catalog at both amplitude and shape

Gates: G9C-0 σ=0 bit-identity + GD = 8X's −1.309 exact; G9C-1
engine probes 1.32e-11; G9C-2 accounting — and the accounting IS a
finding: **R_g solved 0 / fallback 149**. The primary model
(exponential disk anchored at Σ_HI(RHI) = 1) is INFEASIBLE at
SPARC values — the peak of Σ(RHI) over R_g sits below 1 M☉/pc²
for typical (MHI, RHI); my design-time root-existence check was
wrong (owned — I worked an example that had no root and misread
it). The stage ran entirely on the pre-registered fallback R_g =
R_max/2 (outer slope dlnΣ/dlnR = −2 at the edge — a reasonable
outer-disk value); correction fractions at σ=10: GD median 3%,
90th 7%, max 34% of V². 0.7 min.

THE MEASUREMENT: λ̂_GD(σ) = −1.309 / −1.324 / −1.334 / −1.346 /
−1.371 at σ = 0/8/10/12/15 — **monotone in the WRONG direction**;
σ* NOT REACHED for any threshold; DD control −0.079 (no bluntness
flag); 100-rep bootstrap at σ=10: −2.000/−1.385/−0.327. Verdict
OUT-OF-REACH by the letter, a fortiori robust to the R_g
feasibility miss: ANY positive radius-weighted correction moves
λ̂ the measured wrong way. Physics: raising the deep outer points
tilts the GD arm toward even LESS boost — the dial is NOT the
signature of under-corrected pressure support. Combined with 8Y
(amplitude) the pressure hypothesis is dead at catalog grade;
remaining: genuine physics in slow gas-rich dwarfs, or a
systematic outside the pressure class (external Σ_gas(R) data
would now need the sign to flip, which no Σ shape provides).

Plain verdict: SUCCESS as a kill — the last cheap mundane
explanation of the GD dial is excluded in both its forms; the
model-feasibility miss is owned (fallback carried the stage,
pre-registered as such).

ELI12: We gave the gas-pressure idea its best fair shot: the
correction now grows toward each galaxy's edge, where it should
matter most. Result: applying it makes the disagreement slightly
WORSE, not better — the data bend the opposite way from what the
pressure fix predicts. So pressure isn't secretly faking the
dwarfs' answer; something about slow gas-rich dwarfs is genuinely
different. (Also owned: our fancy formula for each galaxy's gas
size had no solution — a design bug — so the simple backup rule
ran the show; the wrong-way result doesn't depend on that
detail.)

## Stage 9D EXECUTED (2026-08-05, first run, ALL GATES PASS): Q4-CARRIED 4/4 — at block grade the α localization lives in the worst-RUWE quartile; 9A's "de-fanged" reading is DOWNGRADED at its own grade

Gates: G9D-0 cube bit-identity 0.00e+00 4/4; G9D-1 free_std α̂
reproduces 9A exactly (0.446/0.522/0.472/0.420). 28.3 min.

Results: α̂ std/drop/ext = 0.446/0.000/0.442, 0.522/0.000/0.467,
0.472/0.000/0.385, 0.420/0.000/0.000 — **d(drop) = −0.42…−0.52 in
4/4 (D2 fired)**; the decensor read keeps α̂ in 3/4 but Q4's fpm
rides the NEW 4.2 edge in 4/4 (the pathological stratum's noise
demand chases any grid top; sq̂ settles at 0.3 interior once fpm
is freer). The clean 75% of the catalog (Q1–Q3), under honest
per-stratum noise, puts block-grade α̂ AT ZERO in every law-seed.
Instrument-grade caveats carried: frozen η/wr, profile+LNPI, the
drop-profile SHARPNESS not printed (lean-output miss — how many
lnL α=0 wins by is unrecorded here; the successor must print
profiles).

Consequence chain (booked): (1) 9A's concluding reading ("the
width-model threat is de-fanged at block grade") is SUPERSEDED by
this stage's D2 — the 9A measurements stand, the reading does
not; ledger bin-9a-stratalpha → CO-QUOTED with pointer. (2) The
7I-S strict-cut collapse, the 8W width address, and this stage
are ONE OBJECT — THE QUALITY-CONCENTRATION OBJECT: fitted-α
collapses under quality restriction while (as of 7I) model-light
statistics survived it. (3) The decisive next instrument is
therefore the MODEL-LIGHT quality cross-check (9E, pre-registered
below): does the fit-free median boost concentrate in Q4? If the
core is quality-blind, the fitted-α concentration reads as an
absorber-field/power effect at block grade and the full marginal
stratified re-run (now REQUIRED per the fired clause) arbitrates;
if the core is ALSO Q4-carried, the anomaly itself is in mortal
danger and the next review round's credence map handles it.

Plain verdict: NEEDS REFINEMENT (the reading, not the gates) —
the sharpest internal threat in the program is now open and
named; no credence movement this round (pre-stated), the movement
belongs to the pre-signed successors.

ELI12: Bad news round. Take away the quarter of star pairs with
the wobbliest position-tracking, and — at this quick-look grade —
the gravity-boost knob falls to zero on the clean three-quarters.
So either the boost evidence really lives in the flagged pairs
(very bad for the anomaly), or this quick instrument loses its
power without them (the full careful machine must re-decide).
Before panic: the anomaly's simplest number — the raw median
speed excess, no model at all — survived exactly this kind of cut
before. Checking THAT against the quality grades is the next,
nearly-free test, and we pre-registered it before running.

## Stage 9E pre-registration (2026-08-05, committed before the run): THE MODEL-LIGHT QUALITY CROSS-CHECK

Measurement round; no credence movement (any move belongs to the
next review round's pre-signed map). Ledger row: bin-9e-mlquality.
([calcs/stage9e_mlquality.py](calcs/stage9e_mlquality.py)) The
fit-free core statistic per RUWE quartile: B(q) = median(ṽ | wide,
q) / median(ṽ | narrow, q) with wide = s ∈ [6, 50) kAU, narrow =
s ∈ [0.2, 2) kAU (the narrow arm ≈ Newtonian internal control;
quartile-level common modes cancel to first order in the ratio);
corrected-kernel ṽ (the in-script vt_d). Report B(all), B(Q1..Q4),
Δ = B(Q4) − B(Q123 pooled), 2000-rep pair-bootstrap CIs; also the
raw wide-arm medians per quartile (no ratio) as the co-read.
Gates: G9E-0 count conservation (quartile × arm cells sum);
G9E-1 rng logged; G9E-2 the quartile edges match 8Z/9A (1.051/
1.118/1.231). Bars: E1 CORE-BLIND iff |Δ| ≤ 0.05 OR P(sign flip)
≥ 0.20 (the model-light core does not concentrate in Q4 — the 9D
collapse reads as an absorber/power effect at block grade; the
marginal stratified re-run arbitrates the fitted channel). E2
CORE-CARRIED iff Δ ≥ 0.10 AND P(Δ ≤ 0) ≤ 0.05 (the anomaly's core
is quality-carried — mortal-danger flag; next review round's map
decides). E3 GRAY-CARRIED else.

## Stage 9E EXECUTED (2026-08-05, first run, ALL GATES PASS): GRAY-CARRIED by the letter — and the core is ALIVE IN THE CLEAN QUARTILE: B(Q1) = 1.122; the 9D threat is not corroborated model-light

Gates: G9E-0/1/2 all PASS (quartile edges = 8Z/9A exact). 0.0 min
(instant statistic). A consistency dividend: the in-script
wide/narrow double-ratio lands B(all) = 1.0779 — the corrected
anchor 1.078 reproduced by construction-independent arithmetic.

Results: B(Q1..Q4) = 1.1216 / 1.0355 / 1.0700 / 1.1498; **the
CLEANEST quartile carries the boost at full strength** (second-
highest); Q4 elevation D = B(Q4) − B(Q123) = +0.077, NOT
significant (P(D ≤ 0) = 0.153; 5/50/95 = −0.043/+0.067/+0.187);
pattern non-monotone (Q2 lowest) — not a quality-dose shape. Q4's
narrow-arm median is 13% elevated (its noise floor inflates both
arms — the double-ratio design controls exactly this). Letter:
E1 missed by margins (|D| 0.077 vs 0.05; P(flip) 0.153 vs 0.20),
E2 not met ⇒ GRAY.

Standing after the 9D/9E pair: the fitted-α block-grade collapse
(9D) coexists with a model-light core alive in the clean quartile
(9E) — the split REPRODUCES the 7I-S pattern (fitted collapses,
model-light survives) now at quartile resolution. The
quality-concentration object is therefore about the FITTED
channel's absorber field, not (on current evidence) about the
anomaly's core — but the letter is GRAY, not acquittal. THE
REQUIRED ARBITER (fired by 9D's D2): the full marginal re-run
with stratified noise, profiles printed. Queued as the top binary
item; heavier than a block round (cube-grade).

Plain verdict: SUCCESS as a measurement (the core's quartile map;
the anchor reproduced); the fitted-channel question stays OPEN
pending the marginal arbiter.

ELI12: Emergency check after the scare: forget all fancy fitting
— just the raw speed excess, quality grade by quality grade. The
cleanest quarter of the sky shows the full effect (bigger than
average, even). So the anomaly isn't made of bad measurements.
The scare stage's zero was about how the fancy fitter distributes
blame when you remove its noisiest students — the big careful
re-run will settle that — but the thing itself is visible in the
best data we have.

## Stage 9F PRE-REGISTERED (2026-08-06, THE MARGINAL STRATIFIED RE-RUN — the arbiter 9D fired; DECIDER round, pre-signed credence map)

Question: does the fitted-α channel survive when the stratified-noise
model is read at MARGINAL grade (posterior weights, not profile max)?
9A/9D profiled; the max operator is exactly where a pathological
stratum (Q4 grid-top riding) can carry a spurious optimum. The
marginal integrates over the noise nuisances and pays for volume —
the honest read.

Machinery (bit-verbatim 9D): build_pop / e_of_x / vp_c / project /
eval_block_g; RUWE-quartile STRATA; (eta, wr) frozen at each
law-seed's cube MAP (block grade — disclosed; lineage-gated);
A_GRID = 0..2 step 0.5; SEEDS (31, 101) x laws (simple, BE);
N = 500k. Per α: per-stratum std tables (FPM_GRID, SQ_GRID) + Q4
extended table (FPMX <= 4.2, SQX <= 0.5); ALL tables archived to
data/stage9f_tables_{seed}_{law}.npz; the combiner is a PURE READER
of the npz (round-trip gated).

Combiner: uniform cell priors on noise axes (normalized -ln Ncells),
LNPI on shared fcomp, uniform on ffly/kw/ws and on the 5-point α
grid. Per-stratum (fpm, sq) marginalized INDEPENDENTLY (per-stratum
noise); shared axes marginalized jointly. Configs: M-STD (4 strata,
std grids), M-DROP (Q1-Q3 only), M-EXT (Q4 on extended grids),
M-TIED (single shared (fpm, sq) — the aggregation world). Report
per law-seed: full weight vector w(α), α_marg = grid mean, P(α=0),
dN = [lse(lnZ(α>0)) - ln 4] - lnZ(0); per-stratum noise posteriors
E[fpm_q], P(fpm_q = top), E[sq_q], P(sq_q = top) at M-STD (the 9D
lean-output REPAIR); d(drop) = α_marg(M-DROP) - α_marg(M-STD),
d(ext) likewise.

Gates (any FAIL => STOP, amendment pre-quote): G9F-0 cube identity
at α_MAP (bit, 4/4, verbatim G9D-0); G9F-1 profile-mode regression —
the combiner in max mode reproduces 9D's printed a_hat std/drop/ext
(12 values, bar 0.002); G9F-2 combiner analytic unit check
(synthetic table, hand-computed logsumexp, 1e-12); G9F-3 npz
round-trip bit-identity before any combiner read.

Bars (locked, ORDERED, exhaustive; rows = 4 law-seeds at M-STD with
αm, p0 = P(α=0), dN, αd = α_marg(M-DROP)):
  1. F-COLLAPSED   iff >= 3/4 rows: (αm <= 0.15 OR p0 >= 0.50).
  2. else F-SURVIVES iff >= 3/4 rows: (αm >= 0.40 AND dN >= +8
     AND αd >= αm - 0.15).
  3. else F-Q4-CARRIED iff >= 3/4 rows: (αd - αm <= -0.25).
  4. else F-GRAY-CARRIED (rows stand as measurements).

PRE-SIGNED CREDENCE MAP (decider round — movement by map ONLY):
F-COLLAPSED: anomaly-real 58 -> 48 (floor: 9E model-light core +
the (band=9, cliff=2) census + untouched galaxy legs). F-SURVIVES:
58 -> 63 (cap: AMBIGUOUS-CARRIED sky band + open width-shape
systematic). F-Q4-CARRIED: 58 -> 53. F-GRAY-CARRIED: HOLD 58.
Disclosed limits: frozen (eta, wr); 5-point α grid (weights printed
in full; the grid mean is the summary); M-EXT volume normalized by
uniform cell priors (the 7J-z6 volume lesson); 2 realization seeds.
Output: data/stage9f_stratmarg.txt; script
calcs/stage9f_stratmarg.py.

## Stage 9G PRE-REGISTERED (2026-08-06, the ambient control on the GD dial — measurement round, NO credence movement)

Question: is the GD dial-tension carried by measured ENVIRONMENT?
The 8X type-split + 8Y/9C pressure kill leave "genuine slow-dwarf
physics or an untested-class systematic"; the program's own AMB
function makes environment the loaded axis. Chae+21 Table 3
per-galaxy ambients (data/chae2021_table3.csv, 109 rows; PRIMARY =
maxclust, noclust co-read): name-match gives GD 21 / DD 70 matched
(names-only feasibility count; ambient VALUES unseen at pre-reg).

Engine: verbatim 8V lift (loader, m2ll_vert, m2ll_fast bincount
engine, lam_hat_fast on LGB -2.0..1.5 step 0.25). Reads: A
(between-type): Δmed = med(log eN | GD) - med(log eN | DD),
permutation test 20000 draws rng 7. B (within-GD tracking):
GD-matched split at its FROZEN point median log eN -> λ̂(hi) vs
λ̂(lo) point fits + paired galaxy bootstrap NBOOT = 200 rng 71
(replicates split at the frozen threshold; halves < 5 => skip+count;
point halves < 8 => read B DESCRIPTIVE-ONLY, its bar cannot fire).
C (DD control): same within DD-matched, NBOOT = 100.

Gates: G9G-0 verbatim-lift (GD/DD counts 38/111 + OFF-probe equals
8S-c printed at 1e-6); G9G-1 fast-vs-verbatim 6 probes <= 1e-6;
G9G-2 match audit (>= 60 total AND >= 12 GD matched, else ABORT
UNDERPOWERED); G9G-3 rng fingerprints.

Bars (locked, ordered): A1 ENV-SPECIAL iff |Δmed| >= 0.30 dex AND
perm P <= 0.05 (maxclust decides; noclust printed). B1 ENV-TRACKS
iff read B non-descriptive AND |D_point| >= 0.50 AND
min(P(D>=0), P(D<=0)) <= 0.05. SEVENTH-CONTROL-PASSED iff neither
fires and read B non-descriptive. GRAY-CARRIED else. Both A1 and B1
may fire (print both; verdict names all fired). NO credence
movement (pre-stated). Output: data/stage9g_gdambient.txt; script
calcs/stage9g_gdambient.py.

## Stage 9G EXECUTED (2026-08-06, all gates PASS): SEVENTH-CONTROL-PASSED — the GD dial is not carried by measured environment

Gates: G9G-0 counts 38/111 + OFF probe d = −3.05e-10; G9G-1 fast
cross-gate 3.64e-12; G9G-2 match 21 GD / 70 DD / 91 total.

Read A (between-type): GD and DD occupy the SAME measured ambient
field — Δmed(log e_N) = +0.020 dex (maxclust PRIMARY, perm
P = 0.813) / +0.040 dex (noclust, P = 0.500). The slow-gas-dwarf
population is environmentally field-typical: whatever makes GD
special, it is not where they sit in the Chae ambient map.

Read B (within-GD tracking): frozen-median split hi 10 / lo 11,
λ̂(hi) = −0.568 vs λ̂(lo) = −0.744, D = +0.176 — under the 0.50
bar; bootstrap 5/50/95 = −2.34/+0.27/+2.64, P(D ≥ 0) = 0.57 =
no ordering (replicate halves are small; hi-set replicate edge
activity 48 lo / 21 hi of 200 — disclosed, descriptive noise).

Read C (DD control, descriptive): a weak within-DD lean — hi-eN
λ̂ = +0.201 vs lo-eN λ̂ = +1.500 (the +1.5 GRID TOP; edge-riding
flagged), D = −1.30, P(D ≥ 0) = 0.10. Not bar-grade, direction =
higher-ambient DD sits less BE-ward; logged as a lead for the
ambient-gated function's galaxy leg, not pursued here.

Verdict by the letter: neither A1 nor B1 fires, read B
non-descriptive ⇒ SEVENTH-CONTROL-PASSED. The GD dial-tension has
now survived: gas budget (8S-b), vertical (8S-c), correlation
(8V), correction-lever (8Y), quality flag (9B), radius-weighted
correction sign (9C), and measured environment (9G). Remaining:
genuine slow-gas-dwarf physics, or a systematic outside every
tested class. NO credence movement (measurement round).

Plain verdict: SUCCESS as a control (the environmental escape for
the GD tension is closed in-catalog; the AMB-adjacent reading of
the type-split loses its cheapest support).

ELI12: Maybe the odd gas-rich galaxies live in odd neighborhoods,
and the neighborhood (not the galaxy) causes the odd reading? We
looked up each galaxy's measured neighborhood strength. The odd
galaxies live in perfectly ordinary neighborhoods, and among them
the odd reading doesn't care how strong the neighborhood is. So
the neighborhood excuse is dead too — seven excuses tested, seven
dead.

## Stage 9F EXECUTED (2026-08-06, all gates PASS): F-Q4-CARRIED 4/4 — THE MAP EXECUTES: anomaly-real 58 → 53

Gates: G9F-0 cube bit-identity 0.00e+00 4/4; G9F-1 profile-mode
regression reproduces 9D's printed a_hat std/drop/ext EXACTLY
12/12; G9F-2 analytic 1e-12; G9F-3 npz round-trip 4/4. Tables
archived: data/stage9f_tables_{seed}_{law}.npz (future readers
need no GPU).

The arbiter's answer, per law-seed (M-STD → M-DROP):
  simple 31 : 0.500 (dN +8.0) → 0.000 (P(α=0) = 1.000, dN −14.6)
  BE 31     : 0.505 (dN +3.8) → 0.000 (1.000, −13.4)
  simple 101: 0.500 (dN +5.0) → 0.000 (1.000, −19.7)
  BE 101    : 0.499 (dN +4.7) → 0.000 (1.000, −14.4)
With all four strata the marginal keeps α at the 0.5 grid point
with small-but-positive Newton contrast; WITHOUT the worst-RUWE
quartile the posterior is delta-at-Newton and the contrast turns
NEGATIVE — Q1-Q3 kinematics prefer Newton over the α ≥ 0.5 grid
by 13-20 lnL. M-EXT: three rows hold ≈ 0.5; BE-101 sags to 0.322
(its block ext also collapsed); Q4-extended rides fpm = 4.2 at
P = 0.90-0.97 with sq = 0.3 INTERIOR. M-TIED (aggregation world):
dN ≈ 0 — the stratified model's Newton contrast lives in the
strata structure.

Noise posteriors (the 9D lean-output repair, M-STD): Q1 E[fpm] =
1.97-2.25 (P(fpm=3.0) up to 0.37) — EVEN THE CLEANEST QUARTILE
wants ~2× formal errors; Q2/Q3 fpm 1.5-2.0 with sq 0.06-0.11; Q4
fpm → whatever top exists (3.0 std, 4.2 ext) with P(sq=0.3) = 1.00.
The width object, stratified: a floor-level ~2× hunger everywhere
+ an unbounded Q4 chase.

Letter: COLLAPSED 0/4 (M-STD holds α), SURVIVES 0/4 (every row
fails the drop condition; 3/4 also miss dN ≥ +8), Q4-CARRIED 4/4
⇒ F-Q4-CARRIED. THE PRE-SIGNED MAP EXECUTES: anomaly-real
58 → 53.

Scope/honesty: within-instrument statements only — the operative
+14.5-23.8 band is the full-cube instrument (eta/wr marginalized,
unsplit); no direct subtraction. GRID CAVEAT (named): the 5-point
α grid cannot resolve α < 0.5 — "Q1-Q3 prefer Newton" means
"reject α ≥ 0.5"; a small α (the model-light core's 9E strength
suggests roughly α ~ 0.2-0.4 territory) is UNTESTED in the drop
world. SUCCESSOR NAMED: the fine-α drop scan (α grid 0..0.6 step
0.1 on M-DROP, same machinery, ~30 min) — does clean-quartile
kinematics allow the small α the model-light core implies, or
exclude it? That question is now the entire fitted-channel case.

Plain verdict: NEEDS REFINEMENT (the honest read of a decider
that fired against us: the fitted channel is quality-carried at
marginal grade; the map moved credence down 5; the fine-α scan
and DR4 are the named next instruments — and the model-light
core's 9E result stands untouched on the other side of the
ledger).

ELI12: The big careful re-run agreed with the scare, not the
relief: when each quality group gets its own honest noise knobs,
the gravity signal in the fit is carried ENTIRELY by the
worst-measured quarter of the pairs — remove them and the fit
actually votes for plain Newton, at least against the
half-strength-or-more versions the grid can see. Meanwhile the
raw speed excess is still right there in the cleanest data
(last stage). Both can't be the full story: either the model is
missing something about how clean pairs scatter (the fit can't
see a smaller boost between its grid steps — next test), or the
anomaly is partly an artifact of bad measurements. We moved our
honesty number down five points, exactly as we promised we would
before seeing the answer.

## Stage 9H PRE-REGISTERED (2026-08-06, THE DATA-PROVENANCE MANIFEST — the author's "are you sure about the data?" answered as an instrument, not a reassurance)

SHA256 every load-bearing input + recompute the historical
invariants LIVE with verbatim loader expressions, compared to
stage-of-record printed values. Files hashed: edr3_binaries.fits.gz,
SPARC mrt, rotmod census (count + aggregate hash), chae2021_table3
.csv, efe_boost_{simple,be}_g1p2.npy, the 4 lker cubes,
stage7jz_prior.npz, the 4 stage9f table npz. Invariants: mask N =
14071 (8Z quartile sum); RUWE edges 1.051/1.118/1.231 (8Z/9A);
plxsn edges 317.311/407.883/570.585 (8Z); kept galaxies 149 with
GD/DD = 38/111 (8S-c); B(all)/B(Q1) = 1.0779/1.1216 (9E, verbatim
statistic); Chae rows 109, matched 21/70 (9G). Output:
data/stage9h_manifest.txt + data/MANIFEST.sha256 (committed — any
future input drift is one re-run away from detection). Verdict
grammar: DATA-VERIFIED iff every invariant matches its record;
DRIFT-DETECTED else (loud, per-line). NO credence movement.

## Stage 9I PRE-REGISTERED (2026-08-06, THE FINE-α DROP SCAN — the 9F successor; decider-adjacent, pre-signed map)

Question: 9F's drop world rejects α ≥ 0.5 — but the 5-pt grid
cannot see below 0.5, and the model-light core (9E) is alive in
Q1. Does Q1-Q3 kinematics ALLOW a small boost, or exclude it?

Machinery bit-verbatim 9F; A_FINE = 0.0..0.6 step 0.1 (7 pts);
per α: 4 stratum STD tables (no extended grids this stage);
(eta, wr) frozen at MAP; tables archived
data/stage9i_tables_{seed}_{law}.npz. PRIMARY = M-DROP fine
marginal (Q1-Q3, per-stratum noise, LNPI on shared fcomp);
co-reads: M-STD fine, per-quartile single-stratum posteriors
(diagnostic), profile-mode a_hat, and the full lnZ(α) − lnZ(0)
curve per law-seed (ship the risk axis).

POWER GATE (the 7J fullpow standard — a null letter without
demonstrated sensitivity is not evidence): per law-seed, draw ONE
synthetic Q1-Q3 sky at truth α = 0.3, fcomp = 0.20, ffly = 0.05,
kw = 1.0, ws = 0, per-stratum (fpm, sq) = Q1 (2.1, 0.0), Q2
(1.8, 0.1), Q3 (1.8, 0.1) [nearest-grid to the 9F posteriors];
multinomial per (stratum, s-bin) at the real ND counts, rng
default_rng(9); run the identical M-DROP marginal on it. G9I-P
PASS iff recovered α_marg >= 0.15 AND P(α=0) <= 0.35.
Single-injection grade, disclosed.

Gates: G9I-0 lineage BIT-identity — STD tables at shared α ∈
{0.0, 0.5} equal the 9F npz exactly; G9I-1 analytic combiner
check (1e-12); G9I-2 synthetic census (counts == ND); G9I-3 cube
identity at α = 0.5 (unsplit block vs stored cube, bit).

Bars (locked, ORDERED; rows = 4 law-seeds, M-DROP fine: αm,
p0 = P(α=0), P01 = P(α <= 0.1)):
  1. I-POWER-FAIL iff G9I-P fails in >= 2/4 law-seeds — null
     letters BLOCKED; verdict POWER-FAIL-CARRIED (measurements
     quoted).
  2. else I-SMALL-ALPHA iff >= 3/4 rows: αm ∈ [0.10, 0.45] AND
     p0 <= 0.35.
  3. else I-NEWTON-FLAT iff >= 3/4 rows: P01 >= 0.60 AND that
     row's own G9I-P passed.
  4. else I-GRAY-CARRIED.
PRE-SIGNED MAP: I-SMALL-ALPHA → anomaly-real 53 → 56 (channels
reconcile at deflated amplitude); I-NEWTON-FLAT → 53 → 50 (clean
kinematics exclude even small boosts at fitted grade; model-light
core still stands); POWER-FAIL / GRAY → HOLD 53.
Output: data/stage9i_finealpha.txt; script
calcs/stage9i_finealpha.py.

## Stage 9H EXECUTED (2026-08-06, amended run): DATA-VERIFIED 15/15 — the author's data question answered as an instrument

SHA256 manifest of 15 load-bearing entries (incl. the 175-file
rotmod aggregate) written to data/MANIFEST.sha256 (committed); all
15 live-recomputed invariants match their stage-of-record printed
values: mask N = 14071; RUWE edges 1.051/1.118/1.231; plxsn edges
317.311/407.883/570.585; fit galaxies 149 (GD/DD 38/111); B(all) =
1.0779, B(Q1) = 1.1216; Chae 109/21/70. AMENDMENT owned (pre-quote,
logged in-script): first firing compared the rotmod FILE count
passing meta cuts (153) to the record 149 — but the program's 149
is galaxies with surviving DATA POINTS (4 galaxies lose every point
to per-line quality cuts and never enter a fit); the gate caught
the manifest's own wiring, which is the gate philosophy working.
Any future input drift is one re-run away from detection.

Plain verdict: SUCCESS (the data chain is now hash-pinned AND
invariant-pinned; the answer to "are you sure about the data" is
an auditable yes).

ELI12: We took a fingerprint of every data file we use and
re-derived every headline count from scratch, comparing against
what our old logs printed. All 15 match. If any file ever changes
under us, one command will catch it.

## Stage 9I EXECUTED (2026-08-06, all gates PASS incl. power 4/4): I-SMALL-ALPHA — THE MAP EXECUTES: anomaly-real 53 → 56

Gates: G9I-0 lineage BIT-identity to the 9F npz at shared α (8/8
zeros); G9I-1 analytic; G9I-2 synthetic census 4/4; G9I-3 cube
identity 4/4. POWER 4/4: injected α = 0.3 recovered at 0.300/
0.299/0.220/0.298 with P(α=0) ≤ 0.003 — the drop world's null
direction had demonstrated sensitivity.

M-DROP fine (Q1-Q3, the clean world), per law-seed:
  simple 31 : peak α = 0.1 (w 0.958), α_marg 0.104, ΔlnZ(0.1−0) = +6.5
  BE 31     : peak α = 0.3 (w 0.789), α_marg 0.274, +4.1 over Newton
  simple 101: peak α = 0.2 (w 0.475), α_marg 0.186, +1.5 over Newton
  BE 101    : Newton-flat (w(0) 0.404, spread to 0.4), α_marg 0.138
The clean kinematics are not anti-boost — they are anti-BIG-boost:
interior small-α optimum in 3/4 with mild preference (+1.5..+6.5
over Newton), and the α ≥ 0.5 rejection (−13..−25 at the 0.5/0.6
points) is exactly 9F's finding restated. 9F's scoped claim
("reject α ≥ 0.5") stands; its Newton-preference phrasing is
REFINED by the finer instrument.

Co-reads: M-STD fine α_marg 0.33-0.48 at dN +10..+14 — but
simple-31 puts 0.69 at the NEW grid top 0.6 (edge-riding,
disclosed; the fine grid should have extended further — co-read
limited, primary unaffected since M-DROP dies by −13 before 0.5).
Q-alone posteriors (diagnostic): Q1 0.12-0.26, Q2 0.05-0.11
(most Newton-ish — matches 9E's B(Q2) lowest), Q3 0.31-0.45, Q4
0.49-0.59 riding the 0.6 top 3/4 (the chase again). Fitted
quartile ordering Q2 < Q1 < Q3 < Q4 vs model-light Q2 < Q3 < Q1
< Q4 — Q1 SWAPS (cleanest data: fitted small, model-light high)
= the quality-concentration object's remaining core, named.

Letter: SMALL-ALPHA rows 3/4 (BE-101 fails on P(0) = 0.404),
NEWTON-FLAT rows 2/4 — the pre-registered ORDER decides
(SMALL-ALPHA evaluated first). GRAMMAR WART owned: a row peaked
at α = 0.1 satisfies P(α ≤ 0.1) ≥ 0.6 by construction (simple-31
counted for both bars); the order made the verdict well-defined,
but the NEWTON-FLAT letter should have excluded interior-peak
rows — logged for successor design.

THE MAP EXECUTES: anomaly-real 53 → 56. The arc 58 → 53 → 56
nets: the fitted channel survived its arbiter DEFLATED (clean-
world amplitude 0.1-0.3 at stratified block grade, vs 0.42-0.52
with Q4 in) and quality-ordered; reconciliation-grade, NOT
detection-grade (per-row leans +1.5..+6.5). Queue: fine-STD
extension past 0.6 (the edge co-read); the two-channel Q1
contrast instrument; DR4 = the external arbiter.

Plain verdict: SUCCESS as a measurement (the arbiter's question
answered: small boosts are allowed and mildly preferred in the
clean world; the big-boost reading was the artifact). The
fitted-vs-model-light tension narrows to one named residual: Q1's
channels disagree.

ELI12: Last round's scare said "remove the messy quarter and the
fit wants zero gravity boost." But the fit could only choose
between zero and half-strength — nothing in between. We gave it
the in-between dial. Answer: the clean data actually LIKE a small
boost (a tenth to a third of full strength), just not the big one
the messy quarter was pushing. And when we secretly injected a
small boost into fake data, the machine found it every time — so
the answer isn't blindness. The two ways of looking (raw medians
vs careful fit) now roughly agree: something small and real, plus
a messy quarter that exaggerates. Honesty number back up a bit,
exactly per the promise we signed before running.

## Stage 9J PRE-REGISTERED (2026-08-06, THE FINE-STD EXTENSION — the 9I co-read repair; measurement round, NO credence movement)

9I's M-STD fine co-read edge-rode its 0.6 grid top (simple-31
w(0.6) = 0.69; Q4-alone at top 3/4).  Extend: A_EXT = 0.7..1.2
step 0.1 (6 new α values, identical machinery); REUSE the 9I
tables for 0..0.6 (bit-lineage from data/stage9i_tables_*.npz);
combined 13-point grid 0..1.2.  Reads: M-STD full-grid marginal
(localize?), M-DROP full-grid (confirm unchanged), Q4-alone
full-grid (where does the chase stop?).  New tables archived
data/stage9j_tables_{seed}_{law}.npz.

Gates: G9J-0 9I-npz reload + new-npz round-trip bit-identity;
G9J-1 analytic combiner (1e-12); G9J-2 lineage regression — on
the 0..0.6 SUBSET the combiner reproduces 9I's printed
M-STD/M-DROP α_marg to 0.002 (8 values).

Bars (locked, ordered; rows = 4 law-seeds, M-STD full-grid):
  1. J-INTERIOR iff >= 3/4 rows: α_marg interior AND
     P(α = 1.2 top) <= 0.10.
  2. else J-STILL-EDGE iff >= 3/4 rows: P(top) >= 0.30 (the
     chase has no ceiling on this axis — width masquerading as
     boost).
  3. else J-GRAY-CARRIED.
NO credence movement (measurement round; the 9I primary already
decided by map).  Output: data/stage9j_stdext.txt; script
calcs/stage9j_stdext.py.

CONCURRENT (same round, not a stage): REVIEW ROUND 13 — the Opus
reviewer (an Opus-model agent, per the author's 2026-08-06
correction: I spawn and run him myself; no relay) is briefed via
OPUS-NOTE.md (uncommitted) on the full quality arc + the four
asks; read-only; his review books as ROUND 13 with point-by-point
adopt/rebut.

## REVIEW ROUND 13 (2026-08-06, the Opus agent's return — full text in REVIEW-ROUND13-OPUS.md, uncommitted): the quality arc audited; math CLEAN; one framing CORRECTION; the fpm-ceiling risk axis demanded

His verdict on process: bars/maps/warts match scripts exactly;
"I could not make the combiner lie about what's in the tables";
combiner math SOUND (no invalid step, no verdict-biasing leak;
M-EXT volume handling = the 7J-z6 lesson applied). ADOPTED
point-by-point:

**E2 = CORRECTION (to be numbered at the next paper thaw): the 9E
double-ratio does NOT control differential noise.** σ(ṽ) rises
with s (same physical velocity error over a shrinking Newtonian
scale), so under Newton+underestimated errors the wide arm's
median inflates MORE than the narrow arm's — B rises above 1 with
the boost's sign. Proof from 9E's own medians: narrow arm ×1.126
Q1→Q4 (a pure non-boost pedestal on the Newtonian arm), wide arm
×1.155. The ledger note "noise floor inflates both arms — the
double-ratio controls it" is WRONG AS STATED (controls common-mode
only); B(Q1) = 1.122 = the un-decomposed sum (boost + differential
noise + pedestal), NOT independent corroboration of a boost. The
9E VALUES stand; the "counterweight/survives" READING is
withdrawn. Consequence: the two channels are not in tension — the
forward fit measures boost after subtracting modeled differential
noise; the gap between 1.122 and Q1-α̂ ≈ 0.15 is the noise term
the fit removes and B keeps. Ledger bin-9e + bin-9i notes amended.

**C ADOPTED — the operative headline is the two-sided upper-limit
form:** on clean strata α ≥ 0.5 is excluded (−13..−25); α is
bounded above at a deflated value (≈0.3-0.5 under the
fpm-ceiling/fcomp-conversion treatments); small-real-boost vs
Newton-plus-width UNRESOLVED. "Reconcile at 0.1-0.3" was
over-claimed as a positive statement (his D-iii; also his A2:
the single self-consistent injection showed −0.08 recovery bias
in the ideal case → the point estimates are biased LOW; noted,
not applied, single-injection grade).

**A2 LNPI audit adopted:** effective fcomp weights [0, .214,
.365, .418, .004, 0] — zero mass at fcomp = 0 (data-consistent:
7J rejected fcomp=0 at −660..−811, but the ABSOLUTE α inherits
the conversion; the profile-max over the band is the permissive
direction). Conversion-band sensitivity row QUEUED.

**Credence: FROZEN at 56 pending D1/D2(/D3)** (cadence rule; his
D-ii "HOLD at least as defensible" on the record; the 9I +3
stands by pre-signed process — no retroactive un-moves; the next
pre-signed map carries any movement).

**Instruments adopted into the queue:** D1 = the fpm-ceiling
curve (9K, THIS ROUND — reader-only); D2 = the narrow-pair fpm
meter (next binary item; the mechanism-level adjudicator of
B-vs-fit); D3 = the degraded (misspecified) injection
(conditional on D1/D2); A4 = the rotation-curve
convergence/V_flat-definition split on the GD dial (next galaxy
item — a candidate eighth control targeting the V-ordering
itself; his premise matches 8Y's VMET fallback implementation).
E1 confirmed as the disclosed 9I grammar wart (fix = P(α=0)
letters in successors).

## Stage 9K PRE-REGISTERED (2026-08-06, THE FPM-CEILING CURVE — round-13 D1; reader-only; measurement round, NO credence movement)

From the archived data/stage9i_tables_*.npz (no GPU): re-marginalize
M-DROP fine with the fpm axis truncated at caps = grid prefixes
{1.2, 1.5, 1.8, 2.4, 3.0-full} (Lindegren ≈1.4 → nearest grid
1.5); per law-seed report α_marg(cap), P(α=0)(cap), and the
max-mode peak lnL(cap). Co-reads: M-STD capped likewise;
per-stratum ΔlnZ(≤1.5 vs full) — who pays the cap. Gates: G9K-0
lineage — at the full cap the combiner reproduces 9I's printed
M-DROP/M-STD α_marg (8 values, 0.002); G9K-1 analytic (1e-12);
G9K-2 monotonicity — max-mode peak lnL non-increasing as the cap
tightens (a violation = reader bug).

Bars (locked, ordered; rows = 4 law-seeds, M-DROP):
  1. K-ROBUST iff >= 3/4 rows: |α_marg(≤1.5) − α_marg(full)| < 0.05.
  2. else K-CAP-FRAGILE iff >= 3/4 rows: α_marg(≤1.5) −
     α_marg(full) >= +0.15 (the deflation was noise-absorption;
     direction-1 pressure — the boost was being eaten by fpm).
  3. else K-GRAY-CARRIED.
NO credence movement (pre-stated; the movement map belongs to the
post-D1/D2 decider). Output: data/stage9k_fpmcap.txt; script
calcs/stage9k_fpmcap.py.

## Stage 9J EXECUTED (2026-08-06, all gates PASS): J-INTERIOR 4/4 — the extended axis localizes everything; the α ladder is a RUWE dose curve

Gates: G9J-0 9I-reload + round-trip 4/4; G9J-1 analytic; G9J-2
subset regression to 9I printed EXACT 8/8 (M-DROP unchanged to
the digit — pure-extension lineage).

M-STD full (13-pt grid 0..1.2): α_marg = 0.785/0.723/0.385/0.399
(peaks 0.8/0.9/0.4/0.4), P(top) = 0.000 in 4/4, dN = +14.3/
+12.8/+13.2/+9.4 — the 9I co-read edge was a grid artifact; the
full-strata numbers are quotable. SEED-SPLIT noted: seed 31 wants
0.72-0.79, seed 101 wants 0.38-0.40 — realization scatter at the
α grade (the 3A object), disclosed.

THE LADDER (Q-alone, full grid): Q2 = 0.05-0.11 < Q1 = 0.12-0.26
< Q3 = 0.40-0.79 < Q4 = 0.90-1.05 (Q4 peak 0.9-1.2; one seed
still presses the top at P = 0.45). The quality-concentration
object is now a monotone-in-RUWE DOSE CURVE (Q1/Q2 swap noted vs
naive ordering): fitted α rises with astrometric dirtiness —
under round-13's frame, exactly the α↔fpm shape-split expected
if dirtier strata carry more α-shaped (s-growing) excess; under
the boost frame, the clean strata bound the physics and the
dirty strata exaggerate. Same numbers, both readings — D2/D3
decide.

Plain verdict: SUCCESS as a measurement (the axis repair landed;
the co-read is un-flagged; the ladder is the object's cleanest
portrait yet).

ELI12: We gave the fit a longer dial and the "wants more than
the dial goes" problem vanished — every group now picks a real
setting. The picture: the cleaner the measurements, the smaller
the gravity knob they choose, in a smooth staircase from ~0.1
(cleanest) to ~1.0 (messiest). Either noise dresses up as
gravity in messy data, or gravity is there and messy data
overstate it. Next instrument decides.

## Stage 9K EXECUTED (2026-08-06, round-13 D1, all gates PASS; one pre-verdict wiring amendment owned): K-GRAY by the letter — the cap direction is REAL at half the fragile bar

Reader-only from the 9I npz. Amendment (pre-quote): a 5-dim
LNPI reshape against the 6-dim sliced table crashed the
per-stratum co-read on first firing (after gates, before any
verdict); fixed, re-run clean.

The curve (M-DROP α_marg at fpm caps 1.2/1.5/1.8/2.4/3.0):
  simple 31 : 0.211/0.155/0.137/0.114/0.104   d(1.5) = +0.051
  BE 31     : 0.340/0.355/0.307/0.296/0.274   d(1.5) = +0.081
  simple 101: 0.274/0.231/0.217/0.202/0.186   d(1.5) = +0.045
  BE 101    : 0.128/0.195/0.197/0.220/0.138   d(1.5) = +0.057
ALL FOUR positive at the Lindegren-adjacent cap — the reviewer's
direction-1 pressure (fpm eats boost) is measured and
sign-consistent, at ~half his +0.15 fragile bar ⇒ K-GRAY (1
ROBUST row, 0 FRAGILE rows). THE CO-FACT: under cap 1.5 every
law-seed prefers small α over Newton — P(α=0) = 0.080/0.006/
0.000/0.026 (BE-101's Newton-flatness, 0.404 at full freedom,
COLLAPSES to 0.026: the drop world's Newton-compatibility lives
in the fpm > 1.5 cells). Evidence cost of the cap: peak lnL
1.3-10.2 (full freedom buys ≤10.2); per-stratum-alone ≤3.4
(Q1 pays most, simple-31). G9K-2 cap-monotonicity PASS 4/4.

Standing: the deflated drop-world α is mildly noise-absorption-
suppressed (+0.05-0.08 recoverable at a physical ceiling) but
does NOT jump to the 0.3+ that would have fired CAP-FRAGILE.
The D2 narrow-pair fpm meter (which decides whether fpm ≈ 2 is
legitimate) is now the single decisive instrument, exactly as
round 13 ranked it. NO credence movement (pre-stated; FROZEN 56
per the cadence rule).

Plain verdict: SUCCESS as a measurement (the missing risk axis
shipped same-day; both directions bounded).

ELI12: The referee asked: "if you forbid the fit from calling
more than 1.5× the official error bars 'noise', does the
gravity knob jump up?" Answer: it rises a little in all four
tests (so some gravity WAS being filed under noise), but not a
lot — and with that noise leash on, even the one holdout test
now prefers a small boost to plain Newton. The tiebreaker test
(measuring the true noise on close pairs where gravity can't
matter) is queued next.

## Stage 9L PRE-REGISTERED (2026-08-06, THE NARROW-PAIR FPM METER — round-13 D2, THE DECIDER the freeze waits on; pre-signed credence map)

Question: is Q1's fpm ≈ 2 legitimate noise or eaten boost? The
narrow s-bins (0.2-2, 2-6 kAU) are near-boost-free — whatever fpm
they demand is honest noise; if they demand less than the joint
fit uses, the wide bins are inflating fpm with relabeled boost.

Machinery: 9F loader/machinery bit-verbatim; α FIXED at 0 (at α=0
tab_a ≡ 1 for BOTH laws — the meter is law-blind; rows = 2 seeds);
eval restricted to s-bins 0-1 (narrow-only tables per stratum,
axes (fc, fy, fpm, kw, sq, ws)); per-stratum (fpm, sq) posterior
with the usual uniform cell priors + LNPI on shared fc; the FULL
fpm marginal vector shipped per stratum (risk-axis rule). Tables
archived data/stage9l_tables_{seed}.npz.

Gates: G9L-1 analytic combiner (1e-12); G9L-2 the boost-free
premise MEASURED — model narrow-bin ṽ-median ratio (α=1 vs α=0)
≤ 1.05 per s-bin (needs one α=1 orbit run per seed·law — 4 runs,
co-computed); G9L-3 injected-fpm recovery — synthetic narrow sky
at the Q1 truth cell (fpm = 2.1, sq = 0.0, fc = 0.20, kw = 1.0,
ws = 0; multinomial at real narrow counts, rng 9): recovered
E[fpm(Q1)] within ±0.25 of 2.1; G9L-4 counts census. (No lineage
anchor exists for a first-of-kind eval — the injection gate
carries that burden, disclosed.)

Bars (locked, ordered; the precondition — 9F joint E[fpm(Q1)] ≥
1.9 — holds at 1.97/1.99/2.22/2.25):
  1. L-NOISE-REAL iff BOTH seeds: E[fpm_narrow(Q1)] ∈ [1.8, 2.3].
  2. else L-BOOST-EATEN iff BOTH seeds: E[fpm_narrow(Q1)] ≤ 1.5.
  3. else L-GRAY-CARRIED.
AMENDMENT 1 (2026-08-06, post-G9L-2-FAIL, pre-quote, verdict
withheld by the stop rule on run 1): the registered 0.2-6 kAU
window FAILS its own boost-premise gate at the 2-6 bin (measured
α=1/α=0 ṽ-median ratios 1.052-1.080 > the 1.05 bar; the 0.2-2
bin passes at 1.003-1.008) — boost contamination there would
bias the meter TOWARD noise-real. The meter window narrows to
the 0.2-2 kAU bin ALONE; bars, gates, map unchanged.
PRE-SIGNED MAP (the frozen-56 decider): L-NOISE-REAL →
anomaly-real 56 → 53 (the fit's noise subtraction is legitimate;
the upper-limit reading stands with honest noise; the model-light
pedestal is real). L-BOOST-EATEN → 56 → 60 (α under-reported;
the D3 degraded injection becomes MANDATORY next). L-GRAY → HOLD
56. Co-reads: Q2-Q4 narrow fpm posteriors; sq_narrow per stratum.
Output: data/stage9l_fpmmeter.txt; script
calcs/stage9l_fpmmeter.py.

## Stage 9M PRE-REGISTERED (2026-08-06, THE CONVERGENCE SPLIT on the GD dial — round-13 A4, candidate EIGHTH control; measurement round, NO credence movement)

Round-13's hypothesis: slow gas dwarfs disproportionately have
still-rising outer rotation curves; a rising curve read at its
last point under-states the asymptote → outer points sit low →
the dial pulled negative — a type-correlated systematic that
would masquerade as the 8Y V-ordering.

Design: verbatim 8V engine lift; per GD galaxy the outer slope
s_out = (V_N − V_{N−2})/V_N over the LAST 3 kept rotmod points
(kept = the standard quality cuts); flags FROZEN before any fit:
RISING iff s_out ≥ 0.05; CONVERGED iff s_out ≤ 0.02; AMBIG else
(excluded from the split, counted). λ̂(GD-conv) vs λ̂(GD-rising)
point fits + paired galaxy bootstrap NBOOT = 200 rng 71 (draw
from GD-38, split by frozen flags; halves < 5 → skip+count).
SMALL-N rule: point subset < 8 → DESCRIPTIVE-ONLY, bars cannot
fire. Co-reads: DD control split (point fits only); the flag
census; edge fractions.

Bars (locked, ordered; D = λ̂(rising) − λ̂(conv)):
  1. M-ARTIFACT iff non-descriptive AND D ≤ −0.5 with
     P(D ≥ 0) ≤ 0.05 AND λ̂(conv) ≥ −0.6 (the negative dial
     concentrates in rising curves; the converged half releases
     toward the dial — the V-ordering is a V_flat-definition
     artifact).
  2. else M-HARDENED iff non-descriptive AND λ̂(conv) ≤ −0.8 AND
     the conv bootstrap 95th pct ≤ −0.3 (the converged subset
     alone carries the dial robustly — the EIGHTH control
     passes; physics/selection reading hardens).
  3. else M-GRAY-CARRIED (incl. all descriptive cases).
NO credence movement (control round). Output:
data/stage9m_convsplit.txt; script calcs/stage9m_convsplit.py.

## Stage 9L EXECUTED (2026-08-06, amended run, all gates PASS): L-NOISE-REAL — THE DECIDER FIRES: anomaly-real 56 → 53

Run 1 story: G9L-2 caught the registered window (the 2-6 kAU bin
carries a 5-8% model boost > the 1.05 bar — contamination in the
dangerous, toward-noise-real direction); verdict WITHHELD by the
stop rule; AMENDMENT 1 (pre-quote, committed e63a3f3) narrowed
the meter to the 0.2-2 kAU bin alone (premise 1.003-1.008).

Amended run: G9L-2 4/4 PASS; G9L-3 injection recovered 2.12/2.06
vs truth 2.1; counts = the 9E narrow censuses exactly. THE
MEASUREMENT: **E[fpm_narrow(Q1)] = 2.21 / 2.13 — both seeds
inside the [1.8, 2.3] noise-real band** (joint-fit precondition
1.97-2.25); η = 1.3 co-read 2.28/2.21 (insensitive — risk axis
shipped); and the meter is ~FLAT across quality: Q1-Q4 narrow
E[fpm] = 2.12-2.30 in both seeds (sq_narrow 0.03-0.13). The
boost-free arm demands the same ~2× the joint fit uses ⇒ the
fit's noise subtraction is LEGITIMATE; the round-13 fork
resolves in the direction he called mechanically forced: B's
pedestal is real differential noise, the drop-world small-α /
upper-limit reading stands, and 9K's cap-1.5 world was the wrong
world (the true multiplicative error is ≈2, above Lindegren's
single-star ceiling — consistent with his point that RUWE
certifies 5-parameter single-star fits, not 2-body velocity
budgets). STANDING MEASUREMENT logged: EDR3 wide-pair velocity
errors ≈ 2× formal at 0.2-2 kAU, quality-quartile-independent.

THE MAP EXECUTES: anomaly-real 56 → 53 (the freeze releases;
trail 50 → 58 → 53 → 56 → 53 — the reviewer's contested +3 at 9I
is unwound by his own instrument, which is the process working).
Consequences: D3 (degraded injection) DEMOTED to optional (the
fpm amplitude is now measured-honest; the width-SHAPE question
on the wide arm remains, DR4-facing); the two-channel Q1
residual is RESOLVED (E2 + 9L: B keeps the noise term the fit
removes — no contradiction remains).

Plain verdict: SUCCESS — the decider decided. The binary
forward-fit leg's operative statement: α ≥ 0.5 excluded on clean
strata; α ≲ 0.3-0.5 upper limit with a mild interior preference
(0.1-0.3) under measured-legitimate ≈2× noise; α = 0 not
excluded (one seed-law Newton-flat). The anomaly case now rests
more heavily on the galaxy legs + the census pair.

ELI12: The tiebreaker: close star pairs, where extra gravity
can't matter, need the SAME "double the error bars" as everything
else. So the doubled error bars are real measurement noise, not
hidden gravity — Gaia's official errors are just optimistic for
pairs. That means our careful fit was right to subtract that
noise, the raw-median excess really does carry a noise pedestal,
and the honest star-pair verdict is "a small boost allowed, a
big one ruled out." We moved the honesty number back down 3, by
the rule we signed before running.

## Stage 9M EXECUTED (2026-08-06, all gates PASS): M-GRAY by the letter — the V_flat-definition artifact is POINT-DEAD (wrong sign), the eighth control effectively holds

Flags (frozen): GD-38 → CONVERGED 17 / RISING 11 / AMBIG 10.
Point fits: λ̂(CONV) = −1.395 vs λ̂(RISING) = −1.172, D = +0.223
— the ARTIFACT direction needed rising ≤ conv − 0.5; measured
sign is BACKWARDS (converged is MORE negative), and the
converged subset ALONE reproduces the full GD dial (−1.395 vs
−1.31). Bootstrap: P(D ≥ 0) = 0.62 (no ordering); the M-HARDENED
letter missed only on the small-sample bootstrap spread (conv
95th pct +0.32 vs the −0.3 bar; 17 galaxies) ⇒ GRAY-CARRIED.
DD control (descriptive): the convergence axis DOES pull
negative in DD (conv +1.09 vs rising −0.48, D = −1.57, 16
rising) — the round-13 mechanism exists in the disk-dominated
sample but is NOT what carries the GD dial. Logged for the
reviewer.

Plain verdict: SUCCESS as a control at point grade (the
artifact hypothesis killed by sign; formal hardening blocked
only by subset size). The GD tension keeps its shape: genuine
slow-gas-dwarf physics or a systematic outside eight tested
classes; external resolved kinematics = the clean arbiter.

ELI12: The referee's last galaxy excuse: maybe the odd reading
comes from galaxies whose rotation curves are still climbing at
the last measured point (so we underestimate their top speed).
Checked: the odd reading is actually STRONGEST in the galaxies
whose curves have flattened — the excuse points the wrong way.
Eighth excuse dead in practice (formally "gray" only because 17
galaxies make noisy error bars).

## ROUND 13 ADDENDUM (2026-08-06, the reviewer's reaction to his own instruments' results — same agent, context retained; verified against outputs on disk): L-NOISE-REAL CERTIFIED

(1) The 9L amendment certified as the CONSERVATIVE move (stripping
the boost-contaminated bin made the meter harder to pass; it still
landed noise-real with clean premise and unbiased injection). His
direction-1 worry ("a Lindegren cap would blow α up") formally
RETIRED by 9K+9L. WORDING ADOPTED: Q1's fpm marginal carries 0.31
mass at the 3.0 grid ceiling ⇒ E[fpm] = 2.21 is mildly
ceiling-censored — quote the standing measurement as **"pair
velocity errors ≳ 2.1× formal" (a LOWER BOUND)**, which only
strengthens noise-real; extend the grid only if the central value
is ever needed. Framing kept: 9L pins the WIDTH at ~2×; whether
that width is error-underestimation or non-Gaussian shape is the
DR4 chase; the direction conclusion (the wide excess is width,
not hidden boost) holds either way.
(2) D3 optional CONFIRMED — and BANKED AS THE PROMOTION GATE
(standing rule): before any paper language stronger than "α ≲ 0.3
upper limit, α ≥ 0.5 excluded, wide excess is width", the D3
no-manufacture arm (truth-Newton sky under a wrong error-SHAPE at
the measured 2× scale must NOT manufacture the small-α lean) must
fire. Until then the operative claim stays at upper-limit /
direction-settled grade.
(3) The DD convergence lean: FOOTNOTE in the GD writeup (artifact
point-dead there by sign) + a SHORT CHECK queued (not a full
stage): does the headline deep-regime dial / c₁ zero-point move
when rising-curve galaxies are flagged? (Type-correlated
deep-regime systematic — feeds the function program.)
NO credence movement (certification of the already-executed map).

## Stage 9N+9O PRE-REGISTRATION (2026-08-06, the paper-prep robustness pair; committed BEFORE any run; measurement rounds — NO credence movement, pre-stated)

Context: the author has called the paper thaw. Two queued cheap rows
must land first so the papers quote them: the round-13-addendum
rising-flag short check (galaxy) and the round-13 LNPI
conversion-band row (binary). Both are readers on frozen machinery.

**9N — THE RISING-FLAG DIAL EXPOSURE (galaxy;
calcs/stage9n_risingdial.py → data/stage9n_risingdial.txt).**
Question: does the HEADLINE deep-regime dial (8S-c FULL vertON
lam_hat = 0.960 ⇒ ĉ₁ = 0.480, D1 0.848–1.074) move when
rising-curve galaxies (the 9M frozen flag, s_out ≥ 0.05 on the last
3 kept points) are excluded? Machinery: the 9M verbatim-lift engine
(m2ll_fast + lam_hat_fast, LGB −2.0..1.5 step 0.25). Frozen flag
rule applied to ALL 149 fit galaxies. Instruments: point dials
lam_hat(FULL), lam_hat(FULL−RISING), lam_hat(CONV-only);
GD/DD rows (full vs −RISING) + the GD−DD gap co-read (no letter —
the GD tension has its own arc); paired galaxy bootstrap
NBOOT = 200 rng 83 on D = lam_hat(rep−RISING) − lam_hat(rep), skip
rule: either side < 30 galaxies. Gates: G9N-0 verbatim-lift (GD/DD
counts 38/111 + OFF probe vs 8S-c printed −647.874131278 at 1e-6);
G9N-1 fast-vs-verbatim 3 probes ≤ 1e-6 (rng 7 sequence, 9M
verbatim); G9N-2 lineage regressions — same-engine 9M printed
values at ≤ 0.002 (GD CONV −1.395 / RISING −1.172, DD CONV +1.087
/ RISING −0.478, censuses 17/11/10 and 69/16) and cross-instrument
8S-c printed dials at ≤ 0.05 (FULL 0.960, GD −1.315, DD 1.331 —
looser bar DISCLOSED: same objective + grid, independent optimizer
path). Bars (locked, ORDERED, exhaustive): **N-MOVED** iff
|D_point(FULL)| ≥ 0.25 (one grid step) AND the bootstrap 5–95% of
D excludes 0; **N-ROBUST** iff |D_point(FULL)| ≤ 0.125 AND
bootstrap P(|D| ≥ 0.25) ≤ 0.10; **N-GRAY-CARRIED** otherwise (rows
stand as measurements). Consequence pre-stated: N-MOVED ⇒ the
Paper-2 dial quote gains a rising-flag conditional column (no
credence content); N-ROBUST ⇒ one sentence + this row's pointer.

**9O — THE LNPI CONVERSION-BAND ROW (binary;
calcs/stage9o_lnpiband.py → data/stage9o_lnpiband.txt).** Question:
do the round-13 operative letters depend on the companion-fraction
prior's conversion treatment? Pure reader over the archived
stage9i/9j npz tables (13-pt fine grid, M-DROP = Q1–Q3 primary,
M-STD co-read). Prior variants: OPER (the shipped conversion-band
ENVELOPE, rebuilt and bit-compared to the stored npz LNPI), G-LOW
/ G-MID / G-HIGH (conversion pinned at the band's low edge / mean
/ high edge), FLAT (uniform over all 6 fcomp cells incl. 0.0 — the
strongest stress: it revives the companion-free world the envelope
never priced). Statistics per variant × law × seed: a_marg,
P(a=0), P(a≤0.1), P(a≥0.5), EX05 = lnZ_drop(0.5) − max lnZ_drop,
dN_fine (7-pt sub-grid, regression-anchored). Gates: G9O-0 rebuilt
OPER LNPI bit-equal to stored npz (9I + 9J, 4/4); G9O-1 analytic
lse unit check 1e-12; G9O-2 regression under OPER to 9J printed
M-STD/M-DROP full a_marg + M-DROP P(a=0) (12 values, 0.002) and to
the 9I printed M-DROP curve at a = 0.5 (4 values, 0.06 — 1-decimal
print). Letters under test (round-13 operative): L1 "α ≥ 0.5
excluded on clean strata" ⇔ EX05 ≤ −8; L2 "no Newton verdict
either way" ⇔ −8 ≤ dN_fine ≤ +8. Bars (locked, ORDERED,
exhaustive): **O-FRAGILE** iff any variant breaks L1 in ≥ 2 of 4
law×seed rows OR breaks L2 in ≥ 2 rows; **O-ROBUST** iff L1 AND L2
hold in all 4 rows under ALL variants; **O-GRAY-CARRIED**
otherwise. Consequence pre-stated: O-ROBUST ⇒ the paper's
upper-limit language carries a prior-robustness pointer; O-FRAGILE
⇒ the α band gains an explicit prior-conditional column and the
paper quotes the spread.

Both stages: NO credence movement regardless of letter
(measurement rounds; the operative anomaly-real 53 is untouched by
pre-statement). Wall-clock: 9O seconds; 9N bootstrap ~0.5–2 h CPU.

## Stage 9N+9O EXECUTED (2026-08-06, same day): N-ROBUST + O-FRAGILE-MAPPED — the last pre-paper rows land

**9N N-ROBUST (all gates first-run: G9N-0 counts + OFF probe
3.05e-10; G9N-1 3.64e-12; G9N-2a censuses exact; G9N-2b
same-engine 4/4 ≤ 0.001; G9N-2c cross-instrument FULL +0.969 vs
8S-c 0.960 / GD −1.309 vs −1.315 / DD +1.337 vs 1.331 — the fast
engine reproduces the headline dial).** The deep-regime dial does
NOT depend on rising-curve galaxies: lam_hat(FULL-149) = +0.969 →
without the 27 RISING (n = 122) +1.030, D = +0.061 = a quarter of
one grid step (ĉ₁ 0.485 → 0.515); paired bootstrap (200 reps, 0
skips) D pct 5/50/95 = −0.143/+0.014/+0.252, P(|D| ≥ 0.25) =
0.065. Co-reads: CONV-only (n = 86) +0.703 — dropping AMBIG too
is a DIFFERENT contrast (reported, not lettered); the GD−DD gap
−2.647 → −2.831 without rising galaxies: the GD tension is not
rising-carried (it strengthens slightly; consistent with 9M's
point-dead artifact). Paper-2 consequence: one-sentence
robustness row + this pointer. Wall-clock 3.0 min.

**9O O-FRAGILE by the letter — and the MAP is the deliverable
(all gates first-run: G9O-0 LNPI bit-identity 4/4 × 2 files;
G9O-1 analytic; G9O-2 twenty regressions exact to print
precision).** Three structural facts: (i) OPER ≡ G-LOW to 3
digits in ALL rows — the clean-strata kinematics put ~zero mass
at fcomp ≥ 0.2 on their own, so the shipped envelope was
operatively its low-conversion edge all along (the 7J-z
D2-inversion re-measured from the prior axis). (ii) WITHIN the
measured band every movement is Newton-ward: the α ≥ 0.5
exclusion only strengthens (EX05 −13..−20 at OPER → −25..−61 at
G-MID/G-HIGH) and at G-HIGH the interior lean dissolves (α_marg →
0.000–0.028, dN to −10.5; L2 broken 2/4 — at high conversion
Newton outright wins the drop world). (iii) Only FLAT — no
companion-rate information at all — revives the boost world:
α_marg 0.30–0.66, P(α ≥ 0.5) up to 0.9992 (simple 31), Newton
re-rejected +11..+16, EX05 −6.4..−9.7 (L1 broken 2/4, marginal).
LETTERS: L1 (α ≥ 0.5 excluded) HOLDS at every measured-prior
anchoring and weakens only priorless; L2 (no Newton verdict)
breaks in BOTH directions at the family's ends. Quotable band:
a_marg(DROP) = [0.000, 0.664] across the family. Paper-1
consequence (pre-stated): the α band gains an explicit
prior-conditional column — operative quote: "on clean strata
α ≥ 0.5 is excluded at every anchoring of the measured
companion-rate prior across its conversion band; the mild
0.1–0.3 interior lean is the low-conversion reading and tightens
toward α = 0 at mid/high conversion; only a flat prior (no
companion-rate information) revives α ≈ 0.5–0.7, by letting the
no-companion world back in." The upper-limit FORM survives
everywhere; the interior-lean clause is now explicitly
conversion-conditional.

PLAIN VERDICT: 9N SUCCESS (the dial is rising-flag-robust — the
row the addendum asked for). 9O SUCCESS as an instrument (the
letter says FRAGILE and that is the finding: the fragility is
mapped, the exclusion is not fragile, the lean is).

ELI12: (9N) We asked "does our galaxy dial change if we throw
out the galaxies whose rotation curves are still climbing at the
last measured point?" No — the needle moves a quarter of a tick,
inside the noise. (9O) We asked "does the binary answer change
if we change how strictly we count hidden companion stars?" The
big NO (no strong boost) stays no matter what; the small MAYBE
(a tiny boost) lives only at one gentle end of the assumptions —
and if we pretend to know nothing about companions, the old
big-boost story comes back. So the paper shows a labeled dial,
not one number.

NO credence movement (pre-stated both stages; anomaly-real 53).

## PAPER v4.0 (2026-08-06): THE 8/9-SERIES ABSORB — the author called the thaw

The frozen v3.9 (Stages 1–7J-z8) absorbed the rivals-and-census
arc, the GD dial-tension arc, the quality arc, rounds 9–13, and
9N/9O. Structural changes: header v4.0 arc block; abstract (3)
extended to the upper-limit endpoint + the ≥ 2.1× standing
measurement; §1 trajectory rewrite; **NEW §4.7** (the GD dial
tension, eight controls, stated as an open anomaly against every
function in the ledger ours included); §6.3 + three blocks
(rivals-and-census / quality-stratification / the noise meter,
closing with the operative upper-limit endpoint, the 9O prior
dial, and the D3 promotion gate as standing rule); §7.2
census-defense addendum; §7.4(d) the 8D settling (BOOST-CARRIED)
+ the compatibility close; §8.1 status re-grade (the upper-limit
sector maps to Q₂ ≲ 1–2× the Cassini comparator — the
binary-calibrated tension no longer independently established);
§9.1 three status notes (binary rows co-quoted with the upper
limit; the GD anomaly row; the 8A/8B rivals ladder); §9.2 item 5
rewritten to the stratified prescription + promotion gate; §10
bullet-1 close + census bullet → the (band, cliff) pair form +
credence paragraph → the full pre-signed trail ending at **~53**
with the banked-regardless list; App A → **NINETEEN** (new #19 =
the 9E double-ratio reading, round-13 E2); App B → LOG/manifest
pointers + the full 8/9-series script map. All edits via the Edit
tool (the PS 5.1 corruption hazard avoided); grep audit: no stale
"eighteen"/"v3.9" outside the header's chronological narrative.
Paper remains not-for-circulation before colleague review.

THE SPLIT assembles from this base next (Paper 1 binaries, Paper
2 galaxies+coefficients). **Paper 3 (mechanism) DEFERRED by the
author (2026-08-06): ship 1–2 first, then "continue a bit with
the third" — more mechanism work (the O5 seam) before it is
written.** Author delegated titles/author-line/credence-placement/
spin-offs to Claude (2026-08-06).

## REVIEW ROUND 14 (2026-08-06, journal-referee pass on Paper 1 draft 0.1/0.2; fresh Opus agent, referee brief incl. readability per the author's request): MAJOR-REVISION (LIGHT) — two real catches, both adopted

Verbatim archive REVIEW-ROUND14-OPUS.md (uncommitted). 38
spot-checks: 30 PASS, 4 PARTIAL, 1 FAIL, 3 unverifiable-in-detail;
endpoint discipline certified ("no endpoint overclaim found");
reconciliation fairness to all five cited groups checked against the
record and confirmed. THE TWO CATCHES (both verified on disk before
adoption, both DRAFT errors — the record was right): (1) the draft
attached the eleven-pair leakage probability (3.8e-9) to the
nine-pair corrected census; P(≥9 | 0.9) = 4.8e-7 (Poisson,
recomputed independently) — FIXED, both conventions now quoted with
their own tails; (2) "monotone dose curve in RUWE" — Q1/Q2 are
mutually inverted (8Z's own verdict: NON-MONOTONE/MIXED-CARRIED);
FIXED in the draft (honest two-sentence form + Table 3) and
CLARIFIED in the two ambiguous PAPER.md clauses (dose LADDER,
labels + inversion note; record numbers unchanged). Also adopted:
wobble-law cost re-quoted at defensible multiplicities (~100–500;
the 486–1078 was the retracted-scale conditional); FLAT-prior
revival quoted as 0.30–0.66 (not "near calibration"); PSS
de-lumped from the Cookson <1σ measurement; ṽ-only/joint agreement
0.08 → ≤0.12; arm-D recovery stated precisely (0.48 discriminator,
one arm 0.27); ρ = 0.47 sliced (+0.24 all-pairs); model-Newton
zero-point ~0.98 clause added to the median section; 8D
removed-subsample numbers added (0.62–0.85 at +21–26 — his
underclaim catch); App A → bulleted with correction #15 added;
"not aware of a published quantification"; "validation arms" →
injection–recovery tests; occupation/y glossed; Tables 2/3/4
added. QUEUED for draft 0.4: Table 1 (sample-cut ladder with
per-cut counts — needs a small counting script), the global
sentence-length pass (mean 27.9 → ~20; the 80–108-word semicolon
chains), the phantom-veto figure (§7 headline has no figure), full
script filenames in App B, and the [verify] references at the
pre-circulation pass. Draft 0.2 → 0.3 committed.

PLAIN VERDICT: SUCCESS — the referee did the job the round was
designed for (one arithmetic error, one mis-description, both
mine, both caught before any reader saw them), and the paper's
spine survived 38 adversarial spot-checks.

ELI12: We asked our toughest critic to read the paper like a
strict journal editor. He checked 38 of our numbers against the
raw files: 30 perfect, a few needing footnotes, and two real
mistakes — I had put the wrong "how unlikely is this" number next
to the star census (the number belonged to the OLD count of 11,
not the new count of 9), and I called a staircase "always going
up" when its first two steps are actually swapped. Both fixed.
He also said: fewer giant sentences, more tables. We added three
tables now and queued the rest.

NO credence movement (paper round; anomaly-real 53).

## PAPER 2 FIGURES + REVIEW ROUND 15 (2026-08-06, journal-referee pass on Paper 2 draft 0.2; fresh Opus agent, same referee brief): MAJOR-REVISION (LIGHT) — two real catches, both verified on-disk and adopted; draft → 0.3

The figure stage first: calcs/paper2_figures.py, seven gates ALL PASS
first run (verbatim SPARC loader counts 153/2700/149/38/111 + 12
fiducial lensing rows asserted; the 4S/4Z profile tables re-derive
their own printed Deltas, 56.30/56.3 and 28.91/28.9; every
paper-quoted number regressed against its stage output before
plotting). Five figures: RAR+family+lensing with a residual panel,
c1 profiles, a0 ladder, scatter+inner-disk decomposition, GD dial
split with the eight controls. Provenance dump data/paper2_figs.txt.

The round: ~40 spot-checks, 38 PASS; figures all PASS; register
certified the cleanest yet (2 structural em-dashes, no codenames,
no mid-prose bold). TWO REAL CATCHES, both mine, both
verified against outputs before adopting:
- (A1) §5 said the pre-vertical hierarchical a0 "drifted HIGH" —
  backwards. On disk (5D/5L): 0.79–0.99e-10, LOW, traded against
  high f_ML 1.40–1.65; only the BINARY translation sat high. I had
  conflated the two rows. Fixed to "drifted low, traded against
  high fitted mass-to-light ratios."
- (A2) §7 said the quarter-coefficient function is "excluded once
  the vertical channel is active" — overstated AND self-
  contradictory with §4 (whose hierarchical profile peaks near ¼).
  On disk (5M): dv-ON boot still sits 8.55 AHEAD of BE; what the
  vertical channel does is collapse its lead 76 → 9. The outright
  rejection was the FENCED binary result (O7), which dissolved into
  the ±8 degeneracy under the landed posterior (7J-d). Fixed to
  collapse language + the fenced-history clause + an explicit
  profile-vs-ledger distinction sentence.
Also adopted: a Methods subsection (new §2.2: the nu_lambda family,
likelihood, three M/L treatments, three interval conventions, the
gate discipline — the referee's "not self-contained" item, correct);
"every function in the ledger" → "all four functions tried" (§5);
alpha defined at first use (§9); "discovery/discovery-grade" softened
to "anomaly" everywhere (his argument: the stages themselves pre-
stated NO-credence-movement — right); the Newton ≥ +7.9 row
qualified as full-sample with the stratified upper limit named
operative; five long sentences split; floor 0.034 → 0.035 (0.0354 on
disk; 0.034 is Desmond's number, now phrased as consistency); §6.2
"outer" → "mid-disk" (the gate computes inner/mid = 2.39); abstract
"tenfold" dropped; data-availability line added. His three §9
SUSPECTED items upgraded by direct reads this session (6W verdict
0/6 + EXCLUDED; 7G 451 orders B2 PASS; 7H interior 6/6).
Report archived verbatim in REVIEW-ROUND15-OPUS.md (uncommitted).

REFERENCE PASS (same day, Haiku scout + primary checks): P1's three
flagged refs verified (El-Badry & Rix MNRAS 480, 4884; Lindegren
A&A 649, A2 = the astrometric-solution paper; Tokovinin AJ 147, 87
= paper II, statistics) and flags cleared. ONE REAL CATCH: our
"Desmond 2023, MNRAS 525, 6130" coordinates belong to Stiskalek &
Desmond 2023; the 0.034-dex intrinsic-scatter paper is Desmond
solo, MNRAS 526, 3342 (sigma_int in its abstract verbatim). Fixed
in P2 + PAPER.md's list (bibliographic fix, not a numbered
correction: right author, right year, right claim, wrong volume;
this NOTES line at ~1831 recorded the wrong coordinates and is
superseded by this entry). Chae et al. 2021 ApJ 921, 104 confirmed
against our own archived PDF (data/chae2021.pdf, Table 3 = the
109-row extraction).

Plain verdict: figures SUCCESS; round SUCCESS (the instrument did
its job — two direction/grade slips caught before circulation);
reference pass SUCCESS.

ROUND 15 CLOSE (same day, the verification pass, same agent
resumed): all ten required items RESOLVED on his re-read of the
committed 0.4; §2.2 methods spot-checked accurate against the stage
scripts; the rewritten §7 attribution checked against the ledger +
7J-d; cross-paper coherence pass (P1 0.5 ↔ P2 0.4) found NO
numerical or terminological mismatch (α convention, upper-limit
language, degeneracy figures, a₀ translation row, credence ~53 in
both). One residual, adopted: "verified symbolically" → "verified
numerically" (the 4S G1 gate is a numerical series check).
**FINAL VERDICT: ACCEPT.** Pre-circulation polish landed in the
same arc: P2 abstract 259→249 words (inside the 250 cap), final
sentence splits (mean 21.2), and P1 gained its reciprocal
companion-paper cross-references (it previously never pointed at
P2). Part 2 of the report appended to REVIEW-ROUND15-OPUS.md
(uncommitted). Both papers now referee-ACCEPTED at draft grade:
P1 0.5, P2 0.4. The next step is the author's: circulation
(CERN colleague → Zenodo → arXiv per the standing path).

## Stage 9P — THE MASS-MODEL INSENSITIVITY CONTROL (2026-08-06, pre-reg 9fef99d + amendment e72817c, the RNAAS-note prerequisite): P-GRAY, and the bound is NOT mass-model-carried

The circulation phase opened (playbook papers/CIRCULATION.md; the
author called the go). The reviewer's exposure (c-iv) — "show the
≥2.1× inflation is insensitive to the M_G–mass relation or
selection" — became stage 9P: the 9L narrow-pair fpm meter re-read
under locked deformations of the only mass-model object in the
chain. Variants: global mass ×0.80 / ×1.25 (generous vs the ~5–10%
photometric M-L systematic), tilt ±0.15 at pivot 0.82 M☉ (±10–14%
differential across the MS), and the MG window tightened to
[3.0, 13.0]. One orbit run per seed reused (the population is
mass-table-independent); pick stream at fiducial pool lengths so
the rng is variant-identical; all expressions bit-verbatim 9L.

AMENDMENT 1 (first-firing catch, pre-quote, logged): G9P-0 FAILED
on the first run — the eval was wired to the full-sample noise
pool where 9L indexes each stratum's own pool. The GB0w precedent
repeats: the identity gate caught the wiring slip before anything
was quoted. Fixed to verbatim; the re-run's G9P-0 is BIT-EXACT
(max|dT| = 0.00e+00, 8/8; E[fpm] to 4 decimals).

RESULT: envelope over V1–V5 × Q1–Q4 × both seeds = [1.77, 2.52];
Q1-only min 1.81; fiducial Q1 = 2.21/2.13 (9L reproduced exactly).
VERDICT BY THE LOCKED BAR: P-GRAY — missed P-ROBUST (≥1.8) by 0.03
in exactly one cell (V2-high, seed 101, Q2 = 1.77); never remotely
approaches P-FRAGILE (1.5). The softest direction is masses-high
(+25%), as the physics predicts (higher v_c → lower ṽ → less noise
demanded). Consequence per the pre-registered grammar: the note
ships the envelope minimum as its bound — "at least 2.1× at the
fiducial mass model; at least 1.8× under a generous mass-model
deformation envelope." The reviewer's objection is closed: the
pair-error inflation is real regardless of the mass model.

Plain verdict: SUCCESS (the control did exactly its job, including
catching its own wiring on the first firing).

ELI12: A skeptic could say "maybe your star-speed errors only look
too big because you guessed the star WEIGHTS wrong." So we redid
the measurement pretending the weights were 20% lighter, 25%
heavier, tilted light-to-heavy, and with the edges of the sample
cut off. The answer barely moved: even in the worst case the
catalog's error bars are still 1.8 times too small, and usually
more than 2 times. The skeptic's escape hatch is closed.

NO credence movement (pre-stated; robustness round).

ELI12: We built the five pictures for the second paper, with the
same rule as always — the computer refuses to draw anything it
cannot re-check against the original result files. Then our strict
editor read the whole paper again. He found two sentences where I
said the opposite of what our own files say: I wrote that a number
had drifted UP when it drifted DOWN, and I said a rival formula was
"knocked out" when it actually just lost its big lead. Both fixed.
He also made us add the "how we did it" section every real paper
needs. And a helper double-checked our reading list: one book
number pointed at the wrong paper by the same author — fixed too.

NO credence movement (paper round; anomaly-real 53).

## Stage 9Q + REVIEW ROUND 16 (2026-08-07) — the three-lemma round: the mechanism seam gets its first derivation-grade tests

Pre-reg 471bef4 (bars, gates, convention scans, credence map, and the
disclosure block locked before execution; the author authorized the
Opus red-team in the same message that called the round). No sky
re-fits anywhere; SPARC entered only as a table of measured orbital
frequencies (loader verbatim, census 153/2700/149 gated).

Executed (all six gates PASS first run):
- L1 ONE-MODE: the dS bath at its own temperature T_dS = 2.76e-30 K
  holds 2.1e-3..7.9e-3 thermal quanta in the horizon volume (central
  4.1e-3; a 300 K 1 m^3 lab cavity holds 5.5e14); the thermal
  wavelength exceeds R_H in all three conventions (2*pi exact / 7.95 /
  13.99); the dressing frequency x*H/2pi sits below every non-soft
  scale for x <= 1 (worst ratio 0.159), and where a propagating mode
  first exists (x* >= 2pi) the boost is already < 0.2%. Bars 5/5.
- L2 FROZEN: sympy-exact SdS theorem — kappa^2 = Omega^2 - Lambda*c^2,
  so STABLE circular orbits need Omega/H > sqrt(3*Om_L) = 1.449: no
  bound system can reach the Markovian regime. Measured floors: SPARC
  min Omega/H = 30.2 (UGC09133 outermost point, 108 kpc), p5 = 73,
  median 240; slowest binary bracket 5.0e3. The largest bound
  structures (clusters 5.7, Local Group 2.1) populate the theorem's
  edge, as they must. Retrodiction: the rate-based analog classes
  6K/6M/6N implemented the regime this forbids — their failures were
  forced, not unlucky.
- L3 CEILING: one-period average = g^2/(2(g^2+delta^2)) exact; = 1/2
  at resonance; monotone down in detuning. Capped CONSISTENT-toy by
  pre-registration, upgrade rule assigned to the red-team.

ROUND 16 (Opus red-team; verbatim in REVIEW-ROUND16-OPUS.md,
uncommitted by standing rule): he re-derived every number
independently — "the arithmetic is clean throughout" — and attacked
the inferences. Adopted in full:
- A2/A3 (the sharp ones): counting BULK thermal quanta does not by
  itself deliver "one BOUNDARY degree of freedom" — but his own patch
  closes the availability question: horizon multipoles l >= 1 sit at
  sqrt(l(l+1))*H, gapped far above omega_dress <= H/2pi, so exactly
  ONE soft collective coordinate is available in the anomaly regime;
  the ~1e122 horizon microstates are the RESERVOIR that sets that
  coordinate's occupation, not 1e122 exchange partners (his image: a
  drumhead has 1e23 atoms and one fundamental). With the patch the
  honest label is: the single-quantum-mode REGIME is derivation-grade;
  M = 1 itself is a geometric CONSISTENCY CONDITION — the measured
  number stays anchored by 6Y's negative-binomial rejection. L1 letter
  accordingly DOWNGRADED from the bars-mechanical DERIVED-toy to
  CONSISTENT-toy (the pre-reg's down-only rule executing as designed).
- A4/A5: L2 verified by his own hand (SdS/Kottler exact, not an
  approximation; Einstein-Straus makes pure-Lambda the correct local
  tide; H0-vs-H_Lambda ~15% bookkeeping noted). His correction
  STRENGTHENS the lemma: the KMS correlation time of a thermal bath is
  tau_B ~ hbar/k_B T = 2pi/H, so our 1/H was conservative — under the
  correct value even the Local-Group edge freezes (tau_B*Omega ~ 13).
  Adopted wording: the theorem supplies the universal O(1) floor; the
  STRONG freezing (30-5000x) is data-carried. L2 DERIVED stands.
- A6/A7: mean-vs-other-functionals selection is patchable (secular
  orbit-averaging + the force entering linearly in n); the 6O/6P
  reconciliation is certified logically sound (frozen = a statement
  about DYNAMICS, not a classical random draw; a thermal state's
  P(n >= L) is a deterministic property of the state), with the
  mean-field-occupation caveat now stated. 7E's 0.480 flagged as a
  driven-dissipative anchor: it confirms the VALUE, not the mechanism.
- A9 (the round's finding, UNPATCHED): nothing establishes that the
  operative soft-mode exchange is RESONANT. Off resonance the weight
  is g^2/(2(g^2+delta^2)) < 1/2, with generic detuning
  delta ~ (x_int - x_amb)*H/2pi — field-configuration-dependent. The
  1/2 in the lending law is inherited from the 6X toy's resonant
  construction (validated by fit, not derived). Per the pre-stated
  rule the L3 upgrade does NOT fire.
- NEW NAMED QUESTION (O5-detuning): the sky MEASURES the 1/2 prefactor
  (c1 = 1/2 at ~7 sigma from zero; 7E 0.480). Either the physical
  exchange is resonant (delta << g — needs a mechanism), or 1/2 is
  over-determined some other way. A detuning-resolved lending law
  beta(delta) is the falsifiable successor instrument — the sharpest
  new mechanism question since 6Y.

Scouts (scout-grade, no priority claims printed): the mode-counting ->
one-mode statement and the frozen-bath-for-bound-orbits argument both
NOT FOUND / NEAR-MISS; the SdS stability bound is known territory
(Nandra, Lasenby & Hobson 2012, MNRAS 422, 2931 — on-topic; formula
not verbatim-verified there; ours is sympy-exact); the T_dS <->
horizon-wavelength folklore is standard with no single canonical
source. Parallel scout: GD replication (TODO 28) is EXECUTABLE —
Oh et al. 2015 (VizieR J/AJ/149/180; 26 dwarfs, ~8-10 slow, full
baryonic decomposition) is the starting catalog; Iorio+2017
supplements.

CREDENCE (map executed mechanically): with L1 at CONSISTENT-toy and L2
at DERIVED, the map's exactly-one-derived cell fires ->
bath-mechanism conditional HOLDS at ~8%. No 12, no 15 — the review
traded the credence bump for label honesty, the right trade.
anomaly-real 53 UNTOUCHED (mechanism-side round, pre-stated).

Ledger: +mech-9q-onemode, +mech-9q-frozen, +mech-9q-ceiling (185
rows); worldtable +4 tokens, six gates PASS.

PLAIN VERDICT: L2 SUCCESS (a real theorem, reviewer-verified, now with
the correct bath time constant); L1 SUCCESS-WITH-RELABEL (the regime
numbers stand; the M=1 label honestly demoted to consistency
condition); L3 NEEDS REFINEMENT (the detuning question is the named
target).

ELI12: We tried to turn three of our working guesses into proofs.
(1) The universe's heat-glow is so cold that its typical "particle of
heat" is bigger than the visible universe — so a galaxy talking to it
talks to ONE big slow thing, not a crowd. Proven at toy level — but
our checker made us admit this shows one partner is AVAILABLE, not
that our measured "one partner" number follows from it. (2) Nothing
that orbits can circle slower than the universe's own beat — a real
theorem, checker-approved, and his correction made it stronger. The
bath can never look like fast noise to an orbit; it always looks
frozen. That is why all our earlier noise-style lab models failed —
they had to. (3) The "half" in our sharing rule is only guaranteed
when the two humming modes are perfectly in tune, and nobody has shown
they are. That gap is now our sharpest open question: the sky keeps
measuring exactly one half — WHY?

## Stage 9S (2026-08-07) — the detuning bound: O5-DETUNING gets its first number

Pre-reg 4b46cd5 (bars + co-read structure locked; the in-session sketch
disclosed). The ROUND-16 A9 hole made quantitative: if the exchange
behind the lending law runs off resonance, the operative weight is
r = (1/2)·g_c²/(g_c²+δ²) < 1/2 and the AMB tail exponent drops as
p = 1/2 + r·G/2 (G = s_amb² = the measured ambient gate; exact by the
5P family theorem + q_loc → 1, sympy-gated; at r = 1/2 this is 6E's
postdiction that matched both systems — regression reproduced
0.6884/0.5293 to ≤ 4e-5).

INVERSION (every input parsed from archived outputs — no new fits, 7I
freeze respected): at the hier anchor (5G p̂ = 0.65, the
declared-primary galaxy treatment; point-grade, no σ printed) the
measured tail requires r ≥ 0.315 across ALL three gate treatments
(fiducial 0.398 / maxclust 0.346 / noclust 0.315), hence
**|δ|/g_c ≤ 0.77: the sky forbids far-off-resonance exchange** — the
operative weight sits within about one linewidth of the resonant
ceiling, or the coupling is strong. VERDICT S-BOUNDED by the locked
grammar. Honest conditionality stated: the 4H flat-treatment co-read is
too soft to bound r (0.16–0.21; unbounded at its −1σ edge) — the bound
is hier-anchored; binaries cannot bound r at all (gate 0.117, tail span
0.500–0.529 over the whole r range — consistent with everything).
DIVIDEND: P1's void asymptote becomes a DIRECT gate-free r-meter
(p_void = 1/2 + r/2; r = 1/2 ⇔ 3/4) — annotated on the P1 row.
Successors named: the r-ladder consistency re-fit (freeze-labeled) and
the void r-meter (external/DR4). Gates 4/4 first run. NO credence
movement (pre-stated).

PLAIN VERDICT: SUCCESS — the round-16 objection now carries a measured
bound instead of an open flank.

ELI12: Our checker's sharpest question was "who says the two humming
modes are in tune? If they are not, the shared fraction drops below
half." We answered with data we already had: the measured galaxy tail
is steep enough that the sharing fraction must be at least 0.32 of its
maximum — the mistuning can be at most about one linewidth. The sky
itself says: nearly in tune. WHY it is in tune stays open — but "badly
out of tune" is now excluded by measurement, and future void galaxies
read the sharing fraction directly.

## Stage 9R (2026-08-07) — the GD replication attempt on LITTLE THINGS: POWER-LIMITED, honestly

Pre-reg 7f32d19 (construction, arm ladder, per-arm injection power
gates, letters all locked; the replication outcome genuinely unknown
at commit time). Data: Oh+2015 VizieR tables, SHA256-pinned manifest;
the baryonic curve reconstructed pointwise as V_bar² = V_tot² − V_DM²
by absolute-radius matching inside their own published mass models
(ring match 0.95; negative-baryon drops 5.5%; the five no-Spitzer
galaxies excluded by their table flag). Census: GD 16/21 by mass
ratio; SPARC overlap 7 split out by name+alias; LT-ONLY verdict set =
14 dwarfs. The instrument port is exact: OFF-branch bit identity
d = −1.1e-13 and the SPARC-GD OFF profile reproduces 8S-b's
λ̂ = −1.542 to the printed digit; the one forced wiring change (f on
total baryon — LT publishes no gas/star split) passed its SPARC-side
license (−1.168 vs archived −1.31, inside the 0.3 bar).

VERDICT: **R-POWER-LIMITED.** Neither arm passed its own injection
test, so BY DESIGN the sky was never read: ARM-V (verbatim 0.10 error
cut; 65 points / 7 galaxies) recovers injected worlds attenuated —
the GD world −1.31 comes back at −0.52/−0.46 (the two worlds still
separate by ~1.1: attenuation, not chaos); ARM-R (relaxed 0.20 cut;
172 points / 12 galaxies) is worse — the headline world +0.97 is
recovered at −0.34 and a −2.00 grid edge (the added low-quality
points wash the dial out; the familiar edge-riding pathology). The λ
dial is not measurable on Oh+15-grade dwarf curves at this sample
size and error level. NO credence movement (pre-stated); the P2
falsifier #3 stays OPEN, now with a validated, reusable data
construction and a measured power wall.

Named successors: (a) 9R-b = the two-world CONTEST form (λ = +0.97 vs
−1.31 likelihood ratio with free nuisances and its own injection
calibration — coarser and possibly powered; the ~1.1 separations
holding in ARM-V say it has a chance); (b) Iorio+2017 3D-Barolo
re-derivations as alternative input; (c) anchored/DR4-era dwarf data.

PLAIN VERDICT: NEEDS DIFFERENT DATA (or the contest-form instrument)
— and the discipline held: the machine refused to quote a sky number
its own injections could not certify.

ELI12: We took someone else's 26 small galaxies, rebuilt their
ingredient lists from the published tables (95% of the rings matched
perfectly), and pointed our dial-measuring tool at them. Before
trusting it, we fed it FAKE skies where we knew the answer — and it
kept getting the answer noticeably wrong, because these galaxies are
few and their speeds are hard to measure. So we did not let it read
the real sky at all. The question "do independent dwarfs show the
same weird dial?" stays open — with a sharper plan for next time.

## Stage 9R-b (2026-08-07) — the contest form: near-powered, and the bars held

Pre-reg 025b135. Estimation replaced by ranking (D = m2ll at the
headline world minus m2ll at the GD world, verbatim nuisances),
calibrated on 20 injections per world per arm; powered iff ≥ 16/20
own-side on BOTH worlds. Engine tied to 9R exactly (the power-cell
regression reproduced λ̂ = −0.524 to the digit). RESULT:
**B-POWER-LIMITED — by ONE injection**: ARM-R classified the GD world
20/20 and the headline world 15/20 against the locked 16/20 bar
(5/40 total misclassified; ARM-V 9/40). The bars were locked before
the run and they hold; no post-hoc softening; the sky was never read
(the code path reads the sky only after power passes, same as 9R).
Per the pre-registered clause TODO 28 CLOSES at the public-Oh+15
input grade: three instruments (profile, relaxed profile, contest)
each failed their own injection tests on this input. Live successors:
Iorio+2017 3D-Barolo re-derivations as input (the one-injection miss
says a modestly better error budget likely tips the contest into
powered) and anchored/DR4-era dwarf data. Gates 4/4. NO credence
movement (pre-stated).

PLAIN VERDICT: NEEDS DIFFERENT DATA — with the instrument now
measured to be one hair from sufficient.

ELI12: We built a simpler question for the little galaxies — not
"what is your dial?" but "which of two known worlds do you look more
like?" We rehearsed on fakes again: this time the tool got 35 of 40
rehearsals right — but our rule, written down in advance, demanded
slightly better, and rules don't bend after the fact. So the real sky
stays unread, the question stays open, and we now know exactly how
close the tool is: one rehearsal short. Better data from another team
is the named next try.

## Stage 9T + REVIEW ROUND 17 (2026-08-07) — the resonance theorem: the detuning flank closes; the averaging flank becomes THE question

Pre-reg c6195fe (bars, scan corners, credence map locked; the full
in-session sketch disclosed). The theorem: the anomaly sector spans
H/2π (L1) and no bound system has more than ~1/H of coherent time
(L2), so the accumulated detuning phase obeys φ = δ·t ≤ η·Δx/2π ≈
0.16 rad — detuning inside the soft sector is UNRESOLVABLE. Exact
legs sympy-clean; scan floor r ≥ 0.454 (sharp cell 0.4876); five
gates PASS first run.

ROUND 17 (fresh Opus red-team; verbatim in REVIEW-ROUND17-OPUS.md,
uncommitted): every number reproduced independently; ruling
**UNPATCHED HOLE: NO for the scoped claim** — with three conditions,
ALL ADOPTED:
(1) The credence move credits ONLY the detuning closure — NOT
"resonance established". His separability check is the round's gem:
δ is small relative to EVERY rate in the problem (coherent
suppression ≤ 1/π²; decoherent ≤ 0.025; finite-time ≤ ~0.01), so
detuning perturbs the operative weight by ≲10% (realistically
≲0.3%) no matter how that weight is extracted.
(2) Rewordings adopted: "parameter-free / must / cannot / every
epoch" retired; the budget is PER-HUBBLE-TIME (cumulative phase
grows as N_efolds: floor r = 0.454 at z_form ≤ 6.4, 0.436 at
z_form = 10); the envelope metric defined (cumulative; his
instantaneous peaks 0.0102 — both ~1%).
(3) T3 DOWNGRADED-and-scoped: the boost witnesses γ ~ √(boost) —
deep-MOND galaxies (the systems that actually carry the tail
measurement) witness γ ≈ 0.76–1+ (r ≥ 0.479, in-band); wide
binaries witness only γ ≈ 0.29 and are EXCLUDED from the tail
claim.
Also adopted: A2 (sweeps through resonance are deeply adiabatic,
2πg²/(dδ/dt) ≈ 39.5, and an adiabatic crossing REINFORCES equal
sharing); A3 (the Rabi picture is exactly as strong as 6Y's measured
M=1 — dependency flagged); the honest A6 form: the measured p = 0.65
IMPLIES r = 0.398, BELOW the 0.454 floor — "consistent at current
resolution" stands only on the grid-step/no-σ facts (nearest grid
point 0.70 is in-band); the pre-registered kill stays (σ_p ≤ 0.02
demanding p < 0.67 kills the resonance reading).

THE NEW LOAD-BEARING QUESTION — O5-AVERAGING (replaces
O5-DETUNING): every coupling in the problem is O(H), so γ = g_c·t ~
1 — the universe has completed only about one radian of the Rabi
cycle, and the UN-AVERAGED transfer at γ = 1 is 0.23, which would
read out as p ≈ 0.587 vs the measured 0.65. The ½ the mechanism
needs is a time-average the universe may not have had time to take.
His closing line, adopted into the program verbatim: "book the
detuning win; do not let it be read as pinning the ½." One honest
γ-measurement now separates vindication from a p ≈ 0.59
falsification — the sharpest, most concrete kill-lever the mechanism
has ever carried.

CREDENCE (map executed): T-THEOREM + UNPATCHED-HOLE-NO ⇒
**bath-mechanism conditional 8 → 12** — the first rise since the
6K/6M/6N strike era — crediting only the detuning closure.
anomaly-real 53 UNTOUCHED.

Ledger mech-9t-resonance (189 rows); worldtable 230 tokens six gates
PASS; PREDICTIONS P1 annotation extended (void saturation band,
conditional on the averaging flank).

PLAIN VERDICT: SUCCESS with the scope his review enforced — the
round-16 hole is closed, and the mechanism's fate now hangs on one
nameable, measurable number (γ, the cosmic swing count).

ELI12: Last time our checker asked "who says the two humming modes
are in tune?" This time we proved: the universe is too YOUNG for any
mistuning inside this soft sector to matter — being out of tune
needs time to show, more time than exists. He checked every number
and agreed. Then he found the next real question: the "half" we keep
measuring is an AVERAGE over many swings — and the universe may have
completed only ONE swing. If it really is one swing, the sky should
read 0.59 where it reads 0.65. Counting the swings is now the whole
game — and either answer teaches us something enormous.

## Stage 9U (2026-08-07) — THE GAMMA METER: the swing-count question meets its error bar, and the error bar wins

**Question.** Where is the galaxy tail exponent p CONTINUOUSLY, and with
what σ? The three pre-registered landing bands (fiducial gate g = 0.7536):
one-swing 0.5866, theorem floor 0.6711, full-averaging 0.6884 — 4σ apart
if σ_p ≤ 0.02 (the 9T kill requirement). Pre-reg 21376aa (bands, gates,
verdict letters, credence map all locked before any run).

**Instrument.** ν_p fine profile (23 points, 0.01 steps through the band
region), 5G plain-hier machinery verbatim + 5M/6J vertical-hardened
machinery verbatim (primary per the 7I freeze); injection power gate
(mock skies at p_true = 0.62/0.70 from the fitted model, rng 202);
40-rep paired galaxy bootstrap (rng 53, 6J warm-lite spec), per-rep
parabolic minimum. Gates: G9U-0a 5/5 chain regression on 5G (d ≤ 0.004),
G9U-0b both 5M dv-ON regressions (d ≤ 0.002), G9U-1 interior minima,
G9U-2 injections recovered (+0.007/−0.036, within bars), G9U-4
arithmetic rows regress to 9S/9T exactly. G9U-3 AMENDED per the
pre-authorized A-widen clause: 31/40 reps hit the ±0.06 per-rep grid
edge; widened to ±0.09, the minima STILL pile at the widened edges.

**Result (pre-signed letter U-POWER).**
- PRIMARY: p̂ = 0.6471 ± 0.0746 (bootstrap SD, CLIPPED at the
  instrumented window — honest quote σ_p ≥ 0.075; curvature σ 0.0312;
  percentiles 16/50/84 = 0.56/0.615/0.74 = edge-piled).
- Plain-hier co-read: p̂ = 0.6927 (split +0.046 < 2σ_p — no flag; noted:
  plain sits ON the full-averaging prediction, vertical below the floor;
  at current power this is noise).
- Band distances: one-swing +0.81σ, floor −0.32σ, full-avg −0.55σ.
  ALL three bands inside 1σ. The swing-count question is UNRESOLVABLE
  at SPARC-hier grade: σ_p ≥ 0.075 vs the 0.02 requirement (3.7×).
- Mechanism of the power loss (the finding): the per-rep tail wells are
  ~7–10 lnL deep over ±0.09 — galaxy resampling tilts them by more —
  the tail LOCATION is realization-dominated. The 3A lesson (realization
  scatter > model gaps), recurring on the galaxy side.
- γ-dial central inversions (lean-grade only): fid r̂ = 0.390 ± 0.198 →
  γ_inst = 1.35 rad, γ_window = 2.55 rad. Consistency: the 9T boost-
  witness band (γ ≈ 0.76–1+) overlaps the tail-γ̂ — two channels, one
  swing count, consistent at lean grade.
- Caveat disclosed: warm-lite location noise may contribute to the
  scatter; even a 2× overestimate leaves σ_p ~2× the bar. U-POWER robust.

**Error-calibration dividends (honest updating, both directions):**
1. The 9T "honest tension" (r̂ = 0.398 < floor 0.454) DISSOLVES —
   at measured error the floor is −0.32σ away. It was a point-artifact
   of the 5G grid. The one armed kill against the resonance reading
   evaporates at measured resolution (this is NOT vindication — one-swing
   is equally alive at +0.81σ). mech-9t-resonance row annotated.
2. The 9S bound "r ≥ 0.315 in every gate treatment" demotes to
   point-conditional: error-calibrated r(noclust) = 0.309 ± 0.157.
   r > 0 survives at ~2.0σ (p = ½ is 1.97σ below p̂) — the
   far-off-resonance disfavor softens to lean-grade. mech-9s-detuning
   → CO-QUOTED (pointer gal-9u-ptail).
3. The 5G coarse point (0.65) was fine — now it has its error bar,
   and the error bar is the result.

**Credence (pre-signed cell):** U-POWER → bath-mechanism conditional
HELD 12. anomaly-real 53 UNTOUCHED. No review round (no credence move,
nothing verdict-bearing to red-team; pre-signed map executed
mechanically).

**What the instrument spec buys.** σ_p ≥ 0.075 is now the number every
successor must beat by ~4×: (i) the P1 VOID ASYMPTOTE is promoted to
the LIVE r-meter (p_void = ½ + r/2 is gate-free — no g-division of the
signal — and the 9T kill-band [0.727, 0.750] sits at the ceiling);
(ii) anchored-subsample or DR4-era tails; (iii) the finite-γ p(γ) dial
is built and waiting in the stage (inversions are one line once σ_p
drops). O5-AVERAGING stands OPEN with its instrument built and its
data-grade measured.

Ledger: +gal-9u-ptail (190 rows); worldtable 232 tokens six gates PASS.

PLAIN VERDICT: NEEDS DIFFERENT DATA — the γ-meter is built, gated, and
validated; SPARC-hier cannot power it. The honest products: the first
tail error bar (σ_p ≥ 0.075), one tension dissolved (9T's), one bound
demoted (9S's), and the void channel promoted to the sharpest live
successor.

ELI12: We built the swing-counter and pointed it at the sky. The dial
works — we tested it by hiding fake skies inside it and it found them.
But when we shook the galaxy list (leave some out, count others twice —
the honest way to ask "how sure are you?"), the needle swung across the
WHOLE dial. The sky's answer is 0.65 ± 0.07 — and the three candidate
answers (0.59, 0.67, 0.69) all fit inside that wobble. Nobody gets to
win yet. The good news: the wobble number itself tells us exactly what
data CAN answer the question — emptier patches of sky, where the signal
is biggest and the dial doesn't divide it down. And two of our older
claims that leaned on the needle's exact position got honestly softened
today, because now we know how much the needle wobbles.

## Stage 9V (2026-08-07, second close) — THE R-LADDER + the void-channel feasibility: the in-catalog averaging arc closes at measured power

**Two jobs, one stage (pre-reg 662410c).** (1) The 9S-queued consistency
instrument: fit the exchange weight r DIRECTLY — p_i = ½ + r·g_i/2 with
per-galaxy measured Chae gates inside the likelihood (ν_p family, no new
form). (2) Book the void-tail feasibility verdict from today's primary
read (trap #6 executed on the scout's lead).

**Feasibility first (the primary read decided before any fit):**
Pustilnik+2020 Table 1, read from the PDF: the 8 Lynx-Cancer void
dwarfs have V_max = 31.5–80.3 km/s. SPARC analogs at that mass
(V_max ≤ 85, N = 53) have per-galaxy max y_bar percentiles
0.038/0.069/0.174 — the tail exponent lives at y ≳ 1, which these
objects never reach. **The void-TAIL r-meter is DATA-BLOCKED at public
grade** (VGS is the same mass class). Reopen: a void-crossmatched
sample WITH Newtonian-arm coverage (WALLABY DR2 + environment tags,
DR4-era), or the unified LV database gaining environment columns. The
void channel's PREDICTION (P1 saturation band) is untouched — it waits
for data, not for us.

**The r-ladder (verdict V-POWER, pre-signed):**
- PRIMARY (vertical + per-galaxy maxclust): r̂ = 0.3365 ± 0.1869
  (bootstrap, CLIPPED — 25/40 reps still edge at ±0.225; curvature
  0.0726). Co-reads: fid 0.3895, noclust 0.3085.
- **The 9U cross-tie (pinned in advance by G9V-4): fid direct fit
  0.3895 vs the 9U arithmetic conversion 0.3904 — d = 0.001.** Two
  estimators, one answer: instrument coherence.
- r = 0 disfavored by Δ(−2lnL) ≈ 31 in ALL THREE gate treatments —
  but only +1.80σ in honest bootstrap units (converging with 9U's ~2σ
  demotion of the 9S bound; the nominal-vs-bootstrap gap is the same
  realization-scatter wall wearing likelihood units).
- Band distances: one-swing +0.57σ, floor −0.63σ, full −0.88σ.
  UNRESOLVED, same as 9U.
- G9V-2 FAIL (the honest letter-locker): the r_true = 0.50 injection
  recovered at 0.407 (err −0.093 > bar 0.090) — the instrument is
  biased LOW at high r at single-injection grade. Direction disclosed:
  the full-averaging band is NOT more disfavored than face value.
- Gates: G9V-0 code-path identity 0.00e+00; G9V-1 r=0 ≡ BE d = −0.002;
  G9V-3/4 PASS; G9V-5 AMENDED (pre-authorized widen; clipped).

**Credence (pre-signed cell):** V-POWER → bath-mechanism conditional
HELD 12; anomaly-real 53 untouched. r̂ ± σ_r (lean-labeled) supersedes
the 9S numeric bound as the operative averaging-flank statement.

**Arc close.** 9U (tail exponent) and 9V (direct r) hit the same wall
from opposite parameterizations: SPARC-hier realization scatter. THE
IN-CATALOG O5-AVERAGING PROGRAM IS CLOSED AT MEASURED POWER. What
survives in-catalog: r > 0 at ~1.8–2σ lean, the ~31-unit nominal
preference, and the cross-instrument tie. The question itself — one
swing or many — is intact, parked with two instrument specs and named
external successors.

Ledger: +gal-9v-rladder (191 rows); worldtable 234 tokens six gates
PASS.

PLAIN VERDICT: NEEDS DIFFERENT DATA — same letter as 9U, now proven
from both sides; the honest content is the convergence itself (0.3895
vs 0.3904) and the closed feasibility question.

ELI12: We tried a second way to read the swing-counter — instead of
asking "where is the tail?" we asked the sky directly "how much
exchange is happening?", letting every galaxy speak at its own
loudness. The answer agreed with yesterday's to three decimal places —
our two rulers measure the same thing, which is how you know a ruler
works. But the wobble is just as big this way, and when we hid a fake
sky with maximum exchange inside the machine, it read it ~20% low — so
the machine itself says "don't trust me past this precision." And the
eight void galaxies we found? We read their paper: they're all gentle
slow spinners with no fast inner region — the wrong shape for this
particular question. So the big question stays open, honestly parked,
with a signed note saying exactly what kind of sky data will answer it.

## Stage 9W + REVIEW ROUND 18 (2026-08-08) — THE MULTIMODE REDUCTION: the M=1 seam closes at theorem grade, with the reviewer's scope discipline adopted whole

**The seam** (TODO 30a): 6X derived the lending gate with ONE ambient
mode; 6Y selected M=1 from counting statistics; ROUND 16 downgraded
M=1 to a consistency condition. The real environment has many soft
modes — either the s^L gate breaks at K ≥ 2, or the reduction is
forced. Pre-reg 161de43 (lemmas, gates, bars, credence map with
ROUND-18 conditionality — the 9T pattern).

**Lemma A — the bright-mode reduction (CONFIRMED, every gate):** for
linear exchange coupling, H_int = λ̄(a†A + A†a) EXACTLY with A the
coupling-weighted collective mode — dark modes decouple by
construction (including with the internal Kerr present); a product of
independent thermal modes gives a bright marginal EXACTLY thermal at
n̄_A = Σw_k n_k ⇒ P(n_A ≥ L) = s̄^L. Gates: 6X K=1 port regression
exact to 1e-5 (all six values); sympy set exact; bright statistics to
1.4e-4 (reviewer's independent Fock construction: 4.3e-8); dynamics in
the original basis lands on the weighted-mean form at max 0.88% with
the rival forms at 13.9× and 71.5× the residual; K-invariance 0.001.

**ROUND 18 (Opus, fresh brief; report in REVIEW-ROUND18-OPUS.md,
uncommitted): UNPATCHED HOLE — NO for the scoped theorem**, plus ten
conditions adopted verbatim. The scope discipline he enforced:
- Lemma A's operative role is DEFENSIVE: multimode linear baths do NOT
  force NB tails — the 6Y counting flank closes at theorem grade; the
  weighted-mean content is non-trivial only OFF-degeneracy, exactly
  where the theorem's dynamical assumption fails and (his simulation:
  P2 0.262 → 0.016 as splitting runs 0 → 20λ) lending dies. Kernel
  scans = ROBUSTNESS DIAGNOSTICS; the local 6E anchors stay primary.
- Wording: M=1 is "forced by linear SINGLE-PORT coupling into a single
  final channel" (his counterexample: a linear-per-mode multi-CHANNEL
  absorber with distinguishable final states gives which-path-
  distinguishable = number-additive NB gating — the honest dividing
  line, better than bare "linearity").
- New named assumptions: single dress-ward channel; degeneracy
  tolerance δ ≲ λ (order-of-limits stated); independence of ambient
  modes (correlated/squeezed ambients break the thermal marginal).
- G9W-3 rescoped: since H_int = λ̄(a†A + h.c.) is an identity, the
  dynamics test is a consistency/code check discriminating the correct
  formula from heuristics, not independent mechanism evidence.
- L=1 dynamics only; the L=2 gate is statistics-grade (dynamical
  inheritance argued via single-port, not simulated).

**Lemma B — PARTIAL-SOUND (his grade):** the flat-weight exclusion
rests SOLELY on e_N-blindness — his catch: at cutoff 1e-3 flat gives
p = 0.6874 ≈ the galaxy anchor 0.6884 (magnitude-degenerate,
disclosed; the split — flat misses the binary by 0.16 — is the whole
exclusion). The pre-registered G9W-5 bar FAILED as fired on the binary
leg by 0.0004; diagnosis CONFIRMED as a domain-truncation artifact
(x_bin = 1.095 sat outside my (0,1] integration domain — the kernel
sampled only the soft side of center). Addendum (conditions 5/7/8
executed, [calcs/stage9w_addendum.py](calcs/stage9w_addendum.py),
disclosed post-hoc, letter unchanged): on (0,2] BOTH legs are
convexity-only (binary d = +0.0006 at Γ/x₀ = 0.1; reviewer +0.0007);
integrator-convergence gate PASS. The binary anchor's edge-of-sector
sensitivity (x_amb ≈ 1) is a real named asymmetry the galaxy lacks.
The Lorentzian kernel-tail bound graded "sound but weak."

**The budget — his ruling ILL-POSED, re-scoped:** c₄ = s²/192 − 1/720
is a cancellation amplifying relative shifts ×1.55 (he predicted the
exact factor; addendum confirms 1.55 at every width); the pre-committed
rule fired "multimode-soft" at an unphysical scan edge. OPERATIVE
QUOTE (adopted): c₄(L=2) = 0.002536 single-mode, width-conditional
band ≤ 0.5% at the physical Rabi-grade width; future bars on Δ(s²);
"multimode-soft" RETIRED. **The 30c ordering answer: the c₄ rung is
NOT hindered — its single-mode target stands with a sub-percent band.**

**Credence (the pre-signed cell fires):** W-THEOREM-ONLY +
no-unpatched-hole → **bath-mechanism conditional 12 → 15**, credited
ONLY to the defensive theorem (the 6Y/NB flank closed; M=1 forced
within the named scope) — exactly the firming his summary licenses
("no credence move... beyond firming, at theorem grade, that linear
multimode ambients stay single-mode/geometric"). anomaly-real 53
UNTOUCHED. O5-seam state after 9W: multimode half CLOSED (defensive
theorem); surviving constructive seam = WHY-locality (the 9T resonant
selection, reading-grade) + the single-channel assumption.

Ledger: +mech-9w-multimode (192 rows); worldtable 236 tokens six
gates PASS.

PLAIN VERDICT: SUCCESS at the scope the review enforced — the
multimode seam half is closed as a defensive theorem with every number
reproduced independently; the constructive superstructure was
correctly stripped to robustness-diagnostic grade; the c₄ target is
clean for 30c.

ELI12: We worried that our story only worked if the environment hums
with exactly ONE note, when real environments are full orchestras. The
theorem we proved — and our checker re-derived from scratch — says:
when the coupling is a single doorway, an orchestra IS one note (the
doorway only admits one combined wave). So the "one mode" we kept
assuming was never an assumption — it's what doorways do. Our checker
then made us be honest about two things: the theorem defends the story
(nobody can kill it with "but real baths have many modes") without
proving the fancier bits we'd stacked on top; and one of our test
dials had its ruler end mid-measurement (fixed, disclosed). Net: the
foundation got stronger, the decorations got correctly labeled as
decorations, and the next brick (c₄) has a clean target.

## Stages 9X + 9Y + REVIEW ROUND 19 (2026-08-08) — the kernel and the state meter: both computations verified, both conclusions broken, and the break names the program's next real question

**The hammering directive** ("make sure we don't walk past a diamond")
sent two quantum stages at the remaining seams. Both were
pre-registered (9X: 24245a4; 9Y: f81f212), both ran gates-green as
fired (9X X-CLOSED 5/5; 9Y Y-THEOREM 6/6, one wiring fix logged), and
both went to ROUND 19 as a package. **The reviewer reproduced every
number and ruled UNPATCHED HOLE = YES on both conclusions.** The
pre-signed any-hole cell executed: bath-mechanism conditional HELD 15.
Nine conditions adopted verbatim; conditions 1/2/3/6/7/8/9 executed
same-day ([calcs/stage9xy_addendum.py](calcs/stage9xy_addendum.py)).

**9X (the participation kernel) — what stands and what broke.**
Stands: the kernel rows and 9T-form core fits (his independent
thermal-mixture-of-Rabi-Lorentzians model reproduces our rows to <1%
in the core and *explains the width*: g_c = 2√2·λ√k Rabi factor ×
mixture broadening — the "2–3.3× above λ√(2n)" was never anomalous);
the composition measurements themselves. Broke: **(1) "tail VIRTUAL"
RETRACTED** — my λ-ratio null was wrong: a *pure lending*
saturated-Lorentzian core gives R = 3.81–3.88 at 16λ under the
fixed-absolute-δ test (the width halves with λ), above my bar of 3;
measured 3.37 sits BELOW the pure-lending null. **(2) The closure
claim itself — the sharpest catch of the round: I measured the
occupation-sensitivity at the one detuning where it saturates.** At
δ = 10λ the response saturates by n₂ ≈ 12 (hence S-ratio 0.070); at
the physically relevant soft detuning (δ ≈ the full transition
frequency ≈ 290λ) the response is LINEAR in n₂ — confirmed in the
addendum to his spec: slopes 0.976/0.991 across n₂ = 2–100 at
δ = 125λ/250λ, matching P2 = 4λ²n₂/δ². A single soft mode at
n₂ = 1000 contributes ~22% of the binary gate. **The soft-sector leak
is NOT closed; 9X is relabeled X-OPEN.** Locality reverts to
DATA-SUPPORTED (6Y), mechanism-consistent-in-tested-regime. The l-gap
single-channel corollary is PENDING the same missing scale. The
cumulative dispersive shift (Σρnλ²/δ) needs its own budget — 9T does
not bound it (different object).

**9Y (the squeezing bound → the non-thermality meter).** Stands
unconditionally: the covariance mathematics (Gaussian-classical
correlations keep the bright mode thermal — absorbable; TMS squeezing
makes it anisotropic — not), the exact factorization identity, the
dial arithmetic (his reproduction to 4 digits), and **the P1
hardening line: the void kill-band [0.727, 0.750] is
STATE-INDEPENDENT** (squeezed and thermal converge at the void; he
reproduced 0.9678 vs 0.9675). Broke: the *discrimination* claim — a
classical non-Gaussian mixture of thermals shifts the gate the same
direction (+0.0138 at the binary anchor = the r_sq ≈ 0.2–0.3 signal)
and is degenerate with squeezing on the sky's scalar; parity would
separate them but is lab-only. **Renamed: NON-THERMALITY METER**, with
squeezing as the quantum benchmark. And the off-thermal gate
functional is ambiguous — P(n≥2) and [P(n≥1)]² coincide thermally
(both s²) but split ~50% under squeezing, and the coherent-pair
channel ⟨a²⟩ is unmodeled — so **Y-THEOREM relabels to Y-DIAL
(conditional on the P(n≥2) reading)**; the registered DR4-era binary
lever (Δp_bin = +0.013 at r_sq = 0.3) inherits the conditional; the
frozen-window consistency condition (squeezing persistence vs
averaging convergence vs 6N dephasing) is required and unestablished.

**THE NAMED SUCCESSOR — O5-LEAK** (the round's real product; both
open items collapse onto it): (i) an explicit horizon-side mode
density ρ(ω) for the soft sector; (ii) the single missing scale
g_c^sky/H; (iii) the continuum lending-leak integral
∫ρ(ω)n(ω)(g_c/δ)²dω and the dispersive-shift budget, evaluated
against P2_sys. The reviewer's unification: the soft-leak closure AND
the l-gap single-channel corollary are the same computation. This is
now the sharpest well-posed mechanism question the program has.

Credence: bath-mechanism conditional HELD 15 (pre-signed any-hole
cell); anomaly-real 53 UNTOUCHED. Ledger: +mech-9x-kernel (X-OPEN
operative) +mech-9y-statemeter (Y-DIAL conditional) → 194 rows;
worldtable 240 tokens six gates PASS.

PLAIN VERDICT: NEEDS REFINEMENT, and honestly earned — both
instruments computed correctly and both conclusions overclaimed; the
review converted a false closure into the program's next well-posed
question (O5-LEAK), kept the void-band state-independence hardening,
and left a conditional DR4 dial. The 9W defensive theorem is
untouched by any of this.

ELI12: We built two new instruments and pointed them at the last soft
spots in the story. Both instruments work — the checker rebuilt them
from scratch and got our numbers. But both of OUR conclusions were
too big. The first instrument said "faraway notes can't sneak into
the doorway" — the checker showed we'd tested the sneakiness exactly
where sneaking is hardest, and at the real distance the loudest faraway
notes DO push on the doorway a little each, so whether they add up to
trouble depends on one number we haven't computed yet (how many notes
there are at each pitch — now the named next question). The second
instrument we called a "quantumness meter" — the checker showed a
plain classical crowd of mixed temperatures moves the needle the same
way, so it's really a "not-one-temperature" meter, still useful,
smaller title. Nothing about the sky changed today; the story's
foundation (last week's theorem) stands; and we know exactly which
stone to turn next.

## Stage 9R-c — THE IORIO REOPEN, EXECUTED AND SPENT (2026-08-08; pre-reg 050717d)

The TODO-28 reopen clause ("Iorio+2017 curves as input; the near-miss
says it likely tips") executed as written. Data: NO VizieR catalog
exists for Iorio+17 (TAP zero rows) — the canonical release is
finalrot.zip on the second author's site (fetched + SHA256-pinned,
calcs/fetch_iorio17.py; 17 galaxies + a DDO216 scenario2 variant,
excluded W3). Design: SINGLE-AXIS change — the observed curve swaps to
Iorio's 3D-Barolo drift-corrected Vc; Oh+15 V_bar reconstruction,
vertical channel, contest form, seeds, and bars all verbatim 9R-b.
New wiring gated: V_bar² interpolation (extrap/gap guards), distance
harmonization (DDO87 rescaled ×0.961), mass-model availability
(CVIdwA + DDO47 = no usable Oh model, consistent with 9R).

VERDICT: **C-POWER-LIMITED** (locked grammar; the sky was NOT read —
fourth instrument, same discipline). Analysis-integrity gates 4/4:
the engine regression reproduced the 9R-b archive to the digit
(LT-only 14, t* +0.40, means −1.25/+2.04, counts 17/14, power cell
−0.524); census 15 galaxies, 211/249 rings, interp drops 6.6% (bar
25%). But ARM-V = 26 pts / 3 gal (below the 40/5 floor) and ARM-R =
47 pts / 6 gal calibrated 14/20 + 17/20 (bar 16/20 each; mis 9/40).

THE AUTOPSY (the round's actual product): the 9R-b one-injection
near-miss was read as an ERROR-SIZE problem; it was a SAMPLE-SIZE
problem. Iorio's per-ring errors are genuinely cleaner, but he
released 17 of Oh's 26 dwarfs, his best-measured galaxies (DDO154,
WLM, NGC2366, DDO50...) are SPARC-overlap — excluded from the
LT-only verdict set by construction — and his ring grids are coarser:
the verdict design shrank 4× (47/6 vs 172/12) and the contest got
WORSE (9/40 vs 5/40 misclassified). Standing lesson (DIARY): when a
power gate misses by one, autopsy WHICH axis the miss lives on before
predicting that a data upgrade tips it; "better data" that shrinks
the design is worse data for a population contest.

The Iorio clause of TODO 28 is SPENT; the reopen re-closes at
public-curve grade; the only remaining successor is DR4-era/anchored
dwarf data (a joint Oh+Iorio hybrid = a NEW pre-reg with its own
multiplicity accounting, not this clause). CREDENCE: no movement
(pre-stated measurement round); anomaly-real 53 untouched. Ledger
+gal-9rc-contest (195 rows); worldtable 242 tokens six gates PASS.

PLAIN VERDICT: DIFFERENT PHYSICS is not implicated and nothing broke
— the instrument worked, the data honestly cannot answer the
question, and the last public-data door on TODO 28 is now cleanly
closed. SUCCESS as a closure, NULL as a discovery.

ELI12: We had a tie-breaker question about slow, fuzzy little
galaxies, and last month's referee said "not enough clean
measurements to call it." A better telescope-analysis of the same
galaxies came out, so we re-ran the whole contest with those cleaner
numbers — but the cleaner set covers fewer of the galaxies we
actually needed (the well-measured ones are already in our main
catalog, so they don't count as independent), and a sharper look at
fewer players is still fewer players. The contest still can't be
called. We wrote down the lesson — "sharper" isn't the same as
"more" — locked the question until the next big data release, and
moved on without peeking at the answer sheet.

## Stage 9Z + ROUND 20 — O5-LEAK: THE CONTINUUM LEAK INTEGRAL (2026-08-08; pre-reg 9735d0e)

The round-19 mandate executed as written: horizon-side physics input,
not another toy scan. The de Sitter static patch, system at the
origin; the l=0 radial problem for the minimally-coupled (=
graviton-like, TT-proxy) massless scalar turns out to be the λ=1
Pöschl–Teller well — EXACTLY SOLVABLE. Everything downstream is
elementary and sympy-exact: the Dirichlet solution v = tanh(Hx)cos(ωx)
+ (ω/H)sin(ωx); the near-origin coupling-density ratio D_min/D_conf =
1 + H²/ω² EXACT (the dS infrared enhancement, localized); the ungapped
graviton-like leak integral DIVERGES linearly (integrand →
[1/(8π³Ω²)]/ω², exact coefficient) while the conformal sector
converges. Gates 7/7 first-run. Letter as fired: Z-AMBIG (the
dispersive bar failed; leak legs clean).

ROUND 20 (REVIEW-ROUND20-OPUS.md, uncommitted): the reviewer
reproduced EVERY number from scratch — his own Ricci computation
(R = 12H²), his own quadrature — and his one adversarial swing (an
η² tidal enhancement at l=2) FAILED to land: careful multipole
accounting confirms η⁴, so the physical graviton channel (spin-2
starts at l=2) is SAFER than the l=0 proxy we computed. Ruling:
UNPATCHED HOLE YES (qualified, non-fatal); Z-AMBIG AFFIRMED; the
|Δ|/Ω dispersive bar ruled MIS-POSED (the shift is a ~100% dS-vacuum
common-mode Lamb shift, −0.023/−0.025 H UV-completed, |Δ|/Ω =
0.134/0.712 — excluded-if-real by the deep RAR itself, so it must
renormalize; only the x-DIFFERENTIAL residual distorts observables,
~few-% at the deep end, uncomputed vs the c₁ band). His other
catches: the small-g scan near-miss is PERTURBATIVELY INVALID
(per-mode admixture 10–21) and ambient-fragile (x_amb = 0.10 peaks at
1.94× the smallest gate in the same invalid region — neither bound
nor breach); the fiducial g = H safety is WINDOW-carried, so
gap-SUFFICIENCY is not established.

BANKED (reviewer-reproduced, derivation grade): the PT kernel; the
exact IR enhancement; **gap-NECESSITY — the 6Y ambient gap is
REQUIRED by the graviton-like sector's exact IR structure** (a
reading promoted to a theorem-grade necessity); the fiducial
four-order leak safety at the physical coupling (9T's g ~ H);
l-channel locality at derivation grade, spin-2-safe. Eight conditions
adopted; 1/4/5/6/7 executed same-day (calcs/stage9z_addendum.py);
2/3/8 = the NAMED SUCCESSORS, led by **O5-DIFF: the x-differential
dispersive distortion budget, UV-completed, propagated to the
measured c₁ band — THE actual closure test the mechanism now waits
on.** CREDENCE: pre-signed map → bath-mechanism conditional HOLDS 15
(reviewer-affirmed); anomaly-real 53 untouched. Ledger
+mech-9z-leak (196 rows); worldtable 244 tokens six gates PASS.

PLAIN VERDICT: NEEDS REFINEMENT — with the best half banked. The
question "does the horizon actually have the mode budget the grammar
borrows from" now has a real, exactly-solved spectral answer: the
soft sector is there, its infrared is dangerous exactly as feared,
the ambient gap the data had already voted for (6Y/6G) is now
REQUIRED by the geometry, and at the physical coupling the leak is
negligible. What remains is one honest unfinished ledger line (the
differential Lamb-shift budget), not a hole in the sky.

ELI12: We finally checked the concert hall itself. The universe's
horizon really does have a choir of very deep notes, and our math
says an un-protected instrument would be drowned by them — UNLESS
each system's own surroundings hum a bass note that pushes the
deepest choir voices out of reach. That protective hum is exactly
the thing our earlier sky-data rounds said must be there; now the
hall's blueprints say it MUST be there. With protection on, the
choir barely touches the melody. One thing left: the choir also
slightly re-tunes every instrument by the same tiny amount — we
showed the re-tuning is the same for everyone (so it mostly cancels
out), but we still owe the exact arithmetic for the small
instrument-to-instrument differences. The referee re-derived every
formula himself, tried hard to break the safest part, and instead
proved it was even safer than we claimed.

## Stage 9Z-b + ROUND 21 — O5-DIFF: THE DIFFERENTIAL BUDGET (2026-08-08; pre-reg af00120)

Round-20 conditions 2+3 executed as one stage. The renormalization
scheme DECLARED (condition 3 discharged): exactly two universal
absorbers — one bare-frequency constant (the common-mode Lamb shift
renormalizes into the unobservable bare scale) and one a₀ rescale
(a linear-in-x pull is exactly degenerate with a₀) — affine in x,
nothing more. The fiducial pull has an EXACT closed form, Δ_vac =
−(g²/4π²)·ln(1+Ω/g)/Ω (UV-convergent — trap #10 dies analytically;
five-digit anchor regression to the 9Z addendum). Residual after
affine absorption propagated exactly through ν = 1+n_BE: fiducial
dc₁_eff = +0.0049 (bar 0.05), binary split distortion 0.0004 (bar
0.25), 0.0013 dex (vs 0.13 scatter). Fired DIFF-CLOSED — with my own
stability gate GDB-3 RED (fit-range sensitivity).

ROUND 21 (REVIEW-ROUND21-OPUS.md, uncommitted): every digit
reproduced ("arithmetic integrity is total"), then the letter
DOWNGRADED: **DIFF-GRAY adopted.** His grounds, all now booked: (i)
a red validation gate cannot coexist with a clean letter (the gate
was ALSO mis-posed — tail-only extrapolation the scheme never
performs; redesigned deep-window-containing version PASSES, spread
0.0049 < 0.010 — but the fix exposes the statistic as a RANGE:
[0.0001, 0.0049] in-scheme, 0.023 at the 1/x²-amplified deep END =
2.2× inside the bar, not 10×); (ii) the closure is LOAD-BEARING on
the a₀-rescale — S1-only gives 0.131, near-FATAL; the absorber is
legitimate but SPENDS temperature-lock precision (implied a₀ shift
−2.09% = 0.22σ of the measured band — passes, now disclosed); (iii)
the margin is order-unity on two axes the program has not pinned:
coupling (comfortable g ≳ 0.3H, PT-marginal below) and the
D-NORMALIZATION (dc₁_eff ∝ D; the joint pessimal corner g = 0.3H ×
D×3 = 0.053 mean / 0.267 deep-end = past FATAL). His sharpest
structural point: the η⁴ "real-graviton-safer" buy-back CANNOT be
invoked for the leak while the grammar stays O(1) on the same
coupling — the D-provenance and the l=2 story are one seam. Genuine
win he verified: ambient-independence at fiducial is PROVEN
(identical fit for every x_amb — round-20 condition 5 moot at
fiducial). Conditions 1/2/3/4 executed same-day
(calcs/stage9zb_addendum.py); **THE LOAD-BEARING SUCCESSOR:
discharge the D-normalization provenance (round-20 condition 8) —
"pin down D and the closure becomes real."** CREDENCE: HOLD 15
(rise cell needed CLOSED + no-hole; neither; no strike — nothing
fatal at the physical fiducial). anomaly-real 53 untouched. Ledger
+mech-9zb-diff (197 rows); worldtable 246 tokens six gates PASS.

PLAIN VERDICT: NEEDS REFINEMENT — and the refinement is now a
single named object. The dispersive flank did not close, but it
shrank from "an uncomputed budget" to "one normalization constant
away from closed": at the physical point every number is inside
every band, and the whole remaining question is the absolute
strength of the horizon-mode coupling (the D-provenance). Two
rounds, two downgrades, zero strikes — the mechanism's soft-sector
story is now honest to the digit.

ELI12: The choir's common re-tuning turned out to be two effects we
were always allowed to absorb — a shared reference pitch and a tiny
overall tempo change. After absorbing exactly those two and nothing
else, the leftover wobble is small everywhere we can measure. But
the referee showed our "small" was a bit too proudly stated: it's
small-with-an-asterisk (bigger at the very deepest notes, and it
grows if the choir sings louder than our standard assumption — and
we haven't yet measured the choir's true loudness). So: not closed,
but reduced to one question — how loud is the choir, exactly. That's
the next stone.

## Stage 10A + ROUND 22 — O5-NORM: THE D-PROVENANCE HUNT (2026-08-08; pre-reg fc33426, amendments 1db67d5)

The commissioned derivation ("hunt down the derivation, give it all
you got") executed as one stage, four instruments, every bar and
letter wired pre-run; amendments A1-A3 disclosed pre-quote after a
first firing with two instrument-grade red gates (as-fired output
archived: data/stage10a_dprov_run1.txt).

- T1 (exact, sympy): the SdS CONSTRAINT CREDENTIAL — Birkhoff steps
  (R_tr/(dB/dt) = 1/(rB), derivative-free ⇒ vacuum+Λ forces static;
  the SdS family solves Einstein+Λ exactly) and the first law
  dM = −T_c dS_c EXACT all orders along the family; conjugacy
  |dS_c/dM| = 2πc²/(ħH) = 4.7e99 nats/M_sun; horizon multipole tower
  admixture 1.081e-3 × g_c².
- T2: the l-general radial problem V_l = H²[l(l+1)/sinh² − 2sech²]
  (metric-derived, sympy) + amplitude machinery regressed EXACTLY on
  the l=0 banked ratio 1+1/ω² (1e-13). THE E-LADDER: E_l1 = 1+4/ω²,
  E_l2 = (1+1/ω²)(1+9/ω²) — locked at 1e-12 across four frequencies
  each.
- T3: the per-quantum radiative density (cascade-correspondence ħ
  audit; |dQ̄/ds|² = (2/3)μ²s² exact; 32/5 tensor regression) with
  the T2 enhancement included (conservative): J_rad = 1.4e-46 (bin)
  / 1.2e-35 (gal) — UNDER-SUPPLY U = 7.2e-45 / 1.23e-34 vs the
  weakest bridged normalization. Crossover masses (where a radiative
  channel WOULD suffice): 4.8e-15 / 1.6e-4 kg — 44/34 orders below
  any coordinate that moves an orbit's mass.
- T4: THE NORMALIZATION BRIDGE — under 9W's collectivization the
  resonant coupling is the window integral of the same density,
  g_c² = ∫_window J, so the D-normalization CANCELS out of every
  budget observable; the all-radiative budget then rides the
  measured γ alone (9U dial verbatim; t̂ ∈ {1/H, 2π/H} = an
  irreducible convention envelope). Kill-window: mean-FATAL
  {0.050..0.435} H, deep-end {0.050..1.274} H, SAT rows
  past-FATAL-by-construction; 2/4 convention centrals FATAL; the
  leak co-read breaches the smallest gate up to 8.7× at the small-g_c
  end. B-PHYS at the physical normalization: |dc1| ≤ 3.7e-37,
  eps ≤ 6.2e-40.

Fired N-SPLIT-CLOSED by the locked grammar (13/13 gates green).
**ROUND 22 (REVIEW-ROUND22-OPUS.md, uncommitted): every load-bearing
digit reproduced from scratch (his own Christoffels, his own
integrators — the E-ladder to 1e-13, exceeding the stage's claim),
then DOWNGRADE N-SPLIT-CLOSED → N-GRAY, ADOPTED**, on the program's
own rules read back: the letter grammar carried no gate on its
POSITIVE clause (trap #11's cousin, the R21 defect again); a static
first-law conjugacy is not a derived Rabi coupling; and the named
constraint carrier is in tension with the mechanism's own measured
real-exchange character (6X: virtual = the null; 9U: a Rabi phase).
His verdict on the halves: clause (1) SOLID (the radiative
free-graviton continuum cannot supply the measured exchange, M_eff ≥
μ justified, no lever-arm escape, E_l2-conservative), clause (3)'s
dividend FOLLOWS (**the 9Z-b joint pessimal FATAL corner DISSOLVES —
it priced a radiative normalization physics does not supply** — "a
real advance... removes the sharpest thing round 21 left on the
table"), clauses (2)/(4) NOT EARNED (the D-provenance is DISPLACED,
not discharged; the cavity-QED analogy decorative as used). His gem:
he PROVED E_l1 = 1+4/ω² analytically (the doubled-argument collapse
2/sinh²x − 2/cosh²x = 8/sinh²(2x), executing the stage's own hint)
and REFUTED the stage's parity-ladder gloss at l=3 (set {1,2,4}).

ADDENDUM (calcs/stage10a_addendum.py; conditions 1/2/3 + the
computable half of 5; every reviewer number independently re-verified
per the standing rule): his corrected U CONFIRMED (7.2e-45/1.23e-34;
his A=0.1966 catch CONFIRMED — root cause round(gc,6) beating a 1e-9
guard); his l=3 set CONFIRMED at 1e-12 with this program's own
integrator; **the CORRECTED conjecture — E_l = Π_k(1+k²/ω²), k ∈
{1..l+1} \ {l} — retrodicts all four known members and its fresh l=4
prediction {1,2,3,5} LOCKED at 7e-12 (predicted before computation)**;
the confinement bound closes the collective-radiative rescue
numerically (maximal cavity enhancement (c/(HL))³ still leaves
5.3/11.0 orders short). Also booked per his Q5a, correcting my own
worry: **gap-necessity (9Z) SURVIVES the split** — a divergent
integral times any nonzero prefactor diverges; only the gapped leak's
magnitude demotes (η⁴-trivial rather than marginally safe).

END STATE: the D-normalization question is TRANSFORMED — no longer "a
free constant that could be 3× and kill us" (that corner is dead) but
"which physical channel carries the measured O(H) exchange." The
radiative sector is excluded at every defensible normalization
including maximal confinement; the constraint sector holds the only
O(1) credential (exact, static); the NAMED SUCCESSOR (R22 cond 4 +
5-residual) = derive the near-field/longitudinal exchange coupling as
a Rabi-capable quantized channel compatible with 6X/9U, and bound the
6Y cloud's own internal continuum. CREDENCE (pre-signed map,
mechanical): bath-mechanism conditional HOLDS 15 (rise cell needed
no-hole AND SPLIT-CLOSED — neither; NO strike — nothing fatal fired);
anomaly-real 53 untouched (no sky number moved). Ledger
+mech-10a-dprov +mech-10a-elladder (199 rows); worldtable 250 tokens
six gates PASS.

PLAIN VERDICT: NEEDS REFINEMENT — with the biggest half banked. The
hunt did not capture the derivation, but it killed the monster that
was guarding it: the FATAL corner hanging over the mechanism since
round 21 is gone (the reviewer affirmed it), the coupling's home is
narrowed to one sector, a new exact five-member spectral ladder came
out of the chase with one member proven, and the remaining question
is sharper than it has ever been.

ELI12: We tried to find out how loud the universe's choir really is.
We proved it is NOT the ordinary kind of loudness — gravity waves
from the choir would be weaker than our measured hum by a number
with 40 zeros, and no clever echo chamber fixes that (we checked the
best possible echo chamber: still short by a million-fold or more).
That kills the scary scenario where the choir was loud enough to
wreck our melody — that worry is over, permanently. What we have NOT
yet done is show how the OTHER kind of connection (the quiet,
always-on pull that binds things — the same kind that holds the moon)
can carry the beat our measurements hear. Our referee checked every
single number, proved one of our new formulas himself, corrected our
guess about the formula family (and our corrected guess then
predicted the next family member perfectly), and told us honestly:
you found where the answer is NOT, which is real progress — now go
find where it IS.

## Stage 10B + ROUND 23 — O5-CARRIER: THE CARRIER LEDGER (2026-08-08; pre-reg 8b46d1c)

The chase's second leg ("Let's chase it. All you got"), and the hunt
changed shape mid-stalk: constructing the near-field vertex R22
demanded revealed it is DISPERSIVE-class, not exchange-class — so the
stage built the COMPLETE carrier ledger instead and let the letters
fork. Instruments: T1 the epicyclic dichotomy (κ² = Ω_K² + 4T exact;
secular splitting −(3/2)T/Ω_K exact ⇒ matter modes are either elastic
≥ the measured 30.2H/5e3H floors or SECULAR; the 6Y cloud's internal
continuum is doubly dead — Wien e^−190 AND ≥29H detuned: admixture ≤
10^−84.8·g², R22 condition 5-residual CLOSED). T2 the secular carrier
(the ONE matter family in the soft band — the 10-kAU Galactic-tide
apsidal rate is 2.49H — excluded AS carrier by the dictionary:
ω_sec/ω_dict runs s^{5/2}, a 6-order sweep ⇒ 12 orders of implied a₀
drift vs the measured ~10% constancy). T3 the constraint vertex: an
exterior-sourced l=0 potential is interior-CONSTANT (Laplace) ⇒ it
couples only to total energy, H_int = (Ĥ_sys/c²)δΦ, charge
equivalence-principle-forced — and [H_int, N_s] = 0 EXACTLY: the
vertex is dispersive, gravity's constraint sector has no l=0
off-diagonal soft vertex; single-port ⇒ 9W applies exactly (the S
microstates present ONE collective mode — the R16 patch mechanized);
under hypothesis H-P (per-microstate Planck amplitude) the bright
amplitude is √S·cHℓ_P = √π c².

Fired V-DISPERSIVE-FORCED 8/8 gates. **ROUND 22... ROUND 23
(REVIEW-ROUND23-OPUS.md, uncommitted): every digit reproduced (his
own Fock matrices for the selection rule), then DOWNGRADE →
V-DISPERSIVE-PERMITTED, ADOPTED**, on three catches that all stood
re-verification (calcs/stage10b_addendum.py): (1) the
MOMENTUM-CONSTRAINT vertex p̂_s(b+b†) IS off-diagonal — clause (2)'s
"no off-diagonal vertex, ever" refuted in principle (his ceiling: the
l=1 channel is still ≥15 orders dead as a carrier — wrong DOF
[free-fall-removable dipole] + gapped/empty/geometry-suppressed
partner — CONFIRMED: n(2π√2) = 1.38e-4, admixture 3.5e-4·g², geometric
1.3e-22/5.5e-12); (2) THE HOLE: my amplitude rows √(π(2n+1)) =
2.5–8.1 were the RMS JITTER (mean-zero), mislabeled as the JC
(2n+1)-pull; his zero-detuning mean pull π(2n+1) = 6.3/44.8/66
CONFIRMED; the √S bridge is H-P (an assumption), not 9W (which
collectivizes the COUPLING, not the amplitude) — the entire 61-order
bridge is the one hypothesis; (3) the 6X demotion was asserted — the
program's one real-vs-virtual discriminator (which favored REAL) was
never re-run under the dispersive hypothesis. Also adopted: the
corotation correction (a real slow disk family; excluded by
measure-zero + dictionary), the √π-identity demotion (a tautology —
S's definition rearranged), and the disclosure that V-EMPTY was
unreachable at physical occupations (the fork was cosmetic; trap-#12
failed substantively a THIRD round running — the letter was wired to
the half that could not miss).

THE EARNED STATE: real-exchange carriers are EXCLUDED everywhere in
the enumeration (C1–C7 + the momentum loophole) — "no real-exchange
channel carries the grammar, anywhere" (his words) — and any
microphysical carrier must be DISPERSIVE; whether the dispersive
channel SUPPLIES the measured grammar is NAMED-NOT-DERIVED (H-P
unproven; jitter-not-pull; 6X un-rerun). THE PREVIEW banked in the
addendum (labeled, not a verdict): restoring the frequency-ratio
structure his zero-detuning bound omits gives pull/ω_s ~
x(2n_BE(x)+1)/2 = (x/2)coth(x/2) = 1.098/1.002 at the anchors, exact
series **1 + x²/12 − x⁴/720 — the BE ladder's own Bernoulli rungs
appearing in the mean dressing**. His condition 8 (the real amplitude
bar: derive the honest mean pull to the measured c₁ = ½) is exactly
where that points — THE 10C OPENING MOVE, with his conditions 6 (the
6X dispersive control re-run) and 7 (justify or drop H-P) as the
package. Condition 9 executed: the r = ½ kill-test REGISTERED in
PREDICTIONS.md (dispersive reading ⇒ no γ-running, p_void = ¾
EXACTLY; r < ½ at >2σ kills it while off-resonance exchange
survives). CREDENCE (map, mechanical): bath-mechanism conditional
HOLDS 15 (hole blocks the rise; nothing fatal — the pre-signed
V-EMPTY strike cell did NOT fire and he declined to trigger it);
anomaly-real 53 untouched. Ledger +mech-10b-carrier (200 rows);
worldtable 252 tokens six gates PASS. NO row annotations executed
(the interpretive-layer demotion list is DEFERRED pending the 6X
re-run — his ruling that clause 5 was under-demonstrated).

PLAIN VERDICT: NEEDS REFINEMENT — with the map redrawn. The chase did
not catch the derivation, but it cornered it: every way gravity could
carry the measured coupling by real exchange is now excluded at
derivation grade, the one vertex that remains is dispersive by an
exact selection rule, and the first honest look at its mean pull
lands order-unity with the ladder's own Bernoulli coefficients in it.
The question is no longer "which channel" — it is "does the
dispersive channel's exact second-order structure reproduce the
measured law," and that is a single computable object.

ELI12: We finished searching every room in the house for the thing
that could be trading energy with our melody — and every room is
empty. Real trading is impossible: the notes are either too high, too
out-of-tune, or the doorway itself (we proved this exactly) only lets
the melody's LOUDNESS be felt, never a note passed through. So if the
universe's hum shapes our melody, it must do it the quiet way — by
leaning on it, not by trading with it. Our referee agreed the rooms
are empty, then caught us calling the leaning "proven" when we had
only shown it's allowed — and he found we had measured the WOBBLE of
the leaning where we needed its steady PUSH. When we redid the push
arithmetic honestly, the first numbers came out almost exactly 1,
with the very same magic fractions (1/12, 1/720) our sky formula
carries. That is either a beautiful coincidence or the beginning of
the real derivation — next hunt decides.

## Stage 10C (2026-08-09): O5-PULL — the dispersive closure (R23 conditions 8+6+7). P-TRANSMITTED, GATE-FLAT+ANTI

Pre-reg 13c3662; amendment A1 (bug-class, pre-quote) 3d81fdd; script
[calcs/stage10c_pull.py](calcs/stage10c_pull.py) →
data/stage10c_pull.txt (run-1 as-fired archived _run1.txt).

**T1 THE SEPARATION THEOREM (exact; the honest Delta-structure).** The
10B vertex H = ω(N+σ) + Ω b†b + g(N+σ)(b+b†) is polaron-solvable:
E(n,m) = ω(n+σ) + Ωm − (g²/Ω)(n+σ)² (G1 numeric 9e-15). The two
time-orderings cancel the bath occupation EXACTLY ((m+1)/(−Ω) + m/(+Ω);
G2 symbolic 0) ⇒ **the honest dispersive mean pull carries NO ambient
occupation**; R23's over-production π(2n+1) = 6.3/44.8/66 was the
JC/EXCHANGE template at the dispersive vertex — the wrong template (his
rows reproduced to his digits as the labeled objects: jitter
√(π(2n+1)), zero-detuning exchange-template; G3). The pull reads the
SYSTEM's occupation, linearly, all orders in g: dω(n) =
−(g²/Ω)(2n+1+2σ). T1d (own overclaim struck): σ (EP vs number charge)
enters only as an n-independent constant → absorbed by the 9Z-b
renormalization → **c₁ = ½ is NOT a vertex-selection result**; the 10B
preview's coth-vs-ν gap = x/2 exact = the same absorbed constant.

**T2 THE TRANSMISSION MAP + κ CLOSURE (condition 8 answered).** Thermal
mean → anomalous softening ∝ n_BE(x) → under the standing response
dictionary (compliance-additive, used-not-derived, disclosed):
**ν = 1 + κ·n_BE(x), κ = 4g²/(Ωω)**. At κ = 1 the ENTIRE measured
ladder is automatic (G4 exact: ½, 1/12, 0, −1/720 — all Bernoulli
rungs at once); the polaron pull's strict linearity matches the
measured law's exact linearity (the multiplicative alternative
diverges at x = ln 2 mid-transition — grossly non-RAR). κ METERS
(archived, no new sky fits; G5): a₀-lock κ = 1.00 ± 0.05 (4H/4V/5M:
fitted a₀/horizon = κ²); flat-M/L c₁-meter κ = 1.100 (0.962–1.230);
TENSION disclosed: 4Z hier profile ⇒ κ = 1.48 (1.38–1.58) excludes 1
(the ¼-vs-½ OPEN contest; 5T deep arm + 5K binaries vote ½). Closure
constant g = √(Ωω)/2: binary anchor lands INSIDE the 10A γ-band
(0.083–0.118 H), galaxy legs sit below (0.009–0.030 H) — order-grade
observation, no bar.

**T3 THE REQUIREMENT CURVE; H-P DROPPED (condition 7).** κ = 1 forces
EXACTLY φ_req = ½√(Ω/ω) = ½√(x_amb/x_loc) ∈ 0.19–0.74 at the anchors,
RUNNING ×2–3.2 across the measured windows (zero-point-like 1/√ω at
fixed environment; √Ω across environments). Exact restatement
(identity): φ_req = √(ħΩ/2E_c) with **E_c = 2ħω — the bright mode must
couple with the zero-point of an effective inertia equal to TWO QUANTA
of the dictionary mode it dresses** (named requirement, not a
derivation). H-P supplies the CONSTANT √π = 1.772: level ×2.4–9.4 high
AND structurally wrong running ⇒ **DROP-AND-REPLACE** (G6b criterion;
branch algebra-determined, disclosed). Single-horizon semiclassical
φ₁ = Hℓ_P/c = 1.22e-61 → 60.2 orders under the smallest requirement.
The 61-order H-P bridge is replaced by one exact curve.

**T4 THE 6X DISPERSIVE CONTROL (condition 6): GATE-FLAT+ANTI.** G7
regression: the 6X resonant lending channel reproduces verbatim (ratio
slope +0.989, rms 0.0046, 6×/29× margins). Amendment A1: run-1's
Ω_b = 0.8 = CHI was a parameter collision (resonant sidebands
|1,m⟩~|0,m−1⟩, |2,m⟩~|3,m−1⟩ drained the pair manifold; fingerprint
P̄₂ → ¼ = four-level equipartition; fixed Ω_b → 0.53, gd → 0.15, bars
unchanged). The two dispersive arms then land the signed predictions
exactly: **ARM-A (linear vertex, polaron-compensated = the 9Z-b
renormalization mechanized): saturated weight ½ FLAT (spread 0.0033;
ε-independent) — the dispersive channel is occupation-BLIND at
saturation (existence-counting is Franck–Condon-blind). ARM-B
(cross-Kerr): the exact per-sector law Σ pₘ·4ε²/(8ε²+g²m²) ≈
½/(1+n̄) — THE ANTI-LENDING LAW (gates on the ambient being EMPTY),
opposite monotonicity to the measured gate; GB-exact ≤ 3e-4; anti
slope +0.961 rms 0.0042 vs ratio rms 0.287.** ⇒ a purely dispersive
ambient coupling does NOT reproduce e^{−x} in either tier; R23's
proposed completion of the 6X demotion FAILS; **6X's real-exchange
reading of the GATE channel STANDS**. Pre-signed two-channel
synthesis: dispersive amplitude (l=0 constraint) + real-exchange gate
whose only ledger-surviving partner is the LOCAL 6Y cloud via the
near-field l=2 vertex (the 10B horizon-partner exclusions do not apply
to the local cloud) — the R22-condition-4 derivation returns as the
required next leg. Annotations per R23 cond 5: 6X KEEPS real-exchange
(partner re-identified); 9T NOT mooted.

LETTER: **P-TRANSMITTED** (tokens THM-SEP:PASS, LADDER:PASS,
KAPPA:CONTAINED-WITH-TENSION, REQ:CURVE-EXACT+HP-DROPPED,
GATE-FLAT+ANTI). NOT claimed in any cell: mechanism closed; κ derived;
the response dictionary derived; priority on ν = 1+n (C&T 2019).
CREDENCE: rise cell 15→18 CONDITIONAL on ROUND 24 ruling no unpatched
hole (pre-signed); anomaly-real 53 untouched (no sky fits).

PLAIN VERDICT: SUCCESS at transmission grade — the honest pull is the
right object (bath-free, system-linear, whole ladder at κ=1, κ
measured ≈1), the H-P assumption is replaced by one exact requirement
curve, and the gate channel is proven non-dispersive in the toy class;
closure still awaits κ = 1's derivation (the requirement curve) and
the l=2 cloud vertex.

### ROUND 24 (same day): P-TRANSMITTED AT CONSISTENCY GRADE — hole YES (mild), HOLD 15, two clauses retracted

The round-24 referee reproduced EVERY load-bearing number (his own Fock
matrices, own sympy, own Rabi law for ARM-B) and **independently
validated amendment A1 with the single-fix combinations the stage never
printed** (Ω_b-fix alone spread 0.0269 = PASSES the locked bars;
gd-fix alone 0.1400 = FAILS ⇒ the collision was the culprit, the
gd-halving is FC insurance — amendment certified bug-class). His three
affirmations: T1's separation theorem is "solid and correct — R23's
objection was genuinely mis-templated" (the π(2n+1) row was the
cross-Kerr/JC shift reading the BATH occupation; the longitudinal
polaron pull reads the SYSTEM occupation and is bath-free); T3's H-P
drop = "honestly discharged, the round's best work"; T4's
GATE-FLAT+ANTI negative = "robust" (holds via the anti-law, the flat
arm, AND run-1's decline — none is the growing lending gate). His
circularity probe ANSWERED-and-adopted: the renormalization split is
affine-in-x ONLY under the closure φ ∝ 1/√ω (under constant φ the
"constant" runs x² and is NOT absorbable) — licensed-not-circular but
closure-conditional (condition 3, booked). RULING: unpatched hole YES,
mild ("true-but-inflated," materially narrower than rounds 19–23):
**(1) "κ MEASURED ~1" retracted → the two-pole form (deep-
normalization-consistent 1.00±0.05; transition-shape-contested — 4Z
hier 1.48 excludes 1 at point grade, bootstrap band 0.6–1.8 contains
it; ¼-vs-½ OPEN); (2) "transmits the ENTIRE measured ladder" retracted
→ the polaron transmits the LINEARITY in the local occupation (which
uniquely selects the additive C&T law over the divergent multiplicative
form); the coefficients are C&T's own Taylor series, and only c₁
(weakly c₂) is measured.** Conditions 1–5 executed in
[calcs/stage10c_addendum.py](calcs/stage10c_addendum.py) (every
reviewer number re-verified per the standing rule: population table
exact, both scans exact, κ=1.5 rungs 1/8 and −1/480 exact, bootstrap
band exact); g_close row down-labeled (4/6 anchors below-band);
T4 annotations labeled CONDITIONAL on the l=2 successor. SUCCESSORS
(R24 conds 6/7/8/9): derive κ = 1 (target: is the bright mode's
effective inertia FORCED to E_c = 2ħω?); the l=2 near-field
system↔cloud exchange vertex (R22-cond-4, STRIKE-BEARING: if it cannot
be Rabi-capable the pre-registered 6X/9U strike fires); **the joint
single-κ deep+transition fit (in-catalog, no new data — resolves
κ~1 vs κ~1.5 vs one-κ-form-rejected, the sharpest meter upgrade)**;
lead with linearity. CREDENCE (mechanical): hole blocks the rise →
mech conditional HOLDS 15; anomaly-real 53 untouched. Five O5 rounds,
zero strikes, zero rises — his line: "the derivation is being cornered
honestly, not carried."

PLAIN VERDICT (round-adjusted): the transmission STRUCTURE is
established (linearity → local-occupation argument, bath-free); the
normalization is one contested, underived number with a named exact
requirement; the gate needs a real-exchange partner whose derivation
is now the strike-bearing next leg.

ELI12: We finally computed the push honestly. Three things came out.
(1) The universe's hum can only LEAN on our melody — and the lean
reads the melody's own warmth, not the room's noise; the scary big
numbers from last round were a formula borrowed from the wrong kind
of touching. (2) If the lean has exactly unit strength, our entire
sky formula — every magic fraction, 1/2, 1/12, 1/720 — pops out at
once; and the sky already measured the strength: it is 1 within 5%.
(3) But leaning can never pass the secret handshake we measured (the
e^{−x} gate) — leaning actually works best in an EMPTY room, backwards
from the sky — so the handshake must be a real trade with the system's
own local cloud, and deriving that trade is the last missing piece.

### Stage 10E — O5-VERTEX: the l=2 near-field exchange channel (2026-08-09; pre-reg 244828a; the strike-bearing R22-cond-4 leg)

**Question.** 10B: the constraint sector's soft vertex is dispersive at
l=0 (interior-constancy ⇒ [H_int, N_s] = 0). 10C: a purely dispersive
ambient coupling cannot reproduce the measured gate e^{−x}. The gate's
only surviving partner is the LOCAL 6Y cloud through the near-field
l=2 sector — derive that channel as a quantized Rabi-capable exchange
vertex, or the pre-registered 6X/9U strike fires.

**Verdict: V-GRAY, gates 5/5 first-run**
([calcs/stage10e_vertex.py](calcs/stage10e_vertex.py)).

STRUCTURE EARNED (exact): the l=0 selection rule NEGATES at l=2 (the
interior-regular Laplace solution is r²Y₂ — non-constant — so an
exterior l=2 fluctuation couples to the internal quadrupole, not the
total energy; the 10B dispersive-forcing rule is l=0-specific by its
own scoping). Both partners are dressing modes (matter modes stay
10B-dead); an aspherical source imprints l=2 on its cloud linearly in
the mode coordinate ⇒ H_int ∝ (a+a†)(c+c†), [H_int, N_s] ≠ 0
(Fock-exact; the l=0 contrast commutes exactly); the resonant part
a†c + ac† IS the 6X exchange operator. The 10B horizon exclusions are
inapplicable axis-by-axis: not gapped (Ω = x_amb·T, soft), not
lever-arm-suppressed (co-located at r_M), not empty (n_amb measured
0.502/6.63). THE INERTIA-CANCELLATION LEMMA (sympy exact): λ =
(η₂/4)·e_s·e_a·√(ωΩ) — the mode inertias cancel, so the magnitude is
closure-INDEPENDENT (10C's E_c = 2ħω is not load-bearing here). THE
CEILING IDENTITY (6/6 digit regression vs the archived 10C G5 table):
λ_ceil = ½√(ωΩ) IS g_close — at perfect geometry the vertex supplies
exactly the κ=1 closure coupling; the measured geometry IS the gap.

MAGNITUDE SHORT: central λ = 1.8e-4..2.2e-3 H (e_a = archived 4K/5S
EFE quadrupole 0.086; e_s binary = 1/(4x_loc²) point-pair multipole,
galaxy 0.1–0.5 banded; η₂ = 0.53–1.00 computed overlaps); even
max-favorable edges reach only 1.8e-2 H vs the demanded [0.072,
1.571] H; Rabi angle ≪ π per Hubble (no saturation); at the galaxy
anchors even PERFECT geometry undershoots the band (ceiling
0.009–0.030 H). Sun contrast row: planet-driven ε₂ ≈ 1.05e-9 ⇒ the
channel is SOURCE-ASPHERICITY-GATED, 8.4 orders dead for the Sun
(the trajectory formulation's Cassini quiet gains a vertex-level
echo; reading grade).

V-WEAK did NOT fire (max-favorable > the 1e-3 H bar); V-DIAGONAL did
NOT fire. PRE-SIGNED ROUND-25 ADJUDICATION: rule the shortfall
STRUCTURAL → the strike fires in adoption (mech 15→8); identify a
live escape axis ({e_a re-identification: static-q vs
fluctuating-sector amplitude; coherent/collective amplification;
coherence-time re-scoping}) → HOLD 15, successor O5-GEOMETRY.
anomaly-real 53 untouched.

PLAIN VERDICT: NEEDS REFINEMENT at the adjudication grade — the
carrier EXISTS as a class (structure exact, first constructive vertex
since 6X) but undersupplies the demanded coupling by 1.5–2.5 orders
at measured geometry; whether that shortfall is fatal or escapable is
now the round's question, pre-signed both ways.

ELI12: We went looking for the doorway the secret handshake could
walk through. We found it — a real door, exactly where the earlier
proofs said no door could be (those proofs only covered round rooms;
this door needs corners, and binaries and galaxy disks have corners
while the Sun doesn't — which is WHY the Sun stays quiet). But the
door is heavy: pushed as hard as the measured numbers allow, it opens
about a hundred times too slowly to carry the handshake. Either
someone is allowed to push harder than we assumed, or the handshake
doesn't go through this door.

### Stage 10D — O5-KAPPA: the joint single-κ fit (2026-08-09; pre-reg d2d9440 + amendments A1 14a64b0 / A2 pre-quote; R24 cond 8)

**Question.** 10C's map ν = 1 + κ·n_BE has ONE free constant, read
two-pole by the archived meters (a₀-lock 1.00±0.05; 4Z hier
translation 1.48). First DIRECT κ-family fit: deep amplitude κ²a₀
and transition constant c₁ = 1−κ/2 in one likelihood (4Z hier
apparatus + lensing; 5M vertical co-read; F1 free / F2
horizon-locked / F3 κ=1 / F4 split over data-defined bins).

**Verdict: K-SPLIT, all gates PASS at run-3**
([calcs/stage10d_kappa.py](calcs/stage10d_kappa.py); amendments A1 =
tol-looped convergence + ridge starts + guard widening after run-1
under-convergence, A2 = seed-budget injection gate after two
converged single-draw misses exposed a bar narrower than the
estimator's own realization scatter — the CORRECTION-#10 class,
caught in-stage; both logged pre-quote, runs 1–2 archived).

THE MEASUREMENT: single κ̂ = 1.503 (boot 16–84: 1.32–1.67; implied
c₁ = 0.248 — the direct fit REPRODUCES the 4Z hier translation) —
but the one-κ form is REJECTED: Δ(F1−F4) = +21.7 (boot median 20.0;
36/40 ≥ 6.2, 33/40 ≥ 11.8; split-null max +8.5). THE DECOMPOSITION:
κ_deep = 1.317, κ_trans = 1.036, κ_tail → 0 (proxy at the guard;
the tail wants faster-than-BE die-off = the 5G sharp-tail vote in
κ-language). READING: the archived "hier κ = 1.48" was the
COMPROMISE of a running κ(x), not a transition property — the
transition's own constant sits at c₁ = 1 − κ_t/2 = 0.48 ≈ ½, i.e.
AT the closure value; the running lives in the deep amplitude and
the tail exponent. Injection budget: truth κ=1 → mean 1.019 (draws
0.87–1.33); truth 1.5 → mean 1.500; the estimator is unbiased with
±0.25 single-draw scatter (bootstrap primary, as pre-registered).

MANDATORY CO-READS: lock row (a₀ = horizon prior) κ = 0.925 Planck /
0.888 SH0ES at +73 lnL vs free — under the temperature lock the
compromise sits BELOW 1; vertical treatment (F3v regression exact
−12152.49): F1v κ = 1.394, split still preferred (+14.4) but the
deep/transition assignment FLIPS (κ_d 1.48 / κ_t 1.61) — the
deep-vs-transition split is treatment-dependent (the 4W
identifiability boundary), while κ_mid ≫ κ_tail → 0 is
treatment-STABLE. Boundary variants: pattern monotone, Δ grows
(+36.5/+48.6). Per the pre-signed map: K-SPLIT → mech HOLDS 15 (the
constant-κ closure is killed at face value; the κ=1 target
re-scopes); anomaly-real 53 untouched.

PLAIN VERDICT: SUCCESS as a measurement (the sharpest in-catalog κ
statement the program owns: κ RUNS — mid-transition at the closure
value, tail to zero, deep contested between treatments); DIFFERENT
PHYSICS for the constant-κ closure reading (one number cannot carry
the family; the running IS the 5G/5T structure seen through the κ
window).

ELI12: We asked the sky to read the lean-strength dial once, with
everything else honest. The sky refused ONE number: in the middle of
the transition it reads almost exactly 1 (the derivation's favorite
value!), but the faint outskirts want the dial near 1.3 and the
bright end wants it near 0 — the dial TURNS as you walk across the
sky. The old "1.48" was the average of a turning dial. A lean whose
strength changes with the room is not a single lean — either the
formula's strength genuinely runs (the gate physics already predicts
the bright-end fall!), or two different knobs are being read as one.

### ROUND 25 (2026-08-09, the arc adjudication; report archived REVIEW-ROUND25-OPUS.md, uncommitted)

The referee reproduced EVERY load-bearing number from an independent
driver (F1/F4/Δ and the bootstrap to the digit; the 10E lemma,
anchors, η₂, Sun row exact), then ruled:

**10D K-SPLIT UPHELD, UNPATCHED HOLE: NO** — with his own 10-draw
split-null (mean Δ = 2.00, max 4.08, never near the 11.8 bar; our
+21.74 = a genuine ~5–10× excess) and a 6-seed × 6-draw injection
replication that VINDICATED amendment A2 (all means in [0.964,
1.057]; run-2's 1.262 was ordinary single-draw noise). His no-lens
refit: the deep arm is SPARC-driven (κ_d 1.306 without lensing).
Conditions adopted: **C1 the decomposition is RETIRED** — the
running DIRECTION reverses between treatments (plain deep-boostier
κ_d 1.317 > κ_t 1.036; vertical transition-boostier 1.479 < 1.614);
only "one κ rejected" is treatment-stable; the "κ_mid ≫ κ_tail"
sentence dropped as near-vacuous (guard-pinned proxy). **C2 the
ridge caveat, mandatory wherever K-SPLIT is cited**: the split lives
on the low-a₀/high-f_ML ridge (a₀ = 0.37× horizon); it re-scopes
the κ=1 target to the deep limit and is NOT evidence against the
closure — the temperature-locked world sits at κ = 0.925 (+73 lnL).
**C3**: the estimator's mild +0.05–0.07 up-bias at the operative
regime annotated (1.503 reads deconvolved ≈ 1.43; no letter change).

**10E V-GRAY UPHELD, UNPATCHED HOLE: YES (moderate)** — two
trap-#12 clauses, BOTH cutting against the mechanism: (a) the
"ceiling identity" was a definitional tautology — the lemma's true
perfect-geometry vertex is (1/4)√(ωΩ) = HALF of g_close (the T2b
"supplies exactly the closure coupling" clause is RETRACTED; the
vertex is ×2 weaker than advertised); (b) the verdict band was
mis-specified — galaxy g_close (0.009–0.030 H) sits BELOW the 10A
band floor by 10C's own printed table, so a requirement the
galaxies clear dispersively cannot be the galaxy exchange bar.
Scored against g_close: λ_max reaches 0.88 of closure at the
deep-galaxy anchor; e_a ≈ 1 closes it centrally (needed: 1.02).

**THE PRE-SIGNED ADJUDICATION: (B) ESCAPE-AXIS LIVE — NO STRIKE.**
The successor **O5-GEOMETRY** (strike-bearing at the next round),
two axes: (1) derive the ambient cloud's FLUCTUATING l=2 amplitude
— the exchange-active object; the static EFE |q| = 0.086 is a
different object — and score λ against g_close; (2) derive the
exchange leg's actual requirement (6X dynamical saturation ~H vs
the 6U KMS state-statistic reading, which needs only
coupling-existence + equilibration). **PRE-STATED KILL (his
words): if the fluctuating e_a stays ~0.086 AND the requirement is
band-like saturation, the R22-cond-4 strike fires (15→8) at that
round.**

Verification addendum
([calcs/stage10de_addendum.py](calcs/stage10de_addendum.py), GA-1..
GA-7 ALL CONFIRMED): his ceiling algebra exact; his λ/g_close and
e_a-needed tables to ±0.002; his seed-42 injection row matched
DRAW-BY-DRAW to three decimals (max |diff| = 0.000); his null
characterization reproduced on 6 independent draws (mean 2.47, max
5.06); his no-lens numbers to 0.02. First round of the arc where
every reviewer number came back exact on first verification.
CREDENCE (mechanical, pre-signed): bath-mechanism conditional
HOLDS 15 (SIX O5 rounds, zero strikes, zero rises); anomaly-real
53 untouched all arc.

PLAIN VERDICT: the round did its job in both directions — 10D
hardened (letter clean, decomposition honestly demoted), 10E
reframed from "structurally short" to "input-limited on an
un-derived amplitude and a mis-specified requirement," and the
strike now hangs on ONE named computation with a pre-stated kill.

ELI12: The referee checked every number and agreed with both
results — then caught two of our own labels. The "door opens
exactly as hard as needed" claim was us comparing a formula to
itself (the door at best pushes HALF as hard as needed), and we had
graded the galaxy door against the wrong difficulty chart (the
galaxies' own chart is easier — and on THAT chart the door gets to
88% with honest numbers, 100% if the wobbling part of the room's
bulge is as big as bulges usually are). So: no execution today.
One more calculation — how big the wobble really is, and what the
handshake actually requires — decides it, and the kill condition
is signed in advance.

### Stage 10F — O5-GEOMETRY: the fluctuating amplitude and the exchange requirement (2026-08-09; pre-reg f78ff64; the ROUND-25 strike-bearing successor)

**Question.** The R25 kill, pre-stated: fluctuating e_a stays ~0.086
AND the exchange leg requires band-like saturation ⇒ strike 15→8.
Plus the hardened symmetric form (no reading works ⇒ strike).

**Verdict: G-CLOSED, gates 6/6 first-run**
([calcs/stage10f_geometry.py](calcs/stage10f_geometry.py)) — **both
kill clauses FAIL; the strike does not fire at stage level.**

AXIS 1 (the amplitude): THE SECTOR-COORDINATE THEOREM (GF-1, exact) —
the axisymmetric external drive pins ONLY the m=0 component of the
cloud's five l=2 coordinates; the measured static q = 0.086 is that
one DC amplitude; the four m≠0 components have zero static value and
are free dynamical coordinates (same license as the occupied l=0
mode). THE BRIGHT-MODE PROJECTION (GF-2, exact ≤ 5e-16): the
system's quadrupole tensor selects ONE combination (9W); four dark;
M=1 preserved in the exchange channel. ⇒ the exchange partner is the
l=2 bright mode, amplitude ratio 1 BY CONSTRUCTION — the static |q|
is a different object (the R25 Axis-1 question answered
structurally). Magnitude closed form (GF-3, sympy exact):
λ/g_close = (η₂/2)·e_s·√(Ω₂/Ω) — anchor-frequency-independent; Ω₂
banded [√q, 1]·Ω (disclosed): deep-galaxy anchor reaches 0.53–0.98
of g_close (enhanced band), binaries 0.2–0.4.

AXIS 2 (the requirement): THE MAP (6X Hamiltonian verbatim — GF-4
regression reproduced all six archived 10C G7 values to 1e-5 — +
local thermal dissipators, 60 sparse-Liouvillian steady states over
λ/γ ∈ [0.03, 30] × δ ∈ [0, 10λ] at n̄ = 0.502/2.0 + 6.63 extension):
**the ambient KMS gate ratios and the lending-line detailed balance
hold to ≤ 0.8% at EVERY cell — 32/32 R-EQ-region cells clean (GF-5
PASS)** — the borrowed-configuration weight is a STATE property (6U
detailed balance), not a transfer outcome. Requirement table:
R-DYN (saturation) needs λ ≥ 0.5 H, derived λ = [6e-4, 2.6e-2] H —
short 1.3–2.9 orders, NOT satisfiable; R-ADIA fails a fortiori;
**R-EQ (bath-maintained statistics) = λ ≠ 0 + the standing
thermalization license + soft-band δ — SATISFIED at every anchor.**
THE RECONCILIATION: the 10E "shortfall" scored the exchange leg
against the DISPERSIVE channel's requirement (g_close); each channel
owes its own — dispersive amplitude at g_close (10C, κ=1-conditional),
exchange gate at λ = derived/nonzero/soft-band.

**O3 SELF-CORRECTION (caught pre-book, pre-round; the O3 row was
report-grade, no letter clause rests on it):** the printed exchange
coherence |⟨a†b⟩| ~ 1e-17 is NUMERICAL ZERO — an equal-temperature
steady state carries no standing first-order exchange coherence
(detailed balance; no flux between equal-T reservoirs); the "linear
scaling ratio 2.00" was λ-proportional roundoff, not a witness. The
corrected equilibrium statement: first-order coherence vanishes; the
vertex's role lives in the LOOP AMPLITUDE — which is exactly the
O5-INTERPLAY successor object, not this stage's claim. Flagged to
ROUND 26 explicitly, with the companion honesty question: is the map
PROTECTED-BY-CONSTRUCTION at equal temperatures (GF-5 shows R-EQ
sufficiency, not discrimination)?

GF-8 THE INTERPLAY ROW (report-grade): the exchange back-shift on
the dictionary line λ²/(ωδ) reaches 2–22% at the binary x_loc=0.5
band edge (galaxy ≤ 1%) — potentially measurable back-reaction, NOT
budgeted here ⇒ **THE NAMED SUCCESSOR O5-INTERPLAY: one arithmetic
carrying the 6H dispersive shares + the exchange loop's gate factor
+ this back-shift against the measured c₁/κ bands.**

KILL EVALUATION (GF-7): clause (a) FAILS (bright mode, ratio 1);
clause (b) FAILS (R-EQ carries the structure at λ ≪ γ); hardened
clause NOT MET. CREDENCE (pre-signed): G-CLOSED + ROUND-26 no-hole
→ mech 15→18; G-KILLED reviewer-affirmed → 15→8; anomaly-real 53
untouched (no sky fits).

PLAIN VERDICT: SUCCESS at the adjudication the round demanded — both
strike clauses fail by computation (the static-q was the wrong
object; the requirement is not saturation) — with the honest seams
named: the interplay arithmetic is the one remaining pillar, and the
map demonstrates sufficiency, not discrimination.

ELI12: Two questions decided whether the door dies. First: is the
room's bulge stuck (only its frozen shape matters) or can it wobble?
We proved the stuck part is just ONE of five directions — the drive
can only pin one — so four wobble freely, and the door's handle
grabs exactly one combination of them at full strength. Second: does
the handshake need to be PUMPED hard (it can't be — the push is 100×
too slow), or is it a property the room's warmth maintains for free?
We built the room in a computer with real thermal walls: the
handshake's bookkeeping holds at every pump strength, including
nearly zero — it's warmth-maintained. So the door lives. What's left
is one piece of arithmetic: adding the door's small back-push on the
melody to the full song and checking the total against what the sky
plays.

### ROUND 26 (2026-08-09, the decisive adjudication; report archived REVIEW-ROUND26-OPUS.md, uncommitted) — **10F G-CLOSED RETRACTED → G-OPEN (adopted)**

The referee reproduced every algebraic number (GF-1/2/3/4/8 + R-DYN
exact, his own tensor basis and Lindblad solver), then ran the four
discriminating computations the stage OMITTED and **inverted Axis 2**:

**D1 (the decisive one):** the map steady state's ABSOLUTE P(sys=2)
sits at the GIBBS value (0.0753 at n̄ = 0.502; 0.1846 at 2.0) — a
factor ~2 below the lending law (0.1671/0.3333) — **at every λ/γ from
0.03 to 300, never moving**. My o2 "lending-line balance" was the
trivial degeneracy ratio; the actual weight is Gibbs. **D2:** o1 = o2
= 0 at λ → 1e-6 — the GF-5 gates pass AT ZERO COUPLING: a positive
letter clause was gated by a test that cannot fail (null-power;
trap-#12's sharpest instance; "its own O3 self-correction was the
tell"). **D3:** removing the system dissipator changes nothing — the
exchange itself thermalizes the system: ANY steady state = Gibbs; the
lending law is intrinsically NON-EQUILIBRIUM (the 6X |1⟩ preparation
+ dephasing). **D4:** the perturbative exchange dressing weight is
LINEAR in m (|me_up|² = 2m, |me_dn|² = m+1 — raw-n, the weighting 10C
G7 already rejected); the gate e^{−x} = P(n≥1) is a TAIL probability
appearing ONLY in the saturated diagonal ensemble. **His crux:
thermalizing the MARGINALS does not confer the JOINT
borrowed-configuration weight — the R-EQ license was over-extended.**
The equal-temperature system dissipator was the load-bearing modeling
choice: it makes the joint Gibbs state the Liouvillian fixed point —
the one substitution that guarantees the 6X result cannot reappear.

Axis 1 ruling: the static-q ≠ fluctuation-amplitude distinction is
VALID (kill clause (a) not cleanly met) but **e_a = 1 is asserted,
not derived** (PERMITTED-grade): GF-1/2 prove the m=0 pinning and
the bright-mode projection, not the mode's spectral
weight/participation; the constraint-slaving question (vs 10B-C6) is
untouched. The Axis-1 residual = the SAME object as the collective-
amplification escape.

**ADJUDICATION (B): HOLD 15 — no rise (hole YES, serious, two), and
NO STRIKE** ("a powerless instrument supports neither the closure nor
the kill"; the collectivity axis is uncomputed — the 9W λ̄ = √(Σλ_k²)
√K-question could lift λ toward saturation; a 7-point strike requires
closing it). His summary of the state: the R-EQ escape as computed is
REFUTED; the requirement side is closed against equilibrium by D1–D4;
**the mechanism holds only on the uncomputed amplitude/collectivity
axis — "input-limited on the l=2 coupling normalization, not shown
structurally dead," on a narrower ledge than the stage claimed.**

**THE SHARPENED SUCCESSOR (O5-COLLECTIVE; strike-bearing at its
round, his condition 6 verbatim):** the mechanism closes iff EITHER
(a) collective amplification is derived and lifts λ to saturation
(λ ≳ H), OR (b) a non-saturated real-exchange process is shown to
imprint e^{−Lx} on the system's DRESSING. **IF BOTH FAIL, THE STRIKE
FIRES (15→8).** Also adopted: GF-8's label upgraded to "unbudgeted
consistency tension at the binary anchor" (the 22% band-edge
back-shift exceeds the ±5% a₀-lock κ band; galaxy-led adopted bands
and the 10D K-SPLIT untouched); any future requirement map must be
NON-EQUILIBRIUM at the physical λ.

Verification addendum
([calcs/stage10f_addendum.py](calcs/stage10f_addendum.py)): GA-1..
GA-6 ALL CONFIRMED — his Gibbs values re-derived exactly (the
0.0744-vs-0.0753 gap = infinite-vs-4-level geometric, explained
digit-exact), the null-power demonstration reproduced (o1 = o2 = 0
at λ/γ = 1e-6), the raw-n matrix elements exact, the 39-Hubble-
period Rabi cycle and the 22%-vs-5% arithmetic confirmed. CREDENCE
(mechanical): bath-mechanism conditional HOLDS 15 (SEVEN O5 rounds,
zero strikes, zero rises); anomaly-real 53 untouched all arc.

PLAIN VERDICT: NEEDS REFINEMENT, honestly earned — the stage's
structural half survives (the static-q two-object distinction; the
bright-mode projection; R-DYN closed as unsatisfiable at single-mode
λ), its closure half is retracted (the map had no power; equilibrium
cannot carry the lending law), and the mechanism question is now
compressed to ONE derivable object — the collectivity of the l=2
bright mode — with the kill signed in advance on both of its exits.

ELI12: The referee agreed the room's bulge can wobble — but caught
that our "thermal walls" test could never have failed: warm walls
keep everything at room temperature BY DEFINITION, so of course the
bookkeeping looked fine; when he measured the actual amount of
lending in that warm room, it was exactly the plain room-temperature
amount — none. A handshake needs someone to actually reach out, and
reaching out takes the strong push we don't have — UNLESS many hands
push together (the one thing nobody has computed yet). That
many-hands calculation is now the whole question, and both of its
possible answers have signed consequences.

### Stage 10G — O5-COLLECTIVE: the participation ceiling (2026-08-09; pre-reg 6f76672; the ROUND-26 strike-bearing successor)

The R26 kill, verbatim: the mechanism closes iff EITHER (a) collective
amplification is derived and lifts λ to saturation (λ ≳ H; the 9W
λ̄ = √(Σλ_k²) √K-question), OR (b) a non-saturated real-exchange
process is shown to imprint e^{−Lx} on the system's DRESSING. Both
fail → 15→8, upon ROUND-27 affirmation.

**FIRED: C-STRIKE-CANDIDATE — all gates green, both kill clauses
evaluated FAILS
([calcs/stage10g_collective.py](calcs/stage10g_collective.py),
data/stage10g_collective.txt).**

**Clause (a) — THE √K QUESTION IS CLOSED BY A SUM-RULE CEILING.** The
chain, each step gated: (1) any partner set couples through the ONE
l=2 port (GF-2 rank-1, banked) → bright-mode λ̄² = Σλ_k² (G1b: the
one-excitation eigenvalue theorem, 1e-16); (2) per-mode 10E lemma
λ_k = (η₂/4)e_s c_k√(ωΩ_k) with c_k the mode's zero-point amplitude
of the ONE fractional-quadrupole field (G3 regression reproduces the
archived 10F table to ≤ 5e-4); (3) the modes superpose in one
observable ⇒ THE BUDGET ⟨δq²⟩ = Σc_k²(2n_k+1) ≤ q_phys² (G1a exact;
q_phys ~ 1 = the harmonic-chart premise, disclosed load-bearing);
(4) THE KMS PIN: n_k = n_BE(Ω_k/T_dS) is measurement-anchored (6U
detailed balance IS the measured gate; the a₀ temperature lock) —
not a modeling choice; (5) ⇒ the budget-to-rate trade has an EXACT
optimum: weight per unit budget = Ωn/(2n+1) = **Ω/(e^{Ω/T}+1) — the
FERMI function** — maximized at x* − 1 = e^{−x*} (x* = 1.27846),
f* = x* − 1 EXACTLY (G2; 20-spectrum LP check, no excess at 1e-9).
**THE CEILING: rate_max = (η₂/4)e_s√(ω(n_s+1)f*T_dS)·q_phys·√N_p —
K NEVER APPEARS.** Participation redistributes the fixed budget;
√K growth would need ⟨δq²⟩ ~ K(2n+1), which the budget forbids. The
thermal occupation CANCELS: every quantum of Rabi enhancement is
bought by budget suppression ("the bath cannot fund both the loan
and the collateral").

Numbers (q_phys = 1, requirement 0.5 H): honest (η₂ central, N_p = 1
= GF-2's own rank-1): binary x=0.5 shortfall 27.6×, deep galaxy
101×; MAX-GRANT (η₂ = 1, all five sector components granted as
ports despite GF-2): tightest anchor 9.4× short, every anchor ≥ 5×.
Robustness: loosest budget (normal-ordered + ZP-cap — licenses
unbounded zero-point amplitude at empty modes, chart-dead a
fortiori) still 4.99×/18.3× (strike-leg bar ≥ 3×); non-KMS partner
rows 2.8–4.0×/10–15× printed REPORT-grade, excluded from the strike
by SCOPING (a non-KMS partner imprints a non-KMS gate and breaks
the temperature lock — even saturated it delivers a gate the sky
did not measure; scoping offered to the round). THE CLOSURE SURFACE:
q_phys ≥ 9.4 (binary) / 34.7 (galaxy) at maximal granting = RMS
≥ 945%/3467% fractional distortion — **THE DILEMMA: there the 10E
linear-order lemma un-derives the vertex itself; clause (a) demands
"derived AND lifts"; beyond the chart nothing is derived. Either
way (a) fails.** Report-grade dividend: the 10F e_a = 1 corner
already implied ⟨δq²⟩ = 2n_amb+1 = 2.0/14.3 (RMS 142%/378%) — the
R25 "closes centrally" corner sat past the chart premise
(O5-INTERPLAY consistency family).

**Clause (b) — the dynamical leg on the 6X engine itself** (verbatim
Hamiltonian; closed |1⟩⊗thermal evolution = non-equilibrium at
physical λ, trap-#17 compliant): the γ=0 machinery reproduces all
six archived 10C G7 lending values (SB-4a) and the TIME-AVERAGED
evolved P2 at strong coupling matches the saturated ensemble to
0.000 rel (SB-4b — the dynamical route does saturate; the
instrument has power both directions, π-row prints P2/sat 1.08–1.33).
At the ceiling's own Rabi angles (θ_c = 2π·rate·t_univ = 0.030–0.32
rad): **P2 reached = 0.07–1.5% of saturated** (bars 1%/5% honest/
max-grant), and the structure probe confirms RAW-N (w(m)/w(1) =
1.000/1.996/2.988/3.977 vs 1/2/3/4 — the D4 weighting 10C G7
rejected), NOT the gate P(n≥1). A non-saturated approach cannot
imprint e^{−Lx}. Sector sweep (G7): (ii) free radiation 34–44
orders dead (10A archived U's); (iii) matter kinematics thermally
EMPTY at 10³–10⁶ decimal orders (galactic orbits sit at 361 H,
x = 2275; the one soft-band matter family was 10B
dictionary-killed); (iv) coherent partners ride the MEASURED static
0.086 (mean-vs-variance split) = the original 2.4–3.4-order
shortfall, Poisson-gate postdiction lean +1.37σ vs thermal +0.56σ
(9U), pure drives carry no occupation statistics (the measured
two-anchor tail split unproducible); (v) no partner = the SB-6 leg.
The partition (bound-cloud/free-field/matter × state character) is
offered to the round for attack.

Pre-signed: C-STRIKE-CANDIDATE + ROUND-27 affirms → mech 15→8
(REQUIRED, the R26 signature); + unpatched hole → HOLD 15, the hole
= successor; C-SURFACE-GRAY → HOLD 15 (successor O5-ANHARM).
anomaly-real 53 UNTOUCHED every branch (no sky fits). ROUND 27 =
the adjudication; the credence does NOT move at stage level.

PLAIN VERDICT (stage-level, pre-round): SUCCESS as an instrument —
the commissioned question ("does collectivity rescue λ?") is
answered with an exact, K-free ceiling and a powered dynamical
test; the answer is NO by 9.4× at maximal granting, and the
non-saturated route is closed at 0.07–1.5% with the wrong
occupation structure. Whether this beheads the mechanism is now
ROUND 27's call — the strike is his to affirm or break.

ELI12: We asked whether many hands pushing together could supply
the strong push one hand lacks. The accounting says no: all the
hands share ONE body — the cloud can only wobble so much in total,
and every bit of "many hands" strength has to be paid out of that
same fixed wobble budget. Written out exactly, the best possible
deal the cloud can strike is a push still ten times too weak — and
pretending the cloud could wobble ten times more than 100% would
break the very math that gave us the push in the first place. We
also ran the actual handshake movie at the best allowed push: after
a whole age of the universe the hands have barely twitched (under
2% of a real handshake), and the twitch has the wrong fingerprint.
Now the referee gets to try to break the accounting.

### ROUND 27 (2026-08-09, THE STRIKE ADJUDICATION; report archived REVIEW-ROUND27-OPUS.md, uncommitted) — **C-STRIKE-CANDIDATE UPHELD; NO HOLE; THE STRIKE FIRES: bath-mechanism conditional 15 → 8**

The referee (fresh, five files + his own code, nothing touched in
the repo) reproduced EVERY stage number from scratch — the Fermi
optimum to 1e-16, the full ceiling table to the digit, the 6X
dynamical rows and the raw-n structure probe, the implied-budget
dividend — **zero mismatches** — then attacked the kill on every
named surface and reported it STRENGTHENED:

**His attack (i) = the round's best work, ADOPTED AS THE PRIMARY
clause-(a) statement (C1):** grant the mechanism its OWN
fluctuation–dissipation budget — q_phys = √(2n_amb+1), the FD value
at the measured T_dS, i.e. the entire fluctuation the e_a = 1
reading ever assumed — stacked on max-grant: **the shortfall is
still 6.67× (binary) / 9.18× (galaxy); closure needs e_a = 6.7/9.2,
impossible for a ratio bounded ≤ 1, or ⟨δq²⟩ = 89/1202 = 45×/84×
the FD value.** The strike therefore does NOT rest on the disputed
q_phys ≲ 1 chart premise — it rests on fluctuation–dissipation at
the measured temperature plus the ratio bound, both solid. The
three closure knobs are each independently pinned: K by the
shared-field budget (his explicit demo: budget-saturated rate
EXACTLY flat over K = 1/8/10000), e_a by the ratio bound, the
single-mode occupation by the a₀ temperature lock (which IS the
measured gate).

**His attack (v) = the sharpest clause-(b) escape, self-built and
self-closed:** dephasing-assisted equilibration (pure N_a dephasing
+ resonant exchange — the one route that could have imprinted the
gate ½P(n≥1) without coherent saturation). His Lindblad table: the
gate is reached only near θ ≈ π (saturation) at EVERY dephasing
rate; stronger dephasing is SLOWER (Zeno); at the ceiling angles
P2/gate = 0.002–0.015. The dephasing route needs the same λ ≳ H.

Also ruled: the non-KMS scoping is a REDUNDANT SECOND DEFENSE, not
load-bearing (C2 — the Λ = 1 H row is itself short of 0.5 H at
0.178 H, and 8.22× short at honest N_p = 1; and for resonant
lending the partner's occupation sets both the rate and the gate —
the same number — so a super-thermal rate boost destroys the
measured e^{−Lx}); the requirement is not load-bearing either (even
a 10%-dressing needs λ = 0.166 H = 3.1× the max-grant rate);
N_p ≤ 5 is the true max; the (n+1) channel is the leak, correctly
excluded; formation-history accumulation = continuous driving =
Gibbs (his D1 re-run); the sector partition holds (constraint-
sector and horizon-microstructure candidates are 10B-excluded as
real-exchange carriers; multi-system Dicke enhances the collective
observable, not the per-system dressing, and is budget-bounded).
Trap audit CLEAN: no null-power gate (the q = 120 counterfactual
and the π-row both flip), no damaging self-regression (G3 relabeled
transcription/convention pin per C3), and the one ungated premise
(q_phys) does not carry the verdict under the FD-primary form.

**RULING: C-STRIKE-CANDIDATE UPHELD. UNPATCHED HOLE IN THE KILL
CASE: NO. ADJUDICATION: AFFIRM — the R26 signature executes.**
His scope clause (C5, adopted): the strike kills the REAL-EXCHANGE
leg — the 6X lending gate, the derived s^L screening, the
two-system tail split have NO microphysical realization at the
derived coupling — while the DISPERSIVE leg (the 10C polaron
theorem, the additive C&T-law selection, the κ work) SURVIVES;
hence 15→8, not a total kill; the RAR fit and every sky measurement
are untouched — the gate becomes MEASURED-BUT-NOT-MICROPHYSICALLY-
DERIVED. The sole surviving conditionality (C4, named): the
affirmation is conditional on the soft-sector ⟨δq²⟩ ≲ O(2n+1);
**O5-ANHARM (derive ⟨δq²⟩ horizon-side) is the only computation
that could ever reopen clause (a), with an ADVERSE PRIOR** —
reopening needs 45–84× FD (RMS ≳ 950%), which would overturn the
10A cloud-continuum-leak closure AND the 9Z gapped-η⁴ result AND
break the temperature lock. "A wish against three existing
results, not a live hole."

Verification addendum
([calcs/stage10g_addendum.py](calcs/stage10g_addendum.py)): GA-1..
GA-7 ALL CONFIRMED — the FD shortfalls/closure budgets digit-exact,
K-flatness exact (his Σc² = 0.49879 traced to n re-derived from
x_amb = 1.0954; his allocation = 0.993 of the Fermi optimum), the
non-KMS restatement, the relaxed-requirement arithmetic, **the full
dephasing Lindblad table reproduced to dmax ≤ 0.005 with the Zeno
ordering**, the R26-D1 Gibbs spot (0.07531 exact), the gate-law
track + coherent z row. One verification-side correction disclosed:
the GA-7 first-run bar (0.006 abs) was tighter than the banked 6X
grade (0.2–2.2%); corrected to rel ≤ 0.025 — no reviewer number
involved.

**CREDENCE (mechanical, the pre-signed cell): bath-mechanism
conditional 15 → 8.** The program's first strike executed through
the full pre-signed machinery since 6N (history: ~20-25 → 15 at 6K,
15 → 8 at 6N, recovery 8 → 12 → 15 at 9T/9W, now 15 → 8 — eight O5
rounds, ONE strike, zero unearned rises). anomaly-real 53 UNTOUCHED
all arc (no sky fits anywhere in 10G or the round).

PLAIN VERDICT: DIFFERENT PHYSICS — the honest kind. The
commissioned hunt ("give it all you got") ran eight derivation
stages and eight adversarial rounds, and the answer is that the
real-exchange leg of the bath mechanism cannot be realized at the
derived coupling: collective participation is exactly
budget-cancelled, and every non-saturated route imprints the wrong
statistics. The dispersive leg — the part that selects the additive
C&T law and carries the κ/c₁ phenomenology — stands. The measured
gate stays measured; what died is one candidate explanation of it.

ELI12: We asked if many hands could push the swing hard enough.
The referee tried every trick on our behalf — he even let the cloud
wobble as much as warm physics ever allows (three-and-a-bit times
MORE than 100%), gave it all five handles and perfect grip — and
the push still came up seven-to-nine times too weak. He also tried
the cleverest cheat: constantly bumping the swing so it forgets its
rhythm, hoping it would drift into the right pattern by itself — it
drifts there only as slowly as the honest push, and bumping harder
makes it slower. So he signed the verdict he himself had drafted
last round: this particular engine cannot be what pushes the swing.
The swing still swings — our measurements of it are untouched — and
the OTHER half of the machine (the part that bends the rulebook
rather than trading pushes) survives intact. Honest score: the
mechanism's chance drops from 15 to 8 out of 100.

## ROUND 28 (2026-08-09): the first Paper-3 referee round — MAJOR REVISION, adopted in full the same day (draft 0.2 → 0.3); zero computed errors; the 1.449-vs-√3 "discrepancy" resolved as one theorem in two unit conventions

The round was run as a fresh journal-referee session (OJAp frame, no
program memory, readability an explicit dimension per the author's
standing request), against draft 0.2 + STYLE.md + the companions.
Report archived REVIEW-ROUND28-OPUS.md (uncommitted); every referee
number independently re-verified in calcs/round28_addendum.py
(GA-1..GA-8 ALL CONFIRMED, data/round28_addendum.txt) per the
standing verify-reviewer-math rule.

THE HEADLINE FINDING (his M1): he could not reproduce the abstract's
"every bound orbit satisfies Ω/H > 1.449" — his own SdS computation
put the outermost stable circular orbit at M/r³ = 4H² ⇒ Ω/H = √3 =
1.732. VERIFICATION RESOLVED IT AS A CONVENTION SPLIT, both numbers
correct: the 9Q chain (κ_r² = GM/r³ − (4/3)Λc² ⇒ Ω² > Λc² exact)
quotes the bound at the PHYSICAL Λ = 3Ω_Λ H₀²/c², giving √(3Ω_Λ)H₀ =
1.449 H₀; his H is the pure de Sitter rate (Λ = 3H²), giving √3 H_Λ;
√3·√0.7 = 1.449 exactly, and his OSCO surface is the κ_r² = 0 surface.
The paper's defect was real anyway: "every bound orbit" overstated the
proven class (stable CIRCULAR orbits) and no convention was stated —
§3.2 is now derivation-visible with Ω_orb and both unit forms.

ADOPTED IN FULL (sixteen findings): M2 the weight-spine unification —
the new §2 definitions paragraph writes p = ½ + r·s²/2 once, labels
the tail postdictions as the r = ½ (pure-dispersive, no-averaging)
slice, quotes the 9V fit r̂ = 0.34 ± 0.19 beside it (p_gal 0.688 →
0.628, both inside the band; his arithmetic confirmed), and maps the
three distinct ½/¾ constants (lending prefactor / amplitude ceiling /
tail ceiling ¾ = the r = ½, s² → 1 corner); Figures 1/4 captions now
say the sky points are input occupations on the model curve, NOT
independent weight measurements (m3). M3 the overclaim — "uniquely
selects"/"no longer a choice" now carries the in-reading antecedent,
the characterized two-member family, and the κ-freedom clause;
exhaustiveness language softened. M4 — the acknowledgments now state
plainly that the adversarial reviews are AI referee sessions, not
human peer review. M5 — Ω_orb/Ω(ambient gap)/x_amb disambiguated.
Minors: four orphan references cited in place (McGaugh+16, Milgrom
1983, Chae+21, Desmond+24); §1 band phrasing aligned to 0.26–0.45;
κ = 1 ⟺ E_c = 2ħω labeled a program result; "within 0.001"; abstract
248 words; multipole l-signposting; six jargon glosses; reviewer
attribution dropped; "equal prominence" tempered; the 0.15 exclusion
grade anchored. His numerical-check table matches the banked values
digit for digit (Fermi 1.278465; 0.6885/0.5279; ln 2; 0.0049/0.0525;
6.7²/9.2² = 44.9/84.6 ≈ the banked 45/84 variance factors — his
consistency identity, now GA-6).

CREDENCE: none moved (paper round; anomaly-real 53, mech conditional
8 untouched). P3 status: draft 0.3, first referee round applied; the
next P3 gate is the author's (second round, or the circulation queue).

PLAIN VERDICT: SUCCESS — the round found zero wrong numbers, one
genuinely unclear headline number (now derivation-visible), and the
paper's biggest legibility defect (four weights, one relation, never
written); all sixteen findings adopted same-day with independent
verification of every referee claim.

ELI12: We sent our physics paper to a fresh referee who tried to break
every number in it. He couldn't break any — but he caught us using two
different names for the same speed limit without saying which clock we
used (both were right; we now show the clock), caught us juggling four
look-alike "probability weights" without writing the one formula that
connects them (now written), and made us say in plain letters that our
tough referees so far have been AIs, not humans. The paper got harder
to misread and nothing in it got weaker.

## 10H — THE Z-LADDER CONTEST (2026-08-09, pre-reg 546641f + A1 9f9b4ee + A2 b6452dd + A3 82ad2e1)

The booked MUSE-DARK successor executed: independent reconstruction of the
Ciocan+ 2026 RAR-evolution measurement from its own public release
(dark-matter.osu-lyon.fr; 9 catalogs + per-galaxy galpak DC14 products
fetched, calcs/fetch_musedark3_release.py), aiming leg A (lock-vs-linear on
their tracks) and leg B (lens-bias injections: does their fixed-form
pipeline inflate locked-pair truth 1.02 → their measured 1.20–1.59?). Both
legs sat behind wiring gates set at their published numbers, with the
recipe-ladder firewall (ladder keys only on wiring targets; no contest
quantity computed before a primary recipe locks).

**Verdict: H-FEASIBILITY-LIMITED (the pre-registered letter). No recipe
passed wiring; no contest fit and no injection ever ran; the firewall held
through all five runs.** What the attempt measured:

- **Census:** their N = 79 is unreproducible from the release — the
  quantitative cuts (z, M*, DC14 membership, evidence) give N = 85 under
  every released convention; the 6-galaxy gap is Paper I's *visual*
  "regular" flag, which no released column encodes (A1).
- **Level (ballpark-grade — relabeled by ROUND 29):** the elevated z∼1
  scale reproduces — after the amendment-4 integrator fix the
  Paper-I-pinned recipes PIN20-E1/E2 sit IN the gate band (a₀ =
  2.34/2.51, σ_int 0.218/0.230; run 6), alongside the guessed-recipe
  rows (2.62/0.220 in-band; 2.15/0.190 at 0.01 below — superseded
  recipes, provenance flagged) vs their 2.38 +0.12/−0.10. The
  offset-vs-f_DM match is DC14 self-consistency, not independent
  corroboration (ROUND-29 Finding 5). The headline OFFSET is in the
  released data.
- **Slope: does not survive to pipeline grade (ROUND-29 softened
  form).** Across the recipe space (16 guessed + 6 Paper-I-pinned,
  re-run corrected in run 6), no recipe passed the z-fit gate; the
  level-passing recipes fit a₁ ≤ 0 or a sign that FLIPS with the
  error model (corrected pinned: E1 −0.90/−1.00/−0.60 vs E2
  +1.03/+0.31/+1.84) — the per-point-errors blocker demonstrated
  in-ladder. The offset material is weak and recipe-dependent
  (bin4−bin1 = +0.01…+0.14 dex across recipes/statistics vs the
  ~+0.067 their trend needs), and the ROUND-29 per-galaxy route (the
  paper's own App-E cross-check) shows a positive NON-significant
  lean (referee +0.65, 90% CI [−0.75, +2.68]; our implementation
  +0.22 [−0.93, +1.34]; both CIs contain 0 and their 1.42).
- **The blockers are all z-leveraged and all unreleased:** per-point
  MCMC error profiles (weighting is slope-relevant — the in-ladder
  sign flip above is the demonstration); M_HI / gas-disk extent /
  profile (the exact Casertano solve of a constant-Σ disk nearly
  cancels interior — force 0.6–25% of Mestel, verified — vs the
  v ∝ √(Σr) scaling; ROUND-29 reframe: this ×2/curvature space is
  the STAGE'S own convention variants, not a contradiction in the
  paper, whose methods name the solve); a release version seam
  (folder Vvir 160.6 vs catalog 126.3 at ID3); the 4-galaxy bulge
  list. [Correction: run-5's "pinned recipes drive a₀ to 8–10" was
  OUR ring-integrator π-bug (amendment 4, ROUND-29 Finding 1), not
  the release's gas physics.]
- **Lock-side dividend (no reconstruction needed; ROUND-29 phrasing):**
  at their own Planck-2015 cosmology the lock intercept cH₀/2π =
  1.047e-10 lands ON their MOND-framework row (1.03±0.05, 0.7σ) and
  their per-galaxy best-evidence row (1.05±0.05, 0.1σ), and 2.3σ —
  just outside the 95% CI — from their DM headline row (1.00±0.04;
  all 95% CIs). The lock slope linearized over their window (1.024)
  sits 3.4σ below their MOND-row a₁ = 1.20±0.10 (95%) at face value —
  and the form-coupling bias that would arbitrate that comparison
  (leg B) is exactly what the release cannot support.
- **Trap #19 (standing):** run 3's frozen optimizer "reproduced" their
  Eq. 4 exactly because the fit starts WERE the published target values.
  Wiring-gate fits must never seed at the wiring target.

Credence: map cell FEASIBILITY → anomaly-real HOLD 53; mech 8 untouched.
Successors: the author-contact route (per-point errors + M_HI + flags +
MOND-track products close every blocker; Desmond is a co-author and
already on the outreach list — circulation-adjacent), galpak_dark source
read, Vărăşteanu+25 MIGHTEE (z<0.08, a₁ = 4.47±1.88), DR4-era. Ledger
ext-10h-musedark-recon (206 rows, six gates PASS, 266 tokens).

**Plain verdict: NEEDS DIFFERENT DATA.** The instrument was built, gated,
and honestly refused to fire: the public release reproduces their
headline *offset* but underdetermines their headline *slope* — the
deciding columns exist only on the authors' disks. The lock-relevant
finding that needed no reconstruction at all: their own MOND-framework
intercept lands on cH₀/2π to 0.7σ.

*ELI12: We tried to rebuild their "the gravity knob grows with cosmic
time" measurement using only the parts they published. We could rebuild
the knob's SIZE (it matches), but not its GROWTH — three small
ingredients that control the growth (error bars per point, how far the
invisible gas spreads, which galaxies they hand-picked as "clean") were
never published. So the referee gates said "stop, don't guess." One free
gift: at their own choice of cosmic clock, the knob's starting value in
their MOND fit lands almost exactly where our formula says it must.*

## ROUND 29 (2026-08-09, the 10H referee round — adopted in full)

Fresh adversarial session against the 10H letter and all five claims;
he rebuilt the tracks from the release himself and re-ran the
estimator, the wiring fit, the gas physics, and the lock arithmetic.
Ruling: **UNPATCHED HOLE NO — H-FEASIBILITY-LIMITED AFFIRMED**, with
one major catch of ours and two claim softenings, all verified in
calcs/stage10h_addendum.py (GA-1/3/4/5/6) + run 6 (= GA-2) before
adoption:

- **His Finding 1 (the catch): a factor-of-π normalization bug in the
  stage's ring-force integrator** — the shipped prefactor
  double-divided by π; old/exact v² = 1/π against the exact Freeman
  disk (we verified: 0.316 at every clean radius). It crippled the six
  Paper-I-pinned recipes; run-5's "a₀ → 8–10" was OUR numerics, not
  the release's gas physics. Amendment 4 (38c7907): prefactor fixed,
  NEW gate G10H-9 (integrator vs exact Freeman ≤ 5%; passes at 0.008)
  — the gate nobody had thought to write — and the trap-#19 seed
  removed from the code.
- **His Finding 2 (exculpatory, we reproduce): the letter is robust to
  the fix.** Run 6, corrected ladder, unchanged gates: PIN20-E1/E2 now
  pass the LEVEL gate (a₀ = 2.34/2.51) and every recipe still fails
  the SLOPE gate, with the fitted sign flipping between error models
  (E1 −0.90/−1.00/−0.60; E2 +1.03/+0.31/+1.84) — the unreleased
  per-point errors are slope-deciding, demonstrated in-ladder.
- **His Finding 3 (the route we never tried): the per-galaxy
  regression** (their own App-E cross-check) gives a positive
  NON-significant lean (his Theil-Sen +0.65, 90% CI [−0.75, +2.68];
  our implementation +0.22 [−0.93, +1.34]) — consistent with their
  own weak per-galaxy slope 1.42 +0.94/−0.89 and with zero. No
  usable-grade recovery; "never materializes" softened to "does not
  survive to pipeline grade."
- Claim corrections adopted: C2 relabeled (level = ballpark-grade;
  the fDM-offset match is DC14 self-consistency); C3's "+0.095"
  replaced by the honest range (+0.01…+0.14, recipe-dependent); C4's
  "inside all three rows" corrected (ON the MOND and per-galaxy rows;
  2.3σ — outside 95% — from the DM row); the gas ×2/curvature clause
  reframed as our convention-variant space, not a paper-internal
  contradiction; the dump-artifact count corrected 27 → 11 (we
  conflated window-outliers with dumps).
- Amendment-integrity audit (his A1): git chain monotonic, amendments
  not outcome-tuned, the A2 window cut removes zero points beyond the
  v-floor in clean recipes; trap #19 verified real (the frozen run-3
  fits returned their seed = the published target) and now enforced
  in code.

No credence at stake in any branch (feasibility cell HOLD 53 / mech 8
already executed). Report archived REVIEW-ROUND29-OPUS.md (never
commit). Ledger row ext-10h-musedark-recon patched in place with the
adoption note; six gates re-audited PASS.

**Plain verdict: SUCCESS** — the round did exactly what rounds are
for: the stage's letter survived an independent rebuild, and the one
real defect it found was ours, found fixable, and now gated forever.

*ELI12: We asked a fresh robot referee to try to break our "we
couldn't rebuild their growth measurement" conclusion. He rebuilt
everything his own way and agreed — but he also caught a real bug in
our gravity calculator (we divided by π twice), which had made one
family of our attempts look worse than they were. We fixed it, added
a permanent test so it can never happen again, re-ran everything, and
the conclusion came out the same — now for the right reasons.*

## Stage 10I (2026-08-09): THE ISOLATION LADDER — the 9V void-channel reopen at LV grade. I-FEASIBILITY-LIMITED (run 1 = the letter)

Pre-reg 630516e (bars, letters, credence map, census bars set BLIND,
sky-blind injection firewall); script
[calcs/stage10i_isoladder.py](calcs/stage10i_isoladder.py) →
data/stage10i_isoladder.txt. Data: SPARC × UNGC (Karachentsev+13,
VizieR J/AJ/145/101, fetched by
[calcs/fetch_flynn_corpus.py](calcs/fetch_flynn_corpus.py) alongside
the Flynn 2026 HI corpus, whose primary read confirmed it adds no
decomposed kinematics beyond SPARC — LITTLE THINGS rows are V_obs-only,
WALLABY has no decomposition/errors; trap-#6 pattern, routed around).

The stage STOPPED at the pre-registered G10I-2b gate, first run — the
gate wired to fail on a null axis (trap-#16 compliance) fired:

- CENSUS (G10I-0 PASS at the blind bars 20/6/3; ROUND-30: quote as
  UNDERCOUNTS, not counts): 31/149 SPARC matched (name + distance);
  iso (Ti1 < 0) 13; **iso transition-crossers (≥3 pts y ≥ 0.8) = 3**;
  deep-iso (Ti1 < −1) 4. The matcher misses NGC-catalogued galaxies
  carried under UGC designations in SPARC — exemplar **UGC05721 =
  NGC3274 (Ti1 = −1.2, deep-isolated, in Chae, distance gate
  passes)**; 22/28 unmatched at D ≤ 12 Mpc are UGC-designated. The
  referee force-fixed the alias and the calibration gate fails HARDER
  (ρ 0.000 → −0.126, addendum-verified) — the letter is
  alias-robust. The void-channel sample at public LV grade is THIN
  regardless of any gate question; a real alias table (SIMBAD/NED
  cross-IDs) is the free fix for every future LV cross-match.
- G10I-2a PASS: our Θ₁ re-derivation from UNGC positions+K-masses
  regresses on the published Ti1 at slope 0.956 (offset −10.45,
  N = 31) — parse and 3D geometry validated.
- **THE NULL-AXIS FINDING (G10I-2b STOP; ROUND-30 rescope adopted):
  the point-source neighbor sum is floor-dominated exactly at the
  isolated end** — iso-stratum neighbor terms ≤ 1.7% of the 0.01
  floor (e_N percentiles 16/50/84 = 0.0100/0.0101/0.0102 — the floor
  IS the number); **the meter is ASYMMETRIC: at the group end it has
  real signal (NGC2976, Ti1 = +2.9, neighbor term 109.9% of the
  floor)** — null precisely where the science needs it, powered where
  it doesn't. The iso gate range collapses to [0.8177, 0.8186]; the
  r = ½ predicted stratum contrast Δp = 0.0003. The "open-gate
  extension" this instrument was built to measure would have been
  100% floor-assumption — the pre-registered STOP is the honest
  letter.
- THE DESIGN-LEVEL RESULT (ROUND-30 corrected wording — two facts,
  stated separately, neither being "isolation does not open the
  gate," which the stage did NOT establish): (i) the LV
  tidal-isolation axis, read through a point-source neighbor sum, has
  NO dynamic range at the operative e_N level (floor/LSS-dominated);
  (ii) LV galaxies — isolated and not — sit near-open at Chae's
  measured e_N ≈ 0.005–0.009 a₀ (iso-overlap gates s² = 0.825–0.869,
  the same band as non-iso), so LV tidal isolation does not reach the
  void regime (g → 1 needs e_N ≪ 0.005). The void CHANNEL is
  untested, not shown infeasible: a density-field-selected sample
  could still show low e_N. A real void instrument needs
  density-field selection (2M++/Cosmicflows-grade e_N, not a
  neighbor sum) WITH transition-crossing kinematics (CAVITY-class
  when kinematic products release; WALLABY × void catalogs; DR4-era).
  The Chae-vs-UNGC meter cross-validation did NOT complete (no
  variance on the UNGC side) — do not cite it as a two-meter
  calibration.

Firewall record: no iso-gate likelihood, no injection, no sky read was
ever computed. Credence: pre-signed cell — HOLD mech-conditional 8,
HOLD anomaly-real 53. PREDICTIONS P1/P2 annotated (no status flip).

**Plain verdict: NEEDS DIFFERENT DATA** — with the ROUND-30 meaning
relabel: the letter documents INSTRUMENT infeasibility (this
point-source meter, at this catalog grade), not channel death — the
void channel is untested. The 9V reopen clause is executed and closed
for THIS route with the blocking physics quantified (the gate
variable lives in the density field, which no public kinematic sample
tags), and the census (≥3 isolated transition-crossers, an
undercount-floor) prices the sample any successor must beat. R30's
inverted dividend: the meter's powered end is the GROUP end — a
"does the gate CLOSE (p → ½) toward dense environments" test is the
axis this data actually supports (named successor).

*ELI12: Our theory says lonely galaxies should feel less of the
universe's background hum. We tried to build a lonely-galaxy detector
using a catalog of which galaxies have close neighbors. The check we
wrote in advance caught the problem: our ruler couldn't see any
difference among the loners, because for them the hum comes from huge
far-away stuff (superclusters), not close neighbors. Our referee
confirmed the ruler DOES work near big crowds (one galaxy next to M81
reads double!) — just not where we needed it. So "no neighbors"
catalogs can't find the quiet galaxies; we need maps of the
universe's giant empty regions with the right rotation measurements —
and nobody has published those yet. We wrote down exactly what's
missing.*

## Stage 10J (2026-08-09): O5-KAPPA-DEEP — the dictionary coefficient at the deep limit. K-CONDITIONAL, gates 4/4 + A1

Pre-reg 865c01e + amendment A1 1176e7c (post-run-1 self-catch, logged
pre-quote, run 1 preserved as data/stage10j_kappadeep_run1.txt);
script [calcs/stage10j_kappadeep.py](calcs/stage10j_kappadeep.py) →
data/stage10j_kappadeep.txt. The R24-cond-6 / 10D-C2 target: derive
(or bound) κ = 4g²/(Ωω) at the deep limit, where the temperature-
locked world reads κ = 0.925. Dispersive side only (10G strike scope
respected — no lending-law input anywhere).

**EARNED UNCONDITIONALLY (sympy-exact + numerics):**
- **THE THRESHOLD CATALOGUE** (general σ, A1): the dressed ladder's
  spacing zero sits at κ_fold(n) = 4/(2n+1+2σ). At σ = ½ — and σ = ½
  is DERIVED, not chosen: the 10B vertex couples the system
  HAMILTONIAN (H_int = (Ĥ/c²)δΦ), so the zero-point gravitates
  (premise P-half) — the thresholds are {2, 1, 2/3, …} and **κ = 1
  is EXACTLY the two-rung degeneracy E(2) = E(1)**: at the closure
  value the dressed mode's second quantum becomes free. ROUND-30 J-3
  scope adopted: this is exact and not a σ-convention AT THE VERTEX;
  the identification of the threshold-κ with the RAR-MEASURED κ
  remains seam-conditional (the declared dictionary).
- **THE SUSCEPTIBILITY RESTATEMENT** (ROUND-30 J-2 relabel adopted —
  the same two-rung degeneracy seen as a pole, NOT an independent
  result): the static susceptibility of the system coordinate in the
  one-quantum dressed state diverges as κ → 1⁻ with finite
  Franck–Condon prefactor — lim (1−κ)·χ(1) = 4e^(−d²)/ω exactly
  (χ(1) = 10.4/30.2/293/2927 at κ = 0.5/0.9/0.99/0.999; limit
  verified 0.003%; R30 caveat: the κ = 0.5 value is inflated by an
  unrelated |1,0⟩↔|0,1⟩ resonance at exactly κ = 0.4 where
  −D(0)+Ω = 0 — the κ → 1 limit itself is clean). Divergence
  persists under a thermal ambient (diagonal-FC weight 0.354 > 0 at
  n̄ = 2) and the spacing is bath-occupation-free (10C T1). Its
  identification with the additive dictionary's compliance is a
  READING, not derived.
- **THE IDENTITY WEB** (six exact equivalences): κ = 1 ⟺ g = ½√(Ωω)
  ⟺ **g² = E_zp,s·E_zp,b (the coupling energy is the geometric mean
  of the two zero-point energies)** ⟺ displacement energy per unit
  charge² = ω/4 ⟺ φ_req = ½√(Ω/ω) ⟺ E_c = 2ħω ⟺ Δ(1) = 0. The 10C
  "TWO QUANTA" identity resurfaces as the capacity-2 ladder.
- **THE 9W FUNNEL + THE D STATEMENT**: κ_coll = 4λ̄²/(Ωω) with λ̄² =
  Σλ_k² — K-invariant exactly; every multimode route reduces to one
  collective number, and λ̄² is exactly where the 10A-displaced
  D-normalization lives. κ = 1 is derivable WITHOUT D-provenance iff
  a structural principle ties λ̄² to Ωω/4 directly.

**CONDITIONAL (the premise-audit letter, trap-#12 grammar working as
designed — the letter could not print stronger):**
- **THE TWO-RUNG CONDITION (ROUND-30 J-1 CORRECTION ADOPTED — this is
  NOT a bound):** at the two-rung ordering condition D(1) ≥ 0, κ ≤ 1
  — **a rung-count-dependent LOCAL statement**: the same logic at
  three rungs (D(2) ≥ 0) gives κ ≤ 2/3, which would EXCLUDE the
  locked world (at κ = 0.925, D(2) = −0.3875 — the ladder is already
  inverted at n = 2→3, the stage's own fold row's content), and the
  literal-ladder premise fails for the deep sky (⟨n⟩ ~ 1/x ≫ 1) at
  every κ. **The original consistency-row sentence "the spectral
  bound and the temperature lock pick the SAME world" is RETRACTED:
  the coincidence κ_lock = 0.925 < 1 is consistency at the chosen
  truncation, not a spectral bound picking the locked world.**
  P-lit relabeled in the audit: ASSUMED (two-rung only; refuted
  beyond two rungs / in the deep limit by the fold row). The choice
  of exactly two rungs was outcome-determining — it lands on the
  identity-web value by construction of stopping there.
- Saturation κ = 1 exactly rides a second assumption, P-crit (the
  physical coupling saturates the two-rung/compliance condition).
- FOLD OBSERVATION (report-grade): at any fixed κ > 0 the literal
  ladder is non-monotone at n ≥ 2/κ while the deep sky occupies
  n ~ 1/x ≫ 1 — the fixed-coefficient dictionary must be the leading
  form of a self-consistent structure (no direction claim; the 10D
  running direction is treatment-unstable per R25-C1).

LETTER: **K-CONDITIONAL** (mechanical: the chain carries P-lit
ASSUMED; K-BOUNDED requires zero — and R30 J-1 makes P-lit WEAKER,
so no upgrade path exists from this construction). Credence:
pre-signed HOLD mech-conditional 8; anomaly-real 53 untouched (no
sky fits; every κ quoted is an archived meter). SUCCESSORS (R30
section 5 adopted): (1) **THE SELF-CONSISTENT OCCUPATION BOUND** —
replace the rung-count heuristic with the real object: solve the
self-consistent fixed point (ν = 1 + κn_BE feeding back into the
pull; n̄ set by the dressed spacing at the physical T_dS) and ask
whether κ → 1 is where self-consistency first fails — the successor
that could move K-CONDITIONAL toward K-BOUNDED; (2) the level
crossing E(2) = E(1) at κ = 1 predicts anomalously large
low-frequency response near the deep limit under static tidal
perturbation (hybridization) — a potential DR4-era signature,
reading-grade; (3) pin the 10D meters' κ CONVENTION (bare-ω vs
renormalized-ω) before any ridge-vs-lock adjudication is called
well-posed; (4) the P-scale computation (9Z machinery — the known D
wall).

**Plain verdict: SUCCESS at the conditional grade it was designed
for, with the round's correction absorbed** — the κ = 1 target now
has an exact spectral meaning at the vertex (the EP-anchored
two-rung degeneracy), a candidate physical principle (saturation of
its own softness), an honest statement that the "bound" was
rung-scope-conditional (the same-world sentence retracted), and a
sharpened successor (the self-consistent fixed point) that would
make it a real bound.

*ELI12: Our sky formula has one knob, κ, and the sky reads it as
almost exactly 1. We asked WHY 1. We proved: at exactly 1, the mode
being tuned reaches the point where adding a second quantum of wobble
costs nothing — its springiness to a probe becomes infinite. And the
reason the counting starts at ½ (which is what makes 1 the special
value) is Einstein's rule that even vacuum energy gravitates. But our
referee caught the soft spot in the next step: "κ must stay below 1"
depends on how many rungs of the energy ladder you demand stay in
order — two rungs gives 1, three rungs gives ⅔ — so the pretty
agreement with the sky's 0.93 is not yet a real bound, and we
retracted that sentence. The real next move is written down: solve
the full self-adjusting version of the ladder and see where IT
breaks. If that lands at 1, the derivation is real.*

### ROUND 30 (2026-08-09, the 10I+10J referee round; report REVIEW-ROUND30-OPUS.md, uncommitted): both letters AFFIRMED; 10J hole YES scoped to the consistency row; ALL EIGHT conditions adopted

Fresh adversarial session vs both stages. **Every load-bearing number
in both stages reproduced exactly** (his own scripts; zero arithmetic
errors anywhere) — and the verification ran BOTH directions: the
blind half of [calcs/round30_addendum.py](calcs/round30_addendum.py)
was committed BEFORE his report existed (87a4676; χ(1) by full
eigenvectors to 1e-5, thresholds by blind spectral bisection,
Laguerre sum, Θ-slope by independent regression r = 0.991,
neighbor-scale arithmetic), and the post-report half re-verified
every NEW number he produced, digit-exact: D(2) at κ = 0.925 =
−0.3875; the three-rung condition κ ≤ 2/3; the κ = 0.4 downward-term
resonance (−D(0)+Ω = 0.0000 exactly); NGC2976 neighbor term 1.099e−2
= 109.9% of floor; the UGC05721 = NGC3274 alias row (Ti1 = −1.2,
distance gate 1.80 ≤ 2.79 PASS); the forced-alias Spearman −0.126;
the iso-overlap Chae gates 0.825–0.869; deep occupation n_BE(0.1) =
9.508.

RULINGS (adopted in full, wording changes executed in the two stage
entries above + PREDICTIONS + LEDGER, same commit):
- **10J: UNPATCHED HOLE YES — scoped to the R-C consistency row**,
  not the letter, not any exact result, not credence. His J-1 (the
  round's best work): the "bound κ ≤ 1" is the D(1) ≥ 0 condition —
  a rung-count-arbitrary truncation whose three-rung sibling
  (κ ≤ 2/3) would EXCLUDE the locked world the row claimed to
  contain, and whose required deep extension is refuted by the
  stage's own fold row. **The "same world" sentence is RETRACTED.**
  J-2: the "compliance-divergence theorem" relabeled a SUSCEPTIBILITY
  RESTATEMENT of the two-rung degeneracy (not independent; plus his
  catch of the unrelated κ = 0.4 resonance inflating the low-κ
  ladder). J-3: "EP-anchored" scoped to the vertex; the measured-κ
  identification stays seam-conditional. Letter AFFIRM
  K-CONDITIONAL (P-lit now WEAKER — no upgrade path from this
  construction); credence HOLD 8 confirmed untouchable.
- **10I: UNPATCHED HOLE NO.** Letter AFFIRM I-FEASIBILITY-LIMITED
  with the meaning relabel (instrument-null, channel untested). His
  I-2: the "≤ ~2% of floor" claim was iso-only — the meter is
  ASYMMETRIC (group end reaches 110%; powered exactly where the
  science doesn't need it). His I-3: the census is a systematic
  UNDERCOUNT (UGC-designation aliases; 22/28 unmatched are UGC) —
  and he force-fixed the exemplar alias and showed the letter is
  ROBUST (the calibration gate fails harder, ρ → −0.126). His I-4:
  the Θ₁ offset −10.45 vs Karachentsev's −10.96 = benign
  mass-convention shift, non-propagating.
- SUCCESSORS ADOPTED AS NAMED: 10J — the SELF-CONSISTENT OCCUPATION
  BOUND (the P-lit replacement; the K-BOUNDED path); the κ = 1 level
  crossing as a DR4-era low-frequency-response reading; the meter
  κ-CONVENTION pin (bare vs renormalized ω) before any ridge-vs-lock
  adjudication. 10I — density-field e_N (2M++/Cosmicflows) × y ≳ 1
  kinematics as the only real void instrument; the free ALIAS-TABLE
  fix for every future LV cross-match; **the GROUP-END gate-closing
  test (p → ½ toward dense environments) — the axis this data
  actually supports.**

NEW TRAP FINGERPRINT #20 (DIARY): when a derived "bound" comes from
truncating a FAMILY of conditions, scan the family before quoting —
if the bound moves with the truncation and the truncation lands
exactly on the target value, the choice was outcome-determining, not
physical (the derivation-side sibling of the round-10
"pass-the-distribution" rule).

CREDENCE: no moves in any cell (anomaly-real 53 / mech-conditional 8;
both rounds' cells pre-signed HOLD; the referee confirmed neither
stage's movable cells could fire). Report archived
REVIEW-ROUND30-OPUS.md (never commit).

**Plain verdict: SUCCESS** — the round did exactly what rounds are
for: every number survived independent rebuild in both directions,
the one live overclaim (the "same world" sentence, the arc's most
quotable line) was caught before it entered any paper, and both
stages' successors came out sharper than the stages left them.

*ELI12: We asked a fresh robot referee to break both of today's
results. He redid every single number and got exactly ours — no
mistakes. But he caught our prettiest SENTENCE overreaching: we said
"the ladder math and the sky's thermometer point at the same
answer," and he showed our ladder rule changes its answer depending
on how many rungs you check — so the agreement was partly our choice
of where to stop counting. We deleted the sentence, wrote down the
rule that would make it honest (check the whole ladder, self-
adjusting), and kept everything that survived: the exact math, the
census, and two sharper next experiments.*

---

### Stage 10K (2026-08-09): THE GROUP-END POWER SPEC — the R30 dividend priced dead; the alias census banked

Pre-reg 008109e (bars locked blind: contest signal S_def 4/9, stratum
contrast 1/2 σ_p). Question: ROUND 30's inverted dividend said the
neighbor meter is powered at the GROUP end — is a gate-closing test
(p → ½ toward dense environments) powered for the LIKELIHOOD at
SPARC × LV grade? Executed R30 condition 7 on the way (the SIMBAD
alias table, counts only). **Letter K-POWER-DEAD, all gates green,
first run. The sky's A-vs-B preference was never read (firewall
held).**

- **The alias census** (counts only; the 10I letter stays closed):
  SIMBAD TAP resolved 25/28 of the unmatched-at-D≤12 names (465
  alias rows; D564-8, D631-7, KK98-251 honestly unresolved, logged in
  the CSV header; fetch calcs/fetch_sparc_aliases.py). Census under
  the fixed matcher: **matched 48 (was 31); iso 20 (was 13); iso
  transition-crossers STILL 3**; group Ti1>0 = 27 (4 crossers),
  Ti1≥1 = 13 (2). The alias fix finds seven more isolated dwarfs
  (UGC05721=NGC3274, UGC05764=DDO083, UGC05829=DDO084,
  UGC05918=DDO087, UGC08286=NGC5023, UGC08490=NGC5204,
  UGC12632=DDO217) — none of them transition-crossers, so the
  R30-flagged undercount resolves UPWARD in matches but the
  science-limiting count (iso crossers) stands at 3. The
  UNGC-vs-Chae Spearman under m1: ρ = 0.283 at N = 24 (recorded;
  still nothing like a calibration; the 10I letter is closed and
  does not reopen).
- **The group end is one galaxy, and it is data-orphaned:** the only
  neighbor-powered row (term ≥ 50% of the 0.01 floor) is NGC2976
  (109.9%) — and NGC2976 has NO Chae e_N row, so the lone carrier of
  the group-end signal cannot even enter the contest. The matched-set
  median neighbor term is 1.3% of the floor. "Powered at the group
  end" was, at SPARC × LV grade, a statement about one galaxy the
  gate world cannot reach.
- **The power spec** (Fisher with exact per-galaxy dv-profiling,
  calibrated on the archived 9V r-axis; reader identity |d| = 0.024
  at the archived grid point, r=0 anchor −12152.49 reproduced;
  calibration C_fish/C_arch = 2.69 → measured deflation D_f = 0.371,
  a reusable number): contest signal **S_def ≈ 0.000 at both r =
  0.3365 and r = 0.5** (bars 4/9 — dead by ≳3 orders); stratum
  contrast **Δp(Ti1>0) = −0.00059 = 0.008 σ_p** (bars 1/2). Both
  legs dead, and B was constructed as an UPPER BOUND (the neighbor
  sum stacked on Chae's 2M++ value, double-count direction
  disclosed) — a dead upper bound is decisive.
- Credence: pre-signed HOLD mech 8 / HOLD anomaly-real 53 (all
  cells). PREDICTIONS P2 annotated (the in-catalog group-end axis is
  priced; the ordering prediction stays fully out-of-sample).
  Successor: the gate-closing test needs genuinely dense
  environments — cluster-infall or group-catalog samples with
  decomposed kinematics, or DR4-era — not the Local Volume.

**Plain verdict: SUCCESS** (as a pricing instrument — the R30
dividend is now a number and the number is dead; the census upgrade
is banked and every future LV cross-match inherits the alias table).

*ELI12: We asked: can we hear the universe's hum get muffled where
galaxies crowd together? First we fixed the address book — 17
galaxies were listed under two names, and matching them properly
found 7 more loners (but no new usable ones). Then we priced the
muffling test: our map of crowded places has exactly ONE crowded
galaxy in it, that galaxy is missing the other measurement we'd
need, and even pretending everything lined up perfectly the
predicted signal is hundreds of times smaller than our noise. So we
wrote the price tag down and moved on — this test needs a crowd,
and our neighborhood doesn't have one.*

---

### Stage 10L (2026-08-09): THE SELF-CONSISTENT OCCUPATION BOUND — the K-BOUNDED path closes NEGATIVE; the occupation is bare-pinned; the meter convention pinned

Pre-reg 008109e + amendment A1 b9cfe04 (post-run-1 self-catch logged
pre-quote, run 1 preserved as data/stage10l_occbound_run1.txt: the
G10L-2 κ→0 regression bar was mis-set ARITHMETIC — at κ = 10⁻⁹ the
first-order response of the fixed point is ~2×10⁻⁷, five orders
above the 10⁻¹² bar, so the gate tested the finite-κ response, not
the limit; replaced by a linear-extrapolation limit test with a sign
check; no physics row changed between runs). **Letter L-PINNED,
gates 5/5, both verdict clauses firing family-stable** (trap-#20
discipline: four spacing offsets + two truncated-Gibbs conventions +
two independent existence detectors — the verdict is not allowed to
inherit a construction choice).

- **THE R30 QUESTION ANSWERED.** The round-30 repair of the two-rung
  bound, verbatim: solve the self-consistent occupation (n̄ set by
  the dressed spacing at the physical T_dS) to fixed point and ask
  "whether κ → 1 is where self-consistency first fails." Answer:
  **NO.** The fed-back occupation n̄ = n_BE(x·s(n̄)), s = 1 −
  κ(n̄+c)/2, has no fixed point at sky occupations: κ_max(x) =
  0.015–0.14 over the deep window x ∈ [0.03, 0.3] — **6.6–62×
  below the locked κ = 0.925** — and at κ_lock the first fixed
  point appears only at x ≈ 3.6, the Newtonian arm, where ν ≈ 1 and
  nothing needs explaining. Failure sits at κ ~ O(x), far below 1;
  κ = 1 is not the failure edge, so no bound lives there.
- **PRODUCTION:** no member of the family carries the deep RAR. The
  SCHA members have no fixed point at all at (κ_lock, x = 0.05);
  the truncated-Gibbs members (Gibbs on the rising branch of the
  literal ladder, both kept-set conventions) saturate at ν ≤ 1.93
  against the measured demand 19.0 — under-production ≥ 9.9×/13×.
- **THE CONSEQUENCE (the honest end-state of the R-C axis):** the
  dictionary's occupation is **BARE-PINNED** — structurally (this
  stage: the softening feedback, the vertex's own polaron direction,
  cannot exist at sky occupations) and by measurement (the 5P
  stiffening family solves fine at every β — negative feedback —
  and the sky sits at its no-feedback point: deep/binary β → 0, 5R
  bound β < 0.03; boot = the β = 1 member, sky-dead 5M; the β = 1
  member regresses to the 4F form ν = 1 + n_BE(νy) at 5×10⁻¹⁵).
  P-lit has no self-consistent completion in the scanned family ⇒
  **no occupation-consistency bound on κ exists. The K-BOUNDED path
  is CLOSED, negative.** (Scope, named: mean-occupation feedback
  constructions; joint coherent-state/Hartree treatments out of
  scope.) This is also the fold-row resolution: the "self-consistent
  structure" the fold demanded is NOT a re-thermalized occupation —
  consistent with the 10C separation theorem (the bath occupation
  cancels between time-orderings; the drive is external). The
  occupation is bare-pinned twice over: measured and structural.
- **THE CONVENTION PIN (R30 successor 3, exact):** bare form ν =
  κ_b/x + (1−κ_b/2) + κ_b x/12 ⇒ the deep-amplitude meter reads
  a₀_fit/a₀_h = κ_b² and the c₁ meter reads κ_b = 2(1−c₁); the
  first-transition (renormalized-ω) convention maps κ_r =
  κ_b(1−κ_b/2) with maximum EXACTLY ½ at κ_b = 1 — nearly κ-BLIND
  over the whole measured band (0.494–0.500 for κ_b ∈ [0.888,
  1.10]). All archived meters (4S/4Z/10D) parameterize x =
  √(g_N/a₀) — bare by construction, which IS the vertex convention
  κ = 4g²/(Ωω_bare). ⇒ **the two-pole tension (lock-meter κ_b =
  0.925–1.00 vs c₁-meter 1.10–1.48) is CONVENTION-ROBUST** — a
  physical tension (κ runs; the one-κ form is 10D-rejected), not
  bookkeeping. Never print "κ measured ~1" (R24 discipline).
  Curiosity recorded: at closure a first-transition meter would
  read κ_r = ½ — the same number as c₁; both are faces of (1−κ/2).
- The κ = 1 level-crossing low-frequency response (R30 successor 2)
  booked as a READING only — no numeric kill condition derived, so
  no PREDICTIONS row (discipline held).
- Credence: pre-signed L-PINNED = **HOLD mech-conditional 8**;
  anomaly-real 53 untouched (no sky fits; every quoted sky number is
  an archived meter).

**Plain verdict: SUCCESS** — the round-30 successor executed to its
honest end. The bound question is now CLOSED (negative) rather than
open; the derivation program knows one more thing it is NOT allowed
to assume (re-thermalization on the dressed ladder), one more thing
it must explain some other way (why the occupation is bare-pinned —
the 10C separation theorem is the standing microscopic account), and
the meter bookkeeping the ridge-vs-lock discussion needed is pinned.

*ELI12: Our formula says gravity's "warmth" fills energy levels like
a thermometer. The referee asked: if filling the levels also
squeezes them together (our math says it does, slightly), maybe
demanding perfect self-consistency forces the magic knob to be
exactly 1? We checked every version of that idea — and
self-consistency doesn't fail AT 1, it fails almost immediately, for
any knob setting big enough to matter. That means the levels cannot
be re-filling themselves on their own squeezed spacing at all: the
filling follows the ORIGINAL spacing — which is exactly what the
data had separately told us, and what our bath-cancellation theorem
said microscopically. No bound from this road, and now we know the
road is closed rather than unexplored. Bonus: we proved two of our
rulers genuinely disagree about the knob (it's physics, not
mislabeled units), and found a cute fact — measured with the
"squeezed" ruler, the knob at its special setting would read
exactly ½.*

## Honest credences (2026-07-21, end of Stage 2)

Wide-binary velocity excess is real physics (not systematics): ~65% (up from agnostic; our
pipeline's tails+RUWE+noise accounting did the moving; Banik conflict keeps it under 80%).
Low-acceleration anomaly pattern (galaxies+binaries) is one law: ~55%. Horizon/entropic
microphysics specifically: ~15%. Superfluid-DM hybrid: ~20%. The program remains underpriced.

## Historical credences (2026-07-21, pre-Stage-1)

Anomaly pattern is real physics needing explanation: ~90%. Specifically horizon-inertia:
~10–15%. Condensate/hybrid DM with emergent a₀: ~25%. Program underpriced relative to its
falsifiability density: ~certain.
