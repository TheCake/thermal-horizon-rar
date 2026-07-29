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

## Honest credences (2026-07-21, end of Stage 2)

Wide-binary velocity excess is real physics (not systematics): ~65% (up from agnostic; our
pipeline's tails+RUWE+noise accounting did the moving; Banik conflict keeps it under 80%).
Low-acceleration anomaly pattern (galaxies+binaries) is one law: ~55%. Horizon/entropic
microphysics specifically: ~15%. Superfluid-DM hybrid: ~20%. The program remains underpriced.

## Historical credences (2026-07-21, pre-Stage-1)

Anomaly pattern is real physics needing explanation: ~90%. Specifically horizon-inertia:
~10–15%. Condensate/hybrid DM with emergent a₀: ~25%. Program underpriced relative to its
falsifiability density: ~certain.
