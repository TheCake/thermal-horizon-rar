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

## Honest credences (2026-07-21, end of Stage 2)

Wide-binary velocity excess is real physics (not systematics): ~65% (up from agnostic; our
pipeline's tails+RUWE+noise accounting did the moving; Banik conflict keeps it under 80%).
Low-acceleration anomaly pattern (galaxies+binaries) is one law: ~55%. Horizon/entropic
microphysics specifically: ~15%. Superfluid-DM hybrid: ~20%. The program remains underpriced.

## Historical credences (2026-07-21, pre-Stage-1)

Anomaly pattern is real physics needing explanation: ~90%. Specifically horizon-inertia:
~10–15%. Condensate/hybrid DM with emergent a₀: ~25%. Program underpriced relative to its
falsifiability density: ~certain.
