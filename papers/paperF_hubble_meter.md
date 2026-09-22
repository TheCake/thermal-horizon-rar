# A Hubble meter from galaxy rotation curves: first execution of the a0-to-H0 inversion

**Draft 0.1 — SKELETON (not for circulation; prose incomplete; every
[TODO] blocks a section until discharged. Style contract: papers/STYLE.md
is binding. Numbers below are the banked operative values; provenance in
the reproducibility appendix at assembly time.)**

Author: Filip Hajek (independent researcher).
Acknowledgments will carry the repository collaboration statement
(Claude, Anthropic), as in Papers 1-3.

---

## Abstract (seed; rewrite at assembly)

The radial acceleration relation of galaxies is described by an
acceleration scale a0 whose measured value coincides with cH0/(2 pi).
If that identity is physical, rotation curves are a Hubble meter: a0,
measured from galaxy dynamics, returns H0 with no supernova rung and no
period-luminosity rung. We report the first execution of that
inversion. Because 55 percent of SPARC distances are Hubble-flow
distances that assume H0 = 73, the naive inversion is circular; the
instrument is therefore built on the anchored subsample alone (78
galaxies with tip-of-the-red-giant-branch, Cepheid, supernova, surface
brightness fluctuation, or maser distances after a CosmicFlows-4
per-method reclassification), with the flow subsample retained only as
a self-consistency leg. The anchored leg reads H0 = 65.4 +/- 5.0
(statistical) +/- 0.5 (convergence), with membership variants spanning
64.8-65.4. The dominant systematic is not statistical: the choice of
the relation's functional form moves the reading across 65.4-70.4, and
we measure the levers of that band directly (a per-galaxy
form-preference map; a single galaxy, UGC03580, carries the flow
subsample's contrary form lean; an inner-angular-resolution cut
resolves nothing and excludes nothing while costing the anchored
acceleration scale 0.8 km/s/Mpc). The meter
shares its calibration pegs with the ladder but with roughly twice the
gearing, so at a precision of 2-3 km/s/Mpc it separates a
pegs-systematic world from a new-physics world; at the current
precision it separates nothing and we say so. The method is
timestamped for the ~4000-galaxy BIG-SPARC era, where the statistical
term shrinks toward the deciding range. [TODO: word count; 3-5
headline numbers only]

---

## 1. Introduction

- The identity: the measured galaxy acceleration scale sits on
  cH0/(2 pi) (the temperature lock; P2 sections 3 and 5). Cadoni &
  Tuveri derived the relation's form from horizon thermodynamics; the
  identity makes a0 a cosmological dial regardless of mechanism.
- The proposition: run the identity BACKWARD. Nobody has (three
  independent literature sweeps found no published a0-to-H0 inversion
  as a measurement; the forward direction exists: van Putten).
- The nearest neighbor and the differentiation [SML20 primary read
  done 2026-09-22, verify quotes at print]: Schombert, McGaugh &
  Lelli (2020) measured H0 = 75.1 +/- 2.3 (stat) +/- 1.5 (sys) from
  the baryonic Tully-Fisher relation, calibrated on 50 Cepheid/TRGB
  galaxies (30 of them SPARC) and applied OUT to 95 SPARC flow
  galaxies through CosmicFlows-3 velocities. Their dominant
  systematic is the flow-model choice (72.8-77.5 across four models);
  no acceleration scale, no radial acceleration relation, and no lock
  appears anywhere in their method (confirmed at primary grade). The
  differentiation sentence: same catalog, opposite direction — they
  ladder out through a flow model, we invert the acceleration scale
  on the anchored subsample with no flow model; the two instruments
  share calibration pegs and fail in different places.
- What this paper is: the instrument's design (the circularity and its
  two-leg resolution), its first light, and the measured levers of its
  systematics. What this paper is not: an H0 arbiter at current
  precision.

## 2. The circularity and the two-leg design

- The census (from the SPARC distance-method table, re-parsed twice):
  97/175 galaxies carry flow distances assuming H0 = 73; 50 carry
  H0-independent anchors; 28 share the Ursa Major cluster distance.
  A naive whole-catalog inversion inherits the assumed 73 and is not a
  measurement.
- Leg A: the anchored subsample, direct inversion. Leg B: the flow
  subsample, where the fitted a0 depends on the assumed H0, solved
  self-consistently for a0(H0) = cH0/(2 pi). Agreement between legs is
  the design's internal test; the identity has a reason to pass it,
  numerology does not.
- Pre-registration discipline: legs, bars, and letters were frozen
  before any sky fit, across five stages (10S design + first fire; 10T
  anchored regrow; 10U form contest; 10V stratification; 10W map and
  resolution cut), each with an adversarial review round adopted in
  full; the record is public.

## 3. Data and sample

- SPARC 175 (quality-cut 153 with inclination >= 30 and Q <= 2); the
  rotation-curve mass models with fixed stellar mass-to-light
  conventions (disk 0.5, bulge 0.7 at 3.6 micron) and per-galaxy
  velocity-error propagation.
- The CosmicFlows-4 per-method reclassification (leg A 67 -> 78): nine
  TRGB plus SBF plus SNII per-method upgrades with distance errors
  falling from 15-30 percent to 4-9 percent; the Ursa Major block kept
  with its two direct-Cepheid members as a frozen co-read; the LITTLE
  THINGS candidates excluded by a pre-registered baryon-model
  compatibility clause (their velocities agree to 0.05 dex, their
  published baryon models sit one-sidedly 0.13 dex low against SPARC
  overlaps).
- Census-pending objects named ([TODO: D564-8, D631-7 status at print
  time]).

## 4. Method

- The hierarchical likelihood: per-galaxy distance and inclination
  nuisances (the SPARC error model), a per-galaxy mass-to-light
  scatter prior, a global disk mass-to-light ratio, an intrinsic
  scatter term, and a shared Ursa Major distance nuisance; the
  operative estimator is the deep-refit form with engine-matched
  uncertainties (a flat-treatment variant is quoted only as an
  envelope, never as the statistic).
- The functional-form family: the Bose-Einstein form (the derived
  identity), and three measured alternatives spanning the plausible
  transition-sharpness range (p = 0.65 tail, the geometric-mean
  construction, the self-consistent bootstrap form). H0 is read at
  each; the family span IS the dominant systematic and is reported as
  a band, not folded into the error bar.
- Reproducibility: every fit reproduces bit-identically under
  parallel execution (a certified deterministic worker pool); all
  pre-registrations, gates, and reviews are in the public repository.

## 5. First light (the anchored leg)

- H0(leg A) = 65.4 +/- 5.0 (stat) +/- 0.5 (convergence plateau across
  three independent optimizer implementations).
- Membership variants: 64.8-65.4 (anchored-only, no-UMa, UMa-only
  63.5). Composition slices stable within 1 sigma.
- The flow leg is power-limited by its own structure: the
  self-consistent solve sits 0.14 from the estimator's pole (a x7.1
  error amplifier), and the bootstrap crosses the pole in roughly a
  quarter of resamples, so no finite variance exists at SPARC depth;
  the leg-B agreement test is deferred to gas-dominated flow samples
  (WALLABY-class), where the mass-to-light pole is absent by
  construction. [Never quote a leg-B number from the archive; the
  retracted strings are listed in the verdict files.]
- Validity domain: 62-75 (an injection arm measured a real -3.0 +/-
  0.7 percent one-sided edge bias at truth 85; inside the domain the
  meter is unbiased at the stated grade).

## 6. Systematics: the measured levers

- (a) THE FORM BAND, 65.4-70.4: the four-member contest on the
  anchored leg demotes the bootstrap form at report grade (~2 sigma;
  paired-bootstrap and permutation constructions agree) and cannot
  separate the band's edges at N = 78 under any defensible
  calibration; the geometric-mean member survives on a bootstrap
  co-requirement the data refuse at 1.26 sigma. The band is quoted in
  full.
- (b) THE STRATUM QUESTION, DISSOLVED: the apparent cross-subsample
  form inversion (the flow subsample preferring the bootstrap form)
  measures 0.71 sigma against its own realized galaxy-block scatter
  and is carried by one galaxy: UGC03580 contributes -47.8 of the
  flow subsample's -48.7 total (52 percent of the block variance);
  with it removed the remaining 70 galaxies prefer the Bose-Einstein
  form. Composition (gas fraction) carries no information beyond
  acceleration coverage (partial rank correlation p = 0.54 after
  conditioning); coverage itself reaches only point grade.
- (c) THE RESOLUTION LEVER (adopted Round-49 wording; the operative
  record is the stage verdict file): the inner-angular-resolution
  lever moves the anchored subsample's form preference by at most
  0.10 of its own scatter and the flow subsample's by at most 0.30 of
  its own, in a direction that reverses between 20 and 30 arcsec; it
  resolves nothing and excludes nothing, and it costs the anchored
  acceleration scale 0.8 km/s/Mpc in H0. The 30-arcsec reversal is a
  point-cut effect, not a population change (holding the population
  fixed, the seven ejected low-surface-brightness dwarfs carry -0.6
  of the -24.5 move). A paired galaxy bootstrap gives P = 0.62 (bar
  0.95) that the cut moves the flow lean toward the anchored
  subsample's preferred form: inner resolution is neither established
  nor excluded as the carrier of a residue that is itself 0.71 sigma.
  The carrier structure is unstable under the cut while the total is
  not (UGC03580's variance share runs 52 -> 46 -> 14 -> 0 percent
  across the cut grid) -- the signature of noise, not of a carrier.
  [Never quote cut-flow a0 values: their implied H0 (77-84) sits
  outside the meter's validity domain.]
- (d) PEG SHARING AND GEARING: the meter hangs from the same
  calibration pegs as the ladder (MW parallaxes / LMC / N4258 feeding
  TRGB and Cepheid zero points) with measured gearing 1.3-1.6 times
  the ladder's response to a coherent peg rescale, same direction.
  The three-world separation table [Figure], scored at current
  precision: L 0.3 / N 1.5 / P 1.9-2.7 sigma (statistical) — leaning
  against the shared-pegs-short world, deciding nothing until the
  form band is pinned.
- (e) Conventions disclosed: the frozen bulge mass-to-light (0.7,
  never varied by any gate; freeing it moves the flow contest by -22
  percent), the prior-free global disk mass-to-light, the shared UMa
  distance treatment.

## 7. Forecast and the BIG-SPARC timestamp

- The statistical term at N = 78 is 5.0; the honest requirement for
  form-band arbitration is N ~ 100-1000 (measured from the contest's
  own separation calibration; smaller forecasts from an earlier
  injection model were retracted on review and are not quoted).
- TRGB-era regrow of leg A (CF4 census in hand) buys roughly 9-10
  from 12.5 at the previous census; BIG-SPARC (~4000 homogeneous
  galaxies, in preparation) is the refire target, and its distance
  choices will embed an H0 — this paper's provenance-split design is
  the instrument for it, timestamped first.
- World V (the differential-expansion / local-void class): under the
  lock it transduces into a radial a0 gradient plus a dipole — a
  registrable discriminant against the uniform new-physics world; the
  hemisphere instrument (SPARC north x WALLABY south) is future work.
  [Cite Haslbauer/Kroupa/Banik with care; the KBC void becomes
  lock-testable in rotation curves for the first time.]

## 8. The redshift axis (synthesis)

- If a0 = cH(z)/(2 pi), a0 RISES with redshift. Both measured epochs
  lean that way (MIGHTEE-HI z ~ 0.09; MUSE-DARK III z ~ 1 lock
  comparisons: the lock intercept sits inside all three of their a0(0)
  rows; slope comparisons are release-limited). Four rival scalings
  (van Putten, Gillot, Escala variants, frozen-a0) are sign-opposite
  or flat: the z-axis is the cleanest external discriminant and it is
  live now. [Cross-cite P2 section 5.1; keep this section short — it
  is a synthesis, not a new measurement.]

## 9. Discussion

- What the meter is today: the first executed inversion, reading
  65.4 +/- 5.0 +/- 0.5 with a 65.4-70.4 form band, Planck-adjacent,
  separating none of the live Hubble-tension worlds at this precision.
- What it becomes at sigma ~ 2-3: a referee with different systematics
  from every existing route (no SN rung, no period-luminosity rung,
  gearing x2 on the shared pegs).
- Credence paragraph (house honesty rule): the program's credence
  that the underlying anomaly is real physics stands at ~53 percent;
  this paper's H0 claim is conditional on the lock and is presented
  as an instrument, not as evidence for the lock.
- Honest nulls: no age-of-universe consequence; no relativistic
  completion supplied here; solar-system screening is exponential in
  the measured class.

## 10. Conclusions

[TODO after Round 49 adoption + figure gates.]

## Appendices

- A. Transparency and corrections: the five review rounds (45-49)
  with the adopted downgrades named (the as-fired agreement letter
  downgraded to gray; the sigma engine-mismatch corrected; the
  injection-calibration hole; the stratification relabel; [TODO R49]);
  the retracted strings listed by name so they cannot be re-quoted.
- B. Reproducibility table: stage -> script -> output -> commit, one
  row per quoted number.

## Figures (each gated by calcs/paperf_figures.py before assembly)

1. The census and the two-leg design (distance-method breakdown; the
   circularity arrow; the leg architecture).
2. Leg-A first light: the H0 posterior with membership variants and
   the form band overlaid.
3. The per-galaxy form-preference map (data/stage10w_map.csv): Delta
   per galaxy vs acceleration coverage, colored by leg, carriers
   labeled — the systematics figure.
4. The three-world separation forecast: meter reading vs sigma with
   the L/N/P/V worlds and the gearing band.
5. The redshift axis: a0(z) at the two measured epochs against the
   lock line and the four rival scalings.

## Must-cite list (verify-before-print; differentiation one-liners)

- Cadoni & Tuveri 2019 (the derivation; the identity's origin) + their
  anisotropic-fluid cosmology paper (2002.06988).
- Schombert, McGaugh & Lelli 2020, AJ 160, 71 (the nearest neighbor;
  ladder OUT vs inversion IN; same catalog, no lock) [primary read in
  progress].
- Lelli, McGaugh & Schombert 2016 (SPARC); Tully et al. CosmicFlows-4;
  Anand et al. EDD/TRGB.
- van Putten (forward direction only); Haslbauer, Kroupa & Banik
  2008.07524 (world V; care); Migkas et al., Secrest et al. (the
  anisotropy watch-list); the high-z a0(z) constraint literature
  [verify scout IDs before citing]; MIGHTEE-HI (2608.03576); Ciocan
  et al. (MUSE-DARK III).
- Schombert-McGaugh-Lelli M/L conventions; the DR4 protocol paper
  [Proposal A] if it exists by submission time.

## Guardrails (bound at assembly; from the stage verdict files)

- NEVER QUOTE: the retracted leg-B numbers; the flat-engine sigma
  strings; the N ~ 25 forecast; the flow gas-rich stratum a0 and its
  implied H0; the UMa-26 screening a0; windowed-leg a0s as meter
  readings; injection-bar significances; bootstrap probabilities
  dressed as sigmas.
- The band 65.4-70.4 and stat 5.0 are the operative pair; no
  headline may compress them into one number.
