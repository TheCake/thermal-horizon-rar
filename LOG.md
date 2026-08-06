# LOG

One line per stage: what we did -> what came out. Plain language, newest at
the bottom. Full detail lives in NOTES-horizon-inertia.md; every headline
number with provenance is in LEDGER.csv; signed predictions in
PREDICTIONS.md. PAPER.md is frozen — it gets rewritten only when Filip
calls it, not per stage.

## Galaxies + the function (2026-07)

- 1: Fit the RAR screening index on SPARC -> p = 0.44 fixed-M/L; M/L-marginalized (4H) p = 0.58 +/- 0.12; p = 1/2 sits comfortably inside.
- BE identity: found RAR nu = 1 + n_BE(x) exactly -> already published (Cadoni & Tuveri 2019); priority retracted, our role = first empirical tests.
- EFE solver: built + validated the QUMOND external-field solver -> BE-EFE ~4% weaker than simple-nu; boost tables feed the binary pipeline.
- 4A/4B: First test of the C&T structure -> the 1/2-branch beats standard-mu 198-200/200 raw; deflates to a strong sign-robust lean under honest likelihood.
- 4E: Lensing RAR leg -> honest null (resolving power ~0.1 sigma); Newton still loses by thousands; SPARC declared agnostic between BE and simple (earlier lean retracted).
- 4F: Bath matrix -> simple-nu IS the exact classical self-consistent thermal bath; BE-vs-simple = quantum-vs-classical bath; the 1/4 "boot" cell disfavored.
- 4S: Measured c1 on SPARC+lensing -> 0.45 (0.39-0.52); 1/2 interior; c1 = 0 excluded ~7 sigma.
- 4T: RAR scatter is acceleration-dependent -> oscillator+floor reads N ~ 21 modes over a floor; constraint, not detection.
- 4U/4V: M/L second-moment control -> thermal trend survives; scorecard: galaxy a0 sits on the horizon value, binary a0 +2.5 sigma = sharpest internal tension.
- 4W: Full hierarchy -> identifiability boundary (channels degenerate on SPARC); the x~1 scatter bump survives every marginalization.
- 4X: Binary c1 = 0.37-0.50 -> two disconnected systems read the same dial (later erased at the landed posterior; see 7J-d).
- 4Y/4Z/5A: Hierarchical c1 = 0.26-0.45, 0 excluded everywhere; anchored subsample; the bump is point-level, not per-galaxy.
- 5B-5H: Chae environmental control passed (thermal survives, EFE-collinearity excluded); boot flips ahead under hier M/L, then the binaries veto boot; decomposition: hier galaxies want sharper screening (p ~ 0.65).
- 5I-5M: Quadrupole is amplitude-locked (no interpolating-function escape); geometric-mean bath built (zero params) leads galaxies; binaries prefer p = 1/2; the vertical control returns every a0 to the horizon and kills boot.
- 5N-5P: The mixing family -> one parameter beta spans BE(0)/gm(1/2)/boot(1); galaxies read beta ~ 1/2, binaries beta ~ 0: the two-system split in one dial.
- 5Q-5T: beta = 0 is the family's unique Bernoulli zero; binaries bound beta < 0.03; the galaxy beta is a compromise (deep votes 0, tail votes 1/2-3/4); Saturn's veto is beta-blind.
- 5U-5X: Spontaneous-fraction bath (zero params, running beta) -> leads all galaxy treatments; binaries improved-not-accepted; Saturn lock holds.
- 5Y-6D: Two-leg variants -> biggest galaxy leads yet, but binaries reject every sharpened function under the dominant ambient field; the split is SYSTEM-level; pointwise drive-weighting excluded.
- 6E-6G: THE AMBIENT-GATED BATH -> beta = (ambient classicality)^2 x local quantumness, zero params, postdicts both measured tails; first single function to pass both systems.
- 6H-6J: Gate grammar derived (JC pull, leg count); measured ambients improve the fits; bootstrap 37/40 = top grade, binary-compatible where F4 is not.
- 6K/6L: Lab-analog exclusion round (rate-based realizations die; pre-committed credence strike taken); leg-count claim deflated to a lean (correction #13).
- 6M/6N: Structured-bath + full Lindblad -> the lab leg CLOSED (flat bath = exact source-locking; a buildable fridge measures zero); bath-microphysics conditional cut to ~8%.
- 6O/6P: Frozen-bath tests -> no support either side; real dividend: an s-flat per-system scatter worth ~+37 lnL exists in v7 (the width object's first sighting).
- 6Q: Built LEDGER.csv + the world table with six audit gates; standing rule: every headline number ships its ledger row in the same commit.
- 6R-6T: Resolution bath -> largest plain galaxy lead on record, but binaries shape-reject it; triangulation complete: the ambient gate is REQUIRED.
- 6U/6V: The gate derived from KMS detailed balance (s^L = Boltzmann cost of L borrowed quanta); untied contest tolerates the tied form.
- 6X: First constructive mechanism result -> the lending probability IS the KMS ratio (toy grade).
- 6W: My scalar-EFE escape hypothesis -> excluded by the data in ninety minutes; dividend: first measurement of the EFE composition character; Saturn sharpened, not solved.
- 6Y/6Z: Reservoir identification -> one collective ambient mode selected (M=1 theorem); the in-sample ordering claim withdrawn (correction #14, shuffle control).
- 7A: PREDICTIONS.md signed ledger opened; Einstein fluctuation test -> shot bath excluded; quantum-vs-classical blocked by the x~1 bump.
- 7B/7C: Bump hunt -> it is inner-disk astrophysics (R < 1.5 R_d), not the law; gamma is not measurable on SPARC at any grade tried; honest close.
- 7D/7E: The vacuum +1 is load-bearing (classical share dies by +556); the lending gate translated into a buildable cQED bench experiment.
- 7F: Mixing-mean uniqueness -> only the geometric mean reproduces the measured two-system pair (p = 0 interior).
- 7G/7H: Saturn in the trajectory formulation -> 451 orders below Cassini; the binary tie RESTORED (the 7G gap was proxy discretization); the quadrupole tension is a formulation property, not a function property.

## Binaries: build-up and the anomaly (2026-07)

- 2: Built the EDR3 wide-binary pipeline (14,071 pairs) -> median velocity boost 1.086 (CI 1.064-1.110); RUWE-stable; triples capped ~5% by tail shape.
- 3A: Found the realization systematic -> orbit-population draws shift lnL more than the Newton-vs-MOND gap; apparently unpublished.
- 3B-3E: Hierarchical v3 -> alpha corner-seeking traced to missing width; a multiplicative smear ~0.2 fixes it; the Newton rejection honestly deflates +264 -> +63.
- 3F-3I: v4 -> alpha interior both laws; budget alpha = 0.98+/-0.23 (simple) / 1.21+/-0.30 (BE); Newton loses all 1000 bootstraps.
- 3J: Measured the photometric mass error -> 12x too small to be the smear (hypothesis refuted); found 12.3% overluminous components.
- 3K: v5 physical-companion model -> companions exonerated as the boost; the fit prefers the smear by 420; f_comp fenced at 0.1 (revisited in 7J).
- 3L + 4G: Velocity angles are U-shaped -> Hwang+22 already published the superthermal e-distribution; our radial excess = confirmation, not discovery.
- 3M/3N: Joint 2D (velocity x angle) fits -> circular reading vetoed; w_rad = 0.20 near-parabolic sub-population identified as the broadening carrier.
- 3O-3Q: v7 closed the catalog-selection systematic; budgets alpha ~ 1.5 (superseded by the convention fix below).
- 3S/3T: Caught our own external-field convention bug -> at the correct field BOTH laws land at the parameter-free alpha = 1.
- 3U/3V: Final measurement -> alpha = 1.18+/-0.11 (simple) / 1.13+/-0.13 (BE); Newton excluded in all 2000 bootstraps; binaries lean simple.
- 4I: Chance-alignment robustness (20x tightening) -> boost stable.
- 4J: Resolved the gamma = 82 islands; DISCOVERED the perpendicular-velocity ceiling: 11 pairs in the Newton-forbidden band, cliff at the boosted escape edge (leakage null 3.8e-9).
- 4K: Solar quadrupole -> both families exceed Cassini ~4.3x; binary-calibrated reproduction of the Desmond-Hees-Famaey tension; escapes = modified inertia or EFE screening.
- 4L: MI-vs-MG on our data -> EFE-respecting MI ties MG; no-EFE MI dead; the data demand the EFE amplitude.
- 4M: Flyby template -> the ceiling residual closes as model-shape; boost untouched.
- 4N: Reconciliation with Banik-style pipelines -> freeing companions absorbs ~60% of the Newton deficit but the detection never flips; the gamma channel is the parameter-protector.
- 4O: PAPER v1.1 verified by a nine-agent primary-source pass; ceiling founding prior art found (correction #11).
- 4P/4Q: The Cookson read caught our missing perspective correction -> present, 5x too small to matter; corrected anchor 1.078 (CI 1.052-1.103).
- 4R: Corrected-velocity budget -> alpha shift zero; numbers stand.

## The review-and-repair arc (2026-07-25..28)

- 7I: External review adopted -> the freeze (function search closed at 13); ablation ladder; the strict-multiplicity cut fired MATERIAL: the tension opened (credence 70 -> ~60-65).
- 7J: The decisive instrument (completeness + marginalized posterior) fired against us -> companion-win at the measured prior; then corrections #15-18 dismantled our own prior scale (common-mode over-attributed as companions); the fenced +99-110 Newton rejection is DEAD; verdict AMBIGUOUS at the literature anchor (credence to ~45).
- 7J-z/z2b/z2c: The 2D photometric mixture measured properly -> host rate ~0.3/component, twin-heavy q-shape (the x3 literature tension dissolves); the width channel enters the model: sq = 0.2 demanded interior (P = 1.00).
- 7J-g: The gamma decider -> velocity-only pipelines report a PHANTOM alpha ~ 0.5 that the angle data veto; the field-reconciliation centerpiece.
- 7J-z3-z8 (review rounds 9-12): anchor-strength curve (sigma* = 0.02); the forced-multiplicity rejection is wobble-binding; twin-force and containment: the (band 9, cliff 2) census pair holds both collapse flanks.
- 7J-z5: The arm suite validates the machinery (null arm manufactures nothing; injections recovered); the sky's noise hunger flagged as real width-shape incompleteness (credence ~45 -> ~50).
- 7J-z6: Width-shape contest UNRESOLVED; product = the noise-ceiling profile (operative alpha 0.68-0.74 is noise-DILUTED; physical envelope reads 0.80).
- 7K-a/b: Forward accounting -> the landed cell produces about half the median excess (GRAY); the census count retires to the self-defending (9,2) pair.
- 7J-d: The unsuspension -> NO function discrimination at the landed anchor (all 16 laws within +/-8 lnL); the binary c1 claim erased; Newton function-robustly dead (dN >= +8 on all 32 rows).
- 7L: Cookson quantified -> their flatness replicates on our catalog; their cleaning discards 92% of the joint sample including the whole deep anchor; their null and our detection are compatible.

## The rivals-and-census arc (2026-07-29..08-05)

- 8A: Rivals' ladder -> catalog kills + the pincer: no standard-family function survives both the measured window and its own published Cassini verdict.
- 8B: LambdaCDM attractor -> does not print the zero-point digit (off-family).
- 8C: The p <= 3/4 ceiling test -> ambiguous; break side clean, per-galaxy clause power-limited.
- 8D: The removed-92% settling instrument -> BOOST-CARRIED 4/4; the Cookson fork resolves toward data, not contamination.
- 8F/-b/-c/-d: Fat-tail manufacture excluded 12/12; the census band is unproducible by the error-tail class; the bias-curve turn-on located (10-20% severity).
- 8G: Freed the eccentricity sector -> the width object SURVIVES (against expectation); the noise chase partially releases.
- 8H: Census shape-and-attribution -> all single-dial repairs fail; the census pair attributes to the companion sector (0.81/0.78).
- 8I/8J: Wobble survival-cap and saturation repairs -> both refused by the sky.
- 8K/8K-b: The collapse-world companion population is absent at object level; NSS crossmatch NINE-CLEAN -> the quiet-faker account loses its last in-catalog hiding place (credence 55 -> 57 by pre-signed map).
- 8L/8L-b: Response model licensed; the derived leakage kernel replaces the legacy curve -> KERNEL-PAID, the calibration debt retires; new MAP cells shipped.
- 8M: Joint-coherence scan -> census-coherent worlds exist only in prior-forbidden no-companion corners at ~5,000 lnL price (CLASS-CONTAINS).
- 8N: Flood anatomy -> two culprits: band flank = the sq tail, cliff flank = long-period companions.
- 8O: T2 external-RV archive pass -> RV-CLEAN at thin coverage (credence 57 -> 58 by pre-signed map).
- 8P/8Q: Both census flank repairs EXCLUDED (band flood is bulk-variance, not tail-shape; the measured subsystem period law is kinematically refused); caught a 1-ulp reader bug with the program's first lnL-grade identity gate -> new standing rule.
- 8R/8R-b: The pincer closes -> the width demand is uniform in gamma; NO multiplicative smear (any shape, any localization) reconciles kinematics and the census.
- 8S/8S-b: Gas-dominated c1 (referee T5) -> NOT the hoped M/L-immune confirmation: a real dial-tension (GD galaxies interior at c1 ~ -0.8); gas budget cleared.
- 8T: The additive noise floor is real (45 m/s, interior) and IDENTIFIES the long-standing fpm chase; it does not substitute the smear and worsens the census.
- 8S-c: The measured distance/inclination channel barely moves GD -> tension survives.
- 8U: The mixture instrument -> no zero-inflated smear reconciles kinematics and census either (locus demands monotone-opposed); minority-bars grammar miss owned; DATA-SIDE is the last reconciliation class.
- 8V: Galaxy-level bootstrap of the GD tension -> 0/300 replicates reach the dial; correlation suspect cleared; the GD-DD split is hardened (suspects left: selection, real deep-regime shape).
- 8W: Quality-strata tracking of the width object -> it has an ADDRESS: concentrates in high-RUWE / low-parallax-S/N pairs (clean half nearly zeroes it); physics-blind axis flat; readings = astrometric noise vs astrometry-channel multiplicity; DR4 arbiter.
- 8X: Regime-vs-composition on the GD tension -> GD is 100% deep (total confounding), but DD's own deep points NEVER vote like GD's (0/200): the split is galaxy TYPE, not regime; new named suspect = pressure-support corrections in slow-spinning gas dwarfs (8Y candidate: V_flat split).
- 8Y: Pressure-support dose test -> the GD dial is strongly V-ordered (slow rotators maximally negative, direction P=0.033) BUT the actual correction lever is out of reach by ~10x (k-curve; even sigma_eff=17 km/s moves lam only -1.31 -> -0.93): pressure support cannot resolve the tension.
- 8Z: RUWE dose-response -> the monolithic 0.2 smear DISSOLVES: cleanest quartile needs ZERO width, upper three plateau at ~0.1, and pooling strata manufactures ~2x width (one noise dial for mixed-quality data fakes blur); floor separate (present everywhere); successor = stratified-noise cube re-fit, alpha exposure = the number.
- 9A: Stratified-noise alpha re-fit -> the alpha measurement SURVIVES the width repair: honest per-stratum model preferred by ~90 lnL, alpha moves +0.0..+0.2 (up, never down), Newton contrast strengthens; worst-RUWE quartile pinned at grid-top noise (flagged).
- 9B: Quality-flag control on GD -> QUALITY-BLIND: Q1 and Q2 gas dwarfs give the same dial (D_Q median 0.000), and the V-ordering replicates inside both strata; fifth control cleared - suspects down to real radius-dependent pressure corrections (external data) or genuine slow-dwarf physics.
- 9C: Radius-dependent pressure correction -> SIGN-INVERTED: the proper outward-growing correction makes the GD dial WORSE (-1.31 -> -1.37 at sigma=15); pressure story dead in-catalog at amplitude AND shape (Rg model infeasible, fallback carried the stage - owned; verdict a-fortiori robust).
- 9D: Q4 robustness pair -> Q4-CARRIED 4/4: dropping the worst-RUWE quartile zeroes block-grade alpha in every law-seed; 9A "de-fanged" reading downgraded; 7I-S + 8W + 9D named as ONE object (fitted-alpha collapses under quality cuts); marginal stratified re-run now REQUIRED.
- 9E: Model-light quality cross-check -> the core is ALIVE in the cleanest quartile (B(Q1) = 1.12, second-highest; anchor 1.078 reproduced as B(all) = 1.0779); Q4 elevation +0.077 not significant; the 9D collapse is a fitted-channel effect on current evidence (GRAY by the letter; marginal re-run arbitrates).
- 9F: THE MARGINAL ARBITER -> Q4-CARRIED 4/4: with all strata alpha_marg = 0.5 (modest Newton contrast), but Q1-Q3 alone vote NEWTON (+13..+20 against alpha >= 0.5); even the cleanest quartile wants ~2x formal errors; pre-signed map EXECUTES: anomaly-real 58 -> 53; grid can't see alpha < 0.5 -> successor = fine-alpha drop scan.
- 9G: Ambient control on GD -> SEVENTH control passed: GD neighborhoods are field-typical (dmed = 0.02 dex, perm P = 0.81) and the dial doesn't order on ambient within GD; the environmental escape is closed in-catalog.
- 9H: Data-provenance manifest (the author's "are you sure?" as an instrument) -> DATA-VERIFIED 15/15 invariants + SHA256 manifest committed; first firing caught the manifest's own file-count-vs-fit-count wiring (owned).
- 9I: Fine-alpha drop scan (power-gated, injected 0.3 recovered 4/4) -> I-SMALL-ALPHA: clean quartiles ALLOW and mildly prefer alpha 0.1-0.3 (+1.5..+6.5 over Newton, 3/4); they reject only the big boost; channels reconcile at deflated amplitude; map executes 53 -> 56; residual = Q1's two channels disagree.
- ROUND 13 (Opus agent, run by us - no relay): combiner math CLEAN ("could not make the combiner lie"); E2 CORRECTION - the 9E double-ratio does NOT control differential noise (wide arm inflates x1.155 vs narrow x1.126; B = un-decomposed sum, counterweight reading withdrawn, values stand); headline re-graded to upper-limit form; credence FROZEN 56; instruments D1/D2/D3 + the GD convergence split adopted.
- 9J: Fine-STD extension (0.7-1.2, 9I tables reused bit-exact) -> J-INTERIOR 4/4: full-strata posterior localizes (0.39-0.79, P(top)=0); the per-stratum alpha ladder is a RUWE dose curve (Q2 0.05-0.11 < Q1 0.12-0.26 < Q3 0.40-0.79 < Q4 0.90-1.05); seed-split 31-vs-101 disclosed.
- 9K: fpm-ceiling curve (round-13 D1, reader-only, same day) -> K-GRAY by letter, direction REAL: capping fpm at 1.5 raises drop-world alpha +0.045..+0.081 in 4/4 (half the fragile bar), and under the cap ALL law-seeds prefer small alpha (BE-101 P(a=0) 0.404 -> 0.026); D2 narrow-pair fpm meter = the decider.
- 9L: THE NARROW-PAIR FPM METER (round-13 D2, THE decider; run-1 window caught by its own boost gate, amendment to 0.2-2 kAU) -> L-NOISE-REAL: the boost-free arm demands E[fpm] = 2.21/2.13 = the same ~2x the joint fit uses, flat across quality quartiles; the noise subtraction is LEGITIMATE, B's pedestal real, the upper-limit alpha reading stands; MAP EXECUTES 56 -> 53; standing measurement: EDR3 pair velocity errors ~2x formal at 0.2-2 kAU.
- 9M: Convergence split on GD (round-13 A4) -> GRAY by letter, artifact POINT-DEAD: converged GD alone reproduces the dial (-1.395) and the sign runs backwards (rising -1.172); DD control shows the mechanism exists elsewhere; eighth control effectively holds (formal hardening blocked by 17-galaxy bootstrap noise).
