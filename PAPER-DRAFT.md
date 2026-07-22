# The Radial Acceleration Relation as a Thermal Occupation Law
## (working draft — priority check pending; do not circulate)

**Candidate titles:**
1. *The radial acceleration relation as a Bose–Einstein occupation law: a joint
   galactic, solar-system, and wide-binary test*
2. *A measured screening index for low-acceleration gravity, and its thermal reading*

---

### Abstract (draft)

The interpolating function preferred by the radial acceleration relation (RAR),
ν(y) = [1−e^(−√y)]⁻¹, is algebraically identical to 1 + n_BE, where n_BE is a
Bose–Einstein occupation number in the variable x = √(g_N/a₀) = r_M/r. In this
reading, the deep-MOND regime is the Rayleigh–Jeans (equipartition) limit of a
thermally occupied mode and the Newtonian regime is its Wien tail; the RAR is a
Planck curve. We test this reading three ways. (1) We introduce a one-parameter
screening family ν_p ≡ (1−e^(−y^p))^(−1/2p) that contains the occupation law at
p = ½, and measure p = 0.443 (+0.063/−0.050) from 2,700 SPARC points, with
Cassini radio tracking independently requiring p > 0.234: two datasets eight
decades apart in acceleration select overlapping regions of one parameter, and
the thermal value p = ½ lies within 1.1σ. (2) We measure the wide-binary
velocity excess in the Gaia EDR3 catalog (14,071 pairs after cuts): a separation-
dependent boost reaching 1.202 ± 0.058 at 20–50 kAU, excluding Newtonian
dynamics at ≈7σ in our pipeline, with hierarchical-triple contamination bounded
≲5% by distribution-tail shape and RUWE insensitivity. (3) The occupation
reading fixes how an external field must enter (inside the mode energy); this
naive embedding is disfavored (Δχ² ≈ 5) relative to the full AQUAL external-
field solution, indicating the reading requires a two-field completion. The
framework makes falsifiable commitments: an NLO RAR coefficient of exactly ½,
and a₀ ∝ H(z) if the occupying bath is the instantaneous horizon.

### Sections → existing assets

1. **Introduction** — anomaly census framing; the a₀–horizon coincidences
   ([calcs/coincidences.py](calcs/coincidences.py)).
2. **The identity** — exact algebra + Planck-curve reading; NLO series
   (verified symbolically). *Priority sweep (2 agents, 2026-07-22): the identity is
   not explicitly published. Milgrom corpus: phenomenological only. McGaugh+16:
   empirical fit, no thermal reading. Desmond+23 symbolic regression: no BE remark.
   Nearest neighbors, must cite: Pazy & Argaman (PRD 85, 104021; PRD 87, 084063) —
   quantum statistics with freeze-out on holographic screens, different function;
   q-deformed heat-capacity derivations (arXiv:2010.03530); Smolin (PRD 96, 083523);
   Ho–Minic–Ng. CAVEAT: Haiku-grade sweep, not exhaustive — a citation-graph walk
   of McGaugh+16 is due before submission; frame as "we note that..." regardless.*

#### Section 2 draft text

> The function selected by the radial acceleration relation admits an exact
> rewriting. With x ≡ √(g_bar/a₀), the McGaugh–Lelli–Schombert form
> g_obs = g_bar/(1−e^(−x²)^(1/2))... [notation: ν(y) = (1−e^(−√y))⁻¹] satisfies
>
>   ν = 1 + 1/(e^x − 1) ≡ 1 + n_BE(x),
>
> where n_BE is the Bose–Einstein occupation of a mode of energy ε at temperature
> T with ε/k_BT = x. The observed acceleration is therefore
>
>   g_obs = g_bar·[1 + n_BE(x)],  x = √(g_bar/a₀) = r_M/r,
>
> with r_M ≡ √(GM/a₀) the MOND radius. Three structural consequences follow.
> (i) The deep-MOND limit is the Rayleigh–Jeans regime: n_BE → 1/x yields
> g_obs → √(g_bar·a₀), i.e. the scale invariance of the deep limit is the
> classical equipartition of the occupied mode. (ii) The Newtonian limit is the
> Wien regime: the anomalous component freezes out as e^(−x), consistent with
> the solar-system bounds that exclude power-law tails. (iii) The expansion
> g_obs = √(g_bar·a₀) + g_bar/2 + O(x·g_bar) fixes the next-to-leading
> coefficient at exactly ½ — a parameter-free target for precision RAR fits.
> The mode variable has a natural reading: ε(r) is the Unruh quantum of the
> deep-MOND acceleration √(g_bar·a₀), and the bath temperature corresponds to
> the acceleration a₀ ≈ cH₀/2π, i.e. the de Sitter horizon. If the bath is the
> instantaneous horizon, a₀ ∝ H(z); if it is set by Λ alone, a₀ is constant —
> a dichotomy resolvable with high-redshift rotation curves. We emphasize what
> is and is not claimed: the identity is exact algebra about an empirical fit;
> whether nature computes n_BE is the hypothesis the remainder of this paper
> constrains — the measured screening index p = 0.443(+0.063/−0.050) is
> consistent at 1.1σ with the thermal value p = ½, and the naive embedding of
> an external field inside the mode energy is disfavored (Δχ² ≈ 5) relative to
> the full AQUAL solution, indicating any completion requires genuine two-field
> structure rather than a scalar shortcut.
3. **The screening index** — ν_p family; SPARC fit + bootstrap; Cassini bound;
   joint constraint figure ([calcs/sparc_rar_fit.py](calcs/sparc_rar_fit.py)).
4. **Wide binaries** — EDR3 pipeline, cuts, ṽ statistic; GPU forward model
   (conservative law, self-consistent ICs — [calcs/stage2b_population.py](calcs/stage2b_population.py),
   [calcs/stage2e_refined.py](calcs/stage2e_refined.py), [calcs/stage2f_nup_test.py](calcs/stage2f_nup_test.py));
   triple bounds (tails + RUWE — [calcs/stage2d_ruwe_variant.py](calcs/stage2d_ruwe_variant.py));
   verdict table (Newton χ²=55; AQUAL 2.0; ν_p 3.6).
5. **The EFE-embedding test** — BE-naive vs AQUAL (χ² 7.2 vs 2.0) — the data
   discriminate embeddings ([calcs/stage2g_be_efe.py](calcs/stage2g_be_efe.py)).
6. **Predictions & falsifiers** — NLO ½; a₀(z) ∝ H(z) (high-z Tully–Fisher);
   occupation-fluctuation noise (to be computed); DR4 family discrimination.
7. **Systematics & bias autopsy** — the three self-caught errors (axial folding
   2/π fingerprint; non-conservative EFE recipe energy pumping; IC and a₀-
   convention biases) as a transparency appendix.
8. **The Banik tension** — stated openly; reconciliation-on-one-sample as the
   necessary next work.

### To do before any submission
- [ ] Priority sweep verdict (two agents out)
- [ ] M/L marginalization on the p measurement (widens error bar — quantify)
- [ ] Proper two-field BE theory or drop to "phenomenological reading"
- [ ] Chance-alignment robustness in the 20–50 kAU bin (N=214)
- [ ] Reproduce Banik-style statistic on our sample
- [ ] Occupation-fluctuation noise estimate
