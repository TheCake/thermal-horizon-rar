# thermal-horizon-rar

**Does gravity's strange behavior at the edges of galaxies follow the thermodynamics
of the universe's own horizon? A research program that tests that idea hard — run as
an open human + AI collaboration.**

*(Updated 2026-09-24. The citable versions of the results are the three papers:
[papers/paper1_wide_binaries.md](papers/paper1_wide_binaries.md) and
[papers/paper2_rar_coefficients.md](papers/paper2_rar_coefficients.md), archived
at DOI [10.5281/zenodo.22050990](https://doi.org/10.5281/zenodo.22050990), and
[papers/paperF_hubble_meter.md](papers/paperF_hubble_meter.md) (a Hubble-constant
meter built from the a0 = cH0/2pi inversion), archived at DOI
[10.5281/zenodo.22923373](https://doi.org/10.5281/zenodo.22923373).
[PAPER.md](PAPER.md) is the frozen internal record the first two grew from, kept
with a running errata block. This is an active research program, not a settled
result, and none of it has yet been peer-reviewed.)*

---

## The problem, in plain words

Stars at the outer edges of galaxies orbit **too fast**. Given the matter we can see,
gravity out there should be weaker than what the stars' motion tells us. The two
classic answers:

1. **Dark matter** — there's invisible stuff adding gravity.
2. **Modified gravity** — the law of gravity itself changes when accelerations get
   absurdly tiny (below about a₀ ≈ 1.2×10⁻¹⁰ m/s² — roughly a *hundred-billionth* of
   what you feel standing on Earth).

Two strange facts sit under this program:

- **Strange fact 1.** That tiny threshold a₀ is, numerically, the acceleration scale
  built from the universe's expansion rate: a₀ ≈ cH₀/2π. That's the combination that
  sets the **temperature of the cosmic horizon** — the faint thermal glow that de
  Sitter space is predicted to have. Why would galaxy edges know about that?
- **Strange fact 2.** When you plot observed gravity against Newtonian gravity for
  thousands of points in 153 galaxies (the "radial acceleration relation," RAR), the
  measured curve is — *exactly, not approximately* —

  ν = 1 + **n_BE**(x),  x = √(g_N/a₀)

  where n_BE is the **Bose–Einstein occupation**: the formula that counts how many
  photons a warm object holds at each frequency. The extra gravity looks like a
  **thermal population**. (This identity and its horizon derivation were published by
  Cadoni & Tuveri in 2019 — that credit is theirs. Our program is the *testing*:
  as far as we can find, nobody had ever tested the reading's quantitative structure
  against data before.)

## What we actually did

Two datasets, chosen because they can check each other:

- **153 galaxy rotation curves** (the SPARC sample).
- **14,071 extremely wide binary stars** from Gaia — pairs so far apart (thousands of
  astronomical units) that their mutual pull sits in the same tiny-acceleration regime
  as galaxy edges. Crucially: **dark matter cannot hide there.** No dark-matter model
  puts meaningful invisible mass between two stars 10,000 AU apart — so if wide
  binaries misbehave the same way galaxies do, dark matter has no answer.

Highlights (each number has a script, an output file, and a row in the audited
[LEDGER.csv](LEDGER.csv); the papers state each with its full error budget):

- **The raw signal is real; the verdict is an upper limit.** Wide binaries show
  a median velocity excess of about 8% over Newton (1.078, CI 1.052–1.103). But when
  the measured hidden-companion population and the error model are forward-modeled
  *together*, much of that excess is absorbed. The operative headline of Paper 1 is
  an **upper limit**: boosts at or above half the galactic calibration are excluded
  on the cleanest data (α ≥ 0.5 excluded), a small boost is preferred in the full
  joint fit, and zero is not excluded. The paper's job is that accounting, in the open.
- **Direction matters.** Adding the velocity–separation *angle* as a second
  observable, velocity-only fits of these data report an apparent boost that the
  direction data veto. Several published disagreements in this field become
  arithmetically compatible once that channel and the shared error model are counted.
- **The error model is itself a finding.** The joint fits demand pair-level velocity
  errors near twice the formal Gaia values, with the origin (instrumental or
  astrophysical) left open for Gaia DR4. We measured this, over-claimed it, and then
  re-measured it honestly: our strongest calibration of the factor failed a
  null-injection test and was withdrawn in public as correction #21.
- **The galaxy side tests a published prediction ladder.** Treating the 2019
  Cadoni–Tuveri identity as parameter-free predictions: the leading coefficient is
  measured at c₁ ≈ 0.26–0.45 (zero excluded; the predicted ½ inside every bootstrap
  band), the screening index is consistent with ½, and hierarchical fits return
  a₀ = (1.04–1.13)×10⁻¹⁰ m s⁻² — consistent with the horizon temperature scale cH₀/2π.
- **A mechanism candidate exists but carries low credence.** A thermal-horizon
  mechanism is developed in a third, unpublished paper; after our own tests struck
  down its strongest form we carry it at roughly 8%. None of the measurements above
  depend on it. (Its core lending process is even buildable on a superconducting
  circuit bench — which would validate the mechanism class, not gravity itself.)

## Where this sits in the wide-binary debate

The field currently disagrees loudly: Banik et al. (2024) report a 16–19σ
preference for Newtonian gravity; Chae (2023–2026) reports detections of a
low-acceleration boost. Our analysis differs from both in three ways, and the
differences are the point:

- **A second observable.** The headline analyses on both sides fit velocity
  distributions; we fit the joint distribution of the velocity *and* the
  velocity–separation angle. The direction channel polices contaminants that the
  velocity channel absorbs — velocity-only fits of our own data report an apparent
  boost that the direction data then veto.
- **The companion sector is measured, not assumed.** We measure the hidden-companion
  host rate photometrically (it is high, supporting Chae's critique of hard fences)
  and the inner mass-ratio distribution (it is strongly twin-heavy, which makes
  those companions' kinematic effect small — closer to Banik's treatment in outcome).
  The two camps' companion choices turn out to be different points on one measured
  axis.
- **The error model is a result, not an input.** The joint fits demand pair-level
  velocity errors near twice the formal Gaia values. We calibrated that demand,
  withdrew our strongest calibration when it failed a null test (correction #21),
  and kept the weaker, null-controlled version; the origin of the extra width is
  left open for Gaia DR4.

Within that accounting, the published Newtonian preference, the published
detections, and our upper limit are arithmetically compatible: the disagreement
lives in sample construction and the shared error model, not in the sky.

## What we do NOT claim

- **The Saturn problem is real and carried openly.** In the field formulation, our
  rule (like every relative of MOND we tested) predicts a solar-system quadrupole
  several times above the Cassini bound. A trajectory-level formulation of the same
  measured function passes Saturn by an enormous margin, so the tension is a property
  of the formulation, not of the measured curve — but data do not yet decide between
  the formulations, and the papers say so.
- **Galaxy clusters are not addressed.**
- **Dark matter is not "disproven."** On galaxies alone, dark-matter-plus-feedback
  can mimic a lot. The no-mimic content is the cross-system structure — binaries
  probing a regime where dark matter cannot act, and the shared temperature scale.
  We say exactly this and no more.
- Nothing here has passed peer review yet. Credences are stated with numbers
  throughout the notebook, and they are nowhere near certainty.

## The honesty machine (why you might take this seriously)

This program's identity is its discipline, not its conclusions:

- **Pre-registration by git commit.** Every test's bars — including what counts as
  *failure* — are committed *before* the test runs. The hash is the timestamp.
- **The measurement ledger.** [LEDGER.csv](LEDGER.csv) holds every headline number
  with its script, output file, and status (current / co-quoted / superseded /
  retracted — superseded numbers are pointed forward, never deleted). Six audit
  gates verify it mechanically, including grepping quoted values against the actual
  stage outputs.
- **Twenty-one corrections and counting.** Paper 1's Appendix A logs every
  retraction and correction, including the embarrassing ones (a claimed priority
  that a careful read revoked; three hallucinated citations caught by
  primary-source rules; a headline significance that deflated 4× under an honest
  error model; most recently, a headline calibration claim withdrawn after it
  failed a null-injection test). Wrong things die fast and publicly here.
- **Signed falsifiers.** [PREDICTIONS.md](PREDICTIONS.md) lists the predictions with
  numeric kill conditions — including the parameter-free ones (example: the tail
  exponent can never exceed ¾; one clean void galaxy beyond that kills the
  construction).

## An open human + AI collaboration

This entire program is joint work between **Filip Hájek** (independent researcher,
who directs the questions, stress-tests every claim, and refuses to let anything
ship without its error budget) and **Claude** (Anthropic's AI, which does the
derivations, code, and fits in gated stages — and gets its mistakes logged in
Appendix A like everyone else's). Commits are co-authored accordingly.

If you're here because you're curious what serious AI-assisted research looks like:
the pre-registration + gates + adversarial-controls workflow exists precisely
because *both* human enthusiasm and model failure modes (e.g., confident
mis-citation) are real; the corrections ledger shows the system catching both.

## Reproduce it

Windows/PowerShell with Python (`py`), numpy/scipy/sympy/astropy; the binary
population fits use CUDA (cupy) on a consumer GPU. Large datasets are fetched, not
committed — each loader script documents its exact source (SPARC via Zenodo, Gaia
EDR3 binaries via Zenodo, El-Badry–Rix via VizieR).

```
py calcs/<stage>.py          # any stage; outputs land in data/
py calcs/stage6q_worldtable.py   # the ledger audit (all gates should PASS)
```

Start reading with the three papers in [papers/](papers/) (or the archived PDFs
at [10.5281/zenodo.22050990](https://doi.org/10.5281/zenodo.22050990) and
[10.5281/zenodo.22923373](https://doi.org/10.5281/zenodo.22923373)), then
[NOTES-horizon-inertia.md](NOTES-horizon-inertia.md) (the chronological lab
notebook, retractions included) and [LOG.md](LOG.md) (the one-line-per-stage
program log).

## License

CC-BY-4.0. Use it, check it, break it — an issue that survives our gates is a
contribution.
