# thermal-horizon-rar

**Does gravity's strange behavior at the edges of galaxies follow the thermodynamics
of the universe's own horizon? A research program that tests that idea hard — run as
an open human + AI collaboration.**

*(Numbers below are current as of the notebook's 2026-07-24 state; [PAPER.md](PAPER.md)
is the source of truth. This is an active research program, not a settled result, and
it has not yet been peer-reviewed.)*

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
[LEDGER.csv](LEDGER.csv)):

- Wide binaries orbit **~8% faster than Newton predicts** (velocity boost 1.078,
  CI 1.052–1.103, after a perspective correction we adopted from engaging a critical
  paper). Newton loses in **all 2000** bootstrap contests we ran.
- The binary data **demand the "external field effect"** — the suppression of the
  anomaly by the galaxy's ambient field. That effect violates the strong equivalence
  principle, which means **no dark-matter model can produce it even in principle.**
- The deep-regime expansion of the measured curve carries a first coefficient
  c₁ ≈ ½ — the vacuum's half-quantum — and **two completely disconnected systems
  (galaxies and binaries) measure the same value.** A classical (non-quantum) thermal
  bath is excluded at the level of hundreds of log-likelihood units; the vacuum "+1"
  in the response structure is measurably load-bearing.
- The one function that fits *both* systems — the "ambient-gated bath" — was not
  fitted into existence. Its pieces were **derived**: the gate is the Boltzmann cost
  of borrowing quanta from the environment's thermal cloud, the borrow count comes
  from the interaction's own loop order, and exact quantum mechanics reproduces the
  lending probability to a fraction of a percent. It postdicted both systems' tail
  exponents before they were measured as a pair.
- The core mechanism is even **buildable on a lab bench** (superconducting-circuit
  parameters worked out; it would validate the mechanism class, not gravity itself —
  we're explicit about that distinction).

## What we do NOT claim

- **The Saturn problem is real and carried openly.** In the field formulation, our
  rule (like every relative of MOND we tested) predicts a solar-system quadrupole
  ~4× above the Cassini bound. The derivation itself points at the one known escape
  (a trajectory-level coupling), and making that quantitative is the current work —
  but as of today this is an open tension, stated in the paper's own abstract.
- **Galaxy clusters are not addressed.**
- **Dark matter is not "disproven."** On galaxies alone, dark-matter-plus-feedback
  can mimic a lot. The no-mimic content is the *cross-system lock* — the binaries,
  the shared coefficients, the shared temperature scale. We say exactly this and no
  more.
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
- **Fourteen corrections and counting.** [PAPER.md](PAPER.md) Appendix A logs every
  retraction and correction, including the embarrassing ones (a claimed priority
  that a careful read revoked; three hallucinated citations caught by
  primary-source rules; a headline significance that deflated 4× under an honest
  error model). Wrong things die fast and publicly here.
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

Start reading with [PAPER.md](PAPER.md) (the manuscript),
[NOTES-horizon-inertia.md](NOTES-horizon-inertia.md) (the chronological lab
notebook, retractions included), and [TODO.md](TODO.md) (the live priority queue).

## License

CC-BY-4.0. Use it, check it, break it — an issue that survives our gates is a
contribution.
