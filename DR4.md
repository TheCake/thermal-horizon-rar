# DR4.md — the day-one plan (frozen 2026-08-06, before the data exists)

Gaia DR4 is expected around the end of 2026. This document freezes WHAT
gets tested, in WHAT order, with WHICH existing instruments, so the
first-week analysis is mechanical. The kill conditions live in
[PREDICTIONS.md](PREDICTIONS.md) (registered entries, never edited);
this file only points at them and sequences the work. Nothing here is a
new claim.

## The frozen test set (pointers, not restatements)

| # | Registered entry | What DR4 provides | The instrument that runs |
|---|---|---|---|
| 1 | P9 (width-object discriminator; the first read) | improved epoch astrometry → honest per-pair errors | the 9F/9L fpm meter re-run with DR4 errors: if the noise chase quenches AND sq dies together = error tail; if sq ≈ 0.2 persists while fpm → ≤1.4 = astrophysical width |
| 2 | The pair-error bound itself (the RNAAS note number) | independent re-measurement of the ≥2.1× at 0.2–2 kAU | same meter; the note's number is confirmed, revised, or retired in one run |
| 3 | The operative upper limit (Paper 1 §6) | cleaner strata, more pairs | the marginalized posterior (7J-z machinery) on the DR4 catalog; the D3 promotion gate still bounds any stronger-than-limit language |
| 4 | The census prediction (Paper 1 §8) | ~10× effective clean sample | the perpendicular census re-count; band/cliff structure predicted to sharpen if real, dissolve if noise |
| 5 | P7 (MI eccentricity signature) | eccentricity-resolved orbits | e-resolved fits: MI-vs-MG separation where 4L tied |
| 6 | P4 (weak-ambient pair) | weak-ambient binaries resolvable | tail-index read at e_N ≈ 0.4 (Δp ≈ 0.025 convention split) |
| 7 | P2 (environmental ordering, out-of-sample form) | environment-resolved samples | the ordering test; failure = REAL strike, pre-stated |
| 8 | P6 (the γ discriminator, currently unmeasurable) | DR4-era distance anchors | the 7A instrument on anchored subsamples, if depth allows |

Order rationale: 1–2 are instrument-level and gravity-free (run first;
they recalibrate everything downstream); 3–4 are the paper-level
claims; 5–8 are the theory-facing legs.

## Build list (do before DR4; none needs the data)

- [ ] DR4 loader: the archive schema will differ from the EDR3 binary
  catalog (El-Badry-style pair construction may need re-running or a
  new published catalog adopted; scout for community DR4 wide-binary
  catalogs at release). Build against an EDR3-shaped mock first; the
  9F loader's cut ladder is the spec (Paper 1 Table 1).
- [ ] Error-model adapter: DR4 epoch astrometry may publish per-epoch
  residuals; decide the formal-error convention BEFORE seeing the
  kinematics (pre-commit the convention choice when the data model is
  published — that document precedes the data itself).
- [ ] The P9 read script: a thin reader over the 9L meter with the
  DR4 error columns; bars pre-registered in P9 already.
- [ ] The census re-count script: 4J/7K-b machinery with the DR4
  sample; the (band, cliff) pair convention frozen (corrected
  convention primary, both reported).
- [ ] Freeze the strata definition: RUWE's DR4 analogue (or RUWE
  itself if published) — decide and pre-commit when the data model
  documentation lands, before kinematics.

## The decision tree (committed in advance)

- P9 → error tail: the width object was noise; the operative band's
  noise-diluted reading was RIGHT to be conservative; the PHYS-envelope
  conditional (α = 0.80) dies; the upper limit tightens.
- P9 → astrophysical: the width object is real population structure;
  the sq = 0.2 sector gets a physical identity hunt; the noise ceiling
  reverts toward Lindegren and the PHYS conditional gains standing
  (promotion still gated by D3).
- Census sharpens (band repopulates at ~10× with the cliff intact):
  the strongest single object in the program; the D3 arm becomes
  worth running.
- Census dissolves: the band was tail noise; the census retires with
  a public note (the 7K-b self-refuting-cliff logic already bounds
  this direction).
- Any P2/P4/P7 outcome: moves by the registered kill conditions only.

## Standing rules that survive into the DR4 era

Pre-reg before any run; bars locked before data; the 6-gate audit
after any refetch (stage9h manifest); no credence movement outside
pre-signed maps; NOTES/LOG/LEDGER discipline unchanged. The in-sample
function freeze (PREDICTIONS §0) remains in force: DR4 is
out-of-sample territory, where function-level claims are allowed to
move again.
