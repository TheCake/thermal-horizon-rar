# REPRO — clone → verify → data → recompute

This is the outsider's path into the repository: how to check that the
record is intact, how to get the data and prove it is byte-for-byte what
the papers used, and how to re-derive the numbers. The repository's rule,
in one line: **every claimed number has a named script under calcs/, a
committed output under data/, and a row in [LEDGER.csv](LEDGER.csv)** —
and every analysis stage was pre-registered by a git commit made before
the run, so the commit hash is the timestamp.

## 0. Environment

- Python 3.13 (3.13.2 is the version of record; any recent 3.x should
  work), then `pip install -r requirements.txt`.
- Commands below are for Windows PowerShell, where Python is `py`. On
  Linux/macOS substitute `python3`. The scripts are OS-agnostic.
- Run everything **from the repository root** (scripts use relative paths).
- GPU: only the wide-binary *population fits* need CUDA (cupy; the runs
  of record used an RTX 5090, where a 10^6-binary population takes
  10–70 s). Everything else — the entire galaxy side, all audits, all
  figures, the Hubble-meter stages — is plain CPU. Without a GPU you can
  still verify every committed output and re-derive every galaxy-side
  number; the binary fits' posterior cubes are committed, so downstream
  reads re-run without refitting.
- Disk: about 1.6 GB after the data fetch.

## 1. Verify the skeleton (no data, under a second)

    py calcs/stage6q_worldtable.py

Six audit gates over the measurement ledger: every ledger row's script
and output exist on disk, supersessions resolve, quoted values are
grepped against the actual committed stage outputs, and the world table
regenerates. All six print PASS on a fresh clone.

## 2. Fetch the data and verify provenance (one command)

    py calcs/fetch_inputs.py

Downloads the two core datasets — the Gaia EDR3 wide-binary catalog
(El-Badry, Rix & Heintz 2021; Zenodo 4435257; 1.4 GB) and the SPARC
release (Lelli, McGaugh & Schombert 2016; Zenodo 16284118) — and then
verifies **every** entry of [data/MANIFEST.sha256](data/MANIFEST.sha256),
the pre-registered provenance manifest, including an aggregate hash over
the 175 SPARC rotation-curve files. `--verify` re-checks without network;
`--selftest` tests the download path without touching data/. Exit code 0
means your inputs are byte-identical to the papers'.

Deeper than hashes:

    py calcs/stage9h_manifest.py

recomputes fifteen historical invariants live from the data — sample
sizes after the exact selection cuts, quality quartiles, the median-boost
statistic, galaxy census counts — and compares each to its
stage-of-record value.

## 3. Regenerate the papers' figures (CPU, with hard gates)

    py calcs/paper1_figures.py
    py calcs/paper2_figures.py
    py calcs/paper3_figures.py
    py calcs/paperf_figures.py

Each rebuilds its paper's figures from the data and the committed stage
outputs, and each carries asserts (gates) that fail loudly if any plotted
number drifts from the record.

## 4. Re-derive headline numbers

The general recipe: find the number in [LEDGER.csv](LEDGER.csv); its row
names the stage script and the output file; run the script and compare
against the committed output. (Population fits draw random realizations;
seeds are fixed and printed, and realization scatter is part of every
quoted error budget.) Worked examples:

- **The Hubble meter (paper: H0 = 65.4 ± 5.0 stat, form band 65.4–70.4):**
  `py calcs/stage10s_h0meter.py` (first light, both arms) and
  `py calcs/stage10t_legregrow.py` (the anchored-arm regrow that sets the
  operative numbers). CPU. The operative verdicts are committed at
  [data/stage10s_verdict.txt](data/stage10s_verdict.txt) and
  [data/stage10t_verdict.txt](data/stage10t_verdict.txt); the ladder
  peg-gearing measurement is `py calcs/round51_gearing.py`.
- **The wide-binary joint fits (paper 1):** the operative marginal lives
  in `calcs/stage7j_marginal.py` and its committed cubes
  (data/stage7j_cube_*.npy, hash-pinned in the manifest). Refitting needs
  cupy/CUDA; reading and re-analyzing the cubes does not.
- Stages that need further small datasets (KiDS lensing, LITTLE THINGS,
  X-COP, ...) document their sources in their own docstrings and fetchers
  (`calcs/fetch_*.py`, `calcs/stage4e_lensing_rar.py`).

## 5. Orientation

- [papers/](papers/) — the manuscripts (the citable PDFs are archived at
  [10.5281/zenodo.22050990](https://doi.org/10.5281/zenodo.22050990) and
  [10.5281/zenodo.22923373](https://doi.org/10.5281/zenodo.22923373)).
- [LEDGER.csv](LEDGER.csv) + [LEDGER.md](LEDGER.md) — every headline
  number: value, stage, script, output, status (current / co-quoted /
  superseded / retracted; superseded rows point forward, never deleted).
- [PREDICTIONS.md](PREDICTIONS.md) — predictions registered before test,
  with numeric kill conditions.
- [NOTES-horizon-inertia.md](NOTES-horizon-inertia.md) — the chronological
  lab notebook, corrections and retractions included;
  [LOG.md](LOG.md) — one line per stage.

If a gate fails on your machine, that is exactly what the gates are for:
please open an issue with the printed diff — an issue that survives our
gates is a contribution.
