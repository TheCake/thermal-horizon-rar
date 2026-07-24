# The measurement ledger

`LEDGER.csv` is the program's curated table of every headline number: one row per
quoted quantity, with provenance (stage, script, output file, data dependence) and
**supersession discipline** — superseded and retracted values stay in the table,
visibly marked, so they can never be silently resurrected by a later reader (or a
later fit). This file documents the format; `calcs/stage6q_worldtable.py` audits the
ledger mechanically and assembles the world table (every candidate law scored
against every ledger test).

## Rules

1. **Hand-curated only.** Rows are written by hand from NOTES, never auto-scraped
   from stage outputs — auto-scraping is exactly the mechanism that would resurrect
   a superseded number (e.g. the wrong-convention alpha = 1.54).
2. **The value string preserves the NOTES convention** (sign conventions differ
   between galaxy objectives and binary lnL tables), so every value carries
   direction words ("leads BE by...", "trails BE by...") rather than bare signs.
3. **Status semantics:**
   - `CURRENT` — quote freely; this is the number.
   - `CO-QUOTED` — valid but must be quoted together with its partner row
     (`superseded_by` points at the partner): raw/corrected boost, raw/honest NLO,
     flat/hier c1.
   - `SUPERSEDED` — kept for the record; `superseded_by` points at the replacement.
   - `RETRACTED` — the claim itself was withdrawn (Appendix A of PAPER.md maps
     the corrections); the row records what was wrong and where that was logged.
4. **Every row cites its script**; the audit verifies the script (and output file,
   where given) exists on disk, that ids are unique, and that supersession pointers
   resolve to live rows.
5. **Adding a row that replaces another** means flipping the old row's status and
   pointing it here — never deleting it.

## Columns

`id` — stable slug; `status` — see above; `system` — SPARC / EDR3 / lensing /
solar / Chae21 / meta; `quantity` — what was measured; `value` — the number(s),
NOTES convention, direction words; `stage` — NOTES stage tag(s); `script` —
primary script; `output` — data/ output file; `data_deps` — datasets the number
rests on (rows sharing a dataset are NOT independent votes — the world table
groups them); `superseded_by` — target id; `note` — one line.

## Regenerating the world table

```
py calcs\stage6q_worldtable.py
```

writes `data/stage6q_worldtable.txt`: audit gates first (G1 ids/status, G2
provenance-on-disk, G3 supersession resolution, G4 current-quantity collisions,
G5 spot-check greps of ledger values against stage outputs), then the
function-by-test grade matrix with per-cell ledger provenance, the viability
summary, and the data-overlap groups.
