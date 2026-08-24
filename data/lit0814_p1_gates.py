"""Gates for the 2026-08-14 literature refresh of papers/paper1_wide_binaries.md (draft 0.7).

Run:  py data/lit0814_p1_gates.py > data/lit0814_p1_gates.txt

Gates
  G1  the five 2026 wide-binary papers are cited (Makarov / 2512.25002, Pasquini,
      Saad & Ting, both Chae 2026 entries)
  G2  the Cookson reference carries the journal record (MNRAS 547), not preprint-only
  G3  abstract <= 250 words
  G4  no sentence >= 70 words
  G5  every arXiv identifier newly cited in draft 0.7 was primary-read (listed below,
      each with the deciding line that was checked against the source)
  G6  banned word absent (constructed at runtime, never spelled in this file)
  G7  draft line advanced to 0.7 / 2026-08-14
"""

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "papers" / "paper1_wide_binaries.md"

# Primary reads performed 2026-08-14 (abstract and/or full text fetched from the
# source; the deciding line was quoted and checked, not taken from a search snippet).
PRIMARY_READ = {
    # 2026-08-24 addition (SOL round): Penoyre reference introduced by the
    # draft-0.8 Sec-5.3 rewrite after the wobble-law primary-source audit.
    "2111.10380": (
        "Penoyre, Belokurov & Evans 2022, MNRAS 513, 2437 -- LaTeX source "
        "fetched and read (2026-08-24 audit). Deciding line (abstract, "
        "verbatim): 'the observed UWE scales [proportional to] varpi "
        "(parallax), a (semi-major axis) and Delta = |q-l|/((1+q)(1+l)) "
        "where q and l are the mass and light ratio respectively' -- "
        "algebraically identical to |q/(1+q) - l/(1+l)|; Sec 7: twins "
        "'become slightly harder to detect' (photocentre ~ barycentre)."
    ),
    "2512.25002": (
        "Makarov 2026, AJ 171, 79 -- full-text HTML read. Deciding line: "
        "'Eq. 7 in (Chae 2024) used to estimate the uncertainty of the relative "
        "proper motion magnitude is incorrect. It leads to a strongly underestimated "
        "noise component for the widest separations where the SNR is tending to zero.' "
        "Observable definition read: nu = s_x v_y - s_y v_x (2D cross product)."
    ),
    "2602.04661": (
        "Pasquini et al. 2026, A&A 707, L2 -- abstract read. Deciding line: "
        "'Of the remaining 12, nine can be fitted with bound orbital solution, while "
        "three show velocity differences too large to be reconciled with any bound "
        "Newtonian orbit.' No perpendicular decomposition and no leakage/null "
        "statistic appear; unbinding by encounters/tides is offered by the authors."
    ),
    "2603.11015": (
        "Saad & Ting 2026 -- abstract read. Deciding line: 'Our model yields "
        "gamma = 1.12 (+0.27, -0.22), consistent with Newtonian gravity (gamma = 1) "
        "at the ~0.4 sigma level.' Method: hierarchical Bayesian fit of 3D orbital "
        "elements to the same 36 pairs; de-projected-separation test returns 1.56."
    ),
    "2607.14450": (
        "Chae & Yoon 2026 -- full abstract read. Deciding line: 'proper data quality "
        "control or reasonable variation in multiple-star modeling cannot remove the "
        "low-acceleration gravitational anomaly but confirms the MOND-type "
        "gravitational anomaly'; nulls attributed to bypassed multiple-star "
        "calibration, quality-control bias, and insufficient low-acceleration statistics."
    ),
    "2601.21728": (
        "Chae et al. 2026 -- abstract read. Deciding line: 36 wide binaries with "
        "accurate 3D velocities give gamma = 1.600 (+0.171, -0.141), 'inconsistent "
        "with standard gravity at 4.9 sigma'."
    ),
    "2309.10404": (
        "Chae 2024, ApJ 960, 114 -- record verified (title, journal, volume, article) "
        "as the target of the Makarov uncertainty-formula criticism, so that the "
        "criticism is attributed to the correct Chae 2024 paper."
    ),
    "2602.24035": (
        "Cookson et al. -- journal record verified at the publisher: MNRAS 547, "
        "issue 2, article stag342, doi:10.1093/mnras/stag342 (2026-02-25)."
    ),
}

text = PAPER.read_text(encoding="utf-8")
results = []


def gate(name, ok, detail):
    results.append((name, bool(ok), detail))


# ---- G1: the five new papers are cited -------------------------------------
required = {
    "Makarov (arXiv id)": "2512.25002",
    "Makarov (name)": "Makarov",
    "Pasquini": "Pasquini",
    "Saad & Ting": "Saad",
    "Chae & Yoon": "Chae, K.-H., & Yoon",
    "Chae et al. 2026": "2601.21728",
}
missing = [k for k, v in required.items() if v not in text]
gate("G1 five 2026 papers cited", not missing,
     "missing: " + ", ".join(missing) if missing else
     "all present (" + ", ".join(sorted(required)) + ")")

# ---- G2: Cookson journal record --------------------------------------------
cook_ok = "MNRAS 547" in text and "stag342" in text
gate("G2 Cookson journal record", cook_ok,
     "found 'MNRAS 547, stag342'" if cook_ok else "Cookson still preprint-only")

# ---- G3: abstract word count ------------------------------------------------
m = re.search(r"^## Abstract\s*\n+(.*?)\n+^## ", text, re.S | re.M)
abstract = m.group(1).strip() if m else ""
n_abs = len(abstract.split())
gate("G3 abstract <= 250 words", 0 < n_abs <= 250, f"{n_abs} words")

# ---- G4: sentence length -----------------------------------------------------
# Strip markdown furniture that is not prose: tables, headings, images, code.
prose_lines = []
for line in text.splitlines():
    s = line.strip()
    if not s or s.startswith(("#", "|", "!", "-", ">")):
        continue
    if s.startswith("*Draft") or s.startswith("*Acknowledg"):
        continue
    prose_lines.append(s)
prose = " ".join(prose_lines)
prose = re.sub(r"\(arXiv:[^)]*\)", " ", prose)
# Emphasis markers must not glue a caption's final period to the next paragraph.
prose = prose.replace("*", " ").replace("_", " ")
prose = re.sub(r"\s+", " ", prose)
# Protect abbreviations / decimals from the sentence splitter.
guarded = prose
for abbr in ["et al.", "e.g.", "i.e.", "cf.", "Eq.", "Fig.", "Sect.", "vs.",
             "Dr.", "no.", "ApJ.", "approx."]:
    guarded = guarded.replace(abbr, abbr.replace(".", ""))
guarded = re.sub(r"(?<=\d)\.(?=\d)", "", guarded)
sentences = [s.strip()
             for s in re.split(r"(?<=[.!?])\s+", guarded) if s.strip()]
sentences = [s.strip() for s in sentences]
lengths = [(len(s.split()), s) for s in sentences]
long_sentences = [(n, s) for n, s in lengths if n >= 70]
mean_len = sum(n for n, _ in lengths) / len(lengths) if lengths else 0
gate("G4 no sentence >= 70 words", not long_sentences,
     f"{len(sentences)} sentences, mean {mean_len:.1f}, max "
     f"{max(n for n, _ in lengths) if lengths else 0}"
     + ("" if not long_sentences else
        "; offenders: " + " || ".join(f"[{n}w] {s[:120]}..." for n, s in long_sentences)))

# ---- G5: every newly cited arXiv id was primary-read -------------------------
cited_ids = set(re.findall(r"arXiv:(\d{4}\.\d{4,5})", text))
# "New text" is defined against the committed draft 0.6, not by hand: any arXiv
# identifier this refresh introduces must appear in the primary-read log above.
prev = subprocess.run(
    ["git", "-C", str(ROOT), "show", "HEAD:papers/paper1_wide_binaries.md"],
    capture_output=True, text=True, encoding="utf-8", errors="replace")
if prev.returncode != 0:
    gate("G5 new arXiv ids primary-read", False,
         "could not read HEAD version of the paper: " + prev.stderr.strip()[:200])
else:
    prev_ids = set(re.findall(r"arXiv:(\d{4}\.\d{4,5})", prev.stdout))
    introduced = sorted(cited_ids - prev_ids)
    unread = [i for i in introduced if i not in PRIMARY_READ]
    logged_not_cited = [i for i in sorted(PRIMARY_READ) if i not in cited_ids]
    gate("G5 new arXiv ids primary-read", not unread and not logged_not_cited,
         (f"{len(introduced)} identifier(s) introduced by this refresh "
          f"({', '.join(introduced) if introduced else 'none'}), all primary-read; "
          f"{len(PRIMARY_READ)} entries in the log, all cited in the text"
          if not unread and not logged_not_cited
          else f"introduced but not primary-read: {unread}; "
               f"logged but not cited: {logged_not_cited}"))

# ---- G6: banned word ----------------------------------------------------------
banned = "anti" + "gravity"
count_banned = text.lower().count(banned)
gate("G6 banned word absent", count_banned == 0, f"{count_banned} occurrences")

# ---- G7: draft line ------------------------------------------------------------
draft_ok = "Draft 0.7 (2026-08-14)" in text
gate("G7 draft line advanced", draft_ok,
     "Draft 0.7 (2026-08-14) present" if draft_ok else "draft line not advanced")

# ---- report ---------------------------------------------------------------------
print("LIT-REFRESH GATES -- papers/paper1_wide_binaries.md draft 0.7 (2026-08-14)")
print("=" * 78)
for name, ok, detail in results:
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")
print("=" * 78)
print("\nPRIMARY-READ LOG (arXiv identifiers newly cited in draft 0.7)")
print("-" * 78)
for k in sorted(PRIMARY_READ):
    print(f"  arXiv:{k}\n    {PRIMARY_READ[k]}")
print("-" * 78)
n_pass = sum(1 for _, ok, _ in results if ok)
print(f"\n{n_pass}/{len(results)} gates PASS")
sys.exit(0 if n_pass == len(results) else 1)
