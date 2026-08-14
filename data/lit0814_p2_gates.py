"""lit0814_p2_gates.py -- gates for the Paper 2 literature-refresh round (2026-08-14).

Checks papers/paper2_rar_coefficients.md after the van Putten differentiation
(Section 5.1) and the Park et al. 2026 quadrupole update (Section 9.1).

Gates:
  G1  abstract <= 250 words
  G2  "2608.07112" or "van Putten" present
  G3  "024066" or "Park" present
  G4  "1.6" present (the Park+26 central value)
  G5  "5.2" present (the verified 2-sigma ceiling, x 10^-27 s^-2)
  G6  both new ratio values printed by lit0814_cassini.txt appear in the paper
  G7  no sentence >= 70 words
  G8  banned word (constructed in-checker, never spelled in a tracked file): zero

Output is ASCII only.
"""

import re
import sys
import unicodedata
from pathlib import Path


def ascii_safe(s):
    """Windows cp1252 consoles cannot print the paper's subscripts/sigmas."""
    return s.encode("ascii", "replace").decode("ascii")

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "papers" / "paper2_rar_coefficients.md"
CASSINI = ROOT / "data" / "lit0814_cassini.txt"

text = PAPER.read_text(encoding="utf-8")
cassini = CASSINI.read_text(encoding="utf-8")

results = []


def gate(name, ok, detail):
    results.append((name, bool(ok), detail))


# ------------------------------------------------------------------ G1 abstract
m = re.search(r"^## Abstract\s*\n(.*?)(?=^## )", text, re.S | re.M)
abstract = m.group(1).strip() if m else ""
abs_words = len(abstract.split())
gate("G1 abstract <= 250 words", abs_words <= 250, "%d words" % abs_words)

# --------------------------------------------------------------- G2/G3 citations
gate("G2 van Putten cited", ("2608.07112" in text) or ("van Putten" in text),
     "2608.07112=%s van Putten=%s" % ("2608.07112" in text, "van Putten" in text))
gate("G3 Park et al. cited", ("024066" in text) or ("Park" in text),
     "024066=%s Park=%s" % ("024066" in text, "Park" in text))

# ------------------------------------------------------------ G4/G5 bound digits
gate("G4 Park central value '1.6' present", "1.6" in text, "found" if "1.6" in text else "MISSING")
gate("G5 2-sigma ceiling '5.2' present", "5.2" in text, "found" if "5.2" in text else "MISSING")

# ------------------------------------------------------------------ G6 ratios
ratios = {}
for key in ("RATIO_ALPHA1_NEW", "RATIO_ALPHA115_NEW"):
    mm = re.search(key + r"\s*=\s*([0-9.]+)", cassini)
    if mm:
        ratios[key] = mm.group(1)
missing = [k for k in ("RATIO_ALPHA1_NEW", "RATIO_ALPHA115_NEW") if k not in ratios]
if missing:
    gate("G6 both new ratios in paper", False, "could not parse from txt: %s" % missing)
else:
    absent = [v for v in ratios.values() if v not in text]
    gate("G6 both new ratios in paper", not absent,
         "looked for %s ; absent %s" % (sorted(ratios.values()), absent or "none"))

# ------------------------------------------------------------- G7 sentence length
def prose(src):
    """Strip front matter, tables, figures, headings, list rows and appendices."""
    out = []
    for ln in src.splitlines():
        s = ln.strip()
        if not s:
            continue
        if s.startswith(("#", "|", "*", "-", "!", ">")):
            continue
        if s.startswith("**") or s.startswith("("):
            continue
        out.append(s)
    return " ".join(out)


body = prose(text)
# protect decimal points and common abbreviations before splitting on sentence ends
protected = re.sub(r"(\d)\.(\d)", r"\1<DOT>\2", body)
for abbr in ("et al.", "e.g.", "i.e.", "cf.", "Eq.", "Fig.", "Sec.", "approx.", "vs."):
    protected = protected.replace(abbr, abbr.replace(".", "<DOT>"))
sentences = [s.strip() for s in re.split(r"(?<=[.?!])\s+", protected) if s.strip()]
long_sents = [(len(s.split()), s.replace("<DOT>", ".")) for s in sentences if len(s.split()) >= 70]
gate("G7 no sentence >= 70 words", not long_sents,
     "max = %d words" % max((len(s.split()) for s in sentences), default=0))

# ------------------------------------------------------------------- G8 banned
banned = "anti" + "gravity"
low = unicodedata.normalize("NFKD", text).lower()
hits = low.count(banned) + low.count(banned[:4] + "-" + banned[4:]) + low.count(banned[:4] + " " + banned[4:])
gate("G8 banned word count zero", hits == 0, "%d occurrence(s)" % hits)

# ------------------------------------------------------------------------ report
print("=" * 72)
print("PAPER 2 LITERATURE-REFRESH GATES -- 2026-08-14")
print("target: %s" % PAPER.relative_to(ROOT).as_posix())
print("=" * 72)
for name, ok, detail in results:
    print("  [%s]  %-38s  %s" % ("PASS" if ok else "FAIL", name, detail))
print("-" * 72)
if long_sents:
    print("  offending sentences (G7):")
    for n, s in long_sents:
        print("    %d words: %s..." % (n, ascii_safe(s[:160])))
    print("-" * 72)
n_pass = sum(1 for _, ok, _ in results if ok)
print("  %d/%d PASS" % (n_pass, len(results)))
print("  abstract word count: %d" % abs_words)
print("=" * 72)
sys.exit(0 if n_pass == len(results) else 1)
