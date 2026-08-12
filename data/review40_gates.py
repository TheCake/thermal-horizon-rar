"""ROUND 40 adoption gates for papers/paper1_wide_binaries.md.

Run:  py data/review40_gates.py > data/review40_gates.txt
All gates must PASS before the draft is considered current.
"""
import re
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

PAPER = Path(__file__).resolve().parents[1] / "papers" / "paper1_wide_binaries.md"
TXT = PAPER.read_text(encoding="utf-8")
LINES = TXT.split("\n")

results = []


def gate(name, ok, detail=""):
    results.append((name, bool(ok), detail))
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  |  {detail}" if detail else ""))


# ---------------------------------------------------------------- structure --
abstract = TXT.split("## Abstract")[1].split("## 1.")[0].strip()
n_abs = len(abstract.split())
gate("G1 abstract <= 250 words", n_abs <= 250, f"{n_abs} words")


def prose_lines():
    """Main-body prose: drop headings, tables, images, italic-only notes."""
    out, in_refs = [], False
    for ln in LINES:
        s = ln.strip()
        if s.startswith("## References"):
            in_refs = True
        if in_refs or not s:
            continue
        if s.startswith(("#", "|", "!", "- ", "*")):
            continue
        out.append(s)
    return out


def sentences(text):
    text = re.sub(r"\((?:[^()]*)\)", lambda m: m.group(0).replace(".", "․"), text)
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z“W])", text)
    return [p.replace("․", ".") for p in parts if p.strip()]


sents = []
for ln in prose_lines():
    sents.extend(sentences(ln))
lens = [len(s.split()) for s in sents]
worst = max(range(len(lens)), key=lambda i: lens[i])
mean = sum(lens) / len(lens)
gate("G2 no prose sentence >= 70 words", max(lens) < 70,
     f"n={len(lens)} mean={mean:.1f} max={max(lens)}: '{sents[worst][:70]}...'")

# ------------------------------------------------------------ number checks --
checks = [
    ("G3 fitted-cell band expectation (M1)", ["26.2"]),
    ("G4 model-light Q1 boost (M3)", ["1.1216"]),
    ("G5 model population size (M4)", ["5 × 10⁵"]),
    ("G6 stratified noise posterior (M2)", ["1.97–2.25"]),
    ("G7 band count P(<=9) (M1)", ["9.8 × 10⁻⁵"]),
    ("G8 differential-noise decomposition (M3)", ["1.155", "1.126"]),
    ("G9 mass-model floor wording (m6)", ["never below 1.77"]),
    ("G10 R_chance span (m2)", ["7 × 10⁻⁶", "6 × 10⁻³"]),
    ("G11 overshoot attribution (m3)", ["17 pairs there at the fitted cell"]),
    ("G12 AQUAL field ratio (m7)", ["roughly 65%"]),
    ("G13 sensitivity band endpoint (m8)", ["0.67–0.74"]),
    ("G14 Q1 width value (m10)", ["0.03"]),
    ("G15 heavy-tail class bound (m21)", ["2.7 against the nine"]),
    ("G16 spectroscopic coverage (m20)", ["three of the twenty-two"]),
    ("G17 NSS coverage disclosure (m19)", ["Seven of the nine"]),
    ("G18 secular-mixing staleness item (m22)", ["−0.90", "stage10o_p15erode.py"]),
    ("G19 repository status (M7a)", ["private at the time of writing"]),
    ("G20 companion-note declaration (M7b)", ["observational note submitted separately"]),
]
for name, needles in checks:
    missing = [n for n in needles if n not in TXT]
    gate(name, not missing, "missing: " + ", ".join(missing) if missing else "ok")

# census median RUWE in its own context (m1)
m = re.search(r"median maximum RUWE ([\d.]+) in the corrected nine-pair census, ([\d.]+) in the uncorrected", TXT)
gate("G21 census median max RUWE (m1)", bool(m) and m.group(1) == "1.12" and m.group(2) == "1.06",
     m.group(0) if m else "context sentence not found")

# Lindegren cited in body, not only in the reference list (M5)
body = TXT.split("## References")[0]
lind = [ln for ln in body.split("\n") if "Lindegren" in ln and not ln.strip().startswith("*Draft")]
gate("G22 Lindegren cited in body (M5/m13)", len(lind) >= 2, f"{len(lind)} body mentions")
gate("G23 error envelope defined (M5)",
     "By that envelope we mean a proper-motion inflation of at most 1.8" in TXT,
     "definition at first use")

# other orphan references (m13)
gate("G24 Bekenstein & Milgrom cited (m13)", "Bekenstein & Milgrom 1984)" in body, "")
gate("G25 Pittordis & Sutherland 2023 cited (m13)", "Pittordis & Sutherland 2018, 2023" in body, "")

# ------------------------------------------------------- alpha-sector sites --
SITES = {
    "abstract": abstract,
    "S1 introduction": TXT.split("## 1. Introduction")[1].split("## 2.")[0],
    "S6.3 fine scan + closing": TXT.split("### 6.3")[1].split("### 6.4")[0],
    "Table 4 row": [ln for ln in LINES if ln.startswith("| clean strata (Q1–Q3)")][0],
    "Conclusions 5": [ln for ln in LINES if ln.startswith("5. After quality stratification")][0],
}
for label, chunk in SITES.items():
    has_excl = ("α ≥ 0.5" in chunk) or ("amplitudes α ≥ 0.5" in chunk)
    has_allow = "0.3–0.5" in chunk
    has_interior = "0.1–0.3" in chunk
    conditioned = any(w in chunk for w in ("low-conversion", "dissolves", "prior"))
    gate(f"G26 alpha sector consistent [{label}]",
         has_excl and has_allow and has_interior and conditioned,
         f"excl={has_excl} allow={has_allow} interior={has_interior} conditioned={conditioned}")

gate("G27 retired flat wording absent",
     "allow, and mildly prefer" not in TXT, "the round-13 re-grade wording")

# ------------------------------------------------------------ register bans --
banned = "anti" + "gravity"
n_banned = len(re.findall(banned, TXT, flags=re.I))
gate("G28 banned term absent", n_banned == 0, f"{n_banned} occurrences")
for code in ("companion-win", "boost-win", "the phantom veto", "Not for circulation"):
    gate(f"G29 codename absent [{code}]", code.lower() not in TXT.lower(), "")

# figure appearance order (m12)
order = [int(m.group(1)) for m in re.finditer(r"^!\[Figure (\d)\]", TXT, flags=re.M)]
gate("G30 figure numbering in appearance order (m12)", order == sorted(order), str(order))

# ------------------------------------------------------------------ summary --
n_fail = sum(1 for _, ok, _ in results if not ok)
print()
print(f"GATES: {len(results) - n_fail}/{len(results)} PASS")
print(f"abstract words = {n_abs}; prose sentences = {len(lens)}, mean {mean:.1f}, max {max(lens)}")
sys.exit(1 if n_fail else 0)
