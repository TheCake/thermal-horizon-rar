"""ROUND 42 adoption gates for papers/paper3_mechanism.md (draft 0.8).

Checks the register contract (papers/STYLE.md) and the presence of the
round-42 adopted disclosures.  Run:  py data/review42_gates.py
"""
import re
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

PAPER = Path(__file__).resolve().parents[1] / "papers" / "paper3_mechanism.md"
TEXT = PAPER.read_text(encoding="utf-8")

BANNED = "anti" + "gravity"          # never written literally in this repo
results = []


def gate(name, ok, detail):
    results.append((name, bool(ok), detail))


# ---------------------------------------------------------------- helpers
def section_blocks(text):
    """Return {heading: body} for level-2 headings."""
    out, cur, buf = {}, None, []
    for line in text.splitlines():
        if line.startswith("## "):
            if cur is not None:
                out[cur] = "\n".join(buf)
            cur, buf = line[3:].strip(), []
        else:
            buf.append(line)
    if cur is not None:
        out[cur] = "\n".join(buf)
    return out


SECTIONS = section_blocks(TEXT)
ABSTRACT = SECTIONS.get("Abstract", "").strip()


def words(s):
    return [w for w in re.split(r"\s+", s.strip()) if w]


def paragraphs(body):
    # drop sub-headings so they do not glue neighbouring sentences together
    body = "\n\n".join(
        "" if ln.lstrip().startswith("#") else ln for ln in body.splitlines()
    )
    return [p.strip() for p in re.split(r"\n\s*\n", body) if p.strip()]


def prose_paragraphs(body):
    """Paragraphs excluding tables, figure captions, images, list rows."""
    out = []
    for p in paragraphs(body):
        if p.startswith("|") or p.startswith("!") or p.startswith("*Figure"):
            continue
        if p.startswith("- ") or p.startswith("#"):
            continue
        out.append(p)
    return out


def sentences(body):
    out = []
    for p in prose_paragraphs(body):
        blob = re.sub(r"\s+", " ", p)
        parts = re.split(r"(?<=[.?!])\s+(?=[A-Z(\"“])", blob)
        out.extend(s for s in parts if s.strip())
    return out


# ---------------------------------------------------------------- G1 abstract length
n_abs = len(words(ABSTRACT))
gate("G1 abstract word count <= 250", n_abs <= 250, f"{n_abs} words")

# ---------------------------------------------------------------- G2 abstract numerals
# a numeric quantity = a token containing a digit (years, percents, decimals)
num_tokens = re.findall(r"\d[\d.,]*\s*%?", ABSTRACT)
num_tokens = [t.strip() for t in num_tokens]
gate("G2 abstract numeric quantities in [3,5]",
     3 <= len(num_tokens) <= 5,
     f"{len(num_tokens)}: {num_tokens}")

# ---------------------------------------------------------------- G3 sentence length
body_all = "\n\n".join(v for k, v in SECTIONS.items() if k != "References")
longest = (0, "")
for s in sentences(body_all):
    n = len(words(s))
    if n > longest[0]:
        longest = (n, s[:90])
gate("G3 no sentence >= 70 words", longest[0] < 70,
     f"max {longest[0]} words: {longest[1]}...")

# ---------------------------------------------------------------- G4 paragraph length in section 6
sec6_key = [k for k in SECTIONS if k.startswith("6.")]
sec6 = SECTIONS[sec6_key[0]] if sec6_key else ""
p6 = [(len(words(p)), p[:70]) for p in prose_paragraphs(sec6)]
max6 = max(p6) if p6 else (0, "")
gate("G4 max paragraph in section 6 < 250 words", max6[0] < 250,
     f"max {max6[0]} words: {max6[1]}...")
# reported for information: the other long-paragraph sections
allp = [(len(words(p)), k) for k, v in SECTIONS.items()
        for p in prose_paragraphs(v)]
worst = max(allp) if allp else (0, "")

# ---------------------------------------------------------------- G5 required disclosures
required = {
    "sigma_p lower bound": r"≥\s*0\.075",
    "dispersive gate monotonicity": r"opposite monotonicity",
    "threshold family with sigma": r"4/\(2n\s*\+\s*1\s*\+\s*2σ\)",
    "P10 background fraction": r"40\s*(?:–|-|to)\s*68\s*%",
    "search multiplicity": r"multiplicity",
    "post-hoc flag": r"post-hoc flag",
    "clipped bootstrap disclosure": r"clipped lower bound",
    "correction count": r"nineteen corrections",
}
for name, pat in required.items():
    m = re.search(pat, TEXT)
    gate(f"G5 present: {name}", m is not None,
         (m.group(0) if m else "MISSING"))

# ---------------------------------------------------------------- G6 repository not called public
bad = []
for s in sentences(TEXT):
    if re.search(r"repositor", s, re.I) and re.search(r"\bpublic\b", s, re.I):
        bad.append(s[:100])
gate("G6 no sentence calls the repository public", not bad, bad[:2] or "clean")

# ---------------------------------------------------------------- G7 banned word
cnt = len(re.findall(BANNED, TEXT, re.I))
gate("G7 banned word absent", cnt == 0, f"{cnt} occurrences")

# ---------------------------------------------------------------- G8 draft line
# 2026-09-01: later drafts (0.9-0.11) now head the line; history presence check.
m = re.search(r"Draft 0\.8 \(2026-08-12\)", TEXT)
gate("G8 draft history updated to 0.8 / 2026-08-12", m is not None,
     m.group(0) if m else "MISSING")

# ---------------------------------------------------------------- report
print("ROUND 42 adoption gates -- papers/paper3_mechanism.md")
print("=" * 68)
for name, ok, detail in results:
    print(f"{'PASS' if ok else 'FAIL'}  {name}  [{detail}]")
print("=" * 68)
print(f"abstract: {n_abs} words, {len(num_tokens)} numeric quantities")
print(f"longest sentence: {longest[0]} words")
print(f"longest paragraph anywhere: {worst[0]} words (section {worst[1]})")
n_fail = sum(1 for _, ok, _ in results if not ok)
print(f"{len(results) - n_fail}/{len(results)} gates PASS")
sys.exit(1 if n_fail else 0)
