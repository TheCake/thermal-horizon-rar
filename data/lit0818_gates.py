"""
lit0818_gates.py -- gates for the 2026-08-18 a0(z) engagement patch.

Checks:
  G1  "2608.03576" present in paper2, paper3, and the PREDICTIONS annotation
  G2  "2504.20857" present in paper2 (at least)
  G3  every numeric slope / sigma value quoted in the NEW paper text appears
      verbatim in data/lit0818_a0z.txt
  G4  both draft lines advanced (paper2 -> 0.7, paper3 -> 0.10, both dated
      2026-08-18 with the engagement note)
  G5  abstracts byte-identical to their git HEAD versions
  G6  no sentence of 70 words or more in the new text
  G7  banned word (checked as a concatenation, never spelled in this file)
      occurs zero times in the touched files

Run:  py data/lit0818_gates.py > data/lit0818_gates.txt
"""

import io
import os
import re
import subprocess
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P2 = os.path.join(ROOT, "papers", "paper2_rar_coefficients.md")
P3 = os.path.join(ROOT, "papers", "paper3_mechanism.md")
PRED = os.path.join(ROOT, "PREDICTIONS.md")
CALC = os.path.join(ROOT, "data", "lit0818_a0z.txt")

results = []


def gate(name, ok, detail=""):
    results.append((name, bool(ok), detail))
    print("[%s] %s%s" % ("PASS" if ok else "FAIL", name,
                         ("  -- " + detail) if detail else ""))


def read(path):
    with io.open(path, encoding="utf-8") as fh:
        return fh.read()


def git_show(relpath):
    out = subprocess.run(["git", "show", "HEAD:" + relpath],
                         cwd=ROOT, capture_output=True)
    if out.returncode != 0:
        return None
    return out.stdout.decode("utf-8")


p2 = read(P2)
p3 = read(P3)
pred = read(PRED)
calc = read(CALC)

print("lit0818_gates.py -- gates for the a0(z) engagement patch (2026-08-18)")
print("-" * 78)

# ------------------------------------------------------------------ G1 / G2
gate("G1a 2608.03576 in paper2", "2608.03576" in p2)
gate("G1b 2608.03576 in paper3", "2608.03576" in p3)
ann_start = pred.find("[2026-08-18 annotation")
ann_end = pred.find("### P4", ann_start) if ann_start >= 0 else -1
ann = pred[ann_start:ann_end] if ann_start >= 0 else ""
gate("G1c 2608.03576 in the PREDICTIONS 2026-08-18 annotation",
     "2608.03576" in ann, "annotation length %d chars" % len(ann))
gate("G2 2504.20857 in paper2", "2504.20857" in p2)
gate("G2b 2504.20857 in paper3 (not required)", True,
     "present: %s" % ("2504.20857" in p3))

# --------------------------------------------------------------------- G3
# Every numeric quantity the new paper text asserts about the comparison.
# Each must appear as a substring of the computation's output file.
QUOTED = [
    "0.52",    # a1_lock, secant over their window
    "0.49",    # da0/dz at z = 0
    "1.60",    # their pure-sample a1 central
    "2.33",    # their pure-sample a1 error
    "0.69",    # sigma distance of zero from their band
    "0.91",    # sigma distance of the lock slope from their band
    "0.22",    # separation of the two hypotheses in units of their sigma
    "5.23",    # their anchored a1
    "1.05",    # their anchored a1 error
    "4.5",     # sigma distance of their anchored a1 above the lock
    "0.47",    # anchored slope x window
    "0.39",    # their own intercept gap
    "1.54",    # their sample-alone intercept
    "1.15",    # their anchored intercept
    "4.33",    # slope reproducing the intercept gap alone
    "0.85",    # sigma distance of their anchored a1 from the pure-offset slope
    "21",      # residual as a percentage of the gap
]
missing = [q for q in QUOTED if q not in calc]
gate("G3 all quoted numerics present verbatim in lit0818_a0z.txt",
     not missing, "missing: %s" % (missing if missing else "none"))

# Reverse guard: no stray numeric of the comparison form in the new paragraphs
# that is absent from the computation output.
NEW_BLOCKS = []
for text, marker_start, marker_end in [
    (p2, "A first direct fit of a₀(z) has since appeared at the other end",
     "The binary translation's high-side pull"),
    # 2026-09-01: end marker updated -- the rival-family sentence now reads
    # "at least four" (Escala addition, draft 0.9); scope of this block is
    # unchanged (it ends where the rival-family paragraph begins).
    (p2, "The low-redshift fit of Vărăşteanu et al. (2026) prefers",
     "Stated plainly, at least four"),
    (p2, "A first direct slope now exists as well", "\n"),
    (p3, "A first direct low-redshift fit of the same quantity", "\n"),
]:
    i = text.find(marker_start)
    j = text.find(marker_end, i) if i >= 0 else -1
    if i >= 0 and j > i:
        NEW_BLOCKS.append(text[i:j])
    elif i >= 0:
        NEW_BLOCKS.append(text[i:i + 3000])

gate("G3b new-text blocks located", len(NEW_BLOCKS) == 4,
     "found %d of 4" % len(NEW_BLOCKS))

ALLOWED_CONTEXT = {
    "130", "0.09", "2026", "2025", "19", "10", "70", "0.3", "5.0", "2.5",
    "2608.03576", "2504.20857", "0", "1", "2", "5", "1.02", "3.4", "2.3",
    "0.7", "0.1", "1.05", "5.1",
}
# tokens above are section numbers, sample sizes, years, the target paper's own
# labels, and figures already published elsewhere in the drafts; every remaining
# numeric in the new text must trace to the computation output.
stray = []
for blk in NEW_BLOCKS:
    for tok in re.findall(r"\d+\.\d+|\d+", blk):
        if tok in ALLOWED_CONTEXT:
            continue
        if tok not in calc:
            stray.append(tok)
gate("G3c no comparison numeric in new text absent from the computation output",
     not stray, "stray: %s" % (sorted(set(stray)) if stray else "none"))

# --------------------------------------------------------------------- G4
# 2026-08-24: draft 0.8 (SOL round) now heads the line; the 0.7 record must
# remain in the header history -- check presence, not line-leading position.
gate("G4a paper2 draft advanced to 0.7 dated 2026-08-18",
     "Draft 0.7 (2026-08-18)" in p2)
gate("G4b paper2 draft line carries the engagement note",
     "a₀(z) measurement engaged (MIGHTEE-HI/LADUMA)" in p2.split("\n")[4])
# 2026-09-01: draft 0.11 now heads the line; presence check like G4a.
gate("G4c paper3 draft advanced to 0.10 dated 2026-08-18",
     "Draft 0.10 (2026-08-18)" in p3)
gate("G4d paper3 draft line carries the engagement note",
     "a₀(z) measurement engaged (MIGHTEE-HI/LADUMA)" in p3.split("\n")[4])

# --------------------------------------------------------------------- G5


def abstract_of(text):
    i = text.find("## Abstract")
    j = text.find("\n## ", i + 5)
    return text[i:j]


for label, cur, rel in [("paper2", p2, "papers/paper2_rar_coefficients.md"),
                        ("paper3", p3, "papers/paper3_mechanism.md")]:
    head = git_show(rel)
    if head is None:
        gate("G5 %s abstract vs HEAD" % label, False, "git show failed")
        continue
    a_new = abstract_of(cur)
    a_old = abstract_of(head)
    gate("G5 %s abstract byte-identical to HEAD" % label, a_new == a_old,
         "%d chars" % len(a_new))

# --------------------------------------------------------------------- G6
longest = 0
longest_txt = ""
for blk in NEW_BLOCKS + [ann]:
    flat = re.sub(r"\s+", " ", blk)
    for sent in re.split(r"(?<=[.;:])\s+(?=[A-ZΩ(])", flat):
        n = len(sent.split())
        if n > longest:
            longest, longest_txt = n, sent
gate("G6 no sentence of 70 words or more in new text", longest < 70,
     "longest = %d words: %s" % (longest, longest_txt[:110]))

# --------------------------------------------------------------------- G7
banned = ("anti" + "gravity")
counts = {}
for label, text in [("paper2", p2), ("paper3", p3), ("PREDICTIONS", pred),
                    ("lit0818_a0z.txt", calc)]:
    counts[label] = text.lower().count(banned)
gate("G7 banned word zero occurrences in touched files",
     sum(counts.values()) == 0, str(counts))

print("-" * 78)
n_fail = sum(1 for _, ok, _ in results if not ok)
print("GATES: %d/%d PASS%s" % (len(results) - n_fail, len(results),
                               "" if n_fail == 0 else "  (%d FAIL)" % n_fail))
print("END lit0818_gates.py")
sys.exit(1 if n_fail else 0)
