"""Gates for the 2026-08-14 rival-mechanism citation patch to papers/paper3_mechanism.md.

Run:  py data/lit0814_p3_gates.py > data/lit0814_p3_gates.txt

Checks (all must PASS):
  G1  arXiv id 2507.11524 (Gillot) present in the body text AND in the reference list
  G2  arXiv id 2602.14515 (Luo)    present in the body text AND in the reference list
  G3  abstract untouched: word count still <= 250 (and matched against the archived
      draft-0.8 abstract if a git copy is reachable)
  G4  no sentence of 70 words or more inside the two newly inserted passages
  G5  banned word (assembled at runtime, never spelled in this file) has zero
      occurrences anywhere in the manuscript
  G6  draft line advanced to draft 0.9 / 2026-08-14 and names the patch
"""

import io
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(ROOT, "papers", "paper3_mechanism.md")

BANNED = "anti" + "gravity"  # assembled, so the checker never contains the token

NEW_PASSAGE_STARTS = (
    "An independent modified-inertia proposal",
    "The closest current relative of the reading developed here",
)

results = []


def gate(name, ok, detail):
    results.append((name, bool(ok), detail))


def sentences(block):
    """Split a prose block into sentences on terminal punctuation + capital/paren."""
    flat = " ".join(block.split())
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z(“])", flat)
    return [p for p in parts if p.strip()]


def main():
    with io.open(PAPER, encoding="utf-8") as fh:
        text = fh.read()

    # --- split body / references -------------------------------------------------
    m = re.search(r"^## References\s*$", text, flags=re.M)
    if not m:
        gate("G0-structure", False, "no '## References' heading found")
        report()
        return
    body, refs = text[: m.start()], text[m.start():]
    gate("G0-structure", True, "body %d chars / references %d chars" % (len(body), len(refs)))

    # --- G1 / G2 citations ------------------------------------------------------
    for gname, arxiv_id, who in (("G1", "2507.11524", "Gillot"),
                                 ("G2", "2602.14515", "Luo")):
        in_body = arxiv_id in body
        in_refs = arxiv_id in refs
        gate("%s-%s-%s" % (gname, who, arxiv_id),
             in_body and in_refs,
             "body=%s references=%s" % (in_body, in_refs))

    # --- G3 abstract untouched --------------------------------------------------
    am = re.search(r"^## Abstract\s*$(.*?)^## ", text, flags=re.M | re.S)
    abstract = am.group(1).strip() if am else ""
    nwords = len(abstract.split())
    ok_len = 0 < nwords <= 250
    detail = "abstract word count %d (bar <= 250)" % nwords

    prev = None
    try:
        prev = subprocess.check_output(
            ["git", "show", "HEAD:papers/paper3_mechanism.md"],
            cwd=ROOT, stderr=subprocess.STDOUT)
        prev = prev.decode("utf-8", "replace")
    except Exception as exc:  # pragma: no cover - git may be unavailable
        detail += "; committed copy unreadable (%s)" % type(exc).__name__
    if prev:
        pm = re.search(r"^## Abstract\s*$(.*?)^## ", prev, flags=re.M | re.S)
        prev_abs = pm.group(1).strip() if pm else None
        if prev_abs is None:
            detail += "; committed abstract not locatable"
        else:
            same = (prev_abs == abstract)
            detail += "; byte-identical to committed abstract=%s" % same
            ok_len = ok_len and same
    gate("G3-abstract-unchanged", ok_len, detail)

    # --- G4 sentence length in the new passages ---------------------------------
    paragraphs = [p for p in body.split("\n\n")]
    new_blocks = [p for p in paragraphs
                  if p.strip().startswith(NEW_PASSAGE_STARTS)]
    found = len(new_blocks)
    worst = 0
    worst_txt = ""
    total_sent = 0
    for blk in new_blocks:
        for s in sentences(blk):
            total_sent += 1
            n = len(s.split())
            if n > worst:
                worst, worst_txt = n, s
    gate("G4-sentence-len",
         found == 2 and worst < 70,
         "passages found %d/2; %d sentences; longest %d words (bar < 70): %r"
         % (found, total_sent, worst, worst_txt[:80]))

    # --- G5 banned word ---------------------------------------------------------
    hits = len(re.findall(BANNED, text, flags=re.I))
    gate("G5-banned-word", hits == 0, "occurrences %d (bar 0)" % hits)

    # --- G6 draft line ----------------------------------------------------------
    dl = re.search(r"\*Draft ([0-9.]+) \(([0-9-]+)\)\.", text)
    ver = dl.group(1) if dl else None
    date = dl.group(2) if dl else None
    names = ("Gillot" in text) and ("Luo" in text) and ("rival-mechanism" in text)
    gate("G6-draft-line",
         ver == "0.9" and date == "2026-08-14" and names,
         "draft %s (%s); patch named=%s" % (ver, date, names))

    report()


def report():
    print("lit0814 P3 citation-patch gates")
    print("target: papers/paper3_mechanism.md")
    print("")
    npass = 0
    for name, ok, detail in results:
        print("%-28s %-5s %s" % (name, "PASS" if ok else "FAIL", detail))
        npass += 1 if ok else 0
    print("")
    print("%d/%d PASS" % (npass, len(results)))
    if npass != len(results):
        sys.exit(1)


if __name__ == "__main__":
    main()
