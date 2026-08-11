# -*- coding: utf-8 -*-
"""ROUND 37 verification addendum (P3 referee round 3, draft 0.6).

Protocol (the 87a4676 standard, paper-round form as in round32_addendum):
  GA half = BLIND, written and committed BEFORE the referee report exists.
    Verifies every load-bearing number the draft-0.6 absorb ADDED to
    papers/paper3_mechanism.md against the archived stage outputs and
    program record, plus register/repro checks.
  GB half = appended AFTER the report: independently re-compute every
    load-bearing number the referee's ruling rests on (memory rule
    feedback-verify-reviewer-math).

Run: py calcs/round37_addendum.py  (writes data/round37_addendum.txt)
"""
import io
import os
import re

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def rd(path):
    with io.open(os.path.join(ROOT, path), encoding="utf-8") as f:
        return f.read()


OUT = []


def log(s):
    OUT.append(s)
    print(s)


PAPER = rd("papers/paper3_mechanism.md")
S10M = rd("data/stage10m_hartree.txt")
S10N = rd("data/stage10n_crossing.txt")
S10Q = rd("data/stage10q_o5amplitude.txt")
R33 = rd("data/round33_addendum.txt")
PRED = rd("PREDICTIONS.md")
LEDG = rd("LEDGER.csv")

npass = nfail = 0


def check(name, ok, detail=""):
    global npass, nfail
    tag = "PASS" if ok else "FAIL"
    if ok:
        npass += 1
    else:
        nfail += 1
    log("GA %-58s %s  %s" % (name, tag, detail))


log("=" * 72)
log("ROUND 37 ADDENDUM -- GA half (blind, pre-report)")
log("=" * 72)

# GA-1 abstract word count (STYLE bar <= 250)
m = re.search(r"## Abstract\n\n(.+?)\n\n## 1\.", PAPER, re.S)
words = len(m.group(1).split())
check("1 abstract word count <= 250", words <= 250, "n=%d" % words)

# GA-2 token-vs-archive greps (new text only)
tok = [
    ("0.234 central e_a", "0.234", S10Q),
    ("0.159 envelope lo", "0.159", S10Q),
    ("0.388 envelope hi", "0.388", S10Q),
    ("coth spread 0.089 (quoted 9 percent)", "spread 0.089", S10Q),
    ("0.268 matched co-quote (program record)", "0.268", PRED),
    ("0.0111 hi-envelope S_max", "0.0111", S10Q),
    ("0.0056 central S_max", "0.0056", S10Q),
    ("0.00389 frozen-draw avg (quoted 0.0039)", "0.00389", S10Q),
    ("0.662 revival amplitude (program record)", "0.662", PRED),
    ("1.43 FD ceiling", "1.43", S10Q),
    ("refresh 83 Gyr", "83", S10Q),
    ("refresh 282 Gyr", "282", S10Q),
    ("x33 over wedge", "33", S10Q),
    ("x29 short of closure", "29", S10Q),
    ("2.48 Gyr wedge boundary (ledger)", "2.48", LEDG),
    ("10M max edge 0.700", "0.700", S10M),
    ("10M lock floor 0.888", "0.888", S10M),
    ("10M edge ratio 4.03", "4.03", S10M),
    ("10M edge ratio 4.27", "4.27", S10M),
    ("10M naive ratio -> 2", "ratio (n -> oo):        2", S10M),
    ("10M midpoint theorem exact", "om_MF(n+1/2) == D(n):        residue = 0", S10M),
    ("10M every-moment conservation", "every moment of N_s is\n   conserved", S10M),
    ("10N deep-gal x20 suppression", "suppressed x20", S10N),
    ("10N peak kappa* = 0.8918", "0.8918", S10N),
    ("R33 plateau 1.088 (quoted 1.09)", "1.088", R33),
    ("anharmonic 8-15% (program record)", "8–15%", PRED),
    ("virial factor 2 disclosure (program record)", "virial factor of 2", PRED),
]
for name, needle, hay in tok:
    check("2 " + name, needle in hay)

# paper-side: the quoted forms actually appear in the draft
# (pre-report blind form preserved at commit 3300f88 checked "factor of 29",
#  faithful to the then-draft; the R37-adopted draft quotes the corrected
#  cross-convention chain -- tokens updated post-report, disclosed here)
quoted = ["0.23 ", "0.16 to 0.39", "a factor of about 40", "6.1 times below",
          "83 to 282 Gyr",
          "factor of 33", "about 2.5 Gyr", "0.70 against a lock-band floor of 0.89",
          "4.0 to 4.3", "ratio of 1.09", "0.011 where a survey precision of 0.02",
          "0.0056 at the central amplitude", "0.0039 averaged over frozen draws",
          "1.7 times the high envelope", "16 to 39 percent", "11 to 27 percent",
          "8 to 15 percent", "9 percent correction", "0.27 under a matched-treatment",
          "about 0.002", "ratio 1.84", "in the frozen-bath treatment",
          "neither self-consistency route bounds"]
for q in quoted:
    check("2p in-draft: '%s'" % q[:40], q in PAPER)

# GA-3 arithmetic re-derivations of the rounded claims
# (the blind form at 3300f88 carried "factor-29 = 6.7/0.234" -- a FAITHFUL
#  check of the then-draft's cross-convention slip; R37 M1 corrected the
#  object: 6.7 is FD-relative, the absolute closure surface is q_phys = 9.4)
ar = [
    ("factor-40 = 9.45/0.2338 (max-grant, same units)", 9.45 / 0.2338, 40.4, 0.5),
    ("factor-40 alt = 6.67*qFD/0.2338", 6.67 * 1.4156 / 0.2338, 40.4, 0.5),
    ("6.1-below-FD = 1.4156/0.2338", 1.4156 / 0.2338, 6.05, 0.05),
    ("factor-33 lo = 83/2.48", 83 / 2.48, 33.5, 0.2),
    ("revival x1.7 = 0.662/0.388", 0.662 / 0.388, 1.706, 0.01),
    ("FD-frac lo = 0.159/1.43", 0.159 / 1.43, 0.111, 0.003),
    ("FD-frac hi = 0.388/1.43", 0.388 / 1.43, 0.271, 0.003),
    ("x20 = (4/27)/0.00751/ (1/0.051)", (4.0 / 27.0) / 0.00751, 19.7, 0.5),
    ("plateau 1.09 quoted vs archived 1.088", 1.088, 1.09, 0.005),
]
for name, got, want, tol in ar:
    check("3 " + name, abs(got - want) <= tol, "got %.4g want %.4g" % (got, want))

# GA-4 Appendix B script/output existence (new rows)
for p in ["calcs/stage10m_hartree.py", "calcs/stage10n_crossing.py",
          "calcs/stage10q_o5amplitude.py", "calcs/stage10o_p15erode.py",
          "calcs/round33_addendum.py", "calcs/round34_addendum.py",
          "calcs/round36_addendum.py", "data/stage10m_hartree.txt",
          "data/stage10n_crossing.txt", "data/stage10q_o5amplitude.txt",
          "data/round33_addendum.txt", "data/stage10o_run6.log"]:
    check("4 exists: " + p, os.path.exists(os.path.join(ROOT, p)))

# GA-5 register checks
body = PAPER.split("## Appendix B")[0]
check("5 no stage codenames in body", re.search(r"[Ss]tage\s*10[a-qA-Q]\b", body) is None)
check("5 banned word absent", ("anti" + "gravity") not in PAPER.lower())
check("5 em-dash count small", PAPER.count("—") <= 3, "n=%d" % PAPER.count("—"))
check("5 draft line bumped", "Draft 0.7 (2026-08-11)" in PAPER)

# GA-6 sentence stats of the new paragraphs (bar: none >= 70 words)
new_txt = []
for anchor in ["Since the exclusion, the one environment-side unknown",
               "The named remaining route, a joint coherent-state",
               "The ambient-amplitude test, stated with its grade",
               "The degeneracy's response signature has since been quantified"]:
    i = PAPER.find(anchor)
    check("6 anchor present: %s..." % anchor[:34], i >= 0)
    if i >= 0:
        new_txt.append(PAPER[i:PAPER.find("\n", i)])
sent = []
for t in new_txt:
    sent += [s for s in re.split(r"(?<=[.!?])\s+", t) if len(s.split()) > 2]
mx = max(len(s.split()) for s in sent)
mean = sum(len(s.split()) for s in sent) / float(len(sent))
check("6 max sentence < 70 words", mx < 70, "max=%d mean=%.1f n=%d" % (mx, mean, len(sent)))

log("-" * 72)
log("GA SUMMARY: %d PASS / %d FAIL" % (npass, nfail))
log("=" * 72)

# ---------------------------------------------------------------------------
# GB half (post-report): independently re-compute every load-bearing number
# in the ROUND 37 report before adoption (memory rule).
# ---------------------------------------------------------------------------
log("")
log("=" * 72)
log("ROUND 37 ADDENDUM -- GB half (post-report reviewer-math verification)")
log("=" * 72)

S10G = rd("data/stage10g_collective.txt")

gb_pass = gb_fail = 0


def gb(name, ok, detail=""):
    global gb_pass, gb_fail
    tag = "OK" if ok else "MISMATCH"
    if ok:
        gb_pass += 1
    else:
        gb_fail += 1
    log("GB %-58s %s  %s" % (name, tag, detail))


# GB-1  M1: the closure-factor units chain, from the 10G archive itself
gb("1a 10G table ran at q_phys = 1", "q_phys = 1; requirement 0.5 H" in S10G)
gb("1b 10G closure surface q_phys >= 9.4",
   "closure surface needs q_phys >= 9.4" in S10G)
gb("1c max-grant shortfall 9.45x archived", "9.45x" in S10G)
qfd = (2 * 0.502 + 1) ** 0.5
gb("1d q_FD = sqrt(2*0.502+1) = 1.416", abs(qfd - 1.4156) < 2e-3,
   "%.4f" % qfd)
gb("1e FD-relative 6.67 = 9.45/q_FD", abs(9.45 / qfd - 6.67) < 0.03,
   "%.3f" % (9.45 / qfd))
gb("1f consistent factor 9.45/0.2338 = 40.4",
   abs(9.45 / 0.2338 - 40.4) < 0.3, "%.1f" % (9.45 / 0.2338))
gb("1g his alt route 6.7*(qFD/0.234) = 40.6",
   abs(6.7 * (qfd / 0.234) - 40.6) < 0.4, "%.1f" % (6.7 * (qfd / 0.234)))
gb("1h variance cross-check 6.67^2 = 44.5 (~ the 45 in Sec 7)",
   abs(6.67 ** 2 - 44.5) < 0.3, "%.1f" % 6.67 ** 2)
gb("1i variance cross-check 9.18^2 = 84.3 (~ the 84 in Sec 7)",
   abs(9.18 ** 2 - 84.3) < 0.4, "%.1f" % 9.18 ** 2)
gb("1j the retired 29 = FD-relative/absolute mix 6.7/0.234",
   abs(6.7 / 0.234 - 28.6) < 0.2, "%.1f (units mismatch, conservative)" % (6.7 / 0.234))
gb("1k his honest-surface 118 = 27.6/0.234",
   abs(27.6 / 0.234 - 118.0) < 1.0, "%.0f" % (27.6 / 0.234))

# GB-2  m1: convention-family structure of the fine-grade straddle
gb("2a conv=1.00 central S = 0.00218 archived",
   "conv1.00: Sc=0.00218" in S10Q)
gb("2b S_max(central) 0.00558 = the conv-sqrt2 cell",
   "conv1.41: Sc=0.00558" in S10Q and "S_max(central) = 0.00558" in S10Q)
gb("2c paper quotes ~0.002 nominal", "about 0.002" in PAPER)

# GB-3  m2: static-ratio corner
gb("3a corner ratio 0.159/0.086 = 1.85", abs(0.159 / 0.086 - 1.849) < 0.01,
   "%.3f" % (0.159 / 0.086))
gb("3b central ratio 0.2338/0.086 = 2.72", abs(0.2338 / 0.086 - 2.72) < 0.02,
   "%.2f" % (0.2338 / 0.086))
gb("3c stage flagged the corner", "1.84" in S10Q)

# GB-4  m3: 0.005 provenance = registered target, not demonstrated feasibility
gb("4a 0.005 grade in program record", "0.005" in PRED)
gb("4b paper softened to 'of order 0.005'", "of order 0.005" in PAPER)

# GB-5  his reproduction-table spot values
gb("5a revival 0.6615 -> x1.70", abs(0.6615 / 0.3883 - 1.704) < 0.01,
   "%.3f" % (0.6615 / 0.3883))
gb("5b plateau 1.088", "1.088" in R33)
gb("5c r-sigma note: (0.5-0.34)/0.19 = 0.84",
   abs((0.5 - 0.34) / 0.19 - 0.842) < 0.01, "%.3f" % ((0.5 - 0.34) / 0.19))
gb("5d r-sigma unrounded (0.5-0.3365)/0.1869 = 0.875",
   abs((0.5 - 0.3365) / 0.1869 - 0.875) < 0.01,
   "%.3f" % ((0.5 - 0.3365) / 0.1869))

# GB-6  abstract still within budget after adoption
m2 = re.search(r"## Abstract\n\n(.+?)\n\n## 1\.", PAPER, re.S)
w2 = len(m2.group(1).split())
gb("6 abstract <= 250 post-adoption", w2 <= 250, "n=%d" % w2)

# GB-7  round-36 output now archived (his item iv.3)
gb("7 data/round36_addendum.txt archived",
   os.path.exists(os.path.join(ROOT, "data/round36_addendum.txt")))

log("-" * 72)
log("GB SUMMARY: %d OK / %d MISMATCH" % (gb_pass, gb_fail))
verdict = ("ALL REVIEWER NUMBERS CONFIRMED; M1 units chain independently "
           "re-derived from the 10G archive (closure q_phys = 9.4 absolute; "
           "6.7 FD-relative; consistent factor ~40; the 29 mixed conventions, "
           "conservative direction)" if gb_fail == 0 else
           "MISMATCHES PRESENT -- do not adopt without resolution")
log("GB VERDICT: " + verdict)
log("=" * 72)

with io.open(os.path.join(ROOT, "data/round37_addendum.txt"), "w",
             encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
