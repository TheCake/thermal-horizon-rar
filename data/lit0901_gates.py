# -*- coding: utf-8 -*-
"""lit0901_gates.py -- gates for the 2026-09-01 pre-send delta-scan engagement.

  G1  P1 cites Boufourou 2608.24556 with the Sec-9.4 engagement          -> vetoes the P1 engagement
  G2  P1 engagement carries the NON-CONVERSION sentence (moment rule)    -> vetoes the estimator-convention discipline
  G3  P1 engagement states the error-inflation disagreement + DR4 meet   -> vetoes the collision statement
  G4  P2 rival family = FOUR incl. Escala 2608.10073, sign stated        -> vetoes the P2 engagement
  G5  P3 cites Rostami 2511.05632 with the no-occupation differentiation -> vetoes the P3 engagement
  G6  draft lines advanced: P1 0.9 / P2 0.9 / P3 0.11, dated 2026-09-01  -> vetoes the version record
  G7  reference entries present in all three papers                      -> vetoes the citation booking
  G8  P2 abstract unchanged (<= 252 words; no engagement leak)           -> vetoes the register contract
"""
import io, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
P1 = io.open(os.path.join(ROOT, "papers", "paper1_wide_binaries.md"), encoding="utf-8").read()
P2 = io.open(os.path.join(ROOT, "papers", "paper2_rar_coefficients.md"), encoding="utf-8").read()
P3 = io.open(os.path.join(ROOT, "papers", "paper3_mechanism.md"), encoding="utf-8").read()

fails = []
def gate(name, cond, msg=""):
    print(f"{name}: {'PASS' if cond else 'FAIL'} {msg}")
    if not cond: fails.append(name)

gate("G1", "Boufourou (2026)" in P1 and "arXiv:2608.24556" in P1
     and "recovers boost factors of 1.08–1.13" in P1
     and "rejecting the claimed 1.4 at about 16σ" in P1)
gate("G2", "We do not convert his fitted boost parameter into our median-ratio convention" in P1
     and "differently weighted moments" in P1)
gate("G3", "applies no inflation to the formal Gaia errors" in P1
     and "his time-stamped DR4 protocol and ours" in P1)
gate("G4", "at least four distinguishable cosmological predictions" in P2
     and "Escala" in P2 and "arXiv:2608.10073" in P2
     and "rising with redshift like the lock, from a different normalization" in P2)
gate("G5", "Rostami, Rezazadeh & Rostampour" in P3 and "arXiv:2511.05632" in P3
     and "no occupation number, no derived horizon lock" in P3)
gate("G6", "*Draft 0.9 (2026-09-01)" in P1 and "*Draft 0.9 (2026-09-01)" in P2
     and "*Draft 0.11 (2026-09-01)" in P3)
gate("G7", "- Boufourou, H. 2026, arXiv:2608.24556" in P1
     and "- Escala, A. 2026, arXiv:2608.10073" in P2
     and "- Rostami, A., Rezazadeh, K., & Rostampour, A. 2026" in P3)
i, j = P2.find("## Abstract"), P2.find("## 1.")
nw = len(P2[i:j].replace("## Abstract", "").split())
gate("G8", nw <= 252 and "Boufourou" not in P2[i:j] and "Escala" not in P2[i:j], f"({nw} words)")

print()
print(f"LIT0901 GATES: {8-len(fails)}/8 PASS" + (f"  FAILED: {fails}" if fails else ""))
sys.exit(1 if fails else 0)
