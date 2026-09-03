# -*- coding: utf-8 -*-
"""SOL-round gates: verify the 2026-08-24 external-review adoption (papers r0.8).

Checks (each names what a failure vetoes):
  G-SOL-1  P1 no longer contains the withdrawn field-wide attribution      -> vetoes the Sec 5.3 rewrite claim
  G-SOL-2  P1 Sec 5.3 engagement numbers match data/reviewsol_wobble.txt   -> vetoes the population comparison
  G-SOL-3  P1 carries the Penoyre reference + App A withdrawal bullet      -> vetoes the correction booking
  G-SOL-4  P1 'effective' scoping present at abstract + conclusion 1       -> vetoes the 2.1x rescope
  G-SOL-5  P2 two-grade zero-exclusion wording at all four sites           -> vetoes the c1 fix
  G-SOL-6  P2 horizon-choice scope statement present in falsifier 2        -> vetoes the FRW clarification
  G-SOL-7  correction count 'twenty' consistent (P1 App A, P2 App A+intro) -> vetoes the count sync
  G-SOL-8  both papers carry the bar-locking != preregistration sentence   -> vetoes the methods scope fix
  G-SOL-9  P2 abstract <= 252 words                                        -> vetoes the register contract
"""
import io, re, sys

P1 = io.open(r"papers/paper1_wide_binaries.md" if len(sys.argv) < 2 else sys.argv[1], encoding="utf-8").read() \
     if False else io.open(__file__.replace("data\\reviewsol_gates.py", "papers\\paper1_wide_binaries.md").replace("data/reviewsol_gates.py", "papers/paper1_wide_binaries.md"), encoding="utf-8").read()
P2 = io.open(__file__.replace("data\\reviewsol_gates.py", "papers\\paper2_rar_coefficients.md").replace("data/reviewsol_gates.py", "papers/paper2_rar_coefficients.md"), encoding="utf-8").read()
W  = io.open(__file__.replace("reviewsol_gates.py", "reviewsol_wobble.txt"), encoding="utf-8").read()

fails = []
def gate(name, cond, msg=""):
    print(f"{name}: {'PASS' if cond else 'FAIL'} {msg}")
    if not cond: fails.append(name)

# G-SOL-1
gate("G-SOL-1", "used in the published wide-binary analyses scales as q/(1+q)" not in P1
     and "law in field use omits photocenter cancellation" not in P1)

# G-SOL-2: engagement numbers 74% / 49-63% / 1.3x must be supported by the pinned output
w_meas = re.search(r"measured-twin-t5 survival\s+= ([0-9.]+)\s+\(removes ([0-9.]+)%\)", W)
w_ps   = re.search(r"P&S-population survival R_var\s+= ([0-9.]+)\s+\(cancellation removes ([0-9.]+)%", W)
w_bk   = re.search(r"Banik\+24-population survival\s+= ([0-9.]+)\s+\(removes ([0-9.]+)%\)", W)
ok2 = (w_meas and abs(float(w_meas.group(2)) - 73.8) < 0.3
       and w_ps and abs(float(w_ps.group(2)) - 48.8) < 0.3
       and w_bk and abs(float(w_bk.group(2)) - 63.4) < 0.3
       and "removes 74% of the naive per-companion wobble variance" in P1
       and "49–63% for the published population inputs" in P1
       and "about 1.3 times the per-companion wobble variance" in P1)
gate("G-SOL-2", ok2, f"(txt: {w_meas.group(2) if w_meas else '?'} / {w_ps.group(2) if w_ps else '?'}-{w_bk.group(2) if w_bk else '?'})")

# G-SOL-3
gate("G-SOL-3", "Penoyre, Z., Belokurov, V., & Evans, N. W. 2022, MNRAS 513, 2437" in P1
     and "withdrawn after a primary-source audit" in P1)

# G-SOL-4 baseline superseded 2026-09-03 (correction #21): the 2.1x claim
# these scoping strings guarded is WITHDRAWN; the gate now pins the
# withdrawal at the abstract + conclusion sites and the surviving
# joint-fit scoping in Sec 4.1.
# needle updated 2026-09-03 (R43 adoption pass): abstract wording tightened
# to "failed its null injection" for the 250-word contract.
gate("G-SOL-4", "failed its null injection" in P1
     and "the narrow-bin calibrator that earlier versions quoted" in P1
     and "effective inflation the joint forward model requires" in P1
     and "the effective pair-level velocity errors are at least 2.1 times" not in P1)

# G-SOL-5: all four P2 sites
gate("G-SOL-5", "zero excluded at profile grade; bootstrap positivity 89–96%" in P2      # abstract
     and "excluded at profile grade (see text)" in P2                                    # table 2
     and "at the conservative grade zero is disfavored, not excluded" in P2              # prose
     and "Zero is excluded at profile grade in both treatments (bootstrap positivity 95.5% and 89.0%)" in P2  # fig 2
     and "excluded in every treatment" not in P2
     and "excluded in each row" not in P2
     and "zero excluded everywhere" not in P2)

# G-SOL-6
gate("G-SOL-6", "requires a choice of horizon, and we use the instantaneous Hubble rate" in P2
     and "not a theorem of the original construction" in P2)

# G-SOL-7 baseline advanced 2026-09-03 (correction #21): count = twenty-one.
gate("G-SOL-7", "twenty-one are on record" in P1
     and "The program logged twenty-one corrections" in P2
     and "Twenty-one corrections logged during the program" in P2
     and "nineteen" not in P1.lower() and "nineteen" not in P2.lower())

# G-SOL-8
gate("G-SOL-8", "not conventional preregistration" in P1 and "not conventional preregistration" in P2)

# G-SOL-9
i, j = P2.find("## Abstract"), P2.find("## 1.")
nw = len(P2[i:j].replace("## Abstract", "").split())
gate("G-SOL-9", nw <= 252, f"({nw} words)")

print()
print(f"SOL GATES: {9-len(fails)}/9 PASS" + (f"  FAILED: {fails}" if fails else ""))
sys.exit(1 if fails else 0)
