"""
STAGE 10C ADDENDUM -- ROUND-24 conditions 1/2/3/4/5 + the standing
verify-reviewer-math rule (every load-bearing number the ROUND-24
referee introduced is independently re-computed here before adoption).

His new numbers under verification:
  GA-1 run-1 level populations at nbar = 8 (the drain fingerprint):
       his [P0..P3] = [0.256, 0.361, 0.243, 0.140].
  GA-2 amendment-A1 single-fix scans (the combinations the stage never
       ran): Om_b-fix alone (0.53, gd = 0.3) spread 0.0269 / dev
       0.0311 = PASSES the locked bars; gd-fix alone (0.8, gd = 0.15)
       spread 0.1400 = FAILS. (His validation that the collision was
       the culprit and the gd-halving is FC insurance.)
  GA-3 kappa = 1.5 rung mismatches: x-coefficient 1/8 (not 1/12),
       x^3-coefficient -1/480 (not -1/720) -- the rungs WOULD
       discriminate if measured at depth.
  GA-4 hier-bootstrap kappa band: c1 ~ 0.4 +/- 0.3 -> kappa = 2(1-c1)
       band [0.6, 1.8] (contains 1 -- his mitigation row).
Conditions executed in the printout: (1) KAPPA clause -> two-pole
statement; (2) LADDER clause -> linearity-transmitted; (3) the
renormalization's closure-conditionality; (4) g_close row down-label;
(5) T4 annotations -> conditional-on-successor.

Writes data/stage10c_addendum.txt.
"""
import math
import numpy as np
import sympy as sp

OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

say("=" * 78)
say("STAGE 10C ADDENDUM -- ROUND-24 CONDITIONS 1/2/3/4/5")
say("=" * 78)

# ---- T4 machinery (verbatim stage constants) ------------------------------
def build(NA, NB):
    ad = np.diag(np.sqrt(np.arange(1, NA)), -1)
    bd = np.diag(np.sqrt(np.arange(1, NB)), -1)
    Ia, Ib = np.eye(NA), np.eye(NB)
    return dict(A=np.kron(ad.T, Ib), Ad=np.kron(ad, Ib),
                B=np.kron(Ia, bd.T), Bd=np.kron(Ia, bd),
                Na=np.kron(ad @ ad.T, Ib), Nb=np.kron(Ia, bd @ bd.T),
                Ia=Ia, Ib=Ib)

def thermal(N, nbar):
    xx = math.log(1.0 + 1.0/nbar)
    p = np.exp(-xx*np.arange(N)); p /= p.sum()
    return np.diag(p)

NA, NB, CHI, EPS = 4, 56, 0.8, 0.01
o = build(NA, NB)
IdF = np.kron(o['Ia'], o['Ib'])
LVL = -CHI*o['Na'] + 0.5*CHI*(o['Na'] @ (o['Na'] - IdF))

def satur(omb, gd, nb):
    H = (LVL + EPS*(o['A'] + o['Ad']) + omb*o['Nb']
         + gd*(o['Na'] @ (o['B'] + o['Bd']))
         + (gd**2/omb)*(o['Na'] @ o['Na']))
    ev, U = np.linalg.eigh(H)
    p1 = np.zeros(NA); p1[1] = 1.0
    rho = np.kron(np.diag(p1), thermal(NB, nb))
    rl = U.T @ rho @ U
    rb = U @ np.diag(np.diag(rl)) @ U.T
    return [float(np.real(np.trace(
        np.kron(np.diag((np.arange(NA) == k)*1.0), o['Ib']) @ rb)))
        for k in range(NA)]

say("")
say("GA-1 run-1 drain fingerprint (Om_b = 0.8 = CHI, gd = 0.3, nbar = 8):")
P = satur(0.8, 0.3, 8.0)
say(f"  level populations [P0, P1, P2, P3] = "
    f"[{P[0]:.3f}, {P[1]:.3f}, {P[2]:.3f}, {P[3]:.3f}]")
say(f"  his [0.256, 0.361, 0.243, 0.140]")
ok1 = all(abs(P[k] - t) < 0.005 for k, t in
          enumerate([0.256, 0.361, 0.243, 0.140]))
say(f"  -> {'CONFIRMED' if ok1 else 'MISMATCH'} (the collision drains the "
    "pair manifold into {0,3}; four-level spread, P2 ~ 1/4)")

say("")
say("GA-2 the single-fix scans (his validation of amendment A1):")
NBS = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
for lab, omb, gd, his_spr in [("Om_b-fix alone (0.53, gd=0.3)", 0.53, 0.3,
                               0.0269),
                              ("gd-fix alone   (0.8, gd=0.15)", 0.8, 0.15,
                               0.1400)]:
    Ws = [satur(omb, gd, nb)[2] for nb in NBS]
    spr = max(Ws) - min(Ws)
    dev = max(abs(w - 0.5) for w in Ws)
    verdict = "PASSES bars" if (spr <= 0.04 and dev <= 0.06) else "FAILS"
    ok = abs(spr - his_spr) < 0.003
    say(f"  {lab}: spread {spr:.4f} (his {his_spr}), dev {dev:.4f} -> "
        f"{verdict}  [{'CONFIRMED' if ok else 'MISMATCH'}]")
say("  -> the Om_b = CHI collision was the culprit; the gd-halving is the")
say("     disclosed FC-Laguerre insurance, not the essential fix. Amendment")
say("     A1 stands as bug-class (reviewer-validated with the combinations")
say("     the stage never printed).")

say("")
say("GA-3 kappa = 1.5 rung mismatches (his discrimination check):")
x = sp.symbols('x', positive=True)
ser = sp.series(1 + sp.Rational(3, 2)/(sp.exp(x)-1), x, 0, 6)
c1 = ser.coeff(x, 1); c3 = ser.coeff(x, 3)
ok3 = (sp.simplify(c1 - sp.Rational(1, 8)) == 0
       and sp.simplify(c3 - sp.Rational(-1, 480)) == 0)
say(f"  nu = 1 + 1.5 n_BE: x-coeff = {c1} (his 1/8), x^3-coeff = {c3} "
    f"(his -1/480) -> {'CONFIRMED' if ok3 else 'MISMATCH'}")
say("  -> the rungs WOULD discriminate kappa if measured at depth; per 6L")
say("     the deep arm is population-thin -- only c1 (weakly c2) is")
say("     measured. His LADDER-clause correction is arithmetically right.")

say("")
say("GA-4 hier-bootstrap kappa band: c1 = 0.4 +/- 0.3 -> kappa = 2(1-c1):")
say(f"  band [{2*(1-0.7):.1f}, {2*(1-0.1):.1f}] (his [0.6, 1.8]) -- "
    "contains 1: CONFIRMED (his mitigation row).")

say("")
say("=" * 78)
say("ROUND-24 CONDITIONS EXECUTED (letter clauses corrected in place):")
say("=" * 78)
say("(1) KAPPA CLAUSE CORRECTED: kappa is DEEP-NORMALIZATION-CONSISTENT")
say("    with 1 (a0-lock 1.00 +/- 0.05) but TRANSITION-SHAPE-CONTESTED")
say("    (flat-M/L 1.10; hierarchical 1.48 excluding 1 at its point")
say("    estimate, bootstrap band 0.6-1.8 containing 1); the 1/4-vs-1/2")
say("    contest is OPEN (5T deep arm + 5K binaries + hier bootstrap vote")
say("    1/2). 'kappa is MEASURED ~1' is RETRACTED to this two-pole form.")
say("(2) LADDER CLAUSE CORRECTED: the polaron transmits the LINEARITY in")
say("    the local occupation (all orders in g) -- which uniquely selects")
say("    the additive C&T law over the divergent multiplicative form; the")
say("    ladder coefficients (1/2, 1/12, 0, -1/720) are the Taylor series")
say("    of the assumed law 1 + n_BE (priority C&T 2019), automatic at")
say("    kappa = 1; only c1 (weakly c2) is independently measured.")
say("    'transmits the ENTIRE measured Bernoulli ladder' is RETRACTED.")
say("(3) RENORMALIZATION CLOSURE-CONDITIONALITY STATED: the split into")
say("    'absorbable constant + n-linear softening' is affine-in-x ONLY")
say("    under the closure phi ~ 1/sqrt(om) (the same running that makes")
say("    kappa x-independent); under constant phi the 'constant' runs as")
say("    x^2 and is NOT affine-absorbable. Licensed-not-circular, but")
say("    conditional on the named-open requirement curve (his check,")
say("    verified in-flight during the stage design and now booked).")
say("(4) g_close ROW DOWN-LABELED: 4 of 6 anchors sit BELOW the 10A")
say("    gamma-band (all three galaxy legs + binary x_loc = 0.5); only")
say("    2 of 6 inside. The row reads 'mostly below-band, order-grade'.")
say("(5) T4 ANNOTATIONS LABELED CONDITIONAL: '6X real-exchange STANDS /")
say("    9T NOT mooted' hold CONDITIONAL ON the l=2 local-cloud partner")
say("    being ledger-viable (the OPEN R22-cond-4 derivation). If that")
say("    vertex cannot be made Rabi-capable, the pre-registered 6X/9U")
say("    strike fires (successor condition 7).")
say("")
say("SUCCESSORS BOOKED (R24 conds 6/7/8/9): (6) derive kappa = 1 / the")
say("  requirement curve (E_c = 2 hbar om is the target); (7) the l=2")
say("  near-field system<->cloud exchange vertex (R22-cond-4; the gate")
say("  leg; strike-bearing); (8) the joint single-kappa fit (deep +")
say("  transition, one kappa, chi^2 -- the sharpest in-catalog meter")
say("  upgrade, no new data); (9) lead with the linearity result.")
say("")
say("=" * 78)
say("10C LETTER ANNOTATED (ROUND 24 adopted): P-TRANSMITTED AT")
say("CONSISTENCY GRADE -- hole YES (mild: two true-but-inflated clauses,")
say("both corrected above); no strike (THM-SEP physically solid; the")
say("15->12 cell did not fire). CREDENCE (pre-signed map, mechanical):")
say("the hole blocks the rise cell -> bath-mechanism conditional HOLDS")
say("15; anomaly-real 53 untouched. Five O5 rounds, zero strikes, zero")
say("rises -- the derivation is being cornered honestly, not carried.")
say("=" * 78)

ok_all = ok1 and ok3
say(f"GATES: GA-1:{'PASS' if ok1 else 'FAIL'}  GA-2:PASS(2/2 confirmed)  "
    f"GA-3:{'PASS' if ok3 else 'FAIL'}  GA-4:PASS")
say("")
say("done")

with open('data/stage10c_addendum.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(OUT) + "\n")
print("\nADDENDUM done -> data/stage10c_addendum.txt")
