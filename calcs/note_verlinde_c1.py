"""Mechanism-consistency note (post-round-12): VERLINDE'S EMERGENT-GRAVITY
INTERPOLATION vs THE MEASURED COEFFICIENT LADDER.

Freeze-compliant: a consistency row, not a new fit — the only new content
is exact series algebra (sympy-gated); every quoted number is an existing
ledger measurement.

Verlinde (2016), in the additive apparent-DM reading tested by
Brouwer+17: g_obs = g_B + g_D with g_D = sqrt(a_M * g_B), a_M = c H0 / 6.
In this program's units (x = sqrt(g_N/a0), a0 = c H0 / 2pi):
    nu_V(x) = 1 + sqrt(r)/x,   r = a_M/a0 = (2 pi)/6 = pi/3.
Ladder read: deep amplitude sqrt(pi/3) = 1.023 (a 4.7% a0 shift — inside
the galaxy a0 error, so the AMPLITUDE is viable), constant term c1 = 1
EXACTLY, c2 = 0 EXACTLY.  The measured ladder: c1 profile 0.385-0.519
(flat M/L, 4S) / 0.208-0.309 (hierarchical, 4Z), bootstrap ~0.4 +/- 0.3,
c1 = 0 excluded at Delta(-2lnL) = 56.3 (7.5 sigma profile / 95.5%
bootstrap).  c1 = 1 sits FARTHER outside every measured interval than
the excluded 0 does.
Output: data/note_verlinde_c1.txt
"""
import sympy as sp

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

x = sp.symbols('x', positive=True)

# G1 — the BE ladder regression (the program's canonical form, exact):
nuBE = 1 + 1/(sp.exp(x)-1)
sBE = sp.series(nuBE, x, 0, 2).removeO().expand()
c_m1 = sBE.coeff(x, -1); c_0 = sBE.coeff(x, 0); c_1 = sBE.coeff(x, 1)
assert (c_m1, c_0, c_1) == (1, sp.Rational(1, 2), sp.Rational(1, 12)), \
    (c_m1, c_0, c_1)
P("G1 BE ladder: nu = 1/x + 1/2 + x/12 (exact) -> PASS")

# G2 — the simple-nu rung regression (c1 = 1/2, c2 = 1/8, exact):
nuS = sp.Rational(1, 2) + sp.sqrt(sp.Rational(1, 4) + 1/x**2)
sS = sp.series(nuS, x, 0, 2).removeO().expand()
assert (sS.coeff(x, -1), sS.coeff(x, 0), sS.coeff(x, 1)) == \
    (1, sp.Rational(1, 2), sp.Rational(1, 8)), sS
P("G2 simple-nu ladder: nu = 1/x + 1/2 + x/8 (exact) -> PASS")

# G3 — the Verlinde read (exact):
r = sp.pi/3
nuV = 1 + sp.sqrt(r)/x
sV = sp.series(nuV, x, 0, 3).removeO().expand()
cV_m1 = sV.coeff(x, -1); cV_0 = sV.coeff(x, 0)
cV_1 = sV.coeff(x, 1); cV_2 = sV.coeff(x, 2)
assert cV_m1 == sp.sqrt(sp.pi/3) and cV_0 == 1 \
    and cV_1 == 0 and cV_2 == 0, sV
P("G3 Verlinde ladder: nu = sqrt(pi/3)/x + 1 + 0*x + 0*x^2 (exact) "
  "-> PASS")
P("")
P(f"deep amplitude: sqrt(a_M/a0) = sqrt(pi/3) = "
  f"{float(sp.sqrt(sp.pi/3)):.4f} -> a0 shift "
  f"{100*(float(sp.pi/3)-1):.1f}% (a_M = cH0/6 vs a0 = cH0/2pi) — "
  f"INSIDE the galaxy a0 error (1.05 +/- 0.10 e-10): amplitude VIABLE")
P("constant term: c1(Verlinde) = 1 EXACTLY; c2 = 0 EXACTLY")
P("")
P("CONTEST vs the measured ladder (existing rows, no new fit):")
P("  c1 profile (flat M/L, 4S):     0.385 - 0.519  -> |1 - 0.45| = 0.55")
P("  c1 profile (hierarchical, 4Z): 0.208 - 0.309  -> |1 - 0.26| = 0.74")
P("  c1 bootstrap (both):           ~0.4 +/- 0.3")
P("  c1 = 0 exclusion grade: Delta(-2lnL) = 56.3 (7.5 sigma profile /")
P("    95.5% bootstrap) at distance 0.45 from the flat profile peak")
P("  => c1 = 1 sits FARTHER from every measured interval than the")
P("     excluded c1 = 0 (0.55-0.74 vs 0.26-0.45 away) — under the")
P("     family's profile curvature this is >= the c1 = 0 grade on the")
P("     opposite flank.")
P("")
P("VERDICT (consistency row): the additive emergent-gravity")
P("interpolation is EXCLUDED at LADDER level by the existing c1")
P("measurement, with its amplitude viable — the exclusion is about the")
P("digits, not the scale (the 1/6-vs-1/2pi near-coincidence, 4.7%, is")
P("why amplitude-level tests could not see this).")
P("")
P("Caveats (carried): Verlinde's derivation is spherical/quasi-static;")
P("the additive composition is the commonly tested reading (Brouwer+17")
P("lensing); disk-geometry corrections exist and are not computed here;")
P("the measured lambda grid ended at c1 = 0.625, so the exact Delta at")
P("c1 = 1 is curvature-extrapolation grade (stated, not overquoted);")
P("scout on prior coefficient-level Verlinde exclusions PENDING —")
P("no novelty claim printed.")

with open('data/note_verlinde_c1.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/note_verlinde_c1.txt")
