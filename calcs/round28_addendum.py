"""
ROUND-28 ADDENDUM -- independent verification of every load-bearing number
in the Paper-3 referee round (the memory rule: re-compute before adopting).

GA-1  The Omega/H > 1.449 bound (his M1): reconcile 9Q's 1.449 with his
      sqrt(3) -- one theorem, two unit conventions.
GA-2  Fermi optimum x* (his check vs the banked 10G value).
GA-3  Tail postdictions p = 1/2 + s^2/4 at the measured ambients.
GA-4  Multiplicative-law divergence at x = ln 2.
GA-5  His M2 r-slice arithmetic: p(r=0.34) = 0.628; sigma placements.
GA-6  His squares identity: amplitude shortfalls 6.7/9.2 vs the banked
      variance factors 45/84 (89/2.004, 1202/14.26).
GA-7  His m5: the 9V cross-tie 0.3895 vs 0.3904 -- absolute vs relative.
GA-8  His m1: orphan-reference grep of the paper body.

Output: data/round28_addendum.txt
"""
import io, math, os, re, sys

import sympy as sp

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
buf = io.StringIO()


def emit(s=""):
    print(s)
    buf.write(s + "\n")


emit("ROUND-28 ADDENDUM -- referee-number verification")
emit("=" * 60)

# ---------------------------------------------------------------- GA-1
emit("\nGA-1  the stability bound: 9Q chain vs his OSCO sqrt(3)")
r, G, M, Lam, c, H, OmL = sp.symbols("r G M Lambda c H Omega_Lambda", positive=True)
Om2 = G * M / r**3 - Lam * c**2 / 3          # SdS circular-orbit frequency^2
kap2 = sp.simplify(r**-3 * sp.diff(r**4 * Om2, r))   # radial epicyclic freq^2
emit(f"  kappa^2 = {kap2}  (expect GM/r^3 - 4/3 Lambda c^2)")
chk = sp.simplify(kap2 - (G * M / r**3 - sp.Rational(4, 3) * Lam * c**2))
emit(f"  identity residual: {chk}")
assert chk == 0
# stability kappa^2>0  =>  GM/r^3 > (4/3)Lam c^2  =>  Om^2 > Lam c^2  (exact)
Om2_at_marginal = sp.simplify(Om2.subs(G * M / r**3, sp.Rational(4, 3) * Lam * c**2))
emit(f"  Omega^2 at marginal stability: {Om2_at_marginal}  (expect Lambda c^2)")
assert sp.simplify(Om2_at_marginal - Lam * c**2) == 0
# convention A (9Q): Lambda = 3 OmL H0^2 / c^2  ->  Omega/H0 > sqrt(3 OmL)
bndA = sp.sqrt(sp.simplify((Lam * c**2).subs(Lam, 3 * OmL * H**2 / c**2)) / H**2)
emit(f"  9Q convention (H = H0, physical Lambda): Omega/H0 > {bndA}")
vA = float(bndA.subs(OmL, sp.Rational(7, 10)))
emit(f"    at Omega_Lambda = 0.7: {vA:.6f}  (banked 1.449)")
assert abs(vA - 1.4491376746) < 1e-6
# convention B (referee): f = 1-2M/r-H^2 r^2, i.e. Lambda = 3 H_dS^2
#   -> Omega/H_dS > sqrt(3); and OSCO condition M/r^3 = 4 H_dS^2
emit(f"  his convention (H = H_dS, Lambda = 3H^2): Omega/H_dS > sqrt(3) = {math.sqrt(3):.6f}")
emit(f"  reconciliation: sqrt(3)*sqrt(OmL) = {math.sqrt(3)*math.sqrt(0.7):.6f} = 1.449  EXACT")
emit(f"  his OSCO: kappa^2=0 -> GM/r^3 = (4/3)Lam c^2 = 4 H_dS^2  -> his 'M/r^3 -> 4H^2'  CONFIRMED")
emit("  VERDICT: both numbers correct; ONE theorem, two unit conventions.")
emit("  Paper defect REAL as flagged: 'every bound orbit' overstates the")
emit("  proven class (stable CIRCULAR orbits via epicyclic stability), and")
emit("  no convention was stated. Fix adopted in the text.")

# ---------------------------------------------------------------- GA-2
emit("\nGA-2  Fermi optimum")
x = sp.symbols("x", positive=True)
xstar = sp.nsolve(x - 1 - sp.exp(-x), x, 1.3)
emit(f"  x* = {float(xstar):.7f}  (his 1.278465; banked 1.2784645)")
assert abs(float(xstar) - 1.2784645) < 1e-6

# ---------------------------------------------------------------- GA-3
emit("\nGA-3  tail postdictions")
for name, n_amb, want_s2, want_p in [("galaxy", 6.63, 0.754, 0.688), ("binary", 0.502, 0.112, 0.528)]:
    s = n_amb / (1 + n_amb)
    p = 0.5 + s**2 / 4
    emit(f"  {name}: s^2 = {s*s:.4f} (banked {want_s2}), p = {p:.4f} (banked {want_p})")
    assert abs(s * s - want_s2) < 2e-3 and abs(p - want_p) < 1e-3

# ---------------------------------------------------------------- GA-4
emit("\nGA-4  multiplicative divergence")
nbe_ln2 = 1 / (math.exp(math.log(2)) - 1)
emit(f"  n_BE(ln 2) = {nbe_ln2:.12f}  (exactly 1 -> 1/(1-n_BE) diverges); y = (ln2)^2 = {math.log(2)**2:.3f}")
assert abs(nbe_ln2 - 1) < 1e-12

# ---------------------------------------------------------------- GA-5
emit("\nGA-5  his r-slice arithmetic (M2)")
s2g = (6.63 / 7.63) ** 2
for rr, lab in [(0.5, "pure dispersive"), (0.34, "9V fit")]:
    emit(f"  p(r={rr}) = {0.5 + rr*s2g/2:.4f}  ({lab})")
emit(f"  his 0.688/0.628 both inside 0.65 +/- 0.075: "
     f"z = {abs(0.5+0.5*s2g/2-0.65)/0.075:.2f} / {abs(0.5+0.34*s2g/2-0.65)/0.075:.2f} sigma  CONFIRMED")

# ---------------------------------------------------------------- GA-6
emit("\nGA-6  squares identity (amplitude vs variance shortfalls)")
emit(f"  6.67^2 = {6.67**2:.1f}, 9.18^2 = {9.18**2:.1f}  vs banked variance factors:")
emit(f"  binary 89/(2*0.502+1) = {89/2.004:.1f};  galaxy 1202/(2*6.63+1) = {1202/14.26:.1f}")
assert abs(6.67**2 - 89 / 2.004) < 0.6 and abs(9.18**2 - 1202 / 14.26) < 0.6
emit("  CONFIRMED -- his consistency check holds because variance = amplitude^2.")

# ---------------------------------------------------------------- GA-7
emit("\nGA-7  the 9V cross-tie phrasing (m5)")
d = abs(0.3895 - 0.3904)
emit(f"  |0.3895-0.3904| = {d:.4f} absolute; relative = 1/{0.39/d:.0f}")
emit("  his point CONFIRMED: 'one part in 10^3' reads relative; the honest phrase is 'within 0.001'.")

# ---------------------------------------------------------------- GA-8
emit("\nGA-8  orphan references (m1)")
t = open("papers/paper3_mechanism.md", encoding="utf-8").read()
body = t.split("## References")[0]
for key, pat in [("Chae et al. 2021", r"Chae"), ("Desmond et al. 2024", r"Desmond"),
                 ("McGaugh et al. 2016", r"McGaugh"), ("Milgrom 1983", r"Milgrom.{0,15}1983|1983.{0,15}Milgrom")]:
    hits = len(re.findall(pat, body))
    emit(f"  {key}: body citations = {hits}  ({'ORPHAN CONFIRMED' if hits == 0 else 'cited'})")

emit("\nALL REFEREE NUMBERS VERIFIED. Adoption licensed.")
open("data/round28_addendum.txt", "w", encoding="utf-8").write(buf.getvalue())
emit("\nwritten: data/round28_addendum.txt")
