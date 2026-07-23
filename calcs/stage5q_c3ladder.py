"""
STAGE 5Q (O13c): the c3(beta) NNLO ladder -- closed-form deep-series
coefficients of the mixing family, one more rung.

Family (5P): nu = 1 + n_BE(u), u = y^((1+beta)/2) * nu^beta. Writing
x = sqrt(y) and S = x*nu, the implicit equation becomes exactly
    S = x/2 + S^(-beta) + x^2 S^beta / 12 - x^4 S^(3 beta) / 720 + O(x^6)
(from n_BE(u) = 1/u - 1/2 + u/12 - u^3/720 + ..., u = x S^beta).
Known rungs: c1 = 1/(2(1+beta)), c2 = 1/(12(1+beta)) + beta/(8(1+beta)^2),
p_tail = (1+beta)/2, family relation c1*p_tail = 1/4.

This stage: solve the series to x^4 symbolically -> c3(beta), c4(beta)
closed forms; verify the known rungs EXACTLY; evaluate the ladder at the
member points; hunt exact relations; gate numerically with mpmath
(50-digit implicit roots, residual scaling exponent ~5).
Writes data/stage5q_c3ladder.txt.
"""
import sympy as sp
from mpmath import mp, mpf, exp, findroot

L = []

x = sp.symbols('x', positive=True)
b = sp.symbols('beta', nonnegative=True)
c1, c2, c3, c4 = sp.symbols('c1 c2 c3 c4')
S = 1 + c1*x + c2*x**2 + c3*x**3 + c4*x**4
expr = x/2 + S**(-b) + x**2*S**b/sp.Integer(12) \
       - x**4*S**(3*b)/sp.Integer(720) - S
ser = sp.series(expr, x, 0, 5).removeO().expand()

sol = {}
for k, ck in ((1, c1), (2, c2), (3, c3), (4, c4)):
    coef = sp.expand(ser.coeff(x, k).subs(sol))
    s = sp.solve(sp.Eq(coef, 0), ck)
    assert len(s) == 1, f"order x^{k}: non-unique solve"
    sol[ck] = sp.simplify(sp.together(s[0].subs(sol)))

C1_REF = 1/(2*(1+b))
C2_REF = 1/(12*(1+b)) + b/(8*(1+b)**2)
ok1 = sp.simplify(sol[c1] - C1_REF) == 0
ok2 = sp.simplify(sol[c2] - C2_REF) == 0
L.append("STAGE 5Q: the c3(beta) NNLO ladder (symbolic series of the "
         "mixing family)")
L.append("")
L.append(f"G1 known rungs reproduced exactly: c1 {'PASS' if ok1 else 'FAIL'}"
         f", c2 {'PASS' if ok2 else 'FAIL'}")
assert ok1 and ok2

c3f = sp.factor(sp.simplify(sol[c3]))
c4f = sp.factor(sp.simplify(sol[c4]))
L.append("")
L.append(f"c3(beta) = {sp.sstr(c3f)}")
L.append(f"c4(beta) = {sp.sstr(c4f)}")
L.append("")
L.append("ladder at the member points (c1, c2, c3, c4):")
for bv in (0, sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 4), 1):
    vals = [sp.nsimplify(sp.simplify(sol[ck].subs(b, bv)))
            for ck in (c1, c2, c3, c4)]
    L.append(f"  beta={sp.sstr(bv):>4}: " +
             ", ".join(sp.sstr(v) for v in vals) +
             "   (" + ", ".join(f"{float(v):+.6f}" for v in vals) + ")")

# structural checks
z3 = sp.solve(sp.Eq(sol[c3], 0), b)
L.append("")
L.append(f"zeros of c3(beta): {sp.sstr(z3)}  "
         f"(c3(0) = 0 is the pure-BE Bernoulli zero)")
p_tail = (1+b)/2
cands = {
    "c1*p_tail": sp.simplify(sol[c1]*p_tail),
    "c2*p_tail^2": sp.simplify(sp.factor(sol[c2]*p_tail**2)),
    "c3*p_tail^3": sp.simplify(sp.factor(sol[c3]*p_tail**3)),
    "c3/(beta)": sp.simplify(sol[c3]/b) if sp.simplify(sol[c3]/b) else None,
    "c2 - c1^2/2": sp.simplify(sp.factor(sol[c2] - sol[c1]**2/2)),
    "c3 - c1*c2 + c1^3/3": sp.simplify(sp.factor(
        sol[c3] - sol[c1]*sol[c2] + sol[c1]**3/3)),
}
L.append("relation hunt (simplified combos):")
for k, v in cands.items():
    if v is not None:
        L.append(f"  {k} = {sp.sstr(v)}")

# ---------------- numeric gates (mpmath, 50 digits) ----------------
mp.dps = 50

def nu_root(xv, bv):
    xv, bv = mpf(xv), mpf(bv)
    def F(nu):
        u = xv**(1+bv) * nu**bv
        return nu - 1 - 1/(exp(u) - 1)
    return findroot(F, 1/xv + mpf('0.5'))

L.append("")
L.append("G2 numeric gates (50-digit implicit roots):")
ok_all = True
for bv in ('0.25', '0.75'):
    cn = [float(sol[ck].subs(b, sp.Rational(bv))) for ck in (c1, c2, c3, c4)]
    res = {}
    for xv in ('0.02', '0.01'):
        nu = nu_root(xv, bv)
        Sn = mpf(xv)*nu
        Sp = 1 + sum(mpf(c)*mpf(xv)**(i+1) for i, c in enumerate(cn))
        res[xv] = abs(Sn - Sp)
    import math
    scal = math.log2(float(res['0.02']/res['0.01']))
    ok = 4.3 < scal < 5.7
    ok_all &= ok
    L.append(f"  beta={bv}: |resid|(x=0.02) = {float(res['0.02']):.2e}, "
             f"scaling exponent {scal:.2f} (expect ~5) "
             f"{'PASS' if ok else 'FAIL'}")
L.append(f"G2 -> {'PASS' if ok_all else 'FAIL'}")
assert ok_all

# observational size of the c3 rung
L.append("")
c3_gap = abs(float(sol[c3].subs(b, sp.Rational(1, 2))) -
             float(sol[c3].subs(b, 0)))
L.append(f"observational note: the c3 rung separates beta=0 from beta=1/2 "
         f"by |Dc3| = {c3_gap:.4f}; at x = 0.1 that is a relative shift "
         f"Dnu/nu ~ |Dc3| x^3 = {c3_gap*1e-3:.1e} -- structural algebra, "
         f"not an instrument. The measurable discriminators stay c1 "
         f"(NLO zero-point), p_tail (screening sharpness), and the "
         f"transition value nu(1).")

out = "\n".join(L)
print(out)
with open('data/stage5q_c3ladder.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage5q_c3ladder.txt")
