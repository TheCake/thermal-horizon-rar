"""
STAGE 5Z (O15b): the two-leg refinement of the spontaneous-fraction bath
-- derivation + exact algebra.

The 5W residual is localized at the transition: the binaries want LESS
admixture there than the one-leg share beta = 1/(2 nu) leaves. The
derived next step (a reading, like 5U's): the mode's SELF-frequency can
only be expressed if a quantum completes the round trip -- emitted AND
reabsorbed through the spontaneous channel. Two legs, each with
spontaneous share 1/(1+n) = 1/nu, so the response admixture SQUARES:
  F3 (two-leg rate share):   beta = (1/2) / nu^2
  F4 (two-leg energy share): beta = (1/2) / (2 nu - 1)^2
Limits unchanged: n >> 1 -> beta -> 0 (classical/source-driven);
n -> 0 -> beta -> 1/2 (quantum endpoint) => tail p = 3/4 preserved.
NEW exact consequence (pre-registered): the faster die-off now protects
BOTH Bernoulli rungs -- c1 = 1/2 AND c2 = 1/12 exactly (the one-leg
functions bent c2 to -1/6 / -1/24; the two-leg ones converge to the pure
occupation series through NNLO). Transition: nu(1) rises toward the
binaries (F3 ~ 1.50, F4 ~ 1.54 vs one-leg 1.470/1.494, BE 1.582).

Gates as 5U: G1 sympy series (c1 = 1/2 AND c2 = 1/12 exact, both);
G2 solver residual < 2e-12 + root-uniqueness scans; G3 mpmath 50-digit
series residual scaling ~ x^4; G4 numeric tail exponent = 3/4; nu(1).
Writes data/stage5z_twoleg.txt.
"""
import numpy as np
import sympy as sp
from mpmath import mp, mpf, exp as mexp, findroot

def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))

def make_twoleg(kind):
    """kind='f3': beta = 1/(2 nu^2); kind='f4': beta = 1/(2(2 nu-1)^2)."""
    def nu_run(y):
        y = np.clip(np.asarray(y, float), 1e-14, None)
        ly = np.log(y)
        nu = nu_simple(y)
        for _ in range(80):
            if kind == 'f3':
                b = 0.5/(nu*nu)
                db = -1.0/(nu**3)
            else:
                b = 0.5/((2.0*nu - 1.0)**2)
                db = -2.0/((2.0*nu - 1.0)**3)
            u = np.exp(np.minimum(0.5*(1.0+b)*ly + b*np.log(nu), 60.0))
            eu = np.exp(np.minimum(u, 60.0))
            em1 = np.maximum(eu - 1.0, 1e-300)
            n = np.where(u < 60.0, 1.0/em1, 0.0)
            F = nu - 1.0 - n
            dudnu = u*(db*(0.5*ly + np.log(nu)) + b/nu)
            dF = 1.0 + (eu/(em1*em1))*dudnu
            step = F/np.where(np.abs(dF) > 1e-12, dF, 1e-12)
            nu = np.maximum(nu - np.clip(step, -0.5*nu, 0.5*nu),
                            1.0 + 1e-15)
        return nu
    return nu_run

nu_f3 = make_twoleg('f3')
nu_f4 = make_twoleg('f4')

L = ["STAGE 5Z: the two-leg spontaneous-fraction bath -- derivation + "
     "exact algebra", ""]

# G1: symbolic series
x = sp.symbols('x', positive=True)
c1, c2, c3 = sp.symbols('c1 c2 c3')
S = 1 + c1*x + c2*x**2 + c3*x**3
LS = sp.log(S)
CO = {}
for name, bexpr in (("F3", x**2/(2*S**2)), ("F4", x**2/(2*(2*S - x)**2))):
    eq = x/2 + sp.exp(-bexpr*LS) + x**2*sp.exp(bexpr*LS)/sp.Integer(12) \
         - x**4*sp.exp(3*bexpr*LS)/sp.Integer(720) - S
    ser = sp.series(eq, x, 0, 4).removeO().expand()
    sol = {}
    for k, ck in ((1, c1), (2, c2), (3, c3)):
        coef = sp.expand(ser.coeff(x, k).subs(sol))
        s = sp.solve(sp.Eq(coef, 0), ck)
        assert len(s) == 1
        sol[ck] = sp.simplify(sp.together(s[0].subs(sol)))
    ok1 = (sol[c1] == sp.Rational(1, 2)) and (sol[c2] == sp.Rational(1, 12))
    L.append(f"G1 {name}: c1 = {sp.sstr(sol[c1])}, c2 = {sp.sstr(sol[c2])}, "
             f"c3 = {sp.sstr(sol[c3])}   (c1 = 1/2 AND c2 = 1/12 exact: "
             f"{'PASS' if ok1 else 'FAIL'})")
    assert ok1
    CO[name] = [sp.Rational(1, 2), sp.Rational(1, 12), sol[c3]]
L.append("   (one-leg comparators: c2 = -1/6 (F1) / -1/24 (F2); BE c3 = 0)")
L.append("")

# G2: solver residual + uniqueness
def resid(fn, kind, y):
    nu = fn(y)
    if kind == 'f3': b = 0.5/(nu*nu)
    else:            b = 0.5/((2.0*nu - 1.0)**2)
    u = np.exp(np.minimum(0.5*(1.0+b)*np.log(y) + b*np.log(nu), 60.0))
    n = np.where(u < 60.0, 1.0/np.expm1(np.minimum(u, 60.0)), 0.0)
    return np.max(np.abs(nu - 1.0 - n)/(1.0 + n))

yv = np.logspace(-8, 4, 400)
r3 = resid(nu_f3, 'f3', yv)
r4 = resid(nu_f4, 'f4', yv)
ok2 = max(r3, r4) < 2e-12
L.append(f"G2 solver residual: F3 {r3:.1e}, F4 {r4:.1e} -> "
         f"{'PASS' if ok2 else 'FAIL'}")
uniq = True
for yq in (1e-6, 1e-2, 0.5, 1.0, 3.0, 50.0):
    for kind in ('f3', 'f4'):
        nug = np.linspace(1.0 + 1e-9,
                          float(nu_simple(np.array([yq]))[0])*3, 20000)
        if kind == 'f3': b = 0.5/(nug*nug)
        else:            b = 0.5/((2.0*nug - 1.0)**2)
        u = np.exp(np.minimum(0.5*(1.0+b)*np.log(yq) + b*np.log(nug), 60.0))
        n = np.where(u < 60.0, 1.0/np.expm1(np.minimum(u, 60.0)), 0.0)
        F = nug - 1.0 - n
        uniq &= int(np.sum(np.diff(np.sign(F)) != 0)) == 1
L.append(f"G2b root-uniqueness scan -> {'PASS' if uniq else 'FAIL'}")
assert ok2 and uniq

# G3: mpmath 50-digit
mp.dps = 50
L.append("")
import math
for name, kind, C in (("F3", 'f3', CO["F3"]), ("F4", 'f4', CO["F4"])):
    cn = [float(c) for c in C]
    res = {}
    for xv in ('0.02', '0.01'):
        xm = mpf(xv)
        def Fm(nu):
            if kind == 'f3': b = 1/(2*nu**2)
            else:            b = 1/(2*(2*nu - 1)**2)
            u = xm**(1+b) * nu**b
            return nu - 1 - 1/(mexp(u) - 1)
        nu = findroot(Fm, 1/xm + mpf('0.5'))
        Sn = xm*nu
        Sp = 1 + sum(mpf(c)*xm**(i+1) for i, c in enumerate(cn))
        res[xv] = abs(Sn - Sp)
    scal = math.log2(float(res['0.02']/res['0.01']))
    ok3 = 3.3 < scal < 4.7
    L.append(f"G3 {name}: |resid|(x=0.02) = {float(res['0.02']):.2e}, "
             f"scaling {scal:.2f} (expect ~4) {'PASS' if ok3 else 'FAIL'}")
    assert ok3

# G4: tail + nu(1)
L.append("")
for name, fn in (("F3", nu_f3), ("F4", nu_f4)):
    yw = np.array([20.0, 60.0])
    nm = fn(yw) - 1.0
    p_hat = np.diff(np.log(-np.log(nm)))[0]/np.diff(np.log(yw))[0]
    ok4 = abs(p_hat - 0.75) < 0.03
    n1 = float(fn(np.array([1.0]))[0])
    L.append(f"G4 {name}: tail p_hat = {p_hat:.4f} (pred 0.7500) "
             f"{'PASS' if ok4 else 'FAIL'};  nu(1) = {n1:.4f}  "
             f"(one-leg 1.4702/1.4943, BE 1.5820)")
    assert ok4

L.append("")
L.append("status: two-leg spontaneity (round-trip emission+reabsorption) "
         "squares the response share; BOTH Bernoulli rungs now exact, "
         "tail p = 3/4 preserved, transition raised toward the binaries. "
         "Galaxy leg = 6A, binary leg = 6B.")

out = "\n".join(L)
print(out)
with open('data/stage5z_twoleg.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage5z_twoleg.txt")
