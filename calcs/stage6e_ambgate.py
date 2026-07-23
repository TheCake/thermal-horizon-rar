"""
STAGE 6E (O5/O16 theory): the ambient-gated bath -- the system-level
rule derived, and its exact algebra.

The 6B/6D end state to explain: sparse ambient (binaries, n_amb ~ 0.5)
-> NO admixture (pure occupation law at every separation); occupied
ambient (galaxy environments, n_amb ~ 7) -> FULL two-leg admixture; the
gate is system-level, not pointwise. The naive rule (ambient traffic
DILUTES the spontaneous share) fails on SIGN: weak fields mean deep,
highly occupied ambients. The reading that carries the right sign, in
the same Einstein-coefficient grammar as 5U/5Z: a mode can only dress
itself to its self-consistent frequency if the ambient reservoir can
assist the exchange, and reservoir assistance is a STIMULATED process:

    beta = (1/2) * [local zero-point share]^2 * [ambient stimulated share]^2
         = (1/2) * [ (1/2)/(n_loc + 1/2) ]^2 * [ n_amb/(1 + n_amb) ]^2
         = g_amb * (1/2)/(2 nu - 1)^2 ,   g_amb = [n_amb/(1+n_amb)]^2

    n_amb = n_BE(sqrt(e_N/a0))  (source-driven ambient occupation)

Admixture = (local quantumness) x (ambient classicality). Zero
parameters. POST-HOC STATUS FLAGGED: constructed after the 6B/6D
pattern; its falsifiable content = the exact consequences below + the
galaxy leg surviving the g_amb = 0.754 dilution (6F) + the binary
acceptance numbers (6G) + out-of-sample Chae/DR4 legs.

Exact consequences (gated below, fixed before the fits):
  - c1 = 1/2 and c2 = 1/12 EXACT for every g_amb (the constant prefactor
    only enters at the c3 rung: c3 = c3_BE-ladder - g_amb/16)
  - tail p(n_amb) = 1/2 + (1/4) g_amb:
      galaxies (e_N = 0.02 a0, n_amb = 6.58, g = 0.754) -> p = 0.688
        [the 5G/5T measured band 0.65-0.75]
      binaries (e_N = 1.15 a0, n_amb = 0.520, g = 0.117) -> p = 0.529
        [their held 1/2]
  - nu(1): 1.548 (galaxy gate) / 1.577 (binary gate; 0.3% from BE)
Gates: G1 sympy series (c1, c2 exact at symbolic g); G2 solver residual
+ uniqueness; G3 mpmath 50-digit; G4 numeric tail exponents vs the
formula at both gates; G5 limits (g=1 -> F4, g=0 -> BE).
Writes data/stage6e_ambgate.txt.
"""
import math
import numpy as np
import sympy as sp
from mpmath import mp, mpf, exp as mexp, findroot

def n_amb_of(eN_over_a0):
    x = math.sqrt(eN_over_a0)
    return 1.0/(math.exp(x) - 1.0)

def g_of(n_amb):
    return (n_amb/(1.0 + n_amb))**2

N_GAL = n_amb_of(0.02)
N_BIN = n_amb_of(1.15)
G_GAL = g_of(N_GAL)
G_BIN = g_of(N_BIN)

def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))

def make_amb(g):
    def nu_run(y):
        y = np.clip(np.asarray(y, float), 1e-14, None)
        ly = np.log(y)
        nu = nu_simple(y)
        for _ in range(80):
            b = g*0.5/((2.0*nu - 1.0)**2)
            db = g*(-2.0)/((2.0*nu - 1.0)**3)
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

L = ["STAGE 6E: the ambient-gated bath -- beta = g_amb * (1/2)/(2nu-1)^2,"
     " g_amb = [n_amb/(1+n_amb)]^2, n_amb = n_BE(sqrt(e_N/a0))", "",
     f"gates of record: n_amb(gal, e=0.02) = {N_GAL:.3f} -> g = {G_GAL:.4f}"
     f";  n_amb(bin, e=1.15) = {N_BIN:.4f} -> g = {G_BIN:.4f}", ""]

# G1: symbolic series at symbolic gate g
x = sp.symbols('x', positive=True)
gs = sp.symbols('g', positive=True)
c1, c2, c3 = sp.symbols('c1 c2 c3')
S = 1 + c1*x + c2*x**2 + c3*x**3
LS = sp.log(S)
bexpr = gs*x**2/(2*(2*S - x)**2)
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
L.append(f"G1 series at symbolic g: c1 = {sp.sstr(sol[c1])}, "
         f"c2 = {sp.sstr(sol[c2])}, c3 = {sp.sstr(sol[c3])} -> "
         f"{'PASS' if ok1 else 'FAIL'} (c1, c2 g-independent)")
assert ok1

# G2: solver residual + uniqueness at both gates
def resid(fn, g, y):
    nu = fn(y)
    b = g*0.5/((2.0*nu - 1.0)**2)
    u = np.exp(np.minimum(0.5*(1.0+b)*np.log(y) + b*np.log(nu), 60.0))
    n = np.where(u < 60.0, 1.0/np.expm1(np.minimum(u, 60.0)), 0.0)
    return np.max(np.abs(nu - 1.0 - n)/(1.0 + n))

yv = np.logspace(-8, 4, 400)
ok2 = True
for tag, g in (("gal", G_GAL), ("bin", G_BIN)):
    fn = make_amb(g)
    rr = resid(fn, g, yv)
    ok2 &= rr < 2e-12
    L.append(f"G2 solver residual ({tag}, g={g:.3f}): {rr:.1e}")
uniq = True
for yq in (1e-6, 1e-2, 1.0, 50.0):
    for g in (G_GAL, G_BIN):
        nug = np.linspace(1.0 + 1e-9,
                          float(nu_simple(np.array([yq]))[0])*3, 20000)
        b = g*0.5/((2.0*nug - 1.0)**2)
        u = np.exp(np.minimum(0.5*(1.0+b)*np.log(yq) + b*np.log(nug), 60.0))
        n = np.where(u < 60.0, 1.0/np.expm1(np.minimum(u, 60.0)), 0.0)
        uniq &= int(np.sum(np.diff(np.sign(nug - 1.0 - n)) != 0)) == 1
L.append(f"G2b uniqueness -> {'PASS' if (ok2 and uniq) else 'FAIL'}")
assert ok2 and uniq

# G3: mpmath at the galaxy gate
mp.dps = 50
cn = [0.5, 1.0/12.0, float(sol[c3].subs(gs, G_GAL))]
res = {}
for xv in ('0.02', '0.01'):
    xm = mpf(xv)
    gg = mpf(G_GAL)
    def Fm(nu):
        b = gg/(2*(2*nu - 1)**2)
        u = xm**(1+b) * nu**b
        return nu - 1 - 1/(mexp(u) - 1)
    nu = findroot(Fm, 1/xm + mpf('0.5'))
    Sn = xm*nu
    Sp = 1 + sum(mpf(c)*xm**(i+1) for i, c in enumerate(cn))
    res[xv] = abs(Sn - Sp)
scal = math.log2(float(res['0.02']/res['0.01']))
ok3 = 3.3 < scal < 4.7
L.append(f"G3 mpmath (gal gate): |resid|(x=0.02) = {float(res['0.02']):.2e},"
         f" scaling {scal:.2f} -> {'PASS' if ok3 else 'FAIL'}")
assert ok3

# G4: numeric tail exponents vs p = 1/2 + g/4
L.append("")
for tag, g in (("gal", G_GAL), ("bin", G_BIN)):
    fn = make_amb(g)
    yw = np.array([20.0, 60.0])
    nm = fn(yw) - 1.0
    p_hat = np.diff(np.log(-np.log(nm)))[0]/np.diff(np.log(yw))[0]
    p_pred = 0.5 + 0.25*g
    ok4 = abs(p_hat - p_pred) < 0.03
    n1 = float(fn(np.array([1.0]))[0])
    L.append(f"G4 {tag}: tail p_hat = {p_hat:.4f} (pred {p_pred:.4f}) "
             f"{'PASS' if ok4 else 'FAIL'};  nu(1) = {n1:.4f}")
    assert ok4

# G5: limits (g=1 is F4 by definition; g=0 must collapse to pure BE)
f0 = make_amb(0.0)
dv_be = float(np.max(np.abs(f0(yv)/nu_be(yv) - 1.0)))
ok5 = dv_be < 1e-9
L.append(f"G5 limits: g=0 vs BE {dv_be:.1e} -> {'PASS' if ok5 else 'FAIL'}"
         f" (g=1 is F4 by definition)")
assert ok5

L.append("")
L.append("predictions of record (pre-registered for 6F/6G): galaxy leg "
         "must hold >= F3-grade under the 0.754 dilution; binary leg "
         "lands within ~ -2 of p050 with BE-grade a0 (~ +2 sigma); "
         "out-of-sample: Chae high-e_N outskirts soften with g(n_amb(e)), "
         "DR4 weak-ambient binaries sharpen toward p = 0.69.")

out = "\n".join(L)
print(out)
with open('data/stage6e_ambgate.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage6e_ambgate.txt")
