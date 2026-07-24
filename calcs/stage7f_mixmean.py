"""
STAGE 7F: the mixing-mean uniqueness contest -- is the geometric mean
SELECTED by the measured tails, or merely chosen?

The 5P family interpolates between the two natural occupation arguments
of a self-dressing mode: A = sqrt(y) (source frequency; C&T/BE) and
B = nu*y (total/boot; the 6N dichotomy's other pole). The program has
always used the GEOMETRIC mean, u = A^(1-beta) B^beta -- equivalent to
the solver form y^((1+beta)/2) nu^beta (identity checked below). But the
geometric mean is one member of the power-mean family

    u_p = [ (1-beta) A^p + beta B^p ]^(1/p),   p -> 0 = geometric.

If other p work too, the mixing algebra is a modeling choice; if the
data select p = 0, seam (iv) of the derivation closes to measurement
grade. Two exact facts do the selecting:

  TAIL (the discriminator): Newtonian-ward, B/A = nu*sqrt(y) -> infinity
    -- the two arguments DIVERGE, and the mean choice is maximally
    consequential. Exact endpoint theorem (GF1): for any fixed
    beta in (0,1), p > 0 snaps the tail exponent to 1 (boot grade --
    binary-vetoed program-wide), p < 0 snaps it to 1/2 (ungated BE --
    the two-system tail SPLIT becomes impossible), and ONLY p = 0
    interpolates continuously: p_tail = (1+beta)/2 = 1/2 + g/4, the
    measured pair (0.689 gal / 0.529 bin).
  DEEP (the protection): MOND-ward, B/A -> 1 -- all power means
    degenerate to first order, and with the RUNNING gate beta ~ g x^2/8
    the p-dependence enters the ladder only at the c4 rung (GF-deep:
    a1..a3 p-independent, a4 linear in p). This is WHY the Bernoulli
    ladder never saw the mixing choice -- and why the tail must do the
    selecting.

Because the p != 0 snap happens beyond a finite crossover, the honest
in-window statement is a BOUND: the effective tail index p_eff(y; p)
over the measured window [3, 30], compared to the measured bands
(galaxy [0.65, 0.75] at g = 0.754; binary [0.45, 0.60] at g = 0.112),
gives a two-sided exclusion on the mixing index p_mix.

Gates: GF1 sympy exact endpoint limits (p in {1, 1/2, -1/2, -1} -> 1 or
1/2; p = 0 -> (1+b)/2, symbolic in b). GF2 continuity + solver
regression (p = +-1e-4 vs p = 0 within 1e-5; p = 0 vs the 6G-form
Newton AMB solver within 1e-9). GF3 constant-beta member regressions at
p = 0 (beta = 0/0.5/1 -> BE/gm/boot closed solvers, 1e-9). GF-deep
order-by-order ladder with symbolic p (a1..a3 p-free; da4/dp printed).
GF4 the geometric case's far-tail p_eff must regress to 1/2 + g/4.

PRE-REGISTERED BARS: UNIQUENESS SUPPORTED if GF1 endpoints hold AND the
window bound excludes |p_mix| >= 0.5 in at least one system with the
other consistent and p = 0 inside both bands. Deliverable = the
two-sided bound. HONEST FALLBACK: if some |p| >= 1 stays inside BOTH
bands, uniqueness is NOT establishable at current tail precision --
reported as such.
Writes data/stage7f_mixmean.txt.
"""
import math
import numpy as np
import sympy as sp
from scipy.optimize import brentq

L = []
def say(s=''):
    L.append(s); print(s, flush=True)
def save():
    with open('data/stage7f_mixmean.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

say("STAGE 7F: the mixing-mean uniqueness contest")
say("=" * 72)

G_GAL, G_BIN = 0.7540, 0.1118          # fiducial gates (6F / 6G tables)
BAND_GAL, BAND_BIN = (0.65, 0.75), (0.45, 0.60)

# ---------------- GF1: exact endpoint theorem (sympy) -----------------------
y_s, b_s = sp.symbols('y b', positive=True)
A_s, B_s = sp.sqrt(y_s), y_s            # tail: nu -> 1 so B = nu*y -> y
say("GF1 exact tail endpoints (beta = b symbolic in (0,1)):")
ok1 = True
for p in (sp.Integer(1), sp.Rational(1, 2), sp.Rational(-1, 2),
          sp.Integer(-1)):
    u = ((1 - b_s)*A_s**p + b_s*B_s**p)**(1/p)
    ex = sp.limit(sp.log(u)/sp.log(y_s), y_s, sp.oo)
    want = 1 if p > 0 else sp.Rational(1, 2)
    ok1 &= sp.simplify(ex - want) == 0
    say(f"  p = {p}: tail exponent = {ex}  (expect {want})")
u0 = A_s**(1 - b_s)*B_s**b_s
ex0 = sp.limit(sp.log(u0)/sp.log(y_s), y_s, sp.oo)
ok1 &= sp.simplify(ex0 - (1 + b_s)/2) == 0
say(f"  p = 0: tail exponent = {ex0}  (expect (1+b)/2 -> 1/2 + g/4 at the "
    f"gate's tail beta = g/2)")
say(f"GF1: {'PASS' if ok1 else 'FAIL'}")
save()
assert ok1

# ---------------- solvers ---------------------------------------------------
def n_be_scalar(u):
    if u > 40: return 0.0
    if u < 1e-12: return 1.0/u - 0.5
    return 1.0/math.expm1(u)

# AMENDMENT (post-commit a4696c1, PRE-RESULTS -- GF4 fired before the
# window table existed): the first p_eff estimator read -ln(nu-1) after
# nu-1 UNDERFLOWED at far-tail y (u ~ 2760 at y = 1e5), flatlining the
# gradient. The index is now read from the solution's own argument u,
# via the exact stable identity -ln(nu-1) = ln(e^u - 1) (= u - e^(-u)
# corrections); no underflow at any scale. Estimator fix only.
def solve_nu(y, g, p, bcap=None, want_u=False):
    """nu = 1 + n_BE(u_p(nu)) with the AMB gate beta = g/2 (2nu-1)^-2."""
    lA = 0.5*math.log(y)
    def uofnu(nu):
        beta = 0.5*g/((2.0*nu - 1.0)**2)
        lB = math.log(y) + math.log(nu)
        if p == 0.0:
            lu = (1 - beta)*lA + beta*lB
        else:
            m = max(p*lA, p*lB)
            s = (1 - beta)*math.exp(p*lA - m) + beta*math.exp(p*lB - m)
            lu = (m + math.log(s))/p
        return math.exp(min(lu, 700.0))
    def F(nu):
        return nu - 1.0 - n_be_scalar(uofnu(nu))
    lo, hi = 1.0 + 1e-14, 2.0/math.sqrt(min(y, 1.0)) + 10.0
    if F(lo)*F(hi) > 0:
        nu = 1.0 + n_be_scalar(uofnu(1.0 + 1e-14))
    else:
        nu = brentq(F, lo, hi, xtol=1e-14, rtol=1e-13)
    return (nu, uofnu(nu)) if want_u else nu

def solve_grid(ys, g, p):
    return np.array([solve_nu(y, g, p) for y in ys])

def w_grid(ys, g, p):
    """-ln(nu-1) computed stably from the solution's argument u."""
    out = []
    for y in ys:
        _, u = solve_nu(y, g, p, want_u=True)
        out.append(u if u > 30.0 else math.log(math.expm1(u)))
    return np.array(out)

# GF2 continuity + regression vs the 6G-form Newton solver
YT = np.geomspace(1e-3, 300.0, 120)
def make_amb_newton(g):
    def nu_run(yv):
        yv = np.clip(np.asarray(yv, float), 1e-14, None)
        ly = np.log(yv)
        nu = 0.5 + np.sqrt(0.25 + 1.0/yv)
        for _ in range(120):
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
            nu = np.maximum(nu - np.clip(step, -0.5*nu, 0.5*nu), 1.0 + 1e-15)
        return nu
    return nu_run
nu0 = solve_grid(YT, G_GAL, 0.0)
nuN = make_amb_newton(G_GAL)(YT)
d_newton = float(np.max(np.abs(nu0/nuN - 1.0)))
nup = solve_grid(YT, G_GAL, 1e-4)
num = solve_grid(YT, G_GAL, -1e-4)
d_cont = float(max(np.max(np.abs(nup/nu0 - 1.0)),
                   np.max(np.abs(num/nu0 - 1.0))))
gf2 = d_newton < 1e-9 and d_cont < 1e-5
say(f"GF2 solver: p=0 vs Newton-AMB max dev {d_newton:.2e}; "
    f"p=+-1e-4 continuity {d_cont:.2e} -> {'PASS' if gf2 else 'FAIL'}")
save()
assert gf2

# GF3 constant-beta members at p=0 vs closed solvers
def solve_member(y, kind):
    if kind == 'BE':
        x = math.sqrt(y)
        return 1.0 + n_be_scalar(x)
    def F(nu):
        u = y**0.75*math.sqrt(nu) if kind == 'gm' else y*nu
        return nu - 1.0 - n_be_scalar(u)
    return brentq(F, 1.0 + 1e-14, 2.0/math.sqrt(min(y, 1.0)) + 10.0,
                  xtol=1e-14)
def solve_constbeta(y, beta, p):
    lA = 0.5*math.log(y)
    def F(nu):
        lB = math.log(y) + math.log(nu)
        lu = (1 - beta)*lA + beta*lB if p == 0.0 else None
        u = math.exp(min(lu, 700.0))
        return nu - 1.0 - n_be_scalar(u)
    return brentq(F, 1.0 + 1e-14, 2.0/math.sqrt(min(y, 1.0)) + 10.0,
                  xtol=1e-14)
ok3 = True
for beta, kind in ((0.0, 'BE'), (0.5, 'gm'), (1.0, 'boot')):
    dev = max(abs(solve_constbeta(y, beta, 0.0)/solve_member(y, kind) - 1.0)
              for y in (0.03, 0.3, 1.0, 3.0, 30.0))
    ok3 &= dev < 1e-9
    say(f"GF3 member p=0 beta={beta:g} vs {kind}: max dev {dev:.1e}")
say(f"GF3: {'PASS' if ok3 else 'FAIL'}")
save()
assert ok3

# ---------------- GF-deep: the ladder protection (symbolic p) ---------------
say('')
say("GF-deep: order-by-order deep ladder with symbolic p (phi = nu*sqrt(y) "
    "= 1 + a1 x + a2 x^2 + a3 x^3 + a4 x^4):")
x, p_s, g_s = sp.symbols('x p g', positive=True)
a1, a2, a3, a4 = sp.symbols('a1 a2 a3 a4')
phi = 1 + a1*x + a2*x**2 + a3*x**3 + a4*x**4
beta = (g_s/2)*x**2/(2*phi - x)**2
lnU = (1/p_s)*sp.log((1 - beta) + beta*sp.exp(p_s*sp.log(phi)))
U = sp.exp(lnU)
# self-consistency: phi = 1/U + x/2 + (x^2/12) U - (x^4/720) U^3
rhs = 1/U + x/2 + (x**2/12)*U - (x**4/720)*U**3
eq = sp.series(phi - rhs, x, 0, 5).removeO()
sol = {}
for k in (1, 2, 3, 4):
    c = sp.expand(eq.coeff(x, k)).subs(sol)
    var = {1: a1, 2: a2, 3: a3, 4: a4}[k]
    s_ = sp.solve(sp.Eq(c, 0), var)
    sol[var] = sp.simplify(s_[0])
say(f"  a1 = {sol[a1]}")
say(f"  a2 = {sol[a2]}")
say(f"  a3 = {sp.simplify(sol[a3])}")
a4s = sp.simplify(sol[a4])
say(f"  a4 = {a4s}")
d1, d2, d3 = (sp.simplify(sp.diff(sol[v], p_s)) for v in (a1, a2, a3))
d4 = sp.simplify(sp.diff(a4s, p_s))
okd = (d1 == 0) and (d2 == 0) and (d3 == 0) and (d4 != 0)
say(f"  d a1..a3 / dp = {d1}, {d2}, {d3} (expect 0, 0, 0)")
say(f"  d a4 / dp = {d4}  (nonzero: the mixing choice first enters at c4)")
say(f"GF-deep: {'PASS' if okd else 'FAIL'} -- the Bernoulli ladder is "
    f"p-blind through c3; the tail must do the selecting")
save()
assert okd

# ---------------- GF4 + the windowed bound ----------------------------------
say('')
def p_eff(ys, g, p):
    """effective nu_p-family tail index: dln(-ln(nu-1))/dln y (stable)."""
    w = w_grid(ys, g, p)
    return np.gradient(np.log(w), np.log(ys))
YF = np.geomspace(1e3, 1e5, 40)
pf_gal = p_eff(YF, G_GAL, 0.0)[-1]
pf_bin = p_eff(YF, G_BIN, 0.0)[-1]
t_gal, t_bin = 0.5 + G_GAL/4, 0.5 + G_BIN/4
gf4 = abs(pf_gal - t_gal) < 0.01 and abs(pf_bin - t_bin) < 0.01
say(f"GF4 far-tail regression (p=0): gal {pf_gal:.4f} vs 1/2+g/4 = "
    f"{t_gal:.4f}; bin {pf_bin:.4f} vs {t_bin:.4f} -> "
    f"{'PASS' if gf4 else 'FAIL'}")
save()
assert gf4

YW = np.geomspace(3.0, 30.0, 30)
PGRID = [-2.0, -1.0, -0.5, -0.25, -0.1, 0.0, 0.1, 0.25, 0.5, 1.0, 2.0]
say('')
say("windowed effective tail index (mean over y in [3, 30]):")
say("  p_mix |  galaxy (band 0.65-0.75) |  binary (band 0.45-0.60)")
rows = {}
for p in PGRID:
    pg = float(np.mean(p_eff(YW, G_GAL, p)))
    pb = float(np.mean(p_eff(YW, G_BIN, p)))
    ing = BAND_GAL[0] <= pg <= BAND_GAL[1]
    inb = BAND_BIN[0] <= pb <= BAND_BIN[1]
    rows[p] = (pg, pb, ing, inb)
    say(f"  {p:+5.2f} |  {pg:.4f} {'IN ' if ing else 'OUT'}          |  "
        f"{pb:.4f} {'IN ' if inb else 'OUT'}")
# bounds: largest |p| range where BOTH systems stay in-band
ok_ps = [p for p in PGRID if rows[p][2] and rows[p][3]]
lo_b, hi_b = (min(ok_ps), max(ok_ps)) if ok_ps else (np.nan, np.nan)
say(f"  both-in-band range: [{lo_b:+.2f}, {hi_b:+.2f}]")
save()

zero_in = rows[0.0][2] and rows[0.0][3]
half_out = not (rows[0.5][2] and rows[0.5][3]) and \
           not (rows[-0.5][2] and rows[-0.5][3])
one_in = (rows[1.0][2] and rows[1.0][3]) or (rows[-1.0][2] and rows[-1.0][3])
if ok1 and okd and zero_in and half_out and not one_in:
    verd = ("UNIQUENESS SUPPORTED: p = 0 (geometric) inside both measured "
            "bands; |p_mix| >= 0.5 excluded both signs; exact endpoints "
            "1 / one-half proven; deep ladder p-blind through c3")
elif one_in:
    verd = ("FALLBACK: |p| >= 1 survives both bands -- uniqueness NOT "
            "establishable at current tail precision")
else:
    verd = "PARTIAL -- see table"
say('')
say(f"VERDICT (pre-registered): {verd}")
say("caveats: bands are the nu_p-family measurements (approximate mapping, "
    "same estimator family as GF4); gates g fixed at the fiducial table "
    "values; window [3, 30] pre-registered.")
save()
print("\nSTAGE 7F done -> data/stage7f_mixmean.txt")
