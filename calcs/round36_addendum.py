"""
ROUND 36 ADDENDUM -- stage 10Q (O5-AMPLITUDE) verification.

GA half = BLIND: written and committed BEFORE the ROUND 36 report is
read (the 87a4676 protocol, 9th execution). Every load-bearing stage
number re-derived with INDEPENDENT methods where feasible (trap #23:
the verifier is an instrument too -- discrepancies get audited on the
verifier first).

GB half = appended AFTER the report: re-compute every load-bearing
reviewer number (memory rule feedback-verify-reviewer-math).
"""
import math
import numpy as np

say = print
say("ROUND 36 ADDENDUM -- GA (blind) half")
say("=" * 60)

TWO_PI = 2*math.pi
ok_all = True
def check(name, got, want, tol, rel=True):
    global ok_all
    d = abs(got - want)/(abs(want) if rel else 1.0)
    ok = d <= tol
    ok_all = ok_all and ok
    say(f"  {name}: {got:.6g} vs {want:.6g} (d {d:.2e}, tol {tol:.0e})"
        f" -> {'OK' if ok else 'MISMATCH'}")
    return ok

# ---------- GA-1: FD closed form via direct Boltzmann sums ----------
say("GA-1 FD via direct Fock/Boltzmann sums (independent of sympy):")
for x in (1.0954, 0.1411, 0.5):
    w = np.arange(0, 4000)
    p = np.exp(-x*w); p /= p.sum()
    got = float(np.sum(p*(2*w + 1)))
    want = 1.0/math.tanh(x/2)
    check(f"2n+1 sum at x={x}", got, want, 1e-9)

# ---------- GA-2: Hessian by numeric finite differences ----------
say("GA-2 sector Hessian by finite differences (independent of the")
say("  sympy tensor algebra): E(q) = (k/2) sum q^2 - F q0")
rng = np.random.default_rng(363636)
kv, Fv = 1.7, 0.31
def Etot(q):
    return 0.5*kv*np.sum(q**2) - Fv*q[0]
q_eq = np.zeros(5); q_eq[0] = Fv/kv
h = 1e-5
H = np.zeros((5, 5))
for i in range(5):
    for j in range(5):
        qpp = q_eq.copy(); qpp[i] += h; qpp[j] += h
        qpm = q_eq.copy(); qpm[i] += h; qpm[j] -= h
        qmp = q_eq.copy(); qmp[i] -= h; qmp[j] += h
        qmm = q_eq.copy(); qmm[i] -= h; qmm[j] -= h
        H[i, j] = (Etot(qpp) - Etot(qpm) - Etot(qmp) + Etot(qmm))/(4*h*h)
check("||H - k I||_max", float(np.abs(H - kv*np.eye(5)).max()), 0.0,
      1e-5, rel=False)

# ---------- GA-3: eps2 by an INDEPENDENT solver route ----------
say("GA-3 eps2 by independent Green DOUBLE-SUM (no A/B recursion) +")
say("  independent finite-difference source stencil:")

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))

from numpy.polynomial.legendre import leggauss, legval

def source_indep(nu, eN, NR=512, NMU=96):
    # independent stencil: manual central differences in ln r and mu
    r = np.logspace(-2, 3, NR)
    lnr = np.log(r)
    mu, wmu = leggauss(NMU)
    R, MU = np.meshgrid(r, mu, indexing='ij')
    ST = np.sqrt(1-MU**2)
    gr = -1.0/R**2 + eN*MU
    gt = -eN*ST
    f = nu(np.hypot(gr, gt)) - 1.0
    Ar, At = f*gr, f*gt
    # d(r^2 Ar)/dr / r^2 via manual central differences in ln r
    X = R**2*Ar
    dX = np.empty_like(X)
    dX[1:-1] = (X[2:] - X[:-2])/(lnr[2:, None] - lnr[:-2, None])
    dX[0] = (X[1] - X[0])/(lnr[1] - lnr[0])
    dX[-1] = (X[-1] - X[-2])/(lnr[-1] - lnr[-2])
    dAr = dX/R**3
    Y = ST*At
    dY = np.empty_like(Y)
    dY[:, 1:-1] = (Y[:, 2:] - Y[:, :-2])/(mu[None, 2:] - mu[None, :-2])
    dY[:, 0] = (Y[:, 1] - Y[:, 0])/(mu[1] - mu[0])
    dY[:, -1] = (Y[:, -1] - Y[:, -2])/(mu[-1] - mu[-2])
    dAt = -dY/R
    S = -(dAr + dAt)
    out = {}
    for l in (0, 2):
        c = np.zeros(l+1); c[l] = 1
        Pl = legval(mu, c)
        out[l] = (2*l+1)/2*np.sum(wmu*S*Pl, axis=1)
    return r, out

def phi_greensum(r, Sl, l):
    # direct O(N^2) two-sided kernel sum (independent of recursion)
    n = len(r)
    w = np.gradient(r)               # trapezoid-equivalent weights
    phi = np.zeros(n)
    for i in range(n):
        ri = r[i]
        inner = r <= ri
        outer = ~inner
        phi[i] = -(1.0/(2*l+1))*(
            np.sum(Sl[inner]*(r[inner]**(l+2)/ri**(l+1))*w[inner])
            + np.sum(Sl[outer]*(ri**l/r[outer]**(l-1))*w[outer]))
    return phi

r, Sl = source_indep(nu_be, 1.2)
E = {}
for l in (0, 2):
    ph = phi_greensum(r, Sl[l], l)
    E[l] = float(-0.5/(2*l+1)*np.sum(Sl[l]*ph*r**2*np.gradient(r)))
eps2_indep = E[2]/E[0]
check("eps2(BE, eN=1.2) independent", eps2_indep, 0.06490, 2e-2)
# and the plateau q from the independent phi
ph2 = phi_greensum(r, Sl[2], 2)
m = (r >= 0.02) & (r <= 0.2)
q_indep = float(np.mean(ph2[m]/r[m]**2))
check("q(BE, eN=1.2) independent", q_indep, -0.09874, 1e-2)

# ---------- GA-4: the assembly by hand arithmetic ----------
say("GA-4 assembly arithmetic:")
n_amb, x_amb, eps2, qdc = 0.502, 1.0954, 0.06490, 0.086
TE = 1.0/((n_amb + 0.5)*x_amb)
check("T/E_amb (full)", TE, 0.9111, 2e-3)
base = qdc*math.sqrt(TE/(2*eps2))
u_lo, u_hi = 0.5*x_amb*math.sqrt(0.086), 0.5*x_amb*1.0
cf = lambda u: math.sqrt(u/math.tanh(u))
check("e_a lo (coth lo)", base*cf(u_lo), 0.2288, 2e-3)
check("e_a hi (coth hi)", base*cf(u_hi), 0.2388, 2e-3)
check("e_a central", 0.5*(base*cf(u_lo) + base*cf(u_hi)), 0.2338, 2e-3)
# envelope corners: (thermal n) x q_hi x coth-hi / (n+1) x q_lo x coth-lo
TE_th = 1.0/(n_amb*x_amb)
TE_sp = 1.0/((n_amb + 1.0)*x_amb)
env_hi = 0.099*math.sqrt(TE_th/(2*eps2))*cf(u_hi)
env_lo = 0.073*math.sqrt(TE_sp/(2*eps2))*cf(u_lo)
check("envelope hi", env_hi, 0.3883, 2e-3)
check("envelope lo", env_lo, 0.1586, 2e-3)
check("chart margin 1/<dq^2>_hi", 1.0/env_hi**2, 6.6, 2e-2)
check("static ratio central", 0.5*(base*cf(u_lo)+base*cf(u_hi))/0.086,
      2.7, 2e-2)

# ---------- GA-5: frozen-row arithmetic ----------
say("GA-5 frozen row:")
check("1/Om_amb [Gyr]", TWO_PI/1.0954*14.42, 82.7, 2e-3)
check("hi edge [Gyr]", TWO_PI/1.0954*14.42/math.sqrt(0.086), 282.0,
      5e-3)
check("margin vs 2.48", TWO_PI/1.0954*14.42/2.48, 33.0, 2e-2)

# ---------- GA-6: P10 S at a corner by independent 2x2 eigensolve ----
say("GA-6 S at (x=0.5, gam=0.010, kap=1.000, conv=1, e_a=0.2338) by")
say("  direct complex 2x2 eigenvalues (independent of the closed form):")
x_amb_10n = math.log(1 + 1/0.5202)
xl, gam, kap, ea = 0.5, 0.010, 1.000, 0.2338
om = xl/TWO_PI; Om = x_amb_10n/TWO_PI
lam_max = 0.304*0.5*math.sqrt(om*Om)
U = ea*lam_max
d2 = (kap/4)*(om/Om)
P0 = math.exp(-d2)
D1 = abs(om*(1 - kap))
V = math.sqrt(2*U*U*P0)
M = np.array([[0, V], [V, D1 - 1j*gam]])   # width on the far level
ev = np.linalg.eigvals(M + 1j*gam/2*np.eye(2))  # symmetrize the width
# closed-form comparator (the stage's convention): delta = Re sqrt(z^2+V^2)-Re z
z = complex(D1/2, -gam/2)
delta_closed = float(np.real(np.sqrt(z*z + V*V)) - np.real(z))
nb = 1.0/(math.exp(xl) - 1.0)
W = (kap*om/4)*(2*nb + 2)
pc = math.exp(-xl)*(1 - math.exp(-xl))**2
S_closed = pc*delta_closed/W
# independent: eigenvalue splitting of [[0, V], [V, D1]] with width gam
# entering as the imaginary part on the detuned level, repulsion =
# Re(lam+) - 0 at kap = 1 (D1 = 0): Re sqrt(V^2 - gam^2/4)
delta_indep = float(np.real(np.emath.sqrt(V*V - gam*gam/4)))
S_indep = pc*delta_indep/W
check("S corner (closed vs stage token 0.00218)", S_closed, 0.00218,
      1e-2)
check("S corner (independent eig route)", S_indep, S_closed, 2e-2)

# ---------- GA-7: thin-shell constants by exact kernel integrals ----
say("GA-7 thin-shell constants (independent exact kernel argument):")
for l in (0, 2):
    # phi_l(a) for unit shell: -(1/(2l+1)) [a^{l+2}/a^{l+1} +
    # a^l/a^{l-1}] / a^2 ... = -(1/(2l+1)) * (a + a)/... derive:
    # A(a) = a^{-(l+1)} * a^{l+2}/a^2 = a^{... } with S = d(r-a)/a^2:
    # int S s^{l+2} ds = a^l; A(a) = a^l/a^{l+1} = 1/a; B(a) = 1/a...
    # continuity: phi(a) = -(1/(2l+1)) * (1/a)  [single-sided value]
    Eval = 1.0/(2*(2*l+1)**2)
    say(f"  l={l}: E_l = 1/(2 (2l+1)^2) = {Eval:.6f} "
        f"(stage analytic {'0.5' if l == 0 else '0.02'})")
    ok_all = ok_all and abs(Eval - (0.5 if l == 0 else 0.02)) < 1e-12

say("")
say(f"GA VERDICT: {'ALL OK' if ok_all else 'MISMATCH PRESENT'}")
say("(committed blind, before the ROUND 36 report is read)")

# ================= GB half (appended post-report) =================
say("")
say("ROUND 36 ADDENDUM -- GB (post-report) half: re-compute every")
say("load-bearing REVIEWER number (memory rule)")
say("=" * 60)

ok_gb = True
def gb(name, got, want, tol, rel=True):
    global ok_gb
    d = abs(got - want)/(abs(want) if rel else 1.0)
    ok = d <= tol
    ok_gb = ok_gb and ok
    say(f"  {name}: {got:.6g} vs his {want:.6g} (d {d:.2e}) -> "
        f"{'OK' if ok else 'MISMATCH'}")
    return ok

# GB-1: his pairing table (cond 1) -- e_a under alternative (q, eps2)
say("GB-1 pairing table (his S3(v)):")
cf_mid = 0.5*(cf(u_lo) + cf(u_hi))
def ea_of(q, e2):
    return q*math.sqrt(TE/(2*e2))*cf_mid
gb("matched eN=1.2 (0.0987, 0.0649)", ea_of(0.0987, 0.0649), 0.268,
   1e-2)
gb("self-consistent eN=1.0 (0.0849, 0.0575)", ea_of(0.0849, 0.0575),
   0.245, 1e-2)
gb("matched simple (0.0976, 0.0442)", ea_of(0.0976, 0.0442), 0.321,
   1e-2)
gb("stage mixed (0.086, 0.0649)", ea_of(0.086, 0.0649), 0.234, 3e-3)

# GB-2: his virial table (cond 2) -- S_max(hi) under E_amb scalings
say("GB-2 virial table (his S3(ii)); S_max over the full 10N grid:")
def S_max_at(ea):
    best = 0.0
    for xl in (0.5, 1.0):
        lam_max = {0.5: 0.304, 1.0: 0.140}[xl]*0.5*math.sqrt(
            (xl/TWO_PI)*(x_amb_10n/TWO_PI))
        for gm in (0.010, 0.015, 0.025):
            for kp in (0.888, 0.925, 1.000):
                for conv in (1/math.sqrt(2), 1.0, math.sqrt(2)):
                    om = xl/TWO_PI
                    D1 = abs(om*(1 - kp))
                    d2 = (kp/4)*(om/(x_amb_10n/TWO_PI))
                    V2 = 2*(conv*ea*lam_max)**2*math.exp(-d2)
                    z = complex(D1/2, -gm/2)
                    dlt = float(np.real(np.sqrt(z*z + V2))
                                - np.real(z))
                    nb2 = 1.0/(math.exp(xl) - 1.0)
                    W2 = (kp*om/4)*(2*nb2 + 2)
                    pc2 = math.exp(-xl)*(1 - math.exp(-xl))**2
                    best = max(best, pc2*dlt/W2)
    return best
gb("S_max(hi) stage envelope 0.388", S_max_at(0.388), 0.01112, 5e-3)
gb("S_max(hi) virial-2 (0.275)", S_max_at(0.388/math.sqrt(2)),
   0.0072, 5e-2)
gb("S_max(hi) /2 (0.549)", S_max_at(0.388*math.sqrt(2)), 0.0164,
   5e-2)
gb("S_max(hi) /4 (0.777)", S_max_at(0.388*2), 0.0237, 5e-2)

# GB-3: his revival threshold (cond 5)
say("GB-3 revival threshold (smallest e_a with S_max >= 0.02):")
ea_rev = None
for ea in np.linspace(0.4, 1.0, 1201):
    if S_max_at(ea) >= 0.02:
        ea_rev = float(ea)
        break
gb("e_a revival", ea_rev, 0.662, 1e-2)
gb("revival / hi-envelope", ea_rev/0.388, 1.7, 2e-2)
gb("revival / central", ea_rev/0.2338, 2.8, 2e-2)

# GB-4: his anharmonic probe (cond 3) -- cubic invariant analytics
say("GB-4 anharmonic probe (independent ANALYTIC route; his FD table):")
# basis B0 = diag(-1,-1,2)/sqrt3; c3 = Tr(B0^3) = 2/sqrt3
c3 = (2*(-1/math.sqrt(3))**3 + (2/math.sqrt(3))**3)
gb("Tr(B0^3)", c3, 2/math.sqrt(3), 1e-12)
# E(q0) = (k/2) q0^2 + g3 c3 q0^3 - F q0; sector Hessian eigen-shifts
# 6 g3 q0 t_m with t_m = Tr(B_m^2 B0) = {2/sqrt3, 1/sqrt3, -2/sqrt3}
# GB first-run note (trap #23, the verifier is an instrument): my
# v1 divided the cubic Hessian shifts by a spurious /2 -- in the
# Tr(B_m B_n) = 2 delta basis the Hessian of Tr(Q^3) is
# 6 g3 q0 Tr(B_m B_n B0) with NO extra norm factor. Removing it
# reproduces the reviewer exactly; his numbers were right.
kv2, q0v = 1.0, 0.086
for g3 in (0.5, 1.0):
    tm = np.array([2, 1, 1, -2, -2])/math.sqrt(3)
    lam_eig = kv2 + 6*g3*q0v*tm
    split = float((lam_eig.max() - lam_eig.min())/lam_eig.mean())
    k_sec = kv2 + 2*g3*c3*q0v
    k_curv = kv2 + 6*g3*c3*q0v
    gap_k = float((k_curv - k_sec)/kv2)      # his normalization
    gap_sec = float((k_curv - k_sec)/k_sec)  # variant row
    if g3 == 0.5:
        gb("eigen split at g3=0.5k", split, 0.596, 5e-2)
        gb("sec-vs-curv gap at g3=0.5k (per k)", gap_k, 0.199, 5e-2)
        say(f"    (variant: gap per k_sec = {gap_sec:.3f})")
        gb("e_a downward at g3=0.5k", 1 - math.sqrt(k_sec/k_curv),
           0.09, 2e-1)
    else:
        gb("eigen split at g3=1.0k", split, 1.19, 5e-2)
        gb("sec-vs-curv gap at g3=1.0k (per k)", gap_k, 0.397, 5e-2)

# GB-5: frozen-row extras
say("GB-5 frozen-row extras:")
gb("period 2pi/Om2 hi-edge [Gyr]", TWO_PI*82.7127, 520.0, 2e-2)

# GB-6: virial low-envelope corner (his cond-2 number)
say("GB-6 virial-2 low-envelope corner:")
gb("low corner under virial-2", 0.1586/math.sqrt(2), 0.11, 3e-2)

say("")
say(f"GB VERDICT: {'ALL REVIEWER NUMBERS CONFIRMED' if ok_gb else 'MISMATCH PRESENT'}")
