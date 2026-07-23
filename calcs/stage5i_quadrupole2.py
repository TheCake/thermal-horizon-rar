"""
STAGE 5I (O10a): the solar quadrupole for the sharpened-screening functions.

4K found Q2 = 3.9e-26 (alpha/1.15) s^-2 for BOTH 1/2-branch families
(transition-sourced) -- 4.3x over the Cassini cap 9e-27 (2 sigma, Hees+14).
5G/5H promoted nu_p(0.65) (hier galaxies +56 over BE, binaries accept,
alpha-hat 1.45 +/- ~0.1 across 2 seeds) and O5 proposes the geometric-mean
bootstrap nu_gm: nu = 1 + n_BE(y^(3/4) sqrt(nu)) (omega = sqrt(omega_source
* omega_total); series nu*sqrt(y) = 1 + x/3 + x^2/12 -- c1=1/3, c2=1/12).
Both have SHARPER Newtonian screening and WEAKER transition boosts
(nu(1): 1.423 / 1.433 vs BE 1.582) -- the quadrupole is transition-sourced,
so both directions matter. This computes their q via the 4K multipole
machinery (verbatim) and reports Q2 at alpha=1 and at the fitted amplitudes.

Gates: G1 regression -- simple at eN=1.2 must reproduce 4K's stored
q = -0.09788 within 2%; G2 Newton control; G3 plateau <5%; G4 resolution
doubling <5% (both new functions); G5 series gate for nu_gm (c1=1/3,
c2=1/12) and its solver residual < 1e-9.
Writes data/stage5i_quadrupole2.txt.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss, legval

A0_SI = 1.2e-10
R_M_SI = 1.052e15
UNIT_Q = A0_SI/R_M_SI
CASSINI_CAP = 9e-27          # 2-sigma, Hees+14 (3+/-3)e-27

def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))
def nu_newton(y):
    return np.ones_like(np.asarray(y, float))
def nu_p065(y):
    yc = np.clip(np.asarray(y, float), 1e-14, None)
    ex = np.exp(-np.minimum(yc**0.65, 60.0))
    return (1.0-ex)**(-1.0/1.3)
def nu_gm(y):
    """geometric-mean bootstrap: nu = 1 + n_BE(y^(3/4) sqrt(nu)).
    Newton on H(w) = w^2 - 1 - n_BE(y^(3/4) w), w = sqrt(nu); H monotone
    increasing in w > 0 => unique root; seed w0 = sqrt(nu_simple)."""
    y = np.clip(np.asarray(y, float), 1e-14, None)
    a = y**0.75
    w = np.sqrt(nu_simple(y))
    for _ in range(30):
        u = np.minimum(a*w, 60.0)
        eu = np.exp(u)
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        H = w*w - 1.0 - n
        dH = 2.0*w + a*eu/(em1*em1)
        w = np.maximum(w - H/dH, 1e-8)
    return w*w

def solve_multipole(nu, eN, NR=512, NMU=96, LMAX=16):
    r = np.logspace(-2, 3, NR)
    mu, wmu = leggauss(NMU)
    R, MU = np.meshgrid(r, mu, indexing='ij')
    ST = np.sqrt(1-MU**2)
    Pl = np.zeros((LMAX+1, NMU))
    for l in range(LMAX+1):
        c = np.zeros(l+1); c[l] = 1
        Pl[l] = legval(mu, c)
    gr = -1.0/R**2 + eN*MU
    gt = -eN*ST
    f = nu(np.hypot(gr, gt)) - 1.0
    Ar, At = f*gr, f*gt
    dAr = np.gradient(R**2*Ar, np.log(r), axis=0)/R**3
    dAt = -np.gradient(ST*At, mu, axis=1)/R
    S = -(dAr + dAt)
    Sl = np.array([(2*l+1)/2*np.sum(wmu*S*Pl[l], axis=1)
                   for l in range(LMAX+1)])
    A = np.zeros((LMAX+1, NR)); B = np.zeros((LMAX+1, NR))
    dr = np.diff(r)
    for l in range(LMAX+1):
        g_in = Sl[l]*r
        for i in range(NR-1):
            qf = (r[i]/r[i+1])**(l+1)
            A[l, i+1] = A[l, i]*qf + 0.5*dr[i]*(g_in[i]*qf + g_in[i+1])
        g_out = Sl[l]*r
        for i in range(NR-2, -1, -1):
            qf = (r[i]/r[i+1])**l
            B[l, i] = B[l, i+1]*qf + 0.5*dr[i]*(g_out[i] + g_out[i+1]*qf)
    lv = np.arange(LMAX+1)[:, None]
    phi_l = -(A+B)/(2*lv+1)
    return r, phi_l

L = ["STAGE 5I: solar quadrupole for nu_p(0.65) and the geometric-mean "
     "bootstrap (4K machinery verbatim)",
     f"units: q unit a0/r_M = {UNIT_Q:.3e} s^-2; Q2 = 3|q|*unit*alpha; "
     f"Cassini 2-sigma cap {CASSINI_CAP:.0e}"]

# G5: nu_gm series + solver residual
xg = np.linspace(0.005, 0.05, 40)
yg = xg*xg
h = nu_gm(yg)*xg - 1.0
cfit = np.polyfit(xg, h/xg, 1)
c1n, c2n = cfit[1], cfit[0]
yv = np.logspace(-8, 2, 300)
w_ = np.sqrt(nu_gm(yv))
u_ = np.minimum(yv**0.75*w_, 60.0)
n_ = np.where(u_ < 60, 1.0/np.expm1(u_), 0.0)
res = np.max(np.abs(w_*w_ - 1.0 - n_)/(1.0 + n_))   # relative: deep n ~ 1e8
n1 = float(nu_gm(np.array([1.0]))[0])
g5ok = abs(c1n-1/3) < 2e-3 and abs(c2n-1/12) < 4e-3 and res < 2e-12
L.append(f"G5 nu_gm: c1 = {c1n:.4f} (pred 0.3333), c2 = {c2n:.4f} "
         f"(pred {1/12:.4f}); max|H| = {res:.1e}; nu(1) = {n1:.3f} -> "
         f"{'PASS' if g5ok else 'FAIL'}")

# G2 Newton control
r, phi_l = solve_multipole(nu_newton, 1.2)
g2 = np.max(np.abs(phi_l[2]))
L.append(f"G2 Newton control: max|phi_2| = {g2:.2e} -> "
         f"{'OK' if g2 < 1e-12 else 'FAIL'}")

# G1 regression: simple at eN=1.2 vs 4K stored value
r, phi_l = solve_multipole(nu_simple, 1.2)
q_s = float(np.median((phi_l[2]/r**2)[(r > 0.02) & (r < 0.2)]))
Q4K_SIMPLE = -0.09788
g1 = abs(q_s/Q4K_SIMPLE - 1)
L.append(f"G1 regression simple eN=1.2: q = {q_s:+.5f} (4K {Q4K_SIMPLE}) "
         f"{100*g1:.1f}% -> {'OK' if g1 < 0.02 else 'FAIL'}")
L.append("")

results = {}
for name, nu in (("p065", nu_p065), ("gm", nu_gm)):
    for eN in (1.0, 1.2, 1.4):
        r, phi_l = solve_multipole(nu, eN)
        q_r = phi_l[2]/r**2
        band = (r > 0.02) & (r < 0.2)
        q = float(np.median(q_r[band]))
        flat = float(np.max(np.abs(q_r[band]/q - 1)))
        row = f"  {name:>5} eN={eN}: q = {q:+.5f}  plateau dev {100*flat:.1f}%"
        if eN == 1.2:
            _, phi_hi = solve_multipole(nu, eN, NR=768, NMU=128, LMAX=24)
            r_hi = np.logspace(-2, 3, 768)
            q_hi = float(np.median((phi_hi[2]/r_hi**2)[(r_hi > 0.02)
                                                       & (r_hi < 0.2)]))
            row += (f"  [G4 hi-res {q_hi:+.5f}, {100*abs(q_hi/q-1):.1f}% "
                    f"{'OK' if abs(q_hi/q-1) < 0.05 else 'FAIL'}]")
        results[(name, eN)] = q
        L.append(row)

L.append("")
L.append("Q2 [s^-2] = 3|q|*unit*alpha  (4K comparators: simple/BE q ~ 0.098"
         " -> Q2(alpha=1) = 3.35e-26):")
ALPHAS = {"p065": (1.45, 0.11), "gm": (None, None)}
for name in ("p065", "gm"):
    q12 = results[(name, 1.2)]
    dq_gext = abs(results[(name, 1.4)]-results[(name, 1.0)])/0.4*0.05
    Q2_1 = 3*abs(q12)*UNIT_Q
    a, da = ALPHAS[name]
    L.append(f"  {name:>5}: q(1.2) = {q12:+.5f}; Q2(alpha=1) = {Q2_1:.2e} "
             f"(vs 1/2-branch 3.35e-26: x{Q2_1/3.35e-26:.2f}); "
             f"g_ext +/-0.05a0 -> +/-{3*dq_gext*UNIT_Q:.1e}")
    L.append(f"         vs Cassini cap: {Q2_1/CASSINI_CAP:.1f}x at alpha=1"
             + (f"; at alpha={a}+/-{da}: Q2 = {Q2_1*a:.2e} = "
                f"{Q2_1*a/CASSINI_CAP:.1f}x cap"
                if a else "  (alpha from the 5K binary fit; Q2 scales linearly)"))

out = "\n".join(L)
print(out)
with open('data/stage5i_quadrupole2.txt', 'w') as f_:
    f_.write(out+"\n")
print("\nsaved: data/stage5i_quadrupole2.txt")
