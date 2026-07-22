"""
QUMOND field solver v2: point mass in uniform external field, axisymmetric,
numerically STABLE multipole solve (scaled recurrences, no raw powers).
Gates: (G1) e_N=0 must reproduce g = nu(y) g_N exactly;
       (G2) simple-nu at e_N=1.9 should track the C&M fitting formula.
Then the Bose-Einstein nu row is the theory's true (sphericalized) EFE.
Units: GM=1, a0=1, y = 1/r^2.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss, legval, legder

NR, NMU, LMAX = 512, 96, 16
r  = np.logspace(-2, 3, NR)
mu, wmu = leggauss(NMU)
R, MU = np.meshgrid(r, mu, indexing='ij')
ST = np.sqrt(1-MU**2)

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))

Pl  = np.zeros((LMAX+1, NMU)); dPl = np.zeros((LMAX+1, NMU))
for l in range(LMAX+1):
    c = np.zeros(l+1); c[l] = 1
    Pl[l] = legval(mu, c)
    dPl[l] = legval(mu, legder(c)) if l > 0 else np.zeros(NMU)

def solve(nu, eN):
    gr = -1.0/R**2 + eN*MU
    gt = -eN*ST
    f = nu(np.hypot(gr, gt)) - 1.0
    Ar, At = f*gr, f*gt
    dAr = np.gradient(R**2*Ar, np.log(r), axis=0)/R**3
    dAt = -np.gradient(ST*At, mu, axis=1)/R
    S = -(dAr + dAt)     # source is div[(nu-1) grad(phi_N)] = -div[(nu-1) g_N]
    Sl = np.array([(2*l+1)/2*np.sum(wmu*S*Pl[l], axis=1) for l in range(LMAX+1)])
    # stable scaled integrals:
    # A_l(r) = r^-(l+1) int_0^r s^(l+2) Sl ds ;  B_l(r) = r^l int_r^inf s^(1-l) Sl ds
    A = np.zeros((LMAX+1, NR)); B = np.zeros((LMAX+1, NR))
    dr = np.diff(r)
    for l in range(LMAX+1):
        g_in = Sl[l]*r          # integrand s^(l+2) Sl / s^(l+1) evaluated s=r
        for i in range(NR-1):
            q = (r[i]/r[i+1])**(l+1)
            A[l, i+1] = A[l, i]*q + 0.5*dr[i]*(g_in[i]*q + g_in[i+1])
        g_out = Sl[l]*r
        for i in range(NR-2, -1, -1):
            q = (r[i]/r[i+1])**l
            B[l, i] = B[l, i+1]*q + 0.5*dr[i]*(g_out[i] + g_out[i+1]*q)
    phi_l = -(A+B)/(2*np.arange(LMAX+1)[:,None]+1)
    # forces (per-l analytic derivatives, boundary terms cancel)
    lv = np.arange(LMAX+1)[:,None]
    gph_r = -(( (lv+1)*A - lv*B )/(r[None,:]*(2*lv+1)))
    g_r = gr + np.tensordot(gph_r, Pl, axes=(0,0))
    g_t = gt - (ST/R)*np.tensordot(-(A+B)/(2*lv+1), dPl, axes=(0,0))
    # subtract the uniform barycenter response
    nuE = nu(np.array([max(eN, 1e-14)]))[0] if eN > 0 else 0.0
    if eN > 0:
        g_r -= nuE*eN*MU
        g_t -= -nuE*eN*ST
    boost = np.sum(wmu*(-g_r), axis=1)/2/(1.0/r**2)
    return boost

y = 1.0/r**2
win = (y > 0.05) & (y < 100)

# --- Gate 1: eN = 0, spherical identity ---
for name, nu in (("simple", nu_simple), ("BE", nu_be)):
    b0 = solve(nu, 0.0)
    exact = nu(y)
    err = np.max(np.abs(b0[win]/exact[win]-1))
    print(f"G1 ({name:>6}): max |boost/nu - 1| in window = {100*err:.2f}%")

# --- Gate 2 + production: eN = 1.9 ---
def cm_formula(yq):
    be = 1.1*1.9; yb = np.sqrt(yq*yq+be*be)
    sq = np.sqrt(0.25+1.0/yb); nus = 0.5+sq
    nuhat = (1.0/yb)/(2.0*nus*sq)
    return nus*(1.0+np.tanh((be/yq)**1.2)*nuhat/3.0)

print(f"\n{'y':>6} {'simple(solver)':>15} {'C&M formula':>12} {'BE(solver)':>11}")
b_s = solve(nu_simple, 1.9)
b_b = solve(nu_be, 1.9)
for yq in (0.3, 1.0, 3.0, 10.0, 30.0):
    i = np.argmin(np.abs(y-yq))
    print(f"{yq:>6} {b_s[i]:>15.3f} {cm_formula(yq):>12.3f} {b_b[i]:>11.3f}")
np.save('data/efe_boost_be.npy', np.stack([y, b_b]))
np.save('data/efe_boost_simple.npy', np.stack([y, b_s]))
print("\nsaved tabulated boosts -> data/efe_boost_*.npy")
