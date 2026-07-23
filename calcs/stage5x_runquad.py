"""
STAGE 5X (O14 Saturn leg): the derived running-beta functions against the
Cassini quadrupole -- the 4K/5I/5S machinery, one more scan.

Every function this program has fielded lands at 4.0-5.8x the Cassini cap
once its binary amplitude is applied (amplitude lock, 5I/5S). The derived
F1/F2 have gm-grade screening (p = 3/4) but stronger transitions
(nu(1) = 1.470/1.494), so their raw q should sit between gm's -0.0833 and
BE's -0.0988; the lock decides as always. Gates: Newton control, simple
regression vs 4K, hi-res doubling. alpha-hats joined from the 5W summary.
Writes data/stage5x_runquad.txt.
"""
import re
import numpy as np
from numpy.polynomial.legendre import leggauss, legval

A0_SI = 1.2e-10
R_M_SI = 1.052e15
UNIT_Q = A0_SI/R_M_SI
CASSINI_CAP = 9e-27

def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))
def nu_newton(y):
    return np.ones_like(np.asarray(y, float))

def make_runbeta(kind):
    def nu_run(y):
        y = np.clip(np.asarray(y, float), 1e-14, None)
        ly = np.log(y)
        nu = nu_simple(y)
        for _ in range(80):
            if kind == 'f1':
                b = 1.0/(2.0*nu)
                db = -1.0/(2.0*nu*nu)
            else:
                b = 1.0/(2.0*(2.0*nu - 1.0))
                db = -1.0/((2.0*nu - 1.0)**2)
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

def qval(nu, eN, NR=512, NMU=96, LMAX=16):
    r, phi_l = solve_multipole(nu, eN, NR, NMU, LMAX)
    band = (r > 0.02) & (r < 0.2)
    return float(np.median((phi_l[2]/r**2)[band]))

L = ["STAGE 5X: running-beta functions vs the Cassini quadrupole "
     "(4K/5I/5S machinery)", ""]

r, phi_l = solve_multipole(nu_newton, 1.2)
g2 = np.max(np.abs(phi_l[2]))
L.append(f"G2 Newton control: max|phi_2| = {g2:.2e} -> "
         f"{'OK' if g2 < 1e-12 else 'FAIL'}")
assert g2 < 1e-12
q_s = qval(nu_simple, 1.2)
g1 = abs(q_s/(-0.09788) - 1)
L.append(f"G1 regression simple: q = {q_s:+.5f} (4K -0.09788) "
         f"{100*g1:.1f}% -> {'OK' if g1 < 0.02 else 'FAIL'}")
assert g1 < 0.02

AH = {}
txt = open('data/stage5w_summary.txt').read()
for law in ('rb1', 'rb2'):
    ah = [float(m) for m in re.findall(rf'seed \d+ {law}: a_hat=([0-9.]+)',
                                       txt)]
    AH[law] = (float(np.mean(ah)),
               float(np.std(ah, ddof=1)/np.sqrt(len(ah))))

L.append("")
for name, kind, law in (("F1", 'f1', 'rb1'), ("F2", 'f2', 'rb2')):
    nu = make_runbeta(kind)
    q12 = qval(nu, 1.2)
    q_hi = qval(nu, 1.2, NR=768, NMU=128, LMAX=24)
    g4 = abs(q_hi/q12 - 1)
    a, da = AH[law]
    Q2_1 = 3*abs(q12)*UNIT_Q
    lock = Q2_1*a
    L.append(f"  {name}: q(1.2) = {q12:+.5f} [hi-res {100*g4:.1f}% "
             f"{'OK' if g4 < 0.05 else 'FAIL'}]; Q2(alpha=1) = {Q2_1:.2e}; "
             f"alpha-hat = {a:.3f}+-{da:.3f} -> Q2 = {lock:.2e} = "
             f"{lock/CASSINI_CAP:.1f}x Cassini")
    assert g4 < 0.05
L.append("")
L.append("(comparators: BE 4.0x, gm 4.3x, family scan 4.0-5.8x -- 5S)")

out = "\n".join(L)
print(out)
with open('data/stage5x_runquad.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage5x_runquad.txt")
