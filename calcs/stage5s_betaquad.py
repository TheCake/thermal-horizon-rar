"""
STAGE 5S (O13d): the beta-family quadrupole scan + the amplitude-lock test.

4K found the solar quadrupole Q2 ~ 4.3x the Cassini cap for both 1/2-branch
families; 5I found the same for nu_p(0.65) and the geometric mean AFTER the
binary alpha-hat rescaling -- the amplitude lock: |q| falls with sharper
screening but alpha-hat rises in proportion, so Q2*alpha stays put. This
stage runs the scan over the MIXING DIAL: q(beta) for beta = 0, 1/4, 3/4, 1
(1/2 = gm regression-gated vs 5I), then the lock row Q2(alpha-hat(beta))
using the 5R binary amplitudes. If the lock holds across the whole family,
no choice of beta escapes Cassini either -- the MI/EFE-screening doors
remain the only ones (4K/5I conclusion, family-completed).

Gates: G2 Newton control; G1 simple regression (4K q = -0.09788, <2%);
G-member: nu_beta(0.5) reproduces the 5I gm q = -0.08327 (<0.5%); G4
resolution doubling (<5%) for each new function. Multipole machinery =
4K/5I verbatim. Solver results cached in data/stage5s_q.npy so the lock
join can re-run after the 5R fits land.
Writes data/stage5s_betaquad.txt.
"""
import os, re
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

def nu_beta(y, beta):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    A = y**(0.5*(1.0+beta))
    nu = nu_simple(y)
    for _ in range(40):
        u = np.minimum(A*nu**beta, 60.0)
        eu = np.exp(u)
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        F = nu - 1.0 - n
        dF = 1.0 + (eu/(em1*em1))*u*beta/nu
        nu = np.maximum(nu - F/dF, 1.0 + 1e-15)
    return nu

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

CACHE = 'data/stage5s_q.npy'
BETAS = [0.0, 0.25, 0.75, 1.0]
L = ["STAGE 5S: beta-family quadrupole scan (4K/5I machinery verbatim)",
     f"units: q unit a0/r_M = {UNIT_Q:.3e} s^-2; Q2 = 3|q|*unit*alpha; "
     f"Cassini 2-sigma cap {CASSINI_CAP:.0e}", ""]

if os.path.exists(CACHE):
    Q = np.load(CACHE, allow_pickle=True).item()
    L.append("(solver results loaded from cache)")
else:
    Q = {}
    # G2 Newton control
    r, phi_l = solve_multipole(nu_newton, 1.2)
    g2 = np.max(np.abs(phi_l[2]))
    L.append(f"G2 Newton control: max|phi_2| = {g2:.2e} -> "
             f"{'OK' if g2 < 1e-12 else 'FAIL'}")
    assert g2 < 1e-12
    # G1 simple regression vs 4K
    q_s = qval(nu_simple, 1.2)
    g1 = abs(q_s/(-0.09788) - 1)
    L.append(f"G1 regression simple eN=1.2: q = {q_s:+.5f} (4K -0.09788) "
             f"{100*g1:.1f}% -> {'OK' if g1 < 0.02 else 'FAIL'}")
    assert g1 < 0.02
    # member regression: beta=0.5 vs 5I gm
    q_gm = qval(lambda y: nu_beta(y, 0.5), 1.2)
    gm5i = abs(q_gm/(-0.08327) - 1)
    L.append(f"G-member beta=0.5 vs 5I gm: q = {q_gm:+.5f} (5I -0.08327) "
             f"{100*gm5i:.2f}% -> {'OK' if gm5i < 0.005 else 'FAIL'}")
    assert gm5i < 0.005
    Q[0.5] = {1.2: q_gm}
    for b in BETAS:
        nu = lambda y, b=b: nu_beta(y, b)
        Q[b] = {eN: qval(nu, eN) for eN in (1.0, 1.2, 1.4)}
        q_hi = qval(nu, 1.2, NR=768, NMU=128, LMAX=24)
        g4 = abs(q_hi/Q[b][1.2] - 1)
        L.append(f"  beta={b:4.2f}: q(1.0/1.2/1.4) = "
                 f"{Q[b][1.0]:+.5f} {Q[b][1.2]:+.5f} {Q[b][1.4]:+.5f}  "
                 f"[G4 hi-res {q_hi:+.5f}, {100*g4:.1f}% "
                 f"{'OK' if g4 < 0.05 else 'FAIL'}]")
        assert g4 < 0.05
    np.save(CACHE, Q)

# ---------------- the amplitude-lock join ----------------
# binary alpha-hats: members from 5O/5K six-seed means; b025/b075 from 5R
AHAT = {0.0: (1.078, 0.023, '6 seeds, 5O'),
        0.5: (1.36, 0.07, '6 seeds, 5K/5O')}
if os.path.exists('data/stage5r_summary.txt'):
    txt = open('data/stage5r_summary.txt').read()
    for law, b in (('b025', 0.25), ('b075', 0.75), ('boot', 1.0)):
        ah = [float(m) for m in re.findall(
            rf'seed \d+ {law}: a_hat=([0-9.]+)', txt)]
        ed = re.findall(rf'seed \d+ {law}: a_hat=[0-9.]+ \(grid [0-9.]+, '
                        rf'interior=(\w+)\)', txt)
        if law == 'boot':      # 5F ran 31/101; 5R adds the rest
            ah = [2.00, 2.00] + ah
            ed = ['False', 'False'] + ed
        if ah:
            n_int = sum(e == 'True' for e in ed)
            AHAT[b] = (float(np.mean(ah)),
                       float(np.std(ah, ddof=1)/np.sqrt(len(ah)))
                       if len(ah) > 1 else 0.0,
                       f'{len(ah)} seeds, interior {n_int}/{len(ah)}')

L.append("")
L.append("the amplitude lock across the mixing dial "
         "(Q2 = 3|q(1.2)|*unit*alpha-hat):")
L.append("   beta    q(1.2)    Q2(alpha=1)   alpha-hat        "
         "Q2(alpha-hat)   x Cassini")
for b in (0.0, 0.25, 0.5, 0.75, 1.0):
    q12 = Q[b][1.2]
    Q2_1 = 3*abs(q12)*UNIT_Q
    if b in AHAT:
        a, da, note = AHAT[b]
        lock = Q2_1*a
        edge = ('EDGE-RIDDEN (shape-rejected); lock n/a as measurement'
                if (b == 1.0 and a >= 1.99) else '')
        L.append(f"  {b:5.2f}  {q12:+.5f}    {Q2_1:.2e}   "
                 f"{a:.3f}+-{da:.3f} ({note})   {lock:.2e}   "
                 f"{lock/CASSINI_CAP:.1f}x {edge}")
    else:
        L.append(f"  {b:5.2f}  {q12:+.5f}    {Q2_1:.2e}   [5R pending]")
L.append("")
L.append("g_ext sensitivity: |q(1.4)-q(1.0)|/0.4 * 0.05 a0 shifts Q2 by:")
for b in BETAS:
    dq = abs(Q[b][1.4]-Q[b][1.0])/0.4*0.05
    L.append(f"  beta={b:4.2f}: +-{3*dq*UNIT_Q:.1e}")

out = "\n".join(L)
print(out)
with open('data/stage5s_betaquad.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage5s_betaquad.txt")
