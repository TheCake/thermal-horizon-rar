"""
ROUND 33 ADDENDUM (2026-08-10): verification for the 10M+10N pair.
GA half (BLIND, committed BEFORE the round report exists -- 4th
execution of the 87a4676 protocol): re-derive every load-bearing
number the two stages print, via INDEPENDENT routes wherever one
exists.  GB half (post-report): re-compute every number the ROUND-33
report introduces.  Output: data/round33_addendum.txt.
"""
import math
import numpy as np
import sympy as sp

OUT = "data/round33_addendum.txt"
L = []
def emit(s=""):
    print(s)
    L.append(s)

ok_all = True
def gate(name, cond, detail=""):
    global ok_all
    ok_all = ok_all and bool(cond)
    emit("%s: %s%s" % (name, "PASS" if cond else "FAIL",
                       ("  " + detail) if detail else ""))

emit("ROUND 33 ADDENDUM -- GA half (blind, pre-report)")
emit("")

# ---------- GA-1: 10M operator identities, independent route ----------
# Completing the square: H_n = Om (B'+d)(B+d) + [om(n+sig) - g^2(n+sig)^2/Om]
# as a MATRIX identity, plus midpoint theorem by central difference.
rng = np.random.default_rng(33)
Mb = 120
bd = np.sqrt(np.arange(1, Mb))
ab = np.diag(bd, 1)
qb = ab + ab.T
nbm = np.diag(np.arange(Mb, dtype=float))
worst = 0.0
for _ in range(5):
    omv = float(rng.uniform(0.5, 2.0))
    Omv = float(rng.uniform(0.5, 2.0))
    gv = float(rng.uniform(0.1, 0.7))
    nn = int(rng.integers(0, 4))
    ch = nn + 0.5
    Hn = omv*ch*np.eye(Mb) + Omv*nbm + gv*ch*qb
    d = gv*ch/Omv
    lhs = Omv*(ab.T + d*np.eye(Mb)) @ (ab + d*np.eye(Mb)) \
        + (omv*ch - gv*gv*ch*ch/Omv)*np.eye(Mb)
    worst = max(worst, float(np.max(np.abs(Hn - lhs))))
gate("GA-1a completing-the-square matrix identity", worst <= 1e-10,
     "worst %.1e" % worst)
# midpoint theorem via central difference at random draws
worst = 0.0
for _ in range(5):
    omv = float(rng.uniform(0.5, 2.0))
    Omv = float(rng.uniform(0.5, 2.0))
    gv = float(rng.uniform(0.1, 0.7))
    nn = int(rng.integers(0, 5))
    E = lambda x: omv*(x + 0.5) - gv*gv*(x + 0.5)**2/Omv
    D = E(nn + 1) - E(nn)
    om_mf_mid = omv - 2*gv*gv*(nn + 0.5 + 0.5)/Omv
    worst = max(worst, abs(D - om_mf_mid))
gate("GA-1b midpoint theorem (5 random draws)", worst <= 1e-12,
     "worst %.1e" % worst)
# naive offset and ratio, independent algebra at numbers
omv, Omv, gv, nn = 1.3, 0.9, 0.45, 3
kapv = 4*gv*gv/(Omv*omv)
D = (omv*(nn + 1.5) - gv*gv*(nn + 1.5)**2/Omv) \
    - (omv*(nn + 0.5) - gv*gv*(nn + 0.5)**2/Omv)
om_nv = omv - 4*gv*gv*(nn + 0.5)/Omv
gate("GA-1c naive offset == kap om n/2",
     abs((D - om_nv) - kapv*omv*nn/2) <= 1e-12,
     "d %.1e" % abs((D - om_nv) - kapv*omv*nn/2))

# ---------- GA-2: 10M edge map, independent ----------
XS = [0.03, 0.05, 0.1, 0.2, 0.3]
L10 = [0.0149, 0.0247, 0.0488, 0.0952, 0.1395]
kc = [2.0/(1.0/(math.exp(x) - 1) + 0.5) for x in XS]
rats = [k/l for k, l in zip(kc, L10)]
# endpoint log-slope (independent of polyfit)
s_h = math.log(kc[-1]/kc[0])/math.log(XS[-1]/XS[0])
s_l = math.log(L10[-1]/L10[0])/math.log(XS[-1]/XS[0])
gate("GA-2a edge ratios in [4.0, 4.3]",
     all(4.0 <= r <= 4.3 for r in rats),
     "band [%.2f, %.2f]" % (min(rats), max(rats)))
gate("GA-2b endpoint slopes ~ (1.00, 0.97)",
     abs(s_h - 0.997) <= 0.02 and abs(s_l - 0.972) <= 0.02,
     "%.3f / %.3f" % (s_h, s_l))
kmax_members = max([2.0/(1.0/(math.exp(x) - 1) + 0.5) for x in XS]
                   + [2.0*(math.exp(x) - 1) for x in XS]
                   + [1.0/(1.0/(math.exp(x) - 1) + 0.5) for x in XS])
gate("GA-2c below-lock max == 0.700", abs(kmax_members - 0.700) <= 0.002,
     "%.3f" % kmax_members)

# ---------- GA-3: 10M reachability, independent propagator ----------
from scipy.linalg import expm
Ms, Mb2 = 6, 18
Ns = np.diag(np.arange(Ms, dtype=float))
bd2 = np.sqrt(np.arange(1, Mb2))
ab2 = np.diag(bd2, 1)
qb2 = ab2 + ab2.T
nb2 = np.diag(np.arange(Mb2, dtype=float))
kv, Omv, omv = 1.2, 0.8, 1.0
gv = math.sqrt(kv*Omv*omv/4)
H = (omv*np.kron(Ns, np.eye(Mb2)) + Omv*np.kron(np.eye(Ms), nb2)
     + gv*np.kron(Ns + 0.5*np.eye(Ms), qb2))
U = expm(-1j*H*7.3)
psi = np.zeros(Ms*Mb2, dtype=complex)
cohs = np.array([1.2**k/math.sqrt(math.factorial(k)) for k in range(Ms)])
cohs /= np.linalg.norm(cohs)
coh2 = np.array([0.7**k/math.sqrt(math.factorial(k))
                 for k in range(Mb2)])
coh2 /= np.linalg.norm(coh2)
psi = np.kron(cohs, coh2)
Nop = np.kron(Ns, np.eye(Mb2))
e0 = float(np.real(psi @ (Nop @ psi)))
pt = U @ psi
e1 = float(np.real(np.conj(pt) @ (Nop @ pt)))
gate("GA-3 <N_s> drift under scipy expm propagator (kappa=1.2)",
     abs(e1 - e0) <= 1e-9, "drift %.1e" % abs(e1 - e0))

# ---------- GA-4: 10N response, independent implementation ----------
TWO_PI = 2*math.pi
x_amb_bin = math.log(1 + 1/0.5202)
def S_indep(x_loc, kap, gam, e_a, conv=1.0):
    om = x_loc/TWO_PI
    Om = x_amb_bin/TWO_PI
    gcl = 0.5*math.sqrt(om*Om)
    lam = {0.5: 0.304, 1.0: 0.140}[x_loc]*gcl
    U = conv*e_a*lam
    d2 = (kap/4)*(om/Om)
    P0 = math.exp(-d2)
    D1 = abs(om*(1 - kap))
    # independent form: delta = (sqrt(Dc^2 + 8 U^2 P0) - Dc)/2 with
    # Dc = D1 - i gam  (equivalent algebra, different code path)
    Dc = complex(D1, -gam)
    delta = float(np.real((np.sqrt(Dc*Dc + 8*U*U*P0) - Dc)/2))
    e = math.exp(-x_loc)
    pcon = e*(1 - e)**2
    nbar = 1.0/(math.exp(x_loc) - 1)
    W = (kap*om/4)*(2*nbar + 2)
    return pcon*delta/W
# stage values to reproduce (from data/stage10n_crossing.txt)
CHECKS = [
    (0.5, 0.010, 1.000, 1.0, 0.02145),
    (0.5, 0.010, 0.925, 1.0, 0.02050),
    (0.5, 0.015, 1.000, 1.0, 0.02081),
    (1.0, 0.015, 0.925, 1.0, 0.01049),
    (1.0, 0.025, 1.000, 1.0, 0.00863),
    (0.5, 0.010, 1.000, 0.3, 0.00467),
]
worst = 0.0
for xl, gm, kp, ea, tok in CHECKS:
    v = S_indep(xl, kp, gm, ea)
    worst = max(worst, abs(v - tok))
gate("GA-4a kill-table spot checks (6 cells, independent code path)",
     worst <= 6e-5, "worst |d| %.1e" % worst)
# FC factor via numeric displacement operator
d2 = (1.0/4)*((0.5/TWO_PI)/(x_amb_bin/TWO_PI))
Mq = 80
bdq = np.sqrt(np.arange(1, Mq))
aq = np.diag(bdq, 1)
Dop = expm(math.sqrt(d2)*(aq.T - aq))
fc_num = float(np.abs(Dop[0, 0])**2)
gate("GA-4b FC factor numeric vs e^{-d^2}",
     abs(fc_num - math.exp(-d2)) <= 1e-10,
     "d %.1e" % abs(fc_num - math.exp(-d2)))
# e_a* at the earned corner by bisection (vs the stage's grid)
lo, hi = 0.01, 1.0
for _ in range(60):
    mid = 0.5*(lo + hi)
    if S_indep(0.5, 1.000, 0.010, mid) >= 0.02:
        hi = mid
    else:
        lo = mid
gate("GA-4c e_a*(2%) at the earned corner by bisection",
     abs(hi - 0.936) <= 0.002, "%.4f vs 0.936" % hi)

# ---------- GA-5: 10N wiring, from first principles ----------
gc_gal = 0.5*math.sqrt((0.0953/TWO_PI)*(math.log(1 + 1/6.583)/TWO_PI))
gate("GA-5a g_close(deep-gal) = 0.00924 H",
     abs(gc_gal - 0.00924) <= 5e-5, "%.5f" % gc_gal)
lam05 = 0.304*0.5*math.sqrt((0.5/TWO_PI)*(x_amb_bin/TWO_PI))
lam10 = 0.140*0.5*math.sqrt((1.0/TWO_PI)*(x_amb_bin/TWO_PI))
gate("GA-5b lambda_max table (0.01771, 0.01154) H",
     abs(lam05 - 0.01771) <= 5e-5 and abs(lam10 - 0.01154) <= 5e-5,
     "%.5f / %.5f" % (lam05, lam10))
g_gal = (6.583/(1 + 6.583))**2
g_bin = (0.5202/(1 + 0.5202))**2
gate("GA-5c 6E gate tokens (0.7536 / 0.1171)",
     abs(g_gal - 0.7536) <= 1e-3 and abs(g_bin - 0.1171) <= 1e-3,
     "%.4f / %.4f" % (g_gal, g_bin))

# ---------- GA-6: 10N window function, exact ----------
x = sp.Symbol('x', positive=True)
w = sp.exp(-x)*(1 - sp.exp(-x))**2
xstar = sp.solve(sp.diff(w, x), x)
xs = [s for s in xstar if s.is_real and s > 0]
gate("GA-6 window argmax == ln 3, value 4/27",
     len(xs) == 1 and sp.simplify(xs[0] - sp.log(3)) == 0
     and sp.simplify(w.subs(x, sp.log(3)) - sp.Rational(4, 27)) == 0,
     "x* = %s" % xs)

# ---------- GA-7: 10N structural rows ----------
pc_deep = math.exp(-0.0953)*(1 - math.exp(-0.0953))**2
gate("GA-7a deep-gal contrast suppression x20",
     abs(pc_deep/(4/27) - 0.051) <= 0.002,
     "%.4f of window max" % (pc_deep/(4/27)))
Om_gal = math.log(1 + 1/6.583)/TWO_PI
gate("GA-7b gal-transition comb: gamma/Om spans 1",
     0.010/Om_gal < 1.0 <= 0.025/Om_gal,
     "gamma/Om = %.2f-%.2f" % (0.010/Om_gal, 0.025/Om_gal))
# peak location, independent: maximize S over kappa by golden section
def S_pk(kp):
    return S_indep(1.0, kp, 0.015, 0.086)
a, b = 0.80, 0.999
gr = (math.sqrt(5) - 1)/2
c, dd = b - gr*(b - a), a + gr*(b - a)
for _ in range(80):
    if S_pk(c) < S_pk(dd):
        a = c
    else:
        b = dd
    c, dd = b - gr*(b - a), a + gr*(b - a)
kstar = 0.5*(a + b)
gate("GA-7c peak kappa* = 0.8918 (golden-section, independent)",
     abs(kstar - 0.8918) <= 0.003, "%.4f" % kstar)

# ---------- GA-8: hygiene ----------
import io, os, re
ok8 = True
for fn, must in (("data/stage10m_hartree.txt", "M-FROZEN"),
                 ("data/stage10n_crossing.txt", "N-ROW-EARNED")):
    t = io.open(fn, encoding="utf-8", errors="replace").read()
    ok8 = ok8 and (must in t)
    for banned in ("kappa measured", "same world", "antigravity"):
        ok8 = ok8 and (banned not in t.lower())
    ok8 = ok8 and ("FAIL" not in t)
gate("GA-8 outputs exist, letters present, no banned strings, no "
     "FAIL rows", ok8)

emit("")
emit("GA half: %s" % ("ALL PASS" if ok_all else "FAILURES PRESENT"))
emit("")
emit("(GB half appended post-report: every ROUND-33 number re-computed)")

with open(OUT, "w") as f:
    f.write("\n".join(L) + "\n")
