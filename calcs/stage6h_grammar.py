"""
STAGE 6H (O5, second pass -- the grammar formalized): beta = beta_max *
[q_loc * s_amb]^L, with the leg count L promoted to a MEASURABLE.

6E's rule beta = (1/2)*[(1/2)/(n_loc+1/2)]^2*[n_amb/(1+n_amb)]^2 has
three structural choices. This stage formalizes what selects each and
turns the exponent into a measurement:

(i) BRANCH SELECTION (derivation-grade anchor): the dispersive frequency
pull of a system coupled to a mode of occupation n is proportional to
(2n+1) = 2(n+1/2) -- vacuum 1/2 included -- exactly (JC block, G1).
A frequency-dressing weight therefore carries the ZERO-POINT share
1/(2n+1) = 1/(2nu-1), NOT the emission-rate share 1/(n+1) = 1/nu.
That is the a-posteriori reason the F2/F4 (energy-share) lineage beat
F1/F3 (rate-share) on the controlled galaxy treatment, and it fixes the
local factor of the 6E rule. Reading-grade remains: that horizon-mode
dressing is dispersive.

(ii) the ambient share n/(1+n) = the stimulated fraction of reservoir
contact (6E's sign lesson, unchanged; reading-grade).

(iii) THE LEG COUNT: write beta = (1/2)*[q_loc*s_amb]^L. L=2 is the 6E
round trip. The deep ladder DERIVES a distinct signature per L (G2,
symbolic s): the leg count is read off the deepest SURVIVING Bernoulli
rung --
    L=1 breaks c2:  c2 = 1/12 - s/8          (Bernoulli c2 KILLED)
    L=2 breaks c3:  c2 = 1/12 survives, c3 = -s^2/16   (= 6E)
    L=3 breaks c4:  c2 AND c3 = 0 survive; break first at c4
and the tail runs p = 1/2 + s^L/4 (G3). The 5T decomposition (the
ultra-deep arm votes FOR the Bernoulli deep ladder) is the instrument
that can measure L -- executed in 6I as a three-way hier contest at the
fiducial gate. NEW exact rung for the record: c4(L=2, s) (gate:
c4(s->0) = -1/720 = 5Q's BE value, all L).

Also computed: per-galaxy Chae+21 gates for the 6I measured-ambient leg
(G5) and the PRE-REGISTERED 6I bars, printed here before any 6I fit.
Writes data/stage6h_grammar.txt.
"""
import csv, math
import numpy as np
import sympy as sp

def n_amb_of(e):
    x = math.sqrt(e)
    return 1.0/(math.exp(x) - 1.0)
def g_of(n):
    return (n/(1.0 + n))**2

OUT = ["STAGE 6H: the grammar formalized -- beta = (1/2)*[q_loc*s_amb]^L;"
       " the leg count as a measurable", ""]

# ---------- G1: branch selection -- the dispersive pull counts (2n+1)
lam, Dl, nn, wc = sp.symbols('lambda Delta n omega_c', positive=True)
Ee = (nn + sp.Rational(1, 2))*wc + sp.sqrt(Dl**2/4 + lam**2*(nn + 1))
Eg = (nn - sp.Rational(1, 2))*wc - sp.sqrt(Dl**2/4 + lam**2*nn)
Om = sp.expand(Ee - Eg)
pull = sp.simplify(sp.series(Om, lam, 0, 3).removeO() - (wc + Dl))
coef = sp.simplify(pull*Dl/lam**2)
ok1 = sp.simplify(coef - (2*nn + 1)) == 0
OUT.append(f"G1 dispersive pull (exact JC manifold block, O(lambda^2)): "
           f"delta_omega_q = lambda^2*({sp.sstr(coef)})/Delta -> "
           f"{'PASS' if ok1 else 'FAIL'}")
assert ok1
OUT.append("   = 2(n + 1/2): vacuum 1/2 included -- an ENERGY-type "
           "weight, not a rate-type (n+1) weight. The local gate of a "
           "frequency-dressing process is the zero-point share "
           "1/(2n+1) = 1/(2nu-1): selects the F2/F4 lineage "
           "(derivation-grade); 'horizon dressing is dispersive' stays "
           "reading-grade.")
OUT.append("")

# ---------- G2: the L-ladders at symbolic s
x = sp.symbols('x', positive=True)
ss = sp.symbols('s', positive=True)
c1, c2, c3, c4 = sp.symbols('c1 c2 c3 c4')
lad = {}
for Lc in (1, 2, 3):
    S = 1 + c1*x + c2*x**2 + c3*x**3 + c4*x**4
    LS = sp.log(S)
    bexpr = ss**Lc * x**Lc / (2*(2*S - x)**Lc)
    eq = x/2 + sp.exp(-bexpr*LS) + x**2*sp.exp(bexpr*LS)/sp.Integer(12) \
         - x**4*sp.exp(3*bexpr*LS)/sp.Integer(720) - S
    ser = sp.series(eq, x, 0, 5).removeO().expand()
    sol = {}
    for k, ck in ((1, c1), (2, c2), (3, c3), (4, c4)):
        coefk = sp.expand(ser.coeff(x, k).subs(sol))
        s_ = sp.solve(sp.Eq(coefk, 0), ck)
        assert len(s_) == 1
        sol[ck] = sp.simplify(sp.together(s_[0].subs(sol)))
    lad[Lc] = sol
    OUT.append(f"G2 L={Lc}: c1 = {sp.sstr(sol[c1])}, c2 = "
               f"{sp.sstr(sol[c2])}, c3 = {sp.sstr(sol[c3])}, "
               f"c4 = {sp.sstr(sol[c4])}")

okA = all(sol[c1] == sp.Rational(1, 2) for sol in lad.values())
okB = sp.simplify(lad[1][c2] - (sp.Rational(1, 12) - ss/8)) == 0
okC = (lad[2][c2] == sp.Rational(1, 12)) and \
      sp.simplify(lad[2][c3] + ss**2/16) == 0
okD = (lad[3][c2] == sp.Rational(1, 12)) and sp.simplify(lad[3][c3]) == 0
c40 = [sp.simplify(sol[c4].subs(ss, 0)) for sol in lad.values()]
okE = all(v == sp.Rational(-1, 720) for v in c40)
OUT.append(f"G2 gates: c1 = 1/2 all L {'PASS' if okA else 'FAIL'}; "
           f"L=1 c2-break = 1/12 - s/8 {'PASS' if okB else 'FAIL'}; "
           f"L=2 = 6E (c2 Bernoulli, c3 = -s^2/16) "
           f"{'PASS' if okC else 'FAIL'}; "
           f"L=3 keeps c2 AND c3 {'PASS' if okD else 'FAIL'}; "
           f"c4(s->0) = -1/720 all L (5Q) {'PASS' if okE else 'FAIL'}")
assert okA and okB and okC and okD and okE
OUT.append("   => THE BERNOULLI-BREAK RUNG IS THE LEG COUNT: the first "
           "deep coefficient the admixture touches is c_{L+1}. NEW "
           f"rung of record (L=2): c4 = {sp.sstr(sp.simplify(lad[2][c4]))}")
OUT.append("")

# ---------- G3/G4: numeric solvers, tails, nu(1), residuals
def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))
def nu_be(y):
    xv = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(xv > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(xv, 40))))

def make_ambL(s, Lc):
    sL = float(s)**Lc
    def nu_run(y):
        y = np.clip(np.asarray(y, float), 1e-14, None)
        ly = np.log(y)
        nu = nu_simple(y)
        for _ in range(80):
            d1 = 2.0*nu - 1.0
            b = 0.5*sL/d1**Lc
            db = -Lc*sL/d1**(Lc + 1)
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

S_FID = math.sqrt(g_of(n_amb_of(0.02)))
S_BIN = math.sqrt(g_of(n_amb_of(1.15)))
yv = np.logspace(-8, 4, 400)
OUT.append(f"G3/G4 numeric (fiducial galaxy gate s = {S_FID:.4f}, "
           f"binary s = {S_BIN:.4f}):")
for Lc in (1, 2, 3):
    fn = make_ambL(S_FID, Lc)
    nu = fn(yv)
    sL = S_FID**Lc
    b = 0.5*sL/(2.0*nu - 1.0)**Lc
    u = np.exp(np.minimum(0.5*(1.0+b)*np.log(yv) + b*np.log(nu), 60.0))
    n = np.where(u < 60.0, 1.0/np.expm1(np.minimum(u, 60.0)), 0.0)
    rr = float(np.max(np.abs(nu - 1.0 - n)/(1.0 + n)))
    yw = np.array([20.0, 60.0])
    nm = fn(yw) - 1.0
    p_hat = float(np.diff(np.log(-np.log(nm)))[0]/np.diff(np.log(yw))[0])
    p_pred = 0.5 + 0.25*S_FID**Lc
    okT = abs(p_hat - p_pred) < 0.03 and rr < 2e-12
    uniq = True
    for yq in (1e-6, 1e-2, 1.0, 50.0):
        nug = np.linspace(1.0 + 1e-9,
                          float(nu_simple(np.array([yq]))[0])*3, 20000)
        bq = 0.5*sL/(2.0*nug - 1.0)**Lc
        uq = np.exp(np.minimum(0.5*(1.0+bq)*np.log(yq)
                               + bq*np.log(nug), 60.0))
        nq = np.where(uq < 60.0, 1.0/np.expm1(np.minimum(uq, 60.0)), 0.0)
        uniq &= int(np.sum(np.diff(np.sign(nug - 1.0 - nq)) != 0)) == 1
    n1 = float(fn(np.array([1.0]))[0])
    OUT.append(f"  L={Lc}: resid {rr:.1e}, unique {uniq}, tail p_hat = "
               f"{p_hat:.4f} (pred {p_pred:.4f}), nu(1) = {n1:.4f} -> "
               f"{'PASS' if (okT and uniq) else 'FAIL'}")
    assert okT and uniq
OUT.append(f"  binary-side tails (e = 1.15): p(L=1) = "
           f"{0.5 + 0.25*S_BIN:.3f}, p(L=2) = {0.5 + 0.25*S_BIN**2:.3f}, "
           f"p(L=3) = {0.5 + 0.25*S_BIN**3:.3f} (measured: they hold 1/2;"
           f" L=1's 0.586 sits between the 5K-tested 0.5 and 0.65 -- the"
           f" binary side is NOT the L discriminator, the deep rung is)")
OUT.append("")

# ---------- G5: per-galaxy Chae gates + the 6I bars
emax, eno = [], []
with open('data/chae2021_table3.csv') as f:
    for row in csv.DictReader(l for l in f if not l.startswith('#')):
        emax.append(10.0**float(row['log_eN_maxclust']))
        eno.append(10.0**float(row['log_eN_noclust']))
emax, eno = np.array(emax), np.array(eno)
def gate_stats(ev):
    n = np.array([n_amb_of(e) for e in ev])
    g = np.array([g_of(v) for v in n])
    p = 0.5 + g/4
    return (float(np.median(ev)), float(np.median(g)),
            float(np.median(p)), float(p.min()), float(p.max()))
em, gm_, pm, plo, phi = gate_stats(emax)
en_, gn_, pn, pnlo, pnhi = gate_stats(eno)
OUT.append(f"G5 Chae+21 Table 3 gates ({len(emax)} galaxies):")
OUT.append(f"  maxclust: median e_N = {em:.5f} -> median g = {gm_:.3f}, "
           f"median p = {pm:.3f} (span {plo:.3f}-{phi:.3f})")
OUT.append(f"  noclust:  median e_N = {en_:.5f} -> median g = {gn_:.3f}, "
           f"median p = {pn:.3f} (span {pnlo:.3f}-{pnhi:.3f})")
OUT.append(f"  fiducial (6E/6F, e = 0.02): g = {g_of(n_amb_of(0.02)):.3f}"
           f", p = {0.5 + g_of(n_amb_of(0.02))/4:.3f}")
OUT.append("  => the MEASURED ambients are sharper than the 0.02 "
           "fiducial (both columns); the measured-ambient tail "
           "postdiction moves to p ~ 0.71-0.74 -- still inside the "
           "5G/5T band 0.65-0.75, now nearer its upper edge (stated "
           "plainly, no spin).")
OUT.append("")
OUT.append("PRE-REGISTERED 6I BARS (fixed before any 6I fit runs):")
OUT.append("  A (L contest, vertical-hardened, global fiducial s = "
           f"{S_FID:.4f}): the grammar survives iff L=2 ranks BEST of "
           "{L1, L2, L3}. L=1 carries the deep Bernoulli break (c2 = "
           f"{1.0/12.0 - S_FID/8.0:+.4f}) and the strongest tail "
           "(p = 0.717): an L=1 WIN would mean the tail out-votes the "
           "deep = grammar strike. Whether the ultra-deep arm resolves "
           "the c2 rung is exactly what the contest measures.")
OUT.append("  B (measured-ambient leg, L=2): per-galaxy Chae gates "
           "(maxclust fiducial, noclust variant; unmatched + lensing "
           "carry the matched-median gate, disclosed). Pre-register: "
           "per-galaxy >= global on both treatments (the measured "
           "ambients must not hurt); plain reaching <= -100 RESOLVES "
           "the 6F disclosed partial (the miss was the fiducial-e "
           "artifact). STRIKE if per-galaxy < global.")
OUT.append("  C (quadrupole record): amb at the solar ambient (per-eN "
           "gates) expected within ~2% of BE's q; lock join with the 6G "
           "alpha-hat = 1.060 +/- 0.024 expected ~4x Cassini -- no "
           "rescue, recorded.")

out = "\n".join(OUT)
print(out)
with open('data/stage6h_grammar.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage6h_grammar.txt")
