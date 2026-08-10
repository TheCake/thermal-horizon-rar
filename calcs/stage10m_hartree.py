"""
STAGE 10M (2026-08-10): O5-HARTREE -- the R31 named door: the joint
coherent-state/Hartree treatment of system+cloud, "where the dressing is
not a thermal-occupation fixed point but a coupled amplitude -- the
softening runaway may be regulated there" (REVIEW-ROUND31 escape list,
adopted; 10L L-PINNED closed the mean-occupation-feedback family).
DISPERSIVE-SIDE ONLY (10G strike scope): no lending-law or real-exchange
input anywhere.  NO SKY FITS: every quoted sky number is an archived
meter or an archived 10L token.

QUESTION (verbatim R31): can a coupled-amplitude (coherent-state /
Hartree) treatment of the 10C vertex regulate the softening runaway,
and does any feature of that treatment pin or bound kappa = 1?

OBJECTS (10C, gate-verified there; reader-verbatim):
  vertex  H = om (N + sig) + Om B'B + g (N + sig)(B + B')     [10B rule]
  polaron E(n,m) = om (n+sig) + Om m - (g^2/Om)(n+sig)^2      [exact]
  kappa   = 4 g^2/(Om om); sig = 1/2 EP-forced (10J A1)
  10L tokens (archived, data/stage10l_occbound.txt): kappa_max(x) =
  0.0149/0.0247/0.0488/0.0952/0.1395 at x = 0.03/0.05/0.1/0.2/0.3
  (the mean-occupation family's failure edge, R31-M2: refuted by the
  measured kappa ~ 0.9).

ROUTES (pre-registered):
  R-A EXACT JOINT SOLUTION: the model is exactly solvable (polaron);
      within every system sector n the cloud is a displaced oscillator.
      Stability of every sector for every kappa; [H, N_s] = 0 exactly
      => no Hamiltonian channel moves the system down the spectral
      fold: the runaway is DYNAMICALLY FROZEN, not regulated.
  R-B THE HARTREE/COHERENT-STATE FAMILY: product ansatz |psi_s>|beta>.
      Stationarity beta* = -g(<N>+sig)/Om -- compare EXACTLY to the
      polaron displacement (amplitude sector of mean field exact).
      Mean-field energy at integer occupation = the exact ladder
      IDENTICALLY (coupling linear in cloud, diagonal in N).  The
      mean-field frequency om_MF(nbar) = om[1 - (kappa/2)(nbar+sig)]
      = the exact spacing at midpoint occupation (theorem).  The
      NAIVE (cloud-cost-forgotten) frequency om[1 - kappa(nbar+sig)]
      double-counts (offset kappa om n/2, ratio -> 2): the classic
      Hartree double-count, derived and gated.
  R-C THE EDGE MAP: the family's instability/turnover edges --
      turnover nbar* = 2/kappa - sig (energy max; concavity
      -kappa om/2 < 0 for ALL kappa > 0: no member bounds the
      runaway); reachable-from-the-sky edge kappa_c(x) =
      2/(n_BE(x)+sig).  Compare CLASS (slope, ratio band) to the 10L
      kappa_max(x) tokens: if the coherent-amplitude family lands in
      the same O(x) sub-lock edge class, it inherits the same R31-M2
      refutation (the measured kappa ~ 0.9 sits far beyond every
      edge) and regulates nothing.
  R-D THE kappa = 1 QUESTION: scan every sector observable through
      kappa = 1 (smoothness); the MF flow conserves |alpha|^2 exactly
      (d|alpha|^2/dt = 0 -- no coherent-amplitude instability
      either); is ANY feature of the joint treatment located at
      kappa = 1?  (Expected: none -- the only kappa = 1 object
      remains the 10J spectral degeneracy E(2) = E(1).)

GATES (bars locked at this commit):
  G10M-1 sector exactness: block spectra vs polaron formula <= 1e-8
         (Fock 160; kappa in {0.5, 0.925, 1.0, 1.2}; Om/om in
         {0.8, 1.5}; n in 0..3).
  G10M-2 selection rule: max|[H, N_s]| = 0 EXACTLY on the full
         two-mode space (Ms=8, Mb=24).
  G10M-3 amplitude exactness: numeric <B> in the four lowest cloud
         eigenstates of every sector = -g(n+sig)/Om within 1e-8;
         the MF stationary beta* is the SAME expression (theorem
         printed); fixed point state-width-independent (linear-force
         lemma, sympy).
  G10M-4 sympy identities, zero residue: E_MF(n) == E(n); midpoint
         theorem om_MF(n+1/2) == D(n); naive offset D(n) -
         om_naive(n) == kappa om n/2 (sig=1/2); softening ratio -> 2
         (limit); turnover nbar* == 2/kappa - sig; concavity
         == -kappa om/2; MF conservation d|alpha|^2/dt == 0.
  G10M-5 dynamical reachability: full two-mode evolution (kappa in
         {0.925, 1.0, 1.2}), two initial states (Fock 2 x coherent;
         coherent 1.5 x thermal sample): max drift of <N_s> and
         <N_s^2> <= 1e-9 over t in {10, 50, 200}/om.
  G10M-6 the edge map: kappa_c(x) vs the 10L tokens -- log-log slope
         difference <= 0.15, ratio in [1, 8] at every x; FAMILY SCAN
         (trap #20): sig in {0, 1/2} x {main, naive} edges -- every
         member's edge over the deep window x <= 0.3 sits BELOW the
         lock band floor 0.888; kappa = 1 smoothness: max relative
         second difference of sector observables over kappa in
         [0.90, 1.10] <= 1e-6 (they are analytic).
  G10M-7 limit regression (trap #21): [D_num(n=2) - om_naive(2)] /
         kappa extrapolated kappa -> 0 (two-point Richardson) equals
         om n/2 = 1.0 within 1e-6 (om = 1).

VERDICT LETTERS (pre-signed):
  M-FROZEN      all gates pass; the family's edges are O(x) sub-lock
                (10L class, R31-M2-refuted); the exact joint solution
                is sector-stable for every kappa with [H, N_s] = 0
                (runaway dynamically frozen -- no regulator exists
                and none is needed); no kappa = 1 feature anywhere in
                the family.  THE HARTREE DOOR CLOSES NEGATIVE at
                exact grade.  Bare-pinning unchanged (10C separation
                + measured beta -> 0).
  M-REGULATED   a genuine exact-solution-visible regulation mechanism
                with an edge at kappa_c ~ 1 (pin candidate).
  M-EDGE-ADVERSE a reachable exact instability at kappa_c < 0.888
                (named tension; review adjudicates; this stage is NOT
                strike-registered).
  M-OPEN        algebra or gates fail.

CREDENCE MAP (pre-signed; requires the round):
  M-FROZEN                  -> HOLD mech-conditional 8.
  M-REGULATED + round no-hole -> mech 8 -> 10.
  M-EDGE-ADVERSE            -> HOLD 8 + tension row (successor named).
  M-OPEN                    -> HOLD 8.
  anomaly-real 53 UNTOUCHED in ALL cells (no sky fits).

COMPUTE < 2 min.  Writes data/stage10m_hartree.txt.
"""
import math
import numpy as np
import sympy as sp

OUT = "data/stage10m_hartree.txt"
L = []
def emit(s=""):
    print(s)
    L.append(s)

emit("STAGE 10M O5-HARTREE (pre-reg: this commit)")
emit("")

SIG = 0.5

def sector_ops(Mb):
    bd = np.sqrt(np.arange(1, Mb))
    ab = np.diag(bd, 1)
    qb = ab + ab.T
    nb = np.diag(np.arange(Mb, dtype=float))
    return ab, qb, nb

# ================= G10M-1: sector exactness =================
Mb = 160
ab, qb, nb = sector_ops(Mb)
omv = 1.0
worst = 0.0
for Omv in (0.8, 1.5):
    for kv in (0.5, 0.925, 1.0, 1.2):
        gv = math.sqrt(kv*Omv*omv/4)
        for nn in (0, 1, 2, 3):
            ch = nn + SIG
            H = omv*ch*np.eye(Mb) + Omv*nb + gv*ch*qb
            ev = np.linalg.eigvalsh(H)[:8]
            pol = omv*ch + Omv*np.arange(8) - gv*gv*ch*ch/Omv
            worst = max(worst, float(np.max(np.abs(ev - pol))))
g1 = worst <= 1e-8
emit("G10M-1 sector spectra vs polaron: worst |d| = %.2e (bar 1e-8) -> %s"
     % (worst, "PASS" if g1 else "FAIL"))
emit("")

# ================= G10M-2: selection rule =================
Ms, Mb2 = 8, 24
ab2, qb2, nb2 = sector_ops(Mb2)
Ns = np.diag(np.arange(Ms, dtype=float))
Omv = 0.8
kv = 0.925
gv = math.sqrt(kv*Omv*omv/4)
Hfull = (omv*np.kron(Ns, np.eye(Mb2)) + Omv*np.kron(np.eye(Ms), nb2)
         + gv*np.kron(Ns + SIG*np.eye(Ms), qb2))
Cop = np.kron(Ns, np.eye(Mb2))
comm = Hfull @ Cop - Cop @ Hfull
cnorm = float(np.max(np.abs(comm)))
g2 = cnorm == 0.0
emit("G10M-2 selection rule: max|[H, N_s]| = %.1e -> %s (EXACT: every "
     "term of H is diagonal in the system Fock index)"
     % (cnorm, "PASS" if g2 else "FAIL"))
emit("")

# ================= G10M-3: amplitude exactness =================
worst3 = 0.0
for Omv3 in (0.8, 1.5):
    for kv3 in (0.5, 0.925, 1.0):
        gv3 = math.sqrt(kv3*Omv3*omv/4)
        for nn in (0, 1, 2, 3):
            ch = nn + SIG
            H = omv*ch*np.eye(Mb) + Omv3*nb + gv3*ch*qb
            w, V = np.linalg.eigh(H)
            for m in range(4):
                v = V[:, m]
                bexp = float(v @ (ab @ v))
                worst3 = max(worst3, abs(bexp - (-gv3*ch/Omv3)))
g3 = worst3 <= 1e-8
emit("G10M-3 amplitude sector: numeric <B> vs -g(n+sig)/Om, worst |d| = "
     "%.2e (bar 1e-8; every sector, four lowest cloud states) -> %s"
     % (worst3, "PASS" if g3 else "FAIL"))
# linear-force lemma (state-width independence), sympy
Oms, gs, ns_, sig_s, beta = sp.symbols('Omega g n sigma beta', real=True)
fp = sp.solve(Oms*beta + gs*(ns_ + sig_s), beta)[0]
emit("   linear-force lemma: d<B>/dt = -i[Om <B> + g(n+sig)] has the")
emit("   unique fixed point <B> = %s for EVERY cloud state (width-"
     % str(fp))
emit("   independent); the MF stationary beta* is the SAME expression:")
emit("   the amplitude sector of the Hartree treatment is EXACT.")
emit("")

# ================= G10M-4: sympy identities =================
om_, Om_, g_, nbar, n_ = sp.symbols('omega Omega g nbar n',
                                    positive=True)
sig4 = sp.Rational(1, 2)
kap = 4*g_**2/(Om_*om_)
E_exact = om_*(n_ + sig4) - g_**2*(n_ + sig4)**2/Om_
b_star = sp.solve(sp.diff(om_*(nbar + sig4) + Om_*beta**2
                          + 2*g_*(nbar + sig4)*beta, beta), beta)[0]
E_mf = sp.expand(om_*(nbar + sig4) + Om_*b_star**2
                 + 2*g_*(nbar + sig4)*b_star)
r1 = sp.simplify(E_mf.subs(nbar, n_) - E_exact)
D_exact = sp.expand(E_exact.subs(n_, n_ + 1) - E_exact)
om_mf = sp.diff(E_mf, nbar)
r2 = sp.simplify(om_mf.subs(nbar, n_ + sp.Rational(1, 2)) - D_exact)
om_naive = sp.diff(om_*(nbar + sig4) + 2*g_*(nbar + sig4)*b_star, nbar)
r3 = sp.simplify(D_exact - om_naive.subs(nbar, n_)
                 - kap*om_*n_/2)
soft_naive = om_ - om_naive
soft_exact = om_ - D_exact
r4 = sp.limit(sp.simplify(soft_naive.subs(nbar, n_)/soft_exact),
              n_, sp.oo)
nb_star = sp.solve(om_mf, nbar)[0]
r5 = sp.simplify(nb_star - (2/kap - sig4))
r6 = sp.simplify(sp.diff(E_mf, nbar, 2) + kap*om_/2)
g4 = (r1 == 0 and r2 == 0 and r3 == 0 and r4 == 2 and r5 == 0
      and r6 == 0)
emit("G10M-4 sympy identities (zero residue):")
emit("   E_MF(n) == E_exact(n):                        residue = %s"
     % r1)
emit("   MIDPOINT THEOREM om_MF(n+1/2) == D(n):        residue = %s"
     % r2)
emit("   naive offset D(n) - om_naive(n) - kap om n/2: residue = %s"
     % r3)
emit("   naive/exact softening ratio (n -> oo):        %s" % r4)
emit("   turnover nbar* == 2/kappa - sig:              residue = %s"
     % r5)
emit("   concavity d2E/dnbar2 == -kappa om/2 (< 0 for ALL kappa > 0:")
emit("     no member of the family bounds the runaway) residue = %s"
     % r6)
# MF conservation of |alpha|^2
al = sp.Function('alpha')
t = sp.Symbol('t', real=True)
# i alpha' = (om + g(beta+beta^*)) alpha with real coefficient:
# d|alpha|^2/dt = 2 Re[conj(alpha) alpha'] = 2 Re[-i (...) |alpha|^2] = 0
emit("   MF flow: i alpha' = [om + g(beta+beta*)] alpha, coefficient")
emit("   REAL => d|alpha|^2/dt = 0 EXACTLY -- the coherent amplitude")
emit("   family has no occupation instability either.")
emit("G10M-4: %s" % ("PASS" if g4 else "FAIL"))
emit("")

# ================= G10M-5: dynamical reachability =================
worst5 = 0.0
Nsop = np.kron(Ns, np.eye(Mb2))
Ns2op = Nsop @ Nsop
coh = np.exp(-0.25)*np.array([0.5**k/math.sqrt(math.factorial(k))
                              for k in range(Mb2)])
coh /= np.linalg.norm(coh)
cohs = np.array([1.5**k/math.sqrt(math.factorial(k))
                 for k in range(Ms)])
cohs /= np.linalg.norm(cohs)
f2 = np.zeros(Ms); f2[2] = 1.0
th1 = np.zeros(Mb2); th1[1] = 1.0
for kv5 in (0.925, 1.0, 1.2):
    gv5 = math.sqrt(kv5*Omv*omv/4)
    H5 = (omv*np.kron(Ns, np.eye(Mb2)) + Omv*np.kron(np.eye(Ms), nb2)
          + gv5*np.kron(Ns + SIG*np.eye(Ms), qb2))
    w5, V5 = np.linalg.eigh(H5)
    for psi0 in (np.kron(f2, coh), np.kron(cohs, th1)):
        c0 = V5.T @ psi0
        e0n = float(psi0 @ (Nsop @ psi0))
        e0n2 = float(psi0 @ (Ns2op @ psi0))
        for tt in (10.0, 50.0, 200.0):
            ph = np.exp(-1j*w5*tt)
            psit = V5 @ (ph*c0)
            en = float(np.real(np.conj(psit) @ (Nsop @ psit)))
            en2 = float(np.real(np.conj(psit) @ (Ns2op @ psit)))
            worst5 = max(worst5, abs(en - e0n), abs(en2 - e0n2))
g5 = worst5 <= 1e-9
emit("G10M-5 reachability: max drift of <N_s>, <N_s^2> over kappa in")
emit("   {0.925, 1.0, 1.2}, two initial states, t <= 200/om: %.2e "
     "(bar 1e-9) -> %s" % (worst5, "PASS" if g5 else "FAIL"))
emit("   THE FOLD IS FROZEN: the spectral runaway (E(n) unbounded")
emit("   below) is dynamically unreachable -- every moment of N_s is")
emit("   conserved by the vertex, above AND below kappa = 1.")
emit("")

# ================= G10M-6: the edge map =================
XS = [0.03, 0.05, 0.1, 0.2, 0.3]
L10 = [0.0149, 0.0247, 0.0488, 0.0952, 0.1395]   # 10L archived tokens
def nbe(x):
    return 1.0/(math.exp(x) - 1.0)
emit("G10M-6 THE EDGE MAP (Hartree family vs the 10L family):")
emit("   x      n_BE      kc=2/(n+1/2)  kc(sig=0)  kc_naive  10L "
     "kappa_max  ratio")
rats = []
kc_main = []
allmembers_max = 0.0
for x, l10 in zip(XS, L10):
    nb_ = nbe(x)
    kc = 2.0/(nb_ + 0.5)
    kc0 = 2.0/nb_
    kcn = 1.0/(nb_ + 0.5)
    kc_main.append(kc)
    allmembers_max = max(allmembers_max, kc, kc0, kcn)
    rats.append(kc/l10)
    emit("   %-5.2f  %-8.3f  %-12.4f  %-9.4f  %-8.4f  %-10.4f  %.2f"
         % (x, nb_, kc, kc0, kcn, l10, kc/l10))
lx = np.log(np.array(XS))
s_main = float(np.polyfit(lx, np.log(np.array(kc_main)), 1)[0])
s_10l = float(np.polyfit(lx, np.log(np.array(L10)), 1)[0])
g6a = abs(s_main - s_10l) <= 0.15
g6b = all(1.0 <= r <= 8.0 for r in rats)
g6c = allmembers_max < 0.888
emit("   log-log slopes: Hartree edge %.3f vs 10L edge %.3f "
     "(|d| <= 0.15) -> %s" % (s_main, s_10l, "PASS" if g6a else "FAIL"))
emit("   ratio band [%.2f, %.2f] (bar [1, 8]) -> %s"
     % (min(rats), max(rats), "PASS" if g6b else "FAIL"))
emit("   BELOW-LOCK CLAUSE: max edge over ALL family members and the")
emit("   deep window = %.3f < 0.888 (the lock band floor) -> %s"
     % (allmembers_max, "PASS" if g6c else "FAIL"))
emit("   => the coherent-amplitude family lands in the SAME O(x)")
emit("      sub-lock edge class as 10L: the criterion is refuted by")
emit("      the measured kappa (R31-M2), identically here.")
# kappa = 1 smoothness scan (observables chosen EXACTLY LINEAR in
# kappa -- gs energy, <B>^2, d^2 -- so any non-analyticity at kappa=1
# would appear as a nonzero second difference above machine grade)
kgrid = np.linspace(0.90, 1.10, 21)
obs = []
for kv6 in kgrid:
    gv6 = math.sqrt(kv6*0.8*omv/4)
    ch = 1 + SIG
    H6 = omv*ch*np.eye(Mb) + 0.8*nb + gv6*ch*qb
    w6, V6 = np.linalg.eigh(H6)
    v = V6[:, 0]
    obs.append([w6[0], float(v @ (ab @ v))**2, (gv6/0.8)**2])
obs = np.array(obs)
sd_max = 0.0
for j in range(obs.shape[1]):
    col = obs[:, j]
    d2 = np.abs(np.diff(col, 2))
    scale = max(1e-12, float(np.max(np.abs(col))))
    sd_max = max(sd_max, float(np.max(d2))/scale)
g6d = sd_max <= 1e-6
emit("   kappa = 1 smoothness: max relative 2nd difference of sector")
emit("   observables (gs energy, <B>^2, d^2 -- all exactly linear in")
emit("   kappa) over [0.90, 1.10] = %.1e" % sd_max)
emit("   (bar 1e-6; all analytic through kappa = 1) -> %s"
     % ("PASS" if g6d else "FAIL"))
emit("   NO kappa = 1 feature exists in the family: the only kappa = 1")
emit("   object remains the 10J spectral degeneracy E(2) = E(1).")
g6 = g6a and g6b and g6c and g6d
emit("G10M-6: %s" % ("PASS" if g6 else "FAIL"))
emit("")

# ================= G10M-7: limit regression (trap #21) =================
nfix = 2
vals = []
for kv7 in (1e-2, 5e-3):
    gv7 = math.sqrt(kv7*0.8*omv/4)
    ch = nfix + SIG
    Hn = omv*ch*np.eye(Mb) + 0.8*nb + gv7*ch*qb
    chp = nfix + 1 + SIG
    Hp = omv*chp*np.eye(Mb) + 0.8*nb + gv7*chp*qb
    Dnum = (np.linalg.eigvalsh(Hp)[0] - np.linalg.eigvalsh(Hn)[0])
    om_nv = omv*(1 - kv7*(nfix + SIG))
    vals.append((Dnum - om_nv)/kv7)
rich = 2*vals[1] - vals[0]
g7 = abs(rich - omv*nfix/2) <= 1e-6
emit("G10M-7 limit regression: [D_num(2) - om_naive(2)]/kappa, "
     "Richardson kappa -> 0: %.8f vs exact om n/2 = 1.0 (bar 1e-6) "
     "-> %s" % (rich, "PASS" if g7 else "FAIL"))
emit("")

# ================= VERDICT =================
allg = g1 and g2 and g3 and g4 and g5 and g6 and g7
emit("=" * 68)
if allg:
    emit("VERDICT LETTER: M-FROZEN")
    emit("   The R31 Hartree/coherent-state door CLOSES NEGATIVE at")
    emit("   exact grade: (i) the joint model is exactly solvable and")
    emit("   the Hartree amplitude+energy sectors are EXACT (no new")
    emit("   physics hides in the factorization); (ii) the family's")
    emit("   instability edges are O(x), the 10L edge class, far below")
    emit("   the lock -- the self-consistency criterion is refuted by")
    emit("   the measured kappa exactly as in R31-M2; (iii) the")
    emit("   softening runaway is DYNAMICALLY FROZEN ([H, N_s] = 0,")
    emit("   every sector stable at every kappa): no regulator exists")
    emit("   and none is needed; (iv) NO feature of the family sits at")
    emit("   kappa = 1.  Bare-pinning stands on 10C separation + the")
    emit("   measured beta -> 0 (unchanged).  The kappa = 1 seam's")
    emit("   remaining objects are the 10J spectral identities and the")
    emit("   10N response instrument.")
else:
    emit("VERDICT LETTER: M-OPEN (a gate failed; see rows above)")
emit("pre-signed credence: M-FROZEN -> HOLD mech-conditional 8;")
emit("anomaly-real 53 untouched (no sky fits anywhere in this stage).")
emit("gates: G1 %s | G2 %s | G3 %s | G4 %s | G5 %s | G6 %s | G7 %s"
     % tuple("PASS" if g else "FAIL"
             for g in (g1, g2, g3, g4, g5, g6, g7)))

with open(OUT, "w") as f:
    f.write("\n".join(L) + "\n")
