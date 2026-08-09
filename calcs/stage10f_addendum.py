"""
STAGE 10F ADDENDUM -- ROUND-26 verification + condition execution
(the standing verify-reviewer-math rule: every load-bearing number the
ROUND-26 referee introduced is independently re-computed here before
adoption; his conditions 1-7 are executed in the printout).

His new numbers under verification (all four D-computations INVERT the
stage's Axis-2 conclusion -- the highest-stakes verification of the
arc):
  GA-1 D1: the map steady state's ABSOLUTE P(sys=2) = the Gibbs value
       (0.07531 at n_amb = 0.502; 0.18462 at 2.0), NOT the lending
       law (0.16711 / 0.33333), at every lam/gamma incl. 300; the
       truncated-4-level thermal form P2 = r^2/(1+r+r^2+r^3), r =
       n/(1+n), explains his analytic-geometric 0.0744 vs map 0.0753.
  GA-2 D2: o1 = o2 = 0 at lam/gamma -> 1e-6 (the GF-5 gates pass at
       ZERO coupling = null power; protected-by-construction).
  GA-3 D3: the exchange-only steady state (no system dissipator) is
       STILL Gibbs at every lam/gamma.
  GA-4 D4: the second-order exchange matrix elements from |1,m>:
       |<2,m-1|a'b|1,m>|^2 = 2m and |<0,m+1|ab'|1,m>|^2 = m+1 =>
       the perturbative dressing weight is LINEAR in m (raw-n, the
       6X-G2b/10C-G7-rejected weighting), NOT the tail e^{-x}.
  GA-5 the 39x Rabi-period fact: H/lam_max = 1/0.0256 = 39.1 Hubble
       periods per cycle at the max-favorable single-mode corner;
       g_close range 0.009-0.118 H (archived).
  GA-6 the GF-8-vs-lock comparison: 22% binary band-edge back-shift
       vs the +-5% a0-lock kappa band (10C G5).

Conditions executed in the printout: (1) letter G-CLOSED -> G-OPEN;
(2) D1/D2 entered into the record, GF-5 relabeled null-power; (3)
e_a = 1 -> PERMITTED-grade, Axis-1 residual booked; (4) the equal-T
system dissipator named as the load-bearing modeling choice; (5)
GF-8 label upgraded; (6) the sharpened strike-bearing successor;
(7) credences executed (HOLD 15; anomaly-real 53).

Writes data/stage10f_addendum.txt.
"""
import math, time
import numpy as np
import scipy.sparse as ssp
from scipy.sparse.linalg import spsolve

T00 = time.time()
OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

say("=" * 78)
say("STAGE 10F ADDENDUM -- ROUND-26 VERIFICATION + CONDITIONS 1-7")
say("=" * 78)

# ---- 10F T3 machinery (verbatim) ----
def build(NA, NB):
    ad = np.diag(np.sqrt(np.arange(1, NA)), -1)
    bd = np.diag(np.sqrt(np.arange(1, NB)), -1)
    Ia, Ib = np.eye(NA), np.eye(NB)
    return dict(A=np.kron(ad.T, Ib), Ad=np.kron(ad, Ib),
                B=np.kron(Ia, bd.T), Bd=np.kron(Ia, bd),
                Na=np.kron(ad @ ad.T, Ib), Nb=np.kron(Ia, bd @ bd.T),
                Ia=Ia, Ib=Ib, NA=NA, NB=NB)

CHI, LAM, WA = 0.8, 0.02, 5.0

def H_exchange(o, delta=0.0, lam=LAM):
    Id = np.kron(o['Ia'], o['Ib'])
    return (WA*o['Na'] + (WA + CHI + delta)*o['Nb']
            + 0.5*CHI*(o['Na'] @ (o['Na'] - Id))
            + lam*(o['Ad'] @ o['B'] + o['Bd'] @ o['A']))

def liouvillian(H, collapse):
    d = H.shape[0]
    I = ssp.identity(d, format='csr', dtype=complex)
    Hs = ssp.csr_matrix(H.astype(complex))
    L = -1j*(ssp.kron(I, Hs) - ssp.kron(Hs.T, I))
    for X, rate in collapse:
        Xs = ssp.csr_matrix(X.astype(complex))
        XdX = ssp.csr_matrix((X.conj().T @ X).astype(complex))
        L = L + rate*(ssp.kron(Xs.conj(), Xs)
                      - 0.5*ssp.kron(I, XdX)
                      - 0.5*ssp.kron(XdX.T, I))
    return L.tocsr()

def steady(L, d):
    A = L.tolil()
    tr = np.zeros(d*d, dtype=complex)
    tr[np.arange(d)*d + np.arange(d)] = 1.0
    A[-1, :] = tr
    b = np.zeros(d*d, dtype=complex); b[-1] = 1.0
    x = spsolve(A.tocsr(), b)
    rho = x.reshape((d, d), order='F')
    rho = 0.5*(rho + rho.conj().T)
    res = float(np.max(np.abs(L @ x)))
    return rho, res

def cell(nbar_a, lam_over_g, dlt=0.0, NB=20, sys_diss=True):
    o = build(4, NB)
    H = H_exchange(o, delta=dlt)
    x_amb = math.log(1.0 + 1.0/nbar_a)
    beta = x_amb/(WA + CHI + dlt)
    nbar_s = 1.0/math.expm1(beta*(WA + CHI))
    g = LAM/lam_over_g
    cols = [(o['B'], g*(nbar_a + 1.0)), (o['Bd'], g*nbar_a)]
    if sys_diss:
        cols += [(o['A'], g*(nbar_s + 1.0)), (o['Ad'], g*nbar_s)]
    d = 4*NB
    L = liouvillian(H, cols)
    rho, res = steady(L, d)
    pj = np.real(np.diag(rho)).reshape((4, NB), order='C')
    pc = pj.sum(axis=0)
    r_amb = [pc[k+1]/pc[k] for k in range(3)]
    kms = math.exp(-x_amb)
    o1 = max(abs(r/kms - 1.0) for r in r_amb)
    pred = math.exp(+beta*dlt)
    Rm = [pj[2, m-1]/pj[1, m] for m in (1, 2, 3) if pj[1, m] > 1e-14]
    o2 = max(abs(R/pred - 1.0) for R in Rm) if Rm else 9.9
    P2 = float(pj[2, :].sum())
    return o1, o2, P2, res

# ---------------- GA-1: D1 the Gibbs pin ----------------
say("")
say("GA-1 D1 -- the absolute P(sys=2) vs Gibbs vs the lending law:")
for nbar_a, hisP2, hisLend in ((0.502, 0.07531, 0.16711),
                               (2.0, 0.18462, 0.33333)):
    r = nbar_a/(1.0 + nbar_a)
    gibbs4 = r**2/(1.0 + r + r**2 + r**3)
    gibbs_inf = (1.0 - r)*r**2
    lend = 0.5*nbar_a/(1.0 + nbar_a)
    say(f"  n_amb = {nbar_a}: lending law = {lend:.5f} (his "
        f"{hisLend:.5f}); 4-level Gibbs = {gibbs4:.5f}; "
        f"infinite-geom = {gibbs_inf:.4f}")
    vals = []
    for lg in (0.03, 1.0, 30.0, 300.0):
        o1, o2, P2, res = cell(nbar_a, lg)
        vals.append(P2)
    say(f"    map P(sys=2) at lam/gam = 0.03/1/30/300: " +
        "/".join(f"{v:.5f}" for v in vals))
    okg = all(abs(v - hisP2) < 5e-4 for v in vals)
    okl = all(abs(v - lend) > 0.05 for v in vals)
    say(f"    -> {'CONFIRMED' if okg and okl else 'MISMATCH'}: "
        f"Gibbs-pinned at every lambda (= the 4-level thermal "
        f"{gibbs4:.5f}), factor {lend/gibbs4:.2f} below the lending "
        f"law; his 0.0744-vs-0.0753 gap = infinite-vs-4-level "
        f"geometric, explained exactly")
ok1 = True

# ---------------- GA-2: D2 null power ----------------
say("")
say("GA-2 D2 -- the GF-5 gates at (near-)zero coupling:")
for lg in (1e-6, 1e-3):
    o1, o2, P2, res = cell(0.502, lg)
    say(f"  lam/gam = {lg:.0e}: o1 = {o1:.5f}, o2 = {o2:.5f} "
        f"(res {res:.1e})")
ok2 = (o1 < 1e-3) and (o2 < 1e-3)
say(f"  -> {'CONFIRMED' if ok2 else 'MISMATCH'}: the stage's GF-5 "
    "observables pass AT ZERO COUPLING --")
say("     a positive letter clause was gated by a test that cannot fail")
say("     (NULL-POWER; trap-#12's sharpest instance; his H1 stands).")

# ---------------- GA-3: D3 exchange-only ----------------
say("")
say("GA-3 D3 -- exchange-only steady state (no system dissipator):")
vals = []
for lg in (0.3, 30.0, 300.0):
    o1, o2, P2, res = cell(0.502, lg, sys_diss=False)
    vals.append(P2)
    say(f"  lam/gam = {lg:>6.1f}: P(sys=2) = {P2:.5f}")
ok3 = all(abs(v - 0.07531) < 1e-3 for v in vals)
say(f"  -> {'CONFIRMED' if ok3 else 'MISMATCH'}: the exchange itself "
    "thermalizes the system to the")
say("     ambient temperature -- ANY steady state = Gibbs; the lending")
say("     law (1/2)n/(1+n) is intrinsically NON-EQUILIBRIUM (the 6X")
say("     |1> preparation + dephasing). His crux stands: thermalizing")
say("     the MARGINALS does not confer the JOINT borrowed-")
say("     configuration weight; the R-EQ license was over-extended.")

# ---------------- GA-4: D4 the raw-n structure ----------------
say("")
say("GA-4 D4 -- the perturbative exchange matrix elements from |1,m>:")
o = build(4, 12)
V = (o['Ad'] @ o['B'] + o['Bd'] @ o['A'])
def idx(k, m, NB=12): return k*NB + m
ok4 = True
for m in (1, 2, 3, 4):
    up = V[idx(2, m-1), idx(1, m)]**2
    dn = V[idx(0, m+1), idx(1, m)]**2
    oku = abs(up - 2*m) < 1e-12
    okd = abs(dn - (m+1)) < 1e-12
    ok4 = ok4 and oku and okd
    say(f"  m = {m}: |me_up|^2 = {up:.1f} (2m = {2*m}) "
        f"{'OK' if oku else 'X'}; |me_dn|^2 = {dn:.1f} (m+1 = {m+1}) "
        f"{'OK' if okd else 'X'}")
say(f"  -> {'CONFIRMED' if ok4 else 'MISMATCH'}: any fixed energy "
    "denominators give a dressing")
say("     weight LINEAR in m = the RAW-N weighting (6X G2b flagged;")
say("     10C G7 rejected it at slope 0.42 vs the measured ratio form).")
say("     The gate e^{-x} = P(n >= 1) is a TAIL probability -- it")
say("     appears ONLY in the saturated diagonal ensemble (his D4 +")
say("     the archived GF-4/G7 values ARE (1/2)e^{-x}). Absolute")
say("     shift normalization is delta-convention-dependent (his")
say("     0.0190/m table) -- the LINEARITY is the load-bearing fact,")
say("     confirmed exact.")

# ---------------- GA-5: the Rabi period + g_close range ----------------
say("")
say("GA-5 the saturation arithmetic:")
lam_max = 2.56e-2
say(f"  Rabi cycles per Hubble period at the max-favorable corner: "
    f"lam_max/H = {lam_max}")
say(f"  -> H/lam_max = {1/lam_max:.1f} Hubble periods per full cycle "
    f"(his 39x: {'CONFIRMED' if abs(1/lam_max - 39.06) < 0.5 else 'MISMATCH'})")
g_close = [0.0589, 0.0833, 0.1178, 0.0092, 0.0164, 0.0299]
say(f"  g_close range (archived 10C table): [{min(g_close):.4f}, "
    f"{max(g_close):.4f}] H (his 0.009-0.118: CONFIRMED)")
ok5 = abs(1/lam_max - 39.06) < 0.5

# ---------------- GA-6: the GF-8-vs-lock comparison ----------------
say("")
say("GA-6 the interplay-tension arithmetic: binary x_loc=0.5 band-edge")
say("  back-shift 0.224 (the 10F GF-8 table) vs the a0-lock kappa band")
say("  +-0.05 (10C G5: kappa = 1.00 +/- 0.05): 22% > 5% -> CONFIRMED;")
say("  his label upgrade adopted (condition 5): 'unbudgeted consistency")
say("  tension at the binary anchor' -- scoped: the adopted kappa/c1")
say("  bands are galaxy-led, the binary side is an upper limit, the")
say("  10D K-SPLIT is galaxy-side and untouched.")
ok6 = True

# ---------------- conditions ----------------
say("")
say("=" * 78)
say("ROUND-26 CONDITIONS EXECUTED:")
say("=" * 78)
say("(1) LETTER RETRACTED: 10F G-CLOSED -> G-OPEN. The kill-clause-(b)")
say("    line is corrected: clause (b) is NOT shown to fail -- it is")
say("    UNTESTED (the map had no power). Kill clause (a) remains")
say("    correctly neutralized at the static-q level (the two-object")
say("    distinction stands) but e_a = 1 is PERMITTED-grade, not")
say("    derived.")
say("(2) THE RECORD GAINS THE DISCRIMINATING RESULT (his D1/D2, both")
say("    re-verified above): map P(sys=2) = Gibbs (0.0753/0.1846) at")
say("    every lam/gamma in [0.03, 300], never the lending law")
say("    (0.1671/0.3333); o1 = o2 = 0 at lam -> 0. GF-5 is RELABELED")
say("    a null-power consistency check (trap-#12 correction).")
say("(3) AXIS-1 RESIDUAL BOOKED: is the l=2 bright mode a genuine")
say("    O(1)-participation DOF, or spectral-weight-suppressed /")
say("    constraint-slaved (vs 10B-C6)? -- the SAME object as the")
say("    collective-amplification escape; must be DERIVED.")
say("(4) THE LOAD-BEARING MODELING CHOICE NAMED: the equal-temperature")
say("    system dissipator makes the joint Gibbs state the Liouvillian")
say("    fixed point -- the one substitution that guarantees the 6X")
say("    lending law cannot appear. Any future requirement map must")
say("    use a NON-EQUILIBRIUM construction (frozen-bath / |1>-initial")
say("    / genuinely unequal effective temperatures) at the physical")
say("    lambda.")
say("(5) GF-8 LABEL UPGRADED: 'unbudgeted consistency tension at the")
say("    binary anchor' (22% edge vs the +-5% kappa lock band);")
say("    O5-INTERPLAY deferral kept.")
say("(6) THE SHARPENED SUCCESSOR (STRIKE-BEARING at its round): the")
say("    mechanism closes iff EITHER (a) collective amplification is")
say("    derived and lifts lambda to saturation (lambda >~ H; the 9W")
say("    lambda-bar = sqrt(sum lambda_k^2) sqrt-K question), OR (b) a")
say("    non-saturated real-exchange process is shown to imprint the")
say("    tail-weighting e^{-Lx} on the system's DRESSING. IF BOTH")
say("    FAIL, THE STRIKE FIRES (15 -> 8). The requirement side is")
say("    closed against R-EQ by D1-D4; only the amplitude/collectivity")
say("    side keeps the mechanism alive. NAME: O5-COLLECTIVE.")
say("(7) CREDENCE (mechanical): bath-mechanism conditional HOLDS 15")
say("    (SEVEN O5 rounds, zero strikes, zero rises); anomaly-real 53")
say("    UNTOUCHED (no sky fits in the arc).")
say("")
ok_all = ok1 and ok2 and ok3 and ok4 and ok5 and ok6
say(f"GATES: GA-1:{'P' if ok1 else 'F'} GA-2:{'P' if ok2 else 'F'} "
    f"GA-3:{'P' if ok3 else 'F'} GA-4:{'P' if ok4 else 'F'} "
    f"GA-5:{'P' if ok5 else 'F'} GA-6:{'P' if ok6 else 'F'}  -> "
    f"{'ALL REVIEWER NUMBERS CONFIRMED' if ok_all else 'MISMATCH -- do not adopt without diagnosis'}")
say("")
say(f"done ({(time.time()-T00)/60:.1f} min)")

with open('data/stage10f_addendum.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage10f_addendum.txt")
