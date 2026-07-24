"""
STAGE 6M (= 6K-v2, the structured-bath leg of the desktop analog).
PRE-REGISTRATION COMMITTED BEFORE EXECUTION.

WHAT 6K-v1 LEFT OPEN: v1 put the two endpoints in BY HAND as two
competing jump channels (flat "source" channel vs dressed-resolving
channel) and found the mixing is rate-ratio-weighted — not
share-gated. But standard theory contains a PHYSICAL interpolation
between the same endpoints with NO free rate ratio: a single thermal
bath of finite spectral resolution. If the bath's bandwidth b cannot
resolve the mode's self-shift (b >> K*n: the singular-coupling /
white limit), it populates the mode at the bare (source) frequency —
the beta=0 endpoint; if it resolves the dressed comb (b << K: the
Davies / weak-coupling limit), it populates at the dressed spacings —
Gibbs of the full Hamiltonian, the beta=1 endpoint. The interpolation
parameter is SPECTRAL (b vs the Kerr splitting), and the overall rate
kappa cancels from the steady state exactly — so this realization is
automatically rate-strength-free, removing v1's kappa objection. The
question v2 asks: does the resulting admixture organize by the
grammar's OCCUPATION SHARES (tanh^2(x/2) local structure), or by
RESOLUTION (b/K and the occupation spread)?

THE MODEL (exact birth-death; golden-rule with a smeared kernel):
Kerr ladder E_n = n*w0 + (K/2)n(n-1), link spacing D_n = w0 + K*n;
one thermal bath at temperature T (x0 = w0/T sets the source
occupation n0 = n_BE(x0)); Ohmic-weighted Gaussian resolution kernel
of width b: the bath occupation entering link n is

    n_eff(n) = Int G_b(w - D_n) * w * n_BE(w/T) dw
             / Int G_b(w - D_n) * w dw     (w > 0)

(the Ohmic factor w is the standard IR-regular spectral density; the
per-link Einstein pair up = n_eff*(n+1), down = (n_eff+1)*(n+1)
preserves the pair structure, and kappa cancels in the NESS product).
ENDPOD GATES: b -> 0 gives n_eff = n_BE(D_n) and the NESS must equal
the Davies-Gibbs of the Kerr mode (1e-6); K*nmax << b << w0 gives
n_eff ~ n_BE(w0-grade) and the NESS must approach the source value n0
(2% at the widest gate point); truncation doubling; kernel
normalization 1e-12.

ESTIMATOR (fixed): lam = (nbar - n0)/(n1 - n0), n0 = n_BE(x0),
n1 = Davies-Gibbs mean (exact partition sum). Scan: x0 in
{0.2, 0.5, 1.0, 2.0, 3.5} (occupied/classical -> sparse/quantum),
K/w0 in {0.005, 0.02}, r = b/K log-spaced (b capped at 0.4*w0).

PRE-REGISTERED QUESTIONS AND BANDS:
 (i) SHARE TEST: at fixed r inside the transition window
     (0.2 < lam < 0.8), the grammar's local factor predicts
     lam ratios across x0 tracking tanh^2(x0/2) (ratio band +-30%).
 (ii) RESOLUTION TEST: lam collapses onto a function of the
     resolution variable r_eff = b / (K * spread), spread =
     sqrt(nbar(nbar+1)) (collapse band: transition midpoints within
     a factor 2 across the x0 grid when expressed in r_eff).
 (iii) The pull-statistics lemma (exact, for the record): the
     emission-comb centroid pull of a diagonal state is
     K*(<n^2>-<n>)/<n> (thermal: 2K*nbar, gate 1e-10) — the state's
     STATISTICS, not a share; the vacuum share tanh(x/2) lives in the
     (n+1/2) pull bookkeeping, not in the population dynamics.
OUTCOME TREE (pre-committed):
  PASS-SHARE: (i) holds somewhere -> the structured bath carries the
     grammar's local gating; mechanism re-grounded (credence partial
     restore toward ~25%).
  ALT-RES: (ii) holds and (i) fails -> standard structured-bath
     physics interpolates by RESOLUTION, not shares. Then the gravity
     translation is fixed and parameter-poor: the self-shift over
     bath width is rho(y) ~ x*(sqrt(nu)-1)/zeta (zeta = b_grav/T =
     O(1), assumption flagged), which VANISHES both deep (sqrt(x))
     and in the tail (x*e^-x/2) and peaks ~0.26 at x ~ 1. The sky
     (5T) demands beta = 0 deep AND transition with beta ~ 1/2-3/4
     in the TAIL alone — the resolution shape has the WRONG SIGN in
     the tail for every zeta. Pre-commitment: ALT-RES + that tail
     mismatch = the second standard realization class fails the sky;
     bath-microphysics conditional ~15% -> ~8-10%; survivors = the
     non-Markovian/horizon-specific corner (with the noted irony
     that the dS bath width ~ T makes zeta O(1) natural) and the MI
     trajectory reading.
  AMBIG: neither band -> kernel-detail sensitivity; a pseudomode
     v3 would be required before any lab claim.
Writes data/stage6m_analog2.txt.
"""
import math
import numpy as np

W0 = 1.0
L = ["STAGE 6M (6K-v2): structured-bath admixture -- shares vs "
     "resolution, pre-registered", ""]

def n_be(x):
    return 1.0/math.expm1(x) if x < 700 else 0.0

# ---------- kernel machinery (Ohmic-weighted Gaussian, w > 0)
WGRID = np.concatenate([np.linspace(1e-6, 0.2, 400),
                        np.linspace(0.2, 6.0, 4000)])
def n_eff_links(T, K, b, N):
    D = W0 + K*np.arange(N)
    G = np.exp(-0.5*((WGRID[None, :] - D[:, None])/b)**2)
    wgt = G*WGRID[None, :]
    nbe = 1.0/np.expm1(np.clip(WGRID/T, 1e-9, 700.0))
    num = np.trapz(wgt*nbe[None, :], WGRID, axis=1)
    den = np.trapz(wgt, WGRID, axis=1)
    return num/den

def ness_nbar(T, K, b, N):
    ne = n_eff_links(T, K, b, N)
    lr = np.log(ne) - np.log(ne + 1.0)
    lp = np.concatenate([[0.0], np.cumsum(lr)])
    lp -= lp.max()
    p = np.exp(lp); p /= p.sum()
    return float(np.sum(np.arange(N + 1)*p)), p

def gibbs_mean(T, K, N):
    n = np.arange(N + 1)
    lp = -(n*W0 + 0.5*K*n*(n - 1.0))/T
    lp -= lp.max()
    p = np.exp(lp); p /= p.sum()
    return float(np.sum(n*p))

# ---------- gates
N = 400
T = 1.0/0.5                       # x0 = 0.5
K = 0.02
g_davies = ness_nbar(T, K, 1e-4*K, N)[0]
g_ref = gibbs_mean(T, K, N)
g1 = abs(g_davies - g_ref)
L.append(f"G1 b->0 (Davies): nbar = {g_davies:.8f} vs Gibbs {g_ref:.8f} "
         f"-> d = {g1:.1e} {'PASS' if g1 < 1e-6 else 'FAIL'}")
assert g1 < 1e-6
n0_ref = n_be(0.5)
g_wide = ness_nbar(T, 0.0005, 0.4*W0, N)[0]   # K*nmax tiny vs b
g2 = abs(g_wide - n0_ref)/n0_ref
L.append(f"G2 wide-b (source): nbar = {g_wide:.6f} vs n0 = {n0_ref:.6f} "
         f"-> rel d = {g2:.2%} {'PASS' if g2 < 0.02 else 'FAIL'}")
assert g2 < 0.02
a1 = ness_nbar(T, K, 0.05, N)[0]
a2 = ness_nbar(T, K, 0.05, 2*N)[0]
g3 = abs(a1 - a2)/a2
L.append(f"G3 truncation doubling: rel d = {g3:.1e} "
         f"{'PASS' if g3 < 1e-8 else 'FAIL'}")
assert g3 < 1e-8
# pull-statistics lemma gate: thermal comb centroid = 2K nbar
nb_t = 1.7
pth = (nb_t/(1 + nb_t))**np.arange(N + 1)
pth /= pth.sum()
nn = np.arange(N + 1)
m1, m2 = np.sum(nn*pth), np.sum(nn*nn*pth)
pull = (m2 - m1)/m1
g4 = abs(pull - 2*nb_t)/(2*nb_t)
L.append(f"G4 pull lemma (thermal): (<n^2>-<n>)/<n> = {pull:.6f} vs "
         f"2*nbar = {2*nb_t:.6f} -> {'PASS' if g4 < 1e-6 else 'FAIL'} "
         f"(the coherent pull weighs the state's STATISTICS; the "
         f"vacuum share tanh(x/2) sits in the (n+1/2) bookkeeping, "
         f"not the dynamics)")
assert g4 < 1e-6
L.append("")

# ---------- the scan
X0S = [0.2, 0.5, 1.0, 2.0, 3.5]
KS = [0.005, 0.02]
RS = 10.0**np.linspace(-1.0, 2.2, 12)
L.append("lam(r) curves (r = b/K; cap b <= 0.4 w0; '--' = capped):")
CURVES = {}
for K in KS:
    for x0 in X0S:
        T = W0/x0
        n0 = n_be(x0)
        n1 = gibbs_mean(T, K, N)
        if abs(n1 - n0) < 1e-9:
            continue
        row = []
        for r in RS:
            b = r*K
            if b > 0.4*W0:
                row.append(None)
                continue
            nb_, _ = ness_nbar(T, K, b, N)
            row.append((nb_ - n0)/(n1 - n0))
        CURVES[(K, x0)] = (n0, n1, row)
        cells = "  ".join("  -- " if v is None else f"{v:5.3f}" for v in row)
        L.append(f"  K={K:.3f} x0={x0:.1f} (n0={n0:6.3f}, n1={n1:6.3f}): "
                 f"{cells}")
L.append(f"  r grid: {['%.2g' % r for r in RS]}")
L.append("")

# ---------- (i) SHARE TEST at fixed r in the transition window
share_pass = False
for K in KS:
    for ri, r in enumerate(RS):
        vals, tanhs = [], []
        for x0 in X0S:
            cv = CURVES.get((K, x0))
            if cv is None or cv[2][ri] is None: break
            vals.append(cv[2][ri])
            tanhs.append(math.tanh(x0/2.0)**2)
        else:
            if all(0.2 < v < 0.8 for v in vals):
                ratios = [v/t for v, t in zip(vals, tanhs)]
                sprd = max(ratios)/min(ratios)
                if sprd < 1.3**2:
                    share_pass = True
                    L.append(f"  share-test PASS candidate at K={K}, "
                             f"r={r:.2g}: lam/tanh^2 spread {sprd:.2f}")
if not share_pass:
    # report the actual spread at the most-transition-like column
    best = None
    for K in KS:
        for ri, r in enumerate(RS):
            vals, tanhs = [], []
            for x0 in X0S:
                cv = CURVES.get((K, x0))
                if cv is None or cv[2][ri] is None: break
                vals.append(cv[2][ri]); tanhs.append(math.tanh(x0/2.0)**2)
            else:
                mid = sum(abs(v - 0.5) for v in vals)
                if best is None or mid < best[0]:
                    rat = [v/t for v, t in zip(vals, tanhs)]
                    best = (mid, K, r, vals, max(rat)/min(rat))
    if best:
        L.append(f"  share test: most-transitional column K={best[1]}, "
                 f"r={best[2]:.2g}: lam = "
                 f"{['%.3f' % v for v in best[3]]}, lam/tanh^2 spread = "
                 f"{best[4]:.1f} (band was < 1.69) -> FAIL")
L.append(f"(i) SHARE TEST: {'PASS' if share_pass else 'FAIL'}")

# ---------- (ii) RESOLUTION TEST: collapse in r_eff = b/(K*spread)
mids = {}
for (K, x0), (n0, n1, row) in CURVES.items():
    xs = [(math.log10(RS[i]), row[i]) for i in range(len(RS))
          if row[i] is not None]
    prev = None
    m = None
    for lg, v in xs:
        if prev is not None:
            (lg0, v0) = prev
            if (v0 - 0.5)*(v - 0.5) <= 0 and v != v0:
                m = lg0 + (0.5 - v0)*(lg - lg0)/(v - v0)
                break
        prev = (lg, v)
    if m is not None:
        T = W0/x0
        nbm, _ = ness_nbar(T, K, (10**m)*K, N)
        spread = math.sqrt(max(nbm*(nbm + 1.0), 1e-12))
        mids[(K, x0)] = (10**m, 10**m/spread)
if mids:
    raw = [v[0] for v in mids.values()]
    eff = [v[1] for v in mids.values()]
    c_raw = max(raw)/min(raw)
    c_eff = max(eff)/min(eff)
    res_pass = c_eff < 2.0
    L.append(f"(ii) RESOLUTION TEST: transition midpoints r50 spread "
             f"raw {c_raw:.1f}x -> in r_eff {c_eff:.1f}x "
             f"(band < 2) -> {'PASS' if res_pass else 'FAIL'}")
    for (K, x0), (r50, re50) in sorted(mids.items()):
        L.append(f"    K={K:.3f} x0={x0:.1f}: r50 = {r50:7.2f}, "
                 f"r_eff50 = {re50:6.2f}")
else:
    res_pass = False
    L.append("(ii) RESOLUTION TEST: no midpoints found -> n/a")
L.append("")

# ---------- (iii) the gravity translation of the resolution shape
L.append("(iii) resolution shape translated to gravity (rho = "
         "x*(sqrt(nu)-1), the self-shift over the bath scale T):")
for xg in (0.05, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
    nu = 1.0 + n_be(xg)
    rho = xg*(math.sqrt(nu) - 1.0)
    L.append(f"    x = {xg:4.2f}: nu = {nu:7.3f}, rho = {rho:.4f}")
L.append("    -> rho vanishes deep (~sqrt(x)) AND in the tail "
         "(~x e^-x/2), peaking ~0.26 at x ~ 1-1.5: a resolution-beta "
         "is zero deep, MAXIMAL at the transition, zero in the tail. "
         "The sky (5T) measures beta = 0 deep AND transition, beta ~ "
         "1/2-3/4 in the TAIL alone -- opposite in both non-deep "
         "regimes, for every O(1) zeta.")
L.append("")

# ---------- verdict
if share_pass:
    call = "PASS-SHARE"
elif res_pass:
    call = "ALT-RES (sky-shape mismatch in the tail: pre-committed " \
           "second-class exclusion applies)"
else:
    call = "AMBIG"
L.append(f"PRE-REGISTERED VERDICT: {call}")

out = "\n".join(L)
print(out)
with open('data/stage6m_analog2.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage6m_analog2.txt")
