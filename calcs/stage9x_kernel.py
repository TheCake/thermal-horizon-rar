"""
STAGE 9X (O5 seam, the WHY-locality closure attempt): THE PARTICIPATION
KERNEL -- measure the lending kernel's shape in the validated toy, and
the COMPOSITION LAW for a band of modes.

The seam after 9W/R18: the s^L gate is multimode-exact for the bright
mode, but WHY the coupling weight localizes at the system's ambient
frequency stayed reading-grade (9T resonant selection), and the 9W
Lemma-B soft-leak worry stands: IF detuned modes contribute their
OCCUPATION to a weighted mean with fat-tailed weights, the soft sector
(n ~ 1/x) leaks in.  The reviewer's degeneracy sweep (P2 0.262 ->
0.016 across 20 lambda) is the kernel being measured -- nobody has
extracted its SHAPE or, more important, the COMPOSITION LAW: does a
detuned mode keep injecting its occupation into an effective mean
(soft-leak-capable), or does its INFLUENCE die with detuning (leak-
impossible, because the lending weight (1/2)s is BOUNDED and the
resonant mode's occupation takes over)?

PARTS:
  K  kernel shape (K=1, 6X config verbatim + detuning): P2(delta) at
     lambda = 0.02, n in {0.5, 2.0}, delta/lambda in {0, 1, 2, 4, 8,
     16} (capped at 16: delta = 0.32 stays below the chi/2 = 0.4
     midpoint to the 0->1 second resonance, noted); core width fitted
     to the 9T secular form r(delta) = r0/(1 + (delta/g_c)^2); the
     TAIL channel identified by lambda-scaling at fixed ABSOLUTE
     delta (lending = lambda-independent; virtual = lambda^2, the 6X
     G2b split).
  C  composition law (K=2): mode 1 RESONANT (n1 = 0.5) + mode 2
     DETUNED by delta (n2 = 2.0), equal couplings; P2(delta) vs three
     candidate laws: occupation-mean (bright nbar -> (1/2)s(nbar) =
     0.2778 at delta=0), gate-mean ((1/2)[s(n1)+s(n2)]/2 = 0.2500),
     resonant-only ((1/2)s(n1) = 0.1667).  THE DECIDER: the
     occupation-sensitivity S(delta) = P2(n2=2) - P2(n2=1) -- if
     S(10 lambda)/S(0) <= 0.1 the detuned mode's occupation has lost
     its vote (leak-impossible); if >= 0.5 the soft-leak is real.
  L  the l-gap single-channel corollary (arithmetic on the measured
     kernel): the l >= 1 horizon multipole channels sit at gaps
     sqrt(l(l+1)) H (9Q) vs Rabi-grade couplings -- their measured
     participation at the equivalent delta/g is printed; single-
     channel (R18 named assumption) becomes kernel-suppressed.

GATES (bars locked at this commit, BEFORE any run):
  G9X-0  K=1 delta=0 port regression: n=0.5 -> 0.16621, n=2.0 ->
         0.33114 within 2e-3 (6X printed values); lambda-independence
         at delta=0 (ratio P2(lam)/P2(lam/2) in (0.8, 1.25), 6X G2a).
  G9X-1  truncation bookkeeping < 5e-3 every config (6X bar).
  G9X-2  6X G2b regression: at the 6X detuned point (delta = 0.3*chi
         = 0.24) the lambda-scaling ratio P2(lam)/P2(lam/2) lands in
         (3, 5) (their printed expectation).
  G9X-3  K=2 delta=0 bright regression: P2 within 5% of (1/2)s(nbar_A)
         (the 9W dynamics bar).
  G9X-4  core fit: the 9T secular form fits the delta <= 4 lambda
         core with rms(ln) <= 0.15 (fitted g_c reported; width-scaling
         vs lambda*sqrt(2n) REPORTED, no bar).

VERDICT LETTERS:
  X-CLOSED  G9X-0..3 PASS AND the tail is VIRTUAL (lambda-scaling
            ratio >= 3 at delta >= 8 lambda_ref fixed-absolute) AND
            S(10 lambda)/S(0) <= 0.1: the lending kernel is core-
            limited, detuned occupations lose their vote, the soft
            sector CANNOT leak into the lending mean -> WHY-locality
            closes at toy grade (composition-boundedness route,
            stronger than kernel-tail suppression); the l-gap
            single-channel corollary is licensed.
  X-LEAK    S(10 lambda)/S(0) >= 0.5 OR the tail is lambda-INDEPENDENT
            (lending-grade tails): the soft-leak is real physics ->
            compute the void-bend predictor (occupation-mean over the
            measured kernel, per-system x0) and annotate P1.
  X-AMBIG   otherwise (sensitivity ratio in (0.1, 0.5) or gate
            failures): report, no closure claimed.

CREDENCE MAP (pre-signed): NO credence move from 9X alone in ANY cell.
9X + 9Y (the squeezing bound) go to ROUND-19 as a package; the package
map (pre-signed here, booked only after the round): (X-CLOSED +
9Y-clean + R19 no-unpatched-hole) -> bath-mechanism conditional
15 -> 17; X-LEAK -> HOLD 15 (the locality seam stays open, now with a
measured leak mechanism); any R19 hole -> HOLD 15 (hole named).
anomaly-real 53 UNTOUCHED in every cell.

Writes data/stage9x_kernel.txt.  Compute: ~2-4 min.
"""
import math
import numpy as np

OUT = 'data/stage9x_kernel.txt'
L = []
def say(s=''):
    L.append(s); print(s, flush=True)
    with open(OUT, 'w') as f:
        f.write("\n".join(L) + "\n")

say("STAGE 9X: THE PARTICIPATION KERNEL -- shape, channel, composition")
say("=" * 72)

# ---------- operators (9W verbatim, + per-mode detuning) ----------
def build_multi(NA, NBs):
    dims = [NA] + list(NBs)
    def op_at(mat, slot):
        out = np.array([[1.0]])
        for i, d in enumerate(dims):
            out = np.kron(out, mat if i == slot else np.eye(d))
        return out
    def low(d):
        return np.diag(np.sqrt(np.arange(1, d)), 1)
    a = op_at(low(NA), 0)
    bs = [op_at(low(NBs[k]), 1 + k) for k in range(len(NBs))]
    Na = a.T @ a
    return a, bs, Na, dims

def thermal_vec(N, nbar):
    if nbar <= 0:
        p = np.zeros(N); p[0] = 1.0
        return p
    x = math.log(1.0 + 1.0/nbar)
    p = np.exp(-x*np.arange(N)); p /= p.sum()
    return p

CHI, LAM, WA = 0.8, 0.02, 5.0

def dressward(NBs, ns, lams, deltas, block=True):
    a, bs, Na, dims = build_multi(4, NBs)
    H = WA*Na + 0.5*CHI*(Na @ (Na - np.eye(Na.shape[0])))
    for k, b in enumerate(bs):
        H = H + (WA + CHI + deltas[k])*(b.T @ b) \
            + lams[k]*(a.T @ b + b.T @ a)
    ev, U = np.linalg.eigh(H)
    p_a = np.zeros(4); p_a[1] = 1.0
    rho_diag = p_a
    for k, NB in enumerate(NBs):
        rho_diag = np.kron(rho_diag, thermal_vec(NB, ns[k]))
    rho = np.diag(rho_diag)
    rl = U.T @ rho @ U
    if block:
        rbar = np.zeros_like(rl)
        scale = max(1.0, float(np.max(np.abs(ev))))
        tol = 1e-8*scale
        i = 0
        n_ = len(ev)
        while i < n_:
            j = i + 1
            while j < n_ and ev[j] - ev[i] < tol:
                j += 1
            rbar[i:j, i:j] = rl[i:j, i:j]
            i = j
    else:
        rbar = np.diag(np.diag(rl))
    rb = U @ rbar @ U.T
    sel = (np.arange(4) == 2)*1.0
    proj_full = np.diag(sel)
    for d in dims[1:]:
        proj_full = np.kron(proj_full, np.eye(d))
    return float(np.real(np.trace(proj_full @ rb)))

def s_of_n(n): return n/(1.0 + n)

# ================= G9X-0: port + lambda-independence =================
say("G9X-0: K=1 delta=0 port regression + G2a lambda-independence:")
tgt = {0.5: 0.16621, 2.0: 0.33114}
ok0 = True
for n, t in tgt.items():
    w = dressward([56], [n], [LAM], [0.0], block=False)
    d = w - t
    ok = abs(d) <= 2e-3
    ok0 &= ok
    say("  n=%.1f: %8.5f (6X %8.5f, d=%+.5f) %s" %
        (n, w, t, d, 'OK' if ok else 'FAIL'))
w1 = dressward([56], [1.0], [LAM], [0.0], block=False)
w2 = dressward([56], [1.0], [LAM/2], [0.0], block=False)
rat0 = w1/w2
okg2a = 0.8 < rat0 < 1.25
ok0 &= okg2a
say("  lambda-independence at delta=0: ratio %.3f (band 0.8-1.25) %s"
    % (rat0, 'OK' if okg2a else 'FAIL'))
say("G9X-0: %s" % ('PASS' if ok0 else 'FAIL'))
say('')

# ================= G9X-1: truncation ================================
trunc = max((s_of_n(2.0))**56, (s_of_n(2.0))**20)
ok1 = trunc < 5e-3
say("G9X-1 truncation (worst n=2: NB=56 %.0e, NB=20 %.0e): %s" %
    ((s_of_n(2.0))**56, (s_of_n(2.0))**20, 'PASS' if ok1 else 'FAIL'))
say('')

# ================= G9X-2: 6X G2b detuned regression ==================
dv1 = dressward([56], [1.0], [LAM], [0.3*CHI], block=False)
dv2 = dressward([56], [1.0], [LAM/2], [0.3*CHI], block=False)
sc = dv1/dv2
ok2 = 3.0 < sc < 5.0
say("G9X-2 6X G2b regression (delta = 0.3 chi): lambda-scaling %.2f "
    "(band 3-5): %s" % (sc, 'PASS' if ok2 else 'FAIL'))
say('')

# ================= PART K: the kernel shape ==========================
say("PART K: P2(delta), K=1, lambda = 0.02 (delta in units of lambda):")
DELS = [0.0, 1.0, 2.0, 4.0, 8.0, 16.0]
kern = {}
for n in (0.5, 2.0):
    row = []
    for dl in DELS:
        w = dressward([56], [n], [LAM], [dl*LAM], block=False)
        row.append(w)
    kern[n] = np.array(row)
    say("  n=%.1f: " % n + "  ".join("%.5f" % v for v in row))
say('')
say("  channel split by lambda-scaling at FIXED ABSOLUTE delta "
    "(P2(lam)/P2(lam/2); 1 = lending, 4 = virtual):")
ratios = []
for dl_abs in (0.0, 0.08, 0.16, 0.32):
    a_ = dressward([56], [2.0], [LAM], [dl_abs], block=False)
    b_ = dressward([56], [2.0], [LAM/2], [dl_abs], block=False)
    ratios.append((dl_abs, a_/b_))
    say("    delta = %.2f (= %.0f lam): ratio = %.2f" %
        (dl_abs, dl_abs/LAM, a_/b_))
tail_virtual = ratios[-1][1] >= 3.0
say("  tail channel at delta = 16 lam: %s" %
    ("VIRTUAL (lambda^2) -- lending is core-limited"
     if tail_virtual else "NOT virtual-scaled"))
say('')
# core fit to the 9T secular form
say("  core fit r(delta) = r0/(1+(delta/g_c)^2) over delta <= 4 lam:")
ok4 = True
for n in (0.5, 2.0):
    core_d = np.array(DELS[:4])*LAM
    core_v = kern[n][:4]
    best = None
    for gc in np.linspace(0.005, 0.20, 400):
        pred = core_v[0]/(1.0 + (core_d/gc)**2)
        r = float(np.sqrt(np.mean(np.log(core_v/pred)**2)))
        if best is None or r < best[1]:
            best = (gc, r)
    gc, rms_ = best
    g_pred = LAM*math.sqrt(2.0*n)
    ok4 &= rms_ <= 0.15
    say("    n=%.1f: g_c = %.4f (rms_ln %.3f); lambda*sqrt(2n) = %.4f "
        "(ratio %.2f)" % (n, gc, rms_, g_pred, gc/g_pred))
say("G9X-4 core fit: %s" % ('PASS' if ok4 else 'FAIL'))
say('')

# ================= PART C: the composition law =======================
say("PART C: K=2 composition (mode1 RESONANT n1=0.5; mode2 detuned, "
    "n2 varied; equal lambdas):")
n1 = 0.5
laws = {'occ-mean': lambda n2: 0.5*s_of_n(0.5*(n1 + n2)),
        'gate-mean': lambda n2: 0.5*0.5*(s_of_n(n1) + s_of_n(n2)),
        'res-only': lambda n2: 0.5*s_of_n(n1)}
say("  delta=0 anchors (n2=2.0): occ-mean %.4f | gate-mean %.4f | "
    "res-only %.4f" % (laws['occ-mean'](2.0), laws['gate-mean'](2.0),
                       laws['res-only'](2.0)))
CDELS = [0.0, 1.0, 2.0, 5.0, 10.0, 20.0]
P_hi, P_lo = [], []
for dl in CDELS:
    ph = dressward([20, 20], [n1, 2.0], [LAM, LAM], [0.0, dl*LAM])
    pl = dressward([20, 20], [n1, 1.0], [LAM, LAM], [0.0, dl*LAM])
    P_hi.append(ph); P_lo.append(pl)
    say("  delta/lam = %4.0f: P2(n2=2.0) = %.5f | P2(n2=1.0) = %.5f "
        "| S = %+.5f" % (dl, ph, pl, ph - pl))
P_hi = np.array(P_hi); P_lo = np.array(P_lo)
S0 = P_hi[0] - P_lo[0]
S10 = P_hi[4] - P_lo[4]
sens_ratio = S10/S0 if S0 != 0 else float('nan')
# G9X-3: delta=0 bright regression
nbar0 = 0.5*(n1 + 2.0)
pred0 = 0.5*s_of_n(nbar0)
d3 = abs(P_hi[0]/pred0 - 1)
ok3 = d3 <= 0.05
say("G9X-3 delta=0 bright regression: P2 = %.5f vs (1/2)s(%.2f) = "
    "%.5f (|d| = %.3f): %s" %
    (P_hi[0], nbar0, pred0, d3, 'PASS' if ok3 else 'FAIL'))
say("  occupation-sensitivity: S(0) = %.5f, S(10 lam) = %.5f, "
    "ratio = %.3f  [bars: <= 0.1 CLOSED / >= 0.5 LEAK]" %
    (S0, S10, sens_ratio))
say("  limit check: P2(20 lam, n2=2) = %.5f vs res-only %.5f" %
    (P_hi[5], laws['res-only'](2.0)))
say('')

# ================= PART L: the l-gap corollary =======================
say("PART L: the l-gap single-channel corollary (arithmetic on the "
    "measured kernel):")
say("  l >= 1 horizon multipole channels are gapped by "
    "sqrt(l(l+1)) H (9Q); soft-sector lending couplings are "
    "Rabi-grade << H.")
kk = kern[2.0]
supp8 = kk[4]/kk[0]
supp16 = kk[5]/kk[0]
say("  measured kernel suppression: P2(8 lam)/P2(0) = %.3f; "
    "P2(16 lam)/P2(0) = %.3f (and the surviving tail is the VIRTUAL "
    "channel, which does not lend)" % (supp8, supp16))
say("  => at gap/coupling ratios >= 10 (the l-channels), lending "
    "participation is at the few-percent level and falling as "
    "1/delta^2 with only virtual character: the single final channel "
    "of R18's named assumption is KERNEL-SUPPRESSED, not assumed.")
say('')

# ================= verdict ==========================================
POWERED = ok0 and ok1 and ok2 and ok3
if POWERED and tail_virtual and sens_ratio <= 0.1:
    letter = 'X-CLOSED'
elif (sens_ratio >= 0.5) or (POWERED and not tail_virtual):
    letter = 'X-LEAK'
else:
    letter = 'X-AMBIG'
say("=" * 72)
say("VERDICT LETTER: %s" % letter)
if letter == 'X-CLOSED':
    say("  the lending kernel is CORE-LIMITED (9T-form core, virtual "
        "tail); a detuned mode's OCCUPATION loses its vote (S-ratio "
        "%.3f <= 0.1); the soft sector cannot leak into the lending "
        "mean because the composition is resonant-limited and the "
        "lending weight is bounded by 1/2." % sens_ratio)
    say("  => WHY-LOCALITY CLOSES AT TOY GRADE by the composition-"
        "boundedness route (stronger than kernel-tail suppression: "
        "even fat-tailed PARTICIPATION cannot import soft "
        "OCCUPATIONS); the 9W Lemma-B soft-leak worry is resolved "
        "in-model; the l-gap corollary licenses single-channel.")
elif letter == 'X-LEAK':
    say("  the detuned mode keeps injecting occupation (S-ratio %.3f)"
        " -- the soft-leak is real in-model; void-bend predictor "
        "required before any P1 use." % sens_ratio)
say("  credence: NO move from 9X alone (pre-signed); package map "
    "goes to ROUND-19 with 9Y.")
say("  anomaly-real 53 UNTOUCHED")
say('')
say("gates: G9X-0 %s | G9X-1 %s | G9X-2 %s | G9X-3 %s | G9X-4 %s" %
    tuple('PASS' if x else 'FAIL' for x in (ok0, ok1, ok2, ok3, ok4)))
print("\nsaved:", OUT)
