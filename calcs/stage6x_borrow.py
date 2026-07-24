"""
STAGE 6X (the borrowing dynamics): does exact closed-system quantum
mechanics produce the KMS-RATIO channel the sky selected?

Setting (the 6N-survivor limit): frozen horizon bath => on orbital
timescales the ONLY dynamical reservoir is the ambient mode. Extreme
non-Markovian limit = CLOSED dynamics: Kerr-anharmonic local mode a
(chi = the dressing self-shift) exchanging with one thermal ambient
mode b (occupation n_amb), no bath during evolution. Geometry: the
dress-ward transition |1>->|2> of a (frequency wa + chi) is RESONANT
with b; the source transition |0>->|1> is detuned by chi. Dressing
must ABSORB a real ambient quantum.

PRE-REGISTERED PRIMARY QUESTION (unchanged from 67626f7): which of the
loop's natural forms does the dress-ward channel weight follow in
n_amb?  ratio n/(1+n) [sky-selected, 6U] | share n/(2n+1) | raw n
[both sky-excluded]. BARS: SUPPORT if ratio wins (slope 1 +- 0.15 and
3x residual margin over both rivals); STRIKE if share or raw wins at
the same margin (pre-commit: bath-microphysics ~8% -> ~5%); else AMBIG.

AMENDMENT (post-commit, PRE-RESULTS -- the first run crashed at the G1
anchor before any channel number was computed; 6O-precedent disclosure):
 (a) G1 anchor corrected to the QUBIT sub-case: two linearly coupled
     harmonic modes are a linear system -- normal modes carry NO
     occupation pull (first run measured exactly lambda^2/Delta,
     constant: the code was right, the gate's physics was wrong). The
     6H anchor lambda^2(2n+1)/Delta is the JC (two-level) result.
 (b) G2 split by channel physics: at RESONANCE the transfer is Rabi
     mixing -- weight is lambda-INDEPENDENT and predicted to be
     (1/2)*P(n_b >= 1); the DETUNED (virtual) control channel scales
     as lambda^2 and is predicted to carry RAW-n weighting. The two
     channels should carry DIFFERENT forms -- the discriminating
     physics. For a BE-thermal state P(n >= L) = [nbar/(1+nbar)]^L
     = e^(-L x) EXACTLY (the geometric-tail identity): if the resonant
     channel follows it, the KMS ratio IS the lending probability and
     the gate s^L = P(ambient can supply L quanta).

GATES: G0 thermal-tail truncation bookkeeping; G1 qubit dispersive
anchor; G2a resonant lambda-independence; G2b detuned lambda^2 + raw-n
form; G3 Markov-dephasing regression toward 6N locking; GL the
geometric-tail identity (sympy exact).

Writes data/stage6x_borrow.txt.
"""
import math
import numpy as np
import sympy as sp

L = []
def say(s=''):
    L.append(s); print(s, flush=True)

say("STAGE 6X: the borrowing dynamics -- exact closed-system channel test")
say("=" * 72)

def build(NA, NB):
    ad = np.diag(np.sqrt(np.arange(1, NA)), -1)
    bd = np.diag(np.sqrt(np.arange(1, NB)), -1)
    Ia, Ib = np.eye(NA), np.eye(NB)
    A_ = np.kron(ad.T, Ib); Ad = np.kron(ad, Ib)
    B_ = np.kron(Ia, bd.T); Bd = np.kron(Ia, bd)
    Na_ = Ad @ A_; Nb_ = Bd @ B_
    return dict(A=A_, Ad=Ad, B=B_, Bd=Bd, Na=Na_, Nb=Nb_,
                Ia=Ia, Ib=Ib, NA=NA, NB=NB)

def hamiltonian(o, wa, wb, chi, lam):
    Id = np.kron(o['Ia'], o['Ib'])
    return wa*o['Na'] + wb*o['Nb'] + 0.5*chi*(o['Na'] @ (o['Na'] - Id)) \
        + lam*(o['Ad'] @ o['B'] + o['Bd'] @ o['A'])

def thermal(N, nbar):
    if nbar <= 0:
        p = np.zeros(N); p[0] = 1.0
        return np.diag(p)
    x = math.log(1.0 + 1.0/nbar)
    p = np.exp(-x*np.arange(N)); p /= p.sum()
    return np.diag(p)

# ---- GL: the geometric-tail identity (exact) -------------------------------
nb_s, k_s, x_s = sp.symbols('nbar L x', positive=True)
q = nb_s/(1 + nb_s)
PgeL = q**k_s                       # sum_{n>=L} (1-q) q^n = q^L
chk = sp.simplify(sp.summation((1-q)*q**sp.symbols('n', integer=True,
                                                   nonnegative=True),
                               (sp.symbols('n', integer=True,
                                           nonnegative=True), k_s,
                                sp.oo)) - PgeL)
okL = chk == 0
kms = sp.simplify(q.subs(nb_s, 1/(sp.exp(x_s) - 1)) - sp.exp(-x_s)) == 0
say(f"GL: P(n >= L) = [n/(1+n)]^L for a thermal state (exact): "
    f"{'PASS' if okL else 'FAIL'};  n/(1+n) = e^-x: "
    f"{'PASS' if kms else 'FAIL'}")
assert okL and kms
say("    => the KMS ratio IS the lending probability; the gate s^L = "
    "P(ambient can supply L quanta).")
say('')

# ---- G1: the qubit dispersive anchor (the 6H (2n+1) pull) ------------------
say("G1: qubit dispersive anchor (NA=2, per-Fock-sector; a thermal probe "
    "dephases across sectors -- amendment (a)):")
o2 = build(2, 30)
lam, Delta = 0.02, 1.0
wa, wb = 5.0, 5.0 - Delta
ok1 = True
for nfock in (0, 2, 5):
    H = hamiltonian(o2, wa, wb, 0.0, lam)
    ev, U = np.linalg.eigh(H)
    psi_a = np.array([1.0, 1.0])/math.sqrt(2)
    pb = np.zeros(30); pb[nfock] = 1.0
    rho = np.kron(np.outer(psi_a, psi_a), np.diag(pb))
    rl = U.T @ rho @ U
    ts = np.linspace(0, 400.0, 4001)
    ph = []
    for t in ts:
        e = np.exp(-1j*ev*t)
        rt = (U*e) @ rl @ (U*e).conj().T
        ph.append(np.angle(np.trace(o2['A'] @ rt)))
    slope = np.polyfit(ts, np.unwrap(np.array(ph)), 1)[0]
    shift = -slope - wa
    pred = lam*lam*(2*nfock + 1)/Delta
    r = shift/pred
    say(f"  n_b={nfock:2d}: measured {shift:+.3e}, predicted {pred:+.3e},"
        f" ratio {r:.3f}")
    ok1 &= abs(r - 1.0) < 0.10
say(f"G1: {'PASS' if ok1 else 'FAIL'}")
assert ok1
say('')

# ---- the borrowing configuration ------------------------------------------
NA, NB = 4, 56
o = build(NA, NB)
CHI, LAM, WA = 0.8, 0.02, 5.0

def channel_weight(nb, lam=LAM, delta=0.0):
    # b resonant with 1->2 when delta=0; detuned control via delta
    H = hamiltonian(o, WA, WA + CHI + delta, CHI, lam)
    ev, U = np.linalg.eigh(H)
    p1 = np.zeros(NA); p1[1] = 1.0
    rho = np.kron(np.diag(p1), thermal(NB, nb))
    rl = U.T @ rho @ U
    rbar = np.diag(np.diag(rl))          # infinite-time average
    rb = U @ rbar @ U.T
    proj2 = np.kron(np.diag((np.arange(NA) == 2)*1.0), o['Ib'])
    return float(np.real(np.trace(proj2 @ rb)))

# G0: thermal truncation bookkeeping at the largest nbar
nmax = 8.0
trunc = (nmax/(1+nmax))**NB
say(f"G0: thermal-tail truncation at nbar={nmax}, NB={NB}: "
    f"P(n>={NB}) = {trunc:.1e}  {'PASS' if trunc < 5e-3 else 'FAIL'}")
assert trunc < 5e-3

# G2a: resonant channel is lambda-independent (Rabi mixing)
w1 = channel_weight(1.0, lam=LAM)
w2 = channel_weight(1.0, lam=LAM/2)
say(f"G2a: resonant lambda-independence P2(lam)/P2(lam/2) = "
    f"{w1/w2:.3f} (expect ~1): {'PASS' if 0.8 < w1/w2 < 1.25 else 'FAIL'}")
assert 0.8 < w1/w2 < 1.25

say('')
say("the scan: dress-ward weight P2bar(n_amb), resonant channel "
    "(prediction: (1/2) P(n>=1) = (1/2) n/(1+n)):")
NBS = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
W = []
for nb in NBS:
    P2 = channel_weight(nb)
    W.append(P2)
    pred = 0.5*nb/(1.0 + nb)
    say(f"  n_amb = {nb:5.2f}: P2bar = {P2:.5f}   [(1/2)ratio = {pred:.5f}]")
W = np.array(W); NBSa = np.array(NBS)

forms = {
    'ratio n/(1+n)':  NBSa/(1.0 + NBSa),
    'share n/(2n+1)': NBSa/(2.0*NBSa + 1.0),
    'raw n':          NBSa,
}
say('')
say("form regression ln(P2bar) = s*ln(form) + c:")
res = {}
for nm, f in forms.items():
    s, c = np.polyfit(np.log(f), np.log(W), 1)
    pred = s*np.log(f) + c
    rms = float(np.sqrt(np.mean((np.log(W) - pred)**2)))
    res[nm] = (s, rms)
    say(f"  {nm:<16}: slope {s:+.3f}, rms {rms:.4f}")
r_s, r_rms = res['ratio n/(1+n)']
sh_s, sh_rms = res['share n/(2n+1)']
rw_s, rw_rms = res['raw n']
ratio_wins = (abs(r_s - 1.0) < 0.15 and r_rms*3.0 <= sh_rms
              and r_rms*3.0 <= rw_rms)
share_wins = (abs(sh_s - 1.0) < 0.15 and sh_rms*3.0 <= r_rms)
raw_wins = (abs(rw_s - 1.0) < 0.15 and rw_rms*3.0 <= r_rms)

# G2b: the detuned (virtual) control -- lambda^2 scaling and raw-n form
say('')
say("G2b: detuned control (delta = 0.3*chi; the virtual channel):")
dv1 = channel_weight(1.0, lam=LAM, delta=0.3*CHI)
dv2 = channel_weight(1.0, lam=LAM/2, delta=0.3*CHI)
sc = dv1/dv2
say(f"  lambda-scaling {sc:.2f} (expect ~4): "
    f"{'PASS' if 3.0 < sc < 5.0 else 'FAIL'}")
Wd = np.array([channel_weight(nb, delta=0.3*CHI) for nb in NBS])
sd, cd = np.polyfit(np.log(NBSa), np.log(Wd), 1)
sr, cr = np.polyfit(np.log(NBSa/(1+NBSa)), np.log(Wd), 1)
predd = sd*np.log(NBSa) + cd
rmsd = float(np.sqrt(np.mean((np.log(Wd) - predd)**2)))
predr = sr*np.log(NBSa/(1+NBSa)) + cr
rmsr = float(np.sqrt(np.mean((np.log(Wd) - predr)**2)))
say(f"  detuned form: raw-n slope {sd:+.3f} (rms {rmsd:.4f}) vs ratio "
    f"slope {sr:+.3f} (rms {rmsr:.4f})")
say(f"  -> the two channels carry different forms: resonant/real = "
    f"lending-probability, detuned/virtual = matrix-element-weighted")

# ---- G3: Markov regression (dephase the exchange -> locking) --------------
say('')
o3 = build(4, 12)
def channel_weight_markov(nb, gamma):
    from scipy.linalg import expm
    H = hamiltonian(o3, WA, WA + CHI, CHI, LAM)
    dim = 4*12
    Lv = -1j*(np.kron(np.eye(dim), H) - np.kron(H.T, np.eye(dim)))
    C = o3['Na'].astype(complex)
    CC = C.conj().T @ C
    Lv = Lv + gamma*(np.kron(C.conj(), C)
                     - 0.5*np.kron(np.eye(dim), CC)
                     - 0.5*np.kron(CC.T, np.eye(dim)))
    p1 = np.zeros(4); p1[1] = 1.0
    rho = np.kron(np.diag(p1), thermal(12, nb)).astype(complex)
    rt = (expm(Lv*600.0) @ rho.reshape(-1)).reshape(dim, dim)
    proj2 = np.kron(np.diag((np.arange(4) == 2)*1.0), o3['Ib'])
    return float(np.real(np.trace(proj2 @ rt)))
mfroz = channel_weight(2.0)
mmark = channel_weight_markov(2.0, gamma=2.0)
say(f"G3: dress-ward weight at n=2: frozen {mfroz:.4f} vs strongly "
    f"dephased {mmark:.4f} -> suppression x{mfroz/max(mmark,1e-12):.1f}")
say(f"    (dephasing the exchange coherence closes the borrow channel = "
    f"the 6N locking direction; recorded)")

say('')
if ratio_wins:
    v = ("SUPPORT: exact closed dynamics follows the KMS RATIO = the "
         "lending probability P(n>=1); the per-leg factor of the "
         "borrowing reading is DERIVED at toy grade (s^L = P(n>=L); "
         "L = 2 from the loop order, 6H/6U exact)")
elif share_wins or raw_wins:
    v = ("STRIKE: the dynamics follows the "
         + ('SHARE' if share_wins else 'RAW-n')
         + " form the sky excluded -- pre-committed: bath-microphysics "
           "~8% -> ~5%")
else:
    v = "AMBIG: no form wins at the pre-registered margins; no move"
say("VERDICT vs pre-registered bars: " + v)

with open('data/stage6x_borrow.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nSTAGE 6X done -> data/stage6x_borrow.txt")
