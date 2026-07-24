"""
STAGE 7E: the platform translation -- the borrowing/lending gate written
in circuit-QED language, with real numbers.

HONEST FRAME (stated first, per 6N): nothing here tests gravity. Our own
cancellation theorem (6N) says buildable flat-bath configurations measure
ZERO occupation pull, and the sky's tail-only beta is anti-standard. What
IS testable in a lab is the MECHANISM CLASS the derivation lives in:
  (i)  the LENDING-PROBABILITY gate (6X): a Kerr mode beam-splitter-
       coupled to ONE thermal ancilla transfers dress-ward with weight
       (1/2) P(n>=1) = (1/2) nbar/(1+nbar) -- saturating, not linear;
       and the two-quantum channel carries P(n>=2) = [nbar/(1+nbar)]^2 --
       the geometric-tail rung (the sky's s^L) as a measurable curve;
  (ii) the RESOLUTION crossover (6N): source-locking under a flat bath
       vs KMS-structured occupation under a resolved (filtered) bath.
A positive lab result validates the grammar's building blocks as real
quantum optics; it does NOT confirm the horizon reading. A negative
result at clean parameters WOULD strike the 6X mechanics.

PLATFORM (primary: 3D-cavity circuit QED, standard values):
  Kerr mode  = transmon:  omega_a/2pi = 5.000 GHz, chi/2pi = 250 MHz
                          T1 = 200 us, Tphi = 100 us
  Ancilla    = 3D cavity: omega_b tuned; kappa/2pi = 5 kHz (Q = 1e6)
  Coupling   = parametric beam-splitter conversion: lambda/2pi = 2 MHz
  Thermal    = injected calibrated noise: nbar in [0.25, 8]
               (photon-number-splitting calibration on the transmon)
  Protocol   : thermalize cavity (~5/kappa ~ 160 us) -> reset + pi-pulse
               transmon to |1> -> pump beam-splitter for T = 2 us ->
               measure P(|2>) (dispersive readout); repeat vs nbar.
  L=1 config : omega_b = omega_a + chi     (|1,n> -> |2,n-1> resonant)
  L=2 config : omega_b = omega_a + 1.5 chi (|1,n> -> |3,n-2> two-quantum
               resonant; single-quantum channels detuned by chi/2 >> lam)

GATES:
  GP0 scale separations printed (chi >> lam >> kappa(1+2nbar), Gam1*T).
  GP1 regression: closed-system channel weight at PLATFORM parameters
      reproduces (1/2) nbar/(1+nbar) to < 2.5% over the nbar grid
      (6X band was 0.2-2.2% at a coarser lam/chi).
  GP2 dissipation: Lindblad (cavity kappa up/down at nbar, transmon
      Gamma1 + dephasing) time-averaged weight within 10% of the
      closed-system value at nbar = 1 and 4 -- else report the required
      kappa/lambda. Sparse expm_multiply.
  GP3 the discriminators: (a) saturation visibility -- ratio-law vs
      linear extrapolation from the lowest-nbar point, >= 3x at nbar = 8;
      (b) the L=2 geometric rung: level-3 weight regresses on
      [n/(1+n)]^2 with slope 1 +- 0.2 over the grid.
  GP4 ion-trap + optomechanics comparison rows (order-of-magnitude,
      labeled as such) -- cQED is the platform.
Writes data/stage7e_platform.txt.
"""
import math
import numpy as np
from scipy.sparse import kron as skron, identity as sid, csc_matrix
from scipy.sparse.linalg import expm_multiply

L = []
def say(s=''):
    L.append(s); print(s, flush=True)
def save():
    with open('data/stage7e_platform.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

say("STAGE 7E: the lending gate in circuit-QED language")
say("=" * 72)

# ---- parameters (2pi*MHz angular units; time in us) ------------------------
WA = 2*math.pi*5000.0
CHI = 2*math.pi*250.0
LAM = 2*math.pi*2.0
KAP = 2*math.pi*0.005
G1 = 1.0/200.0
GPHI = 1.0/100.0
TAVG = 2.0
NBARS = [0.25, 0.5, 1.0, 2.0, 4.0, 8.0]

say("platform card (3D circuit QED):")
say(f"  transmon omega_a/2pi = 5.000 GHz, chi/2pi = 250 MHz, "
    f"T1 = 200 us, Tphi = 100 us")
say(f"  cavity kappa/2pi = 5 kHz (Q = 1e6); beam-splitter lambda/2pi = 2 MHz")
say(f"  averaging window T = {TAVG} us; nbar grid {NBARS}")
say("GP0 scale separations:")
say(f"  chi/lambda = {CHI/LAM:.0f} (>> 1: channel selectivity)")
say(f"  lambda/kappa = {LAM/KAP:.0f} (>> 1: exchange beats cavity decay)")
for nb in (1.0, 8.0):
    say(f"  kappa(1+2nbar)*T at nbar={nb:g}: {KAP*(1+2*nb)*TAVG:.3f}")
say(f"  Gamma1*T = {G1*TAVG:.3f}, Gphi*T = {GPHI*TAVG:.3f}")
say('')
save()

# ---- operators -------------------------------------------------------------
def build(NA, NB):
    ad = np.diag(np.sqrt(np.arange(1, NA)), -1)
    bd = np.diag(np.sqrt(np.arange(1, NB)), -1)
    Ia, Ib = np.eye(NA), np.eye(NB)
    A_ = np.kron(ad.T, Ib); Ad = np.kron(ad, Ib)
    B_ = np.kron(Ia, bd.T); Bd = np.kron(Ia, bd)
    return dict(A=A_, Ad=Ad, B=B_, Bd=Bd, Na=Ad @ A_, Nb=Bd @ B_,
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

# ---- GP1: closed-system regression at platform parameters ------------------
NA, NB = 4, 56
o = build(NA, NB)
def weight(nb, wb_off, level, lam=LAM):
    H = hamiltonian(o, WA, WA + wb_off, CHI, lam)
    ev, U = np.linalg.eigh(H)
    p1 = np.zeros(NA); p1[1] = 1.0
    rho = np.kron(np.diag(p1), thermal(NB, nb))
    rl = U.T @ rho @ U
    rb = U @ np.diag(np.diag(rl)) @ U.T
    proj = np.kron(np.diag((np.arange(NA) == level)*1.0), o['Ib'])
    return float(np.real(np.trace(proj @ rb)))

say("GP1 closed-system regression, L=1 config (omega_b = omega_a + chi):")
ok1 = True
W1 = []
for nb in NBARS:
    w = weight(nb, CHI, 2)
    pred = 0.5*nb/(1.0 + nb)
    r = w/pred
    W1.append(w)
    ok1 &= abs(r - 1.0) < 0.025
    say(f"  nbar = {nb:5.2f}: P2bar = {w:.5f}  vs (1/2)n/(1+n) = "
        f"{pred:.5f}  ratio {r:.4f}")
say(f"GP1: {'PASS' if ok1 else 'FAIL'} (band < 2.5%)")
save()
assert ok1
say('')

# ---- GP2: dissipative run (sparse Lindblad, NB=12) -------------------------
say("GP2 dissipation at platform rates (kappa up/down, Gamma1, Gphi):")
NA2, NB2 = 4, 12
o2 = build(NA2, NB2)
H2 = hamiltonian(o2, WA, WA + CHI, CHI, LAM)
dim = NA2*NB2
Isp = sid(dim, format='csc')
def lind(Cd):
    C = csc_matrix(Cd)
    CC = (C.conj().T @ C).tocsc()
    return skron(C.conjugate(), C, format='csc') \
        - 0.5*skron(Isp, CC, format='csc') \
        - 0.5*skron(CC.transpose(), Isp, format='csc')
ok2 = True
for nb in (1.0, 4.0):
    Lv = -1j*(skron(Isp, csc_matrix(H2), format='csc')
              - skron(csc_matrix(H2).transpose(), Isp, format='csc'))
    Lv = Lv + KAP*(1+nb)*lind(o2['B']) + KAP*nb*lind(o2['Bd'])
    Lv = Lv + G1*lind(o2['A']) + GPHI*lind(o2['Na'])
    p1 = np.zeros(NA2); p1[1] = 1.0
    rho0 = np.kron(np.diag(p1), thermal(NB2, nb)).astype(complex).reshape(-1)
    ts = np.linspace(0.0, TAVG, 26)[1:]
    rts = expm_multiply(Lv, rho0, start=ts[0], stop=ts[-1], num=len(ts))
    proj = np.kron(np.diag((np.arange(NA2) == 2)*1.0), o2['Ib'])
    p2s = [float(np.real(np.trace(proj @ rt.reshape(dim, dim))))
           for rt in rts]
    wdis = float(np.mean(p2s))
    # closed-system reference at the SAME truncation NB2
    o2c = o2
    Hc = H2
    ev, U = np.linalg.eigh(Hc)
    rl = U.T @ np.kron(np.diag(p1), thermal(NB2, nb)) @ U
    rb = U @ np.diag(np.diag(rl)) @ U.T
    wcl = float(np.real(np.trace(proj @ rb)))
    shift = wdis/wcl - 1.0
    ok2 &= abs(shift) < 0.10
    say(f"  nbar = {nb:g}: dissipative time-avg = {wdis:.5f} vs closed "
        f"(same NB) = {wcl:.5f}  shift = {100*shift:+.1f}%")
say(f"GP2: {'PASS' if ok2 else 'FAIL'} (|shift| < 10%)")
save()
say('')

# ---- GP3a: saturation visibility ------------------------------------------
say("GP3a saturation discriminator (lending law vs linear response):")
w0 = W1[0]; nb0 = NBARS[0]
lin8 = w0*(NBARS[-1]/nb0)
vis = lin8/W1[-1]
say(f"  linear extrapolation from nbar={nb0}: predicts {lin8:.3f} at "
    f"nbar=8; lending law measures {W1[-1]:.3f} -> separation x{vis:.1f}")
ok3a = vis >= 3.0
say(f"GP3a: {'PASS' if ok3a else 'FAIL'} (>= 3x)")
say('')

# ---- GP3b: the L=2 geometric rung -----------------------------------------
say("GP3b the two-quantum rung, L=2 config (omega_b = omega_a + 1.5 chi):")
W2 = []
for nb in NBARS:
    w = weight(nb, 1.5*CHI, 3)
    W2.append(w)
    say(f"  nbar = {nb:5.2f}: P3bar = {w:.6f}   [n/(1+n)]^2 = "
        f"{(nb/(1+nb))**2:.4f}")
W2 = np.array(W2); NBSa = np.array(NBARS)
q2 = (NBSa/(1+NBSa))**2
s2, c2 = np.polyfit(np.log(q2), np.log(W2), 1)
rms2 = float(np.sqrt(np.mean((np.log(W2) - (s2*np.log(q2)+c2))**2)))
sr, cr = np.polyfit(np.log(NBSa), np.log(W2), 1)
rmsr = float(np.sqrt(np.mean((np.log(W2) - (sr*np.log(NBSa)+cr))**2)))
ok3b = abs(s2 - 1.0) < 0.2 and rms2 < rmsr
say(f"  regression on [n/(1+n)]^2: slope {s2:+.3f} rms {rms2:.4f}  "
    f"(raw-n^1: slope {sr:+.3f} rms {rmsr:.4f})")
say(f"  prefactor e^c = {math.exp(c2):.4f} (recorded; the FORM is the claim)")
say(f"GP3b: {'PASS' if ok3b else 'FAIL'} (slope 1 +- 0.2 on the "
    f"geometric-tail rung)")
save()
say('')

# ---- GP4: platform comparison + the resolution-crossover design ------------
say("GP4 platform comparison (order-of-magnitude, stated as such):")
say("  cQED (3D):   chi/2pi 200-300 MHz, lam/2pi 1-10 MHz, kappa/2pi")
say("               1-50 kHz, nbar 0.01-8 calibrated -> ALL separations")
say("               met with standard hardware; THE platform.")
say("  trapped ion: effective Kerr (laser-induced) chi/2pi ~ 1-10 kHz,")
say("               exchange g ~ 0.1-1 kHz, heating 10-100 quanta/s;")
say("               engineered thermal reservoirs excellent, but")
say("               chi/g ~ 10 marginal for channel selectivity.")
say("  optomech:    single-photon Kerr negligible; not suitable.")
say('')
say("the 6N resolution-crossover on the same chip (design note):")
say(f"  flat regime:     broadband 50-Ohm environment, bandwidth >> chi")
say(f"                   (ratio > 10) -> theorem predicts SOURCE-LOCKING")
say(f"                   (zero occupation pull, 6N max|lambda| grade)")
say(f"  resolved regime: Purcell-style filter, linewidth/2pi = 1-10 MHz")
say(f"                   << chi/2pi = 250 MHz (ratio 25-250) -> the")
say(f"                   KMS-contrast carrier appears (6N)")
say("  both regimes are switchable on one device (tunable coupler to")
say("  either port); the crossover parameter chi/bandwidth spans ~0.1")
say("  to ~250 across settings.")
say('')

verd = "FEASIBLE" if (ok1 and ok2 and ok3a and ok3b) else "STRAINED"
say(f"VERDICT: the platform translation is {verd} at standard cQED")
say("parameters: the lending law (1/2) n/(1+n), its saturation, and the")
say("L=2 geometric rung [n/(1+n)]^2 are all measurable curves; the")
say("dissipative shift at the chosen rates is quantified above. Honest")
say("frame: validates the MECHANISM CLASS (6X grammar), not the sky rule")
say("(6N cancellation theorem); a clean negative WOULD strike 6X.")
save()
print("\nSTAGE 7E done -> data/stage7e_platform.txt")
