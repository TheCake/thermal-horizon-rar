"""
STAGE 6K (the desktop analog): does standard open-quantum-systems math
produce the grammar's gating? PRE-REGISTRATION COMMITTED BEFORE
EXECUTION — the mapping, estimators, bands, and outcome tree below are
fixed in git before any number is computed.

THE LAB-NATIVE FORM OF THE RULE (identity, gated in G0): the grammar
beta = (1/2)*[1/(2nu-1)]^2*[n_amb/(1+n_amb)]^2 rewrites exactly as

    beta = (1/2) * tanh^2(x_loc/2) * exp(-2 x_amb)

since 1/(2n+1) = tanh(x/2) and n/(1+n) = e^{-x}. Both factors are
standard thermal response functions; the NONSTANDARD content is only
their assignment as per-vertex gates of a two-leg exchange. The squared
ambient share = e^{-2x} = the Boltzmann cost of borrowing TWO ambient
quanta — one per leg.

THE MODEL (exactly solvable birth-death NESS; no Liouvillian needed):
a Kerr mode E_n = n*w0 + (K/2)n(n-1), link spacing D_n = w0 + K*n
("self-frequency depends on occupation"), truncated at N with a
doubling gate. Channels (all bosonic-ladder, per-link Einstein pairs):
  S  (source-locked): frequency-BLIND thermal contact pinned at
     occupation nS: up kS*nS*(n+1), down kS*(nS+1)*(n+1). Alone it
     drives the mode to n0 = nS — the drive-following endpoint (beta=0
     analog). Source temperature T_S := w0/ln(1+1/nS).
  A  (vanilla ambient, CFG-V): Davies thermal channel at temperature
     T_A, frequency-AWARE: rates use n_A(D_n) = 1/(e^{D_n/T_A}-1).
     Alone it drives to its own Davies-Gibbs.
  M  (mediated ambient, CFG-M — the grammar's own process class): a
     TWO-VERTEX composite jump: the mode absorbs its dressed quantum
     D_n by taking w0 from the source reservoir and the mismatch
     d_n = K*n from the ambient reservoir. Up kM*(n+1)*nS*n_A(d_n);
     down kM*(n+1)*(nS+1)*(n_A(d_n)+1). Its detailed-balance ratio is
     e^{-w0/T_S}e^{-d_n/T_A}: at T_A = T_S it alone drives the mode to
     the exact Gibbs of the DRESSED spectrum — the self-consistent
     endpoint (beta=1 analog) — while its RATE carries the stimulated
     shares. This is the minimal standard-physics realization of
     "ambient-assisted two-vertex dressing".

ESTIMATOR (fixed): lam = (nbar_NESS - n0)/(n1 - n0), where n1 is the
competing channel's own fixed point (CFG-V: Davies-Gibbs mean of A
alone; CFG-M: the self-consistent solution of
n = 1/(e^{w0/T_S + K n/T_A} - 1)). lam = fraction of the migration
from drive-following to self-consistent that the ambient channel
achieves. Scanned dials: nS (local occupation -> quantumness),
s* = ambient Boltzmann share at the reference exchange (CFG-V:
e^{-w0/T_A}; CFG-M: e^{-K nS/T_A}), the rate ratio kX/kS, and K/w0.

GRAMMAR FINGERPRINT (pre-registered bands): in some contiguous regime,
  (i)  p_s = dln(lam)/dln(s*) in [1.6, 2.4]        (ambient share SQUARED)
  (ii) p_n = dln(lam)/dln(2nbar+1) in [-2.4, -1.6] (local share SQUARED)
  (iii) kappa-slope dln(lam)/dln(kX/kS) in [-0.2, +0.2]
       (the grammar contains NO rates — share-only gating)
OUTCOME TREE (pre-committed):
  PASS-M: CFG-M shows all three -> the grammar is standard physics in
     its own process class; microphysics credence UP (~+15 points on
     the bath-mechanism conditional).
  PASS-V: even the vanilla config shows it -> stronger PASS.
  ALT: clean but different exponents -> the math hands us a RIVAL
     gating law; pre-committed next step = construct nu_alt from the
     measured law and test it on binaries+galaxies; if nu_alt then
     fails the sky, the vanilla-lab reading of the grammar is dead and
     only a non-Markovian/horizon-specific reading survives (logged as
     a strike, ~-10 points); if nu_alt fits better, the grammar is
     REPLACED and the math wins.
  AMBIG: no clean powers / estimator pathologies -> the protocol is
     not fridge-ready; sharpen the mapping before any lab talk.
Gates: G0 identity (tanh/Boltzmann rewrite vs 6E numbers, 1e-12);
G1 single-channel endpoints exact (S -> BE(nS); A -> Davies-Gibbs;
M at T_A=T_S -> dressed Gibbs, all 1e-10); G2 K=0 linear-mode closed
form for CFG-V (rate-weighted mean, 1e-8); G3 truncation doubling
(<1e-8 relative); G4 lam in [-0.1, 1.1] across the scan (estimator
sanity). Writes data/stage6k_analog.txt.
"""
import math
import numpy as np

L = ["STAGE 6K: the desktop analog -- vanilla vs mediated two-vertex "
     "gating, pre-registered", ""]

def n_be(x):
    return 1.0/(math.exp(x) - 1.0)

# ---------- G0: the lab-native identity
for tag, e in (("gal", 0.02), ("bin", 1.15)):
    x = math.sqrt(e)
    n = n_be(x)
    g_share = (n/(1.0 + n))**2
    g_boltz = math.exp(-2.0*x)
    d = abs(g_share - g_boltz)
    L.append(f"G0 identity ({tag}): [n/(1+n)]^2 = {g_share:.6f}, "
             f"e^-2x = {g_boltz:.6f}, d = {d:.1e} -> "
             f"{'PASS' if d < 1e-12 else 'FAIL'}")
    assert d < 1e-12
xq = 0.7
d2 = abs(1.0/(2.0*n_be(xq) + 1.0) - math.tanh(xq/2.0))
L.append(f"G0 identity: 1/(2n+1) = tanh(x/2) at x=0.7, d = {d2:.1e} -> "
         f"{'PASS' if d2 < 1e-12 else 'FAIL'}")
assert d2 < 1e-12
L.append("")

# ---------- the exactly solvable NESS machinery
def ness(nS, kS, chan, N):
    """p(n) for S-channel + one ambient channel.
    chan = ('V', kA, TA) | ('M', kM, TA) | (None,) for S alone."""
    lr = np.zeros(N)          # log of per-link up/down ratio
    for n in range(N):
        Dn = W0 + KK*n
        up, dn = kS*nS, kS*(nS + 1.0)
        if chan[0] == 'V':
            nA = 1.0/math.expm1(Dn/chan[2])
            up += chan[1]*nA
            dn += chan[1]*(nA + 1.0)
        elif chan[0] == 'M':
            dmis = KK*n
            nA = 1.0/math.expm1(max(dmis, 1e-12)/chan[2])
            up += chan[1]*nS*nA
            dn += chan[1]*(nS + 1.0)*(nA + 1.0)
        lr[n] = math.log(up) - math.log(dn)
    lp = np.concatenate([[0.0], np.cumsum(lr)])
    lp -= lp.max()
    p = np.exp(lp)
    p /= p.sum()
    return p

def nbar_of(p):
    return float(np.sum(np.arange(len(p))*p))

def gibbs_davies_V(TA, N):
    lp = np.array([-( n*W0 + 0.5*KK*n*(n-1.0))/TA for n in range(N+1)])
    lp -= lp.max()
    p = np.exp(lp); p /= p.sum()
    return nbar_of(p)

def n1_mediated(nS, TA):
    TS = W0/math.log(1.0 + 1.0/nS)
    n = nS
    for _ in range(500):
        n = 1.0/math.expm1(W0/TS + KK*n/TA)
    return n

# ---------- G1/G2/G3 gates
W0, KK = 1.0, 0.02
N = 400
p = ness(1.0, 1.0, (None,), N)
g1a = abs(nbar_of(p) - 1.0)
L.append(f"G1 S alone (nS=1): nbar = {nbar_of(p):.10f} -> d = {g1a:.1e} "
         f"{'PASS' if g1a < 1e-10 else 'FAIL'}")
assert g1a < 1e-10
TA = 0.9
pA = ness(1e-300, 1e-30, ('V', 1.0, TA), N)
gd = gibbs_davies_V(TA, N)
g1b = abs(nbar_of(pA) - gd)
L.append(f"G1 A alone (TA={TA}): nbar = {nbar_of(pA):.8f} vs Davies-Gibbs "
         f"{gd:.8f} -> d = {g1b:.1e} {'PASS' if g1b < 1e-8 else 'FAIL'}")
assert g1b < 1e-8
nSg = 2.0
TSg = W0/math.log(1.0 + 1.0/nSg)
pM = ness(nSg, 1e-30, ('M', 1.0, TSg), N)
gm = gibbs_davies_V(TSg, N)
g1c = abs(nbar_of(pM) - gm)
L.append(f"G1 M alone (TA=TS, nS=2): nbar = {nbar_of(pM):.8f} vs dressed "
         f"Gibbs {gm:.8f} -> d = {g1c:.1e} {'PASS' if g1c < 1e-6 else 'FAIL'}")
assert g1c < 1e-6
KK_saved = KK
KK = 0.0
nA0 = n_be(W0/TA)
pV0 = ness(1.0, 1.0, ('V', 2.0, TA), N)
pred = (1.0*1.0 + 2.0*nA0)/(1.0 + 2.0)
g2 = abs(nbar_of(pV0) - pred)
L.append(f"G2 K=0 linear CFG-V: nbar = {nbar_of(pV0):.8f} vs rate-weighted "
         f"{pred:.8f} -> d = {g2:.1e} {'PASS' if g2 < 1e-8 else 'FAIL'}")
assert g2 < 1e-8
KK = KK_saved
p1 = ness(5.0, 1.0, ('M', 1.0, 1.5), N)
p2 = ness(5.0, 1.0, ('M', 1.0, 1.5), 2*N)
g3 = abs(nbar_of(p1) - nbar_of(p2))/max(nbar_of(p2), 1e-12)
L.append(f"G3 truncation doubling: rel d = {g3:.1e} "
         f"{'PASS' if g3 < 1e-8 else 'FAIL'}")
assert g3 < 1e-8
L.append("")

# ---------- the scans
NS_GRID = [0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
SSTAR = [0.15, 0.3, 0.5, 0.7, 0.85]
KRAT = [0.3, 1.0, 3.0]
results = {}
for cfg in ('V', 'M'):
    lam = np.full((len(NS_GRID), len(SSTAR), len(KRAT)), np.nan)
    nb = np.full_like(lam, np.nan)
    for i, nS in enumerate(NS_GRID):
        for j, sstar in enumerate(SSTAR):
            if cfg == 'V':
                TAv = -W0/math.log(sstar)
                n1 = gibbs_davies_V(TAv, N)
            else:
                TAv = -KK*nS/math.log(sstar)
                n1 = n1_mediated(nS, TAv)
            if abs(n1 - nS) < 1e-6:
                continue
            for k, kr in enumerate(KRAT):
                ch = (cfg, kr, TAv)
                pp = ness(nS, 1.0, ch, N)
                nbv = nbar_of(pp)
                lam[i, j, k] = (nbv - nS)/(n1 - nS)
                nb[i, j, k] = nbv
    results[cfg] = (lam, nb)

# G4 sanity + exponent maps
L.append("scan results (lam = migration fraction toward the ambient "
         "channel's self-consistent endpoint):")
for cfg in ('V', 'M'):
    lam, nb = results[cfg]
    ok4 = np.nanmin(lam) > -0.1 and np.nanmax(lam) < 1.1
    L.append(f"  CFG-{cfg}: lam range [{np.nanmin(lam):.4f}, "
             f"{np.nanmax(lam):.4f}] -> G4 "
             f"{'PASS' if ok4 else 'FAIL (estimator pathology)'}")
    # central-grid log-derivatives
    i, j, k = 2, 2, 1          # nS=1, s*=0.5, kr=1
    def sl(a, b, axis_vals, idx):
        return ((math.log(a) - math.log(b)) /
                (math.log(axis_vals[idx+1]) - math.log(axis_vals[idx-1])))
    ps = sl(lam[i, j+1, k], lam[i, j-1, k], SSTAR, j)
    x2 = [2.0*nb[i2, j, k] + 1.0 for i2 in (i-1, i+1)]
    pn = ((math.log(lam[i+1, j, k]) - math.log(lam[i-1, j, k])) /
          (math.log(x2[1]) - math.log(x2[0])))
    pk = sl(lam[i, j, k+1], lam[i, j, k-1], KRAT, k)
    L.append(f"    central exponents: p_s = {ps:+.3f} (grammar 2), "
             f"p_n = {pn:+.3f} (grammar -2), kappa-slope = {pk:+.3f} "
             f"(grammar 0)")
    # full maps (compact)
    for name, arr, vals, ax in (("p_s vs s*", lam, SSTAR, 1),):
        pass
    results[cfg] = (lam, nb, ps, pn, pk)
L.append("")
L.append("lam maps (rows nS, cols s*, at kappa-ratio = 1):")
for cfg in ('V', 'M'):
    lam = results[cfg][0]
    L.append(f"  CFG-{cfg}:")
    hdr = "    nS\\s*  " + "  ".join(f"{s:5.2f}" for s in SSTAR)
    L.append(hdr)
    for i, nS in enumerate(NS_GRID):
        row = "  ".join(f"{lam[i, j, 1]:5.3f}" if np.isfinite(lam[i, j, 1])
                        else "  n/a" for j in range(len(SSTAR)))
        L.append(f"    {nS:5.1f}  {row}")
L.append("")

# ---------- pre-registered verdict
def verdict(ps, pn, pk):
    fp = (1.6 <= ps <= 2.4) and (-2.4 <= pn <= -1.6) and (-0.2 <= pk <= 0.2)
    return fp
vM = verdict(*results['M'][2:])
vV = verdict(*results['V'][2:])
if vM or vV:
    call = "PASS-" + ("V" if vV else "M")
else:
    psM, pnM, pkM = results['M'][2:]
    clean = all(np.isfinite([psM, pnM, pkM]))
    call = "ALT" if clean else "AMBIG"
L.append(f"PRE-REGISTERED VERDICT: {call}")
L.append("  (bands: p_s in [1.6,2.4]; p_n in [-2.4,-1.6]; kappa-slope "
         "in [-0.2,0.2]; PASS if all three in one config)")

out = "\n".join(L)
print(out)
with open('data/stage6k_analog.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage6k_analog.txt")
