# -*- coding: utf-8 -*-
"""SOL-review pinned comparison: photocentre-cancellation ENGAGEMENT under
rival close-companion populations.

Context (2026-08-24, the GPT-Sol review round): the primary-source audit
found the published analyses (P&S 2023 f_pb two-branch; Banik+24 per the
audit interim) already carry the CANCELLED photocentre wobble law
a_c(q,l) = |q/(1+q) - l/(1+l)|.  Paper 1 Sec 5.3's field-wide uncancelled
attribution is being rewritten.  The surviving, sharper claim is about
POPULATION INPUTS: P&S draw three masses from a broken IMF and sort
(twin-POOR induced q), while our measured subsystem law (7J-z2c GV7,
q-density 1 + t*1[q>=0.9] with t=5 on q in [0.1,1]) is twin-HEAVY
(f_twin = 0.6/1.4 = 42.9%).  Twins are exactly where the cancellation
engages (a_c -> 0).  This script quantifies that: the fraction of naive
(uncancelled) wobble VARIANCE that survives cancellation, per population,
per main-sequence luminosity-ratio exponent l = q^s.

Populations:
  A. P&S 2023 induced q (their Sec 3.2 recipe, verbatim from the audited
     LaTeX source): three masses, p(m) flat to 0.7 Msun then m^-2.35
     (continuous at 0.7), relabel so M2 > M3 among the binary pair,
     M1, M2 >= 0.5, M3 >= 0.01; q = M3/M2.  (Unresolved branch only —
     the resolved branch is uncancelled by construction in their model;
     resolved fraction is separation-dependent and out of scope here.)
  B. Flat q on [0.1, 1] (our null family).
  C. Measured twin-heavy: density (1 + 5*1[q>=0.9])/1.4 on [0.1, 1]
     (7J-z2c winner t=5; beats flat +162 lnL).
  D. El-Badry & Rix 2018-style mild twin excess: flat + 10% extra mass
     at q >= 0.95 (comparison row).

L(M): l = q^s with s in {3.0, 4.0, 5.4}; 5.4 = the exponent implied by
P&S's own G-band relation m = 10^(0.074*(4.69 - M_G))  =>  L ~ m^5.4.

Metrics per (population, s):
  R_var = E[a_c^2] / E[a_u^2]   (survival of naive wobble variance)
  1-R_var = the ENGAGEMENT of the cancellation
  R_amp = E[a_c] / E[a_u]
  f_twin = P(q >= 0.9)
"""
import numpy as np

rng = np.random.default_rng(7250)
N = 2_000_000

def au(q):          # uncancelled mass-only amplitude
    return q/(1.0+q)

def ac(q, s):       # cancelled photocentre amplitude, l = q^s
    l = q**s
    return np.abs(q/(1.0+q) - l/(1.0+l))

# --- population A: P&S 2023 induced q ------------------------------------
def draw_ps_mass(n):
    # p(m) ~ flat on [mlo, 0.7], ~ (m/0.7)^-2.35 above, continuous at 0.7.
    # Upper cut 2.0 Msun (their samples are ~solar; result is cut-insensitive,
    # checked below at 5.0).
    mlo, mk, mhi, alpha = 0.01, 0.7, 2.0, 2.35
    w_flat = mk - mlo
    # integral of (m/mk)^-a dm from mk to mhi
    w_pl = mk*(1.0 - (mhi/mk)**(1.0-alpha))/(alpha-1.0)
    u = rng.random(n)
    m = np.empty(n)
    fl = u < w_flat/(w_flat+w_pl)
    m[fl] = mlo + rng.random(fl.sum())*(mk-mlo)
    v = rng.random((~fl).sum())
    m[~fl] = mk*(1.0 - v*(1.0 - (mhi/mk)**(1.0-alpha)))**(1.0/(1.0-alpha))
    return m

def draw_ps_q(n):
    got = np.empty(0)
    while got.size < n:
        m = draw_ps_mass(3*n).reshape(-1, 3)
        m1, ma, mb = m[:, 0], m[:, 1], m[:, 2]
        m2 = np.maximum(ma, mb); m3 = np.minimum(ma, mb)
        ok = (m1 >= 0.5) & (m2 >= 0.5) & (m3 >= 0.01)
        got = np.concatenate([got, (m3/m2)[ok]])
    return got[:n]

qA = draw_ps_q(N)

# --- population B/C/D on [0.1, 1] ----------------------------------------
qB = 0.1 + 0.9*rng.random(N)

u = rng.random(N)                       # twin t=5: mass 0.9/1.4 flat, 0.5/1.4 in [0.9,1]
tw = u < (0.5/1.4)
qC = np.where(tw, 0.9 + 0.1*rng.random(N), 0.1 + 0.8*rng.random(N))

u = rng.random(N)                       # ER18-style: 10% extra at q>=0.95
tw = u < 0.10
qD = np.where(tw, 0.95 + 0.05*rng.random(N), 0.1 + 0.9*rng.random(N))

# Banik+24 nominal (audit-verified): p(q) ~ q^0.4 smooth + P_eqm = 0.04
# delta-function at q = 1 (their Sec 3.2.4 / App B).  Same [0.1,1] support
# for comparability with B/C/D.
u = rng.random(N)
tw = u < 0.04
v = rng.random(N)                       # inverse-CDF of q^0.4 on [0.1,1]
qE_s = (v*(1.0 - 0.1**1.4) + 0.1**1.4)**(1.0/1.4)
qE = np.where(tw, 1.0, qE_s)

pops = [("P&S sorted-IMF (induced)", qA),
        ("flat [0.1,1]", qB),
        ("measured twin t=5 (7J-z2c)", qC),
        ("ER18-style 10% twin excess", qD),
        ("Banik+24 q^0.4 + 4% twin", qE)]

print("population                      f_twin(q>=0.9)   median q")
for name, q in pops:
    print(f"{name:30s}  {np.mean(q >= 0.9):14.3f}  {np.median(q):9.3f}")
print()
hdr = "population                      s      E[au^2]    E[ac^2]    R_var   engage   R_amp"
print(hdr)
for name, q in pops:
    for s in (3.0, 4.0, 5.4):
        Eu2 = np.mean(au(q)**2); Ec2 = np.mean(ac(q, s)**2)
        Rv = Ec2/Eu2
        Ra = np.mean(ac(q, s))/np.mean(au(q))
        print(f"{name:30s}  {s:3.1f}  {Eu2:9.5f}  {Ec2:9.5f}  {Rv:6.3f}  {1-Rv:6.3f}  {Ra:6.3f}")
    print()

# cross-population contrast at the P&S-implied exponent
s = 5.4
EA = np.mean(ac(qA, s)**2); EC = np.mean(ac(qC, s)**2); EE = np.mean(ac(qE, s)**2)
EuA = np.mean(au(qA)**2);   EuC = np.mean(au(qC)**2);   EuE = np.mean(au(qE)**2)
print(f"At s = 5.4 (P&S's own G-band L(M)):")
print(f"  P&S-population survival R_var    = {EA/EuA:.3f}  (cancellation removes {100*(1-EA/EuA):.1f}% of naive variance)")
print(f"  Banik+24-population survival     = {EE/EuE:.3f}  (removes {100*(1-EE/EuE):.1f}%)")
print(f"  measured-twin-t5 survival        = {EC/EuC:.3f}  (removes {100*(1-EC/EuC):.1f}%)")
print(f"  per-companion wobble variance under the CANCELLED law:")
print(f"  E[ac^2] P&S-pop  / E[ac^2] twin-t5-pop = {EA/EC:.2f}x")
print(f"  E[ac^2] Banik-pop/ E[ac^2] twin-t5-pop = {EE/EC:.2f}x")
print(f"  => the field's simulated companion populations carry {EE/EC:.2f}-{EA/EC:.2f}x the")
print(f"     per-companion wobble variance of the measured twin-heavy population,")
print(f"     even with ALL parties using the correct cancelled law.")

# sensitivity: upper mass cut
def with_cut(mhi):
    global rng
    rng2 = np.random.default_rng(99)
    mlo, mk, alpha = 0.01, 0.7, 2.35
    w_flat = mk - mlo
    w_pl = mk*(1.0 - (mhi/mk)**(1.0-alpha))/(alpha-1.0)
    n = 500_000; got = np.empty(0)
    while got.size < n:
        u = rng2.random(3*n)
        m = np.empty(3*n)
        fl = u < w_flat/(w_flat+w_pl)
        m[fl] = mlo + rng2.random(fl.sum())*(mk-mlo)
        v = rng2.random((~fl).sum())
        m[~fl] = mk*(1.0 - v*(1.0 - (mhi/mk)**(1.0-alpha)))**(1.0/(1.0-alpha))
        m = m.reshape(-1, 3)
        m1, ma, mb = m[:, 0], m[:, 1], m[:, 2]
        m2 = np.maximum(ma, mb); m3 = np.minimum(ma, mb)
        ok = (m1 >= 0.5) & (m2 >= 0.5) & (m3 >= 0.01)
        got = np.concatenate([got, (m3/m2)[ok]])
    q = got[:n]
    return np.mean(ac(q, 5.4)**2)/np.mean(au(q)**2)

print(f"\nsensitivity, P&S R_var at s=5.4 vs upper mass cut: "
      f"mhi=2.0: {with_cut(2.0):.3f}   mhi=5.0: {with_cut(5.0):.3f}")
print("\nGATE GSOL-W1 (sanity): twin q=1, any s: ac = |1/2 - 1/2| =",
      ac(np.array([1.0]), 5.4)[0], "(must be 0.0)")
print("GATE GSOL-W2 (limit): q=0.02, s=5.4: ac/au =",
      f"{(ac(np.array([0.02]),5.4)/au(np.array([0.02])))[0]:.4f}", "(must be ~1: tiny companion, no light)")
