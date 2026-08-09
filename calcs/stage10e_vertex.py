"""
STAGE 10E -- O5-VERTEX: THE l=2 NEAR-FIELD EXCHANGE CHANNEL (ROUND-22
condition 4 = ROUND-24 condition 7, THE STRIKE-BEARING LEG; author call
2026-08-09 "Go for the next arc close... Hunt continues").

QUESTION. 10B proved gravity's constraint sector gives a bound system
NO off-diagonal soft vertex at l=0 (exterior monopole potential is
interior-constant => H_int ~ H_sys => [H_int, N_s] = 0: dispersive
only), and 10C proved a purely dispersive ambient coupling CANNOT
reproduce the measured lending gate e^{-x} (occupation-blind at
saturation; anti-lending 1/(1+n) perturbatively). The gate's only
ledger-surviving partner is the LOCAL 6Y ambient cloud through the
near-field l=2 (tidal) sector. R22-cond-4: derive that channel as a
quantized Rabi-capable exchange vertex -- "if it comes out O(H) the
mechanism closes for real; if it cannot be Rabi-capable, the 6X/9U
strike fires."

THE DERIVATION (T1 structure exact; T2 magnitude banded):
  T1a  The l=0 selection rule NEGATES at l=2: the interior-regular
       solution of Laplace's equation at l is r^l Y_l -- constant for
       l=0, r^2-STRUCTURED for l=2. An exterior l=2 fluctuation
       couples to the system's internal quadrupole, not its total
       energy. (GE-1, sympy exact.)
  T1b  Both partners are DRESSING modes (matter modes stay dead by
       10B T1/T2): the system's own dressing cloud (mode a, freq om =
       x_loc*T) and the ambient 6Y cloud (mode c, freq Om = x_amb*T;
       T = H/2pi). An aspherical source imprints an l=2 component on
       its dressing cloud PROPORTIONAL to the cloud amplitude =>
       both l=2 amplitudes are LINEAR in their mode coordinates =>
       H_int ~ (a + a')(c + c'): [H_int, N_s] != 0 (GE-2, Fock-exact),
       resonant part a'c + ac' = EXACTLY the 6X toy's exchange
       operator (the toy validated the lending law under this form).
  T1c  The 10B horizon-partner exclusions DO NOT APPLY to the local
       cloud, axis by axis: (gap) its frequency is the ambient
       dressing frequency Om = O(H) soft, not sqrt(l(l+1))H horizon-
       gapped; (geometry) it is CO-LOCATED at r ~ r_M -- no (L/r_c)^l
       lever-arm suppression; (emptiness) its occupation n_amb =
       nu(e_N) - 1 is MEASURED (binary 0.502, galaxy 6.63) -- the
       lending gate P(n >= L) = s^L has support.
  T2   THE INERTIA-CANCELLATION LEMMA (GE-3, sympy exact): writing
       the classical l=2 cross-energy E_12 = eta2*e_s*e_a*
       sqrt(E_c,s*E_c,a) (eta2 = radial overlap <= 1 by Cauchy-
       Schwarz in the l=2 Coulomb inner product; e_i = the l=2/total
       amplitude ratio of cloud i) and quantizing each mode with
       stiffness energy E_c,i (zero-point/classical amplitude ratio
       sqrt(hbar*om_i/(4E_c,i))), the single-quantum exchange
       coupling is
           lam = (eta2/4) * e_s * e_a * sqrt(om*Om)
       -- the E_c CANCEL: the vertex magnitude is INERTIA-BLIND,
       hence INDEPENDENT of the 10C closure assumption (E_c = 2
       hbar om is NOT load-bearing here; no kappa rescaling under
       any 10D outcome).
  T2b  THE CEILING IDENTITY (GE-4b, regression vs the archived 10C G5
       table): at PERFECT geometry (eta2 = e_s = e_a = 1) the vertex
       supplies lam_ceil = (1/2)sqrt(om*Om) = EXACTLY g_close, the
       kappa=1 closure coupling of 10C (same six numbers, digit
       regression). The exchange vertex at unit geometry delivers
       precisely the coupling the closure demands; the measured
       geometric factors ARE the gap.
  T2c  AMPLITUDE INPUTS (declared, banded): e_a = the ambient cloud's
       l=2/l=0 amplitude = the archived 4K/5S EFE quadrupole
       coefficient |q| in [0.073, 0.099], center 0.086, band
       [x1/2, x2] (MG-formulation-derived object; 7G/7H caveat both
       ways -- disclosed ii). e_s binary = point-pair multipole
       ratio eps2 = (s/r_M)^2/4 = 1/(4 x_loc^2) (EXACT within the
       model; capped at 1; the s-running of the channel is DERIVED:
       it dies exactly where the Newtonian regime kills the signal
       anyway). e_s galaxy = disk-asphericity band [0.1, 0.5]
       (banded-not-derived -- disclosed i). eta2 computed for
       declared Gaussian-shell profile families (self-test = 1
       exact; cross-width matrix; convergence gate at 2x nodes).
       Occupation enhancement sqrt((n_s+1) n_a) (thermal means --
       disclosed iv). Cauchy-Schwarz slack x2 carried in lam_max
       (disclosed iii). Near-resonance per 9T unresolvability
       (disclosed v). Single-quantum frame; collective multi-quantum
       amplification NOT assumed (disclosed vii -- it is the GRAY
       successor's question).

GATES (trap #12 -- every letter clause wired):
  GE-1 (STOP-class) interior Laplace: l=0 regular solution constant,
       l=2 regular solution ~ r^2 (sympy exact).
  GE-2 (STOP-class) Fock algebra (6x6 x 6x6): ||[H_int, N_s]|| > 0;
       l=0 contrast commutes exactly; H_int = (a'c + ac') + (a'c' +
       ac) identity exact.
  GE-3 (STOP-class) the cancellation lemma: sympy residual == 0.
  GE-4 anchor table: lam_central, lam_max (all favorable edges),
       Rabi angle per Hubble time, at 6 anchors (binary x_loc = 0.5/
       1.0/2.0 with x_amb = 1.0954, n_a = 0.502; galaxy x_loc =
       0.0953/0.3/1.0 with x_amb = 0.1411, n_a = 6.63) + the Sun
       CONTRAST row (planet-driven eps2 ~ 1e-9: the channel is
       asphericity-gated -- the Cassini-quiet dividend, reading
       grade).
  GE-4b ceiling identity regression: lam_ceil/H reproduces the six
       archived 10C g_close values to 4 digits.
  GE-5 verdict bands (pre-set): V-RABI iff lam_central/H in [0.072,
       1.571] (the 10A gamma-band) at >= 2 non-Sun anchors. V-GRAY
       iff not RABI and lam_max/H >= 1e-3 at >= 1 anchor. V-WEAK iff
       lam_max/H < 1e-3 at ALL non-Sun anchors. V-DIAGONAL iff GE-2
       off-diagonality fails.
  GE-6 Sun suppression factor printed (e_s,sun/e_s,binary).
  GE-7 eta2 convergence: N -> 2N drift < 1% per profile pair
       (one gate per integral family, trap #10).

LETTER GRAMMAR (pre-registered):
  STOP-class: GE-1/2/3 or eta2 self-test fail => abort (bug).
  V-RABI: the carrier exists at the demanded scale => the mechanism's
      exchange leg closes; rise cell (map below).
  V-GRAY: structure EARNED (off-diagonal, Rabi-FORM, exclusions
      inapplicable, ceiling = g_close identity) but magnitude short
      of the band by the measured geometric factors => NO strike at
      stage level; ROUND-25 ADJUDICATION CLAUSE (pre-signed): if the
      round rules the shortfall STRUCTURAL (no live escape axis among
      {e_a re-identification: static-q vs fluctuating-sector
      amplitude; coherent/collective amplification; coherence-time
      re-scoping}), the R22-cond-4 strike FIRES IN ADOPTION (15->8);
      if a live escape axis is identified, HOLD 15 with the successor
      named (O5-GEOMETRY).
  V-WEAK: demonstrated cannot-at-scale even at every favorable edge
      => the strike cell direct (reviewer affirmation required).
  V-DIAGONAL: structural incapacity => strike cell direct.

PRE-SIGNED CREDENCE MAP (arc-level, shared with 10D -- mechanical):
  10E V-RABI + clean round => mech 15 -> 18; with 10D K-ONE => 20.
  10E V-GRAY => HOLD 15 at stage level; round adjudication clause
      above governs adoption.
  10E V-WEAK or V-DIAGONAL (reviewer-affirmed) => 15 -> 8 REGARDLESS
      of 10D. Any unpatched hole blocks that stage's cell.
  anomaly-real 53 UNTOUCHED (no binary-anomaly evidence in play).

Writes data/stage10e_vertex.txt.
"""
import math
import numpy as np
import sympy as sp

OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

say("="*78)
say("STAGE 10E -- O5-VERTEX: the l=2 near-field exchange channel")
say("="*78)

# ---------------- GE-1 interior Laplace ----------------
say("")
say("GE-1 interior Laplace (the l=0 rule and its l=2 negation):")
r = sp.symbols('r', positive=True)
R = sp.Function('R')
sol0 = sp.dsolve(sp.Eq(R(r).diff(r, 2) + 2*R(r).diff(r)/r, 0), R(r)).rhs
sol2 = sp.dsolve(sp.Eq(R(r).diff(r, 2) + 2*R(r).diff(r)/r
                       - 6*R(r)/r**2, 0), R(r)).rhs
say(f"  l=0 general solution: {sol0}  -> regular branch CONSTANT")
say(f"  l=2 general solution: {sol2}  -> regular branch ~ r^2 (NON-constant)")
# explicit check: r^2 solves the l=2 radial equation; const solves l=0
chk0 = sp.simplify((sp.Integer(1)).diff(r, 2) + 2*(sp.Integer(1)).diff(r)/r)
f2 = r**2
chk2 = sp.simplify(f2.diff(r, 2) + 2*f2.diff(r)/r - 6*f2/r**2)
ok1 = (chk0 == 0) and (chk2 == 0)
say(f"  const solves l=0: residual {chk0}; r^2 solves l=2: residual {chk2}"
    f"  -> {'PASS' if ok1 else 'FAIL'}")
say("  => an exterior l=2 fluctuation has INTERIOR STRUCTURE: it couples")
say("     to the internal quadrupole, not the total energy -- the 10B")
say("     dispersive-forcing rule is l=0-SPECIFIC (its own scoping).")

# ---------------- GE-2 Fock algebra ----------------
say("")
say("GE-2 Fock algebra (6x6 x 6x6, numeric-exact):")
N = 6
ad = np.diag(np.sqrt(np.arange(1, N)), -1)
a = ad.T
I6 = np.eye(N)
A, Ad = np.kron(a, I6), np.kron(ad, I6)
C, Cd = np.kron(I6, a), np.kron(I6, ad)
Na = np.kron(ad @ a, I6)
Hint = (A + Ad) @ (C + Cd)
comm = Hint @ Na - Na @ Hint
n_comm = np.max(np.abs(comm))
Hl0 = Na.copy()
comm0 = Hl0 @ Na - Na @ Hl0
n_comm0 = np.max(np.abs(comm0))
res = Ad @ C + A @ Cd
cnt = Ad @ Cd + A @ C
n_split = np.max(np.abs(Hint - (res + cnt)))
ok2 = (n_comm > 0.5) and (n_comm0 == 0.0) and (n_split < 1e-14)
say(f"  ||[H_int, N_s]||_max = {n_comm:.3f} (> 0: OFF-DIAGONAL)")
say(f"  l=0 contrast ||[H_l0, N_s]|| = {n_comm0:.1e} (exactly 0)")
say(f"  H_int = (a'c + ac') + (a'c' + ac) split residual = {n_split:.1e}")
say(f"  -> {'PASS' if ok2 else 'FAIL'}: the resonant part a'c + ac' is")
say("     EXACTLY the 6X exchange operator (the toy's lending law was")
say("     validated under this form; 10C G7 re-ran it verbatim).")

# ---------------- GE-3 the cancellation lemma ----------------
say("")
say("GE-3 the inertia-cancellation lemma (sympy):")
eta2_s, es_s, ea_s, Ecs, Eca, om_s, Om_s, hb = sp.symbols(
    'eta2 e_s e_a E_cs E_ca omega Omega hbar', positive=True)
E12 = eta2_s*es_s*ea_s*sp.sqrt(Ecs*Eca)
zr_s = sp.sqrt(hb*om_s/(4*Ecs))
zr_a = sp.sqrt(hb*Om_s/(4*Eca))
lam_sym = sp.simplify(E12*zr_s*zr_a/hb)
target = eta2_s*es_s*ea_s*sp.sqrt(om_s*Om_s)/4
ok3 = sp.simplify(lam_sym - target) == 0
say(f"  lam = E_12 * (zp/cl)_s * (zp/cl)_a / hbar = {lam_sym}")
say(f"  residual vs (eta2/4) e_s e_a sqrt(om*Om): "
    f"{sp.simplify(lam_sym - target)} -> {'PASS' if ok3 else 'FAIL'}")
say("  => E_c cancels: the vertex is INERTIA-BLIND -- independent of the")
say("     10C closure assumption (not load-bearing here); no kappa")
say("     rescaling under any 10D outcome.")

# ---------------- eta2 overlap (GE-7 convergence) ----------------
say("")
say("eta2 radial overlap (l=2 Coulomb inner product, Gaussian shells")
say("peaked at r_M = 1, widths w; self-test + cross-width matrix):")
def inner_l2(w1, w2, Nn):
    rr = np.linspace(1e-3, 6.0, Nn)
    dr = rr[1] - rr[0]
    rho1 = np.exp(-((rr - 1.0)/w1)**2)
    rho2 = np.exp(-((rr - 1.0)/w2)**2)
    R1, R2 = np.meshgrid(rr, rr, indexing='ij')
    rlt, rgt = np.minimum(R1, R2), np.maximum(R1, R2)
    K = rlt**2/rgt**3
    f1 = rho1*rr*rr
    f2v = rho2*rr*rr
    return float(f1 @ K @ f2v)*dr*dr
def eta2_of(w1, w2, Nn=400):
    x12 = inner_l2(w1, w2, Nn)
    x11 = inner_l2(w1, w1, Nn)
    x22 = inner_l2(w2, w2, Nn)
    return x12/math.sqrt(x11*x22)
ok7 = True
selfv = eta2_of(1.0, 1.0)
say(f"  self-test eta2(1,1) = {selfv:.6f} (exact 1)")
okself = abs(selfv - 1.0) < 1e-9
WS = (0.5, 1.0, 2.0)
etas = {}
for w1 in WS:
    for w2 in WS:
        if w1 < w2:
            e4 = eta2_of(w1, w2, 400)
            e8 = eta2_of(w1, w2, 800)
            drift = abs(e8 - e4)/e4
            ok7 = ok7 and (drift < 0.01)
            etas[(w1, w2)] = e8
            say(f"  eta2({w1},{w2}) = {e8:.4f}  (N->2N drift {drift:.2e})")
eta_min = min(etas.values())
ETA_C = 0.5*(eta_min + 1.0)
say(f"  band: [{eta_min:.3f}, 1.000] (Cauchy-Schwarz cap exact); "
    f"central {ETA_C:.3f}")
say(f"  GE-7 convergence -> {'PASS' if ok7 else 'FAIL'}; self-test -> "
    f"{'PASS' if okself else 'FAIL (STOP)'}")

# ---------------- GE-4 anchors ----------------
say("")
say("GE-4 THE ANCHOR TABLE (lam/H = (eta2/4) e_s e_a sqrt(x_loc*x_amb)/2pi")
say("x occupation enhancement sqrt((n_s+1) n_a)):")
EA_C, EA_LO, EA_HI = 0.086, 0.037, 0.20
GAMMA_LO, GAMMA_HI = 0.072, 1.571
def nbe(x): return 1.0/math.expm1(x)
ANCH = [
    ('binary  x_loc=0.5', 0.5, 1.0954, 0.502, min(1.0, 1.0/(4*0.5**2)), 'pair'),
    ('binary  x_loc=1.0', 1.0, 1.0954, 0.502, 1.0/(4*1.0**2), 'pair'),
    ('binary  x_loc=2.0', 2.0, 1.0954, 0.502, 1.0/(4*2.0**2), 'pair'),
    ('galaxy  x_loc=0.0953', 0.0953, 0.1411, 6.63, 0.3, 'disk'),
    ('galaxy  x_loc=0.3',    0.3,    0.1411, 6.63, 0.3, 'disk'),
    ('galaxy  x_loc=1.0',    1.0,    0.1411, 6.63, 0.3, 'disk'),
]
G_CLOSE_ARCH = [0.0589, 0.0833, 0.1178, 0.0092, 0.0164, 0.0299]
say(f"  {'anchor':<22} {'e_s':>6} {'ceil/H':>8} {'lam_c/H':>9} "
    f"{'lam_max/H':>9} {'Rabi rad/Hub':>12}")
lam_c_list, lam_max_list, ceil_list = [], [], []
for (nm, xl, xa, na, es_c, kind) in ANCH:
    sq = math.sqrt(xl*xa)/(2*math.pi)
    ns = nbe(xl)
    enh = math.sqrt((ns + 1.0)*na)
    lam_c = (ETA_C/4.0)*es_c*EA_C*sq*enh
    es_max = min(1.0, 3.0*es_c) if kind == 'pair' else 0.5
    enh_max = math.sqrt(max((ns + 1.0)*na, ns*(na + 1.0)))
    lam_max = (2.0/4.0)*1.0*es_max*EA_HI*sq*enh_max
    ceil = 0.5*sq
    lam_c_list.append(lam_c); lam_max_list.append(lam_max)
    ceil_list.append(ceil)
    say(f"  {nm:<22} {es_c:>6.3f} {ceil:>8.4f} {lam_c:>9.2e} "
        f"{lam_max:>9.2e} {2*math.pi*lam_c:>12.2e}")
say(f"  (Rabi rad/Hubble = lam/T = 2pi lam/H; saturation of the 6X")
say(f"   lending needs O(pi) -- printed for the honest coherence row)")

say("")
say("GE-4b CEILING IDENTITY regression (vs the archived 10C G5 g_close):")
ok4b = True
for (nm, _, _, _, _, _), cv, av in zip(ANCH, ceil_list, G_CLOSE_ARCH):
    m = abs(cv - av) < 5e-4
    ok4b = ok4b and m
    say(f"  {nm:<22} lam_ceil/H = {cv:.4f}  vs g_close = {av:.4f}  "
        f"{'OK' if m else 'MISMATCH'}")
say(f"  -> {'PASS' if ok4b else 'FAIL'}: lam_ceil = (1/2)sqrt(om*Om) IS")
say("     g_close -- the vertex at PERFECT geometry supplies exactly the")
say("     kappa=1 closure coupling; the measured geometry IS the gap.")

# ---------------- GE-6 the Sun contrast row ----------------
say("")
say("GE-6 the Sun contrast row (asphericity gating, reading grade):")
GMsun = 1.32712e20
A0 = 1.2e-10
AU = 1.495979e11
rM_sun = math.sqrt(GMsun/A0)
eps_J = 9.545e-4*(5.2028*AU/rM_sun)**2
eps_S = 2.858e-4*(9.5826*AU/rM_sun)**2
J2_at_rM = 2.2e-7*(6.96e8/rM_sun)**2
eps_sun = eps_J + eps_S + J2_at_rM
say(f"  r_M(sun) = {rM_sun/AU:.0f} AU (7G: 7030); planet-driven eps2 = "
    f"{eps_sun:.2e}")
say(f"  suppression vs the 10-kAU binary (eps2 = 0.25): "
    f"{0.25/eps_sun:.1e} = {math.log10(0.25/eps_sun):.1f} orders")
say("  => the exchange channel is SOURCE-ASPHERICITY-GATED: strong for")
say("     q ~ 1 binaries and disks, dead ~8 orders for the point-like")
say("     Sun -- the trajectory formulation's Cassini quiet acquires a")
say("     vertex-level echo (reading row, NOT load-bearing).")

# ---------------- GE-5 verdict ----------------
say("")
say("="*78)
nonsun_c = np.array(lam_c_list)
nonsun_m = np.array(lam_max_list)
n_inband = int(((nonsun_c >= GAMMA_LO) & (nonsun_c <= GAMMA_HI)).sum())
any_gray = bool((nonsun_m >= 1e-3).any())
STOP = not (ok1 and ok2 and ok3 and okself)
if STOP:
    letter = 'STOP (bug-class)'
elif n_comm <= 0.5:
    letter = 'V-DIAGONAL'
elif n_inband >= 2:
    letter = 'V-RABI'
elif any_gray:
    letter = 'V-GRAY'
else:
    letter = 'V-WEAK'
say(f"TOKENS: GE1:{'P' if ok1 else 'F'} GE2:{'P' if ok2 else 'F'} "
    f"GE3:{'P' if ok3 else 'F'} GE4b:{'P' if ok4b else 'F'} "
    f"GE7:{'P' if ok7 else 'F'}  inband:{n_inband}/6  "
    f"lam_c[{nonsun_c.min():.1e},{nonsun_c.max():.1e}]H  "
    f"lam_max<= {nonsun_m.max():.1e}H")
say(f"LETTER: {letter}")
say("")
if letter == 'V-GRAY':
    say("V-GRAY earned/short accounting (every clause gated above):")
    say("  EARNED: the off-diagonal l=2 vertex EXISTS (GE-1/GE-2 exact);")
    say("  its resonant part is the 6X exchange operator; the 10B")
    say("  exclusions are inapplicable axis-by-axis (gap/geometry/")
    say("  emptiness); the magnitude is inertia-blind (GE-3, closure-")
    say("  independent); the ceiling identity lam_ceil = g_close (GE-4b).")
    say("  SHORT: central lam ~ 1e-4..1e-3 H vs the demanded band")
    say(f"  [{GAMMA_LO}, {GAMMA_HI}] H; even max-favorable edges reach "
        f"only {nonsun_m.max():.1e} H;")
    say("  the Rabi angle per Hubble time is << pi (no saturation);")
    say("  at the galaxy anchors even PERFECT geometry undershoots the")
    say("  band (ceiling 0.009-0.030 H) -- the frequency product itself")
    say("  caps the channel there.")
    say("  ROUND-25 ADJUDICATION (pre-signed): rule the shortfall")
    say("  STRUCTURAL (strike fires in adoption, 15->8) or identify a")
    say("  live escape axis among {e_a re-identification (static-q vs")
    say("  fluctuating-sector amplitude), coherent/collective")
    say("  amplification, coherence-time re-scoping} (HOLD 15;")
    say("  successor O5-GEOMETRY).")
say("")
say("CREDENCE (pre-signed arc map): V-RABI + clean round -> 15->18")
say("  (with 10D K-ONE -> 20); V-GRAY -> HOLD 15 at stage level, round")
say("  adjudication governs adoption; V-WEAK/V-DIAGONAL reviewer-")
say("  affirmed -> 15->8. anomaly-real 53 UNTOUCHED.")
say("="*78)
say("done")

with open('data/stage10e_vertex.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage10e_vertex.txt")
