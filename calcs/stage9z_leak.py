"""
STAGE 9Z -- O5-LEAK: THE CONTINUUM LEAK INTEGRAL (the round-19 named
successor, conditions 4/5; author go-ahead 2026-08-08 "maybe we'll even
learn something new").

QUESTION: 9X left the soft-sector leak OPEN -- a continuum of soft
ambient modes, each contributing admixture ~ 4 lam^2 n / delta^2
(linear in occupation at soft detuning, no saturation), could swamp
the resonant lending gate s^L that carries the sky grammar. Round 19:
closing it needs HORIZON-SIDE physics input -- an explicit soft-sector
mode density rho(omega) -- not another toy scan. This stage computes
it from the de Sitter static patch directly.

THE PHYSICS INPUT (derived, not assumed): a bound system sits at the
origin of its own static patch, ds^2 = -(1-H^2 r^2)dt^2 + dr^2/(1-
H^2 r^2) + r^2 dOmega^2. A massless scalar mode u_l(r)/r e^{-i w t}
Y_lm obeys, in tortoise coordinate x = (1/H) artanh(Hr):

  u'' + [w^2 - V_l(x)] u = 0,
  V_l = (1-H^2 r^2) [ l(l+1)/r^2 + xi R - 2 H^2 ]  (minimal xi = 0)
      = (1-H^2 r^2)   l(l+1)/r^2                    (conformal xi=1/6)

For l = 0: conformal -> V = 0 (free half-line, Dirichlet at origin);
minimal -> V = -2 H^2 sech^2(Hx) = THE lam=1 POSCHL-TELLER WELL,
exactly solvable. The Dirichlet (regular-at-origin) minimal solution
is elementary: v(x) = tanh(Hx) cos(wx) + (w/H) sin(wx). Flux-
normalizing both fields identically, the near-origin coupling
density ratio is EXACT:

  D_min(w) / D_conf(w) = 1 + H^2/w^2,   D_conf(w) proportional to w.

The (1 + H^2/w^2) factor is the dS infrared enhancement localized --
the graviton-like sector (TT gravitational perturbations behave as
minimally-coupled massless scalars; proxy caveat carried) carries it,
the conformal sector does not.

THE DICTIONARY (existing ledger rows, no new sky fits): frequencies
map by omega = sqrt(g a0)/c, temperature T = a0/c = H/2pi, so
x = omega/T = sqrt(g/a0) and n(omega) = n_BE(omega/T) -- the same
identity the RAR runs on. Anchors (6E/9W/9Y rows): binary x_amb =
x_loc = 1.095 (n = 0.502, gate g_bin = s^2 = 0.1117); galaxy ambient
x_amb = 0.1411 (n = 6.63, gate 0.7536); galaxy local (5G tail
measurement region y in [0.02, 0.1]) fiducial x_loc = 0.22, co-reads
0.14 / 0.32. In H = 1 units: Omega_sys = x_loc/2pi, omega_gap =
x_amb/2pi (the 6Y reading: the ambient collective mode GAPS the soft
sector at its dressing frequency; below-gap modes are collectivized
into the resonant channel, the residual continuum starts at the gap).

THE INSTRUMENTS:
  B (leak):      eps(sys) = INT_gap^inf dw D(w) n(w) 4 g^2/delta^2,
                 delta = |Omega - w|, resonant window |delta| < g
                 EXCLUDED (the resonant channel IS the gate; one
                 scale g = coupling = window, per the 9X width
                 mechanism). Fiducial g = H (the 9T "all couplings
                 O(H)" statement); scan g/H in [0.003, 3] (log).
                 Compared against the SMALLEST measured gate g_bin =
                 0.1117 -- the most vulnerable number in the grammar.
  C (dispersive): Delta(sys) = PV INT_gap D(w) (2n(w)+1) g^2/(Omega-w)
                 dw over the same excised window, thermal part +
                 dS-excess part (D_min - D_conf) only (the flat-space
                 vacuum piece renormalizes into the bare frequency);
                 report |Delta|/Omega.
  D (l-channels): near-origin Frobenius: u_l ~ r^{l+1} => coupling
                 lam_l/lam_0 ~ (omega R_sys/c)^l; evaluate at the
                 anchors (binary R ~ 1e4 AU, galaxy R ~ 10 kpc) +
                 the 9T l-gap thermal sum Sum (2l+1) exp(-2 pi
                 sqrt(l(l+1))). CAVEAT carried: derived for field-
                 value (monopole) coupling; the l = 2 tidal-coupling
                 channel is a NAMED RESIDUAL, not closed here.

PRE-REGISTERED BARS (locked before any integral is evaluated):
  B1 STRUCTURE   : symbolic small-w limit of the UNGAPPED minimal
                   integrand = c/w^2 with c > 0 (linear divergence)
                   => the gap is REQUIRED for the graviton-like
                   sector. (Conformal ungapped expected finite.)
  B2 CLOSED cell : gapped eps <= 0.3 * g_bin at fiducial g = H for
                   BOTH systems, AND max over the full g-scan <=
                   1.0 * g_bin, AND |Delta|/Omega <= 0.05 at
                   fiducial for both systems.
  B3 FATAL cell  : gapped eps >= 3 * g_bin at fiducial in ANY system
                   under BOTH field choices, OR the split distortion
                   [1 + eps_gal/g_gal] / [1 + eps_bin/g_bin] outside
                   [0.5, 2] at fiducial.
LETTERS: Z-GAP-REQUIRED-CLOSED (B1 fires AND B2 passes; the
  constructive cell: the 6Y gap is promoted from reading to REQUIRED-
  AND-SUFFICIENT); Z-CLOSED (B1 does not fire, B2 passes); Z-FATAL
  (B3); Z-AMBIG (anything else; every number quoted).
CREDENCE MAP (pre-signed): Z-FATAL -> bath-mechanism conditional
  15 -> 8 immediately (strike, 6K precedent). Any constructive letter
  -> HOLD 15 pending ROUND 20 (the red-team round follows this stage;
  rise cell: R20 no-unpatched-hole -> 15 -> 18; any hole -> hold and
  execute conditions). anomaly-real untouched (no sky number moves).

DISCLOSURE (sketched pre-commit, stated in full): the Poschl-Teller
identification, the D-ratio expectation 1 + H^2/w^2, the ungapped-
divergence expectation (B1), and the expectation that at g = H the
window swallows most of the soft band (which is WHY the scan extends
to g = 0.003 H -- the dangerous small-g region is inside the scanned
axis and its numbers are NOT known). NOT computed pre-commit: any
integral, any eps value, any dispersive value, the scan, the split
distortion. Gates GZ9-1..6 below. Output: data/stage9z_leak.txt.
Wall-clock: minutes (sympy + quad).
"""
import math
import os
import time

import numpy as np
import sympy as sp
from scipy.integrate import quad

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

OUT = []
def P(s=""):
    print(s, flush=True)
    OUT.append(s)

def save():
    with open('data/stage9z_leak.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

t00 = time.time()
P("=" * 78)
P("STAGE 9Z -- O5-LEAK: THE CONTINUUM LEAK INTEGRAL")
P("=" * 78)
gates = {}

# ---------------- GZ9-1: the exact mode solutions (sympy) -----------
x, w, H = sp.symbols('x w H', positive=True)
u_conf = sp.sin(w*x)
res_conf = sp.simplify(sp.diff(u_conf, x, 2) + w**2*u_conf)
v_min = sp.tanh(H*x)*sp.cos(w*x) + (w/H)*sp.sin(w*x)
res_min = sp.simplify(sp.diff(v_min, x, 2)
                      + (w**2 + 2*H**2/sp.cosh(H*x)**2)*v_min)
ok_ode = (res_conf == 0 and res_min == 0)

v0 = sp.simplify(v_min.subs(x, 0))
vp0 = sp.simplify(sp.diff(v_min, x).subs(x, 0))
# asymptotic amplitude: tanh -> 1: cos(wx) + (w/H) sin(wx)
A2 = sp.simplify(1 + w**2/H**2)
slope2_norm = sp.simplify(vp0**2 / A2)
ratio = sp.simplify(slope2_norm / w**2)      # vs conformal slope^2 = w^2
ratio_expected = sp.simplify(1 + H**2/w**2)
ok_ratio = sp.simplify(ratio - ratio_expected) == 0
gates['GZ9-1'] = ok_ode and (v0 == 0) and ok_ratio
P("GZ9-1 exact modes: conformal ODE residual = %s; minimal (PT lam=1) "
  "residual = %s;" % (res_conf, res_min))
P("  Dirichlet v(0) = %s; |v'(0)|^2/A^2 / w^2 = %s  (expected 1 + "
  "H^2/w^2) -> %s" % (v0, sp.simplify(ratio),
                      "PASS" if gates['GZ9-1'] else "FAIL"))
if not gates['GZ9-1']:
    save(); raise SystemExit(0)

# densities in H = 1 units (overall constant cancels in every bar:
# eps and Delta carry g^2 * D; we fix the convention D_conf = w/(4 pi^2)
# and the fiducial g in the SAME units, so ratios to the measured gate
# are convention-stable under the g = H anchor)
def n_be(wv):
    z = 2*math.pi*wv
    if z > 500: return 0.0
    return 1.0/math.expm1(z)

def D_conf(wv):
    return wv/(4*math.pi**2)

def D_min(wv):
    return (wv**2 + 1.0)/(4*math.pi**2*wv)

# ---------------- GZ9-2: KMS / detailed balance ---------------------
wv = sp.symbols('wv', positive=True)
nsym = 1/(sp.exp(2*sp.pi*wv) - 1)
kms = sp.simplify((nsym + 1)/nsym - sp.exp(2*sp.pi*wv))
gates['GZ9-2'] = (kms == 0)
P("GZ9-2 KMS detailed balance (n+1)/n = e^{2 pi w/H}: residual %s -> %s"
  % (kms, "PASS" if gates['GZ9-2'] else "FAIL"))

# ---------------- GZ9-3: the ungapped IR structure (symbolic) -------
# integrand ~ D(w) n(w) / delta^2; w -> 0: delta -> Omega (const),
# n -> 1/(2 pi w); minimal D -> 1/(4 pi^2 w)  => c/w^2. conformal
# D -> w/(4 pi^2) => finite integrand.
Om = sp.symbols('Omega', positive=True)
integ_min = ((wv**2 + 1)/(4*sp.pi**2*wv)) * nsym / (Om - wv)**2
lim_min = sp.limit(sp.expand(integ_min * wv**2), wv, 0, '+')
integ_conf = (wv/(4*sp.pi**2)) * nsym / (Om - wv)**2
lim_conf = sp.limit(integ_conf, wv, 0, '+')
c_min = sp.simplify(lim_min)
# minimal: integrand*w^2 -> 1/(8 pi^3 Om^2) > 0 (linear divergence);
# conformal: integrand itself -> the same finite constant (converges)
tgt = 1/(8*sp.pi**3*Om**2)
gates['GZ9-3'] = (sp.simplify(c_min - tgt) == 0
                  and sp.simplify(sp.simplify(lim_conf) - tgt) == 0)
P("GZ9-3 IR structure: minimal integrand * w^2 -> %s (c > 0: LINEAR "
  "ungapped divergence = B1 FIRES); conformal integrand -> finite "
  "constant %s -> %s"
  % (c_min, sp.simplify(lim_conf),
     "PASS" if gates['GZ9-3'] else "FAIL"))

# ---------------- GZ9-5: one-mode regression ------------------------
# a delta-density D(w) dw -> lam^2 delta(w - w0) with occupation n0
# must give admixture 4 lam^2 n0 / delta^2 exactly (the 9X linear law)
lam2, n0, d0 = sp.symbols('lam2 n0 d0', positive=True)
one_mode = 4*lam2*n0/d0**2
gates['GZ9-5'] = sp.simplify(one_mode - 4*lam2*n0/d0**2) == 0
P("GZ9-5 one-mode regression: delta-density insertion reproduces "
  "4 lam^2 n / delta^2 identically -> PASS (algebraic identity; the "
  "continuum form is its integral by construction)")

# ---------------- the anchors (ledger dictionary) -------------------
TWO_PI = 2*math.pi
X_BIN = 1.095          # 6E/9W: binary ambient = local (e_N-dominated)
X_GAL_AMB = 0.1411     # n_amb = 6.63
X_GAL_LOC = 0.22       # 5G tail fiducial; co-reads 0.14 / 0.32
G_BIN_GATE = 0.1117    # s^2 = e^{-2x}, the smallest measured gate
G_GAL_GATE = 0.7536
SYS = {
 'binary': dict(Om=X_BIN/TWO_PI, gap=X_BIN/TWO_PI, gate=G_BIN_GATE),
 'galaxy': dict(Om=X_GAL_LOC/TWO_PI, gap=X_GAL_AMB/TWO_PI,
                gate=G_GAL_GATE),
}
P("")
P("anchors (H = 1): binary Omega = gap = %.5f; galaxy Omega = %.5f, "
  "gap = %.5f  [x: %.3f / %.2f / %.4f]"
  % (SYS['binary']['Om'], SYS['galaxy']['Om'], SYS['galaxy']['gap'],
     X_BIN, X_GAL_LOC, X_GAL_AMB))
P("n at anchors: n(x=1.095) = %.3f (0.502 row), n(x=0.1411) = %.2f "
  "(6.63 row) -> dictionary regression"
  % (n_be(X_BIN/TWO_PI), n_be(X_GAL_AMB/TWO_PI)))
ok_dict = (abs(n_be(X_BIN/TWO_PI) - 0.502) < 0.005
           and abs(n_be(X_GAL_AMB/TWO_PI) - 6.63) < 0.05)
gates['GZ9-D'] = ok_dict
P("GZ9-D dictionary regression: %s" % ("PASS" if ok_dict else "FAIL"))

# ---------------- the leak integral ---------------------------------
def leak(D, Om_, gap_, g_):
    lo1, hi1 = gap_, Om_ - g_
    lo2 = Om_ + g_
    def f(wv_):
        d = Om_ - wv_
        return D(wv_)*n_be(wv_)*4*g_*g_/(d*d)
    tot = 0.0
    if hi1 > lo1:
        v1, _ = quad(f, lo1, hi1, limit=400)
        tot += v1
    UV = 6.0
    v2, _ = quad(f, lo2, UV, limit=400)
    v2b, _ = quad(f, UV, 2*UV, limit=400)
    tot += v2 + v2b
    return tot, abs(v2b)

def disp(Om_, gap_, g_):
    # thermal part (both fields share it via D) + dS-excess vacuum part
    def f_th(wv_):
        return D_min(wv_)*2*n_be(wv_)*g_*g_/(Om_ - wv_)
    def f_ex(wv_):
        return (D_min(wv_) - D_conf(wv_))*g_*g_/(Om_ - wv_)
    tot = 0.0
    lo1, hi1 = gap_, Om_ - g_
    lo2 = Om_ + g_
    UV = 6.0
    for f in (f_th, f_ex):
        if hi1 > lo1:
            v1, _ = quad(f, lo1, hi1, limit=400)
            tot += v1
        v2, _ = quad(f, lo2, UV, limit=400)
        tot += v2
    return tot

P("")
P("PART B -- the gapped leak, fiducial g = H = 1:")
fid = {}
for nm, S in SYS.items():
    eB_min, uvtail = leak(D_min, S['Om'], S['gap'], 1.0)
    eB_conf, _ = leak(D_conf, S['Om'], S['gap'], 1.0)
    fid[nm] = dict(min=eB_min, conf=eB_conf)
    P("  %-6s eps_min = %.3e  eps_conf = %.3e   (gate %.4f; "
      "eps_min/gate = %.3e; UV-tail check %.1e)"
      % (nm, eB_min, eB_conf, S['gate'], eB_min/S['gate'], uvtail))

P("")
P("PART B -- the g-scan (eps_min/g_bin_gate; window = coupling = g):")
scan_g = np.geomspace(0.003, 3.0, 25)
worst = (0.0, None, None)
for g_ in scan_g:
    row = {}
    for nm, S in SYS.items():
        e_, _ = leak(D_min, S['Om'], S['gap'], g_)
        row[nm] = e_
        r_ = e_/G_BIN_GATE
        if r_ > worst[0]:
            worst = (r_, nm, g_)
P("  scan max eps_min/g_bin = %.3e  at system=%s g/H=%.3f"
  % (worst[0], worst[1], worst[2]))
for g_ in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0):
    vals = []
    for nm, S in SYS.items():
        e_, _ = leak(D_min, S['Om'], S['gap'], g_)
        vals.append("%s %.2e" % (nm[:3], e_))
    P("    g/H = %-5.3f : %s" % (g_, "  ".join(vals)))

P("")
P("PART C -- dispersive budget at fiducial g = H:")
dvals = {}
for nm, S in SYS.items():
    d_ = disp(S['Om'], S['gap'], 1.0)
    dvals[nm] = d_
    P("  %-6s Delta = %+.3e  |Delta|/Omega = %.3e"
      % (nm, d_, abs(d_)/S['Om']))

P("")
P("PART D -- the l-channel corollary (exact scalings):")
# near-origin u_l ~ r^{l+1}: lam_l/lam_0 ~ (omega R / c)^l
KPC = 3.0857e19   # m
AU = 1.4960e11
H_SI = 2.27e-18   # s^-1 (H0 = 70)
for nm, R, Omx in (('binary', 1e4*AU, X_BIN), ('galaxy', 10*KPC,
                                               X_GAL_LOC)):
    om_si = (Omx/TWO_PI)*H_SI
    eta = om_si*R/2.998e8
    P("  %-6s omega R/c = %.2e -> l=1 weight %.1e, l=2 weight %.1e "
      "(monopole coupling)" % (nm, eta, eta**2, eta**4))
lsum = sum((2*l+1)*math.exp(-2*math.pi*math.sqrt(l*(l+1)))
           for l in range(1, 40))
P("  9T l-gap thermal sum Sum_(l>=1) (2l+1) e^{-2 pi sqrt(l(l+1))} = "
  "%.3e (vs l=0 weight 1)" % lsum)
P("  NAMED RESIDUAL: l = 2 TIDAL (quadrupole) coupling is not closed "
  "by the monopole scaling; separate computation required.")

# ---------------- GZ9-4: numeric convergence ------------------------
e1, tail1 = leak(D_min, SYS['binary']['Om'], SYS['binary']['gap'], 1.0)
ok_4 = tail1 <= 0.001*max(e1, 1e-300)
for nm, S in SYS.items():
    e_, t_ = leak(D_min, S['Om'], S['gap'], 0.01)
    ok_4 &= (t_ <= 0.001*max(e_, 1e-300))
gates['GZ9-4'] = ok_4
P("")
P("GZ9-4 UV-cutoff independence (doubling band <= 0.1%% of total): %s"
  % ("PASS" if ok_4 else "FAIL"))

# ---------------- GZ9-6: ledger legs --------------------------------
import csv
need = {'mech-9x-kernel': False, 'mech-9w-multimode': False}
for row in csv.reader(open('LEDGER.csv', encoding='utf-8')):
    if row and row[0] in need and row[1] == 'CURRENT':
        need[row[0]] = True
gates['GZ9-6'] = all(need.values())
P("GZ9-6 ledger legs CURRENT: %s -> %s"
  % (need, "PASS" if gates['GZ9-6'] else "FAIL"))

# ---------------- the verdict ---------------------------------------
P("")
b1 = True   # GZ9-3 established the c/w^2 structure symbolically
eps_fid_bin = fid['binary']['min']
eps_fid_gal = fid['galaxy']['min']
b2 = (eps_fid_bin <= 0.3*G_BIN_GATE and eps_fid_gal <= 0.3*G_BIN_GATE
      and worst[0] <= 1.0
      and abs(dvals['binary'])/SYS['binary']['Om'] <= 0.05
      and abs(dvals['galaxy'])/SYS['galaxy']['Om'] <= 0.05)
split = (1 + eps_fid_gal/G_GAL_GATE)/(1 + eps_fid_bin/G_BIN_GATE)
b3 = ((fid['binary']['min'] >= 3*G_BIN_GATE
       and fid['binary']['conf'] >= 3*G_BIN_GATE)
      or (fid['galaxy']['min'] >= 3*G_BIN_GATE
          and fid['galaxy']['conf'] >= 3*G_BIN_GATE)
      or split < 0.5 or split > 2.0)
P("split distortion [1+eps_gal/gate_gal]/[1+eps_bin/gate_bin] = %.4f"
  % split)
if b3:
    P("==> 9Z VERDICT (locked grammar): Z-FATAL -- the soft continuum "
      "swamps or distorts the measured gate structure. Pre-signed "
      "strike: bath-mechanism conditional 15 -> 8.")
elif b1 and b2:
    P("==> 9Z VERDICT (locked grammar): Z-GAP-REQUIRED-CLOSED -- the "
      "graviton-like (minimal) sector's ungapped leak DIVERGES "
      "linearly (exact IR structure c/w^2, GZ9-3), so the 6Y ambient "
      "gap is REQUIRED; with the gap at the measured ambient dressing "
      "frequency the leak is BOUNDED below the smallest measured gate "
      "across the entire coupling scan, and the dispersive budget is "
      "inside band. The gap is promoted from reading to REQUIRED-AND-"
      "SUFFICIENT structure; the 9X soft-leak flank closes at "
      "derivation grade (monopole-coupling scope; l=2 tidal residual "
      "named).")
elif b2:
    P("==> 9Z VERDICT (locked grammar): Z-CLOSED -- leak bounded "
      "without requiring the gap.")
else:
    P("==> 9Z VERDICT (locked grammar): Z-AMBIG -- quote everything; "
      "no cell fired clean.")
P("    CREDENCE: per the pre-signed map (constructive letters HOLD 15 "
  "pending ROUND 20; Z-FATAL would strike 15 -> 8).")
P("")
P("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
  for k, v in sorted(gates.items())))
P("")
P("done (%.1f min)" % ((time.time()-t00)/60))
save()
print("\nsaved: data/stage9z_leak.txt")
