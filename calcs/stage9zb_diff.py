"""
STAGE 9Z-b -- O5-DIFF: THE DIFFERENTIAL DISPERSIVE BUDGET (round-20
condition 2 + the condition-3 renormalization specification; the
closure test the reviewer named "the actual O5-LEAK closure test").

QUESTION: 9Z/R20 established that the soft continuum exerts a
dispersive pull Delta ~ -0.02H on every system's mode -- a common-mode
dS-vacuum Lamb shift, excluded-if-real by the deep RAR, so it must
renormalize. What is NOT automatically absorbable is the pull's
x-DEPENDENCE. This stage specifies the renormalization scheme, removes
what the scheme absorbs, propagates the residual through the exact
BE ladder, and tests it against the program's measured bands.

THE RENORMALIZATION SCHEME (condition 3, DECLARED):
  Absorbable transformations -- exactly two, both universal (no
  system tag, no x tag beyond the stated form):
  (S1) one additive bare-frequency constant: the grammar is
       dressed-frequency-native (the 6H JC pull IS the dressing);
       a constant pull redefines where zero sits in the unobservable
       bare scale.
  (S2) one multiplicative x-rescale == an a0 recalibration (x =
       sqrt(g/a0); the measured a0 = 1.05e-10 is the DRESSED scale --
       the 5M temperature lock already lives downstream of S2).
  Everything beyond affine-in-x is OBSERVABLE: it distorts the
  measured ladder shape and/or the two-system split.

THE OBJECT: delta_x(x) = 2 pi Delta(x) with Delta(x) the gapped,
window-excised, UV-complete soft-continuum pull on a system mode at
Omega = x/(2 pi) (H = 1 units; T = 1/(2 pi); dictionary x =
omega/T = sqrt(g_N/a0)). At fiducial coupling g = H the below-window
band is empty (Omega < g) and the vacuum-excess piece has the EXACT
closed form (partial fractions, UV-convergent, no cutoff):

  Delta_vac(x; g) = -(g^2/4 pi^2) * ln(1 + Omega/g)/Omega .

The thermal piece (D_min * 2 n / (Omega - omega), exponentially
convergent) is added numerically with its own convergence gate
(trap #10: one gate per integral family). Small-g co-reads include
the below-window band analytically (partial fractions again).

INSTRUMENTS AND BARS (locked before any number is computed):
  D2 GALAXY LADDER BUDGET: fit delta_x(x) affine (a + b x, = S1 + S2)
     over the measured SPARC range x in [0.14, 1.73]; residual
     delta_x_r(x); propagate EXACTLY: delta_nu(x) = nu(x +
     delta_x_r(x)) - nu(x), nu = 1 + n_BE(x). Deep-window statistic
     dc1_eff = mean of delta_nu over x in [0.14, 0.45] (the deep arm
     reads an additive nu-shift as a c1-shift). BARS: |dc1_eff| <=
     0.05 (the tightest measured c1 band width, the 4Z hier profile
     0.208-0.309) = PASS-grade; |dc1_eff| >= 0.15 (3x, beyond every
     measured band -- the data would have seen it) = FAIL-grade;
     between = GRAY. Co-read (reported, not a bar): max |Delta log10
     g_obs| over the range.
  D3 SPLIT BUDGET: the same universal affine (fit on the galaxy
     range) subtracted at the binary anchor x_bin = 1.095 (binary
     Delta uses the same fiducial-g geometry); the binary gate
     distortion |d gate/gate| = 2 |delta_x_r(x_bin)|. BARS: <= 0.25
     (the 9U/9V gate-precision grade) = PASS; >= 0.75 = FAIL;
     between = GRAY.
  D4 COUPLING CO-READ: dc1_eff(g) over g/H in [0.1, 3] (letter reads
     at fiducial g = H; the scan is a robustness report).
LETTERS: DIFF-CLOSED (D2 PASS and D3 PASS at fiducial) -- the
  dispersive flank of O5-LEAK closes: pull common-mode-absorbable,
  residual inside measurement precision. DIFF-FATAL (D2 FAIL or D3
  FAIL) -- the mechanism predicts a ladder/split distortion the data
  exclude. DIFF-GRAY (anything else).
CREDENCE MAP (pre-signed): DIFF-FATAL -> bath-mechanism conditional
  15 -> 8 immediately. DIFF-CLOSED or DIFF-GRAY -> HOLD 15 pending
  ROUND 21 (rise cell: R21 no-unpatched-hole AND DIFF-CLOSED ->
  15 -> 18; any hole -> hold and execute conditions). anomaly-real
  untouched (no sky number moves).

GATES:
  GDB-0 closed form: sympy verifies d/dOmega and the partial-fraction
        antiderivative reproduce the integrand (exactness of
        Delta_vac); series of nu = 1 + 1/(e^x - 1) reproduces
        1/x + 1/2 + x/12 - x^3/720 (the ladder used for propagation).
  GDB-1 regression to 9Z-addendum: closed-form + thermal at the two
        anchors reproduces the UV-completed Delta = -2.336e-2 (bin) /
        -2.492e-2 (gal) within 0.5%.
  GDB-2 thermal-family convergence: UV 30 vs 60 moves the thermal
        piece < 0.5% of |Delta| (its own gate, per trap #10).
  GDB-3 absorption-fit stability: refitting the affine on the half-
        ranges [0.14, 0.9] / [0.9, 1.73] moves dc1_eff by < 50% of
        its value or < 0.01 absolute (whichever is larger bound).
  GDB-4 ledger leg: mech-9z-leak CURRENT.
DISCLOSURE (sketched pre-commit, stated): the closed form and its
  five-digit match to the addendum anchors (verified by hand
  arithmetic pre-commit); a Taylor sketch suggesting the curvature is
  small at x -> 0 -- but the operative statistic is the affine-fit
  RESIDUAL AT THE RANGE ENDS, amplified by 1/x^2 in the deep arm,
  which the sketch does NOT determine: at delta_x_r ~ 1e-3 the deep
  distortion is ~0.05 = AT the bar. The letter is genuinely open
  between CLOSED and GRAY, with FATAL live if the exact residual is
  larger than the sketch. NOT computed pre-commit: the fit, any
  residual, dc1_eff, the split leg, the scan.
Output: data/stage9zb_diff.txt. Wall-clock: ~1 min.
"""
import csv
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
    with open('data/stage9zb_diff.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

t00 = time.time()
P("=" * 78)
P("STAGE 9Z-b -- O5-DIFF: THE DIFFERENTIAL DISPERSIVE BUDGET")
P("=" * 78)
gates = {}

TWO_PI = 2*math.pi

# ---------------- GDB-0: exactness -----------------------------------
w, Omg, gs, xs = sp.symbols('w Omega g x', positive=True)
anti = (sp.log(w) - sp.log(w - Omg))/Omg          # antiderivative of
integrand = 1/(w*(w - Omg))                        # 1/(w(w-Omega))
ok_pf = sp.simplify(sp.diff(-anti, w) - integrand) == 0
# definite integral over (Omega+g, oo): value = ln(1+Omega/g)/Omega
val = sp.limit(-anti, w, sp.oo) - (-anti.subs(w, Omg + gs))
ok_val = sp.simplify(val - sp.log(1 + Omg/gs)/Omg) == 0
nu_sym = 1 + 1/(sp.exp(xs) - 1)
ser = sp.series(nu_sym, xs, 0, 4).removeO()
ser_ok = sp.simplify(ser - (1/xs + sp.Rational(1, 2) + xs/12
                            - xs**3/720)) == 0
gates['GDB-0'] = ok_pf and ok_val and ser_ok
P("GDB-0 exactness: partial-fraction antiderivative %s; definite "
  "integral = ln(1+Omega/g)/Omega %s; ladder series 1/x + 1/2 + x/12 "
  "- x^3/720 %s -> %s"
  % (ok_pf, ok_val, ser_ok, "PASS" if gates['GDB-0'] else "FAIL"))
if not gates['GDB-0']:
    save(); raise SystemExit(0)

def n_be(wv):
    z = TWO_PI*wv
    if z > 500: return 0.0
    return 1.0/math.expm1(z)

def D_min(wv):
    return (wv**2 + 1.0)/(4*math.pi**2*wv)

def delta_vac(Om_, g_):
    # exact closed form, above-window part; below-window band analytic
    up = -(g_*g_/(4*math.pi**2))*math.log(1 + Om_/g_)/Om_
    lo = 0.0
    # below-window band exists only if Omega - g > gap; at fiducial
    # g = 1 it never does (Omega < 0.3). Included for the g-scan:
    return up, lo

def delta_vac_band(Om_, g_, gap_):
    if Om_ - g_ <= gap_:
        return 0.0
    # + (g^2/4pi^2) * (1/Omega) [ln(w/(Omega-w))] from gap to Omega-g
    f = lambda wv: math.log(wv/(Om_ - wv))
    return (g_*g_/(4*math.pi**2))*(f(Om_ - g_) - f(gap_))/Om_

def delta_thermal(Om_, g_, gap_, UV):
    def f(wv):
        return D_min(wv)*2*n_be(wv)/(Om_ - wv)
    tot = 0.0
    if Om_ - g_ > gap_:
        v1, _ = quad(f, gap_, Om_ - g_, limit=400)
        tot += v1
    v2, _ = quad(f, Om_ + g_, UV, limit=400)
    return (tot + v2)*g_*g_

GAP_GAL = 0.1411/TWO_PI
GAP_BIN = 1.095/TWO_PI
X_BIN = 1.095

def Delta(x_, g_, gap_, UV=30.0):
    Om_ = x_/TWO_PI
    up, _ = delta_vac(Om_, g_)
    band = delta_vac_band(Om_, g_, gap_)
    th = delta_thermal(Om_, g_, gap_, UV)
    return up + band + th

# ---------------- GDB-1 / GDB-2: regression + thermal gate -----------
d_bin = Delta(X_BIN, 1.0, GAP_BIN)
d_gal = Delta(0.22, 1.0, GAP_GAL)
ok_1 = (abs(d_bin - (-2.336e-2)) <= 0.005*2.336e-2
        and abs(d_gal - (-2.492e-2)) <= 0.005*2.492e-2)
gates['GDB-1'] = ok_1
P("GDB-1 anchor regression: Delta_bin = %+.4e (addendum -2.336e-2), "
  "Delta_gal = %+.4e (addendum -2.492e-2) -> %s"
  % (d_bin, d_gal, "PASS" if ok_1 else "FAIL"))
t30 = delta_thermal(0.22/TWO_PI, 1.0, GAP_GAL, 30.0)
t60 = delta_thermal(0.22/TWO_PI, 1.0, GAP_GAL, 60.0)
ok_2 = abs(t60 - t30) <= 0.005*abs(d_gal)
gates['GDB-2'] = ok_2
P("GDB-2 thermal-family convergence (UV 30 vs 60): d = %.1e vs "
  "0.5%% of |Delta| -> %s" % (abs(t60 - t30),
                              "PASS" if ok_2 else "FAIL"))
if not (ok_1 and ok_2):
    P("STOP: regression/convergence failed; nothing downstream quoted")
    save(); raise SystemExit(0)

# ---------------- D2: the galaxy ladder budget ----------------------
def nu(x_):
    if x_ <= 0: return float('inf')
    if x_ > 500: return 1.0
    return 1.0 + 1.0/math.expm1(x_)

XLO, XHI = 0.14, 1.73
DEEP_LO, DEEP_HI = 0.14, 0.45
xg = np.linspace(XLO, XHI, 160)
dx = np.array([TWO_PI*Delta(float(x_), 1.0, GAP_GAL) for x_ in xg])
A = np.vstack([np.ones_like(xg), xg]).T
coef, *_ = np.linalg.lstsq(A, dx, rcond=None)
resid = dx - A @ coef
P("")
P("D2 -- the affine absorption (S1 + S2) on x in [%.2f, %.2f]:"
  % (XLO, XHI))
P("  delta_x(x) = 2 pi Delta(x): range [%+.4f, %+.4f]; affine fit "
  "a = %+.5f, b = %+.5f" % (dx.min(), dx.max(), coef[0], coef[1]))
P("  residual delta_x_r: max |.| = %.2e at x = %.2f; ends %+.2e / "
  "%+.2e" % (np.max(np.abs(resid)), xg[int(np.argmax(np.abs(resid)))],
             resid[0], resid[-1]))
dnu = np.array([nu(float(x_) + float(r_)) - nu(float(x_))
                for x_, r_ in zip(xg, resid)])
mdeep = (xg >= DEEP_LO) & (xg <= DEEP_HI)
dc1_eff = float(np.mean(dnu[mdeep]))
nu_arr = np.array([nu(float(x_)) for x_ in xg])
dex = np.max(np.abs(np.log10((nu_arr + dnu)/nu_arr)))
P("  dc1_eff (deep-window mean of delta_nu, x in [%.2f, %.2f]) = "
  "%+.5f   [bars: PASS <= 0.05, FAIL >= 0.15]"
  % (DEEP_LO, DEEP_HI, dc1_eff))
P("  co-read: max |Delta log10 g_obs| over the range = %.5f dex "
  "(vs 0.13-dex scatter)" % dex)

# GDB-3 stability
d2h = {}
for tag, lo_, hi_ in (('lo-half', 0.14, 0.90), ('hi-half', 0.90, 1.73)):
    m_ = (xg >= lo_) & (xg <= hi_)
    c_, *_ = np.linalg.lstsq(A[m_], dx[m_], rcond=None)
    r_ = dx - A @ c_
    dn_ = np.array([nu(float(x_) + float(rr)) - nu(float(x_))
                    for x_, rr in zip(xg, r_)])
    d2h[tag] = float(np.mean(dn_[mdeep]))
bound = max(0.5*abs(dc1_eff), 0.01)
ok_3 = all(abs(v - dc1_eff) <= bound for v in d2h.values())
gates['GDB-3'] = ok_3
P("  GDB-3 fit stability: half-range dc1_eff = %+.5f / %+.5f "
  "(bound %.3f) -> %s" % (d2h['lo-half'], d2h['hi-half'], bound,
                          "PASS" if ok_3 else "FAIL"))

# ---------------- D3: the split budget ------------------------------
dxb = TWO_PI*Delta(X_BIN, 1.0, GAP_BIN)
rb = dxb - (coef[0] + coef[1]*X_BIN)
gate_dist = 2*abs(rb)
P("")
P("D3 -- the split budget at the binary anchor x = %.3f:" % X_BIN)
P("  delta_x(bin) = %+.5f; universal affine there %+.5f; residual "
  "%+.2e -> |d gate/gate| = %.4f   [bars: PASS <= 0.25, FAIL >= "
  "0.75]" % (dxb, coef[0] + coef[1]*X_BIN, rb, gate_dist))

# ---------------- D4: coupling co-read ------------------------------
P("")
P("D4 -- coupling co-read dc1_eff(g):")
for g_ in (0.1, 0.3, 1.0, 3.0):
    dxg = np.array([TWO_PI*Delta(float(x_), g_, GAP_GAL) for x_ in xg])
    c_, *_ = np.linalg.lstsq(A, dxg, rcond=None)
    r_ = dxg - A @ c_
    dn_ = np.array([nu(float(x_) + float(rr)) - nu(float(x_))
                    for x_, rr in zip(xg, r_)])
    P("  g/H = %-4.1f : dc1_eff = %+.5f   (max |delta_x_r| %.1e)"
      % (g_, float(np.mean(dn_[mdeep])), np.max(np.abs(r_))))

# ---------------- GDB-4 + verdict -----------------------------------
ok_4 = False
for row in csv.reader(open('LEDGER.csv', encoding='utf-8')):
    if row and row[0] == 'mech-9z-leak' and row[1] == 'CURRENT':
        ok_4 = True
gates['GDB-4'] = ok_4
P("")
P("GDB-4 ledger leg mech-9z-leak CURRENT: %s"
  % ("PASS" if ok_4 else "FAIL"))

P("")
d2_pass = abs(dc1_eff) <= 0.05
d2_fail = abs(dc1_eff) >= 0.15
d3_pass = gate_dist <= 0.25
d3_fail = gate_dist >= 0.75
if d2_fail or d3_fail:
    P("==> 9Z-b VERDICT (locked grammar): DIFF-FATAL -- the residual "
      "distortion exceeds what the measured bands allow. Pre-signed "
      "strike: bath-mechanism conditional 15 -> 8.")
elif d2_pass and d3_pass:
    P("==> 9Z-b VERDICT (locked grammar): DIFF-CLOSED -- under the "
      "declared two-parameter universal renormalization (bare "
      "constant + a0 rescale), the soft-continuum dispersive residual "
      "is inside every measured band: the deep-arm distortion "
      "dc1_eff = %+.5f (bar 0.05), the binary gate distortion %.4f "
      "(bar 0.25). The dispersive flank of O5-LEAK closes at "
      "fiducial coupling; with 9Z's leak bound and gap-necessity "
      "this completes the soft-sector accounting at the stated "
      "scope (pending ROUND 21)." % (dc1_eff, gate_dist))
else:
    P("==> 9Z-b VERDICT (locked grammar): DIFF-GRAY -- residuals "
      "between the PASS and FAIL bars; quote everything, no closure "
      "claimed.")
P("    CREDENCE: per the pre-signed map (FATAL -> 15 to 8; otherwise "
  "HOLD 15 pending ROUND 21).")
P("")
P("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
  for k, v in sorted(gates.items())))
P("")
P("done (%.1f min)" % ((time.time()-t00)/60))
save()
print("\nsaved: data/stage9zb_diff.txt")
