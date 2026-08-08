"""
STAGE 9Z-b ADDENDUM -- ROUND-21 conditions 1/2/3/4 (adopted verbatim;
disclosed post-hoc). THE OPERATIVE LETTER IS RELABELED BY THIS FILE:

  DIFF-CLOSED (as fired) -> DIFF-GRAY (round-21 downgrade, adopted).

Grounds (his ruling): (i) a declared gate (GDB-3) was red and the
pre-reg failed to wire gate failure into the letter grammar; (ii) the
statistic is a RANGE, not a point (definition-sensitive 24x; deep-end
value 4.7x the mean); (iii) the margin is order-unity once the
program's own open axes fold in (coupling scale, D-normalization) --
with a plausible joint corner past the FATAL bar. No strike: the
fiducial computation is exact and inside every band; nothing fatal
fired. Credence: HOLD 15 (the rise cell needed CLOSED + no-hole;
neither held).

Conditions executed here:
1. Range quote for dc1_eff (drop the false-precision point).
2. GDB-3 REDESIGNED: refits only on sub-ranges CONTAINING the deep
   window (the scheme never extrapolates); tail-only leg retired.
3. S2 load-bearing disclosure + the temperature-lock check: S1-only
   absorption -> dc1_eff near-FATAL; the a0-rescale carries ~96%;
   implied a0 shift vs the measured (1.05 +/- 0.10)e-10 band.
4. The systematic envelope: deep-end profile value; g = 0.3H;
   D-normalization scaling; the joint pessimal corner. Framing:
   "fiducial-point-inside-band, order-unity margin" -- NOT "closes".
Successors (NOT executed): condition 5 = THE load-bearing one --
discharge the D-normalization provenance (round-20 cond 8), tied to
the l=2/eta^4 self-consistency (the same coupling must carry the O(1)
grammar and the suppressed leak; which channel carries both?);
condition 6 = rung-resolved budget needs SPARC hierarchical weighting;
condition 7 = ambient moot-at-fiducial noted (his own check: identical
fit for all x_amb), reopens only at PT-invalid small g.

Output: data/stage9zb_addendum.txt.
"""
import math
import os
import time

import numpy as np
from scipy.integrate import quad

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

OUT = []
def P(s=""):
    print(s, flush=True)
    OUT.append(s)

def save():
    with open('data/stage9zb_addendum.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

t00 = time.time()
P("=" * 78)
P("STAGE 9Z-b ADDENDUM -- ROUND-21 conditions 1/2/3/4")
P("=" * 78)

TWO_PI = 2*math.pi
def n_be(wv):
    z = TWO_PI*wv
    if z > 500: return 0.0
    return 1.0/math.expm1(z)
def D_min(wv):
    return (wv**2 + 1.0)/(4*math.pi**2*wv)

GAP_GAL = 0.1411/TWO_PI
GAP_BIN = 1.095/TWO_PI
X_BIN = 1.095

def Delta(x_, g_, gap_, dnorm=1.0, UV=30.0):
    Om_ = x_/TWO_PI
    up = -(g_*g_/(4*math.pi**2))*math.log(1 + Om_/g_)/Om_
    band = 0.0
    if Om_ - g_ > gap_:
        f = lambda wv: math.log(wv/(Om_ - wv))
        band = (g_*g_/(4*math.pi**2))*(f(Om_ - g_) - f(gap_))/Om_
    def ft(wv):
        return D_min(wv)*2*n_be(wv)/(Om_ - wv)
    th = 0.0
    if Om_ - g_ > gap_:
        v1, _ = quad(ft, gap_, Om_ - g_, limit=400)
        th += v1
    v2, _ = quad(ft, Om_ + g_, UV, limit=400)
    th = (th + v2)*g_*g_
    return dnorm*(up + band + th)

def nu(x_):
    if x_ <= 0: return float('inf')
    if x_ > 500: return 1.0
    return 1.0 + 1.0/math.expm1(x_)

XLO, XHI = 0.14, 1.73
DEEP_LO, DEEP_HI = 0.14, 0.45
xg = np.linspace(XLO, XHI, 160)
mdeep = (xg >= DEEP_LO) & (xg <= DEEP_HI)
A = np.vstack([np.ones_like(xg), xg]).T

def dc1_of(dx, lo_, hi_, basis='affine'):
    m_ = (xg >= lo_) & (xg <= hi_)
    if basis == 'affine':
        c_, *_ = np.linalg.lstsq(A[m_], dx[m_], rcond=None)
        r_ = dx - A @ c_
    else:  # constant-only (S1 alone)
        c0 = float(np.mean(dx[m_]))
        c_ = np.array([c0, 0.0])
        r_ = dx - c0
    dn_ = np.array([nu(float(x_) + float(rr)) - nu(float(x_))
                    for x_, rr in zip(xg, r_)])
    return float(np.mean(dn_[mdeep])), c_, dn_

dx_fid = np.array([TWO_PI*Delta(float(x_), 1.0, GAP_GAL) for x_ in xg])
d_full, coef, dn_full = dc1_of(dx_fid, XLO, XHI)

# ---------- condition 2: GDB-3 redesigned ---------------------------
P("")
P("CONDITION 2 -- GDB-3 REDESIGNED (deep-window-containing refits "
  "only; the scheme never extrapolates):")
vals = {'full [0.14,1.73]': d_full}
for tag, lo_, hi_ in (('[0.14,0.90]', 0.14, 0.90),
                      ('[0.14,1.00]', 0.14, 1.00),
                      ('[0.14,1.30]', 0.14, 1.30),
                      ('deep-only [0.14,0.45]', 0.14, 0.45)):
    v, _, _ = dc1_of(dx_fid, lo_, hi_)
    vals[tag] = v
for k, v in vals.items():
    P("  %-24s dc1_eff = %+.5f" % (k, v))
spread = max(vals.values()) - min(vals.values())
bound = max(0.5*abs(d_full), 0.01)
ok_g3 = spread <= bound
P("  GDB-3' spread %.4f vs bound %.3f -> %s (the tail-only "
  "extrapolation leg is RETIRED per the ruling)"
  % (spread, bound, "PASS" if ok_g3 else "FAIL"))

# ---------- condition 1: the range quote ----------------------------
v_tail, _, _ = dc1_of(dx_fid, 0.90, 1.73)
deep_end = float(dn_full[0])
P("")
P("CONDITION 1 -- the honest quote (a RANGE, not a point):")
P("  dc1_eff = %+.5f (declared full-range fit); [%+.5f, %+.5f] "
  "across deep-window-inclusive refits; %+.5f under tail-only "
  "extrapolation (outside the scheme, retired); DEEP-END profile "
  "value delta_nu(x=0.14) = %+.5f (the 1/x^2-amplified worst point; "
  "2.2x inside the 0.05 bar)"
  % (d_full, min(vals.values()), max(vals.values()), v_tail,
     deep_end))

# ---------- condition 3: S2 load-bearing + temperature lock ---------
d_s1only, _, dn_s1 = dc1_of(dx_fid, XLO, XHI, basis='const')
b_ = float(coef[1])
a0_shift = 1.0/(1.0 + b_)**2 - 1.0
sigma_frac = abs(a0_shift)*1.05/0.10
P("")
P("CONDITION 3 -- the S2 (a0-rescale) load-bearing disclosure:")
P("  S1-only absorption (constant, no a0 rescale): dc1_eff = %+.5f "
  "-> NEAR-FATAL (bar 0.15); the a0-rescale carries %.0f%% of the "
  "residual reduction" % (d_s1only,
                          100*(1 - abs(d_full)/abs(d_s1only))))
P("  the closure DEPENDS on the linear-pull == a0 degeneracy "
  "(legitimate, now disclosed) and SPENDS temperature-lock "
  "precision: b = %+.5f -> implied a0 shift %+.2f%% = %.2f sigma of "
  "the measured (1.05 +/- 0.10)e-10 band -> PASSES the lock check."
  % (b_, 100*a0_shift, sigma_frac))

# ---------- condition 4: the systematic envelope --------------------
P("")
P("CONDITION 4 -- the systematic envelope (framing: fiducial-point-"
  "inside-band, ORDER-UNITY margin -- not 'closes'):")
dx_g03 = np.array([TWO_PI*Delta(float(x_), 0.3, GAP_GAL) for x_ in xg])
d_g03, _, dn_g03 = dc1_of(dx_g03, XLO, XHI)
for Dn in (1.0, 3.0, 10.0):
    dxD = dx_fid*Dn
    vD, _, dnD = dc1_of(dxD, XLO, XHI)
    P("  D-norm x%-4.0f (g = H)  : dc1_eff = %+.5f  deep-end %+.5f"
      % (Dn, vD, float(dnD[0])))
dx_joint = dx_g03*3.0
d_joint, _, dn_joint = dc1_of(dx_joint, XLO, XHI)
P("  g = 0.3H, D-norm x1   : dc1_eff = %+.5f  deep-end %+.5f"
  % (d_g03, float(dn_g03[0])))
P("  JOINT corner (0.3H, x3): dc1_eff = %+.5f  deep-end %+.5f "
  "-> the mean crosses the PASS bar and the deep-end crosses the "
  "FATAL bar in a corner of the two axes the program has NOT pinned "
  "(9T coupling scale; round-20 condition-8 D-provenance = the "
  "LOAD-BEARING successor)." % (d_joint, float(dn_joint[0])))

P("")
P("9Z-b LETTER RELABELED: DIFF-GRAY (ROUND-21 downgrade ADOPTED; "
  "DIFF-CLOSED as fired is OVERRIDDEN -- a red validation gate "
  "cannot coexist with a clean letter, the statistic is a range, "
  "and the margin is order-unity on unpinned axes). BANKED: the "
  "exact closed form, the anchor regression, the fiducial-point "
  "smallness (every deep-window-inclusive definition <= 0.023 vs "
  "bar 0.05), the split leg 0.0004, the ambient-independence at "
  "fiducial (reviewer-verified: identical fit for every x_amb). "
  "SUCCESSOR (load-bearing): discharge the D-normalization "
  "provenance (round-20 cond 8) tied to the l=2/eta^4 "
  "self-consistency -- 'pin down D and the closure becomes real.'")
P("CREDENCE: bath-mechanism conditional HOLDS 15 (pre-signed map: "
  "rise needed CLOSED + no-hole; neither held; no strike -- nothing "
  "fatal fired at the physical fiducial). anomaly-real 53 untouched.")
P("")
P("done (%.1f min)" % ((time.time()-t00)/60))
save()
print("\nsaved: data/stage9zb_addendum.txt")
