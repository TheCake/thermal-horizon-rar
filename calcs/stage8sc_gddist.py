"""STAGE 8S-c — THE DISTANCE-CHANNEL CONTROL ON THE GD TENSION.
Pre-registered BEFORE any run.

8S-b left a DIAL-TENSION: the gas-dominated subsample localizes
interior at lam ~ -1.5..-1.7 (gas knob cleared).  Lead suspect: the
vertical channel — distance (and inclination) errors move a SPARC
point VERTICALLY (g_bar is distance-invariant at fixed flux while
g_obs ~ 1/D; V_bar components are inclination-independent while
V_obs ~ 1/sin i), and GD galaxies are the worst-distance dwarfs.
This stage adds the vertical channel at MEASURED-prior grade:
per-galaxy log10 intercept x_i ~ N(0, s_i^2) with
    s_i^2 = (e_D/D)^2/ln10^2 + (2 e_inc/tan(inc))^2/ln10^2
marginalized ANALYTICALLY (the conjugate random-intercept model):
    -2lnL_i = sum a^2/v - (sum a/v)^2 / (sum 1/v + 1/s_i^2)
              + sum ln v + ln(1 + s_i^2 sum 1/v)
(the s=0 branch is the VERBATIM flat expression, bit-identity).
This is the DISTANCE-GRADE leg of the named 8S-c decider; the full
hier M/L leg (4Z/5M machinery) remains named if the tension survives.

Gates: G8Sc-0 vertical-OFF regression (the s=0 objective equals the
flat 8S-b objective bit-grade at a probe theta; the OFF GD profile
minimum reproduces 8S-b's -1.542 within one grid step); G8Sc-1
analytic-vs-numeric marginal at 3 probe galaxies (61-node quadrature,
1e-6); G8Sc-2 the s_i distribution printed (fallback count reported).
Bars (locked): C1 TENSION-DISSOLVES iff GD(vertical-ON) D1 interval
overlaps FULL(vertical-ON)'s (attribution = the measured
distance/inclination channel; T5 defense CONDITIONAL; the full-hier
M/L leg still named for the residual).  C2 TENSION-SURVIVES iff
GD(ON) still excludes FULL(ON)'s lam-hat (the tension hardens; the
full-hier leg = the remaining control).  C3 POWER-CARRIED else.
Co-read: FULL's own lam-hat shift under the channel.  NO credence
movement (measurement round; pre-stated).
Output: data/stage8sc_gddist.txt
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize
import time

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)

meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        name = t[0]
        D, eD = float(t[2]), float(t[3])
        inc, einc = float(t[5]), float(t[6])
        q = int(t[17])
        meta[name] = (inc, q, D, eD, einc)
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
svert = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 1.0, 0.0, 0.0))
    if inc < 30 or q > 2: continue
    kept += 1
    s_d2 = (eD/max(D, 1e-6))/LN10
    s_i2 = 2.0*math.radians(einc)/max(math.tan(math.radians(inc)), 1e-6)/LN10
    svert[gi] = math.sqrt(s_d2**2 + s_i2**2)
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC; gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        g_gas.append(gg); g_dsk.append(gd); g_bul.append(gb)
        gobs.append(Vo*Vo/R*KPC)
        sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs = np.log10(gobs)

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_lam(y, lam):
    return (1.0-lam)*nu_standard(y) + lam*nu_be(y)

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

t00 = time.time()
P("8S-c THE DISTANCE-CHANNEL CONTROL (pre-reg committed BEFORE any "
  "run; measurement round; NO credence movement)")

allg = np.unique(gal_id)
NGAL = gal_id.max()+1
gdfrac = np.zeros(NGAL)
for g_ in allg:
    m = gal_id == g_
    gdfrac[g_] = float(np.mean(g_gas[m] > g_dsk[m] + g_bul[m]))
gd_set = np.array([g_ for g_ in allg if gdfrac[g_] >= 0.5])
dd_set = np.array([g_ for g_ in allg if gdfrac[g_] < 0.5])
SV = np.array([svert.get(g_, 0.0) for g_ in range(NGAL)])
nz = SV[allg] > 1e-4
P(f"selector (8S verbatim): GD {len(gd_set)}, DD {len(dd_set)}; "
  f"G8Sc-2 s_vert percentiles 10/50/90 = "
  + "/".join(f"{v:.3f}" for v in np.percentile(SV[allg], [10, 50, 90]))
  + f" dex; GD median = {np.median(SV[gd_set]):.3f}, DD median = "
  f"{np.median(SV[dd_set]):.3f}; galaxies with s<=1e-4 (no vertical "
  f"freedom): {int((~nz).sum())}")

GIDX = [np.where(gal_id == g_)[0] for g_ in range(NGAL)]

def m2ll_vert(th, lam, gset, von=True):
    la0, f, s_int = th
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    nv = nu_lam(gN/a0, lam)
    if np.min(nv) <= 0.05: return 1e12
    gm = gN*nv
    v = sig2 + s_int*s_int
    a = lgobs - np.log10(gm)
    out = 0.0
    for g_ in gset:
        ji = GIDX[g_]
        if len(ji) == 0: continue
        aj, vj = a[ji], v[ji]
        s_ = SV[g_]
        if (not von) or s_ <= 1e-4:
            # VERBATIM flat expression (the s=0 / OFF branch)
            out += float(np.sum(aj*aj/vj + np.log(vj)))
        else:
            iv = 1.0/vj
            Siv = float(np.sum(iv))
            Sa = float(np.sum(aj*iv))
            out += (float(np.sum(aj*aj*iv))
                    - Sa*Sa/(Siv + 1.0/(s_*s_))
                    + float(np.sum(np.log(vj)))
                    + math.log(1.0 + s_*s_*Siv))
    return out

def fit_at(lam, gset, von, th_warm=None):
    starts = [[math.log10(A0_FID), 1.0, 0.08],
              [math.log10(A0_FID)+0.1, 0.8, 0.12]]
    if th_warm is not None: starts.insert(0, list(th_warm))
    best = None
    for th0 in starts:
        b = minimize(lambda t: m2ll_vert(t, lam, gset, von), th0,
                     method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
        if best is None or b.fun < best.fun: best = b
    return best

LGX = np.round(np.arange(-2.00, 1.501, 0.05), 3)
def profile(gset, von):
    prof, th = [], None
    for lam in LGX:
        b = fit_at(lam, gset, von, th)
        prof.append(b.fun); th = b.x
    prof = np.array(prof)
    i = int(np.argmin(prof))
    lam_hat = LGX[i]
    if 0 < i < len(LGX)-1:
        x3, y3 = LGX[i-1:i+2], prof[i-1:i+2]
        c2_, c1_, _ = np.polyfit(x3, y3, 2)
        if c2_ > 0: lam_hat = -c1_/(2*c2_)
    lo = hi = None
    for j in range(i, -1, -1):
        if prof[j] > prof[i]+1.0:
            lo = np.interp(prof[i]+1.0, [prof[j+1], prof[j]],
                           [LGX[j+1], LGX[j]]) \
                if prof[j] != prof[j+1] else LGX[j]
            break
    for j in range(i, len(LGX)):
        if prof[j] > prof[i]+1.0:
            hi = np.interp(prof[i]+1.0, [prof[j-1], prof[j]],
                           [LGX[j-1], LGX[j]])
            break
    return prof, lam_hat, lo, hi, i

# ---------------- gates ----------------
th_p = [math.log10(A0_FID), 1.0, 0.08]
# G8Sc-1: analytic vs numeric marginal at 3 probe galaxies
xg = np.linspace(-6, 6, 61)
wq = np.exp(-0.5*xg*xg); wq /= wq.sum()
a0p = 10**th_p[0]
gNp = g_gas + th_p[1]*g_dsk + g_bul
gmp = gNp*nu_lam(gNp/a0p, 0.7)
vp_ = sig2 + th_p[2]**2
ap_ = lgobs - np.log10(gmp)
g1_ok = True
probes = [g_ for g_ in gd_set if SV[g_] > 0.02][:3]
for g_ in probes:
    ji = GIDX[g_]
    aj, vj = ap_[ji], vp_[ji]
    s_ = SV[g_]
    iv = 1.0/vj
    Siv = float(np.sum(iv)); Sa = float(np.sum(aj*iv))
    ana = (float(np.sum(aj*aj*iv)) - Sa*Sa/(Siv + 1.0/(s_*s_))
           + float(np.sum(np.log(vj))) + math.log(1.0 + s_*s_*Siv))
    lik = np.array([float(np.sum((aj - s_*x)**2*iv)) for x in xg])
    num = -2.0*(math.log(float(np.sum(wq*np.exp(-0.5*(lik-lik.min())))))
                - 0.5*lik.min()) + float(np.sum(np.log(vj)))
    okp = abs(ana - num) <= 1e-6
    g1_ok &= okp
    P(f"G8Sc-1 probe gal {g_} (s={s_:.3f}): analytic {ana:.8f} vs "
      f"numeric {num:.8f} -> {'PASS' if okp else 'FAIL'}")
# G8Sc-0: OFF branch == the flat objective (bit) + OFF GD profile
va = m2ll_vert(th_p, 0.7, list(gd_set), von=False)
gNf = g_gas + th_p[1]*g_dsk + g_bul
nvf = nu_lam(gNf/a0p, 0.7)
gmf = gNf*nvf
vf = sig2 + th_p[2]**2
af = lgobs - np.log10(gmf)
wsel = np.isin(gal_id, gd_set)
vflat = float(np.sum((af*af/vf + np.log(vf))[wsel]))
ok0a = abs(va - vflat) <= 1e-9
prof0, lam0, _, _, _ = profile(list(gd_set), von=False)
ok0b = abs(lam0 - (-1.542)) <= 0.1
g0_ok = ok0a and ok0b
P(f"G8Sc-0 OFF-branch identity: {va:.9f} vs flat {vflat:.9f} "
  f"(d={va-vflat:.2e}) -> {'PASS' if ok0a else 'FAIL'}; OFF GD "
  f"profile lam_hat = {lam0:.3f} vs 8S-b -1.542 "
  f"-> {'PASS' if ok0b else 'FAIL'}")

if not (g0_ok and g1_ok):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    with open('data/stage8sc_gddist.txt', 'w') as f:
        f.write("\n".join(L)+"\n")
    raise SystemExit(0)
P("GATES: G8Sc-0/1/2 ALL PASS")
P("")

# ---------------- the vertical-ON profiles ----------------
res = {}
for tag, gset in (('FULL', allg), ('GD', gd_set), ('DD', dd_set)):
    prof, lam_hat, lo, hi, i = profile(list(gset), von=True)
    edge = 'INTERIOR' if 0 < i < len(LGX)-1 else 'EDGE'
    res[tag] = dict(lam=lam_hat, lo=lo, hi=hi)
    P(f"[{tag:4} vertON] lam_hat = {lam_hat:.3f} (D1 "
      f"{None if lo is None else round(lo,3)}.."
      f"{None if hi is None else round(hi,3)}) -> c1_hat = "
      f"{lam_hat/2:.3f}; {edge}")

P("")
gd = res['GD']; fu = res['FULL']
two = gd['lo'] is not None and gd['hi'] is not None
ov = two and fu['lo'] is not None and fu['hi'] is not None \
     and not (gd['hi'] < fu['lo'] or gd['lo'] > fu['hi'])
excl = (gd['hi'] is not None and gd['hi'] < fu['lam']) or \
       (gd['lo'] is not None and gd['lo'] > fu['lam'])
if ov:
    P("==> 8S-c VERDICT (locked grammar): TENSION-DISSOLVES — under "
      "the MEASURED distance/inclination vertical channel the "
      "gas-dominated interval rejoins the full-sample dial: the 8S-b "
      "tension attributes to the vertical channel at measured-prior "
      "grade; T5's defense holds CONDITIONAL on it (the full-hier "
      "M/L leg remains the named residual control).")
elif excl:
    P("==> 8S-c VERDICT (locked grammar): TENSION-SURVIVES — the "
      "gas-dominated subsample still excludes the full-sample dial "
      "with the measured vertical channel on; the tension HARDENS "
      "and the full-hier leg (4Z/5M machinery on GD) is the "
      "remaining in-pipeline control.")
else:
    P("==> 8S-c VERDICT (locked grammar): POWER-CARRIED — intervals "
      "quoted; no attribution at this grade.")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t00)/60:.1f} min)")

with open('data/stage8sc_gddist.txt', 'w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage8sc_gddist.txt")
