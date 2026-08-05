"""STAGE 9E — THE MODEL-LIGHT QUALITY CROSS-CHECK.  Pre-registered
BEFORE the run (same commit as the 9D booking).

Does the anomaly's FIT-FREE core concentrate in the worst-RUWE
quartile (as the fitted alpha does at block grade, 9D), or is it
quality-blind (as the 7I-S median was under the strict cut)?
Statistic: B(q) = median(vt | wide, q) / median(vt | narrow, q),
wide = s in [6, 50) kAU, narrow = s in [0.2, 2) kAU (the narrow
arm ~ Newtonian internal control; quartile common modes cancel to
first order in the ratio); corrected-kernel vt.

Gates: G9E-0 count conservation; G9E-1 rng logged; G9E-2 quartile
edges match 8Z/9A (1.051/1.118/1.231).
Bars (locked): E1 CORE-BLIND iff |D| <= 0.05 OR P(sign flip) >=
0.20, D = B(Q4) - B(Q123).  E2 CORE-CARRIED iff D >= 0.10 AND
P(D <= 0) <= 0.05 (mortal-danger flag; the next review round's
pre-signed map decides).  E3 GRAY-CARRIED else.
NO credence movement (measurement round; pre-stated).
Output: data/stage9e_mlquality.txt
"""
import numpy as np, time
from astropy.io import fits

L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1, plx2 = d['parallax1'], d['parallax2']
eplx1, eplx2 = d['parallax_error1'], d['parallax_error2']
sep, Rch = d['sep_AU'], d['R_chance_align']
G1m, G2m = d['phot_g_mean_mag1'], d['phot_g_mean_mag2']
plx = 0.5*(plx1+plx2)
MG1 = G1m+5*np.log10(np.maximum(plx1, 1e-6))-10
MG2 = G2m+5*np.log10(np.maximum(plx2, 1e-6))-10
sigv = 4.74047/plx*np.sqrt(d['pmra_error1']**2+d['pmdec_error1']**2
                           + d['pmra_error2']**2+d['pmdec_error2']**2)
ok = (Rch < 0.01) & (plx1 > 5) & (plx2 > 5) \
   & (plx1/np.maximum(eplx1, 1e-6) > 20) & (plx2/np.maximum(eplx2, 1e-6) > 20) \
   & (np.abs(plx1-plx2) < 3*np.hypot(eplx1, eplx2)) \
   & (sep > 200) & (sep < 50000) \
   & (MG1 > 2.6) & (MG1 < 14.2) & (MG2 > 2.6) & (MG2 < 14.2) & (sigv < 0.03)
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot_d = np.interp(MG1[ok],MG_T,MS_T)+np.interp(MG2[ok],MG_T,MS_T)
s_d = sep[ok]/1e3
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
sn_ = np.maximum(np.hypot(sx_, sy_), 1e-12)
sux_, suy_ = sx_/sn_, sy_/sn_
_cn = list(d.columns.names)
def _pick(*names):
    for n in names:
        if n in _cn: return np.asarray(d[n], dtype=np.float64)
    raise KeyError(names)
r1_ = _pick('radial_velocity1', 'dr2_radial_velocity1', 'rv1')
r2_ = _pick('radial_velocity2', 'dr2_radial_velocity2', 'rv2')
try:
    er1_ = _pick('radial_velocity_error1', 'dr2_radial_velocity_error1')
    er2_ = _pick('radial_velocity_error2', 'dr2_radial_velocity_error2')
except KeyError:
    er1_ = np.full(len(r1_), 2.0); er2_ = np.full(len(r2_), 2.0)
h1_, h2_ = np.isfinite(r1_), np.isfinite(r2_)
w1_ = np.where(h1_, 1.0/np.maximum(er1_, 0.5)**2, 0.0)
w2_ = np.where(h2_, 1.0/np.maximum(er2_, 0.5)**2, 0.0)
rvs_ = (np.where(h1_, r1_, 0.0)*w1_ + np.where(h2_, r2_, 0.0)*w2_) \
       / np.maximum(w1_+w2_, 1e-12)
th_ = sep/(2.06265e8/plx)
pmcor_ = rvs_*th_*plx/4.74047
vx_ = d['pmra2']-d['pmra1'] + pmcor_*sux_
vy_ = d['pmdec2']-d['pmdec1'] + pmcor_*suy_
vt_d = (4.74047/plx[ok]*np.hypot(vx_[ok], vy_[ok]))/(0.9417*np.sqrt(Mtot_d/s_d))
ruwe = np.maximum(_pick('ruwe1', 'RUWE1'),
                  _pick('ruwe2', 'RUWE2'))[ok]

t0 = time.time()
P("9E THE MODEL-LIGHT QUALITY CROSS-CHECK (pre-reg committed "
  "BEFORE the run; measurement round; NO credence movement)")
qs_ = np.percentile(ruwe, [25, 50, 75])
ok2 = np.allclose(qs_, [1.051, 1.118, 1.231], atol=0.001)
P(f"G9E-2 quartile edges = " + "/".join(f"{v:.3f}" for v in qs_)
  + f" vs 8Z/9A 1.051/1.118/1.231 -> {'PASS' if ok2 else 'FAIL'}")
QI = np.digitize(ruwe, qs_)   # 0..3
WIDE = (s_d >= 6.0) & (s_d < 50.0)
NARR = (s_d >= 0.2) & (s_d < 2.0)
n_tot = 0
for q in range(4):
    nw = int((WIDE & (QI == q)).sum()); nn = int((NARR & (QI == q)).sum())
    n_tot += nw + nn
    P(f"  Q{q+1}: wide {nw}, narrow {nn}")
ok0 = n_tot == int(WIDE.sum() + NARR.sum())
P(f"G9E-0 count conservation -> {'PASS' if ok0 else 'FAIL'}; "
  f"rng = default_rng(71) (G9E-1)")
if not (ok0 and ok2):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    with open('data/stage9e_mlquality.txt', 'w') as f:
        f.write("\n".join(L_)+"\n")
    raise SystemExit(0)
P("GATES: G9E-0/1/2 ALL PASS")
P("")

def B_of(mask_q):
    mw = np.median(vt_d[WIDE & mask_q])
    mn = np.median(vt_d[NARR & mask_q])
    return mw/mn, mw, mn

ball, wall, nall = B_of(np.ones(len(vt_d), dtype=bool))
P(f"[all] B = {ball:.4f} (wide med {wall:.4f} / narrow med "
  f"{nall:.4f})")
bq = []
for q in range(4):
    b_, w_, n_ = B_of(QI == q)
    bq.append(b_)
    P(f"[Q{q+1}] B = {b_:.4f} (wide med {w_:.4f} / narrow med "
      f"{n_:.4f})")
b123, w123, n123 = B_of(QI <= 2)
b4, w4, n4 = B_of(QI == 3)
D0 = b4 - b123
P(f"[Q123] B = {b123:.4f}; [Q4] B = {b4:.4f}; D = B(Q4)-B(Q123) = "
  f"{D0:+.4f}")
P("")

NB = 2000
rng = np.random.default_rng(71)
n_all = len(vt_d)
Ds = np.empty(NB)
for r in range(NB):
    ix = rng.integers(0, n_all, size=n_all)
    v_, s_, q_ = vt_d[ix], s_d[ix], QI[ix]
    w_ = (s_ >= 6.0) & (s_ < 50.0)
    nn_ = (s_ >= 0.2) & (s_ < 2.0)
    m123 = q_ <= 2; m4 = q_ == 3
    try:
        b123_ = (np.median(v_[w_ & m123])/np.median(v_[nn_ & m123]))
        b4_ = (np.median(v_[w_ & m4])/np.median(v_[nn_ & m4]))
        Ds[r] = b4_ - b123_
    except Exception:
        Ds[r] = np.nan
Ds = Ds[np.isfinite(Ds)]
qd = np.percentile(Ds, [5, 50, 95])
p_le0 = float(np.mean(Ds <= 0))
p_flip = float(np.mean(np.sign(Ds) != np.sign(D0))) if D0 != 0 else 1.0
P(f"bootstrap x{len(Ds)}: D pct 5/50/95 = {qd[0]:+.4f}/{qd[1]:+.4f}/"
  f"{qd[2]:+.4f}; P(D <= 0) = {p_le0:.4f}; P(sign flip) = "
  f"{p_flip:.4f}")
P("")
e1 = (abs(D0) <= 0.05) or (p_flip >= 0.20)
e2 = (D0 >= 0.10) and (p_le0 <= 0.05)
if e1:
    P("==> 9E VERDICT (locked grammar): CORE-BLIND - the "
      "model-light boost does not concentrate in the worst-RUWE "
      "quartile; the 9D block-grade collapse reads as an "
      "absorber/power effect; the marginal stratified re-run "
      "arbitrates the fitted channel.")
elif e2:
    P("==> 9E VERDICT (locked grammar): CORE-CARRIED - the "
      "anomaly's model-light core is quality-carried; "
      "MORTAL-DANGER flag; the next review round's pre-signed map "
      "decides.")
else:
    P("==> 9E VERDICT (locked grammar): GRAY-CARRIED - quoted as "
      "the measurement.")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage9e_mlquality.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage9e_mlquality.txt")
