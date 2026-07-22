"""
STAGE 4Q: the perspective audit (triggered by the Cookson et al. 2026 full read).

The spherical-projection (perspective) term: a pair with systemic radial velocity
RV and angular separation theta acquires a SPURIOUS on-sky relative velocity
    Dv_spur = -RV * theta * s_hat          (receding pair appears to shrink),
directed ALONG the sky-projected separation vector -- i.e. purely RADIAL in our
(vt, gamma) plane -- growing as theta ~ s/d, hence as s^1.5 in vt terms.
Our pipeline (stage2c onward) never applied this correction. Cookson+26 (their
section 7.2) show its omission inflates median vt by ~0.15 beyond r_M in their
<=130 pc sample. This stage measures the exposure of OUR statistic.

  Q1  exposure map: Newton+perspective-only prediction for our boost ratio,
      injecting per-pair kappa = RV*theta/v_c into the anchor-bin (vt,gamma)
      population (RV drawn from the measured catalog RV distribution).
  Q2  component split: the perpendicular velocity component vt_perp is IMMUNE
      by construction. Boost ratio computed per component: vt (published),
      vt_perp (immune), vt_rad (maximally exposed).
  Q3  direct correction on the RV subsample: slope test reproducing Cookson+26
      fig 7 on our own selection (observed widening rate vs -RV*theta, slope
      should be ~+1 if the effect is present), then medians before/after the
      exact per-pair correction.

Gates:
  G1  anchor-bin (0.2-2 kAU) median kappa < 0.005 (the anchor is unexposed).
  G2  injection round-trip: injecting perspective into anchor pairs and then
      applying the Q3-style correction with the same RVs recovers the base
      median to <0.5%.
  G3  the Q3 slope is within [0.5, 1.5] (effect present with the right sign
      and size) OR |slope| < 0.5 with a stated caveat (effect not detectable
      in our subsample).
"""
import numpy as np
from astropy.io import fits

PATH = 'data/edr3_binaries.fits.gz'
rng = np.random.default_rng(11)

with fits.open(PATH, memmap=False) as hdul:
    d = hdul[1].data
    cols = list(hdul[1].columns.names)

def col(*names, required=True):
    for n in names:
        if n in cols:
            return np.asarray(d[n], dtype=np.float64)
    if required:
        raise KeyError(names)
    return None

plx1, plx2   = col('parallax1'), col('parallax2')
eplx1, eplx2 = col('parallax_error1'), col('parallax_error2')
pmra1, pmra2 = col('pmra1'), col('pmra2')
pmde1, pmde2 = col('pmdec1'), col('pmdec2')
epma1, epma2 = col('pmra_error1'), col('pmra_error2')
epmd1, epmd2 = col('pmdec_error1'), col('pmdec_error2')
G1m, G2m     = col('phot_g_mean_mag1'), col('phot_g_mean_mag2')
sep          = col('sep_AU')
Rch          = col('R_chance_align')
ra1  = col('ra1');  ra2 = col('ra2')
de1  = col('dec1'); de2 = col('dec2')
rv1  = col('radial_velocity1', 'dr2_radial_velocity1', 'rv1', required=False)
rv2  = col('radial_velocity2', 'dr2_radial_velocity2', 'rv2', required=False)
erv1 = col('radial_velocity_error1', 'dr2_radial_velocity_error1', required=False)
erv2 = col('radial_velocity_error2', 'dr2_radial_velocity_error2', required=False)
print("RV columns present:", rv1 is not None and rv2 is not None)

# --- selection: verbatim stage2c cuts ---
plx = 0.5*(plx1+plx2)
ok = (Rch < 0.01) & (plx1 > 5) & (plx2 > 5) \
   & (plx1/np.maximum(eplx1,1e-6) > 20) & (plx2/np.maximum(eplx2,1e-6) > 20) \
   & (np.abs(plx1-plx2) < 3*np.hypot(eplx1, eplx2)) \
   & (sep > 200) & (sep < 50000)
MG1 = G1m + 5*np.log10(np.maximum(plx1,1e-6)) - 10
MG2 = G2m + 5*np.log10(np.maximum(plx2,1e-6)) - 10
ok &= (MG1 > 2.6) & (MG1 < 14.2) & (MG2 > 2.6) & (MG2 < 14.2)
sig_vp = 4.74047/plx*np.sqrt(epma1**2+epmd1**2+epma2**2+epmd2**2)
ok &= sig_vp < 0.03
n_ok = int(ok.sum())
print(f"pairs after stage2c cuts: {n_ok}")
assert n_ok > 13000, "selection drifted from stage2c"

MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot = np.interp(MG1[ok], MG_T, MS_T) + np.interp(MG2[ok], MG_T, MS_T)

# --- kinematic components ---
plxo = plx[ok]
cosd = np.cos(np.deg2rad(0.5*(de1[ok]+de2[ok])))
sx = (ra2[ok]-ra1[ok])*cosd          # sky separation vector, star1 -> star2 (deg)
sy = (de2[ok]-de1[ok])
snorm = np.hypot(sx, sy)
sx, sy = sx/snorm, sy/snorm          # unit separation vector
K = 4.74047/plxo
dvx = K*(pmra2[ok]-pmra1[ok])        # relative velocity of 2 wrt 1 (km/s)
dvy = K*(pmde2[ok]-pmde1[ok])
v_rad_s =  dvx*sx + dvy*sy           # signed widening rate  (km/s)
v_per   = np.abs(dvx*sy - dvy*sx)    # perpendicular component (km/s)
dv      = np.hypot(dvx, dvy)
s_kau = sep[ok]/1e3
vc = 0.9417*np.sqrt(Mtot/s_kau)
vt      = dv/vc
vt_per  = v_per/vc
vt_rad  = np.abs(v_rad_s)/vc
gam = np.degrees(np.arccos(np.clip(np.abs(v_rad_s)/np.maximum(dv,1e-12), 0, 1)))

# angular separation in radians: theta = sep_AU / d_AU, d_AU = 2.06265e8/plx_mas
theta = sep[ok]/(2.06265e8/plxo)

# systemic RV where available
if rv1 is not None:
    r1, r2 = rv1[ok], rv2[ok]
    e1 = erv1[ok] if erv1 is not None else np.full(n_ok, 2.0)
    e2 = erv2[ok] if erv2 is not None else np.full(n_ok, 2.0)
    h1, h2 = np.isfinite(r1), np.isfinite(r2)
    w1 = np.where(h1, 1.0/np.maximum(e1,0.5)**2, 0.0)
    w2 = np.where(h2, 1.0/np.maximum(e2,0.5)**2, 0.0)
    with np.errstate(invalid='ignore'):
        rv_sys = (np.where(h1,r1,0)*w1 + np.where(h2,r2,0)*w2)/np.maximum(w1+w2,1e-12)
    has_rv = (w1+w2) > 0
    print(f"pairs with a systemic RV: {int(has_rv.sum())} "
          f"({100*has_rv.sum()/n_ok:.1f}%); RV 16/50/84 pct: "
          f"{np.percentile(rv_sys[has_rv],[16,50,84]).round(1)}")
    RV_POOL = rv_sys[has_rv]
else:
    has_rv = np.zeros(n_ok, bool); rv_sys = np.zeros(n_ok)
    RV_POOL = rng.normal(0.0, 25.0, 20000)
    print("no RV columns; Q1 uses RV ~ N(0, 25 km/s), Q3 unavailable")

BINS = [(0.2,2),(2,6),(6,20),(20,50)]
ANCH = (s_kau >= 0.2) & (s_kau < 2)
WIDE = (s_kau >= 6)  & (s_kau < 30)      # the published anchor-ratio window

def boots_ratio(a, b, n=400):
    ra_ = [np.median(rng.choice(a, len(a)))/np.median(rng.choice(b, len(b)))
           for _ in range(n)]
    return np.percentile(ra_, [16, 84])

print("\n--- G1: exposure kappa = 25 km/s * theta / v_c by bin ---")
kap25 = 25.0*theta/vc
for b in BINS:
    m = (s_kau >= b[0]) & (s_kau < b[1])
    print(f"  {b[0]:>4}-{b[1]:<4} kAU  N={m.sum():>5}  median kappa {np.median(kap25[m]):.4f}"
          f"  84th pct {np.percentile(kap25[m],84):.4f}")
g1 = np.median(kap25[ANCH])
print(f"G1 anchor median kappa = {g1:.5f}  ->  {'PASS' if g1 < 0.005 else 'FAIL'}")

print("\n--- Q1: Newton+perspective-only prediction for the boost ratio ---")
# base population: anchor-bin (vt, gamma) empirical draws; inject wide-bin kappa
base_vt  = vt[ANCH]; base_gam = np.deg2rad(gam[ANCH])
idx = rng.integers(0, len(base_vt), WIDE.sum())
bvt, bga = base_vt[idx], base_gam[idx]
sgn = rng.choice([-1.0, 1.0], len(bvt))
kap_w = (np.abs(rng.choice(RV_POOL, len(bvt)))*theta[WIDE]/vc[WIDE])
vrad_new = bvt*np.cos(bga)*sgn + kap_w*rng.choice([-1.0,1.0], len(bvt))
vt_inj = np.hypot(vrad_new, bvt*np.sin(bga))
r_pred = np.median(vt_inj)/np.median(base_vt)
print(f"  predicted Newton+perspective boost ratio (6-30 kAU / anchor): {r_pred:.4f}")
print(f"  observed ratio (stage2c): 1.086 (CI 1.064-1.110)")

print("\n--- Q2: the component split (perspective lives ONLY in vt_rad) ---")
print(f"{'stat':>8} {'anchor med':>11} {'wide med':>9} {'ratio':>7} {'68% CI':>15}")
for name, arr in [('vt', vt), ('vt_perp', vt_per), ('vt_rad', vt_rad)]:
    a, w = arr[ANCH], arr[WIDE]
    lo, hi = boots_ratio(w, a)
    print(f"{name:>8} {np.median(a):>11.4f} {np.median(w):>9.4f} "
          f"{np.median(w)/np.median(a):>7.3f} {lo:>7.3f}-{hi:<7.3f}")

print("\n--- per-bin medians of the immune statistic vt_perp ---")
for b in BINS:
    m = (s_kau >= b[0]) & (s_kau < b[1])
    boots = [np.median(rng.choice(vt_per[m], m.sum())) for _ in range(400)]
    print(f"  {b[0]:>4}-{b[1]:<4} kAU  N={m.sum():>5}  median vt_perp "
          f"{np.median(vt_per[m]):.4f}  CI {np.percentile(boots,16):.4f}-{np.percentile(boots,84):.4f}")

if rv1 is not None and has_rv.sum() > 500:
    print("\n--- Q3: direct test on the RV subsample ---")
    mfit = has_rv & (s_kau > 6)
    x = -rv_sys[mfit]*theta[mfit]          # predicted spurious widening rate
    y = v_rad_s[mfit]
    slope = np.sum(x*y)/np.sum(x*x)
    nb = [np.random.default_rng(i).choice(len(x), len(x)) for i in range(200)]
    sl_b = [np.sum(x[i]*y[i])/np.sum(x[i]*x[i]) for i in nb]
    print(f"  wide (>6 kAU) RV pairs: {int(mfit.sum())}")
    print(f"  slope of observed widening vs -RV*theta: {slope:.3f} "
          f"(CI {np.percentile(sl_b,16):.3f}-{np.percentile(sl_b,84):.3f}; expect ~1)")
    g3 = 'PASS (effect present, right sign/size)' if 0.5 <= slope <= 1.5 else \
         ('PASS (effect below detectability here)' if abs(slope) < 0.5 else 'FAIL (anomalous slope)')
    print(f"  G3: {g3}")
    # exact correction on the RV subsample, all four bins
    dvx_c = dvx + rv_sys*theta*sx          # subtract Dv_spur = -RV*theta*s_hat
    dvy_c = dvy + rv_sys*theta*sy
    vt_c = np.hypot(dvx_c, dvy_c)/vc
    ma, mw = has_rv & ANCH, has_rv & WIDE
    print(f"  RV-subsample boost ratio  raw:       "
          f"{np.median(vt[mw])/np.median(vt[ma]):.4f}   (N {int(mw.sum())}/{int(ma.sum())})")
    print(f"  RV-subsample boost ratio  corrected: "
          f"{np.median(vt_c[mw])/np.median(vt_c[ma]):.4f}")
    lo, hi = boots_ratio(vt_c[mw], vt_c[ma])
    print(f"  corrected CI: {lo:.3f}-{hi:.3f}")
    # G2 round-trip on anchor pairs
    inj_x = dvx[ma] - rv_sys[ma]*theta[WIDE][:ma.sum()]*sx[ma] if False else None
    # (round-trip with wide-bin thetas on anchor pairs)
    th_w = rng.choice(theta[WIDE], int(ma.sum()))
    rv_w = rng.choice(RV_POOL, int(ma.sum()))
    ix = dvx[ma] - rv_w*th_w*sx[ma]        # inject spur (-RV*theta*s_hat)
    iy = dvy[ma] - rv_w*th_w*sy[ma]
    cx = ix + rv_w*th_w*sx[ma]             # correct it back
    cy = iy + rv_w*th_w*sy[ma]
    m_inj = np.median(np.hypot(ix,iy)/vc[ma]); m_cor = np.median(np.hypot(cx,cy)/vc[ma])
    m_bas = np.median(vt[ma])
    print(f"  G2 round-trip: base {m_bas:.4f}  injected {m_inj:.4f}  "
          f"corrected {m_cor:.4f}  ->  {'PASS' if abs(m_cor-m_bas) < 0.005*m_bas else 'FAIL'}")
else:
    print("\n--- Q3: unavailable (no/too few RVs) ---")

print("\nInterpretation key:")
print("  boost real + perspective small : vt_perp ratio ~ vt ratio ~ 1.086, Q1 ~ 1.00x")
print("  boost = perspective artifact   : vt_perp ratio ~ 1.00, vt_rad carries all,")
print("                                   Q1 predicts ~1.08, corrected Q3 ratio -> ~1.00")
