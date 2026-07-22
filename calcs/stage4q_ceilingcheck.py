"""4Q addendum: apply the exact perspective correction to the ceiling-band pairs
(gamma >= 75 deg, s >= 6 kAU) and report band membership before/after."""
import numpy as np
from astropy.io import fits

with fits.open('data/edr3_binaries.fits.gz', memmap=False) as hdul:
    d = hdul[1].data; cols = list(hdul[1].columns.names)
def col(*names, required=True):
    for n in names:
        if n in cols: return np.asarray(d[n], dtype=np.float64)
    if required: raise KeyError(names)
    return None
plx1, plx2 = col('parallax1'), col('parallax2')
eplx1, eplx2 = col('parallax_error1'), col('parallax_error2')
pmra1, pmra2 = col('pmra1'), col('pmra2'); pmde1, pmde2 = col('pmdec1'), col('pmdec2')
epma1, epma2 = col('pmra_error1'), col('pmra_error2')
epmd1, epmd2 = col('pmdec_error1'), col('pmdec_error2')
G1m, G2m = col('phot_g_mean_mag1'), col('phot_g_mean_mag2')
sep = col('sep_AU'); Rch = col('R_chance_align')
ra1, ra2, de1, de2 = col('ra1'), col('ra2'), col('dec1'), col('dec2')
rv1 = col('radial_velocity1','dr2_radial_velocity1', required=False)
rv2 = col('radial_velocity2','dr2_radial_velocity2', required=False)
erv1 = col('radial_velocity_error1','dr2_radial_velocity_error1', required=False)
erv2 = col('radial_velocity_error2','dr2_radial_velocity_error2', required=False)

plx = 0.5*(plx1+plx2)
ok = (Rch < 0.01) & (plx1 > 5) & (plx2 > 5) \
   & (plx1/np.maximum(eplx1,1e-6) > 20) & (plx2/np.maximum(eplx2,1e-6) > 20) \
   & (np.abs(plx1-plx2) < 3*np.hypot(eplx1, eplx2)) & (sep > 200) & (sep < 50000)
MG1 = G1m + 5*np.log10(np.maximum(plx1,1e-6)) - 10
MG2 = G2m + 5*np.log10(np.maximum(plx2,1e-6)) - 10
ok &= (MG1 > 2.6) & (MG1 < 14.2) & (MG2 > 2.6) & (MG2 < 14.2)
ok &= 4.74047/plx*np.sqrt(epma1**2+epmd1**2+epma2**2+epmd2**2) < 0.03

MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,0.37,0.23,0.162,0.102])
Mtot = np.interp(MG1[ok], MG_T, MS_T) + np.interp(MG2[ok], MG_T, MS_T)
plxo = plx[ok]
cosd = np.cos(np.deg2rad(0.5*(de1[ok]+de2[ok])))
sx = (ra2[ok]-ra1[ok])*cosd; sy = de2[ok]-de1[ok]
nrm = np.hypot(sx, sy); sx, sy = sx/nrm, sy/nrm
K = 4.74047/plxo
dvx = K*(pmra2[ok]-pmra1[ok]); dvy = K*(pmde2[ok]-pmde1[ok])
s_kau = sep[ok]/1e3
vc = 0.9417*np.sqrt(Mtot/s_kau)
theta = sep[ok]/(2.06265e8/plxo)
r1, r2 = rv1[ok], rv2[ok]
e1 = erv1[ok] if erv1 is not None else np.full(ok.sum(), 2.0)
e2 = erv2[ok] if erv2 is not None else np.full(ok.sum(), 2.0)
h1, h2 = np.isfinite(r1) & (r1 != 0), np.isfinite(r2) & (r2 != 0)
w1 = np.where(h1, 1/np.maximum(e1,0.5)**2, 0); w2 = np.where(h2, 1/np.maximum(e2,0.5)**2, 0)
rv_sys = (np.where(h1,r1,0)*w1 + np.where(h2,r2,0)*w2)/np.maximum(w1+w2,1e-12)
has = (w1+w2) > 0
print(f"true RV coverage (nonzero, finite): {has.sum()}/{ok.sum()} ({100*has.mean():.0f}%)")

def vtgam(x, y):
    dv = np.hypot(x, y)
    vr = x*sx + y*sy
    gam = np.degrees(np.arccos(np.clip(np.abs(vr)/np.maximum(dv,1e-12), 0, 1)))
    return dv/vc, gam
vt0, g0 = vtgam(dvx, dvy)
dvx_c = dvx + rv_sys*theta*sx; dvy_c = dvy + rv_sys*theta*sy
vt1, g1_ = vtgam(dvx_c, dvy_c)

band0 = (g0 >= 75) & (s_kau >= 6) & (vt0 >= np.sqrt(2)) & (vt0 < 1.67)
band1 = (g1_ >= 75) & (s_kau >= 6) & (vt1 >= np.sqrt(2)) & (vt1 < 1.67)
print(f"ceiling band pairs raw: {band0.sum()}   after correction: {band1.sum()}")
idx = np.where(band0 | band1)[0]
print(f"{'s_kAU':>7} {'RV':>7} {'kappa':>6} {'vt raw':>7} {'vt corr':>8} {'gam raw':>8} {'gam corr':>9} {'in band':>13}")
for i in idx:
    print(f"{s_kau[i]:>7.1f} {rv_sys[i]:>7.1f} {rv_sys[i]*theta[i]/vc[i]:>6.3f} "
          f"{vt0[i]:>7.3f} {vt1[i]:>8.3f} {g0[i]:>8.1f} {g1_[i]:>9.1f} "
          f"{str(bool(band0[i]))+'->'+str(bool(band1[i])):>13}")
above0 = (g0 >= 75) & (s_kau >= 6) & (vt0 >= 1.67) & (vt0 < 2.2)
above1 = (g1_ >= 75) & (s_kau >= 6) & (vt1 >= 1.67) & (vt1 < 2.2)
print(f"above-cliff [1.67,2.2) raw: {above0.sum()}  corrected: {above1.sum()}")
