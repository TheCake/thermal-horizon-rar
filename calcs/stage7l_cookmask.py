"""Stage 7L — the Cookson-selection PROXY MASK (pre-reg in NOTES
2026-07-27; proxy-grade declared).

Their cuts on our EDR3 catalog: d < 130 pc, plx S/N > 40 both,
RUWE < 1.25 both, both RVs finite with |dRV| < 10 km/s, s in 1-30 kAU,
velocity uncertainty < 0.1 in vtilde units, corrected vtilde < 2.5,
both colors 0.5 < BP-RP < 3.5, no overluminous component (3J ridge
criterion, delta >= -0.4 both).  Saves a full-catalog-length boolean
data/stage7l_cookmask.npy and prints N vs their 1,421.
"""
import numpy as np
from astropy.io import fits

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
   & (sep > 200) & (sep < 50000) & (MG1 > 2.6) & (MG1 < 14.2) \
   & (MG2 > 2.6) & (MG2 < 14.2) & (sigv < 0.03)
print(f"baseline ok sample: {int(ok.sum())}")

# corrected vtilde (4R convention)
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot = np.interp(MG1, MG_T, MS_T)+np.interp(MG2, MG_T, MS_T)
s_kau = sep/1e3
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
sn_ = np.maximum(np.hypot(sx_, sy_), 1e-12)
r1 = np.asarray(d['dr2_radial_velocity1'], dtype=np.float64)
r2 = np.asarray(d['dr2_radial_velocity2'], dtype=np.float64)
er1 = np.asarray(d['dr2_radial_velocity_error1'], dtype=np.float64)
er2 = np.asarray(d['dr2_radial_velocity_error2'], dtype=np.float64)
h1, h2 = np.isfinite(r1), np.isfinite(r2)
w1 = np.where(h1, 1.0/np.maximum(er1, 0.5)**2, 0.0)
w2 = np.where(h2, 1.0/np.maximum(er2, 0.5)**2, 0.0)
rvs = (np.where(h1, r1, 0.0)*w1 + np.where(h2, r2, 0.0)*w2) \
      / np.maximum(w1+w2, 1e-12)
th = sep/(2.06265e8/plx)
pmcor = rvs*th*plx/4.74047
vx = d['pmra2']-d['pmra1'] + pmcor*sx_/sn_
vy = d['pmdec2']-d['pmdec1'] + pmcor*sy_/sn_
vc = 0.9417*np.sqrt(Mtot/np.maximum(s_kau, 1e-9))
vt = (4.74047/plx*np.hypot(vx, vy))/vc

# 3J overluminosity ridge on the baseline sample's components
c1, c2 = np.asarray(d['bp_rp1'], np.float64), np.asarray(d['bp_rp2'], np.float64)
colb = np.concatenate([c1[ok], c2[ok]])
mgb = np.concatenate([MG1[ok], MG2[ok]])
goodb = np.isfinite(colb) & (colb > -0.5) & (colb < 6)
cc = np.linspace(colb[goodb].min(), colb[goodb].max(), 42)
cc = 0.5*(cc[:-1]+cc[1:])
ridge = np.full(len(cc), np.nan)
for i in range(len(cc)):
    m = goodb & (np.abs(colb-cc[i]) <= (cc[1]-cc[0]))
    if m.sum() >= 25:
        ridge[i] = np.median(mgb[m])
vr = np.isfinite(ridge)
d1 = MG1 - np.interp(c1, cc[vr], ridge[vr])
d2 = MG2 - np.interp(c2, cc[vr], ridge[vr])

cook = ok.copy()
steps = [
    ('d < 130 pc', plx > 1000.0/130.0),
    ('plx S/N > 40 both', (plx1/np.maximum(eplx1, 1e-6) > 40)
     & (plx2/np.maximum(eplx2, 1e-6) > 40)),
    ('RUWE < 1.25 both', (np.asarray(d['ruwe1']) < 1.25)
     & (np.asarray(d['ruwe2']) < 1.25)),
    ('both RVs finite', h1 & h2),
    ('|dRV| < 10 km/s', np.abs(np.where(h1, r1, 0)-np.where(h2, r2, 0)) < 10),
    ('s in 1-30 kAU', (s_kau >= 1.0) & (s_kau < 30.0)),
    ('sig_vt < 0.1', (sigv/np.maximum(vc, 1e-9)) < 0.1),
    ('vt_corr < 2.5', vt < 2.5),
    ('0.5 < BP-RP < 3.5 both', (c1 > 0.5) & (c1 < 3.5) & (c2 > 0.5)
     & (c2 < 3.5)),
    ('no overluminous (3J)', (d1 >= -0.4) & (d2 >= -0.4)),
]
for name, m in steps:
    cook = cook & m
    print(f"  + {name:<26} -> {int(cook.sum())}")
np.save('data/stage7l_cookmask.npy', cook)
print(f"COOK proxy sample: {int(cook.sum())} (Cookson et al.: 1,421); "
      f"mask saved (full-catalog length)")
