"""
STAGE 2D-a: empirical triple diagnostic. RUWE flags unresolved inner companions
(astrometric wobble). If the v-tilde boost is driven by triples, tightening
RUWE should shrink it; if it's gravity, it should persist.
"""
import numpy as np
from astropy.io import fits

d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1, plx2 = d['parallax1'], d['parallax2']
eplx1, eplx2 = d['parallax_error1'], d['parallax_error2']
sep, Rch = d['sep_AU'], d['R_chance_align']
G1, G2 = d['phot_g_mean_mag1'], d['phot_g_mean_mag2']
ruwe1, ruwe2 = d['ruwe1'], d['ruwe2']
plx = 0.5*(plx1+plx2)
MG1 = G1 + 5*np.log10(np.maximum(plx1, 1e-6)) - 10
MG2 = G2 + 5*np.log10(np.maximum(plx2, 1e-6)) - 10
sig_vp = 4.74047/plx*np.sqrt(d['pmra_error1']**2 + d['pmdec_error1']**2
                             + d['pmra_error2']**2 + d['pmdec_error2']**2)
base = (Rch < 0.01) & (plx1 > 5) & (plx2 > 5) \
     & (plx1/np.maximum(eplx1, 1e-6) > 20) & (plx2/np.maximum(eplx2, 1e-6) > 20) \
     & (np.abs(plx1-plx2) < 3*np.hypot(eplx1, eplx2)) \
     & (sep > 200) & (sep < 50000) \
     & (MG1 > 2.6) & (MG1 < 14.2) & (MG2 > 2.6) & (MG2 < 14.2) & (sig_vp < 0.03)

MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
rng = np.random.default_rng(5)

def boost(mask, label):
    Mtot = np.interp(MG1[mask], MG_T, MS_T) + np.interp(MG2[mask], MG_T, MS_T)
    s_kau = sep[mask]/1e3
    dv = 4.74047/plx[mask]*np.hypot(d['pmra1'][mask]-d['pmra2'][mask],
                                    d['pmdec1'][mask]-d['pmdec2'][mask])
    vt = dv/(0.9417*np.sqrt(Mtot/s_kau))
    mn = (s_kau >= 0.2) & (s_kau < 2); mw = (s_kau >= 6) & (s_kau < 30)
    vn, vw = vt[mn], vt[mw]
    rs = [np.median(rng.choice(vw, len(vw)))/np.median(rng.choice(vn, len(vn)))
          for _ in range(400)]
    print(f"{label:>22}: N={mask.sum():>6} (wide {mw.sum():>5})  "
          f"boost = {np.median(vw)/np.median(vn):.3f} "
          f"(68% {np.percentile(rs,16):.3f}-{np.percentile(rs,84):.3f})")

boost(base, "baseline (no RUWE cut)")
boost(base & (ruwe1 < 1.4) & (ruwe2 < 1.4), "RUWE < 1.4")
boost(base & (ruwe1 < 1.2) & (ruwe2 < 1.2), "RUWE < 1.2")
boost(base & (ruwe1 < 1.1) & (ruwe2 < 1.1), "RUWE < 1.1 (strict)")
