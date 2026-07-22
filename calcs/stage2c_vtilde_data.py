"""
STAGE 2C: the confrontation. EDR3 wide-binary data vs our GPU forward model.
Statistic: v-tilde = v_perp / v_c(s),  v_c = sqrt(G M_tot / s).
Newtonian prediction: distribution of v-tilde independent of separation.
Our modified-law GPU population: median rises to ~1.14 at wide s (boost 1.36).
"""
import gzip, math, sys
import numpy as np
from astropy.io import fits

PATH = 'data/edr3_binaries.fits.gz'

with fits.open(PATH, memmap=False) as hdul:
    d = hdul[1].data
    cols = hdul[1].columns.names
print(f"rows: {len(d)}")
print("columns:", ', '.join(cols))

def col(*names):
    for n in names:
        if n in cols: return np.asarray(d[n], dtype=np.float64)
    raise KeyError(names)

plx1, plx2   = col('parallax1'), col('parallax2')
eplx1, eplx2 = col('parallax_error1'), col('parallax_error2')
pmra1, pmra2 = col('pmra1'), col('pmra2')
pmde1, pmde2 = col('pmdec1'), col('pmdec2')
epma1, epma2 = col('pmra_error1'), col('pmra_error2')
epmd1, epmd2 = col('pmdec_error1'), col('pmdec_error2')
G1, G2       = col('phot_g_mean_mag1'), col('phot_g_mean_mag2')
sep          = col('sep_AU')
Rch          = col('R_chance_align')

# --- quality cuts (Chae-grade, first pass) ---
plx = 0.5*(plx1+plx2)
ok = (Rch < 0.01) & (plx1 > 5) & (plx2 > 5) \
   & (plx1/np.maximum(eplx1,1e-6) > 20) & (plx2/np.maximum(eplx2,1e-6) > 20) \
   & (np.abs(plx1-plx2) < 3*np.hypot(eplx1, eplx2)) \
   & (sep > 200) & (sep < 50000)
# main sequence band via absolute mag
MG1 = G1 + 5*np.log10(np.maximum(plx1,1e-6)) - 10
MG2 = G2 + 5*np.log10(np.maximum(plx2,1e-6)) - 10
ok &= (MG1 > 2.6) & (MG1 < 14.2) & (MG2 > 2.6) & (MG2 < 14.2)
# velocity precision: sigma_vp per pair < 0.03 km/s
sig_vp = 4.74047/plx*np.sqrt(epma1**2+epmd1**2+epma2**2+epmd2**2)
ok &= sig_vp < 0.03
print(f"pairs after cuts: {ok.sum()}")

# --- masses from M_G (Pecaut-Mamajek anchors + bright extension) ---
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
def mass(MG): return np.interp(MG, MG_T, MS_T)
Mtot = mass(MG1[ok]) + mass(MG2[ok])

# --- v-tilde ---
dv = 4.74047/plx[ok]*np.hypot(pmra1[ok]-pmra2[ok], pmde1[ok]-pmde2[ok])  # km/s
s_kau = sep[ok]/1e3
vc = 0.9417*np.sqrt(Mtot/s_kau)                                          # km/s
vt = dv/vc

# --- binned medians + bootstrap ---
bins = [(0.2,0.6),(0.6,2),(2,6),(6,20),(20,50)]
rng = np.random.default_rng(3)
print(f"\n{'s [kAU]':>12} {'N':>7} {'median v~':>10} {'68% CI':>16}")
med0 = None
for b in bins:
    m = (s_kau >= b[0]) & (s_kau < b[1])
    v = vt[m]; n = m.sum()
    if n < 50:
        print(f"{b[0]:>5}-{b[1]:<6} {n:>7} {'--':>10}"); continue
    boots = [np.median(rng.choice(v, n)) for _ in range(400)]
    lo_, hi_ = np.percentile(boots, [16, 84])
    med = np.median(v)
    if med0 is None and b[1] <= 2: med0 = med
    print(f"{b[0]:>5}-{b[1]:<6} {n:>7} {med:>10.3f} {lo_:>7.3f}-{hi_:<7.3f}")

# --- the verdict numbers ---
mnarrow = (s_kau >= 0.2) & (s_kau < 2)
mwide   = (s_kau >= 6)  & (s_kau < 30)
vn, vw = vt[mnarrow], vt[mwide]
bn = [np.median(rng.choice(vn, len(vn))) for _ in range(400)]
bw = [np.median(rng.choice(vw, len(vw))) for _ in range(400)]
ratio = np.median(vw)/np.median(vn)
rs = np.array(bw[:400])/np.array(bn[:400])
print(f"\nboost ratio  median_vt(6-30 kAU) / median_vt(0.2-2 kAU):")
print(f"  DATA:  {ratio:.3f}  (68% CI {np.percentile(rs,16):.3f}-{np.percentile(rs,84):.3f})")
print(f"  Newton prediction:        1.000")
print(f"  our GPU modified law:     ~1.25-1.36 (quadrature-EFE, deep bins)")
print(f"  Chae 2023 claim (AQUAL):  ~1.20 velocity-equivalent")
