"""
STAGE 7I-MATERIAL: the pre-declared MATERIAL-branch investigation of the
strict-multiplicity collapse (7I part S: alpha_hat 1.1-1.26 -> 0.4-0.75,
Newton +100 -> +4..+14, w_rad 0.2 -> 0.3 in every strict fit).

Competing hypotheses:
  H1 (companion-carried signal): the boost lives disproportionately in the
     flagged pairs -> the model-light median boost should ALSO collapse
     under the cut.
  H2 (instrument mismatch): the v7 model's FIXED noise inflation FPM = 1.5
     (calibrated on the full sample, 3N/3O) overestimates noise on the
     cleaner strict subsample -> fitted alpha absorbs the difference; the
     model-light median boost should SURVIVE the cut.

Instrument (descriptive grade, model-light): the 2C/2D anchor-relative
median statistic — median(vt | wide)/median(vt | anchor) — which contains
no forward model and no FPM. Computed per cut component (RUWE<1.2 both /
overluminous-free / combined) on raw and perspective-corrected velocities,
with 1000-rep pair bootstrap CIs, per-bin retention fractions, and per-bin
median sigma_c (to quantify how much cleaner the strict sample is = the
size of the H2 mismatch). 2D precedent: the median boost was RUWE-stable.
Writes data/stage7i_material.txt.
"""
import numpy as np
from astropy.io import fits

d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1, plx2 = d['parallax1'], d['parallax2']
eplx1, eplx2 = d['parallax_error1'], d['parallax_error2']
sep, Rch = d['sep_AU'], d['R_chance_align']
G1m, G2m = d['phot_g_mean_mag1'], d['phot_g_mean_mag2']
ruwe1, ruwe2 = d['ruwe1'], d['ruwe2']
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
Mtot = np.interp(MG1, MG_T, MS_T)+np.interp(MG2, MG_T, MS_T)
s_kau = sep/1e3
vc = 0.9417*np.sqrt(Mtot/np.maximum(s_kau, 1e-9))

dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx = dra*np.cos(dec_m); sy = d['dec2']-d['dec1']
vx = d['pmra2']-d['pmra1']; vy = d['pmdec2']-d['pmdec1']
vt_raw = (4.74047/plx*np.hypot(vx, vy))/vc
# corrected (6G formulas)
sn_u = np.maximum(np.hypot(sx, sy), 1e-12)
sux, suy = sx/sn_u, sy/sn_u
_cn = list(d.columns.names)
def _pick(*names):
    for n in names:
        if n in _cn:
            return np.asarray(d[n], dtype=np.float64)
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
rvs = (np.where(h1_, r1_, 0.0)*w1_ + np.where(h2_, r2_, 0.0)*w2_) \
      / np.maximum(w1_+w2_, 1e-12)
th = sep/(2.06265e8/plx)
pmcor = rvs*th*plx/4.74047
vt_cor = (4.74047/plx*np.hypot(vx + pmcor*sux, vy + pmcor*suy))/vc
# perpendicular component (4Q mass-immune statistic), raw velocities
vperp = np.abs(-vx*suy + vy*sux)
vt_perp = (4.74047/plx*vperp)/vc

# cut components (identical construction to stage7i_ablations precheck)
c1, c2 = d['bp_rp1'], d['bp_rp2']
colb = np.concatenate([c1[ok], c2[ok]])
mgb = np.concatenate([MG1[ok], MG2[ok]])
goodb = np.isfinite(colb) & np.isfinite(mgb) & (colb > -0.5) & (colb < 6.0)
CB = np.linspace(np.nanpercentile(colb[goodb], 0.5),
                 np.nanpercentile(colb[goodb], 99.5), 41)
cc = 0.5*(CB[:-1]+CB[1:])
ridge = np.full(len(cc), np.nan)
for i in range(len(cc)):
    m = goodb & (colb >= CB[i]) & (colb < CB[i+1])
    if m.sum() > 30:
        ridge[i] = np.median(mgb[m])
vr = np.isfinite(ridge)
d1 = MG1 - np.interp(c1, cc[vr], ridge[vr])
d2 = MG2 - np.interp(c2, cc[vr], ridge[vr])
okc = np.isfinite(c1) & np.isfinite(c2) & (c1 > -0.5) & (c1 < 6.0) \
    & (c2 > -0.5) & (c2 < 6.0)
RU = (ruwe1 < 1.2) & (ruwe2 < 1.2)
OV = okc & (d1 >= -0.4) & (d2 >= -0.4)

BINS = {'anchor 0.2-2': (0.2, 2.0), '2-6': (2.0, 6.0),
        'wide 6-20': (6.0, 20.0), 'wide 20-50': (20.0, 50.0)}
rng = np.random.default_rng(20260725)

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

def ratio(vt, m, lo=6.0, hi=30.0):
    a = vt[m & (s_kau >= 0.2) & (s_kau < 2.0)]
    w = vt[m & (s_kau >= lo) & (s_kau < hi)]
    return np.median(w)/np.median(a), len(a), len(w)

def boot_ratio(vt, m, n=1000, lo=6.0, hi=30.0):
    ia = np.where(m & (s_kau >= 0.2) & (s_kau < 2.0))[0]
    iw = np.where(m & (s_kau >= lo) & (s_kau < hi))[0]
    rs = []
    for _ in range(n):
        rs.append(np.median(vt[rng.choice(iw, len(iw))])
                  / np.median(vt[rng.choice(ia, len(ia))]))
    return np.percentile(rs, [16, 50, 84])

VARIANTS = [('baseline (14071)', ok),
            ('RUWE<1.2 both', ok & RU),
            ('overluminous-free', ok & OV),
            ('strict (RUWE+overlum)', ok & RU & OV)]

P("STAGE 7I-MATERIAL: model-light median boost under the strict-cut "
  "components")
P("(H1 companion-carried -> ratio collapses with the cut; H2 noise-model "
  "mismatch -> ratio survives)")
for name, m in VARIANTS:
    P(f"\n--- {name}: {int(m.sum())} pairs ---")
    for vt, tag in ((vt_raw, 'raw'), (vt_cor, 'corr'), (vt_perp, 'perp')):
        r, na, nw = ratio(vt, m)
        lo_, md_, hi_ = boot_ratio(vt, m)
        P(f"  {tag:>4}: wide(6-30)/anchor median ratio = {r:.4f}  "
          f"boot [{lo_:.4f}, {md_:.4f}, {hi_:.4f}]  (N {na}/{nw})")
    keep = [f"{b}: {int((m & (s_kau >= lo) & (s_kau < hi)).sum())}"
            f"/{int((ok & (s_kau >= lo) & (s_kau < hi)).sum())}"
            for b, (lo, hi) in BINS.items()]
    P("  retention " + "; ".join(keep))
    sig = [f"{b}: {np.median(sigv[m & (s_kau >= lo) & (s_kau < hi)]):.4f}"
           for b, (lo, hi) in BINS.items()]
    P("  median sigma_v " + "; ".join(sig))

P("\ncomparators: 2C boost 1.086 (CI 1.064-1.110); 4Q corrected 1.078 "
  "(1.052-1.103); 4Q vt_perp 1.151 (1.115-1.197); 2D: median boost was "
  "RUWE-stable at 1.4/1.2/1.1 cuts")
with open('data/stage7i_material.txt', 'w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage7i_material.txt")
