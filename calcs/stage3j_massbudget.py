"""
STAGE 3J: independent mass-error budget for sigma_m (TODO #2b, measurement
half). The v4 fit selects sigma_m = 0.30 (= 0.5 * sigma_lnMtot -> ~60%
effective total-mass error). Here we measure how much of that photometric
masses can actually supply, from the sample's own photometry:
  1. Main-sequence ridge M_G(BP-RP) from the 2x14,071 component stars;
     per-star offsets delta_i = M_G - ridge(color).
  2. Robust MS width sigma(M_G | color) -> sigma_lnM via the M_G->mass table
     slope (mass is ASSIGNED from M_G, so the M_G spread at fixed color is the
     inversion error).
  3. The binary trick: components share distance/age/metallicity/extinction,
     so corr(delta_1, delta_2) across pairs splits the scatter into a shared
     (fully correlated in M_tot) and a per-star part.
  4. sigma_m(mass) = 0.5 * sigma_lnMtot, compared against the fitted 0.30.
CPU-only. Writes data/stage3j_summary.txt.
"""
import numpy as np
from astropy.io import fits

d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1,plx2 = d['parallax1'],d['parallax2']
eplx1,eplx2 = d['parallax_error1'],d['parallax_error2']
sep,Rch = d['sep_AU'],d['R_chance_align']
G1,G2 = d['phot_g_mean_mag1'],d['phot_g_mean_mag2']
plx = 0.5*(plx1+plx2)
MG1 = G1+5*np.log10(np.maximum(plx1,1e-6))-10
MG2 = G2+5*np.log10(np.maximum(plx2,1e-6))-10
sigv = 4.74047/plx*np.sqrt(d['pmra_error1']**2+d['pmdec_error1']**2
                           +d['pmra_error2']**2+d['pmdec_error2']**2)
ok = (Rch<0.01)&(plx1>5)&(plx2>5)&(plx1/np.maximum(eplx1,1e-6)>20) \
   &(plx2/np.maximum(eplx2,1e-6)>20)&(np.abs(plx1-plx2)<3*np.hypot(eplx1,eplx2)) \
   &(sep>200)&(sep<50000)&(MG1>2.6)&(MG1<14.2)&(MG2>2.6)&(MG2<14.2)&(sigv<0.03)

MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
M1 = np.interp(MG1[ok], MG_T, MS_T)
M2 = np.interp(MG2[ok], MG_T, MS_T)

c1 = d['bp_rp1'][ok]; c2 = d['bp_rp2'][ok]
mg1 = MG1[ok]; mg2 = MG2[ok]

col = np.concatenate([c1, c2])
mg  = np.concatenate([mg1, mg2])
# catalog encodes missing photometry as 1e20 sentinels, not NaN
good = np.isfinite(col) & np.isfinite(mg) & (col > -0.5) & (col < 6.0)
L = []
def P(s):
    print(s); L.append(s)

P(f"STAGE 3J mass budget: {ok.sum()} pairs, {good.sum()}/{len(col)} stars "
  f"with BP-RP")

# 1. ridge: running median M_G in color bins
CB = np.linspace(np.nanpercentile(col[good],0.5),
                 np.nanpercentile(col[good],99.5), 41)
cc = 0.5*(CB[:-1]+CB[1:])
ridge = np.full(len(cc), np.nan)
for i in range(len(cc)):
    m = good & (col>=CB[i]) & (col<CB[i+1])
    if m.sum() > 30: ridge[i] = np.median(mg[m])
vr = np.isfinite(ridge)
def ridge_f(c): return np.interp(c, cc[vr], ridge[vr])

delta = mg - ridge_f(col)   # <0 = overluminous

# 2. robust width per color bin (MAD; clips the binary sequence)
widths, wcounts = [], []
for i in range(len(cc)):
    m = good & (col>=CB[i]) & (col<CB[i+1])
    if m.sum() > 100:
        w = 1.4826*np.median(np.abs(delta[m]-np.median(delta[m])))
        widths.append(w); wcounts.append(m.sum())
widths = np.array(widths); wcounts = np.array(wcounts)
sig_MG = np.average(widths, weights=wcounts)
P(f"MS robust width sigma(M_G|color): {sig_MG:.3f} mag "
  f"(bin range {widths.min():.3f}-{widths.max():.3f})")
f_over = np.mean(delta[good] < -0.4)
P(f"overluminous fraction (delta < -0.4 mag, unresolved-companion proxy): "
  f"{f_over:.3f}")

# 3. binary trick: shared vs per-star decomposition
d1 = mg1 - ridge_f(c1); d2 = mg2 - ridge_f(c2)
gp = np.isfinite(d1) & np.isfinite(d2) \
   & (c1 > -0.5) & (c1 < 6.0) & (c2 > -0.5) & (c2 < 6.0)
# robust correlation: clip to +/-3*sigma to keep the binary sequence from
# dominating the covariance
cl = 3*sig_MG
m3 = gp & (np.abs(d1)<cl) & (np.abs(d2)<cl)
rho = np.corrcoef(d1[m3], d2[m3])[0,1]
P(f"pair correlation corr(delta1, delta2): {rho:.3f} "
  f"(clipped +/-{cl:.2f} mag, n={m3.sum()})")
P(f"  -> shared (correlated-in-pair) variance fraction: {max(rho,0):.2f}; "
  f"per-star fraction: {1-max(rho,0):.2f}")

# 4. convert to sigma_lnM and sigma_m
lnM_slope = np.gradient(np.log(MS_T), MG_T)     # dlnM/dM_G (negative)
sl1 = np.abs(np.interp(mg1, MG_T, lnM_slope))
sl2 = np.abs(np.interp(mg2, MG_T, lnM_slope))
sig_lnM1 = sig_MG*sl1; sig_lnM2 = sig_MG*sl2
w1 = M1/(M1+M2); w2 = M2/(M1+M2)
rr = max(rho, 0)
var_tot = (w1*sig_lnM1)**2 + (w2*sig_lnM2)**2 + 2*rr*w1*w2*sig_lnM1*sig_lnM2
sig_lnMtot = np.sqrt(var_tot)
sm_mass = 0.5*sig_lnMtot
P(f"per-star sigma_lnM: median {np.median(np.r_[sig_lnM1,sig_lnM2]):.3f} "
  f"(16-84%: [{np.percentile(np.r_[sig_lnM1,sig_lnM2],16):.3f}, "
  f"{np.percentile(np.r_[sig_lnM1,sig_lnM2],84):.3f}])")
P(f"sigma_lnMtot per pair: median {np.median(sig_lnMtot):.3f}")
P(f"==> sigma_m(mass) = 0.5*sigma_lnMtot: median {np.median(sm_mass):.3f} "
  f"(16-84%: [{np.percentile(sm_mass,16):.3f}, {np.percentile(sm_mass,84):.3f}])")
fit_sm = 0.30
med = np.median(sm_mass)
rest = np.sqrt(max(fit_sm**2 - med**2, 0))
P(f"fitted sigma_m = {fit_sm}; mass errors supply {med:.3f} "
  f"-> unexplained residual sqrt({fit_sm}^2 - {med:.3f}^2) = {rest:.3f}")
P(f"(residual candidates: unresolved companions beyond the {f_over:.1%} "
  f"overluminous tail, e-model mismatch, intra-bin selection; extinction is "
  f"small at plx>5 mas but included in the shared part)")

with open('data/stage3j_summary.txt','w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage3j_summary.txt")
