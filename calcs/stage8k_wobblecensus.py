"""STAGE 8K — THE WOBBLE CENSUS (pre-reg 370eed4, committed BEFORE
this run; the reviewer's object-level tie-breaker, fully local).

Does the boost-carrying population show the astrometric-activity
signature the collapse world requires (fcomp 0.35-0.50 per component
of active-but-saturated companions)?  HOT = ruwe > 1.25 OR
aen_sig > 3 (pre-registered constants).  S1 wide rate (validity:
no catalog RUWE ceiling), S2 vt-activity correlation (asymmetric
interpretation pre-stated), S3 the nine census pairs; S4/S5
descriptive.  Verdict OBJECT-LEVEL-ABSENT / WOBBLE-RICH / MIXED;
this stage carries its own pre-registered credence map.
Output: data/stage8k_wobblecensus.txt
"""
import csv
import numpy as np
from astropy.io import fits
from scipy.stats import spearmanr

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

P("8K THE WOBBLE CENSUS (pre-reg 370eed4; bars + credence map "
  "locked before any run)")
P("")

# ---------- data + masks: verbatim 4J construction ----------
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
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx = dra*np.cos(dec_m); sy = d['dec2']-d['dec1']
vx = d['pmra2']-d['pmra1']; vy = d['pmdec2']-d['pmdec1']
e_vx = np.sqrt(d['pmra_error1']**2+d['pmra_error2']**2)
e_vy = np.sqrt(d['pmdec_error1']**2+d['pmdec_error2']**2)
vmag = np.hypot(vx, vy)
sn = vmag/np.sqrt(0.5*(e_vx**2+e_vy**2))
smag = np.hypot(sx, sy)
cosg = np.abs(sx*vx+sy*vy)/np.maximum(smag*vmag, 1e-12)
gam = np.degrees(np.arccos(np.clip(cosg, 0, 1)))
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot = np.interp(MG1, MG_T, MS_T)+np.interp(MG2, MG_T, MS_T)
s_kau = sep/1e3
vc = 0.9417*np.sqrt(Mtot/np.maximum(s_kau, 1e-9))
dv = 4.74047/plx*vmag
vt = dv/vc

ru1, ru2 = d['ruwe1'], d['ruwe2']
ae1, ae2 = d['astrometric_excess_noise1'], d['astrometric_excess_noise2']
as1, as2 = (d['astrometric_excess_noise_sig1'],
            d['astrometric_excess_noise_sig2'])
dru1, dru2 = d['dr2_ruwe1'], d['dr2_ruwe2']
rv1, rv2 = d['dr2_radial_velocity1'], d['dr2_radial_velocity2']
re1, re2 = d['dr2_radial_velocity_error1'], d['dr2_radial_velocity_error2']

# ---------- gates ----------
SBINS = [(0.2, 2), (2, 6), (6, 20), (20, 50)]
nb = [int(((s_kau >= a) & (s_kau < b) & ok).sum()) for a, b in SBINS]
g0 = int(ok.sum()) == 14071 and nb == [9950, 2684, 1223, 214]
P(f"G8K-0 loader identity: N = {int(ok.sum())} (ref 14071), per-bin "
  f"{nb} (ref [9950, 2684, 1223, 214]) -> {'PASS' if g0 else 'FAIL'}")
assert g0

cens = [r for r in csv.DictReader(open('data/ceiling_pairs.csv'))
        if r['census_corr'] == 'True']
cidx = []
oki = np.where(ok)[0]
for r in cens:
    m = (np.abs(s_kau[oki]-float(r['s_kAU'])) < 0.01) \
      & (np.abs(Mtot[oki]-float(r['Mtot_Msun'])) < 0.01) \
      & (np.abs(vc[oki]-float(r['vc_kms'])) < 0.001)
    j = oki[m]
    assert len(j) == 1, f"census match ambiguity: {len(j)} rows"
    cidx.append(int(j[0]))
P(f"G8K-0b census match: {len(cidx)}/9 matched uniquely -> "
  f"{'PASS' if len(cidx) == 9 else 'FAIL'}")
assert len(cidx) == 9

nan_fr = float(np.mean(~np.isfinite(np.concatenate(
    [ru1[ok], ru2[ok], ae1[ok], ae2[ok]]))))
P(f"G8K-1 flag completeness: NaN fraction = {nan_fr:.4f} -> "
  f"{'PASS' if nan_fr < 0.05 else 'DISCLOSED (>5%)'}")

ruall = np.concatenate([ru1[ok], ru2[ok]])
ruall = ruall[np.isfinite(ruall)]
q = np.percentile(ruall, [50, 90, 99])
P(f"RUWE distribution (all ok components): P50 = {q[0]:.3f}, "
  f"P90 = {q[1]:.3f}, P99 = {q[2]:.3f}, max = {ruall.max():.2f}")
s1_valid = q[2] > 1.6
P(f"S1 VALIDITY (pre-stated: P99 > 1.6, no hard ceiling): "
  f"{'VALID' if s1_valid else 'DISCLOSED-INVALID (ceiling detected)'}")
P("")

# ---------- HOT flags ----------
# AMENDMENT 1 (pre-quote, logged in NOTES; run-1 preserved as
# _run1.txt): the aen_sig>3 arm fired on 0.96-0.99 of G<12
# components (a brightness flag, caught by G8K-2's own mag table).
# HOT = ruwe > 1.25 ONLY (Belokurov-class standard; RUWE is
# brightness-normalized).  aen_sig retained descriptively below.
hot1 = (ru1 > 1.25)
hot2 = (ru2 > 1.25)
pair_hot = hot1 | hot2
WIDE = ok & (s_kau >= 6)
NARR = ok & (s_kau < 2)

def wilson(k, n):
    if n == 0:
        return (0.0, 0.0)
    p = k/n; z = 1.96
    den = 1+z*z/n
    c = (p+z*z/(2*n))/den
    h = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return (c-h, c+h)

# G8K-2 mag-dependence
P("G8K-2 f_hot vs G-mag (1-mag bins, all ok components):")
gall = np.concatenate([G1m[ok], G2m[ok]])
hall = np.concatenate([hot1[ok], hot2[ok]])
fr_mag = []
for glo in range(8, 17):
    m = (gall >= glo) & (gall < glo+1) & np.isfinite(gall)
    if m.sum() > 100:
        fr = float(hall[m].mean())
        fr_mag.append(fr)
        P(f"  G {glo}-{glo+1}: f_hot = {fr:.3f} (n = {int(m.sum())})")
magrun = max(fr_mag)/max(min(fr_mag), 1e-9) if fr_mag else 1.0
P(f"  mag run (max/min) = {magrun:.2f} "
  + ("-> STRONG (mag-binned percentile variant primary for S2)"
     if magrun > 2 else "-> mild"))
P("")

# ---------- S1: the rate ----------
for nm, msk in (('WIDE (s>=6)', WIDE), ('NARROW control (<2)', NARR)):
    kc = int(hot1[msk].sum()+hot2[msk].sum())
    nc = int(2*msk.sum())
    lo, hi = wilson(kc, nc)
    P(f"S1 [{nm}]: f_hot per component = {kc/nc:.3f} "
      f"[{lo:.3f}, {hi:.3f}] (k = {kc}, n = {nc}); pair-hot = "
      f"{float(pair_hot[msk].mean()):.3f}")
kW = int(hot1[WIDE].sum()+hot2[WIDE].sum()); nW = int(2*WIDE.sum())
f1 = kW/nW
s1cat = ('RATE-SHORTFALL' if f1 < 0.20 else
         'RATE-MET' if f1 >= 0.30 else 'GRAY')
P(f"S1 verdict (bars < 0.20 / >= 0.30): {s1cat}"
  + ("" if s1_valid else " [DISCLOSED-INVALID - fallback to S2+S3]"))
P("")

# ---------- S2: the correlation (WIDE) ----------
vtW = vt[WIDE]; phW = pair_hot[WIDE].astype(float)
d_obs = float(vtW[phW > 0.5].mean() - vtW[phW < 0.5].mean())
rng = np.random.default_rng(53)
null = np.empty(10000)
for i in range(10000):
    sh = rng.permutation(phW)
    null[i] = vtW[sh > 0.5].mean() - vtW[sh < 0.5].mean()
p_perm = float(np.mean(np.abs(null) >= abs(d_obs)))
# PRIMARY variant (G8K-2 trigger fired): per-component RUWE
# percentile within its 1-mag bin; pair score = max of the two
allg = np.concatenate([G1m[ok], G2m[ok]])
allr = np.concatenate([ru1[ok], ru2[ok]])
def magpct(gv, rv):
    out = np.zeros(len(gv))
    for glo in range(6, 19):
        m = (allg >= glo) & (allg < glo+1) & np.isfinite(allr)
        mm = (gv >= glo) & (gv < glo+1)
        if m.sum() > 50 and mm.sum() > 0:
            ref = np.sort(allr[m])
            out[mm] = np.searchsorted(ref, rv[mm])/len(ref)
    return out
pc1 = magpct(G1m[WIDE], ru1[WIDE])
pc2 = magpct(G2m[WIDE], ru2[WIDE])
score = np.maximum(pc1, pc2)
rho, p_sp = spearmanr(score, vtW)
P(f"S2 [WIDE]: mean vt(hot) - mean vt(cold) = {d_obs:+.4f} "
  f"(perm p = {p_perm:.4f}); PRIMARY mag-binned: spearman("
  f"max-RUWE-pctile, vt) rho = {rho:+.4f} (p = {p_sp:.2e})")
P("   (asymmetric reading pre-stated: positive = necessary for "
  "collapse but error-degenerate; null = collapse-killing for the "
  "vt-excess-as-wobble reading)")
P("")

# ---------- S3: the nine ----------
alla = np.concatenate([as1[ok], as2[ok]])
def aenpct(g, a):
    m = (allg >= g-0.5) & (allg < g+0.5) & np.isfinite(alla)
    if m.sum() < 50:
        return float('nan')
    ref = np.sort(alla[m])
    return float(np.searchsorted(ref, a)/len(ref))
nhot = 0
P("S3 the nine census pairs (individually; aen_sig percentile = "
  "vs mag-matched ok components, DESCRIPTIVE):")
for j, r in zip(cidx, cens):
    h = bool(pair_hot[j])
    nhot += h
    ap1 = aenpct(G1m[j], as1[j]); ap2 = aenpct(G2m[j], as2[j])
    P(f"  s = {s_kau[j]:6.2f} kAU, vt = {vt[j]:.3f}: ruwe = "
      f"({ru1[j]:.3f}, {ru2[j]:.3f}), aen_sig = ({as1[j]:.1f}, "
      f"{as2[j]:.1f}) [mag-pctile ({ap1:.2f}, {ap2:.2f})], "
      f"dr2_ruwe = ({dru1[j]:.3f}, {dru2[j]:.3f}) "
      f"-> {'HOT' if h else 'clean'}")
s3cat = ('CENSUS-CLEAN' if nhot <= 3 else
         'CENSUS-WOBBLE-SUSPECT' if nhot >= 6 else 'MIXED')
P(f"S3 verdict: {nhot}/9 hot (bars <= 3 / >= 6) -> {s3cat}")
P("")

# ---------- S4/S5 descriptive ----------
rv_good = (np.isfinite(rv1) & np.isfinite(rv2)
           & np.isfinite(re1) & np.isfinite(re2)
           & (re1 > 0) & (re1 < 1e3) & (re2 > 0) & (re2 < 1e3))
both_rv = WIDE & rv_good
if both_rv.sum() > 10:
    drv = np.abs(rv1[both_rv]-rv2[both_rv])
    P(f"S4 [desc] WIDE pairs with both DR2 RVs (sentinels filtered): "
      f"n = {int(both_rv.sum())}; median |dRV| = "
      f"{float(np.median(drv)):.2f} km/s; median rv_err hot/cold = "
      f"{float(np.median(np.concatenate([re1[both_rv & pair_hot], re2[both_rv & pair_hot]]))):.2f}"
      f"/{float(np.median(np.concatenate([re1[both_rv & ~pair_hot], re2[both_rv & ~pair_hot]]))):.2f} km/s")
dr_ok = WIDE & np.isfinite(dru1) & np.isfinite(dru2)
ddr = np.concatenate([(ru1-dru1)[dr_ok & pair_hot],
                      (ru2-dru2)[dr_ok & pair_hot]])
ddc = np.concatenate([(ru1-dru1)[dr_ok & ~pair_hot],
                      (ru2-dru2)[dr_ok & ~pair_hot]])
P(f"S5 [desc] DR2->EDR3 RUWE delta (WIDE): median hot = "
  f"{float(np.median(ddr)):+.3f}, cold = {float(np.median(ddc)):+.3f}")
P("")

# ---------- verdict + map ----------
if s1_valid:
    absent = (s1cat == 'RATE-SHORTFALL') and (s3cat == 'CENSUS-CLEAN')
else:
    absent = (s3cat == 'CENSUS-CLEAN') and (p_perm > 0.05)
rich = (s1cat == 'RATE-MET') and (nhot >= 6)
if absent:
    v = ("OBJECT-LEVEL-ABSENT: the collapse world's required "
         "active-companion population is not in the catalog at the "
         "required rate, and the census pairs are astrometrically "
         "clean - the 4J/8F-b defense UPGRADES (the band pairs are "
         "neither error tails nor wobblers); per the pre-registered "
         "map: anomaly-real ~50% -> ~55%; the 24-seed budget is "
         "DE-PRIORITIZED; reviewer to be notified")
elif rich:
    v = ("WOBBLE-RICH: the required population exists where the "
         "collapse world needs it; per the map: anomaly-real ~50% -> "
         "~40%; the census's gravity-evidence reading is SUSPENDED "
         "pending the NSS leg; the budget proceeds under the "
         "registered map")
else:
    v = ("MIXED-REPORTED: neither grammar fires clean; per the map: "
         "~50% HELD; the budget proceeds; decomposition above is "
         "the product")
P(f"==> 8K VERDICT (locked bars): {v}")

with open('data/stage8k_wobblecensus.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8k_wobblecensus.txt")
