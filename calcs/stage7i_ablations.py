"""
STAGE 7I: external-review response — the reviewer-specified robustness rows
+ the ceiling-pair release (review solicited by the author 2026-07-24; the
adopted points are logged in NOTES; the function FREEZE that accompanies this
stage is registered in PREDICTIONS.md in the same commit).

Three parts, all pre-registered here BEFORE execution:

C. THE CEILING-PAIR RELEASE (CPU, deterministic). Reproduces the 4J census
   construction verbatim (wide s >= 6 kAU, top-gamma >= 75 deg, band
   [1.414, 1.67)) plus the 4Q perspective-corrected recount, and writes
   data/ceiling_pairs.csv: every wide top-gamma pair with vt_raw in
   [1.20, 2.6) (band + below-edge + above-cliff context) or in the corrected
   band, with source_id1/2, positions, s, distance, masses, vt/gamma under
   BOTH velocity conventions, systemic RV + perspective exposure kappa,
   S/N, per-pair vt noise, RUWE both, R_chance, |dplx|/sig — enough to
   re-count under any convention (sqrt2-vs-1.414 edge, S/N cut, gamma cut).
   GATE GC: 4J-construction raw band count == 11 under BOTH S/N conventions
   (4J: identical); corrected-recount in [9, 12] (4Q found the census
   +-1-sensitive at the edges: their sqrt2-edge construction gave 10 -> 9).

S-PRE. THE STRICT-MULTIPLICITY MASK (CPU). RUWE < 1.2 both components AND
   no detectable overluminosity: both components have valid BP-RP in
   (-0.5, 6) and delta >= -0.4 mag off the 41-bin running-median color
   ridge (the 3J criterion, ridge built on the BASELINE 14,071-pair
   sample's 2N stars). Saves data/stage7i_strictmask.npy (full-catalog
   boolean). GATE GS: the recomputed overluminous star fraction must land
   in [0.11, 0.14] (3J measured 0.123) — else ABORT part S.

W. w_rad FROZEN TO THE EXTERNAL VALUE (GPU). The one reviewer-contested
   internal nuisance (eccentricity mixture fitted jointly with the boost)
   replaced by the external prior: WR_GRID = [0.21], the midpoint of the
   Hwang+22-implied near-parabolic fraction on our range (4G: 20-22% at
   e > 0.9). Everything else = the 4R corrected-velocity baseline
   (stage3p_v7budget.py + 6G NEW_DATA patch, g1p2 tables, standard laws).
   EXPECTATION (pre-registered): near-null — the free fit already chose
   0.20 in 12/12 fits; freezing at the external value within grid
   resolution should move nothing.

S. THE STRICT-MULTIPLICITY FIT (GPU). Same baseline machinery, data pairs
   additionally required to pass the strict mask. The model keeps its
   companion grid (fcomp in {0, 0.1}) — under the cut the fit is free to
   prefer 0. EXPECTATION: alpha stable (2D found the median boost
   RUWE-stable; 3K exonerated companions); Newton margin shrinks roughly
   with N — the honest comparison is per-pair rate, reported.

PAIRING + BARS (pre-registered). Baseline = data/stage4r_summary.txt
per-seed rows (corrected velocities, g1p2, seeds 31/101/202/303/404/505).
Run seeds [31, 101] first; per law D = mean_s[ alpha_hat_variant(s) -
alpha_hat_baseline(s) ].
  BW (w_rad):  CLOSED if max-law |D| <= 0.10 and max-law |mean dNewton|
               <= 15.  EXTEND to [202,303,404,505] if 0.10 < |D| <= 0.20.
               MATERIAL if |D| > 0.20 (extend too, then disclose).
  BS (strict): STABLE if max-law |D| <= 0.15 and alpha interior in every
               run and per-law mean Newton margin >= +30 absolute.
               EXTEND if 0.15 < |D| <= 0.25 or any non-interior.
               MATERIAL if |D| > 0.25 (extend too, then disclose).
Final verdict at whatever seed count the rule reaches, same thresholds on
the full-mean. MATERIAL = the fence was load-bearing -> disclosure +
investigation, not suppression.

Outputs: data/stage7i_precheck.txt, data/ceiling_pairs.csv,
data/stage7i_w.txt, data/stage7i_s.txt, data/stage7i_verdict.txt.
"""
import math
import os
import re
import sys
import time

import numpy as np
from astropy.io import fits

t00 = time.time()
LOG = []
def P(s):
    print(s, flush=True)
    LOG.append(s)

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
P(f"STAGE 7I: {int(ok.sum())} baseline pairs")

# ---------- Part C: the ceiling-pair release (4J construction verbatim) ----
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
vt = (4.74047/plx*vmag)/vc
sig_vt = sigv/np.sqrt(2)/vc

# corrected velocities (6G NEW_DATA formulas, gate-validated in 4R)
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
nrv = h1_.astype(int) + h2_.astype(int)
th = sep/(2.06265e8/plx)
pmcor = rvs*th*plx/4.74047
vxc = vx + pmcor*sux
vyc = vy + pmcor*suy
vmagc = np.hypot(vxc, vyc)
vtc = (4.74047/plx*vmagc)/vc
cosgc = np.abs(sx*vxc+sy*vyc)/np.maximum(smag*vmagc, 1e-12)
gamc = np.degrees(np.arccos(np.clip(cosgc, 0, 1)))
kap = rvs*th/vc     # perspective exposure (4Q convention, km/s / km/s)

WIDE = ok & (s_kau >= 6)
def census(v, g, sncut):
    return WIDE & (sn > sncut) & (g >= 75) & (v >= 1.414) & (v < 1.67)
n_nocut = int(census(vt, gam, 0.0).sum())
n_sn3 = int(census(vt, gam, 3.0).sum())
n_corr = int(census(vtc, gamc, 0.0).sum())
gc_pass = (n_nocut == 11) and (n_sn3 == 11) and (9 <= n_corr <= 12)
P(f"census [1.414,1.67), gamma>=75, s>=6: raw no-cut {n_nocut}, raw S/N>3 "
  f"{n_sn3} (4J: 11/11), corrected {n_corr} (gate [9,12])")
P(f"GATE GC: {'PASS' if gc_pass else 'FAIL'}")

sel = WIDE & ((gam >= 75) & (vt >= 1.20) & (vt < 2.6)
              | (gamc >= 75) & (vtc >= 1.414) & (vtc < 2.6))
idx = np.where(sel)[0]
idx = idx[np.argsort(-vt[idx])]
def band_of(v):
    return ('below' if v < 1.414 else 'band' if v < 1.67
            else 'above' if v < 2.2 else 'far')
rows = ["source_id1,source_id2,ra1_deg,dec1_deg,s_kAU,dist_pc,Mtot_Msun,"
        "vc_kms,vt_raw,gamma_raw_deg,vt_corr,gamma_corr_deg,n_rv,"
        "rv_sys_kms,kappa_persp,sn_direction,sig_vt,ruwe1,ruwe2,"
        "Rchance_x1e3,dplx_over_sig,band_raw,band_corr,"
        "census_raw_sn3,census_corr"]
for i in idx:
    rows.append(
        f"{int(d['source_id1'][i])},{int(d['source_id2'][i])},"
        f"{d['ra1'][i]:.6f},{d['dec1'][i]:.6f},{s_kau[i]:.2f},"
        f"{1000/plx[i]:.1f},{Mtot[i]:.3f},{vc[i]:.4f},{vt[i]:.4f},"
        f"{gam[i]:.2f},{vtc[i]:.4f},{gamc[i]:.2f},{nrv[i]},{rvs[i]:.2f},"
        f"{kap[i]:.4f},{sn[i]:.2f},{sig_vt[i]:.4f},{ruwe1[i]:.3f},"
        f"{ruwe2[i]:.3f},{1e3*Rch[i]:.4f},"
        f"{abs(plx1[i]-plx2[i])/math.hypot(eplx1[i], eplx2[i]):.2f},"
        f"{band_of(vt[i])},{band_of(vtc[i])},"
        f"{bool(census(vt, gam, 3.0)[i])},{bool(census(vtc, gamc, 0.0)[i])}")
with open('data/ceiling_pairs.csv', 'w') as f:
    f.write("\n".join(rows)+"\n")
P(f"data/ceiling_pairs.csv: {len(rows)-1} rows (band + [1.20,1.414) context "
  f"+ above-cliff, both conventions)")

# ---------- Part S-pre: the strict-multiplicity mask ----------------------
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
delta_b = mgb - np.interp(colb, cc[vr], ridge[vr])
f_over = float(np.mean(delta_b[goodb] < -0.4))
gs_pass = 0.11 <= f_over <= 0.14
P(f"overluminous star fraction (baseline sample, 3J criterion): "
  f"{f_over:.3f}  GATE GS (3J 0.123, band [0.11,0.14]): "
  f"{'PASS' if gs_pass else 'FAIL -> ABORT part S'}")
d1 = MG1 - np.interp(c1, cc[vr], ridge[vr])
d2 = MG2 - np.interp(c2, cc[vr], ridge[vr])
okc = np.isfinite(c1) & np.isfinite(c2) & (c1 > -0.5) & (c1 < 6.0) \
    & (c2 > -0.5) & (c2 < 6.0)
mask = (ruwe1 < 1.2) & (ruwe2 < 1.2) & okc & (d1 >= -0.4) & (d2 >= -0.4)
np.save('data/stage7i_strictmask.npy', mask)
n_strict = int((ok & mask).sum())
P(f"strict cut: RUWE<1.2 both {int((ok & (ruwe1 < 1.2) & (ruwe2 < 1.2)).sum())}"
  f"/{int(ok.sum())}; + no-overluminous + valid color -> {n_strict}"
  f"/{int(ok.sum())} pairs ({n_strict/max(int(ok.sum()), 1):.3f})")

with open('data/stage7i_precheck.txt', 'w') as f:
    f.write("\n".join(LOG)+"\n")
if not gc_pass:
    print("GC FAIL — census not reproduced; investigate before GPU parts.")
    sys.exit(1)

# ---------- GPU variants: patched stage3p_v7budget --------------------------
src0 = open('calcs/stage3p_v7budget.py', encoding='utf-8-sig').read()

OLD_DATA = r"""vt_d = (4.74047/plx[ok]*np.hypot(d['pmra1'][ok]-d['pmra2'][ok],
        d['pmdec1'][ok]-d['pmdec2'][ok]))/(0.9417*np.sqrt(Mtot_d/s_d))
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
vx_ = d['pmra2']-d['pmra1']; vy_ = d['pmdec2']-d['pmdec1']"""
NEW_DATA = r"""dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
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
vt_d = (4.74047/plx[ok]*np.hypot(vx_[ok], vy_[ok]))/(0.9417*np.sqrt(Mtot_d/s_d))"""

OK_OLD = r"""ok = (Rch<0.01)&(plx1>5)&(plx2>5)&(plx1/np.maximum(eplx1,1e-6)>20) \
   &(plx2/np.maximum(eplx2,1e-6)>20)&(np.abs(plx1-plx2)<3*np.hypot(eplx1,eplx2)) \
   &(sep>200)&(sep<50000)&(MG1>2.6)&(MG1<14.2)&(MG2>2.6)&(MG2<14.2)&(sigv<0.03)"""
OK_NEW = OK_OLD + """
ok0_ = int(ok.sum())
ok = ok & np.load('data/stage7i_strictmask.npy')
print(f"7I strict cut active: {int(ok.sum())}/{ok0_} pairs", flush=True)"""

WR_OLD = "WR_GRID = np.array([0.10, 0.20, 0.30])"
WR_NEW = "WR_GRID = np.array([0.21])"
OUT_OLD = "with open('data/stage3u_summary.txt', 'a') as f:"

def make_src(tag, strict, wfroze):
    s = src0
    assert s.count(OLD_DATA) == 1
    s = s.replace(OLD_DATA, NEW_DATA)
    if strict:
        assert s.count(OK_OLD) == 1
        s = s.replace(OK_OLD, OK_NEW)
    if wfroze:
        assert s.count(WR_OLD) == 1
        s = s.replace(WR_OLD, WR_NEW)
    assert s.count(OUT_OLD) == 1
    s = s.replace(OUT_OLD, f"with open('data/stage7i_{tag}.txt', 'a') as f:")
    return s

ROW = re.compile(r"seed (\d+) (simple|BE): a_hat=([0-9.]+) \(grid ([0-9.]+), "
                 r"interior=(\w+)\), dlnL\(Newton\)=([+-][0-9.]+), "
                 r"wr=([0-9.]+)")
def parse(path):
    out = {}
    if os.path.exists(path):
        for m in ROW.finditer(open(path).read()):
            s_, law, ah, _, inter, dn, wr = m.groups()
            out[(law, int(s_))] = dict(a=float(ah), interior=(inter == 'True'),
                                       dn=float(dn), wr=float(wr))
    return out

BASE = parse('data/stage4r_summary.txt')
for s_ in (31, 101, 202, 303, 404, 505):
    assert ('simple', s_) in BASE and ('BE', s_) in BASE, s_

def run_variant(tag, strict, wfroze, seeds):
    have = parse(f'data/stage7i_{tag}.txt')
    todo = [s_ for s_ in seeds
            if ('simple', s_) not in have or ('BE', s_) not in have]
    if todo:
        print(f"===== 7I variant {tag}: seeds {todo} =====", flush=True)
        src_v = make_src(tag, strict, wfroze)
        sys.argv = ['stage7i', '1p2'] + [str(s_) for s_ in todo]
        exec(compile(src_v, f'stage3p_patched_7i_{tag}', 'exec'),
             {'__name__': '__main__'})
    return parse(f'data/stage7i_{tag}.txt')

def deltas(rec, seeds):
    out = {}
    for law in ('simple', 'BE'):
        da = [rec[(law, s_)]['a'] - BASE[(law, s_)]['a'] for s_ in seeds]
        dn = [rec[(law, s_)]['dn'] for s_ in seeds]
        dnb = [BASE[(law, s_)]['dn'] for s_ in seeds]
        ints = [rec[(law, s_)]['interior'] for s_ in seeds]
        out[law] = dict(D=float(np.mean(da)),
                        se=(float(np.std(da, ddof=1)/np.sqrt(len(da)))
                            if len(da) > 1 else 0.0),
                        dn=float(np.mean(dn)), dnb=float(np.mean(dnb)),
                        n_int=sum(ints), n=len(seeds))
    return out

SEED0 = [31, 101]
SEEDX = [202, 303, 404, 505]
V = ["STAGE 7I verdicts (bars pre-registered in the module docstring)", ""]

# --- variant W ---
rec = run_variant('w', strict=False, wfroze=True, seeds=SEED0)
dl = deltas(rec, SEED0)
mx = max(abs(dl['simple']['D']), abs(dl['BE']['D']))
mxn = max(abs(dl['simple']['dn']-dl['simple']['dnb']),
          abs(dl['BE']['dn']-dl['BE']['dnb']))
extw = mx > 0.10
if extw:
    rec = run_variant('w', strict=False, wfroze=True, seeds=SEED0+SEEDX)
    dl = deltas(rec, SEED0+SEEDX)
    mx = max(abs(dl['simple']['D']), abs(dl['BE']['D']))
    mxn = max(abs(dl['simple']['dn']-dl['simple']['dnb']),
              abs(dl['BE']['dn']-dl['BE']['dnb']))
vw = ('CLOSED' if (mx <= 0.10 and mxn <= 15) else
      'SHIFT-DISCLOSED' if mx <= 0.20 else 'MATERIAL')
V.append(f"W (w_rad frozen 0.21, external Hwang value; {dl['simple']['n']} "
         f"seeds{', EXTENDED' if extw else ''}):")
for law in ('simple', 'BE'):
    x = dl[law]
    V.append(f"  {law:>6}: d_alpha = {x['D']:+.3f} +- {x['se']:.3f}; "
             f"Newton {x['dn']:+.1f} (baseline {x['dnb']:+.1f}); "
             f"interior {x['n_int']}/{x['n']}")
V.append(f"  BW verdict: {vw} (max|D|={mx:.3f}, max|dNewton|={mxn:.1f})")
V.append("")

# --- variant S ---
if gs_pass:
    rec = run_variant('s', strict=True, wfroze=False, seeds=SEED0)
    dl = deltas(rec, SEED0)
    mx = max(abs(dl['simple']['D']), abs(dl['BE']['D']))
    all_int = dl['simple']['n_int'] == dl['simple']['n'] and \
        dl['BE']['n_int'] == dl['BE']['n']
    exts = (mx > 0.15) or (not all_int)
    if exts:
        rec = run_variant('s', strict=True, wfroze=False, seeds=SEED0+SEEDX)
        dl = deltas(rec, SEED0+SEEDX)
        mx = max(abs(dl['simple']['D']), abs(dl['BE']['D']))
        all_int = dl['simple']['n_int'] == dl['simple']['n'] and \
            dl['BE']['n_int'] == dl['BE']['n']
    minn = min(dl['simple']['dn'], dl['BE']['dn'])
    vs = ('STABLE' if (mx <= 0.15 and all_int and minn >= 30) else
          'SHIFT-DISCLOSED' if mx <= 0.25 else 'MATERIAL')
    fkeep = n_strict/14071
    V.append(f"S (strict multiplicity: RUWE<1.2 both + no overluminous; "
             f"{n_strict}/14071 pairs = {fkeep:.3f}; {dl['simple']['n']} "
             f"seeds{', EXTENDED' if exts else ''}):")
    for law in ('simple', 'BE'):
        x = dl[law]
        V.append(f"  {law:>6}: d_alpha = {x['D']:+.3f} +- {x['se']:.3f}; "
                 f"Newton {x['dn']:+.1f} (baseline {x['dnb']:+.1f}; "
                 f"per-pair ratio {x['dn']/max(x['dnb']*fkeep, 1e-9):.2f} "
                 f"of N-scaled expectation); interior {x['n_int']}/{x['n']}")
    V.append(f"  BS verdict: {vs} (max|D|={mx:.3f}, min Newton={minn:+.1f})")
else:
    V.append("S: ABORTED (GS gate failed)")
V.append("")
V.append(f"census: raw 11/11 reproduced, corrected {n_corr}; "
         f"data/ceiling_pairs.csv released ({len(rows)-1} rows)")
V.append(f"({(time.time()-t00)/60:.1f} min total)")

out = "\n".join(V)
print("\n"+out)
with open('data/stage7i_verdict.txt', 'w') as f:
    f.write(out+"\n")
print("\nSTAGE 7I done -> data/stage7i_verdict.txt")
