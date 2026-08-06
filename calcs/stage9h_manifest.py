"""STAGE 9H — THE DATA-PROVENANCE MANIFEST.  Pre-registered.

SHA256 every load-bearing input + recompute the historical
invariants LIVE (verbatim loader expressions) vs stage-of-record
printed values.  Verdict: DATA-VERIFIED iff every invariant
matches; DRIFT-DETECTED else.  NO credence movement.
Output: data/stage9h_manifest.txt + data/MANIFEST.sha256
"""
import csv as _csv
import glob, hashlib, math, os
import numpy as np
from astropy.io import fits

L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

P("9H THE DATA-PROVENANCE MANIFEST (pre-reg committed BEFORE the "
  "run; NO credence movement)")
P("")

# ---------------- hashes ----------------
def sha(path, bufsz=1 << 22):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            b = f.read(bufsz)
            if not b: break
            h.update(b)
    return h.hexdigest()

FILES = [
    'data/edr3_binaries.fits.gz',
    'data/sparc/SPARC_Lelli2016c.mrt',
    'data/chae2021_table3.csv',
    'data/efe_boost_simple_g1p2.npy',
    'data/efe_boost_be_g1p2.npy',
    'data/stage7j_cube_full_lker_31_simple.npy',
    'data/stage7j_cube_full_lker_31_BE.npy',
    'data/stage7j_cube_full_lker_101_simple.npy',
    'data/stage7j_cube_full_lker_101_BE.npy',
    'data/stage7jz_prior.npz',
    'data/stage9f_tables_31_simple.npz',
    'data/stage9f_tables_31_BE.npz',
    'data/stage9f_tables_101_simple.npz',
    'data/stage9f_tables_101_BE.npz',
]
man = []
ok_files = True
for p in FILES:
    if not os.path.exists(p):
        P(f"  MISSING: {p}"); ok_files = False; continue
    man.append((p, os.path.getsize(p), sha(p)))
rotm = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                        recursive=True))
hr = hashlib.sha256()
for p in rotm:
    hr.update(sha(p).encode())
man.append((f'data/sparc/rotmod/** ({len(rotm)} files)',
            sum(os.path.getsize(p) for p in rotm), hr.hexdigest()))
P(f"hashed {len(man)} entries ({len(rotm)} rotmod files "
  f"aggregate); manifest -> data/MANIFEST.sha256")
with open('data/MANIFEST.sha256', 'w') as f:
    for p, sz, hx in man:
        f.write(f"{hx}  {sz:>12}  {p}\n")
P("")

# ---------------- invariants: binary side (verbatim 9F loader) ----
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
checks = []
def chk(name, got, want, tol=0.0):
    okc = (abs(got - want) <= tol) if tol else (got == want)
    checks.append(okc)
    P(f"  [{name:24}] live = {got}  record = {want}  -> "
      f"{'PASS' if okc else 'FAIL'}")
    return okc

P("invariants (live recompute vs stage-of-record):")
chk('mask N (8Z)', int(ok.sum()), 14071)
_cn = list(d.columns.names)
def _pick(*names):
    for n in names:
        if n in _cn: return np.asarray(d[n], dtype=np.float64)
    raise KeyError(names)
ruwe = np.maximum(_pick('ruwe1', 'RUWE1'), _pick('ruwe2', 'RUWE2'))[ok]
qs_ = np.percentile(ruwe, [25, 50, 75])
chk('ruwe q25 (8Z/9A)', round(float(qs_[0]), 3), 1.051)
chk('ruwe q50 (8Z/9A)', round(float(qs_[1]), 3), 1.118)
chk('ruwe q75 (8Z/9A)', round(float(qs_[2]), 3), 1.231)
plxsn = np.minimum(plx1/np.maximum(eplx1, 1e-6),
                   plx2/np.maximum(eplx2, 1e-6))[ok]
qp_ = np.percentile(plxsn, [25, 50, 75])
chk('plxsn q25 (8Z)', round(float(qp_[0]), 3), 317.311)
chk('plxsn q50 (8Z)', round(float(qp_[1]), 3), 407.883)
chk('plxsn q75 (8Z)', round(float(qp_[2]), 3), 570.585)

# 9E statistic verbatim
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
WIDE = (s_d >= 6) & (s_d < 50)
NARR = (s_d >= 0.2) & (s_d < 2)
def B_of(m):
    return float(np.median(vt_d[WIDE & m])/np.median(vt_d[NARR & m]))
ALLm = np.ones(len(s_d), dtype=bool)
Q1m = ruwe <= qs_[0]
chk('B(all) (9E)', round(B_of(ALLm), 4), 1.0779)
chk('B(Q1) (9E)', round(B_of(Q1m), 4), 1.1216)

# ---------------- invariants: galaxy side (verbatim 8V loader) ----
KPC = 3.24078e-14
UD, UB = 0.5, 0.7
meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        name = t[0]
        inc, q = float(t[5]), int(t[17])
        meta[name] = (inc, q)
    except ValueError:
        continue
g_gas, g_dsk, g_bul, gal_id = [], [], [], []
gname = {}
kept = 0
for gi, path in enumerate(rotm):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q = meta.get(name, (0, 3))
    if inc < 30 or q > 2: continue
    kept += 1
    gname[gi] = name
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC
        gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        g_gas.append(gg); g_dsk.append(gd); g_bul.append(gb)
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gal_id = map(np.array,
                                  (g_gas, g_dsk, g_bul, gal_id))
allg = np.unique(gal_id)
# AMENDMENT (pre-quote, logged): the program's 149 = galaxies with
# surviving DATA POINTS (= len(allg) = the fit population), not the
# meta-cut file count (153; 4 galaxies lose every point to the
# per-line eV/Vo and g>0 cuts and never enter any fit).  First
# manifest firing compared the wrong live quantity to the right
# record - the gate caught the manifest's own wiring.
P(f"  (info: rotmod files passing meta cuts = {kept}; galaxies "
  f"with surviving data points = {len(allg)})")
chk('fit galaxies (8V allg)', int(len(allg)), 149)
NGAL = gal_id.max()+1
gdfrac = np.zeros(NGAL)
for g_ in allg:
    m = gal_id == g_
    gdfrac[g_] = float(np.mean(g_gas[m] > g_dsk[m] + g_bul[m]))
gd_set = [g_ for g_ in allg if gdfrac[g_] >= 0.5]
dd_set = [g_ for g_ in allg if gdfrac[g_] < 0.5]
chk('GD count (8S-c)', len(gd_set), 38)
chk('DD count (8S-c)', len(dd_set), 111)
chae = {}
with open('data/chae2021_table3.csv') as f:
    for row in _csv.DictReader(l for l in f if not l.startswith('#')):
        chae[row['galaxy']] = 1
chk('Chae rows (9G)', len(chae), 109)
chk('Chae GD matched (9G)',
    sum(1 for g_ in gd_set if gname.get(g_) in chae), 21)
chk('Chae DD matched (9G)',
    sum(1 for g_ in dd_set if gname.get(g_) in chae), 70)

P("")
if ok_files and all(checks):
    P(f"==> 9H VERDICT: DATA-VERIFIED - {len(checks)} invariants "
      f"match their stage-of-record values; {len(man)} manifest "
      f"entries written.  Any future input drift is one re-run "
      f"away from detection.")
else:
    P("==> 9H VERDICT: DRIFT-DETECTED - STOP; do not run further "
      "stages until resolved.")
P("    NO credence movement (pre-stated).")

with open('data/stage9h_manifest.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage9h_manifest.txt")
