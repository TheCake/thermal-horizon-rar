"""
STAGE 4R: the corrected-velocity v7 budget (TODO #2i - the Stage 4Q spherical-
projection correction executed at the fit level, not footnoted).

Patch-runs stage3p_v7budget.py with the perspective correction applied in the
data build: Dv_corr = Dv + RV_sys * theta * s_hat (catalog RVs, inverse-variance
systemic mean; conventions identical to calcs/stage4q_perspective.py, whose G2
round-trip validates the sign). Both vt_d and gam_d inherit the correction.
Physical-field tables (g=1p2), both laws, the six budget seeds.
Output: data/stage4r_summary.txt (the stored stage3u baseline is untouched).

Expected from 4Q: alpha-hat shift within the 1.6%-of-ratio exposure bound;
Newton dlnL ~unchanged; the per-seed BE-minus-simple lean re-measured on
corrected data (the bath-matrix question rides on it).
"""
import sys

src = open('calcs/stage3p_v7budget.py', encoding='utf-8-sig').read()

OLD_DATA = """vt_d = (4.74047/plx[ok]*np.hypot(d['pmra1'][ok]-d['pmra2'][ok],
        d['pmdec1'][ok]-d['pmdec2'][ok]))/(0.9417*np.sqrt(Mtot_d/s_d))
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
vx_ = d['pmra2']-d['pmra1']; vy_ = d['pmdec2']-d['pmdec1']"""

NEW_DATA = """dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
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
rvs_ = (np.where(h1_, r1_, 0.0)*w1_ + np.where(h2_, r2_, 0.0)*w2_) \\
       / np.maximum(w1_+w2_, 1e-12)
th_ = sep/(2.06265e8/plx)
pmcor_ = rvs_*th_*plx/4.74047
vx_ = d['pmra2']-d['pmra1'] + pmcor_*sux_
vy_ = d['pmdec2']-d['pmdec1'] + pmcor_*suy_
vt_d = (4.74047/plx[ok]*np.hypot(vx_[ok], vy_[ok]))/(0.9417*np.sqrt(Mtot_d/s_d))"""

assert src.count(OLD_DATA) == 1, "stage3p data block not found verbatim"
src = src.replace(OLD_DATA, NEW_DATA)

OLD_OUT = "with open('data/stage3u_summary.txt', 'a') as f:"
NEW_OUT = "with open('data/stage4r_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1, "summary path not found"
src = src.replace(OLD_OUT, NEW_OUT)

sys.argv = ['stage4r'] + (sys.argv[1:] if len(sys.argv) > 1
                          else ['1p2', '31', '101', '202', '303', '404', '505'])
ns = {'__name__': '__main__'}
exec(compile(src, 'stage3p_patched_4r', 'exec'), ns)
