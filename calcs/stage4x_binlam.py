"""
STAGE 4X (O2): the binary-side c1 dial. The Stage 4S lambda-family
(nu_lam = (1-lam)*standard + lam*BE, c1 = lam/2) run through the QUMOND EFE
solver at the physical field, then the v7 binary likelihood (perspective-
corrected data, per Stage 4R) profiled over lambda with alpha and all
nuisances refit per node. An independent, x~1-regime reading of the same
dial the galaxies measured at c1_hat = 0.45.

Part 1 gates: T1 the lam=1 table reproduces the stored BE g1p2 table to <1%
in the fit window; T2 the lam=0.5 table at e_N=0 reproduces the spherical
identity nu_lam(y) to <2%.
Part 2: seeds 31+101 per lambda; absolute best lnL per (lambda, seed)
written to data/stage4x_summary.txt (same machinery across lambda => direct
lnL(lambda) profile).
"""
import sys
import numpy as np

# ---------- part 1: tables ----------
src_solver = open('calcs/qumond_efe_solver.py', encoding='utf-8-sig').read()
ns = {}
exec(src_solver.split('# --- Gate 1')[0], ns)
solve, r = ns['solve'], ns['r']
nu_be_ = ns['nu_be']
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
y = 1.0/r**2
EN = 1.2
LAMS = [0.0, 0.25, 0.5, 0.75, 1.0, 1.25]
paths = {}
for lam in LAMS:
    def nul(yy, l=lam):
        return (1.0-l)*nu_standard(yy) + l*nu_be_(yy)
    b = solve(nul, EN)
    p = f'data/efe_boost_lam{int(round(lam*100)):03d}_g1p2.npy'
    np.save(p, np.stack([y, b]))
    paths[lam] = p
    print(f"table lam={lam:.2f} -> {p}")

win = (y > 0.05) & (y < 30)
ref = np.load('data/efe_boost_be_g1p2.npy')
b1 = np.load(paths[1.0])[1]
t1 = np.max(np.abs(b1[win]/ref[1][win]-1.0))
print(f"T1 lam=1 vs stored BE table: max dev {100*t1:.2f}% -> "
      f"{'PASS' if t1 < 0.01 else 'FAIL'}")
def nul05(yy):
    return 0.5*nu_standard(yy) + 0.5*nu_be_(yy)
b05_iso = solve(nul05, 0.0)
t2 = np.max(np.abs(b05_iso[win]/nul05(y)[win]-1.0))
print(f"T2 lam=0.5 spherical identity: max dev {100*t2:.2f}% -> "
      f"{'PASS' if t2 < 0.02 else 'FAIL'}")
assert t1 < 0.01 and t2 < 0.02, "table gates failed"

# ---------- part 2: binary fits ----------
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
assert src.count(OLD_DATA) == 1
src = src.replace(OLD_DATA, NEW_DATA)

OLD_TABS = """if GTAG == '1p9':
    PS, PB = 'data/efe_boost_simple.npy', 'data/efe_boost_be.npy'
else:
    PS = f'data/efe_boost_simple_g{GTAG}.npy'
    PB = f'data/efe_boost_be_g{GTAG}.npy'
LNY_U, TAB_S = load_tab(PS)
_,     TAB_B = load_tab(PB)"""
NEW_TABS = """LNY_U, TAB_S = load_tab(LAMPATH)
TAB_B = TAB_S"""
assert src.count(OLD_TABS) == 1
src = src.replace(OLD_TABS, NEW_TABS)

OLD_LOOP = 'for law, TAB in (("simple", TAB_S), ("BE", TAB_B)):'
NEW_LOOP = 'for law, TAB in ((LAMLAW, TAB_S),):'
assert src.count(OLD_LOOP) == 1
src = src.replace(OLD_LOOP, NEW_LOOP)

OLD_SUM = """    P(f"  seed {seed}: BE-minus-simple best lnL = "
      f"{best_lnl['BE']-best_lnl['simple']:+.1f}  "
      f"({(time.time()-t0)/60:.1f} min)")"""
NEW_SUM = """    P(f"  seed {seed}: best lnL = {best_lnl[LAMLAW]:+.3f}  "
      f"({(time.time()-t0)/60:.1f} min)")"""
assert src.count(OLD_SUM) == 1
src = src.replace(OLD_SUM, NEW_SUM)

OLD_OUT = "with open('data/stage3u_summary.txt', 'a') as f:"
NEW_OUT = "with open('data/stage4x_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1
src = src.replace(OLD_OUT, NEW_OUT)

SEEDS = ['31', '101']
for lam in LAMS:
    ns2 = {'__name__': '__main__',
           'LAMPATH': paths[lam],
           'LAMLAW': f'lam{lam:.2f}'}
    sys.argv = ['stage4x', '1p2'] + SEEDS
    print(f"\n===== lambda = {lam:.2f} =====")
    exec(compile(src, f'stage3p_patched_4x_lam{lam:.2f}', 'exec'), ns2)
print("\nSTAGE 4X done -> data/stage4x_summary.txt")
