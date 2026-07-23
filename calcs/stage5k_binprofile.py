"""
STAGE 5K (O10c + O5 binary side): the binary p-profile and the
geometric-mean bootstrap, four seeds.

Functions through the v7 perspective-corrected 2D likelihood (4X/5F/5H
patch set), seeds 31/101/202/303, all on one machinery footing:
  p050  = the occupation law (existing BE table)   nu(1)=1.582
  p0578 = the sec-3 flat screening optimum (new)   nu(1)=1.487
  p065  = the 5G hier screening optimum (existing) nu(1)=1.423
  gm    = geometric-mean bootstrap (new; c1=1/3, c2=1/12, tail p=3/4,
          nu(1)=1.433) -- the O5 construction
Also generates the kappa-variant EFE tables (e_N = 1.0, 1.4) for p065 and
gm, used by the a0-ladder assembly (O10b): kappa = dln(B-1)/dln a0 needs
d/d e_N across tables.

Gates: spherical identity at e_N=0 (<2%) for each new function table;
nu(1) checks; the p050 seeds 31/101 must land within ~0.3 of the stored
4X lam=1.00 values (-56338.377 / -56352.062) -- same-machinery regression.
Writes data/stage5k_summary.txt and the new tables
data/efe_boost_{p0578,gm}_g1p2.npy, data/efe_boost_{p065,gm}_g{1p0,1p4}.npy.
"""
import sys, time
import numpy as np

src_solver = open('calcs/qumond_efe_solver.py', encoding='utf-8-sig').read()
ns = {}
exec(src_solver.split('# --- Gate 1')[0], ns)
solve, r = ns['solve'], ns['r']

def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))
def nu_p(y, p):
    yc = np.clip(np.asarray(y, float), 1e-14, None)
    ex = np.exp(-np.minimum(yc**p, 60.0))
    return (1.0-ex)**(-1.0/(2.0*p))
def nu_p0578(y): return nu_p(y, 0.578)
def nu_p065(y): return nu_p(y, 0.65)
def nu_gm(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    a = y**0.75
    w = np.sqrt(nu_simple(y))
    for _ in range(30):
        u = np.minimum(a*w, 60.0)
        eu = np.exp(u)
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        H = w*w - 1.0 - n
        dH = 2.0*w + a*eu/(em1*em1)
        w = np.maximum(w - H/dH, 1e-8)
    return w*w

y = 1.0/r**2
win = (y > 0.05) & (y < 30)

def make_table(nu, eN, path, gate=False):
    b = solve(nu, eN)
    np.save(path, np.stack([y, b]))
    msg = f"table {path} (eN={eN})"
    if gate:
        b_iso = solve(nu, 0.0)
        t = np.max(np.abs(b_iso[win]/nu(y)[win]-1.0))
        msg += f"  [spherical identity {100*t:.2f}% " \
               f"{'PASS' if t < 0.02 else 'FAIL'}]"
        assert t < 0.02, f"identity gate failed for {path}"
    print(msg, flush=True)
    return path

t0 = time.time()
P0578 = make_table(nu_p0578, 1.2, 'data/efe_boost_p0578_g1p2.npy', gate=True)
PGM = make_table(nu_gm, 1.2, 'data/efe_boost_gm_g1p2.npy', gate=True)
make_table(nu_p065, 1.0, 'data/efe_boost_p065_g1p0.npy')
make_table(nu_p065, 1.4, 'data/efe_boost_p065_g1p4.npy')
make_table(nu_gm, 1.0, 'data/efe_boost_gm_g1p0.npy')
make_table(nu_gm, 1.4, 'data/efe_boost_gm_g1p4.npy')
print(f"nu(1): p0578 {float(nu_p0578(np.array([1.0]))[0]):.3f}, "
      f"gm {float(nu_gm(np.array([1.0]))[0]):.3f}  "
      f"({(time.time()-t0)/60:.1f} min tables)", flush=True)

# ---------------- binary fits ----------------
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
NEW_OUT = "with open('data/stage5k_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1
src = src.replace(OLD_OUT, NEW_OUT)

SEEDS = ['31', '101', '202', '303']
RUNS = [('p050', 'data/efe_boost_be_g1p2.npy'),
        ('p0578', P0578),
        ('p065', 'data/efe_boost_p065_g1p2.npy'),
        ('gm', PGM)]
for law, path in RUNS:
    ns2 = {'__name__': '__main__', 'LAMPATH': path, 'LAMLAW': law}
    sys.argv = ['stage5k', '1p2'] + SEEDS
    print(f"\n===== {law} on binaries =====", flush=True)
    exec(compile(src, f'stage3p_patched_5k_{law}', 'exec'), ns2)

print("\nregression: p050 seeds 31/101 vs 4X lam1.00 "
      "(-56338.377 / -56352.062)")
print("STAGE 5K done -> data/stage5k_summary.txt")
