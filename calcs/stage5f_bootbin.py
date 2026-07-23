"""
STAGE 5F (O7): the quantum-bootstrap function on the BINARIES.

5C/5D found the bath matrix REORDERS under hierarchical M/L: the implicit
quantum self-consistent cell  nu = 1 + n_BE(nu*y)  (c1=1/4, c2=7/96, e^-y
screening; frequency set by the TOTAL acceleration -- the Unruh-ratio
reading of paper sec. 2.4) beats BE by 76 lnL on SPARC+lensing (86% of
galaxy bootstraps). The 4X binary dial read c1 = 0.37-0.50 through the
LAMBDA-MIXTURE family, whose low-c1 members carry standard-mu's shape --
boot is a different function at the same c1. Direct test: solve the QUMOND
EFE for nu_boot at the physical field (g_N,ext = 1.2 a0), run the v7
perspective-corrected 2D binary likelihood (stage3p machinery, seeds
31+101), and put boot's absolute lnL on the 4X grid footing.

Gates: T2-analog spherical identity at e_N=0 (<2%); transition sanity
(boot nu(1) = 1.35 < BE 1.58).
Comparators (data/stage4x_summary.txt): seed 31 lam1.00 -56338.377,
lam0.75 -56340.783; seed 101 lam1.00 -56352.062, lam0.75 -56344.988.
Writes data/stage5f_summary.txt (+ boot table data/efe_boost_boot_g1p2.npy).
"""
import sys
import numpy as np

# ---------- part 1: boot table ----------
src_solver = open('calcs/qumond_efe_solver.py', encoding='utf-8-sig').read()
ns = {}
exec(src_solver.split('# --- Gate 1')[0], ns)
solve, r = ns['solve'], ns['r']
nu_be_ = ns['nu_be']

def nu_boot(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    u = 0.5*(y + np.sqrt(y*y + 4.0*y))
    for _ in range(14):
        uc = np.minimum(u, 45.0)
        eu = np.exp(uc)
        em1 = eu - 1.0
        F = u - y - y/em1
        dF = 1.0 + y*eu/(em1*em1)
        u = np.maximum(u - F/dF, 1e-13)
    nu = u/y
    big = y > 45.0
    if np.any(big):
        yb = np.minimum(y, 700.0)
        nu = np.where(big, 1.0 + 1.0/np.expm1(yb), nu)
    return nu

y = 1.0/r**2
EN = 1.2
b = solve(nu_boot, EN)
PBOOT = 'data/efe_boost_boot_g1p2.npy'
np.save(PBOOT, np.stack([y, b]))
print(f"table boot -> {PBOOT}")

win = (y > 0.05) & (y < 30)
b_iso = solve(nu_boot, 0.0)
t2 = np.max(np.abs(b_iso[win]/nu_boot(y)[win]-1.0))
print(f"T2 boot spherical identity: max dev {100*t2:.2f}% -> "
      f"{'PASS' if t2 < 0.02 else 'FAIL'}")
nb1 = float(nu_boot(np.array([1.0]))[0])
print(f"sanity: boot nu(1) = {nb1:.3f} (BE {float(nu_be_(np.array([1.0]))[0]):.3f})")
assert t2 < 0.02 and abs(nb1-1.35) < 0.02, "boot table gates failed"

# ---------- part 2: binary fits (4X patch set, boot table) ----------
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
NEW_OUT = "with open('data/stage5f_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1
src = src.replace(OLD_OUT, NEW_OUT)

ns2 = {'__name__': '__main__', 'LAMPATH': PBOOT, 'LAMLAW': 'boot'}
sys.argv = ['stage5f', '1p2', '31', '101']
print("\n===== boot on binaries =====")
exec(compile(src, 'stage3p_patched_5f_boot', 'exec'), ns2)

COMP = {'31': {'lam1.00': -56338.377, 'lam0.75': -56340.783},
        '101': {'lam1.00': -56352.062, 'lam0.75': -56344.988}}
print("\n4X comparators (same machinery): "
      + " | ".join(f"seed {s}: BE {COMP[s]['lam1.00']:.3f}, "
                   f"lam0.75 {COMP[s]['lam0.75']:.3f}" for s in COMP))
print("STAGE 5F done -> data/stage5f_summary.txt")
