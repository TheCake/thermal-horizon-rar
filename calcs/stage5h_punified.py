"""
STAGE 5H (O9): the unification check -- nu_p(0.65) on the binaries.

5G decomposed the hierarchical bath flip: 75% of boot's gain over BE is the
Newtonian-tail sharpness dial (nu_p family, hier optimum p=0.65, gain -56
of boot's -76), not the deep 1/4 coefficient. The binaries (5F) rejected
boot's WEAK TRANSITION (nu(1)=1.35, alpha edge-ride, +17..+24 lnL) while
holding the 1/2 branch. nu_p(0.65) keeps a 1/2-branch-grade transition
(nu(1)=1.423) with the sharper tail the hierarchical galaxies demand.
If the binaries accept it (alpha interior, lnL on the BE/lam-grid level),
a SINGLE measured function -- the screening family at p ~ 0.6-0.65 --
fits the hierarchical galaxy RAR (+56 over the occupation law), the wide
binaries, and the Cassini sharpness direction simultaneously.

Table: QUMOND EFE solve for nu_p(0.65) at g_N,ext = 1.2 a0, spherical-
identity gate. Fit: v7 perspective-corrected 2D likelihood (4X/5F patch
set), seeds 31+101.
Comparators: seed 31 BE -56338.377 (alpha-hat 1.10 interior), boot
-56360.641 (edge); seed 101 BE -56352.062 (1.13 interior), lam0.75
-56344.988, boot -56369.132 (edge).
Writes data/stage5h_summary.txt (+ data/efe_boost_p065_g1p2.npy).
"""
import sys
import numpy as np

src_solver = open('calcs/qumond_efe_solver.py', encoding='utf-8-sig').read()
ns = {}
exec(src_solver.split('# --- Gate 1')[0], ns)
solve, r = ns['solve'], ns['r']

P_STAR = 0.65
def nu_p65(y):
    yc = np.clip(np.asarray(y, float), 1e-14, None)
    ex = np.exp(-np.minimum(yc**P_STAR, 60.0))
    return (1.0-ex)**(-1.0/(2.0*P_STAR))

y = 1.0/r**2
b = solve(nu_p65, 1.2)
PTAB = 'data/efe_boost_p065_g1p2.npy'
np.save(PTAB, np.stack([y, b]))
print(f"table nu_p(0.65) -> {PTAB}")

win = (y > 0.05) & (y < 30)
b_iso = solve(nu_p65, 0.0)
t2 = np.max(np.abs(b_iso[win]/nu_p65(y)[win]-1.0))
n1 = float(nu_p65(np.array([1.0]))[0])
print(f"T2 spherical identity: max dev {100*t2:.2f}% -> "
      f"{'PASS' if t2 < 0.02 else 'FAIL'}; nu(1) = {n1:.3f}")
assert t2 < 0.02 and abs(n1-1.423) < 0.01

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
NEW_OUT = "with open('data/stage5h_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1
src = src.replace(OLD_OUT, NEW_OUT)

ns2 = {'__name__': '__main__', 'LAMPATH': PTAB, 'LAMLAW': 'p065'}
sys.argv = ['stage5h', '1p2', '31', '101']
print("\n===== nu_p(0.65) on binaries =====")
exec(compile(src, 'stage3p_patched_5h_p065', 'exec'), ns2)

print("\ncomparators: seed 31 BE -56338.377 (interior), boot -56360.641 "
      "(edge) | seed 101 BE -56352.062, lam0.75 -56344.988, boot -56369.132")
print("STAGE 5H done -> data/stage5h_summary.txt")
