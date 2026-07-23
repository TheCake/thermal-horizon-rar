"""
STAGE 5R (O13a): the binary beta-profile -- dedicated beta = 1/4, 3/4 EFE
tables + six seeds, completing the 6x5 lnL(beta) matrix.

The mixing family (5P): nu = 1 + n_BE(y^((1+beta)/2) nu^beta), exact
c1 = 1/(2(1+beta)), p_tail = (1+beta)/2, c1*p_tail = 1/4. Galaxies read
beta-hat = 0.45-0.64 (Delta -43..-86 in the -2lnL objective vs beta=0).
Binaries so far hold beta ~ 0 on members only: gm (beta=1/2) +8.5+-2.1
behind per seed (6 seeds, 5K/5O), boot (beta=1) +17/+22 behind (2 seeds,
5F). This stage measures HOW SHARP the binary beta=0 preference is:
  - new tables nu_beta(0.25), nu_beta(0.75) at e_N = 1.2 a0
  - fits on all six seeds (31/101/202/303/404/505)
  - boot completed on the four seeds 5F did not run
  - assembled per-seed profiles + upper bounds + the two-system table.

Gates: G1 end-to-end member regression (nu_beta(b=0) table through the
solver vs the stored BE table); G2 spherical identity (<2%) per new table;
G3 boost-vs-beta monotonicity report; the 4X/5K patch asserts (exact-count
string replacement) guarantee machinery identity with 5F/5K/5O.
Writes data/stage5r_summary.txt (fits) + data/stage5r_profile.txt (matrix).
"""
import os, re, sys, time
import numpy as np

src_solver = open('calcs/qumond_efe_solver.py', encoding='utf-8-sig').read()
ns = {}
exec(src_solver.split('# --- Gate 1')[0], ns)
solve, r = ns['solve'], ns['r']

def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))

def nu_beta(y, beta):
    """nu = 1 + n_BE(y^((1+beta)/2) * nu^beta), Newton on nu (F monotone)."""
    y = np.clip(np.asarray(y, float), 1e-14, None)
    A = y**(0.5*(1.0+beta))
    nu = nu_simple(y)
    for _ in range(40):
        u = np.minimum(A*nu**beta, 60.0)
        eu = np.exp(u)
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        F = nu - 1.0 - n
        dF = 1.0 + (eu/(em1*em1))*u*beta/nu
        nu = np.maximum(nu - F/dF, 1.0 + 1e-15)
    return nu

y = 1.0/r**2
win = (y > 0.05) & (y < 30)

# G1: end-to-end member regression, beta=0 table vs stored BE table
t0 = time.time()
b0 = solve(lambda yy: nu_beta(yy, 0.0), 1.2)
yb, bb = np.load('data/efe_boost_be_g1p2.npy')
g1 = float(np.max(np.abs(b0[win]/bb[win] - 1.0)))
print(f"G1 member regression (beta=0 table vs stored BE): max rel "
      f"{g1:.2e} -> {'PASS' if g1 < 5e-4 else 'FAIL'}", flush=True)
assert g1 < 5e-4, "G1 member regression failed"

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

P025 = make_table(lambda yy: nu_beta(yy, 0.25), 1.2,
                  'data/efe_boost_b025_g1p2.npy', gate=True)
P075 = make_table(lambda yy: nu_beta(yy, 0.75), 1.2,
                  'data/efe_boost_b075_g1p2.npy', gate=True)
print(f"({(time.time()-t0)/60:.1f} min tables)", flush=True)

# G3 monotonicity report across the family at eN=1.2
tabs = {0.0: 'data/efe_boost_be_g1p2.npy', 0.25: P025,
        0.5: 'data/efe_boost_gm_g1p2.npy', 0.75: P075,
        1.0: 'data/efe_boost_boot_g1p2.npy'}
print("G3 boost B(y) vs beta at eN=1.2:", flush=True)
for yq in (0.1, 0.3, 1.0, 3.0):
    i = int(np.argmin(np.abs(y - yq)))
    row = [f"{np.load(p)[1][i]:.4f}" for p in tabs.values()]
    print(f"  y={yq:4.1f}: " + " ".join(row), flush=True)

# ---------------- binary fits (4X/5K patch set, verbatim) ----------------
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
NEW_OUT = "with open('data/stage5r_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1
src = src.replace(OLD_OUT, NEW_OUT)

if os.path.exists('data/stage5r_summary.txt'):
    os.remove('data/stage5r_summary.txt')

SEEDS6 = ['31', '101', '202', '303', '404', '505']
RUNS = [('b025', P025, SEEDS6),
        ('b075', P075, SEEDS6),
        ('boot', 'data/efe_boost_boot_g1p2.npy', ['202', '303', '404', '505'])]
for law, path, seeds in RUNS:
    ns2 = {'__name__': '__main__', 'LAMPATH': path, 'LAMLAW': law}
    sys.argv = ['stage5r', '1p2'] + seeds
    print(f"\n===== {law} on binaries (seeds {seeds}) =====", flush=True)
    exec(compile(src, f'stage3p_patched_5r_{law}', 'exec'), ns2)

# ---------------- profile assembly ----------------
def parse(path):
    txt = open(path).read()
    out = {}
    for m in re.finditer(
            r'seed (\d+) (\w+): a_hat=([0-9.]+) \(grid [0-9.]+, '
            r'interior=(\w+)\), dlnL\(Newton\)=\+([0-9.]+), wr=([0-9.]+),'
            r'.*?\n\s*seed \1: best lnL = (-?[0-9.]+)', txt):
        s, law, ah, inter, dn, wr, lnl = m.groups()
        out[(law, int(s))] = dict(a=float(ah), interior=(inter == 'True'),
                                  dN=float(dn), wr=float(wr),
                                  lnl=float(lnl))
    return out

REC = {}
for p in ('data/stage5k_summary.txt', 'data/stage5o_summary.txt',
          'data/stage5f_summary.txt', 'data/stage5r_summary.txt'):
    REC.update(parse(p))

BMAP = [(0.0, 'p050'), (0.25, 'b025'), (0.5, 'gm'), (0.75, 'b075'),
        (1.0, 'boot')]
SEEDS = [31, 101, 202, 303, 404, 505]
L = ["STAGE 5R: the binary beta-profile (6 seeds x 5 beta; members from "
     "5K/5O/5F, b025/b075/boot-completion new)", ""]
L.append("per-seed lnL(beta) - lnL(0)  [negative = worse than beta=0]:")
D = np.zeros((len(SEEDS), len(BMAP)))
A = np.zeros((len(SEEDS), len(BMAP)))
for i, s in enumerate(SEEDS):
    row = []
    for j, (b, law) in enumerate(BMAP):
        rec = REC[(law, s)]
        D[i, j] = rec['lnl'] - REC[('p050', s)]['lnl']
        A[i, j] = rec['a']
        row.append(f"{D[i, j]:+7.2f}")
    L.append(f"  seed {s:>3}: " + " ".join(row))
mn, sd = D.mean(0), D.std(0, ddof=1)
se = sd/np.sqrt(len(SEEDS))
L.append("  " + "-"*58)
L.append("  mean    : " + " ".join(f"{v:+7.2f}" for v in mn))
L.append("  seed SE : " + " ".join(f"{v:7.2f}" for v in se))
L.append("  sign<0  : " + " ".join(f"{int((D[:, j] < 0).sum())}/6"
                                   for j in range(len(BMAP))))
L.append("  a_hat   : " + " ".join(f"{A[:, j].mean():7.3f}"
                                   for j in range(len(BMAP))))
L.append("  interior: " + " ".join(
    f"{sum(REC[(law, s)]['interior'] for s in SEEDS)}/6"
    for _, law in BMAP))
L.append("")

# upper bounds from the mean profile (linear crossing of -0.5 / -2.0)
bg = np.array([b for b, _ in BMAP])
mprof = mn.copy()   # 0 at beta=0 by construction
def crossing(thr):
    for j in range(1, len(bg)):
        if mprof[j] < -thr:
            b1, b2 = bg[j-1], bg[j]
            v1, v2 = mprof[j-1], mprof[j]
            return b1 + (b2-b1)*(-thr - v1)/(v2 - v1)
    return None
c1s, c2s = crossing(0.5), crossing(2.0)
L.append(f"binary upper bounds (mean profile, DlnL crossings): "
         f"beta < {c1s:.3f} (1 sigma, DlnL=0.5), "
         f"beta < {c2s:.3f} (2 sigma, DlnL=2.0)")
L.append("NOTE: realization scatter dominates -- the seed SE column is the "
         "honest per-beta uncertainty; the crossing uses the seed-mean.")
L.append("")

# the two-system table (galaxy grids from 5P, -2lnL objective, lower=better)
GAL_PLAIN = {0.0: 0.0, 0.25: -62.63, 0.5: -84.79, 0.75: -85.59, 1.0: -75.62}
GAL_VERT = {0.0: 0.0, 0.25: -38.11, 0.5: -42.69, 0.75: -30.91, 1.0: -11.08}
L.append("two-system profile [Delta(-2lnL) vs beta=0; negative = preferred]:")
L.append("   beta   gal plain   gal vert   binaries(-2*DlnL, +-2SE)")
for j, b in enumerate(bg):
    L.append(f"  {b:5.2f}   {GAL_PLAIN[b]:+9.2f}  {GAL_VERT[b]:+9.2f}   "
             f"{-2*mn[j]:+8.2f} +- {2*2*se[j]:.2f}")
L.append("")
L.append("joint (vertical-hardened galaxies + binaries, additive -2lnL):")
for j, b in enumerate(bg):
    L.append(f"  beta={b:4.2f}: {GAL_VERT[b] + (-2*mn[j]):+8.2f}")
L.append("CAVEAT: the galaxy yardstick under galaxy-bootstrap (5N) is "
         "+-~50 on this contrast (population variance) -- the joint row is "
         "indicative, not a likelihood-ratio verdict.")

out = "\n".join(L)
print("\n" + out)
with open('data/stage5r_profile.txt', 'w') as f:
    f.write(out + "\n")
print("\nSTAGE 5R done -> data/stage5r_profile.txt")
