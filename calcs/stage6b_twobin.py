"""
STAGE 6B (O15b binary leg): the two-leg functions on the wide binaries,
six seeds, the 4X/5K machinery.

5Z's F3 (beta = 1/(2 nu^2), nu(1) = 1.503) and F4 (beta = 1/(2(2 nu-1)^2),
nu(1) = 1.537) raise the transition toward the binaries' held value
(BE 1.582) while keeping the p = 3/4 tail and now both Bernoulli rungs.
The one-leg functions paid -5.8/-5.5 at the transition (5W); gm paid
-8.5. Pre-stated read: mean within -3 of p050 = the transition ACCEPTS
the two-leg function; the refinement then closes 5W's residual.

Gates: G1 spherical identity (<2%) per table; exact-count patch asserts.
Writes data/stage6b_summary.txt + data/stage6b_verdict.txt.
"""
import os, re, sys, time
import numpy as np

src_solver = open('calcs/qumond_efe_solver.py', encoding='utf-8-sig').read()
ns = {}
exec(src_solver.split('# --- Gate 1')[0], ns)
solve, r = ns['solve'], ns['r']

def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))

def make_twoleg(kind):
    def nu_run(y):
        y = np.clip(np.asarray(y, float), 1e-14, None)
        ly = np.log(y)
        nu = nu_simple(y)
        for _ in range(80):
            if kind == 'f3':
                b = 0.5/(nu*nu)
                db = -1.0/(nu**3)
            else:
                b = 0.5/((2.0*nu - 1.0)**2)
                db = -2.0/((2.0*nu - 1.0)**3)
            u = np.exp(np.minimum(0.5*(1.0+b)*ly + b*np.log(nu), 60.0))
            eu = np.exp(np.minimum(u, 60.0))
            em1 = np.maximum(eu - 1.0, 1e-300)
            n = np.where(u < 60.0, 1.0/em1, 0.0)
            F = nu - 1.0 - n
            dudnu = u*(db*(0.5*ly + np.log(nu)) + b/nu)
            dF = 1.0 + (eu/(em1*em1))*dudnu
            step = F/np.where(np.abs(dF) > 1e-12, dF, 1e-12)
            nu = np.maximum(nu - np.clip(step, -0.5*nu, 0.5*nu),
                            1.0 + 1e-15)
        return nu
    return nu_run

nu_f3 = make_twoleg('f3')
nu_f4 = make_twoleg('f4')

y = 1.0/r**2
win = (y > 0.05) & (y < 30)

def make_table(nu, eN, path):
    b = solve(nu, eN)
    np.save(path, np.stack([y, b]))
    b_iso = solve(nu, 0.0)
    t = np.max(np.abs(b_iso[win]/nu(y)[win]-1.0))
    print(f"table {path} (eN={eN})  [spherical identity {100*t:.2f}% "
          f"{'PASS' if t < 0.02 else 'FAIL'}]", flush=True)
    assert t < 0.02
    return path

t0 = time.time()
P3 = make_table(nu_f3, 1.2, 'data/efe_boost_rb3_g1p2.npy')
P4 = make_table(nu_f4, 1.2, 'data/efe_boost_rb4_g1p2.npy')
print(f"({(time.time()-t0)/60:.1f} min tables)", flush=True)

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
NEW_OUT = "with open('data/stage6b_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1
src = src.replace(OLD_OUT, NEW_OUT)

if os.path.exists('data/stage6b_summary.txt'):
    os.remove('data/stage6b_summary.txt')

SEEDS6 = ['31', '101', '202', '303', '404', '505']
for law, path in (('rb3', P3), ('rb4', P4)):
    ns2 = {'__name__': '__main__', 'LAMPATH': path, 'LAMLAW': law}
    sys.argv = ['stage6b', '1p2'] + SEEDS6
    print(f"\n===== {law} on binaries =====", flush=True)
    exec(compile(src, f'stage3p_patched_6b_{law}', 'exec'), ns2)

# ---------------- verdict assembly ----------------
def parse(path):
    txt = open(path).read()
    out = {}
    for m in re.finditer(
            r'seed (\d+) (\w+): a_hat=([0-9.]+) \(grid [0-9.]+, '
            r'interior=(\w+)\), dlnL\(Newton\)=\+([0-9.]+), wr=([0-9.]+),'
            r'.*?\n\s*seed \1: best lnL = (-?[0-9.]+)', txt):
        s, law, ah, inter, dn, wr, lnl = m.groups()
        out[(law, int(s))] = dict(a=float(ah), interior=(inter == 'True'),
                                  lnl=float(lnl))
    return out

REC = {}
for p in ('data/stage5k_summary.txt', 'data/stage5o_summary.txt',
          'data/stage5w_summary.txt', 'data/stage6b_summary.txt'):
    REC.update(parse(p))

SEEDS = [31, 101, 202, 303, 404, 505]
LAWS = ('gm', 'rb1', 'rb2', 'rb3', 'rb4')
L = ["STAGE 6B: two-leg functions on the binaries (vs stored "
     "p050/gm/rb1/rb2)", ""]
L.append("per-seed lnL - lnL(p050):")
L.append("  seed      gm      rb1      rb2      rb3      rb4")
stats = {k: [] for k in LAWS}
for s in SEEDS:
    base = REC[('p050', s)]['lnl']
    row = []
    for law in LAWS:
        d = REC[(law, s)]['lnl'] - base
        stats[law].append(d)
        row.append(f"{d:+7.2f}")
    L.append(f"  {s:>4}: " + "  ".join(row))
L.append("")
for law in LAWS:
    a = np.array(stats[law])
    n_int = sum(REC[(law, s)]['interior'] for s in SEEDS)
    ah = np.array([REC[(law, s)]['a'] for s in SEEDS])
    L.append(f"  {law:>3}: mean {a.mean():+.2f} +- "
             f"{a.std(ddof=1)/np.sqrt(6):.2f} SE  (worse in "
             f"{int((a < 0).sum())}/6; interior {n_int}/6; a_hat "
             f"{ah.mean():.3f} +- {ah.std(ddof=1)/np.sqrt(6):.3f})")
L.append("")
L.append("pre-stated read: mean within -3 of p050 = transition ACCEPTS.")

out = "\n".join(L)
print("\n" + out)
with open('data/stage6b_verdict.txt', 'w') as f:
    f.write(out + "\n")
print("\nSTAGE 6B done -> data/stage6b_verdict.txt")
