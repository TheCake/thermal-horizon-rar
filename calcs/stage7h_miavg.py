"""
STAGE 7H: the trajectory-averaged prescription -- was 7G's -8.4 the
proxy's Jensen gap?

7G measured the mi_t proxy (per-orbit boost AT the Kepler time-averaged
acceleration, B(<g>)) trailing the field formulation by -8.43 with
alpha-hat inflated to 1.55/2.00 (field fit: 1.06). SIGNED HYPOTHESIS,
stated before the run: the proxy under-boosts by Jensen's gap --
eccentric orbits spend their time near apoapsis where g is LOW and the
boost is HIGH, so the honest adiabatic trajectory functional is the
time-average OF the boost,

    <B> = (1/2pi) INT B(g_N(r(E))) (1 - e cos E) dE,   r = a(1 - e cos E),

which is strictly larger than B(<g>) for a decreasing convex boost curve
sampled over the population's wide-e orbits (w_rad = 0.20 near-parabolic:
apo/peri g excursions of 100x+). Prediction: mi_avg RECOVERS boost,
alpha-hat falls toward the field value, and the lnL gap closes partly or
fully. The honest alternative: the gap persists -> the trajectory class
pays a REAL binary cost (the door narrows, quantified).

Machinery: 7G's patch pipeline verbatim (corrected velocities, gated AMB
table, single-law mode, exact-Newton engine with G_eff) with the engine
block computing <B> by eccentric-anomaly quadrature (NE = 48 nodes,
vectorized) instead of the single-point interp. Comparators: same-seed
6G MG-AMB rows; 7G mi_t rows quoted alongside.

Gates:
  G0 alpha = 0 bit-identity vs the 6G Newton rows (< 0.5), as 7G.
  G1 circular identity (standalone, pre-fit): at e = 1e-6, <B>(a) must
     equal B(g(a)) to 1e-8 across the population's a-range.
  G2 quadrature convergence: NE = 48 vs 96 max |d<B>| < 1e-6 at the
     (a, e) extremes (e up to 0.999).
PRE-REGISTERED BARS (extension rule fixed NOW -- the 7G lesson):
  seeds 31/101 first. TIE-RESTORED if mean D(mi_avg-AMB minus MG-AMB)
  in [-5, +5] with alpha-hat interior both -> 7G's gap was proxy
  discretization (Jensen). If the 2-seed result is AMBIGUOUS (mean in
  [-10, -5), or interior 1/2, or |seed spread| > 8), EXTEND to seeds
  202/303/404/505 and judge on all six: TIE-RESTORED as above on the
  6-seed mean +- SE; CLASS-COST if 6-seed mean <= -8 sign-consistent
  (>= 5/6) -> the trajectory reading pays a real binary price.
  Also recorded: alpha-hat movement (prediction: toward ~1.06).
Writes data/stage7h_miavg.txt (+ data/stage7h_summary.txt runs).
"""
import math, os, re, sys, time, gc
import numpy as np

L = []
def say(s=''):
    L.append(s); print(s, flush=True)
def save():
    with open('data/stage7h_miavg.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

say("STAGE 7H: the trajectory-averaged prescription (Jensen-gap test)")
say("=" * 72)

# ---------------- G1/G2 standalone gates on the <B> quadrature --------------
# reproduce the fit's table pipeline: load the gated AMB table, alpha-scale
tab = np.load('data/efe_boost_amb_g1p2.npy')
y_t, b_t = tab[0][::-1], tab[1][::-1]
LNY_U = np.log(y_t)
A0_CAN = 1.0                      # canonical units as in the fit (y already
GM = 1.0                          # in a0 units on the table axis)

def bavg_quad(lny_a, e, tab_a, NE=48):
    """time-average of B over the orbit; lny_a = ln(g_N(a)/a0)."""
    E = (np.arange(NE) + 0.5)*(2*np.pi/NE)
    ome = 1.0 - np.outer(e, np.cos(E))          # (N, NE): r/a
    lny = lny_a[:, None] - 2.0*np.log(np.maximum(ome, 1e-6))
    B = np.interp(lny, LNY_U, tab_a, right=1.0)
    w = ome/np.sum(ome, axis=1, keepdims=True)  # dt weight prop (1 - e cosE)
    return np.sum(B*w, axis=1)

tab_test = b_t                                   # alpha = 1 table values
lny_grid = np.linspace(np.log(0.05), np.log(30.0), 40)
e0 = np.full(len(lny_grid), 1e-6)
bav0 = bavg_quad(lny_grid, e0, tab_test)
bpt0 = np.interp(lny_grid, LNY_U, tab_test, right=1.0)
g1 = float(np.max(np.abs(bav0 - bpt0)))
ee = np.array([0.0, 0.5, 0.9, 0.99, 0.999])
lyx = np.full(len(ee), math.log(0.05))
d48_96 = float(np.max(np.abs(bavg_quad(lyx, ee, tab_test, 48)
                             - bavg_quad(lyx, ee, tab_test, 96))))
lyx2 = np.full(len(ee), math.log(30.0))
d48_96 = max(d48_96, float(np.max(np.abs(
    bavg_quad(lyx2, ee, tab_test, 48) - bavg_quad(lyx2, ee, tab_test, 96)))))
say(f"G1 circular identity (e=1e-6): max |<B> - B| = {g1:.2e} -> "
    f"{'PASS' if g1 < 1e-8 else 'FAIL'}")
say(f"G2 quadrature NE 48 vs 96 (e to 0.999): max d = {d48_96:.2e} -> "
    f"{'PASS' if d48_96 < 1e-6 else 'FAIL'}")
jens = bavg_quad(np.full(1, math.log(1.0)), np.array([0.9]), tab_test)[0] \
    - np.interp(math.log(1.0/ (1-0.81)**0.5), LNY_U, tab_test)
say(f"  (Jensen direction at y=1, e=0.9: <B> - B(<g>) = {jens:+.4f} -- "
    f"positive = boost recovered)")
say('')
save()
assert g1 < 1e-8 and d48_96 < 1e-6

# ---------------- the patched run -------------------------------------------
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
NEW_OUT = "with open('data/stage7h_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1
src = src.replace(OLD_OUT, NEW_OUT)

OLD_ENG = """                    e_s = e_of(p, eta, wr)
                    vp = vp_c(p, e_s, tab_a) if al > 0 else None
                    mode = 5 if al > 0 else 1
                    kw = dict(a0=A0_CAN, tab=tab_a, lny0=LNY0, dlny=DLNY,
                              vp=vp) if al > 0 else {}
                    o = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'],
                            p['uph'], 8, 2500, mode, **kw)"""
NEW_ENG = """                    e_s = e_of(p, eta, wr)
                    lny_a = np.log(np.maximum(GM*p['M_s']/p['a_s']**2
                                              / A0_CAN, 1e-12))
                    NEq = 48
                    Eq = (np.arange(NEq) + 0.5)*(2*np.pi/NEq)
                    ome_ = 1.0 - np.outer(e_s, np.cos(Eq))
                    lny_o = lny_a[:, None] - 2.0*np.log(
                        np.maximum(ome_, 1e-6))
                    B_o = np.interp(lny_o, LNY_U, tab_a, right=1.0)
                    w_o = ome_/np.sum(ome_, axis=1, keepdims=True)
                    Bp = np.sum(B_o*w_o, axis=1)
                    o = run(p['a_s'], e_s, p['psi0'], p['f_ip'],
                            p['M_s']*Bp, p['uph'], 8, 2500, 1)"""
assert src.count(OLD_ENG) == 1
src = src.replace(OLD_ENG, NEW_ENG)

def have(k):
    return os.path.exists('data/stage7h_summary.txt') and \
        open('data/stage7h_summary.txt').read().count(' ambav: ') >= k

SEEDS1 = ['31', '101']
SEEDS2 = ['202', '303', '404', '505']
if not have(len(SEEDS1)):
    ns2 = {'__name__': '__main__',
           'LAMPATH': 'data/efe_boost_amb_g1p2.npy', 'LAMLAW': 'ambav'}
    sys.argv = ['stage7h', '1p2'] + SEEDS1
    print("\n===== mi_avg-AMB (trajectory-averaged) on binaries =====",
          flush=True)
    exec(compile(src, 'stage3p_patched_7h_ambav', 'exec'), ns2)
    del ns2; gc.collect()

def parse(path):
    txt = open(path).read()
    out = {}
    for m in re.finditer(
            r'seed (\d+) (\w+): a_hat=([0-9.]+) \(grid ([0-9.]+), '
            r'interior=(\w+)\), dlnL\(Newton\)=\+([0-9.]+), wr=([0-9.]+),'
            r'.*?\n\s*seed \1: best lnL = (-?[0-9.]+)', txt):
        s, law, ah, agrid, inter, dn, wr, lnl = m.groups()
        out[(law, int(s))] = dict(a=float(ah), interior=(inter == 'True'),
                                  lnl=float(lnl), dn=float(dn))
    return out

REC = {}
for p in ('data/stage6g_summary.txt', 'data/stage7g_summary.txt',
          'data/stage7h_summary.txt'):
    if os.path.exists(p): REC.update(parse(p))

def rowset(seeds):
    ds, ahs, ints, spread, g0ok = [], [], 0, 0.0, True
    for s in seeds:
        if ('ambav', s) not in REC or ('amb', s) not in REC: continue
        mi, mg = REC[('ambav', s)], REC[('amb', s)]
        g0ok &= abs((mi['lnl']-mi['dn']) - (mg['lnl']-mg['dn'])) < 0.5
        ds.append(mi['lnl'] - mg['lnl']); ahs.append(mi['a'])
        ints += mi['interior']
    return np.array(ds), np.array(ahs), ints, g0ok

ds2, ah2, int2, g0a = rowset([31, 101])
ambig = (len(ds2) == 2) and ((-10.0 <= ds2.mean() < -5.0) or int2 == 1
                             or abs(ds2[0]-ds2[1]) > 8.0)
if ambig and not have(6):
    say("2-seed result AMBIGUOUS per pre-registered rule -> extending to 6")
    save()
    ns2 = {'__name__': '__main__',
           'LAMPATH': 'data/efe_boost_amb_g1p2.npy', 'LAMLAW': 'ambav'}
    sys.argv = ['stage7h', '1p2'] + SEEDS2
    exec(compile(src, 'stage3p_patched_7h_ambav2', 'exec'), ns2)
    del ns2; gc.collect()
    REC.update(parse('data/stage7h_summary.txt'))

SE_ALL = [31, 101] + ([202, 303, 404, 505] if have(6) else [])
ds, ahs, ints, g0ok = rowset(SE_ALL)
say('')
say("mi_avg-AMB vs same-seed MG-AMB (corrected velocities):")
for i, s in enumerate(SE_ALL):
    if i >= len(ds): break
    mi = REC[('ambav', s)]
    mt = REC.get(('ambmi', s))
    say(f"  seed {s}: D = {ds[i]:+.2f}; a_hat = {ahs[i]:.2f} "
        f"(interior={mi['interior']}); dlnL(Newton) = +{mi['dn']:.1f}"
        + (f"   [7G mi_t D was {mt['lnl']-REC[('amb', s)]['lnl']:+.2f}, "
           f"a_hat {mt['a']:.2f}]" if mt else ""))
se = ds.std(ddof=1)/math.sqrt(len(ds)) if len(ds) > 1 else float('nan')
say(f"  mean D = {ds.mean():+.2f} +- {se:.2f} SE ({len(ds)} seeds); "
    f"interior {ints}/{len(ds)}; a_hat {ahs.mean():.3f}"
    + (f"+-{ahs.std(ddof=1)/math.sqrt(len(ahs)):.3f}" if len(ahs) > 1 else ""))
say(f"G0 bit-identity: {'PASS' if g0ok else 'FAIL'}")

tie = (len(ds) >= 2) and (-5.0 <= ds.mean() <= 5.0) and ints == len(ds)
cost = (len(ds) == 6) and ds.mean() <= -8.0 and (ds < 0).sum() >= 5
if tie:
    verd = ("TIE-RESTORED: the 7G gap was the proxy's Jensen "
            "discretization -- the trajectory formulation fits the "
            "binaries at field grade (and Saturn at 10^451 margin)")
elif cost:
    verd = "CLASS-COST: the trajectory reading pays a real binary price"
else:
    verd = "PARTIAL/AMBIGUOUS -- see rows"
say('')
say(f"VERDICT (pre-registered): {verd}")
say("carried: adiabatic-average = the next MI bracket, still not a full")
say("nonlocal theory; MG-formulation tension stands in the world table.")
save()
print("\nSTAGE 7H done -> data/stage7h_miavg.txt")
