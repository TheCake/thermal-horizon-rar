"""
ROUND-50 GB (post-report verification; the standing verify-reviewer-math
rule). Re-computes every load-bearing Round-50 number in the main
session's own code before adoption. Output: data/round50_gb.txt.

GB-1  P(A1 fires) arithmetic + max-order-statistic null (no fits)
GB-2  the full LOO list re-run (78 fits): sigma_jack, RMS, N_eff,
      max-vs-simple-null percentile
GB-3  the A4 point-cut ladder (y <= 10/3/1/0.5, all 78 kept, UB frozen)
GB-4  the Upsilon_bul ladder (0.55/0.60/0.70/0.80 warm + 0.55 cold):
      slope, Delta(-2lnL), fire-crossing
GB-5  the primary prior-strain census (SD(dv/sv), >2/>3 counts,
      SD(dml/S_ML), pinned)
GB-6  the sv ladder spots (k = 0.5, 2.0, 2.5) + local gradient
      adjudication
GB-7  the A5 correct-prior mock null (3 seeds)
GB-8  random-subset sigma spot checks (N = 66: 12 draws; N = 46: 10)
GB-9  census: the PGC 64824 row, SPARC f_D of the D-pair, 10T line-848
      text
GB-10 the coherent prior scenario (sv x1.5 AND upri x1.5)
"""
import json, math, os, re, sys, time, glob
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

t0 = time.time()
P("ROUND-50 GB (post-report)")

src = open('calcs/stage10t_legregrow.py', encoding='utf-8').read()
header = src[:src.index("# ---------------- G10T-1")]
ns = {'__name__': 'stage10t_inherited',
      '__file__': os.path.abspath('calcs/stage10t_legregrow.py')}
_argv = sys.argv
sys.argv = ['stage10t_legregrow.py', 'gates']
exec(compile(header, 'stage10t_legregrow.py[header]', 'exec'), ns)
sys.argv = _argv
make_world = ns['make_world']; build_sub = ns['build_sub']
fit_deep = ns['fit_deep']; fit_hier = ns['fit_hier']
nu_be = ns['nu_be']; h0_of_a0 = ns['h0_of_a0']
LN10 = ns['LN10']; U_PRIOR = ns['U_PRIOR']; S_ML = ns['S_ML']

W = make_world(lt_slot=None)
legA = [int(g) for g in W['legA']]
subT = build_sub(W, legA)
b_pl = fit_hier(subT, nu_be, use_u=False)
b_lo = fit_hier(subT, nu_be, use_u=True, th0=list(b_pl.x) + [0.0])
bT = fit_deep(subT, nu_be, th0=list(b_lo.x))
a0T = 10**bT.x[0]; fML = bT.x[1]
H0_T = h0_of_a0(a0T)
P(f"primary: H0 = {H0_T:.4f} (rev 65.3592), th = {np.round(bT.x,5)}")

P("")
P("-- GB-1 A1 bar arithmetic + order-statistic null --")
p2 = 2 * (1 - 0.5 * (1 + math.erf(2 / math.sqrt(2))))
pfire = 1 - (1 - p2)**78
P(f"  P(one draw > 2 sigma) = {p2:.4f}; P(any of 78 fires) = "
  f"{pfire:.4f} (rev 0.973)")
rng = np.random.default_rng(4242)
mx = np.max(np.abs(rng.normal(0, 0.5687, size=(20000, 78))), axis=1)
P(f"  simulated E[max|d|] = {np.mean(mx):.3f} (rev 1.515); 95th pct = "
  f"{np.percentile(mx, 95):.3f} (rev 1.935)")

P("")
P("-- GB-2 full LOO re-run (my own list) --")
t1 = time.time()
ds = []
for g in legA:
    bg = fit_deep(build_sub(W, [x for x in legA if x != g]), nu_be,
                  th0=list(bT.x))
    ds.append(h0_of_a0(10**bg.x[0]) - H0_T)
ds = np.array(ds)
sd = float(np.std(ds, ddof=1)); rms = float(np.sqrt(np.mean(ds**2)))
sig_jack = math.sqrt(77) * sd
neff = float((np.sum(ds**2))**2 / np.sum(ds**4))
mxobs = float(np.max(np.abs(ds)))
order = np.argsort(-np.abs(ds))
P(f"  max |d| = {mxobs:.3f} ({W['names'][legA[order[0]]]}); RMS = "
  f"{rms:.4f} (rev 0.6028); SD = {sd:.4f}")
P(f"  sigma_jack = sqrt(77)*SD = {sig_jack:.3f} (rev 5.290; bootstrap "
  f"4.99 -> ratio {sig_jack/4.99:.2f})")
P(f"  N_eff realized = {neff:.1f} (rev 15.1); top-1 share of Sum d^2 = "
  f"{ds[order[0]]**2/np.sum(ds**2)*100:.1f}% (rev 17.9)")
mx2 = np.max(np.abs(rng.normal(0, rms, size=(20000, 78))), axis=1)
P(f"  max percentile under homogeneous null at observed RMS: "
  f"{100*np.mean(mx2 < mxobs):.1f}% (his heterogeneous null: 87.5%)")
P(f"  [{(time.time()-t1)/60:.1f} min]")

P("")
P("-- GB-3 A4 point-cut ladder (y at the primary optimum, UB frozen) --")
gN_fix = W['gg'] + fML * W['gd'] + W['gb']
y_fix = gN_fix / a0T
for ycut, rev in ((10.0, -0.273), (3.0, -4.050), (1.0, -6.074),
                  (0.5, -6.745)):
    keep = y_fix <= ycut
    Wc = dict(W)
    Wc['gpts'] = {g: W['gpts'][g][keep[W['gpts'][g]]] for g in legA}
    mem = [g for g in legA if len(Wc['gpts'][g]) >= 1]
    b = fit_deep(build_sub(Wc, mem), nu_be, th0=list(bT.x))
    d = h0_of_a0(10**b.x[0]) - H0_T
    P(f"  y <= {ycut:5.2f}: pts {int(keep.sum())}/1263, gal {len(mem)}, "
      f"d = {d:+.3f} (rev {rev:+.3f}) {'OK' if abs(d-rev) < 0.05 else 'MISMATCH'}")

P("")
P("-- GB-4 Upsilon_bul ladder --")
gb0 = W['gb'].copy()
tab = {}
for ub in (0.55, 0.60, 0.70, 0.80):
    W['gb'] = gb0 * (ub / 0.7)
    b = fit_deep(build_sub(W, legA), nu_be, th0=list(bT.x))
    tab[ub] = (b.fun, h0_of_a0(10**b.x[0]))
    P(f"  UB {ub:.2f}: -2lnL = {b.fun:9.2f}, H0 = {tab[ub][1]:.3f}")
W['gb'] = gb0 * (0.55 / 0.7)
sc = build_sub(W, legA)
bc_pl = fit_hier(sc, nu_be, use_u=False)
bc_lo = fit_hier(sc, nu_be, use_u=True, th0=list(bc_pl.x) + [0.0])
bc = fit_deep(sc, nu_be, th0=list(bc_lo.x))
P(f"  UB 0.55 COLD chain: H0 = {h0_of_a0(10**bc.x[0]):.3f}, -2lnL = "
  f"{bc.fun:.2f} (warm {tab[0.55][0]:.2f})")
W['gb'] = gb0
slope = (tab[0.80][1] - tab[0.60][1]) / 0.20
dlnl = tab[0.70][0] - tab[0.55][0]
P(f"  slope (0.60-0.80) = {slope:+.1f} /unit (rev +12.5); "
  f"Delta(-2lnL) 0.70 vs 0.55 = {dlnl:+.2f} (rev +23.61)")
P(f"  H0(0.55) = {tab[0.55][1]:.2f} (rev 63.487; delta to primary "
  f"{tab[0.55][1]-H0_T:+.2f}, rev -1.87)")
# fire-crossing: A4 delta at UB = 0.55 vs the 4.2555 bar
bf66 = [g for g in legA if bool(np.all(gb0[W['gpts'][g]] == 0.0))]
W['gb'] = gb0 * (0.55 / 0.7)
b66 = fit_deep(build_sub(W, sorted(bf66)), nu_be, th0=list(bc.x))
d55 = h0_of_a0(10**b66.x[0]) - tab[0.55][1]
P(f"  A4 delta at UB 0.55 = {d55:+.3f} vs bar 4.2555 -> "
  f"{'no fire' if abs(d55) <= 4.2555 else 'fires'} (rev -3.997, no fire)")
W['gb'] = gb0

P("")
P("-- GB-5 primary prior-strain census --")
from scipy.optimize import minimize_scalar
la0, f, s_int, u = bT.x
n = subT['n']
dml = np.zeros(n); dv = np.zeros(n)
for _ in range(4):
    fac = f * np.exp(dml[subT['gidx']])
    gN = subT['gg'] + fac * subT['gd'] + subT['gb']
    r0 = (subT['lg'] - np.log10(gN * nu_be(gN / 10**la0))
          - u * subT['isuma_pt'])
    for gi in range(n):
        mm = subT['gpts'][gi]
        w = 1.0 / (subT['s2'][mm] + s_int**2)
        dv[gi] = np.sum(w * r0[mm]) / (np.sum(w) + 1.0 / subT['sv'][gi]**2)
    for gi in range(n):
        mm = subT['gpts'][gi]
        uoff = u * subT['isuma_g'][gi]
        def od(dl):
            fc = f * math.exp(dl)
            gN2 = subT['gg'][mm] + fc * subT['gd'][mm] + subT['gb'][mm]
            rr = (subT['lg'][mm] - np.log10(gN2 * nu_be(gN2 / 10**la0))
                  - dv[gi] - uoff)
            return (np.sum(rr * rr / (subT['s2'][mm] + s_int**2))
                    + dl * dl / (S_ML * S_ML))
        dml[gi] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                  method='bounded').x
z = dv / subT['sv']
P(f"  SD(dv/sv) = {np.std(z, ddof=1):.3f} (rev 1.696); |z|>2: "
  f"{int(np.sum(np.abs(z) > 2))}/78 (rev 19); |z|>3: "
  f"{int(np.sum(np.abs(z) > 3))}/78 (rev 9)")
P(f"  SD(dml/S_ML) = {np.std(dml/S_ML, ddof=1):.3f} (rev 1.194); "
  f"pinned |dml|~0.7: {int(np.sum(np.abs(np.abs(dml)-0.7) < 0.01))}/78 "
  f"(rev 3)")

P("")
P("-- GB-6 sv ladder spots + gradient --")
h0k = {1.0: H0_T}
for k, rev in ((0.5, 64.157), (1.5, 66.540), (2.0, 68.032),
               (2.5, 69.487)):
    s = build_sub(W, legA)
    s['sv'] = s['sv'] * k
    b = fit_deep(s, nu_be, th0=list(bT.x))
    h0k[k] = h0_of_a0(10**b.x[0])
    P(f"  k = {k}: H0 = {h0k[k]:.3f} (rev {rev})")
g_central = (h0k[1.5] - h0k[0.5]) / math.log(3)
g_fwd = (h0k[1.5] - h0k[1.0]) / math.log(1.5)
g_2 = (h0k[2.0] - h0k[1.0]) / math.log(2.0)
P(f"  gradient dH0/dln sv: central(0.5-1.5) = {g_central:+.2f}, "
  f"fwd(1-1.5) = {g_fwd:+.2f}, (1-2) = {g_2:+.2f}  [rev quoted +1.81 -- "
  f"ADJUDICATE]")

P("")
P("-- GB-7 A5 correct-prior mock null (3 seeds) --")
resp = []
for seed in (11, 12, 13):
    r = np.random.default_rng(seed)
    docc = [-math.log10(max(min(1.0 + r.normal(0, W['rel'][g]), 1.5),
                            0.5)) for g in legA]
    subm = build_sub(W, legA, dlg_per_occ=docc)
    bm = fit_deep(subm, nu_be, th0=list(bT.x))
    subm2 = build_sub(W, legA, dlg_per_occ=docc)
    subm2['sv'] = subm2['sv'] * 1.5
    bm2 = fit_deep(subm2, nu_be, th0=list(bm.x))
    resp.append(h0_of_a0(10**bm2.x[0]) - h0_of_a0(10**bm.x[0]))
    P(f"  seed {seed}: response = {resp[-1]:+.3f}")
P(f"  mean {np.mean(resp):+.3f} (rev +0.30 +- 0.37, 6/6 positive; "
  f"sign check: {sum(r_ > 0 for r_ in resp)}/3 positive)")

P("")
P("-- GB-8 random-subset sigma spot checks --")
for m, nd, rev in ((66, 12, 2.097), (46, 10, 5.236)):
    r = np.random.default_rng(500 + m)
    dd = []
    for _ in range(nd):
        pick = sorted(r.choice(legA, size=m, replace=False))
        b = fit_deep(build_sub(W, list(pick)), nu_be, th0=list(bT.x))
        dd.append(h0_of_a0(10**b.x[0]) - H0_T)
    P(f"  N = {m}: my {nd}-draw SD = {np.std(dd, ddof=1):.2f} (rev "
      f"{rev} at more draws); max |d| = {np.max(np.abs(dd)):.2f}")

P("")
P("-- GB-9 census verifications --")
for l in open('data/cf4/table2.dat'):
    try:
        pgc = int(l[0:7])
    except ValueError:
        continue
    if pgc == 64824:
        dm, edm = l[102:107].strip(), l[108:112].strip()
        nest = int(l[8:15])
        ra, dec = float(l[137:145]), float(l[146:154])
        pos = json.load(open('data/stage10x_pgc.json'))['KK98-251']
        cd = math.cos(math.radians(pos['dec']))
        sep = math.hypot((ra - pos['ra']) * cd, dec - pos['dec']) * 60
        D = 10**((float(dm) - 25) / 5)
        P(f"  PGC 64824: DMtrgb = {dm} +- {edm} -> D = {D:.2f} Mpc "
          f"(rev 6.98); 1PGC nest = {nest} (rev 65001); sep = "
          f"{sep:.2f} arcmin (rev 5.60)")
        break
for nm in ('D564-8', 'D631-7'):
    fdv = ns['meta'][nm][5]
    P(f"  SPARC f_D({nm}) = {fdv} (rev: 2 = TRGB for both)")
l848 = open('calcs/stage10t_legregrow.py', encoding='utf-8')\
    .read().split('\n')[847]
P(f"  stage10t line 848: {l848.strip()[:90]}")

P("")
P("-- GB-10 coherent prior scenario (sv AND upri x1.5) --")
s = build_sub(W, legA)
s['sv'] = s['sv'] * 1.5
bco = fit_hier(s, nu_be, use_u=True, upri=U_PRIOR * 1.5,
               th0=list(bT.x), tol=5e-4, max_rounds=60)
P(f"  sv x1.5 + upri x1.5: d = {h0_of_a0(10**bco.x[0]) - H0_T:+.3f} "
  f"(rev +1.364)")

P("")
P(f"GB total: {(time.time()-t0)/60:.1f} min")
with open('data/round50_gb.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(L) + "\n")
