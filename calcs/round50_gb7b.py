"""GB-7b: the A5 correct-prior null done RIGHT (model-generated mock
skies, offsets drawn at exactly the prior widths; the first GB-7
construction added offsets on top of the real sky and is void for this
purpose -- disclosed). 4 seeds x (base fit + 1.5sv fit).
Appends to data/round50_gb.txt."""
import math, os, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

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
S_ML = ns['S_ML']; U_PRIOR = ns['U_PRIOR']

W = make_world(lt_slot=None)
legA = [int(g) for g in W['legA']]
subT = build_sub(W, legA)
b_pl = fit_hier(subT, nu_be, use_u=False)
b_lo = fit_hier(subT, nu_be, use_u=True, th0=list(b_pl.x) + [0.0])
bT = fit_deep(subT, nu_be, th0=list(b_lo.x))
la0, f, s_int, u = bT.x

P("GB-7b: A5 null on MODEL-GENERATED mocks (truth = the primary optimum)")
resp = []
for seed in (21, 22, 23, 24):
    r = np.random.default_rng(seed)
    dmlT = r.normal(0, S_ML, subT['n'])
    dvT = np.array([r.normal(0, s) for s in subT['sv']])
    uT = r.normal(0, U_PRIOR)
    fac = f * np.exp(dmlT[subT['gidx']])
    gN = subT['gg'] + fac * subT['gd'] + subT['gb']
    mu = (np.log10(gN * nu_be(gN / 10**la0)) + dvT[subT['gidx']]
          + uT * subT['isuma_pt'])
    lg_mock = mu + r.normal(0, np.sqrt(subT['s2'] + s_int**2))
    sm = build_sub(W, legA)
    sm['lg'] = lg_mock
    bb = fit_deep(sm, nu_be, th0=list(bT.x))
    sm2 = build_sub(W, legA)
    sm2['lg'] = lg_mock
    sm2['sv'] = sm2['sv'] * 1.5
    bw = fit_deep(sm2, nu_be, th0=list(bb.x))
    resp.append(h0_of_a0(10**bw.x[0]) - h0_of_a0(10**bb.x[0]))
    P(f"  seed {seed}: base {h0_of_a0(10**bb.x[0]):.2f}, response = "
      f"{resp[-1]:+.3f}")
P(f"  mean {np.mean(resp):+.3f} +- {np.std(resp, ddof=1)/2:.3f} "
  f"(rev +0.30 +- 0.15, SD 0.37); positive {sum(x > 0 for x in resp)}/4")
with open('data/round50_gb.txt', 'a', encoding='utf-8') as fo:
    fo.write("\n" + "\n".join(L) + "\n")
