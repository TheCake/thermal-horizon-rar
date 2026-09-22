"""
ROUND-51 M1 FIX: the peg-gearing measurement the record never made --
the response of the ANCHORED arm's PRIMARY (hier deep) estimator to a
coherent distance rescale (executed as a round condition, addendum
precedent; the record's only prior responses are flow-flat -1.585 and
flow-hier 0.86, neither on the anchored arm).

Construction (the 10S G3 convention, pin-3): a coherent peg rescale
D -> sD shifts every member's lgobs by -log10(s) (g_bar invariant);
the priors are NOT rescaled (a peg error is unknown to the fitter).
gamma = -d ln a0 / d ln s; the ladder's response is 1 by construction.
Fits: hier deep warm at s = 0.92/1.00/1.08 + one cold chain at 1.08
(warm-start control) + the flat co-read at the same s (treatment
contrast). Output: data/round51_gearing.txt
"""
import math, os, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

t0 = time.time()
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
fit_flat = ns['fit_flat_pts']
nu_be = ns['nu_be']; h0_of_a0 = ns['h0_of_a0']

W = make_world(lt_slot=None)
legA = [int(g) for g in W['legA']]
subT = build_sub(W, legA)
b_pl = fit_hier(subT, nu_be, use_u=False)
b_lo = fit_hier(subT, nu_be, use_u=True, th0=list(b_pl.x) + [0.0])
bT = fit_deep(subT, nu_be, th0=list(b_lo.x))
P("ROUND-51 GEARING: anchored arm (78), coherent rescale, hier PRIMARY")
P(f"  s = 1.00: a0 = {10**bT.x[0]:.4e}, H0 = {h0_of_a0(10**bT.x[0]):.3f}, "
  f"f_ML = {bT.x[1]:.4f}")

res = {1.00: (10**bT.x[0], bT.x[1])}
for s in (0.92, 1.08):
    dlg = [-math.log10(s)] * len(legA)
    sub = build_sub(W, legA, dlg_per_occ=dlg)
    b = fit_deep(sub, nu_be, th0=list(bT.x))
    res[s] = (10**b.x[0], b.x[1])
    P(f"  s = {s:.2f}: a0 = {res[s][0]:.4e}, H0 = "
      f"{h0_of_a0(res[s][0]):.3f}, f_ML = {res[s][1]:.4f}")
# cold-chain control at s = 1.08
dlg = [-math.log10(1.08)] * len(legA)
subc = build_sub(W, legA, dlg_per_occ=dlg)
bc_pl = fit_hier(subc, nu_be, use_u=False)
bc_lo = fit_hier(subc, nu_be, use_u=True, th0=list(bc_pl.x) + [0.0])
bc = fit_deep(subc, nu_be, th0=list(bc_lo.x))
P(f"  s = 1.08 COLD chain: a0 = {10**bc.x[0]:.4e} (warm "
  f"{res[1.08][0]:.4e}; rel d = {abs(10**bc.x[0]/res[1.08][0]-1):.2e})")

g_lo = -math.log(res[0.92][0] / res[1.00][0]) / math.log(0.92)
g_hi = -math.log(res[1.08][0] / res[1.00][0]) / math.log(1.08)
g_ch = -math.log(res[1.08][0] / res[0.92][0]) / math.log(1.08 / 0.92)
P(f"  GAMMA (hier, anchored): chord 0.92-1.00 = {g_lo:.3f}; "
  f"1.00-1.08 = {g_hi:.3f}; two-sided = {g_ch:.3f}")
gf = {}
for s in (0.92, 1.00, 1.08):
    dlg_pt = -math.log10(s)
    lgs = subT['lg'] + dlg_pt
    la0f, _fdd = fit_flat(lgs, subT['gg'], subT['gd'], subT['gb'], nu_be)
    gf[s] = 10**la0f
gflat = -math.log(gf[1.08] / gf[0.92]) / math.log(1.08 / 0.92)
P(f"  flat co-read (same arm, same s): a0 ratio -> gamma = {gflat:.3f}")
P(f"  f_ML absorption: {res[0.92][1]:.4f} / {res[1.00][1]:.4f} / "
  f"{res[1.08][1]:.4f} (Newtonian-arm share of the shift)")
P(f"  record context: flow-flat gamma = 1.585 (10S G3); flow-hier "
  f"0.86 (10S leg-B grid); ladder = 1 by construction")
P(f"total: {(time.time()-t0)/60:.1f} min")
with open('data/round51_gearing.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(L) + "\n")
