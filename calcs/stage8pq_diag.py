"""8P/8Q AMENDMENT-1 DIAGNOSTIC (pre-quote; runs preserved).

G8P-1/G8Q-0 cube identity FAILED at 2.6-10.1 lnL (law-seed-dependent,
deterministic across both scripts) while census identity PASSES exactly
(4/4) — the first-ever direct reader-vs-cube lnL comparison in the
program.  This script isolates the layer: it drives stage7j_marginal's
OWN build_pop + lnL_point at the (simple, 31) MAP config and compares
  (a) marginal l2 vs the stored lker cube slice   (code+data drift?)
  (b) marginal population arrays vs the 8N-lineage reader population
      (stream divergence?)
  (c) reader eval_block ON the marginal population + SAME orbits vs
      marginal l2                                  (assembly parity?)
  (d) data-side objects (data_2d / noise_pool / templates) marginal vs
      reader                                       (data-block parity?)
No stage quantities are produced; the record feeds the amendment.
Output: data/stage8pq_diag.txt
"""
import os, sys, time
import numpy as np

t0 = time.time()
L_ = []
def OUT(s):
    print(s, flush=True)
    L_.append(s)

os.environ['FPME'] = '1'
os.environ['LKER'] = '1'
sys.argv = ['stage7j_marginal.py', 'full', 'photow', '31']
srcM = open('calcs/stage7j_marginal.py').read()
headM = srcM.split('def run_seed(')[0]
nsM = {}
exec(headM, nsM)
OUT(f"[{(time.time()-t0)/60:.1f} min] marginal head loaded "
    f"(SAMPLE={nsM['SAMPLE']}, AMP={nsM['AMP']}, FPME={nsM['FPME']}, "
    f"LKER={nsM['LKER']}, N={nsM['N']}, "
    f"FPM_GRID={nsM['FPM_GRID'].tolist()}, "
    f"SQ_GRID={nsM['SQ_GRID'].tolist()})")

sys.argv = ['stage8q_pprior.py']
srcR = open('calcs/stage8q_pprior.py').read()
headR = srcR.split('ship = {}')[0]
nsR = {}
exec(headR, nsR)
OUT(f"[{(time.time()-t0)/60:.1f} min] reader head loaded")

# ---- (d) data-side parity ------------------------------------------
for bi in range(4):
    OUT(f"bin{bi}: data_2d identical = "
        f"{np.array_equal(nsM['data_2d'][bi], nsR['data_2d'][bi])}; "
        f"noise_pool identical = "
        f"{np.array_equal(nsM['noise_pool'][bi], nsR['noise_pool'][bi])} "
        f"(len M={len(nsM['noise_pool'][bi])}, R={len(nsR['noise_pool'][bi])}); "
        f"UNI_B identical = "
        f"{np.array_equal(nsM['UNI_B'][bi], nsR['UNI_B'][bi])}; "
        f"FLY_B identical = "
        f"{np.array_equal(nsM['FLY_B'][bi], nsR['FLY_B'][bi])}")
OUT(f"SC2 identical = {np.array_equal(nsM['SC2'], nsR['SC2'])}; "
    f"VE identical = {np.array_equal(nsM['VE'], nsR['VE'])}; "
    f"GE identical = {np.array_equal(nsM['GE'], nsR['GE'])}")

# ---- (b) population parity -----------------------------------------
seed, law = 31, 'simple'
pM = nsM['build_pop'](seed)
pR = nsR['build_pop'](seed, 'raghavan')
OUT(f"[{(time.time()-t0)/60:.1f} min] populations built")
for k_ in ['a_s','u_e','psi0','f_ip','M_s','uph','u_mix','gn1','gn2','gs']:
    OUT(f"pop[{k_:5}]: "
        f"{'IDENTICAL' if np.array_equal(pM[k_], pR[k_]) else 'DIFFERS'}")
for kk in (1, 2):
    for f_ in ['w', 'uc', 'mh']:
        OUT(f"comp{kk}[{f_:2}]: "
            f"{'IDENTICAL' if np.array_equal(pM['comp'][kk][f_], pR['comp'][kk][f_]) else 'DIFFERS'}")
OUT(f"pick: "
    + str(['IDENTICAL' if np.array_equal(pM['pick'][b], pR['pick'][b])
           else 'DIFFERS' for b in range(4)]))

# ---- (a) marginal-driven l2 vs the stored cube ---------------------
c9 = np.load(f'data/stage7j_cube_full_lker_{seed}_{law}.npy')
al, eta, wr = 0.5, 1.05, 0.30
ai, ei, wi = 1, 0, 2
tab_a = 1.0 + al*(nsM['TAB_S']-1.0)
e_s = nsM['e_of_x'](pM, eta, wr, 1.0, 0.95)
vp = nsM['vp_c'](pM, e_s, tab_a)
o = nsM['run'](pM['a_s'], e_s, pM['psi0'], pM['f_ip'], pM['M_s'],
               pM['uph'], 8, 2500, 5, a0=nsM['A0_CAN'], tab=tab_a,
               lny0=nsM['LNY0'], dlny=nsM['DLNY'], vp=vp)
OUT(f"[{(time.time()-t0)/60:.1f} min] orbits done; running lnL_point...")
l2, lv = nsM['lnL_point'](pM, o)
l2 = l2[:, :, :, :, :, 0]     # wcut singleton
l2 = l2[..., 0]               # ws singleton
sl = c9[ai, ei, wi]           # (6,1,2,6,3,4)
dA = float(np.max(np.abs(l2 - sl)))
OUT(f"[{(time.time()-t0)/60:.1f} min] (a) marginal-driven l2 vs stored "
    f"cube: max|d| = {dA:.3e}"
    + ("  == CODE+DATA UNCHANGED (cube reproducible)" if dA <= 1e-9
       else "  == CUBE NOT REPRODUCIBLE BY TODAY'S CODE+DATA"))

# ---- (c) reader assembly on the marginal pop + SAME orbits ---------
blk = nsR['eval_block'](pM, o, {'logn': pM['gs']})['logn']
dC = float(np.max(np.abs(blk[..., :4] - np.squeeze(sl, axis=1))))
OUT(f"(c) reader eval_block on marginal pop + same orbits vs cube "
    f"slice: max|d| = {dC:.3e}")
dC2 = float(np.max(np.abs(blk[..., :4] - np.squeeze(l2, axis=1))))
OUT(f"(c2) reader eval_block vs marginal l2 (same pop, same orbits): "
    f"max|d| = {dC2:.3e}"
    + ("  == ASSEMBLY PARITY EXACT" if dC2 <= 1e-9
       else "  == ASSEMBLY DIVERGES"))

# ---- orbit determinism probe: second run() call --------------------
o2 = nsM['run'](pM['a_s'], e_s, pM['psi0'], pM['f_ip'], pM['M_s'],
                pM['uph'], 8, 2500, 5, a0=nsM['A0_CAN'], tab=tab_a,
                lny0=nsM['LNY0'], dlny=nsM['DLNY'], vp=vp)
dO = float(np.max(np.abs(np.asarray(o) - np.asarray(o2))))
OUT(f"orbit determinism (two identical run() calls): max|d o| = "
    f"{dO:.3e}"
    + ("  == DETERMINISTIC" if dO == 0.0 else "  == NON-DETERMINISTIC"))
if dO > 0.0:
    l2b, _ = nsM['lnL_point'](pM, o2)
    l2b = l2b[:, :, :, :, :, 0][..., 0]
    dR = float(np.max(np.abs(l2b - l2)))
    OUT(f"lnL spread from orbit non-determinism alone (same pop, two "
        f"run() calls): max|d lnL| = {dR:.3e}")

OUT(f"\ndone ({(time.time()-t0)/60:.1f} min)")
with open('data/stage8pq_diag.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage8pq_diag.txt")
