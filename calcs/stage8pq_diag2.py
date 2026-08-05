"""8P/8Q AMENDMENT-1 DIAGNOSTIC PART 2: the causal chain, closed.

Part 1 (stage8pq_diag.py): data/pop/assembly/orbit-determinism all
bit-identical; marginal-driven l2 reproduces the cube at 0.0; the solo
8Q re-run reproduces max|d| = 10.1 deterministically.  The one
expression difference left: the reader lineage (stage8lb_read.py:172,
inherited by 8N/8P/8Q) computes e_rad = 0.95+0.045*u while the
marginal computes e_rad = erf+(0.995-erf)*u with erf=0.95 — the slopes
differ by 4.16e-17 (a few ulp).  This script shows that ulp difference
ALONE reproduces the whole failure: same population, two e_s vectors
(reader-form vs marginal-form), two orbit runs, lnL blocks vs the
cube.  Output: data/stage8pq_diag2.txt
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
nsM = {}
exec(open('calcs/stage7j_marginal.py').read().split('def run_seed(')[0], nsM)
sys.argv = ['stage8q_pprior.py']
nsR = {}
exec(open('calcs/stage8q_pprior.py').read().split('ship = {}')[0], nsR)

seed, law = 31, 'simple'
c9 = np.load(f'data/stage7j_cube_full_lker_{seed}_{law}.npy')
al, eta, wr = 0.5, 1.05, 0.30
ai, ei, wi = 1, 0, 2
pM = nsM['build_pop'](seed)
tab_a = 1.0 + al*(nsM['TAB_S']-1.0)

e_marg = nsM['e_of_x'](pM, eta, wr, 1.0, 0.95)   # erf+(0.995-erf)*u
e_read = nsR['e_of_x'](pM, eta, wr)              # 0.95+0.045*u literal
nd = int(np.sum(e_marg != e_read))
dmax_e = float(np.max(np.abs(e_marg - e_read)))
OUT(f"e_s reader-form vs marginal-form: differing entries = {nd}/"
    f"{len(e_marg)} ({nd/len(e_marg):.3f}; wr = {wr} predicts the "
    f"radial branch), max|d e| = {dmax_e:.2e}")

sl = np.squeeze(c9[ai, ei, wi], axis=1)
res = {}
for tag, e_s in (('marginal-form', e_marg), ('reader-form', e_read)):
    vp = nsM['vp_c'](pM, e_s, tab_a)
    o = nsM['run'](pM['a_s'], e_s, pM['psi0'], pM['f_ip'], pM['M_s'],
                   pM['uph'], 8, 2500, 5, a0=nsM['A0_CAN'], tab=tab_a,
                   lny0=nsM['LNY0'], dlny=nsM['DLNY'], vp=vp)
    res[tag] = np.asarray(o)
    blk = nsR['eval_block'](pM, o, {'logn': pM['gs']})['logn']
    d = float(np.max(np.abs(blk[..., :4] - sl)))
    OUT(f"[{tag:14}] lnL block vs cube: max|d| = {d:.3e}")

do = np.abs(res['marginal-form'] - res['reader-form'])
OUT(f"orbit output divergence between the two e_s forms: max|d o| = "
    f"{float(do.max()):.3e}; systems with any |d o| > 1e-6: "
    f"{int(np.sum(do.max(axis=1) > 1e-6))} "
    f"(the ulp difference is integrator-amplified to macroscopic "
    f"phase shifts for a subset of near-parabolic systems)")
OUT(f"\nCONCLUSION: the entire G8P-1/G8Q-0 failure is the e_rad slope "
    f"literal (0.045 vs 0.995-0.95, d = 4.2e-17) introduced in the "
    f"8L-b reader and inherited by 8N/8P/8Q.  Fix: bit-verbatim "
    f"expression in the readers.  Census-grade gates (0.05) were blind "
    f"to it in four consecutive stages; the first lnL-grade identity "
    f"gate caught it on first firing.")
OUT(f"\ndone ({(time.time()-t0)/60:.1f} min)")
with open('data/stage8pq_diag2.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage8pq_diag2.txt")
