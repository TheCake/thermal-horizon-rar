# -*- coding: utf-8 -*-
"""Stage 7J-e: low-end host-prior stress (external-review request).

The reviewer's point: f_host enters the verdict through the prior, and
alpha_marg is pinned at a bound - if the prior peak is high by ~20%, the
bound is where alpha goes regardless of the data.  This reader re-runs
the marginalization from the CACHED photo cubes with the measured
envelope prior recentred at its 1-sigma LOW end (full: peak 0.51 ->
0.42; strict: 0.47 -> 0.37), shape preserved (the stored envelope
evaluated at f + shift).  Read-only; the fired verdict is not touched;
this is the pre-registered sensitivity companion to it.

Report per (sample, law, seed): a_marg, dN_marg, posterior fcomp cell,
and the alpha-marginal gap lm(0.5) - lm(0) (the bound-free statistic).
Output: data/stage7j_lowend.txt
"""
import numpy as np

OUT = []
def P(s):
    print(s); OUT.append(s)

E_GRID = np.array([1.05, 1.3])
A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
FCOMP = np.array([0.0, 0.1, 0.2, 0.35, 0.5, 0.7])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
pr = np.load('data/stage7j_prior.npz')
SHIFT = {'full': 0.09, 'strict': 0.10}   # peak 0.51->0.42 / 0.47->0.37

def lnpi_for(sample, shift):
    xg = pr['f_grid'] if sample == 'full' else pr['r_grid']
    lp = pr['lnpi_full'] if sample == 'full' else pr['lnpi_strict']
    ln = np.full(len(FCOMP), -1e9)
    xs = FCOMP + shift
    inr = (xs >= xg.min()) & (xs <= xg.max())
    ln[inr] = np.interp(xs[inr], xg, lp)
    return ln

for sample in ('full', 'strict'):
    for tag, sh in (('peak', 0.0), ('low', SHIFT[sample])):
        lnpi = lnpi_for(sample, sh)
        P(f"[{sample} {tag}] ln pi(fcomp) = {np.round(lnpi, 2).tolist()}")
        for law in ('simple', 'BE'):
            for seed in (31, 101):
                cube = np.load(
                    f'data/stage7j_cube_{sample}_photo_{seed}_{law}.npy')
                cb = (cube
                      + prior_eta[None, :, None, None, None, None, None, None]
                      + lnpi[None, None, None, :, None, None, None, None])
                m0 = np.nanmax(cb)
                ex = np.exp(np.nan_to_num(cb - m0, nan=-np.inf))
                lm = np.log(np.maximum(
                    ex.sum(axis=(1, 2, 3, 4, 5, 6, 7)), 1e-300)) + m0
                ima = int(np.argmax(lm)); am = A_GRID[ima]
                if 0 < ima < 4:
                    x = A_GRID[ima-1:ima+2]; y = lm[ima-1:ima+2]
                    c2, c1, _ = np.polyfit(x, y, 2)
                    if c2 < 0: am = -c1/(2*c2)
                dnm = float(lm.max() - lm[0])
                fpost = ex.sum(axis=(0, 1, 2, 4, 5, 6, 7))
                fpost /= fpost.sum()
                gap = float(lm[1] - lm[0])
                P(f"[{sample} {tag} {seed} {law}] a_marg={am:.2f} "
                  f"dN_marg={dnm:+.1f} gap(0.5-0)={gap:+.1f} "
                  f"P(fcomp)={np.round(fpost, 2).tolist()}")
P("")
P("Bar reference (unchanged, from the 7J pre-reg): full COMPANION-WIN if")
P("a_marg <= 0.7 or dN <= +15 (either law, seed mean).")

with open('data/stage7j_lowend.txt', 'w') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage7j_lowend.txt")
