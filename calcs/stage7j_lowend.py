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

# --- 7J-e2: the LITERATURE-ANCHORED conditional (scout follow-up) ---------
# Published subsystem rates for wide-binary components (scout-grade,
# primary-source verification pending: Tokovinin 2014 AJ 147, 86-87:
# 10.0%/7.3% per component; Tokovinin 2010: 12+-4%; Hwang 2022 field
# tertiary baseline ~5%) sit a factor ~2-3 BELOW part A's 0.24-0.34 per
# component (0.42-0.57 per pair).  If the literature is right, the
# per-pair host rate is ~0.2-0.25 and the relevant posterior cell is
# fcomp = 0.2.  Report that cell's conditional alpha and Newton margin
# (profile over all nuisances at fixed fcomp), and a marginal under a
# literature-centred prior (Gaussian on the fcomp axis, peak 0.22,
# sigma 0.08 - assumption-light stand-in pending the verified requote).
P("")
P("7J-e2: literature-anchored conditionals (scout-grade external rates)")
# round-4 recentring: Tokovinin per-component 10.0%/7.3% combines to a
# per-pair ~0.166, not the 0.22 first used - both centres reported
# ('lit' = 0.22 kept for continuity/regression; 'lit16' = 0.16 primary)
for sample in ('full', 'strict'):
    for cen, tag in ((0.22, 'lit'), (0.16, 'lit16')):
        for law in ('simple', 'BE'):
            for seed in (31, 101):
                cube = np.load(
                    f'data/stage7j_cube_{sample}_photo_{seed}_{law}.npy')
                cb = cube + prior_eta[None, :, None, None, None, None,
                                      None, None]
                sub = cb[:, :, :, 2:3]        # the fcomp = 0.2 cell
                prof = np.nanmax(sub, axis=(1, 2, 3, 4, 5, 6, 7))
                ima = int(np.nanargmax(prof)); ah = A_GRID[ima]
                if 0 < ima < 4:
                    x = A_GRID[ima-1:ima+2]; y = prof[ima-1:ima+2]
                    c2, c1, _ = np.polyfit(x, y, 2)
                    if c2 < 0: ah = -c1/(2*c2)
                dn = float(np.nanmax(prof) - prof[0])
                lit = -0.5*((FCOMP - cen)/0.08)**2
                cbl = cb + lit[None, None, None, :, None, None, None, None]
                m0 = np.nanmax(cbl)
                ex = np.exp(np.nan_to_num(cbl - m0, nan=-np.inf))
                lm = np.log(np.maximum(
                    ex.sum(axis=(1, 2, 3, 4, 5, 6, 7)), 1e-300)) + m0
                iml = int(np.argmax(lm)); am = A_GRID[iml]
                if 0 < iml < 4:
                    x = A_GRID[iml-1:iml+2]; y = lm[iml-1:iml+2]
                    c2, c1, _ = np.polyfit(x, y, 2)
                    if c2 < 0: am = -c1/(2*c2)
                dnm = float(lm.max() - lm[0])
                fp = ex.sum(axis=(0, 1, 2, 4, 5, 6, 7)); fp /= fp.sum()
                P(f"[{sample} {tag} {seed} {law}] cond(fcomp=0.2): "
                  f"a_prof={ah:.2f} "
                  f"dN={dn:+.1f} | lit-prior marg: a_marg={am:.2f} "
                  f"dN_marg={dnm:+.1f} P(fcomp)={np.round(fp, 2).tolist()}")

with open('data/stage7j_lowend.txt', 'w') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage7j_lowend.txt")
