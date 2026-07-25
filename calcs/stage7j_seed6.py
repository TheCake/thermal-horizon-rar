# -*- coding: utf-8 -*-
"""Stage 7J-s6: the six-seed budget read (measurements; cadence rule).

Amendment-6 standing rule (0045eea): the measured-prior verdict stands
unless any NEW seed (202/303/404/505) shows the boundary-free break -
marginal gap lm(0.5) - lm(0) > -5 lnL, or an interior a_marg > 0 - at
the peak prior.  This reader evaluates that rule and reports the
six-seed statistics at BOTH anchorings (the retracted measured prior,
kept as the continuity row, and the literature prior 0.16 = the
operative one).  Per the cadence rule these are measurements; the
verdict labels stay frozen until 7J-z / 7J-g.

Output: data/stage7j_seed6.txt
"""
import numpy as np

OUT = []
def P(s):
    print(s); OUT.append(s)

E_GRID = np.array([1.05, 1.3])
A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
FCOMP = np.array([0.0, 0.1, 0.2, 0.35, 0.5, 0.7])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
LNPI = {'full':   np.array([-1584.82, -768.17, -268.10, -7.78, -1.06, -42.57]),
        'strict': np.array([-1584.82, -608.39, -150.50, -2.45, -3.38, -68.11])}
LIT16 = -0.5*((FCOMP - 0.16)/0.08)**2
SEEDS = (31, 101, 202, 303, 404, 505)

def read(cb, lnpi):
    cbp = cb + lnpi[None, None, None, :, None, None, None, None]
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp - m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=(1, 2, 3, 4, 5, 6, 7)), 1e-300)) + m0
    ima = int(np.argmax(lm)); am = A_GRID[ima]
    if 0 < ima < 4:
        x = A_GRID[ima-1:ima+2]; y = lm[ima-1:ima+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0: am = -c1/(2*c2)
    return am, float(lm[1] - lm[0]), float(lm.max() - lm[0])

breaks = []
for sample in ('full', 'strict'):
    for law in ('simple', 'BE'):
        amp, gp, aml, gl, dnl = [], [], [], [], []
        for seed in SEEDS:
            cube = np.load(
                f'data/stage7j_cube_{sample}_photo_{seed}_{law}.npy')
            cb = cube + prior_eta[None, :, None, None, None, None, None, None]
            a1, g1, _ = read(cb, LNPI[sample])
            a2, g2, d2 = read(cb, LIT16)
            amp.append(a1); gp.append(g1)
            aml.append(a2); gl.append(g2); dnl.append(d2)
            if seed >= 202 and (g1 > -5 or a1 > 0):
                breaks.append((sample, law, seed, a1, g1))
        P(f"[{sample} {law}] measured prior: a_marg per seed "
          f"{[f'{v:.2f}' for v in amp]}  gap(0.5-0) "
          f"{[f'{v:+.1f}' for v in gp]}")
        P(f"[{sample} {law}] lit-0.16:      a_marg per seed "
          f"{[f'{v:.2f}' for v in aml]}  gap "
          f"{[f'{v:+.1f}' for v in gl]}  -> mean a_marg = "
          f"{np.mean(aml):.2f} +- {np.std(aml)/np.sqrt(6):.2f} SE, "
          f"mean dN = {np.mean(dnl):+.1f}")
P("")
if breaks:
    P("AMENDMENT-6 SEED RULE: BREAK -> AMBIGUOUS at the measured prior:")
    for b in breaks:
        P(f"  {b[0]} {b[1]} seed {b[2]}: a_marg={b[3]:.2f} gap={b[4]:+.1f}")
else:
    P("AMENDMENT-6 SEED RULE: no break in 16 new-seed reads (all gaps")
    P("<= -5 and a_marg = 0 at the measured prior) -> the measured-prior")
    P("verdict STANDS on the seed axis (its prior remains retracted; the")
    P("operative statement is the lit-anchor row above + the flat curve).")

with open('data/stage7j_seed6.txt', 'w') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage7j_seed6.txt")
