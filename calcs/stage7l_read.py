# -*- coding: utf-8 -*-
"""Stage 7L L1/L2 reader — the cook-sample 2D likelihood at LANDED-CONV
(bars pre-registered in NOTES before any cook cube existed).

L1 CONSISTENCY: d = lm_max - max(lm at the alpha = 0.5 and 1.0 cells)
   <= 4 -> CONSISTENT; > 8 -> TENSION-NAMED; else GRAY.  dN(Newton)
   reported with the N-scaling expectation (~1194/14071 x the band =
   +1.2-2.0); NO detection claim on ~1.2k pairs.
L2 (descriptive): vt-only vs 2D on the cook sample.
Output: data/stage7l_read.txt
"""
import os

import numpy as np

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID = np.array([1.05, 1.3])
FCOMP = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2

pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
gband = (pz['conv_band']/0.30 if 'conv_band' in pz.files
         else np.array([0.33, 1.30]))
GS = np.linspace(float(gband[0]), float(gband[1]), 25)
LNPI = np.full(len(FCOMP), -1e9)
for gi in GS:
    fh_eq = FCOMP/gi
    m = (fh_eq >= fg[sup].min()) & (fh_eq <= fg[sup].max())
    cand = np.full(len(FCOMP), -1e9)
    cand[m] = np.interp(fh_eq[m], fg, lp)
    LNPI = np.maximum(LNPI, cand)

def read(cw):
    cb = cw + prior_eta.reshape((1, 2) + (1,)*(cw.ndim-2)) \
         + LNPI.reshape((1, 1, 1, -1) + (1,)*(cw.ndim-4))
    m0 = np.nanmax(cb)
    ex = np.exp(np.nan_to_num(cb - m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=tuple(range(1, cw.ndim))),
                           1e-300)) + m0
    ima = int(np.argmax(lm)); am = A_GRID[ima]
    if 0 < ima < 4:
        x = A_GRID[ima-1:ima+2]; y = lm[ima-1:ima+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0: am = -c1/(2*c2)
    return am, lm

for seed in (31, 101):
    for law in ('simple', 'BE'):
        f2 = f'data/stage7j_cube_cook_photow3_{seed}_{law}.npy'
        fv = f'data/stage7j_cubevt_cook_photow3_{seed}_{law}.npy'
        if not os.path.exists(f2):
            continue
        am, lm = read(np.load(f2))
        dn = float(lm.max()-lm[0])
        d_op = float(lm.max() - max(lm[1], lm[2]))
        l1 = ('CONSISTENT' if d_op <= 4.0 else
              'TENSION-NAMED' if d_op > 8.0 else 'GRAY')
        P(f"[cook {law} {seed}] 2D: a_marg={am:.2f} dN={dn:+.1f} "
          f"lm(alpha)-max = "
          f"{np.round(lm-lm.max(), 1).tolist()}; L1 d_op={d_op:.1f} "
          f"-> {l1}")
        if os.path.exists(fv):
            amv, lmv = read(np.load(fv))
            P(f"[cook {law} {seed}] vt-only: a_marg={amv:.2f} "
              f"dN={float(lmv.max()-lmv[0]):+.1f} (L2 descriptive)")

with open('data/stage7l_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("saved: data/stage7l_read.txt")
