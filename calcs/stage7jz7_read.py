# -*- coding: utf-8 -*-
"""Stage 7J-z7 reader — THE TWIN-Q FORCED SCAN (bars pre-registered in
NOTES round 12, BEFORE any qt5 cube existed).

The round-12 question: is the FIFTH-MOVE collapse (forced fcomp >=
0.35 -> alpha 0.00, measured on flat-q cubes in 7J-z3) dead in the
CODE once the solver draws the measured twin-t5 q-law — or alive
regardless of the q-table arithmetic?

Per law x seed on the qt5 cubes (LANDED-CONV anchor):
  FORCED read: marginal restricted to fcomp >= 0.35 (the S3 grammar)
    -> alpha_forced, dN_forced;
  FREE read (descriptive): the unrestricted marginal on qt5;
  COST (descriptive): profile max(free) - max(forced) — the S3
    135-153 comparator under flat-q.
BARS: FIFTH-MOVE-DEAD-IN-CODE = alpha_forced >= 0.5 AND dN_forced >=
+10 in >= 3/4 reads; FIFTH-MOVE-ALIVE = alpha_forced <= 0.3 in >= 3/4
(-> MATERIAL, routed to the next decider); else PARTIAL.
EXPECTATION (stated, not a bar): host 0.35 under twin-t5 ~ kinematic
0.12-0.18 (the conversion band) -> alpha ~ 0.65-0.75 -> DEAD.
Output: data/stage7jz7_read.txt
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

def marg(cb, lnpi, f0):
    """Marginal over axes 1.. with the fcomp axis restricted to
    indices >= f0; returns (alpha_marg, dN)."""
    cbp = (cb[:, :, :, f0:] +
           lnpi[f0:].reshape((1, 1, 1, -1, 1, 1, 1, 1, 1)))
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp - m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=tuple(range(1, 9))), 1e-300)) + m0
    ima = int(np.argmax(lm)); am = A_GRID[ima]
    if 0 < ima < 4:
        x = A_GRID[ima-1:ima+2]; y = lm[ima-1:ima+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0: am = -c1/(2*c2)
    return am, float(lm.max()-lm[0])

res = []
for seed in (31, 101):
    for law in ('simple', 'BE'):
        f = f'data/stage7j_cube_full_qt5_{seed}_{law}.npy'
        if not os.path.exists(f):
            continue
        cw = np.load(f)
        assert cw.ndim == 9 and np.isfinite(cw).all(), f
        cb = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
        af, dnf = marg(cb, LNPI, 3)          # forced fcomp >= 0.35
        a0, dn0 = marg(cb, LNPI, 0)          # free
        cost = float(np.nanmax(cb) - np.nanmax(cb[:, :, :, 3:]))
        res.append((law, seed, af, dnf))
        P(f"[qt5 {law} {seed}] FORCED(fcomp>=0.35): alpha={af:.2f} "
          f"dN={dnf:+.1f} | FREE: alpha={a0:.2f} dN={dn0:+.1f} | "
          f"cost-to-force = {cost:.1f} (S3 flat-q comparator 135-153)")

if len(res) == 4:
    dead = sum(1 for _, _, a, d in res if a >= 0.5 and d >= 10.0)
    alive = sum(1 for _, _, a, _ in res if a <= 0.3)
    v = ('FIFTH-MOVE-DEAD-IN-CODE' if dead >= 3 else
         'FIFTH-MOVE-ALIVE (MATERIAL)' if alive >= 3 else 'PARTIAL')
    P(f"\n==> 7J-z7 VERDICT: {v} (dead-count {dead}/4, "
      f"alive-count {alive}/4; bars pre-registered round 12)")

with open('data/stage7jz7_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("saved: data/stage7jz7_read.txt")
