# -*- coding: utf-8 -*-
"""Stage 7J-z3: adequacy decomposition of the width-complete model
(read-only cube diagnostic; 7J-c precedent; NO new fits, NO verdict,
NO credence move -- interim-round cadence rule applies).

CONTEXT (review round 9).  The reviewer asks: "does the -60..-116
residual misfit close with sq active?  If a substantial misfit remains,
there's another channel out there and alpha will move a fifth time."
The statistic he names ALREADY EXISTS as D2 in the shipped 7J-z part-2
read (pre-registered bars ratio<=0.5 CLOSED / >=0.8 NOT-CLOSED): it
fired NOT-CLOSED at ratio 1.37-1.42 (multiplicity cost 95-112 photo ->
135-153 photow).  What has NOT been computed is the decomposition that
makes the adequacy question decision-relevant:

  S1  anchor strain: G_free - lnL(posterior-mode cell under LIT16),
      both on the likelihood+eta surface.  If ~0, the quoted
      alpha_marg = 0.74/0.70 is the width-complete model's OWN
      unconstrained optimum -- the kinematic face of adequacy is
      closed and the entire residual misfit is the multiplicity
      tension, which is named and carried.
  S2  the free-optimum cell itself (alpha, fcomp, sq, fpm, wr at the
      argmax of likelihood+eta): what the model wants with no prior.
  S3  THE FIFTH-MOVE EXPOSURE: force fcomp >= 0.35 (the measured-host
      grid cell, D2's own slice) and re-profile alpha there.
      Pre-stated reading bands (bands chosen BEFORE the numbers were
      seen; readings, not verdict bars):
        alpha_hat(forced) >= 0.5 AND dN(forced) >= 15
            -> FIFTH-MOVE-SAFE: a future certificated host rate ~0.3
               would leave alpha standing; the D2 tension lives in the
               noise/nuisance sector, not under alpha.
        alpha_hat(forced) <= 0.25 OR dN(forced) <= 5
            -> FIFTH-MOVE-LIVE: v2c + the arm suite decide alpha's
               fate; printed as the standing risk in NOTES + PAPER.
        between -> PARTIAL (report both faces).

GATES (first run):
  G0  D2 regression: cost_p / cost_w reproduce the shipped values
      (109.1/153.3, 95.4/135.2, 111.6/153.4, 95.4/134.8) to 0.05 --
      the identity-point regression rule (GB0w precedent).
  G1  a_marg regression: the LIT16 marginal alpha reproduces the
      shipped 0.74/0.74/0.75/0.65 to 0.02.
  Any gate fail -> inspect, do not quote.

Also recorded here (round-9 phrasing corrections, no new numbers):
  - WR_GRID step is 0.10 and SQ_GRID step is 0.1: every shipped
    "SD = 0.000" means ALL POSTERIOR MASS ON ONE GRID NODE, i.e.
    SD below one grid step -- not infinite precision.
  - The constant SD(sq) = 0.927 in the sq0 rows of stage7jg_read.txt
    is sqrt(0.86): a broadcasting artifact of the pinned length-1 sq
    axis against the full length-4 grid (mu spuriously = sum(grid) =
    0.6).  Cosmetic only: the sq axis is PINNED in that config and no
    verdict metric (M1-M4) reads SD(sq); the sqfree SDs use
    length-matched axes and are sound.

Output: data/stage7j_sqclose.txt
"""
import numpy as np

OUT = []
def P(s):
    print(s); OUT.append(s)

A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FPM = np.array([1.2, 1.5, 1.8, 2.1, 2.4])
KW = np.array([0.7, 1.0, 1.4])
SQ = np.array([0.0, 0.1, 0.2, 0.3])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
lnpi_lit = -0.5*((FCOMP-0.16)/0.08)**2

SHIPPED = {('simple', 31): (109.1, 153.3, 0.74),
           ('simple', 101): (95.4, 135.2, 0.74),
           ('BE', 31): (111.6, 153.4, 0.75),
           ('BE', 101): (95.4, 134.8, 0.65)}

def refine(x, y):
    i = int(np.nanargmax(y)); v = x[i]
    if 0 < i < len(x)-1:
        c2, c1, _ = np.polyfit(x[i-1:i+2], y[i-1:i+2], 2)
        if c2 < 0:
            v = -c1/(2*c2)
    return v, i

g_ok = True
rows = []
for law in ('simple', 'BE'):
    for seed in (31, 101):
        cw = np.load(f'data/stage7j_cube_full_photow_{seed}_{law}.npy')
        cp = np.load(f'data/stage7j_cube_full_photo_{seed}_{law}.npy')
        cb9 = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
        cb8 = cp + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1))
        # G0: D2 regression (identical definition to stage7jz_read)
        cost_p = float(np.nanmax(cb8) - np.nanmax(cb8[:, :, :, 3:]))
        cost_w = float(np.nanmax(cb9) - np.nanmax(cb9[:, :, :, 3:]))
        sp, sw, sam = SHIPPED[(law, seed)]
        ok0 = abs(cost_p-sp) <= 0.05 and abs(cost_w-sw) <= 0.05
        # G1: LIT16 marginal alpha regression
        cbp = cb9 + lnpi_lit.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1))
        m0 = np.nanmax(cbp)
        ex = np.exp(np.nan_to_num(cbp-m0, nan=-np.inf))
        lm = np.log(np.maximum(ex.sum(axis=(1, 2, 3, 4, 5, 6, 7, 8)),
                               1e-300)) + m0
        am, _ = refine(A_GRID, lm)
        ok1 = abs(am-sam) <= 0.02
        g_ok = g_ok and ok0 and ok1
        P(f"[{law} {seed}] G0 D2 regression: photo {cost_p:.1f} "
          f"(shipped {sp}) photow {cost_w:.1f} (shipped {sw}) "
          f"{'PASS' if ok0 else 'FAIL'} | G1 a_marg {am:.2f} "
          f"(shipped {sam}) {'PASS' if ok1 else 'FAIL'}")
        # S1: anchor strain (posterior-mode cell vs free optimum,
        #     both evaluated on the likelihood+eta surface)
        idx = np.unravel_index(int(np.nanargmax(cbp)), cbp.shape)
        gap1 = float(np.nanmax(cb9) - cb9[idx])
        # S2: the free-optimum cell
        ifr = np.unravel_index(int(np.nanargmax(cb9)), cb9.shape)
        # S3: fifth-move exposure -- alpha profiled at fcomp >= 0.35
        sl = cb9[:, :, :, 3:]
        prof = np.nanmax(sl, axis=(1, 2, 3, 4, 5, 6, 7, 8))
        ah_f, i_f = refine(A_GRID, prof)
        dn_f = float(np.nanmax(prof) - prof[0])
        # and at the free-preferred fcomp for contrast
        pro0 = np.nanmax(cb9, axis=(1, 2, 3, 4, 5, 6, 7, 8))
        ah0, _ = refine(A_GRID, pro0)
        dn0 = float(np.nanmax(pro0) - pro0[0])
        rows.append((law, seed, gap1, ah_f, dn_f, i_f))
        P(f"[{law} {seed}] S1 anchor strain = {gap1:.1f} | "
          f"S2 free cell: a={A_GRID[ifr[0]]:.1f} wr={WR_GRID[ifr[2]]:.2f} "
          f"fcomp={FCOMP[ifr[3]]:.2f} fpm={FPM[ifr[6]]:.1f} "
          f"sq={SQ[ifr[8]]:.1f} | free a_hat={ah0:.2f} dN={dn0:+.1f}")
        P(f"[{law} {seed}] S3 FORCED fcomp>=0.35: a_hat={ah_f:.2f} "
          f"(interior={0 < i_f < 4}) dN={dn_f:+.1f}")

P("")
P(f"GATES: {'ALL PASS' if g_ok else 'FAIL -- do not quote'}")
if g_ok:
    g1m = np.mean([r[2] for r in rows])
    ahm = {law: np.mean([r[3] for r in rows if r[0] == law])
           for law in ('simple', 'BE')}
    dnm = {law: np.mean([r[4] for r in rows if r[0] == law])
           for law in ('simple', 'BE')}
    safe = all(np.mean([r[3] for r in rows if r[0] == law]) >= 0.5 and
               np.mean([r[4] for r in rows if r[0] == law]) >= 15
               for law in ('simple', 'BE'))
    live = any(np.mean([r[3] for r in rows if r[0] == law]) <= 0.25 or
               np.mean([r[4] for r in rows if r[0] == law]) <= 5
               for law in ('simple', 'BE'))
    tag = ('FIFTH-MOVE-SAFE' if safe else
           'FIFTH-MOVE-LIVE' if live else 'PARTIAL')
    P(f"S1 mean anchor strain = {g1m:.1f} (near 0 -> the quoted "
      f"alpha is the width-complete model's own optimum)")
    P(f"S3 seed-means: simple a_hat(forced)={ahm['simple']:.2f} "
      f"dN={dnm['simple']:+.1f}; BE a_hat(forced)={ahm['BE']:.2f} "
      f"dN={dnm['BE']:+.1f}")
    P(f"==> READING (pre-stated bands): {tag}")

with open('data/stage7j_sqclose.txt', 'w') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage7j_sqclose.txt")
