# -*- coding: utf-8 -*-
"""Stage 7J-z6 Part B reader — the width-shape contest at the LANDED-CONV
anchor (bars pre-registered at 66a8045, BEFORE any WSHAPE cube existed).

Reads whatever exists, in this order:
  1. Contest cubes  data/stage7j_cube_full_w{floor|tail}_31_{law}.npy
     (10-dim; ws axis last).  Per law:
       B1 SHAPE-WIN        profile-max(full) - profile-max(ws=0 slice)
                           >= +8.0 lnL, BOTH laws (one law = PARTIAL).
       B2 CHASE-DISSOLVED  P(fpm=3.0) <= 0.10 under the shape marginal
                           at LANDED-CONV.
       B3 alpha-STABILITY  alpha_marg in [0.55, 0.90] AND dN >= +10 ->
                           band-stable; < 0.55 or dN < 10 -> MATERIAL-LOW;
                           > 0.90 -> MATERIAL-HIGH (named, not protected).
  2. B4 arm cubes  data/stage7j_cube_fullarmw_photow3_31_{law}.npy
     (9-dim; the SHAPE-truth sky fit by the global-fpm machinery):
       P(fpm=3.0) >= 0.4 REPRODUCED / 0.1-0.4 PARTIAL / < 0.1 FAILED.
  3. GW1 arm cubes data/stage7j_cube_fullarmw_w{shape}_31_{law}.npy
     (10-dim; own-truth fit): shape param recovered within 1 grid step.
VERDICT (printed only when contest + B4 exist): RESOLVED = B1(both) AND
B2 AND B4 >= 0.4; PARTIAL = exactly one of B1/B2/B4 misses;
UNRESOLVED-CARRIED = neither shape clears B1 (the co-quoted band
0.68-0.74 / +14.5-23.8 stands; the item closes).
Output: data/stage7jz6_read.txt
"""
import os

import numpy as np

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FPM6 = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW = np.array([0.7, 1.0, 1.4])
SQ = np.array([0.0, 0.1, 0.2, 0.3])
WS_OF = {'floor': np.array([0.0, 0.015, 0.030, 0.045]),
         'tail':  np.array([0.0, 0.03, 0.08, 0.15])}
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2

# --- LANDED-CONV anchor (amendment 11 construction, verbatim) -------------
pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
assert str(pz['version']).startswith('v2c'), 'certificate npz expected'
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
gband = (pz['conv_band']/0.30 if 'conv_band' in pz.files
         else np.array([0.33, 1.30]))
GS = np.linspace(float(gband[0]), float(gband[1]), 25)
fhmin, fhmax = fg[sup].min(), fg[sup].max()
LNPI = np.full(len(FCOMP), -1e9)
for gi in GS:
    fh_eq = FCOMP/gi
    m = (fh_eq >= fhmin) & (fh_eq <= fhmax)
    cand = np.full(len(FCOMP), -1e9)
    cand[m] = np.interp(fh_eq[m], fg, lp)
    LNPI = np.maximum(LNPI, cand)
P(f"LANDED-CONV anchor: ln pi(fcomp) = {np.round(LNPI, 2).tolist()}")

def read_marg(cb, nax):
    """Marginal read with the completeness prior; cb has eta prior added.
    nax = cb.ndim.  Returns a_marg, dN, posteriors dict."""
    cbp = cb + LNPI.reshape((1, 1, 1, -1) + (1,)*(nax-4))
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp - m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=tuple(range(1, nax))), 1e-300)) + m0
    ima = int(np.argmax(lm)); am = A_GRID[ima]
    if 0 < ima < 4:
        x = A_GRID[ima-1:ima+2]; y = lm[ima-1:ima+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0: am = -c1/(2*c2)
    post = {}
    axmap = {'wr': 2, 'fcomp': 3, 'fpm': 6, 'kw': 7, 'sq': 8}
    if nax == 10:
        axmap['ws'] = 9
    for name, ax in axmap.items():
        m = ex.sum(axis=tuple(i for i in range(nax) if i != ax))
        post[name] = m/max(m.sum(), 1e-300)
    return am, float(lm.max()-lm[0]), post

def eta9(c):
    return c + prior_eta.reshape((1, 2) + (1,)*(c.ndim-2))

res = {}
for shape in ('floor', 'tail'):
    WS = WS_OF[shape]
    for law in ('simple', 'BE'):
        f = f'data/stage7j_cube_full_w{shape}_31_{law}.npy'
        if not os.path.exists(f):
            continue
        cw = np.load(f)
        assert cw.ndim == 10 and np.isfinite(cw).all(), f
        cb = eta9(cw)
        b1 = float(np.nanmax(cb) - np.nanmax(cb[..., 0]))
        am, dn, post = read_marg(cb, 10)
        p30 = float(post['fpm'][-1])
        wsm = WS[int(np.argmax(post['ws']))]
        res[(shape, law)] = (b1, am, dn, p30, post)
        P(f"[{shape} {law} 31] B1 shape-gain = {b1:+.2f} "
          f"({'>=+8 PASS' if b1 >= 8.0 else 'MISS'}); "
          f"a_marg={am:.2f} dN={dn:+.1f}; B2 P(fpm=3.0)={p30:.2f} "
          f"({'DISSOLVED' if p30 <= 0.10 else 'STILL-RIDING'}); "
          f"P(ws)={np.round(post['ws'], 2).tolist()} (mode ws={wsm}); "
          f"P(sq)={np.round(post['sq'], 2).tolist()} "
          f"P(fcomp)={np.round(post['fcomp'], 2).tolist()}")
        b3 = ('band-stable' if (0.55 <= am <= 0.90 and dn >= 10) else
              'MATERIAL-LOW' if (am < 0.55 or dn < 10) else 'MATERIAL-HIGH')
        P(f"[{shape} {law} 31] B3: {b3}")
        # seed 101 stability (winner only; reported, not a bar)
        f101 = f'data/stage7j_cube_full_w{shape}_101_{law}.npy'
        if os.path.exists(f101):
            cb1 = eta9(np.load(f101))
            b1s = float(np.nanmax(cb1) - np.nanmax(cb1[..., 0]))
            am1, dn1, post1 = read_marg(cb1, 10)
            P(f"[{shape} {law} 101] stability: gain {b1s:+.2f}, "
              f"a_marg={am1:.2f} dN={dn1:+.1f} "
              f"P(fpm=3.0)={post1['fpm'][-1]:.2f}")

# --- B4 arm (shape-truth sky, global-fpm fitter) --------------------------
b4 = {}
for law in ('simple', 'BE'):
    f = f'data/stage7j_cube_fullarmw_photow3_31_{law}.npy'
    if not os.path.exists(f):
        continue
    cb = eta9(np.load(f))
    am, dn, post = read_marg(cb, 9)
    p30 = float(post['fpm'][-1])
    b4[law] = p30
    grade = ('REPRODUCED' if p30 >= 0.4 else
             'PARTIAL' if p30 >= 0.1 else 'FAILED')
    P(f"[B4 arm {law}] global-fpm fit of shape-truth sky: "
      f"P(fpm=3.0)={p30:.2f} -> {grade}; a_marg={am:.2f} dN={dn:+.1f} "
      f"P(sq)={np.round(post['sq'], 2).tolist()}")

# --- exchangeability arm (amendment 2, bars pre-stated at ec7514c) --------
# The 7J-z5 arm-B injection (simple alpha=0.74 truth, twin-t5, sq=0.2,
# fpm=2.1, NO floor in truth) read by the floor fitter:
#   INFORMATIVE  a_marg(simple) >= 0.5 AND P(ws=0.045) <= 0.5
#   EATER        a_marg <= 0.3 AND P(ws=0.045) >= 0.5
#   else AMBIG.  Interprets B3; cannot flip B1/B2.
for law in ('simple', 'BE'):
    f = f'data/stage7j_cube_fullarmb_wfloor_31_{law}.npy'
    if not os.path.exists(f):
        continue
    cb = eta9(np.load(f))
    am, dn, post = read_marg(cb, 10)
    pwse = float(post['ws'][-1])
    p30 = float(post['fpm'][-1])
    line = (f"[XCHG arm {law}] boost-truth under floor fitter: "
            f"a_marg={am:.2f} dN={dn:+.1f} P(ws)="
            f"{np.round(post['ws'], 2).tolist()} P(fpm=3.0)={p30:.2f}")
    if law == 'simple':
        x = ('INFORMATIVE' if (am >= 0.5 and pwse <= 0.5) else
             'EATER' if (am <= 0.3 and pwse >= 0.5) else 'AMBIG')
        line += f" ==> {x}"
    P(line)

# --- GW1 arm (own-truth shape recovery) -----------------------------------
for shape in ('floor', 'tail'):
    WS = WS_OF[shape]
    for law in ('simple', 'BE'):
        f = f'data/stage7j_cube_fullarmw_w{shape}_31_{law}.npy'
        if not os.path.exists(f):
            continue
        cb = eta9(np.load(f))
        am, dn, post = read_marg(cb, 10)
        wsm = int(np.argmax(post['ws']))
        P(f"[GW1 arm {shape} {law}] own-truth read: P(ws)="
          f"{np.round(post['ws'], 2).tolist()} mode ws={WS[wsm]}; "
          f"a_marg={am:.2f}")

# --- PHYS conditional (amendment 3, logged BEFORE tail/arm results) -------
# The physical noise envelope: Lindegren+21 inflation allows fpm <= ~1.4
# (grid: <= 1.8) and an angular-covariance systematic floor <= ~0.025-0.03
# mas/yr = ws <= ~0.015 km/s at <= 200 pc.  Reported CONDITIONAL, not
# operative (changing the operative requires more than an interim round).
P("")
for law in ('simple', 'BE'):
    f3 = f'data/stage7j_cube_full_photow3_31_{law}.npy'
    if os.path.exists(f3):
        cb = eta9(np.load(f3))[:, :, :, :, :, :, :3]
        am, dn, post = read_marg(cb, 9)
        P(f"[PHYS {law}] photow3 | fpm<=1.8: a_marg={am:.2f} dN={dn:+.1f} "
          f"P(fpm)={np.round(post['fpm'], 2).tolist()}")
for shape in ('floor', 'tail'):
    for law in ('simple', 'BE'):
        f = f'data/stage7j_cube_full_w{shape}_31_{law}.npy'
        if not os.path.exists(f):
            continue
        cb = eta9(np.load(f))[:, :, :, :, :, :, :3, :, :, :2]
        am, dn, post = read_marg(cb, 10)
        P(f"[PHYS {law}] {shape} | fpm<=1.8, ws<=step1: a_marg={am:.2f} "
          f"dN={dn:+.1f} P(ws)={np.round(post['ws'], 2).tolist()}")

# --- verdict (contest + B4 present) ---------------------------------------
shapes_read = sorted({s for (s, _) in res})
if shapes_read:
    P("")
    for shape in shapes_read:
        if (shape, 'simple') in res and (shape, 'BE') in res:
            b1s = [res[(shape, law)][0] for law in ('simple', 'BE')]
            p30s = [res[(shape, law)][3] for law in ('simple', 'BE')]
            npass = sum(b >= 8.0 for b in b1s)
            b1v = 'PASS-BOTH' if npass == 2 else ('ONE-LAW' if npass == 1
                                                  else 'FAIL')
            b2v = 'PASS' if all(p <= 0.10 for p in p30s) else 'MISS'
            P(f"{shape}: B1 {b1v} ({b1s[0]:+.1f}/{b1s[1]:+.1f}), "
              f"B2 {b2v} (P30 {p30s[0]:.2f}/{p30s[1]:.2f})"
              + (f", B4 {min(b4.values()):.2f}" if len(b4) == 2 else
                 ", B4 pending"))
    winner = None
    for shape in shapes_read:
        if all((shape, law) in res and res[(shape, law)][0] >= 8.0
               for law in ('simple', 'BE')):
            winner = shape
            break
    if winner is None and all(
            all((s, law) in res for law in ('simple', 'BE'))
            for s in ('floor', 'tail')):
        P("==> 7J-z6 VERDICT: UNRESOLVED-CARRIED — neither shape clears "
          "B1; the co-quoted band 0.68-0.74 / +14.5-23.8 stands; the "
          "item closes (pre-registered closure rule)")
    elif winner and len(b4) == 2:
        b2ok = all(res[(winner, law)][3] <= 0.10 for law in ('simple', 'BE'))
        b4min = min(b4.values())
        misses = int(not b2ok) + int(b4min < 0.4)
        v = ('RESOLVED' if misses == 0 else
             'PARTIAL' if misses == 1 else 'UNRESOLVED-CARRIED')
        P(f"==> 7J-z6 VERDICT: {v} (winner {winner}; B2 "
          f"{'ok' if b2ok else 'miss'}; B4 min {b4min:.2f})")
    elif winner:
        P(f"==> 7J-z6: {winner} clears B1 — B4 arm pending before the "
          f"verdict")

with open('data/stage7jz6_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("saved: data/stage7jz6_read.txt")
