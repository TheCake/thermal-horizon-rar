# -*- coding: utf-8 -*-
"""Stage 7J-d reader — the function contest under the landed posterior
(bars pre-registered in NOTES before any function cube existed).

Per function (photow3 9-dim cube, LANDED-CONV anchor):
  alpha_marg + dN(Newton);  PROF a_hat + interiority (edge-riding =
  shape rejection);  VOTE = ln-evidence(fn) - ln-evidence(BE)
  (logsumexp over cube + eta prior + LANDED-CONV fcomp prior + flat
  alpha grid — identical priors for every function); P(fpm) posterior.

BARS (locked in the pre-reg):
  VETO grammar: a previously vetoed member (rb4, boot, resn, dwf)
  STAYS vetoed iff PROF edge-riding OR vote <= -8; flips are named.
  LEAN grammar: |vote| < 5 tie-grade; 5-15 lean; > 15 strong lean.
  Lambda-family: landed c1_hat = lambda_hat/2 from the evidence
  profile, quoted with the lam000 rejection margin + per-seed spread.
  UNSUSPENSION: full set at seed 31 AND {p065, gm, rb1, rb2, rb3,
  amb} + lambda-family sign-stable at seed 101.
Output: data/stage7jd_read.txt
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

FNS = ['simple', 'BE', 'p065', 'gm', 'rb1', 'rb2', 'rb3', 'rb4',
       'boot', 'amb', 'resn', 'dwf',
       'lam000', 'lam025', 'lam050', 'lam075', 'lam100', 'lam125']
LAM = {'lam000': 0.0, 'lam025': 0.25, 'lam050': 0.5, 'lam075': 0.75,
       'lam100': 1.0, 'lam125': 1.25}
VETOED = ('rb4', 'boot', 'resn', 'dwf')
STAB = ('p065', 'gm', 'rb1', 'rb2', 'rb3', 'amb') + tuple(LAM)

# LANDED-CONV anchor (amendment-11 construction, verbatim)
pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
assert str(pz['version']).startswith('v2c')
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

def read(fn, seed):
    f = f'data/stage7j_cube_full_photow3_{seed}_{fn}.npy'
    if not os.path.exists(f):
        return None
    cw = np.load(f)
    assert cw.ndim == 9 and np.isfinite(cw).all(), f
    cb = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
    # PROF
    prof = np.nanmax(cb, axis=tuple(range(1, 9)))
    imax = int(np.nanargmax(prof)); ah = A_GRID[imax]
    if 0 < imax < 4:
        x = A_GRID[imax-1:imax+2]; y = prof[imax-1:imax+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0: ah = -c1/(2*c2)
    interior = 0 < imax < 4
    # MARG + evidence at LANDED-CONV
    cbp = cb + LNPI.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp - m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=tuple(range(1, 9))), 1e-300)) + m0
    ima = int(np.argmax(lm)); am = A_GRID[ima]
    if 0 < ima < 4:
        x = A_GRID[ima-1:ima+2]; y = lm[ima-1:ima+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0: am = -c1/(2*c2)
    ev = float(np.log(max(ex.sum(), 1e-300)) + m0)
    mfpm = ex.sum(axis=(0, 1, 2, 3, 4, 5, 7, 8))
    p30 = float(mfpm[-1]/max(mfpm.sum(), 1e-300))
    return dict(ah=ah, interior=interior, am=am,
                dn=float(lm.max()-lm[0]), ev=ev, p30=p30)

for seed in (31, 101):
    rows = {fn: read(fn, seed) for fn in FNS}
    if all(v is None for v in rows.values()):
        continue
    if rows.get('BE') is None:
        P(f"seed {seed}: BE reference cube missing — no votes")
        continue
    evbe = rows['BE']['ev']
    P(f"\n== 7J-d seed {seed} (LANDED-CONV anchor; vote = ev - ev_BE) ==")
    have = [(fn, r) for fn, r in rows.items() if r is not None]
    for fn, r in sorted(have, key=lambda t: -t[1]['ev']):
        vote = r['ev'] - evbe
        lean = ('tie' if abs(vote) < 5 else
                'lean' if abs(vote) < 15 else 'STRONG')
        edge = '' if r['interior'] else ' EDGE'
        P(f"  {fn:<7} vote={vote:+7.2f} ({lean})  a_marg={r['am']:.2f} "
          f"dN={r['dn']:+6.1f}  PROF a_hat={r['ah']:.2f}{edge}  "
          f"P(fpm=3.0)={r['p30']:.2f}")
    # veto stability
    for fn in VETOED:
        r = rows.get(fn)
        if r is None:
            continue
        vote = r['ev'] - evbe
        stays = (not r['interior']) or (vote <= -8.0)
        P(f"  VETO {fn}: {'STAYS' if stays else 'FLIP (named)'} "
          f"(edge={'yes' if not r['interior'] else 'no'}, "
          f"vote={vote:+.1f})")
    # lambda family -> landed c1
    lam_ev = {LAM[fn]: rows[fn]['ev'] for fn in LAM if rows.get(fn)}
    if len(lam_ev) == 6:
        lams = sorted(lam_ev)
        evs = np.array([lam_ev[l] for l in lams])
        pk = lams[int(np.argmax(evs))]
        P(f"  LAMBDA seed {seed}: ev-lam000 = "
          + " ".join(f"{l}:{lam_ev[l]-lam_ev[0.0]:+.1f}" for l in lams)
          + f" -> peak lambda={pk}, c1_hat={pk/2:.3f}; "
          f"c1=0 rejection = {evs.max()-lam_ev[0.0]:+.1f}")

with open('data/stage7jd_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage7jd_read.txt")
