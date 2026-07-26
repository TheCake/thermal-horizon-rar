# -*- coding: utf-8 -*-
"""STAGE 7J-z part 2 reader: the width-channel verdict at the landed
anchor (pre-registered before any photow cube existed).

INPUTS: photow cubes (full sample, seeds 31/101, both laws; 9-dim with
the SQ axis) + the part-1 landed anchor (data/stage7jz_prior.npz). If
part 1 did NOT ship (gates failed), the pre-registered fallback is the
literature anchor alone and the verdict is labeled LIT-CONDITIONAL.

ANCHORS READ (all three, landed = operative):
  LANDED   part-1 host-axis envelope, interpolated onto FCOMP_GRID
           exactly as 7J-B interpolated part A's (outside support ->
           ln pi = -1e9);
  LIT16    Gaussian -0.5*((f-0.16)/0.08)^2 (continuity; carries the
           per-component/per-pair convention note from part 1);
  MEAS     the RETRACTED part-A measured prior (continuity only).

CHANNEL DIAGNOSTICS (seed-law means, likelihood + eta prior only):
  D1 usage: P(sq > 0) >= 0.7 under the landed marginal -> CHANNEL-USED;
     P(sq = 0.3) >= 0.5 -> EDGE-FLAG (extension queued, correction-#4
     standard; the verdict is then sq-edge-conditional).
  D2 multiplicity-cost closure (the -60..-116 object, cube-readable):
     cost = max lnL - max lnL over fcomp >= 0.35 cells, computed
     identically on the photo (8-dim) and photow (sq-profiled) cubes;
     ratio photow/photo <= 0.5 -> CLOSED; >= 0.8 -> NOT-CLOSED.
  D3 fpm edge release: P(fpm = 2.4) < 0.5 under the landed marginal.
  Also reported: the PROF gain of sq-free over sq = 0 (the 6P +37 lnL
  comparator) and the joint (sq, fpm) posterior (the two noise channels
  may trade).

VERDICT BARS at the LANDED anchor (seed means, evaluated per law):
  BOOST-REVIVES  a_marg >= 0.5 AND dN_marg >= +25 (either law).
  NO-DETECTION   a_marg <= 0.3 AND dN_marg <= +10 (both laws).
  else AMBIGUOUS-CARRIED.
  Seed-extension rule (amendment 7d): |a_marg(31) - a_marg(101)| > 0.25
  at the landed anchor (either law) -> verdict deferred, seeds 202/303
  appended first.

AMENDMENT 9 (2026-07-26, review round 8, logged before any photow cube
existed to read): THE ANCHOR CURVE RE-READ. The 7J-e3 flatness
(a_marg 0.18-0.31, Newton +2..+4 over smooth anchors 0.06-0.30) was
measured under the OLD absorber configuration; the width channel
competes for the same width budget, so the flatness the landed-anchor
logic relies on is re-measured on the new cubes (same centers, step
0.02, extended to 0.34 to cover the interim MAP region; sigma 0.05 and
0.03; seed means per law). Reading rules (pre-registered):
  KNEE-REAPPEARS   any anchor in 0.06-0.30 gives a_marg >= 0.5 AND
                   dN_marg >= +25 (a detection-grade point exists on
                   the smooth curve under the new model);
  FLAT-PRESERVED   sigma=0.03 span (max-min a_marg, centers 0.06-0.30)
                   <= 0.25 AND all dN_marg <= +10, both laws;
  else INTERMEDIATE (reported as-is).
This tells us whether the flatness is a property of the data or of
the model that was running when it was measured (reviewer's phrasing,
adopted).
Output: data/stage7jz_read.txt
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
FPM = np.array([1.2, 1.5, 1.8, 2.1, 2.4])
KW = np.array([0.7, 1.0, 1.4])
SQ = np.array([0.0, 0.1, 0.2, 0.3])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
SEEDS = (31, 101)

anchors = {}
LANDED_OK = os.path.exists('data/stage7jz_prior.npz')
if LANDED_OK:
    pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
    lnp = np.full(len(FCOMP), -1e9)
    fg, lp = pz['fh_grid'], pz['lnpi_host']
    sup = lp > -1e8
    inr = (FCOMP >= fg[sup].min()) & (FCOMP <= fg[sup].max())
    lnp[inr] = np.interp(FCOMP[inr], fg, lp)
    anchors['LANDED'] = lnp
    P(f"landed anchor loaded: ln pi(fcomp) = {np.round(lnp, 2).tolist()}")
else:
    P("part-1 prior ABSENT -> LIT-CONDITIONAL mode (pre-registered "
      "fallback)")
anchors['LIT16'] = -0.5*((FCOMP-0.16)/0.08)**2
pa = np.load('data/stage7j_prior.npz')
lnm = np.full(len(FCOMP), -1e9)
inr = (FCOMP >= pa['f_grid'].min()) & (FCOMP <= pa['f_grid'].max())
lnm[inr] = np.interp(FCOMP[inr], pa['f_grid'], pa['lnpi_full'])
anchors['MEAS(retracted)'] = lnm
OPER = 'LANDED' if LANDED_OK else 'LIT16'

def read(cb9, lnpi):
    cbp = cb9 + lnpi.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp - m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=tuple(range(1, 9))), 1e-300)) + m0
    ima = int(np.argmax(lm)); am = A_GRID[ima]
    if 0 < ima < 4:
        x = A_GRID[ima-1:ima+2]; y = lm[ima-1:ima+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0: am = -c1/(2*c2)
    post = {}
    for name, ax, grid in (('wr', 2, WR_GRID), ('fcomp', 3, FCOMP),
                           ('fpm', 6, FPM), ('kw', 7, KW), ('sq', 8, SQ)):
        m = ex.sum(axis=tuple(i for i in range(9) if i != ax))
        post[name] = m/max(m.sum(), 1e-300)
    return am, float(lm.max()-lm[0]), lm, post

res = {}
for law in ('simple', 'BE'):
    for seed in SEEDS:
        cw = np.load(f'data/stage7j_cube_full_photow_{seed}_{law}.npy')
        cp = np.load(f'data/stage7j_cube_full_photo_{seed}_{law}.npy')
        assert np.isfinite(cw).all() and cw.ndim == 9, "photow cube bad"
        cb9 = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
        cb8 = cp + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1))
        # D2: multiplicity cost, photo vs photow (likelihood + eta only)
        cost_p = float(np.nanmax(cb8) - np.nanmax(cb8[:, :, :, 3:]))
        cost_w = float(np.nanmax(cb9) - np.nanmax(cb9[:, :, :, 3:]))
        gain_sq = float(np.nanmax(cb9) - np.nanmax(cb9[..., 0]))
        P(f"[{law} {seed}] D2 multiplicity cost: photo {cost_p:.1f} -> "
          f"photow {cost_w:.1f} (ratio {cost_w/max(cost_p,1e-9):.2f}); "
          f"PROF sq-gain {gain_sq:+.1f}")
        for aname, lnpi in anchors.items():
            am, dn, lm, post = read(cb9, lnpi)
            res[(law, seed, aname)] = (am, dn, post)
            P(f"[{law} {seed}] {aname}: a_marg={am:.2f} dN={dn:+.1f} "
              f"P(sq)={np.round(post['sq'],2).tolist()} "
              f"P(fpm)={np.round(post['fpm'],2).tolist()} "
              f"P(wr)={np.round(post['wr'],2).tolist()} "
              f"P(fcomp)={np.round(post['fcomp'],2).tolist()}")
        res[(law, seed, 'D2')] = (cost_p, cost_w, gain_sq)

P("")
ext_needed = False
for law in ('simple', 'BE'):
    d = abs(res[(law, 31, OPER)][0] - res[(law, 101, OPER)][0])
    if d > 0.25:
        ext_needed = True
        P(f"EXTENSION RULE FIRES ({law}): |a_marg(31)-a_marg(101)| = "
          f"{d:.2f} > 0.25 -> append seeds 202/303 before the verdict")
verdict_lines = []
for law in ('simple', 'BE'):
    ams = [res[(law, s, OPER)][0] for s in SEEDS]
    dns = [res[(law, s, OPER)][1] for s in SEEDS]
    psq = np.mean([res[(law, s, OPER)][2]['sq'][1:].sum() for s in SEEDS])
    pse = np.mean([res[(law, s, OPER)][2]['sq'][-1] for s in SEEDS])
    pfe = np.mean([res[(law, s, OPER)][2]['fpm'][-1] for s in SEEDS])
    cr = np.mean([res[(law, s, 'D2')][1] /
                  max(res[(law, s, 'D2')][0], 1e-9) for s in SEEDS])
    verdict_lines.append((law, np.mean(ams), np.mean(dns), psq, pse, pfe,
                          cr))
    P(f"{law} seed-mean @ {OPER}: a_marg={np.mean(ams):.2f} "
      f"dN={np.mean(dns):+.1f} | D1 P(sq>0)={psq:.2f} "
      f"edge P(sq=0.3)={pse:.2f} | D3 P(fpm=2.4)={pfe:.2f} | "
      f"D2 ratio={cr:.2f}")
if not ext_needed:
    rev = any(v[1] >= 0.5 and v[2] >= 25 for v in verdict_lines)
    nod = all(v[1] <= 0.3 and v[2] <= 10 for v in verdict_lines)
    v = ('BOOST-REVIVES' if rev else
         'NO-DETECTION' if nod else 'AMBIGUOUS-CARRIED')
    if not LANDED_OK:
        v += ' (LIT-CONDITIONAL)'
    d1 = 'CHANNEL-USED' if all(v_[3] >= 0.7 for v_ in verdict_lines) \
         else 'CHANNEL-QUIET'
    ef = any(v_[4] >= 0.5 for v_ in verdict_lines)
    d2 = ('CLOSED' if all(v_[6] <= 0.5 for v_ in verdict_lines) else
          'NOT-CLOSED' if any(v_[6] >= 0.8 for v_ in verdict_lines)
          else 'PARTIAL')
    d3 = all(v_[5] < 0.5 for v_ in verdict_lines)
    P(f"\n==> 7J-z part 2 VERDICT @ {OPER}: {v}  [D1 {d1}"
      f"{' +SQ-EDGE-FLAG' if ef else ''}; D2 {d2}; "
      f"D3 fpm-edge {'RELEASED' if d3 else 'STILL RIDING'}]")

# --- AMENDMENT 9: the anchor curve re-read on the photow cubes -----------
P("")
P("ANCHOR CURVE under the width channel (7J-e3 reference: a_marg "
  "0.18-0.31, dN +2..+4 over 0.06-0.30):")
CENS = np.arange(0.06, 0.3401, 0.02)
knee, spans, dmaxs = False, [], []
cubes = {(law, s): np.load(f'data/stage7j_cube_full_photow_{s}_{law}.npy')
         + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
         for law in ('simple', 'BE') for s in SEEDS}
for sg in (0.05, 0.03):
    for law in ('simple', 'BE'):
        row = []
        for cen in CENS:
            lnpi_c = -0.5*((FCOMP-cen)/sg)**2
            ams, dns = [], []
            for s in SEEDS:
                am, dn, _, = read(cubes[(law, s)], lnpi_c)[:3]
                ams.append(am); dns.append(dn)
            am_, dn_ = float(np.mean(ams)), float(np.mean(dns))
            row.append((cen, am_, dn_))
            if cen <= 0.301 and am_ >= 0.5 and dn_ >= 25:
                knee = True
        if sg == 0.03:
            in_ = [r for r in row if r[0] <= 0.301]
            spans.append(max(r[1] for r in in_) - min(r[1] for r in in_))
            dmaxs.append(max(r[2] for r in in_))
        P(f"  sg={sg} {law}: " + " ".join(
            f"{c:.2f}:{a:.2f}/{d:+.0f}" for c, a, d in row))
flat = (max(spans) <= 0.25) and (max(dmaxs) <= 10)
cv = ('KNEE-REAPPEARS' if knee else
      'FLAT-PRESERVED' if flat else 'INTERMEDIATE')
P(f"==> ANCHOR-CURVE READING: {cv} (sigma=0.03 spans "
  f"{[round(s,2) for s in spans]}, max dN {[round(d,1) for d in dmaxs]})")
with open('data/stage7jz_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("saved: data/stage7jz_read.txt")
