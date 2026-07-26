# -*- coding: utf-8 -*-
"""STAGE 7J-g: the gamma-collapse separability test — THE PAPER-1
DECIDER (bars locked before any vt cube is read).

QUESTION (the four-absorber diagnosis, review rounds 4-6): boost /
near-parabolic orbits (wr) / hidden companions (fcomp) / noise
inflation (fpm, sq) share ONE width budget in the scalar-velocity
channel; only the DIRECTION channel (gamma) can separate them. Paper
1's methodological claim is that the 2D (vtilde x gamma) likelihood
does that separation. Here the claim is measured: the photow cubes and
their gamma-collapsed twins (same model, same populations, gamma summed
out of data and model) are read at the same anchor and the information
attributable to gamma is quantified.

CONFIGURATION: full sample, seeds 31/101, both laws, photow cubes +
cubevt twins, eta prior, operative anchor = LANDED (part 1's measured
prior; pre-registered fallback LIT16 if part 1 did not ship).

PER (law, seed, channel): a_marg, dN_marg (internal to each channel —
absolute lnL levels differ by binning and are never compared), the
alpha 1-sigma set width W_a (piecewise-linear interp of the alpha
marginal on a 0.01 grid, measure of {lm >= max - 0.5}; crude on a
5-point grid but IDENTICAL procedure in both channels), posterior SDs
of (wr, fcomp, sq, fpm), and posterior corr(alpha, theta).

PRE-REGISTERED METRICS (seed means; "either law" = fires if any law):
  M1 alpha information:  shrink = 1 - W_a(2D)/W_a(vt) >= 0.30
  M2 gamma evidence:     dN_marg(2D) - dN_marg(vt) >= +15
  M3 answer shift:       |a_marg(2D) - a_marg(vt)| >= 0.30
  M4 absorber grade:     posterior SD of wr or fcomp shrinks >= 30%
                         from vt to 2D
THREE-TIER VERDICT:
  SEPARATION-CONFIRMED (alpha grade)  if M1 or M2 or M3 fires.
  ABSORBER-LEVEL SEPARATION           if only M4 fires: gamma pins the
      absorbers while alpha-inference is gamma-flat at this anchor —
      the 2D claim operates at nuisance grade (stated so in Paper 1).
  SEPARATION-ABSENT  if none fires AND all alpha metrics are quiet
      (shrink < 0.15, |shift| < 0.15, dN diff < +7, both laws) AND M4
      quiet (< 15% both absorbers both laws).
  else AMBIGUOUS (reported as-is).
PAPER-1 RULE (pre-committed): SEPARATION-ABSENT -> the 2D
methodological claim does not operate on this data at the landed
configuration; Paper 1 is reframed around the model-light channels
(census, coherence, completeness measurement); alpha rows gain a
gamma-inert caveat. CONFIRMED -> the claim operates; Paper 1 proceeds
with the 7J-z verdict as its binary bottom line. ABSORBER-LEVEL ->
operates at nuisance grade; the alpha bottom line is quoted as
channel-independent.
CREDENCE MAP (releases the 7I/7J freeze; keyed on the 7J-z part-2
verdict x this tier; pre-committed):
  BOOST-REVIVES  + alpha/absorber grade -> anomaly-real ~55-60%
  BOOST-REVIVES  + ABSENT              -> ~50-55% (method demoted)
  NO-DETECTION   + alpha/absorber grade -> ~35-40% (method works; the
                                          binaries decline the boost)
  NO-DETECTION   + ABSENT              -> ~30-35% + Paper-1 reframe
  AMBIGUOUS anywhere                   -> hold ~45%, carried
AMENDMENT 8 (2026-07-26, review-prompted, logged BEFORE any vt cube
existed to read): the width channel sq is gamma-blind BY PHYSICS
(normalization errors do not rotate velocities) and by construction —
so including it in the absorber set dilutes the gamma-carried share of
total information and pulls the M-metrics toward "no separation" for a
design reason as well as a physical one. BOTH configurations are
therefore read from the SAME cubes: (i) sq FREE (the pre-registered
primary — separability of the full absorber set as fitted), and
(ii) sq PINNED to 0 (the sq=0 slice; the original FOUR-absorber
diagnosis of review rounds 4-6). The primary verdict tier comes from
(i) as pre-registered; (ii) is the interpretation guard — if the two
tiers differ, both are reported and the difference is itself the
finding (it localizes how much gamma-nullness the width channel
contributes by dilution).
GATES: both cubes finite everywhere; alpha=0 cross-law agreement
reported for both channels (determinism grade, diagnostic — the
Newton row is law-blind physics); the vt construction is a collapse of
the same model (disclosed, amendment 7c), so M-metrics measure the
gamma CHANNEL, not implementation differences.
Also reported (design check, not a finding): the wr posterior under vt
is expected to flatten — the near-parabolic population is
gamma-defined.
Output: data/stage7jg_read.txt
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
SQ = np.array([0.0, 0.1, 0.2, 0.3])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
SEEDS = (31, 101)

LANDED_OK = os.path.exists('data/stage7jz_prior.npz')
lnp = np.full(len(FCOMP), -1e9)
if LANDED_OK:
    pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
    fg, lp = pz['fh_grid'], pz['lnpi_host']
    inr = (FCOMP >= fg.min()) & (FCOMP <= fg.max())
    lnp[inr] = np.interp(FCOMP[inr], fg, lp)
    OPER = 'LANDED'
else:
    lnp = -0.5*((FCOMP-0.16)/0.08)**2
    OPER = 'LIT16 (fallback)'
P(f"7J-g operative anchor: {OPER}")

AF = np.arange(0.0, 2.0001, 0.01)
def read(cb9):
    cbp = cb9 + lnp.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp - m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=tuple(range(1, 9))), 1e-300)) + m0
    ima = int(np.argmax(lm)); am = A_GRID[ima]
    if 0 < ima < 4:
        x = A_GRID[ima-1:ima+2]; y = lm[ima-1:ima+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0: am = -c1/(2*c2)
    lmf = np.interp(AF, A_GRID, lm)
    Wa = float(0.01*np.sum(lmf >= lmf.max() - 0.5))
    sds, cors, posts = {}, {}, {}
    for name, ax, grid in (('wr', 2, WR_GRID), ('fcomp', 3, FCOMP),
                           ('sq', 8, SQ), ('fpm', 6, FPM)):
        m = ex.sum(axis=tuple(i for i in range(9) if i != ax))
        m = m/max(m.sum(), 1e-300)
        mu = float((m*grid).sum())
        sds[name] = float(np.sqrt(max(((grid-mu)**2*m).sum(), 0)))
        posts[name] = m
        j = ex.sum(axis=tuple(i for i in range(9) if i not in (0, ax)))
        j = j/max(j.sum(), 1e-300)
        av = A_GRID[:, None]*np.ones_like(j)
        tv = grid[None, :]*np.ones_like(j)
        ma_, mt_ = float((j*av).sum()), float((j*tv).sum())
        cv = float((j*(av-ma_)*(tv-mt_)).sum())
        sd_ = np.sqrt(max(float((j*(av-ma_)**2).sum()), 1e-24)
                      * max(float((j*(tv-mt_)**2).sum()), 1e-24))
        cors[name] = cv/max(sd_, 1e-12)
    return am, float(lm.max()-lm[0]), Wa, sds, cors, posts

R = {}
CONFIGS = ('sqfree', 'sq0')
for law in ('simple', 'BE'):
    for seed in SEEDS:
        c2d = np.load(f'data/stage7j_cube_full_photow_{seed}_{law}.npy')
        cvt = np.load(f'data/stage7j_cubevt_full_photow_{seed}_{law}.npy')
        assert np.isfinite(c2d).all() and np.isfinite(cvt).all(), "NaN cube"
        for ch, cc in (('2D', c2d), ('vt', cvt)):
            for cfg in CONFIGS:
                cs = cc if cfg == 'sqfree' else cc[..., :1]
                cb9 = cs + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
                R[(law, seed, ch, cfg)] = read(cb9)
                am, dn, Wa, sds, cors, posts = R[(law, seed, ch, cfg)]
                P(f"[{law} {seed} {ch} {cfg}] a_marg={am:.2f} "
                  f"dN={dn:+.1f} W_a={Wa:.2f} | SD(wr)={sds['wr']:.3f} "
                  f"SD(fcomp)={sds['fcomp']:.3f} SD(sq)={sds['sq']:.3f} "
                  f"SD(fpm)={sds['fpm']:.3f} | "
                  f"corr(a,wr)={cors['wr']:+.2f} "
                  f"corr(a,fcomp)={cors['fcomp']:+.2f} "
                  f"corr(a,sq)={cors['sq']:+.2f}")
                P(f"   P(wr)={np.round(posts['wr'],2).tolist()} "
                  f"P(fcomp)={np.round(posts['fcomp'],2).tolist()} "
                  f"P(sq)={np.round(posts['sq'],2).tolist()}")

# alpha=0 cross-law agreement (diagnostic, determinism grade)
for tagf, nm in (('stage7j_cube', '2D'), ('stage7j_cubevt', 'vt')):
    d0 = []
    for seed in SEEDS:
        a = np.load(f'data/{tagf}_full_photow_{seed}_simple.npy')[0]
        b = np.load(f'data/{tagf}_full_photow_{seed}_BE.npy')[0]
        d0.append(float(np.nanmax(np.abs(a-b))))
    P(f"alpha=0 cross-law max|diff| ({nm}): "
      f"{['%.2e' % v for v in d0]} (diagnostic)")

P("")
tiers = {}
for cfg in CONFIGS:
    fired, quiet_a, m4_fired, m4_quiet = [], True, False, True
    for law in ('simple', 'BE'):
        a2 = np.mean([R[(law, s, '2D', cfg)][0] for s in SEEDS])
        av = np.mean([R[(law, s, 'vt', cfg)][0] for s in SEEDS])
        d2 = np.mean([R[(law, s, '2D', cfg)][1] for s in SEEDS])
        dv = np.mean([R[(law, s, 'vt', cfg)][1] for s in SEEDS])
        W2 = np.mean([R[(law, s, '2D', cfg)][2] for s in SEEDS])
        Wv = np.mean([R[(law, s, 'vt', cfg)][2] for s in SEEDS])
        shrink = 1 - W2/max(Wv, 1e-9)
        shift = abs(a2-av)
        dnd = d2-dv
        P(f"[{cfg}] {law} seed-mean: a_marg 2D {a2:.2f} vs vt {av:.2f} "
          f"(shift {shift:.2f}); dN 2D {d2:+.1f} vs vt {dv:+.1f} "
          f"(diff {dnd:+.1f}); W_a 2D {W2:.2f} vs vt {Wv:.2f} "
          f"(shrink {shrink:+.2f})")
        if shrink >= 0.30: fired.append(f"M1({law})")
        if dnd >= 15: fired.append(f"M2({law})")
        if shift >= 0.30: fired.append(f"M3({law})")
        if not (shrink < 0.15 and shift < 0.15 and dnd < 7):
            quiet_a = False
        for ab in ('wr', 'fcomp'):
            s2 = np.mean([R[(law, s, '2D', cfg)][3][ab] for s in SEEDS])
            sv = np.mean([R[(law, s, 'vt', cfg)][3][ab] for s in SEEDS])
            sh = 1 - s2/max(sv, 1e-9)
            P(f"  [{cfg}] M4 {law} SD({ab}): vt {sv:.3f} -> 2D {s2:.3f} "
              f"(shrink {sh:+.2f})")
            if sh >= 0.30: m4_fired = True
            if sh >= 0.15: m4_quiet = False
    if fired:
        tiers[cfg] = (f"SEPARATION-CONFIRMED (alpha grade; "
                      f"{', '.join(fired)})")
    elif m4_fired:
        tiers[cfg] = "ABSORBER-LEVEL SEPARATION (M4 only)"
    elif quiet_a and m4_quiet:
        tiers[cfg] = "SEPARATION-ABSENT"
    else:
        tiers[cfg] = "AMBIGUOUS (reported as-is)"
P(f"\n==> 7J-g VERDICT (PRIMARY, sq free): {tiers['sqfree']}")
P(f"==> 7J-g co-read (sq pinned 0, the four-absorber configuration): "
  f"{tiers['sq0']}")
if tiers['sqfree'] != tiers['sq0']:
    P("(tiers differ -> the difference is itself the finding: it "
      "localizes the width channel's dilution share, amendment 8)")
P("(Paper-1 rule and the credence map are pre-committed in the "
  "docstring; applied in NOTES at booking.)")
with open('data/stage7jg_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("saved: data/stage7jg_read.txt")
