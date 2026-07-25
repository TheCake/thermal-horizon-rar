# -*- coding: utf-8 -*-
"""Stage 7J-w: arm-recovery diagnosis (review round 5; measurements
only, no verdict bars - the cadence rule is in force).

The round-5 question: the fullpow simple arm recovered 1.48 against a
1.18 truth (+0.3, a 25% bias on the headline parameter) - is it
prior-induced (should move at the literature anchor) or in the
likelihood (propagates into every alpha this machinery quotes)?

The cached PROF rows already answer the dichotomy's third way: the
bias is PRIOR-FREE (PROF 1.47 vs MARG 1.48 under the measured prior;
fullpowbe: PROF = MARG = 0.73).  Both arms' own-truth anomalies
(simple +0.29, BE -0.40) sit in the likelihood/realization layer, at
single-injection grade where the Stage-3A realization systematic
covers part of the amplitude; the multi-truth-seed recovery map is
queued with the post-7J-z arm suite (reviewer ordering: 7J-z -> 7J-g
-> arms at the landed anchor; the BE alpha=0.4 additive-vs-
multiplicative discriminator runs there too).

This script adds the number the cubes give for free: the SAME injected
data re-marginalized under the literature prior (centres 0.16 / 0.22).
The injected truth carries companions at the (retracted-measurement)
0.51 rate, so the literature prior is deliberately MISMATCHED to the
truth - the read measures the direction and size of prior-
misspecification bias on alpha: if a prior that understates the true
companion rate inflates alpha-hat, then the real-data interior 0.4 at
the literature anchor is an UPPER bound under the branch where the sky
truly hosts ~0.5.  Cube caveat made explicit (round-5 point): cubes
are prior-independent, so re-marginalizing is free on the SEED axis -
but the power ARMS at a new anchor need new injected data (different
companion population), which is why no arm is validated at the
operative prior yet.

Output: data/stage7j_armdiag.txt
"""
import numpy as np

OUT = []
def P(s):
    print(s); OUT.append(s)

E_GRID = np.array([1.05, 1.3])
A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
FCOMP = np.array([0.0, 0.1, 0.2, 0.35, 0.5, 0.7])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
LNPI_MEAS = np.array([-1584.82, -768.17, -268.10, -7.78, -1.06, -42.57])

def amax(lm):
    ima = int(np.argmax(lm)); am = A_GRID[ima]
    if 0 < ima < 4:
        x = A_GRID[ima-1:ima+2]; y = lm[ima-1:ima+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0: am = -c1/(2*c2)
    return am, float(lm.max() - lm[0])

def marg(cb, lnpi):
    cbp = cb + lnpi[None, None, None, :, None, None, None, None]
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp - m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=(1, 2, 3, 4, 5, 6, 7)), 1e-300)) + m0
    fp = ex.sum(axis=(0, 1, 2, 4, 5, 6, 7)); fp /= fp.sum()
    return amax(lm), fp

for run, truth in (('fullpow', 'simple 1.18'), ('fullpowbe', 'BE 1.13')):
    for law in ('simple', 'BE'):
        cube = np.load(f'data/stage7j_cube_{run}_photo_31_{law}.npy')
        cb = cube + prior_eta[None, :, None, None, None, None, None, None]
        prof = np.nanmax(cb, axis=(1, 2, 3, 4, 5, 6, 7))
        ap, dp = amax(prof)
        (am_m, dm_m), _ = marg(cb, LNPI_MEAS)
        rows = [f"PROF={ap:.2f}/{dp:+.1f}", f"meas={am_m:.2f}/{dm_m:+.1f}"]
        for cen in (0.22, 0.16):
            lit = -0.5*((FCOMP - cen)/0.08)**2
            (am, dm), fp = marg(cb, lit)
            rows.append(f"lit{cen:.2f}={am:.2f}/{dm:+.1f} "
                        f"fpost@{FCOMP[int(np.argmax(fp))]:.2f}")
        P(f"[{run} truth={truth} arm={law}] " + "  ".join(rows))

P("")
P("READINGS (measurements, per the cadence rule):")
P(" R1 the +0.3/-0.40 own-truth anomalies are PRIOR-FREE (PROF equals")
P("    the matched-prior MARG on both arms): likelihood/realization")
P("    layer, single-injection grade; recovery map queued post-7J-z.")
P(" R2 the lit-prior rows read the prior-MISMATCH direction: truth")
P("    companions at 0.51 marginalized under a ~0.16-0.22 prior - if")
P("    alpha-hat inflates, an understated companion prior biases alpha")
P("    HIGH, making the real-data interior ~0.4 at the literature")
P("    anchor an upper bound under the high-multiplicity branch.")
P(" R3 explicit scope note: cube reuse is free on the SEED axis (prior")
P("    enters at marginalization); power ARMS at a new anchor require")
P("    NEW injected data - no arm is validated at the operative prior;")
P("    fullpowlit runs at whatever anchor 7J-z lands (reviewer order).")

# --- 7J-w2: injection informativeness vs the sky (round 6) ----------------
# The reviewer's asymmetry: on injected data the likelihood pins
# fcomp = 0.50 through a mismatched prior; on real data alpha moves as
# the anchor moves - the injected likelihood is MORE informative about
# fcomp than the sky's, so power validated on injections overstates
# power on data.  Quantify from cached cubes at matched N: the
# fcomp lnL profile (max over all other axes), real vs injected.
P("")
P("7J-w2: fcomp informativeness, real vs injected (full sample, lnL")
P("profile relative to its own maximum; more negative = sharper)")
for tag, path in (('real 31', 'data/stage7j_cube_full_photo_31_simple.npy'),
                  ('real 101', 'data/stage7j_cube_full_photo_101_simple.npy'),
                  ('inj-fullpow', 'data/stage7j_cube_fullpow_photo_31_simple.npy'),
                  ('inj-fullpowbe(BE arm)',
                   'data/stage7j_cube_fullpowbe_photo_31_BE.npy')):
    c = np.load(path)
    cb = c + prior_eta[None, :, None, None, None, None, None, None]
    prof = np.nanmax(cb, axis=(0, 1, 2, 4, 5, 6, 7))
    prof = prof - prof.max()
    P(f"  [{tag}] lnL(fcomp)-max = "
      f"{np.round(prof, 1).tolist()}  (cells {FCOMP.tolist()})")
P(" reading: the real profiles are shallow across 0.1-0.35 (tens of")
P(" lnL) where the injected profiles cliff by hundreds around their")
P(" truth cell - if so, the sky's companion signature is weaker than")
P(" the model's own, and every arm validation inherits that optimism;")
P(" quantified here so the post-7J-z arm suite states it.")

with open('data/stage7j_armdiag.txt', 'w') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage7j_armdiag.txt")
