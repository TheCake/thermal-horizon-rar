# -*- coding: utf-8 -*-
"""Stage 7J-c: post-verdict cube diagnostics (read-only; no refitting).

Runs AFTER the 7J-B photo verdicts fired.  Pre-registered pipeline rule
being applied (7I, commit a17d3a5): "FPM/nuisance grids are sample
properties - profile them; watch grid-edge riding."  The photo PROF rows
rode fpm to the NEW grid edge (1.8 -> 2.4 after the correction-#4
extension), so before the verdict is interpreted this reader asks, from
the cached cubes:

  (a) where the marginalized posterior sits (fcomp, fpm, kw);
  (b) whether COMPANION-WIN survives with fpm restricted to <= 1.8 and
      to the program-standard 1.5 (external anchor: Gaia EDR3 formal-
      error underestimation is bounded at ~1.1-1.4 for this regime,
      Lindegren+21; f_pm = 2.4 is outside anything published);
  (c) the alpha-marginal gap: how much the marginal disfavors alpha = 1
      relative to alpha = 0;
  (d) the conditional lnL(fcomp) rows in BOTH amplitude modes - the
      photo-mode analogue of the amendment-5 as-published conditional
      (raw: alpha = 0 at the measured host fraction at a cost of
      -505..-1078 lnL) - so both conditionals live on disk.

The fired verdicts are NOT touched; bars stand as written.
Output: data/stage7j_diag.txt
"""
import numpy as np

OUT = []
def P(s):
    print(s); OUT.append(s)

E_GRID = np.array([1.05, 1.3])
A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP = np.array([0.0, 0.1, 0.2, 0.35, 0.5, 0.7])
FPM_P = np.array([1.2, 1.5, 1.8, 2.1, 2.4])   # photo cubes
FPM_R = np.array([1.2, 1.5, 1.8])             # raw cubes
KW = np.array([0.7, 1.0, 1.4])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
# discrete-cell completeness priors as printed in the 7J-B headers
LNPI = {'full':   np.array([-1584.82, -768.17, -268.10, -7.78, -1.06, -42.57]),
        'strict': np.array([-1584.82, -608.39, -150.50, -2.45, -3.38, -68.11])}

def prof_a(cb):
    prof = np.nanmax(cb, axis=(1, 2, 3, 4, 5, 6, 7))
    imax = int(np.nanargmax(prof)); ahat = A_GRID[imax]
    if 0 < imax < len(A_GRID)-1:
        x = A_GRID[imax-1:imax+2]; y = prof[imax-1:imax+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0: ahat = -c1/(2*c2)
    return prof, ahat, imax

def marg_lm(cb, lnpi):
    cbp = cb + lnpi[None, None, None, :, None, None, None, None]
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp - m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=(1, 2, 3, 4, 5, 6, 7)), 1e-300)) + m0
    return lm, ex

for sample in ('full', 'strict'):
    lnpi = LNPI[sample]
    for law in ('simple', 'BE'):
        for seed in (31, 101):
            tag = f"{sample} {seed} {law}"
            # ---------------- photo cube ----------------
            cube = np.load(f'data/stage7j_cube_{sample}_photo_{seed}_{law}.npy')
            cb = cube + prior_eta[None, :, None, None, None, None, None, None]
            lm, ex = marg_lm(cb, lnpi)
            fpost = ex.sum(axis=(0, 1, 2, 4, 5, 6, 7)); fpost /= fpost.sum()
            ppost = ex.sum(axis=(0, 1, 2, 3, 4, 5, 7)); ppost /= ppost.sum()
            P(f"[photo {tag}] alpha-marginal lm-lm[0]="
              f"{np.round(lm-lm[0],1).tolist()}")
            P(f"[photo {tag}] P(fcomp)={np.round(fpost,2).tolist()} "
              f"P(fpm)={np.round(ppost,2).tolist()}")
            for lab, sl in (('fpm<=1.8', slice(0, 3)),
                            ('fpm=1.5', slice(1, 2))):
                cbs = cb[:, :, :, :, :, :, sl, :]
                prof, ah2, im2 = prof_a(cbs)
                dn2 = float(np.nanmax(prof) - prof[0])
                lm2, ex2 = marg_lm(cbs, lnpi)
                ima = int(np.argmax(lm2)); am2 = A_GRID[ima]
                dnm2 = float(lm2.max() - lm2[0])
                f2 = ex2.sum(axis=(0, 1, 2, 4, 5, 6, 7)); f2 /= f2.sum()
                k2 = ex2.sum(axis=(0, 1, 2, 3, 4, 5, 6)); k2 /= k2.sum()
                P(f"[photo {tag}] {lab}: PROF a_hat={ah2:.2f} "
                  f"(interior={0<im2<4}) dN={dn2:+.1f} | MARG a_marg={am2:.2f} "
                  f"dN_marg={dnm2:+.1f} P(fcomp)={np.round(f2,2).tolist()} "
                  f"P(kw)={np.round(k2,2).tolist()}")
            rows = np.nanmax(cb, axis=(0, 1, 2, 4, 5, 6, 7)); rows -= rows[1]
            a_at = [round(float(prof_a(cb[:, :, :, fi:fi+1])[1]), 2)
                    for fi in range(len(FCOMP))]
            P(f"[photo {tag}] cond: lnL(fcomp)-lnL(0.1)="
              f"{np.round(rows,1).tolist()}  a_prof(fcomp)={a_at}")
            # ---------------- raw cube (as-published amplitude law) ----
            rw = np.load(f'data/stage7j_cube_{sample}_{seed}_{law}.npy')
            if rw.ndim == 7:
                rw = rw[..., None]
            rb = rw + prior_eta[None, :, None, None, None, None, None, None]
            rrows = np.nanmax(rb, axis=(0, 1, 2, 4, 5, 6, 7)); rrows -= rrows[1]
            ra_at = [round(float(prof_a(rb[:, :, :, fi:fi+1])[1]), 2)
                     for fi in range(len(FCOMP))]
            P(f"[raw   {tag}] cond: lnL(fcomp)-lnL(0.1)="
              f"{np.round(rrows,1).tolist()}  a_prof(fcomp)={ra_at}")

P("")
P("READINGS (7J-c):")
P(" R1 the COMPANION-WIN marginal is NOT carried by the fpm grid edge:")
P("    restricting to fpm<=1.8 or fpm=1.5 leaves a_marg=0.00 dN_marg=+0.0")
P("    and the posterior cell (fcomp=0.35 full / 0.20 strict, kw=0.7)")
P("    unchanged in all 8 sample x law x seed reads.")
P(" R2 mechanism: the likelihood alone still prefers (alpha~1, fcomp~0.1)")
P("    but the measured-multiplicity prior prices that cell at -768/-608;")
P("    at the prior-allowed fcomp the fitted alpha is 0.  The posterior")
P("    chooses 'Newton + measured companions, kinematic misfit -60..-116'")
P("    over 'alpha~1 + forbidden multiplicity'.")
P(" R3 the winning cell is itself a poor absolute fit (-60..-116 vs the")
P("    unconstrained optimum): under this forward model NO cell fits both")
P("    the kinematics and the measured multiplicity.  The photocenter")
P("    correction softened the multiplicity cost ~5x (raw -505..-525 at")
P("    fcomp=0.35; photo -67..-116) without closing it.")
P(" R4 kw pins to 0.7 (lowest wobble amplitude) everywhere - soft edge,")
P("    second-order for the verdict (fpm slices leave it unchanged).")

with open('data/stage7j_diag.txt', 'w') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage7j_diag.txt")
