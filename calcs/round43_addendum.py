"""ROUND 43 addendum -- the correction-arc verification.

GA half (BLIND, written and committed BEFORE the ROUND 43 report arrives;
the 87a4676 protocol's 12th execution): every load-bearing 9L-b number
re-derived by an INDEPENDENT method (plain probability sums, no logsumexp;
fresh marginalization code; scipy-free statistics).
GB half (post-report): re-compute every number the reviewer's report
quotes; appended after adoption.

GA checks:
  GA-1 flat-prior mean of FPM_GRID = 2.00 exactly (the nullinj diagnosis).
  GA-2 calibration means/SDs from data/stage9lb_calib.npz raw allS arrays
       (independent of the stage's own printed summary): reproduce the
       printed per-truth S means/SDs, D, sigma_pool, and the GB9LB-3
       verdicts for both seeds.
  GA-3 sky S(Q1..Q4) + pooled from data/stage9l_tables_{seed}.npz with
       FRESH marginalization code (max-subtraction exp sums, not scipy),
       vs the stage prints (+0.534/+0.481/+0.472/+0.341 pooled +2.204;
       +0.372/+0.640/+0.670/+0.629 pooled +0.921).
  GA-4 archived E[fpm] regression from the same fresh code
       (2.21/2.17/2.26/2.12; 2.13/2.27/2.25/2.30).
  GA-5 monotonicity re-check: adjacent steps vs -1 pooled SE, both seeds
       (PASS seed 31 / FAIL seed 101 expected).
  GA-6 letter selection: apply the pre-signed grammar mechanically to the
       GA-recomputed gate outcomes; must yield B-POWER-DEAD.
Output: data/round43_addendum.txt (GA block first; GB appended later).
"""
import io
import numpy as np

OUT = []
def P(s):
    print(s, flush=True)
    OUT.append(s)

FPM_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
ok_all = True

def check(name, got, want, tol):
    global ok_all
    ok = abs(got - want) <= tol
    ok_all &= ok
    P(f"  {name}: {got:+.4f} vs {want:+.4f} (tol {tol}) -> "
      f"{'OK' if ok else 'MISMATCH'}")
    return ok

P("ROUND 43 ADDENDUM -- GA half (BLIND, pre-report)")
P("")
P("GA-1 flat-prior mean of FPM_GRID:")
check("mean", float(FPM_GRID.mean()), 2.00, 1e-12)
P("")

def marg_fpm_fresh(T, LNPI):
    # fresh code: max-subtraction, plain exps, loops over kept axes
    A = T + LNPI.reshape(6, 1, 1, 1, 1, 1)
    m = A.max()
    W = np.exp(A - m)
    w2 = W.sum(axis=(0, 1, 3, 5))          # (fpm, sq)
    w2 = w2 / w2.sum()
    return w2.sum(axis=1)                   # fpm marginal

def S_fresh(T, LNPI):
    mf = marg_fpm_fresh(T, LNPI)
    return float(np.log(mf[3:].sum()) - np.log(mf[:3].sum()))

PRINTED = {
    31: dict(
        means=[-0.180, -0.224, -0.199, -0.116, -0.104, -0.025],
        sds=[0.432, 0.424, 0.439, 0.460, 0.369, 0.336],
        D=0.155, twosig=0.824,
        Sq=[+0.534, +0.481, +0.472, +0.341], Sp=+2.204,
        efpm=[2.21, 2.17, 2.26, 2.12], g2='PASS', g3='FAIL'),
    101: dict(
        means=[-0.212, -0.316, -0.217, +0.014, -0.005, -0.030],
        sds=[0.426, 0.466, 0.381, 0.350, 0.435, 0.393],
        D=0.182, twosig=0.821,
        Sq=[+0.372, +0.640, +0.670, +0.629], Sp=+0.921,
        efpm=[2.13, 2.27, 2.25, 2.30], g2='FAIL', g3='FAIL'),
}

cal = np.load('data/stage9lb_calib.npz')
TRUTHS = cal['truths']; R = int(cal['R'])
P(f"calibration archive: truths {list(TRUTHS)}, R = {R}")
P("")

letter_gates = {}
for seed in (31, 101):
    P(f"GA-2 [seed {seed}] calibration recomputation from raw allS:")
    allS = cal[f's{seed}_allS']            # (6 truths, R)
    means = allS.mean(axis=1)
    # independent SD: two-pass definition, explicit ddof=1 arithmetic
    sds = np.sqrt(((allS - means[:, None])**2).sum(axis=1) / (R - 1))
    for ti in range(len(TRUTHS)):
        check(f"mean(t={TRUTHS[ti]:.1f})", means[ti],
              PRINTED[seed]['means'][ti], 0.0006)
        check(f"sd(t={TRUTHS[ti]:.1f})", sds[ti],
              PRINTED[seed]['sds'][ti], 0.0006)
    D = means[-1] - means[0]
    sig_pool = float(np.sqrt(np.mean(sds**2)))
    check("D", D, PRINTED[seed]['D'], 0.0006)
    check("2*sigma_pool", 2*sig_pool, PRINTED[seed]['twosig'], 0.0006)
    g3 = 'PASS' if D >= 2*sig_pool else 'FAIL'
    okg = g3 == PRINTED[seed]['g3']
    ok_all &= okg
    P(f"  GB9LB-3 verdict: {g3} vs printed {PRINTED[seed]['g3']} -> "
      f"{'OK' if okg else 'MISMATCH'}")
    P("")
    P(f"GA-5 [seed {seed}] monotonicity:")
    se_pool = sig_pool/np.sqrt(R)
    inv = np.diff(means)
    g2 = 'PASS' if bool(np.all(inv > -se_pool)) else 'FAIL'
    okg = g2 == PRINTED[seed]['g2']
    ok_all &= okg
    P(f"  min step {inv.min():+.4f} vs -1 SE {-se_pool:.4f}: {g2} vs "
      f"printed {PRINTED[seed]['g2']} -> {'OK' if okg else 'MISMATCH'}")
    letter_gates[seed] = (g2, g3)
    P("")
    P(f"GA-3/GA-4 [seed {seed}] sky reads with fresh marginalization:")
    npz = np.load(f'data/stage9l_tables_{seed}.npz')
    T, LNPI = npz['T'], npz['LNPI']
    for qi in range(4):
        check(f"S(Q{qi+1})", S_fresh(T[qi], LNPI),
              PRINTED[seed]['Sq'][qi], 0.0006)
        check(f"E[fpm](Q{qi+1})",
              float(np.sum(marg_fpm_fresh(T[qi], LNPI)*FPM_GRID)),
              PRINTED[seed]['efpm'][qi], 0.005)
    check("S(pooled)", S_fresh(T.sum(axis=0), LNPI),
          PRINTED[seed]['Sp'], 0.0006)
    P("")

P("GA-6 mechanical letter selection from GA-recomputed gates:")
power_dead = any(g2 == 'FAIL' or g3 == 'FAIL'
                 for g2, g3 in letter_gates.values())
letter = 'B-POWER-DEAD' if power_dead else '(exclusion grammar)'
okL = letter == 'B-POWER-DEAD'
ok_all &= okL
P(f"  gates {letter_gates} => {letter} -> {'OK' if okL else 'MISMATCH'}")
P("")
P("GA VERDICT: " + ("ALL OK" if ok_all else "MISMATCH PRESENT -- STOP"))

with io.open('data/round43_addendum.txt', 'a', encoding='utf-8',
             newline='\n') as f:
    f.write("\n".join(OUT) + "\n")
print("\nappended: data/round43_addendum.txt")
