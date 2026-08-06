"""STAGE 9K — THE FPM-CEILING CURVE (round-13 D1; reader-only).
Pre-registered BEFORE the run (measurement round; NO credence
movement).

From data/stage9i_tables_*.npz: re-marginalize M-DROP fine with
the fpm axis truncated at grid-prefix caps {1.2, 1.5, 1.8, 2.4,
3.0}.  Reads: alpha_marg(cap), P(alpha=0)(cap), max-mode peak
lnL(cap); M-STD capped co-read; per-stratum-alone DlnZ(<=1.5 vs
full).  Gates: G9K-0 full-cap lineage to 9I printed (8 values,
0.002); G9K-1 analytic (1e-12); G9K-2 cap-monotone peak lnL.
Bars (ordered): K-ROBUST / K-CAP-FRAGILE / K-GRAY-CARRIED.
Output: data/stage9k_fpmcap.txt
"""
import math, re
import numpy as np
from scipy.special import logsumexp

L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

A_FINE = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
NA = 7
FPM_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
SQ_N = 4
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
SEEDS = (31, 101)
NC_SHR = 2*3*2
CAPS = [(1.2, 1), (1.5, 2), (1.8, 3), (2.4, 5), (3.0, 6)]

pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
gband = pz['conv_band']/0.30
GS = np.linspace(float(gband[0]), float(gband[1]), 25)
fhmin, fhmax = fg[sup].min(), fg[sup].max()
LNPI = np.full(len(FCOMP_GRID), -1e9)
for gi in GS:
    fh_eq = FCOMP_GRID/gi
    m = (fh_eq >= fhmin) & (fh_eq <= fhmax)
    cand = np.full(len(FCOMP_GRID), -1e9)
    cand[m] = np.interp(fh_eq[m], fg, lp)
    LNPI = np.maximum(LNPI, cand)

g9i = {}
for ln in open('data/stage9i_finealpha.txt').read().splitlines():
    mm = re.match(r"\[(simple|BE) (\d+)\] M-(DROP|STD) {0,2}fine: .*"
                  r"a_marg = ([\d.]+)", ln)
    if mm:
        g9i[(mm.group(1), int(mm.group(2)), mm.group(3))] = \
            float(mm.group(4))
assert len(g9i) == 8, g9i

P("9K THE FPM-CEILING CURVE (round-13 D1; pre-reg committed BEFORE "
  "the run; reader-only; measurement round; NO credence movement)")
P("")

def lse_cells(T, axes, ncells):
    return logsumexp(T, axis=axes) - math.log(ncells)
_T = np.log(np.array([[1.0, 2.0], [3.0, 4.0]]))
_v = float(lse_cells(_T, (0, 1), 4))
g1_ok = abs(_v - math.log(2.5)) < 1e-12
P(f"G9K-1 analytic: {_v:.12f} vs ln(2.5) -> "
  f"{'PASS' if g1_ok else 'FAIL'}")
P("")

def marg_cap(STDz, strata, ktop):
    lnZ = np.zeros(NA)
    for ai in range(NA):
        S = LNPI[:, None, None, None].copy()
        for qi in strata:
            S = S + lse_cells(STDz[ai, qi, :, :, :ktop], (2, 4),
                              ktop*SQ_N)
        lnZ[ai] = logsumexp(S) - math.log(NC_SHR)
    return lnZ

def peak_cap(STDz, strata, ktop):
    best = -np.inf
    for ai in range(NA):
        S = LNPI[:, None, None, None].copy()
        for qi in strata:
            S = S + STDz[ai, qi, :, :, :ktop].max(axis=(2, 4))
        best = max(best, float(np.max(S)))
    return best

def am_of(lnZ):
    w = np.exp(lnZ - logsumexp(lnZ))
    return float(np.sum(w*A_FINE)), float(w[0])

g0_ok = g2_ok = True
rows = []
for seed in SEEDS:
    for law in ('simple', 'BE'):
        z = np.load(f'data/stage9i_tables_{seed}_{law}.npz')
        STD = z['STD']
        am_full_d, _ = am_of(marg_cap(STD, [0, 1, 2], 6))
        am_full_s, _ = am_of(marg_cap(STD, [0, 1, 2, 3], 6))
        ok0 = (abs(am_full_d - g9i[(law, seed, 'DROP')]) <= 0.002
               and abs(am_full_s - g9i[(law, seed, 'STD')]) <= 0.002)
        g0_ok &= ok0
        P(f"[{law} {seed}] G9K-0 full-cap lineage: DROP "
          f"{am_full_d:.3f} vs 9I {g9i[(law, seed, 'DROP')]:.3f}, "
          f"STD {am_full_s:.3f} vs {g9i[(law, seed, 'STD')]:.3f} "
          f"-> {'PASS' if ok0 else 'FAIL'}")
        pks = []
        for cap, k in CAPS:
            lnZd = marg_cap(STD, [0, 1, 2], k)
            amd, p0d = am_of(lnZd)
            lnZs = marg_cap(STD, [0, 1, 2, 3], k)
            ams, _ = am_of(lnZs)
            pk = peak_cap(STD, [0, 1, 2], k)
            pks.append(pk)
            P(f"[{law} {seed}] cap {cap:.1f}: DROP a_marg = "
              f"{amd:.3f}, P(a=0) = {p0d:.3f}, peak lnL = "
              f"{pk:+.1f}; STD a_marg = {ams:.3f}")
            if abs(cap - 1.5) < 1e-9:
                am15 = amd
            if k == 6:
                am30 = amd
        mono = all(pks[i] <= pks[i+1] + 1e-9
                   for i in range(len(pks)-1))
        g2_ok &= mono
        P(f"[{law} {seed}] G9K-2 cap-monotone peak lnL -> "
          f"{'PASS' if mono else 'FAIL'}")
        for qi in range(4):
            z15 = np.zeros(NA); zfu = np.zeros(NA)
            for ai in range(NA):
                z15[ai] = float(logsumexp(
                    STD[ai, qi, :, :, :2]
                    + LNPI.reshape(6, 1, 1, 1, 1, 1))
                    - math.log(2*SQ_N*NC_SHR))
                zfu[ai] = float(logsumexp(
                    STD[ai, qi]
                    + LNPI.reshape(6, 1, 1, 1, 1, 1))
                    - math.log(6*SQ_N*NC_SHR))
            d = float(logsumexp(z15) - math.log(NA)
                      - (logsumexp(zfu) - math.log(NA)))
            P(f"    Q{qi+1}-alone DlnZ(cap 1.5 - full) = {d:+.1f}")
        rows.append(dict(law=law, seed=seed, d=am15 - am30))
        P("")

if not (g0_ok and g1_ok and g2_ok):
    P("GATES FAILED (G9K-0 " + ('PASS' if g0_ok else 'FAIL')
      + ", G9K-1 " + ('PASS' if g1_ok else 'FAIL')
      + ", G9K-2 " + ('PASS' if g2_ok else 'FAIL')
      + ") - STOP; DO NOT QUOTE; verdict WITHHELD")
else:
    P("GATES: G9K-0 4/4, G9K-1, G9K-2 4/4 - ALL PASS")
    P("")
    n1 = sum(1 for r in rows if abs(r['d']) < 0.05)
    n2 = sum(1 for r in rows if r['d'] >= 0.15)
    P(f"K-bars: ROBUST-rows {n1}/4 (|d| < 0.05 at cap 1.5); "
      f"CAP-FRAGILE-rows {n2}/4 (d >= +0.15)")
    if n1 >= 3:
        P("==> 9K VERDICT (locked grammar): K-ROBUST - the "
          "drop-world small-alpha statement does not depend on "
          "fpm freedom above the Lindegren-adjacent cap; the "
          "deflation is not noise-absorption at the cap grade.")
    elif n2 >= 3:
        P("==> 9K VERDICT (locked grammar): K-CAP-FRAGILE - "
          "capping fpm at 1.5 raises alpha materially: the "
          "deflated alpha was partly noise-absorption "
          "(direction-1 pressure; the boost was being eaten); "
          "the D2 narrow-fpm meter decides legitimacy.")
    else:
        P("==> 9K VERDICT (locked grammar): K-GRAY-CARRIED - "
          "rows stand as measurements.")
    P("    NO credence movement (pre-stated; measurement round).")

with open('data/stage9k_fpmcap.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage9k_fpmcap.txt")
