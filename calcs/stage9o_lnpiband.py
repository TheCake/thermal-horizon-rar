"""STAGE 9O — THE LNPI CONVERSION-BAND ROW (round-13 queue item).
Pre-registered BEFORE any run (measurement round; NO credence
movement).  PURE READER over the archived 9I/9J npz tables — do
the round-13 operative letters depend on the companion-fraction
prior's conversion treatment?

Variants: OPER (shipped envelope, bit-compared to stored npz),
G-LOW / G-MID / G-HIGH (conversion pinned at band edges / mean),
FLAT (uniform incl. fcomp = 0).  Letters: L1 "alpha >= 0.5
excluded on clean strata" <=> EX05 <= -8; L2 "no Newton verdict
either way" <=> -8 <= dN_fine <= +8.  Bars (ordered): O-FRAGILE /
O-ROBUST / O-GRAY-CARRIED.
Gates: G9O-0 LNPI bit-identity to stored npz (4/4 x 2 files);
G9O-1 analytic lse 1e-12; G9O-2 regression under OPER to 9J/9I
printed values.  Output: data/stage9o_lnpiband.txt
"""
import math, re
import numpy as np
from scipy.special import logsumexp

L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

P("9O THE LNPI CONVERSION-BAND ROW (pre-reg committed BEFORE any "
  "run; measurement round; NO credence movement)")
P("")

FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
A_FULL = np.concatenate([np.arange(0, 0.61, 0.1),
                         np.array([0.7, 0.8, 0.9, 1.0, 1.1, 1.2])])
SEEDS = (31, 101)
NC_STD = 6*4
NC_SHR = 2*3*2

pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
gband = pz['conv_band']/0.30
GS = np.linspace(float(gband[0]), float(gband[1]), 25)
fhmin, fhmax = fg[sup].min(), fg[sup].max()

def lnpi_at(gi_list):
    out = np.full(len(FCOMP_GRID), -1e9)
    for gi in gi_list:
        fh_eq = FCOMP_GRID/gi
        m = (fh_eq >= fhmin) & (fh_eq <= fhmax)
        cand = np.full(len(FCOMP_GRID), -1e9)
        cand[m] = np.interp(fh_eq[m], fg, lp)
        out = np.maximum(out, cand)
    return out

VARIANTS = [
    ('OPER',   lnpi_at(GS)),
    ('G-LOW',  lnpi_at([float(gband[0])])),
    ('G-MID',  lnpi_at([0.5*(float(gband[0])+float(gband[1]))])),
    ('G-HIGH', lnpi_at([float(gband[1])])),
    ('FLAT',   np.zeros(len(FCOMP_GRID))),
]
P(f"conversion band g = [{float(gband[0]):.3f}, "
  f"{float(gband[1]):.3f}]; host support fh = [{fhmin:.3f}, "
  f"{fhmax:.3f}]")
for nm, v in VARIANTS:
    P(f"LNPI[{nm:6s}] = "
      + "/".join(("-inf" if x < -1e8 else f"{x:+.2f}") for x in v))
P("")

def lse_cells(T, axes, ncells):
    return logsumexp(T, axis=axes) - math.log(ncells)
_T = np.log(np.array([[1.0, 2.0], [3.0, 4.0]]))
_v = float(lse_cells(_T, (0, 1), 4))
g1_ok = abs(_v - math.log(2.5)) < 1e-12
P(f"G9O-1 combiner analytic unit check: lse = {_v:.12f} vs "
  f"ln(2.5) = {math.log(2.5):.12f} -> {'PASS' if g1_ok else 'FAIL'}")

# printed regression targets
g9j = {}
for ln in open('data/stage9j_stdext.txt').read().splitlines():
    mm = re.match(r"\[(simple|BE) (\d+)\] M-STD  full: .*"
                  r"a_marg = ([\d.]+)", ln)
    if mm:
        g9j[(mm.group(1), int(mm.group(2)), 'STD')] = \
            float(mm.group(3))
    mm = re.match(r"\[(simple|BE) (\d+)\] M-DROP full: a_marg = "
                  r"([\d.]+), P\(a=0\) = ([\d.]+)", ln)
    if mm:
        g9j[(mm.group(1), int(mm.group(2)), 'DROP')] = \
            (float(mm.group(3)), float(mm.group(4)))
assert len(g9j) == 8, g9j
g9i = {}
for ln in open('data/stage9i_finealpha.txt').read().splitlines():
    mm = re.match(r"\[(simple|BE) (\d+)\] M-DROP lnZ\(a\)-lnZ\(0\)"
                  r" = (.+)$", ln)
    if mm:
        g9i[(mm.group(1), int(mm.group(2)), 'CURVE')] = \
            [float(x) for x in mm.group(3).split("/")]
    mm = re.match(r"\[(simple|BE) (\d+)\] M-DROP fine: .*"
                  r"dN = ([-+\d.]+)$", ln)
    if mm:
        g9i[(mm.group(1), int(mm.group(2)), 'DN')] = \
            float(mm.group(3))
assert len(g9i) == 8, g9i

def marg_config(STDF, LNPIv, strata):
    lnZ = np.zeros(13)
    for ai in range(13):
        S = LNPIv[:, None, None, None].copy()
        for qi in strata:
            S = S + lse_cells(STDF[ai, qi], (2, 4), NC_STD)
        lnZ[ai] = logsumexp(S) - math.log(NC_SHR)
    return lnZ

g0_ok = g2_ok = True
res = {}
for seed in SEEDS:
    for law in ('simple', 'BE'):
        z9i = np.load(f'data/stage9i_tables_{seed}_{law}.npz')
        z9j = np.load(f'data/stage9j_tables_{seed}_{law}.npz')
        ok0 = (np.array_equal(z9i['LNPI'], VARIANTS[0][1])
               and np.array_equal(z9j['LNPI'], VARIANTS[0][1]))
        g0_ok &= ok0
        P(f"[{law} {seed}] G9O-0 LNPI bit-identity (9I + 9J npz) "
          f"-> {'PASS' if ok0 else 'FAIL'}")
        STDF = np.concatenate([z9i['STD'], z9j['STDX']], axis=0)
        for nm, LNPIv in VARIANTS:
            zd = marg_config(STDF, LNPIv, [0, 1, 2])
            zs = marg_config(STDF, LNPIv, [0, 1, 2, 3])
            wd = np.exp(zd - logsumexp(zd))
            ws = np.exp(zs - logsumexp(zs))
            am_d = float(np.sum(wd*A_FULL))
            am_s = float(np.sum(ws*A_FULL))
            p0 = float(wd[0]); p01 = float(wd[0]+wd[1])
            p05 = float(np.sum(wd[5:]))
            ex05 = float(zd[5] - np.max(zd))
            z7 = zd[:7]
            dn = float((logsumexp(z7[1:]) - math.log(6)) - z7[0])
            res[(law, seed, nm)] = dict(am_d=am_d, am_s=am_s, p0=p0,
                                        p01=p01, p05=p05, ex05=ex05,
                                        dn=dn, zd=zd)
            P(f"[{law} {seed}] {nm:6s}: a_marg(DROP) = {am_d:.3f}, "
              f"P(a=0) = {p0:.3f}, P(a<=0.1) = {p01:.3f}, "
              f"P(a>=0.5) = {p05:.4f}, EX05 = {ex05:+.1f}, "
              f"dN_fine = {dn:+.1f}; a_marg(STD) = {am_s:.3f}")
        # G9O-2 regression under OPER
        r = res[(law, seed, 'OPER')]
        tj_s = g9j[(law, seed, 'STD')]
        tj_d, tj_p0 = g9j[(law, seed, 'DROP')]
        curve = g9i[(law, seed, 'CURVE')]
        tdn = g9i[(law, seed, 'DN')]
        c5 = float(r['zd'][5] - r['zd'][0])
        ok2 = (abs(r['am_s'] - tj_s) <= 0.002
               and abs(r['am_d'] - tj_d) <= 0.002
               and abs(r['p0'] - tj_p0) <= 0.002
               and abs(c5 - curve[5]) <= 0.06
               and abs(r['dn'] - tdn) <= 0.06)
        g2_ok &= ok2
        P(f"[{law} {seed}] G9O-2 regression: STD {r['am_s']:.3f} vs "
          f"{tj_s:.3f}; DROP {r['am_d']:.3f} vs {tj_d:.3f}; P0 "
          f"{r['p0']:.3f} vs {tj_p0:.3f}; curve(0.5) {c5:+.1f} vs "
          f"{curve[5]:+.1f}; dN {r['dn']:+.1f} vs {tdn:+.1f} -> "
          f"{'PASS' if ok2 else 'FAIL'}")
        P("")

if not (g0_ok and g1_ok and g2_ok):
    P("GATES FAILED (G9O-0 " + ('PASS' if g0_ok else 'FAIL')
      + ", G9O-1 " + ('PASS' if g1_ok else 'FAIL')
      + ", G9O-2 " + ('PASS' if g2_ok else 'FAIL')
      + ") - STOP; DO NOT QUOTE; verdict WITHHELD")
else:
    P("GATES: G9O-0 4/4, G9O-1, G9O-2 4/4 - ALL PASS")
    P("")
    fragile = False
    all_hold = True
    for nm, _ in VARIANTS:
        v1 = sum(1 for s in SEEDS for lw in ('simple', 'BE')
                 if res[(lw, s, nm)]['ex05'] > -8.0)
        v2 = sum(1 for s in SEEDS for lw in ('simple', 'BE')
                 if not (-8.0 <= res[(lw, s, nm)]['dn'] <= 8.0))
        P(f"O-bars [{nm:6s}]: L1-broken rows {v1}/4; L2-broken "
          f"rows {v2}/4")
        if v1 >= 2 or v2 >= 2: fragile = True
        if v1 > 0 or v2 > 0: all_hold = False
    ams = [res[(lw, s, nm)]['am_d'] for s in SEEDS
           for lw in ('simple', 'BE') for nm, _ in VARIANTS]
    P(f"quotable band: a_marg(DROP) across variants = "
      f"[{min(ams):.3f}, {max(ams):.3f}]")
    P("")
    if fragile:
        P("==> 9O VERDICT (locked grammar): O-FRAGILE - a prior "
          "variant flips an operative letter; the alpha band gains "
          "an explicit prior-conditional column and the paper "
          "quotes the spread.")
    elif all_hold:
        P("==> 9O VERDICT (locked grammar): O-ROBUST - both "
          "operative letters hold under every conversion-band and "
          "flat-prior treatment; the upper-limit language carries "
          "a prior-robustness pointer.")
    else:
        P("==> 9O VERDICT (locked grammar): O-GRAY-CARRIED - rows "
          "stand as measurements.")
    P("    NO credence movement (pre-stated; measurement round).")

with open('data/stage9o_lnpiband.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage9o_lnpiband.txt")
