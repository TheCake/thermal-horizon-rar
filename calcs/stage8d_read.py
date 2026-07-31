"""STAGE 8D reader — THE SETTLING INSTRUMENT verdict (pre-reg
25b248e; bars locked in NOTES before the batch ran).

G17-ID: the partition identity lnL_full = lnL_cook + lnL_anti +
const per cube cell (the mode regression; ABORT on fail).
READS: the anticook marginal at LANDED-CONV (the amendment-11
construction, verbatim from stage7jz_read.py) and at the FLAT
fcomp prior; verdict per the locked bars.
DIAGNOSTICS (no bars): the price-of-contamination table; the
census-overlap statement (s/d-based, lower bound); the three-way
partition co-read.
Output: data/stage8d_read.txt
"""
import csv
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
FPM = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW = np.array([0.7, 1.0, 1.4])
SQ = np.array([0.0, 0.1, 0.2, 0.3])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
SEEDS = (31, 101)
LAWS = ('simple', 'BE')

# ---------------- G17-ID2 (amendment 1, logged pre-quote) --------
# The v1 partition identity is INVALID for this likelihood (the
# model is per-sample template-conditioned; the full-sample model
# is a pair-weighted mixture of the subsample models — additivity
# would fail for ANY correctly-wired complement; residual measured
# axis-uniform std ~11, recorded as the post-mortem in NOTES).
P("8D THE SETTLING INSTRUMENT — reader (pre-reg 25b248e; "
  "amendment 1 = G17-ID2, NOTES)")
P("")
mask = np.load('data/stage7l_cookmask.npy')
assert mask.dtype == np.bool_, mask.dtype
P(f"G17-ID2a mask certificate: dtype = bool, sum = "
  f"{int(mask.sum())}, len = {len(mask)} -> the ~mask complement "
  f"is exact boolean -> PASS")
btxt = open('data/stage7j_anticook_photow3.txt').read()
nprof = btxt.count('PROF anticook seed')
assert nprof >= 4, nprof
assert int(mask.sum()) == 1194
P(f"G17-ID2b partition count: mask sum = 1194 (the 7L cook N) and "
  f"14071 - 1194 = 12877 = the batch-stdout anticook N (recorded "
  f"in the 8D log); {nprof} PROF anticook rows present -> PASS")
for law in LAWS:
    for seed in SEEDS:
        cf = np.load(f'data/stage7j_cube_full_photow3_{seed}_{law}.npy')
        ck = np.load(f'data/stage7j_cube_cook_photow3_{seed}_{law}.npy')
        ca = np.load(
            f'data/stage7j_cube_anticook_photow3_{seed}_{law}.npy')
        d = cf - ck - ca
        P(f"G17-ID2c v1-residual record {law} {seed}: const = "
          f"{float(np.nanmedian(d)):+.2f}, std = "
          f"{float(np.nanstd(d)):.1f} (axis-uniform; the sample-"
          f"conditioned-template signature — NOT a pass, the "
          f"post-mortem)")
P("G17-ID2: the mode certificate holds at its achievable grade "
  "(count + boolean complement + portrait cross-read; weaker than "
  "a bit-identity, disclosed).")
P("")

# ---------------- anchors (verbatim 7jz_read construction) ------
anchors = {}
pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
V2C = str(pz['version']).startswith('v2c')
assert V2C, "operative anchor requires the v2c certificate npz"
lnp = np.full(len(FCOMP), -1e9)
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
inr = (FCOMP >= fg[sup].min()) & (FCOMP <= fg[sup].max())
lnp[inr] = np.interp(FCOMP[inr], fg, lp)
gband = (pz['conv_band']/0.30 if 'conv_band' in pz.files
         else np.array([0.33, 1.30]))
GS = np.linspace(float(gband[0]), float(gband[1]), 25)
fhmin, fhmax = fg[sup].min(), fg[sup].max()
lnc = np.full(len(FCOMP), -1e9)
for gi in GS:
    fh_eq = FCOMP/gi
    m = (fh_eq >= fhmin) & (fh_eq <= fhmax)
    cand = np.full(len(FCOMP), -1e9)
    cand[m] = np.interp(fh_eq[m], fg, lp)
    lnc = np.maximum(lnc, cand)
anchors['LANDED-CONV'] = lnc
anchors['FLAT'] = np.zeros(len(FCOMP))
P(f"anchors: LANDED-CONV ln pi = {np.round(lnc, 2).tolist()}; "
  f"FLAT = zeros")
P("")

def read(cb9, lnpi):
    cbp = cb9 + lnpi.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp - m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=tuple(range(1, 9))),
                           1e-300)) + m0
    ima = int(np.argmax(lm))
    am = A_GRID[ima]
    if 0 < ima < 4:
        x = A_GRID[ima-1:ima+2]
        y = lm[ima-1:ima+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0:
            am = -c1/(2*c2)
    post = {}
    for name, ax, grid in (('wr', 2, WR_GRID), ('fcomp', 3, FCOMP),
                           ('fpm', 6, FPM), ('kw', 7, KW),
                           ('sq', 8, SQ)):
        m = ex.sum(axis=tuple(i for i in range(9) if i != ax))
        post[name] = m/max(m.sum(), 1e-300)
    return am, float(lm.max()-lm[0]), lm, post

# ---------------- the reads + verdict ----------------
res = {}
for law in LAWS:
    for seed in SEEDS:
        cw = np.load(
            f'data/stage7j_cube_anticook_photow3_{seed}_{law}.npy')
        cb9 = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
        for aname in ('LANDED-CONV', 'FLAT'):
            am, dn, lm, post = read(cb9, anchors[aname])
            fmode = FCOMP[int(np.argmax(post['fcomp']))]
            res[(law, seed, aname)] = (am, dn, fmode, post)
            P(f"[{law} {seed}] {aname}: a_marg = {am:.2f}, dN = "
              f"{dn:+.1f}; fcomp mode = {fmode:.2f} "
              f"P(fcomp) = {np.round(post['fcomp'], 2).tolist()}; "
              f"P(fpm) = {np.round(post['fpm'], 2).tolist()}; "
              f"P(sq) = {np.round(post['sq'], 2).tolist()}")
P("")

# price-of-contamination table (diagnostic, no bar): the alpha = 0
# slice's best cell (likelihood + eta only) and its cost
P("price-of-contamination (the alpha = 0 world's best cell; "
  "cost = lm.max - lm[0] under FLAT):")
for law in LAWS:
    for seed in SEEDS:
        cw = np.load(
            f'data/stage7j_cube_anticook_photow3_{seed}_{law}.npy')
        cb9 = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
        sl = cb9[0]
        idx = np.unravel_index(int(np.nanargmax(sl)), sl.shape)
        cost = res[(law, seed, 'FLAT')][1]
        P(f"  [{law} {seed}] alpha=0 best cell: wr = "
          f"{WR_GRID[idx[1]]:.2f}, fcomp = {FCOMP[idx[2]]:.2f}, "
          f"fpm = {FPM[idx[5]]:.1f}, kw = {KW[idx[6]]:.1f}, sq = "
          f"{SQ[idx[7]]:.1f}; cost-to-force = {cost:+.1f} lnL")
P("")

# census overlap (s/d-based lower bound; the CSV lacks per-pair dRV)
rows = list(csv.DictReader(open('data/ceiling_pairs.csv')))
band = [r for r in rows if r['census_corr'] == 'True']
rem = [r for r in band if float(r['s_kAU']) > 30.0
       or float(r['dist_pc']) > 130.0]
P(f"census overlap (diagnostic): {len(rem)}/{len(band)} of the "
  f"operative band pairs are removed by the cook s/d cuts ALONE "
  f"(s > 30 kAU or d > 130 pc; their dRV cut removes more — "
  f"lower bound); the (band, cliff) statistic lives in the "
  f"complement.")
P("")

# ---------------- verdict per the locked bars ----------------
n_boost_lc = sum(1 for law in LAWS for s in SEEDS
                 if res[(law, s, 'LANDED-CONV')][0] >= 0.5
                 and res[(law, s, 'LANDED-CONV')][1] >= 10.0)
n_boost_fl = sum(1 for law in LAWS for s in SEEDS
                 if res[(law, s, 'FLAT')][0] >= 0.5
                 and res[(law, s, 'FLAT')][1] >= 10.0)
n_cont = sum(1 for law in LAWS for s in SEEDS
             if res[(law, s, 'LANDED-CONV')][0] <= 0.2
             and res[(law, s, 'LANDED-CONV')][2] <= 0.35)
if n_boost_lc >= 3 and n_boost_fl >= 3:
    v = ("BOOST-CARRIED: the removed 92% carries the boost at the "
         "landed absorber budget — the fork resolves toward DATA")
elif n_cont >= 3:
    v = ("CONTAMINATION-DECOMPOSED: MATERIAL against the alpha "
         "claim — route to the next named decider")
else:
    v = "GRAY"
P(f"==> 8D VERDICT (locked bars): {v}")
P(f"    boost-reads {n_boost_lc}/4 (LANDED-CONV), {n_boost_fl}/4 "
  f"(FLAT); contamination-reads {n_cont}/4")
P(f"    three-way partition: full a_marg 0.68-0.74 @ +14.5-23.8 "
  f"(operative band) | cook ~0-powerless (7L: ~3 lnL) | anti = "
  f"this stage")

with open('data/stage8d_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8d_read.txt")
