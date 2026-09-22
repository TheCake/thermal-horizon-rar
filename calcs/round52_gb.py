"""ROUND 52 verification -- GB half (post-report): re-compute every
load-bearing REVIEWER number in own implementations before adoption
(the standing verify-reviewer-math rule). Writes data/round52_gb.txt.

Legs:
  GB-1  exact noise-free world signatures incl. the reviewer-added
        WORLD-D calibration (his d_Ii/d_Oi/S1 table).
  GB-2  z arithmetic + Oi-axis mixture fractions.
  GB-3  both-sides-16 restricted d_Oi + its z row.
  GB-4  S1 construction grid: 7 spot cells + the 18-cell survive
        count + galaxy-weighted + cross-galaxy-only variants
        (parametric BRUTE matcher, own code).
  GB-5  shift-null SDs + the Ii floor, own stream (4000 draws).
  GB-6  pool-rebuilt galaxy-block bootstrap for S1, own seed (4000).
  GB-7  S1 carriers (ESO563-G021 jackknife-variance share, drop
        value) + d_Ii breadth (sign census, worst drop).
  GB-8  S1 y-truncation census + sub-20-arcsec R/Rd geometry.
  GB-9  signature specificity under WORLD-R / WORLD-D / WORLD-Y
        (own seeds, 400 draws each).
  GB-10 dispersion audit: point SDs; realized S1 SD under the
        registered and sky-matched-noise worlds; power at the
        sky's own fixed bars.
  GB-11 d_Ii vs the published C17 per-galaxy contrast (r).
  GB-12 P(S1 <= -2 SE | WORLD-R) at the -2 SE reading of the
        section-5 clause.
  GB-13 the single anchored overlap-match qualifier + value.
  (GB-14 the git-ordering audit was done directly in-session:
   prereg 21:38:45 -> A1 21:49:57 -> A2 22:02:16 -> A3 22:09:27 ->
   gates 22:15:59 -> sky 22:18:44 -> GA 22:19:46, matching the
   reviewer's table.)
"""
import math, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = open(os.path.join(HERE, 'stage10y_organizer.py'),
           encoding='utf-8').read()
cut = SRC.index('# ====')
NS = {'__name__': 'stage10y_trunc',
      '__file__': os.path.join(HERE, 'stage10y_organizer.py')}
_argv = sys.argv
sys.argv = ['stage10y_organizer.py', 'gates']
exec(compile(SRC[:cut], 'stage10y_trunc', 'exec'), NS)
sys.argv = _argv

GAL_A, GAL_F = NS['GAL_A'], NS['GAL_F']
QUAL_A = NS['QUAL_A']
TH, RDF = NS['TH_STAR'], NS['RD_FAC']
s1_matched, M_S1_A = NS['s1_matched'], NS['M_S1_A']
world_noise, world_y, world_r = (NS['world_noise'], NS['world_y'],
                                 NS['world_r'])
SIG_PT, RHO_AR = NS['SIG_PT'], NS['RHO_AR']

L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

t0 = time.time()
P("")
P("=" * 68)
P("ROUND 52 GB (post-report; reviewer numbers re-computed)")
P("")

RES_BY = [g_['res'] for g_ in GAL_A]

# ---------------- shared helpers ------------------------------------
def cells(gal, res_by):
    dii, doi = [], []
    for gd_, r in zip(gal, res_by):
        if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
            continue
        it = gd_['th'] < TH
        ir = gd_['rk'] < RDF*gd_['rd']
        Oo = ~it & ~ir
        if not Oo.any():
            continue
        base = float(r[Oo].mean())
        if (it & ir).any():
            dii.append(float(r[it & ir].mean()) - base)
        if (~it & ir).any():
            doi.append(float(r[~it & ir].mean()) - base)
    return np.array(dii), np.array(doi)

def pk(v):
    return (float(v.mean()), float(v.std(ddof=1)/math.sqrt(len(v))),
            len(v))

def brute(gal, qual, res_by, dy, minw, skip=None, weight='point',
          cross_only=False):
    A, Bly, Bres, Bj = [], [], [], []
    for j in qual:
        if skip is not None and j == skip:
            continue
        gd_ = gal[j]
        ly = np.log10(gd_['y'])
        for i in np.where(gd_['th'] < TH)[0]:
            A.append((j, float(res_by[j][i]), ly[i]))
        for i in np.where(gd_['th'] >= TH)[0]:
            Bly.append(ly[i]); Bres.append(float(res_by[j][i]))
            Bj.append(j)
    Bly = np.array(Bly); Bres = np.array(Bres); Bj = np.array(Bj)
    pts, gset = {}, set()
    for (j, ra, la) in A:
        k = np.abs(Bly - la) <= dy
        if cross_only:
            k &= (Bj != j)
        if k.sum() < minw:
            continue
        pts.setdefault(j, []).append(ra - float(Bres[k].mean()))
        gset.add(j)
    if not pts:
        return float('nan'), 0, 0
    if weight == 'point':
        allv = [v for vs in pts.values() for v in vs]
        return float(np.mean(allv)), len(allv), len(gset)
    gm = [float(np.mean(vs)) for vs in pts.values()]
    n = sum(len(vs) for vs in pts.values())
    return float(np.mean(gm)), n, len(gset)

def brute_jack(gal, qual, res_by, dy, minw, weight='point',
               cross_only=False):
    v0, m, ng = brute(gal, qual, res_by, dy, minw, weight=weight,
                      cross_only=cross_only)
    gals = [j for j in qual
            if any(True for i in np.where(gal[j]['th'] < TH)[0])]
    jk, tags = [], []
    for j in sorted(set(gals)):
        v, _, _ = brute(gal, qual, res_by, dy, minw, skip=j,
                        weight=weight, cross_only=cross_only)
        if np.isfinite(v):
            jk.append(v); tags.append(j)
    jk = np.array(jk)
    n = len(jk)
    se = math.sqrt((n - 1)/n*float(np.sum((jk - jk.mean())**2)))
    return v0, se, m, ng, jk, tags

# ---------------- GB-1 exact world signatures -----------------------
P("-- GB-1 exact noise-free world signatures --")
sigY = [NS['B_Y']*np.log10(g_['y']) for g_ in GAL_A]
sigR = [(-0.0857)*(g_['th'] < TH).astype(float) for g_ in GAL_A]
# WORLD-D calibration: step -delta inside 1.5 Rd, delta solved so the
# plain both-sides contrast = -0.0857
num = []
for j in QUAL_A:
    gd_ = GAL_A[j]
    it = gd_['th'] < TH
    ir = (gd_['rk'] < RDF*gd_['rd']) if np.isfinite(gd_['rd']) \
        and gd_['rd'] > 0 else np.zeros(len(gd_['th']), bool)
    num.append(float(ir[it].mean() - ir[~it].mean()))
DELTA_D = 0.0857/float(np.mean(num))
sigD = []
for g_ in GAL_A:
    ir = (g_['rk'] < RDF*g_['rd']) if np.isfinite(g_['rd']) \
        and g_['rd'] > 0 else np.zeros(len(g_['th']), bool)
    sigD.append(-DELTA_D*ir.astype(float))
P(f"  WORLD-D delta = {DELTA_D:.5f} (rev 0.11333)")
for nm, sig in (('WORLD-Y', sigY), ('WORLD-R', sigR),
                ('WORLD-D', sigD)):
    dii, doi = cells(GAL_A, sig)
    s1v, s1m, s1g = brute(GAL_A, QUAL_A, sig, 0.15, 5)
    pl = []
    for j in QUAL_A:
        gd_, r = GAL_A[j], sig[j]
        it = gd_['th'] < TH
        pl.append(float(r[it].mean() - r[~it].mean()))
    P(f"  {nm}: plain {np.mean(pl):+.5f}  d_Ii {dii.mean():+.5f}  "
      f"d_Oi {doi.mean():+.5f}  S1 {s1v:+.5f}")
P("  [rev: Y -0.08570/-0.09554/-0.06991/-0.00006; "
  "R -0.08570/-0.08570/+0.00000/-0.08570; "
  "D -0.08570/-0.11333/-0.11333/-0.08508]")

# ---------------- GB-2 z + mixtures ---------------------------------
P("")
P("-- GB-2 z rows + Oi mixtures --")
dii_s, doi_s = cells(GAL_A, RES_BY)
mOi, sOi, nOi = pk(doi_s)
zY = (mOi - (-0.06991))/sOi
zR = (mOi - 0.0)/sOi
zD = (mOi - (-DELTA_D))/sOi
P(f"  sky d_Oi {mOi:+.6f} +/- {sOi:.6f} (n {nOi}); z(Y) {zY:+.2f} "
  f"[rev +4.02], z(R) {zR:+.2f} [rev -1.37], z(D) {zD:+.2f} "
  f"[rev +7.36]")
fy, fy_se = mOi/(-0.06991), sOi/0.06991
fd, fd_se = mOi/(-DELTA_D), sOi/DELTA_D
P(f"  f_y = {fy:.3f} +/- {fy_se:.3f} [rev 0.254 +/- 0.186]; "
  f"f_disk = {fd:.3f} +/- {fd_se:.3f} [rev 0.157 +/- 0.115]")

# ---------------- GB-3 both-sides-16 d_Oi ---------------------------
P("")
P("-- GB-3 both-sides-restricted d_Oi --")
doi16 = []
for j in QUAL_A:
    gd_, r = GAL_A[j], RES_BY[j]
    if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
        continue
    it = gd_['th'] < TH
    ir = gd_['rk'] < RDF*gd_['rd']
    Oo = ~it & ~ir
    Oi = ~it & ir
    if Oo.any() and Oi.any():
        doi16.append(float(r[Oi].mean() - r[Oo].mean()))
m16, s16, n16 = pk(np.array(doi16))
P(f"  d_Oi (both-sides gal): {m16:+.5f} +/- {s16:.5f} (n {n16}) "
  f"[rev -0.03574 +/- 0.02397]; z(Y) {(m16+0.06991)/s16:+.2f} "
  f"[+1.83], z(R) {m16/s16:+.2f} [-1.49], "
  f"z(D) {(m16+DELTA_D)/s16:+.2f} [+3.24]")

# ---------------- GB-4 construction grid ----------------------------
P("")
P("-- GB-4 S1 construction grid (own brute matcher) --")
spots = [(0.10, 5, -2.744), (0.15, 3, -2.061), (0.15, 5, -1.922),
         (0.20, 5, -2.292), (0.30, 10, -2.282), (0.50, 3, -2.146)]
for dy, mw, tgt in spots:
    v, se, m, ng, _, _ = brute_jack(GAL_A, QUAL_A, RES_BY, dy, mw)
    P(f"  dy {dy:.2f}/minw {mw}: {v:+.5f} +/- {se:.5f} = "
      f"{v/se:+.3f} SE [rev {tgt:+.3f}]")
nsur = 0
for dy in (0.10, 0.15, 0.20, 0.25, 0.30, 0.50):
    for mw in (3, 5, 10):
        v, se, _, _, _, _ = brute_jack(GAL_A, QUAL_A, RES_BY, dy, mw)
        if v <= -2*se:
            nsur += 1
P(f"  18-cell survive count = {nsur} [rev 16]")
v, se, _, _, _, _ = brute_jack(GAL_A, QUAL_A, RES_BY, 0.15, 5,
                               weight='galaxy')
P(f"  galaxy-weighted: {v:+.5f} = {v/se:+.3f} SE [rev -0.08497 = "
  f"-1.692]")
v, se, _, _, _, _ = brute_jack(GAL_A, QUAL_A, RES_BY, 0.15, 5,
                               cross_only=True)
P(f"  cross-galaxy-only: {v:+.5f} = {v/se:+.3f} SE [rev -0.09657 = "
  f"-2.006]")

# ---------------- GB-5 shift-null SDs + floor -----------------------
P("")
P("-- GB-5 shift-null SDs, own stream (4000 draws) --")
rng = np.random.default_rng(41905)
nii, noi = [], []
exc = 0
for _ in range(4000):
    sh = []
    for gd_ in GAL_A:
        s = int(rng.integers(0, len(gd_['res'])))
        sh.append(np.roll(gd_['res'], s))
    dii, doi = cells(GAL_A, sh)
    nii.append(dii.mean()); noi.append(doi.mean())
    if dii.mean() <= dii_s.mean():
        exc += 1
nii = np.array(nii); noi = np.array(noi)
P(f"  null SD(d_Ii) = {nii.std(ddof=1):.6f} [rev 0.017083]; "
  f"SD(d_Oi) = {noi.std(ddof=1):.6f} [rev 0.008232]")
P(f"  sky in shift units: Ii {dii_s.mean()/nii.std(ddof=1):+.2f} "
  f"[rev -5.23], Oi {doi_s.mean()/noi.std(ddof=1):+.2f} [rev -2.16]")
P(f"  Ii exceedances in 4000 = {exc} (floor confirmed) "
  f"[rev 0 of 20000]")

# ---------------- GB-6 pool-rebuilt block bootstrap -----------------
P("")
P("-- GB-6 S1 pool-rebuilt galaxy bootstrap (4000, own seed) --")
rng = np.random.default_rng(52061)
bs = []
for _ in range(4000):
    pick = list(rng.choice(QUAL_A, size=len(QUAL_A), replace=True))
    idx = list(range(len(pick)))
    galB = [GAL_A[j] for j in pick]
    resB = [RES_BY[j] for j in pick]
    v, _, _ = brute(galB, idx, resB, 0.15, 5)
    if np.isfinite(v):
        bs.append(v)
bs = np.array(bs)
q = np.quantile(bs, [0.025, 0.5, 0.975])
P(f"  SD = {bs.std(ddof=1):.6f} [rev 0.038776]; "
  f"P(S1 > 0) = {float(np.mean(bs > 0)):.4f} [rev 0.0097]; "
  f"95% ({q[0]:+.4f}, {q[2]:+.4f}) [rev (-0.1685, -0.0174)]")

# ---------------- GB-7 carriers + d_Ii breadth ----------------------
P("")
P("-- GB-7 carriers --")
v0, se0, m0, ng0, jk, tags = brute_jack(GAL_A, QUAL_A, RES_BY,
                                        0.15, 5)
terms = (jk - jk.mean())**2
shares = terms/terms.sum()
order = np.argsort(shares)[::-1]
nm1 = NS['NAMEOF'][GAL_A[tags[order[0]]]['g']]
P(f"  top jackknife-variance carrier: {nm1} share "
  f"{shares[order[0]]:.3f} [rev ESO563-G021 0.378]; top-two "
  f"{shares[order[0]]+shares[order[1]]:.3f} [rev 0.533]")
jdrop = [j for j in tags
         if NS['NAMEOF'][GAL_A[j]['g']] == 'ESO563-G021']
v, _, _ = brute(GAL_A, QUAL_A, RES_BY, 0.15, 5, skip=jdrop[0])
P(f"  drop ESO563-G021: S1 = {v:+.6f} [rev -0.065180]")
neg = int(np.sum(dii_s < 0))
from math import comb
p_sign = sum(comb(len(dii_s), k) for k in range(neg, len(dii_s)+1)) \
    / 2**len(dii_s)
worst = 0.0
for k in range(len(dii_s)):
    mm2 = np.delete(dii_s, k).mean()
    if abs(mm2 - dii_s.mean()) > abs(worst - dii_s.mean()):
        worst = mm2
P(f"  d_Ii breadth: {neg}/{len(dii_s)} negative, sign p = "
  f"{p_sign:.4f} [rev 14/19, 0.0318]; worst single drop -> "
  f"{worst:+.5f} [rev -0.07024]")

# ---------------- GB-8 y-truncation + geometry ----------------------
P("")
P("-- GB-8 S1 y-truncation + inner geometry --")
outly = []
for j in QUAL_A:
    gd_ = GAL_A[j]
    outly.extend(np.log10(gd_['y'])[gd_['th'] >= TH])
ceil = max(outly)
dropped = []
for (j, i, k) in M_S1_A['wins']:
    if len(k) < NS['MINW']:
        dropped.append((NS['NAMEOF'][GAL_A[j]['g']],
                        float(np.log10(GAL_A[j]['y'][i])), len(k)))
P(f"  outer-pool log y ceiling = {ceil:.3f} [rev 1.113]; dropped "
  f"inner points ({len(dropped)}):")
for nm2, ly, w in dropped:
    P(f"    {nm2} log y {ly:+.3f} (window {w}) [rev list incl. "
      f"NGC0891 1.25, NGC5005 1.345/1.152, NGC6946 1.812, "
      f"NGC7814 1.637]")
rr = []
for gd_ in GAL_A:
    if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
        continue
    it = gd_['th'] < TH
    rr.extend((gd_['rk'][it]/gd_['rd']))
rr = np.array(rr)
P(f"  sub-20-arcsec R/Rd: range {rr.min():.3f}-{rr.max():.3f}, "
  f"median {np.median(rr):.3f} [rev 0.078-1.089, 0.485]")

# ---------------- GB-9 signature specificity ------------------------
P("")
P("-- GB-9 signature specificity (400 draws/world, own seeds) --")
def world_d(gal, rng2):
    return [s + e for s, e in zip(sigD, world_noise(gal, rng2))]
rng = np.random.default_rng(63017)
for nm, wf, tgt in (('WORLD-R', world_r, 'ang 0.618 dsk 0.022'),
                    ('WORLD-D', world_d, 'ang 0.000 dsk 0.990'),
                    ('WORLD-Y', world_y, 'ang 0.000 dsk 0.874')):
    ca = cd = cid = 0
    for _ in range(400):
        rb = wf(GAL_A, rng)
        dii, doi = cells(GAL_A, rb)
        mii, sii, _ = pk(dii); moi, soi, _ = pk(doi)
        deep_i = mii <= -2*sii
        cid += deep_i
        if deep_i and abs(moi) < soi:
            ca += 1
        if deep_i and moi <= -2*soi:
            cd += 1
    P(f"  {nm}: P(d_Ii deep) {cid/400:.3f}, sig_ang {ca/400:.3f}, "
      f"sig_dsk {cd/400:.3f}  [rev {tgt}]")

# ---------------- GB-10 dispersion audit ----------------------------
P("")
P("-- GB-10 dispersion audit (trap #28) --")
allr = np.concatenate(RES_BY)
innr = np.concatenate([g_['res'][g_['th'] < TH] for g_ in GAL_A])
P(f"  anchored point SD = {allr.std(ddof=1):.4f} [rev 0.0570]; "
  f"inner-{len(innr)} SD = {innr.std(ddof=1):.4f} [rev 0.1596]")
rng = np.random.default_rng(70443)
sreg = []
for _ in range(400):
    v, _, _, _ = s1_matched(world_r(GAL_A, rng), M_S1_A)
    sreg.append(v)
sreg = np.array(sreg)
sig_g = [g_['res'].std(ddof=1) if len(g_['res']) > 1 else SIG_PT
         for g_ in GAL_A]
def noise_sky(gal, rng2):
    out = []
    for g_, sg in zip(gal, sig_g):
        n = len(g_['res'])
        e = np.empty(n)
        e[0] = rng2.normal(0, sg)
        for i in range(1, n):
            e[i] = 0.6117*e[i-1] + math.sqrt(1-0.6117**2) * \
                rng2.normal(0, sg)
        out.append(e)
    return out
ssky = []
for _ in range(400):
    rb = [(-0.0857)*(g_['th'] < TH) + e
          for g_, e in zip(GAL_A, noise_sky(GAL_A, rng))]
    v, _, _, _ = s1_matched(rb, M_S1_A)
    ssky.append(v)
ssky = np.array(ssky)
P(f"  realized S1 SD: registered world {sreg.std(ddof=1):.4f} "
  f"[rev 0.0272]; sky-matched noise {ssky.std(ddof=1):.4f} "
  f"[rev 0.0257]")
P(f"  power at fixed sky bars (sky-matched noise): "
  f"P(S1 <= -2x0.0493) = {float(np.mean(ssky <= -0.098606)):.3f}; "
  f"P(S1 <= -2x0.0388) = {float(np.mean(ssky <= -0.077552)):.3f} "
  f"[rev ~0.40 / ~0.58 ballpark]")

# ---------------- GB-11 d_Ii vs C17 ---------------------------------
P("")
P("-- GB-11 d_Ii vs the published C17 contrast --")
c17, dii_m = [], []
for gd_, r in zip(GAL_A, RES_BY):
    if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
        continue
    it = gd_['th'] < TH
    ir = gd_['rk'] < RDF*gd_['rd']
    Oo = ~it & ~ir
    if not (Oo.any() and (it & ir).any() and (~it).any() and it.any()):
        continue
    dii_m.append(float(r[it & ir].mean() - r[Oo].mean()))
    c17.append(float(r[it].mean() - r[~it].mean()))
c17 = np.array(c17); dii_m = np.array(dii_m)
r_p = float(np.corrcoef(c17, dii_m)[0, 1])
P(f"  Pearson r = {r_p:.4f} over {len(c17)} galaxies [rev 0.9819 "
  f"over 19]; means {dii_m.mean():+.5f} vs {c17.mean():+.5f}")

# ---------------- GB-12 WORLD-R -2SE clause -------------------------
P("")
rng = np.random.default_rng(81559)
c2 = 0
for _ in range(400):
    v, se, mp_, ng_ = s1_matched(world_r(GAL_A, rng), M_S1_A)
    if np.isfinite(v) and v <= -2*se:
        c2 += 1
P(f"-- GB-12 P(S1 <= -2 SE_draw | WORLD-R, registered noise) = "
  f"{c2/400:.3f} [rev 0.817-0.863]")

# ---------------- GB-13 the overlap-match qualifier -----------------
P("")
oc = []
for gd_, r in zip(GAL_A, RES_BY):
    mi = gd_['th'] < TH; mo = ~mi
    if not (mi.any() and mo.any()):
        continue
    yi, yo = gd_['y'][mi], gd_['y'][mo]
    lo = max(yi.min(), yo.min()); hi = min(yi.max(), yo.max())
    if not (hi > lo):
        continue
    ki = mi & (gd_['y'] >= lo) & (gd_['y'] <= hi)
    ko = mo & (gd_['y'] >= lo) & (gd_['y'] <= hi)
    if ki.any() and ko.any():
        oc.append((NS['NAMEOF'][gd_['g']],
                   float(r[ki].mean() - r[ko].mean())))
P(f"-- GB-13 anchored overlap-match qualifiers: {oc} "
  f"[rev ESO563-G021 -0.13788]")

P("")
P(f"GB complete; wall-clock {(time.time()-t0)/60:.1f} min")
with open('data/round52_gb.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(L) + "\n")
print("\nwrote: data/round52_gb.txt")
