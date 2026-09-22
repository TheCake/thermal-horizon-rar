"""STAGE 10Y GA (blind half; committed BEFORE the review round's
report is read -- the 87a4676 protocol, 18th execution). Re-computes
the sky read's load-bearing numbers by INDEPENDENT routes:
  GA-1 S1 matched contrast by brute-force double loop (no
       precomputed matcher), value + jackknife -> exact.
  GA-2 d_Ii/d_Oi by global cell labeling + per-galaxy grouping
       -> exact.
  GA-3 circular-shift p (Ii, Oi) with an OWN rng stream -> same
       order (binomial MC band printed).
  GA-4 Spearman rho via scipy.stats.spearmanr -> exact rho.
  GA-5 grammar re-derivation from the GA numbers -> B4.
  GA-6 WORLD-Y bias replication with an own seed -> |bias| <= 0.010.
Appends to data/stage10y_ga.txt.
"""
import math, os, sys, time
import numpy as np
from scipy.stats import spearmanr

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

np_ = np
GAL_A, GAL_F = NS['GAL_A'], NS['GAL_F']
QUAL_A = NS['QUAL_A']
TH_STAR, RD_FAC = NS['TH_STAR'], NS['RD_FAC']
DY, MINW = NS['DY_MATCH'], NS['MINW']

L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

t0 = time.time()
P("")
P("=" * 68)
P("STAGE 10Y GA (blind; independent routes)")
P("")

# ---------------- GA-1 brute-force matched contrast -----------------
def brute_matched(gal, qual, skip=None):
    """No matcher structure: rescan windows from scratch."""
    A = []
    for j in qual:
        if skip is not None and j == skip:
            continue
        gd_ = gal[j]
        ly = np.log10(gd_['y'])
        for i in np.where(gd_['th'] < TH_STAR)[0]:
            A.append((j, i, ly[i]))
    vals, gals = [], set()
    for (j, i, la) in A:
        pool = []
        for j2 in qual:
            if skip is not None and j2 == skip:
                continue
            gd2 = gal[j2]
            ly2 = np.log10(gd2['y'])
            for i2 in np.where(gd2['th'] >= TH_STAR)[0]:
                if abs(ly2[i2] - la) <= DY:
                    pool.append(gd2['res'][i2])
        if len(pool) >= MINW:
            vals.append(gal[j]['res'][i] - float(np.mean(pool)))
            gals.add(j)
    if not vals:
        return float('nan'), 0, set()
    return float(np.mean(vals)), len(vals), gals

vA, mA, gsA = brute_matched(GAL_A, QUAL_A)
jk = []
for j in sorted(gsA):
    v, _, _ = brute_matched(GAL_A, QUAL_A, skip=j)
    if np.isfinite(v):
        jk.append(v)
jk = np.array(jk)
n = len(jk)
seA = math.sqrt((n - 1)/n*float(np.sum((jk - jk.mean())**2)))
P(f"GA-1 anch S1 matched (brute): {vA:+.4f} +/- {seA:.4f} "
  f"({mA} pts, {len(gsA)} gal)  [sky -0.0947 +/- 0.0493 / 22 / 16]")

# ---------------- GA-2 cells by global labeling ---------------------
def cells_indep(gal):
    rows = {}
    for gd_ in gal:
        if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
            continue
        lab = np.where(gd_['th'] < TH_STAR,
                       np.where(gd_['rk'] < RD_FAC*gd_['rd'], 'Ii', 'Io'),
                       np.where(gd_['rk'] < RD_FAC*gd_['rd'], 'Oi', 'Oo'))
        if not (lab == 'Oo').any():
            continue
        base = float(gd_['res'][lab == 'Oo'].mean())
        for c in ('Ii', 'Oi'):
            if (lab == c).any():
                rows.setdefault(c, []).append(
                    float(gd_['res'][lab == c].mean()) - base)
    out = {}
    for c in ('Ii', 'Oi'):
        v = np.array(rows.get(c, []))
        out[c] = (float(v.mean()), float(v.std(ddof=1)/math.sqrt(len(v))),
                  len(v)) if len(v) >= 2 else (float('nan'),)*2 + (len(v),)
    return out

cA = cells_indep(GAL_A)
P(f"GA-2 anch d_Ii: {cA['Ii'][0]:+.4f} +/- {cA['Ii'][1]:.4f} "
  f"(n {cA['Ii'][2]})  [sky -0.0893 +/- 0.0349 / 19]")
P(f"GA-2 anch d_Oi: {cA['Oi'][0]:+.4f} +/- {cA['Oi'][1]:.4f} "
  f"(n {cA['Oi'][2]})  [sky -0.0177 +/- 0.0130 / 55]")
cF = cells_indep(GAL_F)
P(f"GA-2 flow d_Ii: {cF['Ii'][0]:+.4f} +/- {cF['Ii'][1]:.4f} "
  f"(n {cF['Ii'][2]})  [sky -0.0021 +/- 0.0204 / 43]")

# ---------------- GA-3 shift p, own stream --------------------------
rng = np.random.default_rng(77120)
obs_ii, obs_oi = cA['Ii'][0], cA['Oi'][0]
cii = coi = 0
NDR = 2000
for _ in range(NDR):
    sh = []
    for gd_ in GAL_A:
        s = int(rng.integers(0, len(gd_['res'])))
        sh.append(np.roll(gd_['res'], s))
    rows_ii, rows_oi = [], []
    for gd_, r in zip(GAL_A, sh):
        if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
            continue
        inn_th = gd_['th'] < TH_STAR
        inn_rd = gd_['rk'] < RD_FAC*gd_['rd']
        Oo = ~inn_th & ~inn_rd
        if not Oo.any():
            continue
        base = float(r[Oo].mean())
        Ii = inn_th & inn_rd
        Oi = ~inn_th & inn_rd
        if Ii.any():
            rows_ii.append(float(r[Ii].mean()) - base)
        if Oi.any():
            rows_oi.append(float(r[Oi].mean()) - base)
    if len(rows_ii) >= 2 and float(np.mean(rows_ii)) <= obs_ii:
        cii += 1
    if len(rows_oi) >= 2 and float(np.mean(rows_oi)) <= obs_oi:
        coi += 1
pii = (cii + 1)/(NDR + 1); poi = (coi + 1)/(NDR + 1)
P(f"GA-3 own-stream shift p: Ii {pii:.4f} [sky 0.0005], "
  f"Oi {poi:.4f} [sky 0.0190] (binomial MC scatter expected)")

# ---------------- GA-4 spearman cross-check -------------------------
d, thmin, dist = [], [], []
for gd_ in GAL_A:
    inn = gd_['res'][gd_['th'] < TH_STAR]
    out = gd_['res'][gd_['th'] >= TH_STAR]
    if len(inn) and len(out):
        d.append(float(inn.mean() - out.mean()))
        thmin.append(float(gd_['th'].min()))
        dist.append(gd_['D'])
r1 = spearmanr(d, thmin)
r2 = spearmanr(d, dist)
P(f"GA-4 anch spearman: thmin rho {r1.statistic:+.3f} "
  f"[sky +0.019], dist rho {r2.statistic:+.3f} [sky -0.318] "
  f"(scipy analytic p {r1.pvalue:.3f}/{r2.pvalue:.3f}; perm p stays "
  f"operative)")

# ---------------- GA-5 grammar re-derivation ------------------------
s1_pop = (mA >= NS['NMIN_MPT']) and (len(gsA) >= NS['NMIN_S1'])
survives = s1_pop and vA <= -2*seA
collapses = s1_pop and abs(vA) < seA
ii_deep = cA['Ii'][2] >= NS['NMIN_CELL'] and cA['Ii'][0] <= -2*cA['Ii'][1]
oi_deep = cA['Oi'][2] >= NS['NMIN_CELL'] and cA['Oi'][0] <= -2*cA['Oi'][1]
oi_flat = cA['Oi'][2] >= NS['NMIN_CELL'] and abs(cA['Oi'][0]) < cA['Oi'][1]
sig_ang = ii_deep and oi_flat
sig_dsk = ii_deep and oi_deep
br = 'B1' if (s1_pop and collapses) else 'B4'
if survives and sig_ang:
    br = 'B2?'
if survives and sig_dsk:
    br = 'B3?'
P(f"GA-5 grammar: populated {s1_pop}, survives {survives}, "
  f"collapses {collapses}, sig_ang {sig_ang}, sig_dsk {sig_dsk} "
  f"-> BRANCH {br}  [sky B4]")

# ---------------- GA-6 WORLD-Y bias, own seed -----------------------
rng2 = np.random.default_rng(88231)
vy = []
for _ in range(500):
    rb = NS['world_y'](GAL_A, rng2)
    mm, sm, mp_, ng_ = NS['s1_matched'](rb, NS['M_S1_A'])
    if np.isfinite(mm):
        vy.append(mm)
vy = np.array(vy)
P(f"GA-6 WORLD-Y bias (own seed): {vy.mean():+.4f} +/- "
  f"{vy.std(ddof=1)/math.sqrt(len(vy)):.4f} "
  f"[gates -0.0017 +/- 0.0013; bar |bias| <= 0.010]")

P("")
P(f"GA complete; wall-clock {(time.time()-t0)/60:.1f} min")
with open('data/stage10y_ga.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(L) + "\n")
print("\nwrote: data/stage10y_ga.txt")
