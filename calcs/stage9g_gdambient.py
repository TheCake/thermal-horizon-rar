"""STAGE 9G — THE AMBIENT CONTROL ON THE GD DIAL (galaxy).
Pre-registered BEFORE any run (measurement round; NO credence
movement).

Is the GD dial-tension carried by measured ENVIRONMENT?  Chae+21
Table 3 per-galaxy ambients (PRIMARY = maxclust; noclust co-read).
Reads: A between-type (median log eN, permutation test 20000 rng 7);
B within-GD tracking (frozen point-median split, lam_hat point fits
+ paired galaxy bootstrap NBOOT=200 rng 71; halves < 5 => skip;
point halves < 8 => DESCRIPTIVE-ONLY); C DD control (NBOOT=100).
Gates: G9G-0 verbatim-lift (counts 38/111 + OFF probe 1e-6);
G9G-1 fast-vs-verbatim 6 probes <= 1e-6; G9G-2 match audit (>= 60
total, >= 12 GD else ABORT UNDERPOWERED); G9G-3 rng fingerprints.
Bars (locked, ordered): A1 ENV-SPECIAL (|dmed| >= 0.30 dex AND
perm P <= 0.05, maxclust); B1 ENV-TRACKS (non-descriptive AND
|D_point| >= 0.50 AND min flip P <= 0.05); SEVENTH-CONTROL-PASSED
(neither, B non-descriptive); GRAY-CARRIED else.
Output: data/stage9g_gdambient.txt
"""
import csv as _csv
import glob, math, os, re, time
import numpy as np
from scipy.optimize import minimize

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)

meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        name = t[0]
        D, eD = float(t[2]), float(t[3])
        inc, einc = float(t[5]), float(t[6])
        q = int(t[17])
        meta[name] = (inc, q, D, eD, einc)
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
svert = {}
gname = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 1.0, 0.0, 0.0))
    if inc < 30 or q > 2: continue
    kept += 1
    s_d2 = (eD/max(D, 1e-6))/LN10
    s_i2 = 2.0*math.radians(einc)/max(math.tan(math.radians(inc)), 1e-6)/LN10
    svert[gi] = math.sqrt(s_d2**2 + s_i2**2)
    gname[gi] = name
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC; gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        g_gas.append(gg); g_dsk.append(gd); g_bul.append(gb)
        gobs.append(Vo*Vo/R*KPC)
        sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs = np.log10(gobs)

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_lam(y, lam):
    return (1.0-lam)*nu_standard(y) + lam*nu_be(y)

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

t00 = time.time()
P("9G THE AMBIENT CONTROL ON THE GD DIAL (pre-reg committed BEFORE "
  "any run; measurement round; NO credence movement)")

allg = np.unique(gal_id)
NGAL = gal_id.max()+1
gdfrac = np.zeros(NGAL)
for g_ in allg:
    m = gal_id == g_
    gdfrac[g_] = float(np.mean(g_gas[m] > g_dsk[m] + g_bul[m]))
gd_set = np.array([g_ for g_ in allg if gdfrac[g_] >= 0.5])
dd_set = np.array([g_ for g_ in allg if gdfrac[g_] < 0.5])
SV = np.array([svert.get(g_, 0.0) for g_ in range(NGAL)])
GIDX = [np.where(gal_id == g_)[0] for g_ in range(NGAL)]
ALL_LIST = [int(g_) for g_ in allg]

txt8sc = open('data/stage8sc_gddist.txt').read()
OFFV = float(re.search(r"G8Sc-0 OFF-branch identity: ([-\d.]+) vs",
                       txt8sc).group(1))

def m2ll_vert(th, lam, gset, von=True):
    la0, f, s_int = th
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    nv = nu_lam(gN/a0, lam)
    if np.min(nv) <= 0.05: return 1e12
    gm = gN*nv
    v = sig2 + s_int*s_int
    a = lgobs - np.log10(gm)
    out = 0.0
    for g_ in gset:
        ji = GIDX[g_]
        if len(ji) == 0: continue
        aj, vj = a[ji], v[ji]
        s_ = SV[g_]
        if (not von) or s_ <= 1e-4:
            out += float(np.sum(aj*aj/vj + np.log(vj)))
        else:
            iv = 1.0/vj
            Siv = float(np.sum(iv))
            Sa = float(np.sum(aj*iv))
            out += (float(np.sum(aj*aj*iv))
                    - Sa*Sa/(Siv + 1.0/(s_*s_))
                    + float(np.sum(np.log(vj)))
                    + math.log(1.0 + s_*s_*Siv))
    return out

def make_instance(gset_list):
    cats, labs, svals = [], [], []
    li = 0
    for g_ in gset_list:
        ji = GIDX[g_]
        if len(ji) == 0: continue
        cats.append(ji)
        labs.append(np.full(len(ji), li, dtype=np.int64))
        svals.append(SV[g_])
        li += 1
    return (np.concatenate(cats), np.concatenate(labs),
            np.array(svals), li)

def m2ll_fast(th, lam, cat, lab, svec, ninst):
    la0, f, s_int = th
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    nv = nu_lam(gN/a0, lam)
    if np.min(nv) <= 0.05: return 1e12
    gm = gN*nv
    v = sig2 + s_int*s_int
    a = lgobs - np.log10(gm)
    ac, vc_ = a[cat], v[cat]
    iv = 1.0/vc_
    Siv = np.bincount(lab, weights=iv, minlength=ninst)
    Sa  = np.bincount(lab, weights=ac*iv, minlength=ninst)
    Saa = np.bincount(lab, weights=ac*ac*iv, minlength=ninst)
    Slv = np.bincount(lab, weights=np.log(vc_), minlength=ninst)
    s2 = svec*svec
    on = svec > 1e-4
    out = np.where(on,
                   Saa - Sa*Sa/(Siv + 1.0/np.maximum(s2, 1e-30))
                   + Slv + np.log(1.0 + s2*Siv),
                   Saa + Slv)
    return float(np.sum(out))

LGB = np.round(np.arange(-2.0, 1.501, 0.25), 3)   # 15 pts
def lam_hat_fast(cat, lab, svec, ninst):
    prof = []; th = None
    for lam in LGB:
        starts = (([list(th)] if th is not None else [])
                  + [[math.log10(A0_FID), 1.0, 0.08]])
        best = None
        for th0 in starts:
            b = minimize(lambda t: m2ll_fast(t, lam, cat, lab, svec,
                                             ninst), th0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6,
                                      fatol=1e-7))
            if best is None or b.fun < best.fun: best = b
        prof.append(best.fun); th = best.x
    prof = np.array(prof)
    i = int(np.argmin(prof))
    lam_hat = float(LGB[i])
    lo_edge = (i == 0); hi_edge = (i == len(LGB)-1)
    if not (lo_edge or hi_edge):
        x3, y3 = LGB[i-1:i+2], prof[i-1:i+2]
        c2_, c1_, _ = np.polyfit(x3, y3, 2)
        if c2_ > 0:
            lh = -c1_/(2*c2_)
            if LGB[i-1] <= lh <= LGB[i+1]: lam_hat = float(lh)
    return lam_hat, lo_edge, hi_edge

# ---------------- gates ----------------
th_p = [math.log10(A0_FID), 1.0, 0.08]
ok_n = (len(gd_set) == 38) and (len(dd_set) == 111)
va = m2ll_vert(th_p, 0.7, list(gd_set), von=False)
ok_off = abs(va - OFFV) <= 1e-6
g0_ok = ok_n and ok_off
P(f"G9G-0 verbatim-lift regression: GD {len(gd_set)} / DD "
  f"{len(dd_set)} (want 38/111) -> {'PASS' if ok_n else 'FAIL'}; "
  f"OFF probe {va:.9f} vs 8S-c {OFFV:.9f} (d={va-OFFV:.2e}) -> "
  f"{'PASS' if ok_off else 'FAIL'}")

rngp = np.random.default_rng(7)
dup = [ALL_LIST[j] for j in rngp.integers(0, len(ALL_LIST),
                                          size=len(ALL_LIST))]
probe_sets = [ALL_LIST, [int(g_) for g_ in gd_set], dup]
g1_ok = True; dmx = 0.0
for pi in range(6):
    thp = [-10.5+rngp.random(), 0.4+2.0*rngp.random(),
           0.01+0.29*rngp.random()]
    lam = -2.0+3.5*rngp.random()
    gl = probe_sets[pi % 3]
    slow = m2ll_vert(thp, lam, gl, von=True)
    cat, lab, svec, ninst = make_instance(gl)
    fast = m2ll_fast(thp, lam, cat, lab, svec, ninst)
    if slow < 1e11:
        dmx = max(dmx, abs(slow-fast))
g1_ok = dmx <= 1e-6
P(f"G9G-1 fast-vs-verbatim cross-gate (6 probes): max|d| = "
  f"{dmx:.2e} -> {'PASS' if g1_ok else 'FAIL'}")

# ---------------- the Chae match ----------------
chae = {}
with open('data/chae2021_table3.csv') as f:
    for row in _csv.DictReader(l for l in f if not l.startswith('#')):
        chae[row['galaxy']] = (float(row['log_eN_maxclust']),
                               float(row['log_eN_noclust']))
gd_m = [int(g_) for g_ in gd_set if gname.get(g_) in chae]
dd_m = [int(g_) for g_ in dd_set if gname.get(g_) in chae]
nm = len(gd_m) + len(dd_m)
g2_ok = (nm >= 60) and (len(gd_m) >= 12)
P(f"G9G-2 match audit: GD matched {len(gd_m)}, DD matched "
  f"{len(dd_m)}, total {nm} -> "
  f"{'PASS' if g2_ok else 'FAIL (UNDERPOWERED)'}")

if not (g0_ok and g1_ok and g2_ok):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    with open('data/stage9g_gdambient.txt', 'w') as f:
        f.write("\n".join(L)+"\n")
    raise SystemExit(0)
P("GATES: G9G-0/1/2 ALL PASS")
P("")

eN_max = {g_: chae[gname[g_]][0] for g_ in gd_m + dd_m}
eN_noc = {g_: chae[gname[g_]][1] for g_ in gd_m + dd_m}

# ---------------- read A: between-type ----------------
for tag, eN in (('maxclust', eN_max), ('noclust', eN_noc)):
    vg = np.array([eN[g_] for g_ in gd_m])
    vd = np.array([eN[g_] for g_ in dd_m])
    dmed = float(np.median(vg) - np.median(vd))
    pool = np.concatenate([vg, vd])
    rngA = np.random.default_rng(7)
    cnt = 0; NPERM = 20000
    for _ in range(NPERM):
        pi_ = rngA.permutation(len(pool))
        dm_ = np.median(pool[pi_[:len(vg)]]) - np.median(pool[pi_[len(vg):]])
        if abs(dm_) >= abs(dmed): cnt += 1
    pperm = cnt/NPERM
    P(f"[A {tag:8}] med(GD) = {np.median(vg):+.3f}, med(DD) = "
      f"{np.median(vd):+.3f}, dmed = {dmed:+.3f} dex; perm P = "
      f"{pperm:.4f}" + ("  <- PRIMARY" if tag == 'maxclust' else ""))
    if tag == 'maxclust':
        dmed_pri, pperm_pri = dmed, pperm

# ---------------- read B: within-GD tracking ----------------
thrB = float(np.median([eN_max[g_] for g_ in gd_m]))
hiB = [g_ for g_ in gd_m if eN_max[g_] > thrB]
loB = [g_ for g_ in gd_m if eN_max[g_] <= thrB]
descriptive = min(len(hiB), len(loB)) < 8
P("")
P(f"[B] GD split at frozen median log eN(maxclust) = {thrB:+.3f}: "
  f"hi {len(hiB)} / lo {len(loB)}"
  + ("  DESCRIPTIVE-ONLY (half < 8)" if descriptive else ""))
lh_hi, elo1, ehi1 = lam_hat_fast(*make_instance(hiB))
lh_lo, elo2, ehi2 = lam_hat_fast(*make_instance(loB))
Dpoint = lh_hi - lh_lo
P(f"[B] point fits: lam_hat(hi eN) = {lh_hi:+.3f} "
  f"(edge lo/hi {int(elo1)}/{int(ehi1)}), lam_hat(lo eN) = "
  f"{lh_lo:+.3f} (edge lo/hi {int(elo2)}/{int(ehi2)}); D = "
  f"{Dpoint:+.3f}")

NBOOT_B = 200
rng = np.random.default_rng(71)
Ds = []; skipsB = 0
eloH = ehiH = eloL = ehiL = 0
for r in range(NBOOT_B):
    draw = rng.integers(0, len(gd_m), size=len(gd_m))
    rep = [gd_m[j] for j in draw]
    if r == 0:
        P(f"G9G-3 first-replicate fingerprint (B): sum(draw)="
          f"{int(np.sum(draw))}, unique={len(np.unique(draw))}")
    hi_ = [g_ for g_ in rep if eN_max[g_] > thrB]
    lo_ = [g_ for g_ in rep if eN_max[g_] <= thrB]
    if min(len(hi_), len(lo_)) < 5:
        skipsB += 1; continue
    lh1, e1, e2 = lam_hat_fast(*make_instance(hi_))
    lh2, e3, e4 = lam_hat_fast(*make_instance(lo_))
    Ds.append(lh1-lh2)
    eloH += int(e1); ehiH += int(e2); eloL += int(e3); ehiL += int(e4)
    if (r+1) % 50 == 0:
        P(f"  ... B {r+1}/{NBOOT_B} ({(time.time()-t00)/60:.1f} min)")
Ds = np.array(Ds)
pge = float(np.mean(Ds >= 0)); ple = float(np.mean(Ds <= 0))
qb = np.percentile(Ds, [5, 50, 95])
P(f"[B] bootstrap: reps {len(Ds)} (skips {skipsB}); D pct 5/50/95 = "
  f"{qb[0]:+.3f}/{qb[1]:+.3f}/{qb[2]:+.3f}; P(D >= 0) = {pge:.3f}, "
  f"P(D <= 0) = {ple:.3f}; edges hi-set lo/hi = {eloH}/{ehiH}, "
  f"lo-set lo/hi = {eloL}/{ehiL}")

# ---------------- read C: DD control ----------------
thrC = float(np.median([eN_max[g_] for g_ in dd_m]))
hiC = [g_ for g_ in dd_m if eN_max[g_] > thrC]
loC = [g_ for g_ in dd_m if eN_max[g_] <= thrC]
P("")
P(f"[C] DD split at frozen median log eN(maxclust) = {thrC:+.3f}: "
  f"hi {len(hiC)} / lo {len(loC)}")
lh_hiC, _, _ = lam_hat_fast(*make_instance(hiC))
lh_loC, _, _ = lam_hat_fast(*make_instance(loC))
P(f"[C] point fits: lam_hat(hi eN) = {lh_hiC:+.3f}, lam_hat(lo eN) "
  f"= {lh_loC:+.3f}; D = {lh_hiC-lh_loC:+.3f}")
NBOOT_C = 100
DsC = []; skipsC = 0
for r in range(NBOOT_C):
    draw = rng.integers(0, len(dd_m), size=len(dd_m))
    rep = [dd_m[j] for j in draw]
    if r == 0:
        P(f"G9G-3 first-replicate fingerprint (C): sum(draw)="
          f"{int(np.sum(draw))}, unique={len(np.unique(draw))}")
    hi_ = [g_ for g_ in rep if eN_max[g_] > thrC]
    lo_ = [g_ for g_ in rep if eN_max[g_] <= thrC]
    if min(len(hi_), len(lo_)) < 5:
        skipsC += 1; continue
    lh1, _, _ = lam_hat_fast(*make_instance(hi_))
    lh2, _, _ = lam_hat_fast(*make_instance(lo_))
    DsC.append(lh1-lh2)
    if (r+1) % 50 == 0:
        P(f"  ... C {r+1}/{NBOOT_C} ({(time.time()-t00)/60:.1f} min)")
DsC = np.array(DsC)
qc = np.percentile(DsC, [5, 50, 95])
P(f"[C] bootstrap: reps {len(DsC)} (skips {skipsC}); D pct 5/50/95 "
  f"= {qc[0]:+.3f}/{qc[1]:+.3f}/{qc[2]:+.3f}; P(D >= 0) = "
  f"{float(np.mean(DsC >= 0)):.3f}")
P("")

# ---------------- verdict (locked grammar) ----------------
a1 = (abs(dmed_pri) >= 0.30) and (pperm_pri <= 0.05)
b1 = (not descriptive) and (abs(Dpoint) >= 0.50) \
     and (min(pge, ple) <= 0.05)
fired = []
if a1: fired.append("A1 ENV-SPECIAL")
if b1: fired.append("B1 ENV-TRACKS")
if fired:
    P("==> 9G VERDICT (locked grammar): " + " + ".join(fired)
      + " - the GD dial has a measured environmental axis; the "
        "environmental reading of the type-split is LIVE (AMB-"
        "adjacent); interpretation deferred to the booked entry.")
elif not descriptive:
    P("==> 9G VERDICT (locked grammar): SEVENTH-CONTROL-PASSED - "
      "measured environment neither distinguishes GD from DD nor "
      "orders the dial within GD; the tension is NOT carried by "
      "the Chae ambient axis.")
else:
    P("==> 9G VERDICT (locked grammar): GRAY-CARRIED - read B "
      "descriptive; rows stand as measurements.")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t00)/60:.1f} min)")

with open('data/stage9g_gdambient.txt', 'w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage9g_gdambient.txt")
