"""STAGE 8V — THE GALAXY-LEVEL BOOTSTRAP OF THE GD DIAL-TENSION
(the correlation control).  Pre-registered BEFORE any run.

The 8S-c suspect ledger ranks the within-curve correlation object
(7B lag-1 rho ~ 0.87) first: point-level D1 intervals overcount.
The honest resampling unit is the GALAXY.  PAIRED design: NBOOT = 300
replicates, rng 71; each replicate draws 149 galaxies with
replacement from the kept 149 (GD 38 + DD 111); per replicate fit
lam-hat on the full draw AND on its GD-member instances AND on its
DD-member instances (duplicated galaxies enter the marginal once per
instance — independent-intercept bootstrap semantics).  Statistic:
D_rep = lam_GD_rep − lam_FULL_rep (paired).  Estimator: coarse lam
grid −2.0..1.5 step 0.25 (15 pts), warm-chained Nelder-Mead (warm +
one cold start per lam), parabolic refine at interior minima,
clipped to the bracket; a FAST bincount-vectorized m2ll is the
bootstrap engine, cross-gated against the VERBATIM 8S-c expression
(the GB0w new-MODE precedent).  Skip rule: subset < 10 instances =>
NaN + count.

Gates (any FAIL => STOP; amendment pre-quote):
  G8V-0 verbatim-lift regression — the OFF-branch probe equals
        8S-c's printed value at 1e-6; GD/DD counts 38/111 exact.
  G8V-1 fast-vs-verbatim cross-gate at 20 random theta probes incl.
        a duplicated-instance list, max|d| <= 1e-6.
  G8V-2 identity-replicate estimator check — coarse lam-hat within
        0.10 of 8S-c's fine values for FULL/GD/DD.
  G8V-3 rng bookkeeping (first-replicate fingerprint printed).
Bars (locked): B1 TENSION-ROBUST iff P(D >= 0) <= 0.05 AND the 95th
percentile of lam_GD_rep < the 8S-c FULL point (ONE-SIDED,
explicitly — the 8S-b grammar lesson).  B2 CORRELATION-ABSORBED iff
P(D >= 0) >= 0.20.  B3 GRAY-CARRIED else.  Edge rule (pre-stated
asymmetry): GD low-edge (−2.0) censoring RAISES lam_GD toward the
dial, i.e. biases AGAINST B1 — a B1 pass is conservative; report
edge fractions per set; if GD edge fraction > 30%, append
ESTIMATE-CENSORED.  Co-reads: percentiles (5/50/95) of
lam_FULL/lam_GD/lam_DD/D; P(lam_GD >= lam_DD); FULL high-edge
fraction.  NO credence movement (measurement round; pre-stated).
Output: data/stage8v_gdboot.txt
"""
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
P("8V THE GALAXY-LEVEL BOOTSTRAP OF THE GD DIAL-TENSION (pre-reg "
  "committed BEFORE any run; measurement round; NO credence "
  "movement)")

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

# ---- 8S-c targets (lnL-grade + estimator regression anchors) ----
txt8sc = open('data/stage8sc_gddist.txt').read()
TGT = {m.group(1): float(m.group(2)) for m in
       re.finditer(r"\[(FULL|GD|DD)\s+vertON\] lam_hat = ([-\d.]+)",
                   txt8sc)}
OFFV = float(re.search(r"G8Sc-0 OFF-branch identity: ([-\d.]+) vs",
                       txt8sc).group(1))
assert set(TGT) == {'FULL', 'GD', 'DD'}, TGT

# ---- VERBATIM 8S-c objective (the slow reference) ----
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
            # VERBATIM flat expression (the s=0 / OFF branch)
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

# ---- FAST bincount engine (new MODE; cross-gated G8V-1) ----
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
P(f"G8V-0 verbatim-lift regression: GD {len(gd_set)} / DD "
  f"{len(dd_set)} (want 38/111) -> {'PASS' if ok_n else 'FAIL'}; "
  f"OFF probe {va:.9f} vs 8S-c {OFFV:.9f} (d={va-OFFV:.2e}) -> "
  f"{'PASS' if ok_off else 'FAIL'}")

rngp = np.random.default_rng(7)
dup = [ALL_LIST[j] for j in rngp.integers(0, len(ALL_LIST),
                                          size=len(ALL_LIST))]
probe_sets = [ALL_LIST, [int(g_) for g_ in gd_set], dup]
g1_ok = True; dmx = 0.0
for pi in range(20):
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
P(f"G8V-1 fast-vs-verbatim cross-gate (20 probes incl. duplicated "
  f"list): max|d| = {dmx:.2e} -> {'PASS' if g1_ok else 'FAIL'}")

g2_ok = True
id_hats = {}
for tag, gl in (('FULL', ALL_LIST),
                ('GD', [int(g_) for g_ in gd_set]),
                ('DD', [int(g_) for g_ in dd_set])):
    cat, lab, svec, ninst = make_instance(gl)
    lh, lo_e, hi_e = lam_hat_fast(cat, lab, svec, ninst)
    okt = abs(lh - TGT[tag]) <= 0.10
    g2_ok &= okt
    id_hats[tag] = lh
    P(f"G8V-2 identity replicate [{tag}]: coarse lam_hat = {lh:.3f} "
      f"vs 8S-c fine {TGT[tag]:.3f} -> {'PASS' if okt else 'FAIL'}")

if not (g0_ok and g1_ok and g2_ok):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    with open('data/stage8v_gdboot.txt', 'w') as f:
        f.write("\n".join(L)+"\n")
    raise SystemExit(0)
P("GATES: G8V-0/1/2 ALL PASS")
P("")

# ---------------- the bootstrap ----------------
NBOOT = 300
rng = np.random.default_rng(71)
res = {'FULL': [], 'GD': [], 'DD': []}
edges_lo = {'FULL': 0, 'GD': 0, 'DD': 0}
edges_hi = {'FULL': 0, 'GD': 0, 'DD': 0}
skips = 0
for r in range(NBOOT):
    draw = rng.integers(0, len(ALL_LIST), size=len(ALL_LIST))
    rep = [ALL_LIST[j] for j in draw]
    if r == 0:
        P(f"G8V-3 first-replicate fingerprint: sum(draw)="
          f"{int(np.sum(draw))}, unique={len(np.unique(draw))}, "
          f"first12(sorted)={sorted(rep)[:12]}")
    subs = {'FULL': rep,
            'GD': [g for g in rep if gdfrac[g] >= 0.5],
            'DD': [g for g in rep if gdfrac[g] < 0.5]}
    for tag in ('FULL', 'GD', 'DD'):
        ls = subs[tag]
        if len(ls) < 10:
            res[tag].append(np.nan); skips += 1; continue
        cat, lab, svec, ninst = make_instance(ls)
        lh, lo_e, hi_e = lam_hat_fast(cat, lab, svec, ninst)
        res[tag].append(lh)
        edges_lo[tag] += int(lo_e); edges_hi[tag] += int(hi_e)
    if (r+1) % 50 == 0:
        P(f"  ... {r+1}/{NBOOT} replicates "
          f"({(time.time()-t00)/60:.1f} min)")

fullv = np.array(res['FULL']); gdv = np.array(res['GD'])
ddv = np.array(res['DD'])
val = np.isfinite(fullv) & np.isfinite(gdv)
dva = gdv[val] - fullv[val]
pge = float(np.mean(dva >= 0))
p95_gd = float(np.percentile(gdv[np.isfinite(gdv)], 95))
pgd_dd = float(np.mean(gdv[np.isfinite(gdv) & np.isfinite(ddv)]
                       >= ddv[np.isfinite(gdv) & np.isfinite(ddv)]))
P("")
P(f"replicates: {NBOOT}; skips {skips}; valid pairs {int(val.sum())}")
for tag, arr in (('FULL', fullv), ('GD', gdv), ('DD', ddv)):
    q = np.percentile(arr[np.isfinite(arr)], [5, 50, 95])
    P(f"[{tag:4}] lam_hat pct 5/50/95 = {q[0]:+.3f}/{q[1]:+.3f}/"
      f"{q[2]:+.3f}; edge lo/hi = {edges_lo[tag]}/{edges_hi[tag]} "
      f"({100*edges_lo[tag]/NBOOT:.0f}%/{100*edges_hi[tag]/NBOOT:.0f}%)")
qd = np.percentile(dva, [5, 50, 95])
P(f"[D = GD-FULL] pct 5/50/95 = {qd[0]:+.3f}/{qd[1]:+.3f}/"
  f"{qd[2]:+.3f}; P(D >= 0) = {pge:.4f}")
P(f"P(lam_GD >= lam_DD) = {pgd_dd:.4f}; 95th pct lam_GD = "
  f"{p95_gd:+.3f} vs FULL point {TGT['FULL']:+.3f}")
P("")
gd_edge_frac = edges_lo['GD']/NBOOT
cens = "; ESTIMATE-CENSORED (GD low-edge fraction > 30%)" \
       if gd_edge_frac > 0.30 else ""
b1 = (pge <= 0.05) and (p95_gd < TGT['FULL'])
b2 = pge >= 0.20
if b1:
    P("==> 8V VERDICT (locked grammar): TENSION-ROBUST - the "
      "galaxy-level (correlation-honest) resampling keeps the "
      "gas-dominated optimum below the full-sample dial (one-sided, "
      "pre-stated; GD low-edge censoring biases AGAINST this bar, "
      "so the pass is conservative)" + cens + ".")
elif b2:
    P("==> 8V VERDICT (locked grammar): CORRELATION-ABSORBED - the "
      "point-level D1 overcounted; at the honest galaxy resampling "
      "unit the GD optimum reaches the full-sample dial; the "
      "dial-tension downgrades to gray and TODO 26b re-ranks"
      + cens + ".")
else:
    P("==> 8V VERDICT (locked grammar): GRAY-CARRIED - "
      f"P(D>=0) = {pge:.3f} sits between the locked bars; "
      "percentiles quoted as the measurement" + cens + ".")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t00)/60:.1f} min)")

with open('data/stage8v_gdboot.txt', 'w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage8v_gdboot.txt")
