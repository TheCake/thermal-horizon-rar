"""STAGE 8Y — THE PRESSURE-SUPPORT DOSE TEST (the GD tension).
Pre-registered BEFORE any run.

Pressure-support (asymmetric-drift) bias scales as (sigma/V)^2: if
it drives the GD dial, the slow rotators must carry it.
(a) V-split: GD halved at the identity-sample median of V (Vflat
from the SPARC table; fallback max V_obs where Vflat = 0); paired
bootstrap 300 reps (rng 71) over the 38 GD galaxies, threshold
FROZEN, halves fit when >= 6 galaxies; D_fs = lam_fast - lam_slow.
(b) The correction lever: refit lam_GD with g_obs -> (V^2 +
k*sigma^2)/R at sigma = 10 km/s, k in {0, 0.5, 1, 1.5, 2, 3}
(flat-sigma bound; real AD corrections are radius-dependent — this
brackets the magnitude); DD at k in {0, 2} = the bluntness control.

Gates: G8Y-0 the k=0 GD fit reproduces 8X's printed lam (bar
0.002); G8Y-1 fast-vs-slow engine cross-probes <= 1e-6 (the
lg-argument variant of the gated 8V/8X engine); G8Y-2 V-metric
accounting printed.
Bars (locked): Y1 DIRECTION iff P(D_fs <= 0) <= 0.05; sub-clause
Y1a REJOINS iff p95(lam_fast) >= 0.960 (the 8S-c FULL dial).  Y2
FLAT iff P(D_fs <= 0) >= 0.20 AND |median D_fs| <= 0.5.  Y3
GRAY-CARRIED else.  k-curve = co-read, pre-stated reading: lam_GD
crossing 0 by k <= 2 => plausible-magnitude corrections neutralize
the dial; no crossing by k = 3 => out of reach.
NO credence movement (measurement round; pre-stated).
Output: data/stage8y_pressure.txt
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
        vflat = float(t[15])
        q = int(t[17])
        meta[name] = (inc, q, D, eD, einc, vflat)
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
vobs2, rkpc = [], []
svert = {}
VFLAT = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc, vfl = meta.get(name, (0, 3, 1.0, 0.0, 0.0, 0.0))
    if inc < 30 or q > 2: continue
    s_d2 = (eD/max(D, 1e-6))/LN10
    s_i2 = 2.0*math.radians(einc)/max(math.tan(math.radians(inc)), 1e-6)/LN10
    svert[gi] = math.sqrt(s_d2**2 + s_i2**2)
    VFLAT[gi] = vfl
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
        vobs2.append(Vo*Vo); rkpc.append(R)
        sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id, vobs2, rkpc = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id, vobs2, rkpc))
sig2 = sig*sig
LGOBS0 = np.log10(gobs)

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
P("8Y THE PRESSURE-SUPPORT DOSE TEST (pre-reg committed BEFORE any "
  "run; measurement round; NO credence movement)")

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
GD_LIST = [int(g_) for g_ in gd_set]
DD_LIST = [int(g_) for g_ in dd_set]

# V metric: Vflat, fallback max Vobs
VMET = {}
n_vf, n_fb = 0, 0
for g_ in allg:
    vf = VFLAT.get(int(g_), 0.0)
    if vf > 0:
        VMET[int(g_)] = vf; n_vf += 1
    else:
        VMET[int(g_)] = float(np.sqrt(vobs2[GIDX[g_]].max())); n_fb += 1
VTHR = float(np.median([VMET[g] for g in GD_LIST]))
gd_slow = [g for g in GD_LIST if VMET[g] <= VTHR]
gd_fast = [g for g in GD_LIST if VMET[g] > VTHR]
P(f"G8Y-2 V metric: Vflat available {n_vf}, fallback max-Vobs "
  f"{n_fb}; GD threshold (frozen) = {VTHR:.1f} km/s; "
  f"slow/fast = {len(gd_slow)}/{len(gd_fast)}; GD V range "
  f"{min(VMET[g] for g in GD_LIST):.0f}-"
  f"{max(VMET[g] for g in GD_LIST):.0f}")

txt8x = open('data/stage8x_regime.txt').read()
TGT_GD = float(re.search(r"\[GD-all\s*\] lam_hat = ([+-][\d.]+)",
                         txt8x).group(1))
DIAL = 0.960   # 8S-c FULL vertON

def m2ll_slow(th, lam, gset, lg):
    la0, f, s_int = th
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    nv = nu_lam(gN/a0, lam)
    if np.min(nv) <= 0.05: return 1e12
    gm = gN*nv
    v = sig2 + s_int*s_int
    a = lg - np.log10(gm)
    out = 0.0
    for g_ in gset:
        ji = GIDX[g_]
        if len(ji) == 0: continue
        aj, vj = a[ji], v[ji]
        s_ = SV[g_]
        if s_ <= 1e-4:
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
    if li == 0:
        return None
    return (np.concatenate(cats), np.concatenate(labs),
            np.array(svals), li)

def m2ll_fast(th, lam, cat, lab, svec, ninst, lg):
    la0, f, s_int = th
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    nv = nu_lam(gN/a0, lam)
    if np.min(nv) <= 0.05: return 1e12
    gm = gN*nv
    v = sig2 + s_int*s_int
    a = lg - np.log10(gm)
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

LGB = np.round(np.arange(-2.0, 1.501, 0.25), 3)
def lam_hat_fast(inst, lg):
    cat, lab, svec, ninst = inst
    prof = []; th = None
    for lam in LGB:
        starts = (([list(th)] if th is not None else [])
                  + [[math.log10(A0_FID), 1.0, 0.08]])
        best = None
        for th0 in starts:
            b = minimize(lambda t: m2ll_fast(t, lam, cat, lab, svec,
                                             ninst, lg), th0,
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

def lg_of_k(k):
    return np.log10((vobs2 + 100.0*k)/rkpc*KPC)

# ---------------- gates ----------------
inst_gd = make_instance(GD_LIST)
lh0, _, _ = lam_hat_fast(inst_gd, LGOBS0)
ok0 = abs(lh0 - TGT_GD) <= 0.002
P(f"G8Y-0 k=0 GD fit: lam_hat = {lh0:+.3f} vs 8X {TGT_GD:+.3f} -> "
  f"{'PASS' if ok0 else 'FAIL'}")
lgk0 = lg_of_k(0.0)
dk0 = float(np.max(np.abs(lgk0 - LGOBS0)))
P(f"G8Y-0b lg_of_k(0) vs LGOBS0: max|d| = {dk0:.2e} -> "
  f"{'PASS' if dk0 <= 1e-12 else 'FAIL'}")
rngp = np.random.default_rng(7)
dmx = 0.0
psets = [GD_LIST, gd_slow, DD_LIST]
for pi in range(12):
    thp = [-10.5+rngp.random(), 0.4+2.0*rngp.random(),
           0.01+0.29*rngp.random()]
    lam = -2.0+3.5*rngp.random()
    kk = [0.0, 1.0, 2.0][pi % 3]
    gl = psets[pi % 3]
    lg = lg_of_k(kk)
    slow = m2ll_slow(thp, lam, gl, lg)
    inst = make_instance(gl)
    fast = m2ll_fast(thp, lam, *inst, lg)
    if slow < 1e11:
        dmx = max(dmx, abs(slow-fast))
ok1 = dmx <= 1e-6
P(f"G8Y-1 fast-vs-slow (12 probes, k in 0/1/2): max|d| = {dmx:.2e} "
  f"-> {'PASS' if ok1 else 'FAIL'}")
if not (ok0 and dk0 <= 1e-12 and ok1):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    with open('data/stage8y_pressure.txt', 'w') as f:
        f.write("\n".join(L)+"\n")
    raise SystemExit(0)
P("GATES: G8Y-0/0b/1/2 ALL PASS")
P("")

# ---------------- identity halves + k-curve ----------------
for tag, gl in (('GD-slow', gd_slow), ('GD-fast', gd_fast)):
    lh, lo_e, hi_e = lam_hat_fast(make_instance(gl), LGOBS0)
    ed = ' LO-EDGE' if lo_e else (' HI-EDGE' if hi_e else '')
    P(f"[{tag:7}] lam_hat = {lh:+.3f} (V "
      f"{'<=' if tag == 'GD-slow' else '>'} {VTHR:.1f}){ed}")
KGRID = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]
kcur = []
for k in KGRID:
    lh, lo_e, hi_e = lam_hat_fast(inst_gd, lg_of_k(k))
    kcur.append(lh)
    P(f"[GD k-curve] k={k:.1f} (sigma_eff={10*math.sqrt(k) if k > 0 else 0:.1f} "
      f"km/s): lam_hat = {lh:+.3f}"
      + (' LO-EDGE' if lo_e else (' HI-EDGE' if hi_e else '')))
inst_dd = make_instance(DD_LIST)
for k in (0.0, 2.0):
    lh, _, _ = lam_hat_fast(inst_dd, lg_of_k(k))
    P(f"[DD control] k={k:.1f}: lam_hat = {lh:+.3f}")
kcur = np.array(kcur)
for thr_, nm_ in ((-0.5, 'lam=-0.5'), (0.0, 'lam=0'), (DIAL, 'dial')):
    if kcur[0] >= thr_:
        P(f"k*({nm_}): already at k=0")
    elif np.all(kcur < thr_):
        P(f"k*({nm_}): NOT REACHED by k=3")
    else:
        j = int(np.argmax(kcur >= thr_))
        kstar = np.interp(thr_, [kcur[j-1], kcur[j]],
                          [KGRID[j-1], KGRID[j]])
        P(f"k*({nm_}) = {kstar:.2f} (sigma_eff = "
          f"{10*math.sqrt(kstar):.1f} km/s)")
P("")

# ---------------- the paired V-split bootstrap ----------------
NBOOT = 300
rng = np.random.default_rng(71)
fs, sl = [], []
skips = 0
for r in range(NBOOT):
    draw = rng.integers(0, len(GD_LIST), size=len(GD_LIST))
    rep = [GD_LIST[j] for j in draw]
    r_sl = [g for g in rep if VMET[g] <= VTHR]
    r_fs = [g for g in rep if VMET[g] > VTHR]
    if len(r_sl) < 6 or len(r_fs) < 6:
        fs.append(np.nan); sl.append(np.nan); skips += 1; continue
    lh_s, _, _ = lam_hat_fast(make_instance(r_sl), LGOBS0)
    lh_f, _, _ = lam_hat_fast(make_instance(r_fs), LGOBS0)
    sl.append(lh_s); fs.append(lh_f)
    if (r+1) % 100 == 0:
        P(f"  ... {r+1}/{NBOOT} replicates "
          f"({(time.time()-t00)/60:.1f} min)")
fs = np.array(fs); sl = np.array(sl)
val = np.isfinite(fs) & np.isfinite(sl)
dfs = fs[val] - sl[val]
p_le0 = float(np.mean(dfs <= 0))
qf = np.percentile(fs[val], [5, 50, 95])
qs = np.percentile(sl[val], [5, 50, 95])
qd = np.percentile(dfs, [5, 50, 95])
P("")
P(f"replicates {NBOOT}, skips {skips}")
P(f"[GD-fast] pct 5/50/95 = {qf[0]:+.3f}/{qf[1]:+.3f}/{qf[2]:+.3f}")
P(f"[GD-slow] pct 5/50/95 = {qs[0]:+.3f}/{qs[1]:+.3f}/{qs[2]:+.3f}")
P(f"[D_fs = fast-slow] pct 5/50/95 = {qd[0]:+.3f}/{qd[1]:+.3f}/"
  f"{qd[2]:+.3f}; P(D_fs <= 0) = {p_le0:.4f}")
P("")
y1 = p_le0 <= 0.05
y1a = y1 and (qf[2] >= DIAL)
y2 = (p_le0 >= 0.20) and (abs(qd[1]) <= 0.5)
if y1:
    P("==> 8Y VERDICT (locked grammar): DIRECTION - the slow "
      "rotators carry the negative dial (pressure-support "
      "direction confirmed)"
      + ("; Y1a REJOINS - the fast half reaches the dial"
         if y1a else "; fast half does NOT reach the dial (partial)")
      + ".")
elif y2:
    P("==> 8Y VERDICT (locked grammar): FLAT - no V-dependence; "
      "pressure support disfavored as the carrier; genuine "
      "composition shape strengthens.")
else:
    P("==> 8Y VERDICT (locked grammar): GRAY-CARRIED - "
      "percentiles quoted as the measurement.")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t00)/60:.1f} min)")

with open('data/stage8y_pressure.txt', 'w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage8y_pressure.txt")
