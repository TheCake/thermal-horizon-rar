"""STAGE 9B — THE QUALITY-FLAG CONTROL (GD).  Pre-registered
BEFORE any run.

Is the GD dial carried by the SPARC quality flag?  GD split Q=1
(high) vs Q=2 (medium); identity fits + paired bootstrap 300 reps
(rng 71); D_Q = lam(Q1) - lam(Q2).  Co-read: the 8Y V-split
repeated within each Q stratum (identity fits only,
power-limited).  Sizes printed; < 8 galaxies either side =>
POWER-FLAG (pre-stated).

Gates: G9B-0 GD-all reproduces 8X's lam (bar 0.002); G9B-1
accounting incl. the 8Y V-threshold regression (66.2).
Bars (locked): B1 QUALITY-CARRIED iff P(D_Q <= 0) <= 0.05 AND
median D_Q >= 1.0.  B2 QUALITY-BLIND iff P(D_Q <= 0) >= 0.20 AND
|median D_Q| <= 0.5.  B3 GRAY-CARRIED else.
NO credence movement (measurement round; pre-stated).
Output: data/stage9b_qflag.txt
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
vobs2 = []
svert = {}
VFLAT = {}
QFLAG = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc, vfl = meta.get(name, (0, 3, 1.0, 0.0, 0.0, 0.0))
    if inc < 30 or q > 2: continue
    s_d2 = (eD/max(D, 1e-6))/LN10
    s_i2 = 2.0*math.radians(einc)/max(math.tan(math.radians(inc)), 1e-6)/LN10
    svert[gi] = math.sqrt(s_d2**2 + s_i2**2)
    VFLAT[gi] = vfl
    QFLAG[gi] = q
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
        vobs2.append(Vo*Vo)
        sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id, vobs2 = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id, vobs2))
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
P("9B THE QUALITY-FLAG CONTROL (pre-reg committed BEFORE any run; "
  "measurement round; NO credence movement)")

allg = np.unique(gal_id)
NGAL = gal_id.max()+1
gdfrac = np.zeros(NGAL)
for g_ in allg:
    m = gal_id == g_
    gdfrac[g_] = float(np.mean(g_gas[m] > g_dsk[m] + g_bul[m]))
gd_set = np.array([g_ for g_ in allg if gdfrac[g_] >= 0.5])
SV = np.array([svert.get(g_, 0.0) for g_ in range(NGAL)])
GIDX = [np.where(gal_id == g_)[0] for g_ in range(NGAL)]
GD_LIST = [int(g_) for g_ in gd_set]

VMET = {}
for g_ in allg:
    vf = VFLAT.get(int(g_), 0.0)
    VMET[int(g_)] = vf if vf > 0 else \
        float(np.sqrt(vobs2[GIDX[g_]].max()))
VTHR = float(np.median([VMET[g] for g in GD_LIST]))
GDQ1 = [g for g in GD_LIST if QFLAG[g] == 1]
GDQ2 = [g for g in GD_LIST if QFLAG[g] == 2]
pflag = (len(GDQ1) < 8) or (len(GDQ2) < 8)
P(f"G9B-1 accounting: GD Q1/Q2 = {len(GDQ1)}/{len(GDQ2)}"
  + ("  POWER-FLAG" if pflag else "")
  + f"; V threshold = {VTHR:.1f} (8Y: 66.2) -> "
  f"{'PASS' if abs(VTHR-66.2) <= 0.1 else 'FAIL'}")

txt8x = open('data/stage8x_regime.txt').read()
TGT_GD = float(re.search(r"\[GD-all\s*\] lam_hat = ([+-][\d.]+)",
                         txt8x).group(1))

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

LGB = np.round(np.arange(-2.0, 1.501, 0.25), 3)
def lam_hat_fast(inst):
    cat, lab, svec, ninst = inst
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

lh0, _, _ = lam_hat_fast(make_instance(GD_LIST))
ok0 = abs(lh0 - TGT_GD) <= 0.002
P(f"G9B-0 GD-all: lam_hat = {lh0:+.3f} vs 8X {TGT_GD:+.3f} -> "
  f"{'PASS' if ok0 else 'FAIL'}")
if not (ok0 and abs(VTHR-66.2) <= 0.1):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    with open('data/stage9b_qflag.txt', 'w') as f:
        f.write("\n".join(L)+"\n")
    raise SystemExit(0)
P("GATES: G9B-0/1 ALL PASS")
P("")

for tag, gl in (('GD-Q1', GDQ1), ('GD-Q2', GDQ2)):
    inst = make_instance(gl)
    if inst is None:
        P(f"[{tag}] EMPTY"); continue
    lh, lo_e, hi_e = lam_hat_fast(inst)
    ed = ' LO-EDGE' if lo_e else (' HI-EDGE' if hi_e else '')
    P(f"[{tag}] lam_hat = {lh:+.3f} (n = {len(gl)}){ed}")
for qt, gl in (('Q1', GDQ1), ('Q2', GDQ2)):
    for vt, sub in (('slow', [g for g in gl if VMET[g] <= VTHR]),
                    ('fast', [g for g in gl if VMET[g] > VTHR])):
        inst = make_instance(sub)
        if inst is None or len(sub) < 4:
            P(f"[co-read {qt}-{vt}] n = {len(sub)} (too thin)")
            continue
        lh, lo_e, hi_e = lam_hat_fast(inst)
        ed = ' LO-EDGE' if lo_e else (' HI-EDGE' if hi_e else '')
        P(f"[co-read {qt}-{vt}] lam_hat = {lh:+.3f} "
          f"(n = {len(sub)}){ed}")
P("")

NBOOT = 300
rng = np.random.default_rng(71)
q1v, q2v = [], []
skips = 0
for r in range(NBOOT):
    draw = rng.integers(0, len(GD_LIST), size=len(GD_LIST))
    rep = [GD_LIST[j] for j in draw]
    r1 = [g for g in rep if QFLAG[g] == 1]
    r2 = [g for g in rep if QFLAG[g] == 2]
    if len(r1) < 5 or len(r2) < 5:
        q1v.append(np.nan); q2v.append(np.nan); skips += 1; continue
    lh1, _, _ = lam_hat_fast(make_instance(r1))
    lh2, _, _ = lam_hat_fast(make_instance(r2))
    q1v.append(lh1); q2v.append(lh2)
    if (r+1) % 100 == 0:
        P(f"  ... {r+1}/{NBOOT} replicates "
          f"({(time.time()-t00)/60:.1f} min)")
q1v = np.array(q1v); q2v = np.array(q2v)
val = np.isfinite(q1v) & np.isfinite(q2v)
dq = q1v[val] - q2v[val]
p_le0 = float(np.mean(dq <= 0))
qq1 = np.percentile(q1v[val], [5, 50, 95])
qq2 = np.percentile(q2v[val], [5, 50, 95])
qqd = np.percentile(dq, [5, 50, 95])
P("")
P(f"replicates {NBOOT}, skips {skips}")
P(f"[GD-Q1] pct 5/50/95 = {qq1[0]:+.3f}/{qq1[1]:+.3f}/{qq1[2]:+.3f}")
P(f"[GD-Q2] pct 5/50/95 = {qq2[0]:+.3f}/{qq2[1]:+.3f}/{qq2[2]:+.3f}")
P(f"[D_Q = Q1-Q2] pct 5/50/95 = {qqd[0]:+.3f}/{qqd[1]:+.3f}/"
  f"{qqd[2]:+.3f}; P(D_Q <= 0) = {p_le0:.4f}")
P("")
b1 = (p_le0 <= 0.05) and (qqd[1] >= 1.0)
b2 = (p_le0 >= 0.20) and (abs(qqd[1]) <= 0.5)
pf = "  [POWER-FLAG]" if pflag else ""
if b1:
    P("==> 9B VERDICT (locked grammar): QUALITY-CARRIED - the "
      "high-quality GD subset carries materially less dial; the "
      "data-quality suspect strengthens" + pf + ".")
elif b2:
    P("==> 9B VERDICT (locked grammar): QUALITY-BLIND - the dial "
      "ignores the SPARC quality flag; the data-quality suspect "
      "weakens" + pf + ".")
else:
    P("==> 9B VERDICT (locked grammar): GRAY-CARRIED - "
      "percentiles quoted as the measurement" + pf + ".")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t00)/60:.1f} min)")

with open('data/stage9b_qflag.txt', 'w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage9b_qflag.txt")
