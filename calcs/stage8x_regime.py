"""STAGE 8X — THE REGIME-vs-COMPOSITION DECOMPOSITION (the GD
tension).  Pre-registered BEFORE any run.

Is the GD-DD dial split carried by COMPOSITION (gas-dominated
galaxies) or by REGIME (deep points)?  GD galaxies are
deep-dominated; DD spirals carry the transition/tail.  Decider:
DD's own DEEP points.  Deep := gN(f=1) < 1.2e-10 (fixed pre-fit;
cut on the model-independent baryonic gN, not on gobs).  Six
identity fits (GD, DD) x (all, deep, nondeep), vertical-ON, the 8V
fast engine with point masks; paired galaxy bootstrap NBOOT = 200
(rng 71) over the three deciding subsets (DD-deep, GD-deep,
DD-nondeep — pre-stated cost cut).

Gates (any FAIL => STOP; amendment pre-quote): G8X-0 unmasked
GD/FULL coarse lam-hat reproduce 8V's printed values (bar 0.002);
G8X-1 masked fast-vs-slow cross-probes <= 1e-6 (incl. a duplicated
list); G8X-2 mask accounting exact (deep + nondeep = all, per set).
Bars (locked): X1 REGIME iff p95(DD-deep) < 0 AND
P(DD-deep >= DD-nondeep) <= 0.05 — the split follows the deep
regime.  X2 COMPOSITION iff p5(DD-deep) > 0 AND
P(DD-deep <= GD-deep) <= 0.05 — DD's deep points vote with DD.
X3 GRAY-CARRIED else.  Co-reads: GD-nondeep point fit, counts,
s_int per subset, edge fractions.  NO credence movement
(measurement round; pre-stated).
Output: data/stage8x_regime.txt
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
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 1.0, 0.0, 0.0))
    if inc < 30 or q > 2: continue
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
P("8X THE REGIME-vs-COMPOSITION DECOMPOSITION (pre-reg committed "
  "BEFORE any run; measurement round; NO credence movement)")

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
GD_LIST = [int(g_) for g_ in gd_set]
DD_LIST = [int(g_) for g_ in dd_set]

gN1 = g_gas + g_dsk + g_bul
DEEP = gN1 < A0_FID
MASKS = {'all': np.ones(len(gN1), dtype=bool),
         'deep': DEEP, 'nondeep': ~DEEP}

txt8v = open('data/stage8v_gdboot.txt').read()
TGT = {m.group(1): float(m.group(2)) for m in
       re.finditer(r"G8V-2 identity replicate \[(FULL|GD|DD)\]: "
                   r"coarse lam_hat = ([-\d.]+)", txt8v)}
assert {'FULL', 'GD'} <= set(TGT), TGT

def m2ll_vert_m(th, lam, gset, mask):
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
        ji = ji[mask[ji]]
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

def make_instance(gset_list, mask):
    cats, labs, svals = [], [], []
    li = 0
    for g_ in gset_list:
        ji = GIDX[g_]
        ji = ji[mask[ji]]
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
    return lam_hat, lo_edge, hi_edge, th

# ---------------- gates ----------------
g2_ok = True
for tag, gl in (('GD', GD_LIST), ('DD', DD_LIST)):
    n_a = sum(int(MASKS['all'][GIDX[g_]].sum()) for g_ in gl)
    n_d = sum(int(MASKS['deep'][GIDX[g_]].sum()) for g_ in gl)
    n_n = sum(int(MASKS['nondeep'][GIDX[g_]].sum()) for g_ in gl)
    okc = (n_d + n_n == n_a)
    g2_ok &= okc
    gd_ = sum(1 for g_ in gl if MASKS['deep'][GIDX[g_]].any())
    gn_ = sum(1 for g_ in gl if MASKS['nondeep'][GIDX[g_]].any())
    P(f"G8X-2 [{tag}] points all/deep/nondeep = {n_a}/{n_d}/{n_n} "
      f"(galaxies with >=1 pt: {gd_}/{gn_}) -> "
      f"{'PASS' if okc else 'FAIL'}")

g0_ok = True
for tag, gl in (('FULL', ALL_LIST), ('GD', GD_LIST)):
    lh, _, _, _ = lam_hat_fast(make_instance(gl, MASKS['all']))
    okt = abs(lh - TGT[tag]) <= 0.002
    g0_ok &= okt
    P(f"G8X-0 unmasked [{tag}]: lam_hat = {lh:.3f} vs 8V "
      f"{TGT[tag]:.3f} -> {'PASS' if okt else 'FAIL'}")

rngp = np.random.default_rng(7)
dup = [ALL_LIST[j] for j in rngp.integers(0, len(ALL_LIST),
                                          size=len(ALL_LIST))]
probe_sets = [(DD_LIST, 'deep'), (GD_LIST, 'deep'),
              (dup, 'nondeep')]
dmx = 0.0
for pi in range(15):
    thp = [-10.5+rngp.random(), 0.4+2.0*rngp.random(),
           0.01+0.29*rngp.random()]
    lam = -2.0+3.5*rngp.random()
    gl, mk = probe_sets[pi % 3]
    slow = m2ll_vert_m(thp, lam, gl, MASKS[mk])
    inst = make_instance(gl, MASKS[mk])
    fast = m2ll_fast(thp, lam, *inst)
    if slow < 1e11:
        dmx = max(dmx, abs(slow-fast))
g1_ok = dmx <= 1e-6
P(f"G8X-1 masked fast-vs-slow (15 probes incl. duplicated list): "
  f"max|d| = {dmx:.2e} -> {'PASS' if g1_ok else 'FAIL'}")

if not (g0_ok and g1_ok and g2_ok):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    with open('data/stage8x_regime.txt', 'w') as f:
        f.write("\n".join(L)+"\n")
    raise SystemExit(0)
P("GATES: G8X-0/1/2 ALL PASS")
P("")

# ---------------- identity fits (six) ----------------
idfit = {}
for tag, gl in (('GD', GD_LIST), ('DD', DD_LIST)):
    for mk in ('all', 'deep', 'nondeep'):
        inst = make_instance(gl, MASKS[mk])
        if inst is None:
            # amendment 1 (wiring, pre-quote): empty subset guard —
            # GD-nondeep has ZERO points (GD is 100% deep; the
            # confounding is total, which is why DD-deep decides)
            idfit[(tag, mk)] = np.nan
            P(f"[{tag}-{mk:7}] EMPTY (0 points)")
            continue
        lh, lo_e, hi_e, th = lam_hat_fast(inst)
        idfit[(tag, mk)] = lh
        ed = ' LO-EDGE' if lo_e else (' HI-EDGE' if hi_e else '')
        P(f"[{tag}-{mk:7}] lam_hat = {lh:+.3f} -> c1 = {lh/2:+.3f}; "
          f"s_int = {th[2]:.3f}{ed}")
P("")

# ---------------- paired bootstrap ----------------
NBOOT = 200
rng = np.random.default_rng(71)
keys = [('DD', 'deep'), ('GD', 'deep'), ('DD', 'nondeep')]
res = {k: [] for k in keys}
edges_lo = {k: 0 for k in keys}; edges_hi = {k: 0 for k in keys}
skips = 0
for r in range(NBOOT):
    draw = rng.integers(0, len(ALL_LIST), size=len(ALL_LIST))
    rep = [ALL_LIST[j] for j in draw]
    gd_rep = [g for g in rep if gdfrac[g] >= 0.5]
    dd_rep = [g for g in rep if gdfrac[g] < 0.5]
    for (tag, mk) in keys:
        gl = gd_rep if tag == 'GD' else dd_rep
        inst = make_instance(gl, MASKS[mk])
        if inst is None or inst[3] < 8:
            res[(tag, mk)].append(np.nan); skips += 1; continue
        lh, lo_e, hi_e, _ = lam_hat_fast(inst)
        res[(tag, mk)].append(lh)
        edges_lo[(tag, mk)] += int(lo_e)
        edges_hi[(tag, mk)] += int(hi_e)
    if (r+1) % 50 == 0:
        P(f"  ... {r+1}/{NBOOT} replicates "
          f"({(time.time()-t00)/60:.1f} min)")

ddd = np.array(res[('DD', 'deep')])
gdd = np.array(res[('GD', 'deep')])
ddn = np.array(res[('DD', 'nondeep')])
P("")
P(f"replicates: {NBOOT}; skips {skips}")
for k, arr in (('DD-deep', ddd), ('GD-deep', gdd),
               ('DD-nondeep', ddn)):
    q = np.percentile(arr[np.isfinite(arr)], [5, 50, 95])
    kk = ('DD', 'deep') if k == 'DD-deep' else \
         (('GD', 'deep') if k == 'GD-deep' else ('DD', 'nondeep'))
    P(f"[{k:10}] pct 5/50/95 = {q[0]:+.3f}/{q[1]:+.3f}/{q[2]:+.3f}; "
      f"edge lo/hi = {edges_lo[kk]}/{edges_hi[kk]}")
v1 = np.isfinite(ddd) & np.isfinite(ddn)
v2 = np.isfinite(ddd) & np.isfinite(gdd)
p_dd_ge_ddn = float(np.mean(ddd[v1] >= ddn[v1]))
p_dd_le_gdd = float(np.mean(ddd[v2] <= gdd[v2]))
p95_ddd = float(np.percentile(ddd[np.isfinite(ddd)], 95))
p5_ddd = float(np.percentile(ddd[np.isfinite(ddd)], 5))
P(f"P(DD-deep >= DD-nondeep) = {p_dd_ge_ddn:.4f}; "
  f"P(DD-deep <= GD-deep) = {p_dd_le_gdd:.4f}")
P("")
x1 = (p95_ddd < 0.0) and (p_dd_ge_ddn <= 0.05)
x2 = (p5_ddd > 0.0) and (p_dd_le_gdd <= 0.05)
if x1:
    P("==> 8X VERDICT (locked grammar): REGIME - DD's own deep "
      "points vote negative and differ from DD's nondeep arm: the "
      "GD-DD split follows the DEEP REGIME, not galaxy composition "
      "(connects to 5T's deep-arm vote; the selection story "
      "weakens).")
elif x2:
    P("==> 8X VERDICT (locked grammar): COMPOSITION - DD's deep "
      "points vote with DD and differ from GD's deep points: the "
      "split is carried by galaxy type, not by regime; selection "
      "and genuine composition-shape remain the suspects.")
else:
    P("==> 8X VERDICT (locked grammar): GRAY-CARRIED - percentiles "
      "quoted as the measurement.")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t00)/60:.1f} min)")

with open('data/stage8x_regime.txt', 'w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage8x_regime.txt")
