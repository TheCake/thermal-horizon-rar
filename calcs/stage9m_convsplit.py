"""STAGE 9M — THE CONVERGENCE SPLIT on the GD dial (round-13 A4).
Pre-registered BEFORE any run (measurement round; NO credence
movement).

Frozen flags from the outer slope s_out = (V_N - V_{N-2})/V_N over
the last 3 KEPT rotmod points: RISING >= 0.05, CONVERGED <= 0.02,
AMBIG else (excluded, counted; < 3 kept points => AMBIG).
lam_hat(GD-conv) vs lam_hat(GD-rising) point fits + paired galaxy
bootstrap NBOOT=200 rng 71; DD control point fits.  SMALL-N: point
subset < 8 => DESCRIPTIVE-ONLY.  Bars (ordered): M-ARTIFACT /
M-HARDENED / M-GRAY-CARRIED.
Output: data/stage9m_convsplit.txt
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
vkept = {}
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
        vkept.setdefault(gi, []).append(Vo)
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
P("9M THE CONVERGENCE SPLIT ON THE GD DIAL (round-13 A4; pre-reg "
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

LGB = np.round(np.arange(-2.0, 1.501, 0.25), 3)
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
P(f"G9M-0 verbatim-lift: GD {len(gd_set)} / DD {len(dd_set)} -> "
  f"{'PASS' if ok_n else 'FAIL'}; OFF probe d = {va-OFFV:.2e} -> "
  f"{'PASS' if ok_off else 'FAIL'}")

rngp = np.random.default_rng(7)
g1_ok = True; dmx = 0.0
probe_sets = [[int(g_) for g_ in gd_set],
              [int(g_) for g_ in dd_set],
              [int(g_) for g_ in allg]]
for pi in range(3):
    thp = [-10.5+rngp.random(), 0.4+2.0*rngp.random(),
           0.01+0.29*rngp.random()]
    lam = -2.0+3.5*rngp.random()
    gl = probe_sets[pi]
    slow = m2ll_vert(thp, lam, gl, von=True)
    cat, lab, svec, ninst = make_instance(gl)
    fast = m2ll_fast(thp, lam, cat, lab, svec, ninst)
    if slow < 1e11:
        dmx = max(dmx, abs(slow-fast))
g1_ok = dmx <= 1e-6
P(f"G9M-1 fast-vs-verbatim (3 probes): max|d| = {dmx:.2e} -> "
  f"{'PASS' if g1_ok else 'FAIL'}")

# ---------------- frozen flags ----------------
flags = {}
for g_ in gd_set:
    V = vkept.get(int(g_), [])
    if len(V) < 3:
        flags[int(g_)] = 'AMBIG'; continue
    s_out = (V[-1] - V[-3])/max(V[-1], 1e-6)
    flags[int(g_)] = ('RISING' if s_out >= 0.05 else
                      'CONVERGED' if s_out <= 0.02 else 'AMBIG')
conv = [g_ for g_, f_ in flags.items() if f_ == 'CONVERGED']
rise = [g_ for g_, f_ in flags.items() if f_ == 'RISING']
namb = sum(1 for f_ in flags.values() if f_ == 'AMBIG')
P(f"G9M-2 flag census (GD-38): CONVERGED {len(conv)}, RISING "
  f"{len(rise)}, AMBIG {namb} (excluded)")
descriptive = min(len(conv), len(rise)) < 8

if not (g0_ok and g1_ok):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    with open('data/stage9m_convsplit.txt', 'w') as f:
        f.write("\n".join(L)+"\n")
    raise SystemExit(0)
P("GATES: G9M-0/1 PASS; census printed")
P("")

lh_c, ec1, ec2 = lam_hat_fast(*make_instance(conv))
lh_r, er1, er2 = lam_hat_fast(*make_instance(rise))
Dpoint = lh_r - lh_c
P(f"[point] lam_hat(CONVERGED) = {lh_c:+.3f} (edges {int(ec1)}/"
  f"{int(ec2)}), lam_hat(RISING) = {lh_r:+.3f} (edges {int(er1)}/"
  f"{int(er2)}); D = {Dpoint:+.3f}"
  + ("  DESCRIPTIVE-ONLY (subset < 8)" if descriptive else ""))

NBOOT = 200
rng = np.random.default_rng(71)
Ds = []; skips = 0
cvals = []
gd_list = [int(g_) for g_ in gd_set]
for r in range(NBOOT):
    draw = rng.integers(0, len(gd_list), size=len(gd_list))
    rep = [gd_list[j] for j in draw]
    if r == 0:
        P(f"G9M-3 first-replicate fingerprint: sum(draw) = "
          f"{int(np.sum(draw))}, unique = {len(np.unique(draw))}")
    c_ = [g_ for g_ in rep if flags[g_] == 'CONVERGED']
    r_ = [g_ for g_ in rep if flags[g_] == 'RISING']
    if min(len(c_), len(r_)) < 5:
        skips += 1; continue
    lc, _, _ = lam_hat_fast(*make_instance(c_))
    lr, _, _ = lam_hat_fast(*make_instance(r_))
    Ds.append(lr - lc); cvals.append(lc)
    if (r+1) % 50 == 0:
        P(f"  ... {r+1}/{NBOOT} ({(time.time()-t00)/60:.1f} min)")
Ds = np.array(Ds); cvals = np.array(cvals)
qd = np.percentile(Ds, [5, 50, 95])
qc = np.percentile(cvals, [5, 50, 95])
pge = float(np.mean(Ds >= 0))
P(f"[boot] reps {len(Ds)} (skips {skips}); D pct 5/50/95 = "
  f"{qd[0]:+.3f}/{qd[1]:+.3f}/{qd[2]:+.3f}; P(D >= 0) = {pge:.3f}")
P(f"[boot] lam(CONV) pct 5/50/95 = {qc[0]:+.3f}/{qc[1]:+.3f}/"
  f"{qc[2]:+.3f}")

flagsD = {}
for g_ in dd_set:
    V = vkept.get(int(g_), [])
    if len(V) < 3:
        flagsD[int(g_)] = 'AMBIG'; continue
    s_out = (V[-1] - V[-3])/max(V[-1], 1e-6)
    flagsD[int(g_)] = ('RISING' if s_out >= 0.05 else
                       'CONVERGED' if s_out <= 0.02 else 'AMBIG')
convD = [g_ for g_, f_ in flagsD.items() if f_ == 'CONVERGED']
riseD = [g_ for g_, f_ in flagsD.items() if f_ == 'RISING']
P("")
P(f"[DD control] CONVERGED {len(convD)}, RISING {len(riseD)}")
lhcD, _, _ = lam_hat_fast(*make_instance(convD))
lhrD, _, _ = lam_hat_fast(*make_instance(riseD))
P(f"[DD control] lam_hat(CONV) = {lhcD:+.3f}, lam_hat(RISING) = "
  f"{lhrD:+.3f}; D = {lhrD-lhcD:+.3f}")
P("")

m_art = ((not descriptive) and (Dpoint <= -0.5) and (pge <= 0.05)
         and (lh_c >= -0.6))
m_hard = ((not descriptive) and (lh_c <= -0.8)
          and (float(np.percentile(cvals, 95)) <= -0.3))
if m_art:
    P("==> 9M VERDICT (locked grammar): M-ARTIFACT - the negative "
      "dial concentrates in rising-curve galaxies and the "
      "converged half releases toward the dial: the V-ordering "
      "reads as a V_flat-definition artifact; the GD tension "
      "re-ranks.")
elif m_hard:
    P("==> 9M VERDICT (locked grammar): M-HARDENED - the "
      "converged GD subset alone carries the dial robustly; the "
      "EIGHTH control passes and the physics/selection reading "
      "hardens; external resolved kinematics = the clean arbiter.")
else:
    P("==> 9M VERDICT (locked grammar): M-GRAY-CARRIED - rows "
      "stand as measurements.")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t00)/60:.1f} min)")

with open('data/stage9m_convsplit.txt', 'w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage9m_convsplit.txt")
