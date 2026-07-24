"""
STAGE 6Z (the 6Y sky check): the gate-assignment ORDERING test.

The reservoir identification predicts the per-galaxy gate ORDERING is
physics: galaxies in weaker ambients carry sharper tails (p monotone
decreasing in e_N). 6I measured that substituting the per-galaxy Chae
gates for the global one IMPROVES the vertical treatment by ~2.6
(-61.68 vs -59.05) -- but improvement alone does not isolate ordering.
The clean control: SHUFFLE the matched galaxies' gates among
themselves. A permutation preserves the gate DISTRIBUTION exactly
(gates are inputs, not fitted -- zero flexibility difference), so any
advantage of the TRUE assignment over the shuffle null is pure ordering
information = the identification's signature in the sky.

PRE-REGISTERED BARS (committed before execution; K = 8 shuffles, rng
71, vertical-hardened treatment, warm-started from the true fit):
  SUPPORT if the true assignment beats ALL K shuffles (empirical
          p <= 1/(K+1) ~ 0.11);
  LEAN    if it beats >= 6/8;
  NULL    if in the bulk (ordering unresolved at 149-galaxy grade --
          the expected outcome if the ~2.6-point signal is marginal);
  ANTI    if it loses to ALL K (a strike against the gate ordering).
Regression gates: BE vs 5P (-12152.49, d = -0.00 grade); TRUE delta vs
the 6I PGmax record (-61.68, |dd| < 2).

Writes data/stage6z_ordering.txt.
"""
import csv, glob, math, os, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25
LN10 = math.log(10)
S_ML = 0.1*LN10

def s_of(e):
    nn = 1.0/(math.exp(math.sqrt(e)) - 1.0)
    return nn/(1.0 + nn)

meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta[t[0]] = (float(t[5]), int(t[17]), float(t[2]), float(t[3]),
                      float(t[6]))
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
sigv_g_map = {}
gal_name = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 10.0, 1.0, 3.0))
    if inc < 30 or q > 2: continue
    kept += 1
    gal_name[gi] = name
    irad = math.radians(inc)
    sv = math.hypot((eD/max(D, 1e-3))/LN10,
                    2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10)
    sigv_g_map[gi] = max(sv, 0.01)
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
ug = np.unique(gal_id)
NGal = len(ug)
gmap = {g: i for i, g in enumerate(ug)}
gidx = np.array([gmap[g] for g in gal_id])
SIGV = np.array([sigv_g_map[g] for g in ug])
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

chae = {}
with open('data/chae2021_table3.csv') as f:
    for row in csv.DictReader(l for l in f if not l.startswith('#')):
        chae[row['galaxy']] = 10.0**float(row['log_eN_maxclust'])
s_raw = np.full(NGal, np.nan)
for i in range(NGal):
    nm = gal_name[ug[i]]
    if nm in chae:
        s_raw[i] = s_of(chae[nm])
matched = np.isfinite(s_raw)
NM = int(matched.sum())
S_MED = float(np.median(s_raw[matched]))
S_TRUE = np.where(matched, s_raw, S_MED)

def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))

def nu_pg(y, spt):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    spt = np.broadcast_to(np.asarray(spt, float), y.shape)
    sL = spt**2
    ly = np.log(y)
    nu = nu_simple(y)
    for _ in range(80):
        d1 = 2.0*nu - 1.0
        b = 0.5*sL/d1**2
        db = -2.0*sL/d1**3
        u = np.exp(np.minimum(0.5*(1.0+b)*ly + b*np.log(nu), 60.0))
        eu = np.exp(np.minimum(u, 60.0))
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        F = nu - 1.0 - n
        dudnu = u*(db*(0.5*ly + np.log(nu)) + b/nu)
        dF = 1.0 + (eu/(em1*em1))*dudnu
        step = F/np.where(np.abs(dF) > 1e-12, dF, 1e-12)
        nu = np.maximum(nu - np.clip(step, -0.5*nu, 0.5*nu), 1.0 + 1e-15)
    return nu

def m2hv(th, dml, dv, spt, sl):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*nu_pg(gN/a0, spt)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm_) - dv[gidx]
    out = np.sum(r*r/se2 + np.log(se2))
    lg = l_gbar[lmask] + dlt
    rl = l_gobs[lmask] - (lg + np.log10(nu_pg(10**lg/a0, sl)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(dml*dml)/(S_ML*S_ML)
    out += np.sum(dv*dv/(SIGV*SIGV))
    return out

def fit_conv(spt_g, sl, th0=None, dml0=None, dv0=None, tol=0.05,
             max_rounds=15):
    spt = spt_g[gidx]
    dml = np.zeros(NGal) if dml0 is None else dml0.copy()
    dv = np.zeros(NGal) if dv0 is None else dv0.copy()
    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08, 0.0]]
                  if rd == 0 and th0 is None else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2hv(t, dml, dv, spt, sl),
                         t0, method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int, dlt = best.x
        se2c = s_int*s_int
        for _ in range(3):
            fac = f*np.exp(dml[gidx])
            gN = g_gas + fac*g_dsk + g_bul
            r0_ = lgobs - np.log10(gN*nu_pg(gN/10**la0, spt))
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                w = 1.0/(sig2[mm] + se2c)
                dv[gi2] = np.sum(w*r0_[mm])/(np.sum(w) + 1.0/SIGV[gi2]**2)
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = (lgobs[mm]
                          - np.log10(gN2*nu_pg(gN2/10**la0, spt[mm]))
                          - dv[gi2])
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2hv(best.x, dml, dv, spt, sl)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2hv(t, dml, dv, spt, sl),
                 list(best.x), method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml, dv

REF_BE = -12152.49
REF_TRUE_DELTA = -61.68        # 6I PGmax vertical delta

L = [f"STAGE 6Z: the gate-assignment ordering test -- {kept} galaxies "
     f"({NM} Chae-matched), vertical-hardened, K = 8 shuffles (rng 71)",
     "pre-registered: SUPPORT beats 8/8; LEAN >= 6/8; NULL bulk; "
     "ANTI 0/8",
     ""]
def save():
    with open('data/stage6z_ordering.txt', 'w') as f:
        f.write("\n".join(L) + "\n")
save()

t0 = time.time()
bb, dm, dvv = fit_conv(np.zeros(NGal), 0.0)
d = bb.fun - REF_BE
L.append(f"  BE: {bb.fun:10.2f}  [5P {REF_BE:.2f}, d={d:+.2f} "
         f"{'OK' if abs(d) < 2.0 else 'FAIL'}]  "
         f"({(time.time()-t0)/60:.1f} min)")
print(L[-1], flush=True)
save()
assert abs(d) < 2.0
OBJ_BE = bb.fun
th = bb.x

t0 = time.time()
bb, dm, dvv = fit_conv(S_TRUE, S_MED, th0=th, dml0=dm, dv0=dvv)
d_true = bb.fun - OBJ_BE
dd = d_true - REF_TRUE_DELTA
L.append(f"  TRUE: {bb.fun:10.2f}  delta {d_true:+.2f}  [6I {REF_TRUE_DELTA:+.2f}, "
         f"dd={dd:+.2f} {'OK' if abs(dd) < 2.0 else 'FAIL'}]  "
         f"({(time.time()-t0)/60:.1f} min)")
print(L[-1], flush=True)
save()
assert abs(dd) < 2.0
th_true, dm_true, dv_true = bb.x, dm.copy(), dvv.copy()

rng = np.random.default_rng(71)
idx_m = np.where(matched)[0]
d_shuf = []
for k in range(8):
    perm = rng.permutation(len(idx_m))
    S_k = S_TRUE.copy()
    S_k[idx_m] = s_raw[idx_m][perm]
    t0 = time.time()
    bb, _, _ = fit_conv(S_k, S_MED, th0=th_true, dml0=dm_true,
                        dv0=dv_true)
    dk = bb.fun - OBJ_BE
    d_shuf.append(dk)
    L.append(f"  shuffle {k}: delta {dk:+.2f}  "
             f"({(time.time()-t0)/60:.1f} min)")
    print(L[-1], flush=True)
    save()

d_shuf = np.array(d_shuf)
nbeat = int((d_true < d_shuf).sum())
L.append("")
L.append(f"true delta {d_true:+.2f} vs shuffle null: mean "
         f"{d_shuf.mean():+.2f} +- {d_shuf.std(ddof=1):.2f}, "
         f"range [{d_shuf.min():+.2f}, {d_shuf.max():+.2f}]; true beats "
         f"{nbeat}/8")
if nbeat == 8: v = "SUPPORT (ordering signal; empirical p <= 0.11)"
elif nbeat >= 6: v = "LEAN (ordering favored, not resolved)"
elif nbeat == 0: v = "ANTI (strike against the gate ordering)"
else: v = "NULL (ordering unresolved at this sample grade)"
L.append(f"VERDICT vs pre-registered bars: {v}")
print("\n".join(L[-2:]), flush=True)
save()
print("\nSTAGE 6Z done -> data/stage6z_ordering.txt")
