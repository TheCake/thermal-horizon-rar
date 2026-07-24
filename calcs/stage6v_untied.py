"""
STAGE 6V (the 6U falsifier): the untied-exponent contest.

The 6U loop derivation REQUIRES the grammar's two exponents to be tied:
each vertex of the one dressing loop carries the local zero-point share
AND the ambient KMS ratio, so beta = (1/2) * q^L * s^L with a single L.
An untied form beta = (1/2) * q^L1 * s^L2 (L1 != L2) has no per-vertex
pairing. This stage contests the off-diagonal cells directly on the
vertical-hardened hier ladder at the fiducial gate (6I comparators:
tied L1 -52.76 / L2 -59.05 / L3 -54.52 vs BE).

Deep-series instrument split (exact consequence of the untied form):
deep q ~ x/2 so beta ~ (1/2) s^L2 (x/2)^L1 -- the Bernoulli break rung
is c_(L1+1) (reads L1); the tail has q -> 1 so p = 1/2 + s^L2/4 (reads
L2). The vertical fit is tail-dominated (6L: the deep arm of 153
galaxies cannot read the rung at population grade), so the derivation
PREDICTS: varying L1 at fixed L2 moves the fit much less than varying
L2 at fixed L1.

PRE-REGISTERED BARS (committed before execution; all point-grade, with
the 6L population caveat carried verbatim):
  B1 (tied survival): no untied cell beats the tied (2,2) cell by more
      than +5 lnL at point grade -- else STRIKE against the per-vertex
      pairing. SUPPORT if (2,2) remains the best cell overall.
  B2 (instrument split): mean |Delta vs (2,2)| over {(1,2),(3,2)}
      (L1-varied) < mean |Delta vs (2,2)| over {(2,1),(2,3)}
      (L2-varied) -- PASS = the derivation's exponent-instrument
      mapping confirmed at point grade.
Caveat carried: point preferences at this instrument deflate under the
40-rep bootstrap (correction #13); B1/B2 are structural-direction
tests, not population-grade measurements.

Cells run fresh: BE (regression gate d = -0.00 vs 5P), U22 (regression
vs 6I AMBg delta -59.05, |dd| < 2), U12, U32, U21, U23. Tied (1,1) and
(3,3) quoted from the 6I file (identical machinery).

Writes data/stage6v_untied.txt.
"""
import glob, math, os, time
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
    n = 1.0/(math.exp(math.sqrt(e)) - 1.0)
    return n/(1.0 + n)
S_GLOB = s_of(0.02)          # 0.8681 -- the 6E/6F/6I fiducial gate

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
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 10.0, 1.0, 3.0))
    if inc < 30 or q > 2: continue
    kept += 1
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

def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))

def nu_be(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    x = np.sqrt(y)
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))

def make_untied(s, L1, L2):
    sL = s**L2
    def nu_run(y):
        y = np.clip(np.asarray(y, float), 1e-14, None)
        ly = np.log(y)
        nu = nu_simple(y)
        for _ in range(80):
            d1 = 2.0*nu - 1.0
            b = 0.5*sL/d1**L1
            db = -L1*sL/d1**(L1 + 1)
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
    return nu_run

def m2hv(th, dml, dv, use_v, nu_fn):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*nu_fn(gN/a0)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm_) - (dv[gidx] if use_v else 0.0)
    out = np.sum(r*r/se2 + np.log(se2))
    lg = l_gbar[lmask] + dlt
    rl = l_gobs[lmask] - (lg + np.log10(nu_fn(10**lg/a0)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(dml*dml)/(S_ML*S_ML)
    if use_v: out += np.sum(dv*dv/(SIGV*SIGV))
    return out

def fit_conv(use_v, nu_fn, th0=None, dml0=None, dv0=None, tol=0.05,
             max_rounds=15):
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
            b = minimize(lambda t: m2hv(t, dml, dv, use_v, nu_fn),
                         t0, method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int, dlt = best.x
        se2c = s_int*s_int
        for _ in range(3):
            if use_v:
                fac = f*np.exp(dml[gidx])
                gN = g_gas + fac*g_dsk + g_bul
                r0_ = lgobs - np.log10(gN*nu_fn(gN/10**la0))
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
                          - np.log10(gN2*nu_fn(gN2/10**la0))
                          - (dv[gi2] if use_v else 0.0))
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2hv(best.x, dml, dv, use_v, nu_fn)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2hv(t, dml, dv, use_v, nu_fn),
                 list(best.x), method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml, dv

REF_BE = -12152.49
REF_U22_DELTA = -59.05          # 6I AMBg vertical delta (tied L=2)
CMP_TIED = {1: -52.76, 3: -54.52}   # 6I L1g/L3g vertical deltas

L = [f"STAGE 6V: the untied-exponent contest -- {kept} galaxies, "
     f"{len(gobs)} points + {int(lmask.sum())} lensing; vertical-hardened; "
     f"fiducial gate s = {S_GLOB:.4f}",
     "pre-registered: B1 no untied cell beats (2,2) by > +5 (point grade); "
     "B2 |d| over {(1,2),(3,2)} < |d| over {(2,1),(2,3)}; 6L population "
     "caveat carried",
     ""]

def save():
    with open('data/stage6v_untied.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

save()
CELLS = [('BE', None), ('U22', (2, 2)), ('U12', (1, 2)), ('U32', (3, 2)),
         ('U21', (2, 1)), ('U23', (2, 3))]
VALS = {}
th, dm, dvv = None, None, None
for name, cell in CELLS:
    nu_fn = nu_be if cell is None else make_untied(S_GLOB, *cell)
    t0 = time.time()
    bb, dm, dvv = fit_conv(True, nu_fn, th0=th, dml0=dm, dv0=dvv)
    th = bb.x
    VALS[name] = bb.fun
    extra = ""
    if name == 'BE':
        d = bb.fun - REF_BE
        extra = f"  [5P {REF_BE:.2f}, d={d:+.2f} " \
                f"{'OK' if abs(d) < 2.0 else 'FAIL'}]"
        assert abs(d) < 2.0
    if name == 'U22':
        dd = (bb.fun - VALS['BE']) - REF_U22_DELTA
        extra = f"  [6I delta {REF_U22_DELTA:+.2f}, dd={dd:+.2f} " \
                f"{'OK' if abs(dd) < 2.0 else 'FAIL'}]"
        assert abs(dd) < 2.0
    L.append(f"  {name}: {bb.fun:10.2f}  la0={bb.x[0]:+.3f} f={bb.x[1]:.3f} "
             f"s_int={bb.x[2]:.3f}{extra}  ({(time.time()-t0)/60:.1f} min)")
    print(L[-1], flush=True)
    save()

be = VALS['BE']
d = {n: VALS[n] - be for n, c in CELLS if c is not None}
L.append("")
L.append("Delta vs BE (vertical, point grade): " +
         ", ".join(f"{n} {d[n]:+.2f}" for n in ('U22', 'U12', 'U32',
                                                'U21', 'U23')))
L.append(f"6I tied comparators: (1,1) {CMP_TIED[1]:+.2f}, "
         f"(2,2) {REF_U22_DELTA:+.2f}, (3,3) {CMP_TIED[3]:+.2f}")
best_untied_gain = max(d['U22'] - d[n] for n in ('U12', 'U32', 'U21', 'U23'))
# gain of untied over tied = d[U22] - d[untied] ... negative d = better;
# untied beats tied when d[untied] < d[U22]; the beat margin:
beats = {n: d['U22'] - d[n] for n in ('U12', 'U32', 'U21', 'U23')}
worst = max(beats.values())
mL1 = 0.5*(abs(d['U12'] - d['U22']) + abs(d['U32'] - d['U22']))
mL2 = 0.5*(abs(d['U21'] - d['U22']) + abs(d['U23'] - d['U22']))
L.append("")
L.append("untied-beats-(2,2) margins (positive = untied better): " +
         ", ".join(f"{n} {beats[n]:+.2f}" for n in beats))
b1 = "STRIKE against the per-vertex pairing" if worst > 5.0 else \
     ("SUPPORT ((2,2) best overall)" if worst <= 0.0 else
      "TOLERATED (within point-grade noise)")
b2 = "PASS" if mL1 < mL2 else "FAIL"
L.append(f"B1 (tied survival): worst margin {worst:+.2f} -> {b1}")
L.append(f"B2 (instrument split): mean|d| L1-varied {mL1:.2f} vs "
         f"L2-varied {mL2:.2f} -> {b2} "
         f"(prediction: tail reads L2, rung reads L1)")
L.append("[point grade; population-grade caveat per correction #13]")
print("\n".join(L[-6:]), flush=True)
save()
print("\nSTAGE 6V done -> data/stage6v_untied.txt")
