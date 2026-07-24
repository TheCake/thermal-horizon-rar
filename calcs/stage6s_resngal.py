"""
STAGE 6S (6R test, galaxy leg): the resolution bath on the hier ladder.

nu_R = the 6R zero-parameter function (beta(y) = (1/2)R^2/(1+R^2),
R = nu*y - sqrt(y), frequency mixing x_eff = y^((1+beta)/2) nu^beta).
Machinery = the 6I/6F hier fit verbatim (m2hv objective, NM + closed-form
dv + per-galaxy dml sweeps), laws BE (regression gate vs 5P refs) and
RESN, both treatments (vertical-hardened + plain hier).

PRE-REGISTERED (6R, committed before execution): PASS if vertical
Delta vs BE <= -40; STRONG if <= -55; STRIKE against the resolution
reading if > -20; else partial. Expected band: gm-grade or better
(deep = BE exact, tail = gm exact).

Writes data/stage6s_resngal.txt.
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

def nu_resn(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    sy = np.sqrt(y); ly = np.log(y)
    nu = nu_simple(y)
    for _ in range(120):
        Rv = nu*y - sy
        R2 = Rv*Rv
        b = 0.5*R2/(1.0 + R2)
        db = 0.5*(2.0*Rv/((1.0 + R2)**2))*y
        lnu = np.log(nu)
        u = np.exp(np.minimum(0.5*(1.0 + b)*ly + b*lnu, 60.0))
        eu = np.exp(np.minimum(u, 60.0))
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        F = nu - 1.0 - n
        dudnu = u*(db*(0.5*ly + lnu) + b/nu)
        dF = 1.0 + (eu/(em1*em1))*dudnu
        step = F/np.where(np.abs(dF) > 1e-12, dF, 1e-12)
        nu = np.maximum(nu - np.clip(step, -0.5*nu, 0.5*nu), 1.0 + 1e-15)
    return nu

NU = {'BE': nu_be, 'RESN': nu_resn}

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

REF_BE = {False: -10435.00, True: -12152.49}

L = [f"STAGE 6S: the resolution bath on the hier galaxy ladder -- "
     f"{kept} galaxies, {len(gobs)} points + {int(lmask.sum())} lensing",
     "pre-registered (6R): PASS <= -40 vertical; STRONG <= -55; "
     "STRIKE > -20; comparators vertical: gm -42.7, AMB -59.05/-61.68, "
     "F4 -64.2",
     ""]

def save():
    with open('data/stage6s_resngal.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

save()
VALS = {}
for use_v, tag in ((True, "vertical-hardened"), (False, "plain hier")):
    L.append(f"{tag}:")
    th, dm, dvv = None, None, None
    for name in ('BE', 'RESN'):
        t0 = time.time()
        bb, dm, dvv = fit_conv(use_v, NU[name], th0=th, dml0=dm, dv0=dvv)
        th = bb.x
        VALS[(use_v, name)] = bb.fun
        extra = ""
        if name == 'BE':
            d = bb.fun - REF_BE[use_v]
            extra = f"  [5P {REF_BE[use_v]:.2f}, d={d:+.2f} " \
                    f"{'OK' if abs(d) < 2.0 else 'FAIL'}]"
            assert abs(d) < 2.0
        L.append(f"  {name}: {bb.fun:10.2f}  la0={bb.x[0]:+.3f} "
                 f"f={bb.x[1]:.3f} s_int={bb.x[2]:.3f}{extra}  "
                 f"({(time.time()-t0)/60:.1f} min)")
        print(L[-1], flush=True)
        save()
    d = VALS[(use_v, 'RESN')] - VALS[(use_v, 'BE')]
    L.append(f"  Delta RESN vs BE: {d:+.2f}")
    L.append("")
    print(L[-2], flush=True)
    save()

dv_ = VALS[(True, 'RESN')] - VALS[(True, 'BE')]
dp_ = VALS[(False, 'RESN')] - VALS[(False, 'BE')]
if dv_ <= -55: v = "STRONG (F4/AMB grade)"
elif dv_ <= -40: v = "PASS"
elif dv_ > -20: v = "STRIKE against the resolution reading"
else: v = "PARTIAL"
L.append(f"VERDICT vs pre-registered bars: vertical {dv_:+.2f}, plain "
         f"{dp_:+.2f} -> {v}")
print(L[-1], flush=True)
save()
print("\nSTAGE 6S done -> data/stage6s_resngal.txt")
