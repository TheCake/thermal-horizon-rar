"""
STAGE 5M (O11): does the hierarchical tail vote survive the vertical channel?

The night's central ambiguity: sharp-tail functions (gm, boot, p065) win the
hierarchical-disk-M/L galaxy contest by 56-85 over the occupation law -- but
the binaries counter-lean 6-10 the other way (5K, sign-consistent 12/12),
their a0 translation under sharp functions runs +4.9 sigma off the horizon
temperature (5L) while the FLAT galaxy fits sit ON it and flat M/L shows no
tail preference at all (5J: BE flat-best by 2). The suspicious coupling: the
hier tail vote lives exactly where a0 runs low (0.84-0.89e-10) and f_ML runs
high (1.5+). The 4W lesson says per-galaxy VERTICAL structure (distance/
inclination, measured priors) is the great absorber at SPARC depth.

Test: add the 4W vertical channel -- per-galaxy offset dv_g with MEASURED
prior sigma_v = hypot((e_D/D)/ln10, 2 cot(i) e_i,rad/ln10) from the SPARC
table, closed-form profiled -- to the 5D converged hierarchical machinery
(delta_d + joint lensing), and re-run the function ladder {BE, p065, gm,
boot}. If Delta(gm - BE) collapses, the tail vote was vertical structure in
disguise (the sober reading: occupation law + horizon a0 coherent across
flat galaxies AND binaries stands); if it survives, the sharp-tail
preference is robust to every per-galaxy channel with measured priors.

Gates: G1 sigma_v -> 1e-4 regression on BE reproduces the 5D value
(-10435.00) within 1.0; G2 nesting -- every dv-ON fit must be <= its dv-OFF
comparator (5D/5G/5J values); G3 trajectory convergence (tol 0.05).
Writes data/stage5m_hierv.txt.
"""
import glob, math, os
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

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def nu_p065(y):
    yc = np.clip(np.asarray(y, float), 1e-14, None)
    ex = np.exp(-np.minimum(yc**0.65, 60.0))
    return (1.0-ex)**(-1.0/1.3)
def nu_gm(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    a = y**0.75
    w = np.sqrt(nu_simple(y))
    for _ in range(30):
        u = np.minimum(a*w, 60.0)
        eu = np.exp(u)
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        H = w*w - 1.0 - n
        dH = 2.0*w + a*eu/(em1*em1)
        w = np.maximum(w - H/dH, 1e-8)
    return w*w
def nu_boot(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    u = 0.5*(y + np.sqrt(y*y + 4.0*y))
    for _ in range(14):
        uc = np.minimum(u, 45.0)
        eu = np.exp(uc)
        em1 = eu - 1.0
        F = u - y - y/em1
        dF = 1.0 + y*eu/(em1*em1)
        u = np.maximum(u - F/dF, 1e-13)
    nu = u/y
    big = y > 45.0
    if np.any(big):
        yb = np.minimum(y, 700.0)
        nu = np.where(big, 1.0 + 1.0/np.expm1(yb), nu)
    return nu

FAMS = {'BE': nu_be, 'p065': nu_p065, 'gm': nu_gm, 'boot': nu_boot}
COMP_OFF = {'BE': -10435.00, 'p065': -10491.38, 'gm': -10519.78,
            'boot': -10510.60}

def m2hv(th, nu, dml, dv, sv, w_g):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*nu(gN/a0)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm_) - dv[gidx]
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    rl = l_gobs[lmask] - (lg + np.log10(nu(10**lg/a0)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    out += np.sum(w_g*dv*dv/(sv*sv))
    return out

def fit_v(nu, sv, tol=0.05, max_rounds=15, th0=None, trace=None):
    w_g = np.ones(NGal)
    dml = np.zeros(NGal)
    dv = np.zeros(NGal)
    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08, 0.0]] if rd == 0 else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2hv(t, nu, dml, dv, sv, w_g), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int, dlt = best.x
        se2c = s_int*s_int
        for _ in range(3):
            # (a) closed-form vertical offsets at current dml
            fac = f*np.exp(dml[gidx])
            gN = g_gas + fac*g_dsk + g_bul
            r0 = lgobs - np.log10(gN*nu(gN/10**la0))
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                w = 1.0/(sig2[mm] + se2c)
                dv[gi2] = np.sum(w*r0[mm])/(np.sum(w) + 1.0/sv[gi2]**2)
            # (b) disk-M/L offsets at current dv
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = (lgobs[mm] - np.log10(gN2*nu(gN2/10**la0))
                          - dv[gi2])
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2hv(best.x, nu, dml, dv, sv, w_g)
        if trace is not None: trace.append(cur)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2hv(t, nu, dml, dv, sv, w_g), list(best.x),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml, dv

L = [f"STAGE 5M vertical-channel disambiguation: {kept} galaxies, "
     f"{len(gobs)} points + {int(lmask.sum())} lensing; sigma_v measured "
     f"16/50/84 = {np.percentile(SIGV,[16,50,84]).round(3).tolist()} dex", ""]

# G1: sigma_v -> 0 regression on BE
bg, _, _ = fit_v(nu_be, np.full(NGal, 1e-4), tol=0.1, max_rounds=12)
ok1 = abs(bg.fun - COMP_OFF['BE']) < 1.0
L.append(f"G1 sigma_v->0 (BE): {bg.fun:.2f} (5D -10435.00) -> "
         f"{'PASS' if ok1 else 'FAIL'}")
L.append("")

res = {}
L.append("Function ladder WITH the measured vertical channel (dv-ON):")
for name in ('BE', 'p065', 'gm', 'boot'):
    tr = []
    b, dml, dv = fit_v(FAMS[name], SIGV, trace=tr)
    res[name] = b
    la0, f, s_int, dlt = b.x
    nest = b.fun <= COMP_OFF[name] + 0.5
    L.append(f"  {name:>5}: {b.fun:10.2f}  ({len(tr)} rd)  a0={10**la0:.3e}"
             f"  f_ML={f:.2f}  s_int={s_int:.3f}  std(dv)={np.std(dv):.3f}"
             f"  [G2 nesting vs {COMP_OFF[name]:.2f}: "
             f"{'OK' if nest else 'FAIL'}]")
L.append("")
d_gm = res['gm'].fun - res['BE'].fun
d_p = res['p065'].fun - res['BE'].fun
d_bo = res['boot'].fun - res['BE'].fun
L.append(f"Delta vs BE, dv-ON:  gm {d_gm:+.2f} | p065 {d_p:+.2f} | "
         f"boot {d_bo:+.2f}")
L.append(f"Delta vs BE, dv-OFF: gm -84.78 | p065 -56.38 | boot -75.59")
L.append("")
if d_gm > -15:
    L.append("VERDICT: the tail vote COLLAPSES under the measured vertical "
             "channel -> it was per-galaxy distance/inclination structure; "
             "the occupation law + horizon a0 reading stands coherent "
             "(flat galaxies + binaries).")
elif d_gm < -50:
    L.append("VERDICT: the tail vote SURVIVES the vertical channel -> "
             "robust to every measured per-galaxy nuisance; the two-system "
             "tension is physical.")
else:
    L.append("VERDICT: partial absorption -- the tail vote is "
             "vertical-degenerate at SPARC depth; decision needs anchors "
             "(O1b).")

out = "\n".join(L)
print(out)
with open('data/stage5m_hierv.txt', 'w') as f:
    f.write(out+"\n")
print("\nsaved: data/stage5m_hierv.txt")
