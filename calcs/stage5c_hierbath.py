"""
STAGE 5C (O6): the bath matrix under hierarchical M/L -- the live 1/4-vs-1/2
rung.

Stage 4F ranked the four bath cells (+ the dead-branch control) under FLAT
M/L: BE(1/2) -8397.72 > simple -8379.01 > boot(1/4) -8370.46 > standard
-8341.95 > CSD -8169.32. Stage 4Z then showed the continuous-lambda profile
RELOCATES from lam=0.90 to lam=0.52 (c1 0.45 -> 0.26) once per-galaxy disk
M/L is profiled -- peaking AT the 1/4 value. So the 4F ordering is
flat-M/L-conditional, and the open question is whether the quantum-bootstrap
1/4-function (its own c2=7/96, e^-y screening) overtakes BE once the same
hierarchical freedom is granted to every family.

This script re-ranks all five functions under the 4Z treatment: per-galaxy
disk-M/L offsets (0.1-dex lognormal prior) profiled at every fit, joint
SPARC+lensing scatter-marginalized objective (4E fiducial config).
100-rep galaxy bootstrap (offsets refit per rep) for the sign stability of
boot-BE, boot-simple, BE-simple.

Gates:
  G1 prior->0 regression: with sigma_ML = 1e-4 each family must reproduce
     its Stage 4F flat -2lnL within 1.0 (all five).
  G2 injection: synthetic data at BE truth + drawn offsets: hierarchical
     BE must beat hierarchical boot on the injected set (sign check that
     the machinery can still tell the functions apart with offsets free).
Writes data/stage5c_hierbath.txt.
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
        meta[t[0]] = (float(t[5]), int(t[17]))
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q = meta.get(name, (0, 3))
    if inc < 30 or q > 2: continue
    kept += 1
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
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

# ---------------- the five bath functions (verbatim 4F) ----------------
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_csd(y):
    return 1.0+1.0/np.sqrt(np.clip(y, 1e-14, None))
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

FAMS = {'BE': nu_be, 'simple': nu_simple, 'boot': nu_boot,
        'standard': nu_standard, 'CSD': nu_csd}
REF4F = {'BE': -8397.72, 'simple': -8379.01, 'boot': -8370.46,
         'standard': -8341.95, 'CSD': -8169.32}
PAIR_NOTE = {'BE': 'quantum/source (c1=1/2, c2=1/12)',
             'simple': 'classical/self-consistent (c1=1/2, c2=1/8)',
             'boot': 'quantum/self-consistent (c1=1/4, c2=7/96)',
             'standard': 'dead-branch control (c1=0)',
             'CSD': 'classical/source (c1=1; Cassini-dead a priori)'}

# ---------------- 4Z objective with nu passed ----------------
def m2h(th, nu, dml, w_g, pd):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm = gN*nu(gN/a0)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm)
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    yl = 10**lg/a0
    rl = l_gobs[lmask] - (lg + np.log10(nu(yl)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(pd*pd)
    return out

def fit_node(nu, w_g, pd=S_ML, th0=None, dml0=None, rounds=3,
             lg_override=None):
    global lgobs
    lg_save = lgobs
    if lg_override is not None: lgobs = lg_override
    dml = np.zeros(NGal) if dml0 is None else dml0.copy()
    best = None
    try:
        for rd in range(rounds):
            starts = ([list(best.x)] if best is not None else []) + \
                     ([list(th0)] if th0 is not None else []) + \
                     [[math.log10(A0_FID), 1.0, 0.08, 0.0]]
            bb = None
            for t0 in starts:
                b = minimize(lambda t: m2h(t, nu, dml, w_g, pd), t0,
                             method='Nelder-Mead',
                             options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
                if bb is None or b.fun < bb.fun: bb = b
            best = bb
            la0, f, s_int, dlt = best.x
            if pd < 1e-3: continue
            se2c = s_int*s_int
            for gi2 in range(NGal):
                if w_g[gi2] == 0: dml[gi2] = 0.0; continue
                mm = GIDXS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = lgobs[mm] - np.log10(gN2*nu(gN2/10**la0))
                    s2 = sig2[mm] + se2c
                    return (w_g[gi2]*np.sum(rr*rr/s2)
                            + w_g[gi2]*dl*dl/(pd*pd))
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
    finally:
        lgobs = lg_save
    return best, dml

ones = np.ones(NGal)
L = [f"STAGE 5C bath matrix under hierarchical M/L: {kept} galaxies, "
     f"{len(gobs)} points + {int(lmask.sum())} lensing; delta prior 0.1 dex",
     ""]
for n_, note in PAIR_NOTE.items(): L.append(f"  {n_:>8}: {note}")
L.append("")

# ---------------- G1: prior->0 endpoints vs 4F flat ----------------
g1ok = True
L.append("G1 prior->0 regression vs 4F flat table:")
for name, nu in FAMS.items():
    b, _ = fit_node(nu, ones, pd=1e-4, rounds=1)
    ok = abs(b.fun - REF4F[name]) < 1.0
    g1ok &= ok
    L.append(f"  {name:>8}: {b.fun:10.2f}  (4F {REF4F[name]:10.2f})  "
             f"{'OK' if ok else 'FAIL'}")
L.append(f"G1 -> {'PASS' if g1ok else 'FAIL'}")
L.append("")

# ---------------- hierarchical contest ----------------
L.append("Hierarchical contest (delta_d profiled, 0.1-dex prior):")
hier = {}
dms = {}
for name, nu in FAMS.items():
    b, dm = fit_node(nu, ones)
    hier[name] = b
    dms[name] = dm
    la0, f, s_int, dlt = b.x
    L.append(f"  {name:>8}: -2lnL+prior = {b.fun:10.2f}  a0={10**la0:.3e}  "
             f"f_ML={f:.2f}  s_int={s_int:.3f}  dlt={dlt:+.3f}  "
             f"std(dml)={np.std(dm)/LN10:.3f} dex")
hb, he, hs, hstd = (hier[k].fun for k in ('boot', 'BE', 'simple', 'standard'))
L.append(f"  Delta: boot - BE = {hb-he:+.2f} | boot - simple = {hb-hs:+.2f} | "
         f"BE - simple = {he-hs:+.2f} | boot - standard = {hb-hstd:+.2f}")
L.append(f"  [flat 4F comparators: boot-BE +27.27, boot-simple +8.55, "
         f"BE-simple -18.71, boot-standard -28.51]")
L.append("")

# ---------------- G2: injection at BE truth ----------------
rng = np.random.default_rng(11)
d_t = rng.normal(0, S_ML, NGal)
gN_t = g_gas + 1.0*np.exp(d_t[gidx])*g_dsk + g_bul
lg_inj = (np.log10(gN_t*nu_be(gN_t/A0_FID))
          + rng.normal(0, np.sqrt(sig2 + 0.08**2)))
bi_be, _ = fit_node(nu_be, ones, rounds=2, lg_override=lg_inj)
bi_bo, _ = fit_node(nu_boot, ones, rounds=2, lg_override=lg_inj)
ok2 = bi_be.fun < bi_bo.fun
L.append(f"G2 injection (BE truth + offsets): hier BE {bi_be.fun:.2f} vs "
         f"hier boot {bi_bo.fun:.2f} (BE - boot = {bi_be.fun-bi_bo.fun:+.2f})"
         f" -> {'PASS' if ok2 else 'FAIL'}")
L.append("")

# ---------------- bootstrap: sign stability ----------------
rng2 = np.random.default_rng(47)
allg = np.arange(NGal)
d_bb, d_bs, d_es = [], [], []
th_w = {k: hier[k].x for k in ('boot', 'BE', 'simple')}
for k in range(100):
    pick = rng2.choice(allg, NGal, replace=True)
    w = np.zeros(NGal)
    for g_ in pick: w[g_] += 1
    lo = l_gobs + rng2.normal(0, np.sqrt(l_sig2))
    lgobs_save2 = l_gobs.copy()
    vals = {}
    for name in ('boot', 'BE', 'simple'):
        nu = FAMS[name]
        dml_b = np.zeros(NGal)
        best = None
        for rd in range(2):
            bb = None
            for t0 in ([best.x.tolist()] if best is not None else []) + \
                      [list(th_w[name])]:
                def ob(t):
                    la0, f, s_int, dlt = t
                    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5):
                        return 1e12
                    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8:
                        return 1e12
                    a0 = 10**la0
                    fac = f*np.exp(dml_b[gidx])
                    gN = g_gas + fac*g_dsk + g_bul
                    gm = gN*nu(gN/a0)
                    se2 = sig2 + s_int*s_int
                    r = lgobs - np.log10(gm)
                    out = np.sum(w[gidx]*(r*r/se2 + np.log(se2)))
                    lg = l_gbar[lmask] + dlt
                    rl = lo[lmask] - (lg + np.log10(nu(10**lg/a0)))
                    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
                    out += (dlt/DELTA_PRIOR)**2
                    out += np.sum(w*dml_b*dml_b)/(S_ML*S_ML)
                    return out
                b = minimize(ob, t0, method='Nelder-Mead',
                             options=dict(maxiter=3000, xatol=1e-6, fatol=1e-6))
                if bb is None or b.fun < bb.fun: bb = b
            best = bb
            la0, f, s_int, dlt = best.x
            for gi2 in range(NGal):
                if w[gi2] == 0: dml_b[gi2] = 0.0; continue
                mm = GIDXS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = lgobs[mm] - np.log10(gN2*nu(gN2/10**la0))
                    s2 = sig2[mm] + s_int*s_int
                    return w[gi2]*(np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML))
                dml_b[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                             method='bounded').x
        vals[name] = best.fun
    d_bb.append(vals['boot'] - vals['BE'])
    d_bs.append(vals['boot'] - vals['simple'])
    d_es.append(vals['BE'] - vals['simple'])
d_bb, d_bs, d_es = map(np.array, (d_bb, d_bs, d_es))
L.append(f"bootstrap (100 reps, hierarchical, offsets refit per rep):")
L.append(f"  boot - BE     = {d_bb.mean():+7.2f} +/- {d_bb.std(ddof=1):6.2f}"
         f"  (boot better in {int((d_bb<0).sum())}/100)")
L.append(f"  boot - simple = {d_bs.mean():+7.2f} +/- {d_bs.std(ddof=1):6.2f}"
         f"  (boot better in {int((d_bs<0).sum())}/100)")
L.append(f"  BE - simple   = {d_es.mean():+7.2f} +/- {d_es.std(ddof=1):6.2f}"
         f"  (BE better in {int((d_es<0).sum())}/100)")

out = "\n".join(L)
print(out)
with open('data/stage5c_hierbath.txt', 'w') as f:
    f.write(out+"\n")
print("\nsaved: data/stage5c_hierbath.txt")
