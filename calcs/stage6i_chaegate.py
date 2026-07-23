"""
STAGE 6I (O16a + the 6H L-contest): measured ambients + the leg count.

Part A -- the L contest (vertical-hardened, global fiducial gate
s = 0.8681): BE vs L = 1/2/3 of the grammar beta = (1/2)*(q_loc*s_amb)^L.
PRE-REGISTERED (6H): the grammar survives iff L=2 ranks best of
{L1, L2, L3}. L=1 carries the deep Bernoulli break (c2 = -0.0252) and
the strongest tail (p = 0.717); L=3 keeps c2 AND c3 (first break at c4,
= -s^3/32 - 1/720). The Bernoulli-break rung IS the leg count, so this
is a direct measurement of L.

Part B -- the measured-ambient leg (O16a): per-galaxy Chae+21 Table 3
gates s_i = n_i/(1+n_i), n_i = n_BE(sqrt(e_i/a0)); maxclust fiducial,
noclust variant; unmatched galaxies + lensing leg carry the
matched-median gate (disclosed). 6H: measured medians are SHARPER than
the 0.02 fiducial (median s_i = 0.932/0.976 vs 0.868). Pre-registered:
per-galaxy >= global on both treatments; plain reaching <= -100
RESOLVES the 6F disclosed partial. STRIKE if per-galaxy < global.
NOTE a value coincidence to avoid confusion: median per-galaxy g
(maxclust) = 0.868 happens to equal the fiducial s = 0.8681 -- the
solver parameter is s, and the per-galaxy s median is 0.932.

Part C -- the AMB quadrupole record (O16c): 5S machinery verbatim,
per-eN solar gates g(e), lock join with the 6G alpha-hat 1.060+-0.024.

Gates: BE regressions vs 5P refs (|d| < 2); AMBg deltas vs 6F (-59.05
vertical / -92.40 plain, |dd| < 2); quadrupole G2 Newton control, G1
simple vs 4K, BE cross-check vs the 5S cache, G4 hi-res.
Writes data/stage6i_chaegate.txt (+ cache data/stage6i_q.npy).
"""
import csv, glob, math, os, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from numpy.polynomial.legendre import leggauss, legval

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25
LN10 = math.log(10)
S_ML = 0.1*LN10

def n_amb_of(e):
    x = math.sqrt(e)
    return 1.0/(math.exp(x) - 1.0)
def g_of(n):
    return (n/(1.0 + n))**2
def s_of(e):
    n = n_amb_of(e)
    return n/(1.0 + n)

S_GLOB = s_of(0.02)          # 0.8681 -- the 6E/6F fiducial gate
G_FID = S_GLOB*S_GLOB        # 0.7536

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

# ---- Chae join
chae = {}
with open('data/chae2021_table3.csv') as f:
    for row in csv.DictReader(l for l in f if not l.startswith('#')):
        chae[row['galaxy']] = (10.0**float(row['log_eN_maxclust']),
                               10.0**float(row['log_eN_noclust']))
sm_raw = np.full(NGal, np.nan)
sn_raw = np.full(NGal, np.nan)
for i in range(NGal):
    nm = gal_name[ug[i]]
    if nm in chae:
        sm_raw[i] = s_of(chae[nm][0])
        sn_raw[i] = s_of(chae[nm][1])
matched = np.isfinite(sm_raw)
NM = int(matched.sum())
S_MED_MX = float(np.median(sm_raw[matched]))
S_MED_NO = float(np.median(sn_raw[matched]))
S_GAL_MX = np.where(matched, sm_raw, S_MED_MX)
S_GAL_NO = np.where(matched, sn_raw, S_MED_NO)
SPT_MAX = S_GAL_MX[gidx]
SPT_NO = S_GAL_NO[gidx]

def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))

def nu_pg(y, spt, Lc):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    spt = np.broadcast_to(np.asarray(spt, float), y.shape)
    sL = spt**Lc
    ly = np.log(y)
    nu = nu_simple(y)
    for _ in range(80):
        d1 = 2.0*nu - 1.0
        b = 0.5*sL/d1**Lc
        db = -Lc*sL/d1**(Lc + 1)
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

def m2hv(th, dml, dv, use_v, spt, sl, Lc):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*nu_pg(gN/a0, spt, Lc)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm_) - (dv[gidx] if use_v else 0.0)
    out = np.sum(r*r/se2 + np.log(se2))
    lg = l_gbar[lmask] + dlt
    rl = l_gobs[lmask] - (lg + np.log10(nu_pg(10**lg/a0, sl, Lc)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(dml*dml)/(S_ML*S_ML)
    if use_v: out += np.sum(dv*dv/(SIGV*SIGV))
    return out

def fit_conv(use_v, spt, sl, Lc, th0=None, dml0=None, dv0=None, tol=0.05,
             max_rounds=15):
    dml = np.zeros(NGal) if dml0 is None else dml0.copy()
    dv = np.zeros(NGal) if dv0 is None else dv0.copy()
    spt_a = np.broadcast_to(np.asarray(spt, float), g_gas.shape)
    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08, 0.0]]
                  if rd == 0 and th0 is None else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2hv(t, dml, dv, use_v, spt_a, sl, Lc),
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
                r0_ = lgobs - np.log10(gN*nu_pg(gN/10**la0, spt_a, Lc))
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
                          - np.log10(gN2*nu_pg(gN2/10**la0, spt_a[mm], Lc))
                          - (dv[gi2] if use_v else 0.0))
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2hv(best.x, dml, dv, use_v, spt_a, sl, Lc)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2hv(t, dml, dv, use_v, spt_a, sl, Lc),
                 list(best.x), method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml, dv

REF_BE = {False: -10435.00, True: -12152.49}
REF_AMBG = {False: -92.40, True: -59.05}

L = [f"STAGE 6I: measured ambients + the leg count -- {kept} galaxies, "
     f"{len(gobs)} points + {int(lmask.sum())} lensing; Chae matched "
     f"{NM}/{NGal} (unmatched + lensing carry the matched median, "
     f"disclosed)",
     f"gates: fiducial s = {S_GLOB:.4f} (g = {G_FID:.4f}); measured "
     f"medians s(maxclust) = {S_MED_MX:.4f}, s(noclust) = {S_MED_NO:.4f}",
     ""]

def save():
    with open('data/stage6i_chaegate.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

save()
VALS = {}
for use_v, tag in ((True, "vertical-hardened"), (False, "plain hier")):
    L.append(f"{tag}:")
    runs = [('BE', 0.0, 0.0, 2), ('AMBg', S_GLOB, S_GLOB, 2),
            ('PGmax', SPT_MAX, S_MED_MX, 2)]
    if use_v:
        runs += [('L3g', S_GLOB, S_GLOB, 3), ('L1g', S_GLOB, S_GLOB, 1)]
    else:
        runs += [('PGno', SPT_NO, S_MED_NO, 2)]
    th, dm, dvv = None, None, None
    for name, spt, sl, Lc in runs:
        t0 = time.time()
        # per-galaxy runs pass the per-POINT gate array; scalars broadcast
        bb, dm, dvv = fit_conv(use_v, spt, sl, Lc, th0=th, dml0=dm,
                               dv0=dvv)
        th = bb.x
        VALS[(use_v, name)] = bb.fun
        extra = ""
        if name == 'BE':
            d = bb.fun - REF_BE[use_v]
            extra = f"  [5P {REF_BE[use_v]:.2f}, d={d:+.2f} " \
                    f"{'OK' if abs(d) < 2.0 else 'FAIL'}]"
            assert abs(d) < 2.0
        if name == 'AMBg':
            dd = (bb.fun - VALS[(use_v, 'BE')]) - REF_AMBG[use_v]
            extra = f"  [6F delta {REF_AMBG[use_v]:+.2f}, dd={dd:+.2f} " \
                    f"{'OK' if abs(dd) < 2.0 else 'FAIL'}]"
            assert abs(dd) < 2.0
        L.append(f"  {name}: {bb.fun:10.2f}  la0={bb.x[0]:+.3f} "
                 f"f={bb.x[1]:.3f} s_int={bb.x[2]:.3f}{extra}  "
                 f"({(time.time()-t0)/60:.1f} min)")
        print(L[-1], flush=True)
        save()
    be = VALS[(use_v, 'BE')]
    row = "  Delta vs BE: " + ", ".join(
        f"{n} {VALS[(use_v, n)]-be:+.2f}"
        for n, _, _, _ in runs if n != 'BE')
    L.append(row)
    L.append("")
    print(row, flush=True)
    save()

# ---------------- Part C: the AMB quadrupole record ----------------
A0_SI = 1.2e-10
R_M_SI = 1.052e15
UNIT_Q = A0_SI/R_M_SI
CASSINI_CAP = 9e-27

def nu_newton(y):
    return np.ones_like(np.asarray(y, float))
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def make_amb_scalar(g):
    s = math.sqrt(g)
    return lambda y: nu_pg(y, s, 2)

def solve_multipole(nu, eN, NR=512, NMU=96, LMAX=16):
    r = np.logspace(-2, 3, NR)
    mu, wmu = leggauss(NMU)
    R, MU = np.meshgrid(r, mu, indexing='ij')
    ST = np.sqrt(1-MU**2)
    Pl = np.zeros((LMAX+1, NMU))
    for l in range(LMAX+1):
        c = np.zeros(l+1); c[l] = 1
        Pl[l] = legval(mu, c)
    gr = -1.0/R**2 + eN*MU
    gt = -eN*ST
    f = nu(np.hypot(gr, gt)) - 1.0
    Ar, At = f*gr, f*gt
    dAr = np.gradient(R**2*Ar, np.log(r), axis=0)/R**3
    dAt = -np.gradient(ST*At, mu, axis=1)/R
    S = -(dAr + dAt)
    Sl = np.array([(2*l+1)/2*np.sum(wmu*S*Pl[l], axis=1)
                   for l in range(LMAX+1)])
    A = np.zeros((LMAX+1, NR)); B = np.zeros((LMAX+1, NR))
    dr = np.diff(r)
    for l in range(LMAX+1):
        g_in = Sl[l]*r
        for i in range(NR-1):
            qf = (r[i]/r[i+1])**(l+1)
            A[l, i+1] = A[l, i]*qf + 0.5*dr[i]*(g_in[i]*qf + g_in[i+1])
        g_out = Sl[l]*r
        for i in range(NR-2, -1, -1):
            qf = (r[i]/r[i+1])**l
            B[l, i] = B[l, i+1]*qf + 0.5*dr[i]*(g_out[i] + g_out[i+1]*qf)
    lv = np.arange(LMAX+1)[:, None]
    phi_l = -(A+B)/(2*lv+1)
    return r, phi_l

def qval(nu, eN, NR=512, NMU=96, LMAX=16):
    r, phi_l = solve_multipole(nu, eN, NR, NMU, LMAX)
    band = (r > 0.02) & (r < 0.2)
    return float(np.median((phi_l[2]/r**2)[band]))

try:
    CACHE = 'data/stage6i_q.npy'
    L.append("Part C: AMB quadrupole record (per-eN solar gates):")
    if os.path.exists(CACHE):
        QA = np.load(CACHE, allow_pickle=True).item()
        L.append("  (solver results loaded from cache)")
    else:
        QA = {}
        r_, phi_l = solve_multipole(nu_newton, 1.2)
        g2 = np.max(np.abs(phi_l[2]))
        L.append(f"  G2 Newton control: max|phi_2| = {g2:.2e} -> "
                 f"{'OK' if g2 < 1e-12 else 'FAIL'}")
        assert g2 < 1e-12
        q_s = qval(nu_simple, 1.2)
        g1 = abs(q_s/(-0.09788) - 1)
        L.append(f"  G1 regression simple eN=1.2: q = {q_s:+.5f} "
                 f"(4K -0.09788) {100*g1:.1f}% -> "
                 f"{'OK' if g1 < 0.02 else 'FAIL'}")
        assert g1 < 0.02
        q_be = qval(nu_be, 1.2)
        if os.path.exists('data/stage5s_q.npy'):
            Q5 = np.load('data/stage5s_q.npy', allow_pickle=True).item()
            gbe = abs(q_be/Q5[0.0][1.2] - 1)
            L.append(f"  G-member BE vs 5S cache: q = {q_be:+.5f} "
                     f"(5S {Q5[0.0][1.2]:+.5f}) {100*gbe:.2f}% -> "
                     f"{'OK' if gbe < 0.005 else 'FAIL'}")
            assert gbe < 0.005
        QA['BE'] = {1.2: q_be}
        QA['amb'] = {}
        for eN in (1.0, 1.2, 1.4):
            gA = g_of(n_amb_of(eN))
            QA['amb'][eN] = qval(make_amb_scalar(gA), eN)
        q_hi = qval(make_amb_scalar(g_of(n_amb_of(1.2))), 1.2,
                    NR=768, NMU=128, LMAX=24)
        g4 = abs(q_hi/QA['amb'][1.2] - 1)
        L.append(f"  G4 hi-res: {q_hi:+.5f} vs {QA['amb'][1.2]:+.5f} "
                 f"({100*g4:.1f}%) -> {'OK' if g4 < 0.05 else 'FAIL'}")
        assert g4 < 0.05
        np.save(CACHE, QA)
    ga = {eN: g_of(n_amb_of(eN)) for eN in (1.0, 1.2, 1.4)}
    L.append(f"  amb gates g(1.0/1.2/1.4) = {ga[1.0]:.4f} {ga[1.2]:.4f} "
             f"{ga[1.4]:.4f}; q = {QA['amb'][1.0]:+.5f} "
             f"{QA['amb'][1.2]:+.5f} {QA['amb'][1.4]:+.5f} "
             f"(BE q(1.2) = {QA['BE'][1.2]:+.5f}, amb/BE - 1 = "
             f"{QA['amb'][1.2]/QA['BE'][1.2]-1:+.2%})")
    AH, dAH = 1.060, 0.024        # 6G six-seed alpha-hat
    Q2 = 3*abs(QA['amb'][1.2])*UNIT_Q*AH
    dq = abs(QA['amb'][1.4]-QA['amb'][1.0])/0.4*0.05
    L.append(f"  lock join: Q2 = 3|q|*unit*alpha-hat = {Q2:.2e} s^-2 = "
             f"{Q2/CASSINI_CAP:.1f}x Cassini (alpha-hat {AH}+-{dAH}; "
             f"g_ext +-0.05a0 -> +-{3*dq*UNIT_Q:.1e})")
except Exception as ex:
    L.append(f"  Part C INCOMPLETE: {type(ex).__name__}: {ex}")
L.append("")

# ---------------- verdict assembly (guarded) ----------------
try:
    vb, va, vp = (VALS[(True, 'BE')], VALS[(True, 'AMBg')],
                  VALS[(True, 'PGmax')])
    v1, v3 = VALS[(True, 'L1g')], VALS[(True, 'L3g')]
    pb, pa, pp, pn = (VALS[(False, 'BE')], VALS[(False, 'AMBg')],
                      VALS[(False, 'PGmax')], VALS[(False, 'PGno')])
    L.append("VERDICT BLOCK:")
    order = sorted([('L1', v1), ('L2', va), ('L3', v3)], key=lambda t: t[1])
    winner = order[0][0]
    L.append(f"  A (leg count, vertical): L1 {v1-vb:+.2f}, L2 {va-vb:+.2f}"
             f", L3 {v3-vb:+.2f} vs BE -> ranking "
             f"{' < '.join(n for n, _ in order)} (lower = better); "
             f"winner {winner} -> "
             f"{'GRAMMAR PASSES (L=2)' if winner == 'L2' else 'GRAMMAR STRIKE'}")
    L.append(f"  B (measured ambients): vertical PGmax {vp-vb:+.2f} vs "
             f"AMBg {va-vb:+.2f} ({'improves' if vp < va else 'HURTS'}); "
             f"plain PGmax {pp-pb:+.2f}, PGno {pn-pb:+.2f} vs AMBg "
             f"{pa-pb:+.2f}; plain <= -100 bar: "
             f"{'RESOLVED' if min(pp, pn)-pb <= -100 else 'still short'}")
except Exception as ex:
    L.append(f"  verdict INCOMPLETE: {type(ex).__name__}: {ex}")

save()
print("\n".join(L[-8:]))
print("\nSTAGE 6I done -> data/stage6i_chaegate.txt")
