"""
STAGE 7D: the VACUUM-SHARE contest -- is the "+1" in the gate's local
factor demanded by the sky?

The AMB local factor q = 1/(2nu-1) = 1/(2n+1) is the vacuum share of the
JC dispersive pull lambda^2(2n+1)/Delta (6H/6U, derivation-grade). The
(2n+1) is the QUANTUM symmetrized noise (n + 1/2 per mode); a CLASSICAL
bath's pull carries 2n only -- no zero-point, no spontaneous channel.
6U's rival table contested the AMBIENT factor; the local factor was never
contested. Two classical replacements, fit at full grade on both systems:

  QCL (the vacuum-share contest, PRIMARY): q_cl = 1/(2(nu-1)) = 1/(2n) --
      the same admixture grammar built on the classical pull. beta_cl =
      g/(8(nu-1)^2) diverges toward the Newtonian end (no vacuum floor);
      capped at BETA_CAP = 8 (cap binds only beyond x ~ 2.3 at galaxy g;
      cap-16 galaxy robustness run; cap range disclosed).
      Deep limit: beta_cl -> 0, BE ladder preserved -- the contest is
      pure transition+tail shape, exactly where the data vote (5T).
  RJA (ambient statistics, SECONDARY): n_amb classical-equipartition
      n_RJ = 1/x_amb instead of n_BE (s_amb = 1/(1+x)). Galaxies sit
      deep-ambient (x ~ 0.14: g 0.754 -> 0.768, near-null by
      construction -- BE->RJ converge in occupied ambients); the rung
      reads through the BINARIES (x ~ 1.1: g 0.112 -> 0.228).

Also QUOTED (banked rungs, no new fits): R0 "no vacuum -> zero share ->
beta = 0 -> pure BE" is already excluded by the galaxy gate credit
(-59 vertical / -100 plain); R1 the ladder's c1 = 0 branch (no zero-point
rung) is excluded in BOTH systems (4S galaxies 7.5-sigma profile / 95.5%
bootstrap; 4X binaries ~20 lnL/seed) -- with the open nuance that the
full classical self-consistent FUNCTION (simple-nu) is hier-rejected on
galaxies (-99) while binaries retain a residual vtilde-shape lean toward
it (-12.3 +- 2.2, quoted openly).

Solver gates: GQ1 g->0 recovers BE (<1e-9); GQ2 Newton residual < 1e-8
on the y-grid; GQ3 cap-binding range disclosed; spherical identity per
binary table (<2%); galaxy-side BE regression vs 5P refs (assert d<2);
AMB node cross-check vs 6F -59.05 (disclose).

PRE-REGISTERED BARS (both directions):
  QCL galaxies (vertical): REJECTED if Delta(QCL-AMB) >= +10; KILL-grade
    >= +25; UNRESOLVED |Delta| < 10; QCL BETTER <= -10 = FLAG (the
    vacuum-share reading of 6H takes a real seam strike).
  QCL binaries (2 seeds, 31/101): REJECTED if mean penalty vs amb
    >= +5/seed sign-consistent or alpha-hat edge-rides both seeds;
    KILL-grade >= +10/seed.
  VERDICT "VACUUM RUNG MEASURED" only if QCL rejected in BOTH systems
    (conditional-on-grammar caveat always stated).
  RJA binaries (4 seeds): REJECTED if >= +5/seed sign-consistent;
    UNRESOLVED otherwise (expected: the weak rung); galaxies reported.
Writes data/stage7d_vacuumshare.txt (incremental) +
data/stage7d_summary.txt (binary runs).
"""
import glob, math, os, re, sys, time, gc
import numpy as np
from scipy.optimize import minimize, minimize_scalar

BETA_CAP = 8.0

# ---------------------------------------------------------------- solvers
def nu_simple_arr(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))

def nu_be(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    x = np.sqrt(y)
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))

def make_amb(g):
    def nu_run(y):
        y = np.clip(np.asarray(y, float), 1e-14, None)
        ly = np.log(y)
        nu = nu_simple_arr(y)
        for _ in range(120):
            b = g*0.5/((2.0*nu - 1.0)**2)
            db = g*(-2.0)/((2.0*nu - 1.0)**3)
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

def make_qcl(g, bcap=BETA_CAP):
    def nu_run(y):
        y = np.clip(np.asarray(y, float), 1e-14, None)
        ly = np.log(y)
        nu = nu_simple_arr(y)
        for _ in range(160):
            d = np.maximum(2.0*(nu - 1.0), 1e-12)
            braw = g*0.5/(d*d)
            b = np.minimum(braw, bcap)
            db = np.where(braw < bcap, -2.0*g/(d*d*d), 0.0)
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

def resid(nu_fn, g, classical, y, bcap=BETA_CAP):
    nu = nu_fn(y)
    if classical:
        d = np.maximum(2.0*(nu - 1.0), 1e-12)
        b = np.minimum(g*0.5/(d*d), bcap)
    else:
        b = g*0.5/((2.0*nu - 1.0)**2)
    u = np.exp(np.minimum(0.5*(1.0+b)*np.log(y) + b*np.log(nu), 60.0))
    n = np.where(u < 60.0, 1.0/np.maximum(np.exp(np.minimum(u, 60.0))-1.0,
                                          1e-300), 0.0)
    return np.max(np.abs(nu - 1.0 - n))

def n_amb_of(e):
    x = math.sqrt(e)
    return 1.0/(math.exp(x) - 1.0)
def n_rj_of(e):
    return 1.0/math.sqrt(e)
def g_of(n):
    return (n/(1.0 + n))**2

G_GAL = g_of(n_amb_of(0.02))          # 0.754 (6F fiducial)
G_GAL_RJ = g_of(n_rj_of(0.02))        # 0.768

L = ["STAGE 7D: the vacuum-share contest (bars in script header)", ""]
def save():
    with open('data/stage7d_vacuumshare.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

# AMENDMENT (post-commit c53f2e2, PRE-RESULTS -- the gate fired before
# any fit number existed; 6O/6X precedent), two parts:
# (1) the g->0 limit is POINTWISE, not uniform: with a fixed absolute
#     cap the far tail (nu-1 <= sqrt(g/64)) is capped for ANY g, so the
#     original y<=1e3 grid read the cap, not the solver;
# (2) even cap-free, convergence carries the O(g/(nu-1)) first-order
#     amplification (measured 6.35e-8 at g=1e-12, y<=100 -- exactly the
#     analytic response, not an error). GQ1 now: BE-recovery on the
#     physics window y <= 30 at g = 1e-12 (threshold 1e-9) PLUS the
#     real regression that the deviation scales linearly in g
#     (dev(1e-12)/dev(1e-13) in [7, 13]).
YG = np.geomspace(1e-3, 30.0, 400)
q0 = make_qcl(1e-12)(YG)
gq1 = float(np.max(np.abs(q0/nu_be(YG) - 1.0)))
q0b = make_qcl(1e-13)(YG)
gq1b = float(np.max(np.abs(q0b/nu_be(YG) - 1.0)))
gq1_ratio = gq1/max(gq1b, 1e-300)
YGW = np.geomspace(1e-3, 1e3, 400)
r_amb = resid(make_amb(G_GAL), G_GAL, False, YGW)
r_qcl = resid(make_qcl(G_GAL), G_GAL, True, YGW)
nu_q = make_qcl(G_GAL)(YGW)
d_ = np.maximum(2.0*(nu_q-1.0), 1e-12)
capm = (G_GAL*0.5/(d_*d_)) >= BETA_CAP
ycap = YGW[capm]
gq1_ok = (gq1 < 1e-9) and (7.0 < gq1_ratio < 13.0)
L.append(f"GQ1 g->0 recovers BE on y<=30: dev(1e-12) = {gq1:.2e}, "
         f"dev ratio 1e-12/1e-13 = {gq1_ratio:.1f} (linear-in-g) -> "
         f"{'PASS' if gq1_ok else 'FAIL'}")
L.append(f"GQ2 Newton residual: amb {r_amb:.2e}, qcl {r_qcl:.2e} -> "
         f"{'PASS' if max(r_amb, r_qcl) < 1e-8 else 'FAIL'}")
L.append(f"GQ3 cap (galaxy g={G_GAL:.3f}): binds for y >= "
         f"{ycap.min():.2f}" if len(ycap) else "GQ3: cap never binds (gal)")
nu1 = {'AMB': float(make_amb(G_GAL)(np.array([1.0]))[0]),
       'QCL': float(make_qcl(G_GAL)(np.array([1.0]))[0]),
       'RJA': float(make_amb(G_GAL_RJ)(np.array([1.0]))[0])}
L.append(f"nu(1) galaxy-gate: AMB {nu1['AMB']:.4f} | QCL {nu1['QCL']:.4f} "
         f"| RJA {nu1['RJA']:.4f}")
L.append("")
save()
assert gq1_ok and max(r_amb, r_qcl) < 1e-8

# ------------------------------------------------- binary tables (QUMOND)
# AMENDMENT 2 (post-commit, PRE-RESULTS): the capped-classical tail is a
# CLIFF (nu-1 ~ exp(-y^4.5)); at the solver's stock NR=512 the spherical
# identity failed at 2%+ (adjacent-grid nu-1 jumps ~1.9x at the cliff).
# 7D tables run at NR=2048 (same solver, 4x radial resolution); identity
# gate graded: <2% PASS, 2-5% DISCLOSED (QCL binary lnL carries a table
# systematic, stated), >=5% abort the binary leg (galaxy leg -- pointwise
# exact, no PDE -- carries the contest alone).
src_solver = open('calcs/qumond_efe_solver.py', encoding='utf-8-sig').read()
assert src_solver.count("NR, NMU, LMAX = 512, 96, 16") == 1
src_solver = src_solver.replace("NR, NMU, LMAX = 512, 96, 16",
                                "NR, NMU, LMAX = 2048, 96, 16", 1)
ns = {}
exec(src_solver.split('# --- Gate 1')[0], ns)
solve, r = ns['solve'], ns['r']
y = 1.0/r**2
win = (y > 0.05) & (y < 30)

BIN_OK = {'qcl': True, 'rja': True}
t0 = time.time()
for fam, mk, gfun in (('qcl', make_qcl, lambda e: g_of(n_amb_of(e))),
                      ('rja', make_amb, lambda e: g_of(n_rj_of(e)))):
    for eN, tag in ((1.0, 'g1p0'), (1.2, 'g1p2'), (1.4, 'g1p4')):
        pth = f'data/efe_boost_{fam}_{tag}.npy'
        g = gfun(eN)
        if os.path.exists(pth):
            msg = f"table {pth} exists (g={g:.4f})"
            L.append(msg); print(msg, flush=True); save(); continue
        nu = mk(g)
        b = solve(nu, eN)
        np.save(pth, np.stack([y, b]))
        msg = f"table {pth} (eN={eN}, g={g:.4f})"
        if eN == 1.2:
            b_iso = solve(nu, 0.0)
            t = np.max(np.abs(b_iso[win]/nu(y)[win]-1.0))
            grade = ('PASS' if t < 0.02 else
                     'DISCLOSED-SYSTEMATIC' if t < 0.05 else 'ABORT-LEG')
            msg += f"  [spherical identity {100*t:.2f}% {grade}]"
            if t >= 0.05: BIN_OK[fam] = False
        L.append(msg); print(msg, flush=True); save()
L.append(f"({(time.time()-t0)/60:.1f} min tables, NR=2048)"); save()

# ------------------------------------------------- binary patch-runs
src = open('calcs/stage3p_v7budget.py', encoding='utf-8-sig').read()
OLD_DATA = """vt_d = (4.74047/plx[ok]*np.hypot(d['pmra1'][ok]-d['pmra2'][ok],
        d['pmdec1'][ok]-d['pmdec2'][ok]))/(0.9417*np.sqrt(Mtot_d/s_d))
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
vx_ = d['pmra2']-d['pmra1']; vy_ = d['pmdec2']-d['pmdec1']"""
NEW_DATA = """dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
sn_ = np.maximum(np.hypot(sx_, sy_), 1e-12)
sux_, suy_ = sx_/sn_, sy_/sn_
_cn = list(d.columns.names)
def _pick(*names):
    for n in names:
        if n in _cn: return np.asarray(d[n], dtype=np.float64)
    raise KeyError(names)
r1_ = _pick('radial_velocity1', 'dr2_radial_velocity1', 'rv1')
r2_ = _pick('radial_velocity2', 'dr2_radial_velocity2', 'rv2')
try:
    er1_ = _pick('radial_velocity_error1', 'dr2_radial_velocity_error1')
    er2_ = _pick('radial_velocity_error2', 'dr2_radial_velocity_error2')
except KeyError:
    er1_ = np.full(len(r1_), 2.0); er2_ = np.full(len(r2_), 2.0)
h1_, h2_ = np.isfinite(r1_), np.isfinite(r2_)
w1_ = np.where(h1_, 1.0/np.maximum(er1_, 0.5)**2, 0.0)
w2_ = np.where(h2_, 1.0/np.maximum(er2_, 0.5)**2, 0.0)
rvs_ = (np.where(h1_, r1_, 0.0)*w1_ + np.where(h2_, r2_, 0.0)*w2_) \\
       / np.maximum(w1_+w2_, 1e-12)
th_ = sep/(2.06265e8/plx)
pmcor_ = rvs_*th_*plx/4.74047
vx_ = d['pmra2']-d['pmra1'] + pmcor_*sux_
vy_ = d['pmdec2']-d['pmdec1'] + pmcor_*suy_
vt_d = (4.74047/plx[ok]*np.hypot(vx_[ok], vy_[ok]))/(0.9417*np.sqrt(Mtot_d/s_d))"""
assert src.count(OLD_DATA) == 1
src = src.replace(OLD_DATA, NEW_DATA)
OLD_TABS = """if GTAG == '1p9':
    PS, PB = 'data/efe_boost_simple.npy', 'data/efe_boost_be.npy'
else:
    PS = f'data/efe_boost_simple_g{GTAG}.npy'
    PB = f'data/efe_boost_be_g{GTAG}.npy'
LNY_U, TAB_S = load_tab(PS)
_,     TAB_B = load_tab(PB)"""
NEW_TABS = """LNY_U, TAB_S = load_tab(LAMPATH)
TAB_B = TAB_S"""
assert src.count(OLD_TABS) == 1
src = src.replace(OLD_TABS, NEW_TABS)
OLD_LOOP = 'for law, TAB in (("simple", TAB_S), ("BE", TAB_B)):'
NEW_LOOP = 'for law, TAB in ((LAMLAW, TAB_S),):'
assert src.count(OLD_LOOP) == 1
src = src.replace(OLD_LOOP, NEW_LOOP)
OLD_SUM = """    P(f"  seed {seed}: BE-minus-simple best lnL = "
      f"{best_lnl['BE']-best_lnl['simple']:+.1f}  "
      f"({(time.time()-t0)/60:.1f} min)")"""
NEW_SUM = """    P(f"  seed {seed}: best lnL = {best_lnl[LAMLAW]:+.3f}  "
      f"({(time.time()-t0)/60:.1f} min)")"""
assert src.count(OLD_SUM) == 1
src = src.replace(OLD_SUM, NEW_SUM)
OLD_OUT = "with open('data/stage3u_summary.txt', 'a') as f:"
NEW_OUT = "with open('data/stage7d_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1
src = src.replace(OLD_OUT, NEW_OUT)

def have(law, k):
    return os.path.exists('data/stage7d_summary.txt') and \
        open('data/stage7d_summary.txt').read().count(f' {law}: ') >= k

for law, seeds in (('qcl', ['31', '101']),
                   ('rja', ['31', '101', '202', '303'])):
    if not BIN_OK[law]:
        L.append(f"binary leg SKIPPED for {law} (identity >= 5%)"); save()
        continue
    if have(law, len(seeds)): continue
    ns2 = {'__name__': '__main__',
           'LAMPATH': f'data/efe_boost_{law}_g1p2.npy', 'LAMLAW': law}
    sys.argv = ['stage7d', '1p2'] + seeds
    print(f"\n===== {law} on binaries =====", flush=True)
    exec(compile(src, f'stage3p_patched_7d_{law}', 'exec'), ns2)
    del ns2; gc.collect()
L.append("binary runs done"); save()

# ------------------------------------------------- galaxy nodes (6S clone)
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
gmap = {g_: i for i, g_ in enumerate(ug)}
gidx = np.array([gmap[g_] for g_ in gal_id])
SIGV = np.array([sigv_g_map[g_] for g_ in ug])
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]
ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

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
        for t0_ in starts:
            b = minimize(lambda t: m2hv(t, dml, dv, use_v, nu_fn),
                         t0_, method='Nelder-Mead',
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
NU = {'BE': nu_be, 'AMB': make_amb(G_GAL), 'QCL': make_qcl(G_GAL),
      'QCL16': make_qcl(G_GAL, 16.0), 'RJA': make_amb(G_GAL_RJ)}
VALS = {}
for use_v, tag, names in ((True, "vertical-hardened",
                           ('BE', 'AMB', 'QCL', 'QCL16', 'RJA')),
                          (False, "plain hier", ('BE', 'AMB', 'QCL'))):
    L.append(f"{tag}:")
    th, dm, dvv = None, None, None
    for name in names:
        t0_ = time.time()
        bb, dm, dvv = fit_conv(use_v, NU[name], th0=th, dml0=dm, dv0=dvv)
        th = bb.x
        VALS[(use_v, name)] = bb.fun
        extra = ""
        if name == 'BE':
            d = bb.fun - REF_BE[use_v]
            extra = f"  [5P {REF_BE[use_v]:.2f}, d={d:+.2f} " \
                    f"{'OK' if abs(d) < 2.0 else 'FAIL'}]"
            assert abs(d) < 2.0
        if name == 'AMB' and use_v:
            extra = f"  [6F fiducial delta was -59.05; here " \
                    f"{bb.fun - VALS[(use_v, 'BE')]:+.2f}]"
        L.append(f"  {name}: {bb.fun:10.2f}  la0={bb.x[0]:+.3f} "
                 f"f={bb.x[1]:.3f} s_int={bb.x[2]:.3f}{extra}  "
                 f"({(time.time()-t0_)/60:.1f} min)")
        print(L[-1], flush=True)
        save()
    for name in names[1:]:
        L.append(f"  Delta {name} vs AMB: "
                 f"{VALS[(use_v, name)] - VALS[(use_v, 'AMB')]:+.2f}   "
                 f"(vs BE: {VALS[(use_v, name)] - VALS[(use_v, 'BE')]:+.2f})")
    L.append("")
    save()

# ------------------------------------------------- parse + verdict
def parse(path):
    txt = open(path).read()
    out = {}
    for m in re.finditer(
            r'seed (\d+) (\w+): a_hat=([0-9.]+) \(grid ([0-9.]+), '
            r'interior=(\w+)\), dlnL\(Newton\)=\+([0-9.]+), wr=([0-9.]+),'
            r'.*?\n\s*seed \1: best lnL = (-?[0-9.]+)', txt):
        s, law, ah, agrid, inter, dn, wr, lnl = m.groups()
        out[(law, int(s))] = dict(a=float(ah), interior=(inter == 'True'),
                                  lnl=float(lnl), dn=float(dn))
    return out

REC = {}
for p in ('data/stage5k_summary.txt', 'data/stage5o_summary.txt',
          'data/stage6g_summary.txt', 'data/stage7d_summary.txt'):
    if os.path.exists(p): REC.update(parse(p))

L.append("binaries (corrected velocities; vs amb and p050 comparators):")
res_bin = {}
for law, seeds in (('qcl', [31, 101]), ('rja', [31, 101, 202, 303])):
    ds_amb, ds_p, ahs, ints = [], [], [], 0
    for s in seeds:
        if (law, s) not in REC: continue
        rr = REC[(law, s)]
        ds_amb.append(rr['lnl'] - REC[('amb', s)]['lnl'])
        ds_p.append(rr['lnl'] - REC[('p050', s)]['lnl'])
        ahs.append(rr['a']); ints += rr['interior']
    if not ds_amb:
        L.append(f"  {law}: NO RUNS PARSED"); continue
    ds_amb, ds_p, ahs = map(np.array, (ds_amb, ds_p, ahs))
    res_bin[law] = dict(damb=ds_amb, ah=ahs, ints=ints, n=len(ds_amb))
    L.append(f"  {law}: lnL-amb per seed " +
             " ".join(f"{d:+.2f}" for d in ds_amb) +
             f"  (mean {ds_amb.mean():+.2f}); vs p050 mean {ds_p.mean():+.2f}; "
             f"a_hat {ahs.mean():.3f}+-{ahs.std(ddof=1)/max(np.sqrt(len(ahs)),1):.3f}; "
             f"interior {ints}/{len(ds_amb)}")
L.append("")

dq_v = VALS.get((True, 'QCL'), np.nan) - VALS.get((True, 'AMB'), np.nan)
dq16 = VALS.get((True, 'QCL16'), np.nan) - VALS.get((True, 'AMB'), np.nan)
dq_p = VALS.get((False, 'QCL'), np.nan) - VALS.get((False, 'AMB'), np.nan)
dr_v = VALS.get((True, 'RJA'), np.nan) - VALS.get((True, 'AMB'), np.nan)
gal_rej = dq_v >= 10.0
gal_kill = dq_v >= 25.0
bq = res_bin.get('qcl')
bin_rej = bq is not None and bq['n'] >= 2 and \
    ((np.all(bq['damb'] < 0) and abs(bq['damb'].mean()) >= 5.0) or
     bq['ints'] == 0)
bin_kill = bq is not None and bq['n'] >= 2 and \
    np.all(bq['damb'] < 0) and abs(bq['damb'].mean()) >= 10.0
if gal_rej and bin_rej:
    verd = "VACUUM RUNG MEASURED (QCL rejected both systems" + \
           (", KILL-grade" if gal_kill and bin_kill else "") + \
           "; conditional on the AMB grammar)"
elif dq_v <= -10.0:
    verd = "FLAG: QCL BETTER on galaxies -- the vacuum-share reading " \
           "takes a seam strike"
else:
    verd = "PARTIAL/UNRESOLVED -- see components"
L.append(f"QCL: galaxy vertical {dq_v:+.2f} (cap16 {dq16:+.2f}; plain "
         f"{dq_p:+.2f}); binary see above")
L.append(f"RJA: galaxy vertical {dr_v:+.2f} (near-null expected)")
L.append(f"VERDICT: {verd}")
L.append("caveats: conditional on the admixture grammar (like all 6U")
L.append("rivals); beta cap disclosed; binary comparators share seeds.")
save()
print("\nSTAGE 7D done -> data/stage7d_vacuumshare.txt")
