"""STAGE 10Y -- THE RESIDUAL-CONTRAST ORGANIZER.
Pre-registration: PREREG-CONTRAST-10Y.md (committed pre-sky).

What organizes the anchored inner-minus-outer residual depression
(R49 C17: -0.0857 +/- 0.0301 dex, 20 galaxies, 20-arcsec split)?
Candidates: acceleration y (form misfit), angular radius (the 10V
resolution rival), disk-scaled radius R/1.5Rd (the 7B/7C inner-disk
object). No new fits: all statistics are functions of residuals at
the published parent optima (round49_gb.py construction, verbatim).

Modes:
  py calcs/stage10y_organizer.py gates   -> data/stage10y_gates.txt
  py calcs/stage10y_organizer.py sky     -> data/stage10y_skyread.txt

Gate/sky separation (prereg section 2): gates touch published
regression targets, the DESIGN (y/theta/R/Rd censuses at published
optima) and synthetic residuals only; no new statistic of the real
residuals is computed or printed in gates mode (the identity gate
re-prints the published GB-9 numbers).
"""
import math, os, sys, time
import numpy as np
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)

MODE = sys.argv[1] if len(sys.argv) > 1 else 'gates'
assert MODE in ('gates', 'sky')

# ---------------- inherit the 10V header (round49_gb.py pattern) ----
SRCV = open(os.path.join(HERE, 'stage10v_strat.py'),
            encoding='utf-8').read()
cutv = SRCV.index("# ---------------- baselines (both modes)")
NSV = {'__name__': 'stage10v_trunc',
       '__file__': os.path.join(HERE, 'stage10v_strat.py')}
_argv = sys.argv
sys.argv = ['stage10v_strat.py', 'gates']
exec(compile(SRCV[:cutv], 'stage10v_trunc', 'exec'), NSV)
sys.argv = _argv

W78, LEGA, FLOW = NSV['W78'], NSV['LEGA'], NSV['FLOW']
FAMS, MEMBERS, ARC = NSV['FAMS'], NSV['MEMBERS'], NSV['ARC']
build_sub, fit_hier, fit_deep = (NSV['build_sub'], NSV['fit_hier'],
                                 NSV['fit_deep'])
nu_be = NSV['nu_be']
S_ML, U_PRIOR, LN10 = NSV['S_ML'], NSV['U_PRIOR'], NSV['LN10']
NAMEOF = NSV['NAMEOF']
sub_fields = NSV['sub_fields']
NS10U = NSV['NS']
SP, UD, UB, KPC = NS10U['SP'], NS10U['UD'], NS10U['UB'], NS10U['KPC']

from scipy.optimize import minimize_scalar as msc  # noqa: E402

L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

t0 = time.time()
RNG = np.random.default_rng(20260923)
TH_STAR = 20.0          # arcsec (prereg locked)
RD_FAC = 1.5            # x Rdisk (prereg locked, 7B/7C boundary)
RHO_AR = 0.66
SIG_PT = 0.10
N_PERM = 4000
N_SHIFT_SIZE = 2000
N_SHIFT_BR = 1000
N_WORLD = 1000
N_WIRE = 500
NMIN_S1, NMIN_CELL, NMIN_S2 = 8, 5, 15

# ---------------- baselines (verbatim round49_gb.py) ----------------
subP = build_sub(W78, LEGA)
b_plain = fit_hier(subP, nu_be, use_u=False)
b_loose = fit_hier(subP, nu_be, use_u=True, th0=list(b_plain.x)+[0.0])
arch = {'BE': fit_deep(subP, nu_be, th0=list(b_loose.x))}
for nm in ('p065', 'gm', 'boot'):
    arch[nm] = fit_deep(subP, FAMS[nm], th0=list(arch['BE'].x))
subF = build_sub(W78, FLOW)
fw = {m: fit_deep(subF, FAMS[m],
                  th0=[math.log10(ARC[m]), 1.0, 0.08, 0.0])
      for m in MEMBERS}
P(f"baselines  [{time.time()-t0:.0f}s]")

# ---------------- per-point theta + R_kpc (GB-9 parse, extended) ----
import glob as _glob
ROTPATH = {os.path.basename(p).replace('_rotmod.dat', ''): p
           for p in _glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                               recursive=True)}
PT = np.full(len(W78['lgobs']), np.nan)
RK = np.full(len(W78['lgobs']), np.nan)
for g in LEGA + FLOW:
    nm2 = NAMEOF[int(g)]
    D0 = SP['dist'][int(g)][0]
    ths, rks = [], []
    for l in open(ROTPATH[nm2]):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC
        gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        ths.append(R/D0*206.265)
        rks.append(R)
    PT[W78['gpts'][int(g)]] = ths
    RK[W78['gpts'][int(g)]] = rks

# ---------------- Rdisk from the SPARC master table -----------------
RDISK = {}
for l in open('data/sparc/SPARC_Lelli2016c.mrt', encoding='utf-8'):
    if len(l) < 66: continue
    nm2 = l[0:11].strip()
    try:
        rd = float(l[61:66])
    except ValueError:
        continue
    if nm2:
        RDISK[nm2] = rd

# ---------------- fields (sub_fields verbatim + y + nuisances) ------
def fields_ext(sub, nm, th):
    la0, f, s_int, u = th
    a0 = 10**la0
    nu = FAMS[nm]
    npt = len(sub['lg'])
    MU = np.empty(npt); RW = np.empty(npt); RES = np.empty(npt)
    Y = np.empty(npt)
    DML = np.empty(sub['n']); DV = np.empty(sub['n'])
    for k in range(sub['n']):
        pts = sub['gpts'][k]
        lg = sub['lg'][pts]
        gg = sub['gg'][pts]; gd = sub['gd'][pts]; gb = sub['gb'][pts]
        s2 = sub['s2'][pts]
        sv = sub['sv'][k]
        isu = sub['isuma_g'][k]
        se2 = s2 + s_int*s_int
        dml, dv = 0.0, 0.0
        for _ in range(200):
            gN = gg + f*math.exp(dml)*gd + gb
            r0 = lg - np.log10(gN*nu(gN/a0)) - u*isu
            w = 1.0/se2
            dvn = float(np.sum(w*r0)/(np.sum(w) + 1.0/sv**2))

            def od(dl):
                gN2 = gg + f*math.exp(dl)*gd + gb
                rr = lg - np.log10(gN2*nu(gN2/a0)) - dvn - u*isu
                return float(np.sum(rr*rr/se2) + dl*dl/(S_ML*S_ML))
            dmln = msc(od, bounds=(-0.7, 0.7), method='bounded').x
            if abs(dmln - dml) < 1e-12 and abs(dvn - dv) < 1e-12:
                dml, dv = dmln, dvn
                break
            dml, dv = dmln, dvn
        gN = gg + f*math.exp(dml)*gd + gb
        mu = np.log10(gN*nu(gN/a0)) + dv + u*isu
        MU[pts] = mu; RW[pts] = 1.0/se2; RES[pts] = lg - mu
        Y[pts] = gN/a0
        DML[k] = dml; DV[k] = dv
    return MU, RW, RES, Y, DML, DV

def lag1_of(sub, res):
    num = den = 0.0
    for k in range(sub['n']):
        r = res[sub['gpts'][k]]
        if len(r) < 3: continue
        r = r - r.mean()
        num += float(np.sum(r[:-1]*r[1:]))
        den += float(np.sum(r*r))
    return num/den if den > 0 else float('nan')

# ---------------- leg design objects --------------------------------
def leg_design(sub, seq, opt):
    """A2 best form (GB-9 rule), fields at its optimum, per-galaxy
    point arrays in radius order (rotmod order = increasing R)."""
    A2 = 'BE' if opt['BE'].fun <= opt['boot'].fun else 'boot'
    MU, RW, RES, Y, DML, DV = fields_ext(sub, A2, list(opt[A2].x))
    gal = []
    for k, g in enumerate(seq):
        pts = sub['gpts'][k]
        th_g = PT[W78['gpts'][int(g)]]
        rk_g = RK[W78['gpts'][int(g)]]
        rd = RDISK.get(NAMEOF[int(g)], np.nan)
        gal.append(dict(k=k, g=int(g), name=NAMEOF[int(g)],
                        pts=pts, th=np.asarray(th_g),
                        rk=np.asarray(rk_g), rd=rd,
                        D=SP['dist'][int(g)][0],
                        y=Y[pts], res=RES[pts]))
    return A2, gal, RES, Y

A2_A, GAL_A, RES_A, Y_A = leg_design(subP, LEGA, arch)
A2_F, GAL_F, RES_F, Y_F = leg_design(subF, FLOW, fw)
P(f"designs: anch best form {A2_A}, flow best form {A2_F}  "
  f"[{time.time()-t0:.0f}s]")

# ---------------- the statistics (residual vector -> everything) ----
def contrast_plain(gal, res_by_gal):
    diffs, names = [], []
    for gd_, r in zip(gal, res_by_gal):
        inn = r[gd_['th'] < TH_STAR]; out = r[gd_['th'] >= TH_STAR]
        if len(inn) and len(out):
            diffs.append(float(inn.mean() - out.mean()))
            names.append(gd_['name'])
    d = np.array(diffs)
    if len(d) < 2:
        return d, names, float('nan'), float('nan')
    return d, names, float(d.mean()), float(d.std(ddof=1)/math.sqrt(len(d)))

def contrast_matched(gal, res_by_gal):
    diffs = []
    for gd_, r in zip(gal, res_by_gal):
        mi = gd_['th'] < TH_STAR; mo = ~mi
        if not (mi.any() and mo.any()):
            continue
        yi, yo = gd_['y'][mi], gd_['y'][mo]
        lo = max(yi.min(), yo.min()); hi = min(yi.max(), yo.max())
        if not (hi > lo):
            continue
        ki = mi & (gd_['y'] >= lo) & (gd_['y'] <= hi)
        ko = mo & (gd_['y'] >= lo) & (gd_['y'] <= hi)
        if ki.any() and ko.any():
            diffs.append(float(r[ki].mean() - r[ko].mean()))
    d = np.array(diffs)
    if len(d) < 2:
        return d, float('nan'), float('nan')
    return d, float(d.mean()), float(d.std(ddof=1)/math.sqrt(len(d)))

def detrend_contrast(gal, res_by_gal, deg):
    ly = np.concatenate([np.log10(gd_['y']) for gd_ in gal])
    rr = np.concatenate(res_by_gal)
    cf = np.polyfit(ly, rr, deg)
    out = []
    for gd_, r in zip(gal, res_by_gal):
        out.append(r - np.polyval(cf, np.log10(gd_['y'])))
    return contrast_plain(gal, out)

def spearman(a, b):
    ra, rb = rankdata(a), rankdata(b)
    ra -= ra.mean(); rb -= rb.mean()
    den = math.sqrt(float(np.sum(ra*ra)*np.sum(rb*rb)))
    return float(np.sum(ra*rb)/den) if den > 0 else float('nan')

def s2_tests(gal, res_by_gal, rng, nperm=N_PERM):
    d, names, _, _ = contrast_plain(gal, res_by_gal)
    sel = []
    for gd_, r in zip(gal, res_by_gal):
        inn = r[gd_['th'] < TH_STAR]; out = r[gd_['th'] >= TH_STAR]
        if len(inn) and len(out):
            sel.append(gd_)
    thmin = np.array([float(g_['th'].min()) for g_ in sel])
    dist = np.array([g_['D'] for g_ in sel])
    outp = {}
    for tag, cov, sgn in (('S2a_thmin', thmin, +1),
                          ('S2b_dist', dist, -1)):
        if len(d) < NMIN_S2:
            outp[tag] = (float('nan'), float('nan'), sgn, len(d))
            continue
        rho = spearman(d, cov)
        cnt = 0
        for _ in range(nperm):
            pr = spearman(d, rng.permutation(cov))
            if abs(pr) >= abs(rho):
                cnt += 1
        outp[tag] = (rho, (cnt + 1)/(nperm + 1), sgn, len(d))
    return outp

def s3_cells(gal, res_by_gal):
    dio, doi, n_io, n_oi = [], [], 0, 0
    for gd_, r in zip(gal, res_by_gal):
        if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
            continue
        inn_th = gd_['th'] < TH_STAR
        inn_rd = gd_['rk'] < RD_FAC*gd_['rd']
        Io = inn_th & ~inn_rd
        Oi = ~inn_th & inn_rd
        Oo = ~inn_th & ~inn_rd
        if not Oo.any():
            continue
        base = float(r[Oo].mean())
        if Io.any():
            dio.append(float(r[Io].mean()) - base)
        if Oi.any():
            doi.append(float(r[Oi].mean()) - base)
    def pack(v):
        v = np.array(v)
        if len(v) < 2:
            return v, float('nan'), float('nan')
        return v, float(v.mean()), float(v.std(ddof=1)/math.sqrt(len(v)))
    return pack(dio), pack(doi)

def s3_shift_p(gal, res_by_gal, obs_io, obs_oi, rng,
               nshift=N_SHIFT_SIZE):
    cio = coi = nio = noi = 0
    for _ in range(nshift):
        shifted = []
        for gd_, r in zip(gal, res_by_gal):
            s = int(rng.integers(0, len(r))) if len(r) else 0
            shifted.append(np.roll(r, s))
        (dio, mio, _), (doi, moi, _) = s3_cells(gal, shifted)
        if np.isfinite(mio) and np.isfinite(obs_io):
            nio += 1
            if mio <= obs_io:
                cio += 1
        if np.isfinite(moi) and np.isfinite(obs_oi):
            noi += 1
            if moi <= obs_oi:
                coi += 1
    pio = (cio + 1)/(nio + 1) if nio else float('nan')
    poi = (coi + 1)/(noi + 1) if noi else float('nan')
    return pio, poi

def grammar(gal, res_by_gal, rng, nperm=N_PERM, nshift=N_SHIFT_SIZE,
            want_p=True):
    """Full branch evaluation (prereg section 4). Returns branch +
    the component record."""
    dm, mm, sm = contrast_matched(gal, res_by_gal)
    rec = dict(s1_n=len(dm), s1_mean=mm, s1_se=sm)
    s1_pop = len(dm) >= NMIN_S1
    survives = s1_pop and np.isfinite(mm) and mm <= -2*sm
    collapses = s1_pop and np.isfinite(mm) and abs(mm) < 1*sm
    (dio, mio, sio), (doi, moi, soi) = s3_cells(gal, res_by_gal)
    rec.update(io_n=len(dio), io_mean=mio, io_se=sio,
               oi_n=len(doi), oi_mean=moi, oi_se=soi)
    io_pop = len(dio) >= NMIN_CELL
    oi_pop = len(doi) >= NMIN_CELL
    sig_ang = (io_pop and oi_pop and np.isfinite(mio)
               and mio <= -2*sio and abs(moi) < 1*soi)
    sig_dsk = (io_pop and oi_pop and np.isfinite(moi)
               and moi <= -2*soi and abs(mio) < 1*sio)
    s2 = s2_tests(gal, res_by_gal, rng, nperm)
    rec['s2'] = s2
    ps = sorted([(s2['S2a_thmin'][1], 'S2a_thmin'),
                 (s2['S2b_dist'][1], 'S2b_dist')])
    s2_fire = False
    if np.isfinite(ps[0][0]):
        holm = [(ps[0][0]*2, ps[0][1]), (ps[1][0]*1, ps[1][1])]
        for padj, tag in holm:
            rho, _, sgn, _ = s2[tag]
            if padj <= 0.05 and np.isfinite(rho) and rho*sgn > 0:
                s2_fire = True
            else:
                break
    rec.update(sig_ang=sig_ang, sig_dsk=sig_dsk, s2_fire=s2_fire,
               s1_pop=s1_pop, io_pop=io_pop, oi_pop=oi_pop)
    if want_p:
        rec['io_p'], rec['oi_p'] = s3_shift_p(
            gal, res_by_gal, mio, moi, rng, nshift)
    if s1_pop and collapses:
        br = 'B1'
    elif survives and sig_ang and s2_fire:
        br = 'B2'
    elif survives and sig_dsk and not s2_fire:
        br = 'B3'
    else:
        br = 'B4'
    rec['branch'] = br
    return br, rec

# ---------------- synthetic worlds (design real, residuals not) -----
def ar1(rng, n):
    e = np.empty(n)
    e[0] = rng.normal(0, SIG_PT)
    for i in range(1, n):
        e[i] = RHO_AR*e[i-1] + math.sqrt(1-RHO_AR**2)*rng.normal(0, SIG_PT)
    return e

def world_noise(gal, rng):
    return [ar1(rng, len(gd_['res'])) for gd_ in gal]

def yslope_calibration(gal):
    comps = []
    for gd_ in gal:
        mi = gd_['th'] < TH_STAR; mo = ~mi
        if mi.any() and mo.any():
            ly = np.log10(gd_['y'])
            comps.append(float(ly[mi].mean() - ly[mo].mean()))
    C = float(np.mean(comps))
    return -0.0857/C, C, len(comps)

B_Y, C_Y, N_CAL = yslope_calibration(GAL_A)

def world_y(gal, rng):
    return [B_Y*np.log10(gd_['y']) + e
            for gd_, e in zip(gal, world_noise(gal, rng))]

def world_r(gal, rng):
    return [(-0.0857)*(gd_['th'] < TH_STAR) + e
            for gd_, e in zip(gal, world_noise(gal, rng))]

# ====================================================================
if MODE == 'gates':
    P("")
    P("=" * 68)
    P("STAGE 10Y GATES (pre-sky; PREREG-CONTRAST-10Y.md)")
    P("")

    # -- G10Y-1 identity: the GB-9 regression (published targets) ----
    P("-- G10Y-1 identity (GB-9 regression) --")
    ok1 = True
    for legnm, gal, sub, opt in (('anch', GAL_A, subP, arch),
                                 ('flow', GAL_F, subF, fw)):
        A2 = 'BE' if opt['BE'].fun <= opt['boot'].fun else 'boot'
        fld = sub_fields(sub, A2, list(opt[A2].x))
        res_inh = fld[2]
        res_by_gal = [gd_['res'] for gd_ in gal]
        dmax = max(float(np.max(np.abs(res_inh[gd_['pts']] - r)))
                   for gd_, r in zip(gal, res_by_gal))
        d, names, m, s = contrast_plain(gal, res_by_gal)
        tgt = {'anch': (-0.0857, 0.0301, 20),
               'flow': (-0.0055, 0.0176, 45)}[legnm]
        p1 = (abs(m - tgt[0]) <= 1.5e-4 and abs(s - tgt[1]) <= 1.5e-4
              and len(d) == tgt[2] and dmax == 0.0)
        ok1 &= p1
        P(f"  {legnm}: mean {m:+.4f} +/- {s:.4f} over {len(d)} gal "
          f"(target {tgt[0]:+.4f} +/- {tgt[1]:.4f} / {tgt[2]}); "
          f"fields_ext vs inherited max|d| = {dmax:.1e}  "
          f"{'PASS' if p1 else 'FAIL'}")
    P(f"  G10Y-1 {'PASS' if ok1 else 'FAIL'}")

    # -- G10Y-2 Rdisk census + cell-population feasibility -----------
    P("")
    P("-- G10Y-2 Rdisk census --")
    miss = [NAMEOF[int(g)] for g in LEGA + FLOW
            if RDISK.get(NAMEOF[int(g)], 0) <= 0]
    ok2 = (len(miss) == 0)
    P(f"  matched with Rdisk > 0: {len(LEGA)+len(FLOW)-len(miss)}/"
      f"{len(LEGA)+len(FLOW)}; missing: {miss if miss else 'none'}")
    for nm2 in ('IC2574', 'DDO154', 'NGC2403'):
        P(f"  spot: Rdisk({nm2}) = {RDISK.get(nm2, float('nan')):.2f} kpc")
    ncells = dict(Io=0, Oi=0)
    for gd_ in GAL_A:
        if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
            continue
        inn_th = gd_['th'] < TH_STAR
        inn_rd = gd_['rk'] < RD_FAC*gd_['rd']
        Oo = ~inn_th & ~inn_rd
        if not Oo.any():
            continue
        if (inn_th & ~inn_rd).any(): ncells['Io'] += 1
        if (~inn_th & inn_rd).any(): ncells['Oi'] += 1
    P(f"  anchored cell-population census (design): Io {ncells['Io']} "
      f"gal, Oi {ncells['Oi']} gal (N_min {NMIN_CELL} each)")
    P(f"  WORLD-Y calibration: b = {B_Y:+.4f} per dex(y) "
      f"(composition term C = {C_Y:+.4f} over {N_CAL} gal)")
    P(f"  G10Y-2 {'PASS' if ok2 else 'FAIL'}")

    # -- G10Y-3 null machinery: size + branch-rate tables ------------
    P("")
    P("-- G10Y-3 size (noise-only, 2000 draws; inner perms 500 for")
    P("   the null tables -- the sky grammar uses the full 4000) --")
    rng = np.random.default_rng(20260923)
    fires = dict(sig_ang=0, sig_dsk=0, s2=0, s1_surv=0, s1_coll=0)
    nsz = 2000
    for _ in range(nsz):
        rb = world_noise(GAL_A, rng)
        br, rec = grammar(GAL_A, rb, rng, nperm=500, nshift=0,
                          want_p=False)
        fires['sig_ang'] += rec['sig_ang']
        fires['sig_dsk'] += rec['sig_dsk']
        fires['s2'] += rec['s2_fire']
        dm, mm, sm = rec['s1_n'], rec['s1_mean'], rec['s1_se']
        if rec['s1_pop'] and np.isfinite(mm):
            fires['s1_surv'] += (mm <= -2*sm)
            fires['s1_coll'] += (abs(mm) < sm)
    for k in fires:
        fires[k] /= nsz
    okA = 0.01 <= fires['sig_ang'] <= 0.10
    okD = 0.01 <= fires['sig_dsk'] <= 0.10
    okS = 0.01 <= fires['s2'] <= 0.10
    P(f"  P(sig_ang) = {fires['sig_ang']:.3f}  "
      f"{'PASS' if okA else 'FAIL'} [0.01, 0.10]")
    P(f"  P(sig_dsk) = {fires['sig_dsk']:.3f}  "
      f"{'PASS' if okD else 'FAIL'} [0.01, 0.10]")
    P(f"  P(S2 Holm+sign) = {fires['s2']:.3f}  "
      f"{'PASS' if okS else 'FAIL'} [0.01, 0.10]")
    P(f"  (context: P(S1 survive|noise) = {fires['s1_surv']:.3f}, "
      f"P(S1 collapse|noise) = {fires['s1_coll']:.3f})")
    ok3a = okA and okD and okS

    P("")
    P("-- G10Y-3 branch rates (1000 draws/world; trap #32/#33) --")
    tables = {}
    for wnm, wfun in (('WORLD-Y', world_y), ('WORLD-R', world_r)):
        cnt = dict(B1=0, B2=0, B3=0, B4=0, B4_s1surv=0)
        for _ in range(N_WORLD):
            rb = wfun(GAL_A, rng)
            br, rec = grammar(GAL_A, rb, rng, nperm=500, nshift=0,
                              want_p=False)
            cnt[br] += 1
            if br == 'B4' and rec['s1_pop'] and \
               np.isfinite(rec['s1_mean']) and \
               rec['s1_mean'] <= -rec['s1_se']:
                cnt['B4_s1surv'] += 1
        for k in cnt:
            cnt[k] /= N_WORLD
        tables[wnm] = cnt
        P(f"  {wnm}: B1 {cnt['B1']:.3f}  B2 {cnt['B2']:.3f}  "
          f"B3 {cnt['B3']:.3f}  B4 {cnt['B4']:.3f}  "
          f"(B4 w/ S1<=-1SE: {cnt['B4_s1surv']:.3f})")
    pw1 = tables['WORLD-Y']['B1'] >= 0.60
    keep_r = (tables['WORLD-R']['B2'] + tables['WORLD-R']['B3']
              + tables['WORLD-R']['B4_s1surv'])
    pw2 = keep_r >= 0.60
    P(f"  POWER: P(B1|WORLD-Y) = {tables['WORLD-Y']['B1']:.3f} "
      f"{'PASS' if pw1 else 'FAIL'} (>= 0.60); "
      f"P(radius kept|WORLD-R) = {keep_r:.3f} "
      f"{'PASS' if pw2 else 'FAIL'} (>= 0.60)")
    ok3 = ok3a and pw1 and pw2

    # -- G10Y-4 y-match wiring ----------------------------------------
    P("")
    P("-- G10Y-4 y-match wiring (500 draws/world) --")
    czero = ckeep = 0
    for _ in range(N_WIRE):
        dm, mm, sm = contrast_matched(GAL_A, world_y(GAL_A, rng))
        if len(dm) >= NMIN_S1 and np.isfinite(mm) and abs(mm) < sm:
            czero += 1
        dm, mm, sm = contrast_matched(GAL_A, world_r(GAL_A, rng))
        if len(dm) >= NMIN_S1 and np.isfinite(mm) and mm <= -sm:
            ckeep += 1
    czero /= N_WIRE; ckeep /= N_WIRE
    ok4 = (czero >= 0.80) and (ckeep >= 0.60)
    P(f"  P(matched ~0 | WORLD-Y) = {czero:.3f} "
      f"{'PASS' if czero >= 0.80 else 'FAIL'} (>= 0.80)")
    P(f"  P(matched <= -1SE | WORLD-R) = {ckeep:.3f} "
      f"{'PASS' if ckeep >= 0.60 else 'FAIL'} (>= 0.60)")

    P("")
    allok = ok1 and ok2 and ok3 and ok4
    P(f"GATES {'ALL PASS' if allok else 'FAIL'}; "
      f"wall-clock {(time.time()-t0)/60:.1f} min")
    with open('data/stage10y_gates.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(L) + "\n")
    print("\nwrote: data/stage10y_gates.txt")

# ====================================================================
else:  # ---------------------------- sky ----------------------------
    P("")
    P("=" * 68)
    P("STAGE 10Y SKY READ (post-prereg-commit)")
    P("")
    rng = np.random.default_rng(20260923 + 1)

    for legnm, gal, sub in (('ANCHORED', GAL_A, subP),
                            ('FLOW', GAL_F, subF)):
        res_by_gal = [gd_['res'] for gd_ in gal]
        P(f"---- {legnm} leg (best form "
          f"{A2_A if legnm == 'ANCHORED' else A2_F}) ----")
        d, names, m, s = contrast_plain(gal, res_by_gal)
        P(f"  S0 plain contrast: {m:+.4f} +/- {s:.4f} over {len(d)} gal")
        dm, mm, sm = contrast_matched(gal, res_by_gal)
        P(f"  S1 y-matched:      {mm:+.4f} +/- {sm:.4f} over "
          f"{len(dm)} gal (N_min {NMIN_S1})")
        for deg, tag in ((1, 'lin'), (2, 'quad'), (3, 'cub')):
            _, _, md, sd_ = detrend_contrast(gal, res_by_gal, deg)
            P(f"  S1 co-read detrend-{tag}: {md:+.4f} +/- {sd_:.4f} "
              f"(biased toward collapse; no branch weight)")
        s2 = s2_tests(gal, res_by_gal, rng)
        for tag in ('S2a_thmin', 'S2b_dist'):
            rho, p, sgn, n = s2[tag]
            P(f"  {tag}: rho = {rho:+.3f}, perm p = {p:.4f} "
            f"(smearing sign {'+' if sgn > 0 else '-'}; n = {n})")
        (dio, mio, sio), (doi, moi, soi) = s3_cells(gal, res_by_gal)
        P(f"  S3 d_Io (inner-angle only vs Oo): {mio:+.4f} +/- "
          f"{sio:.4f} over {len(dio)} gal")
        P(f"  S3 d_Oi (inner-disk only vs Oo):  {moi:+.4f} +/- "
          f"{soi:.4f} over {len(doi)} gal")
        pio, poi = s3_shift_p(gal, res_by_gal, mio, moi, rng)
        P(f"  S3 circular-shift p: Io {pio:.4f}, Oi {poi:.4f} "
          f"(one-sided toward depression; {N_SHIFT_SIZE} draws)")
        legres = RES_A if legnm == 'ANCHORED' else RES_F
        P(f"  rho_lag1 = {lag1_of(sub, legres):.3f}")
        P("")

    # banded profile (descriptive, anchored)
    P("---- S4 banded residual profile (ANCHORED; descriptive) ----")
    TB = [(0, 10), (10, 20), (20, 40), (40, 1e9)]
    RB = [(0, 1.0), (1.0, 1.5), (1.5, 3.0), (3.0, 1e9)]
    yall = np.concatenate([gd_['y'] for gd_ in GAL_A])
    yt = np.quantile(np.log10(yall), [1/3, 2/3])
    for it, (t1, t2) in enumerate(TB):
        for ir, (r1, r2) in enumerate(RB):
            for iy in range(3):
                vals = []
                for gd_ in GAL_A:
                    if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
                        continue
                    ly = np.log10(gd_['y'])
                    msk = ((gd_['th'] >= t1) & (gd_['th'] < t2)
                           & (gd_['rk']/gd_['rd'] >= r1)
                           & (gd_['rk']/gd_['rd'] < r2))
                    if iy == 0:
                        msk &= ly < yt[0]
                    elif iy == 1:
                        msk &= (ly >= yt[0]) & (ly < yt[1])
                    else:
                        msk &= ly >= yt[1]
                    if msk.any():
                        vals.extend(gd_['res'][msk])
                if len(vals) >= 10:
                    v = np.array(vals)
                    P(f"  th[{t1:g},{t2:g}) R/Rd[{r1:g},{r2:g}) "
                      f"y-ter{iy}: {v.mean():+.4f} +/- "
                      f"{v.std(ddof=1)/math.sqrt(len(v)):.4f} "
                      f"(n = {len(v)})")
    P("")

    # the letter
    P("---- LETTER (prereg section 4 grammar; ANCHORED) ----")
    rngL = np.random.default_rng(20260923 + 2)
    br, rec = grammar(GAL_A, [gd_['res'] for gd_ in GAL_A], rngL)
    P(f"  S1: n {rec['s1_n']}, mean {rec['s1_mean']:+.4f} +/- "
      f"{rec['s1_se']:.4f}; populated {rec['s1_pop']}")
    P(f"  S3: Io {rec['io_mean']:+.4f} +/- {rec['io_se']:.4f} "
      f"(n {rec['io_n']}, p {rec.get('io_p', float('nan')):.4f}); "
      f"Oi {rec['oi_mean']:+.4f} +/- {rec['oi_se']:.4f} "
      f"(n {rec['oi_n']}, p {rec.get('oi_p', float('nan')):.4f})")
    P(f"  signatures: angular {rec['sig_ang']}, disk {rec['sig_dsk']}; "
      f"S2 fire {rec['s2_fire']}")
    P(f"  BRANCH: {br}")
    P("")
    P(f"carrier note: innermost angles UGC03580 "
      f"{min(g_['th'].min() for g_ in GAL_F if g_['name'] == 'UGC03580'):.1f} "
      f"arcsec (flow), NGC5371 "
      f"{min(g_['th'].min() for g_ in GAL_F if g_['name'] == 'NGC5371'):.1f} "
      f"arcsec (flow)")
    P(f"sky complete; wall-clock {(time.time()-t0)/60:.1f} min")
    with open('data/stage10y_skyread.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(L) + "\n")
    print("\nwrote: data/stage10y_skyread.txt")
