"""STAGE 10Y -- THE RESIDUAL-CONTRAST ORGANIZER (A1 form).
Pre-registration: PREREG-CONTRAST-10Y.md + amendment A1 (gates run 1
preserved as data/stage10y_gates_r1.txt; A1 committed pre-quote).

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
TH_STAR = 20.0          # arcsec (prereg locked)
RD_FAC = 1.5            # x Rdisk (prereg locked, 7B/7C boundary)
RHO_AR = 0.66
SIG_PT = 0.10
N_PERM = 4000
N_SHIFT_SIZE = 2000
N_WORLD = 1000
N_WIRE = 500
NMIN_S1, NMIN_CELL, NMIN_S2 = 8, 5, 15
NMIN_BAND_PTS, NMIN_BANDS = 5, 3

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

# ---------------- Rdisk (A1-i: token parse; index 11) ---------------
RDISK = {}
for l in open('data/sparc/SPARC_Lelli2016c.mrt', encoding='utf-8'):
    t = l.split()
    if len(t) < 18:
        continue
    try:
        int(t[1]); rd = float(t[11])
    except ValueError:
        continue
    RDISK[t[0]] = rd

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
    return A2, gal, RES

A2_A, GAL_A, RES_A = leg_design(subP, LEGA, arch)
A2_F, GAL_F, RES_F = leg_design(subF, FLOW, fw)
P(f"designs: anch best form {A2_A}, flow best form {A2_F}  "
  f"[{time.time()-t0:.0f}s]")

def qual_and_edges(gal):
    """Both-sides galaxies + y-band quintile edges (DESIGN, fixed)."""
    qual = [j for j, gd_ in enumerate(gal)
            if (gd_['th'] < TH_STAR).any()
            and (gd_['th'] >= TH_STAR).any()]
    ly = np.concatenate([np.log10(gal[j]['y']) for j in qual])
    e = np.quantile(ly, [0, .2, .4, .6, .8, 1.0])
    e[0] -= 1e-9; e[-1] += 1e-9
    return qual, e

QUAL_A, EDGES_A = qual_and_edges(GAL_A)
QUAL_F, EDGES_F = qual_and_edges(GAL_F)

# ---------------- the statistics ------------------------------------
def contrast_plain(gal, res_by_gal):
    diffs = []
    for gd_, r in zip(gal, res_by_gal):
        inn = r[gd_['th'] < TH_STAR]; out = r[gd_['th'] >= TH_STAR]
        if len(inn) and len(out):
            diffs.append(float(inn.mean() - out.mean()))
    d = np.array(diffs)
    if len(d) < 2:
        return d, float('nan'), float('nan')
    return d, float(d.mean()), float(d.std(ddof=1)/math.sqrt(len(d)))

def banded_value(gal, res_by_gal, edges, qual, skip=None):
    """A1-ii: count-weighted mean of per-band inner-minus-outer
    contrasts over the pooled both-sides galaxies."""
    nb = len(edges) - 1
    si = np.zeros(nb); ni = np.zeros(nb)
    so = np.zeros(nb); no = np.zeros(nb)
    for j in qual:
        if skip is not None and j == skip:
            continue
        gd_, r = gal[j], res_by_gal[j]
        ly = np.log10(gd_['y'])
        bi = np.clip(np.searchsorted(edges, ly, side='right') - 1,
                     0, nb - 1)
        mi = gd_['th'] < TH_STAR
        np.add.at(si, bi[mi], r[mi]); np.add.at(ni, bi[mi], 1)
        np.add.at(so, bi[~mi], r[~mi]); np.add.at(no, bi[~mi], 1)
    use = (ni >= NMIN_BAND_PTS) & (no >= NMIN_BAND_PTS)
    if use.sum() == 0:
        return float('nan'), 0, np.full(nb, np.nan), ni, no
    cb = si[use]/ni[use] - so[use]/no[use]
    wb = ni[use]*no[use]/(ni[use] + no[use])
    val = float(np.sum(wb*cb)/np.sum(wb))
    cb_full = np.full(nb, np.nan)
    cb_full[use] = cb
    return val, int(use.sum()), cb_full, ni, no

def s1_banded(gal, res_by_gal, edges, qual):
    val, nbands, cb, ni, no = banded_value(gal, res_by_gal, edges, qual)
    if not np.isfinite(val) or len(qual) < 2:
        return float('nan'), float('nan'), nbands, len(qual), cb, ni, no
    jk = []
    for j in qual:
        v, nb_, _, _, _ = banded_value(gal, res_by_gal, edges, qual,
                                       skip=j)
        if np.isfinite(v):
            jk.append(v)
    jk = np.array(jk)
    n = len(jk)
    se = math.sqrt((n - 1)/n*float(np.sum((jk - jk.mean())**2))) \
        if n >= 2 else float('nan')
    return val, se, nbands, len(qual), cb, ni, no

DY_MATCH, MINW, NMIN_MPT = 0.15, 5, 15

def build_matcher(gal, quse, inner_of, outer_of):
    """A2-i: precomputed nearest-neighbor y-match design. inner_of/
    outer_of: fn(gd_) -> point mask. Windows are design (y fixed)."""
    wins = []
    B_j, B_i, B_ly = [], [], []
    for j in quse:
        gd_ = gal[j]
        ly = np.log10(gd_['y'])
        for i in np.where(outer_of(gd_))[0]:
            B_j.append(j); B_i.append(i); B_ly.append(ly[i])
    B_j = np.array(B_j, dtype=int)
    B_i = np.array(B_i, dtype=int)
    B_ly = np.array(B_ly)
    for j in quse:
        gd_ = gal[j]
        ly = np.log10(gd_['y'])
        for i in np.where(inner_of(gd_))[0]:
            k = np.where(np.abs(B_ly - ly[i]) <= DY_MATCH)[0]
            wins.append((j, int(i), k))
    return dict(wins=wins, B_j=B_j, B_i=B_i)

def matcher_census(M):
    q = [w for w in M['wins'] if len(w[2]) >= MINW]
    return len(M['wins']), len(q), len({w[0] for w in q})

def s1_matched(res_by_gal, M):
    """Value + leave-one-galaxy-out jackknife SE (A2-i)."""
    Bres = np.array([res_by_gal[j][i]
                     for j, i in zip(M['B_j'], M['B_i'])])
    B_j = M['B_j']

    def value(skip=None):
        vals, gals = [], set()
        for (j, i, k) in M['wins']:
            if skip is not None and j == skip:
                continue
            kk = k if skip is None else k[B_j[k] != skip]
            if len(kk) < MINW:
                continue
            vals.append(res_by_gal[j][i] - float(Bres[kk].mean()))
            gals.add(j)
        if not vals:
            return float('nan'), 0, set()
        return float(np.mean(vals)), len(vals), gals

    val, m, gals = value()
    if not np.isfinite(val) or len(gals) < 2:
        return val, float('nan'), m, len(gals)
    jk = []
    for j in sorted(gals):
        v, _, _ = value(skip=j)
        if np.isfinite(v):
            jk.append(v)
    jk = np.array(jk)
    n = len(jk)
    se = math.sqrt((n - 1)/n*float(np.sum((jk - jk.mean())**2))) \
        if n >= 2 else float('nan')
    return val, se, m, len(gals)

M_S1_A = build_matcher(GAL_A, QUAL_A,
                       lambda gd_: gd_['th'] < TH_STAR,
                       lambda gd_: gd_['th'] >= TH_STAR)
M_S1_F = build_matcher(GAL_F, QUAL_F,
                       lambda gd_: gd_['th'] < TH_STAR,
                       lambda gd_: gd_['th'] >= TH_STAR)

def _oi_mask(gd_):
    if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
        return np.zeros(len(gd_['th']), dtype=bool)
    return (gd_['th'] >= TH_STAR) & (gd_['rk'] < RD_FAC*gd_['rd'])

def _oo_mask(gd_):
    if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
        return np.zeros(len(gd_['th']), dtype=bool)
    return (gd_['th'] >= TH_STAR) & (gd_['rk'] >= RD_FAC*gd_['rd'])

M_OI_A = build_matcher(GAL_A, list(range(len(GAL_A))),
                       _oi_mask, _oo_mask)

def overlap_census(gal, res_by_gal):
    """The retired per-galaxy y-overlap match (co-read census)."""
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
    return diffs

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
    d, sel = [], []
    for gd_, r in zip(gal, res_by_gal):
        inn = r[gd_['th'] < TH_STAR]; out = r[gd_['th'] >= TH_STAR]
        if len(inn) and len(out):
            d.append(float(inn.mean() - out.mean()))
            sel.append(gd_)
    d = np.array(d)
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
    """A1-iii: d_Ii and d_Oi vs the Oo baseline; Io census only."""
    dii, doi, n_io = [], [], 0
    for gd_, r in zip(gal, res_by_gal):
        if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
            continue
        inn_th = gd_['th'] < TH_STAR
        inn_rd = gd_['rk'] < RD_FAC*gd_['rd']
        Ii = inn_th & inn_rd
        Io = inn_th & ~inn_rd
        Oi = ~inn_th & inn_rd
        Oo = ~inn_th & ~inn_rd
        if not Oo.any():
            continue
        if Io.any():
            n_io += 1
        base = float(r[Oo].mean())
        if Ii.any():
            dii.append(float(r[Ii].mean()) - base)
        if Oi.any():
            doi.append(float(r[Oi].mean()) - base)
    def pack(v):
        v = np.array(v)
        if len(v) < 2:
            return v, float('nan'), float('nan')
        return v, float(v.mean()), float(v.std(ddof=1)/math.sqrt(len(v)))
    return pack(dii), pack(doi), n_io

def s3_shift_p(gal, res_by_gal, obs_ii, obs_oi, rng,
               nshift=N_SHIFT_SIZE):
    cii = coi = nii = noi = 0
    for _ in range(nshift):
        shifted = []
        for gd_, r in zip(gal, res_by_gal):
            s = int(rng.integers(0, len(r))) if len(r) else 0
            shifted.append(np.roll(r, s))
        (dii, mii, _), (doi, moi, _), _ = s3_cells(gal, shifted)
        if np.isfinite(mii) and np.isfinite(obs_ii):
            nii += 1
            if mii <= obs_ii:
                cii += 1
        if np.isfinite(moi) and np.isfinite(obs_oi):
            noi += 1
            if moi <= obs_oi:
                coi += 1
    pii = (cii + 1)/(nii + 1) if nii else float('nan')
    poi = (coi + 1)/(noi + 1) if noi else float('nan')
    return pii, poi

def grammar(gal, res_by_gal, rng, M, nperm=N_PERM, want_p=False):
    mm, sm, mpts, ng = s1_matched(res_by_gal, M)
    rec = dict(s1_mean=mm, s1_se=sm, s1_mpts=mpts, s1_ngal=ng)
    s1_pop = (mpts >= NMIN_MPT) and (ng >= NMIN_S1) \
        and np.isfinite(mm) and np.isfinite(sm)
    survives = s1_pop and mm <= -2*sm
    collapses = s1_pop and abs(mm) < 1*sm
    (dii, mii, sii), (doi, moi, soi), n_io = s3_cells(gal, res_by_gal)
    rec.update(ii_n=len(dii), ii_mean=mii, ii_se=sii,
               oi_n=len(doi), oi_mean=moi, oi_se=soi, io_gal=n_io)
    ii_pop = len(dii) >= NMIN_CELL
    oi_pop = len(doi) >= NMIN_CELL
    ii_deep = ii_pop and np.isfinite(mii) and mii <= -2*sii
    oi_deep = oi_pop and np.isfinite(moi) and moi <= -2*soi
    oi_flat = oi_pop and np.isfinite(moi) and abs(moi) < 1*soi
    sig_ang = ii_deep and oi_pop and oi_flat
    sig_dsk = ii_deep and oi_deep
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
               s1_pop=s1_pop, ii_deep=ii_deep, oi_deep=oi_deep,
               oi_flat=oi_flat, survives=survives, collapses=collapses)
    if want_p:
        rec['ii_p'], rec['oi_p'] = s3_shift_p(
            gal, res_by_gal, mii, moi, rng)
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
    P("STAGE 10Y GATES (pre-sky; PREREG-CONTRAST-10Y.md + A1)")
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
        d, m, s = contrast_plain(gal, res_by_gal)
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

    # -- G10Y-2 Rdisk census + design feasibility (A1-i/iii) ---------
    P("")
    P("-- G10Y-2 Rdisk census (A1-i token parse) --")
    miss = [NAMEOF[int(g)] for g in LEGA + FLOW
            if RDISK.get(NAMEOF[int(g)], 0) <= 0]
    ok2 = (len(miss) == 0)
    P(f"  matched with Rdisk > 0: {len(LEGA)+len(FLOW)-len(miss)}/"
      f"{len(LEGA)+len(FLOW)}; missing: {miss if miss else 'none'}")
    for nm2, truth in (('IC2574', 2.78), ('DDO154', 0.37),
                       ('NGC2403', 1.39)):
        got = RDISK.get(nm2, float('nan'))
        oks = abs(got - truth) < 0.005
        ok2 &= oks
        P(f"  spot: Rdisk({nm2}) = {got:.2f} kpc (truth {truth:.2f}) "
          f"{'PASS' if oks else 'FAIL'}")
    cens = dict(Ii=0, Io=0, Oi=0, Oo=0)
    for gd_ in GAL_A:
        if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
            continue
        inn_th = gd_['th'] < TH_STAR
        inn_rd = gd_['rk'] < RD_FAC*gd_['rd']
        if not (~inn_th & ~inn_rd).any():
            continue
        if (inn_th & inn_rd).any(): cens['Ii'] += 1
        if (inn_th & ~inn_rd).any(): cens['Io'] += 1
        if (~inn_th & inn_rd).any(): cens['Oi'] += 1
        cens['Oo'] += 1
    P(f"  anchored cell census (design, galaxies with Oo): "
      f"Ii {cens['Ii']}, Io {cens['Io']}, Oi {cens['Oi']}, "
      f"Oo {cens['Oo']} (N_min {NMIN_CELL} for Ii/Oi)")
    P(f"  S1 design: both-sides galaxies {len(QUAL_A)} "
      f"(N_min {NMIN_S1}); y-band edges (log10 y) "
      f"{np.array2string(EDGES_A, precision=2)}")
    ta_, qa_, ga_ = matcher_census(M_S1_A)
    P(f"  S1 matcher census (A2): inner points {ta_}, "
      f"window>={MINW} qualifying {qa_}, from {ga_} galaxies "
      f"(floors {NMIN_MPT} pts / {NMIN_S1} gal)")
    to_, qo_, go_ = matcher_census(M_OI_A)
    P(f"  Oi matcher census (A2-ii): Oi points {to_}, qualifying "
      f"{qo_}, from {go_} galaxies (descriptive co-read)")
    P(f"  WORLD-Y calibration: b = {B_Y:+.4f} per dex(y) "
      f"(composition term C = {C_Y:+.4f} over {N_CAL} gal)")
    P(f"  G10Y-2 {'PASS' if ok2 else 'FAIL'}")

    # -- G10Y-3 null machinery: size + branch-rate tables ------------
    P("")
    P("-- G10Y-3 size (noise-only, 2000 draws; inner perms 500 for")
    P("   the null tables -- the sky grammar uses the full 4000) --")
    rng = np.random.default_rng(20260923)
    fires = dict(ii_deep=0, oi_deep=0, s2=0, sig_ang=0, sig_dsk=0,
                 s1_surv=0, s1_coll=0, s1_pop=0)
    nsz = 2000
    for _ in range(nsz):
        rb = world_noise(GAL_A, rng)
        br, rec = grammar(GAL_A, rb, rng, M_S1_A, nperm=500)
        fires['ii_deep'] += rec['ii_deep']
        fires['oi_deep'] += rec['oi_deep']
        fires['s2'] += rec['s2_fire']
        fires['sig_ang'] += rec['sig_ang']
        fires['sig_dsk'] += rec['sig_dsk']
        fires['s1_surv'] += rec['survives']
        fires['s1_coll'] += rec['collapses']
        fires['s1_pop'] += rec['s1_pop']
    for k in fires:
        fires[k] /= nsz
    okI = 0.01 <= fires['ii_deep'] <= 0.10
    okO = 0.01 <= fires['oi_deep'] <= 0.10
    okS = 0.01 <= fires['s2'] <= 0.10
    P(f"  constituent P(d_Ii <= -2SE) = {fires['ii_deep']:.3f}  "
      f"{'PASS' if okI else 'FAIL'} [0.01, 0.10]")
    P(f"  constituent P(d_Oi <= -2SE) = {fires['oi_deep']:.3f}  "
      f"{'PASS' if okO else 'FAIL'} [0.01, 0.10]")
    P(f"  constituent P(S2 Holm+sign) = {fires['s2']:.3f}  "
      f"{'PASS' if okS else 'FAIL'} [0.01, 0.10]")
    P(f"  compound rates (printed, no bar; A1-iv): "
      f"sig_ang {fires['sig_ang']:.3f}, sig_dsk {fires['sig_dsk']:.3f}")
    P(f"  context: P(S1 populated) = {fires['s1_pop']:.3f}, "
      f"P(S1 survive|noise) = {fires['s1_surv']:.3f}, "
      f"P(S1 collapse|noise) = {fires['s1_coll']:.3f}")
    ok3a = okI and okO and okS

    P("")
    P("-- G10Y-3 branch rates (1000 draws/world; trap #32/#33) --")
    tables = {}
    for wnm, wfun in (('WORLD-Y', world_y), ('WORLD-R', world_r)):
        cnt = dict(B1=0, B2=0, B3=0, B4=0, B4_s1surv=0)
        for _ in range(N_WORLD):
            rb = wfun(GAL_A, rng)
            br, rec = grammar(GAL_A, rb, rng, M_S1_A, nperm=500)
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

    # -- G10Y-4 y-match wiring (banded; A1-v) -------------------------
    P("")
    P("-- G10Y-4 y-match wiring (500 draws/world; matched S1) --")
    czero = ckeep = 0
    for _ in range(N_WIRE):
        mm, sm, mp_, ng_ = s1_matched(world_y(GAL_A, rng), M_S1_A)
        if mp_ >= NMIN_MPT and ng_ >= NMIN_S1 and np.isfinite(mm) \
           and abs(mm) < sm:
            czero += 1
        mm, sm, mp_, ng_ = s1_matched(world_r(GAL_A, rng), M_S1_A)
        if mp_ >= NMIN_MPT and ng_ >= NMIN_S1 and np.isfinite(mm) \
           and mm <= -sm:
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
    P("STAGE 10Y SKY READ (post-prereg-commit; A1 form)")
    P("")
    rng = np.random.default_rng(20260923 + 1)

    for legnm, gal, sub, edges, qual, M in (
            ('ANCHORED', GAL_A, subP, EDGES_A, QUAL_A, M_S1_A),
            ('FLOW', GAL_F, subF, EDGES_F, QUAL_F, M_S1_F)):
        res_by_gal = [gd_['res'] for gd_ in gal]
        P(f"---- {legnm} leg (best form "
          f"{A2_A if legnm == 'ANCHORED' else A2_F}) ----")
        d, m, s = contrast_plain(gal, res_by_gal)
        P(f"  S0 plain contrast: {m:+.4f} +/- {s:.4f} over {len(d)} gal")
        mm, sm, mp_, ng_ = s1_matched(res_by_gal, M)
        P(f"  S1 y-matched (A2 primary): {mm:+.4f} +/- {sm:.4f} "
          f"(jackknife; {mp_} matched inner pts, {ng_} gal)")
        mb, sb, nb_, nq, cb, ni, no = s1_banded(gal, res_by_gal,
                                                edges, qual)
        P(f"  S1 co-read banded (A1-ii, demoted): {mb:+.4f} +/- "
          f"{sb:.4f} ({nb_} usable bands, {nq} both-sides gal)")
        oc = overlap_census(gal, res_by_gal)
        P(f"  S1 co-read per-galaxy overlap match: {len(oc)} gal "
          f"qualify" + (f", mean {np.mean(oc):+.4f}" if len(oc) >= 2
                        else " (unpopulated, as the A1 census found)"))
        for deg, tag in ((1, 'lin'), (2, 'quad'), (3, 'cub')):
            _, md, sd_ = detrend_contrast(gal, res_by_gal, deg)
            P(f"  S1 co-read detrend-{tag}: {md:+.4f} +/- {sd_:.4f} "
              f"(biased toward collapse; no branch weight)")
        s2 = s2_tests(gal, res_by_gal, rng)
        for tag in ('S2a_thmin', 'S2b_dist'):
            rho, p, sgn, n = s2[tag]
            P(f"  {tag}: rho = {rho:+.3f}, perm p = {p:.4f} "
              f"(smearing sign {'+' if sgn > 0 else '-'}; n = {n})")
        (dii, mii, sii), (doi, moi, soi), n_io = s3_cells(gal,
                                                          res_by_gal)
        P(f"  S3 d_Ii (inner-both vs Oo):      {mii:+.4f} +/- "
          f"{sii:.4f} over {len(dii)} gal")
        P(f"  S3 d_Oi (outer-angle, inner-disk): {moi:+.4f} +/- "
          f"{soi:.4f} over {len(doi)} gal")
        P(f"  S3 Io census: {n_io} galaxies carry any "
          f"inner-angle-outer-disk point")
        pii, poi = s3_shift_p(gal, res_by_gal, mii, moi, rng)
        P(f"  S3 circular-shift p: Ii {pii:.4f}, Oi {poi:.4f} "
          f"(one-sided toward depression; {N_SHIFT_SIZE} draws)")
        if legnm == 'ANCHORED':
            vo, so_, mo_, ngo = s1_matched(res_by_gal, M_OI_A)
            P(f"  S3 co-read y-matched d_Oi (A2-ii): {vo:+.4f} +/- "
              f"{so_:.4f} ({mo_} Oi pts, {ngo} gal; descriptive)")
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
    P("---- LETTER (prereg section 4 + A1 grammar; ANCHORED) ----")
    rngL = np.random.default_rng(20260923 + 2)
    br, rec = grammar(GAL_A, [gd_['res'] for gd_ in GAL_A], rngL,
                      M_S1_A, nperm=N_PERM, want_p=True)
    P(f"  S1: matched pts {rec['s1_mpts']}, gal {rec['s1_ngal']}, "
      f"mean {rec['s1_mean']:+.4f} +/- {rec['s1_se']:.4f}; "
      f"populated {rec['s1_pop']}; survives {rec['survives']}, "
      f"collapses {rec['collapses']}")
    P(f"  S3: d_Ii {rec['ii_mean']:+.4f} +/- {rec['ii_se']:.4f} "
      f"(n {rec['ii_n']}, shift-p {rec.get('ii_p', float('nan')):.4f}); "
      f"d_Oi {rec['oi_mean']:+.4f} +/- {rec['oi_se']:.4f} "
      f"(n {rec['oi_n']}, shift-p {rec.get('oi_p', float('nan')):.4f})")
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
