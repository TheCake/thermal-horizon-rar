"""
STAGE 5B (O4): the environmental control -- real per-galaxy g_N,ext vs the
thermal deep trend.

4T/4U found the RAR intrinsic scatter x-dependent (deep bins broader), read
thermally as oscillator shot noise (N ~ 20-60 modes). 5A showed a deep-
weighted EFE template with FREE per-galaxy amplitudes can absorb the deep
bins -- but free amplitudes prove nothing about environment. The decisive
version needs the real fields. Chae et al. 2021 (ApJ 921, 104), Table 3
publishes log10 e_N,env = log10(g_N,env/a0) for 109 SPARC galaxies in the
SDSS footprint, for a "max clustering" and a "no clustering" mass model
(the two columns are nearly a constant 0.9-dex offset apart, so only the
overall amplitude distinguishes them -- we fit a global scale beta times
the max-clustering pattern).
  PDF: https://astroweb.case.edu/ssm/papers/Chae_2021_ApJ_921_104.pdf
  (downloaded to data/chae2021.pdf; text-extracted to data/chae2021_text.txt;
  Table 3 parsed below).

EFE model: the exact 1D QUMOND/AQUAL collinear ratio, Chae+21 Eq. (2)
(simple-IF family):
    nu_e(y) = 1/2 + [ D - C ] / (2y),
    D = (y+e)*sqrt(1+4/(y+e)),  C = e*sqrt(1+4/e).
This is the formula their Table-2 e-tilde fits are calibrated to; it is the
MAXIMAL (collinear) geometry, so beta absorbs both the clustering amplitude
and the orientation average (beta ~ 0.5 expected even if max clustering is
exactly right; beta ~ 0 = no environmental signature). We apply its
suppression RATIO as a template on the BE mean (family-consistent with
4T/4U/5A): dlg(y, e) = lg nu_e^simple(y, e) - lg nu^simple(y).
NOTE: 5A used a tanh-form approximation (Chae-Milgrom style) whose
curvature-term sign was ambiguous in our transcription; 5A is unaffected
(its per-galaxy template coefficients were free-sign), but here the sign
is physical, so we use the exact Eq. (2) and gate it.

Contest on the matched subsample S (hierarchical per-galaxy disk-M/L
offsets, 0.1-dex prior, BE mean, SPARC only -- the 4U treatment):
  E0: const scatter                     (thermal off, env off)
  E1: oscillator+floor                  (thermal on,  env off)
  E2: const scatter + beta*e_g EFE      (thermal off, env on)
  E3: osc+floor     + beta*e_g EFE      (both)
plus the free 6-bin scatter profile (frozen 4T sextile edges) at beta=0
and beta=beta_hat, to see WHICH bins the real-pattern EFE absorbs.

Direct channels:
  beta profile (E2 grid) with D1 interval;
  scramble gate: 10 permutations of e_g across galaxies -- the EFE credit
    must be pattern-specific, not generic-deep-freedom;
  correlation: per-galaxy deep-bin mean residual (from E0) vs log e_g
    (Spearman), and observed-vs-predicted slope at beta=1.

Gates:
  G0 formula: e->0 identity to 1e-12; y->0 cap matches the analytic limit;
     suppression monotone in e.
  G1 full-sample regression: the machinery with no EFE reproduces 4U
     M0h/M1bh (-10371.19 / -10393.96) within 1.5.
  G2 injection on S: beta_true=3 (with drawn M/L offsets + noise)
     recovered within a factor ~1.5; beta_true=0 recovers ~0.
Writes data/stage5b_envtest.txt.
"""
import glob, math, os, re
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import spearmanr

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)
SIG_ML = 0.1*LN10

# ---------------- SPARC loader (verbatim 4U) ----------------
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

g_gas, g_dsk, g_bul, gobs, sig, gal_id, gal_name = [], [], [], [], [], [], {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q = meta.get(name, (0, 3))
    if inc < 30 or q > 2: continue
    kept += 1
    gal_name[gi] = name
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

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def osc_shape(x):
    xc = np.minimum(np.clip(x, 1e-9, None), 80.0)
    n = 1.0/np.expm1(xc)
    return np.sqrt(n/(1.0+n))

# frozen full-sample sextile edges (4T convention)
x_fid = np.sqrt((g_gas + g_dsk + g_bul)/A0_FID)
QE = np.quantile(x_fid, np.linspace(0, 1, 7))
QE[0], QE[-1] = 0.0, np.inf
BIN_FID = np.clip(np.searchsorted(QE, x_fid, side='right')-1, 0, 5)

# ---------------- Chae+21 Table 3 parse ----------------
rows = {}
pat = re.compile(r'^(.+?)\s+(-\d+\.\d+)\s*\+/-\s*(\d+\.\d+)\s+'
                 r'(-\d+\.\d+)\s*\+/-\s*(\d+\.\d+)\s*$')
with open('data/chae2021_text.txt', encoding='utf-8') as f:
    for l in f:
        l = l.replace('−', '-').replace('±', '+/-').strip()
        m = pat.match(l)
        if not m: continue
        nm = m.group(1).replace(' ', '')
        rows[nm] = (float(m.group(2)), float(m.group(3)),
                    float(m.group(4)), float(m.group(5)))
assert len(rows) == 109, f"Table 3 parse: got {len(rows)} rows, expected 109"

# match to our kept galaxies
gal_e, gal_eerr = {}, {}
for gi, nm in gal_name.items():
    if nm in rows:
        gal_e[gi] = 10.0**rows[nm][0]        # max clustering, linear e_N
        gal_eerr[gi] = rows[nm][1]
S_gals = sorted(gal_e)
inS = np.isin(gal_id, S_gals)
sub = np.where(inS)[0]
ugS = np.array(S_gals)
NGalS = len(ugS)
gmapS = {g: i for i, g in enumerate(ugS)}
gidxS = np.array([gmapS[g] for g in gal_id[sub]])
GIDXS_S = [np.where(gidxS == i)[0] for i in range(NGalS)]
e_gal = np.array([gal_e[g] for g in ugS])
gS_gas, gS_dsk, gS_bul = g_gas[sub], g_dsk[sub], g_bul[sub]
lgobsS, sigS2 = lgobs[sub], sig2[sub]
BIN_S = BIN_FID[sub]

L = [f"STAGE 5B environmental control: Chae+21 Table 3 parsed 109 rows; "
     f"matched {len(gal_e)}/{kept} kept galaxies -> subsample S: "
     f"{NGalS} galaxies, {len(sub)} points",
     f"e_N (max clustering): median {np.median(e_gal):.5f}, "
     f"16/84 pct {np.percentile(e_gal,16):.5f}/{np.percentile(e_gal,84):.5f}; "
     f"log-err median {np.median([gal_eerr[g] for g in ugS]):.2f} dex",
     ""]

# ---------------- EFE template: exact 1D QUMOND Eq.(2) ----------------
def nu_e_simple(y, e):
    """Chae+21 Eq.(2): exact 1D collinear QUMOND/AQUAL ratio, simple IF.
    Vectorized over y and e (e >= 0)."""
    y = np.clip(np.asarray(y, float), 1e-14, None)
    e = np.clip(np.asarray(e, float), 0.0, None)
    ye = y + e
    D = ye*np.sqrt(1.0 + 4.0/ye)
    C = np.where(e > 0, e*np.sqrt(1.0 + 4.0/np.clip(e, 1e-300, None)), 0.0)
    return 0.5 + (D - C)/(2.0*y)

def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))

def dlg_efe(y, e):
    """log10 suppression template rel. isolated simple-IF; 0 where e=0."""
    out = np.log10(nu_e_simple(y, e)) - np.log10(nu_simple(y))
    return np.where(np.asarray(e) > 0, out, 0.0)

# G0 formula gates
yt = np.logspace(-3, 1, 50)
g0a = np.max(np.abs(nu_e_simple(yt, 0.0) - nu_simple(yt)))
ee = 0.005
cap_num = nu_e_simple(np.array([1e-9]), ee)[0]
sq = math.sqrt(0.25 + 1.0/ee)
cap_ana = 0.5 + sq - 1.0/(2.0*ee*sq)
dsup = dlg_efe(np.full(3, 0.02), np.array([0.001, 0.005, 0.02]))
g0ok = (g0a < 1e-12 and abs(cap_num-cap_ana)/cap_ana < 1e-6
        and dsup[0] > dsup[1] > dsup[2] and dsup[0] < 0)
L.append(f"G0 formula: e->0 max|diff| = {g0a:.2e}; y->0 cap {cap_num:.4f} "
         f"(analytic {cap_ana:.4f}); suppression at y=0.02 for e=(.001,.005,.02): "
         f"{np.round(dsup,4).tolist()} -> {'PASS' if g0ok else 'FAIL'}")
L.append(f"  template size at (y,e)=(0.01, 0.0045): "
         f"{dlg_efe(np.array([0.01]), np.array([0.0045]))[0]:+.4f} dex")
L.append("")

# ---------------- hierarchical fit machinery ----------------
S_BOUNDS = {'M0h': [(1e-4, 0.5)],
            'M1bh': [(-2, 12), (0.0, 0.5)],
            'M2h': [(1e-4, 0.5)]*6}
S_START = {'M0h': [[0.10], [0.06]],
           'M1bh': [[math.log(30.0), 0.05], [math.log(8.0), 0.10]],
           'M2h': [[0.10]*6, [0.06]*6]}

def s_model_S(model, sth, x):
    if model == 'M0h':
        return np.full(len(x), sth[0])
    if model == 'M1bh':
        lnN, sf = sth
        return np.sqrt((osc_shape(x)/(math.sqrt(math.exp(lnN))*LN10))**2 + sf*sf)
    if model == 'M2h':
        return np.asarray(sth)[BIN_S]
    raise ValueError(model)

def m2S(la0, f, beta, dml, model, sth, prior_sig, e_vec):
    fac = f*np.exp(dml[gidxS])
    gN = gS_gas + fac*gS_dsk + gS_bul
    y = gN/10**la0
    x = np.sqrt(y)
    lgm = np.log10(gN*nu_be(y)) + dlg_efe(y, beta*e_vec[gidxS])
    s = s_model_S(model, sth, x)
    se2 = sigS2 + s*s
    r = lgobsS - lgm
    return (np.sum(r*r/se2 + np.log(se2))
            + np.sum(dml*dml)/(prior_sig*prior_sig))

def fit_S(model, use_beta=False, prior_sig=SIG_ML, rounds=3, e_vec=None,
          beta_fix=None, th0=None, lg_override=None):
    """Hierarchical fit on S. Globals: [la0, f] (+ [beta]) + sth."""
    global lgobsS
    lg_save = lgobsS
    if lg_override is not None: lgobsS = lg_override
    if e_vec is None: e_vec = e_gal
    dml = np.zeros(NGalS)
    best = None
    nb = 1 if use_beta else 0
    try:
        for rd in range(rounds):
            def gobj(th):
                la0, f = th[0], th[1]
                if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
                if use_beta:
                    beta = th[2]
                    if not (0.0 <= beta <= 30.0): return 1e12
                else:
                    beta = 0.0 if beta_fix is None else beta_fix
                sth = th[2+nb:]
                for v, (lo, hi) in zip(sth, S_BOUNDS[model]):
                    if not (lo <= v <= hi): return 1e12
                return m2S(la0, f, beta, dml, model, sth, prior_sig, e_vec)
            starts = ([list(best.x)] if best is not None else []) + \
                     ([list(th0)] if th0 is not None else [])
            for s0 in S_START[model]:
                base = [math.log10(A0_FID), 1.0]
                if use_beta:
                    starts.append(base + [1.0] + s0)
                    starts.append(base + [6.0] + s0)
                else:
                    starts.append(base + s0)
            gb = None
            for t0 in starts:
                b = minimize(gobj, t0, method='Nelder-Mead',
                             options=dict(maxiter=6000, xatol=1e-6, fatol=1e-6))
                if gb is None or b.fun < gb.fun: gb = b
            best = gb
            la0, f = best.x[0], best.x[1]
            beta = best.x[2] if use_beta else (0.0 if beta_fix is None
                                               else beta_fix)
            sth = best.x[2+nb:]
            if prior_sig < 1e-3: continue
            for gi in range(NGalS):
                mm = GIDXS_S[gi]
                ev = beta*e_vec[gi]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN = gS_gas[mm] + fc*gS_dsk[mm] + gS_bul[mm]
                    yv = gN/10**la0
                    lgm = (np.log10(gN*nu_be(yv))
                           + dlg_efe(yv, np.full(len(mm), ev)))
                    s = s_model_S(model, sth, np.sqrt(yv))[0:len(mm)] \
                        if model != 'M2h' else np.asarray(sth)[BIN_S[mm]]
                    if model == 'M0h': s = np.full(len(mm), sth[0])
                    elif model == 'M1bh':
                        lnN, sf = sth
                        s = np.sqrt((osc_shape(np.sqrt(yv))
                                     / (math.sqrt(math.exp(lnN))*LN10))**2
                                    + sf*sf)
                    se2 = sigS2[mm] + s*s
                    r = lgobsS[mm] - lgm
                    return (np.sum(r*r/se2 + np.log(se2))
                            + dl*dl/(prior_sig*prior_sig))
                dml[gi] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                          method='bounded').x
    finally:
        lgobsS = lg_save
    return best, dml

# ---------------- G1: full-sample regression vs 4U ----------------
ug = np.unique(gal_id)
NGal = len(ug)
gmap = {g: i for i, g in enumerate(ug)}
gidx = np.array([gmap[g] for g in gal_id])
GIDXS_F = [np.where(gidx == i)[0] for i in range(NGal)]

def fit_full(model, rounds=4):
    dml = np.zeros(NGal)
    best = None
    for rd in range(rounds):
        def gobj(th):
            la0, f = th[0], th[1]
            if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
            sth = th[2:]
            for v, (lo, hi) in zip(sth, S_BOUNDS[model]):
                if not (lo <= v <= hi): return 1e12
            fac = f*np.exp(dml[gidx])
            gN = g_gas + fac*g_dsk + g_bul
            x = np.sqrt(gN/10**la0)
            gm = gN*nu_be(gN/10**la0)
            if model == 'M0h': s = np.full(len(x), sth[0])
            else:
                lnN, sf = sth
                s = np.sqrt((osc_shape(x)/(math.sqrt(math.exp(lnN))*LN10))**2
                            + sf*sf)
            se2 = sig2 + s*s
            r = lgobs - np.log10(gm)
            return np.sum(r*r/se2 + np.log(se2)) + np.sum(dml*dml)/(SIG_ML**2)
        starts = ([list(best.x)] if best is not None else []) + \
                 [[math.log10(A0_FID), 1.0] + s0 for s0 in S_START[model]]
        gb = None
        for t0 in starts:
            b = minimize(gobj, t0, method='Nelder-Mead',
                         options=dict(maxiter=6000, xatol=1e-6, fatol=1e-6))
            if gb is None or b.fun < gb.fun: gb = b
        best = gb
        la0, f, sth = gb.x[0], gb.x[1], gb.x[2:]
        for gi in range(NGal):
            mm = GIDXS_F[gi]
            def od(dl):
                fc = f*math.exp(dl)
                gN = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                x = np.sqrt(gN/10**la0)
                gm = gN*nu_be(gN/10**la0)
                if model == 'M0h': s = np.full(len(mm), sth[0])
                else:
                    lnN, sf = sth
                    s = np.sqrt((osc_shape(x)/(math.sqrt(math.exp(lnN))*LN10))**2
                                + sf*sf)
                se2 = sig2[mm] + s*s
                r = lgobs[mm] - np.log10(gm)
                return np.sum(r*r/se2 + np.log(se2)) + dl*dl/(SIG_ML**2)
            dml[gi] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                      method='bounded').x
    return best

gf0 = fit_full('M0h')
gf1 = fit_full('M1bh')
ok1 = abs(gf0.fun-(-10371.19)) < 1.5 and abs(gf1.fun-(-10393.96)) < 1.5
L.append(f"G1 full-sample regression: M0h {gf0.fun:.2f} (4U -10371.19), "
         f"M1bh {gf1.fun:.2f} (4U -10393.96) -> {'PASS' if ok1 else 'FAIL'}")
L.append("")

# ---------------- E-contest on S ----------------
E0, dml0 = fit_S('M0h')
E1, _ = fit_S('M1bh')
E2, dml2 = fit_S('M0h', use_beta=True)
E3, _ = fit_S('M1bh', use_beta=True)
beta2 = E2.x[2]; beta3 = E3.x[2]
L.append("E-contest on S (hierarchical delta_d, BE mean, SPARC only):")
L.append(f"  E0 const:            {E0.fun:10.2f}  s0={E0.x[2]:.4f}")
L.append(f"  E1 osc+floor:        {E1.fun:10.2f}  N-hat={math.exp(E1.x[2]):8.1f}"
         f"  floor={E1.x[3]:.4f}")
L.append(f"  E2 const+EFE:        {E2.fun:10.2f}  beta={beta2:.2f}  "
         f"s0={E2.x[3]:.4f}")
L.append(f"  E3 osc+floor+EFE:    {E3.fun:10.2f}  beta={beta3:.2f}  "
         f"N-hat={math.exp(E3.x[3]):8.1f}  floor={E3.x[4]:.4f}")
L.append(f"  D: thermal E1-E0 = {E1.fun-E0.fun:+.2f} | env E2-E0 = "
         f"{E2.fun-E0.fun:+.2f} | both E3-E0 = {E3.fun-E0.fun:+.2f} | "
         f"E3-E1 (env credit after thermal) = {E3.fun-E1.fun:+.2f} | "
         f"E3-E2 (thermal credit after env) = {E3.fun-E2.fun:+.2f}")
L.append("")

# free-bin profile at beta=0 and beta_hat
B0, _ = fit_S('M2h')
Bb, _ = fit_S('M2h', use_beta=False, beta_fix=beta2)
L.append(f"free 6-bin scatter (frozen 4T edges), beta=0:      "
         f"{np.round(B0.x[2:8],4).tolist()}  ({B0.fun:.2f})")
L.append(f"free 6-bin scatter,                beta={beta2:.2f}: "
         f"{np.round(Bb.x[2:8],4).tolist()}  ({Bb.fun:.2f})")
L.append("  [bins deep->Newtonian; bin4 = the x~1 bump bin on full-sample "
         "edges]")
L.append("")

# ---------------- beta profile (E2 scatter model) ----------------
BGRID = [0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0]
prof = []
thw = None
for bt in BGRID:
    b, _ = fit_S('M0h', use_beta=False, beta_fix=bt, rounds=2, th0=thw)
    thw = b.x
    prof.append(b.fun)
prof = np.array(prof)
ib = int(np.argmin(prof))
dz = prof - prof.min()
lo = hi = None
for j in range(ib, -1, -1):
    if dz[j] > 1.0:
        lo = np.interp(1.0, [dz[j+1], dz[j]], [BGRID[j+1], BGRID[j]]); break
for j in range(ib, len(BGRID)):
    if dz[j] > 1.0:
        hi = np.interp(1.0, [dz[j-1], dz[j]], [BGRID[j-1], BGRID[j]]); break
L.append("beta profile (const-scatter treatment):")
L.append("  " + "  ".join(f"{bt:g}:{p:.1f}" for bt, p in zip(BGRID, prof)))
L.append(f"  min at beta = {BGRID[ib]:g}; D1 interval "
         f"[{'<0' if lo is None else f'{lo:.2f}'}, "
         f"{'>20' if hi is None else f'{hi:.2f}'}]; "
         f"D(-2lnL) beta=0 vs min = {dz[0]:+.2f}")
L.append("")

# ---------------- scramble gate ----------------
rng = np.random.default_rng(19)
cred_real = E0.fun - E2.fun
creds = []
for k in range(10):
    ep = e_gal[rng.permutation(NGalS)]
    bs, _ = fit_S('M0h', use_beta=True, rounds=2, e_vec=ep, th0=E2.x)
    creds.append(E0.fun - bs.fun)
creds = np.array(creds)
L.append(f"scramble gate (10 perms of e_g): real EFE credit = {cred_real:.2f}"
         f" vs scrambled {np.round(creds,2).tolist()}")
L.append(f"  scrambles beaten: {int((creds < cred_real).sum())}/10; "
         f"scramble mean {creds.mean():.2f} +/- {creds.std(ddof=1):.2f}")
L.append("")

# ---------------- correlation channel ----------------
la0_0, f_0 = E0.x[0], E0.x[1]
fac = f_0*np.exp(dml0[gidxS])
gN = gS_gas + fac*gS_dsk + gS_bul
y0 = gN/10**la0_0
r0 = lgobsS - np.log10(gN*nu_be(y0))
pred1 = dlg_efe(y0, e_gal[gidxS])           # beta=1 prediction
deep = BIN_S <= 1
rbar, pbar, elog = [], [], []
for gi in range(NGalS):
    mm = GIDXS_S[gi]
    md = mm[deep[mm]]
    if len(md) < 2: continue
    rbar.append(np.mean(r0[md])); pbar.append(np.mean(pred1[md]))
    elog.append(math.log10(e_gal[gi]))
rbar, pbar, elog = map(np.array, (rbar, pbar, elog))
rho, pval = spearmanr(rbar, elog)
sl = (np.sum((pbar-pbar.mean())*(rbar-rbar.mean()))
      / np.sum((pbar-pbar.mean())**2)) if len(rbar) > 3 else float('nan')
L.append(f"correlation channel ({len(rbar)} galaxies with >=2 deep-bin "
         f"points):")
L.append(f"  Spearman(deep residual, log e_N) rho = {rho:+.3f} (p = {pval:.3f})"
         f"   [EFE predicts NEGATIVE]")
L.append(f"  slope of observed on predicted(beta=1) = {sl:+.2f}   "
         f"[~beta if Eq.(2) pattern right]")
L.append("")

# ---------------- G2: injections on S ----------------
d_t = rng.normal(0, SIG_ML, NGalS)
for btru in (3.0, 0.0):
    fac_t = 1.0*np.exp(d_t[gidxS])
    gN_t = gS_gas + fac_t*gS_dsk + gS_bul
    y_t = gN_t/A0_FID
    lg_t = (np.log10(gN_t*nu_be(y_t)) + dlg_efe(y_t, btru*e_gal[gidxS])
            + rng.normal(0, np.sqrt(sigS2 + 0.07**2)))
    bi, _ = fit_S('M0h', use_beta=True, rounds=2, lg_override=lg_t)
    L.append(f"G2 injection beta_true={btru:g}: beta_rec = {bi.x[2]:.2f}"
             + ("  -> " + ("PASS" if 2.0 <= bi.x[2] <= 4.5 else "FAIL")
                if btru == 3.0 else
                "  -> " + ("PASS" if bi.x[2] <= 0.8 else "FAIL")))
L.append("")

# ---------------- verdict ----------------
L.append("VERDICT LOGIC: env real => E2 credit ~ E1 credit, pattern-specific "
         "(beats scrambles), rho<0, slope~beta; thermal survives control => "
         "E2 credit small/unspecific, E3 keeps the osc term.")

out = "\n".join(L)
print(out)
with open('data/stage5b_envtest.txt', 'w') as f:
    f.write(out+"\n")
print("\nsaved: data/stage5b_envtest.txt")
