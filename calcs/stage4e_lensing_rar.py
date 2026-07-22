"""
STAGE 4E: the rung-2 (NNLO) test with the weak-lensing RAR anchor (TODO #16).
Within the c1=1/2 branch the next Bernoulli rung separates the laws:
    BE:     nu = 1/(1-exp(-x)), x=sqrt(y):   c2 = 1/12 = 0.0833
    simple: nu = 1/2+sqrt(1/4+1/y):          c2 = 1/8  = 0.1250
    standard (branch control, killed in 4B): c1 = 0,  c2 = 1/4
4B could not separate 1/12 from 1/8 on SPARC alone: (a0, f_ML) refitting
absorbs the shape difference, and its raw-chi2 likelihood (chi2/dof ~ 57, no
intrinsic scatter) makes deltas uninterpretable. Two upgrades here:
  (1) proper likelihood: -2lnL = sum[ r^2/(sig^2+s_int^2) + ln(sig^2+s_int^2) ]
      with the intrinsic RAR scatter s_int profiled per family;
  (2) the weak-lensing RAR extends the relation ~2 dex below SPARC, anchoring
      the deep-MOND normalization sqrt(g_bar*a0) so the SPARC mid-range and
      Wien tail must carry the shape alone.

Data:
  A) Mistele, McGaugh, Lelli, Schombert & Li 2024 (JCAP 04, 020;
     arXiv:2310.15248) Table 1: 15 stacked exactly-deprojected points,
     isolated KiDS-bright lenses, log g_bar in [-14.86, -11.41]; per-point
     stat (+) syst in quadrature; GLOBAL 0.2-dex stellar-mass systematic
     modeled as one log-g_bar offset nuisance `dlt` with N(0, 0.2) prior.
     File: data/lensing_rar/mistele2024_table1.txt (from arXiv source
     plots/RAR-table.tex, fetched 2026-07-22).
  B) Brouwer et al. 2021 (A&A 650, A113) KiDS-1000 release,
     Fig-4-5-C1_RAR-KiDS-isolated (SIS-approx conversion g_obs = 4G*ESD,
     their Eq. 7) + full covariance: independent-pipeline GLS cross-check.
     https://kids.strw.leidenuniv.nl/sci_data/brouwer2021_rar.tar
  C) SPARC rotmod, cuts identical to 4A/4B (inc>30, Q<=2, eV/V<0.10).

Gates: G1 a0 sane + consistent with 4B; G2 regression vs stored 4B chi2;
G3 Newton control on lensing; G4 branch control at depth; G5 verdict-sign
stability across windows / EFE / mass-offset treatments.
Writes data/stage4e_lensing.txt.
"""
import glob, math, os, re
import numpy as np
from scipy.optimize import minimize

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2          # dex, global stellar-mass systematic (lensing g_bar)
LENS_CUT_FID = -14.25      # fiducial: drop the 2 deepest points (syst 0.25/0.67)

# ---------------- SPARC (identical to 4A/4B) ----------------
meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        name, inc, q = t[0], float(t[5]), int(t[17])
        meta[name] = (inc, q)
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

# ---------------- lensing A: Mistele+24 ----------------
ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2

# ---------------- lensing B: Brouwer+21 (cross-check) ----------------
G_PC, PC_M = 4.52e-30, 3.086e16
B = np.loadtxt('data/lensing_rar/Fig-4-5-C1_RAR-KiDS-isolated_Nobins.txt')
b_gbar = B[:, 0]
b_gobs = 4*G_PC*(B[:, 1]/B[:, 4])*PC_M
CV = np.loadtxt('data/lensing_rar/Fig-4-5-C1_RAR-KiDS-isolated_covmatrix.txt')
nb = len(b_gbar)
b_cov = np.zeros((nb, nb))
for row in CV:
    i = int(np.argmin(np.abs(np.log10(row[2]) - np.log10(b_gbar))))
    j = int(np.argmin(np.abs(np.log10(row[3]) - np.log10(b_gbar))))
    b_cov[i, j] = row[4]/row[6]          # bias-corrected ESD covariance
b_cov *= (4*G_PC*PC_M)**2                # -> covariance of g_obs (linear)
b_icov = np.linalg.inv(b_cov)

# ---------------- nu families ----------------
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_newton(y):
    return np.ones_like(np.asarray(y, dtype=float))
FAMS = {'BE': nu_be, 'simple': nu_simple, 'standard': nu_standard}

LN10 = math.log(10)

# ---------------- joint likelihood ----------------
def m2ll(th, nu, sel, w_pt, lmask, lens_obs, use_sparc=True, use_lens=True,
         efe=0.0, dlt_fixed=None):
    la0, f, s_int, dlt = th
    if dlt_fixed is not None: dlt = dlt_fixed
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    out = 0.0
    if use_sparc:
        gN = (g_gas + f*g_dsk + g_bul)[sel]
        gm = gN*nu(gN/a0)
        se2 = sig2[sel] + s_int*s_int
        r = lgobs[sel] - np.log10(gm)
        out += np.sum(w_pt[gal_id][sel]*(r*r/se2 + np.log(se2)))
    if use_lens:
        lg = l_gbar[lmask] + dlt
        yl = 10**lg/a0
        if efe > 0: yl = np.sqrt(yl*yl + efe*efe)
        lgm = lg + np.log10(nu(yl))
        rl = lens_obs[lmask] - lgm
        out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
        if dlt_fixed is None:
            out += (dlt/DELTA_PRIOR)**2
    return out

def fit(nu, sel, w_pt, lmask, lens_obs, **kw):
    best = None
    for th0 in ([math.log10(A0_FID), 1.0, 0.08, 0.0],
                [math.log10(A0_FID)+0.1, 0.8, 0.12, -0.1]):
        b = minimize(lambda t: m2ll(t, nu, sel, w_pt, lmask, lens_obs, **kw),
                     th0, method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
        if best is None or b.fun < best.fun: best = b
    return best

L = [f"STAGE 4E rung-2 test with lensing anchor: {kept} galaxies, "
     f"{len(gobs)} SPARC points + {len(l_gbar)} Mistele+24 lensing points "
     f"(fiducial keeps log g_bar >= {LENS_CUT_FID}: "
     f"{int((l_gbar >= LENS_CUT_FID).sum())})",
     "-2lnL = SPARC[r^2/(sig^2+s_int^2)+ln(.)] + lens[r^2/sig_L^2+ln(.)] "
     "+ (dlt/0.2)^2; free (la0, f_ML, s_int, dlt) per family",
     "predictions within 1/2-branch: BE c2=1/12=0.0833 | simple c2=1/8=0.1250"]

allsel = np.ones(len(gobs), bool)
ones = np.ones(gal_id.max()+1)
y_fid = (g_gas + g_dsk + g_bul)/A0_FID

# ---------------- G2: regression vs stored 4B (raw chi2, no lens) ----------------
def chi2_4b(nu, sel):
    def c2(th):
        la0, f = th
        gN = (g_gas + f*g_dsk + g_bul)[sel]
        r = (lgobs[sel] - np.log10(gN*nu(gN/10**la0)))/sig[sel]
        return np.sum(r*r)
    return minimize(c2, [math.log10(A0_FID), 1.0], method='Nelder-Mead',
                    options=dict(maxiter=2000, xatol=1e-7, fatol=1e-7)).fun
stored = {}
if os.path.exists('data/stage4b_branch.txt'):
    for m in re.finditer(r'\[all y\]\s+(\w+): chi2=\s*([\d.]+)',
                         open('data/stage4b_branch.txt').read()):
        stored[m.group(1)] = float(m.group(2))
g2 = []
for name, nu in FAMS.items():
    c = chi2_4b(nu, y_fid < 30.0)   # 4B's "all y" was y_fid<30 (n=2693)
    ok = name in stored and abs(c - stored[name]) < 1.0
    g2.append(f"  G2 {name}: recomputed 4B chi2 = {c:.1f} vs stored "
              f"{stored.get(name, float('nan')):.1f} -> {'OK' if ok else 'FAIL'}")
L += ["", "G2 regression (4B objective, y<30 as in 4B 'all y'):"] + g2

# ---------------- lensing-only sanity + Newton control ----------------
lmask_fid = l_gbar >= LENS_CUT_FID
L.append("")
L.append("Lensing-only fits (la0, dlt free; f_ML/s_int irrelevant):")
for name, nu in list(FAMS.items()) + [('Newton', nu_newton)]:
    b = fit(nu, allsel, ones, lmask_fid, l_gobs, use_sparc=False)
    L.append(f"  {name:>8}: -2lnL = {b.fun:9.2f}  a0 = {10**b.x[0]:.3e}  "
             f"dlt = {b.x[3]:+.3f}")
L.append("  (G3: Newton must be catastrophic; a0-dlt near-degenerate at LO "
         "-> lensing-only a0 is prior-limited by design)")

# ---------------- fiducial joint fit ----------------
L.append("")
L.append("JOINT fiducial (all SPARC + lensing >= -14.25, dlt ~ N(0,0.2)):")
fid = {}
for name, nu in FAMS.items():
    b = fit(nu, allsel, ones, lmask_fid, l_gobs)
    fid[name] = b
    la0, f, s_int, dlt = b.x
    # decompose
    m_sp = m2ll(b.x, nu, allsel, ones, lmask_fid, l_gobs, use_lens=False)
    m_le = b.fun - m_sp
    L.append(f"  {name:>8}: -2lnL = {b.fun:10.2f} (SPARC {m_sp:10.2f} | lens "
             f"{m_le:7.2f})  a0={10**la0:.3e}  f_ML={f:.2f}  "
             f"s_int={s_int:.3f}  dlt={dlt:+.3f}")
d_bs = fid['BE'].fun - fid['simple'].fun
d_bst = fid['BE'].fun - fid['standard'].fun
L.append(f"  Delta(-2lnL) BE - simple   = {d_bs:+.2f}  (joint shape statistic;"
         f" NOT a pure rung-2 readout -- see stage4e_diag for the regime"
         f" decomposition and galaxy jackknife)")
L.append(f"  Delta(-2lnL) BE - standard = {d_bst:+.2f}  (G4 branch control)")

# ---------------- variants: windows / EFE / mass offset ----------------
L.append("")
L.append("G5 verdict stability [Delta(-2lnL) BE - simple]:")
variants = []
for tag, sel in (('SPARC y<0.5', y_fid < 0.5), ('SPARC y<1', y_fid < 1.0),
                 ('SPARC all', allsel)):
    for ltag, lm in (('lens all15', np.ones(15, bool)),
                     ('lens>=-14.25', lmask_fid),
                     ('lens>=-13.25', l_gbar >= -13.25)):
        r = {n: fit(nu, sel, ones, lm, l_gobs).fun for n, nu in FAMS.items()}
        variants.append((f"{tag:12} {ltag:13}", r['BE']-r['simple'],
                         r['BE']-r['standard']))
r0 = {n: fit(nu, allsel, ones, lmask_fid, l_gobs, dlt_fixed=0.0).fun
      for n, nu in FAMS.items()}
variants.append(("SPARC all    dlt=0 (face)", r0['BE']-r0['simple'],
                 r0['BE']-r0['standard']))
for e in (0.01, 0.02, 0.03, 0.05):
    r = {n: fit(nu, allsel, ones, np.ones(15, bool), l_gobs, efe=e).fun
         for n, nu in FAMS.items()}
    variants.append((f"SPARC all    EFE e={e:.2f} L15", r['BE']-r['simple'],
                     r['BE']-r['standard']))
for tag, dbs, dbst in variants:
    L.append(f"  {tag:28}: BE-simple {dbs:+9.2f} | BE-standard {dbst:+10.2f}")

# ---------------- Brouwer GLS cross-check ----------------
L.append("")
L.append("Brouwer+21 release cross-check (GLS, full covariance, face-value "
         "masses; their SIS-approx conversion):")
for name, nu in list(FAMS.items()) + [('Newton', nu_newton)]:
    def gls(th, nu=nu):
        la0 = th[0]
        if not (-10.6 < la0 < -9.4): return 1e12
        r = b_gobs - b_gbar*nu(b_gbar/10**la0)
        return float(r @ b_icov @ r)
    b = minimize(gls, [math.log10(A0_FID)], method='Nelder-Mead',
                 options=dict(maxiter=1000, xatol=1e-7))
    L.append(f"  {name:>8}: chi2 = {b.fun:9.2f} (n=15)  a0 = {10**b.x[0]:.3e}")

# ---------------- interpretable c2 estimator ----------------
Y_DEEP = 0.05
F_GRID = np.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])

def a0_from_lensing(lens_obs, dlt):
    lg = l_gbar + dlt
    m = lg >= LENS_CUT_FID
    w = 1.0/l_sig2[m]
    la0 = math.log10(A0_FID)
    for _ in range(3):
        x = np.sqrt(10**lg[m]/10**la0)
        corr = np.log10(1 + 0.5*x + 0.1*x*x)
        la0 = np.sum(w*(2*(lens_obs[m] - corr) - lg[m]))/np.sum(w)
    return 10**la0

def c2_step(gal_w, fml, a0, ymax, fix_c1):
    gN = g_gas + fml*g_dsk + g_bul
    y = gN/a0
    w_pt = gal_w[gal_id]
    win = (y >= Y_DEEP) & (y < ymax) & (w_pt > 0)
    if win.sum() < 50: return None
    x = np.sqrt(y[win])
    u = gobs[win]/np.sqrt(gN[win]*a0) - 1.0
    su = LN10*np.sqrt(sig2[win] + 0.06**2)*np.maximum(1+u, 0.1)
    ww = w_pt[win]/su**2
    if fix_c1:
        t = u - 0.5*x
        c2 = np.sum(ww*x*x*t)/np.sum(ww*x**4)
        chi2 = np.sum(ww*(t - c2*x*x)**2)
        return (0.5, c2, chi2)
    X = np.stack([x, x*x], axis=1)
    XtW = X.T*ww
    beta = np.linalg.solve(XtW@X, XtW@u)
    chi2 = np.sum(ww*(u - X@beta)**2)
    return (beta[0], beta[1], chi2)

L.append("")
L.append(f"c2 estimator: a0 anchored on lensing (c1=1/2; window >= "
         f"{LENS_CUT_FID}), then SPARC WLS on {Y_DEEP} <= y < ymax "
         f"(s_floor 0.06 dex; f_ML profiled):")
NG = gal_id.max()+1
allg = np.unique(gal_id)
rng = np.random.default_rng(17)
for ymax in (0.5, 1.0):
    for fix_c1 in (True, False):
        a0L = a0_from_lensing(l_gobs, 0.0)
        cand = [(c2_step(ones, f, a0L, ymax, fix_c1), f) for f in F_GRID]
        cand = [(r, f) for r, f in cand if r is not None]
        (c1, c2, _), fbest = min(cand, key=lambda t: t[0][2])
        c2b = []
        for k in range(400):
            pick = rng.choice(allg, len(allg), replace=True)
            wg = np.zeros(NG)
            for g_ in pick: wg[g_] += 1
            lo = l_gobs + rng.normal(0, np.sqrt(l_sig2))
            a0k = a0_from_lensing(lo, 0.0)
            cnd = [(c2_step(wg, f, a0k, ymax, fix_c1), f) for f in F_GRID]
            cnd = [(r, f) for r, f in cnd if r is not None]
            if cnd: c2b.append(min(cnd, key=lambda t: t[0][2])[0][1])
        c2b = np.array(c2b)
        # mass systematic: dlt = +/-0.2 dex shifts the anchor
        c2_hi = min([(c2_step(ones, f, a0_from_lensing(l_gobs, +0.2), ymax,
                              fix_c1), f) for f in F_GRID],
                    key=lambda t: t[0][2])[0][1]
        c2_lo = min([(c2_step(ones, f, a0_from_lensing(l_gobs, -0.2), ymax,
                              fix_c1), f) for f in F_GRID],
                    key=lambda t: t[0][2])[0][1]
        sy = max(abs(c2_hi-c2), abs(c2_lo-c2))
        st = c2b.std(ddof=1)
        tag = "c1=1/2 fixed" if fix_c1 else f"c1 free ({c1:+.3f})"
        sep = abs(1/8 - 1/12)/math.sqrt(st*st + sy*sy)
        L.append(f"  ymax={ymax} {tag:16}: c2 = {c2:+.4f} +/- {st:.4f} (stat) "
                 f"+/- {sy:.4f} (mass syst)  a0_lens={a0L:.3e} f_ML={fbest:.1f}")
        L.append(f"      vs BE 1/12: {(c2-1/12)/math.sqrt(st*st+sy*sy):+.2f} "
                 f"sig | vs simple 1/8: {(c2-1/8)/math.sqrt(st*st+sy*sy):+.2f} "
                 f"sig | rung-2 resolving power: {sep:.2f} sig")

# ---------------- bootstrap sign stability (fiducial joint) ----------------
L.append("")
NBOOT = 200
wins_s, wins_st, d_list, a0_list = 0, 0, [], []
rng2 = np.random.default_rng(23)
for k in range(NBOOT):
    pick = rng2.choice(allg, len(allg), replace=True)
    wg = np.zeros(NG)
    for g_ in pick: wg[g_] += 1
    lo = l_gobs + rng2.normal(0, np.sqrt(l_sig2))
    r = {}
    for name, nu in FAMS.items():
        r[name] = fit(nu, allsel, wg, lmask_fid, lo)
    d = r['BE'].fun - r['simple'].fun
    d_list.append(d)
    a0_list.append(10**r['BE'].x[0])
    wins_s += int(d < 0)
    wins_st += int(r['BE'].fun < r['standard'].fun)
d_arr = np.array(d_list)
L.append(f"Bootstrap ({NBOOT} reps: SPARC galaxies x lensing perturbation):")
L.append(f"  Delta(-2lnL) BE - simple = {d_arr.mean():+.2f} +/- "
         f"{d_arr.std(ddof=1):.2f}; BE better in {wins_s}/{NBOOT}")
L.append(f"  BE beats standard in {wins_st}/{NBOOT} (G4 at depth)")
L.append(f"  joint a0 (BE) = {np.mean(a0_list):.3e} +/- {np.std(a0_list):.2e} "
         f"(G1: 4B all-y gave 1.206e-10)")

out = "\n".join(L)
print(out)
with open('data/stage4e_lensing.txt', 'w') as f_:
    f_.write(out+"\n")
print("\nsaved: data/stage4e_lensing.txt")
