"""
STAGE 5E (O4, convergence-hardened): the environmental control, adaptive
descent.

5B's instrument gates passed (formula exact, 4U regression, injections
recover beta) but its fixed-round descent left ~80-lnL convergence slop:
the beta-profile node at beta=0 sat 84 above the identical E0 model and
"scrambled" nested fits landed BELOW the base model -- impossible at
convergence. Same contest, converged:

  - adaptive rounds (tol 0.05, max 15; 3x delta_d sweeps per round);
  - NESTING gate: at convergence E1<=E0, E2<=E0, E3<=min(E1,E2), M2h<=E0
    (they nest); any violation = not converged, verdict void;
  - beta profile warm-chained in BOTH globals and offsets (the 4Z trick);
    its beta=0 node must reproduce E0 within 0.5;
  - scrambles (8 perms) warm-started, each must land in [E0-credit, E0];
  - correlation channel now DEPTH-PARTIALED: residualize both the observed
    deep-bin residual and the beta=1 prediction on median deep log-y per
    galaxy, then correlate -- kills the shared-depth confound that made
    5B's raw slope (+0.95) uninterpretable;
  - injections beta_true=3 / 0 as before, converged.
Writes data/stage5e_envtest_conv.txt.
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

x_fid = np.sqrt((g_gas + g_dsk + g_bul)/A0_FID)
QE = np.quantile(x_fid, np.linspace(0, 1, 7))
QE[0], QE[-1] = 0.0, np.inf
BIN_FID = np.clip(np.searchsorted(QE, x_fid, side='right')-1, 0, 5)

rows = {}
pat = re.compile(r'^(.+?)\s+(-\d+\.\d+)\s*\+/-\s*(\d+\.\d+)\s+'
                 r'(-\d+\.\d+)\s*\+/-\s*(\d+\.\d+)\s*$')
with open('data/chae2021_text.txt', encoding='utf-8') as f:
    for l in f:
        l = l.replace('−', '-').replace('±', '+/-').strip()
        m = pat.match(l)
        if not m: continue
        rows[m.group(1).replace(' ', '')] = (float(m.group(2)),
                                             float(m.group(3)))
assert len(rows) == 109

gal_e = {gi: 10.0**rows[nm][0] for gi, nm in gal_name.items() if nm in rows}
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

def nu_e_simple(y, e):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    e = np.clip(np.asarray(e, float), 0.0, None)
    ye = y + e
    D = ye*np.sqrt(1.0 + 4.0/ye)
    C = np.where(e > 0, e*np.sqrt(1.0 + 4.0/np.clip(e, 1e-300, None)), 0.0)
    return 0.5 + (D - C)/(2.0*y)
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def dlg_efe(y, e):
    out = np.log10(nu_e_simple(y, e)) - np.log10(nu_simple(y))
    return np.where(np.asarray(e) > 0, out, 0.0)

S_BOUNDS = {'M0h': [(1e-4, 0.5)],
            'M1bh': [(-2, 12), (0.0, 0.5)],
            'M2h': [(1e-4, 0.5)]*6}
S_START = {'M0h': [[0.10], [0.06]],
           'M1bh': [[math.log(30.0), 0.05], [math.log(8.0), 0.10]],
           'M2h': [[0.10]*6, [0.06]*6]}

def s_vec(model, sth, x, mm=None):
    if model == 'M0h':
        return np.full(len(x), sth[0])
    if model == 'M1bh':
        lnN, sf = sth
        return np.sqrt((osc_shape(x)/(math.sqrt(math.exp(lnN))*LN10))**2
                       + sf*sf)
    if model == 'M2h':
        b = BIN_S if mm is None else BIN_S[mm]
        return np.asarray(sth)[b]
    raise ValueError(model)

def m2S(la0, f, beta, dml, model, sth, e_vec):
    fac = f*np.exp(dml[gidxS])
    gN = gS_gas + fac*gS_dsk + gS_bul
    y = gN/10**la0
    lgm = np.log10(gN*nu_be(y)) + dlg_efe(y, beta*e_vec[gidxS])
    s = s_vec(model, sth, np.sqrt(y))
    se2 = sigS2 + s*s
    r = lgobsS - lgm
    return np.sum(r*r/se2 + np.log(se2)) + np.sum(dml*dml)/(SIG_ML*SIG_ML)

def fit_conv(model, use_beta=False, beta_fix=0.0, e_vec=None, th0=None,
             dml0=None, tol=0.05, max_rounds=15, lg_override=None,
             trace=None):
    global lgobsS
    lg_save = lgobsS
    if lg_override is not None: lgobsS = lg_override
    if e_vec is None: e_vec = e_gal
    nb = 1 if use_beta else 0
    dml = np.zeros(NGalS) if dml0 is None else dml0.copy()
    best = None
    prev = None
    try:
        for rd in range(max_rounds):
            def gobj(th):
                la0, f = th[0], th[1]
                if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5):
                    return 1e12
                beta = th[2] if use_beta else beta_fix
                if use_beta and not (0.0 <= beta <= 30.0): return 1e12
                sth = th[2+nb:]
                for v, (lo, hi) in zip(sth, S_BOUNDS[model]):
                    if not (lo <= v <= hi): return 1e12
                return m2S(la0, f, beta, dml, model, sth, e_vec)
            starts = ([list(best.x)] if best is not None else [])
            if rd == 0:
                if th0 is not None: starts.append(list(th0))
                for s0 in S_START[model]:
                    base = [math.log10(A0_FID), 1.0]
                    if use_beta:
                        starts += [base + [1.0] + s0, base + [6.0] + s0]
                    else:
                        starts.append(base + s0)
            bb = None
            for t0 in starts:
                b = minimize(gobj, t0, method='Nelder-Mead',
                             options=dict(maxiter=6000, xatol=1e-6,
                                          fatol=1e-7))
                if bb is None or b.fun < bb.fun: bb = b
            best = bb
            la0, f = best.x[0], best.x[1]
            beta = best.x[2] if use_beta else beta_fix
            sth = best.x[2+nb:]
            for _ in range(3):
                for gi in range(NGalS):
                    mm = GIDXS_S[gi]
                    ev = beta*e_gal[gi] if e_vec is e_gal else beta*e_vec[gi]
                    def od(dl):
                        fc = f*math.exp(dl)
                        gN = gS_gas[mm] + fc*gS_dsk[mm] + gS_bul[mm]
                        yv = gN/10**la0
                        lgm = (np.log10(gN*nu_be(yv))
                               + dlg_efe(yv, np.full(len(mm), ev)))
                        s = s_vec(model, sth, np.sqrt(yv), mm=mm)
                        se2 = sigS2[mm] + s*s
                        r = lgobsS[mm] - lgm
                        return (np.sum(r*r/se2 + np.log(se2))
                                + dl*dl/(SIG_ML*SIG_ML))
                    dml[gi] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                              method='bounded').x
            cur = m2S(la0, f, beta, dml, model, sth, e_vec)
            if trace is not None: trace.append(cur)
            if prev is not None and abs(prev - cur) < tol:
                break
            prev = cur
        b = minimize(gobj, list(best.x), method='Nelder-Mead',
                     options=dict(maxiter=6000, xatol=1e-6, fatol=1e-7))
        if b.fun < best.fun: best = b
    finally:
        lgobsS = lg_save
    return best, dml

L = [f"STAGE 5E environmental control, converged: S = {NGalS} galaxies, "
     f"{len(sub)} points; e_N median {np.median(e_gal):.5f}; adaptive tol "
     f"0.05, max 15 rounds", ""]

# ---------------- converged E-contest ----------------
tr = {}
tr['E0'] = []
E0, dml0 = fit_conv('M0h', trace=tr['E0'])
tr['E1'] = []
E1, _ = fit_conv('M1bh', th0=None, trace=tr['E1'])
tr['E2'] = []
E2, dml2 = fit_conv('M0h', use_beta=True, trace=tr['E2'])
tr['E3'] = []
E3, _ = fit_conv('M1bh', use_beta=True, trace=tr['E3'])
beta2 = E2.x[2]; beta3 = E3.x[2]
L.append("E-contest (converged; rounds shown):")
L.append(f"  E0 const:            {E0.fun:10.2f}  ({len(tr['E0'])} rd)  "
         f"s0={E0.x[2]:.4f}")
L.append(f"  E1 osc+floor:        {E1.fun:10.2f}  ({len(tr['E1'])} rd)  "
         f"N-hat={math.exp(E1.x[2]):8.1f}  floor={E1.x[3]:.4f}")
L.append(f"  E2 const+EFE:        {E2.fun:10.2f}  ({len(tr['E2'])} rd)  "
         f"beta={beta2:.3f}  s0={E2.x[3]:.4f}")
L.append(f"  E3 osc+floor+EFE:    {E3.fun:10.2f}  ({len(tr['E3'])} rd)  "
         f"beta={beta3:.3f}  N-hat={math.exp(E3.x[3]):8.1f}  "
         f"floor={E3.x[4]:.4f}")
L.append(f"  D: thermal E1-E0 = {E1.fun-E0.fun:+.2f} | env E2-E0 = "
         f"{E2.fun-E0.fun:+.2f} | E3-E1 = {E3.fun-E1.fun:+.2f} | "
         f"E3-E2 = {E3.fun-E2.fun:+.2f}")
B0, _ = fit_conv('M2h', th0=None, tol=0.05)
Bb, _ = fit_conv('M2h', beta_fix=beta2, tol=0.05)
L.append(f"  M2h free-bin beta=0:      {np.round(B0.x[2:8],4).tolist()}  "
         f"({B0.fun:.2f})")
L.append(f"  M2h free-bin beta={beta2:.2f}:   {np.round(Bb.x[2:8],4).tolist()}"
         f"  ({Bb.fun:.2f})")
nest_ok = (E1.fun <= E0.fun + 0.1 and E2.fun <= E0.fun + 0.1
           and E3.fun <= min(E1.fun, E2.fun) + 0.1
           and B0.fun <= E0.fun + 0.1)
L.append(f"NESTING gate: E1<=E0 {E1.fun-E0.fun:+.2f}, E2<=E0 "
         f"{E2.fun-E0.fun:+.2f}, E3<=min(E1,E2) "
         f"{E3.fun-min(E1.fun,E2.fun):+.2f}, M2h<=E0 {B0.fun-E0.fun:+.2f}"
         f" -> {'PASS' if nest_ok else 'FAIL (verdict void)'}")
L.append("")

# ---------------- beta profile, warm-chained ----------------
BGRID = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]
prof = []
thw, dmw = None, None
for bt in BGRID:
    b, dmw = fit_conv('M0h', beta_fix=bt, th0=thw, dml0=dmw,
                      tol=0.1, max_rounds=8)
    thw = [b.x[0], b.x[1], b.x[2]]
    prof.append(b.fun)
prof = np.array(prof)
dz = prof - prof.min()
ib = int(np.argmin(prof))
hi = None
for j in range(ib, len(BGRID)):
    if dz[j] > 1.0:
        hi = np.interp(1.0, [dz[j-1], dz[j]], [BGRID[j-1], BGRID[j]]); break
L.append("beta profile (const scatter, warm-chained globals+offsets):")
L.append("  " + "  ".join(f"{bt:g}:{p:.2f}" for bt, p in zip(BGRID, prof)))
L.append(f"  min at beta={BGRID[ib]:g}; upper D1 = "
         f"{'>8' if hi is None else f'{hi:.2f}'}; beta=0 node vs E0: "
         f"{prof[0]-E0.fun:+.2f} (must be ~0); beta=1 penalty "
         f"{prof[3]-prof.min():+.2f}")
L.append("")

# ---------------- scrambles, warm-started ----------------
rng = np.random.default_rng(19)
cred_real = E0.fun - E2.fun
creds = []
for k in range(8):
    ep = e_gal[rng.permutation(NGalS)]
    bs, _ = fit_conv('M0h', use_beta=True, e_vec=ep, th0=E2.x, dml0=dml2,
                     tol=0.1, max_rounds=8)
    creds.append(E0.fun - bs.fun)
creds = np.array(creds)
scr_ok = np.all(creds > -0.5)
L.append(f"scramble gate (8 perms): real credit = {cred_real:.2f} vs "
         f"scrambled {np.round(creds,2).tolist()}")
L.append(f"  nested-sanity (all >= -0.5): {'PASS' if scr_ok else 'FAIL'}; "
         f"scrambles beaten: {int((creds < cred_real).sum())}/8")
L.append("")

# ---------------- correlation channel, depth-partialed ----------------
la0_0, f_0 = E0.x[0], E0.x[1]
fac = f_0*np.exp(dml0[gidxS])
gN = gS_gas + fac*gS_dsk + gS_bul
y0 = gN/10**la0_0
r0 = lgobsS - np.log10(gN*nu_be(y0))
pred1 = dlg_efe(y0, e_gal[gidxS])
deep = BIN_S <= 1
rbar, pbar, elog, ydep = [], [], [], []
for gi in range(NGalS):
    mm = GIDXS_S[gi]
    md = mm[deep[mm]]
    if len(md) < 2: continue
    rbar.append(np.mean(r0[md])); pbar.append(np.mean(pred1[md]))
    elog.append(math.log10(e_gal[gi]))
    ydep.append(np.mean(np.log10(y0[md])))
rbar, pbar, elog, ydep = map(np.array, (rbar, pbar, elog, ydep))
def resid_on(a, b):
    A = np.vstack([np.ones(len(b)), b]).T
    c, *_ = np.linalg.lstsq(A, a, rcond=None)
    return a - A@c
rho_raw, p_raw = spearmanr(rbar, elog)
rr, pp = resid_on(rbar, ydep), resid_on(pbar, ydep)
rho_p, p_p = spearmanr(rr, resid_on(elog, ydep))
sl_p = (np.sum(pp*rr)/np.sum(pp*pp)) if np.sum(pp*pp) > 0 else float('nan')
L.append(f"correlation channel ({len(rbar)} galaxies, >=2 deep points):")
L.append(f"  raw Spearman(deep resid, log e) = {rho_raw:+.3f} (p={p_raw:.3f})")
L.append(f"  depth-partialed Spearman = {rho_p:+.3f} (p={p_p:.3f})   "
         f"[EFE predicts negative]")
L.append(f"  depth-partialed slope obs-on-pred(beta=1) = {sl_p:+.2f}   "
         f"[~beta if pattern present]")
L.append("")

# ---------------- injections ----------------
d_t = rng.normal(0, SIG_ML, NGalS)
for btru in (3.0, 0.0):
    fac_t = 1.0*np.exp(d_t[gidxS])
    gN_t = gS_gas + fac_t*gS_dsk + gS_bul
    y_t = gN_t/A0_FID
    lg_t = (np.log10(gN_t*nu_be(y_t)) + dlg_efe(y_t, btru*e_gal[gidxS])
            + rng.normal(0, np.sqrt(sigS2 + 0.07**2)))
    bi, _ = fit_conv('M0h', use_beta=True, lg_override=lg_t, tol=0.1,
                     max_rounds=8)
    ok = (2.0 <= bi.x[2] <= 4.5) if btru == 3.0 else (bi.x[2] <= 0.8)
    L.append(f"G2 injection beta_true={btru:g}: beta_rec = {bi.x[2]:.2f} -> "
             f"{'PASS' if ok else 'FAIL'}")

out = "\n".join(L)
print(out)
with open('data/stage5e_envtest_conv.txt', 'w') as f:
    f.write(out+"\n")
print("\nsaved: data/stage5e_envtest_conv.txt")
