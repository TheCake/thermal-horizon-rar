"""
STAGE 7B: the bump hunt -- what is the point-level scatter excess at x ~ 1?

The 4W bump (survives vertical + disk/bulge M/L + sharpness + environment
marginalization; point-level within curves) now ALSO blocks the 7A gamma
contest (it occupies the quantum-vs-classical discriminant window x 0.4-1.2).
This stage runs the pre-registered DISCRIMINATION MATRIX on candidates:

  mundane candidates (all baryonic-structured):
    - spiral-arm/bar streaming + Renzo-rule kinematic features (coherent
      wiggles at 1-3 disk scale lengths, strongest in star-dominated disks)
    - radial M/L gradients / local population structure (star-dominated only)
    - bulge-disk decomposition shape errors (bulge galaxies only)
  law-like alternative: transition-localized scatter of the relation itself
    (composition-BLIND, white, organized by x not by R/R_d).

Three orthogonal discriminators, all at FIXED x (the bump window), so the
mechanical x-composition correlation is frozen out:

  B1 COMPOSITION: within W = [0.8,1.4] (x_fid, the 4W window, pre-registered),
     split points by fiducial star fraction f_* = (g_disk+g_bul)/g_N terciles.
     Excess e_S = s^2_W(S) - mean[s^2_C1(S), s^2_C2(S)] per tercile class
     (controls C1 = [0.4,0.8), C2 = (1.4,2.2], same class). D = e_top - e_bot.
     Baryonic candidates: D > 0. Bath/law-like: D = 0 (composition-blind).
     Galaxy-bootstrap (250) for the CI.
  B2 COHERENCE: within-curve lag-1 correlation of adjacent-R residual pairs,
     both points in W vs both in C1 u C2 (per-galaxy offsets removed first).
     Streaming/wiggles: rho_W > 0 (coherent over adjacent beams), rho_C
     smaller. Per-point draw / error underestimate: white everywhere.
     Null: 500 within-curve permutations.
  B3 GEOMETRY: does the excess organize by x or by R/R_d? Renzo structure
     lives at R/R_d ~ 1.5-3 at ANY x; the law's window is x-organized. Probe:
     the mid-R/R_d excess inside the LOW-x control slice C1, vs the window
     excess.

Gates:
  G0 regression: the frozen mean law (EQ refit) reproduces 7A EQ within 1.0.
  G2 calibrated splits (3 seeds each on the real design, floor 0.101, bump
     amplitude 0.083 = 7A's own EGB fit):
     (a) composition-BLIND bump injected -> |D| must NOT fake a signal
         (|mean D| < 0.5 x 0.083^2 = 0.0034) and e_all recovers within x1.5;
     (b) star-linked bump injected (top tercile only) -> D recovers > 0.0034;
     (c) white injection -> no seed gives rho_W p < 0.02 (coherence null).
Pre-registered verdicts:
  MUNDANE-BARYONIC: B1 fires star-linked (P_boot(D<=0) < 0.023) OR B2 fires
    coherent (p_W < 0.023 with p_C > 0.16).
  LAW-LIKE LEAN: B1 bounds |D| (CI inside +-0.0034) AND B2 null in both
    windows AND B3 x-organized (C1 mid-R/Rd excess < 0.25 x window excess).
  MIXED / UNRESOLVED otherwise. Honest note: LAW-LIKE here stays a LEAN --
  a composition-blind white x-organized excess still needs a positive
  identification before any bath-fluctuation language.
Writes data/stage7b_bumphunt.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)
REF_7A_EQ = -8363.17
BUMP2_REF = 0.083**2          # 7A EGB in-window added variance (raw grade)

# ---------------- loader (4T lineage + T, Rdisk, R per point) ---------------
meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta[t[0]] = (float(t[5]), int(t[17]), int(t[1]), float(t[11]))
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc, rdisk_a, type_a = \
    [], [], [], [], [], [], [], [], []
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, Ty, Rd = meta.get(name, (0, 3, -1, 0.0))
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
        gobs.append(Vo*Vo/R*KPC); sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi); R_kpc.append(R); rdisk_a.append(Rd); type_a.append(Ty)
(g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc, rdisk_a, type_a) = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc, rdisk_a, type_a))
sig2 = sig*sig
lgobs = np.log10(gobs)
NP = len(gobs)

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))

x_fid = np.sqrt((g_gas + g_dsk + g_bul)/A0_FID)
fstar = (g_dsk + g_bul)/(g_gas + g_dsk + g_bul)
rrd = np.where(rdisk_a > 0, R_kpc/np.where(rdisk_a > 0, rdisk_a, 1.0), -1.0)

W  = (x_fid >= 0.8) & (x_fid <= 1.4)
C1 = (x_fid >= 0.4) & (x_fid < 0.8)
C2 = (x_fid > 1.4) & (x_fid <= 2.2)
CC = C1 | C2

# ---------------- G0: frozen mean law (EQ refit, must match 7A) -------------
def shape_EQ(x):
    return np.exp(-np.minimum(np.clip(x, 1e-9, None), 80.0))
def m2ll_eq(th):
    la0, f, lnN, sf = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5) or \
       not (-2 < lnN < 14) or not (0 <= sf < 0.5): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    x = np.sqrt(gN/a0)
    s2 = shape_EQ(x)/(math.exp(lnN)*LN10*LN10) + sf*sf
    se2 = sig2 + s2
    r = lgobs - np.log10(gN*nu_be(gN/a0))
    return np.sum(r*r/se2 + np.log(se2))
best = None
for t0 in ([math.log10(A0_FID), 1.0, math.log(30.0), 0.05],
           [math.log10(A0_FID), 1.0, math.log(8.0), 0.10]):
    b = minimize(m2ll_eq, t0, method='Nelder-Mead',
                 options=dict(maxiter=12000, xatol=1e-6, fatol=1e-6))
    if best is None or b.fun < best.fun: best = b
LA0, FML = best.x[0], best.x[1]
g0 = abs(best.fun - REF_7A_EQ) < 1.0
gN0 = g_gas + FML*g_dsk + g_bul
RES = lgobs - np.log10(gN0*nu_be(gN0/10**LA0))

# per-galaxy offset-subtracted residuals (weights 1/(sig2+0.01))
RESO = RES.copy()
wts = 1.0/(sig2 + 0.01)
for gi in np.unique(gal_id):
    m = gal_id == gi
    RESO[m] = RES[m] - np.sum(RES[m]*wts[m])/np.sum(wts[m])

def s_ml(res, s2obs):
    """1-param ML intrinsic scatter for a cell."""
    if len(res) < 8: return np.nan
    def nll(s):
        se2 = s2obs + s*s
        return np.sum(res*res/se2 + np.log(se2))
    return minimize_scalar(nll, bounds=(1e-4, 0.5), method='bounded').x

L = [f"STAGE 7B bump hunt: {kept} galaxies, {NP} points",
     f"G0 frozen mean: EQ refit {best.fun:.2f} vs 7A {REF_7A_EQ:.2f} -> "
     f"{'PASS' if g0 else 'FAIL'}  (la0 {LA0:.4f}, f {FML:.3f})",
     f"windows (x_fid): W [0.8,1.4] n={W.sum()}  C1 [0.4,0.8) n={C1.sum()}  "
     f"C2 (1.4,2.2] n={C2.sum()}", ""]

# fine profile (deciles), raw + offset-subtracted
QD = np.quantile(x_fid, np.linspace(0, 1, 11)); QD[0], QD[-1] = 0, np.inf
L.append("fine profile (deciles): median-x | s_raw | s_offsub")
for b_ in range(10):
    m = (x_fid >= QD[b_]) & (x_fid < QD[b_+1])
    L.append(f"  {np.median(x_fid[m]):6.3f} | {s_ml(RES[m], sig2[m]):.4f} | "
             f"{s_ml(RESO[m], sig2[m]):.4f}")
L.append("")

# ---------------- B1 composition ----------------
TQ = np.quantile(fstar[W], [1/3, 2/3])
cls = np.digitize(fstar, TQ)          # 0 = gas-rich, 2 = star-dominated
def excess(mask_cls, res):
    sw = s_ml(res[W & mask_cls], sig2[W & mask_cls])
    c1 = s_ml(res[C1 & mask_cls], sig2[C1 & mask_cls])
    c2 = s_ml(res[C2 & mask_cls], sig2[C2 & mask_cls])
    base = np.nanmean([None if np.isnan(c1) else c1*c1,
                       None if np.isnan(c2) else c2*c2][0:2] if False else
                      [v*v for v in (c1, c2) if not np.isnan(v)])
    return sw*sw - base, sw, c1, c2
def bstat(res, idx=None):
    """D = e_top - e_bot on a galaxy subset (None = all)."""
    if idx is None:
        sel = np.ones(NP, bool)
    else:
        sel = np.isin(gal_id, idx)
    out = []
    for c in (2, 0):
        m = sel & (cls == c)
        sw = s_ml(res[W & m], sig2[W & m])
        cs = [s_ml(res[cw & m], sig2[cw & m]) for cw in (C1, C2)]
        cs = [v*v for v in cs if not np.isnan(v)]
        if np.isnan(sw) or not cs: return np.nan
        out.append(sw*sw - np.mean(cs))
    return out[0] - out[1]

e_top, sw_t, c1_t, c2_t = excess(cls == 2, RES)
e_bot, sw_b, c1_b, c2_b = excess(cls == 0, RES)
e_mid, sw_m, _, _ = excess(cls == 1, RES)
e_all, sw_a, c1_a, c2_a = excess(np.ones(NP, bool), RES)
D_obs = e_top - e_bot
L.append(f"B1 COMPOSITION (f_* terciles at fixed x; thresholds "
         f"{TQ[0]:.3f}/{TQ[1]:.3f}):")
L.append(f"  all:  s_W {sw_a:.4f} vs C {c1_a:.4f}/{c2_a:.4f}  "
         f"excess e = {e_all:+.5f}")
L.append(f"  star-dom (top): s_W {sw_t:.4f} (C {c1_t:.4f}/{c2_t:.4f}) "
         f"e = {e_top:+.5f}   n_W = {int((W & (cls==2)).sum())}")
L.append(f"  mid:            s_W {sw_m:.4f}                    "
         f"e = {e_mid:+.5f}   n_W = {int((W & (cls==1)).sum())}")
L.append(f"  gas-rich (bot): s_W {sw_b:.4f} (C {c1_b:.4f}/{c2_b:.4f}) "
         f"e = {e_bot:+.5f}   n_W = {int((W & (cls==0)).sum())}")
rng = np.random.default_rng(77)
ugal = np.unique(gal_id)
Db = []
for _ in range(250):
    Db.append(bstat(RES, rng.choice(ugal, len(ugal), replace=True)))
Db = np.array([d for d in Db if not np.isnan(d)])
pD = np.mean(Db <= 0)
L.append(f"  D = e_top - e_bot = {D_obs:+.5f}  boot CI "
         f"[{np.percentile(Db,2.5):+.5f}, {np.percentile(Db,97.5):+.5f}]  "
         f"P(D<=0) = {pD:.3f}  ({len(Db)} reps)")
L.append("")

# ---------------- B2 coherence ----------------
sbar2 = np.nanmean([v*v for v in (sw_a,) if not np.isnan(v)])
Z = RESO/np.sqrt(sig2 + sbar2)
order = np.lexsort((R_kpc, gal_id))
def rho_pairs(z, win):
    num = den1 = den2 = 0.0; n = 0
    oz = z[order]; og = gal_id[order]; ow = win[order]
    for i in range(len(oz)-1):
        if og[i] != og[i+1]: continue
        if not (ow[i] and ow[i+1]): continue
        num += oz[i]*oz[i+1]; den1 += oz[i]**2; den2 += oz[i+1]**2; n += 1
    return (num/math.sqrt(den1*den2) if n > 10 else np.nan), n
rho_W, nW = rho_pairs(Z, W)
rho_C, nC = rho_pairs(Z, CC)
null_W, null_C = [], []
rng2 = np.random.default_rng(78)
Zs = Z.copy()
for _ in range(500):
    for gi in ugal:
        m = gal_id == gi
        Zs[m] = rng2.permutation(Z[m])
    null_W.append(rho_pairs(Zs, W)[0]); null_C.append(rho_pairs(Zs, CC)[0])
null_W = np.array(null_W); null_C = np.array(null_C)
p_W = np.mean(null_W >= rho_W); p_C = np.mean(null_C >= rho_C)
L.append(f"B2 COHERENCE (adjacent-R pairs, offset-subtracted):")
L.append(f"  rho_W = {rho_W:+.4f} ({nW} pairs)  perm p = {p_W:.4f}   "
         f"null [{np.percentile(null_W,2.5):+.4f},{np.percentile(null_W,97.5):+.4f}]")
L.append(f"  rho_C = {rho_C:+.4f} ({nC} pairs)  perm p = {p_C:.4f}")
L.append("")

# ---------------- B3 geometry ----------------
rd_ok = rrd > 0
def s2cell(res, m):
    v = s_ml(res[m], sig2[m]); return v*v if not np.isnan(v) else np.nan
L.append("B3 GEOMETRY (excess by R/R_d within x-slices; s^2 cells):")
for nmx, msk in (('C1 [0.4,0.8)', C1), ('W  [0.8,1.4]', W)):
    row = []
    for lo, hi in ((0, 1.5), (1.5, 3.0), (3.0, 99.0)):
        m = msk & rd_ok & (rrd >= lo) & (rrd < hi)
        row.append(f"R/Rd {lo:.0f}-{hi:.0f}: {s2cell(RES, m):.5f} "
                   f"(n={int(m.sum())})")
    L.append(f"  {nmx}: " + " | ".join(row))
mid_c1 = C1 & rd_ok & (rrd >= 1.5) & (rrd < 3.0)
out_c1 = C1 & rd_ok & ((rrd < 1.5) | (rrd >= 3.0))
exc_c1_mid = s2cell(RES, mid_c1) - s2cell(RES, out_c1)
ratio = exc_c1_mid/e_all if e_all and not np.isnan(exc_c1_mid) else np.nan
L.append(f"  C1 mid-R/Rd excess = {exc_c1_mid:+.5f} = {ratio:.2f} x window "
         f"excess ({e_all:+.5f})")
L.append("")

# descriptive splits
tsp = type_a <= 6
L.append("descriptive: in-window excess by class (e = s2_W - mean s2_C):")
for nm, msk in (('type<=6 (Sa-Scd)', tsp), ('type>=7 (Sd-BCD)', ~tsp),
                ('bulge gal', g_bul > 0), ('bulgeless', ~(g_bul > 0))):
    gm = np.ones(NP, bool) & msk
    sw = s_ml(RES[W & gm], sig2[W & gm])
    cs = [s_ml(RES[cw & gm], sig2[cw & gm]) for cw in (C1, C2)]
    cs = [v*v for v in cs if not np.isnan(v)]
    e = sw*sw - np.mean(cs) if cs and not np.isnan(sw) else np.nan
    L.append(f"  {nm:18s}: e = {e:+.5f}  (n_W = {int((W & gm).sum())})")
L.append("")

# ---------------- G2 calibrated split instrument ----------------
lg_true = np.log10(gN0*nu_be(gN0/10**LA0))
lg_save = lgobs.copy(); res_save = RES.copy()
g2a_D, g2a_e, g2b_D, g2c_p = [], [], [], []
for k in range(3):
    rngi = np.random.default_rng(801 + k)
    for mode in ('blind', 'star', 'white'):
        s2i = np.full(NP, 0.101**2)
        if mode == 'blind': s2i[W] += BUMP2_REF
        if mode == 'star':  s2i[W & (cls == 2)] += BUMP2_REF
        lgobs = lg_true + rngi.normal(0, np.sqrt(sig2 + s2i))
        RESi = lgobs - lg_true
        if mode == 'blind':
            g2a_D.append(bstat(RESi));
            swl = s_ml(RESi[W], sig2[W])
            csl = [s_ml(RESi[cw], sig2[cw]) for cw in (C1, C2)]
            g2a_e.append(swl*swl - np.mean([v*v for v in csl]))
        elif mode == 'star':
            g2b_D.append(bstat(RESi))
        else:
            RESio = RESi.copy()
            for gi in ugal:
                m = gal_id == gi
                RESio[m] -= np.sum(RESi[m]*wts[m])/np.sum(wts[m])
            Zi = RESio/np.sqrt(sig2 + 0.101**2)
            rW = rho_pairs(Zi, W)[0]
            g2c_p.append(np.mean(null_W >= rW))
lgobs = lg_save; RES = res_save
g2a = abs(np.mean(g2a_D)) < 0.5*BUMP2_REF and \
      (1/1.5 < np.mean(g2a_e)/BUMP2_REF < 1.5)
g2b = np.mean(g2b_D) > 0.5*BUMP2_REF
g2c = min(g2c_p) > 0.02
L.append(f"G2a blind-bump: D = {np.round(g2a_D,5).tolist()} "
         f"(|mean| < {0.5*BUMP2_REF:.5f}?), e recov = "
         f"{np.round(g2a_e,5).tolist()} vs {BUMP2_REF:.5f} -> "
         f"{'PASS' if g2a else 'FAIL'}")
L.append(f"G2b star-bump:  D = {np.round(g2b_D,5).tolist()} -> "
         f"{'PASS' if g2b else 'FAIL'}")
L.append(f"G2c white rho null: p = {np.round(g2c_p,3).tolist()} -> "
         f"{'PASS' if g2c else 'FAIL'}")
L.append("")

# ---------------- verdict ----------------
b1_star = (pD < 0.023) and g2a and g2b
b2_coh = (p_W < 0.023) and (p_C > 0.16) and g2c
b1_blind = (abs(D_obs) < 0.5*BUMP2_REF and
            np.percentile(Db, 2.5) > -BUMP2_REF and
            np.percentile(Db, 97.5) < BUMP2_REF) and g2a and g2b
b2_null = (p_W > 0.16) and g2c
b3_xorg = (not np.isnan(ratio)) and (ratio < 0.25)
if b1_star or b2_coh:
    verdict = "MUNDANE-BARYONIC (source: " + \
        ("composition-linked " if b1_star else "") + \
        ("coherent-wiggles" if b2_coh else "") + ")"
elif b1_blind and b2_null and b3_xorg:
    verdict = "LAW-LIKE LEAN (composition-blind + white + x-organized)"
else:
    verdict = "MIXED/UNRESOLVED"
L.append(f"VERDICT (pre-registered): {verdict}")
L.append("caveats: offset subtraction deflates all cells ~(1-1/n_g) "
         "uniformly; thin cells (n<8) = NaN; LAW-LIKE stays a LEAN pending "
         "positive identification.")

out = "\n".join(L)
print(out)
with open('data/stage7b_bumphunt.txt', 'w') as f:
    f.write(out+"\n")
