"""
STAGE 7C: bump identification confirmed? -- the stratified test + the
UNBLOCKED gamma contest.

7B's descriptive cells found (post-hoc, hence this pre-registered
confirmatory stage): at FIXED x, inner-disk points (R < 1.5 R_d) carry
~2.4x the scatter of outer points in BOTH x-slices, and the window-vs-C1
excess vanishes within R/R_d strata. Hypothesis: the x~1 bump is a
RADIUS-MIX artifact -- inner-disk structure (beam smearing, bars,
non-circular motions, decomposition) entering the window, which samples
49% inner points vs 16% in the deep control.

Pre-registered bars:
  C1 IDENTIFICATION: stratified excess table (W vs controls within each
     R/R_d stratum; if a C2 cell has n < 30 its stratum baseline falls back
     to C1-only, disclosed). EXPLAINED = 1 - sum_str w_str*max(e_str,0)/e_all,
     w_str = window composition. Bar: EXPLAINED >= 0.75 -> the bump is
     IDENTIFIED as inner-disk radius mix (the 4W "point-level mystery"
     resolves mundane). 0.4-0.75 partial; < 0.4 not explained.
  C2 CLEAN-SUBSET BUMP: EGB refit (7A model) on R/Rd >= 1.5 points only:
     bump amplitude b_clean < 0.04 (half of 7A's 0.083) -> bump gone.
  C3 THE GAMMA CONTEST, two variants:
     variant B (PRIMARY, full data): all scatter rivals gain one shared
       inner-disk variance term c_in^2 on R < 1.5 R_d points -- EQ_i/EC_i/
       EG_i + free-gamma profile. Same 7A bars: exclusion at D(-2lnL) >= 9,
       gamma_hat interior, G2 injections (now with inner-term truth
       c_in = 0.114) must separate quantum from classical.
     variant A (robustness, clean subset R/Rd >= 1.5): same models, no
       inner term; own G2.
  G5 CORRELATION HONESTY: 7B measured lag-1 rho ~ 0.87 (smooth within-curve
     misfit; independence bars are nominally calibrated only). Control:
     every-3rd-point thinning per curve, EG_i profile D(gamma=1 vs 2)
     re-quoted. Verdict language carries "nominal-independence" always.
Outcome space pre-stated: bump identified + gamma resolves (best);
bump identified + gamma still flat (legitimate -- window cleared but
discriminant power spent); bump not explained (the mix reading dies).
Writes data/stage7c_gammaclean.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)
REF_7A = {'EQ': -8363.17, 'EC': -8362.74, 'EG': -8366.36}

meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta[t[0]] = (float(t[5]), int(t[17]), float(t[11]))
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc, rdisk_a = [], [], [], [], [], [], [], []
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, Rd = meta.get(name, (0, 3, 0.0))
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
        gal_id.append(gi); R_kpc.append(R); rdisk_a.append(Rd)
(g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc, rdisk_a) = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id, R_kpc, rdisk_a))
sig2 = sig*sig
lgobs = np.log10(gobs)
NP = len(gobs)

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))

x_fid = np.sqrt((g_gas + g_dsk + g_bul)/A0_FID)
rd_ok = rdisk_a > 0
rrd = np.where(rd_ok, R_kpc/np.where(rd_ok, rdisk_a, 1.0), -1.0)
INNER = rd_ok & (rrd < 1.5)
CLEAN = rd_ok & (rrd >= 1.5)
W  = (x_fid >= 0.8) & (x_fid <= 1.4)
C1 = (x_fid >= 0.4) & (x_fid < 0.8)
C2 = (x_fid > 1.4) & (x_fid <= 2.2)

def s_ml(res, s2obs):
    if len(res) < 8: return np.nan
    def nll(s):
        se2 = s2obs + s*s
        return np.sum(res*res/se2 + np.log(se2))
    return minimize_scalar(nll, bounds=(1e-4, 0.5), method='bounded').x

# frozen mean from the EQ fit (7A values; refit for the residuals)
def m2ll_eq(th):
    la0, f, lnN, sf = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5) or \
       not (-2 < lnN < 14) or not (0 <= sf < 0.5): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    x = np.sqrt(gN/a0)
    s2 = np.exp(-np.minimum(x, 80))/(math.exp(lnN)*LN10*LN10) + sf*sf
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
gN0 = g_gas + FML*g_dsk + g_bul
RES = lgobs - np.log10(gN0*nu_be(gN0/10**LA0))
g0 = abs(best.fun - REF_7A['EQ']) < 1.0

L = [f"STAGE 7C: {kept} galaxies, {NP} points; INNER (R<1.5Rd) n={INNER.sum()}"
     f", CLEAN n={CLEAN.sum()}, no-Rd n={(~rd_ok).sum()}",
     f"G0 frozen mean: {best.fun:.2f} vs 7A EQ {REF_7A['EQ']:.2f} -> "
     f"{'PASS' if g0 else 'FAIL'}", ""]

# ---------------- C1: stratified excess ----------------
STR = [(0.0, 1.5), (1.5, 3.0), (3.0, 99.0)]
def s2c(m):
    v = s_ml(RES[m], sig2[m]); return (v*v if not np.isnan(v) else np.nan), int(m.sum())
e_all_W, _ = s2c(W); e_c1, _ = s2c(C1); e_c2, _ = s2c(C2)
e_all = e_all_W - 0.5*(e_c1 + e_c2)
L.append("C1 stratified excess (s^2 cells; e_str = W - baseline per stratum):")
tot_resid = 0.0
for lo, hi in STR:
    st = rd_ok & (rrd >= lo) & (rrd < hi)
    sw, nw = s2c(W & st)
    sc1, n1 = s2c(C1 & st)
    sc2, n2 = s2c(C2 & st)
    if n2 >= 30 and not np.isnan(sc2):
        base = 0.5*(sc1 + sc2); btxt = 'C1+C2'
    else:
        base = sc1; btxt = 'C1-only'
    e_str = sw - base if not (np.isnan(sw) or np.isnan(base)) else np.nan
    wgt = nw/W.sum()
    if not np.isnan(e_str): tot_resid += wgt*max(e_str, 0.0)
    L.append(f"  R/Rd {lo:g}-{hi:g}: W {sw:.5f} (n={nw}) vs {btxt} {base:.5f} "
             f"(n={n1}/{n2}) -> e_str = {e_str:+.5f}  w = {wgt:.2f}")
explained = 1.0 - tot_resid/e_all if e_all > 0 else np.nan
L.append(f"  unstratified e_all = {e_all:+.5f}; residual (weighted, positive "
         f"parts) = {tot_resid:+.5f}; EXPLAINED = {explained:.2f}")
c1_ok = explained >= 0.75
L.append(f"  C1 bar (>=0.75): {'IDENTIFIED' if c1_ok else ('PARTIAL' if explained>=0.4 else 'NOT EXPLAINED')}")
L.append("")

# ---------------- shared model machinery (7A + inner term) ----------------
def fit(model, sel, use_inner, x0):
    gfix = None
    if model.startswith('EGFIX:'):
        gfix = float(model.split(':')[1]); model = 'EGFIX'
    inn = INNER[sel] if use_inner else None
    gg, gd, gb = g_gas[sel], g_dsk[sel], g_bul[sel]
    s2o, lo_ = sig2[sel], lgobs[sel]
    def obj(th):
        la0, f = th[0], th[1]
        if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
        gN = gg + f*gd + gb
        x = np.sqrt(gN/10**la0)
        k = 2
        if use_inner:
            cin = th[2]
            if not (0 <= cin < 0.4): return 1e12
            k = 3
        if model == 'E0':
            s0 = th[k]
            if not (1e-4 <= s0 < 0.5): return 1e12
            s2 = np.full(len(gN), s0*s0)
        elif model in ('EQ', 'EC'):
            lnN, sf = th[k], th[k+1]
            if not (-2 < lnN < 14) or not (0 <= sf < 0.5): return 1e12
            fac = 1.0 if model == 'EQ' else 2.0
            s2 = np.exp(-fac*np.minimum(x, 80))/(math.exp(lnN)*LN10*LN10) + sf*sf
        elif model == 'EG':
            lnS0, gam, sf = th[k], th[k+1], th[k+2]
            if not (-19 < lnS0 < -1.4) or not (0.05 < gam < 6.0) \
               or not (0 <= sf < 0.5): return 1e12
            s2 = math.exp(lnS0)*np.exp(-gam*np.minimum(x, 80)) + sf*sf
        elif model == 'EGFIX':
            lnS0, sf = th[k], th[k+1]
            if not (-19 < lnS0 < -1.4) or not (0 <= sf < 0.5): return 1e12
            s2 = math.exp(lnS0)*np.exp(-gfix*np.minimum(x, 80)) + sf*sf
        else:  # EGB
            lnS0, gam, bb, sf = th[k], th[k+1], th[k+2], th[k+3]
            if not (-19 < lnS0 < -1.4) or not (0.05 < gam < 6.0) \
               or not (0 <= bb < 0.3) or not (0 <= sf < 0.5): return 1e12
            s2 = (math.exp(lnS0)*np.exp(-gam*np.minimum(x, 80))
                  + bb*bb*np.exp(-(x-1.1)**2/(2*0.25**2)) + sf*sf)
        if use_inner:
            s2 = s2 + np.where(inn, cin*cin, 0.0)
        se2 = s2o + s2
        gm = gN*nu_be(gN/10**la0)
        r = lo_ - np.log10(gm)
        return np.sum(r*r/se2 + np.log(se2))
    bst = None
    for t0 in x0:
        b = minimize(obj, t0, method='Nelder-Mead',
                     options=dict(maxiter=15000, xatol=1e-6, fatol=1e-6))
        if bst is None or b.fun < bst.fun: bst = b
    return bst

LNS0_G = math.log(1.0/(21.0*LN10*LN10))
LA, FM = math.log10(A0_FID), 1.0

def contest(tag, sel, use_inner):
    ci = [0.11] if use_inner else []
    ki = 1 if use_inner else 0
    r = {}
    r['E0'] = fit('E0', sel, use_inner, [[LA, FM]+ci+[0.10]])
    r['EQ'] = fit('EQ', sel, use_inner, [[LA, FM]+ci+[math.log(30.0), 0.05],
                                         [LA, FM]+ci+[math.log(8.0), 0.10]])
    r['EC'] = fit('EC', sel, use_inner, [[LA, FM]+ci+[math.log(30.0), 0.05],
                                         [LA, FM]+ci+[math.log(8.0), 0.10]])
    r['EG'] = fit('EG', sel, use_inner, [[LA, FM]+ci+[LNS0_G, 1.0, 0.05],
                                         [LA, FM]+ci+[LNS0_G, 2.0, 0.05],
                                         [LA, FM]+ci+[LNS0_G+1, 0.5, 0.08]])
    gam = r['EG'].x[4] if use_inner else r['EG'].x[3]
    L.append(f"[{tag}] E0 {r['E0'].fun:.2f} | EQ {r['EQ'].fun:.2f} | "
             f"EC {r['EC'].fun:.2f} | EG {r['EG'].fun:.2f} "
             f"gamma_hat = {gam:.3f}" +
             (f" c_in = {r['EG'].x[2]:.4f}" if use_inner else ""))
    L.append(f"    EQ-EC = {r['EQ'].fun-r['EC'].fun:+.2f}")
    # profile
    GRID = np.arange(0.25, 3.501, 0.125)
    prof = []
    warm = [LA, FM]+ci+[LNS0_G, 0.05]
    for g in GRID:
        b = fit(f'EGFIX:{g}', sel, use_inner, [warm, [LA, FM]+ci+[LNS0_G, 0.10]])
        prof.append(b.fun); warm = list(b.x)
    prof = np.array(prof)
    ib = int(prof.argmin()); pmin = prof[ib]
    d1 = prof[np.argmin(np.abs(GRID-1.0))] - pmin
    d2 = prof[np.argmin(np.abs(GRID-2.0))] - pmin
    edge = ib in (0, len(GRID)-1)
    L.append(f"    profile min {GRID[ib]:.3f}{' EDGE' if edge else ''}; "
             f"D(g=1) {d1:+.2f}  D(g=2) {d2:+.2f}")
    L.append("    " + " ".join(f"{g:.2f}:{p-pmin:+.1f}"
                               for g, p in zip(GRID[::3], prof[::3])))
    return r, GRID, prof, gam, edge, d1, d2

# ---------------- C2: clean-subset bump ----------------
rEGB = fit('EGB', CLEAN, False, [[LA, FM, LNS0_G, 1.0, 0.05, 0.05],
                                 [LA, FM, LNS0_G, 2.0, 0.05, 0.05]])
b_clean = rEGB.x[4]
c2_ok = b_clean < 0.04
L.append(f"C2 clean-subset EGB: b_clean = {b_clean:.4f} (7A full: 0.083) -> "
         f"{'BUMP GONE' if c2_ok else 'BUMP PERSISTS'}")
L.append("")

# ---------------- C3 contests ----------------
L.append("C3 variant B (PRIMARY: full data + shared inner term):")
rB, GRID, profB, gamB, edgeB, d1B, d2B = contest('full+inner', np.ones(NP, bool), True)
L.append("")
L.append("C3 variant A (robustness: clean subset, no inner term):")
rA, _, profA, gamA, edgeA, d1A, d2A = contest('clean', CLEAN, False)
L.append("")

# ---------------- G2 injections (variant B design) ----------------
lg_true = np.log10(gN0*nu_be(gN0/10**LA0))
xt = np.sqrt(gN0/10**LA0)
lg_save = lgobs.copy()
recs = {'Q': [], 'C': []}
for kind in ('Q', 'C'):
    fac = 1.0 if kind == 'Q' else 2.0
    s2t = (np.exp(-fac*np.minimum(xt, 80))/(21.0*LN10*LN10) + 0.101**2
           + np.where(INNER, 0.114**2, 0.0))
    for k in range(3):
        rng = np.random.default_rng(901 + 10*k + (0 if kind == 'Q' else 5))
        lgobs = lg_true + rng.normal(0, np.sqrt(sig2 + s2t))
        binj = fit('EG', np.ones(NP, bool), True,
                   [[LA, FM, 0.11, LNS0_G, 1.0, 0.10],
                    [LA, FM, 0.11, LNS0_G, 2.0, 0.10]])
        recs[kind].append(binj.x[4])
lgobs = lg_save
mQ, mC = np.mean(recs['Q']), np.mean(recs['C'])
sep = min(recs['C']) - max(recs['Q'])
g2 = (abs(mQ-1.0) < 0.35) and (abs(mC-2.0) < 0.5) and (sep > 0.3)
L.append(f"G2 (variant B design): gammaQ = {np.round(recs['Q'],3).tolist()} "
         f"(mean {mQ:.3f}); gammaC = {np.round(recs['C'],3).tolist()} "
         f"(mean {mC:.3f}); gap = {sep:.3f} -> {'PASS' if g2 else 'FAIL'}")

# ---------------- G5 thinning control ----------------
order = np.lexsort((R_kpc, gal_id))
pos = np.empty(NP, int)
for gi in np.unique(gal_id):
    m = order[gal_id[order] == gi]
    pos[m] = np.arange(len(m))
THIN = pos % 3 == 0
rT = {}
for mdl in ('EQ', 'EC'):
    rT[mdl] = fit(mdl, THIN, True, [[LA, FM, 0.11, math.log(30.0), 0.05],
                                    [LA, FM, 0.11, math.log(8.0), 0.10]])
L.append(f"G5 thinning (every 3rd, n={THIN.sum()}): EQ-EC = "
         f"{rT['EQ'].fun-rT['EC'].fun:+.2f} (nominal-independence caveat "
         f"carried; 7B measured lag-1 rho ~ 0.87)")
L.append("")

# ---------------- verdicts ----------------
interior_B = not edgeB
verd = "UNRESOLVED"
if g2 and interior_B:
    if d2B >= 9.0 and gamB < 1.5:
        verd = "SUPPORT (quantum: classical excluded at nominal >= 9)"
    elif d1B >= 9.0 and gamB > 1.5:
        verd = "ANTI (classical: gamma=1 excluded -> pre-committed strike)"
if not g2:
    verd = "UNRESOLVED (instrument fails on this design)"
L.append(f"C3 VERDICT (variant B, pre-registered): {verd}")
L.append(f"C1/C2 IDENTIFICATION: explained = {explained:.2f}; "
         f"b_clean = {b_clean:.4f} -> " +
         ("bump IDENTIFIED as inner-disk radius mix" if (c1_ok and c2_ok)
          else ("PARTIAL identification" if (explained >= 0.4 or c2_ok)
                else "NOT explained")))
L.append("caveats: nominal-independence -2lnL everywhere (rho ~ 0.87 measured")
L.append("in 7B); inner term is x-independent first order; stratified cells")
L.append("share galaxies between W and controls (partial pairing only).")

out = "\n".join(L)
print(out)
with open('data/stage7c_gammaclean.txt', 'w') as f:
    f.write(out+"\n")
