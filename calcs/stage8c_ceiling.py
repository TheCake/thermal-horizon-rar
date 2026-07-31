"""STAGE 8C — THE CEILING TEST: p <= 3/4, parameter-free (P1's
in-sample clause).  Pre-reg a5b4816 (bars locked before this script
was written).

The 5M vertical-hardened hier machinery (global (a0, f_ML, s_int,
dlt) + per-galaxy dml (0.1 dex prior) + per-galaxy dv (measured
sigma_v priors) + Mistele lensing leg) with the nu_p tail dial
(5G convention, == BE at p = 1/2).
LAYER I: pooled shared-p profile + 25-rep galaxy bootstrap + the
Chae e_N median-split arms (ceiling-only language, correction #14).
LAYER II: per-galaxy p_g census with local (dml, dv) re-profiling,
quotability and the calibrated exceedance criterion.
G-INJ: two near-ceiling null skies (p_true = 0.72) + one power sky
(p_true = 0.90).
Output: data/stage8c_ceiling.txt
"""
import csv
import glob
import math
import os
import time
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25
LN10 = math.log(10)
S_ML = 0.1*LN10
COMP_5M_BE = -12152.49          # 5M dv-ON BE (G1 regression target)

T0 = time.time()
L = []
def P(s):
    print(f"[{(time.time()-T0)/60:6.1f}m] {s}", flush=True)
    L.append(s)

# ---------------- data (5M loader + names) ----------------
meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines)
            if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18:
        continue
    try:
        meta[t[0]] = (float(t[5]), int(t[17]), float(t[2]), float(t[3]),
                      float(t[6]))
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
sigv_g_map, name_map = {}, {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                        recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 10.0, 1.0, 3.0))
    if inc < 30 or q > 2:
        continue
    kept += 1
    irad = math.radians(inc)
    sv = math.hypot((eD/max(D, 1e-3))/LN10,
                    2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10)
    sigv_g_map[gi] = max(sv, 0.01)
    name_map[gi] = name
    for l in open(path):
        if l.startswith('#'):
            continue
        t = l.split()
        if len(t) < 6:
            continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10:
            continue
        gg = Vg*abs(Vg)/R*KPC
        gd = UD*Vd*abs(Vd)/R*KPC
        gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0:
            continue
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
GNAME = [name_map[g] for g in ug]
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

# Chae join (6I convention)
chae = {}
with open('data/chae2021_table3.csv') as f:
    for row in csv.DictReader(l for l in f if not l.startswith('#')):
        chae[row['galaxy']] = float(row['log_eN_maxclust'])
EN = np.full(NGal, np.nan)
for i in range(NGal):
    if GNAME[i] in chae:
        EN[i] = chae[GNAME[i]]
matched = np.isfinite(EN)
en_med = np.median(EN[matched])
ARM_LO = matched & (EN <= en_med)
ARM_HI = matched & (EN > en_med)
ARM_UN = ~matched

# ---------------- the tail dial ----------------
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))

def nu_p_pt(y, pexp):
    yc = np.clip(np.asarray(y, float), 1e-14, None)
    u = np.minimum(yc**pexp, 60.0)
    return (1.0 - np.exp(-u))**(-0.5/pexp)

yt = np.logspace(-6, 2, 300)
assert np.max(np.abs(nu_p_pt(yt, np.full(300, 0.5)) - nu_be(yt))) < 1e-10
P("G1a identity nu_p(1/2) == BE on y in [1e-6, 100] -> PASS")

# ---------------- objective ----------------
def m2c(th, Ppt, dml, dv, sv, w_g, p_len):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5):
        return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8:
        return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*nu_p_pt(gN/a0, Ppt)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm_) - dv[gidx]
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    rl = l_gobs[lmask] - (lg + np.log10(
        nu_p_pt(10**lg/a0, np.full(int(lmask.sum()), p_len))))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    out += np.sum(w_g*dv*dv/(sv*sv))
    return out

def fit_c(P_g, sv, w_g=None, free_p=False, p_len=0.65, tol=0.1,
          max_rounds=12, th0=None, lite=False, lgo_use=None):
    """5M's fit_v generalized: P_g = per-galaxy tail array (modified
    in place if free_p); returns (best, dml, dv, P_g)."""
    global lgobs
    lgobs_save = lgobs
    if lgo_use is not None:
        lgobs = lgo_use
    try:
        if w_g is None:
            w_g = np.ones(NGal)
        dml = np.zeros(NGal)
        dv = np.zeros(NGal)
        best, prev = None, None
        mi = 2000 if lite else 4000
        for rd in range(max_rounds):
            starts = ([list(best.x)] if best is not None else []) + \
                     ([list(th0)] if th0 is not None and best is None
                      else []) + \
                     ([[math.log10(A0_FID), 1.0, 0.08, 0.0]]
                      if rd == 0 else [])
            bb = None
            for t0 in starts:
                b = minimize(lambda t: m2c(t, P_g[gidx], dml, dv, sv,
                                           w_g, p_len), t0,
                             method='Nelder-Mead',
                             options=dict(maxiter=mi, xatol=1e-6,
                                          fatol=1e-7))
                if bb is None or b.fun < bb.fun:
                    bb = b
            best = bb
            la0, f, s_int, dlt = best.x
            se2c = s_int*s_int
            a0 = 10**la0
            for _ in range(2 if lite else 3):
                fac = f*np.exp(dml[gidx])
                gN = g_gas + fac*g_dsk + g_bul
                r0 = lgobs - np.log10(gN*nu_p_pt(gN/a0, P_g[gidx]))
                for gi2 in range(NGal):
                    mm = GIDXS[gi2]
                    w = 1.0/(sig2[mm] + se2c)
                    dv[gi2] = np.sum(w*r0[mm])/(np.sum(w)
                                                + 1.0/sv[gi2]**2)
                for gi2 in range(NGal):
                    mm = GIDXS[gi2]
                    def od(dl):
                        fc = f*math.exp(dl)
                        gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                        rr = (lgobs[mm]
                              - np.log10(gN2*nu_p_pt(gN2/a0,
                                                     np.full(len(mm),
                                                             P_g[gi2])))
                              - dv[gi2])
                        s2 = sig2[mm] + se2c
                        return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                    dml[gi2] = minimize_scalar(
                        od, bounds=(-0.7, 0.7), method='bounded').x
                if free_p:
                    for gi2 in range(NGal):
                        mm = GIDXS[gi2]
                        fc = f*math.exp(dml[gi2])
                        gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                        s2 = sig2[mm] + se2c
                        def op(pv):
                            rr = (lgobs[mm]
                                  - np.log10(gN2*nu_p_pt(gN2/a0,
                                             np.full(len(mm), pv)))
                                  - dv[gi2])
                            return np.sum(rr*rr/s2)
                        P_g[gi2] = minimize_scalar(
                            op, bounds=(0.30, 1.50), method='bounded').x
            cur = m2c(best.x, P_g[gidx], dml, dv, sv, w_g, p_len)
            if prev is not None and abs(prev - cur) < tol:
                break
            prev = cur
        b = minimize(lambda t: m2c(t, P_g[gidx], dml, dv, sv, w_g,
                                   p_len), list(best.x),
                     method='Nelder-Mead',
                     options=dict(maxiter=mi, xatol=1e-6, fatol=1e-7))
        if b.fun < best.fun:
            best = b
        return best, dml, dv, P_g
    finally:
        lgobs = lgobs_save

P(f"data: {kept} galaxies, {len(gobs)} points + {int(lmask.sum())} "
  f"lensing; Chae matched {int(matched.sum())}/{NGal} (arms "
  f"{int(ARM_LO.sum())}/{int(ARM_HI.sum())}/{int(ARM_UN.sum())} at "
  f"median log eN = {en_med:.3f})")

# ---------------- G1b: pinned-1/2 regression vs 5M ----------------
bg, _, _, _ = fit_c(np.full(NGal, 0.5), SIGV, p_len=0.5, tol=0.1,
                    max_rounds=15)
ok1 = abs(bg.fun - COMP_5M_BE) < 1.5
P(f"G1b pinned p = 1/2 everywhere: {bg.fun:.2f} (5M dv-ON BE "
  f"{COMP_5M_BE}) -> {'PASS' if ok1 else 'FAIL'}")
assert ok1

# ---------------- LAYER I: pooled profile ----------------
PGRID = np.round(np.arange(0.40, 1.101, 0.05), 3)
prof, th_w = [], None
for pnode in PGRID:
    b, _, _, _ = fit_c(np.full(NGal, pnode), SIGV, p_len=pnode,
                       tol=0.15, max_rounds=8, th0=th_w, lite=True)
    prof.append(b.fun)
    th_w = b.x
prof = np.array(prof)
i0 = int(np.argmin(prof))
p_hat = PGRID[i0]
if 0 < i0 < len(PGRID)-1:
    x3, y3 = PGRID[i0-1:i0+2], prof[i0-1:i0+2]
    c2_, c1_, _ = np.polyfit(x3, y3, 2)
    if c2_ > 0:
        p_hat = float(-c1_/(2*c2_))
lo = hi = None
for j in range(i0, -1, -1):
    if prof[j] > prof[i0]+1.0:
        lo = float(np.interp(prof[i0]+1.0, [prof[j+1], prof[j]],
                             [PGRID[j+1], PGRID[j]]))
        break
for j in range(i0, len(PGRID)):
    if prof[j] > prof[i0]+1.0:
        hi = float(np.interp(prof[i0]+1.0, [prof[j-1], prof[j]],
                             [PGRID[j-1], PGRID[j]]))
        break
g2_pool = prof[i0] <= bg.fun + 0.5
P(f"LAYER I pooled: p_hat = {p_hat:.3f} (D1 {lo}..{hi}; "
  f"{'INTERIOR' if 0 < i0 < len(PGRID)-1 else 'EDGE'}); "
  f"min -2lnL = {prof[i0]:.2f} [G2 nesting vs pinned-1/2 "
  f"{bg.fun:.2f}: {'OK' if g2_pool else 'FAIL'}]")
P("  profile: " + " ".join(f"{p:.2f}:{v-prof[i0]:+.1f}"
                           for p, v in zip(PGRID, prof)))

# ---------------- LAYER I: galaxy bootstrap (25 reps) ----------------
rng = np.random.default_rng(31)
COARSE = [0.5, 0.6, 0.7, 0.8, 0.9]
boot_p = []
for k in range(25):
    pick = rng.choice(NGal, NGal, replace=True)
    wg = np.bincount(pick, minlength=NGal).astype(float)
    vals = {}
    thw = list(th_w)
    for pn in COARSE:
        b, _, _, _ = fit_c(np.full(NGal, pn), SIGV, w_g=wg, p_len=pn,
                           tol=0.3, max_rounds=5, th0=thw, lite=True)
        vals[pn] = b.fun
        thw = b.x
    pc = min(vals, key=vals.get)
    for pn in (pc-0.05, pc-0.025, pc+0.025, pc+0.05):
        pn = round(min(max(pn, 0.40), 1.10), 3)
        if pn not in vals:
            b, _, _, _ = fit_c(np.full(NGal, pn), SIGV, w_g=wg,
                               p_len=pn, tol=0.3, max_rounds=5,
                               th0=thw, lite=True)
            vals[pn] = b.fun
    ks = sorted(vals)
    vs = [vals[k2] for k2 in ks]
    im = int(np.argmin(vs))
    ph = ks[im]
    if 0 < im < len(ks)-1:
        c2_, c1_, _ = np.polyfit(ks[im-1:im+2], vs[im-1:im+2], 2)
        if c2_ > 0:
            ph = float(-c1_/(2*c2_))
    boot_p.append(ph)
    if (k+1) % 5 == 0:
        P(f"  bootstrap {k+1}/25 done")
boot_p = np.array(boot_p)
b16, b50, b84 = np.percentile(boot_p, (16, 50, 84))
sig_boot = 0.5*(b84-b16)
P(f"LAYER I bootstrap: p 16/50/84 = {b16:.3f}/{b50:.3f}/{b84:.3f} "
  f"(sigma_boot = {sig_boot:.3f}); P(p > 0.75) = "
  f"{(boot_p > 0.75).mean():.3f}")

# ---------------- LAYER I: arms ----------------
arm_out = {}
for tag, msk in (('eN-LOW', ARM_LO), ('eN-HIGH', ARM_HI),
                 ('unmatched', ARM_UN)):
    gsel = np.where(msk)[0]
    psel = np.isin(gidx, gsel)
    # subset views via weights (galaxies outside get weight 0; the
    # lensing leg is EXCLUDED by zeroing dlt freedom is kept but the
    # leg is population-level -- disclosed in pre-reg; implement by
    # weighting only: lensing stays but is shared across arms)
    wg = np.where(msk, 1.0, 0.0)
    vals = {}
    thw = list(th_w)
    for pn in PGRID[::2]:
        b, _, _, _ = fit_c(np.full(NGal, pn), SIGV, w_g=wg, p_len=pn,
                           tol=0.3, max_rounds=5, th0=thw, lite=True)
        vals[pn] = b.fun
        thw = b.x
    ks = sorted(vals)
    vs = [vals[k2] for k2 in ks]
    im = int(np.argmin(vs))
    ph = ks[im]
    if 0 < im < len(ks)-1:
        c2_, c1_, _ = np.polyfit(ks[im-1:im+2], vs[im-1:im+2], 2)
        if c2_ > 0:
            ph = float(-c1_/(2*c2_))
    # D1 interval on the coarse grid
    alo = ahi = None
    for j in range(im, -1, -1):
        if vs[j] > vs[im]+1.0:
            alo = float(np.interp(vs[im]+1.0, [vs[j+1], vs[j]],
                                  [ks[j+1], ks[j]]))
            break
    for j in range(im, len(ks)):
        if vs[j] > vs[im]+1.0:
            ahi = float(np.interp(vs[im]+1.0, [vs[j-1], vs[j]],
                                  [ks[j-1], ks[j]]))
            break
    arm_out[tag] = (ph, alo, ahi)
    P(f"LAYER I arm {tag} ({int(msk.sum())} gal): p_hat = {ph:.3f} "
      f"(D1 {alo}..{ahi}) [ceiling-only read; no ordering claim]")

# ---------------- LAYER II machinery ----------------
SCAN = np.round(np.arange(0.30, 1.501, 0.05), 3)

def census(lgo_use, tag):
    """Joint free-p fit + per-galaxy scans; returns per-galaxy
    (p_hat, halfwidth, quotable, exceed, d75)."""
    Pg = np.full(NGal, 0.65)
    b, dml, dv, Pg = fit_c(Pg, SIGV, free_p=True, p_len=p_hat,
                           tol=0.15, max_rounds=8, th0=th_w,
                           lite=True, lgo_use=lgo_use)
    la0, f, s_int, dlt = b.x
    a0 = 10**la0
    se2c = s_int*s_int
    out = []
    lg_loc = lgobs if lgo_use is None else lgo_use
    for gi2 in range(NGal):
        mm = GIDXS[gi2]
        s2 = sig2[mm] + se2c
        sv2 = SIGV[gi2]**2
        curve = []
        for pv in SCAN:
            dvv, dl = dv[gi2], dml[gi2]
            for _ in range(2):
                fc = f*math.exp(dl)
                gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                r0 = lg_loc[mm] - np.log10(
                    gN2*nu_p_pt(gN2/a0, np.full(len(mm), pv)))
                w = 1.0/s2
                dvv = np.sum(w*r0)/(np.sum(w) + 1.0/sv2)
                def od(dl2):
                    fc2 = f*math.exp(dl2)
                    gN3 = g_gas[mm] + fc2*g_dsk[mm] + g_bul[mm]
                    rr = (lg_loc[mm] - np.log10(
                        gN3*nu_p_pt(gN3/a0, np.full(len(mm), pv)))
                        - dvv)
                    return (np.sum(rr*rr/s2) + dl2*dl2/(S_ML*S_ML)
                            + dvv*dvv/sv2)
                dl = minimize_scalar(od, bounds=(-0.7, 0.7),
                                     method='bounded').x
            curve.append(od(dl))
        curve = np.array(curve)
        im = int(np.argmin(curve))
        ph = SCAN[im]
        if 0 < im < len(SCAN)-1:
            c2_, c1_, _ = np.polyfit(SCAN[im-1:im+2],
                                     curve[im-1:im+2], 2)
            if c2_ > 0:
                ph = float(np.clip(-c1_/(2*c2_), SCAN[0], SCAN[-1]))
        d = curve - curve[im]
        wlo = whi = None
        for j in range(im, -1, -1):
            if d[j] > 1.0:
                wlo = float(np.interp(1.0, [d[j+1], d[j]],
                                      [SCAN[j+1], SCAN[j]]))
                break
        for j in range(im, len(SCAN)):
            if d[j] > 1.0:
                whi = float(np.interp(1.0, [d[j-1], d[j]],
                                      [SCAN[j-1], SCAN[j]]))
                break
        interior = wlo is not None and whi is not None
        hw = 0.5*(whi-wlo) if interior else np.inf
        quot = interior and hw <= 0.30
        d75 = float(np.interp(0.75, SCAN, d))
        exceed = quot and ph > 0.75 and d75 >= 9.0
        out.append((ph, hw, quot, exceed, d75))
    nq = sum(1 for o in out if o[2])
    ne = sum(1 for o in out if o[3])
    P(f"  census[{tag}]: quotable {nq}/{NGal}, exceedances {ne}")
    return out, b

P("LAYER II: the per-galaxy census (real sky)")
cen_real, b_real = census(None, 'real')

# ---------------- G-INJ ----------------
la0, f, s_int, dlt = b_real.x
a0 = 10**la0

def synth(p_true, seed):
    r2 = np.random.default_rng(seed)
    dml_t = r2.normal(0, S_ML, NGal)
    dv_t = r2.normal(0, SIGV)
    gN_t = g_gas + f*np.exp(dml_t[gidx])*g_dsk + g_bul
    mu = np.log10(gN_t*nu_p_pt(gN_t/a0, np.full(len(gobs), p_true)))
    eps = r2.normal(0, np.sqrt(sig2 + s_int*s_int))
    return mu + dv_t[gidx] + eps

inj = {}
for p_true, seed, tag in ((0.72, 101, 'null-A'), (0.72, 202, 'null-B'),
                          (0.90, 303, 'power')):
    P(f"G-INJ sky {tag} (p_true = {p_true})")
    cen, _ = census(synth(p_true, seed), tag)
    nq = sum(1 for o in cen if o[2])
    ne = sum(1 for o in cen if o[3])
    med = np.median([o[0] for o in cen if o[2]]) if nq else np.nan
    inj[tag] = (nq, ne, med)
    P(f"  -> quotable {nq}, exceedances {ne}, median p_hat "
      f"{med:.3f}")

false_ok = inj['null-A'][1] <= 1 and inj['null-B'][1] <= 1
nq_pow = inj['power'][0]
powered = nq_pow > 0 and inj['power'][1] >= max(1, nq_pow//3)
P(f"G-INJ: false-exceedance {'PASS' if false_ok else 'FAIL'} "
  f"(null skies {inj['null-A'][1]}/{inj['null-B'][1]}); power "
  f"{'PASS' if powered else 'FAIL -> census POWER-LIMITED, no bar'} "
  f"({inj['power'][1]}/{nq_pow} fired at p_true = 0.90)")

# ---------------- verdict ----------------
nq_r = sum(1 for o in cen_real if o[2])
exc = [(GNAME[i], o[0], o[4]) for i, o in enumerate(cen_real) if o[3]]
arm_break = any(ph - 0.75 > 2*sig_boot for ph, _, _ in
                arm_out.values())
pool_break = (b50 - 0.75) > 2*sig_boot
cen_break = powered and false_ok and (
    len(exc) >= 2 or any(d >= 16.0 for _, _, d in exc))
holds_arms = all(ph <= 0.75 + sig_boot for ph, _, _ in
                 arm_out.values()) and b50 <= 0.75 + sig_boot
cen_clean = (not powered) or (false_ok and len(exc) == 0)

P("")
for i, o in enumerate(cen_real):
    if o[3]:
        P(f"  EXCEEDANCE: {GNAME[i]} p_hat = {o[0]:.3f}, "
          f"D(0.75) = {o[4]:.1f}")
if pool_break or arm_break or cen_break:
    v = "CEILING-BROKEN"
elif holds_arms and cen_clean and nq_r >= 30:
    v = "CEILING-HOLDS"
else:
    v = "AMBIG"
P(f"==> 8C VERDICT (pre-registered bars): {v}")
P(f"    pooled p = {b50:.3f} +/- {sig_boot:.3f} (profile "
  f"{p_hat:.3f}, D1 {lo}..{hi}); arms " +
  ", ".join(f"{t}:{v_[0]:.3f}" for t, v_ in arm_out.items()) +
  f"; census quotable {nq_r}, exceedances {len(exc)}, powered = "
  f"{powered}")

with open('data/stage8c_ceiling.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8c_ceiling.txt")
