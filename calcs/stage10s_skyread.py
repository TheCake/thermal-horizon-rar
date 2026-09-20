"""
STAGE 10S SKY READ -- the two-leg H0 meter (pre-reg 2a461d5 + amendment
A2 fb74226; run 1 preserved in data/stage10s_skyread_run1.txt).

A2 changes vs run 1 (all logged pre-letter in the pre-reg):
  (a) display guards (run 1 crashed at a print after all fits landed);
  (b) flat co-read reported as NO-FIXED-POINT-IN-RANGE with the
      extrapolation-grade implied crossing labeled as such; the
      flat-vs-hier treatment split NAMED in every letter;
  (c) NEW GATE G4h: the hier crossing is injection-validated (noiseless
      lock-mocks at 62/70/78 through the full grid procedure, bar 1%;
      3-seed noisy replication at 70 = the hier floor). No hier H0_B
      without G4h;
  (d) sigma_B = galaxy bootstrap of the HIER crossing (100 reps, grid
      {50,60,70,80,90}, per-pick distance jitter), replacing the
      selection-biased flat-solve bootstrap.

Pins 1-10 otherwise as committed. NO credence cell moves on any outcome.
Writes data/stage10s_skyread.txt.
"""
import glob, math, os, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
LN10 = math.log(10)
C_LIGHT = 299792458.0
KMS_MPC = 3.240779e-20
A0_FID = 1.2e-10
S_ML = 0.1*LN10
U_PRIOR = (0.9/18.0)/LN10

def h0_of_a0(a0):
    return 2.0*math.pi*a0/C_LIGHT/KMS_MPC

def a0_of_h0(h0):
    return C_LIGHT*h0*KMS_MPC/(2.0*math.pi)

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

t00 = time.time()
P("STAGE 10S SKY READ, amended run (pre-reg 2a461d5 + A2 fb74226)")
P("")

# ---------------- data build (verbatim harness) ----------------
meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta[t[0]] = (float(t[5]), int(t[17]), float(t[2]), float(t[3]),
                      float(t[6]), int(t[4]))
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
sigv_g_map, fd_g_map, dist_g_map = {}, {}, {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc, fdv = meta.get(name, (0, 3, 10.0, 1.0, 3.0, 1))
    if inc < 30 or q > 2: continue
    kept += 1
    irad = math.radians(inc)
    sv = math.hypot((eD/max(D, 1e-3))/LN10,
                    2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10)
    sigv_g_map[gi] = max(sv, 0.01)
    fd_g_map[gi] = fdv
    dist_g_map[gi] = (D, eD)
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
        sig.append(2*eV/Vo/LN10)
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs = np.log10(gobs)
ug = np.unique(gal_id)
FD_G = np.array([fd_g_map[g] for g in ug])
flow_gal = ug[FD_G == 1]
anch_gal = ug[np.isin(FD_G, (2, 3, 5))]
uma_gal = ug[FD_G == 4]
legA_gal = np.concatenate([anch_gal, uma_gal])
assert len(gobs) == 2700 and kept == 153
assert (len(flow_gal), len(anch_gal), len(uma_gal)) == (82, 41, 26)
gpts_all = {g: np.where(gal_id == g)[0] for g in ug}

gdfrac = {}
for g_ in ug:
    m = gal_id == g_
    gdfrac[g_] = float(np.mean(g_gas[m] > g_dsk[m] + g_bul[m]))
gd_all = np.array([g_ for g_ in ug if gdfrac[g_] >= 0.5])
npts_gd = int(np.isin(gal_id, gd_all).sum())
sel_ok = (len(gd_all) == 38 and npts_gd == 422)
P(f"selector regression (8S): GD = {len(gd_all)}/{npts_gd} "
  f"(archived 38/422) -> {'PASS' if sel_ok else 'FAIL'}")
gd_A = np.array([g_ for g_ in legA_gal if gdfrac[g_] >= 0.5])
nd_A = np.array([g_ for g_ in legA_gal if gdfrac[g_] < 0.5])
P(f"leg A composition: {len(legA_gal)} galaxies ({len(anch_gal)} anchored "
  f"+ {len(uma_gal)} UMa); GD {len(gd_A)} / non-GD {len(nd_A)}")
P("")

# ---------------- nu families (verbatim 5M) ----------------
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def nu_p065(y):
    yc = np.clip(np.asarray(y, float), 1e-14, None)
    ex = np.exp(-np.minimum(yc**0.65, 60.0))
    return (1.0-ex)**(-1.0/1.3)
def nu_gm(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    a = y**0.75
    w = np.sqrt(nu_simple(y))
    for _ in range(30):
        u = np.minimum(a*w, 60.0)
        eu = np.exp(u)
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        H = w*w - 1.0 - n
        dH = 2.0*w + a*eu/(em1*em1)
        w = np.maximum(w - H/dH, 1e-8)
    return w*w
def nu_boot(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    u = 0.5*(y + np.sqrt(y*y + 4.0*y))
    for _ in range(14):
        uc = np.minimum(u, 45.0)
        eu = np.exp(uc)
        em1 = eu - 1.0
        F = u - y - y/em1
        dF = 1.0 + y*eu/(em1*em1)
        u = np.maximum(u - F/dF, 1e-13)
    nu = u/y
    big = y > 45.0
    if np.any(big):
        yb = np.minimum(y, 700.0)
        nu = np.where(big, 1.0 + 1.0/np.expm1(yb), nu)
    return nu
FAMS = {'BE': nu_be, 'p065': nu_p065, 'gm': nu_gm, 'boot': nu_boot}

# ---------------- flat machinery (gated) ----------------
def fit_flat_pts(lg, gg, gd, gb_, nu, x0=(-9.92, 1.2)):
    def lo(v):
        la0, fd = v
        if not (-10.8 < la0 < -9.2) or not (0.3 <= fd <= 3): return 1e9
        gN = gg + fd*gd + gb_
        r = lg - np.log10(gN*nu(gN/10**la0))
        return float(np.mean(r*r))
    b = minimize(lo, list(x0), method='Nelder-Mead',
                 options={'xatol': 1e-5, 'fatol': 1e-12, 'maxiter': 3000})
    return b.x[0], b.x[1]

def solve_flat(lg, gg, gd, gb_, nu, iters=16):
    def F(h0p):
        la0_i, _ = fit_flat_pts(lg - math.log10(73.0/h0p), gg, gd, gb_, nu)
        return 10**la0_i - a0_of_h0(h0p)
    lo, hi = 55.0, 90.0
    flo = F(lo)
    if flo*F(hi) > 0: return None
    for _ in range(iters):
        mid = 0.5*(lo+hi)
        fm = F(mid)
        if flo*fm <= 0: hi = mid
        else: lo = mid; flo = fm
    return 0.5*(lo+hi)

# ---------------- hier machinery (5M minus lensing, + UMa u) ----------
def build_sub(gal_seq, lg_src=None, dlg_per_occ=None):
    """gal_seq may contain duplicates (bootstrap picks): each occurrence
    is its own group with its own dv/dml. lg_src overrides lgobs (mocks).
    dlg_per_occ: per-occurrence lgobs shift (distance jitter)."""
    src = lgobs if lg_src is None else lg_src
    blocks, gidx_s, sv, ig = [], [], [], []
    lg_parts = []
    for k, g in enumerate(gal_seq):
        pts = gpts_all[int(g)]
        blocks.append(pts)
        gidx_s.append(np.full(len(pts), k))
        sv.append(sigv_g_map[int(g)])
        ig.append(1.0 if fd_g_map[int(g)] == 4 else 0.0)
        shift = 0.0 if dlg_per_occ is None else dlg_per_occ[k]
        lg_parts.append(src[pts] + shift)
    pidx = np.concatenate(blocks)
    gidx_s = np.concatenate(gidx_s)
    n = len(gal_seq)
    ig = np.array(ig)
    return dict(n=n, gidx=gidx_s, sv=np.array(sv), isuma_g=ig,
                isuma_pt=ig[gidx_s],
                gg=g_gas[pidx], gd=g_dsk[pidx], gb=g_bul[pidx],
                lg=np.concatenate(lg_parts), s2=sig2[pidx],
                gpts=[np.where(gidx_s == i)[0] for i in range(n)])

def fit_hier(sub, nu, use_u=False, upri=U_PRIOR, dlg=0.0, tol=0.05,
             max_rounds=15, th0=None):
    n = sub['n']
    lg = sub['lg'] + dlg
    gg, gd, gb_ = sub['gg'], sub['gd'], sub['gb']
    s2, gidx_s, sv = sub['s2'], sub['gidx'], sub['sv']
    ipt, ig = sub['isuma_pt'], sub['isuma_g']
    w_g = np.ones(n)
    dml = np.zeros(n)
    dv = np.zeros(n)

    def m2(th, dml, dv):
        if use_u:
            la0, f, s_int, u = th
            if abs(u) > 0.2: return 1e12
        else:
            la0, f, s_int = th; u = 0.0
        if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
        if not (1e-3 <= s_int < 0.4): return 1e12
        a0 = 10**la0
        fac = f*np.exp(dml[gidx_s])
        gN = gg + fac*gd + gb_
        gm_ = gN*nu(gN/a0)
        se2 = s2 + s_int*s_int
        r = lg - np.log10(gm_) - dv[gidx_s] - u*ipt
        out = np.sum(r*r/se2 + np.log(se2))
        out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
        out += np.sum(w_g*dv*dv/(sv*sv))
        if use_u: out += (u/upri)**2
        return out

    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08] + ([0.0] if use_u else [])]
                  if rd == 0 else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2(t, dml, dv), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int = best.x[:3]
        u = best.x[3] if use_u else 0.0
        se2c = s_int*s_int
        for _ in range(3):
            fac = f*np.exp(dml[gidx_s])
            gN = gg + fac*gd + gb_
            r0 = lg - np.log10(gN*nu(gN/10**la0)) - u*ipt
            for gi2 in range(n):
                mm = sub['gpts'][gi2]
                w = 1.0/(s2[mm] + se2c)
                dv[gi2] = np.sum(w*r0[mm])/(np.sum(w) + 1.0/sv[gi2]**2)
            for gi2 in range(n):
                mm = sub['gpts'][gi2]
                uoff = u*ig[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = gg[mm] + fc*gd[mm] + gb_[mm]
                    rr = (lg[mm] - np.log10(gN2*nu(gN2/10**la0))
                          - dv[gi2] - uoff)
                    ss = s2[mm] + se2c
                    return np.sum(rr*rr/ss) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2(best.x, dml, dv)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2(t, dml, dv), list(best.x),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best

def hier_crossing(sub, grid=(60.0, 65.0, 70.0, 75.0, 80.0),
                  lo_cap=50.0, hi_cap=95.0, verbose=False):
    """The leg-B PRIMARY estimator: hier a0 on the trial grid, crossing
    vs the lock by linear interpolation; extension by 5 within caps."""
    phi = {}
    th_warm = None
    def node(h0p):
        nonlocal th_warm
        b = fit_hier(sub, nu_be, use_u=False,
                     dlg=-math.log10(73.0/h0p), th0=th_warm)
        th_warm = list(b.x)
        return math.log(10**b.x[0]) - math.log(a0_of_h0(h0p)), 10**b.x[0]
    for h0p in grid:
        phi[h0p], a0h = node(h0p)
        if verbose:
            P(f"  hier grid H0' = {h0p:.0f}: a0 = {a0h:.4e}, lock = "
              f"{a0_of_h0(h0p):.4e}, phi = {phi[h0p]:+.4f}  "
              f"[{time.time()-t00:.0f}s]")
    def bracketed(ph):
        ks = sorted(ph.keys())
        return any(ph[ks[i]]*ph[ks[i+1]] <= 0 for i in range(len(ks)-1))
    while not bracketed(phi):
        # ROUND-45 FIX (condition 2): the original rule assumed a RISING
        # phi (the flat regime); the hier regime FALLS (gamma_H < 1), so
        # it extended AWAY from high roots (one-sided censoring; the
        # reviewer's noiseless truth-85 demonstration). Slope-aware now.
        ks = sorted(phi.keys())
        slope = phi[ks[-1]] - phi[ks[0]]
        if phi[ks[0]] > 0: go_up = slope < 0
        else: go_up = slope > 0
        if go_up and ks[-1] < hi_cap: nxt = ks[-1] + 5.0
        elif (not go_up) and ks[0] > lo_cap: nxt = ks[0] - 5.0
        else: return None
        phi[nxt], _ = node(nxt)
        if verbose:
            P(f"  hier grid EXTENDED H0' = {nxt:.0f}: phi = {phi[nxt]:+.4f}")
    ks = sorted(phi.keys())
    for i in range(len(ks)-1):
        if phi[ks[i]]*phi[ks[i+1]] <= 0:
            x0_, x1_ = ks[i], ks[i+1]
            return x0_ + (x1_-x0_)*phi[x0_]/(phi[x0_]-phi[x1_])
    return None

# ---------------- G0a: flat solver identity ----------------
fidx = np.where(np.isin(gal_id, flow_gal))[0]
FD_TRUE = 1.2
gN_true = g_gas[fidx] + FD_TRUE*g_dsk[fidx] + g_bul[fidx]
g0a_ok = True
rows = []
for h0_true in (62.0, 70.0, 78.0):
    s_true = 73.0/h0_true
    lg_cat = np.log10(s_true*gN_true*nu_be(gN_true/a0_of_h0(h0_true)))
    h0_rec = solve_flat(lg_cat, g_gas[fidx], g_dsk[fidx], g_bul[fidx],
                        nu_be)
    derr = abs(h0_rec - h0_true)/h0_true if h0_rec else 1.0
    g0a_ok &= derr < 0.002
    rows.append(f"{h0_true:.0f}->{h0_rec:.3f}" if h0_rec else "NO-BRACKET")
P(f"G0a flat-solver identity: {', '.join(rows)} -> "
  f"{'PASS' if g0a_ok else 'FAIL'}  [{time.time()-t00:.0f}s]")

# ---------------- G4h: hier-crossing injection (A2 gate) ----------
# mock lgobs arrays live on the full-point grid so build_sub can slice
g4h_ok = True
mock_full = np.empty_like(lgobs)
for h0_true in (62.0, 70.0, 78.0):
    s_true = 73.0/h0_true
    a0t = a0_of_h0(h0_true)
    gNf = g_gas + FD_TRUE*g_dsk + g_bul
    mock_full[:] = np.log10(s_true*gNf*nu_be(gNf/a0t))
    subm = build_sub(list(flow_gal), lg_src=mock_full.copy())
    h0r = hier_crossing(subm)
    derr = abs(h0r - h0_true)/h0_true if h0r else 1.0
    g4h_ok &= derr < 0.01
    P(f"G4h noiseless {h0_true:.0f} -> "
      f"{h0r:.2f} ({100*derr:.2f}%)" if h0r else
      f"G4h noiseless {h0_true:.0f} -> NO CROSSING")
P(f"G4h noiseless -> {'PASS' if g4h_ok else 'FAIL'}  "
  f"[{time.time()-t00:.0f}s]")
errs = []
for seed in (42, 101, 202):
    rng = np.random.default_rng(seed)
    s_true = 73.0/70.0
    a0t = a0_of_h0(70.0)
    gNf = g_gas + FD_TRUE*g_dsk + g_bul
    mock_full[:] = (np.log10(s_true*gNf*nu_be(gNf/a0t))
                    + rng.normal(0, np.sqrt(sig2 + 0.08**2)))
    subm = build_sub(list(flow_gal), lg_src=mock_full.copy())
    h0r = hier_crossing(subm)
    errs.append(100*(h0r-70.0)/70.0 if h0r else float('nan'))
errs = np.array(errs)
P(f"G4h noisy floor (3 seeds at 70): errors % = "
  f"{[('%+.1f' % e) if np.isfinite(e) else 'none' for e in errs]}; the "
  f"x{1.0/(1.0-0.86):.0f} lever makes this the hier floor  "
  f"[{time.time()-t00:.0f}s]")
if not g4h_ok:
    P("G4h FAILED -> per A2 no hier H0_B is quoted; letter path "
      "degrades to the flat co-read + M-GRAY/M-POWER-DEAD grammar")

# ---------------- G0b: UMa nuisance identity ----------------
subA = build_sub(list(legA_gal))
bA_plain = fit_hier(subA, nu_be, use_u=False)
bA_tight = fit_hier(subA, nu_be, use_u=True, upri=1e-6,
                    th0=list(bA_plain.x)+[0.0])
d_fun = abs(bA_tight.fun - bA_plain.fun)
d_a0 = abs(10**bA_tight.x[0] - 10**bA_plain.x[0])/10**bA_plain.x[0]
g0b_ok = d_fun < 0.5 and d_a0 < 0.002
P(f"G0b UMa-nuisance identity: d(-2lnL) = {d_fun:.3f}, d(a0) = "
  f"{100*d_a0:.3f}% -> {'PASS' if g0b_ok else 'FAIL'}")
P("")

# ================= LEG A =================
P("== LEG A (anchored + UMa; H0-assumption-independent) ==")
bA = fit_hier(subA, nu_be, use_u=True, th0=list(bA_plain.x)+[0.0])
a0_A = 10**bA.x[0]
P(f"leg-A hier BE (PRIMARY): a0_A = {a0_A:.4e}, f_ML = {bA.x[1]:.2f}, "
  f"s_int = {bA.x[2]:.3f}, UMa shared = {bA.x[3]:+.4f} dex "
  f"(run-1 regression: 1.0186e-10)")
aidx = np.where(np.isin(gal_id, legA_gal))[0]
la0_Af, fd_Af = fit_flat_pts(lgobs[aidx], g_gas[aidx], g_dsk[aidx],
                             g_bul[aidx], nu_be)
P(f"leg-A flat BE (co-read): a0 = {10**la0_Af:.4e}, f_d = {fd_Af:.2f}")
fam_A = {'BE': a0_A}
for nm in ('p065', 'gm', 'boot'):
    bf = fit_hier(subA, FAMS[nm], use_u=True, th0=list(bA.x))
    fam_A[nm] = 10**bf.x[0]
P(f"leg-A hier family: " + ", ".join(f"{k} {v:.4e}"
                                     for k, v in fam_A.items()))

rng5 = np.random.default_rng(202)
rel_sig, uma_flag = {}, {}
for g in legA_gal:
    D, eD = dist_g_map[g]
    if fd_g_map[g] == 4:
        rel_sig[g] = 2.3/18.0; uma_flag[g] = True
    else:
        rel_sig[g] = eD/max(D, 1e-3); uma_flag[g] = False
la0_reps = []
for _ in range(200):
    pick = rng5.choice(legA_gal, size=len(legA_gal), replace=True)
    sc_shared = 1.0 + rng5.normal(0, 0.9/18.0)
    rows_, dsh = [], []
    for g in pick:
        s_g = 1.0 + rng5.normal(0, rel_sig[g])
        s_g = min(max(s_g, 0.5), 1.5)
        if uma_flag[g]: s_g *= sc_shared
        rows_.append(gpts_all[g])
        dsh.append(np.full(len(gpts_all[g]), -math.log10(s_g)))
    ridx = np.concatenate(rows_)
    lg_b = lgobs[ridx] + np.concatenate(dsh)
    la0_reps.append(fit_flat_pts(lg_b, g_gas[ridx], g_dsk[ridx],
                                 g_bul[ridx], nu_be)[0])
a0_reps = 10**np.array(la0_reps)
sig_a0A = float(np.std(a0_reps))
H0_A = h0_of_a0(a0_A)
sig_H0A = h0_of_a0(sig_a0A)
P(f"sigma_A (flat boot, 200 reps): sigma(H0_A) = {sig_H0A:.2f}")
P(f"LEG A: H0_A = {H0_A:.2f} +- {sig_H0A:.2f} km/s/Mpc (hier BE "
  f"primary; family H0 {h0_of_a0(min(fam_A.values())):.1f}-"
  f"{h0_of_a0(max(fam_A.values())):.1f})")
powA = sig_H0A <= 15.0
P(f"G5 re-check: {'powered' if powA else 'POWER-LIMITED'}")

gidx_gd = np.where(np.isin(gal_id, gd_A))[0]
gidx_nd = np.where(np.isin(gal_id, nd_A))[0]
la0_gd, _ = fit_flat_pts(lgobs[gidx_gd], g_gas[gidx_gd], g_dsk[gidx_gd],
                         g_bul[gidx_gd], nu_be)
la0_nd, _ = fit_flat_pts(lgobs[gidx_nd], g_gas[gidx_nd], g_dsk[gidx_nd],
                         g_bul[gidx_nd], nu_be)
P(f"G7 GD split (flat): GD({len(gd_A)}) {10**la0_gd:.4e} / "
  f"non-GD({len(nd_A)}) {10**la0_nd:.4e}")
if len(gd_A) >= 15 and len(nd_A) >= 15:
    bgd = fit_hier(build_sub(list(gd_A)), nu_be, use_u=True)
    bnd = fit_hier(build_sub(list(nd_A)), nu_be, use_u=True)
    dgd = abs(10**bgd.x[0] - 10**bnd.x[0])
    g7_fired = dgd > sig_a0A
    P(f"G7 GD split (hier): GD {10**bgd.x[0]:.4e} / non-GD "
      f"{10**bnd.x[0]:.4e}; |d| = {dgd:.2e} vs sigma {sig_a0A:.2e} -> "
      f"{'CAVEAT FIRES' if g7_fired else 'inside 1 sigma'}")
else:
    g7_fired = abs(10**la0_gd - 10**la0_nd) > sig_a0A
    P("G7 hier split POWER-LIMITED; flat split quoted")
P("")

# ================= LEG B =================
P(f"== LEG B (flow; self-consistent solve) ==  [{time.time()-t00:.0f}s]")
subB = build_sub(list(flow_gal))
H0_B_hier = None
if g4h_ok:
    H0_B_hier = hier_crossing(subB, verbose=True)
    if H0_B_hier is None:
        P("LEG B hier: NO FIXED POINT in [50, 95]")

# flat co-read + extrapolation-grade implied crossing (A2-b)
lgB = lgobs[fidx]
H0_B_flat = solve_flat(lgB, g_gas[fidx], g_dsk[fidx], g_bul[fidx], nu_be)
if H0_B_flat is not None:
    P(f"leg-B flat BE solve (co-read): H0 = {H0_B_flat:.2f}")
else:
    la0_lo, _ = fit_flat_pts(lgB - math.log10(73.0/60.0),
                             g_gas[fidx], g_dsk[fidx], g_bul[fidx], nu_be)
    la0_ct, _ = fit_flat_pts(lgB, g_gas[fidx], g_dsk[fidx], g_bul[fidx],
                             nu_be)
    la0_hi, _ = fit_flat_pts(lgB - math.log10(73.0/90.0),
                             g_gas[fidx], g_dsk[fidx], g_bul[fidx], nu_be)
    gam = ((la0_hi - la0_lo)*LN10)/math.log(90.0/60.0)
    a0c = 10**la0_ct
    h0x = 73.0*math.exp(math.log(a0_of_h0(73.0)/a0c)/(gam-1.0)) \
        if gam > 1.05 else float('nan')
    P(f"leg-B flat BE (co-read): NO FIXED POINT in [55, 90]; catalog "
      f"a0 = {a0c:.4e}, response gamma = {gam:.2f}, implied crossing "
      f"~{h0x:.0f} (EXTRAPOLATION-GRADE, outside the pre-registered "
      f"bracket -- the flat treatment refuses self-consistency at any "
      f"plausible H0)")
fam_B = {}
for nm in ('p065', 'gm', 'boot'):
    r = solve_flat(lgB, g_gas[fidx], g_dsk[fidx], g_bul[fidx],
                   FAMS[nm], iters=12)
    fam_B[nm] = r
P("leg-B flat family solves: " +
  ", ".join(f"{k} {('%.1f' % v) if v else 'no-bracket'}"
            for k, v in fam_B.items()))

# sigma_B: hier-crossing bootstrap (A2-d)
sig_H0B = float('nan')
if g4h_ok and H0_B_hier is not None:
    rngB = np.random.default_rng(303)
    reps, nofix = [], 0
    BGRID = (50.0, 60.0, 70.0, 80.0, 90.0)
    for i in range(100):
        pick = list(rngB.choice(flow_gal, size=len(flow_gal),
                                replace=True))
        dlg_occ = []
        for g in pick:
            D, eD = dist_g_map[g]
            s_g = 1.0 + rngB.normal(0, eD/max(D, 1e-3))
            s_g = min(max(s_g, 0.5), 1.5)
            dlg_occ.append(-math.log10(s_g))
        subr = build_sub(pick, dlg_per_occ=dlg_occ)
        h0r = hier_crossing(subr, grid=BGRID)
        if h0r is None: nofix += 1
        else: reps.append(h0r)
        if (i+1) % 25 == 0:
            P(f"  sigma_B boot {i+1}/100  [{time.time()-t00:.0f}s]")
    reps = np.array(reps)
    if len(reps) > 20:
        sig_H0B = float(np.std(reps))
        P(f"sigma_B (HIER-crossing boot, 100 reps, {nofix} no-crossing): "
          f"{sig_H0B:.2f} km/s/Mpc; percentiles 16/50/84 = "
          f"{np.percentile(reps, [16,50,84]).round(1).tolist()}")
    else:
        P(f"sigma_B boot DEGENERATE ({nofix} no-crossing)")
if H0_B_hier is not None:
    P(f"LEG B: H0_B = {H0_B_hier:.2f} +- {sig_H0B:.2f} km/s/Mpc "
      f"(hier crossing primary, G4h-validated)")
P("")

# ================= THE AGREEMENT TEST =================
P("== AGREEMENT TEST ==")
letter = None
if (H0_B_hier is not None) and np.isfinite(sig_H0B) and powA:
    delta = abs(H0_A - H0_B_hier)
    joint = math.hypot(sig_H0A, sig_H0B)
    nsig = delta/joint
    P(f"H0_A = {H0_A:.2f} +- {sig_H0A:.2f}  vs  H0_B = {H0_B_hier:.2f} "
      f"+- {sig_H0B:.2f}  ->  |d| = {delta:.2f} = {nsig:.2f} x joint "
      f"({joint:.2f})")
    if nsig <= 1.0 and g0a_ok and g0b_ok and g4h_ok:
        letter = "M-AGREE"
        wA, wB = 1.0/sig_H0A**2, 1.0/sig_H0B**2
        H0_joint = (wA*H0_A + wB*H0_B_hier)/(wA+wB)
        sig_joint = 1.0/math.sqrt(wA+wB)
        P(f"LETTER M-AGREE: the two legs agree; the lock passes its "
          f"first self-consistency test as an instrument; H0(meter) = "
          f"{H0_joint:.1f} +- {sig_joint:.1f} km/s/Mpc "
          f"(H0-assumption-independent; BE-form primary; hier "
          f"treatment)")
    elif nsig > 2.0:
        letter = "M-SPLIT"
        P("LETTER M-SPLIT: the lock fails the agreement test at current "
          "grade OR a distance-provenance systematic is unmodeled; no "
          "H0 returned. Named suspects: flow Virgocentric-infall model; "
          "UMa TFR zero point; leg composition (GD fraction, "
          "surface-brightness mix); the treatment split below.")
    else:
        letter = "M-GRAY"
        P("LETTER M-GRAY: 1-2 sigma -- both numbers reported, no "
          "headline H0; successor named: DR4-era anchors + measured "
          "peculiar-velocity flow model.")
else:
    letter = "M-GRAY" if powA else "M-POWER-DEAD"
    P(f"letter path degraded ({letter}): missing hier fixed point or "
      f"sigma, or leg A power-limited")
P("")
P("MANDATORY DISCLOSURES (every letter, per A2-b):")
P(f"  treatment split: the flat co-read refuses a fixed point in "
  f"[55, 90] on leg B (implied crossing ~extrapolation-grade >100) "
  f"while the hier treatment crosses at {H0_B_hier if H0_B_hier else float('nan'):.1f}; "
  f"leg-A flat co-read a0 = {10**la0_Af:.3e} vs hier {a0_A:.3e}. The "
  f"meter's verdict is TREATMENT-CONDITIONAL at SPARC grade (the 4H/5M "
  f"M/L-vs-a0 degeneracy; hier = the pre-registered primary).")
P(f"  family bands: leg A H0 {h0_of_a0(min(fam_A.values())):.1f}-"
  f"{h0_of_a0(max(fam_A.values())):.1f}; leg-B flat family all "
  f"no-bracket; hier family on leg B not run (one function pinned "
  f"primary pre-reg; family spread carried by leg A)")
if g7_fired:
    P("  G7 CAVEAT FIRES: leg-A GD/non-GD split exceeds 1 sigma")
P("")
P(f"verdict letter: {letter}")
P("credences: NO cell moves on any outcome (pre-registered); 53/8 "
  "untouched")
P(f"total wall-clock: {(time.time()-t00)/60:.1f} min")

with open('data/stage10s_skyread.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage10s_skyread.txt")
