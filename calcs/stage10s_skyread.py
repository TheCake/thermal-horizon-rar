"""
STAGE 10S SKY READ -- the two-leg H0 meter fires (licensed by the
pre-registration commit 2a461d5, author's word 2026-09-20 "let's run it").

Everything here follows PREREG-H0METER-DRAFT.md pins 1-10 + letters.
- LEG A (anchored + UMa, hier no-lensing PRIMARY with the UMa
  shared-distance nuisance; flat BE co-read; sigma_A = the G5 bootstrap
  machinery unmasked; hier family {p065, gm, boot} band; G7 GD split).
- LEG B (flow; PRIMARY central = hier no-lensing a0 on the trial grid
  H0' in {60,65,70,75,80}, crossing vs cH0'/2pi by linear interpolation,
  extension by 5 to at most [50,95]; engine + sigma_B = the gated flat
  BE fixed-point solve, 150-rep galaxy bootstrap with distance jitter;
  flat family band).
- AGREEMENT TEST: |H0_A - H0_B| vs the joint quadrature sigma; letters
  M-AGREE / M-GRAY / M-SPLIT / M-POWER-DEAD per the pre-registered
  grammar; prohibitions: no single-leg headline, no H0 without the
  family band, "H0-assumption-independent" never "ladder-independent".

Wiring identities run FIRST (G0a: the G4a noiseless mocks reproduce
62/70/78 to 0.2%; G0b: the UMa nuisance at prior width 1e-6 reproduces
the plain no-nuisance leg-A fit). GD selector inherited verbatim from
stage 8S (regression: 38 galaxies / 422 points). Implementation notes:
the .mrt e_D = 2.5 Mpc for UMa members already includes the 0.9 shared
part in quadrature (sqrt(2.3^2+0.9^2)); pin 7 keeps it as-is, so the
shared component is very mildly double-counted in the per-galaxy prior
-- conservative (wider), disclosed. Hier G7 split quoted at hier grade
only where both sides have >= 15 galaxies; flat split always printed.

NO credence cell moves on any outcome (pre-registered).
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
U_PRIOR = (0.9/18.0)/LN10        # UMa shared-distance prior, dex

def h0_of_a0(a0):
    return 2.0*math.pi*a0/C_LIGHT/KMS_MPC

def a0_of_h0(h0):
    return C_LIGHT*h0*KMS_MPC/(2.0*math.pi)

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

t00 = time.time()
P("STAGE 10S SKY READ (pre-reg 2a461d5; two legs, one agreement test)")
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
assert len(gobs) == 2700 and kept == 153, "data-build regression FAIL"
assert (len(flow_gal), len(anch_gal), len(uma_gal)) == (82, 41, 26), \
    "leg census regression FAIL"

# GD selector (verbatim 8S)
gdfrac = {}
for g_ in ug:
    m = gal_id == g_
    gdfrac[g_] = float(np.mean(g_gas[m] > g_dsk[m] + g_bul[m]))
gd_all = np.array([g_ for g_ in ug if gdfrac[g_] >= 0.5])
npts_gd = int(np.isin(gal_id, gd_all).sum())
sel_ok = (len(gd_all) == 38 and npts_gd == 422)
P(f"selector regression (8S): GD = {len(gd_all)} galaxies / {npts_gd} "
  f"points (archived 38/422) -> {'PASS' if sel_ok else 'FAIL'}")
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

# ---------------- flat machinery (gated in 10S harness) ----------------
def fit_flat(pidx, nu, dlg=0.0, x0=(-9.92, 1.2)):
    lg = lgobs[pidx] + dlg
    gg, gd, gb_ = g_gas[pidx], g_dsk[pidx], g_bul[pidx]
    def lo(v):
        la0, fd = v
        if not (-10.8 < la0 < -9.2) or not (0.3 <= fd <= 3): return 1e9
        gN = gg + fd*gd + gb_
        r = lg - np.log10(gN*nu(gN/10**la0))
        return float(np.mean(r*r))
    b = minimize(lo, list(x0), method='Nelder-Mead',
                 options={'xatol': 1e-5, 'fatol': 1e-12, 'maxiter': 3000})
    return b.x[0], b.x[1]

def fit_flat_pts(lg, gg, gd, gb_, nu, x0=(-9.92, 1.2)):
    def lo(v):
        la0, fd = v
        if not (-10.8 < la0 < -9.2) or not (0.3 <= fd <= 3): return 1e9
        gN = gg + fd*gd + gb_
        r = lg - np.log10(gN*nu(gN/10**la0))
        return float(np.mean(r*r))
    b = minimize(lo, list(x0), method='Nelder-Mead',
                 options={'xatol': 1e-5, 'fatol': 1e-12, 'maxiter': 3000})
    return b.x[0]

def solve_flat(lg, gg, gd, gb_, nu, iters=16):
    """Fixed-point solve on given point arrays (lg at the catalog frame)."""
    def F(h0p):
        la0_i = fit_flat_pts(lg - math.log10(73.0/h0p), gg, gd, gb_, nu)
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
def build_sub(gal_list):
    gset = set(int(g) for g in gal_list)
    pidx = np.where(np.isin(gal_id, list(gset)))[0]
    glist = sorted(gset)
    gmap = {g: i for i, g in enumerate(glist)}
    gidx_s = np.array([gmap[g] for g in gal_id[pidx]])
    sv = np.array([sigv_g_map[g] for g in glist])
    isuma_g = np.array([1.0 if fd_g_map[g] == 4 else 0.0 for g in glist])
    return dict(pidx=pidx, n=len(glist), gidx=gidx_s, sv=sv,
                isuma_g=isuma_g, isuma_pt=isuma_g[gidx_s],
                gg=g_gas[pidx], gd=g_dsk[pidx], gb=g_bul[pidx],
                lg=lgobs[pidx].copy(), s2=sig2[pidx],
                gpts=[np.where(gidx_s == i)[0] for i in range(len(glist))])

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
    ndim = 4 if use_u else 3

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

# ---------------- G0a: solver wiring identity ----------------
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
P(f"G0a solver identity (noiseless mocks): {', '.join(rows)} -> "
  f"{'PASS' if g0a_ok else 'FAIL'}  [{time.time()-t00:.0f}s]")

# ---------------- G0b: UMa nuisance identity ----------------
subA = build_sub(legA_gal)
bA_plain = fit_hier(subA, nu_be, use_u=False)
bA_tight = fit_hier(subA, nu_be, use_u=True, upri=1e-6,
                    th0=list(bA_plain.x)+[0.0])
d_fun = abs(bA_tight.fun - bA_plain.fun)
d_a0 = abs(10**bA_tight.x[0] - 10**bA_plain.x[0])/10**bA_plain.x[0]
g0b_ok = d_fun < 0.5 and d_a0 < 0.002
P(f"G0b UMa-nuisance identity (width 1e-6 vs plain): d(-2lnL) = "
  f"{d_fun:.3f}, d(a0) = {100*d_a0:.3f}% -> "
  f"{'PASS' if g0b_ok else 'FAIL'}  [{time.time()-t00:.0f}s]")
P("")

# ================= LEG A =================
P("== LEG A (anchored + UMa; H0-assumption-independent distances) ==")
bA = fit_hier(subA, nu_be, use_u=True, th0=list(bA_plain.x)+[0.0])
a0_A = 10**bA.x[0]
u_A = bA.x[3]
P(f"leg-A hier BE (PRIMARY): a0_A = {a0_A:.4e} m/s^2, f_ML = "
  f"{bA.x[1]:.2f}, s_int = {bA.x[2]:.3f}, UMa shared offset = "
  f"{u_A:+.4f} dex (prior width {U_PRIOR:.4f})")
la0_Af, fd_Af = fit_flat(np.where(np.isin(gal_id, legA_gal))[0], nu_be)
P(f"leg-A flat BE (co-read): a0 = {10**la0_Af:.4e}, f_d = {fd_Af:.2f}")

# family band (hier, with nuisance)
fam_A = {'BE': a0_A}
for nm in ('p065', 'gm', 'boot'):
    bf = fit_hier(subA, FAMS[nm], use_u=True, th0=list(bA.x))
    fam_A[nm] = 10**bf.x[0]
    P(f"leg-A hier {nm}: a0 = {fam_A[nm]:.4e}  "
      f"[{time.time()-t00:.0f}s]")

# sigma_A (G5 machinery, unmasked)
rng5 = np.random.default_rng(202)
rel_sig, uma_flag = {}, {}
for g in legA_gal:
    D, eD = dist_g_map[g]
    if fd_g_map[g] == 4:
        rel_sig[g] = 2.3/18.0; uma_flag[g] = True
    else:
        rel_sig[g] = eD/max(D, 1e-3); uma_flag[g] = False
gpts_all = {g: np.where(gal_id == g)[0] for g in ug}
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
                                 g_bul[ridx], nu_be))
a0_reps = 10**np.array(la0_reps)
sig_a0A = float(np.std(a0_reps))
H0_A = h0_of_a0(a0_A)
sig_H0A = h0_of_a0(sig_a0A)
P(f"sigma_A (flat boot, 200 reps): {sig_a0A:.3e} -> sigma(H0_A) = "
  f"{sig_H0A:.2f} km/s/Mpc; boot median a0 = "
  f"{np.median(a0_reps):.4e} (central grade: hier)")
P(f"LEG A: H0_A = {H0_A:.2f} +- {sig_H0A:.2f} km/s/Mpc  (hier BE "
  f"primary; family band a0 {min(fam_A.values()):.3e}-"
  f"{max(fam_A.values()):.3e} -> H0 "
  f"{h0_of_a0(min(fam_A.values())):.1f}-"
  f"{h0_of_a0(max(fam_A.values())):.1f})")
powA = sig_H0A <= 15.0
P(f"G5 re-check at sky grade: {'powered' if powA else 'POWER-LIMITED'}")

# G7: GD split
la0_gd, _ = fit_flat(np.where(np.isin(gal_id, gd_A))[0], nu_be)
la0_nd, _ = fit_flat(np.where(np.isin(gal_id, nd_A))[0], nu_be)
P(f"G7 GD split (flat): a0(GD, {len(gd_A)} gal) = {10**la0_gd:.4e}; "
  f"a0(non-GD, {len(nd_A)} gal) = {10**la0_nd:.4e}")
if len(gd_A) >= 15 and len(nd_A) >= 15:
    bgd = fit_hier(build_sub(gd_A), nu_be, use_u=True)
    bnd = fit_hier(build_sub(nd_A), nu_be, use_u=True)
    a0gd, a0nd = 10**bgd.x[0], 10**bnd.x[0]
    P(f"G7 GD split (hier): a0(GD) = {a0gd:.4e}; a0(non-GD) = "
      f"{a0nd:.4e}; |d| = {abs(a0gd-a0nd):.2e} vs sigma_A "
      f"{sig_a0A:.2e} -> {'CAVEAT FIRES' if abs(a0gd-a0nd) > sig_a0A else 'inside 1 sigma'}")
    g7_fired = abs(a0gd - a0nd) > sig_a0A
else:
    P("G7 hier split: POWER-LIMITED (side < 15 galaxies); flat split "
      "quoted with that caveat")
    g7_fired = abs(10**la0_gd - 10**la0_nd) > sig_a0A
P("")

# ================= LEG B =================
P(f"== LEG B (flow; self-consistent solve) ==  [{time.time()-t00:.0f}s]")
subB = build_sub(flow_gal)
GRID = [60.0, 65.0, 70.0, 75.0, 80.0]
phi = {}
th_warm = None
def leg_b_node(h0p, th_warm):
    b = fit_hier(subB, nu_be, use_u=False,
                 dlg=-math.log10(73.0/h0p), th0=th_warm)
    return b
for h0p in GRID:
    b = leg_b_node(h0p, th_warm)
    th_warm = list(b.x)
    a0h = 10**b.x[0]
    phi[h0p] = math.log(a0h) - math.log(a0_of_h0(h0p))
    P(f"  hier grid H0' = {h0p:.0f}: a0 = {a0h:.4e}, lock = "
      f"{a0_of_h0(h0p):.4e}, phi = {phi[h0p]:+.4f}  "
      f"[{time.time()-t00:.0f}s]")
# extend if unbracketed
gvals = sorted(phi.keys())
def bracketed(ph):
    ks = sorted(ph.keys())
    return any(ph[ks[i]]*ph[ks[i+1]] <= 0 for i in range(len(ks)-1))
while not bracketed(phi):
    ks = sorted(phi.keys())
    if phi[ks[0]] > 0 and ks[0] > 50.0: nxt = ks[0] - 5.0
    elif phi[ks[-1]] < 0 and ks[-1] < 95.0: nxt = ks[-1] + 5.0
    else: break
    b = leg_b_node(nxt, th_warm)
    phi[nxt] = math.log(10**b.x[0]) - math.log(a0_of_h0(nxt))
    P(f"  hier grid EXTENDED H0' = {nxt:.0f}: phi = {phi[nxt]:+.4f}")
H0_B_hier = None
ks = sorted(phi.keys())
for i in range(len(ks)-1):
    if phi[ks[i]]*phi[ks[i+1]] <= 0:
        x0_, x1_ = ks[i], ks[i+1]
        H0_B_hier = x0_ + (x1_-x0_)*phi[x0_]/(phi[x0_]-phi[x1_])
        break
if H0_B_hier is None:
    P("  LEG B hier: NO FIXED POINT in [50, 95] -- reported per grammar")

# flat engine solve + family band
lgB = lgobs[fidx]
H0_B_flat = solve_flat(lgB, g_gas[fidx], g_dsk[fidx], g_bul[fidx], nu_be)
P(f"leg-B flat BE solve (engine): H0 = "
  f"{H0_B_flat:.2f}" if H0_B_flat else "leg-B flat BE solve: NO BRACKET")
fam_B = {'BE': H0_B_flat}
for nm in ('p065', 'gm', 'boot'):
    fam_B[nm] = solve_flat(lgB, g_gas[fidx], g_dsk[fidx], g_bul[fidx],
                           FAMS[nm], iters=12)
    P(f"leg-B flat {nm} solve: "
      f"{'H0 = %.2f' % fam_B[nm] if fam_B[nm] else 'NO BRACKET'}  "
      f"[{time.time()-t00:.0f}s]")

# sigma_B: 150-rep galaxy bootstrap of the flat solve
rngB = np.random.default_rng(303)
h0_reps, nofix = [], 0
for _ in range(150):
    pick = rngB.choice(flow_gal, size=len(flow_gal), replace=True)
    rows_, dsh = [], []
    for g in pick:
        D, eD = dist_g_map[g]
        s_g = 1.0 + rngB.normal(0, eD/max(D, 1e-3))
        s_g = min(max(s_g, 0.5), 1.5)
        rows_.append(gpts_all[g])
        dsh.append(np.full(len(gpts_all[g]), -math.log10(s_g)))
    ridx = np.concatenate(rows_)
    lg_b = lgobs[ridx] + np.concatenate(dsh)
    h0r = solve_flat(lg_b, g_gas[ridx], g_dsk[ridx], g_bul[ridx],
                     nu_be, iters=12)
    if h0r is None: nofix += 1
    else: h0_reps.append(h0r)
h0_reps = np.array(h0_reps)
sig_H0B = float(np.std(h0_reps)) if len(h0_reps) > 10 else float('nan')
P(f"sigma_B (flat-solve boot, 150 reps, {nofix} no-bracket): "
  f"{sig_H0B:.2f} km/s/Mpc; boot median {np.median(h0_reps):.2f}; "
  f"G4b floor cross-check ~1.7  [{time.time()-t00:.0f}s]")
if H0_B_hier is not None:
    P(f"LEG B: H0_B = {H0_B_hier:.2f} +- {sig_H0B:.2f} km/s/Mpc (hier "
      f"crossing primary; flat engine {H0_B_flat if H0_B_flat else float('nan'):.2f}; family "
      f"{min(v for v in fam_B.values() if v):.1f}-"
      f"{max(v for v in fam_B.values() if v):.1f})")
P("")

# ================= THE AGREEMENT TEST =================
P("== AGREEMENT TEST ==")
if H0_B_hier is None or not powA:
    P("letter path: a leg is power-limited or fixed-point-less")
if H0_B_hier is not None:
    delta = abs(H0_A - H0_B_hier)
    joint = math.hypot(sig_H0A, sig_H0B)
    nsig = delta/joint
    P(f"H0_A = {H0_A:.2f} +- {sig_H0A:.2f}  vs  H0_B = "
      f"{H0_B_hier:.2f} +- {sig_H0B:.2f}  ->  |d| = {delta:.2f} = "
      f"{nsig:.2f} x joint sigma ({joint:.2f})")
    if nsig <= 1.0 and powA and g0a_ok and g0b_ok:
        letter = "M-AGREE"
        wA, wB = 1.0/sig_H0A**2, 1.0/sig_H0B**2
        H0_joint = (wA*H0_A + wB*H0_B_hier)/(wA+wB)
        sig_joint = 1.0/math.sqrt(wA+wB)
        P(f"LETTER M-AGREE: the two legs agree; the lock passes its "
          f"first self-consistency test as an instrument; H0(meter) = "
          f"{H0_joint:.1f} +- {sig_joint:.1f} km/s/Mpc "
          f"(H0-assumption-independent; BE-form primary)")
    elif nsig > 2.0:
        letter = "M-SPLIT"
        P("LETTER M-SPLIT: the lock fails the agreement test at current "
          "grade OR a distance-provenance systematic is unmodeled; the "
          "meter does not return an H0. Named suspects: flow "
          "Virgocentric-infall model; UMa TFR zero point; anchored-leg "
          "GD composition (G7); the function-family choice; SPARC "
          "quality-cut selection between legs.")
    else:
        letter = "M-GRAY"
        P("LETTER M-GRAY: 1-2 sigma or one leg marginal -- both numbers "
          "reported, no headline H0; successor = better anchors "
          "(DR4-era TRGB) and a flow model with measured peculiar "
          "velocities.")
    fam_span_A = (h0_of_a0(min(fam_A.values())),
                  h0_of_a0(max(fam_A.values())))
    P(f"family bands (mandatory): leg A H0 {fam_span_A[0]:.1f}-"
      f"{fam_span_A[1]:.1f}; leg B H0 "
      f"{min(v for v in fam_B.values() if v):.1f}-"
      f"{max(v for v in fam_B.values() if v):.1f} km/s/Mpc")
    if g7_fired:
        P("G7 CAVEAT (fires): leg-A GD/non-GD split exceeds 1 sigma -- "
          "no headline H0_A without this line")
else:
    letter = "M-POWER-DEAD" if not powA else "M-SPLIT"
P("")
P(f"verdict letter: {letter}")
P("credences: NO cell moves on any outcome (pre-registered); "
  "anomaly-real 53 / mech 8 untouched")
P(f"total wall-clock: {(time.time()-t00)/60:.1f} min")

with open('data/stage10s_skyread.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage10s_skyread.txt")
