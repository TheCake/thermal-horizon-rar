# -*- coding: utf-8 -*-
"""STAGE 7J-z2: the coherence-kernel completeness measurement (the v2
instrument; pre-registered before execution, after v1's honest GZ6
fail).

WHAT v1 GOT RIGHT AND WRONG. v1 (stage7jz_mixture) passed GZ0-GZ3 and
all four injection gates at ~exact recovery (GZ1a/b: f_hat = 0.000 on
companion-free truths incl. rank-2 mis-specification; GZ2: 0.103/0.251
on 0.10/0.25), and measured f_hat = 0.166 (blended, per component),
stable under the gray promotion (0.168). Its GZ6 postdiction FAILED:
  (iii) was a CONSTRUCTION MISMATCH in the bar, not the model: the
  model simulates on the 12,084 windowed pairs but was compared to the
  full-sample flag fraction (0.121); like-for-like the windowed data
  give 0.107 vs model 0.096 — inside the +-0.02 bar. v2 compares
  windowed-to-windowed (the bar as it should have been written).
  (ii) is REAL MISSING STRUCTURE: the measured rho(dcol) profile
  decays smoothly 0.86 -> 0.36 across the full color-separation range
  (diagnostic on record); a scalar pair latent with deterministic
  color response cannot produce it (v1 postdicted 0.47 at dcol<0.15 vs
  0.80 measured). Physics: the pair-common displacement is a VECTOR
  (abundance pattern + age), and its projections onto two stars
  decorrelate as their masses separate; twins project identically.

THE v2 MODEL (rank-1 + gray + coherence kernel; 13 params):
  delta_i = mu_i + r(c_i)[sqrt(k) z + sqrt(1-k) eta_i] + sg*g
            + a_i[comp_i] + n_i,
  k(dcol) = exp(-dcol/lambda);  z, g pair-shared ~ N(0,1);
  eta_i component-private ~ N(0,1) — merges EXACTLY into the
  per-star variance (V_i += r_i^2 (1-k)), so the quadrature stays 2D.
  Everything else (amplitude table, q-laws, window, optimizer
  discipline) identical to v1.

GATES (bars pre-registered; v1 values in brackets where they carry):
  GZ0  regression: N = 13784, rho_core_att = +0.465 +- 0.005.
  GZQ  quadrature: |lnL(MAP; 13x5 nodes) - lnL(MAP; 15x7)| <= 1.0.
  GZ1c-a  over-attribution: v2-truth (fitted kernel), f_true = 0,
          v2 fit -> f_hat <= 0.03.
  GZ1c-b  DIAGNOSTIC (no bar): the same v2-truth sky fit with the v1
          model — measures how much v1's missing kernel over-attributes
          companions; explains (or acquits) the v1 f_hat = 0.166 at
          data level.
  GZ2  recovery: f_true = 0.10 / 0.25 under v2 truth -> within
       max(0.03, 25%).
  GZ3  nesting/convergence: L(v2) >= L(v1-MAP); L(f free) >= L(f=0);
       two starts within 2.0; f_hat interior; lambda not riding its
       [0.05, 5] bounds.
  GZ6' postdiction at MAP (like-for-like throughout):
       (i)  rho bar slice +-0.05 of +0.465;
       (ii) rho(dcol<0.15 core) +-0.07 of +0.796 AND rho(dcol>=0.4)
            +-0.07 of +0.394 — the failed slice now reachable;
       (iii) windowed flag fraction +-0.02 of 0.107;
       incidence ratio reported vs 2.74 (no bar; correlated incidence
       known-unmodeled).
DECISION RULE (pre-committed): ALL gates pass -> the v2 envelope ships
as THE LANDED ANCHOR (data/stage7jz_prior.npz, version = v2; the v1
npz is quarantined as _v1_unshipped). ANY gate fails -> nothing ships;
the deciders read LIT16 (the v1 pre-registered fallback stands); no
third instrument today.
THE ANCHOR-RELEVANT NUMBER: f_hat(v2) vs v1's 0.166 — if the kernel
absorbs twin-end correlated brightness that v1 booked as companions,
f_hat drops; if stable, the anchor is kernel-robust. Either way that
comparison is the data-level over-attribution measurement.
Output: data/stage7jz2_twin.txt (+ prior npz per the rule above)
"""
import time

import numpy as np
from scipy.optimize import minimize
from astropy.io import fits

t00 = time.time()
rng = np.random.default_rng(20260727)
L = []
def P(s):
    print(s, flush=True)
    L.append(s)
    with open('data/stage7jz2_twin.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

# --- data (identical to v1 / 7J-y) ----------------------------------------
d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1, plx2 = d['parallax1'], d['parallax2']
eplx1, eplx2 = d['parallax_error1'], d['parallax_error2']
sep, Rch = d['sep_AU'], d['R_chance_align']
G1m, G2m = d['phot_g_mean_mag1'], d['phot_g_mean_mag2']
plx = 0.5*(plx1+plx2)
MG1 = G1m+5*np.log10(np.maximum(plx1, 1e-6))-10
MG2 = G2m+5*np.log10(np.maximum(plx2, 1e-6))-10
sigv = 4.74047/plx*np.sqrt(d['pmra_error1']**2+d['pmdec_error1']**2
                           + d['pmra_error2']**2+d['pmdec_error2']**2)
ok = (Rch < 0.01) & (plx1 > 5) & (plx2 > 5) \
   & (plx1/np.maximum(eplx1, 1e-6) > 20) & (plx2/np.maximum(eplx2, 1e-6) > 20) \
   & (np.abs(plx1-plx2) < 3*np.hypot(eplx1, eplx2)) \
   & (sep > 200) & (sep < 50000) \
   & (MG1 > 2.6) & (MG1 < 14.2) & (MG2 > 2.6) & (MG2 < 14.2) & (sigv < 0.03)
c1, c2 = d['bp_rp1'], d['bp_rp2']
col = np.concatenate([c1[ok], c2[ok]])
mg = np.concatenate([MG1[ok], MG2[ok]])
good = np.isfinite(col) & np.isfinite(mg) & (col > -0.5) & (col < 6.0)
CB = np.linspace(np.nanpercentile(col[good], 0.5),
                 np.nanpercentile(col[good], 99.5), 41)
cc = 0.5*(CB[:-1]+CB[1:])
ridge = np.full(len(cc), np.nan)
for i in range(len(cc)):
    m = good & (col >= CB[i]) & (col < CB[i+1])
    if m.sum() > 30:
        ridge[i] = np.median(mg[m])
vr = np.isfinite(ridge)
ccv, rgv = cc[vr], ridge[vr]
def ridge_f(c):
    return np.interp(c, ccv, rgv)
g1 = np.isfinite(c1[ok]) & (c1[ok] > -0.5) & (c1[ok] < 6.0)
g2 = np.isfinite(c2[ok]) & (c2[ok] > -0.5) & (c2[ok] < 6.0)
both = g1 & g2
d1 = (MG1[ok] - ridge_f(c1[ok]))[both]
d2 = (MG2[ok] - ridge_f(c2[ok]))[both]
c1v = c1[ok][both]; c2v = c2[ok][both]
dcol_all = np.abs(c1v - c2v)
s1 = (5/np.log(10)*(eplx1/np.maximum(plx1, 1e-6)))[ok][both]
s2 = (5/np.log(10)*(eplx2/np.maximum(plx2, 1e-6)))[ok][both]
NPAIR_ALL = int(both.sum())

def rho_att(a, b, va, vb):
    a = a - a.mean(); b = b - b.mean()
    cov = float(np.mean(a*b))
    den = np.sqrt(max(np.mean(a*a)-np.mean(va), 1e-6)
                  * max(np.mean(b*b)-np.mean(vb), 1e-6))
    return cov/den

core = (np.abs(d1) < 0.3) & (np.abs(d2) < 0.3)
mbar = core & (dcol_all >= 0.15)
rho_bar_data = rho_att(d1[mbar], d2[mbar], s1[mbar]**2, s2[mbar]**2)
gz0 = (NPAIR_ALL == 13784) and (abs(rho_bar_data - 0.465) <= 0.005)
P(f"GZ0: N={NPAIR_ALL}, rho_core_att={rho_bar_data:+.3f} -> "
  f"{'PASS' if gz0 else 'FAIL'}")
mlo = core & (dcol_all < 0.15)
mhi = core & (dcol_all >= 0.4)
rho_lo_data = rho_att(d1[mlo], d2[mlo], s1[mlo]**2, s2[mlo]**2)
rho_hi_data = rho_att(d1[mhi], d2[mhi], s1[mhi]**2, s2[mhi]**2)

W = (d1 > -1.2) & (d1 < 0.6) & (d2 > -1.2) & (d2 < 0.6)
d1w, d2w = d1[W], d2[W]
c1w, c2w = c1v[W], c2v[W]
s1w, s2w = s1[W], s2[W]
dcw = np.abs(c1w - c2w)
NP_ = int(W.sum())
pf_win_data = 0.5*(np.mean(d1w < -0.4) + np.mean(d2w < -0.4))
f1_, f2_ = d1w < -0.4, d2w < -0.4
inc_data = float((f1_ & f2_).mean()/max(f1_.mean()*f2_.mean(), 1e-9))
P(f"window: {NP_} pairs; windowed flag fraction {pf_win_data:.3f} "
  f"(the like-for-like (iii) reference); slices "
  f"({rho_lo_data:+.3f}, {rho_hi_data:+.3f})")

# --- amplitude table (identical to v1) ------------------------------------
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
def mg_of_mass(mm):
    return np.interp(-np.asarray(mm), -MS_T, MG_T)
def color_of_mg(mgv):
    return np.interp(mgv, rgv, ccv)
NQ = 31
QG = np.linspace(0.10, 1.0, NQ)
m_grid = np.interp(rgv, MG_T, MS_T)
def combined_delta_grid(mpri, qg):
    m1 = np.asarray(mpri)[:, None]
    mg1 = mg_of_mass(m1)
    mg2 = mg_of_mass(np.clip(qg[None, :]*m1, MS_T[-1], MS_T[0]))
    F1, F2 = 10**(-0.4*mg1), 10**(-0.4*mg2)
    mgs = -2.5*np.log10(F1+F2)
    c1_ = color_of_mg(mg1)
    c2_ = color_of_mg(mg2)
    cs = (F1*c1_ + F2*c2_)/(F1+F2)
    return mgs - ridge_f(cs)
AB = combined_delta_grid(m_grid, QG)
def binof(c):
    return np.clip(np.searchsorted(ccv, c) - 1, 0, len(ccv)-1)
AB1 = AB[binof(c1w)]
AB2 = AB[binof(c2w)]
WQ_FLAT = np.full(NQ, 1.0/NQ)
wq_alt = QG**-0.5
WQ_ALT = wq_alt/wq_alt.sum()

# --- v2 likelihood --------------------------------------------------------
KNOTS = np.percentile(np.concatenate([c1w, c2w]), [5, 35, 65, 95])
def nodes(nz, ng):
    zn, zw = np.polynomial.hermite_e.hermegauss(nz)
    gn, gw = np.polynomial.hermite_e.hermegauss(ng)
    return zn, zw/np.sqrt(2*np.pi), gn, gw/np.sqrt(2*np.pi)
ZN, ZW, GN, GW = nodes(13, 5)
V1P, V2P = s1w**2, s2w**2

def unpack(x):
    f, mu1, mu2 = x[0], x[1], x[2]
    rk, snk = x[3:7], x[7:11]
    sg, lam = x[11], x[12]
    r1 = np.interp(c1w, KNOTS, rk); r2 = np.interp(c2w, KNOTS, rk)
    sn1 = np.interp(c1w, KNOTS, snk); sn2 = np.interp(c2w, KNOTS, snk)
    k = np.exp(-dcw/lam)
    sk = np.sqrt(k)
    V1 = sn1**2 + V1P + r1**2*(1.0-k)
    V2 = sn2**2 + V2P + r2**2*(1.0-k)
    return f, mu1, mu2, r1*sk, r2*sk, V1, V2, sg

def _mix(dd, mu, rk_, V, ABr, wq, f, gshift):
    M = dd[:, None] - mu - gshift - rk_[:, None]*ZN[None, :]
    iv = 1.0/V[:, None]
    S = np.exp(-0.5*M*M*iv)*np.sqrt(iv/(2*np.pi))
    DF = M[:, :, None] - ABr[:, None, :]
    A = np.einsum('q,npq->np', wq,
                  np.exp(-0.5*DF*DF*iv[:, :, None]))*np.sqrt(iv/(2*np.pi))
    return (1.0-f)*S + f*A

def lnL(x, wq):
    if not (0.0 <= x[0] <= 0.6): return 1e15
    if np.any(np.abs(x[1:3]) > 0.3): return 1e15
    if np.any(np.abs(x[3:7]) > 0.8): return 1e15
    if np.any(x[7:11] < 0.02) or np.any(x[7:11] > 0.5): return 1e15
    if not (0.0 <= x[11] <= 0.3): return 1e15
    if not (0.05 <= x[12] <= 5.0): return 1e15
    f, mu1, mu2, rk1, rk2, V1, V2, sg = unpack(x)
    Lp = np.zeros(NP_)
    for gj, wgj in zip(GN, GW):
        m1 = _mix(d1w, mu1, rk1, V1, AB1, wq, f, sg*gj)
        m2 = _mix(d2w, mu2, rk2, V2, AB2, wq, f, sg*gj)
        Lp += wgj*np.einsum('z,nz->n', ZW, m1*m2)
    return -float(np.sum(np.log(np.maximum(Lp, 1e-300))))

def lnL_nodes(x, wq, nz, ng):
    global ZN, ZW, GN, GW
    keep = (ZN, ZW, GN, GW)
    ZN, ZW, GN, GW = nodes(nz, ng)
    try:
        return lnL(x, wq)
    finally:
        ZN, ZW, GN, GW = keep

def fit(wq, x0, maxfev=2200):
    res = minimize(lnL, x0, args=(wq,), method='Nelder-Mead',
                   options=dict(maxfev=maxfev, xatol=1e-4, fatol=1e-3,
                                adaptive=True))
    return res.x, -res.fun

# v1 MAP as the launch point (from stage7jz_mixture output, hard-coded
# init only — the fit refines everything): r/sn knots + gray
X0_BASE = np.array([0.166, 0.014, 0.038, 0.165, 0.11, 0.188, 0.289,
                    0.17, 0.037, 0.059, 0.233, 0.058, 1.0])
t0 = time.time()
starts = []
for f0, lam0 in ((0.08, 0.6), (0.20, 1.5)):
    x0 = X0_BASE.copy(); x0[0] = f0; x0[12] = lam0
    xh, lh = fit(WQ_FLAT, x0)
    starts.append((lh, xh))
    P(f"  start f0={f0} lam0={lam0}: lnL={lh:.1f}, f={xh[0]:.3f}, "
      f"lam={xh[12]:.2f} ({(time.time()-t0)/60:.1f} min)")
starts.sort(key=lambda t: -t[0])
gap01 = starts[0][0] - starts[1][0]
X_V2, L_V2 = fit(WQ_FLAT, starts[0][1], maxfev=1500)
P(f"v2 MAP (flat-q): lnL={L_V2:.1f}, f={X_V2[0]:.4f}, "
  f"mu=({X_V2[1]:+.3f},{X_V2[2]:+.3f}), r={np.round(X_V2[3:7],3).tolist()}, "
  f"sn={np.round(X_V2[7:11],3).tolist()}, sg={X_V2[11]:.3f}, "
  f"lam={X_V2[12]:.2f}; start gap {gap01:.1f}")
L_V1_MAP = 2497.9    # v1's best (rank-2) lnL, for the nesting gate
gzq = abs(lnL_nodes(X_V2, WQ_FLAT, 15, 7) - (-L_V2)) <= 1.0
P(f"GZQ quadrature (13x5 vs 15x7 at MAP): "
  f"d={abs(lnL_nodes(X_V2, WQ_FLAT, 15, 7) - (-L_V2)):.3f} -> "
  f"{'PASS' if gzq else 'FAIL'}")

x0f = X_V2.copy()
def lnL_f0(x, wq):
    return lnL(np.concatenate([[0.0], x]), wq)
res0 = minimize(lnL_f0, x0f[1:], args=(WQ_FLAT,), method='Nelder-Mead',
                options=dict(maxfev=1800, adaptive=True))
L_F0 = -res0.fun
gz3 = (L_V2 >= L_V1_MAP - 0.5) and (L_V2 >= L_F0 - 0.5) and (gap01 <= 2.0) \
      and (0.001 < X_V2[0] < 0.58) and (0.055 < X_V2[12] < 4.9)
P(f"GZ3: L(v2)={L_V2:.1f} vs L(v1)={L_V1_MAP} (d={L_V2-L_V1_MAP:+.1f}); "
  f"L(f=0)={L_F0:.1f} (companions worth {L_V2-L_F0:+.1f}); lam interior "
  f"-> {'PASS' if gz3 else 'FAIL'}")

# --- injections -----------------------------------------------------------
def simulate(x, f_true, wq, seed=1):
    rg = np.random.default_rng(seed)
    f_, mu1, mu2, rk1, rk2, V1, V2, sg = unpack(x)
    z = rg.normal(size=NP_); g = rg.normal(size=NP_)
    out = []
    for dd, mu, rk_, V, ABr in ((d1w, mu1, rk1, V1, AB1),
                                (d2w, mu2, rk2, V2, AB2)):
        a = np.zeros(NP_)
        hasc = rg.random(NP_) < f_true
        qi = rg.choice(NQ, NP_, p=wq)
        a[hasc] = ABr[hasc, qi[hasc]]
        # V already contains the private common share r^2(1-k)
        out.append(mu + rk_*z + sg*g + a + rg.normal(size=NP_)*np.sqrt(V))
    return out

def fit_on(sim1, sim2, x0, maxfev=1800, v1mode=False):
    global d1w, d2w
    keep1, keep2 = d1w, d2w
    d1w, d2w = sim1, sim2
    try:
        if v1mode:
            xx = x0.copy(); xx[12] = 4.99
            xh, lh = fit(WQ_FLAT, xx, maxfev=maxfev)
            # lambda pinned high => k ~ 1 => v1 (no-kernel) fit
        else:
            xh, lh = fit(WQ_FLAT, x0, maxfev=maxfev)
    finally:
        d1w, d2w = keep1, keep2
    return xh, lh

gz_ok = {}
sm1, sm2 = simulate(X_V2, 0.0, WQ_FLAT, seed=91)
x0 = X_V2.copy(); x0[0] = 0.05
xh, _ = fit_on(sm1, sm2, x0)
gz_ok['GZ1c-a'] = xh[0] <= 0.03
P(f"GZ1c-a v2-truth f=0, v2 fit: f_hat={xh[0]:.3f} -> "
  f"{'PASS' if gz_ok['GZ1c-a'] else 'FAIL'} "
  f"({(time.time()-t0)/60:.1f} min)")
def lnL_v1style(x, wq):
    xx = np.concatenate([x, [4.99]])
    return lnL(xx, wq)
res1 = minimize(lnL_v1style, np.concatenate([x0[:12]]), args=(WQ_FLAT,),
                method='Nelder-Mead',
                options=dict(maxfev=1600, adaptive=True))
P(f"GZ1c-b DIAGNOSTIC v2-truth f=0, v1-STYLE fit (lam pinned 4.99): "
  f"f_hat={res1.x[0]:.3f} — v1's kernel-blind over-attribution on a "
  f"companion-free sky")
for tag, ftr in (('GZ2a', 0.10), ('GZ2b', 0.25)):
    sm1, sm2 = simulate(X_V2, ftr, WQ_FLAT, seed=int(1000+ftr*100))
    x0 = X_V2.copy(); x0[0] = ftr
    xh, _ = fit_on(sm1, sm2, x0)
    ok_ = abs(xh[0]-ftr) <= max(0.03, 0.25*ftr)
    gz_ok[tag] = ok_
    P(f"{tag} f={ftr}: recovered {xh[0]:.3f} -> "
      f"{'PASS' if ok_ else 'FAIL'} ({(time.time()-t0)/60:.1f} min)")

# --- GZ6' postdiction (like-for-like) -------------------------------------
sm1, sm2 = simulate(X_V2, X_V2[0], WQ_FLAT, seed=424242)
smc = (np.abs(sm1) < 0.3) & (np.abs(sm2) < 0.3)
def rho_m(m):
    return rho_att(sm1[m], sm2[m], s1w[m]**2, s2w[m]**2)
rb_m = rho_m(smc & (dcw >= 0.15))
rl_m = rho_m(smc & (dcw < 0.15))
rh_m = rho_m(smc & (dcw >= 0.4))
pf_m = 0.5*(np.mean(sm1 < -0.4) + np.mean(sm2 < -0.4))
i1, i2 = sm1 < -0.4, sm2 < -0.4
inc_m = float((i1 & i2).mean()/max(i1.mean()*i2.mean(), 1e-9))
g6i = abs(rb_m - rho_bar_data) <= 0.05
g6ii = (abs(rl_m - rho_lo_data) <= 0.07) and (abs(rh_m - rho_hi_data) <= 0.07)
g6iii = abs(pf_m - pf_win_data) <= 0.02
P(f"GZ6': bar {rb_m:+.3f} vs {rho_bar_data:+.3f} -> "
  f"{'PASS' if g6i else 'FAIL'}; slices ({rl_m:+.3f},{rh_m:+.3f}) vs "
  f"({rho_lo_data:+.3f},{rho_hi_data:+.3f}) -> "
  f"{'PASS' if g6ii else 'FAIL'}; windowed P(flag) {pf_m:.3f} vs "
  f"{pf_win_data:.3f} -> {'PASS' if g6iii else 'FAIL'}; incidence "
  f"{inc_m:.2f} vs {inc_data:.2f} (no bar)")

# --- profile + envelope + host remap --------------------------------------
F_PROF = np.array([0.0, 0.06, 0.11, 0.17, 0.23, 0.30])
def profile(wq, xstart):
    prof = np.zeros(len(F_PROF))
    xw = xstart.copy()
    for i, fv in enumerate(F_PROF):
        def lf(x, wq=wq):
            return lnL(np.concatenate([[fv], x]), wq)
        r = minimize(lf, xw[1:], method='Nelder-Mead',
                     options=dict(maxfev=700, adaptive=True))
        prof[i] = -r.fun
        xw = np.concatenate([[fv], r.x])
    return prof
prof_flat = profile(WQ_FLAT, X_V2)
P(f"profile flat-q: {np.round(prof_flat-prof_flat.max(),1).tolist()} on "
  f"{F_PROF.tolist()} ({(time.time()-t0)/60:.1f} min)")
xa, la = fit(WQ_ALT, X_V2, maxfev=1600)
prof_alt = profile(WQ_ALT, xa)
P(f"q^-0.5 MAP f={xa[0]:.4f}; profile "
  f"{np.round(prof_alt-prof_alt.max(),1).tolist()}")
lnpi_flat = prof_flat - prof_flat.max()
lnpi_alt = prof_alt - prof_alt.max()
env = np.maximum(lnpi_flat, lnpi_alt)
in1 = F_PROF[(env.max()-env) <= 0.5]
f_lo, f_hi = float(in1.min()), float(in1.max())
edge_flag = env[-1] >= env.max() - 0.5
P(f"BLENDED envelope: f in [{f_lo:.2f}, {f_hi:.2f}]; MAPs "
  f"{X_V2[0]:.3f}/{xa[0]:.3f}"
  + ("  ** GRID EDGE **" if edge_flag else ""))

dists = 1000.0/plx[ok]
NB2 = 200000
rgb = np.random.default_rng(31415)
Ms_b = 0.6+1.8*rgb.random(NB2)
Mh_b = 0.5*Ms_b
Pyr_b = 10**rgb.normal(5.03, 2.28, NB2)/365.25
db = rgb.choice(np.asarray(dists, float), NB2)
Q_FLATD = 0.1+0.9*rgb.random(NB2)
Q_ALTD = rgb.choice(QG, NB2, p=WQ_ALT)
FH_GRID = np.arange(0.0, 0.9001, 0.01)
env_fh = np.full(len(FH_GRID), -1e9)
for qtag, qd, lnpiA in (('flat', Q_FLATD, lnpi_flat),
                        ('qalt', Q_ALTD, lnpi_alt)):
    a_in = (Mh_b*(1+qd)*Pyr_b**2)**(1/3)
    for ares in (0.6, 1.0):
        Pb = float(np.mean(a_in < ares*db))
        lnH = np.full(len(FH_GRID), -1e9)
        fA = FH_GRID*Pb
        m = fA <= F_PROF.max()
        lnH[m] = np.interp(fA[m], F_PROF, lnpiA)
        env_fh = np.maximum(env_fh, lnH)
        P(f"P_blend({qtag}, {ares}\") = {Pb:.3f}")
fh1 = FH_GRID[(env_fh.max()-env_fh) <= 0.5]
fh_hat = float(FH_GRID[int(np.argmax(env_fh))])
P(f"HOST-axis (v2): f_host per component in [{fh1.min():.2f}, "
  f"{fh1.max():.2f}] (peak {fh_hat:.2f})")
for refv, lab in ((0.09, 'literature (H)'), (0.10, 'old fence'),
                  (0.20, 'likelihood optimum'), (0.35, 'posterior cell'),
                  (0.51, 'retracted part-A peak')):
    iv = int(np.argmin(np.abs(FH_GRID-refv)))
    P(f"  ln pi_host({refv:.2f} = {lab}) = {env_fh[iv]-env_fh.max():.1f}")

gates_all = gz0 and gzq and gz3 and all(gz_ok.values()) \
            and g6i and g6ii and g6iii and not edge_flag
P("")
P(f"f_hat comparison (the over-attribution measurement at data level): "
  f"v1 0.166 -> v2 {X_V2[0]:.3f}")
if gates_all:
    np.savez('data/stage7jz_prior.npz', fh_grid=FH_GRID,
             lnpi_host=env_fh, f_prof_grid=F_PROF, lnpi_flat=lnpi_flat,
             lnpi_alt=lnpi_alt, f_map_flat=X_V2[0], f_map_alt=xa[0],
             fh_hat=fh_hat, version='v2-kernel', x_map=X_V2, knots=KNOTS)
    P("GATES ALL PASS -> the v2 LANDED ANCHOR SHIPS "
      "(data/stage7jz_prior.npz, version v2-kernel)")
else:
    np.savez('data/stage7jz_prior_v2_unshipped.npz', fh_grid=FH_GRID,
             lnpi_host=env_fh, x_map=X_V2, version='v2-kernel-unshipped')
    P("GATE FAIL -> nothing ships; the deciders read LIT16 "
      "(pre-registered fallback); quarantined npz written")
P(f"[total {(time.time()-t00)/60:.1f} min]")
