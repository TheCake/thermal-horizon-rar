# -*- coding: utf-8 -*-
"""STAGE 7J-z2b: the coherence-kernel instrument with EXACT integration
(v2b; pre-registered before execution; SAME model and bars as v2).

WHY v2b EXISTS (ordering disclosed): v2's own pre-registered gates
killed its integration mid-run, before any gate verdict printed:
  GZQ FAILED at d = 18.7 (bar 1.0): with the kernel fitted, twin pairs'
  effective per-star variance collapses (private share + tiny sn), the
  z-integrand becomes a spike narrower than the Gauss-Hermite node
  spacing, and the quadrature carries O(10-20) lnL of error — v1's
  quadrature is retro-flagged approximate by the same finding (its
  estimator was injection-validated and self-consistent; its lnL scale
  carries an uncertainty now disclosed on its ledger row).
  Start gap 43.1 (bar 2.0): a short-lambda local basin traps one start.
The run was ABORTED (the GB0w precedent: never complete a doomed run).
v2's interim MAP for the record only: lnL 2847.9 (approx), f = 0.162,
lam = 2.95 — nothing read beyond the monitor lines, nothing shipped.

THE FIX IS EXACT, NOT FINER: every companion state's (z, gray)
integral is a pure Gaussian marginal, so the likelihood is a
closed-form bivariate normal with q-sums — no nodes, no quadrature
question, and ~10x faster:
  Sigma = [[R1^2+sg^2+V1, R1 R2+sg^2], [R1 R2+sg^2, R2^2+sg^2+V2]],
  R_i = r(c_i) sqrt(k), V_i = sn(c_i)^2 + s_plx,i^2 + r(c_i)^2 (1-k),
  k = exp(-dcol/lambda);
  L = (1-f)^2 P00 + f(1-f)(P10+P01) + f^2 P11, the P's BVN densities
  with companion shifts summed over the q-law.
State-11's q x q sum is thinned 31 -> 21 nodes per axis; gate GZS11
(replaces GZQ): |lnL(MAP; 21x21) - lnL(MAP; 31x31)| <= 0.5.

Everything else IS v2: model parameters (13), bounds, data, window,
amplitude table, bars, and the decision rule — all gates pass -> the
envelope ships as THE LANDED ANCHOR (stage7jz_prior.npz, v2b-exact);
any fail -> LIT16 stands and no further instrument today. Start
dispersion: FIVE starts (f0, lam0) = (0.08,0.5), (0.16,1.0),
(0.16,3.0), (0.25,2.0), (0.10,4.0); GZ3 requires the top two within
2.0 lnL. GZ3's nesting floor L(v2b) >= L(v1) = 2497.9 keeps v1's
(approximate-quadrature) value as the floor, disclosed.
Output: data/stage7jz2b_exact.txt (+ prior npz per the rule)
"""
import time

import numpy as np
from scipy.optimize import minimize
from astropy.io import fits

t00 = time.time()
L = []
def P(s):
    print(s, flush=True)
    L.append(s)
    with open('data/stage7jz2b_exact.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

# --- data (identical to v1/v2/7J-y) ---------------------------------------
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
P(f"window: {NP_} pairs; windowed flag {pf_win_data:.3f}; slices "
  f"({rho_lo_data:+.3f}, {rho_hi_data:+.3f})")

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
IDX11 = np.unique(np.linspace(0, NQ-1, 21).astype(int))

# --- exact likelihood -----------------------------------------------------
KNOTS = np.percentile(np.concatenate([c1w, c2w]), [5, 35, 65, 95])
V1P, V2P = s1w**2, s2w**2

def unpack(x):
    f, mu1, mu2 = x[0], x[1], x[2]
    rk, snk = x[3:7], x[7:11]
    sg, lam = x[11], x[12]
    r1 = np.interp(c1w, KNOTS, rk); r2 = np.interp(c2w, KNOTS, rk)
    sn1 = np.interp(c1w, KNOTS, snk); sn2 = np.interp(c2w, KNOTS, snk)
    k = np.exp(-dcw/lam)
    R1, R2 = r1*np.sqrt(k), r2*np.sqrt(k)
    V1 = sn1**2 + V1P + r1**2*(1.0-k)
    V2 = sn2**2 + V2P + r2**2*(1.0-k)
    S11 = R1*R1 + sg*sg + V1
    S22 = R2*R2 + sg*sg + V2
    S12 = R1*R2 + sg*sg
    det = np.maximum(S11*S22 - S12*S12, 1e-12)
    return f, mu1, mu2, S22/det, S11/det, S12/det, det

def lnL(x, wq, idx11=IDX11):
    if not (0.0 <= x[0] <= 0.6): return 1e15
    if np.any(np.abs(x[1:3]) > 0.3): return 1e15
    if np.any(np.abs(x[3:7]) > 0.8): return 1e15
    if np.any(x[7:11] < 0.02) or np.any(x[7:11] > 0.5): return 1e15
    if not (0.0 <= x[11] <= 0.3): return 1e15
    if not (0.05 <= x[12] <= 5.0): return 1e15
    f, mu1, mu2, a11, a22, a12, det = unpack(x)
    nrm = 1.0/(2*np.pi*np.sqrt(det))
    x1 = d1w - mu1
    x2 = d2w - mu2
    P00 = nrm*np.exp(-0.5*(a11*x1*x1 - 2*a12*x1*x2 + a22*x2*x2))
    X1 = x1[:, None] - AB1                       # (N, nq)
    Q10 = (a11[:, None]*X1*X1 - 2*a12[:, None]*X1*x2[:, None]
           + (a22*x2*x2)[:, None])
    P10 = nrm*np.einsum('q,nq->n', wq, np.exp(-0.5*Q10))
    X2 = x2[:, None] - AB2
    Q01 = ((a11*x1*x1)[:, None] - 2*a12[:, None]*x1[:, None]*X2
           + a22[:, None]*X2*X2)
    P01 = nrm*np.einsum('q,nq->n', wq, np.exp(-0.5*Q01))
    w11 = wq[idx11]/wq[idx11].sum()
    A1, A2 = X1[:, idx11], X2[:, idx11]
    Q11 = (a11[:, None, None]*(A1*A1)[:, :, None]
           - 2*a12[:, None, None]*A1[:, :, None]*A2[:, None, :]
           + a22[:, None, None]*(A2*A2)[:, None, :])
    P11 = nrm*np.einsum('p,q,npq->n', w11, w11, np.exp(-0.5*Q11))
    Lp = ((1-f)**2*P00 + f*(1-f)*(P10+P01) + f*f*P11)
    return -float(np.sum(np.log(np.maximum(Lp, 1e-300))))

def fit(wq, x0, maxfev=2200):
    res = minimize(lnL, x0, args=(wq,), method='Nelder-Mead',
                   options=dict(maxfev=maxfev, xatol=1e-4, fatol=1e-3,
                                adaptive=True))
    return res.x, -res.fun

X0_BASE = np.array([0.166, 0.014, 0.038, 0.165, 0.11, 0.188, 0.289,
                    0.17, 0.037, 0.059, 0.233, 0.058, 1.0])
t0 = time.time()
starts = []
for f0, lam0 in ((0.08, 0.5), (0.16, 1.0), (0.16, 3.0), (0.25, 2.0),
                 (0.10, 4.0)):
    x0 = X0_BASE.copy(); x0[0] = f0; x0[12] = lam0
    xh, lh = fit(WQ_FLAT, x0)
    starts.append((lh, xh))
    P(f"  start f0={f0} lam0={lam0}: lnL={lh:.1f}, f={xh[0]:.3f}, "
      f"lam={xh[12]:.2f} ({(time.time()-t0)/60:.1f} min)")
starts.sort(key=lambda t: -t[0])
gap01 = starts[0][0] - starts[1][0]
X_V2, L_V2 = fit(WQ_FLAT, starts[0][1], maxfev=1500)
P(f"v2b MAP (flat-q, EXACT): lnL={L_V2:.1f}, f={X_V2[0]:.4f}, "
  f"mu=({X_V2[1]:+.3f},{X_V2[2]:+.3f}), r={np.round(X_V2[3:7],3).tolist()}, "
  f"sn={np.round(X_V2[7:11],3).tolist()}, sg={X_V2[11]:.3f}, "
  f"lam={X_V2[12]:.2f}; top-two start gap {gap01:.1f}")
d11 = abs(lnL(X_V2, WQ_FLAT, np.arange(NQ)) - lnL(X_V2, WQ_FLAT))
gzs11 = d11 <= 0.5
P(f"GZS11 state-11 thinning (21^2 vs 31^2 at MAP): d={d11:.3f} -> "
  f"{'PASS' if gzs11 else 'FAIL'}")
L_V1_MAP = 2497.9
x0f = X_V2.copy()
def lnL_f0(x, wq):
    return lnL(np.concatenate([[0.0], x]), wq)
res0 = minimize(lnL_f0, x0f[1:], args=(WQ_FLAT,), method='Nelder-Mead',
                options=dict(maxfev=1800, adaptive=True))
L_F0 = -res0.fun
gz3 = (L_V2 >= L_V1_MAP - 0.5) and (L_V2 >= L_F0 - 0.5) and (gap01 <= 2.0) \
      and (0.001 < X_V2[0] < 0.58) and (0.055 < X_V2[12] < 4.9)
P(f"GZ3: L(v2b)={L_V2:.1f} vs v1 floor {L_V1_MAP} "
  f"(d={L_V2-L_V1_MAP:+.1f}, floor approx-quadrature, disclosed); "
  f"L(f=0)={L_F0:.1f} (companions {L_V2-L_F0:+.1f}); gap {gap01:.1f}; "
  f"lam interior -> {'PASS' if gz3 else 'FAIL'}")

# --- injections -----------------------------------------------------------
def simulate(x, f_true, wq, seed=1):
    rg = np.random.default_rng(seed)
    f_, mu1, mu2 = x[0], x[1], x[2]
    rk, snk = x[3:7], x[7:11]
    sg, lam = x[11], x[12]
    r1 = np.interp(c1w, KNOTS, rk); r2 = np.interp(c2w, KNOTS, rk)
    sn1 = np.interp(c1w, KNOTS, snk); sn2 = np.interp(c2w, KNOTS, snk)
    k = np.exp(-dcw/lam)
    z = rg.normal(size=NP_); g = rg.normal(size=NP_)
    out = []
    for mu, r, sn, sp, ABr in ((mu1, r1, sn1, s1w, AB1),
                               (mu2, r2, sn2, s2w, AB2)):
        eta = rg.normal(size=NP_)
        a = np.zeros(NP_)
        hasc = rg.random(NP_) < f_true
        qi = rg.choice(NQ, NP_, p=wq)
        a[hasc] = ABr[hasc, qi[hasc]]
        out.append(mu + r*(np.sqrt(k)*z + np.sqrt(1-k)*eta) + sg*g + a
                   + rg.normal(size=NP_)*np.sqrt(sn**2 + sp**2))
    return out

def fit_on(sim1, sim2, x0, maxfev=1800, pin_lam=None):
    global d1w, d2w
    keep1, keep2 = d1w, d2w
    d1w, d2w = sim1, sim2
    try:
        if pin_lam is None:
            xh, lh = fit(WQ_FLAT, x0, maxfev=maxfev)
        else:
            def lf(x, wq=WQ_FLAT):
                return lnL(np.concatenate([x, [pin_lam]]), wq)
            r = minimize(lf, x0[:12], method='Nelder-Mead',
                         options=dict(maxfev=maxfev, adaptive=True))
            xh, lh = np.concatenate([r.x, [pin_lam]]), -r.fun
    finally:
        d1w, d2w = keep1, keep2
    return xh, lh

gz_ok = {}
sm1, sm2 = simulate(X_V2, 0.0, WQ_FLAT, seed=911)
x0 = X_V2.copy(); x0[0] = 0.05
xh, _ = fit_on(sm1, sm2, x0)
gz_ok['GZ1c-a'] = xh[0] <= 0.03
P(f"GZ1c-a kernel-truth f=0, v2b fit: f_hat={xh[0]:.3f} -> "
  f"{'PASS' if gz_ok['GZ1c-a'] else 'FAIL'} ({(time.time()-t0)/60:.1f} min)")
xh1, _ = fit_on(sm1, sm2, x0, pin_lam=4.99)
P(f"GZ1c-b DIAGNOSTIC same sky, v1-style fit (lam pinned): "
  f"f_hat={xh1[0]:.3f} — the kernel-blind over-attribution")
for tag, ftr in (('GZ2a', 0.10), ('GZ2b', 0.25)):
    sm1, sm2 = simulate(X_V2, ftr, WQ_FLAT, seed=int(2000+ftr*100))
    x0 = X_V2.copy(); x0[0] = ftr
    xh, _ = fit_on(sm1, sm2, x0)
    ok_ = abs(xh[0]-ftr) <= max(0.03, 0.25*ftr)
    gz_ok[tag] = ok_
    P(f"{tag} f={ftr}: recovered {xh[0]:.3f} -> "
      f"{'PASS' if ok_ else 'FAIL'} ({(time.time()-t0)/60:.1f} min)")

sm1, sm2 = simulate(X_V2, X_V2[0], WQ_FLAT, seed=424243)
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
  f"{X_V2[0]:.3f}/{xa[0]:.3f}" + ("  ** GRID EDGE **" if edge_flag else ""))

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
P(f"HOST-axis (v2b): f_host per component in [{fh1.min():.2f}, "
  f"{fh1.max():.2f}] (peak {fh_hat:.2f})")
for refv, lab in ((0.09, 'literature (H)'), (0.10, 'old fence'),
                  (0.20, 'likelihood optimum'), (0.35, 'posterior cell'),
                  (0.51, 'retracted part-A peak')):
    iv = int(np.argmin(np.abs(FH_GRID-refv)))
    P(f"  ln pi_host({refv:.2f} = {lab}) = {env_fh[iv]-env_fh.max():.1f}")

gates_all = gz0 and gzs11 and gz3 and all(gz_ok.values()) \
            and g6i and g6ii and g6iii and not edge_flag
P("")
P(f"f_hat chain: v1 0.166 -> v2 (approx) 0.162 -> v2b (exact) "
  f"{X_V2[0]:.3f}")
if gates_all:
    np.savez('data/stage7jz_prior.npz', fh_grid=FH_GRID,
             lnpi_host=env_fh, f_prof_grid=F_PROF, lnpi_flat=lnpi_flat,
             lnpi_alt=lnpi_alt, f_map_flat=X_V2[0], f_map_alt=xa[0],
             fh_hat=fh_hat, version='v2b-exact', x_map=X_V2, knots=KNOTS)
    P("GATES ALL PASS -> the v2b LANDED ANCHOR SHIPS "
      "(data/stage7jz_prior.npz, version v2b-exact)")
else:
    np.savez('data/stage7jz_prior_v2b_unshipped.npz', fh_grid=FH_GRID,
             lnpi_host=env_fh, x_map=X_V2, version='v2b-unshipped')
    P("GATE FAIL -> nothing ships; the deciders read LIT16 "
      "(pre-registered fallback); quarantined npz written")
P(f"[total {(time.time()-t00)/60:.1f} min]")
