# -*- coding: utf-8 -*-
"""STAGE 7J-z part 1: the repaired completeness measurement — a joint 2D
(delta1, delta2) mixture with PER-PAIR common-mode response slopes.
(The gating item fixed by review rounds 3-6; replaces the RETRACTED
part-A scale, correction #18.)

WHY 2D. Part A fit the per-star delta distribution: companions and
common-mode astrophysics (metallicity / age / extinction) both put
stars above the ridge, and the 1D fit attributed the whole bright wing
to companions (7J-y measured rho_core_att = +0.465 => C-FAIL). The 2D
geometry separates them WITHOUT external calibration:
  - common-mode = a PAIR-level latent z displacing BOTH components,
    with color-dependent gain r(c) (an F star and an M star respond to
    the same metallicity offset differently — the per-pair "diagonal"
    slope is r(c2)/r(c1), NOT 1; a global slope would recreate the
    over-attribution, review round 5);
  - companion = a ONE-component displacement, negative only (flux sum),
    with the amplitude law inherited unchanged from part A
    (combined magnitude AND combined color vs the same 41-bin ridge);
  - noise = per-component, symmetric, color-dependent width sn(c) plus
    the known parallax-magnitude term.

MODEL (rank-1 primary). Per pair j, latent z ~ N(0,1):
  delta_i = mu_i + r(c_i) z + a_i [comp_i] + n_i,  n_i ~ N(0, sn(c_i)^2
  + s_plx,i^2);  comp_i ~ Bernoulli(f) independent per component;
  a_i ~ the part-A combined-delta table at the component's color (q-law
  flat(0.1,1) primary, q^-0.5 envelope variant). r and sn are free
  curves on 4 color knots (r in mag units = the common-mode width
  profile; z-scale absorbed). Params: f, mu1, mu2, r x4, sn x4 (11).
  RANK-2 variant: + a second pair-common latent with CONSTANT (gray)
  response sigma_g (differential extinction / age direction) — spans
  the {r(c), 1} response plane.
  Likelihood: Gauss-Hermite over z (15 nodes; x7 gray nodes rank-2);
  given z the components factor: L = sum_z w_z prod_i [(1-f) S_i(z)
  + f A_i(z)] — all four companion states exact.
  Fit window delta in [-1.2, +0.6] both components (part A's window).
  Optimizer: Nelder-Mead, multi-start, warm polish (5B convergence
  discipline: nesting gates below).

GATES (pre-registered, committed before execution):
  GZ0 regression: loading identical to 7J-y — N pairs = 13784 exact;
      rho_core_att(|dcol|>=0.15 slice) = +0.465 within +-0.005.
  GZ1 OVER-ATTRIBUTION CONTROL (the #18 gate):
      GZ1a: inject companion-FREE sky (rank-1 common at fitted
            amplitudes + noise, f_true = 0), fit rank-1: f_hat <= 0.03.
      GZ1b: inject rank-2 truth (add gray sigma_g = 0.05, a chosen
            differential-extinction scale, flagged as chosen), fit
            rank-1: f_hat <= 0.05.
      DECISION RULE: GZ1b fail -> the rank-2 fit becomes PRIMARY and
      GZ1a/GZ1b are re-run against it (one promotion, pre-registered).
  GZ2 recovery: inject f_true = 0.10 and 0.25 (common-mode at fitted
      amplitudes): recover within max(+-0.03, 25% rel).
  GZ3 nesting/convergence: L(f free) >= L(f=0); L(rank-2) >= L(rank-1);
      best two starts within 2.0 lnL; f_hat not riding the 0.6 bound;
      profile grid edge (0.30) -> extend before quoting.
  GZ6 postdiction (model-class check, simulate at MAP on the data
      covariates): (i) rho_core_att within +-0.05 of +0.465;
      (ii) rho in the two dcol slices (<0.15: +0.796, >=0.4: +0.394)
      within +-0.07 — the dcol TREND is the rank-1 geometry's own test;
      (iii) per-star P(delta < -0.4) within +-0.02 of the observed;
      (iv) tail-incidence ratio reported vs 2.74 (NO bar: correlated
      subsystem incidence, Tokovinin, is known-unmodeled; independent-f
      underpredicts it; reported, carried).
  Any GZ0/GZ2/GZ3 fail, or GZ6 fail after the promotion branch ->
  NEEDS REFINEMENT: no prior is shipped; part 2 reads at the
  literature anchor only.

OUTPUT (the LANDED ANCHOR): profile lnL over f (7-point grid, envelope
over both q-laws), pushed through part A's HOST-fraction remap verbatim
(P_blend from the model's own (q, logP) laws at the sample's distances,
a_res in {0.6", 1.0"} carried) -> data/stage7jz_prior.npz
(fh_grid, lnpi_host envelope + the measured common-mode field for the
post-7J-z arm suite). All axes PER COMPONENT (the v7 fcomp convention;
stated to close the per-pair/per-component ambiguity in the lit-0.16
anchor, noted this stage). Comparators printed: retracted part-A
f_photo 0.22-0.30 / f_host 0.42-0.57; literature per-component ~0.10;
this is also the first subsystem-rate MEASUREMENT at 0.2-50 kAU outer
separations (novelty scout running; claim held at scout-grade).
Output: data/stage7jz_mixture.txt
"""
import sys
import time

import numpy as np
from scipy.optimize import minimize
from astropy.io import fits

t00 = time.time()
rng = np.random.default_rng(20260726)
L = []
def P(s):
    print(s, flush=True)
    L.append(s)
    with open('data/stage7jz_mixture.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

# --- data (identical to 7J-y / part A) ------------------------------------
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
dcol = np.abs(c1v - c2v)
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
mbar = core & (dcol >= 0.15)
rho_bar_data = rho_att(d1[mbar], d2[mbar], s1[mbar]**2, s2[mbar]**2)
gz0 = (NPAIR_ALL == 13784) and (abs(rho_bar_data - 0.465) <= 0.005)
P(f"GZ0 regression: N={NPAIR_ALL} (ref 13784), rho_core_att="
  f"{rho_bar_data:+.3f} (ref +0.465) -> {'PASS' if gz0 else 'FAIL'}")
mlo = core & (dcol < 0.15)
mhi = core & (dcol >= 0.4)
rho_lo_data = rho_att(d1[mlo], d2[mlo], s1[mlo]**2, s2[mlo]**2)
rho_hi_data = rho_att(d1[mhi], d2[mhi], s1[mhi]**2, s2[mhi]**2)
pflag_data = 0.5*(np.mean(d1 < -0.4) + np.mean(d2 < -0.4))
f1_, f2_ = d1 < -0.4, d2 < -0.4
inc_data = float((f1_ & f2_).mean()/max(f1_.mean()*f2_.mean(), 1e-9))

# fit window (part A's histogram window, both components)
W = (d1 > -1.2) & (d1 < 0.6) & (d2 > -1.2) & (d2 < 0.6)
d1w, d2w = d1[W], d2[W]
c1w, c2w = c1v[W], c2v[W]
s1w, s2w = s1[W], s2[W]
NP_ = int(W.sum())
P(f"fit window [-1.2, +0.6] both: {NP_} pairs ({NPAIR_ALL - NP_} dropped)")

# --- the part-A companion amplitude table, per ridge color bin ------------
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
m_grid = np.interp(rgv, MG_T, MS_T)          # mass at each ridge color
def combined_delta_grid(mpri, qg):
    m1 = np.asarray(mpri)[:, None]
    mg1 = mg_of_mass(m1)
    mg2 = mg_of_mass(np.clip(qg[None, :]*m1, MS_T[-1], MS_T[0]))
    F1, F2 = 10**(-0.4*mg1), 10**(-0.4*mg2)
    mgs = -2.5*np.log10(F1+F2)
    c1_ = color_of_mg(mg1)
    c2_ = color_of_mg(mg2)
    cs = (F1*c1_ + F2*c2_)/(F1+F2)
    return mgs - ridge_f(cs)                  # (nbin, nq), <= 0
AB = combined_delta_grid(m_grid, QG)          # amplitude table per color bin
def binof(c):
    return np.clip(np.searchsorted(ccv, c) - 1, 0, len(ccv)-1)
AB1 = AB[binof(c1w)]                          # (NP_, NQ)
AB2 = AB[binof(c2w)]
WQ_FLAT = np.full(NQ, 1.0/NQ)
wq_alt = QG**-0.5
WQ_ALT = wq_alt/wq_alt.sum()

# --- likelihood machinery -------------------------------------------------
KNOTS = np.percentile(np.concatenate([c1w, c2w]), [5, 35, 65, 95])
ZN, ZW = np.polynomial.hermite_e.hermegauss(15)
ZW = ZW/np.sqrt(2*np.pi)
GN, GW = np.polynomial.hermite_e.hermegauss(7)
GW = GW/np.sqrt(2*np.pi)
V1P, V2P = s1w**2, s2w**2

def curves(x):
    f, mu1, mu2 = x[0], x[1], x[2]
    rk, snk = x[3:7], x[7:11]
    r1 = np.interp(c1w, KNOTS, rk); r2 = np.interp(c2w, KNOTS, rk)
    sn1 = np.interp(c1w, KNOTS, snk); sn2 = np.interp(c2w, KNOTS, snk)
    V1 = sn1**2 + V1P; V2 = sn2**2 + V2P
    return f, mu1, mu2, r1, r2, V1, V2

def _mix(dd, mu, r, V, ABr, wq, f, gshift=0.0):
    # (NP_, nz): (1-f) S + f A at each z node
    M = dd[:, None] - mu - gshift - r[:, None]*ZN[None, :]
    iv = 1.0/V[:, None]
    S = np.exp(-0.5*M*M*iv)*np.sqrt(iv/(2*np.pi))
    DF = M[:, :, None] - ABr[:, None, :]
    A = np.einsum('q,npq->np', wq,
                  np.exp(-0.5*DF*DF*iv[:, :, None]))*np.sqrt(iv/(2*np.pi))
    return (1.0-f)*S + f*A

def lnL(x, wq, rank2=False):
    if not (0.0 <= x[0] <= 0.6): return 1e15
    if np.any(np.abs(x[1:3]) > 0.3): return 1e15
    if np.any(np.abs(x[3:7]) > 0.8): return 1e15
    if np.any(x[7:11] < 0.03) or np.any(x[7:11] > 0.5): return 1e15
    sg = x[11] if rank2 else 0.0
    if rank2 and not (0.0 <= sg <= 0.3): return 1e15
    f, mu1, mu2, r1, r2, V1, V2 = curves(x)
    if rank2:
        Lp = np.zeros(NP_)
        for gj, wgj in zip(GN, GW):
            m1 = _mix(d1w, mu1, r1, V1, AB1, wq, f, sg*gj)
            m2 = _mix(d2w, mu2, r2, V2, AB2, wq, f, sg*gj)
            Lp += wgj*np.einsum('z,nz->n', ZW, m1*m2)
    else:
        m1 = _mix(d1w, mu1, r1, V1, AB1, wq, f)
        m2 = _mix(d2w, mu2, r2, V2, AB2, wq, f)
        Lp = np.einsum('z,nz->n', ZW, m1*m2)
    return -float(np.sum(np.log(np.maximum(Lp, 1e-300))))

# --- ALS prefit of r(c) from the cross-component covariance field ---------
# E[cov(d1,d2) | colorbins a,b] = r_a r_b: companion terms vanish under
# centering + independence, so the prefit is companion-immune (initializer
# + rank-1 adequacy read; the ML fit is the measurement).
NCB = 12
QB = np.percentile(np.concatenate([c1w, c2w]), np.linspace(0, 100, NCB+1))
QB[0] -= 1e-6; QB[-1] += 1e-6
b1 = np.clip(np.searchsorted(QB, c1w)-1, 0, NCB-1)
b2 = np.clip(np.searchsorted(QB, c2w)-1, 0, NCB-1)
Cm = np.zeros((NCB, NCB)); Nm = np.zeros((NCB, NCB))
mu_b = np.zeros(NCB)
allb = np.concatenate([b1, b2]); alld = np.concatenate([d1w, d2w])
for a in range(NCB):
    mu_b[a] = np.mean(alld[allb == a])
pr = (d1w - mu_b[b1])*(d2w - mu_b[b2])
for a in range(NCB):
    for b in range(NCB):
        m = ((b1 == a) & (b2 == b)) | ((b1 == b) & (b2 == a))
        if m.sum() >= 40:
            Cm[a, b] = np.mean(pr[m]); Nm[a, b] = m.sum()
rv = np.full(NCB, 0.15)
for _ in range(200):
    for a in range(NCB):
        num = np.sum(Nm[a]*Cm[a]*rv); den = np.sum(Nm[a]*rv**2)
        if den > 0: rv[a] = num/den
resid = np.sqrt(np.sum(Nm*(Cm - rv[:, None]*rv[None, :])**2)/max(Nm.sum(), 1))
scale = np.sqrt(np.sum(Nm*Cm**2)/max(Nm.sum(), 1))
qcen = 0.5*(QB[:-1]+QB[1:])
r_init = np.interp(KNOTS, qcen, rv)
P(f"ALS prefit r(c) at knots {np.round(KNOTS,2).tolist()}: "
  f"{np.round(r_init,3).tolist()}; rank-1 residual/scale = "
  f"{resid:.4f}/{scale:.4f} = {resid/max(scale,1e-9):.2f}")

# --- fits -----------------------------------------------------------------
def fit(wq, x0, rank2=False, maxfev=2500):
    res = minimize(lnL, x0, args=(wq, rank2), method='Nelder-Mead',
                   options=dict(maxfev=maxfev, xatol=1e-4, fatol=1e-3,
                                adaptive=True))
    return res.x, -res.fun

def x0_of(f0):
    return np.concatenate([[f0, 0.0, 0.05], r_init, [0.15]*4])

t0 = time.time()
starts = []
for f0 in (0.05, 0.20):
    xh, lh = fit(WQ_FLAT, x0_of(f0))
    starts.append((lh, xh))
    P(f"  start f0={f0}: lnL={lh:.1f}, f_hat={xh[0]:.3f} "
      f"({(time.time()-t0)/60:.1f} min)")
starts.sort(key=lambda t: -t[0])
gap01 = starts[0][0] - starts[1][0]
xb, lb = fit(WQ_FLAT, starts[0][1], maxfev=1500)   # polish
X_R1, L_R1 = xb, lb
P(f"rank-1 MAP (flat-q): lnL={L_R1:.1f}, f={X_R1[0]:.4f}, "
  f"mu=({X_R1[1]:+.3f},{X_R1[2]:+.3f}), r={np.round(X_R1[3:7],3).tolist()}, "
  f"sn={np.round(X_R1[7:11],3).tolist()}; start agreement {gap01:.1f} lnL")

x0f = X_R1.copy(); x0f[0] = 0.0
def lnL_f0(x, wq):
    xx = np.concatenate([[0.0], x])
    return lnL(xx, wq)
res0 = minimize(lnL_f0, x0f[1:], args=(WQ_FLAT,), method='Nelder-Mead',
                options=dict(maxfev=2000, adaptive=True))
L_F0 = -res0.fun
X_R2, L_R2 = fit(WQ_FLAT, np.concatenate([X_R1, [0.05]]), rank2=True,
                 maxfev=1800)
P(f"nested: lnL(f=0)={L_F0:.1f} (dL vs rank-1 = {L_R1-L_F0:+.1f}); "
  f"rank-2 lnL={L_R2:.1f} (dL = {L_R2-L_R1:+.1f}), sigma_g={X_R2[11]:.3f}, "
  f"f_rank2={X_R2[0]:.4f}")
gz3 = (L_R1 >= L_F0 - 0.5) and (L_R2 >= L_R1 - 0.5) and (gap01 <= 2.0) \
      and (X_R1[0] < 0.58)
P(f"GZ3 nesting/convergence -> {'PASS' if gz3 else 'FAIL'}")

# --- injections (GZ1 over-attribution, GZ2 recovery) ----------------------
def simulate(x, f_true, wq, sg_true=0.0, seed=1):
    rg = np.random.default_rng(seed)
    f_, mu1, mu2, r1, r2, V1, V2 = curves(x)
    z = rg.normal(size=NP_); g = rg.normal(size=NP_)
    out = []
    for dd, mu, r, V, ABr in ((d1w, mu1, r1, V1, AB1),
                              (d2w, mu2, r2, V2, AB2)):
        a = np.zeros(NP_)
        hasc = rg.random(NP_) < f_true
        qi = rg.choice(NQ, NP_, p=wq)
        a[hasc] = ABr[hasc, qi[hasc]]
        out.append(mu + r*z + sg_true*g + a + rg.normal(size=NP_)*np.sqrt(V))
    return out

def fit_on(sim1, sim2, x0, maxfev=1800):
    global d1w, d2w
    keep1, keep2 = d1w, d2w
    d1w, d2w = sim1, sim2
    try:
        xh, lh = fit(WQ_FLAT, x0, maxfev=maxfev)
    finally:
        d1w, d2w = keep1, keep2
    return xh, lh

gz_ok = {}
for tag, ftr, sgt, bar in (('GZ1a f=0 rank1-truth', 0.0, 0.0, 0.03),
                           ('GZ1b f=0 rank2-truth', 0.0, 0.05, 0.05),
                           ('GZ2a f=0.10', 0.10, 0.0, None),
                           ('GZ2b f=0.25', 0.25, 0.0, None)):
    sm1, sm2 = simulate(X_R1, ftr, WQ_FLAT, sg_true=sgt,
                        seed=int(1000+100*ftr*100+sgt*1000))
    x0 = X_R1.copy(); x0[0] = max(ftr, 0.05)
    xh, _ = fit_on(sm1, sm2, x0)
    if bar is not None:
        ok_ = xh[0] <= bar
    else:
        ok_ = abs(xh[0]-ftr) <= max(0.03, 0.25*ftr)
    gz_ok[tag] = ok_
    P(f"{tag}: recovered f_hat={xh[0]:.3f} (truth {ftr}) -> "
      f"{'PASS' if ok_ else 'FAIL'}  ({(time.time()-t0)/60:.1f} min)")
PRIMARY = 'rank1'
if not gz_ok['GZ1b f=0 rank2-truth']:
    P("DECISION RULE FIRES: rank-1 over-attributes under rank-2 truth -> "
      "rank-2 PROMOTED to primary; re-running GZ1 against it")
    PRIMARY = 'rank2'
    def fit_on2(sim1, sim2, x0, maxfev=1800):
        global d1w, d2w
        keep1, keep2 = d1w, d2w
        d1w, d2w = sim1, sim2
        try:
            xh, lh = fit(WQ_FLAT, x0, rank2=True, maxfev=maxfev)
        finally:
            d1w, d2w = keep1, keep2
        return xh, lh
    for tag, sgt in (('GZ1a-r2', 0.0), ('GZ1b-r2', 0.05)):
        sm1, sm2 = simulate(X_R1, 0.0, WQ_FLAT, sg_true=sgt, seed=77)
        x0 = np.concatenate([X_R2[:11], [max(sgt, 0.02)]])
        x0[0] = 0.05
        xh, _ = fit_on2(sm1, sm2, x0)
        gz_ok[tag] = xh[0] <= (0.03 if sgt == 0 else 0.05)
        P(f"{tag}: f_hat={xh[0]:.3f} -> {'PASS' if gz_ok[tag] else 'FAIL'}")

X_P = X_R1 if PRIMARY == 'rank1' else X_R2
RANK2 = (PRIMARY == 'rank2')

# --- GZ6 postdiction at the primary MAP -----------------------------------
sgp = X_P[11] if RANK2 else 0.0
sm1, sm2 = simulate(X_P, X_P[0], WQ_FLAT, sg_true=sgp, seed=424242)
smc = (np.abs(sm1) < 0.3) & (np.abs(sm2) < 0.3)
def rho_m(m):
    return rho_att(sm1[m], sm2[m], s1w[m]**2, s2w[m]**2)
dcw = np.abs(c1w - c2w)
rb_m = rho_m(smc & (dcw >= 0.15))
rl_m = rho_m(smc & (dcw < 0.15))
rh_m = rho_m(smc & (dcw >= 0.4))
pf_m = 0.5*(np.mean(sm1 < -0.4) + np.mean(sm2 < -0.4))
i1, i2 = sm1 < -0.4, sm2 < -0.4
inc_m = float((i1 & i2).mean()/max(i1.mean()*i2.mean(), 1e-9))
g6i = abs(rb_m - rho_bar_data) <= 0.05
g6ii = (abs(rl_m - rho_lo_data) <= 0.07) and (abs(rh_m - rho_hi_data) <= 0.07)
g6iii = abs(pf_m - pflag_data) <= 0.02
P(f"GZ6 postdiction: rho_bar model {rb_m:+.3f} vs data "
  f"{rho_bar_data:+.3f} -> {'PASS' if g6i else 'FAIL'}; "
  f"slices ({rl_m:+.3f},{rh_m:+.3f}) vs ({rho_lo_data:+.3f},"
  f"{rho_hi_data:+.3f}) -> {'PASS' if g6ii else 'FAIL'}; "
  f"P(flag) {pf_m:.3f} vs {pflag_data:.3f} -> "
  f"{'PASS' if g6iii else 'FAIL'}; incidence ratio model {inc_m:.2f} vs "
  f"data {inc_data:.2f} (reported, no bar)")

# --- profile over f + q-law envelope --------------------------------------
F_PROF = np.array([0.0, 0.03, 0.06, 0.10, 0.15, 0.22, 0.30])
def profile(wq, xstart):
    prof = np.zeros(len(F_PROF))
    xw = xstart.copy()
    for i, fv in enumerate(F_PROF):
        def lnL_fix(x, wq=wq):
            xx = np.concatenate([[fv], x])
            if RANK2:
                return lnL(xx, wq, rank2=True)
            return lnL(xx, wq)
        r = minimize(lnL_fix, xw[1:], method='Nelder-Mead',
                     options=dict(maxfev=900, adaptive=True))
        prof[i] = -r.fun
        xw = np.concatenate([[fv], r.x])
    return prof
prof_flat = profile(WQ_FLAT, X_P)
P(f"profile (flat-q): {np.round(prof_flat - prof_flat.max(), 1).tolist()} "
  f"on f = {F_PROF.tolist()}  ({(time.time()-t0)/60:.1f} min)")
x0a = X_P.copy()
if RANK2:
    xa, la = fit(WQ_ALT, x0a, rank2=True, maxfev=1800)
else:
    xa, la = fit(WQ_ALT, x0a, maxfev=1800)
prof_alt = profile(WQ_ALT, xa)
P(f"q^-0.5 MAP: f={xa[0]:.4f}; profile: "
  f"{np.round(prof_alt - prof_alt.max(), 1).tolist()}")
lnpi_flat = prof_flat - prof_flat.max()
lnpi_alt = prof_alt - prof_alt.max()
env = np.maximum(lnpi_flat, lnpi_alt)
in1 = F_PROF[(env.max()-env) <= 0.5]
f_lo, f_hi = float(in1.min()), float(in1.max())
edge_flag = env[-1] >= env.max() - 0.5
P(f"BLENDED-axis envelope: f_photo(per component) in [{f_lo:.2f}, "
  f"{f_hi:.2f}] 1-sigma; MAPs {X_P[0]:.3f}/{xa[0]:.3f}"
  + ("  ** GRID-EDGE at 0.30 - extend before quoting **" if edge_flag
     else ""))

# --- part-A host remap, verbatim machinery --------------------------------
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
P(f"HOST-axis prior (LANDED ANCHOR): f_host(per component) in "
  f"[{fh1.min():.2f}, {fh1.max():.2f}] (peak {fh_hat:.2f}) at 1-sigma")
for refv, lab in ((0.10, 'old fence'), (0.16, 'lit-16 anchor'),
                  (0.35, 'posterior cell'), (0.51, 'retracted part-A peak'),
                  (0.69, 'Banik free fit')):
    iv = int(np.argmin(np.abs(FH_GRID-refv)))
    P(f"  ln pi_host({refv:.2f} = {lab}) = {env_fh[iv]-env_fh.max():.1f}")
np.savez('data/stage7jz_prior.npz', fh_grid=FH_GRID, lnpi_host=env_fh,
         f_prof_grid=F_PROF, lnpi_flat=lnpi_flat, lnpi_alt=lnpi_alt,
         f_map_flat=X_P[0], f_map_alt=xa[0], fh_hat=fh_hat,
         primary=PRIMARY, x_map=X_P, knots=KNOTS,
         rho_bar_model=rb_m, rho_bar_data=rho_bar_data)

gates_all = gz0 and gz3 and all(gz_ok.values()) and g6i and g6ii and g6iii \
            and not edge_flag
P("")
P(f"COMPARATORS: retracted part-A f_photo 0.22-0.30 / f_host peak 0.51; "
  f"literature per-component ~0.10 (per-pair 0.166); measured here: "
  f"blended [{f_lo:.2f}, {f_hi:.2f}], host peak {fh_hat:.2f}")
P(f"GATES {'ALL PASS -> the landed anchor SHIPS' if gates_all else 'FAIL '
  'branch -> NEEDS REFINEMENT: anchor NOT shipped (part 2 reads lit only)'}"
  f"  [total {(time.time()-t00)/60:.1f} min]")
