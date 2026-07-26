# -*- coding: utf-8 -*-
"""STAGE 7J-z2c: THE CERTIFICATE RUN (v2c-plus; pre-registered before
execution; sanction chain: v2b's own pre-commit queued "the certificate
run" as the sanctioned next instrument, the round-10 review upgraded it
to v2c-plus — ship a (q, P)-resolved rate, never a scalar — and the
user ordered the queue executed 2026-07-26).

WHAT v2b LEFT OPEN (both numerical-protocol, neither physics): GZS11
d = 2.080 (bar 0.5): the state-11 21^2 q-thinning is under-resolved
where twin pairs' conditional kernels are ~0.04 mag wide.  Start gap
2.8 (bar 2.0): a short-lambda basin.  The no-third-iteration rule held
that day; f_hat was stable 0.166 -> 0.162 -> 0.159 across every class
tried and GZ1c-b acquitted it — what was missing is a CERTIFICATE, not
a number.

THE FIX IS AGAIN EXACT, NOT FINER (the v2b move, applied to q):
  Conditional decomposition of every companion state.  With
  Sigma = [[S11, S12], [S12, S22]]:
    P10 = N(x2; sqrt(S22)) * I1(x1 - (S12/S22) x2),
      I1(t) = int pi(q) N(t - A1(q); sc1) dq,  sc1^2 = S11 - S12^2/S22.
  Under the model's amplitude interpolation (A piecewise-LINEAR on the
  31-knot q-grid) and a piecewise-constant q-density, I1 is a CLOSED
  FORM — a sum of 30 segment terms (w_s/|dA_s|)[Phi((t-A_lo)/sc) -
  Phi((t-A_hi)/sc)] — with NO resolution parameter: the spiked (twin)
  direction that killed v2's GH quadrature and under-resolved v2b's
  state-11 grid is integrated exactly.  P01 symmetric.  P11 keeps ONE
  numeric axis (the outer q2, whose integrand is smooth with width
  >= sqrt(S22) ~ 0.1 mag): M-point Gauss-Legendre, resolution-GATED by
  node doubling (GV1).  The q-density model: flat = uniform on
  [0.1, 1]; alt ~ q^-0.5 (segment masses exact: 2(sqrt(qb)-sqrt(qa)));
  twin-tilt laws ~ 1 + t*1[q >= 0.9] for t in {2, 5} (the GV7
  q-shape deliverable).  Discrete-atom mode (v2b verbatim) is kept as
  a code path for the identity-point regression only.

PROFILED LAMBDA (kills the start-gap direction): lam on the grid
[0.5, 1.0, 1.64, 2.5, 4.0] (1.64 = the v2b MAP), inner Nelder-Mead
over the remaining 12 params, warm-chained outward from 1.64 plus one
cold start per node; the certificate MAP is the profile peak.

GATES (all pre-registered; ANY fail -> quarantine npz, LIT16 stands,
STOP - no further completeness instrument without the user):
  GV0 identity-point regression (GB0w rule): the discrete-mode lnL at
      the v2b MAP reproduces the printed 2973.8 to 0.1.
  GV1 resolution: |lnL(M=31) - lnL(M=62)| <= 0.5 at the cont MAP
      (the only numeric axis left).
  GV2 protocol: top-two start gap at lam* <= 2.0 (5 starts); lam*
      interior in the grid (one pre-authorized extension if edge).
  GV3 nesting: L_cont(MAP) >= L_cont(v2b params) - 0.01.
  GV4 injections (cont mode, kernel-truth skies): GZ1c-a f = 0 ->
      f_hat <= 0.03; GZ2a/b f = 0.10/0.25 recovered within
      max(0.03, 25%); GZ1c-b kernel-blind diagnostic reported.
  GV6 postdiction (v2b bars verbatim): rho bar +-0.05, slices +-0.07,
      windowed flag +-0.02; incidence no-bar.
  GV7 q-SHAPE DELIVERABLE (round-10 requirement; no bar, shipped):
      lnL profile across the four q-laws (flat / q^-0.5 / twin t=2 /
      twin t=5) each at its own f MAP; the host conversion PER LAW;
      and the flat-q-equivalent statement carried as the [0.10, 0.39]
      conversion band (bin-7j-qmoments) stored in the npz - the
      shipped anchor DECLARES that imposing it on the flat-q fcomp
      axis at face precision is invalid (sigma* = 0.02 measured).
  SHIP RULE: GV0-GV6 all pass -> stage7jz_prior.npz version
      'v2c-cert' (fields backward-compatible with the readers + the
      GV7 block); any fail -> stage7jz_prior_v2c_unshipped.npz.

Output: data/stage7jz2c_cert.txt
"""
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr
from astropy.io import fits

t00 = time.time()
L = []
def P(s):
    print(s, flush=True)
    L.append(s)
    with open('data/stage7jz2c_cert.txt', 'w') as f:
        f.write("\n".join(L) + "\n")

# --- data (v2b verbatim) --------------------------------------------------
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
P(f"window: {NP_} pairs; windowed flag {pf_win_data:.3f}")

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
KNOTS = np.percentile(np.concatenate([c1w, c2w]), [5, 35, 65, 95])
V1P, V2P = s1w**2, s2w**2

# --- q-density laws (piecewise; segment masses exact) ---------------------
def seg_masses(law, t=0.0):
    qa, qb = QG[:-1], QG[1:]
    if law == 'flat':
        w = qb - qa
    elif law == 'alt':
        w = 2*(np.sqrt(qb) - np.sqrt(qa))
    elif law == 'twin':
        w = (qb - qa)*(1.0 + t*((qa + qb)/2 >= 0.9))
    return w/w.sum()

def gl_nodes(M, law, t=0.0):
    xg, wg = np.polynomial.legendre.leggauss(M)
    q = 0.5*(xg+1)*(1.0-0.10)+0.10
    w = wg*0.5*(1.0-0.10)
    if law == 'flat':
        dens = np.ones(M)/0.9
    elif law == 'alt':
        dens = q**-0.5/(2*(1-np.sqrt(0.1)))
    elif law == 'twin':
        nrm = 0.9 + t*0.1
        dens = (1.0 + t*(q >= 0.9))/nrm
    return q, w*dens

# --- exact-cont likelihood ------------------------------------------------
def unpack2(x):
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
    return f, mu1, mu2, S11, S22, S12

def seg_int(t, Anod, sc, wseg):
    # exact integral of the piecewise model: sum_s (w_s/|dA_s|)
    #   [Phi((t-A_lo)/sc) - Phi((t-A_hi)/sc)]; degenerate segs -> point
    Alo = np.minimum(Anod[..., :-1], Anod[..., 1:])
    Ahi = np.maximum(Anod[..., :-1], Anod[..., 1:])
    dA = Ahi - Alo
    scb = sc[..., None]
    tb = t[..., None]
    flat = dA < 1e-9
    Phi = ndtr((tb - Alo)/scb) - ndtr((tb - Ahi)/scb)
    dens = np.where(flat, 0.0, wseg/np.maximum(dA, 1e-9))
    pt = wseg/np.sqrt(2*np.pi)/scb*np.exp(
        -0.5*((tb - 0.5*(Alo+Ahi))/scb)**2)
    return np.sum(np.where(flat, pt, dens*Phi), axis=-1)

BOUNDS_MSG = 1e15
def lnL_cont(x, law='flat', t=0.0, M=31):
    if not (0.0 <= x[0] <= 0.6): return BOUNDS_MSG
    if np.any(np.abs(x[1:3]) > 0.3): return BOUNDS_MSG
    if np.any(np.abs(x[3:7]) > 0.8): return BOUNDS_MSG
    if np.any(x[7:11] < 0.02) or np.any(x[7:11] > 0.5): return BOUNDS_MSG
    if not (0.0 <= x[11] <= 0.3): return BOUNDS_MSG
    if not (0.05 <= x[12] <= 5.0): return BOUNDS_MSG
    f, mu1, mu2, S11, S22, S12 = unpack2(x)
    x1 = d1w - mu1
    x2 = d2w - mu2
    det = np.maximum(S11*S22 - S12*S12, 1e-12)
    P00 = 1.0/(2*np.pi*np.sqrt(det))*np.exp(
        -0.5*(S22*x1*x1 - 2*S12*x1*x2 + S11*x2*x2)/det)
    wseg = seg_masses(law, t)
    sc1 = np.sqrt(np.maximum(S11 - S12*S12/S22, 1e-12))
    sc2 = np.sqrt(np.maximum(S22 - S12*S12/S11, 1e-12))
    n2 = 1.0/np.sqrt(2*np.pi*S22)*np.exp(-0.5*x2*x2/S22)
    n1 = 1.0/np.sqrt(2*np.pi*S11)*np.exp(-0.5*x1*x1/S11)
    P10 = n2*seg_int(x1 - (S12/S22)*x2, AB1, sc1, wseg)
    P01 = n1*seg_int(x2 - (S12/S11)*x1, AB2, sc2, wseg)
    qn, qw = gl_nodes(M, law, t)
    A2n = _A2GL[(law, t, M)]
    x2b = x2[:, None] - A2n
    n2q = 1.0/np.sqrt(2*np.pi*S22)[:, None]*np.exp(
        -0.5*x2b*x2b/S22[:, None])
    tq = x1[:, None] - (S12/S22)[:, None]*x2b
    I1q = seg_int(tq, AB1[:, None, :], sc1[:, None], wseg)
    P11 = np.einsum('m,nm->n', qw, n2q*I1q)
    Lp = ((1-f)**2*P00 + f*(1-f)*(P10+P01) + f*f*P11)
    return -float(np.sum(np.log(np.maximum(Lp, 1e-300))))

_A2GL = {}
def prep_gl(law, t, M):
    if (law, t, M) not in _A2GL:
        qn, _ = gl_nodes(M, law, t)
        _A2GL[(law, t, M)] = np.stack(
            [np.interp(qn, QG, AB2[i]) for i in range(NP_)])
    return _A2GL[(law, t, M)]

# --- discrete mode (v2b verbatim, for GV0 only) ---------------------------
def lnL_disc(x, wq, idx11=IDX11):
    f, mu1, mu2, S11, S22, S12 = unpack2(x)
    det = np.maximum(S11*S22 - S12*S12, 1e-12)
    a11, a22, a12 = S22/det, S11/det, S12/det
    nrm = 1.0/(2*np.pi*np.sqrt(det))
    x1 = d1w - mu1
    x2 = d2w - mu2
    P00 = nrm*np.exp(-0.5*(a11*x1*x1 - 2*a12*x1*x2 + a22*x2*x2))
    X1 = x1[:, None] - AB1
    Q10 = (a11[:, None]*X1*X1 - 2*a12[:, None]*X1*x2[:, None]
           + (a22*x2*x2)[:, None])
    P10 = nrm*np.einsum('q,nq->n', wq, np.exp(-0.5*Q10))
    X2 = x2[:, None] - AB2
    Q01 = ((a11*x1*x1)[:, None] - 2*a12[:, None]*x1[:, None]*X2
           + a22[:, None]*X2*X2)
    P01 = nrm*np.einsum('q,nq->n', wq, np.exp(-0.5*Q01))
    w11 = wq[idx11]/wq[idx11].sum()
    A1s, A2s = X1[:, idx11], X2[:, idx11]
    Q11 = (a11[:, None, None]*(A1s*A1s)[:, :, None]
           - 2*a12[:, None, None]*A1s[:, :, None]*A2s[:, None, :]
           + a22[:, None, None]*(A2s*A2s)[:, None, :])
    P11 = nrm*np.einsum('p,q,npq->n', w11, w11, np.exp(-0.5*Q11))
    Lp = ((1-f)**2*P00 + f*(1-f)*(P10+P01) + f*f*P11)
    return -float(np.sum(np.log(np.maximum(Lp, 1e-300))))

# AMENDMENT A3 (logged after the first launch died pre-fit; nothing
# quoted): the original GV0 evaluated lnL at the PRINT-ROUNDED v2b
# vector against the print-rounded reference with a 0.1 bar - a
# bar-design miss (rounding alone moves lnL by O(0.1)).  The v2b
# quarantined npz stores the full-precision x_map; the regression now
# evaluates THAT, bar 0.15 (= the reference's print-rounding envelope).
# Same launch: np.trapz -> np.trapezoid (NumPy 2.x API).
pz_v2b = np.load('data/stage7jz_prior_v2b_unshipped.npz',
                 allow_pickle=True)
X_V2B = np.asarray(pz_v2b['x_map'], float)
P(f"v2b full-precision MAP loaded: f={X_V2B[0]:.4f}, lam={X_V2B[12]:.3f}")
l_disc = -lnL_disc(X_V2B, WQ_FLAT)
gv0 = abs(l_disc - 2973.8) <= 0.15
P(f"GV0 identity-point regression: lnL_disc(v2b x_map) = {l_disc:.2f} "
  f"(printed 2973.8, bar 0.15) -> {'PASS' if gv0 else 'FAIL'}")
# GV0b: f = 0 companion-free identity — cont and disc likelihoods are
# the SAME closed form there (P00 only); catches unpack/P00 wiring.
prep_gl('flat', 0.0, 31)
x_f0 = X_V2B.copy(); x_f0[0] = 0.0
d_f0 = abs(lnL_cont(x_f0, 'flat', 0.0, 31) - lnL_disc(x_f0, WQ_FLAT))
gv0b = d_f0 <= 1e-6
P(f"GV0b f=0 cont/disc identity: d = {d_f0:.2e} -> "
  f"{'PASS' if gv0b else 'FAIL'}")
# GV0c: the segment-erf core vs a brute-force dense Riemann sum of the
# SAME piecewise model (3001 nodes), spot-checked on 200 pairs.
_, _, _, S11c, S22c, S12c = unpack2(X_V2B)
sc1c = np.sqrt(np.maximum(S11c - S12c*S12c/S22c, 1e-12))
tt = (d1w - X_V2B[1]) - (S12c/S22c)*(d2w - X_V2B[2])
sub = slice(0, 200)
ws = seg_masses('flat')
I_erf = seg_int(tt[sub], AB1[sub], sc1c[sub], ws)
qd_ = np.linspace(0.10, 1.0, 3001)
I_rie = np.zeros(200)
for i in range(200):
    Ad = np.interp(qd_, QG, AB1[sub][i])
    dens = np.ones(len(qd_))/0.9
    I_rie[i] = np.trapezoid(dens/np.sqrt(2*np.pi)/sc1c[sub][i]
                            * np.exp(-0.5*((tt[sub][i]-Ad)/sc1c[sub][i])**2),
                            qd_)
rel = float(np.max(np.abs(I_erf - I_rie)/np.maximum(I_rie, 1e-30)))
gv0c = rel <= 1e-3
P(f"GV0c segment-erf vs dense Riemann (200 pairs): max rel d = "
  f"{rel:.2e} -> {'PASS' if gv0c else 'FAIL'}")
gv0 = gv0 and gv0b and gv0c
if not gv0:
    P("GV0 FAIL -> STOP (pre-registered)")
    raise SystemExit(1)

# --- profiled-lambda fit (cont, flat law) ---------------------------------
LAM_GRID = [0.5, 1.0, 1.64, 2.5, 4.0]
prep_gl('flat', 0.0, 31)

def fit12(lam, x0_12, maxfev=700):
    def lf(x):
        return lnL_cont(np.concatenate([x, [lam]]), 'flat', 0.0, 31)
    r = minimize(lf, x0_12, method='Nelder-Mead',
                 options=dict(maxfev=maxfev, xatol=1e-4, fatol=1e-3,
                              adaptive=True))
    return r.x, -r.fun

t0 = time.time()
prof = {}
x_chain = X_V2B[:12].copy()
i164 = LAM_GRID.index(1.64)
order = [i164] + [i for i in range(len(LAM_GRID)) if i != i164]
xs_at = {}
for i in order:
    lam = LAM_GRID[i]
    xw, lw = fit12(lam, x_chain, maxfev=700)
    xc = X_V2B[:12].copy(); xc[0] = 0.08
    xcold, lcold = fit12(lam, xc, maxfev=700)
    if lcold > lw:
        xw, lw = xcold, lcold
    prof[lam] = lw
    xs_at[lam] = xw
    if lam == 1.64:
        x_chain = xw.copy()
    P(f"  lam={lam}: lnL={lw:.1f}, f={xw[0]:.4f} "
      f"({(time.time()-t0)/60:.1f} min)")
lam_star = max(prof, key=prof.get)
interior = LAM_GRID.index(lam_star) not in (0, len(LAM_GRID)-1)
P(f"profile peak lam* = {lam_star} "
  f"({'interior' if interior else 'EDGE - one extension authorized'})")
if not interior:
    ext = 0.25 if lam_star == LAM_GRID[0] else 5.0
    xw, lw = fit12(ext, xs_at[lam_star], maxfev=700)
    prof[ext] = lw; xs_at[ext] = xw
    P(f"  extension lam={ext}: lnL={lw:.1f}, f={xw[0]:.4f}")
    lam_star = max(prof, key=prof.get)
    lams = sorted(prof)
    interior = lams.index(lam_star) not in (0, len(lams)-1)

# start dispersion at lam* (GV2)
best = []
for f0 in (0.08, 0.16, 0.25, 0.10):
    x0 = X_V2B[:12].copy(); x0[0] = f0
    xh, lh = fit12(lam_star, x0, maxfev=900)
    best.append((lh, xh))
xh, lh = fit12(lam_star, xs_at[lam_star], maxfev=900)
best.append((lh, xh))
best.sort(key=lambda z: -z[0])
gap = best[0][0] - best[1][0]
X_C = np.concatenate([best[0][1], [lam_star]])
L_C = best[0][0]
gv2 = (gap <= 2.0) and interior
P(f"v2c MAP (cont, profiled-lam): lnL={L_C:.1f}, f={X_C[0]:.4f}, "
  f"lam*={lam_star}; top-two gap {gap:.2f} -> "
  f"{'PASS' if gv2 else 'FAIL'} ({(time.time()-t0)/60:.1f} min)")

# GV1 resolution (outer doubling), GV3 nesting
prep_gl('flat', 0.0, 62)
d31_62 = abs(-lnL_cont(X_C, 'flat', 0.0, 31) + lnL_cont(X_C, 'flat', 0.0, 62))
gv1 = d31_62 <= 0.5
P(f"GV1 outer-node doubling at MAP: d = {d31_62:.3f} -> "
  f"{'PASS' if gv1 else 'FAIL'}")
l_v2b_cont = -lnL_cont(X_V2B, 'flat', 0.0, 31)
gv3 = L_C >= l_v2b_cont - 0.01
P(f"GV3 nesting: L_cont(MAP)={L_C:.1f} >= L_cont(v2b params)="
  f"{l_v2b_cont:.1f} -> {'PASS' if gv3 else 'FAIL'}")
P(f"f_hat chain: v1 0.166 -> v2b 0.159 -> v2c (exact-cont) {X_C[0]:.3f}")

# --- injections (cont mode; GV4) ------------------------------------------
def simulate_c(x, f_true, seed=1):
    rg = np.random.default_rng(seed)
    _, mu1, mu2 = x[0], x[1], x[2]
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
        qd = 0.1 + 0.9*rg.random(NP_)
        idx = np.where(hasc)[0]
        a[idx] = np.array([np.interp(qd[i], QG, ABr[i]) for i in idx])
        out.append(mu + r*(np.sqrt(k)*z + np.sqrt(1-k)*eta) + sg*g + a
                   + rg.normal(size=NP_)*np.sqrt(sn**2 + sp**2))
    return out

def fit_on_c(sim1, sim2, x0, lam, maxfev=900, pin_lam=None):
    global d1w, d2w
    keep1, keep2 = d1w, d2w
    d1w, d2w = sim1, sim2
    try:
        lamv = pin_lam if pin_lam is not None else lam
        xh, lh = fit12(lamv, x0[:12], maxfev=maxfev)
    finally:
        d1w, d2w = keep1, keep2
    return np.concatenate([xh, [lamv]]), lh

gv4 = {}
sm1, sm2 = simulate_c(X_C, 0.0, seed=911)
x0 = X_C.copy(); x0[0] = 0.05
xh, _ = fit_on_c(sm1, sm2, x0, lam_star)
gv4['GZ1c-a'] = xh[0] <= 0.03
P(f"GZ1c-a kernel-truth f=0: f_hat={xh[0]:.3f} -> "
  f"{'PASS' if gv4['GZ1c-a'] else 'FAIL'} ({(time.time()-t0)/60:.1f} min)")
xh1, _ = fit_on_c(sm1, sm2, x0, lam_star, pin_lam=4.99)
P(f"GZ1c-b kernel-blind diagnostic: f_hat={xh1[0]:.3f}")
for tag, ftr in (('GZ2a', 0.10), ('GZ2b', 0.25)):
    sm1, sm2 = simulate_c(X_C, ftr, seed=int(2000+ftr*100))
    x0 = X_C.copy(); x0[0] = ftr
    xh, _ = fit_on_c(sm1, sm2, x0, lam_star)
    ok_ = abs(xh[0]-ftr) <= max(0.03, 0.25*ftr)
    gv4[tag] = ok_
    P(f"{tag} f={ftr}: recovered {xh[0]:.3f} -> "
      f"{'PASS' if ok_ else 'FAIL'} ({(time.time()-t0)/60:.1f} min)")

# --- GV6 postdiction ------------------------------------------------------
sm1, sm2 = simulate_c(X_C, X_C[0], seed=424243)
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
P(f"GV6: bar {rb_m:+.3f} vs {rho_bar_data:+.3f} -> "
  f"{'PASS' if g6i else 'FAIL'}; slices ({rl_m:+.3f},{rh_m:+.3f}) vs "
  f"({rho_lo_data:+.3f},{rho_hi_data:+.3f}) -> "
  f"{'PASS' if g6ii else 'FAIL'}; flag {pf_m:.3f} vs {pf_win_data:.3f} "
  f"-> {'PASS' if g6iii else 'FAIL'}; incidence {inc_m:.2f} vs "
  f"{inc_data:.2f} (no bar)")

# --- GV7: the q-shape deliverable + profiles ------------------------------
F_PROF = np.array([0.0, 0.06, 0.11, 0.17, 0.23, 0.30])
def profile_c(law, t, xstart, lam):
    prep_gl(law, t, 31)
    prof = np.zeros(len(F_PROF))
    xw = xstart.copy()
    for i, fv in enumerate(F_PROF):
        def lf(x):
            return lnL_cont(np.concatenate([[fv], x, [lam]]), law, t, 31)
        r = minimize(lf, xw[1:12], method='Nelder-Mead',
                     options=dict(maxfev=500, adaptive=True))
        prof[i] = -r.fun
        xw = np.concatenate([[fv], r.x])
    return prof

qlaw_rows = []
for law, t, tag in (('flat', 0.0, 'flat'), ('alt', 0.0, 'q^-0.5'),
                    ('twin', 2.0, 'twin t=2'), ('twin', 5.0, 'twin t=5')):
    prep_gl(law, t, 31)
    def lf(x):
        return lnL_cont(np.concatenate([x, [lam_star]]), law, t, 31)
    r = minimize(lf, X_C[:12], method='Nelder-Mead',
                 options=dict(maxfev=900, adaptive=True))
    qlaw_rows.append((tag, law, t, -r.fun, r.x[0]))
    P(f"q-law {tag}: lnL={-r.fun:.1f}, f_hat={r.x[0]:.4f} "
      f"({(time.time()-t0)/60:.1f} min)")
lbest = max(r[3] for r in qlaw_rows)
for tag, law, t, lv, fv in qlaw_rows:
    P(f"  q-shape table: {tag}: dlnL={lv-lbest:+.1f}, f={fv:.3f}")

prof_flat = profile_c('flat', 0.0, X_C, lam_star)
lnpi_flat = prof_flat - prof_flat.max()
P(f"profile flat-q: {np.round(lnpi_flat,1).tolist()} on {F_PROF.tolist()}")
xalt = [r for r in qlaw_rows if r[0] == 'q^-0.5'][0]
xa0 = X_C.copy(); xa0[0] = xalt[4]
prof_alt = profile_c('alt', 0.0, xa0, lam_star)
lnpi_alt = prof_alt - prof_alt.max()
P(f"profile q^-0.5: {np.round(lnpi_alt,1).tolist()}")
env = np.maximum(lnpi_flat, lnpi_alt)
in1 = F_PROF[(env.max()-env) <= 0.5]
f_lo, f_hi = float(in1.min()), float(in1.max())
edge_flag = env[-1] >= env.max() - 0.5
P(f"BLENDED envelope: f in [{f_lo:.2f}, {f_hi:.2f}]"
  + ("  ** GRID EDGE **" if edge_flag else ""))

# host conversion (v2b verbatim MC) per q-law
dists = 1000.0/plx[ok]
NB2 = 200000
rgb = np.random.default_rng(31415)
Ms_b = 0.6+1.8*rgb.random(NB2)
Mh_b = 0.5*Ms_b
Pyr_b = 10**rgb.normal(5.03, 2.28, NB2)/365.25
db = rgb.choice(np.asarray(dists, float), NB2)
FH_GRID = np.arange(0.0, 0.9001, 0.01)
env_fh = np.full(len(FH_GRID), -1e9)
qdraws = {'flat': 0.1+0.9*rgb.random(NB2),
          'qalt': rgb.choice(QG, NB2, p=WQ_ALT)}
for qtag, qd, lnpiA in (('flat', qdraws['flat'], lnpi_flat),
                        ('qalt', qdraws['qalt'], lnpi_alt)):
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
P(f"HOST-axis (v2c): f_host per component in [{fh1.min():.2f}, "
  f"{fh1.max():.2f}] (peak {fh_hat:.2f})")

gates_all = gz0 and gv0 and gv1 and gv2 and gv3 and all(gv4.values()) \
            and g6i and g6ii and g6iii and not edge_flag
P("")
if gates_all:
    np.savez('data/stage7jz_prior.npz', fh_grid=FH_GRID,
             lnpi_host=env_fh, f_prof_grid=F_PROF, lnpi_flat=lnpi_flat,
             lnpi_alt=lnpi_alt, f_map_flat=X_C[0], fh_hat=fh_hat,
             version='v2c-cert', x_map=X_C, knots=KNOTS,
             qlaw_tags=np.array([r[0] for r in qlaw_rows]),
             qlaw_lnl=np.array([r[3] for r in qlaw_rows]),
             qlaw_f=np.array([r[4] for r in qlaw_rows]),
             conv_band=np.array([0.10, 0.39]),
             conv_note=np.array(
                 'imposing lnpi_host on the flat-q fcomp axis at face '
                 'precision is INVALID (sigma*=0.02, bin-7j-qmoments); '
                 'operative use = conversion-widened profile, factors '
                 'g in [0.33, 1.30]'))
    P("GATES ALL PASS -> THE CERTIFICATE SHIPS "
      "(data/stage7jz_prior.npz, version v2c-cert, q-resolved block)")
else:
    np.savez('data/stage7jz_prior_v2c_unshipped.npz', fh_grid=FH_GRID,
             lnpi_host=env_fh, x_map=X_C, version='v2c-unshipped')
    P("GATE FAIL -> nothing ships; LIT16 stands; quarantined npz "
      "written; STOP per pre-registration")
P(f"[total {(time.time()-t00)/60:.1f} min]")
