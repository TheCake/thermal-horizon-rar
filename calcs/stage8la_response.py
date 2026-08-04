"""STAGE 8L-a — THE RESPONSE MODEL, decisive half (pre-reg 5b1ee6f,
committed BEFORE this run).

(a) Analytic leakage/residual curves for a circular photocenter
orbit against a 5-parameter fit over baseline T (1-D projection,
continuous uniform sampling, phase-averaged):
    L(u) = 3 (sin u - u cos u) / u^3        (PM leakage; = 3 j1(u)/u)
    R^2(u) = (1/2)[1 - sinc^2(u) - 3 (sin u - u cos u)^2 / u^4]
    with u = pi T / P.  Limits: L->1, R->0 (curvature) as P>>T;
    L->0, R->1/sqrt2 as P<<T.  Gated numerically (GL0).
(b) THE SELF-CONSISTENCY NUMBER: the seed-31 collapse world
    (fcomp 0.35/0.50) through the RUWE forward -> predicted f_hot
    vs the measured 0.090 (bars + map in the pre-reg).
(c) The S2 forward: predicted Delta-vt(hot-cold) at fcomp = 0.1
    vs the measured +0.17 (descriptive consistency).
GL1: the derived L vs the legacy S = min(1, P/17.8 yr) — REPORTED,
interpretation deferred to 8L-b.  GL2 validity: the operative-world
forward must postdict the sky's 0.090 within [0.05, 0.15].
Output: data/stage8la_response.txt
"""
import numpy as np
from astropy.io import fits

L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

P("8L-a THE RESPONSE MODEL (pre-reg 5b1ee6f; bars + map locked "
  "before any run)")
P("")

# ---------- (a) the curves + GL0 ----------
def Lk(u):
    u = np.asarray(u, dtype=float)
    small = u < 1e-4
    out = np.where(small, 1.0 - u*u/10.0,
                   3.0*(np.sin(u) - u*np.cos(u))/np.maximum(u, 1e-300)**3)
    return out

def R2(u):
    # amendment 1(i): series branch widened to u < 0.05 + core
    # clipped at 0 (cancellation NaNs at extreme P, run-1)
    u = np.asarray(u, dtype=float)
    small = u < 0.05
    s = np.sinc(u/np.pi)          # sin(u)/u
    core = 0.5*(1.0 - s*s
                - 3.0*(np.sin(u) - u*np.cos(u))**2
                / np.maximum(u, 1e-300)**4)
    return np.where(small, u**4/90.0, np.maximum(core, 0.0))

# GL0: numeric least-squares regression at sample u
rng = np.random.default_rng(7)
gl0_max = 0.0
for u in (0.3, 1.0, 2.0, 5.0, 12.0):
    # midpoint grid (endpoint-inclusive linspace has an O(2/N)
    # second-moment offset vs continuous - caught by GL0 first run)
    t = np.linspace(-0.5, 0.5, 20000, endpoint=False) + 0.5/20000
    om = 2.0*u
    Ls, Rs = [], []
    for ph in np.linspace(0, 2*np.pi, 48, endpoint=False):
        x = np.cos(om*t + ph)
        A = np.vstack([np.ones_like(t), t]).T
        coef, *_ = np.linalg.lstsq(A, x, rcond=None)
        res = x - A@coef
        vinst = -om*np.sin(ph)
        if abs(vinst) > 1e-3:
            Ls.append(coef[1]/vinst)
        Rs.append(np.mean(res**2))
    dL = abs(np.mean(Ls) - Lk(u))
    dR = abs(np.mean(Rs) - R2(u))
    gl0_max = max(gl0_max, float(dL), float(dR))
lim_ok = (abs(Lk(1e-6)-1) < 1e-9 and Lk(200.0) < 1e-3
          and R2(1e-6) < 1e-12 and abs(R2(300.0)-0.5) < 1e-3)
P(f"GL0 (numeric-vs-analytic, 5 u-nodes x 48 phases): max|d| = "
  f"{gl0_max:.2e}; limits {'OK' if lim_ok else 'BAD'} -> "
  f"{'PASS' if gl0_max <= 1e-6 and lim_ok else 'FAIL'}")
assert gl0_max <= 1e-6 and lim_ok

# ---------- constants (scout-verified where possible) ----------
# T: EDR3 baseline 34 months = 1037.9 d (Lindegren+21, A&A 649, A2).
# Period-scaling shape (excess peaks at P ~ T; P^(2/3) rise below;
# strong suppression above) confirmed by Belokurov+20 (MNRAS 496,
# 1922).  sigma_AL per-OBSERVATION values below are standard-ballpark
# (consistent with end-of-mission 0.02-0.03 mas at G 9-14 over
# N ~ 43-249 obs); the scout could not extract granular per-G values,
# so THE ABSOLUTE SCALE IS JUDGED BY GL2 (the sky-postdiction gate),
# per the pre-reg.  RUWE = UWE/u0(G, color) (GAIA-C3-TN-LU-LL-124-01).
# The legacy S = min(1, P/17.8 yr): NOT FOUND in the literature at
# scout level - provenance likely internal; GL1 reports the physics.
T_BASE = 2.83                      # yr
SIG_G = np.array([8.0, 9.0, 11.0, 13.0, 15.0, 17.0])
SIG_AL = np.array([0.25, 0.25, 0.30, 0.45, 0.90, 2.30])
def sig_al(G):
    return np.interp(np.clip(G, SIG_G[0], SIG_G[-1]), SIG_G, SIG_AL)
NOBS_EFF = 1.0   # per-observation convention; absolute scale judged by GL2
P(f"constants: T = {T_BASE} yr; sigma_AL nodes {SIG_AL.tolist()} mas "
  f"at G {SIG_G.tolist()} (scout-sourced; GL2 judges the absolute "
  f"scale)")
P("")

# ---------- GL1: L vs the legacy S (REPORTED) ----------
P("GL1 (reported; interpretation deferred to 8L-b): derived L(P/T) "
  "vs legacy S = min(1, P/17.8):")
for Pyr in (1.0, 2.0, 3.0, 6.0, 10.0, 17.8, 30.0):
    u = np.pi*T_BASE/Pyr
    P(f"  P = {Pyr:5.1f} yr: L = {float(Lk(u)):.3f}, "
      f"S_legacy = {min(1.0, Pyr/17.8):.3f}, R = "
      f"{float(np.sqrt(R2(u))):.3f}")
P("")

# ---------- catalog hosts ----------
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
vx = d['pmra2']-d['pmra1']; vy = d['pmdec2']-d['pmdec1']
vmag = np.hypot(vx, vy)
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot = np.interp(MG1, MG_T, MS_T)+np.interp(MG2, MG_T, MS_T)
s_kau = sep/1e3
vc = 0.9417*np.sqrt(Mtot/np.maximum(s_kau, 1e-9))
vt = (4.74047/plx*vmag)/vc
WIDE = ok & (s_kau >= 6)
iW = np.where(WIDE)[0]
assert int(ok.sum()) == 14071 and len(iW) == 1437
P(f"hosts: {int(ok.sum())} pairs, {len(iW)} WIDE")
P("")

# ---------- the RUWE forward ----------
def fwd_fhot(fcomp, nmc=40, seed=11, base='unity'):
    """Predicted per-component hot fraction in WIDE from the model's
    own companion law (q flat 0.1-1; logP N(5.03, 2.28) days; valid
    a_in < 130 AU and < a_s/5; wfac = |q/(1+q) - l/(1+l)|), photo-
    center amplitude a_phot = wfac*a_in*plx [mas], RUWE inflation
    sqrt(base^2 + (R(u)*a_phot/sigma_AL)^2) > 1.25.  base='unity'
    (primary, conservative-low) or 'p50' (1.059, reported)."""
    rng = np.random.default_rng(seed)
    b0 = 1.0 if base == 'unity' else 1.059
    thr2 = 1.25**2 - b0**2
    if thr2 <= 0:
        return 0.0
    hot, tot = 0, 0
    for comp, (Gm, Mg) in enumerate(((G1m, MG1), (G2m, MG2))):
        Gc = Gm[iW]; Mc = np.interp(Mg[iW], MG_T, MS_T)
        pc = plx[iW]; a_s = s_kau[iW]*1e3
        for _ in range(nmc):
            act = rng.random(len(iW)) < fcomp
            n = int(act.sum())
            if n == 0:
                continue
            q = 0.1+0.9*rng.random(n)
            logP = rng.normal(5.03, 2.28, n)
            P_yr = 10**logP/365.25
            M_h = Mc[act]
            a_in = (M_h*(1+q)*P_yr**2)**(1/3)
            valid = (a_in < 130.0) & (a_in < a_s[act]/5.0)
            MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
            MGs = np.interp(-np.clip(q*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
            l_ = 10**(-0.4*(MGs-MGp))
            wfac = np.abs(q/(1+q) - l_/(1+l_))
            a_phot = wfac*a_in*pc[act]          # mas
            u = np.pi*T_BASE/np.maximum(P_yr, 1e-6)
            excess = np.sqrt(R2(u))*a_phot
            hot += int(np.sum(valid & (excess**2
                                       > thr2*sig_al(Gc[act])**2)))
            tot += len(iW)
    return hot/max(tot, 1)

f_op = fwd_fhot(0.10)
P(f"[REPORTED CONTEXT - the absolute forward, DISCLOSED-INVALID "
  f"per amendment 1: singles-base-dominated observable] "
  f"f_hot(fcomp = 0.10) = {f_op:.3f} vs measured 0.090 "
  f"(run-1 GL2 abort record)")
f_op50 = fwd_fhot(0.10, base='p50')
P(f"  (base = P50 variant: {f_op50:.3f}; primary = unity base, "
  f"conservative-low)")
P("")

f35 = fwd_fhot(0.35)
f50 = fwd_fhot(0.50)
fcol = 0.5*(f35+f50)     # seed-31 posterior mass ~equal at the two
P(f"THE SELF-CONSISTENCY NUMBER: f_hot(collapse) = {f35:.3f} "
  f"(fcomp 0.35) / {f50:.3f} (0.50); mass-weighted = {fcol:.3f} "
  f"vs measured 0.090")
f35p, f50p = fwd_fhot(0.35, base='p50'), fwd_fhot(0.50, base='p50')
P(f"  (base = P50 variant: {f35p:.3f}/{f50p:.3f})")
P("")

# ---------- GL2' license: the catalog-cut-repaired S2 forward ----
sigvW = sigv[WIDE]
ceilW = 2.978/np.sqrt(s_kau[WIDE]) + 2.8284*sigvW   # km/s
rng = np.random.default_rng(23)
vt_cold_pool = vt[WIDE]
nmc = 40
dvts = []
for _ in range(nmc):
    act1 = rng.random(len(iW)) < 0.10
    act2 = rng.random(len(iW)) < 0.10
    leak = np.zeros(len(iW))
    ishot = np.zeros(len(iW), dtype=bool)
    for act, Gm, Mg in ((act1, G1m, MG1), (act2, G2m, MG2)):
        n = int(act.sum())
        if n == 0:
            continue
        q = 0.1+0.9*rng.random(n)
        logP = rng.normal(5.03, 2.28, n)
        P_yr = 10**logP/365.25
        M_h = np.interp(Mg[iW][act], MG_T, MS_T)
        a_in = (M_h*(1+q)*P_yr**2)**(1/3)
        valid = (a_in < 130.0) & (a_in < (s_kau[iW][act]*1e3)/5.0)
        MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        MGs = np.interp(-np.clip(q*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        l_ = 10**(-0.4*(MGs-MGp))
        wfac = np.abs(q/(1+q) - l_/(1+l_))
        v_orb = 29.78*np.sqrt(M_h*(1+q)/np.maximum(a_in, 1e-3))
        u = np.pi*T_BASE/np.maximum(P_yr, 1e-6)
        lv = Lk(u)*wfac*v_orb*valid          # leaked speed, km/s
        a_phot = wfac*a_in*plx[iW][act]
        ex = np.sqrt(R2(u))*a_phot
        hotc = valid & (ex**2 > (1.25**2-1.0)*sig_al(Gm[iW][act])**2)
        li = np.zeros(len(iW)); li[act] = lv
        hi = np.zeros(len(iW), dtype=bool); hi[act] = hotc
        leak = np.hypot(leak, li)
        ishot |= hi
    proj = np.sqrt(1.0 - rng.uniform(-1, 1, len(iW))**2)
    v_new = np.hypot(vt_cold_pool*vc[WIDE], proj*leak)   # km/s
    keep = v_new <= ceilW           # amendment 1(ii): catalog survival
    vt_new = v_new/vc[WIDE]
    hk, ck = ishot & keep, (~ishot) & keep
    if hk.sum() > 3 and ck.sum() > 3:
        dvts.append(float(vt_new[hk].mean() - vt_new[ck].mean()))
d_pred = float(np.mean(dvts))
gl2p = 0.5*0.174 <= d_pred <= 2.0*0.174
P(f"GL2' LICENSE (catalog-cut-repaired S2 forward, fcomp = 0.10, "
  f"{nmc} reps): predicted dvt(hot-cold) = {d_pred:+.3f} vs "
  f"measured +0.174; [0.5x, 2x] -> "
  f"{'PASS - mapping sky-calibrated' if gl2p else 'FAIL - full stop'}")
P("")

# ---------- 8L-a2: THE FAKER-CONDITIONAL (pre-reg 3d98c7a) --------
import csv
cens = [r for r in csv.DictReader(open('data/ceiling_pairs.csv'))
        if r['census_corr'] == 'True']
okidx = np.where(ok)[0]
cidx = []
for r in cens:
    m = (np.abs(s_kau[okidx]-float(r['s_kAU'])) < 0.01) \
      & (np.abs(Mtot[okidx]-float(r['Mtot_Msun'])) < 0.01) \
      & (np.abs(vc[okidx]-float(r['vc_kms'])) < 0.001)
    j = okidx[m]
    assert len(j) == 1
    cidx.append(int(j[0]))
OBS_HOT = 2                     # the 8K measurement (ruwe > 1.25)
rng = np.random.default_rng(97)
NF = 400_000
pj, nfk = [], []
P("8L-a2 THE FAKER-CONDITIONAL (per census pair; the model's own "
  "companion law restricted to draws that could fake the pair):")
for j in cidx:
    M_h = Mtot[j]/2.0
    q = 0.1+0.9*rng.random(NF)
    logP = rng.normal(5.03, 2.28, NF)
    P_yr = 10**logP/365.25
    a_in = (M_h*(1+q)*P_yr**2)**(1/3)
    valid = (a_in < 130.0) & (a_in < sep[j]/5.0)
    MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
    MGs = np.interp(-np.clip(q*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
    l_ = 10**(-0.4*(MGs-MGp))
    wfac = np.abs(q/(1+q) - l_/(1+l_))
    v_orb = 29.78*np.sqrt(M_h*(1+q)/np.maximum(a_in, 1e-3))
    u = np.pi*T_BASE/np.maximum(P_yr, 1e-6)
    proj = np.sqrt(1.0 - rng.uniform(-1, 1, NF)**2)
    leak = Lk(u)*wfac*v_orb*proj
    vobs = vt[j]*vc[j]
    ceilj = 2.978/np.sqrt(s_kau[j]) + 2.8284*sigv[j]
    fak = valid & (leak >= 0.5*vobs) & (leak <= ceilj)
    Gh = np.where(rng.random(NF) < 0.5, G1m[j], G2m[j])
    a_phot = wfac*a_in*plx[j]
    hot = (np.sqrt(R2(u))*a_phot)**2 > (1.25**2-1.0)*sig_al(Gh)**2
    n = int(fak.sum())
    p = float(np.mean(hot[fak])) if n > 0 else float('nan')
    pj.append(p); nfk.append(n)
    P(f"  s = {s_kau[j]:6.2f} kAU, vt = {vt[j]:.3f}, vobs = "
      f"{vobs:.3f} km/s: n_faker = {n}, P(hot | faker) = {p:.3f}")
pbar = float(np.nanmean(pj))
lown = min(nfk)
# exact Poisson-binomial P(<= OBS_HOT hot | all nine are fakers)
dp = np.zeros(11); dp[0] = 1.0
for p in pj:
    dp[1:] = dp[1:]*(1-p) + dp[:-1]*p
    dp[0] *= (1-p)
p_le2 = float(dp[:OBS_HOT+1].sum())
P(f"mean P(hot | faker) over the nine = {pbar:.3f} (min n_faker = "
  f"{lown}); exact Poisson-binomial P(<= {OBS_HOT} of 9 hot | all "
  f"fakers) = {p_le2:.2e}")
P("")

# ---------- verdict per the 8L-a2 map ----------
if not gl2p:
    P("==> 8L-a2 VERDICT: GL2' FAIL - full stop; hold ~55%; the "
      "NSS leg decides (pre-stated).")
else:
    if pbar >= 0.6 and p_le2 <= 0.05:
        v = ("FAKER-LOUD: companions capable of faking the census "
             "pairs would have lit the RUWE flags (mean P(hot|faker)"
             f" = {pbar:.2f}; observing <= 2 of 9 hot has P = "
             f"{p_le2:.1e}) - the collapse account of the nine is "
             "object-level DEAD; the S1 conditional is DISCHARGED "
             "at the census; per the map: anomaly-real ~55% -> ~57%")
    elif pbar < 0.3:
        v = (f"FAKER-QUIET (mean P = {pbar:.2f}): the fakers can "
             "hide from RUWE; the census blind spot is REAL at "
             "correction grade; per the map: anomaly-real ~55% -> "
             "~50%; the NSS leg decides")
    else:
        v = (f"GRAY (mean P = {pbar:.2f}, P(<=2) = {p_le2:.1e}): "
             "hold ~55%; the NSS leg decides")
    P(f"==> 8L-a2 VERDICT (locked bars + map): {v}")

with open('data/stage8la_response.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage8la_response.txt")
