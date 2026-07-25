"""
STAGE 7J part B: the marginalized-companion v7 fit — the decisive
instrument for 7I's COMPANION-DIRECTION verdict (design from the Opus
review exchange, adopted before the result made it urgent).

Everything the reviews contested is FREED JOINTLY, and the companion
fraction is constrained by the part-A completeness measurement instead
of a fence:
  WR_GRID    = [0.10, 0.20, 0.30, 0.40, 0.50]   (7I rode the 0.3 edge)
  FCOMP_GRID = [0.0, 0.10, 0.20, 0.35, 0.50, 0.70]  (spans Banik's 0.69)
  FPM_GRID   = [1.2, 1.5, 1.8]
Two samples: 'full' (14,071) and 'strict' (the 7I mask, 8,047). Two
result grades per seed-law: PROF (free profile, no prior — the
Banik-style full-freedom row, 4N continuity) and MARG (marginal over
eta/wr/fcomp/ffly/fpm with the part-A completeness prior on fcomp —
full sample: pi(f_true); strict sample: pi(residual rate r(f),
photometric-only = upper bound; grid values outside the prior's
support get ln pi = -1e9). Discrete-cell marginalization, stated.

GATE GB0 (regression, per seed-law, before anything is quoted): the
sub-cube at (wr <= 0.3, fcomp <= 0.1, fpm = 1.5) must reproduce the
stored baseline rows — data/stage4r_summary.txt (full) /
data/stage7i_s.txt (strict) — within |d a_hat| <= 0.02 and
|d dlnL(Newton)| <= 2.0. FAIL -> abort that sample's verdict.

BARS (pre-registered; evaluated on seed means; base seeds 31/101,
auto-extension to 202/303/404/505 on the ambiguous branch):
  FULL sample:
    ANCHOR-CONFIRMED  if a_marg >= 0.9 AND dNewton_marg >= +30, both laws.
    COMPANION-WIN     if a_marg <= 0.7 OR dNewton_marg <= +15, either law.
    else AMBIGUOUS -> extend, re-evaluate, then report as-is.
  STRICT sample:
    RESOLVED-INSTRUMENTAL if a_marg >= 0.9 AND dNewton_marg >= +25, both.
    COMPANION-CONFIRMED   if a_marg <= 0.7, both laws.
    else MIDDLE -> extend, re-evaluate, then report as-is.
  wr watch: flag EXTERNAL-INCONSISTENT if the marginal posterior puts
  >= 50% of its wr mass at wr >= 0.4 (beyond the Hwang-implied 20-22%).

SYNTHESIS QUADRANTS (credence moves pre-committed, NOTES-grade):
  (A) FULL confirmed & STRICT resolved      -> the 7I collapse was
      instrumentation (fence+grids); anomaly-real back to ~70%.
  (B) FULL confirmed & STRICT companion     -> unexplained subsample
      split; ~55-60%; carried openly; next = 7K/7L/DR4.
  (C) FULL companion-win (either strict)    -> the completeness-honest
      fit no longer needs the boost at detection grade; ~35-45%; the
      binary section is rebuilt around the model-light channels and the
      census; the mimicry ledger is reframed.
  (D) ambiguous after extension             -> reported as-is, carried.

PRE-VERDICT ADDITIONS (review exchange; logged before any MARG verdict
is read):
  PRIOR AXIS (amendment 3): the part-A prior is remapped to the HOST
  fraction — the axis the v7 fcomp parameter actually lives on (the
  blended-axis prior was restrictive-misspecified, biasing alpha HIGH).
  MARG rows produced under the earlier blended-axis prior are STALE and
  are purged before the verdict is read; cubes are prior-independent.
  STRICTPOW (the power gate; run as sample 'strictpow'): inject an
  alpha = 1.18 truth at the measured residual host rate into the strict
  selection (truth pop seed 777, count draws 888, multinomial at the
  strict bins' observed totals). GP bar: a_marg in [0.9, 1.5] and
  dN_marg >= +25 -> POWER-OK; otherwise POWER-FAIL and the
  strict-companion quadrant carries no weight (pre-stated).
  RIDGE + STRESS per seed-law: the (wr, fcomp) joint-posterior
  correlation (does the fit trade the two absorbers?) and the stress
  conditionals a_prof at fcomp = 0.5 and 0.7 (the conservative
  high-end quote a referee will ask for).

AMENDMENT 4 (logged pre-verdict; the third axis): the part-A host
inversion gives f_host = 0.42-0.57 (peak 0.51) — reproducing the known
solar-type multiplicity (~46%, Raghavan 2010) — which exposes the
model-side half of the units problem: 3K's fcomp <= 0.1 preference
under a TRUE ~0.5 host rate means the v7 companion amplitude law is
too strong per companion. The missing physics is identified: Gaia
astrometry tracks the PHOTOCENTER, whose wobble is
|q/(1+q) - l/(1+l)| * v_orb (l = L2/L1 from the same mass-magnitude
table) — the luminosity term cancels the mass wobble, exactly for
twins (the as-published law gives twins the maximum). The corrected
instrument ('photo' mode, the DECISIVE one for the bars):
  - photocenter amplitude law in the companion sector (mh unchanged —
    dynamical mass is real regardless of light);
  - amplitude-scale nuisance KW_GRID = [0.7, 1.0, 1.4] profiled
    jointly (residual amplitude-model error absorbed in BOTH
    directions; the S period-smearing factor is kept as published and
    NOT relitigated);
  - GATE GB0p (cross-mode regression): the photo-mode cube at
    fcomp = 0 must equal the as-published cube at fcomp = 0 to 1e-6
    (companions off => identical models);
  - the as-published ('raw') cubes are kept and reported as the
    continuity row (PROF only; the host prior does not apply to the
    raw fcomp axis);
  - photo-mode extension on ambiguity is +2 seeds (202/303), not +4
    (cost honesty: ~45-55 min per seed-law).

argv: <sample: full|strict|strictpow> [photo] <seeds...>
Appends rows to data/stage7j_<sample>[_photo].txt; verdict appended to
data/stage7j_verdict.txt once both laws of all requested seeds exist.
"""
import os
import re
import sys
import time

import numpy as np
from astropy.io import fits

SAMPLE = sys.argv[1] if len(sys.argv) > 1 else 'full'
assert SAMPLE in ('full', 'strict', 'strictpow'), SAMPLE
_rest = sys.argv[2:]
AMP = 'photo' if (_rest and _rest[0] == 'photo') else 'raw'
if AMP == 'photo':
    _rest = _rest[1:]
SEEDS_REQ = [int(x) for x in _rest] or [31, 101]
TAG = '_photo' if AMP == 'photo' else ''
OUT = f'data/stage7j_{SAMPLE}{TAG}.txt'
KW_GRID = np.array([0.7, 1.0, 1.4]) if AMP == 'photo' else np.array([1.0])

src = open('calcs/stage2b_population.py').read()
ns = {}
exec(src.split('# ---------- validation')[0], ns)
run = ns['run']
GM = 4*np.pi**2; A0_CAN = 7.99e-7

def load_tab(path):
    t = np.load(path); y, b = t[0][::-1], t[1][::-1]
    lny = np.log(y); lny_u = np.linspace(lny[0], lny[-1], 512)
    return lny_u, np.interp(lny_u, lny, b)
LNY_U, TAB_S = load_tab('data/efe_boost_simple_g1p2.npy')
_,     TAB_B = load_tab('data/efe_boost_be_g1p2.npy')
LNY0, DLNY = LNY_U[0], LNY_U[1]-LNY_U[0]

d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1,plx2 = d['parallax1'],d['parallax2']
eplx1,eplx2 = d['parallax_error1'],d['parallax_error2']
sep,Rch = d['sep_AU'],d['R_chance_align']
G1m,G2m = d['phot_g_mean_mag1'],d['phot_g_mean_mag2']
plx = 0.5*(plx1+plx2)
MG1 = G1m+5*np.log10(np.maximum(plx1,1e-6))-10
MG2 = G2m+5*np.log10(np.maximum(plx2,1e-6))-10
sigv = 4.74047/plx*np.sqrt(d['pmra_error1']**2+d['pmdec_error1']**2
                           +d['pmra_error2']**2+d['pmdec_error2']**2)
ok = (Rch<0.01)&(plx1>5)&(plx2>5)&(plx1/np.maximum(eplx1,1e-6)>20) \
   &(plx2/np.maximum(eplx2,1e-6)>20)&(np.abs(plx1-plx2)<3*np.hypot(eplx1,eplx2)) \
   &(sep>200)&(sep<50000)&(MG1>2.6)&(MG1<14.2)&(MG2>2.6)&(MG2<14.2)&(sigv<0.03)
if SAMPLE in ('strict', 'strictpow'):
    ok0 = int(ok.sum())
    ok = ok & np.load('data/stage7i_strictmask.npy')
    print(f"strict sample: {int(ok.sum())}/{ok0} pairs", flush=True)
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot_d = np.interp(MG1[ok],MG_T,MS_T)+np.interp(MG2[ok],MG_T,MS_T)
s_d = sep[ok]/1e3
# corrected velocities (4R convention, 6G formulas)
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
sn_ = np.maximum(np.hypot(sx_, sy_), 1e-12)
sux_, suy_ = sx_/sn_, sy_/sn_
_cn = list(d.columns.names)
def _pick(*names):
    for n in names:
        if n in _cn: return np.asarray(d[n], dtype=np.float64)
    raise KeyError(names)
r1_ = _pick('radial_velocity1', 'dr2_radial_velocity1', 'rv1')
r2_ = _pick('radial_velocity2', 'dr2_radial_velocity2', 'rv2')
try:
    er1_ = _pick('radial_velocity_error1', 'dr2_radial_velocity_error1')
    er2_ = _pick('radial_velocity_error2', 'dr2_radial_velocity_error2')
except KeyError:
    er1_ = np.full(len(r1_), 2.0); er2_ = np.full(len(r2_), 2.0)
h1_, h2_ = np.isfinite(r1_), np.isfinite(r2_)
w1_ = np.where(h1_, 1.0/np.maximum(er1_, 0.5)**2, 0.0)
w2_ = np.where(h2_, 1.0/np.maximum(er2_, 0.5)**2, 0.0)
rvs_ = (np.where(h1_, r1_, 0.0)*w1_ + np.where(h2_, r2_, 0.0)*w2_) \
       / np.maximum(w1_+w2_, 1e-12)
th_ = sep/(2.06265e8/plx)
pmcor_ = rvs_*th_*plx/4.74047
vx_ = d['pmra2']-d['pmra1'] + pmcor_*sux_
vy_ = d['pmdec2']-d['pmdec1'] + pmcor_*suy_
vt_d = (4.74047/plx[ok]*np.hypot(vx_[ok], vy_[ok]))/(0.9417*np.sqrt(Mtot_d/s_d))
e_vx = np.sqrt(d['pmra_error1']**2+d['pmra_error2']**2)
e_vy = np.sqrt(d['pmdec_error1']**2+d['pmdec_error2']**2)
cosg = np.abs(sx_*vx_+sy_*vy_)/np.maximum(np.hypot(sx_,sy_)*np.hypot(vx_,vy_),
                                          1e-12)
gam_d = np.degrees(np.arccos(np.clip(cosg, 0, 1)))[ok]
sig_c = ((4.74047/plx)*np.sqrt(0.5*(e_vx**2+e_vy**2)))
SBINS = [(0.2,2),(2,6),(6,20),(20,50)]
SC2 = (np.array([0.63, 3.46, 11.0, 31.6])/31.6)**2
VE = np.logspace(np.log10(0.02), np.log10(6.0), 21)
GE = np.linspace(0, 90, 7)
NV, NG = 20, 6
data_2d, noise_pool = [], []
for b in SBINS:
    m = (s_d>=b[0])&(s_d<b[1])
    h,_,_ = np.histogram2d(np.clip(vt_d[m],0.021,5.9), gam_d[m], bins=[VE, GE])
    data_2d.append(h.astype(float))
    noise_pool.append(sig_c[ok][m])
vcen = np.sqrt(VE[:-1]*VE[1:]); gcen = 0.5*(GE[:-1]+GE[1:])
FLY = np.zeros((NV, NG))
for i in range(NV):
    for j in range(NG):
        FLY[i,j] = (1.0 if (0.5 <= vcen[i] <= 3.0) else 0.0) \
                   * max(30.0-gcen[j], 0.0)
FLY /= max(FLY.sum(), 1e-12)
UNI = np.ones((NV, NG))/(NV*NG)
sig_ok = sig_c[ok]
vc_ok = 0.9417*np.sqrt(Mtot_d/s_d)
UNI_B, FLY_B = [], []
for b in SBINS:
    m = (s_d>=b[0])&(s_d<b[1])
    cutp = 2.978/np.sqrt(s_d[m]) + 2.8284*sig_ok[m]
    acc = np.array([(vcen[i]*vc_ok[m] <= cutp).mean() for i in range(NV)])
    for tpl, store in ((UNI, UNI_B), (FLY, FLY_B)):
        t = tpl*acc[:,None]
        store.append(t/max(t.sum(), 1e-12))

N = 500_000
A_GRID  = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID  = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FC0_GRID = np.array([0.10])
FFLY_GRID = np.array([0.05, 0.10])
FPM_GRID = np.array([1.2, 1.5, 1.8])

def build_pop(seed):
    rng = np.random.default_rng(seed)
    p = {}
    u = rng.random(N); lo, hi, g = 0.15, 60.0, -0.6
    p['a_s']  = ((lo**g+u*(hi**g-lo**g))**(1/g))*1e3
    p['u_e']  = rng.random(N)
    p['psi0'] = rng.random(N)*2*np.pi
    nrm = rng.normal(size=(N,3)); nrm /= np.linalg.norm(nrm,axis=1,keepdims=True)
    p['f_ip'] = np.sqrt(np.clip(1-nrm[:,0]**2,0,1))
    p['M_s']  = 0.6+1.8*rng.random(N)
    p['uph']  = rng.random(N)
    xhat = np.zeros((N,3)); xhat[:,0]=1
    ef = xhat-nrm*nrm[:,[0]]
    ef /= np.maximum(np.linalg.norm(ef,axis=1,keepdims=True),1e-12)
    p['ef'] = ef; p['e2'] = np.cross(nrm,ef)
    los = rng.normal(size=(N,3)); los /= np.linalg.norm(los,axis=1,keepdims=True)
    p['los'] = los
    p['u_mix'] = rng.random(N)
    p['pick'] = [rng.integers(0, max(len(q_),1), N) for q_ in noise_pool]
    p['gn1'], p['gn2'] = rng.normal(size=N), rng.normal(size=N)
    M_h = 0.5*p['M_s']
    p['comp'] = {}
    for k in (1, 2):
        q = 0.1+0.9*rng.random(N)
        logP = rng.normal(5.03, 2.28, N)
        P_yr = 10**logP/365.25
        a_in = (M_h*(1+q)*P_yr**2)**(1/3)
        valid = (a_in < 130.0) & (a_in < p['a_s']/5.0)
        v_orb = 29.78*np.sqrt(M_h*(1+q)/np.maximum(a_in,1e-3))
        S = np.minimum(1.0, P_yr/17.8)
        if AMP == 'photo':
            # photocenter wobble: astrometry follows the light, not the
            # mass — |q/(1+q) - l/(1+l)|, l = L2/L1 from the mass-mag
            # table; exactly zero for twins
            MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
            MGs = np.interp(-np.clip(q*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
            l_ = 10**(-0.4*(MGs-MGp))
            wfac = np.abs(q/(1+q) - l_/(1+l_))
        else:
            wfac = q/(1+q)
        w = wfac*v_orb*S/4.74047*valid
        wd = rng.normal(size=(N,3)); wd /= np.linalg.norm(wd,axis=1,keepdims=True)
        p['comp'][k] = dict(w=w, wd=wd, uc=rng.random(N), mh=q*M_h*valid)
    return p

def e_of(p, eta, wr):
    al = np.interp(np.log10(p['a_s']), np.log10([100,500,1000,50000]),
                   [0.6, 1.0, eta, eta])
    e_pow = 0.95*p['u_e']**(1/(1+al))
    e_rad = 0.9+0.095*p['u_e']
    return np.where(p['u_mix'] < wr, e_rad, e_pow)

def vp_c(p, e_s, tab_a):
    a_s, M_s = p['a_s'], p['M_s']
    rp, ra = a_s*(1-e_s), a_s*(1+e_s)
    xg, wg = np.polynomial.legendre.leggauss(32)
    lo_, hi_ = np.log(rp), np.log(ra)
    lr = 0.5*(hi_-lo_)[:,None]*xg[None,:] + 0.5*(hi_+lo_)[:,None]
    r = np.exp(lr)
    gN = GM*M_s[:,None]/r**2
    bst = np.interp(np.log(gN/A0_CAN), LNY_U, tab_a, right=1.0)
    dPhi = np.sum(wg[None,:]*bst*gN*r, axis=1)*0.5*(hi_-lo_)
    return np.sqrt(np.maximum(2*dPhi/(1-(rp/ra)**2), 0))

def lnL_point(p, o):
    ef, e2, los = p['ef'], p['e2'], p['los']
    s3 = o[:,0,None]*ef+o[:,1,None]*e2
    v3 = o[:,2,None]*ef+o[:,3,None]*e2
    ssky = s3-los*np.sum(s3*los,axis=1,keepdims=True)
    vsky = v3-los*np.sum(v3*los,axis=1,keepdims=True)
    smag = np.linalg.norm(ssky,axis=1)
    b1 = ssky/np.maximum(smag[:,None],1e-12)
    b2 = np.cross(los, b1)
    b2 /= np.maximum(np.linalg.norm(b2,axis=1,keepdims=True),1e-12)
    vpar = np.sum(vsky*b1,axis=1)
    vper = np.sum(vsky*b2,axis=1)
    s_kau = smag/1e3
    out = np.zeros((len(FCOMP_GRID), len(FC0_GRID), len(FFLY_GRID),
                    len(FPM_GRID), len(KW_GRID)))
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500 or len(noise_pool[bi]) == 0: continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = noise_pool[bi][p['pick'][bi][idx] % len(noise_pool[bi])]/4.74047
        for fi, fcm in enumerate(FCOMP_GRID):
            cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
            mh_tot = np.zeros(len(idx))
            for k in (1, 2):
                c = p['comp'][k]
                act = c['uc'][idx] < fcm
                mh_tot += act*c['mh'][idx]
                cvp += act*c['w'][idx]*c['wd'][idx,0]
                cvq += act*c['w'][idx]*c['wd'][idx,1]
            boost = np.sqrt(1+mh_tot/p['M_s'][idx])
            for ki, kwv in enumerate(KW_GRID):
                vp_b = vpar[idx] + kwv*cvp
                vq_b = vper[idx] + kwv*cvq
                for pi, fpm in enumerate(FPM_GRID):
                    vp_n = vp_b*boost + p['gn1'][idx]*sg0*fpm
                    vq_n = vq_b*boost + p['gn2'][idx]*sg0*fpm
                    vmag = np.hypot(vp_n, vq_n)
                    keep = vmag*4.74047 <= (2.978/np.sqrt(s_kau[idx])
                                            + 2.8284*sg0*4.74047)
                    vtn = (vmag/vc)[keep]
                    gmn = np.degrees(np.arccos(np.clip(
                        np.abs(vp_n[keep])/np.maximum(vmag[keep],1e-12),
                        0, 1)))
                    h,_,_ = np.histogram2d(np.clip(vtn,0.021,5.9), gmn,
                                           bins=[VE, GE])
                    p0 = np.maximum(h/max(h.sum(),1), 1e-5); p0 /= p0.sum()
                    for ci, fc in enumerate(FC0_GRID):
                        for yi, ff in enumerate(FFLY_GRID):
                            wch = min(fc*SC2[bi], 0.5)
                            wfl = min(ff*SC2[bi], 0.5)
                            wtot = min(wch+wfl, 0.6)
                            mixc = (wch*UNI_B[bi] + wfl*FLY_B[bi])/(wch+wfl)
                            pp = (1-wtot)*p0 + wtot*mixc
                            out[fi, ci, yi, pi, ki] += \
                                np.sum(data_2d[bi]*np.log(pp))
    return out

def forward_pp(p, o, fcm, fpm, fc, ff):
    ef, e2, los = p['ef'], p['e2'], p['los']
    s3 = o[:,0,None]*ef+o[:,1,None]*e2
    v3 = o[:,2,None]*ef+o[:,3,None]*e2
    ssky = s3-los*np.sum(s3*los,axis=1,keepdims=True)
    vsky = v3-los*np.sum(v3*los,axis=1,keepdims=True)
    smag = np.linalg.norm(ssky,axis=1)
    b1 = ssky/np.maximum(smag[:,None],1e-12)
    b2 = np.cross(los, b1)
    b2 /= np.maximum(np.linalg.norm(b2,axis=1,keepdims=True),1e-12)
    vpar = np.sum(vsky*b1,axis=1)
    vper = np.sum(vsky*b2,axis=1)
    s_kau = smag/1e3
    pps = []
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = noise_pool[bi][p['pick'][bi][idx] % len(noise_pool[bi])]/4.74047
        vp_b = vpar[idx].copy(); vq_b = vper[idx].copy()
        mh_tot = np.zeros(len(idx))
        for k in (1, 2):
            c = p['comp'][k]
            act = c['uc'][idx] < fcm
            mh_tot += act*c['mh'][idx]
            vp_b += act*c['w'][idx]*c['wd'][idx,0]
            vq_b += act*c['w'][idx]*c['wd'][idx,1]
        boost = np.sqrt(1+mh_tot/p['M_s'][idx])
        vp_n = vp_b*boost + p['gn1'][idx]*sg0*fpm
        vq_n = vq_b*boost + p['gn2'][idx]*sg0*fpm
        vmag = np.hypot(vp_n, vq_n)
        keep = vmag*4.74047 <= (2.978/np.sqrt(s_kau[idx])
                                + 2.8284*sg0*4.74047)
        vtn = (vmag/vc)[keep]
        gmn = np.degrees(np.arccos(np.clip(
            np.abs(vp_n[keep])/np.maximum(vmag[keep],1e-12), 0, 1)))
        h,_,_ = np.histogram2d(np.clip(vtn,0.021,5.9), gmn, bins=[VE, GE])
        p0 = np.maximum(h/max(h.sum(),1), 1e-5); p0 /= p0.sum()
        wch = min(fc*SC2[bi], 0.5); wfl = min(ff*SC2[bi], 0.5)
        wtot = min(wch+wfl, 0.6)
        mixc = (wch*UNI_B[bi] + wfl*FLY_B[bi])/(wch+wfl)
        pps.append((1-wtot)*p0 + wtot*mixc)
    return pps

def P(s):
    print(s, flush=True)
    with open(OUT, 'a') as f:
        f.write(s+"\n")

# --- the completeness prior on FCOMP_GRID ---------------------------------
pr = np.load('data/stage7j_prior.npz')
if SAMPLE == 'full':
    xg, lp = pr['f_grid'], pr['lnpi_full']
else:
    xg, lp = pr['r_grid'], pr['lnpi_strict']
LNPI = np.full(len(FCOMP_GRID), -1e9)
inr = (FCOMP_GRID >= xg.min()) & (FCOMP_GRID <= xg.max())
LNPI[inr] = np.interp(FCOMP_GRID[inr], xg, lp)
P(f"STAGE 7J-B {SAMPLE} seeds {SEEDS_REQ}: wr={WR_GRID.tolist()}, "
  f"fcomp={FCOMP_GRID.tolist()}, fpm={FPM_GRID.tolist()}; "
  f"ln pi(fcomp) = {np.round(LNPI, 2).tolist()}")

if SAMPLE == 'strictpow':
    R_TR = float(np.load('data/stage7j_prior.npz')['r_host_hat'])
    pt = build_pop(777)
    e_t = e_of(pt, 1.3, 0.2)
    tab_t = 1.0 + 1.18*(TAB_S-1.0)
    vp_t = vp_c(pt, e_t, tab_t)
    ot = run(pt['a_s'], e_t, pt['psi0'], pt['f_ip'], pt['M_s'],
             pt['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_t, lny0=LNY0,
             dlny=DLNY, vp=vp_t)
    pps = forward_pp(pt, ot, R_TR, 1.5, 0.10, 0.05)
    rgs = np.random.default_rng(888)
    for bi in range(len(SBINS)):
        n_obs = int(data_2d[bi].sum())
        cnt = rgs.multinomial(n_obs, (pps[bi]/pps[bi].sum()).ravel())
        data_2d[bi] = cnt.reshape(NV, NG).astype(float)
    P(f"strictpow: injected alpha=1.18 truth at r_host={R_TR:.2f} "
      f"(truth pop 777, count draws 888, strict-bin totals)")

# --- GB0 references -------------------------------------------------------
ROWRE = re.compile(r"seed (\d+) (simple|BE): a_hat=([0-9.]+) \(grid [0-9.]+, "
                   r"interior=(\w+)\), dlnL\(Newton\)=([+-][0-9.]+), "
                   r"wr=([0-9.]+)")
REF = {}
if SAMPLE != 'strictpow':
    REFF = ('data/stage4r_summary.txt' if SAMPLE == 'full'
            else 'data/stage7i_s.txt')
    for m in ROWRE.finditer(open(REFF).read()):
        s_, law, ah, _, dn, _wr = m.groups()
        REF[(law, int(s_))] = (float(ah), float(dn))

prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
def profile_of(cb):
    prof = np.nanmax(cb, axis=(1,2,3,4,5,6,7))
    imax = int(np.nanargmax(prof))
    ahat = A_GRID[imax]
    if 0 < imax < len(A_GRID)-1:
        x = A_GRID[imax-1:imax+2]; y = prof[imax-1:imax+2]
        c2_, c1_, _ = np.polyfit(x, y, 2)
        if c2_ < 0: ahat = -c1_/(2*c2_)
    return prof, ahat, imax

MROWRE = re.compile(r"MARG " + SAMPLE + r" seed (\d+) (simple|BE): "
                    r"a_marg=([0-9.]+), dN_marg=([+-][0-9.]+)")
def have_rows():
    out = {}
    if os.path.exists(OUT):
        for m in MROWRE.finditer(open(OUT).read()):
            out[(m.group(2), int(m.group(1)))] = (float(m.group(3)),
                                                  float(m.group(4)))
    return out

def run_seed(seed):
    t0 = time.time()
    p = build_pop(seed)
    for law, TAB in (("simple", TAB_S), ("BE", TAB_B)):
        cpath = f'data/stage7j_cube_{SAMPLE}{TAG}_{seed}_{law}.npy'
        if os.path.exists(cpath):
            cube = np.load(cpath)
            if cube.ndim == 7:      # legacy raw cube without the K axis
                cube = cube[..., None]
        else:
            cube = np.full((len(A_GRID), len(E_GRID), len(WR_GRID),
                            len(FCOMP_GRID), len(FC0_GRID), len(FFLY_GRID),
                            len(FPM_GRID), len(KW_GRID)), np.nan)
            newt_cache = {}
            for ai, al in enumerate(A_GRID):
                tab_a = 1.0 + al*(TAB-1.0)
                for ei, eta in enumerate(E_GRID):
                    for wi, wr in enumerate(WR_GRID):
                        if al == 0.0 and law == "BE" and (ei, wi) in newt_cache:
                            cube[ai, ei, wi] = newt_cache[(ei, wi)]
                            continue
                        e_s = e_of(p, eta, wr)
                        vp = vp_c(p, e_s, tab_a) if al > 0 else None
                        mode = 5 if al > 0 else 1
                        kw = dict(a0=A0_CAN, tab=tab_a, lny0=LNY0, dlny=DLNY,
                                  vp=vp) if al > 0 else {}
                        o = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'],
                                p['uph'], 8, 2500, mode, **kw)
                        cube[ai, ei, wi] = lnL_point(p, o)
                        if al == 0.0:
                            newt_cache[(ei, wi)] = cube[ai, ei, wi]
            np.save(cpath, cube)
        # Newton cache trick above only fills BE from simple within one
        # process run; cube files are per-law so reruns are consistent.
        cb = cube + prior_eta[None,:,None,None,None,None,None,None]
        if AMP == 'photo':
            # GB0p: cross-mode regression at fcomp = 0 (companions off
            # => the two amplitude laws are the same model)
            rawp = f'data/stage7j_cube_{SAMPLE}_{seed}_{law}.npy'
            if os.path.exists(rawp):
                rw = np.load(rawp)
                if rw.ndim == 7:
                    rw = rw[..., None]
                dmax = float(np.nanmax(np.abs(
                    cube[:, :, :, 0:1, :, :, :, 0:1]
                    - rw[:, :, :, 0:1, :, :, :, 0:1])))
                P(f"GB0p {SAMPLE} seed {seed} {law}: max|photo-raw| at "
                  f"fcomp=0 = {dmax:.2e} -> "
                  f"{'PASS' if dmax <= 1e-3 else 'FAIL'}")
                if dmax > 1e-3:
                    P(f"ABORT {SAMPLE} seed {seed} {law}: GB0p failed")
                    continue
            else:
                P(f"GB0p {SAMPLE} seed {seed} {law}: raw cube absent - "
                  f"SKIPPED (disclosed)")
        if SAMPLE != 'strictpow' and AMP != 'photo':
            # GB0: legacy sub-cube (wr<=0.3, fcomp<=0.1, fpm=1.5)
            sub = cb[:, :, :3, :2, :, :, 1:2, 0:1]
            _, ah_sub, _ = profile_of(sub)
            prof_sub = np.nanmax(sub, axis=(1,2,3,4,5,6,7))
            dn_sub = float(np.nanmax(prof_sub) - prof_sub[0])
            ra, rd = REF[(law, seed)]
            gb0 = (abs(ah_sub-ra) <= 0.02) and (abs(dn_sub-rd) <= 2.0)
            P(f"GB0 {SAMPLE} seed {seed} {law}: sub-cube a_hat={ah_sub:.2f} "
              f"(ref {ra:.2f}), dN={dn_sub:+.1f} (ref {rd:+.1f}) -> "
              f"{'PASS' if gb0 else 'FAIL'}")
            if not gb0:
                P(f"ABORT {SAMPLE} seed {seed} {law}: GB0 failed")
                continue
        # PROF: full freedom, no prior
        prof, ahat, imax = profile_of(cb)
        best = np.unravel_index(np.nanargmax(cb), cb.shape)
        P(f"PROF {SAMPLE} seed {seed} {law}: a_hat={ahat:.2f} "
          f"(interior={0<imax<len(A_GRID)-1}), "
          f"dN={float(np.nanmax(prof)-prof[0]):+.1f}, "
          f"wr={WR_GRID[best[2]]}, fcomp={FCOMP_GRID[best[3]]}, "
          f"fpm={FPM_GRID[best[6]]}")
        # MARG: discrete-cell marginal with the completeness prior
        cbp = cb + LNPI[None, None, None, :, None, None, None, None]
        m0 = np.nanmax(cbp)
        with np.errstate(over='ignore'):
            ex = np.exp(np.nan_to_num(cbp - m0, nan=-np.inf))
        lm = np.log(np.maximum(ex.sum(axis=(1,2,3,4,5,6,7)), 1e-300)) + m0
        ima = int(np.argmax(lm))
        am = A_GRID[ima]
        if 0 < ima < len(A_GRID)-1:
            x = A_GRID[ima-1:ima+2]; y = lm[ima-1:ima+2]
            c2_, c1_, _ = np.polyfit(x, y, 2)
            if c2_ < 0: am = -c1_/(2*c2_)
        dnm = float(lm.max() - lm[0])
        # wr posterior mass beyond the Hwang band
        wmass = ex.sum(axis=(0,1,3,4,5,6,7))
        whi = float(wmass[3:].sum()/wmass.sum())
        kmass = ex.sum(axis=(0,1,2,3,4,5,6))
        kpost = np.round(kmass/max(kmass.sum(), 1e-300), 2).tolist()
        P(f"MARG {SAMPLE} seed {seed} {law}: a_marg={am:.2f}, "
          f"dN_marg={dnm:+.1f}, P(wr>=0.4)={whi:.2f}, P(kw)={kpost}")
        # ridge + stress diagnostics (pre-verdict additions)
        mw = ex.sum(axis=(0,1,4,5,6,7))
        ww = mw/np.maximum(mw.sum(), 1e-300)
        wrv = WR_GRID[:, None]*np.ones_like(ww)
        fcv = FCOMP_GRID[None, :]*np.ones_like(ww)
        mwr = float((ww*wrv).sum()); mfc = float((ww*fcv).sum())
        cov = float((ww*(wrv-mwr)*(fcv-mfc)).sum())
        sd = np.sqrt(max(float((ww*(wrv-mwr)**2).sum()), 1e-24)
                     * max(float((ww*(fcv-mfc)**2).sum()), 1e-24))
        rho = cov/max(sd, 1e-12)
        st = {}
        for fv in (0.5, 0.7):
            fi_ = int(np.argmin(np.abs(FCOMP_GRID-fv)))
            _, a_c, _ = profile_of(cb[:, :, :, fi_:fi_+1])
            st[fv] = a_c
        P(f"RIDGE {SAMPLE} seed {seed} {law}: corr(wr,fcomp)={rho:+.2f}; "
          f"stress a_prof(fcomp=0.5)={st[0.5]:.2f}, (0.7)={st[0.7]:.2f}")
    P(f"  seed {seed} done ({(time.time()-t0)/60:.1f} min)")

for seed in SEEDS_REQ:
    have = have_rows()
    if ('simple', seed) in have and ('BE', seed) in have:
        continue
    run_seed(seed)

# --- verdict on whatever seeds now exist ----------------------------------
have = have_rows()
seeds = sorted({s for (_, s) in have})
if all(('simple', s) in have and ('BE', s) in have for s in seeds) and seeds:
    am_ = {law: np.mean([have[(law, s)][0] for s in seeds])
           for law in ('simple', 'BE')}
    dn_ = {law: np.mean([have[(law, s)][1] for s in seeds])
           for law in ('simple', 'BE')}
    amin, dmin = min(am_.values()), min(dn_.values())
    if SAMPLE == 'strictpow':
        v = ('POWER-OK' if (amin >= 0.9 and max(am_.values()) <= 1.5
                            and dmin >= 25) else
             'POWER-FAIL (the strict-companion quadrant carries no weight)')
    elif SAMPLE == 'full':
        v = ('ANCHOR-CONFIRMED' if (amin >= 0.9 and dmin >= 30) else
             'COMPANION-WIN' if (amin <= 0.7 or dmin <= 15) else 'AMBIGUOUS')
    else:
        v = ('RESOLVED-INSTRUMENTAL' if (amin >= 0.9 and dmin >= 25) else
             'COMPANION-CONFIRMED' if max(am_.values()) <= 0.7 else 'MIDDLE')
    lines = [f"", f"7J {SAMPLE} [{AMP}] verdict ({len(seeds)} seeds): "
             f"simple a_marg={am_['simple']:.2f} dN={dn_['simple']:+.1f}; "
             f"BE a_marg={am_['BE']:.2f} dN={dn_['BE']:+.1f}  ==> {v}"]
    print("\n".join(lines))
    with open('data/stage7j_verdict.txt', 'a') as f:
        f.write("\n".join(lines)+"\n")
print(f"\nSTAGE 7J-B {SAMPLE} pass complete")
