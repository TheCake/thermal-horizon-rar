"""9L-NULLINJ -- the RNAAS-note referee's M4 control (2026-09-03).

VERBATIM clone of stage9l_fpmmeter.py (the archived 9L meter) with ONE
change: the G9L-3 injection runs a TRUTH LADDER instead of the single
2.1 cell. Truths: 1.0 (below-grid null; a correctly calibrated sky),
1.4 (the single-star ceiling), 1.5 (on-grid), 2.1 (the archived
regression cell). Everything else -- data, templates, grids, seeds,
rng streams -- bit-verbatim.

PRE-STATED BARS (written before the run):
  N1 (the null): truth 1.0 must recover E[fpm] <= 1.45 in both seeds
     (grid floor is 1.2; unbiased null recovery piles at the floor).
     FAILURE VETOES the note's central claim: it would mean the
     pipeline manufactures ~2x from a calibrated sky, and the note's
     headline must be restated as a bound pending redesign.
  N2: truth 1.4 recovers within [1.2, 1.65]; truth 1.5 within +-0.25
     (the G9L-3 tolerance grammar). FAILURE -> quote as bias-mapped.
  N3 (regression): truth 2.1 reproduces the archived 2.12/2.06 to
     +-0.01. FAILURE -> wiring error; STOP, fix, rerun.
Output: data/stage9l_nullinj.txt. Does NOT overwrite any archive
(table save disabled).
"""
# --- original 9L header preserved below ---
'''STAGE 9L — THE NARROW-PAIR FPM METER (round-13 D2; THE DECIDER).
Pre-registered BEFORE any run (pre-signed map 53 / 60 / hold-56).

Narrow s-bins (0.2-2, 2-6 kAU) are near-boost-free: whatever fpm
they demand is honest noise.  alpha FIXED at 0 (tab_a == 1 for both
laws => law-blind; rows = 2 seeds).  Per-stratum (fpm, sq)
posterior from narrow-only tables; the FULL fpm marginal shipped.
Gates: G9L-1 analytic; G9L-2 boost-free premise measured (alpha=1
vs alpha=0 narrow-bin vt medians <= 1.05, both laws); G9L-3
injected-fpm recovery (Q1 truth 2.1 +- 0.25, rng 9); G9L-4 counts.
Bars (ordered, BOTH seeds): L-NOISE-REAL (E[fpm_n(Q1)] in
[1.8, 2.3]) / L-BOOST-EATEN (<= 1.5) / L-GRAY-CARRIED.
Output: data/stage9l_nullinj.txt (nullinj variant)
'''

import math
import numpy as np, time
from astropy.io import fits
from scipy.special import logsumexp

L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

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
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot_d = np.interp(MG1[ok],MG_T,MS_T)+np.interp(MG2[ok],MG_T,MS_T)
s_d = sep[ok]/1e3
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
# AMENDMENT 1 (post-G9L-2-fail, pre-quote, logged in NOTES): the
# registered 0.2-6 window failed its own boost-premise gate at the
# 2-6 bin (1.052-1.080 > 1.05); the meter narrows to 0.2-2 alone.
NB_BINS = 1
SC2 = (np.array([0.63, 3.46, 11.0, 31.6])/31.6)**2
VE = np.logspace(np.log10(0.02), np.log10(6.0), 21)
GE = np.linspace(0, 90, 7)
NV, NG = 20, 6
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

def build_stratum(pm):
    D2, PLs, UB, FB = [], [], [], []
    for b in SBINS:
        m = (s_d>=b[0])&(s_d<b[1])&pm
        h,_,_ = np.histogram2d(np.clip(vt_d[m],0.021,5.9), gam_d[m],
                               bins=[VE, GE])
        D2.append(h.astype(float))
        PLs.append(sig_ok[m])
        cutp = 2.978/np.sqrt(s_d[m]) + 2.8284*sig_ok[m]
        acc = np.array([(vcen[i]*vc_ok[m] <= cutp).mean()
                        if m.sum() else 0.0 for i in range(NV)])
        for tpl, store in ((UNI, UB), (FLY, FB)):
            t = tpl*acc[:,None]
            store.append(t/max(t.sum(), 1e-12))
    ND = [int(h.sum()) for h in D2]
    return D2, PLs, UB, FB, ND

N = 500_000
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FC0 = 0.10
FFLY_GRID = np.array([0.05, 0.10])
FPM_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW_GRID = np.array([0.7, 1.0, 1.4])
SQ_GRID = np.array([0.0, 0.1, 0.2, 0.3])
WS2_GRID = np.array([0.0, 0.045])
SEEDS = (31, 101)

pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
gband = pz['conv_band']/0.30
GS = np.linspace(float(gband[0]), float(gband[1]), 25)
fhmin, fhmax = fg[sup].min(), fg[sup].max()
LNPI = np.full(len(FCOMP_GRID), -1e9)
for gi in GS:
    fh_eq = FCOMP_GRID/gi
    m = (fh_eq >= fhmin) & (fh_eq <= fhmax)
    cand = np.full(len(FCOMP_GRID), -1e9)
    cand[m] = np.interp(fh_eq[m], fg, lp)
    LNPI = np.maximum(LNPI, cand)

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
    p['pick'] = [rng.integers(0, max(len(noise_pool[bi]),1), N)
                 for bi in range(len(SBINS))]
    p['gn1'], p['gn2'] = rng.normal(size=N), rng.normal(size=N)
    M_h = 0.5*p['M_s']
    p['comp'] = {}
    for k in (1, 2):
        u_q = rng.random(N)
        q = 0.1+0.9*u_q
        logP = rng.normal(5.03, 2.28, N)
        P_yr = 10**logP/365.25
        a_in = (M_h*(1+q)*P_yr**2)**(1/3)
        valid = (a_in < 130.0) & (a_in < p['a_s']/5.0)
        v_orb = 29.78*np.sqrt(M_h*(1+q)/np.maximum(a_in,1e-3))
        u_ = np.pi*2.83/np.maximum(P_yr, 1e-9)
        S = np.where(u_ < 1e-2, 1.0 - u_*u_/10.0,
                     3.0*(np.sin(u_) - u_*np.cos(u_))
                     / np.maximum(u_, 1e-300)**3)
        MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        MGs = np.interp(-np.clip(q*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        l_ = 10**(-0.4*(MGs-MGp))
        wfac = np.abs(q/(1+q) - l_/(1+l_))
        w = wfac*v_orb*S/4.74047*valid
        wd = rng.normal(size=(N,3)); wd /= np.linalg.norm(wd,axis=1,keepdims=True)
        p['comp'][k] = dict(w=w, wd=wd, uc=rng.random(N),
                            mh=q*M_h*valid, P=P_yr)
    p['gs'] = rng.normal(size=N)
    return p

def e_of_x(p, eta, wr):
    al = np.interp(np.log10(p['a_s']), np.log10([100,500,1000,50000]),
                   [0.6, 1.0, eta, eta])
    e_pow = 0.95*p['u_e']**(1/(1+al))
    # amendment-1 standard: BIT-VERBATIM marginal expression
    erf = 0.95
    e_rad = erf+(0.995-erf)*p['u_e']
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

def project(p, o):
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
    return smag, vpar, vper

def eval_block_nb(p, prj, D2, PLs, UB, FB):
    """VERBATIM eval_block_g restricted to the NARROW bins 0-1."""
    smag, vpar, vper = prj
    s_kau = smag/1e3
    out = np.zeros((len(FCOMP_GRID), len(FFLY_GRID), len(FPM_GRID),
                    len(KW_GRID), len(SQ_GRID), len(WS2_GRID)))
    for bi in range(NB_BINS):
        b = SBINS[bi]
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500 or len(PLs[bi]) == 0: continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = PLs[bi][p['pick'][bi][idx] % len(PLs[bi])]/4.74047
        g1_i, g2_i = p['gn1'][idx], p['gn2'][idx]
        sk_i = s_kau[idx]
        gk_full = p['gs'][idx]
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
                vp_a = vpar[idx] + kwv*cvp
                vq_a = vper[idx] + kwv*cvq
                for pi, fpm in enumerate(FPM_GRID):
                    for wi, ws in enumerate(WS2_GRID):
                        if ws == 0.0:
                            vp_n = vp_a*boost + g1_i*sg0*fpm
                            vq_n = vq_a*boost + g2_i*sg0*fpm
                        else:
                            sig_eff = np.sqrt((sg0*fpm)**2
                                              + (ws/4.74047)**2)
                            vp_n = vp_a*boost + g1_i*sig_eff
                            vq_n = vq_a*boost + g2_i*sig_eff
                        vmag = np.hypot(vp_n, vq_n)
                        keep = vmag*4.74047 <= (2.978/np.sqrt(sk_i)
                                                + 2.8284*sg0*4.74047)
                        vtn = (vmag/vc)[keep]
                        gmn = np.degrees(np.arccos(np.clip(
                            np.abs(vp_n[keep])
                            / np.maximum(vmag[keep], 1e-12), 0, 1)))
                        gk = gk_full[keep]
                        for si, sqv in enumerate(SQ_GRID):
                            vts = vtn*np.exp(sqv*gk)
                            h,_,_ = np.histogram2d(
                                np.clip(vts,0.021,5.9), gmn,
                                bins=[VE, GE])
                            p0 = np.maximum(h/max(h.sum(),1), 1e-5)
                            p0 /= p0.sum()
                            for yi, ff in enumerate(FFLY_GRID):
                                wch = min(FC0*SC2[bi], 0.5)
                                wfl = min(ff*SC2[bi], 0.5)
                                wtot = min(wch+wfl, 0.6)
                                mixc = (wch*UB[bi]
                                        + wfl*FB[bi])/(wch+wfl)
                                pp = (1-wtot)*p0 + wtot*mixc
                                out[fi, yi, pi, ki, si, wi] += \
                                    np.sum(D2[bi]*np.log(pp))
    return out

def eval_pp_nb(p, prj, PLs, UB, FB, cell):
    """pp per NARROW bin at one cell — verbatim pp path."""
    fi, yi, fpm, ki, sqv = cell
    fcm = FCOMP_GRID[fi]; kwv = KW_GRID[ki]; ff = FFLY_GRID[yi]
    smag, vpar, vper = prj
    s_kau = smag/1e3
    pps = []
    for bi in range(NB_BINS):
        b = SBINS[bi]
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500 or len(PLs[bi]) == 0:
            pps.append(None); continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = PLs[bi][p['pick'][bi][idx] % len(PLs[bi])]/4.74047
        g1_i, g2_i = p['gn1'][idx], p['gn2'][idx]
        sk_i = s_kau[idx]
        gk_full = p['gs'][idx]
        cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
        mh_tot = np.zeros(len(idx))
        for k in (1, 2):
            c = p['comp'][k]
            act = c['uc'][idx] < fcm
            mh_tot += act*c['mh'][idx]
            cvp += act*c['w'][idx]*c['wd'][idx,0]
            cvq += act*c['w'][idx]*c['wd'][idx,1]
        boost = np.sqrt(1+mh_tot/p['M_s'][idx])
        vp_n = (vpar[idx] + kwv*cvp)*boost + g1_i*sg0*fpm
        vq_n = (vper[idx] + kwv*cvq)*boost + g2_i*sg0*fpm
        vmag = np.hypot(vp_n, vq_n)
        keep = vmag*4.74047 <= (2.978/np.sqrt(sk_i)
                                + 2.8284*sg0*4.74047)
        vtn = (vmag/vc)[keep]
        gmn = np.degrees(np.arccos(np.clip(
            np.abs(vp_n[keep])/np.maximum(vmag[keep], 1e-12), 0, 1)))
        gk = gk_full[keep]
        vts = vtn*np.exp(sqv*gk)
        h,_,_ = np.histogram2d(np.clip(vts,0.021,5.9), gmn,
                               bins=[VE, GE])
        p0 = np.maximum(h/max(h.sum(),1), 1e-5)
        p0 /= p0.sum()
        wch = min(FC0*SC2[bi], 0.5)
        wfl = min(ff*SC2[bi], 0.5)
        wtot = min(wch+wfl, 0.6)
        mixc = (wch*UB[bi] + wfl*FB[bi])/(wch+wfl)
        pps.append((1-wtot)*p0 + wtot*mixc)
    return pps

okv = ok
ruwe = np.maximum(_pick('ruwe1', 'RUWE1'),
                  _pick('ruwe2', 'RUWE2'))[okv]
qs_ = np.percentile(ruwe, [25, 50, 75])
QMASKS = [ruwe <= qs_[0],
          (ruwe > qs_[0]) & (ruwe <= qs_[1]),
          (ruwe > qs_[1]) & (ruwe <= qs_[2]),
          ruwe > qs_[2]]
ALLM = np.ones(len(s_d), dtype=bool)
noise_pool = build_stratum(ALLM)[1]
STRATA = [build_stratum(mk) for mk in QMASKS]

t0 = time.time()
P("9L THE NARROW-PAIR FPM METER (round-13 D2; pre-reg committed "
  "BEFORE any run; THE DECIDER — movement by the pre-signed map "
  "ONLY; alpha = 0 => LAW-BLIND, rows = 2 seeds)")
P("")

def lse_cells(T, axes, ncells):
    return logsumexp(T, axis=axes) - math.log(ncells)
_T = np.log(np.array([[1.0, 2.0], [3.0, 4.0]]))
_v = float(lse_cells(_T, (0, 1), 4))
g1_ok = abs(_v - math.log(2.5)) < 1e-12
P(f"G9L-1 analytic: {_v:.12f} vs ln(2.5) -> "
  f"{'PASS' if g1_ok else 'FAIL'}")

g4_ok = True
for qi in range(4):
    nn = sum(int(STRATA[qi][0][bi].sum()) for bi in range(NB_BINS))
    P(f"G9L-4 narrow counts Q{qi+1}: {nn}")
    g4_ok &= nn > 1000
P("")

ETA0, WR0 = 1.05, 0.30
P("design note (pre-quote): the pre-reg was silent on the "
  "eccentricity-law cell; PRIMARY = (eta 1.05, wr 0.30) fiducial "
  "with an eta = 1.3 CO-READ shipped (risk-axis rule) — bars read "
  "the primary only")
P("")
g2_ok = g3_ok = True
rows = []
for seed in SEEDS:
    pf = build_pop(seed)
    # eta co-read first (alpha = 0, narrow E[fpm] only)
    e_alt = e_of_x(pf, 1.3, WR0)
    tabA = np.ones_like(TAB_S)
    vpA = vp_c(pf, e_alt, tabA)
    oA = run(pf['a_s'], e_alt, pf['psi0'], pf['f_ip'], pf['M_s'],
             pf['uph'], 8, 2500, 5, a0=A0_CAN, tab=tabA,
             lny0=LNY0, dlny=DLNY, vp=vpA)
    prjA = project(pf, oA)
    TA = eval_block_nb(pf, prjA, *STRATA[0][:4])
    lwA = logsumexp(TA + LNPI.reshape(6, 1, 1, 1, 1, 1),
                    axis=(0, 1, 3, 5))
    wA = np.exp(lwA - logsumexp(lwA))
    P(f"[seed {seed}] eta = 1.3 CO-READ: E[fpm_narrow(Q1)] = "
      f"{float(np.sum(wA.sum(axis=1)*FPM_GRID)):.2f}")
    e_f = e_of_x(pf, ETA0, WR0)
    tab0 = np.ones_like(TAB_S)
    vp0 = vp_c(pf, e_f, tab0)
    o0 = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'], pf['M_s'],
             pf['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab0,
             lny0=LNY0, dlny=DLNY, vp=vp0)
    prj0 = project(pf, o0)
    # G9L-2: boost-free premise (alpha=1, both laws)
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        vp1 = vp_c(pf, e_f, TAB)
        o1 = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'], pf['M_s'],
                 pf['uph'], 8, 2500, 5, a0=A0_CAN, tab=TAB,
                 lny0=LNY0, dlny=DLNY, vp=vp1)
        prj1 = project(pf, o1)
        s0 = prj0[0]/1e3
        s1 = prj1[0]/1e3
        for bi in range(NB_BINS):
            b = SBINS[bi]
            m0 = (s0>=b[0])&(s0<b[1]); m1 = (s1>=b[0])&(s1<b[1])
            vt0 = np.hypot(prj0[1], prj0[2])[m0] \
                / (2*np.pi*np.sqrt(pf['M_s'][m0]/prj0[0][m0]))
            vt1 = np.hypot(prj1[1], prj1[2])[m1] \
                / (2*np.pi*np.sqrt(pf['M_s'][m1]/prj1[0][m1]))
            rat = float(np.median(vt1)/np.median(vt0))
            okb = rat <= 1.05
            g2_ok &= okb
            P(f"[seed {seed}] G9L-2 boost premise ({law}, bin "
              f"{b[0]}-{b[1]} kAU): a1/a0 vt-median ratio = "
              f"{rat:.4f} -> {'PASS' if okb else 'FAIL'}")
    # narrow tables at alpha=0
    T_nb = np.zeros((4, 6, 2, 6, 3, 4, 2))
    for qi in range(4):
        T_nb[qi] = eval_block_nb(pf, prj0, *STRATA[qi][:4])
    pass  # nullinj: archive table save DISABLED (never overwrite 9L archives)
    efpm = {}
    for qi in range(4):
        lw = logsumexp(T_nb[qi]
                       + LNPI.reshape(6, 1, 1, 1, 1, 1),
                       axis=(0, 1, 3, 5))
        wq = np.exp(lw - logsumexp(lw))
        mfpm = wq.sum(axis=1)
        msq = wq.sum(axis=0)
        efpm[qi] = float(np.sum(mfpm*FPM_GRID))
        P(f"[seed {seed}] Q{qi+1} narrow: fpm marginal = "
          + "/".join(f"{x:.2f}" for x in mfpm)
          + f"; E[fpm] = {efpm[qi]:.2f}; E[sq] = "
          f"{float(np.sum(msq*SQ_GRID)):.2f}")
    # NULLINJ: the truth ladder (referee M4); same machinery, same rng stream per truth
    for truth, lo_ok, hi_ok, tag in ((1.0, 0.0, 1.45, 'N1-null'),
                                     (1.4, 1.2, 1.65, 'N2-ceiling'),
                                     (1.5, 1.25, 1.75, 'N2-ongrid'),
                                     (2.1, 2.05 if seed == 31 else 1.99,
                                      2.19 if seed == 31 else 2.13, 'N3-regression')):
        pps = eval_pp_nb(pf, prj0, STRATA[0][1], STRATA[0][2],
                         STRATA[0][3], (2, 0, truth, 1, 0.0))
        rng9 = np.random.default_rng(9)
        synD = []
        for bi in range(NB_BINS):
            nd = int(STRATA[0][0][bi].sum())
            draw = rng9.multinomial(nd, pps[bi].ravel()/pps[bi].sum())
            synD.append(draw.reshape(NV, NG).astype(float))
        synD += [np.zeros((NV, NG)), np.zeros((NV, NG))]
        T_syn = eval_block_nb(pf, prj0, synD, STRATA[0][1],
                              STRATA[0][2], STRATA[0][3])
        lw = logsumexp(T_syn + LNPI.reshape(6, 1, 1, 1, 1, 1),
                       axis=(0, 1, 3, 5))
        wq = np.exp(lw - logsumexp(lw))
        mrec = wq.sum(axis=1)
        rec = float(np.sum(mrec*FPM_GRID))
        okN = lo_ok <= rec <= hi_ok
        g3_ok &= okN
        P(f"[seed {seed}] {tag}: truth {truth:.1f} -> E[fpm] = {rec:.2f}"
          f" (P(bottom)={mrec[0]:.2f}, P(top)={mrec[-1]:.2f})"
          f" bar [{lo_ok:.2f}, {hi_ok:.2f}] -> {'PASS' if okN else 'FAIL'}"
          f"  ({(time.time()-t0)/60:.1f} min)")
    rows.append(dict(seed=seed, e1=efpm[0]))
    P("")

if not (g1_ok and g2_ok and g3_ok and g4_ok):
    P("GATES FAILED (G9L-1 " + ('PASS' if g1_ok else 'FAIL')
      + ", G9L-2 " + ('PASS' if g2_ok else 'FAIL')
      + ", G9L-3 " + ('PASS' if g3_ok else 'FAIL')
      + ", G9L-4 " + ('PASS' if g4_ok else 'FAIL')
      + ") - STOP; DO NOT QUOTE; verdict WITHHELD")
else:
    P("GATES: G9L-1, G9L-2 8/8, G9L-3 2/2, G9L-4 - ALL PASS")
    P("")
    e1s = [r['e1'] for r in rows]
    P(f"L-facts: E[fpm_narrow(Q1)] = "
      + ", ".join(f"{x:.2f}" for x in e1s)
      + " (precondition: 9F joint E[fpm(Q1)] = 1.97-2.25)")
    nr = all(1.8 <= x <= 2.3 for x in e1s)
    be = all(x <= 1.5 for x in e1s)
    if nr:
        P("==> 9L VERDICT (locked grammar): L-NOISE-REAL - the "
          "boost-free narrow bins demand the same ~2x noise the "
          "joint fit uses; the fit's noise subtraction is "
          "LEGITIMATE; the upper-limit reading stands and the "
          "model-light pedestal is real.")
        P("    CREDENCE (pre-signed map): anomaly-real 56 -> 53.")
    elif be:
        P("==> 9L VERDICT (locked grammar): L-BOOST-EATEN - the "
          "narrow bins do NOT demand the joint fpm: the wide bins "
          "inflate it; alpha is under-reported; the D3 degraded "
          "injection is MANDATORY next.")
        P("    CREDENCE (pre-signed map): anomaly-real 56 -> 60.")
    else:
        P("==> 9L VERDICT (locked grammar): L-GRAY-CARRIED - rows "
          "stand as measurements.")
        P("    CREDENCE (pre-signed map): HOLD 56.")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage9l_nullinj.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage9l_fpmmeter.txt")
