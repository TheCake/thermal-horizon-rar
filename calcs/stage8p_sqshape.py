"""STAGE 8P — THE SQ-TAIL SHAPE CONTEST (the band-flank repair candidate).

Pre-registered BEFORE any run (bars/gates/stop rules locked in NOTES; this
docstring is the operational summary).  8N localized the census BAND flood
in the WIDTH channel's lognormal tail (the no-companion channel alone
floods the band ~x20 at sq=0.2).  This stage asks the two-channel
question: can a BOUNDED width distribution keep the KINEMATIC likelihood
(which demanded sq>0 four independent ways: 3E / 6P / the -60..-116
misfit / the fpm edge) while un-flooding the band?

Instrument: a single-config kinematic block evaluator (the lnL_point
ws=0 / wcut-off / legacy path lifted VERBATIM from
calcs/stage7j_marginal.py) + the 8N census forward, both at the operative
corrected-kernel MAP (alpha, eta, wr FROZEN at the lker-cube argmax;
conditioning disclosed — no alpha re-derivation this round).  The width
draw m = exp(sq*T(g)) where T is a unit-variance DETERMINISTIC transform
of the SAME standard-normal stream draw (stream-preserving; matched
second moment in ln m):
    logn   T(g) = g                       (operative; the identity)
    clip2  T(g) = c2*clip(g, -2, 2)       (clipped at 2 sigma, rescaled)
    ulog   T(g) = sqrt(3)*(2*Phi(g) - 1)  (hard-bounded uniform in ln m)
    twopt  T(g) = sign(g)                 (two-point; extreme bounded)
    lapl   T(g) = Laplace^-1(Phi(g))      (HEAVIER tail; direction ctrl)
The sq axis is EXTENDED to [0, 0.1, 0.2, 0.3, 0.4, 0.5] for ALL shapes
symmetrically (bounded shapes may want more bulk variance once their
tails are capped; the cube identity gate uses the first four nodes).

Gates (any FAIL => STOP: no shape rows quoted; diagnose; amendment
logged pre-quote; the run file is preserved):
  G8P-0 in-evaluator identity: the logn path uses the verbatim legacy
        arrays (bit-compare probe).
  G8P-1 cube identity: evaluator logn block (sq<=0.3) vs the stored lker
        cube at the MAP (alpha,eta,wr): max|d| <= 0.05 lnL
        (keep-boundary grade; exact max reported).
  G8P-2 census identity: logn census at the shipped MAP nuisances
        reproduces stage8lb_read.txt mu to 0.05 (the G8N-0 bar).
  G8P-3 moment calibration: each T has |sample var - 1| <= 0.02 and
        |sample mean| <= 0.01 on the N draws.
Bars (locked; law-seed majority = >=3/4):
  B1 SHAPE-ARTIFACT-CONFIRMED iff some bounded shape (clip2/ulog/twopt)
     has profiled Dkin >= -2.0 vs logn AND Poisson pmf(9 | mu_band at
     its profiled cell) >= 1e-3 in >=3/4 law-seeds, with no law-seed of
     that shape at Dkin <= -5.
  B3 TAIL-DEMANDED iff EVERY bounded shape has Dkin <= -5.0 in >=3/4
     law-seeds (the kinematics DEMAND the lognormal tail => the two
     data channels want opposite tails = a named inconsistency).
  Else MIXED-CARRIED.
  B2 direction control: lapl mu_band >= logn mu_band at FIXED operative
     nuisances in >=3/4 law-seeds; else DIRECTION-FLAG.
  Cliff cross-check: fixed-nuisance |d mu_cliff|/mu_cliff <= 0.20
     expected (the 8N attribution); breach = flag, not verdict.
NO credence movement (measurement round; pre-stated).  Machine line for
the 8Q combined leg: 'WINNER: <shape>|none' printed at the end.
Output: data/stage8p_sqshape.txt
"""
import re
import numpy as np, time
from astropy.io import fits
from scipy.stats import poisson
from scipy.special import ndtr

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
# corrected velocities (4R convention, 6G formulas) — VERBATIM 7J block
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
NDATA = [int(h.sum()) for h in data_2d]
POOLS = noise_pool
NOBS, NHI = 9, 2

N = 500_000
A_GRID  = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID  = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FC0 = 0.10
FFLY_GRID = np.array([0.05, 0.10])
FPM_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW_GRID = np.array([0.7, 1.0, 1.4])
SQX = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])   # extended; cube = first 4
SEEDS = (31, 101)
SHAPES = ['logn', 'clip2', 'ulog', 'twopt', 'lapl']
BOUNDED = ['clip2', 'ulog', 'twopt']
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2

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

# unit-variance shape transforms of the SAME standard-normal draw
zg = np.linspace(-12.0, 12.0, 2_000_001)
wq = np.exp(-0.5*zg*zg); wq /= wq.sum()
cc_ = np.clip(zg, -2.0, 2.0)
C2S = 1.0/np.sqrt(float(np.sum(wq*cc_*cc_)))

def shape_maps(gs):
    u = ndtr(gs)
    t = {}
    t['logn']  = gs                       # the identity: SAME array
    t['clip2'] = C2S*np.clip(gs, -2.0, 2.0)
    t['ulog']  = np.sqrt(3.0)*(2.0*u - 1.0)
    t['twopt'] = np.sign(gs)
    ul = np.clip(np.abs(2.0*u - 1.0), 0.0, 1.0 - 1e-16)
    t['lapl']  = -np.sign(u - 0.5)*np.log1p(-ul)/np.sqrt(2.0)
    return t

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
    p['pick'] = [rng.integers(0, max(len(POOLS[bi]),1), N)
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
    # amendment 1: BIT-VERBATIM marginal expression — the slope must be
    # written 0.995-erf, not the literal 0.045 (they differ by 4.2e-17
    # and the near-parabolic integrator amplifies it to keep-boundary
    # flips worth up to ~10 lnL at floor-amplified cells; the 8L-b
    # reader lineage carried the literal; see stage8pq_diag2)
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

def eval_block(p, o, tmap):
    """lnL over (fcomp, ffly, fpm, kw, sqx) per shape; the logn path is
    the lnL_point ws=0/wcut-off legacy expression sequence VERBATIM."""
    smag, vpar, vper = project(p, o)
    s_kau = smag/1e3
    out = {nm: np.zeros((len(FCOMP_GRID), len(FFLY_GRID), len(FPM_GRID),
                         len(KW_GRID), len(SQX))) for nm in tmap}
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500 or len(noise_pool[bi]) == 0: continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = noise_pool[bi][p['pick'][bi][idx] % len(noise_pool[bi])]/4.74047
        g1_i, g2_i = p['gn1'][idx], p['gn2'][idx]
        sk_i = s_kau[idx]
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
                    vp_n = vp_a*boost + g1_i*sg0*fpm
                    vq_n = vq_a*boost + g2_i*sg0*fpm
                    vmag = np.hypot(vp_n, vq_n)
                    keep = vmag*4.74047 <= (2.978/np.sqrt(sk_i)
                                            + 2.8284*sg0*4.74047)
                    vtn = (vmag/vc)[keep]
                    gmn = np.degrees(np.arccos(np.clip(
                        np.abs(vp_n[keep])/np.maximum(vmag[keep],1e-12),
                        0, 1)))
                    for nm, tv in tmap.items():
                        tk = tv[idx][keep]
                        for si, sqv in enumerate(SQX):
                            vts = vtn*np.exp(sqv*tk)
                            h,_,_ = np.histogram2d(np.clip(vts,0.021,5.9),
                                                   gmn, bins=[VE, GE])
                            p0 = np.maximum(h/max(h.sum(),1), 1e-5)
                            p0 /= p0.sum()
                            for yi, ff in enumerate(FFLY_GRID):
                                wch = min(FC0*SC2[bi], 0.5)
                                wfl = min(ff*SC2[bi], 0.5)
                                wtot = min(wch+wfl, 0.6)
                                mixc = (wch*UNI_B[bi]
                                        + wfl*FLY_B[bi])/(wch+wfl)
                                pp = (1-wtot)*p0 + wtot*mixc
                                out[nm][fi, yi, pi, ki, si] += \
                                    np.sum(data_2d[bi]*np.log(pp))
    return out

def census_mu(p, o, fcm, kw, fpm, sq, tv):
    smag, vpar, vper = project(p, o)
    s_kau = smag/1e3
    mu, mu_hi = 0.0, 0.0
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500:
            continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = POOLS[bi][p['pick'][bi][idx] % len(POOLS[bi])]/4.74047
        cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
        mh_tot = np.zeros(len(idx))
        for k in (1, 2):
            c = p['comp'][k]
            act = c['uc'][idx] < fcm
            cvp += act*c['w'][idx]*c['wd'][idx,0]
            cvq += act*c['w'][idx]*c['wd'][idx,1]
            mh_tot += act*c['mh'][idx]
        boost = np.sqrt(1+mh_tot/p['M_s'][idx])
        vp_n = (vpar[idx]+kw*cvp)*boost + p['gn1'][idx]*sg0*fpm
        vq_n = (vper[idx]+kw*cvq)*boost + p['gn2'][idx]*sg0*fpm
        vmag = np.hypot(vp_n, vq_n)
        keep = vmag*4.74047 <= (2.978/np.sqrt(s_kau[idx])
                                + 2.8284*sg0*4.74047)
        vtn = (vmag/vc)[keep]
        gmn = np.degrees(np.arccos(np.clip(
            np.abs(vp_n[keep])/np.maximum(vmag[keep],1e-12), 0, 1)))
        vts = vtn*np.exp(sq*tv[idx][keep])
        nk = max(int(keep.sum()), 1)
        mb = (gmn >= 75) & (vts >= 1.414) & (vts < 1.67)
        mc = (gmn >= 75) & (vts >= 1.67) & (vts < 2.2)
        mu += NDATA[bi]*float(np.sum(mb))/nk
        mu_hi += NDATA[bi]*float(np.sum(mc))/nk
    return mu, mu_hi

def pmf(k, mu):
    return float(poisson.pmf(k, max(mu, 1e-12)))

ship = {}
for ln in open('data/stage8lb_read.txt').read().splitlines():
    mm = re.match(r"\[(simple|BE) (\d+)\] MAP cell: alpha=([\d.]+), "
                  r"wr=([\d.]+), fcomp=([\d.]+), fpm=([\d.]+), "
                  r"kw=([\d.]+), sq=([\d.]+); CENSUS corrected-kernel: "
                  r"mu=\(([\d.]+), ([\d.]+)\)", ln)
    if mm:
        ship[(mm.group(1), int(mm.group(2)))] = dict(
            alpha=float(mm.group(3)), wr=float(mm.group(4)),
            fcomp=float(mm.group(5)), fpm=float(mm.group(6)),
            kw=float(mm.group(7)), sq=float(mm.group(8)),
            mub=float(mm.group(9)), muh=float(mm.group(10)))
assert len(ship) == 4, ship

t0 = time.time()
P("8P THE SQ-TAIL SHAPE CONTEST (pre-reg committed BEFORE any run; "
  "measurement round; NO credence movement)")
P(f"data N per bin = {NDATA}; observed pair = (band {NOBS}, cliff {NHI}); "
  f"shapes = {SHAPES}; sq grid = {SQX.tolist()} (cube = first 4); "
  f"clip2 rescale c2 = {C2S:.6f}")
P("")
g0_ok = g1_ok = g2_ok = g3_ok = True
rows = {nm: [] for nm in SHAPES}      # (law,seed) -> dict per shape
b2_dir, cliff_xchk = [], []
for seed in SEEDS:
    pf = build_pop(seed)
    tmap = shape_maps(pf['gs'])
    ok30 = tmap['logn'] is pf['gs']
    g0_ok &= ok30
    for nm in SHAPES:
        tv = tmap[nm]
        mv, vv = float(np.mean(tv)), float(np.var(tv))
        ok3 = (abs(mv) <= 0.01) and (abs(vv-1.0) <= 0.02)
        g3_ok &= ok3
        P(f"[seed {seed}] G8P-3 {nm:<5}: mean={mv:+.4f}, var={vv:.4f} "
          f"-> {'PASS' if ok3 else 'FAIL'}")
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        c9 = np.load(f'data/stage7j_cube_full_lker_{seed}_{law}.npy')
        pri = (prior_eta.reshape((1, 2) + (1,)*7)
               + LNPI.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1)))
        bx = np.unravel_index(np.nanargmax(c9 + pri), c9.shape)
        al = A_GRID[bx[0]]; eta = E_GRID[bx[1]]; wr = WR_GRID[bx[2]]
        sh = ship[(law, seed)]
        assert all(abs(v-sh[k]) < 1e-9 for v, k in
                   ((al, 'alpha'), (wr, 'wr'))), (bx, sh)
        e_f = e_of_x(pf, eta, wr)
        tab_a = 1.0 + al*(TAB-1.0)
        vp_f = vp_c(pf, e_f, tab_a)
        o_f = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'], pf['M_s'],
                  pf['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_a,
                  lny0=LNY0, dlny=DLNY, vp=vp_f)
        blk = eval_block(pf, o_f, tmap)
        # ---- G8P-1 cube identity (logn, sq<=0.3) --------------------
        sl = np.squeeze(c9[bx[0], bx[1], bx[2]], axis=1)   # (6,2,6,3,4)
        dmax = float(np.max(np.abs(blk['logn'][..., :4] - sl)))
        ok1 = dmax <= 0.05
        g1_ok &= ok1
        P(f"[{law} {seed}] G8P-1 cube identity at MAP (alpha={al}, "
          f"eta={eta}, wr={wr}): max|d| = {dmax:.2e} "
          f"-> {'PASS' if ok1 else 'FAIL'}")
        # ---- G8P-2 census identity ----------------------------------
        mu_l, mh_l = census_mu(pf, o_f, sh['fcomp'], sh['kw'], sh['fpm'],
                               sh['sq'], tmap['logn'])
        ok2 = abs(mu_l-sh['mub']) <= 0.05 and abs(mh_l-sh['muh']) <= 0.05
        g2_ok &= ok2
        P(f"[{law} {seed}] G8P-2 census identity: mu=({mu_l:.2f}, "
          f"{mh_l:.2f}) vs shipped ({sh['mub']:.2f}, {sh['muh']:.2f}) "
          f"-> {'PASS' if ok2 else 'FAIL'}")
        # ---- the contest --------------------------------------------
        sc = {nm: blk[nm] + LNPI.reshape(6, 1, 1, 1, 1) for nm in SHAPES}
        mx0 = float(np.max(sc['logn']))
        for nm in SHAPES:
            cix = np.unravel_index(np.argmax(sc[nm]), sc[nm].shape)
            fcm_h = FCOMP_GRID[cix[0]]; ffl_h = FFLY_GRID[cix[1]]
            fpm_h = FPM_GRID[cix[2]]; kw_h = KW_GRID[cix[3]]
            sq_h = SQX[cix[4]]
            dk = float(np.max(sc[nm])) - mx0
            dk_pure = float(np.max(blk[nm])) - float(np.max(blk['logn']))
            w_ = np.exp(sc[nm] - np.max(sc[nm]))
            w_ /= w_.sum()
            pfpm3 = float(w_[:, :, FPM_GRID == 3.0, :, :].sum())
            psqx = float(w_[..., 4:].sum())
            mu_p, mh_p = census_mu(pf, o_f, fcm_h, kw_h, fpm_h, sq_h,
                                   tmap[nm])
            mu_x, mh_x = census_mu(pf, o_f, sh['fcomp'], sh['kw'],
                                   sh['fpm'], sh['sq'], tmap[nm])
            edge = []
            if cix[4] == len(SQX)-1: edge.append('sq=0.5-EDGE')
            if abs(fpm_h-3.0) < 1e-9: edge.append('fpm=3.0-EDGE')
            if abs(kw_h-0.7) < 1e-9: edge.append('kw-floor')
            if cix[0] in (0, len(FCOMP_GRID)-1): edge.append('fcomp-edge')
            P(f"[{law} {seed}] {nm:<5}: Dkin={dk:+8.2f} (pure "
              f"{dk_pure:+8.2f}); cell fcomp={fcm_h:.2f} ffly={ffl_h:.2f} "
              f"fpm={fpm_h:.1f} kw={kw_h:.1f} sq={sq_h:.1f}; "
              f"P(fpm=3.0)={pfpm3:.2f} P(sq>0.3)={psqx:.2f}; "
              f"census@prof mu=({mu_p:.2f}, {mh_p:.2f}) "
              f"pmf9={pmf(NOBS, mu_p):.1e} pmf2={pmf(NHI, mh_p):.1e}; "
              f"census@fixed mu=({mu_x:.2f}, {mh_x:.2f})"
              + (("; " + ",".join(edge)) if edge else ""))
            rows[nm].append(dict(law=law, seed=seed, dk=dk, mu_p=mu_p,
                                 mh_p=mh_p, p9=pmf(NOBS, mu_p),
                                 mu_x=mu_x, mh_x=mh_x))
            if nm == 'lapl':
                b2_dir.append(mu_x >= rows['logn'][-1]['mu_x'])
            if nm in BOUNDED:
                base = rows['logn'][-1]['mh_x']
                cliff_xchk.append(abs(mh_x-base)/max(base, 1e-9) <= 0.20)
        P("")

if not (g0_ok and g1_ok and g2_ok and g3_ok):
    P("GATES FAILED (G8P-0 " + ('PASS' if g0_ok else 'FAIL')
      + ", G8P-1 " + ('PASS' if g1_ok else 'FAIL')
      + ", G8P-2 " + ('PASS' if g2_ok else 'FAIL')
      + ", G8P-3 " + ('PASS' if g3_ok else 'FAIL')
      + ") - STOP; DO NOT QUOTE; verdict WITHHELD")
    P("WINNER: none")
else:
    P("GATES: G8P-0 2/2, G8P-1 4/4, G8P-2 4/4, G8P-3 10/10 - ALL PASS")
    P("")
    passing = []
    for nm in BOUNDED:
        r = rows[nm]
        n_ok = sum(1 for x in r if (x['dk'] >= -2.0 and x['p9'] >= 1e-3))
        n_bad = sum(1 for x in r if x['dk'] <= -5.0)
        P(f"B1 {nm:<5}: pass-rows {n_ok}/4 (Dkin>=-2 AND pmf9>=1e-3); "
          f"kill-rows {n_bad}/4 (Dkin<=-5)")
        if n_ok >= 3 and n_bad == 0:
            passing.append((float(np.mean([x['dk'] for x in r])), nm))
    winner = max(passing)[1] if passing else None
    tail_dem = all(sum(1 for x in rows[nm] if x['dk'] <= -5.0) >= 3
                   for nm in BOUNDED)
    b2ok = sum(b2_dir) >= 3
    cxok = sum(cliff_xchk)
    P(f"B2 direction (lapl band >= logn band at fixed nuisances): "
      f"{sum(b2_dir)}/4 -> {'OK' if b2ok else 'DIRECTION-FLAG'}")
    P(f"Cliff cross-check (bounded |d mu_cliff| <= 20% at fixed): "
      f"{cxok}/{len(cliff_xchk)}"
      + ("" if cxok == len(cliff_xchk) else " -> FLAG (8N attribution "
         "strained)"))
    if winner is not None:
        P(f"==> 8P VERDICT (locked grammar): SHAPE-ARTIFACT-CONFIRMED — "
          f"the bounded shape '{winner}' keeps the kinematic likelihood "
          f"and un-floods the band; the census band flank is attributable "
          f"to the LOGNORMAL TAIL CONVENTION of the width channel, not "
          f"to the width channel's existence. Successor: cube-grade "
          f"re-run of '{winner}' (own pre-reg) = operative-model "
          f"candidate.")
        P(f"WINNER: {winner}")
    elif tail_dem:
        P("==> 8P VERDICT (locked grammar): TAIL-DEMANDED — every "
          "bounded shape loses >= 5 lnL in >= 3/4 law-seeds: the "
          "kinematic likelihood DEMANDS the lognormal's far tail while "
          "the census rejects it — the two data channels want OPPOSITE "
          "tails of the same object = a named model inconsistency for "
          "the population successor.")
        P("WINNER: none")
    else:
        P("==> 8P VERDICT (locked grammar): MIXED-CARRIED — neither the "
          "repair bar nor the tail-demand bar fires; per-shape rows "
          "stand as measurements.")
        P("WINNER: none")
    P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage8p_sqshape.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage8p_sqshape.txt")
