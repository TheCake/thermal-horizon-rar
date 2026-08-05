"""STAGE 8Q — THE SUBSYSTEM P-PRIOR BRACKET (the cliff-flank repair).

Pre-registered BEFORE any run (bars/gates/stop rules locked in NOTES).
8N attributed the census CLIFF flood to long-period companion wobble
spikes (share 0.77-0.81; RUWE-silent locus 0.70 at P >= 10 yr) drawn
from the FIELD-binary Raghavan lognormal logP[d] ~ N(5.03, 2.28).
Tokovinin 2014 (AJ 147:87, Table 2 — scout-primary read; full primary
read owed before paper use) MEASURES the inner-subsystem law for
components of wider pairs: L11 logP[d] ~ N(3.25 +- 0.12, 1.80 +- 0.09),
ML-corrected for detection.  Caveat carried: his sample is not
stratified by outer separation (the shortness is stability-truncation-
dominated for close outers; our kAU outers are unconstrained by
stability), so this stage runs a BRACKET, not a single swap:
    raghavan  (5.03, 2.28)   the operative identity
    tokL11    (3.25, 1.80)   the measured subsystem law
    mid       (4.14, 2.04)   representative midpoint (labeled bracket)
Stream-preserving recast: logP_new = x0 + sg*z with z the SAME
standard-normal draw; the raghavan branch passes the UNTOUCHED legacy
array (bit-exact identity).  Only deterministic values change — every
other stream draw is identical across priors, so the orbit set is
shared (G8Q-4 asserts it).

Instruments (both at the operative lker MAP; alpha/eta/wr FROZEN,
conditioning disclosed): the single-config kinematic block evaluator
(lnL_point legacy path, as 8P) over (fcomp, ffly, fpm, kw, sq<=0.5) +
the 8N census forward with P-locus decomposition.  The LNPI host prior
is kept at the operative convention; its fcomp mapping is
P-prior-dependent in MEANING (detectability mix shifts) — pure-lnL
profiled numbers are co-printed and the caveat is carried.

Gates (any FAIL => STOP; amendment pre-quote; run file preserved):
  G8Q-0 cube identity: raghavan block (sq<=0.3) vs the stored lker cube
        at MAP — max|d| <= 0.05 lnL.
  G8Q-1 census identity: raghavan census at shipped MAP nuisances vs
        stage8lb_read.txt mu to 0.05.
  G8Q-2 recast calibration: per prior, |mean(logP) - x0| <= 0.02 and
        |std(logP) - sg| <= 0.02 on the N draws (per companion slot).
  G8Q-4 shared-orbit wiring: non-companion arrays bit-identical across
        priors (a_s probe).
  G8Q-5 locus completeness: P-locus rows sum to mu (1e-9).
Bars (locked; law-seed majority = >=3/4):
  Q1 CLIFF-REPAIRED iff mu_cliff(tokL11 @ its profiled cell) <=
     0.5 * mu_cliff(raghavan @ its profiled cell) AND pmf(2 | mu_cliff)
     >= 1e-3, in >=3/4 law-seeds.
  Q2 kinematics: tokL11 ACCEPTED iff Dkin >= -2.0 in >=3/4; REJECTED
     iff Dkin <= -5.0 in >=3/4; else CARRIED.
  Q3 (secondary, pre-stated): the noise/wobble edge fingerprints relax
     under tokL11 — P(fpm=3.0) < 0.5 or kw off the 0.7 floor in the
     profiled block (informational, not a verdict axis).
  Verdict = Q1-status x Q2-status; mid row reported as the bracket.
NO credence movement (measurement round; pre-stated).

COMBINED LEG (own invocation: `py calcs/stage8q_pprior.py combined`,
run only after 8P's record exists; rule pre-stated): width shape = 8P's
'WINNER:' line (none -> logn); population prior = tokL11.  Per
law-seed: profile the combined block, census at the profiled cell,
jointP = pmf(9|mu_band)*pmf(2|mu_cliff).  Bar: JOINT-COHERENT iff
jointP >= 1e-3 AND Dkin(combined vs operative raghavan-logn) >= -5.0
in >=3/4 law-seeds; else JOINT-NOT-REACHED.  This re-asks 8M's
admissibility question at the repaired model, POINT-grade (cube-grade
re-run = the named successor if coherent).
Outputs: data/stage8q_pprior.txt / data/stage8q_combined.txt
"""
import re, sys
import numpy as np, time
from astropy.io import fits
from scipy.stats import poisson
from scipy.special import ndtr

MODE = sys.argv[1] if len(sys.argv) > 1 else 'main'
assert MODE in ('main', 'combined'), MODE

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
SQX = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
SEEDS = (31, 101)
PPRIORS = {'raghavan': (5.03, 2.28), 'tokL11': (3.25, 1.80),
           'mid': (4.14, 2.04)}
PORDER = ['raghavan', 'tokL11', 'mid']
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
PBINS = [3.0, 10.0, 30.0, 100.0]
PLBL = ['P<3', '3-10', '10-30', '30-100', 'P>=100', 'nocomp']

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

zg = np.linspace(-12.0, 12.0, 2_000_001)
wq = np.exp(-0.5*zg*zg); wq /= wq.sum()
cc_ = np.clip(zg, -2.0, 2.0)
C2S = 1.0/np.sqrt(float(np.sum(wq*cc_*cc_)))

def shape_maps(gs):
    u = ndtr(gs)
    t = {}
    t['logn']  = gs
    t['clip2'] = C2S*np.clip(gs, -2.0, 2.0)
    t['ulog']  = np.sqrt(3.0)*(2.0*u - 1.0)
    t['twopt'] = np.sign(gs)
    ul = np.clip(np.abs(2.0*u - 1.0), 0.0, 1.0 - 1e-16)
    t['lapl']  = -np.sign(u - 0.5)*np.log1p(-ul)/np.sqrt(2.0)
    return t

def build_pop(seed, pp):
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
    p['_pdiag'] = {}
    for k in (1, 2):
        u_q = rng.random(N)
        q = 0.1+0.9*u_q
        logP = rng.normal(5.03, 2.28, N)
        if pp != 'raghavan':
            # stream-preserving recast: same z-score, new lognormal;
            # the raghavan branch passes the UNTOUCHED legacy array
            x0, sg = PPRIORS[pp]
            logP = x0 + sg*((logP - 5.03)/2.28)
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
        p['_pdiag'][k] = (float(np.mean(logP)), float(np.std(logP)),
                          float(np.mean(valid)),
                          float(np.mean(valid & (P_yr >= 10.0))))
    p['gs'] = rng.normal(size=N)
    return p

def e_of_x(p, eta, wr):
    al = np.interp(np.log10(p['a_s']), np.log10([100,500,1000,50000]),
                   [0.6, 1.0, eta, eta])
    e_pow = 0.95*p['u_e']**(1/(1+al))
    e_rad = 0.95+0.045*p['u_e']
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

def census_locus(p, o, fcm, kw, fpm, sq, tv):
    smag, vpar, vper = project(p, o)
    s_kau = smag/1e3
    mu, mu_hi = 0.0, 0.0
    lb = np.zeros(len(PLBL)); lc = np.zeros(len(PLBL))
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500:
            continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = POOLS[bi][p['pick'][bi][idx] % len(POOLS[bi])]/4.74047
        cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
        mh_tot = np.zeros(len(idx))
        w1 = np.zeros(len(idx)); w2 = np.zeros(len(idx))
        for k in (1, 2):
            c = p['comp'][k]
            act = c['uc'][idx] < fcm
            cvp += act*c['w'][idx]*c['wd'][idx,0]
            cvq += act*c['w'][idx]*c['wd'][idx,1]
            mh_tot += act*c['mh'][idx]
            aw = act*np.abs(c['w'][idx])
            if k == 1:
                w1 = aw
            else:
                w2 = aw
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
        hasact = ((w1 > 0) | (w2 > 0))[keep]
        Pdom = np.where(w1 >= w2, p['comp'][1]['P'][idx],
                        p['comp'][2]['P'][idx])[keep]
        cls = np.digitize(Pdom, PBINS)
        cls = np.where(hasact, cls, len(PLBL)-1)
        for ci in range(len(PLBL)):
            sel = cls == ci
            lb[ci] += NDATA[bi]*float(np.sum(mb & sel))/nk
            lc[ci] += NDATA[bi]*float(np.sum(mc & sel))/nk
    return mu, mu_hi, lb, lc

def jointp(mu_b, mu_c):
    return float(poisson.pmf(NOBS, max(mu_b, 1e-12))
                 * poisson.pmf(NHI, max(mu_c, 1e-12)))

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
if MODE == 'main':
    P("8Q THE SUBSYSTEM P-PRIOR BRACKET (pre-reg committed BEFORE any "
      "run; measurement round; NO credence movement)")
    P(f"priors: raghavan (5.03, 2.28) identity | tokL11 (3.25, 1.80) "
      f"Tokovinin-2014 Table-2 L11 (scout-primary) | mid (4.14, 2.04) "
      f"representative bracket; data N per bin = {NDATA}; observed pair "
      f"= (band {NOBS}, cliff {NHI})")
    P("")
    g0_ok = g1_ok = g2_ok = g4_ok = g5_ok = True
    R = {pp: [] for pp in PORDER}
    for seed in SEEDS:
        pops = {pp: build_pop(seed, pp) for pp in PORDER}
        ok4 = all(np.array_equal(pops[pp]['a_s'], pops['raghavan']['a_s'])
                  and np.array_equal(pops[pp]['gs'], pops['raghavan']['gs'])
                  for pp in PORDER)
        g4_ok &= ok4
        P(f"[seed {seed}] G8Q-4 shared non-companion streams: "
          f"{'PASS' if ok4 else 'FAIL'}")
        for pp in PORDER:
            x0, sg = PPRIORS[pp]
            for k in (1, 2):
                mv, sv, fv, f10 = pops[pp]['_pdiag'][k]
                ok2 = (abs(mv-x0) <= 0.02) and (abs(sv-sg) <= 0.02)
                g2_ok &= ok2
                if k == 1:
                    P(f"[seed {seed}] G8Q-2 {pp:<8}: mean={mv:.3f} "
                      f"(x0={x0}), std={sv:.3f} (sg={sg}) "
                      f"-> {'PASS' if ok2 else 'FAIL'}; P(valid)={fv:.3f}, "
                      f"P(valid & P>=10yr)={f10:.3f}")
        for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
            c9 = np.load(f'data/stage7j_cube_full_lker_{seed}_{law}.npy')
            pri = (prior_eta.reshape((1, 2) + (1,)*7)
                   + LNPI.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1)))
            bx = np.unravel_index(np.nanargmax(c9 + pri), c9.shape)
            al = A_GRID[bx[0]]; eta = E_GRID[bx[1]]; wr = WR_GRID[bx[2]]
            sh = ship[(law, seed)]
            assert all(abs(v-sh[k_]) < 1e-9 for v, k_ in
                       ((al, 'alpha'), (wr, 'wr'))), (bx, sh)
            pf = pops['raghavan']
            e_f = e_of_x(pf, eta, wr)
            tab_a = 1.0 + al*(TAB-1.0)
            vp_f = vp_c(pf, e_f, tab_a)
            o_f = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'], pf['M_s'],
                      pf['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_a,
                      lny0=LNY0, dlny=DLNY, vp=vp_f)
            base_mx = None
            for pp in PORDER:
                pq = pops[pp]
                blk = eval_block(pq, o_f, {'logn': pq['gs']})['logn']
                if pp == 'raghavan':
                    sl = np.squeeze(c9[bx[0], bx[1], bx[2]], axis=1)
                    dmax = float(np.max(np.abs(blk[..., :4] - sl)))
                    ok0 = dmax <= 0.05
                    g0_ok &= ok0
                    P(f"[{law} {seed}] G8Q-0 cube identity at MAP "
                      f"(alpha={al}, eta={eta}, wr={wr}): max|d| = "
                      f"{dmax:.2e} -> {'PASS' if ok0 else 'FAIL'}")
                    mu_l, mh_l, lbl_, lcl_ = census_locus(
                        pq, o_f, sh['fcomp'], sh['kw'], sh['fpm'],
                        sh['sq'], pq['gs'])
                    ok1 = (abs(mu_l-sh['mub']) <= 0.05
                           and abs(mh_l-sh['muh']) <= 0.05)
                    g1_ok &= ok1
                    P(f"[{law} {seed}] G8Q-1 census identity: mu="
                      f"({mu_l:.2f}, {mh_l:.2f}) vs shipped "
                      f"({sh['mub']:.2f}, {sh['muh']:.2f}) "
                      f"-> {'PASS' if ok1 else 'FAIL'}")
                sc = blk + LNPI.reshape(6, 1, 1, 1, 1)
                mx = float(np.max(sc))
                if pp == 'raghavan':
                    base_mx = mx
                    base_pure = float(np.max(blk))
                cix = np.unravel_index(np.argmax(sc), sc.shape)
                fcm_h = FCOMP_GRID[cix[0]]; fpm_h = FPM_GRID[cix[2]]
                kw_h = KW_GRID[cix[3]]; sq_h = SQX[cix[4]]
                dk = mx - base_mx
                dkp = float(np.max(blk)) - base_pure
                w_ = np.exp(sc - np.max(sc)); w_ /= w_.sum()
                pfpm3 = float(w_[:, :, FPM_GRID == 3.0, :, :].sum())
                pkwf = float(w_[:, :, :, KW_GRID == 0.7, :].sum())
                mu_p, mh_p, lb, lc = census_locus(pq, o_f, fcm_h, kw_h,
                                                  fpm_h, sq_h, pq['gs'])
                ok5 = (abs(lb.sum()-mu_p) <= 1e-9
                       and abs(lc.sum()-mh_p) <= 1e-9)
                g5_ok &= ok5
                edge = []
                if cix[4] == len(SQX)-1: edge.append('sq=0.5-EDGE')
                if abs(fpm_h-3.0) < 1e-9: edge.append('fpm=3.0-EDGE')
                if abs(kw_h-0.7) < 1e-9: edge.append('kw-floor')
                if cix[0] in (0, len(FCOMP_GRID)-1): edge.append('fcomp-edge')
                P(f"[{law} {seed}] {pp:<8}: Dkin={dk:+8.2f} (pure "
                  f"{dkp:+8.2f}); cell fcomp={fcm_h:.2f} fpm={fpm_h:.1f} "
                  f"kw={kw_h:.1f} sq={sq_h:.1f}; P(fpm=3.0)={pfpm3:.2f} "
                  f"P(kw=0.7)={pkwf:.2f}; census@prof mu=({mu_p:.2f}, "
                  f"{mh_p:.2f}) pmf9={pmf(NOBS, mu_p):.1e} "
                  f"pmf2={pmf(NHI, mh_p):.1e} jP={jointp(mu_p, mh_p):.1e}"
                  + (("; " + ",".join(edge)) if edge else ""))
                P(f"[{law} {seed}] {pp:<8} locus band : "
                  + "  ".join(f"{n}={v:.2f}" for n, v in zip(PLBL, lb)))
                P(f"[{law} {seed}] {pp:<8} locus cliff: "
                  + "  ".join(f"{n}={v:.2f}" for n, v in zip(PLBL, lc)))
                R[pp].append(dict(law=law, seed=seed, dk=dk, mu_p=mu_p,
                                  mh_p=mh_p, pfpm3=pfpm3, pkwf=pkwf,
                                  p2=pmf(NHI, mh_p)))
            P("")
    if not (g0_ok and g1_ok and g2_ok and g4_ok and g5_ok):
        P("GATES FAILED (G8Q-0 " + ('PASS' if g0_ok else 'FAIL')
          + ", G8Q-1 " + ('PASS' if g1_ok else 'FAIL')
          + ", G8Q-2 " + ('PASS' if g2_ok else 'FAIL')
          + ", G8Q-4 " + ('PASS' if g4_ok else 'FAIL')
          + ", G8Q-5 " + ('PASS' if g5_ok else 'FAIL')
          + ") - STOP; DO NOT QUOTE; verdict WITHHELD")
    else:
        P("GATES: G8Q-0 4/4, G8Q-1 4/4, G8Q-2 12/12, G8Q-4 2/2, "
          "G8Q-5 all - ALL PASS")
        P("")
        n_q1 = sum(1 for rt, rr in zip(R['tokL11'], R['raghavan'])
                   if (rt['mh_p'] <= 0.5*rr['mh_p'] and rt['p2'] >= 1e-3))
        n_acc = sum(1 for x in R['tokL11'] if x['dk'] >= -2.0)
        n_rej = sum(1 for x in R['tokL11'] if x['dk'] <= -5.0)
        n_q3 = sum(1 for x in R['tokL11']
                   if (x['pfpm3'] < 0.5 or x['pkwf'] < 0.5))
        q1 = 'CLIFF-REPAIRED' if n_q1 >= 3 else 'NO-REPAIR'
        q2 = ('ACCEPTED' if n_acc >= 3 else
              ('REJECTED' if n_rej >= 3 else 'CARRIED'))
        P(f"Q1 (tokL11 cliff <= 0.5x raghavan AND pmf2 >= 1e-3): "
          f"{n_q1}/4 -> {q1}")
        P(f"Q2 (tokL11 kinematics): accept-rows {n_acc}/4, reject-rows "
          f"{n_rej}/4 -> {q2}")
        P(f"Q3 (secondary, informational): edge relaxation in "
          f"{n_q3}/4 rows (P(fpm=3.0) < 0.5 or kw off floor)")
        P(f"==> 8Q VERDICT (locked grammar): {q1} x {q2}; the mid "
          f"bracket rows stand as the outer-separation-caveat envelope. "
          f"NO credence movement (pre-stated).")
    P(f"\ndone ({(time.time()-t0)/60:.1f} min)")
    with open('data/stage8q_pprior.txt', 'w') as f:
        f.write("\n".join(L_) + "\n")
    print("\nsaved: data/stage8q_pprior.txt")
else:
    # ------------------- COMBINED LEG ------------------------------
    win = None
    for ln in open('data/stage8p_sqshape.txt').read().splitlines():
        mm = re.match(r"WINNER: (\w+)$", ln)
        if mm:
            win = mm.group(1)
    assert win is not None, 'no WINNER line in the 8P record'
    shp = 'logn' if win == 'none' else win
    P("8Q COMBINED LEG (rule pre-stated: width shape = 8P WINNER, "
      f"population prior = tokL11): shape = '{shp}' (8P line: {win})")
    P("")
    g_ok = True
    rowsC = []
    for seed in SEEDS:
        pf = build_pop(seed, 'raghavan')
        pq = build_pop(seed, 'tokL11')
        tm_f = shape_maps(pf['gs'])
        tm_q = shape_maps(pq['gs'])
        for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
            c9 = np.load(f'data/stage7j_cube_full_lker_{seed}_{law}.npy')
            pri = (prior_eta.reshape((1, 2) + (1,)*7)
                   + LNPI.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1)))
            bx = np.unravel_index(np.nanargmax(c9 + pri), c9.shape)
            al = A_GRID[bx[0]]; eta = E_GRID[bx[1]]; wr = WR_GRID[bx[2]]
            e_f = e_of_x(pf, eta, wr)
            tab_a = 1.0 + al*(TAB-1.0)
            vp_f = vp_c(pf, e_f, tab_a)
            o_f = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'], pf['M_s'],
                      pf['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_a,
                      lny0=LNY0, dlny=DLNY, vp=vp_f)
            blk0 = eval_block(pf, o_f, {'logn': tm_f['logn']})['logn']
            sl = np.squeeze(c9[bx[0], bx[1], bx[2]], axis=1)
            dmax = float(np.max(np.abs(blk0[..., :4] - sl)))
            okg = dmax <= 0.05
            g_ok &= okg
            P(f"[{law} {seed}] G8Qc-0 baseline cube identity: max|d| = "
              f"{dmax:.2e} -> {'PASS' if okg else 'FAIL'}")
            base = float(np.max(blk0 + LNPI.reshape(6, 1, 1, 1, 1)))
            blkC = eval_block(pq, o_f, {shp: tm_q[shp]})[shp]
            scC = blkC + LNPI.reshape(6, 1, 1, 1, 1)
            mxC = float(np.max(scC))
            cix = np.unravel_index(np.argmax(scC), scC.shape)
            fcm_h = FCOMP_GRID[cix[0]]; fpm_h = FPM_GRID[cix[2]]
            kw_h = KW_GRID[cix[3]]; sq_h = SQX[cix[4]]
            dk = mxC - base
            mu_p, mh_p, lb, lc = census_locus(pq, o_f, fcm_h, kw_h,
                                              fpm_h, sq_h, tm_q[shp])
            jp = jointp(mu_p, mh_p)
            P(f"[{law} {seed}] COMBINED ({shp}+tokL11): Dkin={dk:+8.2f} "
              f"vs operative; cell fcomp={fcm_h:.2f} fpm={fpm_h:.1f} "
              f"kw={kw_h:.1f} sq={sq_h:.1f}; census mu=({mu_p:.2f}, "
              f"{mh_p:.2f}) jointP={jp:.1e}")
            P(f"[{law} {seed}] COMBINED locus band : "
              + "  ".join(f"{n}={v:.2f}" for n, v in zip(PLBL, lb)))
            P(f"[{law} {seed}] COMBINED locus cliff: "
              + "  ".join(f"{n}={v:.2f}" for n, v in zip(PLBL, lc)))
            rowsC.append(dict(law=law, seed=seed, dk=dk, jp=jp))
        P("")
    if not g_ok:
        P("GATES FAILED (G8Qc-0) - STOP; DO NOT QUOTE; verdict WITHHELD")
    else:
        n_jc = sum(1 for x in rowsC if (x['jp'] >= 1e-3 and
                                        x['dk'] >= -5.0))
        P(f"COMBINED bar (jointP >= 1e-3 AND Dkin >= -5): {n_jc}/4")
        if n_jc >= 3:
            P("==> 8Q COMBINED VERDICT: JOINT-COHERENT — the repaired "
              "model (measured subsystem P-prior + the 8P width shape) "
              "reproduces the (band, cliff) census pair at the 8M "
              "admissibility grade WITHOUT leaving the kinematic "
              "posterior's neighborhood.  8M's inconsistency is CLOSED "
              "at point grade; the cube-grade re-run (alpha "
              "re-derivation under the repaired model) is the named "
              "successor and the final-stamp decider.")
        else:
            P("==> 8Q COMBINED VERDICT: JOINT-NOT-REACHED — the "
              "repaired model does not reach census coherence inside "
              "the kinematic neighborhood; the 8M price stands and the "
              "residual flank(s) are read from the locus rows.")
        P("    NO credence movement (pre-stated).")
    P(f"\ndone ({(time.time()-t0)/60:.1f} min)")
    with open('data/stage8q_combined.txt', 'w') as f:
        f.write("\n".join(L_) + "\n")
    print("\nsaved: data/stage8q_combined.txt")
