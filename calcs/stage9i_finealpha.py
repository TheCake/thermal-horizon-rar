"""STAGE 9I — THE FINE-ALPHA DROP SCAN (the 9F successor).
Pre-registered BEFORE any run (decider-adjacent; pre-signed map
56 / 50 / hold-53).

9F's drop world rejects alpha >= 0.5 but its 5-pt grid cannot see
below 0.5.  A_FINE = 0..0.6 step 0.1; PRIMARY = M-DROP fine
marginal (Q1-Q3, per-stratum noise, LNPI on shared fcomp).
POWER GATE (7J fullpow standard): synthetic Q1-Q3 sky at truth
alpha = 0.3 (pre-stated nuisances, multinomial at real counts,
rng 9); null letters blocked unless the machinery recovers it.

Gates: G9I-0 lineage BIT-identity to the 9F npz at shared alphas
{0.0, 0.5}; G9I-1 analytic combiner 1e-12; G9I-2 synthetic census;
G9I-3 cube identity at alpha = 0.5 (bit).
Bars (locked, ORDERED): I-POWER-FAIL / I-SMALL-ALPHA /
I-NEWTON-FLAT / I-GRAY-CARRIED.
Output: data/stage9i_finealpha.txt
"""
import math, re
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

ALLM = np.ones(len(s_d), dtype=bool)
data_2d, noise_pool, UNI_B, FLY_B, NDATA = build_stratum(ALLM)

N = 500_000
A_FINE = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
A_GRID5 = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID  = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FC0 = 0.10
FFLY_GRID = np.array([0.05, 0.10])
FPM_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW_GRID = np.array([0.7, 1.0, 1.4])
SQ_GRID = np.array([0.0, 0.1, 0.2, 0.3])
WS2_GRID = np.array([0.0, 0.045])
SEEDS = (31, 101)
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

def eval_block_g(p, prj, D2, PLs, UB, FB, FPM_G, SQ_G):
    """generalized grids; VERBATIM 9A/8Z block at the std grids."""
    smag, vpar, vper = prj
    s_kau = smag/1e3
    out = np.zeros((len(FCOMP_GRID), len(FFLY_GRID), len(FPM_G),
                    len(KW_GRID), len(SQ_G), len(WS2_GRID)))
    for bi, b in enumerate(SBINS):
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
                for pi, fpm in enumerate(FPM_G):
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
                        for si, sqv in enumerate(SQ_G):
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

def eval_pp_cell(p, prj, PLs, UB, FB, cell):
    """pp per s-bin at ONE cell (fi, yi, fpm, ki, sqv, ws=0) —
    the VERBATIM pp construction path of eval_block_g."""
    fi, yi, fpm, ki, sqv = cell
    fcm = FCOMP_GRID[fi]; kwv = KW_GRID[ki]; ff = FFLY_GRID[yi]
    smag, vpar, vper = prj
    s_kau = smag/1e3
    pps = []
    for bi, b in enumerate(SBINS):
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
        vp_a = vpar[idx] + kwv*cvp
        vq_a = vper[idx] + kwv*cvq
        vp_n = vp_a*boost + g1_i*sg0*fpm
        vq_n = vq_a*boost + g2_i*sg0*fpm
        vmag = np.hypot(vp_n, vq_n)
        keep = vmag*4.74047 <= (2.978/np.sqrt(sk_i)
                                + 2.8284*sg0*4.74047)
        vtn = (vmag/vc)[keep]
        gmn = np.degrees(np.arccos(np.clip(
            np.abs(vp_n[keep])
            / np.maximum(vmag[keep], 1e-12), 0, 1)))
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
        pp = (1-wtot)*p0 + wtot*mixc
        pps.append(pp)
    return pps

ship = {}
for ln in open('data/stage8lb_read.txt').read().splitlines():
    mm = re.match(r"\[(simple|BE) (\d+)\] MAP cell: alpha=([\d.]+), "
                  r"wr=([\d.]+), fcomp=([\d.]+), fpm=([\d.]+), "
                  r"kw=([\d.]+), sq=([\d.]+); CENSUS corrected-kernel: "
                  r"mu=\(([\d.]+), ([\d.]+)\)", ln)
    if mm:
        ship[(mm.group(1), int(mm.group(2)))] = dict(
            alpha=float(mm.group(3)), wr=float(mm.group(4)))
assert len(ship) == 4, ship

okv = ok
ruwe = np.maximum(_pick('ruwe1', 'RUWE1'),
                  _pick('ruwe2', 'RUWE2'))[okv]
qs_ = np.percentile(ruwe, [25, 50, 75])
QMASKS = [ruwe <= qs_[0],
          (ruwe > qs_[0]) & (ruwe <= qs_[1]),
          (ruwe > qs_[1]) & (ruwe <= qs_[2]),
          ruwe > qs_[2]]
STRATA = [build_stratum(mk) for mk in QMASKS]

t0 = time.time()
P("9I THE FINE-ALPHA DROP SCAN (pre-reg committed BEFORE any run; "
  "decider-adjacent — movement by the pre-signed map ONLY)")
P("")

def a_hat_g(prof, grid):
    i = int(np.argmax(prof))
    ah = float(grid[i]); edge = (i == 0 or i == len(grid)-1)
    if not edge:
        x3, y3 = grid[i-1:i+2], prof[i-1:i+2]
        c2_, c1_, _ = np.polyfit(x3, y3, 2)
        if c2_ < 0:
            aa = -c1_/(2*c2_)
            if grid[i-1] <= aa <= grid[i+1]:
                ah = float(aa)
    return ah, edge, i

def lse_cells(T, axes, ncells):
    return logsumexp(T, axis=axes) - math.log(ncells)
_T = np.log(np.array([[1.0, 2.0], [3.0, 4.0]]))
_v = float(lse_cells(_T, (0, 1), 4))
g1_ok = abs(_v - math.log(2.5)) < 1e-12
P(f"G9I-1 combiner analytic unit check: lse = {_v:.12f} vs "
  f"ln(2.5) = {math.log(2.5):.12f} -> {'PASS' if g1_ok else 'FAIL'}")
P("")

NC_STD = 6*4
NC_SHR = 2*3*2
NA = len(A_FINE)

def marg_config(STDz, strata):
    lnZ = np.zeros(NA)
    for ai in range(NA):
        S = LNPI[:, None, None, None].copy()
        for qi in strata:
            S = S + lse_cells(STDz[ai, qi], (2, 4), NC_STD)
        lnZ[ai] = logsumexp(S) - math.log(NC_SHR)
    return lnZ

def summarize(lnZ):
    w = np.exp(lnZ - logsumexp(lnZ))
    am = float(np.sum(w*A_FINE))
    p0 = float(w[0])
    p01 = float(w[0] + w[1])
    dn = float((logsumexp(lnZ[1:]) - math.log(NA-1)) - lnZ[0])
    return w, am, p0, p01, dn

# synthetic truth cell (pre-registered): alpha = 0.3 (A_FINE idx 3);
# fcomp = 0.20 (fi 2), ffly = 0.05 (yi 0), kw = 1.0 (ki 1), ws = 0;
# per-stratum (fpm, sq): Q1 (2.1, 0.0), Q2 (1.8, 0.1), Q3 (1.8, 0.1)
AI_TRUE = 3
TRUTH = {0: (2, 0, 2.1, 1, 0.0),
         1: (2, 0, 1.8, 1, 0.1),
         2: (2, 0, 1.8, 1, 0.1)}

g0_ok = g2_ok = g3_ok = True
rows = []
power = {}
for seed in SEEDS:
    pf = build_pop(seed)
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        c9 = np.load(f'data/stage7j_cube_full_lker_{seed}_{law}.npy')
        pri = (prior_eta.reshape((1, 2) + (1,)*7)
               + LNPI.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1)))
        bx = np.unravel_index(np.nanargmax(c9 + pri), c9.shape)
        al_map = A_GRID5[bx[0]]; eta = E_GRID[bx[1]]
        wr = WR_GRID[bx[2]]
        sh = ship[(law, seed)]
        assert all(abs(v-sh[k]) < 1e-9 for v, k in
                   ((al_map, 'alpha'), (wr, 'wr'))), (bx, sh)
        e_f = e_of_x(pf, eta, wr)
        z9f = np.load(f'data/stage9f_tables_{seed}_{law}.npz')
        STD9F = z9f['STD']
        STD = np.zeros((NA, 4, 6, 2, 6, 3, 4, 2))
        SYN = np.zeros((NA, 3, 6, 2, 6, 3, 4, 2))
        syn_D2 = None
        prj_by_ai = {}
        for ai, al in enumerate(A_FINE):
            tab_a = 1.0 + al*(TAB-1.0)
            vp_f = vp_c(pf, e_f, tab_a)
            o_f = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'],
                      pf['M_s'], pf['uph'], 8, 2500, 5, a0=A0_CAN,
                      tab=tab_a, lny0=LNY0, dlny=DLNY, vp=vp_f)
            prj = project(pf, o_f)
            for qi in range(4):
                STD[ai, qi] = eval_block_g(pf, prj, *STRATA[qi][:4],
                                           FPM_GRID, SQ_GRID)
            if ai == AI_TRUE:
                rng9 = np.random.default_rng(9)
                syn_D2 = []
                for qi in range(3):
                    pps = eval_pp_cell(pf, prj, STRATA[qi][1],
                                       STRATA[qi][2], STRATA[qi][3],
                                       TRUTH[qi])
                    D2s = []
                    for bi in range(len(SBINS)):
                        nd = int(STRATA[qi][0][bi].sum())
                        pp = pps[bi]
                        draw = rng9.multinomial(nd, pp.ravel()
                                                / pp.sum())
                        D2s.append(draw.reshape(NV, NG)
                                   .astype(float))
                    syn_D2.append(D2s)
                    ok2 = sum(int(h.sum()) for h in D2s) == \
                        sum(int(h.sum()) for h in STRATA[qi][0])
                    g2_ok &= ok2
                P(f"[{law} {seed}] G9I-2 synthetic census "
                  f"(counts == ND, 3 strata) -> "
                  f"{'PASS' if g2_ok else 'FAIL'}")
            if abs(al - al_map) < 1e-9:
                blk_u = eval_block_g(pf, prj, data_2d, noise_pool,
                                     UNI_B, FLY_B, FPM_GRID,
                                     SQ_GRID)
                sl = np.squeeze(c9[bx[0], bx[1], bx[2]], axis=1)
                dmax = float(np.max(np.abs(blk_u[..., 0] - sl)))
                ok3 = dmax <= 1e-9
                g3_ok &= ok3
                P(f"[{law} {seed}] G9I-3 cube identity: max|d| = "
                  f"{dmax:.2e} -> {'PASS' if ok3 else 'FAIL'}")
            prj_by_ai[ai] = prj
        # synthetic tables at every alpha (Q1-Q3)
        for ai in range(NA):
            for qi in range(3):
                SYN[ai, qi] = eval_block_g(
                    pf, prj_by_ai[ai],
                    syn_D2[qi], STRATA[qi][1], STRATA[qi][2],
                    STRATA[qi][3], FPM_GRID, SQ_GRID)
        del prj_by_ai
        # G9I-0 lineage bit-identity at shared alphas
        d00 = float(np.max(np.abs(STD[0] - STD9F[0])))
        d05 = float(np.max(np.abs(STD[5] - STD9F[1])))
        ok0 = (d00 == 0.0) and (d05 == 0.0)
        g0_ok &= ok0
        P(f"[{law} {seed}] G9I-0 lineage bit-identity: d(a=0.0) = "
          f"{d00:.2e}, d(a=0.5) = {d05:.2e} -> "
          f"{'PASS' if ok0 else 'FAIL'}  "
          f"({(time.time()-t0)/60:.1f} min)")
        np.savez(f'data/stage9i_tables_{seed}_{law}.npz',
                 STD=STD, SYN=SYN, A_FINE=A_FINE, LNPI=LNPI,
                 eta=eta, wr=wr,
                 synD2=np.array([[syn_D2[qi][bi] for bi in
                                  range(4)] for qi in range(3)]))
        # reads
        zdrop = marg_config(STD, [0, 1, 2])
        zstd = marg_config(STD, [0, 1, 2, 3])
        zsyn = marg_config(SYN, [0, 1, 2])
        wd_, am_d, p0_d, p01_d, dn_d = summarize(zdrop)
        ws_, am_s, p0_s, p01_s, dn_s = summarize(zstd)
        wy_, am_y, p0_y, p01_y, dn_y = summarize(zsyn)
        curve = zdrop - zdrop[0]
        P(f"[{law} {seed}] M-DROP lnZ(a)-lnZ(0) = "
          + "/".join(f"{x:+.1f}" for x in curve))
        P(f"[{law} {seed}] M-DROP fine: w(a) = "
          + "/".join(f"{x:.3f}" for x in wd_)
          + f"; a_marg = {am_d:.3f}, P(a=0) = {p0_d:.3f}, "
            f"P(a<=0.1) = {p01_d:.3f}, dN = {dn_d:+.1f}")
        P(f"[{law} {seed}] M-STD  fine: w(a) = "
          + "/".join(f"{x:.3f}" for x in ws_)
          + f"; a_marg = {am_s:.3f}, P(a=0) = {p0_s:.3f}, "
            f"dN = {dn_s:+.1f}")
        for qi in range(4):
            lnq = np.zeros(NA)
            for ai in range(NA):
                lnq[ai] = float(logsumexp(
                    STD[ai, qi]
                    + LNPI.reshape(6, 1, 1, 1, 1, 1))
                    - math.log(NC_STD*NC_SHR))
            wq = np.exp(lnq - logsumexp(lnq))
            amq = float(np.sum(wq*A_FINE))
            P(f"    Q{qi+1}-alone: w(a) = "
              + "/".join(f"{x:.2f}" for x in wq)
              + f"; a_marg = {amq:.3f}")
        prof_d = np.array([float(np.max(
            sum(STD[ai, qi].max(axis=(2, 4)) for qi in (0, 1, 2))
            + LNPI.reshape(6, 1, 1, 1))) for ai in range(NA)])
        ah_d, ed_d, _ = a_hat_g(prof_d, A_FINE)
        P(f"[{law} {seed}] profile a_hat(drop, fine) = {ah_d:.3f}"
          + (" EDGE" if ed_d else ""))
        okp = (am_y >= 0.15) and (p0_y <= 0.35)
        power[(law, seed)] = okp
        P(f"[{law} {seed}] G9I-P power (truth a=0.3): w(a) = "
          + "/".join(f"{x:.3f}" for x in wy_)
          + f"; recovered a_marg = {am_y:.3f}, P(a=0) = {p0_y:.3f} "
            f"-> {'PASS' if okp else 'FAIL'}")
        rows.append(dict(law=law, seed=seed, am=am_d, p0=p0_d,
                         p01=p01_d, pw=okp))
        P("")

if not (g0_ok and g1_ok and g2_ok and g3_ok):
    P("GATES FAILED (G9I-0 " + ('PASS' if g0_ok else 'FAIL')
      + ", G9I-1 " + ('PASS' if g1_ok else 'FAIL')
      + ", G9I-2 " + ('PASS' if g2_ok else 'FAIL')
      + ", G9I-3 " + ('PASS' if g3_ok else 'FAIL')
      + ") - STOP; DO NOT QUOTE; verdict WITHHELD")
else:
    P("GATES: G9I-0 4/4, G9I-1, G9I-2 4/4, G9I-3 4/4 - ALL PASS")
    P("")
    npf = sum(1 for r in rows if not r['pw'])
    n2 = sum(1 for r in rows if (0.10 <= r['am'] <= 0.45
                                 and r['p0'] <= 0.35))
    n3 = sum(1 for r in rows if (r['p01'] >= 0.60 and r['pw']))
    P(f"I-bars: POWER-FAIL rows {npf}/4; SMALL-ALPHA rows {n2}/4; "
      f"NEWTON-FLAT rows {n3}/4")
    if npf >= 2:
        P("==> 9I VERDICT (locked grammar): POWER-FAIL-CARRIED - "
          "the machinery cannot recover the injected small alpha; "
          "null letters BLOCKED; rows stand as measurements.")
        P("    CREDENCE (pre-signed map): HOLD 53.")
    elif n2 >= 3:
        P("==> 9I VERDICT (locked grammar): I-SMALL-ALPHA - the "
          "clean-quartile kinematics ALLOW a small interior boost; "
          "the fitted channel and the model-light core reconcile "
          "at deflated amplitude.")
        P("    CREDENCE (pre-signed map): anomaly-real 53 -> 56.")
    elif n3 >= 3:
        P("==> 9I VERDICT (locked grammar): I-NEWTON-FLAT - the "
          "clean-quartile kinematics exclude even small boosts at "
          "fitted grade (power demonstrated); the fitted channel "
          "is fully quality-carried.")
        P("    CREDENCE (pre-signed map): anomaly-real 53 -> 50.")
    else:
        P("==> 9I VERDICT (locked grammar): I-GRAY-CARRIED - rows "
          "stand as measurements.")
        P("    CREDENCE (pre-signed map): HOLD 53.")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage9i_finealpha.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage9i_finealpha.txt")
