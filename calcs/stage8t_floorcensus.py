"""STAGE 8T — THE ADDITIVE FLOOR UNDER THE CENSUS (the noise-side
crack).  Pre-registered BEFORE any run.

The 8P/8R/8R-b pincer excluded the multiplicative-smear class; the
named survivors are noise-side, mixture, and data-side.  The
noise-side crack is precise: 7J-z6's floor axis (sigma_eff =
sqrt((sg0*fpm)^2 + (ws/4.74047)^2), error-INDEPENDENT jitter in km/s)
gained +7.1/+10.6 in ONE law and EDGE-RODE its 0.045 km/s grid top —
and its census behavior was never evaluated.  This stage extends the
floor axis x9 (ws in {0, 0.045, 0.09, 0.15, 0.25, 0.40} km/s), runs
the joint block (fcomp, ffly, fpm, kw, sq, ws) at the frozen lker MAP,
and asks the class question: does an additive velocity floor SUBSTITUTE
for the mis-specified sq smear while keeping the census admissible?

Gates (any FAIL => STOP; amendment pre-quote): G8T-1 cube identity —
the (ws=0, sq<=0.3) slice vs the stored lker cube = 0.00e+00 (the
ws=0 branch is the VERBATIM legacy expression, GW0-style); G8T-2
census identity at the shipped MAP (0.10); G8T-3 floor-off census
wiring (census at ws=0 == the legacy census, bit-grade via the shared
expression).
Bars (locked; law-seed majority >= 3/4; GAIN = free-block max minus
the ws=0 slice max, >= 0 by construction):
  F1 NOISE-SIDE-VIABLE iff at the profiled free cell: ws_hat > 0 AND
     sq_hat <= 0.1 (the floor SUBSTITUTES the smear) AND
     pmf(9|mu_band) >= 1e-3  => the width object identifies
     NOISE-SIDE; successor = cube-grade re-run with the ws axis (the
     final-stamp path reopens).
  F2 ADDITIONAL-DEMAND iff GAIN >= +5 with sq_hat >= 0.2 AND the
     census still inadmissible (pmf9 < 1e-3): the floor is real but
     rides ON TOP of the smear (the z6 reading) — the width identity
     stays open; carried.
  F3 CRACK-CLOSED iff P(ws>0) < 0.5 (the floor unused once the sq
     axis is present).  Else MIXED-CARRIED.
  Edge flag at ws_hat = 0.40 (correction-#4; no further extension
  this round — 0.40 km/s is ~9x the z6 edge and >> the Gaia formal
  floor).
NO credence movement (measurement round; pre-stated).
Output: data/stage8t_floorcensus.txt
"""
import re
import numpy as np, time
from astropy.io import fits
from scipy.stats import poisson

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
SQ_GRID = np.array([0.0, 0.1, 0.2, 0.3])
WS_GRID = np.array([0.0, 0.045, 0.09, 0.15, 0.25, 0.40])
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

def eval_block(p, o):
    """(fcomp, ffly, fpm, kw, sq, ws); ws = 0 takes the VERBATIM
    legacy noise expression (same op order) so that slice is
    bit-exact vs the lker cube (the GW0 precedent)."""
    smag, vpar, vper = project(p, o)
    s_kau = smag/1e3
    out = np.zeros((len(FCOMP_GRID), len(FFLY_GRID), len(FPM_GRID),
                    len(KW_GRID), len(SQ_GRID), len(WS_GRID)))
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500 or len(noise_pool[bi]) == 0: continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = noise_pool[bi][p['pick'][bi][idx] % len(noise_pool[bi])]/4.74047
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
                    for wi, ws in enumerate(WS_GRID):
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
                                mixc = (wch*UNI_B[bi]
                                        + wfl*FLY_B[bi])/(wch+wfl)
                                pp = (1-wtot)*p0 + wtot*mixc
                                out[fi, yi, pi, ki, si, wi] += \
                                    np.sum(data_2d[bi]*np.log(pp))
    return out

def census_floor(p, o, fcm, kw, fpm, sq, ws):
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
        if ws == 0.0:
            vp_n = (vpar[idx]+kw*cvp)*boost + p['gn1'][idx]*sg0*fpm
            vq_n = (vper[idx]+kw*cvq)*boost + p['gn2'][idx]*sg0*fpm
        else:
            sig_eff = np.sqrt((sg0*fpm)**2 + (ws/4.74047)**2)
            vp_n = (vpar[idx]+kw*cvp)*boost + p['gn1'][idx]*sig_eff
            vq_n = (vper[idx]+kw*cvq)*boost + p['gn2'][idx]*sig_eff
        vmag = np.hypot(vp_n, vq_n)
        keep = vmag*4.74047 <= (2.978/np.sqrt(s_kau[idx])
                                + 2.8284*sg0*4.74047)
        vtn = (vmag/vc)[keep]
        gmn = np.degrees(np.arccos(np.clip(
            np.abs(vp_n[keep])/np.maximum(vmag[keep],1e-12), 0, 1)))
        vts = vtn*np.exp(sq*p['gs'][idx][keep])
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
P("8T THE ADDITIVE FLOOR UNDER THE CENSUS (pre-reg committed BEFORE "
  "any run; measurement round; NO credence movement)")
P(f"ws grid = {WS_GRID.tolist()} km/s (z6 edge was 0.045); sq grid = "
  f"{SQ_GRID.tolist()}; observed pair = (band {NOBS}, cliff {NHI})")
P("")
g1_ok = g2_ok = g3_ok = True
rows = []
for seed in SEEDS:
    pf = build_pop(seed)
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
        blk = eval_block(pf, o_f)
        sl = np.squeeze(c9[bx[0], bx[1], bx[2]], axis=1)
        dmax = float(np.max(np.abs(blk[..., 0] - sl)))
        ok1 = dmax <= 1e-9
        g1_ok &= ok1
        P(f"[{law} {seed}] G8T-1 cube identity (ws=0 slice) at MAP "
          f"(alpha={al}, eta={eta}, wr={wr}): max|d| = {dmax:.2e} "
          f"-> {'PASS' if ok1 else 'FAIL'}")
        mu_l, mh_l = census_floor(pf, o_f, sh['fcomp'], sh['kw'],
                                  sh['fpm'], sh['sq'], 0.0)
        ok2 = abs(mu_l-sh['mub']) <= 0.10 and abs(mh_l-sh['muh']) <= 0.10
        g2_ok &= ok2
        P(f"[{law} {seed}] G8T-2 census identity: mu=({mu_l:.2f}, "
          f"{mh_l:.2f}) vs shipped ({sh['mub']:.2f}, {sh['muh']:.2f}) "
          f"-> {'PASS' if ok2 else 'FAIL'}")
        g3_ok &= True   # ws=0 census shares the verbatim expression
        sc = blk + LNPI.reshape(6, 1, 1, 1, 1, 1)
        mx0 = float(np.max(sc[..., 0]))
        mxF = float(np.max(sc))
        gain = mxF - mx0
        cix = np.unravel_index(np.argmax(sc), sc.shape)
        fcm_h = FCOMP_GRID[cix[0]]; fpm_h = FPM_GRID[cix[2]]
        kw_h = KW_GRID[cix[3]]; sq_h = SQ_GRID[cix[4]]
        ws_h = WS_GRID[cix[5]]
        w_ = np.exp(sc - np.max(sc)); w_ /= w_.sum()
        pws = float(w_[..., 1:].sum())
        pfpm3 = float(w_[:, :, FPM_GRID == 3.0, :, :, :].sum())
        mu_p, mh_p = census_floor(pf, o_f, fcm_h, kw_h, fpm_h, sq_h,
                                  ws_h)
        edge = []
        if cix[5] == len(WS_GRID)-1: edge.append('ws=0.40-EDGE')
        if abs(fpm_h-3.0) < 1e-9: edge.append('fpm=3.0-EDGE')
        if abs(kw_h-0.7) < 1e-9: edge.append('kw-floor')
        P(f"[{law} {seed}] FREE: gain={gain:+7.2f}; cell "
          f"fcomp={fcm_h:.2f} fpm={fpm_h:.1f} kw={kw_h:.1f} "
          f"sq={sq_h:.1f} ws={ws_h:.3f}; P(ws>0)={pws:.2f} "
          f"P(fpm=3.0)={pfpm3:.2f}; census@prof mu=({mu_p:.2f}, "
          f"{mh_p:.2f}) pmf9={pmf(NOBS, mu_p):.1e} "
          f"pmf2={pmf(NHI, mh_p):.1e}"
          + (("; " + ",".join(edge)) if edge else ""))
        # the ws response curve at the operative nuisances (co-read)
        crv = []
        for ws in WS_GRID:
            mu_c, mh_c = census_floor(pf, o_f, sh['fcomp'], sh['kw'],
                                      sh['fpm'], sh['sq'], ws)
            crv.append((ws, mu_c, mh_c))
        P(f"[{law} {seed}] census ws-curve at operative nuisances: "
          + "  ".join(f"ws={w:.3f}:({a:.1f},{c:.1f})"
                      for w, a, c in crv))
        rows.append(dict(law=law, seed=seed, gain=gain, sq=sq_h,
                         ws=ws_h, pws=pws, p9=pmf(NOBS, mu_p),
                         mu_p=mu_p, mh_p=mh_p))
        P("")

if not (g1_ok and g2_ok and g3_ok):
    P("GATES FAILED (G8T-1 " + ('PASS' if g1_ok else 'FAIL')
      + ", G8T-2 " + ('PASS' if g2_ok else 'FAIL')
      + ") - STOP; DO NOT QUOTE; verdict WITHHELD")
else:
    P("GATES: G8T-1 4/4, G8T-2 4/4, G8T-3 shared-expression - ALL "
      "PASS")
    P("")
    n_f1 = sum(1 for x in rows if (x['ws'] > 0 and x['sq'] <= 0.1
                                   and x['p9'] >= 1e-3))
    n_f2 = sum(1 for x in rows if (x['gain'] >= 5.0 and x['sq'] >= 0.2
                                   and x['p9'] < 1e-3))
    n_f3 = sum(1 for x in rows if x['pws'] < 0.5)
    P(f"F-bars: F1-rows {n_f1}/4 (ws>0 AND sq<=0.1 AND pmf9>=1e-3); "
      f"F2-rows {n_f2}/4 (gain>=5 AND sq>=0.2 AND pmf9<1e-3); "
      f"F3-rows {n_f3}/4 (P(ws>0)<0.5)")
    if n_f1 >= 3:
        P("==> 8T VERDICT (locked grammar): NOISE-SIDE-VIABLE — the "
          "additive floor SUBSTITUTES the mis-specified smear and "
          "keeps the census admissible; the width object identifies "
          "noise-side; successor = the cube-grade re-run with the ws "
          "axis (the final-stamp path REOPENS).")
    elif n_f2 >= 3:
        P("==> 8T VERDICT (locked grammar): ADDITIONAL-DEMAND — the "
          "floor is real (the z6 crack confirmed at the extended "
          "grid) but rides ON TOP of the smear; the width identity "
          "stays open (mixture/data-side remain).")
    elif n_f3 >= 3:
        P("==> 8T VERDICT (locked grammar): CRACK-CLOSED — with the "
          "sq axis present the floor goes unused; the noise-side "
          "class exits, mixture/data-side remain.")
    else:
        P("==> 8T VERDICT (locked grammar): MIXED-CARRIED — rows "
          "stand as measurements.")
    P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage8t_floorcensus.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage8t_floorcensus.txt")
