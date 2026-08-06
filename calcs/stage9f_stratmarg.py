"""STAGE 9F — THE MARGINAL STRATIFIED RE-RUN (binary; THE ARBITER
9D fired).  Pre-registered BEFORE any run (decider round; movement
by the pre-signed credence map ONLY).

9A/9D read the stratified-noise model at PROFILE grade (max over
per-stratum noise cells).  The max operator is exactly where a
pathological stratum (Q4 grid-top riding) can carry a spurious
optimum.  9F reads the SAME tables at MARGINAL grade: logsumexp
with uniform cell priors on noise axes, LNPI on shared fcomp,
uniform on ffly/kw/ws and on the 5-point alpha grid.  Configs:
M-STD (per-stratum noise), M-DROP (Q1-Q3), M-EXT (Q4 extended),
M-TIED (shared noise = the aggregation world).  Per-stratum noise
posteriors printed (the 9D lean-output repair).

Gates: G9F-0 cube identity at alpha_MAP (bit); G9F-1 profile-mode
regression vs 9D printed a_hat std/drop/ext (12 values, 0.002);
G9F-2 combiner analytic unit check (1e-12); G9F-3 npz round-trip
bit-identity before any combiner read.
Bars (locked, ORDERED): F-COLLAPSED / F-SURVIVES / F-Q4-CARRIED /
F-GRAY-CARRIED.  Credence map: 48 / 63 / 53 / hold-58.
Output: data/stage9f_stratmarg.txt
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
A_GRID  = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID  = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FC0 = 0.10
FFLY_GRID = np.array([0.05, 0.10])
FPM_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW_GRID = np.array([0.7, 1.0, 1.4])
SQ_GRID = np.array([0.0, 0.1, 0.2, 0.3])
FPMX_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0, 3.6, 4.2])
SQX_GRID = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
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
    """generalized grids; at (FPM_GRID, SQ_GRID) this is the
    VERBATIM 9A/8Z block."""
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
g9d = {}
for ln in open('data/stage9d_q4robust.txt').read().splitlines():
    mm = re.match(r"\[(simple|BE) (\d+)\] a_hat std/drop/ext = "
                  r"([\d.]+)/([\d.]+)/([\d.]+)", ln)
    if mm:
        g9d[(mm.group(1), int(mm.group(2)))] = (
            float(mm.group(3)), float(mm.group(4)), float(mm.group(5)))
assert len(g9d) == 4, g9d

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
P("9F THE MARGINAL STRATIFIED RE-RUN (pre-reg committed BEFORE any "
  "run; DECIDER round — movement by the pre-signed map ONLY)")
P("")

def a_hat(prof):
    i = int(np.argmax(prof))
    ah = float(A_GRID[i]); edge = (i == 0 or i == 4)
    if not edge:
        x3, y3 = A_GRID[i-1:i+2], prof[i-1:i+2]
        c2_, c1_, _ = np.polyfit(x3, y3, 2)
        if c2_ < 0:
            aa = -c1_/(2*c2_)
            if A_GRID[i-1] <= aa <= A_GRID[i+1]:
                ah = float(aa)
    return ah, edge, i

# ---------- G9F-2: combiner analytic unit check ----------
def lse_cells(T, axes, ncells):
    return logsumexp(T, axis=axes) - math.log(ncells)
_T = np.log(np.array([[1.0, 2.0], [3.0, 4.0]]))
_v = float(lse_cells(_T, (0, 1), 4))
g2_ok = abs(_v - math.log(2.5)) < 1e-12
P(f"G9F-2 combiner analytic unit check: lse = {_v:.12f} vs "
  f"ln(2.5) = {math.log(2.5):.12f} -> {'PASS' if g2_ok else 'FAIL'}")
P("")

# ---------- Part 1: GPU tables (archived) ----------
g0_ok = True
for seed in SEEDS:
    pf = build_pop(seed)
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        c9 = np.load(f'data/stage7j_cube_full_lker_{seed}_{law}.npy')
        pri = (prior_eta.reshape((1, 2) + (1,)*7)
               + LNPI.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1)))
        bx = np.unravel_index(np.nanargmax(c9 + pri), c9.shape)
        al_map = A_GRID[bx[0]]; eta = E_GRID[bx[1]]
        wr = WR_GRID[bx[2]]
        sh = ship[(law, seed)]
        assert all(abs(v-sh[k]) < 1e-9 for v, k in
                   ((al_map, 'alpha'), (wr, 'wr'))), (bx, sh)
        e_f = e_of_x(pf, eta, wr)
        STD = np.zeros((5, 4, 6, 2, 6, 3, 4, 2))
        EXT = np.zeros((5, 6, 2, 8, 3, 6, 2))
        for ai, al in enumerate(A_GRID):
            tab_a = 1.0 + al*(TAB-1.0)
            vp_f = vp_c(pf, e_f, tab_a)
            o_f = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'],
                      pf['M_s'], pf['uph'], 8, 2500, 5, a0=A0_CAN,
                      tab=tab_a, lny0=LNY0, dlny=DLNY, vp=vp_f)
            prj = project(pf, o_f)
            for qi in range(4):
                STD[ai, qi] = eval_block_g(pf, prj, *STRATA[qi][:4],
                                           FPM_GRID, SQ_GRID)
            EXT[ai] = eval_block_g(pf, prj, *STRATA[3][:4],
                                   FPMX_GRID, SQX_GRID)
            if abs(al - al_map) < 1e-9:
                blk_u = eval_block_g(pf, prj, data_2d, noise_pool,
                                     UNI_B, FLY_B, FPM_GRID,
                                     SQ_GRID)
                sl = np.squeeze(c9[bx[0], bx[1], bx[2]], axis=1)
                dmax = float(np.max(np.abs(blk_u[..., 0] - sl)))
                ok0 = dmax <= 1e-9
                g0_ok &= ok0
                P(f"[{law} {seed}] G9F-0 cube identity: max|d| = "
                  f"{dmax:.2e} -> {'PASS' if ok0 else 'FAIL'}")
        np.savez(f'data/stage9f_tables_{seed}_{law}.npz',
                 STD=STD, EXT=EXT, A_GRID=A_GRID, LNPI=LNPI,
                 FPM_GRID=FPM_GRID, SQ_GRID=SQ_GRID,
                 FPMX_GRID=FPMX_GRID, SQX_GRID=SQX_GRID,
                 eta=eta, wr=wr)
        z = np.load(f'data/stage9f_tables_{seed}_{law}.npz')
        ok3 = (np.array_equal(z['STD'], STD)
               and np.array_equal(z['EXT'], EXT)
               and np.array_equal(z['LNPI'], LNPI))
        if not ok3: g0_ok = False
        P(f"[{law} {seed}] G9F-3 npz round-trip bit-identity -> "
          f"{'PASS' if ok3 else 'FAIL'}  "
          f"({(time.time()-t0)/60:.1f} min)")
P("")

# ---------- Part 2: PURE READER (combiner from npz only) ----------
NC_STD = 6*4     # per-stratum noise cells (fpm x sq), std grids
NC_EXT = 8*6
NC_SHR = 2*3*2   # ffly x kw x ws (uniform)

def marg_config(STDz, EXTz, LNPIz, strata, ext_q4=False):
    """lnZ(alpha) with per-stratum (fpm,sq) marginalized
    independently; shared (fc,fy,kw,ws) marginalized jointly."""
    lnZ = np.zeros(5)
    for ai in range(5):
        S = LNPIz[:, None, None, None].copy()
        for qi in strata:
            if ext_q4 and qi == 3:
                S = S + lse_cells(EXTz[ai], (2, 4), NC_EXT)
            else:
                S = S + lse_cells(STDz[ai, qi], (2, 4), NC_STD)
        lnZ[ai] = logsumexp(S) - math.log(NC_SHR)
    return lnZ

def marg_tied(STDz, LNPIz):
    lnZ = np.zeros(5)
    for ai in range(5):
        Ttot = STDz[ai].sum(axis=0)
        S = (LNPIz[:, None, None, None]
             + lse_cells(Ttot, (2, 4), NC_STD))
        lnZ[ai] = logsumexp(S) - math.log(NC_SHR)
    return lnZ

def summarize(lnZ):
    w = np.exp(lnZ - logsumexp(lnZ))
    am = float(np.sum(w*A_GRID))
    p0 = float(w[0])
    dn = float((logsumexp(lnZ[1:]) - math.log(4)) - lnZ[0])
    return w, am, p0, dn

def stratum_noise_post(STDz, EXTz, LNPIz, qi, ext=False):
    """posterior over stratum qi's (fpm, sq), everything else
    (incl. alpha) marginalized; M-STD/M-EXT structure."""
    acc = None
    for ai in range(5):
        Tq = (EXTz[ai] if ext else STDz[ai, qi])
        B = LNPIz[:, None, None, None]
        for r in range(4):
            if r == qi: continue
            B = B + lse_cells(STDz[ai, r], (2, 4), NC_STD)
        A = Tq + B[:, :, None, :, None, :]
        lw = logsumexp(A, axis=(0, 1, 3, 5))
        acc = lw if acc is None else np.logaddexp(acc, lw)
    w = np.exp(acc - logsumexp(acc))
    return w   # shape (n_fpm, n_sq)

rows = []
g1_ok = True
for seed in SEEDS:
    for law in ('simple', 'BE'):
        z = np.load(f'data/stage9f_tables_{seed}_{law}.npz')
        STDz, EXTz, LNPIz = z['STD'], z['EXT'], z['LNPI']
        # G9F-1 profile-mode regression vs 9D printed
        prof_s, prof_d, prof_e = [], [], []
        for ai in range(5):
            mq = [STDz[ai, qi].max(axis=(2, 4)) for qi in range(4)]
            mq4x = EXTz[ai].max(axis=(2, 4))
            base = mq[0] + mq[1] + mq[2]
            LN4 = LNPIz.reshape(6, 1, 1, 1)
            prof_s.append(float(np.max(base + mq[3] + LN4)))
            prof_d.append(float(np.max(base + LN4)))
            prof_e.append(float(np.max(base + mq4x + LN4)))
        ah_s = a_hat(np.array(prof_s))[0]
        ah_d = a_hat(np.array(prof_d))[0]
        ah_e = a_hat(np.array(prof_e))[0]
        tgt = g9d[(law, seed)]
        ok1 = (abs(ah_s-tgt[0]) <= 0.002 and abs(ah_d-tgt[1]) <= 0.002
               and abs(ah_e-tgt[2]) <= 0.002)
        g1_ok &= ok1
        P(f"[{law} {seed}] G9F-1 profile regression: "
          f"{ah_s:.3f}/{ah_d:.3f}/{ah_e:.3f} vs 9D "
          f"{tgt[0]:.3f}/{tgt[1]:.3f}/{tgt[2]:.3f} -> "
          f"{'PASS' if ok1 else 'FAIL'}")
        # marginal configs
        zs = marg_config(STDz, EXTz, LNPIz, [0, 1, 2, 3])
        zd = marg_config(STDz, EXTz, LNPIz, [0, 1, 2])
        ze = marg_config(STDz, EXTz, LNPIz, [0, 1, 2, 3], ext_q4=True)
        zt = marg_tied(STDz, LNPIz)
        ws_, am_s, p0_s, dn_s = summarize(zs)
        wd_, am_d, p0_d, dn_d = summarize(zd)
        we_, am_e, p0_e, dn_e = summarize(ze)
        wt_, am_t, p0_t, dn_t = summarize(zt)
        P(f"[{law} {seed}] M-STD : w(a) = "
          + "/".join(f"{x:.3f}" for x in ws_)
          + f"; a_marg = {am_s:.3f}, P(a=0) = {p0_s:.3f}, "
            f"dN = {dn_s:+.1f}")
        for qi in range(4):
            wq = stratum_noise_post(STDz, EXTz, LNPIz, qi)
            efpm = float(np.sum(wq.sum(axis=1)*FPM_GRID))
            ptop = float(wq.sum(axis=1)[-1])
            esq = float(np.sum(wq.sum(axis=0)*SQ_GRID))
            pstop = float(wq.sum(axis=0)[-1])
            P(f"    Q{qi+1} noise post: E[fpm] = {efpm:.2f}, "
              f"P(fpm=3.0) = {ptop:.2f}; E[sq] = {esq:.2f}, "
              f"P(sq=0.3) = {pstop:.2f}")
        P(f"[{law} {seed}] M-DROP: w(a) = "
          + "/".join(f"{x:.3f}" for x in wd_)
          + f"; a_marg = {am_d:.3f}, P(a=0) = {p0_d:.3f}, "
            f"dN = {dn_d:+.1f}; d(drop) = {am_d-am_s:+.3f}")
        wq4x = stratum_noise_post(STDz, EXTz, LNPIz, 3, ext=True)
        efpm4 = float(np.sum(wq4x.sum(axis=1)*FPMX_GRID))
        ptop4 = float(wq4x.sum(axis=1)[-1])
        esq4 = float(np.sum(wq4x.sum(axis=0)*SQX_GRID))
        pstop4 = float(wq4x.sum(axis=0)[-1])
        P(f"[{law} {seed}] M-EXT : a_marg = {am_e:.3f}, P(a=0) = "
          f"{p0_e:.3f}, dN = {dn_e:+.1f}; d(ext) = {am_e-am_s:+.3f}; "
          f"Q4x post: E[fpm] = {efpm4:.2f}, P(fpm=4.2) = {ptop4:.2f}, "
          f"E[sq] = {esq4:.2f}, P(sq=0.5) = {pstop4:.2f}")
        P(f"[{law} {seed}] M-TIED: a_marg = {am_t:.3f}, P(a=0) = "
          f"{p0_t:.3f}, dN = {dn_t:+.1f}")
        rows.append(dict(law=law, seed=seed, am=am_s, p0=p0_s,
                         dn=dn_s, ad=am_d))
        P("")

if not (g0_ok and g1_ok and g2_ok):
    P("GATES FAILED (G9F-0 " + ('PASS' if g0_ok else 'FAIL')
      + ", G9F-1 " + ('PASS' if g1_ok else 'FAIL')
      + ", G9F-2 " + ('PASS' if g2_ok else 'FAIL')
      + ") - STOP; DO NOT QUOTE; verdict WITHHELD")
else:
    P("GATES: G9F-0 4/4, G9F-1 4/4, G9F-2, G9F-3 4/4 - ALL PASS")
    P("")
    n1 = sum(1 for r in rows if (r['am'] <= 0.15 or r['p0'] >= 0.50))
    n2 = sum(1 for r in rows if (r['am'] >= 0.40 and r['dn'] >= 8.0
                                 and r['ad'] >= r['am']-0.15))
    n3 = sum(1 for r in rows if (r['ad']-r['am'] <= -0.25))
    P(f"F-bars: COLLAPSED-rows {n1}/4; SURVIVES-rows {n2}/4; "
      f"Q4-CARRIED-rows {n3}/4")
    if n1 >= 3:
        P("==> 9F VERDICT (locked grammar): F-COLLAPSED - at "
          "marginal grade the stratified-noise model does not keep "
          "the fitted alpha; the fitted channel's boost does not "
          "survive the honest noise read.")
        P("    CREDENCE (pre-signed map): anomaly-real 58 -> 48.")
    elif n2 >= 3:
        P("==> 9F VERDICT (locked grammar): F-SURVIVES - the "
          "marginal keeps alpha with Newton contrast and without "
          "leaning on Q4; the 9D block-grade collapse was a "
          "profile artifact.")
        P("    CREDENCE (pre-signed map): anomaly-real 58 -> 63.")
    elif n3 >= 3:
        P("==> 9F VERDICT (locked grammar): F-Q4-CARRIED - the "
          "marginal alpha survives only through the worst-RUWE "
          "quartile; the fitted channel is quality-carried at "
          "marginal grade.")
        P("    CREDENCE (pre-signed map): anomaly-real 58 -> 53.")
    else:
        P("==> 9F VERDICT (locked grammar): F-GRAY-CARRIED - rows "
          "stand as measurements.")
        P("    CREDENCE (pre-signed map): HOLD 58.")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage9f_stratmarg.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage9f_stratmarg.txt")
