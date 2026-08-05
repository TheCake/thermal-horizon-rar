"""STAGE 9A — THE STRATIFIED-NOISE ALPHA RE-FIT (the alpha-exposure
decider).  Pre-registered BEFORE any run.

8Z: the global sq = 0.2 is partly an aggregation artifact (one
error-scale over quality-heterogeneous strata manufactures ~2x
width).  This stage measures what the repair does to alpha.  RUWE
quartiles; per alpha in {0, 0.5, 1, 1.5, 2} (eta, wr frozen at MAP
— BLOCK grade, profile + LNPI; NOT comparable to the 7J marginal
numbers; direction-measurement): per-stratum blocks (fcomp, ffly,
fpm, kw, sq, ws in {0, 0.045}); TIED = Sum_q blocks maximized with
SHARED (fpm, sq); FREE = Sum_q max over per-stratum (fpm_q, sq_q)
— nested, same stratified data; ws and population axes shared (the
floor is global per 8T/8Z).

Gates: G9A-1 unsplit block at MAP alpha, ws=0 vs the lker cube =
0.00e+00 (bit); G9A-2 stratum count conservation; G9A-3 8Z
regression — sqbar recomputed from the alpha=MAP stratum blocks vs
8Z's printed values (16 checks, bar 0.0005).
Bars (locked): A1 EXPOSURE-MATERIAL iff |d alpha| >= 0.25 in >=
3/4 law-seeds (successor = full marginal re-run).  A2
EXPOSURE-CONTAINED iff |d alpha| <= 0.15 in >= 3/4.  A3
GRAY-CARRIED else.  Co-reads: DlnL(free-tied) at alpha-hat,
per-stratum (fpm_q, sq_q), Newton contrast under both models, edge
flags.  NO credence movement (measurement round; pre-stated).
Output: data/stage9a_stratalpha.txt
"""
import re
import numpy as np, time
from astropy.io import fits

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

def eval_block_s(p, prj, D2, PLs, UB, FB):
    smag, vpar, vper = prj
    s_kau = smag/1e3
    out = np.zeros((len(FCOMP_GRID), len(FFLY_GRID), len(FPM_GRID),
                    len(KW_GRID), len(SQ_GRID), len(WS2_GRID)))
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
g8z = {}
for ln in open('data/stage8z_dose.txt').read().splitlines():
    mm = re.match(r"\[(simple|BE) (\d+)\] ruwe-Q(\d): sqbar=([\d.]+)",
                  ln)
    if mm:
        g8z[(mm.group(1), int(mm.group(2)), int(mm.group(3)))] = \
            float(mm.group(4))
assert len(g8z) == 16, len(g8z)

okv = ok
ruwe = np.maximum(_pick('ruwe1', 'RUWE1'),
                  _pick('ruwe2', 'RUWE2'))[okv]
qs_ = np.percentile(ruwe, [25, 50, 75])
QMASKS = [ruwe <= qs_[0],
          (ruwe > qs_[0]) & (ruwe <= qs_[1]),
          (ruwe > qs_[1]) & (ruwe <= qs_[2]),
          ruwe > qs_[2]]

t0 = time.time()
P("9A THE STRATIFIED-NOISE ALPHA RE-FIT (pre-reg committed BEFORE "
  "any run; measurement round; NO credence movement)")
g2_ok = True
for bi, b in enumerate(SBINS):
    m_all = (s_d>=b[0])&(s_d<b[1])
    g2_ok &= (sum(int((m_all&mk).sum()) for mk in QMASKS)
              == int(m_all.sum()))
P(f"ruwe quartile edges = " + "/".join(f"{v:.3f}" for v in qs_)
  + f"; G9A-2 count conservation -> {'PASS' if g2_ok else 'FAIL'}")
STRATA = [build_stratum(mk) for mk in QMASKS]
P("")

g1_ok = g3_ok = True
rows = []
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
        tied_prof, free_prof = [], []
        q_at = {}
        for ai, al in enumerate(A_GRID):
            tab_a = 1.0 + al*(TAB-1.0)
            vp_f = vp_c(pf, e_f, tab_a)
            o_f = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'],
                      pf['M_s'], pf['uph'], 8, 2500, 5, a0=A0_CAN,
                      tab=tab_a, lny0=LNY0, dlny=DLNY, vp=vp_f)
            prj = project(pf, o_f)
            bs_q = [eval_block_s(pf, prj, *st[:4]) for st in STRATA]
            if abs(al - al_map) < 1e-9:
                blk_u = eval_block_s(pf, prj, data_2d, noise_pool,
                                     UNI_B, FLY_B)
                sl = np.squeeze(c9[bx[0], bx[1], bx[2]], axis=1)
                dmax = float(np.max(np.abs(blk_u[..., 0] - sl)))
                ok1 = dmax <= 1e-9
                g1_ok &= ok1
                P(f"[{law} {seed}] G9A-1 cube identity at "
                  f"alpha={al}: max|d| = {dmax:.2e} -> "
                  f"{'PASS' if ok1 else 'FAIL'}")
                for qi, bs in enumerate(bs_q):
                    ss = bs + LNPI.reshape(6, 1, 1, 1, 1, 1)
                    w_ = np.exp(ss - np.max(ss)); w_ /= w_.sum()
                    sqb = float(np.tensordot(
                        w_.sum(axis=(0, 1, 2, 3, 5)), SQ_GRID, 1))
                    ok3 = abs(sqb - g8z[(law, seed, qi+1)]) <= 5e-4
                    g3_ok &= ok3
                    if not ok3:
                        P(f"[{law} {seed}] G9A-3 Q{qi+1} sqbar "
                          f"{sqb:.4f} vs 8Z "
                          f"{g8z[(law, seed, qi+1)]:.3f} -> FAIL")
            # tied: shared (fpm, sq) across strata
            bsum = bs_q[0] + bs_q[1] + bs_q[2] + bs_q[3]
            tied = float(np.max(bsum
                                + LNPI.reshape(6, 1, 1, 1, 1, 1)))
            # free: per-stratum (fpm, sq) profiled, shared cell
            mq = [b_.max(axis=(2, 4)) for b_ in bs_q]
            F = mq[0] + mq[1] + mq[2] + mq[3]
            Fp = F + LNPI.reshape(6, 1, 1, 1)
            free = float(np.max(Fp))
            tied_prof.append(tied); free_prof.append(free)
            cix = np.unravel_index(np.argmax(Fp), Fp.shape)
            q_at[ai] = (cix, [np.unravel_index(
                np.argmax(b_[cix[0], cix[1], :, cix[2], :, cix[3]]),
                (6, 4)) for b_ in bs_q])
        tied_prof = np.array(tied_prof)
        free_prof = np.array(free_prof)

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

        ah_t, ed_t, it_ = a_hat(tied_prof)
        ah_f, ed_f, if_ = a_hat(free_prof)
        dNt = float(tied_prof.max() - tied_prof[0])
        dNf = float(free_prof.max() - free_prof[0])
        dLn = float(free_prof[if_] - tied_prof[if_])
        cix, qarg = q_at[if_]
        sqq = [float(SQ_GRID[j[1]]) for j in qarg]
        fpq = [float(FPM_GRID[j[0]]) for j in qarg]
        P(f"[{law} {seed}] tied prof: "
          + "/".join(f"{v:.1f}" for v in tied_prof - tied_prof.max())
          + f" -> a_hat = {ah_t:.3f}"
          + (" EDGE" if ed_t else ""))
        P(f"[{law} {seed}] free prof: "
          + "/".join(f"{v:.1f}" for v in free_prof - free_prof.max())
          + f" -> a_hat = {ah_f:.3f}"
          + (" EDGE" if ed_f else ""))
        P(f"[{law} {seed}] d_alpha = {ah_f-ah_t:+.3f}; "
          f"DlnL(free-tied)@a_hat = {dLn:+.1f}; Newton DlnL "
          f"tied/free = {dNt:+.1f}/{dNf:+.1f}; per-stratum at free "
          f"a_hat: sq = " + "/".join(f"{v:.1f}" for v in sqq)
          + "; fpm = " + "/".join(f"{v:.1f}" for v in fpq))
        rows.append(dict(law=law, seed=seed, da=ah_f-ah_t,
                         ah_t=ah_t, ah_f=ah_f, dLn=dLn,
                         dNt=dNt, dNf=dNf))
        P("")

if not (g1_ok and g2_ok and g3_ok):
    P("GATES FAILED (G9A-1 " + ('PASS' if g1_ok else 'FAIL')
      + ", G9A-2 " + ('PASS' if g2_ok else 'FAIL')
      + ", G9A-3 " + ('PASS' if g3_ok else 'FAIL')
      + ") - STOP; DO NOT QUOTE; verdict WITHHELD")
else:
    P("GATES: G9A-1 4/4, G9A-2 PASS, G9A-3 16/16 - ALL PASS")
    P("")
    n_a1 = sum(1 for x in rows if abs(x['da']) >= 0.25)
    n_a2 = sum(1 for x in rows if abs(x['da']) <= 0.15)
    P(f"A-bars: A1-rows {n_a1}/4 (|d_alpha| >= 0.25); A2-rows "
      f"{n_a2}/4 (|d_alpha| <= 0.15)")
    if n_a1 >= 3:
        P("==> 9A VERDICT (locked grammar): EXPOSURE-MATERIAL - "
          "the honest stratified width model materially moves "
          "alpha; successor = the full marginal re-run with "
          "stratified noise.")
    elif n_a2 >= 3:
        P("==> 9A VERDICT (locked grammar): EXPOSURE-CONTAINED - "
          "the alpha measurement is robust to the width-model "
          "repair at block grade.")
    else:
        P("==> 9A VERDICT (locked grammar): GRAY-CARRIED - rows "
          "stand as measurements.")
    P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage9a_stratalpha.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage9a_stratalpha.txt")
