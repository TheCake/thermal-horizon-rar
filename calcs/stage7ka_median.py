"""Stage 7K-a — THE FORWARD MEDIAN (pre-reg in NOTES 2026-07-27,
committed before execution; G0 baseline-noise clause amended pre-run).

Question: how much of the model-free 2C anchor (corrected boost 1.078,
CI 1.052-1.103; ratio of median vtilde between 6-30 and 0.2-2 kAU;
Newton predicts 1.000 exactly) does the landed width-complete model's
NEWTON-BEST cell produce through its absorbers (companions incl. kw,
noise fpm, smear sq)?

Legs per seed (31 primary, 101 stability):
  G0-i   noise-off identity        R in [0.995, 1.005]
  G0-ii  baseline noise (fpm 1.2)  R in [0.995, 1.015] (~0.005-0.010
         lift predicted from the sg0*sqrt(s) growth)
  R_N    NEWTON-BEST cell (argmax of the alpha=0 slice of cube+LNPI)
  R_A    sanity: each law's PROF argmax cell (expect [1.04, 1.12])
Bars (s-reweighted R_N, seed 31): >= 1.052 MEDIAN-ABSORBED;
<= 1.030 MEDIAN-SURVIVES; else GRAY.  Seed 101 must not flip the
category, else AMBIG-quoted.
Output: data/stage7ka_median.txt
"""
import numpy as np, os, time
from astropy.io import fits

OUT = 'data/stage7ka_median.txt'
open(OUT, 'w').close()
def P(s):
    print(s, flush=True)
    with open(OUT, 'a') as f:
        f.write(s+"\n")

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

# --- catalog (verbatim marginal-pipeline loading, corrected velocities) ---
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
sig_c = ((4.74047/plx)*np.sqrt(0.5*(e_vx**2+e_vy**2)))
sig_ok = sig_c[ok]

# the 2C statistic on the corrected data (reproduction check)
B2C = [(0.2, 2.0), (6.0, 30.0)]
mn = (s_d >= 0.2) & (s_d < 2.0)
mw = (s_d >= 6.0) & (s_d < 30.0)
R_data = float(np.median(vt_d[mw])/np.median(vt_d[mn]))
P(f"STAGE 7K-a: data corrected-velocity 2C ratio = {R_data:.3f} "
  f"(anchor 1.078, CI 1.052-1.103; N = {mn.sum()}/{mw.sum()})")

# data s-distributions for reweighting (10 log sub-bins per 2C bin)
SUBS = {}
for b in B2C:
    edges = np.logspace(np.log10(b[0]), np.log10(b[1]), 11)
    m = (s_d >= b[0]) & (s_d < b[1])
    h, _ = np.histogram(s_d[m], bins=edges)
    SUBS[b] = (edges, h/max(h.sum(), 1))

# per-2C-bin noise pools (the marginal pipeline's convention, 2C bins)
pool2c = {b: sig_ok[(s_d >= b[0]) & (s_d < b[1])] for b in B2C}

# --- grids + landed anchor (photow3 convention) ---------------------------
N = 500_000
A_GRID  = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID  = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FC0_GRID = np.array([0.10])
FFLY_GRID = np.array([0.05, 0.10])
FPM_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW_GRID = np.array([0.7, 1.0, 1.4])
SQ_GRID = np.array([0.0, 0.1, 0.2, 0.3])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
assert str(pz['version']).startswith('v2c')
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
gband = (pz['conv_band']/0.30 if 'conv_band' in pz.files
         else np.array([0.33, 1.30]))
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
    p['pk2c'] = {b: rng.integers(0, max(len(pool2c[b]), 1), N) for b in B2C}
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
        MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        MGs = np.interp(-np.clip(q*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        l_ = 10**(-0.4*(MGs-MGp))
        wfac = np.abs(q/(1+q) - l_/(1+l_))
        w = wfac*v_orb*S/4.74047*valid
        wd = rng.normal(size=(N,3)); wd /= np.linalg.norm(wd,axis=1,keepdims=True)
        p['comp'][k] = dict(w=w, wd=wd, uc=rng.random(N), mh=q*M_h*valid)
    p['gs'] = rng.normal(size=N)
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

def ratio_of(p, o, fcm, kw, fpm, sq, noise=True):
    """The 2C median ratio for a forward sky at the given absorber cell.
    Returns (natural, s-reweighted)."""
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
    meds_nat, meds_rw = {}, {}
    for b in B2C:
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        pool = pool2c[b]
        sg0 = pool[p['pk2c'][b][idx] % len(pool)]/4.74047
        cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
        mh_tot = np.zeros(len(idx))
        for k in (1, 2):
            c = p['comp'][k]
            act = c['uc'][idx] < fcm
            mh_tot += act*c['mh'][idx]
            cvp += act*c['w'][idx]*c['wd'][idx,0]
            cvq += act*c['w'][idx]*c['wd'][idx,1]
        boost = np.sqrt(1+mh_tot/p['M_s'][idx])
        vp_b = vpar[idx] + kw*cvp
        vq_b = vper[idx] + kw*cvq
        if noise:
            vp_n = vp_b*boost + p['gn1'][idx]*sg0*fpm
            vq_n = vq_b*boost + p['gn2'][idx]*sg0*fpm
        else:
            vp_n = vp_b*boost
            vq_n = vq_b*boost
        vmag = np.hypot(vp_n, vq_n)
        keep = vmag*4.74047 <= (2.978/np.sqrt(s_kau[idx])
                                + 2.8284*sg0*4.74047)
        vtn = (vmag/vc)[keep]
        if sq != 0.0:
            vtn = vtn*np.exp(sq*p['gs'][idx][keep])
        sk = s_kau[idx][keep]
        meds_nat[b] = float(np.median(vtn))
        edges, dfr = SUBS[b]
        mi = np.clip(np.digitize(sk, edges)-1, 0, 9)
        mh, _ = np.histogram(sk, bins=edges)
        wgt = dfr[mi]/np.maximum((mh/max(mh.sum(),1))[mi], 1e-9)
        o_ = np.argsort(vtn)
        cw = np.cumsum(wgt[o_]); cw /= cw[-1]
        meds_rw[b] = float(np.interp(0.5, cw, vtn[o_]))
    return (meds_nat[B2C[1]]/meds_nat[B2C[0]],
            meds_rw[B2C[1]]/meds_rw[B2C[0]])

t0 = time.time()
for seed in (31, 101):
    P(f"\n== seed {seed} ==")
    p = build_pop(seed)
    # cells from the photow3 cubes + LANDED-CONV
    cells = {}
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        cw = np.load(f'data/stage7j_cube_full_photow3_{seed}_{law}.npy')
        cb = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1)) \
             + LNPI.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
        if law == 'simple':
            bx = np.unravel_index(np.nanargmax(cb[0]), cb[0].shape)
            cells['newton'] = (0.0, TAB) + bx
        bxa = np.unravel_index(np.nanargmax(cb), cb.shape)
        cells[f'prof-{law}'] = (A_GRID[bxa[0]], TAB) + bxa[1:]
    for tag, cell in cells.items():
        al, TAB, ei, wi, fi, ci, yi, pi, ki, si = cell
        eta, wr = E_GRID[ei], WR_GRID[wi]
        fcm, fpm = FCOMP_GRID[fi], FPM_GRID[pi]
        kw, sq = KW_GRID[ki], SQ_GRID[si]
        e_s = e_of(p, eta, wr)
        if al > 0:
            tab_a = 1.0 + al*(TAB-1.0)
            vp = vp_c(p, e_s, tab_a)
            o = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'],
                    p['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_a,
                    lny0=LNY0, dlny=DLNY, vp=vp)
        else:
            o = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'],
                    p['uph'], 8, 2500, 1)
        if tag == 'newton':
            g0i = ratio_of(p, o, 0.0, 1.0, 1.2, 0.0, noise=False)
            g0ii = ratio_of(p, o, 0.0, 1.0, 1.2, 0.0, noise=True)
            P(f"G0-i  noise-off identity: R = {g0i[0]:.3f} nat / "
              f"{g0i[1]:.3f} rw -> "
              f"{'PASS' if 0.995 <= g0i[1] <= 1.005 else 'FAIL'} "
              f"[0.995, 1.005]")
            P(f"G0-ii baseline noise:     R = {g0ii[0]:.3f} nat / "
              f"{g0ii[1]:.3f} rw -> "
              f"{'PASS' if 0.995 <= g0ii[1] <= 1.015 else 'FAIL'} "
              f"[0.995, 1.015] (noise contribution "
              f"{g0ii[1]-g0i[1]:+.3f})")
        rn, rw = ratio_of(p, o, fcm, kw, fpm, sq)
        P(f"{tag:<12} alpha={al} eta={eta} wr={wr} fcomp={fcm} "
          f"fpm={fpm} kw={kw} sq={sq}: R = {rn:.3f} nat / "
          f"{rw:.3f} reweighted")
        if tag == 'newton':
            cat = ('MEDIAN-ABSORBED' if rw >= 1.052 else
                   'MEDIAN-SURVIVES' if rw <= 1.030 else 'GRAY')
            P(f"  ==> R_N (reweighted) = {rw:.3f} vs anchor 1.078 "
              f"(CI 1.052-1.103): {cat}"
              + (f"; unabsorbed excess = {1.078-rw:+.3f}"
                 if rw < 1.052 else ""))

P(f"\ndone ({(time.time()-t0)/60:.1f} min)")
