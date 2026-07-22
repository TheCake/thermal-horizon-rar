"""
STAGE 3S: the g_ext scan (TODO #2f) â€” the v7 fit (seed 31) repeated with EFE
boost tables solved at g_ext/a0 in {1.4, 1.6, 1.9, 2.2, 2.4} (tables from
stage3r; the 1.9 tables are the originals). Question: does the v7-model
alpha = 1.54 +/- 0.13 (BE) fall toward the parameter-free 1.0 when the
assumed external field weakens â€” i.e., is the alpha>1 tension an EFE
calibration artifact? argv: g_ext tag ('1p4','1p6','1p9','2p2','2p4').
Appends to data/stage3s_summary.txt.
"""
import sys
import time
import numpy as np
from astropy.io import fits

GTAG = sys.argv[1] if len(sys.argv) > 1 else '1p9'
SEEDS = [31]

src = open('calcs/stage2b_population.py').read()
ns = {}
exec(src.split('# ---------- validation')[0], ns)
run = ns['run']
GM = 4*np.pi**2; A0_CAN = 7.99e-7

def load_tab(path):
    t = np.load(path); y, b = t[0][::-1], t[1][::-1]
    lny = np.log(y); lny_u = np.linspace(lny[0], lny[-1], 512)
    return lny_u, np.interp(lny_u, lny, b)
if GTAG == '1p9':
    PS, PB = 'data/efe_boost_simple.npy', 'data/efe_boost_be.npy'
else:
    PS = f'data/efe_boost_simple_g{GTAG}.npy'
    PB = f'data/efe_boost_be_g{GTAG}.npy'
LNY_U, TAB_S = load_tab(PS)
_,     TAB_B = load_tab(PB)
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
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot_d = np.interp(MG1[ok],MG_T,MS_T)+np.interp(MG2[ok],MG_T,MS_T)
s_d = sep[ok]/1e3
vt_d = (4.74047/plx[ok]*np.hypot(d['pmra1'][ok]-d['pmra2'][ok],
        d['pmdec1'][ok]-d['pmdec2'][ok]))/(0.9417*np.sqrt(Mtot_d/s_d))
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
vx_ = d['pmra2']-d['pmra1']; vy_ = d['pmdec2']-d['pmdec1']
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
WR_GRID = np.array([0.10, 0.20, 0.30])
FPM = 1.5
FCOMP_GRID = np.array([0.0, 0.10])
FC0_GRID = np.array([0.10])
FFLY_GRID = np.array([0.05, 0.10])

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
        w = (q/(1+q))*v_orb*S/4.74047*valid
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
    out = np.zeros((len(FCOMP_GRID), len(FC0_GRID), len(FFLY_GRID)))
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500 or len(noise_pool[bi]) == 0: continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = noise_pool[bi][p['pick'][bi][idx] % len(noise_pool[bi])]/4.74047
        for fi, fcm in enumerate(FCOMP_GRID):
            vp_b = vpar[idx].copy(); vq_b = vper[idx].copy()
            mh_tot = np.zeros(len(idx))
            for k in (1, 2):
                c = p['comp'][k]
                act = c['uc'][idx] < fcm
                mh_tot += act*c['mh'][idx]
                vp_b += act*c['w'][idx]*c['wd'][idx,0]
                vq_b += act*c['w'][idx]*c['wd'][idx,1]
            boost = np.sqrt(1+mh_tot/p['M_s'][idx])
            vp_n = vp_b*boost + p['gn1'][idx]*sg0*FPM
            vq_n = vq_b*boost + p['gn2'][idx]*sg0*FPM
            vmag = np.hypot(vp_n, vq_n)
            keep = vmag*4.74047 <= (2.978/np.sqrt(s_kau[idx])
                                    + 2.8284*sg0*4.74047)
            vtn = (vmag/vc)[keep]
            gmn = np.degrees(np.arccos(np.clip(
                np.abs(vp_n[keep])/np.maximum(vmag[keep],1e-12), 0, 1)))
            h,_,_ = np.histogram2d(np.clip(vtn,0.021,5.9), gmn, bins=[VE, GE])
            p0 = np.maximum(h/max(h.sum(),1), 1e-5); p0 /= p0.sum()
            for ci, fc in enumerate(FC0_GRID):
                for yi, ff in enumerate(FFLY_GRID):
                    wch = min(fc*SC2[bi], 0.5); wfl = min(ff*SC2[bi], 0.5)
                    wtot = min(wch+wfl, 0.6)
                    mixc = (wch*UNI_B[bi] + wfl*FLY_B[bi])/(wch+wfl)
                    pp = (1-wtot)*p0 + wtot*mixc
                    out[fi, ci, yi] += np.sum(data_2d[bi]*np.log(pp))
    return out

def P(s):
    print(s)
    with open('data/stage3s_summary.txt', 'a') as f:
        f.write(s+"\n")

prior = -0.5*((E_GRID-1.3)/0.3)**2
P(f"STAGE 3S g_ext={GTAG} (seed 31): a={A_GRID.tolist()}, eta={E_GRID.tolist()}, "
  f"wr={WR_GRID.tolist()}, fpm={FPM}, fcomp={FCOMP_GRID.tolist()}, "
  f"fc0={FC0_GRID.tolist()}, ffly={FFLY_GRID.tolist()}")
for seed in SEEDS:
    t0 = time.time()
    p = build_pop(seed)
    best_lnl = {}
    for law, TAB in (("simple", TAB_S), ("BE", TAB_B)):
        cube = np.full((len(A_GRID), len(E_GRID), len(WR_GRID),
                        len(FCOMP_GRID), len(FC0_GRID), len(FFLY_GRID)),
                       np.nan)
        for ai, al in enumerate(A_GRID):
            tab_a = 1.0 + al*(TAB-1.0)
            for ei, eta in enumerate(E_GRID):
                for wi, wr in enumerate(WR_GRID):
                    if al == 0.0 and law == "BE" and (ei, wi) in best_lnl.get(
                            '_newt', {}):
                        cube[ai, ei, wi] = best_lnl['_newt'][(ei, wi)]
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
                        best_lnl.setdefault('_newt', {})[(ei, wi)] = \
                            cube[ai, ei, wi]
        cube_p = cube + prior[None,:,None,None,None,None]
        prof = np.nanmax(cube_p, axis=(1,2,3,4,5))
        imax = int(np.nanargmax(prof))
        ahat = A_GRID[imax]
        if 0 < imax < len(A_GRID)-1:
            x = A_GRID[imax-1:imax+2]; y = prof[imax-1:imax+2]
            c2, c1, _ = np.polyfit(x, y, 2)
            if c2 < 0: ahat = -c1/(2*c2)
        best = np.unravel_index(np.nanargmax(cube_p), cube_p.shape)
        best_lnl[law] = np.nanmax(cube_p)
        P(f"g={GTAG} {law}: a_hat={ahat:.2f} "
          f"(grid {A_GRID[imax]}, interior={0<imax<len(A_GRID)-1}), "
          f"dlnL(Newton)={np.nanmax(prof)-prof[0]:+.1f}, "
          f"wr={WR_GRID[best[2]]}, prof={np.round(prof-prof.max(),1).tolist()}")
    P(f"  g={GTAG}: BE-minus-simple best lnL = "
      f"{best_lnl['BE']-best_lnl['simple']:+.1f}  "
      f"({(time.time()-t0)/60:.1f} min)")
print("\nbatch done; appended data/stage3s_summary.txt")

