"""
STAGE 3M: hierarchical fit v6 — the joint (vtilde, gamma) fit. Everything
phenomenological is replaced by measured/physical structure:
  - e-MIXTURE (from 3L's U-shaped gamma): p(e) = w_c*[near-circular, e~U(0,0.2)]
    + (1-w_c)*[power family, index a-interpolated 0.6->1.0->eta]. NO sigma_m —
    the mixture must supply the broadening physically, or alpha will inflate
    and expose it.
  - 2D likelihood: per s-bin, 20(log vtilde) x 6(gamma) histogram, multinomial.
    No independence approximation; (v, angle) correlations carry e-information.
  - Companion sector from 3K (physical continuum, f_comp in {0, 0.1}) — wobble
    also DEcorrelates gamma, so the joint fit disciplines it.
  - Two contaminants with distinct 2D signatures, both weighted by SC2 ~ s^2:
    chance alignments (flat in vtilde x flat in gamma, f_c0) and unbound
    flybys (vtilde ~ U(0.5,3) x gamma triangular on [0,30] deg, f_fly).
  - Corrected per-component noise sigv/sqrt(2) on BOTH observables (3L flag).
Batched by law via argv ('simple' or 'BE') to stay under the shell timeout.
Writes data/hier6_cube_<law>.npy + appends to data/stage3m_summary.txt.
"""
import sys
import time
import numpy as np
from astropy.io import fits

LAW = sys.argv[1] if len(sys.argv) > 1 else "simple"

src = open('calcs/stage2b_population.py').read()
ns = {}
exec(src.split('# ---------- validation')[0], ns)
run = ns['run']
GM = 4*np.pi**2; A0_CAN = 7.99e-7

def load_tab(path):
    t = np.load(path); y, b = t[0][::-1], t[1][::-1]
    lny = np.log(y); lny_u = np.linspace(lny[0], lny[-1], 512)
    return lny_u, np.interp(lny_u, lny, b)
LNY_U, TAB_S = load_tab('data/efe_boost_simple.npy')
_,     TAB_B = load_tab('data/efe_boost_be.npy')
LNY0, DLNY = LNY_U[0], LNY_U[1]-LNY_U[0]
TAB = TAB_S if LAW == "simple" else TAB_B

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
vc_d = 0.9417*np.sqrt(Mtot_d/s_d)
vt_d = (4.74047/plx[ok]*np.hypot(d['pmra1'][ok]-d['pmra2'][ok],
        d['pmdec1'][ok]-d['pmdec2'][ok]))/vc_d
# gamma (data)
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
vx_ = d['pmra2']-d['pmra1']; vy_ = d['pmdec2']-d['pmdec1']
e_vx = np.sqrt(d['pmra_error1']**2+d['pmra_error2']**2)
e_vy = np.sqrt(d['pmdec_error1']**2+d['pmdec_error2']**2)
cosg = np.abs(sx_*vx_+sy_*vy_)/np.maximum(np.hypot(sx_,sy_)*np.hypot(vx_,vy_),
                                          1e-12)
gam_all = np.degrees(np.arccos(np.clip(cosg, 0, 1)))
gam_d = gam_all[ok]
sig_c = (4.74047/plx)*np.sqrt(0.5*(e_vx**2+e_vy**2))   # per-component, km/s
SBINS = [(0.2,2),(2,6),(6,20),(20,50)]
SC2 = (np.array([0.63, 3.46, 11.0, 31.6])/31.6)**2
VE = np.logspace(np.log10(0.02), np.log10(6.0), 21)   # 20 vtilde bins
GE = np.linspace(0, 90, 7)                            # 6 gamma bins
NV, NG = 20, 6
data_2d, noise_pool = [], []
for b in SBINS:
    m = (s_d>=b[0])&(s_d<b[1])
    h,_,_ = np.histogram2d(np.clip(vt_d[m],0.021,5.9), gam_d[m],
                           bins=[VE, GE])
    data_2d.append(h.astype(float))
    noise_pool.append(sig_c[ok][m])

# flyby template in the 2D space (vtilde ~ U(0.5,3), gamma triangular [0,30])
FLY = np.zeros((NV, NG))
vcen = np.sqrt(VE[:-1]*VE[1:]); gcen = 0.5*(GE[:-1]+GE[1:])
for i in range(NV):
    for j in range(NG):
        pv = 1.0 if (0.5 <= vcen[i] <= 3.0) else 0.0
        pg = max(30.0-gcen[j], 0.0)
        FLY[i,j] = pv*pg
FLY /= max(FLY.sum(), 1e-12)
UNI = np.ones((NV, NG))/(NV*NG)

SEED = 31
rng = np.random.default_rng(SEED)
N = 500_000
u = rng.random(N); lo, hi, g = 0.15, 60.0, -0.6
a_s  = ((lo**g+u*(hi**g-lo**g))**(1/g))*1e3
u_e  = rng.random(N)
psi0 = rng.random(N)*2*np.pi
nrm  = rng.normal(size=(N,3)); nrm /= np.linalg.norm(nrm,axis=1,keepdims=True)
f_ip = np.sqrt(np.clip(1-nrm[:,0]**2,0,1))
M_s  = 0.6+1.8*rng.random(N)
uph  = rng.random(N)
xhat = np.zeros((N,3)); xhat[:,0]=1
ef = xhat-nrm*nrm[:,[0]]
ef /= np.maximum(np.linalg.norm(ef,axis=1,keepdims=True),1e-12)
e2 = np.cross(nrm,ef)
los = rng.normal(size=(N,3)); los /= np.linalg.norm(los,axis=1,keepdims=True)
u_mix = rng.random(N)                      # mixture branch selector
pick = [rng.integers(0, max(len(p_),1), N) for p_ in noise_pool]
gn1, gn2 = rng.normal(size=N), rng.normal(size=N)
# companion sector (3K model)
T_SUPP = 17.8; A_RES = 130.0
M_h = 0.5*M_s
comp = {}
for k in (1, 2):
    q = 0.1+0.9*rng.random(N)
    logP = rng.normal(5.03, 2.28, N)
    P_yr = 10**logP/365.25
    a_in = (M_h*(1+q)*P_yr**2)**(1/3)
    valid = (a_in < A_RES) & (a_in < a_s/5.0)
    v_orb = 29.78*np.sqrt(M_h*(1+q)/np.maximum(a_in,1e-3))
    S = np.minimum(1.0, P_yr/T_SUPP)
    w = (q/(1+q))*v_orb*S/4.74047*valid
    wd = rng.normal(size=(N,3)); wd /= np.linalg.norm(wd,axis=1,keepdims=True)
    uc = rng.random(N)
    comp[k] = dict(w=w, wd=wd, uc=uc, mh=q*M_h*valid)

A_GRID  = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID  = np.array([1.05, 1.3])
WC_GRID = np.array([0.0, 0.10, 0.20, 0.30])
FCOMP_GRID = np.array([0.0, 0.10])
FC0_GRID = np.array([0.0, 0.02, 0.10])
FFLY_GRID = np.array([0.0, 0.02, 0.05])

def e_of(eta, wc):
    al = np.interp(np.log10(a_s), np.log10([100,500,1000,50000]),
                   [0.6, 1.0, eta, eta])
    e_pow = 0.95*u_e**(1/(1+al))
    e_cir = 0.2*u_e
    return np.where(u_mix < wc, e_cir, e_pow)

def vp_c(e_s, tab_a):
    rp, ra = a_s*(1-e_s), a_s*(1+e_s)
    xg, wg = np.polynomial.legendre.leggauss(32)
    lo_, hi_ = np.log(rp), np.log(ra)
    lr = 0.5*(hi_-lo_)[:,None]*xg[None,:] + 0.5*(hi_+lo_)[:,None]
    r = np.exp(lr)
    gN = GM*M_s[:,None]/r**2
    bst = np.interp(np.log(gN/A0_CAN), LNY_U, tab_a, right=1.0)
    dPhi = np.sum(wg[None,:]*bst*gN*r, axis=1)*0.5*(hi_-lo_)
    return np.sqrt(np.maximum(2*dPhi/(1-(rp/ra)**2), 0))

def lnL_point(o):
    """returns (FCOMP, FC0, FFLY) lnL array from the 2D histograms"""
    s3 = o[:,0,None]*ef+o[:,1,None]*e2
    v3 = o[:,2,None]*ef+o[:,3,None]*e2
    ssky = s3-los*np.sum(s3*los,axis=1,keepdims=True)
    vsky = v3-los*np.sum(v3*los,axis=1,keepdims=True)
    smag = np.linalg.norm(ssky,axis=1)
    b1 = ssky/np.maximum(smag[:,None],1e-12)
    b2 = np.cross(los, b1)
    b2 /= np.maximum(np.linalg.norm(b2,axis=1,keepdims=True),1e-12)
    vpar = np.sum(vsky*b1,axis=1)   # AU/yr
    vper = np.sum(vsky*b2,axis=1)
    s_kau = smag/1e3
    out = np.zeros((len(FCOMP_GRID), len(FC0_GRID), len(FFLY_GRID)))
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500 or len(noise_pool[bi]) == 0: continue
        vc = 2*np.pi*np.sqrt(M_s[idx]/smag[idx])
        sg = noise_pool[bi][pick[bi][idx] % len(noise_pool[bi])]/4.74047
        for fi, fcm in enumerate(FCOMP_GRID):
            vp_n = vpar[idx].copy(); vq_n = vper[idx].copy()
            mh_tot = np.zeros(len(idx))
            for k in (1, 2):
                c = comp[k]
                act = c['uc'][idx] < fcm
                mh_tot += act*c['mh'][idx]
                vp_n += act*c['w'][idx]*c['wd'][idx,0]
                vq_n += act*c['w'][idx]*c['wd'][idx,1]
            boost = np.sqrt(1+mh_tot/M_s[idx])
            vp_n = vp_n*boost + gn1[idx]*sg
            vq_n = vq_n*boost + gn2[idx]*sg
            vmag = np.hypot(vp_n, vq_n)
            vtn = vmag/vc
            cg = np.abs(vp_n)/np.maximum(vmag, 1e-12)
            gmn = np.degrees(np.arccos(np.clip(cg, 0, 1)))
            h,_,_ = np.histogram2d(np.clip(vtn,0.021,5.9), gmn, bins=[VE, GE])
            p0 = np.maximum(h/max(h.sum(),1), 1e-5); p0 /= p0.sum()
            for ci, fc in enumerate(FC0_GRID):
                for yi, ff in enumerate(FFLY_GRID):
                    wch = min(fc*SC2[bi], 0.5); wfl = min(ff*SC2[bi], 0.5)
                    wtot = min(wch+wfl, 0.6)
                    if wch+wfl > 0:
                        mixc = (wch*UNI + wfl*FLY)/(wch+wfl)
                    else:
                        mixc = UNI
                    p = (1-wtot)*p0 + wtot*mixc
                    out[fi, ci, yi] += np.sum(data_2d[bi]*np.log(p))
    return out

cube = np.full((len(A_GRID), len(E_GRID), len(WC_GRID),
                len(FCOMP_GRID), len(FC0_GRID), len(FFLY_GRID)), np.nan)
t0 = time.time()
for ai, al in enumerate(A_GRID):
    tab_a = 1.0 + al*(TAB-1.0)
    for ei, eta in enumerate(E_GRID):
        for wi, wc in enumerate(WC_GRID):
            e_s = e_of(eta, wc)
            vp = vp_c(e_s, tab_a) if al > 0 else None
            mode = 5 if al > 0 else 1
            kw = dict(a0=A0_CAN, tab=tab_a, lny0=LNY0, dlny=DLNY, vp=vp) \
                 if al > 0 else {}
            o = run(a_s, e_s, psi0, f_ip, M_s, uph, 8, 2500, mode, **kw)
            cube[ai, ei, wi] = lnL_point(o)
    print(f"alpha={al} done, {(time.time()-t0)/60:.1f} min")
prior = -0.5*((E_GRID-1.3)/0.3)**2
cube_p = cube + prior[None,:,None,None,None,None]
np.save(f'data/hier6_cube_{LAW}.npy', cube_p)

best = np.unravel_index(np.nanargmax(cube_p), cube_p.shape)
prof = np.nanmax(cube_p, axis=(1,2,3,4,5))
rel = prof - np.nanmax(prof)
imax = int(np.nanargmax(prof))
wc_prof = np.nanmax(cube_p, axis=(0,1,3,4,5))
lines = [
    f"{LAW}: {(time.time()-t0)/60:.1f} min (joint 2D vtilde x gamma fit)",
    f"  best: alpha={A_GRID[best[0]]}, eta={E_GRID[best[1]]}, "
    f"w_circ={WC_GRID[best[2]]}, f_comp={FCOMP_GRID[best[3]]}, "
    f"f_c0={FC0_GRID[best[4]]}, f_fly={FFLY_GRID[best[5]]}",
    f"  best total lnL: {np.nanmax(cube_p):.1f}",
    f"  alpha profile: {np.round(rel,1).tolist()}  "
    f"(interior: {0 < imax < len(A_GRID)-1})",
    f"  Newton dlnL: {np.nanmax(prof)-prof[0]:+.1f}",
    f"  w_circ profile: {np.round(wc_prof-np.nanmax(wc_prof),1).tolist()}",
]
for l in lines: print(l)
with open('data/stage3m_summary.txt', 'a') as f:
    f.write("\n".join(lines)+"\n")
print(f"saved: data/hier6_cube_{LAW}.npy; appended data/stage3m_summary.txt")
