"""
STAGE 3C: hierarchical fit v3 — the ridge-collapse run.
Adds to v2: (1) chance-alignment/contamination component, flat in vtilde-bins,
weight w_b = f_c0*(s_b/31.6 kAU)^2 (physically: alignments scale with sky area,
so they live in the wide bins); (2) Gaussian prior eta ~ N(1.3, 0.3) from
Hwang+22; (3) denser eta grid near the prior. Populations re-realized per
(alpha, eta) as before. Writes cubes + a human-readable summary to data/.
"""
import time
import numpy as np
from astropy.io import fits

src = open('calcs/stage2b_population.py').read()
ns = {}
exec(src.split('# ---------- validation')[0], ns)
run = ns['run']
GM = 4*np.pi**2; A0_CAN = 7.99e-7; GEXT = 1.9*A0_CAN

def load_tab(path):
    t = np.load(path); y, b = t[0][::-1], t[1][::-1]
    lny = np.log(y); lny_u = np.linspace(lny[0], lny[-1], 512)
    return lny_u, np.interp(lny_u, lny, b)
LNY_U, TAB_S = load_tab('data/efe_boost_simple.npy')
_,     TAB_B = load_tab('data/efe_boost_be.npy')
LNY0, DLNY = LNY_U[0], LNY_U[1]-LNY_U[0]

d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1,plx2 = d['parallax1'],d['parallax2']
eplx1,eplx2 = d['parallax_error1'],d['parallax_error2']
sep,Rch = d['sep_AU'],d['R_chance_align']
G1,G2 = d['phot_g_mean_mag1'],d['phot_g_mean_mag2']
plx = 0.5*(plx1+plx2)
MG1 = G1+5*np.log10(np.maximum(plx1,1e-6))-10
MG2 = G2+5*np.log10(np.maximum(plx2,1e-6))-10
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
srel_d = sigv[ok]/vc_d
SBINS = [(0.2,2),(2,6),(6,20),(20,50)]
SC2 = (np.array([0.63, 3.46, 11.0, 31.6])/31.6)**2     # contamination shape
VEDGE = np.logspace(np.log10(0.02), np.log10(6.0), 41)
data_counts, data_srel = [], []
for b in SBINS:
    m = (s_d>=b[0])&(s_d<b[1])
    c,_ = np.histogram(np.clip(vt_d[m], 0.021, 5.9), bins=VEDGE)
    data_counts.append(c.astype(float)); data_srel.append(srel_d[m])

rng = np.random.default_rng(31)
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
def sky(v3): return v3-los*np.sum(v3*los,axis=1,keepdims=True)
q_in  = 0.1+0.9*rng.random(N)
a_in  = 10**(np.log10(2)+rng.random(N)*np.log10(50))
wdir  = rng.normal(size=(N,3)); wdir /= np.linalg.norm(wdir,axis=1,keepdims=True)
v_wob = (q_in/(1+q_in))*29.78*np.sqrt(0.5*M_s*(1+q_in)/a_in)/4.74047
m_hid = q_in*0.5*M_s
is3_u = rng.random(N)
noise_pick = [rng.integers(0, max(len(s),1), N) for s in data_srel]
gauss1, gauss2 = rng.normal(size=N), rng.normal(size=N)

def e_of(eta):
    al = np.interp(np.log10(a_s), np.log10([100,500,1000,50000]),
                   [0.6, 1.0, eta, eta])
    return 0.95*u_e**(1/(1+al))

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

A_GRID  = np.array([0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0])
E_GRID  = np.array([0.8,1.0,1.15,1.3,1.45,1.6,1.9])
FT_GRID = np.array([0.0,0.05,0.10,0.15,0.20])
FC_GRID = np.array([0.0,0.01,0.02,0.05,0.10])

def lnL_point(o):
    s3 = o[:,0,None]*ef+o[:,1,None]*e2
    v3 = o[:,2,None]*ef+o[:,3,None]*e2
    sp = np.linalg.norm(sky(s3),axis=1)
    vt0 = np.linalg.norm(sky(v3),axis=1)
    s_kau = sp/1e3
    out = np.zeros((len(FT_GRID), len(FC_GRID)))
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500 or len(data_srel[bi]) == 0: continue
        vc = 2*np.pi*np.sqrt(M_s[idx]/sp[idx])
        sr = data_srel[bi][noise_pick[bi][idx] % len(data_srel[bi])]
        for ti, ft in enumerate(FT_GRID):
            trip = is3_u[idx] < ft
            vv = vt0[idx]*np.where(trip, np.sqrt(1+m_hid[idx]/M_s[idx]), 1.0)
            wob = np.where(trip, v_wob[idx], 0.0)
            vtn = np.hypot(vv/vc + gauss1[idx]*sr + wob*np.abs(wdir[idx,0])/vc,
                           gauss2[idx]*sr + wob*np.abs(wdir[idx,1])/vc)
            h,_ = np.histogram(np.clip(vtn, 0.021, 5.9), bins=VEDGE)
            p_orb = np.maximum(h/max(h.sum(),1), 1e-5); p_orb /= p_orb.sum()
            for ci, fc in enumerate(FC_GRID):
                w = min(fc*SC2[bi], 0.5)
                p = (1-w)*p_orb + w/40.0
                out[ti, ci] += np.sum(data_counts[bi]*np.log(p))
    return out

summary = []
results = {}
newton_cache = {}
for law, TAB in (("simple", TAB_S), ("BE", TAB_B)):
    cube = np.full((len(A_GRID), len(E_GRID), len(FT_GRID), len(FC_GRID)), np.nan)
    t0 = time.time()
    for ai, al in enumerate(A_GRID):
        tab_a = 1.0 + al*(TAB-1.0)
        for ei, eta in enumerate(E_GRID):
            if al == 0.0:
                if ei in newton_cache: cube[ai, ei] = newton_cache[ei]; continue
            e_s = e_of(eta)
            vp = vp_c(e_s, tab_a) if al > 0 else None
            mode = 5 if al > 0 else 1
            kw = dict(a0=A0_CAN, tab=tab_a, lny0=LNY0, dlny=DLNY, vp=vp) \
                 if al > 0 else {}
            o = run(a_s, e_s, psi0, f_ip, M_s, uph, 8, 2500, mode, **kw)
            val = lnL_point(o)
            cube[ai, ei] = val
            if al == 0.0: newton_cache[ei] = val
    # eta prior
    prior = -0.5*((E_GRID-1.3)/0.3)**2
    cube_p = cube + prior[None,:,None,None]
    results[law] = cube_p
    prof = np.nanmax(cube_p, axis=(1,2,3))
    best = np.unravel_index(np.nanargmax(cube_p), cube_p.shape)
    mx = np.nanmax(prof)
    # 1-sigma interval on alpha via profile crossing (dlnL = -0.5), linear interp
    rel = prof - mx
    inside = A_GRID[rel >= -0.5]
    lines = [
        f"{law}: {(time.time()-t0)/60:.1f} min",
        f"  best: alpha={A_GRID[best[0]]}, eta={E_GRID[best[1]]}, "
        f"f_t={FT_GRID[best[2]]}, f_c0={FC_GRID[best[3]]}",
        f"  alpha profile (lnL-max): "
        f"{np.round(rel,1).tolist()}",
        f"  alpha 1-sigma grid range: [{inside.min()}, {inside.max()}]"
        f"  (interior max: {A_GRID[best[0]] < A_GRID[-1]})",
        f"  Newton rejection dlnL: {mx - prof[0]:+.1f}",
    ]
    for l in lines: print(l)
    summary += lines
np.save('data/hier3_cube_simple.npy', results['simple'])
np.save('data/hier3_cube_BE.npy', results['BE'])
with open('data/stage3c_summary.txt', 'w') as f:
    f.write("STAGE 3C v3 fit (contamination + eta prior), grids:\n"
            f"alpha={A_GRID.tolist()}\neta={E_GRID.tolist()}\n"
            f"f_t={FT_GRID.tolist()}\nf_c0={FC_GRID.tolist()}\n\n"
            + "\n".join(summary) + "\n")
print("\nsaved: data/hier3_cube_*.npy, data/stage3c_summary.txt")
