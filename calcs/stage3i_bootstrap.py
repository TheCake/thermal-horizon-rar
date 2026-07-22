"""
STAGE 3I: data-bootstrap half of the error budget (TODO #3). One reduced-grid
model evaluation (seed 31) with the per-grid-point histogram probabilities
SAVED, then 1000 bootstrap replicates of the 14,071 data pairs re-scored
against the saved model — pure arithmetic, no orbit re-runs. Approximation:
the model's noise draws stay tied to the original per-bin noise distributions
(second-order). Output: bootstrap scatter of alpha-hat and of the Newton dlnL,
to be combined in quadrature-of-judgment with the Stage-3H realization scatter.
Writes data/stage3i_summary.txt.
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
SC2 = (np.array([0.63, 3.46, 11.0, 31.6])/31.6)**2
VEDGE = np.logspace(np.log10(0.02), np.log10(6.0), 41)
data_counts, data_srel = [], []
for b in SBINS:
    m = (s_d>=b[0])&(s_d<b[1])
    c,_ = np.histogram(np.clip(vt_d[m], 0.021, 5.9), bins=VEDGE)
    data_counts.append(c.astype(float)); data_srel.append(srel_d[m])

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
def sky(v3): return v3-los*np.sum(v3*los,axis=1,keepdims=True)
q_in  = 0.1+0.9*rng.random(N)
a_in  = 10**(np.log10(2)+rng.random(N)*np.log10(50))
wdir  = rng.normal(size=(N,3)); wdir /= np.linalg.norm(wdir,axis=1,keepdims=True)
v_wob = (q_in/(1+q_in))*29.78*np.sqrt(0.5*M_s*(1+q_in)/a_in)/4.74047
m_hid = q_in*0.5*M_s
is3_u = rng.random(N)
noise_pick = [rng.integers(0, max(len(s),1), N) for s in data_srel]
gauss1, gauss2 = rng.normal(size=N), rng.normal(size=N)
gauss3 = rng.normal(size=N)

A_GRID  = np.array([0.0,0.5,0.75,1.0,1.25,1.5,2.0])
E_GRID  = np.array([0.8,1.05,1.3])
SM_GRID = np.array([0.20,0.25,0.30,0.35])
FT_GRID = np.array([0.0,0.05])
FC_GRID = np.array([0.0,0.02,0.10])

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

def model_p(o):
    """per-grid-point probability vectors: (SM,FT,FC,4,40)"""
    s3 = o[:,0,None]*ef+o[:,1,None]*e2
    v3 = o[:,2,None]*ef+o[:,3,None]*e2
    sp = np.linalg.norm(sky(s3),axis=1)
    vt0 = np.linalg.norm(sky(v3),axis=1)
    s_kau = sp/1e3
    P = np.full((len(SM_GRID),len(FT_GRID),len(FC_GRID),4,40), 1/40.0)
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500 or len(data_srel[bi]) == 0: continue
        vc = 2*np.pi*np.sqrt(M_s[idx]/sp[idx])
        sr = data_srel[bi][noise_pick[bi][idx] % len(data_srel[bi])]
        for si, sm in enumerate(SM_GRID):
            smear = np.exp(sm*gauss3[idx])
            for ti, ft in enumerate(FT_GRID):
                trip = is3_u[idx] < ft
                vv = vt0[idx]*np.where(trip, np.sqrt(1+m_hid[idx]/M_s[idx]), 1.0)
                wob = np.where(trip, v_wob[idx], 0.0)
                vtn = np.hypot((vv/vc)*smear + gauss1[idx]*sr
                               + wob*np.abs(wdir[idx,0])/vc,
                               gauss2[idx]*sr + wob*np.abs(wdir[idx,1])/vc)
                h,_ = np.histogram(np.clip(vtn, 0.021, 5.9), bins=VEDGE)
                p_orb = np.maximum(h/max(h.sum(),1), 1e-5); p_orb /= p_orb.sum()
                for ci, fc in enumerate(FC_GRID):
                    w = min(fc*SC2[bi], 0.5)
                    P[si, ti, ci, bi] = (1-w)*p_orb + w/40.0
    return P

t0 = time.time()
logp = {}
for law, TAB in (("simple", TAB_S), ("BE", TAB_B)):
    lp = np.full((len(A_GRID), len(E_GRID), len(SM_GRID),
                  len(FT_GRID), len(FC_GRID), 4, 40), np.nan)
    for ai, al in enumerate(A_GRID):
        tab_a = 1.0 + al*(TAB-1.0)
        for ei, eta in enumerate(E_GRID):
            if al == 0.0 and law == "BE":
                lp[ai, ei] = logp["simple"][ai, ei]; continue
            e_s = e_of(eta)
            vp = vp_c(e_s, tab_a) if al > 0 else None
            mode = 5 if al > 0 else 1
            kw = dict(a0=A0_CAN, tab=tab_a, lny0=LNY0, dlny=DLNY, vp=vp) \
                 if al > 0 else {}
            o = run(a_s, e_s, psi0, f_ip, M_s, uph, 8, 2500, mode, **kw)
            lp[ai, ei] = np.log(model_p(o))
    logp[law] = lp
    print(f"{law} model grid done, {(time.time()-t0)/60:.1f} min")

prior = -0.5*((E_GRID-1.3)/0.3)**2

def alpha_hat(prof):
    imax = int(np.nanargmax(prof))
    ahat = A_GRID[imax]
    if 0 < imax < len(A_GRID)-1:
        x = A_GRID[imax-1:imax+2]; y = prof[imax-1:imax+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0: ahat = -c1/(2*c2)
    return ahat, (0 < imax < len(A_GRID)-1)

NBOOT = 1000
brng = np.random.default_rng(7)
npairs = len(s_d)
res = {law: {'ahat': [], 'dnewt': [], 'interior': 0} for law in logp}
for k in range(NBOOT):
    pick = brng.integers(0, npairs, npairs)
    sb, vb = s_d[pick], vt_d[pick]
    counts = np.zeros((4, 40))
    for bi, b in enumerate(SBINS):
        m = (sb>=b[0])&(sb<b[1])
        counts[bi],_ = np.histogram(np.clip(vb[m], 0.021, 5.9), bins=VEDGE)
    for law in logp:
        lnl = np.tensordot(logp[law], counts, axes=([5,6],[0,1]))
        lnl += prior[None,:,None,None,None]
        prof = np.nanmax(lnl, axis=(1,2,3,4))
        ah, interior = alpha_hat(prof)
        res[law]['ahat'].append(ah)
        res[law]['dnewt'].append(np.nanmax(prof)-prof[0])
        res[law]['interior'] += int(interior)

L = [f"STAGE 3I bootstrap: seed {SEED} model grid, {NBOOT} data replicates",
     f"grids: alpha={A_GRID.tolist()}, eta={E_GRID.tolist()}, "
     f"sm={SM_GRID.tolist()}, ft={FT_GRID.tolist()}, fc={FC_GRID.tolist()}"]
REAL_SCATTER = {'simple': 0.111, 'BE': 0.154}   # Stage 3H
for law in res:
    ah = np.array(res[law]['ahat']); dn = np.array(res[law]['dnewt'])
    tot = np.sqrt(ah.std(ddof=1)**2 + REAL_SCATTER[law]**2)
    L += [f"{law}: bootstrap alpha_hat = {ah.mean():.3f} +/- {ah.std(ddof=1):.3f} "
          f"(16-84%: [{np.percentile(ah,16):.2f}, {np.percentile(ah,84):.2f}]); "
          f"interior {res[law]['interior']}/{NBOOT}",
          f"{law}: bootstrap Newton dlnL = {dn.mean():.1f} +/- {dn.std(ddof=1):.1f} "
          f"(min {dn.min():.1f})",
          f"{law}: COMBINED alpha sigma (boot + realization {REAL_SCATTER[law]}) "
          f"= {tot:.3f}"]
for l in L: print(l)
with open('data/stage3i_summary.txt','w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage3i_summary.txt")
