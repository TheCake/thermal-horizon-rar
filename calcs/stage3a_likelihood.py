"""
STAGE 3A: distribution-level likelihood engine.
Per separation bin, evaluate the log-likelihood of the DATA's full v-tilde
distribution under each model population (noise-convolved with the data's own
per-bin error distribution). Each modified model runs under TWO IC conventions:
  A: (rp,ra)-preserving self-consistent ICs   B: Keplerian ICs (orbit relaxes)
The convention spread = IC systematic band. Ranking is trusted only where
Delta(lnL) between models exceeds the band.
"""
import time
import numpy as np
from astropy.io import fits

src = open('calcs/stage2b_population.py').read()
ns = {}
exec(src.split('# ---------- validation')[0], ns)
run = ns['run']
GM = 4*np.pi**2; A0_CAN = 7.99e-7; GEXT = 1.9*A0_CAN

tabf = np.load('data/efe_boost_be.npy')
y_t, b_t = tabf[0][::-1], tabf[1][::-1]
lny = np.log(y_t); lny_u = np.linspace(lny[0], lny[-1], 512)
tab_be = np.interp(lny_u, lny, b_t)
lny0, dlny = lny_u[0], lny_u[1]-lny_u[0]
def boost_be(y):
    return np.interp(np.log(np.clip(y,1e-12,None)), lny_u, tab_be, right=1.0)
def boost_cm(y):
    eN = GEXT/A0_CAN; be = 1.1*eN; yb = np.sqrt(y*y+be*be)
    sq = np.sqrt(0.25+1.0/yb); nus = 0.5+sq
    nuhat = (1.0/yb)/(2.0*nus*sq)
    return nus*(1.0+np.tanh((be/np.maximum(y,1e-12))**1.2)*nuhat/3.0)

def vp_consistent(a, e, M, boost):
    rp, ra = a*(1-e), a*(1+e)
    xg, wg = np.polynomial.legendre.leggauss(48)
    lo_, hi_ = np.log(rp), np.log(ra)
    lr = 0.5*(hi_-lo_)[:,None]*xg[None,:] + 0.5*(hi_+lo_)[:,None]
    r = np.exp(lr)
    gN = GM*M[:,None]/r**2
    gg = boost(gN/A0_CAN)*gN
    dPhi = np.sum(wg[None,:]*gg*r, axis=1)*0.5*(hi_-lo_)
    return np.sqrt(np.maximum(2*dPhi/(1-(rp/ra)**2), 0))

rng = np.random.default_rng(31)
N = 1_000_000
u = rng.random(N); lo, hi, g = 0.15, 60.0, -0.6
a_s  = ((lo**g+u*(hi**g-lo**g))**(1/g))*1e3
alpha = np.interp(np.log10(a_s), np.log10([100,500,1000,50000]),
                  [0.6,1.0,1.2,1.3])
e_s  = 0.95*rng.random(N)**(1/(1+alpha))
psi0 = rng.random(N)*2*np.pi
nrm  = rng.normal(size=(N,3)); nrm /= np.linalg.norm(nrm,axis=1,keepdims=True)
f_ip = np.sqrt(np.clip(1-nrm[:,0]**2,0,1))
M_s  = 0.6+1.8*rng.random(N)
uph  = rng.random(N)

runs = {
 'newton':    dict(mode=1, kw={}),
 'aqual-A':   dict(mode=2, kw=dict(a0=A0_CAN, vp=vp_consistent(a_s,e_s,M_s,boost_cm))),
 'aqual-B':   dict(mode=2, kw=dict(a0=A0_CAN)),
 'BE-A':      dict(mode=5, kw=dict(a0=A0_CAN, vp=vp_consistent(a_s,e_s,M_s,boost_be),
                                   tab=tab_be, lny0=lny0, dlny=dlny)),
 'BE-B':      dict(mode=5, kw=dict(a0=A0_CAN, tab=tab_be, lny0=lny0, dlny=dlny)),
}
pops = {}
for name, cfg in runs.items():
    t0 = time.time()
    pops[name] = run(a_s,e_s,psi0,f_ip,M_s,uph,10,3000,cfg['mode'],**cfg['kw'])
    print(f"{name:>9}: {time.time()-t0:.0f} s")

# geometry/projection (shared)
xhat = np.zeros((N,3)); xhat[:,0]=1
ef = xhat-nrm*nrm[:,[0]]
ef /= np.maximum(np.linalg.norm(ef,axis=1,keepdims=True),1e-12)
e2 = np.cross(nrm,ef)
los = rng.normal(size=(N,3)); los /= np.linalg.norm(los,axis=1,keepdims=True)
def sky(v3): return v3-los*np.sum(v3*los,axis=1,keepdims=True)

# ---- data: v-tilde and per-pair relative noise, per bin ----
d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1,plx2 = d['parallax1'],d['parallax2']
eplx1,eplx2 = d['parallax_error1'],d['parallax_error2']
sep,Rch = d['sep_AU'],d['R_chance_align']
G1,G2 = d['phot_g_mean_mag1'],d['phot_g_mean_mag2']
plx = 0.5*(plx1+plx2)
MG1 = G1+5*np.log10(np.maximum(plx1,1e-6))-10
MG2 = G2+5*np.log10(np.maximum(plx2,1e-6))-10
sig = 4.74047/plx*np.sqrt(d['pmra_error1']**2+d['pmdec_error1']**2
                          +d['pmra_error2']**2+d['pmdec_error2']**2)
ok = (Rch<0.01)&(plx1>5)&(plx2>5)&(plx1/np.maximum(eplx1,1e-6)>20) \
   &(plx2/np.maximum(eplx2,1e-6)>20)&(np.abs(plx1-plx2)<3*np.hypot(eplx1,eplx2)) \
   &(sep>200)&(sep<50000)&(MG1>2.6)&(MG1<14.2)&(MG2>2.6)&(MG2<14.2)&(sig<0.03)
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot = np.interp(MG1[ok],MG_T,MS_T)+np.interp(MG2[ok],MG_T,MS_T)
s_d = sep[ok]/1e3
vc_d = 0.9417*np.sqrt(Mtot/s_d)
dv = 4.74047/plx[ok]*np.hypot(d['pmra1'][ok]-d['pmra2'][ok],
                              d['pmdec1'][ok]-d['pmdec2'][ok])
vt_d = dv/vc_d
srel_d = sig[ok]/vc_d

BINS = [(0.2,2),(2,6),(6,20),(20,50)]
EDGES = np.logspace(np.log10(0.02), np.log10(6.0), 41)
def lnL(o):
    s3 = o[:,0,None]*ef+o[:,1,None]*e2
    v3 = o[:,2,None]*ef+o[:,3,None]*e2
    sp = np.linalg.norm(sky(s3),axis=1)
    vt = np.linalg.norm(sky(v3),axis=1)/(2*np.pi*np.sqrt(M_s/sp))
    s_kau = sp/1e3
    tot = 0.0; per = []
    for b in BINS:
        mm = (s_kau>=b[0])&(s_kau<b[1])
        md = (s_d>=b[0])&(s_d<b[1])
        v = vt[mm]
        # convolve model with the data's own noise distribution in this bin
        sr = rng.choice(srel_d[md], size=len(v), replace=True)
        v_noisy = np.hypot(v + rng.normal(size=len(v))*sr,
                           rng.normal(size=len(v))*sr)
        h, _ = np.histogram(v_noisy, bins=EDGES, density=True)
        h = np.maximum(h, 1e-4)
        idx = np.clip(np.searchsorted(EDGES, vt_d[md])-1, 0, len(h)-1)
        ll = np.sum(np.log(h[idx]))
        per.append(ll); tot += ll
    return tot, per

print(f"\n{'model':>9} {'lnL total':>11} " + " ".join(f"{str(b):>12}" for b in BINS))
base = None
for name in runs:
    tot, per = lnL(pops[name])
    if base is None: base = tot
    print(f"{name:>9} {tot-base:>+11.1f} " + " ".join(f"{p:>12.1f}" for p in per))
print("\n(lnL relative to 'newton'; higher = better. Compare A-vs-B spread per")
print(" model against the A/B-to-model differences to judge ranking robustness.)")
