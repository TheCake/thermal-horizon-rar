"""Control: simple-nu through the SAME QUMOND solver/table path as BE-true,
separating interpolating-function choice from AQUAL-vs-QUMOND formulation."""
import numpy as np

src = open('calcs/stage2b_population.py').read()
ns = {}
exec(src.split('# ---------- validation')[0], ns)
run = ns['run']
GM = 4*np.pi**2
A0_CAN = 7.99e-7; GEXT = 1.9*A0_CAN

tabf = np.load('data/efe_boost_simple.npy')
y_t, b_t = tabf[0][::-1], tabf[1][::-1]
lny = np.log(y_t)
lny_u = np.linspace(lny[0], lny[-1], 512)
tab = np.interp(lny_u, lny, b_t)
lny0, dlny = lny_u[0], lny_u[1]-lny_u[0]
def boost_tab(y, eN=None):
    return np.interp(np.log(np.clip(y,1e-12,None)), lny_u, tab, right=1.0)

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

o = run(a_s,e_s,psi0,f_ip,M_s,uph,10,3000,5,a0=A0_CAN,
        vp=vp_consistent(a_s,e_s,M_s,boost_tab), tab=tab, lny0=lny0, dlny=dlny)

xhat = np.zeros((N,3)); xhat[:,0]=1
ef = xhat-nrm*nrm[:,[0]]
ef /= np.maximum(np.linalg.norm(ef,axis=1,keepdims=True),1e-12)
e2 = np.cross(nrm,ef)
los = rng.normal(size=(N,3)); los /= np.linalg.norm(los,axis=1,keepdims=True)
def sky(v3): return v3-los*np.sum(v3*los,axis=1,keepdims=True)
BINS = [(0.2,2),(2,6),(6,20),(20,50)]
s3 = o[:,0,None]*ef+o[:,1,None]*e2
v3 = o[:,2,None]*ef+o[:,3,None]*e2
sp = np.linalg.norm(sky(s3),axis=1)
vt = np.linalg.norm(sky(v3),axis=1)/(2*np.pi*np.sqrt(M_s/sp))
s_kau = sp/1e3
med = [np.median(vt[(s_kau>=b[0])&(s_kau<b[1])]) for b in BINS]
p = np.array(med)/med[0]
DATA = np.array([1.000,1.045,1.073,1.202]); SIG = np.array([0.001,0.018,0.028,0.058])
print("simple-nu via QUMOND table:", np.round(p,3))
print(f"chi2 (3 wide bins): {np.sum(((DATA[1:]-p[1:])/SIG[1:])**2):.1f}")
# --- localization diagnostics ---
print("ABS medians per bin:", np.round(med,3))
wide = (s_kau>=6)&(s_kau<30)
r_now = np.linalg.norm(s3,axis=1)
yy = (4*np.pi**2*M_s/r_now**2)/A0_CAN
print("wide: N=", wide.sum(), " median y=", np.round(np.median(yy[wide]),2))
vk = 2*np.pi*np.sqrt(M_s/r_now)
vnow = np.linalg.norm(v3,axis=1)
print("wide median v/vK(r_now):", np.round(np.median((vnow/vk)[wide]),3))
ctrl = (s_kau>=0.2)&(s_kau<2)
print("ctrl median v/vK(r_now):", np.round(np.median((vnow/vk)[ctrl]),3))
