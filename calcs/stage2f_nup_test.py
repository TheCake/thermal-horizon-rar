"""
STAGE 2F: the consilience test.
Question: does the SPARC+Cassini-measured screening family (nu_p, p=0.443,
a0=1.03e-10), inserted into the Chae-Milgrom EFE structure, PREDICT the
6-20 kAU dip that the standard simple-nu AQUAL leaves unexplained?
Three model populations (self-consistent ICs each), one data profile with
bootstrap errors, one table.
"""
import time
import numpy as np
from astropy.io import fits

src = open('calcs/stage2b_population.py').read()
ns = {}
exec(src.split('# ---------- validation')[0], ns)
run = ns['run']
GM = 4*np.pi**2
A0_OURS = 1.03e-10*6.658e3
A0_CAN  = 7.99e-7
GEXT    = 1.9*A0_CAN
P_OURS  = 0.443

def boost_cm(y, eN):
    be = 1.1*eN; yb = np.sqrt(y*y + be*be)
    sq = np.sqrt(0.25 + 1.0/yb); nus = 0.5 + sq
    nuhat = (1.0/yb)/(2.0*nus*sq)
    return nus*(1.0 + np.tanh((be/np.maximum(y,1e-12))**1.2)*nuhat/3.0)

def boost_np(y, eN, p=P_OURS):
    be = 1.1*eN; yb = np.sqrt(y*y + be*be)
    t = yb**p; et = np.exp(-t)
    nus = (1.0-et)**(-0.5/p)
    nuhat = t*et/(2.0*(1.0-et))
    return nus*(1.0 + np.tanh((be/np.maximum(y,1e-12))**1.2)*nuhat/3.0)

def vp_consistent(a, e, M, boost, a0):
    rp, ra = a*(1-e), a*(1+e)
    xg, wg = np.polynomial.legendre.leggauss(48)
    lo, hi = np.log(rp), np.log(ra)
    lr = 0.5*(hi-lo)[:,None]*xg[None,:] + 0.5*(hi+lo)[:,None]
    r = np.exp(lr)
    gN = GM*M[:,None]/r**2
    g = boost(gN/a0, GEXT/a0)*gN
    dPhi = np.sum(wg[None,:]*g*r, axis=1)*0.5*(hi-lo)
    return np.sqrt(np.maximum(2*dPhi/(1-(rp/ra)**2), 0))

rng = np.random.default_rng(31)
N = 1_000_000
u = rng.random(N); lo, hi, g = 0.15, 60.0, -0.6
a_s  = ((lo**g + u*(hi**g - lo**g))**(1/g))*1e3
alpha = np.interp(np.log10(a_s), np.log10([100,500,1000,50000]),
                  [0.6,1.0,1.2,1.3])
e_s  = 0.95*rng.random(N)**(1/(1+alpha))
psi0 = rng.random(N)*2*np.pi
nrm  = rng.normal(size=(N,3)); nrm /= np.linalg.norm(nrm,axis=1,keepdims=True)
f_ip = np.sqrt(np.clip(1-nrm[:,0]**2,0,1))
M_s  = 0.6 + 1.8*rng.random(N)
uph  = rng.random(N)

pops = {}
t0=time.time(); pops['newton'] = run(a_s,e_s,psi0,f_ip,M_s,uph,10,3000,1)
print(f"newton: {time.time()-t0:.0f} s")
vpc = vp_consistent(a_s,e_s,M_s,boost_cm,A0_CAN)
t0=time.time()
pops['aqual-CM'] = run(a_s,e_s,psi0,f_ip,M_s,uph,10,3000,2,a0=A0_CAN,vp=vpc)
print(f"aqual-CM (canonical a0, simple nu): {time.time()-t0:.0f} s")
vpp = vp_consistent(a_s,e_s,M_s,boost_np,A0_OURS)
t0=time.time()
pops['nu_p-EFE'] = run(a_s,e_s,psi0,f_ip,M_s,uph,10,3000,3,a0=A0_OURS,
                       pscr=P_OURS,vp=vpp)
print(f"nu_p-EFE (our measured law): {time.time()-t0:.0f} s")

# projection
xhat = np.zeros((N,3)); xhat[:,0]=1
ef = xhat - nrm*nrm[:,[0]]
ef /= np.maximum(np.linalg.norm(ef,axis=1,keepdims=True),1e-12)
e2 = np.cross(nrm, ef)
los = rng.normal(size=(N,3)); los /= np.linalg.norm(los,axis=1,keepdims=True)
def sky(v3): return v3 - los*np.sum(v3*los,axis=1,keepdims=True)
BINS = [(0.2,2),(2,6),(6,20),(20,50)]
def prof(o):
    s3 = o[:,0,None]*ef + o[:,1,None]*e2
    v3 = o[:,2,None]*ef + o[:,3,None]*e2
    sp = np.linalg.norm(sky(s3),axis=1)
    vt = np.linalg.norm(sky(v3),axis=1)/(2*np.pi*np.sqrt(M_s/sp))
    s_kau = sp/1e3
    med = [np.median(vt[(s_kau>=b[0])&(s_kau<b[1])]) for b in BINS]
    return np.array(med)/med[0]

# data profile with bootstrap
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
s_kau = sep[ok]/1e3
dv = 4.74047/plx[ok]*np.hypot(d['pmra1'][ok]-d['pmra2'][ok],
                              d['pmdec1'][ok]-d['pmdec2'][ok])
vt_d = dv/(0.9417*np.sqrt(Mtot/s_kau))
idx = [np.where((s_kau>=b[0])&(s_kau<b[1]))[0] for b in BINS]
data_med = np.array([np.median(vt_d[i]) for i in idx])
data_prof = data_med/data_med[0]
boots = []
for _ in range(500):
    m = [np.median(vt_d[rng.choice(i,len(i))]) for i in idx]
    boots.append(np.array(m)/m[0])
boots = np.array(boots)
lo_ci, hi_ci = np.percentile(boots,16,axis=0), np.percentile(boots,84,axis=0)

profs = {k: prof(v) for k,v in pops.items()}
print(f"\n{'bin [kAU]':>10} {'DATA (68% CI)':>22} {'Newton':>8} "
      f"{'AQUAL-CM':>9} {'nu_p-EFE':>9}")
for i,b in enumerate(BINS):
    print(f"{str(b):>10} {data_prof[i]:>8.3f} ({lo_ci[i]:.3f}-{hi_ci[i]:.3f})"
          f" {profs['newton'][i]:>8.3f} {profs['aqual-CM'][i]:>9.3f}"
          f" {profs['nu_p-EFE'][i]:>9.3f}")
# chi-square-ish per model over the three wide bins
sig_b = 0.5*(hi_ci-lo_ci)
for k in profs:
    chi2 = np.sum(((data_prof[1:]-profs[k][1:])/sig_b[1:])**2)
    print(f"chi2 (3 wide bins) {k:>9}: {chi2:6.1f}")
