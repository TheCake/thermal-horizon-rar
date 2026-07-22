"""
STAGE 2E: hardened comparison.
Upgrades over 2D: (1) empirical eccentricity law alpha(s) from Hwang+22;
(2) SELF-CONSISTENT initial conditions in the modified potential (perihelion
speed set by energy+angular-momentum in the actual force law, fixing the
circularization bias); (3) distribution-SHAPE statistics per bin: a gravity
boost lifts all quantiles; triples fatten only the tail (P90/P50).
"""
import time
import numpy as np
from astropy.io import fits

# ---------------- model machinery ----------------
src = open('calcs/stage2b_population.py').read()
ns = {}
exec(src.split('# ---------- validation')[0], ns)
run, A0, GEXT = ns['run'], ns['A0'], ns['GEXT']
GM = 4*np.pi**2

def cm_boost(y, eN):
    be = 1.1*eN
    yb = np.sqrt(y*y + be*be)
    sq = np.sqrt(0.25 + 1.0/yb)
    nus = 0.5 + sq
    nuhat = (1.0/yb)/(2.0*nus*sq)
    return nus*(1.0 + np.tanh((be/np.maximum(y,1e-12))**1.2)*nuhat/3.0)

def vp_consistent(a, e, M):
    """perihelion speed in the C&M potential with apoapsis a(1+e)."""
    rp, ra = a*(1-e), a*(1+e)
    # Gauss-Legendre in log r for dPhi = int g dr
    xg, wg = np.polynomial.legendre.leggauss(48)
    lo, hi = np.log(rp), np.log(ra)
    lr = 0.5*(hi-lo)[:,None]*xg[None,:] + 0.5*(hi+lo)[:,None]
    r = np.exp(lr)
    gN = GM*M[:,None]/r**2
    g = cm_boost(gN/A0, GEXT/A0)*gN
    dPhi = np.sum(wg[None,:]*g*r, axis=1)*0.5*(hi-lo)   # int g dr = int g r dlnr
    v2 = 2*dPhi/(1-(rp/ra)**2)
    return np.sqrt(np.maximum(v2, 0))

rng = np.random.default_rng(31)
N = 1_000_000
u = rng.random(N); lo, hi, g = 0.15, 60.0, -0.6
a_s  = ((lo**g + u*(hi**g - lo**g))**(1/g))*1e3
# Hwang+22 alpha(s): uniform->superthermal with separation
alpha = np.interp(np.log10(a_s), np.log10([100,500,1000,50000]),
                  [0.6, 1.0, 1.2, 1.3])
e_s  = 0.95*rng.random(N)**(1/(1+alpha))
psi0 = rng.random(N)*2*np.pi
nrm  = rng.normal(size=(N,3)); nrm /= np.linalg.norm(nrm,axis=1,keepdims=True)
f_ip = np.sqrt(np.clip(1-nrm[:,0]**2, 0, 1))
M_s  = 0.6 + 1.8*rng.random(N)
uph  = rng.random(N)

t0 = time.time()
vpc = vp_consistent(a_s, e_s, M_s)
print(f"self-consistent ICs computed in {time.time()-t0:.0f} s")
pops = {}
t0 = time.time(); pops['newton'] = run(a_s,e_s,psi0,f_ip,M_s,uph,10,3000,1)
print(f"newton pop: {time.time()-t0:.0f} s")
t0 = time.time(); pops['aqual'] = run(a_s,e_s,psi0,f_ip,M_s,uph,10,3000,2,vp=vpc)
print(f"aqual pop (consistent ICs): {time.time()-t0:.0f} s")
# sanity: apoapsis fidelity of the aqual population
ra_target = a_s*(1+e_s)
print(f"aqual rmax/target-apoapsis median: "
      f"{np.median(pops['aqual'][:,6]/ra_target):.3f} (want ~1)")

# ---------------- shared projection + triples ----------------
xhat = np.zeros((N,3)); xhat[:,0] = 1
ef = xhat - nrm*nrm[:,[0]]
ef /= np.maximum(np.linalg.norm(ef,axis=1,keepdims=True),1e-12)
e2 = np.cross(nrm, ef)
los = rng.normal(size=(N,3)); los /= np.linalg.norm(los,axis=1,keepdims=True)
def sky(v3): return v3 - los*np.sum(v3*los,axis=1,keepdims=True)
is3_u = rng.random(N)
q_in  = 0.1 + 0.9*rng.random(N)
a_in  = 10**(np.log10(2)+rng.random(N)*(np.log10(100)-np.log10(2)))
wdir  = rng.normal(size=(N,3)); wdir /= np.linalg.norm(wdir,axis=1,keepdims=True)
v_wob = (q_in/(1+q_in))*29.78*np.sqrt(0.5*M_s*(1+q_in)/a_in)/4.74047
m_hid = q_in*0.5*M_s

BINS = [(0.2,2),(2,6),(6,20),(20,50)]
def profile(o, f_t):
    s3 = o[:,0,None]*ef + o[:,1,None]*e2
    v3 = o[:,2,None]*ef + o[:,3,None]*e2
    trip = is3_u < f_t
    v3 = v3*np.where(trip, np.sqrt(1+m_hid/M_s), 1.0)[:,None]
    vp = sky(v3) + np.where(trip,1.0,0.0)[:,None]*sky(wdir*v_wob[:,None])
    sp = np.linalg.norm(sky(s3), axis=1)
    vt = np.linalg.norm(vp, axis=1)/(2*np.pi*np.sqrt(M_s/sp))
    s_kau = sp/1e3
    out = []
    for b in BINS:
        m = (s_kau>=b[0])&(s_kau<b[1])
        v = vt[m]
        out.append((np.median(v), np.percentile(v,90)/np.median(v)))
    return out

# ---------------- data profile ----------------
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
data_prof = []
for b in BINS:
    m = (s_kau>=b[0])&(s_kau<b[1])
    v = vt_d[m]
    data_prof.append((np.median(v), np.percentile(v,90)/np.median(v), m.sum()))

# ---------------- verdict tables ----------------
print(f"\nPER-BIN MEDIAN v-tilde, normalized to the 0.2-2 kAU bin")
hdr = f"{'bin [kAU]':>10} {'DATA':>7}"
for f_t in (0.0, 0.1, 0.2): hdr += f" {'Nw f='+str(f_t):>9}"
for f_t in (0.0, 0.1, 0.2): hdr += f" {'AQ f='+str(f_t):>9}"
print(hdr)
profs = {(lab,f): profile(pops[lab], f)
         for lab in ('newton','aqual') for f in (0.0,0.1,0.2)}
for i,b in enumerate(BINS):
    row = f"{str(b):>10} {data_prof[i][0]/data_prof[0][0]:>7.3f}"
    for lab in ('newton','aqual'):
        for f in (0.0,0.1,0.2):
            p = profs[(lab,f)]
            row += f" {p[i][0]/p[0][0]:>9.3f}"
    print(row)

print(f"\nTAIL SHAPE P90/P50 per bin (triples fatten tails; boosts do not)")
print(f"{'bin [kAU]':>10} {'DATA':>7} {'Nw f=0':>8} {'Nw f=0.2':>9} "
      f"{'AQ f=0':>8} {'AQ f=0.2':>9}")
for i,b in enumerate(BINS):
    print(f"{str(b):>10} {data_prof[i][1]:>7.2f} "
          f"{profs[('newton',0.0)][i][1]:>8.2f} {profs[('newton',0.2)][i][1]:>9.2f} "
          f"{profs[('aqual',0.0)][i][1]:>8.2f} {profs[('aqual',0.2)][i][1]:>9.2f}")
print(f"\ndata bin counts: {[p[2] for p in data_prof]}")
