"""
STAGE 2D-b: triple-star forward model. GPU populations under Newton and
Chae-Milgrom AQUAL; inject RUWE-invisible triples (inner companions at
2-100 AU); grid the triple fraction; compare boost ratios with the data.
DATA targets: boost = 1.086 (baseline) / 1.069 (strict RUWE).
"""
import time
import numpy as np
import importlib.util

spec = importlib.util.spec_from_file_location(
    's2b', 'calcs/stage2b_population.py')
# import only the kernel machinery: execute the file up to the validation block
src = open('calcs/stage2b_population.py').read()
head = src.split('# ---------- validation')[0]
ns = {}
exec(head, ns)
run = ns['run']

rng = np.random.default_rng(21)
N = 1_000_000
u = rng.random(N); lo, hi, g = 0.2, 50.0, -0.6
a_s  = ((lo**g + u*(hi**g - lo**g))**(1/g))*1e3
e_s  = 0.95*rng.random(N)**(1/2.25)
psi0 = rng.random(N)*2*np.pi
nrm  = rng.normal(size=(N,3)); nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
f_ip = np.sqrt(np.clip(1-nrm[:,0]**2, 0, 1))
M_s  = 0.6 + 1.8*rng.random(N)
uph  = rng.random(N)

pops = {}
for label, mode in (("newton", 1), ("aqual", 2)):
    t0 = time.time()
    pops[label] = run(a_s, e_s, psi0, f_ip, M_s, uph, 10, 3000, mode)
    print(f"{label} population: {time.time()-t0:.0f} s")

# --- 3D reconstruction and sky projection (shared geometry) ---
xhat = np.zeros((N,3)); xhat[:,0] = 1
ef = xhat - nrm*nrm[:,[0]]
ef /= np.maximum(np.linalg.norm(ef, axis=1, keepdims=True), 1e-12)
e2 = np.cross(nrm, ef)
los = rng.normal(size=(N,3)); los /= np.linalg.norm(los, axis=1, keepdims=True)
def sky(v3):
    return v3 - los*np.sum(v3*los, axis=1, keepdims=True)

# --- triple wobble ingredients (same draw reused across scenarios) ---
is3_u   = rng.random(N)
q_in    = 0.1 + 0.9*rng.random(N)
a_in    = 10**(np.log10(2) + rng.random(N)*(np.log10(100)-np.log10(2)))
wdir    = rng.normal(size=(N,3)); wdir /= np.linalg.norm(wdir, axis=1, keepdims=True)
m_host  = 0.5*M_s
v_wob   = (q_in/(1+q_in))*29.78*np.sqrt(m_host*(1+q_in)/a_in)/4.74047  # AU/yr
m_hidden = q_in*m_host

def boost_ratio(o, f_t):
    s3 = o[:,0,None]*ef + o[:,1,None]*e2
    v3 = o[:,2,None]*ef + o[:,3,None]*e2
    trip = is3_u < f_t
    # hidden-mass speedup of the outer orbit
    v3 = v3*np.where(trip, np.sqrt(1+m_hidden/M_s), 1.0)[:,None]
    vp = sky(v3) + np.where(trip, 1.0, 0.0)[:,None]*sky(wdir*v_wob[:,None])
    sp = np.linalg.norm(sky(s3), axis=1)
    vt = np.linalg.norm(vp, axis=1)/(2*np.pi*np.sqrt(M_s/sp))
    s_kau = sp/1e3
    mn = (s_kau >= 0.2) & (s_kau < 2); mw = (s_kau >= 6) & (s_kau < 30)
    return np.median(vt[mw])/np.median(vt[mn])

print(f"\n{'f_triple':>9} {'Newton boost':>13} {'AQUAL boost':>12}")
for f_t in (0.0, 0.05, 0.10, 0.15, 0.20, 0.30):
    bn = boost_ratio(pops['newton'], f_t)
    ba = boost_ratio(pops['aqual'], f_t)
    print(f"{f_t:>9.2f} {bn:>13.3f} {ba:>12.3f}")
print("\nDATA: 1.086 (68% 1.065-1.105); strict-RUWE 1.069 (1.036-1.110)")
