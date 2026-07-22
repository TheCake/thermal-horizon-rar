"""
STAGE 3E: mass-error smearing test. The 3D autopsy shows the data vtilde
distribution is BROADER than every model (low-vtilde excess + 0.8-1.5 excess,
0.3-0.8 deficit, in all s-bins incl. the Newtonian control bin). Data vtilde is
normalized by photometric masses (error ~delta M/M), so it is multiplicatively
smeared by ~0.5*dM/M relative to truth; the model is not. Test: smear model
vtilde by exp(sigma_m * g) over a sigma_m grid. Questions: (1) does sigma_m>0
improve ALL models (missing ingredient)? (2) does the alpha=1 -> alpha=2 lnL
gain collapse (corner-seeking was misfit absorption)? (3) control-bin chi2?
Same realization as v3 (seed 31), eta=1.3, f_t=0.05, f_c0=0.1.
Writes data/stage3e_smear.txt.
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
gauss3 = rng.normal(size=N)          # NEW draw: mass-error smear (after v3 draws)

ETA = 1.3; FT = 0.05; FC = 0.10
al = np.interp(np.log10(a_s), np.log10([100,500,1000,50000]),
               [0.6, 1.0, ETA, ETA])
e_s = 0.95*u_e**(1/(1+al))

def vp_c(tab_a):
    rp, ra = a_s*(1-e_s), a_s*(1+e_s)
    xg, wg = np.polynomial.legendre.leggauss(32)
    lo_, hi_ = np.log(rp), np.log(ra)
    lr = 0.5*(hi_-lo_)[:,None]*xg[None,:] + 0.5*(hi_+lo_)[:,None]
    r = np.exp(lr)
    gN = GM*M_s[:,None]/r**2
    bst = np.interp(np.log(gN/A0_CAN), LNY_U, tab_a, right=1.0)
    dPhi = np.sum(wg[None,:]*bst*gN*r, axis=1)*0.5*(hi_-lo_)
    return np.sqrt(np.maximum(2*dPhi/(1-(rp/ra)**2), 0))

SM_GRID = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]

def eval_bins(o, sm):
    s3 = o[:,0,None]*ef+o[:,1,None]*e2
    v3 = o[:,2,None]*ef+o[:,3,None]*e2
    sp = np.linalg.norm(sky(s3),axis=1)
    vt0 = np.linalg.norm(sky(v3),axis=1)
    s_kau = sp/1e3
    lnl = np.zeros(4); chi = np.zeros(4)
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500: continue
        vc = 2*np.pi*np.sqrt(M_s[idx]/sp[idx])
        sr = data_srel[bi][noise_pick[bi][idx] % len(data_srel[bi])]
        trip = is3_u[idx] < FT
        vv = vt0[idx]*np.where(trip, np.sqrt(1+m_hid[idx]/M_s[idx]), 1.0)
        smear = np.exp(sm*gauss3[idx])
        wob = np.where(trip, v_wob[idx], 0.0)
        vtn = np.hypot((vv/vc)*smear + gauss1[idx]*sr + wob*np.abs(wdir[idx,0])/vc,
                       gauss2[idx]*sr + wob*np.abs(wdir[idx,1])/vc)
        h,_ = np.histogram(np.clip(vtn, 0.021, 5.9), bins=VEDGE)
        p_orb = np.maximum(h/max(h.sum(),1), 1e-5); p_orb /= p_orb.sum()
        w = min(FC*SC2[bi], 0.5)
        p = (1-w)*p_orb + w/40.0
        lnl[bi] = np.sum(data_counts[bi]*np.log(p))
        n = data_counts[bi]; Nb = n.sum()
        chi[bi] = np.sum((n-Nb*p)**2/(Nb*p))/39
    return lnl, chi

MODELS = [
    ("Newton",   None,  0.0),
    ("simple a1", TAB_S, 1.0),
    ("BE a1",     TAB_B, 1.0),
    ("simple a2", TAB_S, 2.0),
    ("BE a2",     TAB_B, 2.0),
]
L = []
def P(s):
    print(s); L.append(s)

P(f"STAGE 3E: model-vtilde smear test, sigma_m grid {SM_GRID}")
P(f"(sigma_m = 0.5 * fractional total-mass error; eta={ETA}, f_t={FT}, f_c0={FC})")
P("")
tot = {}; ctl_chi = {}
for name, TAB, alpha in MODELS:
    if alpha == 0:
        o = run(a_s, e_s, psi0, f_ip, M_s, uph, 8, 2500, 1)
    else:
        tab_a = 1.0 + alpha*(TAB-1.0)
        o = run(a_s, e_s, psi0, f_ip, M_s, uph, 8, 2500, 5,
                a0=A0_CAN, tab=tab_a, lny0=LNY0, dlny=DLNY, vp=vp_c(tab_a))
    row_l, row_c = [], []
    for sm in SM_GRID:
        lnl, chi = eval_bins(o, sm)
        row_l.append(lnl.sum()); row_c.append(chi[0])
    tot[name] = np.array(row_l); ctl_chi[name] = np.array(row_c)
    print(f"{name} done")

P("Total lnL vs sigma_m (rel. to each model's sigma_m=0):")
P(f"{'model':<10}" + "".join(f"{s:>9}" for s in SM_GRID))
for name in tot:
    r = tot[name] - tot[name][0]
    P(f"{name:<10}" + "".join(f"{v:>+9.1f}" for v in r))
P("")
P("Absolute total lnL rel. to Newton at same sigma_m:")
for name in tot:
    if name == 'Newton': continue
    r = tot[name] - tot['Newton']
    P(f"{name:<10}" + "".join(f"{v:>+9.1f}" for v in r))
P("")
P("Control-bin (0.2-2 kAU) chi2/dof vs sigma_m:")
for name in ctl_chi:
    P(f"{name:<10}" + "".join(f"{v:>9.2f}" for v in ctl_chi[name]))
P("")
best_sm = {n: SM_GRID[int(np.argmax(tot[n]))] for n in tot}
P(f"best sigma_m per model: {best_sm}")
a2a1_0 = (tot['simple a2']-tot['simple a1'])[0]
i_s = int(np.argmax(tot['simple a1']))
a2a1_b = (tot['simple a2']-tot['simple a1'])[i_s]
P(f"simple: (a2 - a1) lnL at sigma_m=0: {a2a1_0:+.1f}; at best-a1 sigma_m: {a2a1_b:+.1f}")
b0 = (tot['BE a2']-tot['BE a1'])[0]
i_b = int(np.argmax(tot['BE a1']))
bb = (tot['BE a2']-tot['BE a1'])[i_b]
P(f"BE:     (a2 - a1) lnL at sigma_m=0: {b0:+.1f}; at best-a1 sigma_m: {bb:+.1f}")

with open('data/stage3e_smear.txt','w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage3e_smear.txt")
