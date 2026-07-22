"""
STAGE 3D: control-bin autopsy (TODO #2). Where in (s-bin, vtilde) space is the
likelihood discrimination generated? Five focused models on the SAME realization
as stage3c (seed 31, eta=1.3): Newton, simple/BE at alpha=1, simple/BE at
alpha=2. Per-bin lnL, per-vtilde-region decomposition of the differences, and
model-vs-data shape fractions in each bin. Writes data/stage3d_autopsy.txt.
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
VCEN = np.sqrt(VEDGE[:-1]*VEDGE[1:])
data_counts, data_srel = [], []
for b in SBINS:
    m = (s_d>=b[0])&(s_d<b[1])
    c,_ = np.histogram(np.clip(vt_d[m], 0.021, 5.9), bins=VEDGE)
    data_counts.append(c.astype(float)); data_srel.append(srel_d[m])

# ---- identical realization to stage3c (seed 31) ----
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

def eval_bins(o, fc):
    """returns (lnl_per_bin[4], p_model[4][40])"""
    s3 = o[:,0,None]*ef+o[:,1,None]*e2
    v3 = o[:,2,None]*ef+o[:,3,None]*e2
    sp = np.linalg.norm(sky(s3),axis=1)
    vt0 = np.linalg.norm(sky(v3),axis=1)
    s_kau = sp/1e3
    lnl = np.zeros(4); pm = [None]*4
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500: continue
        vc = 2*np.pi*np.sqrt(M_s[idx]/sp[idx])
        sr = data_srel[bi][noise_pick[bi][idx] % len(data_srel[bi])]
        trip = is3_u[idx] < FT
        vv = vt0[idx]*np.where(trip, np.sqrt(1+m_hid[idx]/M_s[idx]), 1.0)
        wob = np.where(trip, v_wob[idx], 0.0)
        vtn = np.hypot(vv/vc + gauss1[idx]*sr + wob*np.abs(wdir[idx,0])/vc,
                       gauss2[idx]*sr + wob*np.abs(wdir[idx,1])/vc)
        h,_ = np.histogram(np.clip(vtn, 0.021, 5.9), bins=VEDGE)
        p_orb = np.maximum(h/max(h.sum(),1), 1e-5); p_orb /= p_orb.sum()
        w = min(fc*SC2[bi], 0.5)
        p = (1-w)*p_orb + w/40.0
        lnl[bi] = np.sum(data_counts[bi]*np.log(p))
        pm[bi] = p
    return lnl, pm

MODELS = [
    ("Newton",   None,  0.0),
    ("simple a1", TAB_S, 1.0),
    ("BE a1",     TAB_B, 1.0),
    ("simple a2", TAB_S, 2.0),
    ("BE a2",     TAB_B, 2.0),
]
res = {}
for name, TAB, alpha in MODELS:
    t0 = time.time()
    if alpha == 0:
        o = run(a_s, e_s, psi0, f_ip, M_s, uph, 8, 2500, 1)
    else:
        tab_a = 1.0 + alpha*(TAB-1.0)
        o = run(a_s, e_s, psi0, f_ip, M_s, uph, 8, 2500, 5,
                a0=A0_CAN, tab=tab_a, lny0=LNY0, dlny=DLNY, vp=vp_c(tab_a))
    lnl, pm = eval_bins(o, FC)
    lnl0, pm0 = eval_bins(o, 0.0)
    res[name] = (lnl, pm, lnl0, pm0)
    print(f"{name}: {time.time()-t0:.0f}s  lnL/bin = {np.round(lnl,1).tolist()}")

L = []
def P(s):
    print(s); L.append(s)

P(f"STAGE 3D autopsy — eta={ETA}, f_t={FT}, f_c0={FC}, seed 31 (matches v3)")
P(f"data counts/bin: {[int(c.sum()) for c in data_counts]}   bins: {SBINS}")
P("")
P("Per-bin lnL, difference vs Newton (positive = better than Newton):")
hdr = f"{'model':<10}" + "".join(f"{f'{b[0]}-{b[1]}kAU':>14}" for b in SBINS) + f"{'total':>10}"
P(hdr)
ln_new = res['Newton'][0]
for name in res:
    lnl = res[name][0]
    dif = lnl - ln_new
    if name == 'Newton':
        P(f"{'Newton':<10}" + "".join(f"{v:>14.1f}" for v in lnl) + f"{lnl.sum():>10.1f}")
    else:
        P(f"{name:<10}" + "".join(f"{v:>+14.1f}" for v in dif) + f"{dif.sum():>+10.1f}")
P("")
P("Same, contamination OFF (f_c0=0):")
ln_new0 = res['Newton'][2]
for name in res:
    if name == 'Newton': continue
    dif = res[name][2] - ln_new0
    P(f"{name:<10}" + "".join(f"{v:>+14.1f}" for v in dif) + f"{dif.sum():>+10.1f}")
P("")
REGIONS = [(0.02,0.3),(0.3,0.8),(0.8,1.5),(1.5,6.0)]
P("vtilde-region decomposition of (model - Newton) lnL, per s-bin:")
P("  regions: " + ", ".join(f"R{k}={r[0]}-{r[1]}" for k,r in enumerate(REGIONS)))
for name in res:
    if name == 'Newton': continue
    P(f"  {name}:")
    for bi, b in enumerate(SBINS):
        pN = res['Newton'][1][bi]; pM = res[name][1][bi]
        terms = data_counts[bi]*np.log(pM/pN)
        row = []
        for r in REGIONS:
            m = (VCEN>=r[0])&(VCEN<r[1])
            row.append(terms[m].sum())
        P(f"    {b[0]}-{b[1]}kAU: " + "  ".join(f"R{k}={v:+.1f}" for k,v in enumerate(row)))
P("")
P("Shape check — fraction of counts per vtilde region (data vs models), per s-bin:")
for bi, b in enumerate(SBINS):
    P(f"  {b[0]}-{b[1]}kAU (N={int(data_counts[bi].sum())}):")
    fr_d = []
    for r in REGIONS:
        m = (VCEN>=r[0])&(VCEN<r[1])
        fr_d.append(data_counts[bi][m].sum()/data_counts[bi].sum())
    P("    data:      " + "  ".join(f"R{k}={v:.4f}" for k,v in enumerate(fr_d)))
    for name in res:
        pM = res[name][1][bi]
        fr = [pM[(VCEN>=r[0])&(VCEN<r[1])].sum() for r in REGIONS]
        P(f"    {name:<10}" + "  ".join(f"R{k}={v:.4f}" for k,v in enumerate(fr)))
P("")
P("Pearson chi2/dof per bin (multinomial, dof=39):")
for name in res:
    row = []
    for bi in range(4):
        n = data_counts[bi]; Nb = n.sum(); p = res[name][1][bi]
        row.append(np.sum((n-Nb*p)**2/(Nb*p))/39)
    P(f"  {name:<10}" + "  ".join(f"{b[0]}-{b[1]}: {v:7.2f}" for b,v in zip(SBINS,row)))

with open('data/stage3d_autopsy.txt','w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage3d_autopsy.txt")
