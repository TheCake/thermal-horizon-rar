"""
STAGE 3L: eccentricity measurement via the v-angle method (Hwang+22 style) —
TODO #15, promoted. gamma = angle in the sky plane between the pair's
separation vector and its relative proper-motion vector, folded to [0,90] deg
(orbit-direction sign is unobservable). Circular orbits -> gamma ~ 90 deg;
radial-heavy populations -> small gamma. KEY PROPERTY: gamma uses only
DIRECTIONS, so it is immune to the mass normalization — the clean channel to
test whether the e-distribution shape is the missing sigma_m broadening.

Forward model: our orbit engine (Newton and BE alpha=1), eta grid, per-pair
matched astrometric noise, same S/N cut on data and model.
Gates: G1 circular population concentrates at 90 deg; G2 recover eta=1.3 from
a mock. (Angle statistics are where this project caught its first bug — the
2/pi double-fold. Suspicious round numbers remain confessions.)

Note discovered while matching noise: the main pipeline applies sigv per sky
component; the data's per-component error is sigv/sqrt(2). Conservative
(over-noises the model; cannot create the too-narrow-model problem — it works
against it). Flagged in NOTES; here the correct sigv/sqrt(2) is used.
Writes data/stage3l_summary.txt.
"""
import time
import numpy as np
from astropy.io import fits

src = open('calcs/stage2b_population.py').read()
ns = {}
exec(src.split('# ---------- validation')[0], ns)
run = ns['run']
GM = 4*np.pi**2; A0_CAN = 7.99e-7

def load_tab(path):
    t = np.load(path); y, b = t[0][::-1], t[1][::-1]
    lny = np.log(y); lny_u = np.linspace(lny[0], lny[-1], 512)
    return lny_u, np.interp(lny_u, lny, b)
LNY_U, TAB_B = load_tab('data/efe_boost_be.npy')
LNY0, DLNY = LNY_U[0], LNY_U[1]-LNY_U[0]

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

# --- data-side gamma ---
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx = dra*np.cos(dec_m); sy = d['dec2']-d['dec1']
vx = d['pmra2']-d['pmra1']; vy = d['pmdec2']-d['pmdec1']
e_vx = np.sqrt(d['pmra_error1']**2+d['pmra_error2']**2)
e_vy = np.sqrt(d['pmdec_error1']**2+d['pmdec_error2']**2)
vmag = np.hypot(vx, vy)
sn = vmag/np.sqrt(0.5*(e_vx**2+e_vy**2))       # per-component S/N
smag = np.hypot(sx, sy)
cosg = np.abs(sx*vx+sy*vy)/np.maximum(smag*vmag, 1e-12)
gam_d = np.degrees(np.arccos(np.clip(cosg, 0, 1)))
s_d = sep/1e3
SNCUT = 3.0
usable = ok & (sn > SNCUT)
SBINS = [(0.2,2),(2,6),(6,20),(20,50)]
GEDGE = np.linspace(0, 90, 10)
L = []
def P(s):
    print(s); L.append(s)
P(f"STAGE 3L v-angle: {ok.sum()} pairs; S/N>{SNCUT}: {usable.sum()}")
data_g = []
for b in SBINS:
    m = usable & (s_d>=b[0]) & (s_d<b[1])
    h,_ = np.histogram(gam_d[m], bins=GEDGE)
    data_g.append(h.astype(float))
    P(f"  {b[0]}-{b[1]} kAU: N={int(h.sum())}, mean gamma="
      f"{gam_d[m].mean():.1f} deg")

# --- model side ---
# noise: per-component sigma in km/s, resampled from the data per s-bin
noise_pool = []
for b in SBINS:
    m = ok & (s_d>=b[0]) & (s_d<b[1])
    noise_pool.append((4.74047/plx[m])*np.sqrt(0.5*(e_vx[m]**2+e_vy[m]**2)))

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
pickers = [rng.integers(0, max(len(p_),1), N) for p_ in noise_pool]
gn1, gn2 = rng.normal(size=N), rng.normal(size=N)

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

def gamma_hists(o, noisy=True):
    """returns per-bin gamma histograms with matched noise + S/N cut"""
    s3 = o[:,0,None]*ef+o[:,1,None]*e2
    v3 = o[:,2,None]*ef+o[:,3,None]*e2
    ssky = s3-los*np.sum(s3*los,axis=1,keepdims=True)
    vsky = v3-los*np.sum(v3*los,axis=1,keepdims=True)
    smag_ = np.linalg.norm(ssky,axis=1)
    b1 = ssky/np.maximum(smag_[:,None],1e-12)
    b2 = np.cross(los, b1)
    b2 /= np.maximum(np.linalg.norm(b2,axis=1,keepdims=True),1e-12)
    vpar = np.sum(vsky*b1,axis=1)*4.74047   # km/s (model units AU/yr *4.74)
    vper = np.sum(vsky*b2,axis=1)*4.74047
    s_kau = smag_/1e3
    hists = []
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx)==0 or len(noise_pool[bi])==0:
            hists.append(np.ones(9)); continue
        sg = noise_pool[bi][pickers[bi][idx] % len(noise_pool[bi])]
        if noisy:
            vp_ = vpar[idx]+gn1[idx]*sg; vq_ = vper[idx]+gn2[idx]*sg
            keep = np.hypot(vp_,vq_)/sg > SNCUT
        else:
            vp_, vq_ = vpar[idx], vper[idx]
            keep = np.ones(len(idx), bool)
        cg = np.abs(vp_[keep])/np.maximum(np.hypot(vp_[keep],vq_[keep]),1e-12)
        gm = np.degrees(np.arccos(np.clip(cg,0,1)))
        h,_ = np.histogram(gm, bins=GEDGE)
        hists.append(np.maximum(h,1e-3))
    return hists

def lnl_bins(hm, hd):
    out = np.zeros(4)
    for bi in range(4):
        p = hm[bi]/hm[bi].sum(); p = np.maximum(p, 1e-5); p /= p.sum()
        out[bi] = np.sum(hd[bi]*np.log(p))
    return out

# gate G1: circular orbits, no noise -> gamma at 90
e0 = np.full(N, 0.001)
o = run(a_s, e0, psi0, f_ip, M_s, uph, 8, 2500, 1)
h = gamma_hists(o, noisy=False)
frac_hi = h[1][-2:].sum()/h[1].sum()
P(f"G1 circular gate: fraction of gamma in [70,90] (2-6 kAU, noiseless) = "
  f"{frac_hi:.3f} (expect >0.9)")

ETA_GRID = np.array([0.6, 0.8, 1.05, 1.3, 1.6, 2.0, 2.4])
models = {}
for law in ("newton", "BE1"):
    hh = []
    for eta in ETA_GRID:
        e_s = e_of(eta)
        if law == "newton":
            o = run(a_s, e_s, psi0, f_ip, M_s, uph, 8, 2500, 1)
        else:
            o = run(a_s, e_s, psi0, f_ip, M_s, uph, 8, 2500, 5,
                    a0=A0_CAN, tab=TAB_B, lny0=LNY0, dlny=DLNY,
                    vp=vp_c(e_s, TAB_B))
        hh.append(gamma_hists(o))
    models[law] = hh
    print(f"{law} eta grid done")

# gate G2: mock recovery (BE1 eta=1.3 histograms as fake data)
mock = [models['BE1'][3][bi]*data_g[bi].sum()/models['BE1'][3][bi].sum()
        for bi in range(4)]
for law in ("BE1",):
    tot = np.array([lnl_bins(models[law][k], mock).sum()
                    for k in range(len(ETA_GRID))])
    P(f"G2 mock(eta=1.3) recovery [{law}]: eta profile "
      f"{np.round(tot-tot.max(),1).tolist()} -> "
      f"eta_hat={ETA_GRID[int(np.argmax(tot))]} (expect 1.3)")

# measurement: per s-bin eta from data
BNAMES = [f"{b[0]}-{b[1]}kAU" for b in SBINS]
for law in models:
    P(f"[{law}] per-bin lnL(eta) profiles (rel max), eta grid "
      f"{ETA_GRID.tolist()}:")
    for bi in range(4):
        prof = np.array([lnl_bins(models[law][k], data_g)[bi]
                         for k in range(len(ETA_GRID))])
        rel = prof - prof.max()
        imax = int(np.argmax(prof))
        # parabolic interp on log-spaced-ish grid
        ehat = ETA_GRID[imax]
        if 0 < imax < len(ETA_GRID)-1:
            x = ETA_GRID[imax-1:imax+2]; y = prof[imax-1:imax+2]
            c2, c1, _ = np.polyfit(x, y, 2)
            if c2 < 0: ehat = -c1/(2*c2)
        P(f"  {BNAMES[bi]}: eta_hat={ehat:.2f} "
          f"{'(interior)' if 0<imax<len(ETA_GRID)-1 else '(EDGE)'} "
          f"prof={np.round(rel,1).tolist()}")

with open('data/stage3l_summary.txt','w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage3l_summary.txt")
