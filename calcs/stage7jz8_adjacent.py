"""Stage 7J-z8 — THE ADJACENT STATISTICS vs THE TWIN-FORCED WORLD
(pre-reg in NOTES 2026-07-28, committed before execution).

The 7J-z7 MATERIAL: forced fcomp >= 0.35 under twin-t5 collapses alpha
via the hidden-MASS channel.  That channel is ratio-flat in separation
(cannot produce an s-rising median) and inflates the pericenter pile
(must flood the census overshoot).  This instrument reads the
twin-forced cell (argmax of the forced slice of the qt5 cubes, per
law, seed 31) against both model-light statistics:
  MEDIAN leg: R >= 1.052 reproduces the anchor / <= 1.030 REJECTED /
  GRAY between.
  CENSUS leg: P(<= 2 obs overshoot | mu) < 1e-3 REJECTED / >= 0.01
  consistent / GRAY between.
VERDICT: EXPOSURE-CONTAINED if >= 1 leg REJECTS in both laws;
EXPOSURE-STANDS if both legs consistent in either law.
Output: data/stage7jz8_adjacent.txt
"""
import csv
import numpy as np, time
from astropy.io import fits
from scipy.stats import poisson

OUT = 'data/stage7jz8_adjacent.txt'
open(OUT, 'w').close()
def P(s):
    print(s, flush=True)
    with open(OUT, 'a') as f:
        f.write(s+"\n")

src = open('calcs/stage2b_population.py').read()
ns = {}
exec(src.split('# ---------- validation')[0], ns)
run = ns['run']
GM = 4*np.pi**2; A0_CAN = 7.99e-7

def load_tab(path):
    t = np.load(path); y, b = t[0][::-1], t[1][::-1]
    lny = np.log(y); lny_u = np.linspace(lny[0], lny[-1], 512)
    return lny_u, np.interp(lny_u, lny, b)
LNY_U, TAB_S = load_tab('data/efe_boost_simple_g1p2.npy')

# --- catalog (7K machinery pieces) ----------------------------------------
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
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
s_d = sep[ok]/1e3
e_vx = np.sqrt(d['pmra_error1']**2+d['pmra_error2']**2)
e_vy = np.sqrt(d['pmdec_error1']**2+d['pmdec_error2']**2)
sig_ok = ((4.74047/plx)*np.sqrt(0.5*(e_vx**2+e_vy**2)))[ok]
SBINS = [(0.2,2),(2,6),(6,20),(20,50)]
NDATA, POOLS = [], []
for b in SBINS:
    m = (s_d>=b[0])&(s_d<b[1])
    NDATA.append(int(m.sum())); POOLS.append(sig_ok[m])
B2C = [(0.2, 2.0), (6.0, 30.0)]
pool2c = {b: sig_ok[(s_d >= b[0]) & (s_d < b[1])] for b in B2C}
# data s-distributions for the median reweighting
SUBS = {}
for b in B2C:
    edges = np.logspace(np.log10(b[0]), np.log10(b[1]), 11)
    m = (s_d >= b[0]) & (s_d < b[1])
    h, _ = np.histogram(s_d[m], bins=edges)
    SUBS[b] = (edges, h/max(h.sum(), 1))
NOBS, NHI = 9, 2   # the corrected census counts (7K-b)

# --- grids + anchor -------------------------------------------------------
N = 500_000
A_GRID  = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID  = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FC0_GRID = np.array([0.10])
FFLY_GRID = np.array([0.05, 0.10])
FPM_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW_GRID = np.array([0.7, 1.0, 1.4])
SQ_GRID = np.array([0.0, 0.1, 0.2, 0.3])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
gband = (pz['conv_band']/0.30 if 'conv_band' in pz.files
         else np.array([0.33, 1.30]))
GS = np.linspace(float(gband[0]), float(gband[1]), 25)
LNPI = np.full(len(FCOMP_GRID), -1e9)
for gi in GS:
    fh_eq = FCOMP_GRID/gi
    m = (fh_eq >= fg[sup].min()) & (fh_eq <= fg[sup].max())
    cand = np.full(len(FCOMP_GRID), -1e9)
    cand[m] = np.interp(fh_eq[m], fg, lp)
    LNPI = np.maximum(LNPI, cand)

def build_pop(seed):
    rng = np.random.default_rng(seed)
    p = {}
    u = rng.random(N); lo, hi, g = 0.15, 60.0, -0.6
    p['a_s']  = ((lo**g+u*(hi**g-lo**g))**(1/g))*1e3
    p['u_e']  = rng.random(N)
    p['psi0'] = rng.random(N)*2*np.pi
    nrm = rng.normal(size=(N,3)); nrm /= np.linalg.norm(nrm,axis=1,keepdims=True)
    p['f_ip'] = np.sqrt(np.clip(1-nrm[:,0]**2,0,1))
    p['M_s']  = 0.6+1.8*rng.random(N)
    p['uph']  = rng.random(N)
    xhat = np.zeros((N,3)); xhat[:,0]=1
    ef = xhat-nrm*nrm[:,[0]]
    ef /= np.maximum(np.linalg.norm(ef,axis=1,keepdims=True),1e-12)
    p['ef'] = ef; p['e2'] = np.cross(nrm,ef)
    los = rng.normal(size=(N,3)); los /= np.linalg.norm(los,axis=1,keepdims=True)
    p['los'] = los
    p['u_mix'] = rng.random(N)
    p['pick'] = [rng.integers(0, max(len(POOLS[bi]),1), N)
                 for bi in range(len(SBINS))]
    p['pk2c'] = {b: rng.integers(0, max(len(pool2c[b]), 1), N) for b in B2C}
    p['gn1'], p['gn2'] = rng.normal(size=N), rng.normal(size=N)
    M_h = 0.5*p['M_s']
    p['comp'] = {}
    for k in (1, 2):
        u_q = rng.random(N)
        s_ = 0.9/1.4          # the twin-t5 marginal (z7 convention)
        q = np.where(u_q < s_, 0.1 + 0.8*(u_q/s_),
                     0.9 + 0.1*((u_q-s_)/(1.0-s_)))
        logP = rng.normal(5.03, 2.28, N)
        P_yr = 10**logP/365.25
        a_in = (M_h*(1+q)*P_yr**2)**(1/3)
        valid = (a_in < 130.0) & (a_in < p['a_s']/5.0)
        v_orb = 29.78*np.sqrt(M_h*(1+q)/np.maximum(a_in,1e-3))
        S = np.minimum(1.0, P_yr/17.8)
        MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        MGs = np.interp(-np.clip(q*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        l_ = 10**(-0.4*(MGs-MGp))
        wfac = np.abs(q/(1+q) - l_/(1+l_))
        w = wfac*v_orb*S/4.74047*valid
        wd = rng.normal(size=(N,3)); wd /= np.linalg.norm(wd,axis=1,keepdims=True)
        p['comp'][k] = dict(w=w, wd=wd, uc=rng.random(N), mh=q*M_h*valid)
    p['gs'] = rng.normal(size=N)
    return p

def e_of(p, eta, wr):
    al = np.interp(np.log10(p['a_s']), np.log10([100,500,1000,50000]),
                   [0.6, 1.0, eta, eta])
    e_pow = 0.95*p['u_e']**(1/(1+al))
    e_rad = 0.9+0.095*p['u_e']
    return np.where(p['u_mix'] < wr, e_rad, e_pow)

def project(p, o):
    ef, e2, los = p['ef'], p['e2'], p['los']
    s3 = o[:,0,None]*ef+o[:,1,None]*e2
    v3 = o[:,2,None]*ef+o[:,3,None]*e2
    ssky = s3-los*np.sum(s3*los,axis=1,keepdims=True)
    vsky = v3-los*np.sum(v3*los,axis=1,keepdims=True)
    smag = np.linalg.norm(ssky,axis=1)
    b1 = ssky/np.maximum(smag[:,None],1e-12)
    b2 = np.cross(los, b1)
    b2 /= np.maximum(np.linalg.norm(b2,axis=1,keepdims=True),1e-12)
    return smag, np.sum(vsky*b1,axis=1), np.sum(vsky*b2,axis=1)

def channels(p, idx, fcm, kw):
    cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
    mh_tot = np.zeros(len(idx))
    for k in (1, 2):
        c = p['comp'][k]
        act = c['uc'][idx] < fcm
        mh_tot += act*c['mh'][idx]
        cvp += act*c['w'][idx]*c['wd'][idx,0]
        cvq += act*c['w'][idx]*c['wd'][idx,1]
    return np.sqrt(1+mh_tot/p['M_s'][idx]), kw*cvp, kw*cvq

def median_R(p, smag, vpar, vper, fcm, kw, fpm, sq):
    s_kau = smag/1e3
    meds = {}
    for b in B2C:
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        pool = pool2c[b]
        sg0 = pool[p['pk2c'][b][idx] % len(pool)]/4.74047
        boost, cvp, cvq = channels(p, idx, fcm, kw)
        vp_n = (vpar[idx]+cvp)*boost + p['gn1'][idx]*sg0*fpm
        vq_n = (vper[idx]+cvq)*boost + p['gn2'][idx]*sg0*fpm
        vmag = np.hypot(vp_n, vq_n)
        keep = vmag*4.74047 <= (2.978/np.sqrt(s_kau[idx])
                                + 2.8284*sg0*4.74047)
        vtn = (vmag/vc)[keep]
        if sq != 0.0:
            vtn = vtn*np.exp(sq*p['gs'][idx][keep])
        sk = s_kau[idx][keep]
        edges, dfr = SUBS[b]
        mi = np.clip(np.digitize(sk, edges)-1, 0, 9)
        mh_, _ = np.histogram(sk, bins=edges)
        wgt = dfr[mi]/np.maximum((mh_/max(mh_.sum(),1))[mi], 1e-9)
        o_ = np.argsort(vtn)
        cw = np.cumsum(wgt[o_]); cw /= cw[-1]
        meds[b] = float(np.interp(0.5, cw, vtn[o_]))
    return meds[B2C[1]]/meds[B2C[0]]

def census_mu(p, smag, vpar, vper, fcm, kw, fpm, sq):
    s_kau = smag/1e3
    mu, mu_hi = 0.0, 0.0
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500:
            continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = POOLS[bi][p['pick'][bi][idx] % len(POOLS[bi])]/4.74047
        boost, cvp, cvq = channels(p, idx, fcm, kw)
        vp_n = (vpar[idx]+cvp)*boost + p['gn1'][idx]*sg0*fpm
        vq_n = (vper[idx]+cvq)*boost + p['gn2'][idx]*sg0*fpm
        vmag = np.hypot(vp_n, vq_n)
        keep = vmag*4.74047 <= (2.978/np.sqrt(s_kau[idx])
                                + 2.8284*sg0*4.74047)
        vtn = (vmag/vc)[keep]
        gmn = np.degrees(np.arccos(np.clip(
            np.abs(vp_n[keep])/np.maximum(vmag[keep],1e-12), 0, 1)))
        vts = vtn*np.exp(sq*p['gs'][idx][keep])
        nk = max(int(keep.sum()), 1)
        inb = float(np.sum((gmn >= 75) & (vts >= 1.414) & (vts < 1.67)))
        inh = float(np.sum((gmn >= 75) & (vts >= 1.67) & (vts < 2.2)))
        mu += NDATA[bi]*inb/nk
        mu_hi += NDATA[bi]*inh/nk
    return mu, mu_hi

t0 = time.time()
p = build_pop(31)
verdicts = []
for law in ('simple', 'BE'):
    cw = np.load(f'data/stage7j_cube_full_qt5_31_{law}.npy')
    cb = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1)) \
         + LNPI.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
    sub = cb[:, :, :, 3:]
    bx = np.unravel_index(np.nanargmax(sub), sub.shape)
    ai, ei, wi, fi, ci, yi, pi, ki, si = bx
    al = A_GRID[ai]; eta, wr = E_GRID[ei], WR_GRID[wi]
    fcm = FCOMP_GRID[fi+3]; fpm = FPM_GRID[pi]
    kw, sq = KW_GRID[ki], SQ_GRID[si]
    P(f"\n[{law}] twin-forced cell: alpha={al}, eta={eta}, wr={wr}, "
      f"fcomp={fcm}, fpm={fpm}, kw={kw}, sq={sq}")
    assert al == 0.0, "the forced argmax is expected at alpha=0 (z7)"
    e_s = e_of(p, eta, wr)
    o = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'], p['uph'],
            8, 2500, 1)
    smag, vpar, vper = project(p, o)
    R = median_R(p, smag, vpar, vper, fcm, kw, fpm, sq)
    med = ('REJECTED' if R <= 1.030 else
           'REPRODUCES' if R >= 1.052 else 'GRAY')
    P(f"[{law}] MEDIAN leg: R_forced = {R:.3f} vs anchor 1.078 "
      f"(CI 1.052-1.103) -> {med}")
    mu, mu_hi = census_mu(p, smag, vpar, vper, fcm, kw, fpm, sq)
    pc = float(poisson.cdf(NHI, mu_hi))
    cen = ('REJECTED' if pc < 1e-3 else
           'CONSISTENT' if pc >= 0.01 else 'GRAY')
    P(f"[{law}] CENSUS leg: mu_band = {mu:.1f} (obs {NOBS}), "
      f"mu_overshoot = {mu_hi:.1f} (obs {NHI}), P(<={NHI}) = {pc:.2e} "
      f"-> {cen}")
    verdicts.append((law, med == 'REJECTED', cen == 'REJECTED'))

rej_both = all(m or c for _, m, c in verdicts)
cons_any = any((not m) and (not c) for _, m, c in verdicts)
v = ('EXPOSURE-CONTAINED' if rej_both else
     'EXPOSURE-STANDS' if cons_any else 'PARTIAL')
P(f"\n==> 7J-z8 VERDICT: {v} "
  f"({[(l, 'med-rej' if m else '', 'cen-rej' if c else '') for l, m, c in verdicts]})")
P(f"done ({(time.time()-t0)/60:.1f} min)")
