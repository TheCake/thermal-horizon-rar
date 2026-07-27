"""Stage 7K-b — THE CENSUS LEAKAGE NULL AT THE LANDED CELL (pre-reg in
NOTES 2026-07-27, committed before execution).

The forward Newton-best sky at the landed posterior mode carries the
companion wobble + hidden-mass machinery, the demanded noise, and the
smear — its occupancy of the Newton-forbidden perpendicular band
(vtilde in [1.414, 1.67), gamma >= 75 deg) IS the companion-
marginalized leakage null.  Observed: n_obs = 9 (corrected convention;
raw-sn3 11 reported).  Bars: NULL-INTACT P <= 1e-4 / NULL-BROKEN
P >= 0.01 / GRAY.  Plus the wobble-lite unconditional per-pair bounds
(hidden-mass channel disclosed as non-boundable without photometry).
Output: data/stage7kb_census.txt
"""
import csv
import numpy as np, time
from astropy.io import fits
from scipy.stats import poisson

OUT = 'data/stage7kb_census.txt'
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
_,     TAB_B = load_tab('data/efe_boost_be_g1p2.npy')
LNY0, DLNY = LNY_U[0], LNY_U[1]-LNY_U[0]

# --- catalog: per-bin N + noise pools (marginal-pipeline convention) ------
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
s_d = sep[ok]/1e3
e_vx = np.sqrt(d['pmra_error1']**2+d['pmra_error2']**2)
e_vy = np.sqrt(d['pmdec_error1']**2+d['pmdec_error2']**2)
sig_ok = ((4.74047/plx)*np.sqrt(0.5*(e_vx**2+e_vy**2)))[ok]
SBINS = [(0.2,2),(2,6),(6,20),(20,50)]
NDATA, POOLS = [], []
for b in SBINS:
    m = (s_d>=b[0])&(s_d<b[1])
    NDATA.append(int(m.sum())); POOLS.append(sig_ok[m])
P(f"STAGE 7K-b: data N per bin = {NDATA}")

MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])

# --- the observed census + wobble-lite bounds -----------------------------
cens = [r for r in csv.DictReader(open('data/ceiling_pairs.csv'))
        if r['census_corr'] == 'True']
NOBS = len(cens)
NRAW = sum(1 for r in csv.DictReader(open('data/ceiling_pairs.csv'))
           if r['census_raw_sn3'] == 'True')
# the observed OVERSHOOT (gamma >= 75, vt_corr in [1.67, 2.2)) — the
# cliff's far side; quality-cut caveat shared with the band count
NHI = sum(1 for r in csv.DictReader(open('data/ceiling_pairs.csv'))
          if float(r['gamma_corr_deg']) >= 75
          and 1.67 <= float(r['vt_corr']) < 2.2)
P(f"observed census: corrected n_obs = {NOBS} (raw-sn3 {NRAW}); "
  f"observed overshoot [1.67, 2.2) at gamma>=75: {NHI}")
qg = np.linspace(0.1, 1.0, 91)
nsafe = 0
for r in cens:
    Ms = float(r['Mtot_Msun']); vc = float(r['vc_kms'])
    need = float(r['vt_corr']) - 1.414
    M_h = 0.5*Ms
    MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
    MGs = np.interp(-np.clip(qg*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
    l_ = 10**(-0.4*(MGs-MGp))
    wfac = np.abs(qg/(1+qg) - l_/(1+l_))
    a_in = (M_h*(1+qg)*17.8**2)**(1/3)
    v_orb = 29.78*np.sqrt(M_h*(1+qg)/a_in)
    dvt = 1.4*np.max(wfac*v_orb)/vc          # kw_max = 1.4, S = 1
    safe = dvt < need
    nsafe += safe
    P(f"  pair s={float(r['s_kAU']):6.2f}: need {need:.4f}, "
      f"wobble-max {dvt:.4f} -> {'SAFE' if safe else 'VULNERABLE'}")
P(f"wobble-lite: {nsafe}/{NOBS} SAFE "
  f"({'sentence licensed' if nsafe >= 7 else 'not licensed'}); "
  f"hidden-mass channel NOT boundable without photometry (disclosed)")

# --- model machinery (7K-a copy) ------------------------------------------
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
fhmin, fhmax = fg[sup].min(), fg[sup].max()
LNPI = np.full(len(FCOMP_GRID), -1e9)
for gi in GS:
    fh_eq = FCOMP_GRID/gi
    m = (fh_eq >= fhmin) & (fh_eq <= fhmax)
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
    p['gn1'], p['gn2'] = rng.normal(size=N), rng.normal(size=N)
    M_h = 0.5*p['M_s']
    p['comp'] = {}
    for k in (1, 2):
        q = 0.1+0.9*rng.random(N)
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

def vp_c(p, e_s, tab_a):
    a_s, M_s = p['a_s'], p['M_s']
    rp, ra = a_s*(1-e_s), a_s*(1+e_s)
    xg, wg = np.polynomial.legendre.leggauss(32)
    lo_, hi_ = np.log(rp), np.log(ra)
    lr = 0.5*(hi_-lo_)[:,None]*xg[None,:] + 0.5*(hi_+lo_)[:,None]
    r = np.exp(lr)
    gN = GM*M_s[:,None]/r**2
    bst = np.interp(np.log(gN/A0_CAN), LNY_U, tab_a, right=1.0)
    dPhi = np.sum(wg[None,:]*bst*gN*r, axis=1)*0.5*(hi_-lo_)
    return np.sqrt(np.maximum(2*dPhi/(1-(rp/ra)**2), 0))

def band_mu(p, o, fcm, kw, fpm, sq):
    """Expected census count: sum over bins of N_data * f_band, plus the
    overshoot band [1.67, 2.2) expectation."""
    ef, e2, los = p['ef'], p['e2'], p['los']
    s3 = o[:,0,None]*ef+o[:,1,None]*e2
    v3 = o[:,2,None]*ef+o[:,3,None]*e2
    ssky = s3-los*np.sum(s3*los,axis=1,keepdims=True)
    vsky = v3-los*np.sum(v3*los,axis=1,keepdims=True)
    smag = np.linalg.norm(ssky,axis=1)
    b1 = ssky/np.maximum(smag[:,None],1e-12)
    b2 = np.cross(los, b1)
    b2 /= np.maximum(np.linalg.norm(b2,axis=1,keepdims=True),1e-12)
    vpar = np.sum(vsky*b1,axis=1)
    vper = np.sum(vsky*b2,axis=1)
    s_kau = smag/1e3
    mu, mu_hi = 0.0, 0.0
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500:
            continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = POOLS[bi][p['pick'][bi][idx] % len(POOLS[bi])]/4.74047
        cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
        mh_tot = np.zeros(len(idx))
        for k in (1, 2):
            c = p['comp'][k]
            act = c['uc'][idx] < fcm
            mh_tot += act*c['mh'][idx]
            cvp += act*c['w'][idx]*c['wd'][idx,0]
            cvq += act*c['w'][idx]*c['wd'][idx,1]
        boost = np.sqrt(1+mh_tot/p['M_s'][idx])
        vp_n = (vpar[idx]+kw*cvp)*boost + p['gn1'][idx]*sg0*fpm
        vq_n = (vper[idx]+kw*cvq)*boost + p['gn2'][idx]*sg0*fpm
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
for seed in (31, 101):
    P(f"\n== seed {seed} ==")
    p = build_pop(seed)
    cw = np.load(f'data/stage7j_cube_full_photow3_{seed}_simple.npy')
    cbs = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1)) \
          + LNPI.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
    bx = np.unravel_index(np.nanargmax(cbs[0]), cbs[0].shape)
    ei, wi, fi, ci, yi, pi, ki, si = bx
    eta, wr = E_GRID[ei], WR_GRID[wi]
    fcm, fpm = FCOMP_GRID[fi], FPM_GRID[pi]
    kw, sq = KW_GRID[ki], SQ_GRID[si]
    e_s = e_of(p, eta, wr)
    o_n = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'], p['uph'],
              8, 2500, 1)
    variants = [('newton-best', fcm, kw, fpm), ('fcomp=0.2', 0.20, kw, fpm),
                ('kw=1.4', fcm, 1.4, fpm)]
    cats = []
    for tag, f_, k_, pm_ in variants:
        mu, mu_hi = band_mu(p, o_n, f_, k_, pm_, sq)
        pv = float(poisson.sf(NOBS-1, max(mu, 1e-12)))
        cat = ('NULL-INTACT' if pv <= 1e-4 else
               'NULL-BROKEN' if pv >= 0.01 else 'GRAY')
        cats.append(cat)
        P(f"{tag:<12} (fpm={pm_}, sq={sq}): mu_band = {mu:.2f}, "
          f"P(>={NOBS}) = {pv:.2e} -> {cat}  [overshoot mu = {mu_hi:.2f}, "
          f"P(<={NHI} obs | mu) = {float(poisson.cdf(NHI, mu_hi)):.2e}]")
    P(f"seed {seed} category (weakest of variants): "
      f"{'NULL-BROKEN' if 'NULL-BROKEN' in cats else 'GRAY' if 'GRAY' in cats else 'NULL-INTACT'}")
    # POST-HOC variant (labeled as such — added AFTER the first run
    # showed the fpm=3.0 corner over-leaking band AND cliff): the
    # PHYSICAL noise envelope.  Does the census null return INTACT
    # under physical noise?  Reported conditional, no bar.
    for pm_ in (1.5, 1.8):
        mu, mu_hi = band_mu(p, o_n, fcm, kw, pm_, sq)
        pv = float(poisson.sf(NOBS-1, max(mu, 1e-12)))
        P(f"POST-HOC PHYS fpm={pm_}: mu_band = {mu:.2f}, "
          f"P(>={NOBS}) = {pv:.2e}  [overshoot mu = {mu_hi:.2f}, "
          f"P(<={NHI} obs) = {float(poisson.cdf(NHI, mu_hi)):.2e}]")
    # POST-HOC attribution (labeled): the smear-off legs — the PHYS
    # runs showed mu_band ~ fpm-INDEPENDENT, i.e. the flood is the
    # gamma-blind sq smear lifting the perpendicular 1.1-1.3
    # population.  (a) Newton, sq=0, physical noise: does the census's
    # Newton-rejection power RETURN once the vacuous smear tail is
    # removed?  (b) at the end of the seed block: boost cells at
    # sq=0 — does boost-without-smear reproduce band + cliff?
    for pm_ in (1.5, 3.0):
        mu, mu_hi = band_mu(p, o_n, fcm, kw, pm_, 0.0)
        pv = float(poisson.sf(NOBS-1, max(mu, 1e-12)))
        P(f"POST-HOC Newton sq=0 fpm={pm_}: mu_band = {mu:.2f}, "
          f"P(>={NOBS}) = {pv:.2e}  [overshoot mu = {mu_hi:.2f}]")
    # boost cells (descriptive)
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        cw = np.load(f'data/stage7j_cube_full_photow3_{seed}_{law}.npy')
        cb = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1)) \
             + LNPI.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
        bxa = np.unravel_index(np.nanargmax(cb), cb.shape)
        al = A_GRID[bxa[0]]
        eta_, wr_ = E_GRID[bxa[1]], WR_GRID[bxa[2]]
        fcm_, fpm_ = FCOMP_GRID[bxa[3]], FPM_GRID[bxa[6]]
        kw_, sq_ = KW_GRID[bxa[7]], SQ_GRID[bxa[8]]
        e_b = e_of(p, eta_, wr_)
        tab_a = 1.0 + al*(TAB-1.0)
        vp = vp_c(p, e_b, tab_a)
        o_b = run(p['a_s'], e_b, p['psi0'], p['f_ip'], p['M_s'],
                  p['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_a,
                  lny0=LNY0, dlny=DLNY, vp=vp)
        mu, mu_hi = band_mu(p, o_b, fcm_, kw_, fpm_, sq_)
        P(f"boost {law} (alpha={al}): mu_band = {mu:.2f} vs obs {NOBS}; "
          f"overshoot [1.67,2.2) mu = {mu_hi:.2f}")
        mu0, mu0_hi = band_mu(p, o_b, fcm_, kw_, 1.5, 0.0)
        P(f"POST-HOC boost {law} sq=0 fpm=1.5: mu_band = {mu0:.2f} "
          f"(obs {NOBS}), overshoot mu = {mu0_hi:.2f} (obs {NHI}) "
          f"[P(>={NOBS})={float(poisson.sf(NOBS-1, max(mu0,1e-12))):.2e}, "
          f"P(<={NHI})={float(poisson.cdf(NHI, mu0_hi)):.2e}]")

P(f"\ndone ({(time.time()-t0)/60:.1f} min)")
