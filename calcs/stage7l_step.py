"""Stage 7L L3 — THEIR STATISTIC ON OUR MODEL (pre-reg in NOTES
2026-07-27): the Cookson median-vtilde step, computed on (i) our data
under the cook proxy mask (bootstrap CI), (ii) the landed operative
boost cells' forwards, (iii) the landed Newton-best forwards — all
with their vtilde < 2.5 ceiling, binned in r_sky/r_M (step = median
in [1, 3.1] / median in [0.05, 0.5]; r_M = 7030 sqrt(M_tot) AU).

7K-a lessons built in up front: G0-7L scale-free control (al == 1,
absorbers off) must give 1.000 +/- 0.010; the fitted-e-run
model-Newton baseline is the honest zero point, reported FIRST.
Bars: RECONCILED-BY-DESIGN boost-step <= 1.11 AND data consistent;
GENUINE-TENSION predicted >= 1.15 AND data <= 1.05; else GRAY.
Output: data/stage7l_step.txt
"""
import numpy as np, time
from astropy.io import fits

OUT = 'data/stage7l_step.txt'
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
RM_AU = 7030.0     # sqrt(G Msun / 1.2e-10) in AU

def load_tab(path):
    t = np.load(path); y, b = t[0][::-1], t[1][::-1]
    lny = np.log(y); lny_u = np.linspace(lny[0], lny[-1], 512)
    return lny_u, np.interp(lny_u, lny, b)
LNY_U, TAB_S = load_tab('data/efe_boost_simple_g1p2.npy')
_,     TAB_B = load_tab('data/efe_boost_be_g1p2.npy')
LNY0, DLNY = LNY_U[0], LNY_U[1]-LNY_U[0]

# --- data under the cook mask ---------------------------------------------
d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
cook = np.load('data/stage7l_cookmask.npy')
plx1, plx2 = d['parallax1'], d['parallax2']
sep = d['sep_AU']
G1m, G2m = d['phot_g_mean_mag1'], d['phot_g_mean_mag2']
plx = 0.5*(plx1+plx2)
MG1 = G1m+5*np.log10(np.maximum(plx1, 1e-6))-10
MG2 = G2m+5*np.log10(np.maximum(plx2, 1e-6))-10
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot = np.interp(MG1, MG_T, MS_T)+np.interp(MG2, MG_T, MS_T)
s_kau = sep/1e3
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
sn_ = np.maximum(np.hypot(sx_, sy_), 1e-12)
r1 = np.asarray(d['dr2_radial_velocity1'], np.float64)
r2 = np.asarray(d['dr2_radial_velocity2'], np.float64)
er1 = np.asarray(d['dr2_radial_velocity_error1'], np.float64)
er2 = np.asarray(d['dr2_radial_velocity_error2'], np.float64)
h1, h2 = np.isfinite(r1), np.isfinite(r2)
w1 = np.where(h1, 1.0/np.maximum(er1, 0.5)**2, 0.0)
w2 = np.where(h2, 1.0/np.maximum(er2, 0.5)**2, 0.0)
rvs = (np.where(h1, r1, 0.0)*w1 + np.where(h2, r2, 0.0)*w2) \
      / np.maximum(w1+w2, 1e-12)
th = sep/(2.06265e8/plx)
pmcor = rvs*th*plx/4.74047
vx = d['pmra2']-d['pmra1'] + pmcor*sx_/sn_
vy = d['pmdec2']-d['pmdec1'] + pmcor*sy_/sn_
vc_d = 0.9417*np.sqrt(Mtot/np.maximum(s_kau, 1e-9))
vt_d = (4.74047/plx*np.hypot(vx, vy))/vc_d
rrm_d = sep/(RM_AU*np.sqrt(Mtot))
e_vx = np.sqrt(d['pmra_error1']**2+d['pmra_error2']**2)
e_vy = np.sqrt(d['pmdec_error1']**2+d['pmdec_error2']**2)
sig_c = ((4.74047/plx)*np.sqrt(0.5*(e_vx**2+e_vy**2)))

BA, BB = (0.05, 0.5), (1.0, 3.1)
mA = cook & (rrm_d >= BA[0]) & (rrm_d < BA[1]) & (vt_d < 2.5)
mB = cook & (rrm_d >= BB[0]) & (rrm_d < BB[1]) & (vt_d < 2.5)
step_d = float(np.median(vt_d[mB])/np.median(vt_d[mA]))
rng = np.random.default_rng(7)
va, vb = vt_d[mA], vt_d[mB]
boots = [np.median(rng.choice(vb, len(vb)))/np.median(rng.choice(va, len(va)))
         for _ in range(400)]
lo, hi = np.percentile(boots, [16, 84])
P(f"STAGE 7L L3: data step (cook mask, vt<2.5, r/rM {BB}/{BA}) = "
  f"{step_d:.3f} (68% CI {lo:.3f}-{hi:.3f}; N = {int(mA.sum())}/"
  f"{int(mB.sum())})")

# model noise pools from the cook-selected pairs, per statistic bin
poolA = sig_c[cook & (rrm_d >= BA[0]) & (rrm_d < BA[1])]
poolB = sig_c[cook & (rrm_d >= BB[0]) & (rrm_d < BB[1])]

# --- model machinery ------------------------------------------------------
N = 500_000
A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
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
    u = rng.random(N); lo_, hi_, g = 0.15, 60.0, -0.6
    p['a_s'] = ((lo_**g+u*(hi_**g-lo_**g))**(1/g))*1e3
    p['u_e'] = rng.random(N)
    p['psi0'] = rng.random(N)*2*np.pi
    nrm = rng.normal(size=(N, 3)); nrm /= np.linalg.norm(nrm, axis=1,
                                                        keepdims=True)
    p['f_ip'] = np.sqrt(np.clip(1-nrm[:, 0]**2, 0, 1))
    p['M_s'] = 0.6+1.8*rng.random(N)
    p['uph'] = rng.random(N)
    xhat = np.zeros((N, 3)); xhat[:, 0] = 1
    ef = xhat-nrm*nrm[:, [0]]
    ef /= np.maximum(np.linalg.norm(ef, axis=1, keepdims=True), 1e-12)
    p['ef'] = ef; p['e2'] = np.cross(nrm, ef)
    los = rng.normal(size=(N, 3)); los /= np.linalg.norm(los, axis=1,
                                                        keepdims=True)
    p['los'] = los
    p['u_mix'] = rng.random(N)
    p['pkA'] = rng.integers(0, max(len(poolA), 1), N)
    p['pkB'] = rng.integers(0, max(len(poolB), 1), N)
    p['gn1'], p['gn2'] = rng.normal(size=N), rng.normal(size=N)
    M_h = 0.5*p['M_s']
    p['comp'] = {}
    for k in (1, 2):
        q = 0.1+0.9*rng.random(N)
        logP = rng.normal(5.03, 2.28, N)
        P_yr = 10**logP/365.25
        a_in = (M_h*(1+q)*P_yr**2)**(1/3)
        valid = (a_in < 130.0) & (a_in < p['a_s']/5.0)
        v_orb = 29.78*np.sqrt(M_h*(1+q)/np.maximum(a_in, 1e-3))
        S = np.minimum(1.0, P_yr/17.8)
        MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        MGs = np.interp(-np.clip(q*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        l_ = 10**(-0.4*(MGs-MGp))
        wfac = np.abs(q/(1+q) - l_/(1+l_))
        w = wfac*v_orb*S/4.74047*valid
        wd = rng.normal(size=(N, 3)); wd /= np.linalg.norm(wd, axis=1,
                                                          keepdims=True)
        p['comp'][k] = dict(w=w, wd=wd, uc=rng.random(N), mh=q*M_h*valid)
    p['gs'] = rng.normal(size=N)
    return p

def e_of(p, eta, wr, scale_free=False):
    if scale_free:
        al = np.full(N, 1.0)
    else:
        al = np.interp(np.log10(p['a_s']), np.log10([100, 500, 1000, 50000]),
                       [0.6, 1.0, eta, eta])
    e_pow = 0.95*p['u_e']**(1/(1+al))
    e_rad = 0.9+0.095*p['u_e']
    return np.where(p['u_mix'] < wr, e_rad, e_pow)

def vp_c(p, e_s, tab_a):
    a_s, M_s = p['a_s'], p['M_s']
    rp, ra = a_s*(1-e_s), a_s*(1+e_s)
    xg, wg = np.polynomial.legendre.leggauss(32)
    lo_, hi_ = np.log(rp), np.log(ra)
    lr = 0.5*(hi_-lo_)[:, None]*xg[None, :] + 0.5*(hi_+lo_)[:, None]
    r = np.exp(lr)
    gN = GM*M_s[:, None]/r**2
    bst = np.interp(np.log(gN/A0_CAN), LNY_U, tab_a, right=1.0)
    dPhi = np.sum(wg[None, :]*bst*gN*r, axis=1)*0.5*(hi_-lo_)
    return np.sqrt(np.maximum(2*dPhi/(1-(rp/ra)**2), 0))

def step_of(p, o, fcm, kw, fpm, sq, noise=True):
    ef, e2, los = p['ef'], p['e2'], p['los']
    s3 = o[:, 0, None]*ef+o[:, 1, None]*e2
    v3 = o[:, 2, None]*ef+o[:, 3, None]*e2
    ssky = s3-los*np.sum(s3*los, axis=1, keepdims=True)
    vsky = v3-los*np.sum(v3*los, axis=1, keepdims=True)
    smag = np.linalg.norm(ssky, axis=1)
    b1 = ssky/np.maximum(smag[:, None], 1e-12)
    b2 = np.cross(los, b1)
    b2 /= np.maximum(np.linalg.norm(b2, axis=1, keepdims=True), 1e-12)
    vpar = np.sum(vsky*b1, axis=1)
    vper = np.sum(vsky*b2, axis=1)
    s_kau_m = smag/1e3
    rrm = smag/(RM_AU*np.sqrt(p['M_s']))
    meds = {}
    for tag, (blo, bhi), pool, pick in (('A', BA, poolA, p['pkA']),
                                        ('B', BB, poolB, p['pkB'])):
        idx = np.where((rrm >= blo) & (rrm < bhi)
                       & (s_kau_m >= 1.0) & (s_kau_m < 30.0))[0]
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = pool[pick[idx] % len(pool)]/4.74047
        cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
        mh_tot = np.zeros(len(idx))
        for k in (1, 2):
            c = p['comp'][k]
            act = c['uc'][idx] < fcm
            mh_tot += act*c['mh'][idx]
            cvp += act*c['w'][idx]*c['wd'][idx, 0]
            cvq += act*c['w'][idx]*c['wd'][idx, 1]
        boost = np.sqrt(1+mh_tot/p['M_s'][idx])
        vp_b = vpar[idx] + kw*cvp
        vq_b = vper[idx] + kw*cvq
        if noise:
            vp_n = vp_b*boost + p['gn1'][idx]*sg0*fpm
            vq_n = vq_b*boost + p['gn2'][idx]*sg0*fpm
        else:
            vp_n = vp_b*boost
            vq_n = vq_b*boost
        vmag = np.hypot(vp_n, vq_n)
        keep = vmag*4.74047 <= (2.978/np.sqrt(s_kau_m[idx])
                                + 2.8284*sg0*4.74047)
        vtn = (vmag/vc)[keep]
        if sq != 0.0:
            vtn = vtn*np.exp(sq*p['gs'][idx][keep])
        vtn = vtn[vtn < 2.5]          # their ceiling
        meds[tag] = float(np.median(vtn))
    return meds['B']/meds['A']

t0 = time.time()
for seed in (31, 101):
    P(f"\n== seed {seed} ==")
    p = build_pop(seed)
    # G0-7L: scale-free control, absorbers off
    e_sf = e_of(p, 1.05, 0.2, scale_free=True)
    o_sf = run(p['a_s'], e_sf, p['psi0'], p['f_ip'], p['M_s'], p['uph'],
               8, 2500, 1)
    g0 = step_of(p, o_sf, 0.0, 1.0, 1.2, 0.0, noise=False)
    P(f"G0-7L scale-free control: step = {g0:.3f} -> "
      f"{'PASS' if 0.990 <= g0 <= 1.010 else 'FAIL'} [0.990, 1.010]")
    # cells from the FULL-sample landed cubes (the operative model)
    cells = {}
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        cw = np.load(f'data/stage7j_cube_full_photow3_{seed}_{law}.npy')
        cb = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1)) \
             + LNPI.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
        if law == 'simple':
            bx = np.unravel_index(np.nanargmax(cb[0]), cb[0].shape)
            cells['newton-best'] = (0.0, TAB) + bx
        bxa = np.unravel_index(np.nanargmax(cb), cb.shape)
        cells[f'boost-{law}'] = (A_GRID[bxa[0]], TAB) + bxa[1:]
    for tag, cell in cells.items():
        al, TAB, ei, wi, fi, ci, yi, pi, ki, si = cell
        eta, wr = E_GRID[ei], WR_GRID[wi]
        fcm, fpm = FCOMP_GRID[fi], FPM_GRID[pi]
        kw, sq = KW_GRID[ki], SQ_GRID[si]
        e_s = e_of(p, eta, wr)
        if al > 0:
            tab_a = 1.0 + al*(TAB-1.0)
            vp = vp_c(p, e_s, tab_a)
            o = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'],
                    p['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_a,
                    lny0=LNY0, dlny=DLNY, vp=vp)
        else:
            o = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'],
                    p['uph'], 8, 2500, 1)
        if tag == 'newton-best':
            base = step_of(p, o, 0.0, 1.0, 1.2, 0.0)
            P(f"e-run Newton baseline (absorbers off, fpm=1.2): "
              f"step = {base:.3f} (the honest zero point)")
        st = step_of(p, o, fcm, kw, fpm, sq)
        P(f"{tag:<12} alpha={al}: step = {st:.3f} "
          f"(fcomp={fcm}, fpm={fpm}, kw={kw}, sq={sq})")

P(f"\ndone ({(time.time()-t0)/60:.1f} min)")
