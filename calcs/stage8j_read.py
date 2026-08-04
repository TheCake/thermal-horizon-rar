"""STAGE 8J reader — THE WOBBLE SATURATION INSTRUMENT (pre-reg
2cec321, committed BEFORE any run; credence FROZEN every branch).

Reads the wsat cubes (10-dim), G8J-1 arithmetic identity vs the esec
(1.0, 0.95) slice, posteriors on w0/fcomp/sq/fpm, alpha/dN, the Dwob
analogue, and THE CENSUS FORWARD at the repaired MAP cell with the
SATURATED wobble (mh untouched), 8H convention at (band 9, cliff 2).
Bars: T-DEMANDED / T-REFUSED / T-PARTIAL + CENSUS-REOPENED;
T-REFUSED = the pre-committed stop rule (reviewer round, no 4th
dial).  Output: data/stage8j_read.txt
"""
import os
import numpy as np, time
from astropy.io import fits
from scipy.stats import poisson

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

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
NOBS, NHI = 9, 2

MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])

N = 500_000
A_GRID  = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID  = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FPM_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW_GRID = np.array([0.7, 1.0, 1.4])
W0SAT_GRID = np.array([0.1, 0.2, 0.4, 0.8, 1e9])
SQ_GRID = np.array([0.0, 0.1, 0.2, 0.3])
SEEDS = (31, 101)
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2

P("8J THE WOBBLE SATURATION READ (pre-reg 2cec321; bars locked "
  "before any run; credence FROZEN every branch)")
P("")
pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
assert str(pz['version']).startswith('v2c')
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
gband = pz['conv_band']/0.30
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

def e_of_x(p, eta, wr, ein, erf):
    al = np.interp(np.log10(p['a_s']), np.log10([100,500,1000,50000]),
                   [0.6*ein, 1.0*ein, eta, eta])
    e_pow = 0.95*p['u_e']**(1/(1+al))
    e_rad = erf+(0.995-erf)*p['u_e']
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

def band_mu(p, o, fcm, kw, fpm, sq, w0):
    """8H band_mu + the 8J saturation: w_eff = (w0/4.74047) *
    tanh(w*4.74047/w0) per companion; mh (dynamical mass) untouched;
    w0 >= 1e8 = the untransformed identity."""
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
            wk = c['w'][idx]
            if w0 < 1e8:
                wk = (w0/4.74047)*np.tanh(wk*4.74047/w0)
            cvp += act*wk*c['wd'][idx,0]
            cvq += act*wk*c['wd'][idx,1]
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

def jointp(mu_b, mu_c):
    return float(poisson.pmf(NOBS, max(mu_b, 1e-12))
                 * poisson.pmf(NHI, max(mu_c, 1e-12)))

def refine(x, y):
    i = int(np.nanargmax(y)); v = x[i]
    if 0 < i < len(x)-1:
        c2, c1, _ = np.polyfit(x[i-1:i+2], y[i-1:i+2], 2)
        if c2 < 0:
            v = -c1/(2*c2)
    return v

def read_any(cb, axes_post=()):
    nd = cb.ndim
    m0 = np.nanmax(cb)
    ex = np.exp(np.nan_to_num(cb-m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=tuple(range(1, nd))),
                           1e-300)) + m0
    def post(ax):
        o = ex.sum(axis=tuple(i for i in range(nd) if i != ax))
        return o/o.sum()
    return (refine(A_GRID, lm), float(lm.max()-lm[0]), lm,
            {ax: post(ax) for ax in axes_post})

# ---------------- reads ---------------------------------------------------
t0 = time.time()
res = {}
g1_ok, pend = True, []
for seed in SEEDS:
    pop = None
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        cp = f'data/stage7j_cube_full_wsat_{seed}_{law}.npy'
        ep = f'data/stage7j_cube_full_esec_{seed}_{law}.npy'
        if not os.path.exists(cp):
            pend.append(f'{law}@{seed}')
            continue
        c10 = np.load(cp)
        assert c10.shape == (5, 2, 5, 6, 1, 2, 6, 3, 5, 4), c10.shape
        cb10 = (c10
                + prior_eta.reshape((1, 2) + (1,)*8)
                + LNPI.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1, 1)))
        es = (np.load(ep)[:, :, :, 1, 2]
              + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
              + LNPI.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1)))
        dq = float(np.nanmax(np.abs(cb10[:, :, :, :, :, :, :, :, 4]
                                    - es)))
        ok1 = dq <= 1e-9
        g1_ok &= ok1
        a_id, dn_id, _, po_id = read_any(es, axes_post=(3,))
        a_fr, dn_fr, _, po = read_any(cb10, axes_post=(3, 6, 8, 9))
        cks = []
        for ki in (0, 2):
            sl = cb10[:, :, :, :, :, :, :, ki:ki+1]
            cks.append(float(np.nanmax(sl) - np.nanmax(sl[:, :, :, 3:])))
        dwob = cks[1]-cks[0]
        bx = np.unravel_index(np.nanargmax(cb10), cb10.shape)
        al = A_GRID[bx[0]]; eta = E_GRID[bx[1]]; wr = WR_GRID[bx[2]]
        fcm = FCOMP_GRID[bx[3]]; fpm = FPM_GRID[bx[6]]
        kw = KW_GRID[bx[7]]; w0 = W0SAT_GRID[bx[8]]; sq = SQ_GRID[bx[9]]
        if pop is None:
            pop = build_pop(seed)
        e_f = e_of_x(pop, eta, wr, 1.0, 0.95)
        tab_a = 1.0 + al*(TAB-1.0)
        vp_f = vp_c(pop, e_f, tab_a)
        o_f = run(pop['a_s'], e_f, pop['psi0'], pop['f_ip'], pop['M_s'],
                  pop['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_a,
                  lny0=LNY0, dlny=DLNY, vp=vp_f)
        mu_r, mh_r = band_mu(pop, o_f, fcm, kw, fpm, sq, w0)
        mu_o, mh_o = band_mu(pop, o_f, fcm, kw, fpm, sq, 1e9)
        res[(law, seed)] = dict(
            a_id=a_id, dn_id=dn_id, a_fr=a_fr, dn_fr=dn_fr,
            pw0=po[8], pfc=po[3], pfp=po[6], psq=po[9], dwob=dwob,
            cen=(mu_r, mh_r, jointp(mu_r, mh_r)),
            cen_off=(mu_o, mh_o, jointp(mu_o, mh_o)))
        P(f"[{law} {seed}] G8J-1 max|off-esec| = {dq:.2e} -> "
          f"{'PASS' if ok1 else 'FAIL'}")
        P(f"[{law} {seed}] OFF:  a={a_id:.2f} dN={dn_id:+.1f}")
        P(f"[{law} {seed}] WSAT: a={a_fr:.2f} dN={dn_fr:+.1f} "
          f"P(w0)={np.round(po[8], 2).tolist()} "
          f"P(fcomp)={np.round(po[3], 2).tolist()}")
        P(f"[{law} {seed}]       P(sq)={np.round(po[9], 2).tolist()} "
          f"P(fpm)={np.round(po[6], 2).tolist()} Dwob={dwob:+.1f}")
        P(f"[{law} {seed}] MAP cell: alpha={al}, eta={eta}, wr={wr}, "
          f"fcomp={fcm}, fpm={fpm}, kw={kw}, w0={w0}, sq={sq}")
        P(f"[{law} {seed}] CENSUS repaired: mu=({mu_r:.2f}, {mh_r:.2f}) "
          f"jointP={jointp(mu_r, mh_r):.2e} | saturation-off: "
          f"mu=({mu_o:.2f}, {mh_o:.2f}) jointP={jointp(mu_o, mh_o):.2e}")
        P("")

if pend:
    P(f"INCOMPLETE - pending cubes: {pend}; verdict PENDING")
elif not g1_ok:
    P("G8J-1 FAIL - reader wiring suspect; DO NOT QUOTE")
else:
    P("G8J-1 4/4 PASS (saturation-off slice = esec (1.0,0.95) slice, "
      "arithmetic)")
    P("")
    branch = {}
    for law in ('simple', 'BE'):
        pw0 = np.mean([res[(law, s)]['pw0'] for s in SEEDS], axis=0)
        pfin = float(pw0[:-1].sum())
        mode_i = int(np.argmax(pw0))
        da = float(np.mean([res[(law, s)]['a_fr']-res[(law, s)]['a_id']
                            for s in SEEDS]))
        br = ('T-DEMANDED' if pfin >= 0.90 else
              'T-REFUSED' if pw0[-1] >= 0.50 else 'T-PARTIAL')
        edge = mode_i == 0 and pw0[0] >= 0.5
        branch[law] = (br, pfin, mode_i, da, edge)
        P(f"BAR [{law}]: P(w0<inf)={pfin:.2f}, mode="
          f"{W0SAT_GRID[mode_i]}, d_alpha={da:+.3f}"
          + (' MATERIAL' if abs(da) > 0.11 else '')
          + (' [EDGE at 0.1]' if edge else '') + f" -> {br}")
    reop = sum(1 for r in res.values() if r['cen'][2] >= 1e-3)
    P(f"CENSUS-REOPENED count (jointP >= 1e-3 at the repaired cell): "
      f"{reop}/4 -> {'REOPENED' if reop >= 3 else 'NOT reopened'}")
    P("")
    brs = {branch[l][0] for l in ('simple', 'BE')}
    v = (('T-DEMANDED: the sky accepts the saturation; the tail-shape '
          'repair is DATA-SUPPORTED; the POWERED ROUND revives as the '
          'next decider')
         if brs == {'T-DEMANDED'} else
         ('T-REFUSED: the third distribution-level refusal; the wobble '
          'sector goes to the REVIEWER ROUND with the 8H/8I-a/8J '
          'scorecard (pre-committed stop rule - no fourth dial)')
         if brs == {'T-REFUSED'} else
         'T-PARTIAL/mixed: decomposition reported as-is')
    if reop >= 3:
        v += ('; CENSUS-REOPENED (conditional on the powered round) - '
              'the 8H demotion lifts pending validation')
    fifth = all(res[(l, s)]['a_fr'] <= 0.2 and res[(l, s)]['dn_fr'] <= 5
                for l in ('simple', 'BE') for s in SEEDS)
    if fifth:
        v += ('; fifth-move-shaped movement observed - VERDICT '
              'DEFERRED to the powered round')
    P(f"==> 8J VERDICT (locked bars): {v}")
    P("    CREDENCE FROZEN (pre-stated, every branch).")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage8j_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8j_read.txt")
