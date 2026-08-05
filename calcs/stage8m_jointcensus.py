"""STAGE 8M — THE JOINT-COHERENCE SCAN (the census-closure round;
pre-reg committed BEFORE any run; the reviewer's round-6 carried
item).  For every census-distinct cell of the operative
corrected-kernel (lker) cube, forward the (band, cliff) census pair
with the 8L-b reader-verbatim machinery and score
CLL = lnPois(9|mu_b) + lnPois(2|mu_h).  Products: the admissible
set (8H bar), the kinematic-price frontier, and the joint-diagnostic
posterior (CO-QUOTED-DIAGNOSTIC; double-count disclosed in NOTES).
Primary floor = per-cell MC half-count resolution; legacy 1e-12
co-read.  NO credence movement in any branch.
Output: data/stage8m_jointcensus.txt
"""
import re
import numpy as np, time
from astropy.io import fits
from scipy.stats import poisson

L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

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

# --- catalog (8L-b reader verbatim) ---------------------------------------
d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1, plx2 = d['parallax1'], d['parallax2']
eplx1, eplx2 = d['parallax_error1'], d['parallax_error2']
sep, Rch = d['sep_AU'], d['R_chance_align']
G1m, G2m = d['phot_g_mean_mag1'], d['phot_g_mean_mag2']
plx = 0.5*(plx1+plx2)
MG1 = G1m+5*np.log10(np.maximum(plx1, 1e-6))-10
MG2 = G2m+5*np.log10(np.maximum(plx2, 1e-6))-10
sigv = 4.74047/plx*np.sqrt(d['pmra_error1']**2+d['pmdec_error1']**2
                           + d['pmra_error2']**2+d['pmdec_error2']**2)
ok = (Rch < 0.01) & (plx1 > 5) & (plx2 > 5) \
   & (plx1/np.maximum(eplx1, 1e-6) > 20) & (plx2/np.maximum(eplx2, 1e-6) > 20) \
   & (np.abs(plx1-plx2) < 3*np.hypot(eplx1, eplx2)) \
   & (sep > 200) & (sep < 50000) \
   & (MG1 > 2.6) & (MG1 < 14.2) & (MG2 > 2.6) & (MG2 < 14.2) & (sigv < 0.03)
s_d_all = sep[ok]/1e3
e_vx = np.sqrt(d['pmra_error1']**2+d['pmra_error2']**2)
e_vy = np.sqrt(d['pmdec_error1']**2+d['pmdec_error2']**2)
sig_ok = ((4.74047/plx)*np.sqrt(0.5*(e_vx**2+e_vy**2)))[ok]
SBINS = [(0.2, 2), (2, 6), (6, 20), (20, 50)]
NDATA, POOLS = [], []
for b in SBINS:
    m = (s_d_all >= b[0]) & (s_d_all < b[1])
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
SQ_GRID = np.array([0.0, 0.1, 0.2, 0.3])
SEEDS = (31, 101)
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
FRX = [0, 2, 5, 10, 15, 20, 30, 40, 60, 100]
AXN = ['alpha', 'eta', 'wr', 'fcomp', 'fpm', 'kw', 'sq']
AXG = [A_GRID, E_GRID, WR_GRID, FCOMP_GRID, FPM_GRID, KW_GRID, SQ_GRID]

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

def refine(x, y):
    i = int(np.nanargmax(y)); v = x[i]
    if 0 < i < len(x)-1:
        c2, c1, _ = np.polyfit(x[i-1:i+2], y[i-1:i+2], 2)
        if c2 < 0:
            v = -c1/(2*c2)
    return v

def read9(cb, axes_post=()):
    m0 = np.nanmax(cb)
    ex = np.exp(np.nan_to_num(cb-m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=tuple(range(1, 9))),
                           1e-300)) + m0
    def post(ax):
        o = ex.sum(axis=tuple(i for i in range(9) if i != ax))
        return o/o.sum()
    return (refine(A_GRID, lm), float(lm.max()-lm[0]),
            {ax: post(ax) for ax in axes_post})

def lk_S(P_yr):
    u_ = np.pi*2.83/np.maximum(P_yr, 1e-9)
    return np.where(u_ < 1e-2, 1.0 - u_*u_/10.0,
                    3.0*(np.sin(u_) - u_*np.cos(u_))
                    / np.maximum(u_, 1e-300)**3)

def build_pop(seed, lker=True):
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
        S = lk_S(P_yr) if lker else np.minimum(1.0, P_yr/17.8)
        MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        MGs = np.interp(-np.clip(q*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        l_ = 10**(-0.4*(MGs-MGp))
        wfac = np.abs(q/(1+q) - l_/(1+l_))
        w = wfac*v_orb*S/4.74047*valid
        wd = rng.normal(size=(N,3)); wd /= np.linalg.norm(wd,axis=1,keepdims=True)
        p['comp'][k] = dict(w=w, wd=wd, uc=rng.random(N), mh=q*M_h*valid)
    p['gs'] = rng.normal(size=N)
    return p

def e_of_x(p, eta, wr):
    al = np.interp(np.log10(p['a_s']), np.log10([100,500,1000,50000]),
                   [0.6, 1.0, eta, eta])
    e_pow = 0.95*p['u_e']**(1/(1+al))
    e_rad = 0.95+0.045*p['u_e']
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
    """8L-b reader verbatim — the G8M-0/G8M-3 direct reference."""
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

def scan_block(p, o):
    """All-nuisance census forward for one orbit config.  Returns
    (mub, muh, res) with axes (FCOMP, KW, FPM, SQ) and the per-
    (f,k,p) count-resolution sum; the expression sequence is
    band_mu verbatim so the shipped-cell values agree bit-for-bit
    (G8M-3)."""
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
    NF, NK, NP, NS = (len(FCOMP_GRID), len(KW_GRID), len(FPM_GRID),
                      len(SQ_GRID))
    mub = np.zeros((NF, NK, NP, NS)); muh = np.zeros((NF, NK, NP, NS))
    res = np.zeros((NF, NK, NP))
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500:
            continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = POOLS[bi][p['pick'][bi][idx] % len(POOLS[bi])]/4.74047
        vpar_i, vper_i = vpar[idx], vper[idx]
        g1_i, g2_i = p['gn1'][idx], p['gn2'][idx]
        gs_i, sk_i = p['gs'][idx], s_kau[idx]
        for fi, fcm in enumerate(FCOMP_GRID):
            cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
            mh_tot = np.zeros(len(idx))
            for k in (1, 2):
                c = p['comp'][k]
                act = c['uc'][idx] < fcm
                mh_tot += act*c['mh'][idx]
                cvp += act*c['w'][idx]*c['wd'][idx,0]
                cvq += act*c['w'][idx]*c['wd'][idx,1]
            boost = np.sqrt(1+mh_tot/p['M_s'][idx])
            for ki, kw in enumerate(KW_GRID):
                vp_a = vpar_i+kw*cvp
                vq_a = vper_i+kw*cvq
                for pi, fpm in enumerate(FPM_GRID):
                    vp_n = vp_a*boost + g1_i*sg0*fpm
                    vq_n = vq_a*boost + g2_i*sg0*fpm
                    vmag = np.hypot(vp_n, vq_n)
                    keep = vmag*4.74047 <= (2.978/np.sqrt(sk_i)
                                            + 2.8284*sg0*4.74047)
                    vtn = (vmag/vc)[keep]
                    gmn = np.degrees(np.arccos(np.clip(
                        np.abs(vp_n[keep])/np.maximum(vmag[keep],1e-12),
                        0, 1)))
                    nk = max(int(keep.sum()), 1)
                    res[fi,ki,pi] += NDATA[bi]/nk
                    gsel = gmn >= 75
                    vg = vtn[gsel]
                    gg = gs_i[keep][gsel]
                    for si, sq in enumerate(SQ_GRID):
                        vts = vg*np.exp(sq*gg)
                        mub[fi,ki,pi,si] += NDATA[bi]*float(
                            np.sum((vts >= 1.414) & (vts < 1.67)))/nk
                        muh[fi,ki,pi,si] += NDATA[bi]*float(
                            np.sum((vts >= 1.67) & (vts < 2.2)))/nk
    return mub, muh, res

def jointp(mu_b, mu_c):
    return float(poisson.pmf(NOBS, max(mu_b, 1e-12))
                 * poisson.pmf(NHI, max(mu_c, 1e-12)))

# --- shipped anchors (parsed at runtime, the 8H pattern) ------------------
ship = {}
for ln in open('data/stage8lb_read.txt').read().splitlines():
    mm = re.match(r"\[(simple|BE) (\d+)\] MAP cell: alpha=([\d.]+), "
                  r"wr=([\d.]+), fcomp=([\d.]+), fpm=([\d.]+), "
                  r"kw=([\d.]+), sq=([\d.]+); CENSUS corrected-kernel: "
                  r"mu=\(([\d.]+), ([\d.]+)\)", ln)
    if mm:
        ship[(mm.group(1), int(mm.group(2)))] = dict(
            alpha=float(mm.group(3)), wr=float(mm.group(4)),
            fcomp=float(mm.group(5)), fpm=float(mm.group(6)),
            kw=float(mm.group(7)), sq=float(mm.group(8)),
            mub=float(mm.group(9)), muh=float(mm.group(10)))
assert len(ship) == 4, ship

# --------------------------------------------------------------------------
t0 = time.time()
P("8M THE JOINT-COHERENCE SCAN (pre-reg committed BEFORE any run; "
  "bars locked; NO credence movement any branch)")
P(f"data N per bin = {NDATA}; observed pair = (band {NOBS}, "
  f"cliff {NHI}); best possible CLL = "
  f"{poisson.logpmf(NOBS, NOBS)+poisson.logpmf(NHI, NHI):.2f}")
P("")
g0_ok = g1_ok = g2_ok = g3_ok = True
R = {}
ncache = {}
for seed in SEEDS:
    pop = build_pop(seed)
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        c9 = np.load(f'data/stage7j_cube_full_lker_{seed}_{law}.npy')
        assert c9.shape == (5, 2, 5, 6, 1, 2, 6, 3, 4), c9.shape
        pri = (prior_eta.reshape((1, 2) + (1,)*7)
               + LNPI.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1)))
        cb = c9 + pri
        bx = np.unravel_index(np.nanargmax(cb), cb.shape)
        # ---- G8M-0: direct band_mu at the recomputed MAP vs shipped -----
        al = A_GRID[bx[0]]; eta = E_GRID[bx[1]]; wr = WR_GRID[bx[2]]
        fcm = FCOMP_GRID[bx[3]]; fpm = FPM_GRID[bx[6]]
        kwv = KW_GRID[bx[7]]; sqv = SQ_GRID[bx[8]]
        sh = ship[(law, seed)]
        okp = all(abs(v-sh[k]) < 1e-9 for v, k in
                  ((al, 'alpha'), (wr, 'wr'), (fcm, 'fcomp'),
                   (fpm, 'fpm'), (kwv, 'kw'), (sqv, 'sq')))
        e_f = e_of_x(pop, eta, wr)
        tab_a = 1.0 + al*(TAB-1.0)
        vp_f = vp_c(pop, e_f, tab_a)
        o_f = run(pop['a_s'], e_f, pop['psi0'], pop['f_ip'],
                  pop['M_s'], pop['uph'], 8, 2500, 5, a0=A0_CAN,
                  tab=tab_a, lny0=LNY0, dlny=DLNY, vp=vp_f)
        mu_d, mh_d = band_mu(pop, o_f, fcm, kwv, fpm, sqv)
        ok0 = (okp and abs(mu_d-sh['mub']) <= 0.05
               and abs(mh_d-sh['muh']) <= 0.05)
        g0_ok &= ok0
        P(f"[{law} {seed}] G8M-0 MAP (alpha={al}, eta={eta}, wr={wr}, "
          f"fcomp={fcm}, fpm={fpm}, kw={kwv}, sq={sqv}): direct mu="
          f"({mu_d:.2f}, {mh_d:.2f}) vs shipped ({sh['mub']:.2f}, "
          f"{sh['muh']:.2f}) -> {'PASS' if ok0 else 'FAIL'}")
        # ---- THE SCAN ---------------------------------------------------
        MUB = np.zeros((5, 2, 5, 6, 3, 6, 4))
        MUH = np.zeros_like(MUB)
        RES = np.zeros((5, 2, 5, 6, 3, 6))
        for ai, alv in enumerate(A_GRID):
            tab_i = 1.0 + alv*(TAB-1.0)
            for ei, etav in enumerate(E_GRID):
                for wi, wrv in enumerate(WR_GRID):
                    ck = (seed, ei, wi)
                    if alv == 0.0 and ck in ncache:
                        blk = ncache[ck]
                    else:
                        e_i = e_of_x(pop, etav, wrv)
                        vp_i = vp_c(pop, e_i, tab_i)
                        o_i = run(pop['a_s'], e_i, pop['psi0'],
                                  pop['f_ip'], pop['M_s'], pop['uph'],
                                  8, 2500, 5, a0=A0_CAN, tab=tab_i,
                                  lny0=LNY0, dlny=DLNY, vp=vp_i)
                        blk = scan_block(pop, o_i)
                        if alv == 0.0:
                            ncache[ck] = blk
                    MUB[ai, ei, wi] = blk[0]
                    MUH[ai, ei, wi] = blk[1]
                    RES[ai, ei, wi] = blk[2]
            print(f"  .. {law} {seed} alpha={alv} scanned "
                  f"({(time.time()-t0)/60:.1f} min)", flush=True)
        # to cube axis order (A,E,WR,FCOMP,FPM,KW,SQ)
        MUBc = MUB.transpose(0, 1, 2, 3, 5, 4, 6)
        MUHc = MUH.transpose(0, 1, 2, 3, 5, 4, 6)
        FLR = 0.5*RES.transpose(0, 1, 2, 3, 5, 4)[..., None]
        # ---- G8M-3: scan block vs the direct call at the MAP ------------
        idM = (bx[0], bx[1], bx[2], bx[3], bx[6], bx[7], bx[8])
        ok3 = (abs(MUBc[idM]-mu_d) <= 1e-9
               and abs(MUHc[idM]-mh_d) <= 1e-9)
        g3_ok &= ok3
        P(f"[{law} {seed}] G8M-3 scan-vs-direct at MAP: "
          f"({MUBc[idM]:.6f}, {MUHc[idM]:.6f}) vs ({mu_d:.6f}, "
          f"{mh_d:.6f}) -> {'PASS' if ok3 else 'FAIL'}")
        # ---- G8M-1: bare-Newton row -------------------------------------
        mun = MUBc[0, 1, 2, 0, 1, 0, 0]
        ok1 = mun < 3.0
        g1_ok &= ok1
        P(f"[{law} {seed}] G8M-1 bare-Newton (a=0, eta=1.3, wr=0.3, "
          f"fcomp=0, fpm=1.5, kw=0.7, sq=0): mu_band = {mun:.2f} "
          f"(< 3.0) -> {'PASS' if ok1 else 'FAIL'}")
        # ---- CLL, admissibility, frontier, joint ------------------------
        CLL = (poisson.logpmf(NOBS, np.maximum(MUBc, FLR))
               + poisson.logpmf(NHI, np.maximum(MUHc, FLR)))
        CLLl = (poisson.logpmf(NOBS, np.maximum(MUBc, 1e-12))
                + poisson.logpmf(NHI, np.maximum(MUHc, 1e-12)))
        JP, JPl = np.exp(CLL), np.exp(CLLl)
        cbr = np.nanmax(cb, axis=(4, 5))
        dkin = float(np.nanmax(cb)) - cbr
        supp = dkin < 1e6          # the operative prior support
        admc, admlc = (JP >= 1e-3), (JPl >= 1e-3)
        adm, adml = admc & supp, admlc & supp
        Nadm, Nadml = int(adm.sum()), int(adml.sum())
        Ncen = int(admc.sum())
        P(f"[{law} {seed}] floor stats: median {np.median(FLR):.3f}, "
          f"max {np.max(FLR):.3f} (half-count resolution)")
        price = None
        if Nadm:
            price = float(np.nanmin(dkin[adm]))
            a_adm = sorted(set(A_GRID[i] for i in np.where(adm)[0]))
            fl = np.where(adm.ravel())[0]
            best = fl[int(np.nanargmin(dkin.ravel()[fl]))]
            bc = np.unravel_index(best, adm.shape)
            cell = ", ".join(f"{n}={g[i]:g}" for n, g, i
                             in zip(AXN, AXG, bc))
            P(f"[{law} {seed}] ADMISSIBLE SET (prior-supported): "
              f"N={Nadm} of 21600 (legacy-floor N={Nadml}; "
              f"prior-blind census-only N={Ncen}); alphas={a_adm}; "
              f"cheapest dkin={price:.1f} at ({cell}) mu="
              f"({MUBc[bc]:.2f}, {MUHc[bc]:.2f}) jointP={JP[bc]:.2e}")
        else:
            P(f"[{law} {seed}] ADMISSIBLE SET (prior-supported): N=0 "
              f"of 21600 (legacy-floor N={Nadml}; prior-blind "
              f"census-only N={Ncen}) - no supported cell reaches "
              f"the 8H bar; max supported jointP = "
              f"{np.nanmax(np.where(supp, JP, np.nan)):.2e}")
        fr = []
        for X in FRX:
            msk = dkin <= X
            if msk.any():
                cf = np.where(msk, CLL, -np.inf)
                ic = np.unravel_index(int(np.argmax(cf)), CLL.shape)
                fr.append(f"{X}:{CLL[ic]:.1f}@a={A_GRID[ic[0]]:g}")
            else:
                fr.append(f"{X}:-")
        P(f"[{law} {seed}] FRONTIER dkin<=X -> best CLL@alpha: "
          + "  ".join(fr))
        cbJ = cb + CLL.reshape(5, 2, 5, 6, 1, 1, 6, 3, 4)
        assert cbJ.shape == cb.shape and np.all(np.isfinite(CLL))
        bxJ = np.unravel_index(np.nanargmax(cbJ), cbJ.shape)
        idJ = (bxJ[0], bxJ[1], bxJ[2], bxJ[3], bxJ[6], bxJ[7], bxJ[8])
        ok2 = abs((cbJ[bxJ]-cb[bxJ]) - CLL[idJ]) <= 1e-9
        g2_ok &= ok2
        a_k, dn_k, _ = read9(cb)
        aJ, dnJ, poJ = read9(cbJ, axes_post=(3, 6, 7, 8))
        aJl, dnJl, _ = read9(cb + CLLl.reshape(5, 2, 5, 6, 1, 1,
                                               6, 3, 4))
        emax = [AXN[i] for i, (j, g_) in
                enumerate(zip(idJ, AXG)) if j == len(g_)-1]
        emin = [AXN[i] for i, j in enumerate(idJ) if j == 0]
        cellJ = ", ".join(f"{n}={g[i]:g}" for n, g, i
                          in zip(AXN, AXG, idJ))
        P(f"[{law} {seed}] JOINT MAP: ({cellJ}); mu=({MUBc[idJ]:.2f}, "
          f"{MUHc[idJ]:.2f}) jointP={JP[idJ]:.2e}; dkin="
          f"{dkin[idJ]:.1f}; CLL={CLL[idJ]:.2f}; G8M-2 "
          f"{'PASS' if ok2 else 'FAIL'}; edges max={emax or '-'} "
          f"min={emin or '-'}")
        P(f"[{law} {seed}] JOINT-DIAGNOSTIC: aJ={aJ:.2f} dNJ={dnJ:+.1f} "
          f"(kin read a={a_k:.2f} dN={dn_k:+.1f}); legacy-floor "
          f"co-read aJ={aJl:.2f} dNJ={dnJl:+.1f}")
        P(f"[{law} {seed}]   P(fcomp)={np.round(poJ[3], 2).tolist()} "
          f"P(fpm)={np.round(poJ[6], 2).tolist()}")
        P(f"[{law} {seed}]   P(kw)={np.round(poJ[7], 2).tolist()} "
          f"P(sq)={np.round(poJ[8], 2).tolist()}")
        P("")
        R[(law, seed)] = dict(aJ=aJ, dnJ=dnJ, aJl=aJl, dnJl=dnJl,
                              a_k=a_k, dn_k=dn_k, Nadm=Nadm,
                              Nadml=Nadml, price=price)

# ---- bars ----------------------------------------------------------------
if not (g0_ok and g1_ok and g2_ok and g3_ok):
    P("GATES FAILED (G8M-0 " + ('PASS' if g0_ok else 'FAIL')
      + ", G8M-1 " + ('PASS' if g1_ok else 'FAIL')
      + ", G8M-2 " + ('PASS' if g2_ok else 'FAIL')
      + ", G8M-3 " + ('PASS' if g3_ok else 'FAIL')
      + ") - DO NOT QUOTE; verdict WITHHELD")
else:
    P("GATES: G8M-0 4/4, G8M-1 4/4, G8M-2 4/4, G8M-3 4/4 - ALL PASS")
    P("")
    ncontain = sum(1 for r in R.values() if r['Nadm'] >= 1)
    if ncontain >= 3:
        b1 = 'CLASS-CONTAINS'
    elif ncontain <= 1:
        b1 = 'CLASS-EMPTY'
    else:
        b1 = 'CLASS-SPLIT'
    parts = [f"{b1} ({ncontain}/4 law-seeds admissible)"]
    if b1 == 'CLASS-CONTAINS':
        for law in ('simple', 'BE'):
            ps = [R[(law, s)]['price'] for s in SEEDS
                  if R[(law, s)]['price'] is not None]
            pm = float(np.mean(ps)) if ps else float('nan')
            grade = ('CHEAP' if pm <= 10 else
                     'PRICED' if pm <= 40 else 'SEVERE')
            parts.append(f"{law} admission price {pm:.1f} lnL "
                         f"({len(ps)}/2 seeds) -> {grade}")
    elif b1 == 'CLASS-EMPTY':
        parts.append("no cell of the operative class reproduces "
                     "(9,2) at the 8H bar - the census inconsistency "
                     "is MODEL-CLASS-level at cube grade; the wobble "
                     "tail-SHAPE successor is MANDATORY before the "
                     "alpha machinery is called final")
    for law in ('simple', 'BE'):
        am = float(np.mean([R[(law, s)]['aJ'] for s in SEEDS]))
        dm = float(np.mean([R[(law, s)]['dnJ'] for s in SEEDS]))
        parts.append(f"JOINT-DIAGNOSTIC {law}: alpha_J = {am:.2f}, "
                     f"dN_J = {dm:+.1f} (seed means; CO-QUOTED-"
                     f"DIAGNOSTIC only)")
    viol = sum(1 for r in R.values() if r['dnJ'] < r['dn_k'] - 1e-9)
    parts.append("B4 Newton flank: "
                 + (f"AGAINST-EXPECTATION ({viol}/4 law-seeds have "
                    f"dN_J < dN_kin - the rescuing corner is in the "
                    f"JOINT MAP rows)" if viol >= 3 else
                    f"expectation held ({4-viol}/4 law-seeds "
                    f"dN_J >= dN_kin)"))
    P("==> 8M VERDICT (locked bars): " + "; ".join(parts))
    P("    NO credence movement (pre-stated; external-only per "
      "8K-b/round 6 - internal-coherence accounting only).")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage8m_jointcensus.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage8m_jointcensus.txt")
