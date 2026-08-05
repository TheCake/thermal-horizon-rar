"""STAGE 8N — THE FLOOD ANATOMY (pre-reg committed BEFORE any run).
Part 0 regenerates the design-time kill of the residual-survival
candidate (the flood carriers are RUWE-silent: L->1, R->0 at long P).
Parts a/b/c measure the census flood's anatomy at the operative
(8L-b corrected-kernel) MAP cells: the P-locus of the band/cliff
carriers, the wobble-vs-mass channel attribution, and the paired
flat-q vs twin-heavy (t5, the v2c GV7-measured law) census dial.
NO credence movement in any branch (measurement round).
Output: data/stage8n_floodanatomy.txt
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
PBINS = [3.0, 10.0, 30.0, 100.0]
PLBL = ['P<3', '3-10', '10-30', '30-100', 'P>=100', 'nocomp']

pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
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

def lk_S(P_yr):
    u_ = np.pi*2.83/np.maximum(P_yr, 1e-9)
    return np.where(u_ < 1e-2, 1.0 - u_*u_/10.0,
                    3.0*(np.sin(u_) - u_*np.cos(u_))
                    / np.maximum(u_, 1e-300)**3)

def R2f(u):
    u = np.asarray(u, dtype=float)
    small = u < 0.05
    s = np.sinc(u/np.pi)
    core = 0.5*(1.0 - s*s
                - 3.0*(np.sin(u) - u*np.cos(u))**2
                / np.maximum(u, 1e-300)**4)
    return np.where(small, u**4/90.0, np.maximum(core, 0.0))

def part0():
    rng = np.random.default_rng(31)
    M_s = 0.6+1.8*rng.random(N)
    M_h = 0.5*M_s
    a_s = ((0.15**-0.6+rng.random(N)*(60.0**-0.6-0.15**-0.6))
           ** (1/-0.6))*1e3
    KW, FCM = 0.7, 0.10
    r_tot = np.zeros(N); w_tot = np.zeros(N); nact = np.zeros(N)
    for k in (1, 2):
        q = 0.1+0.9*rng.random(N)
        logP = rng.normal(5.03, 2.28, N)
        P_yr = 10**logP/365.25
        a_in = (M_h*(1+q)*P_yr**2)**(1/3)
        valid = (a_in < 130.0) & (a_in < a_s/5.0)
        v_orb = 29.78*np.sqrt(M_h*(1+q)/np.maximum(a_in, 1e-3))
        MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        MGs = np.interp(-np.clip(q*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        l_ = 10**(-0.4*(MGs-MGp))
        wfac = np.abs(q/(1+q) - l_/(1+l_))
        u_ = np.pi*2.83/np.maximum(P_yr, 1e-9)
        Lk = np.where(u_ < 1e-2, 1.0-u_*u_/10.0,
                      3.0*(np.sin(u_)-u_*np.cos(u_))
                      / np.maximum(u_, 1e-300)**3)
        act = rng.random(N) < FCM
        r_tot += act*wfac*v_orb*np.sqrt(R2f(u_))/4.74047*valid
        w_tot += act*wfac*v_orb*Lk/4.74047*valid
        nact += act*valid
    sg0 = sig_ok[rng.integers(0, len(sig_ok), N)]/4.74047
    ratio = KW*r_tot/np.maximum(sg0, 1e-12)
    has = nact > 0
    spike = has & (KW*w_tot*4.74047 > 0.3)
    rs = ratio[spike]
    P("PART 0 (the SRVR design-time kill, rng 31): spike carriers "
      f"(kw*w_leak > 0.3 km/s) = {spike.mean():.4f} of systems; "
      f"their residual-to-noise ratio q50 = {np.percentile(rs, 50):.2f}, "
      f"q90 = {np.percentile(rs, 90):.1f} -> the flood is RUWE-SILENT "
      f"(L->1, R->0 at long P); a survival cut at wsr=100 removes "
      f"{np.mean(ratio[spike] > 100):.3f} of the flood carriers. "
      f"SRVR NOT BUILT.")

def build_pop(seed, qlaw='flat'):
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
        u_q = rng.random(N)
        if qlaw == 't5':
            # the 7J-z7 stream-preserving twin-t5 map (GV7 winner)
            s_ = 0.9/1.4
            q = np.where(u_q < s_, 0.1 + 0.8*(u_q/s_),
                         0.9 + 0.1*((u_q-s_)/(1.0-s_)))
        else:
            q = 0.1+0.9*u_q
        logP = rng.normal(5.03, 2.28, N)
        P_yr = 10**logP/365.25
        a_in = (M_h*(1+q)*P_yr**2)**(1/3)
        valid = (a_in < 130.0) & (a_in < p['a_s']/5.0)
        v_orb = 29.78*np.sqrt(M_h*(1+q)/np.maximum(a_in,1e-3))
        S = lk_S(P_yr)
        MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        MGs = np.interp(-np.clip(q*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
        l_ = 10**(-0.4*(MGs-MGp))
        wfac = np.abs(q/(1+q) - l_/(1+l_))
        w = wfac*v_orb*S/4.74047*valid
        wd = rng.normal(size=(N,3)); wd /= np.linalg.norm(wd,axis=1,keepdims=True)
        p['comp'][k] = dict(w=w, wd=wd, uc=rng.random(N),
                            mh=q*M_h*valid, P=P_yr)
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

def band_mu_anat(p, o, fcm, kw, fpm, sq, chan='full', locus=False):
    """8L-b band_mu with a channel switch and an optional P-locus
    decomposition; chan='full' is the verbatim expression sequence."""
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
    lb = np.zeros(len(PLBL)); lc = np.zeros(len(PLBL))
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500:
            continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = POOLS[bi][p['pick'][bi][idx] % len(POOLS[bi])]/4.74047
        cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
        mh_tot = np.zeros(len(idx))
        w1 = np.zeros(len(idx)); w2 = np.zeros(len(idx))
        for k in (1, 2):
            c = p['comp'][k]
            act = (c['uc'][idx] < fcm) if chan != 'none' \
                else np.zeros(len(idx), dtype=bool)
            if chan != 'nowob':
                cvp += act*c['w'][idx]*c['wd'][idx,0]
                cvq += act*c['w'][idx]*c['wd'][idx,1]
            if chan != 'nomass':
                mh_tot += act*c['mh'][idx]
            aw = act*np.abs(c['w'][idx])
            if k == 1:
                w1 = aw
            else:
                w2 = aw
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
        mb = (gmn >= 75) & (vts >= 1.414) & (vts < 1.67)
        mc = (gmn >= 75) & (vts >= 1.67) & (vts < 2.2)
        mu += NDATA[bi]*float(np.sum(mb))/nk
        mu_hi += NDATA[bi]*float(np.sum(mc))/nk
        if locus:
            hasact = ((w1 > 0) | (w2 > 0))[keep]
            Pdom = np.where(w1 >= w2, p['comp'][1]['P'][idx],
                            p['comp'][2]['P'][idx])[keep]
            cls = np.digitize(Pdom, PBINS)      # 0..4
            cls = np.where(hasact, cls, len(PLBL)-1)
            for ci in range(len(PLBL)):
                sel = cls == ci
                lb[ci] += NDATA[bi]*float(np.sum(mb & sel))/nk
                lc[ci] += NDATA[bi]*float(np.sum(mc & sel))/nk
    if locus:
        return mu, mu_hi, lb, lc
    return mu, mu_hi

def jointp(mu_b, mu_c):
    return float(poisson.pmf(NOBS, max(mu_b, 1e-12))
                 * poisson.pmf(NHI, max(mu_c, 1e-12)))

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

t0 = time.time()
P("8N THE FLOOD ANATOMY (pre-reg committed BEFORE any run; "
  "measurement round; NO credence movement any branch)")
P(f"data N per bin = {NDATA}; observed pair = (band {NOBS}, "
  f"cliff {NHI})")
P("")
part0()
P("")
g0_ok = g1_ok = g2_ok = g3_ok = True
r1_fr, r3_rows = [], []
for seed in SEEDS:
    pf = build_pop(seed, 'flat')
    pt = build_pop(seed, 't5')
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        c9 = np.load(f'data/stage7j_cube_full_lker_{seed}_{law}.npy')
        pri = (prior_eta.reshape((1, 2) + (1,)*7)
               + LNPI.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1)))
        bx = np.unravel_index(np.nanargmax(c9 + pri), c9.shape)
        al = A_GRID[bx[0]]; eta = E_GRID[bx[1]]; wr = WR_GRID[bx[2]]
        sh = ship[(law, seed)]
        assert all(abs(v-sh[k]) < 1e-9 for v, k in
                   ((al, 'alpha'), (wr, 'wr'))), (al, wr, sh)
        fcm, fpm = sh['fcomp'], sh['fpm']
        kwv, sqv = sh['kw'], sh['sq']
        e_f = e_of_x(pf, eta, wr)
        tab_a = 1.0 + al*(TAB-1.0)
        vp_f = vp_c(pf, e_f, tab_a)
        o_f = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'], pf['M_s'],
                  pf['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_a,
                  lny0=LNY0, dlny=DLNY, vp=vp_f)
        # t5 pop: identical non-companion draws -> same orbits
        o_t = o_f
        # ---- G8N-0 + channels + locus (flat) ------------------------
        mu_f, mh_f, lb, lc = band_mu_anat(pf, o_f, fcm, kwv, fpm,
                                          sqv, 'full', locus=True)
        ok0 = abs(mu_f-sh['mub']) <= 0.05 and abs(mh_f-sh['muh']) <= 0.05
        g0_ok &= ok0
        P(f"[{law} {seed}] G8N-0 MAP (alpha={al}, eta={eta}, wr={wr}, "
          f"fcomp={fcm}, fpm={fpm}, kw={kwv}, sq={sqv}): mu="
          f"({mu_f:.2f}, {mh_f:.2f}) vs shipped ({sh['mub']:.2f}, "
          f"{sh['muh']:.2f}) -> {'PASS' if ok0 else 'FAIL'}")
        ok3 = (abs(lb.sum()-mu_f) <= 1e-9 and abs(lc.sum()-mh_f) <= 1e-9)
        g3_ok &= ok3
        P(f"[{law} {seed}] G8N-3 locus completeness: band "
          f"{lb.sum():.6f}={mu_f:.6f}, cliff {lc.sum():.6f}={mh_f:.6f}"
          f" -> {'PASS' if ok3 else 'FAIL'}")
        P(f"[{law} {seed}] P-LOCUS band : "
          + "  ".join(f"{n}={v:.2f}" for n, v in zip(PLBL, lb)))
        P(f"[{law} {seed}] P-LOCUS cliff: "
          + "  ".join(f"{n}={v:.2f}" for n, v in zip(PLBL, lc)))
        carr_b = lb[:5].sum()
        fr10 = float(lb[2:5].sum()/max(carr_b, 1e-12))
        r1_fr.append(fr10)
        P(f"[{law} {seed}] R1: band companion-carriers = "
          f"{carr_b:.2f} of {mu_f:.2f}; fraction at P >= 10 yr = "
          f"{fr10:.2f} (expectation >= 0.60)")
        for tag, ch in (('nowob (mass only)', 'nowob'),
                        ('nomass (wobble only)', 'nomass'),
                        ('none (both off)', 'none')):
            mu_c, mh_c = band_mu_anat(pf, o_f, fcm, kwv, fpm, sqv, ch)
            P(f"[{law} {seed}] CHANNEL {tag:<22}: mu=({mu_c:.2f}, "
              f"{mh_c:.2f})")
            if ch == 'none':
                mu_0, mh_0 = band_mu_anat(pf, o_f, 0.0, kwv, fpm,
                                          sqv, 'full')
                ok2 = (abs(mu_c-mu_0) <= 1e-9 and abs(mh_c-mh_0) <= 1e-9)
                g2_ok &= ok2
                P(f"[{law} {seed}] G8N-2 none==fcomp0: "
                  f"{'PASS' if ok2 else 'FAIL'}")
        # ---- G8N-1: t5 stream preservation at fcomp=0 ---------------
        a0f = band_mu_anat(pf, o_f, 0.0, kwv, fpm, sqv, 'full')
        a0t = band_mu_anat(pt, o_t, 0.0, kwv, fpm, sqv, 'full')
        ok1 = (abs(a0f[0]-a0t[0]) <= 1e-9 and abs(a0f[1]-a0t[1]) <= 1e-9)
        g1_ok &= ok1
        P(f"[{law} {seed}] G8N-1 t5==flat at fcomp=0: "
          f"({a0t[0]:.6f}, {a0t[1]:.6f}) vs ({a0f[0]:.6f}, "
          f"{a0f[1]:.6f}) -> {'PASS' if ok1 else 'FAIL'}")
        # ---- (c) the q-law dial -------------------------------------
        for fv in (fcm, 0.20, 0.35):
            muF, mhF = band_mu_anat(pf, o_f, fv, kwv, fpm, sqv, 'full')
            muT, mhT = band_mu_anat(pt, o_t, fv, kwv, fpm, sqv, 'full')
            P(f"[{law} {seed}] QDIAL fcomp={fv:.2f}: flat mu="
              f"({muF:.2f}, {mhF:.2f}) jP={jointp(muF, mhF):.1e} | "
              f"t5 mu=({muT:.2f}, {mhT:.2f}) jP={jointp(muT, mhT):.1e}"
              f" | d_cliff={100*(mhT-mhF)/max(mhF,1e-9):+.0f}%, "
              f"d_band={100*(muT-muF)/max(muF,1e-9):+.0f}%")
            if abs(fv-fcm) < 1e-9:
                r3_rows.append(((mhT-mhF)/max(mhF, 1e-9), muT-muF))
        P("")

if not (g0_ok and g1_ok and g2_ok and g3_ok):
    P("GATES FAILED (G8N-0 " + ('PASS' if g0_ok else 'FAIL')
      + ", G8N-1 " + ('PASS' if g1_ok else 'FAIL')
      + ", G8N-2 " + ('PASS' if g2_ok else 'FAIL')
      + ", G8N-3 " + ('PASS' if g3_ok else 'FAIL')
      + ") - DO NOT QUOTE; verdict WITHHELD")
else:
    P("GATES: G8N-0 4/4, G8N-1 4/4, G8N-2 4/4, G8N-3 4/4 - ALL PASS")
    P("")
    r1m = float(np.mean(r1_fr))
    dclm = float(np.mean([t[0] for t in r3_rows]))
    dbdm = float(np.mean([t[1] for t in r3_rows]))
    r1v = ('CONFIRMED' if r1m >= 0.60 else 'AGAINST-EXPECTATION')
    r3v = ('REDUCES-MATERIALLY' if (dclm <= -0.30 and dbdm <= 0.0)
           else 'NEUTRAL-OR-WORSENS')
    P(f"==> 8N READINGS (locked grammar): R1 P-locus: mean band-"
      f"carrier fraction at P >= 10 yr = {r1m:.2f} -> {r1v}; "
      f"R3 q-law dial at MAP fcomp: mean d_cliff = {100*dclm:+.0f}%, "
      f"mean d_band = {dbdm:+.2f} -> {r3v}")
    if r3v == 'REDUCES-MATERIALLY':
        P("    SUCCESSOR NAMED: the combined t5+lker operative mode "
          "(GPU round, own pre-reg).")
    else:
        P("    SUCCESSOR NAMED: the population-prior round (the "
          "subsystem-period literature scout) + the external T2/8O "
          "arbiter (RV channel sees the long-P flood carriers).")
    P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage8n_floodanatomy.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage8n_floodanatomy.txt")
