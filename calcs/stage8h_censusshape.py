"""STAGE 8H — THE CENSUS SHAPE-AND-ATTRIBUTION CONTEST (pre-reg
652cb4c, committed BEFORE any run; the width-object program leg 2).

Forward ladder A/B/C/K2a/K2b/K2n at the FREED-e MAP cells (the esec
cubes), scored by joint Poisson LL at the observed census pair
(band=9, cliff=2); G8H-0 regresses the machinery against the shipped
7K-b rows at the landed cells (parsed at runtime).  Bars locked in
the pre-reg; expected outcome pre-stated ALL-FAIL with companion
attribution; NO credence move any branch.
Output: data/stage8h_censusshape.txt
"""
import re
import numpy as np, time
from astropy.io import fits
from scipy.stats import poisson, norm

OUT = 'data/stage8h_censusshape.txt'
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

# --- catalog (7K-b verbatim) ----------------------------------------------
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
P(f"STAGE 8H: data N per bin = {NDATA}; observed pair = "
  f"(band {NOBS}, cliff {NHI})")

MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])

# --- grids + priors -------------------------------------------------------
N = 500_000
A_GRID  = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID  = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
EIN_GRID = np.array([0.5, 1.0, 1.5, 2.0])
ERF_GRID = np.array([0.80, 0.90, 0.95])
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FPM_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW_GRID = np.array([0.7, 1.0, 1.4])
SQ_GRID = np.array([0.0, 0.1, 0.2, 0.3])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
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

def e_of_x(p, eta, wr, ein, erf):
    if ein == 1.0 and erf == 0.90:
        return e_of(p, eta, wr)
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

def band_mu(p, o, fcm, kw, fpm, sq, ceil=None):
    """7K-b band_mu + the AMENDMENT-1 smear bound: with ceil set, the
    per-system log-normal draw is TRUNCATED at g_max = ln(max(vtn,
    ceil)/vtn)/sq and RENORMALIZED via the uniform map (conditioned-
    on-bound; no pile-up).  The first-run CLAMP form manufactured
    band counts and was ceiling-degenerate — caught by G8H-1, logged
    pre-quote, preserved in ..._run1.txt."""
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
        g_ = p['gs'][idx][keep]
        if ceil is not None and sq > 0.0:
            cap = np.maximum(vtn, ceil)
            gmax = np.log(np.maximum(cap/vtn, 1e-300))/sq
            u = norm.cdf(g_)
            g_ = norm.ppf(np.clip(u*norm.cdf(gmax), 1e-16, 1-1e-16))
        vts = vtn*np.exp(sq*g_)
        nk = max(int(keep.sum()), 1)
        inb = float(np.sum((gmn >= 75) & (vts >= 1.414) & (vts < 1.67)))
        inh = float(np.sum((gmn >= 75) & (vts >= 1.67) & (vts < 2.2)))
        mu += NDATA[bi]*inb/nk
        mu_hi += NDATA[bi]*inh/nk
    return mu, mu_hi

def jointp(mu_b, mu_c):
    return float(poisson.pmf(NOBS, max(mu_b, 1e-12))
                 * poisson.pmf(NHI, max(mu_c, 1e-12)))

def jll(mu_b, mu_c):
    return float(poisson.logpmf(NOBS, max(mu_b, 1e-12))
                 + poisson.logpmf(NHI, max(mu_c, 1e-12)))

# --- G8H-0 reference rows (shipped 7K-b output, parsed at runtime) --------
kb = open('data/stage7kb_census.txt').read()
ship_boost, ship_sq0 = {}, {}
seed_ctx = None
for ln in kb.splitlines():
    ms = re.match(r"== seed (\d+) ==", ln)
    if ms:
        seed_ctx = int(ms.group(1))
    mb = re.match(r"boost (simple|BE) \(alpha=([\d.]+)\): mu_band = "
                  r"([\d.]+) vs obs 9; overshoot \[1.67,2.2\) mu = "
                  r"([\d.]+)", ln)
    if mb:
        ship_boost[(mb.group(1), seed_ctx)] = (
            float(mb.group(2)), float(mb.group(3)), float(mb.group(4)))
    mp = re.match(r"POST-HOC boost (simple|BE) sq=0 fpm=1.5: mu_band = "
                  r"([\d.]+) \(obs 9\), overshoot mu = ([\d.]+)", ln)
    if mp:
        ship_sq0[(mp.group(1), seed_ctx)] = (float(mp.group(2)),
                                             float(mp.group(3)))
assert len(ship_boost) == 4 and len(ship_sq0) == 4, (ship_boost, ship_sq0)

t0 = time.time()
g0_ok, g1_ok = True, True
rows = {}
for seed in (31, 101):
    P(f"\n== seed {seed} ==")
    p = build_pop(seed)
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        # ---- G8H-0: landed-cell regression vs the shipped rows ----------
        cw = np.load(f'data/stage7j_cube_full_photow3_{seed}_{law}.npy')
        cb = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1)) \
             + LNPI.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
        bxa = np.unravel_index(np.nanargmax(cb), cb.shape)
        al_l = A_GRID[bxa[0]]
        e_l = e_of(p, E_GRID[bxa[1]], WR_GRID[bxa[2]])
        tab_l = 1.0 + al_l*(TAB-1.0)
        vp_l = vp_c(p, e_l, tab_l)
        o_l = run(p['a_s'], e_l, p['psi0'], p['f_ip'], p['M_s'],
                  p['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_l,
                  lny0=LNY0, dlny=DLNY, vp=vp_l)
        sa, sb, sc = ship_boost[(law, seed)]
        assert abs(al_l - sa) < 1e-9, (al_l, sa)
        mu_l, mh_l = band_mu(p, o_l, FCOMP_GRID[bxa[3]], KW_GRID[bxa[7]],
                             FPM_GRID[bxa[6]], SQ_GRID[bxa[8]])
        mu_0, mh_0 = band_mu(p, o_l, FCOMP_GRID[bxa[3]], KW_GRID[bxa[7]],
                             1.5, 0.0)
        za, zb = ship_sq0[(law, seed)]
        ok0 = (abs(mu_l-sb) <= 0.05 and abs(mh_l-sc) <= 0.05
               and abs(mu_0-za) <= 0.05 and abs(mh_0-zb) <= 0.05)
        g0_ok &= ok0
        P(f"G8H-0 {law} seed {seed}: landed {mu_l:.2f}/{mh_l:.2f} vs "
          f"shipped {sb:.2f}/{sc:.2f}; sq0 {mu_0:.2f}/{mh_0:.2f} vs "
          f"{za:.2f}/{zb:.2f} -> {'PASS' if ok0 else 'FAIL'}")
        # ---- the freed-e MAP cell ---------------------------------------
        c11 = np.load(f'data/stage7j_cube_full_esec_{seed}_{law}.npy')
        cb11 = (c11
                + prior_eta.reshape((1, 2) + (1,)*9)
                + LNPI.reshape((1, 1, 1, 1, 1, 6, 1, 1, 1, 1, 1)))
        bx = np.unravel_index(np.nanargmax(cb11), cb11.shape)
        al = A_GRID[bx[0]]; eta = E_GRID[bx[1]]; wr = WR_GRID[bx[2]]
        ein = EIN_GRID[bx[3]]; erf = ERF_GRID[bx[4]]
        fcm = FCOMP_GRID[bx[5]]; fpm = FPM_GRID[bx[8]]
        kw = KW_GRID[bx[9]]; sq = SQ_GRID[bx[10]]
        P(f"FREED CELL {law} seed {seed}: alpha={al}, eta={eta}, "
          f"wr={wr}, ein={ein}, erf={erf}, fcomp={fcm}, fpm={fpm}, "
          f"kw={kw}, sq={sq}")
        e_f = e_of_x(p, eta, wr, ein, erf)
        tab_a = 1.0 + al*(TAB-1.0)
        vp_f = vp_c(p, e_f, tab_a)
        o_f = run(p['a_s'], e_f, p['psi0'], p['f_ip'], p['M_s'],
                  p['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_a,
                  lny0=LNY0, dlny=DLNY, vp=vp_f)
        o_n = run(p['a_s'], e_f, p['psi0'], p['f_ip'], p['M_s'],
                  p['uph'], 8, 2500, 1)
        C_al = float(np.sqrt(2*(1+al*0.36)))
        ladder = [
            ('A  fcm=0 sq=0',        o_f, 0.0, kw, fpm, 0.0, None),
            ('B  companions sq=0',   o_f, fcm, kw, fpm, 0.0, None),
            ('C  +smear (K1)',       o_f, fcm, kw, fpm, sq,  None),
            (f'K2a bound C({al})={C_al:.3f}', o_f, fcm, kw, fpm, sq, C_al),
            ('K2b bound C(1)=1.649', o_f, fcm, kw, fpm, sq, 1.649),
            ('K2n bound C(0)=1.414', o_f, fcm, kw, fpm, sq, 1.414),
        ]
        res = {}
        for tag, o_, f_, k_, pm_, sq_, cl_ in ladder:
            mu, mh = band_mu(p, o_, f_, k_, pm_, sq_, ceil=cl_)
            res[tag.split()[0]] = (tag, mu, mh)
            P(f"  {tag:<26} mu_band = {mu:6.2f}, mu_hi = {mh:6.2f}, "
              f"jointP = {jointp(mu, mh):.2e}, LL = {jll(mu, mh):7.2f}")
        # Newton consistency note (not a gate)
        mu_n, mh_n = band_mu(p, o_n, 0.0, kw, fpm, 0.0)
        P(f"  [note] Newton A at freed e: mu_band = {mu_n:.2f}, "
          f"mu_hi = {mh_n:.2f} (8F-b analytic leakage scale ~<= 2.7)")
        # G8H-1 (amendment-1 invariants for the truncation kernel)
        okc = (res['K2a'][2] <= res['C'][2] + 1e-6
               and res['K2b'][2] <= res['C'][2] + 1e-6
               and res['K2n'][2] <= res['C'][2] + 1e-6
               and res['K2n'][2] <= res['B'][2] + 1e-6
               and (abs(res['K2a'][1]-res['K2n'][1])
                    + abs(res['K2b'][1]-res['K2n'][1])) > 1e-9)
        g1_ok &= okc
        P(f"  G8H-1 (cliff-only-removes + K2n<=B cliff + "
          f"non-degeneracy): {'PASS' if okc else 'FAIL'}")
        rows[(law, seed)] = res

P("")
if not (g0_ok and g1_ok):
    P("GATES FAILED (G8H-0 " + ('PASS' if g0_ok else 'FAIL')
      + ", G8H-1 " + ('PASS' if g1_ok else 'FAIL')
      + ") - DO NOT QUOTE; verdict WITHHELD")
else:
    P("GATES: G8H-0 4/4 PASS, G8H-1 4/4 PASS")
    P("")
    # ---- bars -------------------------------------------------------
    keys = ['A', 'B', 'C', 'K2a', 'K2b', 'K2n']
    adm = {k: 0 for k in keys}
    for (law, seed), res in rows.items():
        for k in keys:
            if jointp(res[k][1], res[k][2]) >= 1e-3:
                adm[k] += 1
    P("ADMISSIBILITY (jointP >= 1e-3, of 4 law-seed runs): "
      + ", ".join(f"{k}={adm[k]}" for k in keys))
    winners = [k for k in keys if adm[k] >= 3]
    # attribution: companion share of the overshoot at C
    attr = {}
    for law in ('simple', 'BE'):
        sh = np.mean([(rows[(law, s)]['B'][2]-rows[(law, s)]['A'][2])
                      / max(rows[(law, s)]['C'][2], 1e-9)
                      for s in (31, 101)])
        attr[law] = float(sh)
    P(f"ATTRIBUTION (companion share of mu_hi at C, seed means): "
      f"simple {attr['simple']:.2f}, BE {attr['BE']:.2f}")
    if not winners:
        comp = all(attr[l] >= 0.5 for l in ('simple', 'BE'))
        v = ("ALL-FAIL: no ladder config reproduces (9,2) at the "
             "1e-3 admissibility bar in >= 3/4 runs - the census "
             "reproduction problem STANDS at kernel grade"
             + ("; ATTRIBUTION >= 0.5 both laws -> the census pair "
                "is named a COMPANION-SECTOR diagnostic (feeds the "
                "7J-z4 wobble-binding line); the width-object "
                "reading of the census is DEMOTED" if comp else
                "; attribution below the naming bar - flood shared "
                "across channels, reported as-is"))
    elif len(winners) == 1:
        k = winners[0]
        if k in ('K2a', 'K2b', 'K2n'):
            v = (f"SHAPE-SELECTED: {k} uniquely admissible - the "
                 f"width object has BOUNDED support at the "
                 f"{rows[('simple', 31)][k][0]} ceiling (the ceiling "
                 f"position doubles as an alpha read, mechanical)")
        else:
            v = (f"FREED-MODEL-RESOLVED at {k}: the 8G revision "
                 f"itself repairs the census; the smear does not "
                 f"live in this channel")
    else:
        v = (f"RESOLVED-DEGENERATE: admissible set {winners} - "
             f"reported as-is, no selection sentence")
    P("")
    P(f"==> 8H VERDICT (locked bars): {v}")
    P("    NO credence move (pre-stated in every branch).")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")
