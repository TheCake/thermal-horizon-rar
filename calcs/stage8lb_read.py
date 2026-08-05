"""STAGE 8L-b reader — THE KERNEL ROUND (pre-reg 2eeb161, committed
BEFORE any run).  The derived leakage kernel's corrected reads vs the
esec-slice baseline: the SIGN of d_alpha (the headline), the honest
Newton band, posteriors, Dwob', and the census retest at the
corrected MAP cell.  NO credence movement in any branch.
Output: data/stage8lb_read.txt
"""
import csv
import os
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
SHIP_OFF = {('simple', 31): (0.57, 7.8), ('simple', 101): (0.54, 19.5),
            ('BE', 31): (0.55, 3.4), ('BE', 101): (0.58, 21.9)}

P("8L-b THE KERNEL ROUND READ (pre-reg 2eeb161; bars locked before "
  "any run; NO credence movement any branch)")
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

def dwob(cb):
    cks = []
    for ki in (0, 2):
        sl = cb[:, :, :, :, :, :, :, ki:ki+1]
        cks.append(float(np.nanmax(sl) - np.nanmax(sl[:, :, :, 3:])))
    return cks[1]-cks[0]

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

def jointp(mu_b, mu_c):
    return float(poisson.pmf(NOBS, max(mu_b, 1e-12))
                 * poisson.pmf(NHI, max(mu_c, 1e-12)))

# ---------------- reads ---------------------------------------------------
t0 = time.time()
res = {}
g3_ok, pend = True, []
for seed in SEEDS:
    pop = None
    for law, TAB in (('simple', TAB_S), ('BE', TAB_B)):
        cp = f'data/stage7j_cube_full_lker_{seed}_{law}.npy'
        ep = f'data/stage7j_cube_full_esec_{seed}_{law}.npy'
        if not os.path.exists(cp):
            pend.append(f'{law}@{seed}')
            continue
        c9 = np.load(cp)
        assert c9.shape == (5, 2, 5, 6, 1, 2, 6, 3, 4), c9.shape
        pri = (prior_eta.reshape((1, 2) + (1,)*7)
               + LNPI.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1)))
        cb = c9 + pri
        eb = np.load(ep)[:, :, :, 1, 2] + pri
        a_b, dn_b, _ = read9(eb)
        sa, sd = SHIP_OFF[(law, seed)]
        ok3 = abs(a_b-sa) <= 0.01 and abs(dn_b-sd) <= 0.1
        g3_ok &= ok3
        a_c, dn_c, po = read9(cb, axes_post=(3, 6, 7, 8))
        P(f"[{law} {seed}] G8Lb-3 baseline {a_b:.2f}/{dn_b:+.1f} vs "
          f"shipped {sa:.2f}/{sd:+.1f} -> {'PASS' if ok3 else 'FAIL'}")
        P(f"[{law} {seed}] CORRECTED: a={a_c:.2f} dN={dn_c:+.1f} "
          f"(baseline a={a_b:.2f} dN={dn_b:+.1f}; d_alpha="
          f"{a_c-a_b:+.3f})")
        P(f"[{law} {seed}]   P(fcomp)={np.round(po[3], 2).tolist()} "
          f"P(fpm)={np.round(po[6], 2).tolist()}")
        P(f"[{law} {seed}]   P(kw)={np.round(po[7], 2).tolist()} "
          f"P(sq)={np.round(po[8], 2).tolist()} "
          f"Dwob'={dwob(cb):+.1f} (baseline {dwob(eb):+.1f})")
        bx = np.unravel_index(np.nanargmax(cb), cb.shape)
        al = A_GRID[bx[0]]; eta = E_GRID[bx[1]]; wr = WR_GRID[bx[2]]
        fcm = FCOMP_GRID[bx[3]]; fpm = FPM_GRID[bx[6]]
        kw = KW_GRID[bx[7]]; sq = SQ_GRID[bx[8]]
        if pop is None:
            pop = build_pop(seed, lker=True)
        e_f = e_of_x(pop, eta, wr)
        tab_a = 1.0 + al*(TAB-1.0)
        vp_f = vp_c(pop, e_f, tab_a)
        o_f = run(pop['a_s'], e_f, pop['psi0'], pop['f_ip'],
                  pop['M_s'], pop['uph'], 8, 2500, 5, a0=A0_CAN,
                  tab=tab_a, lny0=LNY0, dlny=DLNY, vp=vp_f)
        mu_r, mh_r = band_mu(pop, o_f, fcm, kw, fpm, sq)
        res[(law, seed)] = dict(a_b=a_b, dn_b=dn_b, a_c=a_c,
                                dn_c=dn_c, cen=jointp(mu_r, mh_r),
                                mu=(mu_r, mh_r))
        P(f"[{law} {seed}] MAP cell: alpha={al}, wr={wr}, "
          f"fcomp={fcm}, fpm={fpm}, kw={kw}, sq={sq}; CENSUS "
          f"corrected-kernel: mu=({mu_r:.2f}, {mh_r:.2f}) "
          f"jointP={jointp(mu_r, mh_r):.2e}")
        P("")

if pend:
    P(f"INCOMPLETE - pending cubes: {pend}; verdict PENDING")
elif not g3_ok:
    P("G8Lb-3 FAIL - reader wiring suspect; DO NOT QUOTE")
else:
    P("G8Lb-3 4/4 PASS")
    P("")
    verdicts = []
    for law in ('simple', 'BE'):
        da = float(np.mean([res[(law, s)]['a_c']-res[(law, s)]['a_b']
                            for s in SEEDS]))
        dn = float(np.mean([res[(law, s)]['dn_c'] for s in SEEDS]))
        cat = ('WITHIN-SYSTEMATIC' if abs(da) <= 0.11 else
               'BOOSTIER-ALPHA' if da > 0.11 else 'QUIETER-ALPHA')
        verdicts.append((law, da, dn, cat))
        P(f"SIGN [{law}]: d_alpha = {da:+.3f} -> {cat}; honest "
          f"Newton band dN = {dn:+.1f} (seed mean)")
    reop = sum(1 for r in res.values() if r['cen'] >= 1e-3)
    P(f"CENSUS-RETEST: jointP >= 1e-3 in {reop}/4 -> "
      f"{'REOPENED' if reop >= 3 else 'NOT reopened'}")
    fifth = all(res[(l, s)]['a_c'] <= 0.2 and res[(l, s)]['dn_c'] <= 5
                for l in ('simple', 'BE') for s in SEEDS)
    P("")
    cats = {v[3] for v in verdicts}
    v = ("KERNEL-PAID: " + "; ".join(
        f"{law} {cat} (d_alpha {da:+.3f}, dN {dn:+.1f})"
        for law, da, dn, cat in verdicts))
    if reop >= 3:
        v += "; CENSUS-REOPENED at the corrected kernel"
    if fifth:
        v += ("; fifth-move-shaped movement observed - VERDICT "
              "DEFERRED to a powered round")
    P(f"==> 8L-b VERDICT (locked bars): {v}")
    P("    NO credence movement (pre-stated; external-only per "
      "8K-b/round 6).")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage8lb_read.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage8lb_read.txt")
