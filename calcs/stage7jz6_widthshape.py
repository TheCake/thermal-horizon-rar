"""Stage 7J-z6 Part A — the width-shape DIAGNOSTIC (pre-reg in NOTES,
committed before execution).

The sky rides fpm = 3.0 (past the Lindegren ceiling) while no injected sky
does.  This diagnostic asks WHERE the demand lives, conditionally at the
sky's PROF mode cells (seed 31, per law):

  A1  error-half split   — within each s-bin, pairs split at the bin-median
      formal error sg0; data histograms AND model noise pools split
      identically.  LOW-half-carried = error-independent (FLOOR family);
      HIGH-half-carried = error-correlated (TAIL family).
  A2  cell attribution   — per-(v-row, gamma-column) decomposition of the
      data*log(pp) gain fpm 2.4 -> 3.0 at the mode cell; per-bin shares.
  A3  magnitude split    — halves by G_faint (the Lindegren axis).

Gate G0-diag: the recomputed full-data conditional lnL at the mode cell
must reproduce the stored cube value (|d| <= 0.02) before any split is
read.  Conditional previews steer the Part B contest; they decide nothing
(the 5B coordinate-descent lesson).

Output: data/stage7jz6_diag.txt
"""
import numpy as np, os, time
from astropy.io import fits

OUT = 'data/stage7jz6_diag.txt'
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

# --- catalog (verbatim from stage7j_marginal.py) --------------------------
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
Mtot_d = np.interp(MG1[ok],MG_T,MS_T)+np.interp(MG2[ok],MG_T,MS_T)
s_d = sep[ok]/1e3
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
sn_ = np.maximum(np.hypot(sx_, sy_), 1e-12)
sux_, suy_ = sx_/sn_, sy_/sn_
_cn = list(d.columns.names)
def _pick(*names):
    for n in names:
        if n in _cn: return np.asarray(d[n], dtype=np.float64)
    raise KeyError(names)
r1_ = _pick('radial_velocity1', 'dr2_radial_velocity1', 'rv1')
r2_ = _pick('radial_velocity2', 'dr2_radial_velocity2', 'rv2')
try:
    er1_ = _pick('radial_velocity_error1', 'dr2_radial_velocity_error1')
    er2_ = _pick('radial_velocity_error2', 'dr2_radial_velocity_error2')
except KeyError:
    er1_ = np.full(len(r1_), 2.0); er2_ = np.full(len(r2_), 2.0)
h1_, h2_ = np.isfinite(r1_), np.isfinite(r2_)
w1_ = np.where(h1_, 1.0/np.maximum(er1_, 0.5)**2, 0.0)
w2_ = np.where(h2_, 1.0/np.maximum(er2_, 0.5)**2, 0.0)
rvs_ = (np.where(h1_, r1_, 0.0)*w1_ + np.where(h2_, r2_, 0.0)*w2_) \
       / np.maximum(w1_+w2_, 1e-12)
th_ = sep/(2.06265e8/plx)
pmcor_ = rvs_*th_*plx/4.74047
vx_ = d['pmra2']-d['pmra1'] + pmcor_*sux_
vy_ = d['pmdec2']-d['pmdec1'] + pmcor_*suy_
vt_d = (4.74047/plx[ok]*np.hypot(vx_[ok], vy_[ok]))/(0.9417*np.sqrt(Mtot_d/s_d))
e_vx = np.sqrt(d['pmra_error1']**2+d['pmra_error2']**2)
e_vy = np.sqrt(d['pmdec_error1']**2+d['pmdec_error2']**2)
cosg = np.abs(sx_*vx_+sy_*vy_)/np.maximum(np.hypot(sx_,sy_)*np.hypot(vx_,vy_),
                                          1e-12)
gam_d = np.degrees(np.arccos(np.clip(cosg, 0, 1)))[ok]
sig_c = ((4.74047/plx)*np.sqrt(0.5*(e_vx**2+e_vy**2)))
SBINS = [(0.2,2),(2,6),(6,20),(20,50)]
SC2 = (np.array([0.63, 3.46, 11.0, 31.6])/31.6)**2
VE = np.logspace(np.log10(0.02), np.log10(6.0), 21)
GE = np.linspace(0, 90, 7)
NV, NG = 20, 6
sig_ok = sig_c[ok]
Gfaint = np.maximum(G1m, G2m)[ok]
vc_ok = 0.9417*np.sqrt(Mtot_d/s_d)
vcen = np.sqrt(VE[:-1]*VE[1:]); gcen = 0.5*(GE[:-1]+GE[1:])
FLY = np.zeros((NV, NG))
for i in range(NV):
    for j in range(NG):
        FLY[i,j] = (1.0 if (0.5 <= vcen[i] <= 3.0) else 0.0) \
                   * max(30.0-gcen[j], 0.0)
FLY /= max(FLY.sum(), 1e-12)
UNI = np.ones((NV, NG))/(NV*NG)

def data_pack(mask_extra=None):
    """Per-bin data histogram + noise pool + mix templates for a subsample."""
    packs = []
    for bi, b in enumerate(SBINS):
        m = (s_d>=b[0])&(s_d<b[1])
        if mask_extra is not None:
            m = m & mask_extra
        h,_,_ = np.histogram2d(np.clip(vt_d[m],0.021,5.9), gam_d[m],
                               bins=[VE, GE])
        pool = sig_ok[m]
        cutp = 2.978/np.sqrt(s_d[m]) + 2.8284*sig_ok[m]
        acc = np.array([(vcen[i]*vc_ok[m] <= cutp).mean() if m.sum() else 0.0
                        for i in range(NV)])
        tpls = []
        for tpl in (UNI, FLY):
            t = tpl*acc[:,None]
            tpls.append(t/max(t.sum(), 1e-12))
        packs.append(dict(h=h.astype(float), pool=pool, uni=tpls[0],
                          fly=tpls[1], n=int(m.sum())))
    return packs

# per-bin split masks (data side)
split_masks = {'full': None}
for tag, arr in (('sg0', sig_ok), ('Gf', Gfaint)):
    lo = np.zeros(len(s_d), bool); hi = np.zeros(len(s_d), bool)
    for b in SBINS:
        m = (s_d>=b[0])&(s_d<b[1])
        med = np.median(arr[m])
        lo |= m & (arr <= med); hi |= m & (arr > med)
    split_masks[f'{tag}-lo'] = lo
    split_masks[f'{tag}-hi'] = hi
packs_of = {k: data_pack(v) for k, v in split_masks.items()}

# --- population (verbatim physics, seed 31) --------------------------------
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
    p['pick'] = [rng.integers(0, max(len(packs_of['full'][bi]['pool']),1), N)
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

def bin_lnL(p, o, packs, mode, fpm, return_cells=False):
    """Conditional lnL over bins at the mode cell, one fpm value.
    mode = dict(fcm, kw, sq, fc, ff).  Model pools come from packs."""
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
    tot = 0.0; nbin = []
    cells = []
    for bi, b in enumerate(SBINS):
        pk = packs[bi]
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500 or len(pk['pool']) == 0:
            nbin.append(0.0); cells.append(None); continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = pk['pool'][p['pick'][bi][idx] % len(pk['pool'])]/4.74047
        cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
        mh_tot = np.zeros(len(idx))
        for k in (1, 2):
            c = p['comp'][k]
            act = c['uc'][idx] < mode['fcm']
            mh_tot += act*c['mh'][idx]
            cvp += act*c['w'][idx]*c['wd'][idx,0]
            cvq += act*c['w'][idx]*c['wd'][idx,1]
        boost = np.sqrt(1+mh_tot/p['M_s'][idx])
        vp_b = vpar[idx] + mode['kw']*cvp
        vq_b = vper[idx] + mode['kw']*cvq
        vp_n = vp_b*boost + p['gn1'][idx]*sg0*fpm
        vq_n = vq_b*boost + p['gn2'][idx]*sg0*fpm
        vmag = np.hypot(vp_n, vq_n)
        keep = vmag*4.74047 <= (2.978/np.sqrt(s_kau[idx])
                                + 2.8284*sg0*4.74047)
        vtn = (vmag/vc)[keep]
        gmn = np.degrees(np.arccos(np.clip(
            np.abs(vp_n[keep])/np.maximum(vmag[keep],1e-12), 0, 1)))
        vts = vtn*np.exp(mode['sq']*p['gs'][idx][keep])
        h,_,_ = np.histogram2d(np.clip(vts,0.021,5.9), gmn, bins=[VE, GE])
        p0 = np.maximum(h/max(h.sum(),1), 1e-5); p0 /= p0.sum()
        wch = min(mode['fc']*SC2[bi], 0.5); wfl = min(mode['ff']*SC2[bi], 0.5)
        wtot = min(wch+wfl, 0.6)
        mixc = (wch*pk['uni'] + wfl*pk['fly'])/(wch+wfl)
        pp = (1-wtot)*p0 + wtot*mixc
        l = float(np.sum(pk['h']*np.log(pp)))
        tot += l; nbin.append(l)
        cells.append(pk['h']*np.log(pp) if return_cells else None)
    return tot, nbin, cells

t0 = time.time()
P("STAGE 7J-z6 PART A: width-shape diagnostic (conditional at PROF mode, "
  "seed 31; pre-reg in NOTES 2026-07-27)")
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
p = build_pop(31)

for law, TAB in (("simple", TAB_S), ("BE", TAB_B)):
    cube = np.load(f'data/stage7j_cube_full_photow3_31_{law}.npy')
    cb = cube + prior_eta.reshape((1, len(E_GRID)) + (1,)*(cube.ndim-2))
    bx = np.unravel_index(np.nanargmax(cb), cb.shape)
    ai, ei, wi, fi, ci, yi, pi, ki, si = bx
    mode = dict(fcm=FCOMP_GRID[fi], kw=KW_GRID[ki], sq=SQ_GRID[si],
                fc=FC0_GRID[ci], ff=FFLY_GRID[yi])
    al, eta, wr = A_GRID[ai], E_GRID[ei], WR_GRID[wi]
    P(f"\n== {law}: mode cell alpha={al}, eta={eta}, wr={wr}, "
      f"fcomp={mode['fcm']}, kw={mode['kw']}, sq={mode['sq']}, "
      f"ff={mode['ff']}, fpm={FPM_GRID[pi]} (cube lnL {cube[bx]:.2f})")
    tab_a = 1.0 + al*(TAB-1.0)
    e_s = e_of(p, eta, wr)
    if al > 0:
        vp = vp_c(p, e_s, tab_a)
        o = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'], p['uph'],
                8, 2500, 5, a0=A0_CAN, tab=tab_a, lny0=LNY0, dlny=DLNY, vp=vp)
    else:
        o = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'], p['uph'],
                8, 2500, 1)

    # G0-diag: recompute the mode cell on the FULL data, compare to cube
    lC, _, _ = bin_lnL(p, o, packs_of['full'], mode, FPM_GRID[pi])
    dG0 = abs(lC - float(cube[bx]))
    P(f"G0-diag {law}: recomputed {lC:.2f} vs cube {cube[bx]:.2f} "
      f"-> |d| = {dG0:.3f} -> {'PASS' if dG0 <= 0.02 else 'FAIL'}")
    if dG0 > 0.02:
        P(f"ABORT {law}: G0-diag failed — wiring mismatch, splits not read")
        continue

    # A2: cell attribution fpm 2.4 -> 3.0 on full data
    l24, nb24, c24 = bin_lnL(p, o, packs_of['full'], mode, 2.4,
                             return_cells=True)
    l30, nb30, c30 = bin_lnL(p, o, packs_of['full'], mode, 3.0,
                             return_cells=True)
    P(f"A2 {law}: lnL(fpm=3.0) - lnL(fpm=2.4) = {l30-l24:+.2f} total; "
      f"per-bin " + ", ".join(f"{b}:{n3-n2:+.1f}"
                              for b, n2, n3 in zip(SBINS, nb24, nb30)))
    gain = np.zeros((NV, NG))
    for cc24, cc30 in zip(c24, c30):
        if cc24 is not None:
            gain += cc30 - cc24
    vrow = gain.sum(axis=1); gcol = gain.sum(axis=0)
    pos = np.maximum(vrow, 0).sum()
    tail_share = np.maximum(vrow[vcen > 1.5], 0).sum()/max(pos, 1e-9)
    core_share = np.maximum(vrow[vcen < 0.5], 0).sum()/max(pos, 1e-9)
    P(f"A2 {law}: positive-gain shares — core(v<0.5) {core_share:.2f}, "
      f"mid(0.5-1.5) {1-core_share-tail_share:.2f}, "
      f"tail(v>1.5) {tail_share:.2f}")
    P(f"A2 {law}: v-row gains " +
      " ".join(f"{v:.2f}:{g:+.1f}" for v, g in zip(vcen, vrow)
               if abs(g) >= 0.5))
    P(f"A2 {law}: gamma-col gains " +
      " ".join(f"{g:.0f}d:{x:+.1f}" for g, x in zip(gcen, gcol)))

    # A1 + A3: split fpm profiles
    for tag, lab in (('sg0', 'A1 error'), ('Gf', 'A3 G_faint')):
        for half in ('lo', 'hi'):
            packs = packs_of[f'{tag}-{half}']
            prof = []
            for fpm in FPM_GRID:
                l, _, _ = bin_lnL(p, o, packs, mode, fpm)
                prof.append(l)
            prof = np.array(prof)
            im = int(np.argmax(prof))
            P(f"{lab} {half} {law}: fpm profile " +
              " ".join(f"{f}:{v-prof.max():+.1f}"
                       for f, v in zip(FPM_GRID, prof)) +
              f" -> peak {FPM_GRID[im]}"
              f"{' (EDGE)' if im == len(FPM_GRID)-1 else ''}, "
              f"d(3.0-2.4) = {prof[-1]-prof[-2]:+.2f}")

P(f"\ndone ({(time.time()-t0)/60:.1f} min)")
