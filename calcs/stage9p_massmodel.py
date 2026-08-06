"""STAGE 9P — THE MASS-MODEL INSENSITIVITY CONTROL ON THE FPM METER
(pre-registered BEFORE any run; robustness round; NO credence movement;
the RNAAS-note prerequisite — round-15 reviewer exposure (c-iv): "show
explicitly the >=2.1x inflation is insensitive to the M_G-mass relation
or selection").

DESIGN (locked at pre-reg):
The 9L narrow-pair meter (alpha = 0, law-blind, 0.2-2 kAU) is re-read
under deformations of the ONLY mass-model object in the chain, the
MG_T->MS_T interpolation table, plus one sample-window deformation.
The model-side orbit population is mass-table-independent (M_s is drawn
uniformly, vt normalized by its own mass), so one orbit run per seed is
reused across all variants; a variant changes (a) the data-side masses
Mtot_d hence vt_d and the contaminant-acceptance templates, and (b) the
companion wobble's mass->magnitude interpolation. All expressions are
BIT-VERBATIM 9L except the table parameterization. The pick stream is
generated at the FIDUCIAL pool length for every variant (9L's own
modulo re-maps it), keeping the rng streams identical across variants.

VARIANTS (amplitudes locked): V0 fiducial (identity gate); V1 global
mass scale x0.80; V2 global x1.25 (generous vs the ~5-10% photometric
M-L systematic for MS dwarfs); V3/V4 tilt MS_T' = MS_T*(MS_T/0.82)^(+-0.15)
(a +-10-14% differential across the MS, pivot 0.82 Msun at MG = 6.0);
V5 main-sequence window tightened to MG in [3.0, 13.0] both components
(data-side sample cut; RUWE quartile masks recomputed on the cut sample).
Parallax-convention variations are OUT OF SCOPE for this stage (named
limitation carried to the note).

GATES: G9P-0 identity — the V0 narrow tables must reproduce the
archived stage9l_tables_{seed}.npz at max|dT| <= 1e-4 lnL AND the
E[fpm] marginals recomputed from both tables must agree to 1e-6
(lnL-grade reader identity per the 8P/8Q rule). G9P-1 counts — every
variant/quartile narrow cell keeps >= 1000 pairs. The 9L G9L-2
boost-premise and G9L-3 injection gates are model-side and
table-independent; they are inherited, not re-run (stated in output).

BARS (ordered, read over V1-V5, all quartiles, both seeds):
  P-ROBUST  iff min E[fpm] >= 1.8  (the meter's noise-real band edge:
            the >=2.1x-class bound survives every deformation)
  P-FRAGILE iff any E[fpm] <= 1.5  (the bound is mass-model-carried;
            the note MUST NOT ship the 2.1x without mass casework)
  P-GRAY    otherwise (quote the envelope; the note ships the minimum
            as its bound)
NO credence movement in any branch (pre-stated; robustness round).
Output: data/stage9p_massmodel.txt + data/stage9p_grid.npz
"""
import math
import numpy as np, time
from astropy.io import fits
from scipy.special import logsumexp

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
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
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
vrel_kms = (4.74047/plx[ok]*np.hypot(vx_[ok], vy_[ok]))
e_vx = np.sqrt(d['pmra_error1']**2+d['pmra_error2']**2)
e_vy = np.sqrt(d['pmdec_error1']**2+d['pmdec_error2']**2)
cosg = np.abs(sx_*vx_+sy_*vy_)/np.maximum(np.hypot(sx_,sy_)*np.hypot(vx_,vy_),
                                          1e-12)
gam_d = np.degrees(np.arccos(np.clip(cosg, 0, 1)))[ok]
sig_c = ((4.74047/plx)*np.sqrt(0.5*(e_vx**2+e_vy**2)))
sig_ok = sig_c[ok]
MG1ok, MG2ok = MG1[ok], MG2[ok]
SBINS = [(0.2,2),(2,6),(6,20),(20,50)]
NB_BINS = 1
SC2 = (np.array([0.63, 3.46, 11.0, 31.6])/31.6)**2
VE = np.logspace(np.log10(0.02), np.log10(6.0), 21)
GE = np.linspace(0, 90, 7)
NV, NG = 20, 6
vcen = np.sqrt(VE[:-1]*VE[1:]); gcen = 0.5*(GE[:-1]+GE[1:])
FLY = np.zeros((NV, NG))
for i in range(NV):
    for j in range(NG):
        FLY[i,j] = (1.0 if (0.5 <= vcen[i] <= 3.0) else 0.0) \
                   * max(30.0-gcen[j], 0.0)
FLY /= max(FLY.sum(), 1e-12)
UNI = np.ones((NV, NG))/(NV*NG)

N = 500_000
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FC0 = 0.10
FFLY_GRID = np.array([0.05, 0.10])
FPM_GRID = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW_GRID = np.array([0.7, 1.0, 1.4])
SQ_GRID = np.array([0.0, 0.1, 0.2, 0.3])
WS2_GRID = np.array([0.0, 0.045])
SEEDS = (31, 101)

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

# ---------------- the variants (amplitudes LOCKED at pre-reg) -------
def table_variant(tag):
    if tag in ('V0-fid', 'V5-window'):
        return MS_T.copy()
    if tag == 'V1-low':
        return MS_T*0.80
    if tag == 'V2-high':
        return MS_T*1.25
    if tag == 'V3-tiltp':
        return MS_T*(MS_T/0.82)**(+0.15)
    if tag == 'V4-tiltm':
        return MS_T*(MS_T/0.82)**(-0.15)
    raise ValueError(tag)

VARIANTS = ['V0-fid', 'V1-low', 'V2-high', 'V3-tiltp', 'V4-tiltm',
            'V5-window']

def data_side(tag):
    """Variant data-side arrays: (Mtot, vt, vc_ok, window-mask)."""
    ms = table_variant(tag)
    Mtot = np.interp(MG1ok, MG_T, ms) + np.interp(MG2ok, MG_T, ms)
    vc = 0.9417*np.sqrt(Mtot/s_d)
    vt = vrel_kms/vc
    if tag == 'V5-window':
        win = (MG1ok > 3.0) & (MG1ok < 13.0) \
            & (MG2ok > 3.0) & (MG2ok < 13.0)
    else:
        win = np.ones(len(s_d), dtype=bool)
    return vt, vc, win

def build_stratum(pm, vt_v, vc_v):
    D2, PLs, UB, FB = [], [], [], []
    for b in SBINS:
        m = (s_d>=b[0])&(s_d<b[1])&pm
        h,_,_ = np.histogram2d(np.clip(vt_v[m],0.021,5.9), gam_d[m],
                               bins=[VE, GE])
        D2.append(h.astype(float))
        PLs.append(sig_ok[m])
        cutp = 2.978/np.sqrt(s_d[m]) + 2.8284*sig_ok[m]
        acc = np.array([(vcen[i]*vc_v[m] <= cutp).mean()
                        if m.sum() else 0.0 for i in range(NV)])
        for tpl, store in ((UNI, UB), (FLY, FB)):
            t = tpl*acc[:,None]
            store.append(t/max(t.sum(), 1e-12))
    ND = [int(h.sum()) for h in D2]
    return D2, PLs, UB, FB, ND

def build_pop(seed, ms_tab, pool_len_fid):
    """VERBATIM 9L build_pop; ms_tab parameterizes the wobble
    interpolation only; pick is drawn at the FIDUCIAL pool lengths
    (the eval modulo re-maps) so rng streams are variant-identical."""
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
    p['pick'] = [rng.integers(0, max(pool_len_fid[bi],1), N)
                 for bi in range(len(SBINS))]
    p['gn1'], p['gn2'] = rng.normal(size=N), rng.normal(size=N)
    M_h = 0.5*p['M_s']
    p['comp'] = {}
    for k in (1, 2):
        u_q = rng.random(N)
        q = 0.1+0.9*u_q
        logP = rng.normal(5.03, 2.28, N)
        P_yr = 10**logP/365.25
        a_in = (M_h*(1+q)*P_yr**2)**(1/3)
        valid = (a_in < 130.0) & (a_in < p['a_s']/5.0)
        v_orb = 29.78*np.sqrt(M_h*(1+q)/np.maximum(a_in,1e-3))
        u_ = np.pi*2.83/np.maximum(P_yr, 1e-9)
        S = np.where(u_ < 1e-2, 1.0 - u_*u_/10.0,
                     3.0*(np.sin(u_) - u_*np.cos(u_))
                     / np.maximum(u_, 1e-300)**3)
        MGp = np.interp(-np.clip(M_h, ms_tab[-1], ms_tab[0]), -ms_tab, MG_T)
        MGs = np.interp(-np.clip(q*M_h, ms_tab[-1], ms_tab[0]), -ms_tab, MG_T)
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
    erf = 0.95
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
    vpar = np.sum(vsky*b1,axis=1)
    vper = np.sum(vsky*b2,axis=1)
    return smag, vpar, vper

def eval_block_nb(p, prj, D2, PLs, UB, FB):
    """VERBATIM 9L eval_block_nb (narrow bins only)."""
    smag, vpar, vper = prj
    s_kau = smag/1e3
    out = np.zeros((len(FCOMP_GRID), len(FFLY_GRID), len(FPM_GRID),
                    len(KW_GRID), len(SQ_GRID), len(WS2_GRID)))
    for bi in range(NB_BINS):
        b = SBINS[bi]
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        if len(idx) < 500 or len(PLs[bi]) == 0: continue
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = PLs[bi][p['pick'][bi][idx] % len(PLs[bi])]/4.74047
        g1_i, g2_i = p['gn1'][idx], p['gn2'][idx]
        sk_i = s_kau[idx]
        gk_full = p['gs'][idx]
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
            for ki, kwv in enumerate(KW_GRID):
                vp_a = vpar[idx] + kwv*cvp
                vq_a = vper[idx] + kwv*cvq
                for pi, fpm in enumerate(FPM_GRID):
                    for wi, ws in enumerate(WS2_GRID):
                        if ws == 0.0:
                            vp_n = vp_a*boost + g1_i*sg0*fpm
                            vq_n = vq_a*boost + g2_i*sg0*fpm
                        else:
                            sig_eff = np.sqrt((sg0*fpm)**2
                                              + (ws/4.74047)**2)
                            vp_n = vp_a*boost + g1_i*sig_eff
                            vq_n = vq_a*boost + g2_i*sig_eff
                        vmag = np.hypot(vp_n, vq_n)
                        keep = vmag*4.74047 <= (2.978/np.sqrt(sk_i)
                                                + 2.8284*sg0*4.74047)
                        vtn = (vmag/vc)[keep]
                        gmn = np.degrees(np.arccos(np.clip(
                            np.abs(vp_n[keep])
                            / np.maximum(vmag[keep], 1e-12), 0, 1)))
                        gk = gk_full[keep]
                        for si, sqv in enumerate(SQ_GRID):
                            vts = vtn*np.exp(sqv*gk)
                            h,_,_ = np.histogram2d(
                                np.clip(vts,0.021,5.9), gmn,
                                bins=[VE, GE])
                            p0 = np.maximum(h/max(h.sum(),1), 1e-5)
                            p0 /= p0.sum()
                            for yi, ff in enumerate(FFLY_GRID):
                                wch = min(FC0*SC2[bi], 0.5)
                                wfl = min(ff*SC2[bi], 0.5)
                                wtot = min(wch+wfl, 0.6)
                                mixc = (wch*UB[bi]
                                        + wfl*FB[bi])/(wch+wfl)
                                pp = (1-wtot)*p0 + wtot*mixc
                                out[fi, yi, pi, ki, si, wi] += \
                                    np.sum(D2[bi]*np.log(pp))
    return out

ruwe = np.maximum(_pick('ruwe1', 'RUWE1'), _pick('ruwe2', 'RUWE2'))[ok]
qs_ = np.percentile(ruwe, [25, 50, 75])
QMASKS = [ruwe <= qs_[0],
          (ruwe > qs_[0]) & (ruwe <= qs_[1]),
          (ruwe > qs_[1]) & (ruwe <= qs_[2]),
          ruwe > qs_[2]]
ALLM = np.ones(len(s_d), dtype=bool)

t0 = time.time()
P("9P THE MASS-MODEL INSENSITIVITY CONTROL (pre-reg committed BEFORE "
  "any run; robustness round; NO credence movement; the RNAAS-note "
  "prerequisite)")
P("")
P("inherited gates (model-side, table-independent, NOT re-run): "
  "G9L-2 boost premise 8/8 PASS, G9L-3 injection 2/2 PASS (9L)")

# fiducial pool lengths (for the variant-identical pick stream)
vt_fid, vc_fid, _ = data_side('V0-fid')
pool_fid = build_stratum(ALLM, vt_fid, vc_fid)[1]
POOL_LEN_FID = [len(x) for x in pool_fid]

g0_ok = True; g1_ok = True
E = np.zeros((len(VARIANTS), 4, len(SEEDS)))
for si_, seed in enumerate(SEEDS):
    arch = np.load(f'data/stage9l_tables_{seed}.npz')['T']
    for vi, tag in enumerate(VARIANTS):
        pf = build_pop(seed, table_variant(tag), POOL_LEN_FID)
        if vi == 0:
            e_f = e_of_x(pf, 1.05, 0.30)
            tab0 = np.ones_like(TAB_S)
            vp0 = vp_c(pf, e_f, tab0)
            o0 = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'],
                     pf['M_s'], pf['uph'], 8, 2500, 5, a0=A0_CAN,
                     tab=tab0, lny0=LNY0, dlny=DLNY, vp=vp0)
            prj0 = project(pf, o0)
        vt_v, vc_v, win = data_side(tag)
        for qi in range(4):
            # AMENDMENT 1 (G9P-0 first-firing catch, pre-quote, logged):
            # the eval was wired to the full-sample noise pool where 9L
            # indexes each stratum's OWN pool (the pick modulo re-maps);
            # fixed to verbatim — the identity gate exists for exactly
            # this (GB0w precedent).
            D2, PLs, UB, FB, ND = build_stratum(QMASKS[qi] & win,
                                                vt_v, vc_v)
            nn = sum(int(D2[bi].sum()) for bi in range(NB_BINS))
            if nn < 1000:
                g1_ok = False
                P(f"G9P-1 FAIL: {tag} Q{qi+1} narrow count {nn}")
            T = eval_block_nb(pf, prj0, D2, PLs, UB, FB)
            if vi == 0:
                dmax = float(np.max(np.abs(T - arch[qi])))
                lw_a = logsumexp(arch[qi]
                                 + LNPI.reshape(6,1,1,1,1,1),
                                 axis=(0,1,3,5))
                w_a = np.exp(lw_a - logsumexp(lw_a))
                e_a = float(np.sum(w_a.sum(axis=1)*FPM_GRID))
                lw_v = logsumexp(T + LNPI.reshape(6,1,1,1,1,1),
                                 axis=(0,1,3,5))
                w_v = np.exp(lw_v - logsumexp(lw_v))
                e_v = float(np.sum(w_v.sum(axis=1)*FPM_GRID))
                okq = dmax <= 1e-4 and abs(e_v - e_a) <= 1e-6
                g0_ok &= okq
                P(f"[seed {seed}] G9P-0 identity Q{qi+1}: max|dT| = "
                  f"{dmax:.2e}, E[fpm] {e_v:.4f} vs arch {e_a:.4f} "
                  f"-> {'PASS' if okq else 'FAIL'}")
                E[vi, qi, si_] = e_v
            else:
                lw = logsumexp(T + LNPI.reshape(6,1,1,1,1,1),
                               axis=(0,1,3,5))
                wq = np.exp(lw - logsumexp(lw))
                E[vi, qi, si_] = float(np.sum(wq.sum(axis=1)
                                              * FPM_GRID))
        if vi > 0:
            P(f"[seed {seed}] {tag}: E[fpm] Q1-Q4 = "
              + "/".join(f"{E[vi,qi,si_]:.2f}" for qi in range(4))
              + f"  ({(time.time()-t0)/60:.1f} min)")
    P("")

if not (g0_ok and g1_ok):
    P("GATES FAILED (G9P-0 " + ('PASS' if g0_ok else 'FAIL')
      + ", G9P-1 " + ('PASS' if g1_ok else 'FAIL')
      + ") - STOP; DO NOT QUOTE; verdict WITHHELD")
else:
    P("GATES: G9P-0 8/8, G9P-1 - ALL PASS")
    P("")
    sub = E[1:, :, :]
    emin, emax = float(sub.min()), float(sub.max())
    q1min = float(E[1:, 0, :].min())
    P(f"P-facts: envelope over V1-V5 x Q1-Q4 x seeds = "
      f"[{emin:.2f}, {emax:.2f}]; Q1-only min = {q1min:.2f}; "
      f"fiducial Q1 = "
      + ", ".join(f"{E[0,0,s]:.2f}" for s in range(len(SEEDS))))
    if emin >= 1.8:
        P("==> 9P VERDICT (locked grammar): P-ROBUST - the >=2x "
          "pair-error bound survives every mass-table and window "
          "deformation at the meter's noise-real band; the "
          "calibration note ships with the mass-model objection "
          "closed (envelope quoted).")
    elif emin <= 1.5:
        P("==> 9P VERDICT (locked grammar): P-FRAGILE - the bound "
          "is mass-model-carried at one or more deformations; the "
          "note MUST NOT ship the 2.1x headline without mass "
          "casework; quote the envelope only.")
    else:
        P("==> 9P VERDICT (locked grammar): P-GRAY - envelope "
          "quoted; the note ships the envelope minimum as its "
          "bound.")
    P("    NO credence movement (pre-stated; robustness round).")
np.savez('data/stage9p_grid.npz', E=E, variants=np.array(VARIANTS),
         seeds=np.array(SEEDS))
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage9p_massmodel.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage9p_massmodel.txt")
