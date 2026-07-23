"""
STAGE 6D (O16, the reconciliation instrument): the DRIVE-WEIGHTED bath.

6A/6B/5Y sharpened the two-system tension to its final form: the
hierarchical galaxies reward the derived sharp-screening functions
monotonically (F4: -64 vertical, the largest controlled lead ever) while
the binaries -- probing the SAME functions through the dominant external
field -- reject every one of them (-5..-8.5 lnL, a0 translation +4..+6
sigma; only pure BE passes). No nu(y) reshaping can satisfy both,
because the binaries measure the function's behavior under the EFE.

The hypothesis with a mechanism (the O13 configuration-dependence, now
concrete): THE DRIVE IS WHICHEVER FIELD DOMINATES. A mode prepared by
the ambient Galactic field is externally driven -- its response
admixture is suppressed regardless of local occupation -- while a
self-sourced galaxy point keeps the spontaneous-share admixture. The
simplest realization: multiply the two-leg share by the local
self-dominance weight,
    beta = w_self * (1/2)/(2 nu - 1)^2,
    w_self = g_int/(g_int + g_ext)   (Newtonian magnitudes, pointwise).
Limits: w_self = 1 -> F4 exactly (all 6A galaxy results carry over --
SPARC field points have e_N ~ 0.01-0.05 a0, w_self ~ 1; the 5B/5E
environmental control bounds the neglected correction at the few-percent
level); w_self = 0 -> the pure occupation law exactly. This is a
NON-LOCAL modification -- nu is no longer a function of |g_N| alone but
of the ambient/local decomposition -- the same character as Milgrom's
modified-inertia EFE and the 4K/4L doors; stated openly.

This stage: the solver variant (one changed line, barycenter response =
BE(e_N) exactly at w_self = 0), gates (G1 isolated identity -> F4; G2
sandwich: the e_N = 1.2 boost must lie between the BE and F4 tables;
G3 wide-limit approach to BE), the binary 6-seed test (expect ~ p050 if
the reconciliation is right), kappa tables + the a0 translation (expect
BE-grade ~ +2 sigma, resolving the temperature strain).
Writes data/stage6d_summary.txt + data/stage6d_verdict.txt.
"""
import math, os, re, sys, time
import numpy as np
from numpy.polynomial.legendre import leggauss, legval, legder

NR, NMU, LMAX = 512, 96, 16
r = np.logspace(-2, 3, NR)
mu, wmu = leggauss(NMU)
R, MU = np.meshgrid(r, mu, indexing='ij')
ST = np.sqrt(1-MU**2)
Pl = np.zeros((LMAX+1, NMU)); dPl = np.zeros((LMAX+1, NMU))
for l in range(LMAX+1):
    c = np.zeros(l+1); c[l] = 1
    Pl[l] = legval(mu, c)
    dPl[l] = legval(mu, legder(c)) if l > 0 else np.zeros(NMU)

def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))

def nu_dw(y, w):
    """drive-weighted two-leg bath: beta = w * 0.5/(2 nu - 1)^2."""
    y = np.clip(np.asarray(y, float), 1e-14, None)
    w = np.broadcast_to(np.asarray(w, float), y.shape)
    ly = np.log(y)
    nu = nu_simple(y)
    for _ in range(80):
        b = w*0.5/((2.0*nu - 1.0)**2)
        db = w*(-2.0)/((2.0*nu - 1.0)**3)
        u = np.exp(np.minimum(0.5*(1.0+b)*ly + b*np.log(nu), 60.0))
        eu = np.exp(np.minimum(u, 60.0))
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        F = nu - 1.0 - n
        dudnu = u*(db*(0.5*ly + np.log(nu)) + b/nu)
        dF = 1.0 + (eu/(em1*em1))*dudnu
        step = F/np.where(np.abs(dF) > 1e-12, dF, 1e-12)
        nu = np.maximum(nu - np.clip(step, -0.5*nu, 0.5*nu), 1.0 + 1e-15)
    return nu

def solve_dw(eN):
    gr = -1.0/R**2 + eN*MU
    gt = -eN*ST
    wself = (1.0/R**2)/(1.0/R**2 + eN)
    f = nu_dw(np.hypot(gr, gt), wself) - 1.0
    Ar, At = f*gr, f*gt
    dAr = np.gradient(R**2*Ar, np.log(r), axis=0)/R**3
    dAt = -np.gradient(ST*At, mu, axis=1)/R
    S = -(dAr + dAt)
    Sl = np.array([(2*l+1)/2*np.sum(wmu*S*Pl[l], axis=1)
                   for l in range(LMAX+1)])
    A = np.zeros((LMAX+1, NR)); B = np.zeros((LMAX+1, NR))
    dr = np.diff(r)
    for l in range(LMAX+1):
        g_in = Sl[l]*r
        for i in range(NR-1):
            q = (r[i]/r[i+1])**(l+1)
            A[l, i+1] = A[l, i]*q + 0.5*dr[i]*(g_in[i]*q + g_in[i+1])
        g_out = Sl[l]*r
        for i in range(NR-2, -1, -1):
            q = (r[i]/r[i+1])**l
            B[l, i] = B[l, i+1]*q + 0.5*dr[i]*(g_out[i] + g_out[i+1]*q)
    lv = np.arange(LMAX+1)[:, None]
    gph_r = -(((lv+1)*A - lv*B)/(r[None, :]*(2*lv+1)))
    g_r = gr + np.tensordot(gph_r, Pl, axes=(0, 0))
    # barycenter response: at infinity w_self = 0 -> the pure occupation law
    nuE = float(nu_be(np.array([max(eN, 1e-14)]))[0]) if eN > 0 else 0.0
    if eN > 0:
        g_r -= nuE*eN*MU
    boost = np.sum(wmu*(-g_r), axis=1)/2/(1.0/r**2)
    return boost

y = 1.0/r**2
win = (y > 0.05) & (y < 30)
L = ["STAGE 6D: the drive-weighted bath -- beta = w_self * (1/2)/(2nu-1)^2,"
     " w_self = g_int/(g_int+g_ext)", ""]

# G1: isolated identity -> F4
b0 = solve_dw(0.0)
f4 = nu_dw(y, 1.0)
g1 = float(np.max(np.abs(b0[win]/f4[win] - 1.0)))
L.append(f"G1 isolated identity (w=1 -> F4): {100*g1:.2f}% -> "
         f"{'PASS' if g1 < 0.02 else 'FAIL'}")
assert g1 < 0.02

# eN = 1.2 table + G2 sandwich + G3 wide-limit
t0 = time.time()
b12 = solve_dw(1.2)
np.save('data/efe_boost_dwf_g1p2.npy', np.stack([y, b12]))
ybe, bbe = np.load('data/efe_boost_be_g1p2.npy')
yrb, brb = np.load('data/efe_boost_rb4_g1p2.npy')
lo = np.minimum(bbe, brb) - 0.015
hi = np.maximum(bbe, brb) + 0.015
g2 = bool(np.all((b12[win] > lo[win]) & (b12[win] < hi[win])))
L.append(f"G2 sandwich (between BE and F4 tables at eN=1.2, 1.5% slack): "
         f"{'PASS' if g2 else 'FAIL'}")
assert g2
iwide = np.where(y < 0.1)[0]
g3 = float(np.max(np.abs(b12[iwide]/bbe[iwide] - 1.0)))
L.append(f"G3 wide-limit approach to BE (y<0.1): max dev {100*g3:.2f}%")
L.append("  B(y) [BE / drive-weighted / F4] at eN=1.2:")
for yq in (0.1, 0.3, 1.0, 3.0):
    i = int(np.argmin(np.abs(y - yq)))
    L.append(f"    y={yq:4.1f}: {bbe[i]:.4f} / {b12[i]:.4f} / {brb[i]:.4f}")
for eN, tag in ((1.0, 'g1p0'), (1.4, 'g1p4')):
    np.save(f'data/efe_boost_dwf_{tag}.npy', np.stack([y, solve_dw(eN)]))
L.append(f"({(time.time()-t0)/60:.1f} min tables)")
print("\n".join(L), flush=True)

# ---------------- binary fits (4X patch set) ----------------
src = open('calcs/stage3p_v7budget.py', encoding='utf-8-sig').read()

OLD_DATA = """vt_d = (4.74047/plx[ok]*np.hypot(d['pmra1'][ok]-d['pmra2'][ok],
        d['pmdec1'][ok]-d['pmdec2'][ok]))/(0.9417*np.sqrt(Mtot_d/s_d))
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
vx_ = d['pmra2']-d['pmra1']; vy_ = d['pmdec2']-d['pmdec1']"""
NEW_DATA = """dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
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
rvs_ = (np.where(h1_, r1_, 0.0)*w1_ + np.where(h2_, r2_, 0.0)*w2_) \\
       / np.maximum(w1_+w2_, 1e-12)
th_ = sep/(2.06265e8/plx)
pmcor_ = rvs_*th_*plx/4.74047
vx_ = d['pmra2']-d['pmra1'] + pmcor_*sux_
vy_ = d['pmdec2']-d['pmdec1'] + pmcor_*suy_
vt_d = (4.74047/plx[ok]*np.hypot(vx_[ok], vy_[ok]))/(0.9417*np.sqrt(Mtot_d/s_d))"""
assert src.count(OLD_DATA) == 1
src = src.replace(OLD_DATA, NEW_DATA)

OLD_TABS = """if GTAG == '1p9':
    PS, PB = 'data/efe_boost_simple.npy', 'data/efe_boost_be.npy'
else:
    PS = f'data/efe_boost_simple_g{GTAG}.npy'
    PB = f'data/efe_boost_be_g{GTAG}.npy'
LNY_U, TAB_S = load_tab(PS)
_,     TAB_B = load_tab(PB)"""
NEW_TABS = """LNY_U, TAB_S = load_tab(LAMPATH)
TAB_B = TAB_S"""
assert src.count(OLD_TABS) == 1
src = src.replace(OLD_TABS, NEW_TABS)

OLD_LOOP = 'for law, TAB in (("simple", TAB_S), ("BE", TAB_B)):'
NEW_LOOP = 'for law, TAB in ((LAMLAW, TAB_S),):'
assert src.count(OLD_LOOP) == 1
src = src.replace(OLD_LOOP, NEW_LOOP)

OLD_SUM = """    P(f"  seed {seed}: BE-minus-simple best lnL = "
      f"{best_lnl['BE']-best_lnl['simple']:+.1f}  "
      f"({(time.time()-t0)/60:.1f} min)")"""
NEW_SUM = """    P(f"  seed {seed}: best lnL = {best_lnl[LAMLAW]:+.3f}  "
      f"({(time.time()-t0)/60:.1f} min)")"""
assert src.count(OLD_SUM) == 1
src = src.replace(OLD_SUM, NEW_SUM)

OLD_OUT = "with open('data/stage3u_summary.txt', 'a') as f:"
NEW_OUT = "with open('data/stage6d_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1
src = src.replace(OLD_OUT, NEW_OUT)

have = os.path.exists('data/stage6d_summary.txt') and \
       open('data/stage6d_summary.txt').read().count(' dwf: ') >= 6
if not have:
    if os.path.exists('data/stage6d_summary.txt'):
        os.remove('data/stage6d_summary.txt')
    SEEDS6 = ['31', '101', '202', '303', '404', '505']
    ns2 = {'__name__': '__main__',
           'LAMPATH': 'data/efe_boost_dwf_g1p2.npy', 'LAMLAW': 'dwf'}
    sys.argv = ['stage6d', '1p2'] + SEEDS6
    print("\n===== dwf (drive-weighted) on binaries =====", flush=True)
    exec(compile(src, 'stage3p_patched_6d_dwf', 'exec'), ns2)
else:
    print("(dwf fits found in data/stage6d_summary.txt -- skipping)",
          flush=True)

# ---------------- verdict: vs p050 per seed + kappa/a0 ----------------
def parse(path):
    txt = open(path).read()
    out = {}
    for m in re.finditer(
            r'seed (\d+) (\w+): a_hat=([0-9.]+) \(grid [0-9.]+, '
            r'interior=(\w+)\), dlnL\(Newton\)=\+([0-9.]+), wr=([0-9.]+),'
            r'.*?\n\s*seed \1: best lnL = (-?[0-9.]+)', txt):
        s, law, ah, inter, dn, wr, lnl = m.groups()
        out[(law, int(s))] = dict(a=float(ah), interior=(inter == 'True'),
                                  lnl=float(lnl))
    return out

REC = {}
for p in ('data/stage5k_summary.txt', 'data/stage5o_summary.txt',
          'data/stage6d_summary.txt'):
    REC.update(parse(p))
SEEDS = [31, 101, 202, 303, 404, 505]
ds = np.array([REC[('dwf', s)]['lnl'] - REC[('p050', s)]['lnl']
               for s in SEEDS])
ah = np.array([REC[('dwf', s)]['a'] for s in SEEDS])
n_int = sum(REC[('dwf', s)]['interior'] for s in SEEDS)
L.append("")
L.append("binary verdict (vs p050 per seed): " +
         " ".join(f"{d:+.2f}" for d in ds))
L.append(f"  mean {ds.mean():+.2f} +- {ds.std(ddof=1)/np.sqrt(6):.2f} SE; "
         f"better in {int((ds > 0).sum())}/6; interior {n_int}/6; "
         f"a_hat {ah.mean():.3f} +- {ah.std(ddof=1)/np.sqrt(6):.3f}")
L.append("  [comparators: gm -8.50, rb2 -5.47, rb3 -6.84; "
         "reconciliation target: within ~ -2 of p050]")

# kappa + a0 (5L machinery)
def load_tab2(path):
    t = np.load(path)
    yy, b = t[0][::-1], t[1][::-1]
    return np.log(yy), np.log(np.maximum(b-1.0, 1e-12))
tabs = {}
for fam in ('be', 'dwf'):
    for e, tag in ((1.0, 'g1p0'), (1.2, 'g1p2'), (1.4, 'g1p4')):
        tabs[(fam, e)] = load_tab2(f'data/efe_boost_{fam}_{tag}.npy')
def lnB1(fam, e, lny):
    lo_, hi_ = (1.0, 1.2) if e <= 1.2 else (1.2, 1.4)
    w = (e-lo_)/(hi_-lo_)
    la = np.interp(lny, *tabs[(fam, lo_)])
    lb = np.interp(lny, *tabs[(fam, hi_)])
    return (1-w)*la + w*lb
def kappa(fam, lny, e=1.184, dl=0.10, de=0.08):
    dy = (lnB1(fam, e, lny+dl) - lnB1(fam, e, lny-dl))/(2*dl)
    dE = (lnB1(fam, e+de, lny) - lnB1(fam, e-de, lny))/(2*de*1.0/e)
    return -(dy + dE)

from astropy.io import fits
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
Mtot = np.interp(MG1[ok], MG_T, MS_T)+np.interp(MG2[ok], MG_T, MS_T)
s_kau = sep[ok]/1e3
deep = (s_kau >= 6) & (s_kau < 30)
yd = (0.8868*Mtot[deep]/s_kau[deep]**2)*6.685e-9/1.2e-10
lnyd = np.log(np.clip(yd, 1e-3, 50.0))
kbe = float(np.mean(kappa('be', lnyd)))
kdw = float(np.mean(kappa('dwf', lnyd)))
okk = abs(kbe - 0.916) < 0.01
L.append(f"kappa: be {kbe:+.3f} (regression vs 4V +0.916 -> "
         f"{'PASS' if okk else 'FAIL'}), dwf {kdw:+.3f}")
m = float(ah.mean())
stot = math.hypot(float(ah.std(ddof=1)/np.sqrt(6)), 0.11)
a0b = 1.2e-10*m**(1.0/kdw)
sa0b = a0b*stot/(m*kdw)
tp = 2.998e8*(67.4*1e3/3.0857e22)/(2*math.pi)
pull = (a0b - tp)/math.hypot(sa0b, tp*0.5/67.4)
L.append(f"a0 translation: alpha {m:.3f}+-{stot:.3f}, kappa {kdw:+.3f} -> "
         f"a0 = {a0b*1e10:.2f}+-{sa0b*1e10:.2f}  pull {pull:+.1f} sigma  "
         f"[BE +1.9; sharp functions +4.3..+6.3]")

out = "\n".join(L)
print("\n" + out)
with open('data/stage6d_verdict.txt', 'w') as f:
    f.write(out + "\n")
print("\nSTAGE 6D done -> data/stage6d_verdict.txt")
