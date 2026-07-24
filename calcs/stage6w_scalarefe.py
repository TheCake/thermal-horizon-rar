"""
STAGE 6W (the Saturn round): SCALAR vs VECTOR EFE on the wide binaries.

The Cassini quadrupole is NOT a prediction of MOND-like dynamics per se;
it is a prediction of the VECTOR composition of internal and external
fields (the transition-shell phantom density acquires a P2(cos theta)
anisotropy from the direction of g_ext). Our own triangulation says the
ambient enters the physics as a SCALAR: 6D excluded the pointwise
field-ratio rule, 6T excluded local-field running, 6G accepted the
system-level scalar gate. If the EFE itself is thermodynamic (a bath
state, magnitude-only) rather than field composition:
  - the phantom shell is SPHERICAL -> by the shell theorem the inner
    solar system feels NOTHING from it: Q2(scalar) = the true galactic
    tide ~ g_gal/R_gal ~ 8e-31 s^-2, ~3.5 orders BELOW the Cassini cap
    (9e-27). The 4K tension is REMOVED at theory-class level, no
    modified inertia required.
  - the binaries then become the decisive instrument: they are the only
    data that probe the EFE's SHAPE at the transition. Scalar tables
    B(y) = nu(y + e) (magnitude composition) vs our vector tables
    B(y) = sphericalized QUMOND solve. Nobody has contested this on
    binary data (scout launched).

PRE-REGISTERED BARS (committed before execution; 6 seeds, corrected
velocities, the 6G patch machinery verbatim):
  For each pairing (sbe vs p050-vector; samb vs amb-vector):
    VIABLE    if mean Delta >= -3 (AMB-acceptance grade)
    PREFERRED if mean Delta > 0 by > 2 SE
    EXCLUDED  if mean Delta <= -5 (the eight-function band)
    interior alpha-hat required; edge = shape rejection.
  a0 row: PASS <= +2.5 sigma.
  CONSEQUENCE MAP (pre-stated): VIABLE or better in a pairing => the
  scalar-EFE theory class stands on all measured data (binaries here;
  galaxy EFE per 5B/5E is magnitude-level; Cassini passed by the shell
  theorem) => the quadrupole problem has a concrete non-MI resolution.
  EXCLUDED in both pairings => the binaries DEMAND the vector
  composition => Saturn stands; logged honestly.

Writes data/stage6w_summary.txt + data/stage6w_verdict.txt.
"""
import math, os, re, sys, time
import numpy as np

def n_amb_of(e):
    x = math.sqrt(e)
    return 1.0/(math.exp(x) - 1.0)
def g_of(n):
    return (n/(1.0 + n))**2

def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))
def nu_be(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    x = np.sqrt(y)
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))

def make_amb(g):
    def nu_run(y):
        y = np.clip(np.asarray(y, float), 1e-14, None)
        ly = np.log(y)
        nu = nu_simple(y)
        for _ in range(80):
            b = g*0.5/((2.0*nu - 1.0)**2)
            db = g*(-2.0)/((2.0*nu - 1.0)**3)
            u = np.exp(np.minimum(0.5*(1.0+b)*ly + b*np.log(nu), 60.0))
            eu = np.exp(np.minimum(u, 60.0))
            em1 = np.maximum(eu - 1.0, 1e-300)
            n = np.where(u < 60.0, 1.0/em1, 0.0)
            F = nu - 1.0 - n
            dudnu = u*(db*(0.5*ly + np.log(nu)) + b/nu)
            dF = 1.0 + (eu/(em1*em1))*dudnu
            step = F/np.where(np.abs(dF) > 1e-12, dF, 1e-12)
            nu = np.maximum(nu - np.clip(step, -0.5*nu, 0.5*nu),
                            1.0 + 1e-15)
        return nu
    return nu_run

# ---- scalar tables (analytic: B(y) = nu(y + e)) + shape diagnostics --------
yg = np.logspace(-4, 4, 1000)[::-1]          # match table convention (y desc)
print("scalar-EFE tables (B(y) = nu(y+e); magnitude composition):",
      flush=True)
for eN, tag in ((1.0, 'g1p0'), (1.2, 'g1p2'), (1.4, 'g1p4')):
    Bs = nu_be(yg + eN)
    np.save(f'data/efe_boost_sbe_{tag}.npy', np.stack([yg, Bs]))
    ga = g_of(n_amb_of(eN))
    Ba = make_amb(ga)(yg + eN)
    np.save(f'data/efe_boost_samb_{tag}.npy', np.stack([yg, Ba]))
    # G1: e -> 0 regression (analytic identity)
    if tag == 'g1p2':
        d0 = np.max(np.abs(nu_be(yg + 1e-12)/nu_be(yg) - 1.0))
        print(f"  G1 (e->0 -> isolated nu; probe eps=1e-12 itself "
              f"contributes ~5e-9 at the deep end): max rel {d0:.1e} "
              f"{'PASS' if d0 < 1e-7 else 'FAIL'}", flush=True)
        assert d0 < 1e-7
        # shape diagnostic vs the vector table in the data window
        t = np.load('data/efe_boost_be_g1p2.npy')
        yv, Bv = t[0], t[1]
        win = (yv > 0.05) & (yv < 30)
        Bs12 = nu_be(yv + 1.2)
        D = np.log(np.maximum(Bs12 - 1, 1e-12)) \
            - np.log(np.maximum(Bv - 1, 1e-12))
        print(f"  shape diff ln(B-1) scalar-vector over y in [0.05,30]: "
              f"min {D[win].min():+.3f}, max {D[win].max():+.3f}, "
              f"mean {D[win].mean():+.3f}  (the discriminating signal)",
              flush=True)
print("tables written.", flush=True)

# ---- the v7 budget, 6G patch verbatim, two scalar laws ---------------------
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
NEW_OUT = "with open('data/stage6w_summary.txt', 'a') as f:"
assert src.count(OLD_OUT) == 1
src = src.replace(OLD_OUT, NEW_OUT)

SEEDS6 = ['31', '101', '202', '303', '404', '505']
def have(lawtag):
    return os.path.exists('data/stage6w_summary.txt') and \
        open('data/stage6w_summary.txt').read().count(f' {lawtag}: ') >= 6
for lawtag in ('sbe', 'samb'):
    if have(lawtag):
        continue
    ns2 = {'__name__': '__main__',
           'LAMPATH': f'data/efe_boost_{lawtag}_g1p2.npy', 'LAMLAW': lawtag}
    sys.argv = ['stage6w', '1p2'] + SEEDS6
    print(f"\n===== {lawtag} (scalar EFE) on binaries =====", flush=True)
    exec(compile(src, f'stage3p_patched_6w_{lawtag}', 'exec'), ns2)

def parse(path):
    txt = open(path).read()
    out = {}
    for m in re.finditer(
            r'seed (\d+) (\w+): a_hat=([0-9.]+) \(grid [0-9.]+, '
            r'interior=(\w+)\), dlnL\(Newton\)=\+([0-9.]+), wr=([0-9.]+),'
            r'.*?\n\s*seed \1: best lnL = (-?[0-9.]+)', txt):
        s, law, ah, inter, dn, wr, lnl = m.groups()
        out[(law, int(s))] = dict(a=float(ah), interior=(inter == 'True'),
                                  dn=float(dn), lnl=float(lnl))
    return out

REC = {}
for p in ('data/stage5k_summary.txt', 'data/stage5o_summary.txt',
          'data/stage6g_summary.txt', 'data/stage6w_summary.txt'):
    REC.update(parse(p))
SEEDS = [31, 101, 202, 303, 404, 505]

L = ["STAGE 6W: scalar vs vector EFE on the binaries", ""]
verdicts = {}
for lawtag, ref in (('sbe', 'p050'), ('samb', 'amb')):
    ds = np.array([REC[(lawtag, s)]['lnl'] - REC[(ref, s)]['lnl']
                   for s in SEEDS])
    ah = np.array([REC[(lawtag, s)]['a'] for s in SEEDS])
    n_int = sum(REC[(lawtag, s)]['interior'] for s in SEEDS)
    md, se = float(ds.mean()), float(ds.std(ddof=1)/math.sqrt(6))
    if n_int < 6:
        v = "SHAPE REJECTION (edge; invalid read)"
    elif md > 0 and md > 2*se: v = "PREFERRED over the vector composition"
    elif md >= -3: v = "VIABLE (scalar-EFE stands)"
    elif md <= -5: v = "EXCLUDED (the data demand the vector composition)"
    else: v = "UNRESOLVED (between bands)"
    verdicts[lawtag] = (md, se, v)
    L.append(f"{lawtag} vs {ref} (vector): per-seed " +
             " ".join(f"{x:+.2f}" for x in ds))
    L.append(f"  mean {md:+.2f} +- {se:.2f} SE; better in "
             f"{int((ds > 0).sum())}/6; interior {n_int}/6; "
             f"a_hat {ah.mean():.3f} +- {ah.std(ddof=1)/math.sqrt(6):.3f}")
    L.append(f"  -> {v}")
    L.append("")

# a0 rows (6G kappa machinery on the scalar tables)
def load_tab2(path):
    t = np.load(path)
    yy, b = t[0][::-1], t[1][::-1]
    return np.log(yy), np.log(np.maximum(b-1.0, 1e-12))
tabs = {}
for fam in ('sbe', 'samb'):
    for e, tag in ((1.0, 'g1p0'), (1.2, 'g1p2'), (1.4, 'g1p4')):
        tabs[(fam, e)] = load_tab2(f'data/efe_boost_{fam}_{tag}.npy')
def lnB1(fam, e, lny):
    lo_, hi_ = (1.0, 1.2) if e <= 1.2 else (1.2, 1.4)
    w = (e-lo_)/(hi_-lo_)
    return (1-w)*np.interp(lny, *tabs[(fam, lo_)]) \
        + w*np.interp(lny, *tabs[(fam, hi_)])
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
tp = 2.998e8*(67.4*1e3/3.0857e22)/(2*math.pi)
for fam in ('sbe', 'samb'):
    ah = np.array([REC[(fam, s)]['a'] for s in SEEDS])
    m = float(ah.mean())
    stot = math.hypot(float(ah.std(ddof=1)/math.sqrt(6)), 0.11)
    kf = float(np.mean(kappa(fam, lnyd)))
    a0b = 1.2e-10*m**(1.0/kf)
    sa0b = a0b*stot/(m*kf)
    pull = (a0b - tp)/math.hypot(sa0b, tp*0.5/67.4)
    L.append(f"{fam}: kappa {kf:+.3f}, alpha {m:.3f}+-{stot:.3f} -> "
             f"a0 = {a0b*1e10:.2f}+-{sa0b*1e10:.2f}  pull {pull:+.1f} sigma"
             f"  [pre-reg PASS <= +2.5; vector rows: BE +1.9, AMB +1.6]")

# solar consequence (analytic; pre-stated in the docstring)
G_GAL = 1.9e-10          # m s^-2 at the solar circle
R_GAL = 2.55e20          # m (8.27 kpc)
Q_TIDE = G_GAL/R_GAL
L.append("")
L.append(f"solar consequence of the scalar class: spherical phantom shell"
         f" => shell theorem => interior quadrupole = the true galactic"
         f" tide ~ {Q_TIDE:.1e} s^-2 = {9e-27/Q_TIDE:.0f}x BELOW the"
         f" Cassini cap (vector class: 3.9e-26 = 4.3x ABOVE, 4K).")
ok_any = any(v[2].startswith(('VIABLE', 'PREFERRED'))
             for v in verdicts.values())
L.append("")
L.append("CONSEQUENCE (pre-stated map): " +
         ("the scalar-EFE class STANDS on all measured data - the "
          "quadrupole problem has a concrete non-MI resolution; the "
          "vector composition is no longer data-forced."
          if ok_any else
          "the binaries DEMAND the vector composition - Saturn stands; "
          "logged honestly."))

out = "\n".join(L)
print("\n" + out)
with open('data/stage6w_verdict.txt', 'w') as f:
    f.write(out + "\n")
print("\nSTAGE 6W done -> data/stage6w_verdict.txt")
