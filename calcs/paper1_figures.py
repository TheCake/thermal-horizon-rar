"""PAPER 1 FIGURES (draft 0.1 -> 0.2).  Every plotted number is either
computed from the catalog in this script, recomputed from archived npz
tables, or parsed from a committed stage output; the provenance dump
(data/paper1_figs.txt) lists all of them.  In-script gates regress the
computed values against the record before any figure is written.
Outputs: papers/figs/fig{1..5}*.png/.pdf
"""
import math, os, re, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import logsumexp

os.makedirs('papers/figs', exist_ok=True)
L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

P("PAPER 1 FIGURES - provenance-gated build")

# ---------------- catalog loader (verbatim 9F construction) ---------
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
sig_ok = sig_c[ok]
vc_ok = 0.9417*np.sqrt(Mtot_d/s_d)
P(f"catalog: N = {int(ok.sum())} pairs (expect 14071)")
assert int(ok.sum()) == 14071

# ---------------- Table 1: the sequential cut ladder ----------------
cuts = [
    ("catalog rows", np.ones(len(sep), dtype=bool)),
    ("chance-alignment R < 0.01", Rch < 0.01),
    ("both parallaxes > 5 mas", (plx1 > 5) & (plx2 > 5)),
    ("parallax S/N > 20 (both)",
     (plx1/np.maximum(eplx1, 1e-6) > 20)
     & (plx2/np.maximum(eplx2, 1e-6) > 20)),
    ("parallax consistency < 3 sigma",
     np.abs(plx1-plx2) < 3*np.hypot(eplx1, eplx2)),
    ("separation 200 AU - 50 kAU", (sep > 200) & (sep < 50000)),
    ("main sequence (both, 2.6 < M_G < 14.2)",
     (MG1 > 2.6) & (MG1 < 14.2) & (MG2 > 2.6) & (MG2 < 14.2)),
    ("velocity precision < 0.03 km/s", sigv < 0.03),
]
mrun = np.ones(len(sep), dtype=bool)
P("cut ladder (sequential):")
for nm, m in cuts:
    mrun &= m
    P(f"  {nm}: {int(mrun.sum())}")
assert int(mrun.sum()) == 14071

SBINS = [(0.2,2),(2,6),(6,20),(20,50)]
VE = np.logspace(np.log10(0.02), np.log10(6.0), 21)
GE = np.linspace(0, 90, 7)
NV, NG = 20, 6
vcen = np.sqrt(VE[:-1]*VE[1:]); gcen = 0.5*(GE[:-1]+GE[1:])

# ---------------- anchor statistic + gate ---------------------------
mask_n = (s_d >= 0.2) & (s_d < 2)
mask_w = (s_d >= 6) & (s_d < 30)
B_all = float(np.median(vt_d[mask_w])/np.median(vt_d[mask_n]))
rngb = np.random.default_rng(11)
reps = []
iw = np.where(mask_w)[0]; inn = np.where(mask_n)[0]
for _ in range(2000):
    rw = vt_d[rngb.choice(iw, len(iw))]
    rn = vt_d[rngb.choice(inn, len(inn))]
    reps.append(np.median(rw)/np.median(rn))
ci = np.percentile(reps, [16, 84])
g1 = abs(B_all - 1.0779) <= 0.002
P(f"G-F1 anchor regression: B(all) = {B_all:.4f} vs 9E printed 1.0779 "
  f"-> {'PASS' if g1 else 'FAIL'}; 68% CI [{ci[0]:.3f}, {ci[1]:.3f}]")

# ---------------- census selection + gate ---------------------------
perp = (s_d >= 6) & (gam_d >= 75)
vtp = vt_d[perp]
n_band = int(np.sum((vtp >= 1.414) & (vtp < 1.67)))
n_over = int(np.sum((vtp >= 1.67) & (vtp < 2.2)))
rows = list(csv.DictReader(open('data/ceiling_pairs.csv')))
c_band = sum(1 for r in rows if r['band_corr'] == 'band')
c_over = sum(1 for r in rows if r['band_corr'] == 'above'
             and 1.67 <= float(r['vt_corr']) < 2.2)
g2 = (n_band == 9 and n_over == 2 and c_band == 9 and c_over == 2)
P(f"G-F2 census regression: loader band/over = {n_band}/{n_over}; "
  f"CSV = {c_band}/{c_over}; record = 9/2 -> {'PASS' if g2 else 'FAIL'}")

# ---------------- parse committed stage outputs ---------------------
t7kb = open('data/stage7kb_census.txt').read()
def kb(pat):
    m = re.search(pat, t7kb)
    return (float(m.group(1)), float(m.group(2)))
mu_nb = kb(r"newton-best  \(fpm=3\.0, sq=0\.2\): mu_band = ([\d.]+), "
           r"P\(>=9\) = [\dEe.+-]+ -> NULL-BROKEN  \[overshoot mu = "
           r"([\d.]+)")
mu_bs = kb(r"POST-HOC boost simple sq=0 fpm=1\.5: mu_band = ([\d.]+) "
           r"\(obs 9\), overshoot mu = ([\d.]+)")
P(f"parsed 7K-b: newton-best mu = {mu_nb}; boost-edge(sq0) mu = {mu_bs}")

t7ka = open('data/stage7ka_median.txt').read()
Rn = [float(x) for x in re.findall(
    r"newton       alpha=0\.0 .*?: R = [\d.]+ nat / ([\d.]+) reweighted",
    t7ka)]
Rs = [float(x) for x in re.findall(
    r"prof-simple  alpha=0\.5 .*?: R = [\d.]+ nat / ([\d.]+) reweighted",
    t7ka)]
Rb = [float(x) for x in re.findall(
    r"prof-BE      alpha=1\.0 .*?: R = [\d.]+ nat / ([\d.]+) reweighted",
    t7ka)]
cellm = re.search(r"newton       alpha=0\.0 eta=([\d.]+) wr=([\d.]+) "
                  r"fcomp=([\d.]+) fpm=([\d.]+) kw=([\d.]+) sq=([\d.]+)",
                  t7ka)
ETA0, WR0, FC_C, FPM_C, KW_C, SQ_C = (float(cellm.group(i))
                                      for i in range(1, 7))
P(f"parsed 7K-a: forward newton R = {Rn}; boost R = {Rs}+{Rb}; "
  f"cell eta={ETA0} wr={WR0} fcomp={FC_C} fpm={FPM_C} kw={KW_C} "
  f"sq={SQ_C}")

t7l = open('data/stage7l_step.txt').read()
mstep = re.search(r"data step .* = ([\d.]+) \(68% CI ([\d.]+)-([\d.]+)",
                  t7l)
step_d = (float(mstep.group(1)), float(mstep.group(2)),
          float(mstep.group(3)))
step_nb = [float(x) for x in re.findall(
    r"newton-best  alpha=0\.0: step = ([\d.]+)", t7l)]
step_bs = [float(x) for x in re.findall(
    r"boost-simple alpha=0\.5: step = ([\d.]+)", t7l)]
step_bb = [float(x) for x in re.findall(
    r"boost-BE     alpha=1\.0: step = ([\d.]+)", t7l)]
P(f"parsed 7L: data step = {step_d}; model steps nb/bs/bb = "
  f"{step_nb}/{step_bs}/{step_bb}")

t9j = open('data/stage9j_stdext.txt').read()
qa = re.findall(r"Q([1-4])-alone full: a_marg = ([\d.]+)", t9j)
dose = np.array([float(v) for _, v in qa]).reshape(4, 4)  # rows =
# (simple 31, BE 31, simple 101, BE 101) in file order; cols = Q1..Q4
P(f"parsed 9J dose rows x Q: {dose.tolist()}")

t9l = open('data/stage9l_fpmmeter.txt').read()
mq = re.findall(r"\[seed (\d+)\] Q([1-4]) narrow: .*E\[fpm\] = ([\d.]+)",
                t9l)
meter = {}
for sd, q, v in mq:
    meter.setdefault(int(sd), {})[int(q)] = float(v)
P(f"parsed 9L meter: {meter}")

# ---------------- Fig 3 recompute (9O combiner, verbatim) -----------
FCOMP_GRID = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
A_FULL = np.concatenate([np.arange(0, 0.61, 0.1),
                         np.array([0.7, 0.8, 0.9, 1.0, 1.1, 1.2])])
NC_STD = 6*4; NC_SHR = 2*3*2
pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
gband = pz['conv_band']/0.30
GS = np.linspace(float(gband[0]), float(gband[1]), 25)
fhmin, fhmax = fg[sup].min(), fg[sup].max()
def lnpi_at(gi_list):
    out = np.full(len(FCOMP_GRID), -1e9)
    for gi in gi_list:
        fh_eq = FCOMP_GRID/gi
        m = (fh_eq >= fhmin) & (fh_eq <= fhmax)
        cand = np.full(len(FCOMP_GRID), -1e9)
        cand[m] = np.interp(fh_eq[m], fg, lp)
        out = np.maximum(out, cand)
    return out
VARIANTS = [('measured prior (envelope)', lnpi_at(GS)),
            ('low-edge conversion', lnpi_at([float(gband[0])])),
            ('mid conversion', lnpi_at([0.5*(float(gband[0])+float(gband[1]))])),
            ('high-edge conversion', lnpi_at([float(gband[1])])),
            ('flat (no measurement)', np.zeros(len(FCOMP_GRID)))]
def lse_cells(T, axes, ncells):
    return logsumexp(T, axis=axes) - math.log(ncells)
def marg13(STDF, LNPIv, strata):
    lnZ = np.zeros(13)
    for ai in range(13):
        S = LNPIv[:, None, None, None].copy()
        for qi in strata:
            S = S + lse_cells(STDF[ai, qi], (2, 4), NC_STD)
        lnZ[ai] = logsumexp(S) - math.log(NC_SHR)
    return lnZ
posts = {}
amarg_tab = {}
for seed in (31, 101):
    for law in ('simple', 'BE'):
        z9i = np.load(f'data/stage9i_tables_{seed}_{law}.npz')
        z9j = np.load(f'data/stage9j_tables_{seed}_{law}.npz')
        STDF = np.concatenate([z9i['STD'], z9j['STDX']], axis=0)
        for nm, LNPIv in VARIANTS:
            zd = marg13(STDF, LNPIv, [0, 1, 2])
            w = np.exp(zd - logsumexp(zd))
            posts[(law, seed, nm)] = w
            amarg_tab[(law, seed, nm)] = float(np.sum(w*A_FULL))
t9o = open('data/stage9o_lnpiband.txt').read()
tgt = {}
for lw, sd, nm2, am in re.findall(
        r"\[(simple|BE) (\d+)\] (OPER|FLAT)\s*: a_marg\(DROP\) = ([\d.]+)",
        t9o):
    tgt[(lw, int(sd), nm2)] = float(am)
g3 = all(abs(amarg_tab[(lw, sd, VARIANTS[0][0])] - tgt[(lw, sd, 'OPER')])
         <= 0.002 and
         abs(amarg_tab[(lw, sd, VARIANTS[4][0])] - tgt[(lw, sd, 'FLAT')])
         <= 0.002 for lw in ('simple', 'BE') for sd in (31, 101))
P(f"G-F3 9O regression (OPER + FLAT a_marg, 8 values, 0.002) -> "
  f"{'PASS' if g3 else 'FAIL'}")

if not (g1 and g2 and g3):
    P("GATES FAILED - STOP, no figures written")
    with open('data/paper1_figs.txt', 'w') as f:
        f.write("\n".join(L_)+"\n")
    raise SystemExit(1)

# ================= FIGURES ==========================================
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9,
                     'figure.dpi': 150})

# ---- Fig 2: the dose curve + the meter -----------------------------
fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.2))
labels = ['simple, r31', 'occupation, r31', 'simple, r101',
          'occupation, r101']
mk = ['o', 's', '^', 'v']
for i in range(4):
    ax[0].plot([1, 2, 3, 4], dose[i], marker=mk[i], lw=1,
               label=labels[i])
ax[0].set_xticks([1, 2, 3, 4])
ax[0].set_xticklabels(['Q1\n(cleanest)', 'Q2', 'Q3', 'Q4\n(worst)'])
ax[0].set_xlabel('RUWE quartile')
ax[0].set_ylabel(r'marginalized boost amplitude $\alpha$')
ax[0].axhline(0.5, color='k', ls=':', lw=0.8)
ax[0].text(1.02, 0.51, 'excluded on clean strata', fontsize=7)
ax[0].legend(fontsize=6.5, loc='upper left')
ax[0].set_title('(a) fitted amplitude by astrometric quality')
for sd, m_ in ((31, 'o'), (101, 's')):
    ax[1].plot([1, 2, 3, 4], [meter[sd][q] for q in (1, 2, 3, 4)],
               marker=m_, lw=1, label=f'realization {sd}')
ax[1].axhline(2.1, color='k', ls=':', lw=0.8)
ax[1].text(1.02, 2.115, 'quoted lower bound (2.1x formal)', fontsize=7)
ax[1].set_xticks([1, 2, 3, 4])
ax[1].set_xticklabels(['Q1', 'Q2', 'Q3', 'Q4'])
ax[1].set_xlabel('RUWE quartile')
ax[1].set_ylabel(r'$E[f_{\rm pm}]$ (boost-free 0.2-2 kAU bin)')
ax[1].set_ylim(1.0, 3.0)
ax[1].legend(fontsize=6.5, loc='lower right')
ax[1].set_title('(b) honest noise inflation by quality')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/fig2_dosecurve.{ext}')
plt.close(fig)
P("fig2 written")

# ---- Fig 3: prior family -------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(8.2, 3.2))
for i, (lw_, sd) in enumerate([('simple', 31), ('BE', 31),
                               ('simple', 101), ('BE', 101)]):
    ax[0].plot(A_FULL, posts[(lw_, sd, VARIANTS[0][0])], marker=mk[i],
               ms=3, lw=1, label=labels[i])
ax[0].axvspan(0.5, 1.2, color='0.9')
ax[0].text(0.52, ax[0].get_ylim()[1]*0.9, 'excluded', fontsize=7)
ax[0].set_xlabel(r'$\alpha$')
ax[0].set_ylabel('posterior weight (clean strata)')
ax[0].legend(fontsize=6.5)
ax[0].set_title('(a) clean-strata posterior, measured prior')
xs = np.arange(5)
for i, (lw_, sd) in enumerate([('simple', 31), ('BE', 31),
                               ('simple', 101), ('BE', 101)]):
    ax[1].plot(xs + (i-1.5)*0.09,
               [amarg_tab[(lw_, sd, nm)] for nm, _ in VARIANTS],
               mk[i], ms=4, label=labels[i])
ax[1].axhline(0.5, color='k', ls=':', lw=0.8)
ax[1].set_xticks(xs)
ax[1].set_xticklabels([nm.replace(' ', '\n', 1) for nm, _ in VARIANTS],
                      fontsize=6.5)
ax[1].set_ylabel(r'marginalized $\alpha$ (clean strata)')
ax[1].set_title('(b) amplitude under companion-prior variants')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/fig3_priorfamily.{ext}')
plt.close(fig)
P("fig3 written")

# ---- Fig 4: the census column --------------------------------------
fig, ax = plt.subplots(figsize=(5.4, 3.4))
be = np.arange(1.0, 2.45, 0.06)
ax.hist(vtp[(vtp >= 1.0) & (vtp < 2.45)], bins=be, color='0.55',
        edgecolor='k', lw=0.4)
ax.axvspan(1.414, 1.67, color='tab:orange', alpha=0.18)
ax.axvspan(1.67, 2.2, color='tab:red', alpha=0.10)
ax.axvline(1.414, color='k', ls='--', lw=0.9)
ax.axvline(1.65, color='k', ls='-', lw=0.9)
ax.text(1.42, ax.get_ylim()[1]*0.92,
        f'forbidden band\nobs {n_band}', fontsize=7)
ax.text(1.70, ax.get_ylim()[1]*0.92,
        f'overshoot\nobs {n_over}', fontsize=7)
txt = (f"expected (band, overshoot):\n"
       f"  Newton, formal errors: (0.9, ~0)\n"
       f"  fitted Newton cell: ({mu_nb[0]:.0f}, {mu_nb[1]:.0f})\n"
       f"  boosted edge, formal: ({mu_bs[0]:.0f}, {mu_bs[1]:.0f})")
ax.text(0.985, 0.55, txt, transform=ax.transAxes, fontsize=6.8,
        ha='right', va='top',
        bbox=dict(fc='white', ec='0.6', lw=0.5))
ax.set_xlabel(r'$\tilde{v}$ (perpendicular column: '
              r'$s \geq 6$ kAU, $\gamma \geq 75^\circ$)')
ax.set_ylabel('pairs')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/fig4_census.{ext}')
plt.close(fig)
P("fig4 written (Newton-formal 0.9 annotation sourced from the Stage "
  "4J leakage null, PAPER.md sec 7.2)")

# ---- Fig 5: the median forest --------------------------------------
fig, ax = plt.subplots(figsize=(6.4, 3.4))
rows5 = [
    ('observed anchor ratio (this work)', B_all, ci[0], ci[1], 'k'),
    ('forward model: fitted Newton cell', np.mean(Rn), min(Rn),
     max(Rn), 'tab:blue'),
    ('forward model: boosted cells', np.mean(Rs+Rb), min(Rs+Rb),
     max(Rs+Rb), 'tab:red'),
    ('Cookson selection: observed step', step_d[0], step_d[1],
     step_d[2], 'k'),
    ('Cookson selection: fitted Newton cell', np.mean(step_nb),
     min(step_nb), max(step_nb), 'tab:blue'),
    ('Cookson selection: boosted cells',
     np.mean(step_bs+step_bb), min(step_bs+step_bb),
     max(step_bs+step_bb), 'tab:red'),
]
ys = np.arange(len(rows5))[::-1]
for y, (lab, v, lo, hi, c) in zip(ys, rows5):
    ax.errorbar([v], [y], xerr=[[v-lo], [hi-v]], fmt='o', ms=4,
                color=c, capsize=2, lw=1)
    ax.text(0.865, y+0.18, lab, fontsize=7)
ax.axvline(1.0, color='0.6', ls=':', lw=0.8)
ax.set_yticks([])
ax.set_xlabel('median velocity ratio (wide / anchor)')
ax.set_xlim(0.86, 1.22)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/fig5_medians.{ext}')
plt.close(fig)
P("fig5 written")

# ---- Fig 1: data vs the two 7K-a worlds (GPU) ----------------------
P("fig1: building model overlays (GPU)...")
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
SC2 = (np.array([0.63, 3.46, 11.0, 31.6])/31.6)**2
FLY = np.zeros((NV, NG))
for i in range(NV):
    for j in range(NG):
        FLY[i,j] = (1.0 if (0.5 <= vcen[i] <= 3.0) else 0.0) \
                   * max(30.0-gcen[j], 0.0)
FLY /= max(FLY.sum(), 1e-12)
UNI = np.ones((NV, NG))/(NV*NG)
def build_stratum(pm):
    D2, PLs, UB, FB = [], [], [], []
    for b in SBINS:
        m = (s_d>=b[0])&(s_d<b[1])&pm
        h,_,_ = np.histogram2d(np.clip(vt_d[m],0.021,5.9), gam_d[m],
                               bins=[VE, GE])
        D2.append(h.astype(float))
        PLs.append(sig_ok[m])
        cutp = 2.978/np.sqrt(s_d[m]) + 2.8284*sig_ok[m]
        acc = np.array([(vcen[i]*vc_ok[m] <= cutp).mean()
                        if m.sum() else 0.0 for i in range(NV)])
        for tpl, store in ((UNI, UB), (FLY, FB)):
            t = tpl*acc[:,None]
            store.append(t/max(t.sum(), 1e-12))
    ND = [int(h.sum()) for h in D2]
    return D2, PLs, UB, FB, ND
data_2d, noise_pool, UNI_B, FLY_B, NDATA = build_stratum(
    np.ones(len(s_d), dtype=bool))
N = 500_000
FFLY_GRID = np.array([0.05, 0.10]); FC0 = 0.10
KW_GRID = np.array([0.7, 1.0, 1.4])
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
    p['pick'] = [rng.integers(0, max(len(noise_pool[bi]),1), N)
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
def pp_at_cell(p, prj, fcm, ff, fpm, kwv, sqv):
    smag, vpar, vper = prj
    s_kau = smag/1e3
    pps = []
    for bi, b in enumerate(SBINS):
        idx = np.where((s_kau>=b[0])&(s_kau<b[1]))[0]
        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
        sg0 = noise_pool[bi][p['pick'][bi][idx] % len(noise_pool[bi])]/4.74047
        g1_i, g2_i = p['gn1'][idx], p['gn2'][idx]
        sk_i = s_kau[idx]
        gk_full = p['gs'][idx]
        cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
        mh_tot = np.zeros(len(idx))
        for k in (1, 2):
            c = p['comp'][k]
            act = c['uc'][idx] < fcm
            mh_tot += act*c['mh'][idx]
            cvp += act*c['w'][idx]*c['wd'][idx,0]
            cvq += act*c['w'][idx]*c['wd'][idx,1]
        boost = np.sqrt(1+mh_tot/p['M_s'][idx])
        vp_a = vpar[idx] + kwv*cvp
        vq_a = vper[idx] + kwv*cvq
        vp_n = vp_a*boost + g1_i*sg0*fpm
        vq_n = vq_a*boost + g2_i*sg0*fpm
        vmag = np.hypot(vp_n, vq_n)
        keep = vmag*4.74047 <= (2.978/np.sqrt(sk_i)
                                + 2.8284*sg0*4.74047)
        vtn = (vmag/vc)[keep]
        gmn = np.degrees(np.arccos(np.clip(
            np.abs(vp_n[keep])/np.maximum(vmag[keep], 1e-12), 0, 1)))
        gk = gk_full[keep]
        vts = vtn*np.exp(sqv*gk)
        h,_,_ = np.histogram2d(np.clip(vts,0.021,5.9), gmn,
                               bins=[VE, GE])
        p0 = np.maximum(h/max(h.sum(),1), 1e-5)
        p0 /= p0.sum()
        wch = min(FC0*SC2[bi], 0.5)
        wfl = min(ff*SC2[bi], 0.5)
        wtot = min(wch+wfl, 0.6)
        mixc = (wch*UNI_B[bi] + wfl*FLY_B[bi])/(wch+wfl)
        pp = (1-wtot)*p0 + wtot*mixc
        pps.append(pp)
    return pps
pf = build_pop(31)
e_f = e_of_x(pf, ETA0, WR0)
pps = {}
for lab, al in (('newton', 0.0), ('boost', 0.5)):
    tab_a = 1.0 + al*(TAB_S-1.0)
    vp_f = vp_c(pf, e_f, tab_a)
    o_f = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'], pf['M_s'],
              pf['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab_a,
              lny0=LNY0, dlny=DLNY, vp=vp_f)
    prj = project(pf, o_f)
    pps[lab] = pp_at_cell(pf, prj, FC_C, 0.05, FPM_C, KW_C, SQ_C)
    P(f"  {lab} overlay done")
fig, ax = plt.subplots(2, 4, figsize=(10.2, 4.6))
for bi, b in enumerate(SBINS):
    nd = data_2d[bi].sum()
    dv = data_2d[bi].sum(axis=1)
    dg = data_2d[bi].sum(axis=0)
    ax[0, bi].step(vcen, dv, where='mid', color='k', lw=1,
                   label='data')
    ax[1, bi].step(gcen, dg, where='mid', color='k', lw=1)
    for lab, c in (('newton', 'tab:blue'), ('boost', 'tab:red')):
        mv = pps[lab][bi].sum(axis=1)*nd
        mg = pps[lab][bi].sum(axis=0)*nd
        ax[0, bi].plot(vcen, mv, color=c, lw=1,
                       label=('fitted Newton cell' if lab == 'newton'
                              else r'boosted cell ($\alpha=0.5$)'))
        ax[1, bi].plot(gcen, mg, color=c, lw=1)
    ax[0, bi].set_xscale('log')
    ax[0, bi].set_title(f'{b[0]}-{b[1]} kAU')
    ax[0, bi].set_xlabel(r'$\tilde{v}$')
    ax[1, bi].set_xlabel(r'$\gamma$ (deg)')
    if bi == 0:
        ax[0, bi].set_ylabel('pairs')
        ax[1, bi].set_ylabel('pairs')
        ax[0, bi].legend(fontsize=6)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/fig1_vgamma.{ext}')
plt.close(fig)
P("fig1 written (overlays at the 7K-a printed cells, realization 31, "
  "simple family; ffly = 0.05 disclosed)")

# ---- Fig 6: the phantom veto ---------------------------------------
t7jg = open('data/stage7jg_read.txt').read()
jg = {}
for lw_, sd, ch, cfg, am, dn in re.findall(
        r"\[(simple|BE) (\d+) (2D|vt) (sqfree|sq0)\] "
        r"a_marg=([\d.]+) dN=\+?([\d.]+)", t7jg):
    jg[(lw_, int(sd), ch, cfg)] = (float(am), float(dn))
assert len(jg) == 16, jg
P(f"parsed 7J-g: 16 rows; vt-sq0 amplitudes = "
  + "/".join(f"{jg[(l_, s_, 'vt', 'sq0')][0]:.2f}"
             for l_ in ('simple', 'BE') for s_ in (31, 101)))
fig, ax = plt.subplots(figsize=(5.6, 3.4))
configs = [('vt', 'sq0', 'velocity-only,\nno width channel'),
           ('2D', 'sq0', 'joint 2D,\nno width channel'),
           ('vt', 'sqfree', 'velocity-only,\nwidth channel on'),
           ('2D', 'sqfree', 'joint 2D,\nwidth channel on')]
for i, (lw_, sd) in enumerate([('simple', 31), ('BE', 31),
                               ('simple', 101), ('BE', 101)]):
    ys = [jg[(lw_, sd, ch, cfg)][0] for ch, cfg, _ in configs]
    ax.plot(np.arange(4) + (i-1.5)*0.07, ys, mk[i], ms=5,
            label=labels[i])
ax.annotate('', xy=(1.1, 0.24), xytext=(0.4, 0.48),
            arrowprops=dict(arrowstyle='->', lw=1, color='0.3'))
ax.text(0.42, 0.30, 'the direction\ndata veto', fontsize=7,
        color='0.25')
ax.set_xticks(np.arange(4))
ax.set_xticklabels([c[2] for c in configs], fontsize=7)
ax.set_ylabel(r'fitted amplitude $\alpha$')
ax.axhline(0, color='0.7', lw=0.6)
ax.legend(fontsize=6.5, loc='upper left')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/fig6_phantom.{ext}')
plt.close(fig)
P("fig6 written")

with open('data/paper1_figs.txt', 'w') as f:
    f.write("\n".join(L_)+"\n")
print("\nsaved: data/paper1_figs.txt")
