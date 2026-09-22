"""PROPOSAL F FIGURES (the Hubble-meter paper).  Every plotted number is
parsed from a committed stage output, loaded from the committed map CSV,
or recomputed here from declared conventions; in-script gates (GF-1..GF-5)
regress every load-bearing value against the record before any figure is
written.  Provenance dump: data/paperf_figs.txt.
Outputs: papers/figs/figf{1..5}.png/.pdf

Sources (all committed):
  data/stage10s_h0meter_gates.txt  -- the distance-method census (G2 line)
  data/stage10t_skyread.txt        -- leg-A first light (R46-corrected run)
  data/stage10t_verdict.txt        -- validity-domain clause
  data/stage10w_verdict.txt        -- map totals / carriers (R49-adopted)
  data/stage10w_map.csv            -- the per-galaxy form-preference map
  data/stage10h_addendum.txt       -- z~1 intercept comparison (R29-verified)
  data/lit0818_a0z.py conventions  -- lock-curve cosmology + MIGHTEE rows
    (verbatim constants restated below; the a1_lock gate recomputes them)
External reference constants (figure landmarks only, cited in captions):
  Planck 2018 H0 = 67.4 +- 0.5; SH0ES (Riess et al. 2022) 73.0 +- 1.0.
NEVER-QUOTE guardrails honored: no leg-B number, no flat-engine sigma as a
statistic (envelope only), no cut-flow a0, no N~25 forecast.
"""
import os, re, csv, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
os.makedirs('papers/figs', exist_ok=True)

L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

P("PROPOSAL F FIGURES - provenance-gated build")

def ftext(path):
    return open(path, encoding='utf-8').read()

def grep1(txt, pat, path='?'):
    m = re.search(pat, txt)
    assert m, f"parse failed in {path}: {pat}"
    return m

# ================= shared style =====================================
plt.rcParams.update({
    'font.size': 9.5, 'axes.titlesize': 10, 'axes.labelsize': 9.5,
    'legend.fontsize': 8.5, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'figure.dpi': 200, 'savefig.dpi': 200, 'axes.linewidth': 0.8,
})
C_ANCH = '#2060a8'   # anchored leg
C_FLOW = '#c8641e'   # flow leg
C_BAND = '#9ec7e8'   # form band fill
C_GRAY = '#666666'

def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(f'papers/figs/{name}.{ext}', bbox_inches='tight')
    plt.close(fig)
    P(f"wrote papers/figs/{name}.png/.pdf")

# ================= parse the record =================================
t10s = ftext('data/stage10s_h0meter_gates.txt')
t10t = ftext('data/stage10t_skyread.txt')
t10tv = ftext('data/stage10t_verdict.txt')
t10w = ftext('data/stage10w_verdict.txt')
t10h = ftext('data/stage10h_addendum.txt')

# --- census (10S G2) ---
m = grep1(t10s, r"counts flow/TRGB/Cep/UMa/SN = (\d+)/(\d+)/(\d+)/(\d+)/(\d+)",
          'stage10s_h0meter_gates')
n_flow0, n_trgb, n_cep, n_uma, n_sn = (int(g) for g in m.groups())
n_anch0 = n_trgb + n_cep + n_sn
P(f"census (10S G2): flow {n_flow0} / anchored {n_anch0} "
  f"(TRGB {n_trgb} + Cep {n_cep} + SN {n_sn}) / UMa {n_uma}")

# --- leg architecture (10T) ---
m = grep1(t10t, r"primary leg A = (\d+) \((\d+) anchored \+ (\d+) UMa \+ "
                r"(\d+) reclassified\)", 'stage10t_skyread')
nA, nA_anch, nA_uma, nA_re = (int(g) for g in m.groups())
n_flowB = int(grep1(t10t, r"flow 82 -> (\d+)", 'stage10t_skyread').group(1))
P(f"legs (10T): leg A = {nA} ({nA_anch}+{nA_uma}+{nA_re}); flow leg = {n_flowB}")

# --- first light (10T sky read, R46-corrected) ---
m = grep1(t10t, r"PRIMARY hier BE \(DEEP, R46-C2\): a0 = ([\d.e+-]+), .*?"
                r"-> H0_A = ([\d.]+)", 'stage10t_skyread')
a0_A, H0_A = float(m.group(1)), float(m.group(2))
m = grep1(t10t, r"hier family \(deep\): BE ([\d.e+-]+) \(H0 ([\d.]+)\), "
                r"p065 ([\d.e+-]+) \(H0 ([\d.]+)\), gm ([\d.e+-]+) \(H0 ([\d.]+)\), "
                r"boot ([\d.e+-]+) \(H0 ([\d.]+)\)", 'stage10t_skyread')
fam = {'BE': (float(m.group(1)), float(m.group(2))),
       'p065': (float(m.group(3)), float(m.group(4))),
       'gm': (float(m.group(5)), float(m.group(6))),
       'boot': (float(m.group(7)), float(m.group(8)))}
m = grep1(t10t, r"sigma_A \(hier/PRIMARY\): ([\d.]+); percentiles 16/50/84 = "
                r"\[([\d.]+), ([\d.]+), ([\d.]+)\]", 'stage10t_skyread')
sigA = float(m.group(1))
pc16, pc50, pc84 = float(m.group(2)), float(m.group(3)), float(m.group(4))
m = grep1(t10t, r"flat co-read \(the pre-registered treatment variant\): "
                r"a0 = ([\d.e+-]+) \(H0 ([\d.]+)\)", 'stage10t_skyread')
H0_flat = float(m.group(2))
sig_flat = float(grep1(t10t, r"sigma_A \(flat engine, co-quote\): ([\d.]+);",
                       'stage10t_skyread').group(1))
H0_var2 = float(grep1(t10t, r"A-var2 .*\(H0 ([\d.]+)\)", 'stage10t_skyread').group(1))
H0_var4 = float(grep1(t10t, r"A-var4 .*\(H0 ([\d.]+)\)", 'stage10t_skyread').group(1))
assert 'validity to readings <= ~75' in t10tv, 'validity clause missing'
VAL_LO, VAL_HI = 62.0, 75.0   # R46-A1 validity domain (injections 62-73 clean;
                              # measured -3.0% one-sided bias at truth 85)
P(f"first light: H0_A = {H0_A} +- {sigA} (pct {pc16}/{pc50}/{pc84}); "
  f"family {fam['BE'][1]}-{fam['p065'][1]}; variants {H0_var4}/{H0_var2}; "
  f"flat {H0_flat} env {sig_flat}")

# gate GF-2 (first light consistency)
band_lo, band_hi = fam['BE'][1], fam['p065'][1]
assert (band_lo, band_hi) == (65.4, 70.4), (band_lo, band_hi)
assert abs(H0_A - 65.36) < 0.01 and abs(sigA - 4.99) < 0.01
assert H0_var4 == 64.8 and H0_var2 == 65.3
assert min(fam[k][1] for k in fam) == band_lo
assert max(fam[k][1] for k in fam) == band_hi
# a0 <-> H0 conversion is linear in the lock: check each family member
for k, (a0k, h0k) in fam.items():
    assert abs(a0k / a0_A * H0_A - h0k) < 0.06, (k, a0k, h0k)
P("GF-2 PASS: first-light values parsed + lock-linearity closed (<0.06)")

# --- the map (10W CSV + verdict cross-check) ---
rows = list(csv.DictReader(open('data/stage10w_map.csv', encoding='utf-8')))
assert len(rows) == 149, len(rows)
name = np.array([r['name'] for r in rows])
leg = np.array([r['leg'] for r in rows])
cov = np.array([float(r['cov']) for r in rows])
d_boot = np.array([float(r['d_boot']) for r in rows])
sd2 = np.array([float(r['sd2_share']) for r in rows])
pend = np.array([r['flag'] == 'census-pending' for r in rows])
mA, mF = leg == 'anch', leg == 'flow'
assert int(mA.sum()) == 78 and int(mF.sum()) == 71, (mA.sum(), mF.sum())

m = grep1(t10w, r"UGC03580 =\s+(-[\d.]+) of the flow total (-[\d.]+) "
                r"\((\d+)% of block variance\); anchored\s+lead broad "
                r"\(top carrier IC2574 \+([\d.]+) = (\d+)%\)", 'stage10w_verdict')
v_ugc, v_flowtot = float(m.group(1)), float(m.group(2))
v_ugcshare, v_ic, v_icshare = int(m.group(3)), float(m.group(4)), int(m.group(5))
anch_tot_v = float(grep1(t10w, r"\+([\d.]+) becomes \+63.5, \+66.0 and \+61.8",
                         'stage10w_verdict').group(1))
flow_tot = float(d_boot[mF].sum())
anch_tot = float(d_boot[mA].sum())
i_ugc = int(np.where(name == 'UGC03580')[0][0])
i_ic = int(np.where(name == 'IC2574')[0][0])
assert abs(flow_tot - v_flowtot) < 0.05, (flow_tot, v_flowtot)
assert abs(anch_tot - anch_tot_v) < 0.05, (anch_tot, anch_tot_v)
assert abs(d_boot[i_ugc] - v_ugc) < 0.05, d_boot[i_ugc]
assert abs(d_boot[i_ic] - v_ic) < 0.05, d_boot[i_ic]
assert abs(round(100 * sd2[i_ugc]) - v_ugcshare) <= 1, sd2[i_ugc]
assert abs(round(100 * sd2[i_ic]) - v_icshare) <= 1, sd2[i_ic]
assert i_ic == int(np.where(mA)[0][np.argmax(d_boot[mA])])
assert set(name[pend]) == {'D564-8', 'D631-7'}
P(f"GF-3 PASS: map closed vs verdict (flow {flow_tot:+.1f} = {v_flowtot:+.1f}; "
  f"anch {anch_tot:+.1f} = +{anch_tot_v:.1f}; UGC03580 {d_boot[i_ugc]:+.1f} "
  f"[{100*sd2[i_ugc]:.0f}%]; IC2574 {d_boot[i_ic]:+.1f} [{100*sd2[i_ic]:.0f}%])")

# --- three worlds (PREDICTIONS F1, scored in the SF annotation) ---
W_L, W_N = 67.0, 73.0            # meter predictions: worlds L and N
W_P_LO, W_P_HI = 75.0, 79.0      # world P registered range (gamma 1.3..2)
sL = abs(H0_A - W_L) / sigA
sN = abs(H0_A - W_N) / sigA
sP_lo = abs(H0_A - W_P_LO) / sigA
sP_hi = abs(H0_A - W_P_HI) / sigA
assert round(sL, 1) == 0.3 and round(sN, 1) == 1.5, (sL, sN)
assert round(sP_lo, 1) == 1.9 and round(sP_hi, 1) == 2.7, (sP_lo, sP_hi)
# the measured peg gearing (R51 M1 fix): parse the measurement record
tgear = ftext('data/round51_gearing.txt')
GAMMA = float(grep1(tgear, r"two-sided = ([\d.]+)", 'round51_gearing')
              .group(1))
GAMMA_FLAT = float(grep1(tgear, r"flat co-read .*gamma = ([\d.]+)",
                         'round51_gearing').group(1))
assert abs(GAMMA - 1.59) < 0.01, GAMMA
W_P_MEAS = W_L * (W_N / W_L) ** GAMMA     # 67 x (73/67)^gamma
sP_meas = abs(H0_A - W_P_MEAS) / sigA
assert abs(W_P_MEAS - 76.8) < 0.1 and abs(sP_meas - 2.29) < 0.02, \
    (W_P_MEAS, sP_meas)
assert W_P_LO <= W_P_MEAS <= W_P_HI
P(f"GF-4 PASS: world scoring reproduced (L {sL:.2f} / N {sN:.2f} / "
  f"P {sP_lo:.2f}-{sP_hi:.2f} sigma at sigma = {sigA}); measured "
  f"gearing {GAMMA:.3f} (flat {GAMMA_FLAT:.2f}) -> world P "
  f"{W_P_MEAS:.1f} at {sP_meas:.2f} sigma")

# --- lock-curve conventions (verbatim data/lit0818_a0z.py) ---
H0_KM_S_MPC, OMEGA_M = 70.0, 0.3
OMEGA_L = 1.0 - OMEGA_M
C_M_S = 2.99792458e8
MPC_M = 3.0856775814913673e22
H0_SI = H0_KM_S_MPC * 1.0e3 / MPC_M
UNIT = 1.0e-10
ZMAX_MIG = 0.09

def Efac(z):
    return math.sqrt(OMEGA_M * (1.0 + z) ** 3 + OMEGA_L)

def a0_lock(z, h0=H0_KM_S_MPC):
    return C_M_S * (h0 * 1.0e3 / MPC_M) * Efac(z) / (2.0 * math.pi)

a0_0 = a0_lock(0.0) / UNIT
a1_lock = (a0_lock(ZMAX_MIG) - a0_lock(0.0)) / ZMAX_MIG / UNIT
assert abs(a0_0 - 1.0826) < 0.001, a0_0
assert abs(a1_lock - 0.52) < 0.01, a1_lock       # the banked +0.52e-10 per z
# MIGHTEE-HI/LADUMA rows (verbatim from data/lit0818_a0z.py, arXiv:2608.03576)
MIG_PURE = (1.54, 0.11, -1.60, 2.33)   # a0, e_a0, a1, e_a1  (sample alone)
MIG_ANCH = (1.15, 0.02, 5.23, 1.05)    # + SPARC as z~0 anchor
# z~1 intercept comparison (MUSE-DARK III; data/stage10h_addendum.txt,
# their Planck-2015 cosmology; lock intercept there = 1.047)
LOCK_Z1 = 1.047
ciocan = {}
for tag in ('DM', 'MOND', 'pergal'):
    m = grep1(t10h, tag + r": \|1\.047 - ([\d.]+)\| / ([\d.]+) = ([\d.]+) sigma",
              'stage10h_addendum')
    ciocan[tag] = (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    d = abs(LOCK_Z1 - ciocan[tag][0]) / ciocan[tag][1]
    assert abs(d - ciocan[tag][2]) < 0.015, (tag, d)
P(f"GF-5 PASS: lock conventions recomputed (a0(0) = {a0_0:.4f}, "
  f"a1_lock = {a1_lock:+.2f}); z~1 rows reproduce "
  + "/".join(f"{ciocan[t][2]:.2f}" for t in ('DM', 'MOND', 'pergal')) + " sigma")

# external landmarks (captions carry the citations)
PLANCK, E_PLANCK = 67.4, 0.5
SH0ES, E_SH0ES = 73.0, 1.0

# --- GF-6: the conversion-branch arithmetic (paper section 8/9) ------
# Asymptotic de Sitter branch: H_inf = H0 sqrt(Omega_Lambda) is z-flat,
# so the same fitted a0 implies H0_asym = H0 / sqrt(0.7) (Omega_Lambda
# = 0.7, the declared comparison cosmology).
fac_asym = 1.0 / math.sqrt(0.7)
H0_asym = H0_A * fac_asym
band_asym = (band_lo * fac_asym, band_hi * fac_asym)
assert abs(fac_asym - 1.195) < 0.001, fac_asym
assert abs(H0_asym - 78.1) < 0.1, H0_asym
assert abs(band_asym[0] - 78.1) < 0.1 and abs(band_asym[1] - 84.1) < 0.1
P(f"GF-6 PASS: conversion-branch arithmetic (factor {fac_asym:.3f}; "
  f"asymptotic-branch reading {H0_asym:.1f}, band "
  f"{band_asym[0]:.1f}-{band_asym[1]:.1f} -- outside the validity "
  f"domain, quoted only as the section-8 branch disclosure)")

# ================= FIGURE 1: census + two-leg design ================
fig, ax = plt.subplots(figsize=(8.0, 3.9))
ax.set_xlim(0, 100); ax.set_ylim(0, 10); ax.axis('off')

def box(x, y, w, h, txt, fc, ec='black', fs=8.5, tc='black', lw=0.8):
    ax.add_patch(plt.Rectangle((x, y), w, h, facecolor=fc, edgecolor=ec,
                               linewidth=lw, zorder=2))
    ax.text(x + w / 2, y + h / 2, txt, ha='center', va='center', fontsize=fs,
            color=tc, zorder=3)

def arrow(x0, y0, x1, y1):
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle='-|>', color=C_GRAY, lw=1.0),
                zorder=1)

# top row: the catalog split (widths proportional to counts)
tot = n_flow0 + n_anch0 + n_uma
x0, W, ytop, hbox = 2, 96, 7.4, 2.0
wf = W * n_flow0 / tot; wa = W * n_anch0 / tot; wu = W * n_uma / tot
box(x0, ytop, wf, hbox, f"flow distances\n{n_flow0} galaxies\n"
    r"(assume $H_0 = 73$)", '#f3d9c4')
box(x0 + wf, ytop, wa, hbox, f"anchored  {n_anch0}\nTRGB {n_trgb}\n"
    f"Cep {n_cep}, SN {n_sn}", '#cfe0f0')
box(x0 + wf + wa, ytop, wu, hbox, f"Ursa Major\n{n_uma}\n(cluster)", '#e2ecf5')
ax.text(x0, ytop + 2.3, f"SPARC distance-method census ({tot} galaxies)",
        fontsize=10, va='bottom')
# leg boxes
yleg, hleg = 0.6, 3.0
box(2, yleg, 44, hleg, f"flow arm ({n_flowB} usable galaxies)\n"
    r"self-consistent solve for $a_0(H_0)$" + "\nconsistency test only:\n"
    "power-limited at SPARC depth (mass-to-light pole)", '#faf0e6',
    ec=C_FLOW, fs=8)
box(52, yleg, 46, hleg, f"anchored arm ({nA} galaxies)\n"
    f"{nA_anch} anchored + {nA_uma} UMa + {nA_re} CF4-reclassified\n"
    r"direct inversion  $H_0 = 2\pi\,a_0/c$" + "\nthe measurement arm",
    '#eef4fa', ec=C_ANCH, lw=1.4, fs=8)
ax.text(75, 4.55, f"curve-quality cuts: {n_anch0} → {nA_anch}, "
        f"{n_uma} → {nA_uma}", fontsize=7.5, ha='center',
        color=C_GRAY)
# arrows (all from box bottom edges to leg-box top edges)
arrow(x0 + wf * 0.30, ytop - 0.05, 20, yleg + hleg + 0.05)
arrow(x0 + wf * 0.80, ytop - 0.05, 62, yleg + hleg + 0.05)
arrow(x0 + wf + wa / 2, ytop - 0.05, 74, yleg + hleg + 0.05)
arrow(x0 + wf + wa + wu / 2, ytop - 0.05, 86, yleg + hleg + 0.05)
# callouts, placed clear of the arrows
ax.text(11, 5.15, "circular for a direct\ninversion: the distance\n"
        r"scale already assumes $H_0$",
        ha='center', va='center', fontsize=8, color='#7a3c00', style='italic')
ax.text(48.5, 5.15, f"CF4 per-method\nreclassification\n({nA_re} galaxies)",
        fontsize=8, ha='center', va='center', color=C_GRAY)
save(fig, 'figf1_census')
P("GF-1 PASS: figure-1 counts are the parsed census/leg values "
  f"({tot} = {n_flow0}+{n_anch0}+{n_uma}; {nA} = {nA_anch}+{nA_uma}+{nA_re})")

# ================= FIGURE 2: leg-A first light ======================
fig, ax = plt.subplots(figsize=(7.2, 3.6))
ax.axvspan(VAL_LO, VAL_HI, color='#f4f4f4', zorder=0)
ax.axvspan(band_lo, band_hi, color=C_BAND, alpha=0.45, zorder=1)
ax.axvline(PLANCK, color='#333333', ls=':', lw=1.0, zorder=2)
ax.axvline(SH0ES, color='#333333', ls='--', lw=1.0, zorder=2)

rows_f2 = [   # top-down reading order
    ("primary (hierarchical, deep)", H0_A, (pc16, pc84), C_ANCH, 'o', 7.0),
    ("membership variants", None, (H0_var4, H0_A), C_ANCH, 's', 6.0),
    ("form family:  BE", fam['BE'][1], None, C_ANCH, 'o', 4.9),
    ("geometric-mean", fam['gm'][1], None, C_GRAY, 'D', 4.2),
    ("bootstrap form", fam['boot'][1], None, C_GRAY, 'D', 3.5),
    ("p = 0.65 tail", fam['p065'][1], None, C_GRAY, 'D', 2.8),
    ("flat treatment (envelope only)", H0_flat, (H0_flat - sig_flat,
                                                 min(H0_flat + sig_flat,
                                                     80.0)),
     '#aaaaaa', 'v', 1.5),
]
for lab, x, span, c, mk, y in rows_f2:
    if span is not None and x is not None:
        ax.errorbar([x], [y], xerr=[[x - span[0]], [span[1] - x]],
                    fmt=mk, color=c, ms=6, capsize=3, lw=1.4, zorder=4)
    elif span is not None:
        ax.plot(span, [y, y], '-', color=c, lw=2.5, zorder=4,
                solid_capstyle='butt')
        ax.plot(span, [y, y], '|', color=c, ms=8, zorder=5)
    else:
        ax.plot([x], [y], mk, color=c, ms=6, zorder=4)
ax.annotate('', xy=(80.4, 1.5), xytext=(79.4, 1.5),
            arrowprops=dict(arrowstyle='-|>', color='#aaaaaa', lw=1.2))
ax.text(79.9, 1.75, f"to {H0_flat + sig_flat:.0f}", fontsize=7.5,
        color='#888888', ha='center')
yt = [r[5] for r in rows_f2]
ax.set_yticks(yt)
ax.set_yticklabels([r[0] for r in rows_f2], fontsize=8.5)
ax.tick_params(axis='y', length=0)

ax.text(PLANCK - 0.3, 8.15, f"Planck {PLANCK}", ha='right', fontsize=8,
        color='#333333')
ax.text(SH0ES + 0.3, 8.15, f"ladder {SH0ES}", ha='left', fontsize=8,
        color='#333333')
ax.text((band_lo + band_hi) / 2, 0.25, "form band 65.4–70.4",
        ha='center', va='bottom', fontsize=8, color='#1c4c78')
ax.text(VAL_HI - 0.3, 0.82, "validity 62–75", ha='right', va='bottom',
        fontsize=8, color='#888888')
ax.set_xlim(56, 80.5); ax.set_ylim(0, 8.8)
ax.set_xlabel(r"$H_0$  [km s$^{-1}$ Mpc$^{-1}$]")
for s in ('left', 'right', 'top'):
    ax.spines[s].set_visible(False)
save(fig, 'figf2_firstlight')

# ================= FIGURE 3: the per-galaxy map =====================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.6, 3.7),
                               gridspec_kw={'width_ratios': [1.25, 1]})
# (a) map: d_boot vs coverage
for mm, c, lab in ((mA, C_ANCH, f"anchored ({int(mA.sum())})"),
                   (mF, C_FLOW, f"flow ({int(mF.sum())})")):
    ax1.scatter(cov[mm & ~pend], d_boot[mm & ~pend], s=16, color=c,
                alpha=0.75, label=lab, zorder=3, linewidths=0)
ax1.scatter(cov[pend], d_boot[pend], s=26, facecolors='none',
            edgecolors=C_ANCH, linewidths=0.9, zorder=4,
            label='census-pending')
ax1.axhline(0, color='#999999', lw=0.7, zorder=1)
for mm in (mA, mF):
    idx = np.where(mm)[0]
    for i in idx[np.argsort(np.abs(d_boot[idx]))[::-1][:3]]:
        ax1.annotate(str(name[i]), (cov[i], d_boot[i]),
                     textcoords='offset points',
                     xytext=(5, -11 if d_boot[i] < 0 else 5), fontsize=7.5,
                     color=C_ANCH if leg[i] == 'anch' else C_FLOW)
ax1.set_ylim(-56, 26)
ax1.set_xlabel(r"acceleration coverage  (median $\log_{10} g_{\rm bar}/a_0$)")
ax1.set_ylabel(r"$\Delta(-2\ln L)$ per galaxy   (bootstrap $-$ BE form)")
ax1.legend(loc='lower left', frameon=False)
ax1.text(0.02, 0.97, "(a)", transform=ax1.transAxes, va='top', fontsize=10)
# (b) cumulative contribution, galaxies ranked by |contribution|
for mm, c, lab, tot_ in ((mA, C_ANCH, 'anchored', anch_tot),
                         (mF, C_FLOW, 'flow', flow_tot)):
    d = d_boot[mm]
    d = d[np.argsort(np.abs(d))[::-1]]     # largest |contribution| first
    ax2.step(np.arange(1, len(d) + 1), np.cumsum(d), where='mid', color=c,
             lw=1.6)
    ax2.axhline(tot_, color=c, lw=0.7, ls=':', zorder=1)
ax2.axhline(0, color='#999999', lw=0.7)
ax2.text(78, anch_tot + 3, f"anchored: {anch_tot:+.1f}", ha='right',
         va='bottom', fontsize=8.5, color=C_ANCH)
ax2.text(78, flow_tot - 3, f"flow: {flow_tot:+.1f}", ha='right', va='top',
         fontsize=8.5, color=C_FLOW)
rest70 = flow_tot - d_boot[i_ugc]
assert abs(rest70 - (-0.9)) < 0.1, rest70
ax2.annotate(f"UGC03580 alone: {d_boot[i_ugc]:+.1f}\n"
             f"({100*sd2[i_ugc]:.0f}% of the flow block variance);\n"
             f"the remaining {int(mF.sum())-1} net to {rest70:+.1f}",
             xy=(1.4, d_boot[i_ugc]), xytext=(26, -22), fontsize=7.5,
             color=C_FLOW,
             arrowprops=dict(arrowstyle='-', color=C_FLOW, lw=0.7))
ax2.set_xlabel(r"galaxies, ranked by $|$contribution$|$")
ax2.set_ylabel(r"cumulative $\Delta(-2\ln L)$")
ax2.set_ylim(-58, 74)
ax2.text(0.02, 0.97, "(b)", transform=ax2.transAxes, va='top', fontsize=10)
fig.tight_layout()
save(fig, 'figf3_map')

# ================= FIGURE 4: the three-world separation =============
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.0, 3.4),
                               gridspec_kw={'width_ratios': [1.6, 1]})
sig_ax = np.linspace(1.0, 6.0, 200)
worlds = [("L  (true 67; ladder error beyond the pegs)", W_L, W_L, '#2a7a2a', '-'),
          ("N  (new physics; true 73)", W_N, W_N, '#7a2a7a', '-'),
          ("P  (shared pegs ~8% short)", W_P_LO, W_P_HI, '#a02020', '-')]
for lab, wlo, whi, c, ls in worlds:
    lo = np.abs(H0_A - wlo) / sig_ax
    hi = np.abs(H0_A - whi) / sig_ax
    if wlo == whi:
        ax1.plot(sig_ax, lo, ls, color=c, lw=1.6, label=lab)
    else:
        ax1.fill_between(sig_ax, lo, hi, color=c, alpha=0.25, lw=0)
        ax1.plot(sig_ax, lo, '-', color=c, lw=1.0, label=lab)
        ax1.plot(sig_ax, hi, '-', color=c, lw=1.0)
ax1.axhline(2.0, color='#888888', ls=':', lw=1.0)
ax1.text(1.05, 2.06, r"$2\sigma$", ha='right', va='bottom', fontsize=8,
         color='#888888')
ax1.axvline(sigA, color=C_ANCH, ls='--', lw=1.2)
ax1.text(sigA - 0.06, 3.55, f"current stat {sigA:.1f}", fontsize=8,
         color=C_ANCH, va='bottom', ha='left')
ax1.axvspan(2.0, 3.0, color='#eeeeee', zorder=0)
ax1.text(2.5, 0.22, "target\n2–3", ha='center', va='bottom', fontsize=8,
         color='#888888')
ax1.set_xlim(6.0, 1.0)         # sharpening to the right
ax1.set_ylim(0, 6)
ax1.set_xlabel(r"meter precision $\sigma(H_0)$  [km s$^{-1}$ Mpc$^{-1}$]")
ax1.set_ylabel(r"separation from the reading  [$\sigma$]")
ax1.legend(loc='upper left', frameon=False, fontsize=8)
ax1.text(0.02, 0.03, "(a)", transform=ax1.transAxes, va='bottom', fontsize=10)
# (b) the gearing: response of each instrument to a coherent peg rescale
s = np.linspace(-0.08, 0.08, 100)      # fractional peg-distance rescale
ax2.plot(100 * s, 100 * (-1.0) * s, color='#333333', lw=1.4,
         label='ladder  ($\\propto s^{-1}$)')
ax2.plot(100 * s, 100 * (-2.0) * s, color='#999999', lw=1.0, ls=':',
         label='deep-limit asymptote (2)')
ax2.plot(100 * s, 100 * (-GAMMA) * s, color=C_ANCH, lw=1.8,
         label=f'meter  (measured {GAMMA:.2f})')
ax2.axhline(0, color='#bbbbbb', lw=0.6); ax2.axvline(0, color='#bbbbbb', lw=0.6)
ax2.set_xlabel("coherent peg-distance rescale  [%]")
ax2.set_ylabel(r"response of $H_0$ reading  [%]")
ax2.legend(loc='upper right', frameon=False, fontsize=7.5)
ax2.text(0.03, 0.03, "(b)", transform=ax2.transAxes, va='bottom', fontsize=10)
fig.tight_layout()
save(fig, 'figf4_worlds')

# ================= FIGURE 5: the redshift axis ======================
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9.6, 3.2),
                                    gridspec_kw={'width_ratios': [1.3, 1, 1]})
# (a) the lock prediction band + the meter's own z = 0 point
zz = np.linspace(0, 1.2, 300)
lo_band = np.array([a0_lock(z, band_lo) for z in zz]) / UNIT
hi_band = np.array([a0_lock(z, band_hi) for z in zz]) / UNIT
ax1.fill_between(zz, lo_band, hi_band, color=C_BAND, alpha=0.5, lw=0,
                 label=(r"lock $a_0(z)=cH(z)/2\pi$" + "\n"
                        "($H_0$ = form band 65.4–70.4)"))
ax1.plot(zz, lo_band, color=C_ANCH, lw=1.0)
ax1.plot(zz, hi_band, color=C_ANCH, lw=1.0)
a0A_u = a0_A / UNIT
ax1.axhline(a0A_u, color='#888888', ls='--', lw=1.0,
            label=("frozen $a_0$ (no evolution;\nalso the "
                   "asymptotic-rate branch)"))
sig_a0 = a0A_u * sigA / H0_A     # linear conversion of the bootstrap sigma
ax1.errorbar([0.0], [a0A_u], yerr=[[sig_a0], [sig_a0]], fmt='o', color=C_ANCH,
             ms=5, capsize=3, zorder=5, label="this work (anchored leg)")
ax1.set_xlabel("redshift  $z$")
ax1.set_ylabel(r"$a_0$  [$10^{-10}$ m s$^{-2}$]")
ax1.set_xlim(-0.04, 1.2); ax1.set_ylim(0.7, 2.3)
ax1.legend(loc='upper left', frameon=False, fontsize=7.5)
ax1.text(0.03, 0.05, "(a)", transform=ax1.transAxes, fontsize=10)
# (b) MIGHTEE window: evolution-rate comparison (slope space, verbatim rows)
labels = ["HI sample\nalone", "+ SPARC\nanchor"]
vals = [MIG_PURE[2], MIG_ANCH[2]]
errs = [MIG_PURE[3], MIG_ANCH[3]]
ypos = [1, 2]
ax2.axvline(0, color='#888888', ls='--', lw=1.0)
ax2.axvline(a1_lock, color=C_ANCH, lw=1.4)
ax2.text(a1_lock + 0.25, 2.62, "lock\n+0.52", fontsize=8, color=C_ANCH)
ax2.text(0 - 0.25, 2.62, "frozen\n0", fontsize=8, color='#888888',
         ha='right')
ax2.errorbar(vals, ypos, xerr=errs, fmt='o', color='#333333', ms=5,
             capsize=3, lw=1.4)
for y, lab in zip(ypos, labels):
    ax2.text(-7.6, y, lab, ha='left', va='center', fontsize=8)
ax2.set_xlim(-8, 8); ax2.set_ylim(0.4, 3.1)
ax2.set_yticks([])
ax2.set_xlabel(r"$da_0/dz$  [$10^{-10}$ m s$^{-2}$ per unit $z$]")
ax2.set_title(r"$z < 0.09$ (MIGHTEE-HI)", fontsize=9)
ax2.text(0.03, 0.05, "(b)", transform=ax2.transAxes, fontsize=10)
# (c) z ~ 1: back-extrapolated a0(0) intercepts vs the lock intercept
labs3 = ["dark-matter\nframework", "MOND-form\nframework", "per-galaxy"]
tags3 = ['DM', 'MOND', 'pergal']
ax3.axvline(LOCK_Z1, color=C_ANCH, lw=1.4)
ax3.text(LOCK_Z1 + 0.004, 3.35, "lock\n1.047", fontsize=8, color=C_ANCH)
for y, (tag, lab) in enumerate(zip(tags3, labs3), start=1):
    v, e, _ = ciocan[tag]
    ax3.errorbar([v], [y], xerr=[[2 * e], [2 * e]], fmt='none', ecolor='#999999',
                 elinewidth=0.9, capsize=0)
    ax3.errorbar([v], [y], xerr=[[e], [e]], fmt='o', color='#333333', ms=5,
                 capsize=3, lw=1.6)
    ax3.text(0.883, y, lab, ha='left', va='center', fontsize=8)
ax3.set_xlim(0.88, 1.13); ax3.set_ylim(0.4, 3.8)
ax3.set_yticks([])
ax3.set_xlabel(r"back-extrapolated $a_0(z\!=\!0)$  [$10^{-10}$ m s$^{-2}$]")
ax3.set_title(r"$z \sim 1$ (MUSE-DARK III)", fontsize=9)
ax3.text(0.03, 0.05, "(c)", transform=ax3.transAxes, fontsize=10)
fig.tight_layout()
save(fig, 'figf5_zaxis')

# ================= provenance dump ==================================
P("ALL GATES PASS (GF-1..GF-6); five figures written")
with open('data/paperf_figs.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(L_) + "\n")
P("provenance dump: data/paperf_figs.txt")
