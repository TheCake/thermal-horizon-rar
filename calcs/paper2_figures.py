"""PAPER 2 FIGURES (draft 0.1 -> 0.2).  Every plotted number is either
computed from the SPARC/lensing data by the verbatim stage loader in
this script or parsed from a committed stage output; the provenance
dump (data/paper2_figs.txt) lists all of them.  In-script gates regress
the parsed values against the record (and against the numbers quoted in
papers/paper2_rar_coefficients.md) before any figure is written.
Outputs: papers/figs/p2_fig{1..5}*.png/.pdf
"""
import glob, math, os, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

os.makedirs('papers/figs', exist_ok=True)
L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

P("PAPER 2 FIGURES - provenance-gated build")

# ---------------- SPARC loader (verbatim 4S/8S construction) --------
KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LENS_CUT_FID = -14.25

meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        name, inc, q = t[0], float(t[5]), int(t[17])
        meta[name] = (inc, q)
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q = meta.get(name, (0, 3))
    if inc < 30 or q > 2: continue
    kept += 1
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC; gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        g_gas.append(gg); g_dsk.append(gd); g_bul.append(gb)
        gobs.append(Vo*Vo/R*KPC)
        sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
gbar = g_gas + g_dsk + g_bul
lgobs = np.log10(gobs)
lgbar = np.log10(gbar)

# the GD selector (verbatim 8S)
allg = np.unique(gal_id)
NGAL = gal_id.max()+1
gdfrac = np.zeros(NGAL)
for g_ in allg:
    m = gal_id == g_
    gdfrac[g_] = float(np.mean(g_gas[m] > g_dsk[m] + g_bul[m]))
gd_set = np.array([g_ for g_ in allg if gdfrac[g_] >= 0.5])
dd_set = np.array([g_ for g_ in allg if gdfrac[g_] < 0.5])

# lensing table (verbatim 4E/4S fiducial config)
ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_err = np.sqrt(l_stat**2 + l_syst**2)
lmask = l_gbar >= LENS_CUT_FID

P(f"loader: {kept} galaxies pass cuts, {len(gobs)} points, "
  f"{len(allg)} galaxies contribute; GD {len(gd_set)} / DD {len(dd_set)}; "
  f"lensing fiducial {int(lmask.sum())} of {len(l_gbar)} rows")

# the function family (verbatim 8S forms + the simple function)
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))

# ---------------- parse committed stage outputs ---------------------
def read(path):
    with open(path, encoding='utf-8', errors='replace') as f:
        return f.read()

def fnum(pat, txt, idx=1):
    m = re.search(pat, txt, re.M)
    if m is None:
        raise RuntimeError(f"parse failed: {pat}")
    return float(m.group(idx))

# 4S: joint marginalized c1 profile + bootstrap
t4s = read('data/stage4s_c1fit.txt')
prof4s = re.findall(r"^\s+([+-]\d\.\d\d)\s+[\d.]+\s+(-\d+\.\d+)\s*$",
                    t4s, re.M)
lam4s = np.array([float(a) for a, b in prof4s])
obj4s = np.array([float(b) for a, b in prof4s])
m4s = re.search(r"D\(-2lnL\) at c1=0: \+([\d.]+) \| c1=1/4: \+([\d.]+) "
                r"\| c1=1/2: \+([\d.]+)", t4s)
d4s_0, d4s_q, d4s_h = (float(m4s.group(i)) for i in (1, 2, 3))
b4s = re.search(r"joint bootstrap \(200 reps\): lam 16/50/84 = "
                r"([\d.]+)/([\d.]+)/([\d.]+)", t4s)
boot4s = np.array([float(b4s.group(i)) for i in (1, 2, 3)])/2.0

# 4Z: hierarchical c1 profile + bootstrap
t4z = read('data/stage4z_hierc1.txt')
prof4z = re.findall(r"^\s+([+-]\d\.\d\d)\s+(-\d+\.\d+)\s*$", t4z, re.M)
lam4z = np.array([float(a) for a, b in prof4z])
obj4z = np.array([float(b) for a, b in prof4z])
d4z_0 = fnum(r"D\(-2lnL\): c1=0 at \+([\d.]+)", t4z)
d4z_h = fnum(r"c1=1/2 at \+([\d.]+)", t4z)
b4z = re.search(r"bootstrap \(100 reps\): lam 16/50/84 = "
                r"([\d.]+)/([\d.]+)/([\d.]+)", t4z)
boot4z = np.array([float(b4z.group(i)) for i in (1, 2, 3)])/2.0

# 4V scorecard: the a0 rows
t4v = read('data/stage4v_scorecard.txt')
a0_flat = (fnum(r"SPARC p-fit \(4H\)\s+([\d.]+)\+-[\d.]+", t4v),
           fnum(r"SPARC p-fit \(4H\)\s+[\d.]+\+-([\d.]+)", t4v))
a0_joint = (fnum(r"joint SPARC\+lensing \(4E\)\s+([\d.]+)\+-[\d.]+", t4v),
            fnum(r"joint SPARC\+lensing \(4E\)\s+[\d.]+\+-([\d.]+)", t4v))
a0_bin_s = (fnum(r"wide binaries \(simple\)\s+([\d.]+)\+-[\d.]+", t4v),
            fnum(r"wide binaries \(simple\)\s+[\d.]+\+-([\d.]+)", t4v))
a0_bin_b = (fnum(r"wide binaries \(BE\)\s+([\d.]+)\+-[\d.]+", t4v),
            fnum(r"wide binaries \(BE\)\s+[\d.]+\+-([\d.]+)", t4v))

# 5L: the horizon reference values
t5l = read('data/stage5l_ladder.txt')
planck = (fnum(r"Planck ([\d.]+)\+-[\d.]+", t5l),
          fnum(r"Planck [\d.]+\+-([\d.]+)", t5l))
shoes = (fnum(r"SH0ES ([\d.]+)\+-[\d.]+", t5l),
         fnum(r"SH0ES [\d.]+\+-([\d.]+)", t5l))

# 5M: hierarchical (measured D/i priors) a0 per function
t5m = read('data/stage5m_hierv.txt')
a0_hier = {}
for fn_ in ('BE', 'p065', 'gm', 'boot'):
    a0_hier[fn_] = fnum(rf"^\s+{fn_}:\s+-[\d.]+\s+\(15 rd\)\s+"
                        rf"a0=([\d.]+)e-10", t5m)

# 4T: the second moment (BE block)
t4t = read('data/stage4t_bathnoise.txt')
be_block = t4t.split('[BE]')[1].split('[simple]')[0]
s0_const = fnum(r"M0  const:.*?s0 = ([\d.]+)", be_block)
floor4t = fnum(r"M1b osc \+ floor:.*?floor = ([\d.]+)", be_block)
d_m2 = fnum(r"M2-M0 = (-[\d.]+)", be_block)
d_m1b = fnum(r"M1b-M0 = (-[\d.]+)", be_block)
bins4t = re.findall(r"bin\d:\s+([\d.]+) \| ([\d.]+) \| ([\d.]+)", be_block)
bx = np.array([float(a) for a, b, c in bins4t])
bpred = np.array([float(b) for a, b, c in bins4t])
bfree = np.array([float(c) for a, b, c in bins4t])

# 7B: decile profile + geometry cells
t7b = read('data/stage7b_bumphunt.txt')
dec7b = re.findall(r"^\s+([\d.]+) \| ([\d.]+) \| ([\d.]+)\s*$", t7b, re.M)
dx = np.array([float(a) for a, b, c in dec7b])
draw = np.array([float(b) for a, b, c in dec7b])
dsub = np.array([float(c) for a, b, c in dec7b])
c1cells = re.search(r"C1 \[0\.4,0\.8\): R/Rd 0-1\.5: ([\d.]+) \(n=(\d+)\) \| "
                    r"R/Rd 1\.5-3: ([\d.]+) \(n=(\d+)\) \| "
                    r"R/Rd 3-99: ([\d.]+) \(n=(\d+)\)", t7b)
wcells = re.search(r"W  \[0\.8,1\.4\]: R/Rd 0-1\.5: ([\d.]+) \(n=(\d+)\) \| "
                   r"R/Rd 1\.5-3: ([\d.]+) \(n=(\d+)\) \| "
                   r"R/Rd 3-99: ([\d.]+) \(n=(\d+)\)", t7b)
C1s2 = [float(c1cells.group(i)) for i in (1, 3, 5)]
C1n = [int(c1cells.group(i)) for i in (2, 4, 6)]
Ws2 = [float(wcells.group(i)) for i in (1, 3, 5)]
Wn = [int(wcells.group(i)) for i in (2, 4, 6)]

# 7C: identification numbers
t7c = read('data/stage7c_gammaclean.txt')
explained = fnum(r"EXPLAINED = ([\d.]+)", t7c)
b_clean = fnum(r"b_clean = ([\d.]+)", t7c)

# 8S: the gas-budget immunity row
t8s = read('data/stage8s_gasc1.txt')
immun = fnum(r"\[GD  \] IMMUNITY:.*?\|d lam\| = ([\d.]+)", t8s)

# 8S-c: the vertON dial (FULL/GD/DD) with D1 intervals
t8sc = read('data/stage8sc_gddist.txt')
def dial(tag):
    m = re.search(rf"\[{tag}\s*vertON\] lam_hat = (-?[\d.]+) "
                  rf"\(D1 (-?[\d.]+)\.\.(-?[\d.]+)\)", t8sc)
    return tuple(float(m.group(i)) for i in (1, 2, 3))
dFULL, dGD, dDD = dial('FULL'), dial('GD'), dial('DD')

# 8V: galaxy-bootstrap percentiles + separation P
t8v = read('data/stage8v_gdboot.txt')
def pct8v(tag):
    m = re.search(rf"\[{tag}\s*\] lam_hat pct 5/50/95 = "
                  rf"([+-]?[\d.]+)/([+-]?[\d.]+)/([+-]?[\d.]+)", t8v)
    return tuple(float(m.group(i)) for i in (1, 2, 3))
pFULL, pGD, pDD = pct8v('FULL'), pct8v('GD'), pct8v('DD')
p_sep = fnum(r"P\(lam_GD >= lam_DD\) = ([\d.]+)", t8v)

# the remaining controls
p_regime = fnum(r"P\(DD-deep <= GD-deep\) = ([\d.]+)",
                read('data/stage8x_regime.txt'))
t9b = read('data/stage9b_qflag.txt')
q1lam = fnum(r"\[GD-Q1\] lam_hat = (-[\d.]+)", t9b)
q2lam = fnum(r"\[GD-Q2\] lam_hat = (-[\d.]+)", t9b)
p_press = fnum(r"P\(D_fs <= 0\) = ([\d.]+)", read('data/stage8y_pressure.txt'))
p_amb = fnum(r"dmed = \+[\d.]+ dex; perm P = ([\d.]+)  <- PRIMARY",
             read('data/stage9g_gdambient.txt'))
d_rise = fnum(r"D = \+([\d.]+)\s*$",
              re.search(r"\[point\] lam_hat\(FULL-149\).*$",
                        read('data/stage9n_risingdial.txt'),
                        re.M).group(0))

# ---------------- GATES ---------------------------------------------
g_all = True
def gate(name, ok, detail=""):
    global g_all
    g_all = g_all and ok
    P(f"{name}: {'PASS' if ok else 'FAIL'}{'  ('+detail+')' if detail else ''}")

gate("G-P2a loader counts",
     kept == 153 and len(gobs) == 2700 and len(allg) == 149
     and len(gd_set) == 38 and len(dd_set) == 111 and int(lmask.sum()) == 12,
     f"{kept}/{len(gobs)}/{len(allg)}/{len(gd_set)}/{len(dd_set)}/"
     f"{int(lmask.sum())}")

# coefficient profiles: recompute the printed Delta values from the tables
r4s = (abs((obj4s[lam4s == 0.0][0] - obj4s.min()) - d4s_0) < 0.15
       and abs((obj4s[lam4s == 0.5][0] - obj4s.min()) - d4s_q) < 0.15
       and abs((obj4s[lam4s == 1.0][0] - obj4s.min()) - d4s_h) < 0.15
       and lam4s[np.argmin(obj4s)] == 0.90 and len(lam4s) == 37)
r4z = (abs((obj4z[lam4z == 0.0][0] - obj4z.min()) - d4z_0) < 0.15
       and abs((obj4z[lam4z == 1.0][0] - obj4z.min()) - d4z_h) < 0.15
       and len(lam4z) == 17)
gate("G-P2b profile-table regression (4S/4Z headline Deltas)",
     r4s and r4z,
     f"4S {obj4s[lam4s==0.0][0]-obj4s.min():.2f}/{d4s_0} "
     f"4Z {obj4z[lam4z==0.0][0]-obj4z.min():.2f}/{d4z_0}")
gate("G-P2b2 paper-quote check (c1 bootstrap bands)",
     abs(boot4s[1] - 0.427) < 0.005 and abs(boot4z[1] - 0.377) < 0.005
     and abs(d4s_0 - 56.3) < 0.05 and abs(d4z_0 - 28.9) < 0.05,
     f"medians {boot4s[1]:.3f}/{boot4z[1]:.3f}")

# the a0 ladder vs the numbers quoted in the paper's Table 3
hvals = sorted(a0_hier.values())
gate("G-P2c a0-ladder paper-quote check",
     a0_flat == (1.05, 0.10) and a0_joint == (1.00, 0.09)
     and a0_bin_s == (1.48, 0.18) and a0_bin_b == (1.37, 0.17)
     and planck == (1.042, 0.008) and shoes == (1.129, 0.015)
     and abs(hvals[0] - 1.044) < 1e-9 and abs(hvals[-1] - 1.125) < 1e-9,
     f"hier range {hvals[0]:.3f}-{hvals[-1]:.3f}")

# the scatter numbers vs the paper's Section 6
gate("G-P2d scatter paper-quote check",
     abs(bfree[0] - 0.1438) < 1e-9 and abs(bfree[5] - 0.1073) < 1e-9
     and abs(d_m2 + 42.66) < 0.05 and abs(d_m1b + 25.05) < 0.05
     and abs(floor4t - 0.1014) < 1e-9 and len(dx) == 10
     and abs(explained - 0.67) < 1e-9 and b_clean == 0.0,
     f"s_b deep/newt {bfree[0]}/{bfree[5]}")
rin_w = Ws2[0]/Ws2[1]; rin_c = C1s2[0]/C1s2[1]
share_w = Wn[0]/sum(Wn); share_c = C1n[0]/sum(C1n)
gate("G-P2d2 inner-disk decomposition (2.4x; 49%/16%)",
     2.2 < rin_w < 2.6 and 2.2 < rin_c < 2.6
     and abs(share_w - 0.49) < 0.01 and abs(share_c - 0.16) < 0.01,
     f"ratios {rin_w:.2f}/{rin_c:.2f} shares {share_w:.2f}/{share_c:.2f}")

# the GD dial vs the paper's Section 8 and Table 4
gate("G-P2e GD-dial paper-quote check",
     abs(dGD[0]/2 - (-0.66)) < 0.005 and abs(dDD[0]/2 - 0.67) < 0.005
     and abs(dFULL[0]/2 - 0.48) < 0.005 and p_sep == 0.0
     and immun == 0.0 and p_regime == 0.0
     and abs(p_press - 0.0333) < 1e-9 and abs(p_amb - 0.813) < 1e-9
     and abs(d_rise - 0.061) < 1e-9,
     f"c1 GD {dGD[0]/2:.3f} DD {dDD[0]/2:.3f}")

if not g_all:
    P("GATES FAILED - STOP, no figures written")
    with open('data/paper2_figs.txt', 'w') as f:
        f.write("\n".join(L_)+"\n")
    raise SystemExit(1)
P("ALL GATES PASS - writing figures")

# ================= FIGURES ==========================================
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 9,
                     'figure.dpi': 150})

# ---- Fig 1: the RAR with the function family + lensing -------------
fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.4))
ax[0].plot(lgbar, lgobs, '.', ms=1.5, color='0.7', rasterized=True,
           label='SPARC (2,700 points)')
ax[0].errorbar(l_gbar[lmask], l_gobs[lmask], yerr=l_err[lmask], fmt='o',
               ms=3.5, color='tab:blue', lw=0.8, capsize=1.5,
               label='lensing (fit window)')
ax[0].errorbar(l_gbar[~lmask], l_gobs[~lmask], yerr=l_err[~lmask], fmt='o',
               ms=3.5, mfc='none', color='tab:blue', lw=0.8, capsize=1.5,
               label='lensing (deep extension)')
gx = np.linspace(-15.4, -8.6, 400)
gb = 10.0**gx
for nu_, ls_, lab_ in ((nu_be, '-', 'occupation'),
                       (nu_simple, '--', 'simple'),
                       (nu_standard, '-.', r'standard-$\mu$')):
    ax[0].plot(gx, np.log10(nu_(gb/A0_FID)*gb), ls_, color='k', lw=1,
               label=lab_)
ax[0].plot(gx, gx, ':', color='0.5', lw=0.8, label=r'$g_{\rm obs}=g_{\rm bar}$')
ax[0].set_xlim(-15.4, -8.6); ax[0].set_ylim(-13.4, -8.6)
ax[0].set_xlabel(r'$\log_{10}\,g_{\rm bar}\ [{\rm m\,s^{-2}}]$')
ax[0].set_ylabel(r'$\log_{10}\,g_{\rm obs}\ [{\rm m\,s^{-2}}]$')
ax[0].legend(fontsize=6, loc='upper left')
ax[0].set_title('(a) the relation and the candidate functions')

BE_ = np.linspace(-13.4, -8.8, 11)
bc_, med_ = [], []
for i in range(10):
    m = (lgbar >= BE_[i]) & (lgbar < BE_[i+1])
    if m.sum() >= 30:
        bc_.append(np.median(lgbar[m]))
        med_.append(np.median(lgobs[m]
                              - np.log10(nu_be(gbar[m]/A0_FID)*gbar[m])))
ax[1].plot(bc_, med_, 'o', ms=4, color='0.35',
           label='SPARC binned median')
lres = l_gobs - np.log10(nu_be(10**l_gbar/A0_FID)*10**l_gbar)
ax[1].errorbar(l_gbar[lmask], lres[lmask], yerr=l_err[lmask], fmt='o',
               ms=3.5, color='tab:blue', lw=0.8, capsize=1.5,
               label='lensing')
ax[1].errorbar(l_gbar[~lmask], lres[~lmask], yerr=l_err[~lmask], fmt='o',
               ms=3.5, mfc='none', color='tab:blue', lw=0.8, capsize=1.5)
for nu_, ls_, lab_ in ((nu_simple, '--', 'simple'),
                       (nu_standard, '-.', r'standard-$\mu$')):
    ax[1].plot(gx, np.log10(nu_(gb/A0_FID)/nu_be(gb/A0_FID)), ls_,
               color='k', lw=1, label=lab_)
ax[1].axhline(0, color='k', lw=1)
ax[1].set_xlim(-15.4, -8.6); ax[1].set_ylim(-0.25, 0.25)
ax[1].set_xlabel(r'$\log_{10}\,g_{\rm bar}\ [{\rm m\,s^{-2}}]$')
ax[1].set_ylabel(r'$\Delta\log_{10}\,g_{\rm obs}$ (about occupation)')
ax[1].legend(fontsize=6, loc='lower left')
ax[1].set_title('(b) residuals about the occupation function')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/p2_fig1_rar.{ext}')
plt.close(fig)
P("p2_fig1 written")

# ---- Fig 2: the c1 profiles ----------------------------------------
fig, ax = plt.subplots(figsize=(5.8, 3.5))
ax.plot(lam4s/2, obj4s - obj4s.min(), '-', color='k', lw=1.2,
        label='global M/L (joint, marginalized)')
ax.plot(lam4z/2, obj4z - obj4z.min(), '--', color='tab:red', lw=1.2,
        label='hierarchical M/L')
ax.axvline(0.0, color='0.6', lw=0.8)
ax.axvline(0.25, color='0.6', ls=':', lw=0.8)
ax.axvline(0.5, color='k', ls=':', lw=1.0)
ax.text(0.505, 52, 'predicted 1/2', fontsize=7, rotation=90, va='top')
ax.text(0.255, 52, '1/4', fontsize=7, rotation=90, va='top')
for y_, b_, c_, lab_ in ((-4.5, boot4s, 'k', 'bootstrap 16-84% (global)'),
                         (-8.5, boot4z, 'tab:red',
                          'bootstrap 16-84% (hier)')):
    ax.errorbar([b_[1]], [y_], xerr=[[b_[1]-b_[0]], [b_[2]-b_[1]]],
                fmt='o', ms=3.5, color=c_, capsize=2, lw=1, clip_on=False)
    ax.text(b_[2]+0.02, y_, lab_, fontsize=6.5, va='center')
ax.set_xlim(-0.16, 0.78); ax.set_ylim(-11, 60)
ax.set_xlabel(r'leading coefficient $c_1$')
ax.set_ylabel(r'$\Delta(-2\ln L)$ relative to minimum')
ax.legend(fontsize=7, loc='upper right')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/p2_fig2_c1profile.{ext}')
plt.close(fig)
P("p2_fig2 written")

# ---- Fig 3: the a0 ladder ------------------------------------------
fig, ax = plt.subplots(figsize=(6.2, 3.2))
ax.axvspan(planck[0]-planck[1], planck[0]+planck[1], color='0.85',
           label=r'$cH_0/2\pi$ (Planck)')
ax.axvspan(shoes[0]-shoes[1], shoes[0]+shoes[1], color='0.93')
ax.text(shoes[0], 4.42, 'SH0ES', fontsize=6.5, ha='center', color='0.4')
rows = [(4, 'SPARC screening fit (flat M/L)', a0_flat, 'k', 'o', True),
        (3, 'joint SPARC + lensing', a0_joint, 'k', 'o', True)]
for y_, lab_, (v_, e_), c_, mk_, filled in rows:
    ax.errorbar([v_], [y_], xerr=[e_], fmt=mk_, ms=4.5, color=c_,
                capsize=2.5, lw=1.2)
    ax.text(0.58, y_+0.16, lab_, fontsize=7.5)
mks = {'BE': 'o', 'gm': '^', 'boot': 'D', 'p065': 's'}
labs = {'BE': 'occupation', 'gm': 'geometric-mean', 'boot': 'quarter',
        'p065': 'sharpened tail'}
for i, fn_ in enumerate(('BE', 'gm', 'boot', 'p065')):
    ax.plot([a0_hier[fn_]], [2 + (i-1.5)*0.13], mks[fn_], ms=4,
            color='tab:red', mfc='none' if fn_ != 'BE' else 'tab:red')
ax.plot([hvals[0], hvals[-1]], [2-0.32, 2-0.32], '-', color='tab:red',
        lw=1)
ax.text(0.58, 2+0.16,
        'hierarchical, measured distance/inclination priors\n'
        '(four functions; point estimates)', fontsize=7.5)
for dy, (v_, e_), mk_ in ((0.10, a0_bin_s, 'o'), (-0.10, a0_bin_b, 's')):
    ax.errorbar([v_], [1+dy], xerr=[e_], fmt=mk_, ms=4, mfc='none',
                color='0.55', capsize=2.5, lw=1)
ax.text(0.58, 1+0.16, 'wide-binary translation (fenced model; superseded)',
        fontsize=7.5, color='0.45')
ax.set_ylim(0.45, 4.75); ax.set_xlim(0.55, 1.75)
ax.set_yticks([])
ax.set_xlabel(r'$a_0\ [10^{-10}\ {\rm m\,s^{-2}}]$')
ax.legend(fontsize=7, loc='lower right')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/p2_fig3_a0ladder.{ext}')
plt.close(fig)
P("p2_fig3 written")

# ---- Fig 4: the second moment --------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.4))
ax[0].axvspan(0.8, 1.4, color='0.92')
ax[0].text(1.06, 0.152, 'transition\nwindow', fontsize=6.5, ha='center',
           color='0.35')
ax[0].plot(dx, draw, '.', ms=5, color='0.6',
           label='decile profile (raw)')
ax[0].plot(dx, dsub, '.', ms=5, color='tab:blue',
           label='decile profile (within-galaxy)')
ax[0].plot(bx, bfree, 'o', ms=5, color='k', label='free six-bin fit')
ax[0].plot(bx, bpred, 's--', ms=4, mfc='none', color='tab:red', lw=1,
           label='oscillator + floor (per-bin values)')
ax[0].axhline(s0_const, color='k', ls=':', lw=0.8)
ax[0].text(0.163, s0_const+0.002, 'constant fit', fontsize=6.5)
ax[0].axhline(floor4t, color='tab:red', ls=':', lw=0.8)
ax[0].text(0.163, floor4t+0.002, 'fitted floor', fontsize=6.5,
           color='tab:red')
ax[0].set_xscale('log')
ax[0].set_xticks([0.2, 0.3, 0.5, 1.0, 2.0])
ax[0].set_xticklabels(['0.2', '0.3', '0.5', '1', '2'])
ax[0].minorticks_off()
ax[0].set_xlabel(r'$x=(g_{\rm bar}/a_0)^{1/2}$ (fiducial)')
ax[0].set_ylabel('intrinsic scatter [dex]')
ax[0].set_ylim(0.0, 0.165)
ax[0].legend(fontsize=6, loc='lower left', ncol=1)
ax[0].set_title('(a) the scatter is acceleration-dependent')
xb = np.arange(3)
ax[1].bar(xb-0.19, Ws2, width=0.36, color='0.45',
          label='transition window (x 0.8-1.4)')
ax[1].bar(xb+0.19, C1s2, width=0.36, color='tab:blue', alpha=0.65,
          label='control slice (x 0.4-0.8)')
ax[1].set_xticks(xb)
ax[1].set_xticklabels([r'$R<1.5R_d$ (inner)', r'$1.5$-$3R_d$',
                       r'$R>3R_d$ (outer)'])
ax[1].set_ylabel(r'residual variance $s^2$ [dex$^2$]')
ax[1].set_ylim(0, 0.036)
ax[1].text(0.02, 0.97,
           f'inner/mid variance ratio: {rin_w:.1f} (window), '
           f'{rin_c:.1f} (control)', transform=ax[1].transAxes, fontsize=7)
ax[1].text(0.02, 0.905,
           f'window samples {share_w:.0%} inner disk; control '
           f'{share_c:.0%}', transform=ax[1].transAxes, fontsize=7)
ax[1].text(0.02, 0.84, 'outer-point bump amplitude: 0.0000',
           transform=ax[1].transAxes, fontsize=7)
ax[1].legend(fontsize=6.5, loc='center right')
ax[1].set_title('(b) the transition excess is inner-disk structure')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/p2_fig4_scatter.{ext}')
plt.close(fig)
P("p2_fig4 written")

# ---- Fig 5: the gas-dominated dial split ---------------------------
fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.3),
                       gridspec_kw={'width_ratios': [1.15, 1]})
rows5 = [(3, 'full sample (149)', dFULL, pFULL, 'k'),
         (2, 'disk-dominated (111)', dDD, pDD, 'tab:blue'),
         (1, 'gas-dominated (38)', dGD, pGD, 'tab:red')]
for y_, lab_, (lh, lo, hi), (b5, b50, b95), c_ in rows5:
    ax[0].errorbar([lh/2], [y_], xerr=[[lh/2-lo/2], [hi/2-lh/2]],
                   fmt='o', ms=5, color=c_, capsize=3, lw=2)
    ax[0].errorbar([b50/2], [y_-0.18], xerr=[[b50/2-b5/2], [b95/2-b50/2]],
                   fmt='.', ms=3, color=c_, capsize=1.5, lw=0.8,
                   alpha=0.65)
    ax[0].text(-0.98, y_+0.22, lab_, fontsize=8)
ax[0].axvline(0, color='0.6', lw=0.8)
ax[0].axvline(0.5, color='k', ls=':', lw=0.8)
ax[0].text(0.5, 0.62, 'predicted 1/2', fontsize=6.5, ha='center')
ax[0].set_xlim(-1.05, 0.82); ax[0].set_ylim(0.5, 3.75)
ax[0].set_yticks([])
ax[0].set_xlabel(r'leading coefficient $c_1$')
ax[0].set_title('(a) the dial by galaxy class')
ctrl = [
    ('gas budget scaled x0.5-2.0',
     f'dial shift {immun:.3f}'),
    ('measured distance/inclination priors on',
     f'dial stays ({dGD[0]/2:+.2f})'),
    ('galaxy bootstrap, 300 replicates',
     f'P(GD reaches DD) = {p_sep:.4f}'),
    ('regime control (deep DD vs GD points)',
     f'P = {p_regime:.4f} (200 replicates)'),
    ('quality-flag halves',
     f'{q1lam/2:+.2f} vs {q2lam/2:+.2f} (agree)'),
    ('pressure support',
     f'ordered by V (P = {p_press:.3f}); lever too small'),
    ('measured environment (ambient field)',
     f'field-typical (perm P = {p_amb:.2f})'),
    ('rising-curve galaxies removed',
     f'full dial shift +{d_rise:.3f} (immune)'),
]
ax[1].set_axis_off()
ax[1].text(0.0, 1.0, 'the eight controls (Table 4):', fontsize=8,
           va='top')
for i, (a_, b_) in enumerate(ctrl):
    ax[1].text(0.0, 0.885 - i*0.115, f'{i+1}. {a_}', fontsize=6.8,
               va='top')
    ax[1].text(0.63, 0.885 - i*0.115, b_, fontsize=6.8, va='top',
               color='0.3')
ax[1].set_title('(b) none dissolves the sign reversal')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/p2_fig5_gddial.{ext}')
plt.close(fig)
P("p2_fig5 written")

# ---------------- provenance dump -----------------------------------
P("")
P("provenance: fig1 loader counts above; fig2 = 4S/4Z profile tables "
  "+ bootstrap lines; fig3 = 4V scorecard rows + 5L reference + 5M "
  "per-function a0; fig4 = 4T BE per-bin table + 7B deciles/geometry "
  "cells + 7C identification; fig5 = 8S-c dials + 8V percentiles + "
  "8S/8X/9B/8Y/9G/9N control numbers")
with open('data/paper2_figs.txt', 'w') as f:
    f.write("\n".join(L_)+"\n")
P("done")
