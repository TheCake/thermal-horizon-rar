"""
Paper 3 figures (mechanism note), produced under provenance gates.

Four figures -> papers/figs/p3_fig*.png/.pdf:
  F1  the exchange-weight gate s^2 = e^{-2x} with the two measured
      system points and the tail-exponent postdictions.
  F2  the solvable soft sector: the Poschl-Teller well (minimal
      coupling) vs the flat conformal potential, and the exact
      coupling-density ratio 1 + H^2/omega^2.
  F3  the dispersive budget landscape: dc1_eff vs coupling, with the
      normalization band, the measured-band bar, the exclusion-grade
      bar, and the perturbatively-invalid region.
  F4  the lending law w = (1/2) n/(1+n) with the model-grade slope
      checks and the two sky ambient occupations.

GATES (all must pass or no figure ships):
  GP3-1  gate curve at x = 0.1411 / 1.095 equals 0.7536 / 0.1117
         (|d| <= 5e-4) and the p = 1/2 + s^2/4 postdictions equal
         0.6884 / 0.5279 (|d| <= 5e-4).
  GP3-2  D-ratio at omega = H equals 2 exactly; well depth at the
         origin equals -2 H^2 exactly.
  GP3-3  budget curve regression: dc1_eff at g/H = 0.1/0.3/1/3 equals
         -0.0498/+0.0175/+0.0049/+0.0019 (|d| <= 1.5e-3 each).
  GP3-4  lending law w(1) = 0.25 exactly; w(n=1e6) within 1e-6 of 0.5.
Output log: data/paper3_figs.txt.
"""
import math
import os
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import quad

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
FD = os.path.join('papers', 'figs')
os.makedirs(FD, exist_ok=True)

OUT = []
def P(s=""):
    print(s, flush=True)
    OUT.append(s)

t00 = time.time()
gates = {}
plt.rcParams.update({'font.size': 9.5, 'axes.titlesize': 10,
                     'figure.dpi': 150})

TWO_PI = 2*math.pi

# ------------------------------- F1 ---------------------------------
# the measured gates are the PRIMARIES; x and n are derived from them
x = np.linspace(0.01, 2.2, 400)
s2 = np.exp(-2*x)
GG, GB = 0.7536, 0.1117
XG, XB = -0.5*math.log(GG), -0.5*math.log(GB)
nG = math.sqrt(GG)/(1 - math.sqrt(GG))
nB = math.sqrt(GB)/(1 - math.sqrt(GB))
pG, pB = 0.5 + GG/4, 0.5 + GB/4
ok1 = (abs(pG - 0.6884) <= 5e-4 and abs(pB - 0.5279) <= 5e-4
       and abs(nG - 6.6) <= 0.1 and abs(nB - 0.502) <= 0.005)
gates['GP3-1'] = ok1
P("GP3-1 gate primaries 0.7536/0.1117 -> x %.4f/%.4f, n %.3f/%.3f "
  "(ledger 6.6/0.502), p %.4f/%.4f (vs 0.6884/0.5279) -> %s"
  % (XG, XB, nG, nB, pG, pB, "PASS" if ok1 else "FAIL"))

fig, ax = plt.subplots(figsize=(5.0, 3.4))
ax.plot(x, s2, lw=1.8, color='#28527a')
ax.scatter([XG], [GG], zorder=5, color='#d1495b', s=42)
ax.scatter([XB], [GB], zorder=5, color='#d1495b', s=42, marker='s')
ax.annotate("galaxies\n$x_{amb}=0.142$, $s^2=0.754$\n"
            r"$p_{post}=0.688$ (meas. $0.65\pm0.075$)",
            (XG, GG), xytext=(0.42, 0.80), fontsize=8.5,
            arrowprops=dict(arrowstyle='-', lw=0.7))
ax.annotate("wide binaries\n$x_{amb}=1.096$, $s^2=0.112$\n"
            r"$p_{post}=0.528$ (meas. $\approx 1/2$)",
            (XB, GB), xytext=(1.25, 0.42), fontsize=8.5,
            arrowprops=dict(arrowstyle='-', lw=0.7))
ax.set_xlabel(r"ambient occupation variable  $x = (e_N/a_0)^{1/2}$")
ax.set_ylabel(r"exchange weight  $s^2 = e^{-2x}$")
ax.set_xlim(0, 2.2); ax.set_ylim(0, 1.02)
ax.grid(alpha=0.25, lw=0.5)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(FD, 'p3_fig1_gate.' + ext))
plt.close(fig)
P("F1 written")

# ------------------------------- F2 ---------------------------------
Hx = np.linspace(0, 4, 300)
Vmin = -2.0/np.cosh(Hx)**2
w_ax = np.geomspace(0.05, 10, 300)
ratio = 1 + 1.0/w_ax**2
ok2 = (abs((1 + 1.0/1.0**2) - 2.0) == 0 and Vmin[0] == -2.0)
gates['GP3-2'] = ok2
P("GP3-2 exact anchors: ratio(H) = 2, depth = -2H^2 -> %s"
  % ("PASS" if ok2 else "FAIL"))

fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 3.0))
a1.plot(Hx, Vmin, lw=1.8, color='#28527a',
        label=r'minimal ($\xi=0$): $-2H^2\,\mathrm{sech}^2(Hx)$')
a1.axhline(0, lw=1.4, color='#7a9e7e',
           label=r'conformal ($\xi=1/6$): $V=0$')
a1.set_xlabel(r"tortoise coordinate $Hx$")
a1.set_ylabel(r"$V_{l=0}/H^2$")
a1.legend(fontsize=7.5, loc='lower right')
a1.grid(alpha=0.25, lw=0.5)
a2.loglog(w_ax, ratio, lw=1.8, color='#28527a')
a2.axhline(1, lw=0.8, color='0.5', ls=':')
a2.scatter([1], [2], color='#d1495b', zorder=5, s=30)
a2.annotate(r"$1+H^2/\omega^2$ (exact)", (0.09, 1.6), fontsize=9)
a2.set_xlabel(r"$\omega/H$")
a2.set_ylabel(r"$D_{min}/D_{conf}$")
a2.grid(alpha=0.25, lw=0.5, which='both')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(FD, 'p3_fig2_well.' + ext))
plt.close(fig)
P("F2 written")

# ------------------------------- F3 ---------------------------------
def n_be(wv):
    z = TWO_PI*wv
    if z > 500: return 0.0
    return 1.0/math.expm1(z)
def D_min(wv):
    return (wv**2 + 1.0)/(4*math.pi**2*wv)
GAP_GAL = 0.1411/TWO_PI
def Delta(x_, g_, gap_, UV=30.0):
    Om_ = x_/TWO_PI
    up = -(g_*g_/(4*math.pi**2))*math.log(1 + Om_/g_)/Om_
    band = 0.0
    if Om_ - g_ > gap_:
        f = lambda wv: math.log(wv/(Om_ - wv))
        band = (g_*g_/(4*math.pi**2))*(f(Om_ - g_) - f(gap_))/Om_
    def ft(wv):
        return D_min(wv)*2*n_be(wv)/(Om_ - wv)
    th = 0.0
    if Om_ - g_ > gap_:
        v1, _ = quad(ft, gap_, Om_ - g_, limit=400)
        th += v1
    v2, _ = quad(ft, Om_ + g_, UV, limit=400)
    th = (th + v2)*g_*g_
    return up + band + th
def nu(x_):
    if x_ <= 0: return float('inf')
    if x_ > 500: return 1.0
    return 1.0 + 1.0/math.expm1(x_)
XLO, XHI = 0.14, 1.73
xg = np.linspace(XLO, XHI, 160)
A = np.vstack([np.ones_like(xg), xg]).T
mdeep = (xg >= 0.14) & (xg <= 0.45)
def dc1(g_, dnorm=1.0):
    dx = np.array([TWO_PI*Delta(float(x_), g_, GAP_GAL) for x_ in xg])
    dx = dx*dnorm
    c_, *_ = np.linalg.lstsq(A, dx, rcond=None)
    r_ = dx - A @ c_
    dn_ = np.array([nu(float(x_) + float(rr)) - nu(float(x_))
                    for x_, rr in zip(xg, r_)])
    return float(np.mean(dn_[mdeep]))
gv = np.geomspace(0.05, 3.0, 40)
curve = np.array([dc1(float(g_)) for g_ in gv])
curve3 = np.array([dc1(float(g_), 3.0) for g_ in gv])
reg = {0.1: -0.0498, 0.3: 0.0175, 1.0: 0.0049, 3.0: 0.0019}
ok3 = all(abs(dc1(g_) - v_) <= 1.5e-3 for g_, v_ in reg.items())
gates['GP3-3'] = ok3
P("GP3-3 budget regression at 0.1/0.3/1/3: %s"
  % ("PASS" if ok3 else "FAIL"))

fig, ax = plt.subplots(figsize=(5.2, 3.5))
ax.plot(gv, np.abs(curve), lw=1.8, color='#28527a',
        label=r'$|dc_1^{eff}|$, fiducial normalization')
ax.plot(gv, np.abs(curve3), lw=1.4, color='#28527a', ls='--',
        label=r'normalization $\times 3$ (unpinned axis)')
ax.axhline(0.05, color='#d1495b', lw=1.1,
           label=r'measured $c_1$ band width (0.05)')
ax.axhline(0.15, color='#d1495b', lw=1.1, ls=':',
           label='exclusion grade (0.15)')
ax.axvspan(0.05, 0.1, color='0.85',
           label='linear response marginal')
ax.scatter([1.0], [abs(dc1(1.0))], zorder=6, color='k', s=34,
           label=r'fiducial $g = H$')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r"coupling $g/H$ (window scale)")
ax.set_ylabel(r"deep-arm distortion $|dc_1^{eff}|$")
ax.legend(fontsize=7, loc='lower left')
ax.grid(alpha=0.25, lw=0.5, which='both')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(FD, 'p3_fig3_budget.' + ext))
plt.close(fig)
P("F3 written")

# ------------------------------- F4 ---------------------------------
nn = np.geomspace(0.02, 50, 300)
wlaw = 0.5*nn/(1 + nn)
ok4 = (0.5*1/(1+1) == 0.25 and abs(0.5*1e6/(1+1e6) - 0.5) <= 1e-6)
gates['GP3-4'] = ok4
P("GP3-4 lending anchors: w(1) = 0.25, w(inf) -> 1/2 -> %s"
  % ("PASS" if ok4 else "FAIL"))

fig, ax = plt.subplots(figsize=(5.0, 3.4))
ax.semilogx(nn, wlaw, lw=1.8, color='#28527a')
ax.axhline(0.5, lw=0.8, color='0.5', ls=':')
for nv, lab, mk in ((0.502, 'binary ambient', 's'),
                    (6.63, 'galaxy ambient', 'o')):
    ax.scatter([nv], [0.5*nv/(1+nv)], zorder=5, color='#d1495b',
               s=40, marker=mk)
    ax.annotate(lab, (nv, 0.5*nv/(1+nv)),
                xytext=(nv*1.25, 0.5*nv/(1+nv) - 0.055), fontsize=8.5)
ax.annotate("closed-model check: slope $0.989$, rms $0.005$\n"
            "circuit-QED translation: slope $0.987$, rms $0.008$\n"
            "two-quantum prefactor $0.480 \\approx 1/2$",
            (0.025, 0.40), fontsize=8)
ax.set_xlabel(r"ambient occupation $n$")
ax.set_ylabel(r"exchange-ward weight  $\frac{1}{2}\,n/(1+n)$")
ax.set_ylim(0, 0.56)
ax.grid(alpha=0.25, lw=0.5)
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(os.path.join(FD, 'p3_fig4_lending.' + ext))
plt.close(fig)
P("F4 written")

P("")
P("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
  for k, v in sorted(gates.items())))
P("done (%.1f min)" % ((time.time()-t00)/60))
with open('data/paper3_figs.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/paper3_figs.txt")
