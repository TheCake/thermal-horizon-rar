"""
STAGE 5L (O10b): the a0 ladder under the sharpened functions.

Assembles the self-consistent a0 comparison for p065 and the geometric-mean
bootstrap, the two functions the 5B-5K cascade promoted:
  - galaxy FLAT a0 (5J joint SPARC+lensing, the scorecard-comparable
    convention),
  - galaxy HIER a0 (disclosed with its f_ML trade -- not a clean a0),
  - binary a0 via the 4V translation a0 = 1.2e-10 * alpha^(1/kappa),
    kappa = -[dln(B-1)/dlny + dln(B-1)/dln e_N] from the e_N = 1.0/1.2/1.4
    EFE tables (5K), averaged over the actual 6-30 kAU deep-pair y
    distribution (4V machinery verbatim);
  - alpha-hat per function = mean over the 5K four-seed fits (realization
    SE disclosed), parsed from data/stage5k_summary.txt.
Gates: G1 alpha=1 -> a0 = 1.2e-10 identically; G2 kappa > 0 and
family-stable; G3 regression -- the BE kappa recomputed here must match
Stage 4V's stored value to 0.01.
Writes data/stage5l_ladder.txt.
"""
import math, re
import numpy as np
from astropy.io import fits

A0_TAB = 1.2e-10
C = 2.998e8

def load_tab(path):
    t = np.load(path)
    y, b = t[0][::-1], t[1][::-1]
    return np.log(y), np.log(np.maximum(b-1.0, 1e-12))

TAGS = {1.0: 'g1p0', 1.2: 'g1p2', 1.4: 'g1p4'}
tabs = {}
for fam in ('be', 'p065', 'gm'):
    for e, tag in TAGS.items():
        tabs[(fam, e)] = load_tab(f'data/efe_boost_{fam}_{tag}.npy')

def lnB1(fam, e, lny):
    lo, hi = (1.0, 1.2) if e <= 1.2 else (1.2, 1.4)
    w = (e-lo)/(hi-lo)
    la = np.interp(lny, *tabs[(fam, lo)])
    lb = np.interp(lny, *tabs[(fam, hi)])
    return (1-w)*la + w*lb

def kappa(fam, lny, e=1.184, dl=0.10, de=0.08):
    dy = (lnB1(fam, e, lny+dl) - lnB1(fam, e, lny-dl))/(2*dl)
    dE = (lnB1(fam, e+de, lny) - lnB1(fam, e-de, lny))/(2*de*1.0/e)
    return -(dy + dE)

# deep-sample y distribution (4V verbatim)
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
GKMS = 0.8868
UNIT = 6.685e-9
yd = (GKMS*Mtot[deep]/s_kau[deep]**2)*UNIT/A0_TAB
lnyd = np.log(np.clip(yd, 1e-3, 50.0))

# alpha-hats from the 5K summary
txt = open('data/stage5k_summary.txt').read()
ah = {}
for law in ('p050', 'p0578', 'p065', 'gm'):
    vals = [float(m) for m in re.findall(
        rf'seed \d+ {law}: a_hat=([0-9.]+)', txt)]
    ah[law] = (np.mean(vals), np.std(vals, ddof=1)/math.sqrt(len(vals)),
               len(vals))

L = [f"STAGE 5L a0 ladder; deep pairs (6-30 kAU): {int(deep.sum())}, "
     f"y 16/50/84 = {np.percentile(yd,[16,50,84]).round(3).tolist()}", ""]
L.append("alpha-hat (5K, mean +- realization SE over seeds):")
for law in ('p050', 'p0578', 'p065', 'gm'):
    m, s, n = ah[law]
    L.append(f"  {law:>6}: {m:.3f} +- {s:.3f}  ({n} seeds)")
L.append("")

kap = {}
for fam in ('be', 'p065', 'gm'):
    kv = kappa(fam, lnyd)
    kap[fam] = float(np.mean(kv))
    L.append(f"kappa({fam:>5}) = {np.mean(kv):+.3f} "
             f"(16/84 {np.percentile(kv,16):+.3f}/{np.percentile(kv,84):+.3f})")
ok3 = abs(kap['be'] - 0.916) < 0.01
L.append(f"G3 BE-kappa regression vs 4V (+0.916): {kap['be']:+.3f} -> "
         f"{'PASS' if ok3 else 'FAIL'}")
g2 = all(k > 0 for k in kap.values())
L.append(f"G2 kappa > 0 all families -> {'PASS' if g2 else 'FAIL'}")

def a0_of(alpha, salpha, fam):
    k = kap[fam]
    a0 = A0_TAB*alpha**(1.0/k)
    return a0, a0*salpha/(alpha*k)

g1v, _ = a0_of(1.0, 0.0, 'p065')
L.append(f"G1 alpha=1 -> {g1v:.3e} (1.200e-10) -> "
         f"{'PASS' if abs(g1v-A0_TAB) < 1e-14 else 'FAIL'}")
L.append("")

def cho2pi(H0, sH0):
    h = H0*1e3/3.0857e22
    v = C*h/(2*math.pi)
    return v, v*sH0/H0
tp, stp = cho2pi(67.4, 0.5)
ts, sts = cho2pi(73.0, 1.0)

L.append("THE p-FAMILY a0 LADDER  [1e-10 m/s^2]")
L.append(f"  cH0/2pi: Planck {tp*1e10:.3f}+-{stp*1e10:.3f}, "
         f"SH0ES {ts*1e10:.3f}+-{sts*1e10:.3f}")
GAL_FLAT = {'p065': (1.084, 0.09), 'gm': (1.023, 0.09)}   # 5J; sigma ~ 4E's
GAL_HIER = {'p065': (0.888, None, 1.53), 'gm': (0.836, None, 1.55)}
for fam in ('p065', 'gm'):
    a, sa = GAL_FLAT[fam]
    pull = (a-tp*1e10)/math.hypot(sa, stp*1e10)
    L.append(f"  {fam:>5} galaxy FLAT (5J):  {a:.2f}+-{sa:.2f}   "
             f"pull vs Planck {pull:+.1f}s")
    ah_key = fam
    m, s, n = ah[ah_key]
    # amplitude error: realization SE (+) the 3V-grade systematic 0.11 in
    # quadrature -- disclosed, dominated by the latter
    stot = math.hypot(s, 0.11)
    a0b, sa0b = a0_of(m, stot, fam)
    pullb = (a0b*1e10-tp*1e10)/math.hypot(sa0b*1e10, stp*1e10)
    L.append(f"  {fam:>5} binary (kappa):    {a0b*1e10:.2f}+-{sa0b*1e10:.2f}"
             f"   pull vs Planck {pullb:+.1f}s   "
             f"(alpha {m:.2f}+-{stot:.2f}, kappa {kap[fam]:+.2f})")
    ha, _, fh = GAL_HIER[fam]
    L.append(f"  {fam:>5} galaxy HIER:       {ha:.2f} (f_ML={fh} -- "
             f"a0/f_ML-traded, not scorecard-grade; disclosed only)")
L.append("")
L.append("  [4V comparators: BE binary a0 = 1.37+-0.17 (+1.9s); "
         "galaxy flat rows 1.05/1.00 (+0.1/-0.5s)]")

out = "\n".join(L)
print(out)
with open('data/stage5l_ladder.txt', 'w') as f:
    f.write(out+"\n")
print("\nsaved: data/stage5l_ladder.txt")
