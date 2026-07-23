"""
STAGE 4V: the oscillator scorecard -- every parameter-free prediction of the
thermal/Planck-oscillator reading against its measured value, plus the
wide-binary alpha translated into an independent a0 determination.

The a0 translation: the binary fit measures alpha multiplying (B-1), the
EFE-suppressed boost. Rescaling a0 -> a0(1+eps) at FIXED physical fields
moves both table arguments: y = g_N/a0 -> y/(1+eps) and e_N = g_ext/a0 ->
e_N/(1+eps). So d ln(B-1)/d ln a0 = kappa = -[dln(B-1)/dlny + dln(B-1)/dln e_N],
evaluated from the solver tables (g1p0/1p2/1p4) and averaged over the actual
deep-sample (6-30 kAU) y-distribution. Then a0_bin = 1.2e-10 * alpha^(1/kappa).

Gates: G1 alpha=1 maps to a0 = 1.2e-10 identically (construction);
G2 kappa > 0 (boost strengthens with a0) and stable between families.
Writes data/stage4v_scorecard.txt.
"""
import math
import numpy as np
from astropy.io import fits

A0_TAB = 1.2e-10

def load_tab(path):
    t = np.load(path)
    y, b = t[0][::-1], t[1][::-1]
    return np.log(y), np.log(np.maximum(b-1.0, 1e-12))

TAGS = {1.0: 'g1p0', 1.2: 'g1p2', 1.4: 'g1p4'}
tabs = {}
for fam in ('simple', 'be'):
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
    dE = (lnB1(fam, e+de, lny) - lnB1(fam, e-de, lny))/(2*de*1.0/e)  # dln e_N
    return -(dy + dE)

# deep-sample y distribution from the actual pairs (stage2c cuts, 6-30 kAU)
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
GKMS = 0.8868  # (km/s)^2 kAU per Msun: g_N = GKMS*M/s^2 in (km/s)^2/kAU
UNIT = 6.685e-9  # 1 (km/s)^2/kAU in m/s^2
yd = (GKMS*Mtot[deep]/s_kau[deep]**2)*UNIT/A0_TAB
lnyd = np.log(np.clip(yd, 1e-3, 50.0))

L = [f"STAGE 4V oscillator scorecard; deep pairs (6-30 kAU): {int(deep.sum())}, "
     f"y 16/50/84 = {np.percentile(yd,[16,50,84]).round(3).tolist()}", ""]

kap = {}
for fam in ('simple', 'be'):
    kv = kappa(fam, lnyd)
    kap[fam] = float(np.mean(kv))
    L.append(f"kappa({fam}) = dln(B-1)/dln a0: mean {np.mean(kv):+.3f} "
             f"(16/84 pct {np.percentile(kv,16):+.3f}/{np.percentile(kv,84):+.3f})")
g2 = all(k > 0 for k in kap.values()) and abs(kap['simple']-kap['be']) < 0.5
L.append(f"G2 kappa positive + family-stable -> {'PASS' if g2 else 'FAIL'}")

def a0_of_alpha(alpha, salpha, fam):
    k = kap[fam]
    a0 = A0_TAB*alpha**(1.0/k)
    return a0, a0*salpha/(alpha*k)

a0s, sa0s = a0_of_alpha(1.18, 0.11, 'simple')
a0b, sa0b = a0_of_alpha(1.13, 0.13, 'be')
g1s, _ = a0_of_alpha(1.00, 0.0, 'simple')
L.append(f"G1 alpha=1 -> a0 = {g1s:.3e} (must be 1.200e-10) -> "
         f"{'PASS' if abs(g1s-A0_TAB) < 1e-14 else 'FAIL'}")
L.append("")
L.append(f"a0(binaries, simple) = ({a0s*1e10:.2f} +- {sa0s*1e10:.2f})e-10")
L.append(f"a0(binaries, BE)     = ({a0b*1e10:.2f} +- {sa0b*1e10:.2f})e-10")

# the temperature prediction
C = 2.998e8
def cho2pi(H0, sH0):
    h = H0*1e3/3.0857e22
    v = C*h/(2*math.pi)
    return v, v*sH0/H0
tp, stp = cho2pi(67.4, 0.5)
ts, sts = cho2pi(73.0, 1.0)
L += ["", f"cH0/2pi (Planck 67.4+-0.5)  = {tp*1e10:.3f}e-10",
      f"cH0/2pi (SH0ES 73.0+-1.0)   = {ts*1e10:.3f}e-10", ""]

rows = [
    ("a0 [1e-10 m/s^2]  SPARC p-fit (4H)", 1.05, 0.10),
    ("a0                joint SPARC+lensing (4E)", 1.00, 0.09),
    ("a0                wide binaries (simple)", a0s*1e10, sa0s*1e10),
    ("a0                wide binaries (BE)", a0b*1e10, sa0b*1e10),
]
L.append("THE OSCILLATOR SCORECARD")
L.append(f"{'quantity':<46}{'measured':>16}{'predicted':>12}{'pull':>7}")
for name, v, s in rows:
    pull = (v-tp*1e10)/math.hypot(s, stp*1e10)
    L.append(f"{name:<46}{v:>10.2f}+-{s:<5.2f}{tp*1e10:>10.2f}{pull:>+7.1f}s")
L.append(f"{'(same rows vs SH0ES cH0/2pi = %.2f)' % (ts*1e10):<46}"
         f"{'pulls:':>16}  "
         + " ".join(f"{(v-ts*1e10)/math.hypot(s, sts*1e10):+.1f}s"
                    for _, v, s in rows))
L.append("")
dials = [
    ("screening index p (Wien tail; 4H)", 0.578, 0.118, 0.5),
    ("zero-point coefficient c1 (4S bootstrap)", 0.427, 0.268, 0.5),
    ("binary amplitude alpha (simple; 3U/3V)", 1.18, 0.11, 1.0),
    ("binary amplitude alpha (BE; 3U/3V)", 1.13, 0.13, 1.0),
]
for name, v, s, pred in dials:
    L.append(f"{name:<46}{v:>10.3f}+-{s:<5.3f}{pred:>10.3f}"
             f"{(v-pred)/s:>+7.1f}s")
L.append(f"{'ceiling edge sqrt(2B) (4J; P of fit)':<46}"
         f"{'1.65 obs':>16}{'1.65':>10}{'P=0.62':>9}")
L.append(f"{'rung 2: c2 (1/12 vs 1/8)':<46}{'unresolved':>16}{'1/12':>10}"
         f"{'0.1s':>8}")
L.append(f"{'mode count N (4T/4U; no prediction yet)':<46}"
         f"{'~21 (constraint)':>16}{'--':>10}{'--':>7}")
L.append("")
L.append("note: a0(binaries) assumes the alpha->a0 mapping through the EFE")
L.append("tables (kappa above); it is a TRANSLATION of the disclosed alpha>1")
L.append("lean, not an independent dataset, and shares its systematics.")

out = "\n".join(L)
print(out)
with open('data/stage4v_scorecard.txt', 'w') as f:
    f.write(out+"\n")
