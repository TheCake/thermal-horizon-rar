"""
STAGE 4J: data-side autopsy of the gamma~82 perpendicular excess (TODO #2e).
The excess: top-gamma cells (gamma >~ 75-90 deg) at TWO vtilde islands
(~0.07 and ~1.66), wide s-bins, z ~ +5 against every model variant.

Kinematic identities that frame the test (vis-viva at gamma_3D = 90):
  apocenter:  vt^2 = 1 - e  -> vt = 0.07  <=>  e = 0.995 (slow apo dwell)
  pericenter: vt^2 = 1 + e  -> Newtonian ceiling sqrt(2) = 1.414;
              boosted-law ceiling sqrt(2*boost) ~ sqrt(2*1.36) = 1.65;
              vt = 1.66 sits AT the boosted ceiling, ABOVE the Newtonian one.
  hyperbolic closest approach: gamma = 90 exactly, vt unbounded above sqrt(2).
So the islands are candidate APO and PERI faces of the near-parabolic
population (e-ceiling in the model: 0.995) -- or, rivals: the apo island as a
NOISE FLOOR ridge (vt_meas ~ sigma_v/v_c, which lands at 0.07-0.09 in the
wide bins; true Delta-v ~ 0 comoving pairs or deep-apo orbits), and the peri
island as flybys caught at closest approach (our v7 flyby model puts gamma
in [0,30] -- the asymptote regime -- so a 90-deg closest-approach ridge is
unmodeled BY CONSTRUCTION).

Tests (data only, no GPU):
 T1 gamma-specificity: are the vt islands concentrated at top gamma, or
    present at all gamma (=> not a gamma structure)?  Noise-floor/comoving
    pairs give a FLAT gamma; true apo orientation peaks at 90.
 T2 the ceiling edge: vt distribution of top-gamma pairs across
    [1.20,1.414,1.67,2.2,3.0] -- a sharp drop above ~1.67 favors BOUND
    boosted pericenters; a smooth tail favors unbound flybys.
 T3 implied-e: vt^2 histogram of the slow top-gamma pairs (e = 1 - vt^2).
 T4 quality slices: island vs control medians of S/N, sigma_v, RUWE_max,
    R_chance, distance, |dplx|/sig -- artifact hunting.
 All with and without the S/N>3 direction cut (the apo island is
 noise-floor-adjacent; the cut itself sculpts it).
Writes data/stage4j_gamma82.txt.
"""
import math
import numpy as np
from astropy.io import fits

d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1, plx2 = d['parallax1'], d['parallax2']
eplx1, eplx2 = d['parallax_error1'], d['parallax_error2']
sep, Rch = d['sep_AU'], d['R_chance_align']
G1m, G2m = d['phot_g_mean_mag1'], d['phot_g_mean_mag2']
ruwe1, ruwe2 = d['ruwe1'], d['ruwe2']
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

# gamma (3L construction) and vtilde (2C construction)
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx = dra*np.cos(dec_m); sy = d['dec2']-d['dec1']
vx = d['pmra2']-d['pmra1']; vy = d['pmdec2']-d['pmdec1']
e_vx = np.sqrt(d['pmra_error1']**2+d['pmra_error2']**2)
e_vy = np.sqrt(d['pmdec_error1']**2+d['pmdec_error2']**2)
vmag = np.hypot(vx, vy)
sn = vmag/np.sqrt(0.5*(e_vx**2+e_vy**2))
smag = np.hypot(sx, sy)
cosg = np.abs(sx*vx+sy*vy)/np.maximum(smag*vmag, 1e-12)
gam = np.degrees(np.arccos(np.clip(cosg, 0, 1)))

MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot = np.interp(MG1, MG_T, MS_T)+np.interp(MG2, MG_T, MS_T)
s_kau = sep/1e3
vc = 0.9417*np.sqrt(Mtot/np.maximum(s_kau, 1e-9))
dv = 4.74047/plx*vmag
vt = dv/vc
noise_floor = sigv/np.sqrt(2)/vc          # per-pair expected vt floor scale
rmax = np.maximum(ruwe1, ruwe2)
dplxsig = np.abs(plx1-plx2)/np.hypot(eplx1, eplx2)

L = [f"STAGE 4J gamma~82 autopsy: {int(ok.sum())} pairs after baseline cuts"]
WIDE = ok & (s_kau >= 6)
APO = (vt >= 0.05) & (vt < 0.10)
PERI = (vt >= 1.35) & (vt < 2.05)
MID = (vt >= 0.30) & (vt < 0.90)

for sncut, tag in ((0.0, "no S/N cut"), (3.0, "S/N>3 (3L/v7 convention)")):
    use = WIDE & (sn > sncut)
    L.append("")
    L.append(f"=== {tag}: wide pairs (s>=6 kAU) N={int(use.sum())} ===")
    # T1 gamma-specificity of the islands
    L.append("T1 gamma-profile of each vt window (fraction of the window's "
             "pairs per gamma column; flat = noise/comoving, 90-peaked = orbital):")
    gcols = [(0, 30), (30, 60), (60, 75), (75, 90)]
    for name, wmask in (("apo  0.05-0.10", APO), ("mid  0.30-0.90", MID),
                        ("peri 1.35-2.05", PERI)):
        m = use & wmask
        n = int(m.sum())
        if n < 10:
            L.append(f"  {name}: N={n} (too few)"); continue
        fr = [float(np.mean((gam[m] >= a) & (gam[m] < b))) for a, b in gcols]
        # uniform expectation = column width/90
        exp = [(b-a)/90 for a, b in gcols]
        z = [(fr[i]-exp[i])/max(math.sqrt(exp[i]*(1-exp[i])/n), 1e-9)
             for i in range(4)]
        L.append(f"  {name}: N={n:5d}  " + "  ".join(
            f"[{a}-{b}) {fr[i]:.3f}(z{z[i]:+.1f})" for i, (a, b) in
            enumerate(gcols)))
    # T2 ceiling edge in the top-gamma column
    top = use & (gam >= 75)
    L.append(f"T2 ceiling edge, top-gamma column (gamma>=75, N={int(top.sum())}):"
             f" counts per vt window")
    edges = [(1.20, 1.414), (1.414, 1.67), (1.67, 2.2), (2.2, 3.0), (3.0, 6.0)]
    cnt = [int((top & (vt >= a) & (vt < b)).sum()) for a, b in edges]
    wid = [b-a for a, b in edges]
    dens = [c/w for c, w in zip(cnt, wid)]
    L.append("  " + "  ".join(f"[{a}-{b}): {c} ({dn:.1f}/unit)"
                              for (a, b), c, dn in zip(edges, cnt, dens)))
    L.append(f"  Newtonian-bound ceiling 1.414 | boosted ceiling ~1.67 | "
             f"beyond = unbound/noise only")
    # same edge for the low-gamma comparison column
    lowg = use & (gam < 30)
    cntl = [int((lowg & (vt >= a) & (vt < b)).sum()) for a, b in edges]
    L.append("  gamma<30 comparison: " + "  ".join(
        f"[{a}-{b}): {c}" for (a, b), c in zip(edges, cntl)))
    # T2b leakage null: can measurement noise on bound pairs below a TRUE
    # edge produce the observed counts above it? True density = flat at the
    # observed sub-edge level, convolved with the median per-pair vt error.
    fast = top & (vt >= 1.2) & (vt < 1.414)
    sig_vt = (sigv/np.sqrt(2)/vc)
    sbar = float(np.median(sig_vt[top & (vt >= 1.0) & (vt < 2.2)])) \
        if int((top & (vt >= 1.0) & (vt < 2.2)).sum()) >= 5 else 0.08
    rho = cnt[0]/wid[0]                     # observed density just below sqrt2
    from math import erf, exp, sqrt, factorial
    def leak(edge, a, b, rho_, s_):
        # integral over true vt in [1.0, edge] of rho * P(obs in [a,b])
        xs = np.linspace(1.0, edge, 400)
        Pab = 0.5*(np.array([erf((b-x)/(s_*sqrt(2)))-erf((a-x)/(s_*sqrt(2)))
                             for x in xs]))
        return float(rho_*np.trapezoid(Pab, xs))
    def pois_p(nobs, mu):
        if mu <= 0: return 0.0 if nobs > 0 else 1.0
        return sum(exp(-mu)*mu**k/factorial(k) for k in range(nobs+1))
    for edge, ename in ((1.414, "Newtonian edge sqrt(2)"),
                        (1.65, "boosted edge sqrt(2*1.36)")):
        mu1 = leak(edge, 1.414, 1.67, rho, sbar)
        mu2 = leak(edge, 1.67, 2.2, rho, sbar)
        p1 = 1-pois_p(cnt[1]-1, mu1)        # P(>= observed | leakage)
        p2 = pois_p(cnt[2], mu2)            # P(<= observed | leakage)
        L.append(f"  T2b {ename}: median sigma_vt={sbar:.3f}; expect "
                 f"{mu1:.1f} in [1.414,1.67) vs {cnt[1]} obs "
                 f"(P(>=obs)={p1:.2e}); expect {mu2:.1f} in [1.67,2.2) vs "
                 f"{cnt[2]} obs (P(<=obs)={p2:.2f})")
    # T3 implied e for slow top-gamma pairs
    slow = top & (vt < 0.35)
    if slow.sum() >= 10:
        ei = 1 - vt[slow]**2
        L.append(f"T3 implied e (apo reading) of slow top-gamma pairs "
                 f"(N={int(slow.sum())}): median {np.median(ei):.4f}, "
                 f"16-84% {np.percentile(ei,16):.4f}-{np.percentile(ei,84):.4f}")
        L.append(f"   median vt/noise_floor = "
                 f"{np.median((vt/np.maximum(noise_floor,1e-9))[slow]):.2f} "
                 f"(~1 => island IS the noise floor)")
    # T4 quality slices
    L.append("T4 medians island vs control (top-gamma):")
    rows = (("apo island", top & APO), ("peri island", top & PERI),
            ("mid control", top & MID),
            ("same-vt lowg apo", lowg & APO), ("same-vt lowg peri", lowg & PERI))
    L.append(f"  {'set':>18} {'N':>5} {'S/N':>6} {'sigv':>7} {'RUWEmax':>8} "
             f"{'Rch*1e3':>8} {'dist[pc]':>9} {'dplx/sig':>9} {'s[kAU]':>7}")
    for name, m in rows:
        n = int(m.sum())
        if n < 5:
            L.append(f"  {name:>18} {n:>5} (too few)"); continue
        L.append(f"  {name:>18} {n:>5} {np.median(sn[m]):>6.1f} "
                 f"{np.median(sigv[m]):>7.4f} {np.median(rmax[m]):>8.2f} "
                 f"{1e3*np.median(Rch[m]):>8.3f} "
                 f"{np.median(1000/plx[m]):>9.0f} "
                 f"{np.median(dplxsig[m]):>9.2f} {np.median(s_kau[m]):>7.1f}")

out = "\n".join(L)
print(out)
with open('data/stage4j_gamma82.txt', 'w') as f_:
    f_.write(out+"\n")
print("\nsaved: data/stage4j_gamma82.txt")
