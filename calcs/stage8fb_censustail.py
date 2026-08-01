"""STAGE 8F-b — THE CENSUS TAIL NULL (pre-reg 4cade0b, committed
BEFORE this run; referee round-3 #1 = the D1 analytic instrument).

Can a fat-error-tail Newton world produce the operative census
pair (band = 9 in [sqrt2, 1.67), cliff = 2 in [1.67, 2.2))?
Per-pair mixture-kernel generalization of the 4J T2b leakage
null: every sub-edge pair contributes its own leakage probability
under sigma_i -> {sigma_i, KT*sigma_i} with weights (1-ftl, ftl).
Estimator direction pre-stated: observed-as-truth + corrected-9
are both EXPOSURE-conservative.  Bars B-EXPOSED / B-IMMUNE / GRAY
locked in NOTES; NO credence move in any branch.
Output: data/stage8fb_censustail.txt
"""
import math
import re
import numpy as np
from astropy.io import fits
from math import erf, sqrt, exp, factorial

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

# ---------- data + masks: verbatim 4J construction ----------
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
sig_vt = sigv/np.sqrt(2)/vc
WIDE = ok & (s_kau >= 6)

BAND = (1.414, 1.67)
CLIFF = (1.67, 2.2)
OBS_B, OBS_C = 9, 2          # operative corrected census pair (7I/7K-b)
FTLS = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50]
KTS = [2.0, 4.0, 8.0]

def pwin(x, s, a, b):
    return 0.5*(erf((b-x)/(s*sqrt(2))) - erf((a-x)/(s*sqrt(2))))

def pois_le(nobs, mu):
    if mu <= 0:
        return 1.0
    return sum(exp(-mu)*mu**k/factorial(k) for k in range(nobs+1))

def pois_ge(nobs, mu):
    return 1.0 - pois_le(nobs-1, mu)

def pb_exact(ps, upto):
    """Exact Poisson-binomial: P(N = 0..upto) via DP convolution
    (amendment 1: the Poisson approximation failed G1 at the
    deciding cell — anti-conservative 5x)."""
    dp = np.zeros(upto+2)
    dp[0] = 1.0
    for p in ps:
        dp[1:] = dp[1:]*(1-p) + dp[:-1]*p
        dp[0] *= (1-p)
    return dp          # dp[k] = P(N = k) for k <= upto; dp[upto+1] = overflow-lumped

P("8F-b THE CENSUS TAIL NULL (pre-reg 4cade0b; bars locked before "
  "this run; analytic per-pair mixture kernel)")
P(f"baseline-cut pairs {int(ok.sum())}; wide (s>=6 kAU) "
  f"{int(WIDE.sum())}")
P("")

# ---------- G0: reproduce the published 4J T2b numbers ----------
P("G0 IDENTITY (the 4J T2b median-sigma path at ftl = 0):")
g0_all_ok = True
for sncut, tag in ((0.0, "no S/N cut"), (3.0, "S/N>3")):
    use = WIDE & (sn > sncut)
    top = use & (gam >= 75)
    edges = [(1.20, 1.414), (1.414, 1.67), (1.67, 2.2)]
    cnt = [int((top & (vt >= a) & (vt < b)).sum()) for a, b in edges]
    sel = top & (vt >= 1.0) & (vt < 2.2)
    sbar = float(np.median(sig_vt[sel])) if int(sel.sum()) >= 5 else 0.08
    rho = cnt[0]/(1.414-1.20)
    def leak(edge, a, b, rho_, s_):
        xs = np.linspace(1.0, edge, 400)
        Pab = 0.5*(np.array([erf((b-x)/(s_*sqrt(2)))-erf((a-x)/(s_*sqrt(2)))
                             for x in xs]))
        return float(rho_*np.trapezoid(Pab, xs))
    for edge, ename in ((1.414, "sqrt2"), (1.65, "boosted")):
        mu1 = leak(edge, *BAND, rho, sbar)
        mu2 = leak(edge, *CLIFF, rho, sbar)
        p1 = pois_ge(cnt[1], mu1)
        p2 = pois_le(cnt[2], mu2)
        P(f"  [{tag}] {ename}: sbar={sbar:.3f} rho={rho:.1f} "
          f"mu_band={mu1:.2f} (obs {cnt[1]}, P>= {p1:.2e}) "
          f"mu_cliff={mu2:.2f} (obs {cnt[2]}, P<= {p2:.2f})")
try:
    ref = open('data/stage4j_gamma82.txt', encoding='utf-8',
               errors='replace').read()
    refln = [ln for ln in ref.splitlines() if 'T2b' in ln]
    P("  reference 4J output lines (verbatim, for the diff):")
    for ln in refln:
        P("   | " + ln.strip())
    P("  G0: compare the recomputed rows against the reference "
      "lines above — the medians/mus/Ps must match at printed "
      "precision (visual certificate; the code path is verbatim).")
except OSError:
    P("  G0: data/stage4j_gamma82.txt absent — recomputed only, "
      "SKIPPED-disclosed")
P("")

# ---------- the per-pair D1 estimator ----------
use = WIDE & (sn > 3.0)                 # primary census convention
top = use & (gam >= 75)
src = top & (vt >= 0.2) & (vt < 1.414)  # sub-edge sources, Newton truth
xv = vt[src]; sv = sig_vt[src]
P(f"per-pair estimator (primary, S/N>3): {int(src.sum())} sub-edge "
  f"source pairs (0.2 <= vt < 1.414, gamma >= 75); observed census "
  f"pair = (band {OBS_B}, cliff {OBS_C})")
P("")
P("KT   ftl   mu_band  mu_cliff  ratio   Pois(Nb>=9)  EXACT(Nb>=9)  "
  "EXACT(Nc<=2)  P_joint(EXACT)")
res = {}
for KT in KTS:
    for ftl in FTLS:
        pb = (1-ftl)*np.array([pwin(x, s, *BAND) for x, s in zip(xv, sv)]) \
           + ftl*np.array([pwin(x, KT*s, *BAND) for x, s in zip(xv, sv)])
        pc = (1-ftl)*np.array([pwin(x, s, *CLIFF) for x, s in zip(xv, sv)]) \
           + ftl*np.array([pwin(x, KT*s, *CLIFF) for x, s in zip(xv, sv)])
        mb, mc = float(pb.sum()), float(pc.sum())
        PbP = pois_ge(OBS_B, mb)                     # failed-gate record
        dpb = pb_exact(pb, OBS_B-1)
        Pb = float(max(1.0 - dpb[:OBS_B].sum(), 0.0))
        dpc = pb_exact(pc, OBS_C)
        Pc = float(min(dpc[:OBS_C+1].sum(), 1.0))
        res[(KT, ftl)] = (mb, mc, Pb, Pc, Pb*Pc, pb, pc)
        P(f"{KT:.0f}  {ftl:5.2f}  {mb:7.3f}  {mc:8.3f}  "
          f"{(mc/max(mb,1e-12)):5.2f}  {PbP:.3e}  {Pb:.3e}  "
          f"{Pc:.3e}  {Pb*Pc:.3e}")
P("")

# ---------- G1: MC bracket of the EXACT estimator (amendment 1) --
P("G1 MC (per-pair Bernoulli draws, 1e6 reps, KT = 4; brackets the "
  "EXACT Poisson-binomial at resolvable cells):")
rng = np.random.default_rng(31)
for ftl in (0.10, 0.50):
    mb, mc, Pb, Pc, Pj, pb, pc = res[(4.0, ftl)]
    hits_b = 0; hits_c = 0
    for _ in range(10):
        u = rng.random((100000, len(pb)))
        hits_b += int(((u < pb).sum(axis=1) >= OBS_B).sum())
        u2 = rng.random((100000, len(pc)))
        hits_c += int(((u2 < pc).sum(axis=1) <= OBS_C).sum())
    mcPb, mcPc = hits_b/1e6, hits_c/1e6
    err = sqrt(max(hits_b, 1))/1e6
    P(f"  ftl={ftl:.2f}: EXACT P(Nb>=9)={Pb:.3e} vs MC {mcPb:.3e} "
      f"(+-{err:.1e}); EXACT P(Nc<=2)={Pc:.3e} vs MC {mcPc:.3e}")
P("")

# ---------- G2: per-pair vs median-sigma consistency at ftl=0 ----
mb0 = res[(4.0, 0.0)][0]
P(f"G2: per-pair mu_band(ftl=0) = {mb0:.3f} vs the 4J median-sigma "
  f"structure (order-level consistency check; the per-pair form is "
  f"D1's refinement and the operative one)")
P("")

# ---------- diagnostics: no-S/N-cut and no-lower-vt-floor --------
useN = WIDE
topN = useN & (gam >= 75)
srcN = topN & (vt >= 0.2) & (vt < 1.414)
xvN = vt[srcN]; svN = sig_vt[srcN]
for ftl in (0.0, 0.10, 0.50):
    pbN = (1-ftl)*np.array([pwin(x, s, *BAND) for x, s in zip(xvN, svN)]) \
        + ftl*np.array([pwin(x, 4.0*s, *BAND) for x, s in zip(xvN, svN)])
    pcN = (1-ftl)*np.array([pwin(x, s, *CLIFF) for x, s in zip(xvN, svN)]) \
        + ftl*np.array([pwin(x, 4.0*s, *CLIFF) for x, s in zip(xvN, svN)])
    P(f"no-S/N-cut diagnostic KT=4 ftl={ftl:.2f}: mu_band="
      f"{float(pbN.sum()):.3f} mu_cliff={float(pcN.sum()):.3f} "
      f"(N_src={len(xvN)})")
P("")
# amendment-1 sensitivity: drop the registered vt >= 0.2 source
# floor (apo-island pairs can reach the band at KT >= 4 in rare
# fluctuations; INCLUDING them increases leakage = the
# exposure-conservative direction). Verdict stays on the
# registered estimator; bar-crossings here are flagged.
srcA = top & (vt < 1.414)
xvA = vt[srcA]; svA = sig_vt[srcA]
P(f"amendment-1 sensitivity (NO lower vt floor; N_src={len(xvA)}, "
  f"exposure-conservative):")
sens_max = 0.0
for ftl in FTLS:
    pbA = (1-ftl)*np.array([pwin(x, s, *BAND) for x, s in zip(xvA, svA)]) \
        + ftl*np.array([pwin(x, 4.0*s, *BAND) for x, s in zip(xvA, svA)])
    pcA = (1-ftl)*np.array([pwin(x, s, *CLIFF) for x, s in zip(xvA, svA)]) \
        + ftl*np.array([pwin(x, 4.0*s, *CLIFF) for x, s in zip(xvA, svA)])
    dpbA = pb_exact(pbA, OBS_B-1)
    PbA = float(max(1.0 - dpbA[:OBS_B].sum(), 0.0))
    dpcA = pb_exact(pcA, OBS_C)
    PcA = float(min(dpcA[:OBS_C+1].sum(), 1.0))
    sens_max = max(sens_max, PbA*PcA)
    P(f"  KT=4 ftl={ftl:.2f}: mu_band={float(pbA.sum()):.3f} "
      f"mu_cliff={float(pcA.sum()):.3f} EXACT P_joint={PbA*PcA:.3e}")
P(f"  sensitivity max P_joint = {sens_max:.3e} "
  f"({'CROSSES' if sens_max >= 1e-3 else 'below'} the 1e-3 bar - "
  f"flagged per the ship-the-risk-axis rule)")
P("")

# ---------- verdict per the locked bars (KT = 4 primary; EXACT
# estimator per amendment 1) ----------
pj = {ftl: res[(4.0, ftl)][4] for ftl in FTLS}
mx = max(pj.values())
mx_ftl = max(pj, key=pj.get)
match = [f for f in FTLS if res[(4.0, f)][0] >= OBS_B]
mech = ""
if match:
    f0 = min(match)
    mech = (f"; mechanism line: smallest band-matching ftl = {f0:.2f} "
            f"has mu_cliff = {res[(4.0, f0)][1]:.1f} "
            f"({'>=' if res[(4.0, f0)][1] >= 6 else '<'} 6 = 3x obs)")
else:
    mech = ("; no ftl reaches mu_band >= 9 in expectation (max "
            f"{max(res[(4.0, f)][0] for f in FTLS):.1f} vs 9: the "
            "band is unproducible by the tail class in expectation "
            "- the mechanism clause resolves vacuously in the "
            "defense's favor, amendment 1)")
P(f"max EXACT P_joint over ftl (KT=4) = {mx:.3e} at ftl = "
  f"{mx_ftl:.2f}{mech}")
edge_note = (" [statistic still rising at the registered class "
             "boundary ftl = 0.50 - class-scoped, disclosed]"
             if mx_ftl >= 0.49 else "")
if mx >= 1e-3:
    v = ("B-EXPOSED: the census is tail-fakeable at the KT = 4 class "
         "- annotation flips EXPOSED; reviewer informed next round")
elif mx <= 1e-5 and (not match or res[(4.0, min(match))][1] >= 6):
    v = ("B-IMMUNE: the census pair (band 9, cliff 2) is NOT "
         "producible by symmetric variance-inflation tails at any "
         "severity to 50%")
else:
    v = ("GRAY-CARRIED: between the locked bars - the band count is "
         "unreachable in expectation (Poisson-luck residual only) "
         "but the residual exceeds the 1e-5 immunity bar at the "
         "class edge")
P(f"==> 8F-b VERDICT (locked bars, EXACT estimator): {v}{edge_note}")
P("    NO credence move (pre-stated in every branch).")

with open('data/stage8fb_censustail.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8fb_censustail.txt")
