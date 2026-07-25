# -*- coding: utf-8 -*-
"""Stage 7J-y: the within-pair overluminosity correlation — the missing
half of Part A's validation (external-review instrument, round 3).

THE POINT (reviewer's design, accepted): the injection gate (GA2) tests
RECOVERY of companions we put in; it does not test OVER-ATTRIBUTION —
whether non-companion overluminosity (metallicity dispersion, age
spread, differential extinction) is absorbed by the companion component
of the delta-mixture. The bias direction is the one that manufactures
the 7J verdict: over-attribution -> C too low -> f_host too high ->
permissive prior -> alpha -> 0.

THE DISCRIMINANT: wide-pair components are coeval and co-chemical.
Metallicity / age / extinction offsets displace BOTH components'
overluminosity delta in the same direction (common-mode); an unresolved
companion displaces ONE (per-component). So the within-pair correlation
rho(delta1, delta2) across the sample reads the common-mode share
directly. N ~ 14k resolves rho ~ 0.03 (1/sqrt(N)); the reviewer's 0.1
easily.

CONSTRUCTION NOTES (this pipeline's specifics, stated before bars):
  - M_G uses each component's OWN parallax (completeness script lines
    63-69), so parallax error is NOT common-mode here — it is
    independent per-component noise that DILUTES rho. We therefore
    correct for attenuation (subtract the per-star parallax-magnitude
    variance from the denominators) rather than partial anything out.
  - Ridge-error confound: two components at SIMILAR color share the
    ridge's local error -> spurious positive rho at small |dcolor|.
    Guard: rho reported in |dcolor| bins; the wide-|dcolor| bins carry
    the metallicity signal cleanly (metallicity shifts delta at any
    color separation).
  - Correlated companion INCIDENCE (Tokovinin 2014: both-component
    subsystem rate ~4% vs ~0.7% independent) also produces positive
    correlation — but concentrated in the flagged TAIL. Guard: the
    common-mode bar is read on the CORE (|delta| < 0.3 both), where
    single stars dominate; the tail incidence ratio
    P(both < -0.4) / P(< -0.4)^2 is reported separately.

PRE-REGISTERED BARS (committed before execution):
  rho_core (attenuation-corrected, |dcolor| > 0.15 slice as primary):
    <= 0.10          -> C-PASS: common-mode subdominant; C = 0.41
                        stands as measured (external-tension flag
                        remains for the absolute scale).
    0.10 - 0.30      -> C-PARTIAL: quantify the common-mode variance
                        share s = rho_core; inflate the Part A envelope
                        by the implied over-attribution band; the
                        literature-anchored prior gains weight.
    >= 0.30          -> C-FAIL: the bright-tail attribution is suspect
                        at factor level; numbered correction; the
                        literature prior becomes primary and the
                        AMBIGUOUS label at that anchoring with it.
  All rho values quoted with the attenuation correction and raw.
Output: data/stage7j_paircorr.txt
"""
import numpy as np
from astropy.io import fits

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

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

c1, c2 = d['bp_rp1'], d['bp_rp2']
col = np.concatenate([c1[ok], c2[ok]])
mg = np.concatenate([MG1[ok], MG2[ok]])
good = np.isfinite(col) & np.isfinite(mg) & (col > -0.5) & (col < 6.0)
CB = np.linspace(np.nanpercentile(col[good], 0.5),
                 np.nanpercentile(col[good], 99.5), 41)
cc = 0.5*(CB[:-1]+CB[1:])
ridge = np.full(len(cc), np.nan)
for i in range(len(cc)):
    m = good & (col >= CB[i]) & (col < CB[i+1])
    if m.sum() > 30:
        ridge[i] = np.median(mg[m])
vr = np.isfinite(ridge)
ccv, rgv = cc[vr], ridge[vr]
def ridge_f(c):
    return np.interp(c, ccv, rgv)

# per-PAIR deltas (identical ridge to Part A)
g1 = np.isfinite(c1[ok]) & (c1[ok] > -0.5) & (c1[ok] < 6.0)
g2 = np.isfinite(c2[ok]) & (c2[ok] > -0.5) & (c2[ok] < 6.0)
both = g1 & g2
d1 = (MG1[ok] - ridge_f(c1[ok]))[both]
d2 = (MG2[ok] - ridge_f(c2[ok]))[both]
dcol = np.abs(c1[ok] - c2[ok])[both]
# per-star parallax-magnitude noise variance (5/ln10 * eplx/plx)^2
s1 = (5/np.log(10)*(eplx1/np.maximum(plx1, 1e-6)))[ok][both]
s2 = (5/np.log(10)*(eplx2/np.maximum(plx2, 1e-6)))[ok][both]
P(f"STAGE 7J-y: {both.sum()} pairs with both components valid; "
  f"median parallax-mag noise {np.median(s1):.3f}/{np.median(s2):.3f} mag")

def rho_att(a, b, va, vb):
    """Pearson rho, raw and attenuation-corrected (subtracting known
    per-star noise variance from the denominators)."""
    a = a - a.mean(); b = b - b.mean()
    cov = float(np.mean(a*b))
    ra = cov/np.sqrt(np.mean(a*a)*np.mean(b*b))
    den = np.sqrt(max(np.mean(a*a)-np.mean(va), 1e-6)
                  * max(np.mean(b*b)-np.mean(vb), 1e-6))
    return ra, cov/den

# global
ra, rc = rho_att(d1, d2, s1**2, s2**2)
P(f"ALL pairs (N={len(d1)}): rho_raw={ra:+.3f}  rho_att={rc:+.3f}  "
  f"(1/sqrt(N) = {1/np.sqrt(len(d1)):.3f})")

# the CORE bar (single-star-dominated), split by color separation
core = (np.abs(d1) < 0.3) & (np.abs(d2) < 0.3)
for lab, m in (('CORE all', core),
               ('CORE |dcol|<0.15', core & (dcol < 0.15)),
               ('CORE |dcol|>=0.15 [BAR]', core & (dcol >= 0.15)),
               ('CORE |dcol|>=0.4', core & (dcol >= 0.4))):
    if m.sum() > 100:
        ra, rc = rho_att(d1[m], d2[m], s1[m]**2, s2[m]**2)
        P(f"{lab} (N={int(m.sum())}): rho_raw={ra:+.3f}  rho_att={rc:+.3f}")

# tail incidence (companion-incidence correlation channel, reported
# separately from the common-mode bar)
f1 = d1 < -0.4; f2 = d2 < -0.4
pf1, pf2 = f1.mean(), f2.mean()
pboth = (f1 & f2).mean()
P(f"TAIL incidence: P(flag1)={pf1:.4f} P(flag2)={pf2:.4f} "
  f"P(both)={pboth:.4f}  ratio P(both)/(P1*P2)={pboth/max(pf1*pf2,1e-9):.2f}")
# cross-channel: one flagged, other's delta (companion in one should not
# brighten the other; common-mode would)
m1 = f1 & ~f2
P(f"CROSS: median delta2 | (flag1, !flag2) = "
  f"{np.median(d2[m1]):+.3f} mag (N={int(m1.sum())}); "
  f"median delta2 | (!flag1, !flag2) = "
  f"{np.median(d2[~f1 & ~f2]):+.3f} (N={int((~f1 & ~f2).sum())})")

# verdict per pre-registered bars (primary = CORE |dcol|>=0.15, attenuated)
m = core & (dcol >= 0.15)
_, rbar = rho_att(d1[m], d2[m], s1[m]**2, s2[m]**2)
v = ('C-PASS' if rbar <= 0.10 else
     'C-PARTIAL' if rbar < 0.30 else 'C-FAIL')
P(f"\n==> {v} (rho_core_att = {rbar:+.3f} against bars 0.10 / 0.30)")

with open('data/stage7j_paircorr.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("saved: data/stage7j_paircorr.txt")
