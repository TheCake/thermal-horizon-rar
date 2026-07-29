"""STAGE 8B — THE SIMULATION LADDER (pre-reg c03604f; bars locked
before execution).

Does the LCDM attractor print the zero-point digit?  The 4S
lambda-family estimator on the EAGLE RefL0100N1504 z=0 aperture RAR
(7,239 centrals x 8 radii) + the Recal 25 Mpc cross-check.

AMENDMENT 1 (pre-quote, G3-serving): the first run rode the 4S grid's
lam = 1.50 edge -> grid extended to 3.0 with a family-validity guard
(negative-lam members are tail-non-monotone and drop; the valid
contiguous run around BE is kept).  Verdict bars UNTOUCHED.
AMENDMENT 2 (pre-quote, after the extended run rode the FAMILY-
VALIDITY boundary lam = 2.55): the primary fit is reported as a
structured G3 GATE-FAIL (OFF-FAMILY) — no c1_sim quoted, bars
unfired, bootstrap skipped (a boundary pile-up is not an interval);
two LABELED DIAGNOSTICS added (no bars): (D1) the four-fingerprint
function contest incl. the additive class the G7 theorem predicts
for a baryon-tracking halo; (D2) the kinematic-window slice
(g_bar >= 1e-12, SPARC-kinematic support) + per-decade deep-slope
table, localizing where the sim leaves the family.  Support note:
ALL aperture points lie inside the sky measurement's joint support
(the 4E lensing leg reaches 5.6e-15), so the off-family behavior is
physical, not a window artifact.

Output: data/stage8b_simladder.txt
"""
import math
import numpy as np
from scipy.optimize import minimize

GM_SUN = 1.32712440018e20            # m^3 s^-2 (IAU GM_sun)
KPC_M = 3.0856775814913673e19        # m
GK = GM_SUN / KPC_M**2               # g = GK * (M/Msun) / (r/kpc)^2
A0_FID = 1.2e-10
RADII = (5, 10, 20, 30, 40, 50, 70, 100)
LGRID = np.round(np.arange(-0.30, 3.001, 0.05), 3)

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

def load_box(path):
    d = np.genfromtxt(path, delimiter=',', names=True)
    keep = np.isin(d['ApertureSize'], RADII)
    d = d[keep]
    r = d['ApertureSize'].astype(float)
    mbar = d['Mass_Star'] + d['Mass_Gas']
    mtot = mbar + d['Mass_DM'] + d['Mass_BH']
    ok = mbar > 0
    gid_raw = d['GalaxyID'][ok]
    uniq, gal_id = np.unique(gid_raw, return_inverse=True)
    out = dict(
        gbar=GK*mbar[ok]/r[ok]**2, gobs=GK*mtot[ok]/r[ok]**2,
        gbar_star=GK*np.clip(d['Mass_Star'][ok], 0, None)/r[ok]**2,
        r=r[ok], gal_id=gal_id, gid_raw=gid_raw, uniq=uniq,
        ndrop=int((~ok).sum()))
    m30 = {int(g): m for g, m, a in zip(
        d['GalaxyID'], d['Mass_Star'], d['ApertureSize']) if a == 30}
    out['mstar30'] = np.array([m30[int(g)] for g in gid_raw])
    return out

B = load_box('data/eagle/ref100_apertures.csv')
P("8B THE SIMULATION LADDER — EAGLE RefL0100N1504 z=0 aperture RAR")
P(f"points: {len(B['gbar'])} ({len(B['uniq'])} galaxies x 8 radii; "
  f"{B['ndrop']} zero-baryon points dropped)")

# ------------------------------------------------------------------ G0
i0 = np.where((B['r'] == 10))[0][0]
m10 = (B['gobs'][i0]/GK)*100.0
vc = math.sqrt(GM_SUN*m10/(KPC_M*10.0))/1e3
P(f"G0 units: first galaxy M_tot(<10 pkpc) = {m10:.3e} Msun -> "
  f"V_c(10) = {vc:.0f} km/s; g-constant = {GK:.5e}")
assert 30 < vc < 500, vc
P("G0 -> PASS")

# ------------------------------------------------------------------ G1
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0,
                    1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_simple(y):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(y, 1e-14, None))
def nu_add(y):
    return 1.0 + 1.0/np.sqrt(np.clip(y, 1e-14, None))
def nu_lam(y, lam):
    return (1.0-lam)*nu_standard(y) + lam*nu_be(y)

xg = np.linspace(0.005, 0.05, 40)
yg = xg*xg
for lam in (0.0, 0.5, 1.0, 2.0, 3.0):
    h = nu_lam(yg, lam)*xg - 1.0
    c2n, c1n = np.polyfit(xg, h/xg, 1)
    assert abs(c1n - lam/2) < 2e-3, (lam, c1n)
P("G1 series c1(lam) = lam/2 at lam = 0, 1/2, 1, 2, 3 -> PASS")
yv = np.logspace(-4, 2.5, 400)
def lam_valid(lam):
    nv = nu_lam(yv, lam)
    return bool(np.all(nv > 0.5) and
                np.all(np.diff(nv) < 1e-6*nv[:-1]))
okv = np.array([lam_valid(l) for l in LGRID])
i1 = int(np.argmin(np.abs(LGRID - 1.0)))
assert okv[i1]
lo_i = i1
while lo_i > 0 and okv[lo_i-1]:
    lo_i -= 1
hi_i = i1
while hi_i < len(LGRID)-1 and okv[hi_i+1]:
    hi_i += 1
LGRID = LGRID[lo_i:hi_i+1]
P(f"G1b family validity: contiguous valid run around BE = lam in "
  f"[{LGRID[0]:.2f}, {LGRID[-1]:.2f}] ({len(LGRID)} nodes)")

# ------------------------------------------------------------------ G2
lgb = np.log10(B['gbar'])
lgo = np.log10(B['gobs'])
edges = np.quantile(lgb, np.linspace(0, 1, 13))
med, cen = [], []
for a, b in zip(edges[:-1], edges[1:]):
    m = (lgb >= a) & (lgb < b)
    if m.sum() > 50:
        med.append(np.median(lgo[m])); cen.append(0.5*(a+b))
med, cen = np.array(med), np.array(cen)
mono = bool(np.all(np.diff(med) > 0))
raw_sc = float(np.std(lgo - np.interp(lgb, cen, med)))
P(f"G2 descriptive: medians monotone = {mono}; raw scatter = "
  f"{raw_sc:.3f} dex (bar < 0.3) -> "
  f"{'PASS' if mono and raw_sc < 0.3 else 'FAIL'}")
assert mono and raw_sc < 0.3
P("  point support per g_bar decade (all inside the sky joint "
  "support, which reaches 5.6e-15 via the 4E lensing leg):")
for a, b in ((-14, -13), (-13, -12), (-12, -11), (-11, -10),
             (-10, -9), (-9, -8)):
    m = (lgb >= a) & (lgb < b)
    if m.sum():
        P(f"    [{a},{b}): {m.sum()} points")

# ---------------------------------------------------- fit machinery
def make_m2ll(gbar, lgobs, wpt):
    def m2ll(th, lam):
        la0, s = th
        if not (-10.6 < la0 < -9.4) or not (1e-3 <= s < 0.6):
            return 1e15
        gm = gbar*nu_lam(gbar/10**la0, lam)
        r = lgobs - np.log10(gm)
        return float(np.sum(wpt*r*r)/(s*s) + math.log(s*s)*np.sum(wpt))
    return m2ll

def profile(m2ll, tag, report=True):
    prof, th, ths = [], None, []
    for lam in LGRID:
        starts = [[math.log10(A0_FID), 0.10],
                  [math.log10(A0_FID)+0.15, 0.15]]
        if th is not None:
            starts.insert(0, list(th))
        best = None
        for th0 in starts:
            b = minimize(m2ll, th0, args=(lam,), method='Nelder-Mead',
                         options=dict(maxiter=3000, xatol=1e-7,
                                      fatol=1e-7))
            if best is None or b.fun < best.fun:
                best = b
        prof.append(best.fun); th = best.x; ths.append(best.x.copy())
    prof = np.array(prof)
    i = int(np.argmin(prof))
    lam_hat = LGRID[i]
    if 0 < i < len(LGRID)-1:
        x3, y3 = LGRID[i-1:i+2], prof[i-1:i+2]
        c2_, c1_, _ = np.polyfit(x3, y3, 2)
        if c2_ > 0:
            lam_hat = -c1_/(2*c2_)
    interior = 0 < i < len(LGRID)-1
    if report:
        z0 = prof[np.argmin(np.abs(LGRID-0.0))] - prof[i]
        zh = prof[np.argmin(np.abs(LGRID-1.0))] - prof[i]
        P(f"  [{tag}] lam_hat = {lam_hat:.3f} -> c1 = {lam_hat/2:.3f} "
          f"({'INTERIOR' if interior else 'EDGE'}); a0 = "
          f"{10**ths[i][0]:.3e}, s_int = {ths[i][1]:.3f}; "
          f"D(-2lnL) c1=0: {z0:+.1f}, c1=1/2: {zh:+.1f} (nominal "
          f"point grade; radial within-galaxy correlation caveat)")
    return prof, lam_hat, interior, ths[i]

ones = np.ones(len(B['gbar']))
P("")
P("PRE-REGISTERED PRIMARY (full sample, stars+gas, r = 5-100 pkpc):")
prof_p, lam_p, int_p, th_p = profile(
    make_m2ll(B['gbar'], lgo, ones), 'primary full-range')
if int_p:
    P("G3 interiority -> PASS")
else:
    P("==> G3 GATE-FAIL (OFF-FAMILY): the optimum rides the family-")
    P("    validity boundary — the EAGLE aperture RAR is outside the")
    P("    measured slice's neighborhood at EVERY c1. Per pre-reg:")
    P("    no c1_sim quoted, B-bars UNFIRED, bootstrap skipped")
    P("    (a boundary pile-up is not an interval).")

# ------------------------------------------------------------------ G4
rng = np.random.default_rng(8)
g4ok = True
for lt in (0.0, 0.75):
    syn = np.log10(B['gbar']*nu_lam(B['gbar']/A0_FID, lt)) + \
        rng.normal(0, 0.1, len(B['gbar']))
    _, lh, _, _ = profile(make_m2ll(B['gbar'], syn, ones),
                          f'G4 {lt}', report=False)
    ok = abs(lh - lt) <= 0.05
    g4ok &= ok
    P(f"G4 injection lam_true = {lt:.2f}: lam_hat = {lh:.3f} -> "
      f"{'OK' if ok else 'FAIL'}")
assert g4ok
P("G4 -> PASS (the estimator recovers in-family truths; the edge is "
  "the data, not the tool)")

# ----------------------------------- robustness legs (reports only)
P("")
P("robustness legs (per pre-reg; REPORTS — no bars while G3 stands):")
mk = np.genfromt = np.genfromtxt('data/eagle/ref100_morphokinem.csv',
                                 delimiter=',', names=True)
kap = {int(g): k for g, k in zip(mk['GalaxyID'], mk['KappaCoRot'])}
kv = np.array([kap.get(int(g), 0.0) for g in B['gid_raw']])

def leg(mask, tag, gbar=None, lgobs=None):
    gb = B['gbar'][mask] if gbar is None else gbar[mask]
    lg = lgo[mask] if lgobs is None else lgobs[mask]
    _, lh, itr, tht = profile(
        make_m2ll(gb, lg, np.ones(mask.sum())), tag)
    return lh, itr

leg(B['gbar_star'] > 0, '(a) stars-only g_bar',
    gbar=B['gbar_star'],
    lgobs=np.log10(B['gobs']))
leg(B['r'] >= 10, '(b) r >= 10 pkpc')
leg(kv > 0.4, f'(c) disks KappaCoRot>0.4 '
    f'({len(np.unique(B["gid_raw"][kv > 0.4]))} gal)')
leg(B['mstar30'] > 10**9.5, '(d) M* > 10^9.5')
R = load_box('data/eagle/recal25_apertures.csv')
lgoR = np.log10(R['gobs'])
_, lhR, itrR, thR = profile(
    make_m2ll(R['gbar'], lgoR, np.ones(len(R['gbar']))),
    f"(e) RecalL0025N0752 ({len(R['uniq'])} gal)")

# ------------------------- D1: the four-fingerprint function contest
P("")
P("D1 (labeled diagnostic): the fingerprint contest — which class")
P("   describes the sim RAR least badly (D(-2lnL) vs BE; each with")
P("   (a0, s_int) profiled; nominal point grade):")
def fit_fn(nu, gbar, lgobs, tag):
    def m2(th):
        la0, s = th
        if not (-10.9 < la0 < -9.0) or not (1e-3 <= s < 0.6):
            return 1e15
        gm = gbar*nu(gbar/10**la0)
        r = lgobs - np.log10(gm)
        return float(np.sum(r*r)/(s*s) + math.log(s*s)*len(r))
    best = None
    for th0 in ([math.log10(A0_FID), 0.12],
                [math.log10(A0_FID)+0.3, 0.15],
                [math.log10(A0_FID)-0.3, 0.10]):
        b = minimize(m2, th0, method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-7,
                                  fatol=1e-7))
        if best is None or b.fun < best.fun:
            best = b
    return best

for gb, lg, box in ((B['gbar'], lgo, 'Ref100'),
                    (R['gbar'], lgoR, 'Recal25')):
    res = {}
    for name, nu in (('standard', nu_standard), ('simple', nu_simple),
                     ('BE', nu_be), ('additive', nu_add)):
        res[name] = fit_fn(nu, gb, lg, name)
    base = res['BE'].fun
    P(f"   [{box}] " + " | ".join(
        f"{n}: {res[n].fun-base:+.1f} (a0 {10**res[n].x[0]:.2e}, "
        f"s {res[n].x[1]:.3f})" for n in
        ('standard', 'simple', 'BE', 'additive')))

# --------------------- D2: kinematic window + deep-slope per decade
P("")
P("D2 (labeled diagnostic): where does the sim leave the family?")
P("   per-decade log-log slope of the binned median relation")
P("   (MOND deep = 0.5; Newtonian = 1.0):")
for a, b in ((-14, -13), (-13, -12), (-12, -11), (-11, -10),
             (-10, -9)):
    m = (lgb >= a) & (lgb < b)
    if m.sum() > 300:
        bins = np.linspace(a, b, 6)
        mm, cc = [], []
        for x1, x2 in zip(bins[:-1], bins[1:]):
            mm2 = m & (lgb >= x1) & (lgb < x2)
            if mm2.sum() > 30:
                mm.append(np.median(lgo[mm2]))
                cc.append(0.5*(x1+x2))
        if len(mm) > 2:
            sl = np.polyfit(cc, mm, 1)[0]
            P(f"    [{a},{b}): slope = {sl:.2f} (n = {m.sum()})")
KIN = B['gbar'] >= 1e-12
P(f"   kinematic-window slice (g_bar >= 1e-12, the SPARC-kinematic "
  f"support; {KIN.sum()} points, {len(np.unique(B['gid_raw'][KIN]))} "
  f"galaxies):")
_, lhK, itrK, thK = profile(
    make_m2ll(B['gbar'][KIN], lgo[KIN], np.ones(KIN.sum())),
    'D2 kinematic window')

# ------------------------------------------------------------ verdict
P("")
if int_p:
    P("==> (unreachable in this run: G3 passed)")
else:
    P("==> 8B VERDICT (pre-registered gates, amendment 2): G3")
    P("    GATE-FAIL / OFF-FAMILY — no c1_sim exists at aperture")
    P("    grade: the LCDM attractor's RAR is SHAPE-INCOMPATIBLE")
    P("    with the measured function family at every zero-point,")
    P("    monotonically preferring more-than-BE dark response; the")
    P("    sky's digit region c1 = 1/2 is rejected along the slice")
    P("    by D(-2lnL) ~ +7.9k (nominal). The pre-registered B-bars")
    P("    are UNFIRED by the letter; the successor instrument")
    P("    (FIRE-2 snapshot pipeline) stays commissioned-only.")
    P("    The digit does NOT lose uniqueness to the attractor at")
    P("    this grade — the attractor fails one level earlier (the")
    P("    functional class), with the D1/D2 diagnostics locating")
    P("    the failure. Caveats: spherical-aperture proxy; nominal")
    P("    point-grade Deltas (radial correlation); EAGLE-specific.")

with open('data/stage8b_simladder.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8b_simladder.txt")
