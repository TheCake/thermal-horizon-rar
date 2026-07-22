"""
STAGE 1: measure the screening index of gravity's transition function.

Data: SPARC (Lelli, McGaugh & Schombert 2016) — 175 disk galaxies, Spitzer
photometry + HI/Halpha rotation curves. We build the radial-acceleration
relation (g_obs vs g_bar) with the standard cuts and mass-to-light ratios
(0.5 disk, 0.7 bulge), then fit:

  One-parameter SCREENING FAMILY (new angle):
     nu_p(y) = (1 - exp(-y^p))^(-1/(2p)),   y = g_bar/a0
  limits:  y->0: nu -> y^(-1/2)  (deep-MOND / Tully-Fisher identity, any p)
           y->inf: nu -> 1 + (1/2p) exp(-y^p)   (Newton, screening sharpness p)
  p = 1/2 is EXACTLY the McGaugh-Lelli-Schombert RAR function.
  Small p = soft tail (Cassini-dangerous); large p = hard cutoff.

Also fitted for comparison: classic nu-families (simple, standard, toy).
Cassini enters as an independent bound on p from Saturn-ephemeris limits.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize_scalar, minimize

KPC = 3.24078e-14      # (km/s)^2/kpc -> m/s^2
UD, UB = 0.5, 0.7      # M/L disk, bulge (McGaugh+16)

# ---- galaxy metadata: inclination + quality ----
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
print(f"metadata rows parsed: {len(meta)} (sample: {list(meta.items())[:2]})")

# ---- build RAR points ----
gbar, gobs, sig, gal_id = [], [], [], []
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept_gal = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q = meta.get(name, (0, 3))
    if inc < 30 or q > 2: continue
    kept_gal += 1
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        vb2 = Vg*abs(Vg) + UD*Vd*abs(Vd) + UB*Vb*Vb
        if vb2 <= 0: continue
        gbar.append(vb2/R*KPC); gobs.append(Vo*Vo/R*KPC)
        sig.append(2*eV/Vo/math.log(10))          # err in dex on log g_obs
        gal_id.append(gi)
gbar, gobs, sig, gal_id = map(np.array, (gbar, gobs, sig, gal_id))
print(f"galaxies kept: {kept_gal}, RAR points: {len(gbar)}")
print(f"g_bar span: {gbar.min():.1e} .. {gbar.max():.1e} m/s^2")

# ---- families ----
def nu_screen(y, p): return (1 - np.exp(-np.clip(y, 1e-12, None)**p))**(-1/(2*p))
def g_model(gb, a0, p): return gb*nu_screen(gb/a0, p)
FAMS = {
  'RAR (p=1/2 fixed)': lambda gb, a0: g_model(gb, a0, 0.5),
  'simple nu':  lambda gb, a0: gb*(0.5+np.sqrt(0.25+a0/gb)),
  'standard nu':lambda gb, a0: gb*np.sqrt((1+np.sqrt(1+4*(a0/gb)**2))/2),
  'toy sqrt':   lambda gb, a0: np.sqrt(gb*gb + a0*gb),
}
def rms_dex(fun, a0):
    r = np.log10(gobs) - np.log10(fun(gbar, a0))
    return math.sqrt(np.mean(r*r))
print(f"\n{'family':>18} {'best a0 [m/s^2]':>16} {'rms scatter':>12} {'dAIC':>8}")
results = {}
for name, fun in FAMS.items():
    f = minimize_scalar(lambda la0: rms_dex(fun, 10**la0), bounds=(-10.5, -9.5),
                        method='bounded')
    a0b, s = 10**f.x, f.fun
    results[name] = (a0b, s)
n = len(gbar)
aic0 = min(n*math.log(results[k][1]**2) + 4 for k in results)
for name in FAMS:
    a0b, s = results[name]
    aic = n*math.log(s*s) + 4
    print(f"{name:>18} {a0b:>16.3e} {s:>10.4f}dx {aic-aic0:>8.1f}")

# ---- the headline: free screening index ----
def loss(v):
    la0, p = v
    if not (0.05 <= p <= 3): return 1e9
    r = np.log10(gobs) - np.log10(g_model(gbar, 10**la0, p))
    return float(np.mean(r*r))
best = minimize(loss, [-9.92, 0.5], method='Nelder-Mead',
                options={'xatol':1e-4,'fatol':1e-10})
la0_b, p_b = best.x
print(f"\nFREE SCREENING INDEX FIT: p = {p_b:.3f}, a0 = {10**la0_b:.3e} m/s^2, "
      f"rms = {math.sqrt(best.fun):.4f} dex")

# galaxy-level bootstrap for p uncertainty
rng = np.random.default_rng(42)
uniq = np.unique(gal_id)
ps, a0s = [], []
for _ in range(200):
    pick = rng.choice(uniq, size=len(uniq), replace=True)
    m = np.isin(gal_id, pick)
    idx = np.where(m)[0]
    gb_, go_ = gbar[idx], gobs[idx]
    def loss_b(v):
        la0, p = v
        if not (0.05 <= p <= 3): return 1e9
        r = np.log10(go_) - np.log10(g_model(gb_, 10**la0, p))
        return float(np.mean(r*r))
    r = minimize(loss_b, [la0_b, p_b], method='Nelder-Mead',
                 options={'xatol':1e-3,'fatol':1e-9,'maxiter':300})
    ps.append(r.x[1]); a0s.append(10**r.x[0])
ps, a0s = np.array(ps), np.array(a0s)
print(f"bootstrap (200 galaxy resamples): p = {np.median(ps):.3f} "
      f"+{np.percentile(ps,84)-np.median(ps):.3f} / "
      f"-{np.median(ps)-np.percentile(ps,16):.3f}")
print(f"                                 a0 = {np.median(a0s):.3e} "
      f"+/- {np.std(a0s):.2e} m/s^2")

# ---- Cassini bound on p (independent) ----
gsat = 6.674e-11*1.989e30/(9.5*1.496e11)**2
a0c = np.median(a0s)
y_sat = gsat/a0c
def dg_sat(p): return gsat*(nu_screen(np.array([y_sat]), p)[0]-1)
CASSINI = 2e-14                       # m/s^2, conservative ephemeris bound
p_min = None
for p in np.arange(0.02, 1.2, 0.001):
    if dg_sat(p) < CASSINI: p_min = p; break
print(f"\nCassini test: y_Saturn = {y_sat:.2e}; anomalous accel < {CASSINI:.0e} m/s^2")
print(f"  requires p > {p_min:.3f}")
print(f"  data give  p = {np.median(ps):.3f} (68%: "
      f"{np.percentile(ps,16):.3f}-{np.percentile(ps,84):.3f})")
lo = np.percentile(ps, 16)
verdict = "CONSISTENT — galaxies choose a Cassini-safe screening index" \
    if lo > p_min else "TENSION — soft tails preferred, Saturn objects"
print(f"  VERDICT: {verdict}")
print(f"  distance of p=1/2 (instanton/RAR form) from best fit: "
      f"{abs(0.5-np.median(ps))/max(np.std(ps),1e-9):.1f} sigma")
