"""
STAGE 4H: M/L marginalization of the screening index p (TODO #4).
Stage 1 measured p = 0.443 +0.063/-0.050 with M/L FIXED at (0.5 disk, 0.7
bulge). This closes the debt three ways:
  (a) regression: reproduce the Stage-1 fit exactly (f_d = 1 fixed, seed 42);
  (b) global disk-M/L scale f_d free (profiled with la0, p) + 200-replicate
      galaxy bootstrap -> the M/L-marginalized p and its stat error;
  (c) per-galaxy M/L jitter: f_d,g = f_d_best * 10^N(0, 0.1 dex) (the
      SPARC-quoted ~0.1 dex photometric M/L uncertainty), 200 realizations,
      refit (la0, p) each -> the M/L-scatter systematic on p.
Cassini re-verdict at the widened band. Writes data/stage4h_p_ml.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize

KPC = 3.24078e-14
UD, UB = 0.5, 0.7

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

g_gas, g_dsk, g_bul, gobs, gal_id = [], [], [], [], []
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
        gobs.append(Vo*Vo/R*KPC); gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, gal_id))
lgobs = np.log10(gobs)
NG = gal_id.max()+1
uniq = np.unique(gal_id)

def nu_screen(y, p):
    return (1 - np.exp(-np.clip(y, 1e-12, None)**p))**(-1/(2*p))

def loss_gb(gb, la0, p):
    if not (0.05 <= p <= 3): return 1e9
    r = lgobs_sub - np.log10(gb*nu_screen(gb/10**la0, p))
    return float(np.mean(r*r))

L = [f"STAGE 4H: M/L-marginalized screening index ({kept} galaxies, "
     f"{len(gobs)} points)"]

def fit_p(idx, fd, fb=1.0, x0=(-9.92, 0.5)):
    global lgobs_sub
    gb = (g_gas + fd*g_dsk + fb*g_bul)[idx]
    lgobs_sub = lgobs[idx]
    b = minimize(lambda v: loss_gb(gb, v[0], v[1]), list(x0),
                 method='Nelder-Mead',
                 options={'xatol': 1e-4, 'fatol': 1e-10})
    return b.x[0], b.x[1], b.fun

def fit_p_fd(idx, x0=(-9.92, 0.5, 1.0)):
    global lgobs_sub
    lgobs_sub = lgobs[idx]
    def lo(v):
        la0, p, fd = v
        if not (0.05 <= p <= 3) or not (0.3 <= fd <= 3): return 1e9
        gb = (g_gas + fd*g_dsk + g_bul)[idx]
        r = lgobs_sub - np.log10(gb*nu_screen(gb/10**la0, p))
        return float(np.mean(r*r))
    b = minimize(lo, list(x0), method='Nelder-Mead',
                 options={'xatol': 1e-4, 'fatol': 1e-10, 'maxiter': 3000})
    return b.x, b.fun

allidx = np.arange(len(gobs))

# (a) regression: Stage-1 fit + bootstrap, seed 42, f_d = 1
la0_b, p_b, f0 = fit_p(allidx, 1.0)
L.append(f"(a) regression (f_d=1): p = {p_b:.3f}, a0 = {10**la0_b:.3e}, "
         f"rms = {math.sqrt(f0):.4f} dex  (Stage-1 point fit to reproduce)")
rng = np.random.default_rng(42)
ps = []
for _ in range(200):
    pick = rng.choice(uniq, size=len(uniq), replace=True)
    idx = np.where(np.isin(gal_id, pick))[0]
    ps.append(fit_p(idx, 1.0, x0=(la0_b, p_b))[1])
ps = np.array(ps)
L.append(f"    bootstrap: p = {np.median(ps):.3f} "
         f"+{np.percentile(ps,84)-np.median(ps):.3f} / "
         f"-{np.median(ps)-np.percentile(ps,16):.3f}  "
         f"(stored Stage 1: 0.443 +0.063/-0.050)")

# (b) global f_d free
(la0_f, p_f, fd_f), lf = fit_p_fd(allidx)
L.append(f"(b) global f_d free: p = {p_f:.3f}, a0 = {10**la0_f:.3e}, "
         f"f_d = {fd_f:.2f} (M/L_disk = {0.5*fd_f:.2f}), "
         f"rms = {math.sqrt(lf):.4f} dex")
rng2 = np.random.default_rng(7)
ps_f, a0s_f, fds_f = [], [], []
for _ in range(200):
    pick = rng2.choice(uniq, size=len(uniq), replace=True)
    idx = np.where(np.isin(gal_id, pick))[0]
    (la0_, p_, fd_), _ = fit_p_fd(idx, x0=(la0_f, p_f, fd_f))
    ps_f.append(p_); a0s_f.append(10**la0_); fds_f.append(fd_)
ps_f, a0s_f, fds_f = map(np.array, (ps_f, a0s_f, fds_f))
L.append(f"    bootstrap: p = {np.median(ps_f):.3f} "
         f"+{np.percentile(ps_f,84)-np.median(ps_f):.3f} / "
         f"-{np.median(ps_f)-np.percentile(ps_f,16):.3f};  "
         f"a0 = {np.median(a0s_f):.3e} +/- {np.std(a0s_f):.2e};  "
         f"f_d = {np.median(fds_f):.2f} +/- {np.std(fds_f):.2f}")

# (c) per-galaxy M/L jitter (0.1 dex), refit (la0, p) at f_d = fd_f
rng3 = np.random.default_rng(13)
ps_j = []
for _ in range(200):
    dj = rng3.normal(0, 0.1, NG)
    fdg = fd_f*10**dj[gal_id]
    gb = g_gas + fdg*g_dsk + g_bul
    global lgobs_sub
    lgobs_sub = lgobs
    b = minimize(lambda v: (1e9 if not (0.05 <= v[1] <= 3) else float(np.mean(
        (lgobs - np.log10(gb*nu_screen(gb/10**v[0], v[1])))**2))),
        [la0_f, p_f], method='Nelder-Mead',
        options={'xatol': 1e-4, 'fatol': 1e-10})
    ps_j.append(b.x[1])
ps_j = np.array(ps_j)
syst_j = ps_j.std(ddof=1)
shift_j = ps_j.mean() - p_f
L.append(f"(c) per-galaxy 0.1-dex M/L jitter (200 realizations): "
         f"p shift = {shift_j:+.3f}, scatter = {syst_j:.3f}")

# combined statement + Cassini
stat_lo = np.median(ps_f) - np.percentile(ps_f, 16)
stat_hi = np.percentile(ps_f, 84) - np.median(ps_f)
tot_lo = math.sqrt(stat_lo**2 + syst_j**2 + shift_j**2)
tot_hi = math.sqrt(stat_hi**2 + syst_j**2 + shift_j**2)
p_med = np.median(ps_f)
L.append(f"FINAL: p = {p_med:.3f} +{tot_hi:.3f} / -{tot_lo:.3f} "
         f"(stat boot + M/L-jitter syst + jitter-shift, quadrature)")

gsat = 6.674e-11*1.989e30/(9.5*1.496e11)**2
a0c = np.median(a0s_f)
y_sat = gsat/a0c
CASSINI = 2e-14
p_min = None
for p in np.arange(0.02, 1.2, 0.001):
    if gsat*(nu_screen(np.array([y_sat]), p)[0]-1) < CASSINI:
        p_min = p; break
lo16 = p_med - tot_lo
L.append(f"Cassini: requires p > {p_min:.3f}; M/L-marginalized 16th pct = "
         f"{lo16:.3f} -> {'CONSISTENT' if lo16 > p_min else 'TENSION'}")
L.append(f"distance of p=1/2 from marginalized fit: "
         f"{abs(0.5-p_med)/max(tot_hi if p_med<0.5 else tot_lo,1e-9):.1f} sigma")

out = "\n".join(L)
print(out)
with open('data/stage4h_p_ml.txt', 'w') as f_:
    f_.write(out+"\n")
print("\nsaved: data/stage4h_p_ml.txt")
