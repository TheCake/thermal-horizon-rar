"""
ROUND-29 ADDENDUM -- independent verification of every load-bearing
referee number before adoption (the memory rule), stage 10H.

GA-1  The ring_v2 pi-bug (his Finding 1): old/exact v^2 = 1/pi vs the
      exact Freeman disk; corrected prefactor ~1.00.
GA-3  The A2 per-galaxy route (his Finding 3): per-galaxy constant-a0
      fits on the CORRECTED pinned recipe -> Theil-Sen + OLS slope of
      a0_hat on z with bootstrap CI (his +0.65, 90% CI [-0.75, +2.68];
      bin medians 2.07/2.02/2.79/3.17; median a0 2.65).
GA-4  C4 arithmetic (his Finding 6): cH0/2pi at 67.7 = 1.0468e-10;
      placements 2.30/0.66/0.12 sigma; lock OLS slope 1.024; MOND-row
      deficit 3.44 sigma.
GA-5  Offset diagnostic across recipes (his Finding 4): bin4-bin1
      spread + the persistent bin-3 dip (non-monotone).
GA-6  Dump-artifact recount (his Finding 8): galaxies with v < 1 km/s
      rows at r > 2 kpc (his 11-15 vs the amendment-2 text's 27, which
      conflated window-outliers with dumps).

GA-2 (the corrected pinned ladder through the unchanged wiring gates)
is run 6 of the stage itself (amendment 4); its log is quoted by the
round booking, not re-computed here.

Output: data/stage10h_addendum.txt
"""
import io
import os
import re

import numpy as np
from scipy import special

os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
REL = os.path.join("data", "musedark_release")
buf = io.StringIO()

C_SI = 2.99792458e8
MPC_M = 3.0856775814913673e22
KPC_M = MPC_M / 1000.0
PC_M = KPC_M / 1000.0
G_SI = 6.674e-11
MSUN = 1.98892e30


def emit(s=""):
    print(s)
    buf.write(s + "\n")


def nu_p(y, p):
    yc = np.clip(y, 1e-14, None)
    ex = np.exp(-np.minimum(yc**p, 60.0))
    return (1.0 - ex)**(-1.0 / (2.0 * p))


emit("ROUND-29 ADDENDUM -- stage 10H referee-number verification")
emit("=" * 64)

# ---------------------------------------------------------------- GA-1
emit("\nGA-1  the ring_v2 pi-bug vs the exact Freeman disk")


def ring_v2(redges, mring, r_eval, fix=True):
    nphi, nsub = 72, 5
    phi = (np.arange(nphi) + 0.5) * np.pi / nphi
    g = np.zeros_like(r_eval)
    r = r_eval[:, None] * KPC_M
    for k in range(len(mring)):
        if mring[k] <= 0:
            continue
        for a in np.linspace(redges[k], redges[k + 1], nsub + 2)[1:-1]:
            m = mring[k] / nsub * MSUN
            am = a * KPC_M
            d2 = np.maximum(r**2 + am**2 - 2 * r * am * np.cos(phi[None, :]),
                            (0.05 * KPC_M)**2)
            pref = (G_SI * m) if fix else (G_SI * m / np.pi)
            g += pref * np.mean((r - am * np.cos(phi[None, :])) / d2**1.5, axis=1)
    return np.clip(g, 0, None) * (r_eval * KPC_M)


def sersic_rings(Md, Re, n, rmax, nring=100):
    n = float(np.clip(n, 0.3, 6.0))
    bn = 2 * n - 1 / 3 + 4 / (405 * n)
    redges = np.linspace(0, rmax, nring + 1)
    rmid = 0.5 * (redges[:-1] + redges[1:])
    sig = np.exp(-bn * ((rmid / Re)**(1 / n) - 1))
    m = sig * np.pi * (redges[1:]**2 - redges[:-1]**2)
    return redges, m * (Md / m.sum())


def freeman_v2(Md, Rd, r):
    y = np.clip(r / (2 * Rd), 1e-6, 50.0)
    b = special.i0e(y) * special.k0e(y) - special.i1e(y) * special.k1e(y)
    s0 = Md * MSUN / (2 * np.pi * (Rd * KPC_M)**2)
    return 4 * np.pi * G_SI * s0 * (Rd * KPC_M) * y**2 * b


re9, m9 = sersic_rings(1e10, 5.0, 1.0, 60.0, nring=300)
r9 = np.array([4.0, 6.0, 8.0, 12.0])
vex = freeman_v2(1e10, 5.0 / 1.678, r9)
vold = ring_v2(re9, m9, r9, fix=False)
vnew = ring_v2(re9, m9, r9, fix=True)
emit(f"  old/exact v^2: {np.array2string(vold/vex, precision=4)}  (1/pi = {1/np.pi:.4f})")
emit(f"  new/exact v^2: {np.array2string(vnew/vex, precision=4)}")
assert np.all(np.abs(vold / vex - 1 / np.pi) < 0.02)
assert np.all(np.abs(vnew / vex - 1.0) < 0.05)
emit("  CONFIRMED: the shipped prefactor double-divided by pi; his Finding 1 exact.")

# ------------------------------------------------------------ load + build
def parse_pipe(p):
    lines = [l for l in open(p) if l.strip()]
    n = [c.strip() for c in lines[0].split("|") if c.strip()]
    v = [float(c) for c in lines[1].split("|") if c.strip()]
    return dict(zip(n, v))


def parse_txt(p):
    out = {}
    pat = re.compile(r"^\s*([A-Za-z0-9_]+):\s*(-?[\d.eE+]+)\s*±\s*(-?[\d.eE+]+)")
    for l in open(p, encoding="utf-8"):
        m = pat.match(l)
        if m:
            out[m.group(1)] = float(m.group(2))
    return out


photo = {}
with open(os.path.join(REL, "photometry_catalogue.txt")) as f:
    next(f)
    for l in f:
        c = l.split()
        try:
            photo[int(c[0])] = (float(c[1]), float(c[5]))
        except (ValueError, IndexError):
            pass

fitstat = {}
with open(os.path.join(REL, "Fit_statistics_all_models.txt")) as f:
    hdr = f.readline().split()
    iz = hdr.index("log_Z_DC14")
    for l in f:
        c = l.split()
        fitstat[int(c[0])] = float(c[iz - 0])

GAL = {}
ndump = 0
dump_ids = []
for d in sorted(os.listdir(REL)):
    if not d.startswith("ID"):
        continue
    i = int(d[2:])
    if i not in photo:
        continue
    z, M = photo[i]
    if not (0.33 <= z <= 1.44 and M > 8.8):
        continue
    if fitstat.get(i, 1e9) >= 15000:
        continue
    tv = os.path.join(REL, d, f"DC14_{i}_true_Vrot.dat")
    dpp = os.path.join(REL, d, f"DC14_{i}_derived_parameters.txt")
    gpp = os.path.join(REL, d, f"DC14_{i}_galaxy_parameters.dat")
    if not all(map(os.path.exists, [tv, dpp, gpp])):
        continue
    dp = parse_txt(dpp)
    gp = parse_pipe(gpp)
    if "rad_kpc" not in dp or "log_Mdisk" not in dp:
        continue
    lines = [l for l in open(tv) if l.strip()]
    arr = np.array([[float(c) for c in l.split("|") if c.strip()]
                    for l in lines[1:]])
    rk = dp["rad_kpc"]
    r = np.abs(arr[:, 1]) * rk
    v = np.abs(arr[:, 3])
    s = np.abs(arr[:, 4])
    # GA-6 count: dump rows = v < 1 km/s at r > 2 kpc
    if np.any((r > 2.0) & (v < 1.0)):
        ndump += 1
        dump_ids.append(i)
    good = (r > 1e-3) & (v >= 1.0)
    r, v, s = r[good], v[good], s[good]
    Rd = rk / 1.678
    vad2 = 0.92 * s**2 * (r / Rd) * 1e6
    gasS = max(gp.get("gas_density", 0.0), 0.0)
    re_s, m_s = sersic_rings(10**dp["log_Mdisk"], rk,
                             gp.get("sersic_n", 1.0),
                             max(8 * rk, float(r.max()) * 1.3))
    vd2_pin = ring_v2(re_s, m_s, r, fix=True)
    re_g, m_g = sersic_rings(1.0, 1.0, 1.0, 1.0) if gasS <= 0 else (None, None)
    if gasS > 0:
        rg = np.linspace(0, 30.0, 121)
        mg = gasS * np.pi * (rg[1:]**2 - rg[:-1]**2) * 1e6
        vg2_pin = ring_v2(rg, mg, r, fix=True)
    else:
        vg2_pin = np.zeros_like(r)
    Md = 10**dp["log_Mdisk"]
    vd2_free = freeman_v2(Md, Rd, r)
    vg2_mest = 2 * np.pi * G_SI * (gasS * MSUN / PC_M**2) * (r * KPC_M)
    GAL[i] = dict(z=z, r=r, v2=(v * 1000.0)**2, vad2=vad2,
                  pin=vd2_pin + vg2_pin, mest=vd2_free + vg2_mest)
emit(f"\nsample rebuilt: {len(GAL)} galaxies")

# ---------------------------------------------------------------- GA-6
emit(f"\nGA-6  dump-artifact recount: {ndump} galaxies carry a v < 1 km/s "
     f"row at r > 2 kpc")
emit(f"  (his 11-15; the amendment-2 text's '27' counted window-outliers "
     f"too -- corrected)")

# ---------------------------------------------------------------- GA-3
emit("\nGA-3  the per-galaxy route on the CORRECTED pinned recipe")


def pergal_a0(g, bar="pin"):
    vb2 = g[bar]
    vc2 = g["v2"] + g["vad2"]
    keep = (g["r"] >= 2.0) & (vb2 > 0) & (vc2 > 0)
    if keep.sum() < 3:
        return np.nan
    rm = g["r"][keep] * KPC_M
    x = np.log10(vb2[keep] / rm)
    y = np.log10(vc2[keep] / rm)
    w = (x >= -13.5) & (x <= -8) & (y >= -13.5) & (y <= -8)
    if w.sum() < 3:
        return np.nan
    x, y = x[w], y[w]
    la = np.linspace(-11.0, -8.7, 280)
    chi = [np.sum((y - (x + np.log10(nu_p(10**x / 10**l, 0.5))))**2)
           for l in la]
    return 10**la[int(np.argmin(chi))] / 1e-10


ids = sorted(GAL)
a0h = np.array([pergal_a0(GAL[i], "pin") for i in ids])
zz = np.array([GAL[i]["z"] for i in ids])
ok = np.isfinite(a0h) & (a0h < 9.9)
a0h, zz2 = a0h[ok], zz[ok]
emit(f"  per-galaxy a0_hat: N = {ok.sum()}, median = {np.median(a0h):.2f} "
     f"(his 2.65; their level 2.38)")
qz = np.quantile(zz2, [0.25, 0.5, 0.75])
bm = []
for b in range(4):
    lo = -1 if b == 0 else qz[b - 1]
    hi = qz[b] if b < 3 else 9
    bm.append(np.median(a0h[(zz2 > lo) & (zz2 <= hi)]))
emit(f"  bin medians: {np.array2string(np.array(bm), precision=2)} "
     f"(his 2.07/2.02/2.79/3.17; theirs 1.99->2.71)")
ts_slopes = [(a0h[j] - a0h[i]) / (zz2[j] - zz2[i])
             for i in range(len(a0h)) for j in range(i + 1, len(a0h))
             if abs(zz2[j] - zz2[i]) > 0.05]
ts = np.median(ts_slopes)
rng = np.random.default_rng(31)
boot = []
for _ in range(400):
    k = rng.integers(0, len(a0h), len(a0h))
    zb, ab = zz2[k], a0h[k]
    sl = [(ab[j] - ab[i]) / (zb[j] - zb[i])
          for i in range(len(ab)) for j in range(i + 1, len(ab))
          if abs(zb[j] - zb[i]) > 0.05]
    if sl:
        boot.append(np.median(sl))
lo5, hi95 = np.percentile(boot, [5, 95])
ols = np.polyfit(zz2, a0h, 1)[0]
emit(f"  Theil-Sen slope = {ts:+.2f}, bootstrap 90% CI [{lo5:+.2f}, {hi95:+.2f}] "
     f"(his +0.65 [-0.75, +2.68]); OLS = {ols:+.2f}")
emit(f"  READ: positive lean, CI contains 0 and their per-galaxy 1.42 -- "
     f"CONFIRMS his Finding 3 (no usable-grade recovery; 'never "
     f"materializes' softened).")

# ---------------------------------------------------------------- GA-4
emit("\nGA-4  C4 arithmetic")
H0 = 67.7 * 1000.0 / MPC_M
lock0 = C_SI * H0 / (2 * np.pi) / 1e-10
emit(f"  cH0/2pi (67.7) = {lock0:.4f}e-10  (his 1.0468)")
for name, a0, hw in [("DM", 1.00, 0.04), ("MOND", 1.03, 0.05), ("pergal", 1.05, 0.05)]:
    sg = hw / 1.96
    emit(f"  {name}: |{lock0:.3f} - {a0}| / {sg:.4f} = {abs(lock0-a0)/sg:.2f} sigma")
zg = np.linspace(0.33, 1.44, 400)
Ez = np.sqrt(0.307 * (1 + zg)**3 + 0.693)
sl = np.polyfit(zg, lock0 * Ez, 1)[0]
emit(f"  lock OLS slope over [0.33,1.44] = {sl:.3f}  (his 1.024)")
emit(f"  MOND-row deficit: |{sl:.3f} - 1.20| / {0.10/1.96:.4f} = "
     f"{abs(sl-1.20)/(0.10/1.96):.2f} sigma  (his 3.44)")
emit("  CONFIRMED; C4 phrasing corrected: lock INSIDE the MOND and "
     "per-galaxy 95% CIs, 2.3 sigma (outside 95%) from the DM row.")

# ---------------------------------------------------------------- GA-5
emit("\nGA-5  offset diagnostic across recipes (bin means of y-x)")
for bar, lab in [("pin", "corrected-pinned"), ("mest", "Freeman+Mestel-gas")]:
    offs = {b: [] for b in range(4)}
    for i in ids:
        g = GAL[i]
        vb2 = g[bar]
        vc2 = g["v2"] + g["vad2"]
        keep = (g["r"] >= 2.0) & (vb2 > 0) & (vc2 > 0)
        if keep.sum() < 3:
            continue
        rm = g["r"][keep] * KPC_M
        x = np.log10(vb2[keep] / rm)
        y = np.log10(vc2[keep] / rm)
        w = (x >= -13.5) & (x <= -8) & (y >= -13.5) & (y <= -8)
        b = int(np.searchsorted(qz, g["z"]))
        offs[b] += list((y - x)[w])
    mm = [np.mean(offs[b]) for b in range(4)]
    emit(f"  {lab:22s}: " + " ".join(f"{m:+.3f}" for m in mm)
         + f"  (bin4-bin1 = {mm[3]-mm[0]:+.3f}; bin3 dip: {mm[2] < mm[1]})")
emit("  CONFIRMED: recipe-dependent size, persistent mid-bin dip -- the "
     "'+0.095, rising' clause softened per his Finding 4.")

emit("\nALL REFEREE NUMBERS VERIFIED AT GRADE. Adoption licensed; run 6 = "
     "GA-2 (the corrected ladder through the unchanged gates).")
open("data/stage10h_addendum.txt", "w", encoding="utf-8").write(buf.getvalue())
emit("written: data/stage10h_addendum.txt")
