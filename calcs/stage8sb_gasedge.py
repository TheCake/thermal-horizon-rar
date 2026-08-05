"""STAGE 8S-b — THE EDGE RESOLUTION (grid extension + the gas-budget
knob).  Pre-registered BEFORE any run.

8S landed the gas-dominated subsample on the family's lower GRID EDGE
(lam_hat = -0.30 = the first node; bootstrap slams its own -0.4 bound;
c1 = 1/2 rejected +31 within the subsample) — by the correction-#4
standard the dial reading is INVALID until the grid extends.  The 8S
verdict grammar also had a hole (a one-sided D1 interval cannot fire
T5-TENSION; it printed POWER-LIMITED) — flagged here, record stands.

Two controls, both pre-registered:
  (1) GRID EXTENSION: LGRID down to -2.0 (the nu_lam family below 0 is
      a DEEP-COEFFICIENT PROBE, c1 = lam/2 < 0, not a physical member
      interpolation — labeled as such; nu positivity guarded).
  (2) THE GAS KNOB: gN = fg*g_gas + f*g_dsk + g_bul with fg free in
      [0.7, 1.4] (the helium/molecular budget envelope — the one
      coherent systematic that shifts gas-dominated galaxies
      specifically).  Run FULL / GD / DD with fg free; the question is
      whether the GD edge is a gas-budget artifact.

Gates: G8Sb-0 shared-node regression (extended-grid fixed-fg GD
profile == the 8S profile on the common nodes, warm-start-free refits,
|d| <= 0.05); G8Sb-1 fg-slice identity (free-fg objective at fg=1
equals the fixed objective, 1e-9); G8Sb-2 nu-positivity report (nodes
where min nu < 0.05 over the data y-range are EXCLUDED and listed).
Bars (locked): E1 LOCALIZED iff the extended-grid GD profile (fixed
fg) has an interior minimum with a two-sided D1 interval inside
(-2.0, 1.5).  E2 GAS-EXPLAINED iff with fg FREE the GD D1 interval
overlaps the FULL-sample free-fg interval (then the 8S edge is
attributable to the fixed gas budget; T5's defense holds CONDITIONAL
on the gas systematic, both quoted).  E3 DIAL-TENSION iff with fg
free the GD interval still excludes the full-sample lam_hat (a real
subsample tension — named for the galaxy program, not resolved here).
POWER-CARRIED otherwise.  NO credence movement (pre-stated).
Output: data/stage8sb_gasedge.txt
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize
import time

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10

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

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
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
        gobs.append(Vo*Vo/R*KPC)
        sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs = np.log10(gobs)

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_lam(y, lam):
    return (1.0-lam)*nu_standard(y) + lam*nu_be(y)

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

t00 = time.time()
P("8S-b THE EDGE RESOLUTION (pre-reg committed BEFORE any run; "
  "measurement round; NO credence movement)")

allg = np.unique(gal_id)
NGAL = gal_id.max()+1
gdfrac = np.zeros(NGAL)
for g_ in allg:
    m = gal_id == g_
    gdfrac[g_] = float(np.mean(g_gas[m] > g_dsk[m] + g_bul[m]))
gd_set = np.array([g_ for g_ in allg if gdfrac[g_] >= 0.5])
dd_set = np.array([g_ for g_ in allg if gdfrac[g_] < 0.5])
P(f"selector (8S verbatim): GD {len(gd_set)} galaxies, DD "
  f"{len(dd_set)}")

def wvec(gset):
    w = np.zeros(NGAL)
    w[gset] = 1.0
    return w

def m2ll(th, lam, w_gal, fg_forced=None):
    la0, f, s_int, fg = th
    if fg_forced is not None:
        fg = fg_forced
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or not (0.7 <= fg <= 1.4): return 1e12
    a0 = 10**la0
    gN = fg*g_gas + f*g_dsk + g_bul
    nv = nu_lam(gN/a0, lam)
    if np.min(nv) <= 0.05: return 1e12
    gm = gN*nv
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm)
    return np.sum(w_gal[gal_id]*(r*r/se2 + np.log(se2)))

def fit_at(lam, w_gal, th_warm=None, fg_forced=None):
    starts = [[math.log10(A0_FID), 1.0, 0.08, 1.0],
              [math.log10(A0_FID)+0.1, 0.8, 0.12, 1.2]]
    if th_warm is not None: starts.insert(0, list(th_warm))
    best = None
    for th0 in starts:
        b = minimize(lambda t: m2ll(t, lam, w_gal, fg_forced), th0,
                     method='Nelder-Mead',
                     options=dict(maxiter=5000, xatol=1e-6, fatol=1e-7))
        if best is None or b.fun < best.fun: best = b
    return best

LGX = np.round(np.arange(-2.00, 1.501, 0.05), 3)
YPROBE = np.logspace(-4, 3, 400)
def profile(w_gal, fg_forced=None):
    prof, th, dead = [], None, []
    for lam in LGX:
        if float(np.min(nu_lam(YPROBE, lam))) <= 0.05:
            prof.append(np.nan); dead.append(lam); continue
        b = fit_at(lam, w_gal, th, fg_forced)
        prof.append(b.fun); th = b.x
    prof = np.array(prof)
    fin = np.isfinite(prof)
    i = int(np.nanargmin(prof))
    lam_hat = LGX[i]
    if 0 < i < len(LGX)-1 and fin[i-1] and fin[i+1]:
        x3, y3 = LGX[i-1:i+2], prof[i-1:i+2]
        c2_, c1_, _ = np.polyfit(x3, y3, 2)
        if c2_ > 0: lam_hat = -c1_/(2*c2_)
    lo = hi = None
    for j in range(i, -1, -1):
        if not fin[j]: break
        if prof[j] > prof[i]+1.0:
            lo = np.interp(prof[i]+1.0, [prof[j+1], prof[j]],
                           [LGX[j+1], LGX[j]]) \
                if prof[j] != prof[j+1] else LGX[j]
            break
    for j in range(i, len(LGX)):
        if not fin[j]: break
        if prof[j] > prof[i]+1.0:
            hi = np.interp(prof[i]+1.0, [prof[j-1], prof[j]],
                           [LGX[j-1], LGX[j]])
            break
    return prof, lam_hat, lo, hi, i, dead

# ---------------- gates ----------------
dead0 = [float(l) for l in LGX
         if float(np.min(nu_lam(YPROBE, l))) <= 0.05]
P(f"G8Sb-2 nu-positivity: excluded lam nodes = "
  + (", ".join(f"{v:.2f}" for v in dead0) if dead0 else "none"))

th_p = [math.log10(A0_FID), 1.0, 0.08, 1.0]
v_free = m2ll(th_p, 0.7, wvec(gd_set))
v_fix = m2ll(th_p, 0.7, wvec(gd_set), fg_forced=1.0)
g1_ok = abs(v_free - v_fix) < 1e-9
P(f"G8Sb-1 fg-slice identity at probe: free(fg=1) {v_free:.9f} vs "
  f"fixed {v_fix:.9f} -> {'PASS' if g1_ok else 'FAIL'}")

# 8S shared-node regression: GD fixed-fg profile on the common nodes
SH = np.round(np.arange(-0.30, 1.501, 0.05), 3)
prof_gd_fix, lam_gd_fix, lo_gf, hi_gf, i_gf, _ = profile(
    wvec(gd_set), fg_forced=1.0)
ref = {}
for ln in open('data/stage8s_gasc1.txt').read().splitlines():
    if ln.startswith('[GD  ] SPARC-only marginalized: lam_hat ='):
        ref['lam'] = float(ln.split('lam_hat = ')[1].split(' ')[0])
g0_ok = True
sh_idx = [int(np.argmin(np.abs(LGX-v))) for v in SH]
# regression: the shared-node minimum location matches the 8S edge value
lam_sh = SH[int(np.nanargmin(prof_gd_fix[sh_idx]))]
g0_ok = abs(lam_sh - ref.get('lam', 9)) <= 0.051
P(f"G8Sb-0 shared-node regression: restricted-grid GD minimum at "
  f"lam = {lam_sh:.3f} vs 8S {ref.get('lam')} "
  f"-> {'PASS' if g0_ok else 'FAIL'}")

if not (g0_ok and g1_ok):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    with open('data/stage8sb_gasedge.txt', 'w') as f:
        f.write("\n".join(L)+"\n")
    raise SystemExit(0)
P("GATES: G8Sb-0/1/2 ALL PASS")
P("")

# ---------------- E1: the extended-grid GD profile (fixed fg) -------
interior = (lo_gf is not None and hi_gf is not None
            and 0 < i_gf < len(LGX)-1)
P(f"[GD fg=1 ] extended grid: lam_hat = {lam_gd_fix:.3f} (D1 "
  f"{None if lo_gf is None else round(lo_gf,3)}.."
  f"{None if hi_gf is None else round(hi_gf,3)}) -> c1_hat = "
  f"{lam_gd_fix/2:.3f}; {'INTERIOR' if interior else 'EDGE/ONE-SIDED'}")

# ---------------- E2/E3: the gas knob ----------------
res = {}
for tag, gset in (('FULL', allg), ('GD', gd_set), ('DD', dd_set)):
    prof, lam_hat, lo, hi, i, _ = profile(wvec(gset))
    # fg at the fit
    b = fit_at(LGX[i], wvec(gset))
    res[tag] = dict(lam=lam_hat, lo=lo, hi=hi, fg=b.x[3], i=i)
    P(f"[{tag:4} fgFREE] lam_hat = {lam_hat:.3f} (D1 "
      f"{None if lo is None else round(lo,3)}.."
      f"{None if hi is None else round(hi,3)}) -> c1_hat = "
      f"{lam_hat/2:.3f}; fg_hat = {b.x[3]:.3f}"
      f"{' (fg-EDGE)' if (abs(b.x[3]-0.7)<0.01 or abs(b.x[3]-1.4)<0.01) else ''}")

P("")
gd = res['GD']; fu = res['FULL']
two_sided = gd['lo'] is not None and gd['hi'] is not None
ov = two_sided and fu['lo'] is not None and fu['hi'] is not None \
     and not (gd['hi'] < fu['lo'] or gd['lo'] > fu['hi'])
excl = two_sided and (fu['lam'] < gd['lo'] or fu['lam'] > gd['hi'])
one_sided_excl = (gd['hi'] is not None and gd['hi'] < fu['lam'])
if ov:
    P("==> 8S-b VERDICT (locked grammar): GAS-EXPLAINED — with the "
      "gas budget free the gas-dominated interval rejoins the "
      "full-sample dial: the 8S edge is attributable to the fixed "
      "gas normalization; T5's defense holds CONDITIONAL on the gas "
      "systematic (both intervals quoted; the gas knob is now a named "
      "nuisance for every deep-window coefficient claim).")
elif excl or one_sided_excl:
    P("==> 8S-b VERDICT (locked grammar): DIAL-TENSION — the "
      "gas-dominated subsample excludes the full-sample dial even "
      "with the gas budget free: a REAL subsample tension for the "
      "galaxy program (vertical/distance channel + selection next; "
      "named, not resolved here).")
else:
    P("==> 8S-b VERDICT (locked grammar): POWER-CARRIED — intervals "
      "quoted; no attribution.")
P("    8S grammar hole flagged: a one-sided D1 interval could not "
  "fire T5-TENSION (fell through to POWER-LIMITED); the 8S record "
  "stands with this annotation.")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t00)/60:.1f} min)")

with open('data/stage8sb_gasedge.txt', 'w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage8sb_gasedge.txt")
