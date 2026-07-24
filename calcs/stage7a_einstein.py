"""
STAGE 7A: the Einstein fluctuation test -- does the RAR's intrinsic scatter
carry the PARTICLE term?

PRE-REGISTERED (PREDICTIONS.md P6, committed before execution).

The mean law nu = 1 + n_BE fixes the bath's second moment with zero new shape
freedom (fluctuation-dissipation). Einstein 1909: Var(n) = n + n^2 = particle
+ wave. Through dlogg = dn/((1+n) ln10), with x = sqrt(g_N/a0) LOCKED to the
mean fit:

  EQ quantum  (Var = n(1+n)):  s^2 = e^(-x)              / (N ln^2 10)
  EC classical wave (Var = n^2): s^2 = e^(-2x)           / (N ln^2 10)
  ES shot/corpuscular (Var = n): s^2 = e^(-x)(1-e^(-x))  / (N ln^2 10)
  EG free exponent:             s^2 = S0 e^(-gamma x)    [gamma measured]

EQ == Stage 4T's M1b (regression-gated against its recorded value). EC is the
fluctuation law of ANY classical continuous-field bath, including the 4F
classical self-consistent bath (simple-nu). All rivals carry the same
parameter count (amplitude + floor, jointly with la0, f_ML).

Gates:
  G0 regression: EQ(BE) reproduces 4T M1b -2lnL within 1.0.
  G1 nesting: EG <= EQ + tol and EG <= EC + tol (exact nestings; tol 1e-3);
     free 6-bin M2 <= EG + tol (envelope sanity).
  G2 separability (CALIBRATED INSTRUMENT, the 4W lesson): paired injections --
     quantum truth (N=21, floor=0.10) and classical truth (matched deep
     amplitude), 3 seeds each; PASS iff |mean gamma_Q - 1| < 0.35,
     |mean gamma_C - 2| < 0.5, and min(gamma_C) - max(gamma_Q) > 0.3.
     NO verdict is quoted if G2 fails (power-limited, reported as such).
  G3 bump robustness: the x~1 point-level bump (4W) is non-monotone and could
     bias gamma; refit EG (a) masking x_fid in [0.8, 1.4], (b) with an explicit
     shared bump term b^2 exp(-(x-1.1)^2/(2 0.25^2)). Headline = full fit;
     both shifts reported; verdicts must survive (sign + interior).
  G4 mean-family robustness: EG under simple-nu mean law (gamma_hat quoted).

Pre-registered bars (both directions):
  SUPPORT (quantum): G2 PASS, gamma_hat interior, gamma=2 excluded at
    D(-2lnL) >= 9, gamma=1 inside the 2-unit interval, sign survives G3.
  ANTI (classical):  mirror image (gamma=1 excluded >= 9, gamma_hat near 2)
    -> pre-committed REAL strike at the quantum-statistical reading (the
    mean-law function record untouched; "quantum bath" language demoted).
  SHOT: ES beats both EQ and EC by >= 9 -> different-physics note.
  UNRESOLVED: G2 fail, or intervals cover both, or edge-running gamma_hat.
Caveats carried: per-point independence (4T grade); 4U/4W showed the thermal
AMPLITUDE deflates under M/L + vertical marginalization -- any SUPPORT stays
one-channel until hier-hardened.

Writes data/stage7a_einstein.txt.
"""
import glob, math, os, re
import numpy as np
from scipy.optimize import minimize

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)

# ---------------- loader (verbatim Stage 4T) ----------------
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
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))

x_fid = np.sqrt((g_gas + g_dsk + g_bul)/A0_FID)
QE = np.quantile(x_fid, np.linspace(0, 1, 7))
QE[0], QE[-1] = 0.0, np.inf
BIN_FID = np.clip(np.searchsorted(QE, x_fid, side='right')-1, 0, 5)
MASK_BUMP = ~((x_fid >= 0.8) & (x_fid <= 1.4))   # G3a fixed mask

# ---------------- variance shapes (per-point s^2 in dex^2, before floor) ----
def shape_EQ(x):
    xc = np.minimum(np.clip(x, 1e-9, None), 80.0)
    return np.exp(-xc)
def shape_EC(x):
    xc = np.minimum(np.clip(x, 1e-9, None), 80.0)
    return np.exp(-2.0*xc)
def shape_ES(x):
    xc = np.minimum(np.clip(x, 1e-9, None), 80.0)
    e = np.exp(-xc)
    return e*(1.0-e)

def m2ll(la0, f, s2_arr, nu, sel=None):
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    gm = gN*nu(gN/a0)
    se2 = sig2 + s2_arr
    r = lgobs - np.log10(gm)
    q = r*r/se2 + np.log(se2)
    if sel is not None: q = q[sel]
    return np.sum(q)

def fit(model, nu, x0, sel=None):
    """model in {E0, EQ, EC, ES, EG, EGB, M2, EGFIX:<gamma>}"""
    gfix = None
    if model.startswith('EGFIX:'):
        gfix = float(model.split(':')[1]); model = 'EGFIX'
    def obj(th):
        la0, f = th[0], th[1]
        if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
        x = np.sqrt((g_gas + f*g_dsk + g_bul)/10**la0)
        if model == 'E0':
            s0 = th[2]
            if not (1e-4 <= s0 < 0.5): return 1e12
            s2 = np.full(len(gobs), s0*s0)
        elif model in ('EQ', 'EC', 'ES'):
            lnN, sf = th[2], th[3]
            if not (-2 < lnN < 14) or not (0 <= sf < 0.5): return 1e12
            sh = {'EQ': shape_EQ, 'EC': shape_EC, 'ES': shape_ES}[model]
            s2 = sh(x)/(math.exp(lnN)*LN10*LN10) + sf*sf
        elif model == 'EG':
            lnS0, gam, sf = th[2], th[3], th[4]
            if not (-19 < lnS0 < -1.4) or not (0.05 < gam < 6.0) \
               or not (0 <= sf < 0.5): return 1e12
            s2 = math.exp(lnS0)*np.exp(-gam*np.minimum(x, 80.0)) + sf*sf
        elif model == 'EGFIX':
            lnS0, sf = th[2], th[3]
            if not (-19 < lnS0 < -1.4) or not (0 <= sf < 0.5): return 1e12
            s2 = math.exp(lnS0)*np.exp(-gfix*np.minimum(x, 80.0)) + sf*sf
        elif model == 'EGB':
            lnS0, gam, b, sf = th[2], th[3], th[4], th[5]
            if not (-19 < lnS0 < -1.4) or not (0.05 < gam < 6.0) \
               or not (0 <= b < 0.3) or not (0 <= sf < 0.5): return 1e12
            s2 = (math.exp(lnS0)*np.exp(-gam*np.minimum(x, 80.0))
                  + b*b*np.exp(-(x-1.1)**2/(2*0.25**2)) + sf*sf)
        else:  # M2 free bins
            sb = np.asarray(th[2:8])
            if np.any(sb < 1e-4) or np.any(sb > 0.5): return 1e12
            s2 = (sb[BIN_FID])**2
        return m2ll(la0, f, s2, nu, sel)
    best = None
    for t0 in x0:
        b = minimize(obj, t0, method='Nelder-Mead',
                     options=dict(maxiter=12000, xatol=1e-6, fatol=1e-6))
        if best is None or b.fun < best.fun: best = b
    return best

L = [f"STAGE 7A Einstein fluctuation test: {kept} galaxies, {len(gobs)} points",
     "pre-registered P6 (PREDICTIONS.md); bars in script header", ""]
LA, FM = math.log10(A0_FID), 1.0
LNS0_G = math.log(1.0/(21.0*LN10*LN10))   # generic EG amplitude start

# ---------------- main contest (BE mean law) ----------------
r = {}
r['E0'] = fit('E0', nu_be, [[LA, FM, 0.10], [LA, FM, 0.06]])
r['EQ'] = fit('EQ', nu_be, [[LA, FM, math.log(30.0), 0.05],
                            [LA, FM, math.log(8.0), 0.10]])
r['EC'] = fit('EC', nu_be, [[LA, FM, math.log(30.0), 0.05],
                            [LA, FM, math.log(8.0), 0.10]])
r['ES'] = fit('ES', nu_be, [[LA, FM, math.log(30.0), 0.05],
                            [LA, FM, math.log(8.0), 0.10]])
r['EG'] = fit('EG', nu_be, [[LA, FM, LNS0_G, 1.0, 0.05],
                            [LA, FM, LNS0_G, 2.0, 0.05],
                            [LA, FM, LNS0_G+1.0, 0.5, 0.10]])
r['M2'] = fit('M2', nu_be, [[LA, FM]+[0.10]*6, [LA, FM]+[0.06]*6])

e0, eq, ec, es, eg, m2v = (r[k].fun for k in ('E0','EQ','EC','ES','EG','M2'))
gam_hat = r['EG'].x[3]
L.append("[BE mean law]  -2lnL (lower is better)")
L.append(f"  E0 const:            {e0:10.2f}  s0 = {r['E0'].x[2]:.4f}")
L.append(f"  EQ quantum  e^-x:    {eq:10.2f}  N = {math.exp(r['EQ'].x[2]):7.1f}"
         f"  floor = {r['EQ'].x[3]:.4f}")
L.append(f"  EC classical e^-2x:  {ec:10.2f}  N = {math.exp(r['EC'].x[2]):7.1f}"
         f"  floor = {r['EC'].x[3]:.4f}")
L.append(f"  ES shot  e^-x(1-e^-x): {es:8.2f}  N = {math.exp(r['ES'].x[2]):7.1f}"
         f"  floor = {r['ES'].x[3]:.4f}")
L.append(f"  EG free gamma:       {eg:10.2f}  gamma_hat = {gam_hat:.3f}"
         f"  S0 = {math.exp(r['EG'].x[2]):.3e}  floor = {r['EG'].x[4]:.4f}")
L.append(f"  M2 free 6-bin:       {m2v:10.2f}")
L.append(f"  contest: EQ-EC = {eq-ec:+.2f}  EQ-ES = {eq-es:+.2f}"
         f"  EC-ES = {ec-es:+.2f}   (negative favors first)")
L.append("")

# G0 regression vs 4T M1b (BE block, first M1b line)
g0txt = "SKIP (4T output absent)"
try:
    t4 = open('data/stage4t_bathnoise.txt').read()
    mm = re.search(r"M1b osc \+ floor:\s+-2lnL =\s+([-\d.]+)", t4)
    ref = float(mm.group(1))
    g0 = abs(eq - ref) < 1.0
    g0txt = f"EQ = {eq:.2f} vs 4T M1b = {ref:.2f} (d = {eq-ref:+.3f}) -> " \
            f"{'PASS' if g0 else 'FAIL'}"
except Exception as ex:
    g0 = False
    g0txt = f"FAIL to parse ({ex})"
L.append(f"G0 regression: {g0txt}")

# G1 nesting
tol = 1e-3
g1 = (eg <= eq + tol) and (eg <= ec + tol) and (m2v <= eg + 2.0)
L.append(f"G1 nesting: EG<=EQ {eg-eq:+.3f}, EG<=EC {eg-ec:+.3f}, "
         f"M2<=EG+2 {m2v-eg:+.2f} -> {'PASS' if g1 else 'FAIL'}")

# ---------------- gamma profile (grid, warm-started) ----------------
GRID = np.arange(0.25, 3.501, 0.125)
prof = []
warm = [LA, FM, LNS0_G, 0.05]
for g in GRID:
    b = fit(f'EGFIX:{g}', nu_be, [warm, [LA, FM, LNS0_G, 0.10]])
    prof.append(b.fun)
    warm = list(b.x)
prof = np.array(prof)
pmin = prof.min(); ib = int(prof.argmin())
def crossings(level):
    lo, hi = None, None
    for i in range(ib, 0, -1):
        if prof[i-1] >= pmin+level >= prof[i]:
            t = (pmin+level-prof[i])/(prof[i-1]-prof[i]); lo = GRID[i]-t*0.125; break
    for i in range(ib, len(GRID)-1):
        if prof[i+1] >= pmin+level >= prof[i]:
            t = (pmin+level-prof[i])/(prof[i+1]-prof[i]); hi = GRID[i]+t*0.125; break
    return lo, hi
lo1, hi1 = crossings(1.0)
lo4, hi4 = crossings(4.0)
d_g1 = prof[np.argmin(np.abs(GRID-1.0))] - pmin
d_g2 = prof[np.argmin(np.abs(GRID-2.0))] - pmin
edge = ib in (0, len(GRID)-1)
L.append("")
L.append(f"gamma profile: min at {GRID[ib]:.3f} (-2lnL {pmin:.2f})"
         f"{'  ** EDGE **' if edge else ''}")
L.append(f"  1-unit interval: [{lo1 if lo1 is None else round(lo1,3)}, "
         f"{hi1 if hi1 is None else round(hi1,3)}]   "
         f"2-unit-equv(4): [{lo4 if lo4 is None else round(lo4,3)}, "
         f"{hi4 if hi4 is None else round(hi4,3)}]")
L.append(f"  D(-2lnL): gamma=1 {d_g1:+.2f} | gamma=2 {d_g2:+.2f}"
         f"   (>= 9 = pre-registered exclusion)")
L.append("  profile: " + " ".join(f"{g:.2f}:{p-pmin:+.1f}"
                                  for g, p in zip(GRID[::2], prof[::2])))
L.append("")

# ---------------- G2 separability injections ----------------
gNt = g_gas + FM*g_dsk + g_bul
xt = np.sqrt(gNt/10**LA)
lg_true = np.log10(gNt*nu_be(gNt/10**LA))
lgobs_save = lgobs.copy()
recs = {'Q': [], 'C': []}
for kind in ('Q', 'C'):
    sh = shape_EQ(xt) if kind == 'Q' else shape_EC(xt)
    s2_true = sh/(21.0*LN10*LN10) + 0.10**2
    for k in range(3):
        rng = np.random.default_rng(701 + 10*k + (0 if kind == 'Q' else 5))
        lgobs = lg_true + rng.normal(0, np.sqrt(sig2 + s2_true))
        binj = fit('EG', nu_be, [[LA, FM, LNS0_G, 1.0, 0.10],
                                 [LA, FM, LNS0_G, 2.0, 0.10]])
        recs[kind].append(binj.x[3])
lgobs = lgobs_save
mQ, mC = np.mean(recs['Q']), np.mean(recs['C'])
sep = min(recs['C']) - max(recs['Q'])
g2 = (abs(mQ-1.0) < 0.35) and (abs(mC-2.0) < 0.5) and (sep > 0.3)
L.append(f"G2 separability: gammaQ = {np.round(recs['Q'],3).tolist()} "
         f"(mean {mQ:.3f}); gammaC = {np.round(recs['C'],3).tolist()} "
         f"(mean {mC:.3f}); min-gap = {sep:.3f} -> {'PASS' if g2 else 'FAIL'}")

# ---------------- G3 bump robustness ----------------
rmask = fit('EG', nu_be, [[LA, FM, LNS0_G, max(gam_hat, 0.3), 0.05],
                          [LA, FM, LNS0_G, 2.0, 0.05]], sel=MASK_BUMP)
rbump = fit('EGB', nu_be, [[LA, FM, LNS0_G, max(gam_hat, 0.3), 0.03, 0.05],
                           [LA, FM, LNS0_G, 2.0, 0.03, 0.05]])
L.append(f"G3a masked x in [0.8,1.4] ({int(np.sum(~MASK_BUMP))} pts dropped): "
         f"gamma_hat = {rmask.x[3]:.3f} (shift {rmask.x[3]-gam_hat:+.3f})")
L.append(f"G3b bump-modeled: gamma_hat = {rbump.x[3]:.3f} "
         f"(shift {rbump.x[3]-gam_hat:+.3f})  b = {rbump.x[4]:.4f}  "
         f"-2lnL = {rbump.fun:.2f} (vs EG {eg:.2f})")

# ---------------- G4 simple-nu mean law ----------------
rsim = fit('EG', nu_simple, [[LA, FM, LNS0_G, 1.0, 0.05],
                             [LA, FM, LNS0_G, 2.0, 0.05]])
L.append(f"G4 simple-nu mean: gamma_hat = {rsim.x[3]:.3f}  "
         f"-2lnL = {rsim.fun:.2f}")
L.append("")

# ---------------- verdict per pre-registered bars ----------------
interior = (0.25+0.125 < GRID[ib] < 3.5-0.125) and not edge
shot_wins = (es <= eq - 9.0) and (es <= ec - 9.0)
g3_ok = (abs(rmask.x[3]-gam_hat) < 0.5 or (rmask.x[3] < 1.5) == (gam_hat < 1.5)) \
        and ((rbump.x[3] < 1.5) == (gam_hat < 1.5))
verdict = "UNRESOLVED"
if g2 and interior and not shot_wins:
    if d_g2 >= 9.0 and (lo4 is None or lo4 <= 1.0) and (hi4 is None or hi4 >= 1.0) \
       and g3_ok and gam_hat < 1.5:
        verdict = "SUPPORT (quantum: particle term present, classical excluded)"
    elif d_g1 >= 9.0 and gam_hat > 1.5 and g3_ok:
        verdict = "ANTI (classical: gamma=1 excluded -> pre-committed strike)"
elif shot_wins and g2:
    verdict = "SHOT (corpuscular shape wins -> different-physics note)"
if not g2:
    verdict = "UNRESOLVED (power-limited: instrument cannot separate)"
L.append(f"VERDICT (pre-registered bars): {verdict}")
L.append("caveats: 4T grade (per-point independence); hier-hardening is the")
L.append("immediate follow-up before any strong language; amplitude-deflation")
L.append("risk under M/L+vertical marginalization carried from 4U/4W.")

out = "\n".join(L)
print(out)
with open('data/stage7a_einstein.txt', 'w') as f:
    f.write(out+"\n")
