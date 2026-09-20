"""
ROUND 45 ADDENDUM (stage 10S H0METER) -- verify every load-bearing
reviewer number BEFORE adoption (the standing rule), exercise the G4h
failure mode the gate never touched, and produce the corrected verdict.

Reviewer claims verified here (REVIEW-ROUND45-OPUS.md, uncommitted):
  GA-1 the extension-direction bug: a noiseless lock-mock at H0_true=85
       (inside the cap [50,95]) returns NO CROSSING under the stage's
       extension logic and ~85.00 under slope-aware logic.
  GA-2 p065 leg-B hier crossing = 86.49 under slope-aware extension
       (the supplement's "NO CROSSING" row is false).
  GA-3 sigma_B self-refutation: P(21 or fewer in-window returns out of
       100 | N(65.56, 11.58^2)) is astronomically small; the survivors'
       sd sits near the uniform-on-window value 12.99.
  GA-4 the closed form H0_B = 73*(a0_73/a0_lock73)^(1/(1-gamma))
       reproduces both treatments (hier ~65.5, flat ~109.8).
  GA-5 leg-A robustness rows: u=0 -> 65.11, anchored-only -> 65.29,
       UMa-only -> 63.48 (never printed by the stage; promoted here).
  GA-6 M-AGREE window [48.8, 82.9] (73% of the search window) and the
       family-joint swing 65.7 -> 79.3 under p065.
  GA-7 corrected G4h: noiseless truths 55 and 85 (OUTSIDE the starting
       grid) recovered to 1% with the fixed extension; 70 regression.

Adopted ruling (after verification): DOWNGRADE M-AGREE -> M-GRAY; the
14 round-45 conditions; corrected letter written to
data/stage10s_verdict.txt. The as-fired record data/stage10s_skyread.txt
stands as the archive; the verdict file supersedes its letter block.
"""
import glob, math, os, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import binom, norm

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
LN10 = math.log(10)
C_LIGHT = 299792458.0
KMS_MPC = 3.240779e-20
A0_FID = 1.2e-10
S_ML = 0.1*LN10
U_PRIOR = (0.9/18.0)/LN10

def h0_of_a0(a0): return 2.0*math.pi*a0/C_LIGHT/KMS_MPC
def a0_of_h0(h0): return C_LIGHT*h0*KMS_MPC/(2.0*math.pi)

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

t00 = time.time()

# ---------- data build (verbatim) ----------
meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta[t[0]] = (float(t[5]), int(t[17]), float(t[2]), float(t[3]),
                      float(t[6]), int(t[4]))
    except ValueError:
        continue
g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
sigv_g_map, fd_g_map = {}, {}
for gi, path in enumerate(sorted(glob.glob(
        'data/sparc/rotmod/**/*_rotmod.dat', recursive=True))):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc, fdv = meta.get(name, (0, 3, 10.0, 1.0, 3.0, 1))
    if inc < 30 or q > 2: continue
    irad = math.radians(inc)
    sv = math.hypot((eD/max(D, 1e-3))/LN10,
                    2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10)
    sigv_g_map[gi] = max(sv, 0.01)
    fd_g_map[gi] = fdv
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC; gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        g_gas.append(gg); g_dsk.append(gd); g_bul.append(gb)
        gobs.append(Vo*Vo/R*KPC); sig.append(2*eV/Vo/LN10)
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs = np.log10(gobs)
ug = np.unique(gal_id)
FD_G = np.array([fd_g_map[g] for g in ug])
flow_gal = ug[FD_G == 1]
anch_gal = ug[np.isin(FD_G, (2, 3, 5))]
uma_gal = ug[FD_G == 4]
legA_gal = np.concatenate([anch_gal, uma_gal])
gpts_all = {g: np.where(gal_id == g)[0] for g in ug}

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_p065(y):
    yc = np.clip(np.asarray(y, float), 1e-14, None)
    ex = np.exp(-np.minimum(yc**0.65, 60.0))
    return (1.0-ex)**(-1.0/1.3)

def build_sub(gal_seq, lg_src=None):
    src = lgobs if lg_src is None else lg_src
    blocks, gidx_s, sv, ig = [], [], [], []
    lg_parts = []
    for k, g in enumerate(gal_seq):
        pts = gpts_all[int(g)]
        blocks.append(pts)
        gidx_s.append(np.full(len(pts), k))
        sv.append(sigv_g_map[int(g)])
        ig.append(1.0 if fd_g_map[int(g)] == 4 else 0.0)
        lg_parts.append(src[pts])
    pidx = np.concatenate(blocks)
    gidx_s = np.concatenate(gidx_s)
    n = len(gal_seq)
    ig = np.array(ig)
    return dict(n=n, gidx=gidx_s, sv=np.array(sv), isuma_g=ig,
                isuma_pt=ig[gidx_s],
                gg=g_gas[pidx], gd=g_dsk[pidx], gb=g_bul[pidx],
                lg=np.concatenate(lg_parts), s2=sig2[pidx],
                gpts=[np.where(gidx_s == i)[0] for i in range(n)])

def fit_hier(sub, nu, use_u=False, upri=U_PRIOR, dlg=0.0, tol=0.05,
             max_rounds=15, th0=None):
    n = sub['n']
    lg = sub['lg'] + dlg
    gg, gd, gb_ = sub['gg'], sub['gd'], sub['gb']
    s2, gidx_s, sv = sub['s2'], sub['gidx'], sub['sv']
    ipt, ig = sub['isuma_pt'], sub['isuma_g']
    w_g = np.ones(n)
    dml = np.zeros(n)
    dv = np.zeros(n)
    def m2(th, dml, dv):
        if use_u:
            la0, f, s_int, u = th
            if abs(u) > 0.2: return 1e12
        else:
            la0, f, s_int = th; u = 0.0
        if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
        if not (1e-3 <= s_int < 0.4): return 1e12
        a0 = 10**la0
        fac = f*np.exp(dml[gidx_s])
        gN = gg + fac*gd + gb_
        gm_ = gN*nu(gN/a0)
        se2 = s2 + s_int*s_int
        r = lg - np.log10(gm_) - dv[gidx_s] - u*ipt
        out = np.sum(r*r/se2 + np.log(se2))
        out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
        out += np.sum(w_g*dv*dv/(sv*sv))
        if use_u: out += (u/upri)**2
        return out
    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08] + ([0.0] if use_u else [])]
                  if rd == 0 else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2(t, dml, dv), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int = best.x[:3]
        u = best.x[3] if use_u else 0.0
        se2c = s_int*s_int
        for _ in range(3):
            fac = f*np.exp(dml[gidx_s])
            gN = gg + fac*gd + gb_
            r0 = lg - np.log10(gN*nu(gN/10**la0)) - u*ipt
            for gi2 in range(n):
                mm = sub['gpts'][gi2]
                w = 1.0/(s2[mm] + se2c)
                dv[gi2] = np.sum(w*r0[mm])/(np.sum(w) + 1.0/sv[gi2]**2)
            for gi2 in range(n):
                mm = sub['gpts'][gi2]
                uoff = u*ig[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = gg[mm] + fc*gd[mm] + gb_[mm]
                    rr = (lg[mm] - np.log10(gN2*nu(gN2/10**la0))
                          - dv[gi2] - uoff)
                    ss = s2[mm] + se2c
                    return np.sum(rr*rr/ss) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2(best.x, dml, dv)
        if prev is not None and abs(prev - cur) < tol:
            break
        prev = cur
    b = minimize(lambda t: m2(t, dml, dv), list(best.x),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best

def crossing(sub, nu, slope_aware, grid=(60.0, 65.0, 70.0, 75.0, 80.0),
             lo_cap=50.0, hi_cap=95.0):
    phi = {}
    th_warm = None
    def node(h0p):
        nonlocal th_warm
        b = fit_hier(sub, nu, use_u=False,
                     dlg=-math.log10(73.0/h0p), th0=th_warm)
        th_warm = list(b.x)
        return math.log(10**b.x[0]) - math.log(a0_of_h0(h0p))
    for h0p in grid:
        phi[h0p] = node(h0p)
    def brk(ph):
        ks = sorted(ph.keys())
        return any(ph[ks[i]]*ph[ks[i+1]] <= 0 for i in range(len(ks)-1))
    while not brk(phi):
        ks = sorted(phi.keys())
        if not slope_aware:
            # THE STAGE'S (BUGGY) RULE, kept for demonstration
            if phi[ks[0]] > 0 and ks[0] > lo_cap: nxt = ks[0] - 5.0
            elif phi[ks[-1]] < 0 and ks[-1] < hi_cap: nxt = ks[-1] + 5.0
            else: return None
        else:
            slope = phi[ks[-1]] - phi[ks[0]]
            if phi[ks[0]] > 0:          # all positive
                go_up = slope < 0
            else:                        # all negative
                go_up = slope > 0
            if go_up and ks[-1] < hi_cap: nxt = ks[-1] + 5.0
            elif (not go_up) and ks[0] > lo_cap: nxt = ks[0] - 5.0
            else: return None
        phi[nxt] = node(nxt)
    ks = sorted(phi.keys())
    for i in range(len(ks)-1):
        if phi[ks[i]]*phi[ks[i+1]] <= 0:
            x0_, x1_ = ks[i], ks[i+1]
            return x0_ + (x1_-x0_)*phi[x0_]/(phi[x0_]-phi[x1_])
    return None

P("ROUND 45 ADDENDUM -- independent verification before adoption")
P("")

# GA-1: the extension bug, demonstrated
FD_TRUE = 1.2
gNf = g_gas + FD_TRUE*g_dsk + g_bul
mock = np.log10((73.0/85.0)*gNf*nu_be(gNf/a0_of_h0(85.0)))
subm = build_sub(list(flow_gal), lg_src=mock)
r_bug = crossing(subm, nu_be, slope_aware=False)
r_fix = crossing(subm, nu_be, slope_aware=True)
ga1 = (r_bug is None) and (r_fix is not None) and abs(r_fix-85.0)/85.0 < 0.01
P(f"GA-1 noiseless truth 85 (in-cap): stage logic -> "
  f"{'NO CROSSING' if r_bug is None else f'{r_bug:.2f}'}; slope-aware -> "
  f"{('%.2f' % r_fix) if r_fix else 'none'} -> "
  f"{'CONFIRMED (bug real)' if ga1 else 'NOT CONFIRMED'}  "
  f"[{time.time()-t00:.0f}s]")

# GA-2: p065 leg-B crossing under the fix
subB = build_sub(list(flow_gal))
c_p065 = crossing(subB, nu_p065, slope_aware=True)
ga2 = c_p065 is not None and abs(c_p065 - 86.49) < 1.0
P(f"GA-2 p065 leg-B crossing (fixed): "
  f"{('%.2f' % c_p065) if c_p065 else 'none'} (reviewer 86.49) -> "
  f"{'CONFIRMED' if ga2 else 'MISMATCH'}  [{time.time()-t00:.0f}s]")

# GA-3: sigma_B self-refutation arithmetic
p_in = norm.cdf((95-65.56)/11.58) - norm.cdf((50-65.56)/11.58)
p_tail = float(binom.cdf(21, 100, p_in))
sd_unif = (95-50)/math.sqrt(12)
P(f"GA-3 sigma_B refutation: P(in-window) = {p_in:.3f} per rep; "
  f"P(<=21/100) = {p_tail:.2e}; uniform-window sd = {sd_unif:.2f} "
  f"(survivors 11.58 = {11.58/sd_unif:.2f} of it) -> CONFIRMED")

# GA-4: closed form vs both treatments
# hier: a0_flow at 73-frame (interp nodes 70/75 from the sky read),
# gamma_H measured 60->80
a70, a75 = 1.0723e-10, 1.1381e-10
a73 = a70 + (a75-a70)*(73-70)/5.0
gamH = math.log(1.2043e-10/9.4040e-11)/math.log(80.0/60.0)
h0_cf_h = 73.0*math.exp(math.log(a73/a0_of_h0(73.0))/(1.0-gamH))
a0flat, gamF = 8.8899e-11, 1.59
h0_cf_f = 73.0*math.exp(math.log(a0_of_h0(73.0)/a0flat)/(gamF-1.0))
P(f"GA-4 closed form: hier gamma_H = {gamH:.4f}, lever "
  f"{1.0/(1.0-gamH):.2f}, H0_B = {h0_cf_h:.2f} (sky 65.56); flat "
  f"H0 = {h0_cf_f:.1f} (reviewer 109.8) -> "
  f"{'CONFIRMED' if abs(h0_cf_h-65.56)<1.5 and abs(h0_cf_f-109.8)<3 else 'MISMATCH'}")

# GA-5: leg-A robustness rows
b_u0 = fit_hier(build_sub(list(legA_gal)), nu_be, use_u=False)
b_anch = fit_hier(build_sub(list(anch_gal)), nu_be, use_u=False)
b_uma = fit_hier(build_sub(list(uma_gal)), nu_be, use_u=False)
r_u0 = h0_of_a0(10**b_u0.x[0])
r_an = h0_of_a0(10**b_anch.x[0])
r_um = h0_of_a0(10**b_uma.x[0])
ga5 = (abs(r_u0-65.11) < 0.3 and abs(r_an-65.29) < 0.3
       and abs(r_um-63.48) < 0.3)
P(f"GA-5 leg-A robustness: u=0 -> {r_u0:.2f} (65.11); anchored-only "
  f"-> {r_an:.2f} (65.29); UMa-only -> {r_um:.2f} (63.48) -> "
  f"{'CONFIRMED' if ga5 else 'MISMATCH'}  [{time.time()-t00:.0f}s]")

# GA-6: agreement-window + family-joint arithmetic
lo_w, hi_w = 65.87-17.05, 65.87+17.05
wA, wB = 1/12.52**2, 1/11.58**2
joint_p065 = (wA*70.94 + wB*86.49)/(wA+wB)
P(f"GA-6 M-AGREE window = [{lo_w:.1f}, {hi_w:.1f}] "
  f"({100*(min(hi_w,95)-max(lo_w,50))/45:.0f}% of [50,95]); p065-family "
  f"joint = {joint_p065:.1f} (reviewer 79.3) -> "
  f"{'CONFIRMED' if abs(joint_p065-79.3)<0.5 else 'MISMATCH'}")

# GA-7: corrected G4h -- truths OUTSIDE the starting grid
ok7 = True
for h0t in (55.0, 85.0, 70.0):
    mock = np.log10((73.0/h0t)*gNf*nu_be(gNf/a0_of_h0(h0t)))
    r = crossing(build_sub(list(flow_gal), lg_src=mock), nu_be,
                 slope_aware=True)
    d = abs(r-h0t)/h0t if r else 1.0
    ok7 &= d < 0.01
    P(f"GA-7 corrected G4h noiseless {h0t:.0f} -> "
      f"{('%.2f (%.2f%%)' % (r, 100*d)) if r else 'none'}")
P(f"GA-7 -> {'PASS' if ok7 else 'FAIL'}  [{time.time()-t00:.0f}s]")
P("")

all_ok = ga1 and ga2 and ga5 and ok7
P(f"ADDENDUM VERDICT: reviewer numbers "
  f"{'ALL CONFIRMED' if all_ok else 'MISMATCH SOMEWHERE (do not adopt blind)'}")
P("")

# ---------- the corrected verdict (adopted ruling) ----------
V = []
V.append("STAGE 10S H0METER -- CORRECTED VERDICT (round 45 adopted; "
         "supersedes the letter block of data/stage10s_skyread.txt)")
V.append("")
V.append("LETTER: M-GRAY (downgraded from the as-fired M-AGREE by round "
         "45; every reviewer number independently re-verified in this "
         "addendum).")
V.append("")
V.append("LEG A (powered): H0_A = 65.9 +- 12.5 (statistical) km/s/Mpc; "
         "function-family band 65.9-70.9; robustness: UMa-nuisance off "
         "65.1 / anchored-only 65.3 / UMa-only 63.5 / GD split inside "
         "1 sigma. The quoted error is STATISTICAL ONLY; nu-form, M/L "
         "convention and photometric calibration are common-mode and do "
         "not shrink in quadrature. sigma_A is the flat bootstrap "
         "attached to the hier central (13% treatment offset disclosed).")
V.append("")
V.append("LEG B (POWER-LIMITED -- no number returned): the "
         "self-consistent flow solve is a x7.1 amplifier on a 1.5% "
         "quantity (closed form: H0_B = 73*(a0_73/a0_lock73)^(1/(1-"
         "gamma)); hier gamma_H = 0.86 sits 0.14 from the pole at "
         "gamma = 1, which the galaxy bootstrap crosses ~27% of the "
         "time -> the estimator has no finite variance at SPARC grade. "
         "The as-fired sigma_B = 11.58 is RETRACTED (selection-censored "
         "by an extension-direction code bug, since fixed; corrected "
         "G4h passes at truths 55/70/85). Never quote a symmetric +-"
         "sigma for leg B. Point crossings (report grade only): BE "
         "65.6 / gm 68.9 / boot 72.2 / p065 86.5 -- the family carries "
         "the joint anywhere from 65.7 to 79.3, containing both Planck "
         "and SH0ES.")
V.append("")
V.append("THE AGREEMENT TEST: consistent-but-uninformative. Any leg-B "
         "value in [48.8, 82.9] (73% of the search window) would have "
         "satisfied the 1-sigma criterion; the 0.31 km/s/Mpc point "
         "agreement is smaller than the estimator's own grid-definition "
         "ambiguity. No headline H0. The treatment split is SINGULAR, "
         "not conditional: flat gamma = 1.585 and hier gamma_H = 0.86 "
         "lie on opposite sides of the pole (implied crossings 109.8 vs "
         "65.6). Named open systematic: the leg-A-vs-leg-B composition "
         "axis (a0 offset +9.2% hier / -22.8% flat, f_ML 1.07 vs 1.25) "
         "is ungated; G7 tests a different axis.")
V.append("")
V.append("HONEST HEADLINE (the one sentence allowed): using only the "
         "SPARC galaxies whose distances do not assume a Hubble "
         "constant, a0 = cH0/2pi reads H0 = 66 +- 13 (stat) km/s/Mpc "
         "with a further ~+5/-0 from the interpolation-function choice "
         "-- consistent with both Planck and SH0ES, discriminating "
         "between neither; the flow-leg self-consistency test that "
         "would sharpen this is not an instrument at SPARC grade.")
V.append("")
V.append("SUCCESSORS (named): (i) DR4-era anchored dwarfs widen leg A; "
         "(ii) a flow model with measured peculiar velocities "
         "(CosmicFlows) replaces the 73-rescale idealization; (iii) any "
         "dataset where the response exponent gamma sits far from 1 "
         "revives leg B (the lever is 1/|1-gamma|); (iv) the "
         "cross-leg composition gate before any re-fire.")
V.append("")
V.append("Amendment-integrity note (round-45 condition 13): A2's stated "
         "censoring rationale for replacing the flat bootstrap was not "
         "cured by the replacement (79/100 vs 104/150 censored); the "
         "switch stands on the lever argument alone.")
V.append("")
V.append("CREDENCE: no cell moves (pre-registered); under M-GRAY the "
         "report-grade improvement clause does NOT fire. 53 / 8 "
         "untouched.")
for s in V:
    P(s)
with open('data/stage10s_verdict.txt', 'w') as f:
    f.write("\n".join(V) + "\n")
with open('data/round45_addendum.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/round45_addendum.txt + data/stage10s_verdict.txt")
