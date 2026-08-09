"""
STAGE 10D/10E ADDENDUM -- ROUND-25 verification + condition execution
(the standing verify-reviewer-math rule: every load-bearing number the
ROUND-25 referee introduced is independently re-computed here before
adoption; his conditions C1-C6 are executed in the printout).

His new numbers under verification:
  GA-1 the ceiling correction: the 10E lemma at perfect geometry
       (eta2 = e_s = e_a = 1) gives (1/4)sqrt(om*Om) = HALF of
       g_close = (1/2)sqrt(om*Om); GE-4b was a definitional identity
       (0.5*sq vs 0.5*sqrt(x_loc x_amb)/2pi -- same formula), not a
       vertex property.  [his §3(ii); C4a]
  GA-2 the lam/g_close table: lam_max/g_close = 0.304/0.140/0.029/
       0.876/0.505/0.324 at the six anchors; e_a needed to reach
       g_close centrally = 1.020/1.722/2.689 (galaxy) and 2.313
       (binary x_loc = 0.5).  [his §4; C4b/C6]
  GA-3 galaxy g_close 0.0092/0.0164/0.0299 < 0.072 = the 10A band
       floor (the 10C table's own "below the 10A band" rows) -- the
       band cannot be the universal closure requirement.  [his §3(iv)]
  GA-4 the split-null characterization: his 10 fresh kappa=1 draws
       gave D(F1-F4) mean 2.00, max 4.08, 0/10 >= 11.8. Verified
       here with 6 INDEPENDENT fresh draws (seed 61) -- the claim
       under test is "the null D never approaches the 11.8 bar".
  GA-5 the A2 seed-replication: his seed-42 row [0.970, 0.999,
       1.001, 0.970, 0.815, 1.026] mean 0.964 -- re-run draw-level
       with the stage machinery.
  GA-6 the off-regime bias: truth a0 = 3.9e-11 gave his 6-draw mean
       1.067 (mild +0.07 up-bias) -- spot-verified with 3 draws
       (seed 77; band [0.90, 1.25]).
  GA-7 the no-lens refits: F1 kappa = 1.496 (a0 3.903e-11); F4
       1.306/1.024/0.020, D = +22.32 -- re-fit with lensing removed.

Conditions executed in the printout: C1 (running direction
UNRESOLVED, decomposition reading retired), C2 (ridge caveat), C3
(bias annotation, deconvolved ~1.43), C4 (ceiling + band clauses
corrected), C5 (O5-GEOMETRY scope + pre-stated kill), C6 (composed
reading).

Writes data/stage10de_addendum.txt.
"""
import glob, math, os, time
import numpy as np
import sympy as sp
from scipy.optimize import minimize, minimize_scalar

T00 = time.time()
KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25
LN10 = math.log(10)
S_ML = 0.1*LN10
XB1, XB2 = 0.7, 2.5

OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

say("=" * 78)
say("STAGE 10D/10E ADDENDUM -- ROUND-25 VERIFICATION + CONDITIONS C1-C6")
say("=" * 78)

# ---------------- GA-1: the ceiling correction (sympy) ----------------
say("")
say("GA-1 the ceiling correction (his §3(ii)):")
eta2_s, es_s, ea_s, om_s, Om_s = sp.symbols(
    'eta2 e_s e_a omega Omega', positive=True)
lam_sym = eta2_s*es_s*ea_s*sp.sqrt(om_s*Om_s)/4
lam_perfect = lam_sym.subs([(eta2_s, 1), (es_s, 1), (ea_s, 1)])
g_close = sp.sqrt(om_s*Om_s)/2
ratio = sp.simplify(lam_perfect/g_close)
say(f"  lemma at eta2 = e_s = e_a = 1: {lam_perfect}")
say(f"  g_close = {g_close};  ratio = {ratio}")
ok1 = ratio == sp.Rational(1, 2)
say(f"  -> {'CONFIRMED' if ok1 else 'MISMATCH'}: the perfect-geometry "
    "vertex is (1/4)sqrt(om*Om) = HALF of g_close; the stage's GE-4b")
say("     'ceiling identity' compared 0.5*sq against 0.5*sq -- a")
say("     DEFINITIONAL identity, not a lemma property (his catch;")
say("     cuts AGAINST the mechanism: the vertex is x2 weaker than")
say("     the T2b clause advertised).")

# ---------------- GA-2 + GA-3: lam/g_close and the band ----------------
say("")
say("GA-2 lam/g_close ratios + e_a-needed (his §4 table):")
ETA_C, EA_C = 0.765, 0.086
GAMMA_LO = 0.072
def nbe(x): return 1.0/math.expm1(x)
ANCH = [
    ('binary  x_loc=0.5', 0.5, 1.0954, 0.502, 1.0),
    ('binary  x_loc=1.0', 1.0, 1.0954, 0.502, 0.25),
    ('binary  x_loc=2.0', 2.0, 1.0954, 0.502, 0.0625),
    ('galaxy  x_loc=0.0953', 0.0953, 0.1411, 6.63, 0.3),
    ('galaxy  x_loc=0.3',    0.3,    0.1411, 6.63, 0.3),
    ('galaxy  x_loc=1.0',    1.0,    0.1411, 6.63, 0.3),
]
HIS_RATIO = [0.304, 0.140, 0.029, 0.876, 0.505, 0.324]
HIS_EA = {0: 2.313, 3: 1.020, 4: 1.722, 5: 2.689}
ok2 = True
ok3 = True
for i, (nm, xl, xa, na, es_c) in enumerate(ANCH):
    sq = math.sqrt(xl*xa)/(2*math.pi)
    gcl = 0.5*sq
    ns = nbe(xl)
    enh = math.sqrt((ns + 1.0)*na)
    lam_c = (ETA_C/4.0)*es_c*EA_C*sq*enh
    es_max = min(1.0, 3.0*es_c) if i < 3 else 0.5
    enh_max = math.sqrt(max((ns + 1.0)*na, ns*(na + 1.0)))
    lam_max = (2.0/4.0)*1.0*es_max*0.20*sq*enh_max
    r = lam_max/gcl
    okr = abs(r - HIS_RATIO[i]) < 0.01
    ok2 = ok2 and okr
    line = (f"  {nm:<22} lam_max/g_close = {r:.3f} (his "
            f"{HIS_RATIO[i]:.3f}) {'OK' if okr else 'MISMATCH'}")
    if i in HIS_EA:
        ea_need = EA_C*gcl/lam_c
        oke = abs(ea_need - HIS_EA[i]) < 0.02
        ok2 = ok2 and oke
        line += (f";  e_a needed = {ea_need:.3f} (his {HIS_EA[i]:.3f}) "
                 f"{'OK' if oke else 'MISMATCH'}")
    say(line)
    if i >= 3:
        ok3 = ok3 and (gcl < GAMMA_LO)
say(f"  GA-2 -> {'CONFIRMED' if ok2 else 'MISMATCH'} (with e_a ~ 1 the "
    "deep-galaxy anchor reaches g_close CENTRALLY; max-favorable 88%)")
say("")
say("GA-3 the band mis-specification (his §3(iv)):")
say(f"  galaxy g_close = 0.0092/0.0164/0.0299 H, ALL below the 10A band "
    f"floor {GAMMA_LO} -> {'CONFIRMED' if ok3 else 'MISMATCH'}")
say("  (the 10C G5 table's own rows say 'below the 10A band'; a")
say("   requirement the galaxies clear dispersively at sub-band g_close")
say("   cannot be the bar the galaxy exchange vertex must clear --")
say("   his structural point, adopted as C4b.)")

# ---------------- 10D machinery (verbatim from the stage) ----------------
meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta[t[0]] = (float(t[5]), int(t[17]))
    except ValueError:
        continue
g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q = meta.get(name, (0, 3))
    if inc < 30 or q > 2: continue
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
LGOBS = np.log10(gobs)
ug = np.unique(gal_id)
NGal = len(ug)
gmap = {g: i for i, g in enumerate(ug)}
gidx = np.array([gmap[g] for g in gal_id])
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]
ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID
gN_fid = g_gas + 1.0*g_dsk + g_bul
x_fid = np.sqrt(gN_fid/A0_FID)
BSP = np.where(x_fid < XB1, 0, np.where(x_fid < XB2, 1, 2))
xl_ = np.sqrt(10**l_gbar[lmask]/A0_FID)
BLE = np.where(xl_ < XB1, 0, np.where(xl_ < XB2, 1, 2))

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def occ_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 0.0, 1.0/np.expm1(np.minimum(x, 40)))
def nu_kap(y, kap):
    return 1.0 + kap*occ_be(y)

def m2core(kap_sp, kap_le, th, dml, w_g, lgv, use_lens=True):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm = gN*(1.0 + kap_sp*occ_be(gN/a0))
    se2 = sig2 + s_int*s_int
    r = lgv - np.log10(gm)
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    if use_lens:
        lg = l_gbar[lmask] + dlt
        yl = 10**lg/a0
        rl = l_gobs[lmask] - (lg + np.log10(1.0 + kap_le*occ_be(yl)))
        out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
        out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    return out

def kaps_from(mode, tk, bsp, ble):
    if mode == 'one':
        return tk[0], tk[0]
    kd, kt, kl = tk
    ksp = np.where(bsp == 0, kd, np.where(bsp == 1, kt, kl))
    kle = np.where(ble == 0, kd, np.where(ble == 1, kt, kl))
    return ksp, kle

NK = {'one': 1, 'three': 3}

def m2k(t, mode, dml, w_g, lgv, use_lens=True):
    nk = NK[mode]
    tk, th = t[:nk], t[nk:]
    if not all(0.02 < k < 3.0 for k in tk): return 1e12
    ksp, kle = kaps_from(mode, tk, BSP, BLE)
    return m2core(ksp, kle, th, dml, w_g, lgv, use_lens)

def fit_plain(mode, starts, rounds=2, max_rounds=8, lgv=None, tol=0.05,
              use_lens=True):
    if lgv is None: lgv = LGOBS
    w_g = np.ones(NGal)
    nk = NK[mode]
    dml = np.zeros(NGal)
    best, prev = None, None
    for rd in range(max_rounds):
        ss = ([list(best.x)] if best is not None else []) + starts
        bb = None
        for t0 in ss:
            b = minimize(lambda t: m2k(t, mode, dml, w_g, lgv, use_lens),
                         t0, method='Nelder-Mead',
                         options=dict(maxiter=4000+1500*nk, xatol=1e-6,
                                      fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        tk, (la0, f, s_int, dlt) = best.x[:nk], best.x[nk:]
        ksp, _ = kaps_from(mode, tk, BSP, BLE)
        se2c = s_int*s_int
        for gi2 in range(NGal):
            mm = GIDXS[gi2]
            kk = ksp if np.isscalar(ksp) else ksp[mm]
            def od(dl):
                fc = f*math.exp(dl)
                gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                rr = lgv[mm] - np.log10(gN2*(1.0 + kk*occ_be(gN2/10**la0)))
                s2 = sig2[mm] + se2c
                return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
            dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                       method='bounded').x
        cur = m2k(best.x, mode, dml, w_g, lgv, use_lens)
        if rd >= rounds - 1 and prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2k(t, mode, dml, w_g, lgv, use_lens),
                 list(best.x), method='Nelder-Mead',
                 options=dict(maxiter=4000+1500*nk, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best

S1 = [[1.0, math.log10(A0_FID), 1.0, 0.08, 0.0],
      [0.7, math.log10(A0_FID), 1.0, 0.08, 0.0],
      [1.0, math.log10(8e-11), 1.4, 0.07, 0.0],
      [0.7, math.log10(8e-11), 1.4, 0.07, 0.0]]
S4B = [1.0, 1.5, 1.0, math.log10(A0_FID), 1.0, 0.08, 0.0]

def inject(nu_truth, a0_truth, rr):
    d_t = rr.normal(0, S_ML, NGal)
    gN_t = g_gas + 1.0*np.exp(d_t[gidx])*g_dsk + g_bul
    return (np.log10(gN_t*nu_truth(gN_t/a0_truth))
            + rr.normal(0, np.sqrt(sig2 + 0.08**2)))

# ---------------- GA-4: split-null extension ----------------
say("")
say("GA-4 split-null extension (6 fresh kappa=1 draws, seed 61; claim")
say("under test: the null D never approaches the 11.8 bar):")
rng4 = np.random.default_rng(61)
Dn = []
for dr in range(6):
    lg_i = inject(nu_be, A0_FID, rng4)
    b1n = fit_plain('one', S1[:2], lgv=lg_i)
    b4n = fit_plain('three', [[b1n.x[0]]*3 + list(b1n.x[1:]), S4B],
                    lgv=lg_i)
    Dn.append(b1n.fun - b4n.fun)
    say(f"    draw {dr+1}: D = {Dn[-1]:+.2f}   "
        f"({(time.time()-T00)/60:.1f} min)")
ok4 = max(Dn) < 11.8
say(f"  6-draw null: mean {np.mean(Dn):.2f}, max {max(Dn):+.2f}, "
    f">= 11.8: {sum(d >= 11.8 for d in Dn)}/6 -> "
    f"{'CONFIRMED' if ok4 else 'MISMATCH'} (his 10-draw mean 2.00 / "
    "max 4.08 / 0/10; stage's 3 draws max +8.53; observed data +21.74)")

# ---------------- GA-5: seed-42 draw-level replication ----------------
say("")
say("GA-5 A2 seed-42 replication (his row: 0.970 0.999 1.001 0.970 "
    "0.815 1.026, mean 0.964):")
rng5 = np.random.default_rng(42)
rec42 = []
for dr in range(6):
    lg_i = inject(nu_be, A0_FID, rng5)
    bi = fit_plain('one', S1, lgv=lg_i)
    rec42.append(bi.x[0])
HIS42 = [0.970, 0.999, 1.001, 0.970, 0.815, 1.026]
dmax = max(abs(a - b) for a, b in zip(rec42, HIS42))
m42 = float(np.mean(rec42))
ok5 = (dmax <= 0.05) and (abs(m42 - 0.964) <= 0.02)
say("  re-run: " + " ".join(f"{r:.3f}" for r in rec42)
    + f"  mean {m42:.3f}")
say(f"  max per-draw |diff| = {dmax:.3f} (bar 0.05); mean diff = "
    f"{abs(m42-0.964):.3f} (bar 0.02) -> "
    f"{'CONFIRMED' if ok5 else 'MISMATCH'}")

# ---------------- GA-6: off-regime spot ----------------
say("")
say("GA-6 off-regime bias spot (truth a0 = 3.9e-11, 3 draws, seed 77;")
say("his 6-draw mean 1.067; spot band [0.90, 1.25]):")
rng6 = np.random.default_rng(77)
rec_off = []
for dr in range(3):
    lg_i = inject(nu_be, 3.9e-11, rng6)
    bi = fit_plain('one', [[1.0, math.log10(3.9e-11), 1.0, 0.08, 0.0],
                           [1.0, math.log10(8e-11), 1.4, 0.07, 0.0]],
                   lgv=lg_i)
    rec_off.append(bi.x[0])
m_off = float(np.mean(rec_off))
ok6 = 0.90 <= m_off <= 1.25
say("  recoveries: " + " ".join(f"{r:.3f}" for r in rec_off)
    + f"  mean {m_off:.3f} -> "
    f"{'CONSISTENT with his +0.07-grade bias' if ok6 else 'OUT OF BAND'}")

# ---------------- GA-7: no-lens refits ----------------
say("")
say("GA-7 no-lens refits (his: F1 1.496 / a0 3.903e-11; F4 "
    "1.306/1.024/0.020, D +22.32):")
b1nl = fit_plain('one', S1, max_rounds=10, use_lens=False)
b4nl = fit_plain('three', [[b1nl.x[0]]*3 + list(b1nl.x[1:]), S4B],
                 max_rounds=10, use_lens=False)
Dnl = b1nl.fun - b4nl.fun
ok7 = (abs(b1nl.x[0] - 1.496) <= 0.02
       and abs(b4nl.x[0] - 1.306) <= 0.02
       and abs(b4nl.x[1] - 1.024) <= 0.02
       and abs(Dnl - 22.32) <= 1.5)
say(f"  F1 kappa = {b1nl.x[0]:.3f}, a0 = {10**b1nl.x[1]:.3e}; F4 = "
    f"{b4nl.x[0]:.3f}/{b4nl.x[1]:.3f}/{b4nl.x[2]:.3f}, D = {Dnl:+.2f}")
say(f"  -> {'CONFIRMED' if ok7 else 'MISMATCH'} (the deep arm is "
    "SPARC-driven; lensing does not carry kappa_d)")

# ---------------- conditions ----------------
say("")
say("=" * 78)
say("ROUND-25 CONDITIONS EXECUTED (letters corrected in place):")
say("=" * 78)
say("(C1) 10D DECOMPOSITION RETIRED: the treatment-stable content of")
say("     K-SPLIT is 'one kappa is rejected' (plain D = +21.7, vertical")
say("     D = +14.4). The RUNNING DIRECTION is UNRESOLVED: plain kappa_d")
say("     = 1.317 > kappa_t = 1.036 REVERSES under the vertical channel")
say("     (F4v 1.479 < 1.614). 'kappa_d/kappa_t' is NOT quotable as a")
say("     decomposition; 'kappa_mid >> kappa_tail' is near-vacuous (the")
say("     tail is the guard-pinned proxy) and is DROPPED.")
say("(C2) RIDGE CAVEAT (print wherever K-SPLIT is cited): the split is")
say("     measured on the low-a0/high-f_ML ridge (a0 = 3.87e-11 = 0.37x")
say("     horizon, f_ML = 1.52). It is NOT evidence against the kappa=1")
say("     closure: the temperature-locked F2 world sits at kappa = 0.925")
say("     (Planck) / 0.888 (SH0ES) at +73 lnL. K-SPLIT re-scopes the")
say("     kappa = 1 target to the deep limit; it does not exclude it.")
say("(C3) A2 BIAS ANNOTATION: the kappa=1 estimator carries a mild")
say("     +0.05..+0.07 upward bias at the operative low-a0 regime")
say("     (his 6-seed table, means 0.964-1.067); observed kappa = 1.503")
say("     reads deconvolved ~1.43. No letter change.")
say("(C4) 10E CLAUSES CORRECTED: (a) the perfect-geometry vertex is")
say("     (1/4)sqrt(om*Om) = HALF of g_close (GA-1 exact); GE-4b was the")
say("     definitional identity g_close = (1/2)sqrt(om*Om), NOT a vertex")
say("     property -- the T2b 'supplies exactly the closure coupling'")
say("     clause is RETRACTED (the vertex is x2 weaker than advertised).")
say("     (b) Galaxy anchors are scored against g_close, not the 10A")
say("     band (GA-3): lam_max/g_close reaches 0.88 at the deep-galaxy")
say("     anchor; e_a ~ 1 closes it centrally (GA-2).")
say("(C5) O5-GEOMETRY SCOPED (the successor, strike-bearing at the")
say("     NEXT round): Axis 1 = derive the ambient cloud's FLUCTUATING")
say("     l=2/l=0 amplitude ratio (the exchange-active object; the")
say("     static EFE |q| = 0.086 is a different object), score against")
say("     g_close. Axis 2 = derive the exchange leg's requirement (6X")
say("     dynamical saturation ~H vs 6U KMS state-statistic, which")
say("     needs only coupling-existence + equilibration). PRE-STATED")
say("     KILL (reviewer's words): fluctuating e_a ~ 0.086 AND the")
say("     requirement band-like => the R22-cond-4 strike fires (15->8).")
say("(C6) COMPOSED READING (the honest form): 10E V-GRAY = 'input-")
say("     limited on an un-derived amplitude and a mis-specified")
say("     requirement', NOT 'structurally short'. No upgrade without")
say("     the O5-GEOMETRY computation.")
say("")
say("CREDENCE (mechanical, pre-signed): 10D K-SPLIT -> HOLD 15; 10E")
say("V-GRAY + ruling (B) ESCAPE-AXIS LIVE -> HOLD 15, NO strike.")
say("bath-mechanism conditional = 15 (six O5 rounds, zero strikes,")
say("zero rises). anomaly-real = 53 UNTOUCHED.")
say("=" * 78)
ok_all = ok1 and ok2 and ok3 and ok4 and ok5 and ok6 and ok7
say(f"GATES: GA-1:{'P' if ok1 else 'F'} GA-2:{'P' if ok2 else 'F'} "
    f"GA-3:{'P' if ok3 else 'F'} GA-4:{'P' if ok4 else 'F'} "
    f"GA-5:{'P' if ok5 else 'F'} GA-6:{'P' if ok6 else 'F'} "
    f"GA-7:{'P' if ok7 else 'F'}  -> "
    f"{'ALL REVIEWER NUMBERS CONFIRMED' if ok_all else 'MISMATCH -- do not adopt without diagnosis'}")
say("")
say(f"done ({(time.time()-T00)/60:.1f} min)")

with open('data/stage10de_addendum.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage10de_addendum.txt")
