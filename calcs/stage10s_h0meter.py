"""
STAGE 10S "H0METER" -- instrument-side gate harness (SKY-BLIND).

Pre-registration: PREREG-H0METER-DRAFT.md (UNCOMMITTED at build time; the
author's word commits it, THEN the sky legs fire). This script is the
firewall's step 1: gates G1-G5 on archived/synthetic inputs only.

FIREWALL DISCIPLINE (what this script never prints):
  - no anchored-subsample a0 central (leg A) -- G5 prints the bootstrap
    SIGMA only;
  - no flow-subsample a0 central (leg B ingredient) -- G3 prints scaling
    RATIOS and the measured exponent only;
  - no H0 of any kind from real V_obs. G4's H0 recoveries are synthetic
    truths (62/70/78) on synthesized observations.

Gates (bars in the pre-reg draft, fixed before this run):
  G2  census: token + fixed-width parsers agree; 97/45/3/28/2; the two
      SPARC reassignments (NGC3992 f_D=5, UGC06446 f_D=1) asserted.
  G1a flat regression: the 4H (b) machinery (nu_screen, la0+p+f_d free,
      full sample) reproduces archived a0 = 1.041e-10 to 0.5%.
  G1b hier regression: the 5M vertical machinery (BE, dv+dml channels,
      + lensing) reproduces archived -2lnL = -12152.49 (+-1.0) and
      a0 = 1.044e-10 to 0.5%.
  G3  scaling injection: flow subsample rescaled by s in {0.85,1.0,1.15}
      (lgobs -> lgobs - log10 s; g_bar invariant -- the pin-3 laws);
      measured d(ln a0)/d(ln s) must land in [-2.5,-1.5] for the analytic
      shortcut to be licensed (else leg B uses the measured curve).
  G4  leg-B solver (AMENDED 2026-09-20 pre-commit, first firing -- the
      original single-noisy-mock 1% bar was mis-posed: the fixed-point
      geometry amplifies fractional a0 offsets by 1/(gamma-1) = x1.71 at
      the measured gamma = 1.585, so one noisy realization sits at 3-5%
      by construction; see the pre-reg amendment note):
      G4a wiring: NOISELESS mocks recover H0_true to 0.2% at 62/70/78.
      G4b calibration: 5 noise seeds at each outer truth (measured
      per-point + 0.08 dex intrinsic); |mean error| <= max(1.5%, 2 SE);
      replication SD reported as the leg-B single-realization floor.
  G5  leg-A power: galaxy bootstrap (200 reps, seed 202) of the anchored
      +UMa subsample with per-galaxy distance jitter (anchored: e_D/D;
      UMa: shared 0.9/18 + per-galaxy depth 2.3/18 per the desk check);
      sigma(H0_A) > 15 km/s/Mpc => POWER-LIMITED.
  G6  literature bar (logged here per the pre-reg): targeted search
      2026-09-20 (session): no published execution of a0 -> H0 as a
      measurement found (scout round NOT-FOUND 2026-09-20; second
      targeted round run at harness time -- result noted in the gates
      file by the session log).
  G7  fires at sky time (GD-composition split of leg A).

UMa desk check (pin 4, executed 2026-09-20 pre-fit): SPARC UMa distance
18 +- 0.9 Mpc = Sorce et al. 2013b (ApJ 765, 94) mid-IR TFR with a
Cepheid/TRGB-calibrated zero point => ANCHOR-BASED => UMa joins LEG A
with a shared-distance nuisance. Leg A raw = 78, leg B raw = 97.

Writes data/stage10s_h0meter_gates.txt. No sky numbers in the file.
"""
import glob, math, os, sys
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
LN10 = math.log(10)
C_LIGHT = 299792458.0            # m/s
KMS_MPC = 3.240779e-20           # 1 km/s/Mpc in s^-1  (pin 5)

def h0_of_a0(a0):                # km/s/Mpc  (pin 5; machine-to-machine)
    return 2.0*math.pi*a0/C_LIGHT/KMS_MPC

def a0_of_h0(h0):
    return C_LIGHT*h0*KMS_MPC/(2.0*math.pi)

L = ["STAGE 10S H0METER -- instrument-side gates (sky-blind)", ""]

# ---------------- G2: census, two parsers ----------------
MRT = 'data/sparc/SPARC_Lelli2016c.mrt'
with open(MRT) as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1

import re
# parser 1: token split; parser 2: anchored regex (the .mrt's true layout
# deviates from its byte-by-byte header spec -- the 2026-09-20 census
# finding, reproduced at harness time; fixed-width slicing is invalid)
RX = re.compile(r'^\s*(\S+)\s+(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d)\s')
tok_fd, rex_fd = {}, {}
meta = {}   # name -> (inc, Q, D, e_D, e_inc, f_D)
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        name = t[0]
        fdv_tok = int(t[4])
        tok_fd[name] = fdv_tok
        meta[name] = (float(t[5]), int(t[17]), float(t[2]), float(t[3]),
                      float(t[6]), fdv_tok)
    except ValueError:
        continue
    m = RX.match(l)
    if m and m.group(1) == name:
        rex_fd[name] = int(m.group(5))

agree = (set(tok_fd) == set(rex_fd) and
         all(tok_fd[n] == rex_fd[n] for n in tok_fd))
cnt = {m: sum(1 for v in tok_fd.values() if v == m) for m in (1, 2, 3, 4, 5)}
ok_counts = (cnt[1], cnt[2], cnt[3], cnt[4], cnt[5]) == (97, 45, 3, 28, 2)
ok_reassign = (tok_fd.get('NGC3992') == 5 and tok_fd.get('UGC06446') == 1)
g2 = agree and ok_counts and ok_reassign
L.append(f"G2 census: parsers {'AGREE' if agree else 'DISAGREE'}; "
         f"counts flow/TRGB/Cep/UMa/SN = {cnt[1]}/{cnt[2]}/{cnt[3]}/"
         f"{cnt[4]}/{cnt[5]} (bar 97/45/3/28/2); reassignments "
         f"NGC3992->SN, UGC06446->flow {'OK' if ok_reassign else 'FAIL'} "
         f"-> {'PASS' if g2 else 'FAIL'}")

# ---------------- point tables (the archived cuts) ----------------
g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
sigv_g_map, fd_g_map, dist_g_map = {}, {}, {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc, fdv = meta.get(name, (0, 3, 10.0, 1.0, 3.0, 1))
    if inc < 30 or q > 2: continue
    kept += 1
    irad = math.radians(inc)
    sv = math.hypot((eD/max(D, 1e-3))/LN10,
                    2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10)
    sigv_g_map[gi] = max(sv, 0.01)
    fd_g_map[gi] = fdv
    dist_g_map[gi] = (D, eD)
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
        sig.append(2*eV/Vo/LN10)
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs = np.log10(gobs)
ug = np.unique(gal_id)
NGal = len(ug)
FD_G = np.array([fd_g_map[g] for g in ug])
n_flow = int((FD_G == 1).sum())
n_anch = int(np.isin(FD_G, (2, 3, 5)).sum())
n_uma = int((FD_G == 4).sum())
L.append(f"quality-cut sample: {kept} galaxies, {len(gobs)} points; "
         f"leg membership after cuts: flow {n_flow} / anchored {n_anch} / "
         f"UMa {n_uma}  (leg A = {n_anch + n_uma})")
L.append("")

# ---------------- G1a: the 4H flat machinery, verbatim math ----------------
def nu_screen(y, p):
    return (1 - np.exp(-np.clip(y, 1e-12, None)**p))**(-1/(2*p))

def fit_flat_screen(idx, x0=(-9.92, 0.5, 1.0)):
    lg = lgobs[idx]
    def lo(v):
        la0, p, fd = v
        if not (0.05 <= p <= 3) or not (0.3 <= fd <= 3): return 1e9
        gb = (g_gas + fd*g_dsk + g_bul)[idx]
        r = lg - np.log10(gb*nu_screen(gb/10**la0, p))
        return float(np.mean(r*r))
    b = minimize(lo, list(x0), method='Nelder-Mead',
                 options={'xatol': 1e-4, 'fatol': 1e-10, 'maxiter': 3000})
    return b.x, b.fun

allidx = np.arange(len(gobs))
(la0_f, p_f, fd_f), rms_f = fit_flat_screen(allidx)
a0_flat = 10**la0_f
A0_4H = 1.041e-10
d_flat = abs(a0_flat - A0_4H)/A0_4H
g1a = d_flat < 0.005
L.append(f"G1a flat regression (4H-b machinery): a0 = {a0_flat:.4e} vs "
         f"archived 1.041e-10 (d = {100*d_flat:.2f}%), p = {p_f:.3f} "
         f"(arch 0.572), f_d = {fd_f:.2f} (arch 1.23) -> "
         f"{'PASS' if g1a else 'FAIL'}")

# ---------------- BE-form flat fit (the leg engine; G3/G4/G5) ----------
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))

def fit_flat_be(idx, dlg=None, x0=(-9.92, 1.2)):
    """Flat BE-form fit, (la0, fd) free; dlg = per-point lgobs shift."""
    lg = lgobs[idx] if dlg is None else lgobs[idx] + dlg
    gg, gd, gb_ = g_gas[idx], g_dsk[idx], g_bul[idx]
    def lo(v):
        la0, fd = v
        if not (-10.6 < la0 < -9.4) or not (0.3 <= fd <= 3): return 1e9
        gN = gg + fd*gd + gb_
        r = lg - np.log10(gN*nu_be(gN/10**la0))
        return float(np.mean(r*r))
    b = minimize(lo, list(x0), method='Nelder-Mead',
                 options={'xatol': 1e-5, 'fatol': 1e-12, 'maxiter': 3000})
    return b.x[0], b.x[1], b.fun

def fit_be_synth(idx, lg_synth, x0=(-9.92, 1.2)):
    gg, gd, gb_ = g_gas[idx], g_dsk[idx], g_bul[idx]
    def lo(v):
        la0, fd = v
        if not (-10.8 < la0 < -9.2) or not (0.3 <= fd <= 3): return 1e9
        gN = gg + fd*gd + gb_
        r = lg_synth - np.log10(gN*nu_be(gN/10**la0))
        return float(np.mean(r*r))
    b = minimize(lo, list(x0), method='Nelder-Mead',
                 options={'xatol': 1e-5, 'fatol': 1e-12, 'maxiter': 3000})
    return b.x[0], b.x[1], b.fun

# ---------------- G3: scaling injection (flow subsample) ----------------
flow_g = set(ug[FD_G == 1])
fidx = np.where(np.isin(gal_id, list(flow_g)))[0]
SVALS = (0.85, 1.0, 1.15)
la0_s = []
for s in SVALS:
    la0_i, _, _ = fit_flat_be(fidx, dlg=-math.log10(s))
    la0_s.append(la0_i)
la0_s = np.array(la0_s)
lns = np.log(np.array(SVALS))
lna = la0_s*LN10
expo = float(np.polyfit(lns, lna, 1)[0])
rat_lo = 10**(la0_s[0]-la0_s[1])
rat_hi = 10**(la0_s[2]-la0_s[1])
g3 = 1.5 <= abs(expo) <= 2.5
L.append(f"G3 scaling injection (flow, {len(flow_g)} gal, "
         f"{len(fidx)} pts): a0(s)/a0(1) = {rat_lo:.4f} @ s=0.85, "
         f"{rat_hi:.4f} @ s=1.15; measured d(ln a0)/d(ln s) = {expo:+.3f} "
         f"(analytic -2; bar |e| in [1.5,2.5]) -> "
         f"{'PASS (analytic shortcut licensed)' if g3 else 'FAIL (leg B uses the measured curve)'}"
         )
# absolutes deliberately not printed (firewall)

# ------- G4a/G4b: leg-B fixed-point solver (synthetic; amended) -------
S_INT_TRUE = 0.08
FD_TRUE = 1.2
gN_true = g_gas[fidx] + FD_TRUE*g_dsk[fidx] + g_bul[fidx]

def solve_h0(lg_cat):
    def F(h0p):
        sp = 73.0/h0p
        la0_i, _, _ = fit_be_synth(fidx, lg_cat - math.log10(sp))
        return 10**la0_i - a0_of_h0(h0p)
    lo, hi = 55.0, 90.0
    flo = F(lo)
    if flo*F(hi) > 0: return None
    for _ in range(16):
        mid = 0.5*(lo+hi)
        fm = F(mid)
        if flo*fm <= 0: hi = mid
        else: lo = mid; flo = fm
    return 0.5*(lo+hi)

g4a_rows, g4a_ok = [], True
for h0_true in (62.0, 70.0, 78.0):
    s_true = 73.0/h0_true
    a0_true = a0_of_h0(h0_true)
    lg_cat = np.log10(s_true*gN_true*nu_be(gN_true/a0_true))
    h0_rec = solve_h0(lg_cat)
    if h0_rec is None:
        g4a_rows.append(f"  {h0_true:.0f}: NO BRACKET"); g4a_ok = False
        continue
    derr = abs(h0_rec - h0_true)/h0_true
    ok = derr < 0.002
    g4a_ok &= ok
    g4a_rows.append(f"  H0_true {h0_true:.0f} -> {h0_rec:.3f} "
                    f"({100*derr:.3f}%) {'OK' if ok else 'FAIL'}")
L.append("G4a leg-B solver wiring (NOISELESS mocks, f_d_true=1.2; "
         "bar 0.2%):")
L.extend(g4a_rows)

g4b_rows, g4b_ok = [], True
sd_floor = []
for h0_true in (62.0, 78.0):
    s_true = 73.0/h0_true
    a0_true = a0_of_h0(h0_true)
    g_true = gN_true*nu_be(gN_true/a0_true)
    errs = []
    for seed in (42, 101, 202, 303, 404):
        rng = np.random.default_rng(seed)
        lg_cat = (np.log10(s_true*g_true)
                  + rng.normal(0, np.sqrt(sig2[fidx] + S_INT_TRUE**2)))
        h0_rec = solve_h0(lg_cat)
        errs.append(100*(h0_rec - h0_true)/h0_true)
    e = np.array(errs)
    se = e.std(ddof=1)/math.sqrt(len(e))
    bar = max(1.5, 2*se)
    ok = abs(e.mean()) <= bar
    g4b_ok &= ok
    sd_floor.append(e.std(ddof=1))
    g4b_rows.append(f"  H0_true {h0_true:.0f}: mean {e.mean():+.2f}% "
                    f"(bar +-{bar:.2f}%), SD {e.std(ddof=1):.2f}% "
                    f"{'OK' if ok else 'FAIL'}")
L.append("G4b leg-B calibration (5 noise seeds x outer truths; "
         "unbiasedness bar |mean| <= max(1.5%, 2 SE)):")
L.extend(g4b_rows)
L.append(f"  leg-B single-realization noise floor (replication SD): "
         f"~{np.mean(sd_floor):.1f}% of H0 (~"
         f"{0.01*np.mean(sd_floor)*70:.1f} km/s/Mpc at H0~70); the "
         f"fixed-point amplification 1/(gamma-1) = "
         f"{1.0/(abs(expo)-1.0):.2f} at measured gamma = {abs(expo):.3f}")
g4_ok = g4a_ok and g4b_ok
L.append(f"G4 -> {'PASS' if g4_ok else 'FAIL'} "
         f"(G4a {'P' if g4a_ok else 'F'} / G4b {'P' if g4b_ok else 'F'})")

# ---------------- G5: leg-A power (MASKED central) ----------------
legA_g = ug[np.isin(FD_G, (2, 3, 4, 5))]
aidx = np.where(np.isin(gal_id, legA_g))[0]
# per-galaxy relative distance sigma
rel_sig = {}
uma_flag = {}
for g in legA_g:
    D, eD = dist_g_map[g]
    if fd_g_map[g] == 4:
        rel_sig[g] = 2.3/18.0          # independent depth part
        uma_flag[g] = True
    else:
        rel_sig[g] = eD/max(D, 1e-3)
        uma_flag[g] = False
UMA_SHARED = 0.9/18.0

rng5 = np.random.default_rng(202)
la0_reps = []
n_reps = 200
# per-galaxy point index cache
gpts = {g: np.where(gal_id == g)[0] for g in legA_g}
for _ in range(n_reps):
    pick = rng5.choice(legA_g, size=len(legA_g), replace=True)
    sc_shared = 1.0 + rng5.normal(0, UMA_SHARED)
    rows, dshift = [], []
    for g in pick:
        s_g = 1.0 + rng5.normal(0, rel_sig[g])
        s_g = min(max(s_g, 0.5), 1.5)
        if uma_flag[g]:
            s_g *= sc_shared
        rows.append(gpts[g])
        dshift.append(np.full(len(gpts[g]), -math.log10(s_g)))
    ridx = np.concatenate(rows)
    dlg = np.concatenate(dshift)
    la0_i, _, _ = fit_flat_be(ridx, dlg=dlg)
    la0_reps.append(la0_i)
la0_reps = np.array(la0_reps)
a0_reps = 10**la0_reps
sig_a0 = float(np.std(a0_reps))
sig_h0 = h0_of_a0(sig_a0)          # sigma converts linearly
g5_power = sig_h0 <= 15.0
L.append("")
L.append(f"G5 leg-A power ({len(legA_g)} galaxies = {n_anch} anchored + "
         f"{n_uma} UMa, {len(aidx)} pts; 200 boot reps, galaxy resample + "
         f"distance jitter incl. UMa shared nuisance; CENTRAL MASKED):")
L.append(f"  sigma(a0_A) = {sig_a0:.3e} m/s^2  ->  sigma(H0_A) = "
         f"{sig_h0:.2f} km/s/Mpc  (bar: >15 = POWER-LIMITED)  -> "
         f"{'PASS (leg A is powered)' if g5_power else 'POWER-LIMITED'}")
L.append("")

# ---------------- G1b: the 5M hier machinery (slow; last) ----------------
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25
S_ML = 0.1*LN10
gmap = {g: i for i, g in enumerate(ug)}
gidx = np.array([gmap[g] for g in gal_id])
SIGV = np.array([sigv_g_map[g] for g in ug])
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]
ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

def m2hv(th, nu, dml, dv, sv, w_g):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    gm_ = gN*nu(gN/a0)
    se2 = sig2 + s_int*s_int
    r = lgobs - np.log10(gm_) - dv[gidx]
    out = np.sum(w_g[gidx]*(r*r/se2 + np.log(se2)))
    lg = l_gbar[lmask] + dlt
    rl = l_gobs[lmask] - (lg + np.log10(nu(10**lg/a0)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
    out += np.sum(w_g*dv*dv/(sv*sv))
    return out

def fit_v(nu, sv, tol=0.05, max_rounds=15, th0=None):
    w_g = np.ones(NGal)
    dml = np.zeros(NGal)
    dv = np.zeros(NGal)
    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08, 0.0]] if rd == 0 else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2hv(t, nu, dml, dv, sv, w_g), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int, dlt = best.x
        se2c = s_int*s_int
        for _ in range(3):
            fac = f*np.exp(dml[gidx])
            gN = g_gas + fac*g_dsk + g_bul
            r0 = lgobs - np.log10(gN*nu(gN/10**la0))
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                w = 1.0/(sig2[mm] + se2c)
                dv[gi2] = np.sum(w*r0[mm])/(np.sum(w) + 1.0/sv[gi2]**2)
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    rr = (lgobs[mm] - np.log10(gN2*nu(gN2/10**la0))
                          - dv[gi2])
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2hv(best.x, nu, dml, dv, sv, w_g)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2hv(t, nu, dml, dv, sv, w_g), list(best.x),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml, dv

print("G1b: running the 5M hier regression (slow)...", flush=True)
bhier, _, _ = fit_v(nu_be, SIGV)
a0_hier = 10**bhier.x[0]
A0_5M = 1.044e-10
FUN_5M = -12152.49
d_hier = abs(a0_hier - A0_5M)/A0_5M
d_fun = abs(bhier.fun - FUN_5M)
g1b = (d_hier < 0.005) and (d_fun < 1.0)
L.append(f"G1b hier regression (5M machinery, BE dv-ON): -2lnL = "
         f"{bhier.fun:.2f} (arch -12152.49, d = {d_fun:.2f}), a0 = "
         f"{a0_hier:.4e} (arch 1.044e-10, d = {100*d_hier:.2f}%), "
         f"f_ML = {bhier.x[1]:.2f} (arch 1.16) -> "
         f"{'PASS' if g1b else 'FAIL'}")
L.append("")

L.append("G6 literature bar (second targeted round, 2026-09-20, session "
         "logs): no published a0->H0 lock-inversion measurement found. "
         "Nearest neighbor PRIMARY-READ at abstract grade: Schombert, "
         "McGaugh & Lelli 2020 (AJ 160, 71; arXiv:2006.08615) measure "
         "H0 = 75.1 +- 2.3 (stat) +- 1.5 (sys) with the bTFR as a "
         "conventional distance-ladder rung on the SAME catalog and the "
         "SAME 50-anchor/95-flow split (CosmicFlows-3 velocities) -- no "
         "lock anywhere in the method; MUST-CITE + census "
         "cross-validation (their '50 galaxies with accurate distances' "
         "= our anchor count). Landscape: the 2pi a0 ~ cH0 coincidence "
         "literature (e.g. arXiv:1110.2580) has no measurement "
         "execution; Zenodo 20377329 'Quantum Horizon Gravity I' = "
         "derivation-side deposit, not a measurement. -> PASS")
L.append("")
allpass = g2 and g1a and g1b and g3 and g4_ok
L.append(f"GATE SUMMARY: G2 {'P' if g2 else 'F'} | G1a "
         f"{'P' if g1a else 'F'} | G1b {'P' if g1b else 'F'} | G3 "
         f"{'P' if g3 else 'F'} | G4a/b {'P' if g4_ok else 'F'} | G5 "
         f"{'powered' if g5_power else 'POWER-LIMITED'} | G6 P")
L.append(f"INSTRUMENT {'READY' if allpass else 'NOT READY'} -- sky legs "
         f"remain LOCKED until the pre-reg is committed (author's word).")

out = "\n".join(L)
print(out)
with open('data/stage10s_h0meter_gates.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage10s_h0meter_gates.txt")
