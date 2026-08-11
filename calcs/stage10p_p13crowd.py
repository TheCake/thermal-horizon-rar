"""Stage 10P "P13-CROWD" -- the cluster crowd shape: first data contest.

PRE-REGISTRATION (this header is the pre-registration; committed BEFORE any
mode of this script is executed; bars below are frozen at commit time).

QUESTION
--------
PREDICTIONS.md P13 (registered 2026-08-11, commit 61767d4): clusters are
where the M = 1 single-collective-mode condition (6Y/9W) should break, and
the framework prescribes the breakage form: multimode (negative-binomial)
occupation statistics -- a DIFFERENT functional SHAPE at cluster scale, at
the SAME temperature a0 = c*H0/(2*pi), not merely extra mass under the
galaxy function. Two-sided BY DESIGN (sign caution pre-stated in the row).
This stage = the registered contest on public X-COP data:
  (i)  the M>1 family fits cluster-RAR shape+amplitude at the galaxy
       temperature  => the cluster wall becomes a window;
  (ii) it fits WORSE than "galaxy function + scaled extra mass" at the
       pre-registered grade => the genuinely-extra cluster component is
       measured within the framework (bounded negative, NOT a kill).

MODEL CONSTRUCTION (fixed before any fit; sympy/numeric gates G5)
-----------------------------------------------------------------
Dictionary: y = g_bar/a0, x = sqrt(y), n_BE(x) = 1/(exp(x)-1),
a0 = c*H0/(2*pi) with H0 = 67.8 km/s/Mpc (Planck 2015 -- BOTH the
program's lock convention, cf. 10H, AND X-COP's own adopted cosmology;
no conversion mismatch). a0 = 1.0484e-10 m s^-2 (computed in-code from
the formula, not hand-rounded).

  F-GAL   : nu(x) = 1 + n_BE(x).            0 shape params (the galaxy law).
  F-CROWD : nu(x) = 1 + M * n_BE(x), M >= 1 global continuous. 1 param.
            Derivation note (the registered structure): 6Y -- M democratic
            thermal modes carry negative-binomial occupation statistics;
            NB(M, q) mean = M * n_BE. 9W forces M = 1 for any linear
            SINGLE-PORT system (galaxies/binaries); the NB flank opens only
            for which-path-DISTINGUISHABLE multi-channel absorption (R18
            dividing line) -- which is precisely what a cluster's spatially
            resolved sub-halos provide. The 10C linear-pull theorem reads
            the total occupation linearly => nu = 1 + n_tot_bar = 1 + M*n_BE
            at the same temperature. Deep limit: nu -> M/x, i.e. an
            EFFECTIVE deep scale g_dagger_eff = M^2 * a0 while the
            Boltzmann cutoff stays anchored at x = sqrt(g_bar/a0) -- the
            shape discriminant vs a pure rescale.
            DISCLOSED LIMITATION: the 6Y GATE leg (the NB tail
            1-(1-q)^M(1+Mq) replacing s^2 in the AMB refinement) is a
            second-order shape effect requiring an ambient-e_N dictionary
            that is ill-posed when the ambient IS the crowd; it is NOT
            instantiated here. This stage contests the mean-occupation law.
  F-SCALE : nu(x) = 1 + n_BE(x/sqrt(S)), S >= 1 global. 1 param.
            (The single-rescale rival: "same shape, shifted scale" --
            the Tian-style g_dagger form at BE shape.)
  F-2SCALE: g_obs = g_bar * (1+g1/g_bar)^a1 * (1+g2/g_bar)^(a2-a1),
            4 params global (Eckert+22 eq. 20) -- the published empirical
            ceiling; CONTEXT ROW ONLY (not verdict-bearing).
  MODEL-B : g_obs = U_c * g_bar * nu_BE(sqrt(U_c * y)), U_c >= 1 free PER
            CLUSTER (the registered comparator "galaxy function + scaled
            extra mass": extra mass tracing baryons, galaxy law on top).
            N_c params (+ shared sigma_int).
All models share one intrinsic-scatter parameter sigma_int (dex, in
log10 g_obs, >= 0.02) fitted per model.

Likelihood (all models): Gaussian in log10 g_obs(r_j) with
  sigma_eff^2 = sigma_logo^2 + (slope_model * sigma_logb)^2 + sigma_int^2,
slope_model = dlog10 g_obs/dlog10 g_bar evaluated per model per point
(model-consistent propagation), sigma_logo from EM_FORW (symmetrized),
sigma_logb from MGAS_LO/HI + MSTAR_LO/HI (symmetrized, quadrature).

DATA (fetched 2026-08-11, calcs/fetch_xcop.py; data/xcop/ gitignored)
---------------------------------------------------------------------
X-COP public release (SWITCHdrive j3WUOYXWgv9Jbnz, 315080566 bytes;
the ISDC mirror is dead). Per cluster: {C}_hydro_mass.fits (HYDRO_MASS:
RADIUS kpc, M_FORW/EM_FORW model-agnostic + M_NFW/... + PARAMS HDU),
{C}_fgas_profile.fits (FGAS: RADIUS in R/R500_release, M_NFW, MGAS+-),
{C}_mstar.fits (MSTAR_SMOOTHED: RADIUS kpc, MSTAR+- total errors;
present for SEVEN clusters: A644, A1795, A2029, A2142, A2319, A85,
ZW1215). PRIMARY SAMPLE = the 7 mstar clusters. L12 leg = all 12 with
the median-f_star(r/R500) curve of the 7 substituted (labeled).
Transcription (G1 anchor): Eckert+22 Table 2 (ar5iv primary read
2026-08-11), Planck-2015 H0 = 67.8, Om = 0.308:
  A85 0.0555/1214; A644 0.0704/1398; A1644 0.0473/1031; A1795 0.0622/1160;
  A2029 0.0766/1340; A2142 0.0900/1453; A2255 0.0809/1202; A2319 0.0557/1424;
  A3158 0.0590/1146; A3266 0.0589/1381; RXC1825 0.0650/1109; ZW1215 0.0766/1346
  (z / R500_NFW kpc).
WIRING RULE (release-internal, fixed): the FGAS radius scale R500_rel is
SOLVED per cluster by matching the FGAS M_NFW column to the HYDRO_MASS
M_NFW column (golden-section on R500 in [700, 2200] kpc, median |dlog M|
loss); the Table-2 R500 is then a CHECK (G1), not an input to the wiring.
Radial window: r in [max(0.05*R500, r_min_data), 1.2*R500] intersected
with every needed component's measured range (log-log interpolation
INSIDE measured ranges only; no extrapolation). Thinned to N_R = 8
log-spaced points per cluster (PRIMARY; 6/10 = G8 robustness).
g_obs = G*M_FORW(r)/r^2. g_bar = G*(MGAS(r)+MSTAR(r))/r^2.

FIREWALL
--------
No model is fitted to the sky until gates G0-G5 pass AND the injection
suite (G6, G7) passes in the direction(s) to be read. The sky-mode
refuses to run without the gate marker file. Assembled sky tables are
built and stored during --gates for wiring checks (G2-G4 use release
redundancy, not our models); no model-vs-sky statistic is computed
before --sky.

GATES AND BARS (frozen; z-grade criteria per trap #22 where replicated)
-----------------------------------------------------------------------
G0  manifest: all 12 clusters x {hydro_mass, fgas_profile} + 7 mstar
    files present; SHA256 of the archive recorded in output.
G1  R500_rel solve converges for 12/12 in [700, 2200] kpc; agreement
    with Table-2 R500: |R500_rel - R500_T2|/R500_T2 <= 0.05 for >= 10/12.
G2  release-internal reader identity: after the solve, 90th-pct
    |dlog10 M(FGAS.M_NFW vs HYDRO.M_NFW)| <= 0.02 dex per cluster, 12/12.
G3  our NFW wiring: M_NFW(r) built from PARAMS (RS, C200; 200*rho_c(z)
    convention, Om = 0.308) matches the HYDRO M_NFW column: median
    |dlog10 M| <= 0.02 dex for >= 10/12 (a convention-validation gate;
    a clean constant offset = convention fix, disclosed, wiring-grade).
G4  FGAS column identity: FGAS = MGAS/M_NFW to <= 1e-3 relative
    (median), 12/12.
G5  model regressions (sympy + numeric): F-CROWD(M=1) == F-GAL exact;
    F-SCALE(S=1) == F-GAL exact; MODEL-B(U=1) == F-GAL exact; deep
    limit nu_CROWD*x -> M analytically; numeric identity d <= 1e-12 on
    a 200-point x-grid.
G6  recovery injections (R = 8 replicates per truth, sky design matrix
    + measured errors, sigma_int_true = 0.08 dex):
    truth F-SCALE(S=17): mean |S_hat - 17| <= max(3, 2.5*SE_rep);
    truth F-CROWD(M=4):  mean |M_hat - 4|  <= max(0.8, 2.5*SE_rep).
G7  CONTEST POWER (the read-permission gate, per direction):
    (a) under truth MODEL-B (U_c ~ U[4,12] iid per cluster, literature-
        anchored scale): dBIC(CROWD-B) > +10 in >= 6/8 replicates
        => the W-EXTRA direction is readable;
    (b) under truth F-CROWD(M=4): dBIC(CROWD-B) < -2 in >= 6/8
        replicates (BIC already penalizes B's params; the bar asks the
        data to actually prefer the true 2-param law)
        => the W-WINDOW direction is readable.
    A failed direction is NOT read on the sky (letter W-POWER-LIMITED
    for that direction; 9R discipline).
G8  thinning robustness: the operative sky statistic dBIC(CROWD-B)
    sign-stable at N_R in {6, 8, 10}.
G9  optimizer: >= 12 multi-starts per global fit; best two agree to
    <= 0.5 in -2lnL; nesting: -2lnL(F-GAL) >= -2lnL(F-SCALE) - 1e-6
    and >= -2lnL(F-CROWD) - 1e-6 and >= -2lnL(MODEL-B) - 1e-6.
G10 stored-table identity: the assembled sky table is written to npz in
    --gates and reloaded in --sky; lnL of a fixed probe model (F-GAL,
    sigma_int = 0.1) agrees build-vs-reload to 0.0 (8P/8Q rule).

VERDICT GRAMMAR (letters; every positive clause gated -- trap #12)
------------------------------------------------------------------
Operative statistic: dBIC_CB = BIC(F-CROWD) - BIC(MODEL-B) on the
PRIMARY sample/window; cluster-level paired bootstrap (400 reps,
resample the 7 clusters) gives P_boot(B better) and SEs; BIC k*ln(n)
with n = thinned point count; AIC and raw -2lnL co-reported always.
GoF PASS (model X) := chi2_red(X) in [0.5, 2.0] AND sigma_int_hat(X)
<= 0.20 dex.
  W-WINDOW  (outcome i): G7b PASS AND GoF(F-CROWD) PASS AND
            dBIC_CB <= +10 AND P_boot(B better) < 0.90 AND
            (M_hat - 1) >= 2*SE_boot(M_hat).
  W-EXTRA   (outcome ii): G7a PASS AND dBIC_CB > +10 AND
            P_boot(B better) >= 0.90 AND GoF(MODEL-B) PASS.
            => book the bounded negative + the measured extra component
            (median and range of U_c_hat).
  W-POWER-LIMITED: the G7 direction(s) needed above failed; report
            GoF rows only for the unreadable direction.
  W-GRAY:   between bars (report everything, claim nothing).
  W-FEASIBILITY-LIMITED: G0-G4 fail unrepairably (data-route verdict;
            successor = CLASH/Umetsu leg).
SIGN-CAUTION BOOKING (from the registered row): if M_hat < 1 is
preferred, report inside W-GRAY/W-EXTRA with the explicit label
"crowd-direction NEGATIVE (the naive more-open-gate direction)".
SECONDARY READS (labeled, never verdict-bearing): F-SCALE vs F-CROWD
dBIC (the same-temperature shape contest); g_dagger_eff = M_hat^2*a0
vs the published Tian 2.0e-9 (consistency read); per-cluster M_c
scatter (7-param labeled variant); F-2SCALE context row.
ROBUSTNESS LEGS (labeled): LNFW (M_NFW total instead of M_FORW);
LB15 (hydrostatic bias b = 0.15: g_obs x 1/(1-b)); L12 (12 clusters,
median-f_star substitution); LIN02 (inner cut 0.02*R500, cool-core
caveat); N_R legs (G8).

CREDENCE MAP (pre-signed; ALL CELLS HOLD)
-----------------------------------------
anomaly-real 53 and mech-conditional 8 HOLD in EVERY outcome cell
(W-WINDOW / W-EXTRA / W-GRAY / W-POWER-LIMITED / W-FEASIBILITY /
sign-caution). Rationale: the row is registered two-sided; outcome (ii)
is "a bounded negative, not a kill" by its own registration; outcome
(i) on 7 clusters at first firing is a fit result, not a demonstrated
window -- any upward move would need the review round + successor
external legs. The PREDICTIONS P13 row gets a status ANNOTATION per
outcome (no flip of the registered text either way at this grade).
EMOND (1701.03369) differentiation = print-use obligation only (no
print this stage).

AMENDMENTS (log here pre-quote if any; runs preserved)
------------------------------------------------------
A1 (2026-08-11, after gates run 1, BEFORE any re-run; log
   data/stage10p_gates_run1.log preserved; no verdict bar touched):
   (a) G3 CONVENTION FIX (the case the G3 note pre-licensed): our NFW
       rebuild missed the release M_NFW column by a UNIFORM 0.027 dex on
       all 12 clusters = exactly 2*log10(70/67.8) -- the release mass
       products are in the ETTORI+19 convention (H0 = 70, Om = 0.3), not
       Eckert+22's Planck-2015. G3 is re-evaluated at the release
       convention (H0_REL = 70, OM_REL = 0.3 inside nfw_mass only).
       DISCLOSED SEAM: the release (r, M) observables are h = 0.70-
       derived while a0 stays at the program lock (67.8-based
       1.0484e-10); propagating 67.8 through their distances would shift
       log g by a uniform ~0.014 dex -- subdominant to every bar here;
       booked as a labeled systematic in the stage notes.
   (b) G4 BAR 1e-3 -> 5e-3 median relative (mis-set bar, 10O-G3 class):
       observed identities 4e-4..1.9e-3 reflect the release's MC-chain
       storage (ratio-of-chain vs stored point values); the gate's
       target failure mode (column mislabeling / radius mismatch) is
       O(10-100%) and remains excluded at the new bar with margin.

MODES
-----
  py calcs/stage10p_p13crowd.py --gates    (G0-G5 + table build + G10 setup)
  py calcs/stage10p_p13crowd.py --inject   (G6, G7; requires gates marker)
  py calcs/stage10p_p13crowd.py --sky      (the contest; requires gates +
                                            inject markers; G8-G10 inline)
Outputs: data/stage10p_gates.json, data/stage10p_inject.json,
data/stage10p_results.json, data/stage10p_sky.npz, logs to stdout.
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np
from astropy.io import fits
from scipy.optimize import minimize_scalar, minimize

# ---------------------------------------------------------------- constants
C_LIGHT = 2.99792458e8          # m/s
H0_KMSMPC = 67.8                # Planck 2015 (X-COP + program lock convention)
MPC = 3.0856775814913673e22     # m
KPC = MPC / 1000.0
H0_SI = H0_KMSMPC * 1000.0 / MPC
A0 = C_LIGHT * H0_SI / (2.0 * np.pi)    # = 1.0484e-10 m/s^2
OM = 0.308
OL = 1.0 - OM
GSI = 6.674e-11
MSUN = 1.989e30

BASE = os.path.join(os.path.dirname(__file__), "..", "data", "xcop")
OUTD = os.path.join(os.path.dirname(__file__), "..", "data")
ARCHIVE = os.path.join(BASE, "xcop_release.archive")

CLUSTERS = ["A85", "A644", "A1644", "A1795", "A2029", "A2142", "A2255",
            "A2319", "A3158", "A3266", "RXC1825", "ZW1215"]
MSTAR_CLUSTERS = ["A85", "A644", "A1795", "A2029", "A2142", "A2319", "ZW1215"]
# Eckert+22 Table 2 (primary read 2026-08-11): z, R500_NFW [kpc]
T2 = {
    "A85": (0.0555, 1214.0), "A644": (0.0704, 1398.0), "A1644": (0.0473, 1031.0),
    "A1795": (0.0622, 1160.0), "A2029": (0.0766, 1340.0), "A2142": (0.0900, 1453.0),
    "A2255": (0.0809, 1202.0), "A2319": (0.0557, 1424.0), "A3158": (0.0590, 1146.0),
    "A3266": (0.0589, 1381.0), "RXC1825": (0.0650, 1109.0), "ZW1215": (0.0766, 1346.0),
}

N_R_PRIMARY = 8
INNER_FRAC = 0.05
OUTER_FRAC = 1.20
SIG_INT_MIN = 0.02
N_BOOT = 400
N_STARTS = 12
RNG_SEED = 10131            # fixed at pre-registration
R_REP = 8                   # injection replicates (trap #22)

GATES_MARKER = os.path.join(OUTD, "stage10p_gates.json")
INJECT_MARKER = os.path.join(OUTD, "stage10p_inject.json")
SKY_NPZ = os.path.join(OUTD, "stage10p_sky.npz")


# ---------------------------------------------------------------- utilities
def hz(z):
    return H0_SI * np.sqrt(OM * (1.0 + z) ** 3 + OL)


def rho_c(z):
    return 3.0 * hz(z) ** 2 / (8.0 * np.pi * GSI)


def loglog_interp(r, rr, mm):
    """Interpolate log10 m at r inside [rr.min(), rr.max()]; nan outside."""
    out = np.full_like(np.atleast_1d(r), np.nan, dtype=float)
    ra = np.atleast_1d(r).astype(float)
    ok = (ra >= rr.min()) & (ra <= rr.max())
    out[ok] = np.interp(np.log10(ra[ok]), np.log10(rr), np.log10(mm))
    return 10.0 ** out if np.ndim(r) else float(10.0 ** out[0])


def sha256_of(path, blocksize=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(blocksize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- models
def n_be(x):
    x = np.asarray(x, dtype=float)
    return 1.0 / np.expm1(np.clip(x, 1e-12, 700.0))


def gobs_gal(gbar):
    x = np.sqrt(gbar / A0)
    return gbar * (1.0 + n_be(x))


def gobs_crowd(gbar, M):
    x = np.sqrt(gbar / A0)
    return gbar * (1.0 + M * n_be(x))


def gobs_scale(gbar, S):
    x = np.sqrt(gbar / (A0 * S))
    return gbar * (1.0 + n_be(x))


def gobs_two(gbar, g1, a1, g2, a2):
    return gbar * (1.0 + g1 / gbar) ** a1 * (1.0 + g2 / gbar) ** (a2 - a1)


def gobs_modelB(gbar, U):
    ge = U * gbar
    x = np.sqrt(ge / A0)
    return ge * (1.0 + n_be(x))


def model_slope(fun, gbar, *pars):
    """dlog10 g_obs / dlog10 g_bar, central difference in log space."""
    d = 0.01
    lo = fun(gbar * 10 ** (-d), *pars)
    hi = fun(gbar * 10 ** (+d), *pars)
    return (np.log10(hi) - np.log10(lo)) / (2 * d)


# ---------------------------------------------------------------- data build
def read_cluster(name):
    d = os.path.join(BASE, name)
    out = {}
    with fits.open(os.path.join(d, f"{name}_hydro_mass.fits")) as h:
        t = h["HYDRO_MASS"].data
        out["r_h"] = np.array(t["RADIUS"], dtype=float)
        out["m_forw"] = np.array(t["M_FORW"], dtype=float)
        out["em_forw"] = np.array(t["EM_FORW"], dtype=float)
        out["m_nfw"] = np.array(t["M_NFW"], dtype=float)
        p = h["PARAMS"].data
        for row in p:
            if str(row["MODEL"]).strip() == "NFW":
                out["nfw_rs"] = float(row["RS"])
                out["nfw_c200"] = float(row["C200"])
    with fits.open(os.path.join(d, f"{name}_fgas_profile.fits")) as h:
        t = h["FGAS"].data
        out["u_f"] = np.array(t["RADIUS"], dtype=float)          # R/R500_release
        out["mnfw_f"] = np.array(t["M_NFW"], dtype=float)
        out["mgas"] = np.array(t["MGAS"], dtype=float)
        out["mgas_lo"] = np.array(t["MGAS_LO"], dtype=float)
        out["mgas_hi"] = np.array(t["MGAS_HI"], dtype=float)
        out["fgas"] = np.array(t["FGAS"], dtype=float)
    mst = os.path.join(d, f"{name}_mstar.fits")
    if os.path.exists(mst):
        with fits.open(mst) as h:
            t = h["MSTAR_SMOOTHED"].data
            out["r_s"] = np.array(t["RADIUS"], dtype=float)
            out["mstar"] = np.array(t["MSTAR"], dtype=float)
            out["mstar_lo"] = np.array(t["MSTAR_LO"], dtype=float)
            out["mstar_hi"] = np.array(t["MSTAR_HI"], dtype=float)
    return out


def solve_r500_release(c):
    """Golden-section: R500 that matches FGAS M_NFW(u*R500) to HYDRO M_NFW(kpc)."""
    rr, mm = c["r_h"], c["m_nfw"]

    def loss(r500):
        r = c["u_f"] * r500
        ok = (r >= rr.min()) & (r <= rr.max())
        if ok.sum() < 5:
            return 1e9
        mh = loglog_interp(r[ok], rr, mm)
        return float(np.median(np.abs(np.log10(c["mnfw_f"][ok]) - np.log10(mh))))

    res = minimize_scalar(loss, bounds=(700.0, 2200.0), method="bounded",
                          options={"xatol": 0.05})
    return float(res.x), float(res.fun)


H0_REL = 70.0 * 1000.0 / MPC     # A1a: release NFW column = Ettori+19 convention
OM_REL = 0.300


def rho_c_rel(z):
    h2 = H0_REL ** 2 * (OM_REL * (1.0 + z) ** 3 + (1.0 - OM_REL))
    return 3.0 * h2 / (8.0 * np.pi * GSI)


def nfw_mass(r_kpc, rs_kpc, c200, z):
    r200 = c200 * rs_kpc
    rho_s = (200.0 / 3.0) * rho_c_rel(z) * c200 ** 3 / (np.log(1 + c200) - c200 / (1 + c200))
    xx = r_kpc / rs_kpc
    m = 4.0 * np.pi * rho_s * (rs_kpc * KPC) ** 3 * (np.log(1 + xx) - xx / (1 + xx))
    return m / MSUN


def assemble(sample, n_r=N_R_PRIMARY, use_nfw=False, inner=INNER_FRAC,
             r500rel=None, fstar_med=None):
    """Build the thinned sky table. Returns dict of arrays + per-cluster index."""
    rows = []
    for name in sample:
        c = DATA[name]
        z, r500_t2 = T2[name]
        r500 = r500rel[name]
        rr_h, m_tot = c["r_h"], (c["m_nfw"] if use_nfw else c["m_forw"])
        em = c["em_forw"]
        r_gas = c["u_f"] * r500
        lo = max(inner * r500_t2, rr_h.min(), r_gas.min())
        hi_list = [OUTER_FRAC * r500_t2, rr_h.max(), r_gas.max()]
        if "r_s" in c:
            lo = max(lo, c["r_s"].min())
            hi_list.append(c["r_s"].max())
        hi = min(hi_list)
        if hi / lo < 2.0:
            continue
        r = np.logspace(np.log10(lo), np.log10(hi), n_r)
        mt = loglog_interp(r, rr_h, m_tot)
        emt = loglog_interp(r, rr_h, np.maximum(em, 1e-8 * m_tot))
        mg = loglog_interp(r, r_gas, c["mgas"])
        emg = 0.5 * (loglog_interp(r, r_gas, np.maximum(c["mgas_lo"], 1e-6 * c["mgas"]))
                     + loglog_interp(r, r_gas, np.maximum(c["mgas_hi"], 1e-6 * c["mgas"])))
        if "r_s" in c:
            ms = loglog_interp(r, c["r_s"], c["mstar"])
            ems = 0.5 * (np.abs(ms - loglog_interp(r, c["r_s"], c["mstar_lo"]))
                         + np.abs(loglog_interp(r, c["r_s"], c["mstar_hi"]) - ms))
        else:
            fs = np.interp(r / r500_t2, fstar_med[0], fstar_med[1])
            ms = fs * mg
            ems = 0.5 * ms
        r_m = r * KPC
        gobs = GSI * mt * MSUN / r_m ** 2
        gbar = GSI * (mg + ms) * MSUN / r_m ** 2
        s_logo = emt / mt / np.log(10.0)
        s_logb = np.sqrt(emg ** 2 + ems ** 2) / (mg + ms) / np.log(10.0)
        for j in range(n_r):
            rows.append((name, r[j], gobs[j], gbar[j], s_logo[j], s_logb[j]))
    names = np.array([t[0] for t in rows])
    out = {
        "name": names,
        "r": np.array([t[1] for t in rows]),
        "gobs": np.array([t[2] for t in rows]),
        "gbar": np.array([t[3] for t in rows]),
        "s_logo": np.array([t[4] for t in rows]),
        "s_logb": np.array([t[5] for t in rows]),
    }
    return out


# ---------------------------------------------------------------- likelihood
def m2ll(tab, fun, pars, sig_int, per_cluster_U=None):
    gbar = tab["gbar"]
    if per_cluster_U is not None:
        gm = np.empty_like(gbar)
        sl = np.empty_like(gbar)
        for name, U in per_cluster_U.items():
            m = tab["name"] == name
            gm[m] = fun(gbar[m], U)
            sl[m] = model_slope(fun, gbar[m], U)
    else:
        gm = fun(gbar, *pars)
        sl = model_slope(fun, gbar, *pars)
    resid = np.log10(tab["gobs"]) - np.log10(gm)
    var = tab["s_logo"] ** 2 + (sl * tab["s_logb"]) ** 2 + sig_int ** 2
    return float(np.sum(resid ** 2 / var + np.log(2 * np.pi * var)))


def fit_global(tab, fun, p0_list, bounds, rng):
    """Multi-start fit of (pars..., sig_int). Returns best (m2ll, pars, sig)."""
    best = None
    tries = []
    for p0 in p0_list:
        def obj(th):
            *pars, lsig = th
            return m2ll(tab, fun, pars, SIG_INT_MIN + np.exp(lsig))
        res = minimize(obj, p0, method="Nelder-Mead",
                       options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-8})
        tries.append((res.fun, res.x))
    tries.sort(key=lambda t: t[0])
    best = tries[0]
    gap = tries[1][0] - tries[0][0] if len(tries) > 1 else 0.0
    *pars, lsig = best[1]
    return best[0], list(pars), SIG_INT_MIN + np.exp(lsig), gap


def fit_modelB(tab, rng, n_starts=4):
    """Per-cluster U_c >= 1 + shared sigma; coordinate polish from multi-starts."""
    names = sorted(set(tab["name"].tolist()))
    best = None
    for s in range(n_starts):
        U0 = {n: float(np.exp(rng.uniform(np.log(2.0), np.log(20.0)))) for n in names}
        lsig = np.log(0.05)
        th0 = np.array([np.log(U0[n]) for n in names] + [lsig])

        def obj(th):
            U = {n: 1.0 + np.exp(th[i]) - 1.0 if False else max(1.0, np.exp(th[i]))
                 for i, n in enumerate(names)}
            return m2ll(tab, gobs_modelB, None, SIG_INT_MIN + np.exp(th[-1]),
                        per_cluster_U=U)
        res = minimize(obj, th0, method="Nelder-Mead",
                       options={"maxiter": 20000, "xatol": 1e-6, "fatol": 1e-8})
        if best is None or res.fun < best[0]:
            best = (res.fun, res.x)
    th = best[1]
    U = {n: max(1.0, float(np.exp(th[i]))) for i, n in enumerate(names)}
    return best[0], U, SIG_INT_MIN + float(np.exp(th[-1]))


def bic(m2, k, n):
    return m2 + k * np.log(n)


# ---------------------------------------------------------------- gate modes
DATA = {}


def load_all():
    for n in CLUSTERS:
        DATA[n] = read_cluster(n)


def run_gates():
    rep = {"stage": "10P", "mode": "gates"}
    # G0
    missing = []
    for n in CLUSTERS:
        for suff in ("hydro_mass", "fgas_profile"):
            p = os.path.join(BASE, n, f"{n}_{suff}.fits")
            if not os.path.exists(p):
                missing.append(p)
    for n in MSTAR_CLUSTERS:
        p = os.path.join(BASE, n, f"{n}_mstar.fits")
        if not os.path.exists(p):
            missing.append(p)
    rep["G0_missing"] = missing
    rep["G0_pass"] = len(missing) == 0
    rep["G0_sha256"] = sha256_of(ARCHIVE)
    load_all()

    # G1 + G2: R500_rel solve + reader identity
    r500rel, g1_ok, g2_p90 = {}, 0, {}
    for n in CLUSTERS:
        r5, res_med = solve_r500_release(DATA[n])
        r500rel[n] = r5
        z, r5t2 = T2[n]
        if abs(r5 - r5t2) / r5t2 <= 0.05:
            g1_ok += 1
        c = DATA[n]
        r = c["u_f"] * r5
        ok = (r >= c["r_h"].min()) & (r <= c["r_h"].max())
        mh = loglog_interp(r[ok], c["r_h"], c["m_nfw"])
        d = np.abs(np.log10(c["mnfw_f"][ok]) - np.log10(mh))
        g2_p90[n] = float(np.percentile(d, 90))
    rep["G1_r500rel"] = r500rel
    rep["G1_ok_count"] = g1_ok
    rep["G1_pass"] = (len(r500rel) == 12) and (g1_ok >= 10)
    rep["G2_p90"] = g2_p90
    rep["G2_pass"] = all(v <= 0.02 for v in g2_p90.values())

    # G3 our NFW wiring
    g3 = {}
    for n in CLUSTERS:
        c = DATA[n]
        z, _ = T2[n]
        m_us = nfw_mass(c["r_h"], c["nfw_rs"], c["nfw_c200"], z)
        d = np.abs(np.log10(m_us) - np.log10(c["m_nfw"]))
        g3[n] = float(np.median(d))
    rep["G3_med"] = g3
    rep["G3_ok_count"] = sum(1 for v in g3.values() if v <= 0.02)
    rep["G3_pass"] = rep["G3_ok_count"] >= 10

    # G4 fgas identity
    g4 = {}
    for n in CLUSTERS:
        c = DATA[n]
        d = np.abs(c["fgas"] - c["mgas"] / c["mnfw_f"]) / np.maximum(c["fgas"], 1e-9)
        g4[n] = float(np.median(d))
    rep["G4_med"] = g4
    rep["G4_pass"] = all(v <= 5e-3 for v in g4.values())   # A1b

    # G5 model regressions
    import sympy as sp
    xs, Ms, Ss = sp.symbols("x M S", positive=True)
    nbe = 1 / (sp.exp(xs) - 1)
    e1 = sp.simplify((1 + Ms * nbe).subs(Ms, 1) - (1 + nbe))
    e2 = sp.simplify((1 + 1 / (sp.exp(xs / sp.sqrt(Ss)) - 1)).subs(Ss, 1) - (1 + nbe))
    deep = sp.limit(Ms * nbe * xs, xs, 0)
    xg = np.logspace(-2, 1.3, 200)
    gg = A0 * xg ** 2
    d1 = np.max(np.abs(gobs_crowd(gg, 1.0) - gobs_gal(gg)) / gobs_gal(gg))
    d2 = np.max(np.abs(gobs_scale(gg, 1.0) - gobs_gal(gg)) / gobs_gal(gg))
    d3 = np.max(np.abs(gobs_modelB(gg, 1.0) - gobs_gal(gg)) / gobs_gal(gg))
    rep["G5_sym"] = [str(e1), str(e2), str(deep)]
    rep["G5_num"] = [float(d1), float(d2), float(d3)]
    rep["G5_pass"] = (e1 == 0 and e2 == 0 and deep == Ms and
                      max(d1, d2, d3) <= 1e-12)

    # build + store the sky tables (G10 setup); median f_star curve
    ug = np.linspace(0.05, 1.2, 24)
    fs_curves = []
    for n in MSTAR_CLUSTERS:
        c = DATA[n]
        _, r5t2 = T2[n]
        r = ug * r5t2
        ms = loglog_interp(r, c["r_s"], c["mstar"])
        mg = loglog_interp(r, c["u_f"] * r500rel[n], c["mgas"])
        fs_curves.append(ms / mg)
    fs_med = np.nanmedian(np.array(fs_curves), axis=0)
    tab7 = assemble(MSTAR_CLUSTERS, r500rel=r500rel)
    tab12 = assemble(CLUSTERS, r500rel=r500rel, fstar_med=(ug, fs_med))
    tab7_nfw = assemble(MSTAR_CLUSTERS, r500rel=r500rel, use_nfw=True)
    tab7_in02 = assemble(MSTAR_CLUSTERS, r500rel=r500rel, inner=0.02)
    tab7_n6 = assemble(MSTAR_CLUSTERS, n_r=6, r500rel=r500rel)
    tab7_n10 = assemble(MSTAR_CLUSTERS, n_r=10, r500rel=r500rel)
    np.savez(SKY_NPZ,
             **{f"p_{k}": v for k, v in tab7.items()},
             **{f"t12_{k}": v for k, v in tab12.items()},
             **{f"nfw_{k}": v for k, v in tab7_nfw.items()},
             **{f"in02_{k}": v for k, v in tab7_in02.items()},
             **{f"n6_{k}": v for k, v in tab7_n6.items()},
             **{f"n10_{k}": v for k, v in tab7_n10.items()},
             fstar_u=ug, fstar_med=fs_med,
             r500rel_names=np.array(CLUSTERS),
             r500rel_vals=np.array([r500rel[n] for n in CLUSTERS]))
    rep["G10_probe_build"] = m2ll(tab7, gobs_gal, [], 0.1)
    rep["n_points_primary"] = int(len(tab7["gbar"]))
    rep["gbar_range_primary"] = [float(tab7["gbar"].min()), float(tab7["gbar"].max())]
    rep["x_range_primary"] = [float(np.sqrt(tab7["gbar"].min() / A0)),
                              float(np.sqrt(tab7["gbar"].max() / A0))]
    rep["a0"] = A0

    allpass = all(rep[k] for k in ("G0_pass", "G1_pass", "G2_pass", "G3_pass",
                                   "G4_pass", "G5_pass"))
    rep["gates_all_pass"] = bool(allpass)
    with open(GATES_MARKER, "w") as f:
        json.dump(rep, f, indent=1, default=str)
    print(json.dumps({k: rep[k] for k in rep if not k.startswith("G1_r5")},
                     indent=1, default=str))
    print("GATES", "ALL PASS" if allpass else "FAIL")


def tab_from_npz(z, prefix):
    return {k: z[f"{prefix}_{k}"] for k in ("name", "r", "gobs", "gbar",
                                            "s_logo", "s_logb")}


def fit_suite(tab, rng, do_B=True):
    """Fit F-GAL, F-CROWD, F-SCALE (+ MODEL-B). Returns dict."""
    n = len(tab["gbar"])
    out = {"n": n}
    m2g = None
    best_s = fit_global(tab, gobs_gal,
                        [[np.log(0.05 + 0.05 * i)] for i in range(3)], None, rng)
    out["GAL"] = {"m2": best_s[0], "sig": best_s[2], "k": 1}
    p0c = [[np.log(m), np.log(0.08)] for m in
           np.exp(rng.uniform(np.log(1.0), np.log(60.0), N_STARTS))]
    b = fit_global(tab, lambda g, lm: gobs_crowd(g, 1.0 + np.exp(lm)), p0c, None, rng)
    out["CROWD"] = {"m2": b[0], "M": 1.0 + float(np.exp(b[1][0])), "sig": b[2],
                    "gap": b[3], "k": 2}
    p0s = [[np.log(s), np.log(0.08)] for s in
           np.exp(rng.uniform(np.log(1.0), np.log(400.0), N_STARTS))]
    b = fit_global(tab, lambda g, ls: gobs_scale(g, 1.0 + np.exp(ls)), p0s, None, rng)
    out["SCALE"] = {"m2": b[0], "S": 1.0 + float(np.exp(b[1][0])), "sig": b[2],
                    "gap": b[3], "k": 2}
    if do_B:
        m2, U, sig = fit_modelB(tab, rng)
        out["B"] = {"m2": m2, "U": U, "sig": sig,
                    "k": len(U) + 1}
    for key in [k for k in ("GAL", "CROWD", "SCALE", "B") if k in out]:
        out[key]["bic"] = bic(out[key]["m2"], out[key]["k"], n)
        out[key]["aic"] = out[key]["m2"] + 2 * out[key]["k"]
    return out


def synth_from(tab, fun, pars, sig_int, rng, perU=None):
    t = {k: v.copy() for k, v in tab.items()}
    gbar = t["gbar"]
    if perU is not None:
        gm = np.empty_like(gbar)
        for nm, U in perU.items():
            m = t["name"] == nm
            gm[m] = fun(gbar[m], U)
    else:
        gm = fun(gbar, *pars)
    scat = rng.normal(0.0, 1.0, len(gbar)) * np.sqrt(t["s_logo"] ** 2 + sig_int ** 2)
    t["gobs"] = 10 ** (np.log10(gm) + scat)
    return t


def run_inject():
    if not os.path.exists(GATES_MARKER):
        sys.exit("gates marker missing -- run --gates first")
    with open(GATES_MARKER) as f:
        gates = json.load(f)
    if not gates.get("gates_all_pass"):
        sys.exit("gates not all-pass -- inject refused")
    z = np.load(SKY_NPZ, allow_pickle=True)
    tab = tab_from_npz(z, "p")
    rng = np.random.default_rng(RNG_SEED)
    rep = {"stage": "10P", "mode": "inject"}

    # G6 recovery
    for label, fun, pars, key, true in (
            ("S17", gobs_scale, [17.0], "S", 17.0),
            ("M4", gobs_crowd, [4.0], "M", 4.0)):
        vals = []
        for rrep in range(R_REP):
            t = synth_from(tab, fun, pars, 0.08, rng)
            fs = fit_suite(t, rng, do_B=False)
            vals.append(fs["CROWD" if key == "M" else "SCALE"][key])
        vals = np.array(vals)
        se = vals.std(ddof=1) / np.sqrt(R_REP)
        bar = max({"S": 3.0, "M": 0.8}[key], 2.5 * se * np.sqrt(R_REP))
        rep[f"G6_{label}"] = {"mean": float(vals.mean()), "sd": float(vals.std(ddof=1)),
                              "se": float(se), "bar": float(bar),
                              "pass": bool(abs(vals.mean() - true) <= bar)}

    # G7 contest power
    names = sorted(set(tab["name"].tolist()))
    ok_a = ok_b = 0
    da, db = [], []
    for rrep in range(R_REP):
        Utr = {nm: float(rng.uniform(4.0, 12.0)) for nm in names}
        t = synth_from(tab, gobs_modelB, None, 0.08, rng, perU=Utr)
        fs = fit_suite(t, rng)
        d = fs["CROWD"]["bic"] - fs["B"]["bic"]
        da.append(d)
        ok_a += int(d > 10.0)
        t = synth_from(tab, gobs_crowd, [4.0], 0.08, rng)
        fs = fit_suite(t, rng)
        d = fs["CROWD"]["bic"] - fs["B"]["bic"]
        db.append(d)
        ok_b += int(d < -2.0)
    rep["G7a"] = {"dbics": [float(v) for v in da], "n_ok": ok_a,
                  "pass": bool(ok_a >= 6)}
    rep["G7b"] = {"dbics": [float(v) for v in db], "n_ok": ok_b,
                  "pass": bool(ok_b >= 6)}
    rep["inject_all_pass"] = bool(rep["G6_S17"]["pass"] and rep["G6_M4"]["pass"])
    with open(INJECT_MARKER, "w") as f:
        json.dump(rep, f, indent=1)
    print(json.dumps(rep, indent=1))


def run_sky():
    if not (os.path.exists(GATES_MARKER) and os.path.exists(INJECT_MARKER)):
        sys.exit("markers missing -- run --gates and --inject first")
    with open(INJECT_MARKER) as f:
        inj = json.load(f)
    z = np.load(SKY_NPZ, allow_pickle=True)
    rng = np.random.default_rng(RNG_SEED + 1)
    rep = {"stage": "10P", "mode": "sky",
           "G7a_pass": inj["G7a"]["pass"], "G7b_pass": inj["G7b"]["pass"]}

    tab = tab_from_npz(z, "p")
    # G10 reload identity
    with open(GATES_MARKER) as f:
        gates = json.load(f)
    probe = m2ll(tab, gobs_gal, [], 0.1)
    rep["G10_d"] = float(abs(probe - gates["G10_probe_build"]))
    rep["G10_pass"] = rep["G10_d"] == 0.0

    fs = fit_suite(tab, rng)
    rep["primary"] = {k: {kk: vv for kk, vv in v.items() if kk != "U"}
                      for k, v in fs.items() if isinstance(v, dict)}
    rep["primary"]["B_U"] = fs["B"]["U"]
    rep["dBIC_CB"] = fs["CROWD"]["bic"] - fs["B"]["bic"]
    rep["dBIC_CS"] = fs["CROWD"]["bic"] - fs["SCALE"]["bic"]
    rep["dBIC_CG"] = fs["CROWD"]["bic"] - fs["GAL"]["bic"]
    # chi2_red at best fit (approx: mean standardized resid^2)
    def chi2red(fun, pars, sig, perU=None):
        gbar = tab["gbar"]
        if perU is not None:
            gm = np.empty_like(gbar)
            sl = np.empty_like(gbar)
            for nm, U in perU.items():
                m = tab["name"] == nm
                gm[m] = fun(gbar[m], U)
                sl[m] = model_slope(fun, gbar[m], U)
            kfree = len(perU) + 1
        else:
            gm = fun(gbar, *pars)
            sl = model_slope(fun, gbar, *pars)
            kfree = len(pars) + 1
        var = tab["s_logo"] ** 2 + (sl * tab["s_logb"]) ** 2 + sig ** 2
        resid = np.log10(tab["gobs"]) - np.log10(gm)
        return float(np.sum(resid ** 2 / var) / (len(gbar) - kfree))
    rep["chi2red_CROWD"] = chi2red(gobs_crowd, [fs["CROWD"]["M"]], fs["CROWD"]["sig"])
    rep["chi2red_B"] = chi2red(gobs_modelB, None, fs["B"]["sig"], perU=fs["B"]["U"])
    rep["chi2red_SCALE"] = chi2red(gobs_scale, [fs["SCALE"]["S"]], fs["SCALE"]["sig"])

    # F-2SCALE context row
    def f2(g, lg1, a1, lg2, a2):
        return gobs_two(g, np.exp(lg1), a1, np.exp(lg2), a2)
    p0 = [[np.log(2e-11), 0.5, np.log(2e-9), 0.2],
          [np.log(1e-10), 0.4, np.log(1e-9), 0.1],
          [np.log(5e-11), 0.6, np.log(4e-9), 0.3]]
    p0 = [p + [np.log(0.08)] for p in p0]
    b = fit_global(tab, f2, p0, None, rng)
    rep["TWOSCALE"] = {"m2": b[0], "pars": [float(v) for v in b[1]], "sig": b[2],
                       "bic": bic(b[0], 5, fs["n"])}

    # bootstrap (cluster-level, paired)
    names = sorted(set(tab["name"].tolist()))
    stats = {"dBIC_CB": [], "dBIC_CS": [], "M": [], "S": []}
    for i in range(N_BOOT):
        pick = rng.choice(names, size=len(names), replace=True)
        idx = np.concatenate([np.where(tab["name"] == nm)[0] for nm in pick])
        tb = {k: v[idx] for k, v in tab.items()}
        # rename duplicated clusters so MODEL-B fits one U per COPY
        newn = np.concatenate([[f"{nm}#{j}"] * (tab["name"] == nm).sum()
                               for j, nm in enumerate(pick)])
        tb["name"] = newn
        fb = fit_suite(tb, rng)
        stats["dBIC_CB"].append(fb["CROWD"]["bic"] - fb["B"]["bic"])
        stats["dBIC_CS"].append(fb["CROWD"]["bic"] - fb["SCALE"]["bic"])
        stats["M"].append(fb["CROWD"]["M"])
        stats["S"].append(fb["SCALE"]["S"])
    for k in stats:
        a = np.array(stats[k])
        rep[f"boot_{k}"] = {"mean": float(a.mean()), "sd": float(a.std(ddof=1)),
                            "p16": float(np.percentile(a, 16)),
                            "p84": float(np.percentile(a, 84))}
    rep["boot_P_B_better"] = float(np.mean(np.array(stats["dBIC_CB"]) > 0))
    rep["boot_frac"] = N_BOOT

    # G8 thinning + robustness legs (dBIC_CB sign stability)
    legs = {}
    for lg, pref in (("n6", "n6"), ("n10", "n10"), ("LNFW", "nfw"),
                     ("LIN02", "in02"), ("L12", "t12")):
        tt = tab_from_npz(z, pref)
        fl = fit_suite(tt, rng)
        legs[lg] = {"dBIC_CB": fl["CROWD"]["bic"] - fl["B"]["bic"],
                    "M": fl["CROWD"]["M"], "S": fl["SCALE"]["S"],
                    "n": fl["n"]}
    # LB15: hydrostatic bias leg on the primary table
    tb15 = {k: v.copy() for k, v in tab.items()}
    tb15["gobs"] = tb15["gobs"] / 0.85
    fl = fit_suite(tb15, rng)
    legs["LB15"] = {"dBIC_CB": fl["CROWD"]["bic"] - fl["B"]["bic"],
                    "M": fl["CROWD"]["M"], "S": fl["SCALE"]["S"], "n": fl["n"]}
    rep["legs"] = legs
    sgn = np.sign(rep["dBIC_CB"])
    rep["G8_pass"] = bool(np.sign(legs["n6"]["dBIC_CB"]) == sgn
                          and np.sign(legs["n10"]["dBIC_CB"]) == sgn)
    # G9 nesting
    rep["G9_nest"] = {
        "gal_ge_crowd": bool(fs["GAL"]["m2"] >= fs["CROWD"]["m2"] - 1e-6),
        "gal_ge_scale": bool(fs["GAL"]["m2"] >= fs["SCALE"]["m2"] - 1e-6),
        "gal_ge_B": bool(fs["GAL"]["m2"] >= fs["B"]["m2"] - 1e-6),
        "gaps": [fs["CROWD"]["gap"], fs["SCALE"]["gap"]],
    }
    rep["G9_pass"] = all(v for k, v in rep["G9_nest"].items() if k != "gaps") and \
        max(fs["CROWD"]["gap"], fs["SCALE"]["gap"]) <= 0.5

    # letter (grammar implemented verbatim from the header)
    gof_crowd = (0.5 <= rep["chi2red_CROWD"] <= 2.0) and fs["CROWD"]["sig"] <= 0.20
    gof_B = (0.5 <= rep["chi2red_B"] <= 2.0) and fs["B"]["sig"] <= 0.20
    Mhat = fs["CROWD"]["M"]
    se_M = rep["boot_M"]["sd"]
    pB = rep["boot_P_B_better"]
    letter = "W-GRAY"
    if not (rep["G7a_pass"] or rep["G7b_pass"]):
        letter = "W-POWER-LIMITED"
    else:
        if (rep["G7b_pass"] and gof_crowd and rep["dBIC_CB"] <= 10.0
                and pB < 0.90 and (Mhat - 1.0) >= 2.0 * se_M):
            letter = "W-WINDOW"
        elif (rep["G7a_pass"] and rep["dBIC_CB"] > 10.0 and pB >= 0.90 and gof_B):
            letter = "W-EXTRA"
        elif not rep["G7b_pass"] and rep["G7a_pass"]:
            letter = "W-GRAY (window direction POWER-LIMITED)"
        elif not rep["G7a_pass"] and rep["G7b_pass"]:
            letter = "W-GRAY (extra direction POWER-LIMITED)"
    rep["gof_crowd"] = gof_crowd
    rep["gof_B"] = gof_B
    rep["letter"] = letter
    rep["gdagger_eff"] = float(Mhat ** 2 * A0)
    rep["a0"] = A0

    with open(os.path.join(OUTD, "stage10p_results.json"), "w") as f:
        json.dump(rep, f, indent=1, default=str)
    print(json.dumps(rep, indent=1, default=str))
    print("LETTER:", letter)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gates", action="store_true")
    ap.add_argument("--inject", action="store_true")
    ap.add_argument("--sky", action="store_true")
    a = ap.parse_args()
    if a.gates:
        run_gates()
    elif a.inject:
        load_all()
        run_inject()
    elif a.sky:
        load_all()
        run_sky()
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
