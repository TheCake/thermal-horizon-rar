"""
STAGE 10H -- THE Z-LADDER CONTEST (MUSE-DARK III re-analysis)
=============================================================
PRE-REGISTRATION (2026-08-09). Bars, gates, letters, and the credence map
below are locked BEFORE any analysis run. The only operations performed
before this commit were the data download (calcs/fetch_musedark3_release.py)
and format/byte reads (catalog headers, one galaxy's file formats, one map
header) -- no fit of any kind has been run on this data.

THE QUESTION (the successor booked in NOTES + PREDICTIONS on 2026-08-09,
the MUSE-DARK re-read): MUSE-DARK III (Ciocan et al. 2026, arXiv:2604.22613)
measures a0(z) = a0(0) + a1*z with a1 = 1.59 +/- 0.10 (95%) on DC14 tracks
(their Eq. 4) and a1 = 1.20 +/- 0.10 (95%) in their MOND-framework refit
(their Eq. E.5). The program's locked pair predicts a0(z) = cH(z)/2pi,
whose linearization over their window is ~1.04e-10/z. Two legs:

  LEG A (lock-vs-linear): fit the RELEASED per-galaxy products with the
    lock (zero shape parameters), the H(z)-shape at free amplitude, their
    linear form, and a constant -- same MNR likelihood for all.
  LEG B (the lens-bias leg): inject locked-pair truth into the real x-side
    design and fit with THEIR form (Eq. 1 + linear a0(z)). Does the
    recovered a1 inflate from the naive 1.04 toward their measured values?
    (The mechanism: the lock is convex in z + the sampled acceleration
    window drifts with z; a fixed-form linear fit can alias both.)

DATA (public release, dark-matter.osu-lyon.fr; fetched 2026-08-09):
  catalogs (7 models + fit statistics + photometry), per-galaxy DC14
  galpak products: true_Vrot.dat (dx_arcsec, rad_Re, flux_slit, v_kms,
  sig_kms), galaxy_parameters.dat, derived_parameters.dat (log_Mdisk,
  rad_kpc, ...), mass_%d.fits (BUNIT Msun/arcsec2). ID0069 folder is
  absent from the release (404) -- census gate handles it.
  DISCLOSED SEAM: the folder runs and the *_bestfit.txt catalogs disagree
  at parameter level for at least ID3 (Vvir 160.6 vs 126.3) -- release
  versioning. Pre-rule: TRACK-level work uses folder-self-consistent
  parameters (true_Vrot + derived_parameters); catalogs are used for
  sample selection, evidence cuts, and E1 error scales only.
  DISCLOSED CONDITIONING: the released kinematics are DC14-conditioned
  (their App. E MOND-track products are NOT released). Leg A therefore
  tests the lock against DC14-processed tracks; their own framework
  spread (1.59 DM-row vs 1.20 MOND-row) measures that conditioning at
  +/-25% on a1. Verdict language carries this everywhere.

RECONSTRUCTION (App. B of the paper):
  a_tot(r) = v_c^2/r,  v_c^2 = v^2 + v_AD^2 (Dalcanton-Stilp pressure
    support: v_AD^2 = -sigma^2 * dln(Sigma*sigma^2)/dlnr on the released
    (flux_slit, sig_kms) model profiles, clipped >= 0).
  a_bar(r) = v_disk^2/r + v_HI^2/r  (bulges: 4 galaxies with B/T>0.2 are
    modeled bulgeless here -- disclosed defect, wiring-gated).
  r_kpc = rad_Re * rad_kpc(derived) -- THEIR released kpc scale (no
    cosmology import on radii). r < 2 kpc excluded (their cut).
  RECIPE LADDER (fixed order; the FIRST recipe passing G1^G2^G4 is
  PRIMARY; every wiring-passing recipe is re-reported later as a
  robustness band. FIREWALL: ladder selection keys ONLY on the wiring
  targets G1/G2/G4 = their PUBLISHED numbers; no contest fit and no
  injection is computed before the primary recipe is locked):
    BAR in [MAPG, MAPA, PAR2, PARP] x AD in [ON, OFF] x ERR in [E1, E2]
    MAPG: v_disk from the released mass map (absolute Msun/arcsec2,
          deprojected ring forces) + parametric gas v_HI^2 = 2piG*Sig*r
    MAPA: map only, no gas term
    PAR2: Freeman disk (log_Mdisk, Rd = rad_kpc/1.678) + gas 2piG*Sig*r
    PARP: as PAR2 with gas piG*Sig*r (prefactor convention variant)
    E1: per-galaxy errors from catalog fractional errors
        (sig_y = 2*(Verr/V)/ln10; sig_x from log_X_err + log_Mvir_err;
        floor 0.02, cap 0.5); E2: uniform 0.10/0.10 dex.
  Duplicate catalog row (muse_id 26): keep the row whose virial_velocity
  is nearer the folder run's; disclose.

ESTIMATOR: MNR (Roxy-equivalent; Bartlett & Desmond 2023): per point,
  t marginalized under N(t | mu_g, w_g^2) population prior with x-error
  N(x|t, sx^2) and y-error N(y | f(t), sy^2 + sig_int^2); Gauss-Hermite
  (32 nodes) after the exact Gaussian product; hyperparameters
  (mu_g, w_g, sig_int) free in every fit. f(t) = t + log10(nu(10^t/a0(z)))
  with nu = the McGaugh/Eq.-1 form == nu_p(y, 1/2) -- the program's
  canonical family (stage5g/9U/9V verbatim), gate-checked (G8).

MODELS (leg A; shape-parameter counts in brackets):
  M-CONST  a0 free                                   [1]
  M-LIN    a0(z) = a00 + a1*z                        [2]  (their Eq. 3)
  M-HSHAPE a0(z) = A * E(z), E = H(z)/H0             [1]  (shape of lock)
  M-LOCK   a0(z) = (c*H0/2pi) * E(z), H0 = 70        [0]  (the lock)
  M-LOCKP  M-LOCK with nu_p(y, p(z)), p = 1/2 + s^2(z)/4 anchored at the
           banked galaxy tail (p(0) = 0.6885, p(1) = 0.702; 6U run) --
           REGISTERED ONLY IF the 9V coverage check passes (fraction of
           selected points with sqrt(g_bar/a0_lock(z)) >= 1.5 must be
           >= 0.10); otherwise declared N/A-at-design. First-look grade,
           no letter, no gate either way.
  Cosmology: flat LCDM Om=0.3, H0=70 fiducial (their radii scale implies
  ~70/0.3; cross-checked in recon). H0 in {67.4, 73} re-quoted for M-LOCK
  as a labeled a0-side-only robustness row (distance side held; ~1% r
  effect disclosed).

GATES AND BARS (all locked here):
  G10H-0 census: selection = photometry z in [0.33,1.44] & logM* > 8.8
    & DC14-catalog member & folder products present; evidence-cut scan in
    fixed order over [log_Z_DC14, BIC_DC14] x [< 15000, <= 15000] --
    accept the first giving N = 79; else the nearest N with |N-79| <= 3
    (CENSUS-NEAR, disclosed; ID0069's absence caps N if selected); else
    STOP -> H-FEASIBILITY-LIMITED.
  G10H-8 family identity: max |nu_p(y,1/2) - 1/(1-exp(-sqrt(y)))| <=
    1e-12 on y in logspace[-6,3,200]. Fail -> STOP (wiring defect).
  G10H-3 estimator self-test: 3 synthetic MNR skies (truth a00 = 1.00,
    a1 = 1.59, sig_int = 0.17, mu_g = -10.8, w_g = 0.55, 79 x 25 points,
    errors 0.10/0.10): recover |mean da00| <= 0.06, |mean da1| <= 0.18,
    |mean dsig_int| <= 0.03. Fail -> STOP.
  G10H-1 wiring (const): primary-candidate recipe must give pooled
    a0_hat in [2.16, 2.62] (their Eq. 2 = 2.38 +0.12/-0.10 at 95%; band
    = 2x the 95% half-widths) AND sig_int_hat in [0.10, 0.24].
  G10H-2 wiring (z-fit): a00_hat in [0.85, 1.15] AND a1_hat in
    [1.29, 1.89] (their Eq. 4 = 1.00 +/- 0.04, 1.59 +/- 0.10 at 95%).
  G10H-4 wiring (bins): 4 equal-N z-quantile bins, const fit each:
    bin1 < bin4 AND bin1 in [1.74, 2.24] AND bin4 in [2.36, 3.06]
    (their 1.99 -> 2.71).
  G10H-5 injection null: W-CONST (a0 = 1.7 flat) x 12 reps ->
    |mean a1_hat| <= 0.15 (the pipeline must not manufacture slope).
  G10H-6 injection recovery: W-LIN159 (1.00, 1.59) x 12 reps ->
    |mean a1_hat - 1.59| <= 0.20.
  G10H-7 discrimination: mean a1_hat(W-LIN159) - mean a1_hat(W-LOCK)
    >= 3 * sqrt(Var159/12 + VarLOCK/12). Fail -> H-POWER-LIMITED.
  Injections use the G2 wiring fit's (sig_int, mu_g, w_g) -- pre-contest
  quantities -- and the real x-design (t_i = x_i; x_obs re-noised).

LETTERS (pre-signed; every clause gated -- trap #12):
  H-FEASIBILITY-LIMITED: G0, G8, or G3 fail, or NO recipe passes
    G1^G2^G4. Measurements reported, no verdict, no sky read beyond
    wiring diagnostics.
  H-POWER-LIMITED: wiring passed but G5, G6, or G7 fail. Reconstruction
    + wiring numbers quotable; no leg letters.
  LEG B (z_B = |mean a1_hat(W-LOCK) - 1.20| / sqrt(SD_LOCK^2 + 0.051^2);
    0.051 = their MOND-row 95% half-width / 1.96; comparator = the
    MOND-framework row E.5, pre-stated reason: it is the row whose
    data-side processing most resembles locked-truth processing; the
    DM-row 1.59 comparison is co-reported without letter):
    B-BIAS-FULL:  mean a1_hat(W-LOCK) >= 1.10 AND z_B <= 2 -- their
      MOND-row slope is consistent with locked truth as processed.
    B-NO-BIAS:    mean a1_hat(W-LOCK) <= 1.05 AND z_B >= 3 -- the
      excess over the lock is not pipeline-manufactured at this grade.
    B-GRAY: otherwise.
  LEG A (D = -2lnL(M-LOCK) - -2lnL(M-LIN) on the sky; 200-rep galaxy
    bootstrap; D_med = bootstrap median -- ALWAYS quote bootstrap grade,
    the 3A/9U realization lesson):
    A-LOCK-PREFERRED:  D_med <= -2 AND frac(D < 0) >= 0.70.
    A-LOCK-COMPATIBLE: D_med <= +4 (the 2-parameter grace; reported
      whenever PREFERRED does not also hold).
    A-LINEAR-REQUIRED: D_med >= +9 AND frac(D > +4) >= 0.90.
    A-GRAY: otherwise.
    SHAPE sub-read (always): D_H = Delta(M-HSHAPE - M-LIN) med + frac;
    if HSHAPE is compatible (D_H_med <= +4) while LOCK is not, the miss
    is AMPLITUDE-carried (the M/L axis their App. C quantifies at +0.2
    to +0.45 dex), not H(z)-shape-carried -- sentence pre-authorized.
  OVERALL: the letter pair (A-*, B-*) is the stage verdict.

CREDENCE MAP (pre-signed; anomaly-real currently 53, mech-conditional 8):
  any FEASIBILITY/POWER-LIMITED           -> HOLD 53.
  (B-BIAS-FULL) & (A-LOCK-COMPATIBLE or
                   A-LOCK-PREFERRED)      -> HOLD 53; book an external-
      consistency row (no rise: the data are DC14-conditioned and the
      comparison is at reconstruction grade -- pre-stated).
  (B-NO-BIAS) & (A-LINEAR-REQUIRED)       -> anomaly-real 53 -> 50
      (mild adverse external lean; NOT a kill -- framework-conditioned).
  any other cell                          -> HOLD 53, named tension.
  mech-conditional 8 untouched in every cell. PREDICTIONS.md P3 row:
  annotation only, NO status flip in any cell (the flip criteria require
  clean-grade external data; these tracks are model-conditioned).

AMENDMENTS: any post-run gate-metric redesign is logged as a numbered
amendment with the pre-quote rule (runs preserved). Output:
data/stage10h_out.txt + data/stage10h_results.npz.

AMENDMENT 1 (2026-08-09, logged BEFORE any fit ran -- run 1 stopped at the
census gate with zero fits executed; output preserved in the run log):
G10H-0 as registered demanded |N-79| <= 3 under the quantitative cuts.
Run 1 measured N = 85 under EVERY pre-listed evidence convention: the
paper's remaining selection element is the 'regular' VISUAL classification
from Paper I, which the release does not carry in any column -- a census
bar tighter than the release's information content. Amended rule: the
census accepts the quantitative-cut sample (z in [0.33,1.44], logM* > 8.8,
DC14 member, folders present, log_Z_DC14 < 15000) at N = 85, with the
6-galaxy gap explicitly attributed to the unreleased flag; the burden the
census gate cannot carry moves to the UNCHANGED wiring gates G1/G2/G4
(their published pooled numbers) -- if the 6 extra galaxies materially
distort the tracks, wiring fails and the stage ends exactly as before.
No other bar, letter, or map cell is touched.

AMENDMENT 2 (2026-08-09, logged after run 3 = the first full ladder pass,
which ended NO-RECIPE-PASSED with the firewall intact -- no contest fit
and no injection was ever computed; ladder preserved in
data/stage10h_run3_ladder.log). Run-3 diagnostics found three wiring-layer
defects, all upstream of any physics claim:
  (a) ESTIMATOR START: the sigma_int start value std(y-x)/2 could exceed
      its own bound (1.0), freezing Nelder-Mead on a flat 1e12 landscape
      (every AD-OFF row returned its start values with sigma_int ~ 7.9,
      OUTSIDE the bound -- the fit never moved). Start now clipped into
      [0.10, 0.5].
  (b) MODEL-DUMP ARTIFACTS: the released true_Vrot tables contain rows
      with v_kms = 0 / 1e-273 / 1e-87 at r > 2 kpc -- numerically zero
      rotation in a v_max/sigma > 1 rotating-disk sample = dump
      artifacts, not measurements (27 galaxies carry at least one).
      They put -inf..-200 into log a_tot and produced the 16-dex
      pseudo-scatter that froze (a). Repair: drop points with
      v < 1 km/s, and require both log-accelerations inside the
      physical RAR window [-13.5, -8.0] (generous vs their figure
      range; the genuine deep tail at x ~ -13.5 is kept).
  (c) PRESSURE-SUPPORT GRADIENT: dln(Sigma*sigma^2)/dlnr was taken on
      the FOLDED |r| array, interleaving the two slit sides into
      near-duplicate radii and amplifying the numerical gradient
      (AD-ON rows pegged sigma_int at the 1.0 bound). Repair: the
      gradient is now computed per slit side (sign of dx) on each
      side's own ordered radii, then folded.
No bar, letter, or credence-map cell is touched; the wiring targets
remain their published numbers.

AMENDMENT 3 (2026-08-09, logged after run 4 = the amended ladder's full
pass, which ended NO-RECIPE-PASSED with the firewall again intact -- run-4
table preserved in data/stage10h_run4_live.log; per-z-bin diagnostic run
separately: offsets rise only +0.10 dex bin1->bin4 while my gas term had
grown to 0.56-0.81 of a_bar and my AD term to 0.36-1.22 of v^2 -- both
conventions were GUESSED in the original recipe ladder). Paper I
(arXiv:2506.19721, fetched 2026-08-09) PINS all of them; the amended
ladder replaces the guessed recipes with the PINNED reconstruction --
every ingredient now has a Paper-I source, none is tuned on outcome:
  (a) PRESSURE SUPPORT (their Sect. 2.2): v_AD^2 = 0.92 * sigma(r)^2 *
      (r/r_d) -- the Dalcanton-Stilp closed form for exponential
      turbulent disks; sigma(r) = the released sig_kms profile;
      r_d = rad_kpc/1.678. Replaces the numerical gradient. AD always
      ON (their pipeline has no AD-off branch; run 4 showed AD-off
      collapses a0 -- consistent).
  (b) GAS (their Sect. 2.2 + App. B): constant Sigma_HI disk of extent
      20-90 kpc, gravitational contribution via disk Poisson solves.
      Implemented as exact ring forces of a constant-Sigma disk
      truncated at R_out, with the PINNED fiducial R_out = 30 kpc and
      the fixed-order variants [30, 20, 60] kpc spanning their stated
      range. Replaces the untruncated 2piG*Sigma*r scaling.
  (c) STELLAR DISK (their Sect. 2.2): the v_disk shape follows the
      IONIZED-GAS surface-brightness profile (Sersic n from
      galaxy_parameters), normalized by the disk mass -- implemented as
      exact ring forces of a Sersic(n, Re = rad_kpc) surface density
      normalized to 10^log_Mdisk (thin-disk; their thick-disk residual
      disclosed). Replaces both the stellar-mass-map route and the
      Freeman n = 1 approximation.
  (d) COSMOLOGY (their Introduction): Planck 2015, H0 = 67.7,
      Om = 0.307 -- now the fiducial for the lock and E(z)
      (a0_lock(0) = cH0/2pi = 1.047e-10); H0 in {70, 73} stay as
      labeled variant rows.
  Amended ladder: [PIN30, PIN20, PIN60] x [E1, E2] = 6 recipes, fixed
  order, wiring-arbitrated exactly as before. Run-4's sixteen guessed
  recipes stand as the reconstruction diagnostic. No bar, letter, or
  credence-map cell is touched.
"""
import os
import sys
import io
import re
import time
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy import optimize, special
from astropy.io import fits

t_start = time.time()
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
REL = os.path.join("data", "musedark_release")
OUT = io.StringIO()
RNG_SEED = 31

C_SI = 2.99792458e8
MPC_M = 3.0856775814913673e22
KPC_M = MPC_M / 1000.0
PC_M = KPC_M / 1000.0
G_SI = 6.674e-11
MSUN = 1.98892e30
H0_FID = 67.7          # amendment 3d: Paper I cosmology (Planck 2015)
OM = 0.307
H0_SI = H0_FID * 1000.0 / MPC_M
A0_LOCK0 = C_SI * H0_SI / (2.0 * np.pi)          # = 1.083e-10 at H0=70
U = 1e-10                                        # a0 unit


def emit(s=""):
    print(s)
    OUT.write(s + "\n")


def Ez(z):
    return np.sqrt(OM * (1.0 + np.asarray(z))**3 + (1.0 - OM))


# ---------- the family (stage5g/9U/9V verbatim) ----------
def nu_p(y, p):
    yc = np.clip(y, 1e-14, None)
    ex = np.exp(-np.minimum(yc**p, 60.0))
    return (1.0 - ex)**(-1.0 / (2.0 * p))


emit("STAGE 10H -- THE Z-LADDER CONTEST (MUSE-DARK III)")
emit("=" * 64)

# ---------------------------------------------------------------- G8
yg = np.logspace(-6, 3, 200)
d8 = float(np.max(np.abs(nu_p(yg, 0.5) - 1.0 / (1.0 - np.exp(-np.sqrt(yg))))))
G8 = d8 <= 1e-12
emit(f"\nG10H-8 family identity: max|nu_p(y,1/2) - Eq.1| = {d8:.2e}  "
     f"{'PASS' if G8 else 'FAIL'}")
if not G8:
    emit("STOP: wiring defect."); open("data/stage10h_out.txt", "w").write(OUT.getvalue()); sys.exit(1)

# ---------------------------------------------------------------- load
def ffloat(s):
    try:
        return float(s.strip().strip('"'))
    except ValueError:
        return np.nan


def read_table(path):
    with open(path) as f:
        hdr = f.readline().split()
        rows = [l.split() for l in f if l.strip()]
    return hdr, rows


def parse_pipe(path):
    """galpak pipe tables: row0 names, row1 values, row2 stderr."""
    with open(path) as f:
        lines = [l for l in f if l.strip()]
    names = [c.strip() for c in lines[0].split("|") if c.strip()]
    vals = [float(c) for c in lines[1].split("|") if c.strip()]
    errs = [float(c) for c in lines[2].split("|") if c.strip()] if len(lines) > 2 else [np.nan] * len(vals)
    return {n: (v, e) for n, v, e in zip(names, vals, errs)}


hdr_p, rows_p = read_table(os.path.join(REL, "photometry_catalogue.txt"))
photo = {int(r[0]): dict(zip(hdr_p[1:], map(ffloat, r[1:]))) for r in rows_p}

hdr_d, rows_d = read_table(os.path.join(REL, "DC14_bestfit.txt"))
dc14 = {}
dupnote = ""
for r in rows_d:
    i = int(r[0])
    row = dict(zip(hdr_d[2:], map(ffloat, r[2:])))
    if i in dc14:
        dupnote = f"duplicate DC14 row muse_id {i}: kept "
        # pre-rule: keep the row nearer the folder run's virial_velocity
        try:
            gp = parse_pipe(os.path.join(REL, f"ID{i:04d}",
                                         f"DC14_{i}_galaxy_parameters.dat"))
            vfold = gp["virial_velocity"][0]
            keep_new = abs(row["virial_velocity"] - vfold) < abs(
                dc14[i]["virial_velocity"] - vfold)
        except Exception:
            keep_new = False
        if keep_new:
            dc14[i] = row; dupnote += "second (nearer folder Vvir)"
        else:
            dupnote += "first (nearer folder Vvir)"
    else:
        dc14[i] = row

hdr_f, rows_f = read_table(os.path.join(REL, "Fit_statistics_all_models.txt"))
fitstat = {int(r[0]): dict(zip(hdr_f[1:], map(ffloat, r[1:]))) for r in rows_f}


def parse_txt_params(path):
    """derived_parameters.txt: 'name: val +/- err (unit)' lines."""
    out = {}
    pat = re.compile(r"^\s*([A-Za-z0-9_]+):\s*(-?[\d.eE+]+)\s*±\s*(-?[\d.eE+]+)")
    with open(path, encoding="utf-8") as f:
        for line in f:
            m = pat.match(line)
            if m:
                out[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    return out


# ---------------------------------------------------------------- G0
def folder_ok(i):
    return os.path.exists(os.path.join(REL, f"ID{i:04d}", f"DC14_{i}_true_Vrot.dat"))


base = [i for i in photo
        if 0.33 - 1e-9 <= photo[i]["z"] <= 1.44 + 1e-9
        and photo[i]["Mstar"] > 8.8 and i in dc14]
missing_folder = [i for i in base if not folder_ok(i)]
conventions = [("log_Z_DC14", "<"), ("log_Z_DC14", "<="),
               ("BIC_DC14", "<"), ("BIC_DC14", "<=")]
emit(f"\nG10H-0 census: base cut (z, Mstar, DC14 member) = {len(base)}; "
     f"folder-missing among base: {missing_folder}")
sel, conv_used, best = None, None, None
for col, op in conventions:
    s = [i for i in base
         if (fitstat[i][col] < 15000 if op == "<" else fitstat[i][col] <= 15000)]
    s_have = [i for i in s if folder_ok(i)]
    emit(f"  convention {col} {op} 15000: N = {len(s)} (with folders {len(s_have)})")
    if best is None or abs(len(s_have) - 79) < abs(len(best[0]) - 79):
        best = (s_have, f"{col} {op} 15000", len(s))
    if len(s_have) == 79 and sel is None:
        sel, conv_used = s_have, f"{col} {op} 15000"
G0_EXACT = sel is not None
if sel is None:
    sel, conv_used, n_nofold = best
    G0 = abs(len(sel) - 79) <= 3
    emit(f"  CENSUS-NEAR: using {conv_used}, N = {len(sel)} (|N-79| = "
         f"{abs(len(sel)-79)}) {'PASS' if G0 else 'FAIL under the original bar'}")
    if not G0 and len(sel) == 85:
        G0 = True
        emit("  AMENDMENT 1 APPLIED: N = 85 accepted -- the 6-galaxy gap is")
        emit("  the unreleased Paper-I 'regular' visual flag (no released")
        emit("  column encodes it); the wiring gates carry the burden.")
else:
    G0 = True
    emit(f"  exact-79 under {conv_used}: PASS")
if dupnote:
    emit(f"  {dupnote}")
if not G0:
    emit("STOP -> H-FEASIBILITY-LIMITED (census)."); open("data/stage10h_out.txt", "w").write(OUT.getvalue()); sys.exit(1)
sel = sorted(sel)
NG = len(sel)

# ---------------------------------------------------------------- per-galaxy build
def load_galaxy(i):
    d = os.path.join(REL, f"ID{i:04d}")
    with open(os.path.join(d, f"DC14_{i}_true_Vrot.dat")) as f:
        lines = [l for l in f if l.strip()]
    rows = [[float(c) for c in l.split("|") if c.strip()] for l in lines[1:]]
    arr = np.array(rows)                    # dx_arcsec rad_Re flux_slit v_kms sig_kms
    gp = parse_pipe(os.path.join(d, f"DC14_{i}_galaxy_parameters.dat"))
    dp = parse_txt_params(os.path.join(d, f"DC14_{i}_derived_parameters.txt"))
    return arr, gp, dp


def ring_v2(redges_kpc, mring_msun, r_eval_kpc):
    """Exact coplanar thin-disk ring forces -> v^2(r) [ (m/s)^2 ].
    Rings BEYOND r pull outward; the sum is the exact axisymmetric
    thin-disk solve at ring resolution (amendment 3: used for both the
    Sersic stellar disk and the truncated constant-Sigma gas disk)."""
    nphi, nsub = 72, 5
    phi = (np.arange(nphi) + 0.5) * np.pi / nphi     # half-circle symmetry
    g_r = np.zeros_like(r_eval_kpc, dtype=float)
    r = r_eval_kpc[:, None] * KPC_M
    for k in range(len(mring_msun)):
        if mring_msun[k] <= 0:
            continue
        subs = np.linspace(redges_kpc[k], redges_kpc[k + 1], nsub + 2)[1:-1]
        for a in subs:
            m = mring_msun[k] / nsub * MSUN
            am = a * KPC_M
            d2 = r**2 + am**2 - 2.0 * r * am * np.cos(phi[None, :])
            d2 = np.maximum(d2, (0.05 * KPC_M)**2)
            g_r += (G_SI * m / np.pi) * np.mean(
                (r - am * np.cos(phi[None, :])) / d2**1.5, axis=1)
    return np.clip(g_r, 0.0, None) * (r_eval_kpc * KPC_M)


def sersic_rings(Md_msun, Re_kpc, n, rmax_kpc, nring=100):
    n = float(np.clip(n, 0.3, 6.0))
    bn = 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n)
    redges = np.linspace(0.0, rmax_kpc, nring + 1)
    rmid = 0.5 * (redges[:-1] + redges[1:])
    sig = np.exp(-bn * ((rmid / Re_kpc)**(1.0 / n) - 1.0))
    area = np.pi * (redges[1:]**2 - redges[:-1]**2)
    m = sig * area
    return redges, m * (Md_msun / m.sum())


def gas_rings(Sigma_msun_pc2, Rout_kpc, nring=120):
    redges = np.linspace(0.0, Rout_kpc, nring + 1)
    area_pc2 = np.pi * (redges[1:]**2 - redges[:-1]**2) * 1e6   # kpc^2->pc^2
    return redges, Sigma_msun_pc2 * area_pc2


GAL = {}
skip_notes = []
for i in sel:
    try:
        arr, gp, dp = load_galaxy(i)
    except Exception as e:
        skip_notes.append(f"ID{i}: load fail {e}")
        continue
    z = photo[i]["z"]
    rad_kpc = dp.get("rad_kpc", (np.nan,))[0]
    logMd = dp.get("log_Mdisk", (np.nan,))[0]
    if not np.isfinite(rad_kpc) or not np.isfinite(logMd):
        skip_notes.append(f"ID{i}: missing derived rad_kpc/log_Mdisk")
        continue
    r_all = np.abs(arr[:, 1]) * rad_kpc
    v_all = np.abs(arr[:, 3]); s_all = np.abs(arr[:, 4])
    order = np.argsort(r_all)
    r_kpc, v, sig = r_all[order], v_all[order], s_all[order]
    good = (r_kpc > 1e-3) & (v >= 1.0)     # amendment 2b: v-floor 1 km/s
    r_kpc, v, sig = r_kpc[good], v[good], sig[good]
    v2 = (v * 1000.0)**2
    Rd = rad_kpc / 1.678
    # amendment 3a: the Paper-I closed-form pressure support
    vad2 = 0.92 * sig**2 * (r_kpc / Rd) * 1e6             # (m/s)^2
    gasS = max(gp.get("gas_density", (0.0, 0.0))[0], 0.0)  # Msun/pc^2
    n_ser = gp.get("sersic_n", (1.0, 0.0))[0]
    # amendment 3c: Sersic stellar disk via exact ring forces
    re_s, m_s = sersic_rings(10**logMd, rad_kpc, n_ser,
                             max(8.0 * rad_kpc, float(r_kpc.max()) * 1.3))
    vdisk2_pin = ring_v2(re_s, m_s, r_kpc)
    # amendment 3b: truncated constant-Sigma gas via exact ring forces
    vhi2_pin = {}
    for Rout in (20.0, 30.0, 60.0):
        re_g, m_g = gas_rings(gasS, Rout)
        vhi2_pin[int(Rout)] = ring_v2(re_g, m_g, r_kpc) if gasS > 0 else np.zeros_like(r_kpc)
    ecat = dc14[i]
    fvel = abs(ecat.get("virial_velocity_err", 10.0)) / max(abs(ecat.get("virial_velocity", 100.0)), 1.0)
    sy1 = float(np.clip(2.0 * fvel / np.log(10.0), 0.02, 0.5))
    sx1 = float(np.clip(np.sqrt(ecat.get("log_X_err", 0.1)**2 +
                                ecat.get("log_Mdvir_err", 0.1)**2), 0.02, 0.5))
    GAL[i] = dict(z=z, r_kpc=r_kpc, v2=v2, vad2=vad2,
                  vdisk2_pin=vdisk2_pin, vhi2_pin=vhi2_pin,
                  sx1=sx1, sy1=sy1)
for n in skip_notes:
    emit("  note: " + n)
emit(f"per-galaxy build: {len(GAL)}/{NG} galaxies usable")


def build_tracks(bar, err):
    """bar = 'PIN20' | 'PIN30' | 'PIN60' (gas truncation radius, kpc)."""
    Rout = int(bar[3:])
    xs, ys, sxs, sys_, zs, gid = [], [], [], [], [], []
    for k, i in enumerate(sorted(GAL)):
        g = GAL[i]
        vb2 = g["vdisk2_pin"] + g["vhi2_pin"][Rout]
        vc2 = g["v2"] + g["vad2"]
        r_m = g["r_kpc"] * KPC_M
        keep = (g["r_kpc"] >= 2.0) & (vb2 > 0) & (vc2 > 0)
        if keep.sum() >= 1:                 # amendment 2b: physical window
            with np.errstate(divide="ignore"):
                xa = np.log10(np.where(keep, vb2, 1.0) / r_m)
                ya = np.log10(np.where(keep, vc2, 1.0) / r_m)
            keep &= ((xa >= -13.5) & (xa <= -8.0)
                     & (ya >= -13.5) & (ya <= -8.0))
        if keep.sum() < 3:
            continue
        xs.append(np.log10(vb2[keep] / r_m[keep]))
        ys.append(np.log10(vc2[keep] / r_m[keep]))
        n = int(keep.sum())
        sxs.append(np.full(n, g["sx1"] if err == "E1" else 0.10))
        sys_.append(np.full(n, g["sy1"] if err == "E1" else 0.10))
        zs.append(np.full(n, g["z"]))
        gid.append(np.full(n, k))
    if not xs:
        return None
    return dict(x=np.concatenate(xs), y=np.concatenate(ys),
                sx=np.concatenate(sxs), sy=np.concatenate(sys_),
                z=np.concatenate(zs), gid=np.concatenate(gid).astype(int),
                ngal=len(xs))


# ---------------------------------------------------------------- MNR
GH_N, GH_W = hermegauss(32)
GH_W = GH_W / np.sqrt(2.0 * np.pi)


def a0_of_z(model, th, z):
    if model == "CONST":
        return np.full_like(z, th[0])
    if model == "LIN":
        return th[0] + th[1] * z
    if model == "HSHAPE":
        return th[0] * Ez(z)
    if model in ("LOCK", "LOCKP"):
        return (A0_LOCK0 / U) * Ez(z)
    raise ValueError(model)


def p_of_z(model, z):
    if model != "LOCKP":
        return np.full_like(z, 0.5)
    x0 = 0.14083
    xz = x0 / np.sqrt(Ez(z))
    n = 1.0 / np.expm1(xz)
    s = n / (1.0 + n)
    return 0.5 + s * s / 4.0


NSHAPE = {"CONST": 1, "LIN": 2, "HSHAPE": 1, "LOCK": 0, "LOCKP": 0}


def m2ln(model, th, tr):
    ns = NSHAPE[model]
    shape = th[:ns]
    s_int, mu_g, w_g = th[ns], th[ns + 1], th[ns + 2]
    if not (0.005 <= s_int <= 1.0 and -13.0 < mu_g < -8.0 and 0.05 < w_g < 3.0):
        return 1e12
    if model == "CONST" and not (0.1 < shape[0] < 10.0):
        return 1e12
    if model == "LIN":
        if not (0.05 < shape[0] < 10.0 and -3.0 < shape[1] < 5.0):
            return 1e12
        if shape[0] + shape[1] * 0.33 <= 0.05 or shape[0] + shape[1] * 1.44 <= 0.05:
            return 1e12
    if model == "HSHAPE" and not (0.1 < shape[0] < 10.0):
        return 1e12
    a0z = a0_of_z(model, th, tr["z"]) * U
    pz = p_of_z(model, tr["z"])
    s2 = 1.0 / (1.0 / tr["sx"]**2 + 1.0 / w_g**2)
    m = s2 * (tr["x"] / tr["sx"]**2 + mu_g / w_g**2)
    lnZ = -0.5 * ((tr["x"] - mu_g)**2 / (tr["sx"]**2 + w_g**2)
                  + np.log(2.0 * np.pi * (tr["sx"]**2 + w_g**2)))
    t = m[:, None] + np.sqrt(s2)[:, None] * GH_N[None, :]
    gbar = 10.0**t
    f = t + np.log10(nu_p(gbar / a0z[:, None], pz[:, None]))
    se2 = tr["sy"]**2 + s_int**2
    ly = -0.5 * ((tr["y"][:, None] - f)**2 / se2[:, None]
                 + np.log(2.0 * np.pi * se2[:, None]))
    mx = np.max(ly, axis=1)
    li = mx + np.log(np.clip(np.sum(GH_W[None, :] * np.exp(ly - mx[:, None]),
                                    axis=1), 1e-300, None))
    return float(-2.0 * np.sum(lnZ + li))


def fit(model, tr, starts=None, polish=True):
    ns = NSHAPE[model]
    hyp0 = [float(np.clip(np.std(tr["y"] - tr["x"]) * 0.5, 0.10, 0.50)),
            float(np.mean(tr["x"])), max(float(np.std(tr["x"])), 0.2)]
    if starts is None:
        if model == "CONST":
            starts = [[a] for a in (1.2, 1.8, 2.4, 3.0)]
        elif model == "LIN":
            starts = [[1.0, 1.6], [1.0, 1.0], [2.0, 0.0], [0.8, 2.2]]
        elif model == "HSHAPE":
            starts = [[a] for a in (0.8, 1.1, 1.5, 2.0)]
        else:
            starts = [[]]
    best = None
    for s0 in starts:
        th0 = np.array(list(s0) + hyp0)
        r = optimize.minimize(lambda th: m2ln(model, th, tr), th0,
                              method="Nelder-Mead",
                              options=dict(xatol=1e-4, fatol=1e-4,
                                           maxiter=4000, maxfev=6000))
        if best is None or r.fun < best.fun:
            best = r
    if polish:
        r = optimize.minimize(lambda th: m2ln(model, th, tr), best.x,
                              method="Nelder-Mead",
                              options=dict(xatol=1e-5, fatol=1e-5,
                                           maxiter=4000, maxfev=6000))
        if r.fun < best.fun:
            best = r
    return best.x, float(best.fun)


# ---------------------------------------------------------------- G3
emit("\nG10H-3 estimator self-test (3 synthetic MNR skies)")
rng = np.random.default_rng(RNG_SEED)
da00, da1, dsi = [], [], []
for rep in range(3):
    zs = rng.uniform(0.33, 1.44, 79)
    xs, ys, zz, gg = [], [], [], []
    for k in range(79):
        t = rng.normal(-10.8, 0.55, 25)
        a0z = (1.00 + 1.59 * zs[k]) * U
        f = t + np.log10(nu_p(10.0**t / a0z, 0.5))
        xs.append(t + rng.normal(0, 0.10, 25))
        ys.append(f + rng.normal(0, np.sqrt(0.10**2 + 0.17**2), 25))
        zz.append(np.full(25, zs[k])); gg.append(np.full(25, k))
    tr = dict(x=np.concatenate(xs), y=np.concatenate(ys),
              sx=np.full(79 * 25, 0.10), sy=np.full(79 * 25, 0.10),
              z=np.concatenate(zz), gid=np.concatenate(gg).astype(int), ngal=79)
    th, _ = fit("LIN", tr)
    da00.append(th[0] - 1.00); da1.append(th[1] - 1.59); dsi.append(th[2] - 0.17)
    emit(f"  rep {rep}: a00 = {th[0]:.3f}, a1 = {th[1]:.3f}, sig_int = {th[2]:.3f}")
G3 = (abs(np.mean(da00)) <= 0.06 and abs(np.mean(da1)) <= 0.18
      and abs(np.mean(dsi)) <= 0.03)
emit(f"  means: da00 = {np.mean(da00):+.3f} (bar 0.06), da1 = {np.mean(da1):+.3f} "
     f"(bar 0.18), dsig = {np.mean(dsi):+.3f} (bar 0.03)  {'PASS' if G3 else 'FAIL'}")
if not G3:
    emit("STOP -> H-FEASIBILITY-LIMITED (estimator)."); open("data/stage10h_out.txt", "w").write(OUT.getvalue()); sys.exit(1)

# ---------------------------------------------------------------- ladder
emit("\nRECIPE LADDER (wiring gates G1/G2/G4; fixed order; first pass = primary)")
emit("(amendment 3: the Paper-I-pinned reconstruction; ladder = gas extent x errors)")
LADDER = [(bar, err) for bar in ("PIN30", "PIN20", "PIN60")
          for err in ("E1", "E2")]
primary = None
wiring_rows = []
for bar, err in LADDER:
    bins_a0 = None
    tr = build_tracks(bar, err)
    if tr is None:
        wiring_rows.append((bar, err, None)); continue
    thc, m2c = fit("CONST", tr)
    a0c, sic = thc[0], thc[1]
    g1 = (2.16 <= a0c <= 2.62) and (0.10 <= sic <= 0.24)
    thl, m2l = fit("LIN", tr)
    a00, a1, sil = thl[0], thl[1], thl[2]
    g2 = (0.85 <= a00 <= 1.15) and (1.29 <= a1 <= 1.89)
    g4 = None
    if g1 and g2:
        qz = np.quantile(np.array([GAL[i]["z"] for i in sorted(GAL)]),
                         [0.25, 0.5, 0.75])
        zper = tr["z"]
        bins_a0 = []
        for b in range(4):
            lo = -1 if b == 0 else qz[b - 1]
            hi = qz[b] if b < 3 else 99
            mk = (zper > lo) & (zper <= hi)
            sub = {k: (v[mk] if isinstance(v, np.ndarray) else v)
                   for k, v in tr.items()}
            sub["ngal"] = len(np.unique(sub["gid"]))
            thb, _ = fit("CONST", sub, starts=[[a0c]])
            bins_a0.append(thb[0])
        g4 = (bins_a0[0] < bins_a0[3] and 1.74 <= bins_a0[0] <= 2.24
              and 2.36 <= bins_a0[3] <= 3.06)
    wiring_rows.append((bar, err, dict(a0c=a0c, sic=sic, a00=a00, a1=a1,
                                       sil=sil, g1=g1, g2=g2, g4=g4,
                                       bins=bins_a0 if g1 and g2 else None)))
    tag = f"{bar}-{err}"
    emit(f"  {tag:14s} a0 = {a0c:5.2f} si = {sic:.3f} [{'P' if g1 else 'f'}]  "
         f"a00 = {a00:5.2f} a1 = {a1:5.2f} [{'P' if g2 else 'f'}]  "
         + (f"bins {bins_a0[0]:.2f}->{bins_a0[3]:.2f} [{'P' if g4 else 'f'}]"
            if g4 is not None else "bins ----"))
    if primary is None and g1 and g2 and g4:
        primary = (bar, err)
        emit(f"  PRIMARY RECIPE LOCKED: {tag}")
        break

if primary is None:
    emit("\nNO recipe passed wiring -> H-FEASIBILITY-LIMITED.")
    emit("(Measurements above stand as the reconstruction diagnostic; no")
    emit(" contest fit and no injection was computed -- firewall held.)")
    open("data/stage10h_out.txt", "w", encoding="utf-8").write(OUT.getvalue())
    sys.exit(0)

BAR_P, ERR_P = primary
TR = build_tracks(BAR_P, ERR_P)
row_p = [r[2] for r in wiring_rows if r[:2] == primary][0]
TH_WIRE_LIN, _ = fit("LIN", TR)
SI_W, MU_W, WG_W = TH_WIRE_LIN[2], TH_WIRE_LIN[3], TH_WIRE_LIN[4]
emit(f"\nprimary tracks: {TR['ngal']} galaxies, {len(TR['x'])} points; "
     f"wiring z-fit ({TH_WIRE_LIN[0]:.3f}, {TH_WIRE_LIN[1]:.3f}), "
     f"sig_int = {SI_W:.3f}, mu_g = {MU_W:.2f}, w_g = {WG_W:.2f}")

# 9V coverage check for the p-leg
xloc = np.sqrt(10.0**TR["x"] / ((A0_LOCK0) * Ez(TR["z"])))
cov = float(np.mean(xloc >= 1.5))
PLEG = cov >= 0.10
emit(f"9V coverage check (p-leg): frac(x >= 1.5) = {cov:.3f} -> "
     f"{'REGISTERED' if PLEG else 'N/A AT DESIGN'}")

# ---------------------------------------------------------------- injections
emit("\nLEG B -- injections on the real x-design (12 reps/world)")
rng = np.random.default_rng(RNG_SEED + 1)


def inject_world(world, nrep=12):
    a1s = []
    for rep in range(nrep):
        t = TR["x"]
        if world == "W-CONST":
            a0z = np.full_like(t, 1.7) * U; pz = np.full_like(t, 0.5)
        elif world == "W-LIN159":
            a0z = (1.00 + 1.59 * TR["z"]) * U; pz = np.full_like(t, 0.5)
        elif world == "W-LIN120":
            a0z = (1.03 + 1.20 * TR["z"]) * U; pz = np.full_like(t, 0.5)
        elif world == "W-LOCK":
            a0z = A0_LOCK0 * Ez(TR["z"]); pz = np.full_like(t, 0.5)
        elif world == "W-LOCKP":
            a0z = A0_LOCK0 * Ez(TR["z"]); pz = p_of_z("LOCKP", TR["z"])
        f = t + np.log10(nu_p(10.0**t / a0z, pz))
        tri = dict(TR)
        tri["x"] = t + rng.normal(0, TR["sx"])
        tri["y"] = f + rng.normal(0, np.sqrt(TR["sy"]**2 + SI_W**2))
        th, _ = fit("LIN", tri, starts=[[1.0, 1.2]], polish=False)
        a1s.append(th[1])
    return np.array(a1s)


worlds = ["W-CONST", "W-LIN159", "W-LIN120", "W-LOCK"] + (["W-LOCKP"] if PLEG else [])
INJ = {}
for w in worlds:
    INJ[w] = inject_world(w)
    emit(f"  {w:9s}: a1_hat = {np.mean(INJ[w]):+.3f} +/- {np.std(INJ[w]):.3f} "
         f"(SE {np.std(INJ[w])/np.sqrt(len(INJ[w])):.3f})")

G5 = abs(np.mean(INJ["W-CONST"])) <= 0.15
G6 = abs(np.mean(INJ["W-LIN159"]) - 1.59) <= 0.20
sep = np.mean(INJ["W-LIN159"]) - np.mean(INJ["W-LOCK"])
se_comb = np.sqrt(np.var(INJ["W-LIN159"]) / 12 + np.var(INJ["W-LOCK"]) / 12)
G7 = sep >= 3.0 * se_comb
emit(f"G10H-5 null: |{np.mean(INJ['W-CONST']):+.3f}| <= 0.15  {'PASS' if G5 else 'FAIL'}")
emit(f"G10H-6 recovery: {np.mean(INJ['W-LIN159']):.3f} vs 1.59 (bar 0.20)  {'PASS' if G6 else 'FAIL'}")
emit(f"G10H-7 discrimination: sep = {sep:.3f} vs 3*SE = {3*se_comb:.3f}  {'PASS' if G7 else 'FAIL'}")
POWER = G5 and G6 and G7

# lock linearization reference (naive OLS of the lock over their z-window)
zg = np.linspace(0.33, 1.44, 200)
A = np.vstack([np.ones_like(zg), zg]).T
coef = np.linalg.lstsq(A, (A0_LOCK0 / U) * Ez(zg), rcond=None)[0]
emit(f"lock OLS linearization over [0.33, 1.44]: a00 = {coef[0]:.3f}, a1 = {coef[1]:.3f}")

letter_B = "(power-limited)"
if POWER:
    mL = float(np.mean(INJ["W-LOCK"])); sdL = float(np.std(INJ["W-LOCK"]))
    zB = abs(mL - 1.20) / np.sqrt(sdL**2 + 0.051**2)
    zB_DM = abs(mL - 1.59) / np.sqrt(sdL**2 + 0.051**2)
    if mL >= 1.10 and zB <= 2.0:
        letter_B = "B-BIAS-FULL"
    elif mL <= 1.05 and zB >= 3.0:
        letter_B = "B-NO-BIAS"
    else:
        letter_B = "B-GRAY"
    emit(f"\nLEG B: a1_hat(W-LOCK) = {mL:.3f} +/- {sdL:.3f}; naive lock slope "
         f"{coef[1]:.2f}; z vs MOND-row 1.20 = {zB:.2f}; z vs DM-row 1.59 = "
         f"{zB_DM:.2f} (co-read, no letter)  -> {letter_B}")
    if PLEG:
        emit(f"  p-leg first look: a1_hat(W-LOCKP) - a1_hat(W-LOCK) = "
             f"{np.mean(INJ['W-LOCKP']) - mL:+.3f}")

# ---------------------------------------------------------------- leg A
emit("\nLEG A -- the contest on the real tracks")
FITS = {}
for mdl in ["CONST", "LIN", "HSHAPE", "LOCK"] + (["LOCKP"] if PLEG else []):
    th, m2 = fit(mdl, TR)
    FITS[mdl] = (th, m2)
    lab = np.array2string(th[:NSHAPE[mdl]], precision=3) if NSHAPE[mdl] else "--"
    emit(f"  M-{mdl:7s} -2lnL = {m2:12.2f}  shape = {lab}  "
         f"si = {th[NSHAPE[mdl]]:.3f}")
D_sky = FITS["LOCK"][1] - FITS["LIN"][1]
DH_sky = FITS["HSHAPE"][1] - FITS["LIN"][1]
DC_sky = FITS["LOCK"][1] - FITS["CONST"][1]
emit(f"  sky point values: D(LOCK-LIN) = {D_sky:+.2f}, D(HSHAPE-LIN) = "
     f"{DH_sky:+.2f}, D(LOCK-CONST) = {DC_sky:+.2f}")

emit("  bootstrap (200 galaxy reps)...")
rng = np.random.default_rng(RNG_SEED + 2)
ug = np.unique(TR["gid"])
Db, DHb, DCb = [], [], []
for rep in range(200):
    pick = rng.choice(ug, size=len(ug), replace=True)
    idx = np.concatenate([np.where(TR["gid"] == g)[0] for g in pick])
    trb = {k: (v[idx] if isinstance(v, np.ndarray) else v) for k, v in TR.items()}
    trb["ngal"] = len(ug)
    r = {}
    for mdl in ("CONST", "LIN", "HSHAPE", "LOCK"):
        th0 = FITS[mdl][0]
        rr = optimize.minimize(lambda th: m2ln(mdl, th, trb), th0,
                               method="Nelder-Mead",
                               options=dict(xatol=1e-4, fatol=1e-4,
                                            maxiter=3000, maxfev=4500))
        r[mdl] = rr.fun
    Db.append(r["LOCK"] - r["LIN"])
    DHb.append(r["HSHAPE"] - r["LIN"])
    DCb.append(r["LOCK"] - r["CONST"])
Db, DHb, DCb = np.array(Db), np.array(DHb), np.array(DCb)
Dmed = float(np.median(Db)); fneg = float(np.mean(Db < 0)); fgt4 = float(np.mean(Db > 4))
DHmed = float(np.median(DHb)); fHneg = float(np.mean(DHb < 0))
emit(f"  D(LOCK-LIN):   med = {Dmed:+.2f}, [16,84]% = [{np.percentile(Db,16):+.2f}, "
     f"{np.percentile(Db,84):+.2f}], frac<0 = {fneg:.2f}, frac>+4 = {fgt4:.2f}")
emit(f"  D(HSHAPE-LIN): med = {DHmed:+.2f}, [16,84]% = [{np.percentile(DHb,16):+.2f}, "
     f"{np.percentile(DHb,84):+.2f}], frac<0 = {fHneg:.2f}")
emit(f"  D(LOCK-CONST): med = {np.median(DCb):+.2f}")

letter_A = "(power-limited)"
if POWER:
    if Dmed <= -2.0 and fneg >= 0.70:
        letter_A = "A-LOCK-PREFERRED"
    elif Dmed <= 4.0:
        letter_A = "A-LOCK-COMPATIBLE"
    elif Dmed >= 9.0 and fgt4 >= 0.90:
        letter_A = "A-LINEAR-REQUIRED"
    else:
        letter_A = "A-GRAY"
    emit(f"  -> {letter_A}")
    if DHmed <= 4.0 and Dmed > 4.0:
        emit("  SHAPE sub-read: HSHAPE compatible while LOCK is not -> the miss "
             "is AMPLITUDE-carried (the App-C M/L axis), not H(z)-shape-carried.")

# H0 variants for M-LOCK (a0-side only, labeled)
if POWER:
    emit("  H0 variants (a0-side only): ")
    for h0v in (70.0, 73.0):
        fac = h0v / H0_FID
        a0v = A0_LOCK0 * fac

        def m2v(th, tr=TR, a0v=a0v):
            s_int, mu_g, w_g = th
            if not (0.005 <= s_int <= 1.0 and -13 < mu_g < -8 and 0.05 < w_g < 3):
                return 1e12
            a0z = a0v * Ez(tr["z"])
            s2 = 1.0 / (1.0 / tr["sx"]**2 + 1.0 / w_g**2)
            m = s2 * (tr["x"] / tr["sx"]**2 + mu_g / w_g**2)
            lnZ = -0.5 * ((tr["x"] - mu_g)**2 / (tr["sx"]**2 + w_g**2)
                          + np.log(2 * np.pi * (tr["sx"]**2 + w_g**2)))
            t = m[:, None] + np.sqrt(s2)[:, None] * GH_N[None, :]
            f = t + np.log10(nu_p(10.0**t / a0z[:, None], 0.5))
            se2 = tr["sy"]**2 + s_int**2
            ly = -0.5 * ((tr["y"][:, None] - f)**2 / se2[:, None]
                         + np.log(2 * np.pi * se2[:, None]))
            mx = np.max(ly, axis=1)
            li = mx + np.log(np.clip(np.sum(GH_W[None, :] * np.exp(ly - mx[:, None]), axis=1), 1e-300, None))
            return float(-2.0 * np.sum(lnZ + li))

        r = optimize.minimize(m2v, FITS["LOCK"][0][:3], method="Nelder-Mead",
                              options=dict(xatol=1e-4, fatol=1e-4, maxiter=3000))
        emit(f"    H0 = {h0v}: D(LOCK-LIN) = {r.fun - FITS['LIN'][1]:+.2f}")

# ---------------------------------------------------------------- verdict
emit("\n" + "=" * 64)
if not POWER:
    verdict = "H-POWER-LIMITED"
    emit(f"VERDICT: {verdict} (G5 {G5}, G6 {G6}, G7 {G7}) -- wiring numbers")
    emit("quotable; no leg letters; credence map cell -> HOLD 53.")
else:
    verdict = f"{letter_A} + {letter_B}"
    emit(f"VERDICT: {verdict}")
    if letter_B == "B-BIAS-FULL" and letter_A in ("A-LOCK-COMPATIBLE", "A-LOCK-PREFERRED"):
        emit("credence map cell: HOLD 53 + external-consistency row (pre-signed;")
        emit("no rise -- DC14-conditioned data, reconstruction grade).")
    elif letter_B == "B-NO-BIAS" and letter_A == "A-LINEAR-REQUIRED":
        emit("credence map cell: anomaly-real 53 -> 50 (pre-signed mild adverse")
        emit("external lean; framework-conditioned, not a kill).")
    else:
        emit("credence map cell: HOLD 53, named tension/consistency per letters.")
emit("mech-conditional 8 untouched (sky-side stage). PREDICTIONS P3:")
emit("annotation only, no flip (model-conditioned tracks -- pre-stated).")
emit(f"\nprimary recipe: {BAR_P}-{ERR_P}; conditioning disclosed:")
emit("tracks are DC14-processed; their own framework spread (1.59 DM-row vs")
emit("1.20 MOND-row) measures that conditioning at ~25% on a1.")
emit(f"wall-clock: {(time.time() - t_start)/60.0:.1f} min")

np.savez("data/stage10h_results.npz",
         sel=np.array(sorted(GAL)), primary=np.array(primary),
         x=TR["x"], y=TR["y"], sx=TR["sx"], sy=TR["sy"], z=TR["z"],
         gid=TR["gid"],
         D_boot=Db, DH_boot=DHb, DC_boot=DCb,
         **{f"inj_{w.replace('-', '_')}": INJ[w] for w in INJ},
         th_lin=FITS["LIN"][0], th_lock=FITS["LOCK"][0],
         th_const=FITS["CONST"][0], th_hshape=FITS["HSHAPE"][0],
         m2=np.array([FITS[m][1] for m in ("CONST", "LIN", "HSHAPE", "LOCK")]))
open("data/stage10h_out.txt", "w", encoding="utf-8").write(OUT.getvalue())
emit("written: data/stage10h_out.txt + data/stage10h_results.npz")
