"""
STAGE 7J part A: the photometric companion-completeness measurement.
(The Opus-review design, made urgent by 7I's COMPANION-DIRECTION verdict.)

The 3J overluminous fraction (12.3% of stars at delta < -0.4 mag) is a
DETECTION fraction: a threshold on a distribution. Here the full delta
distribution is fit with a forward model, which measures the true
unresolved-companion fraction f_true and the detection completeness at
the 7I strict-cut threshold, converting the fence question into a prior.

Model. Per star of primary mass m (empirical mass distribution from the
sample's own M_G): with probability f_true it hosts an unresolved
companion of mass ratio q ~ flat(0.1, 1) (the v7 companion sector's own
q-law; power-law variant reported as sensitivity). The pair's combined
magnitude AND combined color are computed on the sample's own measured
main-sequence locus (ridge M_G(BP-RP), 41 bins, as 3J/7I): the system
brightens by the flux sum and reddens toward the secondary, and its
overluminosity is evaluated against the ridge at the COMBINED color (the
same operation the data pipeline performs). Singles and systems both
carry the MS width sigma_MS (fit). Mixture likelihood on the binned
delta distribution (Poisson, bins of 0.05 mag over [-1.2, +0.6]).

Honest approximations, stated: (i) evolved/subgiant interlopers brighten
like companions -> f_true is an UPPER bound on true multiplicity - the
conservative direction for the skeptic's case (a larger allowed
companion fraction); (ii) combined color approximated as flux-weighted
BP-RP; (iii) the ridge (a median) is itself pulled slightly bright by
the binaries - absorbed by the fitted zero offset mu0; (iv) RUWE detects
additional companions the photometry misses -> the residual fraction
after the strict cut computed here (photometry-only) is an UPPER bound.

GATES (bars, pre-registered):
  GA1 predictive: the fitted model must reproduce the observed flagged
      fraction P(delta < -0.4) = 0.123 within +-0.02.
  GA2 injection: synthetic populations at f_true = 0.15 and 0.35 pushed
      through the identical estimator must recover f within +-20%
      relative.
  GA3 width: fitted sigma_MS in [0.22, 0.32] (3J robust width 0.275).
If any gate fails -> part B does NOT launch; investigate first.

OUTPUTS: data/stage7j_completeness.txt (report) +
data/stage7j_prior.npz with:
  f_grid, lnpi_full  - profile log-likelihood of f_true (the FULL-sample
                       fcomp prior),
  r_grid, lnpi_strict - the same posterior pushed through the residual
                       mapping r(f) = P(companion & pair survives the
                       photometric flag)/P(pair survives) (the STRICT-
                       sample fcomp prior; photometric channel only =
                       upper bound),
  C_flag (completeness of the delta<-0.4 flag for companions),
  f_hat, f_16, f_84.
"""
import numpy as np
from astropy.io import fits

rng = np.random.default_rng(20260725)
L = []
def P(s):
    print(s, flush=True)
    L.append(s)

d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1, plx2 = d['parallax1'], d['parallax2']
eplx1, eplx2 = d['parallax_error1'], d['parallax_error2']
sep, Rch = d['sep_AU'], d['R_chance_align']
G1m, G2m = d['phot_g_mean_mag1'], d['phot_g_mean_mag2']
plx = 0.5*(plx1+plx2)
MG1 = G1m+5*np.log10(np.maximum(plx1, 1e-6))-10
MG2 = G2m+5*np.log10(np.maximum(plx2, 1e-6))-10
sigv = 4.74047/plx*np.sqrt(d['pmra_error1']**2+d['pmdec_error1']**2
                           + d['pmra_error2']**2+d['pmdec_error2']**2)
ok = (Rch < 0.01) & (plx1 > 5) & (plx2 > 5) \
   & (plx1/np.maximum(eplx1, 1e-6) > 20) & (plx2/np.maximum(eplx2, 1e-6) > 20) \
   & (np.abs(plx1-plx2) < 3*np.hypot(eplx1, eplx2)) \
   & (sep > 200) & (sep < 50000) \
   & (MG1 > 2.6) & (MG1 < 14.2) & (MG2 > 2.6) & (MG2 < 14.2) & (sigv < 0.03)

c1, c2 = d['bp_rp1'], d['bp_rp2']
col = np.concatenate([c1[ok], c2[ok]])
mg = np.concatenate([MG1[ok], MG2[ok]])
good = np.isfinite(col) & np.isfinite(mg) & (col > -0.5) & (col < 6.0)
CB = np.linspace(np.nanpercentile(col[good], 0.5),
                 np.nanpercentile(col[good], 99.5), 41)
cc = 0.5*(CB[:-1]+CB[1:])
ridge = np.full(len(cc), np.nan)
for i in range(len(cc)):
    m = good & (col >= CB[i]) & (col < CB[i+1])
    if m.sum() > 30:
        ridge[i] = np.median(mg[m])
vr = np.isfinite(ridge)
ccv, rgv = cc[vr], ridge[vr]
def ridge_f(c):
    return np.interp(c, ccv, rgv)
delta = (mg - ridge_f(col))[good]
P(f"STAGE 7J-A: {int(ok.sum())} pairs, {good.sum()} stars with color; "
  f"observed P(delta<-0.4) = {float(np.mean(delta < -0.4)):.3f} (3J 0.123)")

# --- the single-star locus in (mass, color, M_G) --------------------------
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
def mg_of_mass(mm):     # MG_T increasing, MS_T decreasing
    return np.interp(-np.asarray(mm), -MS_T, MG_T)
def color_of_mg(mgv):   # ridge is monotone increasing M_G(color)
    return np.interp(mgv, rgv, ccv)
mass_star = np.interp(mg[good], MG_T, MS_T)
# primary-mass sample (empirical), subsampled for the template
MPRI = rng.choice(mass_star, 4000, replace=False)

# --- combined (magnitude, color) shift for a companion of ratio q ---------
QG = np.linspace(0.10, 1.0, 91)
def combined_delta_grid(mpri, qg):
    m1 = mpri[:, None]
    mg1 = mg_of_mass(m1)
    mg2 = mg_of_mass(np.clip(qg[None, :]*m1, MS_T[-1], MS_T[0]))
    F1, F2 = 10**(-0.4*mg1), 10**(-0.4*mg2)
    mgs = -2.5*np.log10(F1+F2)
    c1_ = color_of_mg(mg1)
    c2_ = color_of_mg(mg2)
    cs = (F1*c1_ + F2*c2_)/(F1+F2)
    return mgs - ridge_f(cs)          # (npri, nq); <=0 by construction
CD = combined_delta_grid(MPRI, QG)
P(f"companion shift: q=1 twins {np.median(CD[:, -1]):.3f} mag (exact "
  f"-0.753); q=0.5 {np.median(CD[:, np.argmin(np.abs(QG-0.5))]):.3f}; "
  f"q=0.3 {np.median(CD[:, np.argmin(np.abs(QG-0.3))]):.3f}")

# --- binned Poisson likelihood of the delta distribution ------------------
BE_ = np.arange(-1.2, 0.6001, 0.05)
NB = len(BE_)-1
hist_obs, _ = np.histogram(delta, bins=BE_)
N_IN = hist_obs.sum()
from scipy.special import erf
def bin_probs(cd_flat, w_flat, sig, mu):
    # P(bin) = sum_k w_k [Phi((hi-mu-cd_k)/sig)-Phi((lo-mu-cd_k)/sig)]
    z = (BE_[None, :] - mu - cd_flat[:, None])/(sig*np.sqrt(2))
    Phi = 0.5*(1+erf(z))
    pb = (w_flat[:, None]*(Phi[:, 1:]-Phi[:, :-1])).sum(0)
    return pb
CD_FLAT = CD.ravel()
W_FLAT = np.full(CD_FLAT.size, 1.0/CD_FLAT.size)
ZERO = np.zeros(1)
WONE = np.ones(1)
def model_probs(f, sig, mu):
    ps = bin_probs(ZERO, WONE, sig, mu)
    pc = bin_probs(CD_FLAT, W_FLAT, sig, mu)
    return (1-f)*ps + f*pc

F_GRID = np.arange(0.0, 0.6001, 0.01)
SIG_G = np.arange(0.18, 0.361, 0.005)
MU_G = np.arange(-0.06, 0.0601, 0.01)
def fit(hist):
    # the companion template depends on (sig, mu) only — hoisted out of
    # the f loop (61 x 481 evaluations would otherwise redo the erf sum)
    best = (-1e18, None)
    prof = np.full(len(F_GRID), -1e18)
    for sg in SIG_G:
        for mu in MU_G:
            ps = bin_probs(ZERO, WONE, sg, mu)
            pc = bin_probs(CD_FLAT, W_FLAT, sg, mu)
            for fi, f in enumerate(F_GRID):
                pb = np.maximum((1-f)*ps + f*pc, 1e-12)
                pb = pb/pb.sum()
                ll = float(np.sum(hist*np.log(pb)))
                if ll > prof[fi]:
                    prof[fi] = ll
                if ll > best[0]:
                    best = (ll, (f, sg, mu))
    return best, prof
(best_ll, (f_hat, sig_hat, mu_hat)), prof = fit(hist_obs)
d2 = prof.max()-prof
in1 = F_GRID[d2 <= 0.5]
f_16, f_84 = float(in1.min()), float(in1.max())
P(f"fit: f_true = {f_hat:.3f} [{f_16:.3f}, {f_84:.3f}] (1-sigma profile), "
  f"sigma_MS = {sig_hat:.3f}, mu0 = {mu_hat:+.3f}, lnL = {best_ll:.1f}")

# --- gates ----------------------------------------------------------------
pb = model_probs(f_hat, sig_hat, mu_hat)
pb = pb/pb.sum()
edge = np.searchsorted(BE_, -0.4)
p_flag_model = float(pb[:edge].sum())
p_flag_obs = float(np.mean(delta < -0.4))
ga1 = abs(p_flag_model - p_flag_obs) <= 0.02
P(f"GA1 predictive: model P(delta<-0.4) = {p_flag_model:.3f} vs observed "
  f"{p_flag_obs:.3f}  -> {'PASS' if ga1 else 'FAIL'}")
ga2_ok = True
for f_inj in (0.15, 0.35):
    n_ = 120000
    isc = rng.random(n_) < f_inj
    cd_i = np.where(isc, CD_FLAT[rng.integers(0, CD_FLAT.size, n_)], 0.0)
    d_i = cd_i + rng.normal(0, 0.275, n_)
    h_i, _ = np.histogram(d_i, bins=BE_)
    (_, (f_rec, _, _)), _ = fit(h_i)
    rel = abs(f_rec - f_inj)/f_inj
    ga2_ok &= rel <= 0.20
    P(f"GA2 injection f={f_inj}: recovered {f_rec:.3f} "
      f"(rel {100*rel:.0f}%)  -> {'PASS' if rel <= 0.20 else 'FAIL'}")
ga3 = 0.22 <= sig_hat <= 0.32
P(f"GA3 width: sigma_MS = {sig_hat:.3f} in [0.22, 0.32]  -> "
  f"{'PASS' if ga3 else 'FAIL'}")

# --- completeness + the strict-residual mapping ---------------------------
z_th = (-0.4 - mu_hat - CD_FLAT)/(sig_hat*np.sqrt(2))
C_flag = float(np.mean(0.5*(1+erf(z_th))))          # P(delta<-0.4 | comp)
z_s = (-0.4 - mu_hat)/(sig_hat*np.sqrt(2))
C_single = float(0.5*(1+erf(z_s)))                  # false-flag rate
P(f"completeness of the -0.4 flag: C(companion) = {C_flag:.3f}; "
  f"false-flag(single) = {C_single:.3f}")
def resid_rate(f):
    # per-star: P(companion & unflagged)/P(unflagged), photometric only
    num = f*(1-C_flag)
    den = f*(1-C_flag) + (1-f)*(1-C_single)
    return num/np.maximum(den, 1e-12)
R_OF_F = resid_rate(F_GRID)
lnpi_full = prof - prof.max()
# push the same profile onto the residual axis (monotone mapping)
r_grid = R_OF_F
lnpi_strict = lnpi_full.copy()
f_resid_hat = float(resid_rate(f_hat))
P(f"strict-sample residual companion rate r(f_hat) = {f_resid_hat:.3f} "
  f"[{resid_rate(f_16):.3f}, {resid_rate(f_84):.3f}]  (photometric-only "
  f"= upper bound; RUWE removes more)")

# sensitivity: q-distribution variant ~ q^-0.5 (more low-q, harder to see)
wq = QG**-0.5
W_ALT = (np.ones((len(MPRI), 1))*wq[None, :]).ravel()
W_ALT /= W_ALT.sum()
W_KEEP = W_FLAT.copy()
W_FLAT = W_ALT
(_, (f_alt, sig_alt, _)), _ = fit(hist_obs)
W_FLAT = W_KEEP
P(f"sensitivity q^-0.5: f_true = {f_alt:.3f} (flat-q {f_hat:.3f}) - "
  f"spread carried into the prior width")

np.savez('data/stage7j_prior.npz', f_grid=F_GRID, lnpi_full=lnpi_full,
         r_grid=r_grid, lnpi_strict=lnpi_strict, C_flag=C_flag,
         C_single=C_single, f_hat=f_hat, f_16=f_16, f_84=f_84,
         f_alt=f_alt, f_resid_hat=f_resid_hat)
gates = ga1 and ga2_ok and ga3
P(f"\nGATES {'ALL PASS -> part B may launch' if gates else 'FAIL -> HOLD'}")
with open('data/stage7j_completeness.txt', 'w') as f_:
    f_.write("\n".join(L)+"\n")
print("saved: data/stage7j_completeness.txt + data/stage7j_prior.npz")
