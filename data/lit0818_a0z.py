"""
lit0818_a0z.py -- convention-checked comparison of the horizon temperature lock
a0(z) = c*H(z)/(2*pi) against the first direct low-redshift a0(z) fit.

External source (primary-read 2026-08-18):
  Varasteanu, Jarvis, Desmond, Ponomareva, Yasin et al. 2026, arXiv:2608.03576
  "MIGHTEE-HI / LADUMA: Investigating the link between baryons and dynamics with
  130 resolved HI-selected galaxies" (submitted 2026-08-04), 130 HI-selected
  galaxies with resolved kinematics to z ~ 0.09.
  Their linear parametrization (their Eq. 14):  a(z) = a0 + a1 * z.
    sample alone      : a0 = (1.54 +/- 0.11)e-10, a1 = (-1.60 +/- 2.33)e-10 m/s^2
    + SPARC as z~0 anchor: a0 = (1.15 +/- 0.02)e-10, a1 = (5.23 +/- 1.05)e-10,
                          "a formal 5.0 sigma preference".
  Origin of the evolution claim: Varasteanu et al. 2025, arXiv:2504.20857,
  "the first tentative evidence for redshift evolution in the acceleration scale".

THIS SCRIPT IS THE SOURCE FOR EVERY NUMBER QUOTED FROM THIS COMPARISON IN
papers/paper2_rar_coefficients.md, papers/paper3_mechanism.md, and PREDICTIONS.md.
Nothing not printed below may be printed in those files.

Run:  py data/lit0818_a0z.py > data/lit0818_a0z.txt
"""

import math

# ---------------------------------------------------------------- conventions
# Stated cosmology for the lock curve (declared, not fitted):
#   flat LambdaCDM, H0 = 70 km/s/Mpc, Omega_m = 0.3, Omega_Lambda = 0.7.
# This is the cosmology used for the lock throughout; the comparison targets'
# own adopted cosmologies are not re-derived here (see the caveat block).
H0_KM_S_MPC = 70.0
OMEGA_M = 0.3
OMEGA_L = 1.0 - OMEGA_M

C_M_S = 2.99792458e8          # exact, SI
MPC_M = 3.0856775814913673e22  # exact IAU-derived metres per megaparsec

H0_SI = H0_KM_S_MPC * 1.0e3 / MPC_M       # s^-1
UNIT = 1.0e-10                            # all a0 values reported in 1e-10 m/s^2

ZMAX = 0.09                               # their sample's redshift reach


def E(z):
    """H(z)/H0 for flat LambdaCDM."""
    return math.sqrt(OMEGA_M * (1.0 + z) ** 3 + OMEGA_L)


def a0_lock(z):
    """The temperature lock a0(z) = c H(z) / (2 pi), in SI."""
    return C_M_S * H0_SI * E(z) / (2.0 * math.pi)


def sig(value, ref, err):
    """Signed distance of `value` from `ref` in units of `err`."""
    return (value - ref) / err


line = "-" * 78

print("lit0818_a0z.py -- the temperature lock vs the first direct a0(z) fit")
print("source paper: arXiv:2608.03576 (Varasteanu et al. 2026, MIGHTEE-HI/LADUMA)")
print("origin of the evolution claim: arXiv:2504.20857 (Varasteanu et al. 2025)")
print(line)

# ------------------------------------------------------------ 0. the lock curve
a0_0 = a0_lock(0.0)
a0_zmax = a0_lock(ZMAX)

# Analytic derivative at z = 0: dH/dz|0 = H0 * (3/2) Omega_m, so
#   da0/dz|0 = a0(0) * (3/2) Omega_m  (with E(0) = 1).
dlnE_dz_0 = 1.5 * OMEGA_M
deriv0 = a0_0 * dlnE_dz_0

# Numerical check of the same derivative (central difference).
h = 1.0e-6
deriv0_num = (a0_lock(h) - a0_lock(-h)) / (2.0 * h)

# Mean (secant) slope over the sample's own window 0 < z < 0.09.  This is the
# quantity directly comparable to their a1, because their a1 is the coefficient
# of a LINEAR form fitted across that window, not a derivative at a point.
secant = (a0_zmax - a0_0) / ZMAX

# Least-squares linear fit of the lock curve over the same window, as a check
# that the secant is a fair linearization (uniform z weighting, 4001 nodes).
N = 4001
zs = [ZMAX * i / (N - 1) for i in range(N)]
ys = [a0_lock(z) for z in zs]
zbar = sum(zs) / N
ybar = sum(ys) / N
num = sum((z - zbar) * (y - ybar) for z, y in zip(zs, ys))
den = sum((z - zbar) ** 2 for z in zs)
lsq_slope = num / den
lsq_intercept = ybar - lsq_slope * zbar

print("0. THE LOCK CURVE (stated cosmology)")
print("   flat LambdaCDM, H0 = %.1f km/s/Mpc, Omega_m = %.2f, Omega_Lambda = %.2f"
      % (H0_KM_S_MPC, OMEGA_M, OMEGA_L))
print("   c = %.8e m/s, H0 = %.6e s^-1" % (C_M_S, H0_SI))
print("   a0_lock(0)      = %.4f e-10 m/s^2" % (a0_0 / UNIT))
print("   a0_lock(%.2f)   = %.4f e-10 m/s^2" % (ZMAX, a0_zmax / UNIT))
print("   E(%.2f) = H(z)/H0 = %.5f" % (ZMAX, E(ZMAX)))
print("   da0/dz at z = 0 (analytic, = a0(0)*1.5*Omega_m) = %+.4f e-10 per unit z"
      % (deriv0 / UNIT))
print("   da0/dz at z = 0 (numerical check)               = %+.4f e-10 per unit z"
      % (deriv0_num / UNIT))
print("   mean (secant) slope over 0 < z < %.2f           = %+.4f e-10 per unit z"
      % (ZMAX, secant / UNIT))
print("   least-squares slope over the same window        = %+.4f e-10 per unit z"
      % (lsq_slope / UNIT))
print("   least-squares intercept                         = %+.4f e-10" %
      (lsq_intercept / UNIT))
print()
print("   PRIMARY QUANTITY: a1_lock = %+.2f e-10 m/s^2 per unit z" % (secant / UNIT))
print("   (the secant over their own window, the like-for-like partner of their a1;")
print("    the point derivative at z = 0 is %+.2f e-10 and the two bracket the"
      % (deriv0 / UNIT))
print("    linearization ambiguity at the %.3f e-10 level)"
      % (abs(secant - deriv0) / UNIT))
print(line)

a1_lock = secant / UNIT   # in 1e-10 units, used everywhere below

# --------------------------------------------- their measured values (verbatim)
# sample alone
S_A0, S_A0E = 1.54, 0.11
S_A1, S_A1E = -1.60, 2.33
# SPARC-anchored combined fit
C_A0, C_A0E = 1.15, 0.02
C_A1, C_A1E = 5.23, 1.05

print("1. THEIR MEASURED VALUES (verbatim from arXiv:2608.03576, in 1e-10 m/s^2)")
print("   pure HI-selected sample alone : a0 = %.2f +/- %.2f, a1 = %+.2f +/- %.2f"
      % (S_A0, S_A0E, S_A1, S_A1E))
print("   + SPARC as z~0 anchor         : a0 = %.2f +/- %.2f, a1 = %+.2f +/- %.2f"
      % (C_A0, C_A0E, C_A1, C_A1E))
print("   their label on the second row : 'a formal 5.0 sigma preference for an")
print("   acceleration scale that increases with redshift'")
print(line)

# ------------------------------------- (a) lock vs the pure-sample slope band
d_lock_sample = sig(a1_lock, S_A1, S_A1E)
d_zero_sample = sig(0.0, S_A1, S_A1E)
lo, hi = S_A1 - S_A1E, S_A1 + S_A1E
lo2, hi2 = S_A1 - 2 * S_A1E, S_A1 + 2 * S_A1E

print("(a) LOCK SLOPE vs THEIR PURE-SAMPLE BAND  a1 = %+.2f +/- %.2f" % (S_A1, S_A1E))
print("    1-sigma band : [%+.2f, %+.2f]   2-sigma band : [%+.2f, %+.2f]"
      % (lo, hi, lo2, hi2))
print("    a1_lock = %+.2f  ->  distance %+.2f sigma  ->  INSIDE 1 sigma: %s"
      % (a1_lock, d_lock_sample, "YES" if lo <= a1_lock <= hi else "NO"))
print("    a1 = 0  ->  distance %+.2f sigma  ->  INSIDE 1 sigma: %s"
      % (d_zero_sample, "YES" if lo <= 0.0 <= hi else "NO"))
print("    VERDICT: the pure-sample band contains BOTH no evolution and the lock's")
print("    predicted slope; it separates them by %.2f sigma of its own width and so"
      % abs(d_lock_sample - d_zero_sample))
print("    discriminates between them at no useful grade.")
print(line)

# ------------------------------ (b) lock vs the SPARC-anchored slope
d_lock_comb = sig(C_A1, a1_lock, C_A1E)
d_zero_comb = sig(C_A1, 0.0, C_A1E)
ratio_comb = C_A1 / a1_lock

print("(b) LOCK SLOPE vs THEIR SPARC-ANCHORED FIT  a1 = %+.2f +/- %.2f"
      % (C_A1, C_A1E))
print("    their value sits %+.2f sigma ABOVE a1_lock = %+.2f"
      % (d_lock_comb, a1_lock))
print("    their value sits %+.2f sigma above zero on the same error bar" % d_zero_comb)
print("    (their own 5.0 sigma label is computed inside their fit and is not")
print("    reproduced here)")
print("    face-value ratio their_a1 / a1_lock = %.1f" % ratio_comb)
print("    VERDICT: at face value the anchored fit prefers a slope %.1f sigma above"
      % d_lock_comb)
print("    the lock; the sign agrees with the lock and disagrees with any falling")
print("    prediction.")
print(line)

# ---------------------- (c) the zero-point-offset reading, in their own numbers
rise_across_window = C_A1 * ZMAX          # implied rise of the anchored fit
intercept_gap = S_A0 - C_A0               # their own two intercepts
residual = rise_across_window - intercept_gap
slope_from_pure_offset = intercept_gap / ZMAX
d_anchored_from_offset = sig(C_A1, slope_from_pure_offset, C_A1E)

print("(c) THE ZERO-POINT-OFFSET READING (their numbers only)")
print("    anchored slope x window   : %+.2f x %.2f = %+.4f e-10"
      % (C_A1, ZMAX, rise_across_window))
print("    their own intercept gap   : %.2f - %.2f  = %+.4f e-10"
      % (S_A0, C_A0, intercept_gap))
print("    residual                  : %+.4f e-10 (%.0f%% of the gap)"
      % (residual, 100.0 * residual / intercept_gap))
print("    a slope reproducing the intercept gap alone across the window:")
print("      %+.4f / %.2f = %+.2f e-10 per unit z" % (intercept_gap, ZMAX,
                                                      slope_from_pure_offset))
print("    the anchored a1 sits %+.2f sigma from that pure-offset slope"
      % d_anchored_from_offset)
print("    VERDICT: the rise the anchored fit implies across its own redshift range")
print("    is accounted for, to within %.0f%% of the gap, by the offset between the"
      % (100.0 * residual / intercept_gap))
print("    two intercepts the same paper reports; an inter-sample zero-point")
print("    difference and a slope are not separated by this lever arm.")
print(line)

# ------------------------------------------------- supplementary: intercepts
d_lock_S = sig(S_A0, a0_0 / UNIT, S_A0E)
d_lock_C = sig(C_A0, a0_0 / UNIT, C_A0E)

print("SUPPLEMENTARY (intercepts; systematics-dominated, reported for completeness)")
print("    a0_lock(0) = %.2f e-10" % (a0_0 / UNIT))
print("    their sample-alone intercept   %.2f +/- %.2f -> %+.2f sigma from the lock"
      % (S_A0, S_A0E, d_lock_S))
print("    their SPARC-anchored intercept %.2f +/- %.2f -> %+.2f sigma from the lock"
      % (C_A0, C_A0E, d_lock_C))
print("    the gap between their OWN two intercepts (%.2f e-10) exceeds the distance"
      % intercept_gap)
print("    from the lock to the anchored intercept (%.2f e-10) by a factor %.1f,"
      % (abs(C_A0 - a0_0 / UNIT), intercept_gap / abs(C_A0 - a0_0 / UNIT)))
print("    which is the size of the inter-sample calibration unknown they flag.")
print(line)

print("CAVEATS ON THE COMPARISON")
print(" 1. The lock curve is evaluated at the stated cosmology above. Their fits")
print("    adopt their own cosmology; over 0 < z < %.2f the choice moves the lock" % ZMAX)
print("    slope far less than the quoted uncertainties on either of their a1 rows.")
print(" 2. Their a1 is the coefficient of a linear form; the lock is not linear.")
print("    The secant and the point derivative differ by %.3f e-10, which is the"
      % (abs(secant - deriv0) / UNIT))
print("    linearization ambiguity and is negligible against their error bars.")
print(" 3. No fit-direction or selection correction is applied here. Their own text")
print("    flags fit direction, selection effects, and the constant mass-to-light")
print("    assumption in the anchor sample as the limiting systematics.")
print(line)
print("END lit0818_a0z.py")
