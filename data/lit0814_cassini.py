"""lit0814_cassini.py -- pin the "x Cassini" arithmetic against BOTH quadrupole bounds.

Literature refresh 2026-08-14 (Paper 2, Section 9.1).

Primary reads behind the numbers below:
  Hees, A., Folkner, W. M., Jacobson, R. A., & Park, R. S. 2014, PRD 89, 102002
      -- the bound the program's archived stages (4K, 5I, 5S, 5X, 6I) were computed
         against: Q2 = (3 +/- 3) e-27 s^-2.
  Park, R. S., Hees, A., Famaey, B., Desmond, H., & Durakovic, A. 2026,
      PRD 114, 024066 (arXiv:2602.17884), "Improved constraints on modified
      Newtonian gravity from Cassini radio tracking data", published 2026-07-28.
      Abstract, verbatim: "we find Q2 = (1.6 +/- 1.8) x 10^-27 s^-2 (1-sigma),
      representing an improvement of 40% over previous estimates. ... we update
      previously acknowledged tensions with external galaxy rotation curves, now
      leading to discrepancies at the 3-15 sigma level ... The updated Q2 posterior
      finally confirms that Solar System measurements provide stronger constraints
      than current wide-binary data on classical modified gravity versions of MOND."
      (Crossref 10.1103/r7n8-kw38: Phys. Rev. D 114, 024066, 2026-07-28.)

Our Q2 is unchanged: it is the archived output of the program's own axisymmetric
QUMOND external-field solver (calcs/stage4k_quadrupole.py, cross-validated at the
15% level against Blanchet & Novak 2011). Only the denominator moves.

Bound convention: the program's published "x Cassini" ratios use the 2-sigma
ceiling, central + 2*sigma. Reproduced below for Hees+14 as a regression on the
paper's existing 3.8 / 4.3 figures before the new ratios are quoted.

Output is ASCII only (Windows cp1252 console).
"""

# ---------------------------------------------------------------- our numerator
Q2_ALPHA1 = 3.4e-26   # s^-2, parameter-free galactic calibration alpha = 1
Q2_ALPHA115 = 3.9e-26  # s^-2, superseded fenced wide-binary amplitude alpha = 1.15

# the function-family amplitude-locked band, as ratios against Hees+14 (archived)
BAND_HEES = (4.0, 5.8)

# ---------------------------------------------------------------- the two bounds
BOUNDS = {
    "Hees+14  (PRD 89, 102002)": (3.0e-27, 3.0e-27),
    "Park+26  (PRD 114, 024066)": (1.6e-27, 1.8e-27),
}


def ceiling(central, sigma, k=2.0):
    return central + k * sigma


def line(ch="-", n=72):
    return ch * n


print(line("="))
print("Q2 QUADRUPOLE RATIOS -- archived bound vs current bound (2026-08-14)")
print(line("="))
print()
print("Numerator (unchanged, stage4k_quadrupole.py):")
print("  Q2(alpha = 1.00, galactic calibration)      = %.2e s^-2" % Q2_ALPHA1)
print("  Q2(alpha = 1.15, superseded fenced binary)  = %.2e s^-2" % Q2_ALPHA115)
print()

results = {}
for name, (central, sigma) in BOUNDS.items():
    c1 = ceiling(central, sigma, 1.0)
    c2 = ceiling(central, sigma, 2.0)
    r1 = Q2_ALPHA1 / c2
    r115 = Q2_ALPHA115 / c2
    results[name] = (c2, r1, r115)
    print(line())
    print(name)
    print("  bound            : (%.1f +/- %.1f) e-27 s^-2  (1-sigma)" % (central * 1e27, sigma * 1e27))
    print("  1-sigma ceiling  : %.2e s^-2" % c1)
    print("  2-sigma ceiling  : %.2e s^-2   <-- the program's ratio denominator" % c2)
    print("  ratio at alpha=1.00 : %.2f x" % r1)
    print("  ratio at alpha=1.15 : %.2f x" % r115)
print(line())
print()

# ---------------------------------------------------------------- regression
old = results["Hees+14  (PRD 89, 102002)"]
new = results["Park+26  (PRD 114, 024066)"]

print("REGRESSION against the figures already in Paper 2 draft 0.5 (Hees+14):")
print("  paper says 3.8x at alpha=1.00 ; computed %.2f  -> %s"
      % (old[1], "OK" if abs(old[1] - 3.8) < 0.05 else "MISMATCH"))
print("  paper says 4.3x at alpha=1.15 ; computed %.2f  -> %s"
      % (old[2], "OK" if abs(old[2] - 4.3) < 0.05 else "MISMATCH"))
print()

print("THE UPDATE (Park+26 is the current bound):")
print("  alpha = 1.00 : %.1f x  (was %.1f x)" % (new[1], old[1]))
print("  alpha = 1.15 : %.1f x  (was %.1f x)" % (new[2], old[2]))
print("  tightening factor on every archived 'x Cassini' number: %.3f"
      % (old[0] / new[0]))
print()

print("QUOTED RATIOS FOR THE PAPER (one decimal, as printed in Section 9.1):")
print("  RATIO_ALPHA1_NEW   = %.1f" % new[1])
print("  RATIO_ALPHA115_NEW = %.1f" % new[2])
print("  RATIO_ALPHA1_OLD   = %.1f" % old[1])
print("  RATIO_ALPHA115_OLD = %.1f" % old[2])
print()

scale = old[0] / new[0]
band_new = (BAND_HEES[0] * scale, BAND_HEES[1] * scale)
print("AMPLITUDE-LOCKED FUNCTION-FAMILY BAND (Section 9.1, second paragraph):")
print("  against Hees+14 (archived) : %.1f - %.1f x" % BAND_HEES)
print("  against Park+26 (current)  : %.1f - %.1f x" % band_new)
print("  (the band is a ratio to the same denominator, so it rescales exactly)")
print()

print("CONTEXT, quoted from Park et al. (2026) and not recomputed here:")
print("  - rotation-curve tension now 3-15 sigma depending on mass modeling/subset")
print("  - Milky Way MOND boost at the Sun bounded to 2% (95% confidence)")
print("  - Solar System now constrains classical MG-MOND more strongly than")
print("    current wide-binary data")
print()
print(line("="))
print("Direction of the update: the tension roughly doubles. It binds the FIELD")
print("formulation harder. The trajectory-formulation carve-out (Section 9.3) is")
print("untouched -- it is a statement about which accelerations Saturn's worldline")
print("samples, not about the size of the bound.")
print(line("="))
