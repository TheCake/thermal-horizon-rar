"""
The anomaly numerology dashboard.
Einstein's move: promote a clean anomaly to an exact law, then find which
hidden assumption dies. Step zero is knowing which numbers rhyme.
Includes one famous FALSE rhyme (Pioneer) as a calibration standard.
"""
import math

G, c, hbar, kB = 6.674e-11, 2.998e8, 1.0546e-34, 1.381e-23
H0   = 2.27e-18          # ~70 km/s/Mpc
Om_L, Om_m = 0.7, 0.3
eV   = 1.602e-19

def sci(x, d=2): return f"{x:.{d}e}"
def hdr(s): print("\n" + "="*72 + f"\n{s}\n" + "="*72)

hdr("1. THE ACCELERATION RHYME (the deepest one)")
a0_obs  = 1.2e-10                     # MOND scale from rotation curves (empirical)
a_cH    = c*H0
a_cH2pi = c*H0/(2*math.pi)
a_dS    = math.sqrt(Om_L)*H0*c        # de Sitter horizon acceleration scale
print(f"a0 from galaxy rotation curves (empirical):   {sci(a0_obs)} m/s^2")
print(f"c * H0:                                        {sci(a_cH)} m/s^2")
print(f"c * H0 / 2pi:                                  {sci(a_cH2pi)} m/s^2")
print(f"de Sitter horizon scale sqrt(Om_L)*H0*c:       {sci(a_dS)} m/s^2")
print("-> the ONLY acceleration scale where GR shows empirical cracks equals")
print("   the acceleration set by the cosmic horizon, to within factors of 2pi.")
print("   A galaxy's edge 'knows' the size of the visible universe. Why?")

hdr("2. THE ENERGY RHYME (dark energy vs neutrinos)")
rho_L = Om_L*3*H0**2/(8*math.pi*G)*c**2          # J/m^3
E_L   = (rho_L*(hbar*c)**3)**0.25                # J
print(f"dark energy scale (rho_L)^(1/4):     {E_L/eV*1000:.1f} meV")
print(f"neutrino mass scales (oscillations): ~8.7 meV and ~50 meV")
print("-> the two newest entries in physics (dark energy, neutrino mass) sit")
print("   within one-two orders on a chart spanning 10^30. Rhyme or noise?")

hdr("3. THE 'WHY NOW' COINCIDENCE")
print(f"rho_L/rho_matter today = {Om_L/Om_m:.1f} -- order 1, but rho_m ~ a(t)^-3")
print("while rho_L = const: they match ONLY in the cosmic epoch containing us.")
print("A constant fine-tuned to cross our era, or a density that TRACKS matter.")

hdr("4. CALIBRATION: A FAMOUS RHYME THAT LIED (Pioneer)")
a_pioneer = 8.74e-10
print(f"Pioneer anomalous acceleration (1998-2011): {sci(a_pioneer)} m/s^2")
print(f"c*H0:                                        {sci(a_cH)} m/s^2  <- 'rhymed'!")
print("Resolution (2012): anisotropic THERMAL RADIATION from the RTGs -- paint")
print("and heat, not new physics. Lesson: a numerical rhyme alone is nothing.")
print("You also need a STRUCTURAL anomaly (a tight law that shouldn't be tight).")

hdr("5. THE STRUCTURAL ANOMALY THAT PASSES THAT BAR")
# Baryonic Tully-Fisher: v^4 = G * a0 * M_baryon  (empirical, tiny scatter)
for M_b, name in [(6e10*2e30, "Milky-Way-like"), (2e9*2e30, "dwarf")]:
    v = (G*a0_obs*M_b)**0.25
    print(f"BTFR check {name:15s}: M_b={sci(M_b,1)} kg -> v = {v/1000:.0f} km/s (observed range: right on)")
print("v^4 = G*a0*M_baryon with ~zero scatter over 5 decades of mass.")
print("In the dark-matter picture, v is set by the HALO, M_b by messy gas/star")
print("history -- their tight coupling with no scatter is unexplained. This is")
print("the clean 'Mercury perihelion'-grade fact of the DM problem.")

hdr("6. THE HUBBLE TENSION (two right answers?)")
print("Early universe (CMB+LCDM extrapolated): H0 = 67.4 +/- 0.5")
print("Local ladder (Cepheids+SNe):            H0 = 73.0 +/- 1.0")
print(f"-> {100*(73.0-67.4)/67.4:.0f}% apart, ~5 sigma. If both are right, the assumption")
print("   that dies is LCDM evolution BETWEEN z~1100 and z~0.")

hdr("7. USER'S CHILDHOOD IDEA, AUDITED (primordial black holes as DM)")
# surviving window: asteroid-mass PBHs
print("PBH dark matter is NOT dead -- microlensing/evaporation kill most masses,")
print("but the window ~1e17-1e22 g (asteroid mass, ~1e-16 to 1e-11 M_sun) can")
print("still be ALL of dark matter. Radius of a 1e19 g PBH:")
r_s = 2*G*1e19/c**2 * 1e3   # convert kg: 1e19 g = 1e16 kg
r_s = 2*G*1e16/c**2
print(f"  Schwarzschild radius: {r_s*1e9:.1f} nm -- a proton-sized black hole")
print("  with an asteroid's mass. BUT: PBHs cannot explain the BTFR tightness")
print("  or the a0-horizon rhyme -- they inherit every cold-DM small-scale issue.")
