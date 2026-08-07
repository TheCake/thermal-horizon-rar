"""
PREDICTOR BATTERY for horizon-coupled inertia (and its rivals).
Each section = one observable that separates: (A) inertia couples to H(z),
(B) inertia couples to Lambda (Verlinde-style, constant), (C) superfluid
dark matter (condensate halo + phonon force), (D) plain LCDM.
"""
import math
G, c, hbar, kB = 6.674e-11, 2.998e8, 1.0546e-34, 1.381e-23
H0 = 2.27e-18
M_sun, kpc, AU = 1.989e30, 3.086e19, 1.496e11
a0 = 1.2e-10
a_c = a0/2

def mu(x): return (math.sqrt(x*x+1)-1)/x if x>0 else 0.0
def g_eff(gN, gext=0.0):
    lo, hi = gN, gN + 4*math.sqrt(gN*a_c) + 4*a_c
    for _ in range(200):
        m_ = 0.5*(lo+hi)
        if m_*mu((m_+gext)/a_c) < gN: lo = m_
        else: hi = m_
    return 0.5*(lo+hi)

def hdr(s): print("\n" + "="*74 + f"\n{s}\n" + "="*74)

hdr("1. THE TERRESTRIAL WINDOW (Ignatiev-type 'SHLEM' event)")
# Modified-INERTIA versions trigger on total kinematic acceleration.
# Earth's spin-centrifugal and orbital accelerations can momentarily cancel
# at specific high-latitude spots at the equinoxes:
om_spin, om_orb = 7.27e-5, 1.99e-7
a_cf_80 = om_spin**2 * 6.37e6 * math.cos(math.radians(80))
a_orb   = om_orb**2 * 1.496e11
print(f"spin-centrifugal accel at 80deg lat: {a_cf_80:.2e} m/s^2")
print(f"orbital accel around Sun:            {a_orb:.2e} m/s^2   <- same order!")
window_x = 2*a0/om_spin**2
window_t = 2*a0/(a_orb*om_orb)
print(f"spatial size of the |a| < a0 null spot:  ~{100*window_x:.0f} cm")
print(f"time the spot condition holds:           ~{window_t:.1f} s, twice a year")
print("-> modified-inertia theories predict a gravimeter/atom-interferometer")
print("   GLITCH in a ~cm spot at high latitude at equinox. LCDM and modified-")
print("   GRAVITY theories predict exactly nothing. A weekend-scale experiment.")

hdr("2. ENVIRONMENT DEPENDENCE (external field effect = SEP violation)")
M_b = 5e10*M_sun; r = 30*kpc
gN = G*M_b/r**2
for gext_units, env in [(0.02, "deep void"), (1.8, "solar neighborhood"), (5.0, "cluster outskirts")]:
    v = math.sqrt(g_eff(gN, gext_units*a0)*r)/1000
    print(f"  identical galaxy, g_ext = {gext_units:>4} a0 ({env:>18}): v(30kpc) = {v:.0f} km/s")
print("-> same galaxy rotates measurably slower in crowded environments.")
print("   NO dark-matter theory does this (halos don't care about uniform")
print("   external fields -- strong equivalence principle). Chae et al. claim")
print("   this signature exists in SPARC data. Sharpen it: it must CORRELATE")
print("   with large-scale structure maps, galaxy by galaxy.")

hdr("3. TIDAL DWARF GALAXIES (recycled galaxies = DM-free by construction)")
M_tdg = 1e9*M_sun; r = 5*kpc
gN = G*M_tdg/r**2
vN  = math.sqrt(gN*r)/1000
vM  = math.sqrt(g_eff(gN)*r)/1000
print(f"1e9 M_sun tidal dwarf at 5 kpc: LCDM (no halo captured): v = {vN:.0f} km/s")
print(f"                                horizon-inertia:          v = {vM:.0f} km/s")
print("-> galaxies born from tidal debris inherit NO dark halo in LCDM, but")
print("   modified dynamics applies to them identically. Observed TDGs")
print("   (NGC5291 system etc.) show the boosted value. Keep pushing: JWST can")
print("   find more. This is the cleanest existing discriminator.")

hdr("4. LENSING vs KINEMATICS SPLIT (kills or crowns superfluid DM)")
print("Superfluid-DM: the extra force is PHONON-mediated -> pulls matter but")
print("NOT light. Prediction: kinematic RAR boosted, lensing RAR follows only")
print("actual mass -> the two RARs SPLIT inside the superfluid core.")
print("True modified inertia/gravity: both follow the SAME effective potential")
print("-> RARs identical. Weak-lensing RAR (Brouwer+ 2021) already tracks the")
print("kinematic one down to ~1e-12.5 m/s^2 -- superfluid DM is being cornered.")

hdr("5. THE REDSHIFT DICHOTOMY (separates the two horizon couplings)")
print(f"{'z':>3} {'a0 if ~H(z)':>12} {'a0 if ~sqrt(Lambda)':>20}")
for z in (0,1,2,3):
    Hz = math.sqrt(0.3*(1+z)**3+0.7)
    print(f"{z:>3} {Hz:>11.2f}x {'1.00x':>20}")
print("-> coupling to the EPOCH (H) vs to the COSMOLOGICAL CONSTANT (Verlinde)")
print("   give identical local physics today but diverge in lookback time.")
print("   High-z Tully-Fisher zero point separates them. Nothing else does.")

hdr("6. VACUUM-ENGINEERING CROSS-CHECK (ties back to the gravity-modification bounds)")
rho_cas = 4.33e4                      # |Casimir| J/m^3 at 10 nm gaps
frac = rho_cas/ (1150*c**2)           # vs half-volume silicon plates, kg/m^3
print(f"If inertia is a vacuum reaction, Casimir-structured matter (10 nm gaps)")
print(f"has its local vacuum edited by delta-m/m ~ {frac:.0e}.")
print(f"MICROSCOPE already confirms universal free fall at 2.7e-15, and binding-")
print(f"energy differences (1e-3 level) fall universally to 1e-13 -> naive")
print(f"'vacuum-drag inertia' is ALREADY cornered to fine-tuning. Any real")
print(f"mechanism must couple to the HORIZON (global), not the local mode")
print(f"structure -- which is why lab-scale Casimir gravity modification never had a chance,")
print(f"and why the galactic-scale condensate (superfluid DM) is the version")
print(f"of the superconductor dream that is still scientifically alive.")
