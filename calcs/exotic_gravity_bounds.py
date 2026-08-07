"""
First-principles order-of-magnitude bounds on exotic gravity-modification
claims (frame dragging, rotating superconductors, negative-energy
engineering) against the GR baseline.

Sections:
  1. What GR actually predicts for a lab-scale spinning disc (frame dragging / GEM)
  2. Why a rotating superconductor's REAL field is magnetic, not gravitational
     (London moment), and the coupling-ratio suppression
  3. What the Podkletnov / Tajmar claims would require vs. that GR baseline
  4. Negative mass / negative energy engineering budgets (Casimir, quantum bounds)

All numbers SI unless noted. This is order-of-magnitude physics: O(1) factors
from conventions (GEM factor 2s, shell vs disc geometry) are irrelevant at the
20+ orders of magnitude gaps involved.
"""

import math

# --- constants ---
G      = 6.674e-11        # m^3 kg^-1 s^-2
c      = 2.998e8          # m/s
hbar   = 1.0546e-34       # J s
e      = 1.602e-19        # C
m_e    = 9.109e-31        # kg
m_p    = 1.673e-27        # kg
eps0   = 8.854e-12        # F/m
g0     = 9.81             # m/s^2
k_B    = 1.381e-23

def sci(x, digits=2):
    return f"{x:.{digits}e}"

def header(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)

# ===========================================================================
header("1. GR BASELINE: frame dragging from a lab-scale spinning disc")
# Podkletnov-style disc: YBCO, ~270 mm diameter, ~10 mm thick, ~5000 RPM
R_disc   = 0.135                 # m
m_disc   = 5.0                   # kg  (~YBCO density * disc volume, generous)
rpm      = 5000.0
omega    = rpm * 2*math.pi/60    # rad/s
I_disc   = 0.5 * m_disc * R_disc**2
J_disc   = I_disc * omega        # angular momentum

# GEM dipole field on axis at height z above disc (order of magnitude):
#   B_g ~ 2 G J / (c^2 z^3)      [units: 1/s; the "gravitomagnetic" field]
z = 0.10
B_g = 2*G*J_disc/(c**2 * z**3)
print(f"disc: m={m_disc} kg, R={R_disc} m, {rpm:.0f} RPM  ->  J={sci(J_disc)} kg m^2/s")
print(f"gravitomagnetic field 10 cm above disc:  B_g ~ {sci(B_g)} s^-1")

# Force on a test mass moving at v through that field: a ~ 4 v B_g
v_test = 100.0
a_gem = 4*v_test*B_g
print(f"accel on test mass moving {v_test} m/s:   a ~ {sci(a_gem)} m/s^2  "
      f"= {sci(a_gem/g0)} g")

# Frame-dragging ratio inside a rotating shell (Brill-Cohen, weak field):
#   Omega_drag/omega ~ (4/3) GM/(R c^2)
drag_ratio = (4/3)*G*m_disc/(R_disc*c**2)
print(f"frame-drag ratio Omega/omega:            ~ {sci(drag_ratio)}")

# Compare: what fraction of Earth's g would even the disc's ORDINARY Newtonian
# attraction give right above it?
g_newton_disc = G*m_disc/z**2
print(f"disc's own Newtonian pull at 10 cm:      {sci(g_newton_disc)} m/s^2 "
      f"= {sci(g_newton_disc/g0)} g")

# ===========================================================================
header("2. ROTATING SUPERCONDUCTOR: the real effect (London moment) and the "
       "coupling suppression")
# A rotating superconductor spontaneously generates a REAL magnetic field:
#   B_London = -(2 m_e / e) * omega        (measured; used in Gravity Probe B)
B_london = 2*m_e/e * omega
print(f"London moment at {rpm:.0f} RPM:  B = {sci(B_london)} T  "
      f"({B_london*1e9:.1f} nT)  <-- real, measured physics")

# Gravito-electric vs electric coupling of the same Cooper pair current:
# ratio = 4 pi eps0 G m^2 / e^2
ratio_e  = 4*math.pi*eps0*G*m_e**2 / e**2
ratio_p  = 4*math.pi*eps0*G*m_p**2 / e**2
print(f"(grav coupling)/(EM coupling), electron pairs: {sci(ratio_e)}")
print(f"same for proton-mass carriers:                 {sci(ratio_p)}")
print("Any 'gravitational London moment' built from mass currents inherits a")
print("suppression of this order relative to the EM one -- ~40 orders down.")

# ===========================================================================
header("3. WHAT THE CLAIMS REQUIRE vs. THE GR BASELINE")
# Podkletnov claim: 0.3% - 2% weight reduction above the disc.
claim_frac = 0.02
a_claim = claim_frac*g0
print(f"Podkletnov claim: {claim_frac*100:.0f}% of g = {a_claim} m/s^2 repulsive/shielded")
print(f"  vs disc's own total Newtonian gravity:  x{sci(a_claim/g_newton_disc)} "
      f"stronger than ALL the disc's gravity")
print(f"  vs GR gravitomagnetic effect:           x{sci(a_claim/a_gem)} enhancement needed")

# Tajmar claim (2006-2008): accelerometers near a spinning Nb ring read
# ~ 1e-4 g during angular acceleration; GR predicts ~ 1e-31 g scale.
a_tajmar_claim = 1e-4*g0
ring_m, ring_R, ring_omega = 0.4, 0.07, 500.0
J_ring = ring_m*ring_R**2*ring_omega
Bg_ring = 2*G*J_ring/(c**2 * ring_R**3)
print(f"\nTajmar ring GR baseline: B_g ~ {sci(Bg_ring)} s^-1 "
      f"-> accel scale ~ {sci(4*ring_R*ring_omega*Bg_ring/g0)} g")
print(f"Tajmar claimed ~1e-4 g  ->  enhancement over GR ~ "
      f"{sci(a_tajmar_claim/(4*ring_R*ring_omega*Bg_ring*g0)/g0, 1)} (huge; see report)")

# ===========================================================================
header("4. NEGATIVE MASS / NEGATIVE ENERGY BUDGETS")
# To cancel Earth's g at distance r from a point device: M_neg = -g r^2 / G
for r in (1.0, 10.0):
    M_neg = g0*r**2/G
    print(f"negative mass to null g at r={r:>4} m:  {sci(M_neg)} kg "
          f"(~{sci(M_neg/2700)} m^3 of rock, negative)")

# Casimir energy density between plates, gap a:  rho = -pi^2 hbar c / (720 a^4)
for a_gap in (10e-9, 1e-9):
    rho_E = -math.pi**2*hbar*c/(720*a_gap**4)     # J/m^3
    rho_m = rho_E/c**2                             # kg/m^3
    M_needed = g0*1.0**2/G
    V_needed = M_needed/abs(rho_m)
    print(f"Casimir gap {a_gap*1e9:.0f} nm: rho = {sci(rho_E)} J/m^3 "
          f"= {sci(rho_m)} kg/m^3;  volume for r=1m null: {sci(V_needed)} m^3 "
          f"(cube {sci(V_needed**(1/3))} m/side)")
print("...and the plates themselves (~1e3 kg/m^3 positive) outweigh the negative")
print("   energy by ~1e15, so the NET mass of any Casimir device is positive.")

# Ford-Roman quantum inequality (flat space, massless field, sampling time t):
#   <rho> >~ - (3/32 pi^2) hbar c / (c t)^4   -- order hbar c/(ct)^4
for t in (1.0, 1e-6, 1e-9):
    rho_bound = hbar*c/( (c*t)**4 )
    print(f"quantum-inequality bound, sustained {sci(t,0)} s: "
          f"|rho| <~ {sci(rho_bound)} J/m^3")
print("Negative energy is real (Casimir, squeezed light) but QFT taxes it:")
print("magnitude x duration x volume is bounded; macroscopic static negative")
print("mass is ruled out by these bounds as far as anyone can prove.")

# ===========================================================================
header("5. EXTREME-KINEMATICS BUDGET (what 'metric engineering' would demand)")
# Hypothetical benchmark maneuver: 8.5 km displacement in ~0.78 s, stop-to-stop
d_drop, t_drop = 8534.0, 0.78
# accelerate half the distance, decelerate half:
a_req = 4*d_drop/t_drop**2
v_peak = a_req*t_drop/2
print(f"8.5 km in {t_drop} s (stop-to-stop): a ~ {sci(a_req)} m/s^2 "
      f"= {a_req/g0:,.0f} g,  v_peak ~ {v_peak/343:.0f} Mach")

m_craft = 1000.0
KE = 0.5*m_craft*v_peak**2
P_peak = m_craft*a_req*v_peak
print(f"for a {m_craft:.0f} kg craft: peak KE = {sci(KE)} J "
      f"({KE/4.184e9:.0f} tons TNT), peak power ~ {sci(P_peak)} W "
      f"({P_peak/1e9:,.0f} GW)")

# Aerodynamic heating if it actually flew through air at that speed:
M_mach = v_peak/343
T_stag = 288*(1+0.2*M_mach**2)
print(f"stagnation temperature at Mach {M_mach:.0f} in air: ~{T_stag:,.0f} K "
      f"(any material ablates; none observed)")

# Structural: force per kg at that acceleration
print(f"load on structure/occupants: {a_req/g0:,.0f} g  "
      f"(biology fails >~50 g sustained; airframes >~20 g)")
print("=> For such a maneuver, the craft cannot be ACCELERATING through air in")
print("   the Newtonian sense. Everything inside must be in free fall: the")
print("   motion must be geodesic -- i.e., spacetime around it is being shaped.")
print("   That is a statement about the METRIC, and sourcing it needs the")
print("   negative-energy budgets of section 4. This is the whole crux.")
