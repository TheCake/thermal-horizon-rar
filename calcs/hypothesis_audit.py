"""
THE JACOBIAN MOVE, APPLIED TO GRAVITY.
Don't attack the impossibility theorems -- attack their hypotheses.
Every number here is computed from constants; no lookups.

Sections:
  1. The audit: what each no-go theorem actually assumes
  2. The existence proof: repulsive gravity is already real (dark energy)
  3. The engineering gap: vacuum energy as a dial (the 10^122 hole)
  4. The graviton: what spin-2 forces on us, and what remains untested
  5. The residual slivers: graviphoton, the MOND coincidence
"""
import math

G     = 6.674e-11
c     = 2.998e8
hbar  = 1.0546e-34
H0    = 2.27e-18        # ~70 km/s/Mpc
Om_L  = 0.7
g0    = 9.81
GeV   = 1.602e-10       # J
mu0   = 4e-7*math.pi
M_earth = 5.97e24
eV    = 1.602e-19

def sci(x, d=1): return f"{x:.{d}e}"
def hdr(s): print("\n" + "="*76 + f"\n{s}\n" + "="*76)

hdr("1. HYPOTHESIS AUDIT OF THE NO-GO THEOREMS")
rows = [
 ("Positive mass theorem",   "classical matter; DOMINANT ENERGY CONDITION; 3+1D",
  "DEC is violated by real QFT (Casimir). Survives only b/c violations are tiny."),
 ("No-shielding argument",   "gravity = metric theory; one charge sign",
  "airtight IF gravity is pure spin-2 metric. Hypothesis: no second field."),
 ("ANEC 'theorem'",          "FLAT spacetime; proven for free/holographic QFT",
  "NOT proven for interacting fields in curved spacetime. Genuine gap."),
 ("Ford-Roman inequalities", "FREE fields, inertial worldlines",
  "no general interacting-4D proof exists (only 2D CFT). Genuine gap."),
 ("Warp NEC violation (SSV)","the warp ansatz is sourced by MATTER stress-energy",
  "if vacuum energy itself is a dial, 'matter' bookkeeping changes."),
 ("Universal attraction",    "mediator = massless spin-2 coupled to conserved T",
  "THEOREM-level (Weinberg soft limit). But repulsion via rho+3p<0 is ALLOWED."),
]
for name, hyp, gap in rows:
    print(f"* {name}\n    assumes: {hyp}\n    crack:   {gap}")

hdr("2. EXISTENCE PROOF: REPULSIVE GRAVITY IS ALREADY REAL")
# Dark energy: p = -rho. Gravitating mass density is rho+3p = -2rho < 0.
rho_L = Om_L * 3*H0**2/(8*math.pi*G)          # kg/m^3
print(f"dark energy density: {sci(rho_L)} kg/m^3 = {sci(rho_L*c**2)} J/m^3")
print("equation of state p = -rho  ->  rho + 3p = -2 rho < 0  ->  g_tt repels.")
g_L = Om_L * H0**2 * 1.0                       # (8piG/3)rho_L * r at r=1m
print(f"repulsive acceleration it produces at r = 1 m: {sci(g_L)} m/s^2")
r_win = (G*M_earth/(Om_L*H0**2))**(1/3)
print(f"distance from Earth where dark-energy repulsion WINS over Earth's pull:")
print(f"  r = (GM/(Om_L H0^2))^(1/3) = {sci(r_win)} m = {r_win/9.46e15:.1f} light-years")
print("Repulsive gravity is not exotic. It runs the universe. It is just DILUTE.")

hdr("3. THE ENGINEERING GAP: VACUUM ENERGY AS A DIAL")
rho_need = 3*g0/(8*math.pi*G*1.0)              # 1 g repulsion at r = 1 m
print(f"vacuum-energy density needed for 1 g repulsion at 1 m (de Sitter bubble):")
print(f"  rho = 3g/(8 pi G r) = {sci(rho_need)} kg/m^3 = {sci(rho_need*c**2)} J/m^3")
# known condensates of the Standard Model vacuum:
def E4_density(E_GeV):                          # (E)^4 / (hbar c)^3
    return (E_GeV*GeV)**4 / (hbar*c)**3
rho_qcd = E4_density(0.2)
rho_ew  = E4_density(246.0)
rho_pl  = c**7/(hbar*G**2)                      # Planck density (J/m^3 via c^2 factor folded)
print(f"QCD vacuum condensate scale  (0.2 GeV)^4:  {sci(rho_qcd)} J/m^3")
print(f"electroweak condensate scale (246 GeV)^4:  {sci(rho_ew)} J/m^3")
print(f"naive QFT zero-point (Planck cutoff):      {sci(rho_pl)} J/m^3")
print(f"observed vacuum energy:                    {sci(rho_L*c**2)} J/m^3")
print(f"mismatch (the cosmological constant problem): 10^{math.log10(rho_pl/(rho_L*c**2)):.0f}")
print(f"ratio (QCD condensate)/(needed for 1g@1m):    {sci(rho_qcd/(rho_need*c**2),0)}")
print("=> The Standard Model vacuum ALREADY stores ~1e7 times more energy density")
print("   than engineered gravitational repulsion would require. We observe that it does NOT")
print("   gravitate accordingly -- and NOBODY KNOWS WHY. That ignorance is 122")
print("   orders of magnitude deep. This is the hole in the impossibility proof.")
# and the superconducting cameo:
cond_Nb = 0.2**2/(2*mu0)                       # B_c^2/2mu0, Nb thermodynamic critical field
print(f"\nsuperconducting condensation energy (Nb): {sci(cond_Nb)} J/m^3")
print(f"  = a REAL, human-switchable vacuum-energy change; fractional weight")
print(f"  effect if it gravitates normally: {sci(cond_Nb/(6300*c**2))} -- what")
print(f"  Archimedes is built to weigh. First-ever direct test of the dial.")

hdr("4. THE GRAVITON: WHAT IS FORCED, WHAT IS OPEN")
print("Forced (Lorentz inv. + massless + long-range + couples to energy):")
print("  spin 0: attracts, but no light bending -> dead (light DOES bend, 2x Newton)")
print("  spin 1: like charges REPEL -> matter-matter would repel -> dead")
print("  spin 2: universal attraction, bends light correctly -> the only survivor;")
print("  Weinberg soft-theorem: it MUST couple to everything identically (EP).")
print("  => selective gravitational shielding is spin-theoretically dead;")
print("     ONLY route: source with rho+3p<0 (sec. 2/3). Consistent story.")
m_g_bound = 1.2e-22*eV/c**2                    # GW dispersion bound (order of magnitude)
lam_g = hbar/(m_g_bound*c)
print(f"graviton mass bound ~1e-22 eV -> Compton wavelength {sci(lam_g)} m "
      f"(~{lam_g/9.46e15:.2f} ly)")
# is gravity quantum at all? The BMV phase between two mesoscopic masses:
m_bmv, d_bmv = 1e-14, 250e-6
phase_rate = G*m_bmv**2/(hbar*d_bmv)
print(f"entangling phase rate, two 1e-14 kg masses at 250 um: "
      f"{phase_rate:.2f} rad/s -> O(1) in seconds.")
print("Whether that phase exists (gravity quantum?) is UNTESTED. 2602.12266 lives here.")

hdr("5. RESIDUAL SLIVERS (allowed by everything measured)")
print(f"* B-L 'graviphoton' (spin-1, repels matter-matter): allowed below ~1e-9")
print(f"  of gravity's strength at long range (EP tests). Max levitation payoff:")
print(f"  cancels ~{sci(1e-9*g0)} m/s^2 of your weight. A real repulsive force,")
print(f"  legally permitted, utterly useless.")
a_mond = c*H0/(2*math.pi)
print(f"* the acceleration coincidence: BOTH empirical cracks in GR (rotation")
print(f"  curves, dark energy) live at the same scale a0 ~ 1.2e-10 m/s^2, and")
print(f"  c*H0/(2*pi) = {sci(a_mond)} m/s^2. Coincidence or the biggest clue we")
print(f"  have. No lab has ever tested gravity's force law AT this acceleration")
print(f"  in clean conditions -- the regime is masked by Earth's 9.8 m/s^2.")
