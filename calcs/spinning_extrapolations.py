"""
Extrapolations: what happens if you fiddle with a spinning-superconductor rig?
- material spin limits (hoop stress), London moment and frame dragging at those limits
- the absurd upper bound: rim speed = c
- spin-up transients (gravito-Faraday induction)
- counter-rotating pair (multipole cancellation)
- gravitational-wave emission of a maximally asymmetric rotor
"""
import math

G, c = 6.674e-11, 2.998e8
m_e, e = 9.109e-31, 1.602e-19
g0 = 9.81

def sci(x, d=1): return f"{x:.{d}e}"

print("=" * 76)
print("1. HOW FAST CAN YOU ACTUALLY SPIN? (thin-ring hoop stress: v_max = sqrt(s/rho))")
print("=" * 76)
mats = [
    ("sintered YBCO ceramic", 30e6, 6300),
    ("maraging steel",        1.8e9, 8000),
    ("carbon fiber wrap",     4.0e9, 1800),
    ("CNT (theoretical)",     50e9, 1300),
]
R = 0.10   # m
M = 5.0    # kg ring
for name, sig, rho in mats:
    v = math.sqrt(sig/rho)
    om = v/R
    rpm = om*60/(2*math.pi)
    B_london = 2*m_e*om/e
    # frame dragging near the ring: B_g ~ 2 G M v / (c^2 R^2)
    Bg = 2*G*M*v/(c**2 * R**2)
    a_test = 4*300*Bg   # accel on a 300 m/s test mass
    print(f"{name:24s} v_rim={v:7.0f} m/s ({rpm:12,.0f} RPM) "
          f"B_London={B_london*1e6:8.3f} uT   B_g={sci(Bg)} 1/s  "
          f"a_test={sci(a_test)} m/s^2")
print("best detector (Eot-Wash torsion balance): ~1e-15 m/s^2")

print()
print("=" * 76)
print("2. THE ABSURD LIMIT: 1000 kg ring, R = 1 m, rim at light speed")
print("=" * 76)
for v, label in [(2.1e3, "carbon-fiber-limit rim"), (0.1*c, "0.1c rim"), (c, "c rim")]:
    Bg = 2*G*1000*v/(c**2 * 1.0**2)
    a = 4*300*Bg
    print(f"{label:24s} B_g = {sci(Bg)} 1/s   accel on 300 m/s test mass: {sci(a)} m/s^2")
print("Earth's own field for comparison:  B_g ~ 1.7e-14 1/s")
print("=> even a physically impossible light-speed rim barely reaches Earth's")
print("   natural frame-dragging strength. The knob that matters is compactness")
print(f"   GM/Rc^2: lab rig ~ {sci(G*1000/(1*c**2))}, neutron star ~ 0.2.")

print()
print("=" * 76)
print("3. SPIN-UP TRANSIENT (gravito-Faraday: E_g ~ (R/2) dB_g/dt)")
print("=" * 76)
Bg_cf = 2*G*5*1490/(c**2*0.1**2)      # 5 kg ring at carbon-fiber limit
for tau in (1.0, 1e-3):
    Eg = 0.05*Bg_cf/tau
    print(f"0 -> max spin in {tau:5.0e} s:  induced gravitoelectric accel ~ {sci(Eg)} m/s^2")
print("(this induction effect is exactly what Tajmar claimed at ~1e-3 m/s^2;")
print(" GR says it exists -- at this level, ~21 orders lower)")

print()
print("=" * 76)
print("4. COUNTER-ROTATING PAIR (coaxial, separation d, observer at distance r)")
print("=" * 76)
print("GEM is LINEAR: fields superpose, no mixing terms at lab strength.")
print("Equal + opposite spins cancel the gravitomagnetic dipole; what's left is")
print("a quadrupole suppressed by ~d/r and falling as 1/r^4 instead of 1/r^3:")
for r_obs, d in [(1.0, 0.1), (10.0, 0.1)]:
    print(f"  r={r_obs:4.0f} m, d={d} m: residual ~ {d/r_obs:.0e} of the single-disc field")
print("GR nonlinear cross-terms enter at (GM/Rc^2)^2 ~ 1e-47. Nothing to mine.")

print()
print("=" * 76)
print("5. GRAVITATIONAL-WAVE OUTPUT of a maximally asymmetric 5 kg rotor at CF limit")
print("=" * 76)
I = 0.5*5*0.1**2
om = 1490/0.1
P = (32/5)*(G/c**5)*(I**2)*(om**6)     # eps = 1
E_graviton = 1.0546e-34*2*om
print(f"radiated power: {sci(P)} W  (~{P/E_graviton:.0f} gravitons/s at f = {om/math.pi:,.0f} Hz)")
print("absorption cross-section of matter for one graviton: ~1e-68 cm^2.")
