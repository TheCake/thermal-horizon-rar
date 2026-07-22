"""
GRAVITY AT THE LIMIT OF ITS FUNCTION.
Where does the below-a0 regime physically exist, how big, how hard to reach,
and what does the theory say the world looks like there.
"""
import math
G, c = 6.674e-11, 2.998e8
H0 = 2.27e-18
a0 = 1.2e-10
M_sun, M_e, M_m = 1.989e30, 5.972e24, 7.35e22
AU, kpc, pc = 1.496e11, 3.086e19, 3.086e16
d_moon = 3.844e8

def hdr(s): print("\n" + "="*74 + f"\n{s}\n" + "="*74)

hdr("1. WE ALREADY LIVE ON THE EDGE")
v_sun, R_sun = 233e3, 8.2*kpc
a_gal = v_sun**2/R_sun
print(f"the Sun's acceleration around the galaxy: {a_gal:.2e} m/s^2 = {a_gal/a0:.1f} a0")
print("every atom of you is, right now, in the transition zone of the anomaly.")

hdr("2. GALAXIES ARE *BUILT* AT THIS SCALE (Freeman's law)")
sigma_crit = a0/(2*math.pi*G)
sigma_freeman = 140*M_sun/pc**2
print(f"critical surface density a0/(2 pi G):        {sigma_crit:.2f} kg/m^2")
print(f"observed Freeman limit (~140 M_sun/pc^2):    {sigma_freeman:.2f} kg/m^2")
print("-> spiral discs everywhere saturate at exactly the critical density.")
print("   Galaxies aren't just AFFECTED by the a0 scale; they're SHAPED by it.")

hdr("3. THE NEAREST TRUE below-a0 BUBBLES: GRAVITATIONAL SADDLE POINTS")
# Earth-Sun saddle
r_es = AU/(1+math.sqrt(M_sun/M_e))
grad_es = 2*G*M_e/r_es**3 + 2*G*M_sun/AU**3
b_es = a0/grad_es
# Earth-Moon saddle
r_em = d_moon/(1+math.sqrt(M_e/M_m))
grad_em = 2*G*M_m/r_em**3 + 2*G*M_e/(d_moon-r_em)**3
b_em = a0/grad_em
for name, r, b, v in [("Earth-Sun ", r_es, b_es, 300.0), ("Earth-Moon", r_em, b_em, 500.0)]:
    print(f"{name} saddle: {r/1e6:7.1f} thousand km out, bubble radius ~{b:4.1f} m,")
    print(f"             crossing time at {v:.0f} m/s: {2*b/v*1000:.0f} ms "
          f"(transition halo extends ~100s of km -- that's the measurable part)")
print("-> Bekenstein & Magueijo (2006) proposed flying LISA-Pathfinder-class")
print("   accelerometers through it. LPF existed, the flyby was never funded.")
print("   A cubesat with a cold-atom gradiometer could do it. Nobody has.")

hdr("4. THE PLANET-9 ZONE (the anomaly wearing a planet costume?)")
for au in (100, 500, 1000, 7000):
    gN = G*M_sun/(au*AU)**2
    print(f"  {au:>5} AU: g_N = {gN:.1e} m/s^2 = {gN/a0:8.1f} a0   EFE-scale perturbation ~{100*min(a0/gN,1):5.2f}%")
print("-> percent-level secular perturbations at 500-1000 AU over Gyr can")
print("   shepherd TNO orbits. Several groups argue the 'Planet Nine' clustering")
print("   is EXACTLY this signature. If Vera Rubin Observatory finds no planet")
print("   but the clustering persists and aligns with the galactic external")
print("   field direction -- that's the anomaly caught red-handed in our yard.")

hdr("5. WHAT'S ALREADY CONSTRAINED (free-fall triggers)")
print("LISA Pathfinder test masses reached proper accelerations ~1e-14 m/s^2 --")
print("BELOW a0 -- and behaved perfectly conventionally. So 'modified inertia")
print("triggered by PROPER acceleration' is dead. Surviving trigger: acceleration")
print("w.r.t. the cosmic/galactic frame -- unreachable on Earth (we carry 1.8 a0,")
print("sec. 1) EXCEPT at Ignatiev's cancellation spots, where spin+orbit+galaxy")
print("momentarily null out:")
tau = 0.5e-3
print(f"  Ignatiev window ~{tau*1e3:.1f} ms -> displacement ~ a0*tau^2/2 = "
      f"{0.5*a0*tau**2:.1e} m")
print("  (his published number: 2e-17 m -- LIGO-grade, on a glacier. Untouched.)")

hdr("6. STRETCHED TO THE LIMIT: WHAT GRAVITY *BECOMES*")
print("Deep-MOND dynamics (a << a0) is invariant under (t,r) -> (L*t, L*r):")
print("SCALE-FREE. Below a0 the universe loses its acceleration ruler --")
print("dynamics stops knowing 'how big' anything is, only shapes survive.")
print("The two arrows of that statement:")
Hinf = math.sqrt(0.7)
print(f"  PAST:   a0(z) ~ H(z) was larger -> stronger effective gravity ->")
print(f"          structure forms EARLY. (Massive z~10 galaxies were predicted")
print(f"          on these grounds ~15 yr before JWST embarrassed LCDM with them.)")
print(f"  FUTURE: H -> {Hinf:.2f} H0 and freezes -> the scale-free law becomes")
print(f"          exact and eternal. The universe is mid-crystallization.")

hdr("7. BONUS STANDING ANOMALY: BIG-G METROLOGY")
print("Measured values of Newton's G from world-class labs disagree by ~5e-4 --")
print("~10x their quoted uncertainties, unresolved for 40 years. A modified-")
print("inertia paper (arXiv:1901.02604) claims the scatter correlates with")
print("apparatus acceleration regimes. Probably systematics. But note the")
print("pattern: the ONE constant we cannot pin down is gravity's.")
