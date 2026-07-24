"""
STAGE 6R (O5, the horizon-side round): THE RESOLUTION BATH.

The 6N cancellation-theorem corollary: WHICH occupation a self-shifting
mode takes = whether the bath's KMS structure is RESOLVED across the
self-shift. For lab baths the resolution scale is arbitrary bath
geography (6N: carrier works but is bath-geographic; sky shape
anti-standard). For the dS bath there is ONE scale: the structural
width of the KMS occupation equals the temperature itself (n_BE varies
by O(1) over dw ~ T_dS, and T_dS = H/2pi is also the a0 scale). So the
resolution parameter is FORCED:

    R = (w_tot - w_src)/T_dS = nu*y - sqrt(y)     (a0 units)

the gap between the dressed (total) and source frequencies in
bath-temperature units. Deep: R -> 0 like x^2/2 (locking EXACT - the
C&T/BE limit; explains deep+transition beta = 0 in EVERY dataset, 5T).
Tail: R -> y (fully resolved - dressing turns on). Crossover R = 1 at
y ~ 1.6: AT the MOND transition, where 5T measured the beta turn-on.
No lab bath puts the crossover there naturally; the one-scale dS bath
CANNOT put it anywhere else. This is the horizon-specific reading 6N
left standing, made quantitative.

THE FUNCTION (zero fitted parameters): the 5P frequency mixing
w = w_src^(1-beta) * w_tot^beta with RUNNING admixture

    beta(y) = (1/2) * R^2/(1+R^2),   R = nu*y - sqrt(y)
    nu = 1 + n_BE( y^((1+beta)/2) * nu^beta )        [self-consistent]

beta_max = 1/2 is the exchange-symmetric endpoint (5P); the Lorentzian
weight R^2/(1+R^2) is the simplest resolution profile and is FLAGGED as
the one representative choice (alternatives shift the crossover O(1)).

EXACT results derived + gated here:
  (i)  deep ladder: c1 = 1/2, c2 = 1/12, c3 = 0, c4 = -1/720 ALL
       PRESERVED (beta dies QUARTICALLY: beta = x^4/8 + O(x^5));
       the break is at rung FIVE: c5 = -1/16 (BE: 0). Deepest ladder
       preservation of any candidate in the program (AMB breaks at c3).
  (ii) tail: x_eff -> y^(3/4) sqrt(nu) = EXACTLY the gm argument;
       tail exponent p = 3/4.
 (iii) the exact-temperature matrix (Deser-Levin T = sqrt(T_dS^2+T_U^2))
       is CLOSED analytically: source-frequency assignment gives an
       anti-Newtonian catastrophe (nu-1 ~ sqrt(y)/2pi grows); total-
       frequency gives boot deep (dead, 5F/5M) + an invisible constant
       floor n_BE(2pi) = pure G-renormalization; a free-fall detector
       (geodesic: zero proper acceleration) gives T_U = 0 = pure C&T,
       i.e. the binary beta < 0.03 is the free-fall answer.

PRE-REGISTERED BARS (committed before any fit; tests = 6S galaxies,
6T binaries, existing machinery unchanged):
  6S (hier galaxy ladder vs BE): expected gm-grade-or-better (deep = BE
      exact, tail = gm exact). PASS if vertical Delta <= -40; STRONG if
      <= -55 (F4/AMB grade); STRIKE against the resolution reading if
      > -20; between = partial.
  6T (v7 binary budget vs p050, corrected velocities, 6 seeds): ACCEPT
      if mean Delta >= -3 (AMB grade -0.88+-2.66); REJECT if <= -5
      (joins the eight sharpened functions); between = unresolved.
      Interior alpha-hat required for a valid read; edge-ride = shape
      rejection. Temperature row: PASS if a0 pull <= +2.5 sigma.
  HONEST PRIOR RISK, stated now: nu_R(1) ~ 1.54 (F4-grade transition)
      and beta(y=1.2) ~ 0.15 with NO ambient gate, so the prior from
      the eight-function pattern leans REJECT on binaries. The test
      DISCRIMINATES the two surviving carriers of the two-system split:
      REJECT => the ambient gate is REQUIRED (AMB's system-level story
      stands alone); ACCEPT => a pre-hoc-derived rival to AMB with
      cleaner provenance (this file precedes the fits in git).

Writes data/stage6r_resolution.txt.
"""
import math
import numpy as np
import sympy as sp
import mpmath as mp

L = []
def say(s=''):
    L.append(s); print(s, flush=True)

say("STAGE 6R: the resolution bath -- derivation gates")
say("=" * 72)

# ---------- G1-sym: exact deep series (sympy, rational arithmetic) ----------
x = sp.symbols('x', positive=True)
ORD = 8
nser = sp.series(1/(sp.exp(x) - 1), x, 0, ORD).removeO()      # n_BE(x)
phi0 = sp.expand(x + x*nser)                                   # nu*sqrt(y) at BE
R = sp.expand(x*(phi0 - 1))
w = sp.series(R**2/(1 + R**2), x, 0, ORD).removeO()
beta = w/2
dlt = sp.series(beta*sp.log(phi0), x, 0, ORD).removeO()        # ln(x_eff/x)
xeff = sp.series(x*sp.exp(dlt), x, 0, ORD).removeO()
z = sp.symbols('z', positive=True)
Nz = sp.series(1/(sp.exp(z) - 1), z, 0, ORD).removeO()
phiR = sp.series(x*(1 + Nz.subs(z, xeff)), x, 0, ORD).removeO()
diff = sp.expand(phiR - phi0)
c_break = [sp.nsimplify(diff.coeff(x, k)) for k in range(ORD)]
say("G1-sym: phi_R - phi_BE series coefficients (rungs 0..7):")
say("        " + ", ".join(str(c) for c in c_break))
ok1 = all(c_break[k] == 0 for k in range(5)) and c_break[5] == sp.Rational(-1, 16)
say(f"G1-sym (rungs 1-4 preserved; break c5 = -1/16 exact): "
    f"{'PASS' if ok1 else 'FAIL'}")
assert ok1
phi_lad = [sp.nsimplify(sp.expand(phiR).coeff(x, k)) for k in range(7)]
say(f"        full nu_R ladder: c1={phi_lad[1]}, c2={phi_lad[2]}, "
    f"c3={phi_lad[3]}, c4={phi_lad[4]}, c5={phi_lad[5]}")
say(f"        (BE Bernoulli:    c1=1/2, c2=1/12, c3=0, c4=-1/720, c5=0)")

# ---------- the numeric fixed-point solver (mpmath reference) ----------
mp.mp.dps = 50
def nuR_mp(y):
    y = mp.mpf(y); sy = mp.sqrt(y)
    def F(nu):
        Rv = nu*y - sy
        b = mp.mpf('0.5')*Rv*Rv/(1 + Rv*Rv)
        xe = sy**(1 - b) * (nu*y)**b
        return nu - 1 - 1/mp.expm1(xe)
    lo = 1 + 1/mp.expm1(sy)                       # BE seed
    return mp.findroot(F, lo)

# G1-num: the break rung numerically
say("")
vals = []
for xv in ('1e-2', '5e-3'):
    xm = mp.mpf(xv); ym = xm*xm
    phiRn = xm*nuR_mp(ym)
    phiBn = xm*(1 + 1/mp.expm1(xm))
    vals.append(float((phiRn - phiBn)/xm**5))
say(f"G1-num: (phi_R - phi_BE)/x^5 at x=1e-2, 5e-3: "
    f"{vals[0]:+.6f}, {vals[1]:+.6f}  (exact -1/16 = -0.0625)")
ok1n = abs(vals[1] + 1.0/16) < 2e-3
say(f"G1-num: {'PASS' if ok1n else 'FAIL'}")
assert ok1n

# ---------- G2: tail identity (the gm argument) ----------
say("")
tail = []
for yv in (1e4, 1e6):
    nu = nuR_mp(yv)
    b = mp.mpf('0.5')*((nu*yv - mp.sqrt(yv))**2)/(1 + (nu*yv - mp.sqrt(yv))**2)
    xe = mp.sqrt(yv)**(1 - b)*(nu*yv)**b
    tail.append(float(xe/(mp.mpf(yv)**mp.mpf('0.75')*mp.sqrt(nu)) - 1))
say(f"G2: x_eff/(y^0.75 sqrt(nu)) - 1 at y=1e4, 1e6: "
    f"{tail[0]:+.2e}, {tail[1]:+.2e}")
ok2 = abs(tail[1]) < 1e-4
say(f"G2 (tail = gm argument, p = 3/4): {'PASS' if ok2 else 'FAIL'}")
assert ok2

# ---------- the production numpy solver + gates ----------
def nu_simple(yy):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(yy, 1e-14, None))

def make_resn(wfac=1.0):
    def nu_run(yy):
        yy = np.clip(np.asarray(yy, float), 1e-14, None)
        sy = np.sqrt(yy); ly = np.log(yy)
        nu = nu_simple(yy)
        for _ in range(120):
            Rv = nu*yy - sy
            R2 = Rv*Rv
            wv = wfac*R2/(1.0 + R2)
            b = 0.5*wv
            dw = wfac*2.0*Rv/((1.0 + R2)**2)
            db = 0.5*dw*yy
            lnu = np.log(nu)
            u = np.exp(np.minimum(0.5*(1.0 + b)*ly + b*lnu, 60.0))
            eu = np.exp(np.minimum(u, 60.0))
            em1 = np.maximum(eu - 1.0, 1e-300)
            n = np.where(u < 60.0, 1.0/em1, 0.0)
            F = nu - 1.0 - n
            dudnu = u*(db*(0.5*ly + lnu) + b/nu)
            dF = 1.0 + (eu/(em1*em1))*dudnu
            step = F/np.where(np.abs(dF) > 1e-12, dF, 1e-12)
            nu = np.maximum(nu - np.clip(step, -0.5*nu, 0.5*nu), 1.0 + 1e-15)
        return nu
    return nu_run

nu_resn = make_resn(1.0)
say("")
yg = np.logspace(-6, 6, 241)
nug = nu_resn(yg)
sy = np.sqrt(yg); Rg = nug*yg - sy
bg = 0.5*Rg*Rg/(1 + Rg*Rg)
ue = sy**(1 - bg)*(nug*yg)**bg
res = np.abs(nug - 1 - np.where(ue < 60, 1.0/np.expm1(np.minimum(ue, 60)), 0.0))
say(f"G3: production-solver residual |nu-1-n(x_eff)| max over y in "
    f"[1e-6,1e6]: {res.max():.2e}  {'PASS' if res.max() < 1e-9 else 'FAIL'}")
assert res.max() < 1e-9
spot = [1e-4, 1e-2, 0.3, 1.0, 1.2, 3.0, 30.0]
dmp = max(abs(float(nuR_mp(yv)) - float(nu_resn(np.array([yv]))[0]))
          for yv in spot)
say(f"G3b: numpy vs mpmath fixed point (7 spots): max |d| = {dmp:.2e}  "
    f"{'PASS' if dmp < 1e-8 else 'FAIL'}")
assert dmp < 1e-8
nu_be_chk = 1.0/(1.0 - np.exp(-np.sqrt(yg)))
d0 = np.max(np.abs(make_resn(0.0)(yg)/nu_be_chk - 1.0))
say(f"G4: wfac=0 regression to pure BE: max rel {d0:.2e}  "
    f"{'PASS' if d0 < 1e-10 else 'FAIL'}")
assert d0 < 1e-10

# ---------- the facts table ----------
say("")
say("FACTS:")
def bat(yv):
    nu = float(nu_resn(np.array([yv]))[0])
    Rv = nu*yv - math.sqrt(yv)
    return nu, Rv, 0.5*Rv*Rv/(1 + Rv*Rv)
nu1, R1, b1 = bat(1.0)
say(f"  nu_R(1) = {nu1:.4f}  [BE 1.582, gm 1.433, F4 1.537, AMB(bin) 1.577]")
ystar = None
for yv in np.linspace(0.5, 4.0, 3501):
    _, Rv, _ = bat(yv)
    if Rv >= 1.0:
        ystar = yv; break
say(f"  resolution crossover R=1 at y* = {ystar:.2f}  (the 5T beta turn-on "
    f"arm is y > 1; binaries sit at y_tot ~ 1.2)")
for yv in (0.1, 0.3, 1.0, 1.2, 1.9, 3.0, 10.0):
    nu, Rv, bb = bat(yv)
    say(f"  y={yv:5.2f}: nu_R={nu:7.4f}  R={Rv:7.3f}  beta={bb:.3f}  "
        f"(nu_BE={1.0/(1.0-math.exp(-math.sqrt(yv))):7.4f})")

# ---------- (iii) the exact-temperature matrix, closed ----------
say("")
say("THE DESER-LEVIN MATRIX (T_eff = sqrt(T_dS^2 + T_U^2), exact closures):")
yv = 100.0
nu = 1.0
runaway = False
for _ in range(200):
    xe = math.sqrt(yv)/math.sqrt(1.0 + (nu*yv/(2*math.pi))**2)
    nu = 1.0 + 1.0/math.expm1(xe)
    if nu > 1e6:
        runaway = True; break
say(f"  source-frequency branch: NO finite fixed point at y=100 "
    f"(runaway confirmed: {runaway}); small-occupation branch nu-1 ~ "
    f"sqrt(y)/2pi = {math.sqrt(yv)/(2*math.pi):.2f} and growing -- "
    f"anti-Newtonian catastrophe: EXCLUDED analytically")
nfloor = 1.0/math.expm1(2*math.pi)
say(f"  total-frequency branch: deep -> boot (DEAD, 5F/5M); tail -> "
    f"constant floor n_BE(2pi) = {nfloor:.5f} = pure G-renormalization "
    f"(invisible): EXCLUDED by inheritance")
say(f"  free-fall detector (geodesic, zero proper acceleration): T_U = 0 "
    f"exactly -> pure C&T (beta = 0) -- the binary beta < 0.03 (5R) IS "
    f"the free-fall answer; the resolution channel is what remains")

# ---------- pre-registered bars ----------
say("")
say("PRE-REGISTERED BARS (this file precedes 6S/6T execution in git):")
say("  6S galaxies (hier vs BE): PASS <= -40 vertical; STRONG <= -55; "
    "STRIKE > -20; else partial.")
say("  6T binaries (vs p050, 6 seeds): ACCEPT >= -3 mean; REJECT <= -5; "
    "else unresolved. Interior alpha-hat required; edge = shape rejection. "
    "a0 pull PASS <= +2.5 sigma.")
nu12, _, b12 = bat(1.2)
say("  PRIOR RISK stated: nu_R(1) = %.3f (F4-grade), beta(1.2) = %.2f, "
    "no ambient gate -> prior leans binary-REJECT; the test discriminates "
    "ambient-gating vs pure-resolution as the split's carrier." % (nu1, b12))

with open('data/stage6r_resolution.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nSTAGE 6R done -> data/stage6r_resolution.txt")
