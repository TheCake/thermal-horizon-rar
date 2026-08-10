"""
STAGE 10N (2026-08-10): O5-CROSSING -- the kappa = 1 level-crossing
low-frequency response: derive the NUMERIC KILL CONDITION (R30
successor 2, R31: "the one place a positive kappa ~ 1 signal could
still emerge, and it currently carries no numeric kill condition";
the standing rule: deriving it is the PREREQUISITE for any
PREDICTIONS row).  DISPERSIVE-SIDE ONLY (10G strike scope).  NO SKY
FITS: every number below is an archived meter, an archived-stage
token, or a derived quantity of the model.

THE OBJECT.  The dressed ladder E(n) = om(n+1/2) - (kappa om/4)
(n+1/2)^2 has the pair gap D1 = E(2)-E(1) = om(1-kappa) EXACT: at
kappa -> 1 the states |1>,|2> become degenerate and hybridize under
any static Delta-n = +-1 perturbation (R30).  The dispersive vertex
itself CANNOT mix them ([H, N_s] = 0, 10B); the only archived mixing
agent is the l=2 near-field vertex (10E, structure grade), whose
sub-saturation amplitude is licensed -- the 10G strike killed ONLY
the saturated (lambda >~ H) real-exchange lending leg; here lambda
enters as a perturbative drive bounded by the archived FD/
participation ceilings, and the lending law is used NOWHERE.

ARCHIVED ANCHORS (provenance in comments; gates re-derive tokens):
  n_amb(gal) = 6.583, n_amb(bin) = 0.5202; gates 0.7536/0.1171
    [data/stage6e_ambgate.txt line 3]
  lambda_max/g_close (e_a = 1): binary x=0.5: 0.304, x=1.0: 0.140,
    x=2.0: 0.029; deep-galaxy x=0.0953: 0.874; g_close =
    (1/2)sqrt(om Om); galaxy g_close = 0.0092 H at x=0.0953
    [data/stage10de_addendum.txt GA-2/GA-3]
  static EFE amplitude |q| = 0.086 (the DC quadrupole; the
    FLUCTUATION amplitude e_a <= ~1 is the 10F open object)
  soft-sector width band gamma/H in [0.010, 0.025] (R17
    separability rates: decoherent <= 0.025, finite-time <= ~0.01;
    NOTES-horizon-inertia.md:13192); scan family [0.005, 0.1]
  kappa rows (archived meters, 10J G4): lock 0.925 (Planck) / 0.888
    (SH0ES); c1-meters 1.100 (flat) / 1.484 (hier)
  scales: om = x_loc H/(2 pi); Om = x_amb H/(2 pi), x_amb =
    ln(1 + 1/n_amb)  [the 6E convention; g_close gate below]

THE RESPONSE (derived, gated): thermal pull modification from the
pair channel, S = (p1 - p2) * delta_gamma / |W|, with p_n the BARE
thermal weights (10L bare-pinning), |W| = (kappa om/4)(2 nbar + 2)
the dictionary pull, and delta_gamma the width-regularized exact
2x2 level repulsion, delta = Re[sqrt(z^2 + 2 U^2 P0) - z],
z = |D1|/2 - i gamma/2, U = e_a lambda_max (convention family
{1/sqrt2, 1, sqrt2} disclosed), P0 = e^{-d^2} the Franck-Condon
diagonal, d = g/Om.  Perturbative limit delta -> 2 U^2 P0 D1/
(D1^2 + gamma^2) (Lorentzian; peak at 1 - kappa = gamma/om).
Window function (p1-p2)(x) = e^{-x}(1-e^{-x})^2, max 4/27 at
x = ln 3 = 1.0986 -- the TRANSITION window.

GATES (bars locked at this commit):
  G10N-1 sympy, zero residue: D1 = om(1-kappa); 2x2 eigenvalues
         D/2 +- sqrt((D/2)^2+V^2); perturbative series delta =
         V^2/D + O(V^4); Lorentzian peak at u = gamma, max 1/(2
         gamma).
  G10N-2 reader identity vs 10J: the chi1 machinery (verbatim form,
         Om = 0.8 om) reproduces the printed table {10.365, 30.226,
         293.478, 2927.284} at kappa {0.5, 0.9, 0.99, 0.999},
         rel <= 5e-4.
  G10N-3 sideband census over kappa in [0.5, 1.05]: binary anchors
         and deep-galaxy anchor CLEAN (min_m |D0 - m Om|, ||D1| -
         m Om| >= 0.1 om); galaxy-transition anchor honestly mapped
         (expected DISQUALIFIED: sideband inside the width band AND
         comb spacing Om ~ gamma -- print the numbers).
  G10N-4 scale wiring: (n_amb/(1+n_amb))^2 reproduces the 6E gate
         tokens <= 1e-3; g_close(deep-gal) reproduces 0.0092 H
         <= 2%; the lambda_max table printed in H units.
  G10N-5 surface checks: e_a -> 0 log-slope = 2.000 +- 1e-3
         (Richardson two-point at e_a = 1e-3/2e-3); perturbative
         peak location |(1 - kappa*) om - gamma| <= 0.006 at
         e_a = 0.086, gamma = 0.015, x = 1.0.
  G10N-6 the kill table fires the letter MECHANICALLY (grammar
         below) and the operative clause of WHICHEVER letter has a
         numeric check (trap #12: gate every positive clause):
         EARNED -> |S(e_a*) - 0.02| <= 1e-3 at the quoted corner;
         BELOW-FLOOR -> revival bar in [0.002, 0.05] and S
         monotone in e_a on the grid; RETIRED -> best S(e_a=1)
         < 0.002 across the anchored grid.
  G10N-7 scope audit printed: dispersive-side licensing, no
         lending-law input, no sky fits, convention family spread.

VERDICT LETTERS (pre-signed; anchored grid = x in {0.5, 1.0} binary,
gamma in {0.010, 0.015, 0.025}, kappa in {0.888, 0.925, 1.000}):
  N-ROW-EARNED  some anchored corner has e_a*(sigma_S = 0.02) <= 1:
                the kill condition is derived AND a DR4-grade window
                exists inside the permitted e_a band -> register the
                PREDICTIONS row (conditional triple: kappa near
                lock/1 x e_a >= e_a* x dispersive reading) in the
                books commit.
  N-BELOW-FLOOR kill condition derived; best S(e_a = 1) in [0.002,
                0.02) -> NO row (below the pre-set DR4 floor); the
                surface + revival bar sigma_S^rev = best S are the
                standing design spec.
  N-RETIRED     best S(e_a = 1) < 0.002 everywhere anchored: the
                reading retires; no row.
  N-SCOPE-BLOCKED the only mixing agent fails licensing, or the
                two-level scope fails structurally (sideband gates).
  N-AMBIG       gates fail / family-unstable.

CREDENCE MAP (pre-signed): ALL cells -> HOLD mech-conditional 8
(an instrument-derivation stage moves nothing); anomaly-real 53
UNTOUCHED in all cells (no sky fits).

COMPUTE < 2 min.  Writes data/stage10n_crossing.txt.
"""
import math
import numpy as np
import sympy as sp

OUT = "data/stage10n_crossing.txt"
L = []
def emit(s=""):
    print(s)
    L.append(s)

emit("STAGE 10N O5-CROSSING (pre-reg: this commit)")
emit("")

TWO_PI = 2*math.pi
n_amb_gal, n_amb_bin = 6.583, 0.5202
g_tok_gal, g_tok_bin = 0.7536, 0.1171
x_amb_gal = math.log(1 + 1/n_amb_gal)
x_amb_bin = math.log(1 + 1/n_amb_bin)
lam_ratio_bin = {0.5: 0.304, 1.0: 0.140, 2.0: 0.029}
lam_ratio_gal_deep = 0.874
E_A_STATIC = 0.086
GAM_GRID = (0.010, 0.015, 0.025)
KAP_GRID = (0.888, 0.925, 1.000)
X_BIN = (0.5, 1.0)

def om_of(x):
    return x/TWO_PI

def nbe(x):
    return 1.0/(math.exp(x) - 1.0)

def pcontrast(x):
    e = math.exp(-x)
    return e*(1 - e)**2

def g_close(x_loc, x_amb):
    return 0.5*math.sqrt(om_of(x_loc)*om_of(x_amb))

def W_dict(x_loc, kap):
    return (kap*om_of(x_loc)/4)*(2*nbe(x_loc) + 2)

def fc0(x_loc, x_amb, kap):
    om, Om = om_of(x_loc), om_of(x_amb)
    d2 = (kap/4)*(om/Om)
    return math.exp(-d2)

def delta_pair(x_loc, x_amb, kap, gam, U, conv=1.0):
    om = om_of(x_loc)
    D1 = abs(om*(1 - kap))
    Vp2 = 2*(conv*U)**2*fc0(x_loc, x_amb, kap)
    z = complex(D1/2, -gam/2)
    w = np.sqrt(z*z + Vp2)
    return float(np.real(w) - np.real(z))

def S_stat(x_loc, x_amb, kap, gam, e_a, lam_max, conv=1.0):
    U = e_a*lam_max
    return (pcontrast(x_loc)*delta_pair(x_loc, x_amb, kap, gam, U,
                                        conv)/W_dict(x_loc, kap))

# ================= G10N-1: sympy =================
omS, kapS, VS, DS, uS, gamS, nS = sp.symbols(
    'omega kappa V D u gamma n', positive=True)
sig = sp.Rational(1, 2)
E = omS*(nS + sig) - (kapS*omS/4)*(nS + sig)**2
r1 = sp.simplify((E.subs(nS, 2) - E.subs(nS, 1)) - omS*(1 - kapS))
M2 = sp.Matrix([[0, VS], [VS, DS]])
evs = M2.eigenvals()
r2 = sp.simplify(sum(evs.keys()) - DS)
delta_sym = sp.sqrt((DS/2)**2 + VS**2) - DS/2
r3 = sp.simplify(sp.series(delta_sym, VS, 0, 4).removeO()
                 - VS**2/DS)
lor = uS/(uS**2 + gamS**2)
ustar = sp.solve(sp.diff(lor, uS), uS)
r4a = [u for u in ustar if u.is_positive]
r4 = sp.simplify(r4a[0] - gamS) if r4a else sp.Integer(1)
r5 = sp.simplify(lor.subs(uS, gamS) - 1/(2*gamS))
g1 = (r1 == 0 and r2 == 0 and r3 == 0 and r4 == 0 and r5 == 0)
emit("G10N-1 sympy: D1 = om(1-kappa) residue %s; 2x2 trace residue %s;"
     % (r1, r2))
emit("   perturbative delta = V^2/D residue %s; Lorentzian peak at"
     % r3)
emit("   u = gamma (residue %s), max 1/(2 gamma) (residue %s) -> %s"
     % (r4, r5, "PASS" if g1 else "FAIL"))
emit("")

# ================= G10N-2: reader identity vs 10J =================
def chi1_10j(kv, Mq=240):
    omv, Omv = 1.0, 0.8
    gv = math.sqrt(kv*Omv*omv/4)
    d = gv/Omv
    D1 = omv*(1 - kv)
    D0m = omv*(1 - kv/2)
    P = [math.exp(-d*d)*d**(2*m)/math.factorial(m) for m in range(60)]
    up = 2*2*sum(P[m]/(D1 + m*Omv) for m in range(60))
    dn = 2*1*sum(P[m]/(-D0m + m*Omv) for m in range(60)
                 if abs(-D0m + m*Omv) > 1e-12)
    return up + dn

PRINTED = {0.5: 10.365, 0.9: 30.226, 0.99: 293.478, 0.999: 2927.284}
worst2 = 0.0
for kv, tok in PRINTED.items():
    c = chi1_10j(kv)
    worst2 = max(worst2, abs(c - tok)/tok)
g2 = worst2 <= 5e-4
emit("G10N-2 reader identity vs the 10J chi(1) table: worst rel d = "
     "%.1e (bar 5e-4) -> %s" % (worst2, "PASS" if g2 else "FAIL"))
emit("")

# ================= G10N-3: sideband census =================
emit("G10N-3 sideband census (min over m >= 1, kappa in [0.5, 1.05]):")
def census(x_loc, x_amb):
    om, Om = om_of(x_loc), om_of(x_amb)
    gap_min = 1e9
    for kap in np.linspace(0.5, 1.05, 56):
        D0 = om*(1 - kap/2)
        D1 = abs(om*(1 - kap))
        for m in range(1, 9):
            gap_min = min(gap_min, abs(D0 - m*Om)/om,
                          abs(D1 - m*Om)/om)
    return gap_min

rows = [("binary x=0.5", 0.5, x_amb_bin),
        ("binary x=1.0", 1.0, x_amb_bin),
        ("gal-deep x=0.0953", 0.0953, x_amb_gal)]
g3 = True
for name, xl, xa in rows:
    gm = census(xl, xa)
    ok = gm >= 0.1
    g3 = g3 and ok
    emit("   %-18s Om/om = %.4f  min sideband gap = %.3f om (bar "
         "0.1) -> %s" % (name, om_of(xa)/om_of(xl), gm,
                         "PASS" if ok else "FAIL"))
# galaxy-transition: honest disqualification map
xl, xa = 1.0, x_amb_gal
om, Om = om_of(xl), om_of(xa)
gap_gt = census(xl, xa)
comb_ratio_lo = GAM_GRID[0]/Om
comb_ratio_hi = GAM_GRID[-1]/Om
emit("   gal-transition x=1.0: Om/om = %.4f, min sideband gap = "
     "%.4f om = %.4f H" % (Om/om, gap_gt, gap_gt*om))
emit("     vs width band [%.3f, %.3f] H: sideband INSIDE the width; "
     % (GAM_GRID[0], GAM_GRID[-1]))
emit("     comb spacing Om = %.4f H vs gamma: gamma/Om = %.2f-%.2f "
     "-> comb UNRESOLVED" % (Om, comb_ratio_lo, comb_ratio_hi))
gt_disq = (gap_gt*om < GAM_GRID[-1]) and (comb_ratio_hi >= 1.0)
emit("     => galaxy-transition channel DISQUALIFIED (%s)"
     % ("as expected" if gt_disq else "UNEXPECTEDLY CLEAN -- flag"))
# deep-galaxy contrast kill
ratio_deep = pcontrast(0.0953)/(4/27)
emit("   gal-deep contrast: (p1-p2)(0.0953) = %.5f = %.3f of the "
     "window max 4/27" % (pcontrast(0.0953), ratio_deep))
emit("     => deep-galaxy channel suppressed x%.0f by occupation "
     "contrast (classical limit:" % (1/ratio_deep))
emit("        adjacent levels near-equally populated); the crossing "
     "meter is BINARY-NATIVE,")
emit("        window (p1-p2)(x) = e^-x (1-e^-x)^2 max at x = ln 3 = "
     "%.4f (4/27 = %.4f)" % (math.log(3), 4/27))
emit("G10N-3: %s" % ("PASS" if (g3 and gt_disq) else "FAIL"))
emit("")

# ================= G10N-4: scale wiring =================
d_gal = abs((n_amb_gal/(1 + n_amb_gal))**2 - g_tok_gal)
d_bin = abs((n_amb_bin/(1 + n_amb_bin))**2 - g_tok_bin)
gc_gal = g_close(0.0953, x_amb_gal)
d_gc = abs(gc_gal - 0.0092)/0.0092
g4 = d_gal <= 1e-3 and d_bin <= 1e-3 and d_gc <= 0.02
emit("G10N-4 scale wiring: 6E gate tokens reproduce (d = %.1e/%.1e, "
     "bar 1e-3);" % (d_gal, d_bin))
emit("   g_close(deep-gal) = %.5f H vs archived 0.0092 (rel %.3f, "
     "bar 0.02) -> %s" % (gc_gal, d_gc, "PASS" if g4 else "FAIL"))
emit("   the lambda_max table (e_a = 1, H units):")
LAM = {}
for xl in X_BIN:
    gcb = g_close(xl, x_amb_bin)
    LAM[xl] = lam_ratio_bin[xl]*gcb
    emit("     binary x=%.1f: g_close = %.5f H, lambda_max = %.5f H"
         % (xl, gcb, LAM[xl]))
emit("     (deep-gal x=0.0953: lambda_max = %.5f H -- channel "
     "contrast-killed above)" % (lam_ratio_gal_deep*gc_gal))
emit("")

# ================= G10N-5: surface checks =================
s1 = S_stat(1.0, x_amb_bin, 0.925, 0.015, 1e-3, LAM[1.0])
s2 = S_stat(1.0, x_amb_bin, 0.925, 0.015, 2e-3, LAM[1.0])
slope = math.log(s2/s1)/math.log(2.0)
g5a = abs(slope - 2.0) <= 1e-3
kfine = np.linspace(0.80, 0.999, 400)
svals = [S_stat(1.0, x_amb_bin, k, 0.015, E_A_STATIC, LAM[1.0])
         for k in kfine]
kstar = float(kfine[int(np.argmax(svals))])
pk_res = abs((1 - kstar)*om_of(1.0) - 0.015)
g5b = pk_res <= 0.006
g5 = g5a and g5b
emit("G10N-5 surface: e_a log-slope = %.5f (bar 2 +- 1e-3) -> %s;"
     % (slope, "PASS" if g5a else "FAIL"))
emit("   perturbative peak kappa* = %.4f, |(1-k*) om - gamma| = "
     "%.4f (bar 0.006) -> %s" % (kstar, pk_res,
                                 "PASS" if g5b else "FAIL"))
emit("")

# ================= THE KILL TABLE =================
emit("THE KILL TABLE: S = fractional transition-window pull "
     "modification (pair channel),")
emit("S(e_a) at the anchored grid; e_a* = smallest e_a with S >= "
     "sigma_S (None if S(1) < sigma_S).")
emit("")
best = 0.0
best_corner = None
earned_corner = None
mono_ok = True
for xl in X_BIN:
    for gam in GAM_GRID:
        for kap in KAP_GRID:
            svals_ea = []
            for ea in (E_A_STATIC, 0.3, 1.0):
                svals_ea.append(S_stat(xl, x_amb_bin, kap, gam, ea,
                                       LAM[xl]))
            grid_ea = np.linspace(0.01, 1.0, 991)
            sg = np.array([S_stat(xl, x_amb_bin, kap, gam, e,
                                  LAM[xl]) for e in grid_ea])
            if np.any(np.diff(sg) < -1e-15):
                mono_ok = False
            s1v = float(sg[-1])
            if s1v > best:
                best = s1v
                best_corner = (xl, gam, kap)
            estars = {}
            for bar in (0.05, 0.02, 0.01):
                idx = np.where(sg >= bar)[0]
                estars[bar] = float(grid_ea[idx[0]]) if len(idx) \
                    else None
            if estars[0.02] is not None and (
                    earned_corner is None
                    or estars[0.02] < earned_corner[3]):
                earned_corner = (xl, gam, kap, estars[0.02])
            emit("  x=%.1f gam=%.3f kap=%.3f: S(0.086)=%.5f "
                 "S(0.3)=%.5f S(1)=%.5f | e_a*(5%%)=%s "
                 "e_a*(2%%)=%s e_a*(1%%)=%s"
                 % (xl, gam, kap, svals_ea[0], svals_ea[1],
                    svals_ea[2],
                    "%.3f" % estars[0.05] if estars[0.05] else "None",
                    "%.3f" % estars[0.02] if estars[0.02] else "None",
                    "%.3f" % estars[0.01] if estars[0.01] else "None"))
emit("")
emit("  convention family at the best corner (x=%.1f gam=%.3f "
     "kap=%.3f):" % best_corner)
for conv, cname in ((1/math.sqrt(2), "conv=1/sqrt2"),
                    (1.0, "conv=1 (quoted)"),
                    (math.sqrt(2), "conv=sqrt2")):
    sv = S_stat(best_corner[0], x_amb_bin, best_corner[2],
                best_corner[1], 1.0, LAM[best_corner[0]], conv)
    emit("    %-16s S(e_a=1) = %.5f" % (cname, sv))
emit("")
# background context row (perturbative all-channel sum, pair excluded)
def background(x_loc, x_amb, kap, gam, e_a, lam_max):
    om, Om = om_of(x_loc), om_of(x_amb)
    U2 = (e_a*lam_max)**2
    d2 = (kap/4)*(om/Om)
    P = [math.exp(-d2)*d2**m/math.factorial(m) for m in range(41)]
    tot = 0.0
    for n in range(6):
        pn = math.exp(-n*x_loc)*(1 - math.exp(-x_loc))
        for m in range(41):
            Dup = om*(1 - (kap/2)*(n + 1)) + m*Om
            if not (n == 1 and m == 0):
                tot += -pn*(n + 1)*U2*P[m]*Dup/(Dup**2 + gam**2)
            if n >= 1:
                Ddn = om*(1 - (kap/2)*n) - m*Om
                if not (n == 2 and m == 0):
                    tot += pn*n*U2*P[m]*Ddn/(Ddn**2 + gam**2)
    return abs(tot)/W_dict(x_loc, kap)
bg = background(1.0, x_amb_bin, 0.925, 0.015, 1.0, LAM[1.0])
emit("  background (all non-pair channels, x=1.0 gam=0.015 "
     "kap=0.925, e_a=1): |B| = %.5f" % bg)
emit("  (smooth in kappa; the kill statistic is the pair channel = "
     "the kappa-localized excess)")
emit("")

# ================= LETTER (mechanical) =================
g6 = True
if earned_corner is not None:
    xl, gam, kap, ea_star = earned_corner
    s_chk = S_stat(xl, x_amb_bin, kap, gam, ea_star, LAM[xl])
    g6 = abs(s_chk - 0.02) <= 1e-3
    letter = "N-ROW-EARNED"
elif best >= 0.002:
    g6 = (0.002 <= best <= 0.05) and mono_ok
    letter = "N-BELOW-FLOOR"
else:
    g6 = best < 0.002 and mono_ok
    letter = "N-RETIRED"
emit("G10N-6 positive-clause gate (%s): %s"
     % (letter, "PASS" if g6 else "FAIL"))
emit("")
emit("G10N-7 scope audit: dispersive-side only -- the mixing agent is")
emit("   the 10E l=2 vertex at SUB-SATURATION amplitude (the 10G")
emit("   strike killed only lambda >~ H saturation; ceilings used as")
emit("   UPPER bounds); the lending law appears NOWHERE; no sky fits")
emit("   (archived meters only); convention family spread disclosed")
emit("   above.")
emit("")
emit("=" * 68)
emit("VERDICT LETTER: %s" % letter)
if letter == "N-ROW-EARNED":
    xl, gam, kap, ea_star = earned_corner
    emit("   The kill condition is DERIVED and a DR4-grade window")
    emit("   exists inside the permitted amplitude band: at the")
    emit("   anchored corner x=%.1f, gamma=%.3f H, kappa=%.3f the"
         % (xl, gam, kap))
    emit("   crossing signature reaches sigma_S = 0.02 at e_a* = "
         "%.3f <= 1." % ea_star)
    emit("   THE REGISTERED FORM (books commit -> PREDICTIONS P10):")
    emit("   IF the deep-limit kappa sits in the lock band [0.888,")
    emit("   1.0] AND the ambient l=2 fluctuation amplitude e_a >=")
    emit("   e_a* AND the dispersive reading holds, THEN wide-binary")
    emit("   pairs in the transition window x in [0.5, 1.2] carry an")
    emit("   ambient-quadrupole-correlated boost modification of")
    emit("   size S (table above, up to ~%.3f at e_a = 1), absent"
         % best)
    emit("   for kappa away from the crossing.  KILL: a DR4-era")
    emit("   transition-window split at sigma_S <= 0.02 showing no")
    emit("   excess kills the (kappa near 1) x (e_a >= e_a*) cell;")
    emit("   a detected x-localized, tide-correlated excess of the")
    emit("   predicted shape is the crossing's fingerprint.")
    emit("   TWO DR4-EXECUTABLE STATISTICS: (i) the x-window boost")
    emit("   residual split by ambient-quadrupole proxy (R_gal /")
    emit("   Z-height bins); (ii) the eccentricity-resolved P7")
    emit("   channel restricted to the window (distinct prediction:")
    emit("   the excess rides e_a, not e).")
elif letter == "N-BELOW-FLOOR":
    emit("   The kill condition is DERIVED (the prerequisite is")
    emit("   discharged) but the signature tops out at S = %.4f <"
         % best)
    emit("   0.02 at every anchored corner even at e_a = 1: NO")
    emit("   PREDICTIONS row (below the pre-set DR4 floor).  The")
    emit("   surface above is the standing design spec; revival bar")
    emit("   sigma_S^rev = %.4f (the precision at which a" % best)
    emit("   transition-window split statistic would begin to read")
    emit("   the crossing at the most favorable anchored corner")
    emit("   x=%.1f, gamma=%.3f, kappa=%.3f)." % best_corner)
else:
    emit("   The crossing response is below every conceivable floor")
    emit("   (best S(e_a=1) = %.5f < 0.002): the R30 reading" % best)
    emit("   RETIRES; no PREDICTIONS row.")
emit("")
emit("STRUCTURAL FINDINGS (letter-independent, gated above):")
emit("   (1) the crossing meter is BINARY-NATIVE: deep galaxies die")
emit("       by occupation contrast (x%.0f), transition galaxies by"
     % (1/ratio_deep))
emit("       comb smearing (gamma >= Om); the window function is")
emit("       (p1-p2)(x) = e^-x(1-e^-x)^2, max at x = ln 3 -- the")
emit("       wide-binary regime.")
emit("   (2) the response is WIDTH-SATURATED: at the lock kappa the")
emit("       pair gap D1 = om(1-kappa) sits at or below the R17")
emit("       width band, so the kappa -> 1 divergence is cut at")
emit("       ~om/gamma; distinguishing kappa = 0.925 from 1.0 by")
emit("       response size alone is width-limited (the peak sits at")
emit("       1 - kappa* = gamma/om).")
emit("pre-signed credence: ALL cells HOLD mech-conditional 8;")
emit("anomaly-real 53 untouched (no sky fits).")
emit("gates: G1 %s | G2 %s | G3 %s | G4 %s | G5 %s | G6 %s | G7 "
     "stated" % tuple("PASS" if g else "FAIL"
                      for g in (g1, g2, g3 and gt_disq, g4, g5, g6)))

with open(OUT, "w") as f:
    f.write("\n".join(L) + "\n")
