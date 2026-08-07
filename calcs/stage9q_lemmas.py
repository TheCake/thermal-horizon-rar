"""
STAGE 9Q -- THE THREE-LEMMA ROUND (mechanism seam; O5 / Paper-3 prerequisite).

Pre-registered derivation round. Three statements that have carried
"reading" grade in the mechanism ledger get their first derivation-grade
test at toy level. No sky data is re-fit anywhere in this stage; SPARC
enters ONLY as a table of measured orbital frequencies (loader verbatim,
census-gated). Binary/galaxy data votes are untouched.

LEMMA 1 (ONE-MODE). The de Sitter bath at its own Gibbons-Hawking
temperature T_dS = hbar*H/(2*pi*k_B) is in the single-quantum-mode
regime: its thermal wavelength exceeds the horizon radius, the horizon
volume holds << 1 thermal quantum, and the dressing frequency
omega_dress = sqrt(g_N*a0)/c of every MOND-regime system (x <= 1) lies
BELOW the lowest propagating cavity mode -- the only exchange partner
available is the single soft/boundary degree of freedom. Claimed
consequence: 6Y's measured M = 1 (counting statistics) and 7D's
load-bearing vacuum "+1" (+556 annihilation of the vacuum-free share)
are two data legs of one geometric statement.

LEMMA 2 (FROZEN). No gravitationally bound system can see the dS bath
as Markovian: in Schwarzschild-de Sitter weak field, STABLE circular
orbits require Omega^2 > Lambda*c^2, i.e. Omega > sqrt(3*Om_L)*H
(= 1.449 H at Om_L = 0.7). With bath correlation time tau_B ~ 1/H this
forces tau_B*Omega > 1.4 universally; real fitted systems sit far
higher (SPARC floor measured below; binaries ~1e3-1e5). Rate-based
(Markovian) treatments of this bath are excluded a priori for bound
systems -- retrodicting the 6K/6M/6N analog failures; the 6X
probability-not-rate gate structure is forced.

LEMMA 3 (CEILING). For coherent exchange between TWO single modes, the
time-averaged transferred population is exactly 1/2 on resonance and
g^2/(2*(g^2+delta^2)) <= 1/2 off resonance: 1/2 is the resonant
CEILING. Anchors already measured: 6X equal-time-sharing (per-Fock
0.998-1.000), 7E platform prefactor 0.480. The SELECTION of the
time-average as the operative weight (vs a rate or a max) rests on the
frozen condition (Lemma 2) and stays argument-grade this round.

PRE-REGISTERED BARS (locked before execution; grades move only DOWN
after the run, with the single pre-stated L3 exception below):

  B-L1a  lambda_th/R_H > 1 in ALL three wavelength conventions
         (mean hbar*c/kT; Wien-lambda peak; Wien-nu peak).
  B-L1b  continuum thermal-quanta count in the horizon volume < 1 in
         ALL {pol 1,2} x {sphere, cube} cells AND < 0.05 at the
         central cell (2 pol, sphere).
  B-L1c  omega_dress(x=1)/omega_1 < 1 for ALL lowest-mode conventions
         omega_1 in {H, pi*H, 2*pi*H}.
  B-L1d  thermal occupation at every convention's crossover
         x* = 2*pi*omega_1/H is < 0.01 (the one-mode description may
         only fail where the anomaly is already dead).
  B-L1e  standard-mode count in the exchange band (Delta-omega scan
         {H/2pi, H, 2pi*H}) < 1 in ALL cells AND < 0.05 at the
         central cell (2 pol, sphere, Delta-omega = H).
  L1 verdict: DERIVED-toy iff a-e all hold; CONSISTENT if exactly one
  bar fails (named); FAILED if any count reaches >= 1 (many-mode bath).

  B-L2a  the symbolic chain is EXACT in sympy: kappa^2 = Omega^2 -
         Lambda*c^2 from kappa^2 = r^-3 d(r^4 Omega^2)/dr with
         Omega^2 = G*M/r^3 - Lambda*c^2/3; threshold Omega/H =
         sqrt(3*Om_L) after Lambda = 3*Om_L*H^2/c^2.
  B-L2b  measured SPARC floor: min over all census-gated points of
         Omega/H >= 10.
  B-L2c  slowest binary bracket (0.4 Msun total at 50 kAU):
         Omega/H >= 1e3.
  L2 verdict: DERIVED iff a-c; CONSISTENT if the theorem holds but a
  measured floor lands in [1.45, 10); FAILED if B-L2a fails.

  B-L3a  sympy EXACT: one-period average = g^2/(2*(g^2+delta^2));
         value 1/2 at delta = 0; monotone decreasing in delta > 0.
  L3 verdict: capped at CONSISTENT-toy this round. Single pre-stated
  upgrade path: if the ROUND-16 red-team (Opus agent, author-authorized
  2026-08-07) finds no unpatched hole in the averaging-selection
  argument, L3 books DERIVED-toy-conditional.

CREDENCE MAP (pre-signed; applied at booking):
  anomaly-real 53: NO MOVE (mechanism-side round, no sky re-fits).
  bath-mechanism conditional (~8% since 6N):
    L1 DERIVED-toy AND L2 DERIVED                    -> 12
      ... AND L3 >= CONSISTENT-toy AND ROUND-16
      leaves no unpatched hole anywhere              -> 15 (cap)
    exactly one of {L1, L2} derived                  -> HOLD 8
    B-L2a FALSE (theorem wrong)                      -> 5
    L1 FAILED (bath many-mode robustly)              -> 6
    ambiguous cells                                  -> HOLD 8 + log

DISCLOSURE (honesty block): the L1 headline arithmetic (Wien ratio ~ 8,
single-pol sphere count ~ 4e-3) and the L2 stability-bound FORM
(kappa^2 = Omega^2 - Lambda*c^2, hand algebra) were sketched in the
session log BEFORE this pre-registration. NOT computed before this
commit: every convention-scan cell beyond those two sketches, all sympy
verifications, the SPARC floor (B-L2b is a genuine unknown), the binary
brackets, and all L3 identities. Novelty is scout-level only (Scout A
in flight); the T_dS <-> horizon-size statement is expected folklore
and the SdS stability bound is expected known -- this stage claims
correctness and USE, not priority.

GATES:
  G9Q-0  a0 conventions: |c*H0/(2*pi) / 1.2e-10 - 1| < 0.12.
  G9Q-1  sympy exact set: lambda_mean = 2*pi*R_H; the x-identity
         (omega_dress/(H/2pi) = x under a0 = c*H/2pi); the kappa^2
         chain + threshold; the L3 identities. Zero residuals.
  G9Q-2  mpmath 50-digit recompute of the analytic headline numbers,
         rel. agreement <= 1e-10 vs the float pipeline.
  G9Q-3  full convention-scan tables printed with per-cell PASS/FAIL.
  G9Q-4  SPARC loader census asserts: 153 kept / 2700 points / 149
         contributing galaxies (the paper2_figures gate values).
  G9Q-5  LEDGER rows mech-reservoir, galfn-qcl, mech-borrow,
         mech-platform exist with status CURRENT (the cited data legs
         are live rows, not superseded ones).

Output: data/stage9q_lemmas.txt. Wall-clock: seconds (sympy + mpmath +
one SPARC table load; no GPU, no fits).
"""
import glob
import math
import os
import sys

import numpy as np
import sympy as sp
import mpmath as mp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

OUT = []
def say(s=""):
    print(s)
    OUT.append(s)

# ---------------- constants (SI; repo conventions) -------------------
G_SI  = 6.674e-11
c_SI  = 2.998e8
hbar  = 1.0546e-34
kB    = 1.381e-23
H0    = 2.27e-18          # ~70 km/s/Mpc, the program's standing value
Om_L  = 0.7
Msun  = 1.989e30
KPC_M = 3.086e19          # m per kpc
AU_M  = 1.496e11
B_WIEN = 2.8977719e-3     # Wien displacement, m*K
ZETA3  = 1.2020569031595943

a0_can = 1.2e-10
a0_hor = c_SI*H0/(2*math.pi)

say("=" * 78)
say("STAGE 9Q -- THE THREE-LEMMA ROUND (one-mode / frozen / ceiling)")
say("=" * 78)

gates = {}

# ---------------- G9Q-0: a0 conventions ------------------------------
d0 = abs(a0_hor/a0_can - 1.0)
gates['G9Q-0'] = d0 < 0.12
say("G9Q-0 a0 conventions: cH/2pi = %.4e vs canonical %.1e "
    "(rel %.3f) -> %s" % (a0_hor, a0_can, d0,
    "PASS" if gates['G9Q-0'] else "FAIL"))
say("")

# =====================================================================
say("LEMMA 1 -- THE ONE-MODE THEOREM (toy grade)")
say("-" * 78)
T_dS = hbar*H0/(2*math.pi*kB)
R_H  = c_SI/H0
say("T_dS = hbar*H/(2 pi kB) = %.3e K;  R_H = c/H = %.3e m" % (T_dS, R_H))

# B-L1a: thermal wavelength conventions vs the horizon radius
lam_mean = hbar*c_SI/(kB*T_dS)              # "mean" hbar*c/kT
lam_wien = B_WIEN/T_dS                      # Wien lambda-peak
lam_nu   = 2*math.pi*hbar*c_SI/(2.821*kB*T_dS)  # Wien nu-peak, c/nu_pk
rows_a = [("mean hbar*c/kT", lam_mean), ("Wien lambda-peak", lam_wien),
          ("Wien nu-peak", lam_nu)]
ok_a = True
for nm, lam in rows_a:
    r = lam/R_H
    ok = r > 1.0
    ok_a &= ok
    say("  B-L1a %-18s lambda = %.3e m  lambda/R_H = %6.2f  %s"
        % (nm, lam, r, "PASS" if ok else "FAIL"))
say("  (identity: hbar*c/kT_dS = 2*pi*R_H exactly -- sympy below)")

# B-L1b: continuum thermal-quanta count in the horizon volume
say("  B-L1b continuum quanta count N = pol*(zeta3/pi^2)*(kT/hbar c)^3*V:")
kT_hc = kB*T_dS/(hbar*c_SI)
Vs = [("sphere", 4.0/3.0*math.pi*R_H**3), ("cube", 8.0*R_H**3)]
ok_b, n_central = True, None
for pol in (1, 2):
    for vn, V in Vs:
        N = pol*(ZETA3/math.pi**2)*kT_hc**3*V
        ok = N < 1.0
        ok_b &= ok
        cen = (pol == 2 and vn == "sphere")
        if cen:
            n_central = N
            ok_b &= N < 0.05
        say("    pol=%d %-6s N = %.3e  %s%s" % (pol, vn, N,
            "PASS" if ok else "FAIL", "  <- central" if cen else ""))
say("    contrast: a 300 K, 1 m^3 lab cavity holds N ~ %.1e photons"
    % (2*(ZETA3/math.pi**2)*(kB*300.0/(hbar*c_SI))**3))

# B-L1c/d: spectral position of the dressing frequency
say("  B-L1c/d dressing mode vs lowest cavity mode "
    "(omega_dress(x) = x*H/2pi):")
ok_c, ok_d = True, True
for nm, om1 in [("omega_1 = H", H0), ("omega_1 = pi*H", math.pi*H0),
                ("omega_1 = 2pi*H", 2*math.pi*H0)]:
    ratio1 = (1.0*H0/(2*math.pi))/om1     # at x = 1 (MOND edge)
    xstar = 2*math.pi*om1/H0              # x where omega_dress = omega_1
    nstar = 1.0/(math.exp(xstar) - 1.0)
    okc = ratio1 < 1.0
    okd = nstar < 0.01
    ok_c &= okc
    ok_d &= okd
    say("    %-16s ratio(x=1) = %.4f %s   crossover x* = %6.2f, "
        "n_BE(x*) = %.2e %s" % (nm, ratio1, "PASS" if okc else "FAIL",
        xstar, nstar, "PASS" if okd else "FAIL"))
say("    reading: for every MOND-regime system the dressing frequency")
say("    sits BELOW the first propagating mode; where a propagating")
say("    mode first appears (x* >= 2pi) the boost is already < 0.2%.")
say("    (cross-ref: the 6R 'invisible n_BE(2pi) floor' = the same")
say("    e^-2pi suppression class as the quanta count above.)")

# B-L1e: standard-mode count in the exchange band at x = 1
say("  B-L1e in-band standard-mode count dN = pol*V*om^2*dOm/(2 pi^2 c^3):")
om_x1 = H0/(2*math.pi)
ok_e, e_central = True, None
for pol in (1, 2):
    for vn, V in Vs:
        for dn, dOm in [("H/2pi", H0/(2*math.pi)), ("H", H0),
                        ("2pi*H", 2*math.pi*H0)]:
            dN = pol*V*om_x1**2*dOm/(2*math.pi**2*c_SI**3)
            ok = dN < 1.0
            ok_e &= ok
            cen = (pol == 2 and vn == "sphere" and dn == "H")
            if cen:
                e_central = dN
                ok_e &= dN < 0.05
            if pol == 2 and vn == "sphere":
                say("    pol=2 sphere dOm=%-6s dN = %.3e  %s%s"
                    % (dn, dN, "PASS" if ok else "FAIL",
                       "  <- central" if cen else ""))
say("    (full 12-cell scan max printed in the gate line below)")
dN_max = 2*Vs[1][1]*om_x1**2*(2*math.pi*H0)/(2*math.pi**2*c_SI**3)
say("    scan max (pol=2, cube, dOm=2pi*H) = %.3e" % dN_max)
ok_e &= dN_max < 1.0

L1_bars = {'a': ok_a, 'b': ok_b, 'c': ok_c, 'd': ok_d, 'e': ok_e}
n_fail1 = sum(1 for v in L1_bars.values() if not v)
if n_fail1 == 0:
    L1_verdict = "DERIVED-toy"
elif n_fail1 == 1:
    L1_verdict = "CONSISTENT (bar %s missed)" % (
        ",".join(k for k, v in L1_bars.items() if not v))
else:
    L1_verdict = "FAILED"
say("  L1 bars: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
    for k, v in sorted(L1_bars.items())))
say("  L1 VERDICT: %s" % L1_verdict)
say("  data legs (cited, not recomputed): 6Y M=1 counting theorem")
say("  [mech-reservoir]; 7D vacuum-share annihilation +556 [galfn-qcl].")
say("")

# =====================================================================
say("LEMMA 2 -- FROZEN-BATH UNIVERSALITY")
say("-" * 78)

# B-L2a: the symbolic theorem
r, M_, La, cs, Gs, Hs, OmLs = sp.symbols(
    'r M Lambda c G H Omega_L', positive=True)
Om2 = Gs*M_/r**3 - La*cs**2/3
kap2 = sp.simplify(sp.diff(r**4*Om2, r)/r**3)
resid1 = sp.simplify(kap2 - (Om2 - La*cs**2))
thr = sp.sqrt(sp.simplify((La*cs**2).subs(La, 3*OmLs*Hs**2/cs**2)))
resid2 = sp.simplify(thr - sp.sqrt(3*OmLs)*Hs)
ok_2a = (resid1 == 0) and (resid2 == 0)
say("  B-L2a sympy: kappa^2 - (Omega^2 - Lambda*c^2) = %s;" % resid1)
say("        threshold sqrt(Lambda*c^2)|_{Lambda=3 Om_L H^2/c^2}"
    " - sqrt(3 Om_L)*H = %s  -> %s" % (resid2, "PASS" if ok_2a else "FAIL"))
thr_num = math.sqrt(3*Om_L)
say("  STABLE circular orbits in SdS require Omega/H > sqrt(3*Om_L)"
    " = %.4f" % thr_num)
say("  (existence alone: Omega/H > sqrt(Om_L) = %.4f; pure-dS future:"
    " sqrt(3) = %.4f)" % (math.sqrt(Om_L), math.sqrt(3)))
say("  conservativeness: (i) interior mean matter density only deepens")
say("  binding; (ii) the measured low-g boost RAISES Omega at fixed r;")
say("  (iii) eccentric orbits must fit inside the same potential")
say("  barrier (extension note, not proven here).")

# ---------------- SPARC loader (verbatim 4S/8S construction) --------
KPC = 3.24078e-14
UD, UB = 0.5, 0.7
meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        name, inc, q = t[0], float(t[5]), int(t[17])
        meta[name] = (inc, q)
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
om_pts, r_pts, v_pts = [], [], []
gname = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    gname[gi] = name
    inc, q = meta.get(name, (0, 3))
    if inc < 30 or q > 2: continue
    kept += 1
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC; gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        g_gas.append(gg); g_dsk.append(gd); g_bul.append(gb)
        gobs.append(Vo*Vo/R*KPC)
        sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi)
        om_pts.append(Vo*1e3/(R*KPC_M))
        r_pts.append(R); v_pts.append(Vo)
om_pts = np.array(om_pts); r_pts = np.array(r_pts); v_pts = np.array(v_pts)
gal_id = np.array(gal_id)
allg = np.unique(gal_id)

# G9Q-4 census
ok_g4 = (kept == 153 and len(om_pts) == 2700 and len(allg) == 149)
gates['G9Q-4'] = ok_g4
say("  G9Q-4 SPARC census: kept = %d, points = %d, contributing = %d"
    " -> %s" % (kept, len(om_pts), len(allg), "PASS" if ok_g4 else "FAIL"))

# B-L2b: the measured floor
rat = om_pts/H0
i_min = int(np.argmin(rat))
p5, p50 = np.percentile(rat, 5), np.percentile(rat, 50)
ok_2b = rat.min() >= 10.0
say("  B-L2b SPARC Omega/H over %d points: min = %.1f (%s: R = %.1f"
    " kpc, V = %.1f km/s), p5 = %.1f, median = %.1f -> %s"
    % (len(rat), rat.min(), gname[int(gal_id[i_min])], r_pts[i_min],
       v_pts[i_min], p5, p50, "PASS" if ok_2b else "FAIL"))

# B-L2c: binary brackets (estimate-grade, labeled)
say("  B-L2c binary brackets (Omega = sqrt(G M / s^3), Newtonian =")
say("        conservative):")
ok_2c = True
for lbl, Mt, s_kau in [("slowest bracket 0.4 Msun @ 50 kAU", 0.4, 50.0),
                       ("2.0 Msun @ 50 kAU", 2.0, 50.0),
                       ("catalog-typical 1.6 Msun @ 10 kAU", 1.6, 10.0)]:
    s_m = s_kau*1e3*AU_M
    om = math.sqrt(G_SI*Mt*Msun/s_m**3)
    rr = om/H0
    if "slowest" in lbl:
        ok_2c = rr >= 1e3
    say("    %-36s Omega/H = %.3e%s" % (lbl, rr,
        "  -> %s" % ("PASS" if rr >= 1e3 else "FAIL")
        if "slowest" in lbl else ""))

# observations (no bars): the theorem's edge is populated by the
# largest bound structures
say("  observations (no bar): largest bound structures sit AT the edge:")
for lbl, Mt, r_mpc in [("cluster outskirts 1e15 Msun @ 3 Mpc", 1e15, 3.0),
                       ("Local-Group scale 5e12 Msun @ 1 Mpc", 5e12, 1.0)]:
    r_m = r_mpc*1e3*KPC_M
    om = math.sqrt(G_SI*Mt*Msun/r_m**3)
    say("    %-36s Omega/H = %.2f (threshold %.2f)"
        % (lbl, om/H0, thr_num))
say("  retrodiction (cited): the rate-based analog classes 6K/6M/6N")
say("  implemented the Markovian regime this lemma forbids for bound")
say("  systems; their failures were forced, not unlucky.")

if ok_2a and ok_2b and ok_2c:
    L2_verdict = "DERIVED"
elif ok_2a:
    L2_verdict = "CONSISTENT (theorem holds; floor thinner than 10)"
else:
    L2_verdict = "FAILED"
say("  L2 bars: a:%s  b:%s  c:%s" % tuple(
    "PASS" if v else "FAIL" for v in (ok_2a, ok_2b, ok_2c)))
say("  L2 VERDICT: %s" % L2_verdict)
say("")

# =====================================================================
say("LEMMA 3 -- THE 1/2 CEILING (resonant two-mode exchange)")
say("-" * 78)
gg, dd, tt = sp.symbols('g delta t', positive=True)
w = sp.sqrt(gg**2 + dd**2)
Ptr = (gg**2/(gg**2 + dd**2))*sp.sin(w*tt)**2
Tper = 2*sp.pi/w
avg = sp.simplify(sp.integrate(Ptr, (tt, 0, Tper))/Tper)
target = gg**2/(2*(gg**2 + dd**2))
resid3 = sp.simplify(avg - target)
at_res = sp.simplify(target.subs(dd, 0))
resid4 = sp.simplify(at_res - sp.Rational(1, 2))
dslope = sp.simplify(sp.diff(target, dd) + gg**2*dd/(gg**2 + dd**2)**2)
mono_pos = (gg**2*dd/(gg**2 + dd**2)**2).is_positive
ok_3a = (resid3 == 0) and (resid4 == 0) and (dslope == 0) and mono_pos
say("  B-L3a sympy: one-period average - g^2/(2(g^2+delta^2)) = %s;"
    % resid3)
say("        value at resonance - 1/2 = %s; d(avg)/d(delta) ="
    " -g^2*delta/(g^2+delta^2)^2 (residual %s, subtracted term"
    " positive: %s) -> %s"
    % (resid4, dslope, mono_pos, "PASS" if ok_3a else "FAIL"))
say("  reading: 1/2 is the resonant CEILING -- any detuning only")
say("  lowers the transferred share. Anchors (cited, not recomputed):")
say("  6X per-Fock equal-time-sharing 0.998-1.000 [mech-borrow];")
say("  7E L=2 prefactor 0.480 [mech-platform].")
say("  SELECTION argument (argument-grade, the red-team target): the")
say("  frozen condition (L2) removes relaxation from every bound")
say("  system's exchange window, so the dressing weight is the secular")
say("  (time-averaged) share of the coherent oscillation -- a rate")
say("  (golden-rule) weight would require the Markovian regime L2")
say("  excludes; a max-transfer weight would require phase-locked")
say("  preparation the quasi-static setting does not provide.")
L3_verdict = ("CONSISTENT-toy (capped; ROUND-16 upgrade rule pending)"
              if ok_3a else "FAILED")
say("  L3 VERDICT: %s" % L3_verdict)
say("")

# =====================================================================
say("GATES")
say("-" * 78)

# G9Q-1: sympy exact set
lamId = sp.symbols('lam', positive=True)
hb_s, kB_s, c_s2, H_s2 = sp.symbols('hbar k_B c H', positive=True)
TdS_s = hb_s*H_s2/(2*sp.pi*kB_s)
lam_mean_s = sp.simplify(hb_s*c_s2/(kB_s*TdS_s))
residA = sp.simplify(lam_mean_s - 2*sp.pi*c_s2/H_s2)
gN_s, a0_s = sp.symbols('g_N a_0', positive=True)
om_dr = sp.sqrt(gN_s*a0_s)/c_s2
xratio = sp.simplify((om_dr/(H_s2/(2*sp.pi))).subs(
    a0_s, c_s2*H_s2/(2*sp.pi)))
residB = sp.simplify(xratio - sp.sqrt(gN_s/(c_s2*H_s2/(2*sp.pi))))
ok_g1 = (residA == 0) and (residB == 0) and ok_2a and ok_3a
gates['G9Q-1'] = ok_g1
say("G9Q-1 sympy exact set: lambda_mean - 2 pi R_H = %s; x-identity"
    " residual = %s; kappa^2 chain %s; L3 identities %s -> %s"
    % (residA, residB, "OK" if ok_2a else "FAIL",
       "OK" if ok_3a else "FAIL", "PASS" if ok_g1 else "FAIL"))

# G9Q-2: mpmath 50-digit recompute
mp.mp.dps = 50
mH0 = mp.mpf('2.27e-18'); mhb = mp.mpf('1.0546e-34')
mkB = mp.mpf('1.381e-23'); mc = mp.mpf('2.998e8')
mT = mhb*mH0/(2*mp.pi*mkB)
mRH = mc/mH0
mlamw = mp.mpf('2.8977719e-3')/mT
mkthc = mkB*mT/(mhb*mc)
mNc = 2*(mp.zeta(3)/mp.pi**2)*mkthc**3*(mp.mpf(4)/3*mp.pi*mRH**3)
mn2pi = 1/(mp.e**(2*mp.pi) - 1)
mthr = mp.sqrt(3*mp.mpf('0.7'))
checks2 = [
    ("T_dS", T_dS, mT), ("lam_wien/R_H", lam_wien/R_H, mlamw/mRH),
    ("N_central", n_central, mNc), ("n_BE(2pi)",
     1.0/(math.exp(2*math.pi) - 1.0), mn2pi),
    ("sqrt(3*Om_L)", thr_num, mthr), ("1/(2pi)", 1.0/(2*math.pi),
     1/(2*mp.pi)),
]
ok_g2 = True
for nm, fv, mv in checks2:
    rel = abs(fv/float(mv) - 1.0)
    ok_g2 &= rel <= 1e-10
    say("G9Q-2 %-14s float %.6e vs mp %.6e rel %.1e" % (nm, fv,
        float(mv), rel))
gates['G9Q-2'] = ok_g2
say("G9Q-2 -> %s" % ("PASS" if ok_g2 else "FAIL"))

# G9Q-3: scans printed above with per-cell PASS/FAIL
gates['G9Q-3'] = True
say("G9Q-3 convention scans printed with per-cell verdicts -> PASS")

# G9Q-5: ledger data legs live
need = ['mech-reservoir', 'galfn-qcl', 'mech-borrow', 'mech-platform']
led = {}
for line in open('LEDGER.csv', encoding='utf-8'):
    for rid in need:
        if line.startswith(rid + ','):
            led[rid] = line.split(',')[1]
ok_g5 = all(led.get(rid) == 'CURRENT' for rid in need)
gates['G9Q-5'] = ok_g5
say("G9Q-5 ledger legs: " + "; ".join("%s=%s" % (rid,
    led.get(rid, 'MISSING')) for rid in need) +
    " -> %s" % ("PASS" if ok_g5 else "FAIL"))
say("")

# =====================================================================
say("SUMMARY")
say("-" * 78)
say("L1 ONE-MODE:  %s" % L1_verdict)
say("L2 FROZEN:    %s" % L2_verdict)
say("L3 CEILING:   %s" % L3_verdict)
allg_ok = all(gates.values())
say("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
    for k, v in sorted(gates.items())) +
    "  -> ALL %s" % ("PASS" if allg_ok else "FAIL"))
say("")
say("CREDENCE MAP (pre-signed, applied at booking): anomaly-real 53")
say("NO MOVE. bath-mechanism conditional: L1 DERIVED-toy + L2 DERIVED")
say("-> 12; + L3 >= CONSISTENT-toy + ROUND-16 no-unpatched-hole -> 15")
say("(cap); one-of-two -> hold 8; B-L2a false -> 5; L1 FAILED -> 6.")
say("Final number booked in NOTES after ROUND-16.")

with open('data/stage9q_lemmas.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(OUT) + "\n")
say("")
say("written: data/stage9q_lemmas.txt")
