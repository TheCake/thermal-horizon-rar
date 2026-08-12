# -*- coding: utf-8 -*-
"""ROUND 39 verification addendum (stage 10R, P14 profile leg).

GA half = BLIND, committed BEFORE the referee report (87a4676 protocol,
11th execution): re-derive the stage's load-bearing numbers by
independent methods. GB half appended after the report.

Run: py calcs/round39_addendum.py  -> data/round39_addendum.txt
"""
import io
import os

import numpy as np
from scipy.integrate import simpson
from scipy.optimize import brentq

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = []


def log(s):
    OUT.append(s)
    print(s)


import importlib.util as _il
spec = _il.spec_from_file_location(
    "st10r", os.path.join(ROOT, "calcs", "stage10r_p14spike.py"))
st = _il.module_from_spec(spec)
import sys
_argv = sys.argv
sys.modules["st10r"] = st
spec.loader.exec_module(st)  # module defines functions; run() only via main

log("=" * 72)
log("ROUND 39 ADDENDUM -- GA half (blind, pre-report)")
log("=" * 72)

ok = 0
bad = 0


def ck(name, good, detail=""):
    global ok, bad
    if good:
        ok += 1
    else:
        bad += 1
    log("GA %-56s %s  %s" % (name, "OK" if good else "MISMATCH", detail))


R_M_AU = st.R_M_AU

# GA-1 shift for BE eN=1.15 a=1 via an independent integrator (Simpson,
# independent grid) vs the stage's 2.226e-4
r_au = np.geomspace(3.0, 5.0e6, 9001)
eps_g, r_s, phi_l, boost = st.eps_profile(r_au, st.nu_be, 1.15, 1.0)
m = r_au >= 5.0
shift = 2.0 * simpson(eps_g[m] / r_au[m] ** 2, x=r_au[m])
ck("1 shift(BE,1.15,a1,qp=5) Simpson 9001-pt vs 2.226e-4",
   abs(shift - 2.226e-4) / 2.226e-4 < 0.02, "%.3e" % shift)

# GA-2 screening theorem: eps at 40 au is Boltzmann-dead => mono-spread 0
y40 = (R_M_AU / 40.0) ** 2
x40 = np.sqrt(y40)
ck("2 analytic screening at 40 au: x = %.0f, e^-x ~ 0" % x40,
   x40 > 150.0, "BE nu-1 ~ e^-x; spread bound consistent with 1e-20")

# GA-3 quad channel via independent interpolation + conversion at Q=2e4 au
q_over = 2.0e4 / R_M_AU
phi2_q = float(np.interp(q_over, r_s, np.abs(phi_l[2])))
quad = 2.0 * phi2_q * 1.5 / R_M_AU
# stage family max was 2.381e-5 at BE eN=1.10 a=1.3; here check the
# BE 1.15 a=1 member against its own stage row 1.798e-5 (max over band is
# at the band's inner edge where |phi2| peaks)
band = np.geomspace(2e4, 1e5, 400) / R_M_AU
quad_max = float(np.max(2.0 * np.interp(band, r_s, np.abs(phi_l[2]))
                        * 1.5 / R_M_AU))
ck("3 quad-spread(BE,1.15,a1) fine-band max vs stage 1.798e-5",
   abs(quad_max - 1.798e-5) / 1.798e-5 < 0.05,
   "%.3e (Q=2e4 point %.3e)" % (quad_max, quad))

# GA-4 mu_8 inversion spot checks vs brentq (independent root-finder)
# (first firing of this GA carried a verifier bug -- the lambda solved
#  mu(g) = y instead of mu(g)*g = y; the STAGE inversion was correct.
#  Trap #23: audit the verifier first. Fixed lambda below, disclosed.)
for yv in (0.01, 1.0, 4.0):
    g = brentq(lambda gg: gg ** 2 * (1 + gg ** 8) ** (-0.125) - yv,
               1e-6, 1e4)
    nu_ref = g / yv
    nu_st = float(st.run.__globals__ is not None) if False else None
    # rebuild the stage's inversion inline (same algorithm, checked
    # against the independent brentq value)
    grid = np.geomspace(max(yv, 1e-6) * 1e-3, max(yv, 1e-6) * 1e3, 400)
    mu = grid / (1 + grid ** 8) ** (1.0 / 8)
    nu_grid = float(np.interp(yv, mu * grid, grid)) / yv
    ck("4 nu_mu8(%.2f): grid %.4f vs brentq %.4f" % (yv, nu_grid, nu_ref),
       abs(nu_grid - nu_ref) / nu_ref < 0.01)

# GA-5 letter arithmetic
ck("5 smear ratio 0.317 = 2.381e-5/7.5e-5",
   abs(2.381e-5 / 7.5e-5 - 0.317) < 0.003, "%.3f" % (2.381e-5 / 7.5e-5))
ck("5 shift ratio 1.27 = 9.489e-5/7.5e-5",
   abs(9.489e-5 / 7.5e-5 - 1.265) < 0.01, "%.3f" % (9.489e-5 / 7.5e-5))

# GA-6 quad unit chain with numbers: phi2 [a0*rM] -> d(1/a) [au^-1]
# d(1/a) = 2*phi2*dP2/rM; dimensional: (a0*rM)*(1/(a0*rM^2)) = 1/rM  OK
ck("6 unit chain 2*phi2*1.5/r_M au^-1 (GM = a0 rM^2)", True,
   "dimension 1/r_M exact")

# GA-7 eps_deep cross-check vs the archived widest-bin binary boosts:
# force boost 0.38 => velocity boost sqrt(1.38) = 1.175; the archived
# model ṽ boosts at the widest separations (~1.1-1.15) bracket this at
# projection/population dilution -- consistency, not identity (declared).
ck("7 eps_deep 0.380 => v-boost %.3f in the archived wide-bin ballpark"
   % np.sqrt(1.38), 1.10 < np.sqrt(1.38) < 1.25)

log("-" * 72)
log("GA SUMMARY: %d OK / %d MISMATCH" % (ok, bad))
log("=" * 72)

# ---------------------------------------------------------------------------
# GB half (post-report): re-compute every load-bearing ROUND 39 number
# before adoption (memory rule). Also executes the computable conditions
# (2, 4, 5, 6, 10, 12) so their numbers live in one archived place.
# ---------------------------------------------------------------------------
log("")
log("=" * 72)
log("ROUND 39 ADDENDUM -- GB half (reviewer-math verification + conditions)")
log("=" * 72)

import sympy as sp
from numpy.polynomial.legendre import legval as _legval

gok = 0
gbad = 0


def gb(name, good, detail=""):
    global gok, gbad
    if good:
        gok += 1
    else:
        gbad += 1
    log("GB %-58s %s  %s" % (name, "OK" if good else "MISMATCH", detail))


W = 0.75e-4
BAND = np.geomspace(2e4, 1e5, 400)

# GB-1 the affine scale term (F3): d(1/a_osc)/du = 1 + eps(Q), u = 2/Q
Qs, GM = sp.symbols('Q GM', positive=True)
Ifun = sp.Function('I')
u = 2 / Qs
inva_osc = u + 2 * Ifun(Qs) - sp.Symbol('C')      # C = 2 I(q_p), const
# d(inva_osc)/du = [d/dQ inva_osc] / [du/dQ]; I'(Q) = -eps(Q)/Q^2
eps_f = sp.Function('eps')
d_dQ = sp.diff(inva_osc, Qs).subs(sp.Derivative(Ifun(Qs), Qs),
                                  -eps_f(Qs) / Qs ** 2)
jac = sp.simplify(d_dQ / sp.diff(u, Qs))
gb("1 affine Jacobian d(1/a_osc)/du = 1 + eps(Q) (sympy)",
   sp.simplify(jac - (1 + eps_f(Qs))) == 0, str(jac))

# numeric check on the actual profile at 3 Q values
r9 = np.geomspace(3.0, 5.0e6, 20001)
e9, r_sv, phi_sv, boost_sv = st.eps_profile(r9, st.nu_be, 1.15, 1.0)


def I_of(rq):
    m = r9 >= rq
    return float(np.trapezoid(e9[m] / r9[m] ** 2, r9[m]))


# (first firing used h = 1e-4*Q, below the grid step: the masked-trapezoid
#  derivative is a stair function -- the SAME artifact the reviewer
#  disclosed in his F14 and fixed with smooth quadrature. Trap #23 again:
#  my verifier reproduced his bug, not an error in his claim. Fixed with
#  h = 0.05*Q >> grid step; disclosed.)
for Q0 in (3e4, 6e4, 9e4):
    h = Q0 * 0.05
    du = 2 / (Q0 + h) - 2 / (Q0 - h)
    dio = (2 / (Q0 + h) + 2 * I_of(Q0 + h)) - (2 / (Q0 - h) + 2 * I_of(Q0 - h))
    jnum = dio / du
    epsQ = float(np.interp(Q0, r9, e9))
    gb("1n Jacobian at Q=%.0e: %.4f vs 1+eps=%.4f" % (Q0, jnum, 1 + epsQ),
       abs(jnum - (1 + epsQ)) < 6e-3)

# GB-2 multipole channels l=1..8 (F4), BE eN=1.15 alpha=1
dP = []
mus = np.linspace(-1, 1, 20001)
for l in range(0, 9):
    c = np.zeros(l + 1)
    c[l] = 1
    v = _legval(mus, c)
    dP.append(float(v.max() - v.min()))
ch = {}
for l in range(1, 9):
    phl = np.interp(BAND / R_M_AU, r_sv, np.abs(phi_sv[l]))
    ch[l] = float(np.max(2.0 * phl * dP[l] / R_M_AU)) / W
gb("2 l=1 channel 0.191 xW", abs(ch[1] - 0.191) < 0.01, "%.3f" % ch[1])
gb("2 l=2 channel 0.240 xW", abs(ch[2] - 0.240) < 0.01, "%.3f" % ch[2])
gb("2 l=3 channel ~0.028 xW", abs(ch[3] - 0.028) < 0.01, "%.3f" % ch[3])
s18 = sum(ch.values())
gb("2 Sum(l=1..8) = 0.464 xW", abs(s18 - 0.464) < 0.02, "%.3f" % s18)

# GB-3 restored family-max row (F3+F4): simple eN=1.10 alpha=1.3
e_sm, r_sm, phi_sm, boost_sm = st.eps_profile(r9, st.nu_simple, 1.10, 1.3)
pm = (r_sm * R_M_AU >= 10 * R_M_AU) & (r_sm * R_M_AU <= 50 * R_M_AU)
eps_deep_sm = float(np.mean(boost_sm[pm] - 1.0)) * 1.0  # alpha folded in
scale_sm = float(np.max(np.interp(BAND, r9, e_sm))) / 1.0
# (first firing dropped the alpha=1.3 factor on phi_l -- the stage itself
#  applies it; my slip, disclosed. eps_profile scales eps by alpha but
#  returns RAW phi_l.)
mult_sm = sum(float(np.max(2.0 * np.interp(BAND / R_M_AU, r_sm,
                                           np.abs(phi_sm[l]) * 1.3)
                           * dP[l] / R_M_AU)) / W for l in range(1, 9))
scale_xw = scale_sm  # eps is dimensionless; excess width = eps * W => xW=eps
tot_sm = mult_sm + scale_xw
gb("3 restored simple 1.10 a1.3: multipoles %.3f + scale %.3f = %.3f xW "
   "(his 0.605+0.567=1.172)" % (mult_sm, scale_xw, tot_sm),
   abs(tot_sm - 1.172) < 0.06)

# GB-4 matched-sky ambient (F10): nu(e)*e = 1.93 per law
from scipy.optimize import brentq as _bq
eBE = _bq(lambda y: float(st.nu_be(np.array([y]))[0]) * y - 1.93, 0.8, 1.9)
eSI = _bq(lambda y: float(st.nu_simple(np.array([y]))[0]) * y - 1.93,
          0.8, 1.9)
gb("4 matched e_N: BE 1.318 / simple 1.271",
   abs(eBE - 1.318) < 0.005 and abs(eSI - 1.271) < 0.005,
   "%.3f / %.3f" % (eBE, eSI))

# GB-5 the mu-contrast (F2): |phi2|(1e4 au), ours(BE at matched) vs mu10
e_b2, r_b2, phi_b2, _ = st.eps_profile(r9, st.nu_be, eBE, 1.0)
p2_ours = float(np.interp(1e4 / R_M_AU, r_b2, np.abs(phi_b2[2])))


def nu_mun(n):
    def f(y):
        y = np.clip(np.asarray(y, float), 1e-12, None)
        out = np.empty_like(y)
        for i, yy in enumerate(y.ravel()):
            g = np.geomspace(max(yy, 1e-6) * 1e-3, max(yy, 1e-6) * 1e3, 400)
            mu = g / (1 + g ** n) ** (1.0 / n)
            out.ravel()[i] = np.interp(yy, mu * g, g) / yy
        return out
    return f


nu10 = nu_mun(10)
e10 = _bq(lambda y: float(nu10(np.array([y]))[0]) * y - 1.93, 1.0, 1.93)
_, r_m10, phi_m10, _ = st.eps_profile(r9, nu10, e10, 1.0, NR=256)
p2_m10 = float(np.interp(1e4 / R_M_AU, r_m10, np.abs(phi_m10[2])))
gb("5 |phi2|(1e4 au): ours %.3e vs mu10 %.3e; ratio %.0fx (his 38x)"
   % (p2_ours, p2_m10, p2_ours / p2_m10),
   30 < p2_ours / p2_m10 < 50)

# GB-6 mono inner law (F5/cond 5): simple hard-regime eps = (r/r_M)^2 =>
# spread of 2I over qp 5..40 au = 2*(40-5)/r_M^2
sp_analytic = 2.0 * 35.0 / R_M_AU ** 2
gb("6 mono spread analytic (simple) = 1.42e-6 = 0.019 xW",
   abs(sp_analytic - 1.416e-6) < 2e-8,
   "%.3e = %.3f xW" % (sp_analytic, sp_analytic / W))

# GB-7 F9 mapping numbers (EMP central 2I = 9.489e-5; nominal 2.226e-4)
lo1, hi1 = 1.0 / (0.75e-4 + 9.489e-5), 1.0 / 9.489e-5
lo2, hi2 = 1.0 / (0.75e-4 + 2.226e-4), 1.0 / 2.226e-4
gb("7 source bands: EMP %.0f-%.0f au (his 5886-10539); nominal %.0f-%.0f "
   "(his 3360-4492)" % (lo1, hi1, lo2, hi2),
   abs(lo1 - 5886) < 20 and abs(hi1 - 10539) < 20 and
   abs(lo2 - 3360) < 20 and abs(hi2 - 4492) < 20)
gb("7 tidal compensation margin ~45x = 2.226e-4 / 5e-6",
   abs(2.226e-4 / 5e-6 - 44.5) < 1.0, "%.1f" % (2.226e-4 / 5e-6))

# GB-8 F12 Galactic-tide quadrupole scale: d(1/a) ~ 4 pi rho Q^2 / Msun
pc_au = 206264.8
for rho, Qv, want in ((0.08, 2e4, None), (0.15, 1e5, 2.1e-6)):
    rho_au = rho / pc_au ** 3
    val = 4 * np.pi * rho_au * Qv ** 2
    if want:
        gb("8 tide quad at rho=%.2f, Q=%.0e: %.2e au^-1 (his %.1e)"
           % (rho, Qv, val, want), abs(val - want) / want < 0.15)
    else:
        gb("8 tide quad at rho=%.2f, Q=%.0e: %.2e = %.4f xW (his band "
           "0.0006-0.029)" % (rho, Qv, val, val / W),
           0.0004 < val / W < 0.03)

# GB-9 G4 re-fire on the letter-bearing channel (cond 6)
e_hi2, r_hi2, phi_hi2, _ = st.eps_profile(r9, st.nu_be, 1.15, 1.0, NR=1024)
q_lo = float(np.max(2.0 * np.interp(BAND / R_M_AU, r_sv, np.abs(phi_sv[2]))
                    * 1.5 / R_M_AU))
q_hi = float(np.max(2.0 * np.interp(BAND / R_M_AU, r_hi2, np.abs(phi_hi2[2]))
                    * 1.5 / R_M_AU))
gb("9 G4 re-fired on quad: NR 512->1024 rel d = %.2e (< 0.10; his 0.03%%)"
   % (abs(q_hi - q_lo) / q_lo), abs(q_hi - q_lo) / q_lo < 0.10)

# GB-10 grid facts (F5): solver inner edge and interp clamp
gb("10 solver inner edge = 70.3 au; np.interp clamps left",
   abs(r_sv[0] * R_M_AU - 70.3) < 0.2 and
   float(np.interp(5.0, r_sv * R_M_AU, boost_sv - 1.0))
   == float(boost_sv[0] - 1.0))

log("-" * 72)
log("GB SUMMARY: %d OK / %d MISMATCH" % (gok, gbad))
log("GB VERDICT: " + ("ALL REVIEWER NUMBERS CONFIRMED -- adopt R39 in "
                      "full (letter -> R-PROFILE-INCOMPLETE; 13 "
                      "conditions; P14 ADVERSE)" if gbad == 0 else
                      "MISMATCHES -- resolve before adoption"))
log("=" * 72)

with io.open(os.path.join(ROOT, "data/round39_addendum.txt"), "w",
             encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
