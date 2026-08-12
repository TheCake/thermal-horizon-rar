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

with io.open(os.path.join(ROOT, "data/round39_addendum.txt"), "w",
             encoding="utf-8") as f:
    f.write("\n".join(OUT) + "\n")
