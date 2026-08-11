"""ROUND 35 addendum -- verification for stage 10P (P13-CROWD).

GA half (BLIND, written and committed BEFORE the ROUND 35 report is read;
the 87a4676 protocol, 7th execution): independent re-derivations of the
stage's load-bearing numbers with FRESH code (no import from the stage
script). Each check loads the stage outputs only to COMPARE, never to
compute.

GB half (post-report): re-computation of every load-bearing number the
ROUND 35 reviewer produces (memory rule feedback-verify-reviewer-math);
appended after the report lands, results printed with GB- prefixes.

Run: py calcs/round35_addendum.py --ga   (blind half)
     py calcs/round35_addendum.py --gb   (after the report; edits below)
"""
import argparse
import json
import os

import numpy as np
from astropy.io import fits

BASE = os.path.join(os.path.dirname(__file__), "..", "data", "xcop")
OUTD = os.path.join(os.path.dirname(__file__), "..", "data")

# fresh constants (independent arithmetic)
C = 2.99792458e8
MPC = 3.0856775814913673e22
KPC = MPC / 1e3
G = 6.674e-11
MSUN = 1.989e30
A0_CHECK = C * (67.8e3 / MPC) / (2 * np.pi)

T2 = {  # Eckert+22 Table 2 (independent retype from the primary read)
    "A85": (0.0555, 1214.0), "A644": (0.0704, 1398.0), "A1644": (0.0473, 1031.0),
    "A1795": (0.0622, 1160.0), "A2029": (0.0766, 1340.0), "A2142": (0.0900, 1453.0),
    "A2255": (0.0809, 1202.0), "A2319": (0.0557, 1424.0), "A3158": (0.0590, 1146.0),
    "A3266": (0.0589, 1381.0), "RXC1825": (0.0650, 1109.0), "ZW1215": (0.0766, 1346.0),
}


def nbe(x):
    return 1.0 / np.expm1(x)


def fresh_m2ll(tab, model_log10gobs, slope_fun, sig_int):
    """Independent -2lnL implementation (log-space Gaussian, model-slope
    propagated g_bar errors) -- structurally re-derived from the header."""
    lg_model = model_log10gobs(tab["gbar"])
    sl = slope_fun(tab["gbar"])
    var = tab["s_logo"] ** 2 + (sl * tab["s_logb"]) ** 2 + sig_int ** 2
    r = np.log10(tab["gobs"]) - lg_model
    return float(np.sum(r * r / var + np.log(2 * np.pi * var)))


def num_slope(f, g):
    d = 0.01
    return (np.log10(f(g * 10 ** d)) - np.log10(f(g * 10 ** -d))) / (2 * d)


def load_tab(z, prefix):
    return {k: z[f"{prefix}_{k}"] for k in ("name", "r", "gobs", "gbar",
                                            "s_logo", "s_logb")}


def ga():
    print("=== ROUND 35 GA (blind) ===")
    res = json.load(open(os.path.join(OUTD, "stage10p_results.json")))
    gat = json.load(open(os.path.join(OUTD, "stage10p_gates.json")))
    z = np.load(os.path.join(OUTD, "stage10p_sky.npz"), allow_pickle=True)
    tab = load_tab(z, "p")

    # GA-1: a0 arithmetic + deep-limit identity
    print(f"GA-1 a0: fresh {A0_CHECK:.6e} vs stage {res['a0']:.6e} "
          f"d={abs(A0_CHECK - res['a0']):.2e}")
    Mhat = res["primary"]["CROWD"]["M"]
    print(f"GA-1 gdagger: M^2*a0 = {Mhat**2 * A0_CHECK:.4e} vs stage "
          f"{res['gdagger_eff']:.4e}")

    # GA-2: fresh NFW rebuild at the release convention for 3 clusters
    H0R = 70e3 / MPC
    for name in ("A85", "A2142", "RXC1825"):
        zc, _ = T2[name]
        with fits.open(os.path.join(BASE, name, f"{name}_hydro_mass.fits")) as h:
            t = h["HYDRO_MASS"].data
            rr = np.array(t["RADIUS"], float)
            mn = np.array(t["M_NFW"], float)
            p = h["PARAMS"].data
            for row in p:
                if str(row["MODEL"]).strip() == "NFW":
                    rs, c200 = float(row["RS"]), float(row["C200"])
        hz2 = H0R ** 2 * (0.3 * (1 + zc) ** 3 + 0.7)
        rhoc = 3 * hz2 / (8 * np.pi * G)
        rho_s = (200 / 3) * rhoc * c200 ** 3 / (np.log(1 + c200) - c200 / (1 + c200))
        xx = rr / rs
        m_us = 4 * np.pi * rho_s * (rs * KPC) ** 3 * (np.log(1 + xx) - xx / (1 + xx)) / MSUN
        med = np.median(np.abs(np.log10(m_us) - np.log10(mn)))
        print(f"GA-2 {name}: fresh NFW median |dlog M| = {med:.4f} "
              f"(stage G3 bar 0.02) {'OK' if med <= 0.02 else 'FAIL'}")

    # GA-3: fresh -2lnL of the G10 probe + the reported best fits
    probe = fresh_m2ll(tab, lambda g: np.log10(g * (1 + nbe(np.sqrt(g / A0_CHECK)))),
                       lambda g: num_slope(lambda gg: gg * (1 + nbe(np.sqrt(gg / A0_CHECK))), g),
                       0.1)
    print(f"GA-3 probe: fresh {probe:.6f} vs stored {gat['G10_probe_build']:.6f} "
          f"d={abs(probe - gat['G10_probe_build']):.2e}")
    for key, fun in (
            ("CROWD", lambda g, M=res["primary"]["CROWD"]["M"]:
             g * (1 + M * nbe(np.sqrt(g / A0_CHECK)))),
            ("SCALE", lambda g, S=res["primary"]["SCALE"]["S"]:
             g * (1 + nbe(np.sqrt(g / (A0_CHECK * S)))))):
        m2 = fresh_m2ll(tab, lambda g: np.log10(fun(g)),
                        lambda g: num_slope(fun, g),
                        res["primary"][key]["sig"])
        print(f"GA-3 {key}: fresh -2lnL at reported params = {m2:.3f} vs "
              f"reported {res['primary'][key]['m2']:.3f} "
              f"d={m2 - res['primary'][key]['m2']:+.4f}")
    # MODEL-B at reported per-cluster U
    U = res["primary"]["B_U"]

    def bfun_all(g):
        out = np.empty_like(g)
        for nm, u in U.items():
            m = tab["name"] == nm
            ge = u * g[m]
            out[m] = ge * (1 + nbe(np.sqrt(ge / A0_CHECK)))
        return out
    lgm = np.log10(bfun_all(tab["gbar"]))
    sl = np.empty_like(tab["gbar"])
    for nm, u in U.items():
        m = tab["name"] == nm
        f = lambda g, uu=u: (uu * g) * (1 + nbe(np.sqrt(uu * g / A0_CHECK)))
        sl[m] = num_slope(f, tab["gbar"][m])
    var = tab["s_logo"] ** 2 + (sl * tab["s_logb"]) ** 2 + res["primary"]["B"]["sig"] ** 2
    r = np.log10(tab["gobs"]) - lgm
    m2B = float(np.sum(r * r / var + np.log(2 * np.pi * var)))
    print(f"GA-3 B: fresh -2lnL at reported params = {m2B:.3f} vs reported "
          f"{res['primary']['B']['m2']:.3f} d={m2B - res['primary']['B']['m2']:+.4f}")

    # GA-4: hand assembly of A2029 at 3 radii (fresh interpolation code)
    name = "A2029"
    r500rel = dict(zip([str(s) for s in z["r500rel_names"]],
                       [float(v) for v in z["r500rel_vals"]]))[name]
    with fits.open(os.path.join(BASE, name, f"{name}_hydro_mass.fits")) as h:
        t = h["HYDRO_MASS"].data
        rh, mf = np.array(t["RADIUS"], float), np.array(t["M_FORW"], float)
    with fits.open(os.path.join(BASE, name, f"{name}_fgas_profile.fits")) as h:
        t = h["FGAS"].data
        rg = np.array(t["RADIUS"], float) * r500rel
        mg = np.array(t["MGAS"], float)
    with fits.open(os.path.join(BASE, name, f"{name}_mstar.fits")) as h:
        t = h["MSTAR_SMOOTHED"].data
        rs_, ms_ = np.array(t["RADIUS"], float), np.array(t["MSTAR"], float)

    def li(r, rr, mm):
        return 10 ** np.interp(np.log10(r), np.log10(rr), np.log10(mm))
    sel = tab["name"] == name
    for rq in (tab["r"][sel][0], tab["r"][sel][3], tab["r"][sel][-1]):
        go = G * li(rq, rh, mf) * MSUN / (rq * KPC) ** 2
        gb = G * (li(rq, rg, mg) + li(rq, rs_, ms_)) * MSUN / (rq * KPC) ** 2
        j = np.argmin(np.abs(tab["r"][sel] - rq))
        print(f"GA-4 A2029 r={rq:7.1f} kpc: gobs fresh {go:.4e} vs table "
              f"{tab['gobs'][sel][j]:.4e} (d={abs(np.log10(go / tab['gobs'][sel][j])):.1e} dex); "
              f"gbar fresh {gb:.4e} vs {tab['gbar'][sel][j]:.4e} "
              f"(d={abs(np.log10(gb / tab['gbar'][sel][j])):.1e} dex)")

    # GA-5: BIC/AIC arithmetic from reported m2/k/n
    n = res["primary"]["CROWD"].get("n", None) or len(tab["gbar"])
    for key in ("GAL", "CROWD", "SCALE", "B"):
        k = res["primary"][key]["k"]
        m2 = res["primary"][key]["m2"]
        print(f"GA-5 {key}: BIC fresh {m2 + k * np.log(len(tab['gbar'])):.3f} vs "
              f"reported {res['primary'][key]['bic']:.3f}")
    print(f"GA-5 dBIC_CB fresh "
          f"{res['primary']['CROWD']['bic'] - res['primary']['B']['bic']:+.3f} "
          f"vs reported {res['dBIC_CB']:+.3f}")
    print("=== GA done ===")


def gb():
    print("=== ROUND 35 GB (post-report) ===")
    print("(filled after the ROUND 35 report lands)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ga", action="store_true")
    ap.add_argument("--gb", action="store_true")
    a = ap.parse_args()
    if a.ga:
        ga()
    if a.gb:
        gb()
