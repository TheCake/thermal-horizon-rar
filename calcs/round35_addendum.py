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
    """Re-compute every load-bearing number the ROUND 35 reviewer produced
    (REVIEW-ROUND35-OPUS.md; memory rule feedback-verify-reviewer-math).
    His claims checked here:
      GB-1 within-cluster lag-1 autocorrelation of F-CROWD residuals
           rho ~ 0.538 (range 0.42-0.61) => n_eff ~ 17 under AR(1).
      GB-2 in-window curve separation nu_CROWD(M_hat) vs nu_SCALE(S_hat)
           <= 9.6 pct over the observed g_bar range; ~25 pct at x = 2.
      GB-3 per-x decomposition of the SCALE-over-CROWD -2lnL gap:
           +4.69 from x < 0.5, +3.19 from x >= 0.5 (sum 7.88).
      GB-4 raw -2lnL gaps and 6*ln(n) at N_R = 6/8/10:
           +20.77/+28.53/+36.13 vs 22.43/24.15/25.49.
      GB-5 per-N_R cluster bootstrap (his 300 reps, independent seed):
           P(B better) = 0.367/0.880/1.000; means -1.22/+4.84/+11.12;
           P(dBIC > +10) = 0.003/0.110/0.550. Ours: 150 reps, fresh
           seed, B at 2 starts (distribution-grade check, +-0.04 on P).
      GB-6 run-1 G3 offsets: mean 0.0274, range 0.0264-0.0285 vs
           2*log10(70/67.8) = 0.02774; h-seam -0.0139 dex; B U-spread
           0.047 dex half-range.
    """
    print("=== ROUND 35 GB (post-report) ===")
    import sys
    sys.path.insert(0, os.path.dirname(__file__))
    import importlib
    st = importlib.import_module("stage10p_p13crowd")
    res = json.load(open(os.path.join(OUTD, "stage10p_results.json")))
    z = np.load(os.path.join(OUTD, "stage10p_sky.npz"), allow_pickle=True)
    rng = np.random.default_rng(777001)

    # GB-1 residual autocorrelation
    tab = load_tab(z, "p")
    Mhat = res["primary"]["CROWD"]["M"]
    resid = np.log10(tab["gobs"]) - np.log10(st.gobs_crowd(tab["gbar"], Mhat))
    rhos_raw, rhos_dm = [], []
    for nm in sorted(set(tab["name"].tolist())):
        e = resid[tab["name"] == nm]
        rhos_raw.append(float(np.corrcoef(e[:-1], e[1:])[0, 1]))
        ed = e - e.mean()
        rhos_dm.append(float(np.corrcoef(ed[:-1], ed[1:])[0, 1]))
    rho_r, rho_d = float(np.mean(rhos_raw)), float(np.mean(rhos_dm))
    neff_d = 56 * (1 - rho_d) / (1 + rho_d)
    print(f"GB-1 raw rho {rho_r:.3f} [{min(rhos_raw):.2f}, {max(rhos_raw):.2f}]"
          f" (offset-dominated; ~1 eff pt/cluster => n_eff ~ 7);")
    print(f"GB-1 de-meaned rho {rho_d:.3f} [{min(rhos_dm):.2f}, {max(rhos_dm):.2f}]"
          f" (his 0.538 [0.42, 0.61]); AR1 n_eff {neff_d:.1f} (his ~17)")

    # GB-2 curve separation
    gg = np.logspace(np.log10(tab["gbar"].min()), np.log10(tab["gbar"].max()), 400)
    Shat = res["primary"]["SCALE"]["S"]
    r_in = np.abs(st.gobs_crowd(gg, Mhat) / st.gobs_scale(gg, Shat) - 1.0)
    g2 = st.A0 * 4.0   # x = 2
    r_x2 = float(np.abs(st.gobs_crowd(np.array([g2]), Mhat)
                        / st.gobs_scale(np.array([g2]), Shat) - 1.0)[0])
    print(f"GB-2 in-window max sep {100*r_in.max():.2f} pct (his 9.6); "
          f"at x=2: {100*r_x2:.1f} pct (his ~25)")

    # GB-3 per-x decomposition
    def perpoint(fun, pars, sig):
        gm = fun(tab["gbar"], *pars)
        sl = st.model_slope(fun, tab["gbar"], *pars)
        var = tab["s_logo"] ** 2 + (sl * tab["s_logb"]) ** 2 + sig ** 2
        r = np.log10(tab["gobs"]) - np.log10(gm)
        return r * r / var + np.log(2 * np.pi * var)
    pc = perpoint(st.gobs_crowd, [Mhat], res["primary"]["CROWD"]["sig"])
    ps = perpoint(st.gobs_scale, [Shat], res["primary"]["SCALE"]["sig"])
    x = np.sqrt(tab["gbar"] / st.A0)
    lo = float(np.sum((pc - ps)[x < 0.5]))
    hi = float(np.sum((pc - ps)[x >= 0.5]))
    print(f"GB-3 dm2 split: x<0.5 {lo:+.2f} (his +4.69), x>=0.5 {hi:+.2f} "
          f"(his +3.19), sum {lo+hi:+.2f} (his 7.88)")

    # GB-4 raw gaps at N_R legs (refit)
    for pref, nr, his_raw in (("n6", 42, 20.77), ("p", 56, 28.53), ("n10", 70, 36.13)):
        tt = load_tab(z, pref)
        fs = st.fit_suite(tt, np.random.default_rng(5150 + nr))
        raw = fs["CROWD"]["m2"] - fs["B"]["m2"]
        print(f"GB-4 N={nr}: raw {raw:+.2f} (his {his_raw:+.2f}); "
              f"6ln(n) = {6*np.log(fs['n']):.2f}")

    # GB-5 per-N_R bootstrap (150 reps, B at 2 starts)
    for pref, his_p, his_mean in (("n6", 0.367, -1.22), ("p", 0.880, 4.84),
                                  ("n10", 1.000, 11.12)):
        tt = load_tab(z, pref)
        names = sorted(set(tt["name"].tolist()))
        ds = []
        for i in range(150):
            pick = rng.choice(names, size=len(names), replace=True)
            idx = np.concatenate([np.where(tt["name"] == nm)[0] for nm in pick])
            tb = {k: v[idx] for k, v in tt.items()}
            tb["name"] = np.concatenate(
                [[f"{nm}#{j}"] * int((tt["name"] == nm).sum())
                 for j, nm in enumerate(pick)])
            fb = st.fit_suite(tb, rng)
            ds.append(fb["CROWD"]["bic"] - fb["B"]["bic"])
        a = np.array(ds)
        print(f"GB-5 {pref}: P(B better) {np.mean(a > 0):.3f} (his {his_p}); "
              f"mean {a.mean():+.2f} (his {his_mean:+.2f}); "
              f"P(>+10) {np.mean(a > 10):.3f}")

    # GB-6 run-1 G3 offsets + arithmetic
    g1 = json.load(open(os.path.join(OUTD, "stage10p_gates.json")))
    # run-1 values are in the preserved log; re-derive fresh at 67.8:
    offs = []
    H067 = 67.8e3 / MPC
    for name in T2:
        zc, _ = T2[name]
        with fits.open(os.path.join(BASE, name, f"{name}_hydro_mass.fits")) as h:
            t = h["HYDRO_MASS"].data
            rr = np.array(t["RADIUS"], float)
            mn = np.array(t["M_NFW"], float)
            for row in h["PARAMS"].data:
                if str(row["MODEL"]).strip() == "NFW":
                    rs, c200 = float(row["RS"]), float(row["C200"])
        hz2 = H067 ** 2 * (0.308 * (1 + zc) ** 3 + 0.692)
        rhoc = 3 * hz2 / (8 * np.pi * G)
        rho_s = (200 / 3) * rhoc * c200 ** 3 / (np.log(1 + c200) - c200 / (1 + c200))
        xx = rr / rs
        m_us = 4 * np.pi * rho_s * (rs * KPC) ** 3 * (np.log(1 + xx) - xx / (1 + xx)) / MSUN
        offs.append(float(np.median(np.abs(np.log10(m_us) - np.log10(mn)))))
    print(f"GB-6 G3-at-67.8 offsets: mean {np.mean(offs):.4f} range "
          f"[{min(offs):.4f}, {max(offs):.4f}] (his 0.0274 [0.0264, 0.0285]); "
          f"2log10(70/67.8) = {2*np.log10(70/67.8):.5f}")
    U = res["primary"]["B_U"]
    spread = 0.5 * (np.log10(max(U.values())) - np.log10(min(U.values())))
    print(f"GB-6 h-seam {np.log10(67.8/70):.4f} dex; U half-spread "
          f"{spread:.3f} dex (his 0.047)")
    print("=== GB done ===")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--ga", action="store_true")
    ap.add_argument("--gb", action="store_true")
    a = ap.parse_args()
    if a.ga:
        ga()
    if a.gb:
        gb()
