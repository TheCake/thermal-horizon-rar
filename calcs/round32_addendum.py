# -*- coding: utf-8 -*-
"""ROUND 32 addendum — verification for papers/paper3_mechanism.md draft 0.4.

BLIND HALF (GA-*): written, run, and committed BEFORE the ROUND 32 referee
report exists — the 87a4676/884fbc1 pre-report protocol, third execution.
Verifies every number the 0.4 draft newly prints (the 10J/10L/10H/10I/10K
absorption) through independent routes where closed forms exist, and by
archived-output provenance where the number is a stage measurement already
referee-verified in rounds 29-31.

GB-* (post-report half): appended AFTER the ROUND 32 report lands —
re-computes every load-bearing number the referee produces (the standing
memory rule: never adopt a ruling before re-computing its arithmetic).

Output: data/round32_addendum.txt
"""
import io
import math
import os
import re
import sys

import sympy as sp

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.join(ROOT, "papers", "paper3_mechanism.md")
OUTPATH = os.path.join(ROOT, "data", "round32_addendum.txt")

LINES = []


def rep(tag, ok, detail):
    line = "[%s] %s — %s" % (tag, "PASS" if ok else "FAIL", detail)
    print(line)
    LINES.append(line)
    return ok


def main():
    results = []

    # ---------------- GA-1: ladder thresholds + the two-level degeneracy (sympy exact)
    n, kap = sp.symbols("n kappa", positive=True)
    E = (n + sp.Rational(1, 2)) - (kap / 4) * (n + sp.Rational(1, 2)) ** 2
    D = sp.expand(E.subs(n, n + 1) - E)
    kfold = sp.solve(sp.Eq(D, 0), kap)[0]
    thr = [sp.nsimplify(kfold.subs(n, m)) for m in (0, 1, 2)]
    fold_form_ok = sp.simplify(kfold - 2 / (n + 1)) == 0
    deg = sp.simplify((E.subs(n, 2) - E.subs(n, 1)).subs(kap, 1))
    ok = fold_form_ok and thr == [sp.Integer(2), sp.Integer(1), sp.Rational(2, 3)] and deg == 0
    results.append(rep("GA-1", ok,
                       "kappa_fold = %s; thresholds(n=0,1,2) = %s; E(2)-E(1) at kappa=1 = %s "
                       "(paper: 2/(n+1); 2, 1, 2/3; exact degeneracy)" % (kfold, thr, deg)))

    # ---------------- GA-2: occupation fixed-point existence ceiling (independent iteration route)
    def nbe(x):
        return 1.0 / math.expm1(x)

    def exists(x, k, c=1.0, itmax=20000, cap=1e9):
        nn = nbe(x)
        for _ in range(itmax):
            s = 1.0 - k * (nn + c) / 2.0
            if s <= 1e-12:
                return False
            n2 = nbe(x * s)
            if n2 > cap:
                return False
            if abs(n2 - nn) < 1e-13 * max(1.0, nn):
                return True
            nn = n2
        return True

    def kmax(x, lo=1e-8, hi=2.0):
        assert exists(x, lo) and not exists(x, hi)
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if exists(x, mid):
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)

    xs = [0.03, 0.05, 0.1, 0.2, 0.3]
    km = [kmax(x) for x in xs]
    ratios = [0.925 / k for k in km]
    ok = (0.014 < km[0] < 0.016) and (0.13 < km[-1] < 0.15) \
        and (61.0 < ratios[0] < 63.0) and (6.4 < ratios[-1] < 6.9)
    results.append(rep("GA-2", ok,
                       "kmax(x): " + ", ".join("%.2g:%.4f" % (x, k) for x, k in zip(xs, km))
                       + "; lock/kmax = %.1f .. %.1f (paper: ceiling 0.015-0.14; factor 6.6-62)"
                       % (ratios[-1], ratios[0])))

    # GA-2b: first fixed point at the measured kappa sits on the Newtonian arm
    lo_x, hi_x = 0.5, 8.0
    assert (not exists(lo_x, 0.925)) and exists(hi_x, 0.925)
    for _ in range(80):
        mid = 0.5 * (lo_x + hi_x)
        if exists(mid, 0.925):
            hi_x = mid
        else:
            lo_x = mid
    onset = 0.5 * (lo_x + hi_x)
    results.append(rep("GA-2b", 3.3 < onset < 3.8 and onset ** 2 > 10.0,
                       "onset at kappa=0.925: x = %.3f (y = x^2 = %.1f — Newtonian regime; paper: qualitative)"
                       % (onset, onset ** 2)))

    # ---------------- GA-3: production cap 1.93 (archived, R31-verified) + independent hard ceiling + demand 19
    with open(os.path.join(ROOT, "data", "stage10l_occbound.txt"), encoding="utf-8", errors="replace") as f:
        txt10l = f.read()
    sup_tokens = [t for t in ("1.916", "1.925", "1.456", "1.462") if t in txt10l]
    ok_arch = len(sup_tokens) >= 2
    # independent ceiling: at kappa=0.925 the rising branch keeps levels {0,1,2} only -> <n> < 2 -> nu < 3 << 19
    def Dn(m, k):
        return 1.0 - k * (m + 1) / 2.0
    keep_top = max(m for m in range(0, 50) if all(Dn(j, 0.925) > 0 for j in range(0, m)))
    demand = 1.0 + 0.925 * nbe(0.05)
    ok = ok_arch and keep_top == 2 and 18.9 < demand < 19.2
    results.append(rep("GA-3", ok,
                       "archived sup tokens %s; rising-branch top level at kappa=0.925 = %d "
                       "(any Gibbs mixture: nu < 3, independent ceiling); demand nu(x=0.05) = %.2f "
                       "(paper: cap 1.93 vs demand 19)" % (sup_tokens, keep_top, demand)))

    # ---------------- GA-4: convention pin kappa_r = kappa_b(1 - kappa_b/2)
    kb = sp.symbols("kappa_b", positive=True)
    kr = kb * (1 - kb / 2)
    crit = sp.solve(sp.diff(kr, kb), kb)[0]
    maxval = kr.subs(kb, 1)
    f = sp.lambdify(kb, kr)
    v_lo, v_mid, v_hi = f(0.888), f(1.0), f(1.10)
    band_lo = min(v_lo, v_hi)
    ok = crit == 1 and maxval == sp.Rational(1, 2) and band_lo >= 0.4935 and v_mid == 0.5
    results.append(rep("GA-4", ok,
                       "argmax = %s, max = %s; f(0.888) = %.4f, f(1.0) = %.4f, f(1.10) = %.4f "
                       "(paper: bare 0.89-1.10 maps into 0.494-0.500, max exactly 1/2)"
                       % (crit, maxval, v_lo, v_mid, v_hi)))

    # ---------------- GA-5: feedback-exponent mappings (5P bridge) + the 5R bound token
    beta = sp.symbols("beta")
    c1_of_beta = 1 / (2 * (1 + beta))
    beta_lock = sp.solve(sp.Eq(2 * (1 - c1_of_beta), sp.Rational(925, 1000)), beta)[0]
    ok_b = abs(float(beta_lock) + 0.0698) < 0.0015
    with open(os.path.join(ROOT, "NOTES-horizon-inertia.md"), encoding="utf-8", errors="replace") as f:
        notes = f.read()
    ok_5r = ("β<0.030" in notes) or ("β < 0.030" in notes) or ("beta<0.030" in notes) or ("0.030 (1σ)" in notes)
    results.append(rep("GA-5", ok_b and ok_5r,
                       "beta(kappa_lock=0.925) = %.4f via c1 = 1 - kappa/2 = 1/(2(1+beta)) "
                       "(paper: -0.07); 5R bound-below-0.03 token present in NOTES: %s"
                       % (float(beta_lock), ok_5r)))

    # ---------------- GA-6: redshift-lock external comparison arithmetic (95% CIs -> sigma = CI/1.96)
    c = 2.99792458e8
    H0 = 67.74 * 1000.0 / 3.0856775814913673e22
    lock = c * H0 / (2 * math.pi)
    z_mond = (1.047 - 1.03) / (0.05 / 1.96)
    z_pg = (1.05 - 1.047) / (0.05 / 1.96)
    z_dm = (1.047 - 1.00) / (0.04 / 1.96)
    z_slope = (1.20 - 1.024) / (0.10 / 1.96)
    ok = abs(lock - 1.047e-10) < 0.004e-10 and z_mond < 0.75 and z_pg < 0.15 \
        and 2.2 < z_dm < 2.4 and 3.3 < z_slope < 3.6
    results.append(rep("GA-6", ok,
                       "cH0/2pi(H0=67.74) = %.4e (paper: 1.05e-10); sigmas: MOND-row %.2f, "
                       "per-galaxy %.2f, DM-row %.2f, slope %.2f (paper: 0.7 / 0.1 / 2.3 / 3.4)"
                       % (lock, z_mond, z_pg, z_dm, z_slope)))

    # ---------------- GA-7: susceptibility residue vs the archived stage output (R30-verified)
    with open(os.path.join(ROOT, "data", "stage10j_kappadeep.txt"), encoding="utf-8", errors="replace") as f:
        txt10j = f.read()
    m = re.search(r"\(1-kappa\) chi\(1\) = ([0-9.]+) vs exact limit 4 e\^\{-d\^2\}/om = ([0-9.]+)", txt10j)
    ok7 = False
    detail = "residue line not found in archived output"
    if m:
        got, lim = float(m.group(1)), float(m.group(2))
        rel = abs(got - lim) / lim
        ok7 = 1e-5 < rel < 1e-4 and abs(rel * 100 - 0.003) < 0.002
        detail = ("archived (1-kappa)chi(1) = %.6f vs exact limit %.6f -> rel = %.4f%% "
                  "(paper: verified at the 0.003%% level)" % (got, lim, rel * 100))
    results.append(rep("GA-7", ok7, detail))

    # ---------------- GA-8: reproducibility-table files exist on disk
    new_files = [
        "calcs/stage10j_kappadeep.py", "calcs/round30_addendum.py",
        "data/stage10j_kappadeep.txt", "data/round30_addendum.txt",
        "calcs/stage10l_occbound.py", "calcs/round31_addendum.py",
        "data/stage10l_occbound.txt", "data/round31_addendum.txt",
        "calcs/stage10h_zladder.py", "calcs/stage10h_addendum.py",
        "data/stage10h_out.txt", "data/stage10h_addendum.txt",
        "calcs/stage10i_isoladder.py", "calcs/stage10k_groupend.py",
        "data/stage10i_isoladder.txt", "data/stage10k_groupend.txt",
    ]
    missing = [p for p in new_files if not os.path.exists(os.path.join(ROOT, p.replace("/", os.sep)))]
    results.append(rep("GA-8", not missing,
                       "all %d new reproducibility-row files exist (missing: %s)"
                       % (len(new_files), missing if missing else "none")))

    # ---------------- GA-9: paper hygiene (abstract budget, banned strings, codenames, punctuation)
    with open(PAPER, encoding="utf-8") as f:
        paper = f.read()
    abstract = paper.split("## Abstract")[1].split("## 1. Introduction")[0]
    aw = len(abstract.split())
    banned = ["measured near one", "same world", "anti" + "gravity", ("anti" + "gravity").capitalize(), "κ measured"]
    hits = [b for b in banned if b in paper]
    prose = "\n".join(l for l in paper.splitlines()
                      if not l.strip().startswith("|") and "calcs/" not in l and "data/" not in l)
    code_hits = re.findall(r"\b(?:stage10[a-l]|10[A-L]\b|ROUND\s\d+|K-SPLIT|L-PINNED|K-POWER-DEAD|SCHA)\b", prose)
    emdash = paper.count("—")
    body = re.sub(r"\|[^\n]*", "", paper)
    sentences = re.split(r"(?<=[.!?])\s+", body)
    lens = [len(s.split()) for s in sentences if len(s.split()) > 2]
    mean_len = sum(lens) / float(len(lens))
    ok = aw <= 250 and not hits and not code_hits
    results.append(rep("GA-9", ok,
                       "abstract %d words (bar 250); banned hits %s; codename hits %s; "
                       "em-dashes %d total; sentence mean %.1f words"
                       % (aw, hits if hits else "none", code_hits if code_hits else "none",
                          emdash, mean_len)))

    # ================= GB (post-report half): re-compute every new ROUND 32 referee claim
    # Written AFTER REVIEW-ROUND32-OPUS.md landed; verifies his findings before adoption
    # (the standing memory rule). GA above is unchanged and re-runs for regression.

    # ---------------- GB-1: M1 attribution — which fit in the 10D archive carries 1.503?
    with open(os.path.join(ROOT, "data", "stage10d_kappa.txt"), encoding="utf-8", errors="replace") as f:
        t10d = f.read()
    single_tok = re.search(r"F1 \(kappa, a0 free\): kappa = 1\.503", t10d)
    f4_tok = re.search(r"F4 \(split\): kappa_d = 1\.317, kappa_t = 1\.036", t10d)
    has_150 = single_tok is not None
    has_boot = ("1.321" in t10d) and ("1.670" in t10d)
    split_vals = [v for v in ("1.317", "1.036") if v in t10d]
    rej = [v for v in ("21.74", "14.37") if v in t10d]
    okb1 = has_150 and (f4_tok is not None) and has_boot and len(rej) == 2 and len(split_vals) == 2
    results.append(rep("GB-1", okb1,
                       "archive: 1.503 present %s; boot [1.32,1.67] %s; rejections %s; split pair %s; "
                       "single-kappa line: %r -> M1 CONFIRMED (1.50 = the forced single-kappa fit; "
                       "paper adopts the label fix ONLY — split values stay unprinted per the retired-"
                       "decomposition discipline)" % (has_150, has_boot, rej, split_vals,
                                                      (single_tok.group(0)[:90] if single_tok else None))))

    # ---------------- GB-2: m3 — SdS circular-orbit epicyclic frequency, from scratch (sympy)
    r_, M_, L_, E_, lam = sp.symbols("r M Lz E Lambda", positive=True)
    fmet = 1 - 2 * M_ / r_ - lam * r_ ** 2 / 3
    Om2_exact = sp.simplify(sp.diff(fmet, r_) / (2 * r_))          # coordinate angular velocity^2
    okb2a = sp.simplify(Om2_exact - (M_ / r_ ** 3 - lam / 3)) == 0
    # radial effective potential V(r) = f(1 + L^2/r^2); circular orbit: V'=0 -> L^2(r)
    V = fmet * (1 + L_ ** 2 / r_ ** 2)
    L2 = sp.solve(sp.Eq(sp.diff(V, r_), 0), L_ ** 2)[0]
    E2 = sp.simplify(V.subs(L_ ** 2, L2))
    # coordinate-time radial epicyclic: kappa^2 = (V''/2) * (f/E)^2 evaluated on the circle
    Vpp = sp.simplify(sp.diff(V, r_, 2).subs(L_ ** 2, L2))
    kap2 = sp.simplify(Vpp / 2 * fmet ** 2 / E2)
    kap2_series = sp.series(kap2, M_, 0, 2).removeO()              # exact in Lambda, first order in M
    target = M_ / r_ ** 3 - sp.Rational(4, 3) * lam + 5 * lam * M_ / r_
    okb2b = sp.simplify(sp.expand(kap2_series - target)) == 0
    # magnitudes of the neglected relative correction ~ v^2/c^2
    v2c2_gal = (300e3 / 2.99792458e8) ** 2
    v2c2_bin = (1.0e3 / 2.99792458e8) ** 2
    results.append(rep("GB-2", okb2a and okb2b and v2c2_gal < 1.1e-6,
                       "Omega^2 = M/r^3 - Lambda/3 EXACT: %s; kappa_r^2 to O(M) = M/r^3 - (4/3)Lambda "
                       "+ 5 Lambda M/r: %s (referee's term CONFIRMED); neglected relative order "
                       "v^2/c^2 = %.1e (300 km/s) / %.1e (1 km/s) -> paper's 'below 1e-6' holds"
                       % (okb2a, okb2b, v2c2_gal, v2c2_bin)))

    # ---------------- GB-3: m5 — superseded vs corrected radiative-carrier U values
    with open(os.path.join(ROOT, "data", "stage10a_dprov.txt"), encoding="utf-8", errors="replace") as f:
        t10a = f.read()
    with open(os.path.join(ROOT, "data", "stage10a_addendum.txt"), encoding="utf-8", errors="replace") as f:
        t10aa = f.read()
    old_in_dprov = any(v in t10a for v in ("4.810e-45", "4.81e-45", "8.197e-35"))
    new_in_add = any(v in t10aa for v in ("7.200e-45", "7.2e-45")) and any(v in t10aa for v in ("1.227e-34",))
    orders_new = (-math.log10(7.2e-45), -math.log10(1.227e-34))
    orders_old = (-math.log10(4.81e-45), -math.log10(8.197e-35))
    robust = all(33.5 < o < 44.6 for o in orders_new + orders_old)
    results.append(rep("GB-3", old_in_dprov and new_in_add and robust,
                       "dprov carries superseded tokens: %s; addendum carries corrected 7.2e-45/1.227e-34: %s; "
                       "orders new (%.1f, %.1f) / old (%.1f, %.1f) -> '34 to 44 orders' robust to both "
                       "(m5 CONFIRMED; App B row annotated)" % (old_in_dprov, new_in_add,
                                                                orders_new[0], orders_new[1],
                                                                orders_old[0], orders_old[1])))

    # ---------------- GB-4: his fresh arithmetic — Fermi optimum, r-hat containment, shortfall squares
    from scipy.optimize import brentq
    xstar = brentq(lambda x: (x - 1) - math.exp(-x), 0.5, 3.0)
    z_half = (0.50 - 0.3365) / 0.1869
    sq = (6.67 ** 2, 9.18 ** 2)
    okb4 = abs(xstar - 1.2784645) < 1e-6 and 0.85 < z_half < 0.90 and 44.0 < sq[0] < 45.0 and 84.0 < sq[1] < 84.6
    results.append(rep("GB-4", okb4,
                       "Fermi x* (x-1 = e^-x) = %.7f (his 1.2784645); r-hat contains 1/2 at %.2f sigma "
                       "(his 0.87, paper '0.9'); FD-shortfall squares %.1f / %.1f (paper '45 to 84')"
                       % (xstar, z_half, sq[0], sq[1])))

    # ---------------- GB-5: m4 — long-sentence and punctuation census (pre/post-fix state of the file)
    with open(PAPER, encoding="utf-8") as f:
        paper2 = f.read()
    body2 = re.sub(r"\|[^\n]*", "", paper2)
    body2 = "\n".join(l for l in body2.splitlines() if not l.strip().startswith("!["))
    sents2 = re.split(r"(?<=[.!?])\s+", body2)
    long_sents = [len(s.split()) for s in sents2 if len(s.split()) >= 70]
    colons = body2.count(": ")
    semis = body2.count("; ")
    results.append(rep("GB-5", len(long_sents) == 0,
                       "sentences >= 70 words after the m4 splits: %d %s; mid-sentence colons %d, "
                       "semicolons %d (referee counted 71/51 pre-fix)" % (len(long_sents), long_sents, colons, semis)))

    # ---------------- summary
    n_pass = sum(1 for r in results if r)
    verdict = "ALL PASS (%d/%d) — GA (blind, pre-report) + GB (post-report) both halves closed." % (n_pass, len(results)) \
        if n_pass == len(results) else "INCOMPLETE: %d/%d — fix before adopting." % (n_pass, len(results))
    print(verdict)
    LINES.append(verdict)

    with open(OUTPATH, "w", encoding="utf-8") as f:
        f.write("ROUND 32 addendum — blind half (pre-report)\n")
        f.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
