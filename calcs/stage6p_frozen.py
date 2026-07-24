"""
STAGE 6P (frozen-bath Test 2): the s-shaped per-system boost scatter
on the wide binaries. PRE-REGISTRATION COMMITTED BEFORE EXECUTION.

THE PREDICTION UNDER TEST: a frozen horizon bath gives each SYSTEM one
draw: the boost's occupation part scatters as delta_n/n =
sqrt((n+1)/n)/sqrt(N) with ONE standard normal per binary, so the
induced velocity smear GROWS with separation exactly like the signal
(deeper pairs, larger n) — distinguishable in principle from the
s-flat mass smear and the fenced contaminants the v7 model already
carries. Implementation (proxy-grade, disclosed): per-system factor
fz = sqrt[(1 + n_ch*max(1 + eps*sig, 0))/(1 + n_ch)] applied to the
dynamical velocities before companions/noise, with n_ch the ALPHA-
SCALED boost occupation at the system's characteristic field
g_N(a_s) and sig = sqrt((n_ch+1)/n_ch)/SQN, SQN = sqrt(N). The v7
machinery (3P on-disk config: BE law, g = 1.2 a0 tables, full
contaminant fences, catalog acceptance) is patched with exact-count
replacements; SQN = 1e12 recovers the unpatched model to machine
precision (fz = 1) — the paired baseline. Same seed => identical
population, orbits, and draws across the SQN grid: the contrast is
exactly paired.

SCAN: SQN grid = off (1e12), sqrt(60), sqrt(20), sqrt(8); seed 31;
BE law only (runtime; second seed = follow-up if a signal appears).
GATES: G1 baseline sanity — alpha_hat interior in [0.85, 1.35]
(on-disk-velocity config; 3S-era BE alpha_hat = 1.11); G2 patch
counts exact; G3 the SQN = 1e12 run's best lnL is the reference
(paired by construction).

PRE-REGISTERED OUTCOMES (Delta = best lnL(N) - best lnL(off), one
seed, so realization error ~ +-2-4 by 5K-era per-seed scatter):
 PREFERRED: Delta > +5 at some finite N with alpha_hat stable
   interior -> first binary-side evidence of a frozen component;
   check N-compatibility with the galaxies' point-level 20-60.
 STRIKE: Delta < -5 at N = 20 (the galaxy-favored value under the
   simplest same-N frozen reading) -> the binaries REJECT the
   frozen component at the strength the galaxy channel suggests;
   logged as a strike against that reading (assumption named: one
   frequency band per system, same N both systems).
 FLAT: |Delta| < 5 across the grid -> the binaries are insensitive
   at one-seed grade; quote the bound and stop (no credence move).
Writes data/stage6p_summary.txt + data/stage6p_frozen.txt.
"""
import os, re, sys, time
import numpy as np

SRC = open('calcs/stage3p_v7budget.py', encoding='utf-8-sig').read()

def patch(src, old, new, n=1):
    c = src.count(old)
    assert c == n, f"patch count {c} != {n} for: {old[:60]!r}"
    return src.replace(old, new)

src = SRC
src = patch(src, 'for law, TAB in (("simple", TAB_S), ("BE", TAB_B)):',
            'for law, TAB in (("BE", TAB_B),):')
src = patch(src, "p['gn1'], p['gn2'] = rng.normal(size=N), "
            "rng.normal(size=N)",
            "p['gn1'], p['gn2'] = rng.normal(size=N), "
            "rng.normal(size=N)\n    p['gn3'] = rng.normal(size=N)")
src = patch(src, "            tab_a = 1.0 + al*(TAB-1.0)",
            "            tab_a = 1.0 + al*(TAB-1.0)\n"
            "            _G['tab'] = tab_a")
src = patch(src, "        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])",
            "        vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])\n"
            "        gN_ch = GM*p['M_s'][idx]/np.maximum("
            "p['a_s'][idx]**2, 1e-12)\n"
            "        nch = np.interp(np.log(gN_ch/A0_CAN), LNY_U,"
            " _G['tab'], right=1.0) - 1.0\n"
            "        nch = np.maximum(nch, 0.0)\n"
            "        sigf = np.sqrt((nch+1.0)/np.maximum(nch, 1e-9))/SQN\n"
            "        occ = nch*np.maximum(1.0 + p['gn3'][idx]*sigf, 0.0)\n"
            "        fz = np.sqrt((1.0+occ)/(1.0+nch))")
src = patch(src, "            vp_b = vpar[idx].copy(); "
            "vq_b = vper[idx].copy()",
            "            vp_b = vpar[idx]*fz; vq_b = vper[idx]*fz")
src = patch(src, "    P(f\"  seed {seed}: BE-minus-simple best lnL = \"\n"
            "      f\"{best_lnl['BE']-best_lnl['simple']:+.1f}  \"\n"
            "      f\"({(time.time()-t0)/60:.1f} min)\")",
            "    P(f\"  seed {seed}: best lnL = {best_lnl['BE']:+.2f}  \"\n"
            "      f\"({(time.time()-t0)/60:.1f} min)\")")
src = patch(src, "'data/stage3u_summary.txt'",
            "'data/stage6p_summary.txt'", n=1)
src = patch(src, 'print("\\nbatch done; appended data/stage3u_summary.txt")',
            'print("batch done")')

SQNS = [1e12, 60.0**0.5, 20.0**0.5, 8.0**0.5]
TAGS = ['off', 'N=60', 'N=20', 'N=8']
for sqn, tag in zip(SQNS, TAGS):
    with open('data/stage6p_summary.txt', 'a') as f:
        f.write(f"--- 6P run {tag} (SQN={sqn:.4g}) ---\n")
    print(f"=== 6P {tag} ===", flush=True)
    sys.argv = ['stage6p', '1p2', '31']
    ns = {'SQN': sqn, '_G': {}, '__name__': '__main__'}
    t0 = time.time()
    exec(src, ns)
    print(f"=== {tag} done ({(time.time()-t0)/60:.1f} min) ===", flush=True)

# ---------------- verdict assembly (guarded) ----------------
L = ["STAGE 6P: frozen per-system boost scatter on the v7 binaries -- "
     "seed 31, BE law, g=1.2a0 tables, paired SQN grid", ""]
try:
    txt = open('data/stage6p_summary.txt').read()
    blocks = txt.split('--- 6P run ')[1:]
    res = {}
    for b in blocks:
        tag = b.split(' (SQN', 1)[0]
        m = re.search(r'seed 31 BE: a_hat=([0-9.]+) \(grid [0-9.]+, '
                      r'interior=(\w+)\), dlnL\(Newton\)=([+-][0-9.]+), '
                      r'wr=([0-9.]+)', b)
        mb = re.search(r'seed 31: best lnL = ([+-][0-9.]+)', b)
        if m and mb:
            res[tag] = (float(m.group(1)), m.group(2) == 'True',
                        float(m.group(3)), float(m.group(4)),
                        float(mb.group(1)))
    base = res['off'][4]
    g1 = 0.85 <= res['off'][0] <= 1.35 and res['off'][1]
    L.append(f"G1 baseline: a_hat = {res['off'][0]:.3f} interior="
             f"{res['off'][1]}, dlnL(Newton) = {res['off'][2]:+.1f}, "
             f"wr = {res['off'][3]} -> {'PASS' if g1 else 'FAIL'}")
    L.append("")
    L.append("   run    a_hat  interior  dlnL(Newt)   wr    best lnL"
             "    Delta vs off")
    for tag in TAGS:
        a, itr, dn, wr, bl = res[tag]
        L.append(f"  {tag:>5}  {a:5.3f}  {str(itr):>5}    {dn:+7.1f}"
                 f"   {wr:.2f}  {bl:+10.2f}   {bl-base:+7.2f}")
    d20 = res['N=20'][4] - base
    ds = {t: res[t][4] - base for t in TAGS[1:]}
    best_tag = max(ds, key=ds.get)
    if ds[best_tag] > 5.0 and res[best_tag][1]:
        gcomp = 'YES' if best_tag in ('N=20', 'N=60') else 'NO'
        call = (f"PREFERRED: Delta = {ds[best_tag]:+.1f} at {best_tag}"
                f" -- first binary-side frozen evidence; galaxy-N "
                f"compatibility: {gcomp}")
    elif d20 < -5.0:
        call = (f"STRIKE: Delta(N=20) = {d20:+.1f} -- the binaries "
                f"reject the galaxy-strength frozen component "
                f"(same-N reading)")
    else:
        call = (f"FLAT: max |Delta| = "
                f"{max(abs(v) for v in ds.values()):.1f} < 5 -- "
                f"insensitive at one-seed grade; bound quoted, no move")
    L.append("")
    L.append(f"PRE-REGISTERED VERDICT: {call}")
except Exception as ex:
    L.append(f"verdict INCOMPLETE: {type(ex).__name__}: {ex}")

out = "\n".join(L)
print(out)
with open('data/stage6p_frozen.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage6p_frozen.txt")
