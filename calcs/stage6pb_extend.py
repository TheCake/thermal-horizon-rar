"""
STAGE 6P-b (the edge-and-specificity extension). PRE-REGISTRATION
COMMITTED BEFORE EXECUTION.

6P fired its PREFERRED branch at the GRID EDGE (Delta = +11.1 at N=8,
monotone, no interior optimum). Per this program's own edge-flagging
discipline (correction #4 lineage), an edge preference is NOT a
measurement, and the v7 history (the sigma_m saga: fitted broadening
persistently exceeding measured mass errors) supplies an obvious
mundane absorber: ANY added scatter channel may soak residual
broadening. This stage supplies the three missing legs:

 (1) GRID EXTENSION: shaped N = 4, 2 -- find the turnover (or ride
     the edge again, which itself is diagnostic: N < 8 means
     per-system velocity scatter > 10%, into generic-broadening
     territory).
 (2) THE s-FLAT CONTROL (decisive): same per-system Gaussian draw,
     s-INDEPENDENT amplitude matched to the shaped model's mean over
     wide systems (a_s in 6-50 kAU -- the GENEROUS convention: the
     flat control gets MORE scatter than shaped in the close bins,
     so a surviving shaped-over-flat lead is conservative;
     disclosed). Run at matched N=8 and N=20 strengths.
 (3) REALIZATION: seed 101 at off and shaped N=8.

PRE-REGISTERED BANDS:
 SHAPE-SPECIFIC: shaped(N=8) - flat(matched-8) >= +5 on seed 31 AND
   seed-101 shaped Delta sign-consistent -> the s-shape itself is
   preferred: a genuine frozen-signature candidate; N-hat then read
   from the (extended) profile if interior.
 GENERIC: |shaped - flat| < 3 at matched strength -> the improvement
   is s-blind extra broadening = the old v7 residual finding a new
   absorber; NOT frozen evidence; logged as an error-model finding
   (real, useful, not mechanism support). No credence move.
 AMBIG: between bands, or seed-101 sign flip -> one-seed noise;
   park pending more seeds.
Appends data/stage6p_summary.txt; writes data/stage6pb_extend.txt.
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
src = patch(src, "    s_kau = smag/1e3",
            "    s_kau = smag/1e3\n"
            "    gN_ch = GM*p['M_s']/np.maximum(p['a_s']**2, 1e-12)\n"
            "    nch = np.maximum(np.interp(np.log(gN_ch/A0_CAN), LNY_U,"
            " _G['tab'], right=1.0) - 1.0, 0.0)\n"
            "    if MODE == 'flat':\n"
            "        wide = (p['a_s'] >= 6e3) & (p['a_s'] <= 5e4)\n"
            "        sbar = np.mean(0.5*np.sqrt(nch[wide]/(nch[wide]+1.0))"
            ")/SQN\n"
            "        fz_full = np.maximum(1.0 + p['gn3']*sbar, 0.05)\n"
            "    else:\n"
            "        sigf = np.sqrt((nch+1.0)/np.maximum(nch, 1e-9))/SQN\n"
            "        occ = nch*np.maximum(1.0 + p['gn3']*sigf, 0.0)\n"
            "        fz_full = np.sqrt((1.0+occ)/(1.0+nch))")
src = patch(src, "            vp_b = vpar[idx].copy(); "
            "vq_b = vper[idx].copy()",
            "            vp_b = vpar[idx]*fz_full[idx]; "
            "vq_b = vper[idx]*fz_full[idx]")
src = patch(src, "    P(f\"  seed {seed}: BE-minus-simple best lnL = \"\n"
            "      f\"{best_lnl['BE']-best_lnl['simple']:+.1f}  \"\n"
            "      f\"({(time.time()-t0)/60:.1f} min)\")",
            "    P(f\"  seed {seed}: best lnL = {best_lnl['BE']:+.2f}  \"\n"
            "      f\"({(time.time()-t0)/60:.1f} min)\")")
src = patch(src, "'data/stage3u_summary.txt'",
            "'data/stage6p_summary.txt'", n=1)
src = patch(src, 'print("\\nbatch done; appended data/stage3u_summary.txt")',
            'print("batch done")')

RUNS = [('shaped', 4.0**0.5, 'N=4', 31),
        ('shaped', 2.0**0.5, 'N=2', 31),
        ('flat', 8.0**0.5, 'flat-8', 31),
        ('flat', 20.0**0.5, 'flat-20', 31),
        ('shaped', 1e12, 'off-s101', 101),
        ('shaped', 8.0**0.5, 'N=8-s101', 101)]
for mode, sqn, tag, seed in RUNS:
    with open('data/stage6p_summary.txt', 'a') as f:
        f.write(f"--- 6Pb run {tag} (mode={mode}, SQN={sqn:.4g}, "
                f"seed={seed}) ---\n")
    print(f"=== 6Pb {tag} ===", flush=True)
    sys.argv = ['stage6pb', '1p2', str(seed)]
    ns = {'SQN': sqn, 'MODE': mode, '_G': {}, '__name__': '__main__'}
    t0 = time.time()
    exec(src, ns)
    print(f"=== {tag} done ({(time.time()-t0)/60:.1f} min) ===", flush=True)

# ---------------- verdict assembly (guarded) ----------------
L = ["STAGE 6P-b: edge extension + s-flat control + seed 101 -- "
     "paired on the 6P baseline", ""]
try:
    txt = open('data/stage6p_summary.txt').read()
    res = {}
    for hdr, b in re.findall(r'--- 6Pb? run ([^\s]+) [^\n]*---\n'
                             r'(.*?)(?=--- 6Pb? run |\Z)', txt, re.S):
        m = re.search(r'seed (\d+) BE: a_hat=([0-9.]+) \(grid [0-9.]+, '
                      r'interior=(\w+)\)', b)
        mb = re.search(r'seed \d+: best lnL = ([+-][0-9.]+)', b)
        if m and mb:
            res[hdr] = (float(m.group(2)), m.group(3) == 'True',
                        float(mb.group(1)))
    base31 = res['off'][2]
    L.append("   run       a_hat  interior   best lnL    Delta vs off")
    for tag in ('off', 'N=60', 'N=20', 'N=8', 'N=4', 'N=2',
                'flat-8', 'flat-20'):
        if tag not in res: continue
        a, itr, bl = res[tag]
        L.append(f"  {tag:>8}  {a:5.3f}  {str(itr):>5}  {bl:+10.2f}"
                 f"   {bl-base31:+7.2f}")
    b101 = res.get('off-s101', (0, 0, np.nan))[2]
    if 'N=8-s101' in res:
        L.append(f"  seed 101: off {b101:+.2f}; shaped N=8 Delta = "
                 f"{res['N=8-s101'][2]-b101:+.2f}")
    L.append("")
    d8 = res['N=8'][2] - base31
    df8 = res['flat-8'][2] - base31
    spec = d8 - df8
    d101 = res['N=8-s101'][2] - b101 if 'N=8-s101' in res else np.nan
    ds = {t: res[t][2] - base31 for t in ('N=60', 'N=20', 'N=8', 'N=4',
                                          'N=2') if t in res}
    vals = [ds.get(t, np.nan) for t in ('N=8', 'N=4', 'N=2')]
    interior_turn = (not np.isnan(vals[1]) and not np.isnan(vals[2])
                     and vals[2] < vals[1])
    L.append(f"shaped(N=8) - flat(matched-8) = {spec:+.2f}; seed-101 "
             f"shaped Delta = {d101:+.2f}; turnover on extended grid: "
             f"{interior_turn} (profile {['%.1f' % v for v in vals]})")
    if spec >= 5.0 and (np.isnan(d101) or d101 > 0):
        call = ("SHAPE-SPECIFIC: the s-shape itself is preferred -- "
                "frozen-signature candidate"
                + ("" if interior_turn else
                   " (STILL EDGE-RIDING: N-hat unresolved, quote as "
                   "amplitude bound)"))
    elif abs(spec) < 3.0:
        call = ("GENERIC: s-blind broadening absorbs equally -- the "
                "old v7 residual, NOT frozen evidence; error-model "
                "finding only")
    else:
        call = "AMBIG: between bands or seed inconsistency; park"
    L.append(f"PRE-REGISTERED VERDICT: {call}")
except Exception as ex:
    L.append(f"verdict INCOMPLETE: {type(ex).__name__}: {ex}")

out = "\n".join(L)
print(out)
with open('data/stage6pb_extend.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage6pb_extend.txt")
