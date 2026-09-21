"""ROUND-47 VERIFICATION ADDENDUM (stage 10U pin-the-nu-form).

GA = the BLIND half: written and committed BEFORE the round-47 report
is read (the 87a4676 protocol, 14th execution). GB = the post-report
half (re-compute every load-bearing reviewer number); implemented
after the report lands.

Inheritance: the stage source calcs/stage10u_nuform.py is exec'd
TRUNCATED at its G10U-2 marker -- everything up to there (data
loading, world builder, FAMS, fit machinery, contest_fit) runs
bit-identically with zero copy divergence; nothing downstream (gates,
sky) executes. OUTFILE in the exec namespace is rebound to a scratch
path so no stage output can be touched.

GA legs:
  GA-1 independent objective evaluation -- the four member optima are
       refit with the stage engine (deterministic), then the profiled
       -2lnL is re-evaluated at those optima by a FRESH implementation
       (own profiling loop, own model arithmetic); the three sky
       deltas must match data/stage10u_skyread.txt within 0.5.
  GA-2 mpmath 30-digit member identity at the y grid (independent of
       both the stage vectorized forms and the stage's brentq checks);
       bar 1e-9.
  GA-3 pure-arithmetic recomputation of the 8-seed separation matrix
       and the A1 pooled separations FROM THE PRINTED STATS (regex
       parse of the outputs); bar +-0.15 (printed 0.1 rounding
       propagates ~0.1 worst-case).
  GA-4 rule-logic re-application: every demotion decision and the
       letter, re-derived from the printed powered set, bars, P_boot;
       must match exactly.
  GA-5 paired-bootstrap spot reproduction: reps 1-10 of the seed-202
       stream re-run (TH0B rebuilt by the stage's own selection rule);
       win indicators must be binomially consistent with the printed
       fractions (95% windows); rep fingerprints printed for GB/the
       reviewer.
  GA-6 mpmath H0 conversions of the four a0 + survivor-band
       arithmetic; bar 1e-6 relative.
"""
import math, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
MODE = sys.argv[1] if len(sys.argv) > 1 else 'ga'

OUT = []
def Q(s=""):
    print(s, flush=True)
    OUT.append(s)

def flush(path='data/round47_addendum.txt'):
    with open(path, 'a', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

# ---------- truncated-exec inheritance of the stage source ----------
SRC = open(os.path.join(HERE, 'stage10u_nuform.py'), encoding='utf-8').read()
MARK = "# ================= G10U-2: member identity"
cut = SRC.index(MARK)
NS = {'__name__': 'stage10u_trunc', '__file__':
      os.path.join(HERE, 'stage10u_nuform.py')}
_argv = sys.argv
sys.argv = ['stage10u_nuform.py', 'gates']
exec(compile(SRC[:cut], 'stage10u_trunc', 'exec'), NS)
sys.argv = _argv
NS['OUTFILE'] = os.path.join(ROOT, 'data', 'round47_scratch.txt')
W78, LEGA, FAMS = NS['W78'], NS['LEGA'], NS['FAMS']
build_sub, fit_hier, fit_deep = (NS['build_sub'], NS['fit_hier'],
                                 NS['fit_deep'])
nu_be = NS['nu_be']
MEMBERS = NS['MEMBERS']
ARC = NS['ARC']
S_ML, U_PRIOR, A0_FID = NS['S_ML'], NS['U_PRIOR'], NS['A0_FID']
LN10 = NS['LN10']
h0_of_a0 = NS['h0_of_a0']

GATES_TXT = open('data/stage10u_gates.txt', encoding='utf-8').read()
SKY_TXT = open('data/stage10u_skyread.txt', encoding='utf-8').read()

Q("")
Q("=" * 72)
Q(f"ROUND-47 ADDENDUM -- mode = {MODE}")
Q("=" * 72)

if MODE == 'ga':
    ok_all = True

    # ---------------- GA-1: independent objective evaluation --------
    Q("")
    Q("-- GA-1 independent objective at the refit optima --")
    subP = build_sub(W78, LEGA)
    b_plain = fit_hier(subP, nu_be, use_u=False)
    b_loose = fit_hier(subP, nu_be, use_u=True,
                       th0=list(b_plain.x) + [0.0])
    bBE = fit_deep(subP, nu_be, th0=list(b_loose.x))
    opt = {'BE': bBE}
    for nm in ('p065', 'gm', 'boot'):
        opt[nm] = fit_deep(subP, FAMS[nm], th0=list(bBE.x))

    from scipy.optimize import minimize_scalar as msc

    def indep_obj(nm, th):
        """Fresh implementation: profiled -2lnL at fixed (la0,f,s,u),
        own alternation over per-galaxy dml/dv to 1e-10."""
        la0, f, s_int, u = th
        a0 = 10**la0
        nu = FAMS[nm]
        tot_pts = 0.0
        tot_pri = 0.0
        for k, g in enumerate(LEGA):
            pts = W78['gpts'][int(g)]
            lg = W78['lgobs'][pts]
            gg = W78['gg'][pts]; gd = W78['gd'][pts]; gb = W78['gb'][pts]
            s2 = W78['sig2'][pts]
            sv = W78['sigv'][int(g)]
            isu = 1.0 if W78['fd'][int(g)] == 4 else 0.0
            se2 = s2 + s_int*s_int
            dml, dv = 0.0, 0.0
            prev = None
            for _ in range(300):
                gN = gg + f*math.exp(dml)*gd + gb
                r0 = lg - np.log10(gN*nu(gN/a0)) - u*isu
                w = 1.0/se2
                dv = float(np.sum(w*(r0 - 0.0))/(np.sum(w) + 1.0/sv**2))

                def od(dl):
                    gN2 = gg + f*math.exp(dl)*gd + gb
                    rr = lg - np.log10(gN2*nu(gN2/a0)) - dv - u*isu
                    return float(np.sum(rr*rr/se2) + dl*dl/(S_ML*S_ML))
                dml = msc(od, bounds=(-0.7, 0.7), method='bounded').x
                gN = gg + f*math.exp(dml)*gd + gb
                r = lg - np.log10(gN*nu(gN/a0)) - dv - u*isu
                cur = (float(np.sum(r*r/se2 + np.log(se2)))
                       + dml*dml/(S_ML*S_ML) + dv*dv/(sv*sv))
                if prev is not None and abs(prev - cur) < 1e-10:
                    prev = cur; break
                prev = cur
            tot_pts += prev
        return tot_pts + (u/U_PRIOR)**2

    iv = {nm: indep_obj(nm, list(opt[nm].x)) for nm in MEMBERS}
    ref = {'gm': 24.0, 'p065': 25.9, 'boot': 63.3}
    ga1_ok = True
    for nm in ('gm', 'p065', 'boot'):
        d_ind = iv[nm] - iv['BE']
        d_eng = opt[nm].fun - opt['BE'].fun
        ok = abs(d_ind - ref[nm]) <= 0.5
        ga1_ok &= ok
        Q(f"  d({nm} - BE): independent {d_ind:+.2f} | engine "
          f"{d_eng:+.2f} | sky print {ref[nm]:+.1f} -> "
          f"{'OK' if ok else 'MISMATCH'}")
    Q(f"GA-1: {'OK' if ga1_ok else 'MISMATCH'}")
    ok_all &= ga1_ok

    # ---------------- GA-2: mpmath 30-digit identity -----------------
    Q("")
    Q("-- GA-2 mpmath 30-digit member identity --")
    import mpmath as mp
    mp.mp.dps = 30

    def nu_mp(nm, y):
        y = mp.mpf(y)
        if nm == 'BE':
            return 1/(1 - mp.e**(-mp.sqrt(y)))
        if nm == 'p065':
            return (1 - mp.e**(-y**mp.mpf('0.65')))**(-1/mp.mpf('1.3'))
        if nm == 'gm':
            f = lambda v: v - 1 - 1/mp.expm1(y**mp.mpf('0.75')*mp.sqrt(v))
            return mp.findroot(f, mp.mpf(2.0))
        if nm == 'boot':
            f = lambda u: u - y - y/mp.expm1(u)
            u = mp.findroot(f, y + mp.mpf(1.0))
            return u/y

    ga2_ok = True
    for nm in MEMBERS:
        ds = []
        for y in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0):
            v = float(np.atleast_1d(FAMS[nm](np.array([y])))[0])
            vm = float(nu_mp(nm, y))
            ds.append(abs(v - vm)/vm)
        ok = max(ds) <= 1e-9
        ga2_ok &= ok
        Q(f"  {nm:5s}: max rel dev vs mpmath = {max(ds):.2e} "
          f"(bar 1e-9) -> {'OK' if ok else 'MISMATCH'}")
    Q(f"GA-2: {'OK' if ga2_ok else 'MISMATCH'}")
    ok_all &= ga2_ok

    # ---------------- GA-3: separation arithmetic from prints -------
    Q("")
    Q("-- GA-3 separation arithmetic re-derived from printed stats --")
    tr = {}
    for m in re.finditer(
            r"truth (\w+)\s*: d\(T beats R\) = ([^\n]+)", GATES_TXT):
        T = m.group(1)
        row = {}
        for mm in re.finditer(r"(\w+) ([+-][\d.]+)\+-([\d.]+)",
                              m.group(2)):
            row[mm.group(1)] = (float(mm.group(2)), float(mm.group(3)))
        tr[T] = row
    sep_print = {}
    lines = GATES_TXT.splitlines()
    for i, ln in enumerate(lines):
        mm = re.match(r"\s+(\w+)\s+vs\s+BE\s+p065\s+gm\s+boot", ln)
        if mm and i + 1 < len(lines):
            vals = lines[i+1].split()
            A = mm.group(1)
            for B, v in zip(MEMBERS, vals):
                if v != '.':
                    sep_print[(A, B)] = float(v)
    ga3_ok = True
    for A in MEMBERS:
        for B in MEMBERS:
            if A == B or (A, B) not in sep_print: continue
            mA, sA = tr[A][B]           # d_AB under truth A (positive)
            mBn, sB = tr[B][A]          # d_BA under truth B (positive)
            sep = (mA + mBn)/math.sqrt((sA*sA + sB*sB)/2.0)
            ok = abs(sep - sep_print[(A, B)]) <= 0.15
            ga3_ok &= ok
            if not ok:
                Q(f"  sep({A},{B}) recomputed {sep:.2f} vs printed "
                  f"{sep_print[(A, B)]:.2f} -> MISMATCH")
    for mm in re.finditer(
            r"(\w+)-(\w+): sep8 = ([\d.]+) -> sep24 = ([\d.]+) \(dA "
            r"([+-][\d.]+)\+-([\d.]+), dB ([+-][\d.]+)\+-([\d.]+)\)",
            GATES_TXT):
        A, B = mm.group(1), mm.group(2)
        dAm, dAs = float(mm.group(5)), float(mm.group(6))
        dBm, dBs = float(mm.group(7)), float(mm.group(8))
        sep = (dAm - dBm)/math.sqrt((dAs*dAs + dBs*dBs)/2.0)
        ok = abs(sep - float(mm.group(4))) <= 0.15
        ga3_ok &= ok
        Q(f"  A1 {A}-{B}: sep24 recomputed {sep:.2f} vs printed "
          f"{mm.group(4)} -> {'OK' if ok else 'MISMATCH'}")
    Q(f"GA-3: {'OK' if ga3_ok else 'MISMATCH'} (8-seed matrix + A1 "
      f"pooled, from printed stats, bar +-0.15)")
    ok_all &= ga3_ok

    # ---------------- GA-4: rule-logic re-application ---------------
    Q("")
    Q("-- GA-4 demotion logic re-applied from the printed record --")
    powered = set()
    mm = re.search(r"powered ordered pairs \(final, A1-pooled\): "
                   r"\[([^\]]*)\]", SKY_TXT)
    for pair in re.finditer(r"\('(\w+)', '(\w+)'\)", mm.group(1)):
        powered.add((pair.group(1), pair.group(2)))
    dem_lines = re.findall(
        r"(\w+) -> (\w+): powered \(sep ([\d.]+)\); sky d = "
        r"([+-][\d.]+) vs bar ([+-][\d.]+) \(truth-(\w+): "
        r"([+-][\d.]+)\+-([\d.]+)\); P_boot = ([\d.]+) -> (\w+)",
        SKY_TXT)
    ga4_ok = len(dem_lines) == len(powered)
    demoted = {}
    for (A, B, sep, dsky, bar, TB, mB, sB, pb, verdict) in dem_lines:
        bar_re = float(mB) + 2*float(sB) + 0.0
        hit = ((A, B) in powered and float(dsky) > bar_re
               and float(pb) >= 0.95)
        ok = (abs(bar_re - float(bar)) <= 0.15 and
              (verdict == 'DEMOTED') == hit)
        ga4_ok &= ok
        if hit: demoted.setdefault(B, []).append(A)
        if not ok:
            Q(f"  {A}->{B}: recomputed bar {bar_re:+.1f}/hit {hit} vs "
              f"printed {bar}/{verdict} -> MISMATCH")
    survivors = [m for m in MEMBERS if m not in demoted]
    letter = ('U-FORM-PINNED' if demoted and len(survivors) == 1 else
              'U-FORM-NARROWED' if demoted else
              'U-FORM-LEAN' if powered else 'U-POWER-LIMITED')
    lp = re.search(r"verdict letter: (\S+)", SKY_TXT).group(1)
    ga4_ok &= (letter == lp) and (sorted(demoted) == ['boot']) \
        and (survivors == ['BE', 'p065', 'gm'])
    Q(f"  re-derived: demoted = {dict(demoted)}, survivors = "
      f"{survivors}, letter = {letter} (printed {lp})")
    Q(f"GA-4: {'OK' if ga4_ok else 'MISMATCH'}")
    ok_all &= ga4_ok

    # ---------------- GA-5: paired-boot spot reproduction -----------
    Q("")
    Q("-- GA-5 paired bootstrap reps 1-10 (seed-202 stream) --")
    sky_fit = {}
    b_cold_be = fit_deep(subP, nu_be, th0=None)
    sky_fit['BE'] = b_cold_be if b_cold_be.fun < bBE.fun else bBE
    for nm in ('p065', 'gm', 'boot'):
        b_cold = fit_deep(subP, FAMS[nm], th0=None)
        sky_fit[nm] = (b_cold if b_cold.fun < opt[nm].fun
                       else opt[nm])
    TH0B = {m: list(sky_fit[m].x) for m in MEMBERS}
    rng5 = np.random.default_rng(202)
    rel_sig, uma_flag = {}, {}
    for g in W78['legA']:
        if W78['fd'][g] == 4:
            rel_sig[g] = 2.3/18.0; uma_flag[g] = True
        else:
            rel_sig[g] = W78['rel'][g]; uma_flag[g] = False
    wins = {('BE', 'gm'): 0, ('BE', 'boot'): 0, ('p065', 'boot'): 0,
            ('BE', 'p065'): 0}
    fps = []
    for rep in range(10):
        pick = rng5.choice(W78['legA'], size=len(W78['legA']),
                           replace=True)
        sc_shared = 1.0 + rng5.normal(0, 0.9/18.0)
        docc = []
        for g in pick:
            s_g = 1.0 + rng5.normal(0, rel_sig[g])
            s_g = min(max(s_g, 0.5), 1.5)
            if uma_flag[g]: s_g *= sc_shared
            docc.append(-math.log10(s_g))
        fps.append((int(np.sum(pick)), round(float(sc_shared), 6)))
        subr = build_sub(W78, list(pick), dlg_per_occ=docc)
        fr = {m: fit_deep(subr, FAMS[m], th0=TH0B[m]).fun
              for m in MEMBERS}
        for (A, B) in wins:
            if fr[B] - fr[A] > 0: wins[(A, B)] += 1
    Q(f"  rep fingerprints (sum pick ids, shared draw): {fps[:3]} ...")
    pref = {('BE', 'gm'): 0.915, ('BE', 'boot'): 0.990,
            ('p065', 'boot'): 1.000, ('BE', 'p065'): 0.955}
    from math import comb
    ga5_ok = True
    for (A, B), w in wins.items():
        p = pref[(A, B)]
        # exact binomial two-sided 95% acceptance for n=10
        probs = [comb(10, k)*p**k*(1-p)**(10-k) for k in range(11)]
        lo = 0
        while sum(probs[:lo+1]) < 0.025: lo += 1
        hi = 10
        while sum(probs[hi:]) < 0.025: hi -= 1
        ok = lo <= w <= hi
        ga5_ok &= ok
        Q(f"  P({A} beats {B}) printed {p:.3f}: reps 1-10 wins = "
          f"{w}/10 (95% window [{lo},{hi}]) -> "
          f"{'OK' if ok else 'MISMATCH'}")
    Q(f"GA-5: {'OK' if ga5_ok else 'MISMATCH'}")
    ok_all &= ga5_ok

    # ---------------- GA-6: mpmath conversions ----------------------
    Q("")
    Q("-- GA-6 mpmath H0 conversions + band arithmetic --")
    ga6_ok = True
    href = {'BE': 65.4, 'p065': 70.4, 'gm': 66.3, 'boot': 67.6}
    for nm in MEMBERS:
        a0 = 10**sky_fit[nm].x[0]
        hm = float(2*mp.pi*mp.mpf(a0)/mp.mpf(299792458.0) /
                   mp.mpf('3.240779e-20'))
        ok = (abs(hm - h0_of_a0(a0))/hm <= 1e-6 and
              abs(hm - href[nm]) <= 0.06)
        ga6_ok &= ok
        Q(f"  {nm:5s}: H0(mpmath) = {hm:.2f} vs engine "
          f"{h0_of_a0(a0):.2f} vs printed {href[nm]:.1f} -> "
          f"{'OK' if ok else 'MISMATCH'}")
    surv_band = (min(href[m] for m in ('BE', 'p065', 'gm')),
                 max(href[m] for m in ('BE', 'p065', 'gm')))
    ok = surv_band == (65.4, 70.4)
    ga6_ok &= ok
    Q(f"  survivor band = {surv_band} vs printed 65.4-70.4 -> "
      f"{'OK' if ok else 'MISMATCH'}")
    Q(f"GA-6: {'OK' if ga6_ok else 'MISMATCH'}")
    ok_all &= ga6_ok

    Q("")
    Q(f"GA BLIND HALF: {'ALL OK' if ok_all else 'MISMATCHES PRESENT'}")
    flush()

elif MODE == 'gb':
    Q("GB half: implemented after the round-47 report is read.")
    flush()
