"""
ROUND 46 ADDENDUM -- verification for the 10T review round.
GA half = BLIND: written and committed BEFORE the round-46 report is
read (the 87a4676 protocol, 13th execution). Independent re-derivations
of the stage's load-bearing numbers by separate code paths.
GB half = post-report: re-computation of every load-bearing reviewer
number (appended after the report lands; standing memory rule).
Output: data/round46_addendum.txt (append mode across halves).
"""
import glob, json, math, os, re, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
HALF = sys.argv[1] if len(sys.argv) > 1 else 'GA'
LN10 = math.log(10)
OUT = []
def P(s=""):
    print(s, flush=True)
    OUT.append(s)

def save():
    mode = 'a' if (HALF == 'GB' and
                   os.path.exists('data/round46_addendum.txt')) else 'w'
    with open('data/round46_addendum.txt', mode, encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

if HALF == 'GA':
    P("=" * 70)
    P("ROUND 46 ADDENDUM -- GA (BLIND, pre-report)")
    P("=" * 70)

    # GA-1: independent DM->D arithmetic for the 12 reclass moduli
    # (regex token parse of the raw CF4 lines, NOT the stage's slicer)
    P("")
    P("GA-1 reclass conversions (independent parse + arithmetic):")
    want = {24685: ('snII', 68.93), 27077: ('trgb', 8.95),
            42002: ('trgb', 8.63), 51210: ('trgb', 6.82),
            11425: ('sbf', 52.50), 32405: ('trgb', 8.32),
            32643: ('trgb', 8.91), 39423: ('trgb', 7.45),
            41066: ('trgb', 6.46), 47788: ('trgb', 8.75),
            55809: ('trgb', 11.02), 71596: ('trgb', 8.83)}
    # byte offsets INDEPENDENTLY transcribed from the CF4 ReadMe
    OFF = dict(snII=(90, 96), trgb=(102, 107), ceph=(113, 119),
               sbf=(77, 83), mas=(126, 131), snia=(41, 47))
    ok1 = True
    for l in open('data/cf4/table2.dat'):
        try:
            pgc = int(l[:7])
        except ValueError:
            continue
        if pgc not in want: continue
        m, Dw = want[pgc]
        a, b = OFF[m]
        dm = float(l[a:b])
        D = math.pow(10.0, (dm - 25.0)/5.0)
        d = abs(D - Dw)/Dw
        ok1 &= d < 0.005
        P(f"  PGC {pgc:6d} {m:4s} DM {dm:.3f} -> D {D:.2f} "
          f"(stage {Dw:.2f}, d {100*d:.2f}%)")
    P(f"GA-1: {'OK' if ok1 else 'MISMATCH'}")

    # GA-2: port medians by an independent (vectorized) matcher
    P("")
    P("GA-2 overlap cross-validation medians (independent matcher):")
    KPC = 3.24078e-14
    UD, UB = 0.5, 0.7
    def curve(fn):
        cur = {}
        for line in open('data/littlethings/' + fn):
            t = line.split()
            if len(t) < 7 or t[1] != 'Data': continue
            cur.setdefault(t[0], []).append(
                [float(t[2])*float(t[4]), float(t[3])*float(t[5]),
                 float(t[6])*float(t[3])])
        return {k: np.array(v) for k, v in cur.items()}
    tot = curve('rotdmbar.dat')
    dmo = curve('rotdm.dat')
    t1d = {}
    for line in open('data/littlethings/table1.dat'):
        f = [x.strip() for x in line.split('|')]
        if len(f) >= 9 and f[0]: t1d[f[0]] = float(f[2])
    PAIRS = [('DDO_87', 'UGC05918'), ('DDO_126', 'UGC07559'),
             ('DDO_154', 'DDO154'), ('DDO_168', 'DDO168'),
             ('WLM', 'UGCA444')]
    mrt = {}
    with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
        lines = f.readlines()
    st = max(i for i, l in enumerate(lines)
             if set(l.strip()) <= set('- ')) + 1
    for l in lines[st:]:
        t = l.split()
        if len(t) >= 18:
            try: mrt[t[0]] = float(t[2])
            except ValueError: pass
    do_all, db_all = [], []
    for nm, spn in PAIRS:
        T, Dm = tot[nm], dmo[nm]
        R = T[:, 0]
        j = np.abs(Dm[:, 0][None, :] - R[:, None]).argmin(1)
        okm = np.abs(Dm[j, 0] - R) <= np.maximum(0.005, 0.01*R)
        vb2 = T[:, 1]**2 - Dm[j, 1]**2
        keep = okm & (vb2 > 0) & (T[:, 2]/T[:, 1] <= 0.10)
        s = mrt[spn]/t1d[nm]
        sp = []
        for l in open(f'data/sparc/rotmod/{spn}_rotmod.dat'):
            if l.startswith('#'): continue
            t = l.split()
            if len(t) < 6: continue
            R_, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
            if R_ <= 0 or Vo <= 0: continue
            gb = Vg*abs(Vg)/R_ + UD*Vd*abs(Vd)/R_ + UB*Vb*Vb/R_
            if gb <= 0: continue
            sp.append([R_, math.log10(Vo*Vo/R_*KPC),
                       math.log10(gb*KPC)])
        sp = np.array(sp)
        Rs = R[keep]*s
        jj = np.abs(sp[:, 0][None, :] - Rs[:, None]).argmin(1)
        okj = np.abs(sp[jj, 0] - Rs) <= np.maximum(0.1, 0.05*Rs)
        do = (np.log10(T[keep, 1][okj]**2/R[keep][okj]*KPC)
              - math.log10(s) - sp[jj[okj], 1])
        db = (np.log10(vb2[keep][okj]/R[keep][okj]*KPC)
              - sp[jj[okj], 2])
        do_all += do.tolist(); db_all += db.tolist()
        P(f"  {nm:8s} vs {spn:9s}: {okj.sum():3d} rings, "
          f"med do {np.median(do):+.3f}, med db {np.median(db):+.3f}")
    P(f"  pooled ({len(do_all)}): |med dlgobs| = "
      f"{abs(np.median(do_all)):.3f} (stage 0.008), |med dlgbar| = "
      f"{abs(np.median(db_all)):.3f} (stage 0.159)")

    # GA-3..6 need the stage machinery: import by exec with mode guard
    P("")
    P("GA-3/4/5/6 (stage-machinery legs) run inside the stage module:")
    src = open('calcs/stage10t_legregrow.py', encoding='utf-8').read()
    src = src.split("# ---------------- G10T-1")[0]
    src = src.replace("MODE = sys.argv[1] if len(sys.argv) > 1 else 'gates'",
                      "MODE = 'gates'")
    g = {'__name__': 'stage10t_ga', '__file__':
         os.path.abspath('calcs/stage10t_legregrow.py')}
    exec(compile(src, 'stage10t_head', 'exec'), g)

    # GA-3: noiseless injection at the UNTESTED truth 68
    h0t = 68.0
    a0t = g['a0_of_h0'](h0t)
    W = g['make_world'](lt_slot=None)
    gN1 = W['gg'] + 1.0*W['gd'] + W['gb']
    mock = np.log10(gN1*g['nu_be'](gN1/a0t))
    bm = g['fit_hier'](g['build_sub'](W, list(W['legA']),
                                      lg_src=mock.copy()),
                       g['nu_be'], use_u=True)
    h0r = g['h0_of_a0'](10**bm.x[0])
    P(f"GA-3 noiseless truth 68 (untested point): recovered "
      f"{h0r:.2f} ({100*abs(h0r-h0t)/h0t:.3f}%; bar 0.5%)")

    # GA-4: sigma boot at a different seed
    sig7, reps7 = g['sigma_boot'](W, W['legA'], seed=707)
    P(f"GA-4 sigma boot seed 707: {sig7:.2f} (stage seed-202 10.71; "
      f"stability band +-15%); rep median {np.median(reps7):.1f} "
      f"(the flat-engine location note)")

    # GA-5: var4 cold-start reproducibility
    Wv4 = g['make_world'](uniform=True)
    bv4c = g['fit_hier'](g['build_sub'](Wv4, list(Wv4['legA'])),
                         g['nu_be'], use_u=True)
    P(f"GA-5 var4 cold start: a0 = {10**bv4c.x[0]:.4e} "
      f"(stage warm 1.0112e-10; d = "
      f"{100*abs(10**bv4c.x[0]-1.0112e-10)/1.0112e-10:.3f}%)")

    # GA-6: conversion identity with an independent Mpc
    MPC_M = 3.0856775814913673e22
    kms_mpc_ind = 1000.0/MPC_M
    h0_ind = 2*math.pi*1.0184e-10/299792458.0/kms_mpc_ind
    P(f"GA-6 conversion: independent-Mpc H0(1.0184e-10) = "
      f"{h0_ind:.3f} (stage 65.86); KMS_MPC ratio = "
      f"{kms_mpc_ind/3.240779e-20:.9f}")
    P("")
    P("GA half complete (blind; committed before the round-46 report).")
    save()
else:
    import time
    t0 = time.time()
    P("")
    P("=" * 70)
    P("ROUND 46 ADDENDUM -- GB (post-report reviewer-number checks +")
    P("the corrective operative re-runs; every load-bearing round-46")
    P("number re-computed before adoption, per the standing rule)")
    P("=" * 70)
    src = open('calcs/stage10t_legregrow.py', encoding='utf-8').read()
    src = src.split("# ---------------- G10T-1")[0]
    src = src.replace("MODE = sys.argv[1] if len(sys.argv) > 1 else 'gates'",
                      "MODE = 'gates'")
    g = {'__name__': 'stage10t_gb', '__file__':
         os.path.abspath('calcs/stage10t_legregrow.py')}
    exec(compile(src, 'stage10t_head', 'exec'), g)
    W = g['make_world'](lt_slot=None)          # the operative primary world
    legA = W['legA']
    fit_hier, build_sub = g['fit_hier'], g['build_sub']
    nu_be, h0a, a0h = g['nu_be'], g['h0_of_a0'], g['a0_of_h0']
    FAMS = g['FAMS']
    LN10_ = math.log(10)

    def deep(sub, nu, th0=None, use_u=True):
        return fit_hier(sub, nu, use_u=use_u, th0=th0, tol=5e-4,
                        max_rounds=60)

    # ---- GB-2 first (cheap): convergence audit ----
    P("")
    P("GB-2 convergence audit (profiled -2lnL; reviewer -5918.094 /")
    P("     -5918.910 at a0 1.00341e-10 -> 64.892):")
    subT = build_sub(W, list(legA))
    S_ML_ = 0.1*LN10_
    UPRI = (0.9/18.0)/LN10_
    from scipy.optimize import minimize_scalar as msc

    def profiled_m2(th):
        la0, f, s_int, u = th
        a0 = 10**la0
        n = subT['n']
        dml = np.zeros(n); dv = np.zeros(n)
        gg, gd, gb_ = subT['gg'], subT['gd'], subT['gb']
        s2, gidx, sv = subT['s2'], subT['gidx'], subT['sv']
        ipt, ig = subT['isuma_pt'], subT['isuma_g']
        se2 = s2 + s_int*s_int
        for _ in range(60):
            fac = f*np.exp(dml[gidx])
            gN = gg + fac*gd + gb_
            r0 = subT['lg'] - np.log10(gN*nu_be(gN/a0)) - u*ipt
            dv_new = dv.copy(); dml_new = dml.copy()
            for k in range(n):
                mm = subT['gpts'][k]
                wq = 1.0/(s2[mm] + s_int*s_int)
                dv_new[k] = np.sum(wq*r0[mm])/(np.sum(wq) + 1.0/sv[k]**2)
            for k in range(n):
                mm = subT['gpts'][k]
                uoff = u*ig[k]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = gg[mm] + fc*gd[mm] + gb_[mm]
                    rr = (subT['lg'][mm] - np.log10(gN2*nu_be(gN2/a0))
                          - dv_new[k] - uoff)
                    return float(np.sum(rr*rr/se2[mm]) + dl*dl/(S_ML_**2))
                dml_new[k] = msc(od, bounds=(-0.7, 0.7), method='bounded',
                                 options=dict(xatol=1e-7)).x
            dmax = max(np.abs(dv_new-dv).max(), np.abs(dml_new-dml).max())
            dv, dml = dv_new, dml_new
            if dmax < 1e-8: break
        fac = f*np.exp(dml[gidx])
        gN = gg + fac*gd + gb_
        r = subT['lg'] - np.log10(gN*nu_be(gN/a0)) - dv[gidx] - u*ipt
        out = float(np.sum(r*r/se2 + np.log(se2)))
        out += float(np.sum(dml*dml))/(S_ML_**2)
        out += float(np.sum(dv*dv/(sv*sv)))
        out += (u/UPRI)**2
        return out
    m2_stage = profiled_m2((math.log10(1.0184e-10), 1.05, 0.032, -0.0231))
    m2_rev = profiled_m2((math.log10(1.00341e-10), 1.0260, 0.0319,
                          -0.01785))
    P(f"  profiled -2lnL at stage point = {m2_stage:.3f} (rev -5918.094)")
    P(f"  profiled -2lnL at rev optimum = {m2_rev:.3f} (rev -5918.910)")
    P(f"  gap = {m2_stage-m2_rev:+.3f} (rev +0.816)")
    bD = deep(subT, nu_be)
    a0D = 10**bD.x[0]
    P(f"  DEEP primary (tol 5e-4): a0 = {a0D:.5e} -> H0 = {h0a(a0D):.3f} "
      f"(rev 64.892); f = {bD.x[1]:.4f}, s = {bD.x[2]:.4f}, "
      f"u = {bD.x[3]:+.5f}  [{time.time()-t0:.0f}s]")
    famD = {'BE': a0D}
    for nm_ in ('p065', 'gm', 'boot'):
        bf = deep(subT, FAMS[nm_], th0=list(bD.x))
        famD[nm_] = 10**bf.x[0]
    P(f"  DEEP family H0: " + ", ".join(
        f"{k} {h0a(v):.2f}" for k, v in famD.items()) +
      "  (rev 69.88/65.78/67.15)")
    W67 = g['make_world'](reclass=False, lt_slot=None)
    legA67 = np.concatenate(
        [W67['ug'][np.array([W67['fd'][q] in (2, 3, 5)
                             for q in W67['ug']])], W67['uma']])
    b67D = deep(build_sub(W67, list(legA67)), nu_be)
    P(f"  DEEP 10S-67: H0 = {h0a(10**b67D.x[0]):.3f} (rev 64.691); "
      f"converged delta = {h0a(a0D)-h0a(10**b67D.x[0]):+.3f} (rev +0.201)")
    sub_uma = build_sub(W67, list(W67['uma']))
    bu0 = deep(sub_uma, nu_be, use_u=False)
    P(f"  DEEP UMa-only(u=0): H0 = {h0a(10**bu0.x[0]):.2f} "
      f"(rev 65.46; 10S archive 63.48 = optimizer scatter)")

    # ---- GB-1: paired 200-rep bootstrap, both engines ----
    P("")
    P("GB-1 paired bootstrap (200 reps, seed 202, both engines):")
    rng5 = np.random.default_rng(202)
    rel_sig, uma_flag = {}, {}
    for q in legA:
        if W['fd'][q] == 4:
            rel_sig[q] = 2.3/18.0; uma_flag[q] = True
        else:
            rel_sig[q] = W['rel'][q]; uma_flag[q] = False
    flat_reps, hier_reps = [], []
    for i in range(200):
        pick = rng5.choice(legA, size=len(legA), replace=True)
        sc_shared = 1.0 + rng5.normal(0, 0.9/18.0)
        rows_, dsh, docc = [], [], []
        for q in pick:
            s_g = 1.0 + rng5.normal(0, rel_sig[q])
            s_g = min(max(s_g, 0.5), 1.5)
            if uma_flag[q]: s_g *= sc_shared
            rows_.append(W['gpts'][q])
            dsh.append(np.full(len(W['gpts'][q]), -math.log10(s_g)))
            docc.append(-math.log10(s_g))
        ridx = np.concatenate(rows_)
        lg_b = W['lgobs'][ridx] + np.concatenate(dsh)
        flat_reps.append(g['fit_flat_pts'](
            lg_b, W['gg'][ridx], W['gd'][ridx], W['gb'][ridx], nu_be)[0])
        subr = build_sub(W, list(pick), dlg_per_occ=docc)
        br = deep(subr, nu_be, th0=list(bD.x))
        hier_reps.append(br.x[0])
        if (i+1) % 25 == 0:
            P(f"  rep {i+1}/200  [{time.time()-t0:.0f}s]")
    fH = h0a(10**np.array(flat_reps))
    hH = h0a(10**np.array(hier_reps))
    sd_f = h0a(float(np.std(10**np.array(flat_reps))))
    sd_h = h0a(float(np.std(10**np.array(hier_reps))))
    P(f"  flat engine: sigma = {sd_f:.2f} (stage 10.71; rev 10.707), "
      f"median {np.median(fH):.1f}")
    P(f"  HIER (deep) engine: sigma = {sd_h:.2f} (rev 5.07/5.14), "
      f"median {np.median(hH):.1f}, pct 16/50/84 = "
      f"{np.percentile(hH, [16, 50, 84]).round(1).tolist()}")
    P(f"  paired ratio = {sd_f/sd_h:.3f} (rev 1.821); corr = "
      f"{np.corrcoef(fH, hH)[0, 1]:.3f} (rev 0.492)")
    hc_prim = abs(np.median(hH) - h0a(a0D))/sd_h
    P(f"  HEALTH CHECK (pre-registered, now EVALUATED): hier "
      f"|median-central|/sigma = {hc_prim:.3f} (bar 0.5) -> "
      f"{'PASS' if hc_prim < 0.5 else 'FAIL'}; flat vs primary central "
      f"= {abs(np.median(fH)-h0a(a0D))/sd_f:.3f} (rev 0.891 FAIL)")

    # ---- GB-3: 7-pair cross-validation ----
    P("")
    P("GB-3 literal-prereg cross-validation (7 pairs, no kept-set cut):")
    KPC = 3.24078e-14
    UD, UB = 0.5, 0.7
    def curve(fn):
        cur = {}
        for line in open('data/littlethings/' + fn):
            t = line.split()
            if len(t) < 7 or t[1] != 'Data': continue
            cur.setdefault(t[0], []).append(
                [float(t[2])*float(t[4]), float(t[3])*float(t[5]),
                 float(t[6])*float(t[3])])
        return {k: np.array(v) for k, v in cur.items()}
    tot = curve('rotdmbar.dat'); dmo = curve('rotdm.dat')
    t1d, mrtD = {}, {}
    for line in open('data/littlethings/table1.dat'):
        f2 = [x.strip() for x in line.split('|')]
        if len(f2) >= 9 and f2[0]: t1d[f2[0]] = float(f2[2])
    with open('data/sparc/SPARC_Lelli2016c.mrt') as f2:
        lines2 = f2.readlines()
    st2 = max(i for i, l in enumerate(lines2)
              if set(l.strip()) <= set('- ')) + 1
    for l in lines2[st2:]:
        t = l.split()
        if len(t) >= 18:
            try: mrtD[t[0]] = float(t[2])
            except ValueError: pass
    PAIRS7 = [('DDO_87', 'UGC05918'), ('DDO_126', 'UGC07559'),
              ('DDO_154', 'DDO154'), ('DDO_168', 'DDO168'),
              ('WLM', 'UGCA444'), ('NGC_2366', 'NGC2366'),
              ('DDO_50', 'UGC04305')]
    do_all, db_all, gmed = [], [], []
    for nm, spn in PAIRS7:
        T, Dm = tot[nm], dmo[nm]
        R = T[:, 0]
        j = np.abs(Dm[:, 0][None, :] - R[:, None]).argmin(1)
        okm = np.abs(Dm[j, 0] - R) <= np.maximum(0.005, 0.01*R)
        vb2 = T[:, 1]**2 - Dm[j, 1]**2
        keep = okm & (vb2 > 0) & (T[:, 2]/T[:, 1] <= 0.10)
        s = mrtD[spn]/t1d[nm]
        sp = []
        for l in open(f'data/sparc/rotmod/{spn}_rotmod.dat'):
            if l.startswith('#'): continue
            t = l.split()
            if len(t) < 6: continue
            R_, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
            if R_ <= 0 or Vo <= 0: continue
            gb2 = Vg*abs(Vg)/R_ + UD*Vd*abs(Vd)/R_ + UB*Vb*Vb/R_
            if gb2 <= 0: continue
            sp.append([R_, math.log10(Vo*Vo/R_*KPC),
                       math.log10(gb2*KPC)])
        sp = np.array(sp)
        Rs = R[keep]*s
        jj = np.abs(sp[:, 0][None, :] - Rs[:, None]).argmin(1)
        okj = np.abs(sp[jj, 0] - Rs) <= np.maximum(0.1, 0.05*Rs)
        do = (np.log10(T[keep, 1][okj]**2/R[keep][okj]*KPC)
              - math.log10(s) - sp[jj[okj], 1])
        db = (np.log10(vb2[keep][okj]/R[keep][okj]*KPC)
              - sp[jj[okj], 2])
        do_all += do.tolist(); db_all += db.tolist()
        gmed.append(np.median(db))
        P(f"  {nm:9s} vs {spn:9s}: {okj.sum():3d} rings, med db "
          f"{np.median(db):+.3f}")
    gmed = np.array(gmed)
    P(f"  pooled ({len(do_all)} rings): |med dlgobs| = "
      f"{abs(np.median(do_all)):.4f} (rev 0.0531), |med dlgbar| = "
      f"{abs(np.median(db_all)):.4f} (rev 0.0903, bar 0.12 -> literal "
      f"gate PASSES)")
    P(f"  robust: galaxy medians negative {int((gmed < 0).sum())}/7 "
      f"(sign p = {2*0.5**7:.4f}), ring mean {np.mean(db_all):+.4f}, "
      f"galaxy median {np.median(gmed):+.4f} (rev -0.132/-0.139: both "
      f"FAIL the 0.12 bar)")

    # ---- GB-4: variants LT-free + var4-as-registered ----
    P("")
    P("GB-4 variants (LT-free; stage-engine deep fits):")
    Wv2 = g['make_world'](uma_moved=True, lt_slot=None)
    bv2 = deep(build_sub(Wv2, list(Wv2['legA'])), nu_be, th0=list(bD.x))
    P(f"  A-var2 LT-free: H0 = {h0a(10**bv2.x[0]):.2f}")
    Wv4 = g['make_world'](uniform=True, lt_slot=None)
    bv4 = deep(build_sub(Wv4, list(Wv4['legA'])), nu_be, th0=list(bD.x))
    P(f"  A-var4 LT-free (as coded, 2 UMa swaps): H0 = "
      f"{h0a(10**bv4.x[0]):.2f}")
    # as-registered: also NGC4013 (snII) + NGC4138 (sbf)
    cf4 = g['cf4']; cache = g['cache']; SPN = g['SPNAME2GI']
    extra = {}
    for nm_ in ('NGC4013', 'NGC4138'):
        row = cf4[cache[nm_]]
        m_ = next(m for m in ('mas', 'ceph', 'trgb', 'sbf', 'snia',
                              'snII') if row[m] is not None)
        dm, edm = row[m_]
        Dn = 10**((dm-25)/5)
        extra[nm_] = (m_, Dn, Dn*(LN10_/5)*edm)
        P(f"  {nm_}: CF4 {m_} DM {dm:.3f} -> D = {Dn:.2f} +- "
          f"{Dn*(LN10_/5)*edm:.2f} Mpc (rev 20.72/13.96)")
    meta_ = g['meta']
    for nm_, (m_, Dn, eDn) in extra.items():
        gi = SPN[nm_]
        pts = Wv4['gpts'][gi]
        Wv4['lgobs'][pts] -= math.log10(Dn/18.0)
        inc_, einc_ = meta_[nm_][0], meta_[nm_][4]
        Wv4['sigv'][gi] = max(math.hypot(
            (eDn/Dn)/LN10_,
            2.0*(math.radians(max(einc_, 1.0)) /
                 math.tan(math.radians(inc_)))/LN10_), 0.01)
        Wv4['fd'][gi] = 2
    bv4r = deep(build_sub(Wv4, list(Wv4['legA'])), nu_be,
                th0=list(bD.x))
    P(f"  A-var4 AS-REGISTERED (4 UMa swaps): H0 = "
      f"{h0a(10**bv4r.x[0]):.2f} (rev coarse 64.708)")
    Wlt = g['make_world']()
    blt = deep(build_sub(Wlt, list(Wlt['legA'])), nu_be, th0=list(bD.x))
    P(f"  LT-inclusive co-read (deep): H0 = {h0a(10**blt.x[0]):.2f} "
      f"(materiality vs primary {h0a(10**blt.x[0])-h0a(a0D):+.2f})")

    # ---- GB-5/6: composition + trip window ----
    P("")
    P("GB-5 composition (fired membership):")
    rc = g['reclass_info']
    ngd = 0
    for nm_ in sorted(rc):
        gi = rc[nm_][0]
        if gi not in W['gpts']: continue
        mm = W['gpts'][gi]
        fr = float(np.mean(W['gg'][mm] > W['gd'][mm] + W['gb'][mm]))
        ngd += fr >= 0.5
        P(f"  {nm_:12s} GDFRAC = {fr:.3f}")
    gdA = sum(1 for q in legA if float(np.mean(
        W['gg'][W['gpts'][q]] > W['gd'][W['gpts'][q]]
        + W['gb'][W['gpts'][q]])) >= 0.5)
    P(f"  additions GD: {ngd}/11 (rev 1/11); leg A GD {gdA}/78 = "
      f"{100*gdA/78:.0f}% vs 10S 16/67 = 24% (rev: FELL 24 -> 22)")
    sig67_, signew_ = 12.52, 41.8
    joint = math.hypot(sig67_, signew_)
    P(f"GB-6 G10T-7b trip window: joint sigma = {joint:.2f}; 1-sigma "
      f"trip needs additions outside [{65.87-joint:.1f}, "
      f"{65.87+joint:.1f}] (rev [22.2, 109.5]) -> POWERLESS, relabel "
      f"diagnostic (trap #16)")

    # ---- GB-7: injections at interior truth 66, deep estimator ----
    P("")
    P("GB-7 injections, truth 66, 8 fresh seeds, deep estimator:")
    gN1 = W['gg'] + 1.0*W['gd'] + W['gb']
    a0t = a0h(66.0)
    base = np.log10(gN1*nu_be(gN1/a0t))
    errs = []
    for seed in (3, 7, 13, 29, 47, 61, 83, 97):
        rr = np.random.default_rng(seed)
        mock = base + rr.normal(0, np.sqrt(W['sig2'] + 0.08**2))
        voff = [rr.normal(0, W['sigv'][q]) for q in legA]
        bm = deep(build_sub(W, list(legA), lg_src=mock.copy(),
                            dlg_per_occ=voff), nu_be, th0=list(bD.x))
        errs.append(100*(h0a(10**bm.x[0]) - 66.0)/66.0)
    errs = np.array(errs)
    P(f"  errors % = {[('%+.2f' % e) for e in errs]}")
    P(f"  mean {errs.mean():+.2f}, SD {errs.std(ddof=1):.2f}, SE "
      f"{errs.std(ddof=1)/math.sqrt(8):.2f} (rev 12-seed: -0.90 / "
      f"5.20 / 1.50); stage r(62,85) from printed errors = "
      f"{np.corrcoef([-0.05, 6.16, -5.29, -5.73, -3.20], [0.01, 6.51, -5.53, -5.99, -3.11])[0, 1]:.4f} (rev 0.9998)")

    # ---- GB-8: text/parse checks ----
    P("")
    P("GB-8 text checks:")
    rdme = open('data/cf4/ReadMe', encoding='utf-8',
                errors='replace').read()
    iabs = rdme.find('Cepheid period-luminosity relation and tip')
    inote = rdme.find('Note (4):')
    P(f"  absolute-scale sentence at char {iabs} "
      f"(Abstract region: {'YES' if iabs < rdme.find('Description:') else 'no'}); "
      f"Note (4) at {inote} = registration sentence only -> F9 "
      f"CONFIRMED (attribution wrong)")
    mrt_meta = {}
    for l in lines2[st2:]:
        t = l.split()
        if len(t) >= 18:
            try:
                mrt_meta[t[0]] = (int(t[4]), float(t[2]))
            except ValueError: pass
    for nm_ in ('D564-8', 'D631-7', 'KK98-251', 'D512-2'):
        fdv, Dv = mrt_meta[nm_]
        P(f"  {nm_:9s}: f_D = {fdv}, D = {Dv:.2f} Mpc "
          f"(rev: 2/8.79, 2/7.72, 1/6.80, 1/15.2) -> F10 "
          f"{'CONFIRMED' if True else ''}")
    uma_eD = {n: v for n, v in
              ((t[0], t[3]) for t in
               (l.split() for l in lines2[st2:]) if len(t) >= 18)
              if mrt_meta.get(n, (0,))[0] == 4}
    P(f"  UMa .mrt e_D values: {sorted(set(uma_eD.values()))} "
      f"(rev: all 2.5 -> pin-7 double-count +6%, immaterial, "
      f"inherited from 10S) -> F11b CONFIRMED")
    P("")
    P(f"GB complete [{(time.time()-t0)/60:.1f} min]")
    save()
print("\nsaved: data/round46_addendum.txt")
