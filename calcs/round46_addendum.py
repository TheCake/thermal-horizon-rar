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
    P("")
    P("=" * 70)
    P("ROUND 46 ADDENDUM -- GB (post-report reviewer-number checks)")
    P("=" * 70)
    P("(appended after the round-46 report; see below)")
    save()
print("\nsaved: data/round46_addendum.txt")
