"""STAGE 8F reader — THE FAT-TAIL ARMS (W4 manufacture test;
pre-reg in NOTES 2026-07-31, committed BEFORE any run; referee T6).

Reads the ARMTAG'd fullarmw cubes (operative photow3 fitter on
fat-tail-truth injections; the tail is the ONLY truth-model
mismatch — companions model-matched flat-q at the pinned 0.10).

G8F-WT: each config's OUT must carry the WIDTH-SHAPE injection
line matching the registered WTRUTH exactly (float-compared) plus
the injected-hist checksum line.  SYMPTOM-MATCH RULE + BARS as
locked in NOTES.  Incremental: reports whatever exists; the
verdict fires only when the full pre-registered set is on disk.
Output: data/stage8f_read.txt
"""
import os
import re
import numpy as np

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID = np.array([1.05, 1.3])
WR_GRID = np.array([0.10, 0.20, 0.30, 0.40, 0.50])
FCOMP = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FPM = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
KW = np.array([0.7, 1.0, 1.4])
SQ = np.array([0.0, 0.1, 0.2, 0.3])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
LAWS = ('simple', 'BE')

# the registry (pre-reg): armtag -> (WTRUTH fields, expected seeds)
# F0 sweep = seed 31; the calibrated member gains seed 101; fb/fc
# are filled in AFTER calibration (ftl* substituted) — the reader
# derives them from the calibration result so no number is chosen
# post-hoc: truths and bars are fixed, only ftl* is data-selected
# per the locked symptom-match rule.
F0 = {
    't05': (0.00, 0.10, 1.2, 0.0, 0.0, 0.05),
    't10': (0.00, 0.10, 1.2, 0.0, 0.0, 0.10),
    't20': (0.00, 0.10, 1.2, 0.0, 0.0, 0.20),
    't35': (0.00, 0.10, 1.2, 0.0, 0.0, 0.35),
    'tsq': (0.00, 0.10, 1.2, 0.2, 0.0, 0.10),
}
PURE = ('t05', 't10', 't20', 't35')
INJRE = re.compile(
    r"WIDTH-SHAPE injection (\w+) alpha=([\d.]+) at fcomp=([\d.]+), "
    r"fpm=([\d.]+), sq=([\d.]+), flr=([\d.]+), ftl=([\d.]+)")

P("8F THE FAT-TAIL ARMS — reader (pre-reg in NOTES 2026-07-31; "
  "bars locked before any run)")
P("")

def gate_wt(tag, alpha, fcm, fpm, sq, flr, ftl):
    out = f'data/stage7j_fullarmw_photow3_{tag}.txt'
    if not os.path.exists(out):
        return None
    txt = open(out).read()
    ms = INJRE.findall(txt)
    ok = len(ms) >= 1 and all(
        abs(float(m[1])-alpha) < 1e-9 and abs(float(m[2])-fcm) < 1e-9
        and abs(float(m[3])-fpm) < 1e-9 and abs(float(m[4])-sq) < 1e-9
        and abs(float(m[5])-flr) < 1e-9 and abs(float(m[6])-ftl) < 1e-9
        for m in ms)
    nchk = txt.count('injected-hist weighted checksum')
    P(f"G8F-WT [{tag}]: {len(ms)} injection line(s) vs registered "
      f"truth -> {'PASS' if ok else 'FAIL'}; {nchk} checksum "
      f"line(s) {'PASS' if nchk >= len(ms) else 'FAIL'}")
    assert ok and nchk >= len(ms), f'G8F-WT failed for {tag}'
    return True

# ---------------- anchors (verbatim 8D/7jz_read construction) ----
anchors = {}
pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
assert str(pz['version']).startswith('v2c'), \
    "operative anchor requires the v2c certificate npz"
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
gband = (pz['conv_band']/0.30 if 'conv_band' in pz.files
         else np.array([0.33, 1.30]))
GS = np.linspace(float(gband[0]), float(gband[1]), 25)
fhmin, fhmax = fg[sup].min(), fg[sup].max()
lnc = np.full(len(FCOMP), -1e9)
for gi in GS:
    fh_eq = FCOMP/gi
    m = (fh_eq >= fhmin) & (fh_eq <= fhmax)
    cand = np.full(len(FCOMP), -1e9)
    cand[m] = np.interp(fh_eq[m], fg, lp)
    lnc = np.maximum(lnc, cand)
anchors['LANDED-CONV'] = lnc
anchors['FLAT'] = np.zeros(len(FCOMP))
P(f"anchors: LANDED-CONV ln pi = {np.round(lnc, 2).tolist()}; "
  f"FLAT = zeros")
P("")

def read(cb9, lnpi):
    cbp = cb9 + lnpi.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp - m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=tuple(range(1, 9))),
                           1e-300)) + m0
    ima = int(np.argmax(lm))
    am = A_GRID[ima]
    if 0 < ima < 4:
        x = A_GRID[ima-1:ima+2]
        y = lm[ima-1:ima+2]
        c2, c1, _ = np.polyfit(x, y, 2)
        if c2 < 0:
            am = -c1/(2*c2)
    post = {}
    for name, ax in (('fcomp', 3), ('fpm', 6), ('sq', 8)):
        m = ex.sum(axis=tuple(i for i in range(9) if i != ax))
        post[name] = m/max(m.sum(), 1e-300)
    return am, float(lm.max()-lm[0]), post

def read_cfg(tag, seed):
    got = {}
    for law in LAWS:
        cp = f'data/stage7j_cube_fullarmw_photow3_{tag}_{seed}_{law}.npy'
        if not os.path.exists(cp):
            return None
        cw = np.load(cp)
        cb9 = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
        for an in ('LANDED-CONV', 'FLAT'):
            am, dn, post = read(cb9, anchors[an])
            got[(law, an)] = (am, dn, post)
            if an == 'LANDED-CONV':
                P(f"[{tag} {law} {seed}] {an}: a_marg = {am:.2f}, "
                  f"dN = {dn:+.1f}; P(fpm=3.0) = "
                  f"{post['fpm'][-1]:.2f}; fcomp mode = "
                  f"{FCOMP[int(np.argmax(post['fcomp']))]:.2f}; "
                  f"P(sq) = {np.round(post['sq'], 2).tolist()}")
            else:
                P(f"[{tag} {law} {seed}] {an}: a_marg = {am:.2f}, "
                  f"dN = {dn:+.1f}")
    return got

# ---------------- F0: gates + reads + calibration ----------------
f0 = {}
for tag, (al, fcm, fpm, sq, flr, ftl) in F0.items():
    if gate_wt(tag, al, fcm, fpm, sq, flr, ftl):
        r = read_cfg(tag, 31)
        if r is not None:
            f0[tag] = r
P("")
missing = [t for t in F0 if t not in f0]
if missing:
    P(f"F0 INCOMPLETE — missing configs: {missing}; calibration "
      f"and verdict PENDING")
else:
    p3 = {t: 0.5*(f0[t][('simple', 'LANDED-CONV')][2]['fpm'][-1]
                  + f0[t][('BE', 'LANDED-CONV')][2]['fpm'][-1])
          for t in PURE}
    P("SYMPTOM TABLE (P3 = law-mean LANDED-CONV mass on fpm = 3.0; "
      "sky: 0.54/0.97):")
    for t in PURE:
        P(f"  {t} (ftl = {F0[t][5]:.2f}): P3 = {p3[t]:.3f}")
    tstar = min(PURE, key=lambda t: abs(p3[t]-0.75))
    matched = p3[tstar] >= 0.30
    P(f"CALIBRATED MEMBER: {tstar} (ftl* = {F0[tstar][5]:.2f}), "
      f"P3 = {p3[tstar]:.3f} -> "
      f"{'SYMPTOM-MATCHED' if matched else 'NON-MATCHED'}")
    if not matched:
        tmax = 't35'
        P(f"CHASE-UNREPRODUCIBLE-BY-TAILS (KT = 4 class) — verdict "
          f"leg on max-ftl member {tmax}, labeled NON-MATCHED")
        tstar = tmax
    ftl_s = F0[tstar][5]
    P(f"control (tsq, ftl = 0.10 + sq_true = 0.2): P3 = "
      f"{0.5*(f0['tsq'][('simple', 'LANDED-CONV')][2]['fpm'][-1] + f0['tsq'][('BE', 'LANDED-CONV')][2]['fpm'][-1]):.3f}")
    P("")

    # ------------- verdict legs (need tstar@101 + fb + fc) -------
    fbr = {'fb': (0.74, 0.10, 1.2, 0.2, 0.0, ftl_s, 'simple'),
           'fc': (0.70, 0.10, 1.2, 0.2, 0.0, ftl_s, 'BE')}
    legs = {}
    g101 = read_cfg(tstar, 101)
    if g101 is not None:
        legs['t101'] = g101
    for tag, (al, fcm, fpm, sq, flr, ftl, lawn) in fbr.items():
        if os.path.exists(f'data/stage7j_fullarmw_photow3_{tag}.txt'):
            if gate_wt(tag, al, fcm, fpm, sq, flr, ftl):
                r = read_cfg(tag, 31)
                if r is not None:
                    legs[tag] = r
    P("")
    need = [k for k in ('t101', 'fb', 'fc') if k not in legs]
    if need:
        P(f"VERDICT PENDING — awaiting legs: {need} (launch fb/fc "
          f"with WTRUTH ftl* = {ftl_s:.2f} per pre-reg)")
    else:
        fn = {(law, s): f0[tstar][(law, 'LANDED-CONV')] if s == 31
              else legs['t101'][(law, 'LANDED-CONV')]
              for law in LAWS for s in (31, 101)}
        man = any(v[0] >= 0.5 and v[1] >= 10.0 for v in fn.values())
        cln_n = all(v[0] <= 0.3 for v in fn.values())
        am_b = legs['fb'][('simple', 'LANDED-CONV')][0]
        am_c = legs['fc'][('BE', 'LANDED-CONV')][0]
        rec = (abs(am_b-0.74) <= 0.25) and (abs(am_c-0.70) <= 0.25)
        P(f"F-N (calibrated {tstar}, both seeds): " + "; ".join(
            f"{law}@{s} a={v[0]:.2f}/dN={v[1]:+.1f}"
            for (law, s), v in fn.items()))
        P(f"F-B recovery (truth 0.74): a_marg = {am_b:.2f}; "
          f"F-C recovery (truth 0.70): a_marg = {am_c:.2f}")
        if man:
            v = ("B-MAN MANUFACTURE-FIRED: the operative alpha "
                 "instrument is tail-breakable — CREDENCE MAP "
                 "(binary ledger): ~50% -> ~40%; alpha rows "
                 "annotated tail-conditional")
        elif cln_n and rec:
            v = ("B-CLEAN MANUFACTURE-EXCLUDED (KT = 4 class): " +
                 ("symptom-matched -> CREDENCE MAP ~50% -> ~55% "
                  "(capped; T2/T3 remain)" if matched else
                  "NON-MATCHED -> HOLD ~50%; W4 NARROWED (the tail "
                  "class cannot even produce the chase) — a "
                  "finding, not a credence event"))
        else:
            v = "UNRESOLVED-CARRIED: hold ~50%; report per grammar"
        P(f"==> 8F VERDICT (locked bars): {v}")

with open('data/stage8f_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8f_read.txt")
