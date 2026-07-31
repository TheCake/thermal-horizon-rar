"""STAGE 8F-c reader — THE BIAS CURVE (pre-reg in NOTES
2026-07-31, committed BEFORE any run).

Maps the own-law recovery bias b(ftl) = a_hat - a_truth on
boosted-truth fat-tail skies (ftl 0..0.35; the 0.35 end = 8F's
fb/fc, joined not rerun).  Bars C-VANISH / C-PRESENT /
GRAY-CARRIED + edge rule + baseline rule as locked in NOTES.
G8F-WT config certificates reused.  NO credence move in any
branch (pre-stated).
Output: data/stage8fc_read.txt
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
FCOMP = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
LAWS = ('simple', 'BE')

# registry: tag -> (truth law, alpha_truth, ftl); all carry
# fcomp 0.10, fpm 1.2, sq 0.2, flr 0 (the 8F fb/fc convention)
CFG = [
    ('fb00', 'simple', 0.74, 0.00), ('fb05', 'simple', 0.74, 0.05),
    ('fb10', 'simple', 0.74, 0.10), ('fb20', 'simple', 0.74, 0.20),
    ('fb',   'simple', 0.74, 0.35),
    ('fc00', 'BE', 0.70, 0.00), ('fc05', 'BE', 0.70, 0.05),
    ('fc10', 'BE', 0.70, 0.10), ('fc20', 'BE', 0.70, 0.20),
    ('fc',   'BE', 0.70, 0.35),
]
REAL = (0.05, 0.10)          # the realistic-zone members
INJRE = re.compile(
    r"WIDTH-SHAPE injection (\w+) alpha=([\d.]+) at fcomp=([\d.]+), "
    r"fpm=([\d.]+), sq=([\d.]+), flr=([\d.]+), ftl=([\d.]+)")

P("8F-c THE BIAS CURVE — reader (pre-reg in NOTES 2026-07-31; "
  "bars locked before any run; no credence move in any branch)")
P("")

def gate_wt(tag, lawn, alpha, ftl):
    out = f'data/stage7j_fullarmw_photow3_{tag}.txt'
    if not os.path.exists(out):
        return False
    txt = open(out).read()
    ms = INJRE.findall(txt)
    ok = len(ms) >= 1 and all(
        m[0] == lawn and abs(float(m[1])-alpha) < 1e-9
        and abs(float(m[2])-0.10) < 1e-9 and abs(float(m[3])-1.2) < 1e-9
        and abs(float(m[4])-0.2) < 1e-9 and abs(float(m[5])-0.0) < 1e-9
        and abs(float(m[6])-ftl) < 1e-9 for m in ms)
    nchk = txt.count('injected-hist weighted checksum')
    P(f"G8F-WT [{tag}]: {len(ms)} injection line(s) -> "
      f"{'PASS' if ok else 'FAIL'}; checksums "
      f"{'PASS' if nchk >= len(ms) else 'FAIL'}")
    assert ok and nchk >= len(ms), f'G8F-WT failed for {tag}'
    return True

pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
assert str(pz['version']).startswith('v2c')
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
ANCH = {'LANDED-CONV': lnc, 'FLAT': np.zeros(len(FCOMP))}

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
        mm = ex.sum(axis=tuple(i for i in range(9) if i != ax))
        post[name] = mm/max(mm.sum(), 1e-300)
    return am, float(lm.max()-lm[0]), post

def own_read(tag, lawn, seed):
    cp = f'data/stage7j_cube_fullarmw_photow3_{tag}_{seed}_{lawn}.npy'
    if not os.path.exists(cp):
        return None
    cw = np.load(cp)
    cb9 = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
    out = {}
    for an in ('LANDED-CONV', 'FLAT'):
        out[an] = read(cb9, ANCH[an])
    # cross-law diagnostic at LANDED-CONV
    other = 'BE' if lawn == 'simple' else 'simple'
    cpo = f'data/stage7j_cube_fullarmw_photow3_{tag}_{seed}_{other}.npy'
    xam = None
    if os.path.exists(cpo):
        co = np.load(cpo) + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
        xam = read(co, ANCH['LANDED-CONV'])[0]
    am, dn, post = out['LANDED-CONV']
    P(f"[{tag} {lawn}@{seed}] LC: a_hat = {am:.2f} (FLAT "
      f"{out['FLAT'][0]:.2f}), dN = {dn:+.1f}; P(fpm=3.0) = "
      f"{post['fpm'][-1]:.2f}; sq mode = "
      f"{[0.0, 0.1, 0.2, 0.3][int(np.argmax(post['sq']))]:.1f}; "
      f"fcomp mode = {FCOMP[int(np.argmax(post['fcomp']))]:.2f}; "
      f"cross-law a_hat = {('%.2f' % xam) if xam is not None else 'n/a'}")
    return am

# ---------------- the curve ----------------
bias = {}
pending = []
for tag, lawn, atr, ftl in CFG:
    if not gate_wt(tag, lawn, atr, ftl):
        pending.append(tag)
        continue
    ams = []
    for seed in (31, 101):
        am = own_read(tag, lawn, seed)
        if am is not None:
            ams.append(am)
    if not ams:
        pending.append(tag)
        continue
    bias[(lawn, ftl)] = (float(np.mean(ams)) - atr, len(ams))
P("")
if pending:
    P(f"CURVE INCOMPLETE — pending: {pending}; verdict PENDING")
else:
    P("THE BIAS CURVE b(ftl) = own-law a_hat - truth (LANDED-CONV; "
      "mean over available seeds):")
    for lawn in LAWS:
        row = "; ".join(f"ftl={f:.2f}: {bias[(lawn, f)][0]:+.2f}"
                        f"{'*' if bias[(lawn, f)][1] > 1 else ''}"
                        for _, lw, _, f in
                        [(t, lw, a, f) for t, lw, a, f in CFG
                         if lw == lawn])
        P(f"  {lawn}: {row}   (* = two-seed mean)")
    P("")
    # baseline rule
    eff = {}
    for lawn in LAWS:
        b0 = bias[(lawn, 0.00)][0]
        attr = abs(b0) > 0.15
        for f in REAL:
            eff[(lawn, f)] = (bias[(lawn, f)][0] - b0 if attr
                              else bias[(lawn, f)][0])
        P(f"baseline b(0) [{lawn}] = {b0:+.2f}"
          + (" -> |b(0)| > 0.15: ATTRIBUTION MODE (bars on "
             "delta-b, disclosed)" if attr else " (raw-b bars)"))
    # edge rule
    edges = []
    for f in REAL:
        for lawn, thr in (('simple', (0.15,)), ('BE', (0.15, 0.25))):
            d = min(abs(abs(eff[(lawn, f)]) - t) for t in thr)
            if d < 0.05 and bias[(lawn, f)][1] < 2:
                edges.append((lawn, f, d))
    if edges:
        P(f"EDGE RULE FIRED: {[(l, f) for l, f, _ in edges]} within "
          f"0.05 of a threshold at single-seed grade -> seed-101 "
          f"confirmation of those members REQUIRED before the "
          f"verdict is quoted (pre-registered conditional "
          f"extension); verdict WITHHELD this pass")
    else:
        van = all(abs(eff[(lawn, f)]) <= 0.15
                  for lawn in LAWS for f in REAL)
        pres = any(eff[('BE', f)] >= 0.25 for f in REAL)
        if pres:
            v = ("C-PRESENT: the coupling operates at realistic "
                 "severities — ANNOTATION consequence: BE alpha "
                 "rows carry the tail-exposure caveat; the "
                 "simple-law alpha is the better-conditioned "
                 "amplitude. NO credence move (pre-stated).")
        elif van:
            v = ("C-VANISH: the coupling is a heavy-tail-only "
                 "phenomenon — sky-relevant amplitude risk bounded "
                 "at the systematic scale. NO credence move "
                 "(pre-stated).")
        else:
            v = ("GRAY-CARRIED: the curve is the product, quoted "
                 "with both flags. NO credence move (pre-stated).")
        P(f"==> 8F-c VERDICT (locked bars): {v}")

with open('data/stage8fc_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8fc_read.txt")
