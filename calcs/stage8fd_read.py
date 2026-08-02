"""STAGE 8F-d reader — THE PATCH BLOCK (pre-reg in NOTES
2026-08-02, committed BEFORE any run; the reviewer's H1-H3).

Two-seed hardening of the 8F-c bias curve: member b = mean over
available seeds (half-spread reported); the 8F-c bars re-applied
at two-seed grade; H1 baseline diagnosis + H3 onset grammar per
the pre-reg.  GD0: the seed-31 slice must reproduce the shipped
stage8fc_read.txt values to 0.01.  NO credence move any branch.
Output: data/stage8fd_read.txt
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

# tag -> (truth law, alpha_truth, ftl, expected seeds)
CFG = [
    ('fb00', 'simple', 0.74, 0.00, (31, 101)),
    ('fb05', 'simple', 0.74, 0.05, (31, 101)),
    ('fb10', 'simple', 0.74, 0.10, (31, 101)),
    ('fb20', 'simple', 0.74, 0.20, (31,)),
    ('fb',   'simple', 0.74, 0.35, (31,)),
    ('fc00', 'BE', 0.70, 0.00, (31, 101)),
    ('fc05', 'BE', 0.70, 0.05, (31, 101)),
    ('fc10', 'BE', 0.70, 0.10, (31, 101)),
    ('fc15', 'BE', 0.70, 0.15, (31, 101)),
    ('fc20', 'BE', 0.70, 0.20, (31, 101)),
    ('fc',   'BE', 0.70, 0.35, (31,)),
]
REAL = (0.05, 0.10)
INJRE = re.compile(
    r"WIDTH-SHAPE injection (\w+) alpha=([\d.]+) at fcomp=([\d.]+), "
    r"fpm=([\d.]+), sq=([\d.]+), flr=([\d.]+), ftl=([\d.]+)")
FCRE = re.compile(r"\[(\w+) (\w+)@31\] LC: a_hat = ([\d.]+)")

P("8F-d THE PATCH BLOCK — reader (pre-reg in NOTES 2026-08-02; "
  "bars locked before any run; two-seed grade; no credence move)")
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
ANCH = lnc

def read(cb9):
    cbp = cb9 + ANCH.reshape((1, 1, 1, -1, 1, 1, 1, 1, 1))
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
    return am

# ---------------- reads ----------------
seed31 = {}
bias = {}
half = {}
pending = []
for tag, lawn, atr, ftl, seeds in CFG:
    if not gate_wt(tag, lawn, atr, ftl):
        pending.append(tag)
        continue
    ams = {}
    for seed in seeds:
        cp = f'data/stage7j_cube_fullarmw_photow3_{tag}_{seed}_{lawn}.npy'
        if not os.path.exists(cp):
            continue
        cw = np.load(cp)
        cb9 = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
        ams[seed] = read(cb9)
    if set(seeds) - set(ams):
        pending.append(f"{tag}@{sorted(set(seeds)-set(ams))}")
    if not ams:
        continue
    if 31 in ams:
        seed31[(tag, lawn)] = ams[31]
    vals = np.array(list(ams.values()))
    bias[(lawn, ftl)] = (float(vals.mean()) - atr, len(vals))
    half[(lawn, ftl)] = float((vals.max()-vals.min())/2)
    P(f"[{tag} {lawn}] seeds {sorted(ams)}: a_hat = "
      + ", ".join(f"{ams[s]:.2f}" for s in sorted(ams))
      + f" -> b = {bias[(lawn, ftl)][0]:+.2f} "
        f"(half-spread {half[(lawn, ftl)]:.2f})")
P("")

# ---------------- GD0 regression vs the shipped 8F-c read -------
try:
    ref = open('data/stage8fc_read.txt', encoding='utf-8',
               errors='replace').read()
    refv = {(m[0], m[1]): float(m[2]) for m in FCRE.findall(ref)}
    bad = []
    for (tag, lawn), v in seed31.items():
        if (tag, lawn) in refv and abs(v - refv[(tag, lawn)]) > 0.01:
            bad.append(f"{tag}: {v:.2f} vs shipped {refv[(tag, lawn)]:.2f}")
    P(f"GD0 (seed-31 slice vs shipped 8F-c, {len(refv)} reference "
      f"rows): {'PASS' if not bad else 'FAIL ' + '; '.join(bad)}")
    assert not bad, 'GD0 failed'
except OSError:
    P("GD0: data/stage8fc_read.txt absent - SKIPPED (disclosed)")
P("")

if pending:
    P(f"INCOMPLETE - pending: {pending}; verdict PENDING")
else:
    P("THE TWO-SEED BIAS CURVE b(ftl) (LANDED-CONV; mean over "
      "seeds, half-spread in parens):")
    for lawn in ('simple', 'BE'):
        fts = sorted(f for lw, f in bias if lw == lawn)
        P(f"  {lawn}: " + "; ".join(
            f"ftl={f:.2f}: {bias[(lawn, f)][0]:+.2f}"
            f"({half[(lawn, f)]:.2f})"
            f"{'*' if bias[(lawn, f)][1] > 1 else ''}" for f in fts)
          + "   (* = two-seed)")
    P("")
    # H1 baseline rule (two-seed)
    eff = {}
    attr_note = {}
    for lawn in ('simple', 'BE'):
        b0 = bias[(lawn, 0.00)][0]
        attr = abs(b0) > 0.15
        attr_note[lawn] = attr
        for f in REAL:
            eff[(lawn, f)] = (bias[(lawn, f)][0] - b0 if attr
                              else bias[(lawn, f)][0])
        P(f"H1 [{lawn}]: two-seed b(0) = {b0:+.2f} "
          f"(half-spread {half[(lawn, 0.00)]:.2f}) -> "
          + (f"|b(0)| > 0.15: the config baseline is REAL, "
             f"ATTRIBUTION MODE; the +-0.2 quote-precision "
             f"annotation ADOPTED for this config" if attr else
             f"raw-b bars (fluke component confirmed)"))
    P("")
    # H3 onset grammar
    b15 = bias[('BE', 0.15)][0]
    b20 = bias[('BE', 0.20)][0]
    P(f"H3 onset (two-seed): b_BE(0.15) = {b15:+.2f}"
      f"({half[('BE', 0.15)]:.2f}), b_BE(0.20) = {b20:+.2f}"
      f"({half[('BE', 0.20)]:.2f}); the +0.25 crossing "
      + ("sits between 0.15 and 0.20" if b15 < 0.25 <= b20 else
         "is below 0.15" if b15 >= 0.25 else
         "is above 0.20" if b20 < 0.25 else "n/a"))
    P("")
    # bars at two-seed grade
    flags = []
    for f in REAL:
        for lawn, thr in (('simple', (0.15,)), ('BE', (0.15, 0.25))):
            d = min(abs(abs(eff[(lawn, f)]) - t) for t in thr)
            if d < 0.05:
                flags.append((lawn, f))
    van = all(abs(eff[(lawn, f)]) <= 0.15
              for lawn in ('simple', 'BE') for f in REAL)
    pres = any(eff[('BE', f)] >= 0.25 for f in REAL)
    if pres:
        v = ("C-PRESENT at two-seed grade: annotation consequence "
             "only (BE tail-exposure caveat; simple named the "
             "better-conditioned amplitude)")
    elif van:
        v = ("C-VANISH CONFIRMED-HARDENED at two-seed grade: the "
             "coupling is heavy-tail-only; the 8F-c verdict stands "
             "upgraded from single-seed grade")
    else:
        v = "GRAY at two-seed grade: curve is the product"
    fl = (f" [SEED-LIMITED flag on {flags} - within 0.05 of a "
          f"threshold at two seeds; no third seed per pre-reg]"
          if flags else "")
    P(f"==> 8F-d VERDICT (locked bars, two-seed): {v}{fl}")
    P("    NO credence move (pre-stated in every branch).")

with open('data/stage8fd_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8fd_read.txt")
