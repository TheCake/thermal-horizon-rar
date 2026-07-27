# -*- coding: utf-8 -*-
"""Stage 7J-z5 A-D arm reader (bars pre-registered f6916e0; the
chase-reading bands below are pre-stated HERE, committed before the
arm batch runs).

ARMS (twin-mismatch injections, photow3 cubes, seed 31; extension to
101 pre-authorized if any bar lands within 0.1 of its edge):
  fullarma  null:   alpha_true = 0.00, fcomp 0.20, fpm 2.4
  fullarmb  simple: alpha_true = 0.74, fcomp 0.10, fpm 2.1
  fullarmc  BE:     alpha_true = 0.70, fcomp 0.10, fpm 2.4
  fullarmd  BE:     alpha_true = 0.40, fcomp 0.10, fpm 2.4
All: sq_true 0.2, wr 0.2, eta 1.3, twin-t5 companions, fc/ff
0.10/0.05, truth kw = 1.0.  Read at the OPERATIVE (LANDED-CONV)
anchor, mirroring amendment 11.

PRE-REGISTERED BARS (per arm, truth-law row):
  A: alpha_marg <= 0.3 in BOTH laws -> ARM-PASS else ARM-FAIL
     (false-positive direction).
  B: |alpha_marg(simple) - 0.74| <= 0.25 AND dN(simple) >= +10.
  C: |alpha_marg(BE) - 0.70| <= 0.25 AND dN(BE) >= +10.
  D: |alpha_marg(BE) - 0.40| <= 0.25 (PROF-vs-MARG offset reported).

THE CHASE QUESTION (pre-stated bands; the E-read found the sky
chasing the fpm edge to 3.0 at P = 0.54/0.97 despite the physical
ceiling ~1.4):
  CHASE-REPRODUCED  if the truth-law P(fpm = 3.0) >= 0.5 in >= 2 of
      the 4 arms (injected fpm_true was 2.1-2.4; twin-mismatch skies
      reproducing the ride attribute the sky's chase to the twin
      population + width mismatch = the fingerprint).
  CHASE-UNEXPLAINED if P(fpm = 3.0) < 0.5 in ALL arms (the sky's
      chase is then real missing structure beyond the twin
      population; the width-shape refinement is promoted).
  CHASE-PARTIAL     otherwise.
Also reported per arm: P(sq) (does the fitter recover the injected
0.2?), P(fcomp), and the alpha-exposure delta between the 6-node and
5-node (<= 2.4) reads (the arm-level analogue of the E-read's ~8 lnL
cession).

Output: data/stage7jz5_armread.txt
"""
import numpy as np

OUT = []
def P(s):
    print(s); OUT.append(s)

A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID = np.array([1.05, 1.3])
FCOMP = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2

pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
fg, lp = pz['fh_grid'], pz['lnpi_host']
sup = lp > -1e8
gband = pz['conv_band']/0.30
GS = np.linspace(float(gband[0]), float(gband[1]), 25)
fhmin, fhmax = fg[sup].min(), fg[sup].max()
lnc = np.full(len(FCOMP), -1e9)
for gi in GS:
    fh_eq = FCOMP/gi
    m = (fh_eq >= fhmin) & (fh_eq <= fhmax)
    cand = np.full(len(FCOMP), -1e9)
    cand[m] = np.interp(fh_eq[m], fg, lp)
    lnc = np.maximum(lnc, cand)

def refine(x, y):
    i = int(np.nanargmax(y)); v = x[i]
    if 0 < i < len(x)-1:
        c2, c1, _ = np.polyfit(x[i-1:i+2], y[i-1:i+2], 2)
        if c2 < 0:
            v = -c1/(2*c2)
    return v, i

def read(cb9):
    cbp = cb9 + lnc.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1))
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp-m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=(1, 2, 3, 4, 5, 6, 7, 8)),
                           1e-300)) + m0
    am, _ = refine(A_GRID, lm)
    pf = ex.sum(axis=(0, 1, 2, 3, 4, 5, 7, 8))
    ps = ex.sum(axis=(0, 1, 2, 3, 4, 5, 6, 7))
    pc = ex.sum(axis=(0, 1, 2, 4, 5, 6, 7, 8))
    prof = np.nanmax(cb9, axis=(1, 2, 3, 4, 5, 6, 7, 8))
    ah, _ = refine(A_GRID, prof)
    return (am, float(lm.max()-lm[0]), pf/pf.sum(), ps/ps.sum(),
            pc/pc.sum(), ah)

ARMS = {
    'fullarma': ('null', 0.00, None, 0.20),
    'fullarmb': ('simple', 0.74, 'simple', 0.10),
    'fullarmc': ('BE', 0.70, 'BE', 0.10),
    'fullarmd': ('BE', 0.40, 'BE', 0.10),
}
SEED = 31
res = {}
edge_flag = False
for arm, (lawn, atr, tlaw, ftr) in ARMS.items():
    for law in ('simple', 'BE'):
        c6 = np.load(f'data/stage7j_cube_{arm}_photow3_{SEED}_{law}.npy')
        cb9 = c6 + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
        am, dn, pf, ps, pc, ah = read(cb9)
        am5, dn5, *_ = read(cb9[:, :, :, :, :, :, :5])
        res[(arm, law)] = (am, dn, pf, ps, pc, ah, am5, dn5)
        P(f"[{arm} {law}] a_marg={am:.2f} dN={dn:+.1f} (<=2.4: "
          f"{am5:.2f}/{dn5:+.1f}) PROF a_hat={ah:.2f} | "
          f"P(fpm)={np.round(pf, 2).tolist()} "
          f"P(sq)={np.round(ps, 2).tolist()} "
          f"P(fcomp)={np.round(pc, 2).tolist()}")

P("")
verdicts = {}
a_ok = all(res[('fullarma', law)][0] <= 0.3 for law in ('simple', 'BE'))
verdicts['A'] = 'ARM-PASS' if a_ok else 'ARM-FAIL (false-positive)'
margin_a = max(res[('fullarma', law)][0] for law in ('simple', 'BE'))
b = res[('fullarmb', 'simple')]
verdicts['B'] = ('ARM-VALIDATED' if abs(b[0]-0.74) <= 0.25 and
                 b[1] >= 10 else 'ARM-FAIL')
c = res[('fullarmc', 'BE')]
verdicts['C'] = ('ARM-VALIDATED' if abs(c[0]-0.70) <= 0.25 and
                 c[1] >= 10 else 'ARM-FAIL')
d = res[('fullarmd', 'BE')]
verdicts['D'] = ('ARM-VALIDATED' if abs(d[0]-0.40) <= 0.25
                 else 'ARM-FAIL')
for k, v in verdicts.items():
    P(f"ARM {k}: {v}")
edges = [abs(margin_a-0.3) < 0.1, abs(abs(b[0]-0.74)-0.25) < 0.1,
         abs(abs(c[0]-0.70)-0.25) < 0.1, abs(abs(d[0]-0.40)-0.25) < 0.1]
if any(edges):
    P("EXTENSION RULE: a bar landed within 0.1 of its edge -> seed 101"
      " pre-authorized (run before quoting that arm)")

chase = [res[(arm, ARMS[arm][2] or 'simple')][2][-1]
         for arm in ARMS]
n_ride = sum(1 for x in chase if x >= 0.5)
tag = ('CHASE-REPRODUCED' if n_ride >= 2 else
       'CHASE-UNEXPLAINED' if n_ride == 0 else 'CHASE-PARTIAL')
P(f"CHASE: truth-law P(fpm=3.0) per arm = "
  f"{[round(float(x), 2) for x in chase]} -> {tag}")
P(f"D PROF-vs-MARG: a_hat={d[5]:.2f} vs a_marg={d[0]:.2f} "
  f"(offset {d[5]-d[0]:+.2f})")

with open('data/stage7jz5_armread.txt', 'w') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage7jz5_armread.txt")
