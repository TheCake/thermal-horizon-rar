# -*- coding: utf-8 -*-
"""Stage 7J-z5-E reader: the fpm -> 3.0 extension read (pre-reg
f6916e0; cubes built behind GB0w/GB0e identity gates, all 0.00e+00).

QUESTIONS (pre-registered): (1) D3 — does the BE noise-grid edge
release with 3.0 available?  RELEASED if P(fpm = 3.0) < 0.5 under the
operative (LANDED-CONV) marginal and the posterior mode is interior;
STILL-RIDING if P(3.0) >= 0.5 (correction-#4 standard: grid-edge =
artifact until the grid extends — a ride at 3.0 queues one further
extension decision, not an automatic run).  (2) alpha exposure: does
alpha_marg move with the extended grid (report; the operative claim
carries |shift| as a systematic note).  (3) the kw attribution at the
extended grid: Dwob = cost(fcomp >= 0.35 | kw = 1.4) - cost(kw = 0.7)
recomputed on photow3 (comparator: +314/+306 at <= 2.4).

GATE GE0 (slice identity, arithmetic): the <= 2.4 sub-grid operative
marginal must reproduce the shipped anchored read exactly
(alpha 0.74/0.74/0.75/0.65, dN +24.4/+23.2/+24.1/+22.4; tol 0.01/0.1)
— the cubes are bit-identical on the shared grid, so any deviation is
reader wiring.  Fail -> inspect, do not quote.

Output: data/stage7jz5_eread.txt
"""
import numpy as np

OUT = []
def P(s):
    print(s); OUT.append(s)

A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID = np.array([1.05, 1.3])
FCOMP = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FPM6 = np.array([1.2, 1.5, 1.8, 2.1, 2.4, 3.0])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
SHIP = {('simple', 31): (0.74, 24.4), ('simple', 101): (0.74, 23.2),
        ('BE', 31): (0.75, 24.1), ('BE', 101): (0.65, 22.4)}

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
P(f"operative anchor (LANDED-CONV): {np.round(lnc, 2).tolist()}")

def refine(x, y):
    i = int(np.nanargmax(y)); v = x[i]
    if 0 < i < len(x)-1:
        c2, c1, _ = np.polyfit(x[i-1:i+2], y[i-1:i+2], 2)
        if c2 < 0:
            v = -c1/(2*c2)
    return v

def read6(cb9):
    cbp = cb9 + lnc.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1))
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp-m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=(1, 2, 3, 4, 5, 6, 7, 8)),
                           1e-300)) + m0
    pf = ex.sum(axis=(0, 1, 2, 3, 4, 5, 7, 8))
    return refine(A_GRID, lm), float(lm.max()-lm[0]), pf/pf.sum()

g_ok = True
rows = []
for law in ('simple', 'BE'):
    for seed in (31, 101):
        c6 = np.load(f'data/stage7j_cube_full_photow3_{seed}_{law}.npy')
        cb9 = c6 + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
        # GE0: truncated-grid identity vs the shipped anchored read
        am5, dn5, _ = read6(cb9[:, :, :, :, :, :, :5])
        sa, sd = SHIP[(law, seed)]
        ok0 = abs(am5-sa) <= 0.01 and abs(dn5-sd) <= 0.1
        g_ok &= ok0
        # extended-grid operative read
        am6, dn6, pf6 = read6(cb9)
        # kw attribution on the extended grid
        cks = []
        for ki in (0, 2):
            sl = cb9[:, :, :, :, :, :, :, ki:ki+1, :]
            cks.append(float(np.nanmax(sl) - np.nanmax(sl[:, :, :, 3:])))
        rows.append((law, seed, am6, dn6, pf6, cks[1]-cks[0]))
        P(f"[{law} {seed}] GE0 {am5:.2f}/{dn5:+.1f} vs shipped "
          f"{sa}/{sd:+.1f} {'PASS' if ok0 else 'FAIL'} | EXT: "
          f"a_marg={am6:.2f} dN={dn6:+.1f} "
          f"P(fpm)={np.round(pf6, 2).tolist()} | Dwob={cks[1]-cks[0]:+.1f}")

P("")
P(f"GATES: {'GE0 ALL PASS' if g_ok else 'GE0 FAIL -- do not quote'}")
if g_ok:
    for law in ('simple', 'BE'):
        rs = [r for r in rows if r[0] == law]
        p30 = np.mean([r[4][-1] for r in rs])
        mode_int = all(int(np.argmax(r[4])) < 5 for r in rs)
        am = np.mean([r[2] for r in rs])
        dn = np.mean([r[3] for r in rs])
        dw = np.mean([r[5] for r in rs])
        sh = am - np.mean([SHIP[(law, s)][0] for s in (31, 101)])
        tag = ('RELEASED' if p30 < 0.5 and mode_int else 'STILL-RIDING')
        P(f"{law}: P(fpm=3.0)={p30:.2f}, mode "
          f"{'interior' if mode_int else 'AT EDGE'} -> D3-EXT {tag}; "
          f"a_marg={am:.2f} (shift {sh:+.03f}) dN={dn:+.1f}; "
          f"Dwob={dw:+.1f}")

with open('data/stage7jz5_eread.txt', 'w') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage7jz5_eread.txt")
