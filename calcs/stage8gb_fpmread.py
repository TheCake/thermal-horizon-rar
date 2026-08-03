"""STAGE 8G-b — the freed noise-posterior read (read-only diagnostic
on the shipped 8G cubes; frame committed BEFORE the run).

MEASUREMENT ONLY, no verdict, no credence move (the cadence rule):
does the freed e-sector release the fpm = 3.0 noise chase at MARGINAL
level (the 8G PROF rows hinted 3-of-4 off the edge)?  The numbers
feed the 8H design premise; any interpretation happens in the 8H
pre-reg, not here.

Gate G8Gb-0: the identity-slice P(fpm) must reproduce the shipped
E-arm P(fpm) rows (data/stage7jz5_eread.txt, parsed at runtime)
elementwise to 0.01 — the cubes are bit-identical on that slice, so
any deviation is reader wiring.  Fail -> do not quote.

Output: data/stage8gb_fpmread.txt
"""
import re
import numpy as np

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID = np.array([1.05, 1.3])
FCOMP = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
FPM6 = [1.2, 1.5, 1.8, 2.1, 2.4, 3.0]
WR5 = [0.10, 0.20, 0.30, 0.40, 0.50]
SEEDS = (31, 101)
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2

P("8G-b THE FREED NOISE-POSTERIOR READ (read-only diagnostic; "
  "MEASUREMENT ONLY - no verdict, no credence move; feeds the 8H "
  "design premise)")
P("")

pz = np.load('data/stage7jz_prior.npz', allow_pickle=True)
assert str(pz['version']).startswith('v2c')
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

# shipped E-arm P(fpm) reference rows
FPRE = re.compile(r"\[(\w+) (\d+)\] GE0 .*P\(fpm\)=\[([0-9., ]+)\]")
ship = {}
for m in FPRE.finditer(open('data/stage7jz5_eread.txt').read()):
    ship[(m.group(1), int(m.group(2)))] = np.array(
        [float(x) for x in m.group(3).split(',')])
assert len(ship) == 4, ship

def post(ex, nd, ax):
    o = ex.sum(axis=tuple(i for i in range(nd) if i != ax))
    return o/o.sum()

g_ok = True
for law in ('simple', 'BE'):
    for seed in SEEDS:
        c11 = np.load(f'data/stage7j_cube_full_esec_{seed}_{law}.npy')
        cb11 = (c11
                + prior_eta.reshape((1, 2) + (1,)*9)
                + lnc.reshape((1, 1, 1, 1, 1, 6, 1, 1, 1, 1, 1)))
        # identity slice: 9-dim (A,E,WR,FCOMP,FC0,FFLY,FPM,KW,SQ)
        cb9 = cb11[:, :, :, 1, 1]
        ex9 = np.exp(np.nan_to_num(cb9 - np.nanmax(cb9), nan=-np.inf))
        pf_id = post(ex9, 9, 6)
        pw_id = post(ex9, 9, 2)
        d0 = float(np.max(np.abs(pf_id - ship[(law, seed)])))
        ok = d0 <= 0.01
        g_ok &= ok
        P(f"[{law} {seed}] G8Gb-0: max|P(fpm)_id - shipped| = {d0:.3f}"
          f" -> {'PASS' if ok else 'FAIL'}")
        ex11 = np.exp(np.nan_to_num(cb11 - np.nanmax(cb11),
                                    nan=-np.inf))
        pf_fr = post(ex11, 11, 8)
        pw_fr = post(ex11, 11, 2)
        P(f"[{law} {seed}] P(fpm) id:    "
          f"{np.round(pf_id, 2).tolist()}  P(3.0)={pf_id[-1]:.2f}")
        P(f"[{law} {seed}] P(fpm) freed: "
          f"{np.round(pf_fr, 2).tolist()}  P(3.0)={pf_fr[-1]:.2f}")
        P(f"[{law} {seed}] P(wr)  id: {np.round(pw_id, 2).tolist()} "
          f"| freed: {np.round(pw_fr, 2).tolist()}")
P("")
if g_ok:
    P("G8Gb-0 ALL PASS (grids: fpm " + str(FPM6) + ", wr "
      + str(WR5) + ")")
    P("MEASUREMENT ONLY - no verdict line by design; interpretation "
      "belongs to the 8H pre-reg (the width-shape premise).")
else:
    P("G8Gb-0 FAIL - reader wiring suspect; DO NOT QUOTE")

with open('data/stage8gb_fpmread.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8gb_fpmread.txt")
