"""STAGE 8G reader — THE E-SECTOR CONTROL (pre-reg fd33fe0, committed
BEFORE any run; the width-object program, leg 1, mundane-first).

Question: does freeing the 7J-z6-named e-sector (inner power-law
anchors x EIN, radial floor ERF) absorb the sq = 0.2 width demand?

Gates: G8G-1 (identity-slice LANDED-CONV read vs the shipped E-arm
EXT rows, parsed at runtime, |da| <= 0.01, |ddN| <= 0.1); G8G-2
(bookkeeping + correction-#4 edge report on the new axes).  G8G-0
(cube bit-identity) lives in the run itself (abort-grade there).

Bars (locked in the pre-reg; seed means over {31, 101}, both laws):
E-ABSORB P(sq>0) <= 0.50 | E-SURVIVE P(sq>0) >= 0.90 AND mode >= 0.1
| else E-PARTIAL; SEED-SPLIT flagged; |d alpha| > 0.11 -> MATERIAL.
Expected outcome pre-stated: E-ABSORB.  NO credence move any branch.
Output: data/stage8g_read.txt
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
EIN_GRID = np.array([0.5, 1.0, 1.5, 2.0])
ERF_GRID = np.array([0.80, 0.90, 0.95])
FCOMP = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
SQ_GRID = np.array([0.0, 0.1, 0.2, 0.3])
SEEDS = (31, 101)
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2

P("8G THE E-SECTOR CONTROL — reader (pre-reg fd33fe0; bars locked "
  "before any run; expected outcome pre-stated E-ABSORB; no credence "
  "move any branch)")
P("")

# ---- LANDED-CONV anchor (verbatim from the shipped readers) ---------
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
P(f"operative anchor (LANDED-CONV): {np.round(lnc, 2).tolist()}")
P("")

def refine(x, y):
    i = int(np.nanargmax(y)); v = x[i]
    if 0 < i < len(x)-1:
        c2, c1, _ = np.polyfit(x[i-1:i+2], y[i-1:i+2], 2)
        if c2 < 0:
            v = -c1/(2*c2)
    return v

def read_any(cb, sq_ax, extra_ax=()):
    """Marginal read on a prior-augmented cube of any rank.
    Returns a_marg, dN, P(sq), {ax: posterior}, max lnL."""
    nd = cb.ndim
    m0 = np.nanmax(cb)
    ex = np.exp(np.nan_to_num(cb-m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=tuple(range(1, nd))),
                           1e-300)) + m0
    def post(ax):
        o = ex.sum(axis=tuple(i for i in range(nd) if i != ax))
        return o/o.sum()
    return (refine(A_GRID, lm), float(lm.max()-lm[0]), post(sq_ax),
            {ax: post(ax) for ax in extra_ax}, float(m0))

# ---- G8G-1 reference: the shipped E-arm EXT rows --------------------
EXTRE = re.compile(r"\[(\w+) (\d+)\] GE0 .*EXT: a_marg=([\d.]+) "
                   r"dN=([+-][\d.]+) ")
ship = {}
for m in EXTRE.finditer(open('data/stage7jz5_eread.txt').read()):
    ship[(m.group(1), int(m.group(2)))] = (float(m.group(3)),
                                           float(m.group(4)))
assert len(ship) == 4, f'expected 4 shipped EXT rows, got {ship}'

# ---- reads ----------------------------------------------------------
res = {}          # (law, seed) -> dict
g1_ok, pend = True, []
for law in ('simple', 'BE'):
    for seed in SEEDS:
        cp = f'data/stage7j_cube_full_esec_{seed}_{law}.npy'
        if not os.path.exists(cp):
            pend.append(f'{law}@{seed}')
            continue
        c11 = np.load(cp)
        assert c11.shape == (5, 2, 5, 4, 3, 6, 1, 2, 6, 3, 4), c11.shape
        cb11 = (c11
                + prior_eta.reshape((1, 2) + (1,)*9)
                + lnc.reshape((1, 1, 1, 1, 1, 6, 1, 1, 1, 1, 1)))
        # identity slice (ein=1.0 -> index 1, erf=0.90 -> index 1)
        cb9 = cb11[:, :, :, 1, 1]
        a_id, dn_id, psq_id, _, mx_id = read_any(cb9, sq_ax=8)
        sa, sd = ship[(law, seed)]
        ok = abs(a_id-sa) <= 0.01 and abs(dn_id-sd) <= 0.1
        g1_ok &= ok
        P(f"[{law} {seed}] G8G-1: id-slice {a_id:.2f}/{dn_id:+.1f} vs "
          f"shipped {sa:.2f}/{sd:+.1f} -> {'PASS' if ok else 'FAIL'}")
        a_fr, dn_fr, psq_fr, po, mx_fr = read_any(
            cb11, sq_ax=10, extra_ax=(3, 4))
        res[(law, seed)] = dict(
            a_id=a_id, dn_id=dn_id, psq_id=psq_id,
            a_fr=a_fr, dn_fr=dn_fr, psq_fr=psq_fr,
            pein=po[3], perf=po[4], ge=mx_fr-mx_id)
        P(f"[{law} {seed}] ID:    a={a_id:.2f} dN={dn_id:+.1f} "
          f"P(sq)={np.round(psq_id, 2).tolist()}")
        P(f"[{law} {seed}] FREED: a={a_fr:.2f} dN={dn_fr:+.1f} "
          f"P(sq)={np.round(psq_fr, 2).tolist()} "
          f"P(ein)={np.round(po[3], 2).tolist()} "
          f"P(erf)={np.round(po[4], 2).tolist()} "
          f"G_e={mx_fr-mx_id:+.1f}")
P("")

if pend:
    P(f"INCOMPLETE - pending cubes: {pend}; verdict PENDING")
elif not g1_ok:
    P("G8G-1 FAIL - reader wiring suspect; DO NOT QUOTE; verdict "
      "WITHHELD pending inspection")
else:
    # ---- G8G-2: bookkeeping + edge report ---------------------------
    edge = []
    for law in ('simple', 'BE'):
        pein = np.mean([res[(law, s)]['pein'] for s in SEEDS], axis=0)
        perf = np.mean([res[(law, s)]['perf'] for s in SEEDS], axis=0)
        for nm, pp, gr in (('ein', pein, EIN_GRID),
                           ('erf', perf, ERF_GRID)):
            im = int(np.argmax(pp))
            at_edge = im in (0, len(gr)-1) and pp[im] >= 0.5
            if at_edge:
                edge.append(f"{law}:{nm}={gr[im]} (P={pp[im]:.2f})")
        P(f"G8G-2 [{law}]: seed-mean P(ein)="
          f"{np.round(pein, 2).tolist()}, P(erf)="
          f"{np.round(perf, 2).tolist()}")
    P(f"G8G-2 edge report: "
      + (f"EDGE flag on {edge} (correction-#4: extension is a "
         f"decision, not an auto-run)" if edge else
         "no new-axis edge mode with >= 0.5 mass"))
    P("")
    # ---- bars (seed means, both laws) -------------------------------
    branch = {}
    for law in ('simple', 'BE'):
        psq = np.mean([res[(law, s)]['psq_fr'] for s in SEEDS], axis=0)
        pgt0 = float(psq[1:].sum())
        mode = float(SQ_GRID[int(np.argmax(psq))])
        da = float(np.mean([res[(law, s)]['a_fr']-res[(law, s)]['a_id']
                            for s in SEEDS]))
        ge = float(np.mean([res[(law, s)]['ge'] for s in SEEDS]))
        br = ('E-ABSORB' if pgt0 <= 0.50 else
              'E-SURVIVE' if (pgt0 >= 0.90 and mode >= 0.1) else
              'E-PARTIAL')
        # per-seed branch check for the SEED-SPLIT flag
        sbr = []
        for s in SEEDS:
            pq = res[(law, s)]['psq_fr']
            p1 = float(pq[1:].sum())
            m1 = float(SQ_GRID[int(np.argmax(pq))])
            sbr.append('E-ABSORB' if p1 <= 0.50 else
                       'E-SURVIVE' if (p1 >= 0.90 and m1 >= 0.1) else
                       'E-PARTIAL')
        split = len(set(sbr)) > 1
        mat = abs(da) > 0.11
        branch[law] = (br, pgt0, mode, da, ge, split, mat, sbr)
        P(f"BAR [{law}]: P(sq>0)={pgt0:.2f}, sq-mode={mode:.1f}, "
          f"d_alpha={da:+.3f}{' MATERIAL' if mat else ''}, "
          f"G_e={ge:+.1f} -> {br}"
          + (f" [SEED-SPLIT {sbr}]" if split else ""))
    P("")
    brs = {branch[l][0] for l in ('simple', 'BE')}
    if brs == {'E-ABSORB'}:
        v = ("E-ABSORB: the width demand is absorbable by the named "
             "e-sector; sq = e-sector-inadequacy CONFIRMED at control "
             "grade (matches the pre-stated expectation); the freed "
             "alpha/dN are conditional candidates — promotion is a "
             "review-round question, not this stage's verdict")
    elif brs == {'E-SURVIVE'}:
        v = ("E-SURVIVE: the width object survives the freed e-sector "
             "— AGAINST THE PRE-STATED EXPECTATION; leg 2 (8H "
             "boundedness contest) becomes the decider")
    else:
        v = ("E-PARTIAL/mixed: the decomposition is the product; "
             "per-law rows above are the result, no interpretation "
             "beyond the printed numbers")
    if any(branch[l][6] for l in ('simple', 'BE')):
        v += " [alpha MATERIAL flag: the operative quote gains the " \
             "e-sector exposure as a systematic annotation]"
    P(f"==> 8G VERDICT (locked bars, two-seed): {v}")
    P("    NO credence move (pre-stated in every branch).")

with open('data/stage8g_read.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8g_read.txt")
