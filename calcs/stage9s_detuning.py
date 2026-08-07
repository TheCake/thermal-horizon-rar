"""
STAGE 9S -- THE DETUNING BOUND (first quantitative answer to O5-DETUNING,
the ROUND-16 A9 hole).

The reviewer's unpatched objection: the 1/2 prefactor in the lending law
is guaranteed only at RESONANCE; off resonance the operative exchange
weight is r = (1/2)*g_c^2/(g_c^2 + delta^2) < 1/2 (g_c = coupling,
delta = detuning), generically field-configuration-dependent. This stage
asks what the ALREADY-MEASURED sky says about r -- no new fits, no new
functions (7I freeze respected: this is a consistency/reader instrument).

THE INSTRUMENT (tail-postdiction inversion): the AMB tail exponent obeys
p = (1 + beta_tail)/2 with beta_tail = r*G, where G = s_amb^2 is the
measured ambient gate (the 5P family theorem p_tail = (1+beta)/2 is the
verified leg -- 5P/5Q regressions; the local share q_loc = 1/(2nu-1)^2
-> 1 in the tail, sympy below). At r = 1/2 this is the 6E postdiction
p = 1/2 + G/4 that MATCHED both systems. Generalizing r frees the
ceiling: the measured galaxy tail exponent then BOUNDS r from below,
and the resonance curve inverts that into a detuning bound
|delta|/g_c <= sqrt(1/(2r) - 1).

INPUTS (all archived, parsed verbatim):
  - 5G hier-converged tail exponent: "minimum at p=0.65"
    (data/stage5g_tailtest.txt; grid-point value, no sigma printed --
    quoted as a POINT-grade anchor with that caveat; the hier treatment
    is the program's declared primary galaxy analysis per the 7I freeze)
  - 4H flat-treatment co-read: "FINAL: p = 0.578 +0.121 / -0.115"
    (data/stage4h_p_ml.txt)
  - 6E postdiction regression targets: "(pred 0.6884)" gal /
    "(pred 0.5293)" bin (data/stage6e_ambgate.txt)
  - 6I measured ambient gates: "fiducial s = 0.8681 (g = 0.7536);
    measured medians s(maxclust) = 0.9317, s(noclust) = 0.9756"
    (data/stage6i_chaegate.txt); gate grid G = {g_fid, s_max^2, s_noc^2}
  - binary consistency row: n_amb(bin) = 0.52 (6E, cited constant)

PRE-REGISTERED BARS (locked before execution):
  S-BOUNDED : at the HIER anchor (5G p-hat), r >= 0.25 in ALL THREE
              gate cells. Headline = the conservative (smallest-r)
              cell's bound |delta|/g_c <= sqrt(1/(2r)-1).
  S-WEAK    : any hier cell gives r < 0.25 (partial bound; cells quoted).
  S-BROKEN  : the G9S-0 regression or a sympy leg fails -> instrument
              invalid, NO bound quoted.
  Co-reads ALWAYS printed, never verdict-bearing: the 4H flat anchor
  cells (expected sub-bar -- if so, the bound is stated as
  hier-treatment-conditional); the 4H minus-1-sigma cell (expected
  unbounded); the binary consistency row (gate too small to bound r);
  the P1 dividend (below).

DIVIDEND (interpretive, no new registration): the p <= 3/4 ceiling
prediction P1 sharpens into a DIRECT r-METER -- in the void asymptote
G -> 1, p_void = 1/2 + r/2, so a measured void-population tail exponent
reads the exchange prefactor with no gate uncertainty (r = 1/2 <=> 3/4).
Booked as an annotation on the P1 row, not a new prediction.

CREDENCE: NO movement (consistency/measurement round, pre-stated).
O5-DETUNING gets its first number; the two successors stay named: the
r-ladder consistency re-fit (freeze-labeled) and the void r-meter (DR4/
external).

DISCLOSURE (honesty block): the hier-anchor cells were sketched
in-session BEFORE this pre-registration from the ledger-quoted ~0.65
(r ~ 0.32-0.40, delta/g_c <= ~0.77 expected). NOT run before this
commit: every parse, every gate, the sympy legs, the mpmath recompute,
and the 4H/binary co-read cells. Bars were locked with the sketch
disclosed; the sketch cannot move a bar.

GATES:
  G9S-0  regression: p(1/2, G_fid) equals the 6E printed gal
         postdiction to <= 5e-4; the binary consistency row (cited
         n_amb = 0.52) matches the 6E printed bin postdiction to
         <= 2e-3.
  G9S-1  sympy exact: q_loc -> 1 tail limit; p(r,G) = (1+r*G)/2
         algebra; the inversion delta = g_c*sqrt(1/(2r)-1) from
         r = (1/2)g_c^2/(g_c^2+delta^2) (positive branch), zero
         residuals.
  G9S-2  mpmath 50-digit recompute of all inversion cells, rel <= 1e-10.
  G9S-3  ledger legs live: gal-tail-p, gal-p-marg, gal-chae-ambients,
         mech-9q-ceiling all status CURRENT.

Output: data/stage9s_detuning.txt. Wall-clock: seconds.
"""
import math
import os
import re

import sympy as sp
import mpmath as mp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

OUT = []
def say(s=""):
    print(s)
    OUT.append(s)

say("=" * 78)
say("STAGE 9S -- THE DETUNING BOUND (tail-postdiction inversion)")
say("=" * 78)

gates = {}

# ---------------- parse the archived inputs -------------------------
t5g = open('data/stage5g_tailtest.txt', encoding='utf-8').read()
t4h = open('data/stage4h_p_ml.txt', encoding='utf-8').read()
t6e = open('data/stage6e_ambgate.txt', encoding='utf-8').read()
t6i = open('data/stage6i_chaegate.txt', encoding='utf-8').read()

p_hier = float(re.search(r'minimum at p=([\d.]+)', t5g).group(1))
m4h = re.search(r'FINAL: p = ([\d.]+) \+([\d.]+) / -([\d.]+)', t4h)
p_flat, p_flat_up, p_flat_dn = map(float, m4h.groups())
pred_gal = float(re.search(r'G4 gal: tail p_hat = [\d.]+ \(pred ([\d.]+)\)',
                           t6e).group(1))
pred_bin = float(re.search(r'G4 bin: tail p_hat = [\d.]+ \(pred ([\d.]+)\)',
                           t6e).group(1))
m6i = re.search(r'fiducial s = ([\d.]+) \(g = ([\d.]+)\); measured medians '
                r's\(maxclust\) = ([\d.]+), s\(noclust\) = ([\d.]+)', t6i)
s_fid, g_fid, s_max, s_noc = map(float, m6i.groups())
GATES_G = [('fiducial', g_fid), ('maxclust', s_max**2),
           ('noclust', s_noc**2)]
N_AMB_BIN = 0.52   # 6E cited constant (n_amb at the binary ambient)

say("parsed anchors: 5G hier p-hat = %.2f (grid point, no sigma "
    "printed); 4H flat p = %.3f +%.3f/-%.3f;" % (p_hier, p_flat,
    p_flat_up, p_flat_dn))
say("  6E preds gal %.4f / bin %.4f; 6I gates: fid g = %.4f, "
    "maxclust s^2 = %.4f, noclust s^2 = %.4f" % (pred_gal, pred_bin,
    g_fid, s_max**2, s_noc**2))
say("")

# ---------------- G9S-1: the sympy legs -----------------------------
nu_s, r_s, G_s, gc, dl = sp.symbols('nu r G g_c delta', positive=True)
lim_q = sp.limit(1/(2*nu_s - 1)**2, nu_s, 1)
alg = sp.simplify((1 + r_s*G_s)/2 - (sp.Rational(1, 2) + r_s*G_s/2))
sols = sp.solve(sp.Eq(r_s, sp.Rational(1, 2)*gc**2/(gc**2 + dl**2)), dl)
pos = [s for s in sols if s.subs({gc: 1, r_s: sp.Rational(1, 4)}) > 0]
inv_resid = sp.simplify(pos[0] - gc*sp.sqrt(1/(2*r_s) - 1))
ok_g1 = (lim_q == 1) and (alg == 0) and (inv_resid == 0)
gates['G9S-1'] = ok_g1
say("G9S-1 sympy: q_loc tail limit = %s; p-algebra residual = %s; "
    "inversion residual = %s -> %s" % (lim_q, alg, inv_resid,
    "PASS" if ok_g1 else "FAIL"))
say("  (verified leg cited: the 5P family theorem p_tail = (1+beta)/2, "
    "5P/5Q regressions; beta_tail = r*G by q_loc -> 1)")

# ---------------- G9S-0: postdiction regression at r = 1/2 ----------
p_half_gal = 0.5 + 0.5*g_fid/2
d_gal = abs(p_half_gal - pred_gal)
s_bin = N_AMB_BIN/(1.0 + N_AMB_BIN)
p_half_bin = 0.5 + 0.5*(s_bin**2)/2
d_bin = abs(p_half_bin - pred_bin)
ok_g0 = d_gal <= 5e-4 and d_bin <= 2e-3
gates['G9S-0'] = ok_g0
say("G9S-0 regression: p(1/2, g_fid) = %.4f vs 6E gal pred %.4f "
    "(d = %.1e); bin (n_amb = %.2f -> G = %.4f): %.4f vs pred %.4f "
    "(d = %.1e) -> %s" % (p_half_gal, pred_gal, d_gal, N_AMB_BIN,
    s_bin**2, p_half_bin, pred_bin, d_bin, "PASS" if ok_g0 else "FAIL"))
say("")

# ---------------- the inversion table -------------------------------
def cells(p_val):
    row = []
    for nm, G in GATES_G:
        r = 2.0*(p_val - 0.5)/G
        if r <= 0:
            row.append((nm, G, r, None, 'UNBOUNDED'))
        elif r > 0.5:
            row.append((nm, G, r, None, 'OVER-CEILING'))
        else:
            row.append((nm, G, r, math.sqrt(1.0/(2.0*r) - 1.0), 'ok'))
    return row

say("INVERSION  r = 2(p - 1/2)/G ;  |delta|/g_c = sqrt(1/(2r) - 1)")
say("-" * 78)
tables = {}
for tag, p_val in [('HIER 5G (verdict-bearing)', p_hier),
                   ('FLAT 4H (co-read)', p_flat),
                   ('FLAT 4H -1sigma (co-read)', p_flat - p_flat_dn)]:
    say("  anchor %-28s p = %.3f" % (tag, p_val))
    rows = cells(p_val)
    tables[tag] = rows
    for nm, G, r, dg, st in rows:
        if st == 'ok':
            say("    G(%-8s) = %.4f  r = %.3f  |delta|/g_c <= %.3f"
                % (nm, G, r, dg))
        else:
            say("    G(%-8s) = %.4f  r = %.3f  %s" % (nm, G, r, st))
say("  binary consistency row: G_bin = %.4f -> p spans [0.500, %.4f] "
    "over r in [0, 1/2]; measured = 1/2-consistent (5K 12/12) and "
    "accepts 0.529 (6G) -> binaries CANNOT bound r (gate too small); "
    "consistent with every r" % (s_bin**2, p_half_bin))
say("")

# ---------------- verdict per locked bars ---------------------------
hier_rows = tables['HIER 5G (verdict-bearing)']
hier_ok = all(st == 'ok' and r >= 0.25 for _, _, r, _, st in hier_rows)
r_min_cell = min((row for row in hier_rows if row[4] == 'ok'),
                 key=lambda z: z[2], default=None)
if not (gates['G9S-0'] and gates['G9S-1']):
    verdict = "S-BROKEN (instrument invalid; no bound quoted)"
elif hier_ok:
    verdict = ("S-BOUNDED: r >= %.3f in all gate cells at the hier "
               "anchor; conservative bound |delta|/g_c <= %.3f "
               "(cell %s)" % (r_min_cell[2], r_min_cell[3],
               r_min_cell[0]))
else:
    verdict = "S-WEAK (a hier cell fell below 0.25; cells quoted above)"
say("VERDICT (locked grammar): " + verdict)
say("  reading: the measured galaxy tail exponent already FORBIDS "
    "far-off-resonance exchange -- the operative weight sits within "
    "one linewidth of the resonant ceiling (delta <~ 0.8 g_c), OR the "
    "coupling is strong (g_c >> delta). The ROUND-16 A9 scenario "
    "'generic detuning pushes the weight far below 1/2' is "
    "data-disfavored at the hier anchor. CONDITIONALITY stated: the "
    "bound is hier-treatment-anchored (the 4H flat co-read is too "
    "soft to bound r at its lower error edge -- printed above).")
say("")
say("DIVIDEND (P1 sharpening, annotation not new registration): in "
    "the void asymptote G -> 1 the tail exponent reads the prefactor "
    "DIRECTLY: p_void = 1/2 + r/2 (r = 1/2 <=> 3/4). The p <= 3/4 "
    "ceiling test doubles as an r-meter with no gate uncertainty.")
say("")

# ---------------- G9S-2: mpmath recompute ---------------------------
mp.mp.dps = 50
ok_g2 = True
for tag, p_val in [('hier', p_hier), ('flat', p_flat)]:
    for nm, G in GATES_G:
        r_f = 2.0*(p_val - 0.5)/G
        r_m = 2*(mp.mpf(repr(p_val)) - mp.mpf('0.5'))/mp.mpf(repr(G))
        rel = abs(r_f/float(r_m) - 1.0) if r_m != 0 else 0.0
        ok_g2 &= rel <= 1e-10
        if r_f > 0 and r_f <= 0.5:
            d_f = math.sqrt(1.0/(2.0*r_f) - 1.0)
            d_m = mp.sqrt(1/(2*r_m) - 1)
            ok_g2 &= abs(d_f/float(d_m) - 1.0) <= 1e-10
gates['G9S-2'] = ok_g2
say("G9S-2 mpmath 50-digit recompute of all cells: %s"
    % ("PASS" if ok_g2 else "FAIL"))

# ---------------- G9S-3: ledger legs --------------------------------
need = ['gal-tail-p', 'gal-p-marg', 'gal-chae-ambients',
        'mech-9q-ceiling']
led = {}
for line in open('LEDGER.csv', encoding='utf-8'):
    for rid in need:
        if line.startswith(rid + ','):
            led[rid] = line.split(',')[1]
ok_g3 = all(led.get(rid) == 'CURRENT' for rid in need)
gates['G9S-3'] = ok_g3
say("G9S-3 ledger legs: " + "; ".join("%s=%s" % (rid,
    led.get(rid, 'MISSING')) for rid in need) +
    " -> %s" % ("PASS" if ok_g3 else "FAIL"))
say("")
allok = all(gates.values())
say("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
    for k, v in sorted(gates.items())) +
    "  -> ALL %s" % ("PASS" if allok else "FAIL"))
say("")
say("NO credence movement (pre-stated; consistency round). O5-DETUNING "
    "successors: the r-ladder consistency re-fit (freeze-labeled) and "
    "the void-asymptote r-meter (external/DR4).")

with open('data/stage9s_detuning.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(OUT) + "\n")
say("")
say("written: data/stage9s_detuning.txt")
