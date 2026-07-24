"""
STAGE 6Y (O5: the dS-side reservoir identification): WHAT lends the
quanta -- and the exclusion theorem that pins it.

THE IDENTIFICATION. The 6X toy needed a thermal partner mode. On the
sky the candidates are: (a) ONE collective ambient mode per system, or
(b) the ambient field's MANY modes lending democratically. The gate
data decide this EXACTLY:
  - single collective mode (thermal): P(n >= L) = s^L  (GL, 6X) --
    the measured form;
  - M democratic thermal modes: the total quanta are negative-binomial;
    P(N_tot >= 2) = 1 - (1-q)^M (1 + M q), q = e^(-x_amb). Already at
    M = 2 the galaxy gate saturates (0.95 vs the measured 0.75) and
    the e_N-dependence that the measured Chae gates IMPROVE the fits
    with (6I) washes out; the binary gate rises to 0.26 (p-postdiction
    0.565, at the band edge the binaries reject).
  => THE MEASURED GATE EXCLUDES DEMOCRATIC MULTIMODE LENDING: one
     collective ambient mode per system. This is the same system-level
     scalar structure the sky measured dynamically (6D pointwise
     excluded; 6T local excluded; 6G system-level accepted) -- two
     independent routes to one structure.

WHAT the collective mode is (reading, stated plainly): the system's
barycentric coordinate in the ambient field -- one degree of freedom
per system, frequency at the ambient scale (the program's standing
omega = g/c assignment, a postulate not re-derived here), occupation
n_amb = n_BE(x_amb) = nu(e_N) - 1: THE ENVIRONMENT'S OWN DRESSING
CLOUD. Systems borrow from their environment's boost. Continuity:
e_N -> 0 sends the collective mode into the horizon's soft (IR)
sector, n -> infinity, gate -> 1 -- the isolated limit rejoins the
one-bath picture; a strong ambient GAPS the soft sector and closes
the gate.

THE SATURN COROLLARY (reading-grade, three measured legs): a coupling
to the barycentric collective mode is a TRAJECTORY-STATE coupling --
the modified-inertia class -- and that class (i) ties vector-MG on the
binaries (4L: mi_t -3.5+-3.3), (ii) satisfies the 6W composition
demand at current resolving power for the same reason, and (iii)
carries no capped-type solar quadrupole (the 4K escape). The borrowing
mechanics lands the thermal rule in the ONE door Saturn left open --
by derivation, not preference. Kill test: DR4 eccentricity-resolved
boosts (MI: trajectory-dependent; MG: not).

PREDICTIONS extracted (the falsifier ledger):
  P1 the tail-exponent CEILING: p = 1/2 + gate/4 <= 3/4 ALWAYS (gate
     = P(n>=2) <= 1); practical void asymptote p ~ 0.72 at e ~ 0.005.
     One population-grade galaxy tail beyond 3/4 kills the gate.
  P2 the ambient ORDERING: p(e_N) monotone decreasing -- tested NOW
     at full power by the 6Z shuffle test (pre-registered with this
     commit).
  P3 the z-rung (6U): p runs with H(z) at fixed environment.
  P4 DR4: weak-ambient binaries sharpen toward p ~ 0.69; the
     source-vs-dressed convention splits (dp ~ 0.03).
  P5 single-collective-mode structure: no partial-gate dilution in
     systems with resolved substructure (nested ambients gate by the
     LOCAL total field -- operationally what 6G/6I already used).

Gates: GY1 negative-binomial tail identity (sympy exact); GY2 the
M-scan table (numeric); GY3 the ceiling and ordering curves.
Writes data/stage6y_reservoir.txt.
"""
import math
import sympy as sp

L = []
def say(s=''):
    L.append(s); print(s, flush=True)

say("STAGE 6Y: the reservoir identification -- exclusion theorem + "
    "predictions")
say("=" * 72)

# ---- GY1: the multimode tail identity (exact) ------------------------------
q, M = sp.symbols('q M', positive=True)
n = sp.symbols('n', integer=True, nonnegative=True)
# M iid geometric modes: P(N=0) = (1-q)^M; P(N=1) = M q (1-q)^M
P0 = (1 - q)**M
P1 = M*q*(1 - q)**M
Pge2 = sp.simplify(1 - P0 - P1)
# check against the negative-binomial pmf at k=0,1: NB(k; M, q)
nb0 = sp.binomial(M - 1, 0)*(1 - q)**M
nb1 = sp.binomial(M, 1)*(1 - q)**M*q
okA = sp.simplify(P0 - nb0) == 0 and sp.simplify(P1 - nb1) == 0
okB = sp.simplify(Pge2.subs(M, 1) - q**2) == 0
say(f"GY1: P(N_tot >= 2 | M geometric modes) = 1 - (1-q)^M (1+Mq): "
    f"NB check {'PASS' if okA else 'FAIL'}; M=1 reduces to q^2 = s^2: "
    f"{'PASS' if okB else 'FAIL'}")
assert okA and okB

# ---- GY2: the M-scan at the measured ambients ------------------------------
def gate_M(x, Mm):
    qq = math.exp(-x)
    return 1.0 - (1.0 - qq)**Mm * (1.0 + Mm*qq)
X_GAL = math.sqrt(0.02)     # the 6E/6F fiducial galaxy ambient
X_BIN = math.sqrt(1.2)
say('')
say("GY2: the democratic-multimode gate vs M (measured M=1 values: "
    "galaxy 0.754, binary 0.112):")
say(f"    {'M':>3} {'gate_gal':>9} {'p_gal':>7} {'gate_bin':>9} {'p_bin':>7}")
ok2 = True
for Mm in (1, 2, 3, 5, 10):
    gg, gb = gate_M(X_GAL, Mm), gate_M(X_BIN, Mm)
    say(f"    {Mm:3d} {gg:9.4f} {0.5+gg/4:7.4f} {gb:9.4f} {0.5+gb/4:7.4f}")
    if Mm == 1:
        ok2 &= abs(gg - 0.7536) < 1e-3 and abs(gb - 0.1118) < 1e-3
say(f"GY2 (M=1 regression to the measured gate): {'PASS' if ok2 else 'FAIL'}")
assert ok2
say("    => M >= 2: the galaxy gate saturates (0.95+) and the measured")
say("       e_N-dependence (6I: per-galaxy gates IMPROVE every treatment)")
say("       washes out; the binary p-postdiction moves to the band edge")
say("       the binaries reject. THE MEASURED GATE SELECTS M = 1: one")
say("       collective ambient mode per system -- the same system-level")
say("       scalar structure measured dynamically (6D/6T/6G).")

# ---- GY3: the ceiling + ordering curves ------------------------------------
say('')
say("GY3: predictions P1/P2 -- the tail ceiling and the ambient ordering:")
say(f"    p_max = 1/2 + 1/4 = 0.75 EXACT (gate <= 1). The curve p(e_N):")
for e in (0.001, 0.005, 0.02, 0.05, 0.1, 0.4, 1.2):
    x = math.sqrt(e)
    s = math.exp(-x)
    p = 0.5 + 0.25*s*s
    say(f"      e_N = {e:6.3f} a0: gate = {s*s:.4f}, p = {p:.4f}")
say("    -> void asymptote ~0.72; monotone DECREASING in e_N (the 6Z")
say("       shuffle test = this ordering at full sample power); ANY")
say("       population-grade tail beyond 3/4 kills the gate.")

say('')
say("THE SATURN COROLLARY (reading-grade; three measured legs): the")
say("collective-mode coupling is a trajectory-state (modified-inertia-")
say("class) structure; that class ties vector-MG on the binaries (4L),")
say("meets the 6W composition demand at current resolving power, and")
say("carries no capped-type quadrupole (4K). The borrowing mechanics")
say("lands the thermal rule in the one door Saturn left open -- by")
say("derivation, not preference. Kill test: DR4 eccentricity-resolved")
say("boosts (MI trajectory-dependence vs MG).")

with open('data/stage6y_reservoir.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nSTAGE 6Y done -> data/stage6y_reservoir.txt")
