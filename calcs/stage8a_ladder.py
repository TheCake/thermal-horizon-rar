"""STAGE 8A — THE RIVALS' LADDER (the fingerprint table).

The rivals arc, part i: classify the field's interpolating-function
catalog and the rival-theory forms by their coefficient-ladder
fingerprint (deep expansion nu = A/x + c1 + c2*x + ..., x = sqrt(g_N/a0)),
and price each against the EXISTING ladder measurements (4S/4Z/5T rows).
Freeze-compliant: consistency rows only — no new fit anywhere in this
script; every quoted constraint is an existing ledger measurement.

Contents (all sympy-exact, gated):
  G1-G3  regressions: BE, simple-nu, Verlinde rungs (note_verlinde_c1).
  G4     THE MASTER INVERSION: for any AQUAL mu(z) = z + m2 z^2 + m3 z^3
         + ..., the ladder is c1 = -m2/2, c2 = 5 m2^2/8 - m3/2.  The
         measured c1 therefore measures mu''(0): the data demand
         m2 in [-1.04, -0.42]; the catalog's standard member has m2 = 0.
  G5     the n-family (Hees nu_alpha): n = 1 (simple) has c1 = 1/2; ALL
         n >= 2 have c1 = 0 EXACTLY (n = 2 = the "standard" function).
  G6     THE IDENTITY: the exponential mu-function mu(z) = 1 - e^(-z)
         is EXACTLY the 4F quantum-bootstrap bath nu = 1 + n_BE(nu*y)
         (both solve e^u = nu/(nu-1), u = nu*y); fingerprint (1/4, 7/96)
         and nu(1) = 1.3500 match 4F's printed cell.  Scout: this bare
         form is apparently NOWHERE named or tested (attribution NOT
         FOUND) — the program's boot adjudication (4F/5C/5F/5M) is then
         apparently its FIRST data contest (scout-level).
  G7     THE ADDITIVE CLASS THEOREM: any "apparent dark matter" model
         g_obs = g_N + sqrt(a0 g_N) * F(g_N/a0) with F analytic at 0,
         F(0) = 1, has c1 = 1 EXACTLY and no even rungs.  Members:
         Verlinde (exact), the superfluid-DM MOND-limit composition
         (membership flagged: primary-source equation check pending),
         any dark component tracking baryons analytically.
  G8     non-RAR structures excluded upstream: constant-acceleration
         additions (Rindler); screened chameleon f(R) (Naik+19: upturn
         signature, no universal deep-MOND slope).
  G9     THE HEES+16 CATALOG (arXiv:1510.01369, the field's function
         zoo): nu-tilde_a = (1-e^-y)^(-1/2) + a e^-y  ->  c1 = a (the
         family's parameter IS the zero-point coefficient), c2 = 1/4;
         nu-bar_a = (1-e^(-y^a))^(-1/(2a)) + (1-1/(2a)) e^(-y^a)  ->
         c1 = 1 - 1/(2a); nu-hat_a = (1-e^(-y^(a/2)))^(-1/a)  ->
         c1 = 1/2 at a = 1 (== BE/RAR-fit EXACTLY, gated) and c1 = 0
         for a >= 2.
  G10    THE PINCER: Hees+16's Cassini verdicts (scout-quoted) vs the
         measured ladder window have EMPTY INTERSECTION across the
         entire published catalog — the amplitude-locked solar tension
         (4K/5S/5I) rederived in coefficient language, catalog-wide.

Measured anchors quoted (existing rows, no new fits): c1 profile
0.385-0.519 (flat M/L, 4S) / 0.208-0.309 (hierarchical, 4Z);
bootstrap ~0.4 +/- 0.3; c1 = 0 excluded Delta(-2lnL) = 56.3
(7.5 sigma profile / 95.5% bootstrap); binary c1 ERASED (7J-d) and
NOT quoted.  Direct-contest cross-references: standard-mu already
contested (4B: 198-200/200 raw, -56 honest); boot cell already
adjudicated (4F raw dead 4-7/200; 5C hier flip +75.6; 5M vertical
collapse -9 = dead on the primary treatment).

Output: data/stage8a_ladder.txt
"""
import sympy as sp

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

x = sp.symbols('x', positive=True)

# ---------------------------------------------------------------- G1-G3
nuBE = 1 + 1/(sp.exp(x) - 1)
sBE = sp.series(nuBE, x, 0, 2).removeO().expand()
assert (sBE.coeff(x, -1), sBE.coeff(x, 0), sBE.coeff(x, 1)) == \
    (1, sp.Rational(1, 2), sp.Rational(1, 12)), sBE
P("G1 BE ladder regression: nu = 1/x + 1/2 + x/12 (exact) -> PASS")

nuS = sp.Rational(1, 2) + sp.sqrt(x**2 + 4)/(2*x)
sS = sp.series(nuS, x, 0, 2).removeO().expand()
assert (sS.coeff(x, -1), sS.coeff(x, 0), sS.coeff(x, 1)) == \
    (1, sp.Rational(1, 2), sp.Rational(1, 8)), sS
P("G2 simple-nu ladder regression: nu = 1/x + 1/2 + x/8 (exact) -> PASS")

nuV = 1 + sp.sqrt(sp.pi/3)/x
sV = sp.series(nuV, x, 0, 3).removeO().expand()
assert sV.coeff(x, 0) == 1 and sV.coeff(x, 1) == 0, sV
P("G3 Verlinde regression: nu = sqrt(pi/3)/x + 1 (c1 = 1, c2 = 0, "
  "exact) -> PASS")
P("")

# ------------------------------------------------------------------ G4
# THE MASTER INVERSION.  AQUAL: mu(z) z = y with z = g/a0, y = g_N/a0,
# nu = z/y.  mu(z) = z + m2 z^2 + m3 z^3.  Set w = sqrt(y); solve
# z = w + b2 w^2 + b3 w^3 order by order.
w, m2, m3, b2, b3 = sp.symbols('w m2 m3 b2 b3')
z = w + b2*w**2 + b3*w**3
eq = sp.expand(z**2 + m2*z**3 + m3*z**4 - w**2)
b2_sol = sp.solve(eq.coeff(w, 3), b2)[0]
b3_sol = sp.solve(eq.coeff(w, 4).subs(b2, b2_sol), b3)[0]
assert b2_sol == -m2/2, b2_sol
assert sp.expand(b3_sol - (5*m2**2/8 - m3/2)) == 0, b3_sol
# nu = z/y = z/w^2 = 1/w + b2 + b3 w  ->  c1 = b2, c2 = b3.
P("G4 MASTER INVERSION (exact): for mu(z) = z + m2 z^2 + m3 z^3 + ...")
P("     c1 = -m2/2          c2 = 5 m2^2/8 - m3/2")
P("   the measured c1 is a measurement of mu''(0) = 2 m2:")
P("   c1 in [0.21, 0.52]  <=>  m2 in [-1.04, -0.42]")
for name, mm2, mm3, exp_c1, exp_c2 in [
        ('simple   mu = z/(1+z)        ', -1, 1,
         sp.Rational(1, 2), sp.Rational(1, 8)),
        ('standard mu = z/sqrt(1+z^2)  ', 0, -sp.Rational(1, 2),
         0, sp.Rational(1, 4)),
        ('exp-mu   mu = 1 - e^(-z)     ', -sp.Rational(1, 2),
         sp.Rational(1, 6), sp.Rational(1, 4), sp.Rational(7, 96)),
        ('BE (mu-side shadow)          ', -1, sp.Rational(13, 12),
         sp.Rational(1, 2), sp.Rational(1, 12))]:
    c1v = b2_sol.subs(m2, mm2)
    c2v = b3_sol.subs([(m2, mm2), (m3, mm3)])
    assert (c1v, c2v) == (exp_c1, exp_c2), (name, c1v, c2v)
    P(f"   {name}: c1 = {c1v}, c2 = {c2v} -> CHECK")
P("")

# ------------------------------------------------------------------ G5
nu1 = sp.Rational(1, 2) + sp.sqrt(x**2 + 4)/(2*x)
nu2 = sp.sqrt(x**2/2 + sp.sqrt(1 + x**4/4))/x
nu3 = (x**3/2 + sp.sqrt(1 + x**6/4))**sp.Rational(1, 3)/x
s1 = sp.series(nu1, x, 0, 3).removeO().expand()
s2 = sp.series(nu2, x, 0, 4).removeO().expand()
s3 = sp.series(nu3, x, 0, 4).removeO().expand()
assert (s1.coeff(x, 0), s1.coeff(x, 1)) == \
    (sp.Rational(1, 2), sp.Rational(1, 8))
assert (s2.coeff(x, 0), s2.coeff(x, 1), s2.coeff(x, 3)) == \
    (0, sp.Rational(1, 4), sp.Rational(1, 32)), s2
assert (s3.coeff(x, 0), s3.coeff(x, 1), s3.coeff(x, 2)) == \
    (0, 0, sp.Rational(1, 6)), s3
P("G5 n-FAMILY (Hees nu_alpha; exact series, cross-checked vs G4):")
P("   n = 1 (simple):   nu = 1/x + 1/2 + x/8            c1 = 1/2")
P("   n = 2 (standard): nu = 1/x + 0   + x/4 + x^3/32   "
  "standard-mu family: c1 = 0 EXACTLY")
P("   n = 3:            nu = 1/x + 0   + 0   + x^2/6    c1 = c2 = 0")
P("   -> the whole n >= 2 family (incl. the field's historical")
P("      'standard' function) carries NO zero-point term.")
P("")

# ------------------------------------------------------------------ G6
nu_s, E = sp.symbols('nu_s E', positive=True)
solA = sp.solve(nu_s*(1 - 1/E) - 1, E)[0]      # e^u from exp-mu
solB = sp.solve(nu_s - 1 - 1/(E - 1), E)[0]    # e^u from boot
assert sp.simplify(solA - solB) == 0, (solA, solB)
nu1_num = sp.nsolve(nu_s*(1 - sp.exp(-nu_s)) - 1, nu_s, 1.3)
P("G6 THE IDENTITY exp-mu == boot: e^u = nu/(nu-1) from BOTH")
P("   mu(z) = 1 - e^(-z) in AQUAL  <=>  nu = 1 + n_BE(nu*y) (the 4F")
P("   quantum-bootstrap bath) — EXACT equivalence, gated.")
P("   fingerprint via G4: c1 = 1/4, c2 = 7/96 = the 4F printed cell;")
P(f"   nu(1) = {float(nu1_num):.4f} (4F boot transition amplitude).")
P("   SCOUT (2026-07-29): the bare exponential mu is apparently NOWHERE")
P("   named or tested in the MOND literature (attribution NOT FOUND;")
P("   the catalog's exponential families live in nu-space with other")
P("   arguments — G9).  The program's boot adjudication — raw dead")
P("   4-7/200 [4F]; hier-M/L flip +75.6 [5C]; binary veto [5F];")
P("   vertical-channel collapse -9 [5M] = 'the 1/4 cell dead")
P("   everywhere' — is then apparently this function's FIRST data")
P("   contest (scout-level), and it is dead on the primary treatment.")
# companion identity (already 4F): simple-mu == classical bath.
yS = sp.symbols('yS', positive=True)
nu_cl = sp.solve(sp.Eq(nu_s, 1 + 1/(nu_s*yS)), nu_s)
nu_cl_pos = [r for r in nu_cl if r.subs(yS, 1) > 0][0]
ref = sp.Rational(1, 2) + sp.sqrt(sp.Rational(1, 4) + 1/yS)
assert sp.simplify(nu_cl_pos - ref) == 0
P("   companion identity (4F, re-gated): mu = z/(1+z) <=> the classical")
P("   bath nu = 1 + 1/(nu*y) = 1/2 + sqrt(1/4 + 1/y) — exact.")
P("   READING: the mu-side's two natural one-parameter-free members are")
P("   the two SELF-CONSISTENT bath closures (classical & quantum); the")
P("   data kill both cells (5C -99, 5M -9) and keep the SOURCE-DRIVEN")
P("   occupation forms — the same exponential in the temperature")
P("   variable x = sqrt(y) (McGaugh's RAR fit = BE, surviving) vs in")
P("   the acceleration variable z = nu*y (exp-mu = boot, dead): the")
P("   data distinguish the exponential's ARGUMENT.")
P("")

# ------------------------------------------------------------------ G7
f1, f2 = sp.symbols('f1 f2')
F = 1 + f1*x**2 + f2*x**4          # analytic in y = x^2, F(0) = 1
nu_add = 1 + F/x
s_add = sp.expand(nu_add - 1/x)
assert s_add.coeff(x, 0) == 1
assert s_add.coeff(x, 2) == 0
assert s_add.coeff(x, 1) == f1 and s_add.coeff(x, 3) == f2
P("G7 ADDITIVE CLASS THEOREM (exact, symbolic F): any")
P("   g_obs = g_N + sqrt(a0 g_N) * F(g_N/a0), F analytic, F(0) = 1")
P("   has nu = 1/x + 1 + F'(0) x + ...:")
P("   additive-analytic class: c1 = 1 EXACTLY (independent of F), and")
P("   the expansion carries ODD rungs only (no x^0-beyond-1, no x^2).")
P("   MEMBERS: Verlinde emergent gravity (F = 1, exact — note-V); the")
P("   superfluid-DM MOND-limit composition a_N + sqrt(a0 a_N)")
P("   (membership FLAGGED: the primary-source composition equation was")
P("   not retrievable by scout — check pending); any dark component")
P("   tracking baryons analytically in g_N.")
P("   A c1 != 1 requires F non-analytic in g_N (half-integer powers)")
P("   — i.e. a structure carrying sqrt(g_N): the temperature variable.")
P("   PRIOR ART (scout, resolves note-V's pending item): Lelli+17")
P("   (arXiv:1702.04355) excluded Verlinde's EG on SPARC at FIT level")
P("   (M/L amplitude tension + the radius-residual signature); the")
P("   COEFFICIENT-level exclusion (c1 = 1 vs measured) is apparently")
P("   new (scout-level).")
P("")

# ------------------------------------------------------------------ G8
P("G8 NON-RAR STRUCTURES (excluded upstream of the ladder):")
P("   - constant-acceleration additions (Rindler/Grumiller, g_obs =")
P("     g_N + a_R): deep limit g_obs -> const, not sqrt(a0 g_N).")
P("   - screened chameleon f(R) (Naik+19, arXiv:1905.13330): no")
P("     universal RAR — screening-radius 'upturn' signature, location")
P("     environment-dependent; their own SPARC bound |f_R0| < 6e-8.")
P("     Structurally not a RAR rival; no fingerprint assigned.")
P("")

# ------------------------------------------------------------------ G9
# THE HEES+16 CATALOG (arXiv:1510.01369).  Forms scout-quoted from the
# paper; ladders exact here.
a_sym = sp.symbols('a_sym')
# nu-tilde_a = (1-e^-y)^(-1/2) + a e^-y,  y = x^2  (symbolic a):
nu_til = (1 - sp.exp(-x**2))**sp.Rational(-1, 2) + a_sym*sp.exp(-x**2)
s_til = sp.series(nu_til, x, 0, 2).removeO().expand()
assert s_til.coeff(x, 0) == a_sym and s_til.coeff(x, 1) == \
    sp.Rational(1, 4), s_til
P("G9 THE HEES+16 CATALOG (exact ladders; forms scout-quoted):")
P("   nu-tilde_a = (1-e^-y)^(-1/2) + a e^-y:")
P("     c1 = a EXACTLY (the family parameter IS the zero-point")
P("     coefficient), c2 = 1/4 -> ladder-viable iff a in [0.21, 0.52]")
# nu-bar_a members:
nu_bar1 = (1 - sp.exp(-x**2))**sp.Rational(-1, 2) + \
    sp.Rational(1, 2)*sp.exp(-x**2)
nu_bar2 = (1 - sp.exp(-x**4))**sp.Rational(-1, 4) + \
    sp.Rational(3, 4)*sp.exp(-x**4)
s_b1 = sp.series(nu_bar1, x, 0, 2).removeO().expand()
s_b2 = sp.series(nu_bar2, x, 0, 4).removeO().expand()
assert (s_b1.coeff(x, 0), s_b1.coeff(x, 1)) == \
    (sp.Rational(1, 2), sp.Rational(1, 4)), s_b1
assert (s_b2.coeff(x, 0), s_b2.coeff(x, 1), s_b2.coeff(x, 3)) == \
    (sp.Rational(3, 4), 0, sp.Rational(1, 8)), s_b2
P("   nu-bar_a = (1-e^(-y^a))^(-1/(2a)) + (1-1/(2a)) e^(-y^a):")
P("     c1 = 1 - 1/(2a) EXACTLY (members gated: a=1 -> 1/2, a=2 -> 3/4)")
P("     -> ladder-viable iff a in [0.63, 1.04]; a >= 2 has c1 >= 3/4")
# nu-hat_a members; a=1 == BE identity:
nu_hat1 = (1 - sp.exp(-x))**sp.Rational(-1, 1)
assert sp.simplify(nu_hat1 - nuBE) == 0
nu_hat2 = (1 - sp.exp(-x**2))**sp.Rational(-1, 2)
nu_hat3 = (1 - sp.exp(-x**3))**sp.Rational(-1, 3)
s_h2 = sp.series(nu_hat2, x, 0, 2).removeO().expand()
s_h3 = sp.series(nu_hat3, x, 0, 3).removeO().expand()
assert (s_h2.coeff(x, 0), s_h2.coeff(x, 1)) == (0, sp.Rational(1, 4))
assert (s_h3.coeff(x, 0), s_h3.coeff(x, 1), s_h3.coeff(x, 2)) == \
    (0, 0, sp.Rational(1, 6)), s_h3
P("   nu-hat_a = (1-e^(-y^(a/2)))^(-1/a):")
P("     a = 1: nu-hat_1 == BE == the RAR-fit EXACTLY (gated identity —")
P("     McGaugh's function is the a=1 member of Hees's hat family);")
P("     a = 2: c1 = 0, c2 = 1/4 (the standard fingerprint); a = 3:")
P("     c1 = c2 = 0 -> only a = 1 carries a zero-point term.")
P("")

# ----------------------------------------------------------------- G10
P("G10 THE PINCER (their Cassini verdicts x our ladder window):")
P("   Hees+16 published verdicts (scout-quoted; primary read queued")
P("   before any paper use): nu-tilde 'completely excluded' by Cassini")
P("   at ALL a; nu_alpha and nu-hat 'marginally acceptable' only at")
P("   a >= 7-8; nu-bar compatible for a >= 2.")
P("   The ladder window [0.21, 0.52] selects: nu-tilde a in")
P("   [0.21, 0.52] (their Cassini-dead); nu-bar a in [0.63, 1.04]")
P("   (their Cassini-dead, needs a >= 2); nu_alpha only n = 1 and")
P("   nu-hat only a = 1 (their Cassini-dead, need a >= 7-8, whose")
P("   members carry c1 = 0 = our dead rung).")
P("   => EMPTY INTERSECTION: no member of the published catalog")
P("   passes BOTH the measured coefficient window and their")
P("   AQUAL-EFE Cassini bound.  This is the amplitude-locked solar")
P("   tension (4K/5S/5I: Q2 ~ 4-6x Cassini for every RAR-compatible")
P("   function) rederived in coefficient language through the field's")
P("   own catalog — function choice cannot escape it; the doors are")
P("   structural (the 7G trajectory formulation passes by 451 orders;")
P("   EFE-screened composition is 6W-excluded on the binaries).")
P("   Caveats: their bounds are their analysis (MG/AQUAL-formulation-")
P("   conditional, fixed-a0), quoted not re-derived; scout-grade until")
P("   the primary table read.")
P("")

# ------------------------------------------------- the verdict table
P("=" * 72)
P("THE FINGERPRINT TABLE (x = sqrt(g_N/a0); measured anchors: c1")
P("profile 0.385-0.519 flat [4S] / 0.208-0.309 hier [4Z], bootstrap")
P("~0.4 +/- 0.3, c1 = 0 excluded Delta(-2lnL) = 56.3 = 7.5 sigma")
P("profile / 95.5% bootstrap; binary c1 ERASED [7J-d], not quoted)")
P("=" * 72)
P("  function/theory        c1       c2     nu(1)   verdict")
P("  BE == RAR-fit == hat_1 1/2      1/12   1.582   reference; in-window")
P("  simple == class. bath  1/2      1/8    1.618   ladder-viable; hier-")
P("                                                 rejected -99 (5C)")
P("  standard-mu (n = 2)    0        1/4    1.272   EXCLUDED at the")
P("                                                 c1 = 0 grade (56.3);")
P("                                                 direct: 4B -56")
P("  n = 3                  0        0      --      excluded (same rung)")
P("  exp-mu == boot         1/4      7/96   1.349   dead on the primary")
P("                                                 treatment (5M);")
P("                                                 hier-flip 5C printed")
P("  Hees nu-tilde_a        a        1/4    --      viable only a in")
P("                                                 [0.21,0.52]; their")
P("                                                 Cassini kills ALL a")
P("  Hees nu-bar_a          1-1/(2a) --     --      viable a in")
P("                                                 [0.63,1.04]; their")
P("                                                 Cassini needs a>=2")
P("  Hees nu-hat_(a>=2)     0        --     --      excluded (c1 = 0)")
P("  Verlinde (note-V)      1        0      --      EXCLUDED (ladder);")
P("                                                 fit-level: Lelli+17")
P("  additive class (G7)    1        free   --      EXCLUDED (ladder;")
P("                                                 extrapolation grade)")
P("  Rindler const-accel    --       --     --      excluded (deep limit)")
P("  chameleon f(R)         --       --     --      structurally non-RAR")
P("                                                 (Naik+19)")
P("")
P("HONESTY BLOCK: the c1 intervals are FAMILY-PROFILE measurements")
P("(the BE-anchored continuous-lambda family, 4S/4Z); transfers to")
P("functions with foreign higher structure (standard's c2 = 1/4,")
P("additive's c2 = free) are curvature-extrapolation grade — quoted as")
P("such.  Where DIRECT contests exist they carry the sharp verdict and")
P("the ladder supplies the WHY: standard-mu (4B: 198-200/200 raw,")
P("Delta(-2lnL) = -56 honest, sign-robust strong lean <-> missing")
P("zero-point digit), boot cell (5C/5M arc <-> c1 = 1/4 with the")
P("self-consistent argument), Verlinde/additive (note-V <-> c1 = 1).")
P("The hier interval CONTAINS 1/4: the exp-mu/boot kill is the")
P("FUNCTION-level vertical-channel result (5M), not a c1-digit result")
P("— the distinction is load-bearing and stated.")
P("")
P("PRIOR-ART STATUS (scout round 2026-07-29): (i) coefficient-level")
P("expansion-and-constraint of the interpolating function: NOT FOUND")
P("(the ladder instrument apparently unpublished, scout-level);")
P("(ii) bare exp-mu: attribution NOT FOUND, apparently untested;")
P("(iii) Verlinde-vs-RAR prior art = Lelli+17 fit-level (named in G7);")
P("(iv) superfluid composition equation: primary source not retrieved")
P("— membership flagged; (v) analytic closed-form LCDM nu(y): NOT")
P("FOUND (Navarro+17 sim-based; Paranjape-Sheth 'analytic' but")
P("equations inaccessible to scout; Desmond semi-empirical, scatter")
P("overpredicted 3.5 sigma) — the LCDM leg needs the DATA-side sim")
P("ladder (8B), no algebra shortcut exists.")

with open('data/stage8a_ladder.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8a_ladder.txt")
