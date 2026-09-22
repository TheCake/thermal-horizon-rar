# PRE-REGISTRATION — STAGE 10X: THE METER ROBUSTNESS BATTERY

Registered 2026-09-22, BEFORE any variant fit is run. Parent: 10T
(PREREG-H0METER-10T.md; operative verdict data/stage10t_verdict.txt,
R46-corrected). Purpose: the referee-grade robustness axes never run
on the anchored central — leave-one-out influence, inclination,
quality flag, bulge presence, nuisance-prior width — plus a census
leg for the two flag-pending anchored galaxies. DIAGNOSTIC grade
throughout: no variant replaces the primary; the band 65.4-70.4 and
stat 5.0 are untouched by construction; every credence cell is
pre-signed HOLD on every outcome of this stage.

## 1. Estimator and inheritance (frozen)

calcs/stage10x_meterrobust.py exec-inherits calcs/stage10t_legregrow.py
verbatim to the marker "# ---------------- G10T-1" (data load, CF4,
world builder, nu families, fit_hier/fit_deep, sigma machinery). The
primary world is make_world(lt_slot=None) (leg A = 78, LT out, per
G10T-3b). The operative estimator is fit_deep (deep-converged hier,
use_u=True), warm-started from the reproduced primary optimum. All
variant fits are engine-matched to the primary (same estimator, same
convergence schedule). No flow-leg fit of any kind is run.

## 2. Gates (sky-blind: no variant central is computed in gates mode)

- G10X-1 primary reproduction: the inherited chain (plain -> loose ->
  deep) must reproduce the archived primary a0 = 1.0106e-10 to
  relative 1e-3 and H0_A = 65.36 +- 0.02 (parse
  data/stage10t_skyread.txt; STOP on miss).
- G10X-2 variant-harness identity: the A-var2 world (uma_moved,
  LT-free) must reproduce the archived a0 = 1.0098e-10 to relative
  1e-3 (STOP on miss).
- G10X-3 axis census: per-axis membership counts printed; each axis
  subset must be a strict subset of the 78. An axis with N_sub < 20
  is labeled N-THIN, reported as a co-read, and EXCLUDED from letter
  evaluation (pre-registered clause; deterministic from counts).
- G10X-4 warm-start determinism: refitting the full 78 from the
  reproduced optimum returns the same a0 to relative 1e-6 (the LOO
  loop's identity point; STOP on miss).

## 3. The axes (frozen list; nothing added after this commit)

- A1 LEAVE-ONE-OUT: 78 deep refits, each with one galaxy removed.
  Statistic: max_g |Delta H0(g)|; the full influence list is reported
  (top carriers named).
- A2 INCLINATION LADDER: subsets i >= 40, 45, 50 deg (SPARC master
  inclinations; UMa and reclassified members use their SPARC rows).
- A3 QUALITY: Q = 1 only.
- A4 BULGE-FREE: galaxies whose every surviving point has zero bulge
  contribution (tests the frozen bulge mass-to-light on the meter).
- A5 NUISANCE-PRIOR WIDTH x1.5: the per-galaxy distance+inclination
  prior widths (sv) scaled by 1.5, full membership (tests prior
  sensitivity of the central, not membership).

## 4. Bars (pre-registered form; all computable from counts alone)

Nested subsets wander from the full-sample central by galaxy scatter
alone, so flat bars would fire spuriously. Each axis therefore gets

  bar_axis = max( 2 x sigma_exp , 1.0 km/s/Mpc ),
  sigma_exp = sigma_stat x sqrt(N_full/N_sub - 1),

with sigma_stat = 4.99 (the archived primary-engine bootstrap) and
N_full = 78. For A5 (same membership) sigma_exp = 0, so the floor
1.0 applies. For A1 the same formula with N_sub = 77 applies to the
MAXIMUM over the 78 leave-one-out deltas (a conservative carrier
bar: the max of 78 draws is compared against a 2-sigma single-draw
band, biasing the axis TOWARD firing; accepted, since firing is the
informative direction for an influence axis).

## 5. Letter grammar (trap-#31 form: exclusive, exhaustive, in-grid
falsifiable; tie-break unnecessary by construction)

Let E = the set of letter-eligible axes (all of A1-A5 minus any
N-THIN exclusions from G10X-3).

- X-METER-SENSITIVE(list): fires iff at least one axis in E exceeds
  its bar; the letter names every firing axis and its delta. This
  branch is FALSE when no axis fires.
- X-METER-ROBUST: fires iff no axis in E exceeds its bar. This
  branch is FALSE when any axis fires.

The two branches partition the outcome space (any-vs-none over the
same frozen set E); both are reachable in-grid (subsets of N ~ 25-50
have expected null wander of several km/s/Mpc, so bars CAN be
exceeded; the 10T membership variants moved <= 0.6, so bars CAN be
passed). Grade cap (trap-#29 wire, in code): both letters are
DIAGNOSTIC — the stage prints, in the same block as the letter, that
no variant reading replaces the primary and that the band and stat
are uncut-primary properties. A SENSITIVE letter triggers a review
round before anything is quoted from it; a ROBUST letter is quotable
as robustness rows in Proposal F section 5/6e (co-read grade).

## 6. The census leg (informational; NO letter clause — external
resolvability cannot be wired into letter grammar)

- C1 D564-8 and D631-7 (the two flag-pending anchored members): PGC
  resolution (cache first, then external lookup), then the CF4
  cross-check. CLEARING RULE (frozen): the flag is cleared iff a CF4
  per-method anchor-grade modulus exists (mas/ceph/trgb/sbf/snII in
  the stage's PRIO order) AND |D_CF4/D_SPARC - 1| <= max(2 eD/D,
  0.15). Any other outcome: the flag STANDS with the measured gap
  reported. NO membership or distance change in either case
  (membership moves only through a registered census stage).
- C2 KK98-251 (flow; census hygiene only): one more PGC resolution
  attempt by external identifier lookup. Outcome recorded either
  way; no fit touches it.

## 7. Outputs and discipline

Gates mode -> data/stage10x_gates.txt (no variant centrals). Run
mode -> data/stage10x_battery.txt with the letter, the influence
list, every axis table, and the census-leg record. Prereg committed
before the gates run; gates committed before the battery run; any
amendment is a dated pre-quote entry here with the prior run
preserved. Credences: HOLD 53/8 on every cell (pre-signed).
Wall-clock estimate: ~90 deep fits, minutes on the certified pool or
serial.

## Amendments

- A1 (2026-09-22, BEFORE any run; self-caught at design review): the
  G10X-4 identity bar as first registered (relative 1e-6 on a0) is
  tighter than the estimator's own convergence floor — the inner
  optimizer runs at fatol 1e-7 on the objective, which with the
  measured curvature (sigma_la0 ~ 0.033 dex) implies ~1e-5 dex
  path jitter, i.e. ~2e-5 relative on a0, and the alternation stops
  at tol 5e-4 objective units. A perfect replication could therefore
  fail the 1e-6 bar (the 10W-A2 failure class: a bar a correct
  replication cannot meet). G10X-4 bar re-set to relative 1e-4 on
  a0, matching the verification-grade reproduction precision the
  10T/R46 record actually achieved. No run has been executed under
  either bar.

- A2 (2026-09-22, BEFORE any run): the census leg's external PGC
  resolution is performed OUTSIDE the stage (session-documented
  lookups written to data/stage10x_pgc.json with provenance); the
  stage script itself is deterministic and offline — it reads the
  file if present and reports UNRESOLVED for any name it does not
  contain. The clearing rule of section 6 is unchanged.
