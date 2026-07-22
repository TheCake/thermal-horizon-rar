"""
STAGE 4N (#6): the reconciliation experiments -- where is the field's
10-sigma-Newton vs 16-sigma-MOND disagreement manufactured? Ablations of
OUR v7 fit toward the published analyses' choices (patch-runner over
stage3p, MG engine path untouched, stored 3U baselines comparable):

  'vtonly'    : direction channel REMOVED (gamma collapsed to one bin) --
                every published analysis is vtilde-only; how much of our
                Newton rejection and alpha localization did the gamma
                dimension buy?
  'freecomp'  : the companion fraction UNFENCED (grid 0..0.8 instead of the
                photometry-bounded 0/0.1) -- the Banik-style freedom: an
                unbounded close-binary fraction can absorb velocity excess.
  'banikproxy': vtonly + freecomp + their separation window (2-30 kAU,
                DROPPING our 0.2-2 kAU Newtonian anchor bin) -- the closest
                one-patch proxy to the published Newton-favored setups.

Readout per variant vs stored baseline (+108.7/+98.8, seed 31): dlnL(Newton)
and alpha-hat. Labeled honestly: these are BANIK-STYLE ablations of our
pipeline, not a line-by-line reproduction of their code.
argv: <variant> <gtag> <seeds...>. Appends data/stage4n_summary.txt.
"""
import sys

VAR = sys.argv[1] if len(sys.argv) > 1 else 'vtonly'
assert VAR in ('vtonly', 'freecomp', 'banikproxy'), VAR
sys.argv = [sys.argv[0]] + sys.argv[2:]

src = open('calcs/stage3p_v7budget.py', encoding='utf-8-sig').read()

if VAR in ('vtonly', 'banikproxy'):
    old = "GE = np.linspace(0, 90, 7)\nNV, NG = 20, 6"
    new = "GE = np.linspace(0, 90, 2)\nNV, NG = 20, 1"
    assert old in src, "GE/NV block not found"
    src = src.replace(old, new, 1)
    old = """        FLY[i,j] = (1.0 if (0.5 <= vcen[i] <= 3.0) else 0.0) \\
                   * max(30.0-gcen[j], 0.0)"""
    new = """        FLY[i,j] = (1.0 if (0.5 <= vcen[i] <= 3.0) else 0.0) \\
                   * (max(30.0-gcen[j], 0.0) if NG > 1 else 1.0)"""
    assert old in src, "FLY factor not found"
    src = src.replace(old, new, 1)

if VAR in ('freecomp', 'banikproxy'):
    old = "FCOMP_GRID = np.array([0.0, 0.10])"
    new = "FCOMP_GRID = np.array([0.0, 0.2, 0.4, 0.6, 0.8])"
    assert old in src, "FCOMP grid not found"
    src = src.replace(old, new, 1)

if VAR == 'banikproxy':
    old = ("SBINS = [(0.2,2),(2,6),(6,20),(20,50)]\n"
           "SC2 = (np.array([0.63, 3.46, 11.0, 31.6])/31.6)**2")
    new = ("SBINS = [(2,6),(6,20),(20,30)]\n"
           "SC2 = (np.array([3.46, 11.0, 24.5])/31.6)**2")
    assert old in src, "SBINS/SC2 block not found"
    src = src.replace(old, new, 1)

src = src.replace("data/stage3u_summary.txt", "data/stage4n_summary.txt")
src = src.replace('STAGE 3P g={GTAG}', 'STAGE 4N [' + VAR + '] g={GTAG}')

ns = {}
exec(src, ns)
