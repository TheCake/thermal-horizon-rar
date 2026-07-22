"""
STAGE 4M (#2e-b): close the 4J residual cells in the model. Two ingredients
the 4J autopsy identified as missing from v7, patched into stage3p (MG
engine path untouched -> directly comparable to the stored 3U baselines,
same seed/population/grids):
  - the closest-approach flyby arm: v7's FLY template is gamma in [0,30]
    only (asymptote regime); 4J showed the data's unbound/peri structure at
    gamma ~ 90 in the vt ceiling band. Variants: 'fly90' (pure 90-arm,
    vt in [1.2, 2.2] x gamma-weight max(gcen-60,0)), 'flymix' (50/50
    asymptote + 90-arm).
  - the radial e-ceiling: e_rad in [0.9, 0.995] -> [0.9, 0.9995] ('erad'),
    letting the apo face reach the implied e ~ 0.99 of the 4J slow island.
  - 'both' = flymix + erad.
Readout: dlnL(variant - stored 3U baseline) per law + alpha-hat stability
(residual explained AND alpha invariant = the 4J cells were a model-shape
gap, not an alpha bias). argv: <variant> <gtag> <seeds...>
Appends data/stage4m_summary.txt.
"""
import sys

VAR = sys.argv[1] if len(sys.argv) > 1 else 'both'
assert VAR in ('fly90', 'flymix', 'erad', 'both'), VAR
sys.argv = [sys.argv[0]] + sys.argv[2:]

src = open('calcs/stage3p_v7budget.py', encoding='utf-8-sig').read()

old_fly = """FLY = np.zeros((NV, NG))
for i in range(NV):
    for j in range(NG):
        FLY[i,j] = (1.0 if (0.5 <= vcen[i] <= 3.0) else 0.0) \\
                   * max(30.0-gcen[j], 0.0)
FLY /= max(FLY.sum(), 1e-12)"""
new_fly = """FLYA = np.zeros((NV, NG)); FLY90 = np.zeros((NV, NG))
for i in range(NV):
    for j in range(NG):
        FLYA[i,j] = (1.0 if (0.5 <= vcen[i] <= 3.0) else 0.0) \\
                    * max(30.0-gcen[j], 0.0)
        FLY90[i,j] = (1.0 if (1.2 <= vcen[i] <= 2.2) else 0.0) \\
                     * max(gcen[j]-60.0, 0.0)
FLYA /= max(FLYA.sum(), 1e-12); FLY90 /= max(FLY90.sum(), 1e-12)
FLY = FLY90 if FVAR == 'fly90' else \\
      (0.5*FLYA + 0.5*FLY90 if FVAR in ('flymix', 'both') else FLYA)
FLY /= max(FLY.sum(), 1e-12)"""
assert old_fly in src, "FLY block not found"
src = src.replace(old_fly, new_fly, 1)

old_er = "    e_rad = 0.9+0.095*p['u_e']"
new_er = ("    e_rad = 0.9+(0.0995 if FVAR in ('erad', 'both') else 0.095)"
          "*p['u_e']")
assert old_er in src, "e_rad line not found"
src = src.replace(old_er, new_er, 1)

src = src.replace("data/stage3u_summary.txt", "data/stage4m_summary.txt")
src = src.replace('STAGE 3P g={GTAG}', 'STAGE 4M [' + VAR + '] g={GTAG}')

ns = {'FVAR': VAR}
exec(src, ns)
