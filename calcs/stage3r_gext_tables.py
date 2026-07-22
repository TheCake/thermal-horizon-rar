"""
STAGE 3R: EFE boost tables for the g_ext scan (TODO #2f). Reuses the validated
stage-2G solver (gates G1/G2 passed at 0.01%) to produce boost tables at
g_ext/a0 in {1.4, 1.6, 2.2, 2.4} for both nu-families (1.9 already exists as
data/efe_boost_*.npy). Writes data/efe_boost_{law}_g{val}.npy.
"""
import numpy as np

src = open('calcs/qumond_efe_solver.py').read()
ns = {}
exec(src.split("y = 1.0/r**2")[0], ns)   # solver machinery only, no production
solve, nu_simple, nu_be = ns['solve'], ns['nu_simple'], ns['nu_be']
r = ns['r']
y = 1.0/r**2

for eN in (1.4, 1.6, 2.2, 2.4):
    for name, nu in (("simple", nu_simple), ("be", nu_be)):
        b = solve(nu, eN)
        tag = str(eN).replace('.', 'p')
        np.save(f'data/efe_boost_{name}_g{tag}.npy', np.stack([y, b]))
        i = np.argmin(np.abs(y-1.0))
        print(f"eN={eN} {name}: boost(y=1) = {b[i]:.3f}  saved")
print("done")
