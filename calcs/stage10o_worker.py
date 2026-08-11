"""Stage 10O worker module (amendment A3): one-orbit secular integration.

Lives in its own module so ProcessPoolExecutor workers import ONLY this
(never the stage script body). The Hamiltonian is M&H 2023 eq. 30
verbatim (same transcription as the stage; the round-34 GA half carries
an independent transcription).
"""

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

_j, _jz, _w, _G = sp.symbols('j j_z omega Gamma')
_H = (1 / _j**2) * ((_j**2 - 3 * _G * _jz**2) * (5 - 3 * _j**2)
                    - 15 * _G * (_j**2 - _jz**2) * (1 - _j**2)
                    * sp.cos(2 * _w))
H_fn = sp.lambdify((_j, _w, _jz, _G), _H, 'numpy')
dHdj_fn = sp.lambdify((_j, _w, _jz, _G), sp.diff(_H, _j), 'numpy')
dHdw_fn = sp.lambdify((_j, _w, _jz, _G), sp.diff(_H, _w), 'numpy')


def mix_orbit(task):
    """Integrate one orbit; return its two half-window e-sample arrays,
    end state, and relative H drift."""
    j0, w0, jz, gam, tau_max, t_frac, n_eval = task
    jlo = abs(jz) + 1e-9

    def f(t, y):
        jj = min(max(y[0], jlo), 1 - 1e-12)
        return (-dHdw_fn(jj, y[1], jz, gam), dHdj_fn(jj, y[1], jz, gam))

    t_ev = np.linspace(t_frac * tau_max, tau_max, n_eval)
    sol = solve_ivp(f, (0.0, tau_max), [j0, w0], method='DOP853',
                    rtol=1e-8, atol=1e-10, t_eval=t_ev)
    jj = np.clip(sol.y[0], jlo, 1 - 1e-12)
    e_s = np.sqrt(np.clip(1 - jj**2, 0, 1))
    half = n_eval // 2
    H0 = H_fn(j0, w0, jz, gam)
    H1 = H_fn(jj[-1], sol.y[1, -1], jz, gam)
    drift = abs(H1 - H0) / (abs(H0) + 1.0)
    return e_s[:half], e_s[half:], jj[-1], sol.y[1, -1] % (2 * np.pi), drift


if __name__ == '__main__':
    import sys

    z = np.load(sys.argv[1])
    j0a, w0a, jza = z['j0'], z['w0'], z['jz']
    gam = float(z['gam'])
    tau = float(z['tau'])
    frac = float(z['frac'])
    neval = int(z['neval'])
    n = len(j0a)
    half = neval // 2
    e1 = np.empty((n, half))
    e2 = np.empty((n, neval - half))
    jend = np.empty(n)
    wend = np.empty(n)
    drift = np.empty(n)
    for m in range(n):
        a, b, je, we, dr = mix_orbit((j0a[m], w0a[m], jza[m], gam,
                                      tau, frac, neval))
        e1[m], e2[m], jend[m], wend[m], drift[m] = a, b, je, we, dr
    np.savez(sys.argv[2], e1=e1, e2=e2, jend=jend, wend=wend, drift=drift)
