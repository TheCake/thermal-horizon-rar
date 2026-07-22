"""
STAGE 4L: modified-inertia vs modified-gravity on our own data (TODO #18).
Patch-runner over stage3p_v7budget.py (the physical-field v7 fit): replaces
only the engine-invocation block. Modified inertia = the boost is a
per-ORBIT functional, not a per-point field value (Milgrom 2011 gives the
principle, no closed form -- we bracket with two prescriptions):
    mi_a: y_char = g_N(a)/a0                 (semi-major-axis scale)
    mi_t: y_char = g_N(a)/(sqrt(1-e^2) a0)   (exact Kepler time-average <1/r^2>)
A global per-orbit boost B means exact Kepler dynamics with G_eff = B*G:
implemented as engine mode 1 (pure Newton) with M_eff = M*B, while the
vtilde normalization downstream keeps the TRUE mass -- so the entire v7
nuisance machinery (e-mixture, companions, fences, acceptance, noise) is
bit-identical between the MI and MG fits. B interpolates the SAME alpha-
scaled EFE table as MG (conservative EFE bracket; the no-EFE bracket is a
later variant). alpha=0 rows reduce to exact Newton (gate: must reproduce
the MG run's Newton lnL for the same seed).
argv: <presc: mi_a|mi_t> <gtag> <seeds...>   e.g.  mi_t 1p2 31
Appends data/stage4l_summary.txt (same fields as stage3u_summary.txt ->
directly comparable best-lnL against the stored MG rows, same data, same
population seed, same nuisance grids).
"""
import sys

PRESC = sys.argv[1] if len(sys.argv) > 1 else 'mi_t'
assert PRESC in ('mi_a', 'mi_t', 'mi_a_iso', 'mi_t_iso'), PRESC
ISO = PRESC.endswith('_iso')                 # no-EFE bracket (Milgrom's
                                             # frequency argument: internal
                                             # motion decouples from the slow
                                             # external field -> bare nu)
sys.argv = [sys.argv[0]] + sys.argv[2:]      # inner script sees [gtag, seeds]

src = open('calcs/stage3p_v7budget.py', encoding='utf-8-sig').read()

old = """                    e_s = e_of(p, eta, wr)
                    vp = vp_c(p, e_s, tab_a) if al > 0 else None
                    mode = 5 if al > 0 else 1
                    kw = dict(a0=A0_CAN, tab=tab_a, lny0=LNY0, dlny=DLNY,
                              vp=vp) if al > 0 else {}
                    o = run(p['a_s'], e_s, p['psi0'], p['f_ip'], p['M_s'],
                            p['uph'], 8, 2500, mode, **kw)"""
new = """                    e_s = e_of(p, eta, wr)
                    gch = GM*p['M_s']/p['a_s']**2
                    if PRESC.startswith('mi_t'):
                        gch = gch/np.sqrt(np.maximum(1.0-e_s**2, 1e-4))
                    if ISO:
                        yv = np.maximum(gch/A0_CAN, 1e-12)
                        if law == "BE":
                            nub = 1.0/(1.0-np.exp(-np.minimum(np.sqrt(yv),
                                                              40.0)))
                        else:
                            nub = 0.5+np.sqrt(0.25+1.0/yv)
                        Bp = 1.0 + al*(nub-1.0)
                    else:
                        Bp = np.interp(np.log(gch/A0_CAN), LNY_U, tab_a,
                                       right=1.0)
                    o = run(p['a_s'], e_s, p['psi0'], p['f_ip'],
                            p['M_s']*Bp, p['uph'], 8, 2500, 1)"""
assert old in src, "engine block not found -- stage3p changed?"
src = src.replace(old, new, 1)
assert src.count("data/stage3u_summary.txt") >= 1
src = src.replace("data/stage3u_summary.txt", "data/stage4l_summary.txt")
src = src.replace('STAGE 3P g={GTAG}', 'STAGE 4L [' + PRESC + '] g={GTAG}')

ns = {'PRESC': PRESC, 'ISO': ISO}
exec(src, ns)
