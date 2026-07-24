"""
STAGE 6O (frozen-bath Test 1): the galaxy-coherence contest.
PRE-REGISTRATION COMMITTED BEFORE EXECUTION.

THE PREDICTION UNDER TEST: if the horizon bath decorrelates on ~1/H
(Hubble time), the modes are FROZEN on orbital timescales and a mode
wavelength c^2/sqrt(g a0) >> galaxy size makes the draw coherent
across a galaxy: each galaxy g carries ONE standard-normal draw z_g
imprinting delta log g_obs(point) = z_g * t(x)/( sqrt(N)*ln10 ),
t(x) = e^{-x/2} = sqrt(n/(n+1)) — the thermal mode-variance shape 4T
fit at point level, now as a per-galaxy COHERENT mean channel. The
astrometric channels (distance: g_bar-invariant, uniform in
log g_obs; inclination: uniform) are x-FLAT — so the discriminating
leverage is the within-galaxy radial SHAPE, which the 4W generic
vertical channel never used. This is the identifiability unlock 4W
lacked, run as a shape-template contest.

MODEL: 5M/6F vertical machinery + one new per-galaxy nuisance:
r = lgobs - log10(gN*nu_BE(gN/a0)) - dv_g - A*z_g*t(x), priors
dv_g ~ N(0, sigma_v_g) (measured), z_g ~ N(0,1), A >= 0 global
(A = 1/(sqrt(N)*ln10); N = mode count). Inner (dv_g, z_g) exact 2x2
linear solves; dml as before; lensing leg carries NO template
(stacked ensembles average the draw; disclosed).

ESTIMATOR CORRECTION (post-commit, pre-results; disclosed): the
committed MAP treatment of z_g was monotone-degenerate (a new
per-galaxy parameter always lowers a penalized objective; the first
launch crashed on its own G0 gate before any result was read). The
z-channel is now EXACTLY MARGINALIZED: the per-galaxy 2x2 Gaussian
integral over (dv_g, z_g) adds the Occam term
ln[det M_g] - ln[M11_g], M_g = [[Sw+1/sv^2, S_At],[S_At, S_A2t2+1]],
which is IDENTICALLY ZERO at A = 0 — the historic baseline objective
and the G0 regression are unchanged. Bands and outcome tree unchanged.

GATES: G0 A=0 regression = the 5P/6I BE vertical fit (-12152.49,
|d| < 2 — same model). G-INJ-R (recovery): inject A_inj = 0.10 with
fresh draws into the real data -> profiled A_hat in [0.06, 0.14]
(the 4W failure mode would absorb it into dv and return ~0). G-INJ-N
(null): same-amplitude test with the template SHUFFLED within each
galaxy (x-shape broken, amplitude structure kept) -> Delta < 4.

PRE-REGISTERED OUTCOMES:
 DETECT: Delta(-2lnL) <= -9 at A_hat with N_hat in [5, 500] and both
   injection gates pass -> first evidence of a galaxy-coherent
   thermal channel; the frozen reading gains its predicted
   signature. Credence (bath-microphysics conditional, now ~8%):
   pre-committed direction UP; magnitude 12-20% if N_hat is
   compatible with 4T's point-level 20-60, less if not.
 BOUND: Delta > -9 -> quote the A where Delta crosses +4 as the
   ~95% bound; N_coh >= N_95. Constrains the COHERENT limit of the
   frozen reading only; the frequency-local limit is 4T's channel,
   already characterized. No credence move.
 AMBIG: either injection gate fails -> the identifiability boundary
   survives shape leverage; reported as the honest 4W echo.
Writes data/stage6o_cohtest.txt.
"""
import glob, math, os, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
DELTA_PRIOR = 0.2
LENS_CUT_FID = -14.25
LN10 = math.log(10)
S_ML = 0.1*LN10

meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta[t[0]] = (float(t[5]), int(t[17]), float(t[2]), float(t[3]),
                      float(t[6]))
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
sigv_g_map = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 10.0, 1.0, 3.0))
    if inc < 30 or q > 2: continue
    kept += 1
    irad = math.radians(inc)
    sv = math.hypot((eD/max(D, 1e-3))/LN10,
                    2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10)
    sigv_g_map[gi] = max(sv, 0.01)
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC; gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        g_gas.append(gg); g_dsk.append(gd); g_bul.append(gb)
        gobs.append(Vo*Vo/R*KPC)
        sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs0 = np.log10(gobs)
ug = np.unique(gal_id)
NGal = len(ug)
gmap = {g: i for i, g in enumerate(ug)}
gidx = np.array([gmap[g] for g in gal_id])
SIGV = np.array([sigv_g_map[g] for g in ug])
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]

ML = np.loadtxt('data/lensing_rar/mistele2024_table1.txt')
l_gbar, l_gobs, l_stat, l_syst = ML.T
l_sig2 = l_stat**2 + l_syst**2
lmask = l_gbar >= LENS_CUT_FID

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))

def tshape(gN, a0):
    x = np.sqrt(np.clip(gN/a0, 1e-14, None))
    return np.exp(-0.5*x)

def m2(th, dml, dv, zz, A, lgobs, tperm=None):
    la0, f, s_int, dlt = th
    if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4) or abs(dlt) > 0.8: return 1e12
    a0 = 10**la0
    fac = f*np.exp(dml[gidx])
    gN = g_gas + fac*g_dsk + g_bul
    tt = tshape(gN, a0) if tperm is None else tperm
    r = (lgobs - np.log10(gN*nu_be(gN/a0)) - dv[gidx]
         - A*zz[gidx]*tt)
    se2 = sig2 + s_int*s_int
    out = np.sum(r*r/se2 + np.log(se2))
    # exact z-marginalization: the per-galaxy Occam determinant
    # (identically zero at A = 0 -> baseline objective preserved)
    w = 1.0/se2
    At = A*tt
    Sw_g = np.bincount(gidx, w, minlength=NGal)
    St_g = np.bincount(gidx, w*At, minlength=NGal)
    Stt_g = np.bincount(gidx, w*At*At, minlength=NGal)
    M11 = Sw_g + 1.0/(SIGV*SIGV)
    out += np.sum(np.log(np.maximum(M11*(Stt_g + 1.0) - St_g*St_g,
                                    1e-300)) - np.log(M11))
    lg = l_gbar[lmask] + dlt
    rl = l_gobs[lmask] - (lg + np.log10(nu_be(10**lg/a0)))
    out += np.sum(rl*rl/l_sig2[lmask] + np.log(l_sig2[lmask]))
    out += (dlt/DELTA_PRIOR)**2
    out += np.sum(dml*dml)/(S_ML*S_ML)
    out += np.sum(dv*dv/(SIGV*SIGV)) + np.sum(zz*zz)
    return out

def fit_t(A, lgobs, th0=None, dml0=None, dv0=None, zz0=None,
          tol=0.05, max_rounds=15, tperm=None):
    dml = np.zeros(NGal) if dml0 is None else dml0.copy()
    dv = np.zeros(NGal) if dv0 is None else dv0.copy()
    zz = np.zeros(NGal) if zz0 is None else zz0.copy()
    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08, 0.0]]
                  if rd == 0 and th0 is None else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2(t, dml, dv, zz, A, lgobs, tperm), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int, dlt = best.x
        a0 = 10**la0
        se2c = s_int*s_int
        for _ in range(3):
            fac = f*np.exp(dml[gidx])
            gN = g_gas + fac*g_dsk + g_bul
            tt = tshape(gN, a0) if tperm is None else tperm
            r0 = lgobs - np.log10(gN*nu_be(gN/a0))
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                w = 1.0/(sig2[mm] + se2c)
                Sw = np.sum(w)
                St = np.sum(w*A*tt[mm])
                Stt = np.sum(w*(A*tt[mm])**2)
                Sr = np.sum(w*r0[mm])
                Str = np.sum(w*A*tt[mm]*r0[mm])
                M11 = Sw + 1.0/SIGV[gi2]**2
                det = M11*(Stt + 1.0) - St*St
                if abs(det) < 1e-30:
                    dv[gi2], zz[gi2] = Sr/M11, 0.0
                else:
                    dv[gi2] = ((Stt + 1.0)*Sr - St*Str)/det
                    zz[gi2] = (M11*Str - St*Sr)/det
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = g_gas[mm] + fc*g_dsk[mm] + g_bul[mm]
                    tt2 = (tshape(gN2, a0) if tperm is None
                           else tperm[mm])
                    rr = (lgobs[mm] - np.log10(gN2*nu_be(gN2/a0))
                          - dv[gi2] - A*zz[gi2]*tt2)
                    s2 = sig2[mm] + se2c
                    return np.sum(rr*rr/s2) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2(best.x, dml, dv, zz, A, lgobs, tperm)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2(t, dml, dv, zz, A, lgobs, tperm),
                 list(best.x), method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best, dml, dv, zz

def profile(lgobs, agrid, tag, tperm=None):
    out = []
    th, dm, dvv, zz = None, None, None, None
    for A in agrid:
        t0 = time.time()
        bb, dm, dvv, zz = fit_t(A, lgobs, th0=th, dml0=dm, dv0=dvv,
                                zz0=zz, tperm=tperm)
        th = bb.x
        out.append((A, bb.fun))
        Nval = (1.0/(A*LN10))**2 if A > 0 else float('inf')
        line = (f"  {tag} A={A:.3f} (N={Nval:7.1f}): {bb.fun:10.2f}  "
                f"({(time.time()-t0)/60:.1f} min)")
        print(line, flush=True)
        L.append(line)
        save()
    return out

REF_BE_V = -12152.49
L = [f"STAGE 6O galaxy-coherence contest: {kept} galaxies, "
     f"{len(gobs)} points + {int(lmask.sum())} lensing; template "
     f"t(x) = e^(-x/2), one draw per galaxy; A = 1/(sqrt(N) ln10)", ""]
def save():
    with open('data/stage6o_cohtest.txt', 'w') as f:
        f.write("\n".join(L) + "\n")
save()

AGRID = [0.0, 0.02, 0.035, 0.05, 0.065, 0.08, 0.10, 0.125, 0.15]
prof = profile(lgobs0, AGRID, "real")
base = prof[0][1]
d0 = base - REF_BE_V
L.append(f"G0 A=0 regression: {base:.2f} (5P {REF_BE_V:.2f}, d={d0:+.2f})"
         f" -> {'PASS' if abs(d0) < 2.0 else 'FAIL'}")
assert abs(d0) < 2.0
deltas = [(A, v - base) for A, v in prof]
Abest, dbest = min(deltas, key=lambda t: t[1])
L.append(f"profile: best A = {Abest:.3f} "
         f"(N = {(1.0/(Abest*LN10))**2 if Abest > 0 else float('inf'):.1f})"
         f", Delta = {dbest:+.2f}")
# 95% bound: first A where Delta > +4
Abound = None
for A, d in deltas:
    if A > 0 and d > 4.0:
        Abound = A
        break
if Abound:
    L.append(f"bound: Delta > +4 at A = {Abound:.3f} -> N_coh >= "
             f"{(1.0/(Abound*LN10))**2:.1f}")
save()

# ---------------- injection gates ----------------
rng = np.random.default_rng(101)
a0f = 10**prof[0][1] if False else A0_FID
gNf = g_gas + 1.0*g_dsk + g_bul
tt_f = tshape(gNf, a0f)
AINJ = 0.10
zstar = rng.standard_normal(NGal)
lg_inj = lgobs0 + AINJ*zstar[gidx]*tt_f
COARSE = [0.0, 0.05, 0.08, 0.10, 0.13, 0.16]
L.append("")
L.append(f"G-INJ-R: A_inj = {AINJ} with fresh draws:")
save()
profR = profile(lg_inj, COARSE, "injR")
baseR = profR[0][1]
Ah = min(((A, v - baseR) for A, v in profR), key=lambda t: t[1])[0]
okR = 0.06 <= Ah <= 0.14
L.append(f"G-INJ-R: recovered A_hat = {Ah:.3f} -> "
         f"{'PASS' if okR else 'FAIL (absorption = 4W echo)'}")
save()
# null: shuffled template within each galaxy
tperm = tt_f.copy()
for gi2 in range(NGal):
    mm = GIDXS[gi2]
    tperm[mm] = rng.permutation(tperm[mm])
L.append("")
L.append("G-INJ-N: real data, template x-shape shuffled within galaxy:")
save()
profN = profile(lgobs0, [0.0, 0.08, 0.15], "injN", tperm=tperm)
baseN = profN[0][1]
dN = min(v - baseN for _, v in profN)
okN = dN > -4.0
L.append(f"G-INJ-N: best Delta with shuffled template = {dN:+.2f} -> "
         f"{'PASS' if okN else 'FAIL'}")
save()

# ---------------- verdict ----------------
Nh = (1.0/(Abest*LN10))**2 if Abest > 0 else float('inf')
if dbest <= -9.0 and 5 <= Nh <= 500 and okR and okN:
    call = f"DETECT: coherent channel at N_hat = {Nh:.0f}"
elif not (okR and okN):
    call = "AMBIG (injection gate failure -- the 4W identifiability echo)"
else:
    call = (f"BOUND: no coherent detection (best Delta {dbest:+.1f}); "
            + (f"N_coh >= {(1.0/(Abound*LN10))**2:.0f}" if Abound
               else "no crossing on grid"))
L.append("")
L.append(f"PRE-REGISTERED VERDICT: {call}")
save()
print("\n".join(L[-6:]))
print("\nSTAGE 6O done -> data/stage6o_cohtest.txt")
