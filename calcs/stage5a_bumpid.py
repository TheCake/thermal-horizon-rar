"""
STAGE 5A (O4a): the transition bump's identity -- a per-galaxy template
contest. Two closed-form templates, both linear channels in residual space:

  t_p(x) = d log10 nu_p / dp at p = 1/2 (the section-3 screening family) --
           PEAKED AT THE TRANSITION; coefficient = per-galaxy screening-index
           scatter delta_p (prior 0.1).
  t_e(x) = d log10 g_model / d e_N at e ~ 0 (Chae-Milgrom external-field
           formula, simple family) -- DEEP-WEIGHTED; coefficient = per-galaxy
           environmental field delta_e (prior 0.03 a0; one-sided physical but
           fit symmetric, disclosed).

Machinery: per-galaxy VERTICAL offsets (measured priors, as 4W) plus the
template channels; all per-galaxy coefficients solved EXACTLY per galaxy
(linear ridge solve, basis [1, t_p, t_e]); free 6-bin scatter (frozen
full-sample edges) so the bump bin is read directly.

Variants: V0 = vertical only | VP = +t_p | VE = +t_e | VPE = both.
Verdict key: the bump (bin4) drops to neighbor level under VP => sharpness
scatter; under VE => environmental; under neither => structure unexplained.
Deep bins under VE monitor the EFE-vs-thermal rivalry for the monotone term.

Gates: G1 injection of synthetic p-scatter (0.1) produces a bump in V0 and
VP removes it (bin4 drop > 50%); G2 template orthogonality check printed
(corr(t_p, t_e) over points; high overlap would void the contest).
Writes data/stage5a_bumpid.txt.
"""
import glob, math, os
import numpy as np
from scipy.optimize import minimize

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)
P_PRIOR = 0.1
E_PRIOR = 0.03

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

g_gas, g_dsk, g_bul, gobs, sig, gal_id, svl = [], [], [], [], [], [], []
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 10.0, 1.0, 3.0))
    if inc < 30 or q > 2: continue
    kept += 1
    irad = math.radians(inc)
    sv = max(math.hypot((eD/max(D, 1e-3))/LN10,
                        2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10),
             0.01)
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
        gal_id.append(gi); svl.append(sv)
g_gas, g_dsk, g_bul, gobs, sig, gal_id, svl = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id, svl))
sig2 = sig*sig
lgobs = np.log10(gobs)
ug = np.unique(gal_id)
NGal = len(ug)
gmap = {g: i for i, g in enumerate(ug)}
gidx = np.array([gmap[g] for g in gal_id])
GIDXS = [np.where(gidx == i)[0] for i in range(NGal)]
SIGV = np.array([svl[GIDXS[i][0]] for i in range(NGal)])

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_p(y, p):
    yc = np.clip(y, 1e-14, None)
    return (1.0-np.exp(-yc**p))**(-1.0/(2.0*p))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def cm_gmodel(y, e):
    """Chae-Milgrom EFE formula (simple family): returns log10 of the boost."""
    yq = np.clip(y, 1e-14, None)
    if e <= 0: return np.log10(nu_simple(yq))
    be = 1.1*e
    yb = np.sqrt(yq*yq + be*be)
    sq = np.sqrt(0.25+1.0/yb); nus = 0.5+sq
    nuhat = (1.0/yb)/(2.0*nus*sq)
    return np.log10(nus*(1.0+np.tanh((be/yq)**1.2)*nuhat/3.0))

x_fid = np.sqrt((g_gas + g_dsk + g_bul)/A0_FID)
QE = np.quantile(x_fid, np.linspace(0, 1, 7)); QE[0], QE[-1] = 0.0, np.inf
BIN_FID = np.clip(np.searchsorted(QE, x_fid, side='right')-1, 0, 5)

y_fid = x_fid*x_fid
T_P = (np.log10(nu_p(y_fid, 0.55)) - np.log10(nu_p(y_fid, 0.45)))/0.1
T_E = (cm_gmodel(y_fid, 0.05) - cm_gmodel(y_fid, 0.0))/0.05
cc_pe = np.corrcoef(T_P, T_E)[0, 1]

def fit_var(channels, rounds=3, lg_override=None):
    """channels subset of {'p','e'}; vertical always on. Free 6-bin scatter."""
    global lgobs
    lg_save = lgobs
    if lg_override is not None: lgobs = lg_override
    coef = np.zeros((NGal, 3))          # [dv, dp, de]
    use = [True, 'p' in channels, 'e' in channels]
    pri = np.array([0.0, 1.0/P_PRIOR**2, 1.0/E_PRIOR**2])
    gb = None
    try:
        for rd in range(rounds):
            def gobj(th):
                la0, f = th[0], th[1]
                if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
                sb = np.asarray(th[2:8])
                if np.any(sb < 1e-4) or np.any(sb > 0.5): return 1e12
                gN = g_gas + f*g_dsk + g_bul
                x = np.sqrt(gN/10**la0)
                r = lgobs - np.log10(gN*nu_be(gN/10**la0))
                r = r - coef[gidx, 0] - coef[gidx, 1]*T_P - coef[gidx, 2]*T_E
                s = sb[BIN_FID]
                se2 = sig2 + s*s
                out = np.sum(r*r/se2 + np.log(se2))
                out += np.sum(coef[:, 0]**2/np.maximum(SIGV, 1e-3)**2)
                out += pri[1]*np.sum(coef[:, 1]**2) + pri[2]*np.sum(coef[:, 2]**2)
                return out
            starts = ([list(gb.x)] if gb is not None else []) + \
                     [[math.log10(A0_FID), 1.0] + [0.08]*6]
            gbest = None
            for t0 in starts:
                b = minimize(gobj, t0, method='Nelder-Mead',
                             options=dict(maxiter=8000, xatol=1e-6, fatol=1e-6))
                if gbest is None or b.fun < gbest.fun: gbest = b
            gb = gbest
            la0, f, sb = gb.x[0], gb.x[1], np.asarray(gb.x[2:8])
            gN = g_gas + f*g_dsk + g_bul
            r0 = lgobs - np.log10(gN*nu_be(gN/10**la0))
            se2 = sig2 + sb[BIN_FID]**2
            U = np.stack([np.ones(len(r0)), T_P, T_E], axis=1)
            for gi2 in range(NGal):
                mm = GIDXS[gi2]
                w = 1.0/se2[mm]
                act = [0] + ([1] if use[1] else []) + ([2] if use[2] else [])
                Ua = U[mm][:, act]
                A = (Ua*w[:, None]).T @ Ua
                pr = np.array([1.0/max(SIGV[gi2], 1e-3)**2, pri[1], pri[2]])
                A += np.diag(pr[act])
                bvec = (Ua*w[:, None]).T @ r0[mm]
                sol = np.linalg.solve(A, bvec)
                coef[gi2, :] = 0.0
                for kk, a in enumerate(act):
                    coef[gi2, a] = sol[kk]
    finally:
        lgobs = lg_save
    return gb, coef

L = [f"STAGE 5A bump identity: {kept} galaxies, {len(gobs)} points; "
     f"templates t_p (transition-peaked) and t_e (deep-weighted)",
     f"G2 template overlap: corr(t_p, t_e) = {cc_pe:+.3f} over points "
     f"({'contest valid' if abs(cc_pe) < 0.85 else 'WARNING: templates degenerate'})",
     f"template shapes by bin (median t_p | t_e):"]
for b in range(6):
    mb = BIN_FID == b
    L.append(f"  bin{b}: {np.median(T_P[mb]):+7.4f} | {np.median(T_E[mb]):+7.4f}")
L.append("")

VARS = {'V0': (), 'VP': ('p',), 'VE': ('e',), 'VPE': ('p', 'e')}
res = {}
for vn, ch in VARS.items():
    res[vn] = fit_var(ch)
    gbv, cf = res[vn]
    sb = np.round(gbv.x[2:8], 4).tolist()
    L.append(f"{vn:>4} (+{'+'.join(ch) if ch else 'vertical only'}): obj = "
             f"{gbv.fun:9.2f}  s_b = {sb}")
L.append("")
s0 = res['V0'][0].x[2:8]; sp = res['VP'][0].x[2:8]
se_ = res['VE'][0].x[2:8]; spe = res['VPE'][0].x[2:8]
def bump(sb): return sb[4] - 0.5*(sb[3]+sb[5])
L.append(f"bump statistic s_b[4] - mean(s_b[3], s_b[5]): "
         f"V0 {bump(s0):+.4f} | VP {bump(sp):+.4f} | VE {bump(se_):+.4f} | "
         f"VPE {bump(spe):+.4f}")
L.append(f"deep bins (0+1) mean: V0 {(s0[0]+s0[1])/2:.4f} | VE "
         f"{(se_[0]+se_[1])/2:.4f}  (EFE-template absorption of the deep trend)")
cfP = res['VP'][1]
L.append(f"delta_p spread (VP): std = {np.std(cfP[:,1]):.4f} (prior {P_PRIOR})")
cfE = res['VE'][1]
L.append(f"delta_e spread (VE): std = {np.std(cfE[:,2]):.4f} (prior {E_PRIOR}); "
         f"note: fit symmetric, physical e_N is one-sided (disclosed)")
L.append("")

# G1 injection: synthetic p-scatter must produce a V0 bump and VP must kill it
rng = np.random.default_rng(23)
dp_t = rng.normal(0, P_PRIOR, NGal)
gN = g_gas + g_dsk + g_bul
lg_inj = (np.log10(gN*nu_be(gN/A0_FID)) + dp_t[gidx]*T_P
          + rng.normal(0, np.sqrt(sig2 + 0.05**2)))
gi0, _ = fit_var((), lg_override=lg_inj)
giP, cfi = fit_var(('p',), lg_override=lg_inj)
b_v0, b_vp = bump(gi0.x[2:8]), bump(giP.x[2:8])
okG1 = (b_v0 > 0.005) and (b_vp < 0.5*b_v0)
L.append(f"G1 injection (p-scatter 0.1): V0 bump {b_v0:+.4f} -> VP bump "
         f"{b_vp:+.4f}; corr(dp_hat, dp_true) = "
         f"{np.corrcoef(cfi[:,1], dp_t)[0,1]:.3f} -> {'PASS' if okG1 else 'FAIL'}")

out = "\n".join(L)
print(out)
with open('data/stage5a_bumpid.txt', 'w') as f:
    f.write(out+"\n")
