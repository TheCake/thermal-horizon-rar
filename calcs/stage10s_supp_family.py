"""10S SUPPLEMENT (labeled, post-pin, direction-blind): leg-B HIER
family crossings {p065, gm, boot}. Pin 9 pinned flat variants on leg B
(all no-bracket); the joint H0's family band should not ride one
function. Machinery imported verbatim from the sky-read module pattern
(re-declared; regression: BE crossing must reproduce 65.56).
Writes data/stage10s_supp_family.txt."""
import glob, math, os, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
LN10 = math.log(10)
C_LIGHT = 299792458.0
KMS_MPC = 3.240779e-20
A0_FID = 1.2e-10
S_ML = 0.1*LN10

def a0_of_h0(h0):
    return C_LIGHT*h0*KMS_MPC/(2.0*math.pi)

meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta[t[0]] = (float(t[5]), int(t[17]), float(t[2]), float(t[3]),
                      float(t[6]), int(t[4]))
    except ValueError:
        continue
g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
sigv_g_map, fd_g_map = {}, {}
for gi, path in enumerate(sorted(glob.glob(
        'data/sparc/rotmod/**/*_rotmod.dat', recursive=True))):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc, fdv = meta.get(name, (0, 3, 10.0, 1.0, 3.0, 1))
    if inc < 30 or q > 2: continue
    irad = math.radians(inc)
    sv = math.hypot((eD/max(D, 1e-3))/LN10,
                    2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10)
    sigv_g_map[gi] = max(sv, 0.01)
    fd_g_map[gi] = fdv
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC; gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        g_gas.append(gg); g_dsk.append(gd); g_bul.append(gb)
        gobs.append(Vo*Vo/R*KPC); sig.append(2*eV/Vo/LN10)
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs = np.log10(gobs)
ug = np.unique(gal_id)
flow_gal = ug[np.array([fd_g_map[g] for g in ug]) == 1]
gpts_all = {g: np.where(gal_id == g)[0] for g in ug}

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None))
def nu_p065(y):
    yc = np.clip(np.asarray(y, float), 1e-14, None)
    ex = np.exp(-np.minimum(yc**0.65, 60.0))
    return (1.0-ex)**(-1.0/1.3)
def nu_gm(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    a = y**0.75
    w = np.sqrt(nu_simple(y))
    for _ in range(30):
        u = np.minimum(a*w, 60.0)
        eu = np.exp(u)
        em1 = np.maximum(eu - 1.0, 1e-300)
        n = np.where(u < 60.0, 1.0/em1, 0.0)
        H = w*w - 1.0 - n
        dH = 2.0*w + a*eu/(em1*em1)
        w = np.maximum(w - H/dH, 1e-8)
    return w*w
def nu_boot(y):
    y = np.clip(np.asarray(y, float), 1e-14, None)
    u = 0.5*(y + np.sqrt(y*y + 4.0*y))
    for _ in range(14):
        uc = np.minimum(u, 45.0)
        eu = np.exp(uc)
        em1 = eu - 1.0
        F = u - y - y/em1
        dF = 1.0 + y*eu/(em1*em1)
        u = np.maximum(u - F/dF, 1e-13)
    nu = u/y
    big = y > 45.0
    if np.any(big):
        yb = np.minimum(y, 700.0)
        nu = np.where(big, 1.0 + 1.0/np.expm1(yb), nu)
    return nu

def build_sub(gal_seq):
    blocks, gidx_s, sv = [], [], []
    for k, g in enumerate(gal_seq):
        pts = gpts_all[int(g)]
        blocks.append(pts)
        gidx_s.append(np.full(len(pts), k))
        sv.append(sigv_g_map[int(g)])
    pidx = np.concatenate(blocks)
    gidx_s = np.concatenate(gidx_s)
    n = len(gal_seq)
    return dict(n=n, gidx=gidx_s, sv=np.array(sv),
                gg=g_gas[pidx], gd=g_dsk[pidx], gb=g_bul[pidx],
                lg=lgobs[pidx].copy(), s2=sig2[pidx],
                gpts=[np.where(gidx_s == i)[0] for i in range(n)])

def fit_hier(sub, nu, dlg=0.0, tol=0.05, max_rounds=15, th0=None):
    n = sub['n']
    lg = sub['lg'] + dlg
    gg, gd, gb_ = sub['gg'], sub['gd'], sub['gb']
    s2, gidx_s, sv = sub['s2'], sub['gidx'], sub['sv']
    w_g = np.ones(n)
    dml = np.zeros(n)
    dv = np.zeros(n)
    def m2(th, dml, dv):
        la0, f, s_int = th
        if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
        if not (1e-3 <= s_int < 0.4): return 1e12
        a0 = 10**la0
        fac = f*np.exp(dml[gidx_s])
        gN = gg + fac*gd + gb_
        gm_ = gN*nu(gN/a0)
        se2 = s2 + s_int*s_int
        r = lg - np.log10(gm_) - dv[gidx_s]
        out = np.sum(r*r/se2 + np.log(se2))
        out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
        out += np.sum(w_g*dv*dv/(sv*sv))
        return out
    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08]] if rd == 0 else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2(t, dml, dv), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int = best.x
        se2c = s_int*s_int
        for _ in range(3):
            fac = f*np.exp(dml[gidx_s])
            gN = gg + fac*gd + gb_
            r0 = lg - np.log10(gN*nu(gN/10**la0))
            for gi2 in range(n):
                mm = sub['gpts'][gi2]
                w = 1.0/(s2[mm] + se2c)
                dv[gi2] = np.sum(w*r0[mm])/(np.sum(w) + 1.0/sv[gi2]**2)
            for gi2 in range(n):
                mm = sub['gpts'][gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = gg[mm] + fc*gd[mm] + gb_[mm]
                    rr = (lg[mm] - np.log10(gN2*nu(gN2/10**la0))
                          - dv[gi2])
                    ss = s2[mm] + se2c
                    return np.sum(rr*rr/ss) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2(best.x, dml, dv)
        if prev is not None and abs(prev - cur) < tol:
            break
        prev = cur
    b = minimize(lambda t: m2(t, dml, dv), list(best.x),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best

def crossing(sub, nu):
    phi = {}
    th_warm = None
    def node(h0p):
        nonlocal th_warm
        b = fit_hier(sub, nu, dlg=-math.log10(73.0/h0p), th0=th_warm)
        th_warm = list(b.x)
        return math.log(10**b.x[0]) - math.log(a0_of_h0(h0p))
    for h0p in (60.0, 65.0, 70.0, 75.0, 80.0):
        phi[h0p] = node(h0p)
    def brk(ph):
        ks = sorted(ph.keys())
        return any(ph[ks[i]]*ph[ks[i+1]] <= 0 for i in range(len(ks)-1))
    while not brk(phi):
        # ROUND-45 FIX (condition 2): slope-aware extension (the
        # original rule censored high roots in the falling-phi regime;
        # this is why the p065 row printed a false "NO CROSSING").
        ks = sorted(phi.keys())
        slope = phi[ks[-1]] - phi[ks[0]]
        if phi[ks[0]] > 0: go_up = slope < 0
        else: go_up = slope > 0
        if go_up and ks[-1] < 95.0: nxt = ks[-1] + 5.0
        elif (not go_up) and ks[0] > 50.0: nxt = ks[0] - 5.0
        else: return None
        phi[nxt] = node(nxt)
    ks = sorted(phi.keys())
    for i in range(len(ks)-1):
        if phi[ks[i]]*phi[ks[i+1]] <= 0:
            x0_, x1_ = ks[i], ks[i+1]
            return x0_ + (x1_-x0_)*phi[x0_]/(phi[x0_]-phi[x1_])
    return None

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

subB = build_sub(list(flow_gal))
P("10S SUPPLEMENT: leg-B HIER family crossings (labeled, post-pin)")
for nm, nu in (('BE', nu_be), ('p065', nu_p065), ('gm', nu_gm),
               ('boot', nu_boot)):
    t0 = time.time()
    c = crossing(subB, nu)
    P(f"  {nm:>5}: crossing = "
      f"{('%.2f km/s/Mpc' % c) if c else 'NO CROSSING in [50,95]'}"
      f"  [{time.time()-t0:.0f}s]")
P("(BE row = regression vs the sky read's 65.56)")
with open('data/stage10s_supp_family.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("saved: data/stage10s_supp_family.txt")
