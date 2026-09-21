"""STAGE 10U -- PIN THE NU-FORM: the anchored-leg function-family contest.

Pre-registration: PREREG-NUFORM-10U.md (bars written before any contest
number existed). World + engine + family implementations copied
BIT-VERBATIM from calcs/stage10t_legregrow.py (that module exits in
gates mode on import, so verbatim copy + the G10U-1 identity regression
is the inheritance route). LT loader stubbed inert: lt_slot=None always
(G10T-3b standing), the LT block of make_world is never evaluated.

Modes: gates (sky-blind; sky member fits computed only inside the
G10U-1 regression with objective values MASKED -- a0 regressions
printed, funs never printed) | sky (the contest read).

Blind objects (never computed before this stage): the four members'
profiled -2lnL values and deltas on the 78 leg, their injection power
calibration, the paired-bootstrap delta distribution, any demotion
ruling. PRE-KNOWN (public in data/stage10t_skyread.txt): the four
members' deep-fit centrals (a0/H0) and the BE primary boot sigma.

R46 standing rules honored: health checks EVALUATED with printed
PASS/FAIL; sigma/boot engine-matched to its central; deltas printed at
0.1 precision, never quoted below the measured convergence spread.
"""
import glob, json, math, os, re, sys, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar, brentq

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
MODE = sys.argv[1] if len(sys.argv) > 1 else 'gates'
SKY = (MODE == 'sky')
OUTFILE = 'data/stage10u_skyread.txt' if SKY else 'data/stage10u_gates.txt'

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
LN10 = math.log(10)
C_LIGHT = 299792458.0
KMS_MPC = 3.240779e-20
A0_FID = 1.2e-10
S_ML = 0.1*LN10
U_PRIOR = (0.9/18.0)/LN10

def h0_of_a0(a0):
    return 2.0*math.pi*a0/C_LIGHT/KMS_MPC

def a0_of_h0(h0):
    return C_LIGHT*h0*KMS_MPC/(2.0*math.pi)

L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

def save():
    with open(OUTFILE, 'w', encoding='utf-8') as f:
        f.write("\n".join(L) + "\n")

t00 = time.time()
P(f"STAGE 10U PIN-THE-NU-FORM -- mode = {MODE} "
  f"({'sky read' if SKY else 'sky-blind gates'})")
P("")

# ---- the 10T archive (public, data/stage10t_skyread.txt) ----
MEMBERS = ('BE', 'p065', 'gm', 'boot')
MIDX = {m: i for i, m in enumerate(MEMBERS)}
ARC = dict(BE=1.0106e-10, p065=1.0892e-10, gm=1.0250e-10, boot=1.0446e-10)
ARC_BE = dict(f=1.050, s=0.032, u=-0.0211, H0=65.36)
ARC_BOOT = dict(sig=4.99, pct=(60.8, 65.4, 70.5))

# ---------------- frozen census (10T pre-reg R1) ----------------
RECLASS12 = ['ESO563-G021', 'NGC2903', 'NGC4559', 'NGC5585', 'UGC02487',
             'UGC05918', 'UGC05986', 'UGC07323', 'UGC07603', 'UGC08550',
             'UGC09992', 'UGC12632']
ZERO_PT = 'UGC09992'
UMA_CEPH = {'NGC3972': 21.23, 'NGC4051': 16.62}   # frozen co-read, Mpc

# ---------------- SPARC load (verbatim 10S/10T) ----------------
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

SP = dict(g_gas=[], g_dsk=[], g_bul=[], gobs=[], sig=[], gal_id=[],
          name={}, sigv={}, fd={}, dist={})
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
kept = 0
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc, fdv = meta.get(name, (0, 3, 10.0, 1.0, 3.0, 1))
    if inc < 30 or q > 2: continue
    kept += 1
    irad = math.radians(inc)
    sv = math.hypot((eD/max(D, 1e-3))/LN10,
                    2.0*(math.radians(max(einc, 1.0))/math.tan(irad))/LN10)
    SP['name'][gi] = name
    SP['sigv'][gi] = max(sv, 0.01)
    SP['fd'][gi] = fdv
    SP['dist'][gi] = (D, eD)
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC; gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        SP['g_gas'].append(gg); SP['g_dsk'].append(gd); SP['g_bul'].append(gb)
        SP['gobs'].append(Vo*Vo/R*KPC)
        SP['sig'].append(2*eV/Vo/LN10)
        SP['gal_id'].append(gi)
for k in ('g_gas', 'g_dsk', 'g_bul', 'gobs', 'sig', 'gal_id'):
    SP[k] = np.array(SP[k])
assert len(SP['gobs']) == 2700 and kept == 153
SPNAME2GI = {v: k for k, v in SP['name'].items()}

# ---------------- CF4 + cache (verbatim 10T) ----------------
MCOLS = dict(snia=(41, 47, 48, 52), tf=(53, 59, 60, 64),
             fp=(65, 71, 72, 76), sbf=(77, 83, 84, 89),
             snII=(90, 96, 97, 101), trgb=(102, 107, 108, 112),
             ceph=(113, 119, 120, 125), mas=(126, 131, 132, 136))
PRIO = ['mas', 'ceph', 'trgb', 'sbf', 'snia', 'snII']
cf4 = {}
for l in open('data/cf4/table2.dat'):
    try:
        pgc = int(l[0:7])
    except ValueError:
        continue
    row = {}
    for m, (a, b, c, d) in MCOLS.items():
        s1, s2 = l[a:b].strip(), l[c:d].strip()
        row[m] = (float(s1), float(s2)) if s1 and s2 else None
    cf4[pgc] = row
cache = json.load(open('data/scout10t_pgc_cache.json'))

def adopt(pgc):
    row = cf4.get(pgc)
    if not row: return None
    for m in PRIO:
        if row[m] is not None:
            dm, edm = row[m]
            D = 10**((dm - 25)/5)
            return (m, D, D*(LN10/5)*edm)
    return None

m31 = cf4[2557]
n4258 = cf4.get(39600, {})
chk_m31 = (m31['ceph'] == (24.397, 0.066) and m31['trgb'] == (24.53, 0.1)
           and m31['sbf'] == (24.214, 0.08))
chk_mas = (n4258.get('mas') is not None
           and 29.3 <= n4258['mas'][0] <= 29.5)
P(f"CF4 parse checks: M31 ceph/trgb/sbf -> {'PASS' if chk_m31 else 'FAIL'}; "
  f"NGC4258 DMmas -> {'PASS' if chk_mas else 'FAIL'}")

# ---------------- census (reclass identity; LT inert) ----------------
reclass_info = {}
for gi, nm in SP['name'].items():
    if SP['fd'][gi] != 1: continue
    p = cache.get(nm)
    a = adopt(p) if p else None
    if a:
        reclass_info[nm] = (gi, p) + a
rc_names = sorted(reclass_info)
zero_ok = int((SP['gal_id'] == reclass_info[ZERO_PT][0]).sum()) == 0
gcen_ok = chk_m31 and chk_mas and (rc_names == sorted(RECLASS12)) and zero_ok
P(f"census: reclass table = {len(rc_names)} on paper (frozen 12) "
  f"-> {'PASS' if gcen_ok else 'FAIL'}")
if not gcen_ok:
    P("STOP: census drift vs 10T frozen lists"); save(); sys.exit(0)

# LT stubs -- NEVER evaluated (lt_slot=None always in 10U; the LT block
# of make_world is guarded by `if lt_slot is not None`).
t1, lt_rows, lt_in = {}, {}, []

# ---------------- world builder (verbatim 10T) ----------------
def make_world(reclass=True, lt_slot='gas', lt_cut=0.10,
               uma_moved=False, uniform=False):
    """Assemble point arrays + maps. Base = verbatim SPARC arrays;
    distance moves via pin-3 (lgobs -= log10 s, g_bar invariant)."""
    lgobs = np.log10(SP['gobs']).copy()
    gg = SP['g_gas'].copy(); gd = SP['g_dsk'].copy(); gb = SP['g_bul'].copy()
    sig2 = (SP['sig']**2).copy()
    gal_id = SP['gal_id'].copy()
    sigv = dict(SP['sigv']); fd = dict(SP['fd'])
    rel = {}
    names = dict(SP['name'])
    for gi in SP['name']:
        D, eD = SP['dist'][gi]
        if fd[gi] == 4:
            rel[gi] = 2.3/18.0
        else:
            rel[gi] = eD/max(D, 1e-3)

    def move_gal(gi, Dnew, eDnew, Dold, newclass):
        s = Dnew/Dold
        pts = np.where(gal_id == gi)[0]
        lgobs[pts] -= math.log10(s)
        inc, q, D0, eD0, einc, fdv = meta[names[gi]]
        irad = math.radians(inc)
        sv = math.hypot((eDnew/max(Dnew, 1e-3))/LN10,
                        2.0*(math.radians(max(einc, 1.0)) /
                             math.tan(irad))/LN10)
        sigv[gi] = max(sv, 0.01)
        rel[gi] = eDnew/max(Dnew, 1e-3)
        fd[gi] = newclass

    if reclass:
        for nm in rc_names:
            gi, p, m, D, eD = reclass_info[nm]
            D0, _ = SP['dist'][gi]
            move_gal(gi, D, eD, D0, 2)
    if uma_moved or uniform:
        for nm, Dc in UMA_CEPH.items():
            gi = SPNAME2GI[nm]
            row = cf4[cache[nm]]
            dm, edm = row['ceph']
            Dc_ = 10**((dm-25)/5)
            move_gal(gi, Dc_, Dc_*(LN10/5)*edm, 18.0, 2)
    if uniform:
        for gi, nm in list(SP['name'].items()):
            if fd[gi] != 4: continue
            a = adopt(cache.get(nm))
            if a:
                m, D, eD = a
                move_gal(gi, D, eD, 18.0, 2)
        for gi, nm in SP['name'].items():
            if fd[gi] not in (2, 3, 5): continue
            if reclass and nm in reclass_info: continue
            p = cache.get(nm)
            a = adopt(p) if p else None
            if a:
                m, D, eD = a
                D0, _ = SP['dist'][gi]
                if abs(D/D0 - 1) > 1e-6:
                    move_gal(gi, D, eD, D0, fd[gi])

    # LT append (inert in 10U: lt_slot=None always)
    ltmap = {}
    if lt_slot is not None:
        nxt = max(SP['name']) + 1
        add = dict(gg=[], gd=[], gb=[], lg=[], s2=[], gid=[])
        for nm in lt_in:
            p, (m, D, eD), rings = lt_rows[nm]
            s = D/t1[nm]['D']
            gi = nxt; nxt += 1
            ltmap[gi] = nm
            names[gi] = 'LT:' + nm
            irad = math.radians(t1[nm]['i'])
            sv = math.hypot((eD/max(D, 1e-3))/LN10,
                            2.0*(math.radians(max(t1[nm]['ei'], 1.0)) /
                                 math.tan(irad))/LN10)
            sigv[gi] = max(sv, 0.01)
            rel[gi] = eD/max(D, 1e-3)
            fd[gi] = 9
            for (R, V, eV, Vb) in rings:
                if eV/V > lt_cut: continue
                gbar = Vb*Vb/R*KPC
                gobs_ = V*V/R*KPC
                if lt_slot == 'gas':
                    add['gg'].append(gbar); add['gd'].append(0.0)
                else:
                    add['gg'].append(0.0); add['gd'].append(gbar)
                add['gb'].append(0.0)
                add['lg'].append(math.log10(gobs_) - math.log10(s))
                add['s2'].append((2*eV/V/LN10)**2)
                add['gid'].append(gi)
        if add['gid']:
            gg = np.concatenate([gg, add['gg']])
            gd = np.concatenate([gd, add['gd']])
            gb = np.concatenate([gb, add['gb']])
            lgobs = np.concatenate([lgobs, add['lg']])
            sig2 = np.concatenate([sig2, add['s2']])
            gal_id = np.concatenate([gal_id,
                                     np.array(add['gid'], dtype=int)])

    ug = np.unique(gal_id)
    gpts = {g: np.where(gal_id == g)[0] for g in ug}
    FD = np.array([fd[g] for g in ug])
    W = dict(gg=gg, gd=gd, gb=gb, lgobs=lgobs, sig2=sig2, gal_id=gal_id,
             ug=ug, gpts=gpts, sigv=sigv, fd=fd, rel=rel, names=names,
             flow=ug[FD == 1], legA=ug[np.isin(FD, (2, 3, 5, 4, 9))],
             uma=ug[FD == 4], lt=ug[FD == 9])
    return W

W78 = make_world(reclass=True, lt_slot=None)
LEGA = list(W78['legA'])
assert len(W78['legA']) == 78 and len(W78['flow']) == 71
assert len(W78['uma']) == 26 and len(W78['lt']) == 0
P(f"world: primary leg A = {len(W78['legA'])} (41 anchored + 26 UMa + "
  f"11 reclassified; LT inert), flow = {len(W78['flow'])}")

# ---------------- nu families (verbatim 10T) ----------------
def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_simple(y):
    return 0.5+np.sqrt(0.25+1.0/np.clip(y, 1e-14, None)**2)
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
FAMS = {'BE': nu_be, 'p065': nu_p065, 'gm': nu_gm, 'boot': nu_boot}

# ---------------- fit machinery (verbatim 10T) ----------------
def build_sub(W, gal_seq, lg_src=None, dlg_per_occ=None):
    src = W['lgobs'] if lg_src is None else lg_src
    blocks, gidx_s, sv, ig = [], [], [], []
    lg_parts = []
    for k, g in enumerate(gal_seq):
        pts = W['gpts'][int(g)]
        blocks.append(pts)
        gidx_s.append(np.full(len(pts), k))
        sv.append(W['sigv'][int(g)])
        ig.append(1.0 if W['fd'][int(g)] == 4 else 0.0)
        shift = 0.0 if dlg_per_occ is None else dlg_per_occ[k]
        lg_parts.append(src[pts] + shift)
    pidx = np.concatenate(blocks)
    gidx_s = np.concatenate(gidx_s)
    n = len(gal_seq)
    ig = np.array(ig)
    return dict(n=n, gidx=gidx_s, sv=np.array(sv), isuma_g=ig,
                isuma_pt=ig[gidx_s],
                gg=W['gg'][pidx], gd=W['gd'][pidx], gb=W['gb'][pidx],
                lg=np.concatenate(lg_parts), s2=W['sig2'][pidx],
                gpts=[np.where(gidx_s == i)[0] for i in range(n)])

def fit_hier(sub, nu, use_u=False, upri=U_PRIOR, dlg=0.0, tol=0.05,
             max_rounds=15, th0=None):
    n = sub['n']
    lg = sub['lg'] + dlg
    gg, gd, gb_ = sub['gg'], sub['gd'], sub['gb']
    s2, gidx_s, sv = sub['s2'], sub['gidx'], sub['sv']
    ipt, ig = sub['isuma_pt'], sub['isuma_g']
    w_g = np.ones(n)
    dml = np.zeros(n)
    dv = np.zeros(n)

    def m2(th, dml, dv):
        if use_u:
            la0, f, s_int, u = th
            if abs(u) > 0.2: return 1e12
        else:
            la0, f, s_int = th; u = 0.0
        if not (-10.6 < la0 < -9.4) or not (0.3 < f < 2.5): return 1e12
        if not (1e-3 <= s_int < 0.4): return 1e12
        a0 = 10**la0
        fac = f*np.exp(dml[gidx_s])
        gN = gg + fac*gd + gb_
        gm_ = gN*nu(gN/a0)
        se2 = s2 + s_int*s_int
        r = lg - np.log10(gm_) - dv[gidx_s] - u*ipt
        out = np.sum(r*r/se2 + np.log(se2))
        out += np.sum(w_g*dml*dml)/(S_ML*S_ML)
        out += np.sum(w_g*dv*dv/(sv*sv))
        if use_u: out += (u/upri)**2
        return out

    best = None
    prev = None
    for rd in range(max_rounds):
        starts = ([list(best.x)] if best is not None else []) + \
                 ([list(th0)] if th0 is not None and best is None else []) + \
                 ([[math.log10(A0_FID), 1.0, 0.08] + ([0.0] if use_u else [])]
                  if rd == 0 else [])
        bb = None
        for t0 in starts:
            b = minimize(lambda t: m2(t, dml, dv), t0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
            if bb is None or b.fun < bb.fun: bb = b
        best = bb
        la0, f, s_int = best.x[:3]
        u = best.x[3] if use_u else 0.0
        se2c = s_int*s_int
        for _ in range(3):
            fac = f*np.exp(dml[gidx_s])
            gN = gg + fac*gd + gb_
            r0 = lg - np.log10(gN*nu(gN/10**la0)) - u*ipt
            for gi2 in range(n):
                mm = sub['gpts'][gi2]
                w = 1.0/(s2[mm] + se2c)
                dv[gi2] = np.sum(w*r0[mm])/(np.sum(w) + 1.0/sv[gi2]**2)
            for gi2 in range(n):
                mm = sub['gpts'][gi2]
                uoff = u*ig[gi2]
                def od(dl):
                    fc = f*math.exp(dl)
                    gN2 = gg[mm] + fc*gd[mm] + gb_[mm]
                    rr = (lg[mm] - np.log10(gN2*nu(gN2/10**la0))
                          - dv[gi2] - uoff)
                    ss = s2[mm] + se2c
                    return np.sum(rr*rr/ss) + dl*dl/(S_ML*S_ML)
                dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                           method='bounded').x
        cur = m2(best.x, dml, dv)
        if prev is not None and abs(prev - cur) < tol:
            prev = cur
            break
        prev = cur
    b = minimize(lambda t: m2(t, dml, dv), list(best.x),
                 method='Nelder-Mead',
                 options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
    if b.fun < best.fun: best = b
    return best

def fit_deep(sub, nu, th0=None, use_u=True):
    """R46-C2: the OPERATIVE estimator."""
    return fit_hier(sub, nu, use_u=use_u, th0=th0, tol=5e-4,
                    max_rounds=60)

def contest_fit(sub, nm, warm):
    """The pre-registered CONTEST estimator: best-of-two-starts
    fit_deep (cold + warm), keeping the smaller objective. Returns
    (best_fit, |fun_cold - fun_warm|)."""
    b1 = fit_deep(sub, FAMS[nm], th0=None)
    b2 = fit_deep(sub, FAMS[nm], th0=warm)
    gap = abs(b1.fun - b2.fun)
    return (b1 if b1.fun <= b2.fun else b2), gap

# ================= G10U-2: member identity (independent scalar) ======
P("")
P("-- G10U-2 member identity (independent scalar solves) --")
def nu_indep(nm, y):
    if nm == 'BE':
        return 1.0/(1.0 - math.exp(-math.sqrt(y)))
    if nm == 'p065':
        return (1.0 - math.exp(-y**0.65))**(-1.0/1.3)
    if nm == 'gm':
        # defining relation: nu = 1 + n_BE(y^0.75 * sqrt(nu)).
        # (verifier fix pre-run, before any gate result: the raw
        # bracket end overflowed expm1 -- guard the exponent; trap #23
        # class, checker-side only, stage engine untouched)
        def f(v):
            arg = y**0.75*math.sqrt(v)
            n = 1.0/math.expm1(arg) if arg < 700.0 else 0.0
            return v - 1.0 - n
        return brentq(f, 1.0 + 1e-13, 1e6, xtol=1e-14, rtol=1e-15)
    if nm == 'boot':
        # defining relation: u = y + y/(e^u - 1), nu = u/y
        f = lambda u: u - y - y/math.expm1(u)
        return brentq(f, y + 1e-13, y + 60.0, xtol=1e-15, rtol=1e-15)/y

YGRID = (0.03, 0.1, 0.3, 1.0, 3.0, 10.0)
g2u_ok = True
for nm in MEMBERS:
    ds = []
    vals = [float(np.atleast_1d(FAMS[nm](np.array([y])))[0]) for y in YGRID]
    for y, v in zip(YGRID, vals):
        vi = nu_indep(nm, y)
        ds.append(abs(v - vi)/vi)
    mono = all(vals[i] > vals[i+1] for i in range(len(vals)-1))
    ok = max(ds) <= 1e-6 and mono
    g2u_ok &= ok
    P(f"  {nm:5s}: max |d nu|/nu = {max(ds):.2e} (bar 1e-6), "
      f"monotone-decreasing {'yes' if mono else 'NO'} -> "
      f"{'PASS' if ok else 'FAIL'}")
P(f"G10U-2: {'PASS' if g2u_ok else 'FAIL'}")
if not g2u_ok:
    P("STOP: member identity failed"); save(); sys.exit(0)

# ================= G10U-1: regression (archived convention) ==========
P("")
P("-- G10U-1 regression (archived 10T deep fits; funs MASKED) --")
subP = build_sub(W78, LEGA)
b_plain = fit_hier(subP, nu_be, use_u=False)
b_loose = fit_hier(subP, nu_be, use_u=True, th0=list(b_plain.x)+[0.0])
bT_arch = fit_deep(subP, nu_be, th0=list(b_loose.x))
arch_fit = {'BE': bT_arch}
for nm in ('p065', 'gm', 'boot'):
    arch_fit[nm] = fit_deep(subP, FAMS[nm], th0=list(bT_arch.x))
g1u_ok = True
for nm in MEMBERS:
    a0 = 10**arch_fit[nm].x[0]
    d = abs(a0 - ARC[nm])/ARC[nm]
    g1u_ok &= d <= 1e-3
    P(f"  {nm:5s}: a0 = {a0:.4e} vs archived {ARC[nm]:.4e} "
      f"(d = {100*d:.3f}% vs bar 0.1%)")
P(f"  BE params: f_ML = {bT_arch.x[1]:.3f} ({ARC_BE['f']:.3f}), "
  f"s_int = {bT_arch.x[2]:.3f} ({ARC_BE['s']:.3f}), "
  f"u = {bT_arch.x[3]:+.4f} ({ARC_BE['u']:+.4f})")
P(f"G10U-1: {'PASS' if g1u_ok else 'FAIL'}  [{time.time()-t00:.0f}s]")
if not g1u_ok:
    P("STOP: regression vs 10T archive failed"); save(); sys.exit(0)

# ================= G10U-3: the POWER arm ==============================
P("")
P("-- G10U-3 power (4 truths x 8 seeds, G10T-4 generative model) --")
SEEDS8 = (11, 23, 42, 101, 202, 303, 404, 505)
gN1 = W78['gg'] + 1.0*W78['gd'] + W78['gb']
SIGPT = np.sqrt(W78['sig2'] + 0.08**2)

def make_mock(nm, seed):
    """VERBATIM G10T-4 noisy pattern with nu_T, a0_T = the member's own
    public 10T central; rng stream 1000*seed + 100 + member index
    (distinct from the 10T streams)."""
    rng = np.random.default_rng(1000*seed + 100 + MIDX[nm])
    base = np.log10(gN1*FAMS[nm](gN1/ARC[nm]))
    mock = base + rng.normal(0, SIGPT)
    voff = [rng.normal(0, W78['sigv'][g]) for g in LEGA]
    return mock, voff

FUNS = {}      # FUNS[truth][seed][member] = contest fun
GAPS = {}      # start gaps per (truth, seed, member)
for T in MEMBERS:
    FUNS[T] = {}
    warmT = [math.log10(ARC[T]), 1.0, 0.08, 0.0]
    for seed in SEEDS8:
        mock, voff = make_mock(T, seed)
        subI = build_sub(W78, LEGA, lg_src=mock, dlg_per_occ=voff)
        row = {}
        for m in MEMBERS:
            bf, gap = contest_fit(subI, m, warmT)
            row[m] = bf.fun
            GAPS[(T, seed, m)] = gap
        FUNS[T][seed] = row
    P(f"  truth {T:5s} worlds fit  [{time.time()-t00:.0f}s]")

def dstat(A, B, T):
    """d_AB = fun_B - fun_A over the truth-T seeds."""
    return np.array([FUNS[T][s][B] - FUNS[T][s][A] for s in SEEDS8])

# self-recovery (STOP-grade, trap #26)
g3u_self = True
P("  self-recovery (truth never significantly loses; bar mean >= -2SE):")
for T in MEMBERS:
    parts = []
    for R in MEMBERS:
        if R == T: continue
        d = dstat(T, R, T)
        se = d.std(ddof=1)/math.sqrt(len(d))
        ok = d.mean() >= -2*se
        g3u_self &= ok
        parts.append(f"{R} {d.mean():+.1f}+-{d.std(ddof=1):.1f}"
                     f"{'' if ok else ' FAIL'}")
    P(f"    truth {T:5s}: d(T beats R) = " + ", ".join(parts))
P(f"  self-recovery -> {'PASS' if g3u_self else 'FAIL'}")
if not g3u_self:
    P("STOP: contest instrument broken (truth loses to a rival)")
    save(); sys.exit(0)

# separations (R43 rule: NP separation measured before any bar)
P("  separation matrix sep(A,B) (pair powered iff >= 2.0; 8+8 seeds,")
P("  ~0.35*sep sampling noise -- a measured, lean-grade boundary):")
SEP = {}
POWERED = set()
for A in MEMBERS:
    cells = []
    for B in MEMBERS:
        if A == B:
            cells.append("   .  "); continue
        dA = dstat(A, B, A)   # truth A: expect positive
        dB = dstat(A, B, B)   # truth B: expect negative
        sep = (dA.mean() - dB.mean())/math.sqrt(
            (dA.std(ddof=1)**2 + dB.std(ddof=1)**2)/2.0)
        SEP[(A, B)] = sep
        if sep >= 2.0: POWERED.add((A, B))
        cells.append(f"{sep:6.2f}")
    P(f"    {A:5s} vs " + " ".join(f"{B:>6s}" for B in MEMBERS))
    P(f"    {'':5s}    " + " ".join(cells))
P(f"  powered ordered pairs (sep >= 2.0): "
  f"{sorted(POWERED) if POWERED else 'NONE'}")

# truth-side calibration stats for the demotion bar
CAL = {}
for A in MEMBERS:
    for B in MEMBERS:
        if A == B: continue
        dB = dstat(A, B, B)
        CAL[(A, B)] = (dB.mean(), dB.std(ddof=1))

# ---- AMENDMENT A1 (pre-sky, pre-quote; PREREG-NUFORM-10U.md):
# boundary-pair seed extension -- any unordered pair with
# |sep8 - 2.0| <= 0.5 gets +16 seeds/truth side; POOLED 24-seed sep
# and calibration REPLACE the 8-seed values for those pairs.
SEEDS16 = (606, 707, 808, 909, 1010, 1111, 1212, 1313, 1414, 1515,
           1616, 1717, 1818, 1919, 2020, 2121)
EXT_PAIRS = sorted({tuple(sorted((A, B))) for (A, B) in SEP
                    if abs(SEP[(A, B)] - 2.0) <= 0.5})
if EXT_PAIRS:
    P("")
    P("-- A1 boundary-pair seed extension (pre-sky amendment) --")
    P(f"  pairs with |sep8 - 2.0| <= 0.5: {EXT_PAIRS}")
    need = {}
    for (A, B) in EXT_PAIRS:
        for Tt in (A, B):
            need.setdefault(Tt, set()).update((A, B))
    FUNS_E = {T: {} for T in need}
    for T in sorted(need):
        warmT = [math.log10(ARC[T]), 1.0, 0.08, 0.0]
        for seed in SEEDS16:
            mock, voff = make_mock(T, seed)
            subI = build_sub(W78, LEGA, lg_src=mock, dlg_per_occ=voff)
            row = {}
            for m in sorted(need[T]):
                bf, _g = contest_fit(subI, m, warmT)
                row[m] = bf.fun
            FUNS_E[T][seed] = row
        P(f"  truth {T:5s} +16 worlds fit  [{time.time()-t00:.0f}s]")
    def dstat_pool(A, B, T):
        if T not in FUNS_E:
            return dstat(A, B, T)
        dE = np.array([FUNS_E[T][s][B] - FUNS_E[T][s][A]
                       for s in SEEDS16])
        return np.concatenate([dstat(A, B, T), dE])
    for (A, B) in EXT_PAIRS:
        dA = dstat_pool(A, B, A)
        dB = dstat_pool(A, B, B)
        sep = (dA.mean() - dB.mean())/math.sqrt(
            (dA.std(ddof=1)**2 + dB.std(ddof=1)**2)/2.0)
        P(f"  {A}-{B}: sep8 = {SEP[(A, B)]:.2f} -> sep24 = {sep:.2f} "
          f"(dA {dA.mean():+.1f}+-{dA.std(ddof=1):.1f}, dB "
          f"{dB.mean():+.1f}+-{dB.std(ddof=1):.1f}) -> "
          f"{'POWERED' if sep >= 2.0 else 'unpowered'}")
        for (X, Y) in ((A, B), (B, A)):
            SEP[(X, Y)] = sep
            if sep >= 2.0: POWERED.add((X, Y))
            else: POWERED.discard((X, Y))
            dYp = dstat_pool(X, Y, Y)
            CAL[(X, Y)] = (dYp.mean(), dYp.std(ddof=1))
    P(f"  powered ordered pairs (final, A1-pooled): {sorted(POWERED)}")

# ================= G10U-4: convergence repeat =========================
P("")
c_conv = max(GAPS[('BE', 11, m)] for m in MEMBERS)
gap_all = max(GAPS.values())
P(f"-- G10U-4 convergence repeat --")
P(f"  c_conv (truth-BE seed-11 world, max |fun_cold - fun_warm| over "
  f"members) = {c_conv:.2f} (bar 2.0)")
P(f"  cross-world max start gap (diagnostic) = {gap_all:.2f}")
g4u_ok = c_conv <= 2.0
P(f"G10U-4: {'PASS' if g4u_ok else 'FAIL'}")
if not g4u_ok:
    P("STOP: cross-member deltas unquotable at this convergence grade")
    save(); sys.exit(0)

# ================= G10U-6: N-forecast (diagnostic) ====================
P("")
P("-- G10U-6 N-forecast (DIAGNOSTIC/FORECAST; BE vs p065; cold fits) --")
def forecast_sub(nm, seed, nmult):
    """N-fold duplicated 78 world with independent noise per copy."""
    rng = np.random.default_rng(7000*seed + 700 + MIDX[nm] + 10*nmult)
    base = np.log10(gN1*FAMS[nm](gN1/ARC[nm]))
    blocks, gidx_s, sv, ig, lg_parts, s2p = [], [], [], [], [], []
    k = 0
    for c in range(nmult):
        for g in LEGA:
            pts = W78['gpts'][int(g)]
            blocks.append(pts)
            gidx_s.append(np.full(len(pts), k))
            sv.append(W78['sigv'][int(g)])
            ig.append(1.0 if W78['fd'][int(g)] == 4 else 0.0)
            voff = rng.normal(0, W78['sigv'][int(g)])
            lg_parts.append(base[pts] + rng.normal(0, SIGPT[pts]) + voff)
            s2p.append(W78['sig2'][pts])
            k += 1
    pidx = np.concatenate(blocks)
    gidx_s = np.concatenate(gidx_s)
    ig = np.array(ig)
    return dict(n=k, gidx=gidx_s, sv=np.array(sv), isuma_g=ig,
                isuma_pt=ig[gidx_s],
                gg=W78['gg'][pidx], gd=W78['gd'][pidx],
                gb=W78['gb'][pidx],
                lg=np.concatenate(lg_parts), s2=np.concatenate(s2p),
                gpts=[np.where(gidx_s == i)[0] for i in range(k)])

FSEEDS = (11, 23, 42, 101)
fc = {}
for nmult in (1, 2, 4):
    dd = {}
    for T in ('BE', 'p065'):
        ds = []
        for seed in FSEEDS:
            subF = forecast_sub(T, seed, nmult)
            fB = fit_deep(subF, nu_be).fun
            fP = fit_deep(subF, nu_p065).fun
            ds.append(fP - fB)   # d(BE, p065)
        dd[T] = np.array(ds)
    sep = (dd['BE'].mean() - dd['p065'].mean())/math.sqrt(
        (dd['BE'].std(ddof=1)**2 + dd['p065'].std(ddof=1)**2)/2.0)
    fc[nmult] = sep
    P(f"  N x{nmult} (N = {78*nmult}): d|truthBE = "
      f"{dd['BE'].mean():+.1f}+-{dd['BE'].std(ddof=1):.1f}, d|truthP = "
      f"{dd['p065'].mean():+.1f}+-{dd['p065'].std(ddof=1):.1f} -> "
      f"sep = {sep:.2f}  [{time.time()-t00:.0f}s]")
# A1(3): the x4 point is flagged convergence-suspect (4-seed SD
# blow-up + rounds-cap risk on 312-galaxy cold fits); the sqrt(N)
# coefficient uses the x1-x2 points only, x4 printed as-is above.
P("  [A1(3) flag: the x4 row is convergence-suspect -- 4-seed SD "
  "noise + rounds-cap risk on 312-galaxy cold fits; excluded from "
  "the sqrt(N) fit]")
c1 = np.mean([fc[nm]/math.sqrt(nm) for nm in (1, 2)])
if c1 > 0:
    P(f"  sqrt(N) law (x1-x2): sep ~ {c1:.2f} x sqrt(N/78); sep = 2 "
      f"at N ~ {78*(2.0/c1)**2:.0f}, sep = 3 at N ~ "
      f"{78*(3.0/c1)**2:.0f} (FORECAST grade; BIG-SPARC requirement)")
else:
    P("  sqrt(N) law: coefficient non-positive at these seeds -- "
      "forecast not extractable (FORECAST grade)")

# ================= gate record ========================================
P("")
P(f"GATES: G10U-1 PASS  G10U-2 PASS  G10U-3self PASS  G10U-4 PASS "
  f"(c_conv = {c_conv:.2f})  powered pairs = {len(POWERED)}/12")
if not SKY:
    P("")
    P("gates mode complete -- sky member objectives were computed only "
      "inside the G10U-1 regression and are MASKED (a0 regressions "
      "printed, funs never printed; R46-C12 pattern). No contest "
      "delta, no paired bootstrap, no ruling.")
    P(f"wall-clock: {(time.time()-t00)/60:.1f} min")
    save()
    print("\nsaved:", OUTFILE)
    sys.exit(0)

# ================= SKY: the contest read ==============================
P("")
P("== 10U SKY READ: the four-member contest on the anchored 78 ==")
sky_fit, sky_gap = {}, {}
b_cold_be = fit_deep(subP, nu_be, th0=None)
sky_fit['BE'] = b_cold_be if b_cold_be.fun < bT_arch.fun else bT_arch
sky_gap['BE'] = abs(b_cold_be.fun - bT_arch.fun)
for nm in ('p065', 'gm', 'boot'):
    b_warm = arch_fit[nm]
    b_cold = fit_deep(subP, FAMS[nm], th0=None)
    sky_fit[nm] = b_cold if b_cold.fun < b_warm.fun else b_warm
    sky_gap[nm] = abs(b_cold.fun - b_warm.fun)
fmin = min(sky_fit[m].fun for m in MEMBERS)
B_sky = min(MEMBERS, key=lambda m: sky_fit[m].fun)
P("member   a0          H0     fun-min   startgap")
for m in MEMBERS:
    a0 = 10**sky_fit[m].x[0]
    h0 = h0_of_a0(a0)
    drift = abs(h0 - h0_of_a0(ARC[m]))
    note = "  [INSTRUMENT NOTE: contest optimum drifted >0.3 vs archive]" \
        if drift > 0.3 else ""
    P(f"  {m:5s}  {a0:.4e}  {h0:5.1f}  {sky_fit[m].fun-fmin:8.1f}  "
      f"{sky_gap[m]:6.2f}{note}")
P(f"sky-best member: {B_sky}")
d_sky = {}
for A in MEMBERS:
    for B in MEMBERS:
        if A != B:
            d_sky[(A, B)] = sky_fit[B].fun - sky_fit[A].fun

# ---------- G10U-5 paired bootstrap (stream-verbatim sigma_boot_hier) --
P("")
P(f"-- G10U-5 paired bootstrap (200 reps, seed 202, R46 stream) "
  f"[{time.time()-t00:.0f}s] --")
rng5 = np.random.default_rng(202)
rel_sig, uma_flag = {}, {}
for g in W78['legA']:
    if W78['fd'][g] == 4:
        rel_sig[g] = 2.3/18.0; uma_flag[g] = True
    else:
        rel_sig[g] = W78['rel'][g]; uma_flag[g] = False
TH0B = {m: list(sky_fit[m].x) for m in MEMBERS}
la0_be_reps = []
fun_reps = {m: [] for m in MEMBERS}
for _ in range(200):
    pick = rng5.choice(W78['legA'], size=len(W78['legA']), replace=True)
    sc_shared = 1.0 + rng5.normal(0, 0.9/18.0)
    docc = []
    for g in pick:
        s_g = 1.0 + rng5.normal(0, rel_sig[g])
        s_g = min(max(s_g, 0.5), 1.5)
        if uma_flag[g]: s_g *= sc_shared
        docc.append(-math.log10(s_g))
    subr = build_sub(W78, list(pick), dlg_per_occ=docc)
    bb = fit_deep(subr, nu_be, th0=TH0B['BE'])
    la0_be_reps.append(bb.x[0])
    fun_reps['BE'].append(bb.fun)
    for m in ('p065', 'gm', 'boot'):
        fun_reps[m].append(fit_deep(subr, FAMS[m], th0=TH0B[m]).fun)
for m in MEMBERS:
    fun_reps[m] = np.array(fun_reps[m])
h0_be_reps = h0_of_a0(10**np.array(la0_be_reps))
sig_be = float(np.std(h0_be_reps))
pct = np.percentile(h0_be_reps, [16, 50, 84])
d5a_sig = abs(sig_be - ARC_BOOT['sig'])/ARC_BOOT['sig']
d5a_pct = max(abs(pct[i] - ARC_BOOT['pct'][i]) for i in range(3))
g5a_ok = (d5a_sig <= 0.005) and (d5a_pct <= 0.2)
P(f"  G10U-5a BE-marginal stream identity: sigma = {sig_be:.2f} vs "
  f"archived {ARC_BOOT['sig']:.2f} (d = {100*d5a_sig:.2f}% vs 0.5%); "
  f"pct = {pct.round(1).tolist()} vs {list(ARC_BOOT['pct'])} "
  f"(max d = {d5a_pct:.2f} vs 0.2) -> {'PASS' if g5a_ok else 'FAIL'}")
PBOOT = {}
for A in MEMBERS:
    for B in MEMBERS:
        if A == B: continue
        dr = fun_reps[B] - fun_reps[A]
        PBOOT[(A, B)] = float(np.mean(dr > 0))
runner = min((m for m in MEMBERS if m != B_sky),
             key=lambda m: sky_fit[m].fun)
dr_top = fun_reps[runner] - fun_reps[B_sky]
hgap = abs(float(np.median(dr_top)) - d_sky[(B_sky, runner)])
hsd = float(np.std(dr_top))
g5b_ok = hgap <= 1.0*hsd
P(f"  G10U-5b health (EVALUATED; top pair {B_sky} vs {runner}): "
  f"|boot median d - sky d| = {hgap:.1f} vs 1.0 x SD_boot = {hsd:.1f} "
  f"-> {'PASS' if g5b_ok else 'FAIL'}")
P("  paired-bootstrap win fractions P(A beats B):")
for A in MEMBERS:
    row = []
    for B in MEMBERS:
        row.append("  .  " if A == B else f"{PBOOT[(A, B)]:.3f}")
    P(f"    {A:5s}: " + "  ".join(f"{c:>6s}" for c in row))

# ---------- demotion evaluation (pre-registered rule) -----------------
P("")
P("-- demotion evaluation (rule: powered AND sky d > m_B + 2 s_B + "
  f"c_conv [{c_conv:.2f}] AND P_boot >= 0.95, health PASS) --")
demoted = {}
for B in MEMBERS:
    for A in MEMBERS:
        if A == B: continue
        if (A, B) not in POWERED: continue
        mB, sB = CAL[(A, B)]
        bar = mB + 2*sB + c_conv
        hit = (d_sky[(A, B)] > bar) and (PBOOT[(A, B)] >= 0.95) and g5b_ok
        P(f"  {A} -> {B}: powered (sep {SEP[(A, B)]:.2f}); sky d = "
          f"{d_sky[(A, B)]:+.1f} vs bar {bar:+.1f} "
          f"(truth-{B}: {mB:+.1f}+-{sB:.1f}); P_boot = "
          f"{PBOOT[(A, B)]:.3f} -> {'DEMOTED' if hit else 'no'}")
        if hit:
            demoted.setdefault(B, []).append(A)
if not POWERED:
    P("  no powered ordered pair -- demotion machinery NOT applied "
      "(pre-registered U-POWER-LIMITED branch); sky deltas above are "
      "point deltas, no calibrated grade.")
survivors = [m for m in MEMBERS if m not in demoted]
h0s = {m: h0_of_a0(10**sky_fit[m].x[0]) for m in MEMBERS}
band = (min(h0s[m] for m in survivors), max(h0s[m] for m in survivors))

# ---------- letter ----------------------------------------------------
if not POWERED:
    letter = 'U-POWER-LIMITED'
elif demoted and len(survivors) == 1:
    letter = 'U-FORM-PINNED'
elif demoted:
    letter = 'U-FORM-NARROWED'
else:
    letter = 'U-FORM-LEAN'
P("")
P("== LETTER ==")
P(f"{letter}: the anchored-leg four-member contest at N = 78.")
if letter == 'U-POWER-LIMITED':
    P("  No ordered pair reaches sep 2.0 at leg-A noise: the anchored "
      "leg CANNOT distinguish the registered family members at any "
      "calibrated grade. The nu-form band 65.4-70.4 STANDS as the "
      "measured floor of the meter; the deliverables are the power "
      "table and the N-forecast (the BIG-SPARC requirement).")
elif letter == 'U-FORM-LEAN':
    P(f"  Powered pairs exist ({sorted(POWERED)}) but no demotion "
      f"fires: the band 65.4-70.4 STANDS; the leans are quoted at "
      f"bootstrap grade only (win fractions above).")
else:
    P(f"  demoted: " +
      "; ".join(f"{B} (by {'/'.join(As)})" for B, As in demoted.items()))
    P(f"  survivors: {survivors} -> nu-form band {band[0]:.1f}-"
      f"{band[1]:.1f} km/s/Mpc (supersedes the 10T band clause per "
      f"pre-reg pin 7; REPORT grade, multiplicity disclosed)")
P("  In every branch: the 10T stat 5.0, membership band, flat "
  "treatment variant, peg clause, and validity domain are UNTOUCHED.")
P("")
P("MANDATORY DISCLOSURES:")
P(f"  precision rule: deltas printed at 0.1; the measured convergence "
  f"spread c_conv = {c_conv:.2f} (cross-world max {gap_all:.2f}) is "
  f"charged inside every demotion bar; never quote sub-c_conv deltas.")
P(f"  multiplicity: 3 rival tests per member at 2-SD grade; bootstrap "
  f"co-requirement is the family-wise control; any demotion is "
  f"REPORT-grade, sharpenable by seed extension at review.")
P(f"  expectation disclosure (pre-reg SS0): full-catalog leans toward "
  f"sharper members were on record BEFORE this stage; they were not "
  f"imported as priors or bars; this contest is standalone on the "
  f"anchored 78.")
P(f"  sky-best = {B_sky}; per-member H0: " +
  ", ".join(f"{m} {h0s[m]:.1f}" for m in MEMBERS))
P("")
P(f"verdict letter: {letter}")
P("credences: NO cell moves on any outcome (pre-registered); 53/8 "
  "untouched")
P(f"total wall-clock: {(time.time()-t00)/60:.1f} min")
save()
print("\nsaved:", OUTFILE)
