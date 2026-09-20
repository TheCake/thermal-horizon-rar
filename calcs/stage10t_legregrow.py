"""
STAGE 10T -- THE LEG-A REGROW (pre-reg PREREG-H0METER-10T.md; parent
10S 2a461d5/fb74226/98afffb; Round-45 verdict standing).

Grows the H0 meter's anchored leg 67 -> 82 galaxies: 11 SPARC flow
galaxies reclassified on CF4 per-method anchors (R1/R2; UGC09992 is a
12th on paper with zero surviving points), the 26-galaxy UMa block kept
verbatim (R4; NGC3972/4051 Cepheids = frozen co-read), 4 LITTLE THINGS
dwarfs through the 9R-validated port in the gas slot (R5/R7). Leg B is
NOT refired (R3; R45 standing). Modes:
  py calcs/stage10t_legregrow.py gates  -> data/stage10t_gates.txt
     (sky-blind: no extended-set central computed or printed anywhere)
  py calcs/stage10t_legregrow.py sky    -> data/stage10t_skyread.txt
     (re-verifies gates in-process, then fires the extended leg A)

G10T-6 literature log: three independent NOT-FOUNDs on any published
a0->H0 inversion-as-measurement are on file as of 2026-09-20 (the 10S
G6 round + the morning scout + the data-sweep delta scout, all run
2026-09-20; SML20 = nearest neighbor, must-cite; van Putten = forward
direction only; CF4's own H0 = 74.6 is an output of THEIR flow
analysis and is not imported by using their per-method moduli).
"""
import glob, json, math, os, re, sys, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
MODE = sys.argv[1] if len(sys.argv) > 1 else 'gates'
SKY = (MODE == 'sky')
OUTFILE = 'data/stage10t_skyread.txt' if SKY else 'data/stage10t_gates.txt'

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
P(f"STAGE 10T LEG-A REGROW -- mode = {MODE} "
  f"({'sky read' if SKY else 'sky-blind gates'})")
P("")

# ---------------- frozen census (pre-reg R1/R5) ----------------
RECLASS12 = ['ESO563-G021', 'NGC2903', 'NGC4559', 'NGC5585', 'UGC02487',
             'UGC05918', 'UGC05986', 'UGC07323', 'UGC07603', 'UGC08550',
             'UGC09992', 'UGC12632']
ZERO_PT = 'UGC09992'
LT4 = ['CVnIdwA', 'DDO_52', 'DDO_133', 'DDO_210']
UMA_CEPH = {'NGC3972': 21.23, 'NGC4051': 16.62}   # frozen co-read, Mpc

# ---------------- SPARC load (verbatim 10S) ----------------
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

# ---------------- CF4 + cache ----------------
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

# parse checks (G10T-2 part): M31 + the NGC 4258 maser anchor
m31 = cf4[2557]
n4258 = cf4.get(39600, {})
chk_m31 = (m31['ceph'] == (24.397, 0.066) and m31['trgb'] == (24.53, 0.1)
           and m31['sbf'] == (24.214, 0.08))
chk_mas = (n4258.get('mas') is not None
           and 29.3 <= n4258['mas'][0] <= 29.5)
P(f"CF4 parse checks: M31 (PGC 2557) ceph/trgb/sbf = {m31['ceph']}/"
  f"{m31['trgb']}/{m31['sbf']} -> {'PASS' if chk_m31 else 'FAIL'}; "
  f"NGC4258 (PGC 39600) DMmas = {n4258.get('mas')} -> "
  f"{'PASS' if chk_mas else 'FAIL'}")

# ---------------- LITTLE THINGS (verbatim 9R port) ----------------
LT = 'data/littlethings'
def norm9r(s):
    s = s.upper().replace('_', '').replace('-', '').replace(' ', '')
    return re.sub(r'(?<=[A-Z])0+(?=\d)', '', s)
t1, t2 = {}, {}
for line in open(os.path.join(LT, 'table1.dat')):
    f = [x.strip() for x in line.split('|')]
    if len(f) < 9 or not f[0]: continue
    t1[f[0]] = dict(D=float(f[2]), i=float(f[7]), ei=float(f[8]))
for line in open(os.path.join(LT, 'table2.dat')):
    f = [x.strip() for x in line.split('|')]
    if len(f) < 29 or not f[0]: continue
    t2[f[0]] = dict(flag=bool(f[20] or f[23]),
                    Mgas=float(f[24]) if f[24] else None,
                    MstarK=float(f[25]) if f[25] else None,
                    MstarSED=float(f[26]) if f[26] else None)
def read_curve(fn):
    cur = {}
    for line in open(os.path.join(LT, fn)):
        t = line.split()
        if len(t) < 7 or t[1] != 'Data': continue
        R03, V03, Rs, Vs, eVs = map(float, t[2:7])
        cur.setdefault(t[0], []).append((Rs*R03, Vs*V03, eVs*V03))
    return cur
lt_tot = read_curve('rotdmbar.dat')
lt_dmo = read_curve('rotdm.dat')

def lt_port(nm):
    """9R construction: matched-ring V_bar; returns [(R,V,eV,Vb)] at
    the Hunter+12 distance."""
    got = []
    dm_ = lt_dmo.get(nm, [])
    for (R, V, eV) in lt_tot.get(nm, []):
        best, bd = None, 1e9
        for (Rd, Vd, _) in dm_:
            d = abs(Rd - R)
            if d < bd: bd, best = d, Vd
        if best is None or bd > max(0.005, 0.01*R): continue
        vb2 = V*V - best*best
        if vb2 <= 0: continue
        got.append((R, V, eV, math.sqrt(vb2)))
    return got

# ---------------- census (R1-R8) ----------------
P("")
P("-- census (R1-R8; frozen lists asserted) --")
reclass_info = {}
for gi, nm in SP['name'].items():
    if SP['fd'][gi] != 1: continue
    p = cache.get(nm)
    a = adopt(p) if p else None
    if a:
        reclass_info[nm] = (gi, p) + a
rc_names = sorted(reclass_info)
g2_ok = chk_m31 and chk_mas and (rc_names == sorted(RECLASS12))
P(f"R1 reclassification table ({len(rc_names)} on paper):")
for nm in rc_names:
    gi, p, m, D, eD = reclass_info[nm]
    D0, eD0 = SP['dist'][gi]
    npts = int((SP['gal_id'] == gi).sum())
    P(f"  {nm:12s} PGC {p:7d} {m:4s} D {D0:6.2f}->{D:6.2f} "
      f"(s={D/D0:.3f}) e_D {eD0:.2f}->{eD:.2f} Mpc "
      f"({100*eD0/D0:.0f}%->{100*eD/D:.1f}%)  pts={npts}"
      f"{'  [ZERO POINTS - enters no fit]' if npts == 0 else ''}")
zero_ok = int((SP['gal_id'] == reclass_info[ZERO_PT][0]).sum()) == 0

lt_rows = {}
lt_census = []
for nm in sorted(t1):
    key = 'LT:' + nm
    p = cache.get(key)
    ovl = bool(p and p in {v for k, v in cache.items()
                           if v and not k.startswith('LT:')})
    a = adopt(p) if p else None
    fl = t2[nm]['flag']
    iok = t1[nm]['i'] >= 30
    rings = lt_port(nm)
    p10 = sum(1 for (R, V, eV, Vb) in rings if eV/V <= 0.10)
    stat = ('overlap' if ovl else 'flagged' if fl else
            'no-CF4-anchor' if not a else 'i<30' if not iok else
            'thin-points' if p10 < 3 else 'IN')
    if stat == 'IN':
        lt_rows[nm] = (p, a, rings)
    lt_census.append((nm, stat, p10))
lt_in = sorted(lt_rows)
g2_ok &= (lt_in == sorted(LT4)) and zero_ok
P(f"R5 LITTLE THINGS: entering = {lt_in}")
P("   excluded: " + "; ".join(f"{nm}({st})" for nm, st, _ in lt_census
                              if st not in ('IN',)))
for nm in lt_in:
    p, (m, D, eD), rings = lt_rows[nm]
    p10 = sum(1 for r in rings if r[2]/r[1] <= 0.10)
    gd_true = (t2[nm]['Mgas'] or 0) > (t2[nm]['MstarSED']
                                       if t2[nm]['MstarSED'] is not None
                                       else (t2[nm]['MstarK'] or 0))
    P(f"  {nm:10s} PGC {p:7d} {m:4s} D_H12 {t1[nm]['D']:.1f} -> "
      f"D_CF4 {D:.2f} (s={D/t1[nm]['D']:.3f}) e_D {eD:.2f} Mpc "
      f"({100*eD/D:.1f}%)  pts(0.10)={p10}  i={t1[nm]['i']:.0f} "
      f"Mgas>Mstar={'y' if gd_true else 'N'}")
P(f"R4 UMa co-read (frozen): NGC3972 ceph 21.23 Mpc (s=1.180), "
  f"NGC4051 ceph 16.62 Mpc (s=0.923) vs block 18.0 +- 0.9 (depth 2.3)")
P(f"G10T-2 census vs frozen lists: {'PASS' if g2_ok else 'FAIL'}")
if not g2_ok:
    P("STOP: census drift vs pre-registered lists"); save(); sys.exit(0)

# ---------------- world builder ----------------
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

    # LT append
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
                gbar = Vb*Vb/R*KPC          # invariant under s (pin 3)
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

W10S = make_world(reclass=False, lt_slot=None)
W10T = make_world()
assert (len(W10S['legA']), len(W10S['flow'])) == (67, 82)
P("")
P(f"worlds: 10S legA = {len(W10S['legA'])} (regression), 10T legA = "
  f"{len(W10T['legA'])} ({len(W10T['legA'])-len(W10T['uma'])-len(W10T['lt'])-41} "
  f"reclass + 41 anchored + {len(W10T['uma'])} UMa + {len(W10T['lt'])} LT); "
  f"flow {len(W10S['flow'])} -> {len(W10T['flow'])} (leg B not refired)")
assert len(W10T['legA']) == 82 and len(W10T['flow']) == 71

# ---------------- nu families (verbatim) ----------------
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

# ---------------- fit machinery (verbatim 10S bodies) ----------------
def fit_flat_pts(lg, gg, gd, gb_, nu, x0=(-9.92, 1.2)):
    def lo(v):
        la0, fdd = v
        if not (-10.8 < la0 < -9.2) or not (0.3 <= fdd <= 3): return 1e9
        gN = gg + fdd*gd + gb_
        r = lg - np.log10(gN*nu(gN/10**la0))
        return float(np.mean(r*r))
    b = minimize(lo, list(x0), method='Nelder-Mead',
                 options={'xatol': 1e-5, 'fatol': 1e-12, 'maxiter': 3000})
    return b.x[0], b.x[1]

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

def sigma_boot(W, legA, seed=202, reps=200, mask=True):
    """Verbatim 10S sigma_A machinery on world W. Returns (sigma_H0,
    reps_H0). Prints nothing central unless mask=False."""
    rng5 = np.random.default_rng(seed)
    rel_sig, uma_flag = {}, {}
    for g in legA:
        if W['fd'][g] == 4:
            rel_sig[g] = 2.3/18.0; uma_flag[g] = True
        else:
            rel_sig[g] = W['rel'][g]; uma_flag[g] = False
    la0_reps = []
    for _ in range(reps):
        pick = rng5.choice(legA, size=len(legA), replace=True)
        sc_shared = 1.0 + rng5.normal(0, 0.9/18.0)
        rows_, dsh = [], []
        for g in pick:
            s_g = 1.0 + rng5.normal(0, rel_sig[g])
            s_g = min(max(s_g, 0.5), 1.5)
            if uma_flag[g]: s_g *= sc_shared
            rows_.append(W['gpts'][g])
            dsh.append(np.full(len(W['gpts'][g]), -math.log10(s_g)))
        ridx = np.concatenate(rows_)
        lg_b = W['lgobs'][ridx] + np.concatenate(dsh)
        la0_reps.append(fit_flat_pts(lg_b, W['gg'][ridx], W['gd'][ridx],
                                     W['gb'][ridx], nu_be)[0])
    a0_reps = 10**np.array(la0_reps)
    return h0_of_a0(float(np.std(a0_reps))), h0_of_a0(a0_reps)

# ---------------- G10T-1: regression on the 10S world ----------------
P("")
P("-- G10T-1 regression (10S world, archived sky read) --")
ARC = dict(a0=1.0186e-10, f=1.07, s=0.031, u=-0.0258, flat=1.1522e-10,
           p065=1.0969e-10, gm=1.0334e-10, boot=1.0537e-10, sig=12.52,
           gd=1.0056e-10, nd=1.0287e-10)
legA67 = np.concatenate([W10S['ug'][np.array([W10S['fd'][g] in (2, 3, 5)
                                              for g in W10S['ug']])],
                         W10S['uma']])
subA = build_sub(W10S, list(legA67))
bA_plain = fit_hier(subA, nu_be, use_u=False)
bA_tight = fit_hier(subA, nu_be, use_u=True, upri=1e-6,
                    th0=list(bA_plain.x)+[0.0])
d_fun = abs(bA_tight.fun - bA_plain.fun)
d_a0t = abs(10**bA_tight.x[0] - 10**bA_plain.x[0])/10**bA_plain.x[0]
bA67 = fit_hier(subA, nu_be, use_u=True, th0=list(bA_plain.x)+[0.0])
a067 = 10**bA67.x[0]
aidx = np.where(np.isin(W10S['gal_id'], legA67))[0]
la0f, fdf = fit_flat_pts(W10S['lgobs'][aidx], W10S['gg'][aidx],
                         W10S['gd'][aidx], W10S['gb'][aidx], nu_be)
fam67 = {}
for nm_ in ('p065', 'gm', 'boot'):
    bf = fit_hier(subA, FAMS[nm_], use_u=True, th0=list(bA67.x))
    fam67[nm_] = 10**bf.x[0]
sig67, _ = sigma_boot(W10S, legA67)
gdfrac67 = {g: float(np.mean(W10S['gg'][W10S['gpts'][g]] >
                             W10S['gd'][W10S['gpts'][g]] +
                             W10S['gb'][W10S['gpts'][g]]))
            for g in legA67}
gdA67 = [g for g in legA67 if gdfrac67[g] >= 0.5]
ndA67 = [g for g in legA67 if gdfrac67[g] < 0.5]
bgd67 = fit_hier(build_sub(W10S, list(gdA67)), nu_be, use_u=True)
bnd67 = fit_hier(build_sub(W10S, list(ndA67)), nu_be, use_u=True)
checks = [
    ('hier a0', a067, ARC['a0'], 1e-3),
    ('flat a0', 10**la0f, ARC['flat'], 1e-3),
    ('p065', fam67['p065'], ARC['p065'], 1e-3),
    ('gm', fam67['gm'], ARC['gm'], 1e-3),
    ('boot', fam67['boot'], ARC['boot'], 1e-3),
    ('sigma', sig67, ARC['sig'], 5e-3),
    ('GD', 10**bgd67.x[0], ARC['gd'], 5e-3),
    ('nonGD', 10**bnd67.x[0], ARC['nd'], 5e-3),
]
g1_ok = d_fun < 0.5 and d_a0t < 0.002
for lbl, got, want, bar in checks:
    d = abs(got - want)/want
    g1_ok &= d <= bar
    P(f"  {lbl:8s}: {got:.4e} vs archived {want:.4e} "
      f"(d = {100*d:.3f}% vs bar {100*bar:.1f}%)")
P(f"  u-identity: d(-2lnL) = {d_fun:.3f}, d(a0) = {100*d_a0t:.3f}%")
P(f"  params: f_ML = {bA67.x[1]:.2f} (1.07), s_int = {bA67.x[2]:.3f} "
  f"(0.031), u = {bA67.x[3]:+.4f} (-0.0258)")
g1_ok &= (abs(bA67.x[1] - ARC['f']) < 0.02 and
          abs(bA67.x[2] - ARC['s']) < 0.005 and
          abs(bA67.x[3] - ARC['u']) < 0.005)
P(f"G10T-1: {'PASS' if g1_ok else 'FAIL'}  [{time.time()-t00:.0f}s]")
if not g1_ok:
    P("STOP: regression failed; do not fire"); save(); sys.exit(0)

# ---------------- G10T-3: port regression + overlap cross-val --------
P("")
P("-- G10T-3 port (9R identity + overlap cross-validation) --")
ALIAS9R = {
 'WLM': ['UGCA444'], 'DDO50': ['UGC4305', 'HOII'],
 'DDO70': ['UGC5373', 'SEXTANSB'], 'DDO46': ['UGC3966'],
 'DDO52': ['UGC4426'], 'DDO53': ['UGC4459'], 'DDO87': ['UGC5918'],
 'DDO101': ['UGC6900'], 'DDO126': ['UGC7559'], 'DDO133': ['UGC7698'],
 'DDO43': ['UGC3860'], 'HARO29': ['UGCA281'], 'HARO36': ['UGC7950'],
 'NGC3738': ['UGC6565'], 'NGC1569': ['UGC3056'],
 'DDO216': ['PEGDIG', 'UGC12613'], 'DDO69': ['UGC5364', 'LEOA'],
 'DDO75': ['SEXTANSA', 'UGCA205'], 'DDO154': ['NGC4789A'],
 'DDO168': ['UGC8320'], 'DDO155': ['UGC8091'], 'IC10': ['UGC192'],
}
sparc_norm = {norm9r(n) for n in meta}
rows9r = {}
for nm in t1:
    if t2[nm]['flag']: continue
    got = lt_port(nm)
    if got and len(got) >= 0.5*len(lt_tot.get(nm, [])):
        rows9r[norm9r(nm)] = (nm, got)
ltonly9r = []
for key, (nm, got) in rows9r.items():
    cand = [key] + ALIAS9R.get(key, [])
    if not any(c in sparc_norm for c in cand):
        ltonly9r.append(key)
n9r_pts = n9r_gal = 0
for key in ltonly9r:
    nm, got = rows9r[key]
    if t1[nm]['i'] < 30: continue
    pts = sum(1 for (R, V, eV, Vb) in got if eV/V <= 0.10)
    if pts:
        n9r_pts += pts; n9r_gal += 1
g3a_ok = (n9r_pts == 65 and n9r_gal == 7)
P(f"  9R ARM-V census identity: {n9r_pts} pts / {n9r_gal} gal "
  f"(archived 65/7) -> {'PASS' if g3a_ok else 'FAIL'}")

# overlap cross-validation (LT port vs SPARC native, at SPARC distance)
dlg_obs_all, dlg_bar_all = [], []
P(f"  overlap cross-validation (matched rings, LT rescaled to SPARC D):")
for key, (nm, got) in rows9r.items():
    cand = [key] + ALIAS9R.get(key, [])
    hit = [c for c in cand if c in sparc_norm]
    if not hit: continue
    spn = next(n for n in meta if norm9r(n) == hit[0])
    gi = SPNAME2GI.get(spn)
    if gi is None: continue
    D_sp = SP['dist'][gi][0]
    s = D_sp/t1[nm]['D']
    # reread the rotmod for radii (cheap, once per overlap galaxy)
    rr, lo_, lb_ = [], [], []
    for l in open(f"data/sparc/rotmod/{spn}_rotmod.dat"):
        if l.startswith('#'): continue
        tt = l.split()
        if len(tt) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, tt[:6])
        if R <= 0 or Vo <= 0: continue
        gbar = Vg*abs(Vg)/R*KPC + UD*Vd*abs(Vd)/R*KPC + UB*Vb*Vb/R*KPC
        if gbar <= 0: continue
        rr.append(R); lo_.append(math.log10(Vo*Vo/R*KPC))
        lb_.append(math.log10(gbar))
    do_, db_, fb_ = [], [], []
    for (R, V, eV, Vb) in got:
        if eV/V > 0.10: continue
        Rs = R*s
        j = min(range(len(rr)), key=lambda i: abs(rr[i]-Rs)) if rr else None
        if j is None or abs(rr[j]-Rs) > max(0.1, 0.05*Rs): continue
        do_.append(math.log10(V*V/R*KPC) - math.log10(s) - lo_[j])
        db_.append(math.log10(Vb*Vb/R*KPC) - lb_[j])
        fb_.append(Vb*Vb/(V*V))
    if do_:
        dlg_obs_all += do_; dlg_bar_all += db_
        P(f"    {nm:10s} vs {spn:10s}: {len(do_):3d} rings, "
          f"med dlgobs = {np.median(do_):+.3f}, "
          f"med dlgbar = {np.median(db_):+.3f} "
          f"(port Vbar2/Vtot2 med {np.median(fb_):.2f})")
mo, mb = (abs(np.median(dlg_obs_all)), abs(np.median(dlg_bar_all))) \
    if dlg_obs_all else (9, 9)
g3b_ok = mo <= 0.06 and mb <= 0.12
P(f"  pooled ({len(dlg_obs_all)} rings): |med dlgobs| = {mo:.3f} "
  f"(bar 0.06), |med dlgbar| = {mb:.3f} (bar 0.12) -> "
  f"{'PASS' if g3b_ok else 'FAIL (LT excluded from primary)'}")
g3_ok = g3a_ok and g3b_ok
if not g3b_ok:
    W10T = make_world(lt_slot=None)   # LT out of primary per the clause
    P("  LT members EXCLUDED from primary per the pre-registered clause")
    P(f"  primary leg A = {len(W10T['legA'])} (41 anchored + 26 UMa + "
      f"11 reclassified); LT = co-read only")

# ---------------- G10T-4: injections (truths beyond every edge) ------
P("")
P("-- G10T-4 injections (extended design) --")
legA82 = W10T['legA']
mock_full = np.empty_like(W10T['lgobs'])
gN1 = W10T['gg'] + 1.0*W10T['gd'] + W10T['gb']
g4_ok = True
for h0t in (50.0, 62.0, 73.0, 85.0, 100.0):
    a0t = a0_of_h0(h0t)
    mock_full[:] = np.log10(gN1*nu_be(gN1/a0t))
    bm = fit_hier(build_sub(W10T, list(legA82), lg_src=mock_full.copy()),
                  nu_be, use_u=True)
    h0r = h0_of_a0(10**bm.x[0])
    derr = abs(h0r - h0t)/h0t
    g4_ok &= derr <= 0.005
    P(f"  noiseless {h0t:5.0f} -> {h0r:7.2f} ({100*derr:.3f}%)")
P(f"  noiseless -> {'PASS' if g4_ok else 'FAIL'}  [{time.time()-t00:.0f}s]")
for h0t in (62.0, 85.0):
    errs = []
    a0t = a0_of_h0(h0t)
    base = np.log10(gN1*nu_be(gN1/a0t))
    for seed in (11, 23, 42, 101, 202):
        rng = np.random.default_rng(seed)
        mock_full[:] = base + rng.normal(0, np.sqrt(W10T['sig2'] + 0.08**2))
        voff = [rng.normal(0, W10T['sigv'][g]) for g in legA82]
        bm = fit_hier(build_sub(W10T, list(legA82),
                                lg_src=mock_full.copy(), dlg_per_occ=voff),
                      nu_be, use_u=True)
        errs.append(100*(h0_of_a0(10**bm.x[0]) - h0t)/h0t)
    errs = np.array(errs)
    se = errs.std(ddof=1)/math.sqrt(len(errs))
    bias_ok = abs(errs.mean()) <= max(1.0, 2*se)
    g4_ok &= bias_ok
    P(f"  noisy truth {h0t:.0f}: errors % = "
      f"{[('%+.2f' % e) for e in errs]}; mean {errs.mean():+.2f} "
      f"(bar max(1%, 2SE={2*se:.2f}%)), SD {errs.std(ddof=1):.2f}% "
      f"= realization floor -> {'PASS' if bias_ok else 'FAIL'}")
P(f"G10T-4: {'PASS' if g4_ok else 'FAIL'}  [{time.time()-t00:.0f}s]")
if not g4_ok:
    P("STOP: injection gate failed"); save(); sys.exit(0)

# ---------------- G10T-5: masked power ----------------
P("")
sigT, repsT = sigma_boot(W10T, legA82)
P(f"-- G10T-5 power (central masked in gates mode) --")
P(f"  sigma(H0_A, extended) = {sigT:.2f} km/s/Mpc "
  f"(10S: 12.52; bar for T-SHARPENED path: <= 12.0)")
g5_open = sigT <= 12.0
P(f"G10T-5: {'T-SHARPENED path OPEN' if g5_open else 'T-NULL-GROWTH'}")

# ---------------- G10T-6: provenance ----------------
P("")
P("-- G10T-6 provenance --")
P("  CF4 ReadMe Note 4 (verbatim): per-method moduli are given 'after")
P("  registration to a common scale with the MCMC analysis'; absolute")
P("  scale = 'Cepheid period-luminosity relation and tip of the red")
P("  giant branch observations founded on local stellar parallax")
P("  measurements along with the geometric maser distance to NGC 4258'.")
P("  TF/FP columns EXCLUDED (dynamical circularity); combined DM")
P("  excluded (mixes TF/FP). CF4's own H0 = 74.6 NOT imported.")
P("  Lit: three NOT-FOUNDs on a0->H0 inversion (2026-09-20, docstring).")

gates_all = g1_ok and g2_ok and g3a_ok and g4_ok
P("")
P(f"GATES: G1 {'PASS' if g1_ok else 'FAIL'}  G2 {'PASS' if g2_ok else 'FAIL'}  "
  f"G3a {'PASS' if g3a_ok else 'FAIL'}  G3b {'PASS' if g3b_ok else 'LT-EXCLUDED'}  "
  f"G4 {'PASS' if g4_ok else 'FAIL'}  G5 {'open' if g5_open else 'null-growth'}")
if not SKY:
    P("")
    P("gates mode complete -- no extended-set central computed anywhere.")
    P(f"wall-clock: {(time.time()-t00)/60:.1f} min")
    save()
    print("\nsaved:", OUTFILE)
    sys.exit(0)

# ================= SKY (extended leg A) =================
P("")
P("== 10T SKY READ: the extended leg A ==")
subT = build_sub(W10T, list(legA82))
bT_plain = fit_hier(subT, nu_be, use_u=False)
bT = fit_hier(subT, nu_be, use_u=True, th0=list(bT_plain.x)+[0.0])
a0_T = 10**bT.x[0]
H0_T = h0_of_a0(a0_T)
P(f"PRIMARY hier BE: a0 = {a0_T:.4e}, f_ML = {bT.x[1]:.2f}, s_int = "
  f"{bT.x[2]:.3f}, u = {bT.x[3]:+.4f} dex -> H0_A = {H0_T:.2f}")
tidx = np.where(np.isin(W10T['gal_id'], legA82))[0]
la0Tf, fdTf = fit_flat_pts(W10T['lgobs'][tidx], W10T['gg'][tidx],
                           W10T['gd'][tidx], W10T['gb'][tidx], nu_be)
P(f"flat co-read: a0 = {10**la0Tf:.4e} (H0 {h0_of_a0(10**la0Tf):.1f}), "
  f"f_d = {fdTf:.2f}")
famT = {'BE': a0_T}
for nm_ in ('p065', 'gm', 'boot'):
    bf = fit_hier(subT, FAMS[nm_], use_u=True, th0=list(bT.x))
    famT[nm_] = 10**bf.x[0]
P("hier family: " + ", ".join(f"{k} {v:.4e} (H0 {h0_of_a0(v):.1f})"
                              for k, v in famT.items()))
P(f"sigma_A (flat boot, 200 reps, seed 202): {sigT:.2f}; percentiles "
  f"16/50/84 = {np.percentile(repsT, [16, 50, 84]).round(1).tolist()}")
P(f"LEG A (extended): H0_A = {H0_T:.2f} +- {sigT:.2f} km/s/Mpc "
  f"(hier BE primary; family H0 "
  f"{h0_of_a0(min(famT.values())):.1f}-{h0_of_a0(max(famT.values())):.1f})")
P("")

# variants
P("-- variants --")
Wv2 = make_world(uma_moved=True)
bv2 = fit_hier(build_sub(Wv2, list(Wv2['legA'])), nu_be, use_u=True,
               th0=list(bT.x))
P(f"A-var2 (UMa Cepheid pair moved, u-group 24): a0 = "
  f"{10**bv2.x[0]:.4e} (H0 {h0_of_a0(10**bv2.x[0]):.1f}), "
  f"u = {bv2.x[3]:+.4f}")
if not g3b_ok:
    Wlt = make_world()   # gas slot, LT included (the would-be primary)
    blt = fit_hier(build_sub(Wlt, list(Wlt['legA'])), nu_be, use_u=True,
                   th0=list(bT.x))
    P(f"LT-inclusive co-read (gas slot; EXCLUDED from primary by "
      f"G10T-3b): a0 = {10**blt.x[0]:.4e} "
      f"(H0 {h0_of_a0(10**blt.x[0]):.1f})")
Wv3 = make_world(lt_slot='disk')
bv3 = fit_hier(build_sub(Wv3, list(Wv3['legA'])), nu_be, use_u=True,
               th0=list(bT.x))
P(f"A-var3 (LT disk-slot, f+dml on total{'; co-read only' if not g3b_ok else ''}): "
  f"a0 = {10**bv3.x[0]:.4e} (H0 {h0_of_a0(10**bv3.x[0]):.1f})")
Wv4 = make_world(uniform=True)
nswap = sum(1 for gi, nm in SP['name'].items()
            if W10T['fd'].get(gi) in (2, 3, 5) and nm not in reclass_info
            and cache.get(nm) and adopt(cache[nm]))
bv4 = fit_hier(build_sub(Wv4, list(Wv4['legA'])), nu_be, use_u=True,
               th0=list(bT.x))
P(f"A-var4 (CF4-uniform; {nswap} original anchors swapped + UMa pair): "
  f"a0 = {10**bv4.x[0]:.4e} (H0 {h0_of_a0(10**bv4.x[0]):.1f})")
Wv5 = make_world(lt_cut=0.20)
bv5 = fit_hier(build_sub(Wv5, list(Wv5['legA'])), nu_be, use_u=True,
               th0=list(bT.x))
nlt5 = int(np.isin(Wv5['gal_id'], Wv5['lt']).sum())
P(f"A-var5 (LT cut 0.20, {len(Wv5['lt'])} LT gal / {nlt5} pts"
  f"{'; co-read only' if not g3b_ok else ''}): a0 = "
  f"{10**bv5.x[0]:.4e} (H0 {h0_of_a0(10**bv5.x[0]):.1f})")
P(f"10S-67 baseline (G10T-1): a0 = {a067:.4e} (H0 {h0_of_a0(a067):.2f} "
  f"+- {sig67:.2f})")
varH = [h0_of_a0(10**b.x[0]) for b in ((bv2, bv4) if not g3b_ok
                                       else (bv2, bv3, bv4, bv5))] + [H0_T]
P(f"variant spread (primary-eligible variants): "
  f"{min(varH):.1f} - {max(varH):.1f}")
P("")

# G10T-7 composition splits
P("-- G10T-7 composition --")
gdfracT = {g: float(np.mean(W10T['gg'][W10T['gpts'][g]] >
                            W10T['gd'][W10T['gpts'][g]] +
                            W10T['gb'][W10T['gpts'][g]]))
           for g in legA82}
gdT = [g for g in legA82 if gdfracT[g] >= 0.5]
ndT = [g for g in legA82 if gdfracT[g] < 0.5]
bgdT = fit_hier(build_sub(W10T, list(gdT)), nu_be, use_u=True)
bndT = fit_hier(build_sub(W10T, list(ndT)), nu_be, use_u=True)
sig_a0T = a0_of_h0(sigT)
d_gd = abs(10**bgdT.x[0] - 10**bndT.x[0])
r_gd = d_gd/sig_a0T
P(f"GD split (hier): GD({len(gdT)}) {10**bgdT.x[0]:.4e} / "
  f"non-GD({len(ndT)}) {10**bndT.x[0]:.4e}; |d| = {d_gd:.2e} = "
  f"{r_gd:.2f} sigma (LT members GD by construction of the port; all 4 "
  f"genuinely Mgas > Mstar per Oh+15)")
new_g = [g for g in legA82
         if (W10T['fd'][g] == 9) or
            (W10T['names'][g] in reclass_info)]
old_g = [g for g in legA82 if g not in set(new_g)]
bnew = fit_hier(build_sub(W10T, list(new_g)), nu_be, use_u=False)
sig_new, _ = sigma_boot(W10T, np.array(new_g))
a0_new = 10**bnew.x[0]
d_on = abs(a067 - a0_new)
sig_on = math.hypot(a0_of_h0(sig67), a0_of_h0(sig_new))
r_on = d_on/sig_on
P(f"old-vs-new split: old(67) {a067:.4e} (H0 {h0_of_a0(a067):.1f}) vs "
  f"new({len(new_g)}) {a0_new:.4e} (H0 {h0_of_a0(a0_new):.1f} +- "
  f"{sig_new:.1f}); |d| = {r_on:.2f} x joint sigma")
split_max = max(r_gd, r_on)
P("")

# G10T-8 UMa Cepheid co-read
P("-- G10T-8 UMa co-read (annotation) --")
u_hat = bT.x[3]
P(f"fitted shared u = {u_hat:+.4f} dex (distance factor "
  f"{10**(-u_hat):.3f} on 18.0 Mpc -> {18.0*10**(-u_hat):.2f} Mpc); "
  f"frozen Cepheids: NGC3972 21.23 Mpc (+0.072 dex), NGC4051 16.62 Mpc "
  f"(-0.035 dex); N=2 inside a 2.3 Mpc depth -- annotation grade")
P("")

# letter
P("== LETTER ==")
if split_max > 2.0:
    letter = 'T-COMPOSITION-SPLIT'
    P("T-COMPOSITION-SPLIT: a split exceeds 2 sigma; both subset numbers "
      "stand, NO combined headline; the 10S-67 read stays operative. "
      "Suspects: GD composition (dissident-dwarf interaction), the LT "
      "port, the CF4-vs-SPARC zero-point seam.")
elif not g5_open:
    letter = 'T-NULL-GROWTH'
    P("T-NULL-GROWTH: sigma did not improve; the 10S read stands; growth "
      "re-routes to BIG-SPARC/WALLABY triggers.")
else:
    letter = 'T-SHARPENED'
    caveat = " (with the 1-2 sigma composition caveat)" \
        if split_max > 1.0 else ""
    if not g3b_ok:
        caveat += " [port caveat: LT excluded from primary by G10T-3b]"
    P(f"T-SHARPENED{caveat}: the regrown anchored leg reads")
    P(f"  H0(leg A) = {H0_T:.1f} +- {sigT:.1f} (stat) km/s/Mpc,")
    P(f"  function-family band {h0_of_a0(min(famT.values())):.1f}-"
      f"{h0_of_a0(max(famT.values())):.1f}, variant band {min(varH):.1f}-"
      f"{max(varH):.1f} (H0-assumption-independent; hier BE primary).")
    P("  This is a leg-A update of the 10S M-GRAY state: the meter still "
      "returns no joint H0 (leg B power-limited at SPARC grade, R45).")
P("")
P("MANDATORY DISCLOSURES:")
P(f"  treatment split: flat co-read a0 = {10**la0Tf:.3e} vs hier "
  f"{a0_T:.3e} (the 4H/5M M/L-vs-a0 degeneracy; hier = pre-registered "
  f"primary).")
P(f"  composition: additions are GD-heavy by design; GD split "
  f"{r_gd:.2f} sigma, old-vs-new {r_on:.2f} sigma.")
P(f"  flow census: 82 -> 71 with-points (11 reclassified; UGC09992 "
  f"reclassified on paper, zero surviving points); leg B NOT refired "
  f"(R45 standing).")
if not g3b_ok:
    P(f"  LT port: velocity fields agree (pooled med dlgobs 0.008 dex) "
      f"but the baryon models split (med dlgbar -0.16 dex, bar 0.12; "
      f"port Vbar2/Vtot2 ~ 0.1-0.2 = DM-dominated subtraction) -> the 4 "
      f"LT dwarfs are co-reads, never primary members, per the "
      f"pre-registered G10T-3b clause.")
P(f"  three-world read (PREDICTIONS SF, report grade): at sigma ~ "
  f"{sigT:.0f} the table separates nothing; registered for the "
  f"BIG-SPARC-era refire.")
P("")
P(f"verdict letter: {letter}")
P("credences: NO cell moves on any outcome (pre-registered); 53/8 "
  "untouched")
P(f"total wall-clock: {(time.time()-t00)/60:.1f} min")
save()
print("\nsaved:", OUTFILE)
