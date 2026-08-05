"""STAGE 9C — THE RADIUS-DEPENDENT AD CORRECTION (GD).
Pre-registered BEFORE any run.

The 8Y flat-sigma lever lacked radius dependence; the real
asymmetric drift grows outward.  Bounding model (crudeness
pre-stated: outer-disk exponential, HI central holes ignored):
Sigma_HI(R) = S0 exp(-R/Rg), Rg solved per galaxy from the catalog
pair (MHI, RHI) via Sigma(RHI) = 1 Msun/pc^2; correction g_obs ->
(V^2 + sigma^2 * R/Rg)/R (normalization-free slope, constant
sigma).  lam_GD(sigma) at sigma in {0, 8, 10, 12, 15} km/s; DD
control at sigma = 10; 100-rep GD bootstrap at sigma = 10.

Gates: G9C-0 sigma=0 bit-identity to LGOBS0 + GD fit = 8X (0.002);
G9C-1 engine probes <= 1e-6; G9C-2 Rg accounting.
Bars (at sigma = 10): C1 NEUTRALIZED iff lam_GD >= 0.0.  C2
PARTIAL iff -1.0 < lam_GD < 0.0.  C3 OUT-OF-REACH iff lam_GD <=
-1.0 (the pressure story dead at catalog grade).  DD bluntness
flag if |d lam_DD| > 0.3.  Co-read: sigma* crossings.
NO credence movement (measurement round; pre-stated).
Output: data/stage9c_adcorr.txt
"""
import glob, math, os, re, time
import numpy as np
from scipy.optimize import minimize, brentq

KPC = 3.24078e-14
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)

meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines) if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        name = t[0]
        D, eD = float(t[2]), float(t[3])
        inc, einc = float(t[5]), float(t[6])
        mhi, rhi = float(t[13]), float(t[14])
        q = int(t[17])
        meta[name] = (inc, q, D, eD, einc, mhi, rhi)
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
vobs2, rkpc = [], []
svert = {}
MHIRHI = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc, mhi, rhi = meta.get(
        name, (0, 3, 1.0, 0.0, 0.0, 0.0, 0.0))
    if inc < 30 or q > 2: continue
    s_d2 = (eD/max(D, 1e-6))/LN10
    s_i2 = 2.0*math.radians(einc)/max(math.tan(math.radians(inc)), 1e-6)/LN10
    svert[gi] = math.sqrt(s_d2**2 + s_i2**2)
    MHIRHI[gi] = (mhi, rhi)
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
        vobs2.append(Vo*Vo); rkpc.append(R)
        sig.append(2*eV/Vo/math.log(10))
        gal_id.append(gi)
g_gas, g_dsk, g_bul, gobs, sig, gal_id, vobs2, rkpc = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id, vobs2, rkpc))
sig2 = sig*sig
LGOBS0 = np.log10(gobs)

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_lam(y, lam):
    return (1.0-lam)*nu_standard(y) + lam*nu_be(y)

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

t00 = time.time()
P("9C THE RADIUS-DEPENDENT AD CORRECTION (pre-reg committed BEFORE "
  "any run; measurement round; NO credence movement)")

allg = np.unique(gal_id)
NGAL = gal_id.max()+1
gdfrac = np.zeros(NGAL)
for g_ in allg:
    m = gal_id == g_
    gdfrac[g_] = float(np.mean(g_gas[m] > g_dsk[m] + g_bul[m]))
gd_set = np.array([g_ for g_ in allg if gdfrac[g_] >= 0.5])
dd_set = np.array([g_ for g_ in allg if gdfrac[g_] < 0.5])
SV = np.array([svert.get(g_, 0.0) for g_ in range(NGAL)])
GIDX = [np.where(gal_id == g_)[0] for g_ in range(NGAL)]
GD_LIST = [int(g_) for g_ in gd_set]
DD_LIST = [int(g_) for g_ in dd_set]

# ---- Rg per galaxy: Sigma_HI(RHI) = 1 with Sigma0 = MHI/(2 pi Rg^2)
RG = np.zeros(NGAL)
n_solved, n_fb = 0, 0
for g_ in allg:
    mhi, rhi = MHIRHI.get(int(g_), (0.0, 0.0))
    rg = None
    if mhi > 0 and rhi > 0:
        def gfun(rg_):
            s0 = mhi*1e3/(2*np.pi*rg_**2)   # Msun/pc^2, rg in kpc
            if s0 <= 0: return 1e9
            return rhi/rg_ - math.log(s0)
        try:
            lo_, hi_ = rhi/20.0, rhi*3.0
            if gfun(lo_)*gfun(hi_) < 0:
                rg = brentq(gfun, lo_, hi_, xtol=1e-6)
                n_solved += 1
        except Exception:
            rg = None
    if rg is None:
        rmax = float(rkpc[GIDX[g_]].max()) if len(GIDX[g_]) else 1.0
        rg = 0.5*rmax
        n_fb += 1
    RG[g_] = rg
rg_pts = RG[gal_id]
corr_frac10 = (100.0*rkpc/rg_pts)/vobs2   # sigma=10 fraction of V^2
gd_pts = np.isin(gal_id, gd_set)
P(f"G9C-2 Rg accounting: solved {n_solved}, fallback {n_fb}; "
  f"Rg pct 10/50/90 = "
  + "/".join(f"{v:.1f}" for v in np.percentile(RG[allg], [10, 50, 90]))
  + " kpc; GD correction fraction (sigma=10) pct 50/90/max = "
  + "/".join(f"{v:.2f}" for v in
             np.percentile(corr_frac10[gd_pts], [50, 90, 100])))

txt8x = open('data/stage8x_regime.txt').read()
TGT_GD = float(re.search(r"\[GD-all\s*\] lam_hat = ([+-][\d.]+)",
                         txt8x).group(1))

def lg_of_sig(sg):
    return np.log10((vobs2 + sg*sg*rkpc/rg_pts)/rkpc*KPC)

def m2ll_slow(th, lam, gset, lg):
    la0, f, s_int = th
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    nv = nu_lam(gN/a0, lam)
    if np.min(nv) <= 0.05: return 1e12
    gm = gN*nv
    v = sig2 + s_int*s_int
    a = lg - np.log10(gm)
    out = 0.0
    for g_ in gset:
        ji = GIDX[g_]
        if len(ji) == 0: continue
        aj, vj = a[ji], v[ji]
        s_ = SV[g_]
        if s_ <= 1e-4:
            out += float(np.sum(aj*aj/vj + np.log(vj)))
        else:
            iv = 1.0/vj
            Siv = float(np.sum(iv))
            Sa = float(np.sum(aj*iv))
            out += (float(np.sum(aj*aj*iv))
                    - Sa*Sa/(Siv + 1.0/(s_*s_))
                    + float(np.sum(np.log(vj)))
                    + math.log(1.0 + s_*s_*Siv))
    return out

def make_instance(gset_list):
    cats, labs, svals = [], [], []
    li = 0
    for g_ in gset_list:
        ji = GIDX[g_]
        if len(ji) == 0: continue
        cats.append(ji)
        labs.append(np.full(len(ji), li, dtype=np.int64))
        svals.append(SV[g_])
        li += 1
    if li == 0:
        return None
    return (np.concatenate(cats), np.concatenate(labs),
            np.array(svals), li)

def m2ll_fast(th, lam, cat, lab, svec, ninst, lg):
    la0, f, s_int = th
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    nv = nu_lam(gN/a0, lam)
    if np.min(nv) <= 0.05: return 1e12
    gm = gN*nv
    v = sig2 + s_int*s_int
    a = lg - np.log10(gm)
    ac, vc_ = a[cat], v[cat]
    iv = 1.0/vc_
    Siv = np.bincount(lab, weights=iv, minlength=ninst)
    Sa  = np.bincount(lab, weights=ac*iv, minlength=ninst)
    Saa = np.bincount(lab, weights=ac*ac*iv, minlength=ninst)
    Slv = np.bincount(lab, weights=np.log(vc_), minlength=ninst)
    s2 = svec*svec
    on = svec > 1e-4
    out = np.where(on,
                   Saa - Sa*Sa/(Siv + 1.0/np.maximum(s2, 1e-30))
                   + Slv + np.log(1.0 + s2*Siv),
                   Saa + Slv)
    return float(np.sum(out))

LGB = np.round(np.arange(-2.0, 1.501, 0.25), 3)
def lam_hat_fast(inst, lg):
    cat, lab, svec, ninst = inst
    prof = []; th = None
    for lam in LGB:
        starts = (([list(th)] if th is not None else [])
                  + [[math.log10(A0_FID), 1.0, 0.08]])
        best = None
        for th0 in starts:
            b = minimize(lambda t: m2ll_fast(t, lam, cat, lab, svec,
                                             ninst, lg), th0,
                         method='Nelder-Mead',
                         options=dict(maxiter=4000, xatol=1e-6,
                                      fatol=1e-7))
            if best is None or b.fun < best.fun: best = b
        prof.append(best.fun); th = best.x
    prof = np.array(prof)
    i = int(np.argmin(prof))
    lam_hat = float(LGB[i])
    lo_edge = (i == 0); hi_edge = (i == len(LGB)-1)
    if not (lo_edge or hi_edge):
        x3, y3 = LGB[i-1:i+2], prof[i-1:i+2]
        c2_, c1_, _ = np.polyfit(x3, y3, 2)
        if c2_ > 0:
            lh = -c1_/(2*c2_)
            if LGB[i-1] <= lh <= LGB[i+1]: lam_hat = float(lh)
    return lam_hat, lo_edge, hi_edge

# ---------------- gates ----------------
inst_gd = make_instance(GD_LIST)
lg0 = lg_of_sig(0.0)
d0 = float(np.max(np.abs(lg0 - LGOBS0)))
lh0, _, _ = lam_hat_fast(inst_gd, LGOBS0)
ok0 = (d0 <= 1e-12) and (abs(lh0 - TGT_GD) <= 0.002)
P(f"G9C-0 sigma=0: lg identity max|d| = {d0:.2e}; GD lam_hat = "
  f"{lh0:+.3f} vs 8X {TGT_GD:+.3f} -> {'PASS' if ok0 else 'FAIL'}")
rngp = np.random.default_rng(7)
dmx = 0.0
psets = [GD_LIST, DD_LIST, GD_LIST]
for pi in range(12):
    thp = [-10.5+rngp.random(), 0.4+2.0*rngp.random(),
           0.01+0.29*rngp.random()]
    lam = -2.0+3.5*rngp.random()
    sg = [0.0, 10.0, 15.0][pi % 3]
    gl = psets[pi % 3]
    lg = lg_of_sig(sg)
    slow = m2ll_slow(thp, lam, gl, lg)
    inst = make_instance(gl)
    fast = m2ll_fast(thp, lam, *inst, lg)
    if slow < 1e11:
        dmx = max(dmx, abs(slow-fast))
ok1 = dmx <= 1e-6
P(f"G9C-1 fast-vs-slow (12 probes, sigma 0/10/15): max|d| = "
  f"{dmx:.2e} -> {'PASS' if ok1 else 'FAIL'}")
if not (ok0 and ok1):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    with open('data/stage9c_adcorr.txt', 'w') as f:
        f.write("\n".join(L)+"\n")
    raise SystemExit(0)
P("GATES: G9C-0/1/2 ALL PASS")
P("")

# ---------------- the sigma curve ----------------
SGRID = [0.0, 8.0, 10.0, 12.0, 15.0]
scur = []
for sg in SGRID:
    lh, lo_e, hi_e = lam_hat_fast(inst_gd, lg_of_sig(sg))
    scur.append(lh)
    P(f"[GD sigma-curve] sigma={sg:4.1f}: lam_hat = {lh:+.3f}"
      + (' LO-EDGE' if lo_e else (' HI-EDGE' if hi_e else '')))
inst_dd = make_instance(DD_LIST)
lh_dd0, _, _ = lam_hat_fast(inst_dd, lg_of_sig(0.0))
lh_dd10, _, _ = lam_hat_fast(inst_dd, lg_of_sig(10.0))
P(f"[DD control] sigma=0: {lh_dd0:+.3f}; sigma=10: {lh_dd10:+.3f} "
  f"(d = {lh_dd10-lh_dd0:+.3f})")
scur = np.array(scur)
for thr_, nm_ in ((-1.0, 'lam=-1'), (-0.5, 'lam=-0.5'),
                  (0.0, 'lam=0')):
    if scur[0] >= thr_:
        P(f"sigma*({nm_}): already at sigma=0")
    elif np.all(scur < thr_):
        P(f"sigma*({nm_}): NOT REACHED by sigma=15")
    else:
        j = int(np.argmax(scur >= thr_))
        ss = np.interp(thr_, [scur[j-1], scur[j]],
                       [SGRID[j-1], SGRID[j]])
        P(f"sigma*({nm_}) = {ss:.1f} km/s")
P("")

# ---------------- 100-rep GD bootstrap at sigma=10 ----------------
NB = 100
rng = np.random.default_rng(71)
lg10 = lg_of_sig(10.0)
bv = []
for r in range(NB):
    draw = rng.integers(0, len(GD_LIST), size=len(GD_LIST))
    rep = [GD_LIST[j] for j in draw]
    lh, _, _ = lam_hat_fast(make_instance(rep), lg10)
    bv.append(lh)
bv = np.array(bv)
qb = np.percentile(bv, [5, 50, 95])
P(f"[GD sigma=10 bootstrap x{NB}] pct 5/50/95 = "
  f"{qb[0]:+.3f}/{qb[1]:+.3f}/{qb[2]:+.3f}")
P("")
lam10 = scur[SGRID.index(10.0)]
blunt = abs(lh_dd10-lh_dd0) > 0.3
bf = "  [DD-BLUNT-FLAG]" if blunt else ""
if lam10 >= 0.0:
    P("==> 9C VERDICT (locked grammar): NEUTRALIZED - the "
      "radius-dependent AD correction at physical sigma brings the "
      "GD dial to zero: the dial is pressure-support at catalog "
      "grade; T5's defense collapses to conditional" + bf + ".")
elif lam10 > -1.0:
    P("==> 9C VERDICT (locked grammar): PARTIAL - the "
      "radius-dependent correction moves the dial materially but "
      "does not neutralize it at physical sigma" + bf + ".")
else:
    P("==> 9C VERDICT (locked grammar): OUT-OF-REACH - the "
      "radius-dependent version also fails at catalog grade; the "
      "pressure story is dead here; remaining = external "
      "Sigma_gas(R) data or genuine physics" + bf + ".")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t00)/60:.1f} min)")

with open('data/stage9c_adcorr.txt', 'w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage9c_adcorr.txt")
