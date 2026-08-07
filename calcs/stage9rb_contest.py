"""
STAGE 9R-b -- THE TWO-WORLD CONTEST ON LITTLE THINGS (the 9R successor).

9R's profile estimator failed its injection power gates (GD world
-1.31 recovered at ~-0.5; relaxed arm edge-rides), so the sky was
never read. But the injected WORLDS still separated (~1.1 on the
lam-hat scale). This stage replaces estimation with a CONTEST: which
of the two pre-specified worlds fits the independent dwarfs better?

  D = m2ll(lam = +0.97) - m2ll(lam = -1.31),

each side minimized over the verbatim nuisances (log10 a0, f, s_int),
vertON, LT-ONLY set. D > 0 means the GD-dial world fits better. The
contest does NOT estimate lam; it ranks two fixed hypotheses. A sky
that is neither world lands GRAY by construction (stated limitation:
the calibration is same-model; misspecification widens GRAY, it does
not manufacture a side).

CALIBRATION (the power gate, now the measurement of the statistic's
own distributions): inject 20 mock skies per world per arm (the 9R
mock construction verbatim: log-space point noise sqrt(sig^2+0.08^2)
+ per-galaxy vertical draws; truth nuisances f = 1, a0 = A0_FID,
s_int = 0.08; seed scheme: headline truths default_rng(1000+k), GD
truths default_rng(1500+k), k = 0..19). Threshold t* = midpoint of
the two calibration means. An arm is POWERED iff each world lands on
its own side of t* in >= 16/20 injections. The verdict arm = the
powered arm with the lower total misclassification (tie -> ARM-V).

PRE-REGISTERED BARS (locked; read on the verdict arm's sky D):
  B-GD-SIDE       : D_sky > t* AND D_sky > max(all 20 headline-truth
                    calibration draws). The independent dwarfs prefer
                    the GD-dial world over the headline at
                    beyond-calibration-range grade (~p < 0.05
                    empirical).
  B-HEADLINE-SIDE : D_sky < t* AND D_sky < min(all 20 GD-truth
                    draws). The independent dwarfs prefer the
                    headline world (the dissolution direction).
  B-GRAY          : anything else (D_sky quoted against both
                    calibration distributions; no side claimed).
  B-POWER-LIMITED : no powered arm -> TODO 28 closes to "awaits
                    better input" (Iorio+2017 curves / DR4-era) at
                    the public-Oh+15 grade.
CO-READS (never verdict-bearing): the contest D on LT-ALL and on the
SPARC-overlap subset (verdict arm's cut); both calibration
distributions printed in full.

CREDENCE: NO movement any branch (measurement round, pre-stated);
paper consequence booked at the next author-called thaw.

DISCLOSURE: nothing about the contest has been run before this
commit. Motivation disclosed: the 9R ARM-V world separations (~1.1)
suggested a ranking statistic might be powered where the estimator
is not -- that is a hope, not a result; the calibration decides.

GATES:
  G9Rb-M  manifest: the four LT files match the committed SHA256s.
  G9Rb-1  census regression: the arm censuses reproduce 9R's printed
          values exactly (ARM-V 65 pts / 7 gal; ARM-R 172 / 12), and
          the LT-only membership (14) matches.
  G9Rb-0  engine regression: re-running the 9R power-gate cell
          (truth -1.31, seed 31, ARM-V) through this script's copied
          engine reproduces the printed lam_hat = -0.524 within 0.02
          (ties the two stages' engines + rng streams exactly).
  G9Rb-L  ledger leg: gal-9r-ltrepl status CURRENT.

Output: data/stage9rb_contest.txt. Wall-clock: ~3-6 min CPU.
"""
import csv
import glob
import hashlib
import math
import os
import re
import time

import numpy as np
from scipy.optimize import minimize

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

KPC = 3.24078e-14
A0_FID = 1.2e-10
LN10 = math.log(10)
LT = 'data/littlethings'

OUT = []
def P(s=""):
    print(s, flush=True)
    OUT.append(s)

def save():
    with open('data/stage9rb_contest.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

t00 = time.time()
P("=" * 78)
P("STAGE 9R-b -- THE TWO-WORLD CONTEST (headline +0.97 vs GD -1.31)")
P("=" * 78)

gates = {}

# ---------------- G9Rb-M: manifest ----------------------------------
man = {}
for line in open('data/littlethings_manifest.sha256'):
    h, n = line.split()
    man[n.lstrip('*')] = h
ok_m = True
for fn in ['table1.dat', 'table2.dat', 'rotdmbar.dat', 'rotdm.dat']:
    h = hashlib.sha256(open(os.path.join(LT, fn), 'rb').read()).hexdigest()
    ok_m &= (man.get(fn) == h)
gates['G9Rb-M'] = ok_m
P("G9Rb-M manifest: %s" % ("PASS" if ok_m else "FAIL"))
if not ok_m:
    save(); raise SystemExit(0)

# ---------------- LT construction (verbatim 9R) ---------------------
def norm(s):
    s = s.upper().replace('_', '').replace('-', '').replace(' ', '')
    return re.sub(r'(?<=[A-Z])0+(?=\d)', '', s)

t1 = {}
for line in open(os.path.join(LT, 'table1.dat')):
    f = [x.strip() for x in line.split('|')]
    if len(f) < 9 or not f[0]: continue
    t1[norm(f[0])] = dict(name=f[0], D=float(f[2]), i=float(f[7]),
                          ei=float(f[8]))
t2 = {}
for line in open(os.path.join(LT, 'table2.dat')):
    f = [x.strip() for x in line.split('|')]
    if len(f) < 29 or not f[0]: continue
    mg = float(f[24]) if f[24] else None
    mk = float(f[25]) if f[25] else None
    ms = float(f[26]) if f[26] else None
    t2[norm(f[0])] = dict(flag=bool(f[20] or f[23]), Mgas=mg,
                          MstarK=mk, MstarSED=ms)

def read_curve(fn):
    cur = {}
    for line in open(os.path.join(LT, fn)):
        t = line.split()
        if len(t) < 7 or t[1] != 'Data': continue
        nm = norm(t[0])
        R03, V03, Rs, Vs, eVs = map(float, t[2:7])
        cur.setdefault(nm, []).append((Rs*R03, Vs*V03, eVs*V03))
    return cur

tot = read_curve('rotdmbar.dat')
dmo = read_curve('rotdm.dat')
rows = {}
for nm, rings in tot.items():
    if nm not in t2 or t2[nm]['flag'] or nm not in t1: continue
    dm = dmo.get(nm, [])
    got = []
    for (R, V, eV) in rings:
        best, bd = None, 1e9
        for (Rd, Vd, _) in dm:
            d = abs(Rd - R)
            if d < bd: bd, best = d, (Rd, Vd)
        if best is None or bd > max(0.005, 0.01*R): continue
        vb2 = V*V - best[1]*best[1]
        if vb2 <= 0: continue
        got.append((R, V, eV, math.sqrt(vb2)))
    if got and len(got) >= 0.5*len(rings): rows[nm] = got

ALIAS = {
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
sparc_names = set()
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines_ = f.readlines()
start = max(i for i, l in enumerate(lines_)
            if set(l.strip()) <= set('- ')) + 1
for l in lines_[start:]:
    t = l.split()
    if len(t) >= 18:
        sparc_names.add(norm(t[0]))
overlap, ltonly = [], []
for nm in rows:
    cand = [nm] + ALIAS.get(nm, [])
    hit = [c for c in cand if c in sparc_names]
    (overlap if hit else ltonly).append(nm)

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_lam(y, lam):
    return (1.0-lam)*nu_standard(y) + lam*nu_be(y)

def build_ctx_lt(cut, subset):
    gb, go, sg, gid, sv = [], [], [], [], {}
    for k, nm in enumerate(sorted(subset)):
        m1 = t1[nm]
        if m1['i'] < 30: continue
        s_d2 = 0.10/LN10
        s_i2 = 2.0*math.radians(m1['ei']) / \
            max(math.tan(math.radians(m1['i'])), 1e-6)/LN10
        sv[k] = math.sqrt(s_d2**2 + s_i2**2)
        for (R, V, eV, Vb) in rows[nm]:
            if R <= 0 or V <= 0 or eV/V > cut: continue
            gbar = Vb*Vb/R*KPC
            if gbar <= 0: continue
            gb.append(gbar); go.append(V*V/R*KPC)
            sg.append(2*eV/V/math.log(10)); gid.append(k)
    return dict(gbar=np.array(gb), lgobs=np.log10(np.array(go)),
                sig2=np.array(sg)**2, gid=np.array(gid, dtype=int),
                SV=sv, mode='total')

def m2ll(th, lam, ctx, gset, von=True):
    la0, f, s_int = th
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4): return 1e12
    a0 = 10**la0
    gN = f*ctx['gbar']
    nv = nu_lam(gN/a0, lam)
    if np.min(nv) <= 0.05: return 1e12
    gm = gN*nv
    v = ctx['sig2'] + s_int*s_int
    a = ctx['lgobs'] - np.log10(gm)
    out = 0.0
    for g_ in gset:
        ji = np.where(ctx['gid'] == g_)[0]
        if len(ji) == 0: continue
        aj, vj = a[ji], v[ji]
        s_ = ctx['SV'].get(g_, 0.0)
        if (not von) or s_ <= 1e-4:
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

def fit_at(lam, ctx, gset, von, th_warm=None):
    starts = [[math.log10(A0_FID), 1.0, 0.08],
              [math.log10(A0_FID)+0.1, 0.8, 0.12]]
    if th_warm is not None: starts.insert(0, list(th_warm))
    best = None
    for th0 in starts:
        b = minimize(lambda t: m2ll(t, lam, ctx, gset, von), th0,
                     method='Nelder-Mead',
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-7))
        if best is None or b.fun < best.fun: best = b
    return best

LGX = np.round(np.arange(-2.00, 1.501, 0.05), 3)
def profile(ctx, gset, von):
    prof, th = [], None
    for lam in LGX:
        b = fit_at(lam, ctx, gset, von, th)
        prof.append(b.fun); th = b.x
    prof = np.array(prof)
    i = int(np.argmin(prof))
    lam_hat = LGX[i]
    if 0 < i < len(LGX)-1:
        x3, y3 = LGX[i-1:i+2], prof[i-1:i+2]
        c2_, c1_, _ = np.polyfit(x3, y3, 2)
        if c2_ > 0: lam_hat = -c1_/(2*c2_)
    return lam_hat

LAM_H, LAM_G = 0.97, -1.31

def contest(ctx, gset):
    fH = fit_at(LAM_H, ctx, gset, True).fun
    fG = fit_at(LAM_G, ctx, gset, True).fun
    return fH - fG   # > 0: the GD world fits better

def mock(ctx, gset, lam_t, rng):
    gN = ctx['gbar']
    mu = np.log10(gN*nu_lam(gN/A0_FID, lam_t))
    eps = rng.normal(0, np.sqrt(ctx['sig2'] + 0.08**2))
    voff = {g_: rng.normal(0, ctx['SV'].get(g_, 0.0)) for g_ in gset}
    ctx_m = dict(ctx)
    ctx_m['lgobs'] = mu + eps + np.array([voff[g] for g in ctx['gid']])
    return ctx_m

sub_only = [nm for nm in ltonly if nm in rows]

# ---------------- G9Rb-1: census regression -------------------------
cen = {}
for tag, cut in (('ARM-V', 0.10), ('ARM-R', 0.20)):
    ctx = build_ctx_lt(cut, sub_only)
    cen[tag] = (len(ctx['lgobs']), len(set(ctx['gid'].tolist())))
ok_1 = (cen['ARM-V'] == (65, 7) and cen['ARM-R'] == (172, 12)
        and len(sub_only) == 14)
gates['G9Rb-1'] = ok_1
P("G9Rb-1 census regression: ARM-V %s vs (65, 7); ARM-R %s vs "
  "(172, 12); LT-only %d vs 14 -> %s" % (cen['ARM-V'], cen['ARM-R'],
  len(sub_only), "PASS" if ok_1 else "FAIL"))
if not ok_1:
    save(); raise SystemExit(0)

# ---------------- G9Rb-0: engine regression (9R power cell) ---------
ctxV = build_ctx_lt(0.10, sub_only)
gsetV = sorted(set(ctxV['gid'].tolist()))
rng = np.random.default_rng(31)
lam_rep = profile(mock(ctxV, gsetV, -1.31, rng), gsetV, True)
ok_0 = abs(lam_rep - (-0.524)) <= 0.02
gates['G9Rb-0'] = ok_0
P("G9Rb-0 engine regression: 9R power cell (truth -1.31, s31, ARM-V) "
  "-> lam_hat = %+0.3f vs printed -0.524 -> %s"
  % (lam_rep, "PASS" if ok_0 else "FAIL"))
if not ok_0:
    P("STOP: the engines diverge; do not quote"); save()
    raise SystemExit(0)
P("")

# ---------------- calibration ---------------------------------------
def calibrate(tag, cut):
    ctx = build_ctx_lt(cut, sub_only)
    gset = sorted(set(ctx['gid'].tolist()))
    DH, DG = [], []
    for k in range(20):
        rH = np.random.default_rng(1000 + k)
        DH.append(contest(mock(ctx, gset, LAM_H, rH), gset))
        rG = np.random.default_rng(1500 + k)
        DG.append(contest(mock(ctx, gset, LAM_G, rG), gset))
    DH, DG = np.array(DH), np.array(DG)
    tstar = 0.5*(DH.mean() + DG.mean())
    okH = int(np.sum(DH < tstar))
    okG = int(np.sum(DG > tstar))
    powered = okH >= 16 and okG >= 16
    P("[%s cut %.2f] calibration (20 per world):" % (tag, cut))
    P("  headline-truth D: mean %+0.2f SD %.2f range [%+0.2f, %+0.2f] "
      "own-side %d/20" % (DH.mean(), DH.std(), DH.min(), DH.max(), okH))
    P("  GD-truth       D: mean %+0.2f SD %.2f range [%+0.2f, %+0.2f] "
      "own-side %d/20" % (DG.mean(), DG.std(), DG.min(), DG.max(), okG))
    P("  t* = %+0.2f -> %s (misclassified %d/40)"
      % (tstar, "POWERED" if powered else "UNPOWERED",
         (20-okH) + (20-okG)))
    return dict(ctx=ctx, gset=gset, DH=DH, DG=DG, tstar=tstar,
                powered=powered, mis=(20-okH)+(20-okG), cut=cut)

cal = {}
for tag, cut in (('ARM-V', 0.10), ('ARM-R', 0.20)):
    cal[tag] = calibrate(tag, cut)
P("")

arm = None
pw = [(t, c) for t, c in cal.items() if c['powered']]
if pw:
    arm = min(pw, key=lambda z: (z[1]['mis'],
              0 if z[0] == 'ARM-V' else 1))[0]

# ---------------- the sky read + verdict ----------------------------
if arm is None:
    P("==> 9R-b VERDICT (locked grammar): B-POWER-LIMITED -- the "
      "contest is not powered on either arm at Oh+15 grade; TODO 28 "
      "closes at this input grade (successors: Iorio+2017 curves, "
      "anchored/DR4-era dwarfs). The sky was not read.")
else:
    c = cal[arm]
    D_sky = contest(c['ctx'], c['gset'])
    P("[%s] SKY contest: D = m2ll(+0.97) - m2ll(-1.31) = %+0.3f "
      "(t* = %+0.2f; headline-truth range [%+0.2f, %+0.2f]; GD-truth "
      "range [%+0.2f, %+0.2f])" % (arm, D_sky, c['tstar'],
      c['DH'].min(), c['DH'].max(), c['DG'].min(), c['DG'].max()))
    # co-reads
    for tag, subset in (('LT-ALL', [nm for nm in rows]),
                        ('OVERLAP', [nm for nm in overlap
                                     if nm in rows])):
        if len(subset) < 3:
            P("[co-read %s] < 3 galaxies; skipped" % tag); continue
        ctx2 = build_ctx_lt(c['cut'], subset)
        gset2 = sorted(set(ctx2['gid'].tolist()))
        if len(ctx2['lgobs']) < 20:
            P("[co-read %s] < 20 points; skipped" % tag); continue
        P("[co-read %s] D = %+0.3f (%d pts / %d gal)"
          % (tag, contest(ctx2, gset2), len(ctx2['lgobs']),
             len(gset2)))
    P("")
    if D_sky > c['tstar'] and D_sky > c['DH'].max():
        P("==> 9R-b VERDICT (locked grammar): B-GD-SIDE -- the "
          "independent LITTLE THINGS dwarfs prefer the GD-dial world "
          "over the headline world beyond the full headline-truth "
          "calibration range (empirical ~p < 0.05). The SPARC GD "
          "tension direction REPLICATES at contest grade on "
          "drift-corrected independent curves. NOT a lam "
          "measurement; the profile stays power-limited (9R).")
    elif D_sky < c['tstar'] and D_sky < c['DG'].min():
        P("==> 9R-b VERDICT (locked grammar): B-HEADLINE-SIDE -- the "
          "independent dwarfs prefer the headline world beyond the "
          "full GD-truth calibration range: the dissolution "
          "direction at contest grade; the P2 tension section takes "
          "a caution-grade annotation at the next thaw.")
    else:
        P("==> 9R-b VERDICT (locked grammar): B-GRAY -- the sky D "
          "lands inside the calibration overlap; no side claimed at "
          "this grade (an intermediate/misspecified sky lands here "
          "by construction).")
P("    NO credence movement (pre-stated; measurement round).")

# ---------------- G9Rb-L --------------------------------------------
ok_l = False
for row in csv.reader(open('LEDGER.csv', encoding='utf-8')):
    if row and row[0] == 'gal-9r-ltrepl' and row[1] == 'CURRENT':
        ok_l = True
gates['G9Rb-L'] = ok_l
P("")
P("G9Rb-L ledger leg gal-9r-ltrepl CURRENT: %s"
  % ("PASS" if ok_l else "FAIL"))
P("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
  for k, v in sorted(gates.items())))
P("")
P("done (%.1f min)" % ((time.time()-t00)/60))
save()
print("\nsaved: data/stage9rb_contest.txt")
