"""
STAGE 9R-c -- THE TWO-WORLD CONTEST ON IORIO+17 CURVES (the TODO-28
reopen, executed per 9R-b's pre-registered successor clause: "Iorio+2017
curves ... near-miss says it likely tips").

SINGLE-AXIS CHANGE (the design principle): the OBSERVED curve (rings,
velocities, errors) swaps from Oh+15 to Iorio et al. 2017 (MNRAS 466,
4159; 3D-Barolo, asymmetric-drift-corrected circular velocities
Vc +/- err_Vc; canonical release = finalrot.zip from the second
author's site, fetched by calcs/fetch_iorio17.py, NO VizieR catalog
exists). EVERYTHING ELSE VERBATIM from 9R/9R-b: the Oh+15 baryonic
reconstruction (V_bar^2 = V_tot^2 - V_DM^2 at Oh rings), t1/t2
metadata (flags, e_i vertical channel, i >= 30 cut, e_D/D = 0.10),
nu_lam, m2ll, fit_at, the contest D = m2ll(+0.97) - m2ll(-1.31), mock
construction + seed scheme (headline default_rng(1000+k), GD
default_rng(1500+k), k = 0..19), powered bar >= 16/20 per world per
arm, arm cuts ARM-V 0.10 / ARM-R 0.20, viability >= 40 pts and >= 5
LT-only galaxies, verdict bars verbatim.

NEW WIRING (forced by cross-catalog ring grids), disclosed and gated:
  W1 V_bar at Iorio radii by LINEAR interpolation in V_bar^2 over the
     Oh ring set. NO extrapolation (Iorio rings outside Oh's V_bar
     span dropped, counted). Bracket-gap guard: keep a ring only if
     its Oh brackets satisfy (R_hi - R_lo) <= max(0.5 kpc, 0.30*R);
     wider-gap rings dropped, counted.
  W2 Distance harmonization, d = |D_Io - D_Oh|/D_Oh per galaxy:
     d <= 0.02 use as-is; 0.02 < d <= 0.10 rescale Oh radii AND
     V_bar^2 by (D_Io/D_Oh) (M ~ D^2, R ~ D => V^2 ~ D, first
     order), flagged; d > 0.10 galaxy DROPPED (listed).
  W3 DDO216 primary file only (ddo216b = "scenario2, see Iorio+16" =
     an alternative kinematic solution; existence disclosed, excluded
     from every set).
  W4 Iorio galaxies with no usable Oh mass model (absent or
     Spitzer-flagged) dropped (listed): no V_bar source.
  W5 Vertical channel keeps Oh table1 e_i (same physical galaxy;
     Iorio headers carry no inclination error) -- preserves the
     single-axis principle; Iorio's own i is NOT used anywhere.

PRE-REGISTERED BARS (locked; read on the powered arm with lower
misclassification, tie -> ARM-V; LT-ONLY set; verbatim 9R-b grammar):
  C-GD-SIDE       : D_sky > t* AND D_sky > max(headline-truth draws).
  C-HEADLINE-SIDE : D_sky < t* AND D_sky < min(GD-truth draws).
  C-GRAY          : anything else (D_sky quoted against both).
  C-POWER-LIMITED : no powered arm -> the Iorio clause of TODO 28 is
                    SPENT; the reopen re-closes at public-curve grade
                    (remaining successor: DR4-era dwarfs).
CO-READS (never verdict-bearing): LT-ALL and OVERLAP contest D at the
powered arm's cut; the unpowered arm's sky is NOT read.
CREDENCE: NO movement any branch (measurement round, pre-stated);
paper consequence booked at the next author-called thaw.

GATES:
  G9Rc-M  manifests: the four LT files match
          data/littlethings_manifest.sha256 AND the zip + 18 tables
          match data/iorio17_manifest.sha256.
  G9Rc-0  engine regression at archived-print grade (the 8P/8Q rule):
          the copied engine on the OH-INPUT context reproduces the
          9R-b archive EXACTLY -- ARM-V calibration t* rounds to
          +0.40, DH/DG means round to -1.25/+2.04, own-side counts
          17/20 and 14/20 exact; and the 9R power cell (truth -1.31,
          seed 31, ARM-V, profile) gives lam_hat = -0.524 within
          0.02. FAIL -> STOP, nothing downstream is quoted.
  G9Rc-1  census: per-galaxy disposition printed (rings read / kept /
          dropped by reason: no-Oh, distance, extrap, gap, Vbar2<=0,
          i<30). Bars: overall interp-related drops (extrap+gap) <=
          25% of rings of surviving galaxies; a galaxy keeping < 50%
          of its rings is dropped (listed).
  G9Rc-P  the power gate: >= 16/20 own-side per world per arm on the
          IORIO design + the 9R viability floor (>= 40 pts, >= 5
          LT-only galaxies).
  G9Rc-L  ledger legs gal-9r-ltrepl and gal-9rb-contest CURRENT.

DISCLOSURE (run before this commit: NOTHING of this stage). Seen
before commit: the 18-table listing; the ddo154 / ddo216 / ddo216b /
wlm headers (byte-format read; the scenario2 note); the archived 9R-b
output (the regression pins above). The ddo154 errors were SEEN to be
smaller than Oh's -- that is 9R-b's own stated reopen premise, not a
result of this stage. No fit, census, calibration, or contest has
touched Iorio input at commit time.

Output: data/stage9rc_contest.txt. Wall-clock: ~10-20 min CPU.
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
IO = 'data/iorio17'

OUT = []
def P(s=""):
    print(s, flush=True)
    OUT.append(s)

def save():
    with open('data/stage9rc_contest.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

t00 = time.time()
P("=" * 78)
P("STAGE 9R-c -- THE CONTEST ON IORIO+17 CURVES (headline +0.97 vs GD -1.31)")
P("=" * 78)

gates = {}

# ---------------- G9Rc-M: manifests ---------------------------------
ok_m = True
man = {}
for line in open('data/littlethings_manifest.sha256'):
    h, n = line.split()
    man[n.lstrip('*')] = h
for fn in ['table1.dat', 'table2.dat', 'rotdmbar.dat', 'rotdm.dat']:
    h = hashlib.sha256(open(os.path.join(LT, fn), 'rb').read()).hexdigest()
    ok_m &= (man.get(fn) == h)
man_io = {}
for line in open('data/iorio17_manifest.sha256'):
    h, n = line.split()
    man_io[n] = h
n_io = 0
for n, h in man_io.items():
    p = os.path.join('data', n)
    hh = hashlib.sha256(open(p, 'rb').read()).hexdigest()
    ok_m &= (hh == h)
    n_io += 1
ok_m &= (n_io == 19)
gates['G9Rc-M'] = ok_m
P("G9Rc-M manifests (4 LT + %d Iorio entries): %s"
  % (n_io, "PASS" if ok_m else "FAIL"))
if not ok_m:
    save(); raise SystemExit(0)

# ---------------- LT construction (verbatim 9R/9R-b) ----------------
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

def is_overlap(nm):
    cand = [nm] + ALIAS.get(nm, [])
    return any(c in sparc_names for c in cand)

def nu_be(y):
    x = np.sqrt(np.clip(y, 1e-14, None))
    return np.where(x > 40, 1.0, 1.0/(1.0-np.exp(-np.minimum(x, 40))))
def nu_standard(y):
    return np.sqrt((1.0+np.sqrt(1.0+4.0/np.clip(y, 1e-14, None)**2))/2.0)
def nu_lam(y, lam):
    return (1.0-lam)*nu_standard(y) + lam*nu_be(y)

def build_ctx(cut, subset, rows_src):
    gb, go, sg, gid, sv = [], [], [], [], {}
    for k, nm in enumerate(sorted(subset)):
        m1 = t1[nm]
        if m1['i'] < 30: continue
        s_d2 = 0.10/LN10
        s_i2 = 2.0*math.radians(m1['ei']) / \
            max(math.tan(math.radians(m1['i'])), 1e-6)/LN10
        sv[k] = math.sqrt(s_d2**2 + s_i2**2)
        for (R, V, eV, Vb) in rows_src[nm]:
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

def calibrate(tag, cut, subset, rows_src):
    ctx = build_ctx(cut, subset, rows_src)
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

sub_only_oh = sorted([nm for nm in rows if not is_overlap(nm)])

# ---------------- G9Rc-0: engine regression on the OH path ----------
P("")
P("G9Rc-0 engine regression (OH-input, verbatim 9R-b path):")
c0 = calibrate('ARM-V/OH', 0.10, sub_only_oh, rows)
ok_cal = (len(sub_only_oh) == 14
          and round(c0['tstar'], 2) == 0.40
          and round(float(c0['DH'].mean()), 2) == -1.25
          and round(float(c0['DG'].mean()), 2) == 2.04
          and int(np.sum(c0['DH'] < c0['tstar'])) == 17
          and int(np.sum(c0['DG'] > c0['tstar'])) == 14)
ctxV = build_ctx(0.10, sub_only_oh, rows)
gsetV = sorted(set(ctxV['gid'].tolist()))
rng = np.random.default_rng(31)
lam_rep = profile(mock(ctxV, gsetV, -1.31, rng), gsetV, True)
ok_cell = abs(lam_rep - (-0.524)) <= 0.02
gates['G9Rc-0'] = ok_cal and ok_cell
P("  archive match: LT-only %d vs 14, t* %+0.2f vs +0.40, means "
  "%+0.2f/%+0.2f vs -1.25/+2.04, counts %d/%d vs 17/14 -> %s"
  % (len(sub_only_oh), c0['tstar'], c0['DH'].mean(), c0['DG'].mean(),
     int(np.sum(c0['DH'] < c0['tstar'])),
     int(np.sum(c0['DG'] > c0['tstar'])),
     "PASS" if ok_cal else "FAIL"))
P("  9R power cell lam_hat = %+0.3f vs -0.524 -> %s"
  % (lam_rep, "PASS" if ok_cell else "FAIL"))
if not gates['G9Rc-0']:
    P("STOP: the engines diverge; do not quote"); save()
    raise SystemExit(0)

# ---------------- the Iorio loader + W1-W5 --------------------------
P("")
P("IORIO+17 CONSTRUCTION (W1 interpolation / W2 distance / W3 primary "
  "/ W4 mass-model / W5 metadata):")
rows_io = {}
census = []
drop_ext_tot = drop_gap_tot = kept_tot = read_tot = 0
for fp in sorted(glob.glob(os.path.join(IO, 'finalrot', '*_onlinetab.txt'))):
    base = os.path.basename(fp).replace('_onlinetab.txt', '')
    if base == 'ddo216b':
        census.append(('DDO216b', 'W3 scenario2 excluded', 0, 0))
        continue
    hdr, rr = {}, []
    for line in open(fp):
        if line.startswith('#'):
            m = re.match(r'#\s*Galaxy:\s*(\S+)', line)
            if m: hdr['name'] = m.group(1)
            m = re.match(r'#\s*Distance:\s*([\d.]+)', line)
            if m: hdr['D'] = float(m.group(1))
            continue
        t = line.split()
        if len(t) >= 8:
            rr.append((float(t[1]), float(t[6]), float(t[7])))
    nm = norm(hdr['name'])
    read_tot += len(rr)
    if nm not in rows:
        census.append((hdr['name'], 'W4 no usable Oh mass model',
                       len(rr), 0))
        continue
    D_oh = t1[nm]['D']
    drel = abs(hdr['D'] - D_oh)/D_oh
    scale = 1.0
    tag = 'D ok'
    if drel > 0.10:
        census.append((hdr['name'], 'W2 distance mismatch %.0f%%'
                       % (100*drel), len(rr), 0))
        continue
    if drel > 0.02:
        scale = hdr['D']/D_oh
        tag = 'D rescaled x%.3f' % scale
    Roh = np.array([r[0] for r in rows[nm]])*scale
    Vb2 = np.array([r[3]**2 for r in rows[nm]])*scale
    o = np.argsort(Roh); Roh, Vb2 = Roh[o], Vb2[o]
    got = []
    d_ext = d_gap = 0
    for (R, Vc, eVc) in rr:
        if R < Roh[0] or R > Roh[-1]:
            d_ext += 1; continue
        j = int(np.searchsorted(Roh, R))
        j = min(max(j, 1), len(Roh)-1)
        gapw = Roh[j] - Roh[j-1]
        if gapw > max(0.5, 0.30*R):
            d_gap += 1; continue
        vb2 = float(np.interp(R, Roh, Vb2))
        if vb2 <= 0: continue
        got.append((R, Vc, eVc, math.sqrt(vb2)))
    if len(got) < 0.5*len(rr):
        census.append((hdr['name'], '<50%% rings kept (%s)' % tag,
                       len(rr), len(got)))
        continue
    drop_ext_tot += d_ext; drop_gap_tot += d_gap
    rows_io[nm] = got
    kept_tot += len(got)
    census.append((hdr['name'], 'KEPT (%s; ext %d gap %d)'
                   % (tag, d_ext, d_gap), len(rr), len(got)))
for (name, disp, nr, nk) in census:
    P("  %-10s %3d -> %3d  %s" % (name, nr, nk, disp))
surv_read = sum(nr for (n_, d_, nr, nk) in census if nk > 0)
frac_drop = (drop_ext_tot + drop_gap_tot)/max(surv_read, 1)
ok_1 = frac_drop <= 0.25 and len(rows_io) >= 1
gates['G9Rc-1'] = ok_1
P("G9Rc-1 census: %d galaxies kept, %d/%d rings; interp drops "
  "ext %d + gap %d = %.1f%% of surviving-galaxy rings (bar 25%%) -> %s"
  % (len(rows_io), kept_tot, read_tot, drop_ext_tot, drop_gap_tot,
     100*frac_drop, "PASS" if ok_1 else "FAIL"))
if not ok_1:
    save(); raise SystemExit(0)

sub_only = sorted([nm for nm in rows_io if not is_overlap(nm)])
sub_over = sorted([nm for nm in rows_io if is_overlap(nm)])
P("LT-only (verdict set): %d  %s" % (len(sub_only), sub_only))
P("overlap (co-read):     %d  %s" % (len(sub_over), sub_over))
P("")

# ---------------- G9Rc-P: calibration on the Iorio design -----------
cal = {}
for tag, cut in (('ARM-V', 0.10), ('ARM-R', 0.20)):
    ctx = build_ctx(cut, sub_only, rows_io)
    npts = len(ctx['lgobs'])
    ngal = len(set(ctx['gid'].tolist()))
    P("[%s] census: %d pts / %d gal" % (tag, npts, ngal))
    if npts < 40 or ngal < 5:
        P("[%s] below viability floor (40 pts / 5 gal) -> not "
          "calibrated" % tag)
        cal[tag] = dict(powered=False, mis=99, cut=cut)
        continue
    cal[tag] = calibrate(tag, cut, sub_only, rows_io)
P("")
gates['G9Rc-P'] = any(c['powered'] for c in cal.values())

arm = None
pw = [(t, c) for t, c in cal.items() if c['powered']]
if pw:
    arm = min(pw, key=lambda z: (z[1]['mis'],
              0 if z[0] == 'ARM-V' else 1))[0]

# ---------------- the sky read + verdict ----------------------------
if arm is None:
    P("==> 9R-c VERDICT (locked grammar): C-POWER-LIMITED -- the "
      "contest is not powered on either arm at Iorio+17 grade; the "
      "TODO-28 Iorio clause is SPENT; the reopen re-closes at "
      "public-curve grade (remaining successor: DR4-era dwarfs). "
      "The sky was not read.")
else:
    c = cal[arm]
    D_sky = contest(c['ctx'], c['gset'])
    P("[%s] SKY contest: D = m2ll(+0.97) - m2ll(-1.31) = %+0.3f "
      "(t* = %+0.2f; headline-truth range [%+0.2f, %+0.2f]; GD-truth "
      "range [%+0.2f, %+0.2f])" % (arm, D_sky, c['tstar'],
      c['DH'].min(), c['DH'].max(), c['DG'].min(), c['DG'].max()))
    for tag, subset in (('LT-ALL', sorted(rows_io)),
                        ('OVERLAP', sub_over)):
        if len(subset) < 3:
            P("[co-read %s] < 3 galaxies; skipped" % tag); continue
        ctx2 = build_ctx(c['cut'], subset, rows_io)
        gset2 = sorted(set(ctx2['gid'].tolist()))
        if len(ctx2['lgobs']) < 20:
            P("[co-read %s] < 20 points; skipped" % tag); continue
        P("[co-read %s] D = %+0.3f (%d pts / %d gal)"
          % (tag, contest(ctx2, gset2), len(ctx2['lgobs']),
             len(gset2)))
    P("")
    if D_sky > c['tstar'] and D_sky > c['DH'].max():
        P("==> 9R-c VERDICT (locked grammar): C-GD-SIDE -- the "
          "independent Iorio+17 dwarfs prefer the GD-dial world over "
          "the headline world beyond the full headline-truth "
          "calibration range (empirical ~p < 0.05). The SPARC GD "
          "tension direction REPLICATES at contest grade on 3D-Barolo "
          "drift-corrected independent curves. NOT a lam measurement.")
    elif D_sky < c['tstar'] and D_sky < c['DG'].min():
        P("==> 9R-c VERDICT (locked grammar): C-HEADLINE-SIDE -- the "
          "independent dwarfs prefer the headline world beyond the "
          "full GD-truth calibration range: the dissolution direction "
          "at contest grade; the P2 tension section takes a "
          "caution-grade annotation at the next thaw.")
    else:
        P("==> 9R-c VERDICT (locked grammar): C-GRAY -- the sky D "
          "lands inside the calibration overlap; no side claimed at "
          "this grade (an intermediate/misspecified sky lands here "
          "by construction).")
P("    NO credence movement (pre-stated; measurement round).")

# ---------------- G9Rc-L --------------------------------------------
ok_l = {'gal-9r-ltrepl': False, 'gal-9rb-contest': False}
for row in csv.reader(open('LEDGER.csv', encoding='utf-8')):
    if row and row[0] in ok_l and row[1] == 'CURRENT':
        ok_l[row[0]] = True
gates['G9Rc-L'] = all(ok_l.values())
P("")
P("G9Rc-L ledger legs CURRENT: %s -> %s"
  % (ok_l, "PASS" if gates['G9Rc-L'] else "FAIL"))
P("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
  for k, v in sorted(gates.items())))
P("")
P("done (%.1f min)" % ((time.time()-t00)/60))
save()
print("\nsaved: data/stage9rc_contest.txt")
