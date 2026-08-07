"""
STAGE 9R -- THE GD-DIAL REPLICATION ON INDEPENDENT DWARFS (TODO 28; the
named Paper-2 falsifier #3, executable by us).

QUESTION: SPARC's gas-dominated dwarfs pull the c1 dial strongly
negative (8S-b OFF lam = -1.542; 8S-c vertON GD lam ~ -1.31 vs FULL
+0.96), a tension that survived eight in-catalog controls. Does an
INDEPENDENT dwarf data set reproduce it? LITTLE THINGS (Oh et al. 2015,
VizieR J/AJ/149/180): 26 gas-rich dwarfs, VLA HI curves at ~6",
asymmetric-drift CORRECTED by the data providers (the pressure-support
axis 8Y/9C could not reach in-catalog is already applied here --
interpretive sharpening either way).

DATA CONSTRUCTION (inside Oh+'s own published mass models): the VizieR
tables carry the total curve (rotdmbar.dat) and the DM-only curve
(rotdm.dat), each scaled by its own per-row (R0.3, V0.3). Un-scale both
to absolute units; match rings by absolute radius (tolerance
max(0.005 kpc, 1% of R)); reconstruct the baryonic curve pointwise:
V_bar^2 = V_tot^2 - V_DM^2. Rings that fail matching or give
V_bar^2 <= 0 are dropped and counted. Galaxies flagged in table2
(f_alphamin / f_alpha3.6 nonempty = no Spitzer 3.6um; their "DM"
includes the stars) are EXCLUDED entirely. Galaxy-level GD selector:
Mgas > Mstar (SED; kinematic fallback) -- expected nearly all.
Overlap split: normalized-name + alias match against the SPARC master
list; the verdict reads on LT-ONLY (non-SPARC); the overlap subset is
a cross-catalog co-read. Both memberships printed in full.

INSTRUMENT (the 8S-c dial, ported): profile of lam over the verbatim
grid [-2.00, 1.50] step 0.05 with nu_lam = (1-lam)*nu_standard +
lam*nu_be; free (log10 a0, f, s_int) per lam (verbatim priors); the
measured vertical channel (distance/inclination random intercept,
analytic marginal, formula verbatim from 8S-c). ONE wiring difference,
forced by the data: LT has no gas/disk split, so f multiplies the
TOTAL baryonic g_bar (SPARC wiring: f multiplies the disk only). The
wiring is LICENSED by a SPARC-side test (below); if the license bar
fails, the letters carry a wiring caveat.

Vertical-channel inputs for LT: e_i from table1 (measured); distances
have no printed error -> e_D/D = 0.10 FIDUCIAL (typical for the
Hunter+2012 distance methods), with scan co-reads at 0.05 / 0.20.

ARM LADDER (pre-registered fallback order; the SPARC loader's verbatim
point cut eV/V <= 0.10 may decimate LT's low-V dwarfs):
  ARM-V (verbatim cut 0.10) -> ARM-R (relaxed cut 0.20).
  An arm is VIABLE if census >= 40 points AND >= 5 LT-only galaxies
  AND its power gate passes. Letters read on the FIRST viable arm;
  the other arm is a co-read. If neither is viable: R-POWER-LIMITED.

POWER GATE G9R-P (per arm, before any sky read): inject two worlds on
the LT-only design (same radii, errors, s_vert; noise = point sigma
(+) s_int 0.08 (+) per-galaxy vertical draw): lam_true = +0.97 (the
SPARC headline) and -1.31 (the GD dial), seeds 31/101. PASS iff all
four recoveries land within +-0.5 of truth AND the recovered
separation lam(+0.97) - lam(-1.31) >= 1.0 in both seeds.

PRE-REGISTERED BARS (locked; read on the viable arm's LT-ONLY vertON
profile, D1 interval [lo, hi]):
  R-REPLICATES : lam_hat <= -0.5 AND hi < +0.4 (the independent dwarfs
                 reproduce the GD direction and exclude the headline).
  R-DISSOLVES  : lam_hat >= +0.3 AND lo > -0.5 (the independent dwarfs
                 sit at the headline world and exclude the GD dial ->
                 the SPARC-GD tension re-reads as a SPARC-subset
                 systematic; public-retraction-grade consequence for
                 the P2 tension section at the next thaw).
  R-GRAY       : anything else (quote the interval).
  R-POWER-LIMITED : no viable arm.
CO-READS (never verdict-bearing): LT-ALL profile; LT-overlap profile;
e_D/D scan rows; the other arm's profile; the drift-correction note
(LT curves are pressure-corrected at source: REPLICATES would further
weaken the pressure-support suspect; DISSOLVES would revive it).

GATES:
  G9R-M  manifest: SHA256 of the four fetched files match
         data/littlethings_manifest.sha256.
  G9R-0  instrument-port regression: SPARC-GD OFF profile under the
         ORIGINAL wiring reproduces 8S-b/8S-c lam = -1.542 within one
         grid step (0.1; the 8S-c G8Sc-0 bar), and the OFF branch
         equals the directly-summed flat objective at a probe theta
         (<= 1e-9, the bit-identity).
  G9R-W  wiring license: SPARC-GD vertON under the LT-style
         f-on-total wiring lands within 0.3 of the archived -1.31
         (8S-c). PASS -> the LT wiring is licensed; FAIL -> letters
         carry the wiring caveat (disclosed, not blocking).
  G9R-1  LT parse census: 26 rows in table1/table2; flags in [0, 6];
         ring-match fraction >= 0.5 overall; V_bar^2 <= 0 drops <= 15%
         of matched rings; per-galaxy match < 50% -> galaxy dropped
         (listed).
  G9R-P  the power gate (above).
  G9R-L  ledger legs live: at least one CURRENT row each at stages
         8S-b/8S-c and 8V (the dial + its bootstrap are live objects).
CREDENCE: NO movement either way (measurement round, pre-stated); the
paper consequence is booked at the next author-called thaw.

DISCLOSURE: sketched pre-commit: the CVnIdwA absolute-radius alignment
hand-check (0.190 kpc both files, first-ring V_bar/V_tot ~ 0.98) and
the observation that CVnIdwA's rings fail the verbatim 0.10 cut (the
reason the arm ladder exists). NOT run pre-commit: every census count,
every profile, the power gate, both SPARC legs. The replication
outcome is a genuine unknown to us at commit time.

Output: data/stage9r_gdrepl.txt. Wall-clock: ~10-30 min CPU (Nelder-
Mead profiles; no GPU).
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
UD, UB = 0.5, 0.7
A0_FID = 1.2e-10
LN10 = math.log(10)
LT = 'data/littlethings'

OUT = []
def P(s=""):
    print(s, flush=True)
    OUT.append(s)

def save():
    with open('data/stage9r_gdrepl.txt', 'w', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")

t00 = time.time()
P("=" * 78)
P("STAGE 9R -- GD-DIAL REPLICATION ON LITTLE THINGS (Oh+2015)")
P("=" * 78)

gates = {}

# ---------------- G9R-M: manifest -----------------------------------
man = {}
for line in open('data/littlethings_manifest.sha256'):
    h, n = line.split()
    man[n.lstrip('*')] = h
ok_m = True
for fn in ['table1.dat', 'table2.dat', 'rotdmbar.dat', 'rotdm.dat']:
    h = hashlib.sha256(open(os.path.join(LT, fn), 'rb').read()).hexdigest()
    ok_m &= (man.get(fn) == h)
gates['G9R-M'] = ok_m
P("G9R-M manifest (4 files vs committed sha256): %s"
  % ("PASS" if ok_m else "FAIL"))
if not ok_m:
    P("STOP: manifest mismatch; refetch via calcs/fetch_littlethings.py")
    save(); raise SystemExit(0)

# ---------------- LT parse ------------------------------------------
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

ok_1 = (len(t1) == 26 and len(t2) == 26)
nflag = sum(1 for v in t2.values() if v['flag'])
ok_1 &= (0 <= nflag <= 6)
P("LT parse: table1 %d, table2 %d galaxies; flagged (no-3.6um, "
  "excluded): %d [%s]" % (len(t1), len(t2), nflag,
  ", ".join(t1[k]['name'] for k in t2 if t2[k]['flag'])))

# ring matching + baryon reconstruction
rows = {}
n_tot_rings = n_matched = n_negbar = 0
drop_gal = []
for nm, rings in tot.items():
    if nm not in t2 or t2[nm]['flag'] or nm not in t1: continue
    dm = dmo.get(nm, [])
    got = []
    for (R, V, eV) in rings:
        n_tot_rings += 1
        best, bd = None, 1e9
        for (Rd, Vd, _) in dm:
            d = abs(Rd - R)
            if d < bd: bd, best = d, (Rd, Vd)
        if best is None or bd > max(0.005, 0.01*R): continue
        n_matched += 1
        vb2 = V*V - best[1]*best[1]
        if vb2 <= 0:
            n_negbar += 1; continue
        got.append((R, V, eV, math.sqrt(vb2)))
    if len(got) < 0.5*len(rings):
        drop_gal.append((t1[nm]['name'], len(got), len(rings))); continue
    if got: rows[nm] = got
frac_match = n_matched/max(n_tot_rings, 1)
frac_neg = n_negbar/max(n_matched, 1)
ok_1 &= frac_match >= 0.5 and frac_neg <= 0.15
gates['G9R-1'] = ok_1
P("rings: %d total, matched %.2f, V_bar^2<=0 dropped %.3f; per-galaxy "
  "<50%% matched -> dropped: %s" % (n_tot_rings, frac_match, frac_neg,
  [d[0] for d in drop_gal] if drop_gal else "none"))
P("G9R-1 census: %s" % ("PASS" if ok_1 else "FAIL"))
if not ok_1:
    P("STOP: census gate failed; do not quote"); save(); raise SystemExit(0)

# GD selector + overlap split
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

gd_lt, dd_lt, overlap, ltonly = [], [], [], []
for nm in rows:
    m2 = t2[nm]
    mstar = m2['MstarSED'] if m2['MstarSED'] is not None else \
            (m2['MstarK'] if m2['MstarK'] is not None else 0.0)
    isgd = (m2['Mgas'] is not None) and m2['Mgas'] > mstar
    (gd_lt if isgd else dd_lt).append(nm)
    cand = [nm] + ALIAS.get(nm, [])
    hit = [c for c in cand if c in sparc_names]
    if hit: overlap.append((nm, hit[0]))
    else: ltonly.append(nm)
P("GD census (Mgas > Mstar): GD %d / non-GD %d [%s]"
  % (len(gd_lt), len(dd_lt),
     ", ".join(t1[n]['name'] for n in dd_lt) if dd_lt else "-"))
P("overlap vs SPARC: %d [%s]" % (len(overlap),
  ", ".join("%s=%s" % (t1[a]['name'], b) for a, b in overlap)))
P("LT-ONLY (verdict set): %d [%s]" % (len(ltonly),
  ", ".join(t1[n]['name'] for n in ltonly)))
P("")

# ---------------- instrument (8S-c port) ----------------------------
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
    if ctx['mode'] == 'total':
        gN = f*ctx['gbar']
    else:
        gN = ctx['g_gas'] + f*ctx['g_dsk'] + ctx['g_bul']
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
    lo = hi = None
    for j in range(i, -1, -1):
        if prof[j] > prof[i]+1.0:
            lo = np.interp(prof[i]+1.0, [prof[j+1], prof[j]],
                           [LGX[j+1], LGX[j]]) \
                if prof[j] != prof[j+1] else LGX[j]
            break
    for j in range(i, len(LGX)):
        if prof[j] > prof[i]+1.0:
            hi = np.interp(prof[i]+1.0, [prof[j-1], prof[j]],
                           [LGX[j-1], LGX[j]])
            break
    return lam_hat, lo, hi, i, prof

# ---------------- SPARC legs (G9R-0 regression + G9R-W license) -----
meta_s = {}
for l in lines_[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta_s[t[0]] = (float(t[5]), int(t[17]), float(t[2]),
                        float(t[3]), float(t[6]))
    except ValueError:
        continue
g_gas, g_dsk, g_bul, gobs_s, sig_s, gid_s = [], [], [], [], [], []
sv_s = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                        recursive=True))
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta_s.get(name, (0, 3, 1.0, 0.0, 0.0))
    if inc < 30 or q > 2: continue
    s_d2 = (eD/max(D, 1e-6))/LN10
    s_i2 = 2.0*math.radians(einc) / \
        max(math.tan(math.radians(inc)), 1e-6)/LN10
    sv_s[gi] = math.sqrt(s_d2**2 + s_i2**2)
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gg = Vg*abs(Vg)/R*KPC; gd = UD*Vd*abs(Vd)/R*KPC
        gb = UB*Vb*Vb/R*KPC
        if gg+gd+gb <= 0: continue
        g_gas.append(gg); g_dsk.append(gd); g_bul.append(gb)
        gobs_s.append(Vo*Vo/R*KPC); sig_s.append(2*eV/Vo/math.log(10))
        gid_s.append(gi)
g_gas, g_dsk, g_bul = map(np.array, (g_gas, g_dsk, g_bul))
ctx_sp = dict(g_gas=g_gas, g_dsk=g_dsk, g_bul=g_bul,
              gbar=g_gas+g_dsk+g_bul,
              lgobs=np.log10(np.array(gobs_s)),
              sig2=np.array(sig_s)**2,
              gid=np.array(gid_s, dtype=int), SV=sv_s, mode='sparc')
allg_s = np.unique(ctx_sp['gid'])
gdfrac = {}
for g_ in allg_s:
    m = ctx_sp['gid'] == g_
    gdfrac[g_] = float(np.mean(g_gas[m] > g_dsk[m] + g_bul[m]))
gd_sp = [g_ for g_ in allg_s if gdfrac[g_] >= 0.5]
P("SPARC leg: %d galaxies, GD %d (verbatim selector)"
  % (len(allg_s), len(gd_sp)))

th_p = [math.log10(A0_FID), 1.0, 0.08]
va = m2ll(th_p, 0.7, ctx_sp, gd_sp, von=False)
a0p = 10**th_p[0]
gNf = g_gas + th_p[1]*g_dsk + g_bul
af = ctx_sp['lgobs'] - np.log10(gNf*nu_lam(gNf/a0p, 0.7))
vf = ctx_sp['sig2'] + th_p[2]**2
wsel = np.isin(ctx_sp['gid'], gd_sp)
vflat = float(np.sum((af*af/vf + np.log(vf))[wsel]))
ok_bit = abs(va - vflat) <= 1e-9
lam0, _, _, _, _ = profile(ctx_sp, gd_sp, von=False)
ok_reg = abs(lam0 - (-1.542)) <= 0.1
gates['G9R-0'] = ok_bit and ok_reg
P("G9R-0 port regression: OFF-branch bit identity d = %.2e -> %s; "
  "SPARC-GD OFF lam_hat = %.3f vs -1.542 -> %s"
  % (va - vflat, "PASS" if ok_bit else "FAIL", lam0,
     "PASS" if ok_reg else "FAIL"))

ctx_spT = dict(ctx_sp); ctx_spT['mode'] = 'total'
lamW, loW, hiW, _, _ = profile(ctx_spT, gd_sp, von=True)
ok_w = abs(lamW - (-1.31)) <= 0.3
gates['G9R-W'] = ok_w
P("G9R-W wiring license: SPARC-GD vertON under f-on-TOTAL wiring "
  "lam_hat = %.3f (archived disk-wiring -1.31) -> %s"
  % (lamW, "LICENSED" if ok_w else "CAVEAT (|d| > 0.3)"))
P("")
if not gates['G9R-0']:
    P("STOP: the instrument port failed its regression; do not quote")
    save(); raise SystemExit(0)

# ---------------- arms ----------------------------------------------
def power_gate(ctx, gset, tag):
    okp = True
    rec = {}
    for lam_t in (0.97, -1.31):
        for seed in (31, 101):
            rng = np.random.default_rng(seed)
            gN = ctx['gbar']
            mu = np.log10(gN*nu_lam(gN/A0_FID, lam_t))
            eps = rng.normal(0, np.sqrt(ctx['sig2'] + 0.08**2))
            voff = {g_: rng.normal(0, ctx['SV'].get(g_, 0.0))
                    for g_ in gset}
            lg = mu + eps + np.array([voff[g] for g in ctx['gid']])
            ctx_m = dict(ctx); ctx_m['lgobs'] = lg
            lam_r, _, _, _, _ = profile(ctx_m, gset, von=True)
            rec[(lam_t, seed)] = lam_r
            okp &= abs(lam_r - lam_t) <= 0.5
    for seed in (31, 101):
        okp &= (rec[(0.97, seed)] - rec[(-1.31, seed)]) >= 1.0
    P("G9R-P power [%s]: " % tag + "; ".join(
        "true %+0.2f s%d -> %+0.3f" % (lt_, sd, rec[(lt_, sd)])
        for (lt_, sd) in rec) + " -> %s" % ("PASS" if okp else "FAIL"))
    return okp

verdict_arm = None
arm_results = {}
for tag, cut in (('ARM-V', 0.10), ('ARM-R', 0.20)):
    sub_only = [nm for nm in ltonly if nm in rows]
    ctx = build_ctx_lt(cut, sub_only)
    gset = sorted(set(ctx['gid'].tolist()))
    npts, ngal = len(ctx['lgobs']), len(gset)
    P("[%s cut %.2f] LT-only census: %d points / %d galaxies"
      % (tag, cut, npts, ngal))
    viable = npts >= 40 and ngal >= 5
    if viable and verdict_arm is None:
        viable = power_gate(ctx, gset, tag)
    if not viable:
        P("[%s] NOT VIABLE (census or power)" % tag)
        arm_results[tag] = None
        continue
    lam_h, lo, hi, i, prof = profile(ctx, gset, von=True)
    edge = 'INTERIOR' if 0 < i < len(LGX)-1 else 'EDGE'
    arm_results[tag] = (lam_h, lo, hi, edge, npts, ngal)
    P("[%s] LT-ONLY vertON: lam_hat = %+0.3f (D1 %s..%s) %s -> "
      "c1_hat = %+0.3f" % (tag, lam_h,
      "?" if lo is None else "%+0.3f" % lo,
      "?" if hi is None else "%+0.3f" % hi, edge, lam_h/2))
    if verdict_arm is None:
        verdict_arm = tag
P("")

# co-reads on the verdict arm's cut
if verdict_arm is not None:
    cutv = 0.10 if verdict_arm == 'ARM-V' else 0.20
    for tag, subset in (('LT-ALL', [nm for nm in rows]),
                        ('LT-OVERLAP', [a for a, b in overlap
                                        if a in rows])):
        if len(subset) < 3:
            P("[co-read %s] < 3 galaxies; skipped" % tag); continue
        ctx = build_ctx_lt(cutv, subset)
        gset = sorted(set(ctx['gid'].tolist()))
        if len(ctx['lgobs']) < 20:
            P("[co-read %s] < 20 points; skipped" % tag); continue
        lam_h, lo, hi, i, _ = profile(ctx, gset, von=True)
        P("[co-read %s] lam_hat = %+0.3f (D1 %s..%s), %d pts / %d gal"
          % (tag, lam_h, "?" if lo is None else "%+0.3f" % lo,
             "?" if hi is None else "%+0.3f" % hi,
             len(ctx['lgobs']), len(gset)))
    # e_D scan (rebuild with scanned s_d)
    for ed in (0.05, 0.20):
        gb, go, sg, gid, sv = [], [], [], [], {}
        for k, nm in enumerate(sorted([n for n in ltonly if n in rows])):
            m1 = t1[nm]
            if m1['i'] < 30: continue
            s_d2 = ed/LN10
            s_i2 = 2.0*math.radians(m1['ei']) / \
                max(math.tan(math.radians(m1['i'])), 1e-6)/LN10
            sv[k] = math.sqrt(s_d2**2 + s_i2**2)
            for (R, V, eV, Vb) in rows[nm]:
                if R <= 0 or V <= 0 or eV/V > cutv: continue
                gb.append(Vb*Vb/R*KPC); go.append(V*V/R*KPC)
                sg.append(2*eV/V/math.log(10)); gid.append(k)
        ctx = dict(gbar=np.array(gb), lgobs=np.log10(np.array(go)),
                   sig2=np.array(sg)**2,
                   gid=np.array(gid, dtype=int), SV=sv, mode='total')
        gset = sorted(set(gid))
        lam_h, lo, hi, _, _ = profile(ctx, gset, von=True)
        P("[co-read e_D/D = %.2f] LT-only lam_hat = %+0.3f (D1 %s..%s)"
          % (ed, lam_h, "?" if lo is None else "%+0.3f" % lo,
             "?" if hi is None else "%+0.3f" % hi))
P("")

# ---------------- G9R-L + verdict -----------------------------------
have = {'8S-b': False, '8S-c': False, '8V': False}
for row in csv.reader(open('LEDGER.csv', encoding='utf-8')):
    if len(row) > 5 and row[5] in have and row[1] == 'CURRENT':
        have[row[5]] = True
ok_l = have['8S-c'] and have['8V'] and have['8S-b']
gates['G9R-L'] = ok_l
P("G9R-L ledger legs (8S-b/8S-c/8V CURRENT rows): %s"
  % ("PASS" if ok_l else "FAIL"))

allok = all(v for v in gates.values() if v is not None)
P("GATES: " + "  ".join("%s:%s" % (k, "PASS" if v else "FAIL")
  for k, v in sorted(gates.items())))
P("")
if verdict_arm is None:
    P("==> 9R VERDICT (locked grammar): R-POWER-LIMITED -- no arm "
      "passed census + power; no replication read at this data grade.")
else:
    lam_h, lo, hi, edge, npts, ngal = arm_results[verdict_arm]
    wc = "" if gates['G9R-W'] else " [WIRING CAVEAT: license bar failed]"
    if lam_h <= -0.5 and (hi is not None and hi < 0.4):
        P("==> 9R VERDICT (locked grammar): R-REPLICATES%s -- the "
          "independent LITTLE THINGS dwarfs reproduce the GD dial "
          "direction (lam_hat = %+0.3f, D1 %+0.3f..%+0.3f on %d pts / "
          "%d galaxies, %s arm) and exclude the SPARC headline. The "
          "GD tension is NOT a SPARC-subset artifact; and these "
          "curves are drift-corrected at source, further weakening "
          "the pressure-support suspect."
          % (wc, lam_h, lo if lo is not None else float('nan'),
             hi, npts, ngal, verdict_arm))
    elif lam_h >= 0.3 and (lo is not None and lo > -0.5):
        P("==> 9R VERDICT (locked grammar): R-DISSOLVES%s -- the "
          "independent dwarfs sit at the headline world (lam_hat = "
          "%+0.3f, D1 %+0.3f..%+0.3f) and exclude the GD dial: the "
          "SPARC-GD tension re-reads as a SPARC-subset systematic. "
          "P2's tension section takes a correction-grade update at "
          "the next thaw." % (wc, lam_h,
             lo, hi if hi is not None else float('nan')))
    else:
        P("==> 9R VERDICT (locked grammar): R-GRAY%s -- lam_hat = "
          "%+0.3f (D1 %s..%s); neither letter's exclusion fired; "
          "interval quoted, no attribution."
          % (wc, lam_h, "?" if lo is None else "%+0.3f" % lo,
             "?" if hi is None else "%+0.3f" % hi))
P("    NO credence movement (pre-stated; measurement round).")
P("")
P("done (%.1f min)" % ((time.time()-t00)/60))
save()
print("\nsaved: data/stage9r_gdrepl.txt")
