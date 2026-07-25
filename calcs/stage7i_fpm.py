"""
STAGE 7I-FPM: MATERIAL-branch investigation, part 2 (fit-level). The 7I
strict-multiplicity fit collapsed (alpha_hat 0.4-0.75, Newton +4..+20)
while the model-light median boost SURVIVED and slightly strengthened
(1.086 -> 1.103 raw, 1.151 -> 1.185 perpendicular; stage7i_material.py)
— pointing at H2: the fixed proper-motion noise inflation FPM = 1.5,
calibrated on the FULL sample (3N/3O), overestimates the clean strict
subsample's noise (the flagged 43% of pairs — RUWE >= 1.2, overluminous
= companion wobble — are plausibly what the 1.5 was measuring), and the
over-broadened model pre-explains the boost width.

TEST (expectations stated before running): profile FPM in {1.0, 1.2, 1.5}
on the strict sample, seeds 31/101, both laws, absolute best-lnL printed
so the FPM profile is readable across runs (same data -> direct
likelihood comparison; FPM = a model nuisance being profiled).
  H2 CONFIRMED if the strict sample prefers FPM < 1.5 (>= 3/4 seed-law
     rows) AND alpha_hat at the preferred FPM recovers to >= 0.9
     (baselines 1.10-1.26).
  H1 SUPPORTED if alpha_hat stays <= 0.8 at every FPM (the collapse is
     not a noise-model artifact -> the companion fence question reopens
     at full strength and 7J becomes decisive).
  MIXED otherwise -> carried to 7J.
Writes data/stage7i_fpm.txt + appends verdict to data/stage7i_verdict.txt.
"""
import os
import re
import sys

import numpy as np

src0 = open('calcs/stage3p_v7budget.py', encoding='utf-8-sig').read()

OLD_DATA = r"""vt_d = (4.74047/plx[ok]*np.hypot(d['pmra1'][ok]-d['pmra2'][ok],
        d['pmdec1'][ok]-d['pmdec2'][ok]))/(0.9417*np.sqrt(Mtot_d/s_d))
dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
vx_ = d['pmra2']-d['pmra1']; vy_ = d['pmdec2']-d['pmdec1']"""
NEW_DATA = r"""dec_m = np.deg2rad(0.5*(d['dec1']+d['dec2']))
dra = (d['ra2']-d['ra1']+180.0) % 360.0 - 180.0
sx_ = dra*np.cos(dec_m); sy_ = d['dec2']-d['dec1']
sn_ = np.maximum(np.hypot(sx_, sy_), 1e-12)
sux_, suy_ = sx_/sn_, sy_/sn_
_cn = list(d.columns.names)
def _pick(*names):
    for n in names:
        if n in _cn: return np.asarray(d[n], dtype=np.float64)
    raise KeyError(names)
r1_ = _pick('radial_velocity1', 'dr2_radial_velocity1', 'rv1')
r2_ = _pick('radial_velocity2', 'dr2_radial_velocity2', 'rv2')
try:
    er1_ = _pick('radial_velocity_error1', 'dr2_radial_velocity_error1')
    er2_ = _pick('radial_velocity_error2', 'dr2_radial_velocity_error2')
except KeyError:
    er1_ = np.full(len(r1_), 2.0); er2_ = np.full(len(r2_), 2.0)
h1_, h2_ = np.isfinite(r1_), np.isfinite(r2_)
w1_ = np.where(h1_, 1.0/np.maximum(er1_, 0.5)**2, 0.0)
w2_ = np.where(h2_, 1.0/np.maximum(er2_, 0.5)**2, 0.0)
rvs_ = (np.where(h1_, r1_, 0.0)*w1_ + np.where(h2_, r2_, 0.0)*w2_) \
       / np.maximum(w1_+w2_, 1e-12)
th_ = sep/(2.06265e8/plx)
pmcor_ = rvs_*th_*plx/4.74047
vx_ = d['pmra2']-d['pmra1'] + pmcor_*sux_
vy_ = d['pmdec2']-d['pmdec1'] + pmcor_*suy_
vt_d = (4.74047/plx[ok]*np.hypot(vx_[ok], vy_[ok]))/(0.9417*np.sqrt(Mtot_d/s_d))"""

OK_OLD = r"""ok = (Rch<0.01)&(plx1>5)&(plx2>5)&(plx1/np.maximum(eplx1,1e-6)>20) \
   &(plx2/np.maximum(eplx2,1e-6)>20)&(np.abs(plx1-plx2)<3*np.hypot(eplx1,eplx2)) \
   &(sep>200)&(sep<50000)&(MG1>2.6)&(MG1<14.2)&(MG2>2.6)&(MG2<14.2)&(sigv<0.03)"""
OK_NEW = OK_OLD + """
ok0_ = int(ok.sum())
ok = ok & np.load('data/stage7i_strictmask.npy')
print(f"7I strict cut active: {int(ok.sum())}/{ok0_} pairs", flush=True)"""

OLD_SUM = """    P(f"  seed {seed}: BE-minus-simple best lnL = "
      f"{best_lnl['BE']-best_lnl['simple']:+.1f}  "
      f"({(time.time()-t0)/60:.1f} min)")"""
NEW_SUM = """    P(f"  seed {seed}: abs best lnL simple {best_lnl['simple']:+.3f} "
      f"BE {best_lnl['BE']:+.3f}  ({(time.time()-t0)/60:.1f} min)")"""
OUT_OLD = "with open('data/stage3u_summary.txt', 'a') as f:"

def make_src(fpm):
    s = src0
    assert s.count(OLD_DATA) == 1
    s = s.replace(OLD_DATA, NEW_DATA)
    assert s.count(OK_OLD) == 1
    s = s.replace(OK_OLD, OK_NEW)
    assert s.count("FPM = 1.5") == 1
    s = s.replace("FPM = 1.5", f"FPM = {fpm}")
    assert s.count(OLD_SUM) == 1
    s = s.replace(OLD_SUM, NEW_SUM)
    assert s.count(OUT_OLD) == 1
    s = s.replace(OUT_OLD, "with open('data/stage7i_fpm.txt', 'a') as f:")
    # tag rows with the FPM value so the parse can separate variants
    assert s.count('P(f"seed {seed} {law}: ') == 1
    s = s.replace('P(f"seed {seed} {law}: ',
                  'P(f"fpm ' + str(fpm) + ' seed {seed} {law}: ')
    return s

ROW = re.compile(r"fpm ([0-9.]+) seed (\d+) (simple|BE): a_hat=([0-9.]+) "
                 r"\(grid ([0-9.]+), interior=(\w+)\), "
                 r"dlnL\(Newton\)=([+-][0-9.]+), wr=([0-9.]+)")
ABS = re.compile(r"seed (\d+): abs best lnL simple ([+-][0-9.]+) "
                 r"BE ([+-][0-9.]+)")

def parse():
    rows, absl = {}, {}
    if os.path.exists('data/stage7i_fpm.txt'):
        txt = open('data/stage7i_fpm.txt').read()
        for m in ROW.finditer(txt):
            f_, s_, law, ah, _, inter, dn, wr = m.groups()
            rows[(float(f_), int(s_), law)] = dict(
                a=float(ah), interior=(inter == 'True'), dn=float(dn),
                wr=float(wr))
        # abs lines carry no fpm tag; associate by order of appearance
        fpms = [float(m.group(1)) for m in ROW.finditer(txt)][::2]
        for (f_), m in zip(fpms, ABS.finditer(txt)):
            s_ = int(m.group(1))
            absl[(f_, s_, 'simple')] = float(m.group(2))
            absl[(f_, s_, 'BE')] = float(m.group(3))
    return rows, absl

SEEDS = [31, 101]
for fpm in (1.5, 1.2, 1.0):
    rows, _ = parse()
    todo = [s_ for s_ in SEEDS if (fpm, s_, 'simple') not in rows]
    if todo:
        print(f"===== 7I-FPM {fpm}: seeds {todo} =====", flush=True)
        sys.argv = ['stage7i_fpm', '1p2'] + [str(s_) for s_ in todo]
        exec(compile(make_src(fpm), f'stage3p_7i_fpm{fpm}', 'exec'),
             {'__name__': '__main__'})

rows, absl = parse()
V = ["", "FPM profile on the strict sample (7I MATERIAL part 2; "
     "H2 = prefers FPM<1.5 with a_hat recovering >=0.9):"]
n_pref, n_rec = 0, 0
for s_ in SEEDS:
    for law in ('simple', 'BE'):
        prof = {f_: absl.get((f_, s_, law), np.nan) for f_ in (1.0, 1.2, 1.5)}
        fhat = max(prof, key=lambda k: (prof[k] if np.isfinite(prof[k])
                                        else -np.inf))
        r = rows[(fhat, s_, law)]
        V.append(f"  seed {s_} {law}: lnL(1.0/1.2/1.5) = "
                 f"{prof[1.0]:+.1f}/{prof[1.2]:+.1f}/{prof[1.5]:+.1f} -> "
                 f"FPM_hat {fhat}, a_hat {r['a']:.2f} "
                 f"(interior {r['interior']}), Newton {r['dn']:+.1f}, "
                 f"wr {r['wr']}")
        n_pref += fhat < 1.5
        n_rec += (fhat < 1.5) and (r['a'] >= 0.9)
verdict = ('H2 CONFIRMED' if (n_pref >= 3 and n_rec >= 3) else
           'H1 SUPPORTED' if n_rec == 0 and n_pref <= 1 else 'MIXED')
V.append(f"  FPM<1.5 preferred {n_pref}/4; recovered a_hat>=0.9 at FPM_hat "
         f"{n_rec}/4  ==> {verdict} (bars in module docstring)")
out = "\n".join(V)
print(out)
with open('data/stage7i_verdict.txt', 'a') as f:
    f.write(out + "\n")
print("\nSTAGE 7I-FPM done -> appended data/stage7i_verdict.txt")
