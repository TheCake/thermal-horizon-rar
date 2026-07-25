"""
STAGE 7I-SW: MATERIAL-branch investigation, part 3. Part 2 (FPM profile)
REFUTED the noise-inflation hypothesis under its locked bar: the strict
sample prefers FPM = 1.5 in 4/4 rows and alpha_hat stays collapsed at
every FPM (H1 SUPPORTED printed, as pre-committed). The surviving
instrument flag: w_rad rides its GRID EDGE (0.3, the maximum) in all six
strict fits, versus 0.2 interior in 12/12 full-sample fits — the fit
saturates the radial-orbit dial. The alpha <-> w_rad degeneracy is what
the gamma direction channel normally breaks (4N: "the parameter
protector"); the strict cut thins the wide gamma statistics to 776/131
pairs.

TEST (bars locked here before execution): compose the two existing
variants — strict mask + w_rad FROZEN to the externally validated 0.21
(Hwang-implied; the W variant proved freezing it is free on the full
sample: d_alpha <= 0.01). FPM stays 1.5 (the value the strict data
themselves prefer, part 2). Seeds 31/101, both laws; auto-extend to
202/303/404/505 if ambiguous.
  DEGENERACY-RESOLVED if alpha_hat >= 0.9 both laws (interior) and
     Newton margin >= +30 (about half the N-scaled expectation ~ +55):
     the collapse is a nuisance-degeneracy artifact of the thinned
     sample, cured by the external eccentricity prior; the strict row
     then reads "alpha stable under the strictest cleaning once the
     externally validated w_rad prior is applied."
  COMPANION-DIRECTION if alpha_hat <= 0.8 or Newton < +20 with wr
     frozen: the low-boost preference is not the radial dial; the
     hidden-multiplicity reading gains real support and stage 7J
     (joint completeness marginalization) becomes the decisive
     instrument. Stated as the honest outcome, not absorbed.
  MIXED otherwise -> extend seeds; if still mixed, carried to 7J.
Writes data/stage7i_sw.txt + appends verdict to data/stage7i_verdict.txt.
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

WR_OLD = "WR_GRID = np.array([0.10, 0.20, 0.30])"
WR_NEW = "WR_GRID = np.array([0.21])"
OUT_OLD = "with open('data/stage3u_summary.txt', 'a') as f:"

def make_src():
    s = src0
    assert s.count(OLD_DATA) == 1
    s = s.replace(OLD_DATA, NEW_DATA)
    assert s.count(OK_OLD) == 1
    s = s.replace(OK_OLD, OK_NEW)
    assert s.count(WR_OLD) == 1
    s = s.replace(WR_OLD, WR_NEW)
    assert s.count(OUT_OLD) == 1
    s = s.replace(OUT_OLD, "with open('data/stage7i_sw.txt', 'a') as f:")
    return s

ROW = re.compile(r"seed (\d+) (simple|BE): a_hat=([0-9.]+) \(grid ([0-9.]+), "
                 r"interior=(\w+)\), dlnL\(Newton\)=([+-][0-9.]+), "
                 r"wr=([0-9.]+)")
def parse():
    out = {}
    if os.path.exists('data/stage7i_sw.txt'):
        for m in ROW.finditer(open('data/stage7i_sw.txt').read()):
            s_, law, ah, _, inter, dn, _wr = m.groups()
            out[(law, int(s_))] = dict(a=float(ah),
                                       interior=(inter == 'True'),
                                       dn=float(dn))
    return out

def run(seeds):
    have = parse()
    todo = [s_ for s_ in seeds if ('simple', s_) not in have]
    if todo:
        print(f"===== 7I-SW (strict + wr=0.21): seeds {todo} =====",
              flush=True)
        sys.argv = ['stage7i_sw', '1p2'] + [str(s_) for s_ in todo]
        exec(compile(make_src(), 'stage3p_7i_sw', 'exec'),
             {'__name__': '__main__'})
    return parse()

def judge(rec, seeds):
    a = {law: np.mean([rec[(law, s_)]['a'] for s_ in seeds])
         for law in ('simple', 'BE')}
    dn = {law: np.mean([rec[(law, s_)]['dn'] for s_ in seeds])
          for law in ('simple', 'BE')}
    n_int = sum(rec[(law, s_)]['interior'] for law in ('simple', 'BE')
                for s_ in seeds)
    amin, dmin = min(a.values()), min(dn.values())
    if amin >= 0.9 and dmin >= 30 and n_int == 2*len(seeds):
        v = 'DEGENERACY-RESOLVED'
    elif amin <= 0.8 or dmin < 20:
        v = 'COMPANION-DIRECTION'
    else:
        v = 'MIXED'
    return a, dn, n_int, v

SEED0 = [31, 101]
SEEDX = [202, 303, 404, 505]
rec = run(SEED0)
a, dn, n_int, v = judge(rec, SEED0)
seeds = SEED0
if v == 'MIXED':
    rec = run(SEED0+SEEDX)
    seeds = SEED0+SEEDX
    a, dn, n_int, v = judge(rec, seeds)
V = ["", f"SW: strict + wr frozen 0.21, FPM 1.5 ({len(seeds)} seeds):"]
for law in ('simple', 'BE'):
    per = " ".join(f"{rec[(law, s_)]['a']:.2f}" for s_ in seeds)
    V.append(f"  {law:>6}: a_hat mean {a[law]:.2f} (per-seed {per}), "
             f"Newton {dn[law]:+.1f}")
V.append(f"  interior {n_int}/{2*len(seeds)}  ==> {v} "
         f"(bars in stage7i_sw.py docstring)")
out = "\n".join(V)
print(out)
with open('data/stage7i_verdict.txt', 'a') as f:
    f.write(out + "\n")
print("\nSTAGE 7I-SW done -> appended data/stage7i_verdict.txt")
