"""STAGE 9N — THE RISING-FLAG DIAL EXPOSURE (round-13-addendum
short check).  Pre-registered BEFORE any run (measurement round;
NO credence movement).

Does the headline deep-regime dial (8S-c FULL vertON lam_hat =
0.960 => c1_hat = 0.480) move when rising-curve galaxies (9M
frozen flag: s_out >= 0.05 over the last 3 kept points) are
excluded?  Engine = the 9M verbatim lift (m2ll_fast/lam_hat_fast).
Gates: G9N-0 verbatim-lift (counts + OFF probe 1e-6); G9N-1
fast-vs-verbatim 3 probes 1e-6; G9N-2 lineage regressions
(same-engine 9M printed at 0.002; cross-instrument 8S-c printed
dials at 0.05, looser bar disclosed).  Paired bootstrap NBOOT=200
rng 83 on D = lam(rep-RISING) - lam(rep), skip if either side
< 30.  Bars (ordered): N-MOVED / N-ROBUST / N-GRAY-CARRIED.
Output: data/stage9n_risingdial.txt
"""
import glob, math, os, re, time
import numpy as np
from scipy.optimize import minimize

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
        q = int(t[17])
        meta[name] = (inc, q, D, eD, einc)
    except ValueError:
        continue

g_gas, g_dsk, g_bul, gobs, sig, gal_id = [], [], [], [], [], []
svert = {}
vkept = {}
gals = sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat', recursive=True))
for gi, path in enumerate(gals):
    name = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D, eD, einc = meta.get(name, (0, 3, 1.0, 0.0, 0.0))
    if inc < 30 or q > 2: continue
    s_d2 = (eD/max(D, 1e-6))/LN10
    s_i2 = 2.0*math.radians(einc)/max(math.tan(math.radians(inc)), 1e-6)/LN10
    svert[gi] = math.sqrt(s_d2**2 + s_i2**2)
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
        gal_id.append(gi)
        vkept.setdefault(gi, []).append(Vo)
g_gas, g_dsk, g_bul, gobs, sig, gal_id = map(
    np.array, (g_gas, g_dsk, g_bul, gobs, sig, gal_id))
sig2 = sig*sig
lgobs = np.log10(gobs)

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
P("9N THE RISING-FLAG DIAL EXPOSURE (pre-reg committed BEFORE any "
  "run; measurement round; NO credence movement)")

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

txt8sc = open('data/stage8sc_gddist.txt').read()
OFFV = float(re.search(r"G8Sc-0 OFF-branch identity: ([-\d.]+) vs",
                       txt8sc).group(1))
DIAL8 = {}
for tag in ('FULL', 'GD', 'DD'):
    DIAL8[tag] = float(re.search(
        r"\[" + tag + r" +vertON\] lam_hat = ([-\d.]+) ",
        txt8sc).group(1))

txt9m = open('data/stage9m_convsplit.txt').read()
mm = re.search(r"\[point\] lam_hat\(CONVERGED\) = ([-+\d.]+) \(edges"
               r" \d/\d\), lam_hat\(RISING\) = ([-+\d.]+)", txt9m)
GD9M = (float(mm.group(1)), float(mm.group(2)))
mm = re.search(r"G9M-2 flag census \(GD-38\): CONVERGED (\d+), "
               r"RISING (\d+), AMBIG (\d+)", txt9m)
CEN9M = (int(mm.group(1)), int(mm.group(2)), int(mm.group(3)))
mm = re.search(r"\[DD control\] CONVERGED (\d+), RISING (\d+)", txt9m)
CEND = (int(mm.group(1)), int(mm.group(2)))
mm = re.search(r"\[DD control\] lam_hat\(CONV\) = ([-+\d.]+), "
               r"lam_hat\(RISING\) = ([-+\d.]+)", txt9m)
DD9M = (float(mm.group(1)), float(mm.group(2)))

def m2ll_vert(th, lam, gset, von=True):
    la0, f, s_int = th
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    nv = nu_lam(gN/a0, lam)
    if np.min(nv) <= 0.05: return 1e12
    gm = gN*nv
    v = sig2 + s_int*s_int
    a = lgobs - np.log10(gm)
    out = 0.0
    for g_ in gset:
        ji = GIDX[g_]
        if len(ji) == 0: continue
        aj, vj = a[ji], v[ji]
        s_ = SV[g_]
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
    return (np.concatenate(cats), np.concatenate(labs),
            np.array(svals), li)

def m2ll_fast(th, lam, cat, lab, svec, ninst):
    la0, f, s_int = th
    if not (-10.6 < la0 < -9.4) or not (0.3 <= f <= 2.5): return 1e12
    if not (1e-3 <= s_int < 0.4): return 1e12
    a0 = 10**la0
    gN = g_gas + f*g_dsk + g_bul
    nv = nu_lam(gN/a0, lam)
    if np.min(nv) <= 0.05: return 1e12
    gm = gN*nv
    v = sig2 + s_int*s_int
    a = lgobs - np.log10(gm)
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
def lam_hat_fast(cat, lab, svec, ninst):
    prof = []; th = None
    for lam in LGB:
        starts = (([list(th)] if th is not None else [])
                  + [[math.log10(A0_FID), 1.0, 0.08]])
        best = None
        for th0 in starts:
            b = minimize(lambda t: m2ll_fast(t, lam, cat, lab, svec,
                                             ninst), th0,
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
th_p = [math.log10(A0_FID), 1.0, 0.08]
ok_n = (len(gd_set) == 38) and (len(dd_set) == 111)
va = m2ll_vert(th_p, 0.7, list(gd_set), von=False)
ok_off = abs(va - OFFV) <= 1e-6
g0_ok = ok_n and ok_off
P(f"G9N-0 verbatim-lift: GD {len(gd_set)} / DD {len(dd_set)} -> "
  f"{'PASS' if ok_n else 'FAIL'}; OFF probe d = {va-OFFV:.2e} -> "
  f"{'PASS' if ok_off else 'FAIL'}")

rngp = np.random.default_rng(7)
g1_ok = True; dmx = 0.0
probe_sets = [[int(g_) for g_ in gd_set],
              [int(g_) for g_ in dd_set],
              [int(g_) for g_ in allg]]
for pi in range(3):
    thp = [-10.5+rngp.random(), 0.4+2.0*rngp.random(),
           0.01+0.29*rngp.random()]
    lam = -2.0+3.5*rngp.random()
    gl = probe_sets[pi]
    slow = m2ll_vert(thp, lam, gl, von=True)
    cat, lab, svec, ninst = make_instance(gl)
    fast = m2ll_fast(thp, lam, cat, lab, svec, ninst)
    if slow < 1e11:
        dmx = max(dmx, abs(slow-fast))
g1_ok = dmx <= 1e-6
P(f"G9N-1 fast-vs-verbatim (3 probes): max|d| = {dmx:.2e} -> "
  f"{'PASS' if g1_ok else 'FAIL'}")

# ---------------- frozen flags (9M rule, ALL fit galaxies) -------
flags = {}
for g_ in allg:
    V = vkept.get(int(g_), [])
    if len(V) < 3:
        flags[int(g_)] = 'AMBIG'; continue
    s_out = (V[-1] - V[-3])/max(V[-1], 1e-6)
    flags[int(g_)] = ('RISING' if s_out >= 0.05 else
                      'CONVERGED' if s_out <= 0.02 else 'AMBIG')
gd_l = [int(g_) for g_ in gd_set]
dd_l = [int(g_) for g_ in dd_set]
full_l = [int(g_) for g_ in allg]
cen_gd = (sum(1 for g_ in gd_l if flags[g_] == 'CONVERGED'),
          sum(1 for g_ in gd_l if flags[g_] == 'RISING'),
          sum(1 for g_ in gd_l if flags[g_] == 'AMBIG'))
cen_dd = (sum(1 for g_ in dd_l if flags[g_] == 'CONVERGED'),
          sum(1 for g_ in dd_l if flags[g_] == 'RISING'),
          sum(1 for g_ in dd_l if flags[g_] == 'AMBIG'))
ok_cen = (cen_gd == CEN9M) and (cen_dd[0] == CEND[0]
                                and cen_dd[1] == CEND[1])
P(f"G9N-2a census regression: GD {cen_gd} vs 9M {CEN9M}; DD "
  f"({cen_dd[0]}, {cen_dd[1]}) vs 9M {CEND} -> "
  f"{'PASS' if ok_cen else 'FAIL'}")

# same-engine regressions (0.002)
conv_gd = [g_ for g_ in gd_l if flags[g_] == 'CONVERGED']
rise_gd = [g_ for g_ in gd_l if flags[g_] == 'RISING']
conv_dd = [g_ for g_ in dd_l if flags[g_] == 'CONVERGED']
rise_dd = [g_ for g_ in dd_l if flags[g_] == 'RISING']
reg_ok = True
regs = []
for lab_, gl_, tgt in (('GD-CONV', conv_gd, GD9M[0]),
                       ('GD-RISE', rise_gd, GD9M[1]),
                       ('DD-CONV', conv_dd, DD9M[0]),
                       ('DD-RISE', rise_dd, DD9M[1])):
    lh, _, _ = lam_hat_fast(*make_instance(gl_))
    okr = abs(lh - tgt) <= 0.002
    reg_ok &= okr
    regs.append(f"{lab_} {lh:+.3f} vs {tgt:+.3f}")
P("G9N-2b same-engine regressions (0.002): "
  + "; ".join(regs) + f" -> {'PASS' if reg_ok else 'FAIL'}")

# cross-instrument regressions (0.05, disclosed looser bar)
lh_full, ef1, ef2 = lam_hat_fast(*make_instance(full_l))
lh_gd, _, _ = lam_hat_fast(*make_instance(gd_l))
lh_dd, _, _ = lam_hat_fast(*make_instance(dd_l))
x_ok = (abs(lh_full - DIAL8['FULL']) <= 0.05
        and abs(lh_gd - DIAL8['GD']) <= 0.05
        and abs(lh_dd - DIAL8['DD']) <= 0.05)
P(f"G9N-2c cross-instrument regressions (0.05): FULL {lh_full:+.3f}"
  f" vs {DIAL8['FULL']:+.3f}; GD {lh_gd:+.3f} vs {DIAL8['GD']:+.3f};"
  f" DD {lh_dd:+.3f} vs {DIAL8['DD']:+.3f} -> "
  f"{'PASS' if x_ok else 'FAIL'}")
g2_ok = ok_cen and reg_ok and x_ok

if not (g0_ok and g1_ok and g2_ok):
    P("GATES FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
    with open('data/stage9n_risingdial.txt', 'w') as f:
        f.write("\n".join(L)+"\n")
    raise SystemExit(0)
P("GATES: G9N-0/1/2 PASS")
P("")

# ---------------- point dials ----------------
norise = [g_ for g_ in full_l if flags[g_] != 'RISING']
convon = [g_ for g_ in full_l if flags[g_] == 'CONVERGED']
lh_nr, en1, en2 = lam_hat_fast(*make_instance(norise))
lh_co, ec1, ec2 = lam_hat_fast(*make_instance(convon))
Dp = lh_nr - lh_full
P(f"[point] lam_hat(FULL-149) = {lh_full:+.3f} (edges {int(ef1)}/"
  f"{int(ef2)}); lam_hat(FULL-RISING, n={len(norise)}) = "
  f"{lh_nr:+.3f} (edges {int(en1)}/{int(en2)}); D = {Dp:+.3f}")
P(f"[point] lam_hat(CONV-only, n={len(convon)}) = {lh_co:+.3f} "
  f"(edges {int(ec1)}/{int(ec2)}); c1_hat: full {lh_full/2:+.3f} "
  f"-> norise {lh_nr/2:+.3f}")
gd_nr = [g_ for g_ in gd_l if flags[g_] != 'RISING']
dd_nr = [g_ for g_ in dd_l if flags[g_] != 'RISING']
lh_gnr, _, _ = lam_hat_fast(*make_instance(gd_nr))
lh_dnr, _, _ = lam_hat_fast(*make_instance(dd_nr))
P(f"[co-read] GD full {lh_gd:+.3f} -> -RISING {lh_gnr:+.3f}; DD "
  f"full {lh_dd:+.3f} -> -RISING {lh_dnr:+.3f}; GAP full "
  f"{lh_gd-lh_dd:+.3f} -> -RISING {lh_gnr-lh_dnr:+.3f}")
P("")

# ---------------- paired bootstrap ----------------
NBOOT = 200
rng = np.random.default_rng(83)
Ds = []; skips = 0
for r in range(NBOOT):
    draw = rng.integers(0, len(full_l), size=len(full_l))
    rep = [full_l[j] for j in draw]
    if r == 0:
        P(f"G9N-3 first-replicate fingerprint: sum(draw) = "
          f"{int(np.sum(draw))}, unique = {len(np.unique(draw))}")
    nr_ = [g_ for g_ in rep if flags[g_] != 'RISING']
    if min(len(rep), len(nr_)) < 30:
        skips += 1; continue
    lf, _, _ = lam_hat_fast(*make_instance(rep))
    ln_, _, _ = lam_hat_fast(*make_instance(nr_))
    Ds.append(ln_ - lf)
    if (r+1) % 25 == 0:
        P(f"  ... {r+1}/{NBOOT} ({(time.time()-t00)/60:.1f} min)")
Ds = np.array(Ds)
qd = np.percentile(Ds, [5, 50, 95])
pge = float(np.mean(Ds >= 0))
pbig = float(np.mean(np.abs(Ds) >= 0.25))
P(f"[boot] reps {len(Ds)} (skips {skips}); D pct 5/50/95 = "
  f"{qd[0]:+.3f}/{qd[1]:+.3f}/{qd[2]:+.3f}; P(D >= 0) = {pge:.3f}; "
  f"P(|D| >= 0.25) = {pbig:.3f}")
P("")

moved = (abs(Dp) >= 0.25) and (qd[0] > 0 or qd[2] < 0)
robust = (abs(Dp) <= 0.125) and (pbig <= 0.10)
if moved:
    P("==> 9N VERDICT (locked grammar): N-MOVED - the headline "
      "deep-regime dial depends materially on rising-curve "
      "galaxies; the Paper-2 dial quote gains a rising-flag "
      "conditional column.")
elif robust:
    P("==> 9N VERDICT (locked grammar): N-ROBUST - the headline "
      "dial does not depend on rising-curve galaxies; one-sentence "
      "robustness row for the paper.")
else:
    P("==> 9N VERDICT (locked grammar): N-GRAY-CARRIED - rows "
      "stand as measurements.")
P("    NO credence movement (pre-stated; measurement round).")
P(f"\ndone ({(time.time()-t00)/60:.1f} min)")

with open('data/stage9n_risingdial.txt', 'w') as f:
    f.write("\n".join(L)+"\n")
print("\nsaved: data/stage9n_risingdial.txt")
