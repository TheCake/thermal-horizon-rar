"""
ROUND-50 GA (blind half, committed BEFORE the review report exists;
17th execution of the 87a4676 protocol). Verifies the 10X battery's
as-fired record by independent routes and adds the C15 box-proximity
print the battery omitted (standing harness rule: profile-check
box-pinned parameters before quoting any variant).

GA-1  A4 bulge-free subset independently rebuilt from the rotmod files
      (bypasses the stage's world object), fit re-run COLD (no warm
      start) and WARM; full parameter vector + box proximity printed
      (f in (0.3, 2.5), s_int in [1e-3, 0.4), per-galaxy dml +-0.7
      pinned counts, both census conventions).
GA-2  A1 top influence (NGC5907) re-fit COLD; box proximity printed.
GA-3  census leg re-verified with astropy separations + independent
      modulus->distance arithmetic on the named CF4 rows.
GA-4  bar arithmetic recomputed from counts.
GA-5  A5 construction check (sv scaling reached the fit) + box print.
Output: data/round50_ga.txt
"""
import json, math, os, re, sys, time
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

t0 = time.time()
P("ROUND-50 GA (blind, pre-report)")

# inherit the engine exactly as the stage does
src = open('calcs/stage10t_legregrow.py', encoding='utf-8').read()
header = src[:src.index("# ---------------- G10T-1")]
ns = {'__name__': 'stage10t_inherited',
      '__file__': os.path.abspath('calcs/stage10t_legregrow.py')}
_argv = sys.argv
sys.argv = ['stage10t_legregrow.py', 'gates']
exec(compile(header, 'stage10t_legregrow.py[header]', 'exec'), ns)
sys.argv = _argv
make_world = ns['make_world']; build_sub = ns['build_sub']
fit_deep = ns['fit_deep']; fit_hier = ns['fit_hier']
nu_be = ns['nu_be']; h0_of_a0 = ns['h0_of_a0']
meta = ns['meta']; LN10 = ns['LN10']

W = make_world(lt_slot=None)
legA = [int(g) for g in W['legA']]
subT = build_sub(W, legA)
b_plain = fit_hier(subT, nu_be, use_u=False)
b_loose = fit_hier(subT, nu_be, use_u=True, th0=list(b_plain.x) + [0.0])
bT = fit_deep(subT, nu_be, th0=list(b_loose.x))
H0_T = h0_of_a0(10**bT.x[0])
P(f"primary re-derived: H0 = {H0_T:.2f} (battery printed 65.36)")

def boxprint(tag, sub, b):
    """C15: parameter vector + box proximity + dml pinned census."""
    la0, f, s_int, u = b.x
    P(f"  {tag}: la0 = {la0:.4f}, f_ML = {f:.3f} (box 0.3/2.5), "
      f"s_int = {s_int:.4f} (box 1e-3/0.4), u = {u:+.4f} (box +-0.2)")
    near = []
    if min(f - 0.3, 2.5 - f) < 0.05: near.append('f_ML')
    if s_int < 2e-3 or s_int > 0.38: near.append('s_int')
    if 0.2 - abs(u) < 0.01: near.append('u')
    # per-galaxy dml at the optimum: re-derive by one profiling sweep
    n = sub['n']
    lg = sub['lg']; gg, gd, gb_ = sub['gg'], sub['gd'], sub['gb']
    s2 = sub['s2']; sv = sub['sv']
    ig = sub['isuma_g']
    from scipy.optimize import minimize_scalar
    a0 = 10**la0
    se2c = s_int * s_int
    # dv from the closed form given dml = 0 start, then dml profile
    dml = np.zeros(n); dv = np.zeros(n)
    for _ in range(3):
        fac = f * np.exp(dml[sub['gidx']])
        gN = gg + fac * gd + gb_
        r0 = lg - np.log10(gN * nu_be(gN / a0)) - u * sub['isuma_pt']
        for gi2 in range(n):
            mm = sub['gpts'][gi2]
            w = 1.0 / (s2[mm] + se2c)
            dv[gi2] = np.sum(w * r0[mm]) / (np.sum(w) + 1.0 / sv[gi2]**2)
        for gi2 in range(n):
            mm = sub['gpts'][gi2]
            uoff = u * ig[gi2]
            def od(dl):
                fc = f * math.exp(dl)
                gN2 = gg[mm] + fc * gd[mm] + gb_[mm]
                rr = (lg[mm] - np.log10(gN2 * nu_be(gN2 / a0))
                      - dv[gi2] - uoff)
                ss = s2[mm] + se2c
                return np.sum(rr * rr / ss) + dl * dl / (0.1 * LN10)**2
            dml[gi2] = minimize_scalar(od, bounds=(-0.7, 0.7),
                                       method='bounded').x
    pin = np.abs(np.abs(dml) - 0.7) < 0.01
    P(f"    dml pinned (|dml| ~ 0.7): {int(pin.sum())}/{n} per-member; "
      f"largest |dml| = {np.max(np.abs(dml)):.3f}; "
      f"box-adjacent globals: {near if near else 'NONE'}")
    return near

P("")
P("-- GA-1 A4 bulge-free: independent subset + cold/warm refit --")
# independent membership: rotmod files, Vbul column, same survival cuts
KPC = ns['KPC']; UD, UB = ns['UD'], ns['UB']
indep = []
for g in legA:
    nm = W['names'][g]
    path = None
    import glob as _g
    for c in _g.glob(f'data/sparc/rotmod/**/{nm}_rotmod.dat',
                     recursive=True):
        path = c
    assert path, nm
    has_bul = False
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV / Vo > 0.10: continue
        gg_ = Vg * abs(Vg) / R * KPC
        gd_ = UD * Vd * abs(Vd) / R * KPC
        gb_ = UB * Vb * Vb / R * KPC
        if gg_ + gd_ + gb_ <= 0: continue
        if gb_ > 0: has_bul = True
    if not has_bul:
        indep.append(g)
stage_bf = [g for g in legA
            if bool(np.all(W['gb'][W['gpts'][g]] == 0.0))]
P(f"  independent bulge-free N = {len(indep)}; stage construction N = "
  f"{len(stage_bf)}; sets {'IDENTICAL' if set(indep) == set(stage_bf) else 'DIFFER: ' + str(set(indep) ^ set(stage_bf))}")
sub4 = build_sub(W, sorted(stage_bf))
b4w = fit_deep(sub4, nu_be, th0=list(bT.x))
b4c_pl = fit_hier(sub4, nu_be, use_u=False)
b4c_lo = fit_hier(sub4, nu_be, use_u=True, th0=list(b4c_pl.x) + [0.0])
b4c = fit_deep(sub4, nu_be, th0=list(b4c_lo.x))
P(f"  warm H0 = {h0_of_a0(10**b4w.x[0]):.2f}; COLD-chain H0 = "
  f"{h0_of_a0(10**b4c.x[0]):.2f} (battery printed 59.49); "
  f"-2lnL warm {b4w.fun:.2f} vs cold {b4c.fun:.2f}")
best4 = b4c if b4c.fun < b4w.fun else b4w
near4 = boxprint("A4 optimum", sub4, best4)

P("")
P("-- GA-2 A1 NGC5907 leave-out: cold refit --")
gi5907 = ns['SPNAME2GI']['NGC5907']
sub_l = build_sub(W, [x for x in legA if x != gi5907])
blw = fit_deep(sub_l, nu_be, th0=list(bT.x))
blc_pl = fit_hier(sub_l, nu_be, use_u=False)
blc_lo = fit_hier(sub_l, nu_be, use_u=True, th0=list(blc_pl.x) + [0.0])
blc = fit_deep(sub_l, nu_be, th0=list(blc_lo.x))
P(f"  warm dH0 = {h0_of_a0(10**blw.x[0]) - H0_T:+.3f}; cold-chain dH0 = "
  f"{h0_of_a0(10**blc.x[0]) - H0_T:+.3f} (battery printed +2.250); "
  f"-2lnL warm {blw.fun:.2f} vs cold {blc.fun:.2f}")
bestl = blc if blc.fun < blw.fun else blw
boxprint("LOO-5907 optimum", sub_l, bestl)

P("")
P("-- GA-3 census leg: astropy separations + modulus arithmetic --")
from astropy.coordinates import SkyCoord
import astropy.units as uu
pos = json.load(open('data/stage10x_pgc.json'))
rows = {}
for l in open('data/cf4/table2.dat'):
    try:
        pgc = int(l[0:7])
    except ValueError:
        continue
    if pgc in (86668, 22277):
        rows[pgc] = l
for nm, pgc_expect in (('D564-8', 86668), ('D631-7', 22277)):
    l = rows[pgc_expect]
    ra, dec = float(l[137:145]), float(l[146:154])
    c1 = SkyCoord(pos[nm]['ra'] * uu.deg, pos[nm]['dec'] * uu.deg)
    c2 = SkyCoord(ra * uu.deg, dec * uu.deg)
    sep = c1.separation(c2).arcmin
    dm, edm = float(l[102:107]), float(l[108:112])
    D = 10**((dm - 25) / 5)
    eD = D * (LN10 / 5) * edm
    gi = ns['SPNAME2GI'][nm]
    D_sp, eD_sp = ns['SP']['dist'][gi]
    P(f"  {nm} vs PGC {pgc_expect}: sep = {sep:.2f} arcmin (astropy); "
      f"DMtrgb = {dm} +- {edm} -> D = {D:.2f} +- {eD:.2f} Mpc; SPARC "
      f"{D_sp:.2f}; gap {D/D_sp-1:+.1%} (battery: -1.8%/-0.4%)")
# KK98-251: nearest CF4 row distance (full scan)
best, bd = None, 1e9
kra, kdec = pos['KK98-251']['ra'], pos['KK98-251']['dec']
cd = math.cos(math.radians(kdec))
for l in open('data/cf4/table2.dat'):
    try:
        pgc = int(l[0:7])
    except ValueError:
        continue
    s1, s2 = l[137:145].strip(), l[146:154].strip()
    if not (s1 and s2): continue
    dd = math.hypot((float(s1) - kra) * cd, float(s2) - kdec)
    if dd < bd: best, bd = pgc, dd
P(f"  KK98-251: nearest CF4 row anywhere = PGC {best} at "
  f"{bd*60:.1f} arcmin (battery rule: no row within 3 arcmin -> "
  f"{'CONFIRMED' if bd*60 > 3 else 'CONTRADICTED'})")

P("")
P("-- GA-4 bar arithmetic --")
SIG, NF = 4.99, 78
for nm, nsub in (('LOO', 77), ('incl40', 74), ('incl45', 71),
                 ('incl50', 63), ('Q1', 46), ('bulge-free', 66)):
    bar = max(2.0 * SIG * math.sqrt(NF / nsub - 1.0), 1.0)
    P(f"  {nm}: N = {nsub} -> bar = {bar:.2f}")
P("  A5: floor 1.00 (same membership)")

P("")
P("-- GA-5 A5 construction --")
sub5 = build_sub(W, legA)
sv0 = sub5['sv'].copy()
sub5['sv'] = sub5['sv'] * 1.5
P(f"  sv ratio applied: min/max = {np.min(sub5['sv']/sv0):.3f}/"
  f"{np.max(sub5['sv']/sv0):.3f} (expect 1.5/1.5)")
b5 = fit_deep(sub5, nu_be, th0=list(bT.x))
P(f"  refit H0 = {h0_of_a0(10**b5.x[0]):.2f} (battery printed 66.54)")
boxprint("A5 optimum", sub5, b5)

P("")
P(f"GA total: {(time.time()-t0)/60:.1f} min")
with open('data/round50_ga.txt', 'w', encoding='utf-8') as f:
    f.write("\n".join(L) + "\n")
