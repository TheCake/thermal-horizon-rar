"""
ROUND-31 verification addendum -- BLIND HALF, committed BEFORE the
referee's report exists (the 87a4676 protocol, second execution).

Re-derives every load-bearing 10K/10L number through routes INDEPENDENT
of the stage code: different root-finders, different parsers, different
series extractions.  GA-1..7 = 10L; GA-8..11 = 10K.  The post-report
half (his numbers) is appended below the marker after the round lands.
"""
import csv, glob, math, os
import numpy as np
from scipy.optimize import brentq, minimize_scalar
import sympy as sp

OUT = 'data/round31_addendum.txt'
L = []
def emit(s=""):
    L.append(s)
    print(s, flush=True)
    with open(OUT, 'w') as f:
        f.write("\n".join(L) + "\n")

emit("ROUND-31 ADDENDUM, BLIND HALF (pre-report)")
emit("")

def nbe(z):
    if z > 700: return 0.0
    return 1.0/math.expm1(z)

# ---------- GA-1: 10L existence boundary via brentq (3rd detector) ----
def exists_brentq(kv, xv, c=1.0):
    nstar = 2.0/kv - c
    if nstar <= 0: return False
    f = lambda nb: nbe(xv*(1.0 - kv*(nb + c)/2.0)) - nb
    grid = np.linspace(0.0, nstar*(1-1e-9), 200)
    vals = [f(nb) for nb in grid]
    for a, b in zip(range(199), range(1, 200)):
        if vals[a] > 0 and vals[b] <= 0:
            brentq(f, grid[a], grid[b])
            return True
    return False

def kmax_brentq(xv, c=1.0):
    lo, hi = 1e-4, 4.0
    if not exists_brentq(lo, xv, c): return 0.0
    for _ in range(50):
        mid = 0.5*(lo+hi)
        if exists_brentq(mid, xv, c): lo = mid
        else: hi = mid
    return 0.5*(lo+hi)

STAGE_KMAX = {0.03: 0.0149, 0.05: 0.0247, 0.1: 0.0488,
              0.2: 0.0952, 0.3: 0.1395}
ok1 = True
for xv, want in STAGE_KMAX.items():
    got = kmax_brentq(xv)
    ok1 &= abs(got - want) <= 0.005
    emit("GA-1 kappa_max(x=%.2f) brentq: %.4f vs stage %.4f  %s" %
         (xv, got, want, "OK" if abs(got-want) <= 0.005 else "MISMATCH"))

# ---------- GA-2: x-onset at kappa_lock ----------
KLOCK = 0.925
f_on = lambda xv: 1.0 if exists_brentq(KLOCK, xv) else -1.0
lo, hi = 2.0, 5.0
for _ in range(40):
    mid = 0.5*(lo+hi)
    if f_on(mid) > 0: hi = mid
    else: lo = mid
x_on = 0.5*(lo+hi)
ok2 = abs(x_on - 3.56) <= 0.05
emit("GA-2 existence onset at kappa_lock: x = %.3f vs stage ~3.56  %s" %
     (x_on, "OK" if ok2 else "MISMATCH"))

# ---------- GA-3: truncated-Gibbs via log-sum-exp (independent) -------
def tg(kv, xv, keep_peak):
    E = []
    nn = 0
    while True:
        Ev = (nn+0.5) - (kv/4.0)*(nn+0.5)**2
        if nn > 0 and Ev <= E[-1]:
            break
        E.append(Ev)
        nn += 1
    if not keep_peak and len(E) > 1:
        E = E[:-1]
    lw = -xv*np.array(E)
    lw -= lw.max()
    w = np.exp(lw)
    return float((np.arange(len(E))*w).sum()/w.sum())

rows = [(True, 1.916, 1.925), (False, 1.456, 1.462)]
ok3 = True
for kp, want05, wantsup in rows:
    v05 = 1 + KLOCK*tg(KLOCK, 0.05, kp)
    sup = max(1 + KLOCK*tg(KLOCK, float(x), kp)
              for x in np.geomspace(1e-3, 10, 200))
    o = abs(v05 - want05) <= 0.002 and abs(sup - wantsup) <= 0.002
    ok3 &= o
    emit("GA-3 truncated-Gibbs keep_peak=%s: nu(0.05) = %.3f (stage "
         "%.3f), sup = %.3f (stage %.3f)  %s" %
         (kp, v05, want05, sup, wantsup, "OK" if o else "MISMATCH"))

# ---------- GA-4: classical tangency closed form (sympy) ----------
nb_, kap_, x_, s_ = sp.symbols('nbar kappa x s', positive=True)
s_expr = 1 - kap_*(nb_ + 1)/2
G_cl = 1/(x_*s_expr) - sp.Rational(1, 2)
eq1 = sp.Eq(G_cl, nb_)
eq2 = sp.Eq(sp.diff(G_cl, nb_), 1)
sol = sp.solve([eq1, eq2], [nb_, x_], dict=True)
xc_sym = None
for so in sol:
    xc = sp.simplify(so[x_])
    sst = sp.simplify(1 - kap_*(so[nb_] + 1)/2)
    if sp.simplify(sst - (sp.Rational(1, 2) - kap_/8)) == 0:
        xc_sym = xc
ok4 = xc_sym is not None and \
    sp.simplify(xc_sym - kap_/(2*(sp.Rational(1, 2) - kap_/8)**2)) == 0
emit("GA-4 classical tangency: s* = 1/2 - kappa/8 and x_c = "
     "kappa/(2(1/2 - kappa/8)^2) %s (sympy; value at kappa=1: %.3f "
     "vs stage 3.56)" %
     ("CONFIRMED" if ok4 else "NOT CONFIRMED",
      float(xc_sym.subs(kap_, 1)) if xc_sym is not None else float('nan')))

# ---------- GA-5: c1(beta) numerically (independent of the series) ----
def nu_stiff(xv, bv):
    nu = 1.0/xv
    for _ in range(600):
        u = xv**(1+bv)*nu**bv
        nu2 = 1 + nbe(u)
        if abs(nu2-nu) < 1e-15*nu: break
        nu = 0.5*(nu+nu2)
    return nu
ok5 = True
for bv, want in ((0.0, 0.5), (0.5, 1/3.0), (1.0, 0.25)):
    c1a = nu_stiff(1e-4, bv) - 1e4
    c1b = nu_stiff(1e-5, bv) - 1e5
    c1_extr = c1b + (c1b - c1a)/9.0
    o = abs(c1_extr - want) <= 1e-4
    ok5 &= o
    emit("GA-5 c1(beta=%.1f) numeric-extrapolated: %.6f vs %.6f  %s" %
         (bv, c1_extr, want, "OK" if o else "MISMATCH"))

# ---------- GA-6: kappa_r map ----------
ok6 = True
for kb, want in ((0.888, 0.494), (0.925, 0.497), (1.0, 0.5),
                 (1.10, 0.495), (1.48, 0.385), (1.503, 0.373)):
    v = kb*(1-kb/2)
    o = abs(v - want) <= 5e-4
    ok6 &= o
kk = np.linspace(0.2, 2.0, 100001)
kmax_at = kk[np.argmax(kk*(1-kk/2))]
ok6 &= abs(kmax_at - 1.0) <= 1e-4
emit("GA-6 kappa_r table + argmax at kappa_b = %.5f (want 1)  %s" %
     (kmax_at, "OK" if ok6 else "MISMATCH"))

# ---------- GA-7: the deep demand ----------
d1 = nbe(0.05)
ok7 = abs(d1 - 19.5042) <= 1e-3 and abs(1 + KLOCK*d1 - 19.04) <= 0.01
emit("GA-7 n_BE(0.05) = %.4f; demand 1 + 0.925 n_BE = %.2f  %s" %
     (d1, 1 + KLOCK*d1, "OK" if ok7 else "MISMATCH"))
emit("")

# ---------- GA-8: 10K census via an independent parser ----------
def canon2(nm):
    s = "".join(ch for ch in str(nm).upper() if ch not in " -_")
    for p in ('NGC', 'UGC', 'IC', 'DDO', 'ESO', 'PGC', 'UGCA'):
        if s.startswith(p):
            t = s[len(p):].lstrip('0')
            if t and t[0].isdigit() or t:
                return p + t if t else s
    return s

KPC = 3.24078e-14
meta2 = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    ls = f.readlines()
st = max(i for i, l in enumerate(ls) if set(l.strip()) <= set('- ')) + 1
for l in ls[st:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta2[t[0]] = (float(t[5]), int(t[17]), float(t[2]))
    except ValueError:
        continue
snames, sD = [], []
strans = {}
for path in sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                             recursive=True)):
    nm = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D = meta2.get(nm, (0, 3, 10.0))
    if inc < 30 or q > 2: continue
    npts = 0
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gN = (Vg*abs(Vg) + 0.5*Vd*abs(Vd) + 0.7*Vb*Vb)/R*KPC
        if gN <= 0: continue
        if gN/1.2e-10 >= 0.8: npts += 1
    snames.append(nm); sD.append(D); strans[nm] = npts >= 3
ung = []
with open('data/karachentsev_ungc.csv', encoding='utf-8',
          errors='replace') as f:
    for row in csv.DictReader(f):
        try:
            ung.append((row['Name'].strip(),
                        row.get('SimbadName', '').strip(),
                        row.get('NEDname', '').strip(),
                        float(row['RAJ2000']), float(row['DEJ2000']),
                        float(row['Dist']),
                        float(row['Ti1']) if row.get('Ti1', '') not in
                        ('', None) else float('nan'),
                        float(row['KLum']) if row.get('KLum', '') not in
                        ('', None) else float('nan')))
        except (ValueError, KeyError):
            continue
uc = {}
for k, u in enumerate(ung):
    for nm in u[:3]:
        if nm: uc.setdefault(canon2(nm), k)
AL = {}
with open('data/sparc_simbad_aliases.csv') as f:
    rd = csv.reader(l for l in f if not l.startswith('#'))
    next(rd, None)
    for row in rd:
        if len(row) >= 2: AL.setdefault(row[0], []).append(row[1])
nm_, ni_, nt_, ng_, ng1_ = 0, 0, 0, 0, 0
for nm, D in zip(snames, sD):
    hit = -1
    for c in [canon2(nm)] + [canon2(a) for a in AL.get(nm, [])]:
        k = uc.get(c, -1)
        if k >= 0 and abs(D - ung[k][5]) <= max(2.0, 0.35*ung[k][5]):
            hit = k
            break
    if hit >= 0:
        nm_ += 1
        ti = ung[hit][6]
        if ti < 0:
            ni_ += 1
            if strans[nm]: nt_ += 1
        if ti > 0: ng_ += 1
        if ti >= 1: ng1_ += 1
ok8 = (nm_, ni_, nt_, ng_, ng1_) == (48, 20, 3, 27, 13)
emit("GA-8 independent census parse: matched %d iso %d iso-TR %d "
     "group %d group>=1 %d vs stage 48/20/3/27/13  %s" %
     (nm_, ni_, nt_, ng_, ng1_, "OK" if ok8 else "MISMATCH"))

# ---------- GA-9: NGC2976 -- not in Chae; neighbor term ----------
chae_names = set()
with open('data/chae2021_table3.csv') as f:
    for row in csv.DictReader(l for l in f if not l.startswith('#')):
        chae_names.add(row['galaxy'])
in_chae = 'NGC2976' in chae_names
k29 = next(k for k, u in enumerate(ung) if canon2(u[0]) == 'NGC2976')
G_SI, MSUN, MPC = 6.674e-11, 1.989e30, 3.0857e22
def xyz(u):
    ra, de, d = math.radians(u[3]), math.radians(u[4]), u[5]
    return (d*math.cos(de)*math.cos(ra), d*math.cos(de)*math.sin(ra),
            d*math.sin(de))
p0 = xyz(ung[k29])
ssum = 0.0
for k, u in enumerate(ung):
    if k == k29 or not np.isfinite(u[7]): continue
    p = xyz(u)
    d2 = sum((a-b)**2 for a, b in zip(p0, p))
    if d2 <= 9.0:
        ssum += G_SI*(10.0**u[7])*MSUN/(d2*MPC**2)
frac = ssum/1.2e-10/0.01*100
ok9 = (not in_chae) and abs(frac - 109.9) <= 1.0
emit("GA-9 NGC2976: in Chae table = %s (stage: NO); neighbor term = "
     "%.1f%% of the 0.01 floor (stage 109.9)  %s" %
     (in_chae, frac, "OK" if ok9 else "MISMATCH"))

# ---------- GA-10: stratum arithmetic ----------
def s_of(e):
    n = 1.0/math.expm1(math.sqrt(e))
    return n/(1.0+n)
dp = 0.3365*(-0.00352)/2
dpc = 0.50*(-0.00352)/2
ok10 = abs(dp - (-0.000593)) <= 2e-6 and abs(dp/0.075 - (-0.0079)) <= 2e-3 \
    and abs(dpc - (-0.000881)) <= 3e-6
emit("GA-10 stratum arithmetic: Dp = %.6f (stage -0.000593) = %.4f "
     "sigma_p; co-read %.6f (stage -0.000881)  %s" %
     (dp, dp/0.075, dpc, "OK" if ok10 else "MISMATCH"))

# ---------- GA-11: Fisher consistency + nuisance-free upper bound -----
C_arch = 2.0/0.073**2
ok11a = abs(C_arch - 375.3) <= 0.1 and abs(0.371*1010.7 - C_arch) <= 1.0
# nuisance-free upper bound on the contest signal: no dv-profiling,
# sigma^2 >= sig2 only, all 24 affected galaxies at their largest Dp
def nu_p(y, p):
    yc = max(y, 1e-14)
    return (1.0 - math.exp(-min(yc**p, 60.0)))**(-1.0/(2.0*p))
DP_MAX = 0.00341        # largest per-galaxy Dp in the stage table
S_ub = 0.0
npts_aff = 0
for path in sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                             recursive=True)):
    nm = os.path.basename(path).replace('_rotmod.dat', '')
    if nm not in ('UGC05986', 'UGC07323', 'UGC07399', 'UGC08837',
                  'UGC08550', 'UGC07151', 'NGC2903', 'DDO168'):
        continue
    inc, q, D = meta2.get(nm, (0, 3, 10.0))
    if inc < 30 or q > 2: continue
    for l in open(path):
        if l.startswith('#'): continue
        t = l.split()
        if len(t) < 6: continue
        R, Vo, eV, Vg, Vd, Vb = map(float, t[:6])
        if R <= 0 or Vo <= 0 or eV/Vo > 0.10: continue
        gN = (Vg*abs(Vg) + 0.5*Vd*abs(Vd) + 0.7*Vb*Vb)/R*KPC
        if gN <= 0: continue
        y = gN/1.2e-10
        p0v = 0.5 + 0.3365*0.85/2
        dmu = abs(math.log10(nu_p(y, p0v + DP_MAX))
                  - math.log10(nu_p(y, p0v)))
        sg = 2*eV/Vo/math.log(10)
        S_ub += (dmu/sg)**2
        npts_aff += 1
# the 16 remaining affected galaxies have smaller Dp than any of these
# 8; bound their total by 16/8 times this sum (conservative)
S_ub_tot = S_ub*3.0
ok11b = S_ub_tot < 4.0
emit("GA-11 Fisher: C_arch = %.1f, D_f x C_fish = %.1f (consistency "
     "OK=%s); nuisance-FREE upper bound on the contest signal (no "
     "profiling, formal errors only, top-8 galaxies x3): S <= %.3f "
     "(%d pts) < bar 4  %s" %
     (C_arch, 0.371*1010.7, ok11a, S_ub_tot, npts_aff,
      "OK" if ok11b else "MISMATCH"))
emit("")
allok = all((ok1, ok2, ok3, ok4, ok5, ok6, ok7, ok8, ok9, ok10,
             ok11a, ok11b))
emit("BLIND HALF: %s (GA-1..11)" %
     ("ALL CONFIRMED" if allok else "MISMATCHES PRESENT"))
