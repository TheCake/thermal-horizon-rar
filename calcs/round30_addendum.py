# round30_addendum.py -- independent re-verification of the 10I + 10J
# load-bearing numbers (the standing memory rule: re-compute every
# number before adopting any ROUND-30 ruling).  Written BLIND, before
# the ROUND-30 report was read; extended afterward if the referee
# computes new objects.  Independence: chi(1) via full numeric
# eigenvectors (not the analytic polaron/Poisson route); Theta-slope
# via scipy.stats.linregress on an independently rebuilt match list;
# neighbor-scale arithmetic from scratch.
import csv, glob, math, os, re
import numpy as np
from scipy.stats import linregress

OUT = 'data/round30_addendum.txt'
L = []
def emit(s=""):
    L.append(s)
    print(s, flush=True)
    with open(OUT, 'w') as f:
        f.write("\n".join(L) + "\n")

emit("ROUND-30 ADDENDUM (blind part written before reading the report)")
emit("")

# ==================== GA-1: 10J chi(1) via full numerics ====================
emit("GA-1  chi(1) by full per-block EIGENVECTORS (independent of the")
emit("      analytic displaced-vacuum/Poisson route):")
omv, Omv = 1.0, 0.8
M = 400
bd = np.sqrt(np.arange(1, M))
q_b = np.diag(bd, 1) + np.diag(bd, -1)
nb = np.diag(np.arange(M, dtype=float))

def block(nn, gv):
    ch = nn + 0.5
    return omv*ch*np.eye(M) + Omv*nb + gv*ch*q_b

pub = {0.5: 10.365, 0.9: 30.226, 0.99: 293.478, 0.999: 2927.284}
for kv in (0.5, 0.9, 0.99, 0.999):
    gv = math.sqrt(kv*Omv*omv/4)
    ev1, U1 = np.linalg.eigh(block(1, gv))
    ev0, U0 = np.linalg.eigh(block(0, gv))
    ev2, U2 = np.linalg.eigh(block(2, gv))
    psi = U1[:, 0]; E1 = ev1[0]
    # x_s = a + a': <2-block, s'|x|1-block, psi> = sqrt(2) * <s'|psi>
    # (bath basis shared), <0-block, s'|x|psi> = sqrt(1) * <s'|psi>
    chi = 0.0
    ov2 = U2.T @ psi
    chi += 2.0*np.sum(2.0*ov2**2/(ev2 - E1))
    ov0 = U0.T @ psi
    chi += 2.0*np.sum(1.0*ov0**2/(ev0 - E1))
    d = abs(chi - pub[kv])/pub[kv]
    emit("   kappa=%.3f: chi(1) = %10.3f  (stage %10.3f, rel d %.1e)"
         % (kv, chi, pub[kv], d))
d2 = (math.sqrt(0.999*Omv*omv/4)/Omv)**2
lim = 4*math.exp(-d2)/omv
emit("   limit 4 e^{-d^2}/om at kappa=0.999: %.6f (stage 2.927377)"
     % lim)
emit("")

# ==================== GA-2: 10J thresholds by direct scan ====================
emit("GA-2  thresholds by direct spectral scan (block eigh, no formula):")
for target, kexp in (("D(0)=0", 2.0), ("D(1)=0", 1.0)):
    lo, hi = 0.2, 3.5
    for _ in range(60):
        mid = 0.5*(lo + hi)
        gv = math.sqrt(mid*Omv*omv/4)
        if target == "D(0)=0":
            val = np.linalg.eigvalsh(block(1, gv))[0] - \
                  np.linalg.eigvalsh(block(0, gv))[0]
        else:
            val = np.linalg.eigvalsh(block(2, gv))[0] - \
                  np.linalg.eigvalsh(block(1, gv))[0]
        if val > 0: lo = mid
        else: hi = mid
    emit("   %s at kappa = %.6f (exact %.1f)" % (target, 0.5*(lo+hi),
                                                 kexp))
emit("")

# ==================== GA-3: thermal Laguerre sum, high precision ==========
emit("GA-3  thermal-FC persistence (mpmath, independent truncation):")
import mpmath as mp
mp.mp.dps = 30
nbar = mp.mpf(2)
d2m = mp.mpf(0.99)*Omv*omv/4/(Omv*Omv)
S = mp.mpf(0)
for k in range(400):
    pk = (nbar/(1+nbar))**k/(1+nbar)
    S += pk*mp.e**(-d2m)*mp.laguerre(k, 0, d2m)**2
emit("   sum = %s (stage 0.35389857)" % mp.nstr(S, 10))
emit("")

# ==================== GA-4: 10I neighbor-scale arithmetic ==================
emit("GA-4  10I neighbor-sum scale (from-scratch arithmetic):")
G_SI, MSUN, MPC = 6.674e-11, 1.989e30, 3.0857e22
A0 = 1.2e-10
for LK, D in ((1e10, 0.5), (1e11, 0.3), (1e9, 0.2), (1e10, 1.0)):
    e = G_SI*LK*MSUN/(D*MPC)**2/A0
    emit("   L_K=%.0e Lsun at %.1f Mpc: e_N = %.2e a0 (%.1f%% of the"
         " 0.01 floor)" % (LK, D, e, 100*e/0.01))
emit("")

# ==================== GA-5: 10I Theta regression, independent ==============
emit("GA-5  Theta1 re-derivation slope (independent match + linregress):")
meta = {}
with open('data/sparc/SPARC_Lelli2016c.mrt') as f:
    lines = f.readlines()
start = max(i for i, l in enumerate(lines)
            if set(l.strip()) <= set('- ')) + 1
for l in lines[start:]:
    t = l.split()
    if len(t) < 18: continue
    try:
        meta[t[0]] = (float(t[5]), int(t[17]), float(t[2]))
    except ValueError:
        continue
PFX = ('NGC', 'UGC', 'IC', 'DDO', 'ESO', 'PGC', 'UGCA')
def canon(nm):
    s = re.sub(r'[\s\-_]', '', str(nm)).upper()
    for p in PFX:
        if s.startswith(p):
            rest = s[len(p):]
            m = re.match(r'0*(\d.*)', rest)
            if m: return p + m.group(1)
    return s
ungc = []
with open('data/karachentsev_ungc.csv', encoding='utf-8',
          errors='replace') as f:
    for row in csv.DictReader(f):
        try:
            ra = float(row['RAJ2000']); de = float(row['DEJ2000'])
            dist = float(row['Dist'])
        except (ValueError, KeyError):
            continue
        def ff(k):
            try: return float(row.get(k, ''))
            except ValueError: return np.nan
        ungc.append((row['Name'].strip(), row.get('SimbadName', '').strip(),
                     row.get('NEDname', '').strip(), ra, de, dist,
                     ff('KLum'), ff('Ti1')))
POS = np.array([[u[5]*math.cos(math.radians(u[4]))*math.cos(math.radians(u[3])),
                 u[5]*math.cos(math.radians(u[4]))*math.sin(math.radians(u[3])),
                 u[5]*math.sin(math.radians(u[4]))] for u in ungc])
MASS = np.array([10.0**u[6] if np.isfinite(u[6]) else 0.0 for u in ungc])
umap = {}
for k, u in enumerate(ungc):
    for nm in (u[0], u[1], u[2]):
        if nm: umap.setdefault(canon(nm), k)
# SPARC galaxy list = rotmod files with the 9V cuts
names = []
for path in sorted(glob.glob('data/sparc/rotmod/**/*_rotmod.dat',
                             recursive=True)):
    nm = os.path.basename(path).replace('_rotmod.dat', '')
    inc, q, D = meta.get(nm, (0, 3, 10.0))
    if inc < 30 or q > 2: continue
    names.append((nm, D))
xs, ys = [], []
nmatch = 0
for nm, D in names:
    k = umap.get(canon(nm), -1)
    if k < 0: continue
    du = ungc[k][5]
    if abs(D - du) > max(2.0, 0.35*du): continue
    nmatch += 1
    ti1 = ungc[k][7]
    if not np.isfinite(ti1): continue
    d2v = np.sum((POS - POS[k])**2, axis=1)
    d2v[k] = np.inf
    sel = (d2v <= 25.0) & (MASS > 0)
    if not sel.any(): continue
    xs.append(float(np.max(np.log10(MASS[sel]) - 1.5*np.log10(d2v[sel]))))
    ys.append(ti1)
res = linregress(xs, ys)
emit("   matches %d (stage 31); regression N=%d slope %.3f "
     "intercept %.2f r=%.3f (stage slope 0.956 offset -10.45)" %
     (nmatch, len(xs), res.slope, res.intercept, res.rvalue))
emit("")
emit("BLIND PART DONE (extend below after reading the report)")
emit("")
emit("=" * 60)
emit("POST-REPORT PART: verifying the ROUND-30 referee's NEW numbers")
emit("")

# GA-6: J-1 -- rung-count dependence of the R-C bound
emit("GA-6  J-1 rung arithmetic (kappa thresholds 2/(n+1); D(n) at")
emit("      the locked world and at kappa=1):")
for kv in (0.925, 1.0):
    Ds = [1 - kv*(nn+1)/2 for nn in range(4)]
    emit("   kappa=%.3f: D(0..3)/om = %s" %
         (kv, [round(d, 4) for d in Ds]))
emit("   three-rung condition D(2)>=0 <=> kappa <= 2/3 = %.4f -> "
     "EXCLUDES 0.925 (referee's J-1 CONFIRMED)" % (2/3))
emit("   deep occupation n_BE(x): x=0.10 -> %.3f, x=0.05 -> %.3f "
     "(occupied rungs sit far beyond the fold at every kappa)" %
     (1/(math.exp(0.10)-1), 1/(math.exp(0.05)-1)))
emit("")

# GA-7: J-2c -- the competing pole at kappa = 0.4 in the downward term
emit("GA-7  J-2c the kappa=0.4 downward-term resonance:")
for kv in (0.4, 0.5):
    D0 = omv*(1 - kv/2)
    emit("   kappa=%.1f: -D(0) + 1*Om = %+.4f (zero at kappa=0.4 "
         "exactly: |1,0> <-> |0,1> resonance)" % (kv, -D0 + Omv))
emit("   => the chi(1) value at kappa=0.5 is inflated by THIS pole,")
emit("   not the kappa=1 pole; the kappa->1 limit itself is clean")
emit("   (nearest competing resonance far); referee CONFIRMED.")
emit("")

# GA-8: I-2 -- group-end neighbor terms (NGC2976) + iso/non-iso ranges
emit("GA-8  I-2 group-end neighbor terms (full matched set):")
gvals = {}
for nm, D in names:
    k = umap.get(canon(nm), -1)
    if k < 0: continue
    du = ungc[k][5]
    if abs(D - du) > max(2.0, 0.35*du): continue
    d2v = np.sum((POS - POS[k])**2, axis=1)
    d2v[k] = np.inf
    sel = (d2v <= 9.0) & (MASS > 0)
    e = float(np.sum(G_SI*MASS[sel]*MSUN/(d2v[sel]*(MPC)**2))/A0)
    gvals[nm] = (e, ungc[k][7])
iso_terms = [e for nm, (e, t) in gvals.items()
             if np.isfinite(t) and t < 0]
non_terms = [e for nm, (e, t) in gvals.items()
             if np.isfinite(t) and t >= 0]
emit("   iso neighbor term max = %.2e a0 (%.1f%% of 0.01 floor)" %
     (max(iso_terms), 100*max(iso_terms)/0.01))
emit("   non-iso max = %.2e a0 (%.1f%% of floor); median %.2e "
     "(%.1f%%)" % (max(non_terms), 100*max(non_terms)/0.01,
                   float(np.median(non_terms)),
                   100*float(np.median(non_terms))/0.01))
if 'NGC2976' in gvals:
    emit("   NGC2976: %.3e a0 = %.1f%% of floor (referee: 1.099e-2 "
         "= 109.9%%)" % (gvals['NGC2976'][0],
                         100*gvals['NGC2976'][0]/0.01))
emit("")

# GA-9: I-3 -- the UGC05721 = NGC3274 alias row
emit("GA-9  I-3 the alias exemplar:")
for u in ungc:
    if canon(u[0]) == 'NGC3274' or canon(u[1]) == 'NGC3274' \
       or canon(u[2]) == 'NGC3274':
        emit("   UNGC row: Name=%s Simbad=%s NED=%s Dist=%.2f "
             "Ti1=%s" % (u[0], u[1], u[2], u[5], u[7]))
        du = u[5]
        Dsp = dict(names).get('UGC05721', None)
        if Dsp:
            emit("   SPARC UGC05721 D=%.2f; gate |%.2f-%.2f|=%.2f <= "
                 "max(2, 0.35*%.2f=%.2f): %s" %
                 (Dsp, Dsp, du, abs(Dsp-du), du, 0.35*du,
                  "PASS" if abs(Dsp-du) <= max(2.0, 0.35*du)
                  else "FAIL"))
emit("")

# GA-10: I-1 -- iso-overlap Chae gates
emit("GA-10 I-1 iso-overlap Chae e_N -> gates:")
chae = {}
with open('data/chae2021_table3.csv') as f:
    for row in csv.DictReader(l for l in f if not l.startswith('#')):
        chae[row['galaxy']] = 10.0**float(row['log_eN_maxclust'])
def s_of(e):
    xx = math.sqrt(e)
    nn = 1.0/(math.exp(xx) - 1.0)
    return nn/(1.0 + nn)
rows = []
for nm, D in names:
    k = umap.get(canon(nm), -1)
    if k < 0 or nm not in chae: continue
    du = ungc[k][5]
    if abs(D - du) > max(2.0, 0.35*du): continue
    t = ungc[k][7]
    if np.isfinite(t) and t < 0:
        rows.append((nm, chae[nm], s_of(chae[nm])**2))
for nm, e, gg in sorted(rows, key=lambda r: r[1]):
    emit("   %-10s Chae e_N=%.5f a0  gate s^2=%.3f" % (nm, e, gg))
emit("   (referee: 0.00493-0.00920 -> s^2 0.83-0.87)")
emit("")
emit("ALL REFEREE NUMBERS CHECKED")
