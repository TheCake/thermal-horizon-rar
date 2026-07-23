"""
STAGE 5Y (O15a): kappa tables + the a0 ladder for the running-beta
functions -- the temperature-strain adjudication.

The sharp fixed functions bought galaxy likelihood at the price of the
binary a0 translation (+4.3/+5.1 sigma, 5L). The running functions'
alpha-hats (1.395/1.487, and the 6B two-leg values when they land) face
the same ledger: a0_bin = 1.2e-10 * alpha^(1/kappa), kappa from their
own e_N = 1.0/1.2/1.4 EFE tables averaged over the 4V deep-pair y
distribution (5L machinery verbatim).

Gates: G1 alpha=1 -> a0 = 1.2e-10 identically; G2 kappa > 0; G3 the BE
kappa recomputed here matches 4V's +0.916 to 0.01. The e_N = 1.2 tables
are 5W/6B's (waited on, not rebuilt -- no write race); this stage builds
only the 1.0/1.4 kappa variants. Solver results cached; the alpha-join
re-runs cheaply after 6B. Writes data/stage5y_kladder.txt.
"""
import math, os, re, time
import numpy as np
from astropy.io import fits

A0_TAB = 1.2e-10
C = 2.998e8

src_solver = open('calcs/qumond_efe_solver.py', encoding='utf-8-sig').read()
ns = {}
exec(src_solver.split('# --- Gate 1')[0], ns)
solve, r = ns['solve'], ns['r']
y = 1.0/r**2

def nu_simple(yv):
    return 0.5 + np.sqrt(0.25 + 1.0/np.clip(yv, 1e-14, None))

def make_run(kind):
    def nu_run(yv):
        yv = np.clip(np.asarray(yv, float), 1e-14, None)
        ly = np.log(yv)
        nu = nu_simple(yv)
        for _ in range(80):
            if kind == 'rb1':
                b = 1.0/(2.0*nu); db = -1.0/(2.0*nu*nu)
            elif kind == 'rb2':
                b = 1.0/(2.0*(2.0*nu-1.0)); db = -1.0/((2.0*nu-1.0)**2)
            elif kind == 'rb3':
                b = 0.5/(nu*nu); db = -1.0/(nu**3)
            else:
                b = 0.5/((2.0*nu-1.0)**2); db = -2.0/((2.0*nu-1.0)**3)
            u = np.exp(np.minimum(0.5*(1.0+b)*ly + b*np.log(nu), 60.0))
            eu = np.exp(np.minimum(u, 60.0))
            em1 = np.maximum(eu - 1.0, 1e-300)
            n = np.where(u < 60.0, 1.0/em1, 0.0)
            F = nu - 1.0 - n
            dudnu = u*(db*(0.5*ly + np.log(nu)) + b/nu)
            dF = 1.0 + (eu/(em1*em1))*dudnu
            step = F/np.where(np.abs(dF) > 1e-12, dF, 1e-12)
            nu = np.maximum(nu - np.clip(step, -0.5*nu, 0.5*nu),
                            1.0 + 1e-15)
        return nu
    return nu_run

FAMS = ('rb1', 'rb2', 'rb3', 'rb4')
for fam in FAMS:
    for eN, tag in ((1.0, 'g1p0'), (1.4, 'g1p4')):
        path = f'data/efe_boost_{fam}_{tag}.npy'
        if not os.path.exists(path):
            b = solve(make_run(fam), eN)
            np.save(path, np.stack([y, b]))
            print(f"built {path}", flush=True)

# wait (not rebuild) for the 1.2 tables owned by 5W/6B
for fam in FAMS:
    path = f'data/efe_boost_{fam}_g1p2.npy'
    t0 = time.time()
    while not os.path.exists(path):
        if time.time() - t0 > 600:
            raise RuntimeError(f"timeout waiting for {path} (6B builds it)")
        time.sleep(10)

def load_tab(path):
    t = np.load(path)
    yy, b = t[0][::-1], t[1][::-1]
    return np.log(yy), np.log(np.maximum(b-1.0, 1e-12))

TAGS = {1.0: 'g1p0', 1.2: 'g1p2', 1.4: 'g1p4'}
tabs = {}
for fam in ('be',) + FAMS:
    for e, tag in TAGS.items():
        tabs[(fam, e)] = load_tab(f'data/efe_boost_{fam}_{tag}.npy')

def lnB1(fam, e, lny):
    lo, hi = (1.0, 1.2) if e <= 1.2 else (1.2, 1.4)
    w = (e-lo)/(hi-lo)
    la = np.interp(lny, *tabs[(fam, lo)])
    lb = np.interp(lny, *tabs[(fam, hi)])
    return (1-w)*la + w*lb

def kappa(fam, lny, e=1.184, dl=0.10, de=0.08):
    dy = (lnB1(fam, e, lny+dl) - lnB1(fam, e, lny-dl))/(2*dl)
    dE = (lnB1(fam, e+de, lny) - lnB1(fam, e-de, lny))/(2*de*1.0/e)
    return -(dy + dE)

# deep-sample y distribution (4V/5L verbatim)
d = fits.open('data/edr3_binaries.fits.gz', memmap=False)[1].data
plx1, plx2 = d['parallax1'], d['parallax2']
eplx1, eplx2 = d['parallax_error1'], d['parallax_error2']
sep, Rch = d['sep_AU'], d['R_chance_align']
G1m, G2m = d['phot_g_mean_mag1'], d['phot_g_mean_mag2']
plx = 0.5*(plx1+plx2)
MG1 = G1m+5*np.log10(np.maximum(plx1, 1e-6))-10
MG2 = G2m+5*np.log10(np.maximum(plx2, 1e-6))-10
sigv = 4.74047/plx*np.sqrt(d['pmra_error1']**2+d['pmdec_error1']**2
                           + d['pmra_error2']**2+d['pmdec_error2']**2)
ok = (Rch < 0.01) & (plx1 > 5) & (plx2 > 5) \
   & (plx1/np.maximum(eplx1, 1e-6) > 20) & (plx2/np.maximum(eplx2, 1e-6) > 20) \
   & (np.abs(plx1-plx2) < 3*np.hypot(eplx1, eplx2)) \
   & (sep > 200) & (sep < 50000) \
   & (MG1 > 2.6) & (MG1 < 14.2) & (MG2 > 2.6) & (MG2 < 14.2) & (sigv < 0.03)
MG_T = np.array([2.6,3.4,4.2,4.8,5.4,6.0,6.83,7.57,8.16,8.82,9.29,10.05,
                 11.21,12.45,14.26])
MS_T = np.array([1.60,1.33,1.12,1.00,0.90,0.82,0.70,0.64,0.57,0.50,0.44,
                 0.37,0.23,0.162,0.102])
Mtot = np.interp(MG1[ok], MG_T, MS_T)+np.interp(MG2[ok], MG_T, MS_T)
s_kau = sep[ok]/1e3
deep = (s_kau >= 6) & (s_kau < 30)
GKMS = 0.8868
UNIT = 6.685e-9
yd = (GKMS*Mtot[deep]/s_kau[deep]**2)*UNIT/A0_TAB
lnyd = np.log(np.clip(yd, 1e-3, 50.0))

# alpha-hats
def ahat(path, laws):
    if not os.path.exists(path): return {}
    txt = open(path).read()
    out = {}
    for law in laws:
        vals = [float(m) for m in re.findall(
            rf'seed \d+ {law}: a_hat=([0-9.]+)', txt)]
        if vals:
            out[law] = (np.mean(vals),
                        np.std(vals, ddof=1)/math.sqrt(len(vals)), len(vals))
    return out

ah = {}
ah.update(ahat('data/stage5w_summary.txt', ('rb1', 'rb2')))
ah.update(ahat('data/stage6b_summary.txt', ('rb3', 'rb4')))

L = [f"STAGE 5Y kappa/a0 ladder for the running functions; deep pairs "
     f"(6-30 kAU): {int(deep.sum())}", ""]

kap = {}
for fam in ('be',) + FAMS:
    kv = kappa(fam, lnyd)
    kap[fam] = float(np.mean(kv))
    L.append(f"kappa({fam:>4}) = {np.mean(kv):+.3f} "
             f"(16/84 {np.percentile(kv,16):+.3f}/{np.percentile(kv,84):+.3f})")
ok3 = abs(kap['be'] - 0.916) < 0.01
L.append(f"G3 BE-kappa regression vs 4V (+0.916): {kap['be']:+.3f} -> "
         f"{'PASS' if ok3 else 'FAIL'}")
assert ok3
g2 = all(k > 0 for k in kap.values())
L.append(f"G2 kappa > 0 -> {'PASS' if g2 else 'FAIL'}")
assert g2
g1v = A0_TAB*1.0**(1.0/kap['rb1'])
L.append(f"G1 alpha=1 -> {g1v:.3e} -> "
         f"{'PASS' if abs(g1v-A0_TAB) < 1e-14 else 'FAIL'}")
L.append("")

def cho2pi(H0, sH0):
    h = H0*1e3/3.0857e22
    v = C*h/(2*math.pi)
    return v, v*sH0/H0
tp, stp = cho2pi(67.4, 0.5)

L.append("binary a0 = 1.2e-10 * alpha^(1/kappa)  [x1e-10; alpha error = "
         "realization SE (+) 0.11 systematic in quadrature, 5L convention]")
L.append(f"  cH0/2pi (Planck) = {tp*1e10:.3f} +- {stp*1e10:.3f}")
for fam in FAMS:
    if fam not in ah:
        L.append(f"  {fam:>4}: [6B pending]")
        continue
    m, s, n = ah[fam]
    stot = math.hypot(s, 0.11)
    k = kap[fam]
    a0b = A0_TAB*m**(1.0/k)
    sa0b = a0b*stot/(m*k)
    pull = (a0b*1e10-tp*1e10)/math.hypot(sa0b*1e10, stp*1e10)
    L.append(f"  {fam:>4}: alpha {m:.3f}+-{stot:.3f} ({n} seeds), kappa "
             f"{k:+.3f} -> a0 = {a0b*1e10:.2f}+-{sa0b*1e10:.2f}  pull "
             f"{pull:+.1f} sigma")
L.append("")
L.append("  [comparators (5L/4V): BE 1.37+-0.17 (+1.9s); p065 1.59+-0.11 "
         "(+5.1s); gm 1.51+-0.11 (+4.3s); galaxy vertical-hier a0 for "
         "rb1/rb2 = 0.99/0.98 (5V, controlled treatment)]")

out = "\n".join(L)
print(out)
with open('data/stage5y_kladder.txt', 'w') as f:
    f.write(out + "\n")
print("\nsaved: data/stage5y_kladder.txt")
