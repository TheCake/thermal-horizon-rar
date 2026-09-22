"""ROUND 52 GB addendum: two corrections to the GB half itself,
disclosed (trap #23 -- the verifier is an instrument too):
  (a) GB-3's z(Y) used the 55-galaxy WORLD-Y prediction against the
      16-galaxy sky value; the reviewer's subsample-consistent
      construction is correct -- recompute WORLD-Y's d_Oi ON the
      both-sides subsample and the proper z.
  (b) GB-7's worst-single-drop print had an initialization bug
      (worst=0.0 dominates the comparison); recompute honestly.
Appends to data/round52_gb.txt.
"""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = open(os.path.join(HERE, 'stage10y_organizer.py'),
           encoding='utf-8').read()
cut = SRC.index('# ====')
NS = {'__name__': 'stage10y_trunc',
      '__file__': os.path.join(HERE, 'stage10y_organizer.py')}
_argv = sys.argv
sys.argv = ['stage10y_organizer.py', 'gates']
exec(compile(SRC[:cut], 'stage10y_trunc', 'exec'), NS)
sys.argv = _argv

GAL_A, QUAL_A = NS['GAL_A'], NS['QUAL_A']
TH, RDF = NS['TH_STAR'], NS['RD_FAC']

L = []
def P(s=""):
    print(s, flush=True)
    L.append(s)

P("")
P("-- GB addendum (verifier corrections, disclosed) --")

# (a) WORLD-Y d_Oi restricted to the both-sides subsample
sigY = {j: NS['B_Y']*np.log10(GAL_A[j]['y']) for j in range(len(GAL_A))}
doiY, doiS = [], []
for j in QUAL_A:
    gd_ = GAL_A[j]
    if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
        continue
    it = gd_['th'] < TH
    ir = gd_['rk'] < RDF*gd_['rd']
    Oo = ~it & ~ir
    Oi = ~it & ir
    if Oo.any() and Oi.any():
        doiY.append(float(sigY[j][Oi].mean() - sigY[j][Oo].mean()))
        doiS.append(float(gd_['res'][Oi].mean() -
                          gd_['res'][Oo].mean()))
doiY = float(np.mean(doiY))
m16 = float(np.mean(doiS))
s16 = float(np.std(doiS, ddof=1)/math.sqrt(len(doiS)))
P(f"  (a) WORLD-Y d_Oi on the both-sides-{len(doiS)} subsample = "
  f"{doiY:+.5f}; z(Y) = {(m16 - doiY)/s16:+.2f} "
  f"[reviewer +1.83; the GB-3 +1.43 used the 55-galaxy world value "
  f"and is RETIRED]")

# (b) worst single drop of d_Ii, honest recompute
dii = []
for gd_ in GAL_A:
    if not np.isfinite(gd_['rd']) or gd_['rd'] <= 0:
        continue
    it = gd_['th'] < TH
    ir = gd_['rk'] < RDF*gd_['rd']
    Oo = ~it & ~ir
    if not Oo.any() or not (it & ir).any():
        continue
    dii.append(float(gd_['res'][it & ir].mean()) -
               float(gd_['res'][Oo].mean()))
dii = np.array(dii)
m0 = dii.mean()
worst, wk = m0, -1
for k in range(len(dii)):
    mk = float(np.delete(dii, k).mean())
    if abs(mk - m0) > abs(worst - m0):
        worst, wk = mk, k
P(f"  (b) d_Ii worst single drop -> {worst:+.5f} "
  f"({(m0-worst)/m0:.1%} of the mean) [reviewer -0.07024 / 21.3%; "
  f"the GB-7 +0.00000 was an init bug in the verifier, disclosed]")

with open('data/round52_gb.txt', 'a', encoding='utf-8') as f:
    f.write("\n".join(L) + "\n")
print("\nappended: data/round52_gb.txt")
