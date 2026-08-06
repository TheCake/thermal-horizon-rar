"""THE RNAAS NOTE FIGURE (provenance-gated). One panel: the fpm meter
per RUWE quartile (9L fiducial, both seeds) over the 9P mass-model
deformation envelope. Gates regress every plotted number against the
committed stage outputs before drawing.
Output: papers/figs/note_fig1_meter.png/.pdf + provenance lines in
data/note_figure.txt
"""
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

P("NOTE FIGURE - provenance-gated build")

t9l = open('data/stage9l_fpmmeter.txt', encoding='utf-8',
           errors='replace').read()
rows = re.findall(r"\[seed (\d+)\] Q(\d) narrow: fpm marginal = "
                  r"[\d./]+; E\[fpm\] = ([\d.]+)", t9l)
fid = {}
for sd, q, e in rows:
    fid[(int(sd), int(q))] = float(e)
assert len(fid) == 8, fid

gz = np.load('data/stage9p_grid.npz', allow_pickle=True)
E = gz['E']            # (variants, quartile, seed); E[0] = fiducial
t9p = open('data/stage9p_massmodel.txt', encoding='utf-8',
           errors='replace').read()
m = re.search(r"envelope over V1-V5 x Q1-Q4 x seeds = "
              r"\[([\d.]+), ([\d.]+)\]", t9p)
env_lo, env_hi = float(m.group(1)), float(m.group(2))

g1 = all(abs(fid[(sd, q+1)] - E[0, q, si]) <= 0.005
         for si, sd in enumerate((31, 101)) for q in range(4))
P(f"G-N1 9L printed rows vs 9P fiducial table (8 values, 0.005): "
  f"{'PASS' if g1 else 'FAIL'}")
sub = E[1:, :, :]
g2 = (abs(float(sub.min()) - env_lo) <= 0.005
      and abs(float(sub.max()) - env_hi) <= 0.005)
P(f"G-N2 envelope regression ({env_lo}/{env_hi} vs "
  f"{float(sub.min()):.2f}/{float(sub.max()):.2f}): "
  f"{'PASS' if g2 else 'FAIL'}")
if not (g1 and g2):
    P("GATES FAILED - STOP, no figure written")
    with open('data/note_figure.txt', 'w') as f:
        f.write("\n".join(L_)+"\n")
    raise SystemExit(1)

plt.rcParams.update({'font.size': 9, 'figure.dpi': 150})
fig, ax = plt.subplots(figsize=(4.4, 3.2))
qx = np.array([1, 2, 3, 4])
lo_q = sub.min(axis=(0, 2)); hi_q = sub.max(axis=(0, 2))
ax.fill_between(qx, lo_q, hi_q, color='0.88',
                label='mass-model deformation envelope')
for si, (sd, mk) in enumerate(((31, 'o'), (101, 's'))):
    ax.plot(qx, [fid[(sd, q)] for q in qx], mk + '-', ms=4.5, lw=1,
            color='k' if si == 0 else '0.45',
            label=f'fiducial, realization {sd}')
ax.axhline(1.0, color='0.5', ls=':', lw=0.8)
ax.text(3.62, 1.03, 'formal errors', fontsize=7, color='0.35')
ax.axhline(1.4, color='tab:blue', ls='--', lw=0.9)
ax.text(3.35, 1.43, 'single-star ceiling', fontsize=7,
        color='tab:blue')
ax.axhline(env_lo, color='tab:red', ls='--', lw=0.9)
ax.text(1.02, env_lo-0.09, f'deformation floor ({env_lo:.2f})',
        fontsize=7, color='tab:red')
ax.set_xticks(qx)
ax.set_xticklabels(['Q1\n(cleanest)', 'Q2', 'Q3', 'Q4\n(worst)'])
ax.set_xlabel('RUWE quartile')
ax.set_ylabel(r'$E[f_{\rm pm}]$ (0.2-2 kAU, boost-free)')
ax.set_ylim(0.9, 2.75)
ax.legend(fontsize=6.5, loc='upper left')
fig.tight_layout()
for ext in ('png', 'pdf'):
    fig.savefig(f'papers/figs/note_fig1_meter.{ext}')
plt.close(fig)
P("note_fig1 written")
with open('data/note_figure.txt', 'w') as f:
    f.write("\n".join(L_)+"\n")
P("done")
