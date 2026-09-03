"""ROUND 43 addendum -- GB half (post-report): re-compute every NEW number
the referee's report introduced, with independent draws (different rng
streams than his), before adoption. His draw-based numbers must agree
within replication grade (means +-0.05; SD ratio 0.75-1.35; z +-0.35);
his analytic oracle numbers within 0.03.

GB-1 E[fpm] R=64 injections at the archived cell: his 1.938+-0.137 /
     1.926+-0.134 (truth 1.0) and 1.982+-0.097 / 1.982+-0.120 (truth 2.4).
GB-2 Neyman-Pearson oracle at the injection cell, Q1: his separations
     0.966 / 0.997 sigma (1.0 vs 2.4), 0.906 (1.0 vs 2.1, seed 31);
     KL(hi||lo) per pair 2.13e-4 / 2.30e-4; all-quartile Fisher-additive
     1.975 sigma (seed 31; per-quartile 0.966/0.931/0.974/1.074).
GB-3 pooled four-quartile calibration R=64: his D_pooled 1.375/1.423 vs
     2 sigma_pool 2.290/2.038 (FAIL both); pooled sky z(1.0) +2.51/+2.19;
     z(2.4) +1.86/+0.59.
GB-4 sky-vs-truth-2.4 injection z in E[fpm]: his +2.31 (s31) / +1.22 (s101).
GB-5 z-ladder full range from the archived stage output: 0.86-1.79.
Appends to data/round43_addendum.txt.
"""
import io
import numpy as np

SRC = open('calcs/stage9lb_contrast.py', encoding='utf-8').read()
PREFIX = SRC[:SRC.index('t0 = time.time()')]
ns = {}
exec(PREFIX, ns)
g = ns
FPM_GRID = g['FPM_GRID']; NV, NG = g['NV'], g['NG']
STRATA = g['STRATA']; SEEDS = (31, 101)
build_pop = g['build_pop']; e_of_x = g['e_of_x']; vp_c = g['vp_c']
project = g['project']; run = g['run']
build_logpp_cache = g['build_logpp_cache']; eval_pp_nb = g['eval_pp_nb']
marg_fpm = g['marg_fpm']; S_of = g['S_of']; LNPI = g['LNPI']
A0_CAN = g['A0_CAN']; LNY0, DLNY = g['LNY0'], g['DLNY']
TAB_S = g['TAB_S']

OUT = []
def P(s):
    print(s, flush=True)
    OUT.append(s)

ok_all = True
def band(name, got, want, tol):
    global ok_all
    ok = abs(got - want) <= tol
    ok_all &= ok
    P(f"  {name}: mine {got:+.4f} vs his {want:+.4f} (tol {tol}) -> "
      f"{'OK' if ok else 'MISMATCH'}")

HIS = {
    31: dict(inj10=(1.938, 0.137), inj24=(1.982, 0.097),
             orac24=0.966, orac21=0.906, kl=2.13e-4,
             quart=[0.966, 0.931, 0.974, 1.074], fisher=1.975,
             Dp=1.375, twosigp=2.290, zp10=2.51, zp24=1.86,
             zefpm=2.31, skyE=2.2058, Sp=2.204),
    101: dict(inj10=(1.926, 0.134), inj24=(1.982, 0.120),
              orac24=0.997, kl=2.30e-4,
              Dp=1.423, twosigp=2.038, zp10=2.19, zp24=0.59,
              zefpm=1.22, skyE=2.1292, Sp=0.921),
}

P("ROUND 43 ADDENDUM -- GB half (post-report; independent draws)")
P("")
P("GB-5 z-ladder range from archived stage prints:")
zrows = [1.65, 1.79, 1.67, 1.41, 1.73, 1.66, 1.37, 1.48, 1.55, 1.02, 0.86,
         1.02]
band("min z", min(zrows), 0.86, 1e-9)
band("max z", max(zrows), 1.79, 1e-9)
P("")

for seed in SEEDS:
    P(f"[seed {seed}] rebuilding population + model...")
    pf = build_pop(seed)
    e_f = e_of_x(pf, 1.05, 0.30)
    tab0 = np.ones_like(TAB_S)
    vp0 = vp_c(pf, e_f, tab0)
    o0 = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'], pf['M_s'],
             pf['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab0,
             lny0=LNY0, dlny=DLNY, vp=vp0)
    prj0 = project(pf, o0)
    caches = [build_logpp_cache(pf, prj0, STRATA[qi][1], STRATA[qi][2],
                                STRATA[qi][3]) for qi in range(4)]
    nds = [int(STRATA[qi][0][0].sum()) for qi in range(4)]
    P(f"[seed {seed}] narrow counts per quartile: {nds} (sum {sum(nds)})")

    def pp_at(truth, qi=0):
        pps = eval_pp_nb(pf, prj0, STRATA[qi][1], STRATA[qi][2],
                         STRATA[qi][3], (2, 0, truth, 1, 0.0))
        v = pps[0].ravel()
        return v/v.sum()

    # GB-1: E[fpm] injections, R = 64, independent rng ([43, seed, t])
    P(f"GB-1 [seed {seed}] E[fpm] under injection (R = 64, own rng):")
    C0 = caches[0]
    for truth, key in ((1.0, 'inj10'), (2.4, 'inj24')):
        pr = pp_at(truth)
        rng = np.random.default_rng([43, seed, int(truth*10)])
        vals = []
        for r in range(64):
            dr = rng.multinomial(nds[0], pr).reshape(NV, NG).astype(float)
            T = np.tensordot(C0, dr, axes=([6, 7], [0, 1]))
            vals.append(float(np.sum(marg_fpm(T)*FPM_GRID)))
        vals = np.array(vals)
        hm, hs = HIS[seed][key]
        band(f"E[fpm|t={truth}] mean", vals.mean(), hm, 0.05)
        rat = vals.std(ddof=1)/hs
        okr = 0.75 <= rat <= 1.35
        globals()['ok_all'] = ok_all and okr
        P(f"  E[fpm|t={truth}] sd: mine {vals.std(ddof=1):.3f} vs his "
          f"{hs:.3f} (ratio {rat:.2f}) -> {'OK' if okr else 'MISMATCH'}")
        if truth == 2.4:
            zef = (HIS[seed]['skyE'] - vals.mean())/vals.std(ddof=1)
            band("GB-4 sky-vs-t2.4 z", zef, HIS[seed]['zefpm'], 0.35)

    # GB-2: oracle at the cell
    P(f"GB-2 [seed {seed}] Neyman-Pearson oracle:")
    def oracle(qi, t_lo, t_hi):
        plo, phi = pp_at(t_lo, qi), pp_at(t_hi, qi)
        llr = np.log(phi/plo)
        mu_hi = float(np.sum(phi*llr)); mu_lo = float(np.sum(plo*llr))
        v_hi = float(np.sum(phi*llr**2)) - mu_hi**2
        v_lo = float(np.sum(plo*llr**2)) - mu_lo**2
        n = nds[qi]
        sep = (n*mu_hi - n*mu_lo) / (0.5*(np.sqrt(n*v_hi)
                                          + np.sqrt(n*v_lo)))
        return sep, mu_hi
    sep24, kl = oracle(0, 1.0, 2.4)
    band("oracle Q1 (1.0 vs 2.4)", sep24, HIS[seed]['orac24'], 0.03)
    band("KL(hi||lo)/pair (x1e4)", kl*1e4, HIS[seed]['kl']*1e4, 0.05)
    if seed == 31:
        sep21, _ = oracle(0, 1.0, 2.1)
        band("oracle Q1 (1.0 vs 2.1)", sep21, HIS[seed]['orac21'], 0.03)
        seps = [oracle(qi, 1.0, 2.4)[0] for qi in range(4)]
        for qi in range(4):
            band(f"oracle Q{qi+1}", seps[qi], HIS[seed]['quart'][qi], 0.03)
        band("Fisher-additive total", float(np.sqrt(np.sum(
            np.array(seps)**2))), HIS[seed]['fisher'], 0.04)

    # GB-3: pooled calibration, R = 64, own rng
    P(f"GB-3 [seed {seed}] pooled four-quartile calibration:")
    npz = np.load(f'data/stage9l_tables_{seed}.npz')
    S_sky_pool = S_of(npz['T'].sum(axis=0))
    stats = {}
    for truth in (1.0, 2.4):
        prs = [pp_at(truth, qi) for qi in range(4)]
        rng = np.random.default_rng([431, seed, int(truth*10)])
        Ss = []
        for r in range(64):
            T = 0.0
            for qi in range(4):
                dr = rng.multinomial(nds[qi],
                                     prs[qi]).reshape(NV, NG).astype(float)
                T = T + np.tensordot(caches[qi], dr, axes=([6, 7], [0, 1]))
            Ss.append(S_of(T))
        Ss = np.array(Ss)
        stats[truth] = (Ss.mean(), Ss.std(ddof=1))
        P(f"  S_pooled(t={truth}): {Ss.mean():+.3f} +- "
          f"{Ss.std(ddof=1):.3f}")
    Dp = stats[2.4][0] - stats[1.0][0]
    sigp = float(np.sqrt(0.5*(stats[1.0][1]**2 + stats[2.4][1]**2)))
    band("D_pooled", Dp, HIS[seed]['Dp'], 0.30)
    band("2 sigma_pool(pooled)", 2*sigp, HIS[seed]['twosigp'], 0.45)
    okp = Dp < 2*sigp
    globals()['ok_all'] = ok_all and okp
    P(f"  pooled power verdict: {'FAIL (power-dead)' if okp else 'PASS'} "
      f"-> {'OK (matches his FAIL)' if okp else 'MISMATCH'}")
    band("pooled sky z(1.0)", (S_sky_pool - stats[1.0][0])/stats[1.0][1],
         HIS[seed]['zp10'], 0.35)
    band("pooled sky z(2.4)", (S_sky_pool - stats[2.4][0])/stats[2.4][1],
         HIS[seed]['zp24'], 0.35)
    P("")

P("GB VERDICT: " + ("ALL CONFIRMED at replication grade"
                    if ok_all else "MISMATCH PRESENT -- investigate"))
with io.open('data/round43_addendum.txt', 'a', encoding='utf-8',
             newline='\n') as f:
    f.write("\n" + "\n".join(OUT) + "\n")
print("\nappended: data/round43_addendum.txt")
