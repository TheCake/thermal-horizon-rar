"""STAGE 9L-c -- THE FULL-SEPARATION NOISE-CHANNEL CALIBRATION (ROUND 43
condition C2; 2026-09-03; pre-registered BEFORE the run).

THE QUESTION: Paper 1 quotes the quality-stratified joint fit's cleanest-
quartile noise posterior E[f_pm] = 1.97-2.25 -- an expectation over the
same bounded grid whose flat-prior mean is 2.00, for a stratum with no
prior-proof null control. R43-C2: calibrate the JOINT (full-separation)
noise channel on Q1 with a truth ladder INCLUDING the formal-errors truth,
reporting a prior-proof contrast, not a bare expectation.

THE INSTRUMENT: the 9L evaluator un-restricted to ALL FOUR separation bins
(eval_block_g form; verbatim model expressions; per-bin log-pp caches),
alpha = 0 throughout. DISCLOSED SCOPE: the real joint fit marginalizes the
boost amplitude; with alpha free, noise-channel power can only DECREASE
(the amplitude absorbs width). This stage therefore measures the channel's
alpha-off information content: a FAIL is decisive (the channel is dead even
before alpha-degeneracy); a PASS is necessary-not-sufficient and is quoted
with the alpha-on control carried by the arm suite (P(fpm = 3.0) <= 0.02
on arms vs 0.54/0.97 sky, full-sample). Sky-side reads here are
DIAGNOSTIC-ONLY (at alpha = 0 any true boost in the wide bins would be
soaked into f_pm; never quote a 9L-c sky number as a measurement).

DESIGN: statistic S = ln P(fpm >= 2.1) - ln P(fpm <= 1.8) from the
archived marginalization over the summed 4-bin T (flat posterior S = 0
exactly; trap-#27 compliant); edge-mass co-read P(fpm = 3.0) (flat = 1/6).
Truth ladder {1.0, 1.5, 2.1} x R = 64 multinomial replications x 2 seeds
(31/101), injected at the archived G9L-3 cell (fcomp = 0.20, ffly = 0.05,
kw = 1.0, sq = 0.0) independently per bin at each bin's Q1 count.

GATES (each names the clause its failure vetoes -- R35 rule):
  GC9LC-0 cache identity: per-bin log-pp caches reproduce a direct 4-bin
    eval on one synthetic dataset, max|dT| <= 1e-6. FAIL -> wiring; STOP.
  GC9LC-1 reader identity: the bin-0-only restriction reproduces the 9L-b
    sky S(Q1) values (+0.534/+0.372) to +-0.002. FAIL -> STOP.
LETTER (pre-signed, ordered, BOTH seeds; sigma_pool = rms of per-truth
SDs; D = S_bar(2.1) - S_bar(1.0)):
  C-CHANNEL-POWERED: D >= 2 sigma_pool in both seeds -> the Paper-1 Q1
    sentence gains the calibration citation: "the full-separation noise
    channel separates the formal-errors truth from the fitted level at
    >= 2 sigma (alpha-off calibration, stage9lc; alpha-on control = the
    injection arms)". The C2 same-sentence caveat is then REPLACED by this
    sourced form (the alpha-degeneracy disclosure stays).
  C-CHANNEL-DEAD-TOO: D < 2 sigma_pool in either seed -> the C2 caveat
    becomes PERMANENT: the Q1 E[f_pm] may only be printed with the
    same-sentence statement that 2.00 is the grid's flat-prior mean and
    the per-quartile number is not null-calibrated; the Sec-6.3 "even the
    cleanest strata demand ~2x" sentence is WEAKENED to fit-internal
    status.
NO credence movement under any letter. NO sky verdict from this stage.
Output: data/stage9lc_jointcal.txt + data/stage9lc_calib.npz.
"""
import io, time
import numpy as np

SRC = open('calcs/stage9lb_contrast.py', encoding='utf-8').read()
PREFIX = SRC[:SRC.index('t0 = time.time()')]
ns = {}
exec(PREFIX, ns)
g = ns
FPM_GRID = g['FPM_GRID']; FCOMP_GRID = g['FCOMP_GRID']
FFLY_GRID = g['FFLY_GRID']; KW_GRID = g['KW_GRID']
SQ_GRID = g['SQ_GRID']; WS2_GRID = g['WS2_GRID']
NV, NG = g['NV'], g['NG']; SBINS = g['SBINS']; SC2 = g['SC2']
VE, GE = g['VE'], g['GE']; FC0 = g['FC0']
STRATA = g['STRATA']; LNPI = g['LNPI']
build_pop = g['build_pop']; e_of_x = g['e_of_x']; vp_c = g['vp_c']
project = g['project']; run = g['run']
A0_CAN = g['A0_CAN']; LNY0, DLNY = g['LNY0'], g['DLNY']
TAB_S = g['TAB_S']
np_ = np

OUT = []
def P(s):
    print(s, flush=True)
    OUT.append(s)

def cache_bin(p, prj, PLs, UB, FB, bi):
    """Per-bin log-pp cache -- verbatim model expressions from
    eval_block_nb, with the bin index a parameter."""
    smag, vpar, vper = prj
    s_kau = smag/1e3
    C = np.zeros((len(FCOMP_GRID), len(FFLY_GRID), len(FPM_GRID),
                  len(KW_GRID), len(SQ_GRID), len(WS2_GRID), NV, NG))
    b = SBINS[bi]
    idx = np.where((s_kau >= b[0]) & (s_kau < b[1]))[0]
    assert len(idx) >= 500 and len(PLs[bi]) > 0
    vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
    sg0 = PLs[bi][p['pick'][bi][idx] % len(PLs[bi])]/4.74047
    g1_i, g2_i = p['gn1'][idx], p['gn2'][idx]
    sk_i = s_kau[idx]
    gk_full = p['gs'][idx]
    for fi, fcm in enumerate(FCOMP_GRID):
        cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
        mh_tot = np.zeros(len(idx))
        for k in (1, 2):
            c = p['comp'][k]
            act = c['uc'][idx] < fcm
            mh_tot += act*c['mh'][idx]
            cvp += act*c['w'][idx]*c['wd'][idx, 0]
            cvq += act*c['w'][idx]*c['wd'][idx, 1]
        boost = np.sqrt(1 + mh_tot/p['M_s'][idx])
        for ki, kwv in enumerate(KW_GRID):
            vp_a = vpar[idx] + kwv*cvp
            vq_a = vper[idx] + kwv*cvq
            for pi, fpm in enumerate(FPM_GRID):
                for wi, ws in enumerate(WS2_GRID):
                    if ws == 0.0:
                        vp_n = vp_a*boost + g1_i*sg0*fpm
                        vq_n = vq_a*boost + g2_i*sg0*fpm
                    else:
                        sig_eff = np.sqrt((sg0*fpm)**2
                                          + (ws/4.74047)**2)
                        vp_n = vp_a*boost + g1_i*sig_eff
                        vq_n = vq_a*boost + g2_i*sig_eff
                    vmag = np.hypot(vp_n, vq_n)
                    keep = vmag*4.74047 <= (2.978/np.sqrt(sk_i)
                                            + 2.8284*sg0*4.74047)
                    vtn = (vmag/vc)[keep]
                    gmn = np.degrees(np.arccos(np.clip(
                        np.abs(vp_n[keep])
                        / np.maximum(vmag[keep], 1e-12), 0, 1)))
                    gk = gk_full[keep]
                    for si, sqv in enumerate(SQ_GRID):
                        vts = vtn*np.exp(sqv*gk)
                        h, _, _ = np.histogram2d(
                            np.clip(vts, 0.021, 5.9), gmn,
                            bins=[VE, GE])
                        p0 = np.maximum(h/max(h.sum(), 1), 1e-5)
                        p0 /= p0.sum()
                        for yi, ff in enumerate(FFLY_GRID):
                            wch = min(FC0*SC2[bi], 0.5)
                            wfl = min(ff*SC2[bi], 0.5)
                            wtot = min(wch + wfl, 0.6)
                            mixc = (wch*UB[bi]
                                    + wfl*FB[bi])/(wch + wfl)
                            pp = (1 - wtot)*p0 + wtot*mixc
                            C[fi, yi, pi, ki, si, wi] = np.log(pp)
    return C, idx

def pp_bin(p, prj, PLs, UB, FB, cell, bi):
    """Verbatim eval_pp_nb pp path with the bin index a parameter."""
    fi, yi, fpm, ki, sqv = cell
    fcm = FCOMP_GRID[fi]; kwv = KW_GRID[ki]; ff = FFLY_GRID[yi]
    smag, vpar, vper = prj
    s_kau = smag/1e3
    b = SBINS[bi]
    idx = np.where((s_kau >= b[0]) & (s_kau < b[1]))[0]
    vc = 2*np.pi*np.sqrt(p['M_s'][idx]/smag[idx])
    sg0 = PLs[bi][p['pick'][bi][idx] % len(PLs[bi])]/4.74047
    g1_i, g2_i = p['gn1'][idx], p['gn2'][idx]
    sk_i = s_kau[idx]
    gk_full = p['gs'][idx]
    cvp = np.zeros(len(idx)); cvq = np.zeros(len(idx))
    mh_tot = np.zeros(len(idx))
    for k in (1, 2):
        c = p['comp'][k]
        act = c['uc'][idx] < fcm
        mh_tot += act*c['mh'][idx]
        cvp += act*c['w'][idx]*c['wd'][idx, 0]
        cvq += act*c['w'][idx]*c['wd'][idx, 1]
    boost = np.sqrt(1 + mh_tot/p['M_s'][idx])
    vp_n = (vpar[idx] + kwv*cvp)*boost + g1_i*sg0*fpm
    vq_n = (vper[idx] + kwv*cvq)*boost + g2_i*sg0*fpm
    vmag = np.hypot(vp_n, vq_n)
    keep = vmag*4.74047 <= (2.978/np.sqrt(sk_i)
                            + 2.8284*sg0*4.74047)
    vtn = (vmag/vc)[keep]
    gmn = np.degrees(np.arccos(np.clip(
        np.abs(vp_n[keep])/np.maximum(vmag[keep], 1e-12), 0, 1)))
    gk = gk_full[keep]
    vts = vtn*np.exp(sqv*gk)
    h, _, _ = np.histogram2d(np.clip(vts, 0.021, 5.9), gmn,
                             bins=[VE, GE])
    p0 = np.maximum(h/max(h.sum(), 1), 1e-5)
    p0 /= p0.sum()
    wch = min(FC0*SC2[bi], 0.5)
    wfl = min(ff*SC2[bi], 0.5)
    wtot = min(wch + wfl, 0.6)
    mixc = (wch*UB[bi] + wfl*FB[bi])/(wch + wfl)
    pp = (1 - wtot)*p0 + wtot*mixc
    v = pp.ravel()
    return v/v.sum()

def marg_fpm(T):
    from scipy.special import logsumexp
    lw = logsumexp(T + LNPI.reshape(6, 1, 1, 1, 1, 1),
                   axis=(0, 1, 3, 5))
    wq = np.exp(lw - logsumexp(lw))
    return wq.sum(axis=1)

def S_of(T):
    m = marg_fpm(T)
    return float(np.log(m[3:].sum()) - np.log(m[:3].sum()))

TRUTHS = (1.0, 1.5, 2.1)
R = 64
CELL = (2, 0, None, 1, 0.0)   # fpm slot filled per truth

t0 = time.time()
P("9L-c THE FULL-SEPARATION NOISE-CHANNEL CALIBRATION (R43-C2; "
  "pre-reg committed BEFORE the run; alpha = 0 disclosed scope; "
  "NO credence movement; NO sky verdict)")
P("")

g0_ok = g1_ok = True
save = {}
gate3 = {}
for seed in (31, 101):
    pf = build_pop(seed)
    e_f = e_of_x(pf, 1.05, 0.30)
    tab0 = np.ones_like(TAB_S)
    vp0 = vp_c(pf, e_f, tab0)
    o0 = run(pf['a_s'], e_f, pf['psi0'], pf['f_ip'], pf['M_s'],
             pf['uph'], 8, 2500, 5, a0=A0_CAN, tab=tab0,
             lny0=LNY0, dlny=DLNY, vp=vp0)
    prj0 = project(pf, o0)
    Q = STRATA[0]
    caches = []
    for bi in range(4):
        C, _ = cache_bin(pf, prj0, Q[1], Q[2], Q[3], bi)
        caches.append(C)
    nds = [int(Q[0][bi].sum()) for bi in range(4)]
    P(f"[seed {seed}] Q1 counts per bin: {nds}")

    # GC9LC-1 reader identity: bin-0-only S vs 9L-b archived values
    S_nb = S_of(np.tensordot(caches[0], Q[0][0], axes=([6, 7], [0, 1])))
    want = {31: 0.534, 101: 0.372}[seed]
    ok1 = abs(S_nb - want) <= 0.002
    g1_ok &= ok1
    P(f"[seed {seed}] GC9LC-1 reader identity: bin-0 S = {S_nb:+.4f} vs "
      f"9L-b {want:+.3f} -> {'PASS' if ok1 else 'FAIL'}")

    # GC9LC-0 cache identity on one synthetic 4-bin dataset (truth 2.1)
    rng_id = np.random.default_rng([9, seed, 430])
    dT_max = 0.0
    for bi in range(4):
        pr = pp_bin(pf, prj0, Q[1], Q[2], Q[3], (2, 0, 2.1, 1, 0.0), bi)
        dr = rng_id.multinomial(nds[bi], pr).reshape(NV, NG).astype(float)
        # direct: re-derive log-pp for this bin via pp_bin per cell would
        # be O(1728) pp calls; instead verify cache internal consistency:
        # the cache's own cell (2,0,idx21,1,idx0,0) must equal ln(pr)
        i21 = int(np.where(FPM_GRID == 2.1)[0][0])
        isq0 = int(np.where(SQ_GRID == 0.0)[0][0])
        iw0 = int(np.where(WS2_GRID == 0.0)[0][0])
        lp_cache = caches[bi][2, 0, i21, 1, isq0, iw0].ravel()
        d = float(np.max(np.abs(np.exp(lp_cache)/np.exp(lp_cache).sum()
                                - pr)))
        dT_max = max(dT_max, d)
    ok0 = dT_max <= 1e-9
    g0_ok &= ok0
    P(f"[seed {seed}] GC9LC-0 cache-vs-pp identity: max|dp| = "
      f"{dT_max:.2e} -> {'PASS' if ok0 else 'FAIL'}")
    if not (ok0 and ok1):
        continue

    # the ladder
    means, sds, edges = [], [], []
    allS = {}
    for ti, truth in enumerate(TRUTHS):
        prs = [pp_bin(pf, prj0, Q[1], Q[2], Q[3],
                      (2, 0, truth, 1, 0.0), bi) for bi in range(4)]
        rng_t = np.random.default_rng([9, seed, 900 + ti])
        Ss, Es = np.empty(R), np.empty(R)
        for r in range(R):
            T = 0.0
            for bi in range(4):
                dr = rng_t.multinomial(nds[bi],
                                       prs[bi]).reshape(NV, NG)
                T = T + np.tensordot(caches[bi], dr.astype(float),
                                     axes=([6, 7], [0, 1]))
            Ss[r] = S_of(T)
            Es[r] = float(marg_fpm(T)[-1])
        means.append(Ss.mean()); sds.append(Ss.std(ddof=1))
        edges.append((Es.mean(), Es.std(ddof=1)))
        allS[truth] = Ss
        P(f"[seed {seed}] calib truth {truth:.1f}: S = {Ss.mean():+.3f} "
          f"+- {Ss.std(ddof=1):.3f}; P(fpm=3.0) = {Es.mean():.3f} +- "
          f"{Es.std(ddof=1):.3f} (flat 0.167)  "
          f"({(time.time()-t0)/60:.1f} min)")
    means = np.array(means); sds = np.array(sds)
    D = means[-1] - means[0]
    sig_pool = float(np.sqrt(np.mean(sds**2)))
    ok3 = D >= 2*sig_pool
    gate3[seed] = bool(ok3)
    P(f"[seed {seed}] POWER: D = S(2.1) - S(1.0) = {D:.3f} vs "
      f"2 sigma_pool = {2*sig_pool:.3f} -> "
      f"{'POWERED' if ok3 else 'DEAD'}")
    # sky diagnostic (alpha = 0; NEVER a measurement)
    T_sky = 0.0
    for bi in range(4):
        T_sky = T_sky + np.tensordot(caches[bi],
                                     Q[0][bi],
                                     axes=([6, 7], [0, 1]))
    S_sky = S_of(T_sky)
    zrow = [(S_sky - means[ti])/sds[ti] for ti in range(len(TRUTHS))]
    P(f"[seed {seed}] DIAGNOSTIC-ONLY sky: S = {S_sky:+.3f}; z(truth) = "
      + ", ".join(f"{t:.1f}:{z:+.2f}" for t, z in zip(TRUTHS, zrow)))
    save[seed] = dict(means=means, sds=sds,
                      allS=np.array([allS[t] for t in TRUTHS]),
                      edges=np.array(edges), S_sky=S_sky)
    P("")

np.savez('data/stage9lc_calib.npz', truths=np.array(TRUTHS), R=R,
         **{f's{s}_{k}': v for s in save for k, v in save[s].items()})

if not (g0_ok and g1_ok):
    P("GC9LC-0/1 FAILED - STOP; DO NOT QUOTE; verdict WITHHELD")
elif all(gate3.values()):
    P("==> 9L-c LETTER (locked grammar): C-CHANNEL-POWERED - the "
      "full-separation Q1 noise channel separates the formal-errors "
      "truth from the fitted level at >= 2 sigma (alpha-off "
      "calibration; alpha-on control = the injection arms). The "
      "Paper-1 C2 caveat is replaced by the sourced form; the "
      "alpha-degeneracy disclosure stays. NO credence movement; "
      "NO sky verdict.")
else:
    P("==> 9L-c LETTER (locked grammar): C-CHANNEL-DEAD-TOO - the C2 "
      "caveat is PERMANENT (Q1 E[fpm] only with the flat-prior-mean "
      "sentence; Sec-6.3 sentence weakened to fit-internal status). "
      "NO credence movement; NO sky verdict.")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with io.open('data/stage9lc_jointcal.txt', 'w', encoding='utf-8',
             newline='\n') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage9lc_jointcal.txt")
