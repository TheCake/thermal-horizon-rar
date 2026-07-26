# -*- coding: utf-8 -*-
"""Stage 7J-z4: the q-moment conversion audit + the anchor-strength curve
+ the companion-channel attribution (review round 10).

CONTEXT.  Round 10 (the reviewer's close-out) makes two substantive
points.  (1) The knee: the shipped "flat to 0.01-0.03 across every
anchor 0.06-0.34" and the S3 collapse at forced fcomp >= 0.35 coexist
because they live on DIFFERENT AXES - the anchor curve varied the
CENTER at fixed width (sigma = 0.05/0.03), while the collapse is the
sigma -> 0 (hard-forcing) limit.  A smooth sigma = 0.03 anchor centered
at 0.35 would NOT collapse alpha (the likelihood pays the ~30 lnL prior
penalty to keep fcomp = 0.1 and alpha = 0.74; the kinematic gain is
135-153).  The reviewer's letter ("one grid step past the sampled
range") is therefore not the mechanism - but his substance stands: the
flat-curve claim shipped without the axis that carries the risk.  Part
A measures that axis: alpha_marg(sigma) at center = the measured host
rate, locating sigma* = the certificate precision at which the fifth
move actually fires.  (2) The two-moments objection: photometric
detection is high-q weighted (flag |dm| >= 0.4 needs l = q^k >= 0.445,
i.e. q >~ 0.8 before noise smearing) while photocenter wobble
|q/(1+q) - l/(1+l)| VANISHES at q = 1 - so the companions the
photometry counts are the ones that wobble least, and passing host
0.30 -> the model's flat-q fcomp axis as a scalar compares two
differently-weighted moments of pi(q).  TRUE for the wobble channel -
but the reviewer's "kinematically invisible" is INCOMPLETE: the model
carries a second companion channel, the hidden-mass velocity inflation
boost = sqrt(1 + mh/M_s) = sqrt(1 + q/2), which is MAXIMAL at twins
(+22% per-system velocity scale, one-sided, direction-preserving,
s-flat - the sq-shaped channel).  Whether his recount can rescue the
measured rate depends on WHICH channel drives the 135-153 forced-
multiplicity rejection.  Part C attributes it with the cubes' own kw
axis (kw scales ONLY the wobble kicks): cost(fcomp >= 0.35 | kw slice),
within-slice.  Part B computes the moment table with the model's own
verbatim laws (no new data; it prices the CONVERSION, it does not
measure pi_true - the measured-pi object is v2c-plus scope, the joint
(q, P, l) fit the reviewer names as the principled version).

PRE-STATED READINGS (bands fixed before any number was seen):
  Part A: sigma* = largest sigma in the scan at which the seed-law mean
    alpha_marg(sigma) < 0.35 (half-collapse).  Reported per law.
  Part C: Dwob = cost(kw = 1.4) - cost(kw = 0.7), seed-law mean
    (doubling the wobble amplitude, 4x the kick variance).
    Dwob >= +25 -> WOBBLE-BINDING; |Dwob| <= 10 -> COUNT/MASS-BINDING;
    else MIXED.  (Sign matters: negative Dwob = more wobble HELPS the
    forced fit = the strain is not kick-shaped at all.)
  Overall:
    CONVERSION-LIVE = WOBBLE-BINDING and the detected-informed pi
      variants cut the wobble moment to <= 0.5 of flat -> the 3x
      tension is conversion-flexible; the S3 exposure is demoted to
      q-conditional; v2c MUST ship a (q, P)-resolved rate (v2c-plus
      adopted as a requirement, per the review).
    TENSION-ROBUST = COUNT/MASS-BINDING -> the recount cannot rescue
      the measured rate (twins are mass-channel-loud); branch (b)
      stays live as stated in 7J-z3, and the mass channel is named
      the strain driver (adjacent to sq: one-sided, gamma-invisible).
    MIXED = anything else; both carried, the arm suite decides.
  No verdict bar, no credence move (cadence rule; the deciders remain
  v2c + the arm suite - this stage prices which question they answer).

GATES (first run):
  GA0 wide-sigma identity: sigma = 1e3 anchor reproduces the free
      (eta-only) marginal alpha to 0.02 (prior flattens out exactly).
  GA1 anchor-curve regression: sigma = 0.03, center 0.30 reproduces
      the shipped curve values (simple 0.74/+16, BE 0.67/+15) to
      0.02 in alpha / 2.0 in dN (seed means).
  GC0 D2 regression: the no-slice cost reproduces the shipped
      153.3/135.2/153.4/134.8 to 0.05.
  GB0 law identities: l(q=1) = 1 and wfac(q=1) = 0 exactly at every
      M_h; wfac(q -> 0) -> q/(1+q) (dark-companion limit, 5% at
      q = 0.1 where l <= 1e-3).
  Any gate fail -> inspect, do not quote.

AMENDMENTS (logged after the first run FAILED GB0 and before any
number is quoted; both are instrument-design errors of this stage, not
physics changes - no cube or moment value moved):
  A1  GB0's dark-companion clause tested an analytic limit the shared
      implementation INTENTIONALLY lacks: the MS table clips secondary
      masses at its 0.102 M_sun floor, so a q = 0.1 companion keeps
      l ~ 0.02-0.03 and wfac never reaches q/(1+q) (deviation ~26%,
      the printed fail).  The clause is replaced by the law's true
      identities - twin-zero wfac(q=1) = 0 exact, l(q=1) = 1 exact,
      interior maximum of wfac in q (the photocenter law's shape) -
      and the clip floor min-l(q = 0.1) is PRINTED as a measured model
      property instead of gated.
  A2  the fcomp_equiv column as first designed held the host at 0.30
      while re-weighting pi(q), but a pi re-attribution also rescales
      the completeness that produced the 0.30: the decision object is
      the JOINT factor (0.30/completeness-ratio) x moment-ratio, now
      printed as fce_joint.  The pre-stated overall band read the
      moment ratio alone; the FORMAL reading is kept by that letter
      and the joint band is reported alongside as the corrected
      object - both carried to NOTES, no silent re-labeling.

Output: data/stage7j_qmoments.txt
"""
import numpy as np

OUT = []
def P(s):
    print(s); OUT.append(s)

A_GRID = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
E_GRID = np.array([1.05, 1.3])
FCOMP = np.array([0.0, 0.10, 0.20, 0.35, 0.50, 0.70])
KW = np.array([0.7, 1.0, 1.4])
prior_eta = -0.5*((E_GRID-1.3)/0.3)**2
SHIP_D2 = {('simple', 31): 153.3, ('simple', 101): 135.2,
           ('BE', 31): 153.4, ('BE', 101): 134.8}
SHIP_CURVE = {'simple': (0.74, 16.0), 'BE': (0.67, 15.0)}
SIGMAS = [0.05, 0.03, 0.02, 0.015, 0.012, 0.010, 0.0075, 0.005]
CEN = 0.30

def refine(x, y):
    i = int(np.nanargmax(y)); v = x[i]
    if 0 < i < len(x)-1:
        c2, c1, _ = np.polyfit(x[i-1:i+2], y[i-1:i+2], 2)
        if c2 < 0:
            v = -c1/(2*c2)
    return v

def marg_a(cb9, lnpi):
    cbp = cb9 + lnpi.reshape((1, 1, 1, 6, 1, 1, 1, 1, 1))
    m0 = np.nanmax(cbp)
    ex = np.exp(np.nan_to_num(cbp-m0, nan=-np.inf))
    lm = np.log(np.maximum(ex.sum(axis=(1, 2, 3, 4, 5, 6, 7, 8)),
                           1e-300)) + m0
    fpost = ex.sum(axis=(0, 1, 2, 4, 5, 6, 7, 8))
    return refine(A_GRID, lm), float(lm.max()-lm[0]), fpost/fpost.sum()

g_ok = True
P("== PART A: the anchor-strength curve (center %.2f) ==" % CEN)
curves = {}
cost_kw = {}
for law in ('simple', 'BE'):
    for seed in (31, 101):
        cw = np.load(f'data/stage7j_cube_full_photow_{seed}_{law}.npy')
        cb9 = cw + prior_eta.reshape((1, 2, 1, 1, 1, 1, 1, 1, 1))
        # GA0: wide-sigma identity
        am_free, _, _ = marg_a(cb9, np.zeros(6))
        am_wide, _, _ = marg_a(cb9, -0.5*((FCOMP-CEN)/1e3)**2)
        okA0 = abs(am_wide-am_free) <= 0.02
        g_ok &= okA0
        row = []
        for sg in SIGMAS:
            am, dn, fp = marg_a(cb9, -0.5*((FCOMP-CEN)/sg)**2)
            row.append((sg, am, dn, fp))
        curves[(law, seed)] = row
        P(f"[{law} {seed}] GA0 {'PASS' if okA0 else 'FAIL'} | " +
          " ".join(f"s{sg:g}:{am:.2f}/{dn:+.0f}" for sg, am, dn, _ in row))
        P(f"[{law} {seed}]   P(fcomp) mode per sigma: " +
          " ".join(f"{sg:g}:{FCOMP[int(np.argmax(fp))]:.2f}"
                   for sg, _, _, fp in row))
        # GC0 + Part C: kw-sliced forced-multiplicity cost
        cost = float(np.nanmax(cb9) - np.nanmax(cb9[:, :, :, 3:]))
        okC0 = abs(cost - SHIP_D2[(law, seed)]) <= 0.05
        g_ok &= okC0
        cks = []
        for ki in range(3):
            sl = cb9[:, :, :, :, :, :, :, ki:ki+1, :]
            cks.append(float(np.nanmax(sl) - np.nanmax(sl[:, :, :, 3:])))
        cost_kw[(law, seed)] = cks
        P(f"[{law} {seed}] GC0 {'PASS' if okC0 else 'FAIL'} | "
          f"cost(kw=0.7/1.0/1.4) = {cks[0]:.1f}/{cks[1]:.1f}/{cks[2]:.1f}"
          f"  Dwob = {cks[2]-cks[0]:+.1f}")

P("")
for law in ('simple', 'BE'):
    # GA1 vs shipped curve at sigma = 0.03
    am03 = np.mean([[r for r in curves[(law, s)] if r[0] == 0.03][0][1]
                    for s in (31, 101)])
    dn03 = np.mean([[r for r in curves[(law, s)] if r[0] == 0.03][0][2]
                    for s in (31, 101)])
    sa, sd = SHIP_CURVE[law]
    okA1 = abs(am03-sa) <= 0.02 and abs(dn03-sd) <= 2.0
    g_ok &= okA1
    P(f"GA1 {law}: sigma=0.03 mean {am03:.2f}/{dn03:+.1f} "
      f"(shipped {sa}/{sd:+.0f}) {'PASS' if okA1 else 'FAIL'}")
    means = [(sg, np.mean([curves[(law, s)][i][1] for s in (31, 101)]),
              np.mean([curves[(law, s)][i][2] for s in (31, 101)]))
             for i, sg in enumerate(SIGMAS)]
    star = [sg for sg, am, _ in means if am < 0.35]
    P(f"  {law} alpha_marg(sigma): " +
      " ".join(f"{sg:g}:{am:.2f}" for sg, am, _ in means) +
      f"  -> sigma* = {max(star) if star else 'NOT REACHED in scan'}")
    dw = np.mean([cost_kw[(law, s)][2]-cost_kw[(law, s)][0]
                  for s in (31, 101)])
    tag = ('WOBBLE-BINDING' if dw >= 25 else
           'COUNT/MASS-BINDING' if abs(dw) <= 10 else 'MIXED')
    P(f"  {law} Part C: mean Dwob = {dw:+.1f} -> {tag}")

# ---------------- PART B: the moment table (model-verbatim laws) --------
P("")
P("== PART B: q-moment table (model-verbatim; prices the conversion) ==")
MG_T = np.array([2.6, 3.4, 4.2, 4.8, 5.4, 6.0, 6.83, 7.57, 8.16, 8.82,
                 9.29, 10.05, 11.21, 12.45, 14.26])
MS_T = np.array([1.60, 1.33, 1.12, 1.00, 0.90, 0.82, 0.70, 0.64, 0.57,
                 0.50, 0.44, 0.37, 0.23, 0.162, 0.102])
rng = np.random.default_rng(7)
N = 4_000_000
jw = {}
qs = 0.1 + 0.9*rng.random(N)
logP = rng.normal(5.03, 2.28, N)
P_yr = 10**logP/365.25
S = np.minimum(1.0, P_yr/17.8)
okB = True
for M_h in (0.35, 0.50, 0.75):
    a_in = (M_h*(1+qs)*P_yr**2)**(1/3)
    valid = (a_in < 130.0) & (a_in < 10_000.0/5.0)
    v_orb = 29.78*np.sqrt(M_h*(1+qs)/np.maximum(a_in, 1e-3))
    MGp = np.interp(-np.clip(M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
    MGs = np.interp(-np.clip(qs*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
    l_ = 10**(-0.4*(MGs-MGp))
    wfac = np.abs(qs/(1+qs) - l_/(1+l_))
    w = wfac*v_orb*S/4.74047*valid
    dm = 2.5*np.log10(1+l_)
    Wwob = w**2                                # kick variance weight
    dmass = np.sqrt(1+qs/2)-1                  # per-system scale shift
    from scipy.special import erf as verf
    Dsoft = 0.5*(1+verf((dm-0.4)/(0.275*np.sqrt(2))))
    Dhard = (dm >= 0.4).astype(float)
    # GB0 (amendment A1): the law's true identities
    MGs1 = np.interp(-np.clip(1.0*M_h, MS_T[-1], MS_T[0]), -MS_T, MG_T)
    l1 = 10**(-0.4*(MGs1-MGp))
    w1 = abs(1.0/2.0 - l1/(1+l1))
    imax = float(qs[np.argmax(wfac)])
    okB &= abs(l1-1.0) < 1e-12 and w1 < 1e-12 and 0.15 < imax < 0.85
    lo = np.abs(qs-0.1) < 0.005
    P(f"[M_h={M_h:.2f}] clip floor: l(q=0.1) = {l_[lo].mean():.3f} "
      f"(the model's companions are never fully dark); wfac max at "
      f"q = {imax:.2f}; q_min(hard flag) = "
      f"{qs[Dhard > 0].min() if Dhard.sum() else float('nan'):.2f}; "
      f"<Wwob|detected>/<Wwob|undetected> = "
      f"{Wwob[Dsoft > 0.5].mean()/max(Wwob[Dsoft <= 0.5].mean(), 1e-12):.2f}; "
      f"<dmass|detected>/<dmass|undet> = "
      f"{dmass[Dsoft > 0.5].mean()/max(dmass[Dsoft <= 0.5].mean(), 1e-9):.2f}")
    pis = {
        'flat':      np.ones(N),
        'q^-0.5':    qs**-0.5,
        'twin25':    np.where(qs >= 0.9, 1.0 + 0.25*0.9/0.1, 1.0),
        'det-shape': Dsoft.copy(),
    }
    for name, wt in pis.items():
        wt = wt/wt.sum()
        rw = float((wt*Wwob).sum()/(pis['flat']/N*Wwob).sum()*1.0)
        rm = float((wt*dmass**2).sum()/(pis['flat']/N*dmass**2).sum())
        rd = float((wt*Dsoft).sum()/(pis['flat']/N*Dsoft).sum())
        P(f"  pi={name:9s}: completeness x{rd:.2f} | wobble x{rw:.2f} "
          f"(fce_joint {0.30*rw/rd:.3f}) | mass x{rm:.2f} "
          f"(fce_joint {0.30*rm/rd:.3f})")
        jw[(M_h, name)] = (rw, 0.30*rw/rd, 0.30*rm/rd)
g_ok &= okB
P(f"GB0 law identities (amended A1): {'PASS' if okB else 'FAIL'}")

P("")
P(f"GATES: {'ALL PASS' if g_ok else 'FAIL -- do not quote'}")
if g_ok:
    dws = [np.mean([cost_kw[(law, s)][2]-cost_kw[(law, s)][0]
                    for s in (31, 101)]) for law in ('simple', 'BE')]
    wob = all(d >= 25 for d in dws)
    cnt = all(abs(d) <= 10 for d in dws)
    ctag = ('WOBBLE-BINDING' if wob else
            'COUNT/MASS-BINDING' if cnt else 'MIXED')
    # pre-stated letter: the det-shape wobble MOMENT ratio alone
    # (amendment A2 keeps the letter; the joint band is the corrected
    # object, reported alongside)
    rw_det = max(jw[k][0] for k in jw if k[1] == 'det-shape')
    letter = ('CONVERSION-LIVE' if wob and rw_det <= 0.5 else
              'TENSION-ROBUST' if cnt else 'MIXED')
    jb_w = [jw[k][1] for k in jw if k[1] != 'flat']
    jb_m = [jw[k][2] for k in jw if k[1] != 'flat']
    P(f"==> READING: Part C = {ctag} (Dwob {dws[0]:+.1f}/{dws[1]:+.1f});"
      f" FORMAL (pre-stated letter; det-shape moment x{rw_det:.2f} vs"
      f" the 0.5 bar) = {letter}")
    P(f"==> JOINT CONVERSION BAND (amendment A2): fce_joint(wobble) ="
      f" [{min(jb_w):.2f}, {max(jb_w):.2f}], fce_joint(mass) ="
      f" [{min(jb_m):.2f}, {max(jb_m):.2f}] across pi brackets"
      f" (references: kinematic preference 0.10, scalar-passed 0.30).")
    if min(jb_w) <= 0.15:
        s1 = ("the detection-informed bracket REACHES the kinematic "
              "preference (<= 0.15): the wobble-channel tension is "
              "conversion-flexible")
    else:
        s1 = ("no bracket reaches the kinematic preference: the "
              "tension survives every pi tried")
    P(f"    {s1}; the mass-channel joint is pi-stable (completeness "
      f"and mass-weight co-vary) - the certificate must ship "
      f"(q, P)-resolved output either way (v2c-plus, per the "
      f"review's principled version).")

with open('data/stage7j_qmoments.txt', 'w') as f:
    f.write("\n".join(OUT) + "\n")
print("\nsaved: data/stage7j_qmoments.txt")
