# STAGE 6Q: the measurement-ledger audit + world-table consistency assembly.
# Bookkeeping stage - NO fits. Reads LEDGER.csv (hand-curated; see LEDGER.md),
# runs mechanical gates (provenance on disk, supersession resolution, value
# spot-checks against stage outputs), then assembles the world table: every
# candidate law graded against every ledger test, with per-cell provenance
# and data-overlap groups (rows sharing a dataset are not independent votes).
import csv, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LED  = os.path.join(ROOT, 'LEDGER.csv')
OUT  = os.path.join(ROOT, 'data', 'stage6q_worldtable.txt')

L = []
def say(s=''):
    L.append(s); print(s)

rows = list(csv.DictReader(open(LED, encoding='utf-8')))
byid = {r['id']: r for r in rows}
STATUSES = {'CURRENT', 'CO-QUOTED', 'SUPERSEDED', 'RETRACTED'}

say('STAGE 6Q: ledger audit + world table  (%d ledger rows)' % len(rows))
say('=' * 78)

# ---- G1: ids unique, statuses legal, required fields present -----------------
fails = []
if len(byid) != len(rows):
    fails.append('duplicate ids')
for r in rows:
    if r['status'] not in STATUSES: fails.append('bad status %s: %s' % (r['id'], r['status']))
    if not (r['id'] and r['quantity'] and r['value'] and r['stage']):
        fails.append('missing required field on %s' % r['id'])
say('G1 (ids/status/fields): %s' % ('PASS' if not fails else 'FAIL ' + '; '.join(fails)))

# ---- G2: provenance exists on disk ------------------------------------------
miss = []
for r in rows:
    for f in (r['script'], r['output']):
        if f and not os.path.exists(os.path.join(ROOT, f.replace('/', os.sep))):
            miss.append('%s -> %s' % (r['id'], f))
say('G2 (script/output on disk): %s' % ('PASS' if not miss else 'FAIL ' + '; '.join(miss)))

# ---- G3: supersession pointers resolve to live rows -------------------------
bad = []
for r in rows:
    tgt = r['superseded_by']
    if r['status'] == 'SUPERSEDED' and not tgt:
        bad.append('%s superseded with no pointer' % r['id'])
    if tgt:
        if tgt not in byid: bad.append('%s -> missing target %s' % (r['id'], tgt))
        elif byid[tgt]['status'] not in ('CURRENT', 'CO-QUOTED'):
            bad.append('%s -> non-live target %s (%s)' % (r['id'], tgt, byid[tgt]['status']))
say('G3 (supersession resolves): %s' % ('PASS' if not bad else 'FAIL ' + '; '.join(bad)))

# ---- G4: no two CURRENT rows claim the same quantity ------------------------
seen, dup = {}, []
for r in rows:
    if r['status'] == 'CURRENT':
        q = r['quantity'].strip().lower()
        if q in seen: dup.append('%s / %s' % (seen[q], r['id']))
        seen[q] = r['id']
say('G4 (current-quantity collisions): %s' % ('PASS' if not dup else 'FAIL ' + '; '.join(dup)))

# ---- G5: spot-check ledger values against the stage outputs -----------------
CHECKS = [
    ('bin-boost-corr',  '1.078',  'data/stage4q_perspective.txt'),
    ('bin-boost-perp',  '1.151',  'data/stage4q_perspective.txt'),
    ('gal-p-marg',      '0.578',  'data/stage4h_p_ml.txt'),
    ('gal-c1-flat',     '0.450',  'data/stage4s_c1fit.txt'),
    ('gal-c1-hier',     '0.258',  'data/stage4z_hierc1.txt'),
    ('boot-amb',        '-56.71', 'data/stage6j_ambboot.txt'),
    ('boot-amb',        '37/40',  'data/stage6j_ambboot.txt'),
    ('boot-f4',         '-57.4',  'data/stage6c_f4boot.txt'),
    ('boot-gm',         '-29.27', 'data/stage5n_dvboot.txt'),
    ('bin-sflat-fact',  '+37.00', 'data/stage6pb_extend.txt'),
    ('galfn-amb',       '-59.05', 'data/stage6i_chaegate.txt'),
    ('galfn-f3f4',      '-64.2',  'data/stage6a_twogal.txt'),
    ('binfn-amb',       '-0.88',  'data/stage6g_verdict.txt'),
    ('binfn-drive',     '-9.64',  'data/stage6d_verdict.txt'),
    ('gal-legcount',    '10/40',  'data/stage6l_legboot.txt'),
    ('galfn-resn',      '-113.72', 'data/stage6s_resngal.txt'),
    ('binfn-resn',      '-7.38',  'data/stage6t_verdict.txt'),
    ('mech-gate',       '0.6884', 'data/stage6u_gatederiv.txt'),
    ('gal-untied',      '+3.33',  'data/stage6v_untied.txt'),
    ('mech-borrow',     '+0.989', 'data/stage6x_borrow.txt'),
    ('binfn-scalarefe', '-10.35', 'data/stage6w_verdict.txt'),
    ('mech-reservoir',  '0.9524', 'data/stage6y_reservoir.txt'),
    ('gal-ordering',    '-62.91', 'data/stage6z_ordering.txt'),
    ('galfn-shotbath',  '-25.05', 'data/stage7a_einstein.txt'),
    ('gal-curve-coherence', '+0.8761', 'data/stage7b_bumphunt.txt'),
    ('gal-bump-inner',  '0.0000', 'data/stage7c_gammaclean.txt'),
    ('gal-gamma-final', '-2.376', 'data/stage7c_gammaclean.txt'),
    ('mech-platform',   '+0.987', 'data/stage7e_platform.txt'),
    ('galfn-qcl',       '+556.51', 'data/stage7d_vacuumshare.txt'),
    ('binfn-rja',       '+1.68', 'data/stage7d_vacuumshare.txt'),
    ('mech-mixmean',    '0.7445', 'data/stage7f_mixmean.txt'),
    ('sol-trajmargin',  '451', 'data/stage7g_trajsaturn.txt'),
    ('binfn-ambmi',     '-8.43', 'data/stage7g_trajsaturn.txt'),
    ('binfn-ambmiavg',  '-2.92', 'data/stage7h_miavg.txt'),
    ('bin-ablation-wrad',   '-0.010', 'data/stage7i_verdict.txt'),
    ('bin-ablation-strict', '-0.740', 'data/stage7i_verdict.txt'),
    ('bin-ablation-strict', 'COMPANION-DIRECTION', 'data/stage7i_verdict.txt'),
    ('bin-ceilingpairs',    '11/11 reproduced', 'data/stage7i_verdict.txt'),
    ('bin-7j-completeness', 'f_host in [0.42, 0.57] (peak 0.51)',
     'data/stage7j_completeness.txt'),
    ('bin-7j-completeness', 'C(companion) = 0.410',
     'data/stage7j_completeness.txt'),
    ('bin-7j-marginal',     '==> COMPANION-WIN', 'data/stage7j_verdict.txt'),
    ('bin-7j-marginal',     'a_marg=0.00', 'data/stage7j_full_photo.txt'),
    ('bin-7j-rawcond',      '-1034.4', 'data/stage7j_diag.txt'),
    ('bin-7j-lowend',       'a_marg=0.37', 'data/stage7j_lowend.txt'),
    ('bin-7j-lowend',       'lit-prior marg: a_marg=0.39',
     'data/stage7j_lowend.txt'),
    ('ret-73-manufactured', 'f_host in [0.42, 0.57] (peak 0.51)',
     'data/stage7j_completeness.txt'),
    ('bin-7j-paircorr',     'rho_att=+0.465', 'data/stage7j_paircorr.txt'),
    ('bin-7j-paircorr',     '==> C-FAIL', 'data/stage7j_paircorr.txt'),
    ('bin-7j-paircorr',     'C-FAIL STANDS', 'data/stage7j_paircorr.txt'),
    ('bin-7j-armdiag',      'PROF=1.47/+71.1  meas=1.48/+94.5',
     'data/stage7j_armdiag.txt'),
    ('bin-7j-anchorcurve',  '0.06:0.31/+2', 'data/stage7j_lowend.txt'),
    ('bin-7j-seed6', 'mean a_marg = 0.30 +- 0.09 SE', 'data/stage7j_seed6.txt'),
    ('bin-7j-seed6', 'AMENDMENT-6 SEED RULE: no break', 'data/stage7j_seed6.txt'),
    ('bin-7jz-v1', 'anchor NOT shipped', 'data/stage7jz_mixture.txt'),
    ('bin-7jz-v1', '[0.15, 0.22]', 'data/stage7jz_mixture.txt'),
    ('bin-7jz-v2b', '[0.29, 0.32]', 'data/stage7jz2b_exact.txt'),
    ('bin-7jz-v2b', 'f_hat=0.000', 'data/stage7jz2b_exact.txt'),
    ('bin-7jz-width', 'a_marg=0.74', 'data/stage7jz_read_lit.txt'),
    ('bin-7jz-width', 'AMBIGUOUS-CARRIED (LIT-CONDITIONAL)',
     'data/stage7jz_read_lit.txt'),
    ('bin-7jg-gamma', 'ABSORBER-LEVEL SEPARATION', 'data/stage7jg_read.txt'),
    ('bin-7jg-gamma', 'SEPARATION-CONFIRMED (alpha grade; M3(simple))',
     'data/stage7jg_read.txt'),
    ('bin-7j-sqclose', 'S1 mean anchor strain = 0.0',
     'data/stage7j_sqclose.txt'),
    ('bin-7j-sqclose', 'FIFTH-MOVE-LIVE', 'data/stage7j_sqclose.txt'),
    ('bin-7j-qmoments', 'sigma* = 0.02', 'data/stage7j_qmoments.txt'),
    ('bin-7j-qmoments', 'fce_joint(wobble) = [0.10, 0.39]',
     'data/stage7j_qmoments.txt'),
    ('bin-7jz-v2c', 'GATES ALL PASS -> THE CERTIFICATE SHIPS',
     'data/stage7jz2c_cert.txt'),
    ('bin-7jz-v2c', 'f_host per component in [0.29, 0.32] (peak 0.32)',
     'data/stage7jz2c_cert.txt'),
    ('bin-7jz-qshape', 'twin t=5: dlnL=+0.0, f=0.124',
     'data/stage7jz2c_cert.txt'),
    ('bin-7jz-qshape', 'flat: dlnL=-162.0', 'data/stage7jz2c_cert.txt'),
    ('bin-7jz-anchored', 'VERDICT @ LANDED-CONV: AMBIGUOUS-CARRIED',
     'data/stage7jz_read.txt'),
    ('bin-7jz5-earm', 'D3-EXT STILL-RIDING', 'data/stage7jz5_eread.txt'),
    ('bin-7jz5-earm', 'P(fpm=3.0)=0.97', 'data/stage7jz5_eread.txt'),
    ('bin-7j-anchorcurve',  'hard-wall reference (f <= 0.1 fence, 4R): '
     'a_hat = 1.06, dN = +99.5', 'data/stage7j_lowend.txt'),
    ('bin-7j-armdiag',      '[inj-fullpow] lnL(fcomp)-max = [-4726.2',
     'data/stage7j_armdiag.txt'),
    ('ret-7i-median-immune', 'completeness of the -0.4 flag: C(companion) = 0.410',
     'data/stage7j_completeness.txt'),
]
g5bad = []
for rid, tok, f in CHECKS:
    p = os.path.join(ROOT, f.replace('/', os.sep))
    try:
        txt = open(p, encoding='utf-8', errors='replace').read()
    except OSError:
        g5bad.append('%s: no file %s' % (rid, f)); continue
    if tok not in txt:
        g5bad.append('%s: "%s" NOT in %s' % (rid, tok, f))
say('G5 (value spot-checks, %d tokens): %s' % (len(CHECKS),
    'PASS' if not g5bad else 'FAIL  ' + '; '.join(g5bad)))
say('')

# ---- The world table --------------------------------------------------------
# Columns: BIN  = v7 2D binary contest vs BE          (EDR3-14071)
#          GALH = hierarchical galaxy ladder vs BE     (SPARC-153 [+Chae21])
#          GBOOT= paired galaxy bootstrap grade        (same data as GALH)
#          A0   = binary a0 translation vs cH0/2pi     (EDR3-14071)
#          C1   = exact c1 vs the measured dials (gal 0.26-0.45 open to 1/2 at
#                 flat grade; bin 0.37-0.50)
#          TAILP= exact tail p vs measured (gal 0.65-0.75; bin ~0.5)
#          QUAD = solar quadrupole vs Cassini
# Grades: LEAD / PASS / TIE / LEAN- (trails at SE grade) / REJECT (strong lean)
#         / VETO (shape rejection or 0/N) / DEAD / SHARED-FAIL (the 4.0-5.8x
#         Cassini lock common to every MG member) / ESCAPE / OPEN / na
W = [
 ('Newton', {
   'BIN':  ('DEAD',  '+90..+103 lnL behind, min +84; 24/24 contests', 'bin-newton-v7corr'),
   'GALH': ('DEAD',  '+1659/+2777 on 15 lensing points alone; loses every SPARC fit', 'gal-newton-lensing'),
   'GBOOT':('DEAD',  'excluded 2000/2000, min +53', 'bin-newton-final'),
   'A0':   ('na',    '', ''), 'C1': ('na', 'nu=1: c1 undefined; c1=0 branch dead-to-strong-lean', 'gal-nlo-honest'),
   'TAILP':('na', '', ''), 'QUAD': ('PASS', 'no anomalous quadrupole', 'sol-quadrupole')},
  'DEAD on the sky; alive only at Saturn'),
 ('simple-nu (classical bath)', {
   'BIN':  ('LEAD',  'beats BE by 12.6+-2.4 SE', 'bin-lean-simple'),
   'GALH': ('REJECT','falls LAST hierarchically: 99 behind BE (strong lean; no bootstrap)', 'galfn-simple-hier'),
   'GBOOT':('na',    'flat SPARC agnostic (raw-chi2 lean retracted)', 'ret-sparc-simple'),
   'A0':   ('PASS',  '+2.5 sigma', 'bin-a0-half'),
   'C1':   ('OK',    'c1=1/2 exact (c2=1/8)', 'gal-c1-hier'),
   'TAILP':('na',    'power-law tail; outside the p-family', ''),
   'QUAD': ('SHARED-FAIL', '4.3x Cassini', 'sol-quadrupole')},
  'binary-favored, hier-galaxy-rejected: the two-system tension in one row'),
 ('BE (C&T reference)', {
   'BIN':  ('PASS',  'the reference; beats every sharpened function here', 'binfn-pref-half'),
   'GALH': ('REJECT','trailed by every sharpened function (-42..-64 vertical)', 'galfn-f3f4'),
   'GBOOT':('REJECT','behind F4 and AMB 37/40 (lean grade)', 'boot-amb'),
   'A0':   ('PASS',  '+1.9 sigma', 'bin-a0-half'),
   'C1':   ('OK',    'c1=1/2, c2=1/12 exact', 'gal-c1-hier'),
   'TAILP':('MIXED', 'p=1/2: binaries OK, galaxies want 0.65-0.75', 'gal-tail-p'),
   'QUAD': ('SHARED-FAIL', '4.3x Cassini', 'sol-quadrupole')},
  'binary-anchored baseline; galaxies have outgrown it at lean grade'),
 ('nu_p(0.65)', {
   'BIN':  ('LEAN-', 'trails BE by 5.2+-0.9; interior = accepted', 'binfn-p065'),
   'GALH': ('LEAD',  '-56.4 plain (tail dial alone); -32 vertical', 'gal-tail-p'),
   'GBOOT':('na', '', ''),
   'A0':   ('REJECT','+4.9 sigma (sharp translation)', 'bin-a0-sharp'),
   'C1':   ('na',    'p-family axis, not a c1 member', ''),
   'TAILP':('MIXED', 'p=0.65: galaxies OK, binaries want 1/2', 'binfn-pref-half'),
   'QUAD': ('SHARED-FAIL', 'amplitude-locked', 'sol-quadrupole')},
  'first both-systems-viable function; superseded in role by AMB'),
 ('gm (beta=1/2)', {
   'BIN':  ('LEAN-', 'trails BE by 8.5+-2.1', 'binfn-gm'),
   'GALH': ('LEAD',  '-84.8 plain / -42.7 vertical', 'galfn-gm'),
   'GBOOT':('LEAN+', '-29.27+-52.98, 29/40', 'boot-gm'),
   'A0':   ('REJECT','+4.3 sigma', 'bin-a0-sharp'),
   'C1':   ('MIXED', 'c1=1/3: galaxy dial OK; binary dial 0.37-0.50 excludes it narrowly', 'bin-c1'),
   'TAILP':('MIXED', 'p=3/4: galaxy edge, binaries reject', 'binfn-pref-half'),
   'QUAD': ('SHARED-FAIL', '4.0-5.8x (beta-blind)', 'sol-quadrupole')},
  'best zero-param galaxy function of the mixing family; binary-strained'),
 ('F1/F2 (spontaneous fraction)', {
   'BIN':  ('LEAN-', 'trail BE by 5.8/5.5+-2; interior 12/12', 'binfn-f1f2'),
   'GALH': ('LEAD',  'F2 -107.2/-51.9; F1 -101.3/-44.3', 'galfn-f1f2'),
   'GBOOT':('na', '', ''),
   'A0':   ('REJECT','+5.2..+6.3 sigma', 'bin-a0-sharp'),
   'C1':   ('OK',    'c1=1/2 exact; p=3/4 exact', 'gal-c1-hier'),
   'TAILP':('MIXED', 'p=3/4: galaxy edge, binaries reject', 'binfn-pref-half'),
   'QUAD': ('SHARED-FAIL', '4.8-5.3x', 'sol-quadrupole')},
  'derived not fitted; improved-not-accepted on binaries; a0 wall'),
 ('F3 (two-leg)', {
   'BIN':  ('LEAN-', 'trails BE by 6.84; interior', 'binfn-f3f4'),
   'GALH': ('LEAD',  '-111.4 plain / -50.9 vertical', 'galfn-f3f4'),
   'GBOOT':('na', '', ''),
   'A0':   ('REJECT','+6.2 sigma grade (kappa 1.26-1.35)', 'bin-a0-sharp'),
   'C1':   ('OK',    'c1=1/2 AND c2=1/12 exact', 'gal-c1-hier'),
   'TAILP':('na',    'sharp-screening class (the rejected behavior is screening-region)', ''),
   'QUAD': ('SHARED-FAIL', '', 'sol-quadrupole')},
  'galaxy-strong, binary- and a0-strained'),
 ('F4 (two-leg)', {
   'BIN':  ('VETO',  'alpha-EDGE 6/6 = shape rejection (trails 6.96)', 'binfn-f3f4'),
   'GALH': ('LEAD',  '-108.7/-64.2 = biggest controlled galaxy lead', 'galfn-f3f4'),
   'GBOOT':('LEAN+', '-57.4+-38.3, 37/40', 'boot-f4'),
   'A0':   ('na',    'edge-invalidated', 'bin-a0-sharp'),
   'C1':   ('OK',    'c1=1/2 AND c2=1/12 exact', 'gal-c1-hier'),
   'TAILP':('na', '', ''),
   'QUAD': ('SHARED-FAIL', '', 'sol-quadrupole')},
  'best pure-galaxy function; binary-vetoed'),
 ('boot (quarter cell)', {
   'BIN':  ('VETO',  '+17-24 behind half-branch; alpha edge-rides = shape rejection', 'binfn-boot'),
   'GALH': ('VETO',  'hier lead +75.6 COLLAPSES to -9 under the measured vertical channel', 'galfn-boot-hier'),
   'GBOOT':('na',    '43/50 was the no-vertical conditional', 'galfn-boot-hier'),
   'A0':   ('na', '', ''),
   'C1':   ('MIXED', 'c1=1/4: hier profile peaks there but flat-joint disfavors 3.1 sigma; binary dial excludes', 'gal-c1-flat'),
   'TAILP':('na',    'nu(1)=1.35 weak transition = the rejected shape', ''),
   'QUAD': ('SHARED-FAIL', '', 'sol-quadrupole')},
  'DEAD everywhere once vertical freedom and binaries are counted'),
 ('AMB (ambient-gated, L=2)', {
   'BIN':  ('TIE',   '-0.88+-2.66 vs BE; interior 6/6; ACCEPTED', 'binfn-amb'),
   'GALH': ('LEAD',  '-61.68 vertical / -100.5..-105.9 plain (measured ambients)', 'galfn-amb'),
   'GBOOT':('LEAN+', '-56.71+-35.65, 37/40 = top grade, binary-compatible', 'boot-amb'),
   'A0':   ('PASS',  '+1.6 sigma = best temperature row in the program', 'bin-a0-amb'),
   'C1':   ('OK',    'c1=1/2, c2=1/12 gate-independent', 'gal-c1-hier'),
   'TAILP':('OK',    'postdicts BOTH: 0.689 gal / 0.529 bin', 'amb-p-postdict'),
   'QUAD': ('SHARED-FAIL', '3.60e-26 = 4.0x', 'sol-quadrupole')},
  'the unique two-system pass; post-hoc flag carried; DR4 = out-of-sample decider'),
 ('RESN (resolution bath)', {
   'BIN':  ('VETO',  'SHAPE REJECTION: alpha edge 5/6, 0/6 prefer, -7.38+-2.08 = the eight-function band', 'binfn-resn'),
   'GALH': ('LEAD',  '-58.59 vertical (STRONG bar) / -113.72 plain = largest plain lead on record', 'galfn-resn'),
   'GBOOT':('na', '', ''),
   'A0':   ('na',    'edge-invalidated (kappa 1.351; formally +10 sigma at the edge value)', 'binfn-resn'),
   'C1':   ('OK',    'c1=1/2..c4=-1/720 ALL preserved; break c5=-1/16 (deepest preservation)', 'mech-resolution'),
   'TAILP':('MIXED', 'p=3/4 exact (gm argument): galaxy edge, binaries reject', 'mech-resolution'),
   'QUAD': ('na',    'not separately computed (amplitude-lock pattern expects ~4-5x)', 'sol-quadrupole')},
  'derived pre-hoc from the 6N corollary; galaxy-best-in-class, binary-vetoed: proves the ambient gate REQUIRED'),
 # note: the AMB measured-ambient gains are heterogeneity-generic per 6Z
 # (correction #14); its amplitude-level legs are the load-bearing ones.
 ('drive-weighted (pointwise)', {
   'BIN':  ('VETO',  'EXCLUDED -9.64+-1.49, 0/6; predicted mid-separation sag observed', 'binfn-drive'),
   'GALH': ('na',    'isolated limit reduces to F4 (gate-verified)', 'binfn-drive'),
   'GBOOT':('na', '', ''),
   'A0':   ('PASS',  '+2.1 sigma (the two vetoes are independent)', 'binfn-drive'),
   'C1':   ('OK',    'c1=1/2 isolated', ''),
   'TAILP':('na', '', ''), 'QUAD': ('na', 'not separately computed', '')},
  'the local-field version of the gate is excluded; the gate is system-level'),
 ('MI (EFE-respecting, mi_t)', {
   'BIN':  ('TIE',   'ties MG: -3.5+-3.3 / -0.8+-2.5', 'binfn-mi'),
   'GALH': ('na',    'MI = MG on circular orbits (rotation curves blind)', 'binfn-mi'),
   'GBOOT':('na', '', ''), 'A0': ('na', 'shares the MG tables', ''),
   'C1':   ('na', '', ''), 'TAILP': ('na', '', ''),
   'QUAD': ('ESCAPE','no capped-type EFE quadrupole', 'sol-quadrupole')},
  'the open door through Saturn; DR4 eccentricity-resolved data decides'),
 ('MI (no-EFE)', {
   'BIN':  ('VETO',  'dead 12/12: -20..-28 with alpha driven to 0.5-0.6 = the data demand the EFE amplitude', 'binfn-mi'),
   'GALH': ('na', '', ''), 'GBOOT': ('na', '', ''), 'A0': ('na', '', ''),
   'C1':   ('na', '', ''), 'TAILP': ('na', '', ''),
   'QUAD': ('ESCAPE', '', 'sol-quadrupole')},
  'dead'),
]

# G6: every world-table cell cites a live ledger row --------------------------
g6bad = []
for name, cells, verdict in W:
    for col, (grade, num, src) in cells.items():
        if src and src not in byid:
            g6bad.append('%s/%s -> %s' % (name, col, src))
        elif src and byid[src]['status'] == 'RETRACTED' and grade not in ('na',):
            g6bad.append('%s/%s cites RETRACTED %s' % (name, col, src))
say('G6 (world-table provenance): %s' % ('PASS' if not g6bad else 'FAIL ' + '; '.join(g6bad)))
say('')

say('THE WORLD TABLE  (columns: BIN binary contest | GALH hier-galaxy ladder |')
say('  GBOOT galaxy bootstrap | A0 binary temperature | C1 dial | TAILP dial |')
say('  QUAD Cassini)')
say('-' * 78)
for name, cells, verdict in W:
    say(name)
    for col in ('BIN', 'GALH', 'GBOOT', 'A0', 'C1', 'TAILP', 'QUAD'):
        grade, num, src = cells[col]
        if grade == 'na' and not num: continue
        say('  %-6s %-12s %s%s' % (col, grade, num, ('  [%s]' % src) if src else ''))
    say('  => %s' % verdict)
    say('')

# ---- Viability + independence summary ---------------------------------------
say('-' * 78)
vetoed = [n for n, c, v in W if any(g in ('VETO', 'DEAD') for g, _, _ in c.values())]
alive  = [n for n, c, v in W if not any(g in ('VETO', 'DEAD') for g, _, _ in c.values())]
say('Formal sky vetoes (VETO/DEAD in some row): ' + '; '.join(vetoed))
say('No formal sky veto: ' + '; '.join(alive))
say('')
say('Independence: all GALH/GBOOT rows share SPARC-153 (+Chae21 ambients) = ONE')
say('data vote; all BIN/A0 rows share EDR3-14071 = ONE data vote; QUAD (Cassini)')
say('is independent and 4.0-5.8x AGAINST every MG member (escape: EFE-respecting')
say('MI, which ties the binary contest). "Passes both systems" therefore means')
say('exactly two independent data votes, plus one shared solar tension.')
say('')
say('7J CONDITIONALITY BANNER (through the 7J-z/7J-g deciders, 2026-07-26):')
say('every BIN/A0 row above was measured inside the multiplicity fence')
say('(f_comp <= 0.1) and WITHOUT the per-system width channel.  The fenced')
say('+99-110 Newton rejection stays DEAD.  The width-channel round')
say('(bin-7jz-width): once sq is a parameter (sq = 0.2 interior, demanded')
say('at P = 1.00 - the 3E/6P object), alpha_marg = 0.70-0.74 at')
say('dN = +23.2-23.8 with the anchor curve FLAT to 0.01-0.03 across every')
say('companion anchor 0.06-0.34 - the alpha answer is companion-prior-')
say('INDEPENDENT and the old COMPANION-WIN lives in the sq = 0 slice.')
say('Verdict label AMBIGUOUS-CARRIED at the LANDED operative anchor')
say('(missed the pre-registered detection bar by 1.2-1.8 lnL; the v2c')
say('certificate SHIPPED: rate six-way stable 0.159, host [0.29, 0.32]')
say('flat-convention; the GV7 q-shape table measures the subsystems')
say('TWIN-HEAVY (t=5 beats flat +162, f 0.124, host ~0.23) = the')
say('wobble-quiet population - the x3 multiplicity tension DISSOLVING')
say('in the round-10 direction; the FACE flat-q read = the fifth move')
say('fired live, alpha 0.00, quoted as the conditional whose premise')
say('the q-table rejects).  E-read (bin-7jz5-earm): the fpm edge')
say('CHASES to 3.0 (BE P=0.97) = an unphysical noise demand (Lindegren')
say('ceiling ~1.4) = the missing-width-SHAPE signature; no further')
say('extension; operative numbers co-quoted with the band alpha')
say('0.68-0.74 / dN +14.5-23.8; A-D mismatch arms = next.')
say('7J-g: gamma pins the absorbers (SD(w_rad) -> below one 0.10 grid')
say('step, all mass on the 0.20 node); in the four-absorber configuration')
say('vt-only fitting reports a PHANTOM alpha = 0.5 that the direction')
say('data veto (SEPARATION-CONFIRMED, M3).  Round-9 decomposition')
say('(bin-7j-sqclose): anchor strain 0.0 - the quoted alpha is the')
say('width-complete model\'s own optimum - but forcing the measured ~0.3')
say('host rate collapses alpha to 0.00 in 4/4 reads: FIFTH-MOVE-LIVE,')
say('alpha fully exposed to the multiplicity tension.  Round-10')
say('conversion audit (bin-7j-qmoments) RESHAPES the exposure: the')
say('collapse needs anchor width sigma <= 0.02 (sigma* measured; alpha')
say('holds 0.67-0.75 at sigma >= 0.03), the rejection is WOBBLE-BINDING')
say('(+310 per amplitude doubling), and the joint q-conversion band')
say('[0.10, 0.39] brackets the kinematic preference - detection and')
say('wobble are 6-8x anticorrelated in q, so the x3 multiplicity')
say('tension is NOT YET KINEMATIC: the certificate must ship a')
say('(q,P)-resolved rate (v2c-plus, required).  The binary DATA VOTE')
say('stays suspended pending v2c-plus + the anchored re-read + the arm')
say('suite (degraded injections, fpm -> 3.0, kw-residual attribution);')
say('7K and 7J-d queue behind.  Galaxy rows are untouched.')
say('')
say('The mechanical sentence the table proves: AMB is the only candidate that is')
say('(a) un-vetoed on the binaries (tie with BE), (b) at the top galaxy bootstrap')
say('grade (37/40, lean not detection), (c) best on the binary temperature row')
say('(+1.6 sigma), and (d) postdicts both measured tail exponents - carrying the')
say('post-hoc flag and the same Cassini quadrupole as every other MG member.')
say('The ladder digit (1/2 vs 1/3 vs sharper) remains OPEN at population grade.')
say('')
say('The 6R/6T triangulation (three independent exclusions force the two-factor')
say('grammar): local-resolution running alone = binary-dead (RESN); pointwise')
say('ambient weighting alone = binary-dead (drive-weighted); local quantumness x')
say('ambient classicality = the only structure passing both systems (AMB).')

os.makedirs(os.path.dirname(OUT), exist_ok=True)
open(OUT, 'w', encoding='utf-8').write('\n'.join(L) + '\n')
print('\n[saved -> %s]' % OUT)
