"""STAGE 8O — THE T2 EXTERNAL-RV RECONNAISSANCE (pre-reg committed
BEFORE any query).  All 23 ceiling-census pairs (46 components)
crossmatched via CDS XMatch (3") against the scout-verified VizieR
RV tables (APOGEE DR17, GALAH DR3, LAMOST DR7, RAVE DR6) + the Gaia
DR3 RV row.  PRIMARY = the census region (band 9 + above 2 = 11
pairs); the rest are control rows.  Bars + credence map locked in
the pre-reg; classifications from pre-named columns only (G8O-3).
Output: data/stage8o_extrv.txt (pulls cached to data/, committed)
"""
import csv
import io
import os
import time
import urllib.parse
import urllib.request
import uuid

import numpy as np

L_ = []
def P(s):
    print(s, flush=True)
    L_.append(s)

t0 = time.time()
P("8O THE T2 EXTERNAL-RV RECONNAISSANCE (pre-reg committed BEFORE "
  "any query; bars + map locked)")
P("")

# --- G8O-0: targets -------------------------------------------------------
rows = list(csv.DictReader(open('data/ceiling_pairs.csv')))
nb = sum(1 for r in rows if r['band_corr'] == 'band')
na = sum(1 for r in rows if r['band_corr'] == 'above')
ok0 = (len(rows) == 23 and nb == 9 and na == 2)
P(f"G8O-0 targets: {len(rows)} pairs (band {nb}, above {na}) -> "
  f"{'PASS' if ok0 else 'FAIL'}")
assert ok0
targets = []
for i, r in enumerate(rows):
    for cn in (1, 2):
        targets.append(dict(sid=int(r[f'source_id{cn}']), pair=i,
                            comp=cn,
                            primary=r['band_corr'] in ('band', 'above'),
                            cls=r['band_corr'], skau=float(r['s_kAU'])))
sids = [t['sid'] for t in targets]
assert len(set(sids)) == 46 and all(s > 0 for s in sids)

# --- Gaia DR3 positions + RV row (cached) ---------------------------------
def gaia_tap(query):
    data = urllib.parse.urlencode(dict(REQUEST='doQuery', LANG='ADQL',
                                       FORMAT='csv',
                                       QUERY=query)).encode()
    for att in range(3):
        try:
            return urllib.request.urlopen(
                'https://gea.esac.esa.int/tap-server/tap/sync',
                data=data, timeout=180).read().decode()
        except Exception:
            if att == 2:
                raise
            time.sleep(5*(att+1))

TCACHE = 'data/stage8o_targets.csv'
if not os.path.exists(TCACHE):
    q = ("SELECT source_id, ra, dec, phot_g_mean_mag, radial_velocity, "
         "radial_velocity_error FROM gaiadr3.gaia_source "
         "WHERE source_id IN ("
         + ",".join(str(s) for s in sids) + ")")
    open(TCACHE, 'w', newline='').write(gaia_tap(q))
pos, grv = {}, {}
for r in csv.DictReader(open(TCACHE)):
    s = int(r['source_id'])
    pos[s] = (float(r['ra']), float(r['dec']))
    try:
        grv[s] = (float(r['radial_velocity']),
                  float(r['radial_velocity_error']))
    except (ValueError, TypeError):
        grv[s] = None
okp = len(pos) == 46
P(f"G8O-0 Gaia position pull: {len(pos)}/46 -> "
  f"{'PASS' if okp else 'FAIL'}; Gaia RV present for "
  f"{sum(1 for s in sids if grv.get(s))} components")
assert okp

# --- CDS XMatch -----------------------------------------------------------
XM_URL = 'http://cdsxmatch.u-strasbg.fr/xmatch/api/xmatch/sync'
CSV1 = 'sid,ra,dec\n' + '\n'.join(
    f"{s},{pos[s][0]:.8f},{pos[s][1]:.8f}" for s in sids) + '\n'

def xmatch(cat2, dist=3.0):
    bd = uuid.uuid4().hex
    parts = b''
    for k, v in (('request', 'xmatch'), ('distMaxArcsec', str(dist)),
                 ('RESPONSEFORMAT', 'csv'), ('cat2', cat2),
                 ('colRA1', 'ra'), ('colDec1', 'dec')):
        parts += (f'--{bd}\r\nContent-Disposition: form-data; '
                  f'name="{k}"\r\n\r\n{v}\r\n').encode()
    parts += (f'--{bd}\r\nContent-Disposition: form-data; '
              f'name="cat1"; filename="targets.csv"\r\n'
              f'Content-Type: text/csv\r\n\r\n').encode()
    parts += CSV1.encode() + f'\r\n--{bd}--\r\n'.encode()
    req = urllib.request.Request(
        XM_URL, data=parts,
        headers={'Content-Type': f'multipart/form-data; boundary={bd}'})
    for att in range(3):
        try:
            return urllib.request.urlopen(req, timeout=240).read() \
                .decode(errors='replace')
        except Exception:
            if att == 2:
                raise
            time.sleep(8*(att+1))

SURVEYS = [
    ('apogee',  ['vizier:III/286/catalog']),
    ('galah',   ['vizier:J/MNRAS/506/150/stars']),
    ('lamost',  ['vizier:V/156/dr7slrs']),
    ('lamostm', ['vizier:V/156/dr7melrs']),
    ('rave',    ['vizier:III/283/ravedr6', 'vizier:III/283/master']),
]
COLS = {   # G8O-3 pre-named classification columns per survey
    'apogee':  dict(rv='HRV', erv='e_HRV', scat='s_HRV', nv='Nvis'),
    'galah':   dict(rv='RVgalah', erv='e_RVgalah'),
    'lamost':  dict(rv='RV', erv='e_RV'),
    'lamostm': dict(rv='RV', erv='e_RV'),
    'rave':    dict(rv='HRV', erv='e_HRV'),
}
hits = {}      # sid -> list of dicts(survey, ang, rv, erv, scat, nv)
g1_ok, g2_ok = True, True
for name, cands in SURVEYS:
    cache = f'data/stage8o_xm_{name}.csv'
    txt, used = None, None
    if os.path.exists(cache):
        txt = open(cache, encoding='utf-8').read()
        used = '(cache)'
    else:
        for c2 in cands:
            t_ = xmatch(c2)
            if 'angDist' in t_.splitlines()[0] if t_.strip() else False:
                txt, used = t_, c2
                break
            # keep trying candidates; log the miss
            P(f"  [{name}] candidate {c2}: no angDist header - trying "
              f"next" if c2 != cands[-1] else
              f"  [{name}] candidate {c2}: no angDist header - FAILED")
        if txt is not None:
            open(cache, 'w', encoding='utf-8', newline='').write(txt)
    if txt is None:
        P(f"G8O-1 {name}: PULL-FAILED (all candidates) - "
          f"SKIPPED-DISCLOSED")
        g1_ok = False
        continue
    rdr = list(csv.DictReader(io.StringIO(txt)))
    hdr = rdr[0].keys() if rdr else \
        csv.DictReader(io.StringIO(txt)).fieldnames or []
    cols = COLS[name]
    # case-insensitive fallback lookup (renaming logged, not improvised)
    def find(col, hdr=hdr):
        if col in hdr:
            return col
        for h in hdr:
            if h.lower() == col.lower():
                return h
        return None
    cmap = {k: find(v) for k, v in cols.items()}
    missing = [cols[k] for k, v in cmap.items() if v is None
               and k in ('rv', 'erv')]
    P(f"G8O-1 {name} {used}: {len(rdr)} match rows; columns "
      + ", ".join(f"{cols[k]}->{cmap[k]}" for k in cols) )
    if missing:
        P(f"G8O-3 {name}: pre-named column(s) {missing} ABSENT - "
          f"survey SKIPPED-DISCLOSED")
        continue
    for r in rdr:
        try:
            sid = int(float(r['sid']))
            ang = float(r['angDist'])
        except (KeyError, ValueError):
            continue
        if ang > 3.0:
            g2_ok = False
        def fget(key):
            c = cmap.get(key)
            if not c:
                return None
            try:
                return float(r[c])
            except (ValueError, TypeError):
                return None
        rec = dict(survey=name, ang=ang, rv=fget('rv'),
                   erv=fget('erv'), scat=fget('scat'), nv=fget('nv'))
        cur = hits.setdefault(sid, {})
        if name not in cur or ang < cur[name]['ang']:
            cur[name] = rec
P(f"G8O-2 all match distances <= 3 arcsec: "
  f"{'PASS' if g2_ok else 'FAIL'}")
P("")

# --- classification (locked) ----------------------------------------------
def classify(sid):
    hh = hits.get(sid, {})
    if not hh:
        return 'UNCOVERED', []
    tags = []
    for name, rec in hh.items():
        if rec['rv'] is None:
            continue
        if (name == 'apogee' and rec.get('nv') and rec['nv'] >= 2
                and rec.get('scat') is not None):
            if rec['scat'] >= 1.0:
                tags.append(f'VARIABLE({name} s_HRV='
                            f'{rec["scat"]:.2f}, Nvis={rec["nv"]:.0f})')
        if grv.get(sid) and rec['erv'] is not None:
            g_rv, g_erv = grv[sid]
            thr = 3*np.hypot(rec['erv'], g_erv) + 1.0
            dv = abs(rec['rv'] - g_rv)
            if dv > thr:
                tags.append(f'OFFSET({name} |dRV|={dv:.1f} > '
                            f'{thr:.1f} km/s)')
    if tags:
        return 'ACTIVE', tags
    readable = any(rec['rv'] is not None and
                   (grv.get(sid) or (n == 'apogee' and rec.get('nv')
                                     and rec['nv'] >= 2))
                   for n, rec in hh.items())
    return ('QUIET' if readable else 'COVERED-SNAPSHOT'), []

pair_state = {}
P("PER-COMPONENT (primary census region first):")
for t in sorted(targets, key=lambda x: (not x['primary'], x['pair'])):
    st, tags = classify(t['sid'])
    cov = ",".join(sorted(hits.get(t['sid'], {}).keys())) or '-'
    if t['primary'] or st != 'UNCOVERED':
        P(f"  pair{t['pair']:02d}[{t['cls']:>5}] comp{t['comp']} "
          f"sid={t['sid']} s={t['skau']:.1f}kAU: {st:<16} "
          f"surveys=[{cov}]" + ("  " + "; ".join(tags) if tags else ""))
    ps = pair_state.setdefault((t['pair'], t['primary']), [])
    ps.append(st)
P("")

prim_t = [t for t in targets if t['primary']]
C = sum(1 for t in prim_t if t['sid'] in hits)
act_pairs = sorted({pr for (pr, isp), sts in pair_state.items()
                    if isp and 'ACTIVE' in sts})
ctrl_cov = sum(1 for t in targets
               if not t['primary'] and t['sid'] in hits)
ctrl_act = sorted({pr for (pr, isp), sts in pair_state.items()
                   if not isp and 'ACTIVE' in sts})
P(f"COVERAGE: primary components with >= 1 external row: C = {C}/22; "
  f"control components covered: {ctrl_cov}/24")
P(f"ACTIVE pairs: primary {len(act_pairs)} {act_pairs}; control "
  f"{len(ctrl_act)} {ctrl_act}")
P("")
if not g1_ok:
    P("NOTE: one or more survey pulls FAILED - coverage is a lower "
    "bound; disclosed above.")

# --- bars + map (locked) --------------------------------------------------
if C <= 2:
    v = ("T2-ATTEMPTED, NULL-COVERAGE (C <= 2): the public archives "
         "do not cover the census pairs at classification grade; NO "
         "credence movement; T2 remains future-spectrograph work, "
         "now attempted-and-scoped")
elif len(act_pairs) == 0:
    v = (f"RV-CLEAN (C = {C} >= 3, 0 active primary pairs): map +1 "
         f"-> anomaly-real 57 -> 58 (small BY DESIGN: km/s "
         f"thresholds, partial coverage, within-survey-static "
         f"blindness disclosed)")
elif len(act_pairs) == 1:
    v = (f"RV-ACTIVE-MINOR (1 primary pair: {act_pairs}): NO "
         f"credence movement; the pair is annotated in the census "
         f"record")
else:
    v = (f"RV-ACTIVE-MAJOR ({len(act_pairs)} primary pairs: "
         f"{act_pairs}): map -3 -> anomaly-real 57 -> 54; the pairs "
         f"are conditionalized in the census record")
P(f"==> 8O VERDICT (locked bars + map): {v}")
P(f"\ndone ({(time.time()-t0)/60:.1f} min)")

with open('data/stage8o_extrv.txt', 'w') as f:
    f.write("\n".join(L_) + "\n")
print("\nsaved: data/stage8o_extrv.txt")
