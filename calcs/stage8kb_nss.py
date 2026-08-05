"""STAGE 8K-b — THE NSS/ORBIT-CURVATURE CROSSMATCH (pre-reg in
NOTES; committed BEFORE any query; THE decider named by 8L-a2).

FETCH phase (idempotent, cached to data/stage8kb_*.csv): Gaia DR3
TAP sync pulls for three populations — the 9 census pairs (18 ids),
the WIDE population (s >= 6 kAU, 2874 ids), and a size-matched
NARROW control — from gaia_source (ruwe, non_single_star, RV
statistics), nss_two_body_orbit, and nss_acceleration_astro.
READ phase: the selection-aware bars locked in the pre-reg.
Output: data/stage8kb_nss.txt
"""
import csv
import io
import json
import os
import time
import urllib.parse
import urllib.request
import numpy as np
from astropy.io import fits

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

TAP = 'https://gea.esac.esa.int/tap-server/tap/sync'

def tap_query(adql, retries=4):
    data = urllib.parse.urlencode({
        'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv',
        'QUERY': adql}).encode()
    for a in range(retries):
        try:
            req = urllib.request.Request(TAP, data=data)
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode('utf-8', errors='replace')
        except Exception as e:
            if a == retries-1:
                raise
            time.sleep(3*(a+1))

def fetch_table(name, cols, ids, cache, chunk=500):
    if os.path.exists(cache):
        P(f"[fetch] {cache} cached ({sum(1 for _ in open(cache))-1} rows)")
        return
    rows, hdr = [], None
    for i in range(0, len(ids), chunk):
        sub = ids[i:i+chunk]
        adql = (f"SELECT {cols} FROM {name} WHERE source_id IN "
                f"({','.join(str(s) for s in sub)})")
        txt = tap_query(adql)
        lines = [ln for ln in txt.splitlines() if ln.strip()]
        if not lines:
            continue
        if hdr is None:
            hdr = lines[0]
        rows.extend(lines[1:])
    with open(cache, 'w') as f:
        f.write((hdr or cols.replace(' ', '')) + "\n")
        f.write("\n".join(rows) + ("\n" if rows else ""))
    P(f"[fetch] {name}: {len(rows)} rows -> {cache}")

# ---------- populations (loader verbatim) ----------
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
Mtot = np.interp(MG1, MG_T, MS_T)+np.interp(MG2, MG_T, MS_T)
s_kau = sep/1e3
vc = 0.9417*np.sqrt(Mtot/np.maximum(s_kau, 1e-9))
sid1 = d['source_id1'].astype(np.int64)
sid2 = d['source_id2'].astype(np.int64)
WIDE = ok & (s_kau >= 6)
NARR = ok & (s_kau < 2)
iW = np.where(WIDE)[0]
rngc = np.random.default_rng(41)
iN = rngc.choice(np.where(NARR)[0], size=len(iW), replace=False)

cens = [r for r in csv.DictReader(open('data/ceiling_pairs.csv'))
        if r['census_corr'] == 'True']
okidx = np.where(ok)[0]
cidx = []
for r in cens:
    m = (np.abs(s_kau[okidx]-float(r['s_kAU'])) < 0.01) \
      & (np.abs(Mtot[okidx]-float(r['Mtot_Msun'])) < 0.01) \
      & (np.abs(vc[okidx]-float(r['vc_kms'])) < 0.001)
    j = okidx[m]
    assert len(j) == 1
    cidx.append(int(j[0]))
assert len(cidx) == 9
P(f"G8Kb-0: populations - census 9 pairs, WIDE {len(iW)}, "
  f"narrow control {len(iN)} (matched size, rng 41)")

ids_all = sorted(set(
    [int(x) for j in cidx for x in (sid1[j], sid2[j])]
    + [int(x) for j in iW for x in (sid1[j], sid2[j])]
    + [int(x) for j in iN for x in (sid1[j], sid2[j])]))
P(f"G8Kb-0b: {len(ids_all)} unique source_ids to query")

# ---------- FETCH (cached) ----------
fetch_table('gaiadr3.gaia_source',
            'source_id, ruwe, non_single_star, radial_velocity, '
            'radial_velocity_error, rv_nb_transits, '
            'rv_renormalised_gof, rv_chisq_pvalue, '
            'rv_amplitude_robust, rv_template_teff, phot_g_mean_mag',
            ids_all, 'data/stage8kb_src.csv')
fetch_table('gaiadr3.nss_two_body_orbit',
            'source_id, nss_solution_type, period, eccentricity',
            ids_all, 'data/stage8kb_orb.csv')
fetch_table('gaiadr3.nss_acceleration_astro',
            'source_id, nss_solution_type, significance',
            ids_all, 'data/stage8kb_acc.csv')
fetch_table('gaiadr3.nss_non_linear_spectro',
            'source_id, nss_solution_type',
            ids_all, 'data/stage8kb_trend.csv')

# ---------- load pulls ----------
def load_csv(path, key='source_id'):
    out = {}
    with open(path) as f:
        rd = csv.DictReader(f)
        for r in rd:
            try:
                out.setdefault(int(r[key]), []).append(r)
            except (ValueError, KeyError):
                continue
    return out

src = load_csv('data/stage8kb_src.csv')
orb = load_csv('data/stage8kb_orb.csv')
acc = load_csv('data/stage8kb_acc.csv')
trd = load_csv('data/stage8kb_trend.csv')
P(f"pulled: gaia_source {len(src)}, nss_two_body_orbit {len(orb)}, "
  f"nss_acceleration_astro {len(acc)}, nss_non_linear_spectro "
  f"{len(trd)}")
comp = len(src)/len(ids_all)
P(f"G8Kb-1 pull completeness: {len(src)}/{len(ids_all)} = "
  f"{comp:.4f} -> {'PASS' if comp >= 0.99 else 'DISCLOSED (<99%)'}")
P("")

def fnum(r, k):
    try:
        v = float(r[k])
        return v if np.isfinite(v) else np.nan
    except (ValueError, KeyError, TypeError):
        return np.nan

def flags(sid):
    """Per-component channel flags per the pre-reg definitions
    (Katz+23 sect. 3.7 criterion incl. the teff window; Halbwachs+23
    RUWE > 1.4 entry gate for astro-informativeness)."""
    s = src.get(sid, [{}])
    r = s[0] if s else {}
    nss = fnum(r, 'non_single_star')
    ruwe = fnum(r, 'ruwe')
    nbt = fnum(r, 'rv_nb_transits')
    gof = fnum(r, 'rv_renormalised_gof')
    pval = fnum(r, 'rv_chisq_pvalue')
    teff = fnum(r, 'rv_template_teff')
    rv_cov = (np.isfinite(nbt) and nbt >= 10
              and np.isfinite(teff) and 3900 <= teff <= 8000
              and np.isfinite(gof) and np.isfinite(pval))
    rv_var = bool(rv_cov and pval <= 0.01 and gof > 4)
    has_orb = sid in orb
    has_acc = sid in acc
    has_trd = sid in trd
    return dict(nss=nss, ruwe=ruwe, rv_cov=rv_cov, rv_var=rv_var,
                has_orb=has_orb, has_acc=has_acc, has_trd=has_trd,
                astro_inf=bool(np.isfinite(ruwe) and ruwe > 1.4),
                astro_flag=bool(has_orb or has_acc
                                or (np.isfinite(nss)
                                    and int(nss) & 1)))

# ---------- the census-pair leg ----------
P("THE NINE (per pair; COVERED = rv-covered OR astro-informative "
  "(ruwe > 1.4 entry gate); ACTIVE = astro-active OR rv-variable "
  "OR rv-trend):")
n_cov, n_act_cov, n_act = 0, 0, 0
for j in cidx:
    f1, f2 = flags(int(sid1[j])), flags(int(sid2[j]))
    astro = f1['astro_flag'] or f2['astro_flag']
    ainf = f1['astro_inf'] or f2['astro_inf']
    rvcov = f1['rv_cov'] or f2['rv_cov']
    rvv = f1['rv_var'] or f2['rv_var']
    trend = f1['has_trd'] or f2['has_trd']
    covered = rvcov or ainf
    active = astro or rvv or trend
    n_cov += covered
    n_act += active
    n_act_cov += (covered and active)
    P(f"  s = {s_kau[j]:6.2f} kAU: ruwe = ({f1['ruwe']:.2f}, "
      f"{f2['ruwe']:.2f}), nss = ({f1['nss']:.0f}, {f2['nss']:.0f}), "
      f"orb/acc/trend = ({int(f1['has_orb'])+int(f2['has_orb'])},"
      f"{int(f1['has_acc'])+int(f2['has_acc'])},"
      f"{int(f1['has_trd'])+int(f2['has_trd'])}), rv_cov = "
      f"({f1['rv_cov']}, {f2['rv_cov']}), rv_var = ({f1['rv_var']}, "
      f"{f2['rv_var']}) -> {'COVERED' if covered else 'uncovered'}"
      f"{', ACTIVE' if active else ''}")
P(f"census totals: COVERED {n_cov}/9, ACTIVE {n_act}/9, "
  f"active-among-covered {n_act_cov}/{n_cov}")
P("")

# ---------- population rates (descriptive context) ----------
for nm, idxs in (('WIDE', iW), ('NARROW-matched', iN)):
    na = sum(1 for j in idxs
             if flags(int(sid1[j]))['astro_flag']
             or flags(int(sid2[j]))['astro_flag'])
    ncov = sum(1 for j in idxs
               if flags(int(sid1[j]))['rv_cov']
               or flags(int(sid2[j]))['rv_cov'])
    nv = sum(1 for j in idxs
             if flags(int(sid1[j]))['rv_var']
             or flags(int(sid2[j]))['rv_var'])
    P(f"[{nm}] pair rates: astro-active {na}/{len(idxs)} = "
      f"{na/len(idxs):.4f}; RV-covered {ncov}/{len(idxs)} = "
      f"{ncov/len(idxs):.3f}; RV-variable {nv}/{len(idxs)} = "
      f"{nv/len(idxs):.4f}")
P("")

# ---------- verdict per the pre-registered bars (NOTES 8K-b) -----
COV_MIN = 6
if n_cov < COV_MIN:
    v = (f"STRUCTURALLY-BLIND (COVERED {n_cov}/9 < {COV_MIN}): the "
         "in-catalog channels are EXHAUSTED for the census-pair "
         "faker question; hold ~55%; T2 (external RVs of the nine) "
         "is named THE decider; this landing is itself the round's "
         "product")
elif n_act >= 5:
    v = (f"NINE-ACTIVE ({n_act}/9 active at {n_cov}/9 coverage): "
         "the quiet-faker account gains object-level support; per "
         "the map: anomaly-real ~55% -> ~50%; the seed budget "
         "re-prioritizes")
elif n_act_cov <= 2:
    v = (f"NINE-CLEAN ({n_act_cov} active among {n_cov} covered): "
         "the faker account loses its last in-catalog hiding place "
         "at the covered periods; the S3 immunity partially "
         "RESTORES; per the map: anomaly-real ~55% -> ~57% (the "
         "map's ceiling; the P >~ 15-20 yr window disclosure "
         "carried)")
else:
    v = (f"MIXED ({n_act_cov} active among {n_cov} covered): hold "
         "~55%; T2 decides")
P(f"==> 8K-b VERDICT (locked bars + map): {v}")

with open('data/stage8kb_nss.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8kb_nss.txt")
