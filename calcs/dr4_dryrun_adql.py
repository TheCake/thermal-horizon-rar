"""DR4 day-one ADQL dry-run (red-team finding F10, with F7/F8 probes).

The DR4 day-one plan (DR4.md) contains no data-access step: no ADQL, no table
names, no volume estimate.  The realistic day-one state is ~2.8e9 sources, no
community wide-binary pair catalogue (the EDR3 -> El-Badry/Rix/Heintz lag was
months), and a loaded archive.  This script writes the queries NOW and dry-runs
them against the PUBLIC Gaia DR3 TAP service, which is the same grammar modulo
table/column names.  The DR3 row counts become the day-one regression target.

Probes (all sync mode, small TOP limits, no bulk pulls):

  Q1  patch source census            -- grammar + density baseline
  Q2  solution-type heterogeneity    -- F7: astrometric_params_solved strata
  Q3  RUWE availability by stratum   -- F7: RUWE is undefined for 2-par solutions
  Q4  wide-binary pair self-join     -- F10: the pair-construction grammar, TOP 1000
  Q5  pair count, patch A            -- F10: density for the full-sky extrapolation
  Q6  pair count, patch B            -- F10: second latitude, to bracket the density
  Q7  nss_two_body_orbit row count   -- F8: the companion-census leg
  Q8  nss_acceleration_astro count   -- F8
  Q9  NSS column list (TAP_SCHEMA)   -- F8: the day-one column contract
  Q10 NSS x patch cross-match        -- F8: test-0 grammar end to end

Every call is wrapped: on network failure the ADQL text is still emitted to the
output file, marked UNTESTED, so the queries are banked either way.

Pair-selection physics (El-Badry & Rix 2018 style, stated so a referee can check):
  angular separation      theta[arcsec] = s[AU] * parallax[mas] / 1000
  window 0.2-50 kAU at parallax > 5 mas  ->  theta in [1, 250] arcsec
  parallax consistency    |d(parallax)| < 3 sigma
  pm consistency          |d(mu)| < 0.44 * parallax^1.5 / sqrt(theta) + 3 sigma
                          (the 0.44 coefficient is the projected Keplerian
                           ceiling for a 5 Msun total mass; derived inline below)

Usage:  py calcs/dr4_dryrun_adql.py
Output: data/dr4_dryrun_adql.txt
"""
import csv
import io
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
TIMEOUT = 300
OUT = os.path.join(os.path.dirname(__file__), "..", "data", "dr4_dryrun_adql.txt")

# ---------------------------------------------------------------- patches
# Two cones at contrasting Galactic latitude so the full-sky extrapolation is
# reported as a BRACKET, not a single number (the sample is not isotropic).
PATCH_A = dict(name="A", ra=150.0, dec=10.0, rad=3.0, note="l~228 b~+47, high latitude")
PATCH_B = dict(name="B", ra=100.0, dec=10.0, rad=3.0, note="l~204 b~+12, low latitude")

PLX_MIN = 5.0          # mas; the 9F/Paper-1 sigma_v < 0.03 km/s cut forces this
PLX_OVER_ERR = 20.0
S_MIN_AU = 200.0       # 0.2 kAU
S_MAX_AU = 50000.0     # 50 kAU
JOIN_RAD_DEG = 0.0695  # 250 arcsec = 50 kAU at exactly 5 mas (the widest angle needed)

SKY_DEG2 = 41252.96


def patch_area(p):
    return 3.14159265358979 * p["rad"] ** 2


# ---------------------------------------------------------------- ADQL text
def q_census(p):
    return f"""SELECT COUNT(*) AS n_sources
FROM gaiadr3.gaia_source
WHERE 1 = CONTAINS(POINT('ICRS', ra, dec),
                   CIRCLE('ICRS', {p['ra']}, {p['dec']}, {p['rad']}))
  AND parallax > {PLX_MIN}
  AND parallax_over_error > {PLX_OVER_ERR}"""


def q_solution_types(p):
    """F7. In DR3 astrometric_params_solved is 3 (2-par), 31 (5-par), 95 (6-par).
    DR4's gaia_source is model-CONSOLIDATED (single-star vs binary-star model per
    source), so this axis becomes a first-class stratum axis, not a footnote."""
    return f"""SELECT astrometric_params_solved, COUNT(*) AS n,
       AVG(ruwe) AS ruwe_mean, MIN(ruwe) AS ruwe_min, MAX(ruwe) AS ruwe_max
FROM gaiadr3.gaia_source
WHERE 1 = CONTAINS(POINT('ICRS', ra, dec),
                   CIRCLE('ICRS', {p['ra']}, {p['dec']}, {p['rad']}))
  AND parallax > {PLX_MIN}
GROUP BY astrometric_params_solved"""


def q_ruwe_availability(p):
    """F7. RUWE is a SINGLE-STAR fit statistic; it is null where no 5/6-par
    single-star solution exists.  The DR4 content page does not list RUWE."""
    return f"""SELECT COUNT(*) AS n_total,
       SUM(ruwe) AS ruwe_sum_nonnull_only
FROM gaiadr3.gaia_source
WHERE 1 = CONTAINS(POINT('ICRS', ra, dec),
                   CIRCLE('ICRS', {p['ra']}, {p['dec']}, {p['rad']}))
  AND parallax > {PLX_MIN}
  AND ruwe IS NULL"""


def _pair_body(p, select_clause):
    """Shared FROM/WHERE for the pair self-join.  theta is written out in full
    everywhere because ADQL does not allow select-list aliases in WHERE."""
    sub = f"""(SELECT source_id, ra, dec, parallax, parallax_error,
             pmra, pmdec, pmra_error, pmdec_error,
             phot_g_mean_mag, ruwe, astrometric_params_solved
      FROM gaiadr3.gaia_source
      WHERE 1 = CONTAINS(POINT('ICRS', ra, dec),
                         CIRCLE('ICRS', {p['ra']}, {p['dec']}, {p['rad']}))
        AND parallax > {PLX_MIN}
        AND parallax_over_error > {PLX_OVER_ERR}
        AND astrometric_params_solved >= 31)"""
    theta = ("3600.0 * DISTANCE(POINT('ICRS', a.ra, a.dec), "
             "POINT('ICRS', b.ra, b.dec))")
    s_au = f"({theta} * 2000.0 / (a.parallax + b.parallax))"
    return f"""{select_clause}
FROM {sub} AS a,
     {sub} AS b
WHERE a.source_id < b.source_id
  AND 1 = CONTAINS(POINT('ICRS', a.ra, a.dec),
                   CIRCLE('ICRS', b.ra, b.dec, {JOIN_RAD_DEG}))
  AND {s_au} BETWEEN {S_MIN_AU} AND {S_MAX_AU}
  AND ABS(a.parallax - b.parallax) <
      3.0 * SQRT(a.parallax_error * a.parallax_error
                 + b.parallax_error * b.parallax_error)
  AND SQRT(POWER(a.pmra - b.pmra, 2) + POWER(a.pmdec - b.pmdec, 2)) <
      0.44 * POWER((a.parallax + b.parallax) / 2.0, 1.5) / SQRT({theta})
      + 3.0 * SQRT(a.pmra_error * a.pmra_error + b.pmra_error * b.pmra_error
                   + a.pmdec_error * a.pmdec_error
                   + b.pmdec_error * b.pmdec_error)"""


def q_pairs_list(p):
    theta = ("3600.0 * DISTANCE(POINT('ICRS', a.ra, a.dec), "
             "POINT('ICRS', b.ra, b.dec))")
    sel = f"""SELECT TOP 1000
  a.source_id AS source_id1, b.source_id AS source_id2,
  a.parallax AS parallax1, b.parallax AS parallax2,
  a.pmra AS pmra1, a.pmdec AS pmdec1, b.pmra AS pmra2, b.pmdec AS pmdec2,
  a.pmra_error AS pmra_error1, a.pmdec_error AS pmdec_error1,
  b.pmra_error AS pmra_error2, b.pmdec_error AS pmdec_error2,
  a.phot_g_mean_mag AS g1, b.phot_g_mean_mag AS g2,
  a.ruwe AS ruwe1, b.ruwe AS ruwe2,
  a.astrometric_params_solved AS aps1, b.astrometric_params_solved AS aps2,
  {theta} AS theta_arcsec,
  {theta} * 2000.0 / (a.parallax + b.parallax) AS s_au"""
    return _pair_body(p, sel)


def q_pairs_count(p):
    return _pair_body(p, "SELECT COUNT(*) AS n_pairs")


Q_NSS_ORBIT = """SELECT COUNT(*) AS n_rows FROM gaiadr3.nss_two_body_orbit"""

Q_NSS_ACCEL = """SELECT COUNT(*) AS n_rows FROM gaiadr3.nss_acceleration_astro"""

Q_NSS_COLUMNS = """SELECT TOP 200 column_name, datatype, unit
FROM tap_schema.columns
WHERE table_name = 'gaiadr3.nss_two_body_orbit'
ORDER BY column_name"""


def q_nss_crossmatch(p):
    """F8 / test 0.  Fraction of the wide-binary parent population that already
    carries a published NSS solution.  In DR4 this is the companion census."""
    return f"""SELECT COUNT(*) AS n_with_nss
FROM gaiadr3.gaia_source AS g,
     gaiadr3.nss_two_body_orbit AS n
WHERE g.source_id = n.source_id
  AND 1 = CONTAINS(POINT('ICRS', g.ra, g.dec),
                   CIRCLE('ICRS', {p['ra']}, {p['dec']}, {p['rad']}))
  AND g.parallax > {PLX_MIN}"""


# ---------------------------------------------------------------- TAP driver
def tap_query(adql, label):
    """Return (ok, rows, elapsed_s, message).  rows = list of dicts."""
    data = urllib.parse.urlencode({
        "REQUEST": "doQuery",
        "LANG": "ADQL",
        "FORMAT": "csv",
        "QUERY": adql,
    }).encode("ascii")
    req = urllib.request.Request(
        TAP, data=data,
        headers={"User-Agent": "PhysicsResearch DR4 dry-run (research, sync, small TOP)"})
    t0 = time.time()
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
            raw = r.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        body = ""
        try:
            body = e.read().decode("utf-8", errors="replace")[:600]
        except Exception:
            pass
        return False, [], time.time() - t0, f"HTTP {e.code}: {body}"
    except Exception as e:
        return False, [], time.time() - t0, f"{type(e).__name__}: {e}"
    dt = time.time() - t0
    if raw.lstrip().startswith("<"):
        return False, [], dt, "VOTable/XML error payload: " + raw[:600]
    rows = list(csv.DictReader(io.StringIO(raw)))
    return True, rows, dt, "OK"


def fnum(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------- main
def main():
    t_start = time.time()
    probes = [
        ("Q1  census patch A", q_census(PATCH_A)),
        ("Q1b census patch B", q_census(PATCH_B)),
        ("Q2  solution types patch A (F7)", q_solution_types(PATCH_A)),
        ("Q3  RUWE-null count patch A (F7)", q_ruwe_availability(PATCH_A)),
        ("Q4  pair self-join listing patch A (F10)", q_pairs_list(PATCH_A)),
        ("Q5  pair count patch A (F10)", q_pairs_count(PATCH_A)),
        ("Q6  pair count patch B (F10)", q_pairs_count(PATCH_B)),
        ("Q7  nss_two_body_orbit rows (F8)", Q_NSS_ORBIT),
        ("Q8  nss_acceleration_astro rows (F8)", Q_NSS_ACCEL),
        ("Q9  nss_two_body_orbit columns (F8)", Q_NSS_COLUMNS),
        ("Q10 NSS x patch A cross-match (F8)", q_nss_crossmatch(PATCH_A)),
    ]

    results = []
    for label, adql in probes:
        print(f"[{label}] running ...", flush=True)
        ok, rows, dt, msg = tap_query(adql, label)
        print(f"    {'OK ' if ok else 'FAIL'}  {len(rows)} rows  {dt:.1f} s  {msg[:120]}",
              flush=True)
        results.append(dict(label=label, adql=adql, ok=ok, rows=rows, dt=dt, msg=msg))

    by = {r["label"].split()[0]: r for r in results}

    # ---- full-sky pair-volume extrapolation ---------------------------------
    extrap = []
    dens = {}
    for key, p in (("Q5", PATCH_A), ("Q6", PATCH_B)):
        r = by.get(key)
        if r and r["ok"] and r["rows"]:
            n = fnum(list(r["rows"][0].values())[0])
            if n is not None:
                d = n / patch_area(p)
                dens[p["name"]] = (n, d)
                extrap.append(f"  patch {p['name']} ({p['note']}): "
                              f"{int(n)} pairs / {patch_area(p):.1f} deg2 = {d:.2f} pairs/deg2 "
                              f"-> full sky {d * SKY_DEG2:,.0f}")
    if len(dens) == 2:
        lo = min(v[1] for v in dens.values()) * SKY_DEG2
        hi = max(v[1] for v in dens.values()) * SKY_DEG2
        extrap.append(f"  FULL-SKY BRACKET (crude, isotropy assumed within each "
                      f"latitude class): {lo:,.0f} - {hi:,.0f} pairs at parallax > "
                      f"{PLX_MIN:.0f} mas, s in [0.2, 50] kAU, DR3 astrometry")
        extrap.append("  RISK AXIS (ship it with the number, per the standing rule):")
        extrap.append("    - two cones only; the true distribution is strongly "
                      "latitude- and crowding-dependent, so this is a bracket of "
                      "two samples, not a confidence interval;")
        extrap.append("    - no chance-alignment subtraction (the pm+parallax "
                      "consistency cuts are the only contamination control here; "
                      "the 4I R_chance ladder is NOT applied);")
        extrap.append("    - edge truncation: pairs straddling the cone boundary "
                      "are lost (perimeter annulus ~4.6% of area, partial loss "
                      "within it, so the deficit is <~2%);")
        extrap.append("    - DR4 loosens the parallax floor at fixed sigma_v by the "
                      "DR4 pm-improvement factor f: the volume scales ~f^3 and the "
                      "yield ~f^3 only if the magnitude and main-sequence cuts do "
                      "not bind. F13 says compute that, do not assume '~10x'.")

    # ---- write the bank ------------------------------------------------------
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        w = f.write
        w("DR4 DAY-ONE ADQL DRY-RUN -- banked queries and DR3 regression targets\n")
        w("=" * 78 + "\n")
        w(f"generated : {time.strftime('%Y-%m-%d %H:%M:%S')} local\n")
        w(f"endpoint  : {TAP} (sync, CSV)\n")
        w(f"script    : calcs/dr4_dryrun_adql.py\n")
        w("purpose   : red-team findings F7 (solution heterogeneity),\n")
        w("            F8 (NSS companion census = test 0), F10 (data access).\n")
        w("status    : DRY RUN against Gaia DR3. Nothing here is pre-registered.\n")
        w("            DR4 table/column names MUST be re-verified against the\n")
        w("            published DR4 data model before any of this is fired.\n")
        w("\n")
        w("SELECTION CONSTANTS\n")
        w("-" * 78 + "\n")
        w(f"  parallax > {PLX_MIN} mas (200 pc), parallax_over_error > {PLX_OVER_ERR}\n")
        w(f"  projected separation s in [{S_MIN_AU:.0f}, {S_MAX_AU:.0f}] AU\n")
        w(f"  join radius {JOIN_RAD_DEG} deg = {JOIN_RAD_DEG*3600:.0f} arcsec "
          f"= 50 kAU at exactly {PLX_MIN:.0f} mas\n")
        w("  pm-consistency ceiling coefficient 0.44 derived as:\n")
        w("    v_max = 29.8 km/s * sqrt(M/s_AU) with M = 5 Msun -> 66.6/sqrt(s_AU)\n")
        w("    dmu[mas/yr] = v * plx / 4.74 ; s_AU = theta_arcsec * 1000 / plx\n")
        w("    => dmu_max = 0.444 * plx^1.5 / sqrt(theta_arcsec)\n")
        w("    (this reproduces the El-Badry & Rix 2018 coefficient independently)\n")
        w("\n")
        w("PATCHES\n")
        w("-" * 78 + "\n")
        for p in (PATCH_A, PATCH_B):
            w(f"  {p['name']}: ra={p['ra']} dec={p['dec']} radius={p['rad']} deg "
              f"(area {patch_area(p):.1f} deg2) -- {p['note']}\n")
        w("\n")
        w("FULL-SKY PAIR-VOLUME EXTRAPOLATION\n")
        w("-" * 78 + "\n")
        if extrap:
            for line in extrap:
                w(line + "\n")
        else:
            w("  NOT AVAILABLE -- the pair-count probes did not return "
              "(see per-query status below).\n")
        w("\n")
        w("PER-QUERY RESULTS\n")
        w("=" * 78 + "\n")
        for r in results:
            w("\n")
            w(f"### {r['label']}\n")
            w(f"status   : {'OK' if r['ok'] else 'UNTESTED (call failed)'}\n")
            w(f"elapsed  : {r['dt']:.1f} s\n")
            w(f"message  : {r['msg'][:400]}\n")
            w(f"rows     : {len(r['rows'])}\n")
            w("ADQL:\n")
            for line in r["adql"].splitlines():
                w("    " + line + "\n")
            if r["ok"] and r["rows"]:
                # the TAP_SCHEMA probe IS the column contract -- dump it whole
                cap = 250 if "columns" in r["label"] else 40
                w(f"RESULT (first {cap} rows):\n")
                keys = list(r["rows"][0].keys())
                w("    " + " | ".join(keys) + "\n")
                for row in r["rows"][:cap]:
                    w("    " + " | ".join(str(row.get(k, "")) for k in keys) + "\n")
                if len(r["rows"]) > cap:
                    w(f"    ... {len(r['rows']) - cap} more rows suppressed\n")
            elif not r["ok"]:
                w("RESULT   : none -- ADQL banked UNTESTED.\n")
        w("\n")
        w("=" * 78 + "\n")
        w("DAY-ONE DR4 FORMS (UNTESTED -- schema not published at write time)\n")
        w("=" * 78 + "\n")
        w("Replace gaiadr3 -> gaiadr4 throughout.  Known DR4 differences that\n")
        w("break these queries if not handled (ESA DR4 content page):\n")
        w("  1. gaia_source in DR4 is the ~2.0e9 HIGH-QUALITY SUBSET of ~2.8e9\n")
        w("     processed sources.  The complete astrometry lives in\n")
        w("     all_source_astrometry.  Decide which is the parent sample and\n")
        w("     pre-commit it: the subset is quality-selected, and that selection\n")
        w("     is correlated with binarity.\n")
        w("  2. gaia_source parameters are CONSOLIDATED across processing modules\n")
        w("     -- a binary-star model may have been applied instead of the\n")
        w("     single-star model.  astrometric_params_solved is therefore NOT the\n")
        w("     full stratum axis in DR4; the applied-model flag is (see\n")
        w("     all_source_flags).  This is finding F7.\n")
        w("  3. RUWE is not listed on the DR4 content page.  The 9J dose curve must\n")
        w("     be RE-DERIVED on whatever DR4 publishes, never transferred.\n")
        w("  4. NSS tables in DR4 (10 of them, ESA content page): nss_two_body_orbit,\n")
        w("     nss_acceleration_astro, nss_resolved_pair, nss_multiple_orbits,\n")
        w("     nss_masses, nss_non_linear_spectro, nss_multiplicity, nss_vim_fl,\n")
        w("     optical_pair (+1).  Test 0 rides on these.\n")
        w("  5. epoch_astrometry is ~62 TB.  Any per-epoch work must be a\n")
        w("     source_id-list-driven server-side join, never a bulk pull.\n")
        w("\n")
        w("FALLBACK ORDER (pre-committed here, per F10):\n")
        w("  (a) our own ADQL above, run against the DR4 archive;\n")
        w("  (b) a community DR4 wide-binary catalogue IF one exists AND it\n")
        w("      reproduces our cut ladder within tolerance;\n")
        w("  (c) bulk mirror download (AIP gaia.aip.de, ARI/Heidelberg, VizieR/CDS)\n")
        w("      only if (a) is rate-limited past the first week.\n")
        w("\n")
        w(f"TOTAL WALL CLOCK: {time.time() - t_start:.1f} s\n")

    # ---- console summary -----------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    n_ok = sum(1 for r in results if r["ok"])
    for r in results:
        head = ""
        if r["ok"] and r["rows"]:
            first = r["rows"][0]
            head = " | ".join(f"{k}={v}" for k, v in list(first.items())[:3])
        print(f"  {'OK  ' if r['ok'] else 'FAIL'} {r['label']:<42} "
              f"rows={len(r['rows']):<5} {r['dt']:6.1f}s  {head[:60]}")
    print(f"\n  {n_ok}/{len(results)} probes returned")
    for line in extrap:
        print(line)
    print(f"\n  wall clock: {time.time() - t_start:.1f} s")
    print(f"  written   : {os.path.normpath(OUT)}")
    return 0 if n_ok else 1


if __name__ == "__main__":
    sys.exit(main())
