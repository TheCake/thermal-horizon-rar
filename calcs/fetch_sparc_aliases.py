"""
Fetch SIMBAD cross-identifications for SPARC galaxies that failed the
10I name-based UNGC match (R30 condition 7 / successor: "fix the matcher
with a real alias table before banking any census").

Queries the SIMBAD TAP service (ident table self-join) for every SPARC
galaxy with D <= 12 Mpc that the 10I canon-matcher left unmatched, plus
the exemplar UGC05721 (= NGC3274, R30-verified).  Writes
data/sparc_simbad_aliases.csv with rows (sparc_name, alias) -- one row
per catalogue identifier SIMBAD returns.  Names that do not resolve are
logged in the CSV header comment (honest-miss record).  Idempotent:
skips the fetch if the CSV already exists (delete to refetch).

The list of names below is NOT hand-typed: it is the verbatim
"unmatched SPARC with D<=12 Mpc" line from data/stage10i_isoladder.txt
(run 1, committed bf7a3dd).
"""
import csv, os, re, time
import urllib.parse
import urllib.request

OUT = 'data/sparc_simbad_aliases.csv'
TAP = 'https://simbad.cds.unistra.fr/simbad/sim-tap/sync'

UNMATCHED = ['D564-8', 'D631-7', 'ESO444-G084', 'KK98-251', 'NGC1003',
             'UGC00891', 'UGC02259', 'UGC04278', 'UGC04325', 'UGC05414',
             'UGC05721', 'UGC05764', 'UGC05829', 'UGC05918', 'UGC05986',
             'UGC06446', 'UGC07151', 'UGC07232', 'UGC07323', 'UGC07399',
             'UGC07524', 'UGC07603', 'UGC08286', 'UGC08490', 'UGC08550',
             'UGC08837', 'UGC12632', 'UGCA444']


def candidates(nm):
    """SIMBAD query-string candidates for a SPARC designation."""
    c = [nm]
    m = re.match(r'^(NGC|UGC|IC|DDO|ESO|PGC|UGCA)[\s\-_]?0*(\d.*)$', nm,
                 re.I)
    if m:
        c += ['%s %s' % (m.group(1).upper(), m.group(2)),
              '%s%s' % (m.group(1).upper(), m.group(2))]
    m = re.match(r'^KK98[\s\-_]?0*(\d+)$', nm, re.I)
    if m:
        c += ['[KK98] %s' % m.group(1)]
    m = re.match(r'^ESO(\d+)-G?0*(\d+)$', nm, re.I)
    if m:
        c += ['ESO %s-%s' % (m.group(1), m.group(2)),
              'ESO %s-G%s' % (m.group(1), int(m.group(2)))]
    m = re.match(r'^D(\d+)-(\d+)$', nm)
    if m:
        c += ['[SPS97] D%s-%s' % (m.group(1), m.group(2)),
              'D%s-%s' % (m.group(1), m.group(2))]
    seen, out = set(), []
    for x in c:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


def tap_ids(qname):
    """All SIMBAD identifiers of the object named qname, or []."""
    q = ("SELECT id2.id FROM ident AS id1 JOIN ident AS id2 "
         "ON id1.oidref = id2.oidref WHERE id1.id = '%s'"
         % qname.replace("'", "''"))
    try:
        body = urllib.parse.urlencode(dict(REQUEST='doQuery', LANG='ADQL',
                                           FORMAT='csv', QUERY=q)).encode()
        with urllib.request.urlopen(
                urllib.request.Request(TAP, data=body), timeout=30) as r:
            text = r.read().decode('utf-8', errors='replace')
        lines = text.strip().splitlines()
        return [l.strip().strip('"') for l in lines[1:] if l.strip()]
    except Exception as e:
        print('  TAP error for %r: %s' % (qname, e))
        return None            # None = network error, [] = no object


def main():
    if os.path.exists(OUT):
        print('%s exists -- skipping fetch (delete to refetch)' % OUT)
        return
    rows, misses, errors = [], [], []
    for nm in UNMATCHED:
        got = None
        for cand in candidates(nm):
            ids = tap_ids(cand)
            time.sleep(0.4)
            if ids is None:
                errors.append(nm)
                break
            if ids:
                got = ids
                break
        if got:
            for al in got:
                rows.append((nm, al))
            print('%s: %d identifiers (via %r)' % (nm, len(got), cand))
        elif nm not in errors:
            misses.append(nm)
            print('%s: NOT RESOLVED' % nm)
    with open(OUT, 'w', newline='') as f:
        f.write('# SIMBAD TAP cross-IDs fetched %s\n'
                % time.strftime('%Y-%m-%d'))
        f.write('# unresolved: %s\n' % (misses if misses else 'none'))
        f.write('# tap-errors: %s\n' % (errors if errors else 'none'))
        w = csv.writer(f)
        w.writerow(['sparc_name', 'alias'])
        w.writerows(rows)
    print('wrote %s: %d alias rows, %d unresolved, %d errors'
          % (OUT, len(rows), len(misses), len(errors)))


if __name__ == '__main__':
    main()
