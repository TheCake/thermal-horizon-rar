"""8B fetch — the EAGLE aperture pull (pre-reg c03604f; polite: 6
chunked queries, ~4 s pauses, cached one-time to data/eagle/,
gitignored).  Credentials from private/ (never committed/printed).
Output log: data/stage8b_fetch.txt."""
import os
import time
import urllib.request
import urllib.parse

BASE = 'http://virgodb.dur.ac.uk:8080/Eagle'
CRED = open('private/eagle_creds.txt').read().split()

pm = urllib.request.HTTPPasswordMgrWithDefaultRealm()
pm.add_password(None, BASE, CRED[0], CRED[1])
opener = urllib.request.build_opener(
    urllib.request.HTTPBasicAuthHandler(pm))
opener.addheaders = [('User-Agent',
                      'PhysicsResearch-RAR-ladder (polite scripted '
                      'access; hfilip11@gmail.com)')]

os.makedirs('data/eagle', exist_ok=True)
L = []
def P(s):
    print(s, flush=True)
    L.append(s)

def query_rows(sql, tag):
    url = BASE + '?action=doQuery&SQL=' + urllib.parse.quote(sql)
    for attempt in (1, 2):
        try:
            t0 = time.time()
            with opener.open(url, timeout=900) as r:
                txt = r.read().decode('utf-8', errors='replace')
            lines = txt.splitlines()
            meta = [l for l in lines if l.startswith('#')]
            ok = any(l.startswith('#OK') for l in meta)
            data = [l for l in lines if l.strip() and
                    not l.startswith('#')]
            P(f"[{tag}] {time.time()-t0:.1f}s, ok={ok}, "
              f"{len(data)-1 if data else 0} rows")
            if not ok:
                for l in meta[-3:]:
                    P("   " + l)
            time.sleep(4.0)
            return data if ok else []
        except Exception as e:
            P(f"[{tag}] attempt {attempt} FAILED: {e}")
            time.sleep(15)
    return []

def main_pull(sim, out, bins):
    header, rows = None, []
    for lo, hi in bins:
        sql = (f"SELECT ap.GalaxyID, ap.ApertureSize, ap.Mass_Star, "
               f"ap.Mass_Gas, ap.Mass_DM, ap.Mass_BH "
               f"FROM {sim}_Aperture AS ap, {sim}_Subhalo AS sub, "
               f"{sim}_Aperture AS sel "
               f"WHERE ap.GalaxyID = sub.GalaxyID "
               f"AND sel.GalaxyID = sub.GalaxyID "
               f"AND sub.SnapNum = 28 AND sub.SubGroupNumber = 0 "
               f"AND sub.Spurious = 0 AND sel.ApertureSize = 30 "
               f"AND sel.Mass_Star >= {lo} AND sel.Mass_Star < {hi} "
               f"ORDER BY ap.GalaxyID, ap.ApertureSize")
        data = query_rows(sql, f"{sim} M*[{lo:.0e},{hi:.0e})")
        if data:
            if header is None:
                header = data[0]
            rows += data[1:]
    if header:
        with open(out, 'w') as f:
            f.write(header + "\n" + "\n".join(rows) + "\n")
        P(f"  -> {out}: {len(rows)} rows")

# main box, 4 mass chunks (edges match the probe count window)
main_pull('RefL0100N1504', 'data/eagle/ref100_apertures.csv',
          [(1e9, 3e9), (3e9, 1e10), (1e10, 3e10), (3e10, 1.00001e11)])

# morphokinem (disk-subset leg)
mk = query_rows(
    "SELECT mk.GalaxyID, mk.KappaCoRot, mk.DiscToTotal "
    "FROM RefL0100N1504_MorphoKinem AS mk, "
    "RefL0100N1504_Subhalo AS sub, RefL0100N1504_Aperture AS sel "
    "WHERE mk.GalaxyID = sub.GalaxyID AND sel.GalaxyID = sub.GalaxyID "
    "AND sub.SnapNum = 28 AND sub.SubGroupNumber = 0 "
    "AND sub.Spurious = 0 AND sel.ApertureSize = 30 "
    "AND sel.Mass_Star BETWEEN 1e9 AND 1e11", "MorphoKinem")
if mk:
    with open('data/eagle/ref100_morphokinem.csv', 'w') as f:
        f.write("\n".join(mk) + "\n")
    P(f"  -> data/eagle/ref100_morphokinem.csv: {len(mk)-1} rows")

# the Recal 25 Mpc high-res cross-check box (single query)
main_pull('RecalL0025N0752', 'data/eagle/recal25_apertures.csv',
          [(1e9, 1.00001e11)])

with open('data/stage8b_fetch.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8b_fetch.txt")
