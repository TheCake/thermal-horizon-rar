"""STAGE 8B probe — EAGLE public database schema + selection counts.

Instrument-building only (no measurement, no bars): confirm the
Aperture-table route for the simulation ladder — table names, aperture
sizes, column names, and the size of the z = 0 central-galaxy selection
— with TOP-N and COUNT queries (server load: seconds).

Credentials: private/eagle_creds.txt (gitignored; NEVER committed or
printed).  Politeness: single connection, 2.5 s sleep between queries,
<= 2 tries per query, generous timeout, identifying User-Agent.

Output: data/stage8b_probe.txt (no credentials inside).
"""
import time
import urllib.request
import urllib.parse

BASE = 'http://virgodb.dur.ac.uk:8080/Eagle'
CRED = open('private/eagle_creds.txt').read().split()
USER, PW = CRED[0], CRED[1]

pm = urllib.request.HTTPPasswordMgrWithDefaultRealm()
pm.add_password(None, BASE, USER, PW)
opener = urllib.request.build_opener(
    urllib.request.HTTPBasicAuthHandler(pm))
opener.addheaders = [('User-Agent',
                      'PhysicsResearch-RAR-ladder (polite scripted '
                      'access; hfilip11@gmail.com)')]

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

def query(sql, tag, maxrows=40):
    url = BASE + '?action=doQuery&SQL=' + urllib.parse.quote(sql)
    for attempt in (1, 2):
        try:
            t0 = time.time()
            with opener.open(url, timeout=300) as r:
                txt = r.read().decode('utf-8', errors='replace')
            dt = time.time() - t0
            lines = [l for l in txt.splitlines() if l.strip()]
            P(f"[{tag}] {dt:.1f}s, {len(lines)} lines "
              f"(showing <= {maxrows}):")
            for l in lines[:maxrows]:
                P("   " + l)
            time.sleep(2.5)
            return lines
        except Exception as e:
            P(f"[{tag}] attempt {attempt} FAILED: {e}")
            time.sleep(10)
    return []

P("8B PROBE — EAGLE database, Aperture route feasibility")
P("")

# Q1: aperture table shape
query("SELECT TOP 5 GalaxyID, ApertureSize, Mass_Star, Mass_Gas, "
      "Mass_DM, Mass_BH FROM RefL0100N1504..Aperture "
      "ORDER BY GalaxyID, ApertureSize", "Q1 Aperture TOP5")

# Q2: which aperture sizes exist
query("SELECT DISTINCT ApertureSize FROM RefL0100N1504..Aperture "
      "WHERE GalaxyID = (SELECT MIN(GalaxyID) FROM "
      "RefL0100N1504..Aperture)", "Q2 aperture sizes")

# Q3: subhalo columns we need
query("SELECT TOP 5 GalaxyID, SnapNum, Redshift, SubGroupNumber, "
      "Spurious, MassType_Star, StarFormationRate, "
      "HalfMassRad_Star FROM RefL0100N1504..SubHalo "
      "WHERE SnapNum = 28", "Q3 SubHalo TOP5")

# Q4: selection count at z = 0 — centrals, SPARC-like stellar mass
query("SELECT COUNT(*) AS n FROM RefL0100N1504..SubHalo AS sub, "
      "RefL0100N1504..Aperture AS ap "
      "WHERE ap.GalaxyID = sub.GalaxyID AND sub.SnapNum = 28 "
      "AND sub.SubGroupNumber = 0 AND sub.Spurious = 0 "
      "AND ap.ApertureSize = 30 "
      "AND ap.Mass_Star BETWEEN 1e9 AND 1e11",
      "Q4 selection count")

# Q5: full-pull row estimate (that count x ~10 apertures)
query("SELECT COUNT(*) AS n FROM RefL0100N1504..Aperture AS ap, "
      "RefL0100N1504..SubHalo AS sub "
      "WHERE ap.GalaxyID = sub.GalaxyID AND sub.SnapNum = 28 "
      "AND sub.SubGroupNumber = 0 AND sub.Spurious = 0 "
      "AND ap.GalaxyID IN (SELECT GalaxyID FROM "
      "RefL0100N1504..Aperture WHERE ApertureSize = 30 "
      "AND Mass_Star BETWEEN 1e9 AND 1e11)",
      "Q5 pull-size count")

with open('data/stage8b_probe.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8b_probe.txt")
