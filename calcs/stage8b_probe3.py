"""8B probe part 3 — the Aperture route with correct object names.
Tiny queries.  Output: data/stage8b_probe3.txt (no credentials)."""
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

L = []
def P(s):
    print(s, flush=True)
    L.append(s)

def query(sql, tag, maxrows=40):
    url = BASE + '?action=doQuery&SQL=' + urllib.parse.quote(sql)
    for attempt in (1, 2):
        try:
            t0 = time.time()
            with opener.open(url, timeout=600) as r:
                txt = r.read().decode('utf-8', errors='replace')
            dt = time.time() - t0
            lines = [l for l in txt.splitlines() if l.strip()]
            P(f"[{tag}] {dt:.1f}s, {len(lines)} lines:")
            for l in lines[:maxrows]:
                P("   " + l)
            time.sleep(2.5)
            return lines
        except Exception as e:
            P(f"[{tag}] attempt {attempt} FAILED: {e}")
            time.sleep(10)
    return []

P("8B PROBE 3 — RefL0100N1504 Aperture route")
P("")
query("SELECT TOP 12 GalaxyID, ApertureSize, Mass_Star, Mass_Gas, "
      "Mass_DM, Mass_BH FROM RefL0100N1504_Aperture "
      "ORDER BY GalaxyID, ApertureSize", "Q1 Aperture TOP12")
query("SELECT DISTINCT ApertureSize FROM RefL0100N1504_Aperture "
      "WHERE GalaxyID = (SELECT MIN(GalaxyID) FROM "
      "RefL0100N1504_Aperture) ORDER BY ApertureSize",
      "Q2 aperture sizes")
query("SELECT TOP 5 GalaxyID, SnapNum, Redshift, SubGroupNumber, "
      "Spurious, MassType_Star, StarFormationRate, "
      "HalfMassRad_Star FROM RefL0100N1504_Subhalo "
      "WHERE SnapNum = 28", "Q3 Subhalo TOP5")
query("SELECT COUNT(*) AS n FROM RefL0100N1504_Subhalo AS sub, "
      "RefL0100N1504_Aperture AS ap "
      "WHERE ap.GalaxyID = sub.GalaxyID AND sub.SnapNum = 28 "
      "AND sub.SubGroupNumber = 0 AND sub.Spurious = 0 "
      "AND ap.ApertureSize = 30 "
      "AND ap.Mass_Star BETWEEN 1e9 AND 1e11",
      "Q4 selection count (centrals, 1e9-1e11 Msun)")
query("SELECT TOP 5 GalaxyID, KappaCoRot, DiscToTotal "
      "FROM RefL0100N1504_MorphoKinem", "Q5 MorphoKinem probe")

with open('data/stage8b_probe3.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8b_probe3.txt")
