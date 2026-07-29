"""8B probe part 2 — discover the EAGLE database/table naming.
Tiny queries only.  Output: data/stage8b_probe2.txt (no credentials)."""
import re
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

def get(url, tag):
    try:
        t0 = time.time()
        with opener.open(url, timeout=300) as r:
            txt = r.read().decode('utf-8', errors='replace')
        P(f"[{tag}] {time.time()-t0:.1f}s, {len(txt)} chars")
        time.sleep(2.5)
        return txt
    except Exception as e:
        P(f"[{tag}] FAILED: {e}")
        time.sleep(5)
        return ''

def query(sql, tag, maxrows=60):
    txt = get(BASE + '?action=doQuery&SQL=' + urllib.parse.quote(sql),
              tag)
    for l in [l for l in txt.splitlines() if l.strip()][:maxrows]:
        P("   " + l)

# 1) the query page itself — the left-hand database list
page = get(BASE, 'page')
names = sorted(set(re.findall(
    r'([A-Za-z0-9_]*(?:[Ee]agle|Ref[LN]|Recal|AGNdT)[A-Za-z0-9_]*)',
    page)))
P("name-like tokens on the page: " + ", ".join(names[:60]))
P("")

# 2) system catalog of databases (cheap)
query("SELECT name FROM sys.databases", "sys.databases")
P("")

# 3) tables in the default database
query("SELECT TOP 60 TABLE_CATALOG, TABLE_NAME "
      "FROM information_schema.tables ORDER BY TABLE_NAME",
      "info_schema default")

with open('data/stage8b_probe2.txt', 'w') as f:
    f.write("\n".join(L) + "\n")
print("\nsaved: data/stage8b_probe2.txt")
