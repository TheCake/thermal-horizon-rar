"""Fetch the Oh et al. 2015 LITTLE THINGS mass-model tables (VizieR
J/AJ/149/180) into data/littlethings/ (gitignored; manifest committed
at data/littlethings_manifest.sha256).

Files: table1.dat (sample properties incl. D, i, e_i), table2.dat
(masses incl. Mgas, MstarSED + no-3.6um flags), rotdmbar.dat (total
rotation curves, scaled; Data/Model rows), rotdm.dat (DM-only curves,
scaled). The stage-9R loader un-scales via each file's own per-row
(R0.3, V0.3) and reconstructs the baryonic curve by matched-radius
subtraction V_bar^2 = V_tot^2 - V_DM^2 inside Oh+'s own mass models.
"""
import gzip
import hashlib
import os
import urllib.request

BASE = "https://cdsarc.cds.unistra.fr/ftp/J/AJ/149/180/"
FILES = ["table1.dat", "table2.dat", "rotdmbar.dat", "rotdm.dat"]

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DST = os.path.join(ROOT, "data", "littlethings")
os.makedirs(DST, exist_ok=True)

for f in FILES:
    gz = urllib.request.urlopen(BASE + f + ".gz", timeout=120).read()
    raw = gzip.decompress(gz)
    with open(os.path.join(DST, f), "wb") as out:
        out.write(raw)
    h = hashlib.sha256(raw).hexdigest()
    print("%s  %d bytes  sha256 %s" % (f, len(raw), h))
print("done; verify against data/littlethings_manifest.sha256")
