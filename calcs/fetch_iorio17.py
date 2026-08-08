"""
Fetch the Iorio et al. 2017 (MNRAS 466, 4159) LITTLE THINGS in 3D final
rotation/circular-velocity curves -- the PRIMARY data release.

Source: the paper's arXiv comments field ("The final rotation curves can be
downloaded from [author's website]") -> Filippo Fraternali's downloads page
(https://www.filippofraternali.com/downloads/index.html), item "Rotation
curves of dwarf galaxies": complete dataset of circular velocities for 17
LITTLE THINGS galaxies from Iorio et al. 2017.

There is NO VizieR catalog for this paper (TAP probe 2026-08-08: zero rows
for J/MNRAS/466/4159) -- this zip is the canonical machine-readable release.

Output: data/iorio17_finalrot.zip -> unpacked to data/iorio17/
Prints SHA256 of the zip + a file listing for the manifest.
"""
import hashlib
import os
import urllib.request
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

URL = "https://www.dropbox.com/s/t4j8dacmnwgj0yb/finalrot.zip?dl=1"
ZIP = "data/iorio17_finalrot.zip"
DEST = "data/iorio17"

req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
with urllib.request.urlopen(req, timeout=120) as r, open(ZIP, "wb") as f:
    f.write(r.read())

h = hashlib.sha256(open(ZIP, "rb").read()).hexdigest()
print(f"fetched {ZIP}  bytes={os.path.getsize(ZIP)}  sha256={h}")

os.makedirs(DEST, exist_ok=True)
with zipfile.ZipFile(ZIP) as z:
    z.extractall(DEST)
    names = z.namelist()
print(f"unpacked {len(names)} entries to {DEST}/")
for n in sorted(names):
    print("  ", n)
