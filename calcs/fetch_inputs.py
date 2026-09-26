"""Fetch the program's core public inputs and verify the provenance manifest.

The two remote datasets every headline number ultimately loads:

1. The Gaia EDR3 wide-binary catalog -- El-Badry, Rix & Heintz (2021),
   Zenodo record 4435257, file all_columns_catalog.fits.gz
   (1,419,826,533 bytes).  Saved under the program's local name for the
   same bytes: data/edr3_binaries.fits.gz.  Identity verified both ways:
   the local file's md5 equals Zenodo's published md5
   (9644d53f056758fb2ab1d3d04a45e25d) and its sha256 is pinned in
   data/MANIFEST.sha256.

2. The SPARC rotation-curve release -- Lelli, McGaugh & Schombert (2016),
   Zenodo record 16284118:
     SPARC_Lelli2016c.mrt  ->  data/sparc/SPARC_Lelli2016c.mrt
     Rotmod_LTG.zip        ->  data/sparc/Rotmod_LTG.zip, extracted to
                               data/sparc/rotmod/  (175 *_rotmod.dat files)

Every entry of data/MANIFEST.sha256 (written by the pre-registered
calcs/stage9h_manifest.py) is then verified, including the rotmod
aggregate hash recomputed with stage 9H's own recipe (sorted paths,
sha256 of each file, outer sha256 over the hex digests).  The committed
inputs in the manifest (chae2021_table3.csv and the derived .npy/.npz
tables) are verified on the way.

The program's third dataset is NOT fetched here: the El-Badry & Rix
(2018) DR2 pair list (VizieR J/MNRAS/480/4884, via TAP) feeds only the
exploratory calcs/gaia_widebinary_orientation.py, none of the papers'
numbers.  Stages that need further, smaller datasets name their sources
in their own docstrings (calcs/fetch_*.py, calcs/stage4e_lensing_rar.py).

Usage (from anywhere; the script changes to the repository root):
  py calcs/fetch_inputs.py             fetch what is missing, verify all
  py calcs/fetch_inputs.py --verify    verify only, no network
  py calcs/fetch_inputs.py --selftest  download the two small SPARC files
                                       into a temporary directory and
                                       check their hashes (tests the
                                       network path; never touches data/)

Standard library only.  Exit code 0 iff every check passes.
"""
import glob
import hashlib
import io
import os
import re
import shutil
import sys
import tempfile
import urllib.request
import zipfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

UA = {"User-Agent": "thermal-horizon-rar fetch_inputs (stdlib urllib)"}
MANIFEST = "data/MANIFEST.sha256"

EDR3_URL = ("https://zenodo.org/api/records/4435257/files/"
            "all_columns_catalog.fits.gz/content")
EDR3_DST = "data/edr3_binaries.fits.gz"
EDR3_SHA = "6d67d12bf15f6c579e4072aecbb995e8ce13f057ba1dc83115e7046d9c13bacf"
EDR3_SIZE = 1419826533

MRT_URL = ("https://zenodo.org/api/records/16284118/files/"
           "SPARC_Lelli2016c.mrt/content")
MRT_DST = "data/sparc/SPARC_Lelli2016c.mrt"
MRT_SHA = "5aa0501f6b0d881fa579030e315e7b5b6ef561a5bd3a07472f9929c7e5728243"
MRT_SIZE = 28259

ZIP_URL = ("https://zenodo.org/api/records/16284118/files/"
           "Rotmod_LTG.zip/content")
ZIP_DST = "data/sparc/Rotmod_LTG.zip"
ZIP_SHA = "0a80cc90714828cc28b7dd57923576714d209f2490328c087c4a4ad607faf588"
ZIP_SIZE = 110737

ROTMOD_DIR = "data/sparc/rotmod"
ROTMOD_GLOB = "data/sparc/rotmod/**/*_rotmod.dat"
ROTMOD_COUNT = 175


def sha256_file(path, bufsz=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(bufsz)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def download(url, dst, expect_sha, expect_size, quiet=False):
    """Stream url -> dst, hashing on the way; refuse to keep a mismatch."""
    os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
    part = dst + ".part"
    h = hashlib.sha256()
    done = 0
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=300) as r, open(part, "wb") as f:
        while True:
            b = r.read(1 << 22)
            if not b:
                break
            h.update(b)
            f.write(b)
            done += len(b)
            if not quiet and expect_size >= (1 << 28) and done % (1 << 27) < (1 << 22):
                print(f"    ... {done/1e9:.2f} / {expect_size/1e9:.2f} GB",
                      flush=True)
    if done != expect_size or h.hexdigest() != expect_sha:
        os.remove(part)
        raise RuntimeError(
            f"download mismatch for {url}\n"
            f"  got  size {done}  sha256 {h.hexdigest()}\n"
            f"  want size {expect_size}  sha256 {expect_sha}")
    os.replace(part, dst)
    print(f"  fetched {dst}  ({done} bytes, sha256 OK)")


def extract_rotmod(zip_path):
    os.makedirs(ROTMOD_DIR, exist_ok=True)
    n = 0
    with zipfile.ZipFile(zip_path) as z:
        for member in z.namelist():
            base = os.path.basename(member)
            if not base.endswith("_rotmod.dat"):
                continue
            with z.open(member) as src, \
                    open(os.path.join(ROTMOD_DIR, base), "wb") as out:
                shutil.copyfileobj(src, out)
            n += 1
    print(f"  extracted {n} rotmod files -> {ROTMOD_DIR}/")
    return n


def fetch_missing():
    if os.path.exists(EDR3_DST):
        print(f"  present: {EDR3_DST}")
    else:
        print(f"  fetching the EDR3 catalog (1.4 GB) ...")
        download(EDR3_URL, EDR3_DST, EDR3_SHA, EDR3_SIZE)
    if os.path.exists(MRT_DST):
        print(f"  present: {MRT_DST}")
    else:
        download(MRT_URL, MRT_DST, MRT_SHA, MRT_SIZE)
    have = len(glob.glob(ROTMOD_GLOB, recursive=True))
    if have >= ROTMOD_COUNT:
        print(f"  present: {ROTMOD_DIR}/ ({have} files)")
    else:
        if not (os.path.exists(ZIP_DST)
                and sha256_file(ZIP_DST) == ZIP_SHA):
            download(ZIP_URL, ZIP_DST, ZIP_SHA, ZIP_SIZE)
        else:
            print(f"  present: {ZIP_DST}")
        extract_rotmod(ZIP_DST)


def verify_manifest():
    """Check every entry of data/MANIFEST.sha256; return True iff all pass."""
    if not os.path.exists(MANIFEST):
        print(f"MISSING {MANIFEST} -- run from a full clone")
        return False
    entries = []
    with open(MANIFEST) as f:
        for line in f:
            m = re.match(r"^([0-9a-f]{64})\s+(\d+)\s+(.+?)\s*$", line)
            if m:
                entries.append((m.group(1), int(m.group(2)), m.group(3)))
    print(f"verifying {len(entries)} manifest entries:")
    all_ok = True
    for want_sha, want_size, path in entries:
        agg = re.match(r"^data/sparc/rotmod/\*\* \((\d+) files\)$", path)
        if agg:
            # stage 9H's aggregate recipe, replicated verbatim
            rotm = sorted(glob.glob(ROTMOD_GLOB, recursive=True))
            hr = hashlib.sha256()
            for p in rotm:
                hr.update(sha256_file(p).encode())
            got_sha = hr.hexdigest()
            got_size = sum(os.path.getsize(p) for p in rotm)
            ok = (len(rotm) == int(agg.group(1))
                  and got_size == want_size and got_sha == want_sha)
            label = f"{path}  [{len(rotm)} found]"
        elif not os.path.exists(path):
            ok, label = False, f"{path}  [MISSING]"
        else:
            got_size = os.path.getsize(path)
            got_sha = sha256_file(path) if got_size == want_size else "-"
            ok = (got_size == want_size and got_sha == want_sha)
            label = path
        all_ok &= ok
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}")
    return all_ok


def selftest():
    tmp = tempfile.mkdtemp(prefix="fetch_inputs_selftest_")
    try:
        print(f"selftest (temporary dir, data/ untouched): {tmp}")
        download(MRT_URL, os.path.join(tmp, "SPARC_Lelli2016c.mrt"),
                 MRT_SHA, MRT_SIZE)
        zp = os.path.join(tmp, "Rotmod_LTG.zip")
        download(ZIP_URL, zp, ZIP_SHA, ZIP_SIZE)
        with zipfile.ZipFile(zp) as z:
            n = sum(1 for m in z.namelist()
                    if os.path.basename(m).endswith("_rotmod.dat"))
        ok = (n == ROTMOD_COUNT)
        print(f"  zip holds {n} rotmod members "
              f"(expect {ROTMOD_COUNT}) -> {'PASS' if ok else 'FAIL'}")
        print(f"==> SELFTEST {'PASS' if ok else 'FAIL'}")
        return ok
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "--selftest":
        sys.exit(0 if selftest() else 1)
    if mode != "--verify":
        print("fetch:")
        fetch_missing()
    ok = verify_manifest()
    if ok:
        print("==> ALL VERIFIED -- inputs are byte-identical to the papers'")
    else:
        print("==> VERIFICATION FAILED -- do not trust downstream numbers "
              "until resolved")
    sys.exit(0 if ok else 1)
