"""Fetch the X-COP public data release (P13-stage leg A).

Source: the X-COP project data page
  https://dominiqueeckert.wixsite.com/xcop/data
which (as of 2026-08-11) hosts the entire release as a single archive on
SWITCHdrive (the older ISDC mirror www.isdc.unige.ch/~deckert/XCOP/ times
out). Share token j3WUOYXWgv9Jbnz, single file, 315,080,566 bytes.

Contents per the data page: images/exposure/background maps, XMM density
profiles, spectral fits, Planck SZ profiles, thermodynamic profiles
(Ghirardini+19, A&A 621, A41), mass profiles (Ettori+19, A&A 621, A39),
gas mass / gas fraction (Eckert+19, A&A 621, A40), metal abundances
(Ghizzardi+21), stellar mass profiles (added 2020-08-19).
Clusters: A85, A644, A1644, A1795, A2029, A2142, A2255, A2319, A3158,
A3266, RXC1825, ZW1215.

Destination: data/xcop/ (gitignored; re-run this script to re-fetch).
"""
import os
import sys
import urllib.request

URL = "https://drive.switch.ch/index.php/s/j3WUOYXWgv9Jbnz/download"
EXPECTED = 315080566
DEST_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "xcop")
DEST = os.path.join(DEST_DIR, "xcop_release.archive")


def main():
    os.makedirs(DEST_DIR, exist_ok=True)
    if os.path.exists(DEST) and os.path.getsize(DEST) == EXPECTED:
        print(f"already fetched: {DEST} ({EXPECTED} bytes)")
        return
    req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0 (research fetch)"})
    print(f"downloading {URL} -> {DEST}")
    with urllib.request.urlopen(req, timeout=120) as r, open(DEST, "wb") as f:
        total = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            total += len(chunk)
            if total % (20 << 20) < (1 << 20):
                print(f"  {total/1e6:.0f} MB", flush=True)
    size = os.path.getsize(DEST)
    print(f"done: {size} bytes (expected {EXPECTED})")
    if size != EXPECTED:
        print("WARNING: size mismatch vs the probe -- inventory before use")
        sys.exit(1)


if __name__ == "__main__":
    main()
