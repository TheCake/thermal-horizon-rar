# fetch_flynn_corpus.py — idempotent fetch of the Unified HI Rotation Curve Corpus
# (Flynn 2026, arXiv:2604.13489; Zenodo 10.5281/zenodo.19563417, CC BY 4.0)
# plus the Karachentsev Local Volume isolation-index table (VizieR J/AJ/145/101,
# Updated Nearby Galaxy Catalog: tidal index TI = Theta1) for the void/isolated
# cross-match feasibility read (10I candidate).
#
# Targets:
#   data/flynn_corpus/            — all Zenodo record files (byte tables)
#   data/flynn_corpus/zenodo_record.json — the record metadata (file list, sizes)
#   data/flynn2026_text.txt       — paper text (ar5iv HTML stripped; fallback: abs page)
#   data/karachentsev_ungc.csv    — LV catalog with tidal index Theta1 (+Theta5 if present)
#
# Scout-grade inputs only; every downstream use requires the primary read first
# (trap #6: read the byte tables before designing a stage).
import json
import os
import re
import sys
import time
import urllib.request

BASE = os.path.join(os.path.dirname(__file__), "..", "data")
CORP = os.path.join(BASE, "flynn_corpus")
os.makedirs(CORP, exist_ok=True)

UA = {"User-Agent": "Mozilla/5.0 (research fetch; PhysicsResearch repo)"}


def get(url, timeout=120):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read()


def save(path, blob):
    with open(path, "wb") as f:
        f.write(blob)
    print(f"  saved {os.path.relpath(path, BASE)} ({len(blob):,} bytes)")


def fetch_zenodo():
    rec_path = os.path.join(CORP, "zenodo_record.json")
    if os.path.exists(rec_path):
        rec = json.load(open(rec_path, encoding="utf-8"))
        print("[zenodo] record.json cached")
    else:
        print("[zenodo] fetching record 19563417 metadata")
        blob = get("https://zenodo.org/api/records/19563417")
        save(rec_path, blob)
        rec = json.loads(blob)
    files = rec.get("files", [])
    print(f"[zenodo] {len(files)} files in record")
    for f in files:
        name = f.get("key") or f.get("filename")
        link = (f.get("links") or {}).get("self")
        if link is None:
            link = f"https://zenodo.org/records/19563417/files/{name}?download=1"
        dest = os.path.join(CORP, name)
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            print(f"  cached {name}")
            continue
        print(f"  downloading {name} ({f.get('size', '?')} bytes)")
        for attempt in range(3):
            try:
                save(dest, get(link, timeout=600))
                break
            except Exception as e:
                print(f"    retry {attempt+1}: {e}")
                time.sleep(5)
        else:
            print(f"  FAILED {name}")


def fetch_paper():
    dest = os.path.join(BASE, "flynn2026_text.txt")
    if os.path.exists(dest) and os.path.getsize(dest) > 5000:
        print("[paper] text cached")
        return
    txt = None
    for url in ("https://ar5iv.labs.arxiv.org/html/2604.13489",
                "https://arxiv.org/abs/2604.13489"):
        try:
            html = get(url).decode("utf-8", errors="replace")
            body = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", html,
                          flags=re.S | re.I)
            body = re.sub(r"<[^>]+>", " ", body)
            body = re.sub(r"&[a-zA-Z]+;", " ", body)
            body = re.sub(r"[ \t]+", " ", body)
            body = re.sub(r"\n\s*\n+", "\n\n", body)
            txt = f"SOURCE: {url}\n\n" + body
            print(f"[paper] fetched {url} ({len(txt):,} chars)")
            if len(txt) > 20000:
                break
        except Exception as e:
            print(f"[paper] {url} failed: {e}")
    if txt:
        save(dest, txt.encode("utf-8"))
    else:
        print("[paper] FAILED both routes")


def fetch_karachentsev():
    dest = os.path.join(BASE, "karachentsev_ungc.csv")
    if os.path.exists(dest) and os.path.getsize(dest) > 5000:
        print("[ungc] cached")
        return
    # VizieR TAP: J/AJ/145/101/catalog = Updated Nearby Galaxy Catalog
    # (Karachentsev, Makarov & Kaisina 2013). Column TI = tidal index Theta1.
    query = ("SELECT * FROM \"J/AJ/145/101/catalog\"")
    url = ("https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync?"
           "REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=" +
           urllib.request.quote(query))
    try:
        blob = get(url, timeout=300)
        save(dest, blob)
    except Exception as e:
        print(f"[ungc] TAP failed: {e}")


if __name__ == "__main__":
    fetch_zenodo()
    fetch_paper()
    fetch_karachentsev()
    print("DONE")
