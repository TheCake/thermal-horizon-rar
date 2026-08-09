"""
FETCH: MUSE-DARK UDF (Papers I & III) public data release
=========================================================
Source: https://dark-matter.osu-lyon.fr/  (linked from
https://dark.univ-lyon1.fr/data-releases, footnote 3 of Ciocan et al. 2026,
arXiv:2604.22613 = MUSE-DARK III, the RAR at 0.33<z<1.44).

Products fetched (all public, no auth):
  data/musedark_release/
    *_bestfit.txt                    -- 7 disk-halo decomposition catalogs
    Fit_statistics_all_models.txt    -- per-galaxy chi2/BIC/AIC/log_Z, 7 models
    photometry_catalogue.txt         -- muse_id z r_kpc incl PA Mstar eMstar SFR eSFR
    ID%04d/DC14_%d_true_Vrot.dat     -- intrinsic RC: dx_arcsec rad_Re flux_slit v_kms sig_kms
    ID%04d/DC14_%d_galaxy_parameters.{dat,txt}
    ID%04d/DC14_%d_derived_parameters.{dat,txt}
    ID%04d/DC14_%d_model.txt
    ID%04d/mass_%s.fits              -- baryonic mass map (pattern probed at runtime)

Everything under data/ is gitignored (repo convention); this script is the
re-fetch documentation. Idempotent: skips files already present and non-empty.
"""
import os
import sys
import time
import urllib.request

BASE = "https://dark-matter.osu-lyon.fr/data"
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                    "data", "musedark_release")
os.makedirs(ROOT, exist_ok=True)

CATALOGS = ["DC14_bestfit.txt", "NFW_bestfit.txt", "cNFW_bestfit.txt",
            "Einastot_bestfit.txt", "DZF_bestfit.txt", "Burkert_bestfit.txt",
            "baryons_only_bestfit.txt", "Fit_statistics_all_models.txt",
            "photometry_catalogue.txt"]


def get(url, dest, retries=3):
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        return "cached"
    for i in range(retries):
        try:
            urllib.request.urlretrieve(url, dest)
            return "ok"
        except Exception as e:
            err = str(e)
            time.sleep(0.5 * (i + 1))
    if os.path.exists(dest):
        os.remove(dest)
    return f"FAIL ({err})"


def main():
    log = []
    for c in CATALOGS:
        st = get(f"{BASE}/catalogues/{c}", os.path.join(ROOT, c))
        log.append(f"catalog {c}: {st}")

    ids = []
    with open(os.path.join(ROOT, "DC14_bestfit.txt")) as f:
        next(f)
        for line in f:
            i = int(line.split()[0])
            if i not in ids:
                ids.append(i)
    print(f"{len(ids)} unique muse_ids in DC14 catalog")

    # probe the mass-map name pattern on the first galaxy
    mm_pat = None
    for pat in ("mass_{i}.fits", "mass_{i:04d}.fits", "mass_ID{i:04d}.fits"):
        p = pat.format(i=ids[0])
        st = get(f"{BASE}/ID{ids[0]:04d}/mass_maps/{p}",
                 os.path.join(ROOT, f"ID{ids[0]:04d}", p))
        os.makedirs(os.path.join(ROOT, f"ID{ids[0]:04d}"), exist_ok=True)
        st = get(f"{BASE}/ID{ids[0]:04d}/mass_maps/{p}",
                 os.path.join(ROOT, f"ID{ids[0]:04d}", p))
        if st in ("ok", "cached"):
            mm_pat = pat
            break
    log.append(f"mass-map pattern: {mm_pat}")

    nok = nfail = 0
    for i in ids:
        d = os.path.join(ROOT, f"ID{i:04d}")
        os.makedirs(d, exist_ok=True)
        files = [f"galpak_run_DC14/DC14_{i}_true_Vrot.dat",
                 f"galpak_run_DC14/DC14_{i}_galaxy_parameters.dat",
                 f"galpak_run_DC14/DC14_{i}_galaxy_parameters.txt",
                 f"galpak_run_DC14/DC14_{i}_derived_parameters.dat",
                 f"galpak_run_DC14/DC14_{i}_derived_parameters.txt",
                 f"galpak_run_DC14/DC14_{i}_model.txt"]
        if mm_pat:
            files.append(f"mass_maps/{mm_pat.format(i=i)}")
        for f in files:
            st = get(f"{BASE}/ID{i:04d}/{f}", os.path.join(d, os.path.basename(f)))
            if st == "FAIL" or st.startswith("FAIL"):
                nfail += 1
                log.append(f"ID{i:04d}/{os.path.basename(f)}: {st}")
            else:
                nok += 1
    log.append(f"per-galaxy files: {nok} ok/cached, {nfail} failed")

    with open(os.path.join(ROOT, "fetch_log.txt"), "w") as f:
        f.write("\n".join(log) + "\n")
    print("\n".join(log[-25:]))


if __name__ == "__main__":
    sys.exit(main())
