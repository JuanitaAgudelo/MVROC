#!/usr/bin/env python
"""Query JPL SBDB physical properties, including the rotation period.

SBDB answers one object per HTTP request (astroquery SBDB.query(..., phys=True)).
A burst of parallel requests makes the server answer 403 Forbidden, so this
script sends one request at a time and waits between them. It stops if several
requests in a row are refused. The shared astroquery cache is left off.

The default target list is every asteroid provid in df_set, the same filter as
search_spectra_rotation_orbit.ipynb. Writes data/sbdb_physical_properties.csv.
A later run skips provids already saved in that file.
"""

import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
from astroquery.jplsbdb import SBDB

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
OBS_FILE = DATA_DIR / "obs_sbn_X05_full.parquet"
ROTATION_CATALOG = DATA_DIR / "asteroids-with-rotation.csv"
OUT_FILE = DATA_DIR / "sbdb_physical_properties.csv"
ASTEROID_PROVID = re.compile(r"^\d{4} [A-Z]{2}\d+[A-Z]?$")

WORKERS = 1
PAUSE_S = 1.5
RETRIES = 3
FORBIDDEN_PAUSE_S = 30
STOP_AFTER_FORBIDDEN = 8
COLUMNS = ["provid", "fullname", "parameter", "value", "sigma", "reference", "note"]


def format_elapsed(seconds):
    seconds = int(seconds)
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h {minutes:02d}m {secs:02d}s"
    return f"{minutes}m {secs:02d}s"


def targets_from_df_set():
    """Unique asteroid provids in df_set, the notebook's filtered observations."""
    provids = pd.read_parquet(OBS_FILE, columns=["provid"])["provid"]
    provids = provids.dropna().astype(str).str.strip()
    provids = provids[provids.map(lambda pid: bool(ASTEROID_PROVID.match(pid)))]
    provids = provids.drop_duplicates().sort_values()
    print(f"Number of objects with provid: {len(provids)}")
    return provids.tolist()


def targets_from_rotation():
    """Provids in asteroids-with-rotation.csv."""
    df_rot = pd.read_csv(ROTATION_CATALOG)
    provids = df_rot["provid"].dropna().astype(str).str.strip()
    return provids[provids != ""].drop_duplicates().tolist()


def targets_from_list():
    """Edit this list to query a few objects."""
    return [
        "2025 MM82",
        "2025 MA19",
    ]


def query_phys(provid):
    """One SBDB physical-property lookup. Returns rows for the output table."""
    last_error = None
    for attempt in range(1, RETRIES + 1):
        try:
            sbdb = SBDB.query(str(provid), phys=True, cache=False)
            break
        except Exception as exc:
            last_error = exc
            # 403 means JPL is refusing the client. Wait before trying again.
            if "403" in str(exc):
                time.sleep(FORBIDDEN_PAUSE_S)
            else:
                time.sleep(attempt)
    else:
        return [{
            "provid": provid,
            "fullname": "",
            "parameter": "",
            "value": "",
            "sigma": "",
            "reference": "",
            "note": f"query failed: {last_error}",
        }]

    phys = sbdb.get("phys_par") or {}
    obj = sbdb.get("object") or {}
    fullname = obj.get("fullname") or obj.get("des") or provid
    rows = []
    for key, value in phys.items():
        if key.endswith(("_sig", "_ref", "_note")) or value in (None, ""):
            continue
        rows.append({
            "provid": provid,
            "fullname": fullname,
            "parameter": key,
            "value": value,
            "sigma": phys.get(f"{key}_sig"),
            "reference": phys.get(f"{key}_ref"),
            "note": phys.get(f"{key}_note"),
        })
    if not rows:
        rows.append({
            "provid": provid,
            "fullname": fullname,
            "parameter": "",
            "value": "",
            "sigma": "",
            "reference": "",
            "note": "no physical properties",
        })
    return rows


def already_queried():
    """Provids with a saved result. Failed queries are dropped so they run again."""
    if not OUT_FILE.exists():
        return set()
    saved = pd.read_csv(OUT_FILE)
    failed = saved["note"].fillna("").astype(str).str.startswith("query failed")
    if failed.any():
        saved = saved.loc[~failed].copy()
        saved.to_csv(OUT_FILE, index=False)
    return set(saved["provid"].dropna().astype(str))


def append_rows(rows):
    frame = pd.DataFrame(rows, columns=COLUMNS)
    write_header = not OUT_FILE.exists()
    frame.to_csv(OUT_FILE, mode="a", header=write_header, index=False)


def main():
    # targets_from_df_set: every asteroid provid in the notebook's df_set
    # targets_from_rotation: asteroids-with-rotation.csv
    # targets_from_list: the list written inside that function
    targets = targets_from_df_set()
    done = already_queried()
    pending = [provid for provid in targets if provid not in done]
    print(f"{len(targets)} targets, {len(targets) - len(pending)} already saved, {len(pending)} to query")
    print(f"One request every {PAUSE_S:.1f} s ({WORKERS} worker)")
    if not pending:
        print(OUT_FILE)
        return

    t_start = time.perf_counter()
    n_ok = 0
    forbidden_in_a_row = 0
    stopped = False
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futures = {}
        pending_iter = iter(pending)
        for _ in range(WORKERS):
            provid = next(pending_iter, None)
            if provid is None:
                break
            futures[pool.submit(query_phys, provid)] = provid
        while futures:
            done = next(as_completed(futures))
            provid = futures.pop(done)
            rows = done.result()
            append_rows(rows)
            n_ok += 1
            note = rows[0]["note"] or ""
            if "403" in str(note):
                forbidden_in_a_row += 1
            else:
                forbidden_in_a_row = 0
            names = [row["parameter"] for row in rows if row["parameter"]]
            label = ", ".join(names) if names else note
            elapsed = time.perf_counter() - t_start
            estimated = elapsed * len(pending) / n_ok
            print(
                f"{provid}: {label} | "
                f"elapsed: {format_elapsed(elapsed)} ({n_ok}/{len(pending)}), "
                f"estimated total: {format_elapsed(estimated)}"
            )
            if forbidden_in_a_row >= STOP_AFTER_FORBIDDEN:
                print("SBDB is refusing requests (403). Stopping so the block can expire.")
                print("Saved objects are skipped on the next run.")
                stopped = True
                break
            time.sleep(PAUSE_S)
            provid = next(pending_iter, None)
            if provid is not None:
                futures[pool.submit(query_phys, provid)] = provid

    if stopped:
        print(f"Stopped after {n_ok} asteroids. File: {OUT_FILE}")
    else:
        print(f"Wrote {n_ok} asteroids to {OUT_FILE}")


if __name__ == "__main__":
    main()
