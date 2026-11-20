#!/usr/bin/env python3
"""Download public raw inputs into pipeline/inputs/.

Fetches:

* obs_sbn_X05_full.parquet (Asteroid Institute / B612 public bucket)
* mpcorb_extended.dat.gz (Minor Planet Center)

Skips any destination that already exists locally as a non-empty file
(or as a symlink that resolves to a non-empty file).

Run from pipeline/ (or via make):

    cd pipeline
    python3 .bin/get_public_data.py
    make get_public_data
"""

from __future__ import annotations

import os
import sys
import urllib.request

FILES = [
    {
        "local": "inputs/obs_sbn_X05_full.parquet",
        "url": (
            "https://storage.googleapis.com/asteroid-institute-public/"
            "production/rubin/mpc/obs_sbn/full/parquet/"
            "obs_sbn_X05_full.parquet"
        ),
    },
    {
        "local": "inputs/mpcorb_extended.dat.gz",
        "url": (
            "https://minorplanetcenter.net/Extended_Files/"
            "mpcorb_extended.dat.gz"
        ),
    },
]


def _pipeline_dir():
    """pipeline/ directory (parent of .bin/)."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _local_ready(path):
    """True if path exists and points to a non-empty regular file."""
    if not os.path.lexists(path):
        return False
    try:
        return os.path.isfile(path) and os.path.getsize(path) > 0
    except OSError:
        return False


def _download(url, dest):
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    print(f"download: {dest}")
    print(f"  from {url}")
    tmp = dest + ".part"
    try:
        urllib.request.urlretrieve(url, tmp)
        os.replace(tmp, dest)
    except Exception:
        if os.path.isfile(tmp):
            os.remove(tmp)
        raise
    if not _local_ready(dest):
        raise RuntimeError(f"download failed: {url} -> {dest}")
    print(f"  wrote {dest} ({os.path.getsize(dest):,} bytes)")


def main():
    pipe = _pipeline_dir()
    os.chdir(pipe)

    errors = 0
    for item in FILES:
        dest = os.path.join(pipe, item["local"])
        if _local_ready(dest):
            print(f"skip (already local): {item['local']}")
            continue
        try:
            _download(item["url"], dest)
        except Exception as exc:
            print(f"  ERROR: {exc}")
            errors += 1

    if errors:
        raise SystemExit(1)
    print("Done.")


if __name__ == "__main__":
    main()
