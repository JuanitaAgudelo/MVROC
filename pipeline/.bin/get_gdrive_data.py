#!/usr/bin/env python3
"""Download large pipeline data files from Google Drive with gdown.

Skips any destination that already exists locally as a non-empty file
(or as a symlink that resolves to a non-empty file).

Run from pipeline/ (or via make):

    cd pipeline
    python3 .bin/get_gdrive_data.py
    make get_gdrive_data

The Google Drive folder mvroc-pipeline-data must be shared so gdown can
read it. Folder id is PIPELINE_DATA_FOLDER_ID.

https://drive.google.com/drive/folders/1gddbkDYiGWRI4pVhKPfR5AGp5dNpTrBW?usp=drive_link
"""

from __future__ import annotations

import os
import sys

# ---------------------------------------------------------------------------
# Drive layout (mirrors pipeline/ under mvroc-pipeline-data/)
# ---------------------------------------------------------------------------
# https://drive.google.com/drive/folders/1gddbkDYiGWRI4pVhKPfR5AGp5dNpTrBW?usp=drive_link
PIPELINE_DATA_FOLDER_ID = "1gddbkDYiGWRI4pVhKPfR5AGp5dNpTrBW"

# drive_path: relative to the Drive folder root.
# local: relative to pipeline/ (cwd for drivers and this Makefile).
# Optional "id" entries are used when set; otherwise resolved from the folder.
FILES = [
    {
        "drive_path": "inputs/obs_sbn_X05_full.parquet",
        "local": "inputs/obs_sbn_X05_full.parquet",
        "id": None,
    },
    {
        "drive_path": "inputs/mpcorb_extended.dat.gz",
        "local": "inputs/mpcorb_extended.dat.gz",
        "id": None,
    },
    {
        "drive_path": "outputs/mois/mois.parquet",
        "local": "outputs/mois/mois.parquet",
        "id": None,
    },
    {
        "drive_path": "outputs/states/states.parquet",
        "local": "outputs/states/states.parquet",
        "id": None,
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


def _normalize_drive_path(path):
    return path.replace("\\", "/").lstrip("./")


def _list_folder_files(folder_id):
    """Return {relative_path: file_id} for files under the Drive folder."""
    import gdown

    entries = gdown.download_folder(
        id=folder_id,
        quiet=True,
        skip_download=True,
        use_cookies=True,
        remaining_ok=True,
    )
    mapping = {}
    for entry in entries or []:
        # GoogleDriveFileToDownload: .id, .path, .local_path
        rel = _normalize_drive_path(getattr(entry, "path", "") or "")
        file_id = getattr(entry, "id", None)
        if not rel or not file_id:
            continue
        mapping[rel] = file_id
        mapping[os.path.basename(rel)] = file_id
        parts = rel.split("/")
        if len(parts) > 1:
            mapping["/".join(parts[1:])] = file_id
    return mapping


def _download_file(file_id, dest):
    import gdown

    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    url = f"https://drive.google.com/uc?id={file_id}"
    out = gdown.download(url=url, output=dest, quiet=False, resume=True)
    if out is None or not _local_ready(dest):
        raise RuntimeError(f"gdown failed for id={file_id} -> {dest}")
    return dest


def main():
    pipe = _pipeline_dir()
    os.chdir(pipe)

    try:
        import gdown  # noqa: F401
    except ImportError:
        print("gdown is not installed. Run: make install  (from pipeline/)")
        raise SystemExit(1)

    pending = []
    for item in FILES:
        dest = os.path.join(pipe, item["local"])
        if _local_ready(dest):
            print(f"skip (already local): {item['local']}")
            continue
        pending.append(item)

    if not pending:
        print("All pipeline data files are already present.")
        return

    need_listing = any(not item.get("id") for item in pending)
    listing = {}
    if need_listing:
        print(
            f"Listing Google Drive folder {PIPELINE_DATA_FOLDER_ID} ..."
        )
        try:
            listing = _list_folder_files(PIPELINE_DATA_FOLDER_ID)
        except Exception as exc:
            print(f"Could not list Drive folder: {exc}")
            print(
                "Share mvroc-pipeline-data publicly (or with link) and retry, "
                "or set FILES[*]['id'] in .bin/get_gdrive_data.py."
            )
            raise SystemExit(1)
        if not listing:
            print(
                "Drive folder listing was empty. "
                "Files may still be uploading, or the folder is not shared."
            )
            raise SystemExit(1)

    errors = 0
    for item in pending:
        dest = os.path.join(pipe, item["local"])
        file_id = item.get("id")
        drive_path = _normalize_drive_path(item["drive_path"])
        if not file_id:
            file_id = listing.get(drive_path) or listing.get(
                os.path.basename(drive_path)
            )
        if not file_id:
            print(f"missing on Drive: {drive_path}")
            errors += 1
            continue
        print(f"download: {item['local']}  (id={file_id})")
        try:
            _download_file(file_id, dest)
            print(f"  wrote {item['local']} ({os.path.getsize(dest):,} bytes)")
        except Exception as exc:
            print(f"  ERROR: {exc}")
            errors += 1

    if errors:
        raise SystemExit(1)
    print("Done.")


if __name__ == "__main__":
    main()
