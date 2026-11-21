#!/usr/bin/env python3
"""MPC orbital elements and epoch of every cross-matched Rubin object.

Stage: orbits
Upstream: outputs/crossmatches/crossmatches.csv, inputs/mpcorb_extended.dat.gz
Downstream: states

Usage (always from the pipeline/ directory):

    python orbits.py "1981 US22"
    python orbits.py 1981_US22
    python orbits.py all

For every provid_vro of crossmatches.csv takes the osculating elements of its
MPC designation (provid_mpc) from MPCORB: a, e, q, incl, node, peri, M and the
epoch of the elements (epoch_tdb). Objects without an orbit in MPCORB (e.g.
comets) get a row with empty elements. Rows in orbits.csv are always appended
- delete the CSV yourself to rebuild it (e.g. after a new MPCORB download).
"""

import gzip
import os
import sys
import time

import numpy as np
import pandas as pd

import mvroc as mv

# ---------------------------------------------------------------------------
# Orbits constants
# ---------------------------------------------------------------------------
# Fixed-width MPCORB columns as 0-based [start, end) slices
# (https://www.minorplanetcenter.net/iau/info/MPOrbitFormat.html).
MPCORB_COLSPECS = {
    "epoch_packed": (20, 25),  # epoch of the elements, packed (e.g. K2669)
    "M": (26, 35),  # mean anomaly [deg]
    "peri": (37, 46),  # argument of perihelion [deg]
    "node": (48, 57),  # longitude of ascending node [deg]
    "incl": (59, 68),  # inclination [deg]
    "e": (70, 79),  # eccentricity
    "a": (92, 103),  # semimajor axis [au]
    "designation": (166, 194),  # readable designation, e.g. "(1) Ceres"
}
ELEMENTS = ["a", "e", "q", "incl", "node", "peri", "M"]
OUT_DIR = "outputs/orbits"
RESULTS_CSV = f"{OUT_DIR}/orbits.csv"
LOG_FILE = f"{OUT_DIR}/orbits.log"
RESULTS_COLUMNS = [
    "provid_vro",
    "provid_mpc",
    "permid",
    "epoch_tdb",
] + ELEMENTS

# ---------------------------------------------------------------------------
# MPCORB helpers
# ---------------------------------------------------------------------------
def unpack_epoch(packed):
    """Packed MPC epoch -> TDB date string at 0h.

    MPCORB epochs are 0h TT; TT and TDB differ by less than 2 ms, so the
    same date is written as TDB.

    Parameters
    ----------
    packed : str
        Packed epoch, e.g. K2669 (century letter, year, month, day).

    Returns
    -------
    str
        Epoch as YYYY-MM-DD HH:MM:SS (TDB).

    Examples
    --------
    >>> unpack_epoch("K2669")
    '2026-06-09 00:00:00'
    >>> unpack_epoch("J981A")
    '1998-01-10 00:00:00'
    """
    digit = lambda c: int(c) if c.isdigit() else ord(c) - ord("A") + 10
    year = {"I": 1800, "J": 1900, "K": 2000}[packed[0]] + int(packed[1:3])
    return f"{year}-{digit(packed[3]):02d}-{digit(packed[4]):02d} 00:00:00"

def load_mpcorb_elements(path):
    """Orbital elements of every MPCORB object, indexed by designation.

    Parameters
    ----------
    path : str
        MPCORB extended file (gzip text), e.g. inputs/mpcorb_extended.dat.gz.

    Returns
    -------
    DataFrame
        Index designation; columns epoch_tdb, a, e, q, incl, node, peri, M
        (au and degrees).

    Examples
    --------
    >>> mpcorb = load_mpcorb_elements("inputs/mpcorb_extended.dat.gz")
    >>> mpcorb.loc["(1) Ceres", ["epoch_tdb", "a"]].tolist()
    ['2026-06-09 00:00:00', 2.7655...]
    """
    # --> The data start after the header line of dashes.
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if line.startswith("-----"):
                header_lines = i + 1
                break

    mpcorb = pd.read_fwf(
        path,
        compression="gzip",
        skiprows=header_lines,
        header=None,
        colspecs=list(MPCORB_COLSPECS.values()),
        names=list(MPCORB_COLSPECS.keys()),
        dtype={"epoch_packed": str, "designation": str},
        skip_blank_lines=True,
    )
    mpcorb = mpcorb.dropna(subset=["designation"])

    # --> Epoch as a TDB string and perihelion distance q = a (1 - e).
    epochs = {p: unpack_epoch(p) for p in mpcorb["epoch_packed"].dropna().unique()}
    mpcorb["epoch_tdb"] = mpcorb["epoch_packed"].map(epochs)
    mpcorb["q"] = mpcorb["a"] * (1 - mpcorb["e"])
    return mpcorb.set_index("designation")[["epoch_tdb"] + ELEMENTS]

# ---------------------------------------------------------------------------
# Science: one asteroid
# ---------------------------------------------------------------------------
def mvroc_orbit_asteroid(match, mpcorb):
    """Orbits row for one cross-matched Rubin object.

    Parameters
    ----------
    match : dict-like
        One crossmatches.csv row (provid_vro, provid_mpc, permid).
    mpcorb : DataFrame
        Elements from load_mpcorb_elements.

    Returns
    -------
    dict
        provid_vro, provid_mpc, permid, epoch_tdb, a, e, q, incl, node,
        peri, M. Epoch and elements are empty (NaN) if provid_mpc has no
        orbit in MPCORB.

    Examples
    --------
    >>> row = {"provid_vro": "2026 DO29", "provid_mpc": "2015 CY86", "permid": None}
    >>> mvroc_orbit_asteroid(row, mpcorb)["epoch_tdb"]
    '2026-06-09 00:00:00'
    """
    row = {
        "provid_vro": match["provid_vro"],
        "provid_mpc": match["provid_mpc"],
        "permid": match["permid"],
    }
    if match["provid_mpc"] in mpcorb.index:
        row.update(mpcorb.loc[match["provid_mpc"]].to_dict())
    else:
        row.update({col: np.nan for col in ["epoch_tdb"] + ELEMENTS})
    return row

def print_orbit(row):
    """Print the orbital elements of one object.

    Parameters
    ----------
    row : dict-like
        One orbits.csv row (RESULTS_COLUMNS).

    Returns
    -------
    None

    Examples
    --------
    >>> print_orbit(mvroc_orbit_asteroid(match, mpcorb))
    ====================
    Orbit of 1981 US22
    ====================
      provid_mpc: (12215) Jessicalounsbury
      permid: 12215
      epoch: 2026-06-09 00:00:00 TDB
      a    = 2.30064950 au
      ...
    """
    units = {"a": "au", "q": "au", "e": ""}
    print()
    print("====================")
    print(f"Orbit of {row['provid_vro']}")
    print("====================")
    print(f"  provid_mpc: {row['provid_mpc']}")
    print(f"  permid: {row['permid'] if pd.notna(row['permid']) else '-'}")
    if pd.isna(row["a"]):
        print("  No orbit in MPCORB (empty elements)")
        return
    print(f"  epoch: {row['epoch_tdb']} TDB")
    for col in ELEMENTS:
        print(f"  {col:<4} = {float(row[col]):.8f} {units.get(col, 'deg')}")

# ---------------------------------------------------------------------------
# Orbits driver
# ---------------------------------------------------------------------------
def mvroc_orbits():
    """Extract MPC orbital elements from the command line.

    Must be launched from the pipeline/ directory.

    Stage: orbits
    Upstream: crossmatches (outputs/crossmatches/crossmatches.csv),
    inputs/mpcorb_extended.dat.gz
    Downstream: states

    Command-line
    ------------
    Exactly one positional argument:

    * provid - one Rubin designation (provid_vro), e.g. "1981 US22" or
      1981_US22 (underscores are turned into spaces).
    * all - every object of crossmatches.csv (skipping those already in
      orbits.csv).

    Examples
    --------
    python orbits.py "1981 US22"
    python orbits.py all

    Inputs
    ------
    outputs/crossmatches/crossmatches.csv
        provid_vro -> provid_mpc map (columns provid_vro, provid_mpc, permid).
    inputs/mpcorb_extended.dat.gz
        MPCORB: designation, epoch and osculating elements.

    Outputs
    -------
    outputs/orbits/orbits.csv
        One appended row per provid_vro (created with header if missing):
        epoch_tdb and a, e, q, incl, node, peri, M (au, deg). Objects with
        no orbit in MPCORB get empty elements.
    outputs/orbits/orbits.log
        Same text that is printed on screen (append mode).

    Returns
    -------
    None
        With all, prints a summary (with orbit / without orbit / skipped).
        With one provid, prints its epoch and orbital elements.
    """

    # --> Check command-line arguments.
    if len(sys.argv) != 2:
        print('Usage: python orbits.py "<provid>" | python orbits.py all')
        raise SystemExit(2)

    # --> Get the target provid.
    target = sys.argv[1].strip()
    if target.lower() != "all":
        target = target.replace("_", " ")

    # --> Check the upstream files.
    for path in (mv.CROSSMATCHES_CSV, mv.MPCORB_DAT):
        if not os.path.isfile(path):
            raise SystemExit(f"Missing upstream file: {path}")

    # --> Create the output directory and set up logging.
    os.makedirs(OUT_DIR, exist_ok=True)

    # --> Mirror all print() output to the log file and to the terminal.
    log = open(LOG_FILE, "a", encoding="utf-8")
    sys.stdout = mv.LogOutput(sys.__stdout__, log)

    # --> Initialize the results CSV if it doesn't exist.
    if not os.path.isfile(RESULTS_CSV):
        pd.DataFrame(columns=RESULTS_COLUMNS).to_csv(RESULTS_CSV, index=False)
        print(f"Created {RESULTS_CSV}")

    # --> Never append to a CSV written with other columns (e.g. an old file).
    mv.check_results_header(RESULTS_CSV, RESULTS_COLUMNS)
    print(f"Output: {OUT_DIR}")
    print(f"Log: {LOG_FILE}")
    t0 = time.perf_counter()

    # --> Cross-matched objects (one row per provid_vro).
    print(f"Reading {mv.CROSSMATCHES_CSV} ...")
    matches = pd.read_csv(
        mv.CROSSMATCHES_CSV,
        usecols=["provid_vro", "provid_mpc", "permid"],
        dtype={"permid": "Int64"},
    )
    matches = matches.dropna(subset=["provid_vro"])
    print(f"  {len(matches):,} objects")

    # --> If not processing all objects, filter to the target.
    if target.lower() != "all":
        matches = matches[matches["provid_vro"] == target]
        if matches.empty:
            raise SystemExit(f"Provid not found in crossmatches.csv: {target}")

    # --> Skip the objects already in orbits.csv.
    done = mv.provids_in_outputs(RESULTS_CSV, "provid_vro")
    todo = matches[~matches["provid_vro"].isin(done)]
    print(f"  {len(todo):,} to process, {len(matches) - len(todo):,} already in orbits.csv")

    # --> Single object already done: show its row from orbits.csv.
    if todo.empty and target.lower() != "all":
        done_rows = pd.read_csv(RESULTS_CSV, dtype={"permid": "Int64"})
        print_orbit(done_rows[done_rows["provid_vro"] == target].iloc[0])
        print(f"  Already in {RESULTS_CSV}")
        return
    if todo.empty:
        print(f"Done in {(time.perf_counter() - t0) / 60.0:.1f} m -> {OUT_DIR}")
        return

    # --> MPCORB elements indexed by designation.
    print(f"Reading {mv.MPCORB_DAT} ...")
    mpcorb = load_mpcorb_elements(mv.MPCORB_DAT)
    print(f"  {len(mpcorb):,} MPC orbits")

    # --> Orbit of every object.
    rows = [mvroc_orbit_asteroid(match, mpcorb) for match in todo.to_dict("records")]
    table = pd.DataFrame(rows, columns=RESULTS_COLUMNS)
    table["permid"] = table["permid"].astype("Int64")
    table.to_csv(RESULTS_CSV, mode="a", header=False, index=False)

    # --> Single object: show its elements.
    if target.lower() != "all":
        print_orbit(table.iloc[0])
        print(f"  Appended to {RESULTS_CSV}")
        print(f"Done in {(time.perf_counter() - t0) / 60.0:.1f} m -> {OUT_DIR}")
        return

    # --> Catalogue summary.
    no_orbit = table.loc[table["a"].isna(), "provid_vro"].tolist()
    print()
    print("====================")
    print("Orbits summary")
    print("====================")
    print(f"  Objects processed: {len(table):,}")
    print(f"  With MPC orbit: {len(table) - len(no_orbit):,}")
    print(f"  Without MPC orbit (empty elements): {len(no_orbit):,}")
    print(f"  Distinct MPC objects: {table['provid_mpc'].nunique():,}")
    print(f"  Epochs: {table['epoch_tdb'].value_counts().head(3).to_dict()}")
    print(f"  Appended: {len(table):,} rows to {RESULTS_CSV}")
    if no_orbit:
        print(f"  Without orbit (first 20): {no_orbit[:20]}")
    print(f"Done in {(time.perf_counter() - t0) / 60.0:.1f} m -> {OUT_DIR}")

if __name__ == "__main__":
    mvroc_orbits()
