#!/usr/bin/env python3
"""Cross-match Rubin (VRO) designations with the Minor Planet Center.

Stage: crossmatches
Upstream: inputs/obs_sbn_X05_full.parquet, inputs/mpcorb_extended.dat.gz
Downstream: orbits, states (and photometry via mv.load_moi_photometry)

Usage (always from the pipeline/ directory):

    python crossmatches.py "2026 DO29"
    python crossmatches.py 2026_DO29
    python crossmatches.py all

For every Rubin object (provid_vro) finds the designation of the same object
in the MPC orbit catalogue (provid_mpc, e.g. "(1) Ceres", "(800001) 2013 YW88"
or "2021 AQ"). Search order: permid in MPCORB, provid in MPCORB, provid among
the other designations of MPCORB, and finally the MPC Designation Identifier
API. Rows in crossmatches.csv are always appended - delete the CSV yourself
to rebuild the cross-match from scratch (e.g. after a new MPCORB download).
"""

import gzip
import json
import os
import re
import sys
import time
import urllib.request

import pandas as pd

import mvroc as mv

# ---------------------------------------------------------------------------
# Crossmatch constants
# ---------------------------------------------------------------------------
MPC_API_URL = "https://data.minorplanetcenter.net/api/query-identifier"
API_BATCH = 100  # identifiers per MPC API request
API_TIMEOUT = 120  # seconds per MPC API request
DESIG_COLS = (166, 194)  # MPCORB readable designation, e.g. "(1) Ceres"
OTHER_COL = 217  # MPCORB extended: other designations start here
OTHER_WIDTH = 11  # each other designation is a fixed 11-character field
OUT_DIR = "outputs/crossmatches"
RESULTS_CSV = f"{OUT_DIR}/crossmatches.csv"
LOG_FILE = f"{OUT_DIR}/crossmatches.log"
API_CACHE = f"{OUT_DIR}/mpc_query_identifier.json"
RESULTS_COLUMNS = [
    "provid_vro",
    "provid_mpc",
    "permid",
    "match_kind",
    "n_vro_for_mpc",
    "match_by",
]

# ---------------------------------------------------------------------------
# MPCORB / MPC API helpers
# ---------------------------------------------------------------------------
def split_designation(designation):
    """Split an MPCORB designation into permid and provisional designation.

    Parameters
    ----------
    designation : str
        Readable MPCORB designation.

    Returns
    -------
    permid : str or None
        Number of the object, None if unnumbered.
    provid : str or None
        Provisional (or survey) designation, None if the object is named.

    Examples
    --------
    >>> split_designation("(800001) 2013 YW88")
    ('800001', '2013 YW88')
    >>> split_designation("(1) Ceres")
    ('1', None)
    >>> split_designation("2021 AQ")
    (None, '2021 AQ')
    """
    match = re.match(r"^\((\d+)\)\s*(.*)$", designation)
    if match is None:
        return None, designation
    permid, rest = match.groups()
    if re.match(r"^(\d{4} [A-Z]{2}\d*|\d{4} [PT]-[L123])$", rest):
        return permid, rest
    return permid, None

def unique_map(pairs):
    """Dict key -> designation, dropping keys that point to several objects.

    Parameters
    ----------
    pairs : list of (str, str)
        (identifier, MPCORB designation) pairs.

    Returns
    -------
    dict
        Identifier -> designation for the identifiers that appear once.

    Examples
    --------
    >>> unique_map([("2025 MM27", "2007 CL70"), ("A", "x"), ("A", "y")])
    {'2025 MM27': '2007 CL70'}
    """
    table = pd.DataFrame(pairs, columns=["key", "designation"])
    table = table.drop_duplicates("key", keep=False)
    return dict(zip(table["key"], table["designation"]))

def load_mpcorb_keys(path):
    """Lookup tables key -> MPCORB designation (permid, provid, other).

    Ambiguous provisional designations (pointing to more than one object)
    are dropped.

    Parameters
    ----------
    path : str
        MPCORB extended file (gzip text), e.g. inputs/mpcorb_extended.dat.gz.

    Returns
    -------
    dict
        permid, provid and other: dicts from identifier to designation.
        n_objects: number of MPCORB objects read.

    Examples
    --------
    >>> keys = load_mpcorb_keys("inputs/mpcorb_extended.dat.gz")
    >>> keys["permid"]["1"], keys["other"]["2025 MM27"]
    ('(1) Ceres', '2007 CL70')
    """
    by_permid = {}
    provid_pairs = []
    other_pairs = []
    n_objects = 0
    with gzip.open(path, "rt") as f:
        # --> The data start after the header line of dashes.
        for line in f:
            if line.startswith("-----"):
                break
        for line in f:
            designation = line[DESIG_COLS[0] : DESIG_COLS[1]].strip()
            if not designation:
                continue
            n_objects += 1
            permid, provid = split_designation(designation)
            if permid is not None:
                by_permid[permid] = designation
            if provid is not None:
                provid_pairs.append((provid, designation))
            others = line[OTHER_COL:].rstrip("\n")
            for i in range(0, len(others), OTHER_WIDTH):
                other = others[i : i + OTHER_WIDTH].strip()
                if other:
                    other_pairs.append((other, designation))

    # --> A designation must point to one object: drop the ambiguous ones.
    return {
        "permid": by_permid,
        "provid": unique_map(provid_pairs),
        "other": unique_map(other_pairs),
        "n_objects": n_objects,
    }

def mpcorb_lookup(provid, permid, keys):
    """Find an object in MPCORB by permid, provid, then other designations.

    Parameters
    ----------
    provid : str
        Rubin provisional designation.
    permid : str or None
        Rubin permanent number.
    keys : dict
        Lookup tables from load_mpcorb_keys.

    Returns
    -------
    designation : str or None
        MPCORB designation, None if not found.
    match_by : str or None
        permid, provid or other_desig.

    Examples
    --------
    >>> mpcorb_lookup("2025 MM27", None, keys)
    ('2007 CL70', 'other_desig')
    """
    if permid is not None and permid in keys["permid"]:
        return keys["permid"][permid], "permid"
    if provid in keys["provid"]:
        return keys["provid"][provid], "provid"
    if provid in keys["other"]:
        return keys["other"][provid], "other_desig"
    return None, None

def query_mpc_identifier(ids, cache):
    """Ask the MPC Designation Identifier API for the ids not in the cache.

    Answers are added to cache and written to API_CACHE after every batch,
    so each identifier is only queried once. Stops quietly (keeping what
    was cached) if the API cannot be reached.

    Parameters
    ----------
    ids : list of str
        Identifiers to resolve (provisional designations).
    cache : dict
        Identifier -> API answer; updated in place.

    Returns
    -------
    int
        Number of identifiers queried in this call.

    Examples
    --------
    >>> cache = {}
    >>> query_mpc_identifier(["2026 DO29"], cache)
    1
    >>> cache["2026 DO29"]["unpacked_primary_provisional_designation"]
    '2015 CY86'
    """
    to_query = [i for i in ids if i not in cache]
    n_done = 0
    for start in range(0, len(to_query), API_BATCH):
        batch = to_query[start : start + API_BATCH]
        request = urllib.request.Request(
            MPC_API_URL,
            data=json.dumps({"ids": batch}).encode(),
            headers={"Content-Type": "application/json"},
            method="GET",
        )
        try:
            with urllib.request.urlopen(request, timeout=API_TIMEOUT) as response:
                cache.update(json.loads(response.read().decode()))
        except Exception as error:
            print(f"  MPC API unreachable ({error}); continuing with cache only")
            break
        n_done += len(batch)
        with open(API_CACHE, "w", encoding="utf-8") as f:
            json.dump(cache, f)
        print(f"  queried {n_done:,} of {len(to_query):,}")
    return n_done

def api_designation(answer, keys):
    """MPC designation (MPCORB format) from one Identifier API answer.

    If MPCORB has the object under its current number or primary
    designation, that MPCORB designation is used; otherwise it is built in
    the same format. Comets keep the IAU designation given by the MPC.

    Parameters
    ----------
    answer : dict or None
        One entry of the MPC Identifier API response.
    keys : dict
        Lookup tables from load_mpcorb_keys.

    Returns
    -------
    designation : str or None
        MPC designation, None if the MPC did not find the object.
    is_comet : bool
        True if the object is not a minor planet.

    Examples
    --------
    >>> api_designation(cache["2026 DO29"], keys)
    ('2015 CY86', False)
    """
    if not answer or answer.get("found") != 1:
        return None, False
    if answer["object_type"][0] != "Minor Planet":
        return answer["iau_designation"], True
    permid = answer["permid"]
    provid = answer["unpacked_primary_provisional_designation"]
    if permid in keys["permid"]:
        return keys["permid"][permid], False
    if provid in keys["provid"]:
        return keys["provid"][provid], False
    if permid is not None:
        return f"({permid}) {answer['name'] or provid}", False
    return provid, False

def print_crossmatch(provid, table, written=True):
    """Print the cross-match of one Rubin object and its Rubin aliases.

    Parameters
    ----------
    provid : str
        Rubin designation (provid_vro).
    table : DataFrame
        Cross-match of all objects (RESULTS_COLUMNS).
    written : bool, optional
        True if the row was appended now, False if it was already in
        crossmatches.csv.

    Returns
    -------
    None

    Examples
    --------
    >>> print_crossmatch("2025 MN83", table)
    ====================
    Crossmatch summary for 2025 MN83
    ====================
      provid_mpc: (831800) 2009 WC150
      ...
      Rubin ids for this MPC object: 2 (2009 WC150, 2025 MN83)
    """
    print()
    print("====================")
    print(f"Crossmatch summary for {provid}")
    print("====================")
    match = table[table["provid_vro"] == provid]
    if match.empty:
        print("  No MPC designation found (not written to crossmatches.csv)")
        return
    row = match.iloc[0]
    print(f"  provid_mpc: {row['provid_mpc']}")
    print(f"  permid: {row['permid'] if pd.notna(row['permid']) else '-'}")
    print(f"  match_kind: {row['match_kind']}")
    print(f"  match_by: {row['match_by']}")
    aliases = table.loc[table["provid_mpc"] == row["provid_mpc"], "provid_vro"]
    print(f"  Rubin ids for this MPC object: {len(aliases)} ({', '.join(aliases)})")
    if written:
        print(f"  Appended to {RESULTS_CSV}")
    else:
        print(f"  Already in {RESULTS_CSV}")

# ---------------------------------------------------------------------------
# Science: one asteroid
# ---------------------------------------------------------------------------
def mvroc_crossmatch_asteroid(provid, permid, keys, cache):
    """MPC designation of one Rubin object.

    Parameters
    ----------
    provid : str
        Rubin provisional designation (provid_vro).
    permid : str or None
        Rubin permanent number.
    keys : dict
        Lookup tables from load_mpcorb_keys.
    cache : dict
        MPC Identifier API answers (see query_mpc_identifier).

    Returns
    -------
    dict or None
        provid_vro, provid_mpc, permid, match_kind and match_by.
        None if the object has no MPC designation.

    Examples
    --------
    >>> mvroc_crossmatch_asteroid("2025 MM27", None, keys, {})
    {'provid_vro': '2025 MM27', 'provid_mpc': '2007 CL70', 'permid': None,
     'match_kind': 'provisional', 'match_by': 'other_desig'}
    """
    # --> First MPCORB (permid, provid, other designations), then the MPC API.
    designation, match_by = mpcorb_lookup(provid, permid, keys)
    is_comet = False
    if designation is None:
        designation, is_comet = api_designation(cache.get(provid), keys)
        match_by = "mpc_api"
    if designation is None:
        return None

    # --> Number of the MPC object and kind of designation.
    mpc_permid, _ = split_designation(designation)
    if is_comet:
        match_kind = "comet"
    elif mpc_permid is not None:
        match_kind = "numbered"
    else:
        match_kind = "provisional"
    return {
        "provid_vro": provid,
        "provid_mpc": designation,
        "permid": mpc_permid,
        "match_kind": match_kind,
        "match_by": match_by,
    }

# ---------------------------------------------------------------------------
# Crossmatches driver
# ---------------------------------------------------------------------------
def mvroc_crossmatches():
    """Run the Rubin <-> MPC cross-match from the command line.

    Must be launched from the pipeline/ directory.

    Stage: crossmatches
    Upstream: inputs (obs_sbn_X05_full.parquet, mpcorb_extended.dat.gz)
    Downstream: orbits, states, photometry (dual-provid merge)

    Command-line
    ------------
    Exactly one positional argument:

    * provid - one Rubin designation, e.g. "2026 DO29" or 2026_DO29
      (underscores are turned into spaces).
    * all - every Rubin object with a provid (skipping those already in
      crossmatches.csv).

    Examples
    --------
    python crossmatches.py "2026 DO29"
    python crossmatches.py all

    Inputs
    ------
    inputs/obs_sbn_X05_full.parquet
        Rubin detections (columns permid, provid).
    inputs/mpcorb_extended.dat.gz
        MPCORB extended (designation and other designations).
    MPC Designation Identifier API
        Only for objects not found in MPCORB; answers cached in
        outputs/crossmatches/mpc_query_identifier.json.

    Outputs
    -------
    outputs/crossmatches/crossmatches.csv
        One appended row per matched provid_vro (created with header if
        missing). n_vro_for_mpc counts the Rubin provids that share the
        same provid_mpc. Objects without an MPC designation are not written
        (they are retried on the next run).
    outputs/crossmatches/mpc_query_identifier.json
        Cache of the MPC API answers.
    outputs/crossmatches/crossmatches.log
        Same text that is printed on screen (append mode).

    Returns
    -------
    None
        With all, prints the catalogue summary (matched / unmatched /
        aliases / MPCORB vs API). With one provid, prints its MPC
        designation, how it was matched and all the Rubin ids of the same
        MPC object.
    """

    # --> Check command-line arguments.
    if len(sys.argv) != 2:
        print('Usage: python crossmatches.py "<provid>" | python crossmatches.py all')
        raise SystemExit(2)

    # --> Get the target provid.
    target = sys.argv[1].strip()
    if target.lower() != "all":
        target = target.replace("_", " ")

    # --> Check the upstream inputs.
    for path in (mv.OBSERVATIONS_PARQUET, mv.MPCORB_DAT):
        if not os.path.isfile(path):
            raise SystemExit(f"Missing input: {path} (run make get_gdrive_data)")

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

    # --> Distinct Rubin objects: one (provid, permid) per provid.
    print(f"Reading {mv.OBSERVATIONS_PARQUET} ...")
    obs = pd.read_parquet(mv.OBSERVATIONS_PARQUET, columns=["permid", "provid"])
    n_detections = len(obs)
    n_no_provid = int(obs["provid"].isna().sum())
    objects = obs.dropna(subset=["provid"]).drop_duplicates("provid")
    del obs
    objects = objects.astype(object).where(objects.notna(), None)
    print(f"  {n_detections:,} detections, {len(objects):,} objects with provid")
    print(f"  {n_no_provid:,} detections without provid (not cross-matched)")

    # --> If not processing all objects, select the target (all objects are
    #     still cross-matched below so that its aliases are found).
    selected = objects
    if target.lower() != "all":
        selected = objects[objects["provid"] == target]
        if selected.empty:
            raise SystemExit(f"Provid not found in observations: {target}")

    # --> MPCORB lookup tables (permid, provid, other designations).
    print(f"Reading {mv.MPCORB_DAT} ...")
    keys = load_mpcorb_keys(mv.MPCORB_DAT)
    print(
        f"  {keys['n_objects']:,} MPC objects, "
        f"{len(keys['other']):,} other designations"
    )

    # --> Ask the MPC API for the selected objects not found in MPCORB
    #     (the other objects only use the answers already cached).
    cache = {}
    if os.path.isfile(API_CACHE):
        with open(API_CACHE, encoding="utf-8") as f:
            cache = json.load(f)
    missing = [
        provid
        for provid, permid in zip(selected["provid"], selected["permid"])
        if mpcorb_lookup(provid, permid, keys)[0] is None
    ]
    print(f"{len(missing):,} objects not in MPCORB: asking the MPC API")
    print(f"  {sum(i in cache for i in missing):,} already in {API_CACHE}")
    query_mpc_identifier(missing, cache)

    # --> Cross-match every object (needed for n_vro_for_mpc).
    rows = []
    unmatched = []
    for provid, permid in zip(objects["provid"], objects["permid"]):
        row = mvroc_crossmatch_asteroid(provid, permid, keys, cache)
        if row is None:
            unmatched.append(provid)
        else:
            rows.append(row)
    table = pd.DataFrame(rows, columns=RESULTS_COLUMNS)
    table["permid"] = pd.to_numeric(table["permid"]).astype("Int64")

    # --> Several Rubin provids can point to the same MPC object.
    table["n_vro_for_mpc"] = table.groupby("provid_mpc")["provid_vro"].transform(
        "size"
    )

    # --> Append only the selected objects not yet in crossmatches.csv.
    done = mv.provids_in_outputs(RESULTS_CSV, "provid_vro")
    new = table[
        table["provid_vro"].isin(selected["provid"])
        & ~table["provid_vro"].isin(done)
    ]
    new.to_csv(RESULTS_CSV, mode="a", header=False, index=False)

    # --> Single object: show its cross-match.
    if target.lower() != "all":
        print_crossmatch(target, table, written=not new.empty)
        print(f"Done in {(time.perf_counter() - t0) / 60.0:.1f} m -> {OUT_DIR}")
        return

    # --> Catalogue summary.
    by = table["match_by"].value_counts()
    n_mpcorb = int(by.get("permid", 0) + by.get("provid", 0) + by.get("other_desig", 0))
    aliases = table[table["n_vro_for_mpc"] > 1]
    print()
    print("====================")
    print("Crossmatch summary")
    print("====================")
    print(f"  Rubin objects (provid_vro): {len(objects):,}")
    print(f"  Matched with an MPC designation: {len(table):,}")
    print(f"  Without MPC designation: {len(unmatched):,}")
    print(f"  Found in MPCORB: {n_mpcorb:,}")
    print(f"    by permid: {int(by.get('permid', 0)):,}")
    print(f"    by provid: {int(by.get('provid', 0)):,}")
    print(f"    by other designation: {int(by.get('other_desig', 0)):,}")
    print(f"  Found through the MPC API: {int(by.get('mpc_api', 0)):,}")
    print(f"  Distinct MPC objects (provid_mpc): {table['provid_mpc'].nunique():,}")
    print(
        f"  MPC objects with several Rubin provids: "
        f"{aliases['provid_mpc'].nunique():,} "
        f"({len(aliases):,} provid_vro)"
    )
    print(f"  Match kind: {table['match_kind'].value_counts().to_dict()}")
    print(
        f"  Appended: {len(new):,}, already in "
        f"{os.path.basename(RESULTS_CSV)}: {len(table) - len(new):,}"
    )
    if unmatched:
        print(f"  Unmatched (first 20): {unmatched[:20]}")
    print(f"Done in {(time.perf_counter() - t0) / 60.0:.1f} m -> {OUT_DIR}")

if __name__ == "__main__":
    mvroc_crossmatches()
