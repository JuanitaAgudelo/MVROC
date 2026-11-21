#!/usr/bin/env python3
"""Geometry and absolute magnitude of every Rubin observation.

Stage: states
Upstream: outputs/orbits/orbits.csv, inputs/obs_sbn_X05_full.parquet,
inputs/kernels/ (SPICE)
Downstream: mois, photometry (mv.load_moi_photometry)

Usage (always from the pipeline/ directory):

    python states.py "1981 US22"
    python states.py 1981_US22
    python states.py all

For every observation of an object with an orbit, the MPC elements are turned
into a heliocentric state vector at their epoch (spiceypy.conics), propagated
to the observation time with two-body motion (spiceypy.prop2b), corrected for
light time, and used to compute r (Sun distance), delta (Earth distance),
alpha (phase angle) and H (HG system, G = 0.15, V = Rubin mag).

states.parquet is rewritten on every run. The orbits.csv row of each object
(provid_mpc, permid, epoch and elements) is hashed (orbit_hash): objects whose
row did not change keep their rows (flag_updated = 0); new or changed objects
are recomputed (flag_updated = 1).
"""

import hashlib
import os
import sys
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import spiceypy

import mvroc as mv

# ---------------------------------------------------------------------------
# States constants
# ---------------------------------------------------------------------------
KERNEL_FILES = [
    f"{mv.KERNELS}/naif0012.tls",  # leap seconds (UTC -> ET)
    f"{mv.KERNELS}/gm_de440.tpc",  # GM of the Sun (DE440)
    f"{mv.KERNELS}/de440s.bsp",  # planetary ephemeris (Earth position)
]
AU_KM = 149597870.7  # IAU 2012 astronomical unit [km]
G_SLOPE = 0.15  # HG slope parameter used for every object
N_LIGHT_TIME = 3  # light-time iterations
ELEMENTS = ["a", "e", "q", "incl", "node", "peri", "M"]
OUT_DIR = "outputs/states"
RESULTS_PARQUET = f"{OUT_DIR}/states.parquet"
LOG_FILE = f"{OUT_DIR}/states.log"
SCHEMA = pa.schema(
    [
        ("id", pa.int64()),
        ("provid_vro", pa.string()),
        ("provid_mpc", pa.string()),
        ("permid", pa.int64()),
        ("obstime", pa.timestamp("us")),
        ("obs_et", pa.float64()),
        ("deltat", pa.float64()),
        ("r", pa.float64()),
        ("delta", pa.float64()),
        ("alpha", pa.float64()),
        ("mag", pa.float64()),
        ("H", pa.float64()),
        ("flag_updated", pa.int8()),
        ("orbit_hash", pa.string()),
    ]
)
RESULTS_COLUMNS = SCHEMA.names

# ---------------------------------------------------------------------------
# Orbit / geometry helpers
# ---------------------------------------------------------------------------
def orbit_hash(orbit):
    """Short hash of the identity, epoch and elements of one orbit.

    Parameters
    ----------
    orbit : dict-like
        One orbits.csv row (provid_mpc, permid, epoch_tdb, a, e, q, incl,
        node, peri, M).

    Returns
    -------
    str
        16 hexadecimal characters; changes if provid_mpc, permid, the epoch
        or any MPCORB element (a, e, incl, node, peri, M) changes. q is
        derived from a and e and is left out; elements are rounded to 12
        significant digits, so the hash does not depend on how the CSV
        floats were written.

    Examples
    --------
    >>> orbit_hash({"provid_mpc": "(12215) Jessicalounsbury", "permid": 12215,
    ...             "epoch_tdb": "2026-06-09 00:00:00", "a": 2.3, "e": 0.15,
    ...             "q": 1.955, "incl": 6.5, "node": 76.2, "peri": 288.6,
    ...             "M": 309.6})
    '...'
    """
    values = [str(orbit[c]) for c in ["provid_mpc", "permid", "epoch_tdb"]]
    values += [f"{float(orbit[c]):.12g}" for c in ELEMENTS if c != "q"]
    text = "|".join(values)
    return hashlib.md5(text.encode()).hexdigest()[:16]

def epoch_state(orbit, mu):
    """Heliocentric state vector at the epoch of the elements.

    Parameters
    ----------
    orbit : dict-like
        One orbits.csv row: epoch_tdb and q [au], e, incl, node, peri, M [deg].
    mu : float
        GM of the Sun [km^3/s^2].

    Returns
    -------
    state : ndarray, shape (6,)
        x, y, z [km], vx, vy, vz [km/s], ecliptic J2000.
    epoch_et : float
        Epoch as ephemeris time (TDB seconds past J2000).

    Examples
    --------
    >>> state, epoch_et = epoch_state(orbits.iloc[0], mu)
    >>> np.linalg.norm(state[:3]) / AU_KM
    2.01...
    """
    angles = np.radians([orbit[c] for c in ["incl", "node", "peri", "M"]])
    elts = [orbit["q"] * AU_KM, orbit["e"], *angles, 0.0, mu]
    # With t0 = 0 and et = 0 the state is at the epoch of the elements.
    state = np.ascontiguousarray(spiceypy.conics(elts, 0.0))
    epoch_et = spiceypy.str2et(f"{orbit['epoch_tdb']} TDB")
    return state, epoch_et

def observer_geometry(states, r_earth):
    """Sun distance, Earth distance and phase angle with light-time correction.

    The asteroid position when the light left it is approximated as
    r(t - lt) = r(t) - v(t) lt (iterated N_LIGHT_TIME times); the error of
    this linear step is below 1 km.

    Parameters
    ----------
    states : ndarray, shape (n, 6)
        Heliocentric state of the asteroid at the observation times [km, km/s].
    r_earth : ndarray, shape (n, 3)
        Heliocentric position of the Earth at the same times [km].

    Returns
    -------
    r : ndarray
        Sun distance [au].
    delta : ndarray
        Earth distance [au].
    alpha : ndarray
        Phase angle (Sun-asteroid-Earth) [deg].

    Examples
    --------
    >>> states = np.array([[2 * AU_KM, 0, 0, 0, 0, 0]])
    >>> observer_geometry(states, np.array([[AU_KM, 0, 0]]))
    (array([2.]), array([1.]), array([0.]))
    """
    r_ast, v_ast = states[:, :3], states[:, 3:]
    lt = np.zeros(len(states))
    for _ in range(N_LIGHT_TIME):
        r_emit = r_ast - v_ast * lt[:, None]
        rho = r_emit - r_earth  # Earth -> asteroid
        dist = np.linalg.norm(rho, axis=1)
        lt = dist / spiceypy.clight()
    r_sun = np.linalg.norm(r_emit, axis=1)
    cos_alpha = np.einsum("ij,ij->i", r_emit, rho) / (r_sun * dist)
    alpha = np.degrees(np.arccos(np.clip(cos_alpha, -1, 1)))
    return r_sun / AU_KM, dist / AU_KM, alpha

def absolute_magnitude(V, r, delta, alpha, G=G_SLOPE):
    """H from the apparent magnitude with the HG phase function.

    Bowell et al. (1989):
    H = V - 5 log10(r delta) + 2.5 log10[(1 - G) Phi1 + G Phi2].

    Parameters
    ----------
    V : array-like
        Apparent magnitude.
    r, delta : array-like
        Sun and Earth distances [au].
    alpha : array-like
        Phase angle [deg].
    G : float, optional
        Slope parameter.

    Returns
    -------
    ndarray
        Absolute magnitude H.

    Examples
    --------
    >>> absolute_magnitude(np.array([15.0]), np.array([2.0]), np.array([1.0]),
    ...                    np.array([0.0]))
    array([13.49485002])
    """
    tan_half = np.tan(np.radians(alpha) / 2)
    phi1 = np.exp(-3.33 * tan_half**0.63)
    phi2 = np.exp(-1.87 * tan_half**1.22)
    return V - 5 * np.log10(r * delta) + 2.5 * np.log10((1 - G) * phi1 + G * phi2)

def print_state(provid, rows, recomputed):
    """Print the summary of one object's states.

    Parameters
    ----------
    provid : str
        Rubin designation (provid_vro).
    rows : DataFrame
        states.parquet rows of the object.
    recomputed : bool
        True if the rows were computed in this run.

    Returns
    -------
    None

    Examples
    --------
    >>> print_state("1981 US22", rows, True)
    ====================
    States of 1981 US22
    ====================
      provid_mpc: (12215) Jessicalounsbury
      ...
    """
    print()
    print("====================")
    print(f"States of {provid}")
    print("====================")
    print(f"  provid_mpc: {rows['provid_mpc'].iloc[0]}")
    print(f"  N = {len(rows)} observations")
    print(f"  obstime: {rows['obstime'].min()} -> {rows['obstime'].max()}")
    if rows["r"].isna().all():
        print("  No orbit in orbits.csv (empty geometry)")
    else:
        units = [("deltat", "d"), ("r", "au"), ("delta", "au"), ("alpha", "deg")]
        for col, unit in units:
            print(f"  {col}: {rows[col].min():.4f} -> {rows[col].max():.4f} {unit}")
        hs = rows["H"].dropna()
        if len(hs):
            err = hs.std() / np.sqrt(len(hs)) if len(hs) > 1 else np.nan
            print(
                f"  H = {hs.mean():.2f} +/- {err:.2f} "
                f"(std {hs.std():.2f}, N = {len(hs)})"
            )
    if recomputed:
        print(f"  Computed now (flag_updated = 1) -> {RESULTS_PARQUET}")
    else:
        print(f"  Orbit unchanged, kept from {RESULTS_PARQUET}")

# ---------------------------------------------------------------------------
# Science: one asteroid
# ---------------------------------------------------------------------------
def mvroc_state_asteroid(orbit, obs_ets, mags, r_earth, mu):
    """Geometry and H of every observation of one asteroid.

    Parameters
    ----------
    orbit : dict-like
        orbits.csv row of the asteroid (epoch_tdb and elements).
    obs_ets : ndarray
        Observation times as ephemeris time [s].
    mags : ndarray
        Rubin apparent magnitudes (NaN if missing).
    r_earth : ndarray, shape (n, 3)
        Heliocentric Earth position at obs_ets [km].
    mu : float
        GM of the Sun [km^3/s^2].

    Returns
    -------
    dict of ndarray
        deltat [d], r [au], delta [au], alpha [deg], H. All NaN if the
        orbit has no elements.

    Examples
    --------
    >>> out = mvroc_state_asteroid(orbit, obs_ets, mags, r_earth, mu)
    >>> out["r"][:2]
    array([2.0..., 2.0...])
    """
    n = len(obs_ets)
    if pd.isna(orbit["a"]):
        nan = np.full(n, np.nan)
        return {"deltat": nan, "r": nan, "delta": nan, "alpha": nan, "H": nan}

    # --> State at the epoch of the elements, propagated to each observation.
    state0, epoch_et = epoch_state(orbit, mu)
    dts = obs_ets - epoch_et
    states = np.array([spiceypy.prop2b(mu, state0, dt) for dt in dts])

    # --> Light-time corrected geometry and absolute magnitude.
    r, delta, alpha = observer_geometry(states, r_earth)
    return {
        "deltat": dts / 86400.0,
        "r": r,
        "delta": delta,
        "alpha": alpha,
        "H": absolute_magnitude(mags, r, delta, alpha),
    }

# ---------------------------------------------------------------------------
# States driver
# ---------------------------------------------------------------------------
def mvroc_states():
    """Compute the per-observation geometry and H from the command line.

    Must be launched from the pipeline/ directory.

    Stage: states
    Upstream: orbits (outputs/orbits/orbits.csv),
    inputs/obs_sbn_X05_full.parquet, inputs/kernels/
    Downstream: mois, photometry (mv.load_moi_photometry)

    Command-line
    ------------
    Exactly one positional argument:

    * provid - one Rubin designation (provid_vro), e.g. "1981 US22" or
      1981_US22 (underscores are turned into spaces).
    * all - every object of orbits.csv.

    Examples
    --------
    python states.py "1981 US22"
    python states.py all

    Inputs
    ------
    outputs/orbits/orbits.csv
        Epoch and elements of each provid_vro.
    inputs/obs_sbn_X05_full.parquet
        Rubin observations (id, provid, obstime, mag).
    inputs/kernels/naif0012.tls, gm_de440.tpc, de440s.bsp
        Leap seconds, GM of the Sun and Earth ephemeris.

    Outputs
    -------
    outputs/states/states.parquet
        One row per observation id of an object in orbits.csv (rewritten
        every run): A.5 columns plus orbit_hash. Objects with empty
        elements get empty geometry and H.
    outputs/states/states.log
        Same text that is printed on screen (append mode).

    Returns
    -------
    None
        With all, prints a summary (recomputed / unchanged / rows). With
        one provid, prints the ranges of its geometry and its mean H.
    """

    # --> Check command-line arguments.
    if len(sys.argv) != 2:
        print('Usage: python states.py "<provid>" | python states.py all')
        raise SystemExit(2)

    # --> Get the target provid.
    target = sys.argv[1].strip()
    if target.lower() != "all":
        target = target.replace("_", " ")

    # --> Check the upstream files.
    for path in [mv.ORBITS_CSV, mv.OBSERVATIONS_PARQUET] + KERNEL_FILES:
        if not os.path.isfile(path):
            raise SystemExit(f"Missing upstream file: {path}")

    # --> Create the output directory and set up logging.
    os.makedirs(OUT_DIR, exist_ok=True)

    # --> Mirror all print() output to the log file and to the terminal.
    log = open(LOG_FILE, "a", encoding="utf-8")
    sys.stdout = mv.LogOutput(sys.__stdout__, log)
    print(f"Output: {OUT_DIR}")
    print(f"Log: {LOG_FILE}")
    t0 = time.perf_counter()

    # --> SPICE kernels and GM of the Sun.
    for path in KERNEL_FILES:
        spiceypy.furnsh(path)
    mu = spiceypy.bodvrd("SUN", "GM", 1)[1][0]

    # --> Orbits and their hashes.
    print(f"Reading {mv.ORBITS_CSV} ...")
    orbits = pd.read_csv(
        mv.ORBITS_CSV, dtype={"permid": "Int64"}, float_precision="round_trip"
    )
    orbits = orbits.dropna(subset=["provid_vro"]).drop_duplicates("provid_vro")
    orbits["orbit_hash"] = [orbit_hash(o) for o in orbits.to_dict("records")]
    orbits = orbits.set_index("provid_vro")
    n_empty = int(orbits["a"].isna().sum())
    print(f"  {len(orbits):,} objects, {n_empty} without elements")

    if target.lower() != "all" and target not in orbits.index:
        raise SystemExit(f"Provid not found in orbits.csv: {target}")

    # --> Objects whose orbit did not change keep their rows.
    old_hashes = {}
    if os.path.isfile(RESULTS_PARQUET):
        # --> Never mix rows with a file written with other columns.
        old_columns = pq.read_schema(RESULTS_PARQUET).names
        if old_columns != RESULTS_COLUMNS:
            raise SystemExit(
                f"{RESULTS_PARQUET} has columns {old_columns}, expected "
                f"{RESULTS_COLUMNS}. Delete it and run again."
            )
        old = pq.read_table(RESULTS_PARQUET, columns=["provid_vro", "orbit_hash"])
        old = old.group_by(["provid_vro", "orbit_hash"]).aggregate([]).to_pandas()
        old_hashes = dict(zip(old["provid_vro"], old["orbit_hash"]))
    unchanged = {
        p for p, h in orbits["orbit_hash"].items() if old_hashes.get(p) == h
    }
    if target.lower() == "all":
        todo = [p for p in orbits.index if p not in unchanged]
    else:
        todo = [] if target in unchanged else [target]
    print(f"  {len(todo):,} objects to compute")

    # --> Observations of the objects to compute.
    n_obs = n_no_provid = n_no_orbit = 0
    if todo:
        print(f"Reading {mv.OBSERVATIONS_PARQUET} ...")
        obs = pd.read_parquet(
            mv.OBSERVATIONS_PARQUET, columns=["id", "provid", "obstime", "mag"]
        )
        n_obs = len(obs)
        n_no_provid = int(obs["provid"].isna().sum())
        obs = obs.dropna(subset=["provid"])
        n_no_orbit = int((~obs["provid"].isin(orbits.index)).sum())
        obs = obs[obs["provid"].isin(todo)]
        obs = obs.sort_values(["provid", "obstime"]).reset_index(drop=True)
        obs["mag"] = pd.to_numeric(obs["mag"], errors="coerce")
        print(f"  {len(obs):,} observations to compute")

        # --> Observation times (UTC) as ET, and the Earth position, once per
        #     distinct time.
        times, inverse = np.unique(obs["obstime"].to_numpy(), return_inverse=True)
        ets = np.array(
            [
                spiceypy.str2et(f"{pd.Timestamp(t):%Y-%m-%d %H:%M:%S.%f} UTC")
                for t in times
            ]
        )
        earth = np.array(
            [spiceypy.spkpos("EARTH", et, "ECLIPJ2000", "NONE", "SUN")[0] for et in ets]
        )
        obs["obs_et"] = ets[inverse]
        r_earth = earth[inverse]
        print(f"  {len(times):,} distinct observation times")

        # --> Propagate every object to its observation times.
        geometry_cols = ["deltat", "r", "delta", "alpha", "H"]
        out = {col: np.full(len(obs), np.nan) for col in geometry_cols}
        provids = obs["provid"].to_numpy()
        starts = np.flatnonzero(np.r_[True, provids[1:] != provids[:-1]])
        stops = np.r_[starts[1:], len(obs)]
        obs_ets = obs["obs_et"].to_numpy()
        mags = obs["mag"].to_numpy()
        for k, (i, j) in enumerate(zip(starts, stops), start=1):
            res = mvroc_state_asteroid(
                orbits.loc[provids[i]], obs_ets[i:j], mags[i:j], r_earth[i:j], mu
            )
            for col, values in res.items():
                out[col][i:j] = values
            if k % 20000 == 0:
                elapsed = (time.perf_counter() - t0) / 60.0
                print(f"  {k:,}/{len(starts):,} objects ({elapsed:.1f} m)")

        # --> New rows (flag_updated = 1).
        info = orbits.loc[obs["provid"], ["provid_mpc", "permid", "orbit_hash"]]
        new = pd.DataFrame(
            {
                "id": obs["id"].to_numpy(),
                "provid_vro": obs["provid"].to_numpy(),
                "provid_mpc": info["provid_mpc"].to_numpy(),
                "permid": info["permid"].to_numpy(),
                "obstime": obs["obstime"].to_numpy(),
                "obs_et": obs_ets,
                **out,
                "mag": mags,
                "flag_updated": np.ones(len(obs), dtype=np.int8),
                "orbit_hash": info["orbit_hash"].to_numpy(),
            },
            columns=RESULTS_COLUMNS,
        )
        new_table = pa.Table.from_pandas(new, schema=SCHEMA, preserve_index=False)
        del obs, new, r_earth

        # --> Rewrite states.parquet: kept rows (flag_updated = 0) + new rows.
        #     With all, objects no longer in orbits.csv are dropped; with one
        #     provid, every other object is kept as it is.
        if target.lower() == "all":
            keep_set = unchanged
        else:
            keep_set = set(old_hashes) - set(todo)
        tables = []
        if os.path.isfile(RESULTS_PARQUET) and keep_set:
            kept = pq.read_table(RESULTS_PARQUET, schema=SCHEMA)
            mask = pc.is_in(
                kept["provid_vro"], value_set=pa.array(sorted(keep_set), pa.string())
            )
            kept = kept.filter(mask)
            flags = pa.array(np.zeros(kept.num_rows, dtype=np.int8))
            kept = kept.set_column(
                kept.schema.get_field_index("flag_updated"), "flag_updated", flags
            )
            tables.append(kept)
        tables.append(new_table)
        result = pa.concat_tables(tables)
        tmp = f"{RESULTS_PARQUET}.tmp"
        pq.write_table(result, tmp)
        os.replace(tmp, RESULTS_PARQUET)
        print(f"  Wrote {result.num_rows:,} rows to {RESULTS_PARQUET}")

    # --> Single object: summary of its rows.
    if target.lower() != "all":
        rows = pd.read_parquet(RESULTS_PARQUET, filters=[("provid_vro", "==", target)])
        print_state(target, rows, recomputed=bool(todo))
        print(f"Done in {(time.perf_counter() - t0) / 60.0:.1f} m -> {OUT_DIR}")
        return

    # --> Catalogue summary.
    hs = pq.read_table(RESULTS_PARQUET, columns=["H", "r", "flag_updated"])
    n_rows = hs.num_rows
    n_geom = n_rows - hs["r"].null_count - int(pc.sum(pc.is_nan(hs["r"])).as_py() or 0)
    n_h = int(pc.sum(pc.is_finite(hs["H"])).as_py() or 0)
    n_new = int(pc.sum(hs["flag_updated"]).as_py() or 0)
    print()
    print("====================")
    print("States summary")
    print("====================")
    print(f"  Objects in orbits.csv: {len(orbits):,}")
    print(f"    recomputed: {len(todo):,}")
    print(f"    unchanged (kept): {len(unchanged):,}")
    print(f"  Rows in states.parquet: {n_rows:,}")
    print(f"    flag_updated = 1: {n_new:,}")
    print(f"    with geometry: {n_geom:,}")
    print(f"    with H: {n_h:,}")
    if todo:
        print(f"  Detections read: {n_obs:,}")
        print(f"    without provid (no row): {n_no_provid:,}")
        print(f"    provid not in orbits.csv (no row): {n_no_orbit:,}")
    print(f"Done in {(time.perf_counter() - t0) / 60.0:.1f} m -> {OUT_DIR}")

if __name__ == "__main__":
    mvroc_states()
