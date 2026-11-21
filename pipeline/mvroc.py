"""Shared helpers for MVROC pipeline scripts.

Run every script from the pipeline/ directory so that all paths below
are relative (inputs/..., outputs/...).

mvroc.py only holds routines needed by more than one pipeline stage:
loaders, path helpers, and small photometry utilities. Science steps for a
given stage live in that stage's script (e.g. periodograms.py).
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Directories (cwd = pipeline/)
# ---------------------------------------------------------------------------
INPUTS = "inputs"
OUTPUTS = "outputs"

# ---------------------------------------------------------------------------
# Raw inputs (downloaded observation tables)
# ---------------------------------------------------------------------------
OBSERVATIONS_PARQUET = f"{INPUTS}/obs_sbn_X05_full.parquet"
MPCORB_DAT = f"{INPUTS}/mpcorb_extended.dat.gz"
KERNELS = f"{INPUTS}/kernels"  # SPICE kernels (copied from code/kernels)

# ---------------------------------------------------------------------------
# Stage outputs
# ---------------------------------------------------------------------------
CROSSMATCHES_CSV = f"{OUTPUTS}/crossmatches/crossmatches.csv"
ORBITS_CSV = f"{OUTPUTS}/orbits/orbits.csv"
STATES_PARQUET = f"{OUTPUTS}/states/states.parquet"
MOIS_PARQUET = f"{OUTPUTS}/mois/mois.parquet"
MOCS_CSV = f"{OUTPUTS}/mocs/mocs.csv"

# Aliases used by photometry loaders.
STATE_PARQUET = STATES_PARQUET
MOI_OBS_PARQUET = MOIS_PARQUET
VRO_MPC_CSV = CROSSMATCHES_CSV

# ---------------------------------------------------------------------------
# Shared photometry / figure constants
# ---------------------------------------------------------------------------
BAND = "Lr"
MIN_MAG_ERROR = 0.01  # floor / fill value for rmsmag
DPI = 150
WEBP_QUALITY = 80

# ---------------------------------------------------------------------------
# Logging / figure helpers
# ---------------------------------------------------------------------------
# Filesystem-safe asteroid name (spaces -> underscores).
slug = lambda provid: str(provid).replace(" ", "_")

class LogOutput:
    """Write the same text to several streams (screen + log file).

    Typical use in a pipeline main:

        log = open("outputs/periodograms/periodograms.log", "a", encoding="utf-8")
        sys.stdout = LogOutput(sys.__stdout__, log)

    Parameters
    ----------
    *streams
        File-like objects that implement write and flush.

    Examples
    --------
    >>> import sys
    >>> sys.stdout = LogOutput(sys.__stdout__, open("outputs/periodograms/periodograms.log", "a"))
    """

    def __init__(self, *streams):
        self.streams = streams

    def write(self, text):
        for stream in self.streams:
            stream.write(text)
            stream.flush()

    def flush(self):
        for stream in self.streams:
            stream.flush()

def save_figure(fig, path):
    """Write a WebP with an MVROC watermark on every main axes.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure to save (closed after writing).
    path : str
        Output path (.webp).

    Returns
    -------
    str
        Path written.

    Examples
    --------
    >>> fig, ax = plt.subplots()
    >>> ax.plot([0, 1], [0, 1])
    >>> save_figure(fig, "outputs/periodograms/_tmp_demo.webp")
    'outputs/periodograms/_tmp_demo.webp'
    """
    fig.canvas.draw()
    # Child axes (insets) and colorbars should not get the watermark.
    inset_ids = {
        id(child) for ax in fig.axes for child in getattr(ax, "child_axes", ())
    }
    for ax in fig.axes:
        if id(ax) in inset_ids or ax.get_label() == "<colorbar>":
            continue
        ax.text(
            1,
            1,
            "MVROC",
            rotation=270,
            ha="left",
            va="top",
            transform=ax.transAxes,
            color="pink",
            alpha=0.55,
            fontsize=12,
            zorder=100,
        )
    fig.savefig(path, dpi=DPI, pil_kwargs={"quality": WEBP_QUALITY, "method": 6})
    plt.close(fig)
    return path

# ---------------------------------------------------------------------------
# Common driver helpers
# ---------------------------------------------------------------------------
def cadence_min_period(ts):
    """Shortest resolvable photometric period: 2 x median delta t.

    Parameters
    ----------
    ts : array-like
        Observation times in days from an arbitrary zero point.

    Returns
    -------
    float
        Minimum photometric period in hours. 0.0 if fewer than two times.

    Examples
    --------
    >>> cadence_min_period([0.0, 0.01, 0.02, 0.03])
    0.48
    """
    t_sorted = np.sort(np.asarray(ts, dtype=float))
    dt_hours = np.diff(t_sorted) * 24.0
    dt_hours = dt_hours[dt_hours > 0]
    if len(dt_hours) == 0:
        return 0.0
    return 2.0 * float(np.median(dt_hours))

def merge_frequencies(items, tol=0.1):
    """Group (frequency, label) pairs that agree within tol cycles/day.

    Parameters
    ----------
    items : iterable of (float, str)
        Photometric frequencies (cycles/day) with a label each.
    tol : float, optional
        Absolute frequency tolerance in cycles/day.

    Returns
    -------
    list of dict
        Each dict has freq and labels (list of str).

    Examples
    --------
    >>> merge_frequencies([(5.0, "Lr"), (5.05, "Lg"), (7.0, "Li")])
    [{'freq': 5.0, 'labels': ['Lr', 'Lg']}, {'freq': 7.0, 'labels': ['Li']}]
    """
    groups = []
    for freq, label in items:
        freq = float(freq)
        for group in groups:
            if abs(group["freq"] - freq) <= tol:
                if label not in group["labels"]:
                    group["labels"].append(label)
                break
        else:
            groups.append({"freq": freq, "labels": [label]})
    return groups

def already_in_outputs(provid, results_csv):
    """True if provid already has a row in the stage results CSV.

    Parameters
    ----------
    provid : str
        Asteroid name, e.g. "1981 US22".
    results_csv : str
        Path to the stage CSV, e.g. outputs/foldings/foldings.csv.

    Returns
    -------
    bool

    Examples
    --------
    >>> already_in_outputs("1981 US22", "outputs/periodograms/periodograms.csv")
    True
    """
    if not os.path.isfile(results_csv):
        return False
    done = pd.read_csv(results_csv, usecols=["provid"])["provid"].astype(str)
    return str(provid) in set(done)

def provids_in_outputs(results_csv, column="provid"):
    """Set of names already written in a stage results CSV.

    Faster than already_in_outputs when a stage checks many objects: the
    CSV is read once.

    Parameters
    ----------
    results_csv : str
        Path to the stage CSV, e.g. outputs/orbits/orbits.csv.
    column : str, optional
        Name column, e.g. provid or provid_vro.

    Returns
    -------
    set of str
        Empty if the CSV does not exist.

    Examples
    --------
    >>> "2026 DO29" in provids_in_outputs("outputs/crossmatches/crossmatches.csv", "provid_vro")
    True
    """
    if not os.path.isfile(results_csv):
        return set()
    done = pd.read_csv(results_csv, usecols=[column])[column].astype(str)
    return set(done)

def check_results_header(results_csv, columns):
    """Exit if a stage CSV was written with other columns.

    Appending rows under a different header (e.g. an old or simulated file)
    would corrupt the CSV, so the driver stops and asks to delete it.

    Parameters
    ----------
    results_csv : str
        Path to the stage CSV.
    columns : list of str
        Expected columns (the driver RESULTS_COLUMNS).

    Returns
    -------
    None
        Raises SystemExit if the header differs.

    Examples
    --------
    >>> check_results_header("outputs/orbits/orbits.csv", ["provid_vro", "provid_mpc"])
    SystemExit: outputs/orbits/orbits.csv has columns [...], expected [...]. Delete it and run again.
    """
    header = list(pd.read_csv(results_csv, nrows=0).columns)
    if header != list(columns):
        raise SystemExit(
            f"{results_csv} has columns {header}, expected {list(columns)}. "
            "Delete it and run again."
        )

# ---------------------------------------------------------------------------
# Photometry loaders
# ---------------------------------------------------------------------------
def load_vro_mpc_alias_groups(path=None):
    """Groups of VRO provids that share one MPC designation.

    Parameters
    ----------
    path : str, optional
        CSV path (default outputs/crossmatches/crossmatches.csv).

    Returns
    -------
    list of dict
        Each dict has provid_mpc and provids (list of VRO names with
        at least two aliases).

    Examples
    --------
    >>> groups = load_vro_mpc_alias_groups()
    >>> groups[0]["provids"]
    ['...', '...']
    """
    path = path or VRO_MPC_CSV
    table = pd.read_csv(path)
    table["provid_vro"] = table["provid_vro"].astype(str).str.strip()
    table["provid_mpc"] = table["provid_mpc"].astype(str).str.strip()
    groups = []
    for provid_mpc, part in table.groupby("provid_mpc", sort=False):
        provids = [
            p
            for p in part["provid_vro"].tolist()
            if p and p not in ("nan", "None")
        ]
        provids = list(dict.fromkeys(provids))
        if len(provids) >= 2:
            groups.append({"provid_mpc": str(provid_mpc), "provids": provids})
    return groups

def merge_dual_provids(df):
    """Remap dual VRO names onto the alias with the most rows in df.

    Parameters
    ----------
    df : DataFrame
        Photometry with a provid column.

    Returns
    -------
    out : DataFrame
        Copy of df with remapped provid values.
    alias_map : dict
        Mapping from original provid -> canonical provid.

    Examples
    --------
    >>> df = pd.read_parquet("outputs/mois/mois.parquet", columns=["provid"])
    >>> out, alias_map = merge_dual_provids(df)
    """

    # --> Load the dual-provid alias groups from the CSV file.
    print(f"  Reading {VRO_MPC_CSV} ...")
    groups = load_vro_mpc_alias_groups()

    # --> Count the number of observations per provid.
    counts = df.groupby(df["provid"].astype(str).str.strip()).size()

    # --> Create a mapping from original provid to the canonical provid.
    alias_map = {}
    for group in groups:
        present = [p for p in group["provids"] if p in counts.index]
        if len(present) < 2:
            continue
        canonical = max(present, key=lambda p: int(counts.loc[p]))
        for provid in group["provids"]:
            alias_map[provid] = canonical

    # --> If no alias map was created, return the original DataFrame and an empty map.
    if not alias_map:
        return df.copy(), alias_map

    # --> Remap the provid column to the canonical provid.
    out = df.copy()
    out["provid"] = (
        out["provid"].astype(str).str.strip().map(lambda p: alias_map.get(p, p))
    )

    # --> Sort the data and return the remapped DataFrame and the alias map.
    out = out.sort_values(
        by=["provid", "band", "obstime"], ascending=[True, True, True]
    ).reset_index(drop=True)

    print(
        f"  Dual-provid merge: {len(set(alias_map.values()))} objects from "
        f"{VRO_MPC_CSV}"
    )
    return out, alias_map

def load_moi_photometry():
    """MOI light curves with absolute magnitude H from upstream stages.

    Joins:

    * Bands / errors / sky: outputs/mois/mois.parquet
    * Absolute magnitude: outputs/states/states.parquet (on observation id)
    * Dual-provid merge: outputs/crossmatches/crossmatches.csv

    Returns
    -------
    DataFrame
        Photometry with columns including provid, band, mag,
        rmsmag, H, obstime, night.

    Examples
    --------
    >>> df = load_moi_photometry()
    >>> df["provid"].nunique()
    8000
    """
    cols = [
        "id",
        "obstime",
        "permid",
        "provid",
        "ra",
        "dec",
        "band",
        "mag",
        "rmsmag",
    ]

    print("Loading MOI photometry ...")

    # --> Load the MOI photometry key columns from the parquet file.
    print(f"  Reading {MOI_OBS_PARQUET} ...")
    df = pd.read_parquet(MOI_OBS_PARQUET, columns=cols)
    print(f"    {len(df):,} rows, {df['provid'].nunique():,} provids")

    # --> Load the absolute magnitude H from the parquet file.
    print(f"  Reading {STATE_PARQUET} (H, provid_mpc) ...")
    state = pd.read_parquet(STATE_PARQUET, columns=["id", "H", "provid_mpc"])
    state["id"] = pd.to_numeric(state["id"], errors="coerce")
    df["id"] = pd.to_numeric(df["id"], errors="coerce")
    state = state.dropna(subset=["id"]).drop_duplicates("id", keep="first")
    df = df.merge(state, on="id", how="left")
    n_h = int(df["H"].notna().sum())
    print(f"    H matched for {n_h:,} / {len(df):,} rows")

    # --> Clean up the data.
    df["provid"] = df["provid"].astype(str).str.strip()
    df = df[df["provid"].ne("") & df["provid"].ne("nan")].copy()

    # --> Add the MPC designation column.
    df["mpc_designation"] = df["provid_mpc"].astype(str).str.strip()
    df.loc[
        df["mpc_designation"].isin(["", "nan", "None"]), "mpc_designation"
    ] = pd.NA

    # --> VRO sometimes stores a single letter ("r"); catalogue code uses "Lr".
    df["band"] = df["band"].map(
        lambda b: f"L{b}" if pd.notnull(b) and len(str(b)) == 1 else b
    )

    # --> Convert the columns to numeric.
    for col in ("ra", "dec", "mag", "rmsmag", "H"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["rmsmag"] = df["rmsmag"].where(
        df["rmsmag"].notna() & (df["rmsmag"] > 0), MIN_MAG_ERROR
    )

    # --> Convert the time column to datetime and add the night column (date without time).
    df["obstime"] = pd.to_datetime(df["obstime"])
    df["night"] = df["obstime"].dt.floor("D")

    # --> Merge dual-provid aliases.
    df, alias_map = merge_dual_provids(df)

    # --> Sort the data and add the alias map to the DataFrame attributes.
    df = df.sort_values(
        by=["provid", "band", "obstime"], ascending=[True, True, True]
    ).reset_index(drop=True)
    df.attrs["alias_map"] = alias_map

    print(f"  Photometry ready: {len(df):,} rows, {df['provid'].nunique():,} objects")

    return df

def observation_counts(df):
    """One row per provid with n_obs, richest first.

    Parameters
    ----------
    df : DataFrame
        Photometry with a provid column.

    Returns
    -------
    DataFrame
        Columns provid, n_obs.

    Examples
    --------
    >>> df = load_moi_photometry()
    >>> observation_counts(df).head(1)
       provid  n_obs
    0  ...      ...
    """
    return (
        df.groupby("provid", sort=False)
        .size()
        .reset_index(name="n_obs")
        .sort_values("n_obs", ascending=False)
        .reset_index(drop=True)
    )

def observation_window(df_ast):
    """Start, end (UTC strings) and duration in days for one asteroid.

    Parameters
    ----------
    df_ast : DataFrame
        Rows for a single asteroid with an obstime column.

    Returns
    -------
    ini : str
        First epoch as YYYY-MM-DD HH:MM:SS, or "".
    end : str
        Last epoch as YYYY-MM-DD HH:MM:SS, or "".
    duration : float
        Span in days, or nan if times are missing.

    Examples
    --------
    >>> df_obj = prepare_asteroid(load_moi_photometry(), "1981 US22")
    >>> observation_window(df_obj)
    ('2025-04-22 02:38:48', '2025-05-04 04:49:49', 12.09...)
    """
    t_min = df_ast["obstime"].min()
    t_max = df_ast["obstime"].max()
    ini = f"{t_min:%Y-%m-%d %H:%M:%S}" if pd.notna(t_min) else ""
    end = f"{t_max:%Y-%m-%d %H:%M:%S}" if pd.notna(t_max) else ""
    if pd.notna(t_min) and pd.notna(t_max):
        duration = float((t_max - t_min).total_seconds() / 86400.0)
    else:
        duration = float("nan")
    return ini, end, duration

def prepare_asteroid(df, provid):
    """Rows for one asteroid with t_days measured from the first epoch.

    Parameters
    ----------
    df : DataFrame
        Full MOI photometry.
    provid : str
        Asteroid name.

    Returns
    -------
    DataFrame
        Subset sorted by band and time, with a t_days column. Empty if the
        provid is missing.

    Examples
    --------
    >>> df = load_moi_photometry()
    >>> df_obj = prepare_asteroid(df, "1981 US22")
    >>> df_obj["t_days"].iloc[0]
    0.0
    """

    #--> Filter to the target asteroid.
    df_obj = df[df["provid"].astype(str) == str(provid)].copy()
    if df_obj.empty:
        return df_obj

    #--> Sort the data and add the t_days column.
    df_obj = df_obj.sort_values(["band", "obstime"]).reset_index(drop=True)
    t0 = df_obj["obstime"].iloc[0]

    #--> Create the t_days column as the time since the first observation.
    df_obj["t_days"] = (df_obj["obstime"] - t0).dt.total_seconds() / 86400.0

    return df_obj

def band_counts(df_ast):
    """Count points per LSST-like band (n_obs_u ... n_obs_y).

    Parameters
    ----------
    df_ast : DataFrame
        Rows for one asteroid with a band column (Lu/u, ...).

    Returns
    -------
    dict of str to int

    Examples
    --------
    >>> band_counts(pd.DataFrame({"band": ["Lr", "Lr", "Lg"]}))
    {'n_obs_u': 0, 'n_obs_g': 1, 'n_obs_r': 2, 'n_obs_i': 0, 'n_obs_z': 0, 'n_obs_y': 0}
    """
    return {
        f"n_obs_{b}": int((df_ast["band"].isin([f"L{b}", b])).sum())
        for b in ("u", "g", "r", "i", "z", "y")
    }
