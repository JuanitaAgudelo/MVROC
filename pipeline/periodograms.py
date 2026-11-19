#!/usr/bin/env python3
"""Democratic period search for MOI asteroids.

Usage (always from the pipeline/ directory):

    python periodograms.py "1981 US22"
    python periodograms.py 1981_US22
    python periodograms.py all

Upstream stages: crossmatches, orbits, states, mois. Reads photometry via
mv.load_moi_photometry() from those outputs. Writes figures and the results
CSV to outputs/periodograms/. Outputs: periodograms.csv, *_periodogram.webp
and periodograms.log (same text as the terminal). Rows in periodograms.csv
are always appended - delete the CSV yourself if you do not want duplicates.
"""

import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle, LombScargleMultiband

import mvroc as mv

# ---------------------------------------------------------------------------
# Period-search constants
# ---------------------------------------------------------------------------
NB = 3  # number of best-sampled bands that cast democratic votes
EXCLUDE_LX = ("Lu", "Ly")
PERIODS_H = (0.5 / 60, 24)  # photometric period search window (hours)
ALIAS_P_MIN_H = 2.0  # minimum period for campaign aliases
N_FREQ = 100_000  # number of frequencies for periodogram
OUT_DIR = "outputs/periodograms"
RESULTS_CSV = f"{OUT_DIR}/periodograms.csv"
LOG_FILE = f"{OUT_DIR}/periodograms.log"
RESULTS_COLUMNS = [
    "provid",
    "epoch_utc_ini",
    "epoch_utc_end",
    "obs_duration",
    "n_obs",
    "n_obs_u",
    "n_obs_g",
    "n_obs_r",
    "n_obs_i",
    "n_obs_z",
    "n_obs_y",
    "p_rot_dem_h",
    "votes_dem",
    "n_series_dem",
    "bands_dem",
    "periods_dem_h",
]

# ---------------------------------------------------------------------------
# Periodogram helpers
# ---------------------------------------------------------------------------
def frequency_grid(ts):
    """Photometric frequency grid (cycles/day) covering PERIODS_H.

    Parameters
    ----------
    ts : array-like
        Observation times in days. The time span densifies the grid for long
        campaigns so narrow peaks are still resolved.

    Returns
    -------
    ndarray
        Frequencies in cycles per day, from 24/PERIODS_H[1] to
        24/PERIODS_H[0].

    Examples
    --------
    >>> freqs = frequency_grid([0.0, 1.0, 2.0])
    >>> freqs[0], freqs[-1]
    (1.0, 2880.0)
    """
    f_min = 24.0 / PERIODS_H[1]
    f_max = 24.0 / PERIODS_H[0]
    n_freq = N_FREQ
    baseline = float(np.ptp(np.asarray(ts, dtype=float)))
    if baseline > 0:
        # Densify the grid so long campaigns still resolve narrow peaks.
        n_peak = int(np.ceil((f_max - f_min) * 5.0 * baseline)) + 1
        n_freq = max(n_freq, n_peak)
    return np.linspace(f_min, f_max, n_freq)

def rotation_at_max_power(ts, ys, dys, freq_grid, bands=None):
    """Rotation period at the highest Lomb-Scargle power above the cadence cut.

    Photometric frequency <-> rotation: P_rot [h] = 48 / f [cycles/day]
    (two peaks per rotation for a typical asteroid light curve).

    Parameters
    ----------
    ts : array-like
        Times in days.
    ys : array-like
        Magnitudes (usually absolute magnitude H).
    dys : array-like
        Magnitude uncertainties.
    freq_grid : array-like
        Photometric frequencies in cycles/day (see frequency_grid).
    bands : array-like, optional
        Band labels per point. If given, uses LombScargleMultiband.

    Returns
    -------
    dict
        freq (cycles/day), p_rot_h (hours), power, and
        cadence_cut_h (hours).

    Examples
    --------
    >>> ts = np.linspace(0, 5, 200)
    >>> ys = np.sin(2 * np.pi * ts / (3.0 / 24.0))  # 3 h photometric period
    >>> dys = np.full_like(ys, 0.01)
    >>> freqs = frequency_grid(ts)
    >>> out = rotation_at_max_power(ts, ys, dys, freqs)
    >>> out["p_rot_h"]
    6.0
    """
    dys = np.maximum(np.asarray(dys, dtype=float), 1e-4)
    if bands is None:
        power = LombScargle(ts, ys, dys).power(freq_grid)
    else:
        power = LombScargleMultiband(ts, ys, np.asarray(bands), dys).power(freq_grid)
    p_phot_min = mv.cadence_min_period(ts)
    allowed = (24.0 / freq_grid) >= p_phot_min
    if not np.any(allowed):
        allowed = np.ones(len(freq_grid), dtype=bool)
        p_phot_min = 0.0
    best_i = int(np.flatnonzero(allowed)[np.argmax(power[allowed])])
    freq = float(freq_grid[best_i])
    return {
        "freq": freq,
        "p_rot_h": 48.0 / freq,
        "power": float(power[best_i]),
        "cadence_cut_h": 2.0 * p_phot_min,
    }

def daily_alias_cycles(p1, p2, alias_c=0.1):
    """Return 1 or 2 if two rotation periods are daily aliases, else 0.

    Parameters
    ----------
    p1, p2 : float
        Rotation periods in hours.
    alias_c : float, optional
        Tolerance in cycles/day when matching a 1- or 2-cycle offset.

    Returns
    -------
    int
        1, 2, or 0.

    Examples
    --------
    >>> daily_alias_cycles(6.0, 6.0)
    0
    >>> # Photometric freqs differ by ~1 cycle/day
    >>> daily_alias_cycles(48.0 / 5.0, 48.0 / 6.0)
    1
    """
    separation = abs(48.0 / p1 - 48.0 / p2)
    for cycles in (1.0, 2.0):
        if abs(separation - cycles) <= alias_c:
            return int(cycles)
    return 0

def same_voted_period(p1, p2, frac=0.03, alias_c=0.1):
    """True when two rotation periods agree or are 1-2 cycle/day aliases.

    Parameters
    ----------
    p1, p2 : float
        Rotation periods in hours.
    frac : float, optional
        Relative tolerance for a direct period match (default 3%).
    alias_c : float, optional
        Tolerance passed to daily_alias_cycles.

    Returns
    -------
    bool

    Examples
    --------
    >>> same_voted_period(6.00, 6.10)
    True
    >>> same_voted_period(6.0, 10.0)
    False
    """
    if abs(p1 - p2) <= frac * max(p1, p2):
        return True
    return daily_alias_cycles(p1, p2, alias_c) > 0

def multiband_combos(bands):
    """Every pair of bands, then all of them when more than two are present.

    Parameters
    ----------
    bands : sequence of str
        Band labels, e.g. ["Lr", "Lg", "Li"].

    Returns
    -------
    list of list of str
        Pair combinations, plus the full set if len(bands) > 2.

    Examples
    --------
    >>> multiband_combos(["Lr", "Lg", "Li"])
    [['Lr', 'Lg'], ['Lr', 'Li'], ['Lg', 'Li'], ['Lr', 'Lg', 'Li']]
    """
    bands = [str(b) for b in bands]
    combos = [
        [bands[i], bands[j]]
        for i in range(len(bands))
        for j in range(i + 1, len(bands))
    ]
    if len(bands) > 2:
        combos.append(list(bands))
    return combos

def campaign_alias_k_range(span_days, p_min_h=ALIAS_P_MIN_H):
    """Integer daily shifts k for campaign-duration aliases.

    Parameters
    ----------
    span_days : float
        Campaign length in days.
    p_min_h : float, optional
        Shortest rotation period (hours) to mark as an alias.

    Returns
    -------
    range
        Integer k values, or an empty range if the span is invalid.

    Examples
    --------
    >>> list(campaign_alias_k_range(10.0))[:3]
    [...]
    """
    span_days = float(span_days)
    if not np.isfinite(span_days) or span_days <= 0 or p_min_h <= 0:
        return range(0)
    f_long = 48.0 / (2.0 * float(PERIODS_H[1]))
    f_short = 48.0 / float(p_min_h)
    k_min = int(np.ceil(f_long - 2.0 / span_days - 1e-12))
    k_max = int(np.ceil(f_short - 1.0 / span_days - 1e-12))
    if k_max < k_min:
        return range(0)
    return range(k_min, k_max + 1)

def draw_campaign_alias_vlines(ax, span_days, tol=0.1):
    """Shade 2-span campaign aliases on a period axis.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis whose x-scale is rotation period in hours.
    span_days : float
        Campaign length in days.
    tol : float, optional
        Half-width of each shaded band in cycles/day.

    Returns
    -------
    None

    Examples
    --------
    >>> fig, ax = plt.subplots()
    >>> draw_campaign_alias_vlines(ax, span_days=12.0)
    >>> plt.close(fig)
    """
    span = float(span_days)
    if not np.isfinite(span) or span <= 0 or tol <= 0:
        return
    k_values = list(campaign_alias_k_range(span))
    k_first = k_values[0] if k_values else None
    base = 2.0 / span
    for k in k_values:
        freq = base + float(k)
        f_lo, f_hi = freq - float(tol), freq + float(tol)
        if f_lo <= 0:
            continue
        ax.axvspan(
            48.0 / f_hi,
            48.0 / f_lo,
            color="C4",
            alpha=0.25,
            zorder=1,
            linewidth=0,
            label=(
                f"Campaign aliases (+/-{tol:.1f})" if k == k_first else None
            ),
        )

# ---------------------------------------------------------------------------
# Democratic periodogram plot
# ---------------------------------------------------------------------------

def mvroc_periodogram_plot(
    provid,
    ts,
    ys,
    dys,
    freq_grid,
    series,
    adopted,
    span_days,
    band_label,
):
    """Plot a democratic periodogram and save it as a WebP figure.

    Draws the Lomb-Scargle power of one band, shades campaign aliases, and
    marks every voted rotation period (the adopted one in a stronger style).

    Parameters
    ----------
    provid : str
        Asteroid name (title and output filename).
    ts, ys, dys : array-like
        Times (days), magnitudes and errors for the band shown in the curve.
    freq_grid : array-like
        Photometric frequencies in cycles/day.
    series : list of dict
        Voted solutions; each needs p_rot_h and label.
    adopted : float
        Winning rotation period in hours.
    span_days : float
        Campaign length in days (for alias shading).
    band_label : str
        Band name shown in the title, e.g. "Lr".

    Returns
    -------
    str
        Path of the saved WebP figure.

    Examples
    --------
    >>> path = mvroc_periodogram_plot(
    ...     "1981 US22", ts, ys, dys, freqs, series, 2.94, 12.0, "Lr"
    ... )
    >>> path.endswith(".webp")
    True
    """
    power = np.asarray(LombScargle(ts, ys, dys).power(freq_grid), dtype=float)
    p_phot_min = mv.cadence_min_period(ts)
    period_h = 48.0 / freq_grid  # photometric freq -> rotation period

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(period_h, power, color="C0", lw=1.0)
    draw_campaign_alias_vlines(ax, span_days)

    voted_merged = mv.merge_frequencies(
        (48.0 / float(item["p_rot_h"]), item["label"]) for item in series
    )
    for group in voted_merged:
        p_rot = 48.0 / float(group["freq"])
        is_adopted = abs(p_rot - adopted) <= 0.03 * max(p_rot, adopted)
        ax.axvline(
            p_rot,
            color="0.15" if is_adopted else "C2",
            ls=":" if is_adopted else "--",
            lw=1.6 if is_adopted else 1.2,
            zorder=3,
            label=(
                f"vote, {p_rot:.2f} h"
                if is_adopted
                else f"{', '.join(group['labels'])}, {p_rot:.2f} h"
            ),
        )

    x_min = 2.0 * p_phot_min if p_phot_min > 0 else 2.0 * PERIODS_H[0]
    ax.set_xscale("log")
    ax.set_xlim(x_min, 2.0 * PERIODS_H[1])
    ax.set_xlabel("Rotation period (h)")
    ax.set_ylabel("Lomb-Scargle power")
    ax.set_title(f"{provid} democratic periodogram ({band_label})")
    ax.legend(loc="upper left", fontsize=7, frameon=True)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()

    path = f"{OUT_DIR}/{mv.slug(provid)}_periodogram.webp"
    mv.save_figure(fig, path)
    print(f"  democratic periodogram: {path}")
    return path

# ---------------------------------------------------------------------------
# Democratic periodogram for one asteroid
# ---------------------------------------------------------------------------
def mvroc_periodogram_asteroid(df, provid, Nb=NB, exclude=EXCLUDE_LX, qplot=True):
    """Democratic period search for one asteroid; append a CSV result row.

    Takes the Nb best-sampled bands (excluding Lu/Ly by default). Each band
    casts a Lomb-Scargle vote; every pair and all bands together cast more
    votes. Periods within 3%, or 1-2 cycles/day apart, count as the same vote.
    A unanimous vote prefers the combined series; otherwise the best-sampled
    single band among the winners wins.

    Parameters
    ----------
    df : DataFrame
        Full MOI photometry from mvroc.load_moi_photometry.
    provid : str
        Asteroid name, e.g. "1981 US22".
    Nb : int, optional
        How many best-sampled bands cast votes.
    exclude : sequence of str, optional
        Bands to ignore (default Lu, Ly).
    qplot : bool, optional
        If True (default), write the democratic periodogram figure.

    Returns
    -------
    dict or None
        Result row written to periodograms.csv, or None on failure.

    Examples
    --------
    >>> df = mv.load_moi_photometry()
    >>> row = mvroc_periodogram_asteroid(df, "1981 US22")
    >>> row["p_rot_dem_h"]
    2.94...
    """
    print("=" * 20 + "\n" + provid + "\n" + "=" * 20)

    # --> Prepare the asteroid data.
    df_obj = mv.prepare_asteroid(df, provid)
    if df_obj.empty:
        print("  skip: no photometric observations in pipeline products")
        return None

    print(f"  N = {len(df_obj)} observations")
    
    # --> Check if there are any H values.
    if df_obj["H"].notna().sum() == 0:
        print(f"  skip: no H values in the state product for this object")
        return None

    # --> Keep only rows with H, drop Lu/Ly, count points per band, take top Nb.
    with_h = df_obj.dropna(subset=["H"])
    counts = with_h.loc[~with_h["band"].isin(exclude), "band"].value_counts()
    usable = counts[counts >= 5].head(int(Nb))
    if usable.empty:
        print("  skip: no band with at least 5 points for democratic periodogram")
        return None
    print(f"  {len(usable)} bands with at least 5 points ({[band_i for band_i in usable.index]})")

    # --> Join observations from the top bands into a single DataFrame.
    pieces = [
        df_obj[(df_obj["band"] == band_i) & df_obj["H"].notna()]
        for band_i in usable.index
    ]
    df_top = pd.concat(pieces, ignore_index=True)

    # --> Create the frequency grid from the top bands.
    freq_grid = frequency_grid(df_top["t_days"].to_numpy())

    # --> Compute the rotation period for each band.
    series = []
    for band_i, n_points in usable.items():
        df_i = df_obj[(df_obj["band"] == band_i) & df_obj["H"].notna()]
        found = rotation_at_max_power(
            df_i["t_days"].to_numpy(),
            df_i["H"].to_numpy(),
            df_i["rmsmag"].fillna(mv.MIN_MAG_ERROR).to_numpy(),
            freq_grid,
        )
        found["label"] = str(band_i)
        found["n"] = int(n_points)
        series.append(found)
        print(
            f"  {band_i}: N = {int(n_points)}, P_rot = {found['p_rot_h']:.2f} h, "
            f"power = {found['power']:.3f}"
        )

    # --> Compute extra votes from every pair and from all bands together.
    if len(usable) >= 2:
        bands = [str(b) for b in usable.index]
        for combo in multiband_combos(bands):
            label = "+".join(combo)
            df_i = df_obj[df_obj["band"].isin(combo) & df_obj["H"].notna()]
            found = rotation_at_max_power(
                df_i["t_days"].to_numpy(),
                df_i["H"].to_numpy(),
                df_i["rmsmag"].fillna(mv.MIN_MAG_ERROR).to_numpy(),
                freq_grid,
                bands=df_i["band"].to_numpy(),
            )
            found["label"] = label
            found["n"] = int(len(df_i))
            series.append(found)
            print(
                f"  {label}: N = {len(df_i)}, P_rot = {found['p_rot_h']:.2f} h, "
                f"power = {found['power']:.3f}"
            )

    # --> Union-Find: merge series whose periods agree (or are daily aliases).
    periods = np.array([item["p_rot_h"] for item in series], dtype=float)
    parent = list(range(len(periods)))

    # --> Define the find function for the Union-Find algorithm.
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    # --> Merge series whose periods agree (or are daily aliases).
    for i in range(len(periods)):
        for j in range(i + 1, len(periods)):
            if same_voted_period(periods[i], periods[j]):
                parent[find(i)] = find(j)

    # --> Find the largest vote cluster and determine the adopted period.
    groups = {}
    for i in range(len(periods)):
        groups.setdefault(find(i), []).append(i)

    # Largest vote cluster wins; prefer more bands when the vote is unanimous.
    winner_ids = max(groups.values(), key=lambda members: (len(members), -min(members)))
    single_ids = [i for i in winner_ids if "+" not in series[i]["label"]]
    combined_ids = [i for i in winner_ids if "+" in series[i]["label"]]
    unanimous = len(winner_ids) == len(series)

    # --> Choose the best period from the largest vote cluster.
    if unanimous and combined_ids:
        chosen = max(combined_ids, key=lambda i: series[i]["label"].count("+"))
    elif single_ids:
        chosen = max(single_ids, key=lambda i: series[i]["n"])
    else:
        chosen = winner_ids[0]

    # --> Store the adopted period and the labels of the winning series.
    adopted = float(series[chosen]["p_rot_h"])
    labels = [series[i]["label"] for i in winner_ids]
    periods_h = ",".join(f"{item['label']}={item['p_rot_h']:.4f}" for item in series)

    print(
        f"  democratic (Nb = {int(Nb)}): P_rot = {adopted:.4f} h from "
        f"{series[chosen]['label']}, {len(winner_ids)}/{len(series)} votes "
        f"({', '.join(labels)})"
    )

    # --> Plot the democratic periodogram if requested.
    if qplot:
        band_plot = str(usable.index[0])
        df_plot = df_obj[(df_obj["band"] == band_plot) & df_obj["H"].notna()]
        t_campaign = pd.to_numeric(df_obj["t_days"], errors="coerce").to_numpy(float)
        t_campaign = t_campaign[np.isfinite(t_campaign)]
        span_days = float(np.ptp(t_campaign)) if len(t_campaign) else 0.0
        mvroc_periodogram_plot(
            provid,
            df_plot["t_days"].to_numpy(float),
            df_plot["H"].to_numpy(float),
            df_plot["rmsmag"].fillna(mv.MIN_MAG_ERROR).to_numpy(float),
            freq_grid,
            series,
            adopted,
            span_days,
            band_plot,
        )

    # --> Compute the observation window and create the result row.
    epoch_ini, epoch_end, duration = mv.observation_window(df_obj)

    # --> Create the result row.
    row = {
        "provid": provid,
        "epoch_utc_ini": epoch_ini,
        "epoch_utc_end": epoch_end,
        "obs_duration": duration,
        "n_obs": int(len(df_obj)),
        **mv.band_counts(df_obj),
        "p_rot_dem_h": adopted,
        "votes_dem": int(len(winner_ids)),
        "n_series_dem": int(len(series)),
        "bands_dem": "+".join(str(b) for b in usable.index),
        "periods_dem_h": periods_h,
    }
    pd.DataFrame([row]).to_csv(RESULTS_CSV, mode="a", header=False, index=False)
    
    return row

# ---------------------------------------------------------------------------
# Periodogram driver
# ---------------------------------------------------------------------------
def mvroc_periodograms():
    """Run the democratic period search from the command line.

    Must be launched from the pipeline/ directory.

    Command-line
    ------------
    Exactly one positional argument:

    * provid - one asteroid name, e.g. "1981 US22" or 1981_US22
      (underscores are turned into spaces).
    * all - process every MOI asteroid (skipping those already in
      periodograms.csv).

    Examples
    --------
    python periodograms.py "1981 US22"
    python periodograms.py 1981_US22
    python periodograms.py all

    Inputs
    ------
    outputs/mois/mois.parquet
        MOI photometry (bands, times, errors).
    outputs/states/states.parquet
        Absolute magnitude H joined on observation id.
    outputs/crossmatches/crossmatches.csv
        Dual-provid aliases (optional merge).
    outputs/orbits/orbits.csv
        Orbital elements (upstream of states; not read directly here).

    Outputs
    -------
    outputs/periodograms/periodograms.csv
        One appended row per successful asteroid (created with header if
        missing). Asteroids already listed are skipped.
    outputs/periodograms/<provid>_periodogram.webp
        Democratic periodogram figure for each new asteroid.
    outputs/periodograms/periodograms.log
        Same text that is printed on screen (append mode).

    Returns
    -------
    None
        Prints progress and a final summary (succeeded / skipped / failed).
    """

    # --> Check command-line arguments.
    if len(sys.argv) != 2:
        print('Usage: python periodograms.py "<provid>" | python periodograms.py all')
        raise SystemExit(2)

    # --> Get the target provid
    target = sys.argv[1].strip()
    if target.lower() != "all":
        # Allow 1981_US22 as a convenient shell form of "1981 US22".
        target = target.replace("_", " ")

    # --> Create the output directory and set up logging.
    os.makedirs(OUT_DIR, exist_ok=True)

    # --> Mirror all print() output to the log file and to the terminal.
    log = open(LOG_FILE, "a", encoding="utf-8")
    sys.stdout = mv.LogOutput(sys.__stdout__, log)

    # --> Initialize the results CSV if it doesn't exist.
    if not os.path.isfile(RESULTS_CSV):
        pd.DataFrame(columns=RESULTS_COLUMNS).to_csv(RESULTS_CSV, index=False)
        print(f"Created {RESULTS_CSV}")
    print(f"Output: {OUT_DIR}")
    print(f"Log: {LOG_FILE}")

    # --> Load the MOI photometry and count the observations.
    df = mv.load_moi_photometry()
    selected = mv.observation_counts(df)

    # --> If not processing all asteroids, filter to the target.
    if target.lower() != "all":
        selected = selected[selected["provid"] == target].reset_index(drop=True)
        if selected.empty:
            raise SystemExit(f"Provid not found in MOI photometry: {target}")

    print(f"{len(selected)} asteroids to process (Nb = {NB})")
    t0 = time.perf_counter()
    succeeded = []
    failed = []
    skipped = []

    for i, provid in enumerate(selected["provid"], start=1):
        if mv.already_in_outputs(provid, RESULTS_CSV):
            print(f"{provid}: already in periodograms.csv, skip")
            skipped.append(provid)
            continue

        # --> Run the democratic periodogram for the target asteroid.
        row = mvroc_periodogram_asteroid(df, provid, Nb=NB)
        if row is not None:
            succeeded.append(provid)
        else:
            failed.append(provid)

        elapsed = time.perf_counter() - t0
        print(
            f"  elapsed: {elapsed / 60.0:.1f} m ({i}/{len(selected)}), "
            f"ok={len(succeeded)}, fail={len(failed)}"
        )

    print()
    print(f"Succeeded: {len(succeeded)}")
    print(f"Skipped: {len(skipped)}")
    print(f"Failed: {len(failed)}")
    if failed:
        print(failed[:20])
    print(f"Done in {(time.perf_counter() - t0) / 60.0:.1f} m -> {OUT_DIR}")

if __name__ == "__main__":
    mvroc_periodograms()
