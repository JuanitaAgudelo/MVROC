#!/usr/bin/env python
"""Repeat the rotation figures of search_spectra_rotation_orbit.ipynb.

Which asteroids and which figures to write is set in select_provid_all(),
select_provid_list(), select_provid_rotation(), select_provid_medellin(),
or select_provid_nobs().
The active one is the call at the start of main().

Writes, into figall/, as WebP at half the previous resolution:
  <provid>_lightcurve_allnights_H.webp
  <provid>_periodogram.webp
  <provid>_phase_allbands.webp
  <provid>_timemodel_H.webp
  <provid>_phase_model.webp
  <provid>_phase_spheres_fit.webp
  <provid>_phase_ellipsoid_fit.webp
  <provid>_lightcurve_allnights_Lk.webp
  <provid>_phase_model_Lk.webp
  <provid>_phase_shape_fit_Lk.webp
  <provid>_spectrum.webp
  <provid>_phase_model_all_periods.webp  (select_provid_rotation)
  period_compare_<period>_vs_lsm.webp     (select_provid_rotation)
  period_compare_<period>_vs_fourier.webp (select_provid_rotation)
  period_compare_grid.webp                (select_provid_rotation)
"""

import re
import shutil
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from astropy.timeseries import LombScargle, LombScargleMultiband
from scipy.optimize import least_squares
from scipy.signal import find_peaks

from util import ALL_BANDS, BANDS, band_colors

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figcomp"
_LOG_FILE = None
_PREV_STDOUT = None
ROTATION_CATALOG = DATA_DIR / "asteroids-with-rotation.csv"
MEDELLIN_CATALOG = DATA_DIR / "medellin-vro-list.csv"
DESIGNATIONS = DATA_DIR / "mpc_designations.csv"
OBS_FILE = DATA_DIR / "obs_sbn_X05_full.parquet"

BAND = "Lr"
PROP = "H"
PERIODS_H = (0.5 / 60, 24)  # photometric periods in hours; rotational plot spans 1 min to 48 h
# Lomb-Scargle search grid. 20_000 samples over this range miss the real peak.
N_FREQ = 100_000
MIN_N_OBS = 200
DPI = 150  # half of the previous 300 dpi
WEBP_QUALITY = 80
ASTEROID_PROVID = re.compile(r"^\d{4} [A-Z]{2}\d+[A-Z]?$")

# Per-asteroid analysis names used in main().
ANALYSIS = (
    "lightcurve",
    "periodogram",
    "periodogram_democratic",
    "phase",
    "phase_model",
    "shape_fit",
    "timemodel",
    "multiband_photometry",
    "shape_fit_multiband",
    "spectrum",
    "literature",
    "compare",
)
# Lx means "the band with the most observations", except these.
EXCLUDE_LX = ("Lu", "Ly")
FIGURES = ANALYSIS


def set_fig_dir(fig_dir):
    """Send this run's figures to fig_dir. A relative path is taken from ROOT."""
    global FIG_DIR
    fig_dir = Path(fig_dir)
    if not fig_dir.is_absolute():
        fig_dir = ROOT / fig_dir
    FIG_DIR = fig_dir
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    _start_log(FIG_DIR / "processing.log")
    print(f"Figures: {FIG_DIR}")
    print(f"Log: {FIG_DIR / 'processing.log'}")
    return FIG_DIR


class _Tee:
    """Write each print to the terminal and to the run log."""

    def __init__(self, stream, log_file):
        self.stream = stream
        self.log_file = log_file

    def write(self, data):
        self.stream.write(data)
        self.log_file.write(data)
        self.log_file.flush()

    def flush(self):
        self.stream.flush()
        self.log_file.flush()

    def isatty(self):
        return getattr(self.stream, "isatty", lambda: False)()


def _start_log(path):
    """Copy stdout into path. A later call replaces the previous log."""
    global _LOG_FILE, _PREV_STDOUT
    _stop_log()
    _LOG_FILE = open(path, "w", encoding="utf-8")
    _PREV_STDOUT = sys.stdout
    sys.stdout = _Tee(_PREV_STDOUT, _LOG_FILE)


def _stop_log():
    """Close the run log and restore the terminal stdout."""
    global _LOG_FILE, _PREV_STDOUT
    if _LOG_FILE is None:
        return
    sys.stdout = _PREV_STDOUT
    _LOG_FILE.close()
    _LOG_FILE = None
    _PREV_STDOUT = None


def slug(provid):
    return provid.replace(" ", "_")


def output_path(kind, provid):
    stem = {
        "lightcurve": "lightcurve_allnights_H",
        "periodogram": "periodogram",
        "phase": "phase_allbands",
        "timemodel": "timemodel_H",
        "phase_model": "phase_model",
        "spheres": "phase_spheres_fit",
        "ellipsoid": "phase_ellipsoid_fit",
        "multiband_lightcurve": "lightcurve_allnights_Lk",
        "multiband_phase_model": "phase_model_Lk",
        "spheres_multiband": "phase_spheres_fit_Lk",
        "ellipsoid_multiband": "phase_ellipsoid_fit_Lk",
        "shape_multiband": "phase_shape_fit_Lk",
        "spectrum": "spectrum",
        "literature": "phase_model_all_periods",
    }[kind]
    return FIG_DIR / f"{slug(provid)}_{stem}.webp"


def save_figure(fig, path):
    """Write a WebP at DPI. The path suffix is forced to .webp."""
    path = Path(path).with_suffix(".webp")
    fig.savefig(path, dpi=DPI, pil_kwargs={"quality": WEBP_QUALITY, "method": 6})
    plt.close(fig)
    return path


def format_elapsed(seconds):
    seconds = int(seconds)
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h {minutes:02d}m {secs:02d}s"
    return f"{minutes}m {secs:02d}s"


def check_analysis(analysis, where):
    unknown = set(analysis) - set(ANALYSIS)
    if unknown:
        names = ", ".join(sorted(unknown))
        allowed = ", ".join(ANALYSIS)
        raise SystemExit(f"Unknown analysis in {where}: {names}. Use: {allowed}")
    if "shape_fit_multiband" in analysis and "multiband_photometry" not in analysis:
        raise SystemExit(
            f"Error in {where}: 'shape_fit_multiband' cannot be executed without 'multiband_photometry'."
        )
    if "spectrum" in analysis and "multiband_photometry" not in analysis:
        raise SystemExit(
            f"Error in {where}: 'spectrum' cannot be executed without 'multiband_photometry'."
        )


check_figures = check_analysis


def load_df_set():
    """Photometry table, filtered as in the notebook (df_set)."""
    columns = ["obstime", "permid", "provid", "ra", "dec", "band", "mag", "rmsmag"]
    df_full = pd.read_parquet(OBS_FILE, columns=columns)
    print(f"Original number of observations: {len(df_full)}")

    df_set = df_full[
        ~(df_full["provid"].isna() & df_full["permid"].isna())
    ].reset_index(drop=True)
    print(f"Number of observations with provid or permid: {len(df_set)}")

    df_set = df_set[
        df_set["provid"].apply(lambda pid: bool(ASTEROID_PROVID.match(str(pid))))
    ].reset_index(drop=True)
    n_provid = df_set["provid"].dropna().nunique()
    print(f"Number of observations with asteroid provid: {len(df_set)}")
    print(f"Number of objects with provid: {n_provid}")

    df_set["band"] = df_set["band"].apply(
        lambda b: f"L{b}" if pd.notnull(b) and len(str(b)) == 1 else b
    )
    df_set = df_set.sort_values(
        by=["provid", "band", "obstime"],
        ascending=[True, True, True],
    ).reset_index(drop=True)
    df_set[["ra", "dec", "mag", "rmsmag"]] = df_set[
        ["ra", "dec", "mag", "rmsmag"]
    ].apply(pd.to_numeric, errors="coerce")
    df_set["obstime"] = pd.to_datetime(df_set["obstime"])
    df_set["night"] = df_set["obstime"].dt.floor("D")
    return df_set


def observation_counts(df_set):
    """df_obs from the notebook: number of observations per asteroid provid."""
    provids = df_set["provid"].dropna().unique()
    df_obs = (
        df_set[df_set["provid"].isin(provids)]
        .groupby("provid")
        .size()
        .reset_index(name="n_obs")
        .sort_values(by="n_obs", ascending=False)
    )
    return df_obs.reset_index(drop=True)


def load_rotation_catalog():
    """df_rot from the notebook: asteroids-with-rotation.csv."""
    df_rot = pd.read_csv(ROTATION_CATALOG)
    df_rot["provid"] = df_rot["provid"].dropna().astype(str).str.strip()
    df_rot = df_rot[df_rot["provid"] != ""].drop_duplicates("provid")
    return df_rot.reset_index(drop=True)


def load_mpc_designations():
    """provid -> mpc_designation, from data/mpc_designations.csv."""
    df = pd.read_csv(DESIGNATIONS)
    df["provid"] = df["provid"].astype(str).str.strip()
    return df.drop_duplicates("provid").set_index("provid")["mpc_designation"]


def load_geometry(designation):
    if designation is None or pd.isna(designation) or str(designation).strip() == "":
        return pd.DataFrame(columns=["obstime", "mag", "H"])
    table = pq.read_table(
        DATA_DIR / "obs_final_elements_geometry_H.parquet",
        columns=["obstime", "mag", "H"],
        filters=[("mpc_designation", "==", str(designation))],
    )
    df = table.to_pandas()
    df["obstime"] = pd.to_datetime(df["obstime"])
    return df


def prepare_object(df_set, provid, designations):
    df_obj = df_set[df_set["provid"] == provid].copy()
    if df_obj.empty:
        return df_obj
    df_obj = df_obj.sort_values(["band", "obstime"]).reset_index(drop=True)
    t0 = df_obj["obstime"].iloc[0]
    df_obj["t_days"] = (df_obj["obstime"] - t0).dt.total_seconds() / 86400.0
    designation = designations.get(provid, provid)
    df_spy = load_geometry(designation)
    df_obj = pd.merge(df_obj, df_spy, on="obstime", how="left", suffixes=("", "_spy"))
    return df_obj


def plot_lightcurve_h(df_obj, provid):
    path = output_path("lightcurve", provid)
    if path.exists():
        return path
    nights = df_obj.night.unique()
    n_nights = len(nights)
    fig, axs = plt.subplots(n_nights + 1, 1, figsize=(10, 3 * n_nights))
    if n_nights + 1 == 1:
        axs = [axs]

    ax = axs[0]
    for b in BANDS:
        sb = df_obj[df_obj["band"] == b]
        if len(sb) > 0:
            ax.errorbar(
                sb["obstime"],
                sb["H"],
                yerr=sb["rmsmag"],
                fmt="o",
                label=f"{b} (N={len(sb)})",
                color=band_colors.get(b, "k"),
                markersize=4,
                alpha=0.85,
            )
    ax.set_xlabel("Observation Time")
    ax.set_ylabel("H")
    ax.tick_params(axis="x", rotation=30)
    ax.set_title(f"Light Curve of Asteroid {provid} (all Nights)")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1, frameon=True)
    ax.grid(True, alpha=0.3)

    for i, night in enumerate(nights):
        night_str = pd.to_datetime(night).strftime("%Y-%m-%d")
        ax = axs[i + 1]
        t_night = df_obj[df_obj["night"] == night]
        t_night_0 = t_night["obstime"].min()
        hours = (t_night["obstime"] - t_night_0).dt.total_seconds() / 3600.0
        for b in ALL_BANDS:
            if b == "Lu":
                continue
            sel = t_night["band"] == b
            if not sel.any():
                continue
            ax.errorbar(
                hours[sel],
                t_night.loc[sel, "H"],
                yerr=t_night.loc[sel, "rmsmag"],
                fmt="o",
                label=f"{b} (N={sel.sum()})",
                color=band_colors.get(b, "k"),
                markersize=4,
                alpha=0.85,
            )
        ax.set_title(f"Night {night_str}")
        ax.set_xlabel("Time since start [h]")
        ax.set_ylabel("H")
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1, frameon=True)
        ax.tick_params(axis="x", rotation=30)

    fig.tight_layout()
    return save_figure(fig, output_path("lightcurve", provid))


def lightcurve_amplitudes(ls_band, t_band, y_band, frequency):
    """Peak-to-peak amplitudes of the periodogram model and of the photometry.

    amp_lsm is max(model) - min(model) of ls_band.model over one cycle at
    `frequency`. amp_photo is the width of y_band between the 5th and 95th
    percentiles.
    """
    t0 = float(np.min(t_band))
    phase = np.linspace(0, 1, 400, endpoint=False)
    y_model = ls_band.model(t0 + phase / frequency, frequency)
    amp_lsm = float(np.ptp(y_model))
    amp_photo = float(np.percentile(y_band, 95) - np.percentile(y_band, 5))
    return amp_lsm, amp_photo


PHOTO_BANDS = ("Lr", "Lg", "Li", "Lz", "Ly", "Lu")


def precise_band(df_obj, exclude=()):
    """Band with the most measurements, whatever their magnitude error.

    Lr wins a tie, so the reference band is kept when it is as well sampled
    as any other filter. Bands in `exclude` are not candidates.
    """
    candidates = [band for band in PHOTO_BANDS if band not in exclude]
    if BAND in candidates:
        candidates = [band for band in candidates if band != BAND] + [BAND]
    best_band, best_n = None, -1
    for band in candidates:
        n_points = int(((df_obj["band"] == band) & df_obj[PROP].notna()).sum())
        if n_points >= best_n:
            best_band, best_n = band, n_points
    return best_band, best_n


def resolve_band(df_obj, band, exclude_lx=EXCLUDE_LX):
    """Lr when no band is given. Lx is the band with the most observations."""
    requested = BAND if band is None or str(band).strip() == "" else str(band).strip()
    if requested != "Lx":
        return requested
    chosen, n_points = precise_band(df_obj, exclude=exclude_lx)
    print(f"  Lx -> {chosen} ({n_points} points)")
    return chosen


def cadence_min_period(t_days):
    """Shortest photometric period the sampling can resolve, in hours.

    Twice the median spacing between consecutive points. A photometric period
    shorter than that is an alias of the cadence.
    """
    t_sorted = np.sort(np.asarray(t_days, dtype=float))
    dt_hours = np.diff(t_sorted) * 24.0
    dt_hours = dt_hours[dt_hours > 0]
    if len(dt_hours) == 0:
        return 0.0
    return 2.0 * float(np.median(dt_hours))


def frequency_grid(t_band):
    """Photometric frequencies, in cycles per day, over PERIODS_H.

    Always at least N_FREQ samples. A long arc gets more, so there are five
    samples across a peak of width 1/baseline.
    """
    f_min = 24.0 / PERIODS_H[1]
    f_max = 24.0 / PERIODS_H[0]
    n_freq = N_FREQ
    baseline = float(np.ptp(np.asarray(t_band, dtype=float)))
    if baseline > 0:
        n_peak = int(np.ceil((f_max - f_min) * 5.0 * baseline)) + 1
        n_freq = max(n_freq, n_peak)
    return np.linspace(f_min, f_max, n_freq)


def periodogram(df_obj, provid, save=True, band=None, kind="periodogram"):
    if band is None:
        band = BAND
    df_band = df_obj[df_obj["band"] == band].dropna(subset=[PROP]).copy()
    if len(df_band) < 5:
        print(
            f"  skip periodogram for {provid}: {len(df_band)} {band} points with {PROP}"
        )
        return None

    t_band = df_band["t_days"].to_numpy()
    y_band = df_band[PROP].to_numpy()
    dy_band = df_band["rmsmag"].fillna(0.01).to_numpy()
    freq_grid = frequency_grid(t_band)
    ls_band = LombScargle(t_band, y_band, dy_band)
    power = ls_band.power(freq_grid)

    p_phot_min = cadence_min_period(t_band)
    allowed = (24.0 / freq_grid) >= p_phot_min
    if not np.any(allowed):
        allowed = np.ones(len(freq_grid), dtype=bool)
        p_phot_min = 0.0
    main_peak_power = np.max(power[allowed])
    best_freq = freq_grid[np.flatnonzero(allowed)[np.argmax(power[allowed])]]
    best_phot_period = 24.0 / best_freq
    best_rot_period = 2.0 * best_phot_period

    peaks, _ = find_peaks(power)
    high_power_peaks = np.intersect1d(
        peaks, np.where(power >= 0.8 * main_peak_power)[0]
    )
    high_power_peaks = high_power_peaks[
        (24.0 / freq_grid[high_power_peaks]) >= p_phot_min
    ]
    peak_data = []
    for idx in high_power_peaks:
        freq = freq_grid[idx]
        peak_data.append(
            {
                "frequency": freq,
                "rotational_period_h": 2.0 * 24.0 / freq,
                "relative_power": power[idx] / main_peak_power,
            }
        )
    if not peak_data:
        peak_data.append(
            {
                "frequency": best_freq,
                "rotational_period_h": best_rot_period,
                "relative_power": 1.0,
            }
        )
    df_peaks = (
        pd.DataFrame(peak_data)
        .sort_values("relative_power", ascending=False)
        .reset_index(drop=True)
    )
    all_rot_periods = df_peaks["rotational_period_h"].to_numpy()

    path = output_path(kind, provid)
    if save and not path.exists():
        fig, ax = plt.subplots(figsize=(10, 6))
        period_grid_hours = (24.0 / freq_grid) * 2.0
        ax.plot(period_grid_hours, power, color="b", lw=1.5)
        ax.axvline(
            best_rot_period,
            color="r",
            linestyle="--",
            lw=2,
            label=f"Peak Period = {best_rot_period:.2f} h",
        )
        if p_phot_min > 0:
            ax.axvline(
                2.0 * p_phot_min,
                color="0.45",
                linestyle="--",
                lw=1,
                label=f"Cadence cut = {2.0 * p_phot_min:.2f} h",
            )
        for period in all_rot_periods:
            ax.axvline(period, color="g", linestyle=":", alpha=1, lw=1)
        ax.set_xlabel("Rotational Period [hours]")
        ax.set_ylabel("Spectral Power")
        ax.set_title(
            f"Lomb-Scargle Periodogram for Asteroid {provid} ({band} Band, property={PROP})"
        )
        ax.set_xscale("log")
        x_min = 2.0 * p_phot_min if p_phot_min > 0 else 2.0 * PERIODS_H[0]
        ax.set_xlim(x_min, 2.0 * PERIODS_H[1])
        for level in (0.5, 0.7, 0.9):
            ax.axhline(
                level,
                color="0.45",
                linestyle=":",
                lw=1,
                label=f"power = {level:.1f}",
            )
        ax.legend(loc="upper right", frameon=True)
        ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        path = save_figure(fig, path)
    elif not save or not path.exists():
        path = None

    freq_sel = df_peaks.iloc[0]["frequency"]
    per_sel = df_peaks.iloc[0]["rotational_period_h"]
    amp_lsm, amp_photo = lightcurve_amplitudes(ls_band, t_band, y_band, freq_sel)
    print(
        f"  {provid}: P_rot = {per_sel:.2f} h ({band}, {PROP}), "
        f"amp_lsm = {amp_lsm:.3f}, amp_photo = {amp_photo:.3f}, "
        f"power = {main_peak_power:.3f}, cadence cut = {2.0 * p_phot_min:.2f} h"
    )
    return {
        "band": band,
        "df_band": df_band,
        "t_band": t_band,
        "y_band": y_band,
        "dy_band": dy_band,
        "ls_band": ls_band,
        "best_freq": best_freq,
        "best_phot_period": best_phot_period,
        "best_rot_period": best_rot_period,
        "freq_sel": freq_sel,
        "per_sel": per_sel,
        "period_source": f"Periodogram {band}",
        "amp_lsm": amp_lsm,
        "amp_photo": amp_photo,
        "power": float(main_peak_power),
        "df_peaks": df_peaks,
        "path": path,
    }


def _rotation_at_max_power(t, y, dy, freq_grid, bands=None):
    """Rotational period at the highest Lomb-Scargle power above the cadence cut."""
    dy = np.maximum(np.asarray(dy, dtype=float), 1e-4)
    if bands is None:
        power = LombScargle(t, y, dy).power(freq_grid)
    else:
        power = LombScargleMultiband(t, y, np.asarray(bands), dy).power(freq_grid)
    p_phot_min = cadence_min_period(t)
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


def _daily_alias_cycles(p1, p2, alias_c=0.1):
    """1 or 2 when the photometric frequencies differ by that many cycles per day."""
    separation = abs(48.0 / p1 - 48.0 / p2)
    for cycles in (1.0, 2.0):
        if abs(separation - cycles) <= alias_c:
            return int(cycles)
    return 0


def _same_voted_period(p1, p2, frac=0.03, alias_c=0.1):
    """True when two rotation periods agree, or differ by one or two cycles per day."""
    if abs(p1 - p2) <= frac * max(p1, p2):
        return True
    return _daily_alias_cycles(p1, p2, alias_c) > 0


def periodogram_democratic(df_obj, provid, Nb=3, exclude=EXCLUDE_LX):
    """Vote a rotation period from the Nb bands with the most points.

    Each of those bands votes with its own Lomb-Scargle maximum. Consecutive
    bands in that list, and all of them together, each cast one more vote.
    u and y are excluded. Periods within 3 percent, or one or two cycles per
    day apart, count as the same vote. A unanimous vote uses the period of
    all the bands together. Otherwise the reported period is the Lomb-Scargle
    period of the band with the most points among the series that cast the
    winning vote.
    """
    counts = (
        df_obj.dropna(subset=[PROP])
        .loc[~df_obj["band"].isin(exclude), "band"]
        .value_counts()
    )
    usable = counts[counts >= 5].head(int(Nb))
    if usable.empty:
        print(f"  skip democratic periodogram for {provid}: no band with 5 points")
        return None

    pieces = []
    for band_i, n_points in usable.items():
        df_i = df_obj[(df_obj["band"] == band_i) & df_obj[PROP].notna()]
        pieces.append(df_i)
    df_top = pd.concat(pieces, ignore_index=True)
    freq_grid = frequency_grid(df_top["t_days"].to_numpy())

    series = []
    for band_i, n_points in usable.items():
        df_i = df_obj[(df_obj["band"] == band_i) & df_obj[PROP].notna()]
        found = _rotation_at_max_power(
            df_i["t_days"].to_numpy(),
            df_i[PROP].to_numpy(),
            df_i["rmsmag"].fillna(0.01).to_numpy(),
            freq_grid,
        )
        found["label"] = str(band_i)
        found["n"] = int(n_points)
        series.append(found)
        print(
            f"  {band_i}: N = {int(n_points)}, P_rot = {found['p_rot_h']:.2f} h, "
            f"power = {found['power']:.3f}"
        )

    if len(usable) >= 2:
        bands = [str(band) for band in usable.index]
        combos = [bands[i : i + 2] for i in range(len(bands) - 1)]
        if len(bands) > 2:
            combos.append(bands)
        for combo in combos:
            label = "+".join(combo)
            df_i = df_obj[df_obj["band"].isin(combo) & df_obj[PROP].notna()]
            found = _rotation_at_max_power(
                df_i["t_days"].to_numpy(),
                df_i[PROP].to_numpy(),
                df_i["rmsmag"].fillna(0.01).to_numpy(),
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

    periods = np.array([item["p_rot_h"] for item in series], dtype=float)
    parent = list(range(len(periods)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(periods)):
        for j in range(i + 1, len(periods)):
            if _same_voted_period(periods[i], periods[j]):
                parent[find(i)] = find(j)

    groups = {}
    for i in range(len(periods)):
        groups.setdefault(find(i), []).append(i)
    winner_ids = max(groups.values(), key=lambda members: (len(members), -min(members)))
    single_ids = [i for i in winner_ids if "+" not in series[i]["label"]]
    combined_ids = [i for i in winner_ids if "+" in series[i]["label"]]
    unanimous = len(winner_ids) == len(series)
    if unanimous and combined_ids:
        chosen = max(combined_ids, key=lambda i: series[i]["label"].count("+"))
    elif single_ids:
        chosen = max(single_ids, key=lambda i: series[i]["n"])
    else:
        chosen = winner_ids[0]
    adopted = float(series[chosen]["p_rot_h"])
    labels = [series[i]["label"] for i in winner_ids]
    periods_h = ",".join(f"{item['label']}={item['p_rot_h']:.4f}" for item in series)
    print(
        f"  democratic (Nb = {int(Nb)}): P_rot = {adopted:.4f} h from {series[chosen]['label']}, "
        f"{len(winner_ids)}/{len(series)} votes ({', '.join(labels)})"
    )
    return {
        "p_rot_h": adopted,
        "votes": int(len(winner_ids)),
        "n_series": int(len(series)),
        "bands": "+".join(str(band) for band in usable.index),
        "series": ", ".join(labels),
        "source": series[chosen]["label"],
        "periods_h": periods_h,
    }


def plot_phase_allbands(df_obj, provid, fit, kind="phase"):
    path = output_path(kind, provid)
    if path.exists():
        return path
    t0 = fit["t_band"].min()
    freq_sel = fit["freq_sel"]
    per_sel = fit["per_sel"]
    source = fit.get("period_source") or f"Periodogram {fit.get('band', BAND)}"
    present = [
        b for b in ALL_BANDS if b not in ("u", "Lu") and (df_obj["band"] == b).any()
    ]
    m_max = {b: df_obj.loc[df_obj["band"] == b, PROP].min() for b in present}

    fig, ax = plt.subplots(figsize=(8, 5))
    for b in present:
        sb = df_obj[df_obj["band"] == b]
        phase_b = np.mod((sb["t_days"].to_numpy() - t0) * freq_sel / 2.0, 1.0)
        dy = sb["rmsmag"].fillna(0.01)
        marker = "o" if b.startswith("L") else "s"
        m = m_max[b]
        ax.errorbar(
            phase_b,
            sb[PROP] - m,
            yerr=dy,
            fmt=marker,
            ms=4,
            alpha=0.8,
            color=band_colors.get(b, "k"),
            label=f"{b} (N={len(sb)}, m_max={m:.2f})",
        )

    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Phase at the rotation period")
    ax.set_ylabel(r"$m - m_{\mathrm{max}}$")
    ax.set_title(
        f"Phase curve of {provid}, $P_{{rot}}$ = {per_sel:.2f} h ({source})"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)
    fig.tight_layout()
    return save_figure(fig, path)


def plot_timemodel_h(provid, fit):
    path = output_path("timemodel", provid)
    if path.exists():
        return path
    df_band = fit["df_band"]
    nights_band = df_band["night"].drop_duplicates().sort_values()
    n_nights_band = len(nights_band)
    fig, axs = plt.subplots(
        n_nights_band, 1, figsize=(10, 2.6 * n_nights_band), sharey=True
    )
    if n_nights_band == 1:
        axs = [axs]

    for ax, night in zip(axs, nights_band):
        sel = df_band["night"] == night
        t_night = df_band.loc[sel, "t_days"].to_numpy()
        y_night = df_band.loc[sel, PROP].to_numpy()
        dy_night = df_band.loc[sel, "rmsmag"].fillna(0.01).to_numpy()
        t_left, t_right = t_night.min(), t_night.max()
        pad = max(0.01, 0.05 * (t_right - t_left))
        t_grid = np.linspace(t_left - pad, t_right + pad, 400)
        y_grid = fit["ls_band"].model(t_grid, fit["best_freq"])
        hours = (t_night - t_left) * 24.0
        hours_grid = (t_grid - t_left) * 24.0
        ax.errorbar(
            hours,
            y_night,
            yerr=dy_night,
            fmt="o",
            ms=4,
            alpha=0.85,
            color="0.2",
            ecolor="0.6",
            label="Observed magnitude",
            zorder=3,
        )
        ax.plot(
            hours_grid,
            y_grid,
            color="C3",
            lw=1.8,
            label=f"Periodogram model, P = {fit['best_phot_period']:.2f} h",
            zorder=2,
        )
        ax.set_title(pd.to_datetime(night).strftime("%Y-%m-%d"))
        ax.set_ylabel("Magnitude")
        ax.grid(True, alpha=0.3)

    axs[0].invert_yaxis()
    axs[0].legend(loc="best", frameon=True)
    axs[-1].set_xlabel("Hours since the first observation of the night")
    fig.suptitle(
        f"{provid}: observed light curve and periodogram model ({fit.get('band', BAND)})"
    )
    fig.tight_layout()
    return save_figure(fig, output_path("timemodel", provid))


def literature_periods(df_rot, provid):
    """Fourier and LSM periods from df_rot, in the notebook's panel order."""
    row = df_rot.loc[df_rot["provid"] == provid]
    if row.empty:
        return {}
    periods = {}
    for column, name in (
        ("period_fourier_hr", "Literature (Fourier)"),
        ("period_lsm_hr", "Literature (Lomb-Scargle)"),
    ):
        value = pd.to_numeric(row.iloc[0][column], errors="coerce")
        if np.isfinite(value) and value > 0:
            periods[name] = float(value)
    return periods


def Lr_phase_curve(provid, fit, band=None, kind="phase_model"):
    """Folded photometry and the periodogram model, as in shape_fit.ipynb."""
    if band is None:
        band = fit.get("band", BAND)
    path = output_path(kind, provid)
    if path.exists():
        return path

    t_band = fit["t_band"]
    y_band = fit["y_band"]
    dy_band = fit["dy_band"]
    freq_sel = fit["freq_sel"]
    ls_band = fit["ls_band"]
    per_rot = fit["per_sel"]
    t0 = t_band.min()
    phase_rot = np.mod((t_band - t0) * freq_sel / 2.0, 1.0)
    phase_rot_grid = np.linspace(0, 1, 400)
    y_model_rot = ls_band.model(t0 + phase_rot_grid * 2.0 / freq_sel, freq_sel)

    fig, ax = plt.subplots(1, 1, figsize=(10, 4.5))
    curve_color = band_colors.get(band, "k")
    ax.errorbar(
        phase_rot,
        y_band,
        yerr=dy_band,
        fmt="o",
        ms=4,
        alpha=0.75,
        color=curve_color,
        ecolor="0.6",
        label="Magnitude",
    )
    ax.plot(
        phase_rot_grid,
        y_model_rot,
        color=curve_color,
        lw=2,
        label=rf"Periodogram model, $P_{{rot}}$ = {per_rot:.2f} h",
    )
    ax.hlines(
        [y_model_rot.min(), y_model_rot.max()],
        0,
        1,
        colors=curve_color,
        linestyles="--",
        alpha=0.5,
        label="Model min/max",
    )
    ax.hlines(
        [np.percentile(y_band, 1), np.percentile(y_band, 99)],
        0,
        1,
        colors="grey",
        linestyles=":",
        alpha=0.7,
        label="Photometry 1%/99%",
    )
    ax.set_xlabel("Phase at the rotation period ($2P$)")
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", frameon=True)
    ax.set_ylabel("Magnitude")
    ax.set_title(
        f"Folded light curve of {provid} ({PROP},{band}), "
        rf"$P_{{rot}}$ = {per_rot:.2f} h"
    )
    fig.tight_layout()
    return save_figure(fig, path)


def plot_phase_model_all_periods(provid, fit, lit_periods):
    """Folded light curve at our period and at the literature periods."""
    path = output_path("literature", provid)
    if path.exists():
        return path
    periods = {"VRO": fit["best_rot_period"], **lit_periods}
    t_band = fit["t_band"]
    y_band = fit["y_band"]
    dy_band = fit["dy_band"]
    ls_band = fit["ls_band"]
    t0 = t_band.min()
    band = fit.get("band", BAND)
    curve_color = band_colors.get(band, "k")

    fig, axs = plt.subplots(
        len(periods), 1, figsize=(10, 4.5 * len(periods)), squeeze=False
    )
    for ax, (name, per) in zip(axs.ravel(), periods.items()):
        freq_sel = 24.0 / (per / 2.0)
        phase_rot = np.mod((t_band - t0) * freq_sel / 2.0, 1.0)
        phase_rot_grid = np.linspace(0, 1, 400)
        y_model_rot = ls_band.model(t0 + phase_rot_grid * 2.0 / freq_sel, freq_sel)
        ax.errorbar(
            phase_rot,
            y_band,
            yerr=dy_band,
            fmt="o",
            ms=4,
            alpha=0.75,
            color=curve_color,
            ecolor="0.6",
            label="Magnitude",
        )
        ax.plot(
            phase_rot_grid,
            y_model_rot,
            color=curve_color,
            lw=2,
            label=rf"Periodogram model, $P_{{rot}}$ = {per:.2f} h",
        )
        ax.set_xlabel("Phase at the rotation period ($2P$)")
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", frameon=True)
        ax.set_ylabel("Magnitude")
        ax.set_title(
            f"Folded light curve of {provid} ({PROP},{band}), "
            rf"$P_{{rot}}$ = {per:.2f} h [{name}]"
        )
    fig.tight_layout()
    return save_figure(fig, output_path("literature", provid))


# (results column, filename stem, axis label). periods_h is parsed to its first value.
COMPARE_PERIODS = (
    ("p_periodogram_h", "periodogram", "Periodogram"),
    ("p_rot_dem_h", "democratic", "Democratic"),
    ("p_rot_sph_h", "spheres", "Touching spheres"),
    ("p_rot_ell_h", "ellipsoid", "Ellipsoid"),
    ("p_rot_multi_h", "multiband", "Multiband"),
    ("p_rot_sph_h_k", "spheres_multiband", "Touching spheres, multiband"),
    ("p_rot_ell_h_k", "ellipsoid_multiband", "Ellipsoid, multiband"),
)
LITERATURE_PERIODS = (
    ("period_lsm_hr", "lsm", "Literature, Lomb-Scargle"),
    ("period_fourier_hr", "fourier", "Literature, Fourier"),
)


def _comparison_table(df_rot):
    """Calculated periods joined to the literature periods."""
    csv_path = FIG_DIR / "asteroid_analysis_results.csv"
    if not csv_path.exists():
        print("  no asteroid_analysis_results.csv to compare")
        return None
    analysis = pd.read_csv(csv_path)
    if "periods_h" in analysis.columns:
        analysis["p_periodogram_h"] = pd.to_numeric(
            analysis["periods_h"].astype(str).str.split(",").str[0], errors="coerce"
        )
    lit_cols = ["provid"] + [name for name, _, _ in LITERATURE_PERIODS]
    literature = df_rot.loc[:, [c for c in lit_cols if c in df_rot.columns]]
    return analysis.merge(literature, on="provid", how="inner")


def _comparison_limits(table):
    columns = [name for name, _, _ in COMPARE_PERIODS if name in table.columns]
    columns += [name for name, _, _ in LITERATURE_PERIODS if name in table.columns]
    values = pd.to_numeric(table[columns].stack(), errors="coerce")
    values = values[np.isfinite(values) & (values > 0)]
    if values.empty:
        return 0.01, 48.0
    return float(values.min()) * 0.7, float(values.max()) * 1.4


def _comparison_pairs(table, ours, literature):
    if ours not in table.columns or literature not in table.columns:
        return table.iloc[0:0]
    frame = table[["provid", ours, literature]].copy()
    frame[ours] = pd.to_numeric(frame[ours], errors="coerce")
    frame[literature] = pd.to_numeric(frame[literature], errors="coerce")
    ok = np.isfinite(frame[ours]) & np.isfinite(frame[literature])
    ok &= (frame[ours] > 0) & (frame[literature] > 0)
    return frame.loc[ok]


def _draw_period_comparison(ax, pairs, x_col, y_col, x_label, y_label, limits, annotate):
    lo, hi = limits
    x = pairs[x_col].to_numpy()
    y = pairs[y_col].to_numpy()
    ax.scatter(x, y, s=28, color="#1f4e79", zorder=3)
    if annotate:
        for provid, xi, yi in zip(pairs["provid"], x, y):
            ax.annotate(
                provid,
                (xi, yi),
                textcoords="offset points",
                xytext=(4, 3),
                fontsize=6,
                alpha=0.85,
            )
    line = np.array([lo, hi])
    ax.plot(line, line, color="0.45", lw=1, label="1:1")
    ax.plot(line, 2.0 * line, color="C1", lw=1, ls="--", label="2:1")
    ax.plot(line, 0.5 * line, color="C2", lw=1, ls="--", label="1:2")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel(f"{x_label} [h]")
    ax.set_ylabel(f"{y_label} [h]")
    ax.set_title(f"{y_label} vs {x_label}")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", frameon=True, fontsize=8)


def plot_period_comparison(df_rot):
    """Each calculated period against the literature Lomb-Scargle and Fourier periods.

    Writes one figure per pair, with the 1:1, 2:1 and 1:2 lines, and a grid
    of all of them. Periods come from asteroid_analysis_results.csv.
    """
    table = _comparison_table(df_rot)
    if table is None or table.empty:
        print("  no asteroids in common with the literature catalogue")
        return []
    limits = _comparison_limits(table)
    present = [
        spec
        for spec in COMPARE_PERIODS
        if spec[0] in table.columns and table[spec[0]].notna().any()
    ]
    if not present:
        print("  no calculated periods to compare")
        return []

    paths = []
    for y_col, y_stem, y_label in present:
        for x_col, x_stem, x_label in LITERATURE_PERIODS:
            pairs = _comparison_pairs(table, y_col, x_col)
            if pairs.empty:
                print(f"  no pairs for {y_label} vs {x_label}")
                continue
            fig, ax = plt.subplots(figsize=(8, 8))
            _draw_period_comparison(
                ax, pairs, x_col, y_col, x_label, y_label, limits, annotate=True
            )
            fig.tight_layout()
            path = save_figure(
                fig, FIG_DIR / f"period_compare_{y_stem}_vs_{x_stem}.webp"
            )
            paths.append(path)
            print(f"  {path.name} ({len(pairs)} asteroids)")

    n_rows, n_cols = len(present), len(LITERATURE_PERIODS)
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(6.2 * n_cols, 5.4 * n_rows), squeeze=False
    )
    for i, (y_col, _, y_label) in enumerate(present):
        for j, (x_col, _, x_label) in enumerate(LITERATURE_PERIODS):
            ax = axes[i, j]
            pairs = _comparison_pairs(table, y_col, x_col)
            if pairs.empty:
                ax.set_axis_off()
                continue
            _draw_period_comparison(
                ax, pairs, x_col, y_col, x_label, y_label, limits, annotate=False
            )
    fig.tight_layout()
    path = save_figure(fig, FIG_DIR / "period_compare_grid.webp")
    paths.append(path)
    print(f"  {path.name}")
    return paths


def touching_spheres_area(theta, q, phi):
    """Projected area of two touching spheres (r1 = 1, r2 = q)."""
    r1, r2 = 1.0, q
    d = (r1 + r2) * np.sqrt(np.cos(theta) ** 2 + np.sin(theta) ** 2 * np.cos(phi) ** 2)
    d = np.maximum(d, 1e-9)
    xi1 = np.arccos(np.clip((d**2 + r1**2 - r2**2) / (2 * d * r1), -1, 1))
    xi2 = np.arccos(np.clip((d**2 + r2**2 - r1**2) / (2 * d * r2), -1, 1))
    return r1**2 * (np.pi - xi1 + np.sin(xi1) * np.cos(xi1)) + r2**2 * (
        np.pi - xi2 + np.sin(xi2) * np.cos(xi2)
    )


def touching_spheres_mag(phase, m_sep, q, cos2phi, theta0):
    """Magnitude versus rotational phase. m_sep is the fully separated magnitude."""
    phi = np.arccos(np.sqrt(np.clip(cos2phi, 0, 1)))
    area = touching_spheres_area(2 * np.pi * phase + theta0, q, phi)
    return m_sep - 2.5 * np.log10(area / (np.pi * (1.0 + q**2)))


def ellipsoid_area(theta, b, c, phi, a=1.0):
    """Projected area of a triaxial ellipsoid rotating about c."""
    sphi, cphi = np.sin(phi), np.cos(phi)
    return np.pi * np.sqrt(
        (a * c * np.cos(theta) * sphi) ** 2
        + (b * c * np.sin(theta) * sphi) ** 2
        + (a * b * cphi) ** 2
    )


def ellipsoid_mag(phase, m_ab, b, c_over_b, cos2phi, theta0):
    """Magnitude versus rotational phase. m_ab is the a-b face, with a = 1."""
    phi = np.arccos(np.sqrt(np.clip(cos2phi, 0, 1)))
    area = ellipsoid_area(2 * np.pi * phase + theta0, b, c_over_b * b, phi)
    return m_ab - 2.5 * np.log10(area / (np.pi * b))


def _rotational_phase(t, t0, freq):
    return np.mod((t - t0) * freq / 2.0, 1.0)


def _best_least_squares(residuals, lower, upper, starts):
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    best = None
    for p0 in starts:
        p0 = np.clip(np.asarray(p0, dtype=float), lower + 1e-10, upper - 1e-10)
        sol = least_squares(residuals, p0, bounds=(lower, upper))
        if best is None or sol.cost < best.cost:
            best = sol
    return best


def _fit_scale(cost, jac, n_data, n_param):
    ndof = n_data - n_param
    chi_r = 2 * cost / ndof if ndof > 0 else np.nan
    cov = np.linalg.pinv(jac.T @ jac) * chi_r
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    return chi_r, err


def _fit_touching_spheres(t, y, dy, freq0, t0):
    def residuals(p):
        m_sep, q, cos2phi, theta0, freq = p
        model = touching_spheres_mag(
            _rotational_phase(t, t0, freq), m_sep, q, cos2phi, theta0
        )
        return (y - model) / dy

    lower = [y.min() - 2, 0.05, 0.0, -2 * np.pi, freq0 * 0.98]
    upper = [y.max() + 2, 20.0, 1.0, 4 * np.pi, freq0 * 1.02]
    rng = np.random.default_rng(42)
    starts = [
        [
            np.median(y) + rng.normal(0, 0.3),
            rng.uniform(0.2, 3.0),
            rng.uniform(0, 1),
            rng.uniform(0, 2 * np.pi),
            freq0 * (1 + rng.normal(0, 5e-4)),
        ]
        for _ in range(100)
    ]
    best = _best_least_squares(residuals, lower, upper, starts)
    m_sep, q_fit, cos2phi, theta0, freq = best.x
    chi_r, err = _fit_scale(best.cost, best.jac, len(t), len(best.x))
    q_can = min(q_fit, 1 / q_fit)
    q_err = err[1] * (1 if q_fit <= 1 else 1 / q_fit**2)
    phi = np.arccos(np.sqrt(np.clip(cos2phi, 0, 1)))
    phi_err = err[2] / max(2 * np.sqrt(max(cos2phi, 0)) * np.sin(phi), 1e-3)
    phi_is_limit = cos2phi < 2 * err[2]
    ph_grid = np.linspace(0, 1, 2000)
    m_grid = touching_spheres_mag(ph_grid, m_sep, q_fit, cos2phi, theta0)
    m_bright, m_faint = float(m_grid.min()), float(m_grid.max())
    return {
        "m_sep": float(m_sep),
        "q_fit": float(q_fit),
        "q_sph": float(q_can),
        "q_sph_err": float(q_err),
        "cos2phi": float(cos2phi),
        "theta0": float(theta0),
        "freq": float(freq),
        "chi_r_sph": float(chi_r),
        "phi_sph_deg": float(np.rad2deg(phi)),
        "phi_sph_err_deg": np.nan if phi_is_limit else float(np.rad2deg(phi_err)),
        "amp_sph": m_faint - m_bright,
        "m_mid_sph": 0.5 * (m_bright + m_faint),
        "m_bright": m_bright,
        "m_faint": m_faint,
        "p_rot_sph_h": float(48 / freq),
        "ph_grid": ph_grid,
        "m_grid": m_grid,
        "rms": float(np.std(best.fun * dy)),
    }


def _fit_ellipsoid(t, y, dy, freq0, t0):
    def residuals(p):
        m_ab, b, theta0, freq = p
        model = ellipsoid_mag(_rotational_phase(t, t0, freq), m_ab, b, 1.0, 0.0, theta0)
        return (y - model) / dy

    lower = [y.min() - 2, 0.05, -2 * np.pi, freq0 * 0.98]
    upper = [y.max() + 2, 1.0, 4 * np.pi, freq0 * 1.02]
    rng = np.random.default_rng(7)
    starts = [
        [
            np.median(y) + rng.normal(0, 0.2),
            rng.uniform(0.3, 0.95),
            rng.uniform(0, 2 * np.pi),
            freq0 * (1 + rng.normal(0, 5e-4)),
        ]
        for _ in range(60)
    ]
    best = _best_least_squares(residuals, lower, upper, starts)
    m_ab, b, theta0, freq = best.x
    chi_r, err = _fit_scale(best.cost, best.jac, len(t), len(best.x))
    ph_grid = np.linspace(0, 1, 2000)
    m_grid = ellipsoid_mag(ph_grid, m_ab, b, 1.0, 0.0, theta0)
    m_bright, m_faint = float(m_grid.min()), float(m_grid.max())
    return {
        "m_ab": float(m_ab),
        "b_over_a": float(b),
        "b_over_a_err": float(err[1]),
        "theta0": float(theta0),
        "freq": float(freq),
        "chi_r_ell": float(chi_r),
        "amp_ell": m_faint - m_bright,
        "m_mid_ell": 0.5 * (m_bright + m_faint),
        "p_rot_ell_h": float(48 / freq),
        "ph_grid": ph_grid,
        "m_grid": m_grid,
        "rms": float(np.std(best.fun * dy)),
    }


def _plot_spheres_fit(provid, t, y, dy, t0, spheres, band=BAND, path=None):
    if path is None:
        path = output_path("spheres_multiband" if band == "Lk" else "spheres", provid)
    if path.exists():
        return path
    curve_color = band_colors.get(band, "k")
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(10, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    ph_data = _rotational_phase(t, t0, spheres["freq"])
    ax.errorbar(
        ph_data,
        y,
        yerr=dy,
        fmt="o",
        ms=4,
        alpha=0.6,
        color=curve_color,
        ecolor="0.7",
        label="Magnitude",
    )
    ax.plot(
        spheres["ph_grid"],
        spheres["m_grid"],
        "k",
        lw=2.5,
        label=rf"Touching spheres: $q$={spheres['q_sph']:.2f}, "
        rf"$\phi$={spheres['phi_sph_deg']:.0f}$^\circ$, $\Delta m$={spheres['amp_sph']:.2f}",
    )
    ax.axhline(
        spheres["m_mid_sph"],
        color="C3",
        ls="--",
        lw=1.2,
        label=f"Midpoint = {spheres['m_mid_sph']:.3f}",
    )
    ax.axhline(spheres["m_bright"], color="0.5", ls=":", lw=1)
    ax.axhline(spheres["m_faint"], color="0.5", ls=":", lw=1)
    ax.set_ylabel("Magnitude")
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", frameon=True)
    ax.set_title(
        f"Touching-spheres fit to {provid} ({PROP},{band}), "
        rf"$P_{{rot}}$ = {spheres['p_rot_sph_h']:.3f} h, $\chi^2_\nu$ = {spheres['chi_r_sph']:.1f}"
    )
    _draw_sphere_inset(ax, spheres["q_sph"])
    model = touching_spheres_mag(
        ph_data,
        spheres["m_sep"],
        spheres["q_fit"],
        spheres["cos2phi"],
        spheres["theta0"],
    )
    axr.errorbar(
        ph_data,
        y - model,
        yerr=dy,
        fmt="o",
        ms=3,
        alpha=0.6,
        color=curve_color,
        ecolor="0.7",
    )
    axr.axhline(0, color="k", lw=1)
    axr.set_ylabel("Data - model")
    axr.set_xlabel("Phase at the rotation period ($2P$)")
    axr.invert_yaxis()
    axr.set_xlim(0, 1)
    axr.grid(True, alpha=0.3)
    fig.tight_layout()
    return save_figure(fig, path)


def _draw_sphere_inset(ax, q_can, box=None, radii=True):
    """Oblique view of the two fitted lobes, with radii in the ratio q."""
    alpha, beta = np.deg2rad(18), np.deg2rad(12)
    r1_s, r2_s = 1.0, q_can
    view = np.array([np.sin(alpha), np.cos(alpha) * np.cos(beta), np.sin(beta)])
    view /= np.linalg.norm(view)
    right = np.cross(np.array([0.0, 0.0, 1.0]), -view)
    right /= np.linalg.norm(right)
    up = np.cross(-view, right)

    def sphere_face(color, n=200):
        g = np.linspace(-1, 1, n)
        u, v = np.meshgrid(g, g)
        rho2 = u**2 + v**2
        nz = np.sqrt(np.clip(1 - rho2, 0, None))
        light = np.array([-0.30, 0.48, 0.82])
        light /= np.linalg.norm(light)
        shade = np.clip(u * light[0] + v * light[1] + nz * light[2], 0, 1)
        rgb = np.array(plt.matplotlib.colors.to_rgb(color))
        img = np.zeros((n, n, 4))
        intensity = (0.20 + 0.80 * shade) * (0.45 + 0.55 * nz)
        img[..., :3] = np.clip(
            rgb * intensity[..., None] + 0.40 * (shade**24)[..., None], 0, 1
        )
        img[..., 3] = np.clip((1 - rho2) * (n / 3.5), 0, 1)
        return img

    spheres = [
        {"center": np.array([-r2_s, 0.0, 0.0]), "radius": r1_s, "label": r"$r_1$"},
        {"center": np.array([r1_s, 0.0, 0.0]), "radius": r2_s, "label": r"$r_2$"},
    ]
    spheres.sort(key=lambda s: np.dot(s["center"], view))
    if box is None:
        box = [0.705, 0.025, 0.28, 0.345]
    axins = ax.inset_axes(box)
    axins.set_xticks([])
    axins.set_yticks([])
    for spine in axins.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("0.45")
    axins.set_facecolor("white")
    placed = []
    for s in spheres:
        p = np.array([np.dot(s["center"], right), np.dot(s["center"], up)])
        r = s["radius"]
        axins.imshow(
            sphere_face("#8a8a8a"),
            origin="lower",
            aspect="equal",
            interpolation="bilinear",
            extent=[p[0] - r, p[0] + r, p[1] - r, p[1] + r],
            zorder=2,
        )
        axins.add_patch(
            plt.Circle((p[0], p[1]), r, fill=False, ec="0.12", lw=0.6, zorder=3)
        )
        placed.append((p, r, s["label"]))
    y_drop = 0.05 if radii else 0.08
    y_lab = min(p[1] - r for p, r, _ in placed) - y_drop
    if radii:
        for p, r, lab in placed:
            axins.text(p[0], y_lab, lab, ha="center", va="top", fontsize=8, color="0.1")
        axins.text(
            np.mean([p[0] for p, _, _ in placed]),
            y_lab - 0.42,
            rf"$q = r_2/r_1 = {q_can:.2f}$",
            ha="center",
            va="top",
            fontsize=7.5,
            color="0.15",
        )
        y_floor = y_lab - 0.85
        x_pad = 0.28
    else:
        axins.text(
            np.mean([p[0] for p, _, _ in placed]),
            y_lab - 0.28,
            rf"$q = {q_can:.2f}$",
            ha="center",
            va="top",
            fontsize=7,
            color="0.15",
        )
        y_floor = y_lab - 0.7
        x_pad = 0.22
    xs = np.concatenate([[p[0] - r, p[0] + r] for p, r, _ in placed])
    ys = np.concatenate([[p[1] - r, p[1] + r] for p, r, _ in placed])
    axins.set_xlim(xs.min() - x_pad, xs.max() + x_pad)
    axins.set_ylim(y_floor, ys.max() + (0.22 if radii else 0.18))
    axins.set_aspect("equal", adjustable="box")


def _draw_ellipsoid_inset(ax, b_over_a, box=None):
    """Oblique Lambert view of the equator-on ellipsoid, labelled by b/a."""
    a_ax, b_ax, c_ax = 1.0, float(b_over_a), float(b_over_a)
    alpha, beta = np.deg2rad(18), np.deg2rad(12)
    view = np.array([np.sin(alpha), np.cos(alpha) * np.cos(beta), np.sin(beta)])
    view /= np.linalg.norm(view)
    right = np.cross(np.array([0.0, 0.0, 1.0]), -view)
    right /= np.linalg.norm(right)
    up = np.cross(-view, right)
    n = 200
    span = 1.2 * max(a_ax, b_ax, c_ax)
    g = np.linspace(-span, span, n)
    uu, vv = np.meshgrid(g, g)
    origin = uu[..., None] * right + vv[..., None] * up + 4.0 * view
    direction = np.broadcast_to(-view, origin.shape)
    qdiag = np.array([1 / a_ax**2, 1 / b_ax**2, 1 / c_ax**2])
    quad_a = np.sum(direction**2 * qdiag, axis=-1)
    quad_b = 2 * np.sum(origin * direction * qdiag, axis=-1)
    quad_c = np.sum(origin**2 * qdiag, axis=-1) - 1.0
    disc = quad_b**2 - 4 * quad_a * quad_c
    hit = disc > 0
    sqrt_disc = np.sqrt(np.clip(disc, 0, None))
    travel = np.where(hit, (-quad_b - sqrt_disc) / (2 * quad_a), np.nan)
    point = origin + travel[..., None] * direction
    normal = point * qdiag
    normal /= np.linalg.norm(normal, axis=-1, keepdims=True)
    light = np.array([-0.30, 0.48, 0.82])
    light /= np.linalg.norm(light)
    shade = np.clip(np.sum(normal * light, axis=-1), 0, 1)
    facing = np.clip(np.sum(normal * (-view), axis=-1), 0, 1)
    rgb = np.array(plt.matplotlib.colors.to_rgb("#8a8a8a"))
    image = np.zeros((n, n, 4))
    intensity = (0.20 + 0.80 * shade) * (0.45 + 0.55 * facing)
    image[..., :3] = np.clip(
        rgb * intensity[..., None] + 0.40 * (shade**24)[..., None], 0, 1
    )
    image[..., 3] = np.where(hit, 1.0, 0.0)
    camera = np.column_stack([right, up])
    quadric = camera.T @ np.diag([a_ax**2, b_ax**2, c_ax**2]) @ camera
    evals, evecs = np.linalg.eigh(quadric)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    ell_w, ell_h = 2 * np.sqrt(np.clip(evals, 0, None))
    ell_ang = np.degrees(np.arctan2(evecs[1, 0], evecs[0, 0]))
    if box is None:
        box = [0.80, 0.012, 0.195, 0.255]
    axins = ax.inset_axes(box)
    axins.set_xticks([])
    axins.set_yticks([])
    for spine in axins.spines.values():
        spine.set_linewidth(0.6)
        spine.set_color("0.45")
    axins.set_facecolor("white")
    axins.imshow(
        image,
        origin="lower",
        aspect="equal",
        interpolation="bilinear",
        extent=[g.min(), g.max(), g.min(), g.max()],
        zorder=2,
    )
    axins.add_patch(
        plt.matplotlib.patches.Ellipse(
            (0, 0), ell_w, ell_h, angle=ell_ang, fill=False, ec="0.12", lw=0.6, zorder=3
        )
    )
    axins.text(
        0.0,
        -0.62 * ell_h,
        rf"$b/a = {b_over_a:.2f}$",
        ha="center",
        va="top",
        fontsize=7,
        color="0.15",
    )
    axins.set_xlim(-span, span)
    axins.set_ylim(-0.62 * ell_h - 0.28 * span, span * 0.72)
    axins.set_aspect("equal", adjustable="box")


def _plot_ellipsoid_fit(provid, t, y, dy, t0, spheres, ellipsoid, band=BAND, path=None):
    if path is None:
        path = output_path(
            "ellipsoid_multiband" if band == "Lk" else "ellipsoid", provid
        )
    if path.exists():
        return path
    curve_color = band_colors.get(band, "k")
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(10, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    ph_data = _rotational_phase(t, t0, ellipsoid["freq"])
    ax.errorbar(
        ph_data,
        y,
        yerr=dy,
        fmt="o",
        ms=4,
        alpha=0.6,
        color=curve_color,
        ecolor="0.7",
        label="Magnitude",
    )
    ax.plot(
        ellipsoid["ph_grid"],
        ellipsoid["m_grid"],
        "k",
        lw=2.5,
        label=rf"Ellipsoid, equator-on $b/a={ellipsoid['b_over_a']:.2f}$, $\Delta m={ellipsoid['amp_ell']:.2f}$",
    )
    t_grid = t0 + ellipsoid["ph_grid"] * 2.0 / ellipsoid["freq"]
    m_sph = touching_spheres_mag(
        _rotational_phase(t_grid, t0, spheres["freq"]),
        spheres["m_sep"],
        spheres["q_fit"],
        spheres["cos2phi"],
        spheres["theta0"],
    )
    ax.plot(
        ellipsoid["ph_grid"],
        m_sph,
        color="0.45",
        ls="--",
        lw=1.4,
        label=rf"Touching spheres, $q={spheres['q_sph']:.2f}$",
    )
    ax.axhline(
        ellipsoid["m_mid_ell"],
        color="C3",
        ls="--",
        lw=1.2,
        label=f"Midpoint = {ellipsoid['m_mid_ell']:.3f}",
    )
    ax.set_ylabel("Magnitude")
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", frameon=True)
    ax.set_title(
        f"Triaxial-ellipsoid fit to {provid} ({PROP},{band}), "
        rf"$P_{{rot}}$ = {ellipsoid['p_rot_ell_h']:.3f} h, "
        rf"$\chi^2_\nu$ = {ellipsoid['chi_r_ell']:.1f} (spheres {spheres['chi_r_sph']:.1f})"
    )
    model = ellipsoid_mag(
        ph_data, ellipsoid["m_ab"], ellipsoid["b_over_a"], 1.0, 0.0, ellipsoid["theta0"]
    )
    axr.errorbar(
        ph_data,
        y - model,
        yerr=dy,
        fmt="o",
        ms=3,
        alpha=0.6,
        color=curve_color,
        ecolor="0.7",
    )
    axr.axhline(0, color="k", lw=1)
    axr.set_ylabel("Data - model")
    axr.set_xlabel("Phase at the rotation period ($2P$)")
    axr.invert_yaxis()
    axr.set_xlim(0, 1)
    axr.grid(True, alpha=0.3)
    fig.tight_layout()
    return save_figure(fig, path)


def _plot_shape_compare(provid, t, y, dy, t0, freq0, spheres, ellipsoid, band="Lk", path=None):
    """Both shape fits, their previews, and the two residual panels in one figure."""
    if path is None:
        path = output_path("shape_multiband", provid)
    if path.exists():
        return path
    curve_color = band_colors.get(band, "k")
    fig, (ax, ax_sph, ax_ell) = plt.subplots(
        3, 1, figsize=(10, 9), sharex=True, gridspec_kw={"height_ratios": [3.2, 1, 1]}
    )
    ph_show = _rotational_phase(t, t0, freq0)
    ph_grid = spheres["ph_grid"]
    t_show = t0 + ph_grid * 2.0 / freq0
    m_sph_show = touching_spheres_mag(
        _rotational_phase(t_show, t0, spheres["freq"]),
        spheres["m_sep"],
        spheres["q_fit"],
        spheres["cos2phi"],
        spheres["theta0"],
    )
    m_ell_show = ellipsoid_mag(
        _rotational_phase(t_show, t0, ellipsoid["freq"]),
        ellipsoid["m_ab"],
        ellipsoid["b_over_a"],
        1.0,
        0.0,
        ellipsoid["theta0"],
    )
    ax.errorbar(
        ph_show,
        y,
        yerr=dy,
        fmt="o",
        ms=4,
        alpha=0.6,
        color=curve_color,
        ecolor="0.7",
        label="Magnitude",
        zorder=3,
    )
    ax.plot(
        ph_grid,
        m_sph_show,
        "k",
        lw=2.2,
        label=rf"Touching spheres: $q={spheres['q_sph']:.2f}$, "
        rf"$\phi={spheres['phi_sph_deg']:.0f}^\circ$, "
        rf"$\Delta m={spheres['amp_sph']:.2f}$, $P={spheres['p_rot_sph_h']:.2f}$ h",
    )
    ax.plot(
        ph_grid,
        m_ell_show,
        color="C3",
        lw=2.0,
        label=rf"Ellipsoid: $b/a={ellipsoid['b_over_a']:.2f}$, "
        rf"$\Delta m={ellipsoid['amp_ell']:.2f}$, $P={ellipsoid['p_rot_ell_h']:.2f}$ h",
    )
    ax.set_ylabel("Magnitude")
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", frameon=True, fontsize=8)
    ax.set_title(f"Shape fits to {provid} ({PROP},{band})")
    _draw_sphere_inset(ax, spheres["q_sph"], box=[0.012, 0.02, 0.20, 0.30], radii=False)
    _draw_ellipsoid_inset(ax, ellipsoid["b_over_a"], box=[0.788, 0.02, 0.20, 0.30])

    ph_sph = _rotational_phase(t, t0, spheres["freq"])
    res_sph = y - touching_spheres_mag(
        ph_sph, spheres["m_sep"], spheres["q_fit"], spheres["cos2phi"], spheres["theta0"]
    )
    ax_sph.errorbar(
        ph_show, res_sph, yerr=dy, fmt="o", ms=3, alpha=0.6, color=curve_color, ecolor="0.7"
    )
    ax_sph.axhline(0, color="k", lw=1)
    ax_sph.set_ylabel("Data $-$ spheres")
    ax_sph.invert_yaxis()
    ax_sph.grid(True, alpha=0.3)
    ax_sph.text(
        0.98,
        0.92,
        rf"$\chi^2_\nu = {spheres['chi_r_sph']:.1f}$",
        transform=ax_sph.transAxes,
        ha="right",
        va="top",
        fontsize=9,
    )

    ph_ell = _rotational_phase(t, t0, ellipsoid["freq"])
    res_ell = y - ellipsoid_mag(
        ph_ell, ellipsoid["m_ab"], ellipsoid["b_over_a"], 1.0, 0.0, ellipsoid["theta0"]
    )
    ax_ell.errorbar(
        ph_show, res_ell, yerr=dy, fmt="o", ms=3, alpha=0.6, color=curve_color, ecolor="0.7"
    )
    ax_ell.axhline(0, color="k", lw=1)
    ax_ell.set_ylabel("Data $-$ ellipsoid")
    ax_ell.set_xlabel("Phase at the rotation period ($2P$)")
    ax_ell.invert_yaxis()
    ax_ell.set_xlim(0, 1)
    ax_ell.grid(True, alpha=0.3)
    ax_ell.text(
        0.98,
        0.92,
        rf"$\chi^2_\nu = {ellipsoid['chi_r_ell']:.1f}$",
        transform=ax_ell.transAxes,
        ha="right",
        va="top",
        fontsize=9,
    )
    fig.tight_layout()
    return save_figure(fig, path)


SHAPE_COLUMNS = [
    "chi_r_sph",
    "q_sph",
    "q_sph_err",
    "phi_sph_deg",
    "phi_sph_err_deg",
    "amp_sph",
    "m_mid_sph",
    "m_sep",
    "p_rot_sph_h",
    "chi_r_ell",
    "b_over_a",
    "b_over_a_err",
    "amp_ell",
    "m_mid_ell",
    "m_ab",
    "p_rot_ell_h",
]


def shape_fit(provid, fit):
    """Touching-spheres and equator-on ellipsoid fits of the folded light curve.

    Writes phase_spheres_fit and phase_ellipsoid_fit, and returns the columns
    stored in asteroid_analysis_results.csv, including chi_r_sph and chi_r_ell.
    """
    t = np.asarray(fit["t_band"], dtype=float)
    y = np.asarray(fit["y_band"], dtype=float)
    dy = np.maximum(np.asarray(fit["dy_band"], dtype=float), 1e-4)
    freq0 = float(fit["freq_sel"])
    t0 = float(np.min(t))
    spheres = _fit_touching_spheres(t, y, dy, freq0, t0)
    ellipsoid = _fit_ellipsoid(t, y, dy, freq0, t0)
    band = fit.get("band", BAND)
    paths = [
        _plot_spheres_fit(provid, t, y, dy, t0, spheres, band=band),
        _plot_ellipsoid_fit(provid, t, y, dy, t0, spheres, ellipsoid, band=band),
    ]
    print(
        f"  shape: P_sph = {spheres['p_rot_sph_h']:.3f} h, "
        f"P_ell = {ellipsoid['p_rot_ell_h']:.3f} h, "
        f"chi_r_sph = {spheres['chi_r_sph']:.2f}, q = {spheres['q_sph']:.3f}; "
        f"chi_r_ell = {ellipsoid['chi_r_ell']:.2f}, b/a = {ellipsoid['b_over_a']:.3f}"
    )
    row = {column: spheres[column] for column in SHAPE_COLUMNS if column in spheres}
    row.update(
        {column: ellipsoid[column] for column in SHAPE_COLUMNS if column in ellipsoid}
    )
    return row, paths


def _figure_paths(kind, provid):
    if kind == "shape_fit":
        return [output_path("spheres", provid), output_path("ellipsoid", provid)]
    if kind == "multiband_photometry":
        return [
            # output_path("multiband_lightcurve", provid),
            output_path("multiband_phase_model", provid),
            FIG_DIR / f"{slug(provid)}_band_k.csv",
        ]
    if kind == "shape_fit_multiband":
        return [output_path("shape_multiband", provid)]
    if kind == "periodogram_democratic":
        return []
    return [output_path(kind, provid)]


def _design_matrix(t, flt, f_base, K, filters, ref):
    cols = [(flt == b).astype(float) for b in filters if b != ref] + [np.ones(len(t))]
    names = [f"c_{b}" for b in filters if b != ref] + ["m_ref"]
    ph = 2 * np.pi * f_base * t
    for k in range(1, K + 1):
        cols += [np.sin(k * ph), np.cos(k * ph)]
        names += [f"sin{k}", f"cos{k}"]
    return np.column_stack(cols), names


def _wls_fit(X, y, s, mask):
    beta, *_ = np.linalg.lstsq(X[mask] / s[mask, None], y[mask] / s[mask], rcond=None)
    return beta, y - X @ beta


def _robust_joint_fit(
    t_fit, y_fit, s_fit, f_fit, present, ref, f_base, K, mask0=None, nsig=5.0, n_iter=30
):
    X, names = _design_matrix(t_fit, f_fit, f_base, K, present, ref)
    good = np.ones(len(y_fit), bool) if mask0 is None else mask0.copy()
    jitter = 0.02
    for _ in range(n_iter):
        s = np.sqrt(s_fit**2 + jitter**2)
        beta, r = _wls_fit(X, y_fit, s, good)
        scale = 1.4826 * np.median(np.abs(r[good] - np.median(r[good])))
        jitter = np.sqrt(max(scale**2 - np.median(s_fit[good]) ** 2, 0.0))
        new = np.abs(r) < nsig * max(scale, 1e-4)
        if (new == good).all():
            break
        good = new
    s = np.sqrt(s_fit**2 + jitter**2)
    beta, r = _wls_fit(X, y_fit, s, good)
    A = X[good] / s[good, None]
    cov = np.linalg.pinv(A.T @ A)
    chi2 = np.sum((r[good] / s[good]) ** 2)
    return dict(
        beta=beta,
        cov=cov,
        names=names,
        good=good,
        resid=r,
        jitter=jitter,
        chi2=chi2,
        n=int(good.sum()),
        k=X.shape[1],
        rms=np.sqrt(np.mean(r[good] ** 2)),
    )


def _chi2_fixed(t_fit, y_fit, s_fit, f_fit, present, ref, f_base, K, mask, jitter):
    X, _ = _design_matrix(t_fit, f_fit, f_base, K, present, ref)
    s = np.sqrt(s_fit**2 + jitter**2)
    _, r = _wls_fit(X, y_fit, s, mask)
    return np.sum((r[mask] / s[mask]) ** 2), X.shape[1]


def multiband_photometry(df_obj, provid, fit=None, save=True):
    """Multiband joint fit, color indices, and bolometric (Lk) light curve."""
    alias = {"Lg": "g", "Lr": "r", "Li": "i", "Lz": "z", "Ly": "y", "Lu": "u"}
    filt_order = ["u", "g", "r", "i", "z", "y"]
    filt = lambda b: alias.get(b, b)
    EXCLUDE_FILTERS = ["u"]
    REF_FILTER = "r"
    NSIG_CLIP = 5.0
    use_filters = [f for f in filt_order if f not in EXCLUDE_FILTERS]

    df_color = df_obj.dropna(subset=[PROP, "band"]).copy()
    df_color["filt"] = df_color["band"].map(filt)
    rms = df_color["rmsmag"].to_numpy(float)
    df_color["rmsmag"] = np.where(np.isfinite(rms) & (rms > 0), rms, 0.01)

    df_fit = df_color[df_color["filt"].isin(use_filters)]
    present = [f for f in use_filters if (df_fit["filt"] == f).any()]
    if len(present) < 2:
        print(
            f"  skip multiband_photometry for {provid}: only {len(present)} filter(s) in {present}"
        )
        return None

    if REF_FILTER not in present:
        counts = {f: (df_fit["filt"] == f).sum() for f in present}
        REF_FILTER = max(counts, key=counts.get)

    t_fit = df_fit["t_days"].to_numpy(float)
    y_fit = df_fit[PROP].to_numpy(float)
    s_fit = df_fit["rmsmag"].to_numpy(float)
    f_fit = df_fit["filt"].to_numpy()

    if fit is not None and "freq_sel" in fit and np.isfinite(fit["freq_sel"]):
        f_ph0 = float(fit["freq_sel"])
    else:
        f_scan = np.arange(1.0, 30.0, 0.005)
        chi_scan = np.empty(len(f_scan))
        all_good = np.ones(len(y_fit), bool)
        for j, fq in enumerate(f_scan):
            chi_scan[j] = _chi2_fixed(
                t_fit, y_fit, s_fit, f_fit, present, REF_FILTER, fq, 3, all_good, 0.02
            )[0]
        order = np.argsort(chi_scan)
        cands = []
        for j in order:
            if all(abs(f_scan[j] - f_scan[i]) > 0.3 for i in cands):
                cands.append(j)
            if len(cands) == 4:
                break
        f_ph0 = f_scan[cands[0]]

    p_rot_used_h = 48.0 / f_ph0

    fit0 = _robust_joint_fit(
        t_fit, y_fit, s_fit, f_fit, present, REF_FILTER, f_ph0, 3, nsig=NSIG_CLIP
    )
    good0, jit0 = fit0["good"], fit0["jitter"]
    n0 = max(good0.sum(), 1)

    def best_on_grid(f_center, half_width, step, K):
        grid = np.arange(f_center - half_width, f_center + half_width + step / 2, step)
        chis = np.array(
            [
                _chi2_fixed(
                    t_fit, y_fit, s_fit, f_fit, present, REF_FILTER, fq, K, good0, jit0
                )[0]
                for fq in grid
            ]
        )
        return grid, chis

    f_rot_c = f_ph0 / 2.0
    res_models = {}
    for K in (2, 3, 4, 5, 6, 7, 8, 9, 10):
        grid, chis = best_on_grid(f_rot_c, 0.02, 0.0002, K)
        j = int(np.argmin(chis))
        kpar = 2 * K + len(present)
        bic_val = chis[j] + kpar * np.log(n0)
        res_models[K] = (grid[j], bic_val, grid, chis)

    K_harm, (f_rot, bic_multi, grid_best, chis_best) = min(
        res_models.items(), key=lambda kv: kv[1][1]
    )
    P_rot_h = 24.0 / f_rot

    final_fit = _robust_joint_fit(
        t_fit,
        y_fit,
        s_fit,
        f_fit,
        present,
        REF_FILTER,
        f_rot,
        K_harm,
        mask0=good0,
        nsig=NSIG_CLIP,
    )
    good_fit = final_fit["good"]
    beta = dict(zip(final_fit["names"], final_fit["beta"]))
    err = dict(zip(final_fit["names"], np.sqrt(np.diag(final_fit["cov"]))))

    color = {REF_FILTER: 0.0}
    sig_color = {REF_FILTER: 0.0}
    for f in present:
        if f != REF_FILTER:
            color[f] = beta[f"c_{f}"]
            sig_color[f] = err[f"c_{f}"]

    df_obj = df_obj.copy()
    df_obj["outlier"] = False
    df_obj.loc[df_fit.index[~good_fit], "outlier"] = True

    real = df_obj.dropna(subset=[PROP, "band"]).copy()
    real["filt"] = real["band"].map(filt)
    real = real[~real["outlier"] & real["filt"].isin(use_filters)].copy()
    real["estimated"] = False
    real["source_filt"] = pd.NA
    rms = real["rmsmag"].to_numpy(float)
    real["rmsmag"] = np.where(np.isfinite(rms) & (rms > 0), rms, 0.01)

    sources = real[real["filt"].isin(color)].copy()
    pieces = [real]
    sig_mean = dict(sig_color)
    long_name = {short: long for long, short in alias.items()}
    for b, cb in color.items():
        src = sources[sources["filt"] != b].copy()
        if src.empty:
            continue
        sig_a = src["filt"].map(sig_mean).to_numpy(float)
        src[PROP] = (
            src[PROP].to_numpy(float) + cb - src["filt"].map(color).to_numpy(float)
        )
        src["rmsmag"] = np.sqrt(
            src["rmsmag"].to_numpy(float) ** 2 + sig_mean[b] ** 2 + sig_a**2
        )
        src["source_filt"] = src["filt"]
        src["filt"] = b
        src["band"] = long_name.get(b, f"L{b}")
        src["estimated"] = True
        pieces.append(src)

    df_obj_comp = pd.concat(pieces, ignore_index=True)
    df_obj_comp = df_obj_comp.sort_values(["obstime", "filt", "estimated"]).reset_index(
        drop=True
    )

    base = df_obj_comp[df_obj_comp["band"] != "Lk"].copy()
    sig_mu = dict(sig_mean)
    rows_k = []
    for obstime, g in base.groupby("obstime", sort=False):
        g = (
            g.sort_values("estimated")
            .groupby("filt", as_index=False)
            .first()
            .reset_index(drop=True)
        )
        mag = g[PROP].to_numpy(float)
        flux = 10.0 ** (-0.4 * mag)
        ftot = flux.sum()
        w = flux / ftot
        mag_k = -2.5 * np.log10(ftot)

        coef_obs, coef_mu, obs_rms = {}, {}, {}
        for i in range(len(g)):
            row_i = g.iloc[i]
            f, wi = row_i["filt"], w[i]
            if not row_i["estimated"]:
                coef_obs[f] = coef_obs.get(f, 0.0) + wi
                obs_rms[f] = row_i["rmsmag"]
            else:
                a = row_i["source_filt"]
                coef_obs[a] = coef_obs.get(a, 0.0) + wi
                coef_mu[f] = coef_mu.get(f, 0.0) + wi
                coef_mu[a] = coef_mu.get(a, 0.0) - wi
                if a not in obs_rms:
                    obs_match = g.loc[(g["filt"] == a) & (~g["estimated"]), "rmsmag"]
                    obs_rms[a] = obs_match.iloc[0] if len(obs_match) > 0 else 0.01

        var = sum((c * obs_rms.get(f, 0.01)) ** 2 for f, c in coef_obs.items())
        var += sum((c * sig_mu.get(f, 0.01)) ** 2 for f, c in coef_mu.items())

        anchor = (
            g.loc[~g["estimated"]].iloc[0] if (~g["estimated"]).any() else g.iloc[0]
        )
        out_k = anchor.copy()
        out_k["obs_filt"] = anchor["filt"]
        present_k = [f for f in filt_order if f in set(g["filt"])]
        out_k["band"] = "Lk"
        out_k["filt"] = "k"
        out_k[PROP] = mag_k
        out_k["rmsmag"] = np.sqrt(var)
        out_k["estimated"] = True
        out_k["source_filt"] = "+".join(present_k)
        rows_k.append(out_k)

    df_lk = pd.DataFrame(rows_k)
    df_obj_comp = pd.concat([base, df_lk], ignore_index=True)
    df_obj_comp = df_obj_comp.sort_values(["obstime", "filt", "estimated"]).reset_index(
        drop=True
    )

    t_band_k = df_lk["t_days"].values
    y_band_k = df_lk[PROP].values
    dy_band_k = np.sqrt(
        df_lk["rmsmag"].fillna(0.01).values ** 2 + final_fit["jitter"] ** 2
    )
    ls_k = LombScargle(t_band_k, y_band_k, dy_band_k, nterms=K_harm)
    t0_k = float(t_band_k.min())
    phase_eval = np.linspace(0, 1, 1000, endpoint=False)
    y_model_eval = ls_k.model(t0_k + phase_eval / f_rot, f_rot)
    amp_lsm_k = float(np.ptp(y_model_eval))

    resid_k = y_band_k - ls_k.model(t_band_k, f_rot)
    med_res = float(np.median(resid_k))
    mad_res = float(1.4826 * np.median(np.abs(resid_k - med_res)))
    scale_res = max(mad_res, 0.01)
    inliers_k = np.abs(resid_k - med_res) <= 3.0 * scale_res
    y_k_clean = y_band_k[inliers_k] if inliers_k.sum() >= 5 else y_band_k
    amp_photo_k = float(np.percentile(y_k_clean, 99) - np.percentile(y_k_clean, 1))
    df_lk["resid"] = resid_k
    print(
        f"  {provid} (multiband Lk): P_rot used = {p_rot_used_h:.3f} h, "
        f"P_rot = {P_rot_h:.3f} h, "
        f"amp_lsm = {amp_lsm_k:.3f}, amp_photo = {amp_photo_k:.3f}, BIC = {bic_multi:.1f}"
    )

    ADJACENT = {"u": "g", "g": "u", "y": "z", "z": "y"}
    excl_colors = {}
    for fx in EXCLUDE_FILTERS:
        df_x = df_color[df_color["filt"] == fx]
        if len(df_x) == 0:
            continue
        tx = df_x["t_days"].to_numpy(float)
        X_ref, _ = _design_matrix(
            tx, np.full(len(tx), REF_FILTER), f_rot, K_harm, present, REF_FILTER
        )
        m_ref_model = X_ref @ final_fit["beta"]
        d_x = df_x[PROP].to_numpy(float) - m_ref_model
        s_x = df_x["rmsmag"].to_numpy(float)
        good_x = np.ones(len(d_x), bool)
        for _ in range(30):
            scale_x = 1.4826 * np.median(np.abs(d_x[good_x] - np.median(d_x[good_x])))
            w_x = 1.0 / np.maximum(s_x, 1e-4) ** 2
            c_x = np.sum(w_x[good_x] * d_x[good_x]) / np.sum(w_x[good_x])
            new_x = np.abs(d_x - c_x) < NSIG_CLIP * max(scale_x, np.median(s_x))
            if (new_x == good_x).all():
                break
            good_x = new_x
        n_x = int(good_x.sum())
        chi2_red = float(
            np.sum(((d_x[good_x] - c_x) / s_x[good_x]) ** 2) / max(n_x - 1, 1)
        )
        err_x = float(np.sqrt(1.0 / np.sum(w_x[good_x]) * max(chi2_red, 1.0)))
        adj = ADJACENT.get(fx, REF_FILTER)
        c_adj = c_x - color.get(adj, 0.0)
        e_adj = float(np.sqrt(err_x**2 + sig_color.get(adj, 0.0) ** 2))
        excl_colors[fx] = {
            "fx": fx,
            "ref": REF_FILTER,
            "color_ref": c_x,
            "err_ref": err_x,
            "adj": adj,
            "color_adj": c_adj,
            "err_adj": e_adj,
            "n_pts": n_x,
            "chi2_red": chi2_red,
        }
        print(
            f"  {provid} ({fx}): {fx} - {REF_FILTER} = {c_x:.3f} ± {err_x:.3f}, {fx} - {adj} = {c_adj:+.3f} ± {e_adj:.3f}"
        )

    summary_row = {
        "p_rot_multi_used_h": p_rot_used_h,
        "p_rot_multi_h": P_rot_h,
        "f_rot_multi": f_rot,
        "k_harm_multi": K_harm,
        "bic": float(bic_multi),
        "amp_lsm_k": amp_lsm_k,
        "amp_photo_k": amp_photo_k,
        "rms_multi": float(final_fit["rms"]),
        "n_outliers_multi": int(len(df_fit) - final_fit["n"]),
    }
    # Contiguous pairs: u-g, g-r, r-i, i-z, z-y
    if "u" in excl_colors and "g" in present:
        summary_row["u-g"] = excl_colors["u"]["color_adj"]
        summary_row["u-g_err"] = excl_colors["u"]["err_adj"]
    elif "u" in present and "g" in present:
        ug = color["u"] - color["g"]
        ug_err = (
            float(np.sqrt(sig_color.get("u", 0) ** 2 + sig_color.get("g", 0) ** 2))
            if REF_FILTER != "g"
            else sig_color["u"]
        )
        summary_row["u-g"] = ug
        summary_row["u-g_err"] = ug_err

    if "g" in present and "r" in present:
        gr = color["g"] - color["r"]
        gr_err = (
            float(np.sqrt(sig_color.get("g", 0) ** 2 + sig_color.get("r", 0) ** 2))
            if REF_FILTER != "r"
            else sig_color["g"]
        )
        summary_row["g-r"] = gr
        summary_row["g-r_err"] = gr_err

    if "r" in present and "i" in present:
        ri = color["r"] - color["i"]
        ri_err = (
            float(np.sqrt(sig_color.get("r", 0) ** 2 + sig_color.get("i", 0) ** 2))
            if REF_FILTER != "r"
            else sig_color["i"]
        )
        summary_row["r-i"] = ri
        summary_row["r-i_err"] = ri_err

    if "i" in present and "z" in present:
        iz = color["i"] - color["z"]
        c_i_name = "c_i" if REF_FILTER != "i" else None
        c_z_name = "c_z" if REF_FILTER != "z" else None
        if c_i_name in final_fit["names"] and c_z_name in final_fit["names"]:
            i_idx = final_fit["names"].index(c_i_name)
            z_idx = final_fit["names"].index(c_z_name)
            cov_mat = final_fit["cov"]
            iz_err = float(
                np.sqrt(
                    max(
                        cov_mat[i_idx, i_idx]
                        + cov_mat[z_idx, z_idx]
                        - 2 * cov_mat[i_idx, z_idx],
                        0.0,
                    )
                )
            )
        else:
            iz_err = float(
                np.sqrt(sig_color.get("i", 0) ** 2 + sig_color.get("z", 0) ** 2)
            )
        summary_row["i-z"] = iz
        summary_row["i-z_err"] = iz_err

    if "z" in present and "y" in present:
        zy = color["z"] - color["y"]
        c_z_name = "c_z" if REF_FILTER != "z" else None
        c_y_name = "c_y" if REF_FILTER != "y" else None
        if c_z_name in final_fit["names"] and c_y_name in final_fit["names"]:
            z_idx = final_fit["names"].index(c_z_name)
            y_idx = final_fit["names"].index(c_y_name)
            cov_mat = final_fit["cov"]
            zy_err = float(
                np.sqrt(
                    max(
                        cov_mat[z_idx, z_idx]
                        + cov_mat[y_idx, y_idx]
                        - 2 * cov_mat[z_idx, y_idx],
                        0.0,
                    )
                )
            )
        else:
            zy_err = float(
                np.sqrt(sig_color.get("z", 0) ** 2 + sig_color.get("y", 0) ** 2)
            )
        summary_row["z-y"] = zy
        summary_row["z-y_err"] = zy_err

    # Excluded filters vs reference filter (e.g. u-r)
    for fx, ex in excl_colors.items():
        summary_row[f"{fx}-{ex['ref']}"] = ex["color_ref"]
        summary_row[f"{fx}-{ex['ref']}_err"] = ex["err_ref"]

    written_paths = []
    if save:
        path_k_csv = FIG_DIR / f"{slug(provid)}_band_k.csv"
        df_lk["provid"] = provid
        cols_export = [
            c
            for c in [
                "provid",
                "obstime",
                "t_days",
                "night",
                "band",
                "filt",
                "obs_filt",
                "source_filt",
                PROP,
                "rmsmag",
            ]
            if c in df_lk.columns
        ]
        df_lk[cols_export].to_csv(path_k_csv, index=False)
        written_paths.append(path_k_csv)

        # 1. lightcurve_allnights
        # path_lc = output_path("multiband_lightcurve", provid)
        # nights = df_obj_comp["night"].unique()
        # n_nights = len(nights)
        # fig, axs = plt.subplots(n_nights + 1, 1, figsize=(10, 3 * n_nights))
        # if n_nights + 1 == 1:
        #     axs = [axs]
        # ax = axs[0]
        # for b in ["Lk", "Lr", "Lg", "Li", "Lz", "Ly", "Lu"]:
        #     sb = df_obj_comp[df_obj_comp["band"] == b]
        #     if len(sb) > 0:
        #         ax.errorbar(
        #             sb["obstime"],
        #             sb[PROP],
        #             yerr=sb["rmsmag"],
        #             fmt="o",
        #             label=f"{b} (N={len(sb)})",
        #             color=band_colors.get(b, "k"),
        #             markersize=4,
        #             alpha=0.85,
        #         )
        # ax.set_xlabel("Observation Time")
        # ax.set_ylabel(PROP)
        # ax.tick_params(axis="x", rotation=30)
        # ax.set_title(f"Light Curve of Asteroid {provid} (all Nights)")
        # ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1, frameon=True)
        # ax.grid(True, alpha=0.3)
        # for i, night in enumerate(nights):
        #     night_str = pd.to_datetime(night).strftime("%Y-%m-%d")
        #     ax = axs[i + 1]
        #     t_night = df_obj_comp[df_obj_comp["night"] == night].copy()
        #     t_night_0 = min(t_night["obstime"])
        #     for b in ["Lk", "Lr", "Lg", "Li", "Lz", "Ly", "Lu"]:
        #         sb = t_night[t_night["band"] == b]
        #         if len(sb) > 0:
        #             sb_h = (sb["obstime"] - t_night_0).dt.total_seconds() / 3600.0
        #             ax.errorbar(
        #                 sb_h,
        #                 sb[PROP],
        #                 yerr=sb["rmsmag"],
        #                 fmt="o",
        #                 label=f"{b} (N={len(sb)})",
        #                 color=band_colors.get(b, "k"),
        #                 markersize=4,
        #                 alpha=0.85,
        #             )
        #     ax.set_title(f"Night {night_str}")
        #     ax.set_xlabel("Time since start [h]")
        #     ax.set_ylabel(PROP)
        #     ax.legend(
        #         loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1, frameon=True
        #     )
        #     ax.tick_params(axis="x", rotation=30)
        # fig.tight_layout()
        # path_lc = save_figure(fig, path_lc)
        # written_paths.append(path_lc)

        # 2. Folded phase model
        path_pm = output_path("multiband_phase_model", provid)
        df_band_k = df_lk.copy()
        phase_rot = np.mod((t_band_k - t0_k) * f_rot, 1.0)
        phase_rot_grid = np.linspace(0, 1, 400)
        y_model_rot = ls_k.model(t0_k + phase_rot_grid / f_rot, f_rot)

        fig, (ax, axr) = plt.subplots(
            2, 1, figsize=(10, 6.5), sharex=True, height_ratios=[3, 1.2]
        )
        for f in use_filters:
            sel = (df_band_k["obs_filt"] == f).to_numpy()
            if not sel.any():
                continue
            ax.errorbar(
                phase_rot[sel],
                y_band_k[sel],
                yerr=dy_band_k[sel],
                fmt="o",
                ms=4,
                alpha=0.8,
                color=band_colors.get(f"L{f}", "k"),
                ecolor="0.7",
                label=f"observed in {f} (N={sel.sum()})",
            )
            axr.plot(
                phase_rot[sel],
                resid_k[sel],
                "o",
                ms=3,
                alpha=0.8,
                color=band_colors.get(f"L{f}", "k"),
            )
        ax.plot(
            phase_rot_grid,
            y_model_rot,
            color="k",
            lw=2,
            label=f"Periodogram model ({K_harm} harmonics), $P_{{rot}}$ = {P_rot_h:.4f} h",
        )
        ax.set_ylim(y_band_k.max() + 0.08, y_band_k.min() - 0.25)
        ax.set_ylabel("Magnitude")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper center", ncol=2, fontsize=8, frameon=True)
        ax.set_title(
            f"Folded light curve of {provid} ({PROP},Lk), $P_{{rot}}$ = {P_rot_h:.4f} h"
        )
        axr.axhline(0, color="k", lw=1)
        axr.set_ylabel("Data - model")
        axr.set_xlabel("Phase at the rotation period")
        axr.set_xlim(0, 1)
        axr.grid(True, alpha=0.3)
        fig.tight_layout()
        path_pm = save_figure(fig, path_pm)
        written_paths.append(path_pm)

    return {
        "provid": provid,
        "df_band_k": df_lk,
        "df_obj_comp": df_obj_comp,
        "f_rot": f_rot,
        "p_rot_h": P_rot_h,
        "amp_lsm_k": amp_lsm_k,
        "amp_photo_k": amp_photo_k,
        "freq_sel": 2.0 * f_rot,
        "t_band": df_lk["t_days"].values,
        "y_band": df_lk[PROP].values,
        "dy_band": df_lk["rmsmag"].values,
        "K_harm": K_harm,
        "color": color,
        "sig_color": sig_color,
        "excl_colors": excl_colors,
        "paths": written_paths,
        "summary_row": summary_row,
    }


def shape_fit_multiband(provid, multi_fit, save=True):
    """Touching-spheres and equator-on ellipsoid fits on the multiband (Lk) data.

    Writes one figure, <provid>_phase_shape_fit_Lk.webp, with both model curves,
    the sphere and ellipsoid previews, and a residual panel for each fit.
    """
    if multi_fit is None or "df_band_k" not in multi_fit:
        raise RuntimeError(
            f"shape_fit_multiband cannot be executed for {provid} because multiband_photometry has not been run first."
        )
    df_band_k = multi_fit["df_band_k"]
    t = np.asarray(df_band_k["t_days"], dtype=float)
    y = np.asarray(df_band_k[PROP], dtype=float)
    dy = np.maximum(np.asarray(df_band_k["rmsmag"], dtype=float), 1e-4)
    freq0 = float(multi_fit["freq_sel"])
    t0 = float(np.min(t))

    spheres = _fit_touching_spheres(t, y, dy, freq0, t0)
    ellipsoid = _fit_ellipsoid(t, y, dy, freq0, t0)

    paths = []
    if save:
        paths = [
            _plot_shape_compare(
                provid, t, y, dy, t0, freq0, spheres, ellipsoid, band="Lk"
            )
        ]

    row = {
        f"{column}_k": spheres[column] for column in SHAPE_COLUMNS if column in spheres
    }
    row.update(
        {
            f"{column}_k": ellipsoid[column]
            for column in SHAPE_COLUMNS
            if column in ellipsoid
        }
    )

    print(
        f"  shape (multiband Lk): P_sph = {spheres['p_rot_sph_h']:.3f} h, "
        f"P_ell = {ellipsoid['p_rot_ell_h']:.3f} h, "
        f"chi_r_sph = {spheres['chi_r_sph']:.2f}, q = {spheres['q_sph']:.3f}; "
        f"chi_r_ell = {ellipsoid['chi_r_ell']:.2f}, b/a = {ellipsoid['b_over_a']:.3f}"
    )
    return row, paths


SPECTRUM_WAVELENGTH = {"Lg": 483, "Lr": 622, "Li": 755, "Lz": 869, "Ly": 971}
SPECTRUM_SOLAR = {"Lg": 0.45, "Lr": 0.00, "Li": -0.12, "Lz": -0.16, "Ly": -0.18}
SPECTRUM_ORDER = ("Lg", "Lr", "Li", "Lz", "Ly")
SPECTRUM_SHORT = {"Lg": "g", "Lr": "r", "Li": "i", "Lz": "z", "Ly": "y"}
SPECTRUM_COLUMNS = [
    column
    for short in ("g", "i", "z", "y")
    for column in (
        f"{short}-r_spec",
        f"{short}-r_spec_err",
        f"R_{short}",
        f"R_{short}_err",
        f"n_{short}_spec",
    )
]


def _spectrum_color_wrt_r(short, color, sig_color, excl_colors):
    """Color m_short - m_r from the joint light-curve fit."""
    if short == "r":
        return 0.0, 0.0
    if short in color:
        value = color[short] - color.get("r", 0.0)
        error = float(
            np.sqrt(sig_color.get(short, 0.0) ** 2 + sig_color.get("r", 0.0) ** 2)
        )
        return value, error
    excluded = excl_colors.get(short)
    if excluded is None:
        return None
    value = excluded["color_ref"]
    error = excluded["err_ref"]
    if excluded.get("ref", "r") != "r":
        value = value - color.get("r", 0.0)
        error = float(np.sqrt(error**2 + sig_color.get("r", 0.0) ** 2))
    return value, error


def spectrum(df_obj, provid, multi_fit, save=True):
    """Solar-normalized reflectance from the joint light-curve colors.

    The shared rotation model removes the phase, so a color does not require
    the two filters to have been observed on the same night.
    """
    path = output_path("spectrum", provid)
    if multi_fit is None or "color" not in multi_fit:
        print(f"  skip spectrum for {provid}: the joint light curve has not been fit")
        return None
    color = multi_fit["color"]
    sig_color = multi_fit["sig_color"]
    excl_colors = multi_fit.get("excl_colors") or {}

    points = []
    row = {}
    print(f"  spectrum of {provid}")
    for band in SPECTRUM_ORDER:
        short = SPECTRUM_SHORT[band]
        if short == "r" and short not in color:
            continue
        if short != "r" and short not in color and short not in excl_colors:
            row[f"{short}-r_spec"] = np.nan
            row[f"{short}-r_spec_err"] = np.nan
            row[f"R_{short}"] = np.nan
            row[f"R_{short}_err"] = np.nan
            row[f"n_{short}_spec"] = 0
            continue
        got = _spectrum_color_wrt_r(short, color, sig_color, excl_colors)
        if got is None:
            continue
        color_mr, color_err = got
        excess = color_mr - SPECTRUM_SOLAR[band]
        reflectance = 10.0 ** (-0.4 * excess)
        reflectance_err = reflectance * 0.4 * np.log(10.0) * color_err
        n_points = int((df_obj["band"] == band).sum())
        points.append((band, SPECTRUM_WAVELENGTH[band], reflectance, reflectance_err))
        if short != "r":
            row[f"{short}-r_spec"] = color_mr
            row[f"{short}-r_spec_err"] = color_err
            row[f"R_{short}"] = reflectance
            row[f"R_{short}_err"] = reflectance_err
            row[f"n_{short}_spec"] = n_points
        print(
            f"  {provid} {band}: (m-mr) = {color_mr:.3f} ± {color_err:.3f}, "
            f"R = {reflectance:.3f} ± {reflectance_err:.3f}, points = {n_points}"
        )
    if len(points) < 2:
        print(f"  skip spectrum for {provid}: fewer than two filters with a color")
        return None

    if save and not path.exists():
        fig, ax = plt.subplots(figsize=(10, 6))
        wavelength = np.linspace(450, 1050, 200)
        ax.plot(
            wavelength,
            1.0
            + 0.5 * (wavelength - 622) / 300.0
            - 0.3 * np.exp(-((wavelength - 950) / 150.0) ** 2),
            "--",
            color="goldenrod",
            alpha=0.75,
            label="S-type template",
        )
        ax.plot(
            wavelength,
            1.0 + 0.05 * (wavelength - 622) / 400.0,
            "--",
            color="gray",
            alpha=0.75,
            label="C-type template",
        )
        ax.plot(
            wavelength,
            1.0 + 0.9 * (wavelength - 622) / 400.0,
            "--",
            color="firebrick",
            alpha=0.75,
            label="D-type template",
        )
        ax.errorbar(
            [item[1] for item in points],
            [item[2] for item in points],
            yerr=[item[3] for item in points],
            fmt="o-",
            color="navy",
            ecolor="royalblue",
            elinewidth=2,
            capsize=5,
            markersize=8,
            linewidth=2,
            zorder=5,
            label=provid,
        )
        for band, wave, reflectance, _err in points:
            ax.annotate(
                SPECTRUM_SHORT[band],
                (wave, reflectance),
                textcoords="offset points",
                xytext=(0, 12),
                ha="center",
                fontsize=11,
                fontweight="bold",
                color="navy",
            )
        ax.axhline(1.0, color="0.5", lw=0.8, zorder=0)
        ax.set_title(f"Solar-normalized reflectance of {provid}")
        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel(r"$R(\lambda)\,/\,R(r)$")
        ax.set_xlim(440, 1060)
        reflectances = [item[2] for item in points]
        pad = 0.15
        ax.set_ylim(min(0.55, min(reflectances) - pad), max(2.15, max(reflectances) + pad))
        ax.legend(frameon=True, loc="upper left")
        fig.tight_layout()
        path = save_figure(fig, path)
    elif not save or not path.exists():
        path = None

    return {"row": row, "path": path}


def process(
    df_set,
    provid,
    df_rot,
    analysis,
    designations,
    force=False,
    band=None,
    exclude_lx=EXCLUDE_LX,
    Nb=3,
):
    print(f"{provid}")
    kinds = [kind for kind in ANALYSIS if kind in analysis and kind != "compare"]
    expected = [path for kind in kinds for path in _figure_paths(kind, provid)]
    csv_path = FIG_DIR / "asteroid_analysis_results.csv"
    csv_ok = False
    if csv_path.exists():
        try:
            df_csv = pd.read_csv(csv_path)
            csv_ok = (provid in df_csv["provid"].values) and (
                "epoch_utc_ini" in df_csv.columns
            )
            requested = BAND if band is None or str(band).strip() == "" else str(band).strip()
            if csv_ok and requested != "Lx":
                if "band_periodogram" not in df_csv.columns:
                    csv_ok = requested == BAND
                else:
                    stored = df_csv.loc[df_csv["provid"] == provid, "band_periodogram"]
                    if stored.empty or pd.isna(stored.iloc[0]):
                        csv_ok = requested == BAND
                    elif str(stored.iloc[0]) != requested:
                        csv_ok = False
            if csv_ok and "periodogram_democratic" in analysis:
                if "p_rot_dem_h" not in df_csv.columns:
                    csv_ok = False
        except Exception:
            csv_ok = False

    if not force and kinds and all(path.exists() for path in expected) and csv_ok:
        print("  results and figures already exist")
        return "skip", None

    try:
        df_obj = prepare_object(df_set, provid, designations)
        if df_obj.empty:
            print("  no observations")
            return None, "Sin observaciones fotométricas en el catálogo"
        if df_obj[PROP].notna().sum() == 0:
            print(f"  no {PROP} values after merging the geometry table")
            return (
                None,
                f"Sin valores de {PROP} tras cruzar tabla de geometría (designación MPC ausente en el catálogo de designaciones)",
            )

        written = []
        if "lightcurve" in analysis:
            written.append(plot_lightcurve_h(df_obj, provid))

        fit = None
        needs_single = any(
            k in analysis
            for k in (
                "periodogram",
                "periodogram_democratic",
                "phase",
                "phase_model",
                "shape_fit",
                "timemodel",
                "literature",
            )
        )
        if needs_single:
            band_use = resolve_band(df_obj, band, exclude_lx=exclude_lx)
            fit = periodogram(
                df_obj, provid, save="periodogram" in analysis, band=band_use
            )
            if fit is None:
                n_band = int(
                    ((df_obj["band"] == band_use) & df_obj[PROP].notna()).sum()
                )
                return (
                    None,
                    f"Menos de 5 observaciones válidas en banda {band_use} con magnitud {PROP} ({n_band} observadas)",
                )
            if fit["path"] is not None:
                written.append(fit["path"])
            if "periodogram_democratic" in analysis:
                dem = periodogram_democratic(
                    df_obj, provid, Nb=Nb, exclude=exclude_lx
                )
                if dem is not None:
                    p_rot = float(dem["p_rot_h"])
                    freq = 48.0 / p_rot
                    fit["freq_sel"] = freq
                    fit["per_sel"] = p_rot
                    fit["best_freq"] = freq
                    fit["best_phot_period"] = p_rot / 2.0
                    fit["best_rot_period"] = p_rot
                    fit["period_source"] = f"LS democratic {dem['source']}"
                    fit["democratic"] = dem
                    stale = []
                    if "phase" in analysis:
                        stale.append(output_path("phase", provid))
                    if "phase_model" in analysis:
                        stale.append(output_path("phase_model", provid))
                    if "timemodel" in analysis:
                        stale.append(output_path("timemodel", provid))
                    if "shape_fit" in analysis:
                        stale.extend(
                            [
                                output_path("spheres", provid),
                                output_path("ellipsoid", provid),
                            ]
                        )
                    if "literature" in analysis:
                        stale.append(output_path("literature", provid))
                    if "multiband_photometry" in analysis:
                        stale.append(output_path("multiband_phase_model", provid))
                    if "shape_fit_multiband" in analysis:
                        stale.append(output_path("shape_multiband", provid))
                    for stale_path in stale:
                        if stale_path.exists():
                            stale_path.unlink()
            if "phase" in analysis:
                written.append(plot_phase_allbands(df_obj, provid, fit))
            if "phase_model" in analysis:
                written.append(Lr_phase_curve(provid, fit))
            if "shape_fit" in analysis:
                shape, paths = shape_fit(provid, fit)
                fit["shape"] = shape
                written.extend(paths)
            if "timemodel" in analysis:
                written.append(plot_timemodel_h(provid, fit))
            if "literature" in analysis:
                written.append(
                    plot_phase_model_all_periods(
                        provid, fit, literature_periods(df_rot, provid)
                    )
                )

        multi_fit = None
        if "multiband_photometry" in analysis:
            multi_fit = multiband_photometry(df_obj, provid, fit=fit)
            if multi_fit is None:
                return (
                    None,
                    "Datos insuficientes para fotometría multibanda (menos de 2 filtros con datos o ajuste no convergente)",
                )
            if fit is None:
                fit = {}
            fit["multiband"] = multi_fit
            written.extend(multi_fit.get("paths", []))

        if "shape_fit_multiband" in analysis:
            if multi_fit is None:
                return (
                    None,
                    "Falta ajuste de fotometría multibanda para calcular la forma",
                )
            try:
                shape_k, paths_k = shape_fit_multiband(provid, multi_fit)
                if fit is None:
                    fit = {}
                fit["shape_k"] = shape_k
                written.extend(paths_k)
            except Exception as e:
                print(f"  error in shape_fit_multiband: {e}")
                return None, f"Error en ajuste de forma multibanda: {e}"

        if "spectrum" in analysis:
            if multi_fit is None:
                return (
                    None,
                    "Falta el ajuste conjunto de la curva de luz para calcular el espectro",
                )
            spec = spectrum(df_obj, provid, multi_fit, save=True)
            if spec is None:
                if fit is None and not written:
                    return (
                        None,
                        "Sin colores del ajuste conjunto para calcular el espectro",
                    )
            else:
                if fit is None:
                    fit = {}
                fit["spectrum"] = spec["row"]
                if spec["path"] is not None:
                    written.append(spec["path"])

        if fit is None and written:
            fit = {"provid": provid}

        for path in written:
            print(f"  {path.name}")
        return fit, None
    except Exception as e:
        print(f"  error processing {provid}: {e}")
        return None, f"Error inesperado durante el procesamiento: {e}"


def get_analysis_columns(analysis, rows=None):
    """Determine the ordered column schema for asteroid_analysis_results.csv."""
    cols = [
        "provid",
        "epoch_utc_ini",
        "epoch_utc_end",
        "n_obs",
        "n_obs_u",
        "n_obs_g",
        "n_obs_r",
        "n_obs_i",
        "n_obs_z",
        "n_obs_y",
    ]
    if "periodogram" in analysis:
        cols.extend(["band_periodogram", "periods_h", "power", "amp_lsm", "amp_photo"])
    if "periodogram_democratic" in analysis:
        cols.extend(
            ["p_rot_dem_h", "votes_dem", "n_series_dem", "bands_dem", "periods_dem_h"]
        )
    if "shape_fit" in analysis:
        cols.extend(SHAPE_COLUMNS)
    if "multiband_photometry" in analysis:
        cols.extend(
            [
                "p_rot_multi_used_h",
                "p_rot_multi_h",
                "f_rot_multi",
                "k_harm_multi",
                "bic",
                "amp_lsm_k",
                "amp_photo_k",
                "rms_multi",
                "n_outliers_multi",
                "u-g",
                "u-g_err",
                "g-r",
                "g-r_err",
                "r-i",
                "r-i_err",
                "i-z",
                "i-z_err",
                "z-y",
                "z-y_err",
                "u-r",
                "u-r_err",
            ]
        )
    if "shape_fit_multiband" in analysis:
        cols.extend([f"{c}_k" for c in SHAPE_COLUMNS])
    if "spectrum" in analysis:
        cols.extend(SPECTRUM_COLUMNS)

    if rows:
        mandatory_set = {
            "provid",
            "epoch_utc_ini",
            "epoch_utc_end",
            "n_obs",
            "n_obs_u",
            "n_obs_g",
            "n_obs_r",
            "n_obs_i",
            "n_obs_z",
            "n_obs_y",
        }
        seen_keys = {k for r in rows for k in r.keys()}
        cols = [c for c in cols if c in mandatory_set or c in seen_keys]
        extra = [k for r in rows for k in r.keys() if k not in cols]
        cols.extend(dict.fromkeys(extra))

    return cols


def write_analysis_results(rows, replaced, analysis):
    """One row per asteroid in asteroid_analysis_results.csv.

    Columns correspond strictly to the analyses specified in the active 'analysis' set.
    """
    if not rows:
        return None
    columns = get_analysis_columns(analysis, rows)
    new = pd.DataFrame(rows)
    for c in columns:
        if c not in new.columns:
            new[c] = np.nan
    new = new.reindex(columns=columns)
    path = FIG_DIR / "asteroid_analysis_results.csv"
    if path.exists():
        try:
            old = pd.read_csv(path)
            old = old[~old["provid"].isin(replaced)]
            for c in columns:
                if c not in old.columns:
                    old[c] = np.nan
            old = old.reindex(columns=columns)
            catalogue = pd.concat([old, new], ignore_index=True)
        except Exception:
            catalogue = new
    else:
        catalogue = new
    integer_cols = [
        "n_obs",
        "n_obs_u",
        "n_obs_g",
        "n_obs_r",
        "n_obs_i",
        "n_obs_z",
        "n_obs_y",
        "k_harm_multi",
        "n_outliers_multi",
        "n_g_spec",
        "n_i_spec",
        "n_z_spec",
        "n_y_spec",
        "votes_dem",
        "n_series_dem",
    ]
    for ic in integer_cols:
        if ic in catalogue.columns:
            catalogue[ic] = catalogue[ic].astype("Int64")

    catalogue.to_csv(path, index=False)
    print(f"{path.name} ({len(catalogue)} asteroids)")
    copy_analysis_docs()
    return path


def copy_analysis_docs():
    """Copy documentation markdown file into the figures directory FIG_DIR."""
    src = ROOT / "asteroid_analysis_results.md"
    if src.exists():
        dest = FIG_DIR / "asteroid_analysis_results.md"
        try:
            shutil.copy2(src, dest)
        except Exception as e:
            print(f"Warning: could not copy {src.name} to {FIG_DIR}: {e}")


def _provid_and_band(item):
    """A provid, or (provid, band). A missing band means Lr. Lx is resolved later."""
    if isinstance(item, (tuple, list)):
        provid = str(item[0]).strip()
        band = item[1] if len(item) > 1 else None
        if band is None or str(band).strip() == "":
            band = BAND
        else:
            band = str(band).strip()
        return provid, band
    return str(item).strip(), BAND


def _selection_tables(provids, analysis, exclude_lx=EXCLUDE_LX):
    """Photometry and observation counts for an explicit provid list."""
    check_analysis(analysis, "analysis")

    band_of = {}
    wanted = []
    for item in provids:
        provid, band = _provid_and_band(item)
        if not provid or provid in band_of:
            continue
        band_of[provid] = band
        wanted.append(provid)
    df_rot = load_rotation_catalog()
    df_set = load_df_set()
    present = set(df_set["provid"].dropna().astype(str))
    missing = [provid for provid in wanted if provid not in present]
    if missing:
        print(f"No observations for: {', '.join(missing)}")
    wanted = [provid for provid in wanted if provid in present]
    if not wanted:
        raise SystemExit("No selected asteroids have observations.")

    df_set = df_set[df_set["provid"].isin(wanted)].copy()
    counts = observation_counts(df_set).set_index("provid")["n_obs"]
    selected = pd.DataFrame(
        {
            "provid": wanted,
            "n_obs": [int(counts[provid]) for provid in wanted],
            "band": [band_of[provid] for provid in wanted],
        }
    )
    selected.attrs["missing"] = missing
    selected.attrs["exclude_lx"] = tuple(exclude_lx)
    print(f"{len(selected)} asteroids selected")
    return selected, df_set, df_rot, analysis, missing


def select_provid_all(Nb=3):
    """Asteroids with more than MIN_N_OBS observations, and the analyses to run."""
    fig_dir = ROOT / "figcomp"
    set_fig_dir(fig_dir)
    analysis = {
        "lightcurve",
        "periodogram",
        "periodogram_democratic",
        "phase",
        "phase_model",
        "shape_fit",
        "timemodel",
        "multiband_photometry",
        "shape_fit_multiband",
        "spectrum",
    }
    check_analysis(analysis, "analysis")

    df_rot = load_rotation_catalog()
    df_set = load_df_set()
    df_obs = observation_counts(df_set)
    selected = df_obs[df_obs["n_obs"] > MIN_N_OBS].reset_index(drop=True)
    selected["band"] = BAND
    print(f"{len(selected)} asteroids with more than {MIN_N_OBS} observations")
    if selected.empty:
        raise SystemExit(f"No asteroids with more than {MIN_N_OBS} observations.")

    df_set = df_set[df_set["provid"].isin(selected["provid"])].copy()
    missing = []
    selected.attrs["missing"] = missing
    selected.attrs["exclude_lx"] = EXCLUDE_LX
    selected.attrs["Nb"] = int(Nb)
    return selected, df_set, df_rot, analysis, missing


def unpack_selection(ret):
    if len(ret) == 5:
        return ret
    selected, df_set, df_rot, analysis = ret
    missing = getattr(selected, "attrs", {}).get("missing", [])
    return selected, df_set, df_rot, analysis, missing


def select_provid_list(provids=None, Nb=3):
    """An explicit provid list, and the analyses to run."""
    fig_dir = ROOT / "figures/figures_list"
    set_fig_dir(fig_dir)
    analysis = {
        "lightcurve",
        "periodogram",
        "periodogram_democratic",
        "phase",
        "phase_model",
        "shape_fit",
        "timemodel",
        "multiband_photometry",
        "shape_fit_multiband",
        "spectrum",
    }
    if provids is None:
        provids = [
            #("2002 TH52", "Lr"),
            #("2025 MM82", "Lr"),
            #("2025 MA19", "Lr"),
            #("2002 PJ24", "Lx"),
            #("2025 MP71", "Lr"),
            ('2025 ML17', "Lr"),
        ]
    excluded_bands = ["Lu", "Ly"]
    selected, df_set, df_rot, analysis, missing = _selection_tables(
        provids, analysis, exclude_lx=excluded_bands
    )
    selected.attrs["Nb"] = int(Nb)
    return selected, df_set, df_rot, analysis, missing


def select_provid_rotation(Nb=3):
    """Provids in asteroids-with-rotation.csv, and the analyses to run."""
    fig_dir = ROOT / "figures/figures_rotation"
    set_fig_dir(fig_dir)
    analysis = {
        #"lightcurve",
        "periodogram",
        "periodogram_democratic",
        "phase",
        "phase_model",
        #"shape_fit",
        #"timemodel",
        "literature",
        "compare",
        "multiband_photometry",
        "shape_fit_multiband",
        "spectrum",
    }
    df_rot = load_rotation_catalog()
    provids = df_rot["provid"].tolist()
    # provids = provids[: len(provids) // 2]
    print(f"{len(provids)} of {len(df_rot)} asteroids with a literature period")
    selected, df_set, df_rot, analysis, missing = _selection_tables(provids, analysis)
    selected.attrs["Nb"] = int(Nb)
    return selected, df_set, df_rot, analysis, missing


def select_provid_nobs(n_min=50, n_max=100, Nb=3):
    """Asteroids with n_obs between n_min and n_max, inclusive."""
    fig_dir = ROOT / "figures/figures_nobs"
    set_fig_dir(fig_dir)
    analysis = {
        "lightcurve",
        "periodogram",
        "periodogram_democratic",
        "phase",
        "phase_model",
        "shape_fit",
        "timemodel",
        "multiband_photometry",
        "shape_fit_multiband",
        "spectrum",
    }
    check_analysis(analysis, "analysis")

    df_rot = load_rotation_catalog()
    df_set = load_df_set()
    df_obs = observation_counts(df_set)
    selected = df_obs[
        (df_obs["n_obs"] >= n_min) & (df_obs["n_obs"] <= n_max)
    ].reset_index(drop=True)
    selected["band"] = BAND
    print(f"{len(selected)} asteroids with {n_min} to {n_max} observations")
    if selected.empty:
        raise SystemExit(f"No asteroids with {n_min} to {n_max} observations.")

    df_set = df_set[df_set["provid"].isin(selected["provid"])].copy()
    missing = []
    selected.attrs["missing"] = missing
    selected.attrs["exclude_lx"] = EXCLUDE_LX
    selected.attrs["Nb"] = int(Nb)
    return selected, df_set, df_rot, analysis, missing


def select_provid_medellin(Nb=3):
    """Provids in medellin-vro-catalogue.csv, and the analyses to run."""
    fig_dir = ROOT / "figures/figures_medellin"
    set_fig_dir(fig_dir)
    analysis = {
        "lightcurve",
        #"periodogram",
        "periodogram_democratic",
        "phase",
        #"phase_model",
        #"shape_fit",
        #"timemodel",
        "multiband_photometry",
        "shape_fit_multiband",
        "spectrum",
    }
    catalogue = pd.read_csv(MEDELLIN_CATALOG).dropna(subset=["provid"])
    if "band_periodogram" not in catalogue.columns:
        catalogue["band_periodogram"] = BAND
    provids = list(
        zip(
            catalogue["provid"].astype(str).str.strip(),
            catalogue["band_periodogram"].fillna(BAND).astype(str).str.strip(),
        )
    )
    selected, df_set, df_rot, analysis, missing = _selection_tables(provids, analysis)
    selected.attrs["Nb"] = int(Nb)
    return selected, df_set, df_rot, analysis, missing


def print_processing_report(total_requested, succeeded, failed, fig_dir=None):
    """Print and save a summary report showing success and failure due to missing data."""
    n_succeeded = len(succeeded)
    n_failed = len(failed)
    total = total_requested
    pct_ok = (n_succeeded / total * 100.0) if total > 0 else 0.0
    pct_fail = (n_failed / total * 100.0) if total > 0 else 0.0

    lines = []
    lines.append("")
    lines.append("=" * 80)
    lines.append("              INFORME FINAL DE PROCESAMIENTO DE ASTEROIDES")
    lines.append("=" * 80)
    lines.append(f"Total de asteroides considerados : {total}")
    lines.append(f"  ✓ Cálculos realizados con éxito: {n_succeeded} ({pct_ok:.1f}%)")
    lines.append(f"  ✗ Problemas por falta de datos : {n_failed} ({pct_fail:.1f}%)")
    lines.append("-" * 80)

    if failed:
        lines.append(f"DETALLE DE ASTEROIDES CON PROBLEMAS ({n_failed}):")
        reasons_map = {}
        for provid, reason in failed.items():
            reasons_map.setdefault(reason, []).append(provid)

        for reason, provids in sorted(
            reasons_map.items(), key=lambda kv: len(kv[1]), reverse=True
        ):
            lines.append(f"\n  • Motivo: {reason} ({len(provids)} asteroides)")
            lines.append(repr(sorted(provids)))
        lines.append("-" * 80)

    lines.append("=" * 80)
    lines.append("")
    print("\n".join(lines))

    if fig_dir is not None:
        file_lines = lines[:-2]
        file_lines.append(f"ASTEROIDES CON CÁLCULO EXITOSO ({n_succeeded}):")
        file_lines.append(repr(sorted(succeeded)))
        file_lines.append("=" * 80)
        file_lines.append("")
        try:
            report_path = Path(fig_dir) / "asteroid_processing_report.txt"
            report_path.write_text("\n".join(file_lines), encoding="utf-8")
        except Exception as e:
            print(f"Warning: could not write processing report: {e}")


def _cli_provids(argv):
    """Asteroids named with --provid. Underscores in the name become spaces."""
    provids = []
    args = [arg for arg in argv[1:] if arg not in ("--force", "-f")]
    i = 0
    while i < len(args):
        arg = args[i]
        if arg == "--provid":
            if i + 1 >= len(args):
                raise SystemExit('Missing asteroid name after --provid.')
            provids.append(args[i + 1].replace("_", " ").strip())
            i += 2
            continue
        if arg.startswith("--provid="):
            name = arg.split("=", 1)[1].replace("_", " ").strip()
            if not name:
                raise SystemExit('Missing asteroid name after --provid.')
            provids.append(name)
            i += 1
            continue
        if arg.startswith("-"):
            raise SystemExit(f"Unknown option: {arg}")
        provids.append(arg.replace("_", " ").strip())
        i += 1
    return [name for name in provids if name]


def _main(t_start):
    force = "--force" in sys.argv or "-f" in sys.argv
    cli_provids = _cli_provids(sys.argv)

    if cli_provids:
        selected, df_set, df_rot, analysis, missing = unpack_selection(
            select_provid_list(provids=cli_provids)
        )
    else:
        # select_provid_all: more than MIN_N_OBS observations
        # select_provid_list: the list written inside that function
        # select_provid_rotation: asteroids-with-rotation.csv
        # select_provid_medellin: medellin-vro-catalogue.csv
        # select_provid_nobs: 50 to 100 observations
        # selected, df_set, df_rot, analysis, missing = unpack_selection(select_provid_all())
        # selected, df_set, df_rot, analysis, missing = unpack_selection(select_provid_list())
        selected, df_set, df_rot, analysis, missing = unpack_selection(
            # select_provid_rotation()
            # select_provid_list()
            select_provid_medellin()
            # select_provid_nobs()
        )
        # selected, df_set, df_rot, analysis, missing = unpack_selection(select_provid_medellin())
    designations = load_mpc_designations()

    total_requested = len(selected) + len(missing)
    succeeded = []
    failed = {
        p: "Sin observaciones fotométricas en obs_sbn_X05_full.parquet"
        for p in missing
    }

    peak_rows = []
    replaced = set()
    exclude_lx = tuple(selected.attrs.get("exclude_lx", EXCLUDE_LX))
    nb = int(selected.attrs.get("Nb", 2))
    bands = (
        selected["band"]
        if "band" in selected.columns
        else pd.Series(BAND, index=selected.index)
    )
    n_asteroids = len(selected)
    total_obs = int(selected["n_obs"].sum())
    t_loop = time.perf_counter()
    done_obs = 0
    for i, (provid, n_obs, band) in enumerate(
        zip(selected["provid"], selected["n_obs"], bands), start=1
    ):
        print(f"N = {int(n_obs)}, band = {band}")
        fit, err = process(
            df_set,
            provid,
            df_rot,
            analysis,
            designations,
            force=force,
            band=band,
            exclude_lx=exclude_lx,
            Nb=nb,
        )
        if fit == "skip":
            succeeded.append(provid)
        else:
            replaced.add(provid)
            if fit is not None:
                succeeded.append(provid)

                df_ast = df_set[df_set["provid"] == provid]
                t_min = df_ast["obstime"].min()
                t_max = df_ast["obstime"].max()
                epoch_utc_ini = f"{t_min:%Y-%m-%d %H:%M:%S}" if pd.notna(t_min) else ""
                epoch_utc_end = f"{t_max:%Y-%m-%d %H:%M:%S}" if pd.notna(t_max) else ""
                n_obs_ast = int(len(df_ast))
                n_obs_bands = {
                    f"n_obs_{b}": int((df_ast["band"].isin([f"L{b}", b])).sum())
                    for b in ["u", "g", "r", "i", "z", "y"]
                }

                row = {
                    "provid": provid,
                    "epoch_utc_ini": epoch_utc_ini,
                    "epoch_utc_end": epoch_utc_end,
                    "n_obs": n_obs_ast,
                    **n_obs_bands,
                }

                if "periodogram" in analysis:
                    periods_h = (
                        ",".join(
                            f"{period:.6f}"
                            for period in fit.get("df_peaks", {}).get(
                                "rotational_period_h", []
                            )
                        )
                        if "df_peaks" in fit
                        else ""
                    )
                    row["band_periodogram"] = fit.get("band", band)
                    row["periods_h"] = periods_h
                    row["power"] = fit.get("power", np.nan)
                    row["amp_lsm"] = fit.get("amp_lsm", np.nan)
                    row["amp_photo"] = fit.get("amp_photo", np.nan)

                if "periodogram_democratic" in analysis and fit.get("democratic"):
                    dem = fit["democratic"]
                    row["p_rot_dem_h"] = dem["p_rot_h"]
                    row["votes_dem"] = dem["votes"]
                    row["n_series_dem"] = dem["n_series"]
                    row["bands_dem"] = dem["bands"]
                    row["periods_dem_h"] = dem["periods_h"]

                if "shape_fit" in analysis:
                    row.update(fit.get("shape", {}))

                if "multiband_photometry" in analysis:
                    if "multiband" in fit and fit["multiband"] is not None:
                        row.update(fit["multiband"].get("summary_row", {}))

                if "shape_fit_multiband" in analysis:
                    if "shape_k" in fit:
                        row.update(fit["shape_k"])

                if "spectrum" in analysis:
                    row.update(fit.get("spectrum") or {})

                peak_rows.append(row)
            else:
                failed[provid] = err or "Datos insuficientes o error en el cálculo"

        done_obs += int(n_obs)
        elapsed = time.perf_counter() - t_start
        # Scale the asteroid loop by observations; the catalog load is paid once.
        estimated_total = (t_loop - t_start) + (time.perf_counter() - t_loop) * (
            total_obs / done_obs
        )
        remaining = max(0.0, estimated_total - elapsed)
        print(
            f"  elapsed: {format_elapsed(elapsed)} ({i}/{n_asteroids}), "
            f"estimated total: {format_elapsed(estimated_total)}, "
            f"time remaining (est.): {format_elapsed(remaining)}"
        )

    write_analysis_results(peak_rows, replaced, analysis)
    copy_analysis_docs()

    if "compare" in analysis:
        print("Period comparison with the literature")
        plot_period_comparison(df_rot)
        print(f"  elapsed: {format_elapsed(time.perf_counter() - t_start)}")

    print_processing_report(total_requested, succeeded, failed, fig_dir=FIG_DIR)


def main():
    t_start = time.perf_counter()
    try:
        _main(t_start)
    finally:
        _stop_log()


if __name__ == "__main__":
    main()
