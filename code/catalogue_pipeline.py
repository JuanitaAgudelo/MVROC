#!/usr/bin/env python
"""Repeat the rotation figures of search_spectra_rotation_orbit.ipynb.

Which asteroids and which figures to write is set in select_provid_all(),
select_provid_list(), select_provid_rotation(), select_provid_rotfold(),
select_provid_medellin(), select_provid_nobs(), or select_provid_revisit().
The active one is the call at the start of main().

Writes, into figall/, as WebP at half the previous resolution:
  <provid>_lightcurve_allnights_H.webp
  <provid>_periodogram.webp
  <provid>_periodogram_democratic.webp
  <provid>_phase_periodogram.webp
  <provid>_phase_democratic.webp
  <provid>_phase_democratic_weighted.webp
  <provid>_timemodel_H.webp
  <provid>_timemodel_Lk.webp
  <provid>_phase_model.webp
  <provid>_phase_model_compare.webp
  <provid>_multiband_photometry_compare.webp
  <provid>_multiband_photometry_nightalias.webp
  <provid>_justfold.webp
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

import os
import re
import shutil
import sys
import time
from datetime import datetime, timezone
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
REVISIT_CATALOG = DATA_DIR / "revisit.csv"
DESIGNATIONS = DATA_DIR / "mpc_designations.csv"
OBS_FILE = DATA_DIR / "obs_sbn_X05_full_mpc_desig.parquet"

BAND = "Lr"
PROP = "H"
PERIODS_H = (0.5 / 60, 24)  # photometric periods in hours; rotational plot spans 1 min to 48 h
ALIAS_P_MIN_H = 2.0  # shortest rotation period of a campaign-duration alias, in hours
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
    "periodogram_democratic_weighted",
    "phase",
    "phase_model",
    "phase_model_compare",
    "multiband_photometry_compare",
    "multiband_photometry_nightalias",
    "justfold",
    "shape_fit",
    "timemodel",
    "timemodel_Lk",
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
    """Copy stdout into path. If the file already exists, append a new run."""
    global _LOG_FILE, _PREV_STDOUT
    _stop_log()
    path = Path(path)
    append = path.is_file() and path.stat().st_size > 0
    _LOG_FILE = open(path, "a" if append else "w", encoding="utf-8")
    _PREV_STDOUT = sys.stdout
    sys.stdout = _Tee(_PREV_STDOUT, _LOG_FILE)
    if append:
        stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        print(f"\n{'=' * 72}")
        print(f"New run started: {stamp}")
        print(f"{'=' * 72}")


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
        "periodogram_democratic": "periodogram_democratic",
        "phase": "phase_allbands",
        "phase_periodogram": "phase_periodogram",
        "phase_democratic": "phase_democratic",
        "phase_democratic_weighted": "phase_democratic_weighted",
        "timemodel": "timemodel_H",
        "timemodel_Lk": "timemodel_Lk",
        "phase_model": "phase_model",
        "phase_model_compare": "phase_model_compare",
        "multiband_photometry_compare": "multiband_photometry_compare",
        "multiband_photometry_nightalias": "multiband_photometry_nightalias",
        "justfold": "justfold",
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


def mvroc_mark(ax, fontsize=12, alpha=0.55):
    """Add an MVROC watermark at about the size of the figure title."""
    args = dict(
        rotation=270,
        ha="left",
        va="top",
        transform=ax.transAxes,
        color="pink",
        alpha=alpha,
        fontsize=fontsize,
        zorder=100,
    )
    try:
        ax.add_collection3d
        plt_text = ax.text2D
    except Exception:
        plt_text = ax.text
    return plt_text(1, 1, "MVROC", **args)


def add_mvroc_watermarks(fig):
    """Stamp MVROC on every main axes (skip insets and colorbars)."""
    inset_ids = {
        id(child)
        for ax in fig.axes
        for child in getattr(ax, "child_axes", ())
    }
    for ax in fig.axes:
        if id(ax) in inset_ids:
            continue
        if ax.get_label() == "<colorbar>":
            continue
        mvroc_mark(ax)


def save_figure(fig, path):
    """Write a WebP at DPI. The path suffix is forced to .webp."""
    path = Path(path).with_suffix(".webp")
    # Draw once so axis window extents are known before sizing the mark.
    fig.canvas.draw()
    add_mvroc_watermarks(fig)
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
    has_compare_period = (
        "multiband_photometry_compare" in analysis
        or "multiband_photometry_nightalias" in analysis
    )
    if "shape_fit_multiband" in analysis and not (
        "multiband_photometry" in analysis or has_compare_period
    ):
        raise SystemExit(
            f"Error in {where}: 'shape_fit_multiband' needs "
            "'multiband_photometry_compare', 'multiband_photometry_nightalias', "
            "or 'multiband_photometry'."
        )
    if "spectrum" in analysis and not (
        "multiband_photometry" in analysis or has_compare_period
    ):
        raise SystemExit(
            f"Error in {where}: 'spectrum' needs "
            "'multiband_photometry_compare', 'multiband_photometry_nightalias', "
            "or 'multiband_photometry'."
        )
    if "timemodel_Lk" in analysis and not has_compare_period:
        raise SystemExit(
            f"Error in {where}: 'timemodel_Lk' needs 'multiband_photometry_compare' "
            "or 'multiband_photometry_nightalias'."
        )
    if "phase_model_compare" in analysis and not (
        "periodogram_democratic" in analysis
        or "periodogram_democratic_weighted" in analysis
    ):
        raise SystemExit(
            f"Error in {where}: 'phase_model_compare' needs 'periodogram_democratic' "
            "or 'periodogram_democratic_weighted'."
        )
    if has_compare_period and "periodogram_democratic" not in analysis:
        raise SystemExit(
            f"Error in {where}: 'multiband_photometry_compare' and "
            "'multiband_photometry_nightalias' need 'periodogram_democratic'."
        )


check_figures = check_analysis


def load_df_set():
    """Photometry table, filtered as in the notebook (df_set)."""
    columns = [
        "obstime",
        "permid",
        "provid",
        "ra",
        "dec",
        "band",
        "mag",
        "rmsmag",
        "mpc_designation",
    ]
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


def observation_window(df_ast):
    """Start, end, and length of one asteroid's photometric series.

    The length is the time from the first measurement to the last, in days.
    """
    t_min = df_ast["obstime"].min()
    t_max = df_ast["obstime"].max()
    ini = f"{t_min:%Y-%m-%d %H:%M:%S}" if pd.notna(t_min) else ""
    end = f"{t_max:%Y-%m-%d %H:%M:%S}" if pd.notna(t_max) else ""
    if pd.notna(t_min) and pd.notna(t_max):
        duration = float((t_max - t_min).total_seconds() / 86400.0)
    else:
        duration = np.nan
    return ini, end, duration


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


def _object_designations(df_obj, provid, designations):
    """MPC designations for one provid.

    The observation catalogue carries the designation used by the geometry
    table. The designation file is only a fallback when that column is empty.
    """
    found = []
    if "mpc_designation" in df_obj.columns:
        values = df_obj["mpc_designation"].dropna().astype(str).str.strip()
        values = values[~values.str.lower().isin(["", "nan", "none"])]
        found.extend(values.unique().tolist())
    fallback = None if designations is None else designations.get(provid)
    if fallback is not None and not pd.isna(fallback):
        text = str(fallback).strip()
        if text and text.lower() not in ("nan", "none") and text not in found:
            found.append(text)
    return found


def prepare_object(df_set, provid, designations):
    df_obj = df_set[df_set["provid"] == provid].copy()
    if df_obj.empty:
        return df_obj
    df_obj = df_obj.sort_values(["band", "obstime"]).reset_index(drop=True)
    t0 = df_obj["obstime"].iloc[0]
    df_obj["t_days"] = (df_obj["obstime"] - t0).dt.total_seconds() / 86400.0
    frames = [load_geometry(name) for name in _object_designations(df_obj, provid, designations)]
    frames = [frame for frame in frames if not frame.empty]
    if frames:
        df_spy = pd.concat(frames, ignore_index=True).drop_duplicates("obstime")
    else:
        df_spy = pd.DataFrame(columns=["obstime", "mag", "H"])
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
    ax.set_ylabel("Absolute magnitude, H")
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
        ax.set_ylabel("Absolute magnitude, H")
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
        t_campaign = pd.to_numeric(df_obj["t_days"], errors="coerce").to_numpy(float)
        t_campaign = t_campaign[np.isfinite(t_campaign)]
        span = float(np.ptp(t_campaign)) if len(t_campaign) else 0.0
        _draw_campaign_alias_vlines(ax, span)
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


def _multiband_combos(bands):
    """Every pair of bands, then all of them when more than two are present."""
    bands = [str(band) for band in bands]
    combos = [
        [bands[i], bands[j]]
        for i in range(len(bands))
        for j in range(i + 1, len(bands))
    ]
    if len(bands) > 2:
        combos.append(list(bands))
    return combos


def _draw_campaign_alias_vlines(ax, span_days, tol=0.1):
    """Shaded 2-span campaign aliases from the 48 h edge down to ALIAS_P_MIN_H.

    Only the 2/T + k family is drawn. Each band is an axvspan whose width in
    period matches the frequency tolerance used to tag a voted peak as an
    alias: periods between 48/(f+tol) and 48/(f-tol).
    """
    span = float(span_days)
    if not np.isfinite(span) or span <= 0 or tol <= 0:
        return
    k_values = list(_campaign_alias_k_range(span))
    k_first = k_values[0] if k_values else None
    base = 2.0 / span
    for k in k_values:
        freq = base + float(k)
        f_lo = freq - float(tol)
        f_hi = freq + float(tol)
        if f_lo <= 0:
            continue
        # Lower frequency -> longer rotation period.
        p_lo = 48.0 / f_hi
        p_hi = 48.0 / f_lo
        ax.axvspan(
            p_lo,
            p_hi,
            color="C4",
            alpha=0.25,
            zorder=1,
            linewidth=0,
            label=(
                rf"Campaign aliases ($\pm${tol:.1f})"
                if k == k_first
                else None
            ),
        )


def periodogram_democratic(df_obj, provid, Nb=3, exclude=EXCLUDE_LX):
    """Vote a rotation period from the Nb bands with the most points.

    Each of those bands votes with its own Lomb-Scargle maximum. Every pair
    of those bands, and all of them together, each cast one more vote.
    u and y are excluded. Periods within 3 percent, or one or two cycles per
    day apart, count as the same vote. A unanimous vote uses the period of
    all the bands together. Otherwise the reported period is the Lomb-Scargle
    period of the band with the most points among the series that cast the
    winning vote. A periodogram of the best-sampled band is written with the
    campaign-alias lines and the voted periods marked.
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
        combos = _multiband_combos(bands)
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

    # Periodogram of the best-sampled band, with campaign aliases and votes.
    band_plot = str(usable.index[0])
    df_plot = df_obj[(df_obj["band"] == band_plot) & df_obj[PROP].notna()]
    t_plot = df_plot["t_days"].to_numpy(float)
    y_plot = df_plot[PROP].to_numpy(float)
    dy_plot = df_plot["rmsmag"].fillna(0.01).to_numpy(float)
    power_plot = np.asarray(
        LombScargle(t_plot, y_plot, dy_plot).power(freq_grid), dtype=float
    )
    p_phot_min = cadence_min_period(t_plot)
    period_h = 48.0 / freq_grid
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(period_h, power_plot, color="C0", lw=1.0)
    t_campaign = pd.to_numeric(df_obj["t_days"], errors="coerce").to_numpy(float)
    t_campaign = t_campaign[np.isfinite(t_campaign)]
    span = float(np.ptp(t_campaign)) if len(t_campaign) else 0.0
    _draw_campaign_alias_vlines(ax, span)
    # One line per distinct voted frequency (merged within 0.1 cycles/day).
    voted_merged = _merge_frequencies(
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
    ax.set_ylabel("Lomb–Scargle power")
    ax.set_title(f"{provid} democratic periodogram ({band_plot})")
    ax.legend(loc="upper left", fontsize=7, frameon=True)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    path = save_figure(fig, output_path("periodogram_democratic", provid))
    print(f"  democratic periodogram: {path.name}")

    return {
        "p_rot_h": adopted,
        "votes": int(len(winner_ids)),
        "n_series": int(len(series)),
        "bands": "+".join(str(band) for band in usable.index),
        "series": ", ".join(labels),
        "source": series[chosen]["label"],
        "periods_h": periods_h,
        "path": path,
        "candidates": [
            {
                "label": item["label"],
                "p_rot_h": float(item["p_rot_h"]),
                "freq": 48.0 / float(item["p_rot_h"]),
            }
            for item in series
        ],
    }


def _band_power(df_i, freq_grid, bands=None):
    """Lomb-Scargle power of one series on the shared frequency grid.

    A multiband series passes the band label of each row. A single band does not.
    """
    t = df_i["t_days"].to_numpy(float)
    y = df_i[PROP].to_numpy(float)
    dy = np.maximum(df_i["rmsmag"].fillna(0.01).to_numpy(float), 1e-4)
    if bands is None:
        ls = LombScargle(t, y, dy)
    else:
        ls = LombScargleMultiband(t, y, np.asarray(bands), dy)
    power = np.asarray(ls.power(freq_grid), dtype=float)
    p_phot_min = cadence_min_period(t)
    allowed = (24.0 / freq_grid) >= p_phot_min if p_phot_min > 0 else np.ones(len(freq_grid), bool)
    if not np.any(allowed):
        allowed = np.ones(len(freq_grid), dtype=bool)
        p_phot_min = 0.0
    best_i = int(np.flatnonzero(allowed)[np.argmax(power[allowed])])
    return {
        "ls": ls,
        "t": t,
        "y": y,
        "dy": dy,
        "freq": float(freq_grid[best_i]),
        "p_rot_h": 48.0 / float(freq_grid[best_i]),
        "power_max": float(power[best_i]),
        "cadence_h": 2.0 * p_phot_min,
        "n": int(len(df_i)),
    }


def _power_at(ls, freq):
    return float(np.asarray(ls.power(float(freq)), dtype=float).reshape(-1)[0])


def periodogram_democratic_weighted(df_obj, provid, Nb=3, exclude=EXCLUDE_LX):
    """Rotation period from the daily aliases of the best-sampled band.

    The reference band is the one with the most measurements. Candidates are
    its Lomb-Scargle frequency shifted by k cycles per day, for k from -6 to
    6.     The series are the same as in the unweighted vote: each band, every pair,
    and all of the bands together. Each series scores the
    candidates with its relative power, weighted by its number of points. A
    near tie is broken by the multiband periodogram of the single bands.
    """
    counts = (
        df_obj.dropna(subset=[PROP])
        .loc[~df_obj["band"].isin(exclude), "band"]
        .value_counts()
    )
    usable = counts[counts >= 5].head(int(Nb))
    if usable.empty:
        print(f"  skip weighted periodogram for {provid}: no band with 5 points")
        return None

    frames = {}
    for band_i in usable.index:
        frames[str(band_i)] = df_obj[(df_obj["band"] == band_i) & df_obj[PROP].notna()]
    df_top = pd.concat(frames.values(), ignore_index=True)
    freq_grid = frequency_grid(df_top["t_days"].to_numpy())
    scans = {band: _band_power(frame, freq_grid) for band, frame in frames.items()}
    band_names = list(scans)
    if len(band_names) >= 2:
        for combo in _multiband_combos(band_names):
            label = "+".join(combo)
            frame = df_obj[df_obj["band"].isin(combo) & df_obj[PROP].notna()]
            scans[label] = _band_power(
                frame, freq_grid, bands=frame["band"].to_numpy()
            )

    reference = str(usable.index[0])
    ref = scans[reference]
    f_ref = ref["freq"]
    f_min = 24.0 / PERIODS_H[1]
    f_max = 24.0 / PERIODS_H[0]
    p_phot_min = ref["cadence_h"] / 2.0
    n_total = sum(scan["n"] for scan in scans.values())
    print(
        f"  weighted reference {reference}: N = {ref['n']}, "
        f"P_rot = {ref['p_rot_h']:.4f} h, power = {ref['power_max']:.3f}"
    )

    scores = {}
    for k in range(-6, 7):
        f_k = f_ref + k
        if f_k <= 0 or f_k < f_min or f_k > f_max:
            continue
        if p_phot_min > 0 and (24.0 / f_k) < p_phot_min:
            continue
        score = 0.0
        for scan in scans.values():
            if scan["power_max"] <= 0:
                continue
            rho = _power_at(scan["ls"], f_k) / scan["power_max"]
            score += (scan["n"] / n_total) * rho
        scores[k] = score

    if not scores:
        print(f"  skip weighted periodogram for {provid}: no alias inside the search window")
        return None

    best_score = max(scores.values())
    tied = [k for k, score in scores.items() if score >= (1.0 - 0.05) * best_score]
    single = {band: scan for band, scan in scans.items() if "+" not in band}
    if len(tied) > 1 and len(single) > 1:
        bands = np.concatenate(
            [np.full(scan["n"], band) for band, scan in single.items()]
        )
        ls_multi = LombScargleMultiband(
            np.concatenate([scan["t"] for scan in single.values()]),
            np.concatenate([scan["y"] for scan in single.values()]),
            bands,
            np.concatenate([scan["dy"] for scan in single.values()]),
        )
        k_star = max(tied, key=lambda k: _power_at(ls_multi, f_ref + k))
        print(
            f"  weighted tie among k = {tied}; multiband power selects k = {k_star}"
        )
    else:
        k_star = max(scores, key=scores.get)

    f_star = f_ref + k_star
    p_star = 48.0 / f_star
    candidates = []
    for band, scan in scans.items():
        candidates.append(
            {
                "label": band,
                "p_rot_h": float(scan["p_rot_h"]),
                "freq": float(scan["freq"]),
            }
        )
        if scan["power_max"] <= 0:
            continue
        prefer = max(
            scores,
            key=lambda k, scan=scan: _power_at(scan["ls"], f_ref + k) / scan["power_max"],
        )
        rho_star = _power_at(scan["ls"], f_star) / scan["power_max"]
        print(
            f"  {band}: N = {scan['n']}, P_rot = {scan['p_rot_h']:.2f} h, "
            f"power = {scan['power_max']:.3f}, prefers k = {prefer:+d}, "
            f"rho(k*) = {rho_star:.2f}"
        )
    score_text = ",".join(f"k={k:+d}:{scores[k]:.3f}" for k in sorted(scores))
    print(
        f"  weighted (Nb = {int(Nb)}): P_rot = {p_star:.4f} h from {reference}, "
        f"k = {k_star:+d}, S = {scores[k_star]:.3f}"
    )
    print(f"  alias scores: {score_text}")
    return {
        "p_rot_h": float(p_star),
        "k_alias": int(k_star),
        "score": float(scores[k_star]),
        "band_ref": reference,
        "n_bands": int(len(single)),
        "bands": "+".join(single),
        "scores": score_text,
        "source": reference,
        "candidates": candidates,
    }


def _magnitude_inliers(values, n_iter=5):
    """True for magnitudes that belong to the main light curve of one band.

    The cut is the 5th to 95th percentile, extended by a margin of at least
    0.2 mag. That keeps the rotational maxima and minima. Points that sit
    well outside that interval, like the faint measurements that stretch the
    phase plot, are rejected.
    """
    values = np.asarray(values, dtype=float)
    good = np.isfinite(values)
    if int(good.sum()) < 8:
        return good
    for _ in range(n_iter):
        sample = values[good]
        low, high = np.percentile(sample, [5, 95])
        margin = max(0.2, 0.35 * float(high - low))
        updated = np.isfinite(values) & (values >= low - margin) & (values <= high + margin)
        if int(updated.sum()) < 8 or np.array_equal(updated, good):
            break
        good = updated
    return good


def plot_phase_allbands(
    df_obj, provid, fit, kind="phase_periodogram", freq_sel=None, per_sel=None, source=None
):
    path = output_path(kind, provid)
    if path.exists():
        return path
    t0 = fit["t_band"].min()
    if freq_sel is None:
        freq_sel = fit["freq_sel"]
    if per_sel is None:
        per_sel = fit["per_sel"]
    if source is None:
        source = fit.get("period_source") or f"Periodogram {fit.get('band', BAND)}"
    present = [
        b for b in ALL_BANDS if b not in ("u", "Lu") and (df_obj["band"] == b).any()
    ]

    fig, ax = plt.subplots(figsize=(8, 5))
    n_rejected = 0
    for b in present:
        sb = df_obj[df_obj["band"] == b].dropna(subset=[PROP])
        good = _magnitude_inliers(sb[PROP].to_numpy(float))
        n_rejected += int((~good).sum())
        sb = sb.iloc[np.flatnonzero(good)]
        if sb.empty:
            continue
        phase_b = np.mod((sb["t_days"].to_numpy() - t0) * freq_sel / 2.0, 1.0)
        dy = sb["rmsmag"].fillna(0.01)
        marker = "o" if b.startswith("L") else "s"
        m = float(sb[PROP].min())
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
    if n_rejected:
        print(f"  phase: rejected {n_rejected} magnitude outlier(s)")

    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Phase at the rotation period")
    ax.set_ylabel(r"$m - m_{\mathrm{max}}$")
    ax.set_title(
        f"Phase curve of {provid}, $P_{{rot}}$ = {per_sel:.2f} h\n({source})"
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
        ax.set_ylabel("Absolute magnitude, H")
        ax.grid(True, alpha=0.3)

    axs[0].invert_yaxis()
    axs[0].legend(loc="best", frameon=True)
    axs[-1].set_xlabel("Hours since the first observation of the night")
    fig.suptitle(
        f"{provid}: observed light curve and periodogram model ({fit.get('band', BAND)})"
    )
    fig.tight_layout()
    return save_figure(fig, output_path("timemodel", provid))


def plot_timemodel_lk(provid, df_band_k, p_rot_h):
    """Night-by-night Lk light curve and the two-peak model at one rotation period."""
    path = output_path("timemodel_Lk", provid)
    if path.exists():
        return path
    df_band = df_band_k.dropna(subset=[PROP, "t_days", "night"]).copy()
    if df_band.empty:
        print(f"  skip timemodel_Lk for {provid}: no Lk points")
        return None
    nights_band = df_band["night"].drop_duplicates().sort_values()
    n_nights_band = len(nights_band)
    t_all = df_band["t_days"].to_numpy(float)
    y_all = df_band[PROP].to_numpy(float)
    dy_all = np.maximum(
        pd.to_numeric(df_band["rmsmag"], errors="coerce").fillna(0.01).to_numpy(float),
        1e-4,
    )
    f_rot = 24.0 / float(p_rot_h)
    ls_k = LombScargle(t_all, y_all, dy_all, nterms=2)
    fig, axs = plt.subplots(
        n_nights_band, 1, figsize=(10, 2.6 * n_nights_band), sharey=True
    )
    if n_nights_band == 1:
        axs = [axs]

    for ax, night in zip(axs, nights_band):
        sel = (df_band["night"] == night).to_numpy()
        t_night = t_all[sel]
        y_night = y_all[sel]
        dy_night = dy_all[sel]
        t_left, t_right = t_night.min(), t_night.max()
        pad = max(0.01, 0.05 * (t_right - t_left))
        t_grid = np.linspace(t_left - pad, t_right + pad, 400)
        y_grid = ls_k.model(t_grid, f_rot)
        hours = (t_night - t_left) * 24.0
        hours_grid = (t_grid - t_left) * 24.0
        ax.errorbar(
            hours,
            y_night,
            yerr=dy_night,
            fmt="o",
            ms=4,
            alpha=0.85,
            color=band_colors.get("Lk", "0.2"),
            ecolor="0.6",
            label="Lk magnitude",
            zorder=3,
        )
        ax.plot(
            hours_grid,
            y_grid,
            color="C3",
            lw=1.8,
            label=rf"2-peak model, $P_{{rot}}$ = {p_rot_h:.2f} h",
            zorder=2,
        )
        ax.set_title(pd.to_datetime(night).strftime("%Y-%m-%d"))
        ax.set_ylabel("Absolute magnitude, H")
        ax.grid(True, alpha=0.3)

    axs[0].invert_yaxis()
    axs[0].legend(loc="best", frameon=True)
    axs[-1].set_xlabel("Hours since the first observation of the night")
    fig.suptitle(
        f"{provid}: Lk light curve and 2-peak model, "
        rf"$P_{{rot}}$ = {p_rot_h:.2f} h"
    )
    fig.tight_layout()
    return save_figure(fig, path)


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
    ax.set_ylabel("Absolute magnitude, H")
    ax.set_title(
        f"Folded light curve of {provid} ({PROP},{band}), "
        rf"$P_{{rot}}$ = {per_rot:.2f} h"
    )
    fig.tight_layout()
    return save_figure(fig, path)


def _merge_frequencies(items, tol=0.1):
    """Group (frequency, label) pairs that agree within tol cycles/day."""
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


def _vote_frequencies(periodogram_fold, dem, wdem):
    """Photometric frequencies reported by the periodogram and by each vote."""
    items = [
        (periodogram_fold["freq_sel"], periodogram_fold["source"]),
    ]
    for vote in (dem, wdem):
        if not vote:
            continue
        for cand in vote.get("candidates", []):
            items.append((cand["freq"], cand["label"]))
    return _merge_frequencies(items)


def phase_model_compare(provid, fit, candidates):
    """Fold each voted frequency and keep the period whose model has the lowest chi-squared.

    The model is the Lomb-Scargle fit of the periodogram band, evaluated at
    that frequency, the same comparison drawn by the folded light curve.
    """
    t_band = fit["t_band"]
    y_band = fit["y_band"]
    dy_band = np.maximum(np.asarray(fit["dy_band"], dtype=float), 1e-4)
    ls_band = fit["ls_band"]
    band = fit.get("band", BAND)
    nterms = int(getattr(ls_band, "nterms", 1))
    dof = max(int(len(y_band)) - (1 + 2 * nterms), 1)
    t0 = float(np.min(t_band))
    scored = []
    for cand in candidates:
        freq = float(cand["freq"])
        model = np.asarray(ls_band.model(t_band, freq), dtype=float)
        chi2 = float(np.sum(((y_band - model) / dy_band) ** 2))
        scored.append(
            {
                "freq": freq,
                "p_rot_h": 48.0 / freq,
                "labels": list(cand["labels"]),
                "chi2": chi2,
                "chi2_nu": chi2 / dof,
            }
        )
    scored.sort(key=lambda item: (item["chi2_nu"], item["p_rot_h"]))
    best = scored[0]
    periods = ",".join(
        f"{item['p_rot_h']:.4f}:{item['chi2_nu']:.2f}" for item in scored
    )
    for item in scored:
        mark = " selected" if item is best else ""
        print(
            f"  {item['p_rot_h']:.4f} h ({', '.join(item['labels'])}): "
            f"chi2_nu = {item['chi2_nu']:.2f}{mark}"
        )
    print(
        f"  phase-model compare: P_rot = {best['p_rot_h']:.4f} h, "
        f"chi2_nu = {best['chi2_nu']:.2f} (best of {len(scored)})"
    )

    path = output_path("phase_model_compare", provid)
    if not path.exists():
        n = len(scored)
        ncol = 2 if n > 1 else 1
        nrow = int(np.ceil(n / ncol))
        fig, axs = plt.subplots(
            nrow, ncol, figsize=(11, 3.4 * nrow), sharex=True, sharey=True, squeeze=False
        )
        axes = axs.ravel()
        phase_grid = np.linspace(0, 1, 400)
        curve_color = band_colors.get(band, "k")
        for ax, item in zip(axes, scored):
            freq = item["freq"]
            phase = np.mod((t_band - t0) * freq / 2.0, 1.0)
            y_model = ls_band.model(t0 + phase_grid * 2.0 / freq, freq)
            ax.errorbar(
                phase,
                y_band,
                yerr=dy_band,
                fmt="o",
                ms=3,
                alpha=0.7,
                color=curve_color,
                ecolor="0.7",
                label="Magnitude",
            )
            ax.plot(
                phase_grid,
                y_model,
                color="0.15",
                lw=1.8,
                label="Periodogram model",
            )
            mark = ", selected" if item is best else ""
            ax.set_title(
                rf"$P_{{rot}}$ = {item['p_rot_h']:.2f} h, "
                rf"$\chi^2_\nu$ = {item['chi2_nu']:.1f}{mark}"
                "\n"
                + ", ".join(item["labels"]),
                fontsize=9,
            )
            ax.grid(True, alpha=0.3)
            ax.set_xlim(0, 1)
        axes[0].invert_yaxis()
        axes[0].legend(loc="best", frameon=True, fontsize=8)
        for ax in axes[n:]:
            ax.set_visible(False)
        for ax in axes[:n]:
            if ax in axs[-1]:
                ax.set_xlabel("Phase at the rotation period ($2P$)")
        for ax in axs[:, 0]:
            ax.set_ylabel("Absolute magnitude, H")
        fig.suptitle(
            f"Phase-model comparison of {provid} ({PROP},{band})",
            fontsize=12,
        )
        fig.tight_layout()
        save_figure(fig, path)
    return {
        "p_rot_h": float(best["p_rot_h"]),
        "chi2_nu": float(best["chi2_nu"]),
        "periods": periods,
        "path": path,
    }


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
        ax.set_ylabel("Absolute magnitude, H")
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
    ("p_rot_wdem_h", "democratic_weighted", "Democratic, weighted"),
    ("p_rot_h", "multiband_compare", "Multiband 2-peak compare"),
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
    ax.set_ylabel("Absolute magnitude, H")
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
    ax.set_ylabel("Absolute magnitude, H")
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
    ax.set_ylabel("Absolute magnitude, H")
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


def _figure_paths(kind, provid, analysis=None):
    if kind == "shape_fit":
        return [output_path("spheres", provid), output_path("ellipsoid", provid)]
    if kind == "multiband_photometry":
        return [
            # output_path("multiband_lightcurve", provid),
            output_path("multiband_phase_model", provid),
            FIG_DIR / f"{slug(provid)}_observations.csv",
        ]
    if kind == "shape_fit_multiband":
        return [output_path("shape_multiband", provid)]
    if kind in ("periodogram_democratic", "periodogram_democratic_weighted"):
        return []
    if kind == "phase":
        paths = [output_path("phase_periodogram", provid)]
        selected = set(analysis or ())
        if "periodogram_democratic" in selected:
            paths.append(output_path("phase_democratic", provid))
        if "periodogram_democratic_weighted" in selected:
            paths.append(output_path("phase_democratic_weighted", provid))
        return paths
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


def _brightness_peaks(y_model):
    """Number of brightness maxima of one rotation. Magnitude minima are the peaks.

    The phase is circular, so the curve is repeated and only the middle copy is
    counted. A maximum has to rise at least 0.02 mag, and at least 15% of the
    model amplitude, above the fainter side of its neighbors.
    """
    y = np.asarray(y_model, dtype=float)
    span = float(np.ptp(y))
    if not np.isfinite(span) or span < 0.01:
        return 0
    wrapped = np.concatenate([y, y, y])
    prominence = max(0.02, 0.15 * span)
    idx, _ = find_peaks(-wrapped, prominence=prominence)
    n = len(y)
    return int(np.sum((idx >= n) & (idx < 2 * n)))


def _is_baseline_alias(labels):
    """True when a candidate is a daily alias of the campaign duration."""
    return any(str(label).startswith(("span", "2-span")) for label in labels)


def _scatter_about_curve(y_points, y_curve, sigma, trim=0.10):
    """Mean distance of the points to the model curve, errors included.

    The farthest `trim` fraction of the distances is left out. A point that
    lies on the curve still widens the mean when its magnitude error is large:
    each kept point contributes sqrt(distance^2 + error^2).
    """
    distance = np.abs(np.asarray(y_points, dtype=float) - np.asarray(y_curve, dtype=float))
    sigma = np.maximum(np.asarray(sigma, dtype=float), 0.0)
    n = len(distance)
    if n == 0:
        return np.nan
    n_drop = int(np.floor(trim * n))
    n_keep = max(n - n_drop, 1)
    keep = np.argsort(distance)[:n_keep]
    width = np.sqrt(distance[keep] ** 2 + sigma[keep] ** 2)
    return float(np.mean(width))


def _data_amplitude(y, low=2.5, high=97.5):
    """Peak-to-peak of the magnitudes between two percentiles.

    The points are ordered by magnitude. The amplitude is the 97.5th
    percentile minus the 2.5th, so the outer 5% of the points do not set
    the range. The second and third values are those two magnitudes.
    """
    y = np.asarray(y, dtype=float)
    y = y[np.isfinite(y)]
    if len(y) < 2:
        return np.nan, np.nan, np.nan
    faint, bright = np.percentile(y, [high, low])
    return float(faint - bright), float(bright), float(faint)


def _red(text):
    """Paint a line red on a terminal. FORCE_COLOR keeps the color in a capture."""
    forced = os.environ.get("FORCE_COLOR", "") not in ("", "0")
    try:
        tty = bool(sys.stdout.isatty())
    except Exception:
        tty = False
    if not forced and not tty:
        return text
    return f"\033[31m{text}\033[0m"


def _campaign_alias_k_range(span_days, p_min_h=ALIAS_P_MIN_H):
    """Integer shifts from the long edge of the search down to p_min_h.

    The rotation period is 48/f hours. The long edge is 2*PERIODS_H[1] hours,
    so its frequency is 1 cycle/day, and a period of p_min_h hours has
    frequency 48/p_min_h. The aliases are 1/T + k and 2/T + k. k_min is set
    from 2/T, which needs the smaller shift to reach the long edge, and k_max
    from 1/T, which needs the larger shift to reach p_min_h. Both families
    then use every integer from k_min to k_max.
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


def _baseline_alias_items(t_days, p_min_h=ALIAS_P_MIN_H):
    """Daily aliases of the campaign duration, in cycles per day.

    The observations span T days. One cycle over that span has frequency 1/T,
    and one rotation over that span has frequency 2/T. The aliases are those
    frequencies shifted by an integer number of cycles per day, from the long
    edge of the period search down to p_min_h. Folding at one of them piles
    the nights into a few phase clusters. Used to mark frequencies and to draw
    the alias lines. The two-peak comparison does not fit these as candidates.
    """
    t = np.asarray(t_days, dtype=float)
    span = float(np.ptp(t))
    if not np.isfinite(span) or span <= 0:
        return []
    f_min = 24.0 / PERIODS_H[1]
    f_max = 24.0 / PERIODS_H[0]
    p_min = cadence_min_period(t)
    items = []
    for base, tag in ((1.0 / span, "span"), (2.0 / span, "2-span")):
        for k in _campaign_alias_k_range(span, p_min_h):
            freq = base + float(k)
            if freq < f_min or freq > f_max:
                continue
            if p_min > 0 and (24.0 / freq) < p_min:
                continue
            items.append((freq, f"{tag} k={k:+d}"))
    return items


def _tag_campaign_aliases(merged, t_days, tol=0.1, p_min_h=ALIAS_P_MIN_H):
    """Append span / 2-span labels to voted frequencies that land on an alias.

    The chi-squared comparison fits the voted frequencies only. An alias that
    no series voted for is left unmarked here and is not scored.
    """
    alias_items = _baseline_alias_items(t_days, p_min_h)
    if not alias_items:
        return merged
    for group in merged:
        freq = float(group["freq"])
        for alias_freq, label in alias_items:
            if abs(freq - float(alias_freq)) > tol:
                continue
            if label not in group["labels"]:
                group["labels"].append(label)
    return merged


# Bluer filter first. Every pair among the filters that have a zero point is a color.
COLOR_FILTERS = ("u", "g", "r", "i", "z", "y")


def _color_pairs(filters=COLOR_FILTERS):
    """Every color a-b with a bluer than b, in wavelength order."""
    return [
        (filters[i], filters[j])
        for i in range(len(filters))
        for j in range(i + 1, len(filters))
    ]


COMPARE_COLOR_COLUMNS = tuple(
    name
    for a, b in _color_pairs()
    for name in (f"{a}-{b}", f"{a}-{b}_err")
)


def _color_difference(a, b, color, sig_color, cov, names, ref):
    """Magnitude a minus magnitude b, and its uncertainty, from the joint fit."""
    value = float(color[a] - color[b])
    name_a = None if a == ref else f"c_{a}"
    name_b = None if b == ref else f"c_{b}"
    if name_a in names and name_b in names:
        i = names.index(name_a)
        j = names.index(name_b)
        var = cov[i, i] + cov[j, j] - 2.0 * cov[i, j]
        return value, float(np.sqrt(max(var, 0.0)))
    if name_a in names:
        return value, float(sig_color[a])
    if name_b in names:
        return value, float(sig_color[b])
    err = float(np.sqrt(sig_color.get(a, 0.0) ** 2 + sig_color.get(b, 0.0) ** 2))
    return value, err


def _u_minus_reference(df_obj, beta, f_rot, K, present, ref):
    """u minus the reference curve. u stays out of the joint harmonic fit."""
    alias = {"Lg": "g", "Lr": "r", "Li": "i", "Lz": "z", "Ly": "y", "Lu": "u"}
    df_u = df_obj.dropna(subset=[PROP, "band"]).copy()
    df_u["filt"] = df_u["band"].map(lambda b: alias.get(b, b))
    df_u = df_u[df_u["filt"] == "u"]
    if df_u.empty or "t_days" not in df_u.columns:
        return None
    rms = pd.to_numeric(df_u["rmsmag"], errors="coerce").to_numpy(float)
    t_u = df_u["t_days"].to_numpy(float)
    y_u = df_u[PROP].to_numpy(float)
    s_u = np.where(np.isfinite(rms) & (rms > 0), rms, 0.01)
    X_ref, _ = _design_matrix(
        t_u, np.full(len(t_u), ref), f_rot, K, present, ref
    )
    d_u = y_u - X_ref @ beta
    good = np.ones(len(d_u), bool)
    c_u = 0.0
    for _ in range(30):
        if good.sum() < 2:
            break
        scale = 1.4826 * np.median(np.abs(d_u[good] - np.median(d_u[good])))
        w = 1.0 / np.maximum(s_u, 1e-4) ** 2
        c_u = float(np.sum(w[good] * d_u[good]) / np.sum(w[good]))
        new = np.abs(d_u - c_u) < 5.0 * max(scale, np.median(s_u))
        if (new == good).all():
            good = new
            break
        good = new
    if int(good.sum()) < 2:
        return None
    w = 1.0 / np.maximum(s_u, 1e-4) ** 2
    c_u = float(np.sum(w[good] * d_u[good]) / np.sum(w[good]))
    chi2_red = float(np.sum(((d_u[good] - c_u) / s_u[good]) ** 2) / max(int(good.sum()) - 1, 1))
    err_u = float(np.sqrt(1.0 / np.sum(w[good]) * max(chi2_red, 1.0)))
    return c_u, err_u


def _compare_colors(df_obj, best, present, ref):
    """Color indices of the adopted two-peak model."""
    beta = best["beta"]
    names = best["names"]
    cov = best["cov"]
    color = {ref: 0.0}
    sig_color = {ref: 0.0}
    beta_map = dict(zip(names, beta))
    err_diag = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    err_map = dict(zip(names, err_diag))
    for filt in present:
        if filt == ref:
            continue
        color[filt] = float(beta_map[f"c_{filt}"])
        sig_color[filt] = float(err_map[f"c_{filt}"])
    joint = [filt for filt in COLOR_FILTERS if filt != "u" and filt in color]
    colors = {}
    for i, a in enumerate(joint):
        for b in joint[i + 1 :]:
            value, err = _color_difference(a, b, color, sig_color, cov, names, ref)
            colors[f"{a}-{b}"] = value
            colors[f"{a}-{b}_err"] = err
    # u is not in the joint fit. Its zero point is measured against the reference
    # curve, so u minus any fitted filter is that offset minus the filter's zero point.
    u_ref = _u_minus_reference(df_obj, beta, best["f_rot"], 2, present, ref)
    if u_ref is not None:
        c_u, err_u = u_ref
        for filt in joint:
            colors[f"u-{filt}"] = c_u - color[filt]
            colors[f"u-{filt}_err"] = float(
                np.sqrt(err_u**2 + sig_color.get(filt, 0.0) ** 2)
            )
    return colors


def _band_producing_period(df_obj, labels):
    """Single band with the most points among those that voted for a period."""
    single = {"Lg", "Lr", "Li", "Lz", "Ly", "Lu"}
    bands = []
    for label in labels:
        text = str(label)
        if text.startswith(("span", "2-span")):
            continue
        if text in single:
            bands.append(text)
            continue
        if "+" in text:
            for part in text.split("+"):
                if part in single:
                    bands.append(part)
    if not bands:
        return None
    counts = {
        band: int(((df_obj["band"] == band) & df_obj[PROP].notna()).sum())
        for band in set(bands)
    }
    return max(counts, key=counts.get)


def _leave_one_night_out_peak(df_obj, band, freq_ref, tol=0.1):
    """Leave-one-night-out Lomb–Scargle peak for one band.

    For each calendar night left out, the periodogram peak of the remaining
    points is compared with freq_ref. The period is stable when at least half
    of the folds with at least five points keep a peak within tol cycles/day.
    At least three nights and two successful folds are required.
    """
    empty = {
        "stable": False,
        "n_nights": 0,
        "n_tested": 0,
        "n_kept": 0,
        "periods_h": [],
        "folds": [],
        "band": band,
        "reason": "no band",
    }
    if band is None:
        return empty
    df_band = df_obj[(df_obj["band"] == band) & df_obj[PROP].notna()].copy()
    if df_band.empty:
        empty["reason"] = f"{band} has no points"
        return empty
    if "night" not in df_band.columns:
        df_band["night"] = pd.to_datetime(df_band["obstime"]).dt.floor("D")
    nights = sorted(df_band["night"].dropna().unique())
    n_nights = len(nights)
    if n_nights < 3:
        return {
            "stable": False,
            "n_nights": n_nights,
            "n_tested": 0,
            "n_kept": 0,
            "periods_h": [],
            "folds": [],
            "band": band,
            "reason": f"fewer than 3 nights ({n_nights})",
        }
    kept = 0
    tested = 0
    periods_h = []
    folds = []
    for night in nights:
        part = df_band[df_band["night"] != night]
        if len(part) < 5:
            continue
        t = part["t_days"].to_numpy(float)
        y = part[PROP].to_numpy(float)
        dy = pd.to_numeric(part["rmsmag"], errors="coerce").fillna(0.01).to_numpy(float)
        freq_grid = frequency_grid(t)
        power = np.asarray(LombScargle(t, y, dy).power(freq_grid), dtype=float)
        p_phot_min = cadence_min_period(t)
        allowed = (24.0 / freq_grid) >= p_phot_min
        if not np.any(allowed):
            allowed = np.ones(len(freq_grid), dtype=bool)
        index = int(np.flatnonzero(allowed)[np.argmax(power[allowed])])
        f_peak = float(freq_grid[index])
        p_rot = 48.0 / f_peak
        power_peak = float(power[index])
        periods_h.append(p_rot)
        tested += 1
        fold_kept = abs(f_peak - float(freq_ref)) <= tol
        if fold_kept:
            kept += 1
        night_str = pd.to_datetime(night).strftime("%Y-%m-%d")
        folds.append(
            {
                "night": night_str,
                "p_rot_h": float(p_rot),
                "freq": f_peak,
                "power": power_peak,
                "kept": bool(fold_kept),
                "n_points": int(len(part)),
            }
        )
    # At least half the folds must keep the peak (a true rotation survives
    # dropping individual nights more often than a campaign-window artifact).
    stable = tested >= 2 and kept >= 0.5 * tested
    reason = (
        f"kept {kept}/{tested} folds within {tol:.1f} d^-1 of "
        f"{48.0 / float(freq_ref):.2f} h "
        f"({'stable' if stable else 'unstable'}, need >= half)"
    )
    return {
        "stable": bool(stable),
        "n_nights": n_nights,
        "n_tested": tested,
        "n_kept": kept,
        "periods_h": periods_h,
        "folds": folds,
        "band": band,
        "reason": reason,
    }


def multiband_photometry_compare(df_obj, provid, candidates, nightalias=False):
    """Fit a two-peak multiband model at each voted period and keep the best.

    Candidates are the periods from the democratic vote. Each is compared with
    the campaign aliases 1/T+k and 2/T+k; a match within 0.1 cycles/day gets
    the alias label and is not adopted. Aliases that no series voted for are
    drawn on the periodogram and are not fitted. The model is the joint
    harmonic fit of order 2. Filters share the curve and keep their own zero
    point. u is left out. The figure shows the adopted period and the other
    voted periods that have two peaks. If no model has two peaks, or every
    two-peak model is a campaign alias, the asteroid is left out.

    When nightalias is True and the lowest-chi2 two-peak model is a campaign
    alias, a leave-one-night-out periodogram of the best-sampled band that
    voted for that period decides whether the alias is adopted. The peak is
    stable when at least half of the folds keep it; that period is then kept
    with its alias label. An unstable peak falls back to the best non-alias
    two-peak model, or leaves the asteroid out.
    """
    tag = "nightalias" if nightalias else "compare"
    figure_kind = (
        "multiband_photometry_nightalias"
        if nightalias
        else "multiband_photometry_compare"
    )
    alias = {"Lg": "g", "Lr": "r", "Li": "i", "Lz": "z", "Ly": "y", "Lu": "u"}
    filt_order = ["g", "r", "i", "z", "y"]
    ref = "r"
    df_fit = df_obj.dropna(subset=[PROP, "band"]).copy()
    df_fit["filt"] = df_fit["band"].map(lambda b: alias.get(b, b))
    rms = pd.to_numeric(df_fit["rmsmag"], errors="coerce").to_numpy(float)
    df_fit["rmsmag"] = np.where(np.isfinite(rms) & (rms > 0), rms, 0.01)
    df_fit = df_fit[df_fit["filt"].isin(filt_order)]
    kept = []
    n_rejected = 0
    for filt_name, part in df_fit.groupby("filt", sort=False):
        good = _magnitude_inliers(part[PROP].to_numpy())
        n_rejected += int((~good).sum())
        kept.append(part.loc[good])
    if kept:
        df_fit = pd.concat(kept, ignore_index=True)
    if n_rejected:
        print(
            f"  multiband 2-peak {tag}: rejected {n_rejected} magnitude outlier(s)"
        )
    present = [f for f in filt_order if (df_fit["filt"] == f).any()]
    if len(present) < 2:
        print(
            f"  skip multiband_photometry_{tag} for {provid}: "
            f"only {len(present)} filter(s) in {present}"
        )
        return None
    if ref not in present:
        counts = {f: int((df_fit["filt"] == f).sum()) for f in present}
        ref = max(counts, key=counts.get)

    t_fit = df_fit["t_days"].to_numpy(float)
    y_fit = df_fit[PROP].to_numpy(float)
    s_fit = np.maximum(df_fit["rmsmag"].to_numpy(float), 1e-4)
    f_fit = df_fit["filt"].to_numpy()
    mask = np.ones(len(y_fit), bool)
    t0 = float(np.min(t_fit))
    K = 2
    vote_items = [
        (float(cand["freq"]), label)
        for cand in candidates
        for label in cand["labels"]
    ]
    merged = _tag_campaign_aliases(_merge_frequencies(vote_items), t_fit)
    n_alias = sum(1 for cand in merged if _is_baseline_alias(cand["labels"]))
    print(
        f"  multiband 2-peak {tag}: {len(merged)} voted periods, "
        f"{n_alias} land on a campaign alias"
    )
    ls_power = LombScargleMultiband(t_fit, y_fit, f_fit, s_fit)
    scored = []
    for cand in merged:
        freq = float(cand["freq"])
        f_rot = freq / 2.0
        chi2, kpar = _chi2_fixed(
            t_fit, y_fit, s_fit, f_fit, present, ref, f_rot, K, mask, 0.0
        )
        X, names = _design_matrix(t_fit, f_fit, f_rot, K, present, ref)
        beta, _ = _wls_fit(X, y_fit, s_fit, mask)
        phase_grid = np.linspace(0, 1, 400, endpoint=False)
        Xg, _ = _design_matrix(
            t0 + phase_grid / f_rot,
            np.full(len(phase_grid), ref),
            f_rot,
            K,
            present,
            ref,
        )
        y_model = Xg @ beta
        offset = {ref: 0.0}
        for name, value in zip(names, beta):
            if name.startswith("c_"):
                offset[name[2:]] = float(value)
        y_on_ref = y_fit - np.array([offset.get(f, 0.0) for f in f_fit])
        y_curve = _design_matrix(
            t_fit, np.full(len(t_fit), ref), f_rot, K, present, ref
        )[0] @ beta
        amp = float(np.ptp(y_model))
        scatter = _scatter_about_curve(y_on_ref, y_curve, s_fit)
        ratio = amp / scatter if scatter > 0 else np.nan
        power = float(np.asarray(ls_power.power([freq]), dtype=float).ravel()[0])
        scored.append(
            {
                "freq": freq,
                "f_rot": f_rot,
                "p_rot_h": 48.0 / freq,
                "labels": list(cand["labels"]),
                "chi2": float(chi2),
                "chi2_nu": float(chi2) / max(int(mask.sum()) - kpar, 1),
                "n_peaks": _brightness_peaks(y_model),
                "power": power,
                "amp": amp,
                "scatter": scatter,
                "amp_ratio": ratio,
                "span_alias": int(_is_baseline_alias(cand["labels"])),
                "offset": offset,
                "beta": beta,
                "names": names,
                "cov": np.linalg.pinv(
                    (X[mask] / s_fit[mask, None]).T @ (X[mask] / s_fit[mask, None])
                ),
                "y_model": y_model,
                "y_on_ref": y_on_ref,
                "phase_grid": phase_grid,
            }
        )
    scored.sort(key=lambda item: (item["chi2_nu"], item["p_rot_h"]))
    two_peaks = [item for item in scored if item["n_peaks"] == 2]
    # A campaign alias can fit the nights, so it is not a rotation period.
    rotations = [item for item in two_peaks if not item["span_alias"]]
    nightalias_tested = 0
    nightalias_stable = 0
    nightalias_band = ""
    nightalias_n_nights = 0
    nightalias_n_kept = 0
    nightalias_n_tested = 0
    nightalias_periods = ""
    best = None
    if nightalias and two_peaks and two_peaks[0]["span_alias"]:
        candidate = two_peaks[0]
        band_lono = _band_producing_period(df_obj, candidate["labels"])
        lono = _leave_one_night_out_peak(df_obj, band_lono, candidate["freq"])
        nightalias_tested = 1
        nightalias_stable = int(lono["stable"])
        nightalias_band = str(lono.get("band") or "")
        nightalias_n_nights = int(lono["n_nights"])
        nightalias_n_kept = int(lono["n_kept"])
        nightalias_n_tested = int(lono["n_tested"])
        nightalias_periods = ",".join(
            f"{fold['p_rot_h']:.4f}:{fold['power']:.3f}"
            for fold in lono.get("folds", [])
        )
        print(
            f"  multiband 2-peak {tag}: leave-one-night-out on {nightalias_band} "
            f"at {candidate['p_rot_h']:.4f} h: {lono['reason']}"
        )
        for fold in lono.get("folds", []):
            mark = " kept" if fold["kept"] else ""
            print(
                f"    leave out {fold['night']}: P_rot = {fold['p_rot_h']:.4f} h, "
                f"power = {fold['power']:.3f} (N = {fold['n_points']}){mark}"
            )
        if lono["stable"]:
            best = candidate
            print(
                f"  multiband 2-peak {tag}: campaign alias kept after "
                f"leave-one-night-out ({candidate['p_rot_h']:.4f} h)"
            )
        else:
            best = rotations[0] if rotations else None
            if best is None:
                print(
                    f"  multiband 2-peak {tag}: campaign alias rejected by "
                    f"leave-one-night-out; {provid} is left out"
                )
            else:
                print(
                    f"  multiband 2-peak {tag}: campaign alias rejected by "
                    f"leave-one-night-out; falling back to {best['p_rot_h']:.4f} h"
                )
    else:
        best = rotations[0] if rotations else None
    periods = ",".join(
        f"{item['p_rot_h']:.4f}:{item['chi2_nu']:.2f}:{item['n_peaks']}:"
        f"{item['power']:.3f}:{item['amp_ratio']:.3f}:{item['span_alias']}"
        for item in scored
    )
    for item in scored:
        mark = " selected" if item is best else ""
        alias_mark = ", baseline alias" if item["span_alias"] else ""
        line = (
            f"  {item['p_rot_h']:.4f} h ({', '.join(item['labels'])}): "
            f"chi2_nu = {item['chi2_nu']:.2f}, peaks = {item['n_peaks']}, "
            f"power = {item['power']:.3f}, A/s = {item['amp_ratio']:.3f}"
            f"{alias_mark}{mark}"
        )
        print(_red(line) if item["span_alias"] else line)
    path = output_path(figure_kind, provid)
    if best is None:
        if path.exists():
            path.unlink()
        if two_peaks and nightalias_tested and not nightalias_stable:
            print(
                f"  multiband 2-peak {tag}: leave-one-night-out rejected "
                f"the campaign alias; {provid} is left out"
            )
            reason = "Campaign alias failed leave-one-night-out"
        elif two_peaks:
            print(
                f"  multiband 2-peak {tag}: the two-peak models are "
                f"campaign aliases; {provid} is left out"
            )
            reason = "Only campaign aliases have two peaks"
        else:
            print(
                f"  multiband 2-peak {tag}: no period with 2 peaks; "
                f"{provid} is left out"
            )
            reason = "No two-peak period in the multiband comparison"
        return {
            "no_two_peaks": True,
            "reject_reason": reason,
            "path": None,
            "periods": periods,
            "nightalias_tested": nightalias_tested,
            "nightalias_stable": nightalias_stable,
            "nightalias_band": nightalias_band,
            "nightalias_n_nights": nightalias_n_nights,
            "nightalias_n_tested": nightalias_n_tested,
            "nightalias_n_kept": nightalias_n_kept,
            "nightalias_periods": nightalias_periods,
        }
    n_rot_label = len(rotations) if not best["span_alias"] else len(two_peaks)
    summary = (
        f"  multiband 2-peak {tag}: P_rot = {best['p_rot_h']:.4f} h, "
        f"chi2_nu = {best['chi2_nu']:.2f}, power = {best['power']:.3f}, "
        f"A/s = {best['amp_ratio']:.3f}, peaks = 2 "
        f"(best of {n_rot_label} rotations, {len(scored)} scored)"
    )
    if nightalias_tested:
        summary += (
            f", leave-one-night-out stable = {nightalias_stable} "
            f"({nightalias_band})"
        )
    print(summary)

    def voted(item):
        return any(
            not str(label).startswith(("span", "2-span")) for label in item["labels"]
        )

    others = [item for item in two_peaks if item is not best and voted(item)]
    n_other = len(others)
    nrow_other = int(np.ceil(n_other / 2.0)) if n_other else 0
    fig = plt.figure(
        figsize=(10.5, 5.8 + 2.9 * nrow_other),
        layout="constrained",
    )

    def draw_panel(ax, item, selected):
        phase = np.mod((t_fit - t0) * item["f_rot"], 1.0)
        y_on_ref = item["y_on_ref"]
        beta_of = dict(zip(item["names"], item["beta"]))
        err_of = dict(zip(
            item["names"],
            np.sqrt(np.clip(np.diag(item["cov"]), 0.0, None)),
        ))
        for filt in present:
            sel = f_fit == filt
            if filt == ref:
                filt_label = rf"{filt}, $c_{{{filt}}} = 0$"
            else:
                filt_label = (
                    rf"{filt}, $c_{{{filt}}} = {float(beta_of[f'c_{filt}']):+.3f}"
                    rf" \pm {float(err_of[f'c_{filt}']):.3f}$"
                )
            ax.errorbar(
                phase[sel],
                y_on_ref[sel],
                yerr=s_fit[sel],
                fmt="o",
                ms=5.5 if selected else 2.5,
                alpha=0.85 if selected else 0.7,
                color=band_colors.get(f"L{filt}", "k"),
                ecolor="0.75",
                label=filt_label,
            )
        ax.plot(
            item["phase_grid"],
            item["y_model"],
            color="0.1",
            lw=2.2 if selected else 1.3,
            label="2-peak harmonic model",
        )
        if selected:
            _, bright, faint = _data_amplitude(y_on_ref)
            if np.isfinite(bright) and np.isfinite(faint):
                ax.axhline(
                    bright,
                    color="0.25",
                    ls=":",
                    lw=1.2,
                    label="data 2.5–97.5%",
                )
                ax.axhline(faint, color="0.25", ls=":", lw=1.2)
        if selected and item["span_alias"] and nightalias_stable:
            mark = ", selected, leave-one-night-out"
            title_color = "blue"
        elif selected:
            mark = ", selected"
            title_color = "black"
        else:
            mark = ""
            title_color = "red" if item["span_alias"] else "black"
        ax.set_title(
            rf"$P_{{rot}}$ = {item['p_rot_h']:.2f} h, "
            rf"$\chi^2_\nu$ = {item['chi2_nu']:.1f}, "
            f"{item['n_peaks']} peaks{mark}\n"
            rf"power = {item['power']:.3f}, "
            rf"$A$ = {item['amp']:.3f}, $A/s$ = {item['amp_ratio']:.2f}"
            "\n" + ", ".join(item["labels"]),
            fontsize=11 if selected else 8,
            color=title_color,
        )
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        if selected:
            ax.tick_params(labelsize=10)

    if n_other == 0:
        ax_best = fig.add_subplot(1, 1, 1)
        other_axes = []
    else:
        gs = fig.add_gridspec(
            1 + nrow_other,
            2,
            height_ratios=[2.15] + [1.0] * nrow_other,
        )
        ax_best = fig.add_subplot(gs[0, :])
        other_axes = [
            fig.add_subplot(
                gs[1 + i // 2, i % 2],
                sharex=ax_best,
                sharey=ax_best,
            )
            for i in range(n_other)
        ]
    draw_panel(ax_best, best, True)
    for ax, item in zip(other_axes, others):
        draw_panel(ax, item, False)
    photometry = np.concatenate(
        [np.asarray(item["y_on_ref"], dtype=float) for item in [best, *others]]
    )
    photometry = photometry[np.isfinite(photometry)]
    bright = float(np.min(photometry))
    faint = float(np.max(photometry))
    span = max(faint - bright, 0.05)
    pad = max(0.05, 0.15 * span)
    ax_best.set_ylim(faint + pad, bright - pad)
    ax_best.legend(loc="best", frameon=True, fontsize=7)
    ax_best.set_ylabel("Absolute magnitude, H")
    if n_other == 0:
        ax_best.set_xlabel("Phase at the rotation period ($2P$)")
    else:
        bottom = {}
        for i, ax in enumerate(other_axes):
            bottom[i % 2] = ax
            if i % 2 == 0:
                ax.set_ylabel("Absolute magnitude, H")
            else:
                ax.tick_params(labelleft=False)
        for ax in bottom.values():
            ax.set_xlabel("Phase at the rotation period ($2P$)")
    title = "nightalias" if nightalias else "comparison"
    fig.suptitle(
        f"Multiband 2-peak {title} of {provid} ({PROP})",
        fontsize=12,
    )
    save_figure(fig, path)
    colors = _compare_colors(df_obj, best, present, ref)
    if colors:
        shown = ", ".join(
            f"{key} = {colors[key]:+.3f} ± {colors[key + '_err']:.3f}"
            for key in (f"{a}-{b}" for a, b in _color_pairs())
            if key in colors
        )
        print(f"  colors of the adopted period: {shown}")
    return {
        "p_rot_h": float(best["p_rot_h"]),
        "chi2_nu": float(best["chi2_nu"]),
        "n_peaks": int(best["n_peaks"]),
        "power": float(best["power"]),
        "amp_ratio": float(best["amp_ratio"]),
        "amp_data": _data_amplitude(best["y_on_ref"])[0],
        "span_alias": int(best["span_alias"]),
        "periods": periods,
        "colors": colors,
        "path": path,
        "no_two_peaks": False,
        "nightalias_tested": nightalias_tested,
        "nightalias_stable": nightalias_stable,
        "nightalias_band": nightalias_band,
        "nightalias_n_nights": nightalias_n_nights,
        "nightalias_n_tested": nightalias_n_tested,
        "nightalias_n_kept": nightalias_n_kept,
        "nightalias_periods": nightalias_periods,
    }


def multiband_photometry_nightalias(df_obj, provid, candidates):
    """Two-peak compare with leave-one-night-out for campaign aliases."""
    return multiband_photometry_compare(
        df_obj, provid, candidates, nightalias=True
    )


def _given_periods(df_rot, provid):
    """Rotation periods stored for this asteroid.

    Lomb-Scargle and Fourier are both kept when they differ. A single panel
    is enough when the two values agree.
    """
    row = df_rot.loc[df_rot["provid"] == provid]
    if row.empty:
        return []
    found = []
    for column, label, key in (
        ("period_lsm_hr", "Literature, Lomb-Scargle", "lsm"),
        ("period_fourier_hr", "Literature, Fourier", "fourier"),
    ):
        if column not in row.columns:
            continue
        value = pd.to_numeric(row.iloc[0][column], errors="coerce")
        if np.isfinite(value) and float(value) > 0:
            found.append((float(value), label, key))
    if len(found) == 2:
        first, second = found[0][0], found[1][0]
        if abs(first - second) <= max(0.01, 0.005 * first):
            return [(first, "Literature, Lomb-Scargle and Fourier", "lsm+fourier")]
    return found


def justfold(df_obj, provid, df_rot):
    """Fold the catalogue period with the same two-peak multiband model.

    No search and no other candidates. The period is the one stored for this
    asteroid. Filters share one harmonic curve of order 2 and keep a zero
    point each. u is left out.
    """
    given = _given_periods(df_rot, provid)
    if not given:
        print(f"  skip justfold for {provid}: no literature rotation period")
        return None
    alias = {"Lg": "g", "Lr": "r", "Li": "i", "Lz": "z", "Ly": "y", "Lu": "u"}
    filt_order = ["g", "r", "i", "z", "y"]
    ref = "r"
    df_fit = df_obj.dropna(subset=[PROP, "band"]).copy()
    df_fit["filt"] = df_fit["band"].map(lambda b: alias.get(b, b))
    rms = pd.to_numeric(df_fit["rmsmag"], errors="coerce").to_numpy(float)
    df_fit["rmsmag"] = np.where(np.isfinite(rms) & (rms > 0), rms, 0.01)
    df_fit = df_fit[df_fit["filt"].isin(filt_order)]
    kept = []
    n_rejected = 0
    for _, part in df_fit.groupby("filt", sort=False):
        good = _magnitude_inliers(part[PROP].to_numpy())
        n_rejected += int((~good).sum())
        kept.append(part.loc[good])
    if kept:
        df_fit = pd.concat(kept, ignore_index=True)
    if n_rejected:
        print(f"  justfold: rejected {n_rejected} magnitude outlier(s)")
    present = [f for f in filt_order if (df_fit["filt"] == f).any()]
    if len(present) < 2:
        print(
            f"  skip justfold for {provid}: only {len(present)} filter(s) in {present}"
        )
        return None
    if ref not in present:
        counts = {f: int((df_fit["filt"] == f).sum()) for f in present}
        ref = max(counts, key=counts.get)

    t_fit = df_fit["t_days"].to_numpy(float)
    y_fit = df_fit[PROP].to_numpy(float)
    s_fit = np.maximum(df_fit["rmsmag"].to_numpy(float), 1e-4)
    f_fit = df_fit["filt"].to_numpy()
    mask = np.ones(len(y_fit), bool)
    t0 = float(np.min(t_fit))
    K = 2
    ls_power = LombScargleMultiband(t_fit, y_fit, f_fit, s_fit)
    alias_freqs = [freq for freq, _ in _baseline_alias_items(t_fit)]
    scored = []
    for p_rot_h, label, key in given:
        freq = 48.0 / p_rot_h
        f_rot = freq / 2.0
        chi2, kpar = _chi2_fixed(
            t_fit, y_fit, s_fit, f_fit, present, ref, f_rot, K, mask, 0.0
        )
        X, names = _design_matrix(t_fit, f_fit, f_rot, K, present, ref)
        beta, _ = _wls_fit(X, y_fit, s_fit, mask)
        phase_grid = np.linspace(0, 1, 400, endpoint=False)
        Xg, _ = _design_matrix(
            t0 + phase_grid / f_rot,
            np.full(len(phase_grid), ref),
            f_rot,
            K,
            present,
            ref,
        )
        y_model = Xg @ beta
        offset = {ref: 0.0}
        for name, value in zip(names, beta):
            if name.startswith("c_"):
                offset[name[2:]] = float(value)
        y_on_ref = y_fit - np.array([offset.get(f, 0.0) for f in f_fit])
        y_curve = _design_matrix(
            t_fit, np.full(len(t_fit), ref), f_rot, K, present, ref
        )[0] @ beta
        amp = float(np.ptp(y_model))
        scatter = _scatter_about_curve(y_on_ref, y_curve, s_fit)
        ratio = amp / scatter if scatter > 0 else np.nan
        power = float(np.asarray(ls_power.power([freq]), dtype=float).ravel()[0])
        span_alias = int(any(abs(freq - alias_f) <= 0.1 for alias_f in alias_freqs))
        item = {
            "key": key,
            "label": label,
            "p_rot_h": float(p_rot_h),
            "f_rot": f_rot,
            "chi2_nu": float(chi2) / max(int(mask.sum()) - kpar, 1),
            "n_peaks": _brightness_peaks(y_model),
            "power": power,
            "amp_ratio": ratio,
            "span_alias": span_alias,
            "y_model": y_model,
            "y_on_ref": y_on_ref,
            "phase_grid": phase_grid,
        }
        scored.append(item)
        line = (
            f"  justfold {label}: P_rot = {item['p_rot_h']:.4f} h, "
            f"chi2_nu = {item['chi2_nu']:.2f}, peaks = {item['n_peaks']}, "
            f"power = {item['power']:.3f}, A/s = {item['amp_ratio']:.3f}"
        )
        if span_alias:
            line += ", baseline alias"
        print(_red(line) if span_alias else line)

    path = output_path("justfold", provid)
    n = len(scored)
    fig = plt.figure(figsize=(10.5, 5.6 * n), layout="constrained")
    axes = []
    for i in range(n):
        if i == 0:
            ax = fig.add_subplot(n, 1, i + 1)
        else:
            ax = fig.add_subplot(n, 1, i + 1, sharex=axes[0], sharey=axes[0])
        axes.append(ax)
    for ax, item in zip(axes, scored):
        phase = np.mod((t_fit - t0) * item["f_rot"], 1.0)
        for filt in present:
            sel = f_fit == filt
            ax.errorbar(
                phase[sel],
                item["y_on_ref"][sel],
                yerr=s_fit[sel],
                fmt="o",
                ms=5.5,
                alpha=0.85,
                color=band_colors.get(f"L{filt}", "k"),
                ecolor="0.75",
                label=filt,
            )
        ax.plot(item["phase_grid"], item["y_model"], color="0.1", lw=2.2, label="2-peak model")
        ax.set_title(
            rf"$P_{{rot}}$ = {item['p_rot_h']:.2f} h, "
            rf"$\chi^2_\nu$ = {item['chi2_nu']:.1f}, "
            f"{item['n_peaks']} peaks\n"
            rf"power = {item['power']:.3f}, $A/s$ = {item['amp_ratio']:.2f}"
            f"\n{item['label']}",
            fontsize=11,
            color="red" if item["span_alias"] else "black",
        )
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 1)
        ax.tick_params(labelsize=10)
        ax.set_ylabel("Absolute magnitude, H")
    axes[0].legend(loc="best", frameon=True, fontsize=8)
    axes[-1].set_xlabel("Phase at the rotation period ($2P$)")
    photometry = np.concatenate(
        [np.asarray(item["y_on_ref"], dtype=float) for item in scored]
    )
    photometry = photometry[np.isfinite(photometry)]
    bright = float(np.min(photometry))
    faint = float(np.max(photometry))
    span = max(faint - bright, 0.05)
    pad = max(0.05, 0.15 * span)
    axes[0].set_ylim(faint + pad, bright - pad)
    fig.suptitle(f"Fold at the given period of {provid} ({PROP})", fontsize=12)
    save_figure(fig, path)
    first = scored[0]
    periods = ",".join(
        f"{item['key']}:{item['p_rot_h']:.4f}:{item['chi2_nu']:.2f}:"
        f"{item['n_peaks']}:{item['power']:.3f}:{item['amp_ratio']:.3f}:"
        f"{item['span_alias']}"
        for item in scored
    )
    return {
        "p_rot_h": float(first["p_rot_h"]),
        "chi2_nu": float(first["chi2_nu"]),
        "n_peaks": int(first["n_peaks"]),
        "power": float(first["power"]),
        "amp_ratio": float(first["amp_ratio"]),
        "span_alias": int(first["span_alias"]),
        "periods": periods,
        "path": path,
    }


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
        out_k["mag"] = np.nan
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
        path_obs_csv = FIG_DIR / f"{slug(provid)}_observations.csv"
        observed = real.copy()
        observed["provid"] = provid
        observed["obs_filt"] = observed["filt"]
        observed["source_filt"] = observed["filt"]
        df_lk["provid"] = provid
        df_export = pd.concat([observed, df_lk], ignore_index=True)
        df_export["mag"] = pd.to_numeric(df_export["mag"], errors="coerce")
        df_export[PROP] = pd.to_numeric(df_export[PROP], errors="coerce")
        is_lk = df_export["band"].eq("Lk")
        df_export.loc[is_lk, "mag"] = np.nan
        valid_h = df_export[PROP].notna() & (df_export[PROP] != 0)
        valid_mag = df_export["mag"].notna() & (df_export["mag"] != 0)
        df_export = df_export.loc[(is_lk & valid_h) | (~is_lk & valid_mag & valid_h)]
        df_export = df_export.sort_values(
            ["obstime", "band"], kind="mergesort"
        ).reset_index(drop=True)
        cols_export = [
            "provid",
            "obstime",
            "t_days",
            "night",
            "band",
            "obs_filt",
            "source_filt",
            "mag",
            PROP,
            "rmsmag",
        ]
        df_export[cols_export].to_csv(path_obs_csv, index=False, na_rep="NaN")
        stale_k_csv = FIG_DIR / f"{slug(provid)}_band_k.csv"
        if stale_k_csv.exists():
            stale_k_csv.unlink()
        written_paths.append(path_obs_csv)

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
        ax.set_ylabel("Absolute magnitude, H")
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


def _set_fit_period(fit, p_rot, source):
    """Point later fits at one rotation period. p_rot is in hours."""
    freq = 48.0 / float(p_rot)
    fit["freq_sel"] = freq
    fit["per_sel"] = float(p_rot)
    fit["best_freq"] = freq
    fit["best_phot_period"] = float(p_rot) / 2.0
    fit["best_rot_period"] = float(p_rot)
    fit["period_source"] = source


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
    expected = [path for kind in kinds for path in _figure_paths(kind, provid, analysis)]
    csv_path = FIG_DIR / "asteroid_analysis_results.csv"
    csv_ok = False
    if csv_path.exists():
        try:
            df_csv = pd.read_csv(csv_path)
            csv_ok = (provid in df_csv["provid"].values) and (
                "epoch_utc_ini" in df_csv.columns
                and "obs_duration" in df_csv.columns
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
            if csv_ok and "periodogram_democratic_weighted" in analysis:
                if "p_rot_wdem_h" not in df_csv.columns:
                    csv_ok = False
            if csv_ok and "phase_model_compare" in analysis:
                if "p_rot_cmp_h" not in df_csv.columns:
                    csv_ok = False
            if csv_ok and (
                "multiband_photometry_compare" in analysis
                or "multiband_photometry_nightalias" in analysis
            ):
                for column in (
                    "p_rot_h",
                    "power",
                    "amp_ratio",
                    "amp_data",
                    "span_alias",
                    "g-r",
                ):
                    if column not in df_csv.columns:
                        csv_ok = False
            if csv_ok and "multiband_photometry_nightalias" in analysis:
                for column in (
                    "nightalias_tested",
                    "nightalias_stable",
                    "nightalias_band",
                ):
                    if column not in df_csv.columns:
                        csv_ok = False
            if csv_ok and "justfold" in analysis:
                if "p_rot_justfold_h" not in df_csv.columns:
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
            return None, "No photometric observations in the catalogue"
        if df_obj[PROP].notna().sum() == 0:
            print(f"  no {PROP} values after merging the geometry table")
            return (
                None,
                f"No {PROP} values after merging the geometry table (no MPC designation, or no matching epochs)",
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
                "periodogram_democratic_weighted",
                "phase",
                "phase_model",
                "phase_model_compare",
                "multiband_photometry_compare",
                "multiband_photometry_nightalias",
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
                    f"Fewer than 5 valid observations in band {band_use} with magnitude {PROP} ({n_band} observed)",
                )
            if fit["path"] is not None:
                written.append(fit["path"])
            periodogram_fold = {
                "freq_sel": float(fit["freq_sel"]),
                "per_sel": float(fit["per_sel"]),
                "source": fit.get("period_source")
                or f"Periodogram {fit.get('band', BAND)}",
            }
            dem = None
            if "periodogram_democratic" in analysis:
                dem = periodogram_democratic(
                    df_obj, provid, Nb=Nb, exclude=exclude_lx
                )
                if dem is not None:
                    fit["democratic"] = dem
                    if dem.get("path") is not None:
                        written.append(dem["path"])
            wdem = None
            if "periodogram_democratic_weighted" in analysis:
                wdem = periodogram_democratic_weighted(
                    df_obj, provid, Nb=Nb, exclude=exclude_lx
                )
                if wdem is not None:
                    fit["democratic_weighted"] = wdem
            if "phase" in analysis:
                written.append(
                    plot_phase_allbands(
                        df_obj,
                        provid,
                        fit,
                        kind="phase_periodogram",
                        freq_sel=periodogram_fold["freq_sel"],
                        per_sel=periodogram_fold["per_sel"],
                        source=periodogram_fold["source"],
                    )
                )
                if dem is not None:
                    written.append(
                        plot_phase_allbands(
                            df_obj,
                            provid,
                            fit,
                            kind="phase_democratic",
                            freq_sel=48.0 / float(dem["p_rot_h"]),
                            per_sel=float(dem["p_rot_h"]),
                            source=f"Democratic, unweighted, {dem['source']}",
                        )
                    )
                if wdem is not None:
                    written.append(
                        plot_phase_allbands(
                            df_obj,
                            provid,
                            fit,
                            kind="phase_democratic_weighted",
                            freq_sel=48.0 / float(wdem["p_rot_h"]),
                            per_sel=float(wdem["p_rot_h"]),
                            source=(
                                f"Democratic, weighted, {wdem['band_ref']}, "
                                f"k={wdem['k_alias']:+d}"
                            ),
                        )
                    )
                old_phase = FIG_DIR / f"{slug(provid)}_phase_allbands.webp"
                if old_phase.exists():
                    old_phase.unlink()
            adopted = dem if dem is not None else wdem
            compared = None
            if "phase_model_compare" in analysis:
                compared = phase_model_compare(
                    provid,
                    fit,
                    _vote_frequencies(periodogram_fold, dem, wdem),
                )
                fit["phase_model_compare"] = compared
                if compared["path"] is not None:
                    written.append(compared["path"])
            vote_merged = None
            if dem is not None and dem.get("candidates"):
                vote_merged = _merge_frequencies(
                    (cand["freq"], cand["label"]) for cand in dem["candidates"]
                )
            mb_compared = None
            if "multiband_photometry_compare" in analysis:
                if vote_merged is None:
                    print(
                        "  skip multiband_photometry_compare: "
                        "the democratic vote has no candidate periods"
                    )
                else:
                    mb_compared = multiband_photometry_compare(
                        df_obj, provid, vote_merged
                    )
                    if mb_compared is not None and mb_compared.get("no_two_peaks"):
                        if "multiband_photometry_nightalias" not in analysis:
                            return (
                                None,
                                mb_compared.get("reject_reason")
                                or "No two-peak period in the multiband comparison",
                            )
                    if mb_compared is not None and not mb_compared.get("no_two_peaks"):
                        fit["multiband_photometry_compare"] = mb_compared
                        if mb_compared["path"] is not None:
                            written.append(mb_compared["path"])
            mb_nightalias = None
            if "multiband_photometry_nightalias" in analysis:
                if vote_merged is None:
                    print(
                        "  skip multiband_photometry_nightalias: "
                        "the democratic vote has no candidate periods"
                    )
                else:
                    mb_nightalias = multiband_photometry_nightalias(
                        df_obj, provid, vote_merged
                    )
                    if mb_nightalias is not None and mb_nightalias.get("no_two_peaks"):
                        return (
                            None,
                            mb_nightalias.get("reject_reason")
                            or "No two-peak period in the multiband comparison",
                        )
                    if mb_nightalias is not None:
                        fit["multiband_photometry_nightalias"] = mb_nightalias
                        if mb_nightalias["path"] is not None:
                            written.append(mb_nightalias["path"])
            period_fit = None
            period_label = None
            if mb_nightalias is not None and not mb_nightalias.get("no_two_peaks"):
                period_fit = mb_nightalias
                period_label = (
                    f"Multiband 2-peak nightalias, chi2_nu = "
                    f"{mb_nightalias['chi2_nu']:.2f}"
                )
            elif mb_compared is not None and not mb_compared.get("no_two_peaks"):
                period_fit = mb_compared
                period_label = (
                    f"Multiband 2-peak compare, chi2_nu = "
                    f"{mb_compared['chi2_nu']:.2f}"
                )
            if period_fit is not None:
                _set_fit_period(fit, period_fit["p_rot_h"], period_label)
            elif compared is not None:
                _set_fit_period(
                    fit,
                    compared["p_rot_h"],
                    f"Phase-model compare, chi2_nu = {compared['chi2_nu']:.2f}",
                )
            elif adopted is not None:
                if dem is not None:
                    source = f"LS democratic {dem['source']}"
                else:
                    source = f"LS democratic weighted {wdem['band_ref']}"
                _set_fit_period(fit, adopted["p_rot_h"], source)
            if (
                period_fit is not None
                or compared is not None
                or adopted is not None
            ):
                stale = []
                if "phase_model" in analysis:
                    stale.append(output_path("phase_model", provid))
                if "timemodel" in analysis:
                    stale.append(output_path("timemodel", provid))
                if "timemodel_Lk" in analysis:
                    stale.append(output_path("timemodel_Lk", provid))
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

        if "justfold" in analysis:
            folded = justfold(df_obj, provid, df_rot)
            if folded is None:
                return (
                    None,
                    "The literature period could not be folded",
                )
            if fit is None:
                fit = {}
            fit["justfold"] = folded
            if folded["path"] is not None:
                written.append(folded["path"])

        multi_fit = None
        needs_lk = any(
            k in analysis
            for k in (
                "multiband_photometry",
                "shape_fit_multiband",
                "timemodel_Lk",
                "spectrum",
            )
        )
        if needs_lk:
            multi_fit = multiband_photometry(
                df_obj,
                provid,
                fit=fit,
                save="multiband_photometry" in analysis,
            )
            if multi_fit is None:
                return (
                    None,
                    "Not enough data for multiband photometry (fewer than 2 filters with data, or the fit did not converge)",
                )
            if fit is None:
                fit = {}
            fit["multiband"] = multi_fit
            written.extend(multi_fit.get("paths", []))

        if "timemodel_Lk" in analysis:
            compared = None
            if fit is not None:
                compared = fit.get("multiband_photometry_nightalias") or fit.get(
                    "multiband_photometry_compare"
                )
            if multi_fit is None or "df_band_k" not in (multi_fit or {}):
                return (
                    None,
                    "Lk light curve is missing, so the night-by-night model cannot be drawn",
                )
            if not compared or not np.isfinite(compared.get("p_rot_h", np.nan)):
                return (
                    None,
                    "The multiband period comparison did not choose a period",
                )
            path_lk = plot_timemodel_lk(
                provid, multi_fit["df_band_k"], compared["p_rot_h"]
            )
            if path_lk is not None:
                written.append(path_lk)

        if "shape_fit_multiband" in analysis:
            if multi_fit is None:
                return (
                    None,
                    "Multiband photometry fit is missing, so the shape cannot be calculated",
                )
            try:
                shape_k, paths_k = shape_fit_multiband(provid, multi_fit)
                if fit is None:
                    fit = {}
                fit["shape_k"] = shape_k
                written.extend(paths_k)
            except Exception as e:
                print(f"  error in shape_fit_multiband: {e}")
                return None, f"Error in the multiband shape fit: {e}"

        if "spectrum" in analysis:
            if multi_fit is None:
                return (
                    None,
                    "The joint light-curve fit is missing, so the spectrum cannot be calculated",
                )
            spec = spectrum(df_obj, provid, multi_fit, save=True)
            if spec is None:
                if fit is None and not written:
                    return (
                        None,
                        "No colors from the joint fit, so the spectrum cannot be calculated",
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
        return None, f"Unexpected error during processing: {e}"


def get_analysis_columns(analysis, rows=None):
    """Determine the ordered column schema for asteroid_analysis_results.csv."""
    cols = [
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
    ]
    if "periodogram" in analysis:
        cols.extend(["band_periodogram", "periods_h", "power", "amp_lsm", "amp_photo"])
    if "periodogram_democratic" in analysis:
        cols.extend(
            ["p_rot_dem_h", "votes_dem", "n_series_dem", "bands_dem", "periods_dem_h"]
        )
    if "periodogram_democratic_weighted" in analysis:
        cols.extend(
            [
                "p_rot_wdem_h",
                "k_alias_wdem",
                "score_wdem",
                "band_wdem",
                "scores_wdem",
            ]
        )
    if "phase_model_compare" in analysis:
        cols.extend(["p_rot_cmp_h", "chi2_nu_cmp", "periods_cmp"])
    if (
        "multiband_photometry_compare" in analysis
        or "multiband_photometry_nightalias" in analysis
    ):
        cols.extend(
            [
                "p_rot_h",
                "chi2_nu",
                "n_peaks",
                "power",
                "amp_ratio",
                "amp_data",
                "span_alias",
                *COMPARE_COLOR_COLUMNS,
            ]
        )
    if "multiband_photometry_nightalias" in analysis:
        cols.extend(
            [
                "nightalias_tested",
                "nightalias_stable",
                "nightalias_band",
                "nightalias_n_nights",
                "nightalias_n_tested",
                "nightalias_n_kept",
                "nightalias_periods",
            ]
        )
    if "justfold" in analysis:
        cols.extend(
            [
                "p_rot_justfold_h",
                "chi2_nu_justfold",
                "n_peaks_justfold",
                "power_justfold",
                "amp_ratio_justfold",
                "span_alias_justfold",
                "periods_justfold",
            ]
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
                *COMPARE_COLOR_COLUMNS,
            ]
        )
    if "shape_fit_multiband" in analysis:
        cols.extend([f"{c}_k" for c in SHAPE_COLUMNS])
    if "spectrum" in analysis:
        cols.extend(SPECTRUM_COLUMNS)

    cols = list(dict.fromkeys(cols))

    if rows:
        mandatory_set = {
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
        "k_alias_wdem",
        "n_peaks",
        "span_alias",
        "nightalias_tested",
        "nightalias_stable",
        "nightalias_n_nights",
        "nightalias_n_tested",
        "nightalias_n_kept",
        "n_peaks_justfold",
        "span_alias_justfold",
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
    src = ROOT / "asteroid_analysis_columns.md"
    if src.exists():
        dest = FIG_DIR / "asteroid_analysis_columns.md"
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
        "periodogram_democratic_weighted",
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


def exclude_finished(selected):
    """Drop asteroids already written in this run's results table.

    A new launch reads asteroid_analysis_results.csv in the figure folder and
    leaves those provids out of the loop.
    """
    path = FIG_DIR / "asteroid_analysis_results.csv"
    attrs = dict(getattr(selected, "attrs", {}))
    if not path.exists():
        print("Excluding 0 asteroids already present in asteroid_analysis_results.csv")
        return selected
    try:
        done = pd.read_csv(path, usecols=["provid"])
    except Exception as exc:
        print(f"Could not read {path.name}: {exc}")
        print("Excluding 0 asteroids already present in asteroid_analysis_results.csv")
        return selected
    done_ids = set(done["provid"].dropna().astype(str).str.strip())
    names = selected["provid"].astype(str).str.strip()
    n_exclude = int(names.isin(done_ids).sum())
    print(
        f"Excluding {n_exclude} asteroids already present in asteroid_analysis_results.csv"
    )
    kept = selected.loc[~names.isin(done_ids)].reset_index(drop=True)
    kept.attrs.update(attrs)
    return kept


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
        "periodogram_democratic",
        #"multiband_photometry_compare",
        "multiband_photometry_nightalias",
        #"lightcurve",
        #"periodogram",
        #"periodogram_democratic",
        #"periodogram_democratic_weighted",
        #"phase",
        #"phase",
        #"phase_model",
        #"phase_model_compare",
        #"shape_fit",
        #"multiband_photometry",
        #"multiband_photometry_compare",
        #"shape_fit_multiband",
        #"timemodel_Lk",
        #"spectrum",
    }
    if provids is None:
        provids = [
            #("2002 TH52", "Lr"),
            #("2025 MM82", "Lr"),
            #("2025 MA19", "Lr"),
            #("2002 PJ24", "Lx"),
            #("2025 MP71", "Lr"),
            #('2025 ML17', "Lr"),
            #('2000 TB35', "Lr"),
            ("1997 SX3", "Lr"),
            ("1998 BC26", "Lr"),
            ("1998 HH100", "Lr"),
        ]
    excluded_bands = ["Lu", "Ly"]
    selected, df_set, df_rot, analysis, missing = _selection_tables(
        provids, analysis, exclude_lx=excluded_bands
    )
    selected.attrs["Nb"] = int(Nb)
    return selected, df_set, df_rot, analysis, missing


def select_provid_rotation(Nb=3):
    """Provids in asteroids-with-rotation.csv, and the analyses to run."""
    fig_dir = ROOT / "figures/figures_rotation_parallel"
    set_fig_dir(fig_dir)
    analysis = {
        "periodogram_democratic",
        "multiband_photometry_compare",
        #"lightcurve",
        #"periodogram",
        #"periodogram_democratic",
        #"periodogram_democratic_weighted",
        #"phase",
        #"phase",
        #"phase_model",
        #"phase_model_compare",
        #"shape_fit",
        #"multiband_photometry",
        #"multiband_photometry_compare",
        #"shape_fit_multiband",
        #"timemodel_Lk",
        #"spectrum",
        "compare"
    }
    df_rot = load_rotation_catalog()
    provids = df_rot["provid"].tolist()
    # provids = provids[: len(provids) // 2]
    print(f"{len(provids)} of {len(df_rot)} asteroids with a literature period")
    selected, df_set, df_rot, analysis, missing = _selection_tables(provids, analysis)
    selected.attrs["Nb"] = int(Nb)
    return selected, df_set, df_rot, analysis, missing


def select_provid_rotfold():
    """Asteroids with a literature period. The only analysis is justfold."""
    fig_dir = ROOT / "figures/figures_rotfold"
    set_fig_dir(fig_dir)
    analysis = {"justfold"}
    df_rot = load_rotation_catalog()
    provids = df_rot["provid"].tolist()
    print(f"{len(provids)} of {len(df_rot)} asteroids with a literature period")
    return _selection_tables(provids, analysis)


def select_provid_nobs(n_min=100, n_max=100000, Nb=3):
    """Asteroids with n_obs between n_min and n_max, excluding the Medellin list."""
    fig_dir = ROOT / "figures/figures_nobs"
    set_fig_dir(fig_dir)
    analysis = {
        "periodogram_democratic",
        "multiband_photometry_compare",
        #"lightcurve",
        #"periodogram",
        #"periodogram_democratic",
        #"periodogram_democratic_weighted",
        #"phase",
        #"phase_model",
        #"shape_fit",
        #"timemodel",
        #"multiband_photometry",
        #"shape_fit_multiband",
        #"spectrum",
    }
    check_analysis(analysis, "analysis")

    df_rot = load_rotation_catalog()
    df_set = load_df_set()
    df_obs = observation_counts(df_set)
    selected = df_obs[
        (df_obs["n_obs"] >= n_min) & (df_obs["n_obs"] <= n_max)
    ].reset_index(drop=True)
    medellin = (
        pd.read_csv(MEDELLIN_CATALOG)
        .dropna(subset=["provid"])["provid"]
        .astype(str)
        .str.strip()
    )
    n_range = len(selected)
    selected = selected[~selected["provid"].isin(set(medellin))].reset_index(drop=True)
    selected["band"] = BAND
    print(
        f"{len(selected)} asteroids with {n_min} to {n_max} observations, "
        f"outside the Medellin list ({n_range - len(selected)} excluded)"
    )
    if selected.empty:
        raise SystemExit(
            f"No asteroids with {n_min} to {n_max} observations outside the Medellin list."
        )

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
        "periodogram_democratic",
        "multiband_photometry_compare",
        #"multiband_photometry_nightalias",
        #"lightcurve",
        #"periodogram",
        #"phase_model_compare",
        #"phase_model",
        #"shape_fit",
        #"timemodel",
        #"shape_fit_multiband",
        #"spectrum",
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


def select_provid_revisit(Nb=3):
    """Provids in revisit.csv: leave-one-night-out revisit into figures_nobs."""
    fig_dir = ROOT / "figures/figures_nobs"
    set_fig_dir(fig_dir)
    analysis = {
        "periodogram_democratic",
        "multiband_photometry_nightalias",
        #"multiband_photometry_compare",
        #"lightcurve",
        #"periodogram",
        #"phase_model_compare",
        #"phase_model",
        #"shape_fit",
        #"timemodel",
        #"shape_fit_multiband",
        #"spectrum",
    }
    if not REVISIT_CATALOG.exists():
        raise SystemExit(f"Missing revisit list: {REVISIT_CATALOG}")
    catalogue = pd.read_csv(REVISIT_CATALOG).dropna(subset=["provid"])
    if "band_periodogram" not in catalogue.columns:
        catalogue["band_periodogram"] = BAND
    provids = list(
        zip(
            catalogue["provid"].astype(str).str.strip(),
            catalogue["band_periodogram"].fillna(BAND).astype(str).str.strip(),
        )
    )
    print(f"Revisit list: {len(provids)} asteroids from {REVISIT_CATALOG.name}")
    selected, df_set, df_rot, analysis, missing = _selection_tables(provids, analysis)
    selected.attrs["Nb"] = int(Nb)
    return selected, df_set, df_rot, analysis, missing


def print_processing_report(total_requested, succeeded, failed, fig_dir=None, elapsed_s=None):
    """Print and save a summary report showing success and failure due to missing data."""
    n_succeeded = len(succeeded)
    n_failed = len(failed)
    total = total_requested
    pct_ok = (n_succeeded / total * 100.0) if total > 0 else 0.0
    pct_fail = (n_failed / total * 100.0) if total > 0 else 0.0

    lines = []
    lines.append("")
    lines.append("=" * 80)
    lines.append("                   FINAL ASTEROID PROCESSING REPORT")
    lines.append("=" * 80)
    lines.append(f"Asteroids considered : {total}")
    lines.append(f"  ✓ Completed successfully : {n_succeeded} ({pct_ok:.1f}%)")
    lines.append(f"  ✗ Problems from missing data : {n_failed} ({pct_fail:.1f}%)")
    if elapsed_s is not None:
        lines.append(f"  Total time : {format_elapsed(elapsed_s)}")
    lines.append("-" * 80)

    if failed:
        lines.append(f"ASTEROIDS WITH PROBLEMS ({n_failed}):")
        reasons_map = {}
        for provid, reason in failed.items():
            reasons_map.setdefault(reason, []).append(provid)

        for reason, provids in sorted(
            reasons_map.items(), key=lambda kv: len(kv[1]), reverse=True
        ):
            lines.append(f"\n  • Reason: {reason} ({len(provids)} asteroids)")
            lines.append(repr(sorted(provids)))
        lines.append("-" * 80)

    lines.append("=" * 80)
    lines.append("")
    print("\n".join(lines))

    if fig_dir is not None:
        file_lines = lines[:-2]
        file_lines.append(f"ASTEROIDS PROCESSED SUCCESSFULLY ({n_succeeded}):")
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
        # select_provid_rotfold: the same asteroids, folded at the given period
        # select_provid_medellin: medellin-vro-catalogue.csv
        # select_provid_nobs: 50 to 100 observations, outside the Medellin list
        # select_provid_revisit: data/revisit.csv (LONO) into figures_nobs
        # selected, df_set, df_rot, analysis, missing = unpack_selection(select_provid_all())
        # selected, df_set, df_rot, analysis, missing = unpack_selection(select_provid_list())
        selected, df_set, df_rot, analysis, missing = unpack_selection(
            #select_provid_rotation()
            # select_provid_rotfold()
            # select_provid_list()
            # select_provid_medellin()
            # select_provid_nobs()
            select_provid_revisit()
        )
        # selected, df_set, df_rot, analysis, missing = unpack_selection(select_provid_medellin())
    if not force:
        selected = exclude_finished(selected)
    designations = load_mpc_designations()

    total_requested = len(selected) + len(missing)
    succeeded = []
    failed = {
        p: "No photometric observations in obs_sbn_X05_full_mpc_desig.parquet"
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
                epoch_utc_ini, epoch_utc_end, obs_duration = observation_window(df_ast)
                n_obs_ast = int(len(df_ast))
                n_obs_bands = {
                    f"n_obs_{b}": int((df_ast["band"].isin([f"L{b}", b])).sum())
                    for b in ["u", "g", "r", "i", "z", "y"]
                }

                row = {
                    "provid": provid,
                    "epoch_utc_ini": epoch_utc_ini,
                    "epoch_utc_end": epoch_utc_end,
                    "obs_duration": obs_duration,
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

                if "periodogram_democratic_weighted" in analysis and fit.get(
                    "democratic_weighted"
                ):
                    wdem = fit["democratic_weighted"]
                    row["p_rot_wdem_h"] = wdem["p_rot_h"]
                    row["k_alias_wdem"] = wdem["k_alias"]
                    row["score_wdem"] = wdem["score"]
                    row["band_wdem"] = wdem["band_ref"]
                    row["scores_wdem"] = wdem["scores"]

                if "phase_model_compare" in analysis and fit.get(
                    "phase_model_compare"
                ):
                    compared = fit["phase_model_compare"]
                    row["p_rot_cmp_h"] = compared["p_rot_h"]
                    row["chi2_nu_cmp"] = compared["chi2_nu"]
                    row["periods_cmp"] = compared["periods"]

                mb_row = None
                if "multiband_photometry_nightalias" in analysis and fit.get(
                    "multiband_photometry_nightalias"
                ):
                    mb_row = fit["multiband_photometry_nightalias"]
                elif "multiband_photometry_compare" in analysis and fit.get(
                    "multiband_photometry_compare"
                ):
                    mb_row = fit["multiband_photometry_compare"]
                if mb_row is not None:
                    row["p_rot_h"] = mb_row["p_rot_h"]
                    row["chi2_nu"] = mb_row["chi2_nu"]
                    row["n_peaks"] = mb_row["n_peaks"]
                    row["power"] = mb_row["power"]
                    row["amp_ratio"] = mb_row["amp_ratio"]
                    row["amp_data"] = mb_row["amp_data"]
                    row["span_alias"] = mb_row["span_alias"]
                    row.update(mb_row.get("colors", {}))
                    if "multiband_photometry_nightalias" in analysis:
                        row["nightalias_tested"] = mb_row.get("nightalias_tested", 0)
                        row["nightalias_stable"] = mb_row.get("nightalias_stable", 0)
                        row["nightalias_band"] = mb_row.get("nightalias_band", "")
                        row["nightalias_n_nights"] = mb_row.get(
                            "nightalias_n_nights", 0
                        )
                        row["nightalias_n_tested"] = mb_row.get(
                            "nightalias_n_tested", 0
                        )
                        row["nightalias_n_kept"] = mb_row.get("nightalias_n_kept", 0)
                        row["nightalias_periods"] = mb_row.get("nightalias_periods", "")

                if "justfold" in analysis and fit.get("justfold"):
                    folded = fit["justfold"]
                    row["p_rot_justfold_h"] = folded["p_rot_h"]
                    row["chi2_nu_justfold"] = folded["chi2_nu"]
                    row["n_peaks_justfold"] = folded["n_peaks"]
                    row["power_justfold"] = folded["power"]
                    row["amp_ratio_justfold"] = folded["amp_ratio"]
                    row["span_alias_justfold"] = folded["span_alias"]
                    row["periods_justfold"] = folded["periods"]

                if "shape_fit" in analysis:
                    row.update(fit.get("shape", {}))

                if "multiband_photometry" in analysis:
                    if "multiband" in fit and fit["multiband"] is not None:
                        summary = dict(fit["multiband"].get("summary_row", {}))
                        if (
                            "multiband_photometry_compare" in analysis
                            or "multiband_photometry_nightalias" in analysis
                        ):
                            for key in COMPARE_COLOR_COLUMNS:
                                summary.pop(key, None)
                        row.update(summary)

                if "shape_fit_multiband" in analysis:
                    if "shape_k" in fit:
                        row.update(fit["shape_k"])

                if "spectrum" in analysis:
                    row.update(fit.get("spectrum") or {})

                peak_rows.append(row)
            else:
                failed[provid] = err or "Not enough data, or an error in the calculation"

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
