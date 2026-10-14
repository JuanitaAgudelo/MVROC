#!/usr/bin/env python
"""Repeat the rotation figures of search_spectra_rotation_orbit.ipynb.

One set of figures is written for every provid in data/asteroids-with-rotation.csv
(the notebook's df_rot).

Writes, for each designation, into figmass/:
  lightcurve_allnights_H_<provid>.png
  periodogram_<provid>.png
  phase_allbands_<provid>.png
  timemodel_H_<provid>.png
  phase_model_all_periods_<provid>.png

And two sample comparisons of the Lomb-Scargle period with the literature:
  period_compare_lsm.png
  period_compare_fourier.png
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

from util import ALL_BANDS, BANDS, band_colors

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
FIG_DIR = ROOT / "figmass"
ROTATION_CATALOG = DATA_DIR / "asteroids-with-rotation.csv"

BAND = "Lr"
PROP = "H"
PERIODS_H = (2 / 60, 24)  # photometric periods searched, in hours
N_FREQ = 20000


def slug(provid):
    return provid.replace(" ", "_")


def load_rotation_catalog():
    """df_rot from the notebook: asteroids-with-rotation.csv."""
    df_rot = pd.read_csv(ROTATION_CATALOG)
    df_rot["provid"] = df_rot["provid"].dropna().astype(str).str.strip()
    df_rot = df_rot[df_rot["provid"] != ""].drop_duplicates("provid")
    return df_rot.reset_index(drop=True)


def load_observations(provids):
    table = pq.read_table(
        DATA_DIR / "obs_sbn_X05_full.parquet",
        columns=["obstime", "permid", "provid", "ra", "dec", "band", "mag", "rmsmag"],
        filters=[("provid", "in", provids)],
    )
    df = table.to_pandas()
    df["band"] = df["band"].apply(
        lambda b: f"L{b}" if pd.notnull(b) and len(str(b)) == 1 else b
    )
    for col in ["ra", "dec", "mag", "rmsmag"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["obstime"] = pd.to_datetime(df["obstime"])
    df["night"] = df["obstime"].dt.floor("D")
    return df


def load_geometry(provid):
    table = pq.read_table(
        DATA_DIR / "obs_final_elements_geometry_H.parquet",
        columns=["obstime", "mag", "H"],
        filters=[("mpc_designation", "==", provid)],
    )
    df = table.to_pandas()
    df["obstime"] = pd.to_datetime(df["obstime"])
    return df


def prepare_object(df_obs, provid):
    df_obj = df_obs[df_obs["provid"] == provid].copy()
    if df_obj.empty:
        return df_obj
    df_obj = df_obj.sort_values(["band", "obstime"]).reset_index(drop=True)
    t0 = df_obj["obstime"].iloc[0]
    df_obj["t_days"] = (df_obj["obstime"] - t0).dt.total_seconds() / 86400.0
    df_spy = load_geometry(provid)
    df_obj = pd.merge(df_obj, df_spy, on="obstime", how="left", suffixes=("", "_spy"))
    return df_obj


def plot_lightcurve_h(df_obj, provid):
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
                sb["obstime"], sb["H"], yerr=sb["rmsmag"], fmt="o",
                label=f"{b} (N={len(sb)})", color=band_colors.get(b, "k"),
                markersize=4, alpha=0.85,
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
            sel = t_night["band"] == b
            if not sel.any():
                continue
            ax.errorbar(
                hours[sel], t_night.loc[sel, "H"], yerr=t_night.loc[sel, "rmsmag"],
                fmt="o", label=f"{b} (N={sel.sum()})", color=band_colors.get(b, "k"),
                markersize=4, alpha=0.85,
            )
        ax.set_title(f"Night {night_str}")
        ax.set_xlabel("Time since start [h]")
        ax.set_ylabel("H")
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), ncol=1, frameon=True)
        ax.tick_params(axis="x", rotation=30)

    fig.tight_layout()
    path = FIG_DIR / f"lightcurve_allnights_H_{slug(provid)}.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def periodogram(df_obj, provid):
    df_band = df_obj[df_obj["band"] == BAND].dropna(subset=[PROP]).copy()
    if len(df_band) < 5:
        print(f"  skip periodogram for {provid}: {len(df_band)} {BAND} points with {PROP}")
        return None

    t_band = df_band["t_days"].to_numpy()
    y_band = df_band[PROP].to_numpy()
    dy_band = df_band["rmsmag"].fillna(0.01).to_numpy()
    freq_grid = np.linspace(24 / PERIODS_H[1], 24 / PERIODS_H[0], N_FREQ)
    ls_band = LombScargle(t_band, y_band, dy_band)
    power = ls_band.power(freq_grid)

    main_peak_power = np.max(power)
    best_freq = freq_grid[np.argmax(power)]
    best_phot_period = 24.0 / best_freq
    best_rot_period = 2.0 * best_phot_period

    peaks, _ = find_peaks(power)
    high_power_peaks = np.intersect1d(peaks, np.where(power >= 0.8 * main_peak_power)[0])
    peak_data = []
    for idx in high_power_peaks:
        freq = freq_grid[idx]
        peak_data.append({
            "frequency": freq,
            "rotational_period_h": 2.0 * 24.0 / freq,
            "relative_power": power[idx] / main_peak_power,
        })
    if not peak_data:
        peak_data.append({
            "frequency": best_freq,
            "rotational_period_h": best_rot_period,
            "relative_power": 1.0,
        })
    df_peaks = pd.DataFrame(peak_data).sort_values("relative_power", ascending=False).reset_index(drop=True)
    all_rot_periods = df_peaks["rotational_period_h"].to_numpy()

    fig, ax = plt.subplots(figsize=(10, 6))
    period_grid_hours = (24.0 / freq_grid) * 2.0
    ax.plot(period_grid_hours, power, color="b", lw=1.5)
    ax.axvline(best_rot_period, color="r", linestyle="--", lw=2, label=f"Peak Period = {best_rot_period:.2f} h")
    for period in all_rot_periods:
        ax.axvline(period, color="g", linestyle=":", alpha=1, lw=1)
    ax.set_xlabel("Rotational Period [hours]")
    ax.set_ylabel("Spectral Power")
    ax.set_title(f"Lomb-Scargle Periodogram for Asteroid {provid} ({BAND} Band, property={PROP})")
    ax.set_xlim(0, 10)
    ax.legend(loc="upper right", frameon=True)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = FIG_DIR / f"periodogram_{slug(provid)}.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)

    freq_sel = df_peaks.iloc[0]["frequency"]
    per_sel = df_peaks.iloc[0]["rotational_period_h"]
    print(f"  {provid}: P_rot = {per_sel:.2f} h ({BAND}, {PROP})")
    return {
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
        "path": path,
    }


def plot_phase_allbands(df_obj, provid, fit):
    t0 = fit["t_band"].min()
    freq_sel = fit["freq_sel"]
    per_sel = fit["per_sel"]
    present = [b for b in ALL_BANDS if (df_obj["band"] == b).any()]
    m_max = {b: df_obj.loc[df_obj["band"] == b, PROP].min() for b in present}

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for b in present:
        sb = df_obj[df_obj["band"] == b]
        phase_b = np.mod((sb["t_days"].to_numpy() - t0) * freq_sel / 2.0, 1.0)
        dy = sb["rmsmag"].fillna(0.01)
        marker = "o" if b.startswith("L") else "s"
        m = m_max[b]
        style = dict(
            fmt=marker, ms=4, alpha=0.8, color=band_colors.get(b, "k"),
            label=f"{b} (N={len(sb)}, m_max={m:.2f})",
        )
        axes[0].errorbar(phase_b, sb[PROP], yerr=dy, **style)
        axes[1].errorbar(phase_b, sb[PROP] - m, yerr=dy, **style)

    axes[0].set_ylabel("Magnitude")
    axes[0].set_title("Observed")
    axes[1].set_ylabel(r"$m - m_{\mathrm{max}}$")
    axes[1].set_title("Normalized to the brightest point of each filter")
    for ax in axes:
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xlabel("Phase at the rotation period")
        ax.grid(True, alpha=0.3)
    axes[1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)
    fig.suptitle(f"Phase curve of {provid}, $P_{{rot}}$ = {per_sel:.2f} h")
    fig.tight_layout()
    path = FIG_DIR / f"phase_allbands_{slug(provid)}.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def plot_timemodel_h(provid, fit):
    df_band = fit["df_band"]
    nights_band = df_band["night"].drop_duplicates().sort_values()
    n_nights_band = len(nights_band)
    fig, axs = plt.subplots(n_nights_band, 1, figsize=(10, 2.6 * n_nights_band), sharey=True)
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
            hours, y_night, yerr=dy_night, fmt="o", ms=4, alpha=0.85,
            color="0.2", ecolor="0.6", label="Observed magnitude", zorder=3,
        )
        ax.plot(
            hours_grid, y_grid, color="C3", lw=1.8,
            label=f"Periodogram model, P = {fit['best_phot_period']:.2f} h", zorder=2,
        )
        ax.set_title(pd.to_datetime(night).strftime("%Y-%m-%d"))
        ax.set_ylabel("Magnitude")
        ax.grid(True, alpha=0.3)

    axs[0].invert_yaxis()
    axs[0].legend(loc="best", frameon=True)
    axs[-1].set_xlabel("Hours since the first observation of the night")
    fig.suptitle(f"{provid}: observed light curve and periodogram model")
    fig.tight_layout()
    path = FIG_DIR / f"timemodel_H_{slug(provid)}.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


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


def plot_phase_model_all_periods(provid, fit, lit_periods):
    """Folded light curve at our period and at the literature periods."""
    periods = {"VRO": fit["best_rot_period"], **lit_periods}
    t_band = fit["t_band"]
    y_band = fit["y_band"]
    dy_band = fit["dy_band"]
    ls_band = fit["ls_band"]
    t0 = t_band.min()
    curve_color = band_colors.get(BAND, "k")

    fig, axs = plt.subplots(len(periods), 1, figsize=(10, 4.5 * len(periods)), squeeze=False)
    for ax, (name, per) in zip(axs.ravel(), periods.items()):
        freq_sel = 24.0 / (per / 2.0)
        phase_rot = np.mod((t_band - t0) * freq_sel / 2.0, 1.0)
        phase_rot_grid = np.linspace(0, 1, 400)
        y_model_rot = ls_band.model(t0 + phase_rot_grid * 2.0 / freq_sel, freq_sel)
        ax.errorbar(
            phase_rot, y_band, yerr=dy_band, fmt="o", ms=4, alpha=0.75,
            color=curve_color, ecolor="0.6", label="Magnitude",
        )
        ax.plot(
            phase_rot_grid, y_model_rot, color=curve_color, lw=2,
            label=rf"Periodogram model, $P_{{rot}}$ = {per:.2f} h",
        )
        ax.set_xlabel("Phase at the rotation period ($2P$)")
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", frameon=True)
        ax.set_ylabel("Magnitude")
        ax.set_title(
            f"Folded light curve of {provid} ({PROP},{BAND}), "
            rf"$P_{{rot}}$ = {per:.2f} h [{name}]"
        )
    fig.tight_layout()
    path = FIG_DIR / f"phase_model_all_periods_{slug(provid)}.png"
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path


def plot_period_comparison(df_rot, measured):
    """Our Lomb-Scargle rotation period against the two literature periods."""
    specs = [
        ("period_lsm_hr", "LSM", "period_compare_lsm.png"),
        ("period_fourier_hr", "Fourier", "period_compare_fourier.png"),
    ]
    literature = df_rot.set_index("provid")
    paths = []
    for column, label, filename in specs:
        rows = []
        for provid, p_own in measured.items():
            if provid not in literature.index or column not in literature.columns:
                continue
            p_lit = pd.to_numeric(literature.at[provid, column], errors="coerce")
            if not np.isfinite(p_lit) or p_lit <= 0 or not np.isfinite(p_own) or p_own <= 0:
                continue
            rows.append((provid, float(p_lit), float(p_own)))
        if not rows:
            print(f"  no pairs for the {label} comparison")
            continue

        provids, x, y = zip(*rows)
        x = np.asarray(x)
        y = np.asarray(y)
        lo = min(x.min(), y.min()) * 0.7
        hi = max(x.max(), y.max()) * 1.4

        fig, ax = plt.subplots(figsize=(7.5, 7))
        ax.scatter(x, y, s=28, color="#1f4e79", zorder=3)
        for provid, xi, yi in rows:
            ax.annotate(provid, (xi, yi), textcoords="offset points", xytext=(4, 3), fontsize=6, alpha=0.85)
        ax.plot([lo, hi], [lo, hi], color="0.45", lw=1, label="1:1")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel(f"Literature period, {label} [h]")
        ax.set_ylabel("Lomb-Scargle period from our photometry [h]")
        ax.set_title(f"Rotation period: our photometry vs literature ({label})")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(loc="upper left", frameon=True)
        fig.tight_layout()
        path = FIG_DIR / filename
        fig.savefig(path, dpi=300)
        plt.close(fig)
        paths.append(path)
        print(f"  {path.name} ({len(rows)} asteroids)")
    return paths


def process(df_obs, provid, df_rot):
    print(f"{provid}")
    df_obj = prepare_object(df_obs, provid)
    if df_obj.empty:
        print("  no observations")
        return None
    if df_obj[PROP].notna().sum() == 0:
        print(f"  no {PROP} values after merging the geometry table")
        return None

    written = [plot_lightcurve_h(df_obj, provid)]
    fit = periodogram(df_obj, provid)
    if fit is not None:
        written.append(fit["path"])
        written.append(plot_phase_allbands(df_obj, provid, fit))
        written.append(plot_timemodel_h(provid, fit))
        written.append(plot_phase_model_all_periods(provid, fit, literature_periods(df_rot, provid)))
    for path in written:
        print(f"  {path.name}")
    return None if fit is None else fit["per_sel"]


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    df_rot = load_rotation_catalog()
    provids = df_rot["provid"].tolist()
    if not provids:
        raise SystemExit(f"No provids in {ROTATION_CATALOG}.")
    print(f"Loading observations for {len(provids)} provids from {ROTATION_CATALOG.name}")
    df_obs = load_observations(provids)
    measured = {}
    for provid in provids:
        period = process(df_obs, provid, df_rot)
        if period is not None:
            measured[provid] = period
    print("Period comparison with the literature")
    plot_period_comparison(df_rot, measured)


if __name__ == "__main__":
    main()
