#!/usr/bin/env python
"""Compare rotation periods from figures_rotation with the literature catalogue.

Reads asteroid_analysis_results.csv and asteroids-with-rotation.csv, then
writes one log-log comparison per derived period and literature period, plus
a grid of all of them, into figures/figures_rotation/.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
FIG_DIR = ROOT / "figures" / "figures_rotation"
ANALYSIS_CSV = FIG_DIR / "asteroid_analysis_results.csv"
LITERATURE_CSV = ROOT / "data" / "asteroids-with-rotation.csv"

DPI = 150
WEBP_QUALITY = 80

# (column in the merged table, short filename, axis label)
OUR_PERIODS = (
    ("p_periodogram_h", "periodogram", "Periodogram"),
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


def load_periods():
    """Analysis periods joined to the literature periods on provid."""
    analysis = pd.read_csv(ANALYSIS_CSV)
    literature = pd.read_csv(LITERATURE_CSV)
    analysis["p_periodogram_h"] = pd.to_numeric(
        analysis["periods_h"].astype(str).str.split(",").str[0], errors="coerce"
    )
    keep = ["provid", "p_periodogram_h"] + [name for name, _, _ in OUR_PERIODS[1:]]
    lit_cols = ["provid"] + [name for name, _, _ in LITERATURE_PERIODS]
    merged = analysis[keep].merge(literature[lit_cols], on="provid", how="inner")
    print(f"{len(merged)} asteroids in both tables")
    return merged


def axis_limits(table):
    """Shared log limits covering every period drawn, with the 1:1 line diagonal."""
    columns = [name for name, _, _ in OUR_PERIODS] + [
        name for name, _, _ in LITERATURE_PERIODS
    ]
    values = pd.to_numeric(table[columns].stack(), errors="coerce")
    values = values[np.isfinite(values) & (values > 0)]
    lo = float(values.min()) * 0.7
    hi = float(values.max()) * 1.4
    return lo, hi


def _pairs(table, ours, literature):
    frame = table[["provid", ours, literature]].copy()
    frame[ours] = pd.to_numeric(frame[ours], errors="coerce")
    frame[literature] = pd.to_numeric(frame[literature], errors="coerce")
    ok = np.isfinite(frame[ours]) & np.isfinite(frame[literature])
    ok &= (frame[ours] > 0) & (frame[literature] > 0)
    return frame.loc[ok]


def draw_comparison(ax, pairs, x_col, y_col, x_label, y_label, limits, annotate):
    """Literature period on x, one derived period on y."""
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
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="upper left", frameon=True, fontsize=8)
    ratio = y / x
    n_match = int(np.sum((ratio > 0.8) & (ratio < 1.25)))
    n_double = int(np.sum((ratio > 1.6) & (ratio < 2.5)) + np.sum((ratio > 0.4) & (ratio < 0.625)))
    ax.set_title(f"{y_label} vs {x_label}\n{n_match}/{len(pairs)} within 25% of 1:1, {n_double} near 2:1 or 1:2")


def save_figure(fig, path):
    path = Path(path).with_suffix(".webp")
    fig.savefig(path, dpi=DPI, pil_kwargs={"quality": WEBP_QUALITY, "method": 6})
    plt.close(fig)
    print(f"  {path.name}")
    return path


def plot_individuals(table, limits):
    for y_col, y_stem, y_label in OUR_PERIODS:
        for x_col, x_stem, x_label in LITERATURE_PERIODS:
            pairs = _pairs(table, y_col, x_col)
            if pairs.empty:
                print(f"  no pairs for {y_label} vs {x_label}")
                continue
            fig, ax = plt.subplots(figsize=(8, 8))
            draw_comparison(ax, pairs, x_col, y_col, x_label, y_label, limits, annotate=True)
            fig.tight_layout()
            save_figure(fig, FIG_DIR / f"period_compare_{y_stem}_vs_{x_stem}.webp")


def plot_grid(table, limits):
    n_rows = len(OUR_PERIODS)
    n_cols = len(LITERATURE_PERIODS)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6.2 * n_cols, 5.4 * n_rows))
    for i, (y_col, _, y_label) in enumerate(OUR_PERIODS):
        for j, (x_col, _, x_label) in enumerate(LITERATURE_PERIODS):
            ax = axes[i, j]
            pairs = _pairs(table, y_col, x_col)
            if pairs.empty:
                ax.set_axis_off()
                continue
            draw_comparison(
                ax, pairs, x_col, y_col, x_label, y_label, limits, annotate=False
            )
    fig.tight_layout()
    save_figure(fig, FIG_DIR / "period_compare_grid.webp")


def main():
    table = load_periods()
    limits = axis_limits(table)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    plot_individuals(table, limits)
    plot_grid(table, limits)


if __name__ == "__main__":
    main()
