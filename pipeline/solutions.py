#!/usr/bin/env python3
"""Folded light curves for mocs-retained rotation solutions.

Stage: solutions
Upstream: mocs (outputs/mocs/mocs.csv) + photometry (mois, states, crossmatches)
Downstream: shapefit

Usage (always from the pipeline/ directory):

    python solutions.py "1981 US22"
    python solutions.py 1981_US22
    python solutions.py all

Reads retained periods from outputs/mocs/mocs.csv and photometry via
mv.load_moi_photometry(). Writes outputs/solutions/: solutions.csv,
*_solution.webp and solutions.log. One CSV row per MVROC solution.
"""

import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargleMultiband

import foldings as fd
import mvroc as mv

# ---------------------------------------------------------------------------
# Solutions constants
# ---------------------------------------------------------------------------
OUT_DIR = "outputs/solutions"
RESULTS_CSV = f"{OUT_DIR}/solutions.csv"
LOG_FILE = f"{OUT_DIR}/solutions.log"
UPSTREAM_CSV = mv.MOCS_CSV

HARMONIC_ORDER = fd.HARMONIC_ORDER
FILT_ORDER = fd.FILT_ORDER
REF_FILT = fd.REF_FILT
BAND_ALIAS = fd.BAND_ALIAS
BAND_COLORS = fd.BAND_COLORS
COLOR_COLUMNS = fd.COLOR_COLUMNS

N_OBS_COLUMNS = [f"n_obs_{filt}" for filt in ("u", "g", "r", "i", "z", "y")]

RESULTS_COLUMNS = [
    "provid",
    "solution_id",
    "p_rot_h",
    "chi2_nu",
    "n_peaks",
    "power",
    "amp",
    "amp_ratio",
    "amp_data",
    "scatter",
    "span_alias",
    "alias_label",
    "filters",
    "n_filters",
    "ref_filter",
] + N_OBS_COLUMNS + COLOR_COLUMNS

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def solution_periods(solutions, periods_h):
    """Expand slash-separated mocs solution ids and periods.

    Parameters
    ----------
    solutions : str
        Slash-separated solution ids, e.g. "1/2".
    periods_h : str
        Slash-separated rotation periods in hours, e.g. "2.94/2.77".

    Returns
    -------
    list of (str, float)
        Pairs (solution_id, p_rot_h).

    Examples
    --------
    >>> solution_periods("1/2", "2.94/2.77")
    [('1', 2.94), ('2', 2.77)]
    """
    ids = [part.strip() for part in str(solutions).split("/") if part.strip()]
    periods = [
        float(part)
        for part in str(periods_h).split("/")
        if str(part).strip()
    ]
    if len(ids) != len(periods):
        raise ValueError(
            f"solutions/periods length mismatch: {ids} vs {periods}"
        )
    return list(zip(ids, periods))

def campaign_alias_for_period(p_rot_h, ts, tol=0.1):
    """Return (is_alias, label) if P_rot lands on a campaign-duration alias.

    Parameters
    ----------
    p_rot_h : float
        Rotation period in hours.
    ts : array-like
        Observation times in days.
    tol : float, optional
        Absolute frequency match in cycles/day.

    Returns
    -------
    is_alias : bool
    label : str
        Alias label, or "" when not an alias.

    Examples
    --------
    >>> campaign_alias_for_period(6.0, [0.0, 1.0])[0] in (True, False)
    True
    """
    freq = 48.0 / float(p_rot_h)
    for alias_freq, label in fd.baseline_alias_items(ts):
        if abs(freq - float(alias_freq)) <= tol:
            return True, str(label)
    return False, ""

def prepare_twopeak_frame(df_obj):
    """Build the multiband arrays for a two-peak harmonic fold.

    Parameters
    ----------
    df_obj : DataFrame
        Asteroid photometry from mv.prepare_asteroid.

    Returns
    -------
    dict or None
        Fit arrays and reference filter, or None if fewer than two filters.
    """
    # --> Drop Lu, reject magnitude outliers, choose the reference filter.
    df_fit = df_obj.dropna(subset=["H", "band"]).copy()
    df_fit["filt"] = df_fit["band"].map(lambda b: BAND_ALIAS.get(b, b))
    rms = pd.to_numeric(df_fit["rmsmag"], errors="coerce").to_numpy(float)
    df_fit["rmsmag"] = np.where(
        np.isfinite(rms) & (rms > 0), rms, mv.MIN_MAG_ERROR
    )
    df_fit = df_fit[df_fit["filt"].isin(FILT_ORDER)]
    kept = []
    n_rejected = 0
    for _, part in df_fit.groupby("filt", sort=False):
        good = fd.magnitude_inliers(part["H"].to_numpy())
        n_rejected += int((~good).sum())
        kept.append(part.loc[good])
    if kept:
        df_fit = pd.concat(kept, ignore_index=True)
    if n_rejected:
        print(f"  solutions: rejected {n_rejected} magnitude outlier(s)")
    present = [f for f in FILT_ORDER if (df_fit["filt"] == f).any()]
    if len(present) < 2:
        return None
    ref = REF_FILT
    if ref not in present:
        counts = {f: int((df_fit["filt"] == f).sum()) for f in present}
        ref = max(counts, key=counts.get)
    ts = df_fit["t_days"].to_numpy(float)
    ys = df_fit["H"].to_numpy(float)
    ss = np.maximum(df_fit["rmsmag"].to_numpy(float), 1e-4)
    filts = df_fit["filt"].to_numpy()
    return {
        "present": present,
        "ref": ref,
        "ts": ts,
        "ys": ys,
        "ss": ss,
        "filts": filts,
        "t0": float(np.min(ts)),
        "mask": np.ones(len(ys), bool),
        "ls_power": LombScargleMultiband(ts, ys, filts, ss),
    }

def score_twopeak_frequency(frame, freq, labels):
    """Fit the shared two-peak multiband model at one photometric frequency.

    Parameters
    ----------
    frame : dict
        Output of prepare_twopeak_frame.
    freq : float
        Photometric frequency in cycles/day.
    labels : list of str
        Labels attached to this candidate.

    Returns
    -------
    dict
        Scored model (chi2_nu, amp, colors inputs, folded curve, ...).
    """
    ts = frame["ts"]
    ys = frame["ys"]
    ss = frame["ss"]
    filts = frame["filts"]
    present = frame["present"]
    ref = frame["ref"]
    t0 = frame["t0"]
    mask = frame["mask"]
    K = HARMONIC_ORDER
    freq = float(freq)
    f_rot = freq / 2.0
    chi2, kpar = fd.chi2_fixed(
        ts, ys, ss, filts, present, ref, f_rot, K, mask, 0.0
    )
    X, names = fd.design_matrix(ts, filts, f_rot, K, present, ref)
    beta, _ = fd.wls_fit(X, ys, ss, mask)
    phase_grid = np.linspace(0, 1, 400, endpoint=False)
    Xg, _ = fd.design_matrix(
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
    y_on_ref = ys - np.array([offset.get(f, 0.0) for f in filts])
    y_curve = (
        fd.design_matrix(
            ts, np.full(len(ts), ref), f_rot, K, present, ref
        )[0]
        @ beta
    )
    amp = float(np.ptp(y_model))
    scatter = fd.scatter_about_curve(y_on_ref, y_curve, ss)
    ratio = amp / scatter if scatter > 0 else np.nan
    power = float(
        np.asarray(frame["ls_power"].power([freq]), dtype=float).ravel()[0]
    )
    return {
        "freq": freq,
        "f_rot": f_rot,
        "p_rot_h": 48.0 / freq,
        "labels": list(labels),
        "chi2": float(chi2),
        "chi2_nu": float(chi2) / max(int(mask.sum()) - kpar, 1),
        "n_peaks": fd.brightness_peaks(y_model),
        "power": power,
        "amp": amp,
        "scatter": scatter,
        "amp_ratio": ratio,
        "span_alias": int(fd.is_baseline_alias(labels)),
        "offset": offset,
        "beta": beta,
        "names": names,
        "cov": np.linalg.pinv(
            (X[mask] / ss[mask, None]).T @ (X[mask] / ss[mask, None])
        ),
        "y_model": y_model,
        "y_on_ref": y_on_ref,
        "phase_grid": phase_grid,
    }

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def mvroc_solution_plot(provid, frame, panels, filter_counts):
    """Vertical stack of folded two-peak panels for retained solutions.

    Parameters
    ----------
    provid : str
        Asteroid name.
    frame : dict
        Multiband fit frame from prepare_twopeak_frame.
    panels : list of dict
        Scored solutions to draw (one panel each).
    filter_counts : dict
        Number of points per filter.

    Returns
    -------
    str
        Path of the saved WebP figure.
    """
    # --> One full-width panel per retained MVROC solution.
    n_panel = len(panels)
    fig, axes = plt.subplots(
        n_panel,
        1,
        figsize=(10.5, 4.4 * n_panel),
        squeeze=False,
        sharex=True,
        layout="constrained",
    )
    present = frame["present"]
    ref = frame["ref"]
    ts = frame["ts"]
    ss = frame["ss"]
    filts = frame["filts"]
    t0 = frame["t0"]
    filters_text = ", ".join(
        f"{filt} (N={filter_counts[filt]})" for filt in present
    )

    photometry = np.concatenate(
        [np.asarray(item["y_on_ref"], dtype=float) for item in panels]
    )
    photometry = photometry[np.isfinite(photometry)]
    bright = float(np.min(photometry))
    faint = float(np.max(photometry))
    span = max(faint - bright, 0.05)
    pad = max(0.05, 0.15 * span)
    y_lo, y_hi = faint + pad, bright - pad

    for ax, item in zip(axes.ravel(), panels):
        # --> Fold on the rotation period and overlay the order-2 model.
        phase = np.mod((ts - t0) * item["f_rot"], 1.0)
        y_on_ref = item["y_on_ref"]
        beta_of = dict(zip(item["names"], item["beta"]))
        err_of = dict(
            zip(
                item["names"],
                np.sqrt(np.clip(np.diag(item["cov"]), 0.0, None)),
            )
        )
        for filt in present:
            sel = filts == filt
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
                yerr=ss[sel],
                fmt="o",
                ms=5.5,
                alpha=0.85,
                color=BAND_COLORS.get(f"L{filt}", BAND_COLORS.get(filt, "k")),
                ecolor="0.75",
                label=filt_label,
            )
        ax.plot(
            item["phase_grid"],
            item["y_model"],
            color="0.1",
            lw=2.2,
            label="2-peak harmonic model",
        )
        _, lo, hi = fd.data_amplitude(y_on_ref)
        if np.isfinite(lo) and np.isfinite(hi):
            ax.axhline(lo, color="0.25", ls=":", lw=1.2, label="data 2.5-97.5%")
            ax.axhline(hi, color="0.25", ls=":", lw=1.2)

        color_bits = []
        for key in (f"{a}-{b}" for a, b in fd.color_pairs()):
            if key in item["colors"] and f"{key}_err" in item["colors"]:
                color_bits.append(
                    rf"${key}$ = ${item['colors'][key]:+.3f}"
                    rf" \pm {item['colors'][key + '_err']:.3f}$"
                )
        color_line = "; ".join(color_bits[:4])
        alias_txt = (
            f", campaign alias ({item['alias_label']})"
            if item["span_alias"]
            else ""
        )
        title_color = "red" if item["span_alias"] else "black"
        ax.set_title(
            rf"MVROC sol. {item['solution_id']}: {filters_text}"
            "\n"
            rf"$P_{{\mathrm{{rot}}}}$ = {item['p_rot_h']:.4f} h, "
            rf"$\chi^{{2}}_{{\nu}}$ = {item['chi2_nu']:.2f}, "
            rf"$A$ = {item['amp']:.3f}, $A/s$ = {item['amp_ratio']:.2f}"
            f"{alias_txt}"
            + (f"\n{color_line}" if color_line else ""),
            fontsize=10,
            color=title_color,
        )
        ax.set_ylim(y_lo, y_hi)
        ax.set_xlim(0, 1)
        ax.set_ylabel("Absolute magnitude, H")
        ax.legend(loc="best", frameon=True, fontsize=7)
        ax.grid(True, alpha=0.3)

    axes.ravel()[-1].set_xlabel("Phase at the rotation period (2P)")
    fig.suptitle(f"Folding solutions of {provid} (H)", fontsize=12)

    path = f"{OUT_DIR}/{mv.slug(provid)}_solution.webp"
    mv.save_figure(fig, path)
    print(f"  solution: {path}")
    return path

# ---------------------------------------------------------------------------
# Science: one asteroid
# ---------------------------------------------------------------------------
def mvroc_solution_asteroid(df, provid, solutions, periods_h, qplot=True):
    """Fold every mocs-retained period for one asteroid.

    Parameters
    ----------
    df : DataFrame
        Full MOI photometry from mv.load_moi_photometry.
    provid : str
        Asteroid name.
    solutions : str
        Slash-separated solution ids from mocs.csv.
    periods_h : str
        Slash-separated periods in hours from mocs.csv.
    qplot : bool, optional
        If True (default), write the solution figure.

    Returns
    -------
    list of dict or None
        Result rows for solutions.csv, or None on failure.

    Examples
    --------
    >>> df = mv.load_moi_photometry()
    >>> rows = mvroc_solution_asteroid(df, "1981 US22", "1/2", "2.94/2.77")
    """
    print("=" * 20 + "\n" + provid + "\n" + "=" * 20)
    df_obj = mv.prepare_asteroid(df, provid)
    if df_obj.empty:
        print("  skip: no photometric observations in MOI photometry")
        return None
    print(f"  N = {len(df_obj)} observations")

    # --> Expand the retained solution list from mocs.csv.
    try:
        specs = solution_periods(solutions, periods_h)
    except ValueError as exc:
        print(f"  skip: {exc}")
        return None
    if not specs:
        print("  skip: no periods to fold")
        return None

    frame = prepare_twopeak_frame(df_obj)
    if frame is None:
        print("  skip: fewer than 2 filters for multiband fold")
        return None
    present = frame["present"]
    ref = frame["ref"]
    filter_counts = {filt: int((frame["filts"] == filt).sum()) for filt in present}

    # --> Score each retained period with the joint two-peak model.
    panels = []
    for sol_id, p_rot_h in specs:
        freq = 48.0 / float(p_rot_h)
        is_alias, alias_label = campaign_alias_for_period(p_rot_h, frame["ts"])
        labels = [f"MVROC solution {sol_id}"]
        if is_alias and alias_label:
            labels.append(alias_label)
        item = score_twopeak_frequency(frame, freq, labels)
        item["p_rot_h"] = float(p_rot_h)
        item["span_alias"] = int(is_alias)
        item["alias_label"] = alias_label if is_alias else ""
        item["solution_id"] = str(sol_id)
        item["colors"] = fd.compare_colors(df_obj, item, present, ref)
        panels.append(item)
        alias_mark = f", campaign alias ({alias_label})" if is_alias else ""
        print(
            f"  solution {sol_id}: P_rot = {p_rot_h:.6f} h, "
            f"chi2_nu = {item['chi2_nu']:.2f}, A = {item['amp']:.3f}, "
            f"A/s = {item['amp_ratio']:.2f}{alias_mark}"
        )

    if qplot:
        mvroc_solution_plot(provid, frame, panels, filter_counts)

    # --> One CSV row per retained MVROC solution.
    rows = []
    for item in panels:
        amp_data, _lo, _hi = fd.data_amplitude(item["y_on_ref"])
        row = {
            "provid": provid,
            "solution_id": item["solution_id"],
            "p_rot_h": float(item["p_rot_h"]),
            "chi2_nu": float(item["chi2_nu"]),
            "n_peaks": int(item["n_peaks"]),
            "power": float(item["power"]),
            "amp": float(item["amp"]),
            "amp_ratio": float(item["amp_ratio"]),
            "amp_data": float(amp_data),
            "scatter": float(item["scatter"]),
            "span_alias": int(item["span_alias"]),
            "alias_label": item["alias_label"],
            "filters": ",".join(present),
            "n_filters": len(present),
            "ref_filter": ref,
        }
        for filt in ("u", "g", "r", "i", "z", "y"):
            row[f"n_obs_{filt}"] = int(filter_counts.get(filt, 0))
        for col in COLOR_COLUMNS:
            row[col] = item["colors"].get(col, np.nan)
        rows.append(row)
    pd.DataFrame(rows, columns=RESULTS_COLUMNS).to_csv(
        RESULTS_CSV, mode="a", header=False, index=False
    )
    return rows

# ---------------------------------------------------------------------------
# Solutions driver
# ---------------------------------------------------------------------------
def mvroc_solutions():
    """Run folded mocs solutions from the command line.

    Must be launched from the pipeline/ directory.

    Stage: solutions
    Upstream: mocs
    Downstream: shapefit

    Command-line
    ------------
    Exactly one positional argument:

    * provid - one asteroid name, e.g. "1981 US22" or 1981_US22
    * all - process every asteroid listed in mocs.csv

    Examples
    --------
    python solutions.py "1981 US22"
    python solutions.py 1981_US22
    python solutions.py all

    Inputs
    ------
    outputs/mocs/mocs.csv
        Retained solution ids and periods (solutions, periods_h).
    outputs/mois/mois.parquet, outputs/states/states.parquet,
    outputs/crossmatches/crossmatches.csv
        MOI photometry with H (via mv.load_moi_photometry).

    Outputs
    -------
    outputs/solutions/solutions.csv
        One row per retained MVROC solution, with fold metrics and colors.
    outputs/solutions/<provid>_solution.webp
        Vertical stack of folded panels.
    outputs/solutions/solutions.log
        Same text as the terminal.

    Returns
    -------
    None
    """
    # --> Check command-line arguments.
    if len(sys.argv) != 2:
        print('Usage: python solutions.py "<provid>" | python solutions.py all')
        raise SystemExit(2)

    # --> Get the target provid.
    target = sys.argv[1].strip()
    if target.lower() != "all":
        target = target.replace("_", " ")

    # --> Create the output directory and set up logging.
    os.makedirs(OUT_DIR, exist_ok=True)
    log = open(LOG_FILE, "a", encoding="utf-8")
    # --> Mirror all print() output to the log file and to the terminal.
    sys.stdout = mv.LogOutput(sys.__stdout__, log)

    if not os.path.isfile(UPSTREAM_CSV):
        raise SystemExit(
            f"Missing upstream mocs CSV: {UPSTREAM_CSV}\n"
            "Run the mocs stage first (outputs/mocs/mocs.csv)."
        )
    # --> Initialize the results CSV if it doesn't exist.
    if not os.path.isfile(RESULTS_CSV):
        pd.DataFrame(columns=RESULTS_COLUMNS).to_csv(RESULTS_CSV, index=False)
        print(f"Created {RESULTS_CSV}")
    print(f"Output: {OUT_DIR}")
    print(f"Log: {LOG_FILE}")
    print(f"Upstream: {UPSTREAM_CSV}")

    # --> Asteroids to fold = those listed in mocs.csv.
    upstream = pd.read_csv(UPSTREAM_CSV)
    needed = {"provid", "solutions", "periods_h"}
    if not needed.issubset(upstream.columns):
        raise SystemExit(
            f"{UPSTREAM_CSV} must have columns provid, solutions, periods_h"
        )
    upstream["provid"] = upstream["provid"].astype(str).str.strip()
    upstream = upstream.drop_duplicates("provid", keep="last").reset_index(drop=True)

    # --> If not processing all asteroids, filter to the target.
    if target.lower() != "all":
        upstream = upstream[upstream["provid"] == target].reset_index(drop=True)
        if upstream.empty:
            raise SystemExit(f"Provid not found in {UPSTREAM_CSV}: {target}")

    print(f"{len(upstream)} asteroids to fold")
    # --> Load the MOI photometry once for the whole run.
    df = mv.load_moi_photometry()
    t0 = time.perf_counter()
    succeeded = []
    failed = []
    skipped = []

    for i, (_, row_up) in enumerate(upstream.iterrows(), start=1):
        provid = str(row_up["provid"])
        if mv.already_in_outputs(provid, RESULTS_CSV):
            print(f"{provid}: already in solutions.csv, skip")
            skipped.append(provid)
            continue
        # --> Fold every retained period for this asteroid.
        rows = mvroc_solution_asteroid(
            df,
            provid,
            row_up["solutions"],
            row_up["periods_h"],
        )
        if rows is not None:
            succeeded.append(provid)
        else:
            failed.append(provid)
        elapsed = time.perf_counter() - t0
        print(
            f"  elapsed: {elapsed / 60.0:.1f} m ({i}/{len(upstream)}), "
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
    mvroc_solutions()
