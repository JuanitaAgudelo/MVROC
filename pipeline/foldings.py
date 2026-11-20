#!/usr/bin/env python3
"""Multiband two-peak folding of democratic period candidates.

Stage: foldings
Upstream: periodograms (outputs/periodograms/periodograms.csv)
Downstream: shapefit

Usage (always from the pipeline/ directory):

    python foldings.py "1981 US22"
    python foldings.py 1981_US22
    python foldings.py all

Upstream: periodograms (and through it crossmatches, orbits, states, mois).
Reads photometry via mv.load_moi_photometry() and voted periods from
outputs/periodograms/periodograms.csv. Writes outputs/foldings/:
foldings.csv, *_folding.webp and foldings.log. One CSV row per plotted
two-peak solution, with color columns at the end.
"""

import os
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargleMultiband
from scipy.signal import find_peaks

import mvroc as mv

# ---------------------------------------------------------------------------
# Folding constants
# ---------------------------------------------------------------------------
OUT_DIR = "outputs/foldings"
RESULTS_CSV = f"{OUT_DIR}/foldings.csv"
LOG_FILE = f"{OUT_DIR}/foldings.log"
UPSTREAM_CSV = "outputs/periodograms/periodograms.csv"

PERIODS_H = (0.5 / 60, 24)
ALIAS_P_MIN_H = 2.0
HARMONIC_ORDER = 2
FILT_ORDER = ["g", "r", "i", "z", "y"]
REF_FILT = "r"
BAND_ALIAS = {"Lg": "g", "Lr": "r", "Li": "i", "Lz": "z", "Ly": "y", "Lu": "u"}
BAND_COLORS = {
    "Lu": "#0000ff",
    "Lg": "#00bfbf",
    "Lr": "#008000",
    "Li": "#bfbf00",
    "Lz": "#ff0000",
    "Ly": "#bf00bf",
    "u": "#5252ff",
    "g": "#12ffff",
    "r": "#00d200",
    "i": "#ffff12",
    "z": "#ff5252",
    "y": "#ff12ff",
}
COLOR_FILTERS = ("u", "g", "r", "i", "z", "y")
COLOR_COLUMNS = [
    col
    for i in range(len(COLOR_FILTERS))
    for j in range(i + 1, len(COLOR_FILTERS))
    for col in (f"{COLOR_FILTERS[i]}-{COLOR_FILTERS[j]}",
                f"{COLOR_FILTERS[i]}-{COLOR_FILTERS[j]}_err")
]

RESULTS_COLUMNS = [
    "provid",
    "selected",
    "p_rot_fold_h",
    "chi2_nu",
    "n_peaks",
    "power",
    "amp_ratio",
    "amp_data",
    "span_alias",
] + COLOR_COLUMNS

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def candidates_from_periods_dem_h(periods_dem_h):
    """Build merged vote candidates from periodograms.csv periods_dem_h.

    Parameters
    ----------
    periods_dem_h : str
        Comma-separated label=period_h tokens, e.g. "Lr=2.94,Lg=2.77".

    Returns
    -------
    list of dict
        Merged vote frequencies with keys freq and labels.

    Examples
    --------
    >>> cands = candidates_from_periods_dem_h("Lr=2.94,Lg=2.77")
    >>> len(cands), sorted(cands[0]["labels"] + cands[1]["labels"])
    (2, ['Lg', 'Lr'])
    """
    items = []
    for part in str(periods_dem_h).split(","):
        part = part.strip().strip('"')
        if not part or "=" not in part:
            continue
        label, value = part.split("=", 1)
        p_rot_h = float(value)
        if not np.isfinite(p_rot_h) or p_rot_h <= 0:
            continue
        items.append((48.0 / p_rot_h, label.strip()))
    return mv.merge_frequencies(items)

def magnitude_inliers(values, n_iter=5):
    """True for magnitudes that belong to the main light curve of one band.

    Parameters
    ----------
    values : array-like
        Magnitudes of one filter.
    n_iter : int, optional
        Outlier-rejection passes.

    Returns
    -------
    ndarray of bool
        True for inliers.

    Examples
    --------
    >>> magnitude_inliers([12.0, 12.1, np.nan]).tolist()
    [True, True, False]
    """
    values = np.asarray(values, dtype=float)
    good = np.isfinite(values)
    if int(good.sum()) < 8:
        return good
    for _ in range(n_iter):
        sample = values[good]
        low, high = np.percentile(sample, [5, 95])
        margin = max(0.2, 0.35 * float(high - low))
        updated = (
            np.isfinite(values)
            & (values >= low - margin)
            & (values <= high + margin)
        )
        if int(updated.sum()) < 8 or np.array_equal(updated, good):
            break
        good = updated
    return good

def design_matrix(ts, filts, f_base, K, filters, ref):
    """Harmonic design matrix with one zero point per non-reference filter.

    Parameters
    ----------
    ts : array-like
        Times in days.
    filts : array-like
        Filter label per point.
    f_base : float
        Rotational frequency in cycles/day.
    K : int
        Harmonic order.
    filters : list of str
        Filters in the joint fit.
    ref : str
        Reference filter (zero point fixed at 0).

    Returns
    -------
    X : ndarray
        Design matrix.
    names : list of str
        Coefficient names matching the columns of X.

    Examples
    --------
    >>> X, names = design_matrix([0.0, 0.1], ["r", "g"], 8.0, 2, ["g", "r"], "r")
    >>> "m_ref" in names and "c_g" in names
    True
    """
    cols = [(filts == b).astype(float) for b in filters if b != ref] + [
        np.ones(len(ts))
    ]
    names = [f"c_{b}" for b in filters if b != ref] + ["m_ref"]
    ph = 2 * np.pi * f_base * ts
    for k in range(1, K + 1):
        cols += [np.sin(k * ph), np.cos(k * ph)]
        names += [f"sin{k}", f"cos{k}"]
    return np.column_stack(cols), names

def wls_fit(X, ys, ss, mask):
    """Weighted least-squares fit and residual vector.

    Parameters
    ----------
    X : ndarray
        Design matrix.
    ys : array-like
        Magnitudes.
    ss : array-like
        Magnitude uncertainties.
    mask : array-like of bool
        Points used in the fit.

    Returns
    -------
    beta : ndarray
        Fitted coefficients.
    residual : ndarray
        ys - X @ beta for every row.

    Examples
    --------
    >>> X = np.ones((3, 1)); ys = np.array([1.0, 2.0, 3.0])
    >>> beta, r = wls_fit(X, ys, np.ones(3), np.ones(3, bool))
    >>> float(beta[0])
    2.0
    """
    beta, *_ = np.linalg.lstsq(
        X[mask] / ss[mask, None], ys[mask] / ss[mask], rcond=None
    )
    return beta, ys - X @ beta

def chi2_fixed(ts, ys, ss, filts, present, ref, f_base, K, mask, jitter):
    """Chi-squared of a fixed-frequency harmonic model.

    Parameters
    ----------
    ts, ys, ss, filts : array-like
        Times, magnitudes, errors and filter labels.
    present : list of str
        Filters in the joint fit.
    ref : str
        Reference filter.
    f_base : float
        Rotational frequency in cycles/day.
    K : int
        Harmonic order.
    mask : array-like of bool
        Points used in the fit.
    jitter : float
        Extra variance added in quadrature to ss.

    Returns
    -------
    chi2 : float
        Weighted sum of squared residuals.
    n_par : int
        Number of free coefficients.

    Examples
    --------
    >>> chi2, n_par = chi2_fixed(
    ...     [0.0, 0.1], [12.0, 12.1], [0.01, 0.01], ["r", "r"],
    ...     ["r"], "r", 8.0, 2, np.ones(2, bool), 0.0)
    >>> n_par >= 1 and chi2 >= 0
    True
    """
    X, _ = design_matrix(ts, filts, f_base, K, present, ref)
    s = np.sqrt(ss**2 + jitter**2)
    _, r = wls_fit(X, ys, s, mask)
    return np.sum((r[mask] / s[mask]) ** 2), X.shape[1]

def brightness_peaks(y_model):
    """Number of brightness maxima of one rotation (magnitude minima).

    Parameters
    ----------
    y_model : array-like
        Model magnitudes over one rotational cycle.

    Returns
    -------
    int
        Number of brightness peaks (0, 1 or 2 for order-2 models).

    Examples
    --------
    >>> phase = np.linspace(0, 1, 200, endpoint=False)
    >>> brightness_peaks(np.sin(4 * np.pi * phase))
    2
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

def is_baseline_alias(labels):
    """True when a candidate is a daily alias of the campaign duration.

    Parameters
    ----------
    labels : iterable of str
        Vote / alias labels for one frequency.

    Returns
    -------
    bool

    Examples
    --------
    >>> is_baseline_alias(["Lr", "span k=+3"])
    True
    >>> is_baseline_alias(["Lr", "Lg"])
    False
    """
    return any(str(label).startswith(("span", "2-span")) for label in labels)

def scatter_about_curve(y_points, y_curve, sigma, trim=0.10):
    """Mean distance of the points to the model curve, errors included.

    Parameters
    ----------
    y_points : array-like
        Observed magnitudes shifted onto the reference filter.
    y_curve : array-like
        Model magnitudes at the same times.
    sigma : array-like
        Magnitude uncertainties.
    trim : float, optional
        Fraction of farthest points dropped before the mean.

    Returns
    -------
    float
        Robust scatter s used in the A/s ratio.

    Examples
    --------
    >>> scatter_about_curve([0.0, 0.1], [0.0, 0.0], [0.01, 0.01]) > 0
    True
    """
    distance = np.abs(
        np.asarray(y_points, dtype=float) - np.asarray(y_curve, dtype=float)
    )
    sigma = np.maximum(np.asarray(sigma, dtype=float), 0.0)
    n = len(distance)
    if n == 0:
        return np.nan
    n_drop = int(np.floor(trim * n))
    n_keep = max(n - n_drop, 1)
    keep = np.argsort(distance)[:n_keep]
    width = np.sqrt(distance[keep] ** 2 + sigma[keep] ** 2)
    return float(np.mean(width))

def data_amplitude(ys, low=2.5, high=97.5):
    """Peak-to-peak of the magnitudes between two percentiles.

    Parameters
    ----------
    ys : array-like
        Magnitudes.
    low, high : float, optional
        Percentiles that define the bright and faint edges.

    Returns
    -------
    amp, bright, faint : float
        Amplitude and the two percentile levels.

    Examples
    --------
    >>> amp, bright, faint = data_amplitude([10.0, 10.5, 11.0, 11.5, 12.0])
    >>> amp > 0 and bright < faint
    True
    """
    y = np.asarray(ys, dtype=float)
    y = y[np.isfinite(y)]
    if len(y) < 2:
        return np.nan, np.nan, np.nan
    faint, bright = np.percentile(y, [high, low])
    return float(faint - bright), float(bright), float(faint)

def campaign_alias_k_range(span_days, p_min_h=ALIAS_P_MIN_H):
    """Integer daily shifts for campaign-duration aliases.

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
    >>> len(list(campaign_alias_k_range(10.0))) > 0
    True
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

def baseline_alias_items(ts, p_min_h=ALIAS_P_MIN_H):
    """Daily aliases of the campaign duration, in cycles per day.

    Parameters
    ----------
    ts : array-like
        Observation times in days.
    p_min_h : float, optional
        Shortest rotation period (hours) kept as an alias.

    Returns
    -------
    list of (float, str)
        Frequency and label pairs (span / 2-span).

    Examples
    --------
    >>> items = baseline_alias_items(np.linspace(0, 10, 50))
    >>> isinstance(items, list)
    True
    """
    t = np.asarray(ts, dtype=float)
    span = float(np.ptp(t))
    if not np.isfinite(span) or span <= 0:
        return []
    f_min = 24.0 / PERIODS_H[1]
    f_max = 24.0 / PERIODS_H[0]
    p_min = mv.cadence_min_period(t)
    items = []
    for base, tag in ((1.0 / span, "span"), (2.0 / span, "2-span")):
        for k in campaign_alias_k_range(span, p_min_h):
            freq = base + float(k)
            if freq < f_min or freq > f_max:
                continue
            if p_min > 0 and (24.0 / freq) < p_min:
                continue
            items.append((freq, f"{tag} k={k:+d}"))
    return items

def tag_campaign_aliases(merged, ts, tol=0.1, p_min_h=ALIAS_P_MIN_H):
    """Append span / 2-span labels to voted frequencies that land on an alias.

    Parameters
    ----------
    merged : list of dict
        Vote groups with keys freq and labels.
    ts : array-like
        Observation times in days.
    tol : float, optional
        Absolute frequency match in cycles/day.
    p_min_h : float, optional
        Shortest rotation period (hours) for alias tags.

    Returns
    -------
    list of dict
        Same groups, possibly with extra alias labels.

    Examples
    --------
    >>> tag_campaign_aliases([{"freq": 5.0, "labels": ["Lr"]}], [0.0, 1.0])
    [{'freq': 5.0, 'labels': ['Lr']}]
    """
    alias_items = baseline_alias_items(ts, p_min_h)
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

def color_pairs(filters=COLOR_FILTERS):
    """Every color a-b with a bluer than b, in wavelength order.

    Parameters
    ----------
    filters : sequence of str, optional
        Filters from blue to red.

    Returns
    -------
    list of (str, str)
        Color pairs.

    Examples
    --------
    >>> ("g", "r") in color_pairs(("g", "r", "i"))
    True
    """
    return [
        (filters[i], filters[j])
        for i in range(len(filters))
        for j in range(i + 1, len(filters))
    ]

def color_difference(a, b, color, sig_color, cov, names, ref):
    """Magnitude a minus magnitude b, and its uncertainty, from the joint fit.

    Parameters
    ----------
    a, b : str
        Filter names (a bluer than b).
    color : dict
        Zero point per filter relative to ref.
    sig_color : dict
        Formal error of each free zero point.
    cov : ndarray
        Coefficient covariance matrix.
    names : list of str
        Coefficient names matching cov.
    ref : str
        Reference filter.

    Returns
    -------
    value, err : float
        Color and its uncertainty.

    Examples
    --------
    >>> color_difference("g", "r", {"g": 0.5, "r": 0.0}, {"g": 0.01, "r": 0.0},
    ...                  np.eye(1), ["c_g"], "r")
    (0.5, 0.01)
    """
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

def u_minus_reference(df_obj, beta, f_rot, K, present, ref):
    """u minus the reference curve. u stays out of the joint harmonic fit.

    Parameters
    ----------
    df_obj : DataFrame
        Asteroid photometry including band u when present.
    beta : ndarray
        Coefficients of the joint non-u fit.
    f_rot : float
        Rotational frequency in cycles/day.
    K : int
        Harmonic order of that fit.
    present : list of str
        Filters in the joint fit.
    ref : str
        Reference filter.

    Returns
    -------
    (float, float) or None
        Mean u-minus-reference offset and its uncertainty, or None if too
        few u points.

    Examples
    --------
    >>> u_minus_reference(pd.DataFrame(columns=["H", "band", "t_days", "rmsmag"]),
    ...                   np.zeros(5), 8.0, 2, ["g", "r"], "r") is None
    True
    """
    # --> Keep only u-band points and compare them to the reference model.
    df_u = df_obj.dropna(subset=["H", "band"]).copy()
    df_u["filt"] = df_u["band"].map(lambda b: BAND_ALIAS.get(b, b))
    df_u = df_u[df_u["filt"] == "u"]
    if df_u.empty or "t_days" not in df_u.columns:
        return None
    rms = pd.to_numeric(df_u["rmsmag"], errors="coerce").to_numpy(float)
    t_u = df_u["t_days"].to_numpy(float)
    y_u = df_u["H"].to_numpy(float)
    s_u = np.where(np.isfinite(rms) & (rms > 0), rms, mv.MIN_MAG_ERROR)
    X_ref, _ = design_matrix(
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
    chi2_red = float(
        np.sum(((d_u[good] - c_u) / s_u[good]) ** 2) / max(int(good.sum()) - 1, 1)
    )
    err_u = float(np.sqrt(1.0 / np.sum(w[good]) * max(chi2_red, 1.0)))
    return c_u, err_u

def compare_colors(df_obj, best, present, ref):
    """Color indices of one two-peak model, including u when available.

    Parameters
    ----------
    df_obj : DataFrame
        Asteroid photometry (needed for the post-fit u offset).
    best : dict
        Scored solution with beta, names, cov and f_rot.
    present : list of str
        Filters in the joint fit.
    ref : str
        Reference filter.

    Returns
    -------
    dict
        Color values and *_err uncertainties (keys like g-r, g-r_err).

    Examples
    --------
    >>> isinstance(compare_colors(pd.DataFrame(), {"beta": np.zeros(5),
    ...     "names": ["m_ref", "sin1", "cos1", "sin2", "cos2"],
    ...     "cov": np.eye(5), "f_rot": 8.0}, ["r"], "r"), dict)
    True
    """
    # --> Zero points of the joint fit, then every bluer-minus-redder pair.
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
            value, err = color_difference(a, b, color, sig_color, cov, names, ref)
            colors[f"{a}-{b}"] = value
            colors[f"{a}-{b}_err"] = err
    # --> Attach u-minus-filter colors from the post-fit u offset.
    u_ref = u_minus_reference(df_obj, beta, best["f_rot"], 2, present, ref)
    if u_ref is not None:
        c_u, err_u = u_ref
        for filt in joint:
            colors[f"u-{filt}"] = c_u - color[filt]
            colors[f"u-{filt}_err"] = float(
                np.sqrt(err_u**2 + sig_color.get(filt, 0.0) ** 2)
            )
    return colors

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def mvroc_folding_plot(
    provid,
    ts,
    ys_err,
    filts,
    present,
    ref,
    best,
    others,
    n_obs_total,
    alias_remnant,
):
    """Draw adopted and other two-peak folded panels; save WebP.

    Parameters
    ----------
    provid : str
        Asteroid name.
    ts, ys_err, filts : array-like
        Fit times, errors and filter labels.
    present : list of str
        Filters used in the fit.
    ref : str
        Reference filter.
    best : dict
        Adopted scored period.
    others : list of dict
        Other two-peak voted periods.
    n_obs_total : int
        Total observations for the title.
    alias_remnant : bool
        Mark the title when a better alias was rejected.

    Returns
    -------
    str
        Path of the saved figure.
    """
    # --> Layout: selected panel on top, other two-peak solutions below.
    n_other = len(others)
    nrow_other = int(np.ceil(n_other / 2.0)) if n_other else 0
    fig = plt.figure(
        figsize=(10.5, 5.8 + 2.9 * nrow_other),
        layout="constrained",
    )
    t0 = float(np.min(ts))
    s_fit = np.asarray(ys_err, dtype=float)
    f_fit = np.asarray(filts)

    def draw_panel(ax, item, selected):
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
                color=BAND_COLORS.get(f"L{filt}", BAND_COLORS.get(filt, "k")),
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
            _, bright, faint = data_amplitude(y_on_ref)
            if np.isfinite(bright) and np.isfinite(faint):
                ax.axhline(
                    bright, color="0.25", ls=":", lw=1.2, label="data 2.5-97.5%"
                )
                ax.axhline(faint, color="0.25", ls=":", lw=1.2)
        mark = ", selected" if selected else ""
        title_color = "black" if selected else (
            "red" if item["span_alias"] else "black"
        )
        show_remnant = bool(selected and alias_remnant)
        ax.set_title(
            rf"$P_{{rot}}$ = {item['p_rot_h']:.2f} h, "
            rf"$N$ = {n_obs_total}, "
            rf"$\chi^2_\nu$ = {item['chi2_nu']:.1f}, "
            f"{item['n_peaks']} peaks{mark}\n"
            rf"power = {item['power']:.3f}, "
            rf"$A$ = {item['amp']:.3f}, $A/s$ = {item['amp_ratio']:.2f}"
            "\n" + ", ".join(item["labels"]),
            fontsize=11 if selected else 8,
            color=title_color,
            pad=22 if show_remnant else 6,
        )
        if show_remnant:
            ax.text(
                0.5,
                1.0,
                "alias remnant",
                transform=ax.transAxes,
                color="red",
                fontsize=11,
                fontweight="bold",
                ha="center",
                va="bottom",
                clip_on=False,
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
        ax_best.set_xlabel("Phase at the rotation period (2P)")
    else:
        bottom = {}
        for i, ax in enumerate(other_axes):
            bottom[i % 2] = ax
            if i % 2 == 0:
                ax.set_ylabel("Absolute magnitude, H")
            else:
                ax.tick_params(labelleft=False)
        for ax in bottom.values():
            ax.set_xlabel("Phase at the rotation period (2P)")
    fig.suptitle(f"Multiband 2-peak comparison of {provid} (H)", fontsize=12)

    path = f"{OUT_DIR}/{mv.slug(provid)}_folding.webp"
    mv.save_figure(fig, path)
    print(f"  folding: {path}")
    return path

# ---------------------------------------------------------------------------
# Science: one asteroid
# ---------------------------------------------------------------------------
def mvroc_folding_asteroid(df, provid, candidates, qplot=True):
    """Fit a two-peak multiband model at each voted period; keep the best.

    Candidates come from the democratic vote (upstream periodograms). Each is
    tagged if it lands on a campaign alias. Only non-alias two-peak models can
    be adopted. u is left out of the joint fit. One CSV row is written per
    plotted two-peak solution (selected plus the other panels).

    Parameters
    ----------
    df : DataFrame
        Full MOI photometry from mv.load_moi_photometry.
    provid : str
        Asteroid name.
    candidates : list of dict
        Merged vote frequencies with keys freq and labels.
    qplot : bool, optional
        If True (default), write the folding figure.

    Returns
    -------
    list of dict or None
        Result rows for foldings.csv, or None on failure.

    Examples
    --------
    >>> df = mv.load_moi_photometry()
    >>> cands = candidates_from_periods_dem_h("Lr=2.94,Lg=2.77")
    >>> rows = mvroc_folding_asteroid(df, "1981 US22", cands)
    """
    print("=" * 20 + "\n" + provid + "\n" + "=" * 20)
    df_obj = mv.prepare_asteroid(df, provid)
    if df_obj.empty:
        print("  skip: no photometric observations in MOI photometry")
        return None
    print(f"  N = {len(df_obj)} observations")
    if not candidates:
        print("  skip: no democratic vote candidates from periodograms")
        return None

    # --> Prepare multiband frame (drop Lu, reject magnitude outliers).
    n_obs_total = int(len(df_obj))
    df_fit = df_obj.dropna(subset=["H", "band"]).copy()
    df_fit["filt"] = df_fit["band"].map(lambda b: BAND_ALIAS.get(b, b))
    rms = pd.to_numeric(df_fit["rmsmag"], errors="coerce").to_numpy(float)
    df_fit["rmsmag"] = np.where(
        np.isfinite(rms) & (rms > 0), rms, mv.MIN_MAG_ERROR
    )
    df_fit = df_fit[df_fit["filt"].isin(FILT_ORDER)]
    kept = []
    n_rejected = 0
    for filt_name, part in df_fit.groupby("filt", sort=False):
        good = magnitude_inliers(part["H"].to_numpy())
        n_rejected += int((~good).sum())
        kept.append(part.loc[good])
    if kept:
        df_fit = pd.concat(kept, ignore_index=True)
    if n_rejected:
        print(f"  multiband 2-peak compare: rejected {n_rejected} magnitude outlier(s)")
    present = [f for f in FILT_ORDER if (df_fit["filt"] == f).any()]
    if len(present) < 2:
        print(
            f"  skip: only {len(present)} filter(s) in {present} for multiband compare"
        )
        return None
    ref = REF_FILT
    if ref not in present:
        counts = {f: int((df_fit["filt"] == f).sum()) for f in present}
        ref = max(counts, key=counts.get)

    # --> Build the joint-fit arrays and choose the reference filter.
    ts = df_fit["t_days"].to_numpy(float)
    ys = df_fit["H"].to_numpy(float)
    ss = np.maximum(df_fit["rmsmag"].to_numpy(float), 1e-4)
    filts = df_fit["filt"].to_numpy()
    mask = np.ones(len(ys), bool)
    t0 = float(np.min(ts))
    K = HARMONIC_ORDER

    # --> Merge vote labels and tag campaign-duration aliases.
    vote_items = [
        (float(cand["freq"]), label)
        for cand in candidates
        for label in cand["labels"]
    ]
    merged = tag_campaign_aliases(mv.merge_frequencies(vote_items), ts)
    n_alias = sum(1 for cand in merged if is_baseline_alias(cand["labels"]))
    print(
        f"  multiband 2-peak compare: {len(merged)} voted periods, "
        f"{n_alias} land on a campaign alias"
    )

    # --> Score each voted period with a joint 2-peak harmonic fit.
    ls_power = LombScargleMultiband(ts, ys, filts, ss)
    scored = []
    for cand in merged:
        freq = float(cand["freq"])
        f_rot = freq / 2.0
        chi2, kpar = chi2_fixed(
            ts, ys, ss, filts, present, ref, f_rot, K, mask, 0.0
        )
        X, names = design_matrix(ts, filts, f_rot, K, present, ref)
        beta, _ = wls_fit(X, ys, ss, mask)
        phase_grid = np.linspace(0, 1, 400, endpoint=False)
        Xg, _ = design_matrix(
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
            design_matrix(
                ts, np.full(len(ts), ref), f_rot, K, present, ref
            )[0]
            @ beta
        )
        amp = float(np.ptp(y_model))
        scatter = scatter_about_curve(y_on_ref, y_curve, ss)
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
                "n_peaks": brightness_peaks(y_model),
                "power": power,
                "amp": amp,
                "scatter": scatter,
                "amp_ratio": ratio,
                "span_alias": int(is_baseline_alias(cand["labels"])),
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
        )
    # --> Keep the lowest-chi2 two-peak non-alias model as selected.
    scored.sort(key=lambda item: (item["chi2_nu"], item["p_rot_h"]))
    two_peaks = [item for item in scored if item["n_peaks"] == 2]
    rotations = [item for item in two_peaks if not item["span_alias"]]
    best = rotations[0] if rotations else None

    for item in scored:
        mark = " selected" if item is best else ""
        alias_mark = ", baseline alias" if item["span_alias"] else ""
        print(
            f"  {item['p_rot_h']:.4f} h ({', '.join(item['labels'])}): "
            f"chi2_nu = {item['chi2_nu']:.2f}, peaks = {item['n_peaks']}, "
            f"power = {item['power']:.3f}, A/s = {item['amp_ratio']:.3f}"
            f"{alias_mark}{mark}"
        )

    if best is None:
        if two_peaks:
            print(
                "  skip: the two-peak models are campaign aliases; "
                f"{provid} is left out"
            )
        else:
            print(f"  skip: no period with 2 peaks; {provid} is left out")
        return None

    alias_remnant = any(
        item["span_alias"] and item["chi2_nu"] < best["chi2_nu"]
        for item in two_peaks
    )
    print(
        f"  multiband 2-peak compare: P_rot = {best['p_rot_h']:.4f} h, "
        f"chi2_nu = {best['chi2_nu']:.2f}, power = {best['power']:.3f}, "
        f"A/s = {best['amp_ratio']:.3f}, peaks = 2 "
        f"(best of {len(rotations)} rotations, {len(scored)} scored)"
        + (", alias remnant" if alias_remnant else "")
    )

    def voted(item):
        return any(
            not str(label).startswith(("span", "2-span")) for label in item["labels"]
        )

    others = [item for item in two_peaks if item is not best and voted(item)]
    # --> Draw selected panel plus the other two-peak voted solutions.
    if qplot:
        mvroc_folding_plot(
            provid,
            ts,
            ss,
            filts,
            present,
            ref,
            best,
            others,
            n_obs_total,
            alias_remnant,
        )

    # --> One CSV row per plotted solution (selected + other two-peak panels).
    rows = []
    for item in [best, *others]:
        colors = compare_colors(df_obj, item, present, ref)
        amp_data, _lo, _hi = data_amplitude(item["y_on_ref"])
        row = {
            "provid": provid,
            "selected": int(item is best),
            "p_rot_fold_h": float(item["p_rot_h"]),
            "chi2_nu": float(item["chi2_nu"]),
            "n_peaks": int(item["n_peaks"]),
            "power": float(item["power"]),
            "amp_ratio": float(item["amp_ratio"]),
            "amp_data": float(amp_data),
            "span_alias": int(item["span_alias"]),
        }
        for col in COLOR_COLUMNS:
            row[col] = colors.get(col, np.nan)
        if item is best and colors:
            shown = ", ".join(
                f"{key} = {colors[key]:+.3f} +/- {colors[key + '_err']:.3f}"
                for key in (f"{a}-{b}" for a, b in color_pairs())
                if key in colors
            )
            print(f"  colors of the adopted period: {shown}")
        rows.append(row)
    pd.DataFrame(rows, columns=RESULTS_COLUMNS).to_csv(
        RESULTS_CSV, mode="a", header=False, index=False
    )
    return rows

# ---------------------------------------------------------------------------
# Foldings driver
# ---------------------------------------------------------------------------
def mvroc_foldings():
    """Run multiband two-peak folding from the command line.

    Must be launched from the pipeline/ directory.

    Stage: foldings
    Upstream: periodograms
    Downstream: shapefit

    Command-line
    ------------
    Exactly one positional argument:

    * provid - one asteroid name, e.g. "1981 US22" or 1981_US22
    * all - process every asteroid listed in periodograms.csv

    Examples
    --------
    python foldings.py "1981 US22"
    python foldings.py 1981_US22
    python foldings.py all

    Inputs
    ------
    outputs/mois/mois.parquet, outputs/states/states.parquet,
    outputs/crossmatches/crossmatches.csv
        MOI photometry with H (via mv.load_moi_photometry).
    outputs/periodograms/periodograms.csv
        Democratic vote periods (periods_dem_h) from upstream.

    Outputs
    -------
    outputs/foldings/foldings.csv
        One row per plotted two-peak solution, with color columns.
    outputs/foldings/<provid>_folding.webp
        Folded comparison figure (selected panel plus other two-peak votes).
    outputs/foldings/foldings.log
        Same text as the terminal.

    Returns
    -------
    None
    """
    # --> Check command-line arguments.
    if len(sys.argv) != 2:
        print('Usage: python foldings.py "<provid>" | python foldings.py all')
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
            f"Missing upstream periodograms CSV: {UPSTREAM_CSV}\n"
            "Run python periodograms.py first."
        )
    # --> Initialize the results CSV if it doesn't exist.
    if not os.path.isfile(RESULTS_CSV):
        pd.DataFrame(columns=RESULTS_COLUMNS).to_csv(RESULTS_CSV, index=False)
        print(f"Created {RESULTS_CSV}")
    print(f"Output: {OUT_DIR}")
    print(f"Log: {LOG_FILE}")
    print(f"Upstream: {UPSTREAM_CSV}")

    # --> Asteroids to fold = those already in periodograms.csv.
    upstream = pd.read_csv(UPSTREAM_CSV)
    if "provid" not in upstream.columns or "periods_dem_h" not in upstream.columns:
        raise SystemExit(
            f"{UPSTREAM_CSV} must have columns provid and periods_dem_h"
        )
    # --> Keep the last vote row per provid if the CSV was appended more than once.
    upstream["provid"] = upstream["provid"].astype(str).str.strip()
    upstream = upstream.drop_duplicates("provid", keep="last").reset_index(drop=True)

    # --> If not processing all asteroids, filter to the target.
    if target.lower() != "all":
        upstream = upstream[upstream["provid"] == target].reset_index(drop=True)
        if upstream.empty:
            raise SystemExit(
                f"Provid not found in {UPSTREAM_CSV}: {target}"
            )

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
            print(f"{provid}: already in foldings.csv, skip")
            skipped.append(provid)
            continue
        # --> Fold the democratic vote candidates for this asteroid.
        candidates = candidates_from_periods_dem_h(row_up["periods_dem_h"])
        rows = mvroc_folding_asteroid(df, provid, candidates)
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
    mvroc_foldings()
