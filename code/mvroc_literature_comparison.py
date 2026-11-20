#!/usr/bin/env python
"""Compare Greenstreet and SBDB rotation periods with the obs set and the MOI.

Reports how many literature rotators appear in the Rubin observation set, and
how many of those survive in figures/figures_MOI (MVROC objects of interest).
Matching expands each provid with alt_designation (Greenstreet) and the
two-provid pairs in rubin_objects_two_provids.csv.

Also writes a period–period scatter of the MOI folding solutions
(figures_MOI/solutions.csv) against the published value for every MOI
asteroid that has both. One literature period can map to 1–3 (or more)
MVROC solutions; those solutions are drawn at the same literature abscissa
and joined by a thin vertical line. Marker colour encodes the absolute
magnitude H_r.

Run from code/:

    python mvroc_literature_comparison.py
    python mvroc_literature_comparison.py --n-min 50

Writes figures/figures_MOI/literature_comparison_stats.csv with the summary
counts used in the photometry manuscript (Table tab:literature),
code/mvroc_literature_comparison.log with the full text report, and
test_dev/paper_workshop_LSST/figs/mvroc_vs_literature_periods.png for the paper.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import catalogue_pipeline as cp  # noqa: E402

DATA_DIR = ROOT / "data"
MOI_DIR = ROOT / "figures" / "figures_MOI"
MOI_CSV = MOI_DIR / "asteroid_selection.csv"
MOI_SOLUTIONS = MOI_DIR / "solutions.csv"
TWOPROV_CSV = DATA_DIR / "rubin_objects_two_provids.csv"
OUT_CSV = MOI_DIR / "literature_comparison_stats.csv"
OUT_LOG = ROOT / "mvroc_literature_comparison.log"
PAPER_FIGS = ROOT / "test_dev" / "paper_workshop_LSST" / "figs"
PERIOD_FIG = PAPER_FIGS / "mvroc_vs_literature_periods.png"
PERIOD_CSV = MOI_DIR / "mvroc_vs_literature_periods.csv"
DEFAULT_N_MIN = 100
DPI = 150


class _Tee:
    """Write the same text to the console and to a log file."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for stream in self.streams:
            stream.write(data)
            stream.flush()
        return len(data)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def _load_alias_map(path: Path) -> dict[str, set[str]]:
    alias_map: dict[str, set[str]] = {}
    if not path.exists():
        return alias_map
    raw = pd.read_csv(path)
    for _, row in raw.iterrows():
        provids = [p.strip() for p in str(row["provids"]).split(",") if p.strip()]
        for p in provids:
            alias_map.setdefault(p, set()).update(provids)
    return alias_map


def _expand(name: str, alias_map: dict[str, set[str]]) -> set[str]:
    return {name} | set(alias_map.get(name, ()))


def _paper_identity_sets(paper: pd.DataFrame, alias_map: dict[str, set[str]]):
    """Canonical Greenstreet provid -> all names that identify the same object."""
    paper_ids: dict[str, set[str]] = {}
    for _, row in paper.iterrows():
        provid = str(row["provid"]).strip()
        names = {provid}
        alt = row.get("alt_designation")
        if pd.notna(alt) and str(alt).strip():
            names.add(str(alt).strip())
        all_names: set[str] = set()
        for name in names:
            all_names |= _expand(name, alias_map)
        paper_ids[provid] = all_names
    return paper_ids


def _in_obs(names: set[str], obs_counts: pd.Series, n_min: int | None):
    hits = [n for n in names if n in obs_counts.index]
    if not hits:
        return False, 0
    n_obs = int(sum(int(obs_counts[n]) for n in hits))
    if n_min is not None and n_obs < n_min:
        return False, n_obs
    return True, n_obs


def _in_moi(names: set[str], moi_expanded: set[str], moi_provids: set[str]) -> bool:
    return any((n in moi_expanded) or (n in moi_provids) for n in names)


def compute_stats(
    n_min: int = DEFAULT_N_MIN,
    moi_csv: Path = MOI_CSV,
    solutions_csv: Path = MOI_SOLUTIONS,
):
    """Return (summary DataFrame, detail dict) for the literature cross-match."""
    moi_list = solutions_csv if solutions_csv.exists() else moi_csv
    if not moi_list.exists():
        raise SystemExit(f"Missing MOI list: {moi_list}")

    paper = cp.load_rotation_catalog()
    sbdb = cp.load_sbdb_rotation_periods()
    lit = cp.load_literature_rotation()
    alias_map = _load_alias_map(TWOPROV_CSV)
    paper_ids = _paper_identity_sets(paper, alias_map)

    print("Loading observation set...")
    df_set = cp.load_df_set()
    obs_counts = cp.observation_counts(df_set).set_index("provid")["n_obs"]
    obs_counts.index = obs_counts.index.astype(str)

    moi = pd.read_csv(moi_list)
    moi_provids = set(moi["provid"].astype(str).str.strip())
    print(f"MOI list: {moi_list} ({len(moi_provids)} asteroids)")
    moi_expanded: set[str] = set()
    for p in moi_provids:
        moi_expanded |= _expand(p, alias_map)

    # Greenstreet
    gs_in_obs = []
    gs_in_obs_nmin = []
    gs_in_obs_moi = []
    gs_in_obs_nmin_moi = []
    gs_in_obs_not_moi = []
    for provid, names in paper_ids.items():
        ok, n_obs = _in_obs(names, obs_counts, n_min=None)
        ok_n, n_obs_n = _in_obs(names, obs_counts, n_min=n_min)
        moi_hit = _in_moi(names, moi_expanded, moi_provids)
        if ok:
            gs_in_obs.append(provid)
            if moi_hit:
                gs_in_obs_moi.append(provid)
            else:
                gs_in_obs_not_moi.append((provid, n_obs, sorted(names)))
        if ok_n:
            gs_in_obs_nmin.append(provid)
            if moi_hit:
                gs_in_obs_nmin_moi.append(provid)

    # SBDB
    gs_all_names: set[str] = set()
    for names in paper_ids.values():
        gs_all_names |= names

    sb_in_obs = []
    sb_in_obs_nmin = []
    sb_in_obs_moi = []
    sb_in_obs_nmin_moi = []
    sb_also_gs = []
    sb_also_gs_nmin = []
    for provid in sbdb["provid"].astype(str).str.strip():
        names = _expand(provid, alias_map)
        ok, _ = _in_obs(names, obs_counts, n_min=None)
        ok_n, _ = _in_obs(names, obs_counts, n_min=n_min)
        moi_hit = _in_moi(names, moi_expanded, moi_provids)
        is_gs = bool(names & gs_all_names)
        if ok:
            sb_in_obs.append(provid)
            if moi_hit:
                sb_in_obs_moi.append(provid)
            if is_gs:
                sb_also_gs.append(provid)
        if ok_n:
            sb_in_obs_nmin.append(provid)
            if moi_hit:
                sb_in_obs_nmin_moi.append(provid)
            if is_gs:
                sb_also_gs_nmin.append(provid)

    # Pipeline-strict (exact provid string)
    def strict_ok(provid: str, floor: int | None) -> bool:
        if provid not in obs_counts.index:
            return False
        if floor is not None and int(obs_counts[provid]) < floor:
            return False
        return True

    gs_strict = [p for p in paper["provid"].astype(str) if strict_ok(p, None)]
    gs_strict_nmin = [p for p in paper["provid"].astype(str) if strict_ok(p, n_min)]
    sb_strict = [p for p in sbdb["provid"].astype(str) if strict_ok(p, None)]
    sb_strict_nmin = [p for p in sbdb["provid"].astype(str) if strict_ok(p, n_min)]
    lit_nmin = [
        str(r.provid).strip()
        for _, r in lit.iterrows()
        if strict_ok(str(r.provid).strip(), n_min)
    ]

    # MOI composition
    sbdb_expanded: set[str] = set()
    for p in sbdb["provid"].astype(str).str.strip():
        sbdb_expanded |= _expand(p, alias_map)

    moi_gs = [p for p in moi_provids if _expand(p, alias_map) & gs_all_names]
    moi_sb = [p for p in moi_provids if _expand(p, alias_map) & sbdb_expanded]
    moi_lit = sorted(set(moi_gs) | set(moi_sb))

    summary_rows = [
        ("greenstreet_catalogue", len(paper)),
        ("greenstreet_in_obs", len(gs_in_obs)),
        ("greenstreet_in_obs_and_moi", len(gs_in_obs_moi)),
        (f"greenstreet_in_obs_nobs_ge_{n_min}", len(gs_in_obs_nmin)),
        (f"greenstreet_in_obs_nobs_ge_{n_min}_and_moi", len(gs_in_obs_nmin_moi)),
        ("sbdb_rot_per", len(sbdb)),
        ("sbdb_rot_per_in_obs", len(sb_in_obs)),
        ("sbdb_rot_per_in_obs_and_moi", len(sb_in_obs_moi)),
        (f"sbdb_rot_per_in_obs_nobs_ge_{n_min}", len(sb_in_obs_nmin)),
        (f"sbdb_rot_per_in_obs_nobs_ge_{n_min}_and_moi", len(sb_in_obs_nmin_moi)),
        ("literature_union", len(lit)),
        (f"literature_union_in_obs_nobs_ge_{n_min}_strict", len(lit_nmin)),
        ("moi_asteroids", len(moi_provids)),
        ("moi_with_greenstreet", len(moi_gs)),
        ("moi_with_sbdb_rot_per", len(moi_sb)),
        ("moi_with_any_literature", len(moi_lit)),
        ("moi_without_literature", len(moi_provids) - len(moi_lit)),
        ("obs_unique_provids", int(obs_counts.index.nunique())),
        ("greenstreet_sbdb_overlap_in_obs", len(sb_also_gs)),
        (f"greenstreet_sbdb_overlap_in_obs_nobs_ge_{n_min}", len(sb_also_gs_nmin)),
        ("greenstreet_provid_in_obs_strict", len(gs_strict)),
        (f"greenstreet_provid_in_obs_nobs_ge_{n_min}_strict", len(gs_strict_nmin)),
        ("sbdb_provid_in_obs_strict", len(sb_strict)),
        (f"sbdb_provid_in_obs_nobs_ge_{n_min}_strict", len(sb_strict_nmin)),
    ]
    summary = pd.DataFrame(summary_rows, columns=["quantity", "N"])

    detail = {
        "n_min": n_min,
        "greenstreet_in_obs_not_moi": gs_in_obs_not_moi,
        "moi_greenstreet": sorted(moi_gs),
        "moi_sbdb": sorted(moi_sb),
        "moi_literature": moi_lit,
    }
    return summary, detail


def _positive(value) -> float | None:
    number = pd.to_numeric(value, errors="coerce")
    if pd.isna(number) or float(number) <= 0:
        return None
    return float(number)


def _literature_match_for_names(
    names: set[str], lit_entries: list[dict]
) -> tuple[float, str, str] | None:
    """Preferred published period for an identity set, or None."""
    matches = []
    for entry in lit_entries:
        if not (names & entry["names"]):
            continue
        if entry["period_lsm_hr"] is not None:
            matches.append(
                (entry["period_lsm_hr"], "Greenstreet LSM", "Greenstreet")
            )
        elif entry["period_fourier_hr"] is not None:
            matches.append(
                (
                    entry["period_fourier_hr"],
                    "Greenstreet Fourier",
                    "Greenstreet",
                )
            )
        elif entry["period_sbdb_hr"] is not None:
            matches.append((entry["period_sbdb_hr"], "SBDB", "SBDB"))
    if not matches:
        return None
    matches.sort(key=lambda item: 0 if item[2] == "Greenstreet" else 1)
    return matches[0]


def build_period_pairs(
    moi_csv: Path = MOI_CSV,
    solutions_csv: Path = MOI_SOLUTIONS,
) -> pd.DataFrame:
    """MOI folding solutions paired with Greenstreet / SBDB published periods.

    Reads ``figures_MOI/solutions.csv`` (one row per MOI asteroid, with 1–N
    slash-separated periods). Each solution becomes its own pair against the
    same published period. Literature preference is Greenstreet Lomb–Scargle,
    then Greenstreet Fourier, then SBDB ``rot_per``. Matching is alias-aware,
    as in ``compute_stats``.
    """
    if not solutions_csv.exists():
        raise SystemExit(f"Missing MOI folding solutions: {solutions_csv}")

    solutions = cp.load_moi_solutions(solutions_csv)
    alias_map = _load_alias_map(TWOPROV_CSV)
    paper = cp.load_rotation_catalog()
    paper["provid"] = paper["provid"].astype(str).str.strip()
    sbdb = cp.load_sbdb_rotation_periods()
    sbdb["provid"] = sbdb["provid"].astype(str).str.strip()

    source_by_provid: dict[str, str] = {}
    if moi_csv.exists():
        moi_meta = pd.read_csv(moi_csv)
        moi_meta["provid"] = moi_meta["provid"].astype(str).str.strip()
        if "source" in moi_meta.columns:
            source_by_provid = dict(
                zip(moi_meta["provid"], moi_meta["source"].astype(str))
            )

    lit_entries: list[dict] = []
    for _, row in paper.iterrows():
        names = {str(row["provid"]).strip()}
        alt = row.get("alt_designation")
        if pd.notna(alt) and str(alt).strip():
            names.add(str(alt).strip())
        all_names: set[str] = set()
        for name in names:
            all_names |= _expand(name, alias_map)
        lit_entries.append(
            {
                "names": all_names,
                "period_lsm_hr": _positive(row.get("period_lsm_hr")),
                "period_fourier_hr": _positive(row.get("period_fourier_hr")),
                "period_sbdb_hr": None,
                "catalog": "Greenstreet",
            }
        )
    for _, row in sbdb.iterrows():
        lit_entries.append(
            {
                "names": _expand(str(row["provid"]).strip(), alias_map),
                "period_lsm_hr": None,
                "period_fourier_hr": None,
                "period_sbdb_hr": _positive(row.get("period_sbdb_hr")),
                "catalog": "SBDB",
            }
        )

    h_by_provid: dict[str, float] = {}
    colors_path = DATA_DIR / "colors.parquet"
    if colors_path.exists():
        colors = pd.read_parquet(colors_path, columns=["provid", "H_r", "H_g"])
        colors["provid"] = colors["provid"].astype(str).str.strip()
        for _, row in colors.iterrows():
            h = row.get("H_r")
            if pd.isna(h):
                h = row.get("H_g")
            if pd.notna(h):
                h_by_provid[str(row["provid"]).strip()] = float(h)
    # Fallbacks: Greenstreet H_mag, then SBDB H.
    if "H_mag" in paper.columns:
        for _, row in paper.iterrows():
            provid = str(row["provid"]).strip()
            if provid in h_by_provid:
                continue
            h = _positive(row.get("H_mag"))
            if h is not None:
                h_by_provid[provid] = h
                alt = row.get("alt_designation")
                if pd.notna(alt) and str(alt).strip():
                    h_by_provid.setdefault(str(alt).strip(), h)
    sbdb_phys = DATA_DIR / "sbdb_physical_properties.csv"
    if sbdb_phys.exists():
        raw = pd.read_csv(sbdb_phys, usecols=["provid", "parameter", "value"], low_memory=False)
        raw = raw[raw["parameter"].astype(str).str.strip() == "H"]
        for _, row in raw.iterrows():
            provid = str(row["provid"]).strip()
            if provid in h_by_provid:
                continue
            h = pd.to_numeric(row.get("value"), errors="coerce")
            if pd.notna(h):
                h_by_provid[provid] = float(h)

    rows = []
    for _, sol_row in solutions.iterrows():
        provid = str(sol_row["provid"]).strip()
        names = _expand(provid, alias_map)
        matched = _literature_match_for_names(names, lit_entries)
        if matched is None:
            continue
        p_lit, lit_label, catalog = matched
        h_val = h_by_provid.get(provid)
        if h_val is None:
            for name in names:
                if name in h_by_provid:
                    h_val = h_by_provid[name]
                    break
        sol_periods = list(cp._moi_solution_periods(sol_row))
        n_sol = len(sol_periods)
        for sol_id, p_mvroc in sol_periods:
            rows.append(
                {
                    "provid": provid,
                    "solution": str(sol_id),
                    "n_solutions": n_sol,
                    "p_mvroc_h": float(p_mvroc),
                    "p_literature_h": float(p_lit),
                    "literature_label": lit_label,
                    "literature_catalog": catalog,
                    "mvroc_kind": f"solution_{sol_id}",
                    "source": source_by_provid.get(provid, ""),
                    "H_r": h_val if h_val is not None else np.nan,
                }
            )
    return pd.DataFrame(rows)


def plot_mvroc_vs_literature_periods(
    pairs: pd.DataFrame,
    out_path: Path = PERIOD_FIG,
) -> Path | None:
    """Scatter MOI folding solutions against published periods (log–log).

    Two side-by-side panels: Greenstreet (left) and SBDB (right). Multiple
    MVROC solutions for the same asteroid share the literature abscissa and
    are joined by a thin vertical line. Marker colour encodes H_r.
    """
    if pairs is None or pairs.empty:
        print("  no MVROC–literature period pairs to plot")
        return None

    plot = pairs.copy()
    if "H_r" not in plot.columns:
        plot["H_r"] = np.nan
    plot["H_r"] = pd.to_numeric(plot["H_r"], errors="coerce")

    x_all = plot["p_literature_h"].to_numpy(dtype=float)
    y_all = plot["p_mvroc_h"].to_numpy(dtype=float)
    h_all = plot["H_r"].to_numpy(dtype=float)
    n_asteroids = int(plot["provid"].nunique())
    lo = float(min(x_all.min(), y_all.min()) * 0.7)
    hi = float(max(x_all.max(), y_all.max()) * 1.4)
    lo = max(lo, 0.1)
    line = np.array([lo, hi])

    panels = (
        ("Greenstreet", "o", "Greenstreet et al. (2026)"),
        ("SBDB", "s", "SBDB"),
    )
    h_finite = h_all[np.isfinite(h_all)]
    if len(h_finite) == 0:
        h_vmin, h_vmax = 14.0, 20.0
    else:
        h_vmin = float(np.nanpercentile(h_finite, 5))
        h_vmax = float(np.nanpercentile(h_finite, 95))
        if h_vmax <= h_vmin:
            h_vmax = h_vmin + 1.0
    cmap = plt.cm.viridis_r
    marker_size = 42.0
    fig = plt.figure(figsize=(12.4, 6.4))
    axes = [
        fig.add_subplot(1, 2, 1),
        fig.add_subplot(1, 2, 2, sharex=None, sharey=None),
    ]
    axes[1].sharex(axes[0])
    axes[1].sharey(axes[0])
    scatters = []

    for iax, (ax, (catalog, marker, title)) in enumerate(zip(axes, panels)):
        part = plot[plot["literature_catalog"].astype(str) == catalog]
        x = part["p_literature_h"].to_numpy(dtype=float)
        y = part["p_mvroc_h"].to_numpy(dtype=float)
        h = part["H_r"].to_numpy(dtype=float)
        n_obj = int(part["provid"].nunique()) if not part.empty else 0
        n_sol = len(part)

        for _, group in part.groupby("provid", sort=False):
            if len(group) < 2:
                continue
            g = group.sort_values("p_mvroc_h")
            ax.plot(
                g["p_literature_h"].to_numpy(dtype=float),
                g["p_mvroc_h"].to_numpy(dtype=float),
                color="0.55",
                lw=0.8,
                alpha=0.45,
                zorder=2.5,
                solid_capstyle="round",
            )

        has_h = np.isfinite(h)
        if np.any(has_h):
            sc = ax.scatter(
                x[has_h],
                y[has_h],
                s=marker_size,
                c=h[has_h],
                cmap=cmap,
                vmin=h_vmin,
                vmax=h_vmax,
                marker=marker,
                zorder=3,
                edgecolors="0.15",
                linewidths=0.45,
                alpha=0.9,
            )
            scatters.append(sc)
        if np.any(~has_h):
            ax.scatter(
                x[~has_h],
                y[~has_h],
                s=marker_size,
                color="0.65",
                marker=marker,
                zorder=2.8,
                edgecolors="0.15",
                linewidths=0.45,
                alpha=0.85,
            )

        ax.plot(line, line, color="0.45", lw=1.2, label="1:1", zorder=2)
        ax.plot(line, 2.0 * line, color="C1", lw=1, ls="--", label="2:1", zorder=2)
        ax.plot(line, 0.5 * line, color="C2", lw=1, ls="--", label="1:2", zorder=2)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel(r"Literature $P_{\mathrm{rot}}$ [h]")
        ax.grid(True, which="both", alpha=0.3)
        ax.set_box_aspect(1)

        if iax == 0:
            ax.legend(loc="upper left", frameon=True, fontsize=7.5)

        ax.text(
            0.5,
            1.03,
            f"{title} ({n_obj} objects, {n_sol} solutions)",
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=10,
            clip_on=False,
        )

    axes[0].set_ylabel(r"MVROC folding solution $P_{\mathrm{rot}}$ [h]")
    plt.setp(axes[1].get_yticklabels(), visible=False)
    fig.subplots_adjust(left=0.08, right=0.90, top=0.90, bottom=0.12, wspace=0.06)

    if scatters and np.any(np.isfinite(h_all)):
        cax = fig.add_axes([0.92, 0.18, 0.015, 0.62])
        cbar = fig.colorbar(scatters[0], cax=cax)
        cbar.set_label(r"$H_r$ [mag]")

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)
    n_h = int(np.isfinite(h_all).sum())
    print(
        f"  wrote {out_path} ({n_asteroids} asteroids, {len(plot)} solutions, "
        f"{n_h} with H_r)"
    )
    return out_path


def print_report(summary: pd.DataFrame, detail: dict) -> None:
    n_min = detail["n_min"]
    get = dict(zip(summary["quantity"], summary["N"]))

    def pct(num, den):
        return 100.0 * num / den if den else float("nan")

    print()
    print("=" * 70)
    print("CATALOG SIZES")
    print("=" * 70)
    print(f"Greenstreet (asteroids-with-rotation.csv): {get['greenstreet_catalogue']}")
    print(f"SBDB with rot_per:                         {get['sbdb_rot_per']}")
    print(f"Literature union (paper ∪ SBDB):           {get['literature_union']}")
    print(f"Obs-set unique provids:                    {get['obs_unique_provids']}")
    print(f"MOI selection:                             {get['moi_asteroids']}")

    print()
    print("=" * 70)
    print("IN OBSERVATION SET (alias-aware)")
    print("=" * 70)
    g_obs, g_moi = get["greenstreet_in_obs"], get["greenstreet_in_obs_and_moi"]
    s_obs, s_moi = get["sbdb_rot_per_in_obs"], get["sbdb_rot_per_in_obs_and_moi"]
    print(
        f"Greenstreet in obs set:     {g_obs:4d} / {get['greenstreet_catalogue']}"
    )
    print(
        f"  → of which in MOI:        {g_moi:4d} / {g_obs}  "
        f"({pct(g_moi, g_obs):.1f}%)"
    )
    print(f"SBDB rot_per in obs set:    {s_obs:4d} / {get['sbdb_rot_per']}")
    print(
        f"  → of which in MOI:        {s_moi:4d} / {s_obs}  "
        f"({pct(s_moi, s_obs):.1f}%)"
    )
    print(
        f"  (SBDB∩Greenstreet in obs: {get['greenstreet_sbdb_overlap_in_obs']}; "
        f"SBDB-only: {s_obs - get['greenstreet_sbdb_overlap_in_obs']})"
    )

    print()
    print("=" * 70)
    print(f"WITH n_obs >= {n_min} (select_provid_literature threshold)")
    print("=" * 70)
    g_obs_n = get[f"greenstreet_in_obs_nobs_ge_{n_min}"]
    g_moi_n = get[f"greenstreet_in_obs_nobs_ge_{n_min}_and_moi"]
    s_obs_n = get[f"sbdb_rot_per_in_obs_nobs_ge_{n_min}"]
    s_moi_n = get[f"sbdb_rot_per_in_obs_nobs_ge_{n_min}_and_moi"]
    print(f"Greenstreet in obs (n≥{n_min}):   {g_obs_n:4d} / {get['greenstreet_catalogue']}")
    print(
        f"  → of which in MOI:           {g_moi_n:4d} / {g_obs_n}  "
        f"({pct(g_moi_n, g_obs_n):.1f}%)"
    )
    print(f"SBDB rot_per in obs (n≥{n_min}):  {s_obs_n:4d} / {get['sbdb_rot_per']}")
    print(
        f"  → of which in MOI:           {s_moi_n:4d} / {s_obs_n}  "
        f"({pct(s_moi_n, s_obs_n):.1f}%)"
    )

    print()
    print("=" * 70)
    print("PIPELINE-STRICT (exact provid only)")
    print("=" * 70)
    print(
        f"Greenstreet provid in obs:        "
        f"{get['greenstreet_provid_in_obs_strict']} / {get['greenstreet_catalogue']}"
    )
    print(
        f"Greenstreet provid n≥{n_min}:        "
        f"{get[f'greenstreet_provid_in_obs_nobs_ge_{n_min}_strict']} / "
        f"{get['greenstreet_catalogue']}"
    )
    print(
        f"SBDB provid in obs:               "
        f"{get['sbdb_provid_in_obs_strict']} / {get['sbdb_rot_per']}"
    )
    print(
        f"SBDB provid n≥{n_min}:               "
        f"{get[f'sbdb_provid_in_obs_nobs_ge_{n_min}_strict']} / {get['sbdb_rot_per']}"
    )
    print(
        f"Literature union n≥{n_min}:          "
        f"{get[f'literature_union_in_obs_nobs_ge_{n_min}_strict']} / "
        f"{get['literature_union']}"
    )

    print()
    print("=" * 70)
    print("MOI COMPOSITION vs LITERATURE")
    print("=" * 70)
    n_moi = get["moi_asteroids"]
    print(
        f"MOI ∩ Greenstreet:   {get['moi_with_greenstreet']:4d} / {n_moi}  "
        f"({pct(get['moi_with_greenstreet'], n_moi):.1f}%)"
    )
    print(
        f"MOI ∩ SBDB rot_per:  {get['moi_with_sbdb_rot_per']:4d} / {n_moi}  "
        f"({pct(get['moi_with_sbdb_rot_per'], n_moi):.1f}%)"
    )
    print(
        f"MOI ∩ any literature:{get['moi_with_any_literature']:4d} / {n_moi}  "
        f"({pct(get['moi_with_any_literature'], n_moi):.1f}%)"
    )
    print(
        f"MOI without lit. P:  {get['moi_without_literature']:4d} / {n_moi}  "
        f"({pct(get['moi_without_literature'], n_moi):.1f}%)"
    )

    not_moi = detail["greenstreet_in_obs_not_moi"]
    print()
    print("--- Greenstreet in obs but NOT in MOI ---")
    print(f"count: {len(not_moi)}")
    for provid, n_obs, names in not_moi:
        print(f"  {provid:16s} n_obs≈{n_obs:4d}  ids={names}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--n-min",
        type=int,
        default=DEFAULT_N_MIN,
        help=f"Observation floor matching select_provid_literature (default {DEFAULT_N_MIN})",
    )
    parser.add_argument(
        "--moi-csv",
        type=Path,
        default=MOI_CSV,
        help=f"MOI selection CSV for metadata (default {MOI_CSV})",
    )
    parser.add_argument(
        "--solutions-csv",
        type=Path,
        default=MOI_SOLUTIONS,
        help=f"MOI folding solutions CSV (default {MOI_SOLUTIONS})",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=OUT_CSV,
        help=f"Summary counts CSV (default {OUT_CSV})",
    )
    parser.add_argument(
        "--out-log",
        type=Path,
        default=OUT_LOG,
        help=f"Text report log (default {OUT_LOG})",
    )
    parser.add_argument(
        "--quiet-missing",
        action="store_true",
        help="Do not list Greenstreet objects that are in obs but not in MOI",
    )
    parser.add_argument(
        "--period-fig",
        type=Path,
        default=PERIOD_FIG,
        help=f"Period comparison PNG for the paper (default {PERIOD_FIG})",
    )
    parser.add_argument(
        "--period-csv",
        type=Path,
        default=PERIOD_CSV,
        help=f"Period pairs CSV (default {PERIOD_CSV})",
    )
    parser.add_argument(
        "--skip-period-fig",
        action="store_true",
        help="Skip the MVROC vs literature period figure",
    )
    args = parser.parse_args(argv)

    args.out_log.parent.mkdir(parents=True, exist_ok=True)
    log_file = args.out_log.open("w", encoding="utf-8")
    previous_stdout = sys.stdout
    sys.stdout = _Tee(previous_stdout, log_file)
    try:
        started = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        print(f"mvroc_literature_comparison.py  started {started}")
        print(
            f"n_min={args.n_min}  moi_csv={args.moi_csv}  "
            f"solutions_csv={args.solutions_csv}"
        )
        print()

        summary, detail = compute_stats(
            n_min=args.n_min,
            moi_csv=args.moi_csv,
            solutions_csv=args.solutions_csv,
        )
        if args.quiet_missing:
            detail = dict(detail)
            detail["greenstreet_in_obs_not_moi"] = []

        print_report(summary, detail)

        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.out_csv, index=False)
        print()
        print(f"Wrote {args.out_csv}")
        print(f"Wrote {args.out_log}")

        if not args.skip_period_fig:
            print()
            print("=" * 70)
            print("MVROC vs LITERATURE PERIODS (MOI folding solutions)")
            print("=" * 70)
            pairs = build_period_pairs(
                moi_csv=args.moi_csv,
                solutions_csv=args.solutions_csv,
            )
            args.period_csv.parent.mkdir(parents=True, exist_ok=True)
            pairs.to_csv(args.period_csv, index=False)
            n_ast = int(pairs["provid"].nunique()) if not pairs.empty else 0
            print(
                f"Wrote {args.period_csv} "
                f"({n_ast} asteroids, {len(pairs)} solution pairs)"
            )
            if not pairs.empty:
                print(pairs.drop_duplicates("provid")["literature_catalog"]
                      .value_counts().to_string())
                print(
                    "solutions per asteroid:\n"
                    + pairs.groupby("provid")["solution"].count()
                    .value_counts()
                    .sort_index()
                    .to_string()
                )
            plot_mvroc_vs_literature_periods(pairs, out_path=args.period_fig)
    finally:
        sys.stdout = previous_stdout
        log_file.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
