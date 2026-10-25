#!/usr/bin/env python
"""Compare Greenstreet and SBDB rotation periods with the obs set and the MOI.

Reports how many literature rotators appear in the Rubin observation set, and
how many of those survive in figures/figures_MOI (MVROC objects of interest).
Matching expands each provid with alt_designation (Greenstreet) and the
two-provid pairs in rubin_objects_two_provids.csv.

Run from code/:

    python mvroc_literature_comparison.py
    python mvroc_literature_comparison.py --n-min 50

Writes figures/figures_MOI/literature_comparison_stats.csv with the summary
counts used in the photometry manuscript (Table tab:literature), and
code/mvroc_literature_comparison.log with the full text report.
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

import catalogue_pipeline as cp  # noqa: E402

DATA_DIR = ROOT / "data"
MOI_DIR = ROOT / "figures" / "figures_MOI"
MOI_CSV = MOI_DIR / "asteroid_selection.csv"
TWOPROV_CSV = DATA_DIR / "rubin_objects_two_provids.csv"
OUT_CSV = MOI_DIR / "literature_comparison_stats.csv"
OUT_LOG = ROOT / "mvroc_literature_comparison.log"
DEFAULT_N_MIN = 50


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


def compute_stats(n_min: int = DEFAULT_N_MIN, moi_csv: Path = MOI_CSV):
    """Return (summary DataFrame, detail dict) for the literature cross-match."""
    if not moi_csv.exists():
        raise SystemExit(f"Missing MOI selection: {moi_csv}")

    paper = cp.load_rotation_catalog()
    sbdb = cp.load_sbdb_rotation_periods()
    lit = cp.load_literature_rotation()
    alias_map = _load_alias_map(TWOPROV_CSV)
    paper_ids = _paper_identity_sets(paper, alias_map)

    print("Loading observation set...")
    df_set = cp.load_df_set()
    obs_counts = cp.observation_counts(df_set).set_index("provid")["n_obs"]
    obs_counts.index = obs_counts.index.astype(str)

    moi = pd.read_csv(moi_csv)
    moi_provids = set(moi["provid"].astype(str).str.strip())
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
        help=f"MOI selection CSV (default {MOI_CSV})",
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
    args = parser.parse_args(argv)

    args.out_log.parent.mkdir(parents=True, exist_ok=True)
    log_file = args.out_log.open("w", encoding="utf-8")
    previous_stdout = sys.stdout
    sys.stdout = _Tee(previous_stdout, log_file)
    try:
        started = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        print(f"mvroc_literature_comparison.py  started {started}")
        print(f"n_min={args.n_min}  moi_csv={args.moi_csv}")
        print()

        summary, detail = compute_stats(n_min=args.n_min, moi_csv=args.moi_csv)
        if args.quiet_missing:
            detail = dict(detail)
            detail["greenstreet_in_obs_not_moi"] = []

        print_report(summary, detail)

        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.out_csv, index=False)
        print()
        print(f"Wrote {args.out_csv}")
        print(f"Wrote {args.out_log}")
    finally:
        sys.stdout = previous_stdout
        log_file.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
