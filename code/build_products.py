#!/usr/bin/env python
"""Assemble MVROC pipeline products under products/ from pipeline/inputs/
and code/figures (figures_nobs for MOI; figures_MOI / figures_mvroc_sol for
MOC and solutions).

Run from the repository root:
    python code/build_products.py
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "pipeline" / "inputs"
DATA = ROOT / "code" / "data"
MOI = ROOT / "code" / "figures" / "figures_MOI"
FOLD = ROOT / "code" / "figures" / "figures_mvroc_sol"
NOBS = ROOT / "code" / "figures" / "figures_nobs"
OBS_SRC = INPUTS / "obs_sbn_X05_full.parquet"
OUT = ROOT / "products"

# Catalogue column order from CATALOGUE.md (plus solution / classification).
CATALOGUE_COLUMNS = [
    # Identification + solution
    "provid",
    "solution_id",
    "permid",
    "mpc_designation",
    # Observations
    "epoch_utc_ini",
    "epoch_utc_end",
    "n_obs",
    "n_obs_u",
    "n_obs_g",
    "n_obs_r",
    "n_obs_i",
    "n_obs_z",
    "n_obs_y",
    # Orbit
    "epoch_tdb",
    "a",
    "e",
    "q",
    "incl",
    "node",
    "peri",
    "M",
    # Mean absolute magnitudes
    "H_u",
    "H_g",
    "H_r",
    "H_i",
    "H_z",
    "H_y",
    # Phase-corrected colors (joint fit)
    "u-g",
    "u-g_err",
    "u-r",
    "u-r_err",
    "g-r",
    "g-r_err",
    "r-i",
    "r-i_err",
    "i-z",
    "i-z_err",
    "z-y",
    "z-y_err",
    # Photometric taxonomy (gri slope; coarse C/S/D complex)
    "gri_slope",
    "photo_spec_T",
    # Rotation: adopted solution (fold)
    "p_rot_h",
    "chi2_nu",
    # Rotation: multiband Lk amplitudes / residuals (filled when available)
    "amp_lsm_k",
    "amp_photo_k",
    "rms_multi",
    "n_outliers_multi",
    # Shape: spheres (catalogue names use _k)
    "q_sph_k",
    "q_sph_err_k",
    "phi_sph_deg_k",
    "phi_sph_err_deg_k",
    "chi_r_sph_k",
    "amp_sph_k",
    "m_mid_sph_k",
    "m_sep_k",
    "p_rot_sph_h_k",
    # Shape: ellipsoid
    "b_over_a_k",
    "b_over_a_err_k",
    "amp_ell_k",
    "m_mid_ell_k",
    "m_ab_k",
    "p_rot_ell_h_k",
    "chi_r_ell_k",
    # Classification
    "tags",
    "Q",
    "use_band",
    # SBDB
    "H_sbdb",
    "H_sbdb_err",
    "albedo_sbdb",
    "albedo_sbdb_err",
    "diameter_sbdb",
    "diameter_sbdb_err",
    "extent_sbdb",
    "extent_sbdb_err",
    "pole_sbdb",
    "pole_sbdb_err",
    "density_sbdb",
    "density_sbdb_err",
    "rot_per_sbdb",
    "rot_per_sbdb_err",
    "GM_sbdb",
    "GM_sbdb_err",
    "BV_sbdb",
    "BV_sbdb_err",
    "UB_sbdb",
    "UB_sbdb_err",
    "spec_B_sbdb",
    "spec_B_sbdb_err",
    "spec_T_sbdb",
    "spec_T_sbdb_err",
    # Greenstreet
    "rot_per_lsm_gst",
    "rot_per_fr_gst",
    # Dynamics (osculating elements; resonance is a proximity flag)
    "dyn_class",
    "T_earth",
    "T_mars",
    "T_jup",
    "T_nep",
    "flag_resonant",
    "resonant_ratio",
]

SHAPE_RENAME = {
    "q_sph": "q_sph_k",
    "q_sph_err": "q_sph_err_k",
    "phi_sph_deg": "phi_sph_deg_k",
    "phi_sph_err_deg": "phi_sph_err_deg_k",
    "chi_r_sph": "chi_r_sph_k",
    "amp_sph": "amp_sph_k",
    "m_mid_sph": "m_mid_sph_k",
    "m_sep": "m_sep_k",
    "p_rot_sph_h": "p_rot_sph_h_k",
    "b_over_a": "b_over_a_k",
    "b_over_a_err": "b_over_a_err_k",
    "amp_ell": "amp_ell_k",
    "m_mid_ell": "m_mid_ell_k",
    "m_ab": "m_ab_k",
    "p_rot_ell_h": "p_rot_ell_h_k",
    "chi_r_ell": "chi_r_ell_k",
}

SOLAR_COLOR = {"g": 0.45, "i": -0.12, "z": -0.16, "y": -0.18}

# Mean heliocentric distances (au) used in the Tisserand parameter.
A_EARTH = 1.000
A_MARS = 1.524
A_JUP = 5.204
A_NEP = 30.07

# Jupiter MMRs: (label, asteroid revs, Jupiter revs, half-width in au).
# a0 = A_JUP * (n_jup / n_ast) ** (2/3). Wide islands for 3:2, 4:3, and 1:1;
# narrow windows for the Kirkwood gaps.
_JUPITER_MMR = (
    ("4:1", 4, 1, 0.03),
    ("3:1", 3, 1, 0.03),
    ("5:2", 5, 2, 0.03),
    ("7:3", 7, 3, 0.03),
    ("2:1", 2, 1, 0.04),
    ("3:2", 3, 2, 0.15),
    ("4:3", 4, 3, 0.08),
    ("1:1", 1, 1, 0.25),
)

# gri slope (% per 100 nm) of the C, S, and D sketches in catalogue_pipeline.spectrum.
_LSST_G_NM = 483.0
_LSST_I_NM = 755.0
_GRI_SLOPE_C = 1.25
_GRI_SLOPE_S = 14.6
_GRI_SLOPE_D = 22.5


def tisserand(a: np.ndarray, e: np.ndarray, incl_deg: np.ndarray, a_planet: float) -> np.ndarray:
    """Tisserand parameter w.r.t. a planet, from osculating a (au), e, i (deg)."""
    a = np.asarray(a, dtype=float)
    e = np.asarray(e, dtype=float)
    incl = np.deg2rad(np.asarray(incl_deg, dtype=float))
    ok = np.isfinite(a) & (a > 0) & np.isfinite(e) & (e < 1) & np.isfinite(incl)
    out = np.full(a.shape, np.nan)
    cos_i = np.cos(incl[ok])
    out[ok] = a_planet / a[ok] + 2.0 * cos_i * np.sqrt((a[ok] / a_planet) * (1.0 - e[ok] ** 2))
    return out


def _mmr_centres() -> list[tuple[str, float, float]]:
    centres = []
    for label, n_ast, n_jup, half in _JUPITER_MMR:
        a0 = A_JUP * (n_jup / n_ast) ** (2.0 / 3.0)
        centres.append((label, a0, half))
    return centres


def resonant_ratio_of(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nearest Jupiter MMR inside its half-width. Returns (flag 0/1, ratio or '')."""
    a = np.asarray(a, dtype=float)
    flag = np.zeros(a.shape, dtype=int)
    ratio = np.full(a.shape, "", dtype=object)
    centres = _mmr_centres()
    score = np.full(a.shape, np.inf)
    for label, a0, half in centres:
        dist = np.abs(a - a0)
        hit = np.isfinite(a) & (dist <= half) & ((dist / half) < score)
        score[hit] = dist[hit] / half
        flag[hit] = 1
        ratio[hit] = label
    return flag, ratio


def dynamical_class(a: np.ndarray, e: np.ndarray, incl: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Population label from osculating elements. First match wins."""
    a = np.asarray(a, dtype=float)
    incl = np.asarray(incl, dtype=float)
    q = np.asarray(q, dtype=float)
    labels = np.full(a.shape, "", dtype=object)
    finite = np.isfinite(a) & np.isfinite(q) & np.isfinite(incl)
    labels[finite] = "Other"
    rules = (
        ("TNO", finite & (a >= 30.1)),
        ("Trojan", finite & (a >= 4.95) & (a <= 5.45)),
        ("Centaur", finite & (a > 5.45) & (a < 30.1)),
        ("Thule", finite & (a >= 4.20) & (a < 4.50)),
        ("Hilda", finite & (a >= 3.70) & (a < 4.20)),
        ("Cybele", finite & (a >= 3.28) & (a < 3.70)),
        ("NEA", finite & (q < 1.3)),
        ("Hungaria", finite & (a >= 1.78) & (a < 2.06) & (incl >= 16) & (q >= 1.3)),
        ("Mars-crosser", finite & (q < 1.666) & (a < 2.06)),
        ("InnerMB", finite & (a >= 2.06) & (a < 2.50)),
        ("MiddleMB", finite & (a >= 2.50) & (a < 2.82)),
        ("OuterMB", finite & (a >= 2.82) & (a < 3.28)),
    )
    assigned = np.zeros(a.shape, dtype=bool)
    for name, mask in rules:
        take = mask & ~assigned
        labels[take] = name
        assigned |= take
    return labels


def gri_slope_and_class(g_r: np.ndarray, r_i: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Solar-normalized g–i slope (%/100 nm) and nearest C/S/D template."""
    g_r = np.asarray(g_r, dtype=float)
    r_i = np.asarray(r_i, dtype=float)
    ok = np.isfinite(g_r) & np.isfinite(r_i)
    r_g = np.full(g_r.shape, np.nan)
    r_i_ref = np.full(g_r.shape, np.nan)
    # R(b)/R(r) = 10^(-0.4 * ((m_b - m_r) - solar)).
    r_g[ok] = 10.0 ** (-0.4 * (g_r[ok] - SOLAR_COLOR["g"]))
    i_minus_r = -r_i
    r_i_ref[ok] = 10.0 ** (-0.4 * (i_minus_r[ok] - SOLAR_COLOR["i"]))
    slope = np.full(g_r.shape, np.nan)
    slope[ok] = (r_i_ref[ok] - r_g[ok]) / (_LSST_I_NM - _LSST_G_NM) * 1.0e4
    spec = np.full(g_r.shape, "", dtype=object)
    # Midpoints between the C, S, and D template slopes.
    c_s = 0.5 * (_GRI_SLOPE_C + _GRI_SLOPE_S)
    s_d = 0.5 * (_GRI_SLOPE_S + _GRI_SLOPE_D)
    spec[ok & (slope < c_s)] = "C"
    spec[ok & (slope >= c_s) & (slope < s_d)] = "S"
    spec[ok & (slope >= s_d)] = "D"
    return slope, spec


def add_dynamical_columns(cat: pd.DataFrame) -> pd.DataFrame:
    cat = cat.copy()
    a = cat["a"] if "a" in cat.columns else np.nan
    e = cat["e"] if "e" in cat.columns else np.nan
    incl = cat["incl"] if "incl" in cat.columns else np.nan
    q = cat["q"] if "q" in cat.columns else np.nan
    cat["dyn_class"] = dynamical_class(a, e, incl, q)
    cat["T_earth"] = tisserand(a, e, incl, A_EARTH)
    cat["T_mars"] = tisserand(a, e, incl, A_MARS)
    cat["T_jup"] = tisserand(a, e, incl, A_JUP)
    cat["T_nep"] = tisserand(a, e, incl, A_NEP)
    flag, ratio = resonant_ratio_of(cat["a"].to_numpy(dtype=float) if "a" in cat.columns else np.array([]))
    cat["flag_resonant"] = flag
    cat["resonant_ratio"] = ratio
    return cat


def add_taxonomy_columns(cat: pd.DataFrame) -> pd.DataFrame:
    cat = cat.copy()
    g_r = cat["g-r"] if "g-r" in cat.columns else np.nan
    r_i = cat["r-i"] if "r-i" in cat.columns else np.nan
    slope, spec = gri_slope_and_class(g_r, r_i)
    cat["gri_slope"] = slope
    cat["photo_spec_T"] = spec
    return cat


def _print(msg: str) -> None:
    print(msg, flush=True)


def expand_solutions(solutions: pd.DataFrame) -> pd.DataFrame:
    """Expand solutions.csv into one row per (provid, solution_id, p_rot_h)."""
    rows = []
    for _, row in solutions.iterrows():
        ids = [p.strip() for p in str(row["solutions"]).split("/") if p.strip()]
        periods = [float(p) for p in str(row["periods_h"]).split("/") if p.strip()]
        if len(ids) != len(periods):
            raise ValueError(f"solutions/periods mismatch for {row['provid']}")
        for sid, p_rot in zip(ids, periods):
            rows.append(
                {
                    "provid": str(row["provid"]).strip(),
                    "solution_id": str(sid),
                    "p_rot_h": float(p_rot),
                }
            )
    return pd.DataFrame(rows)


def _parse_sbdb_number(value):
    """Number from an SBDB cell, dropping a trailing unit such as 'h' or 'km'."""
    number = pd.to_numeric(value, errors="coerce")
    if pd.notna(number):
        return float(number)
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return np.nan
    text = str(value).strip().lower().replace("hours", "h").replace("hour", "h")
    text = text.replace("h", " ").replace("km", " ").strip()
    if not text:
        return np.nan
    try:
        return float(text.split()[0])
    except (TypeError, ValueError):
        return np.nan


def pivot_sbdb(path: Path) -> pd.DataFrame:
    """Wide SBDB physical parameters keyed by provid."""
    raw = pd.read_csv(path, low_memory=False)
    raw["provid"] = raw["provid"].astype(str).str.strip()
    raw["parameter"] = raw["parameter"].astype(str).str.strip()
    keep = {
        "H": "H_sbdb",
        "albedo": "albedo_sbdb",
        "diameter": "diameter_sbdb",
        "extent": "extent_sbdb",
        "pole": "pole_sbdb",
        "density": "density_sbdb",
        "rot_per": "rot_per_sbdb",
        "GM": "GM_sbdb",
        "BV": "BV_sbdb",
        "UB": "UB_sbdb",
        "spec_B": "spec_B_sbdb",
        "spec_T": "spec_T_sbdb",
    }
    raw = raw[raw["parameter"].isin(keep)].copy()
    raw["col"] = raw["parameter"].map(keep)
    numeric = {
        "H_sbdb",
        "albedo_sbdb",
        "diameter_sbdb",
        "density_sbdb",
        "rot_per_sbdb",
        "GM_sbdb",
        "BV_sbdb",
        "UB_sbdb",
    }
    frames = []
    for col, part in raw.groupby("col", sort=False):
        part = part.drop_duplicates("provid", keep="first")
        if col in numeric:
            value = part["value"].map(_parse_sbdb_number)
            sigma = pd.to_numeric(part["sigma"], errors="coerce")
        else:
            value = part["value"].astype(str).where(part["value"].notna(), np.nan)
            sigma = pd.Series(np.nan, index=part.index)
        frames.append(
            pd.DataFrame(
                {
                    "provid": part["provid"].to_numpy(),
                    col: value.to_numpy(),
                    f"{col}_err": sigma.to_numpy(),
                }
            )
        )
    out = frames[0]
    for frame in frames[1:]:
        out = out.merge(frame, on="provid", how="outer")
    return out


def greenstreet_periods(path: Path) -> pd.DataFrame:
    gst = pd.read_csv(path)
    rows = []
    for _, row in gst.iterrows():
        for key in ("provid", "alt_designation", "designation"):
            name = row.get(key)
            if pd.isna(name) or not str(name).strip():
                continue
            rows.append(
                {
                    "provid": str(name).strip(),
                    "rot_per_lsm_gst": row.get("period_lsm_hr"),
                    "rot_per_fr_gst": row.get("period_fourier_hr"),
                }
            )
    out = pd.DataFrame(rows).drop_duplicates("provid", keep="first")
    return out


def reflectance_from_colors(row: pd.Series) -> dict:
    """Solar-normalized reflectance relative to r (CATALOGUE.md)."""
    out = {}
    n_map = {"g": "n_obs_g", "i": "n_obs_i", "z": "n_obs_z", "y": "n_obs_y"}
    # Prefer fit colors; fall back to mean colors if needed.
    for band, solar in SOLAR_COLOR.items():
        if band == "g":
            color = row.get("g-r")
            err = row.get("g-r_err")
        elif band == "i":
            # i−r = −(r−i)
            if pd.notna(row.get("r-i")):
                color = -float(row["r-i"])
                err = row.get("r-i_err")
            else:
                color = np.nan
                err = np.nan
        elif band == "z":
            # z−r = (i−z) + (r−i) with opposite? Catalogue: z−r from excess.
            # Prefer chain: z-r = -(r-z) if present, else (i-z) - (r-i) wait:
            # r-z = (r-i)+(i-z) => z-r = -(r-z)
            if pd.notna(row.get("r-z")):
                color = -float(row["r-z"])
                err = row.get("r-z_err")
            elif pd.notna(row.get("r-i")) and pd.notna(row.get("i-z")):
                color = -(float(row["r-i"]) + float(row["i-z"]))
                e1 = row.get("r-i_err")
                e2 = row.get("i-z_err")
                if pd.notna(e1) and pd.notna(e2):
                    err = float(np.sqrt(float(e1) ** 2 + float(e2) ** 2))
                else:
                    err = np.nan
            else:
                color = np.nan
                err = np.nan
        else:  # y
            if pd.notna(row.get("r-y")):
                color = -float(row["r-y"])
                err = row.get("r-y_err")
            elif pd.notna(row.get("z-y")) and pd.notna(row.get("r-z")):
                color = -float(row["r-z"]) + float(row["z-y"])
                err = np.nan
            else:
                color = np.nan
                err = np.nan
        out[f"{band}-r_spec"] = color
        out[f"{band}-r_spec_err"] = err
        if pd.notna(color):
            excess = float(color) - solar
            out[f"R_{band}"] = 10 ** (-0.4 * excess)
            if pd.notna(err):
                # dR/dcolor ≈ R * (-0.4 * ln(10))
                out[f"R_{band}_err"] = abs(
                    out[f"R_{band}"] * (-0.4 * np.log(10)) * float(err)
                )
            else:
                out[f"R_{band}_err"] = np.nan
        else:
            out[f"R_{band}"] = np.nan
            out[f"R_{band}_err"] = np.nan
        out[f"n_{band}_spec"] = row.get(n_map[band], np.nan)
    return out


def build_db_vro_mpc() -> pd.DataFrame:
    _print("Building db_vro_mpc.csv …")
    mpc = pd.read_csv(DATA / "mpc_designations.csv")
    mpc["provid_vro"] = mpc["provid"].astype(str).str.strip()
    mpc["provid_mpc"] = mpc["mpc_designation"].astype(str).str.strip()
    mpc = mpc[["provid_vro", "provid_mpc"]].drop_duplicates("provid_vro")

    colors = pd.read_parquet(
        DATA / "colors.parquet", columns=["provid", "mpc_designation", "permid"]
    )
    colors["provid_vro"] = colors["provid"].astype(str).str.strip()
    colors["provid_mpc"] = colors["mpc_designation"].astype(str).str.strip()
    colors = colors[["provid_vro", "provid_mpc", "permid"]].drop_duplicates(
        "provid_vro"
    )

    # Prefer mpc_designations as the pipeline cross-match; fill permid from colors.
    out = mpc.merge(
        colors[["provid_vro", "permid"]], on="provid_vro", how="left"
    )
    # Objects only in colors (should be rare).
    only_colors = colors[~colors["provid_vro"].isin(out["provid_vro"])]
    out = pd.concat([out, only_colors], ignore_index=True)

    # Match class: numbered vs provisional (heuristic for the report).
    numbered = out["provid_mpc"].str.match(r"^\(\d+\)", na=False)
    out["match_kind"] = np.where(numbered, "numbered", "provisional")
    out["n_vro_for_mpc"] = out.groupby("provid_mpc")["provid_vro"].transform("count")
    out = out.sort_values("provid_vro").reset_index(drop=True)
    out.to_csv(OUT / "db_vro_mpc.csv", index=False)
    _print(f"  {len(out):,} rows → db_vro_mpc.csv")
    return out


def build_info_vro_mpc(db_vro_mpc: pd.DataFrame) -> None:
    _print("Building info_vro_mpc.md …")
    two = pd.read_csv(DATA / "rubin_objects_two_provids.csv")
    n = len(db_vro_mpc)
    n_mpc = db_vro_mpc["provid_mpc"].nunique()
    n_multi = int((db_vro_mpc["n_vro_for_mpc"] > 1).sum())
    n_numbered = int((db_vro_mpc["match_kind"] == "numbered").sum())
    lines = [
        "# Rubin–MPC cross-match report",
        "",
        "Built from `code/data/mpc_designations.csv`, `colors.parquet`, and",
        "`rubin_objects_two_provids.csv` (pipeline step 2).",
        "",
        "## Summary",
        "",
        f"- Unique `provid_vro` with an MPC designation: **{n:,}**",
        f"- Unique `provid_mpc`: **{n_mpc:,}**",
        f"- Numbered MPC designations: **{n_numbered:,}**",
        f"- Provisional / unnumbered MPC designations: **{n - n_numbered:,}**",
        f"- `provid_vro` that share an MPC designation with another Rubin id: **{n_multi:,}**",
        f"- MPC objects known to gather several Rubin provids (`rubin_objects_two_provids.csv`): **{len(two):,}**",
        "",
        "## Many-to-one associations",
        "",
        "A single `provid_mpc` can collect several `provid_vro`. The equivalence",
        "table `db_vro_mpc.csv` records every pair; `n_vro_for_mpc` counts how",
        "many Rubin ids map to that MPC designation.",
        "",
        f"Objects with ≥2 Rubin provids (from the curated list): {len(two):.0f}.",
        "",
        "| mpc_designation | provids | n_obs_total |",
        "|---|---|---:|",
    ]
    for _, row in two.head(25).iterrows():
        lines.append(
            f"| {row['mpc_designation']} | {row['provids']} | {int(row['n_obs_total'])} |"
        )
    if len(two) > 25:
        lines.append(f"| … | ({len(two) - 25} more rows) | |")
    lines.append("")
    (OUT / "info_vro_mpc.md").write_text("\n".join(lines), encoding="utf-8")
    _print("  wrote info_vro_mpc.md")


def build_orbits(db_vro_mpc: pd.DataFrame) -> pd.DataFrame:
    """One orbital-element row per provid_mpc from the geometry table."""
    _print("Building db_mvroc_orbits.csv …")
    cols = [
        "provid",
        "mpc_designation",
        "permid",
        "obstime",
        "a",
        "e",
        "q",
        "incl",
        "node",
        "peri",
        "M",
        "dt_from_epoch_days",
        "obs_et",
    ]
    # Stream geometry and keep the first row per provid (elements are epoch-fixed).
    seen = set()
    chunks = []
    pf = pq.ParquetFile(DATA / "obs_final_elements_geometry_H.parquet")
    for batch in pf.iter_batches(batch_size=500_000, columns=cols):
        df = batch.to_pandas()
        df["provid"] = df["provid"].astype(str).str.strip()
        df = df[~df["provid"].isin(seen)]
        if df.empty:
            continue
        first = df.drop_duplicates("provid", keep="first")
        seen.update(first["provid"].tolist())
        chunks.append(first)
        _print(f"  orbits scan: {len(seen):,} unique provids")
    geo = pd.concat(chunks, ignore_index=True)
    # Epoch of elements ≈ obstime − dt_from_epoch_days (TDB days).
    obstime = pd.to_datetime(geo["obstime"], utc=True)
    epoch = obstime - pd.to_timedelta(geo["dt_from_epoch_days"], unit="D")
    orbits = pd.DataFrame(
        {
            "provid_vro": geo["provid"],
            "provid_mpc": geo["mpc_designation"].astype(str).str.strip(),
            "permid": geo["permid"],
            "epoch_tdb": epoch.dt.strftime("%Y-%m-%d %H:%M:%S"),
            "a": geo["a"],
            "e": geo["e"],
            "q": geo["q"],
            "incl": geo["incl"],
            "node": geo["node"],
            "peri": geo["peri"],
            "M": geo["M"],
        }
    )
    # One row per MPC designation (prefer the first VRO id encountered).
    orbits = orbits.drop_duplicates("provid_mpc", keep="first")
    # Restrict to cross-matched MPC ids when possible; keep geometry extras too.
    wanted = set(db_vro_mpc["provid_mpc"].astype(str))
    orbits = orbits[orbits["provid_mpc"].isin(wanted) | orbits["provid_mpc"].notna()]
    orbits = orbits.sort_values("provid_mpc").reset_index(drop=True)
    orbits.to_csv(OUT / "db_mvroc_orbits.csv", index=False)
    _print(f"  {len(orbits):,} rows → db_mvroc_orbits.csv")
    return orbits


def build_state() -> None:
    """Per-observation state table keyed by observation id."""
    _print("Building db_mvroc_state.parquet …")
    geo_cols = [
        "obstime",
        "provid",
        "mpc_designation",
        "permid",
        "obs_et",
        "dt_from_epoch_days",
        "r",
        "delta",
        "phase_angle",
        "mag",
        "H",
    ]
    # Load observation ids for rows that have a provid (needed for the join).
    _print("  loading observation ids …")
    # Parquet nulls are not filtered by != None; drop empty/NaN after load.
    obs = pq.read_table(
        DATA / "obs_sbn_X05_full_mpc_desig.parquet",
        columns=["id", "provid", "obstime"],
        filters=[("provid", "!=", "")],
    ).to_pandas()
    obs["provid"] = obs["provid"].astype(str).str.strip()
    obs = obs[obs["provid"].ne("") & obs["provid"].ne("None") & obs["provid"].ne("nan")]
    obs["obstime"] = pd.to_datetime(obs["obstime"], utc=True)
    obs = obs.dropna(subset=["provid", "obstime"])
    _print(f"  observation ids with provid: {len(obs):,}")

    writer = None
    n_written = 0
    pf = pq.ParquetFile(DATA / "obs_final_elements_geometry_H.parquet")
    for i, batch in enumerate(pf.iter_batches(batch_size=400_000, columns=geo_cols)):
        geo = batch.to_pandas()
        geo["provid"] = geo["provid"].astype(str).str.strip()
        geo["obstime"] = pd.to_datetime(geo["obstime"], utc=True)
        # Restrict obs lookup to provids in this batch.
        provids = geo["provid"].unique()
        obs_part = obs[obs["provid"].isin(provids)]
        merged = geo.merge(obs_part, on=["provid", "obstime"], how="left")
        state = pd.DataFrame(
            {
                "id": pd.to_numeric(merged["id"], errors="coerce"),
                "provid_vro": merged["provid"],
                "provid_mpc": merged["mpc_designation"].astype(str).str.strip(),
                "permid": merged["permid"].astype("string"),
                "obstime": merged["obstime"].dt.strftime("%Y-%m-%d %H:%M:%S.%f").str[:-3],
                "obs_et": merged["obs_et"].astype(float),
                "deltat": merged["dt_from_epoch_days"].astype(float),
                "r": merged["r"].astype(float),
                "delta": merged["delta"].astype(float),
                "alpha": merged["phase_angle"].astype(float),
                "mag": merged["mag"].astype(float),
                "H": merged["H"].astype(float),
                "flag_updated": np.int64(0),
            }
        )
        table = pa.Table.from_pandas(state, preserve_index=False)
        # Force a stable schema (id may be all-null in early batches).
        schema = pa.schema(
            [
                ("id", pa.float64()),
                ("provid_vro", pa.string()),
                ("provid_mpc", pa.string()),
                ("permid", pa.string()),
                ("obstime", pa.string()),
                ("obs_et", pa.float64()),
                ("deltat", pa.float64()),
                ("r", pa.float64()),
                ("delta", pa.float64()),
                ("alpha", pa.float64()),
                ("mag", pa.float64()),
                ("H", pa.float64()),
                ("flag_updated", pa.int64()),
            ]
        )
        table = table.cast(schema)
        if writer is None:
            writer = pq.ParquetWriter(OUT / "db_mvroc_state.parquet", schema)
        writer.write_table(table)
        n_written += len(state)
        _print(f"  state batch {i + 1}: {n_written:,} rows")
    if writer is not None:
        writer.close()
    _print(f"  {n_written:,} rows → db_mvroc_state.parquet")


def _nobs_provids() -> set[str]:
    """Provids listed under figures/figures_nobs (including updateH/)."""
    frames = []
    main = NOBS / "asteroid_analysis_results.csv"
    if not main.exists():
        raise SystemExit(f"Missing MOI list: {main}")
    frames.append(pd.read_csv(main, usecols=["provid"]))
    update_h = NOBS / "updateH" / "asteroid_analysis_results.csv"
    if update_h.exists():
        frames.append(pd.read_csv(update_h, usecols=["provid"]))
    provids: set[str] = set()
    for frame in frames:
        provids.update(frame["provid"].astype(str).str.strip().tolist())
    return provids


def build_moi() -> Path:
    """MOI = figures_nobs objects; rows = their observations from the source parquet."""
    _print("Building db_mvroc_moi.parquet …")
    if not OBS_SRC.exists():
        raise SystemExit(f"Missing source observations: {OBS_SRC}")
    provids = _nobs_provids()
    _print(f"  {len(provids):,} MOI provids from figures_nobs")

    out_path = OUT / "db_mvroc_moi.parquet"
    old_csv = OUT / "db_mvroc_moi.csv"
    pf = pq.ParquetFile(OBS_SRC)
    schema = pf.schema_arrow
    writer = None
    n_rows = 0
    for i, batch in enumerate(pf.iter_batches(batch_size=400_000), start=1):
        df = batch.to_pandas()
        mask = df["provid"].notna() & df["provid"].astype(str).str.strip().isin(
            provids
        )
        part = df.loc[mask]
        if part.empty:
            continue
        table = pa.Table.from_pandas(part, schema=schema, preserve_index=False)
        if writer is None:
            writer = pq.ParquetWriter(out_path, schema, compression="zstd")
        writer.write_table(table)
        n_rows += len(part)
        _print(f"  moi batch {i}: {n_rows:,} rows")
    if writer is not None:
        writer.close()
    if old_csv.exists():
        old_csv.unlink()
    _print(f"  {n_rows:,} rows → db_mvroc_moi.parquet")
    return out_path


def build_moc_sol() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """MOC + SOL from figures_MOI / figures_mvroc_sol.

    Returns (object_meta, moc, sol). object_meta carries per-object fields used
    by mvroc.csv (from the MOC selection, not the observation-level MOI).
    Fold metrics and colours come from select_MOC_solutions (no shape fits).
    """
    _print("Building db_mvroc_moc / sol …")
    selection = pd.read_csv(MOI / "asteroid_selection.csv")
    selection["provid"] = selection["provid"].astype(str).str.strip()
    solutions = expand_solutions(pd.read_csv(MOI / "solutions.csv"))
    results = pd.read_csv(FOLD / "asteroid_analysis_results.csv")
    results["provid"] = results["provid"].astype(str).str.strip()
    results["solution_id"] = results["solution_id"].astype(str).str.strip()

    object_meta = selection.copy()
    n_sol = solutions.groupby("provid").size().rename("n_solutions")
    object_meta = object_meta.merge(n_sol, on="provid", how="left")
    object_meta["n_solutions"] = object_meta["n_solutions"].fillna(0).astype(int)
    periods = solutions.groupby("provid")["p_rot_h"].apply(
        lambda s: "/".join(f"{p:.6f}" for p in s)
    )
    sol_ids = solutions.groupby("provid")["solution_id"].apply(
        lambda s: "/".join(str(x) for x in s)
    )
    object_meta = object_meta.merge(periods.rename("periods_h"), on="provid", how="left")
    object_meta = object_meta.merge(
        sol_ids.rename("solutions"), on="provid", how="left"
    )

    moc = object_meta[
        ["provid", "n_obs", "n_solutions", "solutions", "periods_h", "p_rot_dem_h"]
    ].copy()
    moc["selection"] = "visual_folding"
    moc.to_csv(OUT / "db_mvroc_moc.csv", index=False)
    _print(f"  {len(moc):,} rows → db_mvroc_moc.csv")

    sol = solutions.merge(
        results,
        on=["provid", "solution_id"],
        how="left",
        suffixes=("", "_res"),
    )
    if "p_rot_h_res" in sol.columns:
        sol["p_rot_h"] = sol["p_rot_h"].fillna(sol["p_rot_h_res"])
        sol = sol.drop(columns=["p_rot_h_res"])
    sol["tags"] = ""
    sol["Q"] = np.nan
    sol["use_band"] = ""
    sol.to_csv(OUT / "db_mvroc_sol.csv", index=False)
    _print(f"  {len(sol):,} rows → db_mvroc_sol.csv")
    return object_meta, moc, sol


def build_catalogue(
    object_meta: pd.DataFrame,
    sol: pd.DataFrame,
    orbits: pd.DataFrame,
    db_vro_mpc: pd.DataFrame,
) -> pd.DataFrame:
    _print("Building mvroc.csv …")
    colors = pd.read_parquet(DATA / "colors.parquet")
    colors["provid"] = colors["provid"].astype(str).str.strip()
    colors = colors.drop_duplicates("provid", keep="first")

    sbdb = pivot_sbdb(DATA / "sbdb_physical_properties.csv")
    gst = greenstreet_periods(DATA / "asteroids-with-rotation.csv")

    # Object-level fields from the MOC selection (not observation-level MOI).
    obj = object_meta.copy()
    # Mean H / designation from colors.parquet
    h_cols = ["H_u", "H_g", "H_r", "H_i", "H_z", "H_y"]
    color_obj = colors[["provid", "mpc_designation", "permid", *h_cols]].rename(
        columns={"mpc_designation": "mpc_designation", "permid": "permid"}
    )
    obj = obj.merge(color_obj, on="provid", how="left")
    # Fill designation / permid from cross-match when colors miss a MOI.
    xmatch = db_vro_mpc.rename(
        columns={
            "provid_vro": "provid",
            "provid_mpc": "mpc_designation_x",
            "permid": "permid_x",
        }
    )
    obj = obj.merge(
        xmatch[["provid", "mpc_designation_x", "permid_x"]], on="provid", how="left"
    )
    if "mpc_designation" not in obj.columns:
        obj["mpc_designation"] = np.nan
    if "permid" not in obj.columns:
        obj["permid"] = np.nan
    obj["mpc_designation"] = obj["mpc_designation"].fillna(obj["mpc_designation_x"])
    obj["permid"] = obj["permid"].fillna(obj["permid_x"])
    # Orbits: match by provid_vro first, then by mpc designation.
    orb_vro = orbits.rename(columns={"provid_vro": "provid"}).drop_duplicates(
        "provid", keep="first"
    )
    orb_cols = ["epoch_tdb", "a", "e", "q", "incl", "node", "peri", "M"]
    obj = obj.merge(
        orb_vro[["provid", *orb_cols]],
        on="provid",
        how="left",
        suffixes=("", "_orb"),
    )
    missing = obj["a"].isna()
    if missing.any():
        orb_mpc = orbits.drop_duplicates("provid_mpc", keep="first")
        fill = obj.loc[missing, ["provid", "mpc_designation"]].merge(
            orb_mpc,
            left_on="mpc_designation",
            right_on="provid_mpc",
            how="left",
        )
        for col in orb_cols:
            obj.loc[missing, col] = fill[col].to_numpy()

    # One catalogue row per solution.
    cat = sol.merge(obj, on="provid", how="left", suffixes=("", "_obj"))
    # Prefer solution-level fold metrics over object-level duplicates.
    for col in (
        "p_rot_h",
        "chi2_nu",
        "g-r",
        "g-r_err",
        "r-i",
        "r-i_err",
        "i-z",
        "i-z_err",
    ):
        obj_col = f"{col}_obj"
        if obj_col in cat.columns:
            cat[col] = cat[col].fillna(cat[obj_col])

    # Rename shape columns to catalogue (_k) names.
    cat = cat.rename(columns=SHAPE_RENAME)

    # Literature cross-matches.
    cat = cat.merge(sbdb, on="provid", how="left")
    cat = cat.merge(gst, on="provid", how="left")

    # Multiband Lk amplitudes / residuals (not produced in the MOI folding run).
    for col in ("amp_lsm_k", "amp_photo_k", "rms_multi", "n_outliers_multi"):
        if col not in cat.columns:
            cat[col] = np.nan
    if "tags" not in cat.columns:
        cat["tags"] = ""
    if "Q" not in cat.columns:
        cat["Q"] = np.nan
    if "use_band" not in cat.columns:
        cat["use_band"] = ""

    cat = add_dynamical_columns(cat)
    cat = add_taxonomy_columns(cat)

    # Keep only catalogue columns, in documented order.
    for col in CATALOGUE_COLUMNS:
        if col not in cat.columns:
            cat[col] = np.nan
    cat = cat[CATALOGUE_COLUMNS].copy()
    cat = cat.sort_values(["provid", "solution_id"]).reset_index(drop=True)
    cat.to_csv(OUT / "mvroc.csv", index=False)
    _print(f"  {len(cat):,} rows × {cat.shape[1]} cols → mvroc.csv")
    return cat


def main(argv: list[str] | None = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-state", action="store_true")
    parser.add_argument("--skip-orbits", action="store_true")
    parser.add_argument("--skip-crossmatch", action="store_true")
    parser.add_argument("--skip-moi", action="store_true")
    args = parser.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    if args.skip_crossmatch and (OUT / "db_vro_mpc.csv").exists():
        db_vro_mpc = pd.read_csv(OUT / "db_vro_mpc.csv")
        _print(f"Reusing db_vro_mpc.csv ({len(db_vro_mpc):,} rows)")
    else:
        db_vro_mpc = build_db_vro_mpc()
        build_info_vro_mpc(db_vro_mpc)

    if args.skip_orbits and (OUT / "db_mvroc_orbits.csv").exists():
        orbits = pd.read_csv(OUT / "db_mvroc_orbits.csv")
        _print(f"Reusing db_mvroc_orbits.csv ({len(orbits):,} rows)")
    else:
        orbits = build_orbits(db_vro_mpc)

    if not args.skip_state:
        build_state()
    else:
        _print("Skipping db_mvroc_state.parquet")

    if not args.skip_moi:
        build_moi()
    else:
        _print("Skipping db_mvroc_moi.parquet")

    object_meta, moc, sol = build_moc_sol()
    build_catalogue(object_meta, sol, orbits, db_vro_mpc)
    _print("Done.")


if __name__ == "__main__":
    main()
