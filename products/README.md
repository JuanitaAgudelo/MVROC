# MVROC pipeline products — summary

Result of the MVROC pipeline run of **4 November 2026**.

This directory sits at the **repository root** (`products/`). Raw observation
tables live under [`../pipeline/inputs/`](../pipeline/inputs/). The large parquet tables are in
[this Google Drive folder](https://drive.google.com/drive/folders/1l99VF8mxHfaqc82aR4SSCKunrfRnXzCZ?usp=sharing).

Input observations downloaded **12 October 2026** (8 078 546 rows).

| File | Format | Rows | Role |
|---|---|---:|---|
| `db_vro_mpc.csv` | CSV | 300 952 | Rubin ↔ MPC equivalence |
| `db_mvroc_orbits.csv` | CSV | 298 225 | Orbital elements per MPC object |
| `db_mvroc_state.parquet` | Parquet | 8 021 267 | Per-observation geometry + H. On Google Drive. |
| `db_mvroc_moi.parquet` | Parquet | 1 962 044 obs<br>(8 000 objects) | MOI observations (same columns as the source table). On Google Drive. |
| `db_mvroc_moc.csv` | CSV | 703 | MOC: objects with candidate light curves (visual selection) |
| `db_mvroc_sol.csv` | CSV | 1 092 | Rotation solutions (one row / solution) |
| `mvroc.csv` | CSV | 1 092 | Public catalogue join (one row / solution) |

---

## `db_vro_mpc.csv`

Equivalence table between Rubin designations and MPC designations.

| Column | Meaning |
|---|---|
| `provid_vro` | Rubin provisional designation |
| `provid_mpc` | Matched MPC designation |
| `permid` | Permanent number when available |
| `match_kind` | `numbered` or `provisional` |
| `n_vro_for_mpc` | How many Rubin ids map to this MPC designation |

**Counts:** 300 951 unique `provid_vro`, 297 966 unique `provid_mpc`. About 2 986 MPC objects gather two Rubin provids.

---

## `db_mvroc_orbits.csv`

Orbital elements at the reported epoch for each `provid_mpc`.

| Column | Unit | Meaning |
|---|---|---|
| `provid_vro` | | One Rubin id attached to the orbit row |
| `provid_mpc` | | MPC designation (unique key) |
| `permid` | | Permanent number |
| `epoch_tdb` | TDB | Epoch of the elements |
| `a` | au | Semimajor axis |
| `e` | | Eccentricity |
| `q` | au | Perihelion distance |
| `incl` | deg | Inclination |
| `node` | deg | Longitude of ascending node |
| `peri` | deg | Argument of perihelion |
| `M` | deg | Mean anomaly |

---

## `db_mvroc_state.parquet`

One row per observation with geometry and absolute magnitude.
Primary key: observation `id`. The file is in the [Google Drive folder](https://drive.google.com/drive/folders/1l99VF8mxHfaqc82aR4SSCKunrfRnXzCZ?usp=sharing) (≈ 472 MB).

| Column | Unit | Meaning |
|---|---|---|
| `id` | | Observation id from the Rubin MPC table (~99.3% filled) |
| `provid_vro` | | Rubin designation |
| `provid_mpc` | | Matched MPC designation |
| `permid` | | Permanent number |
| `obstime` | UTC | Observation time |
| `obs_et` | s | Ephemeris time |
| `deltat` | d | Days from orbital epoch to observation (TDB) |
| `r` | au | Heliocentric distance |
| `delta` | au | Geocentric distance |
| `alpha` | deg | Phase angle |
| `mag` | mag | Apparent magnitude |
| `H` | mag | Absolute magnitude (G = 0.15) |
| `flag_updated` | | Recompute flag for incremental updates |

---

## `db_mvroc_moi.parquet`

**MOI — MVROC Objects of Interest.** The N_obs ≥ 100 selection, including the absolute-magnitude update. These are not the visually selected rotation candidates (those are MOC).

One row per observation, with the same 90 columns as the source observation table. The file is in the [Google Drive folder](https://drive.google.com/drive/folders/1l99VF8mxHfaqc82aR4SSCKunrfRnXzCZ?usp=sharing) (≈ 143 MB). This build:

| Quantity | Value |
|---|---|
| Distinct objects (`provid_mpc`) | 8 000 |
| Distinct `provid` | 8 010 |
| Observation rows | 1 962 044 |
| File size | ≈ 143 MB |

Ten MPC designations gather two Rubin `provid` each; counting unique `provid_mpc` gives 8 000 physical objects.

Photometric periods, colours, and shape fits are **not** stored here; they
appear later in `db_mvroc_moc.csv`, `db_mvroc_sol.csv`, and `mvroc.csv`.

---

## `db_mvroc_moc.csv`

**MOC — MVROC Objects with Candidate light curves.** Objects whose folded light
curves were retained after visual review (703 asteroids in this build).

Columns: `provid`, `n_obs`, `n_solutions`, `solutions`, `periods_h`,
`p_rot_dem_h`, `selection`.

**Sample (this build):** N_obs from 100 to 765 (median 351).
Photometry span 2025-04-22 → 2026-04-15.

---

## `db_mvroc_sol.csv`

MVROC Solutions: one row per retained rotation period (1 092 solutions for the
703 MOC objects). Includes fold diagnostics, colours from the joint fit, and
shape-fit parameters (contact binary + ellipsoid).

| Solutions per object | Objects |
|---:|---:|
| 1 | 386 |
| 2 | 248 |
| 3 | 66 |
| 4 | 3 |

**Period range:** 0.063–44.1 h (median ≈ 4.85 h).

Classification columns `tags`, `Q`, and `use_band` are reserved and empty in this build.

---

## `mvroc.csv`

End product: **one row per solution**, with the object information joined onto each solution.

| Quantity | Value |
|---|---|
| Rows (solutions) | 1 092 |
| Columns | 99 |
| Unique objects | 703 |
| Orbits filled (`a`) | 100% |
| H_r filled | 98.8% |
| Shape `q_sph_k` filled | 99.5% |
| Greenstreet period match | 56 objects |
| SBDB `rot_per` match | 35 objects (44 solutions) |

### Column blocks

1. **Identification** — `provid`, `solution_id`, `permid`, `mpc_designation`
2. **Observations** — `epoch_utc_ini`, `epoch_utc_end`, `n_obs`, `n_obs_u`, `n_obs_g`, `n_obs_r`, `n_obs_i`, `n_obs_z`, `n_obs_y`
3. **Orbit** — `epoch_tdb`, `a`, `e`, `q`, `incl`, `node`, `peri`, `M`
4. **Mean absolute magnitudes** — `H_u`, `H_g`, `H_r`, `H_i`, `H_z`, `H_y`
5. **Phase-corrected colours** — `u-g`, `u-g_err`, `u-r`, `u-r_err`, `g-r`, `g-r_err`, `r-i`, `r-i_err`, `i-z`, `i-z_err`, `z-y`, `z-y_err`
6. **Photometric taxonomy** — `gri_slope`, `photo_spec_T`
7. **Adopted solution** — `p_rot_h`, `chi2_nu`
8. **Multiband Lk** — `amp_lsm_k`, `amp_photo_k`, `rms_multi`, `n_outliers_multi`
9. **Shape, spheres** — `q_sph_k`, `q_sph_err_k`, `phi_sph_deg_k`, `phi_sph_err_deg_k`, `chi_r_sph_k`, `amp_sph_k`, `m_mid_sph_k`, `m_sep_k`, `p_rot_sph_h_k`
10. **Shape, ellipsoid** — `b_over_a_k`, `b_over_a_err_k`, `amp_ell_k`, `m_mid_ell_k`, `m_ab_k`, `p_rot_ell_h_k`, `chi_r_ell_k`
11. **Classification** — `tags`, `Q`, `use_band`
12. **SBDB** — `H_sbdb`, `H_sbdb_err`, `albedo_sbdb`, `albedo_sbdb_err`, `diameter_sbdb`, `diameter_sbdb_err`, `extent_sbdb`, `extent_sbdb_err`, `pole_sbdb`, `pole_sbdb_err`, `density_sbdb`, `density_sbdb_err`, `rot_per_sbdb`, `rot_per_sbdb_err`, `GM_sbdb`, `GM_sbdb_err`, `BV_sbdb`, `BV_sbdb_err`, `UB_sbdb`, `UB_sbdb_err`, `spec_B_sbdb`, `spec_B_sbdb_err`, `spec_T_sbdb`, `spec_T_sbdb_err`
13. **Greenstreet** — `rot_per_lsm_gst`, `rot_per_fr_gst`
14. **Dynamics** — `dyn_class`, `T_earth`, `T_mars`, `T_jup`, `T_nep`, `flag_resonant`, `resonant_ratio`

-----
[MVROC Team](https://github.com/seap-udea/MVROC) (2026-present)
