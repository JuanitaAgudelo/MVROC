# MVROC inputs

Raw downloaded observation tables used by the pipeline. Large parquet files
are not committed; place them here (or keep symlinks to `code/data/`).

| File | Role |
|---|---|
| `obs_sbn_X05_full.parquet` | Rubin MPC observations (Asteroid Institute / B612 download) |
| `obs_sbn_X05_full_mpc_desig.parquet` | Same table with MPC designations attached |

Download (example):

```bash
# from the repository root
mkdir -p pipeline/inputs
wget -P pipeline/inputs/ \
  https://storage.googleapis.com/asteroid-institute-public/production/rubin/mpc/obs_sbn/full/parquet/obs_sbn_X05_full.parquet
```

Pipeline drivers run from `pipeline/` and read `inputs/...`.
Derived tables (crossmatches, orbits, states, mois, catalogue) live under
[`../../products/`](../../products/).
