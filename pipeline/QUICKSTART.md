# MVROC pipeline quick start

Initialize a **fresh clone** and run the stage drivers in order. All commands
below use **cwd = `pipeline/`**.

For driver style, paths, and column schemas see [`PATTERN.md`](PATTERN.md).

---

## 1. Clone and enter the pipeline

```bash
git clone git@github.com:seap-udea/MVROC.git
cd MVROC/pipeline
```

Use a Python 3 environment of your choice (venv, conda, …).

---

## 2. Install dependencies

```bash
make install
```

This installs [`.bin/requirements.txt`](.bin/requirements.txt)
(`numpy`, `pandas`, `astropy`, `gdown`, …).

**Without `make`** (e.g. Windows): the Makefile targets only run Python
scripts, so run them directly from `pipeline/`:

| `make` target | Equivalent command |
|---|---|
| `make install` | `python -m pip install -r .bin/requirements.txt` |
| `make get_gdrive_data` | `python .bin/get_gdrive_data.py` |
| `make get_public_data` | `python .bin/get_public_data.py` |

---

## 3. Download data

Large files are gitignored. On a fresh clone, pull the pipeline data set from
Drive:

```bash
make get_gdrive_data
```

That fills `inputs/` (obs + MPCORB) and early products
(`outputs/mois`, `outputs/states`; `states` can also be rebuilt with
`states.py`) from
[mvroc-pipeline-data](https://drive.google.com/drive/folders/1gddbkDYiGWRI4pVhKPfR5AGp5dNpTrBW?usp=drive_link).
The folder must be shared (*Anyone with the link / Viewer*). Existing local
files are **skipped**.

Check:

```bash
ls -lh inputs/
ls -lh outputs/mois/ outputs/states/
```

### Optional: refresh public inputs

Only when you want to **update** raw inputs from their public sources (not
needed for a normal first setup):

```bash
make get_public_data    # re-fetch inputs/obs_sbn… + inputs/mpcorb…
```

Delete the local files first if you need a forced refresh (`get_public_data`
also skips paths that are already present).

### SPICE kernels

`states.py` needs three SPICE kernels that no `make` target downloads. Copy
them from `code/kernels/` (gitignored under `inputs/` too):

```bash
mkdir -p inputs/kernels
cp ../code/kernels/naif0012.tls ../code/kernels/gm_de440.tpc ../code/kernels/de440s.bsp inputs/kernels/
```

| File | Role |
|---|---|
| `inputs/kernels/naif0012.tls` | Leap seconds (UTC -> ET) |
| `inputs/kernels/gm_de440.tpc` | GM of the Sun (DE440) |
| `inputs/kernels/de440s.bsp` | Planetary ephemeris (Earth position) |

---

## 4. Stage order

The pipeline is sequential. Each driver reads earlier outputs only.

```text
inputs
  → crossmatches → orbits → states → mois
  → periodograms → foldings → mocs → solutions → shapefit
```

| Order | Stage | Driver today | Notes |
|------:|-------|--------------|-------|
| 1 | `crossmatches` | `crossmatches.py` | needs obs + MPCORB (+ MPC API) |
| 2 | `orbits` | `orbits.py` | needs `crossmatches.csv` + MPCORB |
| 3 | `states` | `states.py` | needs `orbits.csv`, obs and `inputs/kernels/` (~4 min for all) |
| 4 | `mois` | (not yet) | parquet via `make get_gdrive_data` |
| 5 | `periodograms` | `periodograms.py` | needs photometry (mois+states+crossmatches) |
| 6 | `foldings` | `foldings.py` | needs `periodograms.csv` |
| 7 | `mocs` | (not yet) | visual selection → `mocs.csv` |
| 8 | `solutions` | `solutions.py` | needs `mocs.csv` |
| 9 | shapefit | (future) | |

Until their drivers exist, keep simulated / Drive-provided files under
`outputs/mois` and `outputs/mocs`.

---

## 5. Run the available drivers

One positional argument: a provid or `all`.

Build the early stages for the whole catalogue (in this order):

```bash
python crossmatches.py all
python orbits.py all
python states.py all
```

With a single provid these three print a per-object summary (MPC
designation and Rubin aliases, orbital elements, geometry ranges and mean H).
`crossmatches.csv` and `orbits.csv` are appended like the other CSVs;
`states.parquet` is rewritten on every run and only objects whose orbit
changed are recomputed (`flag_updated`, `orbit_hash`).

Smoke-test one asteroid end-to-end (after data + upstream CSVs are in place):

```bash
python periodograms.py "1981 US22"
python foldings.py "1981 US22"
python solutions.py "1981 US22"
```

Underscores in the name are fine (`1981_US22`). Process every object listed
in the stage upstream CSV:

```bash
python periodograms.py all
python foldings.py all
python solutions.py all
```

Each run **appends** to `outputs/<stage>/<stage>.csv` and writes a log plus
optional WebP figures. Re-running the same provid **skips** if that provid is
already in the CSV (delete the CSV to start clean).

Typical products after the smoke-test:

```text
outputs/periodograms/periodograms.csv
outputs/periodograms/1981_US22_periodogram.webp
outputs/foldings/foldings.csv
outputs/foldings/1981_US22_folding.webp
outputs/solutions/solutions.csv
outputs/solutions/1981_US22_solution.webp
```

`solutions.py` reads retained periods from `outputs/mocs/mocs.csv` (not from
foldings directly).

---

## 6. Minimal checklist

- [ ] `cd pipeline`
- [ ] `make install`
- [ ] `make get_gdrive_data`
- [ ] Kernels copied to `inputs/kernels/`
- [ ] `python crossmatches.py all`, `python orbits.py all`,
      `python states.py all`
- [ ] Upstream files present (`mois`, and for solutions: `mocs`)
- [ ] `python periodograms.py "1981 US22"`
- [ ] `python foldings.py "1981 US22"`
- [ ] `python solutions.py "1981 US22"`

More detail: [`PATTERN.md`](PATTERN.md).

-----
[MVROC Team](https://github.com/seap-udea/MVROC) (2026-present)
