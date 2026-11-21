# MVROC pipeline driver pattern

How to write a stage driver under `pipeline/`. Global variables and helpers of the pipeline: `mvroc.py` . Reference implementations:
`crossmatches.py`, `orbits.py`, `states.py`, `periodograms.py`, `foldings.py`,
`solutions.py`.

---

## 1. What a driver is

One stage = one script = one scientific step.

| Piece | Role |
|-------|------|
| `pipeline/<stage>.py` | Science, stage constants, plot, CSV/log/CLI |
| `pipeline/mvroc.py` | Shared helpers used by **more than one** stage |
| `pipeline/outputs/<stage>/` | Everything this stage writes |

Do not invent a new architecture. Copy the human style of the reference drivers.

---

## 2. Pipeline layout

Always run with **cwd = `pipeline/`**. Paths are relative strings from there.

```text
pipeline/                       # cwd for every driver
  inputs/                       # raw downloads (not committed to git)
    obs_sbn_X05_full.parquet    # Rubin MPC observations
    mpcorb_extended.dat.gz      # MPCORB extended (MPC)
    kernels/                    # SPICE kernels (used by states)
      naif0012.tls              #   leap seconds
      gm_de440.tpc              #   GM of the Sun (DE440)
      de440s.bsp                #   planetary ephemeris (Earth position)
  outputs/                      # stage outputs only
    crossmatches/
    orbits/
    states/
    mois/
    periodograms/
    foldings/
    mocs/
    solutions/
  mvroc.py
  crossmatches.py
  orbits.py
  states.py
  periodograms.py
  foldings.py
  solutions.py
  PATTERN.md
```

| Path from `pipeline/` | Constant in `mvroc.py` | Role |
|-----------------------|------------------------|------|
| `inputs/` | `INPUTS` | Raw downloads (obs + MPCORB + kernels) |
| `inputs/kernels/` | `KERNELS` | SPICE kernels |
| `outputs/` | `OUTPUTS` | Per-stage pipeline outputs |

**Do not use** `pathlib.Path.resolve()`, `__file__`, or path helper wrappers
(`set_fig_dir`, `results_path`, …).

### Getting the data

Large files are **gitignored**. From **`pipeline/`**:

```bash
cd pipeline
make install
make get_gdrive_data
```

**`make get_gdrive_data`** pulls the mirrored tree from
[mvroc-pipeline-data](https://drive.google.com/drive/folders/1gddbkDYiGWRI4pVhKPfR5AGp5dNpTrBW?usp=drive_link)
and **skips** paths that already exist locally:

| Local path (under `pipeline/`) | Role |
|--------------------------------|------|
| `inputs/obs_sbn_X05_full.parquet` | Rubin MPC observations (X05) |
| `inputs/mpcorb_extended.dat.gz` | MPCORB extended orbital elements |
| `outputs/mois/mois.parquet` | MOI-filtered observations |
| `outputs/states/states.parquet` | Geometry + absolute magnitude H |

Share that Drive folder as *Anyone with the link / Viewer* so `gdown` can
list and download it.

**SPICE kernels** are not downloaded by any `make` target. Copy them from
`code/kernels/` into `inputs/kernels/` (also gitignored):

| Local path (under `pipeline/`) | Role |
|--------------------------------|------|
| `inputs/kernels/naif0012.tls` | Leap seconds (UTC -> ET) |
| `inputs/kernels/gm_de440.tpc` | GM of the Sun (DE440) |
| `inputs/kernels/de440s.bsp` | Planetary ephemeris (Earth position) |

Without `make` (e.g. on Windows), run the same scripts from `pipeline/`:
`python -m pip install -r .bin/requirements.txt`,
`python .bin/get_gdrive_data.py`, `python .bin/get_public_data.py`.

**Optional — update public inputs only** (not part of a normal first setup).
Use when refreshing obs / MPCORB from their public URLs:

```bash
make get_public_data
```

| Local path (under `pipeline/`) | Source |
|--------------------------------|--------|
| `inputs/obs_sbn_X05_full.parquet` | Asteroid Institute public bucket |
| `inputs/mpcorb_extended.dat.gz` | Minor Planet Center |

Also skips existing local files; delete them first to force a refresh.

Column schemas: Appendix A.1–A.2b (inputs), A.3–A.6 (crossmatches / orbits /
states / mois).

---

## 3. Upstream and downstream

The pipeline is **sequential**. A driver assumes earlier stages already ran.

Before writing a driver, declare:

1. **Stage name** — short plural or stage noun (`periodograms`, `foldings`, `solutions`).
2. **Upstream** — stage outputs (or `inputs/`) that must exist first.
3. **Downstream** — stages that will read this stage later.

### Photometry (shared join)

Several photometric stages do **not** read `mois.parquet` alone. They call
`mv.load_moi_photometry()`, which builds one light-curve table by joining:

| Source | Path | Columns used |
|--------|------|--------------|
| MOI observations | `outputs/mois/mois.parquet` | `id`, `obstime`, `permid`, `provid`, `ra`, `dec`, `band`, `mag`, `rmsmag` |
| Absolute mag / geometry | `outputs/states/states.parquet` | `id`, `H`, `provid_mpc` (join on `id`) |
| Dual-provid aliases | `outputs/crossmatches/crossmatches.csv` | `provid_vro`, `provid_mpc` (merge aliases) |

Result columns: `id`, `obstime`, `permid`, `provid`, `ra`, `dec`, `band`,
`mag`, `rmsmag`, `H`, `provid_mpc`, `mpc_designation`, `night`.
Bands are normalized to `Lu`/`Lg`/`Lr`/`Li`/`Lz`/`Ly`. See Appendix A.9.

In this document, **photometry** always means that join — never a free-standing
file named `photometry.*`.

### Current stage order

| Order | Stage | Upstream (exact paths) | Writes |
|------:|-------|------------------------|--------|
| 0 | seed | external downloads (see **Getting the data** above) | `inputs/obs_sbn_X05_full.parquet`, `inputs/mpcorb_extended.dat.gz` |
| 1 | `crossmatches` | `inputs/obs_sbn_X05_full.parquet`, `inputs/mpcorb_extended.dat.gz`, MPC Identifier API | `outputs/crossmatches/crossmatches.csv`, `mpc_query_identifier.json`, `crossmatches.log` |
| 2 | `orbits` | `outputs/crossmatches/crossmatches.csv`, `inputs/mpcorb_extended.dat.gz` | `outputs/orbits/orbits.csv`, `orbits.log` |
| 3 | `states` | `outputs/orbits/orbits.csv`, `inputs/obs_sbn_X05_full.parquet`, `inputs/kernels/` | `outputs/states/states.parquet`, `states.log` |
| 4 | `mois` | `outputs/states/states.parquet` (+ MOI filter on obs) | `outputs/mois/mois.parquet` |
| 5 | `periodograms` | photometry = mois + states + crossmatches | `outputs/periodograms/periodograms.csv`, `*_periodogram.webp`, `periodograms.log` |
| 6 | `foldings` | `outputs/periodograms/periodograms.csv` + photometry | `outputs/foldings/foldings.csv`, `*_folding.webp`, `foldings.log` |
| 7 | `mocs` | `outputs/foldings/foldings.csv` (visual / selected periods) | `outputs/mocs/mocs.csv` |
| 8 | `solutions` | `outputs/mocs/mocs.csv` + photometry | `outputs/solutions/solutions.csv`, `*_solution.webp`, `solutions.log` |
| 9 | shapefit | `outputs/solutions/solutions.csv` | (future) |

`crossmatches`, `orbits` and `states` have drivers. Until the remaining
early-stage drivers exist, their outputs may be **simulated** or come from
Drive under `outputs/<stage>/` (currently `mois` and `mocs`). Column schemas
still apply (Appendix A).

### Rules

- Read earlier stage results from the **exact** upstream paths above.
- Read raw downloads from `inputs/` when needed
  (`mv.OBSERVATIONS_PARQUET`, `mv.MPCORB_DAT`, kernels under `mv.KERNELS`).
- Do **not** reimplement an upstream stage inside a later driver.
- If an upstream file is missing, exit with a clear message.
- Document Inputs / Outputs (including upstream paths) in the driver docstring.
- Split monolithic procedures into thin drivers along this order.

### Example chain

```text
inputs → crossmatches → orbits → states → mois
  → periodograms → foldings → mocs → solutions → shapefit
```

**`crossmatches.py`** reads `inputs/obs_sbn_X05_full.parquet` and
`inputs/mpcorb_extended.dat.gz` (plus the MPC Designation Identifier API for
objects not in MPCORB) and writes:

- `outputs/crossmatches/crossmatches.csv` (one row per matched `provid_vro`)
- `outputs/crossmatches/mpc_query_identifier.json` (cache of API answers)
- `outputs/crossmatches/crossmatches.log`

**`orbits.py`** reads `outputs/crossmatches/crossmatches.csv` and
`inputs/mpcorb_extended.dat.gz` and writes:

- `outputs/orbits/orbits.csv` (one row per `provid_vro`)
- `outputs/orbits/orbits.log`

**`states.py`** reads `outputs/orbits/orbits.csv`,
`inputs/obs_sbn_X05_full.parquet` and `inputs/kernels/` and writes:

- `outputs/states/states.parquet` (one row per observation id; rewritten
  every run, see section 9)
- `outputs/states/states.log`

**`periodograms.py`** reads photometry via `mv.load_moi_photometry()` and writes:

- `outputs/periodograms/periodograms.csv` (one row per asteroid)
- `outputs/periodograms/<provid>_periodogram.webp`
- `outputs/periodograms/periodograms.log`

**`foldings.py`** reads:

- `outputs/periodograms/periodograms.csv` (`periods_dem_h` vote tokens)
- photometry via `mv.load_moi_photometry()`

and writes:

- `outputs/foldings/foldings.csv` (one row per plotted two-peak solution)
- `outputs/foldings/<provid>_folding.webp`
- `outputs/foldings/foldings.log`

**`mocs`** (future driver) selects which folded periods are retained after
visual review and writes `outputs/mocs/mocs.csv`.

**`solutions.py`** reads:

- `outputs/mocs/mocs.csv` (`solutions`, `periods_h`)
- photometry via `mv.load_moi_photometry()`

and writes:

- `outputs/solutions/solutions.csv` (one row per retained solution)
- `outputs/solutions/<provid>_solution.webp`
- `outputs/solutions/solutions.log`

---

## 4. Layout and execution

```bash
cd pipeline
python crossmatches.py all
python orbits.py all
python states.py all
python periodograms.py "1981 US22"
python foldings.py all
python solutions.py "1981 US22"
```

Create the output directory once at driver start:

```python
os.makedirs(OUT_DIR, exist_ok=True)
```

Do not re-check paths inside the per-asteroid loop.

### Stage constants (top of the driver)

```python
OUT_DIR = "outputs/foldings"
RESULTS_CSV = f"{OUT_DIR}/foldings.csv"
LOG_FILE = f"{OUT_DIR}/foldings.log"
UPSTREAM_CSV = "outputs/periodograms/periodograms.csv"  # if needed
RESULTS_COLUMNS = [
    "provid",
    "selected",
    "p_rot_fold_h",
    # ...
]
```

Shared paths live in `mvroc.py`:

```python
INPUTS = "inputs"
OUTPUTS = "outputs"

OBSERVATIONS_PARQUET = "inputs/obs_sbn_X05_full.parquet"
MPCORB_DAT = "inputs/mpcorb_extended.dat.gz"
KERNELS = "inputs/kernels"  # SPICE kernels

CROSSMATCHES_CSV = "outputs/crossmatches/crossmatches.csv"
ORBITS_CSV = "outputs/orbits/orbits.csv"
STATES_PARQUET = "outputs/states/states.parquet"
MOIS_PARQUET = "outputs/mois/mois.parquet"
MOCS_CSV = "outputs/mocs/mocs.csv"
```

Photometry for photometric stages is loaded once with:

```python
df = mv.load_moi_photometry()
# joins mois.parquet + states.parquet + crossmatches.csv
# -> id, obstime, permid, provid, ra, dec, band, mag, rmsmag,
#    H, provid_mpc, mpc_designation, night
```

---

## 5. Naming

| Thing | Pattern | Example |
|-------|---------|---------|
| Driver file | `<stage>.py` | `periodograms.py`, `solutions.py` |
| Entry point | `mvroc_<stage>()` | `mvroc_periodograms()`, `mvroc_solutions()` |
| One-object science | `mvroc_<algo>_asteroid(...)` | `mvroc_periodogram_asteroid(...)` |
| Figure | `mvroc_<algo>_plot(...)` | `mvroc_folding_plot(...)`, `mvroc_solution_plot(...)` |
| Skip helper | `mv.already_in_outputs(provid, RESULTS_CSV)` | shared in `mvroc.py` |
| Skip set (many objects) | `mv.provids_in_outputs(RESULTS_CSV, column)` | reads the CSV once |
| Single-object printout | `print_<algo>(...)` | `print_crossmatch`, `print_orbit`, `print_state` |
| Arrays | plural `s` suffix | `ts`, `ys`, `dys` |
| Constants | `UPPER_CASE` | `NB`, `N_FREQ`, `MIN_MAG_ERROR` |

Repeated magic numbers become named constants (in the driver if stage-local;
in `mvroc.py` if shared). Prefer the literal `"H"` over a useless global
`PROP = "H"`.

Figure stems:

| Stage | Figure path |
|-------|-------------|
| periodograms | `outputs/periodograms/<provid>_periodogram.webp` |
| foldings | `outputs/foldings/<provid>_folding.webp` |
| solutions | `outputs/solutions/<provid>_solution.webp` |

---

## 6. `mvroc.py` vs the driver

**Put in `mvroc.py` only what a second stage already needs:**

- path constants (`INPUTS`, `OUTPUTS`, stage CSV/parquet paths)
- `load_moi_photometry`, `prepare_asteroid`, `observation_counts`, …
- `LogOutput`, `save_figure`, `slug`
- common helpers: `cadence_min_period`, `merge_frequencies`, `already_in_outputs`,
  `provids_in_outputs`
- shared constants: `MIN_MAG_ERROR`, `DPI`, …

**Keep in the driver:**

- stage science
- stage constants
- stage plot
- stage CSV / log / CLI

If unsure, leave it in the driver. Move to `mvroc.py` when a second stage
imports it for real. Prefer importing another driver's helpers only as a
temporary bridge; promote shared code into `mvroc.py` soon after.

---

## 7. File skeleton

```python
#!/usr/bin/env python3
"""One-line summary of the stage.

Stage: <stage>
Upstream: <upstream stages / paths>
Downstream: <who reads this next>

Usage (always from the pipeline/ directory):

    python <stage>.py "1981 US22"
    python <stage>.py 1981_US22
    python <stage>.py all
"""

import os
import sys
import time

import numpy as np
import pandas as pd

import mvroc as mv

# ---------------------------------------------------------------------------
# <Stage> constants
# ---------------------------------------------------------------------------
OUT_DIR = "outputs/<stage>"
RESULTS_CSV = f"{OUT_DIR}/<stage>.csv"
LOG_FILE = f"{OUT_DIR}/<stage>.log"
RESULTS_COLUMNS = ["provid", ...]

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def helper(...):
    """Summary.

    Parameters
    ----------
    ...

    Returns
    -------
    ...

    Examples
    --------
    >>> helper(...)
    ...
    """
    ...

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def mvroc_<algo>_plot(...):
    ...

# ---------------------------------------------------------------------------
# Science: one asteroid
# ---------------------------------------------------------------------------
def mvroc_<algo>_asteroid(df, provid, ..., qplot=True):
    ...

# ---------------------------------------------------------------------------
# <Stage> driver
# ---------------------------------------------------------------------------
def mvroc_<stage>():
    ...

if __name__ == "__main__":
    mvroc_<stage>()
```

Section headers sit **immediately** above the first routine (no blank line
after the `# ---` block). Use **one** blank line between routines — never two.

---

## 8. CLI

One positional argument only:

```text
python <stage>.py "<provid>" | python <stage>.py <provid_with_underscores> | python <stage>.py all
```

```python
# --> Check command-line arguments.
if len(sys.argv) != 2:
    print('Usage: python foldings.py "<provid>" | python foldings.py all')
    raise SystemExit(2)

# --> Get the target provid.
target = sys.argv[1].strip()
if target.lower() != "all":
    target = target.replace("_", " ")
```

- No argparse.
- No force / fig-dir / dead flags.
- If a parameter is almost always the same value, hard-code a constant.
- `qplot=True` is fine when it separates science from plotting.
- With a single provid, catalogue-wide stages (`crossmatches`, `orbits`,
  `states`) print a per-object summary (`print_<algo>`) instead of the
  catalogue summary printed with `all`.

---

## 9. Results CSV

- Path: `outputs/<stage>/<stage>.csv`
- At driver start: create with header if missing
- Per asteroid (or per solution): **append** with `header=False`
- Never rewrite or upsert rows — the user deletes the CSV to start clean
- Skip with `mv.already_in_outputs(provid, RESULTS_CSV)`

```python
# --> Initialize the results CSV if it doesn't exist.
if not os.path.isfile(RESULTS_CSV):
    pd.DataFrame(columns=RESULTS_COLUMNS).to_csv(RESULTS_CSV, index=False)
    print(f"Created {RESULTS_CSV}")

...
if mv.already_in_outputs(provid, RESULTS_CSV):
    print(f"{provid}: already in {os.path.basename(RESULTS_CSV)}, skip")
    continue

...
pd.DataFrame([row], columns=RESULTS_COLUMNS).to_csv(
    RESULTS_CSV, mode="a", header=False, index=False
)
```

Do not write incidental `*_observations.csv` files unless the stage product
truly needs them.

Stages keyed by `provid_vro` (`crossmatches`, `orbits`) cannot use
`mv.already_in_outputs` (it reads a `provid` column); they read the done set
once with `mv.provids_in_outputs(RESULTS_CSV, "provid_vro")`.

**Parquet exception (`states`).** A parquet file cannot be appended.
`states.py` rewrites `states.parquet` on every run: the `orbits.csv` row of
each object (`provid_mpc`, `permid`, epoch and MPCORB elements) is hashed
(`orbit_hash`); objects with an unchanged hash keep their rows
(`flag_updated = 0`) and new or changed objects are recomputed
(`flag_updated = 1`). The file is written to `states.parquet.tmp` and then
renamed.

**Schema guard.** Never append to (or merge with) a file written with other
columns: CSV stages call `mv.check_results_header(RESULTS_CSV,
RESULTS_COLUMNS)` at start; `states.py` checks the parquet schema the same
way. Both stop with "Delete it and run again".

---

## 10. Logging and comments

Mirror stdout to the stage log:

```python
# --> Create the output directory and set up logging.
os.makedirs(OUT_DIR, exist_ok=True)
log = open(LOG_FILE, "a", encoding="utf-8")
# --> Mirror all print() output to the log file and to the terminal.
sys.stdout = mv.LogOutput(sys.__stdout__, log)
```

Print style:

```text
====================
1981 US22
====================
  N = 670 observations
  Lr: N = 312, P_rot = 2.9377 h
  skip: no band with at least 5 points
```

- Detail lines indented with two spaces.
- Milestone comments use `# --> ...` (important procedure steps only).
- ASCII only in messages and docstrings (no fancy unicode arrows or dashes).

Example milestones inside science code:

```python
# --> Prepare multiband frame (drop Lu, reject magnitude outliers).
...
# --> Score each voted period with a joint 2-peak harmonic fit.
...
# --> One CSV row per plotted solution (selected + other two-peak panels).
```

---

## 11. Figures

- Save with `mv.save_figure(fig, path)` (WebP + MVROC watermark).
- Path: `f"{OUT_DIR}/{mv.slug(provid)}_<stem>.webp"`
- Keep **calculation** and **plot** in separate routines.
- Plotting is optional via `qplot=True` when the CSV row is the main product.

```python
path = f"{OUT_DIR}/{mv.slug(provid)}_periodogram.webp"
mv.save_figure(fig, path)
print(f"  democratic periodogram: {path}")
```

---

## 12. Docstrings

NumPy / human style. Full `Parameters` / `Returns` / `Examples` on helpers
and science routines. No double backticks. No `# doctest: +SKIP`.

```python
def cadence_min_period(ts):
    """Shortest resolvable photometric period: 2 x median delta t.

    Parameters
    ----------
    ts : array-like
        Observation times in days from an arbitrary zero point.

    Returns
    -------
    float
        Minimum photometric period in hours. 0.0 if fewer than two times.

    Examples
    --------
    >>> cadence_min_period([0.0, 0.01, 0.02, 0.03])
    0.48
    """
```

The entry point `mvroc_<stage>()` docstring must include:

- Stage / Upstream / Downstream (module docstring or driver docstring)
- Command-line
- Examples of invocation
- Inputs (upstream outputs + `inputs/` when needed)
- Outputs (CSV, figures, log)
- Returns

```python
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
    python solutions.py all

    Inputs
    ------
    outputs/mocs/mocs.csv
        Retained solution ids and periods (columns solutions, periods_h).
    Photometry via mv.load_moi_photometry()
        Join of mois.parquet + states.parquet + crossmatches.csv.

    Outputs
    -------
    outputs/solutions/solutions.csv
        One row per retained MVROC solution.
    outputs/solutions/<provid>_solution.webp
        Vertical stack of folded panels.
    outputs/solutions/solutions.log
        Same text as the terminal.

    Returns
    -------
    None
    """
```

---

## 13. Checklist before merge

- [ ] Lives in `pipeline/<stage>.py`
- [ ] Runs from `pipeline/`
- [ ] Upstream / Downstream declared
- [ ] Reads `outputs/<upstream>/` and/or `inputs/` as documented
- [ ] Does not commit large `inputs/` files (see `.gitignore`)
- [ ] Writes only `outputs/<stage>/`
- [ ] Entry point `mvroc_<stage>()`; stage helpers in the same file
- [ ] CLI: `provid | all`
- [ ] CSV append + `mv.already_in_outputs` (or `mv.provids_in_outputs`;
      parquet stages: hash-based rewrite as in `states.py`)
- [ ] Log via `mv.LogOutput`
- [ ] Full Parameters / Returns / Examples docstrings
- [ ] Important steps marked with `# -->`
- [ ] One blank line between routines; no blank after section headers
- [ ] No dead options, no Path resolve
- [ ] Smoke-test OK

---

## Appendix A. Files and columns

Paths relative to `pipeline/`. Schemas match the current driver outputs
(`crossmatches`, `orbits`, `states`) or the simulated / Drive products (the
rest). Color columns use LSST letters `u,g,r,i,z,y` (not `Lu`…).

### A.0 Color-column set (shared)

Used at the end of `foldings.csv` and `solutions.csv`. For every unordered
pair of filters in `(u, g, r, i, z, y)`:

| Column | Meaning |
|--------|---------|
| `<f1>-<f2>` | Mean color in magnitudes (NaN if either band missing) |
| `<f1>-<f2>_err` | Uncertainty on that color |

Pairs in order: `u-g`, `u-r`, `u-i`, `u-z`, `u-y`, `g-r`, `g-i`, `g-z`,
`g-y`, `r-i`, `r-z`, `r-y`, `i-z`, `i-y`, `z-y` (30 columns total).

---

### A.1 `inputs/obs_sbn_X05_full.parquet`

Rubin / Asteroid Institute MPC observations for station X05. Full ADPACS-like
observation schema (~90 columns). Photometric stages care about:

| Column | Meaning |
|--------|---------|
| `id` | Observation id (join key to states) |
| `provid` | Provisional / VRO designation |
| `permid` | Permanent number if numbered |
| `obstime` | Observation epoch |
| `ra`, `dec` | Sky position |
| `mag`, `rmsmag` | Apparent magnitude and error |
| `band` | Filter letter (`u`/`g`/`r`/`i`/`z`/`y` or `Lu`…) |
| `stn` | Observatory code (`X05`) |

Download: **Getting the data** (section 2) — `make get_gdrive_data`
(optional public refresh: `make get_public_data`).

---

### A.2 `inputs/mpcorb_extended.dat.gz`

MPCORB extended orbital-element file from the Minor Planet Center (gzip text).
Used by `crossmatches` (designation; other designations from column 217 in
fixed 11-character fields) and `orbits` (epoch, a, e, incl, node, peri, M);
not by photometric stages. Not a CSV;
parse with the MPCORB format. Download: **Getting the data** (section 2) —
`make get_gdrive_data` (optional: `make get_public_data`).

---

### A.2b `inputs/kernels/`

SPICE kernels used by `states` (copied from `code/kernels/`, gitignored):
`naif0012.tls` (leap seconds), `gm_de440.tpc` (GM of the Sun),
`de440s.bsp` (Earth position, ecliptic J2000).

---

### A.3 `outputs/crossmatches/crossmatches.csv`

VRO ↔ MPC designation map. One row per matched `provid_vro` (objects without
an MPC designation are not written and are retried on the next run). Search
order: `permid` in MPCORB, `provid` in MPCORB, `provid` among the MPCORB
other designations, MPC Designation Identifier API.

| Column | Meaning |
|--------|---------|
| `provid_vro` | VRO / provisional name used in observations |
| `provid_mpc` | MPC designation string (MPCORB format, e.g. `(1) Ceres`) |
| `permid` | MPC number (integer when numbered) |
| `match_kind` | Kind of MPC designation: `numbered`, `provisional` or `comet` |
| `n_vro_for_mpc` | How many VRO names map to this MPC object |
| `match_by` | Where the match came from: `permid`, `provid`, `other_desig` or `mpc_api` |

Also written under `outputs/crossmatches/`:

- `mpc_query_identifier.json` — cache of the MPC API answers
- `crossmatches.log` — mirrored stdout

---

### A.4 `outputs/orbits/orbits.csv`

MPCORB orbital elements of each cross-matched object. One row per
`provid_vro` (Rubin aliases of the same MPC object repeat its elements).
Objects with no orbit in MPCORB (mostly comets) get a row with empty epoch and
elements.

| Column | Meaning |
|--------|---------|
| `provid_vro` | VRO designation |
| `provid_mpc` | MPC designation |
| `permid` | Permanent number |
| `epoch_tdb` | Element epoch, `YYYY-MM-DD 00:00:00` TDB (MPCORB epochs are 0h TT; TT - TDB < 2 ms) |
| `a` | Semi-major axis (au) |
| `e` | Eccentricity |
| `q` | Perihelion distance a (1 - e) (au) |
| `incl` | Inclination (deg) |
| `node` | Longitude of ascending node (deg) |
| `peri` | Argument of perihelion (deg) |
| `M` | Mean anomaly at the epoch (deg) |

Also written: `outputs/orbits/orbits.log` — mirrored stdout.

---

### A.5 `outputs/states/states.parquet`

Per-observation geometric state and absolute magnitude. One row per
observation id whose `provid` is in `orbits.csv` (objects without elements get
empty geometry and H). Elements -> heliocentric state at the epoch
(`spiceypy.conics`) -> two-body propagation to the observation time
(`spiceypy.prop2b`) -> light-time corrected geometry, Earth from `de440s`.

| Column | Meaning |
|--------|---------|
| `id` | Observation id (join to mois / obs) |
| `provid_vro` | VRO designation |
| `provid_mpc` | MPC designation |
| `permid` | Permanent number (nullable) |
| `obstime` | Observation epoch (UTC) |
| `obs_et` | Observation time as ephemeris time (TDB s past J2000) |
| `deltat` | Days from the element epoch to the observation (TDB) |
| `r` | Heliocentric distance (au), light-time corrected |
| `delta` | Geocentric distance (au), light-time corrected |
| `alpha` | Phase angle (deg) |
| `mag` | Apparent magnitude (Rubin, any band) |
| `H` | Absolute magnitude, HG system with G = 0.15 and V = `mag` |
| `flag_updated` | `1` if computed in the run that last wrote the file, `0` if kept (orbit unchanged) |
| `orbit_hash` | Hash of `provid_mpc`, `permid`, the epoch and the MPCORB elements (a, e, incl, node, peri, M) |

`load_moi_photometry` takes `id`, `H`, `provid_mpc` from this file.

Also written: `outputs/states/states.log` — mirrored stdout.

---

### A.6 `outputs/mois/mois.parquet`

MOI-filtered observation subset. Same column layout as
`obs_sbn_X05_full.parquet` (~90 ADPACS columns). Key columns for later stages:

| Column | Meaning |
|--------|---------|
| `id` | Observation id |
| `provid` | Designation used as asteroid key downstream |
| `permid` | Permanent number (nullable) |
| `obstime` | Observation epoch |
| `ra`, `dec` | Sky position |
| `mag`, `rmsmag` | Apparent magnitude and error |
| `band` | Filter |

`load_moi_photometry` reads the key set above, then joins states + crossmatches.

---

### A.7 `outputs/periodograms/periodograms.csv`

Democratic period search. **One row per asteroid.**

| Column | Meaning |
|--------|---------|
| `provid` | Asteroid name |
| `epoch_utc_ini` | First observation UTC |
| `epoch_utc_end` | Last observation UTC |
| `obs_duration` | Span in days |
| `n_obs` | Total points used |
| `n_obs_u` … `n_obs_y` | Counts per LSST band |
| `p_rot_dem_h` | Winning democratic rotation period (h) |
| `votes_dem` | Votes for the winner |
| `n_series_dem` | Number of series that cast votes |
| `bands_dem` | Band / combo labels that voted (`Lr+Lg+…`) |
| `periods_dem_h` | Vote tokens `label=period_h` joined by commas (input to foldings) |

Also written under `outputs/periodograms/`:

- `<provid>_periodogram.webp` — periodogram figure
- `periodograms.log` — mirrored stdout

---

### A.8 `outputs/foldings/foldings.csv`

Multiband two-peak folding of democratic candidates. **One row per plotted
solution** (selected and other two-peak panels).

| Column | Meaning |
|--------|---------|
| `provid` | Asteroid name |
| `selected` | `1` if chosen as best folding, else `0` |
| `p_rot_fold_h` | Folded rotation period (h) |
| `chi2_nu` | Reduced chi-square of the harmonic fit |
| `n_peaks` | Number of peaks in the folded curve |
| `power` | Lomb–Scargle / fit power |
| `amp_ratio` | Amplitude ratio diagnostic |
| `amp_data` | Peak-to-peak amplitude of the data |
| `span_alias` | `1` if period is a campaign-span alias, else `0` |
| + A.0 | Color columns |

Also written under `outputs/foldings/`:

- `<provid>_folding.webp` — folding figure
- `foldings.log` — mirrored stdout

---

### A.9 Photometry frame (`mv.load_moi_photometry()`)

In-memory only (not written as a stage file). Built from A.6 + A.5 + A.3.

| Column | Meaning |
|--------|---------|
| `id` | Observation id |
| `obstime` | datetime |
| `permid` | Permanent number (nullable) |
| `provid` | Canonical designation after dual-provid merge |
| `ra`, `dec` | Sky position (deg) |
| `band` | Normalized `Lu`/`Lg`/`Lr`/`Li`/`Lz`/`Ly` |
| `mag` | Apparent magnitude |
| `rmsmag` | Magnitude error (floored at `MIN_MAG_ERROR`) |
| `H` | Absolute magnitude from states |
| `provid_mpc` | MPC designation from states |
| `mpc_designation` | Cleaned MPC string (nullable) |
| `night` | `obstime` floored to calendar day |

---

### A.10 `outputs/mocs/mocs.csv`

Visually retained folding solutions. **One row per asteroid.**

| Column | Meaning |
|--------|---------|
| `provid` | Asteroid name |
| `n_obs` | Observation count at selection time |
| `n_solutions` | How many solutions were kept |
| `solutions` | Slash-separated solution ids, e.g. `1` or `1/2` |
| `periods_h` | Slash-separated periods (h), same order as `solutions` |
| `p_rot_dem_h` | Democratic period from periodograms (context) |
| `selection` | How it was chosen (e.g. `visual_folding`) |

---

### A.11 `outputs/solutions/solutions.csv`

Final folded solutions for retained mocs periods. **One row per solution.**

| Column | Meaning |
|--------|---------|
| `provid` | Asteroid name |
| `solution_id` | Id from `mocs.solutions` |
| `p_rot_h` | Rotation period (h) |
| `chi2_nu` | Reduced chi-square |
| `n_peaks` | Peaks in the folded curve |
| `power` | Fit / periodogram power |
| `amp` | Model amplitude |
| `amp_ratio` | Amplitude ratio diagnostic |
| `amp_data` | Peak-to-peak data amplitude |
| `scatter` | Residual scatter |
| `span_alias` | Campaign-span alias flag |
| `alias_label` | Alias label if applicable (nullable) |
| `filters` | Comma-separated filters used, e.g. `g,r,i` |
| `n_filters` | Number of filters |
| `ref_filter` | Reference filter for colors (usually `r`) |
| `n_obs_u` … `n_obs_y` | Counts per band in the solution |
| + A.0 | Color columns |

Also written under `outputs/solutions/`:

- `<provid>_solution.webp` — stacked folded panels
- `solutions.log` — mirrored stdout

-----
[MVROC Team](https://github.com/seap-udea/MVROC) (2026-present)
