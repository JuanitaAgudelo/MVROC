# MVROC pipeline driver pattern

How to write a stage driver under `pipeline/`. Reference implementations:
`mvroc.py`, `periodograms.py`, `foldings.py`, `solutions.py`.

When this document and the code disagree, the code wins — then update this file.

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
  inputs/                       # raw downloads (obs tables)
    obs_sbn_X05_full.parquet
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
  periodograms.py
  foldings.py
  solutions.py
  PATTERN.md
```

| Path from `pipeline/` | Constant in `mvroc.py` | Role |
|-----------------------|------------------------|------|
| `inputs/` | `INPUTS` | Raw observation tables |
| `outputs/` | `OUTPUTS` | Per-stage pipeline outputs |

**Do not use** `pathlib.Path.resolve()`, `__file__`, or path helper wrappers
(`set_fig_dir`, `results_path`, …).

---

## 3. Upstream and downstream

The pipeline is **sequential**. A driver assumes earlier stages already ran.

Before writing a driver, declare:

1. **Stage name** — short plural or stage noun (`periodograms`, `foldings`, `solutions`).
2. **Upstream** — stage outputs (or `inputs/`) that must exist first.
3. **Downstream** — stages that will read this stage later.

### Current stage order

| Order | Stage | Writes | Upstream |
|------:|-------|--------|----------|
| 0 | seed | `inputs/` (raw obs downloads) | external downloads |
| 1 | `crossmatches` | `outputs/crossmatches/crossmatches.csv` | `inputs/` + MPC |
| 2 | `orbits` | `outputs/orbits/orbits.csv` | crossmatches |
| 3 | `states` | `outputs/states/states.parquet` | orbits + observations |
| 4 | `mois` | `outputs/mois/mois.parquet` | states (MOI filter) |
| 5 | `periodograms` | `outputs/periodograms/` | mois + states + crossmatches |
| 6 | `foldings` | `outputs/foldings/` | periodograms + photometry |
| 7 | `mocs` | `outputs/mocs/mocs.csv` | foldings (visual / selected periods) |
| 8 | `solutions` | `outputs/solutions/` | mocs + photometry |
| 9 | shapefit | (future) | solutions |

Until early-stage drivers exist, their outputs may be **simulated** under
`outputs/<stage>/` (see current `crossmatches`, `orbits`, `states`, `mois`,
`mocs`).

### Rules

- Read earlier stage results from `outputs/<upstream>/`.
- Read raw observations from `inputs/` (via `mv.OBSERVATIONS_PARQUET` when needed).
- Do **not** reimplement an upstream stage inside a later driver.
- If an upstream file is missing, exit with a clear message.
- Document Inputs / Outputs (including upstream paths) in the driver docstring.
- Split monolithic procedures into thin drivers along this order.

### Example chain

```text
inputs → crossmatches → orbits → states → mois
  → periodograms → foldings → mocs → solutions → shapefit
```

**`foldings.py`** reads:

- `outputs/periodograms/periodograms.csv` (voted periods)
- photometry via `mv.load_moi_photometry()` (mois + states + crossmatches)

and writes only under `outputs/foldings/`.

**`mocs`** (future driver) selects which folded periods are retained after
visual review and writes `outputs/mocs/mocs.csv`.

**`solutions.py`** reads:

- `outputs/mocs/mocs.csv` (retained solution ids and periods)
- photometry via `mv.load_moi_photometry()`

and writes only under `outputs/solutions/`.

---

## 4. Layout and execution

```bash
cd pipeline
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

CROSSMATCHES_CSV = "outputs/crossmatches/crossmatches.csv"
ORBITS_CSV = "outputs/orbits/orbits.csv"
STATES_PARQUET = "outputs/states/states.parquet"
MOIS_PARQUET = "outputs/mois/mois.parquet"
MOCS_CSV = "outputs/mocs/mocs.csv"
```

Photometry for photometric stages is loaded once with:

```python
df = mv.load_moi_photometry()
# joins outputs/mois + outputs/states + outputs/crossmatches
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
- common helpers: `cadence_min_period`, `merge_frequencies`, `already_in_outputs`
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
- No force / fig-dir / dead flags from the monolith.
- If a parameter is almost always the same value, hard-code a constant.
- `qplot=True` is fine when it separates science from plotting.

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
        Retained solution ids and periods.
    outputs/mois/mois.parquet, outputs/states/states.parquet,
    outputs/crossmatches/crossmatches.csv
        MOI photometry with H (via mv.load_moi_photometry).

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

## 13. Porting from `catalogue_pipeline.py`

1. Locate the monolith function and the private helpers only it uses.
2. List real data dependencies (columns, upstream outputs, inputs).
3. Extract the minimum science. Drop:
   - complex argparse / force / skip-by-figure-existence
   - incidental observations CSV writes
   - absolute Path resolve
   - dead or nearly unused parameters
   - neighboring stage code
4. Rewrite in this pattern.
5. Smoke-test:

```bash
cd pipeline
python <stage>.py "1981 US22"
```

6. Check: CSV append, figure (if any), log, and skip on re-run.

---

## 14. Checklist before merge

- [ ] Lives in `pipeline/<stage>.py`
- [ ] Runs from `pipeline/`
- [ ] Upstream / Downstream declared
- [ ] Reads `outputs/<upstream>/` and/or `inputs/` as documented
- [ ] Writes only `outputs/<stage>/`
- [ ] Entry point `mvroc_<stage>()`; stage helpers in the same file
- [ ] CLI: `provid | all`
- [ ] CSV append + `mv.already_in_outputs`
- [ ] Log via `mv.LogOutput`
- [ ] Full Parameters / Returns / Examples docstrings
- [ ] Important steps marked with `# -->`
- [ ] One blank line between routines; no blank after section headers
- [ ] No dead options, no Path resolve
- [ ] Smoke-test OK

---

## 15. Agent request template

```text
Convert procedure <NAME_IN_catalogue_pipeline.py> into a new pipeline driver
following pipeline/PATTERN.md and the style of pipeline/periodograms.py +
pipeline/mvroc.py.

Stage name: <stage>
Upstream: <earlier outputs / inputs>
Downstream: <who consumes this next>

Scope: only <short list of functions / steps>.
Do not port <neighboring steps>.
Smoke-test with provid "1981 US22".
```
