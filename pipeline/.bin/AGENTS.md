# Agent notes: MVROC pipeline drivers

Living document for agents working under `pipeline/`. When conventions
diverge, prefer the code and then update this file:

- `pipeline/PATTERN.md` — layout, stages, data download, column schemas
- `pipeline/mvroc.py` — shared helpers
- `pipeline/crossmatches.py`, `orbits.py`, `states.py`, `periodograms.py`,
  `foldings.py`, `solutions.py` — reference drivers

---

## Role

Port or write scientific stage drivers under `pipeline/`. Do not invent a new
architecture. Copy the human style of the reference drivers and follow
`PATTERN.md`.

---

## Goal of each driver

Produce:

1. An executable `pipeline/<stage>.py` (cwd = `pipeline/`).
2. Stage-local helpers in that same file.
3. New routines in `mvroc.py` only if a **second** stage needs them.
4. Outputs only under `outputs/<stage>/` (CSV, WebP figures, log).

---

## Upstream / Downstream (required)

The pipeline is **sequential**. Each driver assumes earlier stage outputs
already exist.

Before writing or porting, declare:

1. **Stage name** (e.g. `periodograms`).
2. **Upstream** — exact paths under `inputs/` or `outputs/<upstream>/`.
3. **Downstream** — which stages will read this next.

Rules:

- Raw downloads live in `inputs/` (never absolute monolith paths).
- Earlier results live in `outputs/<upstream>/`.
- Photometry for photometric stages = `mv.load_moi_photometry()` (join of
  `mois` + `states` + `crossmatches`), not a file named `photometry.*`.
- Do not reimplement an upstream stage; exit clearly if a required file is
  missing.
- Document Inputs / Outputs (including upstream paths) in the driver
  docstring.

Stage order (see `PATTERN.md` section 3):

```text
inputs → crossmatches → orbits → states → mois
  → periodograms → foldings → mocs → solutions → shapefit
```

If a historical procedure mixes several steps, split into thin drivers along
that order.

---

## Fresh clone / data

From `pipeline/`:

```bash
make install
make get_gdrive_data
```

Optional (only when refreshing public obs / MPCORB): `make get_public_data`.
Without `make`: `python -m pip install -r .bin/requirements.txt`,
`python .bin/get_gdrive_data.py`, `python .bin/get_public_data.py`.

SPICE kernels for `states.py` are not downloaded: copy `naif0012.tls`,
`gm_de440.tpc` and `de440s.bsp` from `code/kernels/` to `inputs/kernels/`.

Details: `PATTERN.md` (**Getting the data**) and `QUICKSTART.md`.

---

## Layout and execution

- Always cwd = `pipeline/`.
- Relative paths only: `inputs/...`, `outputs/...`.
- No `pathlib.Path.resolve()`, `__file__`, or path wrappers
  (`set_fig_dir`, `results_path`, …).
- Create the output directory once at driver start:
  `os.makedirs(OUT_DIR, exist_ok=True)`.
- Do not re-check paths inside the per-asteroid loop.

Typical constants at the top of the driver:

```python
OUT_DIR = "outputs/<stage>"
RESULTS_CSV = f"{OUT_DIR}/<stage>.csv"
LOG_FILE = f"{OUT_DIR}/<stage>.log"
RESULTS_COLUMNS = [...]
```

Shared path constants live in `mvroc.py` (`INPUTS`, `OUTPUTS`,
`OBSERVATIONS_PARQUET`, `MPCORB_DAT`, `KERNELS`, `MOIS_PARQUET`,
`STATES_PARQUET`, `CROSSMATCHES_CSV`, `ORBITS_CSV`, `MOCS_CSV`, …).

---

## Naming

- Driver file: `pipeline/<stage>.py`
- Entry point: `mvroc_<stage>()`
- Per-asteroid science: `mvroc_<algo>_asteroid(...)`
- Figure: `mvroc_<algo>_plot(...)`
- Skip helper: `mv.already_in_outputs(provid, RESULTS_CSV)`; for many objects
  or a `provid_vro` key: `mv.provids_in_outputs(RESULTS_CSV, column)`
- Single-object printout: `print_<algo>(...)` (`print_crossmatch`,
  `print_orbit`, `print_state`)
- Arrays: plural `s` suffix (`ts`, `ys`, `dys`)
- Constants: `UPPER_CASE`
- Prefer the literal `"H"` over a useless global `PROP = "H"`

Figure stems:

| Stage | Path |
|-------|------|
| periodograms | `outputs/periodograms/<provid>_periodogram.webp` |
| foldings | `outputs/foldings/<provid>_folding.webp` |
| solutions | `outputs/solutions/<provid>_solution.webp` |

---

## `mvroc.py` vs the driver

**`mvroc.py` only** what more than one stage needs:

- path constants
- `load_moi_photometry`, `prepare_asteroid`, `observation_counts`, …
- `LogOutput`, `save_figure`, `slug`
- `cadence_min_period`, `merge_frequencies`, `already_in_outputs`,
  `provids_in_outputs`
- shared constants (`MIN_MAG_ERROR`, `DPI`, …)

**Keep in the driver:** stage science, constants, plot, CSV / log / CLI.

If unsure, leave it in the driver.

---

## CLI

One positional argument only:

```text
python <stage>.py "<provid>" | python <stage>.py <provid_with_underscores> | python <stage>.py all
```

- Underscores → spaces (`1981_US22` → `1981 US22`), except `all`.
- No argparse; no force / fig-dir flags.
- Hard-code near-constant parameters.
- `qplot=True` is fine when it separates science from plotting.
- With a single provid, catalogue-wide stages (`crossmatches`, `orbits`,
  `states`) print a per-object summary instead of the catalogue summary.

---

## Results CSV

- Path: `outputs/<stage>/<stage>.csv`
- Create with header if missing; **append** per asteroid / solution
- Never rewrite or upsert; user deletes the CSV to start clean
- Skip with `mv.already_in_outputs(provid, RESULTS_CSV)` (or
  `mv.provids_in_outputs` when the key is `provid_vro`)
- Do not write incidental `*_observations.csv` files
- Exception: `states.parquet` cannot be appended; `states.py` rewrites it
  every run, keeping rows of objects whose `orbit_hash` (provid_mpc, permid,
  epoch, elements) is unchanged (`flag_updated = 0`) and recomputing the rest
  (`flag_updated = 1`)
- Schema guard: call `mv.check_results_header(RESULTS_CSV, RESULTS_COLUMNS)`
  at start so rows are never appended under an old header

---

## Logging and comments

- Mirror stdout with `mv.LogOutput` to `outputs/<stage>/<stage>.log`
- Detail lines indented with two spaces
- Milestones: `# --> ...`
- ASCII only in messages and docstrings
- One blank line between routines; no blank line after section headers

---

## Figures

- `mv.save_figure(fig, path)` (WebP + watermark)
- Keep calculation and plot in separate routines

---

## Docstrings

NumPy / human style. Full `Parameters` / `Returns` / `Examples` on helpers
and science routines. No double backticks. No `# doctest: +SKIP`.

The entry point must document Stage / Upstream / Downstream, command-line,
Examples, Inputs, Outputs, Returns (see `PATTERN.md` section 12).

---

## Porting an existing procedure

1. Locate the source function and the private helpers only it uses.
2. List real data dependencies (columns, upstream outputs, `inputs/`).
3. Extract the minimum science. Drop argparse / force / Path resolve /
   dead parameters / neighboring stage code / incidental CSV writes.
4. Rewrite in this pattern.
5. Smoke-test:

```bash
cd pipeline
python <stage>.py "1981 US22"
```

6. Check: CSV append, figure (if any), log, skip on re-run.

---

## Checklist before done

- [ ] Lives in `pipeline/<stage>.py`
- [ ] Runs from `pipeline/`
- [ ] Upstream / Downstream declared
- [ ] Reads `outputs/<upstream>/` and/or `inputs/` as documented
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

## Request template

```text
Convert procedure <NAME> into a new pipeline driver following
pipeline/PATTERN.md, pipeline/.bin/AGENTS.md, and the style of
pipeline/periodograms.py + pipeline/mvroc.py.

Stage name: <stage>
Upstream: <exact paths from PATTERN.md section 3>
Downstream: <who consumes this next>

Scope: only <short list of functions / steps>.
Do not port <neighboring steps>.
Smoke-test with provid "1981 US22".
```

-----
[MVROC Team](https://github.com/seap-udea/MVROC) (2026-present)
