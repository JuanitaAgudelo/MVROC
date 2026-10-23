#!/usr/bin/env python
"""Run catalogue_pipeline.py with one asteroid per process.

The science code stays in catalogue_pipeline.py. This file only distributes
the same main loop across workers. Each worker keeps a single BLAS thread, so
the processes do not compete for the same cores.

    python catalogue_pipeline_parallel.py --workers 8
    python catalogue_pipeline_parallel.py --workers 4 --provid "2016 GC22"
    python catalogue_pipeline_parallel.py --force --workers 8

--workers defaults to one less than the number of CPUs. --force and --provid
mean the same thing as in catalogue_pipeline.py. The asteroid list is the call
in choose_selection(), kept in step with catalogue_pipeline._main.
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import io
import multiprocessing
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd

import catalogue_pipeline as cp

_WORKER = {}
_FLUSH_EVERY = 20


def choose_selection(cli_provids):
    """Same asteroid list as catalogue_pipeline._main."""
    if cli_provids:
        return cp.unpack_selection(cp.select_provid_list(provids=cli_provids))
    return cp.unpack_selection(
        # cp.select_provid_rotation()
        # cp.select_provid_rotfold()
        # cp.select_provid_list()
        cp.select_provid_medellin()
        # cp.select_provid_nobs()
    )


def _parse_args(argv):
    """--workers N, plus the --force and --provid options of the serial script."""
    force = False
    workers = None
    provids = []
    args = list(argv[1:])
    i = 0
    while i < len(args):
        arg = args[i]
        if arg in ("--force", "-f"):
            force = True
            i += 1
            continue
        if arg == "--workers":
            if i + 1 >= len(args):
                raise SystemExit("Missing worker count after --workers.")
            workers = int(args[i + 1])
            i += 2
            continue
        if arg.startswith("--workers="):
            workers = int(arg.split("=", 1)[1])
            i += 1
            continue
        if arg == "--provid":
            if i + 1 >= len(args):
                raise SystemExit("Missing asteroid name after --provid.")
            provids.append(args[i + 1].replace("_", " ").strip())
            i += 2
            continue
        if arg.startswith("--provid="):
            name = arg.split("=", 1)[1].replace("_", " ").strip()
            if not name:
                raise SystemExit("Missing asteroid name after --provid.")
            provids.append(name)
            i += 1
            continue
        if arg.startswith("-"):
            raise SystemExit(f"Unknown option: {arg}")
        provids.append(arg.replace("_", " ").strip())
        i += 1
    if workers is not None and workers < 1:
        raise SystemExit("--workers must be at least 1.")
    if workers is None:
        workers = max(1, (os.cpu_count() or 2) - 1)
    return force, workers, [name for name in provids if name]


def _paint_span_lines(text):
    """Paint baseline-alias lines red when the worker captured them uncolored."""
    painted = []
    for line in text.splitlines(keepends=True):
        bare = line.rstrip("\r\n")
        ending = line[len(bare):]
        if "baseline alias" in bare and "\033[31m" not in bare:
            painted.append(f"\033[31m{bare}\033[0m{ending}")
        else:
            painted.append(line)
    return "".join(painted)


def _init_worker(fig_dir, designations, df_rot, analysis, exclude_lx, nb, force):
    """Point this process at the figure directory without reopening the log."""
    os.environ["FORCE_COLOR"] = "1"
    cp.FIG_DIR = cp.Path(fig_dir)
    _WORKER.clear()
    _WORKER.update(
        fig_dir=cp.FIG_DIR,
        designations=designations,
        df_rot=df_rot,
        analysis=analysis,
        exclude_lx=exclude_lx,
        nb=nb,
        force=force,
    )


def _result_row(df_ast, provid, band, fit, analysis):
    """One catalogue row. Same columns as the serial loop."""
    epoch_utc_ini, epoch_utc_end, obs_duration = cp.observation_window(df_ast)
    n_obs_bands = {
        f"n_obs_{b}": int((df_ast["band"].isin([f"L{b}", b])).sum())
        for b in ["u", "g", "r", "i", "z", "y"]
    }
    row = {
        "provid": provid,
        "epoch_utc_ini": epoch_utc_ini,
        "epoch_utc_end": epoch_utc_end,
        "obs_duration": obs_duration,
        "n_obs": int(len(df_ast)),
        **n_obs_bands,
    }
    if "periodogram" in analysis:
        periods_h = (
            ",".join(
                f"{period:.6f}"
                for period in fit.get("df_peaks", {}).get("rotational_period_h", [])
            )
            if "df_peaks" in fit
            else ""
        )
        row["band_periodogram"] = fit.get("band", band)
        row["periods_h"] = periods_h
        row["power"] = fit.get("power", np.nan)
        row["amp_lsm"] = fit.get("amp_lsm", np.nan)
        row["amp_photo"] = fit.get("amp_photo", np.nan)
    if "periodogram_democratic" in analysis and fit.get("democratic"):
        dem = fit["democratic"]
        row["p_rot_dem_h"] = dem["p_rot_h"]
        row["votes_dem"] = dem["votes"]
        row["n_series_dem"] = dem["n_series"]
        row["bands_dem"] = dem["bands"]
        row["periods_dem_h"] = dem["periods_h"]
    if "periodogram_democratic_weighted" in analysis and fit.get("democratic_weighted"):
        wdem = fit["democratic_weighted"]
        row["p_rot_wdem_h"] = wdem["p_rot_h"]
        row["k_alias_wdem"] = wdem["k_alias"]
        row["score_wdem"] = wdem["score"]
        row["band_wdem"] = wdem["band_ref"]
        row["scores_wdem"] = wdem["scores"]
    if "phase_model_compare" in analysis and fit.get("phase_model_compare"):
        compared = fit["phase_model_compare"]
        row["p_rot_cmp_h"] = compared["p_rot_h"]
        row["chi2_nu_cmp"] = compared["chi2_nu"]
        row["periods_cmp"] = compared["periods"]
    if "multiband_photometry_compare" in analysis and fit.get(
        "multiband_photometry_compare"
    ):
        mb_compared = fit["multiband_photometry_compare"]
        row["p_rot_h"] = mb_compared["p_rot_h"]
        row["chi2_nu"] = mb_compared["chi2_nu"]
        row["n_peaks"] = mb_compared["n_peaks"]
        row["power"] = mb_compared["power"]
        row["amp_ratio"] = mb_compared["amp_ratio"]
        row["amp_data"] = mb_compared["amp_data"]
        row["span_alias"] = mb_compared["span_alias"]
        # Every pair among the filters with a zero point, bluer minus redder.
        row.update(mb_compared.get("colors", {}))
    if "justfold" in analysis and fit.get("justfold"):
        folded = fit["justfold"]
        row["p_rot_justfold_h"] = folded["p_rot_h"]
        row["chi2_nu_justfold"] = folded["chi2_nu"]
        row["n_peaks_justfold"] = folded["n_peaks"]
        row["power_justfold"] = folded["power"]
        row["amp_ratio_justfold"] = folded["amp_ratio"]
        row["span_alias_justfold"] = folded["span_alias"]
        row["periods_justfold"] = folded["periods"]
    if "shape_fit" in analysis:
        row.update(fit.get("shape", {}))
    if "multiband_photometry" in analysis:
        if "multiband" in fit and fit["multiband"] is not None:
            summary = dict(fit["multiband"].get("summary_row", {}))
            if "multiband_photometry_compare" in analysis:
                for key in cp.COMPARE_COLOR_COLUMNS:
                    summary.pop(key, None)
            row.update(summary)
    if "shape_fit_multiband" in analysis and "shape_k" in fit:
        row.update(fit["shape_k"])
    if "spectrum" in analysis:
        row.update(fit.get("spectrum") or {})
    return row


def _run_asteroid(task):
    """Process one asteroid. The log text goes back to the parent in one piece."""
    provid, band, df_ast = task
    buffer = io.StringIO()
    previous = sys.stdout
    sys.stdout = buffer
    status = "fail"
    err = None
    row = None
    try:
        cp.FIG_DIR = _WORKER["fig_dir"]
        fit, err = cp.process(
            df_ast,
            provid,
            _WORKER["df_rot"],
            _WORKER["analysis"],
            _WORKER["designations"],
            force=_WORKER["force"],
            band=band,
            exclude_lx=_WORKER["exclude_lx"],
            Nb=_WORKER["nb"],
        )
        if fit == "skip":
            status = "skip"
        elif fit is not None:
            status = "ok"
            row = _result_row(df_ast, provid, band, fit, _WORKER["analysis"])
        else:
            status = "fail"
            err = err or "Not enough data, or an error in the calculation"
    except Exception:
        status = "fail"
        err = traceback.format_exc()
        print(err)
    finally:
        sys.stdout = previous
    return {
        "provid": provid,
        "status": status,
        "err": err,
        "row": row,
        "log": buffer.getvalue(),
    }


def _drop_failed_rows(replaced):
    """Remove asteroids that failed this run from an existing results table."""
    path = cp.FIG_DIR / "asteroid_analysis_results.csv"
    if not replaced or not path.exists():
        return
    old = pd.read_csv(path)
    kept = old[~old["provid"].isin(replaced)]
    if len(kept) == len(old):
        return
    temporary = path.with_suffix(".csv.tmp")
    kept.to_csv(temporary, index=False)
    os.replace(temporary, path)


def _flush_results(rows, replaced, analysis):
    """Write finished asteroids. The parent is the only process that does this."""
    if rows:
        cp.write_analysis_results(rows, replaced, analysis)
        return
    _drop_failed_rows(replaced)


def _consume(results, analysis, succeeded, failed, pending_rows, pending_replaced):
    """Record one finished asteroid and flush the table every few of them."""
    status = results["status"]
    provid = results["provid"]
    if results["log"]:
        text = _paint_span_lines(results["log"])
        print(text, end="" if text.endswith("\n") else "\n")
    if status == "skip":
        succeeded.append(provid)
    elif status == "ok":
        succeeded.append(provid)
        pending_rows.append(results["row"])
        pending_replaced.add(provid)
    else:
        failed[provid] = results["err"] or "Not enough data, or an error in the calculation"
        pending_replaced.add(provid)
    if len(pending_replaced) >= _FLUSH_EVERY:
        _flush_results(pending_rows, pending_replaced, analysis)
        pending_rows.clear()
        pending_replaced.clear()


def _progress(done, n_asteroids, done_obs, total_obs, t_start, t_loop):
    elapsed = time.perf_counter() - t_start
    estimated_total = (t_loop - t_start) + (time.perf_counter() - t_loop) * (
        total_obs / done_obs
    )
    remaining = max(0.0, estimated_total - elapsed)
    print(
        f"  elapsed: {cp.format_elapsed(elapsed)} ({done}/{n_asteroids}), "
        f"estimated total: {cp.format_elapsed(estimated_total)}, "
        f"time remaining (est.): {cp.format_elapsed(remaining)}"
    )


def _finish(pending_rows, pending_replaced, analysis, df_rot, total_requested, succeeded, failed, t_start):
    _flush_results(pending_rows, pending_replaced, analysis)
    cp.copy_analysis_docs()
    if "compare" in analysis:
        print("Period comparison with the literature")
        cp.plot_period_comparison(df_rot)
    cp.print_processing_report(
        total_requested,
        succeeded,
        failed,
        fig_dir=cp.FIG_DIR,
        elapsed_s=time.perf_counter() - t_start,
    )


def main():
    t_start = time.perf_counter()
    force, workers, cli_provids = _parse_args(sys.argv)
    try:
        selected, df_set, df_rot, analysis, missing = choose_selection(cli_provids)
        if not force:
            selected = cp.exclude_finished(selected)
        designations = cp.load_mpc_designations()
        total_requested = len(selected) + len(missing)
        succeeded = []
        failed = {
            provid: "No photometric observations in obs_sbn_X05_full_mpc_desig.parquet"
            for provid in missing
        }
        exclude_lx = tuple(selected.attrs.get("exclude_lx", cp.EXCLUDE_LX))
        nb = int(selected.attrs.get("Nb", 2))
        bands = (
            selected["band"]
            if "band" in selected.columns
            else pd.Series(cp.BAND, index=selected.index)
        )
        n_asteroids = len(selected)
        total_obs = int(selected["n_obs"].sum()) if n_asteroids else 0
        print(
            f"Parallel run: {workers} worker{'s' if workers != 1 else ''}, "
            f"{n_asteroids} asteroids"
        )
        grouped = df_set.groupby("provid", sort=False)
        obs_of = {
            provid: int(n_obs)
            for provid, n_obs in zip(selected["provid"], selected["n_obs"])
        }

        def iter_tasks():
            for provid, band in zip(selected["provid"], bands):
                if provid in grouped.groups:
                    part = grouped.get_group(provid)
                else:
                    part = df_set.iloc[0:0]
                yield provid, str(band), part

        pending_rows = []
        pending_replaced = set()
        t_loop = time.perf_counter()
        done_obs = 0
        initargs = (
            str(cp.FIG_DIR),
            designations,
            df_rot,
            analysis,
            exclude_lx,
            nb,
            force,
        )
        if n_asteroids == 0:
            pass
        elif workers == 1:
            _init_worker(*initargs)
            for done, task in enumerate(iter_tasks(), start=1):
                result = _run_asteroid(task)
                _consume(
                    result, analysis, succeeded, failed, pending_rows, pending_replaced
                )
                done_obs += obs_of.get(result["provid"], 0)
                _progress(done, n_asteroids, done_obs, total_obs, t_start, t_loop)
        else:
            context = multiprocessing.get_context("spawn")
            with context.Pool(
                workers,
                initializer=_init_worker,
                initargs=initargs,
                maxtasksperchild=40,
            ) as pool:
                for done, result in enumerate(
                    pool.imap_unordered(_run_asteroid, iter_tasks(), chunksize=1),
                    start=1,
                ):
                    _consume(
                        result,
                        analysis,
                        succeeded,
                        failed,
                        pending_rows,
                        pending_replaced,
                    )
                    done_obs += obs_of.get(result["provid"], 0)
                    _progress(done, n_asteroids, done_obs, total_obs, t_start, t_loop)
        _finish(
            pending_rows,
            pending_replaced,
            analysis,
            df_rot,
            total_requested,
            succeeded,
            failed,
            t_start,
        )
    finally:
        cp._stop_log()


if __name__ == "__main__":
    main()
