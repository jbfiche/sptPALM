"""Text outputs of the analysis (same files and columns as the MATLAB version)."""

from __future__ import annotations

import os
from typing import Dict

import numpy as np
import pandas as pd

from ..pipeline import AnalysisResult


def write_diffusion_coefficients(result: AnalysisResult, path: str) -> None:
    """``Saved_Diffusion_Coeff.txt``: log10(D) and D (um^2/s) of every accepted track."""
    df = pd.DataFrame({"logD": result.log_d, "D": result.d_app})
    df.to_csv(path, index=False, float_format="%.15g")


def write_msd(result: AnalysisResult, path: str) -> None:
    """``Saved_MSD.txt``: mean, median and standard error of the MSD per population."""
    if result.fit is None:
        return
    dt = result.params.dt_s
    n = result.fit.msd_summary[0].shape[0]
    cols = {"Time_s": dt * np.arange(1, n + 1)}
    for k, summ in enumerate(result.fit.msd_summary, start=1):
        cols[f"MSD{k}_Average"] = summ[:, 0]
        cols[f"MSD{k}_Median"] = summ[:, 2]
        cols[f"MSD{k}_Error_Bar"] = summ[:, 1]
    pd.DataFrame(cols).to_csv(path, index=False, float_format="%.15g")


def write_parameters(result: AnalysisResult, path: str) -> None:
    """``Parameters_analysis.txt``: parameters and number of tracks at each step."""
    p = result.params
    lines = [
        "Parameters used for the sptPALM analysis (Python version)",
        "",
        f"Analysis                                          : {result.name}",
        f"Legacy (MATLAB) behaviour                         : {p.legacy}",
        f"Acquisition time (ms)                             : {p.acquisition_time_ms:g}",
        f"Pixel size (nm)                                   : {p.pixel_size_um * 1000:g}",
        f"Maximum number of blinks (frames)                 : {p.max_blinks}",
        f"Maximum step length (um)                          : {p.max_step_length_um}",
        f"Minimum number of detections per track            : {p.min_points}",
        f"Minimum duration (frames, 0 = not used)           : {p.min_duration_frames}",
        f"Minimum fraction of frames with a detection       : {p.min_fraction}",
        f"Minimum number of pairs for a MSD point           : {p.min_points_msd}",
        f"Number of MSD points used to calculate D          : {p.msd_fit_points}",
        f"Method for D                                      : {p.diffusion_method}",
        f"Maximum display time for the MSD (s)              : {p.max_display_time_s:g}",
        "",
        f"Initial number of trajectories                    : {result.n_loaded}",
        f"Number of trajectories after filtering            : {result.n_filtered}",
        f"Number of trajectories after ROI selection        : {result.n_roi}",
        f"Number of trajectories validated for the MSD      : {result.n_msd}",
        f"Number of trajectories validated for D            : {result.n_dapp}",
        f"Trajectories in several ROIs (excluded)           : {result.n_ambiguous_roi}",
        f"Density of tracks (/um^2)                         : {result.density:.4g}",
    ]
    if result.fit is not None:
        f = result.fit
        lines += ["", f"Number of Gaussians                               : {f.n_gaussians}",
                  f"R^2 of the fit (%)                                : {100 * f.r_squared:.2f}"]
        for k, (c, s) in enumerate(zip(f.centers, f.sigmas), start=1):
            label = "slow" if k == 1 else "fast"
            if f.n_gaussians == 1:
                label = "single"
            lines.append(f"Population {k} ({label}): log10(D) = {c:.3f}, D = {10 ** c:.4g} um^2/s, "
                         f"sigma = {s:.3f}")
        if f.fast_fraction is not None:
            lines.append(f"Mobile (fast) fraction (%)                        : {100 * f.fast_fraction:.1f}")
    elif result.fit_error:
        lines += ["", f"Fit not performed: {result.fit_error}"]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")


def write_results(results: Dict[str, AnalysisResult], out_dir: str) -> None:
    """Write the text outputs. With several analyses (one per ROI) each one gets a sub-folder."""
    multiple = len(results) > 1
    for name, res in results.items():
        folder = os.path.join(out_dir, _safe(name)) if multiple else out_dir
        os.makedirs(folder, exist_ok=True)
        write_diffusion_coefficients(res, os.path.join(folder, "Saved_Diffusion_Coeff.txt"))
        write_msd(res, os.path.join(folder, "Saved_MSD.txt"))
        write_parameters(res, os.path.join(folder, "Parameters_analysis.txt"))
    if multiple:
        summary_table(results).to_csv(os.path.join(out_dir, "Summary_ROIs.csv"), index=False)


def summary_table(results: Dict[str, AnalysisResult]) -> pd.DataFrame:
    """One row per analysis (ROI): track numbers, density and fit results."""
    rows = []
    for name, r in results.items():
        row = {"analysis": name, "n_loaded": r.n_loaded, "n_filtered": r.n_filtered,
               "n_roi": r.n_roi, "n_msd": r.n_msd, "n_dapp": r.n_dapp,
               "density_per_um2": r.density}
        if r.fit is not None:
            row["n_gaussians"] = r.fit.n_gaussians
            row["logD_slow"] = r.fit.centers[0]
            row["logD_fast"] = r.fit.centers[1] if r.fit.n_gaussians == 2 else np.nan
            row["fast_fraction"] = r.fit.fast_fraction if r.fit.fast_fraction is not None else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_ ." else "_" for c in name).strip() or "roi"
