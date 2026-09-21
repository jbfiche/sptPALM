"""Figures of the analysis (matplotlib).

Every function draws on an ``Axes`` (created if not given) so they can be
used in scripts, notebooks and, later, embedded in the GUI.
"""

from __future__ import annotations

import os
from typing import Dict, Optional

import numpy as np

from .pipeline import AnalysisResult
from .tracks import TrackSet

SLOW_COLOR = "#0072BD"   # blue
FAST_COLOR = "#D95319"   # orange


def _ax(ax):
    if ax is None:
        import matplotlib.pyplot as plt
        _, ax = plt.subplots(figsize=(5, 5))
    return ax


def plot_step_length_ecdf(tracks: TrackSet, ax=None):
    """Cumulative distribution of the step lengths, to check that the tracking
    parameters do not artificially shorten the tracks."""
    ax = _ax(ax)
    steps = np.sort(tracks.step_lengths())
    n = len(steps)
    f = np.arange(1, n + 1) / n

    def cut(q):
        i = np.argmax(f > q)
        return steps[i]

    ax.plot(steps, f, "-", color="r", lw=1, label=f"All values, Lmax = {steps[-1]:.2f} um")
    for q, color, lab in ((0.995, (1, 0.5, 0), "0.5%"), (0.99, "g", "1%"),
                          (0.985, (0, 0.5, 1), "1.5%"), (0.98, "b", "2%")):
        sub = steps[steps < cut(q)]
        ax.plot(sub, np.arange(1, len(sub) + 1) / len(sub), "-", color=color, lw=1,
                label=f"All values without the {lab} longest, Lmax = {sub[-1]:.2f} um")
    ax.set_xlim(0, steps[steps < cut(0.995)][-1])
    ax.set_ylim(0, 1)
    ax.set_xlabel("Step length (um)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title("Cumulative distribution of the step length")
    ax.legend(loc="lower right", fontsize=7)
    ax.set_box_aspect(1)
    return ax


def plot_track_durations(tracks: TrackSet, acquisition_time_ms: float, ax=None):
    """Distribution of the track durations (binned up to the 99th percentile)."""
    ax = _ax(ax)
    dt = acquisition_time_ms / 1000.0
    length = np.sort(tracks.duration_frames * dt)
    max99 = length[max(int(round(len(length) * 0.99)) - 1, 0)]
    bins = np.arange(0, max99 + dt, dt)
    counts, edges = np.histogram(length, bins=bins)
    counts = counts / counts.sum() if counts.sum() else counts
    centers = 0.5 * (edges[1:] + edges[:-1])
    ax.bar(centers, counts, width=dt, color="0.7")
    cum = np.cumsum(counts)
    if counts.sum():
        t80 = centers[np.argmax(cum > 0.8)]
        t90 = centers[np.argmax(cum > 0.9)]
        ax.axvline(t80, ls="--", color="g", lw=0.8, label="80% limit")
        ax.axvline(t90, ls="--", color="r", lw=0.8, label="90% limit")
        ax.legend(loc="upper right")
    ax.set_xlim(0, max99)
    ax.set_xlabel("Trajectories duration (s)")
    ax.set_ylabel("Fraction of trajectories")
    ax.set_title("Trajectories duration distribution")
    ax.set_box_aspect(1)
    return ax


def plot_diffusion_distribution(result: AnalysisResult, ax=None):
    """Histogram of log10(D) with the Gaussian fit(s)."""
    ax = _ax(ax)
    fit = result.fit
    if fit is None:
        ax.hist(result.log_d, bins=30, color="0.6")
        ax.set_title(result.fit_error or "")
    else:
        width = fit.hist_bin[1] - fit.hist_bin[0]
        if fit.n_gaussians == 2:
            ax.bar(fit.hist_bin, fit.hist_values_pop[:, 0], width=width, color=SLOW_COLOR)
            ax.bar(fit.hist_bin, fit.hist_values_pop[:, 1], width=width, color=FAST_COLOR,
                   bottom=fit.hist_values_pop[:, 0])
            ax.plot(fit.curve_x, fit.curve_pop[:, 0], "-b", lw=1)
            ax.plot(fit.curve_x, fit.curve_pop[:, 1], "-r", lw=1)
            ax.plot(fit.curve_x, fit.curve_total, "--k", lw=1)
        else:
            ax.bar(fit.hist_bin, fit.hist_values, width=width, color="0.6")
            ax.plot(fit.curve_x, fit.curve_total, "--k", lw=1)
        lines = [f"log10(D{k}) = {c:.2f} -- D{k} = {10 ** c:.3f} um^2/s"
                 for k, c in enumerate(fit.centers, start=1)]
        if fit.fast_fraction is not None:
            lines.append(f"Mobile fraction : {fit.mobile_fraction_percent:.0f}%")
        ax.set_title("\n".join(lines), fontsize=9)
    ax.set_xlabel("Log10 of the apparent diffusion coefficient (um^2/s)")
    ax.set_ylabel("Fraction of molecule (%)")
    ax.set_box_aspect(1)
    return ax


def plot_msd(result: AnalysisResult, ax=None):
    """Median MSD (with standard error) of each population."""
    ax = _ax(ax)
    if result.fit is None:
        return ax
    t = result.params.dt_s * np.arange(1, result.fit.msd_summary[0].shape[0] + 1)
    colors = [SLOW_COLOR, FAST_COLOR] if result.fit.n_gaussians == 2 else [SLOW_COLOR]
    for k, (summ, color) in enumerate(zip(result.fit.msd_summary, colors), start=1):
        ax.errorbar(t, summ[:, 2], yerr=summ[:, 1], fmt="-s", color=color, lw=2,
                    label=f"MSD{k} median")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("MSD (um^2)")
    ax.legend(loc="upper left")
    ax.set_box_aspect(1)
    return ax


def plot_trajectories(result: AnalysisResult, ax=None, by_population: bool = True):
    """Trajectories used for the analysis (um), coloured by population."""
    ax = _ax(ax)
    labels = result.population_labels
    for i in range(result.tracks.n_tracks):
        _, x, y = result.tracks.get(i)
        color = "0.3"
        if by_population and result.fit is not None and result.fit.n_gaussians == 2:
            color = SLOW_COLOR if labels[i] == 1 else FAST_COLOR
        ax.plot(x, y, "-", color=color, lw=0.6)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    return ax


def save_figures(results: Dict[str, AnalysisResult], out_dir: str,
                 loaded_tracks: Optional[TrackSet] = None, trajectories: bool = False) -> None:
    """Save the standard figures as png (same names as the MATLAB version)."""
    import matplotlib
    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    first = next(iter(results.values()))
    os.makedirs(out_dir, exist_ok=True)
    if loaded_tracks is not None and loaded_tracks.n_tracks:
        fig, ax = plt.subplots(figsize=(6, 6))
        plot_step_length_ecdf(loaded_tracks, ax)
        fig.savefig(os.path.join(out_dir, "Cumulative_Distribution_LengthStep.png"), dpi=150)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(6, 6))
        plot_track_durations(loaded_tracks, first.params.acquisition_time_ms, ax)
        fig.savefig(os.path.join(out_dir, "Trajectories_duration.png"), dpi=150)
        plt.close(fig)

    multiple = len(results) > 1
    for name, res in results.items():
        folder = os.path.join(out_dir, "".join(c if c.isalnum() or c in "-_ ." else "_" for c in name)) \
            if multiple else out_dir
        os.makedirs(folder, exist_ok=True)
        method = "Fit_Method" if res.params.diffusion_method == "fit" else "Weighted_Average_Method"
        fig, ax = plt.subplots(figsize=(6, 6))
        plot_diffusion_distribution(res, ax)
        fig.tight_layout()
        fig.savefig(os.path.join(folder, f"Diffusion_distribution_{method}.png"), dpi=150)
        plt.close(fig)
        fig, ax = plt.subplots(figsize=(6, 6))
        plot_msd(res, ax)
        fig.savefig(os.path.join(folder, "MSD_Curves.png"), dpi=150)
        plt.close(fig)
        if trajectories:
            fig, ax = plt.subplots(figsize=(8, 8))
            plot_trajectories(res, ax)
            fig.savefig(os.path.join(folder, "Trajectories.png"), dpi=150)
            plt.close(fig)
