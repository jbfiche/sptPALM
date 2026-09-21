"""Full analysis, replacing ``Trajectory_analysis_v4.m``.

    tracks, info = load_trackmate_folder(folder, pixel_size_um)
    results = run_analysis(tracks, params)              # no ROI
    results = run_analysis(tracks, params, rois=[...])  # pooled or separate (params.roi_mode)

``run_analysis`` returns a dict of :class:`AnalysisResult`:

* ``{"all": ...}`` without ROI,
* ``{"pooled": ...}`` with ROIs and ``roi_mode="pooled"``,
* ``{roi.name: ...}`` with ROIs and ``roi_mode="separate"``.

Nothing is asked to the user: the ROIs, the number of Gaussians and the split
threshold are arguments.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .diffusion import apparent_diffusion
from .filtering import filter_tracks
from .fitting import DiffusionFit, fit_distribution
from .msd import MSDResult, compute_msd
from .params import AnalysisParams
from .roi import ROI, assign_tracks
from .tracks import TrackSet


@dataclass
class AnalysisResult:
    name: str
    params: AnalysisParams
    n_loaded: int                   # tracks loaded (after rejection of zero-length steps)
    n_filtered: int                 # after the filtering of the blinks / duration
    n_roi: int                      # after the ROI selection
    n_msd: int                      # validated for the MSD calculation
    n_dapp: int                     # validated for the calculation of D
    tracks: TrackSet                # the n_dapp tracks used for the distribution of D
    msd: np.ndarray                 # (n_dapp, width) MSD (um^2) of these tracks
    weight: np.ndarray              # (n_dapp, width) Saxton weights
    d_app: np.ndarray               # (n_dapp,) apparent D (um^2/s)
    fit: Optional[DiffusionFit]     # None if the fit was not possible (see fit_error)
    fit_error: Optional[str] = None
    density: float = float("nan")   # tracks / um^2
    area_um2: float = float("nan")
    n_ambiguous_roi: int = 0        # tracks in several ROIs (excluded)
    roi_names: Tuple[str, ...] = ()

    @property
    def log_d(self) -> np.ndarray:
        return np.log10(self.d_app)

    @property
    def population_labels(self) -> np.ndarray:
        """1 = slow, 2 = fast (only slow for one population)."""
        return self.fit.labels if self.fit is not None else np.ones(self.n_dapp, dtype=int)


def _analyze(name: str, ts_roi: TrackSet, params: AnalysisParams, n_loaded: int, n_filtered: int,
             area_px2: Optional[float], n_gaussians: int, split_threshold: Optional[float],
             n_ambiguous: int, roi_names: Tuple[str, ...]) -> AnalysisResult:
    lags = params.display_lags
    res: MSDResult = compute_msd(ts_roi, width=max(params.msd_fit_points, lags),
                                 p=params.msd_fit_points, min_pairs=params.min_points_msd,
                                 legacy=params.legacy)
    ts_msd = ts_roi.select(res.accepted)
    msd = res.msd[res.accepted]
    weight = res.weight[res.accepted]
    n_msd = ts_msd.n_tracks

    d, ok = apparent_diffusion(msd, weight, params.msd_fit_points, params.dt_s,
                               params.diffusion_method)
    ts_d = ts_msd.select(ok)
    msd_d, w_d, d_app = msd[ok], weight[ok], d[ok]

    # track density
    n_roi = ts_roi.n_tracks
    if area_px2 is not None and area_px2 > 0:
        area_um2 = area_px2 * params.pixel_size_um ** 2
        density = n_roi / area_um2
        if params.legacy:
            density = round(1000 * density) / 1000
    else:
        area_um2, density = float("nan"), float("nan")

    fit, fit_error = None, None
    if len(d_app) >= 10:
        try:
            fit = fit_distribution(np.log10(d_app), msd_d, n_gaussians, lags,
                                   split_threshold=split_threshold, legacy=params.legacy)
        except (ValueError, RuntimeError) as exc:
            fit_error = str(exc)
    else:
        fit_error = f"only {len(d_app)} tracks with a valid D: not enough to fit a distribution"

    return AnalysisResult(name=name, params=params, n_loaded=n_loaded, n_filtered=n_filtered,
                          n_roi=n_roi, n_msd=n_msd, n_dapp=len(d_app), tracks=ts_d, msd=msd_d,
                          weight=w_d, d_app=d_app, fit=fit, fit_error=fit_error, density=density,
                          area_um2=area_um2, n_ambiguous_roi=n_ambiguous, roi_names=roi_names)


def run_analysis(tracks: TrackSet, params: AnalysisParams, rois: Optional[Sequence[ROI]] = None,
                 image_shape: Optional[Tuple[int, int]] = None,
                 roi_settings: Optional[Dict[str, dict]] = None) -> Dict[str, AnalysisResult]:
    """Run the analysis.

    Parameters
    ----------
    tracks : all the tracks (see :func:`sptpalm.io.load_trackmate_folder`).
    params : analysis parameters.
    rois : optional ROIs (see :mod:`sptpalm.roi`).
    image_shape : (rows, columns) of the movies, used for the track density
        when there is no ROI.
    roi_settings : optional per-ROI overrides for ``roi_mode="separate"``:
        ``{"cell 1": {"n_gaussians": 1, "split_threshold": -1.2}}``.
    """
    n_loaded = tracks.n_tracks
    filtered = filter_tracks(tracks, params)
    n_filtered = filtered.n_tracks
    ng, split = params.n_gaussians, params.split_threshold
    roi_settings = roi_settings or {}

    if not rois:
        area = image_shape[0] * image_shape[1] if image_shape else None
        return {"all": _analyze("all", filtered, params, n_loaded, n_filtered, area, ng, split, 0, ())}

    index, n_amb = assign_tracks(filtered, rois, params.pixel_size_um)
    names = tuple(r.name for r in rois)
    round_area = params.legacy

    if params.roi_mode == "pooled":
        ts_roi = filtered.select(index >= 0)
        area = sum(r.area_px2(round_area) for r in rois)
        return {"pooled": _analyze("pooled", ts_roi, params, n_loaded, n_filtered, area, ng, split,
                                   n_amb, names)}

    out: Dict[str, AnalysisResult] = {}
    for k, roi in enumerate(rois):
        opts = roi_settings.get(roi.name, {})
        out[roi.name] = _analyze(roi.name, filtered.select(index == k), params, n_loaded, n_filtered,
                                 roi.area_px2(round_area), opts.get("n_gaussians", ng),
                                 opts.get("split_threshold", split), n_amb, (roi.name,))
    return out
