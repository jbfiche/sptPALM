"""Mean square displacement (MSD) of individual tracks.

Port of ``MSD_calculation.m``. For every track and every lag time (in
frames) the squared displacements between all the pairs of detections
separated by that lag are averaged. A lag is kept only if at least
``min_points_msd`` pairs are available, and the MSD of a track stops at the
first lag that does not satisfy this. The weight of each MSD point is
1 / variance of the squared displacements (Saxton, Biophys. J. 1997).

A track is accepted if its MSD has at least ``p`` points, if it moved along
both axes, and if none of its weights is infinite (all the squared
displacements of a lag are identical).

Implementation: the tracks are laid out on a dense frame grid (NaN where
there is no detection, which handles the blinks) and processed by chunks of
tracks of similar duration, so that each lag is one vectorised numpy operation.

Differences with MATLAB (``legacy=True`` restores the MATLAB behaviour):

* hard threshold: a lag is valid with ``>= min_points_msd`` pairs (MATLAB
  required strictly more);
* only the lags needed for the analysis are computed (``max_lag``), whereas
  MATLAB computed all of them (this only matters for the test on infinite
  weights, which is applied on all lags in legacy mode).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .tracks import TrackSet


@dataclass
class MSDResult:
    msd: np.ndarray        # (n_tracks, width) um^2, NaN after the last valid lag
    weight: np.ndarray     # (n_tracks, width)
    n_lags: np.ndarray     # (n_tracks,) number of valid lags (contiguous from lag 1)
    accepted: np.ndarray   # (n_tracks,) bool: track usable for the diffusion analysis

    @property
    def width(self) -> int:
        return self.msd.shape[1]


def compute_msd(ts: TrackSet, width: int, p: int, min_pairs: int = 3,
                legacy: bool = False, chunk_size: int = 256) -> MSDResult:
    """Compute the MSD of every track.

    ``width`` is the number of lags kept in the result (at least ``p`` and the
    number of lags displayed in the MSD plot).
    """
    n_tracks = ts.n_tracks
    width = max(int(width), int(p))
    msd_out = np.full((n_tracks, width), np.nan)
    w_out = np.full((n_tracks, width), np.nan)
    n_lags = np.zeros(n_tracks, dtype=np.int64)
    accepted = np.zeros(n_tracks, dtype=bool)
    if n_tracks == 0:
        return MSDResult(msd_out, w_out, n_lags, accepted)

    # lags to compute: all of them in legacy mode (the Inf-weight test used all lags)
    max_lag = None if legacy else width

    # tracks that moved along both axes
    starts = ts.offsets[:-1]
    ext_x = np.maximum.reduceat(ts.x, starts) - np.minimum.reduceat(ts.x, starts)
    ext_y = np.maximum.reduceat(ts.y, starts) - np.minimum.reduceat(ts.y, starts)
    moved = (ext_x > 0) & (ext_y > 0)

    span = ts.duration_frames
    order = np.argsort(span, kind="stable")
    has_inf = np.zeros(n_tracks, dtype=bool)

    for c0 in range(0, n_tracks, chunk_size):
        idx = order[c0:c0 + chunk_size]
        idx = idx[moved[idx]]
        if idx.size == 0:
            continue
        m = idx.size
        L = int(span[idx].max()) + 1
        gx = np.full((m, L), np.nan)
        gy = np.full((m, L), np.nan)
        for r, i in enumerate(idx):
            a, b = ts.offsets[i], ts.offsets[i + 1]
            col = ts.frame[a:b] - ts.frame[a]
            gx[r, col] = ts.x[a:b]
            gy[r, col] = ts.y[a:b]

        alive = np.ones(m, dtype=bool)
        lag_hi = L - 1 if max_lag is None else min(L - 1, max_lag)
        for lag in range(1, lag_hi + 1):
            d = (gx[:, lag:] - gx[:, :-lag]) ** 2 + (gy[:, lag:] - gy[:, :-lag]) ** 2
            valid = ~np.isnan(d)
            n = valid.sum(axis=1)
            ok = (n > min_pairs) if legacy else (n >= min_pairs)
            alive &= ok
            if not alive.any():
                break
            with np.errstate(invalid="ignore", divide="ignore"):
                mean = np.nansum(d, axis=1) / n
                dev = np.where(valid, d - mean[:, None], 0.0)
                var = (dev ** 2).sum(axis=1) / (n - 1)
                weight = 1.0 / var
            rows = idx[alive]
            if lag <= width:
                msd_out[rows, lag - 1] = mean[alive]
                w_out[rows, lag - 1] = weight[alive]
            n_lags[rows] = lag
            has_inf[rows] |= np.isinf(weight[alive])

    accepted = moved & (n_lags >= p) & ~has_inf
    return MSDResult(msd_out, w_out, n_lags, accepted)
