"""Apparent diffusion coefficient of individual tracks.

Port of ``Diff_calculation.m``. Both methods use the first ``p`` points of the
MSD and the relation MSD = 4 D t (2D diffusion):

``average``
    D = (MSD(p) - MSD(1)) / (4 (t_p - t_1)); a fast estimate that ignores
    the localisation error (the offset).

``fit``
    weighted linear fit MSD = a t + b of the first ``p`` points, with the
    Saxton weights and the constraint b <= min(MSD); D = a / 4. Because the
    model is linear, the (constrained) solution is computed in closed form
    for all the tracks at once instead of calling a fit routine per track.

A track is accepted when its ``p`` first MSD points are all positive and D > 0.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np


def dapp_average(msd: np.ndarray, p: int, dt_s: float) -> np.ndarray:
    """Apparent D with the "average" method (NaN when not defined)."""
    m = msd[:, :p]
    lag = dt_s * (p - 1)
    return 0.25 * (m[:, p - 1] - m[:, 0]) / lag


def dapp_fit(msd: np.ndarray, weight: np.ndarray, p: int, dt_s: float) -> np.ndarray:
    """Apparent D with the weighted fit method (closed form, vectorised).

    Minimises sum_i w_i (msd_i - a t_i - b)^2 with b <= min(msd).
    """
    t = dt_s * np.arange(1, p + 1)
    y = msd[:, :p]
    w = weight[:, :p]
    sw = w.sum(axis=1)
    swt = (w * t).sum(axis=1)
    swtt = (w * t * t).sum(axis=1)
    swy = (w * y).sum(axis=1)
    swty = (w * t * y).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        det = sw * swtt - swt ** 2
        a = (sw * swty - swt * swy) / det
        b = (swy - a * swt) / sw
        # constraint violated: the optimum lies on the boundary b = min(msd)
        bmax = y.min(axis=1)
        over = b > bmax
        a_c = (swty - bmax * swt) / swtt
        a = np.where(over, a_c, a)
    return a / 4.0


def apparent_diffusion(msd: np.ndarray, weight: np.ndarray, p: int, dt_s: float,
                       method: str = "average") -> Tuple[np.ndarray, np.ndarray]:
    """Apparent D for every track.

    Returns ``(D, ok)``: ``D`` in um^2/s and the mask of the tracks that give
    a valid (positive) coefficient. ``D`` is NaN where not valid.
    """
    n = msd.shape[0]
    if n == 0:
        return np.empty(0), np.zeros(0, dtype=bool)
    first = msd[:, :p]
    valid_msd = np.all(np.isfinite(first) & (first > 0), axis=1)
    if method == "fit":
        d = dapp_fit(msd, weight, p, dt_s)
    elif method == "average":
        d = dapp_average(msd, p, dt_s)
    else:
        raise ValueError("method must be 'average' or 'fit'")
    ok = valid_msd & np.isfinite(d) & (d > 0)
    return np.where(ok, d, np.nan), ok
