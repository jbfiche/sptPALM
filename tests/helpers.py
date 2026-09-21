"""Helpers shared by the tests."""

import numpy as np

from sptpalm.tracks import TrackSet


def brownian_tracks(n_tracks, d_um2_s, dt_s=0.02, n_points=20, seed=0, sigma_loc=0.0,
                    blink_prob=0.0):
    """Synthetic 2D Brownian tracks (positions in um), optionally with missing detections."""
    rng = np.random.default_rng(seed)
    d_um2_s = np.broadcast_to(np.asarray(d_um2_s, dtype=float), (n_tracks,))
    tracks = []
    for i in range(n_tracks):
        step_sigma = np.sqrt(2 * d_um2_s[i] * dt_s)
        xy = np.cumsum(rng.normal(0, step_sigma, size=(n_points, 2)), axis=0) + rng.uniform(5, 20, 2)
        xy += rng.normal(0, sigma_loc, size=xy.shape) if sigma_loc else 0
        frame = np.arange(n_points) + 10
        keep = np.ones(n_points, dtype=bool)
        if blink_prob:
            keep = rng.random(n_points) > blink_prob
            keep[0] = keep[-1] = True
        tracks.append((frame[keep], xy[keep, 0], xy[keep, 1]))
    return TrackSet.from_list(tracks)


def naive_msd(frame, x, y, min_pairs, strict):
    """Direct transcription of MSD_calculation.m for one track (0-based).

    Returns (msd list, weight list)."""
    msd, weight = [], []
    span = frame[-1] - frame[0]
    for lag in range(1, span - (min_pairs - 1) + 1):
        d = []
        for n in range(len(frame)):
            if frame[n] + lag > frame[-1]:
                break
            j = np.flatnonzero(frame == frame[n] + lag)
            if j.size:
                d.append((x[n] - x[j[0]]) ** 2 + (y[n] - y[j[0]]) ** 2)
        d = np.array(d)
        ok = len(d) > min_pairs if strict else len(d) >= min_pairs
        if not ok:
            break
        msd.append(d.mean())
        weight.append(1.0 / d.var(ddof=1))
    return msd, weight
