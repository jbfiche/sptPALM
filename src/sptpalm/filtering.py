"""Track filtering.

Port of ``Filter_Trajectories.m``. Three steps:

1. (optional, ``max_step_length_um``) removal of mis-connections: steps longer
   than the threshold are either trimmed at the extremities, removed when they
   are due to a single wrong detection, or used to split the track in two.
2. Splitting at blinks: a gap between two detections longer than
   ``max_blinks`` frames cuts the track.
3. Selection of the (sub-)tracks that are long enough: ``min_points``,
   ``min_duration_frames`` and ``min_fraction`` (see :mod:`sptpalm.params`).

Differences with MATLAB (``legacy=True`` restores the MATLAB behaviour):

* a gap larger than ``max_blinks`` right before the *last* detection cuts the
  track (MATLAB glued the last detection to the previous segment);
* the selection criteria are explicit (see :mod:`sptpalm.params`).
"""

from __future__ import annotations

from typing import List, Tuple

import numpy as np

from .params import AnalysisParams
from .tracks import TrackSet

_Track = Tuple[np.ndarray, np.ndarray, np.ndarray]


# ---------------------------------------------------------------------------
# 1. Mis-connections
# ---------------------------------------------------------------------------
def _steps(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.hypot(np.diff(x), np.diff(y))


def _fix_track_misconnections(track: _Track, max_step: float) -> Tuple[_Track, List[_Track]]:
    """Fix one track. Returns the (possibly shortened) track and the new tracks
    created by splits, in the order MATLAB appended them."""
    frame, x, y = (a.copy() for a in track)
    d = _steps(x, y)
    ok = d[d <= max_step]
    mean_d = ok.mean() if ok.size else np.nan
    std_d = ok.std(ddof=1) if ok.size > 1 else 0.0
    extra: List[_Track] = []

    bad = np.flatnonzero(d > max_step)
    while bad.size > 0:
        m = int(bad[0])
        if m == 0:                                   # first step: drop the first point
            frame, x, y = frame[1:], x[1:], y[1:]
        elif m == len(d) - 1:                        # last step: drop the last point
            frame, x, y = frame[:-1], x[:-1], y[:-1]
        else:
            d2 = np.hypot(x[m] - x[m + 2], y[m] - y[m + 2])
            if d2 < mean_d + 3 * std_d:              # a single wrong detection
                keep = np.ones(len(x), dtype=bool)
                keep[m + 1] = False
                frame, x, y = frame[keep], x[keep], y[keep]
            else:                                    # two different tracks: split
                extra.append((frame[m + 1:], x[m + 1:], y[m + 1:]))
                frame, x, y = frame[:m + 1], x[:m + 1], y[:m + 1]
        d = _steps(x, y)
        bad = np.flatnonzero(d > max_step)
    return (frame, x, y), extra


def remove_misconnections(ts: TrackSet, max_step_um: float) -> TrackSet:
    """Apply the mis-connection removal to all the tracks.

    Tracks created by splitting are appended at the end and are checked in
    turn, as in MATLAB.
    """
    if ts.n_tracks == 0:
        return ts
    steps = np.hypot(np.diff(ts.x), np.diff(ts.y))
    same = np.ones(len(steps), dtype=bool)
    same[ts.offsets[1:-1] - 1] = False
    over = np.zeros(ts.n_tracks, dtype=bool)
    over[ts.track_index[1:][same & (steps > max_step_um)]] = True
    if not over.any():
        return ts

    work: List[_Track] = [ts.get(i) for i in range(ts.n_tracks)]
    movie = list(ts.movie)
    tid = list(ts.track_id)
    need = list(over)
    i = 0
    while i < len(work):
        if need[i]:
            work[i], extra = _fix_track_misconnections(work[i], max_step_um)
            for e in extra:
                work.append(e)
                movie.append(movie[i])
                tid.append(tid[i])
                need.append(True)
        i += 1
    return TrackSet.from_list(work, movie, tid)


# ---------------------------------------------------------------------------
# 2 and 3. Blinks and selection
# ---------------------------------------------------------------------------
def split_and_select(ts: TrackSet, params: AnalysisParams) -> TrackSet:
    """Split the tracks at long blinks and keep the long enough sub-tracks."""
    if ts.n_tracks == 0:
        return ts
    n = len(ts.frame)
    track_start = np.zeros(n, dtype=bool)
    track_start[ts.offsets[:-1]] = True
    track_last = np.zeros(n, dtype=bool)
    track_last[ts.offsets[1:] - 1] = True

    # cut[j] is True when a new sub-track starts at point j (j >= 1)
    gap = np.zeros(n, dtype=bool)
    gap[1:] = (np.diff(ts.frame) > params.max_blinks + 1) & ~track_start[1:]
    if params.legacy:
        # MATLAB never cut before the last detection of a track
        gap &= ~track_last
    seg_start = track_start | gap
    starts = np.flatnonzero(seg_start)
    ends = np.concatenate([starts[1:], [n]])           # exclusive
    n_points = ends - starts
    span = ts.frame[ends - 1] - ts.frame[starts]

    keep = n_points >= max(params.min_points, 2)
    if params.min_duration_frames > 0:
        keep &= span >= params.min_duration_frames
    if params.min_fraction > 0:
        keep &= n_points / (span + 1.0) >= params.min_fraction

    seg_track = ts.track_index[starts]
    sel = np.flatnonzero(keep)
    if sel.size == 0:
        return TrackSet.empty()
    lengths = n_points[sel]
    new_offsets = np.concatenate([[0], np.cumsum(lengths)])
    idx = np.repeat(starts[sel] - new_offsets[:-1], lengths) + np.arange(new_offsets[-1])
    return TrackSet(ts.frame[idx], ts.x[idx], ts.y[idx], new_offsets,
                    ts.movie[seg_track[sel]], ts.track_id[seg_track[sel]])


def filter_tracks(ts: TrackSet, params: AnalysisParams) -> TrackSet:
    """Full filtering chain (mis-connections, blinks, selection)."""
    if params.max_step_length_um is not None:
        ts = remove_misconnections(ts, params.max_step_length_um)
    return split_and_select(ts, params)
