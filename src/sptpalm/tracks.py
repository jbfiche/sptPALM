"""Container for a set of tracks.

The tracks are stored in flat, concatenated numpy arrays (one entry per
detection) plus an ``offsets`` array, so that all the heavy operations can be
vectorised. Track *i* owns the detections ``offsets[i]:offsets[i+1]``, sorted
by frame.

Coordinates are in um and follow the TrackMate convention: ``x`` is
POSITION_X and ``y`` is POSITION_Y (image column and row).

Note: the MATLAB version stored X = POSITION_Y and Y = POSITION_X. This is
intentional there and only matters when exporting coordinates; the analysis
(MSD, D) is not affected.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Sequence, Tuple

import numpy as np


@dataclass
class TrackSet:
    frame: np.ndarray            # int64, one per detection
    x: np.ndarray                # float64, um
    y: np.ndarray                # float64, um
    offsets: np.ndarray          # int64, length n_tracks + 1
    movie: np.ndarray = field(default=None)     # int32, one per track
    track_id: np.ndarray = field(default=None)  # int64, TrackMate id inside its movie

    def __post_init__(self) -> None:
        self.frame = np.asarray(self.frame, dtype=np.int64)
        self.x = np.asarray(self.x, dtype=np.float64)
        self.y = np.asarray(self.y, dtype=np.float64)
        self.offsets = np.asarray(self.offsets, dtype=np.int64)
        n = len(self.offsets) - 1
        if self.movie is None:
            self.movie = np.zeros(n, dtype=np.int32)
        if self.track_id is None:
            self.track_id = np.arange(n, dtype=np.int64)
        self.movie = np.asarray(self.movie, dtype=np.int32)
        self.track_id = np.asarray(self.track_id, dtype=np.int64)

    # ------------------------------------------------------------------
    @classmethod
    def empty(cls) -> "TrackSet":
        return cls(np.empty(0, np.int64), np.empty(0), np.empty(0), np.zeros(1, np.int64))

    @classmethod
    def from_list(cls, tracks: Sequence[Tuple[np.ndarray, np.ndarray, np.ndarray]],
                  movie: Sequence[int] = None, track_id: Sequence[int] = None) -> "TrackSet":
        """Build from a list of ``(frame, x, y)`` arrays."""
        if len(tracks) == 0:
            return cls.empty()
        lengths = np.array([len(t[0]) for t in tracks], dtype=np.int64)
        offsets = np.concatenate([[0], np.cumsum(lengths)])
        frame = np.concatenate([np.asarray(t[0]) for t in tracks])
        x = np.concatenate([np.asarray(t[1]) for t in tracks])
        y = np.concatenate([np.asarray(t[2]) for t in tracks])
        return cls(frame, x, y, offsets, movie, track_id)

    @classmethod
    def concatenate(cls, sets: Iterable["TrackSet"]) -> "TrackSet":
        sets = [s for s in sets if s.n_tracks > 0]
        if not sets:
            return cls.empty()
        frame = np.concatenate([s.frame for s in sets])
        x = np.concatenate([s.x for s in sets])
        y = np.concatenate([s.y for s in sets])
        lengths = np.concatenate([s.lengths for s in sets])
        offsets = np.concatenate([[0], np.cumsum(lengths)])
        movie = np.concatenate([s.movie for s in sets])
        track_id = np.concatenate([s.track_id for s in sets])
        return cls(frame, x, y, offsets, movie, track_id)

    # ------------------------------------------------------------------
    @property
    def n_tracks(self) -> int:
        return len(self.offsets) - 1

    def __len__(self) -> int:
        return self.n_tracks

    @property
    def lengths(self) -> np.ndarray:
        """Number of detections of each track."""
        return np.diff(self.offsets)

    @property
    def track_index(self) -> np.ndarray:
        """Track index of every detection."""
        return np.repeat(np.arange(self.n_tracks), self.lengths)

    @property
    def first_frame(self) -> np.ndarray:
        return self.frame[self.offsets[:-1]]

    @property
    def last_frame(self) -> np.ndarray:
        return self.frame[self.offsets[1:] - 1]

    @property
    def duration_frames(self) -> np.ndarray:
        """Last frame minus first frame."""
        return self.last_frame - self.first_frame

    def get(self, i: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        a, b = self.offsets[i], self.offsets[i + 1]
        return self.frame[a:b], self.x[a:b], self.y[a:b]

    def select(self, which) -> "TrackSet":
        """Sub-set of tracks given a boolean mask or an index array."""
        which = np.asarray(which)
        if which.dtype == bool:
            which = np.flatnonzero(which)
        which = which.astype(np.int64)
        if which.size == 0:
            return TrackSet.empty()
        lengths = self.lengths[which]
        starts = self.offsets[which]
        new_offsets = np.concatenate([[0], np.cumsum(lengths)])
        # index of every detection to keep
        idx = np.repeat(starts - new_offsets[:-1], lengths) + np.arange(new_offsets[-1])
        return TrackSet(self.frame[idx], self.x[idx], self.y[idx], new_offsets,
                        self.movie[which], self.track_id[which])

    # ------------------------------------------------------------------
    def step_lengths(self) -> np.ndarray:
        """Distance (um) between consecutive detections, all tracks pooled."""
        d = np.hypot(np.diff(self.x), np.diff(self.y))
        # remove the "steps" that connect two different tracks
        keep = np.ones(len(d), dtype=bool)
        boundaries = self.offsets[1:-1] - 1
        keep[boundaries[(boundaries >= 0) & (boundaries < len(d))]] = False
        return d[keep]

    def mean_positions(self) -> np.ndarray:
        """Mean (x, y) of each track, shape (n_tracks, 2), in um."""
        if self.n_tracks == 0:
            return np.empty((0, 2))
        starts = self.offsets[:-1]
        n = self.lengths
        mx = np.add.reduceat(self.x, starts) / n
        my = np.add.reduceat(self.y, starts) / n
        return np.column_stack([mx, my])

    def to_dataframe(self):
        import pandas as pd
        return pd.DataFrame({
            "track": self.track_index,
            "movie": np.repeat(self.movie, self.lengths),
            "track_id": np.repeat(self.track_id, self.lengths),
            "frame": self.frame, "x_um": self.x, "y_um": self.y,
        })
