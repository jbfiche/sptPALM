"""Reader for the TrackMate "spots" tables (csv).

TrackMate (v7 and later) can export a table with all the detected spots and
the track they belong to (``*-spots.csv`` with the batcher, or the "Spots"
table of the GUI). The file starts with a header line followed by three extra
lines (long names, short names and units) that are skipped. The columns are
selected by name, so the extra feature columns do not matter.

Make sure that in Fiji the image X/Y calibration was in pixels when TrackMate
was run: positions are converted to um here with ``pixel_size_um``.
"""

from __future__ import annotations

import csv
import glob
import os
import warnings
from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np
import pandas as pd

from ..tracks import TrackSet

REQUIRED_COLUMNS = ("TRACK_ID", "POSITION_X", "POSITION_Y", "FRAME")


@dataclass
class LoadInfo:
    files: List[str] = field(default_factory=list)
    n_tracks_in_files: List[int] = field(default_factory=list)  # tracks found in each file
    n_rejected_zero_step: List[int] = field(default_factory=list)
    skipped_files: List[str] = field(default_factory=list)

    @property
    def n_tracks_total(self) -> int:
        return int(sum(self.n_tracks_in_files))

    @property
    def n_rejected_total(self) -> int:
        return int(sum(self.n_rejected_zero_step))


def _read_header(path: str) -> Tuple[List[str], int]:
    """Return the column names and the number of extra header lines."""
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        try:
            id_col = header.index("ID")
        except ValueError:
            id_col = header.index("TRACK_ID")
        n_extra = 0
        for row in reader:
            try:
                float(row[id_col])
                break
            except (ValueError, IndexError):
                n_extra += 1
    return header, n_extra


def read_spots_csv(path: str) -> pd.DataFrame:
    """Read a TrackMate spots table.

    Returns a DataFrame with columns ``track_id`` (int), ``x_px``, ``y_px``
    and ``frame`` (int). Spots that do not belong to a track are dropped.
    """
    header, n_extra = _read_header(path)
    missing = [c for c in REQUIRED_COLUMNS if c not in header]
    if missing:
        raise ValueError(f"{os.path.basename(path)}: missing TrackMate columns {missing}")
    df = pd.read_csv(path, skiprows=range(1, 1 + n_extra), usecols=list(REQUIRED_COLUMNS))
    df = df.dropna(subset=["TRACK_ID"])
    out = pd.DataFrame({
        "track_id": df["TRACK_ID"].astype(np.int64).to_numpy(),
        "x_px": df["POSITION_X"].to_numpy(dtype=np.float64),
        "y_px": df["POSITION_Y"].to_numpy(dtype=np.float64),
        "frame": df["FRAME"].astype(np.int64).to_numpy(),
    })
    return out


def tracks_from_spots(spots: pd.DataFrame, pixel_size_um: float, movie: int = 0
                      ) -> Tuple[TrackSet, int]:
    """Group the spots in tracks and reject the ones with a zero-length step.

    TrackMate sometimes returns identical positions for detections close to
    the image border (or the same localisation several times in a track).
    Such tracks give a step length of exactly 0 and are removed, as in the
    MATLAB version.

    Returns the TrackSet (frames not yet offset) and the number of rejected tracks.
    """
    if len(spots) == 0:
        return TrackSet.empty(), 0
    order = np.lexsort((spots["frame"].to_numpy(), spots["track_id"].to_numpy()))
    tid = spots["track_id"].to_numpy()[order]
    frame = spots["frame"].to_numpy()[order]
    x = spots["x_px"].to_numpy()[order] * pixel_size_um
    y = spots["y_px"].to_numpy()[order] * pixel_size_um

    new_track = np.concatenate([[True], tid[1:] != tid[:-1]])
    starts = np.flatnonzero(new_track)
    offsets = np.concatenate([starts, [len(tid)]])
    n_tracks = len(starts)
    track_index = np.cumsum(new_track) - 1

    d = np.hypot(np.diff(x), np.diff(y))
    same = ~new_track[1:]                       # step j belongs to a track if point j+1 is not a start
    zero_step_tracks = np.unique(track_index[1:][same & (d == 0)])
    good = np.ones(n_tracks, dtype=bool)
    good[zero_step_tracks] = False

    ts = TrackSet(frame, x, y, offsets, np.full(n_tracks, movie, dtype=np.int32), tid[starts])
    return ts.select(good), int(n_tracks - good.sum())


def load_trackmate_folder(folder: str, pixel_size_um: float, pattern: str = "*.csv",
                          legacy_frame_offset: bool = True) -> Tuple[TrackSet, LoadInfo]:
    """Load every TrackMate spots file of ``folder``.

    Files are read in alphabetical order; each file is one movie. To keep the
    frame numbers of different movies apart (needed when the localisations are
    exported), the frames of a movie are shifted by the sum of the durations of
    the previous ones, as in the MATLAB version (``legacy_frame_offset``).
    """
    files = sorted(glob.glob(os.path.join(folder, pattern)))
    info = LoadInfo()
    sets = []
    frame_offset = 0
    movie = 0
    for path in files:
        try:
            spots = read_spots_csv(path)
        except Exception as exc:  # not a TrackMate table (or unreadable)
            warnings.warn(f"skipping {os.path.basename(path)}: {exc}")
            info.skipped_files.append(path)
            continue
        if len(spots) < 2:
            info.skipped_files.append(path)
            continue
        ts, n_rejected = tracks_from_spots(spots, pixel_size_um, movie)
        max_frame = max(1, int(spots["frame"].max()))
        if legacy_frame_offset:
            ts.frame = ts.frame + frame_offset
            frame_offset += max_frame
        sets.append(ts)
        info.files.append(path)
        info.n_tracks_in_files.append(ts.n_tracks + n_rejected)
        info.n_rejected_zero_step.append(n_rejected)
        movie += 1
    return TrackSet.concatenate(sets), info
