"""Regions of interest (ROI).

A ROI is a named polygon drawn on the image (typically one cell). Vertices
are in image pixels: ``x`` is the column and ``y`` the row, i.e. the TrackMate
POSITION_X and POSITION_Y expressed in pixels.

A track belongs to a ROI when its *mean position* is inside the polygon. A
track whose mean position is inside two or more ROIs cannot be attributed to a
single cell (proteins do not cross cell walls, so such tracks are most likely
artefacts) and is excluded, as in the MATLAB version.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np
from matplotlib.path import Path

from .tracks import TrackSet


@dataclass
class ROI:
    name: str
    vertices_px: np.ndarray  # (n, 2): x (column), y (row)

    def __post_init__(self) -> None:
        v = np.asarray(self.vertices_px, dtype=float)
        if v.ndim != 2 or v.shape[1] != 2 or len(v) < 3:
            raise ValueError("a ROI needs at least 3 vertices, as an (n, 2) array")
        self.vertices_px = v

    @property
    def closed_vertices(self) -> np.ndarray:
        v = self.vertices_px
        return v if np.allclose(v[0], v[-1]) else np.vstack([v, v[0]])

    def area_px2(self, rounded: bool = False) -> float:
        """Polygon area in pixel^2 (shoelace formula).

        MATLAB rounded the area to an integer number of pixels (``rounded=True``).
        """
        v = self.closed_vertices
        area = 0.5 * abs(np.sum(v[:-1, 0] * v[1:, 1]) - np.sum(v[:-1, 1] * v[1:, 0]))
        return float(round(area)) if rounded else float(area)

    def contains(self, points_px: np.ndarray) -> np.ndarray:
        """Boolean mask of the (n, 2) points (x, y in pixels) inside the polygon."""
        if len(points_px) == 0:
            return np.zeros(0, dtype=bool)
        return Path(self.closed_vertices).contains_points(np.asarray(points_px))


def assign_tracks(ts: TrackSet, rois: Sequence[ROI], pixel_size_um: float
                  ) -> Tuple[np.ndarray, int]:
    """Attribute each track to a ROI.

    Returns ``(roi_index, n_ambiguous)``: for every track the index of its ROI,
    or -1 when it is outside all the ROIs or inside several of them (the number
    of the latter is ``n_ambiguous``).
    """
    pos_px = ts.mean_positions() / pixel_size_um
    inside = np.column_stack([r.contains(pos_px) for r in rois]) if len(rois) else \
        np.zeros((ts.n_tracks, 0), dtype=bool)
    count = inside.sum(axis=1)
    index = np.where(count == 1, inside.argmax(axis=1) if inside.shape[1] else -1, -1)
    return index.astype(int), int((count > 1).sum())
