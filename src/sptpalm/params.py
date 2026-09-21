"""Analysis parameters.

All the values that the MATLAB GUI read from text fields live here, in one
dataclass that is saved together with the results.

Correspondence with the MATLAB sptPALM_viewer parameters
--------------------------------------------------------
=========================  =========================================
MATLAB                     Python
=========================  =========================================
AcquisitionTime            ``acquisition_time_ms``
PixelSize                  ``pixel_size_um``
MaxBlinks                  ``max_blinks``
MaxStepLength              ``max_step_length_um`` (``None`` = off)
MinTrajLength              ``min_points`` - 1 (see below)
MinNPoint (0.75)           ``min_fraction`` (see below)
MinNPointMSD               ``min_points_msd``
NumberPointsMSDFit         ``msd_fit_points``
MaxDisplayTime             ``max_display_time_s``
DiffusionCalculationMethod ``diffusion_method`` ("average" or "fit")
=========================  =========================================

Track selection criteria
------------------------
In the MATLAB code the duration test compared *point indices* instead of
frame numbers. As a consequence ``MinTrajLength`` = 7 really meant "at least
8 detections", and the ``MinNPoint`` fraction test was always satisfied.
The Python version exposes three explicit criteria:

* ``min_points``: minimum number of detections in a (sub-)track. The default
  (8) reproduces the MATLAB results obtained with ``MinTrajLength`` = 7.
* ``min_duration_frames``: minimum duration, in frames, between the first and
  the last detection (0 = not used).
* ``min_fraction``: minimum fraction of the frames between the first and the
  last detection that contain a detection (0 = not used, 0.75 is the value
  the MATLAB code was meant to apply).
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Optional


@dataclass
class AnalysisParams:
    # Acquisition -----------------------------------------------------------
    acquisition_time_ms: float
    pixel_size_um: float

    # Track filtering -------------------------------------------------------
    max_blinks: int = 3
    max_step_length_um: Optional[float] = None
    min_points: int = 8
    min_duration_frames: int = 0
    min_fraction: float = 0.0

    # MSD and apparent diffusion coefficient --------------------------------
    min_points_msd: int = 3
    msd_fit_points: int = 4
    max_display_time_s: float = 0.5
    diffusion_method: str = "average"  # "average" or "fit"

    # Distribution of log10(D) ----------------------------------------------
    n_gaussians: int = 2
    split_threshold: Optional[float] = None  # log10(D) used to start the 2-Gaussian fit

    # ROI handling ----------------------------------------------------------
    roi_mode: str = "pooled"  # "pooled" or "separate"

    # Reproduce the MATLAB behaviour (see module docstring of each module) ---
    legacy: bool = False

    def __post_init__(self) -> None:
        if self.diffusion_method not in ("average", "fit"):
            raise ValueError("diffusion_method must be 'average' or 'fit'")
        if self.n_gaussians not in (1, 2):
            raise ValueError("n_gaussians must be 1 or 2")
        if self.roi_mode not in ("pooled", "separate"):
            raise ValueError("roi_mode must be 'pooled' or 'separate'")
        if self.msd_fit_points < 2:
            raise ValueError("msd_fit_points must be at least 2")
        if self.acquisition_time_ms <= 0 or self.pixel_size_um <= 0:
            raise ValueError("acquisition_time_ms and pixel_size_um must be positive")

    # ------------------------------------------------------------------
    @property
    def dt_s(self) -> float:
        return self.acquisition_time_ms / 1000.0

    @property
    def display_lags(self) -> int:
        """Number of MSD points shown/summarised (MaxDisplayTime / dt)."""
        return int(round(self.max_display_time_s * 1000.0 / self.acquisition_time_ms))

    @classmethod
    def from_matlab(cls, *, acquisition_time_ms: float, pixel_size_um: float,
                    min_traj_length: int = 7, **kwargs) -> "AnalysisParams":
        """Build parameters from the MATLAB ``MinTrajLength`` convention."""
        return cls(acquisition_time_ms=acquisition_time_ms,
                   pixel_size_um=pixel_size_um,
                   min_points=int(min_traj_length) + 1, **kwargs)

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "AnalysisParams":
        names = {f.name for f in dataclasses.fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in names})
