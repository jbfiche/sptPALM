"""sptPALM: analysis of single-particle tracking data (TrackMate tracks).

Authors: Claude (main author), JB Fiche.

Python conversion of the MATLAB software *sptPALM_viewer*.
"""

from .filtering import filter_tracks
from .params import AnalysisParams
from .pipeline import AnalysisResult, run_analysis
from .roi import ROI
from .tracks import TrackSet

__version__ = "0.1.0.dev0"

__all__ = ["AnalysisParams", "AnalysisResult", "ROI", "TrackSet", "filter_tracks", "run_analysis"]
