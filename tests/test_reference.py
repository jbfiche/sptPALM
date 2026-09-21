"""Comparison with the results of the MATLAB sptPALM_viewer on the example dataset.

The example data are in ``examples/Test_data/TrackMate_v7_batcher_results``
(TrackMate spot tables of 10 movies and the MATLAB results obtained with the
default parameters). Set the environment variable ``SPTPALM_EXAMPLES`` to use
another location. The tests are skipped when the data are not there.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd

from sptpalm import AnalysisParams, run_analysis
from sptpalm.io.export import write_results
from sptpalm.io.trackmate import load_trackmate_folder

try:
    import pytest
except ImportError:  # the tests can also be run without pytest
    pytest = None


def _folder():
    root = Path(os.environ.get("SPTPALM_EXAMPLES", Path(__file__).resolve().parents[1] / "examples"))
    return root / "Test_data" / "TrackMate_v7_batcher_results"


FOLDER = _folder()
HAVE_DATA = FOLDER.is_dir() and (FOLDER / "Saved_Diffusion_Coeff.txt").exists()
skip_without_data = pytest.mark.skipif(not HAVE_DATA, reason="example dataset not found") \
    if pytest else (lambda f: f)

_cache = {}


def _tracks():
    if "tracks" not in _cache:
        _cache["tracks"], _cache["info"] = load_trackmate_folder(str(FOLDER), 0.102, "AC*-spots.csv")
    return _cache["tracks"]


@skip_without_data
def test_track_counts_and_diffusion_coefficients_match_matlab():
    tracks = _tracks()
    params = AnalysisParams(acquisition_time_ms=20, pixel_size_um=0.102, legacy=True)
    res = run_analysis(tracks, params)["all"]
    # numbers of tracks in Parameters_analysis.txt
    assert (res.n_loaded, res.n_filtered, res.n_msd, res.n_dapp) == (22925, 5783, 5681, 5227)
    ref = pd.read_csv(FOLDER / "Saved_Diffusion_Coeff.txt")
    # tiny D values (2e-6) come from a difference of two nearly equal MSD values:
    # floating point noise of about 1e-10 (relative) is expected
    assert np.allclose(res.d_app, ref["D"].to_numpy(), rtol=1e-6, atol=0)
    assert np.allclose(res.log_d, ref["logD"].to_numpy(), rtol=0, atol=1e-8)


@skip_without_data
def test_population_msd_matches_matlab():
    tracks = _tracks()
    params = AnalysisParams(acquisition_time_ms=20, pixel_size_um=0.102, legacy=True)
    res = run_analysis(tracks, params)["all"]
    ref = pd.read_csv(FOLDER / "Saved_MSD.txt")
    assert np.allclose(res.fit.msd_summary[0][:, 0], ref["MSD1_Average"], rtol=1e-9)
    assert np.allclose(res.fit.msd_summary[0][:, 2], ref["MSD1_Median"], rtol=1e-9)
    assert np.allclose(res.fit.msd_summary[0][:, 1], ref["MSD1_Error_Bar"], rtol=1e-8)
    assert np.allclose(res.fit.msd_summary[1][:, 0], ref["MSD2_Average"], rtol=1e-9)
    assert np.allclose(res.fit.msd_summary[1][:, 2], ref["MSD2_Median"], rtol=1e-9)
    assert np.allclose(res.fit.msd_summary[1][:, 1], ref["MSD2_Error_Bar"], rtol=1e-8)


@skip_without_data
def test_hard_threshold_default_gives_more_tracks(tmp_path):
    tracks = _tracks()
    res = run_analysis(tracks, AnalysisParams(acquisition_time_ms=20, pixel_size_um=0.102))["all"]
    assert (res.n_filtered, res.n_msd, res.n_dapp) == (5783, 5765, 5294)
    write_results({"all": res}, str(tmp_path))
    written = pd.read_csv(tmp_path / "Saved_Diffusion_Coeff.txt")
    assert list(written.columns) == ["logD", "D"] and len(written) == 5294
    msd = pd.read_csv(tmp_path / "Saved_MSD.txt")
    assert list(msd.columns)[:4] == ["Time_s", "MSD1_Average", "MSD1_Median", "MSD1_Error_Bar"]
