# sptPALM

Analysis of single-particle tracking (sptPALM) data: track filtering, MSD and
apparent-diffusion-coefficient calculation, and population fitting from
**TrackMate** tracks.

This repository is being converted from the original MATLAB software
(`sptPALM_viewer`, still available in the `legacy_matlab` branch and in the
`sptPALM_viewer/` folder here) to a Python package, `sptpalm`, so the tool can
run without a MATLAB license and stay easy to maintain. The Python version
drops support for the older **MTT** tracker and works from TrackMate tracks
only. See `docs/python_conversion_plan.md` for the full comparison with the
MATLAB code and `DEVLOG.md` for a dated history of the conversion.

Status: the headless analysis core (loading, filtering, MSD, diffusion
coefficients, population fits, text/plot outputs, command line) is done and
validated against the MATLAB outputs. The GUI has not been ported yet — see
`docs/python_conversion_plan.md` for the plan.

## Input data

The package reads TrackMate "all spots" CSV tables (one file per movie),
exported with **File > Export tracks to CSV** or the TrackMate batcher, with
image dimensions calibrated in pixels (not nm/µm). Positions are converted to
µm using the pixel size you provide.

## Installation

Requires Python ≥ 3.9. From the root of the repository, on the
`python_conversion` branch:

```bash
# with uv (recommended)
uv venv .venv --python 3.12
source .venv/bin/activate
uv pip install -e ".[test]"

# or with plain pip, in a virtual environment of your choice
python -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
```

`.[test]` also installs `pytest` so you can run the test suite. Two other
optional extras are available if/when you need them: `.[image]` (`tifffile`,
for reading movies) and `.[gui]` (`PySide6`, for the future GUI).

## Running an analysis

### Command line

```bash
sptpalm analyze examples/Test_data/TrackMate_v7_batcher_results \
    --pattern "AC*-spots.csv" --dt-ms 20 --pixel-size-um 0.102
```

This reads every CSV file matching `--pattern` in the given folder,
reconstructs and filters the tracks, computes the MSD and apparent diffusion
coefficient of each track, fits the log10(D) distribution, and writes the
results to `<folder>/sptpalm_results` (change with `--out`):
`Saved_Diffusion_Coeff.txt`, `Saved_MSD.txt`, `Parameters_analysis.txt`, and
the figures `Cumulative_Distribution_LengthStep.png`,
`Trajectories_duration.png`, `Diffusion_distribution_*.png`,
`MSD_Curves.png` (same names as the MATLAB version).

Add `--legacy` to reproduce the MATLAB numbers exactly (see "Differences with
the MATLAB version" below). Run `sptpalm analyze --help` for the full list of
options (`--max-blinks`, `--min-points`, `--method`, `--gaussians`, ROI and
plotting options, etc.).

### From Python

```python
from sptpalm import AnalysisParams, ROI, run_analysis
from sptpalm.io.trackmate import load_trackmate_folder
from sptpalm.io.export import write_results

tracks, info = load_trackmate_folder(
    "examples/Test_data/TrackMate_v7_batcher_results",
    pixel_size_um=0.102,
    pattern="AC*-spots.csv",
)
params = AnalysisParams(acquisition_time_ms=20, pixel_size_um=0.102)

results = run_analysis(tracks, params)   # {"all": AnalysisResult}
res = results["all"]
res.n_dapp, res.fit.centers, res.fit.fast_fraction

write_results(results, "out/")
```

Multiple ROIs (polygons in image pixels, `x` = column, `y` = row) can be
analysed pooled together or separately:

```python
cells = [ROI("cell 1", [(10, 10), (80, 10), (80, 90), (10, 90)]),
         ROI("cell 2", [(120, 30), (200, 30), (200, 110)])]
pooled   = run_analysis(tracks, params, rois=cells)                       # {"pooled": ...}
separate = run_analysis(tracks, AnalysisParams(20, 0.102, roi_mode="separate"), rois=cells)
write_results(separate, "out/")   # one sub-folder per ROI + Summary_ROIs.csv
```

See `docs/python_usage.md` for the full MATLAB-to-Python parameter table and
a longer description of the differences with the MATLAB version.

## Tests

```bash
pytest
```

`tests/test_core.py` covers the individual pieces (filtering, MSD, fitting,
...) with synthetic data. `tests/test_reference.py` compares the package's
output against the MATLAB results stored in
`examples/Test_data/TrackMate_v7_batcher_results` (track counts at every
stage, all apparent diffusion coefficients, and the MSD curves of both
populations) — it is skipped automatically if that example data is not
present. Set `SPTPALM_EXAMPLES` to point elsewhere if you keep the reference
data in a different location.

## Differences with the MATLAB version

Every corrected behaviour below has a legacy switch
(`legacy=True` in the Python API, `--legacy` on the command line) that
reproduces the MATLAB numbers exactly, so old analyses stay reproducible.

* **Track length.** The MATLAB filter compared point *indices* instead of
  frame numbers, so `MinTrajLength = 7` actually required 8 detections and
  the 75% "populated fraction" rule had no effect. The Python version has
  three explicit, independent criteria (`min_points`, `min_duration_frames`,
  `min_fraction`); the defaults reproduce the MATLAB results.
* A gap longer than `max_blinks` right before the **last** detection of a
  track now cuts the track there (MATLAB glued it to the previous segment).
* The MSD-point threshold is now a hard `>=` (MATLAB used a strict `>`).
* Gaussian fit widths are reported as standard deviations (the MATLAB
  parametrisation was not the standard form); with a single population the
  MSD is computed on tracks within 3σ of the peak (MATLAB: about 2.1σ).
* Overlapping ROIs: as in MATLAB, a track whose mean position falls in
  several ROIs is excluded (now reported as a count).
* Multiple ROIs can be analysed **pooled** or **separately** (new feature,
  not in MATLAB).
* Speed: MSD, diffusion coefficients and fits are vectorised. On the example
  dataset (23,000 tracks, 10 movies) loading takes about 2 s and the analysis
  well under 1 s, versus about a minute (average method) to several minutes
  (fit method) for the MATLAB version.

## Legacy MATLAB version

The original MATLAB software (`sptPALM_viewer`, MATLAB R2019a, tested on
Windows 10 and Ubuntu 18.04) is preserved as-is in the `legacy_matlab` branch
and in the `sptPALM_viewer/` folder. It supports both MTT and TrackMate
tracks and includes a GUI (control panel + display panel), a visualization
tool to overlay tracks on the acquired images, and a Brownian-motion
simulation module — none of which are ported to Python yet (see
`docs/python_conversion_plan.md`, phases 2 and 3). To use it, clone the
repository, check out `legacy_matlab`, and run `sptPALM_viewer.m` in MATLAB.
