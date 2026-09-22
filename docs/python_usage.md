# sptpalm (Python) quick start

Python version of the MATLAB *sptPALM_viewer*. Authors: Claude (main author), JB Fiche.
Status: headless analysis core (no GUI yet). The MATLAB code is untouched in `sptPALM_viewer/`.

## Install (development)

```
pip install -e .            # numpy, scipy, pandas, matplotlib
pip install pytest          # to run the tests
```

## Command line

```
sptpalm analyze examples/Test_data/TrackMate_v7_batcher_results \
    --pattern "AC*-spots.csv" --dt-ms 20 --pixel-size-um 0.102
```

Input: the TrackMate *spots* tables (one csv per movie, with the image calibration in
pixels when TrackMate was run). Output, in `<folder>/sptpalm_results` (or `--out`):
`Saved_Diffusion_Coeff.txt`, `Saved_MSD.txt`, `Parameters_analysis.txt` and the figures
`Cumulative_Distribution_LengthStep.png`, `Trajectories_duration.png`,
`Diffusion_distribution_*.png`, `MSD_Curves.png` (same names as the MATLAB version).
`--legacy` reproduces the MATLAB numbers exactly, `--help` lists all the options.

## From Python

```python
from sptpalm import AnalysisParams, ROI, run_analysis
from sptpalm.io.trackmate import load_trackmate_folder
from sptpalm.io.export import write_results

tracks, info = load_trackmate_folder("data/", pixel_size_um=0.102, pattern="*-spots.csv")
params = AnalysisParams(acquisition_time_ms=20, pixel_size_um=0.102)

results = run_analysis(tracks, params)                 # {"all": AnalysisResult}
res = results["all"]
res.n_dapp, res.fit.centers, res.fit.fast_fraction     # numbers of tracks, log10(D) of the populations, ...

# several ROIs (polygons in image pixels: x = column, y = row)
cells = [ROI("cell 1", [(10, 10), (80, 10), (80, 90), (10, 90)]),
         ROI("cell 2", [(120, 30), (200, 30), (200, 110)])]
pooled   = run_analysis(tracks, params, rois=cells)            # {"pooled": ...}
separate = run_analysis(tracks, AnalysisParams(20, 0.102, roi_mode="separate"), rois=cells)
write_results(separate, "out/")                                # one sub-folder per ROI + Summary_ROIs.csv
```

## Parameters (MATLAB name -> Python)

| MATLAB | Python | Default |
|---|---|---|
| AcquisitionTime | `acquisition_time_ms` | required |
| PixelSize | `pixel_size_um` | required |
| MaxBlinks | `max_blinks` | 3 |
| MaxStepLength | `max_step_length_um` | `None` (off) |
| MinTrajLength (7) | `min_points` (= MinTrajLength + 1) | 8 |
| MinNPoint (0.75, had no effect) | `min_fraction` | 0 (off) |
| (none) | `min_duration_frames` | 0 (off) |
| MinNPointMSD | `min_points_msd` | 3 |
| NumberPointsMSDFit | `msd_fit_points` | 4 |
| MaxDisplayTime | `max_display_time_s` | 0.5 |
| Average / fit method | `diffusion_method` = "average" / "fit" | "average" |
| one or two Gaussians, click | `n_gaussians`, `split_threshold` (automatic if `None`) | 2 |
| (all ROIs together) | `roi_mode` = "pooled" / "separate" | "pooled" |

## Differences with the MATLAB version (`legacy=True` / `--legacy` restores the MATLAB behaviour)

* Track length: the MATLAB test compared point indices instead of frames, so `MinTrajLength` = 7 meant at least 8 detections and the 75 % rule never removed anything. Now `min_points`, `min_duration_frames` and `min_fraction` are separate and explicit; the default reproduces the MATLAB results.
* A gap longer than `max_blinks` right before the last detection of a track now cuts the track.
* MSD points need at least `min_points_msd` pairs (MATLAB: strictly more).
* Gaussian width: the fit is unchanged but widths are reported as standard deviations, and with one population the MSD is computed on the tracks within 3 sigma of the peak (MATLAB: 3|s|, about 2.1 sigma).
* Overlapping ROIs: as in MATLAB, a track whose mean position is inside several ROIs is excluded (the number is reported).
* Speed: the MSD, D and fits are vectorised (the example dataset, 23 000 tracks in 10 movies, is loaded in about 2 s and analysed in well under 1 s; the MATLAB README quotes about a minute for the average method and up to 4 minutes for the fit method).

## Tests

```
pytest
```

`tests/test_reference.py` compares the results with the MATLAB outputs stored in `examples/`
(track numbers, all apparent diffusion coefficients and the MSD curves of both populations).
See `docs/matlab_python_comparison.md` for a figure-by-figure version of that comparison.
