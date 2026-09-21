# sptPALM Python conversion: read-through, decisions and proposed structure

Authors: Claude (main author), JB Fiche (second author)
Branch: `python_conversion`
Status: phase 1 (headless core) implemented in `src/sptpalm/`, validated against the MATLAB results; GUI not started
Last update: 2026-09-21

## 1. What was read

Read in full: `README.md`, `Trajectory_analysis_v4.m` (the pipeline driver), `ReconstructTraj_TrackMate_csv_v7.m`, `Load_TrackMate_Tracking_Files_v1.m`, `importTrackMateTracks_csv_v1.m`, `Filter_Trajectories.m`, `MSD_calculation.m`, `Diff_calculation.m`, `FitGaussianDistribution_v3.m`, `Select_Trajectories_ROI.m`, `Save_For_Tesseler_sptPALM.m`, `ImageJ batch macro/TrackMate_batch.py`.

Not yet read (to do before porting them): the GUI files (`sptPALM_viewer.m`, `sptPALM_components.m`, `sptPALM_initialize.m`), the simulation code (`SimuDiff_v4/v5`, `Simulated_Trajectory_analysis_v3`, `emCCD_*`), plotting/visualization files, the `_batch` and `_parallel_computing` variants, and `importTrackMateTracks_xml_v0.m`.

## 2. How the MATLAB pipeline works

1. **Load.** Read TrackMate "all spots" CSV files from a folder. Keep spot ID, track ID, X, Y, frame. Spots without a track are dropped. Positions are converted from pixels to um. Frame numbers of successive movies are offset so that tracks from different files never share a frame index.
2. **Reconstruct.** Each track becomes a 4 x N array (frame, x, y, step length). Tracks containing any zero-length step are rejected outright (TrackMate artefact near image borders). The pooled step-length ECDF and track-duration histogram are plotted as QC.
3. **Filter** (`Filter_Trajectories`). Split or trim tracks at steps longer than `MaxStepLength`, split at gaps longer than `MaxBlinks`, then keep sub-tracks with duration >= `MinTrajLength` and a fraction of populated frames >= `MinNPoint`.
4. **ROI.** Optional polygon ROIs drawn on an average image. A track is kept if its mean position is inside exactly one ROI. Track density is computed per um^2.
5. **MSD** per track, lags 1..L, with at least `MinNPointMSD` displacements per lag. Weight = 1/var of the squared displacements (Saxton 1997). Tracks with fewer than `p` MSD points, or infinite weights, are dropped.
6. **Apparent D** per track from the first `p` MSD points, either "average" (slope between point 1 and p, divided by 4) or "fit" (weighted linear fit, intercept bounded above by min(MSD), D = slope/4). Non-positive D is rejected.
7. **Distribution fit.** Histogram of log10(D), 0.05 bins, fitted by one or two Gaussians. With two, the user clicks an approximate threshold, and the split point is the crossing of the two fitted curves. Outputs are the mean log10(D) per population, the mobile fraction, and the median/mean/SEM MSD curves per population.
8. **Outputs.** PNGs, `Saved_Diffusion_Coeff.txt`, `Saved_MSD.txt`, `Parameters_analysis.txt`, a `.mat` results file, and the SR-Tesseler exports (`Tesseler_diffusion.txt`, `Tesseler_Instant_Velocity.txt`).

## 3. Decisions taken (2026-09-20, with JB)

Each behaviour below has a corrected default and a **legacy switch** (`legacy=True` in the API, "Legacy MATLAB behaviour" option in the GUI) that reproduces the MATLAB result, so old analyses can be re-run and compared.

| Topic | MATLAB behaviour | Python default | Legacy switch |
|---|---|---|---|
| Gaussian width | `A*exp(-((x-x0)/(2s))^2)`, so `s` is not the standard sigma (sigma = sqrt(2)*s) | Standard form `A*exp(-(x-x0)^2/(2 sigma^2))`, sigma reported | Old form. |
| One-population MSD selection | tracks within `x0 +/- 3|s|` (about 2.1 sigma) | tracks within `x0 +/- 3 sigma` | Old window. |
| Filter, gap before last point | the last detection is glued to the previous sub-track even if the gap exceeds `MaxBlinks` | split there like any other gap | Old behaviour. |
| Filter, duration criterion | **Found on the example data.** `MinTrajLength` is compared to a difference of point *indices*, not frames, so it means "at least MinTrajLength+1 detections" (8 points for the default 7). The fraction criterion `MinNPoint` (0.75) is always satisfied (points / index span + 1 is always 1), so it has no effect. | Three explicit criteria: `min_points` (default 8, reproduces MATLAB), `min_duration_frames` (default off), `min_fraction` (default off, 0.75 is the value MATLAB meant to apply). Decided 2026-09-21. | Old behaviour (same defaults). |
| MSD threshold | `nMSD > MinNPointMSD` (strict) | hard threshold, `nMSD >= MinNPointMSD` | Strict `>`. |
| ROI overlap | a track whose mean position lies in two or more ROIs is excluded | kept as is: excluded (ambiguous cell assignment). Reported in the log as a count. | not needed |
| Population labels | code: 1 = slow, 2 = fast (comments say the opposite) | 1 = slow (blue), 2 = fast (orange); mobile fraction = share of population 2; named `slow`/`fast` in the API | not needed |
| X/Y ordering | X = TrackMate POSITION_Y, Y = POSITION_X (image row/column) | kept, documented; columns selected by header name | not needed |

Other confirmed decisions:

- **Multiple ROIs (new feature).** ROIs are named. A `roi_mode` option selects `pooled` (current behaviour: all ROIs analysed together) or `separate` (each ROI gets its own MSD, D distribution, fit, density, output folder and one row in a summary table). In separate mode the number of Gaussians and the slow/fast split threshold can be set per ROI, since different cells may differ.
- **Efficiency.** MSD is computed with vectorised numpy over dense, NaN-padded frame arrays. The weighted linear fit for D has a closed form (with the bounded intercept handled analytically), so it is evaluated for all tracks at once instead of calling a generic fit per track. Loading uses `pandas.read_csv`. The `_parallel_computing` duplicates disappear. Expected results are identical to MATLAB up to numerical tolerance, and this will be checked in the tests.
- **GUI required**, with the same organisation as the MATLAB one: a control panel with the six sections (loading data, analysis parameters, diffusion analysis, statistics, visualization tool, simulation) and a separate display panel. Built as a thin layer over the headless core.

## 3b. Validation on `examples/Test_data/TrackMate_v7_batcher_results` (2026-09-21)

A short Python re-implementation of the MATLAB chain (loader, zero-step rejection, filter, MSD, "average" Dapp), run on the ten `AC*-spots.csv` files with the MATLAB default parameters (dt = 20 ms, pixel = 102 nm, `MaxBlinks` = 3, `MinTrajLength` = 7, `p` = 4, `MinNPointMSD` = 3, no `MaxStepLength`), reproduces the top-level MATLAB results of that folder exactly:

| Stage | MATLAB (`Parameters_analysis.txt`) | Python |
|---|---|---|
| Initial trajectories | 22925 | 22925 |
| After filter | 5783 | 5783 |
| Validated for MSD | 5681 | 5681 |
| Validated for Dapp | 5227 | 5227 |

The 5227 apparent diffusion coefficients equal the ones in `Saved_Diffusion_Coeff.txt`, in the same order (files read AC0 to AC9), to 1e-14. The sub-folder `sptPALM_viewer results` (22983 initial trajectories) does not correspond to the current CSV files and is not used as reference.

Other results on this dataset:

- The filter frame-span interpretation (`MinTrajLength` = 7 frames plus the 75% rule) would give 6005 tracks, not 5783; the MATLAB numbers are only reproduced with the index-based reading described in section 3.
- The gap-before-last-detection bug changes nothing here (legacy and fixed give the same count).
- With the hard MSD threshold (`>=` instead of `>`), 5765 tracks pass the MSD step and 5294 the Dapp step (+67 tracks, +1.3%).

Phase 1 result: the Python package (`legacy=True`) reproduces the MATLAB outputs of this folder: track counts at every stage, all 5227 apparent diffusion coefficients, and the MSD curves (mean, median, standard error) of both populations of `Saved_MSD.txt`, all to floating-point precision (relative difference below 1e-9). The two-Gaussian fit was started from an automatically estimated split threshold and converged to the same solution as the click-started MATLAB fit. With the corrected defaults the numbers are 5783 / 5765 / 5294 tracks, log10(D) = -1.825 and -1.043, mobile fraction 74.7 %.

## 4. Proposed package layout

```
sptPALM/
  pyproject.toml
  DEVLOG.md
  src/sptpalm/
    io/
      trackmate.py      # CSV (spots table) reader, by header name; XML reader later if wanted
      results.py        # save/load analysis results (HDF5 + JSON parameters); optional legacy .mat import
      export.py         # SR-Tesseler text files, Saved_MSD / Saved_Diffusion_Coeff
    tracks.py           # Tracks container: pandas DataFrame (track_id, frame, x_um, y_um, movie)
    filtering.py        # misconnection removal, blink splitting, duration / density criteria
    roi.py              # named polygon ROIs, pooled/separate modes, density, average image (tifffile)
    msd.py              # per-track MSD + Saxton weights (vectorised)
    diffusion.py        # Dapp: average and weighted-fit methods (closed form, vectorised)
    fitting.py          # 1- and 2-Gaussian fits of log10(D) (scipy.optimize.curve_fit)
    pipeline.py         # run_analysis(params) -> AnalysisResult; replaces Trajectory_analysis_v4
    params.py           # dataclass with all analysis parameters, MATLAB-equivalent defaults, legacy flag
    plotting.py         # step-length ECDF, duration histogram, D distribution, MSD, trajectories
    simulation/         # port of SimuDiff later, once the analysis core is validated
    cli.py              # `sptpalm analyze <folder> --dt 20 --pixel-size 0.16 ...`
    gui/                # PySide6 control + display panels, matplotlib canvases
  tests/
  examples/             # a notebook reproducing the README workflow
```

Design choices:

- **Data model.** One long-format DataFrame instead of cell arrays of 3 x N matrices.
- **No interactive dialogs in the core.** ROI polygons and the two-Gaussian split threshold are function arguments; the GUI supplies them interactively. This removes the `h.batch` special cases.
- **Parameters** live in one dataclass, saved with the results, replacing the values read from GUI text fields.
- **Dependencies:** numpy, scipy, pandas, matplotlib, tifffile, h5py; PySide6 for the GUI.
- **Removed:** all MTT code and the `_parallel_computing` duplicates.

## 5. Phases

1. **Headless core.** TrackMate loader, filtering, MSD, Dapp, Gaussian fits, multi-ROI (pooled/separate), text/plot outputs, CLI, tests. Validate against MATLAB outputs on a small dataset, in legacy mode first. **Done (2026-09-21)**, except saving/loading a full results file (needed for "Load previous analysis" in the GUI) and the Tesseler / localisation exports, which move to phase 2.
2. **GUI.** Control panel and display panel mirroring the MATLAB layout, ROI drawing, previous-analysis loading and plotting, track overlay visualization tool with movie export.
3. **Simulation module** and remaining tools (Tesseler export polish, XML reader if still needed).

## 6. Open items

1. **Track selection defaults.** Implemented option (c): separate `min_points`, `min_duration_frames` and `min_fraction`, with defaults that reproduce the MATLAB results (frame-based rules off). JB to confirm that this is the intended reading of "I like your idea better" and whether the GUI should propose the frame-based rules (e.g. 0.75) by default.
2. The top-level result set of the example folder is used as reference (see 3b).
2. Simulation module: port, and is the XML TrackMate reader still needed? (not urgent)
3. Confirm the MSD-threshold interpretation: "hard threshold" implemented as `nMSD >= MinNPointMSD`.
