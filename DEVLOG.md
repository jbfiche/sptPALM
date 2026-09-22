# Development log: Python conversion

Authors: Claude (main author), JB Fiche (second author)
Branch: `python_conversion`

## 2026-09-20: read-through and proposal

- Read the core of the MATLAB analysis chain (TrackMate loading, filtering, MSD, apparent D, Gaussian fits, ROI selection, Tesseler export, TrackMate batch macro). See `docs/python_conversion_plan.md` for the summary.
- Identified MATLAB behaviours that need a decision before porting (Gaussian width convention, X/Y ordering, filter edge cases, MSD threshold, ROI overlap, population labelling).
- Proposed a `src/sptpalm/` package layout and three phases.
- No code ported yet. No MATLAB file modified.
- Not yet read: GUI files, simulation code, plotting/visualization files, batch and parallel variants, XML importer.

## 2026-09-20: decisions with JB

- Gaussian width: the MATLAB form is an error (`s` is not the standard sigma). Python uses the standard form; a legacy switch keeps the old one. The fitted peak position (and therefore mean log10(D), two-population split and mobile fraction) is not affected. Only the MSD-selection window of the one-population case changes (old `x0 +/- 3|s|` is about 2.1 sigma, new is 3 sigma).
- X/Y ordering (X = TrackMate POSITION_Y) is intentional and is kept.
- Filter: a gap longer than `MaxBlinks` right before the last detection will split the track (MATLAB glues that last point on). Legacy switch available.
- MSD threshold: hard threshold, `nMSD >= MinNPointMSD`. Legacy switch keeps the strict `>`.
- ROI overlap: tracks lying in several ROIs stay excluded (ambiguous cell assignment).
- Population labels: 1 = slow (blue), 2 = fast (orange), mobile fraction = share of 2. The MATLAB code already behaves this way; only its comments were reversed.
- New feature: several ROIs can be analysed either pooled or separately.
- A GUI with the same organisation as the MATLAB one is required (PySide6 planned).
- All corrected behaviours get a legacy switch. Efficiency is a goal (vectorised MSD and closed-form weighted fit).
- Test dataset from JB expected on 2026-09-21.

## 2026-09-21: validation on the example dataset

- JB added `examples/Test_data` (10 raw movies, TrackMate v7 batcher outputs, single-stack TrackMate output, and previous sptPALM_viewer results).
- Scratch Python re-implementation (not committed) of loader + zero-step rejection + filter + MSD + "average" Dapp reproduces the MATLAB numbers of `TrackMate_v7_batcher_results` exactly: 22925 / 5783 / 5681 / 5227 tracks, and all 5227 Dapp values to 1e-14. The `sptPALM_viewer results` sub-folder (22983 initial tracks) is from an older run and does not match the current CSVs.
- Found: in `Filter_Trajectories.m` the duration test uses index differences instead of frame numbers. `MinTrajLength` therefore means at least `MinTrajLength`+1 detections, and the `MinNPoint` fraction criterion is always 1, so it does nothing. Only the index-based reading reproduces 5783 (a frame-span reading gives 6005). Decision on the Python semantics pending (plan, section 6).
- The last-detection gap bug has no effect on this dataset. The hard MSD threshold (`>=`) adds 67 tracks (5294 vs 5227 with Dapp).
- Format confirmed from `AC*-spots.csv`: 4 header lines, columns LABEL, ID, TRACK_ID, QUALITY, POSITION_X, POSITION_Y, POSITION_Z, POSITION_T, FRAME, ... so the loader should read columns by name and skip rows 2 to 4. MATLAB X = POSITION_Y, Y = POSITION_X (intentional).
- Not yet compared: the two-Gaussian fit, MSD curves (`Saved_MSD.txt`), the `.mat` contents, ROI-related outputs.

## 2026-09-21: phase 1, headless core (`src/sptpalm/`)

- Decision (JB): track selection with separate criteria (minimum detections, minimum duration, minimum fraction of populated frames); frame-based rules off by default so that MATLAB results are reproduced. Legacy switch restores every MATLAB behaviour.
- Added `.gitignore` (Python and MATLAB autosave files, heavy example data: raw tif, xml, mat).
- New package `sptpalm`, dataclass-based, no GUI dialogs:
  - `params.py` (all parameters, MATLAB correspondence), `tracks.py` (flat-array track container), `io/trackmate.py` (spots csv reader by column name, zero-step rejection, frame offsets), `filtering.py` (mis-connections, blink split, selection; vectorised), `msd.py` (vectorised MSD with Saxton weights), `diffusion.py` (average and weighted-fit D, closed form), `fitting.py` (1 and 2 Gaussians, populations, MSD summaries), `roi.py` (named polygons, overlap exclusion), `pipeline.py` (`run_analysis`, pooled or separate ROIs), `io/export.py` (Saved_*.txt, Parameters_analysis.txt, summary table per ROI), `plotting.py`, `cli.py` (`sptpalm analyze ...`).
- Fixes with legacy switch: gap before the last detection cuts the track; hard MSD threshold (`>=`); Gaussian width reported as standard sigma and 3-sigma window for the one-population MSD; explicit track-length criteria.
- Fit method: closed-form weighted linear fit with the intercept bounded by min(MSD), verified against a bounded least-squares solver on random data (tests). Not yet compared with MATLAB `fit` output on real data (no reference available: all example results use the average method).
- Tests (`tests/`): 12 unit tests and 3 reference tests against the MATLAB results in `examples/`. All pass. (Run in a scratch environment without pytest using a small runner; `pytest` itself was not available there, so the first run under pytest is up to JB.)
- Validation on `TrackMate_v7_batcher_results` in legacy mode: 22925 / 5783 / 5681 / 5227 tracks, all D values, and the mean/median/SEM MSD curves of both populations match MATLAB (relative difference < 1e-9). Corrected defaults: 5783 / 5765 / 5294 tracks, log10(D) = -1.825 and -1.043, mobile fraction 74.7 %. Loading takes about 2 s and the analysis about 0.1 s.
- Not done yet: saving/loading a full results file, Tesseler and localisation exports, visualization tool, ROI drawing on an image, simulation module, GUI, XML reader, single-stack `export-spots.csv` check.
- Git: no shell is available in this session, so nothing has been committed. Files are written in the working tree of branch `python_conversion`.

## 2026-09-21: branch layout

- JB reworked the branches by hand. `master`, `legacy_matlab` and `python_conversion` all pointed at the same commit (a56b2ea, the MATLAB code, equal to `origin/master`); HEAD is on `python_conversion`. `legacy_matlab` freezes the MATLAB version, `python_conversion` holds this work. Nothing was pushed yet at that point.
- The repository was re-cloned, which dropped the hidden `.gitignore` from the working tree (so `git add .gitignore` failed). It was written again.

## 2026-09-22: README

- `README.md` still described only the MATLAB GUI. Rewritten to lead with the Python package: overview, Installation (uv or pip, `.[test]`/`.[image]`/`.[gui]` extras), Running an analysis (CLI and Python API, using the example dataset), Tests, and the differences-with-MATLAB summary (moved here from `docs/python_usage.md`, which keeps the full parameter table). The original MATLAB walkthrough is kept, trimmed, under "Legacy MATLAB version", pointing to the `legacy_matlab` branch and `sptPALM_viewer/`.

## 2026-09-22: MATLAB vs Python validation report

- Added `docs/matlab_python_comparison.md` and `docs/validation/*.png`: side-by-side comparison of the Python package against the MATLAB results already saved in `examples/Test_data/TrackMate_v7_batcher_results` (`sptpalm analyze --legacy` vs the corrected defaults vs MATLAB's own `Saved_Diffusion_Coeff.txt` / `Saved_MSD.txt` / PNGs).
- Confirms, with figures, what phase 1 validation established numerically: in legacy mode, track counts match at every stage, the 5227 apparent D values match MATLAB's to ~1e-14, the MSD curves of both populations match to ~1e-13, and MATLAB's own two-Gaussian fit figure (log10(D) = -1.81 / -1.04, mobile fraction 75%) matches the Python fit (-1.814 / -1.042, 74.9%) — this is the first time the MATLAB fit numbers themselves (read off its PNG, not just the reproduced D/MSD values) were compared directly, since the `.mat` results file is MATLAB v7.3 (HDF5) and `h5py` is not installable in this sandbox (no network access to PyPI/apt beyond the pre-approved index).
- Also shows the (expected, small) effect of the corrected defaults: 5294 vs 5227 D-valid tracks, same distribution shape.
