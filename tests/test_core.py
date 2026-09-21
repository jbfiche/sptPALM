"""Unit tests of the analysis core (no external data needed)."""

import numpy as np
from scipy.optimize import lsq_linear

from helpers import brownian_tracks, naive_msd
from sptpalm import AnalysisParams, ROI, TrackSet, run_analysis
from sptpalm.diffusion import apparent_diffusion, dapp_fit
from sptpalm.filtering import filter_tracks, remove_misconnections
from sptpalm.fitting import fit_distribution, auto_split_threshold
from sptpalm.io.trackmate import read_spots_csv, tracks_from_spots
from sptpalm.msd import compute_msd
from sptpalm.roi import assign_tracks


# ---------------------------------------------------------------- tracks
def test_trackset_select_and_positions():
    ts = TrackSet.from_list([(np.array([0, 1, 2]), np.array([0., 1, 2]), np.array([0., 0, 0])),
                             (np.array([5, 6]), np.array([10., 12]), np.array([4., 4]))])
    assert ts.n_tracks == 2 and list(ts.lengths) == [3, 2]
    assert np.allclose(ts.mean_positions(), [[1, 0], [11, 4]])
    assert np.allclose(ts.step_lengths(), [1, 1, 2])   # no step between the two tracks
    sub = ts.select(np.array([False, True]))
    assert sub.n_tracks == 1 and list(sub.frame) == [5, 6]


# ---------------------------------------------------------------- io
def test_read_spots_and_zero_step_rejection(tmp_path):
    csv = tmp_path / "x-spots.csv"
    csv.write_text(
        "LABEL,ID,TRACK_ID,QUALITY,POSITION_X,POSITION_Y,POSITION_Z,POSITION_T,FRAME\n"
        "Label,Spot ID,Track ID,Quality,X,Y,Z,T,Frame\n"
        "Label,Spot ID,Track ID,Quality,X,Y,Z,T,Frame\n"
        ",,,(quality),(pixels),(pixels),(pixels),(frames),\n"
        "ID1,1,0,1,10,10,0,0,0\n"
        "ID2,2,0,1,11,10,0,1,1\n"
        "ID3,3,1,1,20,20,0,0,0\n"
        "ID4,4,1,1,20,20,0,1,1\n"          # zero step -> track 1 rejected
        "ID5,5,,1,30,30,0,0,0\n"           # not in a track
        "ID6,6,2,1,5,5,0,0,3\n"
        "ID7,7,2,1,6,5,0,1,1\n")           # unsorted frames
    spots = read_spots_csv(str(csv))
    assert len(spots) == 6
    ts, rejected = tracks_from_spots(spots, 0.1)
    assert rejected == 1 and ts.n_tracks == 2
    assert list(ts.track_id) == [0, 2]
    assert list(ts.get(1)[0]) == [1, 3]                 # sorted by frame
    assert np.allclose(ts.get(0)[1], [1.0, 1.1])        # um


# ---------------------------------------------------------------- filtering
def _one_track(frames, x=None, y=None):
    frames = np.asarray(frames)
    x = np.arange(len(frames), dtype=float) * 0.05 if x is None else np.asarray(x, dtype=float)
    y = np.sin(np.arange(len(frames))) * 0.03 if y is None else np.asarray(y, dtype=float)
    return TrackSet.from_list([(frames, x, y)])


def test_gap_before_last_detection():
    ts = _one_track(list(range(1, 11)) + [30])
    p_new = AnalysisParams(20, 0.1, max_blinks=3, min_points=8)
    p_old = AnalysisParams(20, 0.1, max_blinks=3, min_points=8, legacy=True)
    new = filter_tracks(ts, p_new)
    old = filter_tracks(ts, p_old)
    assert new.n_tracks == 1 and new.lengths[0] == 10 and new.last_frame[0] == 10
    assert old.n_tracks == 1 and old.lengths[0] == 11 and old.last_frame[0] == 30   # MATLAB quirk


def test_split_at_blinks_and_criteria():
    frames = list(range(0, 12)) + list(range(20, 32))      # gap of 8 frames in the middle
    ts = _one_track(frames)
    out = filter_tracks(ts, AnalysisParams(20, 0.1, max_blinks=3, min_points=8))
    assert out.n_tracks == 2 and list(out.lengths) == [12, 12]
    out = filter_tracks(ts, AnalysisParams(20, 0.1, max_blinks=10, min_points=8))
    assert out.n_tracks == 1 and out.lengths[0] == 24
    # duration and fraction criteria
    sparse = _one_track([0, 3, 6, 9, 12, 15, 18, 21])       # 8 points over 22 frames
    base = AnalysisParams(20, 0.1, max_blinks=3, min_points=8)
    assert filter_tracks(sparse, base).n_tracks == 1
    assert filter_tracks(sparse, AnalysisParams(20, 0.1, max_blinks=3, min_points=8,
                                                min_fraction=0.75)).n_tracks == 0
    assert filter_tracks(sparse, AnalysisParams(20, 0.1, max_blinks=3, min_points=8,
                                                min_duration_frames=30)).n_tracks == 0


def test_misconnections():
    n = 12
    rng = np.random.default_rng(11)
    x = rng.normal(0, 0.03, n)             # confined track, mean step ~0.04 um
    y = rng.normal(0, 0.03, n)
    # a single wrong detection far away: the point is removed
    x1, y1 = x.copy(), y.copy()
    y1[6] += 2.0
    out = remove_misconnections(_one_track(range(n), x1, y1), 0.5)
    assert out.n_tracks == 1 and out.lengths[0] == n - 1
    # two different tracks connected by a jump: the track is split
    x2 = np.concatenate([x[:6], x[:6] + 10.0])
    y2 = np.concatenate([y[:6], y[:6]])
    out = remove_misconnections(_one_track(range(n), x2, y2), 0.5)
    assert out.n_tracks == 2 and list(out.lengths) == [6, 6]
    # jump at the first detection: it is dropped
    x3 = x.copy()
    x3[0] += 9.0
    out = remove_misconnections(_one_track(range(n), x3, y), 0.5)
    assert out.n_tracks == 1 and out.lengths[0] == n - 1
    # nothing to fix
    assert remove_misconnections(_one_track(range(n), x, y), 0.5).lengths[0] == n


# ---------------------------------------------------------------- MSD
def test_msd_matches_transcription_of_matlab():
    ts = brownian_tracks(60, 0.05, n_points=40, seed=3, blink_prob=0.2)
    for legacy in (False, True):
        res = compute_msd(ts, width=10, p=4, min_pairs=3, legacy=legacy)
        for i in range(ts.n_tracks):
            f, x, y = ts.get(i)
            m, w = naive_msd(f, x, y, 3, strict=legacy)
            k = min(len(m), 10)
            assert res.n_lags[i] == (len(m) if legacy else min(len(m), 10))
            assert np.allclose(res.msd[i, :k], m[:k], rtol=1e-10)
            assert np.allclose(res.weight[i, :k], w[:k], rtol=1e-8)
            assert np.all(np.isnan(res.msd[i, k:]))


def test_msd_brownian_motion():
    ts = brownian_tracks(400, 0.1, n_points=30, seed=1)
    res = compute_msd(ts, width=6, p=4)
    mean_msd = np.nanmean(res.msd, axis=0)
    expected = 4 * 0.1 * 0.02 * np.arange(1, 7)
    assert np.allclose(mean_msd, expected, rtol=0.08)
    d, ok = apparent_diffusion(res.msd[res.accepted], res.weight[res.accepted], 4, 0.02)
    assert abs(np.median(d) / 0.1 - 1) < 0.25 and ok.mean() > 0.9


def test_hard_threshold_accepts_more():
    ts = brownian_tracks(80, 0.05, n_points=7, seed=4)   # lag 4 has exactly 3 pairs
    new = compute_msd(ts, width=5, p=4, min_pairs=3, legacy=False)
    old = compute_msd(ts, width=5, p=4, min_pairs=3, legacy=True)
    assert new.accepted.sum() > old.accepted.sum()


# ---------------------------------------------------------------- D
def test_fit_closed_form_matches_bounded_least_squares():
    rng = np.random.default_rng(7)
    p, dt = 5, 0.02
    t = dt * np.arange(1, p + 1)
    n_over = 0
    for _ in range(200):
        a_true = rng.uniform(0.002, 0.3)
        b_true = rng.uniform(-0.003, 0.004)            # can be negative: constraint may bind
        y = a_true * t + b_true + rng.normal(0, 0.002, p)
        y = np.abs(y) + 1e-5
        w = rng.uniform(0.5, 2, p)
        got = dapp_fit(y[None, :], w[None, :], p, dt)[0]
        A = np.column_stack([t, np.ones(p)]) * np.sqrt(w)[:, None]
        sol = lsq_linear(A, y * np.sqrt(w), bounds=([-np.inf, -np.inf], [np.inf, y.min()]),
                         tol=1e-14)
        n_over += sol.x[1] >= y.min() - 1e-12
        assert np.isclose(got, sol.x[0] / 4, rtol=1e-6, atol=1e-9)
    assert n_over > 0          # the bounded case was exercised


# ---------------------------------------------------------------- distribution fit
def _mixture(seed=0, n_slow=1500, n_fast=3500):
    rng = np.random.default_rng(seed)
    logd = np.concatenate([rng.normal(-1.9, 0.35, n_slow), rng.normal(-1.0, 0.2, n_fast)])
    msd = np.abs(rng.normal(0.01, 0.003, (len(logd), 25))) + 1e-4
    return logd, msd


def test_two_gaussian_fit_recovers_populations():
    logd, msd = _mixture()
    fit = fit_distribution(logd, msd, 2, 25)                   # automatic split threshold
    assert np.allclose(fit.centers, [-1.9, -1.0], atol=0.06)
    assert np.allclose(fit.sigmas, [0.35, 0.2], atol=0.06)
    assert abs(fit.fast_fraction - 0.7) < 0.04
    assert fit.centers[0] < fit.centers[1] and set(np.unique(fit.labels)) <= {0, 1, 2}
    assert fit.msd_summary[0].shape == (25, 3)
    fit_b = fit_distribution(logd, msd, 2, 25, split_threshold=-1.5)
    assert np.allclose(fit.centers, fit_b.centers, atol=1e-3)
    assert -1.9 < auto_split_threshold(logd) < -1.0


def test_single_gaussian_window_legacy_vs_new():
    rng = np.random.default_rng(2)
    logd = rng.normal(-1.2, 0.3, 4000)
    msd = np.abs(rng.normal(0.01, 0.003, (4000, 25))) + 1e-4
    new = fit_distribution(logd, msd, 1, 25)
    old = fit_distribution(logd, msd, 1, 25, legacy=True)
    assert abs(new.sigmas[0] - 0.3) < 0.03
    assert np.allclose(new.centers, old.centers) and np.allclose(new.curve_total, old.curve_total)
    # 3 sigma keeps ~99.7 % of the tracks, MATLAB's 3|s| only ~96.6 %
    assert new.msd_selected.mean() > 0.99 and 0.95 < old.msd_selected.mean() < 0.98


# ---------------------------------------------------------------- ROI
def _two_cells():
    """Tracks in two square cells of a 60x60 px image (pixel 0.1 um) plus a few outside/overlap."""
    px = 0.1
    rng = np.random.default_rng(5)
    tracks = []
    centers = [(10, 10)] * 150 + [(45, 45)] * 100 + [(30, 30)] * 5 + [(55, 5)] * 5
    for cx, cy in centers:
        step = np.sqrt(2 * 0.08 * 0.02)
        xy = np.cumsum(rng.normal(0, step, (15, 2)), axis=0) + np.array([cx, cy]) * px
        tracks.append((np.arange(15), xy[:, 0], xy[:, 1]))
    return TrackSet.from_list(tracks), px


def test_roi_assignment_overlap_and_modes():
    ts, px = _two_cells()
    r1 = ROI("cell 1", [(0, 0), (20, 0), (20, 20), (0, 20)])
    r2 = ROI("cell 2", [(35, 35), (55, 35), (55, 55), (35, 55)])
    r3 = ROI("big", [(0, 0), (40, 0), (40, 40), (0, 40)])     # overlaps cell 1 (and holds (30,30))
    assert r1.area_px2() == 400.0
    idx, n_amb = assign_tracks(ts, [r1, r2], px)
    assert (idx == 0).sum() == 150 and (idx == 1).sum() == 100 and n_amb == 0
    assert (idx == -1).sum() == 10
    idx, n_amb = assign_tracks(ts, [r1, r3], px)
    assert n_amb == 150                                       # cell-1 tracks lie in both ROIs

    base = dict(min_points=8, n_gaussians=1)
    pooled = run_analysis(ts, AnalysisParams(20, px, **base), rois=[r1, r2])
    assert list(pooled) == ["pooled"] and pooled["pooled"].n_roi == 250
    assert np.isclose(pooled["pooled"].density, 250 / ((400 + 400) * px ** 2))
    sep = run_analysis(ts, AnalysisParams(20, px, roi_mode="separate", **base), rois=[r1, r2])
    assert list(sep) == ["cell 1", "cell 2"]
    assert sep["cell 1"].n_roi == 150 and sep["cell 2"].n_roi == 100
    assert np.isclose(sep["cell 1"].density, 150 / (400 * px ** 2))
    whole = run_analysis(ts, AnalysisParams(20, px, **base), image_shape=(60, 60))["all"]
    assert whole.n_roi == 260 and np.isclose(whole.density, 260 / (3600 * px ** 2))
