"""Gaussian fit of the distribution of log10(D) and population analysis.

Port of ``FitGaussianDistribution_v3.m``.

The histogram of log10(D) (bins of 0.05, expressed in % of the tracks) is
fitted by one or two Gaussians. With two Gaussians, population 1 is the
*slow* one (lowest mean log10(D)) and population 2 the *fast* one; the
tracks are assigned to a population according to the crossing point of the two
fitted curves.

Width convention
----------------
The MATLAB code fitted ``A exp(-((x-x0)/(2 s))^2)``, where ``s`` is *not* the
standard deviation (sigma = sqrt(2) s). The fit itself is unchanged (so the
fitted curves, the peak positions, the split between populations and the mobile
fraction are identical) but the results are reported with the standard sigma
and, for one population, the tracks used for the MSD are those within
3 sigma of the peak. With ``legacy=True`` the MATLAB window of 3 |s| (about 2.1
sigma) is used instead.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import curve_fit

SQRT2 = np.sqrt(2.0)
BIN_WIDTH = 0.05


# ---------------------------------------------------------------------------
# Model functions (MATLAB parametrisation, s = sigma / sqrt(2))
# ---------------------------------------------------------------------------
def gaussian_s(x, s, x0, a):
    return a * np.exp(-((x - x0) / (2.0 * s)) ** 2)


def gaussian_sigma(x, sigma, x0, a):
    """Standard Gaussian, amplitude ``a`` at the peak."""
    return a * np.exp(-0.5 * ((x - x0) / sigma) ** 2)


def _two_gaussians(x, s1, x01, a1, s2, x02, a2):
    return gaussian_s(x, s1, x01, a1) + gaussian_s(x, s2, x02, a2)


# ---------------------------------------------------------------------------
@dataclass
class DiffusionFit:
    n_gaussians: int
    centers: np.ndarray            # mean log10(D) of each population (slow first)
    sigmas: np.ndarray             # standard deviations (log10 units)
    amplitudes: np.ndarray         # amplitudes (% of tracks per bin)
    r_squared: float
    hist_bin: np.ndarray           # bin centres
    hist_values: np.ndarray        # % of tracks per bin (all tracks)
    hist_values_pop: np.ndarray    # (n_bins, n_gaussians) % of tracks per bin, by population
    curve_x: np.ndarray            # fine grid for the fitted curves
    curve_pop: np.ndarray          # (len(curve_x), n_gaussians) fitted curve of each population
    curve_total: np.ndarray        # sum of the fitted curves
    labels: np.ndarray             # per track: 1 = slow, 2 = fast (0 = on the boundary); 1 for one population
    split: Optional[float] = None  # log10(D) where the two fitted curves cross
    fast_fraction: Optional[float] = None
    msd_selected: np.ndarray = None  # tracks used for the MSD curve (one population)
    msd_summary: list = field(default_factory=list)  # per population (n_lags, 3): mean, SEM, median

    @property
    def d_centers(self) -> np.ndarray:
        """D (um^2/s) at the peak of each population."""
        return 10.0 ** self.centers

    @property
    def mobile_fraction_percent(self) -> Optional[float]:
        return None if self.fast_fraction is None else round(100 * self.fast_fraction)


# ---------------------------------------------------------------------------
def log_histogram(log_d: np.ndarray):
    """Histogram of log10(D) with the MATLAB binning.

    Edges are multiples of 0.05 from floor(20 min)/20 to floor(20 max)/20:
    the (few) values above the last edge are not counted, as in MATLAB.
    Returns (edges, centers, percent per bin, total count used).
    """
    kmin = int(np.floor(20 * np.min(log_d)))
    kmax = int(np.floor(20 * np.max(log_d)))
    if kmax <= kmin:
        kmax = kmin + 1
    edges = (kmin + np.arange(kmax - kmin + 1)) / 20.0
    counts, _ = np.histogram(log_d, bins=edges)
    total = counts.sum()
    centers = 0.5 * (edges[:-1] + edges[1:])
    return edges, centers, counts * 100.0 / total, int(total), counts


def auto_split_threshold(log_d: np.ndarray, n_iter: int = 200) -> float:
    """Estimate a threshold between two populations (2-component Gaussian mixture, EM)."""
    x = np.asarray(log_d, dtype=float)
    mu = np.percentile(x, [25, 75])
    sd = np.full(2, max(x.std() / 2, 1e-3))
    pi = np.array([0.5, 0.5])
    for _ in range(n_iter):
        dens = np.stack([pi[k] * np.exp(-0.5 * ((x - mu[k]) / sd[k]) ** 2) / sd[k] for k in range(2)])
        resp = dens / np.maximum(dens.sum(axis=0), 1e-300)
        nk = resp.sum(axis=1)
        mu_new = (resp * x).sum(axis=1) / nk
        sd = np.sqrt((resp * (x - mu_new[:, None]) ** 2).sum(axis=1) / nk).clip(1e-3)
        pi = nk / nk.sum()
        if np.allclose(mu_new, mu, atol=1e-8):
            mu = mu_new
            break
        mu = mu_new
    lo, hi = np.sort(mu)
    grid = np.linspace(lo, hi, 501)
    k_lo, k_hi = np.argsort(mu)
    f = np.abs(pi[k_lo] * np.exp(-0.5 * ((grid - mu[k_lo]) / sd[k_lo]) ** 2) / sd[k_lo]
               - pi[k_hi] * np.exp(-0.5 * ((grid - mu[k_hi]) / sd[k_hi]) ** 2) / sd[k_hi])
    return float(grid[np.argmin(f)])


def _r_squared(y, yfit) -> float:
    ss_res = np.sum((y - yfit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return float(1.0 - ss_res / ss_tot)


def msd_summary(msd: np.ndarray, n_lags: int) -> np.ndarray:
    """Mean, standard error of the mean and median of the MSD at each lag.

    Only positive, finite values are used (as in MATLAB). Result shape (n_lags, 3).
    """
    out = np.full((n_lags, 3), np.nan)
    for k in range(min(n_lags, msd.shape[1])):
        v = msd[:, k]
        v = v[np.isfinite(v) & (v > 0)]
        if v.size:
            sem = v.std(ddof=1) / np.sqrt(v.size) if v.size > 1 else 0.0
            out[k] = (v.mean(), sem, np.median(v))
    return out


def fit_distribution(log_d: np.ndarray, msd: np.ndarray, n_gaussians: int, n_lags: int,
                     split_threshold: Optional[float] = None, legacy: bool = False) -> DiffusionFit:
    """Fit the distribution of log10(D) and summarise the MSD of each population.

    Parameters
    ----------
    log_d : log10 of the apparent diffusion coefficients (one per track).
    msd : MSD of the same tracks (rows), used for the population MSD curves.
    n_gaussians : 1 or 2.
    n_lags : number of MSD lags to summarise (MaxDisplayTime / acquisition time).
    split_threshold : approximate log10(D) separating the two populations, used
        to start the two-Gaussian fit (the click on the graph in MATLAB). If
        None it is estimated automatically.
    """
    log_d = np.asarray(log_d, dtype=float)
    edges, centers, values, total, _ = log_histogram(log_d)

    if n_gaussians == 1:
        i0 = int(np.argmax(values))
        popt, _ = curve_fit(gaussian_s, centers, values, p0=[1.0, centers[i0], values[i0]],
                            maxfev=20000)
        s, x0, a = popt
        sigma = SQRT2 * abs(s)
        fit_at_bins = gaussian_s(centers, *popt)
        curve_x = centers
        curve = fit_at_bins[:, None]
        half_width = 3 * abs(s) if legacy else 3 * sigma
        selected = (log_d > x0 - half_width) & (log_d < x0 + half_width)
        return DiffusionFit(
            n_gaussians=1, centers=np.array([x0]), sigmas=np.array([sigma]),
            amplitudes=np.array([a]), r_squared=_r_squared(values, fit_at_bins),
            hist_bin=centers, hist_values=values, hist_values_pop=values[:, None],
            curve_x=curve_x, curve_pop=curve, curve_total=fit_at_bins,
            labels=np.ones(len(log_d), dtype=int), msd_selected=selected,
            msd_summary=[msd_summary(msd[selected], n_lags)])

    # ---- two populations -------------------------------------------------
    if split_threshold is None:
        split_threshold = auto_split_threshold(log_d)
    d_hi = log_d[log_d > split_threshold]
    d_lo = log_d[log_d <= split_threshold]
    if d_hi.size < 3 or d_lo.size < 3:
        raise ValueError("the split threshold leaves fewer than 3 tracks in one population")

    def start_amplitude(x_mean):
        above = np.flatnonzero(centers > x_mean)
        return values[above[0]] if above.size else values[-1]

    x_hi, x_lo = d_hi.mean(), d_lo.mean()
    p0 = [d_hi.std(ddof=1), x_hi, start_amplitude(x_hi),
          d_lo.std(ddof=1), x_lo, start_amplitude(x_lo)]
    try:
        popt, _ = curve_fit(_two_gaussians, centers, values, p0=p0, maxfev=50000)
    except RuntimeError as exc:
        raise ValueError(f"the two-Gaussian fit did not converge: {exc}") from exc
    s1, x01, a1, s2, x02, a2 = popt

    # population 1 = the one with the lowest mean
    if x01 < x02:
        (s_slow, c_slow, a_slow), (s_fast, c_fast, a_fast) = (s1, x01, a1), (s2, x02, a2)
    else:
        (s_slow, c_slow, a_slow), (s_fast, c_fast, a_fast) = (s2, x02, a2), (s1, x01, a1)

    n_fine = int(np.floor((centers[-1] - centers[0]) / 0.01 + 1e-9))
    grid = centers[0] + 0.01 * np.arange(n_fine + 1)
    g_slow = gaussian_s(grid, s_slow, c_slow, a_slow)
    g_fast = gaussian_s(grid, s_fast, c_fast, a_fast)
    between = np.flatnonzero((grid >= c_slow) & (grid <= c_fast))
    if between.size:
        split = float(grid[between[np.argmin(np.abs(g_slow[between] - g_fast[between]))]])
    else:
        split = 0.5 * (c_slow + c_fast)

    labels = np.zeros(len(log_d), dtype=int)
    labels[log_d < split] = 1
    labels[log_d > split] = 2
    n_slow, n_fast = int((labels == 1).sum()), int((labels == 2).sum())
    fast_fraction = n_fast / (n_slow + n_fast)

    def pop_hist(mask):
        c, _ = np.histogram(log_d[mask], bins=edges)
        return c * 100.0 / total

    hist_pop = np.column_stack([pop_hist(labels == 1), pop_hist(labels == 2)])
    fit_at_bins = _two_gaussians(centers, *popt)
    return DiffusionFit(
        n_gaussians=2, centers=np.array([c_slow, c_fast]),
        sigmas=SQRT2 * np.abs([s_slow, s_fast]), amplitudes=np.array([a_slow, a_fast]),
        r_squared=_r_squared(values, fit_at_bins),
        hist_bin=centers, hist_values=values, hist_values_pop=hist_pop,
        curve_x=grid, curve_pop=np.column_stack([g_slow, g_fast]), curve_total=g_slow + g_fast,
        labels=labels, split=split, fast_fraction=fast_fraction,
        msd_summary=[msd_summary(msd[labels == 1], n_lags), msd_summary(msd[labels == 2], n_lags)])
