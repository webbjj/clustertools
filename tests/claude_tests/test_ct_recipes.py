import numpy as np

from clustertools.util.recipes import (
    nbinmaker,
    binmaker,
    roaming_nbinmaker,
    roaming_binmaker,
    power_law_distribution_function,
    interpolate,
    minimum_distance,
    distance,
    mean_prof,
    smooth,
    dx_function,
    tapered_dx_function,
    x_hist,
    bin_index,
)


def test_nbinmaker_counts_and_coverage():
    rng = np.random.default_rng(1)
    x = rng.uniform(0, 100, 997)

    x_lower, x_mid, x_upper, x_hist = nbinmaker(x, nbin=10)

    # every point should land in exactly one bin
    total = 0
    for lo, hi in zip(x_lower, x_upper):
        if hi == x_upper[-1]:
            indx = (x >= lo) & (x <= hi)
        else:
            indx = (x >= lo) & (x < hi)
        total += np.sum(indx)
    assert total == len(x)
    assert np.sum(x_hist) == len(x)

    # bins roughly equal population (equal-count binning)
    assert np.max(x_hist) - np.min(x_hist) <= 2

    # bin edges span the full data range
    assert x_lower[0] == np.min(x)
    assert x_upper[-1] == np.max(x)


def test_binmaker_matches_independent_histogram():
    rng = np.random.default_rng(2)
    x = rng.uniform(-5, 5, 2000)
    nbin = 8

    x_lower, x_mid, x_upper, x_hist = binmaker(x, nbin=nbin)

    # independently compute expected counts with np.histogram on the same
    # fixed-width edges
    edges = np.linspace(np.min(x), np.max(x), nbin + 1)
    expected_counts, _ = np.histogram(x, bins=edges)

    np.testing.assert_array_equal(x_hist, expected_counts)
    np.testing.assert_allclose(x_lower, edges[:-1])
    np.testing.assert_allclose(x_upper, edges[1:])
    np.testing.assert_allclose(x_mid, (edges[:-1] + edges[1:]) / 2.0)


def test_binmaker_log_steps():
    rng = np.random.default_rng(3)
    x = rng.uniform(1, 1000, 2000)
    nbin = 6

    x_lower, x_mid, x_upper, x_hist = binmaker(x, nbin=nbin, steptype="log")

    edges = np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), nbin + 1)
    expected_counts, _ = np.histogram(x, bins=edges)

    np.testing.assert_array_equal(x_hist, expected_counts)
    np.testing.assert_allclose(x_lower, edges[:-1])
    np.testing.assert_allclose(x_upper, edges[1:])


def test_roaming_nbinmaker_bins_are_overlapping_and_ordered():
    rng = np.random.default_rng(4)
    x = rng.uniform(0, 50, 500)

    x_lower, x_mid, x_upper, x_hist = roaming_nbinmaker(x, nbin=5, ntot=15)

    assert np.all(np.diff(x_lower) >= 0)
    assert np.all(x_upper >= x_lower)
    # every bin should contain a positive number of points (roaming bins by
    # construction always straddle at least one data point)
    assert np.all(x_hist > 0)


def test_roaming_binmaker_linear_matches_nonoverlap_case():
    # when ntot == nbin, roaming_binmaker reduces to fixed-width bins over
    # the same range as binmaker (minus the very last partial bin), so counts
    # should agree with an independent histogram computed the same way.
    rng = np.random.default_rng(5)
    x = rng.uniform(0, 20, 500)
    nbin = 5

    x_lower, x_mid, x_upper, x_hist = roaming_binmaker(x, nbin=nbin, ntot=nbin)

    dx = (np.max(x) - np.min(x)) / nbin
    for lo, hi, cnt in zip(x_lower, x_upper, x_hist):
        expected = np.sum((x >= lo) & (x < hi)) if hi != x_upper[-1] else np.sum(
            (x >= lo) & (x <= hi)
        )
        assert cnt == expected
        assert np.isclose(hi - lo, dx)


def test_power_law_distribution_recovers_slope():
    rng_state = np.random.default_rng(6)
    np.random.seed(6)  # power_law_distribution_function uses the global RNG

    alpha = -2.35  # Salpeter-like slope
    xmin, xmax = 1.0, 100.0
    n = 200000

    x = power_law_distribution_function(n, alpha, xmin, xmax)

    assert len(x) == n
    assert np.all(x >= xmin) and np.all(x <= xmax)

    # independently fit the slope by log-log linear regression on a
    # histogram of the samples (dN/dlogx ~ x^(alpha+1))
    logbins = np.logspace(np.log10(xmin), np.log10(xmax), 40)
    counts, edges = np.histogram(x, bins=logbins)
    widths = np.diff(edges)
    mids = np.sqrt(edges[:-1] * edges[1:])
    dn_dx = counts / widths

    nonzero = dn_dx > 0
    fit_alpha, _ = np.polyfit(np.log10(mids[nonzero]), np.log10(dn_dx[nonzero]), 1)

    assert abs(fit_alpha - alpha) < 0.1


def test_power_law_distribution_alpha_zero_is_uniform():
    np.random.seed(7)
    xmin, xmax = 2.0, 10.0
    x = power_law_distribution_function(50000, 0.0, xmin, xmax)

    assert np.all(x >= xmin) and np.all(x <= xmax)
    # uniform distribution: mean should be close to the midpoint
    assert abs(np.mean(x) - (xmin + xmax) / 2.0) < 0.05 * (xmax - xmin)


def test_distance_hand_computed():
    p1 = np.array([0.0, 0.0, 0.0])
    p2 = np.array([3.0, 4.0, 0.0])
    assert np.isclose(distance(p1, p2), 5.0)

    p3 = np.array([1.0, 2.0, 2.0])
    p4 = np.array([0.0, 0.0, 0.0])
    assert np.isclose(distance(p3, p4), 3.0)


def test_minimum_distance_hand_computed():
    # three points on a line: 0, 1, 10 -- nearest neighbour distances are
    # 1 (for point 0 and point 1) and 9 (for point 10)
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
        ]
    )
    mind = minimum_distance(pts)
    np.testing.assert_allclose(mind, [1.0, 1.0, 9.0])


def test_interpolate_linear():
    # line through (0,0) and (2,4) -> y = 2x
    val = interpolate((0.0, 0.0), (2.0, 4.0), x=1.0)
    assert np.isclose(val, 2.0)

    val_y = interpolate((0.0, 0.0), (2.0, 4.0), y=3.0)
    assert np.isclose(val_y, 1.5)


def test_mean_prof_matches_manual_binned_mean():
    rng = np.random.default_rng(8)
    x = rng.uniform(0, 10, 400)
    y = 2.0 * x + rng.normal(0, 0.01, 400)

    x_bin, y_bin, y_sig = mean_prof(x, y, nbin=5, bintype="fix")

    x_lower, x_mid, x_upper, x_hist = binmaker(x, nbin=5)
    # mean_prof now assigns points to bins via bin_index(), matching
    # binmaker/nbinmaker's own last-bin-inclusive convention: every bin is
    # left-closed/right-open except the last, which also includes the upper
    # edge (so the single max-value point lands in the last bin).
    for i in range(5):
        if i < 4:
            indx = (x >= x_lower[i]) & (x < x_upper[i])
        else:
            indx = (x >= x_lower[i]) & (x <= x_upper[i])
        if np.sum(indx) > 0:
            expected_mean = np.mean(y[indx])
            assert np.isclose(y_bin[i], expected_mean, atol=1e-6)


def test_smooth_returns_ordered_bins():
    rng = np.random.default_rng(9)
    x = np.sort(rng.uniform(0, 10, 300))
    y = np.sin(x)

    x_bin, y_bin, y_sig, y_min, y_max = smooth(x, y, dx=1.0)

    assert np.all(np.array(y_min) <= np.array(y_bin))
    assert np.all(np.array(y_bin) <= np.array(y_max))


def test_bin_index_matches_manual_membership_for_interior_and_last_bin():
    x_lower = np.array([0.0, 2.0, 5.0])
    x_upper = np.array([2.0, 5.0, 10.0])
    x = np.array([0.0, 1.9, 2.0, 4.9, 5.0, 9.9, 10.0])

    idx = bin_index(x, x_lower, x_upper)

    # interior bins are left-closed/right-open; the last bin is closed on
    # both ends (clamped), matching nbinmaker/binmaker's own convention
    expected = np.array([0, 0, 1, 1, 2, 2, 2])
    np.testing.assert_array_equal(idx, expected)


def test_bin_index_clamps_points_beyond_last_edge():
    x_lower = np.array([0.0, 1.0])
    x_upper = np.array([1.0, 2.0])
    # a point exactly at the global max upper edge should clamp to the last
    # bin index, not overflow to len(x_upper)
    idx = bin_index(np.array([2.0]), x_lower, x_upper)
    assert idx[0] == 1


def test_dx_function_recovers_known_power_law_slope():
    np.random.seed(30)
    alpha_true = -2.0
    x = power_law_distribution_function(50000, alpha_true, 1.0, 100.0)

    x_mean, x_hist, dx, alpha, ealpha, yalpha, eyalpha = dx_function(x, nx=15)

    # dx_function fits dN/dx directly, so the recovered slope should match
    # the input power-law exponent alpha_true (not alpha_true+1, since this
    # is the differential distribution, unlike mass_function's dN/dm binning
    # convention which is fit the same way but interpreted per-bin)
    assert abs(alpha - alpha_true) < 0.1


def test_dx_function_custom_bins_matches_default_bins():
    np.random.seed(31)
    x = power_law_distribution_function(20000, -2.0, 1.0, 50.0)

    x_lower, x_mean, x_upper, x_hist_default = nbinmaker(x, 10)

    x_mean1, x_hist1, dx1, alpha1, ea1, ya1, eya1 = dx_function(x, nx=10)
    x_mean2, x_hist2, dx2, alpha2, ea2, ya2, eya2 = dx_function(
        x, x_lower=x_lower, x_mean=x_mean, x_upper=x_upper
    )

    # supplying the exact same bins nbinmaker would have generated
    # internally should reproduce the same fit
    np.testing.assert_allclose(x_hist1, x_hist2)
    assert np.isclose(alpha1, alpha2, rtol=1e-6)


def test_tapered_dx_function_runs_and_fits_shape():
    rng = np.random.default_rng(32)

    def shape(x, alpha=-2.0, xc=3.0, beta=2.0):
        return (x**alpha) * (1.0 - np.exp(-((x / xc) ** beta)))

    grid = np.linspace(0.3, 30.0, 3000)
    f_max = np.max(shape(grid)) * 1.05
    samples = np.empty(0)
    while len(samples) < 15000:
        batch = 20000
        x_trial = rng.uniform(0.3, 30.0, batch)
        u = rng.uniform(0, f_max, batch)
        samples = np.append(samples, x_trial[u <= shape(x_trial)])
    x = samples[:15000]

    x_mean, x_hist, dx, A, eA, alpha, ealpha, xc, exc, beta, ebeta = tapered_dx_function(x, nx=15)

    assert np.isfinite(alpha)
    assert np.isfinite(xc)
    # the fitted curve should reproduce the binned data it was fit to
    from clustertools.util.recipes import tapered_func

    predicted = tapered_func(x_mean, A, alpha, xc, beta)
    valid = x_hist > 0
    residual = np.abs(np.log10(predicted[valid]) - np.log10(dx[valid]))
    assert np.median(residual) < 0.1


def test_x_hist_matches_nbinmaker_counts():
    rng = np.random.default_rng(33)
    x = rng.uniform(0, 100, 5000)

    x_mean_out, x_hist_out = x_hist(x, nx=12)
    x_lower, x_mid, x_upper, expected_hist = nbinmaker(x, 12)

    np.testing.assert_array_equal(x_hist_out, expected_hist)
    np.testing.assert_allclose(x_mean_out, x_mid)


def test_x_hist_custom_bins_matches_bincount():
    rng = np.random.default_rng(34)
    x = rng.uniform(0, 100, 5000)

    x_lower, x_mean_bins, x_upper, hist0 = nbinmaker(x, 8)
    x_mean_out, x_hist_out = x_hist(x, x_lower=x_lower, x_mean=x_mean_bins, x_upper=x_upper)

    idx = bin_index(x, x_lower, x_upper)
    expected = np.bincount(idx, minlength=len(x_lower))

    np.testing.assert_array_equal(x_hist_out, expected)
