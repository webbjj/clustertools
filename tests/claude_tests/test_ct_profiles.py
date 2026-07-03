import numpy as np
import pytest

from clustertools.analysis.profiles import (
    rho_prof,
    m_prof,
    alpha_prof,
    sigv_prof,
    beta_prof,
    v_prof,
    v2_prof,
    eta_prof,
    meq_prof,
    vcirc_prof,
)
from clustertools.analysis.functions import mass_function as mass_function_fn
from clustertools.util.coordinates import sphere_coords
from clustertools.util.recipes import nbinmaker, bin_index
from clustertools.util.constants import _get_grav


def test_rho_prof_and_m_prof_are_mutually_consistent(equal_mass_cluster_factory):
    # rho_prof (density) and m_prof (mass) both bin by equal star-count using
    # the same default binning algorithm -- cross-check them against each
    # other: density * shell volume should recover the binned mass, an
    # invariant that doesn't require knowing either function's "correct"
    # numeric answer in advance.
    cluster = equal_mass_cluster_factory(n=2000, rmax=10.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    nrad = 12
    rprof, pprof, nprof_rho = rho_prof(cluster, nrad=nrad, projected=False)
    rprof_m, mprof, nprof_m = m_prof(cluster, nrad=nrad, projected=False, cumulative=False)

    np.testing.assert_allclose(rprof, rprof_m)
    np.testing.assert_array_equal(nprof_rho, nprof_m)

    # reconstruct bin edges the same way rho_prof/m_prof do internally
    # (equal-count bins on the "in-bounds" subset), then verify pprof*vol == mprof
    from clustertools.util.recipes import nbinmaker

    r_lower, r_mean, r_upper, r_hist = nbinmaker(cluster.r, nrad)
    vol = (4.0 / 3.0) * np.pi * (r_upper**3 - r_lower**3)

    np.testing.assert_allclose(pprof * vol, mprof, rtol=1e-8)


def test_m_prof_cumulative_matches_manual_cumulative_sum(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=1500, rmax=10.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    nrad = 10
    rprof, mprof_noncum, nprof = m_prof(cluster, nrad=nrad, cumulative=False)
    rprof2, mprof_cum, nprof2 = m_prof(cluster, nrad=nrad, cumulative=True)

    # independent check: cumulative mass at bin i should equal the sum of
    # non-cumulative bin masses up to and including bin i (since bins are
    # contiguous, non-overlapping, and cover [rmin, r_upper[i]])
    manual_cumsum = np.cumsum(mprof_noncum)
    np.testing.assert_allclose(mprof_cum, manual_cumsum, rtol=1e-8)

    # cumulative mass should be monotonically non-decreasing
    assert np.all(np.diff(mprof_cum) >= -1e-9)

    # m_prof now assigns points to bins via bin_index(), which clamps the
    # last bin to include the global-max-radius point -- so the cumulative
    # total in the final bin should reach the cluster's full mass exactly,
    # rather than silently excluding the single outermost star.
    assert np.isclose(mprof_cum[-1], cluster.mtot, rtol=1e-6)


def test_rho_prof_normalize_scales_radius_by_rm(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=1000, rmax=10.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()
    cluster.analyze()
    rm = cluster.rm

    rprof_raw, pprof_raw, nprof_raw = rho_prof(cluster, nrad=8, normalize=False)
    rprof_norm, pprof_norm, nprof_norm = rho_prof(cluster, nrad=8, normalize=True)

    np.testing.assert_allclose(rprof_norm, rprof_raw / rm, rtol=1e-6)


def test_rho_prof_density_decreases_outward_for_plummer(plummer_cluster):
    cluster = plummer_cluster
    rprof, pprof, nprof = rho_prof(cluster, nrad=10)

    # a Plummer sphere's 3D density profile is strictly decreasing with
    # radius -- an independent physical expectation of the model, not a
    # property of this specific binning code
    assert np.all(np.diff(pprof) <= 0)


def test_sigv_prof_matches_manual_std_per_bin(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=3000, rmax=10.0, vscale=2.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    nrad = 5
    rprofn, sigvprof = sigv_prof(cluster, nrad=nrad)

    r_lower, r_mean, r_upper, r_hist = nbinmaker(cluster.r, nrad)
    bin_idx = bin_index(cluster.r, r_lower, r_upper)

    for i in range(nrad):
        mask = bin_idx == i
        if np.sum(mask) > 3:
            expected = np.std(cluster.v[mask])
            assert np.isclose(sigvprof[i], expected, rtol=1e-8)


def test_beta_prof_matches_manual_spherical_velocity_dispersions(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=3000, rmax=10.0, vscale=2.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    nrad = 5
    rprofn, betaprof = beta_prof(cluster, nrad=nrad)

    # independent recomputation using the same spherical-coordinate
    # primitive (sphere_coords) but a separately-written binning+formula
    # pass, rather than calling beta_prof's own internals
    r, phi, theta, vr, vphi, vtheta = sphere_coords(cluster)
    r_lower, r_mean, r_upper, r_hist = nbinmaker(cluster.r, nrad)
    bin_idx = bin_index(cluster.r, r_lower, r_upper)

    for i in range(nrad):
        mask = bin_idx == i
        if np.sum(mask) > 3:
            sigr = np.std(vr[mask])
            sigp = np.std(vphi[mask])
            sigt = np.std(vtheta[mask])
            expected_beta = 1.0 - (sigt**2.0 + sigp**2.0) / (2.0 * sigr**2.0)
            assert np.isclose(betaprof[i], expected_beta, rtol=1e-6)


def test_beta_prof_is_near_zero_for_isotropic_velocities():
    # an isotropic (independent-Gaussian-per-Cartesian-component) velocity
    # distribution has no preferred radial/tangential direction, so beta
    # should be close to zero everywhere -- a physical expectation
    # independent of beta_prof's own implementation
    import clustertools as ctools

    rng = np.random.default_rng(77)
    n = 20000
    x = rng.uniform(-10, 10, n)
    y = rng.uniform(-10, 10, n)
    z = rng.uniform(-10, 10, n)
    vx = rng.normal(0, 5, n)
    vy = rng.normal(0, 5, n)
    vz = rng.normal(0, 5, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, sortstars=True)
    cluster.find_centre()
    cluster.to_centre()

    rprofn, betaprof = beta_prof(cluster, nrad=5)

    assert np.all(np.abs(betaprof) < 0.1)


def test_v_prof_matches_manual_mean_speed_per_bin(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=3000, rmax=10.0, vscale=2.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    nrad = 5
    rprofn, vprof = v_prof(cluster, nrad=nrad)

    r_lower, r_mean, r_upper, r_hist = nbinmaker(cluster.r, nrad)
    bin_idx = bin_index(cluster.r, r_lower, r_upper)

    for i in range(nrad):
        mask = bin_idx == i
        if np.sum(mask) > 3:
            expected = np.mean(cluster.v[mask])
            assert np.isclose(vprof[i], expected, rtol=1e-8)


def test_v2_prof_matches_manual_mean_square_speed_per_bin(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=3000, rmax=10.0, vscale=2.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    nrad = 5
    rprofn, v2prof = v2_prof(cluster, nrad=nrad)

    r_lower, r_mean, r_upper, r_hist = nbinmaker(cluster.r, nrad)
    bin_idx = bin_index(cluster.r, r_lower, r_upper)

    for i in range(nrad):
        mask = bin_idx == i
        if np.sum(mask) > 3:
            expected = np.mean(cluster.v[mask] ** 2.0)
            assert np.isclose(v2prof[i], expected, rtol=1e-8)


def test_vcirc_prof_vmax_matches_independent_formula(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=2000, rmax=10.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()
    cluster.analyze()

    rprof, vcprof, rvmax, vmax = vcirc_prof(cluster)

    # independent recomputation of the circular velocity curve from the
    # standard formula v_circ(r) = sqrt(G*M(<r)/r), using cumulative mass
    # along stars sorted by radius, without calling vcirc_prof's internals
    grav = _get_grav(cluster)
    order = cluster.rorder
    r_sorted = cluster.r[order]
    m_sorted = cluster.m[order]
    msum = np.cumsum(m_sorted)
    vcirc = np.sqrt(grav * msum / r_sorted)

    expected_vmax = np.amax(vcirc)
    expected_rvmax = r_sorted[np.argmax(vcirc)]

    assert np.isclose(vmax, expected_vmax, rtol=1e-8)
    assert np.isclose(rvmax, expected_rvmax, rtol=1e-8)


def test_alpha_prof_bin_matches_independent_mass_function_call():
    # alpha_prof fits a local mass-function slope within each radial bin by
    # calling mass_function() internally with a per-bin star mask -- verify
    # one bin's slope by independently reconstructing that same mask and
    # calling mass_function() directly ourselves.
    import clustertools as ctools
    from clustertools.util.recipes import power_law_distribution_function

    np.random.seed(500)
    n = 9000
    m = power_law_distribution_function(n, -2.3, 0.2, 5.0)

    rng = np.random.default_rng(501)
    x = rng.uniform(-10, 10, n)
    y = rng.uniform(-10, 10, n)
    z = rng.uniform(-10, 10, n)
    vx = rng.normal(0, 1, n)
    vy = rng.normal(0, 1, n)
    vz = rng.normal(0, 1, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, m, sortstars=True)
    cluster.find_centre()
    cluster.to_centre()

    nrad = 3
    nmass = 10
    rprofn, aprof, dalpha, edalpha, ydalpha, eydalpha = alpha_prof(cluster, nrad=nrad, nmass=nmass)

    base_indx = cluster.subset(kwmin=0, kwmax=1)
    r_lower, r_mean, r_upper, r_hist = nbinmaker(cluster.r[base_indx], nrad)
    bin_idx = bin_index(cluster.r, r_lower, r_upper)

    i = 1  # middle bin
    rindx = base_indx & (bin_idx == i)
    m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha = mass_function_fn(
        cluster, nmass=nmass, indx=rindx, projected=False, mcorr=None, plot=False
    )

    assert np.isclose(aprof[i], alpha, rtol=1e-6)


def test_eta_prof_and_meq_prof_run_and_return_finite_profiles(equal_mass_cluster_factory):
    from clustertools.util.recipes import power_law_distribution_function

    np.random.seed(502)
    n = 9000
    m = power_law_distribution_function(n, -2.3, 0.2, 5.0)

    import clustertools as ctools

    rng = np.random.default_rng(503)
    x = rng.uniform(-10, 10, n)
    y = rng.uniform(-10, 10, n)
    z = rng.uniform(-10, 10, n)
    vx = rng.normal(0, 1, n)
    vy = rng.normal(0, 1, n)
    vz = rng.normal(0, 1, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, m, sortstars=True)
    cluster.find_centre()
    cluster.to_centre()

    nrad = 3
    rprofn_eta, eprof, deta, edeta, ydeta, eydeta = eta_prof(cluster, nrad=nrad, nmass=10)
    assert len(rprofn_eta) == nrad
    assert np.all(np.isfinite(eprof))
    assert np.all(np.diff(rprofn_eta) > 0)

    rprofn_meq, meqprof, dmeq, edmeq, ydmeq, eymeq = meq_prof(cluster, nrad=nrad, nmass=10)
    assert len(rprofn_meq) == nrad
    assert np.all(np.isfinite(meqprof))
    assert np.all(np.diff(rprofn_meq) > 0)
