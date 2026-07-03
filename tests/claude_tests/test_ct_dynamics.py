import numpy as np
import pytest
from scipy.spatial import cKDTree
from galpy.potential import MWPotential2014, PlummerPotential, rtide

import clustertools as ctools
from clustertools.analysis.functions import (
    find_centre_of_mass,
    find_centre_of_density,
    rtidal,
    rlimiting,
    ckin,
    closest_star,
)


GRAV_PCKMS = 4.302e-3  # G in pc (km/s)^2 / Msun, independently looked up


def test_two_body_potential_energy_matches_closed_form(two_body_cluster):
    cluster = two_body_cluster
    m1, m2 = cluster.m[0], cluster.m[1]
    sep = abs(cluster.x[1] - cluster.x[0])

    kin, pot = cluster.energies(specific=False, full=True)

    # closed-form two-body potential energy per star: -G*m1*m2/r (each star
    # gets the full pairwise term in this code's convention, see
    # analysis/functions.py::_potential_energy)
    expected_pot = -GRAV_PCKMS * m1 * m2 / sep

    np.testing.assert_allclose(pot, [expected_pot, expected_pot], rtol=1e-6)
    # both stars at rest -> zero kinetic energy
    np.testing.assert_allclose(kin, [0.0, 0.0], atol=1e-12)


def test_energies_kinetic_matches_half_mv_squared():
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    x = np.array([0.0, 5.0, -3.0])
    y = np.array([0.0, 0.0, 0.0])
    z = np.array([0.0, 0.0, 0.0])
    vx = np.array([1.0, -2.0, 3.0])
    vy = np.array([0.0, 1.0, 0.0])
    vz = np.array([0.0, 0.0, 2.0])
    m = np.array([1.0, 2.0, 0.5])
    cluster.add_stars(x, y, z, vx, vy, vz, m)

    kin, pot = cluster.energies(specific=False, full=True)

    v2 = vx**2 + vy**2 + vz**2
    expected_kin = 0.5 * m * v2
    np.testing.assert_allclose(kin, expected_kin, rtol=1e-10)


def test_relaxation_time_matches_spitzer_formula_independently(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=400, rmax=8.0, vscale=1.5, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    trh_code = cluster.half_mass_relaxation_time(coulomb=0.4)

    # independently recompute the Spitzer (1987) half-mass relaxation time
    # formula directly from the definition, without calling clustertools'
    # own relaxation_time/half_mass_relaxation_time functions
    mass = np.sum(cluster.m)
    ntot = float(cluster.ntot)
    mbar = mass / ntot
    lnlambda = np.log(0.4 * ntot)
    rm = cluster.rm

    trh_expected = 0.138 * (mass**0.5) * (rm**1.5) / (mbar * np.sqrt(GRAV_PCKMS) * lnlambda)
    trh_expected *= 3.086e13 / (3600.0 * 24.0 * 365.0 * 1e6)  # to Myr

    assert np.isclose(trh_code, trh_expected, rtol=1e-6)


def test_rlagrange_matches_percentile_definition_for_equal_mass(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=1000, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    nlagrange = 10
    rn = cluster.rlagrange(nlagrange=nlagrange)

    # independent computation from the definition: rn[k] is the radius
    # enclosing (k+1)/nlagrange of the total mass, for equal-mass stars this
    # is just the corresponding order statistic of the sorted radii
    r_sorted = np.sort(cluster.r)
    n = len(r_sorted)
    for k in range(nlagrange - 1):
        frac = (k + 1) / nlagrange
        idx = int(np.ceil(frac * n)) - 1
        idx = min(idx, n - 1)
        assert np.isclose(rn[k], r_sorted[idx], rtol=1e-6)

    assert np.isclose(rn[-1], np.max(cluster.r))
    assert np.all(np.diff(rn) >= -1e-9)


def test_virial_radius_is_positive_and_order_of_magnitude_sane(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=500, rmax=8.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    rv = cluster.virial_radius(method="inverse_distance")

    assert rv > 0
    # virial radius should be within an order of magnitude of the half-mass
    # radius for a roughly uniform-ish population -- a coarse sanity bound,
    # not a precise physical prediction
    assert 0.1 * cluster.rm < rv < 10.0 * cluster.rm


def test_virialize_produces_requested_virial_ratio(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=400, rmax=8.0, vscale=3.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    qvir_target = 0.5
    cluster.virialize(qvir=qvir_target)

    # independent cross-check: recompute energies() from scratch (a
    # different code path than virialize's internal qv calculation) and
    # verify the resulting virial ratio matches the target
    kin, pot = cluster.energies(specific=True, full=True)
    ektot = np.sum(kin)
    ptot = np.sum(pot) / 2.0
    qvir_actual = abs(ektot / ptot)

    assert np.isclose(qvir_actual, qvir_target, rtol=1e-6)


def test_rcore_less_than_half_mass_radius(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=800, rmax=8.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    rc = cluster.rcore()
    assert rc > 0
    assert rc < cluster.rm


def test_core_relaxation_time_matches_stone_ostriker_formula_independently(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=800, rmax=8.0, vscale=1.5, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()

    trc_code = cluster.core_relaxation_time(coulomb=0.4)

    # independently recompute the Stone & Ostriker (2015) core relaxation
    # time formula from its definition, using cluster.rcore()/cluster.rm as
    # already-independently-tested primitives rather than re-deriving them
    mtot = np.sum(cluster.m)
    mbar = np.mean(cluster.m)
    lnlambda = np.log(0.4 * cluster.ntot)
    rc = cluster.rcore()
    rh = cluster.rm

    trc_expected = (0.39 / lnlambda) * np.sqrt(rc**3.0 / (GRAV_PCKMS * mtot)) * (mtot / mbar) * np.sqrt(rc * rh) / (rc + rh)
    trc_expected *= 3.086e13 / (3600.0 * 24.0 * 365.0 * 1e6)  # to Myr

    assert np.isclose(trc_code, trc_expected, rtol=1e-6)


def test_find_centre_of_mass_matches_manual_mass_weighted_centroid(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=500, rmax=8.0, mass=1.0)
    rng = np.random.default_rng(600)
    cluster.m = rng.uniform(0.1, 5.0, cluster.ntot)

    xc, yc, zc, vxc, vyc, vzc = find_centre_of_mass(cluster)

    expected_xc = np.sum(cluster.m * cluster.x) / np.sum(cluster.m)
    expected_yc = np.sum(cluster.m * cluster.y) / np.sum(cluster.m)
    expected_zc = np.sum(cluster.m * cluster.z) / np.sum(cluster.m)
    expected_vxc = np.sum(cluster.m * cluster.vx) / np.sum(cluster.m)

    assert np.isclose(xc, expected_xc)
    assert np.isclose(yc, expected_yc)
    assert np.isclose(zc, expected_zc)
    assert np.isclose(vxc, expected_vxc)


def test_find_centre_of_density_recovers_known_offset_gaussian_peak():
    # build a cluster whose density is sharply peaked around a known offset
    # (x0,y0,z0), well-separated from the coordinate origin, so
    # find_centre_of_density's result can be checked against that known
    # peak location rather than against find_centre_of_mass (which would
    # trivially agree for a symmetric distribution and wouldn't test
    # anything density-specific)
    rng = np.random.default_rng(601)
    n = 5000
    x0, y0, z0 = 5.0, -3.0, 2.0
    x = rng.normal(x0, 0.5, n)
    y = rng.normal(y0, 0.5, n)
    z = rng.normal(z0, 0.5, n)
    vx = rng.normal(0, 1, n)
    vy = rng.normal(0, 1, n)
    vz = rng.normal(0, 1, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, sortstars=True)

    xc, yc, zc, vxc, vyc, vzc = find_centre_of_density(cluster)

    assert np.isclose(xc, x0, atol=0.2)
    assert np.isclose(yc, y0, atol=0.2)
    assert np.isclose(zc, z0, atol=0.2)


def test_closest_star_matches_independent_kdtree_query(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=600, rmax=10.0, mass=1.0)

    min_dist = cluster.closest_star()

    # independent nearest-neighbour computation via scipy's cKDTree (a
    # different algorithm/library than clustertools' own numba-based
    # brute-force minimum_distance)
    pos = np.column_stack([cluster.x, cluster.y, cluster.z])
    tree = cKDTree(pos)
    dist, _ = tree.query(pos, k=2)  # k=1 is the point itself (distance 0)
    expected = dist[:, 1]

    np.testing.assert_allclose(np.sort(min_dist), np.sort(expected), rtol=1e-6)


def test_rtidal_matches_independent_galpy_rtide_call(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=500, rmax=5.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()
    cluster.add_orbit(8000.0, 0.0, 0.0, 0.0, 220.0, 0.0)

    pot = MWPotential2014

    rt_code = cluster.rtidal(pot=pot)

    # rtidal() delegates directly to galpy.potential.rtide(pot, R, z, M=...)
    # in galpy natural units -- reproduce that call independently using the
    # cluster's own recorded galactocentric position/mass, converted the
    # same way, without going through clustertools' internal unit-conversion
    # machinery for the comparison value itself
    cluster.save_cluster()
    units0, origin0, rorder0, rorder_origin0 = cluster.units0, cluster.origin0, cluster.rorder0, cluster.rorder_origin0
    cluster.to_galpy()
    R = np.sqrt(cluster.xgc**2 + cluster.ygc**2)
    z = cluster.zgc
    M = np.sum(cluster.m)
    cluster.return_cluster(units0, origin0, rorder0, rorder_origin0)
    rt_expected_galpy_units = rtide(pot, R, z, M=M, use_physical=False)

    ro = cluster._ro
    rt_expected = rt_expected_galpy_units * ro * 1000.0  # kpc -> pc

    assert np.isclose(rt_code, rt_expected, rtol=1e-3)


def test_rlimiting_is_positive_and_less_than_rtidal(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=500, rmax=5.0, mass=1.0)
    cluster.find_centre()
    cluster.to_centre()
    cluster.add_orbit(8000.0, 0.0, 0.0, 0.0, 220.0, 0.0)

    pot = MWPotential2014

    rl = cluster.rlimiting(pot=pot)
    rt = cluster.rtidal(pot=pot)

    assert rl > 0
    # the limiting radius (where cluster density drops to the local
    # background density) should not exceed the tidal radius
    assert rl <= rt * 1.01


def test_ckin_runs_and_returns_finite_value():
    from clustertools.util.recipes import power_law_distribution_function

    np.random.seed(700)
    n = 6000
    m = power_law_distribution_function(n, -2.0, 0.2, 5.0)

    rng = np.random.default_rng(701)
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

    ck = ckin(cluster, nmass=10)
    assert np.isfinite(ck)
