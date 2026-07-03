import numpy as np
import pytest

import clustertools as ctools
from clustertools.cluster.operations import from_radec


def _make_cluster_with_orbit():
    rng = np.random.default_rng(40)
    n = 200
    x = rng.uniform(-5, 5, n)
    y = rng.uniform(-5, 5, n)
    z = rng.uniform(-5, 5, n)
    vx = rng.uniform(-1, 1, n)
    vy = rng.uniform(-1, 1, n)
    vz = rng.uniform(-1, 1, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, sortstars=True)
    cluster.add_orbit(8275.0, 0.0, 20.0, 10.0, 220.0, 5.0)
    cluster.find_centre()
    return cluster


def test_to_centre_and_back_to_cluster_round_trip():
    cluster = _make_cluster_with_orbit()
    x0, y0, z0 = cluster.x.copy(), cluster.y.copy(), cluster.z.copy()

    cluster.to_centre()
    assert cluster.origin == "centre"
    cluster.to_cluster()
    assert cluster.origin == "cluster"

    np.testing.assert_allclose(cluster.x, x0, atol=1e-9)
    np.testing.assert_allclose(cluster.y, y0, atol=1e-9)
    np.testing.assert_allclose(cluster.z, z0, atol=1e-9)


def test_to_centre_shift_matches_independent_offset():
    cluster = _make_cluster_with_orbit()
    x0 = cluster.x.copy()
    xc = cluster.xc

    cluster.to_centre()

    # independent check: to_centre should have subtracted exactly xc (the
    # centre-of-density/mass offset found by find_centre), nothing more
    np.testing.assert_allclose(cluster.x, x0 - xc, atol=1e-9)


def test_to_galaxy_and_back_to_cluster_round_trip():
    cluster = _make_cluster_with_orbit()
    x0, y0, z0 = cluster.x.copy(), cluster.y.copy(), cluster.z.copy()

    cluster.to_galaxy()
    assert cluster.origin == "galaxy"
    # independent check: galaxy origin should be cluster positions + xgc offset
    np.testing.assert_allclose(cluster.x, x0 + cluster.xgc, atol=1e-6)

    cluster.to_cluster()
    assert cluster.origin == "cluster"
    np.testing.assert_allclose(cluster.x, x0, atol=1e-6)
    np.testing.assert_allclose(cluster.y, y0, atol=1e-6)
    np.testing.assert_allclose(cluster.z, z0, atol=1e-6)


def test_to_radec_from_radec_round_trip():
    cluster = _make_cluster_with_orbit()
    cluster.to_galaxy()
    x0, y0, z0 = cluster.x.copy(), cluster.y.copy(), cluster.z.copy()
    vx0, vy0, vz0 = cluster.vx.copy(), cluster.vy.copy(), cluster.vz.copy()

    cluster.to_radec()
    assert cluster.units == "radec"
    from_radec(cluster)
    # from_radec always returns galactocentric positions in kpckms
    # (galpy's native physical units), regardless of the units the cluster
    # started in -- convert back to pckms before comparing.
    assert cluster.units == "kpckms"
    cluster.to_pckms()

    np.testing.assert_allclose(cluster.x, x0, rtol=1e-4, atol=1e-6)
    np.testing.assert_allclose(cluster.y, y0, rtol=1e-4, atol=1e-6)
    np.testing.assert_allclose(cluster.z, z0, rtol=1e-4, atol=1e-6)
    np.testing.assert_allclose(cluster.vx, vx0, rtol=1e-3, atol=1e-6)
    np.testing.assert_allclose(cluster.vy, vy0, rtol=1e-3, atol=1e-6)
    np.testing.assert_allclose(cluster.vz, vz0, rtol=1e-3, atol=1e-6)


def test_save_and_return_cluster_restores_state():
    cluster = _make_cluster_with_orbit()
    cluster.save_cluster()
    units0, origin0, rorder0, rorder_origin0 = (
        cluster.units0, cluster.origin0, cluster.rorder0, cluster.rorder_origin0
    )

    cluster.to_centre()
    cluster.to_kpckms()
    assert cluster.units != units0 or cluster.origin != origin0

    cluster.return_cluster(units0, origin0, rorder0, rorder_origin0)

    assert cluster.units == units0
    assert cluster.origin == origin0


def test_sub_cluster_restores_parent_state():
    cluster = _make_cluster_with_orbit()
    units0, origin0 = cluster.units, cluster.origin

    sub = ctools.sub_cluster(cluster, rmax=3.0)

    # parent cluster should be unchanged after extracting a subset
    assert cluster.units == units0
    assert cluster.origin == origin0

    # the returned subset should independently satisfy the radius constraint
    assert np.all(sub.r <= 3.0)


def test_sub_cluster_mass_constraint_matches_manual_mask():
    cluster = _make_cluster_with_orbit()
    rng = np.random.default_rng(41)
    cluster.m = rng.uniform(0.1, 5.0, cluster.ntot)

    mmin, mmax = 1.0, 3.0
    sub = ctools.sub_cluster(cluster, mmin=mmin, mmax=mmax, sortstars=False)

    assert np.all(sub.m >= mmin)
    assert np.all(sub.m <= mmax)
    # independent count check via a manual boolean mask on the parent arrays
    expected_n = np.sum((cluster.m >= mmin) & (cluster.m <= mmax))
    assert sub.ntot == expected_n
