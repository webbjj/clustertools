import numpy as np
import pytest
from galpy.potential import MWPotential2014

import clustertools as ctools
from clustertools.tidaltail.tails import to_tail, tail_path, tail_path_match


def test_to_tail_preserves_position_and_velocity_norms():
    # to_tail rotates every star's position/velocity into a frame where the
    # cluster's bulk velocity points along +x. Rotations are norm-preserving,
    # so ||position|| and ||velocity|| for every star must be unchanged --
    # an invariant independent of the rotation matrix's own implementation.
    rng = np.random.default_rng(300)
    n = 100
    x = rng.uniform(-5, 5, n)
    y = rng.uniform(-5, 5, n)
    z = rng.uniform(-5, 5, n)
    vx = rng.uniform(-2, 2, n)
    vy = rng.uniform(-2, 2, n)
    vz = rng.uniform(-2, 2, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz)
    cluster.add_orbit(8000.0, 0.0, 0.0, 50.0, 200.0, 10.0)

    r_before = np.sqrt(cluster.x**2 + cluster.y**2 + cluster.z**2)
    v_before = np.sqrt(cluster.vx**2 + cluster.vy**2 + cluster.vz**2)

    x_tail, y_tail, z_tail, vx_tail, vy_tail, vz_tail = to_tail(cluster)

    r_after = np.sqrt(x_tail**2 + y_tail**2 + z_tail**2)
    v_after = np.sqrt(vx_tail**2 + vy_tail**2 + vz_tail**2)

    np.testing.assert_allclose(r_after, r_before, rtol=1e-8)
    np.testing.assert_allclose(v_after, v_before, rtol=1e-8)


def test_to_tail_does_not_mutate_cluster_coordinates():
    rng = np.random.default_rng(301)
    n = 50
    x = rng.uniform(-5, 5, n)
    y = rng.uniform(-5, 5, n)
    z = rng.uniform(-5, 5, n)
    vx = rng.uniform(-2, 2, n)
    vy = rng.uniform(-2, 2, n)
    vz = rng.uniform(-2, 2, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz)
    cluster.add_orbit(8000.0, 0.0, 0.0, 50.0, 200.0, 10.0)

    x0, y0, z0 = cluster.x.copy(), cluster.y.copy(), cluster.z.copy()

    to_tail(cluster)

    np.testing.assert_array_equal(cluster.x, x0)
    np.testing.assert_array_equal(cluster.y, y0)
    np.testing.assert_array_equal(cluster.z, z0)


def _make_orbit_hugging_cluster(n=400, spread=0.02, seed=900):
    # stars scattered tightly around the cluster's own galactocentric
    # position/velocity, so the tail path (which averages nearby stars per
    # orbital-phase bin) should stay close to the underlying orbit itself.
    # origin="galaxy" since star positions below are absolute galactocentric
    # coordinates, not offsets from a cluster centre.
    rng = np.random.default_rng(seed)
    x = 8000.0 + rng.uniform(-spread, spread, n) * 1000.0
    y = rng.uniform(-spread, spread, n) * 1000.0
    z = rng.uniform(-spread, spread, n) * 1000.0
    vx = rng.uniform(-1, 1, n)
    vy = 220.0 + rng.uniform(-1, 1, n)
    vz = rng.uniform(-1, 1, n)

    cluster = ctools.StarCluster(units="pckms", origin="galaxy")
    cluster.add_stars(x, y, z, vx, vy, vz, sortstars=True)
    cluster.add_orbit(8000.0, 0.0, 0.0, 0.0, 220.0, 0.0)
    return cluster


def test_tail_path_near_t_zero_matches_cluster_galactocentric_position():
    cluster = _make_orbit_hugging_cluster()

    ttail, xtail, ytail, ztail, vxtail, vytail, vztail = tail_path(
        cluster, tfinal=0.02, no=2000, nt=40, ntail=40, pot=MWPotential2014, bintype="num"
    )

    assert np.all(np.diff(ttail) >= -1e-9)
    i0 = np.argmin(np.abs(ttail))

    # tail path positions are returned in cluster.units (pckms here); the
    # bin nearest t=0 should be close to the cluster's own galactocentric
    # position, since all stars were scattered tightly around it
    assert np.isclose(xtail[i0], cluster.xgc, atol=50.0)
    assert np.isclose(ytail[i0], cluster.ygc, atol=50.0)
    assert np.isclose(ztail[i0], cluster.zgc, atol=50.0)


def test_tail_path_dmax_filters_distant_stars():
    cluster = _make_orbit_hugging_cluster(n=300, spread=0.02, seed=901)

    # add a handful of stars far off the orbital path -- with a tight dmax
    # they should be excluded from the tail-path averaging, changing the
    # resulting path compared to not filtering
    n_far = 20
    rng = np.random.default_rng(902)
    x_far = rng.uniform(-500, 500, n_far)
    y_far = 5000.0 + rng.uniform(-500, 500, n_far)
    z_far = rng.uniform(-500, 500, n_far)
    vx_far = rng.uniform(-5, 5, n_far)
    vy_far = rng.uniform(-5, 5, n_far)
    vz_far = rng.uniform(-5, 5, n_far)
    cluster.add_stars(x_far, y_far, z_far, vx_far, vy_far, vz_far, sortstars=True)

    ttail_all, *_rest_all = tail_path(cluster, tfinal=0.02, no=2000, nt=40, ntail=40, pot=MWPotential2014, bintype="num")
    ttail_filt, xtail_filt, *_rest_filt = tail_path(
        cluster, tfinal=0.02, no=2000, nt=40, ntail=40, pot=MWPotential2014, bintype="num", dmax=1.0
    )

    # filtering should not increase the number of populated bins
    assert len(ttail_filt) <= len(ttail_all)
    # and every returned bin position should still be finite
    assert np.all(np.isfinite(xtail_filt))


def test_tail_path_match_consistent_with_precomputed_path():
    cluster = _make_orbit_hugging_cluster(n=200, spread=0.02, seed=903)

    # tail_path_match(path=None) internally calls tail_path(...) without
    # forwarding a bintype kwarg, so it always uses tail_path's own default
    # ('fix') regardless of what's passed to tail_path_match -- match that
    # here so both calls are actually directly comparable.
    path = tail_path(cluster, tfinal=0.02, no=2000, nt=40, ntail=40, pot=MWPotential2014)

    tstar1, dprog1, dpath1 = tail_path_match(
        cluster, tfinal=0.02, no=2000, nt=40, ntail=40, pot=MWPotential2014, path=path
    )

    # calling with an explicitly-precomputed path should give the same
    # result as tail_path_match recomputing that identical path internally
    tstar2, dprog2, dpath2 = tail_path_match(
        cluster, tfinal=0.02, no=2000, nt=40, ntail=40, pot=MWPotential2014
    )

    np.testing.assert_allclose(tstar1, tstar2, rtol=1e-6)
    np.testing.assert_allclose(dpath1, dpath2, rtol=1e-6)
