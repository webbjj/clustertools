import numpy as np
import pytest
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
try:
    from galpy.util import coords
except ImportError:
    import galpy.util.bovy_coords as coords

import clustertools as ctools
from clustertools.analysis.orbits import (
    initialize_orbits,
    interpolate_orbit,
    interpolate_orbits,
    orbital_path,
    orbital_path_match,
)


def test_initialize_orbit_matches_independent_galpy_orbit():
    ro, vo = 8.275, 220.0
    xgc, ygc, zgc = 8000.0, 500.0, 100.0  # pc
    vxgc, vygc, vzgc = 10.0, 220.0, 5.0  # km/s

    cluster = ctools.StarCluster(units="pckms", origin="galaxy")
    cluster.add_stars(np.array([0.0]), np.array([0.0]), np.array([0.0]),
                       np.array([0.0]), np.array([0.0]), np.array([0.0]))
    cluster.add_orbit(xgc, ygc, zgc, vxgc, vygc, vzgc)
    cluster._ro, cluster._vo, cluster._zo = ro, vo, 0.0
    cluster._solarmotion = [-11.1, 12.24, 7.25]

    orbit = cluster.initialize_orbit()

    # independent cross-check: build a galpy Orbit directly from the same
    # cylindrical coordinates, without going through clustertools' own
    # initialize_orbit code path, and confirm they agree. galpy's Orbit
    # constructor expects the array input in natural (dimensionless) units
    # when ro/vo are also supplied -- i.e. R/ro, v/vo -- matching what
    # initialize_orbit does internally via cluster.to_galpy() before
    # constructing its Orbit.
    x_kpc, y_kpc, z_kpc = xgc / 1000.0, ygc / 1000.0, zgc / 1000.0
    R, phi, z = coords.rect_to_cyl(x_kpc, y_kpc, z_kpc)
    vR, vT, vz = coords.rect_to_cyl_vec(vxgc, vygc, vzgc, x_kpc, y_kpc, z_kpc)
    expected_orbit = Orbit([R / ro, vR / vo, vT / vo, z / ro, vz / vo, phi], ro=ro, vo=vo)

    assert np.isclose(orbit.x(), expected_orbit.x(), rtol=1e-8)
    assert np.isclose(orbit.y(), expected_orbit.y(), rtol=1e-8)
    assert np.isclose(orbit.z(), expected_orbit.z(), rtol=1e-8)
    assert np.isclose(orbit.vx(), expected_orbit.vx(), rtol=1e-8)
    assert np.isclose(orbit.vy(), expected_orbit.vy(), rtol=1e-8)
    assert np.isclose(orbit.vz(), expected_orbit.vz(), rtol=1e-8)


def test_initialize_orbit_position_matches_cluster_galactocentric_coords():
    cluster = ctools.StarCluster(units="kpckms", origin="galaxy")
    cluster.add_stars(np.array([0.0]), np.array([0.0]), np.array([0.0]),
                       np.array([0.0]), np.array([0.0]), np.array([0.0]))
    cluster.add_orbit(8.0, 0.5, 0.1, 10.0, 220.0, 5.0)

    orbit = cluster.initialize_orbit()

    assert np.isclose(orbit.x(), cluster.xgc, rtol=1e-6)
    assert np.isclose(orbit.y(), cluster.ygc, rtol=1e-6)
    assert np.isclose(orbit.z(), cluster.zgc, rtol=1e-6)
    assert np.isclose(orbit.vx(), cluster.vxgc, rtol=1e-6)
    assert np.isclose(orbit.vy(), cluster.vygc, rtol=1e-6)
    assert np.isclose(orbit.vz(), cluster.vzgc, rtol=1e-6)


def test_initialize_orbits_matches_per_star_galactocentric_positions():
    cluster = ctools.StarCluster(units="kpckms", origin="galaxy")
    rng = np.random.default_rng(800)
    n = 20
    x = 8.0 + rng.uniform(-0.5, 0.5, n)
    y = rng.uniform(-0.5, 0.5, n)
    z = rng.uniform(-0.1, 0.1, n)
    vx = rng.uniform(-20, 20, n)
    vy = 220.0 + rng.uniform(-20, 20, n)
    vz = rng.uniform(-10, 10, n)
    cluster.add_stars(x, y, z, vx, vy, vz, sortstars=True)
    cluster.add_orbit(8.0, 0.0, 0.0, 0.0, 220.0, 0.0)

    orbits = initialize_orbits(cluster)

    # each star's individually-initialized orbit should reproduce that
    # star's own galactocentric position/velocity at t=0 -- an independent
    # sanity check using the cluster's own recorded arrays, not a re-call
    # into clustertools' own orbit-construction internals
    np.testing.assert_allclose(orbits.x(), cluster.x, rtol=1e-6)
    np.testing.assert_allclose(orbits.y(), cluster.y, rtol=1e-6)
    np.testing.assert_allclose(orbits.z(), cluster.z, rtol=1e-6)
    np.testing.assert_allclose(orbits.vx(), cluster.vx, rtol=1e-6)
    np.testing.assert_allclose(orbits.vy(), cluster.vy, rtol=1e-6)
    np.testing.assert_allclose(orbits.vz(), cluster.vz, rtol=1e-6)


def test_interpolate_orbit_at_zero_time_matches_current_position():
    cluster = ctools.StarCluster(units="kpckms", origin="galaxy")
    cluster.add_stars(np.array([0.0]), np.array([0.0]), np.array([0.0]),
                       np.array([0.0]), np.array([0.0]), np.array([0.0]))
    cluster.add_orbit(8.0, 0.0, 0.0, 0.0, 220.0, 0.0)

    xgc, ygc, zgc, vxgc, vygc, vzgc = interpolate_orbit(cluster, pot=MWPotential2014, tfinal=0.0, nt=10)

    assert np.isclose(xgc, cluster.xgc, atol=1e-4)
    assert np.isclose(ygc, cluster.ygc, atol=1e-4)
    assert np.isclose(zgc, cluster.zgc, atol=1e-4)
    assert np.isclose(vxgc, cluster.vxgc, atol=1e-3)
    assert np.isclose(vygc, cluster.vygc, atol=1e-3)


def test_interpolate_orbits_at_zero_time_matches_current_star_positions():
    cluster = ctools.StarCluster(units="kpckms", origin="galaxy")
    rng = np.random.default_rng(801)
    n = 10
    x = 8.0 + rng.uniform(-0.3, 0.3, n)
    y = rng.uniform(-0.3, 0.3, n)
    z = rng.uniform(-0.05, 0.05, n)
    vx = rng.uniform(-10, 10, n)
    vy = 220.0 + rng.uniform(-10, 10, n)
    vz = rng.uniform(-5, 5, n)
    cluster.add_stars(x, y, z, vx, vy, vz, sortstars=True)
    cluster.add_orbit(8.0, 0.0, 0.0, 0.0, 220.0, 0.0)

    xs, ys, zs, vxs, vys, vzs = interpolate_orbits(cluster, pot=MWPotential2014, tfinal=0.0, nt=10)

    np.testing.assert_allclose(xs, cluster.x, atol=1e-4)
    np.testing.assert_allclose(ys, cluster.y, atol=1e-4)
    np.testing.assert_allclose(zs, cluster.z, atol=1e-4)


def test_orbital_path_at_t_zero_matches_cluster_galactocentric_position():
    cluster = ctools.StarCluster(units="kpckms", origin="galaxy")
    cluster.add_stars(np.array([0.0]), np.array([0.0]), np.array([0.0]),
                       np.array([0.0]), np.array([0.0]), np.array([0.0]))
    cluster.add_orbit(8.0, 0.0, 0.0, 0.0, 220.0, 0.0)

    nt = 5000
    t, x, y, z, vx, vy, vz = orbital_path(cluster, tfinal=0.05, nt=nt, pot=MWPotential2014)

    # the path should be time-ordered and bracket t=0 (it covers -tfinal to
    # +tfinal), and the point nearest t=0 should sit close to the cluster's
    # current galactocentric position. The grid is discrete, so the nearest
    # sample to t=0 is at most half a grid step away, not exactly at t=0.
    assert np.all(np.diff(t) >= 0)
    grid_step = 2 * 0.05 / (nt - 1)
    i0 = np.argmin(np.abs(t))
    assert np.abs(t[i0]) <= grid_step
    # position tolerance must account for how far the star can move (at up
    # to ~300 km/s ~ 300 kpc/Gyr) over half a grid step in time
    pos_tol = 300.0 * grid_step
    assert np.isclose(x[i0], cluster.xgc, atol=pos_tol)
    assert np.isclose(y[i0], cluster.ygc, atol=pos_tol)
    assert np.isclose(z[i0], cluster.zgc, atol=pos_tol)


def test_orbital_path_match_star_on_orbit_has_near_zero_path_distance():
    # a star sitting exactly at the cluster's own galactocentric position
    # and velocity is, by definition, on the orbital path -- its matched
    # distance to the path should be close to zero, an independent physical
    # expectation rather than a property of the matching algorithm itself
    cluster = ctools.StarCluster(units="kpckms", origin="galaxy")
    cluster.add_stars(np.array([8.0]), np.array([0.0]), np.array([0.0]),
                       np.array([0.0]), np.array([220.0]), np.array([0.0]))
    cluster.add_orbit(8.0, 0.0, 0.0, 0.0, 220.0, 0.0)

    tstar, dprog, dpath = orbital_path_match(cluster, tfinal=0.05, nt=5000, pot=MWPotential2014)

    assert np.abs(dpath[0]) < 0.05
