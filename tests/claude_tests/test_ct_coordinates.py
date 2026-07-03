import numpy as np

from clustertools.util.coordinates import (
    cart_to_sphere,
    sphere_to_cart,
    cart_to_cyl,
    cyl_to_cart,
    cart_to_sky,
)


def test_cart_sphere_round_trip():
    rng = np.random.default_rng(10)
    n = 500
    x = rng.uniform(-10, 10, n)
    y = rng.uniform(-10, 10, n)
    z = rng.uniform(-10, 10, n)
    vx = rng.uniform(-5, 5, n)
    vy = rng.uniform(-5, 5, n)
    vz = rng.uniform(-5, 5, n)

    r, phi, theta, vr, vphi, vtheta = cart_to_sphere(x, y, z, vx, vy, vz)
    x2, y2, z2, vx2, vy2, vz2 = sphere_to_cart(r, phi, theta, vr, vphi, vtheta)

    np.testing.assert_allclose(x2, x, atol=1e-9)
    np.testing.assert_allclose(y2, y, atol=1e-9)
    np.testing.assert_allclose(z2, z, atol=1e-9)
    np.testing.assert_allclose(vx2, vx, atol=1e-9)
    np.testing.assert_allclose(vy2, vy, atol=1e-9)
    np.testing.assert_allclose(vz2, vz, atol=1e-9)


def test_cart_sphere_known_axis_points():
    # point on the +x axis
    r, phi, theta, vr, vphi, vtheta = cart_to_sphere(
        np.array([1.0]), np.array([0.0]), np.array([0.0]),
        np.array([0.0]), np.array([0.0]), np.array([0.0]),
    )
    assert np.isclose(r[0], 1.0)
    assert np.isclose(phi[0], 0.0, atol=1e-9)
    assert np.isclose(theta[0], np.pi / 2.0)

    # point on the +z axis
    r, phi, theta, vr, vphi, vtheta = cart_to_sphere(
        np.array([0.0]), np.array([0.0]), np.array([2.0]),
        np.array([0.0]), np.array([0.0]), np.array([0.0]),
    )
    assert np.isclose(r[0], 2.0)
    assert np.isclose(theta[0], 0.0, atol=1e-9)


def test_cart_cyl_round_trip():
    rng = np.random.default_rng(11)
    n = 500
    x = rng.uniform(-10, 10, n)
    y = rng.uniform(-10, 10, n)
    z = rng.uniform(-10, 10, n)
    vx = rng.uniform(-5, 5, n)
    vy = rng.uniform(-5, 5, n)
    vz = rng.uniform(-5, 5, n)

    r, theta, zed, vr, vtheta, vzed = cart_to_cyl(x, y, z, vx, vy, vz)
    x2, y2, z2, vx2, vy2, vz2 = cyl_to_cart(r, theta, zed, vr, vtheta, vzed)

    np.testing.assert_allclose(x2, x, atol=1e-9)
    np.testing.assert_allclose(y2, y, atol=1e-9)
    np.testing.assert_allclose(z2, z, atol=1e-9)
    np.testing.assert_allclose(vx2, vx, atol=1e-9)
    np.testing.assert_allclose(vy2, vy, atol=1e-9)
    np.testing.assert_allclose(vz2, vz, atol=1e-9)


def test_cart_cyl_known_point():
    # point on the +x axis at radius 3, moving purely in +y (tangential)
    r, theta, zed, vr, vtheta, vzed = cart_to_cyl(
        np.array([3.0]), np.array([0.0]), np.array([1.0]),
        np.array([0.0]), np.array([2.0]), np.array([0.0]),
    )
    assert np.isclose(r[0], 3.0)
    assert np.isclose(theta[0], 0.0, atol=1e-9)
    assert np.isclose(zed[0], 1.0)
    assert np.isclose(vr[0], 0.0, atol=1e-9)
    assert np.isclose(vtheta[0], 2.0)


def test_cart_to_sky_preserves_distance_scale():
    # A star along the sun-to-galactic-centre-ish direction: rather than
    # asserting a specific ra/dec (which depends on galpy's frame
    # conventions), check the physically-required invariant that the
    # reported distance to the star equals the norm of its position vector
    # relative to the sun's location, computed independently.
    ro, vo = 8.275, 220.0
    solarmotion = [-11.1, 12.24, 7.25]

    x = np.array([1000.0])  # pc, but cart_to_sky/galpy Orbit uses ro in kpc
    y = np.array([0.0])
    z = np.array([0.0])
    vx = np.array([0.0])
    vy = np.array([0.0])
    vz = np.array([0.0])

    # convert to kpc for galpy
    ra, dec, dist, pmra, pmdec, vlos = cart_to_sky(
        x / 1000.0, y / 1000.0, z / 1000.0, vx, vy, vz, ro=ro, vo=vo, solarmotion=solarmotion
    )

    # independent distance-to-sun calculation: sun sits at (ro, 0, zo=0)
    # in the galactocentric frame used here (galpy convention), so distance
    # from a point at (x,y,z) galactocentric (kpc) to the sun is:
    sun_x, sun_y, sun_z = ro, 0.0, 0.0
    expected_dist = np.sqrt(
        (x[0] / 1000.0 - sun_x) ** 2 + (y[0] / 1000.0 - sun_y) ** 2 + (z[0] / 1000.0 - sun_z) ** 2
    )

    assert np.isclose(dist[0], expected_dist, rtol=1e-3)
