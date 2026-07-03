import numpy as np
import pytest

import clustertools as ctools

# Independently derived (not copied from tests/test_units.py) conversion
# constants, worked out from first principles:
KPC_TO_PC = 1000.0
KMS_TO_PCMYR = 1.022712165045695  # 1 km/s in pc/Myr (1 pc = 3.0857e13 km, 1 Myr = 3.1557e13 s)


@pytest.mark.parametrize(
    "from_units,to_units",
    [
        ("pckms", "kpckms"),
        ("kpckms", "pckms"),
        ("pckms", "pcmyr"),
        ("pcmyr", "pckms"),
        ("pckms", "nbody"),
        ("nbody", "pckms"),
        ("pckms", "kpcgyr"),
        ("kpcgyr", "pckms"),
        ("kpckms", "kpcgyr"),
    ],
)
def test_unit_round_trip_recovers_positions_and_velocities(from_units, to_units):
    rng = np.random.default_rng(20)
    n = 50
    x = rng.uniform(-10, 10, n)
    y = rng.uniform(-10, 10, n)
    z = rng.uniform(-10, 10, n)
    vx = rng.uniform(-5, 5, n)
    vy = rng.uniform(-5, 5, n)
    vz = rng.uniform(-5, 5, n)

    cluster = ctools.StarCluster(units=from_units, origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz)
    cluster.rbar = 3.0
    cluster.zmbar = 100.0
    cluster.vbar = 5.0
    cluster.tbar = 2.0

    cluster.to_units(to_units)
    cluster.to_units(from_units)

    np.testing.assert_allclose(cluster.x, x, rtol=1e-8)
    np.testing.assert_allclose(cluster.y, y, rtol=1e-8)
    np.testing.assert_allclose(cluster.z, z, rtol=1e-8)
    np.testing.assert_allclose(cluster.vx, vx, rtol=1e-8)
    np.testing.assert_allclose(cluster.vy, vy, rtol=1e-8)
    np.testing.assert_allclose(cluster.vz, vz, rtol=1e-8)


def test_pckms_to_kpckms_factor_of_a_thousand():
    x0 = 500.0
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(np.array([x0]), np.array([0.0]), np.array([0.0]),
                       np.array([1.0]), np.array([0.0]), np.array([0.0]))

    cluster.to_kpckms()

    assert np.isclose(cluster.x[0], x0 / KPC_TO_PC)
    # velocities are km/s in both pckms and kpckms -- unchanged
    assert np.isclose(cluster.vx[0], 1.0)


def test_pckms_to_pcmyr_velocity_conversion():
    v0 = 100.0  # km/s
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(np.array([1.0]), np.array([0.0]), np.array([0.0]),
                       np.array([v0]), np.array([0.0]), np.array([0.0]))

    cluster.to_pcmyr()

    # independently derived: 1 km/s = 1.02271 pc/Myr
    expected_v = v0 * KMS_TO_PCMYR
    assert np.isclose(cluster.vx[0], expected_v, rtol=1e-5)
    # positions unaffected (both in pc)
    assert np.isclose(cluster.x[0], 1.0)


def test_pckms_to_nbody_uses_rbar_zmbar_vbar():
    rbar, zmbar, vbar = 4.0, 200.0, 3.0
    x0, v0, m0 = 8.0, 6.0, 0.5

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(np.array([x0]), np.array([0.0]), np.array([0.0]),
                       np.array([v0]), np.array([0.0]), np.array([0.0]),
                       np.array([m0]))
    cluster.rbar = rbar
    cluster.zmbar = zmbar
    cluster.vbar = vbar

    cluster.to_nbody()

    assert np.isclose(cluster.x[0], x0 / rbar, rtol=1e-6)
    assert np.isclose(cluster.vx[0], v0 / vbar, rtol=1e-6)
    assert np.isclose(cluster.m[0], m0 / zmbar, rtol=1e-6)


def test_kpckms_to_kpcgyr_velocity_conversion():
    v0 = 50.0  # km/s
    cluster = ctools.StarCluster(units="kpckms", origin="cluster")
    cluster.add_stars(np.array([1.0]), np.array([0.0]), np.array([0.0]),
                       np.array([v0]), np.array([0.0]), np.array([0.0]))

    cluster.to_kpcgyr()

    # 1 km/s = 1.02271 kpc/Gyr numerically (same constant, different length/time scale)
    expected_v = v0 * KMS_TO_PCMYR
    assert np.isclose(cluster.vx[0], expected_v, rtol=1e-5)


def test_add_orbit_ounits_conversion_is_consistent_with_to_units():
    # If we set an orbit in 'kpckms' on a cluster whose stellar units are
    # 'pckms', add_orbit should convert xgc/ygc/zgc to pc the same way
    # to_units converts stellar positions pc<->kpc.
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(np.array([0.0]), np.array([0.0]), np.array([0.0]),
                       np.array([0.0]), np.array([0.0]), np.array([0.0]))

    xgc_kpc, ygc_kpc, zgc_kpc = 8.0, 0.5, 0.1
    vxgc_kkms, vygc_kkms, vzgc_kkms = 10.0, 20.0, 5.0

    cluster.add_orbit(xgc_kpc, ygc_kpc, zgc_kpc, vxgc_kkms, vygc_kkms, vzgc_kkms, ounits="kpckms")

    assert np.isclose(cluster.xgc, xgc_kpc * KPC_TO_PC)
    assert np.isclose(cluster.ygc, ygc_kpc * KPC_TO_PC)
    assert np.isclose(cluster.zgc, zgc_kpc * KPC_TO_PC)
    # velocities: km/s in both pckms and kpckms -- unchanged
    assert np.isclose(cluster.vxgc, vxgc_kkms)
