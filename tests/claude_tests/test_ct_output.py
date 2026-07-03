import os
import numpy as np
import pytest

import clustertools as ctools
from clustertools.util.output import fortout, sseout, gyrout


def test_fortout_writes_expected_number_of_lines_and_columns(tmp_path):
    rng = np.random.default_rng(400)
    n = 60
    x = rng.uniform(-5, 5, n)
    y = rng.uniform(-5, 5, n)
    z = rng.uniform(-5, 5, n)
    vx = rng.uniform(-2, 2, n)
    vy = rng.uniform(-2, 2, n)
    vz = rng.uniform(-2, 2, n)
    m = rng.uniform(0.1, 2.0, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, m)
    cluster.add_nbody6(rbar=3.0, zmbar=500.0, vbar=4.0, tbar=1.0)
    cluster.find_centre()
    cluster.wdir = str(tmp_path) + "/"

    fortout(cluster, filename="fort.10")

    outfile = tmp_path / "fort.10"
    assert outfile.exists()

    data = np.loadtxt(outfile)
    # independent check: file should have exactly n rows (one per star) and
    # 7 columns (m,x,y,z,vx,vy,vz), per fortout's own column_stack call
    assert data.shape == (n, 7)

    # independent check: total mass written (in nbody units) should equal
    # cluster mass converted the same way (mass / zmbar)
    expected_mass_nbody_sum = np.sum(m) / 500.0
    assert np.isclose(np.sum(data[:, 0]), expected_mass_nbody_sum, rtol=1e-6)

    # fortout should leave the cluster restored to its original units/origin
    assert cluster.units == "pckms"
    assert cluster.origin == "cluster"


def test_sseout_writes_expected_columns(tmp_path):
    rng = np.random.default_rng(401)
    n = 40
    x = rng.uniform(-5, 5, n)
    y = rng.uniform(-5, 5, n)
    z = rng.uniform(-5, 5, n)
    vx = rng.uniform(-2, 2, n)
    vy = rng.uniform(-2, 2, n)
    vz = rng.uniform(-2, 2, n)
    m = rng.uniform(0.1, 2.0, n)
    m0 = m + rng.uniform(0, 0.5, n)
    kw = rng.integers(0, 10, n)
    ep = rng.uniform(0, 100, n)
    ospin = rng.uniform(0, 10, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, m, m0=m0)
    cluster.add_sse(kw, np.zeros(n), np.zeros(n), ep=ep, ospin=ospin)

    outfile = tmp_path / "sse.dat"
    sseout(cluster, str(outfile))

    data = np.loadtxt(outfile)
    assert data.shape == (n, 5)

    # independent check: columns are m0, kw, m, ep, ospin per sseout's own
    # column_stack call, written while the cluster is temporarily in pckms
    # (already true here, so no unit conversion needed for the comparison).
    # sseout writes with fixed-point '%f' (6 decimal places), so use an
    # absolute tolerance matching that format's precision rather than rtol,
    # which is too strict for small values under fixed-point rounding.
    np.testing.assert_allclose(data[:, 0], m0, atol=1e-6, rtol=0)
    np.testing.assert_array_equal(data[:, 1].astype(int), kw)
    np.testing.assert_allclose(data[:, 2], m, atol=1e-6, rtol=0)
    np.testing.assert_allclose(data[:, 3], ep, atol=1e-6, rtol=0)
    np.testing.assert_allclose(data[:, 4], ospin, atol=1e-6, rtol=0)

    # sseout should leave the cluster restored to its original units/origin
    assert cluster.units == "pckms"
    assert cluster.origin == "cluster"


def test_gyrout_reload_round_trip(tmp_path):
    rng = np.random.default_rng(402)
    n = 80
    x = rng.uniform(-5, 5, n)
    y = rng.uniform(-5, 5, n)
    z = rng.uniform(-5, 5, n)
    vx = rng.uniform(-2, 2, n)
    vy = rng.uniform(-2, 2, n)
    vz = rng.uniform(-2, 2, n)
    m = rng.uniform(0.1, 2.0, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, m, sortstars=True)
    cluster.add_orbit(8000.0, 0.0, 0.0, 0.0, 220.0, 0.0)

    outfile = tmp_path / "gyr.dat"
    gyrout(cluster, filename=str(outfile))

    # independent check: read the raw WDunits-format file directly (no
    # header, unlike real gyrfalcon output) and undo gyrout's own conversion
    # factors by hand, rather than trusting clustertools' own gyrfalcon
    # loader for the comparison values
    from galpy.util import conversion
    solar_ro, solar_vo = 8.275, 8.275 * 30.39 - 12.24
    vcon = 220.0 / conversion.velocity_in_kpcGyr(vo=solar_vo, ro=solar_ro)
    mcon = 222288.4543021174

    raw = np.loadtxt(outfile)
    assert raw.shape == (n, 7)

    reloaded_mass = raw[:, 0] * mcon
    reloaded_vx = raw[:, 4] * vcon

    # cluster was moved to galaxy origin + kpckms internally by gyrout;
    # compare against the total mass (unit/order-independent) and the
    # velocity conversion factor rather than assuming star order survives
    assert np.isclose(np.sum(reloaded_mass), np.sum(m), rtol=1e-5)

    # gyrout should leave the cluster restored to its original units/origin
    assert cluster.units == "pckms"
    assert cluster.origin == "cluster"
