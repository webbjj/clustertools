import os
import numpy as np
import pytest

import clustertools as ctools
from clustertools.util.output import snapout

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), os.pardir)


def _count_lines(path):
    with open(path) as f:
        return sum(1 for _ in f)


@pytest.mark.parametrize(
    "filename",
    ["ngc6101_pckms_cluster.dat", "g1phi5rh3m10000.dat"],
)
def test_load_real_snapshot_star_count_matches_independent_line_count(filename):
    path = os.path.join(TEST_DATA_DIR, filename)
    expected_n = _count_lines(path)

    cluster = ctools.load_cluster(
        ctype="snapshot", filename=path, units="pckms", origin="cluster"
    )

    assert cluster.ntot == expected_n


def test_load_real_snapshot_mass_sum_matches_independent_column_sum():
    path = os.path.join(TEST_DATA_DIR, "ngc6101_pckms_cluster.dat")
    raw = np.loadtxt(path)
    expected_mass_sum = np.sum(raw[:, 0])

    cluster = ctools.load_cluster(
        ctype="snapshot", filename=path, units="pckms", origin="cluster"
    )

    assert np.isclose(np.sum(cluster.m), expected_mass_sum, rtol=1e-10)
    # positions/velocities should also match the raw file columns exactly
    # (col order is m,x,y,z,vx,vy,vz,id per clustertools/io/snapshot.py)
    np.testing.assert_allclose(cluster.x, raw[:, 1])
    np.testing.assert_allclose(cluster.vz, raw[:, 6])


def test_load_real_snapshot_has_no_nans_or_infs():
    path = os.path.join(TEST_DATA_DIR, "g1phi5rh3m10000.dat")
    cluster = ctools.load_cluster(
        ctype="snapshot", filename=path, units="pckms", origin="cluster"
    )

    for arr in (cluster.x, cluster.y, cluster.z, cluster.vx, cluster.vy, cluster.vz, cluster.m):
        assert np.all(np.isfinite(arr))


def test_snapout_reload_round_trip(tmp_path):
    rng = np.random.default_rng(200)
    n = 100
    x = rng.uniform(-5, 5, n)
    y = rng.uniform(-5, 5, n)
    z = rng.uniform(-5, 5, n)
    vx = rng.uniform(-2, 2, n)
    vy = rng.uniform(-2, 2, n)
    vz = rng.uniform(-2, 2, n)
    m = rng.uniform(0.1, 2.0, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, m)

    outfile = tmp_path / "snap_roundtrip.dat"
    snapout(cluster, str(outfile))

    reloaded = ctools.load_cluster(
        ctype="snapshot",
        filename=str(outfile),
        units="pckms",
        origin="cluster",
        col_names=["m", "x", "y", "z", "vx", "vy", "vz", "id"],
        col_nums=[0, 1, 2, 3, 4, 5, 6, 7],
    )

    assert reloaded.ntot == n
    # snapout writes stars in cluster.id order, not necessarily original
    # array order -- sort both by id before comparing
    orig_order = np.argsort(cluster.id)
    reload_order = np.argsort(reloaded.id)

    np.testing.assert_allclose(reloaded.x[reload_order], cluster.x[orig_order], rtol=1e-10)
    np.testing.assert_allclose(reloaded.y[reload_order], cluster.y[orig_order], rtol=1e-10)
    np.testing.assert_allclose(reloaded.z[reload_order], cluster.z[orig_order], rtol=1e-10)
    np.testing.assert_allclose(reloaded.vx[reload_order], cluster.vx[orig_order], rtol=1e-10)
    np.testing.assert_allclose(reloaded.m[reload_order], cluster.m[orig_order], rtol=1e-10)
