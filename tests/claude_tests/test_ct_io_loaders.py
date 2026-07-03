import os
import numpy as np
import pytest

import clustertools as ctools

NOTEBOOKS_DIR = os.path.join(
    os.path.dirname(__file__), os.pardir, os.pardir, "docs", "source", "notebooks"
)


def test_gyrfalcon_loader_star_count_and_mass_match_independent_parse():
    path = os.path.join(NOTEBOOKS_DIR, "cluster.nemo.dat")
    if not os.path.isfile(path):
        pytest.skip(f"data file not found: {path}")

    # independently confirm the file's structure before trusting the loader:
    # the first snapshot's particle count is stated in its own header line
    with open(path) as f:
        header_lines = [next(f) for _ in range(13)]
    ntot_line = [l for l in header_lines if "Ntot" in l][0]
    expected_ntot = int(ntot_line.split("Ntot:")[1].split(",")[0].strip())

    raw = np.loadtxt(path, skiprows=13, max_rows=expected_ntot)
    expected_mass_sum = np.sum(raw[:, 0])

    cluster = ctools.load_cluster(
        ctype="gyrfalcon", filename=path, units="WDunits", origin="centre"
    )

    assert cluster.ntot == expected_ntot
    assert np.isclose(np.sum(cluster.m), expected_mass_sum, rtol=1e-6)
    np.testing.assert_allclose(np.sort(cluster.x), np.sort(raw[:, 1]), rtol=1e-6)


def test_nbody6pp_loader_runs_and_advances(tmp_path):
    wdir = os.path.join(NOTEBOOKS_DIR, "")
    if not os.path.isfile(os.path.join(wdir, "conf.3_0")):
        pytest.skip(f"data file not found in {wdir}")

    cluster = ctools.load_cluster("nbody6pp", wdir=wdir)

    assert cluster.ntot > 0
    assert cluster.units == "nbody"
    assert np.all(np.isfinite(cluster.x))
    assert np.all(np.isfinite(cluster.m))
    assert np.isclose(cluster.tphys, 0.0, atol=1e-6)

    if os.path.isfile(os.path.join(wdir, "conf.3_20")):
        advanced = ctools.advance_cluster(cluster, dtout=20)
        assert advanced.ntot > 0
        assert advanced.tphys > cluster.tphys


def test_astropy_table_loader_round_trips_known_columns():
    astropy = pytest.importorskip("astropy")
    from astropy.table import QTable

    rng = np.random.default_rng(950)
    n = 200
    m = rng.uniform(0.1, 2.0, n)
    x = rng.uniform(-10, 10, n)
    y = rng.uniform(-10, 10, n)
    z = rng.uniform(-10, 10, n)
    vx = rng.normal(0, 1, n)
    vy = rng.normal(0, 1, n)
    vz = rng.normal(0, 1, n)

    table = QTable()
    table["m"] = m
    table["x"] = x
    table["y"] = y
    table["z"] = z
    table["vx"] = vx
    table["vy"] = vy
    table["vz"] = vz

    cluster = ctools.load_cluster(
        ctype="astropy_table", particles=table, units="pckms", origin="cluster"
    )

    assert cluster.ntot == n
    np.testing.assert_allclose(cluster.m, m)
    np.testing.assert_allclose(cluster.x, x)
    np.testing.assert_allclose(cluster.vz, vz)


def test_astropy_table_loader_with_explicit_column_mapper():
    astropy = pytest.importorskip("astropy")
    from astropy.table import QTable

    rng = np.random.default_rng(951)
    n = 100
    mass = rng.uniform(0.1, 2.0, n)
    px = rng.uniform(-10, 10, n)
    py = rng.uniform(-10, 10, n)
    pz = rng.uniform(-10, 10, n)
    pvx = rng.normal(0, 1, n)
    pvy = rng.normal(0, 1, n)
    pvz = rng.normal(0, 1, n)

    table = QTable()
    table["mass_col"] = mass
    table["px"] = px
    table["py"] = py
    table["pz"] = pz
    table["pvx"] = pvx
    table["pvy"] = pvy
    table["pvz"] = pvz

    cm = {"m": "mass_col", "x": "px", "y": "py", "z": "pz", "vx": "pvx", "vy": "pvy", "vz": "pvz"}

    cluster = ctools.load_cluster(
        ctype="astropy_table", particles=table, units="pckms", origin="cluster", column_mapper=cm
    )

    assert cluster.ntot == n
    np.testing.assert_allclose(cluster.m, mass)
    np.testing.assert_allclose(cluster.x, px)
