import numpy as np
import matplotlib
import pytest

matplotlib.use("Agg")

import clustertools as ctools

RNG_SEED = 42


@pytest.fixture
def rng():
    return np.random.default_rng(RNG_SEED)


def _uniform_cluster(n=500, rmax=10.0, vscale=1.0, rng=None, mass=1.0):
    """Random points inside a sphere of radius rmax, equal mass, isotropic velocities."""
    if rng is None:
        rng = np.random.default_rng(RNG_SEED)

    # sample uniformly within a sphere (not uniform density in r, just a generic
    # scattered population -- fine for structural/consistency tests)
    x = rng.uniform(-rmax, rmax, n)
    y = rng.uniform(-rmax, rmax, n)
    z = rng.uniform(-rmax, rmax, n)
    r = np.sqrt(x**2 + y**2 + z**2)
    keep = r <= rmax
    x, y, z = x[keep], y[keep], z[keep]

    vx = rng.normal(0, vscale, len(x))
    vy = rng.normal(0, vscale, len(x))
    vz = rng.normal(0, vscale, len(x))
    m = np.ones(len(x)) * mass

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, m, sortstars=True)
    return cluster


@pytest.fixture
def equal_mass_cluster(rng):
    """Small cluster, all stars equal mass -- makes half-mass/Lagrange radii
    exactly reproducible from np.percentile-style definitions."""
    return _uniform_cluster(n=300, rmax=10.0, vscale=2.0, rng=rng)


@pytest.fixture
def equal_mass_cluster_factory(rng):
    """Factory so tests can request a specific N."""

    def _make(n=300, rmax=10.0, vscale=2.0, mass=1.0):
        return _uniform_cluster(n=n, rmax=rmax, vscale=vscale, rng=rng, mass=mass)

    return _make


@pytest.fixture
def two_body_cluster():
    """Two stars at rest relative to each other, separated along x, for
    analytic energy checks."""
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    m1, m2 = 3.0, 5.0
    sep = 4.0
    cluster.add_stars(
        np.array([0.0, sep]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        np.array([0.0, 0.0]),
        np.array([m1, m2]),
    )
    return cluster


@pytest.fixture
def plummer_cluster():
    """Synthetic Plummer-model realization via limepy, skipped if limepy is
    not installed or is incompatible with the installed scipy (the
    installed limepy 1.3.0 passes a float where scipy's ode.integrate now
    requires an int for its dopri5 backend on scipy>=1.17 -- an upstream
    limepy/scipy version mismatch, not a clustertools bug)."""
    pytest.importorskip("limepy")
    try:
        cluster = ctools.load_cluster(
            ctype="limepy", units="pckms", origin="cluster", g=1, phi0=5.0, M=1000.0, rh=3.0, N=2000
        )
    except TypeError as e:
        pytest.skip(f"limepy is incompatible with installed scipy: {e}")
    return cluster
