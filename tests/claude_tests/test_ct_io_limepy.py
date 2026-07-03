import numpy as np
import pytest

pytest.importorskip("limepy")

import clustertools as ctools
from clustertools.analysis.profiles import rho_prof


def test_limepy_cluster_half_mass_radius_matches_requested_rh(plummer_cluster):
    cluster = plummer_cluster
    cluster.find_centre()
    cluster.to_centre()
    cluster.analyze()

    # the model was generated with rh=3.0 pc -- clustertools' own computed
    # half-mass radius (from analyze(), independent of the limepy model
    # object itself) should recover that value for an N-body realization
    assert np.isclose(cluster.rm, 3.0, rtol=0.15)


def test_limepy_cluster_total_mass_matches_requested_M(plummer_cluster):
    cluster = plummer_cluster
    assert np.isclose(np.sum(cluster.m), 1000.0, rtol=0.05)


def test_limepy_cluster_density_profile_decreases_outward(plummer_cluster):
    cluster = plummer_cluster
    cluster.find_centre()
    cluster.to_centre()

    rprof, pprof, nprof = rho_prof(cluster, nrad=10)
    assert np.all(np.diff(pprof) <= 0)
