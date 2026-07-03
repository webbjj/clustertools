import numpy as np
import pytest

import clustertools as ctools
from clustertools.cluster.operations import add_rotation, virialize
try:
    from galpy.util import coords
except ImportError:
    import galpy.util.bovy_coords as coords


def test_add_rotation_called_directly_flips_negative_vtheta_stars():
    rng = np.random.default_rng(500)
    n = 2000
    x = rng.uniform(-5, 5, n)
    y = rng.uniform(-5, 5, n)
    z = rng.uniform(-5, 5, n)
    vx = rng.uniform(-2, 2, n)
    vy = rng.uniform(-2, 2, n)
    vz = rng.uniform(-2, 2, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz)

    # independent computation of vtheta before rotation
    r, theta, zed = coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
    vr0, vtheta0, vz0 = coords.rect_to_cyl_vec(cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z)
    n_negative_before = np.sum(vtheta0 < 0)

    # qrot=1.0 -- every star with vtheta<0 should be flipped to vtheta>=0
    x2, y2, z2, vx2, vy2, vz2 = add_rotation(cluster, qrot=1.0)

    r2, theta2, zed2 = coords.rect_to_cyl(x2, y2, z2)
    vr2, vtheta2, vz2c = coords.rect_to_cyl_vec(vx2, vy2, vz2, x2, y2, z2)

    assert n_negative_before > 0  # sanity: the test setup actually has some
    assert np.all(vtheta2 >= -1e-9)

    # add_rotation should not have changed the cluster's actual x/y/z/vx/vy/vz
    # attributes (it returns new arrays rather than mutating in place)
    np.testing.assert_array_equal(cluster.x, x)


def test_virialize_module_function_returns_consistent_scale_factor():
    rng = np.random.default_rng(501)
    n = 300
    x = rng.uniform(-8, 8, n)
    y = rng.uniform(-8, 8, n)
    z = rng.uniform(-8, 8, n)
    vx = rng.uniform(-3, 3, n)
    vy = rng.uniform(-3, 3, n)
    vz = rng.uniform(-3, 3, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, sortstars=True)
    cluster.find_centre()
    cluster.to_centre()
    # to_centre() legitimately shifts velocities by the cluster's
    # centre-of-mass velocity offset -- snapshot *after* that shift as the
    # baseline for confirming virialize() itself leaves velocities untouched.
    vx_before_virialize = cluster.vx.copy()

    qvir_target = 0.4
    qv = virialize(cluster, qvir=qvir_target)

    # virialize (module-level function) only computes and returns the scale
    # factor -- it does NOT rescale velocities itself (that happens in the
    # StarCluster.virialize() wrapper method). Confirm this directly: qv
    # should be sqrt(|qvir_target / cluster.qvir|), computed independently
    # from the same energies() the function used.
    kin, pot = cluster.energies(specific=True, full=True)
    ektot = np.sum(kin)
    ptot = np.sum(pot) / 2.0
    qvir_current = ektot / ptot

    expected_qv = np.sqrt(abs(qvir_target / qvir_current))
    assert np.isclose(qv, expected_qv, rtol=1e-6)

    # and confirm the module-level function indeed left velocities untouched
    np.testing.assert_array_equal(cluster.vx, vx_before_virialize)
