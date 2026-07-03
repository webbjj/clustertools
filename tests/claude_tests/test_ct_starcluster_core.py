import numpy as np
import pytest
from scipy.spatial import cKDTree

import clustertools as ctools
from clustertools.cluster.operations import reset_nbody_scale, to_WDunits
from clustertools.cluster.cluster import overlap_cluster


def test_constructor_defaults():
    cluster = ctools.StarCluster()
    assert cluster.tphys == 0.0
    assert cluster.units is None
    assert cluster.origin is None
    assert cluster.ctype == "snapshot"
    assert cluster.projected is False
    assert cluster.ntot == 0


def test_constructor_custom_args():
    cluster = ctools.StarCluster(5.0, "pckms", "cluster", "nbody6", True)
    assert cluster.tphys == 5.0
    assert cluster.units == "pckms"
    assert cluster.origin == "cluster"
    assert cluster.ctype == "nbody6"
    assert cluster.projected is True


def test_add_stars_default_mass_and_id():
    n = 50
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    x = np.arange(n, dtype=float)
    y = np.zeros(n)
    z = np.zeros(n)
    vx = np.zeros(n)
    vy = np.zeros(n)
    vz = np.zeros(n)
    cluster.add_stars(x, y, z, vx, vy, vz)

    np.testing.assert_array_equal(cluster.x, x)
    np.testing.assert_array_equal(cluster.m, np.ones(n))
    np.testing.assert_array_equal(cluster.id, np.arange(n))
    assert cluster.ntot == n


def test_add_stars_binary_split():
    n = 40
    nb = 6
    rng = np.random.default_rng(30)
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    x = rng.uniform(-5, 5, n)
    y = rng.uniform(-5, 5, n)
    z = rng.uniform(-5, 5, n)
    vx = rng.uniform(-1, 1, n)
    vy = rng.uniform(-1, 1, n)
    vz = rng.uniform(-1, 1, n)

    cluster.add_stars(x, y, z, vx, vy, vz, nb=nb)

    assert len(cluster.x) == n - nb
    assert len(cluster.xb1) == nb
    assert len(cluster.xb2) == nb
    assert cluster.nb == nb

    # centre-of-mass positions of the binaries should independently match
    # the mass-weighted average of the two components (equal mass here,
    # since add_stars defaults m=1 for all input stars)
    arg1 = np.arange(0, 2 * nb - 1, 2)
    arg2 = arg1 + 1
    expected_xcom = 0.5 * (x[arg1] + x[arg2])
    # the first nb entries of cluster.x are the binary centre-of-mass stars
    np.testing.assert_allclose(np.sort(cluster.x[:nb]), np.sort(expected_xcom))


def test_add_binary_stars_directly():
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    xb1, yb1, zb1 = np.array([0.0]), np.array([0.0]), np.array([0.0])
    vxb1, vyb1, vzb1 = np.array([0.0]), np.array([0.0]), np.array([0.0])
    xb2, yb2, zb2 = np.array([2.0]), np.array([0.0]), np.array([0.0])
    vxb2, vyb2, vzb2 = np.array([0.0]), np.array([0.0]), np.array([0.0])
    mb1, mb2 = np.array([1.0]), np.array([3.0])

    xcom, ycom, zcom, vxcom, vycom, vzcom, mcom = cluster.add_binary_stars(
        xb1, yb1, zb1, vxb1, vyb1, vzb1, xb2, yb2, zb2, vxb2, vyb2, vzb2, mb1, mb2, return_com=True
    )

    # independent mass-weighted centre of mass: (1*0 + 3*2)/4 = 1.5
    assert np.isclose(xcom[0], 1.5)
    assert np.isclose(mcom[0], 4.0)
    assert cluster.nb == 1


def test_add_orbit_sets_galactocentric_position():
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(np.array([0.0]), np.array([0.0]), np.array([0.0]),
                       np.array([0.0]), np.array([0.0]), np.array([0.0]))
    cluster.add_orbit(100.0, 200.0, 50.0, 10.0, -20.0, 5.0)

    assert cluster.xgc == 100.0
    assert cluster.ygc == 200.0
    assert cluster.zgc == 50.0
    # independently recompute galactocentric radius
    expected_rgc = np.sqrt(100.0**2 + 200.0**2 + 50.0**2)
    assert np.isclose(cluster.rgc, expected_rgc)


def test_add_nbody6_sets_scaling_params():
    cluster = ctools.StarCluster(units="nbody", origin="cluster")
    cluster.add_nbody6(nc=5, rc=1.5, rbar=3.0, rtide=20.0, zmbar=500.0, vbar=4.0, tbar=2.0)

    assert cluster.nc == 5
    assert cluster.rc == 1.5
    assert cluster.rbar == 3.0
    assert cluster.zmbar == 500.0

    # tbar_days is derived from rbar/zmbar via Kepler's third law; recompute
    # independently and check consistency
    au_to_cm = 1.49597870700e13
    pc_to_cm = 1296000.0 / (2.0 * np.pi) * au_to_cm
    nbody_to_years = (cluster.rbar * 1296000.0 / (2.0 * np.pi)) ** 1.5 / np.sqrt(cluster.zmbar)
    expected_tbar_days = 365.25 * nbody_to_years
    assert np.isclose(cluster.tbar_days, expected_tbar_days)


def test_add_sse_and_bse():
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(np.array([0.0, 1.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]),
                       np.array([0.0, 0.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]))
    cluster.add_sse(np.array([0, 1]), np.array([1.0, 2.0]), np.array([0.1, 0.2]))
    np.testing.assert_array_equal(cluster.kw, np.array([0, 1]))

    cluster.add_bse(
        id1=np.array([0]), id2=np.array([1]), kw1=np.array([0]), kw2=np.array([1]),
        kcm=np.array([0]), ecc=np.array([0.1]), pb=np.array([10.0]), semi=np.array([2.0]),
        m1=np.array([1.0]), m2=np.array([2.0]), logl1=np.array([0.0]), logl2=np.array([0.0]),
        logr1=np.array([0.0]), logr2=np.array([0.0]),
    )
    # binding energy eb = 0.5*m1*m2/semi is computed internally -- verify
    # independently
    expected_eb = 0.5 * 1.0 * 2.0 / 2.0
    assert np.isclose(cluster.eb[0], expected_eb)


def test_add_energies_computes_qvir():
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(np.array([0.0, 1.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]),
                       np.array([0.0, 0.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]))
    kin = np.array([2.0, 3.0])
    pot = np.array([-4.0, -4.0])
    cluster.add_energies(kin, pot)

    expected_ektot = np.sum(kin)
    expected_ptot = np.sum(pot) / 2.0
    assert np.isclose(cluster.ektot, expected_ektot)
    assert np.isclose(cluster.ptot, expected_ptot)
    assert np.isclose(cluster.qvir, expected_ektot / expected_ptot)


def test_analyze_r_and_v_are_correct_norms(equal_mass_cluster):
    cluster = equal_mass_cluster
    expected_r = np.sqrt(cluster.x**2 + cluster.y**2 + cluster.z**2)
    expected_v = np.sqrt(cluster.vx**2 + cluster.vy**2 + cluster.vz**2)
    np.testing.assert_allclose(cluster.r, expected_r)
    np.testing.assert_allclose(cluster.v, expected_v)


def test_half_mass_radius_matches_definition_for_equal_mass(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=999, mass=1.0)

    # independent computation of the half-mass radius from its definition:
    # the radius that encloses half the total mass, using cumulative sum
    # over stars sorted by radius (not calling any clustertools internals)
    order = np.argsort(cluster.r)
    m_sorted = cluster.m[order]
    r_sorted = cluster.r[order]
    msum = np.cumsum(m_sorted)
    idx = np.searchsorted(msum, 0.5 * cluster.mtot)
    expected_rm = r_sorted[idx]

    assert np.isclose(cluster.rm, expected_rm)


def test_sortstars_produces_nondecreasing_order(equal_mass_cluster):
    cluster = equal_mass_cluster
    cluster.sortstars()
    r_ordered = cluster.r[cluster.rorder]
    assert np.all(np.diff(r_ordered) >= 0)


def test_subset_mass_bounds(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=500)
    # give stars varying mass so mmin/mmax is meaningful
    rng = np.random.default_rng(31)
    cluster.m = rng.uniform(0.1, 10.0, cluster.ntot)

    mmin, mmax = 1.0, 5.0
    indx = cluster.subset(mmin=mmin, mmax=mmax)

    assert np.all(cluster.m[indx] >= mmin)
    assert np.all(cluster.m[indx] <= mmax)
    # independent check: every star outside the mask should violate the bound
    outside = ~indx
    assert np.all((cluster.m[outside] < mmin) | (cluster.m[outside] > mmax))


def test_subset_velocity_bounds(equal_mass_cluster):
    cluster = equal_mass_cluster
    vmax = 1.0
    indx = cluster.subset(vmax=vmax)
    assert np.all(cluster.v[indx] <= vmax)


def test_add_action_and_add_actions_store_values_directly():
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(np.array([0.0]), np.array([0.0]), np.array([0.0]),
                       np.array([0.0]), np.array([0.0]), np.array([0.0]))

    cluster.add_action(1.0, 2.0, 3.0, OR=4.0, Ophi=5.0, Oz=6.0, TR=7.0, Tphi=8.0, Tz=9.0)
    assert (cluster.JR, cluster.Jphi, cluster.Jz) == (1.0, 2.0, 3.0)
    assert (cluster.OR, cluster.Ophi, cluster.Oz) == (4.0, 5.0, 6.0)
    assert (cluster.TR, cluster.Tphi, cluster.Tz) == (7.0, 8.0, 9.0)

    JR = np.array([1.0, 2.0])
    Jphi = np.array([3.0, 4.0])
    Jz = np.array([5.0, 6.0])
    cluster.add_actions(JR, Jphi, Jz)
    np.testing.assert_array_equal(cluster.JRs, JR)
    np.testing.assert_array_equal(cluster.Jphis, Jphi)
    np.testing.assert_array_equal(cluster.Jzs, Jz)


def test_reset_nbody_scale_matches_independent_formula(equal_mass_cluster_factory):
    cluster = equal_mass_cluster_factory(n=500, rmax=8.0, mass=1.0)
    rng = np.random.default_rng(1000)
    cluster.m = rng.uniform(0.1, 5.0, cluster.ntot)
    cluster.find_centre()
    cluster.to_centre()

    zmbar, rbar, vbar, tbar = reset_nbody_scale(cluster, rvirial=False)

    # independent recomputation from the definitions in operations.py:
    # zmbar = total mass, rbar = 4/3 * half-mass radius (non-virial mode),
    # vbar/tbar from the standard N-body scaling relations
    expected_zmbar = np.sum(cluster.m)
    expected_rbar = 4.0 * cluster.rm / 3.0
    expected_vbar = 0.06557 * np.sqrt(expected_zmbar / expected_rbar)
    expected_tbar = expected_rbar / (1.023 * expected_vbar)

    assert np.isclose(zmbar, expected_zmbar, rtol=1e-8)
    assert np.isclose(rbar, expected_rbar, rtol=1e-6)
    assert np.isclose(vbar, expected_vbar, rtol=1e-8)
    assert np.isclose(tbar, expected_tbar, rtol=1e-8)


def test_overlap_cluster_matches_independent_kdtree_cross_query():
    rng = np.random.default_rng(1001)
    n1, n2 = 300, 300

    x1 = rng.uniform(-10, 10, n1)
    y1 = rng.uniform(-10, 10, n1)
    z1 = rng.uniform(-10, 10, n1)
    cluster1 = ctools.StarCluster(units="pckms", origin="cluster")
    cluster1.add_stars(x1, y1, z1, np.zeros(n1), np.zeros(n1), np.zeros(n1))

    # cluster2 overlaps cluster1's region only partially (shifted centre)
    x2 = 5.0 + rng.uniform(-10, 10, n2)
    y2 = rng.uniform(-10, 10, n2)
    z2 = rng.uniform(-10, 10, n2)
    cluster2 = ctools.StarCluster(units="pckms", origin="cluster")
    cluster2.add_stars(x2, y2, z2, np.zeros(n2), np.zeros(n2), np.zeros(n2))

    tol = 0.5
    indx = overlap_cluster(cluster1, cluster2, tol=tol, return_cluster=False)

    # independent cross-set nearest-neighbour distance via scipy's cKDTree
    tree2 = cKDTree(np.column_stack([x2, y2, z2]))
    dmin, _ = tree2.query(np.column_stack([x1, y1, z1]))
    expected_indx = dmin < tol

    np.testing.assert_array_equal(indx, expected_indx)
    assert np.sum(indx) > 0  # sanity: the two regions do overlap somewhat


def test_to_WDunits_matches_independent_mass_velocity_conversion():
    from galpy.util import conversion

    rng = np.random.default_rng(1002)
    n = 100
    x = rng.uniform(-10, 10, n)
    y = rng.uniform(-10, 10, n)
    z = rng.uniform(-10, 10, n)
    vx = rng.normal(0, 50, n)
    vy = rng.normal(0, 50, n)
    vz = rng.normal(0, 50, n)
    m = rng.uniform(0.1, 2.0, n)

    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(x, y, z, vx, vy, vz, m)

    m_before_kpckms = m.copy()
    vx_before = vx.copy()

    cluster.to_WDunits()

    ro, vo = cluster._ro, cluster._vo
    vcon = vo / conversion.velocity_in_kpcGyr(vo=vo, ro=ro)
    mcon = 222288.4543021174

    assert cluster.units == "WDunits"
    np.testing.assert_allclose(cluster.m, m_before_kpckms / mcon, rtol=1e-8)
    np.testing.assert_allclose(cluster.vx, vx_before / vcon, rtol=1e-8)


def test_to_audays_and_to_sudays_round_trip():
    cluster = ctools.StarCluster(units="pckms", origin="cluster")
    cluster.add_stars(np.array([0.0, 1.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]),
                       np.array([0.0, 0.0]), np.array([0.0, 0.0]), np.array([0.0, 0.0]))
    cluster.add_nbody6(rbar=3.0, zmbar=500.0, vbar=4.0, tbar=1.0)

    semi0 = np.array([2.0])
    pb0 = np.array([10.0])
    m1_0 = np.array([1.0])
    m2_0 = np.array([2.0])

    cluster.add_bse(
        id1=np.array([0]), id2=np.array([1]), kw1=np.array([0]), kw2=np.array([1]),
        kcm=np.array([0]), ecc=np.array([0.1]), pb=pb0.copy(), semi=semi0.copy(),
        m1=m1_0.copy(), m2=m2_0.copy(), logl1=np.array([0.0]), logl2=np.array([0.0]),
        logr1=np.array([0.0]), logr2=np.array([0.0]),
    )
    cluster.bunits = "nbody"

    rbar_au = cluster.rbar_au
    tbar_days = cluster.tbar_days
    zmbar = cluster.zmbar

    # nbody -> audays: semi *= rbar_au, pb *= tbar_days, m1/m2 *= zmbar
    cluster.to_audays()
    assert cluster.bunits == "audays"
    np.testing.assert_allclose(cluster.semi, semi0 * rbar_au, rtol=1e-8)
    np.testing.assert_allclose(cluster.pb, pb0 * tbar_days, rtol=1e-8)
    np.testing.assert_allclose(cluster.m1, m1_0 * zmbar, rtol=1e-8)

    # to_sudays(), called from bunits=='audays', cascades through both of
    # its own branches in one call (audays->nbody, then immediately
    # nbody->sudays), ending at 'sudays' rather than stopping at the
    # intermediate 'nbody' state -- confirmed by reading
    # cluster/operations.py::to_sudays. Since pb/m1/m2 use the same zmbar
    # and tbar_days conversion regardless of sudays vs audays, only semi
    # (which uses rbar_su vs rbar_au) actually changes value here.
    rbar_su = cluster.rbar_su
    cluster.to_sudays()
    assert cluster.bunits == "sudays"
    np.testing.assert_allclose(cluster.semi, semi0 * rbar_su, rtol=1e-6)
    np.testing.assert_allclose(cluster.pb, pb0 * tbar_days, rtol=1e-6)
    np.testing.assert_allclose(cluster.m1, m1_0 * zmbar, rtol=1e-6)
