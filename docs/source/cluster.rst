.. _cluster:
Initialization
============

The StarCluster Class
---------------------

The ``StarCluster`` class is the foundation on which ``clustertools`` has been built. All functions have been designed to either act on or use elements of `StarCluster``. Its purpose is to be a python representation of a star cluster, based on a snapshot provided by a simulation or other software pacakge. At minimum, a ``StarCluster`` contains the positions and velocities of all stars in a star cluster. However, ``StarCluster`` has been designed to store all information related to the cluster's orbit, Single Star Evolution, Binary Star Evolution, and metadata related to the software used to generate the data and/or the data file it was read from. Once initiated, ``clustertools`` is able to perform a large number of operations, measurements, and calculations on the ``StarCluster`` class. 

To initialize a ``StarCluster``, one can simply start with:

>>> import clustertools as ctools

>>> cluster=ctools.StarCluster()

``StarCluster`` accepts additional optional arguments, each of which have defaults. They are the current time (``tphys=0``), the units (``units=None``) and origin (``origin=None``) of the coordinate system and the name of the code used to genereate the dataset ``ctype``. The units and origin variables only need to be specified if unit or coordinate trannsformations are going to be done (see :ref:`Units and Coordinate Systems <units_and_coordinate_systems>` for more information). The code type ``ctype`` defaults to ``'snapshot'``, but can alternatively be set to ``'astropy_table'``, ``'amuse'``, ``'galpy'``, ``'gyrfalcon'``, ``'limepy'``, ``'nbody6'``, ``'nbody6pp'``. ``ctype`` informs the ``StarCluster`` of the input files format. See :ref:`Loading and Advancing <_load_and_advance>` for more informatoin on how ``ctype`` is used. Other keywords accepted when initializing a ``StarCluster`` can be found in the complete documentation (see StarCluster).

Once a ``StarCluster`` is initialized, there are a large number of arrays and variables that correspond to individual stars, global properties of the cluster, and information related to the software used to generate the data. However, several functions have then been written to more easily populate ``StarCluster`` with information and carry out helpful calculations. They are:

* add_stars
* add_orbit
* add_nbody6
* add_sse
* add_bse
* add_energies
* add_action 
* add_actions 
* analyze 
* sortstars
* subset

.. automodapi:: clustertools.cluster.cluster
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading:
        :skip: sub_cluster, overlap_cluster
        :noindex:

The ``add_stars`` function is the intended way to actually add stars to ``cluster``. Assuming stellar positions (x,y,z) and velocities (vx,vy,vz) have already been read in via a snapshot, one can call:

>>> cluster.add_stars(x,y,z,vx,vy,vz)

Using ``add_stars`` as opposed to setting variables like ``cluster.x`` mantually ensures that ``cluster.analyze`` is automatically called and the three dimensional radius (``cluster.r``) and velocity (``cluster.v``) of each star is calculated. General cluster properties like total mass (``cluster.mtot``), mean radius (``cluster.rmean``) and maximum radius (``cluster.rmax``) are calculated. Projected values are also calculated assuming the x-y plane is the plane of the sky, and can be called by adding ``pro`` the variable name as in ``cluster.rmeanpro``.

It is also possible to include stellar masses ``m`` and ids ``id`` if they are known via: 

>>> cluster.add_stars(x,y,z,vx,vy,vz,m,id)

Otherwise, masses will be set to 1 and ids will simply be set to integer values between 1 and the number of stars in the cluster. ``add_stars`` can also accept the initial masses of each stars ``m0`` and the star's population number ``npop`` if one is studying multiple population within the cluster.

If your dataset contains binary stars, you can specificy the number of binary stars using ``nb`` via:

>>> cluster.add_stars(x,y,z,vx,vy,vz,nb=10)

Note that ``clustertools`` assumes that the arrays are structured such that the individual binary stars make up the first 2 x ``nb`` stars in the arrays. The position and centre of each binary's centre of mass is then determined and put in the main ``cluster`` arrays. The individual binary star positions and velocities are stored separately (e.g the x-position and velocity of the primary stars are saved as ``cluster.xb1`` and ``cluster.vxb1``). Alternatively you could add the binary stars manually if you have their individual position and velocities.

>>> cluster.add_binary_stars(xb1,yb1,zb1,vxb1,vyb1,vzb1,xb2,yb2,zb2,vxb2,vyb2,vzb2)

``add_binary_stars`` has the same features as ``add_stars`` in that masses (``mb1``, ``mb1``), ids (``id1``, ``id2``), initial masses (``m01``, ``m02``), and population numbers (``npop1``, ``npop2``) can be provided.


Finally, using ``add_stars``  and/or ``add_binary_stars`` results in ``cluster.ntot`` and ``cluster.nb`` being calculated. Otherwise they would have to be set manually.

One other features in ``add_stars`` that is by default set to ``True`` is ``sortstars`` . Having ``sortstars=True`` means stars are being sorted based on their distance from the origin (information stored in ``cluster.rorder``) and the half-mass radius ``rm`` and 10\% Lagrange radius ``r10`` will be calculated as well. Alternatively, one can call:

>>> cluster.analyze(sortstars=True)

at a later point in time if the cluster was populated manually. For computational efficiency, ``sortstars`` can be set to ``False`` when sorting is not required. The full list of parameters calculated by ``analyze`` is:

* ``r`` - clustercentric or galactocentric radius
* ``rpro`` - projeced clustercentric or galactocentric radius
* ``v`` - total velocity
* ``vpro`` - projected total velocity
* ``mtot`` - total mass
* ``mmean`` - mean mass
* ``rmean`` - mean radius
* ``rmeanpro`` - mean projected radius
* ``rmax`` - maximum radius
* ``rmaxpro`` - maximum projected radius
* ``rorder`` - indices of stars sorted by radius
* ``rproorder`` - indices of stars sorted by projected radius
* ``rm`` - half-mass radius
* ``rmpro`` - projected half-mass radius
* ``r10`` - 10% lagrange radius
* ``r10pro`` - projected 10% lagrange radius

If luminosities have been provided:

* ``rh`` - half-light radius
* ``rhpro`` - projected half-light radius
* ``rh10`` - radius containing inner 10% of the cluster's light
* ``rh10pro`` - projected radius containing inner 10% of the cluster's light


If the cluster's galactocentric positiion (xgc,ygc,zgc) and velocity (vxgc,vygc,vzgc) are known, then orbital information can be added via:

>>> cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

It is beneficial to use ``add_orbit`` as opposed to setting variables like ``cluster.xgc`` manually because additional arguments that can be passed. These arguments, and their default values, include ``ounits=None``, ``initialize=False``, ``ro=solar_ro``, and ``vo=solar_vo``. ``ounits`` informs ``StarCluster`` of the units of the galactocentric coordinates that are being provided. If they differ from ``cluster.units``, the unit conversion will be handled internally. ``intialize`` is an example of how strongly ``clustertools`` relies on ``galpy``, for if ``intialize=True`` a ``galpy`` orbit is initialized using the galactocentric coordinates provided, ro (distance from vantage point to the cluster (kpc)), vo (the circular velocity at ro (km/s)) and the Sun's motion (``solarmotion``). The ``galpy`` orbit can be accessed via ``cluster.orbit``.

``add_nbody6``, ``add_sse``, and ``add_bse`` are tailored to the standard output from ``NBODY6``. They are simple functions of convenience for adding information that ``NBODY6`` provides regarding the simulations itself, single star evolution, and binary star evolution. They would be called via:

>>> cluster.add_nbody6(nc, rc, rbar, rtide, xc, yc, zc, zmbar, vstar, rscale, nsbnd, nbbnd)
>>> cluster.add_sse(kw, logl, logr, ep, ospin)
>>> cluster.add_bse(id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,logr1,
                    logr2,ep1,ep2,ospin1,ospin2)

For those not familiar with ``NBODY6``, please consult the documention for ``add_nbody6``, ``add_sse``, and ``add_bse`` for the defintion of each variable. It is important to note that each of the above variables are intialized upon the initialization of ``StarCluster``, hence they can be set manually as well if you are using a code other than ``NBODY6`` and would like to define some of these parameters. Note, no units or origin are associated with any of the values provided via ``add_nbody6``, ``add_sse``, and ``add_bse`` such that they are not adjusted when unit and coordinate transformations are performed. 

One variable that is worth expanding on is the ``kw`` parameter. Motivated by NBODY6, each star's stellar evolution type is described by ``kw``. The below table illustrates what each ``kw`` integer represents. If you need to quickly lookup this table, it can be printed to screen using the function ``kwtypes()`` (See :ref:`Utilities <utilities>`).

.. list-table:: Table 1 - Relationship between ``kw`` and stellar evolution type, as per NBODY6
   :widths: 25 25
   :header-rows: 1

   * - KW
     - Stellar Evolution Types
   * - 0 
     - Low main sequence (M < 0.7).
   * - 1
     - Main sequence.
   * - 2
     - Hertzsprung gap (HG).
   * - 3
     - Red giant.
   * - 4
     - Core Helium burning.
   * - 5
     - First AGB.
   * - 6
     - Second AGB.
   * - 7
     - Helium main sequence.
   * - 8
     - Helium HG.
   * - 9
     - Helium GB.
   * - 10
     - Helium white dwarf.
   * - 11
     - Carbon-Oxygen white dwarf.
   * - 12
     - Oxygen-Neon white dwarf.
   * - 13
     - Neutron star.
   * - 14
     - Black hole.
   * - 15
     - Massless supernova remnant.

It is important to note that if stellar luminosities (``logl``) have been provided, calling ``cluster.analyze`` will also calculate the half-light radius of the cluster ``cluster.rh`` and the radius containing 10% of the light ``cluster.rh10``. The projected half-light radius ``cluster.rhpro`` and the projected radius containing 10% of the light ``cluster.rh10pro`` are calculated as well.

``add_energies`` is slightly more than a convenience function, because if the kinetic energy (kin) and potential energy (pot) of each star added to ``StarCluster`` via:

>>> cluster.add_energies(kin, pot)

then the total energy of each star (``cluster.etot``), total kinetic energy (``cluster.ektot``) and total potential energy (``cluster.ptot``) are calculated. Additionally, the virial parameter Qvir is calcualted and can be accessed via ``cluster.qvir``. 

Orbit actions must be added via ``add_actions`` because the associated variables are not created when a ``StarCluster`` is initialized. Hence once the actions JR, Jphi, and Jz have been calculated they can be added to the cluster via:

>>> cluster.add_actions(JR, Jphi, Jz)

Note that it is also possible to add orbital frequencies and periods by using:

>>> cluster.add_actions(JR, Jphi, Jz, OR, Ophi, Oz, TR, Tphi, Tz)

Finally, it is also possible within ``clustertools`` to define a subset of stars within the ``StarCluster`` via the ``subset`` function:

>>>cluster.subset(rmin=0.,rmax=1.,mmin=0.1,mmax=1.0,vmin=-5.,vmax=5.,emin=-100,rmax=0,kwmin=0,kwmax=1.,projected=False
)

where ``projected=True`` would use the projected radii and velocities to define the subset. Once called, ``cluster.indx`` will be a boolean array that is true only for stars that meet the criteria. Other parameters that can be passed to ``subset`` include npop (population number being studied) and indx (a custom boolean array).

Instead of using ``subset``, its also possible to extract a subset of stars from a ``StarCluster`` to form a new ``StarCluster`` using the ``sub_cluster`` function. For example, to extrct only stars within the cluster's half-mass radius one can call:

>>> new_cluster=ctools.sub_cluster(cluster,rmin=0,rmax=cluster.rm)

In this example, ``new_cluster`` will contain all the same information as ``cluster`` but only for stars within ``cluster.rm``. Please consult the ``sub_cluster`` documentation for the complete list of criteria that can be given to ``sub_cluster``.

Finally, it is also possible to determine where two different ``StarCluster``s overlap. This function is usefull when trying to compare two different systems. The comparison can be done via:

>>> overlap_cluster(cluster1,cluster2,tol=0.1,projected=False,return_cluster=True)

Here, ``tol`` represents the distance tolerance. If a star in ``cluster1`` is within ``tol`` of a star in ``cluster2`` then these regions of the cluster are deemed to overlap. If ``return_cluster=True`` then a ``sub_cluster`` of ``cluster1`` that overlaps with ``cluster2`` is returned. Otherwise a boolean array that is ``True`` for overlapping stars is returned.

.. automodapi:: clustertools.cluster.cluster
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading:
        :skip: StarCluster
        :noindex:

Units
-----

When a ``StarCluster`` is initialized, the default value of ``StarCluster.units`` is ``None``. However it is possible for users to specify the units system used by stars in the ``StarCluster``. For most functions, it is necessary to set units in order for calculations to be carried out. At present, ``clustertools`` supports 8 different string inputs for ``StarCluster.units``. The inputs and their meanings are summarized in Table 1 below.

.. list-table:: Table 1 - Units available in ``clustertools``
   :widths: 35 25 25 25 20 20
   :header-rows: 1

   * - Name
     - Definition
     - Distance Units
     - Velocity Units
     - Mass Units
     - Time Units
   * - nbody
     - nbody or Henon units, where G=M=rv=1
     - nbody
     - nbody
     - nbody
     - nbody
   * - pckms
     - units for when working in clustercentric coordinates
     - pc
     - km/s
     - Msun
     - Myr
   * - pcmyr
     - units for when working in clustercentric coordinates
     - pc
     - pc/Myr
     - Msun
     - Myr
   * - kpckms
     - units for working in galactocentric coordinates
     - kpc
     - km/s
     - Msun
     - Gyr
   * - kpcgyr
     - units for working in galactocentric coordinates
     - kpc
     - kpc/Gyr
     - Msun
     - Gyr
   * - radec
     - units for comparison to observations
     - degrees (Ra, Dec), and kpc (distance)
     - mas (proper motions) and kms (radial velocity)
     - Msun
     - Gyr
   * - galpy
     - galpy or natural units, set so the Sun orbits at a distance of 1 with a velocity of 1 (assumes ``ro,vo=8.275,239.2``)
     - kpc/ro
     - kms/vo
     - Msun/110119572536.69392
     - Gyr/0.03382094817762924
   * - WDunits
     - Walter Dehnen units used by NEMO
     - kpc
     - pc/Myr
     - Msun/222288.4543021174
     - Gyr

In order to convert between certain units, it is necessary to know the distance to the Galactic Centre (``ro``), the rotation velocity at ``ro`` (``vo``), the Sun's height above the Galactic disk (``zo``), and the Sun's motion with respect to the local standard of rest (``solarmotion``). When a ``StarCluster`` is initialzied, the default value for ``ro`` is 8.275 kpc (Gravity Collaboration, Abuter, R., Amorim, A., et al. 2020 ,A&A, 647, A59), the [U,V,W] motion of the Sun is set to [-11.1,12.24,7.25] (Schönrich, R., Binney, J., Dehnen, W., 2010, MNRAS, 403, 1829) and the velocity of the local standard of rest is 239.23 km/s. The choice of ``vo`` is such that ``vo`` + V is consistent with current estimates of the proper motion of Sagitarius A* (Reid, M.J. & Brunthaler, A., ApJ, 892, 1). The height of the Sun above the disk is 0.0208 kpc (Bennett, M. & Bovy, J. 2019, MNRAS, 483, 1417). However each of these values can be changed when intializing a ``StarCluster``.

See :ref:`Operations <operations>` for information on operations that convert a ``StarCluster`` from one set of units to another.

Coordinate Systems
------------------

Similar to ``StarCluster.units``, when a ``StarCluster`` is initialized the default value of ``StarCluster.origin`` is ``None``. However it is possible for users to specify the origin of the coordinate system used by stars in the ``StarCluster``. At present, ``clustertools`` supports 4 different string inputs for ``StarCluster.origin``. The inputs and their meanings are summarized in Table 2 below.

.. list-table:: Table 2 - Coordinate systems available in ``clustertools``
   :widths: 25 25
   :header-rows: 1

   * - Name
     - Definition
   * - cluster
     - clustercentric reference frame with origin equal to cluster's orbital position
   * - centre
     - clustercentric reference frame with origin equal to cluster's centre of density or centre of mass
   * - galaxy
     - galactocentric reference frame with origin at the centre of the galaxy
   * - sky
     - sky coordinate system, used when units are set to degrees

For clarity, it is worth expanding on the difference between ``StarCluster.origin='cluster'`` and ``StarCluster.origin='centre'``. The motivation for the two separate reference frames stems from codes like NBODY6 where the orbital evolution of the cluster in the external tidal field is handled separately from the internal cluster evolution. Hence the cluster's orbital position is integrated forwards in time from its initial conditions while the cluster's centre of density (or mass) will wander slightly due to internal cluster evolution. The true centre may in fact wander a lot if tidal tail stars are kept in the simulation or the cluster reaches dissolution. Codes like NBODY6 provide snapshots where the origin is at the cluster's orbital position (``StarCluster.origin='cluster'``) and then provide the location of the cluster's true centre separately to change to ``StarCluster.origin='centre'``.

In most cases, the cluster's centre may not be done. The centre can be determined using the ``find_centre()`` command, where the default option is to find the centre of density using the method of Harfst, S., Gualandris, A., Merritt, D., et al. 2007, NewA, 12, 357. Alternatively it is possible to find the centre of density using Casertano, S., Hut, P. 1985, ApJ, 298, 80 via ``find_centre(method='casertano')``. If the centre of mass is preferred, use ``find_centre(density=False)``

It is also important to note that converting to ``sky`` requires knowlendge of ``ro``, ``vo``, ``zo``, and ``solarmotion``, the default values of which are discussed above.


See :ref:`Operations <operations>` for information on operations that convert a ``StarCluster`` from one coordinate system to another.
