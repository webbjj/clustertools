.. _units_and_coordinate_systems:

Units and Coordinate Systems
==============================

As discussed in :ref:`Cluster <cluster>`, when initializing a ``StarCluster`` the variables ``units`` and ``origin`` are by default set equal to ``None``. This approach means that any calculated values are done so in whatever units and coordinate system the data is provided in. However, if the units and coordinate system are known, specificying these parameters allows for unit and coordinate transformations to be done on the data. 

Units
------------

When a ``StarCluster`` is initialized, the default value of ``StarCluster.units`` is ``None``. However it is possible for users to specify the units system used by stars in the ``StarCluster``. At present, ``clustertools`` supports 6 different string inputs for ``StarCluster.units``. The inputs and their meanings are summarized in Table 1 below.

.. list-table:: Table 1 - Units available in ``clustertools``
   :widths: 25 25 25 25 25
   :header-rows: 1

   * - Name
     - Definition
     - Distance Units
     - Velocity Units
     - Mass Units
   * - nbody
     - nbody or Henon units, where G=M=rv=1
     - nbody
     - nbody
     - nbody
   * - pckms
     - units for when working in clustercentric coordinates
     - pc
     - km/s
     - Msun
   * - kpckms
     - units for working in galactocentric coordinates
     - kpc
     - km/s
     - Msun
   * - degrees
     - units for comparison to observations
     - degrees (Ra, Dec), and kpc (distance)
     - mas (proper motions) and kms (radial velocity)
     - Msun  
   * - galpy
     - galpy or natural units, set so the Sun orbits at a distance of 1 with a velocity of 1
     - kpc/ro
     - kms/vo
     - Msun/90027307126.905106

See :ref:`Operations <operations>` for information on operations that convert a ``StarCluster`` from one set of units to another.

Coordinate Systems
------------------------

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

For clarity, it is worth expanding on the difference between ``StarCluster.origin='cluster'`` and ``StarCluster.origin='centre'``. The motivation for the two separate reference frames stems from codes like NBODY6 where the orbital evolution of the cluster in the external tidal field is handled separately from the internal cluster evolution. Hence the cluster's orbital position is integrated forwards in time from its initial conditions while the cluster's centre of density (or mass) will wander slightly due to internal cluster evolution. The true centre may in fact wander a lot of tidal tail stars are kept in the simulation or the cluster reaches dissolution. Codes like NBODY6 provide snapshots where the origin is at the cluster's orbital position (``StarCluster.origin='cluster'``) and then provide the location of the cluster's true centre separately to change to ``StarCluster.origin='centre'``.


See :ref:`Operations <operations>` for information on operations that convert a ``StarCluster`` from one coordinate system to another.
