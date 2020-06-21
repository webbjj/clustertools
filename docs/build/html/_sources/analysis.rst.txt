Analysis
===============

Operations
----------------------------------

Several operations are available within ``clustertools`` to perform large scale changes to a ``StarCluster``, primarily related to unit and coordinate system changes. Operations are primarily meant to be called as functions internal to the ``StarCluser`` class (``StarCluster.operation_name()``), however it is also possible to call operations as external functions (``operation_name(StarCluster)`` for cases where the operation returns information of interest (i.e. a rotation matrix). No information is returned when an operation is called internally and no variables within cluster are set when an operation is called externally.

An example of the difference between internal and external operation calls would be finding a cluster's centre. One can call the function internally via:

>>> cluster.find_centre()

When called internally, the position (cluster.xc, cluster.yc, cluster.zc) and velocity (cluster.vxc, cluster.vyc, cluster.vzc) of the cluster's centre are set. Unless ``cluster.origin=='galaxy'``, then the cluster's centre is set equal to its galactocentric position and velocity (cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc). 

Alternatively, when called externally via:

>>> xc,yc,zc,vxc,vyc,vzv=find_centre(cluster)

no variables within ``cluster`` are actually set. 

A change of units or coordinate system can be implemented using:

>>> cluster.to_galaxy()
>>> cluster.to_kpckms()

It is important to note that only ``rv3d`` is called after a change in either units or coordinate system. By default, for computational efficiency, stars are not sorted again in the new coordinate system and key parameters are not re-calculated. If sorting and re-calculations are preferred, then be sure to set ``do_order=True`` and ``do_key_params=True`` via:

>>> cluster.to_galaxy(do_order=True)
>>> cluster.to_kpckms(do_key_params=True)

All available operations are listed below.

.. automodapi:: clustertools.analysis.operations
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading:

Functions
----------------------------------

A long list of functions exist within ``clustertools`` that can be called to measure a speciic property of the ``StarCluster``. Similar to operations, functions can be called internally (``StarCluster.function_name()``) or externally (``function_name(StarCluster)``). When called internally, changes are made to variables within ``StarCluster`` with nothing returned. When called externally, some operations have data returned. For example, a cluster's half-mass relaxation time can be called internally via:

>>> cluster.half_mass_relaxation_time()

where the variable cluster.trh now represents the cluster's half-mass relaxation time. Alternatively, the cluster's half-mass relaxation time can be called externally via:

>>> trh=half_mass_relaxation_time(cluster)

When called externally, no variables within `cluster` are set.

Some functions, including ``mass_function`` and ``eta_function`` can be applied to only a sub-population of the cluster. They have an array of input parameters that allow for only certain stars to be used. For example, to measure the mass function of stars between 0.1 and 0.8 Msun within the cluster's half-mass radius one can call:

>>> m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha = mass_function(cluster,
        mmin=0.1,mmax=0.8,rmin=0.,rmax=cluster.rm)

Other constraints include velocity (``vmin``,``vmax``), energy (``emin``,``emax``) and stellar evolution type (``kwmin``,``kwmax``). Alternatively a boolean array can be passed to ``index`` to ensure only a predefined subset of stars is used.

All functions can be called using projected valeus using ``projeted=True``.

The complete list of available functions are tabulated below in their external form.

.. automodapi:: clustertools.analysis.functions
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading:

Profiles
----------------------------------

A common measurement to make when analysing star clusters is to determine how a certain parameter varies with clustercentric distance. Therefore a list of commonly measured profiles can be measured via ``clustertools``. Unlike operations and functions, profiles can only be called externally. For example, the density profile of a cluster can be measured via:

>>> rprof, pprof, nprof= rho_prof(cluster)

where the radial bins ``rprof``, the density in each radial bin ``pprof``, and the number of stars in each bin ``nprof`` will be returned. 

One feature that is available for all profiles is the ability to plot the profile using ``plot=True``. For example, calling

>>> rprof, pprof, nprof= rho_prof(cluster,plot=True)

returns the below figure

.. image:: /images/rho_prof.png

Some editing can be done to the plot via key word arguments (see PLOTTING for information regarding making figures with ``clustertools``)

All profiles can also just be measured for a select subpopulation of stars based on mass (``mmin``,``mmax``), radius (``rmin``,``rmax``), velocity (``vmin``,``vmax``), energy (``emin``,``emax``) and stellar evolution type (``kwmin``,``kwmax``). Alternatively a boolean array can be passed to ``index`` to ensure only a predefined subset of stars is used.

All available profiles are listed below.

.. automodapi:: clustertools.analysis.profiles
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading:

Orbit
----------------------------------

In cases where the ``StarCluster`` does not evolve in isolation, it is possible to specify both the cluster's orbit and the the external tidal field. When orbital and tidal field information are provided, ``clustertools`` makes use of ``galpy`` to calcuate additional cluster and orbital properties.

**Cluster Properties Dependent on the External Tidal Field**

In cases where the ``StarCluster`` does not evolve in isolation, it is possible to specify both the cluster's orbit and the the external tidal field. When orbital and tidal field information are provided, ``clustertools`` makes use of ``galpy`` to calcuate additional cluster properties like its tidal radius (also referred to as the Jacobi radius) or limiting radius (also referred to as the King tidal radius). The defualt external potential is always the MWPotential2014 model from ``galpy``, which is an observationally motived represetnation of the Milky Way. However it is always possible to define a potential using the ``pot`` variable where applicable.

For example, the limiting radius of a cluster (radius where the density falls to background) can be measured via:

>> cluster.rlimiting(plot=True)

Since ``plot=True``, a figure showing the cluster's density profile and the background density in MWPotential at the cluster's position is returned in order to view exactly where the limiting radius is. cluster.rl is now set equal to the the measured limiting radius.

.. image:: /images/rlplot.png

Alternatively, one can define a different potential and make the calculation again:

>> cluster.rlimiting(pot=LogarithmicHaloPotential,plot=True)

.. image:: /images/rlplot_log.png

**Orbit Properties**

The majority of the orbital analysis performed in ``clusertools`` is just a wrapper around a ``galpy`` function that uses the ``StarCluster`` class, especially ``calc_actions``,``get_cluster_orbit``, and ``ttensor``. At the moment, orbital analysis can be done internally using the ``StarCluster`` class. Hence all orbit analysis function should be called externally. For example, one can easily initialize and integrate a cluster's orbit, given its galactocentric coordinates are known, using:

>> orbit=initialize_orbit(cluster)
>> dt,orbit=integrate_orbit(cluster,tfinal=12)

where ``tfinal=12`` means that the orbit is being integrated 12 Gyr into the future. The returned orbit is also stored in ``cluster.orbit`` while the array dt are the times (in galpy units) for which the orbit was integrated for. ``galpy`` orbits can be initialized and integrated for individual stars as well, via ``cluster.initialize_orbits()`` and ``cluster.integrate_orbits``, however this feature is likely only of interest for tail stars that are no longer bound to the cluster as orbit integration does not account for the cluster's potential.

It may also be helpful to have the cluster's orbital path within +/- dt, especially if looking at escaped stars. Arrays for the cluster's past and future coordinates can be generated via:

>>> t, x, y, z, vx, vy, vz = orbit_path(cluster,dt=0.1,plot=True)

Since ``plot=True`` and dt=0.1, a figure showing all stars in the simulation as well as the orbital path +/- 0.1 Gyr in time is returned.

.. image:: /images/rlplot_log.png

Stars can also be matched to a point along the orbital path, where each stars distance from the path and distance along the path to the progenitor is calculated. This function effectively creates a convenient frame of reference for studying escaped stars.

>>> tstar,dprog,dpath = orbit_path_match(cluster,dt=0.1,plot=True)

Since ``plot=True``, a figure showing ``dpath`` vs ``dprog`` is returned.


**Tidal Tails**

Since stars which escape a cluster do so with velocities slightly larger or smaller than the progenitor cluster, tidal tails will not follow the cluster's orbit path. Therefore ``clustertools`` can generate a tail path by first matching stars to the cluster's orbital path (as above) and then binning stars based on their distance from the progenitor along the path. The mean coordinates of stars in each bin marks the tail path. Similar to a cluster's orbital path, ``clustertools`` can also match stars to the tail path as well.

>>> t, x, y, z, vx, vy, vz = tail_path(cluster,dt=0.1,plot=True)

.. image:: /images/rlplot_log.png

>>> tstar,dprog,dpath = tail_path_match(cluster,dt=0.1,plot=True)

.. image:: /images/rlplot_log.png

A complete list of all functions related to a cluster's orbit can be found below.

.. automodapi:: clustertools.analysis.orbit
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading:
