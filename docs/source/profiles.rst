Profiles
===============

A common measurement to make when analysing star clusters is to determine how a certain parameter varies with clustercentric distance. Therefore a list of commonly measured profiles can be measured via ``clustertools``. Unlike :ref:`Operations <operations>` and :ref:`Functions <functions>`, profiles can only be called externally. For example, the density profile of a cluster can be measured via:

>>> rprof, pprof, nprof= rho_prof(cluster)

where the radial bins ``rprof``, the density in each radial bin ``pprof``, and the number of stars in each bin ``nprof`` will be returned. 

One feature that is available for all profiles is the ability to plot the profile using ``plot=True``. For example, calling

>>> rprof, pprof, nprof= rho_prof(cluster,plot=True)

returns the below figure

FIGURE

Some editing can be done to the plot via key word arguments (see PLOTTING for information regarding making figures with ``clustertools``)

All available profiles are listed below.

Profile List
----------------------------------

.. automodapi:: clustertools.main.profiles
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading: