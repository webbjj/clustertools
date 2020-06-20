Functions
===============

A long list of functions exist within ``clustertools`` that can be called to measure a speciic property of the ``StarCluster``. Similar to :ref:`operations <operations>`, functions can be called as either an internal function (``StarCluster.function_name()``) or an external function (``function_name(StarCluster)``). When called internally, changes are made to the ``StarCluster`` with nothing returned. When called externally, some operations have data returned. For example, a cluster's half-mass relaxation time can be called internally via:

>>> cluster.half_mass_relaxation_time()

where the variable cluster.trh now represents the cluster's half-mass relaxation time. Alternatively, the cluster's half-mass relaxation time can be called externally via:

>>> trh=half_mass_relaxation_time(cluster)

When called externally, no variables within `cluster` are set.

All available functions are listed below in their external form.

Function List
-----------------

.. automodapi:: clustertools.main.functions
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading: