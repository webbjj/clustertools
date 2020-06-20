Operations
===============

Several operations are available within ``clustertools`` to perform large scale changes to a ``StarCluster``, primarily related to unit changes and coordinate system changes. Operations can be called as either an internal function (``StarCluster.operation_name()``) or an external function (``operation_name(StarCluster)``). When called internally, changes are made to the ``StarCluster`` with nothing returned. When called externally, some operations have data returned. 

For example, in order to find the cluster's centre of density one can call the function internally via:

>>> cluster.find_centre()

When called internally, the position (cluster.xc, cluster.yc, cluster.zc) and velocity (cluster.vxc, cluster.vyc, cluster.vzc) of the cluster's centre are set. Alternatively, when called externally via:

>>> xc,yc,zc,vxc,vyc,vzv=find_centre(cluster)

When called externally, no variables within ``cluster`` are set.

All available operations are listed below.

Operation List
-----------------

.. automodapi:: clustertools.main.operations
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading:

