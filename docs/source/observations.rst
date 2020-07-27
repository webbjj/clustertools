Observations
===============

Currently, the majority of ``clustertools`` has been written with full N-body simulations of star clusters in mind where at minimum the three dimensional position and velocity of every stars is known. However, the majority of the analysis can be performed on an observed dataset as well where less than 6 dimensions worth of data is known. When setting up an observed ``StarCluster`` using observational data using:

>>> cluster=StarCluster(units='radec',origin='sky',projected=True)
>>> cluster.add_stars(ra,dec,dist,pmra,pmdec,vlos)

simply set any of the variables to zero if they are unknown. Note that even if ``projected=False``, setting ``dist=0`` and ``vlos=0`` will result in three deminstional values being projected.

Furthermore, a few functions have options that have been included with observations in mind. For example, if one wishes to move the above cluster to a clustercentric coordinate system, first define the cluster's orbital centre:

>>> cluster.add_orbit(ra_gc,dec_gc,dist_gc,pmra_gc,pmdec_gc,vlos_gc)

then move to the cluster's centre using an observational method of choice. The default is just to find each stars angular separation from the centre. However its also possible to move the cluster to orthographic coordinates via Helmi et al 2018:

>>> cluster.to_cluster(centre_method='orthographic')

or using the centring method of Van de Ven et al. 2005

>>> cluster.to_cluster(centre_method='VandeVen')

Finally, when measuring the mass function of an observed dataset, there is an additional parameter ``mcorr`` that allows for stellar masses to be corrected for completeness. Simple call:

>> cluster.mass_function(mmin=0.1,mmax=0.8,rmin=0.,rmax=cluster.rm,mcorr=mcorr)

where ``mcorr`` is an array of estimated completeness values of each stars in the cluster. 

More support for working with observed datasets is planned for future releases.


