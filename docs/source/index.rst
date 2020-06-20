``clustertools``
================

``clustertools`` is a Python package with tools for analysing star cluster simulations. The package is built around the StarCluster class, for which all functions are able use or to act on. ``clustertools`` can be used for unit and coordinate transformations, the calculation of key structural and kinematic parameters, analysis of the cluster's orbit and tidal tails (with the help of `galpy
<https://docs.galpy.org/en/v1.6.0/index.html>`_, and measuring common cluster properties like its mass function, density profile, velocity dispersion profile (among others). 

The package contains methods for loading data from commonly used N-body codes (NBODY6, GYRFALCON, AMUSE) and generic snapshots. With the help of `limepy
<https://limepy.readthedocs.io/en/latest/>`_, it is also possible to generate cluster datasets from a distribution function for immediate analysis.

``clustertools`` is developed on Github. Please go to https://github.com/webbjj/clustertools to report issues or contribute to the code. 

Guide
-----------------

.. toctree::
   :maxdepth: 2

   installation.rst

   cluster.rst

   units_and_coordinate_systems.rst

   initialize.rst

   operations.rst

   functions.rst

   profiles.rst

   orbit.rst

   plotting.rst



Example List
-----------------

.. toctree::
   :maxdepth: 1

   notebooks/clustertools_example.ipynb

Planned Future Additions
----------------------------------

* Calculate meq and c_k from Biancinni et al.

* Caclulate tdiss from Baumgardt and Makino 2003

* Allow for the use of astropy units

* Ability to read NBODY6 binary files as opposed to forcing them to be in ascii format.

* Customizable loading of NBODY6 files to reflect the many customized versions of NBODY6 outthere (including NBODY6tt, NBODY6DF, and NBODY6++)

* Be able to initialize orbits and clusters based on Holger Baumgardt's Galactic Globular Clusters Database (https://people.smp.uq.edu.au/HolgerBaumgardt/globular/)

* Incorporate other methods for estimating a cluster's dynamical age (A+, Blue Stragglers) (REF)

* More analysis features for binary stars

* Allow for stars to be tagged as members of different subpopulations for the study of multiple populations

* More analysis features that are comparable to standard techniques used by observers

* Incorporation of a fast evolution code like in order to predict a clusters forward evolution

Library Reference
-----------------

.. toctree::
   :maxdepth: 6

   /reference/clustertools
