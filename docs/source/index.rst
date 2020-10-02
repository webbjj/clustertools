``clustertools``
================

``clustertools`` is a Python package with tools for analysing star cluster simulations. The package is built around the StarCluster class, for which all functions are able to use or act on. ``clustertools`` can be used for unit and coordinate transformations, the calculation of key structural and kinematic parameters, analysis of the cluster's orbit and tidal tails (with the help of `galpy
<https://docs.galpy.org/en/v1.6.0/index.html>`_), and measuring common cluster properties like its mass function, density profile, and velocity dispersion profile (among others). 

The package contains methods for loading data from commonly used N-body codes (NBODY6, GYRFALCON, AMUSE) and generic snapshots. With the help of `limepy
<https://limepy.readthedocs.io/en/latest/>`_, it is also possible to generate cluster datasets from a distribution function for immediate analysis.

``clustertools`` is developed on Github. Please go to https://github.com/webbjj/clustertools to report issues or contribute to the code. 

Guide
-----------------

.. toctree::
   :maxdepth: 2

   installation.rst

   cluster.rst

   initialize.rst

   analysis.rst

   utilities.rst

   observations.rst

   binaries.rst

Example Notebooks
-----------------

.. toctree::
   :maxdepth: 1

   notebooks/loading_and_advancing.ipynb
   notebooks/setup.ipynb
   notebooks/operations.ipynb
   notebooks/functions.ipynb
   notebooks/profiles.ipynb
   notebooks/orbits_and_tails.ipynb

Planned Future Additions
----------------------------------

* Allow for the use of astropy units

* Be able to initialize orbits and clusters based on Holger Baumgardt's Galactic Globular Clusters Database (https://people.smp.uq.edu.au/HolgerBaumgardt/globular/)

* Incorporate other methods for estimating a cluster's dynamical age (A+, Blue Stragglers)

* More analysis features for binary stars

* Allow for stars to be tagged as members of different subpopulations for the study of multiple populations

* More analysis features that are comparable to standard techniques used by observers

References
-----------------

Bertin, G. & Varri, A.L. 2008, ApJ, 689, 1005

Bianchini, P. et al. 2016, MNRAS, 458, 3644

Bianchini et al. 2018, MNRAS, 475, 96

Bovy J., 2015, ApJS, 216, 29

Claydon, I., Gieles, M., Varri, A.L., Heggie, D.C., Zocchi, A. 2019, MNRAS, 487, 147

de Boer, T. J. L., Gieles, M., Balbinot, E., Hénault-Brunet, V., Sollima, A., Watkins, L. L., Claydon, I. 2019, MNRAS, 485, 4906 and

Gieles, M. & Zocchi, A. 2015, MNRAS, 454, 576

Harfst, S., Gualandris, A., Merritt, D., et al. 2007, NewA, 12, 357

Harris, W.E. 1996 (2010 Edition), AJ, 112, 1487

Portegies Zwart S., McMillan S., 2018, Astrophysical Recipes; The art ofAMUSE, doi:10.1088/978-0-7503-1320-9

Spitzer, L. Jr, Hart, M.H. 1971, ApJ, 164, 399

Spitzer, L. 1987, Dynamical evolution of globular clusters

van de Ven, G. 2005, PhD Thesis, Leiden University, https://ui.adsabs.harvard.edu/abs/2005PhDT........12V/abstract

Vasiliev E., 2019, MNRAS, 484,2832

Webb, J.J. & Vesperini, E. 2016, MNRAS, 463, 2383


Library Reference
-----------------

.. toctree::
   :maxdepth: 6

   /reference/clustertools
