``clustertools``
================

``clustertools`` is a Python package with tools for analysing star cluster simulations. The package is built around the StarCluster class, for which all functions are able to use or act on. ``clustertools`` can be used for unit and coordinate transformations, the calculation of key structural and kinematic parameters, analysis of the cluster's orbit and tidal tails (with the help of `galpy
<https://docs.galpy.org/en/v1.6.0/index.html>`_), and measuring common cluster properties like its mass function, density profile, and velocity dispersion profile (among others). 

The package contains methods for loading data from commonly used N-body codes, generic snapshots, and codes for generating intial conditions. 

``clustertools`` is developed on Github. Please go to https://github.com/webbjj/clustertools to report issues or contribute to the code. 

Supported N-Body Simulation Codes, Snapshot Formats, and Packages
-----------------

AMUSE <https://amusecode.github.io/> (Portegies Zwart S., McMillan S., 2018, Astrophysical Recipes; The art ofAMUSE, doi:10.1088/978-0-7503-1320-9)

ASTROPY_TABLE <https://docs.astropy.org/en/latest/> (Astropy Collaboration 2013, A&A, 558, A33)

GALPY <https://docs.galpy.org/en/v1.7.2/index.html> (Bovy J., 2015, ApJS, 216, 29)

LIMEPY <https://limepy.readthedocs.io/en/latest/> (Gieles, M. & Zocchi, A. 2015, MNRAS, 454, 576)

NBODY6 <https://people.ast.cam.ac.uk/~sverre/web/pages/nbody.htm> (Aarseth S. J., 2003, Gravitational N-Body Simulations)

NBODY6++GPU <https://github.com/nbodyx/Nbody6ppGPU> (Wang L., Spurzem R., Aarseth S., Nitadori K., Berczik P., Kouwenhoven M. B. N., Naab T., 2015, MNRAS, 450, 4070)

NEMO/GYRFALCON <https://astronemo.readthedocs.io/en/latest/> (Teuben P., 1995, The Stellar Dynamics Toolbox NEMO. p. 398)

Guide
-----------------

.. toctree::
   :maxdepth: 2

   installation.rst

   cluster.rst

   load.rst

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

Acknowledging clustertools
--------------------------

If you use ``clustertools`` in a publiclication, please cite the following digital object identifier (DOI):

.. image:: https://zenodo.org/badge/272233602.svg
   :target: https://zenodo.org/badge/latestdoi/272233602

and link to ``https://github.com/webbjj/clustertools``. 

If you perform any analysis related to a cluster's orbit and the external tidal field, please cite ``galpy`` (Bovy 2015).

If you use the ``setup_cluster`` command, please cite ``limepy`` (Gieles & Zocchi 2015) or ``galpy`` (Bovy 2015) depending on your choice of  ``ctype``. If you use ``setup_cluster`` to initialize a specific Galactic globular cluster, please cite either de Boer et al. (2019) or Harris (1999, 2010 Edition) depending on the cluster and your choice for the keyword argument ``source``. If you make use of the cluster's orbit, please cite Vasiliev (2019) as orbital properties are taken from that paper.

For specific functions, it is recomended to check that specific functions documentation to check for a specific reference if applicable. The list below will continually be updated when possible:

References
-----------------

``alpha_prof`` - Webb, J.J. & Vesperini, E. 2016, MNRAS, 463, 2383

``core_relaxation_time`` - Stone, N.C. & Ostriker, J.P. 2015, ApJ, 806, 28

``ckin`` - Bianchini et al. 2018, MNRAS, 475, 96

``cyl_coords`` ,``cart_to_cy``,``cyl_to_cart`` - Bovy J., 2015, ApJS, 216, 29

``find_centre(density=True)`` (default) - Harfst, S., Gualandris, A., Merritt, D., et al. 2007, NewA, 12, 357

``half_mass_relaxation_time`` - Spitzer, L. 1987, Dynamical evolution of globular clusters

``initialize_orbit``, ``initialize_orbits``, ``interpolate_orbit``, ``interpolate_orbits`` - Bovy J., 2015, ApJS, 216, 29

``load_cluster(ctype='limepy')`` - Gieles, M. & Zocchi, A. 2015, MNRAS, 454, 576

``load_cluster('limepy','gcname')`` (default) - Gieles, M. & Zocchi, A. 2015, MNRAS, 454, 576 - de Boer, T. J. L., Gieles, M., Balbinot, E., Hénault-Brunet, V., Sollima, A., Watkins, L. L., Claydon, I. 2019, MNRAS, 485, 4906 - Vasiliev E., 2019, MNRAS, 484,2832

``load_cluster('limepy',gcname', source='harris')`` - Harris, W.E. 1996 (2010 Edition), AJ, 112, 1487 - Gieles, M. & Zocchi, A. 2015, MNRAS, 454, 576 - Vasiliev E., 2019, MNRAS, 484,2832

``meq_function``, ``meq_prof`` - Bianchini, P. et al. 2016, MNRAS, 458, 3644 

``orbital_path``, ``orbital_path_match`` - Bovy J., 2015, ApJS, 216, 29 -

``relaxation_time`` - Spitzer, L. Jr, Hart, M.H. 1971, ApJ, 164, 399

``rlimiting`` - Bovy J., 2015, ApJS, 216, 29 

``rtidal`` - Bertin, G. & Varri, A.L. 2008, ApJ, 689, 1005 - Bovy J., 2015, ApJS, 216, 29 - Webb, J.J., Bovy, J., Carlberg, R.G., Gieles, M. 2019, MNRAS, 448, 4

``tail_path``, ``tail_path_match`` - Bovy J., 2015, ApJS, 216, 29

``to_centre(method='orthonormal')`` - GAIA Collaboration, 2018, A&A, 616, A12 

``to_centre(method='VandeVen')`` - van de Ven, G. 2005, PhD Thesis, Leiden University

``to_cluster(method='orthonormal')`` - GAIA Collaboration, 2018, A&A, 616, A12 

``to_cluster(method='VandeVen')`` - van de Ven, G. 2005, PhD Thesis, Leiden University

``to_sky``,``sky_coords``,``cart_to_sky`` - Bovy J., 2015, ApJS, 216, 29 

``virial_radius(method='critical_density')`` - Bovy J., 2015, ApJS, 216, 29 -

``virial_radius(method=inverse_distance)`` - Portegies Zwart S., McMillan S., 2018, Astrophysical Recipes; The art ofAMUSE, doi:10.1088/978-0-7503-1320-9

Library Reference
-----------------

.. toctree::
   :maxdepth: 6

   /reference/clustertools

Planned Future Additions
----------------------------------

* Add PeTar to list of readable codes (in progress)

* Allow for the use of astropy units

* Be able to initialize orbits and clusters based on Holger Baumgardt's Galactic Globular Clusters Database (https://people.smp.uq.edu.au/HolgerBaumgardt/globular/)

* Incorporate other methods for estimating a cluster's dynamical age (A+, Blue Stragglers)

* More analysis features for binary stars

* Allow for stars to be tagged as members of different subpopulations for the study of multiple populations

* More analysis features that are comparable to standard techniques used by observers

* Include options for planetary systems around individual stars

* A simulation class that reads in all timesteps of a simulation
