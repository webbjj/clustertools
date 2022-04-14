``clustertools``
================

``clustertools`` is a Python package with tools for analysing star cluster simulations. The package is built around the StarCluster class, for which all functions are able to use or act on. ``clustertools`` can be used for unit and coordinate transformations, the calculation of key structural and kinematic parameters, analysis of the cluster's orbit and tidal tails (with the help of `galpy
<https://docs.galpy.org/en/v1.6.0/index.html>`_), and measuring common cluster properties like its mass function, density profile, and velocity dispersion profile (among others). 

The package contains methods for loading data from commonly used N-body codes, generic snapshots and codes for generating intial conditions. 

``clustertools`` is developed on Github. Please go to https://github.com/webbjj/clustertools to report issues or contribute to the code. 

Supported N-Body Simulation Codes
-----------------

AMUSE <https://amusecode.github.io/>

NBODY6 <https://people.ast.cam.ac.uk/~sverre/web/pages/nbody.htm>

NBODY6++GPU <https://github.com/nbodyx/Nbody6ppGPU>

NEMO/GYRFALCON <https://astronemo.readthedocs.io/en/latest/>

Supported Initial Condition Generation Codes
-----------------

galpy <https://docs.galpy.org/en/v1.6.0/index.html> - Very simple version, will be updated soon with Galpy's own distribution function generator

limepy <https://limepy.readthedocs.io/en/latest/>

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

Aarseth S. J., 2003, Gravitational N-Body Simulations - ``load_cluster('nbody6')``, ``load_cluster('nbody6se')``

Astropy Collaboration 2013, A&A, 558, A33 - ``load_cluster('astropy_table')``

Bertin, G. & Varri, A.L. 2008, ApJ, 689, 1005 - ``rtidal``

Bianchini, P. et al. 2016, MNRAS, 458, 3644 - ``meq_function``

Bianchini et al. 2018, MNRAS, 475, 96 - ``ckin``

Bovy J., 2015, ApJS, 216, 29 - ``initialize_orbit``, ``initialize_orbits``, ``integrate_orbit``, ``integrate_orbits``, ``orbit_interpolation``, ``orbital_path``, ``orbital_path_match``, ``calc_actions``, ``ttensor``, ``tail_path``, ``tail_path_match``, ``rtidal``, ``rlimiting``,``virial_radius(method='critical_density')``

Claydon, I., Gieles, M., Varri, A.L., Heggie, D.C., Zocchi, A. 2019, MNRAS, 487, 147

de Boer, T. J. L., Gieles, M., Balbinot, E., Hénault-Brunet, V., Sollima, A., Watkins, L. L., Claydon, I. 2019, MNRAS, 485, 4906 - ``setup_cluster('limepy','gcname')`` (default)

GAIA Collaboration, 2018, A&A, 616, A12 - ``to_cluster(method='orthonormal')``

Gieles, M. & Zocchi, A. 2015, MNRAS, 454, 576 - ``setup_cluster(ctype='limepy')``

Harfst, S., Gualandris, A., Merritt, D., et al. 2007, NewA, 12, 357 - ``find_centre(density=True)`` (default)

Harris, W.E. 1996 (2010 Edition), AJ, 112, 1487 - ``setup_cluster('gcname', source='harris')``

Portegies Zwart S., McMillan S., 2018, Astrophysical Recipes; The art ofAMUSE, doi:10.1088/978-0-7503-1320-9 - ``load_cluster('amuse'``,``virial_radius(method=inverse_distance)``

Spitzer, L. Jr, Hart, M.H. 1971, ApJ, 164, 399 - ``relaxation_time``

Spitzer, L. 1987, Dynamical evolution of globular clusters - ``half_mass_relaxation_time``

Stone, N.C. & Ostriker, J.P. 2015, ApJ, 806, 28 - ``core_relaxation_time``

Teuben P., 1995, The Stellar Dynamics Toolbox NEMO. p. 398 = ``load_cluster('nemo')``, ``load_cluster('gyrfalcon')``

van de Ven, G. 2005, PhD Thesis, Leiden University, https://ui.adsabs.harvard.edu/abs/2005PhDT........12V/abstract - ``to_cluster(method='VandeVen')``

Vasiliev E., 2019, MNRAS, 484,2832 - ``setup_cluster('limepy','gcname')``

Wang L., Spurzem R., Aarseth S., Nitadori K., Berczik P., Kouwenhoven M. B. N., Naab T., 2015, MNRAS, 450, 4070 - ``load_cluster('nbody6pp')``, ``load_cluster('nbody6++')``

Webb, J.J. & Vesperini, E. 2016, MNRAS, 463, 2383 - ``alpha_prof``


Library Reference
-----------------

.. toctree::
   :maxdepth: 6

   /reference/clustertools

Planned Future Additions
----------------------------------

* Allow for the use of astropy units

* Be able to initialize orbits and clusters based on Holger Baumgardt's Galactic Globular Clusters Database (https://people.smp.uq.edu.au/HolgerBaumgardt/globular/)

* Incorporate other methods for estimating a cluster's dynamical age (A+, Blue Stragglers)

* More analysis features for binary stars

* Allow for stars to be tagged as members of different subpopulations for the study of multiple populations

* More analysis features that are comparable to standard techniques used by observers
