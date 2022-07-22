``clustertools``
================

``clustertools`` is a Python package for analysing star cluster simulations. The package is built around the ``StarCluster`` class, which will store all the necessary information from a given star cluster simulation to be used for anaylsis. All functions within ``clustertools`` are then designed to act on a ``StarCluster``. ``clustertools`` can be used for unit and coordinate transformations, the calculation of key structural and kinematic parameters, analysis of the cluster's orbit and tidal tails (with the help of `galpy <https://docs.galpy.org/en/v1.6.0/index.html>`_ , and measuring common cluster properties like its mass function, density profile, and velocity dispersion profile (among others). While originally designed with star clusters in mind, ``clustertools`` can be used to study other types of N-body systems, including stellar streams and dark matter sub-halos.  

The package contains functions for loading data from commonly used N-body codes, generic snapshots, and codes for generating initial conditions. 

``clustertools`` is developed on Github. Please go to https://github.com/webbjj/clustertools to report issues or contribute to the code. 

Supported N-Body Simulation Codes, Snapshot Formats, and Packages
-----------------

`AMUSE <https://amusecode.github.io/>`_ (Portegies Zwart S., McMillan S., 2018, Astrophysical Recipes; The art ofAMUSE, doi:10.1088/978-0-7503-1320-9)

`ASTROPY_TABLE <https://docs.astropy.org/en/latest/>`_ (Astropy Collaboration 2013, A&A, 558, A33)

`GALPY <https://docs.galpy.org/en/v1.7.2/index.html>`_ (Bovy J., 2015, ApJS, 216, 29)

`LIMEPY <https://limepy.readthedocs.io/en/latest/>`_ (Gieles, M. & Zocchi, A. 2015, MNRAS, 454, 576)

`NBODY6 <https://people.ast.cam.ac.uk/~sverre/web/pages/nbody.htm>`_ (Aarseth S. J., 2003, Gravitational N-Body Simulations)

`NBODY6ppGPU <https://github.com/nbodyx/Nbody6ppGPU>`_ (Wang L., Spurzem R., Aarseth S., Nitadori K., Berczik P., Kouwenhoven M. B. N., Naab T., 2015, MNRAS, 450, 4070)

`NEMO <https://astronemo.readthedocs.io/en/latest/>`_ (Teuben P., 1995, The Stellar Dynamics Toolbox NEMO. p. 398)

Version 1.0 (May 6, 2022)
-----------------
Thanks to anyone that has been using ``clustertools`` while it was in devlopment. I am happy to announce that Version 1.0 has now officially been released.  A couple of minor changes that are worth noting include:

-- For consistancy purposes, all profile calls now return the true values of the radial bins and the calculated value. Pre-Version 1.0, some profile calls returned the natual logarithm of one (or both) values, and sometimes normalized by the cluster's effective radius. Any other returned values are returned as per their definition. For example, ``delta_alpha`` is still calculated as d(alpha)/d(ln(r/r_m)) even though calling ``alpha_prof`` returns ``r`` and ``alpha``.

-- ``load_cluster`` now requires ``units`` and ``origin`` to be set, with the exception of when ``ctype=nbody6`` or ``ctype=nbody6pp``. No assumptions are made regarding the format of incoming data

-- The ``setup_cluster`` function has been deprecated. You can now initialize a Galactic cluster using ``load_cluster(ctype='limepy',gcname='Pal5')``. If you wish to initialize a cluster with a specific lowered isothermal model, you must do so in ``limepy`` first and then call ``load_cluster`` accordingly.

-- Previously when centreing the cluster to orthonormal coordinates, positions were in radians and velocities were converted to radians/year, but not actually moved to the cluster's centre. Positions are now returned to degrees and velocities remain in mas/yr. Furthermore the velocities are shifted to be relative to the cluster's centre

The above changes were made to ensure consistency accross ``clustertools``. Now that Version 1.0 is official, I will work to ensure that any future changes will be compatabile with past versions of ``clustertools``.

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

For specific functions, please consult the list below for the appropriate reference (if applicable).

References
-----------------

``alpha_prof`` - Webb, J.J. & Vesperini, E. 2016, MNRAS, 463, 2383 (`ADS <https://ui.adsabs.harvard.edu/abs/2016MNRAS.463.2383W/abstract>`_)

``core_relaxation_time`` - Stone, N.C. & Ostriker, J.P. 2015, ApJ, 806, 28 (`ADS <https://ui.adsabs.harvard.edu/abs/2015ApJ...806L..28S/abstract>`_)

``ckin`` - Bianchini, P. et al. 2018, MNRAS, 475, 96 (`ADS <https://ui.adsabs.harvard.edu/abs/2018MNRAS.475L..96B/abstract>`_)

``cyl_coords`` ,``cart_to_cy``,``cyl_to_cart`` - Bovy J., 2015, ApJS, 216, 29 (`ADS <https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract>`_)

``find_centre(density=True)`` (default) - Harfst, S., Gualandris, A., Merritt, D., et al. 2007, NewA, 12, 357 (`ADS <https://ui.adsabs.harvard.edu/abs/2007NewA...12..357H/abstract>`_) or Casertano, S., Hut, P. 1985, ApJ, 298, 80 (`ADS <https://ui.adsabs.harvard.edu/abs/1985ApJ...298...80C/abstract>`_)

``half_mass_relaxation_time`` - Spitzer, L. 1987, Dynamical evolution of globular clusters (`ADS <https://ui.adsabs.harvard.edu/abs/1987degc.book.....S/abstract>`_)

``initialize_orbit``, ``initialize_orbits``, ``interpolate_orbit``, ``interpolate_orbits`` - Bovy J., 2015, ApJS, 216, 29 (`ADS <https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract>`_)

``load_cluster('limepy','gcname')`` (default) - Bovy J., 2015, ApJS, 216, 29 (`ADS <https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract>`_) - Gieles, M. & Zocchi, A. 2015, MNRAS, 454, 576 (`ADS <https://ui.adsabs.harvard.edu/abs/2015MNRAS.454..576G/abstract>`_) - de Boer, T. J. L., Gieles, M., Balbinot, E., Hénault-Brunet, V., Sollima, A., Watkins, L. L., Claydon, I. 2019, MNRAS, 485, 4906 (`ADS <https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.4906D/abstract>`_) - Vasiliev E., 2019, MNRAS, 484,2832 (`ADS <https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.2832V/abstract>`_)

``load_cluster('limepy',gcname', source='harris')`` - Harris, W.E. 1996 (2010 Edition), AJ, 112, 1487 - Bovy J., 2015, ApJS, 216, 29 (`ADS <https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract>`_) - Gieles, M. & Zocchi, A. 2015, MNRAS, 454, 576 (`ADS <https://ui.adsabs.harvard.edu/abs/2015MNRAS.454..576G/abstract>`_) - Vasiliev E., 2019, MNRAS, 484,2832 (`ADS <https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.2832V/abstract>`_)

``meq_function``, ``meq_prof`` - Bianchini, P. et al. 2016, MNRAS, 458, 3644 (`ADS <https://ui.adsabs.harvard.edu/abs/2016MNRAS.458.3644B/abstract>`_)

``orbital_path``, ``orbital_path_match`` - Bovy J., 2015, ApJS, 216, 29 (`ADS <https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract>`_)

``relaxation_time`` - Spitzer, L. Jr, Hart, M.H. 1971, ApJ, 164, 399 (`ADS <https://ui.adsabs.harvard.edu/abs/1971ApJ...164..399S/abstract>`_)

``rcore`` - Casertano, S., Hut, P. 1985, ApJ, 298, 80 (`ADS <https://ui.adsabs.harvard.edu/abs/1985ApJ...298...80C/abstract>`_)

``rlimiting`` - Bovy J., 2015, ApJS, 216, 29 (`ADS <https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract>`_)

``rtidal`` - Bertin, G. & Varri, A.L. 2008, ApJ, 689, 1005 - Bovy J., 2015, ApJS, 216, 29 (`ADS <https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract>`_) - Webb, J.J., Bovy, J., Carlberg, R.G., Gieles, M. 2019, MNRAS, 448, 4 (`ADS <https://ui.adsabs.harvard.edu/abs/2019MNRAS.488.5748W/abstract>`_)

``tail_path``, ``tail_path_match``, ``to_tail`` - Bovy J., 2015, ApJS, 216, 29 (`ADS <https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract>`_)

``to_centre(method='orthonormal')`` - GAIA Collaboration, 2018, A&A, 616, A12 (`ADS <https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..12G/abstract>`_)

``to_centre(method='VandeVen')`` - van de Ven, G. 2005, PhD Thesis, Leiden University (`ADS <https://ui.adsabs.harvard.edu/abs/2005PhDT........12V/abstract>`_)

``to_cluster(method='orthonormal')`` - GAIA Collaboration, 2018, A&A, 616, A12 (`ADS <https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..12G/abstract>`_)

``to_cluster(method='VandeVen')`` - van de Ven, G. 2005, PhD Thesis, Leiden University (`ADS <https://ui.adsabs.harvard.edu/abs/2005PhDT........12V/abstract>`_)

``to_sky``,``sky_coords``,``cart_to_sky`` - Bovy J., 2015, ApJS, 216, 29 (`ADS <https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract>`_)

``virial_radius(method='critical_density')`` - Bovy J., 2015, ApJS, 216, 29 (`ADS <https://ui.adsabs.harvard.edu/abs/2015ApJS..216...29B/abstract>`_)

``virial_radius(method=inverse_distance)`` - Portegies Zwart S., McMillan S., 2018, Astrophysical Recipes; The art ofAMUSE, doi:10.1088/978-0-7503-1320-9 (`ADS <https://ui.adsabs.harvard.edu/abs/2018araa.book.....P/abstract>`_)

Contributers
-----------------
Jo Bovy - https://github.com/jobovy
Erik Gillis - https://github.com/e-gillis
Nathaniel Starkman - https://github.com/nstarman

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
