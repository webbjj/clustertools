Binaries
===============

At the moment, NBODY6 and NBODY6++GPU are the only supported codes that have a separate treatment for binary stars. In both cases, the main array of stars contain the centre of mass position and velocity of each binary. The masses are then the total mass of the binary star. Separate parameters are used that contain properties of the individual binaries themselves. These include, but are not limited to:

* ids ``i_d1`` and ``i_d2``
* kwtypes ``kw1`` and ``kw2``
* eccentricity ``ecc``
* orbital period ``pb``
* semi-major axis ``semi``
* masses ``m1`` and ``m2``
* luminosities ``logl1`` and ``logl2``

It is possible to add binary star properties manually using ``StarCluster.add_bse``, however ``clustertools`` will always assume the binary's centre of mass positions and velocities, total mass, and primary star ID were added via ``StarCluster.add_stars``.

Generally speaking, unit and coordinate transformations do not act on these parameters. They are saved in the event an experienced NBODY6 user would like to view these parameters (see ``StarCluster.add_bse`` for the full list of saved parameters). The only functions designed specifically designed for binaries are ``to_audays`` and ``to_sudays``. These unit conversion function convert orbital period ``pb`` to units days, semi-major axis ``semi`` to units of solar units (su) or astronomical units (au). When these are called, the masses of individual binary stars (``m1`` and ``m2``) are converted to solar masses. When ``to_nbody`` is called, all of these values are returned to nbody units.


More support for binary stars is planned for future releases.
