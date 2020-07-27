Binaries
===============

At the moment, ``clustertools`` has minimal support for analysing binary stars in a ``StarCluster``. The only function designed specifically designed for binaries is ``convert_binary_units``, where one can convert orbital period ``pb`` to units days,years, or Nbody units, semi-major axis ``semi`` to units of pc, solar units (su), astronomical units (au), or Nbody units, and mass ``mass`` to units of solar masses or Nbody units. An example of this step would be:

>>> cluster.convert_binary_units(['pb','semi','mass'],from_units=['nbody','nbody','nbody')],to_units=['years','su','Msun']

Of course, knowing which stars are binary stars in your datasets allows for them to be extracted using ``sub_cluster`` and analysed individually. More support for binary stars is planned for future releases.
