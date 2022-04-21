Utilities
=========

Plotting
--------

There exists several built in plotting tools to help visualize a ``StarCluster``.

.. automodapi:: clustertools.util.plots
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading:

Output
------

There exists several built in tools for creating standardized output files.

.. automodapi:: clustertools.util.output
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading:

It is of course possible to include custom output functions, which for simplicity could be put in the custom folder which is imported. As an example I have included ``custom_outputs.py`` that demonstrates several output functions that I frequently use.

Coordinates
-----------

While coordinate transformations to a ``StarCluster``  are handled via :ref:`Operations <operations>`, users may find the below functions helpful if trying to manually do some coordinate transformations. Note that several of these functions (``cart_to_cyl``,``sky_coords``,``cart_to_sky``) are wrappers around ``galpy`` coordinate transformations.

.. automodapi:: clustertools.util.coordinates
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading:

Recipes
-------

Finally, several functions that are used throughout ``StarCluster`` that users may find helpful are listed below. 

.. automodapi:: clustertools.util.recipes
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading:

Constants
---------

For convenience purposes, a few constants are saved to be used globally throughout ``clustertools``. A function that prints the conversion table between ``kw`` and stellar evolution type, as discussed in :ref:`Cluster <cluster>` is also included here.

.. automodapi:: clustertools.util.constants
        :no-inheritance-diagram:
        :no-main-docstr:
        :no-heading: