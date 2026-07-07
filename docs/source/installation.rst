Installation
===============

github and pip install
-----------------------

To install ``clustertools``, simply clone the GitHub repository and install via pip:

>>> git clone https://github.com/webbjj/clustertools.git
>>> cd clustertools
>>> pip install .

Please note that if you don't have permission to write files to the default install location (often /usr/local/lib), you will either need to run

>>> pip install --user .

or install into a directory you do have permission to write to and is in your PYTHONPATH:

>>> pip install --prefix='PATH' .

where 'PATH' is that directory.

It has recently become possible to install ``clustertools`` using pip directly from PyPI:

>>> pip install clustertools

however please note this version is not updated as frequently as the GitHub repository. Similarly, if permissions are an issue, you can use:

>>> pip install --user clustertools

or

>>> pip install --prefix='PATH' clustertools

Note: For users looking to take advantage of the plotting features available in clustertools, it may be necessary to install the cm-super, dvipng, texlive, and texlive-extra packages if they aren't installed by default.

Requirements
------------

``clustertools`` requires the following python packages:

* galpy (https://docs.galpy.org/en/v1.7.2/index.html)
* matplotlib
* numpy
* scipy
* astropy
* numba

Optional:

* limepy (https://readthedocs.org/projects/limepy/)
* AMUSE (https://amuse.readthedocs.io/en/latest/)
