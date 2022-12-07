Installation
===============

github and pip install
------------

To install ``clustertools``, simply clone the GitHub repository and install via setup tools:

>>> git clone https://github.com/webbjj/clustertools.git
>>> cd clustertools
>>> python setup.py install

Please note that if you don't have permission to write files to the default install location (often /usr/local/lib), you will either need to run

>>> sudo python setup.py install

or

>>> python setup.py install --prefix='PATH'

where 'PATH' is a directory that you do have permission to write in and is in your PYTHONPATH.

It has recently become possible to install ``clustertools`` using pip:

>>> pip install clustertools

however please note this version is not updated as frequently as the GitHub repository. Similarly, if permissions are an issue, you can use:

>>> pip install --user clustertools

or

>>> pip install clustertools --install-option="--prefix=PATH"

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
