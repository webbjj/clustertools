# clustertools

clustertools is a Python package for analysing star cluster simulations. The package is built around the StarCluster class, which will store all the necessary information from a given star cluster simulation to be used for anaylsis. All functions within clustertools are then designed to act on a StarCluster. clustertools can be used for unit and coordinate transformations, the calculation of key structural and kinematic parameters, analysis of the cluster’s orbit and tidal tails (with the help of galpy , and measuring common cluster properties like its mass function, density profile, and velocity dispersion profile (among others). While originally designed with star clusters in mind, clustertools can be used to study other types of N-body systems, including stellar streams and dark matter sub-halos.

The package contains functions for loading data from commonly used N-body codes, generic snapshots, and codes for generating initial conditions.

clustertools is developed on Github. Please go to https://github.com/webbjj/clustertools to report issues or contribute to the code.

Documentation for clustertools can be found at https://clustertools.readthedocs.io/en/latest/

# Installation

clustertools can be installed either directly from the GitHub repository or via pip

To install clustertools from GitHub, simply clone the repository and install via setup tools:

git clone https://github.com/webbjj/clustertools.git
cd clustertools
python setup.py install
Please note that if you don’t have permission to write files to the default install location (often /usr/local/lib), you will either need to run

sudo python setup.py install
or

python setup.py install --prefix='PATH'
where ‘PATH’ is a directory that you do have permission to write in and is in your PYTHONPATH.

It is also possible to install clustertools using pip:

pip install clustertools
however please note this version is not updated as frequently as the GitHub repository. Similarly, if permissions are an issue, you can use:

pip install --user clustertools
or

pip install clustertools --install-option="--prefix=PATH"

Note: For users looking to take advantage of the plotting features available in clustertools, it may be necessary to install the cm-super and dvipng packages if they aren't installed by default.

# Requirements
clustertools requires the following python packages:

galpy (https://docs.galpy.org/en/v1.7.2/index.html)
matplotlib
numpy
scipy
astropy
numba

Optional:
limepy (https://readthedocs.org/projects/limepy/)
AMUSE (https://amuse.readthedocs.io/en/latest/)