# Welcome to nbodypy, a python packaged aimed at reading in and analysing the results of star cluster simulations. 

# The structure of the code is as follows:
# Read in data from a star cluster simulation (various routines are already written in get_nbody.py, and a custom read in function is easily built
# Use data from simulation to define a StarCluster instance (cluster.py)
# -- Instance get can be defined with or without some initial information
# -- Stellar properties (with or without stellar evolution) can be included when the 
#    instance is first defined or added later.
# Once a StarCluster instance has been defined, various operations and functions are in place to manipulate the data and/or calculation key properties (units.py, functions.py)
