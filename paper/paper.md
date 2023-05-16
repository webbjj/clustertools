---
title: 'clustertools: A Python Package for Analyzing Star Cluster Simulations'
tags:
  - Python
  - astronomy
  - star clusters
  - N-body simulations
authors:
  - name: Jeremy J. Webb # note this makes a footnote saying 'Co-first author'
    orcid: 0000-0003-3613-0854
    affiliation: 1 # (Multiple affiliations must be quoted)
affiliations:
 - name: David A. Dunlap Department of Astronomy and Astrophysics, University of Toronto, 50 St. George Street, Toronto, ON, M5S 3H4, Canada
   index: 1
date: 6 May 2022
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

\texttt{clustertools} is a Python package for analyzing star cluster simulations. The package is built around the \texttt{StarCluster} class, which stores all data read in from the snapshot of a given model star cluster. The package contains functions for loading data from commonly used N-body codes, generic snapshots, and software for generating initial conditions. All operations and functions within \texttt{clustertools} are then designed to act on a \texttt{StarCluster}. \texttt{clustertools} can be used for unit and coordinate transformations, the calculation of key structural and kinematic parameters, analysis of the cluster's orbit and tidal tails, and measuring common cluster properties like its mass function, density profile, and velocity dispersion profile (among others). While originally designed with star clusters in mind, \texttt{clustertools} can be used to study other types of $N$-body systems, including stellar streams and dark matter sub-halos.

# Statement of need

Stars do not form alone, but in clustered environments that in some cases can remain gravitationally bound as a star cluster for billions of years. The details of how exactly star clusters form, either at high-redshifts or at present day, remain unknown. After formation, the subsequent evolution of these stellar systems has been shown to be strongly linked to that of their host galaxy. Hence star clusters are often used as tools for studying star formation, galaxy evolution, and galaxy structure. The comparison of simulated star clusters to observations offers the ability to explore what formation conditions reproduce the properties of observed star cluster populations and how the evolution and structure of a galaxy affects star clusters.  

A large number of $N$-body codes and software packages exist to generate star cluster models and simulate their long-term evolution in different environments [e.g. AMUSE, @amuse; NBODY6, @nbody6; NBODY6++, @nbody6pp; NEMO, @nemo; PETAR, @petar]. Additional software exists for generating model star clusters directly from a known distribution function [e.g. GALPY, @galpy; LIMEPY, @limepy; MCLUSTER, @mcluster]. The output from these codes and packages can differ in units, coordinate system, and format. The subsequent analysis of star cluster simulations and models can then be quite inhomogeneous between studies, as various methods exist for determining things like a cluster's centre, relaxation time, or tidal radius. \texttt{clustertools} allows for data from a star cluster simulation or model to be loaded and analyzed homogeneously by first loading a \texttt{StarCluster} instance. Within a \texttt{StarCluster}, stellar masses, positions, velocities, identification numbers, and all relevant meta-data from the source is saved. \texttt{clustertools} then contains an array of operations for converting units and coordinate systems, functions for calculating key cluster parameters, and functions for measuring how certain parameters vary with clustercentric distance. In some cases, multiple different methods are available for calculating a specific parameter. The software is particularly useful to students or users new to star cluster studies, as functions and profile measurements have a plotting feature that helps illustrate how certain parameters are measured when possible. This approach also ensures that the analysis of any star cluster is done in a homogenous way with open-source code, regardless of the simulation or model from which that data was produced. 


# Acknowledgements

JJW would like to thank Jo Bovy and Nathaniel Starkman for helpful discussions and contributions to \texttt{clustertools}


# References
