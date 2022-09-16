"""A list of commonly used constants
"""

try:
    import amuse.units.units as u
except:
    pass

solar_motion=[-11.1,12.24,7.25] #Schönrich, R., Binney, J., Dehnen, W., 2010, MNRAS, 403, 1829
solar_ro=8.275 #Gravity Collaboration, Abuter, R., Amorim, A., et al. 2020 ,A&A, 647, A59
solar_vo=solar_ro*30.39-solar_motion[1]
# Note vo is set so solar_vo+solar_motion[1] = solar_ro*30.39, where 30.39 is the proper motoin of Sag A* (Reid, M.J. & Brunthaler, A., ApJ, 892, 1)
solar_zo=0.0208 #Bennett, M. & Bovy, J. 2019, MNRAS, 483, 1417

kmperpc = 3.086e13
kmperkpc = 3.086e16
spermyr = 3600.0 * 24.0 * 365.0 * 1000000.0

def kwtypes():
    """ Legend for converting kwtype (from NBODY6) to stellar evolution type

    Parameters
    ----------
    None

    Returns
    -------
    None

    History
    -------
    2019 - Written - Webb (UofT)
    """
    print(
            """\
        *       Stellar evolution types
        *       ***********************
        *
        *       ---------------------------------------------------------------------
        *       0       Low main sequence (M < 0.7).
        *       1       Main sequence.
        *       2       Hertzsprung gap (HG).
        *       3       Red giant.
        *       4       Core Helium burning.
        *       5       First AGB.
        *       6       Second AGB.
        *       7       Helium main sequence.
        *       8       Helium HG.
        *       9       Helium GB.
        *      10       Helium white dwarf.
        *      11       Carbon-Oxygen white dwarf.
        *      12       Oxygen-Neon white dwarf.
        *      13       Neutron star.
        *      14       Black hole.
        *      15       Massless supernova remnant.
        *       ---------------------------------------------------------------------
        *
        *       Binary types
        *       ************
        *
        *       ---------------------------------------------------------------------
        *       0       Standard case.
        *      -1       Chaotic (option 27 = 2).
        *      -2       Continuous circularizing (option 27 = 2).
        *       9       Sequential circularization (option 27 = 1).
        *      10       Circularized.
        *      11       First Roche stage (option 34 = 1/2).
        *      12       End of first Roche stage.
        *      13       Start of second Roche stage.
        *      xx       Further Roche stages.
        *       ---------------------------------------------------------------------
        *



        """
    )

def _get_grav(cluster):
    if cluster.units == "nbody":
        grav = 1.0
    elif cluster.units == "pckms":
        # G has units of pc (km/s)^2 / Msun
        grav = 4.302e-3
    elif cluster.units == "amuse":
        # G has units of pc (km/s)^2 / Msun
        grav = 4.302e-3 | u.pc * u.kms * u.kms / u.MSun
    elif cluster.units == "pcmyr":
        # G has units of pc (pc/myr)^2 / Msun
        grav = (4.302e-3)*(1.022712165045695**2.0)
    elif cluster.units == "kpckms":
        # G has units of kpc (km/s)^2 / Msun
        grav = 4.302e-6
    elif cluster.units == "kpcgyr":
        # G has units of pc (pc/myr)^2 / Msun
        grav = (4.302e-6)*(1.022712165045695**2.0)
    else:
        grav = 1.0

    return grav
