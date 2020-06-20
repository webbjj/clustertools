"""
NAME:

constants

PURPOSE:

A list of commonly used constants

INPUT:

None

OUTPUT:

kmperpc - convert km to pc

kmperkpc - convert km to kpc

spermyr - convert seconds to myr

HISTORY:

2018 - Written - Webb (UofT)
"""

kmperpc = 3.086e13
kmperkpc = 3.086e16
spermyr = 3600.0 * 24.0 * 365.0 * 1000000.0

def kwtypes():
    """
    NAME:

       kwtypes

    PURPOSE:

       Print legend for converting kwtype (from NBODY6) to stellar evolution type

    INPUT:

       None

    OUTPUT:

       None

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    print(
        dedent(
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
    )
