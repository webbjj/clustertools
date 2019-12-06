from ..main.cluster import sub_cluster

# Place a mask over the simulated data in order to replicate an observational dataset


def load_mask(cluster, name):

    if name == "m15" or name == "M15":
        # From Emmanulle Dalessandro - reference GO
        cluster.ra_gc = 322.4929716
        cluster.dec_gc = 12.1670539

        cluster.rt = 900.0 / 3600.0
        cluster.rm = 78.0 / 3600.0

        indx = np.ones(cluster.ntot, dtype=bool)

        ocluster = sub_cluster(cluster, indx=indx)

    return ocluster
