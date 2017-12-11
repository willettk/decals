from get_catalogs.get_joint_nsa_decals_catalog import get_nsa_catalog
import numpy as np
import pandas as pd


def apply_selection_cuts(input_catalog):
    """
    Select only galaxies with PETROTHETA > 3 and not within 1e-3 of default value (i.e. bad measurement)

    Args:
        catalog (astropy.Table): Galaxy catalog including NSA information

    Returns:
        (astropy.Table) catalog of galaxies matching selection criteria above

    """

    petrotheta_above_3 = input_catalog['petrotheta'] > 3
    # redshift_below_p05 = input_catalog['z'] < 0.05

    '''
    NSA catalog’s PETROTHETA calculation sometimes fails to a ‘default’ value that is related to the annulus used to
    measure the Petrosian radius.
    The relation is (found via fitting the observed ‘snap to’ values to the annulus values):
    PETROTHETA_snap_to = 0.997 * PROFTHETA ** 0.998
    Any galaxies with PETROTHETA within 1e-3 of the snap_to value likely has the wrong size.
    '''
    # TODO discuss with Coleman: PROFTHETA is 15 values, which value (or all)?
    #
    # print(input_catalog['PROFTHETA'])
    # proftheta = np.array(input_catalog['PROFTHETA'])
    # snap_values = 0.997 * np.power(input_catalog['PROFTHETA'], 0.998)  # ** syntax behaves unexpectedly with Column
    # snap_tolerance = 1e-3
    # snap_lower_limit = snap_values - snap_tolerance
    # snap_upper_limit = snap_values + snap_tolerance
    #
    # above_snap_lower_limit = input_catalog['PETROTHETA'] > snap_lower_limit
    # below_snap_upper_limit = input_catalog['PETROTHETA'] < snap_upper_limit
    # within_snap_window = above_snap_lower_limit & below_snap_upper_limit
    #
    # selected_catalog = input_catalog[petrotheta_above_3 & ~ within_snap_window]
    # print(len(input_catalog), len(selected_catalog))
    # return selected_catalog

    return input_catalog[petrotheta_above_3]
    # return input_catalog[petrotheta_above_3 & redshift_below_p05]


if __name__ == '__main__':

    nsa_version = '1_0_0'
    nsa_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nsa_v1_0_0.fits'
    nsa = get_nsa_catalog(nsa_catalog_loc, nsa_version)
    import matplotlib.pyplot as plt
    nsa_small = nsa[nsa['petrotheta'] < 50]
    plt.hist(nsa_small['petrotheta'], bins=100)
    plt.ylabel('Count')
    plt.xlabel('Petrotheta for NSA galaxies according to catalog')
    plt.savefig('petrotheta_hist.png')

    plt.clf()
    nsa_wrong = nsa[(nsa['petrotheta'] > 26) & (nsa['petrotheta'] < 29)]
    plt.hist(nsa_wrong['petrotheta'], bins=40)
    plt.ylabel('Count')
    plt.xlabel('Petrotheta for NSA galaxies according to catalog')
    plt.savefig('petrotheta_hist_28.png')

    print(pd.Series(nsa_wrong['petrotheta']).value_counts())

    nsa_bad = nsa[nsa['petrotheta'] == 27.653702]