from astroquery.vizier import Vizier

import numpy as np
from astropy.table import Table, join
from astropy.io import fits

import matplotlib.pyplot as plt


def get_expert_catalog(row_limit=-1):
    catalog_list = Vizier.find_catalogs('J/ApJS/186/427')
    Vizier.ROW_LIMIT = row_limit
    catalog = Vizier.get_catalogs(catalog_list.keys())[0]
    # convert column names to be consistent with nsa
    for column_name in catalog.colnames:
        catalog.rename_column(column_name, column_name.lower())
    catalog.rename_column('sdss', 'iauname')

    catalog['iauname_1dp'] = 'J000000.00-000000.0'
    for galaxy_n, galaxy in enumerate(catalog):
        catalog['iauname_1dp'][galaxy_n] = round_iauname_to_1dp(galaxy['iauname'])

    return catalog


def round_iauname_to_1dp(iauname):

    # simple cut
    iauname_1dp = iauname[:-1]

    # rounding
    # final_digits = float(iauname[-4:])
    # rounded_digits = np.around(final_digits, decimals=1)
    #
    # remaining_string = iauname[:-4]
    #
    # iauname_1dp = remaining_string + str(rounded_digits)
    #
    # print('\n')
    # print(iauname)
    # print(len(iauname))
    # print(final_digits)
    # print(remaining_string)
    # print(iauname_1dp)
    # print(len(iauname_1dp))

    return iauname_1dp


if __name__ == '__main__':
    row_limit = -1
    catalog_dir = '/data/galaxy_zoo/decals/catalogs'
    joint_columns = [
        'nsa_id',
        'iauname',
        'ra',
        'dec']

    expert_catalog = get_expert_catalog(row_limit)
    print(expert_catalog.colnames)

    joint_catalog = Table(fits.getdata('{}/nsa_v1_0_0_decals_dr5.fits'.format(catalog_dir)))[joint_columns]
    joint_catalog['iauname_1dp'] = joint_catalog['iauname']

    plt.scatter(joint_catalog['ra'], joint_catalog['dec'])
    plt.scatter(expert_catalog['_ra'], expert_catalog['_de'], c='r')
    plt.legend(['joint', 'expert'])
    plt.savefig('locations.png')

    # print(expert_catalog['iauname_1dp'])
    # print(joint_catalog['iauname_1dp'])

    output_catalog = join(joint_catalog, expert_catalog, keys='iauname_1dp', join_type='inner', table_names=['joint', 'expert'])
    # output_catalog = join(joint_catalog, expert_catalog, keys='iauname_1dp', join_type='outer', table_names=['joint', 'expert'])
    print(output_catalog.colnames)
    output_catalog.sort('iauname_1dp')
    print(len(output_catalog))

    # target_galaxy = 'J113057.91-010851.1' in viz, matches
    # target_galaxy = 'J110122.00-010824.9'  # in viz, doesn't match'
    # galaxy_loc = np.argmax(expert_catalog['iauname'] == target_galaxy)
    # print(galaxy_loc)
    #
    # print(output_catalog[:5])
    # galaxy_loc = np.argmax(output_catalog['iauname_joint'] == target_galaxy)
    # print(galaxy_loc)
    # print(output_catalog[galaxy_loc-2:galaxy_loc+2])

    # output_catalog.write(format='csv')
    print(output_catalog)
    print(len(output_catalog))