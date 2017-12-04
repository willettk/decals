from astroquery.vizier import Vizier

import numpy as np
from astropy.table import Table, join
from astropy.io import fits
import pandas as pd

import matplotlib.pyplot as plt
import datashader as ds
import datashader.transfer_functions as tf
from datashader.utils import export_image
import seaborn as sns


def get_expert_catalog_joined_with_decals():
    catalog_dir = '/data/galaxy_zoo/decals/catalogs'
    joint_columns = [
        'nsa_id',
        'iauname',
        'ra',
        'dec']

    expert_catalog = get_expert_catalog()

    joint_catalog = Table(fits.getdata('{}/nsa_v1_0_0_decals_dr5.fits'.format(catalog_dir)))[joint_columns]
    joint_catalog['iauname_1dp'] = joint_catalog['iauname']

    output_catalog = join(joint_catalog, expert_catalog, keys='iauname_1dp', join_type='inner', table_names=['joint', 'expert'])

    output_catalog['has_bar'] = output_catalog['bar'] != 0
    output_catalog['has_ring'] = output_catalog['ring'] != 0

    output_catalog['bar_types'] = map(decode_bar_ints, output_catalog['bar'])
    output_catalog['ring_types'] = map(decode_ring_ints, output_catalog['ring'])

    output_catalog['iauname'] = output_catalog['iauname_joint']
    del output_catalog['iauname_joint']

    save_output_catalog_statistics(output_catalog)

    return output_catalog


def get_expert_catalog():

    catalog = Table(fits.getdata('/data/galaxy_zoo/decals/catalogs/nair_sdss_catalog.fit'))

    # convert column names to be consistent with nsa
    for column_name in catalog.colnames:
        catalog.rename_column(column_name, column_name.lower())
    catalog.rename_column('sdss', 'iauname')

    # nsa/decals catalog uses 1dp precision in galaxy names - round expert catalog to match this
    catalog['iauname_1dp'] = 'J000000.00-000000.0'
    for galaxy_n, galaxy in enumerate(catalog):
        catalog['iauname_1dp'][galaxy_n] = round_iauname_to_1dp(galaxy['iauname'])

    return catalog


def plot_joint_and_expert_overlap(joint_catalog, expert_catalog):

    joint_to_plot = joint_catalog[['ra', 'dec']].to_pandas()
    joint_to_plot['catalog'] = 'joint'
    expert_to_plot = expert_catalog[['_ra', '_de']].to_pandas()
    expert_to_plot.rename(columns={'_ra': 'ra', '_de': 'dec'}, inplace=True)
    expert_to_plot['catalog'] = 'expert'

    df_to_plot = pd.concat([joint_to_plot] + [expert_to_plot for n in range(int(len(joint_to_plot)/len(expert_to_plot)))])
    df_to_plot['catalog'] = df_to_plot['catalog'].astype('category')

    canvas = ds.Canvas(plot_width=300, plot_height=300)
    aggc = canvas.points(df_to_plot, 'ra', 'dec', ds.count_cat('catalog'))
    img = tf.shade(aggc)
    export_image(img, 'locations_datashader_overlaid')


def decode_bar_ints(encoded_int):
    bar_int_labels = decode_binary_mask(encoded_int)
    return [bar_int_label_to_str(label) for label in bar_int_labels]


def decode_ring_ints(encoded_int):
    ring_int_labels = decode_binary_mask(encoded_int)
    return [ring_int_label_to_str(label) for label in ring_int_labels]


def decode_binary_mask(encoded_int):
    # convert to binary base, represent as string
    binary_int_string = bin(encoded_int)[2:]
    # convert to array
    binary_int_array = np.array([int(x) for x in binary_int_string])
    # reading right to left, count off the index of each 1
    # flip and then get left-to-right index of nonzero elements
    indices_of_nonzero = np.nonzero(np.flip(binary_int_array, axis=0))[0]
    return list(indices_of_nonzero)


def bar_int_label_to_str(bar_int_label):
    mapping = {
        1: 'strong',
        2: 'intermediate',
        3: 'weak',
        4: 'ansae',
        5: 'peanut',
        6: 'nuclear',
        7: 'unsure'
    }
    return mapping[bar_int_label]


def ring_int_label_to_str(ring_int_label):
    mapping = {
        1: 'nuclear',
        2: 'inner',
        3: 'outer'
    }
    return mapping[ring_int_label]


def save_output_catalog_statistics(output_catalog):
    output_catalog = output_catalog.to_pandas()

    output_catalog['tt'].value_counts().sort_index().plot(kind='bar')
    plt.ylabel('Count')
    plt.xlabel('T-Type of Nair galaxy in DECALS DR5')
    plt.tight_layout()
    plt.savefig('ttype_counts.png')
    plt.clf()

    output_catalog['has_bar'].value_counts().sort_index().plot(kind='bar')
    plt.ylabel('Count')
    plt.xlabel('Bar of Nair galaxy in DECALS DR5')
    plt.tight_layout()
    plt.savefig('bar_counts.png')
    plt.clf()

    output_catalog['has_ring'].value_counts().sort_index().plot(kind='bar')
    plt.ylabel('Count')
    plt.xlabel('Ring of Nair galaxy in DECALS DR5')
    plt.tight_layout()
    plt.savefig('ring_counts.png')
    plt.clf()


def column_value_counts_as_dataframe(data, column):
    counts = pd.Series(data[column]).value_counts()
    return pd.DataFrame({column: counts.index, 'count': counts.values})


def round_iauname_to_1dp(iauname, debug=False):

    # simple cut
    iauname_1dp = iauname[:-1]

    # empirically, simple cut recovers more galaxies than rounding - let's not round for now
    # rounding
    # final_digits = float(iauname[-4:])
    # rounded_digits = np.around(final_digits, decimals=1)
    #
    # remaining_string = iauname[:-4]
    #
    # iauname_1dp = remaining_string + str(rounded_digits)
    #
    # if debug:
    #     print('\n')
    #     print(iauname)
    #     print(len(iauname))
    #     print(final_digits)
    #     print(remaining_string)
    #     print(iauname_1dp)
    #     print(len(iauname_1dp))

    return iauname_1dp


if __name__ == '__main__':
    get_expert_catalog_joined_with_decals()
