import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy import units as u
import pandas as pd

import matplotlib.pyplot as plt

from shared_utilities import match_galaxies_to_catalog


def get_expert_catalog_joined_with_decals(joint_catalog, expert_catalog, plot=False):

    output_catalog, _ = match_galaxies_to_catalog(
        galaxies=expert_catalog,
        catalog=joint_catalog,
        matching_radius=5 * u.arcsec,
        galaxy_suffix='_expert')

    output_catalog['has_bar'] = output_catalog['bar'] != 0
    output_catalog['has_ring'] = output_catalog['ring'] != 0

    output_catalog['bar_types'] = list(map(decode_bar_ints, output_catalog['bar']))
    output_catalog['ring_types'] = list(map(decode_ring_ints, output_catalog['ring']))

    if plot:
        save_output_catalog_statistics(output_catalog)

    return output_catalog


def get_expert_catalog(expert_catalog_loc):

    catalog = Table(fits.getdata(expert_catalog_loc))

    # convert column names to be consistent with nsa
    for column_name in catalog.colnames:
        catalog.rename_column(column_name, column_name.lower())
    catalog.rename_column('sdss', 'iauname')

    return catalog


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
    output_catalog = output_catalog[['tt', 'has_bar', 'has_ring']].to_pandas()

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
