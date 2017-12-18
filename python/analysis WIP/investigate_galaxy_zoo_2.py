import pandas as pd
from astropy.io import fits
from astropy.table import Table

from astropy.coordinates import SkyCoord
from astropy import units as u

import datashader as ds
import datashader.transfer_functions as tf
from datashader.utils import export_image

import matplotlib.pyplot as plt

pd.options.display.max_rows = 200
pd.options.display.max_columns = 100


def match_galaxies_to_catalog(galaxies, catalog, matching_radius=10 * u.arcsec):
    # http://docs.astropy.org/en/stable/coordinates/matchsep.html

    galaxies_coord = SkyCoord(ra=galaxies['ra'] * u.degree, dec=galaxies['dec'] * u.degree)
    catalog_coord = SkyCoord(ra=catalog['ra'] * u.degree, dec=catalog['dec'] * u.degree)
    best_match_catalog_index, sky_separation, _ = galaxies_coord.match_to_catalog_sky(catalog_coord)

    galaxies['best_match'] = best_match_catalog_index
    galaxies['sky_separation'] = sky_separation.to(u.arcsec).value
    matched_galaxies = galaxies[galaxies['sky_separation'] < matching_radius.value]

    catalog['best_match'] = catalog.index.values

    matched_catalog = pd.merge(matched_galaxies, catalog, on='best_match', how='inner', suffixes=['_subject', ''])

    return matched_catalog


def plot_catalog_overlap(catalog_a, catalog_b, legend):

    a_coords = catalog_a[['ra', 'dec']]
    a_coords['catalog'] = legend[0]
    b_coords = catalog_b[['ra', 'dec']]
    b_coords['catalog'] = legend[1]

    df_to_plot = pd.concat([a_coords, b_coords])
    df_to_plot['catalog'] = df_to_plot['catalog'].astype('category')

    canvas = ds.Canvas(plot_width=300, plot_height=300)
    aggc = canvas.points(df_to_plot, 'ra', 'dec', ds.count_cat('catalog'))
    img = tf.shade(aggc)
    export_image(img, 'catalog_overlap')



catalog_dir = '/data/galaxy_zoo/decals/catalogs'

gz_catalog = pd.read_csv('{}/nsa_all_raw_gz_counts_10.0_arcsec.csv'.format(catalog_dir))
print('gz nsa galaxies: {}'.format(len(gz_catalog)))

joint_catalog = Table(fits.getdata('{}/nsa_v1_0_0_decals_dr5.fits'.format(catalog_dir)))
joint_catalog = joint_catalog[['nsa_id', 'ra', 'dec', 'petrotheta']].to_pandas()
print('decals nsa galaxies: {}'.format(len(joint_catalog)))

gz_and_decals = pd.merge(gz_catalog, joint_catalog, on='nsa_id', how='inner')
print('gz and decals: {}'.format(len(gz_and_decals)))

old_galaxies = match_galaxies_to_catalog(joint_catalog, gz_and_decals)
new_galaxies = joint_catalog[~joint_catalog['iauname'].isin(set(old_galaxies['iauname']))]

print(len(new_galaxies))
plot_catalog_overlap(old_galaxies, new_galaxies, ['old', 'new'])

new_galaxies['petrotheta'].dist()



    # interesting_cols = ['data_release', 'nsa_id', 'ra', 'dec']
# decals_dr1_dr2 = pd.read_csv(
#         '/data/galaxy_zoo/decals/catalogs/galaxy_zoo_decals_with_nsa_v1_0_0.csv',
#         nrows=None,
#         usecols=interesting_cols)
# print('decals_dr1_dr2 nsa galaxies: {}'.format(len(decals_dr1_dr2)))
#
# decals_all = pd.merge(decals_dr1_dr2, joint_catalog, on='nsa_id', how='outer')
# print('decals all: {}'.format(len(decals_all)))
# decals_repeated = pd.merge(decals_dr1_dr2, joint_catalog, on='nsa_id', how='inner')
# print('decals repeated: {}'.format(len(decals_repeated)))
#
# gz_or_previous_decals = pd.merge(gz_catalog, decals_dr1_dr2, on='nsa_id', how='outer')
# print('gz or previous decals: {}'.format(len(gz_or_previous_decals)))
#
# gz_and_previous_decals = pd.merge(gz_catalog, decals_dr1_dr2, on='nsa_id', how='inner')
# print('gz and previous decals: {}'.format(len(gz_and_previous_decals)))
#
# all_galaxies = pd.merge(gz_or_previous_decals, decals_all, on='nsa_id', how='outer')
# print('all galaxies: {}'.format(len(all_galaxies)))
#
#
# new_dr5_galaxies = joint_catalog[~joint_catalog['nsa_id'].isin(gz_and_previous_decals['nsa_id'])]
# print(len(new_dr5_galaxies))
#
# new_small_galaxies = new_dr5_galaxies[new_dr5_galaxies['petrotheta'] < 50]
# new_small_galaxies['petrotheta'].hist(bins=30)
# plt.xlabel('r petrosian radius')
# plt.ylabel('count')
# plt.tight_layout()
# plt.savefig('new_petrotheta.png')



# new_dr5_galaxies['']

# https://stackoverflow.com/questions/23284409/how-to-subtract-rows-of-one-pandas-data-frame-from-another
# joint_catalog['name'] = 'joint_catalog'
# # gz_and_previous_decals['name'] = 'gz_and_previous_decals'
# new_dr5_galaxies = pd.merge(joint_catalog, gz_and_previous_decals, on='nsa_id', how='left')