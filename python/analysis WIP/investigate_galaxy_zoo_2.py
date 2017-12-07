import pandas as pd
from astropy.io import fits
from astropy.table import Table

import matplotlib.pyplot as plt

catalog_dir = '/data/galaxy_zoo/decals/catalogs'

gz_catalog = pd.read_csv('{}/nsa_all_raw_gz_counts_10.0_arcsec.csv'.format(catalog_dir))
print('gz nsa galaxies: {}'.format(len(gz_catalog)))

joint_catalog = Table(fits.getdata('{}/nsa_v1_0_0_decals_dr5.fits'.format(catalog_dir)))
joint_catalog = joint_catalog[['nsa_id', 'ra', 'dec', 'petrotheta']].to_pandas()
print('decals nsa galaxies: {}'.format(len(joint_catalog)))

gz_and_decals = pd.merge(gz_catalog, joint_catalog, on='nsa_id', how='inner')
print('gz and decals: {}'.format(len(gz_and_decals)))

interesting_cols = ['data_release', 'nsa_id', 'ra', 'dec']
decals_dr1_dr2 = pd.read_csv(
        '/data/galaxy_zoo/decals/catalogs/galaxy_zoo_decals_with_nsa_v1_0_0.csv',
        nrows=None,
        usecols=interesting_cols)
print('decals_dr1_dr2 nsa galaxies: {}'.format(len(decals_dr1_dr2)))

decals_all = pd.merge(decals_dr1_dr2, joint_catalog, on='nsa_id', how='outer')
print('decals all: {}'.format(len(decals_all)))
decals_repeated = pd.merge(decals_dr1_dr2, joint_catalog, on='nsa_id', how='inner')
print('decals repeated: {}'.format(len(decals_repeated)))

gz_or_previous_decals = pd.merge(gz_catalog, decals_dr1_dr2, on='nsa_id', how='outer')
print('gz or previous decals: {}'.format(len(gz_or_previous_decals)))

gz_and_previous_decals = pd.merge(gz_catalog, decals_dr1_dr2, on='nsa_id', how='inner')
print('gz and previous decals: {}'.format(len(gz_and_previous_decals)))

all_galaxies = pd.merge(gz_or_previous_decals, decals_all, on='nsa_id', how='outer')
print('all galaxies: {}'.format(len(all_galaxies)))


new_dr5_galaxies = joint_catalog[~joint_catalog['nsa_id'].isin(gz_and_previous_decals['nsa_id'])]
print(len(new_dr5_galaxies))

new_small_galaxies = new_dr5_galaxies[new_dr5_galaxies['petrotheta'] < 50]
new_small_galaxies['petrotheta'].hist(bins=30)
plt.xlabel('r petrosian radius')
plt.ylabel('count')
plt.tight_layout()
plt.savefig('new_petrotheta.png')



# new_dr5_galaxies['']

# https://stackoverflow.com/questions/23284409/how-to-subtract-rows-of-one-pandas-data-frame-from-another
# joint_catalog['name'] = 'joint_catalog'
# # gz_and_previous_decals['name'] = 'gz_and_previous_decals'
# new_dr5_galaxies = pd.merge(joint_catalog, gz_and_previous_decals, on='nsa_id', how='left')