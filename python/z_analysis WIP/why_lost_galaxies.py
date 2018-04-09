from astropy.table import Table

import shared_utilities
import pandas as pd

from a_download_decals.get_catalogs import selection_cuts
from a_download_decals.get_images import download_images_threaded
from a_download_decals import download_decals_settings
from a_download_decals.get_catalogs import get_joint_nsa_decals_catalog


def image_download_stats(catalog):
    existing_fits = catalog[catalog['fits_ready']]
    print('existing fits: {}'.format(len(existing_fits)))
    good_fits = catalog[catalog['fits_filled']]
    print('good fits: {}'.format(len(good_fits)))
    good_png = catalog[catalog['png_ready']]
    print('good png: {}'.format(len(good_png)))


lost_loc = '/Data/repos/decals/python/b_to_zooniverse/lost_from_gordon.fits'
lost = Table.read(lost_loc)

nsa_loc = '/Data/repos/decals/nsa_1_0_0_basic_cache.fits'
nsa = Table.read(nsa_loc)

bricks_loc = download_decals_settings.bricks_loc
bricks = Table.read(bricks_loc)
#
#
# print(lost[['objID', 'ra', 'dec', 'sky_separation']])
# print('original: {}'.format(len(lost)))
#
#
# in_nsa, _ = shared_utilities.match_galaxies_to_catalog_table(
#     galaxies=lost,
#     catalog=nsa,
#     galaxy_suffix='',
#     catalog_suffix='_nsa'
# )
#
# print(in_nsa[0])
#
# print('in nsa 1_0_0: {}'.format(len(in_nsa)))
# good_petrotheta = selection_cuts.apply_selection_cuts(in_nsa)
# print('good petrotheta: {}'.format(len(good_petrotheta)))
#
# in_decals_bricks = get_joint_nsa_decals_catalog.create_joint_catalog(in_nsa, bricks, '5')  # dont apply selection cuts
#
# # are they already downloaded? No!
# hopefully_downloaded = in_decals_bricks.copy()
# hopefully_downloaded['fits_loc'] = [
#     download_images_threaded.get_loc(download_decals_settings.fits_dir, galaxy, 'fits') for galaxy in hopefully_downloaded]
# hopefully_downloaded['png_loc'] = [
#     download_images_threaded.get_loc(download_decals_settings.fits_dir, galaxy, 'png') for galaxy in hopefully_downloaded]
# hopefully_downloaded = download_images_threaded.check_images_are_downloaded(hopefully_downloaded, n_processes=1)
# image_download_stats(hopefully_downloaded)
#
# # can we download them? Yes!
# png_dir = 'newly_lost'
# fits_dir = 'newly_lost'
# newly_downloaded = download_images_threaded.download_images_multithreaded(in_decals_bricks, '5', fits_dir, png_dir, overwrite_fits=False, overwrite_png=False)
# image_download_stats(newly_downloaded)


# let's redownload all, without filtering - make temp joint catalog

gordon_galaxies = Table.from_pandas(pd.read_csv('/data/galaxy_zoo/decals/catalogs/gordon_sdss_sample.csv'))
in_nsa, _ = shared_utilities.match_galaxies_to_catalog_table(
    galaxies=gordon_galaxies,
    catalog=nsa,
    galaxy_suffix='',
    catalog_suffix='_nsa'
)
print('not in nsa: {}'.format(len(gordon_galaxies) - len(in_nsa)))
in_decals_bricks = get_joint_nsa_decals_catalog.create_joint_catalog(in_nsa, bricks, '5')  # dont apply selection cuts
png_dir = '/Data/repos/decals/python/z_analysis WIP/newly_lost'
fits_dir = '/Data/repos/decals/python/z_analysis WIP/newly_lost'
joint_catalog = download_images_threaded.download_images_multithreaded(in_decals_bricks, '5', fits_dir, png_dir, overwrite_fits=False, overwrite_png=False)
image_download_stats(joint_catalog)
joint_catalog['iauname'] = list(map(lambda x: str(x), joint_catalog['iauname']))
joint_catalog.write('gordon_joint_catalog_temp.fits', overwrite=True)
