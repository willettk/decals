# Find entries in the NASA-Sloan Atlas (NSA) that have a brick with images in Dark Energy Camera Legacy Survey (DECaLS)

import errno
import os
import random

from astropy.table import Table
from astropy.io import fits

from python.get_catalogs.get_joint_nsa_decals_catalog import create_joint_catalog, get_nsa_catalog, get_decals_bricks, apply_selection_cuts
from python.get_images.download_images_threaded import download_images_multithreaded, get_fits_loc


def download_joint_catalog_images(joint_catalog, dr, nsa_version, fits_dir, jpeg_dir, random_sample=False,
                                  overwrite=False):
    '''
    For galaxies with coverage in the DECaLS bricks, download FITS images from cutout service and make JPGs

    Args:
        joint_catalog (astropy.Table): catalog of galaxies in both nsa and decals including RA, Dec, Petro.
        dr (str): data release to download
        nsa_version (str): version of NSA catalog being used. Defines output catalog filenames.
        random_sample (bool): if True, restrict to a random sample of 101 galaxies only
        overwrite (bool): if True, download FITS and remake JPEG even if identically-named file(s) already exist

    Returns:
        None
    '''

    if random_sample:
        N = 101
        galaxies = random.sample(joint_catalog, N)
    else:
        galaxies = joint_catalog

    joint_catalog = download_images_multithreaded(galaxies, dr, fits_dir, jpeg_dir, overwrite=overwrite)

    return joint_catalog


if __name__ == "__main__":
    # Run all steps to create the NSA-DECaLS-GZ catalog

    data_release = '5'

    catalog_dir = '/data/galaxy_zoo/decals/catalogs'

    fits_dir = '/data/galaxy_zoo/decals/fits/dr{}'.format(data_release)
    jpeg_dir = '/data/galaxy_zoo/decals/jpeg/dr{}'.format(data_release)

    nsa_version = '0_1_2'
    # nsa_version = '1_0_0'
    nsa_catalog_loc = '{}/nsa_v{}.fits'.format(catalog_dir, nsa_version)

    joint_catalog_loc = '{0}/nsa_v{1}_decals_dr{2}.fits'.format(catalog_dir, nsa_version, data_release)

    if data_release == '5' or data_release == '3':
        bricks_filename = 'survey-bricks-dr{}-with-coordinates.fits'.format(data_release)
    elif data_release == '2':
        bricks_filename = 'decals-bricks-dr2.fits'
    elif data_release == '1':
        bricks_filename = 'decals-bricks-dr1.fits'
    else:
        raise ValueError('Data Release "{}" not recognised'.format(data_release))
    bricks_loc = '{}/{}'.format(catalog_dir, bricks_filename)

    nsa = get_nsa_catalog(nsa_catalog_loc)
    bricks = get_decals_bricks(bricks_loc, data_release)

    new_catalog = True
    if new_catalog:
        joint_catalog = create_joint_catalog(nsa, bricks, data_release, nsa_version, run_to=None)  # set None not -1
        joint_catalog.write(joint_catalog_loc, overwrite=True)
        # TODO still need to apply broken PETRO check (small number of cases)
    else:
        joint_catalog = Table(fits.getdata(joint_catalog_loc))

    selected_joint_catalog = apply_selection_cuts(joint_catalog)

    new_images = True
    if new_images:
        joint_catalog_after_download = download_joint_catalog_images(
            selected_joint_catalog,
            data_release,
            nsa_version,
            fits_dir,
            jpeg_dir,
            random_sample=False,
            overwrite=True)
        joint_catalog_after_download.write(joint_catalog_loc, overwrite=True)
