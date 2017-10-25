# Find entries in the NASA-Sloan Atlas (NSA) that have a brick with images in Dark Energy Camera Legacy Survey (DECaLS)

import errno
import os
import random

from astropy.table import Table
from astropy.io import fits

from python.get_catalogs.get_joint_nsa_decals_catalog import create_joint_catalog, get_nsa_catalog, get_decals_bricks
from python.get_images.download_images_threaded import download_images_multithreaded


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

    timed_out, good_images = download_images_multithreaded(joint_catalog, dr, fits_dir, jpeg_dir, overwrite=overwrite)

    # Write catalogs of time-outs and good images to file
    galaxies[timed_out].write('../fits/nsa_v{0}_decals_dr{1}_timedout.fits'.format(nsa_version, dr), overwrite=True)
    galaxies[good_images].write('../fits/nsa_v{0}_decals_dr{1}_goodimgs.fits'.format(nsa_version, dr), overwrite=True)


if __name__ == "__main__":
    # Run all steps to create the NSA-DECaLS-GZ catalog

    catalog_dir = '/data/galaxy_zoo/decals/catalogs'
    # TODO move to decals folder once development complete
    fits_dir = '../fits/nsa'
    jpeg_dir = '../jpeg/dr2'

    nsa_version = '0_1_2'
    nsa_catalog_loc = '{}/nsa_v{}.fits'.format(catalog_dir, nsa_version)

    # TODO should document which files are being used and then rename
    data_release = '3'
    if data_release == '3':
        # http: // legacysurvey.org / dr3 / files /
        bricks_filename = 'survey-bricks-dr3-with-coordinates.fits'
    elif data_release == '2':
        # http: // legacysurvey.org / dr2 / files /
        bricks_filename = 'decals-bricks-dr2.fits'
    bricks_loc = '{}/{}'.format(catalog_dir, bricks_filename)

    nsa = get_nsa_catalog(nsa_catalog_loc)
    bricks = get_decals_bricks(bricks_loc, data_release)

    new_catalog = True
    if new_catalog:
        joint_catalog = create_joint_catalog(nsa, bricks, data_release, nsa_version, run_to=100)
    else:
        # TODO still need to apply cuts - happened externally via Topcat
        joint_catalog = Table(fits.getdata('../fits/nsa_v{0}_decals_dr{1}.fits'.format(nsa_version, data_release), 1))

    download_joint_catalog_images(joint_catalog,
                                  data_release,
                                  nsa_version,
                                  fits_dir,
                                  jpeg_dir,
                                  random_sample=False,
                                  overwrite=True)
