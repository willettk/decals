# Find entries in the NASA-Sloan Atlas (NSA) that have a brick with images in Dark Energy Camera Legacy Survey (DECaLS)

import errno
import os
import random

from astropy.table import Table
from astropy.io import fits
import pandas as pd

from get_catalogs.get_joint_nsa_decals_catalog import create_joint_catalog, get_nsa_catalog, get_decals_bricks, apply_selection_cuts
from get_images.download_images_threaded import download_images_multithreaded, get_fits_loc
from get_catalogs.previous_subjects import get_previous_subjects_with_nsa
from get_catalogs.find_new_subjects import find_new_catalog_images


class Settings():
    # Define parameters, file locations, etc.

    def get_nsa_catalog_loc(self):
        return '{}/nsa_v{}.fits'.format(self.catalog_dir, self.nsa_version)

    def get_joint_catalog_loc(self):
        return '{0}/nsa_v{1}_decals_dr{2}.fits'.format(
            self.catalog_dir, self.nsa_version, self.data_release)

    def get_upload_catalog_loc(self):
        return '{}/dr{}_nsa{}_to_upload.csv'.format(self.catalog_dir, self.data_release, self.nsa_version)

    def get_bricks_loc(self):
        data_release = self.data_release
        if data_release == '5' or data_release == '3':
            bricks_filename = 'survey-bricks-dr{}-with-coordinates.fits'.format(data_release)
        elif data_release == '2':
            bricks_filename = 'decals-bricks-dr2.fits'
        elif data_release == '1':
            bricks_filename = 'decals-bricks-dr1.fits'
        else:
            raise ValueError('Data Release "{}" not recognised'.format(data_release))
        return '{}/{}'.format(self.catalog_dir, bricks_filename)

    def derive_file_paths(self):
        self.nsa_catalog_loc = self.get_nsa_catalog_loc()
        self.joint_catalog_loc = self.get_joint_catalog_loc()
        self.bricks_loc = self.get_bricks_loc()
        self.upload_catalog_loc = self.get_upload_catalog_loc()

    def __init__(self,
                 data_release='3',
                 nsa_version='1_0_0',
                 catalog_dir='/data/galaxy_zoo/decals/catalogs',
                 fits_dir='/data/galaxy_zoo/decals/fits/unknown_dr',
                 jpeg_dir='/data/galaxy_zoo/decals/jpeg/unknown_dr',
                 subject_loc='/data/galaxy_zoo/decals/subjects/decals_dr1_and_dr2.csv'):

        self.data_release = data_release
        self.nsa_version = nsa_version
        self.catalog_dir = catalog_dir
        self.fits_dir = fits_dir
        self.jpeg_dir = jpeg_dir
        self.subject_loc = subject_loc

        self.derive_file_paths()


def get_decals(nsa, bricks, previous_subjects, s):

    nsa_after_cuts = apply_selection_cuts(nsa)

    if s.new_catalog:
        joint_catalog = create_joint_catalog(nsa_after_cuts, bricks, s.data_release, s.nsa_version, run_to=s.run_to)
        joint_catalog.write(s.joint_catalog_loc, overwrite=True)
        # TODO still need to apply broken PETRO check (small number of cases)
    else:
        joint_catalog = Table(fits.getdata(s.joint_catalog_loc))

    if s.new_images:
        joint_catalog = download_images_multithreaded(
            joint_catalog[:s.run_to],
            s.data_release,
            s.fits_dir,
            s.jpeg_dir,
            overwrite_fits=s.overwrite_fits,
            overwrite_jpeg=s.overwrite_jpeg)
        joint_catalog.write(s.joint_catalog_loc, overwrite=True)

    # add nsa info to previous gz subjects downloaded from data dump
    previous_galaxy_zoo_with_nsa = get_previous_subjects_with_nsa(previous_subjects, nsa)

    # find_new_catalog_images expects two Table catalogs
    # galaxy zoo catalog should include where the FITS would be
    previous_galaxy_zoo_with_nsa['fits_loc'] = [get_fits_loc(s.catalog_dir, galaxy) for _, galaxy in previous_galaxy_zoo_with_nsa.iterrows()]
    # TODO a bit messy - I need to settle on astropy vs. pandas as soon as I get to metadata stage
    previous_galaxy_zoo_with_nsa = [row.to_dict() for _, row in previous_galaxy_zoo_with_nsa.iterrows()]  # neaten
    previous_galaxy_zoo_with_nsa = Table(previous_galaxy_zoo_with_nsa)

    # compare DR{} with previous subjects to see what's new
    return find_new_catalog_images(old_catalog=previous_galaxy_zoo_with_nsa, new_catalog=joint_catalog)


def main():
    """
    Run all steps to create the NSA-DECaLS-GZ catalog

    Returns:
        None
    """

    data_release = '5'
    fits_dir = '/data/galaxy_zoo/decals/fits/dr{}'.format(data_release)
    jpeg_dir = '/data/galaxy_zoo/decals/jpeg/dr{}'.format(data_release)

    nondefault_params = {
        'data_release': data_release,
        'fits_dir': fits_dir,
        'jpeg_dir': jpeg_dir,
        'nsa_version': '1_0_0'
    }
    s = Settings(**nondefault_params)

    # for safety, these settings must be separately specified
    s.new_catalog = True
    s.new_images = True
    s.overwrite_fits = False
    s.overwrite_jpeg = False
    s.run_to = 100

    nsa = get_nsa_catalog(s.nsa_catalog_loc)
    bricks = get_decals_bricks(s.bricks_loc, s.data_release)
    previous_subjects = pd.read_csv(s.subject_loc)  # previously extracted decals subjects

    catalog_to_upload = get_decals(nsa, bricks, previous_subjects, s)

    output_columns = ['nsa_id',
                      'iauname',
                      'ra',
                      'dec',
                      'galaxy_is_new',
                      'galaxy_is_updated',
                      'fits_ready',
                      'fits_filled',
                      'jpeg_ready']
    catalog_to_upload[output_columns].to_pandas().to_csv(s.upload_catalog_loc, index=False)

if __name__ == '__main__':
    main()
