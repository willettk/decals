# Find entries in the NASA-Sloan Atlas (NSA) that have a brick with images in Dark Energy Camera Legacy Survey (DECaLS)

import warnings

import pandas as pd
from astropy.io import fits
from astropy.table import Table

from get_catalogs.get_joint_nsa_decals_catalog import create_joint_catalog, get_nsa_catalog, get_decals_bricks
from get_catalogs.selection_cuts import apply_selection_cuts
from get_images.download_images_threaded import download_images_multithreaded
from setup.join_brick_tables import merge_bricks_catalogs


class Settings():
    # Define parameters, file locations, etc.

    def get_nsa_catalog_loc(self):
        """
        Returns:
            (str) location of NASA-Sloan Atlas (NSA) catalog. Input file.
        """
        return '{}/nsa_v{}.fits'.format(self.catalog_dir, self.nsa_version)

    def get_joint_catalog_loc(self):
        """
        Returns:
            (str) location of catalog of all NSA galaxies imaged by DECALS DR{}. Output file.
        """
        return '{0}/nsa_v{1}_decals_dr{2}.fits'.format(
            self.catalog_dir, self.nsa_version, self.data_release)

    def get_upload_catalog_loc(self):
        """
        Returns:
            (str) location of catalog of all NSA galaxies imaged by DECALS DR{} and not yet classified. Output file.
        """
        return '{}/dr{}_nsa{}_to_upload.fits'.format(self.catalog_dir, self.data_release, self.nsa_version)

    def get_bricks_loc(self):
        """
        Returns:
            (str) location of bricks catalog, including center and exposure counts.
                Input file for DR1, DR2. Calculated from brick_coordinates_loc and brick_exposures_loc for DR5.
        """
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

    def get_brick_coordinates_loc(self):
        """
        Returns:
            (str) location of bricks coordinate catalog, including center and edges. Input file for DR5.
                Named 'survey-bricks.fits' on DR5 website
        """
        if self.data_release == '5' or self.data_release == '3':
            return '{}/survey-bricks.fits'.format(self.catalog_dir)

    def get_brick_exposures_loc(self):
        """
        Returns:
            (str) location of bricks exposure catalog, including brick centers and exposure counts.
                Input file for DR5. Named 'survey-bricks-dr5.fits' on DR5 website
        """
        if self.data_release == '5' or self.data_release == '3':
            return '{}/survey-bricks-dr5.fits'.format(self.catalog_dir)

    def derive_file_paths(self):
        self.nsa_catalog_loc = self.get_nsa_catalog_loc()
        self.joint_catalog_loc = self.get_joint_catalog_loc()
        self.brick_coordinates_loc = self.get_brick_coordinates_loc()
        self.brick_exposures_loc = self.get_brick_exposures_loc()
        self.bricks_loc = self.get_bricks_loc()
        self.upload_catalog_loc = self.get_upload_catalog_loc()
        self.brick_coordinates_loc = self.get_brick_coordinates_loc()
        self.brick_exposures_loc = self.get_brick_exposures_loc()

    def __init__(self,
                 # must specify these for safety
                 fits_dir,
                 png_dir,
                 data_release,
                 # these rarely change and are set as default
                 nsa_version='1_0_0',
                 catalog_dir='/data/galaxy_zoo/decals/catalogs',
                 subject_loc='/data/galaxy_zoo/decals/subjects/decals_dr1_and_dr2.csv',
                 run_to=-1):

        self.data_release = data_release
        self.nsa_version = nsa_version
        self.catalog_dir = catalog_dir
        self.fits_dir = fits_dir
        self.png_dir = png_dir
        self.subject_loc = subject_loc
        self.run_to = run_to

        self.derive_file_paths()


def setup_tables(s):
    """
    Generate the 'bricks' and 'previous subjects' data tables used in main program.
    Only needs to be run once after downloading the required files.

    Args:
        s (Settings) object wrapper for user-specified settings

    Returns:
        None
    """
    if s.merge_bricks:
        if s.data_release == '3' or s.data_release == '5':
            coordinate_catalog = Table(fits.getdata(s.brick_coordinates_loc, 1))
            exposure_catalog = Table(fits.getdata(s.brick_exposures_loc, 1))
            bricks_catalog = merge_bricks_catalogs(coordinate_catalog, exposure_catalog)
            bricks_catalog.write(s.brick_loc, overwrite=True)

        else:
            warnings.warn('Data release "{}" does not require joining brick tables - skipping'.format(s.data_release))


def get_decals(nsa, bricks, s):
    """
    Find NSA galaxies imaged by DECALS. Download fits. Create RGB pngs. Return catalog of new galaxies.

    Args:
        nsa (astropy.Table): catalog of NSA galaxies
        bricks (astropy.Table): catalog of DECALS bricks, including center, edges and exposure counts
        previous_subjects (pd.DataFrame):
        s (Settings): object wrapper for user-specified behaviour

    Returns:
        (pd.DataFrame) catalog of DECaLS galaxies not previously classified in Galaxy Zoo
    """

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
            s.png_dir,
            overwrite_fits=s.overwrite_fits,
            overwrite_png=s.overwrite_png)
        joint_catalog.write(s.joint_catalog_loc, overwrite=True)

    return joint_catalog


def main():
    """
    Create the NSA-DECaLS-GZ catalog, download fits, produce png, and identify new subjects

    Returns:
        None
    """

    data_release = '5'
    nsa_version = '1_0_0'
    fits_dir = '/Volumes/external/decals/fits/dr{}'.format(data_release)
    png_dir = '/Volumes/external/decals/png/dr{}'.format(data_release)
    new_bricks_table = False

    nondefault_params = {
        'nsa_version': nsa_version,
        'data_release': data_release,
        'fits_dir': fits_dir,
        'png_dir': png_dir,
    }
    s = Settings(**nondefault_params)

    # specify setup options
    s.merge_bricks = False

    # Setup tasks generate the 'bricks' data table used later.
    # They need only be completed once after downloading the required files
    if new_bricks_table:
        setup_tables(s)

    # specify execution options
    s.new_catalog = True
    s.new_images = True
    s.overwrite_fits = False
    s.overwrite_png = False
    s.run_to = -1

    nsa = get_nsa_catalog(s.nsa_catalog_loc, nsa_version)
    bricks = get_decals_bricks(s.bricks_loc, s.data_release)

    joint_catalog = get_decals(nsa, bricks, s)

    joint_catalog.write(s.upload_catalog_loc)


if __name__ == '__main__':
    main()
