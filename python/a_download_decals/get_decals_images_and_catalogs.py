# Find entries in the NASA-Sloan Atlas (NSA) that have a brick with images in Dark Energy Camera Legacy Survey (DECaLS)

import warnings

import pandas as pd
from astropy.io import fits
from astropy.table import Table

from get_catalogs.get_joint_nsa_decals_catalog import create_joint_catalog, get_nsa_catalog, get_decals_bricks
from get_catalogs.selection_cuts import apply_selection_cuts
from get_images.download_images_threaded import download_images_multithreaded
from setup.join_brick_tables import merge_bricks_catalogs

import download_decals_settings


def setup_tables(s):
    """
    Generate the 'bricks' and 'previous subjects' data tables used in main program.
    Only needs to be run once after downloading the required files.

    Args:
        s (Settings) object wrapper for user-specified settings

    Returns:
        None
    """
    if s.data_release == '3' or s.data_release == '5':
        coordinate_catalog = Table(fits.getdata(s.brick_coordinates_loc, 1))
        exposure_catalog = Table(fits.getdata(s.brick_exposures_loc, 1))
        bricks_catalog = merge_bricks_catalogs(coordinate_catalog, exposure_catalog)
        bricks_catalog.write(s.bricks_loc, overwrite=True)

    else:
        warnings.warn('Data release "{}" does not require joining brick tables - skipping'.format(s.data_release))


def get_decals(nsa=None, bricks=None, s=None):
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

    include_names = [
        'iauname',
        'nsa_id',
        'ra',
        'dec',
        'petrotheta',
        'petroth50',
        'petroth90',
        'z'
    ]

    if s.new_catalog:
        print('get new catalog')
        nsa_after_cuts = apply_selection_cuts(nsa)
        joint_catalog = create_joint_catalog(nsa_after_cuts, bricks, s.data_release, s.nsa_version, run_to=s.run_to)
        joint_catalog = joint_catalog[include_names]
        print('writing new catalog')
        joint_catalog.write(s.joint_catalog_loc, overwrite=True)
    else:
        print('get joint catalog')

        joint_catalog = Table(fits.getdata(s.joint_catalog_loc))[include_names]
        # joint_catalog = Table.read(s.joint_catalog_loc, data_end=1000, include_names=include_names)

    if s.new_images:
        print('get new images')
        joint_catalog = download_images_multithreaded(
            joint_catalog[:s.run_to],
            s.data_release,
            s.fits_dir,
            s.png_dir,
            overwrite_fits=s.overwrite_fits,
            overwrite_png=s.overwrite_png)

    return joint_catalog


def main():
    """
    Create the NSA-DECaLS-GZ catalog, download fits, produce png, and identify new subjects

    Returns:
        None
    """

    # specify setup options
    new_bricks_table = False  # if DR3+, run this on first use.

    # Setup tasks generate the 'bricks' data table used later.
    # They need only be completed once after downloading the required files
    if new_bricks_table:
        setup_tables(download_decals_settings)
        print('setup complete')

    # specify execution options
    download_decals_settings.new_catalog = False
    download_decals_settings.new_images = True
    download_decals_settings.overwrite_fits = False
    download_decals_settings.overwrite_png = False
    download_decals_settings.run_to = None

    nsa = get_nsa_catalog(download_decals_settings.nsa_catalog_loc, download_decals_settings.nsa_version)
    print('nsa loaded')
    bricks = get_decals_bricks(download_decals_settings.bricks_loc, download_decals_settings.data_release)
    print('catalogs loaded')

    # nsa = None
    # bricks = None
    joint_catalog = get_decals(nsa, bricks, download_decals_settings)

    joint_catalog.write(download_decals_settings.upload_catalog_loc, overwrite=True)


if __name__ == '__main__':
    main()
