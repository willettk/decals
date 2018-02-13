# TODO checking for updated images has been disabled - see https://github.com/zooniverse/decals/issues/19

import numpy as np
from astropy.io import fits

from shared_utilities import match_galaxies_to_catalog


def find_new_catalog_images(old_catalog, new_catalog):
    """
    Select only images from a new catalog that were a) missing or b) different in an older catalog
    Assumes catalogs names are unique

    Args:
        old_catalog (astropy.Table): old object catalog including fits_loc and IAUNAME, to be compared
        new_catalog (astropy.Table): new object catalog, as above, to be filtered for only new or changed images

    Returns:
        (astropy.Table) new catalog with only images that are new or changed since the old catalog
    """

    _, catalog_of_new_galaxies = match_galaxies_to_catalog(  # unmatched galaxies are new
        galaxies=new_catalog,
        catalog=old_catalog,
        galaxy_suffix='',
        catalog_suffix='_dr1_2')  # if field exists in both catalogs

    return catalog_of_new_galaxies


def fits_are_identical(fits_a_loc, fits_b_loc):
    """
    Given the location of two fits files, do they have identical pixels?
    TODO move elsewhere - to shared utilities?

    Args:
        fits_a_loc (str): location of one fits file
        fits_b_loc (str): location of other fits file

    Returns:
        (bool) True if both fits files have identical pixels (including shape), else False
    """
    pixels_a = fits.open(fits_a_loc)[0].data
    pixels_b = fits.open(fits_b_loc)[0].data
    return np.array_equal(pixels_a, pixels_b)
