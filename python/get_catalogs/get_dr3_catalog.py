
import numpy as np
from astropy.io import fits
from astropy.table import Table

from tqdm import tqdm


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

    galaxy_is_new = np.zeros(len(new_catalog), dtype=bool)
    galaxy_is_updated = np.zeros(len(new_catalog), dtype=bool)

    old_catalog_names = set(old_catalog['IAUNAME'])

    # TODO could make this unordered map to increase speed
    for row_index, new_galaxy in tqdm(enumerate(new_catalog), total=len(new_catalog)):

        # if the name is not in the old catalog, the galaxy must be new
        if new_galaxy['IAUNAME'] not in old_catalog_names:
            galaxy_is_new[row_index] = True

        # if the pixels have changed, the galaxy should be considered new and re-classified
        else:
            old_galaxy = old_catalog[old_catalog['IAUNAME'] == new_galaxy['IAUNAME']][0]
            new_fits_loc = new_galaxy['fits_loc']
            old_fits_loc = old_galaxy['fits_loc']
            try:
                identical = fits_are_identical(new_fits_loc, old_fits_loc)
            except TypeError or KeyError:  # TODO handle the case where either image is bad
                identical = False
            if not identical:
                galaxy_is_updated[row_index] = True

    new_catalog['galaxy_is_new'] = galaxy_is_new
    new_catalog['galaxy_is_updated'] = galaxy_is_updated

    print(new_catalog[['IAUNAME', 'galaxy_is_new', 'galaxy_is_updated']])

    return new_catalog[galaxy_is_new | galaxy_is_updated]


def fits_are_identical(fits_a_loc, fits_b_loc):
    """
    Given the location of two fits files, do they have identical pixels?

    Args:
        fits_a_loc (str): location of one fits file
        fits_b_loc (str): location of other fits file

    Returns:
        (bool) True if both fits files have identical pixels (including shape), else False
    """
    pixels_a = fits.open(fits_a_loc)[0].data
    pixels_b = fits.open(fits_b_loc)[0].data
    return np.array_equal(pixels_a, pixels_b)


if __name__ == '__main__':

    catalog_dir = '/data/galaxy_zoo/decals/catalogs'
    nsa_version = '0_1_2'

    dr2_loc = '{0}/nsa_v{1}_decals_dr{2}.fits'.format(catalog_dir, nsa_version, '2')
    dr3_loc = '{0}/nsa_v{1}_decals_dr{2}.fits'.format(catalog_dir, nsa_version, '3')

    dr2 = Table(fits.getdata(dr2_loc, 1))
    dr3 = Table(fits.getdata(dr3_loc, 1)[:100])

    dr3_only_new = find_new_catalog_images(old_catalog=dr2, new_catalog=dr3)
    print(len(dr3_only_new))
    dr3_only_new.write('{}/galaxy_zoo_upload.fits'.format(catalog_dir), overwrite=True)
