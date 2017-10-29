import pytest

from astropy.table import Table

from python.get_catalogs.get_dr3_catalog import fits_are_identical, find_new_catalog_images


@pytest.fixture()
def nsa_decals_dr3():

    # galaxy a is a new DR3 image, only in this catalog
    gal_a_dr3 = {'IAUNAME': 'gal_a',
                 'RA': 146.714208787,
                 'DEC': -1.04128156958,
                 'PETROTH50': 3.46419,
                 'PETROTH90': 10.4538,
                 'fits_loc': '../test_examples/example_a.fits'}

    # galaxy b is in dr2, but with a different fits file
    # it should be included in the gz dr3 upload
    gal_b_dr3 = {'IAUNAME': 'gal_b',
                 'RA': 146.631735209,
                 'DEC': -0.988354858412,
                 'PETROTH50': 2.28976,
                 'PETROTH90': 5.20297,
                 'fits_loc': '../test_examples/example_b.fits'}

    # galaxy c is in dr2 with an identical fits file
    # it should not be included in the gz dr3 upload (already classified)
    gal_c_dr3 = {'IAUNAME': 'gal_c',
                 'RA': 147.176446949,
                 'DEC': -0.354030416643,
                 'PETROTH50': 7.16148,
                 'PETROTH90': 24.7535,
                 'fits_loc': '../test_examples/example_c.fits'}

    return Table([gal_a_dr3, gal_b_dr3, gal_c_dr3])


@pytest.fixture()
def nsa_decals_dr2():
    # same properties as dr3 except fits file
    gal_b_dr2 = {'IAUNAME': 'gal_b',
                 'RA': 146.631735209,
                 'DEC': -0.988354858412,
                 'PETROTH50': 2.28976,
                 'PETROTH90': 5.20297,
                 'fits_loc': '../test_examples/example_c.fits'}  # points to dif. fits file than dr3

    # same properties as dr3 including fits file
    gal_c_dr2 = {'IAUNAME': 'gal_c',
                 'RA': 147.176446949,
                 'DEC': -0.354030416643,
                 'PETROTH50': 7.16148,
                 'PETROTH90': 24.7535,
                 'fits_loc': '../test_examples/example_c.fits'}  # points to same fits file as dr3

    return Table([gal_b_dr2, gal_c_dr2])


def test_fits_are_identical():
    fits_a_loc = '../test_examples/example_a.fits'
    fits_b_loc = '../test_examples/example_b.fits'
    assert fits_are_identical(fits_a_loc, fits_a_loc)
    assert not fits_are_identical(fits_a_loc, fits_b_loc)


def test_find_new_catalog_images(nsa_decals_dr2, nsa_decals_dr3):
    only_new_galaxies = find_new_catalog_images(old_catalog=nsa_decals_dr2, new_catalog=nsa_decals_dr3)
    new_galaxy_names = set(only_new_galaxies['IAUNAME'])
    assert 'gal_a' in new_galaxy_names
    assert 'gal_b' in new_galaxy_names
    assert 'gal_c' not in new_galaxy_names
