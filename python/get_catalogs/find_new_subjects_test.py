import pytest

from astropy.table import Table

from get_catalogs.find_new_subjects import fits_are_identical, find_new_catalog_images


@pytest.fixture()
def nsa_decals_dr3():
    return Table([
        # galaxy a is a new DR3 image, only in this catalog
        {'iauname': 'gal_a',
         'ra': 146.714208787,
         'dec': -1.04128156958,
         'petroth50': 3.46419,
         'petroth90': 10.4538,
         'fits_loc': '../test_examples/example_a.fits'},

        # galaxy b is in dr2, but with a different fits file
        # it should be included in the gz dr3 upload
        {'iauname': 'gal_b',
         'ra': 146.631735209,
         'dec': -0.988354858412,
         'petroth50': 2.28976,
         'petroth90': 5.20297,
         'fits_loc': '../test_examples/example_b.fits'},

        # galaxy c is in dr2 with an identical fits file
        # it should not be included in the gz dr3 upload (already classified)
        {'iauname': 'gal_c',
         'ra': 147.176446949,
         'dec': -0.354030416643,
         'petroth50': 7.16148,
         'petroth90': 24.7535,
         'fits_loc': '../test_examples/example_c.fits'},

        # same properties as dr3 including fits file
        {'iauname': 'gal_missing_in_dr2',
         'ra': 147.176446949,
         'dec': -0.354030416643,
         'petroth50': 7.16148,
         'petroth90': 24.7535,
         'fits_loc': '../test_examples/example_c.fits'}  # good path, unlike in dr2
    ])


@pytest.fixture()
def nsa_decals_dr2():
    # same properties as dr3 except fits file
    return Table([
        {'iauname': 'gal_b',
         'ra': 146.631735209,
         'dec': -0.988354858412,
         'petroth50': 2.28976,
         'petroth90': 5.20297,
         'fits_loc': '../test_examples/example_c.fits'},  # points to dif. fits file than dr3

        # same properties as dr3 including fits file
        {'iauname': 'gal_c',
         'ra': 147.176446949,
         'dec': -0.354030416643,
         'petroth50': 7.16148,
         'petroth90': 24.7535,
         'fits_loc': '../test_examples/example_c.fits'},  # points to same fits file as dr3

        # same properties as dr3 including fits file
        {'iauname': 'gal_missing_in_dr2',
         'ra': 147.176446949,
         'dec': -0.354030416643,
         'petroth50': 7.16148,
         'petroth90': 24.7535,
         'fits_loc': 'no_image_here'}  # points nowhere
    ])


def test_fits_are_identical():
    fits_a_loc = '../test_examples/example_a.fits'
    fits_b_loc = '../test_examples/example_b.fits'
    assert fits_are_identical(fits_a_loc, fits_a_loc)
    assert not fits_are_identical(fits_a_loc, fits_b_loc)


def test_find_new_catalog_images(nsa_decals_dr2, nsa_decals_dr3):
    only_new_galaxies = find_new_catalog_images(old_catalog=nsa_decals_dr2, new_catalog=nsa_decals_dr3)
    new_galaxy_names = set(only_new_galaxies['iauname'])
    assert 'gal_a' in new_galaxy_names
    assert 'gal_b' in new_galaxy_names
    assert 'gal_missing_in_dr2' in new_galaxy_names
    assert 'gal_c' not in new_galaxy_names
