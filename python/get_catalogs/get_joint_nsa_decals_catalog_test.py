import pytest

from python.get_catalogs.get_joint_nsa_decals_catalog import *


@pytest.fixture()
def catalog():

    # DR5 is from 0-360 in RA and -20 to +30 in Dec
    return Table([
        {'iauname': 'gal_a',
         'ra': 150.,
         'dec': 3.,  # in Dec window
         'petrotheta': 5.},  # above minimum petrotheta

        # adversarial example identical to gal_a
        {'iauname': 'gal_a_duplicate',
         'ra': 150.,
         'dec': 3.,  # in Dec window
         'petrotheta': 5.},  # above minimum petrotheta

        {'iauname': 'gal_a_small',
         'ra': 150.,
         'dec': 3.,  # in Dec window
         'petrotheta': 1.},  # below minimum petrotheta

        {'iauname': 'gal_a_negative_size',
         'ra': 150.,
         'dec': 3.,  # in Dec window
         'petrotheta': -1.},  # below minimum petrotheta, adversarial example

        # matches to multiple fake bricks
        {'iauname': 'gal_multiple_bricks',
         'ra': 280.,
         'dec': -1.,
         'petrotheta': 5.},

        {'iauname': 'gal_below_bricks',
         'ra': 150.,
         'dec': -50.,  # below allowed Dec window
         'petrotheta': 1.},  # below minimum petrotheta

        {'iauname': 'gal_above_bricks',
         'ra': 150.,
         'dec': 50.,  # above allowed Dec window
         'petrotheta': 0.1},

        # in general sky area but not in any DECaLS brick
        {'iauname': 'gal_outside_bricks',
         'ra': 100.,
         'dec': -1.,
         'petrotheta': 0.1},
    ])


@pytest.fixture()
def data_release():
    return '5'


@pytest.fixture()
def nsa_version():
    return '0_1_2'


@pytest.fixture()
def bricks(data_release):
    catalog_dir = '/data/galaxy_zoo/decals/catalogs'
    if data_release == '5' or data_release == '3':
        bricks_filename = 'survey-bricks-dr{}-with-coordinates.fits'.format(data_release)
    elif data_release == '2':
        bricks_filename = 'decals-bricks-dr2.fits'
    elif data_release == '1':
        bricks_filename = 'decals-bricks-dr1.fits'
    else:
        raise ValueError('Data Release "{}" not recognised'.format(data_release))
    bricks_loc = '{}/{}'.format(catalog_dir, bricks_filename)
    return get_decals_bricks(bricks_loc, data_release)


@pytest.fixture()
def bricks(data_release):
    return Table([
        # galaxies RA 0->90 and dec -10->10 should match
        {'ra1': 0.,
         'ra': 45.,
         'ra2': 90.,
         'dec1': -10.,
         'dec': 0.,
         'dec2': 10.},

        # galaxies RA 90->110 should fail to match
        # galaxies RA 110->270 and dec - 10->10 should match
        {'ra1': 110.,
         'ra': 185.,
         'ra2': 360.,
         'dec1': -10.,
         'dec': 0.,
         'dec2': 10.},

        # galaxies RA = 270->360 should multiple-match
        {'ra1': 270.,
         'ra': 315.,
         'ra2': 360.,
         'dec1': -10.,
         'dec': 0.,
         'dec2': 10.},
    ])


def test_filter_catalog_to_approximate_sky_area_dr5(catalog, bricks, data_release):
    nsa_in_decals_area = filter_catalog_to_approximate_sky_area(catalog, bricks, data_release, visualise=False)
    remaining_names = set(nsa_in_decals_area['iauname'])

    assert 'gal_a' in remaining_names
    assert 'gal_a_duplicate' in remaining_names

    assert 'gal_outside_bricks' in remaining_names
    assert 'gal_multiple_bricks' in remaining_names

    assert 'gal_below_bricks' not in remaining_names
    assert 'gal_above_bricks' not in remaining_names


def test_create_joint_catalog(catalog, bricks, data_release, nsa_version):
    joint_catalog = create_joint_catalog(catalog, bricks, data_release, nsa_version, visualise=False)
    remaining_names = set(joint_catalog['iauname'])

    assert 'gal_a' in remaining_names
    assert 'gal_a_duplicate' in remaining_names
    assert 'gal_multiple_bricks' in remaining_names  # for multi-match galaxies, the first match is recorded

    assert 'gal_below_bricks' not in remaining_names
    assert 'gal_above_bricks' not in remaining_names

    assert 'gal_outside_bricks' not in remaining_names


# it's confusing to have this function in the same place as create joint catalog, which doesn't use it
def test_apply_selection_cuts(catalog):
    filtered_catalog = apply_selection_cuts(catalog)
    filtered_catalog_names = set(filtered_catalog['iauname'])

    assert 'gal_a' in filtered_catalog_names
    assert 'gal_a_duplicate' in filtered_catalog_names

    assert 'gal_a_small' not in filtered_catalog_names
    assert 'gal_a_negative_size' not in filtered_catalog_names
