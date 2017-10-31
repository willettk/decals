import pytest

from python.get_catalogs.get_joint_nsa_decals_catalog import *


@pytest.fixture()
def catalog():

    # DR5 is from 0-360 in RA and -20 to +30 in Dec

    gal_a1 = {'iauname': 'gal_a1',
              'ra': 146.714208787,
              'dec': -1.04128156958,  # in Dec window
              'petrotheta': 5.}  # above minimum petrotheta

    # adversarial example identical to gal_a
    gal_a2 = {'iauname': 'gal_a2',
              'ra': 146.714208787,
              'dec': -1.04128156958,  # in Dec window
              'petrotheta': 5.}  # above minimum petrotheta

    gal_a3 = {'iauname': 'gal_a3',
              'ra': 146.714208787,
              'dec': -1.04128156958,  # in Dec window
              'petrotheta': 1.}  # below minimum petrotheta

    gal_b = {'iauname': 'gal_b',
             'ra': 146.631735209,
             'dec': -50.988354858412,  # below allowed Dec window
             'petrotheta': 1.}  # below minimum petrotheta

    gal_c = {'iauname': 'gal_c',
             'ra': 147.176446949,
             'dec': 50.354030416643,  # above allowed Dec window
             'petrotheta': -1}  # below minimum petrotheta, adversarial example

    # in general sky area but not in any DECaLS brick
    gal_d = {'iauname': 'gal_d',
             'ra': 100.,
             'dec': -1.,
             'petrotheta': 5.}  # above minimum petrotheta

    return Table([gal_a1, gal_a2, gal_a3, gal_b, gal_c, gal_d])


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


def test_filter_catalog_to_approximate_sky_area_dr5(catalog, bricks, data_release):
    nsa_in_decals_area = filter_catalog_to_approximate_sky_area(catalog, bricks, data_release, visualise=False)
    nsa_in_decals_area_names = set(nsa_in_decals_area['iauname'])
    # gal a and d are in the correct general sky region, while b and c are not
    assert 'gal_a1' in nsa_in_decals_area_names
    assert 'gal_a2' in nsa_in_decals_area_names
    assert 'gal_a3' in nsa_in_decals_area_names
    assert 'gal_d' in nsa_in_decals_area_names
    assert 'gal_b' not in nsa_in_decals_area_names
    assert 'gal_c' not in nsa_in_decals_area_names


def test_create_joint_catalog(catalog, bricks, data_release, nsa_version):
    joint_catalog = create_joint_catalog(catalog, bricks, data_release, nsa_version, visualise=False)
    joint_catalog_names = set(joint_catalog['iauname'])
    # gal a is in the correct exact sky region (bricks)
    assert 'gal_a1' in joint_catalog_names
    assert 'gal_a2' in joint_catalog_names
    assert 'gal_a3' in joint_catalog_names
    # gal b and c are not in the correct general area
    assert 'gal_b' not in joint_catalog_names
    assert 'gal_c' not in joint_catalog_names
    # gal d is in the correct general area but NOT the exact sky region (i.e. no bricks cover that RA/DEC)
    assert 'gal_d' not in joint_catalog_names


# it's confusing to have this function in the same place as create joint catalog, which doesn't use it
def test_apply_selection_cuts(catalog):
    filtered_catalog = apply_selection_cuts(catalog)
    filtered_catalog_names = set(filtered_catalog['iauname'])
    # a1, a2 and d are above the minimum size while a3, b and c  are not
    assert 'gal_a1' in filtered_catalog_names
    assert 'gal_a2' in filtered_catalog_names
    assert 'gal_d' in filtered_catalog_names

    assert 'gal_a3' not in filtered_catalog_names
    assert 'gal_b' not in filtered_catalog_names
    assert 'gal_c' not in filtered_catalog_names
