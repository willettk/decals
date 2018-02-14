import pytest
from astropy.table import Table

from b_to_zooniverse.do_upload import manifest

TEST_EXAMPLES_DIR = 'python/test_examples'


@pytest.fixture()
def joint_catalog():
    return Table([
        {
            'nsa_id': 'example_nsa_id',
            'iauname': 'example_iauname',
            'ra': 10.0,
            'dec': 12.0,
            'petroth50': 50.,
            'petrotheta': 4.,
            'petroflux': [20., 21., 22., 23., 24., 25., 26.],
            'nsa_version': '1_0_0',
            'z': 0.1,
            'mag': [0., 1., 2., 3., 4., 5., 6.],
            'absmag': [10., 11., 12., 13., 14., 15., 16.],
            'nmgy': [30., 31., 32., 33., 34., 35., 36.],
            'png_loc': 'jpeg_here.png',
            'fits_loc': 'fits_there.fits'
        }
    ])


@pytest.fixture()
def calibration_catalog():
    return Table([
        {
            'nsa_id': 'example_nsa_id',
            'iauname': 'example_iauname',
            'ra': 10.0,
            'dec': 12.0,
            'petroth50': 50.,
            'petrotheta': 4.,
            'petroflux': [20., 21., 22., 23., 24., 25., 26.],
            'nsa_version': '1_0_0',
            'z': 0.1,
            'mag': [0., 1., 2., 3., 4., 5., 6.],
            'absmag': [10., 11., 12., 13., 14., 15., 16.],
            'nmgy': [30., 31., 32., 33., 34., 35., 36.],
            'png_loc': 'do_not_use.png',
            'fits_loc': 'fits_there.fits',
            'dr2_png_loc': 'type_dr2.png',
            'colour_png_loc': 'type_colour.png'
        }
    ])


def test_create_manifest_from_joint_catalog(joint_catalog):
    new_manifest = manifest.create_manifest_from_joint_catalog(joint_catalog)
    assert len(new_manifest) == len(joint_catalog)
    entry = new_manifest[0]
    assert entry['png_loc'] == 'jpeg_here.png'
    assert type(entry['key_data']) == dict
    assert entry['key_data']['ra'] == 10.


def test_create_manifest_from_calibration_catalog(calibration_catalog):
    image_columns = ['dr2_png_loc', 'colour_png_loc']
    new_manifest = manifest.create_manifest_from_calibration_catalog(calibration_catalog, image_columns)
    assert len(new_manifest) == len(calibration_catalog) * 2

    first_entry = new_manifest[0]
    assert first_entry['png_loc'] == 'type_dr2.png'
    assert type(first_entry['key_data']) == dict
    assert first_entry['key_data']['ra'] == 10.

    last_entry = new_manifest[-1]
    assert last_entry['png_loc'] == 'type_colour.png'
    assert type(last_entry['key_data']) == dict
    assert last_entry['key_data']['ra'] == 10.
