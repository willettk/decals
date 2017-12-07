import pytest

from astropy.table import Table

from manifest import *

TEST_EXAMPLES_DIR = 'python/test_examples'


@pytest.fixture()
def joint_catalog():
    return Table([
        {
            'nsa_id': 'example',
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
            'jpeg_loc': 'jpeg_here.jpg',
            'fits_loc': 'fits_there.fits'
        }
    ])


@pytest.fixture()
def calibration_catalog(joint_catalog):
    return Table([
        {
            'nsa_id': 'example',
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
            'jpeg_loc': 'do_not_use.jpg',
            'fits_loc': 'fits_there.fits',
            'dr2_jpeg_loc': 'type_dr2.jpg',
            'colour_jpeg_loc': 'type_colour.jpg'
        }
    ])


def test_create_manifest_from_joint_catalog(joint_catalog):
    manifest = create_manifest_from_joint_catalog(joint_catalog)
    assert len(manifest) == len(joint_catalog)
    entry = manifest[0]
    assert entry[0] == 'jpeg_here.jpg'
    assert type(entry[1]) == dict
    assert entry[1]['ra'] == 10.


def test_create_manifest_from_calibration_catalog(calibration_catalog):
    image_columns = ['dr2_jpeg_loc', 'colour_jpeg_loc']
    manifest = create_manifest_from_calibration_catalog(calibration_catalog, image_columns)
    assert len(manifest) == len(calibration_catalog) * 2

    first_entry = manifest[0]
    assert first_entry[0] == 'type_dr2.jpg'
    assert type(first_entry[1]) == dict
    assert first_entry[1]['ra'] == 10.

    last_entry = manifest[-1]
    assert last_entry[0] == 'type_colour.jpg'
    assert type(last_entry[1]) == dict
    assert last_entry[1]['ra'] == 10.
