import pytest

import pandas as pd
import numpy as np
from astropy.table import Table

from b_to_zooniverse.upload_decals import upload_decals_to_panoptes


@pytest.fixture
def calibration_dir(tmpdir):
    return tmpdir.mkdir('calibration_dir').strpath


@pytest.fixture
def fits_dir(tmpdir):
    return tmpdir.mkdir('fits_dir').strpath


@pytest.fixture
def png_dir(tmpdir):
    return tmpdir.mkdir('png_dir').strpath


@pytest.fixture
def fits_loc(fits_dir):
    return fits_dir + '/' + 'test_image.fits'


@pytest.fixture
def png_loc(png_dir):
    return png_dir + '/' + 'test_image.png'


@pytest.fixture()
def nsa_catalog():
    return Table([
        {'iauname': 'gal_a',
         'ra': 146.,
         'dec': -1.,
         'petroth50': 2.,
         'petroth90': 5.},

        # adversarial example identical to gal_a
        {'iauname': 'gal_dr1',
         'ra': 14.,
         'dec': -1.,
         'petroth50': 2.,
         'petroth90': 5.},

        {'iauname': 'gal_dr2',
         'ra': 1.,
         'dec': -1.,
         'petroth50': 2.,
         'petroth90': 5.}
    ])


@pytest.fixture()
def fake_metadata():
    # manifest expects many columns from joint catalog
    return {
        'petrotheta': 4.,
        'petroflux': [20., 21., 22., 23., 24., 25., 26.],
        'nsa_version': '1_0_0',
        'z': 0.1,
        'mag': [0., 1., 2., 3., 4., 5., 6.],
        'absmag': [10., 11., 12., 13., 14., 15., 16.],
        'nmgy': [30., 31., 32., 33., 34., 35., 36.],
        'another_column': 'sadness'}


@pytest.fixture()
def joint_catalog(fits_dir, png_dir, fake_metadata):
    # saved from downloader, which adds fits_loc, png_loc and png_ready to nsa_catalog + decals bricks
    gal_a = {
        'iauname': 'gal_a',
        'nsa_id': 0,
        'fits_loc': '{}/gal_a.fits'.format(fits_dir),
        'png_loc': '{}/gal_a.png'.format(png_dir),
        'png_ready': True,
        'ra': 146.,
        'dec': -1.,
        'petroth50': 2.,
        'petroth90': 5.}
    gal_a.update(fake_metadata)

    gal_dr1 = {
        'iauname': 'gal_dr1',
        'nsa_id': 1,
        'fits_loc': '{}/gal_b.fits'.format(fits_dir),
        'png_loc': '{}/gal_b.png'.format(png_dir),
        'png_ready': True,
        'ra': 14.,
        'dec': -1.,
        'petroth50': 2.,
        'petroth90': 5.}
    gal_dr1.update(fake_metadata)

    gal_dr2 = {
        'iauname': 'gal_dr2',
        'nsa_id': 2,
        'fits_loc': '{}/gal_c.fits'.format(fits_dir),
        'png_loc': '{}/gal_c.png'.format(png_dir),
        'png_ready': True,
        'ra': 1.,
        'dec': -1.,
        'petroth50': 2.,
        'petroth90': 5.}
    gal_dr2.update(fake_metadata)

    return Table([gal_a, gal_dr1, gal_dr2])


@pytest.fixture()
def previous_subjects():
    # loaded from GZ data dump. Metadata already corrected in setup (name = main stage).
    return Table([

        # DR1 entries have provided_image_id filled with iau name, nsa_id blank, dr blank
        {'_id': 'ObjectId(0)',
         'zooniverse_id': 'gz_dr1',
         'iauname': 'gal_dr1',
         'nsa_id': 1,
         'dr': 'DR1',
         'ra': 14.,
         'dec': -1.,
         'petroth50': 2.,
         'petroth90': 5.
         },


        # DR2 entries have provided_image_id blank, nsa_id filled with NSA_[number], dr filled with 'DR2'
        {'_id': 'ObjectId(1)',
         'zooniverse_id': 'gz_dr2',
         'iauname': 'gal_dr2',
         'nsa_id': 2,
         'dr': 'DR2',
         'ra': 1.,
         'dec': -1.,
         'petroth50': 2.,
         'petroth90': 5.
         }
    ])


@pytest.fixture()
def expert_catalog():
    return Table([
        {
            # gal a is both a bar and ring galaxy, and so should be included in the calibration set
            'iauname': 'gz_a',
            'ra': 146.,
            'dec': -1.,
            'bar': 2 ** 5,
            'ring': 2 ** 3,
        }
    ])


def test_upload_decals_to_panoptes(joint_catalog, previous_subjects, expert_catalog, calibration_dir):
    main_subjects, calibration_subjects = upload_decals_to_panoptes(
        joint_catalog, previous_subjects, expert_catalog, calibration_dir, subject_set_name='TEST')

    print(main_subjects)
    print(calibration_subjects)

    assert len(main_subjects) == 1
    assert len(calibration_subjects) == len(main_subjects) * 2

    first_main_subject = main_subjects[0]
    assert first_main_subject['png_loc'][-17:] == 'png_dir/gal_a.png'
    assert first_main_subject['key_data']['ra'] == 146.0
    assert first_main_subject['key_data']['dec'] == -1.
    assert first_main_subject['key_data']['nsa_id'] == 0
    assert first_main_subject['key_data']['petroth50'] == 2.0
    assert first_main_subject['key_data']['mag_abs_r'] == 14.0

    # wrong, should have 1 of each version not two
    assert calibration_subjects[0]['png_loc'][-29:] == 'calibration_dir/gal_a_dr2.png'
    assert calibration_subjects[0]['key_data']['selected_image'] == 'dr2_png_loc'
    assert calibration_subjects[1]['png_loc'][-32:] == 'calibration_dir/gal_a_colour.png'
