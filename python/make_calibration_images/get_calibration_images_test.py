import pytest

import os

from get_calibration_images import *


TEST_EXAMPLES_DIR = 'python/test_examples'


@pytest.fixture
def jpeg_dir(tmpdir):
    return tmpdir.mkdir('jpeg_dir').strpath


@pytest.fixture
def galaxy(jpeg_dir):
    return {
        'fits_loc': '{}/example_a.fits'.format(TEST_EXAMPLES_DIR),
        'dr2_jpeg_loc': '{}/example_a_dr2.jpeg'.format(jpeg_dir),
        'colour_jpeg_loc': '{}/example_a_colour.jpeg'.format(jpeg_dir)
    }


def test_save_calibration_images_of_galaxy(galaxy):
    assert not os.path.exists(galaxy['dr2_jpeg_loc'])
    assert not os.path.exists(galaxy['colour_jpeg_loc'])
    save_calibration_images_of_galaxy(galaxy)
    assert os.path.exists(galaxy['dr2_jpeg_loc'])
    assert os.path.exists(galaxy['colour_jpeg_loc'])
