import pytest

import os

from get_calibration_images import *


TEST_EXAMPLES_DIR = 'python/test_examples'


@pytest.fixture
def png_dir(tmpdir):
    return tmpdir.mkdir('png_dir').strpath


@pytest.fixture
def galaxy(png_dir):
    return {
        'fits_loc': '{}/example_a.fits'.format(TEST_EXAMPLES_DIR),
        'dr2_png_loc': '{}/example_a_dr2.png'.format(png_dir),
        'colour_png_loc': '{}/example_a_colour.png'.format(png_dir)
    }


def test_save_calibration_images_of_galaxy(galaxy):
    assert not os.path.exists(galaxy['dr2_png_loc'])
    assert not os.path.exists(galaxy['colour_png_loc'])
    save_calibration_images_of_galaxy(galaxy)
    assert os.path.exists(galaxy['dr2_png_loc'])
    assert os.path.exists(galaxy['colour_png_loc'])
