import numpy as np
import pytest
from PIL import Image
from astropy.io import fits as fits

from python.get_images.image_server import download_image


class TestImageServer():
    @pytest.fixture
    def download_path(self, tmpdir, scope='session'):  # tmpdir is built-in fixture
        return tmpdir.mkdir('tmpdir').strpath  # create directory, get string path

    @pytest.fixture
    def download_loc(self, download_path):
        return download_path + '/' + 'test_image'

    @pytest.fixture
    def image_params(self, download_loc):
        return {
            'download_loc': download_loc,
            'ra': 114.5970,
            'dec': 21.5681,
            'pixscale': 0.262,
            'size': 512}

    def test_image_server_retrieves_jpg(self, image_params):
        download_image(type='jpg', **image_params)
        im = Image.open(image_params['download_loc'] + '.jpg')
        matrix = np.asarray(im)
        assert matrix.shape == (image_params['size'], image_params['size'], 3)
        assert np.max(matrix) > 0
        assert np.min(matrix) <= 0

    def test_image_server_retrieves_fits(self, image_params):
        download_image(type='fits', **image_params)
        hdulist = fits.open(image_params['download_loc'] + '.fits')
        matrix = hdulist[0].data
        assert matrix.shape == (3, image_params['size'], image_params['size'])
        assert np.max(matrix) > 0
        assert np.min(matrix) <= 0


pytest.main()
