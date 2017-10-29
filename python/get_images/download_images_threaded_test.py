import pytest

from python.get_images.download_images_threaded import *
from astropy.table import Table


@pytest.fixture
def fits_dir(tmpdir):
    return tmpdir.mkdir('fits_dir').strpath


@pytest.fixture
def jpeg_dir(tmpdir):
    return tmpdir.mkdir('jpeg_dir').strpath


@pytest.fixture
def fits_download_loc(fits_dir):
    return fits_dir + '/' + 'test_image.fits'


@pytest.fixture
def jpeg_loc(jpeg_dir):
    return jpeg_dir + '/' + 'test_image.jpg'


@pytest.fixture
def fits_image_params(fits_download_loc):
    return {
        'download_loc': fits_download_loc,
        'ra': 114.5970,
        'dec': 21.5681,
        'pixscale': 0.262,
        'size': 64,  # small test image for speed
        'data_release': '3'}


@pytest.fixture
def nsa_decals():
    gal_a = {'IAUNAME': 'J094651.40-010228.5',
             'RA': 146.714208787,
             'DEC': -1.04128156958,
             'PETROTH50': 3.46419,
             'PETROTH90': 10.4538}

    gal_b = {'IAUNAME': 'J094631.60-005917.7',
             'RA': 146.631735209,
             'DEC': -0.988354858412,
             'PETROTH50': 2.28976,
             'PETROTH90': 5.20297}

    gal_c = {'IAUNAME': 'J094842.33-002114.4',
             'RA': 147.176446949,
             'DEC': -0.354030416643,
             'PETROTH50': 7.16148,
             'PETROTH90': 24.7535}

    return Table([gal_a, gal_b, gal_c])



def test_download_fits_cutout_retrieves_fits(fits_image_params):
    download_fits_cutout(**fits_image_params)
    hdulist = fits.open(fits_image_params['download_loc'])
    matrix = hdulist[0].data
    assert matrix.shape == (3, fits_image_params['size'], fits_image_params['size'])
    assert np.max(matrix) > 0
    assert np.min(matrix) <= 0


def test_make_jpeg_from_fits_creates_jpeg(jpeg_loc):
    fits_loc = '../test_examples/example.fits'
    make_jpeg_from_fits(fits_loc, jpeg_loc)
    assert os.path.exists(jpeg_loc)
    # TODO Dustin's scaling method should be tested


def test_download_images_creates_fits_and_jpeg(nsa_decals, fits_dir, jpeg_dir):
    galaxy = nsa_decals[0]
    data_release = '3'
    download_path = '../test_examples'
    download_images(galaxy, fits_dir=fits_dir, jpeg_dir=jpeg_dir, data_release=data_release, overwrite=True)
    assert os.path.exists(get_fits_loc(download_path, galaxy))
    assert os.path.exists(get_jpeg_loc(download_path, galaxy))


def test_download_images_multithreaded(nsa_decals, fits_dir, jpeg_dir):
    data_release = '3'
    download_images_multithreaded(nsa_decals, data_release, fits_dir, jpeg_dir, overwrite=True)

# TODO test bad image counters with mock downloader
# TODO test time out counters with mock downloader
