import pytest

from python.get_images.download_images_threaded import *


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
    # TODO create a short reference catalog here or in test_examples
    # a = [1, 4, 5]
    # b = [2.0, 5.0, 8.2]
    # c = ['x', 'y', 'z']
    # t = Table([a, b, c], names=('a', 'b', 'c'), meta={'name': 'first table'})

    dr = '3'
    nsa_version = '0_1_2'
    return Table(fits.getdata('../../fits/nsa_v{0}_decals_dr{1}.fits'.format(nsa_version, dr), 1))[:3]


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
