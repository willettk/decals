import pytest
from astropy.table import Table

from get_images.download_images_threaded import *

TEST_EXAMPLES_DIR = 'python/test_examples'


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


@pytest.fixture
def fits_image_params(fits_loc):
    return {
        'fits_loc': fits_loc,
        'ra': 114.5970,
        'dec': 21.5681,
        'pixscale': 0.262,
        'size': 64,  # small test image for speed
        'data_release': '5'}


@pytest.fixture()
def new_galaxy(fits_loc, png_loc):
    return {
        'iauname': 'new_galaxy',
        'fits_loc': fits_loc,
        'png_loc': png_loc,
        'ra': 146.714208787,
        'dec': -1.04128156958,
        'petroth50': 3.46419,
        'petroth90': 10.4538
    }


@pytest.fixture
def nsa_decals(fits_dir, png_dir):
    # TODO bad practice: fits_loc is overwritten in some tests
    gal_missing = {
        'iauname': 'J0',
        'fits_loc': '{}/example_missing.fits'.format(fits_dir),
        'png_loc': '{}/example_missing.png'.format(png_dir),
        'ra': 146.714208787,
        'dec': -1.04128156958,
        'petroth50': 3.46419,
        'petroth90': 10.4538}

    gal_incomplete = {
        'iauname': 'J1',
        'fits_loc': '{}/example_incomplete.fits'.format(TEST_EXAMPLES_DIR),
        'png_loc': '{}/example_incomplete.png'.format(TEST_EXAMPLES_DIR),
        'ra': 146.714208787,
        'dec': -1.04128156958,
        'petroth50': 3.46419,
        'petroth90': 10.4538}

    gal_bad_pix = {
        'iauname': 'J2',
        'fits_loc': '{}/example_bad_pix.fits'.format(TEST_EXAMPLES_DIR),
        'png_loc': '{}/example_bad_pix.png'.format(TEST_EXAMPLES_DIR),
        'ra': 146.714208787,
        'dec': -1.04128156958,
        'petroth50': 3.46419,
        'petroth90': 10.4538}

    gal_a = {
        'iauname': 'J3',
        'fits_loc': '{}/example_a.fits'.format(TEST_EXAMPLES_DIR),
        'png_loc': '{}/example_a.png'.format(TEST_EXAMPLES_DIR),
        'ra': 146.714208787,
        'dec': -1.04128156958,
        'petroth50': 3.46419,
        'petroth90': 10.4538}

    gal_b = {
        'iauname': 'J4',
        'fits_loc': '{}/example_b.fits'.format(TEST_EXAMPLES_DIR),
        'png_loc': '{}/example_b.png'.format(TEST_EXAMPLES_DIR),
        'ra': 146.631735209,
        'dec': -0.988354858412,
        'petroth50': 2.28976,
        'petroth90': 5.20297}

    gal_c = {
        'iauname': 'J5',
        'fits_loc': '{}/example_c.fits'.format(TEST_EXAMPLES_DIR),
        'png_loc': '{}/example_c.png'.format(TEST_EXAMPLES_DIR),
        # 'ra': 147.176446949,
        # 'dec': -0.354030416643,
        'ra': 146.631735209,
        'dec': -0.988354858412,
        'petroth50': 7.16148,
        'petroth90': 24.7535}

    gal_mystery = {
        'iauname': 'J121457.39-002412.7',
        'fits_loc': '/data/temp/J121457.39-002412.7.fits',
        'png_loc': '/data/temp/J121457.39-002412.7.png',
        'ra': 183.73974251,
        'dec': -0.403807422075,
        'petroth50': 7.16148,
        'petroth90': 24.7535}
        # http: // legacysurvey.org / viewer / fits - cutout?ra =  & dec =  & pixscale = 0.1 & size = 424 & layer = decals - dr5

    return Table([gal_missing, gal_incomplete, gal_bad_pix, gal_a, gal_b, gal_c, gal_mystery])


@pytest.fixture()
def partially_downloaded_galaxy_params(tmpdir):
    # copy the incomplete fits to a temporary directory
    dir = tmpdir.mkdir('incomplete_fits_dir').strpath
    temp_incomplete_example_loc = '{}/temp_incomplete_example.fits'.format(dir)
    incomplete_example_loc = '../test_examples/example_incomplete.fits'
    os.system('cp {} {}'.format(incomplete_example_loc, temp_incomplete_example_loc))
    # arguments required by download_fits_cutout
    return {
        'fits_loc': temp_incomplete_example_loc,
        'ra': 114.5970,
        'dec': 21.5681,
        'pixscale': 0.262,
        'size': 64,  # small test image for speed
        'data_release': '5'}


def test_download_fits_cutout_retrieves_fits(fits_image_params):
    download_fits_cutout(**fits_image_params)
    hdulist = fits.open(fits_image_params['fits_loc'])
    matrix = hdulist[0].data
    assert matrix.shape == (3, fits_image_params['size'], fits_image_params['size'])
    assert np.max(matrix) > 0
    assert np.min(matrix) <= 0


def test_download_fits_cutout_corrects_incomplete_fits(partially_downloaded_galaxy_params):
    assert not fits_downloaded_correctly(partially_downloaded_galaxy_params['fits_loc'])
    download_fits_cutout(**partially_downloaded_galaxy_params)
    assert fits_downloaded_correctly(partially_downloaded_galaxy_params['fits_loc'])


def test_get_download_quality_of_fits_on_missing_fits():
    downloaded, good_pix = get_download_quality_of_fits('nothing here')
    assert not downloaded
    assert not good_pix


def test_get_download_quality_of_fits_on_incomplete_fits():
    incomplete_fits_loc = '{}/example_incomplete.fits'.format(TEST_EXAMPLES_DIR)
    assert os.path.exists(incomplete_fits_loc)  # file exists
    downloaded, good_pix = get_download_quality_of_fits(incomplete_fits_loc)
    assert not downloaded  # file exists but is not completely downloaded
    assert not good_pix


def test_get_download_quality_of_fits_on_bad_pix_fits():
    partial_fits_loc = '{}/example_bad_pix.fits'.format(TEST_EXAMPLES_DIR)
    assert os.path.exists(partial_fits_loc)
    downloaded, good_pix = get_download_quality_of_fits(partial_fits_loc)
    assert downloaded
    assert not good_pix


def test_make_png_from_fits_creates_png(png_loc):
    fits_loc = '{}/example_a.fits'.format(TEST_EXAMPLES_DIR)
    assert os.path.exists(fits_loc)
    make_png_from_fits(fits_loc, png_loc)
    assert os.path.exists(png_loc)


def test_download_images_creates_fits_and_png(new_galaxy, fits_dir, png_dir):
    data_release = '5'
    assert not os.path.exists(new_galaxy['fits_loc'])
    assert not os.path.exists(new_galaxy['png_loc'])
    download_images(new_galaxy, data_release=data_release, overwrite_fits=False, overwrite_png=False)
    assert fits_downloaded_correctly(new_galaxy['fits_loc'])
    assert os.path.exists(new_galaxy['png_loc'])


def test_check_images_are_downloaded(nsa_decals):
    checked_catalog = check_images_are_downloaded(nsa_decals)  # directly on test_examples directory
    print(checked_catalog[['iauname', 'fits_ready', 'fits_filled', 'png_ready']])
    assert np.array_equal(checked_catalog['fits_ready'], np.array([False, False, True, True, True, True]))
    assert np.array_equal(checked_catalog['fits_filled'], np.array([False, False, False, True, True, True]))
    assert np.array_equal(checked_catalog['png_ready'], np.array([False, False, True, True, True, True]))


def test_download_images_multithreaded(nsa_decals, fits_dir, png_dir):
    data_release = '5'

    # download to new temporary directory
    nsa_decals['fits_loc'] = [get_fits_loc(fits_dir, galaxy) for galaxy in nsa_decals]
    nsa_decals['png_loc'] = [get_fits_loc(png_dir, galaxy) for galaxy in nsa_decals]

    print(nsa_decals['fits_loc'])

    # verify files are not already downloaded
    for galaxy in nsa_decals:
        assert not os.path.exists(galaxy['fits_loc'])
        assert not os.path.exists(galaxy['png_loc'])

    # run download
    output_catalog = download_images_multithreaded(nsa_decals, data_release, fits_dir, png_dir, overwrite_fits=False,
                                                   overwrite_png=False)

    # verify files are now downloaded (to some unknown quality)
    for galaxy in output_catalog:
        assert os.path.exists(galaxy['fits_loc'])
        assert os.path.exists(galaxy['png_loc'])

    # verify catalog correctly reports as now downloaded
    np.array_equal(output_catalog['fits_ready'], np.ones(len(output_catalog), dtype=bool))
    np.array_equal(output_catalog['fits_filled'], np.ones(len(output_catalog), dtype=bool))
    np.array_equal(output_catalog['png_ready'], np.ones(len(output_catalog), dtype=bool))
