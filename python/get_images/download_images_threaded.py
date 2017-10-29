import os
import urllib.parse
import urllib.request
from multiprocessing.dummy import Pool as ThreadPool
import functools
import warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table
from matplotlib import pyplot as plt

from python.get_images.image_utils import dstn_rgb
from tqdm import tqdm

min_pixelscale = 0.10


def download_images_multithreaded(catalog, data_release, fits_dir, jpeg_dir, overwrite=False):
    '''
    Multi-threaded analogy to get_decals_images_and_catalogs.py downloads

    Args:
        catalog (astropy.Table): catalog of NSA galaxies and DECALS bricks
        data_release (str): DECALS data release e.g. '3'
        fits_dir (str): directory to save fits files
        jpeg_dir (str): directory to save jpeg files
        overwrite (bool): if [default=False], overwrite existing fits and jpeg. Else, skip.

    Returns:
        (list) catalog index of galaxies which timed out
        (list) catalog index of galaxies with bad pixels
    '''

    # Table does not support .apply - use list comprehension instead
    catalog['fits_loc'] = [get_fits_loc(fits_dir, catalog[index]) for index in range(len(catalog))]
    catalog['jpeg_loc'] = [get_jpeg_loc(jpeg_dir, catalog[index]) for index in range(len(catalog))]

    pbar = tqdm(total=len(catalog), unit=' images created')

    download_params = {
        'data_release': data_release,
        'fits_dir': fits_dir,
        'jpeg_dir': jpeg_dir,
        'overwrite': overwrite,
        'pbar': pbar
    }
    download_images_partial = functools.partial(download_images, **download_params)

    pool = ThreadPool(80)
    results = pool.map(download_images_partial, catalog)
    pbar.close()
    pool.close()
    pool.join()

    results = np.array(results)
    timed_out = results[:, 0]
    good_images = results[:, 1]

    print("\n{0} total galaxies processed".format(len(catalog)))
    print("{0} good images".format(sum(good_images)))
    print("{0} galaxies with bad pixels".format(len(catalog) - sum(good_images)))
    print("{0} galaxies timed out downloading data from Legacy Skyserver".format(sum(timed_out)))

    catalog['good_image'] = good_images
    catalog['timed_out'] = timed_out

    return catalog


def download_images(galaxy, fits_dir, jpeg_dir, data_release='3', overwrite=False, pbar=None, max_attempts=5):
    '''
    Download a multi-plane FITS image from the DECaLS skyserver
    Write multi-plane FITS images to separate files for each band
    Default arguments are used due to pool.map(download_images, nsa), no more args. Can fix.

    Args:
        galaxy (astropy.TableRow): catalog entry of galaxy to download
        fits_dir (str): directory to save downloaded FITS
        data_release (str): DECALS data release e.g. '2'
        remove_multi_fits (bool): if True, delete multi-band FITS after writing to single-band

    Returns:
        (bool) Download timed out?
        (bool) Image has no bad pixels?
    '''
    timed_out = False
    good_image = False

    pixscale = max(min(galaxy['PETROTH50'] * 0.04, galaxy['PETROTH90'] * 0.02), min_pixelscale)

    # Download multi-band fits images
    # TODO refactor so timed-out -> downloaded
    fits_loc = get_fits_loc(fits_dir, galaxy)
    if os.path.exists(fits_loc) is False or overwrite is True:
        attempt = 0
        downloaded = False
        while attempt < max_attempts:
            try:
                download_fits_cutout(fits_loc, data_release, galaxy['RA'], galaxy['DEC'], pixscale, 424)
                downloaded = True
                break
            except Exception as err:
                print(err, 'on galaxy {}, attempt {}'.format(galaxy['IAUNAME'], attempt))
                attempt += 1

        if not downloaded:
            warnings.warn('Failed to download {} after three attempts'.format(galaxy['IAUNAME']))
            timed_out = True
            good_image = False
            return timed_out, good_image

    # Create artistic jpeg for Galaxy Zoo
    jpeg_loc = get_jpeg_loc(jpeg_dir, galaxy)
    if os.path.exists(jpeg_loc) is False or overwrite is True:
        good_image = make_jpeg_from_fits(fits_loc, jpeg_loc)

    if pbar:
        pbar.update()

    return timed_out, good_image


def get_fits_loc(fits_dir, galaxy):
    '''
    Get full path where fits file of galaxy should be saved, given directory
    Defines standard naming convention for fits images

    Args:
        fits_dir (str): target directory for fits files
        galaxy (astropy.TableRecord): row of NSA/DECALS catalog with galaxy info e.g. name

    Returns:
        (str) full path of where galaxy fits should be saved
    '''
    return '{0}/{1}.fits'.format(fits_dir, galaxy['IAUNAME'])


def get_jpeg_loc(jpeg_dir, galaxy):
    '''
    Get full path where jpeg file of galaxy should be saved, given directory
    Defines standard naming convention for jpeg images

    Args:
        jpeg_dir (str): target directory for jpeg files
        galaxy (astropy.TableRecord): row of NSA/DECALS catalog with galaxy info e.g. name

    Returns:
        (str) full path of where galaxy jpeg should be saved
    '''
    return '{0}/{1}.jpeg'.format(jpeg_dir, galaxy['IAUNAME'])


def download_fits_cutout(download_loc, data_release, ra=114.5970, dec=21.5681, pixscale=0.262, size=424):
    '''
    Retrieve fits image from DECALS server and save to disk

    Args:
        download_loc (str): location to save file, excluding type e.g. /data/fits/test_image.fits
        ra (float): right ascension (corner? center?)
        dec (float): declination (corner? center?)
        pixscale (float): proportional to decals pixels vs. image pixels. 0.262 for 1-1 map.
        size (int): image edge length in pixels

    Returns:
        None
    '''
    params = urllib.parse.urlencode({
        'ra': ra,
        'dec': dec,
        'pixscale': pixscale,
        'size': size})
    if data_release == '1':
        url = "http://imagine.legacysurvey.org/fits-cutout-decals-dr1?{0}".format(params)
    elif data_release == '2':
        url = "http://legacysurvey.org/viewer/fits-cutout-decals-dr2?{0}".format(params)
    elif data_release == '3':
        url = 'http://legacysurvey.org/viewer/fits-cutout-decals-dr3?{0}'.format(params)
    else:
        raise ValueError('Data release "{}" not recognised'.format(data_release))
    urllib.request.urlretrieve(url, download_loc)


def make_jpeg_from_fits(fits_loc, jpeg_loc):
    '''
    Create artistically-scaled JPG from multi-band FITS using Dustin Lang's method

    Args:
        fits_loc (str): location of FITS to read
        jpegpath (str): directory to save JPG in

    Returns:
        (bool) Image successfully created with < 20% bad pixels in any band?
    '''

    # Set parameters for RGB image creation
    _scales = dict(
        g=(2, 0.008),
        r=(1, 0.014),
        z=(0, 0.019))
    _mnmx = (-0.5, 300)

    try:
        img, hdr = fits.getdata(fits_loc, 0, header=True)
    except Exception:
        print('Invalid fits at {}'.format(fits_loc))
        return False

    rgbimg = dstn_rgb(
        (img[0, :, :], img[1, :, :], img[2, :, :]),
        'grz',
        mnmx=_mnmx,
        arcsinh=1.,
        scales=_scales,
        desaturate=True)
    plt.imsave(jpeg_loc, rgbimg, origin='lower')

    badmax = 0.
    for j in range(img.shape[0]):
        band = img[j, :, :]
        nbad = (band == 0.).sum() + np.isnan(band).sum()  # count of bad pixels in band
        fracbad = nbad / np.prod(band.shape)  # fraction of bad pixels in band
        badmax = max(badmax, fracbad)  # update worst band fraction

    if badmax < 0.2:  # if worst fraction of bad pixels is < 0.2, consider image as 'good'
        return True
    else:
        return False


if __name__ == '__main__':
    data_release = '3'
    # nsa_version = '1_0_0'
    nsa_version = '0_1_2'
    fits_dir = '../../fits/nsa/dr3'
    jpeg_dir = '../../jpeg/dr3'
    nsa_decals = Table(fits.getdata(
        '/data/galaxy_zoo/decals/catalogs/nsa_v{0}_decals_dr{1}.fits'.format(nsa_version, data_release), 1))
    download_images_multithreaded(nsa_decals, data_release, fits_dir, jpeg_dir, overwrite=True)
