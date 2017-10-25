import os
import urllib.parse
import urllib.request
from multiprocessing.dummy import Pool as ThreadPool
import functools

import numpy as np
from astropy.io import fits
from astropy.table import Table
from matplotlib import pyplot as plt

from python.get_images.image_utils import dstn_rgb
from tqdm import tqdm


min_pixelscale = 0.10


def download_images_multithreaded(nsa_decals, data_release, fits_dir, jpeg_dir, overwrite=False):
    '''
    Multi-threaded analogy to decals_dr2.py downloads

    Args:
        nsa_decals (astropy.Table): catalog of NSA galaxies and DECALS bricks
        data_release (str): DECALS data release e.g. '3'
        fits_dir (str): directory to save fits files
        jpeg_dir (str): directory to save jpeg files
        overwrite (bool): if [default=False], overwrite existing fits and jpeg. Else, skip.

    Returns:
        (list) catalog index of galaxies which timed out
        (list) catalog index of galaxies with bad pixels
    '''

    pbar = tqdm(total=len(nsa_decals), unit=' images created')

    download_params = {
        'data_release': data_release,
        'fits_dir': fits_dir,
        'jpeg_dir': jpeg_dir,
        'overwrite': overwrite,
        'pbar': pbar
    }
    download_images_partial = functools.partial(download_images, **download_params)

    pool = ThreadPool(8)
    results = pool.map(download_images_partial, nsa_decals)
    pbar.close()
    pool.close()
    pool.join()

    results = np.array(results)
    timed_out = results[:, 0]
    good_images = results[:, 1]

    # TODO fix logs
    # timed_out_log_loc = "../failed_fits_downloads.log"
    # timed_out_log = open(timed_out_log_loc, 'w')
    # print(timed_out_log, "\n".join(nsa_decals['IAUNAME'][timed_out]))
    # timed_out_log.close()
    #
    # bad_pix_log_loc = "../bad_pix_images.log"
    # bad_pix_log = open(bad_pix_log_loc, 'w')
    # print(bad_pix_log, "\n".join(nsa_decals['IAUNAME'][~good_images]))
    # bad_pix_log.close()

    print("\n{0} total galaxies processed".format(len(nsa_decals)))
    print("{0} good images".format(sum(good_images)))
    print("{0} galaxies with bad pixels".format(len(nsa_decals) - sum(good_images)))
    print("{0} galaxies timed out downloading data from Legacy Skyserver".format(sum(timed_out)))

    return timed_out, good_images


def download_images(galaxy, fits_dir, jpeg_dir, data_release='3', overwrite=False, pbar=None):
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

    fits_loc = get_fits_loc(fits_dir, galaxy)
    pixscale = max(min(galaxy['PETROTH50'] * 0.04, galaxy['PETROTH90'] * 0.02), min_pixelscale)

    # Download multi-band fits images
    if os.path.exists(fits_loc) is False or overwrite is True:
        try:
            download_fits_cutout(fits_loc, data_release, galaxy['RA'], galaxy['DEC'], pixscale, 424)
        except IOError as io_err:
            print(io_err)
            exit(1)
        except Exception as err:
            print('Assuming the above error is a time-out')
            print(err)
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
    Returns:
        (bool) FITS Download timed out?
        (bool) Image has no bad pixels?
    '''
    good_image = False


    # Set parameters for RGB image creation
    _scales = dict(
        g=(2, 0.008),
        r=(1, 0.014),
        z=(0, 0.019))
    _mnmx = (-0.5, 300)

    img, hdr = fits.getdata(fits_loc, 0, header=True)

    badmax = 0.
    for j in range(img.shape[0]):
        band = img[j, :, :]
        nbad = (band == 0.).sum() + np.isnan(band).sum()
        fracbad = nbad / np.prod(band.shape)
        badmax = max(badmax, fracbad)

    if badmax < 0.2:
        rgbimg = dstn_rgb(
            (img[0, :, :], img[1, :, :], img[2, :, :]),
            'grz',
            mnmx=_mnmx,
            arcsinh=1.,
            scales=_scales,
            desaturate=True)
        plt.imsave(jpeg_loc, rgbimg, origin='lower')
        good_image = True

    return good_image


if __name__ == '__main__':
    data_release = '3'
    # nsa_version = '1_0_0'
    nsa_version = '0_1_2'
    fits_dir = '../../fits/nsa/dr3'
    jpeg_dir = '../../jpeg/dr3'
    # nsa_decals = Table(fits.getdata('../fits/nsa_v{0}_decals_dr{1}_after_cuts.fits'.format(nsa_version, dr), 1))
    nsa_decals = Table(fits.getdata('../../fits/nsa_v{0}_decals_dr{1}.fits'.format(nsa_version, data_release), 1))
    download_images_multithreaded(nsa_decals, data_release, fits_dir, jpeg_dir, overwrite=True)
