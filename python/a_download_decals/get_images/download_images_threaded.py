import functools
import os
import requests
import urllib.parse
import urllib.request
import warnings
import logging
import multiprocessing
from subprocess import call

import threading
from queue import Queue

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from tqdm import tqdm

import time

from get_images.image_utils import dr2_style_rgb

min_pixelscale = 0.10


def download_images_multithreaded(catalog, data_release, fits_dir, png_dir, overwrite_fits, overwrite_png):
    '''
    Download fits of catalog, create png images. Update catalog with fits and png locations. Run multithreaded.

    Args:
        catalog (astropy.Table): catalog of NSA galaxies and DECALS bricks
        data_release (str): DECALS data release e.g. '3'
        fits_dir (str): directory to save fits files
        png_dir (str): directory to save png files
        overwrite_fits (bool): if True, force new fits download.
        overwrite_png (bool): if True, force new png download.
    Returns:
        (astropy.Table) catalog with new columns for fits/png locations and quality checks
    '''

    catalog = catalog[-10000:]

    # Table does not support .apply - use list comprehension instead
    catalog['fits_loc'] = [get_fits_loc(fits_dir, catalog[index]) for index in range(len(catalog))]
    catalog['png_loc'] = [get_png_loc(png_dir, catalog[index]) for index in range(len(catalog))]

    download_params = {
        'data_release': data_release,
        'overwrite_fits': overwrite_fits,
        'overwrite_png': overwrite_png,
        'pbar': None
    }
    download_images_partial = functools.partial(download_images, **download_params)

    queue = Queue()
    [queue.put(galaxy) for galaxy in catalog]

    n_threads = 10
    for i in range(n_threads):
        t = threading.Thread(target=download_worker, args=(queue, download_images_partial))
        t.daemon = True  # kill once program exits
        t.start()

    wait_time = 40
    previous_num_left = queue.qsize() * 2.
    while not queue.empty():
        num_left = queue.qsize()
        rate = (previous_num_left - num_left) / wait_time
        print("Images left: {}, Rate: {}".format(num_left, rate))
        previous_num_left = num_left
        # if rate < 0.001:
        #     print('rate 0, breaking')
        #     break
        time.sleep(wait_time)

    # blocked_galaxies = [queue.get() for n in range(queue.qsize())]
    # import pandas as pd
    # pd.DataFrame(blocked_galaxies).to_csv('blocked_galaxies.csv')
    # exit(0)

    queue.join()

    catalog = check_images_are_downloaded(catalog)

    print("\n{} total galaxies".format(len(catalog)))
    print("{} fits are downloaded".format(np.sum(catalog['fits_ready'])))
    print("{} png generated".format(np.sum(catalog['png_ready'])))
    print("{} fits have many bad pixels".format(len(catalog) - np.sum(catalog['fits_filled'])))

    return catalog


def download_worker(queue, downloader):
    while True:
        galaxy = queue.get(timeout=30)
        downloader(galaxy)
        queue.task_done()


def download_images(galaxy, data_release, overwrite_fits=False, overwrite_png=False, pbar=None, max_attempts=5):
    """
    Download a multi-plane FITS image from the DECaLS skyserver
    Write multi-plane FITS images to separate files for each band
    Default arguments are used due to pool.map(download_images, nsa), no more args. Can fix.
    Args:
        galaxy (astropy.TableRow): catalog entry of galaxy to download
        data_release (str): DECALS data release e.g. '2'
        overwrite_fits (bool): download fits even if fits file already exists in target location
        overwrite_png (bool): download png even if png file already exists in target location
        pbar (tqdm): progress bar shared between processes, to be updated. If None, no progress bar will be shown.
        max_attempts (int): max number of fits download attempts per file

    Returns:
        None
    """
    print('download triggered', flush=True)

    try:

        pixscale = max(min(galaxy['petroth50'] * 0.04, galaxy['petroth50'] * 0.02), min_pixelscale)

        # For convenience
        fits_loc = galaxy['fits_loc']
        png_loc = galaxy['png_loc']

        # Download multi-band fits images
        # if not fits_downloaded_correctly(fits_loc) or overwrite_fits: TODO removed for speed
        print('checking existence {}'.format(fits_loc), flush=True)
        print('path', os.path.exists(fits_loc), flush=True)
        print('overwrite', overwrite_fits, flush=True)
        print('both', os.path.exists(fits_loc) or overwrite_fits, flush=True)
        print('png', os.path.exists(png_loc))
        if (not os.path.exists(fits_loc)) or overwrite_fits:
            print('begin downloading {}'.format(galaxy['iauname']), flush=True)
            attempt = 0
            while attempt < max_attempts:
                # print('making attempt, thread {}'.format(os.getpid()))
                try:
                    # print('doing the download, thread {}'.format(os.getpid()))
                    download_fits_cutout(fits_loc, data_release, galaxy['ra'], galaxy['dec'], pixscale, 424)
                    assert fits_downloaded_correctly(fits_loc)
                    # print('yay, downloaded, thread {}'.format(os.getpid()))
                    break
                except Exception as err:
                    print(err, 'on galaxy {}, attempt {}'.format(galaxy['iauname'], attempt), flush=True)
                    attempt += 1
                if attempt == max_attempts:
                    warnings.warn('Failed to download {} after three attempts. No fits at {}'.format(
                        galaxy['iauname'], galaxy['fits_loc']))

        # print('moving on to png, thread {}'.format(os.getpid()))
        if not os.path.exists(png_loc) or overwrite_png:
            try:
                # Create artistic png for Galaxy Zoo from the new FITS
                make_png_from_fits(fits_loc, png_loc)
            except:
                warnings.warn('Error creating png from {}'.format(fits_loc))

        if pbar:
            pbar.update()

    except Exception:
        warnings.warn('FATAL THREAD ERROR: {}'.format(Exception), flush=True)
        exit(1)


def fits_downloaded_correctly(fits_loc):
    """
    Is there a readable fits image at fits_loc?
    Does NOT check for bad pixels

    Args:
        fits_loc (str): location of fits file to open

    Returns:
        (bool) True if file at fits_loc is readable, else False
    """

    try:
        img, _ = fits.getdata(fits_loc, 0, header=True)
        return True
    except:  # image fails to open
        return False


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
    return '{0}/{1}.fits'.format(fits_dir, galaxy['iauname'])


def get_png_loc(png_dir, galaxy):
    '''
    Get full path where png file of galaxy should be saved, given directory
    Defines standard naming convention for png images

    Args:
        png_dir (str): target directory for png files
        galaxy (astropy.TableRecord): row of NSA/DECALS catalog with galaxy info e.g. name

    Returns:
        (str) full path of where galaxy png should be saved
    '''
    return '{0}/{1}.png'.format(png_dir, galaxy['iauname'])


def download_fits_cutout(fits_loc, data_release, ra=114.5970, dec=21.5681, pixscale=0.262, size=424):
    '''
    Retrieve fits image from DECALS server and save to disk

    Args:
        fits_loc (str): location to save file, excluding type e.g. /data/fits/test_image.fits
        ra (float): right ascension (center)
        dec (float): declination (center)
        pixscale (float): proportional to decals pixels vs. image pixels. 0.262 for 1-1 map.
        size (int): image edge length in pixels. Default 424 to match GZ2, but consider 512.

    Returns:
        None
    '''
    params = urllib.parse.urlencode({
        'ra': ra,
        'dec': dec,
        'pixscale': pixscale,
        'size': size,
        'layer': 'decals-dr{}'.format(data_release)})
    if data_release == '1':
        url = "http://imagine.legacysurvey.org/fits-cutout?{}".format(params)
    elif data_release == '2' or '3' or '5':
        url = "http://legacysurvey.org/viewer/fits-cutout?{}".format(params)
    else:
        raise ValueError('Data release "{}" not recognised'.format(data_release))

    print('urlib request {}'.format(url))

    # urllib.request.urlretrieve(url, fits_loc)
    # result = urllib.request.urlretrieve(url, fits_loc)
    # print(result)
    open(fits_loc, 'wb').write(requests.get(url, stream=True).content)
    # call('wget {} {}'.format(url, fits_loc), shell=True)
    print('done {}'.format(fits_loc))


def make_png_from_fits(fits_loc, png_loc):
    '''
    Create png from multi-band fits

    Args:
        fits_loc (str): location of .fits to create png from
        png_loc (str): location to save png

    Returns:
        None
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
        warnings.warn('Invalid fits at {}'.format(fits_loc))
    else:  # if no exception getting image
        rgbimg = dr2_style_rgb(
            (img[0, :, :], img[1, :, :], img[2, :, :]),
            'grz',
            mnmx=_mnmx,
            arcsinh=1.,
            scales=_scales,
            desaturate=True)
        plt.imsave(png_loc, rgbimg, origin='lower')


def check_images_are_downloaded(catalog):
    """
    Record if images are downloaded. Add 'fits_ready', 'png_ready' and 'fits_complete' columns to catalog.

    Args:
        catalog (astropy.Table): joint NSA/decals catalog

    Returns:
        (astropy.Table) catalog with image quality check columns added
    """
    catalog['fits_ready'] = np.zeros(len(catalog), dtype=bool)
    catalog['fits_filled'] = np.zeros(len(catalog), dtype=bool)
    catalog['png_ready'] = np.zeros(len(catalog), dtype=bool)

    for row_index, galaxy in tqdm(enumerate(catalog), total=len(catalog), unit=' images checked'):
        downloaded, complete = get_download_quality_of_fits(galaxy['fits_loc'])
        catalog['fits_ready'][row_index] = downloaded
        catalog['fits_filled'][row_index] = complete
        if downloaded:  # fits must be ready for png to be ready
            if os.path.exists(galaxy['png_loc']):
                catalog['png_ready'][row_index] = True
    return catalog


def get_download_quality_of_fits(fits_loc, badmax_limit=0.2):
    """
    Find if fits at fits_loc a) opens b) has few bad pixels

    Args:
        fits_loc (str): location of fits file to open
        badmax_limit(float): maximum ratio of empty pixels for image to be considered 'correct'

    Returns:
        (bool) is image downloaded?
        (bool) is empty pixel ratio below the allowed limit?
    """

    try:
        img, _ = fits.getdata(fits_loc, 0, header=True)
        complete = few_missing_pixels(img, badmax_limit)
        return True, complete
    except:  # image fails to open
        return False, False


def few_missing_pixels(img, badmax_limit):
    """
    Find if img contains few NaN pixels e.g. from incomplete imaging

    Args:
        img (np.array): multi-band (i.e. 3-dim) pixel data to check
        badmax_limit(float): maximum ratio of empty pixels for image to be considered 'correct' e.g. 0.2 (20%)

    Returns:
        (bool) True if few NaN pixels, else False

    """
    badmax = 0.
    for j in range(img.shape[0]):
        band = img[j, :, :]
        nbad = (band == 0.).sum() + np.isnan(band).sum()  # count of bad pixels in band
        fracbad = nbad / np.prod(band.shape)  # fraction of bad pixels in band
        badmax = max(badmax, fracbad)  # update worst band fraction

    # if worst fraction of bad pixels is < badmax_limit, consider image as 'good'
    if badmax < badmax_limit:
        return True
    else:
        return False
