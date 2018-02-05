import os
import functools
from multiprocessing.dummy import Pool as ThreadPool

import matplotlib.pyplot as plt
from astropy.io import fits
from tqdm import tqdm

from get_images import image_utils
from a_download_decals.get_images.download_images_threaded import get_loc


def make_calibration_images(calibration_catalog, calibration_dir, new_images=True):
    """
    Create calibration images in two styles: as with DR2 (by Kyle Willet) and with stronger colours (new)
    The impact on classification performance will be tested
    Args:
        calibration_catalog (astropy.Table):
        calibration_dir (str): root directory to save calibration data. Subdirectories will be created
        new_images (bool): if True, create new calibration images. Else, this function does nothing.

    Returns:
        (astropy.Table) calibration_catalog with dr2_png_loc and colour_png_loc columns added
    """
    dr2_calibration_dir = '{}/dr2_style'.format(calibration_dir)
    colour_calibration_dir = '{}/colour'.format(calibration_dir)
    for subdirectory in [dr2_calibration_dir, colour_calibration_dir]:
        if not os.path.exists(subdirectory):
            os.mkdir(subdirectory)

    # make images with top 2 options selected from mosaic experiment
    calibration_catalog['dr2_png_loc'] = [get_loc(dr2_calibration_dir, galaxy, 'png') for galaxy in calibration_catalog]
    calibration_catalog['colour_png_loc'] = [get_loc(colour_calibration_dir, galaxy, 'png') for galaxy in calibration_catalog]

    # TODO extend to checking/overwriting, similar to downloader?
    if new_images:
        pbar = tqdm(total=len(calibration_catalog), unit=' image sets created')

        save_calibration_images_of_galaxy_partial = functools.partial(save_calibration_images_of_galaxy, **{'pbar': pbar})

        pool = ThreadPool(30)
        pool.map(save_calibration_images_of_galaxy_partial, calibration_catalog)
        pbar.close()
        pool.close()
        pool.join()

    return calibration_catalog


def save_calibration_images_of_galaxy(galaxy, pbar=None):

    try:
        img_data = fits.getdata(galaxy['fits_loc'])
    except:
        print('no fits at ' + galaxy['fits_loc'])
        if pbar:
            pbar.update()
        return None

    dr2_img = get_dr2_style_image(img_data)
    colour_img = get_colour_style_image(img_data)

    plt.imsave(galaxy['dr2_png_loc'], dr2_img, origin='lower')
    plt.imsave(galaxy['colour_png_loc'], colour_img, origin='lower')

    if pbar:
        pbar.update()


def get_dr2_style_image(img_data):
    img_bands = (img_data[0, :, :], img_data[1, :, :], img_data[2, :, :])
    kwargs = {
        'scales': dict(
            g=(2, 0.008),
            r=(1, 0.014),
            z=(0, 0.019)),
        'bands': 'grz',
        'mnmx': (-0.5, 300),
        'arcsinh': 1.,
        'desaturate': False
}
    return image_utils.dr2_style_rgb(img_bands, **kwargs)


def get_colour_style_image(img_data):
    img_bands = (img_data[0, :, :], img_data[1, :, :], img_data[2, :, :])
    kwargs = {
        'arcsinh': .3,
        'mn': 0,
        'mx': .4
    }
    return image_utils.lupton_rgb(img_bands, **kwargs)
