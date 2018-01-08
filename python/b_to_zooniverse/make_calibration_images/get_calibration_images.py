import functools
from multiprocessing.dummy import Pool as ThreadPool

import matplotlib.pyplot as plt
from astropy.io import fits
from tqdm import tqdm

from get_images import image_utils


def make_calibration_images(calibration_catalog, calibration_dir, new_images=True):

    # make images with top 2 options selected from mosaic experiment
    calibration_catalog['dr2_png_loc'] = [get_dr2_png_loc(calibration_dir, galaxy) for galaxy in calibration_catalog]
    calibration_catalog['colour_png_loc'] = [get_colour_png_loc(calibration_dir, galaxy) for galaxy in calibration_catalog]

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


def get_dr2_png_loc(png_dir, galaxy):
    return '{}/{}_dr2.png'.format(png_dir, galaxy['iauname'])


def get_colour_png_loc(png_dir, galaxy):
    return '{}/{}_colour.png'.format(png_dir, galaxy['iauname'])
