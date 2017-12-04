from multiprocessing.dummy import Pool as ThreadPool
import functools

from astropy.io import fits
import matplotlib.pyplot as plt
from tqdm import tqdm

from make_calibration_images.get_calibration_catalog import get_expert_catalog_joined_with_decals
from get_images import image_utils, download_images_threaded


def make_calibration_images(calibration_catalog, calibration_dir, new_images=True):

    # make images with top 2 options selected from mosaic experiment
    calibration_catalog['dr2_jpeg_loc'] = [get_dr2_jpeg_loc(calibration_dir, galaxy) for galaxy in calibration_catalog]
    calibration_catalog['colour_jpeg_loc'] = [get_colour_jpeg_loc(calibration_dir, galaxy) for galaxy in calibration_catalog]

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
    # TODO temporary fix, need to rethink the image location labelling to happen outside actual download run
    fits_dir = '/Volumes/external/decals/fits/dr5'
    fits_loc = download_images_threaded.get_fits_loc(fits_dir, galaxy)
    try:
        img_data = fits.getdata(fits_loc)
    except:
        print('no fits at ' + fits_loc)
        if pbar:
            pbar.update()
        return None

    dr2_img = get_dr2_style_image(img_data)
    colour_img = get_colour_style_image(img_data)

    plt.imsave(galaxy['dr2_jpeg_loc'], dr2_img, origin='lower')
    plt.imsave(galaxy['colour_jpeg_loc'], colour_img, origin='lower')

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
    return image_utils.dstn_rgb(img_bands, **kwargs)


def get_colour_style_image(img_data):
    img_bands = (img_data[0, :, :], img_data[1, :, :], img_data[2, :, :])
    kwargs = {
        'arcsinh': .3,
        'mn': 0,
        'mx': .4
    }
    return image_utils.lupton_rgb(img_bands, **kwargs)


def get_dr2_jpeg_loc(jpeg_dir, galaxy):
    return '{}/{}_dr2.jpeg'.format(jpeg_dir, galaxy['iauname'])


def get_colour_jpeg_loc(jpeg_dir, galaxy):
    return '{}/{}_colour.jpeg'.format(jpeg_dir, galaxy['iauname'])

if __name__ == '__main__':
    # TODO should use the catalog including metadata not raw joint
    calibration_catalog = get_expert_catalog_joined_with_decals()
    print(calibration_catalog.columns.values)
    make_calibration_images(calibration_catalog, '/Volumes/external/decals/jpeg/calibration')