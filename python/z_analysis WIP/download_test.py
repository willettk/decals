from astropy.io import fits
from astropy.table import Table
from PIL import Image
import numpy as np

import urllib.parse
import urllib.request

from smoothing import get_lupton


def get_comparison_images(name, ra, dec):
    target_size = 424
    target_pixscale = 0.1
    native_pixscale = 0.262
    download_size = np.ceil(target_size * target_pixscale / native_pixscale).astype(int)

    native_fits_loc = '/Data/repos/decals/python/test_examples/resizing/{}_native_internet.fits'.format(name)
    params = urllib.parse.urlencode({
            'ra': ra,
            'dec': dec,
            'size': target_size,
            'layer': 'decals-dr5'})
    url = "http://legacysurvey.org/viewer/fits-cutout?{}".format(params)
    urllib.request.urlretrieve(url, native_fits_loc)
    previous_image = fits.getdata(native_fits_loc)
    previous_image_bands = (previous_image[0, :, :], previous_image[1, :, :], previous_image[2, :, :])

    previous_image = get_lupton(previous_image_bands, size=424)
    previous_pil_image = Image.fromarray(np.uint8(previous_image * 255.), mode='RGB')
    previous_pil_image.save('/Data/repos/decals/python/test_examples/resizing/J00001_previous.png')

    native_fits_loc = '/Data/repos/decals/python/test_examples/resizing/{}_{}px_native_internet.fits'.format(name, download_size)
    params = urllib.parse.urlencode({
            'ra': ra,
            'dec': dec,
            'size': download_size,
            'layer': 'decals-dr5'})
    url = "http://legacysurvey.org/viewer/fits-cutout?{}".format(params)
    urllib.request.urlretrieve(url, native_fits_loc)
    native_image = fits.getdata(native_fits_loc)
    native_image_bands = (native_image[0, :, :], native_image[1, :, :], native_image[2, :, :])
    native_image = get_lupton(native_image_bands, size=download_size)

    native_pil_image = Image.fromarray(np.uint8(native_image * 255.), mode='RGB')
    # see http://pillow.readthedocs.io/en/3.1.x/reference/Image.html for resize modes
    nearest_image = native_pil_image.resize(size=(target_size, target_size), resample=Image.NEAREST)
    nearest_image.save('/Data/repos/decals/python/test_examples/resizing/{}_from_native_nearest.png'.format(name))
    nearest_image = native_pil_image.resize(size=(target_size, target_size), resample=Image.BILINEAR)
    nearest_image.save('/Data/repos/decals/python/test_examples/resizing/{}_from_native_bilinear.png'.format(name))
    nearest_image = native_pil_image.resize(size=(target_size, target_size), resample=Image.BICUBIC)
    nearest_image.save('/Data/repos/decals/python/test_examples/resizing/{}_from_native_bicubic.png'.format(name))
    nearest_image = native_pil_image.resize(size=(target_size, target_size), resample=Image.LANCZOS)
    nearest_image.save('/Data/repos/decals/python/test_examples/resizing/{}_from_native_lanczos.png'.format(name))


if __name__ == '__main__':
    catalog_loc = '/data/galaxy_zoo/decals/catalogs/dr5_nsa1_0_0_to_upload.fits'
    nsa = Table.read(catalog_loc)
    for galaxy in nsa[:20]:
        get_comparison_images(name=galaxy['iauname'], ra=galaxy['ra'], dec=galaxy['dec'])