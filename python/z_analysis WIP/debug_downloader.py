
from astropy.table import Table
from astropy.io import fits

from get_images.download_images_threaded import download_images_multithreaded


def main():
    """
    Create the NSA-DECaLS-GZ catalog, download fits, produce png, and identify new subjects

    Returns:
        None
    """

    data_release = '5'
    fits_dir = '/data/temp'
    png_dir = '/data/temp'

    overwrite_fits = False
    overwrite_png = False

    joint_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nsa_v1_0_0_decals_dr5_first_1k.fits'
    # joint_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nsa_v1_0_0_decals_dr5_last_1k.fits'
    joint_catalog = Table(fits.getdata(joint_catalog_loc))

    _ = download_images_multithreaded(
        joint_catalog,
        data_release,
        fits_dir,
        png_dir,
        overwrite_fits=overwrite_fits,
        overwrite_png=overwrite_png)

if __name__ == '__main__':
    main()
