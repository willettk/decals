import astropy.table

from b_to_zooniverse.do_upload import manifest
from shared_utilities import current_date
from a_download_decals.get_images import image_utils
from b_to_zooniverse.make_calibration_images import get_calibration_images
from b_to_zooniverse import to_zooniverse_settings


def upload_galaxy_subject_set(catalog, subject_set_name):
    catalog = catalog[catalog['png_ready'] == True]
    name = '{}_{}'.format(current_date(), subject_set_name)
    galaxy_manifest = manifest.create_manifest_from_joint_catalog(catalog)
    return manifest.upload_manifest_to_galaxy_zoo(name, galaxy_manifest)


def upload_color_calibration_subject_set(calibration_catalog, subject_set_name):
    """
    Upload from the calibration catalog up to 1000 bars or rings, and 1000 others
    Color them with both DR2 style and Lupton style
    Args:
        calibration_catalog (astropy.Table): table of expertly classified galaxies, with nsa and decals info
        subject_set_name (str): name to give uploaded subject set on Panoptes. Must not already exist.

    Returns:
        (astropy.Table): calibration subject set manifest, as uploaded to Panoptes
    """
    calibration_catalog = calibration_catalog[calibration_catalog['png_ready'] == True]
    featured_galaxies = calibration_catalog[calibration_catalog['has_bar'] | calibration_catalog['has_ring']][:1000]
    other_galaxies = calibration_catalog[~calibration_catalog['has_bar'] & ~calibration_catalog['has_ring']][:1000]
    calibration_galaxies = astropy.table.vstack([featured_galaxies, other_galaxies])
    print('Calibration bars: {}'.format(featured_galaxies['has_bar'].sum()))
    print('Calibration rings: {}'.format(featured_galaxies['has_ring'].sum()))
    print('Calibration non-bar, non-ring: {}'.format(len(other_galaxies)))

    get_calibration_images.make_catalog_png_images(calibration_catalog, image_utils.get_dr2_style_image, )
    name = '{}_{}_dr2_style'.format(current_date(), subject_set_name)
    calibration_galaxies_manifest = manifest.create_manifest_from_joint_catalog(calibration_galaxies)
    manifest.upload_manifest_to_galaxy_zoo(name, calibration_galaxies_manifest)

    calibration_galaxies_manifest = manifest.create_manifest_from_joint_catalog(calibration_galaxies)
    manifest.upload_manifest_to_galaxy_zoo(name, calibration_galaxies_manifest)
