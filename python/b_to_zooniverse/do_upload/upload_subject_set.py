import astropy.table

from b_to_zooniverse.do_upload import manifest
from shared_utilities import current_date

# TODO this test routine is no longer appropriate because uploader is changing - disable for now
# def create_prototype_subject_set(catalog_of_new, calibration_catalog, subject_set_name):
#     """
#     This is the part that changes depending on what I want to upload
#     Upload 5000 new galaxies
#     Also upload from the calibration set: all bars, all rings, and 500 others
#     Args:
#         catalog_of_new (astropy.Table): table of galaxies not previously uploaded as subjects, with nsa and decals info
#         calibration_catalog (astropy.Table): table of expertly classified galaxies, with
#
#     Returns:
#
#     """
#     catalog_of_new = catalog_of_new[catalog_of_new['png_ready'] == True]
#     calibration_catalog = calibration_catalog[calibration_catalog['png_ready'] == True]
#
#     # get 1000 new galaxies
#     new_galaxies = catalog_of_new[:1000]
#
#     image_columns = ['dr2_png_loc', 'colour_png_loc']
#
#     # get all calibration bar and ring galaxies
#     featured_galaxies = calibration_catalog[calibration_catalog['has_bar'] | calibration_catalog['has_ring']]
#     other_galaxies = calibration_catalog[~calibration_catalog['has_bar'] & ~calibration_catalog['has_ring']][:500]
#
#     calibration_galaxies = astropy.table.vstack([featured_galaxies, other_galaxies])
#
#     new_galaxies_manifest = manifest.create_manifest_from_joint_catalog(new_galaxies)
#     calibration_galaxies_manifest = manifest.create_manifest_from_calibration_catalog(calibration_galaxies, image_columns)
#
#     main_subject_set = manifest.upload_manifest_to_galaxy_zoo(subject_set_name, new_galaxies_manifest)
#     calibration_subject_set = manifest.upload_manifest_to_galaxy_zoo('{}_calibration'.format(subject_set_name), calibration_galaxies_manifest)
#     return main_subject_set, calibration_subject_set  # for debugging


def upload_galaxy_subject_set(catalog, subject_set_name):
    catalog = catalog[catalog['png_ready'] == True]
    name = '{}_{}'.format(current_date(), subject_set_name)
    galaxy_manifest = manifest.create_manifest_from_joint_catalog(catalog)
    return manifest.upload_manifest_to_galaxy_zoo(name, galaxy_manifest)


def upload_calibration_subject_set(calibration_catalog, subject_set_name):
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

    image_columns = ['dr2_png_loc', 'colour_png_loc']
    name = '{}_{}'.format(current_date(), subject_set_name)
    calibration_galaxies_manifest = manifest.create_manifest_from_calibration_catalog(calibration_galaxies, image_columns)
    return manifest.upload_manifest_to_galaxy_zoo(name, calibration_galaxies_manifest)
