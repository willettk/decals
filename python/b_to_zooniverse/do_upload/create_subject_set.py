import astropy.table

from do_upload.manifest import create_manifest_from_joint_catalog, create_manifest_from_calibration_catalog, \
    upload_manifest_to_galaxy_zoo


def create_prototype_subject_set(catalog_of_new, calibration_catalog, subject_set_name):
    """
    This is the part that changes depending on what I want to upload
    Upload 5000 new galaxies
    Also upload from the calibration set: all bars, all rings, and 500 others
    Args:
        catalog_of_new (astropy.Table): table of galaxies not previously uploaded as subjects, with nsa and decals info
        calibration_catalog (astropy.Table): table of expertly classified galaxies, with

    Returns:

    """
    catalog_of_new = catalog_of_new[catalog_of_new['png_ready'] == True]
    calibration_catalog = calibration_catalog[calibration_catalog['png_ready'] == True]

    # get 1000 new galaxies
    new_galaxies = catalog_of_new[:50]

    image_columns = ['dr2_png_loc', 'colour_png_loc']

    # get all calibration bar and ring galaxies
    featured_galaxies = calibration_catalog[calibration_catalog['has_bar'] | calibration_catalog['has_ring']]
    other_galaxies = calibration_catalog[~calibration_catalog['has_bar'] & ~calibration_catalog['has_ring']][:500]

    calibration_galaxies = astropy.table.vstack([featured_galaxies, other_galaxies])

    new_galaxies_manifest = create_manifest_from_joint_catalog(new_galaxies)
    # calibration_galaxies_manifest = create_manifest_from_calibration_catalog(calibration_galaxies, image_columns)

    main_subject_set = upload_manifest_to_galaxy_zoo(subject_set_name, new_galaxies_manifest)
    # calibration_subject_set = upload_manifest_to_galaxy_zoo('{}_calibration'.format(subject_set_name), calibration_galaxies_manifest)
    # return main_subject_set, calibration_subject_set  # for debugging

    return main_subject_set