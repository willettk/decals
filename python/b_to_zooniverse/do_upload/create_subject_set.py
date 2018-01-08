import astropy.table

from do_upload.manifest import create_manifest_from_joint_catalog, create_manifest_from_calibration_catalog, \
    upload_manifest_to_galaxy_zoo


def create_prototype_subject_set(catalog_of_new, calibration_catalog):
    """
    This is the part that changes depending on what I want to upload
    Upload 5000 new galaxies
    Also upload from the calibration set: all bars, all rings, and 500 others
    Args:
        catalog_of_new (?): table of galaxies not previously uploaded as subjects, with nsa and decals info
        calibration_catalog (astropy.Table): table of expertly classified galaxies, with

    Returns:

    """

    # get 1000 new galaxies
    new_galaxies = catalog_of_new[:2500]

    image_columns = ['dr2_jpeg_loc', 'colour_jpeg_loc']

    # get all calibration bar and ring galaxies
    bar_galaxies = calibration_catalog[calibration_catalog['has_bar']]
    ring_galaxies = calibration_catalog[calibration_catalog['has_ring']]
    other_galaxies = calibration_catalog[~calibration_catalog['has_bar'] & ~calibration_catalog['has_ring']][:500]

    calibration_galaxies = astropy.table.vstack([bar_galaxies, ring_galaxies, other_galaxies])
    print(calibration_galaxies.colnames)

    print('{} bars, {} rings'.format(len(bar_galaxies), len(ring_galaxies)))

    new_galaxies_manifest = create_manifest_from_joint_catalog(new_galaxies)
    calibration_galaxies_manifest = create_manifest_from_calibration_catalog(calibration_galaxies, image_columns)

    upload_manifest_to_galaxy_zoo('decals_prototype', new_galaxies_manifest)
    upload_manifest_to_galaxy_zoo('decals_calibration', calibration_galaxies_manifest)
