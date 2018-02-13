import pandas as pd
from astropy.io import fits
from astropy.table import Table

import to_zooniverse_settings as settings
from a_download_decals.get_catalogs.get_joint_nsa_decals_catalog import get_nsa_catalog
from do_upload.create_subject_set import create_prototype_subject_set
from previous_subjects.previous_decals_subjects import get_previous_decals_subjects
from make_calibration_images.get_calibration_catalog import get_expert_catalog, get_expert_catalog_joined_with_decals
from make_calibration_images.get_calibration_images import make_calibration_images
from new_subjects.find_new_subjects import find_new_catalog_images


def upload_decals_to_panoptes(joint_catalog, previous_subjects, expert_catalog, calibration_dir, subject_set_name, new_calibration_images=False):
    """
    Using the DECALS joint catalog created by a_download_decals, upload DECALS sets to Panoptes
    Only upload new galaxies by checking against previous subjects used
    Create calibration images with different rgb conversions to check if classifications are affected

    etc

    Returns:
        None
    """

    catalog_of_new = find_new_catalog_images(new_catalog=joint_catalog, old_catalog=previous_subjects)
    print('New galaxies: {}'.format(len(catalog_of_new)))

    calibration_catalog = get_expert_catalog_joined_with_decals(joint_catalog, expert_catalog)
    calibration_catalog = make_calibration_images(calibration_catalog,
                                                  calibration_dir,
                                                  new_images=new_calibration_images)

    main_subjects, calibration_subjects = create_prototype_subject_set(catalog_of_new,
                                                                       calibration_catalog,
                                                                       subject_set_name)

    return main_subjects, calibration_subjects  # for debugging


if __name__ == '__main__':

    settings.new_previous_subjects = False
    settings.new_calibration_images = False

    joint_catalog = Table(fits.getdata(settings.joint_catalog_loc))
    expert_catalog = get_expert_catalog(settings.expert_catalog_loc)

    if settings.new_previous_subjects:
        raw_previous_subjects = pd.read_csv(settings.previous_subjects_loc)  # TODO make clear how to produce this
        nsa_v1_0_0 = get_nsa_catalog(settings.nsa_v1_0_0_catalog_loc, '1_0_0')
        previous_subjects = get_previous_decals_subjects(raw_previous_subjects, nsa_v1_0_0)
        previous_subjects.write(settings.subject_loc)
    else:
        previous_subjects_df = pd.read_csv(settings.subject_loc)
        previous_subjects = Table.from_pandas(previous_subjects_df)  # previously extracted decals subjects

    subject_set_name = 'aliasing_test'

    upload_decals_to_panoptes(
        joint_catalog,
        previous_subjects,
        expert_catalog,
        settings.calibration_dir,
        subject_set_name,
        settings.new_calibration_images)
