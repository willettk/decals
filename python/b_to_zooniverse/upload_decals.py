import pandas as pd
from astropy.io import fits
from astropy.table import Table

from a_download_decals.get_catalogs.get_joint_nsa_decals_catalog import get_nsa_catalog
import b_to_zooniverse.to_zooniverse_settings as settings
from b_to_zooniverse.do_upload import upload_subject_set
from b_to_zooniverse.previous_subjects.previous_decals_subjects import get_previous_decals_subjects
from b_to_zooniverse.make_calibration_images.get_calibration_catalog import get_expert_catalog, get_expert_catalog_joined_with_decals
from b_to_zooniverse.make_calibration_images.get_calibration_images import make_calibration_images
from b_to_zooniverse.setup.check_joint_catalog import enforce_joint_catalog_columns
from shared_utilities import match_galaxies_to_catalog_table


def upload_decals_to_panoptes(joint_catalog,
                              previous_subjects,
                              expert_catalog,
                              calibration_dir,
                              new_calibration_images=False):
    """
    Using the DECALS joint catalog created by a_download_decals, upload DECALS sets to Panoptes
    Only upload new galaxies by checking against previous subjects used
    Create calibration images with different rgb conversions to check if classifications are affected
    Upload the calibration images

    Args:
        joint_catalog (astropy.Table): NSA subjects imaged by DECALS. Includes png_loc, png_ready columns.
        previous_subjects (astropy.Table):
        expert_catalog (astropy.Table): Nair 2010 (human expert) catalog of rings, bars, etc.
        calibration_dir (str): directory to save calibration images
        subject_set_name (str): name to give subject set on Panoptes. Must not already exist.
        new_calibration_images (bool): if True, remake calibration images. If False, do not.

    Returns:
        None
    """
    dr2_galaxies, dr5_only_galaxies = match_galaxies_to_catalog_table(  # unmatched galaxies are new
        galaxies=joint_catalog,
        catalog=previous_subjects,
        galaxy_suffix='',
        catalog_suffix='_dr1_2')  # if field exists in both catalogs
    print('Previously classified galaxies: {}'.format(len(dr2_galaxies)))
    print('New galaxies: {}'.format(len(dr5_only_galaxies)))

    # use Nair galaxies previously classified in DR2
    calibration_catalog = get_expert_catalog_joined_with_decals(dr2_galaxies, expert_catalog)
    print(calibration_catalog['nsa_version'])
    print(calibration_catalog['nsa_version'][0])
    # TODO move expanding rows by image earlier, to this point. Generalises better.
    # calibration_catalog_with_colours = make_calibration_images(calibration_catalog,
    #                                               calibration_dir,
    #                                               new_images=new_calibration_images)
    # calibration_set_name = 'decals_dr2_nair_calibration'
    # calibration_subjects = upload_subject_set.upload_calibration_subject_set(calibration_catalog, calibration_set_name)

    dr5_only_name = 'decals_dr5_only'
    # dr5_only_subjects = upload_subject_set.upload_galaxy_subject_set(dr5_only_galaxies, dr5_only_name)

    dr2_name = 'decals_dr2_and_nair'
    dr2_subjects = upload_subject_set.upload_galaxy_subject_set(calibration_catalog, dr2_name)

    # return main_subjects, calibration_subjects  # for debugging


if __name__ == '__main__':

    settings.new_previous_subjects = False
    settings.new_calibration_images = False

    joint_catalog = Table(fits.getdata(settings.joint_catalog_loc))
    joint_catalog = enforce_joint_catalog_columns(joint_catalog, overwrite_cache=False)

    expert_catalog = get_expert_catalog(settings.expert_catalog_loc)

    if settings.new_previous_subjects:
        raw_previous_subjects = pd.read_csv(settings.previous_subjects_loc)  # TODO make clear how to produce this
        nsa_v1_0_0 = get_nsa_catalog(settings.nsa_v1_0_0_catalog_loc, '1_0_0')  # takes a while
        previous_subjects = get_previous_decals_subjects(raw_previous_subjects, nsa_v1_0_0)
        previous_subjects.write(settings.subject_loc, overwrite=True)
    else:
        previous_subjects_df = pd.read_csv(settings.subject_loc)
        previous_subjects = Table.from_pandas(previous_subjects_df)  # previously extracted decals subjects

    upload_decals_to_panoptes(
        joint_catalog,
        previous_subjects,
        expert_catalog,
        settings.calibration_dir,
        settings.new_calibration_images)
