import pandas as pd
from astropy.io import fits
from astropy.table import Table

from a_download_decals.get_catalogs.get_joint_nsa_decals_catalog import get_nsa_catalog
from a_download_decals.get_images import image_utils
import b_to_zooniverse.to_zooniverse_settings as settings
from b_to_zooniverse.do_upload import upload_subject_set
from b_to_zooniverse.previous_subjects.previous_decals_subjects import get_previous_decals_subjects
from b_to_zooniverse.make_calibration_images.get_calibration_catalog import get_expert_catalog, get_expert_catalog_joined_with_decals
from b_to_zooniverse.make_calibration_images.get_calibration_images import make_catalog_png_images
from b_to_zooniverse.setup.check_joint_catalog import enforce_joint_catalog_columns
from shared_utilities import match_galaxies_to_catalog_table


def upload_decals_to_panoptes(joint_catalog,
                              previous_subjects,
                              expert_catalog,
                              calibration_dir):
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
    joint_catalog = joint_catalog[joint_catalog['png_ready'] == True]

    dr2_galaxies, dr5_only_galaxies = match_galaxies_to_catalog_table(  # unmatched galaxies are new
        galaxies=joint_catalog,
        catalog=previous_subjects,
        galaxy_suffix='',
        catalog_suffix='_dr1_2')  # if field exists in both catalogs
    print('Previously classified galaxies: {}'.format(len(dr2_galaxies)))
    print('New galaxies: {}'.format(len(dr5_only_galaxies)))

    # use Nair galaxies previously classified in DR2
    calibration_catalog = get_expert_catalog_joined_with_decals(dr2_galaxies, expert_catalog)

    # upload standard calibration set of Nair/DR2 galaxies, coloured by DR1/2 rules
    # calibration_set_name = 'decals_dr2_nair_calibration_dr2_style_250_each'
    # calibration_catalog_dr2_style = make_catalog_png_images(
    #     calibration_catalog,
    #     image_utils.get_dr2_style_image,
    #     '{}/{}'.format(calibration_dir, calibration_set_name),
    #     size=424,
    #     overwrite=overwrite_png)
    # _ = upload_subject_set.upload_nair_calibration_subject_set(
    #     calibration_catalog_dr2_style, calibration_set_name)

    # upload all Nair/DR2 galaxies, coloured by Lupton rules
    calibration_set_name = 'decals_dr2_nair_lupton_style_all'
    calibration_catalog_lupton_style = make_catalog_png_images(
        calibration_catalog,
        image_utils.get_colour_style_image,
        '{}/{}'.format(calibration_dir, calibration_set_name),
        size=424,
        overwrite=True)  # places new png in calibration folder under this name
    upload_subject_set.upload_galaxy_subject_set(
        calibration_catalog_lupton_style, calibration_set_name)

    # upload all Nair/DR2 galaxies, coloured by DR2 rules
    # calibration_set_name = 'decals_dr2_nair_dr2_style_all'
    # _ = upload_subject_set.upload_galaxy_subject_set(calibration_catalog, calibration_set_name)

    # upload first n DR2-only galaxies
    # dr2_only_name = 'first_1k_decals_dr2'
    # _ = upload_subject_set.upload_galaxy_subject_set(dr2_galaxies[:1000], dr2_only_name)

    # upload first n DR5-only galaxies
    # dr5_only_name = 'first_1k_decals_dr5_only'
    # _ = upload_subject_set.upload_galaxy_subject_set(dr5_only_galaxies[:1000], dr5_only_name)
    # dr5_only_name = 'first_10k_decals_dr5_only'
    # _ = upload_subject_set.upload_galaxy_subject_set(dr5_only_galaxies[:10000], dr5_only_name)


if __name__ == '__main__':

    settings.new_previous_subjects = False

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
    )
