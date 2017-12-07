import pandas as pd
from astropy.io import fits
from astropy.table import Table

from download_decals.get_catalogs.get_joint_nsa_decals_catalog import get_nsa_catalog
from make_calibration_images.get_calibration_catalog import get_expert_catalog, get_expert_catalog_joined_with_decals
from make_calibration_images.get_calibration_images import make_calibration_images
from to_zooniverse.create_subject_set import create_prototype_subject_set
from to_zooniverse.identify_new_images import get_new_images


def main():

    # Settings
    subject_loc = '/data/galaxy_zoo/decals/subjects/decals_dr1_and_dr2.csv'
    nsa_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nsa_v1_0_0.fits'
    joint_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nsa_v1_0_0_decals_dr5.fits'
    expert_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nair_sdss_catalog.fit'
    calibration_dir = '/Volumes/external/decals/jpeg/calibration'

    new_images = False

    previous_subjects = pd.read_csv(subject_loc)  # previously extracted decals subjects
    nsa_catalog = get_nsa_catalog(nsa_catalog_loc)
    joint_catalog = Table(fits.getdata(joint_catalog_loc))
    expert_catalog = get_expert_catalog(expert_catalog_loc)

    catalog_of_new = get_new_images(joint_catalog, previous_subjects, nsa_catalog)

    calibration_catalog = get_expert_catalog_joined_with_decals(joint_catalog, expert_catalog)
    calibration_catalog = make_calibration_images(calibration_catalog, calibration_dir, new_images=new_images)

    create_prototype_subject_set(catalog_of_new, calibration_catalog)

if __name__ == '__main__':
    main()
