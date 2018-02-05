# General settings for to_zooniverse

zooniverse_login_loc = '/data/repos/decals/python/b_to_zooniverse/do_upload/zooniverse_login.txt'

subject_loc = '/data/galaxy_zoo/decals/subjects/decals_dr1_and_dr2.csv'
previous_subjects_loc = '/data/galaxy_zoo/decals/subjects/galaxy_zoo_subjects.csv'  # DR1 and DR2 subjects
nsa_v1_0_0_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nsa_v1_0_0.fits'

nsa_version = '1_0_0'  # to select which joint catalog
data_release = '5'
joint_catalog_loc = '/data/galaxy_zoo/decals/catalogs/dr{}_nsa{}_to_upload.fits'.format(data_release, nsa_version)
expert_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nair_sdss_catalog.fit'
calibration_dir = '/Volumes/external/decals/png/calibration'

new_previous_subjects = False
new_calibration_images = True
