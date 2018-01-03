# General settings for to_zooniverse

nsa_version = '1_0_0'
data_release = '5'
subject_loc = '/data/galaxy_zoo/decals/subjects/decals_dr1_and_dr2.csv'
previous_subjects_loc = '/data/galaxy_zoo/decals/subjects/galaxy_zoo_subjects.csv'  # DR1 and DR2 subjects
nsa_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nsa_v{}.fits'.format(nsa_version)

# new: use the (Table) catalog which includes checks for which images are ready
joint_catalog_loc = '/data/galaxy_zoo/decals/catalogs/dr{}_nsa{}_to_upload.fits'.format(data_release, nsa_version)
# joint_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nsa_v{}_decals_dr5.fits'.format(nsa_version)
expert_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nair_sdss_catalog.fit'
calibration_dir = '/Volumes/external/decals/jpeg/calibration'

new_previous_subjects = False
new_calibration_images = False
