

def get_bricks_loc(catalog_dir, data_release):
    """
    Returns:
        (str) location of bricks catalog, including center and exposure counts.
            Input file for DR1, DR2. Calculated from brick_coordinates_loc and brick_exposures_loc for DR5.
    """
    if data_release == '5' or data_release == '3':
        bricks_filename = 'survey-bricks-dr{}-with-coordinates.fits'.format(data_release)
    elif data_release == '2':
        bricks_filename = 'decals-bricks-dr2.fits'
    elif data_release == '1':
        bricks_filename = 'decals-bricks-dr1.fits'
    else:
        raise ValueError('Data Release "{}" not recognised'.format(data_release))
    return '{}/{}'.format(catalog_dir, bricks_filename)


nsa_version = '1_0_0'
data_release = '5'
catalog_dir = '/media/mike/EXTERNAL/decals/catalogs'
# subject_loc = '/data/galaxy_zoo/decals/subjects/decals_dr1_and_dr2.csv'

# fits_dir = '/data/temp'
# png_dir = '/data/temp'
fits_dir = '/media/mike/EXTERNAL/decals/fits/dr{}'.format(data_release)
png_dir = '/media/mike/EXTERNAL/decals/png/dr{}'.format(data_release)
# fits_dir = '/media/mike/EXTERNAL/decals/temp'
# png_dir = '/media/mike/EXTERNAL/decals/temp'

# only needed for dr3+
brick_coordinates_loc = '{}/survey-bricks.fits'.format(catalog_dir)
brick_exposures_loc = '{}/survey-bricks-dr5.fits'.format(catalog_dir)

# if dr3+, brick coordinate-exposure merged catalog will be placed/expected here
bricks_loc = get_bricks_loc(catalog_dir, data_release)

nsa_catalog_loc = '{}/nsa_v{}.fits'.format(catalog_dir, nsa_version)

joint_catalog_loc = '{}/nsa_v{}_decals_dr{}.fits'.format(
            catalog_dir, nsa_version, data_release)

upload_catalog_loc = '{}/dr{}_nsa{}_to_upload.fits'.format(catalog_dir, data_release, nsa_version)
