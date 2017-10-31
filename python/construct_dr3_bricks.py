import astropy.table
from astropy.io import fits
from astropy.table import Table


def merge_bricks_catalogs(data_release):
    '''
    DR3 has put brick exposure counts and brick coordinate edges in different tables (sad face).
    These tables need to be merged so that bricks with both good exposures and matching coordinates can be identified

    Returns:
        None
    '''
    catalog_dir = '/data/galaxy_zoo/decals/catalogs/'
    brick_coordinate_catalog_filename = 'survey-bricks-data_release{}-edges.fits'.format(data_release)
    brick_exposures_catalog_filename = 'survey-bricks-data_release{}.fits'.format(data_release)

    coordinate_catalog = Table(fits.getdata(catalog_dir + brick_coordinate_catalog_filename, 1))
    exposure_catalog = Table(fits.getdata(catalog_dir + brick_exposures_catalog_filename, 1))

    # coordinate catalog has brick listing for whole survey - more bricks than have been imaged in DR3
    assert len(coordinate_catalog) > len(exposure_catalog)

    # Coordinate catalog has uppercase column names. Rename to lowercase match exposure_catalog.
    for colname in coordinate_catalog.colnames:
        coordinate_catalog.rename_column(colname, colname.lower())

    # merge on brickname, ra, dec by default
    bricks_catalog = astropy.table.join(
        coordinate_catalog,
        exposure_catalog,
        keys=['brickname', 'ra', 'dec'],
        join_type='inner')

    # check that all exposed bricks were 1-1 matched
    assert len(bricks_catalog) == len(exposure_catalog)

    bricks_catalog.write(catalog_dir + 'survey-bricks-data_release{}-with-coordinates.fits'.format(data_release), overwrite=True)

if __name__ == '__main__':
    data_release = '5'
    merge_bricks_catalogs(data_release)
