import astropy.table
from astropy.io import fits
from astropy.table import Table


def merge_bricks_catalogs(brick_coordinates_loc, brick_exposures_loc, brick_loc):
    """
    DR3 has put brick coordinate edges and brick exposure counts in different tables (sad face).
    These tables need to be merged so that bricks with both good exposures and matching coordinates can be identified


    Args:
        brick_coordinates_loc (str): absolute file path to fits table of coordinates of all bricks
        brick_exposures_loc (str):  absolute file path to fits table of coordinates of imaged bricks with exposure counts
        brick_loc (str): absolute file path to save fits table of brick coordinates with exposure counts

    Returns:
        None
    """

    coordinate_catalog = Table(fits.getdata(brick_coordinates_loc, 1))
    exposure_catalog = Table(fits.getdata(brick_exposures_loc, 1))

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

    bricks_catalog.write(brick_loc, overwrite=True)

if __name__ == '__main__':
    data_release = '5'
    merge_bricks_catalogs(data_release)
