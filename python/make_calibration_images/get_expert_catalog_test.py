import pytest

from astropy.table import Table

row_limit = 20
catalog_dir = '/data/galaxy_zoo/decals/catalogs'


@pytest.fixture()
def joint_catalog():
    # DR5 is from 0-360 in RA and -20 to +30 in Dec
    return Table([
        {'iauname': 'example_a',
         'ra': 150.,
         'dec': 3.,
         'fits_loc': 'test_examples/example_a',
         'z': 0.01,
         'petrotheta': 5.},
    ])


@pytest.fixture()
def expert_catalog():
    # DR5 is from 0-360 in RA and -20 to +30 in Dec
    return Table([
        {'iauname': 'example_a',
         'ra': 150.,
         'dec': 3.,
         'TT': 3},
    ])

