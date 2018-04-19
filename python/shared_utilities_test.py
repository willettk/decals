import pytest

from astropy import units as u
from astropy.table import Table

from shared_utilities import match_galaxies_to_catalog_table, fits_are_identical


@pytest.fixture()
def galaxies():
    return Table([
        {
            'name': 'a',
            'ra': 10.,
            'dec': 10.,
            'z': 0.05,
            'galaxy_data': 14.
        },

        {
            'name': 'b',
            'ra': 20.,
            'dec': 10.,
            'z': 0.05,
            'galaxy_data': 14.
        }
    ])


@pytest.fixture()
def catalog():
    return Table([
        {
            'name': 'a',
            'ra': 10.,
            'dec': 10.,
            'z': 0.05,
            'table_data': 12.
        },

        {
            'name': 'c',
            'ra': 100.,
            'dec': 80.,
            'z': 0.05,
            'table_data': 12.
        },
    ])


def test_match_galaxies_to_catalog_table(galaxies, catalog):

    matched, unmatched = match_galaxies_to_catalog_table(galaxies, catalog)

    assert matched['name'] == ['a']
    assert unmatched['name'] == ['b']

    assert set(matched.colnames) == {'dec_subject',  'galaxy_data', 'name_subject', 'ra_subject', 'z_subject', 'best_match', 'sky_separation', 'dec', 'name', 'ra', 'table_data', 'z'}
    assert set(unmatched.colnames) == {'dec', 'name', 'ra', 'z', 'best_match', 'sky_separation', 'galaxy_data'}


def test_match_galaxies_to_catalog_table_awkward_units(galaxies, catalog):
    galaxies['ra'] = galaxies['ra'] * u.deg
    catalog['dec'] = catalog['dec'] * u.deg

    matched, unmatched = match_galaxies_to_catalog_table(galaxies, catalog)

    assert matched['name'] == ['a']
    assert unmatched['name'] == ['b']

    assert set(matched.colnames) == {'dec_subject',  'galaxy_data', 'name_subject', 'ra_subject', 'z_subject', 'best_match', 'sky_separation', 'dec', 'name', 'ra', 'table_data', 'z'}
    assert set(unmatched.colnames) == {'dec', 'name', 'ra', 'z', 'best_match', 'sky_separation', 'galaxy_data'}


# @pytest.fixture()
# def nsa():
#     nsa_v1_0_0_catalog_loc = '/data/galaxy_zoo/decals/catalogs/nsa_v1_0_0.fits'
#     return Table(fits.getdata(nsa_v1_0_0_catalog_loc, 1))[:1000]
#
#
# def test_astropy_table_to_pandas(nsa):
#     df = astropy_table_to_pandas(nsa)
#     assert len(df) == len(nsa)
#     assert len(df.columns.values) == len(nsa.colnames)
    # df.to_csv('temp.csv')
    #
    # df = pd.read_csv('temp.csv')
    # value = df['PROFMEAN'].iloc[0]
    # print(value)
    # print(type(value))
    # evaluated = ast.literal_eval(value)
    # print(evaluated)
    # print(type(evaluated))

TEST_EXAMPLES_DIR = 'python/test_examples'


def test_fits_are_identical():
    fits_a_loc = '{}/example_a.fits'.format(TEST_EXAMPLES_DIR)
    fits_b_loc = '{}/example_b.fits'.format(TEST_EXAMPLES_DIR)
    assert fits_are_identical(fits_a_loc, fits_a_loc)
    assert not fits_are_identical(fits_a_loc, fits_b_loc)
