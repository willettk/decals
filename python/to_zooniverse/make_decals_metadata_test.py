import pytest

from make_decals_metadata import *
from astropy.units import Quantity


@pytest.fixture()
def catalog():
    return Table([
        {'nsa_id': 'example',
         'ra': 10.0,
         'dec': 12.0,
         'petroth50': 50.,
         'petrotheta': 4.,
         'petroflux': [20., 21., 22., 23., 24., 25., 26.],
         'nsa_version': '1_0_0',
         'z': 0.1,
         'mag': [0., 1., 2., 3., 4., 5., 6.],
         'absmag': [10., 11., 12., 13., 14., 15., 16.],
         'nmgy': [30., 31., 32., 33., 34., 35., 36.],
         'another_column': 'sadness'}
    ])


def test_get_r_magnitude(catalog):
    mag = get_r_magnitude(catalog)
    assert mag[0] == 14.


def test_get_extinction_corrected_magnitudes(catalog):
    corrected_mag = get_extinction_corrected_magnitudes(catalog)
    assert corrected_mag[0][0] == 22.5 - 2.5 * np.log10(30.).astype(float)


def test_get_r_flux(catalog):
    flux = get_r_flux(catalog)
    assert flux[0] == 24.


def test_calculate_size_in_kpc(catalog):
    size = calculate_sizes_in_kpc(catalog)
    assert size[0] == Quantity(7.4570048368983, u.kpc)


def test_get_key_astrophysical_columns(catalog):
    key_data = get_key_astrophysical_columns(catalog)
    expected_columns = ['ra',
                        'dec',
                        'nsa_id',
                        'petroth50',
                        'petrotheta',
                        'nsa_version',
                        'absolute_size',
                        'mag_r',
                        'petroflux',
                        'redshift']
    actual_columns = set(key_data.colnames)
    for expected_col in expected_columns:
        assert expected_col in actual_columns
    assert 'another column' not in actual_columns


def test_get_key_astrophysical_columns_warns_if_col_missing(catalog):
    del catalog['nsa_version']
    with pytest.warns(UserWarning):
        get_key_astrophysical_columns(catalog)
