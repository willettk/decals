import pytest

from get_calibration_catalog import *


@pytest.fixture()
def expert_catalog():
    # expert catalog uses iau name at 2dp (originally 'sdss')
    return Table([
        {'iauname': 'J100000.00-000000.00',
         'iauname_1dp': 'J100000.00-000000.0',  # loading helper adds this column
         'bar': 2 ** 5,
         'ring': 2 ** 3},

        {'iauname': 'J200000.00-000000.00',
         'iauname_1dp': 'J200000.00-000000.0',  # loading helper adds this column
         'bar': 2 ** 3,
         'ring': 2 ** 1}
    ])


@pytest.fixture()
def joint_catalog():
    # joint catalog uses iau name at 1dp
    return Table([
        {'iauname': 'J100000.00-000000.0',
         'feature': 'value_1'},

        {'iauname': 'J900000.00-000000.0',
         'feature': 'value_2'}
    ])


def test_decode_binary_mask():
    decoded_ints = decode_binary_mask(2 ** 2 + 2 ** 3)
    assert decoded_ints == [2, 3]

    decoded_ints = decode_binary_mask(2 ** 4)
    assert decoded_ints == [4]

    decoded_ints = decode_binary_mask(0)
    assert decoded_ints == []


def test_decode_bar_ints():

    peanut_bar = decode_bar_ints(2 ** 5)
    assert peanut_bar == ['peanut']

    weak_and_nuclear = decode_bar_ints(2 ** 3 + 2 ** 6)
    assert weak_and_nuclear == ['weak', 'nuclear']


def test_decode_ring_ints():

    inner_ring = decode_ring_ints(2 ** 2)
    assert inner_ring == ['inner']

    nuclear_and_outer = decode_ring_ints(2 ** 1 + 2 ** 3)
    assert nuclear_and_outer == ['nuclear', 'outer']


def test_round_iauname_to_1dp():
    initial_name = 'J000000.00-000000.00'
    assert round_iauname_to_1dp(initial_name) == 'J000000.00-000000.0'


def test_get_expert_catalog_joined_with_decals(joint_catalog, expert_catalog):
    joined_catalog = get_expert_catalog_joined_with_decals(joint_catalog, expert_catalog)
    assert len(joined_catalog) == 1

    expected_catalog = Table([{
        'iauname': 'J100000.00-000000.00',
        'iauname_1dp': 'J100000.00-000000.0',
        'feature': 'value_1',
        'bar': 2 ** 5,
        'ring': 2 ** 3,
        'has_bar': True,
        'has_ring': True,
        'bar_types': ['peanut'],
        'ring_types': ['outer']
    }])

    assert joined_catalog == expected_catalog
