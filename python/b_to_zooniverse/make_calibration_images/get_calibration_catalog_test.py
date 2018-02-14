import pytest

from astropy.table import Table

from b_to_zooniverse.make_calibration_images.get_calibration_catalog import decode_binary_mask, decode_bar_ints,\
    decode_ring_ints, get_expert_catalog_joined_with_decals


@pytest.fixture()
def expert_catalog():
    # expert catalog uses iau name at 2dp (originally 'sdss')
    return Table([
        {'iauname': 'gz_a_2dp',
         'ra': 10.,
         'dec': 10.,
         'bar': 2 ** 5,
         'ring': 2 ** 3},

        {'iauname': 'gz_b_2dp',
         'ra': 20.,
         'dec': 20.,
         'bar': 2 ** 3,
         'ring': 2 ** 1}
    ])


@pytest.fixture()
def joint_catalog():
    # joint catalog uses iau name at 1dp
    return Table([
        {'iauname': 'gz_a_1dp',
         'ra': 10.,
         'dec': 10.,
         'feature': 'value_1'},

        {'iauname': 'gz_b_1dp',
         'ra': 90.,
         'dec': 90.,
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


def test_get_expert_catalog_joined_with_decals(joint_catalog, expert_catalog):
    joined_catalog = get_expert_catalog_joined_with_decals(joint_catalog, expert_catalog)
    assert len(joined_catalog) == 1

    expected_catalog = Table([{
        'iauname': 'gz_a_1dp',
        'iauname_expert': 'gz_a_2dp',
        'feature': 'value_1',
        'ra': 10.,
        'ra_expert': 10.,
        'dec': 10.,
        'dec_expert': 10.,
        'bar': 2 ** 5,
        'ring': 2 ** 3,
        'has_bar': True,
        'has_ring': True,
        'bar_types': ['peanut'],
        'ring_types': ['outer'],
        'best_match': 0,
        'sky_separation': 0.0
    }])

    for col in set(expected_catalog.colnames + joined_catalog.colnames):
        # print(col)
        # print(joined_catalog[0][col])
        # print(expected_catalog[0][col])
        assert joined_catalog[col] == expected_catalog[col]  # more descriptive than assert Table equal
