import pytest

from get_calibration_catalog import decode_binary_mask


def test_decode_binary_mask():
    decoded_ints = decode_binary_mask(2 ** 2 + 2 ** 3)
    assert decoded_ints == [2, 3]

    decoded_ints = decode_binary_mask(2 ** 4)
    assert decoded_ints == [4]

    decoded_ints = decode_binary_mask(0)
    assert decoded_ints == []
