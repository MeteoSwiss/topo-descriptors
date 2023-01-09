import numpy as np

import topo_descriptors.helpers as hlp


def test_round_up_to_odd():
    inputs = np.arange(0.1, 10, 0.7)
    outputs = hlp.round_up_to_odd(inputs)
    expected = [1, 1, 1, 3, 3, 3, 5, 5, 5, 7, 7, 7, 9, 9, 9]
    assert outputs.dtype == np.int64
    assert all([a == b for a, b in zip(outputs, expected)])
