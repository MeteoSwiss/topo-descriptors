import numpy as np

from topo_descriptors import topo


def test_sx_distance():

    radius = 150.0
    dx = 50.0
    dy = 40.0

    output = topo._sx_distance(radius, dx, dy)
    expected_first_row = np.array(
        [
            256.1249695,
            219.31712199,
            188.67962264,
            167.63054614,
            160.0,
            167.63054614,
            188.67962264,
            219.31712199,
            256.1249695,
        ]
    )

    assert np.all(np.isclose(output[0, :], expected_first_row))
    assert output.dtype == np.float64


def test_sx_bresenhamlines():

    start = np.array([[8, 9], [17, 22]])
    end = np.array([15, 15])
    output = topo._sx_bresenhamlines(start, end)
    expected = np.array(
        [
            [9, 10],
            [10, 11],
            [11, 12],
            [12, 12],
            [13, 13],
            [14, 14],
            [17, 21],
            [16, 20],
            [16, 19],
            [16, 18],
            [16, 17],
            [15, 16],
        ]
    )

    assert np.all(output == expected)
    assert output.dtype == np.int64


def test_sx_source_idx_delta():

    azimuths = np.array([3.0, 4.0, 5.0, 6.0])
    radius = 500
    dx = 20
    dy = 30
    output = topo._sx_source_idx_delta(azimuths, radius, dx, dy)
    expected = np.array([[17, 1], [17, 2], [17, 2], [17, 3]])

    assert np.all(output == expected)
    assert output.dtype == np.int64
