'''
This modules contains routines to prepare the DEM data before computing the
topographical descriptors.
'''

import numpy as np
import os
import pyproj
import rasterio
from scipy.interpolate import griddata
from tempfile import TemporaryDirectory
import time
from zipfile import ZipFile


def generate_grid(resolution, domain, crs, buffer=0):
    '''
    Generate destination grid based on target resolution, domain, and
    coordinate system.

    Parameters
    ----------
    resolution: int
        Pixel size of the new grid, in the units defined by the *crs* argument.
    domain: str
        Name of the grid domain. Currently the following domaines are
        implemented:
        +-------------+-------------------------------------------------------+
        | Name        | CRS        | Corners (WN, EN, WS, ES)                 |
        +=============+=======================================================+
        | CCS4        | EPSG:21781 | ((255000, 480000), (965000, 480000),     |
        |             |            |  (255000, -160000), (965000, -160000))   |
        +-------------+-------------------------------------------------------+
        | COSMO       | EPSG:4326  |  (0.16, 49.52), (16.75, 49.73),          |
        |             |            |  (1.33, 42.67), (15.94, 42.85))          |
        +-------------+-------------------------------------------------------+
        | Switzerland | EPSG:21781 | ((480000, 302000), (865000, 302000),     |
        |             |            |  (480000, 74000), (865000, -74000))      |
        +-------------+-------------------------------------------------------+
    crs: str
        pyproj-compatible CRS name defining the coordinate system of the new grid.
    buffer: int, optional
        An optional buffer to be added to the grid extent as defined by the
        *domain* argument, in the same units defined by the *crs* argument.

    Returns
    -------
    coords_x, coords_y : ndarray
        Two 2-d arrays representing the coordinates of the pixels centers of the
        new grid.

    Example
    -------
    >>> resolution = 1000
    >>> domain = 'Switzerland'
    >>> crs = 'EPSG:21781'
    >>> coords_x, coords_y = generate_grid(resolution, domain, crs)
    Grid "Switzerland" extent: 480000, 865000, 74000, 302000
    Grid "Switzerland" shape: 228, 385
    >>> coords_x
    array([[480500., 481500., 482500., ..., 862500., 863500., 864500.],
       [480500., 481500., 482500., ..., 862500., 863500., 864500.],
       [480500., 481500., 482500., ..., 862500., 863500., 864500.],
       ...,
       [480500., 481500., 482500., ..., 862500., 863500., 864500.],
       [480500., 481500., 482500., ..., 862500., 863500., 864500.],
       [480500., 481500., 482500., ..., 862500., 863500., 864500.]],
      dtype=float32)
    >>> coords_x.shape
    (228, 385)
    '''

    destination_proj = pyproj.Proj(init=crs)

    # pick a domain
    if domain == 'CCS4':

        domain_corners = { # chx chy
            'NW' : [255000., 480000.],
            'NE' : [965000., 480000.],
            'SW' : [255000., -160000.],
            'SE' : [965000., -160000.]
            }
        domain_proj = pyproj.Proj(init='EPSG:21781') # LV03

    elif domain == 'COSMO':

        domain_corners = { # lon lat
            'NW' : [0.16, 49.52],
            'NE' : [16.75, 49.73],
            'SW' : [1.33, 42.67],
            'SE' : [15.94, 42.85]
            }
        domain_proj = pyproj.Proj(init='EPSG:4326') # WSG84

    elif domain == 'Switzerland':

        domain_corners = { # chx chy
            # original bounds from Swisstopo DHM25 extent:
            # 'NW' : [479975., 302000.],
            # 'NE' : [865000., 302000.],
            # 'SW' : [479975., 73975.],
            # 'SE' : [865000., 73975.]
            # shifted by 25 meters:
            'NW' : [480000., 302000.],
            'NE' : [865000., 302000.],
            'SW' : [480000., 74000.],
            'SE' : [865000., 74000.]
            }
        domain_proj = pyproj.Proj(init='EPSG:21781') # LV03

    else:
        raise ValueError(f'unknown domain {domain}')

    # project domain corners onto target crs
    corners = {}
    for corner in domain_corners.keys():
        corners[corner] = pyproj.transform(domain_proj, destination_proj,
                            domain_corners[corner][0],
                            domain_corners[corner][1]
                            )

    # compute extent (with buffer) assuming the 4 corners describe a regular grid
    extent = [
        int(min(corners['NW'][0], corners['SW'][0]) - buffer), # left
        int(max(corners['NE'][0], corners['SE'][0]) + buffer), # right
        int(min(corners['SW'][1], corners['SE'][1]) - buffer), # bottom
        int(max(corners['NW'][1], corners['NE'][1]) + buffer)  # top
    ]
    print(f'Grid "{domain}" extent: {extent[0]}, {extent[1]}, {extent[2]}, {extent[3]}')

    # shape of destination grid
    ncols = int((extent[1] - extent[0]) / resolution)
    nrows = int((extent[3] - extent[2]) / resolution)
    shape = (nrows, ncols)
    print(f'Grid "{domain}" shape: {shape[0]}, {shape[1]}')

    # grid coordinates at the centre of the pixel
    coords_x = np.linspace(extent[0], extent[1] - resolution, shape[1]) + resolution / 2.
    coords_y = np.linspace(extent[2], extent[3] - resolution, shape[0]) + resolution / 2.

    # meshgrid as float32 TODO: is this a good idea?
    coords_x, coords_y = np.meshgrid(
        coords_x.astype(np.float32),
        coords_y.astype(np.float32)
        )

    return coords_x, coords_y


def reproject(
        source,
        src_grid,
        dst_grid,
        src_crs,
        dst_crs,
        ref_crs='EPSG:4326',
        resampling='nearest',
        nchunks=10,
    ):
    '''
    Reproject and resample a source raster on a destination grid.

    Parameters
    ----------
    source : ndarray
        The input 2-d array to be reprojected. All values are required to be finite.
    src_grid : list
        A two-element list containing the x and y coordinates of the source grid
        as ndarrays of the same shape as *source*.
    dst_grid : list
        A two-element list containing the x and y coordinates of the destination
        grid as 2-d arrays.
    src_crs : str
        A pyproj-compatible CRS name defining the coordinate system of the
        source grid.
    dst_crs : str
        A pyproj-compatible CRS name defining the coordinate system of the
        destination grid.
    ref_crs : str, optional
        An optional pyproj-compatible CRS name defining the coordinate system
        of the reference grid which is returned together with the destination
        grid. Default is 'EPSG:4326' (WGS84).
    resampling : str, optional
        Name of the interpolation routine to be used for the resampling. See the
        documentation for the *method* argument to the scipy.interpolate.griddata
        function. Default is 'nearest'.
    nchunks : int, optional
        Split the destination grid in *nchunks* x *nchunks* pieces to limit the
        memory usage during the resampling. Default is 10.

    Returns
    -------
    out : ndarray
        Output 2-d array with the same shape as the elements of the *dst_grid*
        argument.
    ref_grid : list
        A two-element list containing the x and y coordinates of the destination
        grid as ndarrays of the same shape as *out* in the coordinates system
        specified by *ref_crs*.
    '''

    if np.any(~np.isfinite(source)):
        raise ValueError("source contains non-finite values")

    # initialize CRSs
    src_proj = pyproj.Proj(init=src_crs)
    dst_proj = pyproj.Proj(init=dst_crs)
    ref_proj = pyproj.Proj(init=ref_crs)

    # flatten source arrays
    src_x = np.array(src_grid[0]).flatten()
    src_y = np.array(src_grid[1]).flatten()
    source = np.array(source).flatten()

    # grid resolution (approximately)
    src_res = np.median(np.abs(np.diff(src_x[:100])))

    # resampling routine always returns floats
    src_dtype = source.dtype

    # init
    grid = []
    lon = []
    lat = []

    # split destination grid in n chunks along y axis first
    print('Resampling grid', end='', flush=True)
    start = time.time()
    subygrids = np.array_split(np.stack(dst_grid), nchunks, 1)
    subygrids = [x for x in subygrids if x.size > 0]
    for i, subygrid in enumerate(subygrids):
        # then split x axis, too
        subgrids = np.array_split(subygrid, nchunks, 2)
        subgrids = [x for x in subgrids if x.size > 0]
        tiles = None
        for j, subgrid in enumerate(subgrids):

            # define corners of destination sub grid
            # TODO: this is wrt pixel centres!
            subdst_corners = { # x, y
                'NW' : [subgrid[0].min(), subgrid[1].max()],
                'NE' : [subgrid[0].max(), subgrid[1].max()],
                'SW' : [subgrid[0].min(), subgrid[1].min()],
                'SE' : [subgrid[0].max(), subgrid[1].min()],
            }

            # project corners to src crs
            subsrc_corners = {}
            for corner in subdst_corners.keys():
                subsrc_corners[corner] = pyproj.transform(dst_proj, src_proj,
                                        subdst_corners[corner][0],
                                        subdst_corners[corner][1]
                                        )

            # dirty subset of source data to keep only relevant points
            idin = _dirty_subset(src_x, src_y, subsrc_corners, 2 * src_res)

            # reproject subset of source points onto destination crs
            src_points = pyproj.transform(src_proj, dst_proj,
                            src_x[idin], src_y[idin])

            # resample source points on destination subgrid
            tile = griddata(
                    src_points,
                    source[idin],
                    (subgrid[0], subgrid[1]),
                    method=resampling,
                    )

            # reproject subset of destination grid onto reference crs
            ref_points = pyproj.transform(dst_proj, ref_proj,
                            subgrid[0].flatten(), subgrid[1].flatten())
            ref_lon = ref_points[0].reshape(tile.shape)
            ref_lat = ref_points[1].reshape(tile.shape)

            if tiles is None:
                tiles = tile.astype(src_dtype)
                ref_lons = ref_lon
                ref_lats = ref_lat
            else:
                tiles = np.concatenate((tiles, tile.astype(src_dtype)), axis=1)
                ref_lons = np.concatenate((ref_lons, ref_lon), axis=1)
                ref_lats = np.concatenate((ref_lats, ref_lat), axis=1)

        print('.', end='', flush=True)

        grid.append(tiles)
        lon.append(ref_lons)
        lat.append(ref_lats)

    grid = np.concatenate(grid, axis=0)
    ref_lon = np.concatenate(lon, axis=0)
    ref_lat = np.concatenate(lat, axis=0)
    ref_coords = (ref_lon, ref_lat)

    print(f' done! Elapsed time: {(time.time() - start) / 60:.1f} minutes')

    return grid, ref_coords


def srtm_index_tiles(tiles):
    '''
    Parse SRTM tile names, and return all tile indices included in the bounding
    box defined by the input tiles.

    Parameters
    ----------
    tiles : list
        List of SRTM tile names.

    Returns
    -------
    idcol, idrow: list
        Two lists of integers including the column and row indices of
        all SRTM tiles included in bounding box defined by the input tiles.

    Example
    -------
    >>> tiles = ['srtm_37_03', 'srtm_39_02']
    >>> srtm_index_tiles(tiles)
    ([37, 38, 39], [2, 3])
    '''

    idcol = [int(tile.split('_')[1]) for tile in tiles]
    idrow = [int(tile.split('_')[2]) for tile in tiles]
    idcol = list(range(min(idcol), max(idcol) + 1))
    idrow = list(range(min(idrow), max(idrow) + 1))

    return idcol, idrow


def srtm_import_zip(dir, tile):
    '''
    Read in a zipped SRTM tile in geotiff format.

    Parameters
    ----------
    dir: string
        Path to the directory containing SRTM zip files.
    tile : string
        Name of the SRTM tile.

    Returns
    -------
    out: ndarray
        2-dimensional field of SRTM elevation data in meters above mean sea
        level. Missing values are assigned with -32768.
    extent: list
        List of four floats, indicating the extent of the output field.
        The extent is specified with the [left, right, bottom, top] limits.
    res: tuple
        2-element tuple of float indicating the pixel size in the y and x
        directions.

    Example
    -------
    >>> dir = '/repos/repos_data/climate/cmsaf/topo/SRTM_CGIAR/GeoTiff'
    >>> tile = 'srtm_37_03'
    >>> out, extent, res = srtm_import_zip(dir, tile)
    >>> out
    array([[-32768, -32768, -32768, ...,    435,    435,    436],
       [-32768, -32768, -32768, ...,    422,    421,    422],
       [-32768, -32768, -32768, ...,    416,    417,    418],
       ...,
       [    33,     35,     37, ...,    173,    173,    174],
       [    35,     36,     35, ...,    173,    173,    173],
       [    35,     38,     37, ...,    175,    175,    174]], dtype=int16)
    >>> out.shape
    (6001, 6001)
    >>> extent
    [-0.0004166182681528685, 5.000416715065181, 44.99958369651887, 50.0004170298522]
    >>> res
    (0.0008333333333333334, 0.0008333333333333334)
    '''

    with TemporaryDirectory() as td:
        with ZipFile(os.path.join(dir, tile + '.zip'), 'r') as zipfile:
            zipfile.extractall(path=td)
        dataset = rasterio.open(os.path.join(td, tile + '.tif'))
    out = dataset.read(1)
    # out[out == -32768] = np.nan
    res = dataset.res
    extent = [dataset.bounds.left, dataset.bounds.right,
              dataset.bounds.bottom, dataset.bounds.top]

    return out, extent, res


def srtm_mosaic(dir, tiles):
    '''
    Produce an SRTM mosaic given a list of tile names.

    Parameters
    ----------
    dir: string
        Path to the directory containing SRTM zip files.
    tiles : list
        List of SRTM tile names.

    Returns
    -------
    out: ndarray
        2-dimensional field of SRTM elevation data in meters above mean sea
        level. Missing values are assigned with -32768.
    extent: list
        List of four floats, indicating the extent of the output field.
        The extent is specified with the [left, right, bottom, top] limits.
    res: tuple
        2-element tuple of float indicating the pixel size in the y and x
        directions.

    Example
    -------
    >>> dir = '/repos/repos_data/climate/cmsaf/topo/SRTM_CGIAR/GeoTiff'
    >>> tiles = ['srtm_37_03', 'srtm_38_02']
    out, extent, res = srtm_mosaic(dir, tiles)
    srtm_37_02
    srtm_38_02
    srtm_37_03
    srtm_38_03
    >>> out.shape
    (12001, 12001)
    >>> extent
    [-0.0004166182681528685,
     10.000416836137116,
     45.0004170298522,
     55.000417150924136]
    >>> res
    (0.0008333333333333334, 0.0008333333333333334)
    '''

    # build list of tile indices needed to build the SRTM composite
    idcol, idrow = srtm_index_tiles(tiles)

    srtm_fn = []
    rasters = []
    res = []
    extent = [9999, -9999, 9999, -9999]
    for row in idrow:
        raster = None
        for col in idcol:
            tilename = f'srtm_{col:02d}_{row:02d}'
            print(tilename)
            tile, bounds, tileres = srtm_import_zip(dir, tilename)

            if raster is not None and (prvstile[:, -1] == tile[:, 0]).all():
                # the outermost column of a tile overlaps with the
                # corresponding column of the next tile
                tile = tile[:, 1:]
                bounds[0] += tileres[1]

            if raster is None:
                # initialize row mosaic
                raster = tile
                res = tileres
            else:
                # check bounds (with some tolerance)
                assert prvsbounds[1] // 1e-6 * 1e-6 == bounds[0] // 1e-6 * 1e-6
                assert res[0] == tileres[0]
                assert res[1] == tileres[1]
                # merge tiles
                raster = np.concatenate((raster, tile), axis=1)

            # merge extents (left, right, bottom, top)
            extent = [max(x[0], x[1]) if i % 2 else min(x[0], x[1])
                        for i,x in enumerate(zip(extent, bounds))]

            prvstile = tile.copy()
            prvsbounds = bounds

        if len(rasters) > 0 and (rasters[0][-1, :] == raster[0, :]).all():
            # the outermost row of a tile overlaps with the corresponding row
            # of the next tile
            raster = raster[1:, :]
            extent[2] += tileres[1]

        rasters.append(raster)

    # merge the row mosaics and return the whole mosaic
    return np.concatenate(rasters, axis=0), extent, res


def _dirty_subset(x, y, corners, offset=0):
    '''
    Quick and dirty subset of an array of coordinates according to the
    corners of a subdomain.
    '''
    idx_x = np.logical_and(
        x < max(corners['NE'][0], corners['SE'][0]) + offset,
        x > min(corners['NW'][0], corners['SW'][0]) - offset
        )

    idx_y = np.logical_and(
        y < max(corners['NW'][1], corners['NE'][1]) + offset,
        y > min(corners['SW'][1], corners['SE'][1]) - offset
        )
    return np.logical_and(idx_x, idx_y)
