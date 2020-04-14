'''
This scripts takes the SRTM mosaic and reprojects it on a regular grid, namely
the radar/INCA CCS4 grid, with Swiss LV03 coordinates, and 50m resolution.
The domain is buffered by 30 km to allow the computation of large scale
topographical descriptors.
Additionally, it combines it with higher quality DEM data from Swisstopo for
grid points within Switzerland.
'''

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from Topo_descriptors.preprocessing import generate_grid, reproject

# path to input datasets
fn_srtm = 'resources/srtm_mosaic_WSG84.nc'
fn_swisstopo = 'resources/dhm25_clean.nc'

# the first source dataset is the SRTM90m v4.1 mosaic on the native WSG84 grid
# with 0.0008333° resolution
srtm = xr.open_dataset(fn_srtm)
src_crs = 'EPSG:4326'
src_grid = np.meshgrid(srtm.lon.values, srtm.lat.values)
srtm = srtm.SRTM90.values.astype(np.int16)

# which we need to reproject on the CCS4 domain with Swiss coordinates, 50m
# resolution, and 15 km buffer to allow the computation of large scale topographical
# descriptors
resolution = 50
buffer = 15e3
domain = 'CCS4'
dst_crs = 'EPSG:21781' # LV03
dst_grid = generate_grid(resolution, domain, dst_crs, buffer)

# reproject and resample the SRTM source data on the CCS4 destination grid
destination, ref_coords = reproject(
    srtm,
    src_grid=src_grid,
    dst_grid=dst_grid,
    src_crs=src_crs,
    dst_crs=dst_crs,
    resampling='linear',
    nchunks=10,
    )

# as xarray
srtm_ccs4 = xr.Dataset(
        {'SRTM90' : (['chy', 'chx'], destination.astype(np.int16))},
        coords={
            'chx' : ('chx', dst_grid[0][0, :].astype(np.float32)),
            'chy' : ('chy', dst_grid[1][:, 0].astype(np.float32)),
            'lon' : (['chy', 'chx'], ref_coords[0].astype(np.float32)),
            'lat' : (['chy', 'chx'], ref_coords[1].astype(np.float32))
            }
        # TODO: add attributes
    )
print(srtm_ccs4)

# the second source is the Swisstopo DHM25 on a LV03 grid
swisstopo = xr.open_dataset(fn_swisstopo)
swisstopo = swisstopo.reindex(chy=swisstopo.chy[::-1])

# which we first need to upscale to 50 m
swisstopo = swisstopo.coarsen(chx=2, chy=2, boundary='trim').mean()
swisstopo = swisstopo.where(swisstopo.DHM25 > 0, -9999)
print(swisstopo)

# and then to align with the CCS4 domain
src_crs = 'EPSG:21781'
src_grid = np.meshgrid(swisstopo.chx.values, swisstopo.chy.values)
swisstopo = swisstopo.DHM25.values.astype(np.int16)
resolution = 50
buffer = 0
domain = 'Switzerland'
dst_crs = 'EPSG:21781' # LV03
dst_grid = generate_grid(resolution, domain, dst_crs, buffer)

# reproject and resample Swisstopo on new destination grid
destination, ref_coords = reproject(
    swisstopo,
    src_grid=src_grid,
    dst_grid=dst_grid,
    src_crs=src_crs,
    dst_crs=dst_crs,
    resampling='linear',
    nchunks=10,
    )
destination[destination < 0] = -9999

# as xarray
swisstopo = xr.Dataset(
        {'DHM25' : (['chy', 'chx'], destination.astype(np.int16))},
        coords={
            'chx' : ('chx', dst_grid[0][0, :].astype(np.float32)),
            'chy' : ('chy', dst_grid[1][:, 0].astype(np.float32)),
            'lon' : (['chy', 'chx'], ref_coords[0].astype(np.float32)),
            'lat' : (['chy', 'chx'], ref_coords[1].astype(np.float32))
            }
        # TODO: add attributes
    )
print(swisstopo)

# get a view on SRTM domain to merge it with the reprojected Swisstopo data
srtm_swiss = srtm_ccs4.sel(
    chx=slice(swisstopo.chx.min(), swisstopo.chx.max()),
    chy=slice(swisstopo.chy.min(), swisstopo.chy.max())
)
dems_swiss = srtm_swiss.merge(swisstopo)
print(dems_swiss)

# set NaNs and compare the two DEMs
dems_swiss = dems_swiss.where(dems_swiss.DHM25 > 0).where(dems_swiss.SRTM90 > 0)
squared_errors = (dems_swiss.SRTM90 - dems_swiss.DHM25) ** 2
rmse = np.sqrt(squared_errors.mean().values)
bias = (dems_swiss.SRTM90 - dems_swiss.DHM25).mean().values
plt.close()
xr.plot.scatter(dems_swiss, 'DHM25', 'SRTM90', s=1, alpha=0.3)
plt.plot([0, 4809], [0, 4809], '--k')
plt.xlim([0, 4809])
plt.ylim([0, 4809])
plt.grid(ls=':')
plt.title(f'RMSE = {rmse:.1f}m, BIAS = {bias:.1f}m')
plt.tight_layout()
plt.savefig('scatter_DHM25_SRTM90.png', dpi=300)

# use DHM25 where available to compile a new product
dem_swiss = srtm_swiss.SRTM90.where(swisstopo.DHM25 < 0, swisstopo.DHM25)
print(dem_swiss)

# and finally merge the combined SRTM/Swisstopo from the Swiss domain into
# the larger CCS4 domain
idx_x = np.logical_and(
            srtm_ccs4.chx >= dem_swiss.chx.min(),
            srtm_ccs4.chx <= dem_swiss.chx.max())
idx_y = np.logical_and(
            srtm_ccs4.chy >= dem_swiss.chy.min(),
            srtm_ccs4.chy <= dem_swiss.chy.max())
dem_ccs4 = srtm_ccs4.SRTM90
dem_ccs4[idx_y, idx_x] = dem_swiss

# convert to dataset
dem_ccs4 = dem_ccs4.to_dataset(name='DEM')
print(dem_ccs4)

# plot and save
plt.close()
dem_ccs4.DEM.where(dem_ccs4.DEM > 0, 0).plot.imshow(cmap='terrain', vmin=0, vmax=4800)
plt.tight_layout()
plt.savefig(f'DEM_CCS4_50m_buffer', dpi=300)
print('saving to netcdf...', end=' ', flush=True)
dem_ccs4.to_netcdf('DEM_CCS4_50m_buffer.nc')
print('done!')
