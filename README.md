# topo-descriptors

A python library to compute DEM-based topographical descriptors.

## Usage

Let's install `topo-descriptors` with few additional packages that will help us
to run a simple example (remember to use a virtual environment):


```python
%pip install topo-descriptors elevation rioxarray matplotlib --quiet
```

    Note: you may need to restart the kernel to use updated packages.


The [elevation](https://github.com/bopen/elevation) package is an python library that
provides an easy access to global elevation data. Here we are going to clip the SRTM 30m
DEM around the Basodino region in southern Switzerland, around 46.4N 8.5E:


```python
!eio clip -o Basodino-30m-DEM.tif --bounds 8.2 46.30 8.6 46.55
```

    make: Entering directory '/prod/gve/home/mts/.cache/elevation/SRTM1'
    make: Nothing to be done for 'download'.
    make: Leaving directory '/prod/gve/home/mts/.cache/elevation/SRTM1'
    make: Entering directory '/prod/gve/home/mts/.cache/elevation/SRTM1'
    make: Nothing to be done for 'all'.
    make: Leaving directory '/prod/gve/home/mts/.cache/elevation/SRTM1'
    make: Entering directory '/prod/gve/home/mts/.cache/elevation/SRTM1'
    cp SRTM1.vrt SRTM1.d73260d0233a450ab8ca8ce05b9b46c6.vrt
    make: Leaving directory '/prod/gve/home/mts/.cache/elevation/SRTM1'
    make: Entering directory '/prod/gve/home/mts/.cache/elevation/SRTM1'
    gdal_translate -q -co TILED=YES -co COMPRESS=DEFLATE -co ZLEVEL=9 -co PREDICTOR=2 -projwin 8.2 46.55 8.6 46.3 SRTM1.d73260d0233a450ab8ca8ce05b9b46c6.vrt /prod/gve/home/mts/git/topo-descriptors/Basodino-30m-DEM.tif
    rm -f SRTM1.d73260d0233a450ab8ca8ce05b9b46c6.vrt
    make: Leaving directory '/prod/gve/home/mts/.cache/elevation/SRTM1'



```python
import logging

logger = logging.getLogger()
handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s %(name)-12s %(levelname)-8s %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.INFO)
```

Now in python we can easily import the
`Basodino-30m-DEM.tif` file generated above:


```python
from topo_descriptors.helpers import get_dem_netcdf, scale_to_pixel

dem_ds = get_dem_netcdf("Basodino-30m-DEM.tif")
varname = list(dem_ds)[0]
dem_ds.attrs.update(crs="epsg:4326")
dem_ds[varname].plot(robust=True)
```

    2023-12-20 17:53:55,120 yaconfigobject INFO     Loading /prod/gve/home/mts/git/topo-descriptors/topo_descriptors/config/topo_descriptors.conf.
    2023-12-20 17:53:55,121 yaconfigobject INFO     Loading configuration file: /prod/gve/home/mts/git/topo-descriptors/topo_descriptors/config/topo_descriptors.conf





    <matplotlib.collections.QuadMesh at 0x7f060ab21fd0>




    
![png](README_files/README_6_2.png)
    



```python
from topo_descriptors import topo

scale_meters = 500
scale_pixel, __ = scale_to_pixel(scale_meters, dem_ds)
topo.tpi(dem_ds[varname], scale_pixel).plot(vmin=-100, vmax=100, cmap="bwr")
```

    2023-12-20 17:53:58,736 topo_descriptors.helpers INFO     Computed in 0:00:00 (HH:mm:ss)





    <matplotlib.collections.QuadMesh at 0x7f0601529750>




    
![png](README_files/README_7_2.png)
    


The Sx is used to describe the horizon in a given direction and spatial scale.
In the example below we compute the Sx for a 0° azimuth (i.e., looking North)
and a radius of 500 meters.


```python
import xarray as xr

sx_500m = topo.sx(dem_ds, azimuth=0, radius=500)
xr.DataArray(sx_500m, coords=dem_ds.coords).plot.imshow()
```

    2023-12-20 17:54:06,246 topo_descriptors.helpers INFO     Computed in 0:00:06 (HH:mm:ss)





    <matplotlib.image.AxesImage at 0x7f05f15b8610>




    
![png](README_files/README_9_2.png)
    


Other topographical descriptors are available, such as slope, aspect, derivatives,
and more. As an example, below we show how to compute a range of descriptors for two
distinc spatial scales (200 and 2000 meters).


```python
from pathlib import Path

output_dir = Path("out/")
output_dir.mkdir(exist_ok=True)

scales_meters = [200, 2000]
domain = {"x": slice(8.25, 8.55), "y": slice(46.50, 46.35)}

topo.compute_gradient(dem_ds, scales_meters, sig_ratios=1, crop=domain, outdir=output_dir)
topo.compute_std(dem_ds, scales_meters, crop=domain, outdir=output_dir)
topo.compute_tpi(dem_ds, scales_meters, crop=domain, outdir=output_dir)
topo.compute_sx(dem_ds, azimuth=0, radius=scales_meters[0], crop=domain, outdir=output_dir)
topo.compute_sx(dem_ds, azimuth=0, radius=scales_meters[1], crop=domain, outdir=output_dir)
```

    2023-12-20 17:54:06,586 topo_descriptors.topo INFO     ***Starting gradients computation for scales [200, 2000] meters***
    2023-12-20 17:54:06,870 topo_descriptors.topo INFO     Computing scale 200 meters with sigma ratio 1 ...
    2023-12-20 17:54:06,920 topo_descriptors.helpers INFO     Computed in 0:00:00 (HH:mm:ss)
    2023-12-20 17:54:06,945 topo_descriptors.helpers INFO     saved: out/topo_WE_DERIVATIVE_200M_SIGRATIO1.nc
    2023-12-20 17:54:06,958 topo_descriptors.helpers INFO     saved: out/topo_SN_DERIVATIVE_200M_SIGRATIO1.nc
    2023-12-20 17:54:06,970 topo_descriptors.helpers INFO     saved: out/topo_SLOPE_200M_SIGRATIO1.nc
    2023-12-20 17:54:06,982 topo_descriptors.helpers INFO     saved: out/topo_ASPECT_200M_SIGRATIO1.nc
    2023-12-20 17:54:06,983 topo_descriptors.topo INFO     Computing scale 2000 meters with sigma ratio 1 ...
    2023-12-20 17:54:07,229 topo_descriptors.helpers INFO     Computed in 0:00:00 (HH:mm:ss)
    2023-12-20 17:54:07,242 topo_descriptors.helpers INFO     saved: out/topo_WE_DERIVATIVE_2000M_SIGRATIO1.nc
    2023-12-20 17:54:07,260 topo_descriptors.helpers INFO     saved: out/topo_SN_DERIVATIVE_2000M_SIGRATIO1.nc
    2023-12-20 17:54:07,272 topo_descriptors.helpers INFO     saved: out/topo_SLOPE_2000M_SIGRATIO1.nc
    2023-12-20 17:54:07,284 topo_descriptors.helpers INFO     saved: out/topo_ASPECT_2000M_SIGRATIO1.nc
    2023-12-20 17:54:07,286 topo_descriptors.topo INFO     ***Starting STD computation for scales [200, 2000] meters***
    2023-12-20 17:54:07,569 topo_descriptors.topo INFO     Computing scale 200 meters with smoothing factor None ...
    2023-12-20 17:54:07,696 topo_descriptors.helpers INFO     Computed in 0:00:00 (HH:mm:ss)
    2023-12-20 17:54:07,717 topo_descriptors.helpers INFO     saved: out/topo_STD_200M.nc
    2023-12-20 17:54:07,719 topo_descriptors.topo INFO     Computing scale 2000 meters with smoothing factor None ...
    2023-12-20 17:54:07,853 topo_descriptors.helpers INFO     Computed in 0:00:00 (HH:mm:ss)
    2023-12-20 17:54:07,871 topo_descriptors.helpers INFO     saved: out/topo_STD_2000M.nc
    2023-12-20 17:54:07,873 topo_descriptors.topo INFO     ***Starting TPI computation for scales [200, 2000] meters***
    2023-12-20 17:54:08,141 topo_descriptors.topo INFO     Computing scale 200 meters with smoothing factor None ...
    2023-12-20 17:54:08,182 topo_descriptors.helpers INFO     Computed in 0:00:00 (HH:mm:ss)
    2023-12-20 17:54:08,196 topo_descriptors.helpers INFO     saved: out/topo_TPI_200M.nc
    2023-12-20 17:54:08,197 topo_descriptors.topo INFO     Computing scale 2000 meters with smoothing factor None ...
    2023-12-20 17:54:08,235 topo_descriptors.helpers INFO     Computed in 0:00:00 (HH:mm:ss)
    2023-12-20 17:54:08,255 topo_descriptors.helpers INFO     saved: out/topo_TPI_2000M.nc
    2023-12-20 17:54:08,257 topo_descriptors.topo INFO     ***Starting Sx computation for azimuth 0 meters and radius 200***
    2023-12-20 17:54:09,313 topo_descriptors.helpers INFO     Computed in 0:00:01 (HH:mm:ss)
    2023-12-20 17:54:09,326 topo_descriptors.helpers INFO     saved: out/topo_SX_RADIUS200_AZIMUTH0.nc
    2023-12-20 17:54:09,328 topo_descriptors.topo INFO     ***Starting Sx computation for azimuth 0 meters and radius 2000***
    2023-12-20 17:54:15,750 topo_descriptors.helpers INFO     Computed in 0:00:06 (HH:mm:ss)
    2023-12-20 17:54:15,764 topo_descriptors.helpers INFO     saved: out/topo_SX_RADIUS2000_AZIMUTH0.nc


Above, the output was written directly to disk, while in the cell below we show how 
to easly import the results and visualize them using xarray.


```python
ds = xr.open_mfdataset(str(output_dir / "topo_*.nc"))
min_max = ds.quantile(q=[0.05, 0.95])
ds = (ds - min_max.isel(quantile=0)) / (
    min_max.isel(quantile=1) - min_max.isel(quantile=0)
)
ds.to_array().plot.imshow(
    col="variable",
    col_wrap=len(scales_meters),
    robust=True,
    add_colorbar=False,
    vmin=0,
    vmax=1,
)
ds.close()
```


    
![png](README_files/README_13_0.png)
    


## Build the README

To use this Jupyter Notebook to compile the markdown's version for GitHub, first install
the conda environment using the `environment.yaml` file:

```shell
conda env create -f environment.yaml
conda activate topo
```

Then generate the `README.md` by running:

```shell
jupyter nbconvert --execute --to markdown README.ipynb
```

The associated figures are saved in the `README_files/` folder.
