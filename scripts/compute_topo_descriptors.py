"""
Example script on how to compute spatial descriptors from a DEM file.
"""

import logging

import topo_descriptors.topo as tp
import topo_descriptors.helpers as hlp
from topo_descriptors import CFG

logger = logging.getLogger(__name__)

if __name__ == "__main__":

    logging.basicConfig(level=logging.INFO)
    logging.captureWarnings(True)

    # get the DEM
    path_dem = "DEM.nc"
    dem_da = hlp.get_dem_netcdf(path_dem)
    ind_nans, dem_da = hlp.fill_na(dem_da)

    # define the target domain
    domain = {"x": slice(255000, 965000), "y": slice(-160000, 480000)}

    # define the convolution scales in meters
    scales_meters = [
        100,
        300,
        500,
        1000,
        2000,
        4000,
        6000,
        10000,
        20000,
        30000,
        60000,
        100000,
    ]

    # Launch computations and save output

    # raw TPI
    tp.compute_tpi(
        dem_da, scales_meters, smth_factors=None, ind_nans=ind_nans, crop=domain
    )

    # TPI with prior smoothing
    tp.compute_tpi(
        dem_da, scales_meters, smth_factors=1, ind_nans=ind_nans, crop=domain
    )

    # Gradients with symmetric kernels
    tp.compute_gradient(
        dem_da, scales_meters, sig_ratios=1, ind_nans=ind_nans, crop=domain
    )

    # Valley Index with prior smoothing
    tp.compute_valley_ridge(
        dem_da,
        scales_meters[3:],
        mode="valley",
        flat_list=[0, 0.2, 0.4],
        smth_factors=0.5,
        ind_nans=ind_nans,
        crop=domain,
    )

    # Ridge Index with prior smoothing
    tp.compute_valley_ridge(
        dem_da,
        scales_meters[3:],
        mode="ridge",
        flat_list=[0, 0.15, 0.3],
        smth_factors=0.5,
        ind_nans=ind_nans,
        crop=domain,
    )
