'''Script to plot topographic descriptors maps'''

from env_setting import DATAPATH, RESOURCES, FIGPATH

import copy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import Topo_descriptors.topo_helpers as hlp
import xarray as xr

# optional dependency, plot national borders
try:
    import noborders as nob
    NOBORDERS_IMPORTED = True
except ModuleNotFoundError:
    NOBORDERS_IMPORTED = False

inca_domain = {'chx':slice(255000,965000), 'chy':slice(-160000,480000)}
zoom_domain = {'chx':slice(654000,704000), 'chy':slice(115500,165500)}
zoom = patches.Rectangle(
    (zoom_domain['chx'].start, zoom_domain['chy'].start),
    zoom_domain['chx'].stop - zoom_domain['chx'].start,
    zoom_domain['chy'].stop - zoom_domain['chy'].start,
    linewidth=1.5,
    edgecolor='k',
    facecolor='none'
    )

# ds = xr.open_dataset(str(DATAPATH / 'topo_TPI_60000M.nc'))
ds = xr.open_mfdataset(str(DATAPATH / 'topo_*.nc'), combine='nested', concat_dim=None)

for var in ds.data_vars:

    fig = plt.figure(figsize=(12, 4.5))

    if 'ASPECT_' in var or 'DIR_' in var:
        robust = False
    else:
        robust = True

    ax1 = plt.subplot(1, 2, 1)
    ds[var].plot.imshow(robust=robust, ax=ax1)
    if NOBORDERS_IMPORTED: nob.addBorders(color="black", linewidth=1)
    ax1.add_patch(copy.copy(zoom))
    ax1.set_aspect(1)

    ax2 = plt.subplot(1, 2, 2)
    ds[var].sel(zoom_domain).plot.imshow(robust=robust, ax=ax2)
    if NOBORDERS_IMPORTED: nob.addBorders(color="black")
    ax2.set_xlim([zoom_domain['chx'].start, zoom_domain['chx'].stop])
    ax2.set_ylim([zoom_domain['chy'].start, zoom_domain['chy'].stop])
    ax2.set_aspect(1)

    fig.tight_layout()
    fn = str(FIGPATH / f'{var}.png')
    fig.savefig(fn)
    print(f'saved: {fn}')
    plt.close()
