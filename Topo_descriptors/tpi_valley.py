import numpy as np
from scipy import ndimage, signal
import xarray as xr
import Topo_descriptors.topo_helpers as hlp
from env_setting import DATAPATH


def tpi_netcdf(path_dem, dist_list):
    
    if not hasattr(dist_list, '__iter__'): dist_list = [dist_list]
    
    dem_da, dist_pxl, res = preprocess_dem(path_dem, dist_list)
    tpi = np.empty(dist_pxl.shape +  dem_da.shape, dtype= np.float32)
    
    for idx,dist in enumerate(dist_pxl):
        xx, yy = np.mgrid[:dist, :dist]
        middle = np.floor(dist/2)
        circle = (xx - middle) ** 2 + (yy - middle) ** 2
        kernel = np.asarray(circle <= (middle**2), dtype= np.float32)
        kernel = kernel/np.sum(kernel)
        conv = ndimage.convolve(dem_da, kernel, mode='reflect')
        conv = (dem_da - conv)
        tpi[idx,:,:] = conv
        print('Distance ' + str(dist_list[idx]) + ' finished')
    
    array_to_netcdf(tpi, dem_da, dist_list, f'tpi{str(int(res))}.nc')
    
    
def valley_netcdf(path_dem, dist_list, flat_list):
     
    if not hasattr(dist_list, '__iter__'): dist_list = [dist_list]
    if not hasattr(flat_list, '__iter__'): flat_list = [flat_list]
    dem_da, dist_pxl, res = preprocess_dem(path_dem, dist_list)
    dem_in = (dem_da - dem_da.mean())/ dem_da.std()
    dem_in = dem_in.interpolate_na(dim='chx', method='nearest', fill_value='extrapolate').values

    n_y, n_x = dem_in.shape
    n_dist = len(dist_list)
    valley_index = np.empty((n_dist,n_y, n_x), dtype= np.float32)
    valley_angle = np.empty((n_dist,n_y, n_x), dtype= np.float32)
    n_kernels = len(flat_list) + 1
    angle_list = np.arange(0,180,dtype=np.float32)
    
    for idx,size in enumerate(dist_pxl):
        
        current_max = np.zeros((n_y, n_x),dtype=np.float32) - np.inf
        current_anglemax = np.empty((n_y, n_x),dtype=np.float32)
        
        middle = int(np.floor(size/2))
        kernel_tmp = np.broadcast_to(np.arange(0,middle+1),(size,middle+1)).T
        kernel_tmp = np.concatenate((np.flip(kernel_tmp[1:,:],axis=0),kernel_tmp),axis=0)
        kernel_tmp = np.asarray(kernel_tmp,dtype=np.float32)
        kernel_all = np.broadcast_to(kernel_tmp, (n_kernels,size,size)).copy()
        
        for id2,flat in enumerate(flat_list,1):
             
            halfwidth = int(np.floor(np.floor(size*flat/2)+0.5))
            kernel_all[id2,middle-halfwidth:middle+halfwidth+1,:] = kernel_all[id2,middle-halfwidth,0]
            kernel_all = (kernel_all - np.mean(kernel_all,axis=(1,2),keepdims=True))/np.std(kernel_all,axis=(1,2),keepdims=True)
            
        for angle in angle_list:
            
            kernel_rotated = ndimage.rotate(kernel_all, angle , axes=(1,2), reshape=True, order=2, mode='constant', cval=-999)
            dem_convolved = np.empty((n_kernels,n_y, n_x), dtype=np.float32)
            
            for id3,filters in enumerate(kernel_rotated):
                
                filter_mask = filters == -999
                filter_good = np.logical_not(filter_mask)
                filters = (filters - np.mean(filters[filter_good]))/np.std(filters[filter_good])
                filters[filter_mask] = 0
                
                dem_convolved[id3,:,:] = signal.convolve(dem_in,filters, mode = 'same')
             
            dem_convolved = np.max(dem_convolved, axis=0)
            bool_greater = dem_convolved > current_max
            current_max[bool_greater] = dem_convolved[bool_greater]
            current_anglemax[bool_greater] = angle
            del bool_greater 
        
        valley_index[idx,:,:] = current_max  
        valley_angle[idx,:,:] = current_anglemax
        print('Size ' + str(dist_list[idx]) + ' finished')
            
    valley_index = np.ndarray.clip(valley_index, min = 0)
    valley_index = xr.DataArray(valley_index, dims=('scale','chy','chx'))
    valley_index = valley_index.where(~np.isnan(dem_da))
    valley_index = (valley_index - valley_index.mean(dim=('chy','chx')))/valley_index.std(dim=('chy','chx'))
    valley_angle = valley_angle * np.pi/180 
    valley_cos = valley_index * np.abs(np.cos(valley_angle))
    valley_sin = valley_index * np.sin(valley_angle)
    
    array_to_netcdf(valley_cos, dem_da, dist_list, f'valley_cos{str(res)}.nc')
    array_to_netcdf(valley_sin, dem_da, dist_list, f'valley_sin{str(res)}.nc')


def array_to_netcdf(array,dem_da,dist_list, name):
    
    array = xr.DataArray(array,
                coords=[('scale', dist_list),
                ('chy', dem_da['chy']),
                ('chx', dem_da['chx'])]).to_dataset(name=str.upper(name))
    
    array.to_netcdf(DATAPATH / str.lower(name))


def preprocess_dem(path_dem, dist_list):
    
    dem = hlp.get_dem_netcdf(path_dem)
    dem = dem.where(dem != -9999.)
    res = dem['chx'].diff('chx').mean().values
    dist_pxl = hlp.round_up_to_odd(dist_list / res)
    
    return dem, dist_pxl, res