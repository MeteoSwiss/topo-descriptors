import numpy as np
import numpy.ma as ma
import xarray as xr
from scipy import ndimage, signal
from skimage.transform import rescale, resize, downscale_local_mean
import Topo_descriptors.topo_helpers as hlp


def compute_tpi(dem_da, dist_list):
    
    if not hasattr(dist_list, '__iter__'): dist_list = [dist_list]
    
    dist_pxl = hlp.dist_to_pixel(dist_list, dem_da)
    
    for idx,dist in enumerate(dist_pxl):
        array = tpi(dem_da.values, dist)
        name = 'TPI_' + str(dist_list[idx]) + 'KM'
        hlp.save_topo_netcdf(hlp._to_da(array, dem_da, name))
        print('TPI distance ' + str(dist_list[idx]) + ' finished')
        del array
        
def tpi(array, size):
        
    xx, yy = np.mgrid[:size, :size]
    middle = np.floor(size/2)
    circle = (xx - middle) ** 2 + (yy - middle) ** 2
    kernel = np.asarray(circle <= (middle**2), dtype= np.float32)
    kernel = kernel/np.sum(kernel)
    conv = ndimage.convolve(array, kernel, mode='reflect')
    return (array - conv)
        

def valley_index(array, size, flat_list):
         
    n_y, n_x = array.shape
    n_kernels = len(flat_list) + 1
    array = np.broadcast_to(array,(n_kernels,n_y,n_x))
    angle_list = np.arange(0,180,dtype=np.float32)
    valley_index = np.zeros((n_y, n_x),dtype=np.float32) - np.inf
    valley_angle = np.empty((n_y, n_x),dtype=np.float32)
    
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
        kernel_rotated = ma.masked_array(kernel_rotated, mask=kernel_rotated == -999)
        kernel_rotated = (kernel_rotated - np.mean(kernel_rotated,axis=(1,2),keepdims=True))/np.std(kernel_rotated,axis=(1,2),keepdims=True)
        kernel_rotated = ma.MaskedArray.filled(kernel_rotated,0)
        dem_convolved = signal.convolve(array,kernel_rotated, mode = 'same', method='auto')
        dem_convolved = np.max(dem_convolved, axis=0)
        bool_greater = dem_convolved > valley_index
        valley_index[bool_greater] = dem_convolved[bool_greater]
        valley_angle[bool_greater] = angle
        del bool_greater 

    valley_index = np.ndarray.clip(valley_index, min=0)
    valley_index = valley_index.where(~np.isnan(array))
    valley_index = (valley_index - valley_index.mean())/valley_index.std()
    valley_angle = valley_angle * np.pi/180 
    valley_cos = valley_index * np.abs(np.cos(valley_angle))
    valley_sin = valley_index * np.sin(valley_angle)
    
    return valley_cos, valley_sin
    
def compute_valley_index(dem_da, dist_list, flat_list=[0.15,0.3]):
    
    if not hasattr(dist_list, '__iter__'): dist_list = [dist_list]
    dist_pxl = hlp.dist_to_pixel(dist_list, dem_da)
    
    for idx,dist in enumerate(dist_pxl):
        array = valley_index(dem_da.values, dist)
        name = 'VALLEY_INDEX_' + str(dist_list[idx]) + 'KM'
        hlp.save_topo_netcdf(hlp._to_da(array, dem_da, name))
        print('Valley index distance ' + str(dist_list[idx]) + ' finished')
        del array


# TODO: change valley function to have option 'valley/ridge'
        # TODO: wrapper for upsampling and down sampling


#def valley_netcdf(path_dem, dist_list, flat_list):
#     
#    if not hasattr(dist_list, '__iter__'): dist_list = [dist_list]
#    if not hasattr(flat_list, '__iter__'): flat_list = [flat_list]
#    dem_da, dist_pxl, res = preprocess_dem(path_dem, dist_list)
#    dem_in = (dem_da - dem_da.mean())/ dem_da.std()
#    dem_in = dem_in.interpolate_na(dim='chx', method='nearest', fill_value='extrapolate').values
#
#    n_y, n_x = dem_in.shape
#    n_dist = len(dist_list)
#    valley_index = np.empty((n_dist,n_y, n_x), dtype= np.float32)
#    valley_angle = np.empty((n_dist,n_y, n_x), dtype= np.float32)
#    n_kernels = len(flat_list) + 1
#    dem_in = np.broadcast_to(dem_in,(n_kernels,n_y,n_x))
#    angle_list = np.arange(0,180,dtype=np.float32)
#    
#    for idx,size in enumerate(dist_pxl):
#        
#        current_max = np.zeros((n_y, n_x),dtype=np.float32) - np.inf
#        current_anglemax = np.empty((n_y, n_x),dtype=np.float32)
#        
#        middle = int(np.floor(size/2))
#        kernel_tmp = np.broadcast_to(np.arange(0,middle+1),(size,middle+1)).T
#        kernel_tmp = np.concatenate((np.flip(kernel_tmp[1:,:],axis=0),kernel_tmp),axis=0)
#        kernel_tmp = np.asarray(kernel_tmp,dtype=np.float32)
#        kernel_all = np.broadcast_to(kernel_tmp, (n_kernels,size,size)).copy()
#        
#        for id2,flat in enumerate(flat_list,1):
#             
#            halfwidth = int(np.floor(np.floor(size*flat/2)+0.5))
#            kernel_all[id2,middle-halfwidth:middle+halfwidth+1,:] = kernel_all[id2,middle-halfwidth,0]
#            kernel_all = (kernel_all - np.mean(kernel_all,axis=(1,2),keepdims=True))/np.std(kernel_all,axis=(1,2),keepdims=True)
#            
#        for angle in angle_list:
#            
#            kernel_rotated = ndimage.rotate(kernel_all, angle , axes=(1,2), reshape=True, order=2, mode='constant', cval=-999)
#            dem_convolved = np.empty((n_kernels,n_y, n_x), dtype=np.float32)
#            kernel_rotated = ma.masked_array(kernel_rotated, mask=kernel_rotated == -999)
#            kernel_rotated = (kernel_rotated - np.mean(kernel_rotated,axis=(1,2),keepdims=True))/np.std(kernel_rotated,axis=(1,2),keepdims=True)
#            kernel_rotated = ma.MaskedArray.filled(kernel_rotated,0)
#            dem_convolved = signal.convolve(dem_in,kernel_rotated, mode = 'same', method='auto')
#            dem_convolved = np.max(dem_convolved, axis=0)
#            bool_greater = dem_convolved > current_max
#            current_max[bool_greater] = dem_convolved[bool_greater]
#            current_anglemax[bool_greater] = angle
#            del bool_greater 
#        
#        valley_index[idx,:,:] = current_max  
#        valley_angle[idx,:,:] = current_anglemax
#        print('Size ' + str(dist_list[idx]) + ' finished')
#            
#    valley_index = np.ndarray.clip(valley_index, min = 0)
#    valley_index = xr.DataArray(valley_index, dims=('scale','chy','chx'))
#    valley_index = valley_index.where(~np.isnan(dem_da))
#    valley_index = (valley_index - valley_index.mean(dim=('chy','chx')))/valley_index.std(dim=('chy','chx'))
#    valley_angle = valley_angle * np.pi/180 
#    valley_cos = valley_index * np.abs(np.cos(valley_angle))
#    valley_sin = valley_index * np.sin(valley_angle)
#    
#    array_to_netcdf(valley_cos, dem_da, dist_list, f'valley_cos{str(res)}.nc')
#    array_to_netcdf(valley_sin, dem_da, dist_list, f'valley_sin{str(res)}.nc')


