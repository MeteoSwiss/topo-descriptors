# -*- coding: utf-8 -*-
import numpy as np
from pyproj import Proj, transform
from multiprocessing import Pool
import Topo_descriptors.topo_helpers as hlp


def get_sxrelief_geotiff(theta, alpha, d_min, d_max, height, path_dem, path_out):
    """compute Sx descriptor using fancy indexing and save as Geotiff"""
    
    (dem, lat, lon, geotransform, proj) = hlp.get_dem_all(path_dem)
    
    sizes_dem = dem.shape
   
    lat = np.reshape(lat,(-1))
    lon = np.reshape(lon,(-1))
    dem = np.reshape(dem,(-1))
    
    n_flat_dem = len(dem)
    ind_array = np.asarray(np.reshape(np.arange(n_flat_dem),sizes_dem), dtype=np.int32)
    LV03 = Proj(init='epsg:21781')  # LV03 CRS for distances in meters
    WGS84 = Proj(init='epsg:4326')  # for lat, lon [deg]
    x, y = transform(WGS84, LV03,  lon, lat)
    
    coord = np.asarray([x,y], dtype=np.float32)
    dx1 = np.abs(coord[0,1]-coord[0,0])
    dx2 = np.abs(coord[0,-1]-coord[0,-2])
    dy1 = np.abs(coord[1,ind_array[1,0]]-coord[1,ind_array[0,0]])
    dy2 = np.abs(coord[1,ind_array[-1,0]]-coord[1,ind_array[-2,0]])
    size_pixel = min(dx1,dx2,dy1,dy2) - 1
    max_dmax = np.max(d_max)
    d_max_pixel = int(np.ceil(max_dmax/size_pixel))
    d_max_squared = max_dmax**2
    del x, y
    
    angle_list = np.arange(-90,270,theta)
    angle_list = angle_list * np.pi/180
    alpha_rad = alpha * np.pi/180
    angles_names_sx = [path_out / ('Sx_' + str(ind) + '.tif') for ind in range(0,360,theta)]
    angles_names_relief = [path_out / ('Relief_' + str(ind) + '.tif') for ind in range(0,360,theta)]
    
    pool = Pool(18)
    for id_angle, angle_rad in enumerate(angle_list):
        print('Angle: '+ str(id_angle))
        pool.apply_async(sx_one_rotation, args = (dem, angle_rad, coord, n_flat_dem,
        ind_array, d_max_pixel, sizes_dem, alpha_rad, d_min, d_max, d_max_squared, 
        height, angles_names_sx[id_angle], angles_names_relief[id_angle], geotransform, proj))

    pool.close()
    pool.join()
    
    return 
    

def sx_one_rotation(dem, angle_rad, coord, n_flat_dem, ind_array, d_max_pixel, 
    sizes_dem, alpha_rad, d_min, d_max, d_max_squared, height, angle_name_sx, angle_name_relief, geotransform, proj):
    print('starting process')
    
    n_dist = len(d_min)
    
    # rotation matrice
    cos_angle = np.cos(angle_rad)
    sin_angle = np.sin(angle_rad)
    R = np.array([[cos_angle,-sin_angle],[sin_angle,cos_angle]], dtype=np.float32)
    coord_rotate = R@coord
    del coord
    sx_tmp = np.empty((n_dist, n_flat_dem), dtype=np.float32)
    relief_tmp = np.empty((n_dist, n_flat_dem), dtype=np.float32)
    
    for idx,(x0,y0) in enumerate(coord_rotate.T):
        
        row_idx, col_idx = np.divmod(idx,sizes_dem[1])
        start_rows = max(0,row_idx-d_max_pixel)
        stop_rows = min(sizes_dem[0],row_idx+d_max_pixel)
        start_col = max(0,col_idx-d_max_pixel)
        stop_col = min(sizes_dem[1],col_idx+d_max_pixel)
        ind_dem = ind_array[start_rows:stop_rows,start_col:stop_col].flatten()
        
        diff_x = coord_rotate[0,ind_dem] - x0
        bool_diff = diff_x > 0
        diff_x = diff_x[bool_diff]
        ind_dem = ind_dem[bool_diff]
        diff_y = coord_rotate[1,ind_dem] - y0
        
        tan_abs = np.abs((diff_y)/(diff_x))
        
        bool_tan = tan_abs < np.tan(alpha_rad)
        ind_dem = ind_dem[bool_tan]
        diff_x = diff_x[bool_tan]
        diff_y = diff_y[bool_tan]
        
        dist = diff_x**2 + diff_y**2
        bool_dist = dist < d_max_squared
        ind_dem = ind_dem[bool_dist] 
        dist = dist[bool_dist]
        dist = np.sqrt(dist)
        sx_tmp[:,idx], relief_tmp[:,idx] = compute_sxrelief(dist, dem[ind_dem], dem[idx], height, d_min, d_max)
    
    del dem, ind_array, coord_rotate
    sx_tmp =  np.reshape(sx_tmp, (n_dist,) + sizes_dem)
    hlp.array_to_tif(sx_tmp, angle_name_sx, sizes_dem, geotransform, proj, n_dist)
    del sx_tmp
    relief_tmp =  np.reshape(relief_tmp, (n_dist,) + sizes_dem)
    hlp.array_to_tif(relief_tmp, angle_name_relief, sizes_dem, geotransform, proj, n_dist)

    return


def compute_sxrelief(dist_reg, dem_reg, dem_pixel, height, d_min, d_max):
    
    delta = dem_reg - (dem_pixel + height)
    phi = np.arctan(delta / dist_reg)
    sx_alldist = np.empty(d_min.shape, dtype = np.float32)
    relief_alldist = sx_alldist.copy()
    dist_reg = dist_reg[:,np.newaxis]
    bool_alldist = np.logical_and(dist_reg > d_min, dist_reg < d_max)
    for ind,col in enumerate(bool_alldist.T):
        try:
            sx = np.max(phi[col])
            relief = - np.mean(delta[col])
        except:
            sx = 0
            relief = 0
            
        sx_alldist[ind] = sx
        relief_alldist[ind] = relief 
    
    return (sx_alldist, relief_alldist)
