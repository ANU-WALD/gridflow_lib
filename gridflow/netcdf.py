from osgeo import gdal
import numpy as np
import netCDF4
import json
from datetime import datetime

def pack_fmc(hdf_file, date, mean_arr, std_arr, q_mask, dest):
    
    with netCDF4.Dataset(dest, 'w', format='NETCDF4_CLASSIC') as ds:
        with open('nc_metadata.json') as data_file:
            attrs = json.load(data_file)
            for key in attrs:
                setattr(ds, key, attrs[key])
        setattr(ds, "date_created", datetime.now().strftime("%Y%m%dT%H%M%S"))
        
        rast = gdal.Open('HDF4_EOS:EOS_GRID:"{}":MOD_Grid_BRDF:Nadir_Reflectance_Band1'.format(hdf_file))
        proj_wkt = rast.GetProjection()
        geot = rast.GetGeoTransform()
        
        t_dim = ds.createDimension("time", 1)
        x_dim = ds.createDimension("x", rast.RasterXSize)
        y_dim = ds.createDimension("y", rast.RasterYSize)

        var = ds.createVariable("time", "f8", ("time",))
        var.units = "seconds since 1970-01-01 00:00:00.0"
        var.calendar = "standard"
        var.long_name = "Time, unix time-stamp"
        var.standard_name = "time"
        var[:] = netCDF4.date2num([date], units="seconds since 1970-01-01 00:00:00.0", calendar="standard")

        var = ds.createVariable("x", "f8", ("x",))
        var.units = "m"
        var.long_name = "x coordinate of projection"
        var.standard_name = "projection_x_coordinate"
        var[:] = np.linspace(geot[0], geot[0]+(geot[1]*rast.RasterXSize), rast.RasterXSize)
        
        var = ds.createVariable("y", "f8", ("y",))
        var.units = "m"
        var.long_name = "y coordinate of projection"
        var.standard_name = "projection_y_coordinate"
        var[:] = np.linspace(geot[3], geot[3]+(geot[5]*rast.RasterYSize), rast.RasterYSize)
        
        var = ds.createVariable("lfmc_mean", 'f4', ("time", "y", "x"), fill_value=-9999.9)
        var.long_name = "LFMC Arithmetic Mean"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = mean_arr[None,...]

        var = ds.createVariable("lfmc_stdv", 'f4', ("time", "y", "x"), fill_value=-9999.9)
        var.long_name = "LFMC Standard Deviation"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = std_arr[None,...]
        
        var = ds.createVariable("quality_mask", 'i1', ("time", "y", "x"), fill_value=0)
        var.long_name = "Combined Bands Quality Mask"
        var.units = 'Cat'
        var.grid_mapping = "sinusoidal"
        var[:] = q_mask.astype(np.int8)[None,...]

        var = ds.createVariable("sinusoidal", 'S1', ())
        var.grid_mapping_name = "sinusoidal"
        var.false_easting = 0.0
        var.false_northing = 0.0
        var.longitude_of_central_meridian = 0.0
        var.longitude_of_prime_meridian = 0.0
        var.semi_major_axis = 6371007.181
        var.inverse_flattening = 0.0
        var.spatial_ref = proj_wkt
        var.GeoTransform = "{} {} {} {} {} {} ".format(*[geot[i] for i in range(6)])


def pack_flammability(fmc_file, date, flam, anom, q_mask, dest):
    
    with netCDF4.Dataset(dest, 'w', format='NETCDF4_CLASSIC') as ds:
        with open('nc_metadata.json') as data_file:
            attrs = json.load(data_file)
            for key in attrs:
                setattr(ds, key, attrs[key])
        setattr(ds, "date_created", datetime.now().strftime("%Y%m%dT%H%M%S"))
        
        rast = gdal.Open('NETCDF:"{}":lfmc_mean'.format(fmc_file))
        proj_wkt = rast.GetProjection()
        geot = rast.GetGeoTransform()
        
        t_dim = ds.createDimension("time", 1)
        x_dim = ds.createDimension("x", rast.RasterXSize)
        y_dim = ds.createDimension("y", rast.RasterYSize)

        var = ds.createVariable("time", "f8", ("time",))
        var.units = "seconds since 1970-01-01 00:00:00.0"
        var.calendar = "standard"
        var.long_name = "Time, unix time-stamp"
        var.standard_name = "time"
        var[:] = netCDF4.date2num([date], units="seconds since 1970-01-01 00:00:00.0", calendar="standard")

        var = ds.createVariable("x", "f8", ("x",))
        var.units = "m"
        var.long_name = "x coordinate of projection"
        var.standard_name = "projection_x_coordinate"
        var[:] = np.linspace(geot[0], geot[0]+(geot[1]*rast.RasterXSize), rast.RasterXSize)
        
        var = ds.createVariable("y", "f8", ("y",))
        var.units = "m"
        var.long_name = "y coordinate of projection"
        var.standard_name = "projection_y_coordinate"
        var[:] = np.linspace(geot[3], geot[3]+(geot[5]*rast.RasterYSize), rast.RasterYSize)
        
        var = ds.createVariable("flammability", 'f4', ("time", "y", "x"), fill_value=-9999.9)
        var.long_name = "Flammability Index"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = flam[None,...]
        
        var = ds.createVariable("anomaly", 'f4', ("time", "y", "x"), fill_value=-9999.9)
        var.long_name = "Flammability Anomaly"
        var.units = '%'
        var.grid_mapping = "sinusoidal"
        var[:] = anom[None,...]

        var = ds.createVariable("quality_mask", 'i1', ("time", "y", "x"), fill_value=0)
        var.long_name = "Combined Bands Quality Mask"
        var.units = 'Cat'
        var.grid_mapping = "sinusoidal"
        var[:] = q_mask.astype(np.int8)[None,...]

        var = ds.createVariable("sinusoidal", 'S1', ())
        var.grid_mapping_name = "sinusoidal"
        var.false_easting = 0.0
        var.false_northing = 0.0
        var.longitude_of_central_meridian = 0.0
        var.longitude_of_prime_meridian = 0.0
        var.semi_major_axis = 6371007.181
        var.inverse_flattening = 0.0
        var.spatial_ref = proj_wkt
        var.GeoTransform = "{} {} {} {} {} {} ".format(*[geot[i] for i in range(6)])

wgs84_wkt = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'

def pack_fmc_mosaic(date, fmc_mean, fmc_stdv, q_mask, dest):
    lat0 = -10.
    lat1 = -44.
    lon0 = 113.
    lon1 = 154.
    res = 0.005
    
    x_size = int((lon1 - lon0)/res)
    y_size = int((lat1 - lat0)/(-1*res))
    geot = [lon0, res, 0., lat0, 0., -1*res]
    
    with netCDF4.Dataset(dest, 'w', format='NETCDF4_CLASSIC') as ds:
        with open('nc_metadata.json') as data_file:
            attrs = json.load(data_file)
            for key in attrs:
                setattr(ds, key, attrs[key])
        setattr(ds, "date_created", datetime.now().strftime("%Y%m%dT%H%M%S"))
        
        t_dim = ds.createDimension("time", 1)
        x_dim = ds.createDimension("longitude", fmc_mean.shape[1])
        y_dim = ds.createDimension("latitude", fmc_mean.shape[0])

        var = ds.createVariable("time", "f8", ("time",))
        var.units = "seconds since 1970-01-01 00:00:00.0"
        var.calendar = "standard"
        var.long_name = "Time, unix time-stamp"
        var.standard_name = "time"
        var[:] = netCDF4.date2num([date], units="seconds since 1970-01-01 00:00:00.0", calendar="standard")

        var = ds.createVariable("longitude", "f8", ("longitude",))
        var.units = "degrees"
        var.long_name = "longitude"
        var.standard_name = "longitude"
        var[:] = np.linspace(lon0, lon1-res, num=x_size)
        
        var = ds.createVariable("latitude", "f8", ("latitude",))
        var.units = "degrees"
        var.long_name = "latitude"
        var.standard_name = "latitude"
        var[:] = np.linspace(lat0, lat1+res, num=y_size)
        
        var = ds.createVariable("fmc_mean", 'f4', ("time", "latitude", "longitude"), fill_value=-9999.9)
        var.long_name = "Mean Live Fuel Moisture Content"
        var.units = '%'
        var[:] = fmc_mean[None,...]

        var = ds.createVariable("fmc_stdv", 'f4', ("time", "latitude", "longitude"), fill_value=-9999.9)
        var.long_name = "Standard Deviation Live Fuel Moisture Content"
        var.units = '%'
        var[:] = fmc_stdv[None,...]
        
        var = ds.createVariable("quality_mask", 'i1', ("time", "latitude", "longitude"), fill_value=0)
        var.long_name = "Quality Mask"
        var.units = 'Cat'
        var[:] = q_mask[None,...]


def pack_flammability_mosaic(date, flam, anom, q_mask, dest):
    lat0 = -10.
    lat1 = -44.
    lon0 = 113.
    lon1 = 154.
    res = 0.005
    
    x_size = int((lon1 - lon0)/res)
    y_size = int((lat1 - lat0)/(-1*res))
    geot = [lon0, res, 0., lat0, 0., -1*res]
    
    with netCDF4.Dataset(dest, 'w', format='NETCDF4_CLASSIC') as ds:
        with open('nc_metadata.json') as data_file:
            attrs = json.load(data_file)
            for key in attrs:
                setattr(ds, key, attrs[key])
        setattr(ds, "date_created", datetime.now().strftime("%Y%m%dT%H%M%S"))
        
        t_dim = ds.createDimension("time", 1)
        x_dim = ds.createDimension("longitude", flam.shape[1])
        y_dim = ds.createDimension("latitude", flam.shape[0])

        var = ds.createVariable("time", "f8", ("time",))
        var.units = "seconds since 1970-01-01 00:00:00.0"
        var.calendar = "standard"
        var.long_name = "Time, unix time-stamp"
        var.standard_name = "time"
        var[:] = netCDF4.date2num([date], units="seconds since 1970-01-01 00:00:00.0", calendar="standard")

        var = ds.createVariable("longitude", "f8", ("longitude",))
        var.units = "degrees"
        var.long_name = "longitude"
        var.standard_name = "longitude"
        var[:] = np.linspace(lon0, lon1-res, num=x_size)
        
        var = ds.createVariable("latitude", "f8", ("latitude",))
        var.units = "degrees"
        var.long_name = "latitude"
        var.standard_name = "latitude"
        var[:] = np.linspace(lat0, lat1+res, num=y_size)
        
        var = ds.createVariable("flammability", 'f4', ("time", "latitude", "longitude"), fill_value=-9999.9)
        var.long_name = "Flammability Index"
        var.units = '%'
        var[:] = flam[None,...]

        var = ds.createVariable("anomaly", 'f4', ("time", "latitude", "longitude"), fill_value=-9999.9)
        var.long_name = "FMC Anomaly"
        var.units = '%'
        var[:] = anom[None,...]
        
        var = ds.createVariable("quality_mask", 'i1', ("time", "latitude", "longitude"), fill_value=0)
        var.long_name = "Quality Mask"
        var.units = 'Cat'
        var[:] = q_mask[None,...]
