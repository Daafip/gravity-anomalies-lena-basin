import rasterio
import matplotlib
import os
import matplotlib.pyplot as plt
import contextily as cx
from rasterio.plot import show as rioshow
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.fill import fillnodata
import numpy as np


# setting up some variables for nice plots
cmap = matplotlib.cm.get_cmap("tab10")


def get_background_map(name, bounds, dst_crs='EPSG:4326'):  
    """
    Given a name of area and bounds saves a local map file from Esri.WorldImagery using Contexitly
    """
    if dst_crs == 'EPSG:4326':
        fname = f'Figures\\{name}_wgs84.tif'
    else:
        fname = f'Figures\\{name}_{dst_crs[-4:]}.tif'
    if os.path.exists(fname):
        pass 
    else:
        # unpack the bounds into the wanted corners
        xmin, ymin, xmax, ymax = bounds

        # get the needed background map of the northsea with the bounds we found above and saves to given path
        _ = cx.bounds2raster(xmin, ymin, xmax, ymax,
                            ll=True,
                            path=f"Figures\\{name}.tif",
                            source=cx.providers.Esri.WorldImagery 
                            )

        # defines the wanted crs: we use WGS84                    
        # dst_crs = 'EPSG:4326'
        
        # the map is in RD coordinates, so we transform it here
        # also see https://rasterio.readthedocs.io/en/latest/topics/reproject.html 
        with rasterio.open(f'Figures\\{name}.tif') as src:
            transform, width, height = calculate_default_transform(src.crs, dst_crs, src.width, src.height, *src.bounds)
            kwargs = src.meta.copy()
            # change the dictionary to what we want
            kwargs.update({
                'crs': dst_crs,
                'transform': transform,
                'width': width,
                'height': height
            })

            # output a new .tif immage in the correct transformation        
            with rasterio.open(fname, 'w', **kwargs) as dst:
                # update the image with the projection we want 
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=dst_crs,
                        resampling=Resampling.nearest)
    return fname


def reproject_raster(source_path, destination_crs,ending="tif"):       
        dst_crs = destination_crs
        # the map is in local coordinates, so we transform it here
        # also see https://rasterio.readthedocs.io/en/latest/topics/reproject.html
        n = len(ending) + 1
        destination_transform_path = str(source_path)[:-n] + "_transform." + ending
        if os.path.exists(destination_transform_path):
            pass 
        else:
            with rasterio.open(source_path) as src:
                transform, width, height = calculate_default_transform(src.crs, dst_crs, src.width, src.height, *src.bounds)
                kwargs = src.meta.copy()
                # change the dictionary to what we want
                kwargs.update({
                    'crs': dst_crs,
                    'transform': transform,
                    'width': width,
                    'height': height
                })

                # output a new .tif immage in the correct transformation        
                with rasterio.open(destination_transform_path, 'w', **kwargs) as dst:
                    # update the image with the projection we want 
                    for i in range(1, src.count + 1):
                        reproject(
                            source=rasterio.band(src, i),
                            destination=rasterio.band(dst, i),
                            src_transform=src.transform,
                            src_crs=src.crs,
                            dst_transform=transform,
                            dst_crs=dst_crs,
                            resampling=Resampling.nearest)

        return destination_transform_path


def remove_below_0(path, ending="tif"):
    """removes any values below 0
    """
    source_raster_path = path
    n = len(ending) + 1
    destination_raster_path = path[:-n] + "_output." + ending
    with rasterio.open(source_raster_path, "r+") as src:
        src.nodata = 0 # set the nodata value
        profile = src.profile
        profile.update(
                dtype=rasterio.uint8,
                compress='lzw'
        )

        with rasterio.open(destination_raster_path, 'w',  **profile) as dst:
            for i in range(1, src.count + 1):
                band = src.read(i)
                band = np.where(band<0,0,band) # for completeness
                dst.write(band,i)

    return destination_raster_path
        




