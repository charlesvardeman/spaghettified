"""
.. module:: MLSOperations
    :synopsis: Runs the MLSFitting and converts the output to GeoTIFF
    

.. moduleauthor:: Nathan Smith <nsmith13@nd.edu>

"""

# import MLSFitting

import numpy
from osgeo import gdal
from osgeo import osr
import csv
import os

def getStandardGeoProperties():
    """Returns the standard coordinates of the corners of the MLS output csv files.
    
    Currently the MLS models outputs two csv files, one close around the island of
    Oahu, and one wider around the surrounding area and touching other islands.
    Hsin refers to the close map and has the coordinates of the top-left and
    bottom-right corners.  Hsout refers to the wider map.  Here are the geoproperties
    for Oahu::
    
        geoproperties = {
            'Hsin': {
                'tl': {
                    'lat': 21.108,
                    'long': -158.584,
                },
                'br': {
                    'lat': 21.9,
                    'long': -157.392,
                }
            },
            'Hsout':{
                'tl': {
                    'lat': 20.33,
                    'long': -159.87,
                },
                'br': {
                    'lat': 22.7,
                    'long': -156.3,
                }
            },
        }
    
    Returns:
        geoproperties (dict)
    """
    geoproperties = {
        'Hsin': {
            'tl': {
                'lat': 21.108,
                'long': -158.584,
            },
            'br': {
                'lat': 21.9,
                'long': -157.392,
            }
        },
        'Hsout':{
            'tl': {
                'lat': 20.33,
                'long': -159.87,
            },
            'br': {
                'lat': 22.7,
                'long': -156.3,
            }
        },
    }
    return geoproperties

# def runMLS(num1,num2,num3,num4,num5,dataDir,outputDir,outfile1,outfile2):
#     """Runs MLSFitting from the extension module.
#     
#     Outputs two csv files to the output directory.
#     
#     Args:
#         | num1-5 (float): Input vector for MLS fitting model
#         | dataDir (str): Where the input data files are located
#         | outputDir (str): Where to output the CSV files
#         | outfile1 (str): Name of first (closer) CSV file
#         | outfile2 (str): Name of second (wider) CSV file
#     
#     Returns:
#         None
#     """
#     
#     MLSFitting.MLSFitting(num1,num2,num3,num4,num5,dataDir,outputDir,outfile1,outfile2)

def convertCSVToGeoTIFF(geoproperties, datafile, outputdir, outfile):
    """Converts MLSFitting CSV file to GeoTIFF
    
    Args:
        | geoproperties (dict): see below
        | datafile (str): csv file to be converted
        | outputdir (str): directory to put GeoTIFF
        | outfile (str): name of GeoTIFF
    
    Takes a set of geographic properties about the incoming csv files such as
    either Hsin or Hsout subdicts as described in :py:meth:`MLSOperations.getStandardGeoProperties`
    
    It will import the csv into a numpy array, check if the last row/column are
    empty and slice them off, then output a geotiff in the specified directory
    with the specified file name.
    
    Returns:
        None
    """
    
    data = numpy.genfromtxt(datafile, delimiter=',')
    
    Ny, Nx = data.shape
    
    #print data.size
    #print data.shape
    #print data
    
    # slice out the last column of null values
    if str(data[Ny-1][Nx-1]) == 'nan':
        data = data[:,:-1]
    
    Ny, Nx = data.shape
    
    #print data.size
    #print data.shape
    #print data
    #print Ny, Nx
    
    startPos = [geoproperties['tl']['long'],geoproperties['tl']['lat']]
    d_lat = (geoproperties['br']['lat'] - geoproperties['tl']['lat']) / (Ny - 1)
    d_long = (geoproperties['br']['long'] - geoproperties['tl']['long']) / (Nx - 1)
    
    #print startPos, d_lat, d_long
    
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(os.path.join(outputdir,outfile),Nx,Ny,1,gdal.GDT_Float32)
    #ds = driver.Create('output/output.tif',Nx,Ny,1,gdal.GDT_Byte)
    #ds.SetGeoTransform( [ -158.584, .008, 0, 21.108, 0, .008 ] )
    ds.SetGeoTransform( [ startPos[0], d_long, 0, startPos[1], 0, d_lat ] )
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection( srs.ExportToWkt() )
    ds.GetRasterBand(1).WriteArray(data)
    ds = None