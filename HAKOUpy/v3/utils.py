__author__ = 'nathansmith'

import os, sys
import numpy
from osgeo import gdal, osr, ogr
import json
from path import path

from compareResults import readDBinaryData, getPercentageVector, generateSurgePlotPoints, vectorChopper

try:
    import ipdb as pdb
except:
    import pdb


def convertLineListToShp(data, outfile):
    # set up shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    ds = driver.CreateDataSource(outfile)
    layer = ds.CreateLayer('surge', geom_type=ogr.wkbLineString)

    for i in range(len(data[0])):
        feat = ogr.Feature(layer.GetLayerDefn())
        feat.SetGeometryDirectly(_getLine([
                                            (data[0][i][0],data[1][i][0]), # from
                                            (data[0][i][1],data[1][i][1])  # to
                                          ]))
        layer.CreateFeature(feat)

    ds.Destroy()
    ds = None

def convertBinToShp(dataDir,binfile,outfile):
    """Get coord lists
    create geom
    write shp
    """
    l = generateSurgePlotPoints(dataDir, binfile)
    convertLineListToShp(l, outfile)

def convertProbabilityLevelToShp(dataDir,binfile,outfile,level):
    l = generateSurgePlotPoints(dataDir,-1,binfile,level)
    convertLineListToShp(l,outfile)

def _getLine(coords):
    line = ogr.Geometry(type=ogr.wkbLineString)
    for xy in coords:
        line.AddPoint_2D(xy[0],xy[1])
    return line

def getGeoPropertiesFromJSON(jsonfile):
    """
    """
    jsonconf = None
    with open(jsonfile) as f:
        jsonconf = json.load(f)

    Ny = int(jsonconf['lat_long_size'][0])
    Nx = int(jsonconf['lat_long_size'][1])
    geoproperties = {
        'tl': {
            'lat': max(jsonconf['lat'][0],jsonconf['lat'][-1]),
            'long': min(jsonconf['long'][0],jsonconf['long'][-1]),
        },
        'br': {
            'lat': min(jsonconf['lat'][0],jsonconf['lat'][-1]),
            'long': max(jsonconf['long'][0],jsonconf['long'][-1]),
        },
        'Ny': Ny,
        'Nx': Nx,
    }

    return geoproperties


def convertArrayToGeoTIFF(geoproperties, data, outfile):
    """
    """
    Nx = geoproperties['Nx']
    Ny = geoproperties['Ny']

    startPos = [geoproperties['tl']['long'],geoproperties['tl']['lat']]
    d_lat = (geoproperties['br']['lat'] - geoproperties['tl']['lat']) / (Ny - 1)
    d_long = (geoproperties['br']['long'] - geoproperties['tl']['long']) / (Nx - 1)

    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(outfile,Nx,Ny,1,gdal.GDT_Float32)
    ds.SetGeoTransform( [ startPos[0], d_long, 0, startPos[1], 0, d_lat ] )
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection( srs.ExportToWkt() )
    ds.GetRasterBand(1).WriteArray(data)
    ds = None

def _convertWaveListToNumpy(geoproperties, l):
    """
    """
    Nx = geoproperties['Nx']
    Ny = geoproperties['Ny']

    data = numpy.array(l,dtype=numpy.float)
    data = data.reshape((Nx,Ny))
    data = numpy.rot90(data)

    return data


def convertBinToGeoTIFF(geoproperties, datadir, binfile, outfile):
    """
    """
    rawtuple = vectorChopper(datadir, readDBinaryData(binfile)[0])[1]
    data = _convertWaveListToNumpy(geoproperties, rawtuple)
    convertArrayToGeoTIFF(geoproperties, data, outfile)

def convertProbabilityLevelToGeoTIFF(geoproperties, datadir, binfile, outfile, level):
    """
    """
    rawtuple = vectorChopper(datadir, getPercentageVector(binfile, level)[2])[1]
    data = _convertWaveListToNumpy(geoproperties, rawtuple)
    convertArrayToGeoTIFF(geoproperties, data, outfile)

def exportSurgeMapFile(outFile,shapeFile,projection,extents,mapName,fontLocation="font.list"):
    """Exports a MapFile for MapServer

    Args:
        | outFile (str): File to write to
        | shapeFile (str): ShapeFile with contours in EPSG:26905

    Kwargs:
        fontLocation (str): Location of fonts list for MapServer

    Returns:
        None
    """

    projection = projection.lower()
    abspath = path(outFile).abspath()
    shapefileAbsPath = path(shapeFile).abspath()
    mapFileContents = """
MAP
NAME '"""+mapName+"""'
STATUS ON
# Map image size
SIZE 512 256
MAXSIZE 4096
UNITS METERS
FONTSET '"""+fontLocation+"""'

EXTENT """+' '.join([str(e) for e in extents])+"""
PROJECTION
    #'proj=utm'
#'zone=5'
#'ellps=GRS80'
#'datum=NAD83'
#'units=m'
#'no_defs'
'init="""+projection+"""'
END

# Background color for the map canvas -- change as desired
#IMAGECOLOR 195 220 252
IMAGECOLOR 255 255 255
IMAGEQUALITY 95
IMAGETYPE png24

#CONFIG "MS_ERRORFILE" "/var/log/mapserver/ms_error.txt"
#DEBUG 1

#OUTPUTFORMAT
#  NAME png
#  DRIVER 'GD/PNG'
#  MIMETYPE 'image/png'
#  IMAGEMODE PC256
#  EXTENSION 'png'
#  TRANSPARENT ON
#END

OUTPUTFORMAT
NAME 'png24'
MIMETYPE 'image/png'
DRIVER 'AGG/PNG'
EXTENSION 'png'
IMAGEMODE RGBA
TRANSPARENT ON
END # OUTPUTFORMAT


# Legend
LEGEND
  IMAGECOLOR 255 255 255
STATUS ON
KEYSIZE 18 12
LABEL
  TYPE BITMAP
  SIZE MEDIUM
  COLOR 0 0 89
END
END

# Web interface definition. Only the template parameter
# is required to display a map. See MapServer documentation
WEB
# Set IMAGEPATH to the path where MapServer should
# write its output.
IMAGEPATH '/tmp/'

# Set IMAGEURL to the url that points to IMAGEPATH
# as defined in your web server configuration
IMAGEURL '/tmp/'

# WMS server settings
METADATA
  'wms_title'           '"""+mapName+"""'
  'wms_onlineresource'  'http://cybereye.crc.nd.edu/portal/mapserver/?map="""+abspath+"""'
  'wms_srs'             '"""+projection+""" epsg:4326 epsg:3857'
  'wms_enable_request'  '*'
END

#Scale range at which web interface will operate
# Template and header/footer settings
# Only the template parameter is required to display a map. See MapServer documentation
#    TEMPLATE 'fooOnlyForWMSGetFeatureInfo'
END



LAYER
    NAME '"""+mapName+"""'
    TYPE LINE
    DUMP true
    TEMPLATE fooOnlyForWMSGetFeatureInfo
    EXTENT """+' '.join([str(e) for e in extents])+"""
    DATA '"""+shapefileAbsPath+"""'
    METADATA
      'ows_title' '"""+mapName+"""'
    END
    STATUS ON
    TRANSPARENCY 100
    PROJECTION
    #'proj=utm'
    #'zone=5'
    #'ellps=GRS80'
    #'datum=NAD83'
    #'units=m'
    #'no_defs'
    'init="""+projection+"""'
    #'init=epsg:4326'
    END

    CLASS
        STYLE
            COLOR 0 0 0
        END
    END
END
END
"""
    with open(outFile,"w") as f:
        f.write(mapFileContents)


if __name__ == '__main__':
    datapath = '/vagrant/HAKOUv3/Hawaii_data/newData1.0/data/'
    resultspath = '/vagrant/HAKOUv3/KSM_CalculationZip_cmake/build/'
    # convertBinToGeoTIFF(getGeoPropertiesFromJSON(resultspath+'results.json'), datapath, resultspath+'result_deter.bin', 'test.tiff')
    # convertBinToShp(datapath, resultspath+'result_deter.bin', 'test.shp')
    # convertProbabilityLevelToGeoTIFF(getGeoPropertiesFromJSON(resultspath+'results.json'), datapath, resultspath+'result_prob.bin', 'test_25.tiff', 25)
    # convertProbabilityLevelToShp(datapath, resultspath+'result_prob.bin', 'test_25.shp', 25)