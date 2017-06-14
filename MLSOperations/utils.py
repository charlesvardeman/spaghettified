"""
.. module:: utils
.. moduleauthor:: Nathan Smith <nsmith13@nd.edu>

"""

import subprocess
from osgeo import gdal
from osgeo import ogr

#import pdb

def getShpProjection(input):
    """
    Finds the projection of the supplied shapefile
    """
    #pdb.set_trace()
    
    proj = []
    ds = ogr.Open(input)
    layer = ds.GetLayer()
    srs = layer.GetSpatialRef()
    srs.AutoIdentifyEPSG()
    proj.append(srs.GetAttrValue("PROJCS|GEOGCS|AUTHORITY",0))
    proj.append(srs.GetAttrValue("PROJCS|GEOGCS|AUTHORITY",1))
    return proj

def getShpExtents(input):
    """
    Calculates the extents of a ShapeFile
    """
    ds = ogr.Open(input)
    layer = ds.GetLayer()
    layerExt = layer.GetExtent()
    
    extent = [
        layerExt[0],
        layerExt[2],
        layerExt[1],
        layerExt[3],
    ]
    
    return extent

def getGeoTiffExtents(input):
    """
    Calculates the extents of a GeoTiff file.
    """
    ds = gdal.Open(input, gdal.GA_ReadOnly)
    (Xleft, deltaX, rotationX, Ytop, rotationY, deltaY) = ds.GetGeoTransform()
    Nx = ds.RasterXSize
    Ny = ds.RasterYSize
    
    Xright = Xleft + ( deltaX * Nx )
    Ybottom = Ytop + ( deltaY * Ny )
    
    extents = [0,0,0,0]
    
    if Xleft < Xright:
        extents[0] = Xleft
        extents[2] = Xright
    else:
        extents[0] = Xright
        extents[2] = Xleft
    
    if Ytop < Ybottom:
        extents[1] = Ytop
        extents[3] = Ybottom
    else:
        extents[1] = Ybottom
        extents[3] = Ytop
    
    return extents

def reprojectGeoTiff(fromProj, toProj, input, output, overwrite=False, noDataValue=None):
    """
    Reprojects a geotiff using gdalwarp
    """
    if noDataValue == None:
        noDataValue = ''
    else:
        noDataValue = '-srcnodata "%s" -dstnodata "%s"' % (noDataValue,noDataValue)
    
    param = ['gdalwarp','-overwrite' if overwrite else '','-s_srs',fromProj,'-t_srs',toProj,noDataValue,'"'+input+'"','"'+output+'"']
    cmd = ' '.join(param)
    #print cmd
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = ''.join(process.stdout.readlines())
    stderr = ''.join(process.stderr.readlines())
    
    if len(stderr) > 0:
        raise IOError(stderr)
    
    #print cmd
    #print stdout
    #print stderr

def reprojectShp(fromProj, toProj, input, output, overwrite=False):
    """
    Reprojects a shapefile using ogr2ogr
    """
    
    param = ['ogr2ogr','-overwrite' if overwrite else '','-s_srs',fromProj,'-t_srs',toProj,'"'+output+'"','"'+input+'"']
    cmd = ' '.join(param)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = ''.join(process.stdout.readlines())
    stderr = ''.join(process.stderr.readlines())
    
    if len(stderr) > 0:
        raise IOError(stderr)
        
    #print cmd
    #print stdout
    #print stderr

def tileGeoTiff(input, output):
    """
    Tiles a getTiff using gdal_translate
    """
    
    param = ['gdal_translate','-co TILED=YES','"'+input+'"','"'+output+'"']
    cmd = ' '.join(param)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = ''.join(process.stdout.readlines())
    stderr = ''.join(process.stderr.readlines())
    
    if len(stderr) > 0:
        raise IOError(stderr)
    

def createGeoTiffOverviews(input, levels):
    """
    Tiles a getTiff using gdaladdo
    
    levels should a string of ints ie '2 4 8 16 32'
    """
    
    param = ['gdaladdo','"'+input+'"', levels]
    cmd = ' '.join(param)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = ''.join(process.stdout.readlines())
    stderr = ''.join(process.stderr.readlines())
    
    if len(stderr) > 0:
        raise IOError(stderr)
    

def indexGeoTiff(input, output):
    """
    Indexes a geoTiff using gdaltindex (outputs shapefile)
    """
    
    param = ['gdaltindex','"'+output+'"','"'+input+'"']
    cmd = ' '.join(param)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = ''.join(process.stdout.readlines())
    stderr = ''.join(process.stderr.readlines())
    
    if len(stderr) > 0:
        raise IOError(stderr)
    

def tileShp(input):
    """
    Tiles a shapefile using Mapserver utility shptree if possible.  Outputs .qix file
    """
    
    param = ['shptree','"'+input+'"']
    cmd = ' '.join(param)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = ''.join(process.stdout.readlines())
    stderr = ''.join(process.stderr.readlines())
    
    if len(stderr) > 0:
        raise IOError(stderr)
    
