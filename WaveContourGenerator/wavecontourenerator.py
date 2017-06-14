"""
.. module:: WaveContourGenerator
    :synopsis: Generates a contour map from a GeoTIFF
    

.. moduleauthor:: Nathan Smith <nsmith13@nd.edu>

"""

import os
import numpy as np
from osgeo import gdal
from osgeo import osr
from osgeo import ogr
import csv
from path import path
import colorsys

import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

try:
    import ipdb as pdb
except:
    import pdb


class WaveContourGenerator:
    
    def readGeoTIFF(self, dataFile):
        """Creates the data matrix from a GeoTIFF file.
        
        Args:
            dataFile (str): GeoTIFF file to read data from.
        
        Returns:
            None
        """
        
        ds = gdal.Open(dataFile, gdal.GA_ReadOnly)

        (Xleft, deltaX, rotationX, Ytop, rotationY, deltaY) = ds.GetGeoTransform()
        srs_wkt = ds.GetProjection()
        Nx = ds.RasterXSize
        Ny = ds.RasterYSize

        ## for a single-band GeoTIFF
        self.data = ds.GetRasterBand(1).ReadAsArray()
        
        self.data_lat = np.array( self._generateCoordList(Ytop,Ny,deltaY) )
        self.data_long = np.array( self._generateCoordList(Xleft,Nx,deltaX) )

    def _generateCoordList(self,start,n,delta):
        l = []
        for i in range(n):
            l.append(start + (delta*i))
        return l

    def readCSV(self,dataFile,dataLat,dataLong):
        """Creates the data matrix from a GeoTIFF file.
        
        Args:
            | dataFile (str): File to read data from
            | dataLat (str): File to read latitudes for each row
            | dataLong (str): File to read longitudes for each column
        
        Returns:
            None
        """
        
        self.data = np.genfromtxt(dataFile, delimiter=',')
        self.data_lat = np.genfromtxt(dataLat, delimiter=',')
        self.data_long = np.genfromtxt(dataLong, delimiter=',')

    def showRawData(self):
        """Displays the raw data using matplotlib.
        
        Must not be using Agg as renderer.
        """
        plt.imshow(np.flipud(self.data), extent=[self.data_long[0],self.data_long[len(self.data_long)-1],self.data_lat[0],self.data_lat[len(self.data_lat)-1]])
        plt.show()

    def createPlot(self, numLevels=None, showPlot=False, outputPDF=False, outfile=None, maxLevel=20, stepLevel=1):
        """Creates the contour plot using MatPlotLib's contourf module.
        
        To show a plot, use a renderer other than Agg.
        
        Kwargs:
            | numLevels (int): Number of levels to compute
            | showPlot (bool): Whether to display the plot with mpl
            | outputPDF (bool): Whether to output a pdf of plot
            | outfile (str): Where to save the pdf file
            | maxLevels (int): Max level value
            | stepLevels (int): Level step size
        
        Returns:
            None
        """
        self.numLevels = numLevels
        self.maxLevel = maxLevel
        self.stepLevel = stepLevel
        
        x,y = np.meshgrid(self.data_long, self.data_lat)
        
        #self.levels = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0]
        #self.levels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
        self._calculateLevels()
        
        self.contour = plt.contourf(x,y,self.data,levels=self.levels)
        
        if showPlot == True:
            plt.show()
        if outputPDF == True and outfile != None:
            plt.savefig(outfile)

    def _calculateLevels(self):
        """Calculates a list of the levels to use for contourf
        """
        min = 0.001
        
        if self.numLevels != None:
            max = np.max(self.data)
            min = 0.1
            levels = []
            
            delta = (max - min) / self.numLevels
            
            for i in range(int(self.numLevels)):
                levels.append(min + (delta * i))
            levels.append(max)
            
            self.levels = levels
        
        # 2.5m intervals up to max
        else:
            self.numLevels = int(self.maxLevel/self.stepLevel)
            self.levels = [self.stepLevel*(i+1) for i in range(self.numLevels)]
            self.levels.append(np.max(self.data))
            self.levels.insert(0,min)
            #self.levels.append(100)
        
        #print self.levels

    def _calculateColors(self):
        """
        This will create a rainbow using hsv and convert to rgb for each color
        to go in the mapfile

        0deg is red
        240deg is blue

        Go from 240 -> 0
        """
        colorTable = []

        colormin = 0
        colormax = 240
        
        numLevelsMax = len(self.levels) - 1 # to get rid of unnecessary extra level at end

        delta = (colormax - colormin) / (numLevelsMax - 1.0)

        for i in range(int(numLevelsMax)):
            colorValue = {}
            colorValue['value'] = self.levels[i]
            colorhue = (colormax - (delta * i)) / 360
            r, g, b = colorsys.hsv_to_rgb(colorhue,1,.95)
            colorValue['r'] = int(r * 255)
            colorValue['g'] = int(g * 255)
            colorValue['b'] = int(b * 255)
            colorTable.append(colorValue)

        self.colorTable = colorTable

    def exportShp(self, outFile):
        """Export a shapefile of contour levels to specified location
        
        Args:
            outFile (str): File to write to
        
        Returns:
            None
        """
        # Generate collection for shapefile
        collection = self.contour.collections
        
        label = [[] for ii in range(len(self.levels))]
        for i in range(0, len(self.levels)):
            label[i] = str(self.levels[i]) +'-' + str(self.levels[i]+5)
        else:
            label[i] = str(self.levels[i]) +'+'
        
        driver = ogr.GetDriverByName('ESRI Shapefile')
        
        if os.path.exists(outFile):
            driver.DeleteDataSource(outFile)
        ds = driver.CreateDataSource(outFile)
        self.layer = ds.CreateLayer('test', geom_type=ogr.wkbPolygon)
        
        fieldDefn = ogr.FieldDefn('HEIGHT', ogr.OFTReal)
        self.layer.CreateField(fieldDefn)
        
        fieldDefn = ogr.FieldDefn('HEIGHT_STR', ogr.OFTString)
        self.layer.CreateField(fieldDefn)
        
        fieldDefn    = ogr.FieldDefn('POLY_ID',ogr.OFTInteger)
        self.layer.CreateField(fieldDefn)
        
        #print "In main, range(len(collection))=", range(len(collection))
        cnLevel_list = map(self._proc_cnLevel,zip(collection,range(len(collection))))
        ds.Destroy()

    def exportMapFile(self,outFile,shapeFile,projection,extents,mapName,fontLocation="font.list"):
        """Exports a MapFile for MapServer
        
        Args:
            | outFile (str): File to write to
            | shapeFile (str): ShapeFile with contours in EPSG:26905
        
        Kwargs:
            fontLocation (str): Location of fonts list for MapServer
        
        Returns:
            None
        """
        
        self._calculateColors()
        
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
      'wms_onlineresource'  'https://cybereye.crc.nd.edu/portal/mapserver/?map="""+abspath+"""'
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
    TYPE POLYGON
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
    
    LABELITEM 'HEIGHT_STR'
    #LABELMAXSCALEDENOM 1451339.3664
    
    CLASSITEM 'HEIGHT_STR'
    """

        for i in range(len(self.colorTable)):
            valueclass = """
        CLASS
            EXPRESSION ('[HEIGHT_STR]' eq '"""+str(self.colorTable[i]['value'])+"""')
            LABEL
                ANGLE FOLLOW
                MAXOVERLAPANGLE 20
                FONT 'arial'
                MAXSIZE 10
                MINSIZE 10
                SIZE 10
                BUFFER 2
                COLOR 255 255 255
                MINDISTANCE 75
                MINFEATURESIZE 20
                OFFSET 0 0
                OUTLINECOLOR 0 0 0
                PARTIALS FALSE
                POSITION CC
                #SHADOWSIZE 1 1
                TYPE TRUETYPE
            END # LABEL
            STYLE
                COLOR """+str(self.colorTable[i]['r'])+" "+str(self.colorTable[i]['g'])+" "+str(self.colorTable[i]['b'])+"""
                WIDTH 0.91
                OUTLINECOLOR 0 0 0
            END
        END
        """
            mapFileContents += valueclass
            # end loop

        mapFileContents += """
  END
END
"""
        with open(outFile,"w") as f:
            f.write(mapFileContents)



    def _proc_ring(self, polygon):
        ring = ogr.Geometry(ogr.wkbLinearRing)
        for (x,y) in polygon: ring.AddPoint(x,y)
        #print "NVERT ",ring.GetPointCount()

        # find mean distance between vertices
        #print "Polygon: ", polygon
        dst = np.zeros((ring.GetPointCount()))
        for n in range(0,ring.GetPointCount()-1):
            (x1,y1) = polygon[n]
            (x2,y2) = polygon[n+1]
            dst[n] = np.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
        #print "Dist :", dst

        # Count non-zero distances and compute mean among them
        meandst, n = 0.0, 0.0
        for rho in dst:
            if rho > 0.0:
                meandst, n = meandst+rho, n+1.0
        meandst = meandst/n
        #print "MDIST ", meandst
        return ring

    def _proc_multipolygon(self, multiArgs):
        multipolygon, count, inner_count = multiArgs
        #print "In _proc_multipolygon:         count =", count
        #print "                    :   inner_count =", inner_count
        polygon = ogr.Geometry(ogr.wkbPolygon)
        multipolygon.should_simplify = False
        polyList = multipolygon.to_polygons()
        ringList = map(self._proc_ring,polyList)
        for ring in ringList: polygon.AddGeometry(ring)
        feat = ogr.Feature( self.layer.GetLayerDefn())
        feat.SetGeometry(polygon)
        feat.SetField( 'HEIGHT', self.levels[count] )
        feat.SetField( 'HEIGHT_STR', str(self.levels[count]) )
        #print "                    : contour_level =", self.levels[count], "\n"
        feat.SetField('POLY_ID',str(count))
        if self.layer.CreateFeature(feat) != 0:
            #print "Failed to create feature in shapefile.\n"
            sys.exit( 1 )
        feat.Destroy()

    def _proc_cnLevel(self, cnArgs):
        cnLevel, count = cnArgs
        #print "In _proc_cnLevel: cnLevel, count =",cnLevel," ",count,"\n"
        multipolygons = cnLevel.get_paths()
        clist = [count for i in range(len(multipolygons))]
        return map(self._proc_multipolygon,zip(multipolygons,clist,range(len(multipolygons))))
