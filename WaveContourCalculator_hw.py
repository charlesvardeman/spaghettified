### DKDK path for this class (and functions), WaveContourGenerator: trunk/WaveContourGenerator/WaveContourGenerator/__init__.py
### DKDK below functions with the suffix, _hw, are modified by DK for hurricane wind model
### DKDK thus those functions should be integrated with the corresponding WCC py file: here, to __init__.py

import os
import numpy as np, numpy
#import numpy as numpy

from path import path
import colorsys


from osgeo import gdal
from osgeo import osr
from osgeo import ogr

'''
import csv
from path import path
import colorsys
'''

import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# DKDK add these
import json


class WaveContourGenerator:
    """Defining functions that are utilized in the hurricane wind model such as wind speed contour and wind vector fields.

    * Suffix of the functions::

       _hw()    : hurricane wind speed contour
       _hwvec() : hurricane wind vector fields

    """
    def readGeoTIFF_hw(self, dataFile):
        """Creates the data matrix from a GeoTIFF file (hurricane wind speed contour), but this is not used anymore.

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
        self.data_GTIF_hw = ds.GetRasterBand(1).ReadAsArray()

        self.data_lat_GTIF_hw = np.array( self._generateCoordList(Ytop,Ny,deltaY) )
        self.data_long_GTIF_hw = np.array( self._generateCoordList(Xleft,Nx,deltaX) )


    def createPlot_hw(self, showPlot=False, outputPDF=False, outfile=None, maxLevel=15, stepLevel=2.5):
        """Creates the wind speed contour plot using MatPlotLib's contourf module.

        Note that this is not used anymore in the wind model. Instead, "createPlot_hw_short" is employed.

        Kwargs:
            | showPlot (bool): Whether to display the plot with mpl
            | outputPDF (bool): Whether to output a pdf of plot
            | outfile (str): Where to save the pdf file
            | maxLevels (int): Max level value
            | stepLevels (int): Level step size

        Returns:
            None
        """

        self.stepLevel_hw = stepLevel
        self.maxLevel_hw = np.amax(self.data_GTIF_hw)+2.5
        print np.amax(self.data_GTIF_hw), self.maxLevel_hw



        x,y = np.meshgrid(self.data_long_GTIF_hw, self.data_lat_GTIF_hw)

        #self.levels = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0]
        #self.levels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]

        #self.levels_hw=[5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40]
        #self.levels_hw=[5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30]
        self._calculateLevels_hw()

        print self.levels_hw

        self.contour_hw = plt.contourf(x,y,self.data_GTIF_hw,levels=self.levels_hw)
        #plt.show()

        if showPlot == True:
            plt.show()
        if outputPDF == True and outfile != None:
            plt.savefig(outfile)


    def _calculateLevels_hw(self):
        """Calculates a list of the wind speed levels to use for matplotlib's contourf
        """
        #min = 0.001
        # DKDK 2.5m/s intervals up to max as the max value was set to datamax+2.5

        self.numLevels_hw = int(self.maxLevel_hw/self.stepLevel_hw)
        self.levels_hw = [self.stepLevel_hw*(i+2) for i in range(self.numLevels_hw-1)]
        #self.levels_hw.append(np.max(self.data))
        #self.levels.insert(0,min)
        #self.levels.append(100)


    def exportShp_hw(self, outFile):
        """Export a shapefile of wind speed contour levels to specified location

        Args:
            outFile (str): File to write to

        Returns:
            None
        """
        # Generate collection for shapefile
        collection_hw = self.contour_hw.collections

        label = [[] for ii in range(len(self.levels_hw))]
        for i in range(0, len(self.levels_hw)):
            label[i] = str(self.levels_hw[i]) +'-' + str(self.levels_hw[i]+5)
        else:
            label[i] = str(self.levels_hw[i]) +'+'

        driver = ogr.GetDriverByName('ESRI Shapefile')

        if os.path.exists(outFile):
            driver.DeleteDataSource(outFile)
        ds = driver.CreateDataSource(outFile)
        self.layer_hw = ds.CreateLayer('test_hw', geom_type=ogr.wkbPolygon)

        fieldDefn = ogr.FieldDefn('HWS', ogr.OFTReal)
        self.layer_hw.CreateField(fieldDefn)

        fieldDefn = ogr.FieldDefn('HWS_STR', ogr.OFTString)
        self.layer_hw.CreateField(fieldDefn)

        fieldDefn    = ogr.FieldDefn('HWS_Pol_ID',ogr.OFTInteger)
        self.layer_hw.CreateField(fieldDefn)

        #print "In main, range(len(collection))=", range(len(collection))
        cnLevel_list_hw = map(self._proc_cnLevel_hw,zip(collection_hw,range(len(collection_hw))))
        print cnLevel_list_hw

        ds.Destroy()


    def _proc_ring_hw(self, polygon):
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

    def _proc_multipolygon_hw(self, multiArgs):
        multipolygon, count, inner_count = multiArgs
        #print "In _proc_multipolygon:         count =", count
        #print "                    :   inner_count =", inner_count
        polygon = ogr.Geometry(ogr.wkbPolygon)
        multipolygon.should_simplify = False
        polyList = multipolygon.to_polygons()
        ringList = map(self._proc_ring_hw,polyList)
        for ring in ringList: polygon.AddGeometry(ring)
        feat = ogr.Feature( self.layer_hw.GetLayerDefn())
        feat.SetGeometry(polygon)
        feat.SetField( 'HWS', self.levels_hw[count] )
        feat.SetField( 'HWS_STR', str(self.levels_hw[count]) )
        #print "                    : contour_level =", self.levels[count], "\n"
        feat.SetField('HWS_Pol_ID',str(count))
        if self.layer_hw.CreateFeature(feat) != 0:
            #print "Failed to create feature in shapefile.\n"
            sys.exit( 1 )

        feat.Destroy()


    def _proc_cnLevel_hw(self, cnArgs):
        cnLevel, count = cnArgs
        #print "In _proc_cnLevel: cnLevel, count =",cnLevel," ",count,"\n"
        multipolygons = cnLevel.get_paths()
        clist = [count for i in range(len(multipolygons))]
        return map(self._proc_multipolygon_hw,zip(multipolygons,clist,range(len(multipolygons))))


    ### DKDK add one more parameter, LOCAL_MAPSERVER_URL - nono
    def exportMapFile_hw(self,outFile,shapeFile,projection,extents,mapName,fontLocation):
    #def exportMapFile_hw(self,outFile,shapeFile,projection,extents,mapName,fontLocation, def_mapurl):
    #def exportMapFile_hw(self,outFile,shapeFile,projection,extents,mapName,fontLocation="font.list"):

        """Exports a MapFile for MapServer - hurricane wind speed contour

        Args:
            | outFile (str): File to write to
            | shapeFile (str): ShapeFile with wind speed contour

        Kwargs:
            fontLocation (str): Location of fonts list for MapServer

        Returns:
            None
        """

        self._calculateColors_hw()

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
    ### DKDK RGBA indicates transparent background - http://mapserver.org/output/agg.html
    IMAGEMODE RGBA
    ### DKDK not quite sure whether the option TRANSPARENT works
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
    ### DKDK wms_onlineresource may need to be changed for server type (demo, dev, staging)
    METADATA
      'wms_title'           '"""+mapName+"""'

      ### DKDK change below to use parameter - nono it seens that below is not actually used
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
    TYPE POLYGON
    DUMP true
    TEMPLATE fooOnlyForWMSGetFeatureInfo
    EXTENT """+' '.join([str(e) for e in extents])+"""
    DATA '"""+shapefileAbsPath+"""'
    METADATA
      'ows_title' '"""+mapName+"""'
    END
    STATUS ON
    ### DKDK TRANSPARENCY is being deprecated (0 = transparent) - OPACITY is recommended instead
    ### DKDK 0 = nothing is shown (at least in an image viewer), 1 = looks a bit transparent than other number
    ### DKDK via mapserver, OPACITY effect is more clear than an image viewer
    ###TRANSPARENCY 100
    OPACITY 100
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

    LABELITEM 'HWS_STR'
    #LABELMAXSCALEDENOM 1451339.3664

    CLASSITEM 'HWS_STR'
    """

        for i in range(len(self.colorTable_hw)):
            valueclass = """
        CLASS
            EXPRESSION ('[HWS_STR]' eq '"""+str(self.colorTable_hw[i]['value'])+"""')
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
                COLOR """+str(self.colorTable_hw[i]['r'])+" "+str(self.colorTable_hw[i]['g'])+" "+str(self.colorTable_hw[i]['b'])+"""
                WIDTH 0.91
                ### DKDK if OUTLINECOLOR is blocked or not used/defined in the STYLE, then no outline will be shown - this may be better in contourf
                #OUTLINECOLOR 0 0 0
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


    def _calculateColors_hw(self):
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

        numLevelsMax = len(self.levels_hw) - 1 # to get rid of unnecessary extra level at end

        delta = (colormax - colormin) / (numLevelsMax - 1.0)

        for i in range(int(numLevelsMax)):
            colorValue = {}
            colorValue['value'] = self.levels_hw[i]
            colorhue = (colormax - (delta * i)) / 360
            r, g, b = colorsys.hsv_to_rgb(colorhue,1,.95)
            colorValue['r'] = int(r * 255)
            colorValue['g'] = int(g * 255)
            colorValue['b'] = int(b * 255)
            colorTable.append(colorValue)

        self.colorTable_hw = colorTable

    ### DKDK new function for the test of shortening the procedure is introduced
    ### DKDK this function will do the following things: read json, array conversion, calculate contourf level (level_hw), contourf via matplotlib

    def createPlot_hw_short(self,jsonfile,showPlot=False, outputPDF=False, outfile=None, maxLevel=15, stepLevel=2.5):
        """Read json file generated from the hwindmodel run and Creates the wind speed contour plot using MatPlotLib's contourf module.

        This is employed for the wind speed contour to shorten the overall process.

        Kwargs:
            | showPlot (bool): Whether to display the plot with mpl
            | outputPDF (bool): Whether to output a pdf of plot
            | outfile (str): Where to save the pdf file
            | maxLevels (int): Max level value
            | stepLevels (int): Level step size

        Returns:
            None
        """

        jsonconf=None
        #with open('results_hwind.json') as f:
        with open(jsonfile) as f:
            jsonconf = json.load(f)
        f.close()

        print jsonfile

        hwdata=jsonconf['output_array_hw']
        thetadata=jsonconf['theta_hw']
        size_one_variable_hw=jsonconf['size_each_variable_and_total_number5']
        size_output_array_hw=jsonconf['size_output_array_hw']
        # size_each_variable_and_total_number5 = [37.0, 14.0, 5.0]

        # DKDK single line like below works too
        hwdata_reshape=numpy.transpose(numpy.reshape(hwdata,(int(size_one_variable_hw[1])*int(size_one_variable_hw[2]),int(size_one_variable_hw[0]))))

        #print hwdata_reshape.shape, len(numpy.transpose(hwdata_reshape))

        # here, 5 = total number of variables, x,y,u,v,total
        num_variables=int(size_one_variable_hw[2])
        lat_hw   =hwdata_reshape[:,numpy.arange(0,len(numpy.transpose(hwdata_reshape)),num_variables)]
        long_hw  =hwdata_reshape[:,numpy.arange(1,len(numpy.transpose(hwdata_reshape)),num_variables)]
        u_hw     =hwdata_reshape[:,numpy.arange(2,len(numpy.transpose(hwdata_reshape)),num_variables)]
        v_hw     =hwdata_reshape[:,numpy.arange(3,len(numpy.transpose(hwdata_reshape)),num_variables)]
        wspeed_hw=hwdata_reshape[:,numpy.arange(4,len(numpy.transpose(hwdata_reshape)),num_variables)]

        print lat_hw.shape, long_hw.shape,u_hw.shape,v_hw.shape,wspeed_hw.shape
        print numpy.amax(lat_hw),numpy.amax(long_hw),numpy.amax(wspeed_hw)
        print numpy.amin(lat_hw),numpy.amin(long_hw),numpy.amin(wspeed_hw)

        Numx=int(size_one_variable_hw[0])
        Numy=int(size_one_variable_hw[1])

        ###
        ### DKDK - for test purpose that can be blocked later - start
        #long_hw_reshp  =numpy.reshape(long_hw,Numx*Numy)
        #lat_hw_reshp   =numpy.reshape(lat_hw,Numx*Numy)
        #wspeed_hw_reshp=numpy.reshape(wspeed_hw,Numx*Numy)
        #print len(long_hw_reshp), long_hw_reshp.shape
        #X1, Y1 = numpy.meshgrid(long_hw_reshp, lat_hw_reshp)
        #V1 = griddata((long_hw_reshp, lat_hw_reshp), wspeed_hw_reshp, (X1, Y1), method='nearest')
        #plt.show()
        ### DKDK - for test purpose that can be blocked later  - end
        ###

        xv=numpy.reshape(long_hw,Numx*Numy)
        yv=numpy.reshape(lat_hw,Numx*Numy)
        values=numpy.reshape(wspeed_hw,Numx*Numy)

        #startPos = [geoproperties['tl']['long'],geoproperties['tl']['lat']]
        #d_lat = (geoproperties['br']['lat'] - geoproperties['tl']['lat']) / (Ny - 1)
        #d_long = (geoproperties['br']['long'] - geoproperties['tl']['long']) / (Nx - 1)

        #startPos = [numpy.amin(X1),numpy.amax(Y1)]
        #print startPos, len(Y1)
        #d_lat = (numpy.amin(Y1) - numpy.amax(Y1)) / (len(Y1) - 1)
        #d_long = (numpy.amax(X1) - numpy.amin(X1)) / (len(X1) - 1)

        ### DKDK not quite sure geoproperties_hw is necessary
        geoproperties_hw = {
            'tl': {
                'lat': numpy.amax(yv),
                'long': numpy.amin(xv),
            },
            'br': {
                'lat': numpy.amin(yv),
                'long': numpy.amax(xv),
            },
            'Ny': len(xv),
            'Nx': len(yv),
        }
        print geoproperties_hw

        # #
        # # Setting the default values for GeoTIFF
        # #
        # #geotransform = None
        # #proj='osr.SpatialReference()'
        # #proj=None
        # xSize=100
        # ySize=100
        # power=5
        # smoothing=0
        # #zField='Z'
        # #dataFile=None
        # driverName='GTiff'
        # #outFile=hwindgeotiff  # DKDK this is now received via function call

        # xMin=min(xv)
        # xMax=max(xv)
        # yMin=min(yv)
        # yMax=max(yv)
        # #power=1
        # #smoothing=20

        # print xMin, xMax, yMin, yMax
        # print (xMax-xMin)/xSize
        # #print min(xv)

        # #geotransform=[]
        # #geotransform.append(xMin)
        # #geotransform.append((xMax-xMin)/xSize)
        # #geotransform.append(0)
        # #geotransform.append(yMax)
        # #geotransform.append(0)
        # #geotransform.append((yMin-yMax)/ySize)

        # geotransform=[xMin,(xMax-xMin)/xSize,0,yMax,0,(yMin-yMax)/ySize]
        # #geotransform=[10,0.6,0,90,0,0.8]
        # print geotransform

        # #Creating the interpolation function and populating the output matrix value
        # data_hw = invDist_hw(xv,yv,values,geotransform,xSize,ySize,power,smoothing,driverName,outFile)

        # return (geoproperties_hw, data_hw)

        """
        self.stepLevel_hw = stepLevel
        self.maxLevel_hw = np.amax(self.data_GTIF_hw)+2.5
        print np.amax(self.data_GTIF_hw), self.maxLevel_hw



        x,y = np.meshgrid(self.data_long_GTIF_hw, self.data_lat_GTIF_hw)

        #self.levels = [1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0]
        #self.levels = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]

        #self.levels_hw=[5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5, 40]
        #self.levels_hw=[5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30]
        self._calculateLevels_hw()

        print self.levels_hw

        self.contour_hw = plt.contourf(x,y,self.data_GTIF_hw,levels=self.levels_hw)
        #plt.show()

        if showPlot == True:
            plt.show()
        if outputPDF == True and outfile != None:
            plt.savefig(outfile)

        """

        self.stepLevel_hw = stepLevel
        self.maxLevel_hw = np.amax(values)+2.5
        print np.amax(values), self.maxLevel_hw


        #self.levels_hw=[5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30]
        self._calculateLevels_hw()
        print self.levels_hw


        self.contour_hw = plt.contourf(long_hw,lat_hw,wspeed_hw,levels=self.levels_hw)
        #self.contour_hw = plt.contourf(long_hw,lat_hw,wspeed_hw,levels=self.levels_hw,edgecolors='none')
        #self.contour_hw = plt.contourf(long_hw,lat_hw,wspeed_hw,levels=self.levels_hw,shading='faceted')


        #plt.clabel(self.contour, fmt = '%.1f', colors = 'black', fontsize=14, inline=1)
        #plt.show()
        plt.savefig('test1.png')

        #collection_hw = self.contour_hw.collections
        #print collection_hw, len(collection_hw)
        # DKDK len(collection_hw) = num of self.levels_hw !!!

        #if showPlot == True:
        #    plt.show()
        #if outputPDF == True and outfile != None:
        #    plt.savefig(outfile)

    #
    # DKDK Below Functions from here are the same with pre-existing ones at FSWP (i.e., no changes)
    #      Thus, there is no need to add these again


    #  DKDK_generateCoordList is reusable, so does not need to change
    def _generateCoordList(self,start,n,delta):
        l = []
        for i in range(n):
            l.append(start + (delta*i))
        return l


    ###
    ### DKDK Below functions are for hurricane wind vector fields
    ### DKDK 1. Please note that invDist_hwvec() and pointValue_hwvec() should be defined ahead of FromJsonToGeoTIFF_hwvec()
    ### DKDK 2. Since FromJsonToGeoTIFF_hwvec(), invDist_hwvec(), and pointValue_hwvec() are now inside the class WaveContourGenerator, 'self' needs to be used
    ### DKDK    - FromJsonToGeoTIFF_hwvec(self,...), invDist_hwvec(self,...), and pointValue_hwvec(self,...)
    ### DKDK 3. Similar to the No. 2, when calling invDist_hwvec() and pointValue_hwvec(), it should be called like below inside the parent function
    ### DKDK    - self.invDist_hwvec() and self.pointValue_hwvec()
    ###

    ### DKDK new function for interpolation of vector fields
    def invDist_hwvec(self, xv,yv,u_vec,v_vec,geotransform,xSize,ySize,power,smoothing,driverName,outFile):
        """Used for wind vector field -
        Create two bands of GeoTIFF (for u & v vectors separately) involving rasterization and interpolation

        """
        #Transform geographic coordinates to pixels
        xv1=numpy.zeros(len(xv))
        yv1=numpy.zeros(len(yv))

        for i in range(0,len(xv)):
             xv1[i] = (xv[i]-geotransform[0])/geotransform[1]
        for i in range(0,len(yv)):
             yv1[i] = (yv[i]-geotransform[3])/geotransform[5]

        #Creating the file
        driver = gdal.GetDriverByName( driverName )
        ds = driver.Create( outFile, xSize, ySize, 2, gdal.GDT_Float32)

        #if proj is not None:
        #    ds.SetProjection(proj.ExportToWkt())
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        ds.SetProjection(srs.ExportToWkt())

        ds.SetGeoTransform(geotransform)
        #srs = osr.SpatialReference()
        #srs.ImportFromEPSG(4326)
        #ds.SetProjection(srs.ExportToWkt())

        #valuesGrid = numpy.zeros((ySize,xSize))
        u = numpy.zeros((ySize,xSize))
        v = numpy.zeros((ySize,xSize))

        #Getting the interpolated values
        """
        ### DKDK using pointValue_hw() will lead to poor performance than pointValue_hwvec()
        for x in range(0,xSize):
            for y in range(0,ySize):
                u[y,x] = pointValue_hw(x,y,power,smoothing,xv1,yv1,u_vec)
                v[y,x] = pointValue_hw(x,y,power,smoothing,xv1,yv1,v_vec)

        ### DKDK pointValue_hwvec() is faster
        ut = numpy.zeros((ySize,xSize))
        vt = numpy.zeros((ySize,xSize))
        for x in range(0,xSize):
            for y in range(0,ySize):
                ut[y,x] = pointValue_hwvec(x,y,power,smoothing,xv1,yv1,u_vec)
                vt[y,x] = pointValue_hwvec(x,y,power,smoothing,xv1,yv1,v_vec)
        """

        for x in range(0,xSize):
            for y in range(0,ySize):
                ### DKDK - since pointValue_hwvec() is now a part of WCC_hw() file so "self." is necessary
                u[y,x] = self.pointValue_hwvec(x,y,power,smoothing,xv1,yv1,u_vec)
                v[y,x] = self.pointValue_hwvec(x,y,power,smoothing,xv1,yv1,v_vec)
                #u[y,x] = pointValue_hwvec(x,y,power,smoothing,xv1,yv1,u_vec)
                #v[y,x] = pointValue_hwvec(x,y,power,smoothing,xv1,yv1,v_vec)

        ### DKDK comparing two arrays in both shape and values
        ### DKDK from tests, two are not the same but close enough - guessing numerical error/difference between python and numpy
        ### DKDK numpy.allclose() test was passed
        #print numpy.array_equal(u,ut)
        #print numpy.array_equal(v,vt)
        #print numpy.allclose(u,ut)
        #print numpy.allclose(v,vt)

        ds.GetRasterBand(1).WriteArray(u)
        ds.GetRasterBand(2).WriteArray(v)
        ds = None
        ### DKDK no return ?
        return u,v


    ### DKDK new function for a test
    ### Bascially a sort of vectorization of pointValue_hw(), which results in about 15 times faster
    def pointValue_hwvec(self, x,y,power,smoothing,xv1,yv1,values):
        """Used for wind vector field Interpolation

        """

        nominator=0
        denominator=0

        ### DKDK below is to convert x,y to numpy arrays
        x=numpy.repeat(x,1)
        y=numpy.repeat(y,1)

        dist=numpy.sqrt((x-xv1)*(x-xv1)+(y-yv1)*(y-yv1)+smoothing*smoothing);

        nominator=numpy.sum((values/numpy.power(dist,power)))
        denominator=numpy.sum((1/numpy.power(dist,power)))

        #Return NODATA if the denominator is zero
        if denominator > 0:
            value = nominator/denominator
        ### DKDK below is to avoid singularity
        elif (numpy.where(dist<0.0000000001)):
            value = values[numpy.where(dist<0.0000000001)[0][0]]
        else:
            value = -9999
        return value

        """
        if denominator > 0:
            value = nominator/denominator
        else:
            value = -9999
        return value
        """


    ###
    ### DKDK new functions (modified from above functions) to generate two bands GeoTIFF for hurricane wind vector fields
    ###
    def FromJsonToGeoTIFF_hwvec(self,jsonfile,outFile):
        """Used for wind vector fields -
        Read json file, and Create two bands GeoTIFF (for u & v vectors separately) involving rasterization and interpolation

        """

        jsonconf=None
        #with open('results_hwind.json') as f:
        with open(jsonfile) as f:
            jsonconf = json.load(f)
        f.close()

        print jsonfile

        hwdata=jsonconf['output_array_hw']
        thetadata=jsonconf['theta_hw']
        size_one_variable_hw=jsonconf['size_each_variable_and_total_number5']
        size_output_array_hw=jsonconf['size_output_array_hw']
        # size_each_variable_and_total_number5 = [37.0, 14.0, 5.0]

        # DKDK single line like below works too
        hwdata_reshape=numpy.transpose(numpy.reshape(hwdata,(int(size_one_variable_hw[1])*int(size_one_variable_hw[2]),int(size_one_variable_hw[0]))))

        #print hwdata_reshape.shape, len(numpy.transpose(hwdata_reshape))

        # here, 5 = total number of variables, x,y,u,v,total
        num_variables=int(size_one_variable_hw[2])
        lat_hw   =hwdata_reshape[:,numpy.arange(0,len(numpy.transpose(hwdata_reshape)),num_variables)]
        long_hw  =hwdata_reshape[:,numpy.arange(1,len(numpy.transpose(hwdata_reshape)),num_variables)]
        u_hw     =hwdata_reshape[:,numpy.arange(2,len(numpy.transpose(hwdata_reshape)),num_variables)]
        v_hw     =hwdata_reshape[:,numpy.arange(3,len(numpy.transpose(hwdata_reshape)),num_variables)]
        wspeed_hw=hwdata_reshape[:,numpy.arange(4,len(numpy.transpose(hwdata_reshape)),num_variables)]

        print lat_hw.shape, long_hw.shape,u_hw.shape,v_hw.shape,wspeed_hw.shape
        print numpy.amax(lat_hw),numpy.amax(long_hw),numpy.amax(wspeed_hw)
        print numpy.amin(lat_hw),numpy.amin(long_hw),numpy.amin(wspeed_hw)

        Numx=int(size_one_variable_hw[0])
        Numy=int(size_one_variable_hw[1])

        ###
        ### DKDK - for test purpose that can be blocked later - start
        #long_hw_reshp  =numpy.reshape(long_hw,Numx*Numy)
        #lat_hw_reshp   =numpy.reshape(lat_hw,Numx*Numy)
        #wspeed_hw_reshp=numpy.reshape(wspeed_hw,Numx*Numy)
        #print len(long_hw_reshp), long_hw_reshp.shape
        #X1, Y1 = numpy.meshgrid(long_hw_reshp, lat_hw_reshp)
        #V1 = griddata((long_hw_reshp, lat_hw_reshp), wspeed_hw_reshp, (X1, Y1), method='nearest')
        #plt.show()
        ### DKDK - for test purpose that can be blocked later  - end
        ###

        xv=numpy.reshape(long_hw,Numx*Numy);
        yv=numpy.reshape(lat_hw,Numx*Numy);
        u_vec=numpy.reshape(u_hw,Numx*Numy);
        v_vec=numpy.reshape(v_hw,Numx*Numy);
        values=numpy.reshape(wspeed_hw,Numx*Numy);

        print u_vec.shape, v_vec.shape

        #startPos = [geoproperties['tl']['long'],geoproperties['tl']['lat']]
        #d_lat = (geoproperties['br']['lat'] - geoproperties['tl']['lat']) / (Ny - 1)
        #d_long = (geoproperties['br']['long'] - geoproperties['tl']['long']) / (Nx - 1)

        #startPos = [numpy.amin(X1),numpy.amax(Y1)]
        #print startPos, len(Y1)
        #d_lat = (numpy.amin(Y1) - numpy.amax(Y1)) / (len(Y1) - 1)
        #d_long = (numpy.amax(X1) - numpy.amin(X1)) / (len(X1) - 1)


        ### DKDK just in case, added _hwvec (keys) at the geoproperties_hw dictionaries
        geoproperties_hwvec = {
            'tl_hwvec': {
                'lat_hwvec': numpy.amax(yv),
                'long_hwvec': numpy.amin(xv),
            },
            'br_hwvec': {
                'lat_hwvec': numpy.amin(yv),
                'long_hwvec': numpy.amax(xv),
            },
            'Ny_hwvec': len(xv),
            'Nx_hwvec': len(yv),
        }

        #
        # Setting the default values for GeoTIFF
        #
        #geotransform = None
        #proj='osr.SpatialReference()'
        #proj=None
        xSize=100
        ySize=100
        power=5
        smoothing=0
        #zField='Z'
        #dataFile=None
        driverName='GTiff'
        #outFile=hwindgeotiff  # DKDK this is now received via function call

        xMin=min(xv)
        xMax=max(xv)
        yMin=min(yv)
        yMax=max(yv)
        #power=1
        #smoothing=20

        print xMin, xMax, yMin, yMax
        print (xMax-xMin)/xSize
        #print min(xv)

        #geotransform=[]
        #geotransform.append(xMin)
        #geotransform.append((xMax-xMin)/xSize)
        #geotransform.append(0)
        #geotransform.append(yMax)
        #geotransform.append(0)
        #geotransform.append((yMin-yMax)/ySize)

        geotransform=[xMin,(xMax-xMin)/xSize,0,yMax,0,(yMin-yMax)/ySize]
        #geotransform=[10,0.6,0,90,0,0.8]
        print geotransform


        #start_time = time.time()
        #Creating the interpolation function and populating the output matrix value
        #data_hw = invDist_hw(xv,yv,values,geotransform,xSize,ySize,power,smoothing,driverName,outFile)
        ### DKDK
        #u,v = invDist_hwvec(xv,yv,u_vec,v_vec,geotransform,xSize,ySize,power,smoothing,driverName,outFile)
        u,v = self.invDist_hwvec(xv,yv,u_vec,v_vec,geotransform,xSize,ySize,power,smoothing,driverName,outFile)

        ### DKDK although below works, I don't think u, v return is necessary
        #return (geoproperties_hwvec, u, v)
        return geoproperties_hwvec


    ### DKDK new function for vector fields to generate a new mapfile !!!
    ### DKDK please note that shapeFile for vector fields actually indicate GeoTIFF file !!!

    ### DKDK add one more parameter, LOCAL_MAPSERVER_URL - nono
    def exportMapFile_hwvec(self,outFile,shapeFile,projection,extents,mapName):
    #def exportMapFile_hwvec(self,outFile,shapeFile,projection,extents,mapName, def_mapurl):
    #def exportMapFile_hw(self,outFile,shapeFile,projection,extents,mapName,fontLocation="font.list"):
        #self._calculateColors_hw()
        """Exports a MapFile for MapServer - hurricane wind vector fields

        Args:
            | outFile (str): File to write to
            | shapeFile (str): ShapeFile with wind vector fields

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
  #STATUS ON
  # Map image size
  EXTENT """+' '.join([str(e) for e in extents])+"""
  ###
  ### DKDK how to determine size for hurricane wind ? Temporarily set to 512 512
  ###
  #SIZE 512 256
  SIZE 512 512
  MAXSIZE 4096
  IMAGETYPE png24
  #UNITS METERS

  PROJECTION
    #'proj=utm'
    #'zone=5'
    #'ellps=GRS80'
    #'datum=NAD83'
    #'units=m'
    #'no_defs'
    'init="""+projection+"""'
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
    ### DKDK wms_onlineresource may need to be changed for server type (demo, dev, staging)
    ### DKDK not quite sure whether METADATA is required for vector fields
    METADATA
      'wms_title'           '"""+mapName+"""'

      ### DKDK change below to use parameter - nono. It seems that below is not actually used...
      'wms_onlineresource'  'http://cybereye.crc.nd.edu/portal/mapserver/?map="""+abspath+"""'
      'wms_srs'             '"""+projection+""" epsg:4326 epsg:3857'
      'wms_enable_request'  '*'
    END

    #Scale range at which web interface will operate
    # Template and header/footer settings
    # Only the template parameter is required to display a map. See MapServer documentation
    # TEMPLATE 'fooOnlyForWMSGetFeatureInfo'
  END

### DKDK below is for vector fields"
SYMBOL
  NAME "horizline"
  TYPE VECTOR
  POINTS
     0 0
     1 0
  END # points
END # symbol
SYMBOL
  NAME "arrowhead"
  TYPE vector
  FILLED true
  ANCHORPOINT 0 0.5
  POINTS
    0 2
    4 1
    0 0
  END # points
END # symbol
SYMBOL
  NAME "arrowtail"
  TYPE vector
  FILLED true
  ANCHORPOINT 1 0.5 # to shift the arrowtail
  POINTS
    0 2
    4 1
    0 0
    -99 -99
    0 1
    4 1
  END # points
END # symbol

LAYER
  PROJECTION
    'init="""+projection+"""'
  END
    NAME "hwindvector_layer"
    TYPE POINT
    CONNECTIONTYPE UVRASTER
    STATUS DEFAULT
    DATA '"""+shapefileAbsPath+"""'
    PROCESSING "BANDS=1,2"
    ### DKDK UV_SPACING: The spacing is simply the distance, in pixels, between arrows to be displayed in the vector field. Default is 32.
    PROCESSING "UV_SPACING=32"
    PROCESSING "UV_SIZE_SCALE=0.7"
    #PROCESSING "UV_SIZE_SCALE=1.0"
    ### DKDK TRANSPARENCY is being deprecated (0 = transparent) - OPACITY is recommended instead
    ### DKDK 0 = nothing is shown (at least in an image viewer), 1 = looks a bit transparent than other number
    ### DKDK via mapserver, OPACITY effect is more clear than an image viewer
    ###TRANSPARENCY 100
    OPACITY 100

  CLASS
    STYLE
      SYMBOL "horizline"
      ANGLE [uv_angle]
      SIZE [uv_length]
      WIDTH 2
      ### DKDK this is the color for arrow line (not head)
      #COLOR 255 112 112
      COLOR 0 0 0
    END # style
    STYLE
      SYMBOL "arrowhead"
      ANGLE [uv_angle]
      ### DKDK test for size (default = 7)
      SIZE 7
      ### DKDK this is the color for arrowhead
      #COLOR 255 0 0
      COLOR 0 0 0
      POLAROFFSET [uv_length_2] [uv_angle]
    END # style
  END # CLASS ?

END # LAYER ?


  ###
  ### Below are used for hwindspeed - not quite sure if this is needed
  ###
  # Background color for the map canvas -- change as desired
  #IMAGECOLOR 195 220 252
  IMAGECOLOR 255 255 255
  IMAGEQUALITY 95
  IMAGETYPE png24

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
    ### DKDK RGBA indicates transparent background - http://mapserver.org/output/agg.html
    IMAGEMODE RGBA
    ### DKDK not quite sure whether the option TRANSPARENT works
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

END # Mapfile
"""
        with open(outFile,"w") as f:
            f.write(mapFileContents)