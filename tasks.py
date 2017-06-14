# tasks.py
#
# Nathan Smith, Dae Kun Kwon, Samuel Njoroge
#
# File defines the celery tasks and helper functions/classes to run both
# the wave and wind model. Both tasks were put there to avoid creating
# another app. Using seperate files does not work because of how celery
# handles the jobs. S
# Same result model (HAKOUv2Result) is used and columns were added to
# help accomodate for new parameters in the wind model.
# These new variables are:
#
# dd - Distance from landfall - defaults to 0 and not available in interface.
# dcp - Pressure difference - Minimum of 1 from the interface.
# z0 - Terrain properties.
# related - id if a related model run (i.e. if a wave/wind model were run at the
# same time with the same parameters). In this way, if they are loaded for
# visualization, both models are loaded. Also, this is shown on lists.
# wind - quick way to determine if model result (row) is a wave or wind.

import os, sys, time
import subprocess

from django.conf import settings
from celery.task import Task
from celery.registry import tasks
from celery.contrib import rdb

import MLSOperations
import MLSOperations.utils
import WaveContourGenerator
import WaveContourCalculator_hw
import WaveContourGenerator.utils

# from WaveContourCalculator.models import WaveContourCalculation
from HAKOUpy.v2.utils import (getGeoPropertiesFromJSON,
                              convertBinToGeoTIFF,
                              convertProbabilityLevelToGeoTIFF,
                              convertBinToShp,
                              exportSurgeMapFile,
                              convertProbabilityLevelToShp,
                              convertTrackPlotJSONtoShp,
                              exportTrackPlotMapFile)
import HAKOUpy.v2.utils as HAKOUv2utils
import HAKOUpy.v3.utils as HAKOUv3utils

import json
import urllib
from osgeo import ogr

try:
    import ipdb as pdb
except:
    import pdb

class CalledProcessError(Exception):
    """Exception with the following properties:

    * returncode: Returned status code of the called process
    * cmd: The command that was run
    * output (optional or None): Output of the command
    """
    def __init__(self, returncode, cmd, output=None):
        self.returncode = returncode
        self.cmd = cmd
        self.output = output
        super(CalledProcessError, self).__init__(returncode,cmd,output)

def saveOutputDict(saveFile,dictObj):
    try:
        f = file(saveFile, 'w')
        json.dump(dictObj,f)
        f.close()

        return True
    except IOError:
        return False

class contourModelCalc(Task):
    """
    """

    README_Contents = """Cybereye Wave Contour Calculator Results (cybereye.crc.nd.edu)

The Wave Contour Calculator is an implementation of the HAKOU model from the
Environmental Fluid Dynamics group of the Department of Civil & Environmental
Engineering & Earth Sciences at the University of Notre Dame
(http://ceees.nd.edu/research-facilities/research-focus-areas/environmental-fluid-dynamics-1).
It employs Moving Least Squares analysis to perform hurricane surge risk
assessment in real-time based on high-fidelity simulation data.

In the version of the model that produced these results two areas are output
with average wave height information for each.  There are the same results for
each area with a different identifying name.  All files are in the EPSG:4326
projection.  Contours plots go from 0-20m in 1m steps.  They are as follows:

mlsresults_close.csv 	  - Average wave height in CSV format
mlsresults_close.tif 	  - Average wave height in GEOTIFF format
output_close.shp	 	  - Shapefile of wave height contours
output_close_reproj.shp	  - Shapefile with projection information
output_close.map		  - MapServer mapfile for WMS access to the shapefiles
						    along with styling
output_close_overview.png - Overview image of the results contour plot

settingsDict.json		  - JSON document with settings used to produce these results
outputDict.json			  - JSON document with locations and names of all the
							output file along with urls to access the MapServer
							results.  This is primarily for system use."""

    def updateState(self, status='', extra='', increment=True, failed=False):
        """Updates the state of the celery task with a message and increments
        the count if increment is True.
        """
        state = 'PROGRESS'
        if failed:
            state = 'FAILURE'
            increment = False

        if increment:
            self.current += 1

        self.update_state(state=state, meta={
            'progress': float(self.current)/self.total,
            'status': status,
            'extra': extra,
        })

    def run(self,initial_parameters,fileseed,wcc_pk,outputDir=None):
        try:
            self.total = 11
            self.current = 0

            self.updateState('Establishing settings')

            inputDir = os.path.join(settings.PROJECT_ROOT, 'data/')
            if outputDir == None:
                outputDir = os.path.join(settings.MEDIA_ROOT, 'outputs/')

            mlscsv1 = 'mlsresults_close.csv'
            mlscsv2 = 'mlsresults_wide.csv'
            mlsgeotiff1 = 'mlsresults_close.tif'
            mlsgeotiff2 = 'mlsresults_wide.tif'
            mlsgeotiff1_reprojected = 'mlsresults_close_reproj.tif'
            mlsgeotiff2_reprojected = 'mlsresults_wide_reproj.tif'

            outFile1 = 'output_close.shp'
            reprojectedFile1 = 'output_close_reproj.shp'
            mapFile1 = 'output_close.map'
            overview1 = 'output_close_overview.png'
            outFile2 = 'output_wide.shp'
            reprojectedFile2 = 'output_wide_reproj.shp'
            mapFile2 = 'output_wide.map'
            overview2 = 'output_wide_overview.png'
            fontLocation = os.path.join(inputDir,'fonts/font.list')

            dataDir = inputDir
            result_geoproperties = {
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

            initial_parameters.append(str(dataDir))
            initial_parameters.append(str(outputDir))
            initial_parameters.append(str(mlscsv1))
            initial_parameters.append(str(mlscsv2))

            print initial_parameters

            self.updateState('Running MLS')

            # Run the moving least squares analysis
            MLSOperations.runMLS(*initial_parameters)

            self.updateState('Converting to GeoTiff')

            # convert the resulting csv files to geotiff
            MLSOperations.convertCSVToGeoTIFF(result_geoproperties['Hsin'], outputDir+mlscsv1, outputDir, mlsgeotiff1)
            MLSOperations.convertCSVToGeoTIFF(result_geoproperties['Hsout'], outputDir+mlscsv2, outputDir, mlsgeotiff2)


            # first
            self.updateState('Calculating first contours')

            contour = WaveContourGenerator.WaveContourGenerator()
            contour.readGeoTIFF(outputDir+mlsgeotiff1)
            #contour.read(inputDir, 'wave_data.csv','wave_data_lat.csv','wave_data_long.csv')
            #contour.createPlot(levels)
            contour.createPlot()

            self.updateState('Exporting & reprojecting first shapefile')

            contour.exportShp(os.path.join(outputDir,outFile1))
            #WaveContourGenerator.utils.reprojectShp(os.path.join(outputDir, outFile1), os.path.join(outputDir, reprojectedFile1))
            MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',os.path.join(outputDir, outFile1), os.path.join(outputDir, reprojectedFile1),True)

            self.updateState('Exporting first mapfile')

            extents1 = MLSOperations.utils.getShpExtents(os.path.join(outputDir, reprojectedFile1))
            contour.exportMapFile(os.path.join(outputDir, mapFile1), os.path.join(outputDir, reprojectedFile1),'EPSG:4326', extents1, 'OahuWaveHeight',fontLocation)

            # second
            self.updateState('Calculating second contours')

            contour = WaveContourGenerator.WaveContourGenerator()
            contour.readGeoTIFF(outputDir+mlsgeotiff2)
            #contour.read(inputDir, 'wave_data.csv','wave_data_lat.csv','wave_data_long.csv')
            #contour.createPlot(levels)
            contour.createPlot()

            self.updateState('Exporting & reprojecting second shapefile')

            contour.exportShp(os.path.join(outputDir,outFile2))
            #WaveContourGenerator.utils.reprojectShp(os.path.join(outputDir, outFile2), os.path.join(outputDir, reprojectedFile2))
            MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',os.path.join(outputDir, outFile2), os.path.join(outputDir, reprojectedFile2),True)

            self.updateState('Exporting second mapfile')

            extents2 = MLSOperations.utils.getShpExtents(os.path.join(outputDir, reprojectedFile2))
            contour.exportMapFile(os.path.join(outputDir, mapFile2), os.path.join(outputDir, reprojectedFile2), 'EPSG:4326', extents2, 'OahuWaveHeight',fontLocation)

            self.updateState('Saving settings')

            #### Saving settings ####

            mapFiles = [os.path.join(outputDir, mapFile1), os.path.join(outputDir, mapFile2), extents1, extents2]

            mapURL1 = settings.MAPSERVER_URL + "?map=" + mapFiles[0]
            mapURL2 = settings.MAPSERVER_URL + "?map=" + mapFiles[1]
            extents1 = mapFiles[2]
            extents2 = mapFiles[3]
            imageURL1_QS = "&LAYERS=OahuWaveHeight&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents1])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents1])+"&WIDTH=700&HEIGHT=508"
            imageURL2_QS = "&LAYERS=OahuWaveHeight&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents2])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents2])+"&WIDTH=700&HEIGHT=508"
            imageURL1 = mapURL1 + imageURL1_QS
            imageURL2 = mapURL2 + imageURL2_QS


            wcc_settings = {
                'inputDir': inputDir,
                'outputDir': outputDir,
                'fileseed': fileseed,
                'mlscsv1': mlscsv1,
                'mlscsv2': mlscsv2,
                'mlsgeotiff1': mlsgeotiff1,
                'mlsgeotiff2': mlsgeotiff2,
                'mlsgeotiff1_reprojected': mlsgeotiff1_reprojected,
                'mlsgeotiff2_reprojected': mlsgeotiff2_reprojected,

                'outFile1': outFile1,
                'reprojectedFile1': reprojectedFile1,
                'mapFile1': mapFile1,
                'outFile2': outFile2,
                'reprojectedFile2': reprojectedFile2,
                'mapFile2': mapFile2,
                'fontLocation': fontLocation,
                'dataDir': dataDir,
                'result_geoproperties': result_geoproperties,
            }

            d = {
                'contour1': mapFiles[0],
                'mapURL1': mapURL1,
                'imageURL1': imageURL1,
                'contour2': mapFiles[1],
                'mapURL2': mapURL2,
                'imageURL2': imageURL2,
                'num1': initial_parameters[0],
                'num2': initial_parameters[1],
                'num3': initial_parameters[2],
                'num4': initial_parameters[3],
                'num5': initial_parameters[4],

                'MAPSERVER_BASELAYERS_URL': settings.MAPSERVER_URL,
                'TILECACHE_BASELAYERS_URL': settings.TILECACHE_BASELAYERS_URL,
            }

            saveOutputDict(os.path.join(outputDir,'outputDict.json'),d)
            saveOutputDict(os.path.join(outputDir,'settingsDict.json'),wcc_settings)

            # save README.txt
            with open(os.path.join(outputDir,'README.txt'),'w') as f:
                f.write(self.README_Contents)

            # fetching overview images
            self.updateState('Fetching overview images')
            urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+mapFiles[0]+imageURL1_QS, os.path.join(outputDir, overview1))
            urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+mapFiles[1]+imageURL2_QS, os.path.join(outputDir, overview2))

            #pdb.set_trace()

        except:
            #self.updateState(self, extra=traceback.format_exc(), failed=True)
            raise


class HAKOUModel_v2(contourModelCalc):
    """
    """
    waveprobbin = 'result_wave_prob.bin'
    wavejson = 'result_wave.json'
    wavegeotiff = 'result_wave_prob_%s.tiff'
    waveshape = 'result_wave_prob_%s.shp'
    wavemap = 'result_wave_prob_%s.map'
    waveoverview = 'result_wave_prob_overview_%s.png'
    wavelayername = 'wave'
    wavegeoprop = None

    surgeprobbin = 'result_surge_prob.bin'
    surgejson = 'result_surge.json'
    surgeshape = 'result_surge_prob_%s.shp'
    surgemap = 'result_surge_prob_%s.map'
    surgeoverview = 'result_surge_prob_overview_%s.png'
    surgelayername = 'surge'

    def run(self, dataDir, outputDir, vector, latitude, wave=True, surge=False, wavedetonly=True, surgedetonly=True, track=True):
        self.dataDir = dataDir
        self.resultDir = outputDir
        self.outputDir = outputDir
        self.force = False

        try:
            self.total = 2
            if wave == True:
                self.total += 7
                if wavedetonly == True:
                    self.total += 1
            if surge == True:
                self.total += 4
                if surgedetonly == True:
                    self.total += 1
            if track == True:
                self.total += 3
            self.current = 0

            self.updateState('Establishing settings')

            #rdb.set_trace()

            wcc_settings = {
                'dataDir': dataDir,
                'outputDir': outputDir,
                'fontLocation': settings.MAPSERVER_FONTS_LOCATION,
                'wave': wave,
                'surge': surge,
                'wavedetonly': wavedetonly,
                'surgedetonly': surgedetonly,
            }

            d = {
                'longitude': vector[0],
                'angle': vector[1],
                'cp': vector[2],
                'vf': vector[3],
                'Rmax': vector[4],
                'time': vector[5],

                'MAPSERVER_BASELAYERS_URL': settings.MAPSERVER_URL,
                'TILECACHE_BASELAYERS_URL': settings.TILECACHE_BASELAYERS_URL,
            }

            # Wave portion
            if wave == True:
                self.updateState('Running wave model')
                cmd = [
                    settings.HAKOU_V2_EXEC,
                    '--iswave',
                    '--detonly' if wavedetonly else '',
                    '--',
                    os.path.join(dataDir,'wave'),
                    outputDir,
                    ' '.join([str(i) for i in vector]),
                ]

                process = subprocess.Popen(' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()
                dumpstr = process.stdout.readlines()
                errorstr = process.stderr.readlines()

                print ' '.join(cmd)
                print ''.join(dumpstr)
                print ''.join(errorstr)

                if process.returncode != 0:
                    raise CalledProcessError(process.returncode, ' '.join(cmd), 'OUT:\n'+''.join(dumpstr)+'\n\nERROR:\n'+''.join(errorstr))

            # Surge portion
            if surge == True:
                self.updateState('Running surge model')
                cmd = [
                    settings.HAKOU_V2_EXEC,
                    '--detonly' if surgedetonly else '',
                    '--',
                    os.path.join(dataDir,'surge'),
                    outputDir,
                    ' '.join([str(i) for i in vector]),
                ]

                process = subprocess.Popen(' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()
                dumpstr = process.stdout.readlines()
                errorstr = process.stderr.readlines()

                print ' '.join(cmd)
                print ''.join(dumpstr)
                print ''.join(errorstr)

                if process.returncode != 0:
                    raise CalledProcessError(process.returncode, ' '.join(cmd), 'OUT:\n'+''.join(dumpstr)+'\n\nERROR:\n'+''.join(errorstr))

            # Track portion
            if track == True:
                self.updateState('Running track plotter')
                cmd = [
                    settings.HAKOU_TRACK_EXEC,
                    '--',
                    os.path.join(dataDir,'track'),
                    outputDir,
                    str(latitude),
                    str(vector[0]),
                    str(vector[1]),
                    str(vector[3]),
                    str(vector[5]),
                ]

                process = subprocess.Popen(' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                process.wait()
                dumpstr = process.stdout.readlines()
                errorstr = process.stderr.readlines()

                print ' '.join(cmd)
                print ''.join(dumpstr)
                print ''.join(errorstr)

                if process.returncode != 0:
                    raise CalledProcessError(process.returncode, ' '.join(cmd), 'OUT:\n'+''.join(dumpstr)+'\n\nERROR:\n'+''.join(errorstr))

            # Wave output
            if wave == True:
                wavedetbin = 'result_wave_deter.bin'
                wavejson = 'result_wave.json'
                wavegeotiff = 'result_wave_deter.tiff'
                waveshape = 'result_wave_deter.shp'
                wavemap = 'result_wave_deter.map'
                waveoverview = 'result_wave_deter_overview.png'
                wavelayername = 'wave'

                self.updateState('Getting geographic extents of wave model outputs')
                geoprop = getGeoPropertiesFromJSON(os.path.join(outputDir,wavejson))

                self.updateState('Converting to Geotiff')
                convertBinToGeoTIFF(geoprop,os.path.join(outputDir,wavedetbin),os.path.join(outputDir,wavegeotiff))

                self.updateState('Generating wave contours')
                contour = WaveContourGenerator.WaveContourGenerator()
                contour.readGeoTIFF(os.path.join(outputDir,wavegeotiff))
                contour.createPlot()

                self.updateState('Exporting & reprojecting wave shapefile')
                contour.exportShp(os.path.join(outputDir,waveshape))
                MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',os.path.join(outputDir, waveshape), os.path.join(outputDir, waveshape),True)

                self.updateState('Exporting wave mapfile')
                extents = MLSOperations.utils.getShpExtents(os.path.join(outputDir, waveshape))
                contour.exportMapFile(os.path.join(outputDir, wavemap), os.path.join(outputDir, waveshape),'EPSG:4326', extents, wavelayername,settings.MAPSERVER_FONTS_LOCATION)

                mapURL = settings.MAPSERVER_URL + "?map=" + os.path.join(outputDir,wavemap)
                imageURL_QS = "&LAYERS="+wavelayername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents])+"&WIDTH=700&HEIGHT=508"
                imageURL = mapURL + imageURL_QS

                # fetching wave overview image
                self.updateState('Fetching wave overview image')
                urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+os.path.join(outputDir,wavemap)+imageURL_QS, os.path.join(outputDir, waveoverview))

                if wavedetonly == False:
                    self.wavegeoprop = getGeoPropertiesFromJSON(os.path.join(self.resultDir,self.wavejson))
                    self.processProbWave(10.0)
                    self.processProbWave(50.0)
                    self.processProbWave(90.0)

                wcc_settings.update({
                    'wavebin': wavedetbin,
                    'wavegeotiff': wavegeotiff,
                    'waveshape': waveshape,
                    'wavemap': wavemap,
                    'wavejson': wavejson,
                    'waveoverview': waveoverview,
                    'wave_geoproperties': geoprop,
                })

                d.update({
                    'wave_map': os.path.join(outputDir, wavemap),
                    'wave_mapURL': mapURL,
                    'wave_imageURL': imageURL,
                    'wave_layername': wavelayername,
                    'wave_extents': extents,
                })

            # Surge output
            if surge == True:
                surgedetbin = 'result_surge_deter.bin'
                surgejson = 'result_surge.json'
                surgeshape = 'result_surge_deter.shp'
                surgemap = 'result_surge_deter.map'
                surgeoverview = 'result_surge_deter_overview.png'
                surgelayername = 'surge'

                self.updateState('Exporting & reprojecting surge shapefile')
                convertBinToShp(os.path.join(dataDir,'surge'),os.path.join(outputDir,surgedetbin),os.path.join(outputDir,surgeshape))
                MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',os.path.join(outputDir, surgeshape), os.path.join(outputDir, surgeshape),True)

                self.updateState('Exporting surge mapfile')
                extents = MLSOperations.utils.getShpExtents(os.path.join(outputDir, surgeshape))
                exportSurgeMapFile(os.path.join(outputDir, surgemap), os.path.join(outputDir, surgeshape),'EPSG:4326', extents, surgelayername,settings.MAPSERVER_FONTS_LOCATION)

                mapURL = settings.MAPSERVER_URL + "?map=" + os.path.join(outputDir,surgemap)
                imageURL_QS = "&LAYERS="+surgelayername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents])+"&WIDTH=700&HEIGHT=508"
                imageURL = mapURL + imageURL_QS

                # fetching surge overview image
                self.updateState('Fetching surge overview image')
                urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+os.path.join(outputDir,surgemap)+imageURL_QS, os.path.join(outputDir, surgeoverview))

                if surgedetonly == False:
                    self.processProbSurge(10.0)
                    self.processProbSurge(50.0)
                    self.processProbSurge(90.0)

                wcc_settings.update({
                    'surgebin': surgedetbin,
                    'surgeshape': surgeshape,
                    'surgemap': surgemap,
                    'surgejson': surgejson,
                    'surgeoverview': surgeoverview,
                })

                d.update({
                    'surge_map': os.path.join(outputDir, surgemap),
                    'surge_mapURL': mapURL,
                    'surge_imageURL': imageURL,
                    'surge_layername': surgelayername,
                    'surge_extents': extents,
                })

            # track output
            if track == True:
                self.updateState('Producing track output')

                trackjson = 'track_plot.json'
                trackmap = 'result_trackplot.map'
                tracknames = ['result_trackplot_output', 'result_trackplot_cone', 'result_trackplot_hist']
                trackshapes = ['%s.shp'%x for x in tracknames]
                trackoverview = 'result_trackplot_overview.png'
                tracklayername = 'trackplot'

                convertTrackPlotJSONtoShp(os.path.join(outputDir,trackjson), outputDir, tracknames)
                extents = []

                shapefiles = [os.path.join(outputDir,x) for x in trackshapes]
                for shp in shapefiles:
                    MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',shp,shp,True)
                    extents.append( MLSOperations.utils.getShpExtents(shp) )

                combinedExtents = extents[0]
                for ext in extents:
                    combinedExtents[0] = min(combinedExtents[0],ext[0])
                    combinedExtents[1] = min(combinedExtents[1],ext[1])
                    combinedExtents[2] = max(combinedExtents[2],ext[2])
                    combinedExtents[3] = max(combinedExtents[3],ext[3])

                exportTrackPlotMapFile(os.path.join(outputDir,trackmap), shapefiles, 'EPSG:4326', extents, combinedExtents, tracklayername, settings.MAPSERVER_FONTS_LOCATION)

                mapURL = settings.MAPSERVER_URL + "?map=" + os.path.join(outputDir,trackmap)
                imageURL_QS = "&LAYERS="+tracklayername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in combinedExtents])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in combinedExtents])+"&WIDTH=700&HEIGHT=508"
                imageURL = mapURL + imageURL_QS

                # fetching track overview image
                self.updateState('Fetching track overview image')
                urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+os.path.join(outputDir,trackmap)+imageURL_QS, os.path.join(outputDir, trackoverview))

                wcc_settings.update({
                    'trackshapes': trackshapes,
                    'trackmap': trackmap,
                    'trackjson': trackjson,
                    'trackoverview': trackoverview,
                })

                d.update({
                    'track_map': os.path.join(outputDir, trackmap),
                    'track_mapURL': mapURL,
                    'track_imageURL': imageURL,
                    'track_layername': tracklayername,
                    'track_extents': combinedExtents,
                })

            self.updateState('Saving settings and outputs')
            saveOutputDict(os.path.join(outputDir,'outputDict.json'),d)
            saveOutputDict(os.path.join(outputDir,'settingsDict.json'),wcc_settings)

            # save README.txt
            with open(os.path.join(outputDir,'README.txt'),'w') as f:
                f.write(self.README_Contents)

        except:
            raise

    def processProbWave(self, level):
        """
        :param level:
        :return:
        """
        probdir = os.path.join(self.outputDir,'wave_%s'%level)
        if not os.path.exists(probdir) or self.force == True:
            if not os.path.exists(probdir):
                os.mkdir(probdir)

            self.updateState('Procesing probability level %s wave results'%level)

            # convert to Geotiff
            convertProbabilityLevelToGeoTIFF(self.wavegeoprop, os.path.join(self.resultDir,self.waveprobbin), os.path.join(probdir,self.wavegeotiff%level), level)

            # generate contours
            contour = WaveContourGenerator.WaveContourGenerator()
            contour.readGeoTIFF(os.path.join(probdir,self.wavegeotiff%level))
            contour.createPlot()

            # export to shapefile
            contour.exportShp(os.path.join(probdir,self.waveshape%level))
            MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',os.path.join(probdir, self.waveshape%level), os.path.join(probdir, self.waveshape%level),True)

            # export mapfile
            extents = MLSOperations.utils.getShpExtents(os.path.join(probdir, self.waveshape%level))
            contour.exportMapFile(os.path.join(probdir, self.wavemap%level), os.path.join(probdir, self.waveshape%level),'EPSG:4326', extents, self.wavelayername,settings.MAPSERVER_FONTS_LOCATION)

            mapURL = settings.MAPSERVER_URL + "?map=" + os.path.join(probdir,self.wavemap%level)
            imageURL_QS = "&LAYERS="+self.wavelayername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents])+"&WIDTH=700&HEIGHT=508"
            imageURL = mapURL + imageURL_QS

            # fetching wave overview image
            urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+os.path.join(probdir,self.wavemap%level)+imageURL_QS, os.path.join(probdir, self.waveoverview%level))
        else:
            self.updateState('Probability level %s wave results already exist'%level)

    def processProbSurge(self, level):
        """
        :param level:
        :return:
        """
        probdir = os.path.join(self.outputDir,'surge_%s'%level)
        if not os.path.exists(probdir) or self.force == True:
            if not os.path.exists(probdir):
                os.mkdir(probdir)

            self.updateState('Procesing probability level %s surge results'%level)

            # convert to shapefile
            convertProbabilityLevelToShp(os.path.join(self.dataDir,'surge'),os.path.join(self.resultDir,self.surgeprobbin),os.path.join(probdir,self.surgeshape%level),level)
            MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',os.path.join(probdir, self.surgeshape%level), os.path.join(probdir, self.surgeshape%level),True)

            # export mapfile
            extents = MLSOperations.utils.getShpExtents(os.path.join(probdir, self.surgeshape%level))
            exportSurgeMapFile(os.path.join(probdir, self.surgemap%level), os.path.join(probdir, self.surgeshape%level),'EPSG:4326', extents, self.surgelayername,settings.MAPSERVER_FONTS_LOCATION)

            mapURL = settings.MAPSERVER_URL + "?map=" + os.path.join(probdir,self.surgemap%level)
            imageURL_QS = "&LAYERS="+self.surgelayername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents])+"&WIDTH=700&HEIGHT=508"
            imageURL = mapURL + imageURL_QS

            # fetching surge overview image
            urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+os.path.join(probdir,self.surgemap%level)+imageURL_QS, os.path.join(probdir, self.surgeoverview%level))
        else:
            self.updateState('Probability level %s surge results already exist'%level)


class processProbabilityLevels_v2(HAKOUModel_v2):
    """
    """
    def run(self, dataDir, resultDir, outputDir, levels, wave=True, surge=False, force=False):
        self.dataDir = dataDir
        self.resultDir = resultDir
        self.outputDir = outputDir
        self.force = force

        try:
            # calculate steps
            self.total = len(levels)*2 if wave and surge else len(levels)
            self.current = 0

            if wave == True:
                self.wavegeoprop = getGeoPropertiesFromJSON(os.path.join(resultDir,self.wavejson))

            # process each level
            for level in levels:
                if wave == True:
                    self.processProbWave(level)
                if surge == True:
                    self.processProbSurge(level)
        except:
            raise



class HAKOUModel_v3(contourModelCalc):
    """
    """
    waveprobbin = 'result_prob.bin'
    wavejson = 'results.json'
    wavegeotiff = 'result_wave_prob_%s.tiff'
    waveshape = 'result_wave_prob_%s.shp'
    wavemap = 'result_wave_prob_%s.map'
    waveoverview = 'result_wave_prob_overview_%s.png'
    wavelayername = 'wave'
    wavegeoprop = None

    surgeprobbin = 'result_prob.bin'
    surgejson = 'results.json'
    surgeshape = 'result_surge_prob_%s.shp'
    surgemap = 'result_surge_prob_%s.map'
    surgeoverview = 'result_surge_prob_overview_%s.png'
    surgelayername = 'surge'

    README_Contents = """Cybereye Wave Contour Calculator Results (cybereye.crc.nd.edu)

************* Deterministic Outputs *************

outputDict.json			  - JSON document with locations and names of all the
							output file along with urls to access the MapServer
							results.  This is primarily for system use.
settingsDict.json		  - JSON document with settings used to produce these
							results.
results.json			  - JSON document produced by the HAKOU model with
							information about the wave grid, input parameters, etc.
track_plot.json			  - JSON document produced by the TrackPlot program
							about the historical path, deterministic path, and
							probability cone.  Used to produce the trackplot
							shapefiles.
result_deter.bin		  - Binary file with outputs from the HAKOU model.
 							Contains both wave and surge deterministic results.

result_surge_deter_overview.png - Overview image of the deterministic surge results
result_surge_deter.map	  - MapServer Mapfile for deterministic surge results WMS.
result_surge_deter.shp	  - Shapefile (.shp, .dbf, .shx, & .prj) of the
							deterministic surge results.

result_trackplot_overview.png - Overview image of the storm track.
result_trackplot_cone.shp - Shapefile (.shp, .dbf, .shx, & .prj) of the storm
							probability cone.
result_trackplot_hist	  - Shapefile (.shp, .dbf, .shx, & .prj) of the storm
							historic path.
result_trackplot_output   - Shapefile (.shp, .dbf, .shx, & .prj) of the storm
							deterministic path.

result_wave_deter_overview.png - Overview image of the deterministic wave results
result_wave_deter.map	  - MapServer Mapfile for deterministic wave results WMS.
result_wave_deter.shp	  - Shapefile (.shp, .dbf, .shx, & .prj) of the
							deterministic wave results contour plot.  The wave
							contour plot is divided into levels from 0m-20m in
							1m intervals.
result_wave_deter.tiff	  - GeoTIFF of the deterministic wave results.
							Graphical representation of the gridded wave
							results in result_deter.bin and used to create the
							contour plot.

************* Probabilistic Outputs *************

result_prob.bin			  - Binary file with outputs from the HAKOU model.
							Contains both wave and surge deterministic results
							for all probability levels (1%-100% in steps of
							0.1% for 1000 levels for each).

For each processed probability level there will be a folder labeled
<surge|wave>_<level>.  Inside each directory is an overview image (.png) as
explained above, a Mapfile (.map) for the MapServer WMS, and a shapefile (.shp,
.dbf, .shx, & .prj) of the results.  If the directory is for a wave result
there will also be a GeoTIFF with the gridded results."""

    def run(self, dataDir, outputDir, vector, latitude, detonly=True, track=True):
        self.dataDir = dataDir
        self.resultDir = outputDir
        self.outputDir = outputDir
        self.force = False

        try:
            self.total = 12
            if detonly is False:
                self.total += 4
            if track is True:
                self.total += 3
            self.current = 0

            self.updateState('Establishing settings')

            #rdb.set_trace()

            wcc_settings = {
                'dataDir': dataDir,
                'outputDir': outputDir,
                'fontLocation': settings.MAPSERVER_FONTS_LOCATION,
                'detonly': detonly,
            }

            d = {
                'longitude': vector[0],
                'angle': vector[1],
                'cp': vector[2],
                'vf': vector[3],
                'Rmax': vector[4],
                'time': vector[5],

                'MAPSERVER_BASELAYERS_URL': settings.MAPSERVER_URL,
            }

            # Create a vector for the flood.
            flood_vector = [vector[0],vector[1],vector[2],vector[3],vector[4],vector[5]]

            # Run Model
            self.updateState('Running model')
            cmd = [
                settings.HAKOU_V3_EXEC,
                '--detonly' if detonly else '',
                '--',
                os.path.join(dataDir, 'data'),
                outputDir,
                ' '.join([str(i) for i in flood_vector]),
            ]

            process = subprocess.Popen(' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            process.wait()
            dumpstr = process.stdout.readlines()
            errorstr = process.stderr.readlines()

            print ' '.join(cmd)
            print ''.join(dumpstr)
            print ''.join(errorstr)

            if process.returncode != 0:
                raise CalledProcessError(process.returncode, ' '.join(cmd), 'OUT:\n'+''.join(dumpstr)+'\n\nERROR:\n'+''.join(errorstr))

            # Track portion
            if track == True:
                self.updateState('Running track plotter')
                runTrackPlot(dataDir, outputDir, flood_vector, latitude)

            #### Wave output ####
            wavedetbin = 'result_deter.bin'
            wavejson = 'results.json'
            wavegeotiff = 'result_wave_deter.tiff'
            waveshape = 'result_wave_deter.shp'
            wavemap = 'result_wave_deter.map'
            waveoverview = 'result_wave_deter_overview.png'
            wavelayername = 'wave'

            self.updateState('Getting geographic extents of wave model outputs')
            geoprop = HAKOUv3utils.getGeoPropertiesFromJSON(os.path.join(outputDir,wavejson))

            self.updateState('Converting to Geotiff')
            HAKOUv3utils.convertBinToGeoTIFF(geoprop,os.path.join(dataDir,'data'),os.path.join(outputDir,wavedetbin),os.path.join(outputDir,wavegeotiff))

            self.updateState('Generating wave contours')
            contour = WaveContourGenerator.WaveContourGenerator()
            contour.readGeoTIFF(os.path.join(outputDir,wavegeotiff))
            contour.createPlot()

            self.updateState('Exporting & reprojecting wave shapefile')
            contour.exportShp(os.path.join(outputDir,waveshape))
            MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',os.path.join(outputDir, waveshape), os.path.join(outputDir, waveshape),True)

            self.updateState('Exporting wave mapfile')
            extents = MLSOperations.utils.getShpExtents(os.path.join(outputDir, waveshape))
            contour.exportMapFile(os.path.join(outputDir, wavemap), os.path.join(outputDir, waveshape),'EPSG:4326', extents, wavelayername,settings.MAPSERVER_FONTS_LOCATION)

            mapURL = settings.MAPSERVER_URL + "?map=" + os.path.join(outputDir,wavemap)
            imageURL_QS = "&LAYERS="+wavelayername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents])+"&WIDTH=700&HEIGHT=508"
            imageURL = mapURL + imageURL_QS

            # fetching wave overview image
            self.updateState('Fetching wave overview image')
            urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+os.path.join(outputDir,wavemap)+imageURL_QS, os.path.join(outputDir, waveoverview))

            if detonly == False:
                self.wavegeoprop = geoprop
                self.processProbWave(10.0)
                self.processProbWave(50.0)
                # self.processProbWave(90.0)

            wcc_settings.update({
                'wavebin': wavedetbin,
                'wavegeotiff': wavegeotiff,
                'waveshape': waveshape,
                'wavemap': wavemap,
                'wavejson': wavejson,
                'waveoverview': waveoverview,
                'wave_geoproperties': geoprop,
            })

            d.update({
                'wave_map': os.path.join(outputDir, wavemap),
                'wave_mapURL': mapURL,
                'wave_imageURL': imageURL,
                'wave_layername': wavelayername,
                'wave_extents': extents,
            })

            #### Surge output ####
            surgedetbin = 'result_deter.bin'
            surgejson = 'results.json'
            surgeshape = 'result_surge_deter.shp'
            surgemap = 'result_surge_deter.map'
            surgeoverview = 'result_surge_deter_overview.png'
            surgelayername = 'surge'

            self.updateState('Exporting & reprojecting surge shapefile')
            HAKOUv3utils.convertBinToShp(os.path.join(dataDir,'data'),os.path.join(outputDir,surgedetbin),os.path.join(outputDir,surgeshape))
            MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',os.path.join(outputDir, surgeshape), os.path.join(outputDir, surgeshape),True)

            self.updateState('Exporting surge mapfile')
            extents = MLSOperations.utils.getShpExtents(os.path.join(outputDir, surgeshape))
            HAKOUv3utils.exportSurgeMapFile(os.path.join(outputDir, surgemap), os.path.join(outputDir, surgeshape),'EPSG:4326', extents, surgelayername,settings.MAPSERVER_FONTS_LOCATION)

            mapURL = settings.MAPSERVER_URL + "?map=" + os.path.join(outputDir,surgemap)
            imageURL_QS = "&LAYERS="+surgelayername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents])+"&WIDTH=700&HEIGHT=508"
            imageURL = mapURL + imageURL_QS

            # fetching surge overview image
            self.updateState('Fetching surge overview image')
            urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+os.path.join(outputDir,surgemap)+imageURL_QS, os.path.join(outputDir, surgeoverview))

            if detonly == False:
                self.processProbSurge(10.0)
                self.processProbSurge(50.0)
                # self.processProbSurge(90.0)

            wcc_settings.update({
                'surgebin': surgedetbin,
                'surgeshape': surgeshape,
                'surgemap': surgemap,
                'surgejson': surgejson,
                'surgeoverview': surgeoverview,
            })

            d.update({
                'surge_map': os.path.join(outputDir, surgemap),
                'surge_mapURL': mapURL,
                'surge_imageURL': imageURL,
                'surge_layername': surgelayername,
                'surge_extents': extents,
            })

            #### track output ####
            if track == True:
                self.updateState('Producing track output')

                trackjson = 'track_plot.json'
                trackmap = 'result_trackplot.map'
                tracknames = ['result_trackplot_output', 'result_trackplot_cone', 'result_trackplot_hist']
                trackshapes = ['%s.shp'%x for x in tracknames]
                trackoverview = 'result_trackplot_overview.png'
                tracklayername = 'trackplot'

                HAKOUv2utils.convertTrackPlotJSONtoShp(os.path.join(outputDir,trackjson), outputDir, tracknames)
                extents = []

                shapefiles = [os.path.join(outputDir,x) for x in trackshapes]
                for shp in shapefiles:
                    MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',shp,shp,True)
                    extents.append( MLSOperations.utils.getShpExtents(shp) )

                combinedExtents = extents[0]
                for ext in extents:
                    combinedExtents[0] = min(combinedExtents[0],ext[0])
                    combinedExtents[1] = min(combinedExtents[1],ext[1])
                    combinedExtents[2] = max(combinedExtents[2],ext[2])
                    combinedExtents[3] = max(combinedExtents[3],ext[3])

                HAKOUv2utils.exportTrackPlotMapFile(os.path.join(outputDir,trackmap), shapefiles, 'EPSG:4326', extents, combinedExtents, tracklayername, settings.MAPSERVER_FONTS_LOCATION)

                mapURL = settings.MAPSERVER_URL + "?map=" + os.path.join(outputDir,trackmap)
                imageURL_QS = "&LAYERS="+tracklayername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in combinedExtents])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in combinedExtents])+"&WIDTH=700&HEIGHT=508"
                imageURL = mapURL + imageURL_QS

                # fetching track overview image
                self.updateState('Fetching track overview image')
                urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+os.path.join(outputDir,trackmap)+imageURL_QS, os.path.join(outputDir, trackoverview))

                wcc_settings.update({
                    'trackshapes': trackshapes,
                    'trackmap': trackmap,
                    'trackjson': trackjson,
                    'trackoverview': trackoverview,
                })

                d.update({
                    'track_map': os.path.join(outputDir, trackmap),
                    'track_mapURL': mapURL,
                    'track_imageURL': imageURL,
                    'track_layername': tracklayername,
                    'track_extents': combinedExtents,
                })

            self.updateState('Saving settings and outputs')
            saveOutputDict(os.path.join(outputDir,'outputDict.json'),d)
            saveOutputDict(os.path.join(outputDir,'settingsDict.json'),wcc_settings)

            # save README.txt
            with open(os.path.join(outputDir,'README.txt'),'w') as f:
                f.write(self.README_Contents)

        except:
            raise

    def processProbWave(self, level):
        """
        :param level:
        :return:
        """
        probdir = os.path.join(self.outputDir,'wave_%s'%level)
        if not os.path.exists(probdir) or self.force == True:
            if not os.path.exists(probdir):
                os.mkdir(probdir)

            self.updateState('Processing probability level %s wave results'%level)

            # convert to Geotiff
            HAKOUv3utils.convertProbabilityLevelToGeoTIFF(self.wavegeoprop, os.path.join(self.dataDir,'data'), os.path.join(self.resultDir,self.waveprobbin), os.path.join(probdir,self.wavegeotiff%level), level)

            # generate contours
            contour = WaveContourGenerator.WaveContourGenerator()
            contour.readGeoTIFF(os.path.join(probdir,self.wavegeotiff%level))
            contour.createPlot()

            # export to shapefile
            contour.exportShp(os.path.join(probdir,self.waveshape%level))
            MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',os.path.join(probdir, self.waveshape%level), os.path.join(probdir, self.waveshape%level),True)

            # export mapfile
            extents = MLSOperations.utils.getShpExtents(os.path.join(probdir, self.waveshape%level))
            contour.exportMapFile(os.path.join(probdir, self.wavemap%level), os.path.join(probdir, self.waveshape%level),'EPSG:4326', extents, self.wavelayername,settings.MAPSERVER_FONTS_LOCATION)

            mapURL = settings.MAPSERVER_URL + "?map=" + os.path.join(probdir,self.wavemap%level)
            imageURL_QS = "&LAYERS="+self.wavelayername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents])+"&WIDTH=700&HEIGHT=508"
            imageURL = mapURL + imageURL_QS

            # fetching wave overview image
            urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+os.path.join(probdir,self.wavemap%level)+imageURL_QS, os.path.join(probdir, self.waveoverview%level))
        else:
            self.updateState('Probability level %s wave results already exist'%level)

    def processProbSurge(self, level):
        """
        :param level:
        :return:
        """
        probdir = os.path.join(self.outputDir,'surge_%s'%level)
        if not os.path.exists(probdir) or self.force == True:
            if not os.path.exists(probdir):
                os.mkdir(probdir)

            self.updateState('Processing probability level %s surge results'%level)

            # convert to shapefile
            HAKOUv3utils.convertProbabilityLevelToShp(os.path.join(self.dataDir,'data'),os.path.join(self.resultDir,self.surgeprobbin),os.path.join(probdir,self.surgeshape%level),level)
            MLSOperations.utils.reprojectShp('EPSG:4326','EPSG:4326',os.path.join(probdir, self.surgeshape%level), os.path.join(probdir, self.surgeshape%level),True)

            # export mapfile
            extents = MLSOperations.utils.getShpExtents(os.path.join(probdir, self.surgeshape%level))
            HAKOUv3utils.exportSurgeMapFile(os.path.join(probdir, self.surgemap%level), os.path.join(probdir, self.surgeshape%level),'EPSG:4326', extents, self.surgelayername,settings.MAPSERVER_FONTS_LOCATION)

            mapURL = settings.MAPSERVER_URL + "?map=" + os.path.join(probdir,self.surgemap%level)
            imageURL_QS = "&LAYERS="+self.surgelayername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents])+"&WIDTH=700&HEIGHT=508"
            imageURL = mapURL + imageURL_QS

            # fetching surge overview image
            urllib.urlretrieve(settings.LOCAL_MAPSERVER_URL+"?map="+os.path.join(probdir,self.surgemap%level)+imageURL_QS, os.path.join(probdir, self.surgeoverview%level))
        else:
            self.updateState('Probability level %s surge results already exist'%level)

class processProbabilityLevels_v3(HAKOUModel_v3):
    """
    """
    def run(self, dataDir, resultDir, outputDir, levels, wave=True, surge=False, force=False):
        self.dataDir = dataDir
        self.resultDir = resultDir
        self.outputDir = outputDir
        self.force = force

        try:
            # calculate steps
            self.total = len(levels)*2 if wave and surge else len(levels)
            self.current = 0

            if wave == True:
                self.wavegeoprop = getGeoPropertiesFromJSON(os.path.join(resultDir,self.wavejson))

            # process each level
            for level in levels:
                if wave == True:
                    self.processProbWave(level)
                if surge == True:
                    self.processProbSurge(level)
        except:
            raise

def runTrackPlot(dataDir, outputDir, vector, latitude):
    cmd = [
        settings.HAKOU_TRACK_EXEC,
        '--',
        os.path.join(dataDir,'track'),
        outputDir,
        str(latitude),
        str(vector[0]),
        str(vector[1]),
        str(vector[3]),
        str(vector[5]),
    ]

    process = subprocess.Popen(' '.join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()
    dumpstr = process.stdout.readlines()
    errorstr = process.stderr.readlines()

    print ' '.join(cmd)
    print ''.join(dumpstr)
    print ''.join(errorstr)

    if process.returncode != 0:
        raise CalledProcessError(process.returncode, ' '.join(cmd), 'OUT:\n'+''.join(dumpstr)+'\n\nERROR:\n'+''.join(errorstr))

###############################################################################################
# Tasks for the wind model.
###############################################################################################
# DKDK this WINDModel_reprojectShp_hw() is the same with WINDModel_reprojectShp() - maybe reusable ?
def WINDModel_reprojectShp(fromProj, toProj, input, output, overwrite=False):
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

def WINDModel_getShpExtents(input):
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

def WINDModel_saveOutputDict(saveFile,dictObj):
    try:
        f = file(saveFile, 'w')
        json.dump(dictObj,f)
        f.close()

        return True
    except IOError:
        return False

class WINDModel_CalledProcessError(Exception):
    """Exception with the following properties:

    * returncode: Returned status code of the called process
    * cmd: The command that was run
    * output (optional or None): Output of the command
    """
    def __init__(self, returncode, cmd, output=None):
        self.returncode = returncode
        self.cmd = cmd
        self.output = output
        super(WINDModel_CalledProcessError, self).__init__(returncode,cmd,output)

class WINDModel_v1(Task):
    """
    Hurricane Wind Model to generate wind speed contour and vector fields

    """

    README_Contents = """Cybereye Hurricane Wind Model (Wind Speed Contour and Vector Fields) Results (cybereye.crc.nd.edu)

************* Wind Model Outputs - JSON files *************

results_hwind.json        - JSON document produced by the Hurricane Wind model with
                            information about the wind grid, input parameters, etc.
outputDict_hw.json        - JSON document with locations and names of all the
                            output file along with urls to access the MapServer
                            results. This is primarily for system use.
settingsDict_hw.json      - JSON document with settings used to produce these
                            results.


************* Wind Speed Contour Outputs *************

result_hwindspeed_shape.shp         - Shapefile (.shp, .dbf, .shx, & .prj) of the
                                      wind speed contour.
result_hwindspeed_map.map           - MapServer Mapfile for wind speed contour
result_hwindspeed_overview.png      - Overview image of the hurricane wind speed contour


************* Wind Vector Fields Outputs *************

results_hwindvector_geotiff.tiff    - Two bands (u,v) GeoTIFF of the wind vector fields
                                      Graphical representation of the gridded vector fields
                                      involving rasterization and interpolation
result_hwindvector_map.map          - MapServer Mapfile for wind vector fields.
result_hwindvector_overview.png     - Overview image of the hurricane wind vector fields"""


    def updateState(self, status='', extra='', progress=0.0, failed=False):
        """Updates the state of the celery task with a message and increments
        the count if increment is True.
        """
        state = 'PROGRESS'
        if failed:
            state = 'FAILURE'

        print '---------------------------------------------', status
        print self.update_state(state=state, meta={
            'progress': progress,
            'status': status,
            'extra': extra,
        })
        print '-------------------------------------------////'

    def run(self, dataDir, outputDir, vector, latitude, detonly=True, track=True):
        self.dataDir = dataDir
        self.resultDir = outputDir
        self.outputDir = outputDir
        self.force = False

        print '-------------------------------inside run'

        try:
            ### DKDK temporary check for execution time - no need to be copied
            #start_time = time.time()

            ### DKDK for hwindmodel, outputDir should have " " format, thus I made another variable, outputDir_hw here.
            ### DKDK In addition, slash in the end of outputDir should be included
            #outputDir_hw='"/home/vagrant/CyberEye01/Vagrant/trunk/django-WaveContourCalculator/DKDK/DK_hwind_tasks/"'
            outputDir_hw='"' +outputDir+ '/"'
            #print outputDir_hw

            hwindspeed_json='results_hwind.json'
            hwindspeed_geotiff='results_hwindspeed_geotiff.tiff'
            hwindspeed_shape='result_hwindspeed_shape.shp'
            ###below, hwindspeed_shape1, is for test purpose only - originally it should indicate hwindspeed_shape
            ###hwindspeed_shape1='result_hwindspeed_shape_repro.shp'
            hwindspeed_map = 'result_hwindspeed_map.map'
            hwindspeed_overview = 'result_hwindspeed_overview.png'
            hwindspeed_layername = 'hwindspeed_layer'


            ###
            ### DKDK below are variables for vector fields - start
            ###
            hwindvector_geotiff='results_hwindvector_geotiff.tiff'
            hwindvector_map = 'result_hwindvector_map.map'
            hwindvector_overview = 'result_hwindvector_overview.png'
            hwindvector_layername = 'hwindvector_layer'
            ###
            ### DKDK below are variables for vector fields - end
            ###


            ### DKDK change below for server setting
            DK_MAPSERVER_FONTS_LOCATION=settings.MAPSERVER_FONTS_LOCATION

            DK_MAPSERVER_URL = settings.MAPSERVER_URL
            ### DKDK LOCAL_MAPSERVER_URL used in the end of this file is basically the same with DK_MAPSERVER_URL
            ### DKDK in the original settings.py, LOCAL_MAPSERVER_URL == MAPSERVER_URL
            LOCAL_MAPSERVER_URL = settings.LOCAL_MAPSERVER_URL

            ### DKDK below is a temporary path to indicate hwindmodel executable file with path - need to be changed !!!
            hwind_V1_EXEC=settings.HWIND_V1_EXEC

            ### DKDK just mockup WCC case: actually below two, hwind_settings and d_hw, are not directly used for generating mapfile starting from json
            hwind_settings = {
            #    'dataDir': dataDir,
                'outputDir': outputDir,
            ### DKDK fontLocation defined here needs to be changed when integrating: settings.MAPSERVER_FONTS_LOCATION
                'fontLocation': DK_MAPSERVER_FONTS_LOCATION,
            #    'detonly': detonly,
            }

            ### DKDK below needs to be rearranged as vector[] would be different from WCC only - i.e., hwind would have more parameters
            ### DKDK this means that some WCC part may also need to be changed by accepting
            ### DKDK assuming vector[] with hwind to be vector=[latitude,longitude,angle,cp,dcp,vf,Rmax,time,dd,z0] here

            ### DKDK I think that at the def run(), we had better deliver an additional variable, 'vector_hw',
            ###      so that FSWP can be relying on its original variable, 'vector', and hwind will rely on 'vector_hw'
            ###      In this way, we will not need to change FSWP.

            ### DKDK 'time' variable is not used for the hwindmodel executable
            ### DKDK thus, vector_hw format should be:
            ### DKDK vector_hw=[latitude,longitude,angle,cp,dcp,vf,Rmax,dd,z0]
            # vector_hw=[21.256601342372, -158.12, 205, 950.0, 73.0, 16, 45.0, 0, 0.3]
            vector_hw=[vector[0],vector[1],vector[2],vector[3],vector[4],vector[5],vector[6],vector[7],vector[8]]
            print vector_hw
            ### DKDK "d_hw" is similar to "d" at the WCC
            d_hw = {
            ### DKDK below is hwind case with assumption of vector[] above
            ### DKDK I think that it is better to make a new variable, vector_hw, so that it cannot interrupt FSWP
                'latitude': vector_hw[0],
                'longitude': vector_hw[1],
                'angle': vector_hw[2],
                'cp': vector_hw[3],
                'dcp': vector_hw[4],
                'vf': vector_hw[5],
                'Rmax': vector_hw[6],
            ### DKDK 'time' is not used in the hwindmodel executable
                ###'time': vector_hw[7],
                'dd': vector_hw[7],
                'z0': vector_hw[8],

            ### DKDK MAPSERVER_BASELAYERS_URL defined here needs to be changed when integrating: settings.MAPSERVER_URL
                'MAPSERVER_BASELAYERS_URL': DK_MAPSERVER_URL,
            }


            ### DKDK Run Model for hwind
            self.updateState(status='Running model', progress=0.0)
            cmd_hw = [
            ### DKDK below need to be changed by adding a new variable to indicate the directory of executable file at the settings.py
                hwind_V1_EXEC,
                '--',
            ## DKDK below includes " "
                outputDir_hw,
            ### DKDK be carefule of vector_hw here
                ' '.join([str(i) for i in vector_hw]),
            ]
            print cmd_hw


            ### DKDK let's roll
            process_hw = subprocess.Popen(' '.join(cmd_hw), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            process_hw.wait()
            dumpstr_hw = process_hw.stdout.readlines()
            errorstr_hw = process_hw.stderr.readlines()

            print ' '.join(cmd_hw)
            print ''.join(dumpstr_hw)
            print ''.join(errorstr_hw)

            if process_hw.returncode != 0:
                raise WINDModel_CalledProcessError(process_hw.returncode, ' '.join(cmd_hw), 'OUT:\n'+''.join(dumpstr_hw)+'\n\nERROR:\n'+''.join(errorstr_hw))


            ###
            ### DKDK - Hurricane wind speed contour
            ###

            ### DKDK In this test, below function calls would not be used.
            ### DKDK basically, the new test is to call json output, contourf plot via matplotlib then read info from the contourf, making ESRI shapefile and then mapfile.
            ### DKDK that is, a process, json to geotiff, read geotiff to get data & interpolation, would not be used.

            # ### DKDK below function will combine read Json and create GeoTIFF
            # ### DKDK i.e., A combination of getGeoPropertiesFromJSON() and convertBinToGeoTIFF()/convertArrayToGeoTIFF()
            # ##self.updateState('Getting geographic extents of hurricane wind speed model outputs')
            # geoprop_hw, data_hw=HAKOUv3utils.FromJsonToGeoTIFF_hw(os.path.join(outputDir,hwindspeed_json),os.path.join(outputDir,hwindspeed_geotiff))
            # print geoprop_hw
            # print data_hw.shape
            #
            #
            # ### DKDK due to the combined function, DKfromJsontoGeoTIFF_hw(), below is no longer used
            # ###self.updateState('Converting to Geotiff')
            # ###HAKOUv3utils.convertBinToGeoTIFF(geoprop,os.path.join(dataDir,'data'),os.path.join(outputDir,wavedetbin),os.path.join(outputDir,wavegeotiff))
            # ###HAKOUv3utils.convertArrayToGeoTIFF_hw(geoprop_hw, data_hw, os.path.join(outputDir,hwindgeotiff))

            self.updateState(status='Generating hurricane wind speed contour', progress=0.1)
            contour_class_hw=WaveContourCalculator_hw.WaveContourGenerator()
            ###contour_class_hw.readGeoTIFF_hw(os.path.join(outputDir,hwindspeed_geotiff))
            contour_class_hw.createPlot_hw_short(os.path.join(outputDir,hwindspeed_json))


            self.updateState(status='Exporting & reprojecting hurricane wind speed shapefile', progress=0.2)
            contour_class_hw.exportShp_hw(os.path.join(outputDir,hwindspeed_shape))

            # ### DKDK I don't understand why reprojects the shapefile with the same input and output coord. and with the same input/output file name. It may be a test purpose for different coord.
            # ### DKDK MLSOperationsutils below should be renamed as MLSOperations.utils when combining
            WINDModel_reprojectShp('EPSG:4326','EPSG:4326',os.path.join(outputDir, hwindspeed_shape), os.path.join(outputDir, hwindspeed_shape),True)


            self.updateState('Exporting hurricane wind speed mapfile', progress=0.3)
            extents_hw = WINDModel_getShpExtents(os.path.join(outputDir, hwindspeed_shape))

            ### DKDK add one more parameter, LOCAL_MAPSERVER_URL - nono
            contour_class_hw.exportMapFile_hw(os.path.join(outputDir, hwindspeed_map), os.path.join(outputDir, hwindspeed_shape),'EPSG:4326', extents_hw, hwindspeed_layername,DK_MAPSERVER_FONTS_LOCATION)
            #contour_class_hw.exportMapFile_hw(os.path.join(outputDir, hwindspeed_map), os.path.join(outputDir, hwindspeed_shape),'EPSG:4326', extents_hw, hwindspeed_layername,DK_MAPSERVER_FONTS_LOCATION, LOCAL_MAPSERVER_URL)
            print extents_hw

            mapURL_hw = DK_MAPSERVER_URL + "?map=" + os.path.join(outputDir,hwindspeed_map)
            imageURL_QS_hw = "&LAYERS="+hwindspeed_layername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents_hw])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents_hw])+"&WIDTH=512&HEIGHT=512"
            imageURL_hw = mapURL_hw + imageURL_QS_hw
            print imageURL_hw

            ### DKDK fetching overview image for wind speed contour
            self.updateState(status='Fetching overview image for hurricane wind speed contour', progress=0.4)
            ### DKDK below seems that it generates hwindspeed_overview file(.png)
            urllib.urlretrieve(LOCAL_MAPSERVER_URL+"?map="+os.path.join(outputDir,hwindspeed_map)+imageURL_QS_hw, os.path.join(outputDir, hwindspeed_overview))


            ###
            ### DKDK - Hurricane wind vector fields
            ###
            ###
            ### DKDK hwindspeed does not produce geoprop anymore due to short workflow, but actually below geoprop_hwvec = geoprop_hw
            ###      Also, keys inside geoprop_hwvec dictionaries is changed to have suffix, _hwvec, just in case

            #geoprop_hwvec, u_vec, v_vec=HAKOUv3utils.FromJsonToGeoTIFF_hwvec(os.path.join(outputDir,hwindspeed_json),os.path.join(outputDir,hwindvector_geotiff))
            self.updateState(status='Generating GeoTIFF for hurricane wind vector fields', progress=0.6)
            ### DKDK FromJsonToGeoTIFF_hwvec() is moved to the class WaveContourCalculator_hw, which is referred by contour_class_hw
            geoprop_hwvec=contour_class_hw.FromJsonToGeoTIFF_hwvec(os.path.join(outputDir,hwindspeed_json),os.path.join(outputDir,hwindvector_geotiff))

            #elapsed_time = time.time() - start_time
            #print geoprop_hwvec, "geoprop_hwvec"
            #print u_vec.shape, v_vec.shape
            #print elapsed_time, "sec_hwvec"


            ### DKDK extents_hw would be obtained from hwindspeed so no need to repeat ~!
            #extents_hw = MLSOperationsutils.getShpExtents(os.path.join(outputDir, hwindspeed_shape))
            #print extents_hw

            ### DKDK a) contour_class_hw. in front of the function name can be the same with hwindspeed case as it only indicates CLASS name; b) fontlist is gone;
            ### DKDK c) hwindspeed_shape is used here (not vector) as it will be used to find path inside the function (not to use shapefile itself)
            # #self.updateState('Exporting mapfile for hurricane wind vector fields ')
            #contour_class_hw=WaveContourGenerator()
            self.updateState('Exporting mapfile for hurricane wind vector fields mapfile', progress=0.7)

            ### DKDK add one more parameter, LOCAL_MAPSERVER_URL - nono
            contour_class_hw.exportMapFile_hwvec(os.path.join(outputDir, hwindvector_map), os.path.join(outputDir, hwindvector_geotiff),'EPSG:4326', extents_hw, hwindvector_layername)
            #contour_class_hw.exportMapFile_hwvec(os.path.join(outputDir, hwindvector_map), os.path.join(outputDir, hwindvector_geotiff),'EPSG:4326', extents_hw, hwindvector_layername, LOCAL_MAPSERVER_URL)

            mapURL_hwvec = DK_MAPSERVER_URL + "?map=" + os.path.join(outputDir,hwindvector_map)
            imageURL_QS_hwvec = "&LAYERS="+hwindvector_layername+"&FORMAT=image%2Fpng&TRANSPARENT=true&ISBASELAYER=true&RESOLUTIONS=65536,32768,16384,8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1&MAXEXTENT="+'%2C'.join([str(e) for e in extents_hw])+"&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&SRS=EPSG%3A4326&BBOX="+','.join([str(e) for e in extents_hw])+"&WIDTH=512&HEIGHT=512"
            imageURL_hwvec = mapURL_hwvec + imageURL_QS_hwvec
            print imageURL_hwvec

            ### DKDK fetching wave overview image
            #self.updateState('Fetching overview image of hurricane wind vector fields')
            ### DKDK below seems that it generates hwindvector_overview file(.png)
            self.updateState(status='Fetching overview image of hurricane wind vector fields', progress=0.8)
            urllib.urlretrieve(LOCAL_MAPSERVER_URL+"?map="+os.path.join(outputDir,hwindvector_map)+imageURL_QS_hwvec, os.path.join(outputDir, hwindvector_overview))


            ### DKDK later, below two variables (used for output json files) need to be checked carefully
            ###      whether objects at the d_hw & hwind_settings are used/called in the web interfaces (html etc.)
            hwind_settings.update({
                ### DKDK wavebin is not used in the hwind
                #'wavebin': wavedetbin,
                ### DKDK below is no longer available in this test procedure
                #'hwindspeed_geotiff': hwindspeed_geotiff,
                'hwindspeed_shape': hwindspeed_shape,
                'hwindspeed_map': hwindspeed_map,
                'hwindspeed_json': hwindspeed_json,
                'hwindspeed_overview': hwindspeed_overview,
                ### DKDK below is no longer available in this test procedure, but may be able to be added if necessary - need some change
                #'hwindspeed_geoproperties': geoprop_hw,

                ###
                ### DKDK for hwindvector
                ###
                ### DKDK Note that shapefile is not generated at the hwindvector
                #'hwindvector_shape': hwindvector_shape,
                'hwindvector_map': hwindvector_map,
                ### DKDK json is the same with hwindspeed
                #'hwindspeed_json': hwindspeed_json,
                'hwindvector_overview': hwindvector_overview,
                ### DKDK geoproperties  is the same with hwindspeed
                #'hwindspeed_geoproperties': geoprop_hwvec,
                ### DKDK not quite sure if below is necessary but just in case
                'hwindvector_geotiff': hwindvector_geotiff,
            })

            d_hw.update({
                'hwindspeed_d_hw_map': os.path.join(outputDir, hwindspeed_map),
                'hwindspeed_d_hw_mapURL': mapURL_hw,
                'hwindspeed_d_hw_imageURL': imageURL_hw,
                'hwindspeed_d_hw_layername': hwindspeed_layername,
                'hwindspeed_d_hw_extents': extents_hw,

                ###
                ### DKDK for hwindvector
                ###
                'hwindvector_d_hw_map': os.path.join(outputDir, hwindvector_map),
                'hwindvector_d_hw_mapURL': mapURL_hwvec,
                'hwindvector_d_hw_imageURL': imageURL_hwvec,
                'hwindvector_d_hw_layername': hwindvector_layername,
                ### DKDK below is duplicate with hwindspeed
                #'hwindvector_d_hw_extents': extents_hw,

            })

            ### DKDK save output json files
            self.updateState(status='Saving settings and outputs', progress=0.9)
            WINDModel_saveOutputDict(os.path.join(outputDir,'outputDict_hw.json'),d_hw)
            WINDModel_saveOutputDict(os.path.join(outputDir,'settingsDict_hw.json'),hwind_settings)

            ### DKDK temporary time check - no need to be copied in the integration
            #elapsed_time = time.time() - start_time
            #print elapsed_time, "sec_total"

            with open(os.path.join(outputDir,'README.txt'),'w') as f:
                f.write(self.README_Contents)

            self.updateState(status='Done',progress=1.0)

        except:
            raise


tasks.register(WINDModel_v1)

tasks.register(contourModelCalc)
tasks.register(HAKOUModel_v2)
tasks.register(HAKOUModel_v3)
