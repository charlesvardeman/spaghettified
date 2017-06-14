__author__ = 'nathansmith'

import subprocess, os, sys, time
import urllib, urllib2
import simplejson

import MLSOperations
import MLSOperations.utils
import WaveContourGenerator
import WaveContourGenerator.utils

import utils as HAKOUv3utils
import HAKOUpy.v2.utils as HAKOUv2utils

try:
    import ipdb as pdb
except:
    import pdb


class Settings(object):
    MAPSERVER_FONTS_LOCATION = '/vagrant/djangoproject/djangoproject/data/fonts/font.list'
    MAPSERVER_URL = 'http://localhost:8080/cgi-bin/mapserv'
    MAPSERVER_BASELAYERS_URL = 'http://cybereye.crc.nd.edu/portal/cgi-bin/mapserv'
    HAKOU_EXEC = '/vagrant/HAKOUv3/KSM_CalculationZip_cmake/build/bin/KSM_CalculationZip'
    HAKOU_TRACK_EXEC = '/vagrant/packages/HAKOU/v2/trackplot/build/bin/MLE_TrackPlot'
    LOCAL_MAPSERVER_URL = 'http://localhost:80/cgi-bin/mapserv'

settings = Settings()


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
        simplejson.dump(dictObj,f)
        f.close()

        return True
    except IOError:
        return False

class HAKOUModel_v3():
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

    README_Contents = 'test'

    def updateState(self, status, *args, **kwargs):
        print status

    def run(self, dataDir, outputDir, vector, latitude, detonly=True, track=True):
        self.dataDir = dataDir
        self.resultDir = outputDir
        self.outputDir = outputDir
        self.force = False

        try:
            self.total = 2
            # if wave == True:
            #     self.total += 7
            #     if wavedetonly == True:
            #         self.total += 1
            # if surge == True:
            #     self.total += 4
            #     if surgedetonly == True:
            #         self.total += 1
            # if track == True:
            #     self.total += 3
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

            # Run Model
            self.updateState('Running model')
            cmd = [
                settings.HAKOU_EXEC,
                '--detonly' if detonly else '',
                '--',
                dataDir,
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
            HAKOUv3utils.convertBinToGeoTIFF(geoprop,dataDir,os.path.join(outputDir,wavedetbin),os.path.join(outputDir,wavegeotiff))

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
            HAKOUv3utils.convertBinToShp(dataDir,os.path.join(outputDir,surgedetbin),os.path.join(outputDir,surgeshape))
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
            HAKOUv3utils.convertProbabilityLevelToGeoTIFF(self.wavegeoprop, self.dataDir, os.path.join(self.resultDir,self.waveprobbin), os.path.join(probdir,self.wavegeotiff%level), level)

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
            HAKOUv3utils.convertProbabilityLevelToShp(self.dataDir,os.path.join(self.resultDir,self.surgeprobbin),os.path.join(probdir,self.surgeshape%level),level)
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



if __name__ == '__main__':
    datapath = '/vagrant/HAKOUv3/Hawaii_data/newData1.0/data/'
    resultspath = '/vagrant/HAKOUv3/KSM_CalculationZip_cmake/build/'

    task = HAKOUModel_v3()
    task.run(datapath, './', [-158.0997, 210, 950, 16, 40, 30], 21.256601342372, False, False)