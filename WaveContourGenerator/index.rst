.. WaveContourGenerator documentation master file, created by
   sphinx-quickstart on Wed May  9 14:36:51 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to WaveContourGenerator's documentation!
================================================

.. toctree::
   :maxdepth: 2

The Wave Height Contour Generator is currently implemented as a web app that produces
contour plots from the mls output.  The python module behind the app accepts a
geotiff in EPSG:4326.  It uses the MatPlotLib plotting library
(http://matplotlib.sourceforge.net/) contourf module to generate filled contour
plots of comparable quality to MATLAB.

These plots are exported as shapefiles and reprojected to EPSG:26905.  Along with
the shapefiles a MapServer (http://mapserver.org/) mapfile is generated with the
coloring and labels of the various levels.  The plots are presented using MapServer
and OpenLayers (http://openlayers.org/).

The web app version of this process is only for testing purposes and in the final
version will be run on a periodic basis as new input data is generated from external sources.

-----
Usage
-----

#. Create WaveContourGenerator object::

    contour = WaveContourGenerator.WaveContourGenerator()

#. Get data from GeoTIFF file::

    contour.readGeoTIFF(outputDir+mlsgeotiff1)

#. Create the plot and export the shapefile::

    contour.createPlot(levels)
    contour.exportShp(os.path.join(outputDir,outFile1))

#. Reproject the resulting shapefile and generate a MapFile::

    WaveContourGenerator.utils.reprojectShp(os.path.join(outputDir, outFile1), os.path.join(outputDir, reprojectedFile1))
    contour.exportMapFile(os.path.join(outputDir, mapFile1), os.path.join(outputDir, reprojectedFile1),fontLocation)

-----------------
Class Definitions
-----------------

.. autoclass:: WaveContourGenerator.WaveContourGenerator
    :members:

.. automodule:: WaveContourGenerator.utils
    :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

