.. MLSOperations documentation master file, created by
   sphinx-quickstart on Wed May  9 15:53:25 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MLSOperations's documentation!
=========================================

.. toctree::
   :maxdepth: 2

The Moving Least Squares (mls) model was originally written in MATLAB as a proof
of concept.  For the CyberEye project it was rewritten in C++.  It uses the
Eigen 3.0.3 C++ template library for linear algebra.  The mls C++ code has been
wrapped using Swig and exposed as a Python module for ease of use.  The C++ code
produces two csv files, which are then converted into geotiffs and reprojected
into EPSG:26905 from ESPG:4326 using the GDAL library.  The projection EPSG:26905
is used for the island of Hawaii, which is the area of focus for CyberEye during
development.

This package is currently being used as a part of several others as web apps.  In
the final version of CyberEye it will be run on a periodic basis as new input
data is generated from external sources.

-----
Usage
-----

#. Define the initial_parameters to use for the Moving Least Squares analysis::

    dataDir = '../data/'
    outputDir = '../outputs/'
    initial_parameters = [
        158.4,
        180,
        950,
        12,
        50,
        dataDir,
        outputDir,
        'results_Hsin.csv',
        'results_Hsout.csv',
    ]

#. Define the geoproperties as explained in :py:meth:`MLSOperations.getStandardGeoProperties`

#. Run the moving least squares analysis::

    MLSOperations.runMLS(*initial_parameters)

#. Convert the resulting csv files to geotiff::

    MLSOperations.convertCSVToGeoTIFF(result_geoproperties['Hsin'], outputDir+'results_Hsin.csv', outputDir, 'resultsHsin.tif')
    MLSOperations.convertCSVToGeoTIFF(result_geoproperties['Hsout'], outputDir+'results_Hsout.csv', outputDir, 'resultsHsout.tif')

#. Reproject the geotiffs to the correct projection::

    MLSOperations.utils.reprojectGeoTiff('EPSG:4326','EPSG:26905', outputDir+'resultsHsin.tif', outputDir+'resultsHsin_reprojected.tif')
    MLSOperations.utils.reprojectGeoTiff('EPSG:4326','EPSG:26905', outputDir+'resultsHsout.tif', outputDir+'resultsHsout_reprojected.tif')

-----------------
Class Definitions
-----------------

`Import Graph <importGraph.pdf>`_

.. automodule:: MLSOperations
    :members:

.. automodule:: MLSOperations.utils
    :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

