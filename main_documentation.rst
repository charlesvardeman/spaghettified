.. django-WaveContourCalculator documentation master file, created by
   sphinx-quickstart on Fri May 11 13:47:33 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to django-WaveContourCalculator's documentation!
========================================================

This tool provides estimates for surge and (possibly) waves for different U.S coastal regions. It is based on surrogate
modeling informed by existing databases of ADCIRC (for surge) and SWAN (for waves) high-fidelity simulations. It
provides the opportunity to estimate surge/wave for a storm with specific characteristics (deterministic assessment)
but also allows for a probabilistic assessment for any storm prior to landfall considering statistical errors related
to the storm track predictions. Batch upload of different storm tracks will be enabled in the future.

User should first select a specific coastal region from the respective dropdown menu and then should define the
characteristics of the storm. These characteristics correspond to the following properties defined for the point that
the storm is projected to make (or has made) landfall (when there is no actual landfall point, a conventional one can
be chosen)

# Latitude and longitude (in degrees)
# Angle of final approach (in degrees defined clockwise from south)
# Storm speed (in knots)
# Central pressure (in mbar)
# Radius of maximum winds (in km)

The first two can be defined either though numerical inputs or by using the map to set location/angle. For the
location/angle, the user can first set the landfall location (latitude and longitude) in the map (first point selection)
and then define the direction of the storm (second point selection). The storm speed, radius of maximum winds and
central pressure can be defined either numerically or through the slide bar next to them. Each input is constrained
within the bounds available in each database.

For a probabilistic assessment of storms prior to landfall, the "Time to Landfall" (in hours) needs to be also
numerically specified. This is used to estimate statistical errors related to the predicted storm track.

As soon as the input is defined, the storm track is plotted along with the cone of possible tracks if the time till
landfall is also provided. The output (surge and possibly waves) can be then calculated for a storm with the specified
characteristics. As soon as the calculation is performed, the results can be visualized on screen and also downloaded.
If a probabilistic assessment is requested, a Monte Carlo simulation is performed for different storms within the cone
of possible tracks and then the output with different probabilities of exceedance can be plotted. Two different
probabilities are set by default, 10% and 50%, but the user is also given the opportunity to select different probabilities.

All runs performed can be managed from the manage tab. From there the results can be also directly downloaded.

.. image:: /images/wcc.png
   :align: center
   :width: 700

------------
Dependencies
------------

This app requires the following packages:

* MLSOperations on `redmine <https://redmine.crc.nd.edu/svn/cybereye/packages/trunk/packages/MLSOperations>`_
* WaveContourGenerator on `redmine <https://redmine.crc.nd.edu/svn/cybereye/packages/trunk/packages/WaveContourGenerator>`_
* Celery
* MatPlotLib

Images are generated with MapServer MapFiles, so MapServer is required to display them.


Contents:
=========

.. toctree::
   :maxdepth: 2
   
   api/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

