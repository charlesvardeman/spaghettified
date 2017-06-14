"""
.. module:: utils
.. moduleauthor:: Nathan Smith <nsmith13@nd.edu>

"""

import subprocess
from exceptions import IOError

def reprojectShp(infile,reproj,fromProj="EPSG:4326",toProj="EPSG:26905"):
    """Reprojects a shapefile using ogr2ogr
    
    Args:
        | infile (str): File to reproject
        | reproj (str): Reprojected output file
    
    Kwargs:
        | fromProj (str): Original projection
        | toProj (str): New projection
    """
    
    cmd = "ogr2ogr -overwrite -s_srs %s -t_srs %s %s %s;" % (fromProj, toProj, reproj, infile)
    
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    dumpstr = process.stdout.readlines()
    errorstr = process.stderr.readlines()
    
    #print len(dumpstr)
    #print len(errorstr)
    #
    #print "".join(dumpstr)
    #print "".join(errorstr)
    
    if len(errorstr) > 0:
        raise IOError("".join(errorstr))