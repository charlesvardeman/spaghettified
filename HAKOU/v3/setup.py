import argparse
import os, sys

import subprocess

THISDIR = os.path.abspath(os.path.dirname(__file__))
EXEC_FILENAME = 'KSM_CalculationZip'

try:
    import ipdb as pdb
except:
    import pdb

def main():
    parser = argparse.ArgumentParser(
                    description='Sets up this version of the HAKOU model'
                )
    parser.add_argument('outfile', nargs='?', type=str, help='Location to put the executable')
    parser.add_argument('--rebuildjson', action='store_true', default=False, help='Whether to reguild the jsoncpp library')

    args = parser.parse_args()

    #pdb.set_trace()

    # build jsoncpp
    if not os.listdir(os.path.join(THISDIR,'source','json')) or args.rebuildjson == True:
        cmd = 'cd %s; sh jsoncpp_0.5.0_build.sh' % os.path.join(THISDIR,'source')
        subprocess.call(cmd, shell=True)

    # build executable
    cmd = 'cd %s; sh build.sh' % os.path.join(THISDIR,'source')
    subprocess.call(cmd, shell=True)

    # copy to the desired location if available
    if args.outfile != None:
        cmd = 'cp %s %s' % (os.path.join(THISDIR,'source','build','bin',EXEC_FILENAME), args.outfile)
        subprocess.call(cmd, shell=True)




if __name__ == '__main__':
    main()
