import struct
import numpy as np
import sys

# path to both C++ output and standard Matlab output
# input for C++: par_s << 200, 0.4, 1.2; theta << 157.1002, 220, 930, 16, 40;
# High-fidelity data extracted from sc=1, cs=4. In test/sample folder
outputFilePath = "/Users/Cheng/Documents/CRC projects/Projects/Cyber-Eye/newMatlabCode/MLE_CalculationZip/result_surge_prob.bin"
testFilePath = "/Users/Cheng/Documents/CRC projects/Projects/Cyber-Eye/newMatlabCode/MLE_CalculationZip/test/sample/ResultsP_41.bin"
probabilisticFilePath = "/Users/Cheng/Documents/CRC projects/Projects/Cyber-Eye/newMatlabCodes2.0/KSM_CalculationZip/result_prob.bin"
deterministicFilePath = "/Users/Cheng/Documents/CRC projects/Projects/Cyber-Eye/newMatlabCodes2.0/KSM_CalculationZip/result_deter.bin"
binaryDataFilePath = "/Users/Cheng/Documents/CRC projects/Projects/Cyber-Eye/newMatlabCodes2.0/Hawaii_data/newData1.0/data"


def readDBinaryData(FilePath, unitbytes=4):
    """function to read binary data. First two elements are row(int, should be 1 for deterministic case, 1000 for
    probabilistic) and col(int), then the array of data, whose size is row*col(float=4 bytes) or
    row*col(double= 8 bytes)
    :param FilePath: Path to binary file
    :param unitbytes: Number of bytes per unit
    :return: Tuple of data, number of rows, number of cols
    """
    fin = open(FilePath,"rb")
    row = struct.unpack('i', fin.read(4))[0]
    col = struct.unpack('i', fin.read(4))[0]
    if unitbytes == 4:
        DeterministicData = struct.unpack('<' + str(row*col) + 'f', fin.read(unitbytes*row*col))
    else:
        DeterministicData = struct.unpack('<' + str(row*col) + 'd', fin.read(unitbytes*row*col))
        
    fin.close()
    return DeterministicData, int(row), int(col)


def testResultCompare(outputFilePath, testFilePath, threshold):
    """function to compare output with test data. If two dataset are the same,
    return (number of unmatched points, total points)

    :param outputFilePath: Path of output file
    :param testFilePath: Path of test file
    :param threshold: Threshold to use
    :return: tuple
    """
    # variable to record all the unmatch points 
    totalUnmatched = 0;
    # read C++ binary output file, get data array, row and col
    outputDataStructure = readDBinaryData(outputFilePath, 4)
    # read original matlab binary output file, get data array, row and col
    testDataStructure = readDBinaryData(testFilePath, 4)
    # compare two data
    for i in range(0 , outputDataStructure[1]*outputDataStructure[2]):
        # check if the difference is large than 0.001. If so, return 0 to alarm an error with the problematic index
        if abs(outputDataStructure[0][i] - testDataStructure[0][i]) > threshold*testDataStructure[0][i]:
            totalUnmatched += 1
    
    return totalUnmatched, outputDataStructure[1]*outputDataStructure[2], float(totalUnmatched)/float(outputDataStructure[1]*outputDataStructure[2])

# function for new data structure. First 1:Nsurge components are for surge, the rest is for wave.
# input should be the whole vector with both infomation
# dataFilePath should be the directory for binary data
def vectorChopper(binaryDataFilePath, inputVector):
    Nsurg = readDBinaryData(binaryDataFilePath + "/Nsurg.bin", 8)[0]
    surge = inputVector[0 : int(Nsurg[0])]
    wave = inputVector[int(Nsurg[0]):]
    return surge, wave

# function to get deterministic data for both surge and wave
# first one is surge, second is wave
def getDeterministicData(binaryDataFilePath, deterministicFilePath):
    Ress = np.array(readDBinaryData(deterministicFilePath, 4)[0])
    return vectorChopper(binaryDataFilePath, Ress)

# function to extract the right prob vector output
def getPercentageVector(probabilisticFilePath, percentage):
    percentage = 100 - percentage # get the adjusted percentage due to different order direction
    # read binary data
    fin = open(probabilisticFilePath,"rb")
    row = struct.unpack('i', fin.read(4))[0]
    col = struct.unpack('i', fin.read(4))[0]

    # get the col that corresponding to the percentage
    rowSelected = int(round(percentage * row/100)) - 1
    # make sure col index is not less than 0, just in case
    if rowSelected < 0:
        rowSelected = 0
    
    # calculate the offset. skip rowSelected rows + two ints (col and row)
    selectedVector = []
    for i in range(0, col):
        offset = (rowSelected + i * row )* 4 + 4 * 2 
        fin.seek(offset)
        selectedVector.append( struct.unpack('<1f', fin.read(4))[0] )
    fin.close()
    
    # return [length of the vector, which row is selected, vector] 
    return col, rowSelected, selectedVector


# function to extract corresponding percentage vector from binary output file.
# percentage should be a number [0, 100]
def generateSurgePlotPoints(binaryDataFilePath, deterministicFilePath, probabilisticFilePath=None, percentage=None):
    # load elev, tide, lgrid, nrem, aux_e from surge binary data set
    elev = readDBinaryData(binaryDataFilePath + "/elev.bin", 8)[0]
    tide = readDBinaryData(binaryDataFilePath + "/tide.bin", 8)[0][0]
    lgrid = readDBinaryData(binaryDataFilePath + "/lgrid.bin", 8)[0][0]
    nrem = readDBinaryData(binaryDataFilePath + "/nrem.bin", 8)[0]
    aux_e = readDBinaryData(binaryDataFilePath + "/aux_e.bin", 8)[0][0]
    # load MLE surge results from deterministic or probabilistic binary file
    if percentage is None:
        Ress = np.array(readDBinaryData(deterministicFilePath, 4)[0])
    else:
        Ress = np.array(getPercentageVector(probabilisticFilePath, percentage)[2])
    
    # Get the first part of the output vector as surge data
    Ress = vectorChopper(binaryDataFilePath, Ress)[0]
    
    #post correction of surge results by elev&tide infomation
    for i in range(0, len(Ress)):
        if Ress[i] > elev[i] + aux_e:
            Ress[i] = max(Ress[i], tide)
        else:
            Ress[i] = -999999
    
    # create storm_height array and fill it with surge info
    storm_height = -999999 * np.ones(lgrid)
    for i in range(0, len(nrem)):
        storm_height[nrem[i] - 1] = Ress[i]
    
    # load surge grid information
    node_lat = readDBinaryData(binaryDataFilePath + "/node_lat.bin", 8)[0]
    node_long = readDBinaryData(binaryDataFilePath + "/node_long.bin", 8)[0]
    node_elev = np.array(readDBinaryData(binaryDataFilePath + "/node_elev.bin", 8)[0]) * -1
    elem_nodes = readDBinaryData(binaryDataFilePath + "/elem_nodes.bin", 8)
    n_elem = elem_nodes[1]
    # save list into array 3*rows, int format
    elem_nodes = np.array([elem_nodes[0][0 : n_elem],
                           elem_nodes[0][n_elem : 2*n_elem],
                           elem_nodes[0][2*n_elem : 3*n_elem]], dtype = np.int32)
    n_node = len(node_lat)

    # start calculate contour line
    istart = 1
    long_line = []
    lat_line = []
    for j in range(0, n_elem):
        a1 = max(storm_height[elem_nodes[:, j] - 1] )
        a2 = min(storm_height[elem_nodes[:, j] - 1] )
        if a1 > -0.1 and a2 <= -100:
            temp_sh = storm_height[elem_nodes[:, j] - 1]
            i1 = np.argsort(temp_sh)
            a1 = np.array([temp_sh[i1[0]], temp_sh[i1[1]], temp_sh[i1[2]]])
            imax = elem_nodes[i1[2], j] - 1
            imin = elem_nodes[i1[0], j] - 1
            imed = elem_nodes[i1[1], j] - 1
            
            # calculate frac of the start point and its lat&long
            frac = (a1[2] - node_elev[imax])/(node_elev[imin] - node_elev[imax])
            frac = max(frac, 0)
            frac = min(frac, 1)
            # calculate lat&long of start point
            lat1 = (1-frac)*node_lat[imax] + frac*node_lat[imin]
            long1 = (1-frac)*node_long[imax] + frac*node_long[imin]
            
            # find end point based on med value of the triangle, calculate corresponding lat&long
            if a1[1] >= 0:
                inext = imin
                frac = (a1[1] - node_elev[imed])/(node_elev[imin]-node_elev[imed])
            else:
                inext = imax
                frac = 1 - (a1[2] - node_elev[imax])/(node_elev[imed]-node_elev[imax])
            frac = max(frac, 0)
            frac = min(frac, 1)
            # calculate lat&long of end point
            lat2 = (1-frac)*node_lat[imed] + frac*node_lat[inext]
            long2 = (1-frac)*node_long[imed] + frac*node_long[inext]
            
            # save start and end points into two arraies. one longitude, one latitude
            long_line.append([long1, long2])
            lat_line.append([lat1, lat2]) 
        
    return long_line, lat_line


if __name__ == '__main__':
    # example how to call these functions
    
    #print testResultCompare(outputFilePath, testFilePath, 0.01)
    #print readDBinaryData(testFilePath)[1]
    #print readDBinaryData(testFilePath)[2]
    #print readDBinaryData(testFilePath)[0][1]
    #print readDBinaryData(outputFilePath)[0][1]
    #print getPercentageVector(probabilisticFilePath, 5.0)[2][0:10]
    print generateSurgePlotPoints(binaryDataFilePath, deterministicFilePath, probabilisticFilePath, 5.0)[0][0: 10]
    #print vectorChopper(binaryDataFilePath, 'inputVector')
    #print getDeterministicData(deterministicFilePath)[0]
