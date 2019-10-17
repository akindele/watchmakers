#!/usr/bin/python

from watchmakers.load import *
from watchmakers.io_operations import *
from watchmakers.analysis import *
from watchmakers.sensitivity import *
from watchmakers.data import *

######################## Start of main function ###########################

if __name__ == "__main__":

    print(docstring)
    if arguments['-v']:
    	print(arguments)
#     print defaultValues
    print("")

    if arguments['-m']:
        generateMacrosNew(int(arguments['-N']),int(arguments['-e']))

    if arguments['-j']:
        generateJobsNew(int(arguments['-N']),arguments)

    if arguments['-M']:
        mergeNtupleFilesNew(arguments)

    if arguments['--histograms']:
        sensitivityMapPass2New()

    if arguments['--PMTAccHists']:
        EfficiencyMapInPMTVol()

    if arguments['--evalRate']:
        readEfficiencyHistogram()

    if arguments['--findRate']:
        findRate()

######################## Waba Luba Dub Dub!! ###########################
