#!/usr/bin/env python3
##
# @file         paramOpt.py
#
# @author       Jakub Chlebik \n
#               Faculty of Information Technology \n
#               Brno University of Technology \n
#               xchleb07@stud.fit.vutbr.com
#
# @brief        A python script trying to find the optimal evolution parameters by usign meta-evolution. 
#               UNFINISHED. Experimental and very computationaly expensiver. For use please study the implementation.
#
# @version      0.1
#
# @date         2019-04-01 (created) \n
#

# IMPORTS ######################################################################
from subprocess import Popen, PIPE
import sys
import os
from datetime import datetime
import argparse
import json
import random
import statistics

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.optimize as opt

# END OF IMPORTS ###############################################################

# METHODS DEFINITIONS ##########################################################

##
# @brief: Runs the specified optimizer with given config and returns results from stdout
# 
# @param [in] optimizerName         - name of the build optimizer binary
# @param [in] config                - dictionary of arguments values with which to run the optimizer (e.g. crossoverRate for GA, ...)
# @return string                    - result of the optimizer from stdout
def runOptimizer(optimizerName, config):
    binaryLaunch = ["./{0}".format(optimizerName)]
    optimizerArgs = [str(value) for _, value in config.items()]

    process = Popen(binaryLaunch + optimizerArgs, stdout=PIPE)
    (output, err) = process.communicate()
    exit_code = process.wait()

    if exit_code != 0:
        errorMessage = "Error while executing the program. "
        if (err != None):
            errorMessage += "stderr output : {0} ".format(err.decode("utf-8"))
        if (output != None):
            errorMessage += "stdout output: : {0} ".format(output.decode("utf-8"))
        raise RuntimeError(errorMessage)

    return output.decode("utf-8")
# END OF runOptimizer ############################################################

##
# @brief Creates a folder where all statistics will be saved
# 
# @param [in]   args                - given arguments to the script
# @param [in]   binaryName          - name of a binary thats doing the optimization
# @return String                    - path to the newly created folder
def createOutputFolder(args, binaryName):
    outputPath = args.outputFolder
    now = datetime.now().replace(microsecond=0)
    outputPathPrefix = os.path.join(outputPath, binaryName, now.strftime("%Y%m%d%H%M%S"))
    os.makedirs(outputPathPrefix, exist_ok=True)
    return outputPathPrefix
    #outputPathPrefix = "{}/{}/{}/{}".format(outputPath, now.isoformat(), binaryName, run)
# END OF createRunsBestOutputsPlot ################################################################

##
# @brief: Creates a python dictionary from the output of the optimization algorithm 
# 
# @param [in]       optimizerOutput     - text output from stdout of the optimizer
# @return JSONified data                - python dictionary object (["Launch"]["Run"]["Result"])
def createDictFromOptimzerOutput(optimizerOutput):
    sections = zip(["Launch", "Run", "Result"], optimizerOutput.split(sep="@@@"))
    result = dict()
    for section in sections:
        result[section[0]] = dict()
        sectionLines = section[1].strip().split(sep="\n")
        for sectionPart in sectionLines:
            for keyVal in sectionPart.split(sep="$"):
                if keyVal != '' and keyVal != '\n':
                    keyVal = keyVal.split(sep=":")
                    if section[0] == "Run":
                        if not keyVal[0] in result[section[0]]:
                            result[section[0]][keyVal[0]] = list()
                        split = keyVal[1].split(sep=",")
                        if len(split) > 1:
                            result[section[0]][keyVal[0]].append([abs(float(y)) for y in split])
                        else:
                            result[section[0]][keyVal[0]].append(abs(float(keyVal[1])))
                    elif section[0] == "Result":
                        split = keyVal[1].split(sep=",")
                        if len(split) > 1:
                            result[section[0]][keyVal[0]] = [abs(float(y)) for y in split]
                        else:
                            result[section[0]][keyVal[0]] = abs(float(keyVal[1]))
                    else:
                        result[section[0]][keyVal[0]] = keyVal[1]
    if not result["Run"]:
        result.pop("Run", None)
    if not result["Launch"] or not result["Result"]:
        raise KeyError("Launch or Result missing in optimizer output") 
    return result
# END OF createDictFromOptimzerOutput ############################################

##
# @brief: Loads the given configuration file for the optimizer and parses it into dictionary
# 
# @param [in]   configJsonPath      - path to the configuration file with arguments
# @return dictionary                - arguments to launch the optimizer with
def parseOptimizerConfiguration(configJsonPath):
    with open(os.path.abspath(configJsonPath), 'r') as jsonFile:
        optimizerArguments = json.load(jsonFile)
    return optimizerArguments
# END OF parseOptimizerConfiguration #############################################

##
# @brief: Wrapper around the fitness function. \n
#         Collects the data and runs the optimizer, using the result as a fitness for this chromosome
# 
# @param [in]   x               - the chromosome we wish to test
# @param [in]   indexes         - a list of indexes
# @param [in]   optimizerName   - name of the runned optimizer we are trying to find inputs for
# @param [in]   config          - configuration dictionary for the optimizer
# @return int                   - fitness value for the given input
def fitnessWrapper(x, indexes, optimizerName, config):
    fitnessWrapper.counter += 1
    i = 0
    for index in indexes:
        if x[i] > 1:
            config[index] = round(x[i])
        else:
            config[index] = round(x[i], 2)
        i = i + 1
    results = []
    for _ in range(0, 15):
        configDict["rngSeed"] = random.getrandbits(32)
        optimizerOutput = runOptimizer(binaryName, configDict)
        dataDictionary  = createDictFromOptimzerOutput(optimizerOutput)
        results.append(dataDictionary["Result"]["fitness"])
    res = statistics.median(results)
    print(res)
    return res
fitnessWrapper.counter = 0
# END OF fitnessWrapper #############################################

##
# @brief: Main entry point of the script.
# 
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-x', '--executable-path', help="path to the executable", required=True)
    parser.add_argument('-a', '--arguments', help="path to JSON file specifying the arguments needed to run the optimizer algorithm in correct order", required=True)
    parser.add_argument('-n', '--arguments-indexes', help="which of the parameters to optimize", nargs='+', type=int, required=True)
    parser.add_argument('-u', '--upper-bounds', nargs='+', type=float, required=True)
    parser.add_argument('-l', '--lower-bounds', nargs='+', type=float, required=True)  
 
    args = parser.parse_args()

    try:
        binaryName = args.executable_path
        configDict = parseOptimizerConfiguration(args.arguments)
        random.seed(configDict["rngSeed"])
        #outputPath = createOutputFolder(args, binaryName)

        #runs = int(args.runs)
        runs = 15
        
        xKey = []
        xVal = []
        b = opt.Bounds(args.lower_bounds, args.upper_bounds)
        b_tuples = []
        for i in range(0, len(args.lower_bounds)):
            b_tuples.append((args.lower_bounds[i], args.upper_bounds[i]))

        i = 0
        for key, value in configDict.items():
            if i in args.arguments_indexes:
                xKey.append(key)
                xVal.append(value)

            i = i + 1 
        xVal = np.array(xVal)
        results = []
        for i in range(0, runs):
            # res = opt.minimize(fitnessWrapper, \
            #         xVal, \
            #         args=(xKey, binaryName, configDict), \
            #         method='nelder-mead', \
            #         bounds=b, \
            #         options={'maxiter': 300, 'disp': True}    )

            # res = opt.differential_evolution(\
            #        fitnessWrapper, \
            #         bounds=b, \
            #         args=(xKey, binaryName, configDict), \
            #         strategy='randtobest1bin', \
            #         maxiter=200, \
            #         seed=random.getrandbits(32), \
            #         disp=True, \
            #         #workers=-1, \
            #         atol=1e-3)
            res = opt.dual_annealing(\
                    fitnessWrapper, \
                    bounds=b_tuples, \
                    args=(xKey, binaryName, configDict), \
                    seed=random.getrandbits(32), \
                    )
            results.append(res.x)
        print(results)

    except Exception as err:
        print(err)
# END of MAIN ###################################################################


if __name__ == "__main__":
    main()