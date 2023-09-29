
#!/bin/python

import numpy as np
import pickle
import warnings
import argparse
import os
import json

import sys
sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')
from multiwayExpectedProbs import multiwayEval

def main():
    dataDir = '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/'

    args = parse_args()
    inputPkl = args.inputPkl
    workingDir = args.workingDir
    plotDir = f'{dataDir}{workingDir}{args.plotDir}/'
    outDir = f'{dataDir}{workingDir}{args.outDir}/'
    probHashOutName = f'{outDir}/{args.probHashOutFile}.json'
    seed = args.seed
    toChoose = args.sampleSize
    print("Input args ------ toChoose:",toChoose)
    toPlotRef = args.plotRef
    print("Input args ------ toPlotRef:",toPlotRef)
    toPlotInd = args.plotInd
    print("Input args ------ toPlotInd:",toPlotInd)
    toPlotScatter = args.plotScatter
    print("Input args ------ toPlotScatter:",toPlotScatter)

    print("Loading in hyperedge pickle file")
    with open(f'{dataDir}/{inputPkl}','rb') as f:
        hpEdges = pickle.load(f)

    print("Processing all the hyperedges from pickle file")
    hpKeys_split = [k.split("_") for k in hpEdges.keys()]
    keyCard = [len(item) for item in hpKeys_split]

    evalInstance = multiwayEval(keyCard, hpEdges, hpKeys_split, seed, 
                                toChoose, toPlotRef, toPlotInd, toPlotScatter,
                                plotDir,outDir)
    
    if os.path.exists(probHashOutName):
        print("Expected probabilities file already exists...moving on")
        with open(probHashOutName,'r') as file:
            tmpHash = json.load(file)
            probHash = {key: {int(k): v for k, v in value.items()} 
                        for key, value in tmpHash.items()}
    else:
        print("Expected probabilities file does not exist...creating now")
        probHash = evalInstance.makeAllReferenceHashDicts()
        with open(probHashOutName,'w') as file:
            json.dump(probHash,file)
        print("File created...moving on")
    
    warnings.filterwarnings('ignore')
    evalInstance.statsForAllReads(probHash)

def parse_args():
    parser = argparse.ArgumentParser(
        description="""Takes in a pickle file of hyperedge frequencies and
        outputs .csv files with expected versus observed stats for a user-defined
        subset of reads. Also outputs plots of probability of all subset card reads
        by distance for all observed cards. Optional arguments to plot individual reads
        or summary statistics for all""")

    parser.add_argument("workingDir",type=str, help="Working directory relative to data dir")
    parser.add_argument("plotDir",type=str, help="Plot directory relative to working dir")
    parser.add_argument("outDir",type=str, help="Output directory relative to working dir")
    parser.add_argument("inputPkl", type=str, help="Input pickle file relative to data dir containing reads as hyperedges")
    parser.add_argument("probHashOutFile", type=str, help="Name of probability reference hash")
    parser.add_argument("seed", type=int, help="Seed for read sampling")
    parser.add_argument("sampleSize", type=int, help="Number of reads to sample")
    parser.add_argument("--plotRef", action="store_true", help="Specify is you want a hist of dist stratified probabilities per card and subset")
    parser.add_argument("--plotInd", action="store_true", help="Specificy if you want an exp vs. obs of individual sampled reads")
    parser.add_argument("--plotScatter", action="store_true", help="Specify if you want a scatterplot of avg similarity measures for sampled reads")

    return parser.parse_args()

if __name__ == "__main__":
    main()