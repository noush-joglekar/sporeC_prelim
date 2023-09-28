
#!/bin/python

import numpy as np
import pickle
import warnings
import argparse

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
    seed = args.seed
    toChoose = args.sampleSize
    toPlotRef = args.toPlotRef
    toPlotInd = args.toPlotInd
    toPlotScatter = args.toPlotScatter

    with open(f'{dataDir}/{inputPkl}','rb') as f:
        hpEdges = pickle.load(f)

    print("Processing all the hyperedges from pickle file")
    hpKeys_split = [k.split("_") for k in hpEdges.keys()]
    keyCard = [len(item) for item in hpKeys_split]

    evalInstance = multiwayEval(keyCard, hpEdges, hpKeys_split, seed, 
                                toChoose, toPlotRef, toPlotInd, toPlotScatter,
                                plotDir,outDir)
    probHash = multiwayEval.makeAllReferenceHashDicts()
    
    warnings.filterwarnings('ignore')
    evalInstance.statsForAllReads(self,probHash)

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
    parser.add_argument("seed", type=int, help="Seed for read sampling")
    parser.add_argument("sampleSize", type=int, help="Number of reads to sample")
    parser.add_argument("toPlotRef", type=bool, help="Plot dist stratified probabilities per card and subset",default=True)
    parser.add_argument("toPlotInd", type=bool, help="Plot exp vs. obs of individual sampled reads", default=False)
    parser.add_argument("toPlotScatter", type=bool, help="Plot scatterplot of avg similarity measures for sampled reads",default=False)

    return parser.parse_args()

if __name__ == "__main__":
    main()