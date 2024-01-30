import numpy as np
import pandas as pd
import pickle
import warnings
import argparse
import os
import pathlib
import json
from collections import Counter

import sys
sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')
from promethData_multiwayExpectedProbs import multiwayEval_realData, cutoffEval

def main():
    args = parse_args()

    if os.path.exists(f'{args.dataDir}/{args.workingDir}/{args.inputPkl}'):
        print("Loading in hyperedge pickle file")
        with open(f'{args.dataDir}/{args.workingDir}/{args.inputPkl}','rb') as f:
            hpEdges = pickle.load(f)
    else:
        print("Hyperedge pickle file does not exist, generating now")
        hpEdges = createPklFile(args)

    evaluateInterestingness(args, hpEdges)

    return

def evaluateInterestingness(args, hpEdges):

    plotDir = f'{args.dataDir}{args.workingDir}{args.plotDir}/'
    path1 = pathlib.Path(plotDir)
    path1.mkdir(parents=True, exist_ok=True)
    outDir = f'{args.dataDir}{args.workingDir}{args.outDir}/'
    path2 = pathlib.Path(outDir)
    path2.mkdir(parents=True, exist_ok=True)

    probHashOutName = f'{outDir}/{args.probHashOutFile}.json'

    print("Processing all the hyperedges from pickle file")
    hpKeys = [k for k in hpEdges.keys()]
    keyCard = [len(k.split("_")) for k in hpKeys]
    print("A total of",len(hpKeys),"initial interactions")

    readSupport = [v for v in hpEdges.values()]
#    cE = cutoffEval(keyCard,readSupport)
#    passedReadIx = cE.runForAllCards()
    atLeastTwoChains = [i for i,x in enumerate(readSupport) if x >=2]
    updatedDict = {hpKeys[i]:readSupport[i] for i in atLeastTwoChains} #passedReadIx

    hpKeys = [k for k in updatedDict.keys()]
    hpKeys_split = [k.split("_") for k in updatedDict.keys()]
    keyCard = [len(item) for item in hpKeys_split]

    print("Updated the input: retaining",len(hpKeys),
          "interactions with chain support of at least 2")

    print("Input args ------ toChoose:",args.sampleSize)
    print("Input args ------ toPlotRef:",args.plotRef)
    print("Input args ------ toPlotInd:",args.plotInd)
    print("Input args ------ toPlotScatter:",args.plotScatter)

    evalInstance = multiwayEval_realData(keyCard, updatedDict, hpKeys, hpKeys_split, args.seed,
                                args.sampleSize, args.plotRef, args.plotInd, args.plotScatter,
                                args.quartile, plotDir,outDir)

    if os.path.exists(probHashOutName):
        print("Expected probabilities file already exists...moving on")
        with open(probHashOutName,'r') as file:
            tmpHash = json.load(file)
            probHash = {key: {int(k): v for k, v in value.items()} for key, value in tmpHash.items() if value is not None}
    else:
        print("Expected probabilities file does not exist...creating now")
        probHash = evalInstance.makeAllReferenceHashDicts()
        with open(probHashOutName,'w') as file:
            json.dump(probHash,file)
        print("File created...moving on")

    warnings.filterwarnings('ignore')
    run = evalInstance.statsForAllReads(probHash)


def sort_key(item):
    return int(item.split('Bin')[1])

def createPklFile(args):
    chromSizes = pd.read_csv(f'{args.dataDir}hg38.chromSizes',sep="\t", names = ['chr','size']).set_index('chr')['size'].to_dict()
    readConcatemersWClosestGene = f'{args.dataDir}/{args.fragFile}'
    colnames = ["chr","start","end","readID","readLen","readQual",
                "geneChr","geneStart","geneEnd","strand","geneID","bioType","geneName","dist","ID"]
    fullBed = pd.read_csv(readConcatemersWClosestGene,sep = "\t",names = colnames)

    chrFile = fullBed[fullBed['chr']==args.chrom]
    binSize = 5*10**4 #1*10**6 #5*10**5
    chrBins = [x for x in range(0,chromSizes[args.chrom]+binSize,binSize)]
    chrFile_binned = pd.cut(chrFile['start'],bins = chrBins, 
                           labels = ["Bin"+str(i+1) for i in range(len(chrBins)-1)]).rename("binID")
    chrFile_wBinID = chrFile.merge(chrFile_binned,left_index=True,right_index=True)
    groupedBins = chrFile_wBinID.groupby('ID')['binID'].apply(list).reset_index(name='Bins')

    edges = ["_".join(sorted(list(set(a)), key=sort_key)) for a in groupedBins['Bins'] if 
             len(list(set(a))) > 1]
    
    hpEdges = dict(Counter(edges))

    with open(f'{args.dataDir}/{args.workingDir}/{args.inputPkl}','wb') as f:
        pickle.dump(hpEdges,f)
    
    return(hpEdges)


def parse_args():
    parser = argparse.ArgumentParser(
        description="""Takes in a pickle file of hyperedge frequencies and
        outputs .csv files with expected versus observed stats for a user-defined
        subset of reads. Also outputs plots of probability of all subset card reads
        by distance for all observed cards. Optional arguments to plot individual reads
        or summary statistics for all""")

    parser.add_argument("dataDir",type=str, help="Main data dir")
    parser.add_argument("workingDir",type=str, help="Working directory relative to data dir")
    parser.add_argument("plotDir",type=str, help="Plot directory relative to working dir")
    parser.add_argument("outDir",type=str, help="Output directory relative to working dir")
    parser.add_argument("fragFile",type=str, help="PoreC fragment file relative to working dir")
    parser.add_argument("chrom",type=str, help="PoreC fragment file relative to working dir")
    parser.add_argument("inputPkl", type=str, help="Input pickle file relative to workingDir")
    parser.add_argument("probHashOutFile", type=str, help="Name of probability reference hash")
    parser.add_argument("seed", type=int, help="Seed for read sampling")
    parser.add_argument("sampleSize", type=int, help="Number of reads to sample")
    parser.add_argument("quartile", type=int, help="Which quartile to set as cutoff", default = 25)
    parser.add_argument("--plotRef", action="store_true", help="Specify is you want a hist of dist stratified probabilities per card and subset")
    parser.add_argument("--plotInd", action="store_true", help="Specificy if you want an exp vs. obs of individual sampled reads")
    parser.add_argument("--plotScatter", action="store_true", help="Specify if you want a scatterplot of avg similarity measures for sampled reads")

    return parser.parse_args()

if __name__ == "__main__":
    main()
   
