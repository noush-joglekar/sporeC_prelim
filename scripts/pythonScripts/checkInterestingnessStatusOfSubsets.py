import argparse
import pickle
import numpy as np
import pandas as pd
import os
import csv

def main():
    dataDir = '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/v2.checkSubsetValidity/'
    args = parse_args()
    pklFile = [file for file in os.listdir(f'{dataDir}{args.outDir}') if ".pkl" in file][0]

    roi = [line.strip('\n').split(" ") for line in open(f'{dataDir}{args.outDir}{args.readOfInterest}')][0]
    card = len(roi[0].split("_"))

    edgeIx = preprocess(args, dataDir, pklFile, roi)
    writeOutput(args,dataDir,card,edgeIx,roi)

def preprocess(args, dataDir, pklFile, roi):
    """Pre-process the pickle file containing the subset chains to only include
    reads that were seen in at least 2 chains. Same as in the getMultiwayExpectedProbsByDist.py
    script"""

    with open(f'{dataDir}{args.outDir}/{pklFile}','rb') as f:
        hpEdges = pickle.load(f)

    print("Processing all the hyperedges from pickle file")
    hpKeys = [k for k in hpEdges.keys()]
    print("A total of",len(hpKeys),"initial interactions")

    readSupport = [v[0] for v in hpEdges.values()]
    chainSupport = [v[1] for v in hpEdges.values()]

    atLeastTwoChains = [i for i,x in enumerate(chainSupport) if x >=2]
    updatedDict = {hpKeys[i]:readSupport[i] for i in atLeastTwoChains}

    hpKeys = [k for k in updatedDict.keys()]
    print("A total of",len(hpKeys),"final interactions")

    edgeIx = hpKeys.index(roi[0])
    return(edgeIx)


def writeOutput(args,dataDir,card,edgeIx,roi):
    print("All done! Writing output to file")

    emp, cos = processDistanceMetrics(args,dataDir,card,edgeIx)
    splitSpecs = roi[3].strip("/").split("_")
    res = [emp,cos]
    res.extend(splitSpecs)

    with open(f'{dataDir}/{args.outFileName}','a') as f:
        writer = csv.writer(f, delimiter = "\t")
        writer.writerow(res)


def processDistanceMetrics(args,dataDir,card,edgeIx):
    print("Reading in distance metrics")
    coSimFile = pd.read_csv(f'{dataDir}{args.outDir}/dfs/cosineSim_card{card}.csv',sep = "\t")
    empDistFile = pd.read_csv(f'{dataDir}{args.outDir}/dfs/empDist_card{card}.csv',sep = "\t")
    empDistStatus = fixEmpCutoff(empDistFile)
    statusIx = empDistFile[empDistFile['Edge_ix'] == edgeIx].index[0]
    emp = empDistStatus[statusIx]
    cos = int(coSimFile[coSimFile['Edge_ix'] == edgeIx]['Status'].iloc[0])

    return([emp,cos])
    

def fixEmpCutoff(empDistFile):
    print("Re-calculating empDist status --- this takes a while")
    summary_empDist = empDistFile.filter(like="Sub").apply(lambda row: [np.mean(row), np.std(row)], 
                                                           axis=1, result_type='expand')
    summary_empDist.columns = ['mean','sd']
    empdistCutoff = getCutoff(summary_empDist,75)
    empdistStatus = [1 if x else 0 for x in (summary_empDist['mean'] >= empdistCutoff)]
    return(empdistStatus)

def getCutoff(summaryDF,quartile):
    q = f'{quartile}%'
    cutoff = pd.Series(summaryDF['mean']).describe()[q]
    return(cutoff)

def parse_args():
    parser = argparse.ArgumentParser(
        description="""Takes in a file containing a single line, which was used as input to
        v2.extractSubsetChainsToTestValidity_sample*.sh. It contains the necessary input to check
        against the interestingness output and report back """)

    parser.add_argument("outDir",type=str, help="Output directory relative to working dir")
    parser.add_argument("outFileName",type=str, help="Output file where summary file is written")
    parser.add_argument("--readOfInterest", type=str, help="file containing info about subset read")

    return parser.parse_args()

if __name__ == "__main__":
    main()
