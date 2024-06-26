import importlib
import numpy as np
import argparse
import pickle
import os
import pandas as pd

import sys
sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')

from v1_chains import IncDFCreator
from chains import dfToDict, chain_from_file, interactions_from_chains

def processFile(args, distMat):
    creator = IncDFCreator(args.num_processes, args.prim_cutoff, args.sec_cutoff, args.offDiagLim)

    oneIter = creator.makeIncDF_fromChainDists_mp(distMat)
    exChain = oneIter[0]
    ratios = oneIter[1]

    print("Calculating basic stats")
    numReads = exChain.shape[1]
    card = exChain.sum()
    maxCard = card.max(0)

    exChain.index = ["Bin"+str(i) for i in exChain.index]

    inc_dict = {}
    inc_dict = dfToDict(exChain,inc_dict)

    print("Checking if cell1 spec multiway contacts present")
    

    print("Writing output part 1")
    file_path = f'{args.outDir}incDF_{args.offDiagLim}_{args.prim_cutoff}_{args.sec_cutoff}_{args.file_num}.pkl'
    exChain.to_pickle(file_path)

    file_path = f'{args.outDir}incDict_{args.offDiagLim}_{args.prim_cutoff}_{args.sec_cutoff}_{args.file_num}.pkl'
    with open(file_path,'wb') as pklFile:
        pickle.dump(bInc_dict,pklFile)

    report = [args.file_num, numReads, maxCard,
              args.prim_cutoff, args.sec_cutoff, args.offDiagLim]

def createDistDF(args,dataDir):
    df = pd.read_csv(os.path.join(dataDir, args.chainDir, args.file_num), sep="\s+", comment="#", header=None)
    chain = chain_from_file(df)
    distMat = interactions_from_chains(chain)
    return(distMat)

def parse_args():
    parser = argparse.ArgumentParser(description="""Takes in a dict per chain, processes them
                                     in chunks, and finally generates a dict containing all the chains""")

    #parser.add_argument("numFiles", type=int, help="Total number of files to process")
    parser.add_argument("file_num", type=int, help="Input file ID")
    parser.add_argument("chunk_size", type=int, help="Number of files in a chunk for processing")
    parser.add_argument("chainDir",type=str, help="Input directory containing chain files relative to data dir")
    parser.add_argument("outDir",type=str, help="Directory where you want your output relative to data dir")
    parser.add_argument("--prim_cutoff", type=int, help="Primary distance cutoff for neighbor detection",default=600)
    parser.add_argument("--sec_cutoff", type=int, help="Secondary distance cutoff for neighbor detection", default=750)
    parser.add_argument("--offDiagDist", type=int, help="How many nodes off-diag do you want to start", default=2)
    parser.add_argument("--num_processes", type=int, help="Number of cores for parallel processing", default=8)
    parser.add_argument("--readOfInterest", type=str, help="file containing info about read that you want to subset")
    parser.add_argument("--seed", type=int, help="Seed for random chain sampling", default=123)
    return parser.parse_args()

def main():
    ## Set up
    dataDir = '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/v3.multiwayConstraints/'
    args = parse_args()
    random.seed(args.seed)

    roi = [line.strip('\n').split(" ") for line in open(f'{dataDir}{args.outDir}{args.readOfInterest}')][0]

    print("About to read in",args.numFiles,"files in chunks of",args.chunk_size)
    constructFullDict_pkl_parallel(dataDir, roi, args)

if __name__ == "__main__":
    main() 
