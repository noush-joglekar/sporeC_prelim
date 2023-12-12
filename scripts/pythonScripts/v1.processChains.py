#!/bin/python

import numpy as np
import pandas as pd
import argparse

import sys
sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')

from incidenceToProjection import makeHiC_fromInc
from v1.chains import IncDFCreator, increaseIncDF_binSize, dfToDict

def main():
    ## Set up
    dataDir = '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/'
    args = parse_args()

    outDir = f'{dataDir}{args.outDir}/'

    ## Read in distance file
    exMat = np.loadtxt(f'{dataDir}{args.inputDir}/chain_dist_{args.file_num}.txt')
    print("Read in distance file, generating incidence matrix")

    creator = IncDFCreator(args.num_processes, args.prim_cutoff, args.sec_cutoff, args.offDiagLim)
    
    oneIter = creator.makeIncDF_fromChainDists_mp(exMat)
    exChain = oneIter[0]
    ratios = oneIter[1]

    exChain.index = ["Bin"+str(i) for i in range(exChain.shape[0])]

    numReads = exChain.shape[1]
    print("Binning into sets of 5")
    bInc_by5 = increaseIncDF_binSize(exChain,5)

    print("Making projection matrix")
    hic_mat = makeHiC_fromInc(exChain)

    print("Converting a binned incDF to dict")
    newCols = [str(bInc_by5.columns[i])+"_"+str(round(ratios[i],2)) 
           for i in range(bInc_by5.shape[1])]
    bInc_by5.columns = newCols

    bInc_dict = {}
    bInc_dict = dfToDict(bInc_by5,bInc_dict)

    card = exChain.sum()
    maxCard = card.max(0)

    report = [args.file_num, numReads, maxCard,
              args.prim_cutoff, args.sec_cutoff, args.offDiagLim]
    report_string = '\t'.join(map(str,report))

    print("Writing output")

    file_path = f'{outDir}projMat_{args.offDiagLim}_{args.prim_cutoff}_{args.sec_cutoff}_{fileNum}.txt'
    np.savetxt(file_path, hic_mat, delimiter='\t', fmt='%d')

    file_path = f'{outDir}incDF_{args.offDiagLim}_{args.prim_cutoff}_{args.sec_cutoff}_{fileNum}.pkl'
    exChain.to_pickle(file_path)

    file_path = f'{outDir}binConcatInc_{args.offDiagLim}_{args.prim_cutoff}_{args.sec_cutoff}_{fileNum}.pkl'
    bInc_dict.to_pickle(file_path)

    file_path = f'{outDir}summary.txt'
    with open(file_path, 'a') as file:
        file.write(report_string + '\n')

def parse_args():
    parser = argparse.ArgumentParser(description="""Takes in a distance matrix and generates:
        1. An incidence dataframe,
        2. Binned incidence df saved to dict,
        3. A projection matrix""")

    parser.add_argument("inputDir",type=str, help="Input directory relative to data dir")
    parser.add_argument("outDir",type=str, help="Output directory relative to data dir")
    parser.add_argument("file_num", type=int, help="Input file ID")
    parser.add_argument("prim_cutoff", type=int, help="Distance cutoff value to make sub-matrices")
    parser.add_argument("sec_cutoff", type=int, help="Distance cutoff value to deem adjacency")
    parser.add_argument("num_processes", type=int, help="Number of parallel processes")
    parser.add_argument("--offDiagLim", type=int, help="How many nodes off-diag do you want to start", default=2)

    return parser.parse_args()

if __name__ == "__main__":
    main()