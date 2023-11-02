#!/bin/python

import numpy as np
import pandas as pd
import argparse

import sys
sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')
from incidenceToProjection import makeHiC_fromInc
from chains import IncDFCreator, increaseIncDF_binSize


def main():
    ## Set up
    dataDir = '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/'
    args = parse_args()

    fileNum = args.file_num
    outDir = f'{dataDir}{args.outDir}/'

    ## Read in distance file
    exMat = np.loadtxt(f'{dataDir}{args.inputDir}/chain_dist_{fileNum}.txt')
    print("Read in distance file, generating incidence matrix")

    prim_cutoff = args.prim_cutoff
    sec_cutoff = args.sec_cutoff
    numProcesses = args.num_processes
    offDiagLim = args.offDiagLim

    creator = IncDFCreator(numProcesses, prim_cutoff, sec_cutoff, offDiagLim)
    exChain = creator.makeIncDF_fromChainDists_mp(exMat)

    numReads = exChain.shape[1]
    exChain_by5 = increaseIncDF_binSize(exChain,5)

    print("Making projection matrix")
    hic_mat = makeHiC_fromInc(exChain)

    card = exChain.sum()
    maxCard = card.max(0)

    report = [fileNum,numReads,maxCard,prim_cutoff, sec_cutoff,offDiagLim]
    report_string = '\t'.join(map(str,report))

    print("Writing output")

    file_path = f'{outDir}projMat_{offDiagLim}_{prim_cutoff}_{sec_cutoff}_{fileNum}.txt'
    np.savetxt(file_path, hic_mat, delimiter='\t', fmt='%d')

    file_path = f'{outDir}incDF_{offDiagLim}_{prim_cutoff}_{sec_cutoff}_{fileNum}.pkl'
    exChain.to_pickle(file_path)

    file_path = f'{outDir}binConcatInc_{offDiagLim}_{prim_cutoff}_{sec_cutoff}_{fileNum}.pkl'
    exChain_by5.to_pickle(file_path)

    file_path = f'{outDir}summary.txt'
    with open(file_path, 'a') as file:
        file.write(report_string + '\n')


def parse_args():
    parser = argparse.ArgumentParser(description="""Takes in a distance matrix and gives 3 outputs:
        1. An incidence dataframe,
        2. the incDF but binned to have 5 nodes combined into 1,
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
