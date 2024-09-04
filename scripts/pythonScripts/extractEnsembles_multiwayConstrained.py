import argparse
import pickle
import os
import pandas as pd

import sys
sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')

from v1_chains import IncDFCreator
from chains import dfToDict, chain_from_file, interactions_from_chains

def processFile(args, distMat, int1, int2, dataDir):

    creator = IncDFCreator(args.num_processes, args.prim_cutoff, args.sec_cutoff, args.offDiagDist)

    oneIter = creator.makeIncDF_fromChainDists_mp(distMat)
    exChain = oneIter[0]

    print("Calculating basic stats")
    numReads = exChain.shape[1]
    card = exChain.sum()
    maxCard = card.max(0)

    exChain.index = ["Bin"+str(i) for i in exChain.index]

    inc_dict = {}
    inc_dict = dfToDict(exChain,inc_dict)

    print("Checking if cell1 spec multiway contacts present")
    saveStatus = checkExclusivity(int1,int2,inc_dict)

    if saveStatus == "10":
        print("Yes! Writing output part 1")
        file_path = f'{dataDir}/{args.outDir}incDF_{args.offDiagDist}_{args.prim_cutoff}_{args.sec_cutoff}_{args.file_num}.pkl'
        exChain.to_pickle(file_path)

        file_path = f'{dataDir}/{args.outDir}incDict_{args.offDiagDist}_{args.prim_cutoff}_{args.sec_cutoff}_{args.file_num}.pkl'
        with open(file_path,'wb') as pklFile:
            pickle.dump(inc_dict,pklFile)

    report = [args.file_num, numReads, maxCard,
              args.prim_cutoff, args.sec_cutoff, 
              args.offDiagDist, saveStatus, args.cellType]
    
    report_string = '\t'.join(map(str,report))
    file_path = f'{dataDir}/{args.outDir}summary.txt'
    with open(file_path, 'a') as file:
        file.write(report_string + '\n')
              

def checkExclusivity(int1, int2, inc_dict):
    A = ["_".join(["Bin"+str(j) for j in i]) for i in int1]
    B = ["_".join(["Bin"+str(j) for j in i]) for i in int2]
    a = [x in inc_dict for x in A]
    b = [x in inc_dict for x in B]
    if any(a) and not any(b):
        return("10")
    elif any(a) and any(b):
        return("11")
    elif any(b) and not any(a):
        return("01")
    elif not any(a) and not any(b):
        return("00")
    

def createDistDF(args,dataDir):
    fileName = "%05d%s" % (args.file_num,".pts")
    df = pd.read_csv(os.path.join(dataDir, args.chainDir, fileName), sep="\s+", comment="#", header=None)
    chain = chain_from_file(df)
    distMat = interactions_from_chains(chain)
    return(distMat)

def parse_args():
    parser = argparse.ArgumentParser(description="""Takes in a dict per chain, processes them
                                     in chunks, and finally generates a dict containing all the chains""")

    #parser.add_argument("numFiles", type=int, help="Total number of files to process")
    parser.add_argument("file_num", type=int, help="Input file ID")
    parser.add_argument("chainDir",type=str, help="Input directory containing chain files relative to data dir")
    parser.add_argument("outDir",type=str, help="Directory where you want your output relative to data dir")
    parser.add_argument("--prim_cutoff", type=int, help="Primary distance cutoff for neighbor detection",default=600)
    parser.add_argument("--sec_cutoff", type=int, help="Secondary distance cutoff for neighbor detection", default=750)
    parser.add_argument("--offDiagDist", type=int, help="How many nodes off-diag do you want to start", default=2)
    parser.add_argument("--num_processes", type=int, help="Number of cores for parallel processing", default=8)
    parser.add_argument("--cellType", type=str, help="Cell type ID (A or B)", default='A')
    return parser.parse_args()

def main():
    ## Set up
    dataDir = '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/v3.multiwayConstraints/'
    args = parse_args()
    intList = [[[1,20,30], [5,9,15], [1,22,35,49], [5,11,20,25,40]],
               [[10,20,30], [30,43,50], [10,15,23,35,43,50]]]
    
    if args.cellType == 'A':
        int1 = intList[0]
        int2 = intList[1]
    elif args.cellType == 'B':
        int1 = intList[0]
        int2 = intList[1]

    distMat = createDistDF(args,dataDir)
    processFile(args, distMat, int1, int2, dataDir)
    return 


if __name__ == "__main__":
    main() 
