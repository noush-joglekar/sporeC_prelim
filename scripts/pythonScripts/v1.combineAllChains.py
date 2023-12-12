import numpy as np
import pandas as pd
import os.path
import pickle
from multiprocessing import Pool

import sys
sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')
from v1_chains import appendSingleBIncDict

def process_chunk(dataDir, args,chunk):
    """Process a chunk of files and write temporary output."""
    result_dict = {}

    start_file = (chunk - 1) * args.chunk_size + 1
    end_file = chunk * args.chunk_size

    for ix in range(start_file, end_file + 1):
        filePath = f'{args.dataDir}/{args.inputDir}/binConcatInc_{args.offDiagDist}_{args.prim_cutoff}_{args.sec_cutoff}_{ix}.pkl'
        if os.path.isfile(filePath):
            with open(filePath,'rb') as f:
                bIncDict = pickle.load(f)
                resultDict = appendSingleBIncDict(bIncDict,resultDict)

    # Write temporary output for the current chunk
    temp_output_file = f'{dataDir}/{args.outDir}hyperEdges_{args.offDiagDist}_{args.prim_cutoff}_{args.sec_cutoff}_chunk{chunk}_chains.pkl'
    with open(temp_output_file, 'wb') as f:
        pickle.dump(result_dict, f)

def constructFullDict_pkl_parallel(dataDir, args):
    """Process data in chunks using parallel processing."""
    num_chunks = args.numFiles // args.chunk_size

    # Create a multiprocessing pool with the number of available CPU cores
    with Pool(num_cores) as pool:
        pool.starmap(process_chunk, [(dataDir, args, chunk) 
                                     for chunk in range(1, num_chunks + 1)])

    # After processing all chunks, merge the results into a final output
    final_result_dict = {}

    for chunk in range(1, num_chunks + 1):
        temp_output_file = f'{dataDir}{args.outDir}hyperEdges_{args.offDiagDist}_{args.prim_cutoff}_{args.sec_cutoff}_chunk{chunk}_chains.pkl'
        with open(temp_output_file, 'rb') as f:
            result_dict = pickle.load(f)
            final_result_dict = appendSingleBIncDict(result_dict,final_result_dict)
        # Remove temporary files
        os.remove(temp_output_file)

    ## Have to prune here??

    # Write the final output
    with open(f'{dataDir}{args.outDir}hyperEdges_{args.offDiagDist}_{args.prim_cutoff}_{args.sec_cutoff}_final_chains.pkl', 'wb') as f:
        pickle.dump(final_result_dict, f)


def parse_args():
    parser = argparse.ArgumentParser(description="""Takes in a dict per chain, processes them
                                     in chunks, and finally generates a dict containing all the chains""")

    parser.add_argument("numFiles", type=int, help="Total number of files to process")
    parser.add_argument("chunk_size", type=int, help="Number of files in a chunk for processing")
    parser.add_argument("inputDir",type=str, help="Input directory containing files relative to data dir")
    parser.add_argument("outDir",type=str, help="Directory where you want your output relative to data dir")
    parser.add_argument("--prim_cutoff", type=int, help="Primary distance cutoff for neighbor detection",default=600)
    parser.add_argument("--sec_cutoff", type=int, help="Secondary distance cutoff for neighbor detection", default=750)
    parser.add_argument("--offDiagDist", type=int, help="How many nodes off-diag do you want to start", default=2)
    parser.add_argument("--num_cores", type=int, help="Number of cores for parallel processing", default=8)

    return parser.parse_args()

def main():
    ## Set up
    dataDir = '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/'
    args = parse_args()

    print("About to read in",args.numFiles,"files in chunks of",args.chunk_size)
    constructFullDict_pkl_parallel(dataDir, args)

if __name__ == "__main__":
    main()