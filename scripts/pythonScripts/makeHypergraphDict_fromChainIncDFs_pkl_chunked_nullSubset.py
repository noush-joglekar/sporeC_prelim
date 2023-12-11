import numpy as np
import pandas as pd
import os.path
import pickle
from multiprocessing import Pool

import sys
sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')
from chains import dfToDict

def getReads(inFile):
    """Read in a list of reads that we want to test our workflow on.
    Currently unsure if they are supposed to be random or interesting"""
    with open(inFile,'rb') as f:
        roi = pd.read_pickle(f)
        readsOfInterest = [r.split("_") for r in roi]
    return(readsOfInterest)

def checkIfValid(pickledFile,readsOfInterest,ix):
    """need to check diff between pd.read_pickle and pickle.load"""
    for roi in readsOfInterest:
        rois = pickledFile.loc[roi]
        nonZeroEdges = rois.columns[rois.sum() == rois.shape[0]]
        if len(nonZeroEdges) >0 : 
            print(f"File #{ix} is a valid file")
            return True

def process_chunk(dataDir, inputDir, outDir, offDiagDist, chunk, chunk_size, readsOfInterest):
    """Process a chunk of files and write temporary output."""
    result_dict = {}
    numEdges = []

    start_file = (chunk - 1) * chunk_size + 1
    end_file = chunk * chunk_size

    for ix in range(start_file, end_file + 1):
        print("File #",ix)
        filePath = f'{dataDir}/{inputDir}/binConcatInc_3_600_750_{ix}.pkl'
        if os.path.isfile(filePath):
            c = False
            bIncDF = pd.read_pickle(filePath)
            c = checkIfValid(bIncDF,readsOfInterest,ix)
            if c is True:
                result_dict = dfToDict(bIncDF, result_dict)
                nE = len(result_dict)
                numEdges.append(nE)

    # Write temporary output for the current chunk
    if result_dict:
        temp_output_file = f'{dataDir}/{outDir}hyperEdges_{offDiagDist}_600_750_chunk{chunk}_chains.pkl'
        with open(temp_output_file, 'wb') as f:
            pickle.dump(result_dict, f)

        temp_num_edges_file = f'{dataDir}{outDir}numEdges_{offDiagDist}_600_750_chunk{chunk}_chains.txt'
        np.savetxt(temp_num_edges_file, numEdges, delimiter='\t', fmt='%d')

    print(f"Chunk {chunk} done")

def constructFullDict_pkl_parallel(dataDir, inputDir, outDir, offDiagDist, numFiles, chunk_size, readsOfInterest, num_cores):
    """Process data in chunks using parallel processing."""
    num_chunks = numFiles // chunk_size

    # Create a multiprocessing pool with the number of available CPU cores
    with Pool(num_cores) as pool:
        pool.starmap(process_chunk, [(dataDir, inputDir, outDir, offDiagDist, chunk, chunk_size, readsOfInterest) 
                                     for chunk in range(1, num_chunks + 1)])

    # After processing all chunks, merge the results into a final output
    final_result_dict = {}
    final_num_edges = []

    for chunk in range(1, num_chunks + 1):
        temp_output_file = f'{dataDir}{outDir}hyperEdges_{offDiagDist}_600_750_chunk{chunk}_chains.pkl'
        if os.path.isfile(temp_output_file):
            with open(temp_output_file, 'rb') as f:
                result_dict = pickle.load(f)
                final_result_dict.update(result_dict)

            temp_num_edges_file = f'{dataDir}{outDir}numEdges_{offDiagDist}_600_750_chunk{chunk}_chains.txt'
            numEdges = np.loadtxt(temp_num_edges_file, dtype=int)
            final_num_edges.extend(numEdges)
            # Remove temporary files
            os.remove(temp_output_file)
            os.remove(temp_num_edges_file)

    # Write the final output
    if final_result_dict:
        with open(f'{dataDir}{outDir}hyperEdges_{offDiagDist}_600_750_subsetted_chains.pkl', 'wb') as f:
            pickle.dump(final_result_dict, f)

        final_num_edges_path = f'{dataDir}{outDir}numEdges_{offDiagDist}_600_750_subsetted_chains.txt'
        np.savetxt(final_num_edges_path, final_num_edges, delimiter='\t', fmt='%d')


## Set up.
dataDir = '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/'

numFiles = int(sys.argv[1])
offDiagDist = sys.argv[2]
inputDir = sys.argv[3]
outDir = sys.argv[4]
chunk_size = int(sys.argv[5])
readList = sys.argv[6]

print("Reading in reads/ edges of interest")
readsOfInterest = getReads(f'{dataDir}{readList}')
print("Need to subset files containing any from a set of", len(readsOfInterest), "random reads")

print("About to read in",numFiles,"files in chunks of",chunk_size)

constructFullDict_pkl_parallel(dataDir, inputDir, outDir, offDiagDist, numFiles, chunk_size,readsOfInterest,num_cores=4)
