import numpy as np
import pandas as pd
from itertools import combinations
import random
import os.path
import pickle

import sys
sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')
from chains import dfToDict, dictToDF

## Set up.
dataDir = '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/'

def constructFullDict(numFiles):
    """Takes in a directory of DFs and outputs a dict"""
    result_dict = {}
    numEdges = []
    for ix in range(1,numFiles+1):
        filePath = f'{dataDir}chains_10k_500_projectionMtxOutput/binConcatInc_600_1000_{ix}.pkl'
        if os.path.isfile(filePath):
            bIncDF = pd.read_pickle(filePath)
            result_dict = dfToDict(bIncDF,result_dict)
            nE = len(result_dict)
            numEdges.append(nE)
    return(result_dict,numEdges)

numFiles = int(sys.argv[1])

print("About to read in",numFiles,"files")

hpEdges, numEdges = constructFullDict(numFiles)

print("Constructed hyperedge dict - writing output")

with open(f'{dataDir}hyperEdges_{numFiles}_chains.pkl', 'wb') as f:
    data = pickle.dump(hpEdges,f)

numEdgesPath = f'{dataDir}numEdges_{numFiles}_chains.txt'
np.savetxt(numEdgesPath, numEdges, delimiter='\t', fmt='%d')
