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

numFiles = int(sys.argv[1])
offDiagDist = sys.argv[2]
inputDir = sys.argv[3]

def constructFullDict(numFiles):
    """Takes in a directory of DFs and outputs a dict"""
    result_dict = {}
    numEdges = []
    for ix in range(1,numFiles+1):
        filePath = f'{dataDir}/{inputDir}/binConcatInc_{offDiagDist}_600_750_{ix}.pkl'
        if os.path.isfile(filePath):
            bIncDF = pd.read_pickle(filePath)
            result_dict = dfToDict(bIncDF,result_dict)
            nE = len(result_dict)
            numEdges.append(nE)
    return(result_dict,numEdges)


print("About to read in",numFiles,"files")

hpEdges, numEdges = constructFullDict(numFiles)

print("Constructed hyperedge dict - writing output")

with open(f'{dataDir}hyperEdges_{offDiagDist}_600_750_{numFiles}_chains.pkl', 'wb') as f:
    data = pickle.dump(hpEdges,f)

numEdgesPath = f'{dataDir}numEdges_{offDiagDist}_600_750_{numFiles}_chains.txt'
np.savetxt(numEdgesPath, numEdges, delimiter='\t', fmt='%d')
