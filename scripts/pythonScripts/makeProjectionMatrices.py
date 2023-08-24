import numpy as np
import pandas as pd
from itertools import combinations
import matplotlib.pyplot as plt
import seaborn as sns
import random

import sys
sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')
from incidenceToProjection import makeHiC_fromInc
from chains import makeIncDF_fromChainDists

## Set up
dataDir = '/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/2023_03_01_v0_dataGathering/v0_hypergraphSimulations/getMultiwayInteractions_fromBPChains/'
outDir = f'{dataDir}chains_10k_500_projectionMtxOutput/'

## Read in distance file
fileNum = sys.argv[1]
exMat = np.loadtxt(f'{dataDir}chains_500_10000_1500_1681171613/chain_dist_{fileNum}.txt')

print("Read in distance file, generating incidence matrix")

nrow = exMat.shape[0]
cutoff = int(sys.argv[2])

exChain = makeIncDF_fromChainDists(exMat,cutoff)
numReads = exChain.shape[1]

print("Making projection matrix")
hic_mat = makeHiC_fromInc(exChain)

card = exChain.sum()
maxCard = card.max(0)

report = [fileNum,numReads,maxCard]
report_string = '\t'.join(map(str,report))

print("Writing output")

file_path = f'{outDir}projMat_{fileNum}.txt'
np.savetxt(file_path, hic_mat, delimiter='\t', fmt='%d')

file_path = f'{outDir}summary.txt'
with open(file_path, 'a') as file:
    file.write(report_string + '\n')
