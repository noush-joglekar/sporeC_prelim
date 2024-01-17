from itertools import combinations
import numpy as np
import pandas as pd

def sort_key(item, delimiter):
    return int(item.split(delimiter)[1])

def fillInGaps_realData(incDF,delimiter,stepSize):
    """Sometimes with real data a lot of the bins may be
    missing from the list of interesting reads. To have a
    projection matrix that has gaps and no weird re-indexing,
    this fills the missing Bin IDs"""
    binIDs = list(incDF.index)
    start = int(binIDs[0].split(delimiter)[1])
    end = int(binIDs[-1].split(delimiter)[1])
    allBins = ["Bin"+str(i) for i in range(start,end + 1,stepSize)]
    missingBins = list(set(allBins) - set(binIDs))
    missingDF = pd.DataFrame(0, index=missingBins, 
                             columns=incDF.columns)
    fullDF = pd.concat([incDF,missingDF])
    return fullDF

def makeHiC_fromInc(incDF):
    ## Convert incidence matrix to 2d hiC matrix
    nrow = incDF.shape[0]
    binIDs = list(incDF.index)
    if isinstance(binIDs[0],str) and ":" in binIDs[0]:
        sorted_binIDs = sorted(binIDs, key=lambda x: sort_key(x,':'))
    elif isinstance(binIDs[0],str) and ":" not in binIDs[0]:
        sorted_binIDs = sorted(binIDs, key=lambda x: sort_key(x,'Bin'))
    else:
        sorted_binIDs = binIDs.sort()
    df = pd.DataFrame(np.zeros(shape = (nrow,nrow)), index=sorted_binIDs, 
                      columns=sorted_binIDs)
    for read in incDF.columns:
        if isinstance(read,str):
            support = int(read.split(":")[1])
        else:
            support = 1
        arr = incDF[read][incDF[read] == 1].index
        for a in arr:
            df.loc[a][a] += support
        combs = list(combinations(arr,2))
        for c in combs:
            df.loc[c[0]][c[1]] += support
            df.loc[c[1]][c[0]] += support
    return(df)

def multiresolutionProjMat(df, ratios, report):
    """Given an incidence DF and information about how many (ratio)
    2-way interactions of a multiway interactions are present in the data,
    calculate projection matrices at each resolution"""
    projMat_dict = {}
    for cutoff in [0.5,1]:
        print(f"Making projection matrix for 2-way cutoff = {cutoff}")
        subset = [index for index, value in enumerate(ratios) if value >= cutoff]
        report.append(len(subset))
        subset_df = df.iloc[:,subset]
        hic_subset = makeHiC_fromInc(subset_df)
        projMat_dict[cutoff] = hic_subset
    return [projMat_dict,report]

def makeNorm_HiC_fromInc(incDF,weightList):
    ## 2D HiC matrix with counts normalized by weights
    nrow = incDF.shape[0]
    binIDs = list(incDF.index)
    if isinstance(binIDs[0],str) and ":" in binIDs[0]:
        sorted_binIDs = sorted(binIDs, key=lambda x: sort_key(x,':'))
    elif isinstance(binIDs[0],str) and ":" not in binIDs[0]:
        sorted_binIDs = sorted(binIDs, key=lambda x: sort_key(x,'Bin'))
    else:
        sorted_binIDs = binIDs.sort()
    df = pd.DataFrame(np.zeros(shape = (nrow,nrow)), index=sorted_binIDs, 
                      columns=sorted_binIDs)
    for read in incDF.columns:
        arr = incDF[read][incDF[read] == 1].index
        for a in arr:
            df.loc[a][a] += float(1/len(arr))
        combs = list(combinations(arr,2))
        for c in combs:
            df.loc[c[0]][c[1]] += weightList[read]
            df.loc[c[1]][c[0]] += weightList[read]
    return(df)
