def makeHiC_fromInc(incDF):
    ## Convert incidence matrix to 2d hiC matrix
    nrow = incDF.shape[0]
    ncol = incDF.shape[1]
    binIDs = list(incDF.index)
    df = pd.DataFrame(np.zeros(shape = (nrow,nrow)), index=binIDs, columns=binIDs)
    for read in incDF.columns:
        arr = incDF[read][incDF[read] == 1].index
        for a in arr:
            df.loc[a][a] += 1
        combs = list(combinations(arr,2))
        for c in combs:
            df.loc[c[0]][c[1]] += 1
            df.loc[c[1]][c[0]] += 1
    return(df)

def makeNorm_HiC_fromInc(incDF,weightList):
    ## 2D HiC matrix with counts normalized by weights
    nrow = incDF.shape[0]
    ncol = incDF.shape[1]
    binIDs = list(incDF.index)
    df = pd.DataFrame(np.zeros(shape = (nrow,nrow)), index=binIDs, columns=binIDs)
    for read in incDF.columns:
        arr = incDF[read][incDF[read] == 1].index
        for a in arr:
            df.loc[a][a] += float(1/len(arr))
        combs = list(combinations(arr,2))
        for c in combs:
            df.loc[c[0]][c[1]] += weightList[read]
            df.loc[c[1]][c[0]] += weightList[read]
    return(df)
