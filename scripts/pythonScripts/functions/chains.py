import numpy as np
import pandas as pd
from itertools import combinations

def makeIncDF_fromChainDists(chainMat,cutoff):
    nrow = chainMat.shape[0]
    index_range = range(nrow)
    df = pd.DataFrame(index=index_range)
    columns_to_add = []
    chainMat_triu = np.triu(chainMat, k=1)
    counter = 0
    for i in range(nrow):
        condition1 = (0 < chainMat_triu[i, :])
        condition2 = (chainMat_triu[i, :] < cutoff)
        possNeighbors = list(np.where(condition1 & condition2)[0])
        if(possNeighbors):
            for ix in range(2,len(possNeighbors)):
                for comb in combinations(possNeighbors,ix):
                    d = chainMat_triu[np.ix_(comb,comb)]
                    if d.max() <= cutoff:
                        counter += 1
                        new_column = np.zeros(nrow)
                        new_column[list(comb)] = 1
                        columns_to_add.append(pd.Series(new_column, name=f'Read_{counter}'))

    df = pd.concat([df] + columns_to_add, axis=1)
    return(df)
