import numpy as np
import pandas as pd
from itertools import combinations
import multiprocessing
from statistics import median

from utils import flatten

def increaseIncDF_binSize(df,binSize):
    """This collapses consective bins in an incidence DF
    to reduce dimensions. This also means that some multiway
    interactions are now collapsed and those edges will be
    pruned"""
    result = []
    names = []
    for i in range(0,len(df) - binSize + 1,binSize):
        summed_value = df.loc[i:i+binSize-1,:].sum()
        summed_value[summed_value > 0] = 1
        names.append(f"Bin{i}:{i+binSize-1}")
        result.append(summed_value)
    result_df = pd.DataFrame(result,index = names)
    result_df = result_df.loc[:,result_df.sum() >= 2]
    return(result_df)

def dfToDict(df,result_dict):
    """Takes in an incidence DF and converts to 
    a dictionary of hyperedges"""
    for col in df.columns:
        twoWayRatio = float(col.split("_")[1])
        indices = df.index[df[col] == 1].tolist()
        key = '_'.join(indices)

        tmpList = result_dict.get(key, [])
        if tmpList:
            tmpList[0] += 1
            tmpList[1].append(twoWayRatio)
        else:
            tmpList = [1,[twoWayRatio]]
        result_dict[key] = tmpList
    return [result_dict]

def appendSingleBIncDict(chainDict,resultDict):
    """Takes in a binned inc dictionary which contains the
    read, number of times it occurred in a chain, and the
    % two-way support for that multiway read. It combines 
    each element, i.e., 3 elements per key to create a
    a bigger chunked dict."""
    for key, value in chainDict.items():
        tmpList = resultDict.get(key,[])
        if tmpList:
            tmpList[0] += value[0]
            tmpList[1] += 1
            tmpList[2].append(median(value[1]))
        else:
            tmpList = [value[0],1,[median(value[1])]]
        resultDict[key] = tmpList
    return(resultDict)

def combineChunkedBIncDicts(chainDict,resultDict):
    """Takes in chunked dicts which contain the
    read, number of chains wherein it occured, and the
    % two-way support for that multiway read across all chains. 
    It combines each element, i.e., 3 elements per key to create a
    a bigger chunked dict. 
    Also adding cardinality as the fourth element"""
    for key, value in chainDict.items():
        tmpList = resultDict.get(key,[])
        card = len(key.split("_"))
        if tmpList:
            tmpList[0] += value[0]
            tmpList[1] += value[1]
            tmpList[2].extend(value[2])
        else:
            tmpList = [value[0],value[1],value[2],card]
        resultDict[key] = tmpList
    return(resultDict)


def appendSingleDict(chainDict,resultDict):
    """Same as above when ratio is not being recorded"""
    for key, value in chainDict.items():
        tmpList = resultDict.get(key,[])
        if tmpList:
            tmpList += value
        else:
            tmpList = value
        resultDict[key] = tmpList
    return(resultDict)

def combineChunkedDicts(chainDict,resultDict):
    """Same as above without ratios"""
    for key, value in chainDict.items():
        tmpList = resultDict.get(key,[])
        card = len(key.split("_"))
        if tmpList:
            tmpList[0] += value
        else:
            tmpList = [value,card]
        resultDict[key] = tmpList
    return(resultDict)

class IncDFCreator:
    """Replacing this function: makeIncDF_fromChainDists"""
    def __init__(self, numProcesses, prim_cutoff, sec_cutoff, offDiagLim):
        self.numProcesses = numProcesses
        self.prim_cutoff = prim_cutoff
        self.sec_cutoff = sec_cutoff
        self.offDiagLim = offDiagLim

    def preprocessMat(self, chainMat):
        """Get upper triangular matrix and number of rows"""
        nrow = chainMat.shape[0]
        chainMat_triu = np.triu(chainMat, k=1)
        return chainMat_triu, nrow
    
    def assessMultiway(self, slice):
        """How many 2-way contacts fall below secondary cutoff"""
        total = np.count_nonzero(slice)
        if total > 0:
            passThresh = np.count_nonzero(slice[slice < self.sec_cutoff])
            ratio = passThresh / total
        else:
            ratio = 0
        return(ratio)
        
    def perRow(self, args):
        """Per row of upper tri matrix, get potential list of
        neighbors that fulfill distance criteria. Take iterative
        n choose k subsets of the matrix if the elements fall below 
        primary cutoff. If all elements fall 
        below the secondary distance cutoff needed to make it interact, then
        report that as a hyperedge"""
        chainMat_triu, row_ix = args
        columns_to_add = []
        ratioVec = []
        self.nrow = chainMat_triu.shape[0]
        vec = chainMat_triu[row_ix, row_ix + self.offDiagLim:]
        #print(vec)
        #condition1 = (0 < vec)
        condition2 = (vec < self.prim_cutoff)
        possNeighbors = [row_ix + self.offDiagLim + index for index in np.where(condition2)[0]]
        possNeighbors.insert(0, row_ix)
        #print(possNeighbors)
        if possNeighbors:
            for ix in range(2, len(possNeighbors)):
                for comb in combinations(possNeighbors, ix):
                    d = chainMat_triu[np.ix_(comb, comb)]
                    # if d.shape[0] == 0:
                    #     print(d.shape)
                    ratioUnderThresh = self.assessMultiway(d)
                    if ratioUnderThresh >= 0.5:
                        new_column = np.zeros(self.nrow)
                        new_column[list(comb)] = 1
                        columns_to_add.append(new_column)
                        ratioVec.append(ratioUnderThresh)
        if columns_to_add:
            return [columns_to_add,ratioVec]
        else:
            return None

    def mp(self, chainMat_triu, nrow):
        """Define multiprocessing pool"""
        pool = multiprocessing.Pool(self.numProcesses)
        argument_pairs = [(chainMat_triu, row_ix) for row_ix in range(nrow)]
        results = pool.map(self.perRow, argument_pairs)
        return results

    def makeIncDF_fromChainDists_mp(self, chainMat):
        """Run the dist matrix --> incidence DF in a parallelized
        fashion"""
        chainMat_triu, nrow = self.preprocessMat(chainMat)
        results = self.mp(chainMat_triu, nrow)
        filtered_results = [arr[0] for arr in results if arr is not None]
        ratioVec = flatten([vec[1] for vec in results if vec is not None])
        df = pd.DataFrame(np.concatenate(filtered_results)).T
        return [df,ratioVec]

    def makeIncDF_fromChainDists_single(self, chainMat):
        """Run the dist matrix --> incidence DF in a single-threaded
        fashion"""
        chainMat_triu, nrow = self.preprocessMat(chainMat)
        res = []
        rV = []
        for i in range(nrow):
            args = (chainMat_triu, i)
            row_result = self.perRow(args)
            if row_result is not None:
                res.append(row_result[0])
                rV.append(row_result[1])
        filtered_results = [arr for arr in res if arr is not None]
        df = pd.DataFrame(np.concatenate(filtered_results)).T
        ratioVec = flatten(rV)
        return [df, ratioVec]