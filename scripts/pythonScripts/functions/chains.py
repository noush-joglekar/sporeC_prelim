import numpy as np
import pandas as pd
from itertools import combinations
import multiprocessing
import os

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
        indices = df.index[df[col] == 1].tolist()
        key = '_'.join(indices)

        result_dict[key] = result_dict.get(key, 0) + 1
    return(result_dict)

def dictToDF(hpDict):
    """Finally, takes in a dict and converts to incidence DF 
    for hypergraph construction. We can implement pruning based on
    heuristics later"""
    indices = list(set(flatten([key.split('_') for key in hpDict.keys()])))
    columns = []
    colnames = []
    counter = 0

    for key, value in hpDict.items():
        counter+=1
        col_ix = key.split('_')
        column = pd.Series([0] * len(indices),index = indices)  # Initialize row with zeros
        column[col_ix] = 1
        colName = f"Read{counter}:{value}"
        colnames.append(colName)
        columns.append(column)

    df = pd.concat(columns,axis=1)
    df.columns = colnames
    return(df)

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

    def perRow(self, args):
        """Per row of upper tri matrix, get potential list of
        neighbors that fulfill distance criteria. Take iterative
        n choose k subsets of the matrix if the elements fall below 
        primary cutoff. If all elements fall 
        below the secondary distance cutoff needed to make it interact, then
        report that as a hyperedge"""
        chainMat_triu, row_ix = args
        columns_to_add = []
        self.nrow = chainMat_triu.shape[0]
        vec = chainMat_triu[row_ix, row_ix + self.offDiagLim:]
        #condition1 = (0 < vec)
        condition2 = (vec < self.prim_cutoff)
        possNeighbors = [row_ix + self.offDiagLim + index for index in np.where(condition2)[0]]
        possNeighbors.insert(0, 0)
        if possNeighbors:
            for ix in range(2, len(possNeighbors)):
                for comb in combinations(possNeighbors, ix):
                    d = chainMat_triu[np.ix_(comb, comb)]
                    if d.max() <= self.sec_cutoff:
                        new_column = np.zeros(self.nrow)
                        new_column[list(comb)] = 1
                        columns_to_add.append(new_column)
        if columns_to_add:
            return columns_to_add
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
        filtered_results = [arr for arr in results if arr is not None]
        df = pd.DataFrame(np.concatenate(filtered_results)).T
        return df

    def makeIncDF_fromChainDists_single(self, chainMat):
        """Run the dist matrix --> incidence DF in a single-threaded
        fashion"""
        chainMat_triu, nrow = self.preprocessMat(chainMat)
        res = []
        for i in range(nrow):
            args = (chainMat_triu, i)
            row_result = self.perRow(args)
            if row_result is not None:
                res.append(row_result)
        filtered_results = [arr for arr in res if arr is not None]
        df = pd.DataFrame(np.concatenate(filtered_results)).T
        return df

class RealHiC:
    """Make a traditional HiC matrix from distance matrices"""
    def __init__(self, chain_dir, num_files):
        self.chain_dir = chain_dir
        self.num_files = num_files

    def distMatToBinary(self,ix):
        """If distance falls below a threshold, make entry 1 else 0"""
        file_path = f'{self.chain_dir}chain_dist_{ix}.txt'
        chain_mat = np.loadtxt(file_path)
        binary_dist_mat = np.where(chain_mat <= 500, 1, 0)
        return binary_dist_mat

    def distFilesToRealHiC(self, num_files):
        """Binarize a whole bunch of distance matrices
        specified by numFiles"""
        real_hic = None
        counter = 0
        for i in range(i, num_files):
            file_path = f'{self.chain_dir}chain_dist_{i}.txt'
            if os.path.isfile(file_path):
                counter += 1
                binary_dist_mat = self.distMatToBinary(i)
                if real_hic is None:
                    real_hic = binary_dist_mat
                else:
                    real_hic += binary_dist_mat
        if real_hic is not None:
            real_hic = real_hic.astype('float64')
            real_hic /= counter
        return real_hic
