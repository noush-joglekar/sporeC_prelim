import argparse
import pandas as pd
import numpy as np
import sys

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
from matplotlib import colors

sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')
from chains import dictToDF, dfToDict
from edgeWeightFormulations import finalBounded, finalBounded_fromEdge

def getDifference(f1,f2):
    npm1 = normalizeMats(f1)
    npm2 = normalizeMats(f2)
    diffMat = npm1.sub(npm2, fill_value=0)
    numeric_values = diffMat.index.to_series().str.extract(r'Bin(\d+)', 
                                                           expand=False).astype(int)
    diffMat_sorted = diffMat.loc[numeric_values.sort_values().index, 
                                 numeric_values.sort_values().index]
    return(diffMat_sorted)

def normalizeMats(fPath):
    pm = pd.read_pickle(fPath)
    np.fill_diagonal(np.asmatrix(pm),0)
    npm = pm / np.nanmax(pm)
    return(npm)

def plotDiffMat(args,diffMat):
    plt.figure(figsize=(6, 4))
    im = plt.imshow(diffMat, cmap="BrBG", norm=colors.SymLogNorm(linthresh=0.02, linscale=0.02,
                                                vmin=-0.5, vmax=0.5, base=10))
    plt.colorbar(im, fraction=0.046, pad=0.04, label='balanced');
    plt.title(f"Projection matrix, Card={args.card}, {args.chrom}: {args.ct1} - {args.ct2}")
    plt.savefig()
    return()

def getEnrichedReads(diffMat,d1,d2):
    npm1_enr = diffMat.where(diffMat >= 0.05, 0)
    npm2_enr =  diffMat.where(diffMat <= -0.05, 0)

    positive_rows, positive_columns = (npm1_enr > 0).any(axis=1), (npm1_enr > 0).any()
    negative_rows, negative_columns = (npm2_enr < 0).any(axis=1), (npm2_enr < 0).any()

    id1 = pd.read_pickle(d1)
    id2 = pd.read_pickle(d2)

    multiwayReads_pm1 = id1.loc[positive_rows].sum() > 2
    enriched_pm1_incDF = id1.loc[positive_rows,multiwayReads_pm1]
    print(enriched_pm1_incDF.shape)

    multiwayReads_pm2 = id2.loc[negative_rows].sum() > 2
    enriched_pm2_incDF = id2.loc[negative_rows,multiwayReads_pm2]
    enriched_pm2_incDF.shape
    return()


def parse_args():
    parser = argparse.ArgumentParser(
        description="""""")
    parser.add_argument("dataDir",type=str, help="Main data dir")
    parser.add_argument("matDir",type=str, help="location of pkl files relative to data dir")
    parser.add_argument("chrom",type=str, help="chr ID")
    parser.add_argument("card",type=str, help="read card")
    parser.add_argument("ct1",type=str, help="cell line 1")
    parser.add_argument("ct2",type=str, help="cell line 2")
    return parser.parse_args()
    
def main():
    args = parse_args()
    dfDir = f'{args.dataDir}{args.runDir}dfs_{args.chrom}/'
    outDir = f'{args.dataDir}{args.outputDir}'
    f1 = f'{args.dataDir}{args.matDir}projMat_Int_{args.ct1}_Card{args.card}_{args.chrom}.pkl'
    f2 = f'{args.dataDir}{args.matDir}projMat_Int_{args.ct2}_Card{args.card}_{args.chrom}.pkl'
    d1 = f'{args.dataDir}{args.matDir}IncDF_IntReads_{args.ct1}_Card{args.card}_{args.chromID}.pkl'
    d2 = f'{args.dataDir}{args.matDir}IncDF_IntReads_{args.ct2}_Card{args.card}_{args.chromID}.pkl'





if __name__ == "__main__":
    main()