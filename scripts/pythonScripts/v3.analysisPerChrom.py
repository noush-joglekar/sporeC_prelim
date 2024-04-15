import argparse
import pickle
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
import os

import sys
sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')
from chains import dictToDF, dfToDict
from incidenceToProjection import makeHiC_fromInc, fillInGaps_realData
from edgeWeightFormulations import finalBounded
from promethData_multiwayExpectedProbs import cutoffEval

def main():
    args = parse_args()
    dfDir = f'{args.dataDir}{args.runDir}dfs_{args.chrom}/'
    outDir = f'{args.dataDir}{args.outputDir}'
    matDir = f'{outDir}/matrices/'
    os.makedirs(matDir,exist_ok=True)
    cardList = [int(c) for c in args.card.split(",")]
    dfList = []
    okayCards = []
    for card in cardList:
        card_fullDF = perCard(args,dfDir,
                            outDir,matDir,card)
        if card_fullDF is not None:
            dfList.append(card_fullDF)
            okayCards.append(card)
    if args.combined is True:
        if(len(okayCards) > 1):
            updatedCards =  '-'.join(map(str,okayCards))
            print(f"Combining the interesting reads from {updatedCards}")
            collatedIncDF = pd.concat(dfList,axis = 1)
            print(collatedIncDF.shape)
            plotInterestingProjMat(args,collatedIncDF,outDir,matDir,updatedCards)
        else:
            print("Nothing to combine!")
    return


def perCard(args,dfDir,outDir,matDir,card):
    coSimFile = pd.read_csv(f'{dfDir}cosineSim_card{card}.csv',sep = "\t")
    empDistFile = pd.read_csv(f'{dfDir}empDist_card{card}.csv',sep = "\t")
    pklFile = f'{args.dataDir}{args.runDir}hyperEdges_{args.identifier}_{args.chrom}.pkl'

    if empDistFile.shape[0] > 500:
        print(f"There are a total of {empDistFile.shape[0]} card{card} reads for {args.identifier} celline for {args.chrom}. Proceeding ...")
        hpKeys, updatedDict, hpEdges = extractInterestingEdges(args,pklFile)

        print("Determining which reads are interesting ...")
        empDistStatus = calcDistMetricStatus(empDistFile,args.empCutoff)
        agreement_status = getAgreementStatus(coSimFile,empDistStatus)
        consensusIx = [coSimFile['Edge_ix'][ix] for ix,x in enumerate(agreement_status) 
                       if x == "Agree:Interesting"]
        print(f'There are {len(consensusIx)} interesting reads by cosine similarity and empirical distance calculation')
        intScores, randScores, subset_incDF = getEdgeScores(args, consensusIx,
                            coSimFile,hpKeys,updatedDict,matDir,card)
        print('Generating and plotting all the output now ...')
        print("Distribution of edge scores")
        plotEdgeScores(args,outDir,intScores,randScores,card)
        print("Projection matrix for interesting reads")
        fullDF = plotInterestingProjMat(args,subset_incDF,outDir,matDir,card)
        print("Projection matrix for the entire dataset")
        plotOG_fullDataset(args,updatedDict,outDir)
        print(f"Skipping the projection matrix for all reads for card = {card}")
        #plotAll_perCard(args,hpEdges,outDir,card)
        return(fullDF)
    else:
        print(f"Less than 500 card{card} reads for {args.identifier} celline for {args.chrom}.")


def extractInterestingEdges(args,pklFile):
    with open(pklFile,'rb') as f:
        hpEdges = pickle.load(f)
    print(f"Processing all the hyperedges for {args.identifier} {args.chrom} from pickle file")
    hpKeys = [k for k in hpEdges.keys()]
    keyCard = [len(k.split("_")) for k in hpKeys]
    print("A total of",len(hpKeys),"initial interactions")
    
    readSupport = [v for v in hpEdges.values()]
    cE = cutoffEval(keyCard,readSupport)
    passedReadIx = cE.runForAllCards()
    # atLeastTwoChains = [i for i,x in enumerate(readSupport) if x >=2]
    updatedDict = {hpKeys[i]:readSupport[i] for i in passedReadIx} # atLeastTwoChains
    hpKeys = [k for k in updatedDict.keys()]
    print("A total of",len(hpKeys),"final interactions")
    return([hpKeys,updatedDict,hpEdges])

def getEdgeScores(args,consensusIx,coSimFile,hpKeys,updatedDict,matDir,card):
    randList = random.sample(range(coSimFile.shape[0]),len(consensusIx))
    subsetKeys = [hpKeys[ix] for ix in consensusIx]
    subsetDict = {key: updatedDict[key] for key in subsetKeys}
    subset_incDF = dictToDF(subsetDict)
    pd.to_pickle(subset_incDF,f'{matDir}IncDF_IntReads_{args.identifier}_Card{card}_{args.chrom}.pkl')
    finalBoundedScores = [finalBounded(list(subset_incDF[c])) for c in subset_incDF.columns]

    rsubsetKeys = [hpKeys[ix] for ix in randList]
    rsubsetDict = {key: updatedDict[key] for key in rsubsetKeys}
    rsubset_incDF = dictToDF(rsubsetDict)

    randomBoundedScores = [finalBounded(list(rsubset_incDF[c])) for c in rsubset_incDF.columns]
    return([finalBoundedScores,randomBoundedScores,subset_incDF])

def plotInterestingProjMat(args,subset_incDF,outDir,matDir,card):
    fullDF = fillInGaps_realData(subset_incDF,"Bin",1)
    projMat = makeHiC_fromInc(fullDF)
    pd.to_pickle(projMat,f'{matDir}projMat_Int_{args.identifier}_Card{card}_{args.chrom}.pkl')
    mv = np.nanmax(projMat)
    if mv <= 100:
        vm = 10
    else:
        vm = round(mv / 100) * 100
    plt.figure(figsize=(6, 4))
    im = plt.imshow(projMat, cmap="YlOrRd",norm = LogNorm(vmax = vm))
    plt.colorbar(im, fraction=0.046, pad=0.04, label='balanced');
    plt.title(f"Card={card}, {args.chrom}: {args.identifier}")
    plt.savefig(f'{outDir}ProjMat_Int_{args.identifier}_Card{card}_{args.chrom}.png',
                bbox_inches = 'tight',facecolor = "white")
    plt.close()
    return(fullDF)

def plotOG_fullDataset(args,updatedDict,outDir):
    plotFile = f'{outDir}ProjMat_AllReads_{args.identifier}_{args.chrom}.png'
    if not os.path.exists(plotFile):
        print("File does not exist. Generating matrix")
        ogDF = dictToDF(updatedDict)
        ogProjMat = makeHiC_fromInc(ogDF)
        mv = np.nanmax(ogProjMat)
        if mv <= 2000:
            vm = 250
        else:
            vm = round(mv / 1000) * 1000
        print("Plotting ....")
        plt.figure(figsize=(6, 4))
        im = plt.imshow(ogProjMat, cmap="YlOrRd",norm = LogNorm(vmax = vm))
        plt.colorbar(im, fraction=0.046, pad=0.04, label='balanced');
        plt.title(f"All Reads (>=2), {args.chrom}: {args.identifier}")
        plt.savefig(plotFile, bbox_inches = 'tight',facecolor = "white")
        plt.close()
        return
    else:
        print("File already exists, moving on.")

def plotAll_perCard(args,hpEdges,outDir,card):
    print("Generating matrix ...")
    nWay = [k for k in hpEdges.keys() if len(k.split("_")) == 3]
    nWayDict = {k:v for k,v in hpEdges.items() if k in nWay and v >= 2}
    nWayIncDF = dictToDF(nWayDict)
    nWayProjMat = makeHiC_fromInc(nWayIncDF)
    mv = np.nanmax(nWayProjMat)
    if mv <= 2000:
        vm = 250
    else:
        vm = round(mv / 1000) * 1000
    print("Plotting ....")
    plt.figure(figsize=(6, 4))
    im = plt.imshow(nWayProjMat, cmap="YlOrRd",norm = LogNorm(vmax = vm))
    plt.colorbar(im, fraction=0.046, pad=0.04, label='balanced');
    plt.title(f"All Card={card}, {args.chrom}: {args.identifier}")
    plt.savefig(f'{outDir}ProjMat_AllCard{card}_{args.identifier}_{args.chrom}.png',
                bbox_inches = 'tight',facecolor = "white")
    plt.close()
    return

def plotEdgeScores(args,outDir,intScores,randScores,card):
    plt.hist(intScores,bins = 21,alpha =0.5)
    plt.hist(randScores,bins = 21,alpha =0.5)
    plt.title(f'{args.identifier} {args.chrom}:Card{card}')
    plt.legend(title="Score",loc="upper right",
               labels = ("Interesting reads","Random reads"))
    plt.savefig(f'{outDir}EdgeScores_Card{card}_{args.identifier}_{args.chrom}.png',
                bbox_inches = 'tight',facecolor = "white")
    plt.close()
    return

def getAgreementStatus(coSimFile,empDistStatus):
    agreement_status = ["Agree:Interesting" if v1 == 1 and v2 == 1 else 
                    "CoSim only" if v1 == 1 else "empDist only" if v2 == 1 
                    else "Expected" for 
                    v1, v2 in zip(coSimFile['Status'], empDistStatus)] #empDistStatus CHANGE
    return(agreement_status)

def getCutoff(summaryDF,quartile):
    q = f'{quartile}%'
    cutoff = pd.Series(summaryDF['mean']).describe()[q]
    return(cutoff)

def calcDistMetricStatus(inFile,empCutoff):
    """Repeat of what is in original file. I kept messing up for emp dist.
    But this is good incase we decide to change the 4th quartile stipulation
    Will just have to change to 25 and less-equal for cosine"""
    summary_dist = inFile.filter(like="Sub").apply(
        lambda row: [np.mean(row), np.std(row)], axis=1, result_type='expand')
    summary_dist.columns = ['mean','sd']
    distCutoff = getCutoff(summary_dist,empCutoff)
    distStatus = [1 if x else 0 for x in (summary_dist['mean'] >= distCutoff)]
    return(distStatus)


def parse_args():
    parser = argparse.ArgumentParser(
        description="""""")

    parser.add_argument("dataDir",type=str, help="Main data dir")
    parser.add_argument("runDir",type=str, help="Working directory relative to data dir")
    parser.add_argument("outputDir",type=str, help="Output directory where plots will be placed")
    parser.add_argument("chrom",type=str, help="PoreC fragment file relative to working dir")
    parser.add_argument("identifier",type=str, help="PoreC fragment file relative to working dir")
    parser.add_argument("empCutoff", type =int, help = "cutoff for int reads from empDist distr", default=75)
    parser.add_argument("card", type=str, help="Comma separated list of cardinalities to consider")
    parser.add_argument("--combined", action="store_true", 
                        help="Specify if you want a combined proj mat of interesting reads from all cards")
    ## CARD. Combined?
    
    return parser.parse_args()

if __name__ == "__main__":
    main()
