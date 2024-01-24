import pandas as pd
import numpy as np
import random
import statistics
from itertools import combinations
from collections import defaultdict
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import wasserstein_distance
from scipy.stats import energy_distance
import matplotlib.pyplot as plt
import seaborn as sns

class cutoffEval:
    def __init__(self, keyCard, readSupport):
        self.keyCard = keyCard
        self.readSupport = readSupport

    def getCutoff(self, card):
        cRI = [i for i,x in enumerate(self.keyCard) if x == card]
        cRS = [self.readSupport[i] for i,x in enumerate(self.keyCard) if x == card]
        cardReadSupport = pd.Series(cRS)
        cutOff = cardReadSupport.describe()['75%']
        print(f"Cutoff = {cutOff} for card = {card}")
        cardIx = [cRI[i] for i in range(len(cRI)) if cRS[i] >= cutOff]
        return(cardIx)

    def runForAllCards(self):
        eL = []
        for card in list(set(self.keyCard)):
            aboveCutoffIx = self.getCutoff(card)
            eL.extend(aboveCutoffIx)
        return(eL)
    
class multiwayEval_realData:

    def __init__(self, keyCard, hpEdges, hpKeys, hpKeys_split, seed, 
                 toChoose,toPlotRef, toPlotInd, toPlotScatter,
                 quartile, plotDir,outDir):
        self.keyCard = keyCard
        self.hpEdges = hpEdges
        self.hpKeys = hpKeys
        self.hpKeys_split = hpKeys_split
        self.seed = seed
        self.toPlotRef = toPlotRef
        self.toPlotInd = toPlotInd
        self.toPlotScatter = toPlotScatter
        self.toChoose = toChoose
        self.qt = quartile
        self.plotDir = plotDir
        self.outDir = outDir

    def getReadSupportStatsPerCard(self, ixList, card):
        """Get disribution of read support per card so as to impose cutoffs
        in the end as to the trustworthiness of reads"""
        readSupps = [self.hpEdges[self.hpKeys[ix]] for ix in ixList]
        rsStats = pd.Series(readSupps).describe()
        rsStats.to_pickle(f'{self.outDir}readStats_card{card}.pkl')
        return

    def makeAllReferenceHashDicts(self):
        """For all available cards, create look up table of probability of
        seeing all possible subsets as a function of mean distances
        between subsets"""
        probHash = defaultdict(dict)

        for card in range(max(self.keyCard),2,-1):
            print("Calculating for card=",card)
            ixList = [index for index,element in enumerate(self.keyCard) if element == card]
            print("There are ",len(ixList),"reads")
            self.getReadSupportStatsPerCard(ixList,card)
            for n in range(2,card):
                print(f"Creating a probability hash for {n}-way subsets")
                hashID = f'{card}sub{n}'
                probHash[hashID] = self.getNWayProbsPerCard(card,n)
        return(probHash)
    
    def getNWayProbsPerCard(self,card,n):
        """Get distance-stratified expected probability distribution for high
        cardinality reads and their subsets of cardinality n"""
        ixList = [index for index,element in enumerate(self.keyCard) if element == card]
        mainDict = defaultdict(int)

        for ix in ixList:
            readSupport, nWayDict = self.makeNWayDict(ix,n)
            for key in nWayDict.keys():
                mainDict[key] += nWayDict[key]
        normalized_values = self.getProbabilitiesByDist(mainDict)

        if normalized_values is None: ##### TEST
            return 
        else:
            if self.toPlotRef is True:
                self.plotReadFreqsPerCard(mainDict,normalized_values,card,n)
            
            probabilityDict = {list(mainDict.keys())[ix]:normalized_values[ix] for ix in range(len(normalized_values))}
            return(probabilityDict)
    
    def makeNWayDict(self,ix,n): ### Breaks down when filtered
        """For a high cardinality read, make a dictionary
        that outputs the frequency and mean distance for subsets
        of cardinality n"""
        readSupport = self.hpEdges[self.hpKeys[ix]]
        splitKey = self.hpKeys_split[ix]
        combs = list(combinations(splitKey,n))
        subsetDict = defaultdict(int)
        for comb in combs:
            subsetEdge = '_'.join(map(str, comb))
            try:
                subsetEdgeReads = self.hpEdges[subsetEdge]
            except KeyError:
                subsetEdgeReads = 0
            meanDist = self.getNWayMeanDistPerSubset(comb)
            subsetDict[meanDist] = subsetEdgeReads
        return(readSupport, subsetDict)
    
    def getNWayMeanDistPerSubset(self,comb):
        """Get the mean pairwise distance given a
        high cardinality read"""
        twoWays = list(combinations(comb,2))
        twoWayDist = [self.getPairwiseDist(c) for c in twoWays]
        # meanDist = round(statistics.mean(twoWayDist))
        meanDist = round(statistics.geometric_mean(twoWayDist)) ## geometric mean
        return(meanDist)
    
    def getPairwiseDist(self,combination):
        """Given a two-way interaction, find dist between them"""
        ends = [int(combination[i].split("Bin")[1]) for i in [0,1]]
        diff = ends[1] - ends[0]
        return(diff)
    
    def getProbabilitiesByDist(self,mainDict):
        """Get probabilities of reads occurring by mean dist to
        obtain the expected probability distribution from the data"""
        totalReadsForCard = sum(mainDict.values())
        if totalReadsForCard == 0:
            return None
        else:
            normalized_values = [value / totalReadsForCard for value in mainDict.values()]
            ## return total reads as well
            return(normalized_values)
    

    def plotReadFreqsPerCard(self,mainDict,normalized_values,card,n):
        """Takes in a generated dictionary and plots the
        read frequency as well as probabilities (first 20) by the
        distance"""

        # Step 1: Plot a barplot with dictionary keys on the x-axis and values on the y-axis
        sns.barplot(x=list(mainDict.keys()), y=list(mainDict.values()),color = "grey")
        plt.title(f"Read frequency of {n}-way subsets for Card={card}")
        plt.xlabel(f"Mean (Geom) distance between {n}-way interactions")
        plt.ylabel("Frequency")
        plt.xticks(rotation=90,fontsize=5)
        plt.savefig(f'{self.plotDir}DistStratProbs_Card{card}sub{n}.png',bbox_inches = 'tight',facecolor = "white")
        #plt.show()

        # Step 2: Plot the same barplot with the ratio of each value to the total
        sns.barplot(x=list(mainDict.keys()), y=normalized_values,color = "grey")
        plt.title(f"Probability of occurrence {n}-way distance for Card={card}")
        plt.xlabel(f"Mean (Geom) distance between {n}-way interactions")
        plt.ylabel("Probabilities")
        plt.xlim(0,19)
        plt.xticks(rotation=90)
        plt.savefig(f'{self.plotDir}DistStratProbs_top20_Card{card}sub{n}.png',bbox_inches = 'tight',facecolor = "white")
        #plt.show()
        plt.close()
        return None
    
    def statsForAllReads(self,probHash):
        """Wrapper for all reads"""
        for card in range(max(self.keyCard),3,-1):
            allStats = self.statsForAllCardSubsets(card,probHash)
            if allStats is not None:
                self.processConcatDFs(allStats[0],"wDist",False,card)
                self.processConcatDFs(allStats[1],"cosineSim",True,card)
                self.processConcatDFs(allStats[2],"eDist",False,card)
                self.processConcatDFs(allStats[3],"empDist",True,card)
        for card in [3]:
            allStats = self.statsForAllCardSubsets(card,probHash)
            self.processConcatDFs(allStats[0],"wDist",False,card)
            self.processConcatDFs(allStats[1],"cosineSim",True,card)
            self.processConcatDFs(allStats[2],"eDist",False,card)
            self.processConcatDFs(allStats[3],"empDist",True,card)            
        return None

    def processConcatDFs(self,metricDF,metricName,simiBool,card):
        """summarize and write output for each metric"""
        if simiBool is True:
            qtCutoff = self.qt
            operator = "leq"
        else:
            qtCutoff = 100 - self.qt
            operator = "geq"
        if card == 3:
            metricCutoff = pd.Series(metricDF['3Sub2']).describe()[f'{qtCutoff}%']
            self.plotSimilarityHist(metricDF['3Sub2'],metricName,card)
            if operator == "leq":
                mStatus = [1 if x else 0 for x in (metricDF['3Sub2'] <= metricCutoff)]
            else:
                mStatus = [1 if x else 0 for x in (metricDF['3Sub2'] >= metricCutoff)]
        else:
            summaryMetric = metricDF.filter(like="Sub").agg((np.mean,np.std),axis = 1)        
            metricCutoff = self.getCutoff(summaryMetric,qtCutoff)
            self.plotSimilarityHist(summaryMetric['mean'],metricName,card)
            if operator == "leq":
                mStatus = [1 if x else 0 for x in (summaryMetric['mean'] <= metricCutoff)]
            else:
                mStatus = [1 if x else 0 for x in (summaryMetric['mean'] >= metricCutoff)]
            if self.toPlotScatter is True:
                self.plotScatterWithErrorBars(summaryMetric,metricName,card)
        metricDF['Status'] = mStatus
        print(f"Writing output for {metricName}")
        metricDF.to_csv(f'{self.outDir}/{metricName}_card{card}.csv',sep = "\t",index=False)
        return None
    
    def getCutoff(self,summaryDF,quartile):
        q = f'{quartile}%'
        cutoff = pd.Series(summaryDF['mean']).describe()[q]
        return(cutoff)

    def statsForAllCardSubsets(self,card,probHash):
        """Wrapper for the stats per card. Calculates for all subsets of a card and binds into a df"""
        print("Calculating for card=",card)
        ixList = [index for index,element in enumerate(self.keyCard) if element == card]
        print("There are ",len(ixList),"reads")
        if len(ixList) > 1:
            cardToChoose = min(self.toChoose,len(ixList))
            print("Calculating stats for",cardToChoose,"reads")
            random.seed(self.seed)
            revised_ixes = random.sample(ixList,cardToChoose)
            C1 = []
            C2 = []
            C3 = []
            C4 = []
            expHashNames = []
            for n in range(2,card):
                stats = self.getStatsPerCard(card,n,probHash,revised_ixes)
                if stats is not None:
                    expHashNames.append(str(card)+"Sub"+str(n))
                    C1.append(stats[0])
                    C2.append(stats[1])
                    C3.append(stats[2])
                    C4.append(stats[3])
                if n == 2: ## None clause?
                    C5 = stats[4]
                    filteredIxes = stats[5]
            #cN = [str(card)+"Sub"+str(i) for i in range(2,card)]
            cN = expHashNames
            df1 = pd.DataFrame(C1).T
            df2 = pd.DataFrame(C2).T
            df3 = pd.DataFrame(C3).T
            df4 = pd.DataFrame(C4).T
            df1.columns = cN
            df1['Edge_ix'] = filteredIxes
            df1['ReadSupport'] = C5
            df2.columns = cN
            df2['Edge_ix'] = filteredIxes
            df2['ReadSupport'] = C5
            df3.columns = cN
            df3['Edge_ix'] = filteredIxes
            df3['ReadSupport'] = C5
            df4.columns = cN
            df4['Edge_ix'] = filteredIxes
            df4['ReadSupport'] = C5
            return(df1,df2,df3,df4)

    def getStatsPerCard(self,card,n,probHash,revised_ixes):
        """Per cardinality, get a subset of reads (for computational
        purposes), calculate the similarity of observed n-way interactions
        to the expected value, and output a list"""
        wdistList = []
        cosList = []
        edistList = []
        empDistList = []
        readSuppList = []
        filteredIx = []
        expHash = f'{card}sub{n}'
        for ix in revised_ixes:
            if expHash in probHash:
                stats = self.getReadExpectednessStats(card,ix,n,probHash[expHash])
                if stats is not None:
                    wdistList.append(stats[0])
                    cosList.append(stats[1])
                    edistList.append(stats[2])
                    empDistList.append(stats[3])
                    readSuppList.append(stats[4])
                    filteredIx.append(ix)
        return(wdistList, cosList, edistList, empDistList, readSuppList, filteredIx)
    
    def getReadExpectednessStats(self,card,ix,n,hash):
        """Per read, get the observed distribution of n-way contacts
        and calculate similarity to the expected distribution"""
        readSupport, readDict = self.makeNWayDict(ix,n)
        readPercs= self.getProbabilitiesByDist(readDict)
        if readPercs is None:
            return None
        else:
            probVals = [hash[k] for k in readDict.keys()]
            normProbs = [p / sum(probVals) for p in probVals]
            probVals = normProbs  ## Normalizing mostly for the empirical calculation
            if self.toPlotInd is True:
                if len(readPercs) > 2 and len(probVals) > 2:
                    self.makeSanityCheckPlotsPerRead(readDict,readPercs,probVals,card,n,ix)
            wDist = wasserstein_distance(readPercs, probVals)
            similarity = cosine_similarity(np.array(readPercs).reshape(1,-1), 
                                        np.array(probVals).reshape(1,-1))[0,0]
            eDist = energy_distance(readPercs, probVals)
            empDist = self.calcEmpDist(readPercs, probVals)
            return(wDist, similarity, eDist, empDist, readSupport)
    
    def calcEmpDist(self,readPercs, probVals):
        """Trying a different metric for distance. 
        We will figure out the cutoff later"""
        obs = np.log(readPercs)
        exp = np.log(probVals)
        if 0 in exp:
            meanOE = 1
        else:
            obsOverExp = obs/exp
            meanOE = np.mean(obsOverExp[np.isfinite(obsOverExp)])
        return(meanOE)

    def makeSanityCheckPlotsPerRead(self,readDict, readPercs, probVals, card, n, ix):
        """For specific reads, plot the observed versus expected distributions
        of two-way interactions along with a fitted spline. Additionally outputs
        the slope (useful only if line) and wDist values"""
        Distances = list(readDict.keys())
        print(card,"sub",n)
        print(ix)
        #print(readPercs)
        #print(probVals)
        
        # Create a 2x2 grid of subplots
        fig, axs = plt.subplots(2, 2, figsize=(6, 6))
        # Plot 1: Observed w/ spline
        lowess1 = sns.regplot(x=Distances, y=readPercs, 
                              lowess=True, ci=None, color='red', ax=axs[0, 0])
        axs[0, 0].set_title(f"Observed w/ spline Card{card}sub{n}")
        # Plot 2: Expected w/ spline
        lowess2 = sns.regplot(x=Distances, y=probVals, 
                              lowess=True, ci=None, color='grey', ax=axs[0, 1])
        axs[0, 1].set_title(f"Expected w/ spline Card{card}sub{n}")

        # Calculate smoothed values and slopes
        y_smoothed1 = lowess1.get_lines()[0].get_ydata()
        y_smoothed2 = lowess2.get_lines()[0].get_ydata()
        slope1 = np.gradient(y_smoothed1)
        slope2 = np.gradient(y_smoothed2)

        # Plot 3: Barplot for readPercs
        sns.barplot(x=list(readDict.keys()), y=readPercs, ax=axs[1, 0])
        axs[1, 0].set_title(f"Barplot for readPercs Card{card}sub{n}")
        # Plot 4: Barplot for probVals
        sns.barplot(x=list(readDict.keys()), y=probVals, ax=axs[1, 1])
        axs[1, 1].set_title(f"Barplot for normalized prior probs Card{card}sub{n}")

        #correlation = np.corrcoef(readPercs, probVals)[0, 1]
        wDist = wasserstein_distance(readPercs, probVals)
        similarity = cosine_similarity(np.array(readPercs).reshape(1,-1), 
                    np.array(probVals).reshape(1,-1))[0,0]
        eDist = energy_distance(readPercs, probVals)
        empDist = self.calcEmpDist(readPercs, probVals)
        
        # Print the slopes
        print("Comparing -------------")
        print(f"Slope for observed: {slope1.mean()}")
        print(f"Slope for expected: {slope2.mean()}")
        print(f"Wasserstein distance: {wDist}")
        print(f"Energy distance: {eDist}")
        print(f"Empirical distance: {empDist}")
        print(f"Cosine similarity: {similarity}")

        # Adjust layout and show the subplots
        plt.tight_layout()
        plt.savefig(f'{self.plotDir}/SingleReadPlot_Card{card}sub{n}_{ix}.png',bbox_inches = 'tight',facecolor = "white")
        plt.show()
        plt.close()
        return None

    def plotSimilarityHist(self,metricDistr,metricName,card):
        """Given distribution for a cardinality, plot the wass dist
        and cosine similarity values so that we can settle on a heuristic"""
        plt.figure(figsize=(4,4))
        plt.hist(metricDistr,color='blue', alpha=0.7, width=0.3,bins = 201)
        plt.xlim(-1,1.01)
        plt.title(f"Mean {metricName} distance for card={card}")
        plt.close()
        return None

    def plotScatterWithErrorBars(self,summaryDF,metric,card):
        """Given summary stats for a certain number of reads, plot scatter plot with error bars 
        representing change over all subsets"""
        sortedSumm = summaryDF.sort_values(by='mean')
        x = [i+1 for i in range(summaryDF.shape[0])] # Use the index as x-axis
        y = sortedSumm['mean']
        error = sortedSumm['std']

        plt.scatter(x, y, label=f'mean {metric}', marker='o')
        plt.errorbar(x, y, yerr=error, linestyle='None', color='grey', capsize=3)
        plt.xlabel('ReadID')
        plt.ylabel(f'mean {metric}')
        plt.ylim(0,2)
        plt.title(f'Distribution for card={card}')
        plt.legend()
        plt.savefig(f'{self.plotDir}/ScatterPlot_{metric}_Card{card}_max{self.toChoose}reads.png',
                    bbox_inches = 'tight',facecolor = "white")
        #plt.show()
        plt.close()
        return None
