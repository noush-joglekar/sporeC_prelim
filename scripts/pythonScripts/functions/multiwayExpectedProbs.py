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

class multiwayEval:

    def __init__(self, keyCard, hpEdges, hpKeys_split, seed, 
                 toChoose,toPlotRef, toPlotInd, toPlotScatter,
                 quartile, plotDir,outDir):
        self.keyCard = keyCard
        self.hpEdges = hpEdges
        self.hpKeys_split = hpKeys_split
        self.seed = seed
        self.toPlotRef = toPlotRef
        self.toPlotInd = toPlotInd
        self.toPlotScatter = toPlotScatter
        self.toChoose = toChoose
        self.qt = quartile
        self.plotDir = plotDir
        self.outDir = outDir

    def makeAllReferenceHashDicts(self):
        """For all available cards, create look up table of probability of
        seeing all possible subsets as a function of mean distances
        between subsets"""
        probHash = defaultdict(dict)

        for card in range(max(self.keyCard),2,-1):
            print("Calculating for card=",card)
            ixList = [index for index,element in enumerate(self.keyCard) if element == card]
            print("There are ",len(ixList),"reads")
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
            nWayDict = self.makeNWayDict(ix,n)
            for key in nWayDict.keys():
                mainDict[key] += nWayDict[key]

        normalized_values = self.getProbabilitiesByDist(mainDict)
        if self.toPlotRef is True:
            self.plotReadFreqsPerCard(mainDict,normalized_values,card,n)
        
        probabilityDict = {list(mainDict.keys())[ix]:normalized_values[ix] for ix in range(len(normalized_values))}
        return(probabilityDict)
    
    def makeNWayDict(self,ix,n):
        """For a high cardinality read, make a dictionary
        that outputs the frequency and mean distance for subsets
        of cardinality n"""
        splitKey = self.hpKeys_split[ix]
        combs = list(combinations(splitKey,n))
        subsetDict = defaultdict(int)
        for comb in combs:
            subsetEdge = '_'.join(map(str, comb))
            subsetEdgeReads = self.hpEdges[subsetEdge]
            meanDist = self.getNWayMeanDistPerSubset(comb)
            subsetDict[meanDist] = subsetEdgeReads
        return(subsetDict)
    
    def getNWayMeanDistPerSubset(self,comb):
        """Get the mean pairwise distance given a
        high cardinality read"""
        twoWays = list(combinations(comb,2))
        twoWayDist = [self.getPairwiseDist(c) for c in twoWays]
        meanDist = round(statistics.mean(twoWayDist))
        return(meanDist)
    
    def getPairwiseDist(self,combination):
        """Given a two-way interaction, find dist between them"""
        ends = [int(combination[i].split(":")[1]) for i in [0,1]]
        diff = (ends[1] - ends[0]) / 5
        return(diff)
    
    def getProbabilitiesByDist(self,mainDict):
        """Get probabilities of reads occurring by mean dist to
        obtain the expected probability distribution from the data"""
        totalReadsForCard = sum(mainDict.values())
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
        plt.xlabel(f"Mean distance between {n}-way interactions")
        plt.ylabel("Frequency")
        plt.xticks(rotation=90,fontsize=5)
        plt.savefig(f'{self.plotDir}DistStratProbs_Card{card}sub{n}.png',bbox_inches = 'tight',facecolor = "white")
        #plt.show()

        # Step 2: Plot the same barplot with the ratio of each value to the total
        sns.barplot(x=list(mainDict.keys()), y=normalized_values,color = "grey")
        plt.title(f"Probability of occurrence {n}-way distance for Card={card}")
        plt.xlabel(f"Mean distance between {n}-way interactions")
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
            wDistStats, cosStats = self.statsForAllCardSubsets(card,probHash)
            summaryCos = cosStats.filter(like="Sub").agg((np.mean,np.std),axis = 1)
            summaryWDist = wDistStats.filter(like="Sub").agg((np.mean,np.std),axis = 1)
            cosCutoff = self.getCutoff(summaryCos,self.qt)
            cStatus = [1 if x else 0 for x in (summaryCos['mean'] <= cosCutoff)]
            cosStats['Status'] = cStatus
            wdistCutoff = self.getCutoff(summaryWDist,100 - self.qt)
            wStatus = [1 if x else 0 for x in (summaryWDist['mean'] >= wdistCutoff)]
            wDistStats['Status'] = wStatus
            if self.toPlotScatter is True:
                self.plotScatterWithErrorBars(summaryCos,"CosineSim",card)
                self.plotScatterWithErrorBars(summaryWDist,"WassDist",card)
            self.plotSimilarityHist(summaryWDist['mean'],summaryCos['mean'],card)
            print("Writing output")
            wDistStats.to_csv(f'{self.outDir}/wDist_card{card}.csv',sep = "\t",index=False)
            cosStats.to_csv(f'{self.outDir}/cosineSim_card{card}.csv',sep = "\t",index=False)
        for card in [3]:
            wDistStats, cosStats = self.statsForAllCardSubsets(card,probHash)
            cosCutoff = pd.Series(cosStats['3Sub2']).describe()[f'{self.qt}%']
            status = [1 if x else 0 for x in (summaryCos['mean'] <= cosCutoff)]
            cosStats['Status'] = status
            wdistCutoff = pd.Series(wDistStats['3Sub2']).describe()[f'{self.qt}%']
            wStatus = [1 if x else 0 for x in (summaryWDist['mean'] <= wdistCutoff)]
            wDistStats['Status'] = wStatus
            self.plotSimilarityHist(wDistStats['3Sub2'],cosStats['3Sub2'],card)
            print("Writing output")
            wDistStats.to_csv(f'{self.outDir}/wDist_card{card}.csv',sep = "\t",index=False)
            cosStats.to_csv(f'{self.outDir}/cosineSim_card{card}.csv',sep = "\t",index=False)
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
        cardToChoose = min(self.toChoose,len(ixList))
        print("Calculating stats for",cardToChoose,"reads")
        random.seed(self.seed)
        revised_ixes = random.sample(ixList,cardToChoose)
        C1 = []
        C2 = []
        C3 = []
        C4 = []
        for n in range(2,card):
            stats = self.getStatsPerCard(card,n,probHash,revised_ixes)
            C1.append(stats[0])
            C2.append(stats[1])
            C3.append(stats[0])
            C4.append(stats[1])
        cN = [str(card)+"Sub"+str(i) for i in range(2,card)]
        df1 = pd.DataFrame(C1).T
        df2 = pd.DataFrame(C2).T
        df3 = pd.DataFrame(C3).T
        df4 = pd.DataFrame(C4).T
        df1.columns = cN
        df1['Edge_ix'] = revised_ixes
        df2.columns = cN
        df2['Edge_ix'] = revised_ixes
        df3.columns = cN
        df3['Edge_ix'] = revised_ixes
        df4.columns = cN
        df4['Edge_ix'] = revised_ixes
        return(df1,df2,df3,df4)

    def getStatsPerCard(self,card,n,probHash,revised_ixes):
        """Per cardinality, get a subset of reads (for computational
        purposes), calculate the similarity of observed n-way interactions
        to the expected value, and output a list"""
        wdistList = []
        cosList = []
        edistList = []
        empDistList = []
        expHash = f'{card}sub{n}'
        for ix in revised_ixes:
            #print(ix)
            wDist, cos, eDist, empDist = self.getReadExpectednessStats(card,ix,n,probHash[expHash])
            wdistList.append(wDist)
            cosList.append(cos)
            edistList.append(eDist)
            empDistList.append(empDist)
        return(wdistList, cosList, edistList, empDistList)
    
    def getReadExpectednessStats(self,card,ix,n,hash):
        """Per read, get the observed distribution of n-way contacts
        and calculate similarity to the expected distribution"""
        readDict = self.makeNWayDict(ix,n)
        readPercs = self.getProbabilitiesByDist(readDict)
        probVals = [hash[k] for k in readDict.keys()]
        normProbs = [p / sum(probVals) for p in probVals]
        if self.toPlotInd is True:
            if len(readPercs) > 2 and len(probVals) > 2:
                self.makeSanityCheckPlotsPerRead(readDict,readPercs,probVals,card,n,ix)
        wDist = wasserstein_distance(readPercs, probVals)
        similarity = cosine_similarity(np.array(readPercs).reshape(1,-1), 
                                    np.array(probVals).reshape(1,-1))[0,0]
        eDist = energy_distance(readPercs, probVals)
        empDist = self.calcEmpDist(readPercs, probVals)
        return(wDist, similarity, eDist, empDist)
    
    def calcEmpDist(self,readPercs, probVals):
        """Trying a different metric for distance. 
        We will figure out the cutoff later"""
        obs = np.log(readPercs)
        exp = np.log(probVals)
        obsOverExp = obs/exp
        meanOE = obsOverExp.mean()
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
        axs[1, 1].set_title(f"Barplot for probVals Card{card}sub{n}")

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

    def plotSimilarityHist(self,wdistList,cosList,card):
        """Given distribution for a cardinality, plot the wass dist
        and cosine similarity values so that we can settle on a heuristic"""
        # Create a 2x2 grid of subplots
        fig, axs = plt.subplots(1, 2, figsize=(8, 4))

        # Plot 1: Wasserstein distance histogram
        axs[0].hist(wdistList, color='blue', alpha=0.7, width=0.3, bins=201)
        axs[0].set_xlim(-1, 1.01)
        axs[0].set_title(f"Mean wasserstein distance for card={card}")

        # Plot 2: Cosine similarity histogram
        axs[1].hist(cosList, color='green', alpha=0.7, width=0.3, bins=101)
        axs[1].set_xlim(0, 1.01)
        axs[1].set_title(f"Mean cosine similarity for card={card}")
        # Adjust layout
        plt.tight_layout()
        plt.savefig(f'{self.plotDir}/Histogram_MeanForCard{card}.png',
                    bbox_inches = 'tight',facecolor = "white")
        #plt.show()
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
        plt.ylim(-1,1)
        plt.title(f'Distribution for card={card}')
        plt.legend()
        plt.savefig(f'{self.plotDir}/ScatterPlot_{metric}_Card{card}_max{self.toChoose}reads.png',
                    bbox_inches = 'tight',facecolor = "white")
        #plt.show()
        plt.close()
        return None
