def getCardProbs(incDF):
    cardVec = incDF.sum()
    cardFreq = Counter(cardVec)
    total = sum(cardFreq.values())
    probCard = {key: (value/total) for key, value in cardFreq.items()}
    return(probCard)

def getLenDistrPerCard(incDF,card):
    readIx = incDF.columns[incDF.sum() == card]
    concatemerLen = []
    for ix in readIx:
        binList = list(incDF[ix])
        ixFirst = binList.index(1)
        ixLast = (len(binList) - 1) - binList[::-1].index(1)    ## Index base 0
        cL = ixLast - ixFirst + 1
        concatemerLen.append(cL)   ## Difference w.r.t. base 0 gives 1 less so add
    return(concatemerLen)

def getFreqPerCard(lenDist):
    lenCounts=Counter(lenDist)
    #plotHistOfLenGivenCard(lenDist,card)
    total = sum(lenCounts.values())
    freqLen = {key: (value/total) for key, value in lenCounts.items()}
    return(freqLen)

def plotHistOfLenGivenCard(lenDist,card):
    fig, ax = plt.subplots(figsize =(5, 4))
    ax.hist(lenDist,density = True,bins = binSpecs)
    plt.title(f"Length | Card for c= {card}")
    plt.show()
    return()

