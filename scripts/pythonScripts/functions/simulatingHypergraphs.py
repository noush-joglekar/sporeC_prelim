def fullyInsulatedBinIDs(ix,readCards,groundTruth_clusters):
    ## Get bin IDs assuming the matrix is fully insulated
    rcard = readCards[ix]
    clustID = groundTruth_clusters[ix]
    whichBins = [bID for bID,cID in binClust.items() if cID == clustID]
    binsSatisfyingCard = list(set(random.choices(whichBins,k = rcard)))
    return(binsSatisfyingCard)

def partiallyInsulatedBinIDs(ix,readCards,groundTruth_clusters,cutoff):
    ## If cardinality is higher than insulated clusters, then in x% of cases, simulate random pairs of bins being part of reads
    rcard = readCards[ix]
    clustID = groundTruth_clusters[ix]
    whichBins = [bID for bID,cID in binClust.items() if cID == clustID]
    binsSatisfyingCard = list(set(random.choices(whichBins,k = rcard)))
    randomProb = random.random()
    if rcard - len(binsSatisfyingCard) >= 2 and randomProb <= cutoff:
        numToSample = math.floor((rcard - len(binsSatisfyingCard)) / 2)
        addedBinIx = random.choices(range(len(binIDs)-1),k = numToSample)
        addedBinIDs = list(set(flatten([binIDs[x:x+2] for x in addedBinIx])))
        binsSatisfyingCard.extend(addedBinIDs)
    elif rcard - len(binsSatisfyingCard) >= 2 and randomProb > cutoff:
        numToSample = math.floor((rcard - len(binsSatisfyingCard)) / 2)
        addedBinIx = random.choices(range(len(whichBins)-1),k = numToSample)
        addedBinIDs = list(set(flatten([whichBins[x:x+2] for x in addedBinIx])))
        binsSatisfyingCard.extend(addedBinIDs)
    binsSatisfyingCard = list(set(binsSatisfyingCard))
    return(binsSatisfyingCard)

def realisticNestedTADs(ix,readCards,groundTruth_clusters,cutoff):
    ## If cardinality is higher than insulated clusters, then in x% of cases, simulate random pairs of bins being part of reads
    rcard = readCards[ix]
    clustID = groundTruth_clusters[ix]
    whichBins = [bID for bID,cID in binClust.items() if cID == clustID]
    binsSatisfyingCard = list(set(random.choices(whichBins,k = rcard)))
    if rcard - len(binsSatisfyingCard) >= 2:
        randomProb = random.random()
        if set(binsSatisfyingCard).issubset(secBinClust.keys()):
            sClustID = secBinClust[binsSatisfyingCard[0]]
            whichSecBins = [bID for bID,cID in secBinClust.items() if 
            cID == sClustID]
            numToSample = (rcard - len(binsSatisfyingCard) -1)
            addedBinIx = random.choices(range(len(whichSecBins)),k = numToSample)
            addedBinIDs = list(set([whichSecBins[x] for x in addedBinIx]))
            binsSatisfyingCard.extend(addedBinIDs)
        elif not set(binsSatisfyingCard).issubset(secBinClust.keys()) and randomProb <= cutoff:
            numToSample = math.floor((rcard - len(binsSatisfyingCard)) / 2)
            addedBinIx = random.choices(range(len(binIDs)-1),k = numToSample)
            addedBinIDs = list(set(flatten([binIDs[x:x+2] for x in addedBinIx])))
            binsSatisfyingCard.extend(addedBinIDs)
        binsSatisfyingCard = list(set(binsSatisfyingCard))
    return(binsSatisfyingCard)

