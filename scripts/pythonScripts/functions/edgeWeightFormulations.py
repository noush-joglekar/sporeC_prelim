def unbounded_withLength(binList):
    ## Calculate unbounded weight per read
    ## mathematically, W = (skip + 1) * cLength/(card)
    rCard = binList.count(1)
    ixFirst = binList.index(1)
    ixLast = (len(binList) - 1) - binList[::-1].index(1)    ## Index base 0
    concatemerLen = ixLast - ixFirst + 1    ## Difference w.r.t. base 0 gives 1 less so add
    consecBinCounts = sum([binList[i:i+2]==[1,1] for i in range(len(binList)-1)])
    skipLen = (rCard - 1 - consecBinCounts)
    score = (skipLen + 1) * concatemerLen / rCard
    return(score)


def bounded_noLength(binList):
    ## Calculate bounded weight per read
    ## mathematically, W = (skip + 1) /((consec + 1) * card)
    rCard = binList.count(1)
    consecBinCounts = sum([binList[i:i+2]==[1,1] for i in range(len(binList)-1)])
    skip = (rCard - consecBinCounts - 1)
    score = (skip + 1) / ((consecBinCounts + 1) * rCard)
    return(score)


def bounded_lengthFreq(binList):
    ## Calculate bounded weight per read
    ## mathematically, W = (skip + 1) * (1-lengthFreq)/((consec + 1) * card)
    rCard = binList.count(1)
    ixFirst = binList.index(1)
    ixLast = (len(binList) - 1) - binList[::-1].index(1)    ## Index base 0
    concatemerLen = ixLast - ixFirst + 1 
    try:
        lGC = freqLenPerCard[rCard][concatemerLen]
    except KeyError:
        lGC = 0
    consecBinCounts = sum([binList[i:i+2]==[1,1] for i in range(len(binList)-1)])
    skip = (rCard - consecBinCounts - 1)
    score = ((skip + 1) * (1-lGC))/ ((consecBinCounts + 1) * rCard)
    return(score)

def finalBounded(binList):
    ## Calculate unbounded weight per read
    ## mathematically, W = (skip + 1) * cLength/(card * maxPossLength)
    maxPossLen = len(binList)
    rCard = binList.count(1)
    ixFirst = binList.index(1)
    ixLast = (len(binList) - 1) - binList[::-1].index(1)    ## Index base 0
    concatemerLen = ixLast - ixFirst + 1    ## Difference w.r.t. base 0 gives 1 less so add
    consecBinCounts = sum([binList[i:i+2]==[1,1] for i in range(len(binList)-1)])
    skipLen = (rCard - 1 - consecBinCounts)
    score = (skipLen + 1) * concatemerLen / (rCard * maxPossLen)
    return(score)


def finalBounded_fromEdge(edge,maxPossLen):
    """Same calculation as above except from edge IDs"""
    split_edge = edge.split("_")
    nonZeroBins = [(int(e.split(":")[1])+1)//5 for e in  split_edge]

    rCard = len(split_edge)
    ixFirst = nonZeroBins[0]
    ixLast = nonZeroBins[-1]

    concatemerLen = ixLast - ixFirst + 1
    consecBinCounts = [i - j for i,j in 
                    zip(nonZeroBins[:0:-1],nonZeroBins[-2::-1])].count(1)
    skipLen = (rCard - 1 - consecBinCounts)
    score = (skipLen + 1) * concatemerLen / (rCard * maxPossLen)
    return(score)