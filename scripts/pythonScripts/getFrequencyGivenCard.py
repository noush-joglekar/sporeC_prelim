
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
    

