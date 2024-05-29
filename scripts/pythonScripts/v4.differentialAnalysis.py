import argparse
import pandas as pd
import numpy as np
import csv
import random
import sys
import os

import hypernetx as hnx
import networkx as nx
import community
import seaborn as sns

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
from matplotlib import colors

sys.path.append('/gpfs/commons/groups/gursoy_lab/ajoglekar/Projects/2023_03_01_multiwayInteractions/v0.analysis/scripts/pythonScripts/functions/')
from chains import dictToDF, dfToDict
from edgeWeightFormulations import finalBounded, finalBounded_fromEdge

def getDifference(args,f1,f2,outDir):
    npm1 = normalizeMats(args,f1,args.ct1,outDir)
    npm2 = normalizeMats(args,f2,args.ct2,outDir)
    diffMat = npm1.sub(npm2, fill_value=0)
    numeric_values = diffMat.index.to_series().str.extract(r'Bin(\d+)', 
                                                           expand=False).astype(int)
    diffMat_sorted = diffMat.loc[numeric_values.sort_values().index, 
                                 numeric_values.sort_values().index]
    print("Obtained sorted difference matrix")
    return(diffMat_sorted)

def normalizeMats(args,fPath,cT,outDir):
    print("Reading in file and normalizing matrix")
    if os.path.exists(fPath):
        pm = pd.read_pickle(fPath)
        np.fill_diagonal(np.asmatrix(pm),0)
        npm = pm / np.nanmax(pm)
        print("Plotting normalized matrix to .pdf")
        plotNormMatToPDF(args,npm,cT,outDir)
        return(npm)
    else:
        print("Main file does not exist - quitting program")
        print(fPath)
        sys.exit()

def plotNormMatToPDF(args,normMat,cT,outDir):
    outname = f'{outDir}/normProjMat_card{args.card}_{args.chrom}_{cT}.pdf'
    plt.figure(figsize=(6, 4))
    im = plt.imshow(normMat, cmap="YlOrRd",norm = LogNorm(vmax = 100))
    plt.colorbar(im, fraction=0.046, pad=0.04, label='balanced');
    plt.title(f"Projection matrix, Card={args.card}, {args.chrom}: {cT}")
    plt.savefig(outname,bbox_inches = 'tight',facecolor = "white")
    return()

def plotDiffMat(args,diffMat,outDir):
    print("Plotting and saving sorted diffMat")
    outname = f'{outDir}/diffProjMat_card{args.card}_{args.chrom}_{args.ct1}-{args.ct2}.pdf'
    plt.figure(figsize=(6, 4))
    im = plt.imshow(diffMat, cmap="BrBG", norm=colors.SymLogNorm(linthresh=0.02, linscale=0.02,
                                                vmin=-0.5, vmax=0.5, base=10))
    plt.colorbar(im, fraction=0.046, pad=0.04, label='balanced');
    plt.title(f"Diff-Proj matrix, Card={args.card}, {args.chrom}: {args.ct1} - {args.ct2}")
    plt.savefig(outname,bbox_inches = 'tight',facecolor = "white")
    return()

def getEnrichedReads(args,diffMat,d1,d2,outDir):
    print("Getting enriched reads per sample")
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
    print(enriched_pm2_incDF.shape)

    print(f"Plotting all the graphs for {args.ct1}")
    plotsPerCT(args,enriched_pm1_incDF,args.ct1,id1,outDir)
    print(f"Plotting all the graphs for {args.ct2}")
    plotsPerCT(args,enriched_pm2_incDF,args.ct2,id2,outDir)
    return([enriched_pm1_incDF,enriched_pm2_incDF])

def sort_key(item, delimiter):
    return int(item.split(delimiter)[1])

def plotRawLineGraph(args,epm,cT,incDF,outDir):
    print("Getting hypergraph")
    HIS_f1 = hnx.Hypergraph.from_incidence_dataframe(epm)
    node_names = sorted(list(list(incDF.index)), key=lambda x: sort_key(x,'Bin'))
    node_colors = sns.color_palette("flare", n_colors=len(node_names))
    color_mapping = dict(zip(node_names, node_colors))

    HIS_f1D = HIS_f1.dual()
    l_f1d = HIS_f1D.get_linegraph(s = 1)
    print("Raw line graph from hypergraph")
    outname = f'{outDir}/rawLinegraph_{cT}_readByBin_card{args.card}_{args.chrom}.pdf'
    plt.figure(figsize=(6, 4))
    a = nx.draw(l_f1d,node_size = 30, 
            with_labels = False, 
            node_color=[color_mapping[node] for node in l_f1d.nodes],
            width = 0.1,
            font_size=5)
    plt.savefig(outname,bbox_inches = 'tight',facecolor = "white")
    return(color_mapping)
    
def getNodeImportance(epm):
    readSupport = [int(rID.split(":")[1]) for rID in epm.columns]
    # finalBoundedScores = [finalBounded(list(epm[c])) for c in epm.columns]
    # node_importance = [readSupport[i]*finalBoundedScores[i] for i in range(len(readSupport))]
    node_importance = readSupport
    return(node_importance)

def addEdgeWeights(epm, node_importance):
    print("Weighted graph from incDF and read support")
    df_transposed = epm.transpose()
    # Initialize a graph
    G = nx.Graph()
    # Iterate through each row in the transposed DataFrame
    for rName, row in df_transposed.iterrows():
        # Get the nodes that are connected in this row
        connected_nodes = [column for column, value in row.items() if value == 1]
        impValue = node_importance[list(df_transposed.index).index(rName)]
        # Add edges between the connected nodes
        for i in range(len(connected_nodes)):
            for j in range(i+1, len(connected_nodes)):
                node1 = connected_nodes[i]
                node2 = connected_nodes[j]
                # Increment edge weight if edge already exists, otherwise add the edge with weight 1
                if G.has_edge(node1, node2):
                    G[node1][node2]['weight'] += 1*impValue
                else:
                    G.add_edge(node1, node2, weight=1*impValue)
    return(G)

def plotWeightedGraph(args,G,cT,outDir,color_mapping):
    print("Plotting weighted graph")
    node_names = sorted(list(list(G.nodes)), key=lambda x: sort_key(x,'Bin'))
    # node_colors = sns.color_palette("flare", n_colors=len(node_names))
    # color_mapping = dict(zip(node_names, node_colors))

    # Get edge weights
    edge_weights = [data['weight'] for u, v, data in G.edges(data=True)]
    # Normalize edge weights
    min_weight = min(edge_weights)
    max_weight = max(edge_weights)
    normalized_weights = [(weight - min_weight) / (max_weight - min_weight) for weight in edge_weights]

    outname = f'{outDir}/weightedLineGraph_{cT}_readByBin_card{args.card}_{args.chrom}.pdf'
    plt.figure(figsize=(6, 4))
    nx.draw(G,node_size = 30, 
            #with_labels = True, 
            node_color=[color_mapping[node] for node in node_names],
            width = normalized_weights,
            # pos = nx.spring_layout(G),
            font_size=5)
    plt.title('Weighted Graph of Connected Nodes')
    plt.savefig(outname,bbox_inches = 'tight',facecolor = "white")

    print("Writing weighted graph to file.")
    gName = f'{outDir}/weightedLineGraph_{cT}_readByBin_card{args.card}_{args.chrom}.graphml'
    nx.write_graphml(G, gName)
    return(normalized_weights)


def writeOutClusters(args,cT,partition,outDir):
    outName = f'{outDir}/binCommunities_{cT}_card{args.card}_{args.chrom}.csv'
    # Create communities dictionary
    communities = {}
    for node, community_id in partition.items():
        if community_id not in communities:
            communities[community_id] = [node]
        else:
            communities[community_id].append(node)

    # Write communities to a CSV file
    print("Writing out louvain clusters to csv")
    with open(outName, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Community ID', 'Nodes'])
        for community_id, nodes in communities.items():
            writer.writerow([community_id, ','.join(nodes)])
    return()

def clusterToGetCommunities(args,G,cT,normalized_weights,outDir):
    print("Clustering w/ Louvain to get communities")
    # Apply the Louvain method for community detection
    partition = community.best_partition(G)
    writeOutClusters(args,cT,partition,outDir)

    # Generate random colors for each community
    color_map = {}
    for node, community_id in partition.items():
        if community_id not in color_map:
            color_map[community_id] = (random.random(), random.random(), random.random())  # RGB color tuple

    print("Plotting louvain")
    outname = f'{outDir}/coloredByCommunity_{cT}_weightedReadByBin_card{args.card}_{args.chrom}.pdf'
    # Plot the graph with nodes colored by community
    plt.figure(figsize=(8, 6))
    pos = nx.spring_layout(G)  # Positions of nodes
    for community_id, color in color_map.items():
        nodes_in_community = [node for node, cid in partition.items() if cid == community_id]
        nx.draw_networkx_nodes(G, pos, nodelist=nodes_in_community, node_color=color, node_size=30, alpha=0.8)
    nx.draw_networkx_edges(G, pos, alpha=0.5,width = normalized_weights)
    plt.title('Graph with Communities')
    plt.savefig(outname,bbox_inches = 'tight',facecolor = "white")
    return()

def plotsPerCT(args,epm,cT,incDF,outDir):
    color_mapping = plotRawLineGraph(args,epm,cT,incDF,outDir)
    nodeImp = getNodeImportance(epm)
    weightedGraph = addEdgeWeights(epm, nodeImp)
    normalized_weights = plotWeightedGraph(args,weightedGraph,cT,outDir,color_mapping)
    clusterToGetCommunities(args,weightedGraph,cT,normalized_weights,outDir)
    return()

def getImpReadsPerCT(epm):
    print("Getting important reads in dict format")
    impReads = {}
    impReads = dfToDict(epm,impReads)
    return(impReads)

def formatReadFile_wBinID(args,bedFile,chromSizes):
    print("Formatting the BED file to assign bin IDs")
    chrFile = bedFile[bedFile['chr']==args.chrom]
    binSize = 5*10**5 #1*10**6 #HARDCODED
    chrBins = [x for x in range(0,chromSizes[args.chrom]+binSize,binSize)]
    chrFile_binned = pd.cut(chrFile['start'],bins = chrBins, labels = ["Bin"+str(i+1) for i in range(len(chrBins)-1)]).rename("binID")
    binStart = [chrBins[i] + 1 for i in chrFile_binned.cat.codes]
    binEnd = [chrBins[i+1] for i in chrFile_binned.cat.codes]
    chrFile_wBinID = chrFile.assign(binID=chrFile_binned, binStart=binStart, binEnd=binEnd)
    return(chrFile_wBinID)

def reformatFileAndSubsetImpReads(args,chrFile_wBinID,impReads,outDir,cT):
    print("Getting read IDs from BED file.... this step takes a long while....")
    groupedBins = chrFile_wBinID.groupby('ID')['binID'].apply(list).reset_index(name='Bins')
    edges = ["_".join(sorted(list(set(a)), key=lambda x: sort_key(x,'Bin'))) for a in groupedBins['Bins'] if len(list(set(a))) > 1]
    readIDs = [groupedBins.iloc[ix][0] for ix in range(len(groupedBins)) if len(list(set(groupedBins.iloc[ix][1]))) > 1]
    print("Read IDs obtained - now subsetting important reads and writing to file")
    impReads_ct = [readIDs[i] for i,x in enumerate(edges) if x in impReads]
    impRecords_ct = chrFile_wBinID[chrFile_wBinID['ID'].isin(impReads_ct)]
    outName = f'{outDir}diffContactReads_card{args.card}_{args.chrom}_{cT}.tab.gz'
    impRecords_ct.to_csv(path_or_buf=outName,index = False,sep = "\t",compression="gzip")
    return(impRecords_ct)

def finalBounded_fromEdge(edge,maxPossLen):
    """Same calculation as above except from edge IDs"""
    split_edge = edge.split("_")
    nonZeroBins = [(int(e.split(":")[1])+1)//5 if ":" in e else (int(e.split("Bin")[1])) for e in  split_edge]
    
    rCard = len(split_edge)
    ixFirst = nonZeroBins[0]
    ixLast = nonZeroBins[-1]

    concatemerLen = ixLast - ixFirst + 1
    consecBinCounts = [i - j for i,j in 
                    zip(nonZeroBins[:0:-1],nonZeroBins[-2::-1])].count(1)
    skipLen = (rCard - 1 - consecBinCounts)
    score = (skipLen + 1) * concatemerLen / (rCard * maxPossLen)
    return(score)

def readBed(args, cT):
    print("Reading in BED file, takes a while depending on size")
    readConcatemersWClosestGene = f'{args.dataDir}/NlaIII_{cT}_output_byChr/NlaIII_{cT}_{args.chrom}.gz'
    colnames = ["chr","start","end","readID","readLen","readQual",
"geneChr","geneStart","geneEnd","strand","geneID","bioType","geneName","dist","ID"]
    fullBed = pd.read_csv(readConcatemersWClosestGene,sep = "\t",names = colnames)
    return(fullBed)

def runCycle(args,f1,f2,d1,d2,chromSizes,outDir):
    print("Running the full cycle on the two datasets")
    diffMat_sorted = getDifference(args,f1,f2,outDir)
    plotDiffMat(args,diffMat_sorted,outDir)
    epm1, epm2 = getEnrichedReads(args,diffMat_sorted,d1,d2,outDir)
    for i in [0,1]:
        cT = [args.ct1,args.ct2][i]
        print(f"Running sample-specific part of pipeline for {cT}")
        epm = [epm1, epm2][i]
        bedFile = readBed(args, cT)
        chrFile_wBinID = formatReadFile_wBinID(args,bedFile,chromSizes)
        impReads = getImpReadsPerCT(epm)
        reformatFileAndSubsetImpReads(args,chrFile_wBinID,impReads,outDir,cT)
    return()

def parse_args():
    print("Parsing arguments")
    parser = argparse.ArgumentParser(
        description="""Takes in two cell types and cardinality and compares interesting reads""")
    parser.add_argument("dataDir",type=str, help="Main data dir")
    parser.add_argument('outDir',type=str, help="output directory relative to data dir")
    parser.add_argument("matDir",type=str, help="location of pkl files relative to data dir")
    parser.add_argument("chrom",type=str, help="chr ID")
    parser.add_argument("card",type=int, help="read card")
    parser.add_argument("ct1",type=str, help="cell line 1")
    parser.add_argument("ct2",type=str, help="cell line 2")
    return parser.parse_args()
 
def main():
    args = parse_args()
    #dfDir = f'{args.dataDir}{args.runDir}dfs_{args.chrom}/' ### don't need???
    outDir = f'{args.dataDir}{args.outDir}'
    f1 = f'{args.dataDir}{args.matDir}projMat_Int_{args.ct1}_Card{args.card}_{args.chrom}.pkl'
    f2 = f'{args.dataDir}{args.matDir}projMat_Int_{args.ct2}_Card{args.card}_{args.chrom}.pkl'
    d1 = f'{args.dataDir}{args.matDir}IncDF_IntReads_{args.ct1}_Card{args.card}_{args.chrom}.pkl'
    d2 = f'{args.dataDir}{args.matDir}IncDF_IntReads_{args.ct2}_Card{args.card}_{args.chrom}.pkl'
    chromSizes = pd.read_csv(f'{args.dataDir}/hg38.chromSizes',sep="\t", names = ['chr','size']).set_index('chr')['size'].to_dict()
    print(f"Starting run: {args.ct1} versus {args.ct2}")
    runCycle(args,f1,f2,d1,d2,chromSizes,outDir)
    print("All done!")
    return()


if __name__ == "__main__":
    main()
