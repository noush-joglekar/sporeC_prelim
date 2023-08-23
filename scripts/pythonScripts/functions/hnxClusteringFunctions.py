## Main
def compareClusteringMethods(hg,K):
    HG = hmod.precompute_attributes(hg)
    rand_qh = randomPartition(hg,HG,K)
    kumarClusters = kumarAlg(HG,K)
    print(kumarClusters)
    louvainClusters = twoSectionAndLouvain(HG,K)
    print(louvainClusters)
    return([rand_qh,kumarClusters,louvainClusters])

def randomPartition(hg,HG,K):
    ## random
    V = list(hg.nodes)
    p = np.random.choice(K, size=len(V))
    RandPart = hmod.dict2part({V[i]:p[i] for i in range(len(V))})
    rand_qh = hmod.modularity(HG, RandPart)
    print(K," random clusters, qH = ",rand_qh)
    return(rand_qh)

def kumarAlg(HG,K):
    ## kumar
    kumarPart = hmod.kumar(HG)
    print("Num kumar clusters = ",len(kumarPart))
    print('Kumar qH =',hmod.modularity(HG, kumarPart))
    return(kumarPart)

def twoSectionAndLouvain(HG,K):
    ## 2-section and louvain
    G = hmod.two_section(HG)
    ## Louvain algorithm
    G.vs['louvain'] = G.community_multilevel(weights='weight').membership
    ml = hmod.dict2part({v['name']:v['louvain'] for v in G.vs})
    print("Num 2-section Louvain clusters = ",len(ml))
    print('Louvain qH =',hmod.modularity(HG, ml))
    return(ml)


def plotAlgoComparison(algo1, algo2, name1, name2):
    # Convert the clusterings to sets
    algo1_sets = [set(bin_set) for bin_set in algo1]
    algo2_sets = [set(bin_set) for bin_set in algo2]

    # Create a contingency matrix
    num_clusters_algo1 = len(algo1_sets)
    num_clusters_algo2 = len(algo2_sets)
    contingency_matrix = np.zeros((num_clusters_algo1, num_clusters_algo2), dtype=int)

    for i, cluster1 in enumerate(algo1_sets):
        for j, cluster2 in enumerate(algo2_sets):
            intersection = len(cluster1 & cluster2)
            contingency_matrix[i, j] = intersection

    # Create a heatmap
    plt.figure(figsize=(6, 4))
    sns.heatmap(contingency_matrix, annot=True, fmt="d", cmap="YlGnBu",
                xticklabels=[f"C {i}" for i in range(1, num_clusters_algo2 + 1)],
                yticklabels=[f"C {i}" for i in range(1, num_clusters_algo1 + 1)])
    plt.title("Contingency Matrix Heatmap")
    plt.xlabel(f"{name1} Clusters")
    plt.ylabel(f"{name2} Clusters")
    plt.show()

    return
