# Multiway contacts influence splicing

## Preliminary analysis for the pore-c project

Order of operations for synthetic chains
- Evaluate parameters
- Make 'reads' based off distance matrix
- Make projection matrix + incidence DFs for each of the 10k chains
- Make a hypergraph dict and combine
- Evaluate interesting reads

Evaluate if interesting reads are truly so
- Arbitrarily choose reads classified as int versus control
- Generate ensembles of chains containing these reads at various percentages (10,20,50 etc)
- Re-run pipeline


Order of operations for real PoreC data
- Process input csv and convert into hypergraph
- Evaluate interesting reads
- Compare across samples
- Isolate differential reads and generate clusters
- Visualize the clusters
