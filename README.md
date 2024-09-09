# Multiway contacts influence splicing

## Preliminary analysis for the pore-c project

### Order of operations for synthetic chains
Workflow:
- Evaluate parameters
- Make 'reads' based off distance matrix
- Make projection matrix + incidence DFs for each of the 10k chains
- Make a hypergraph dict and combine
- Evaluate interesting reads

Side note 1: Evaluate if interesting reads are truly so before moving to real data
- Arbitrarily choose reads classified as int versus control
- Generate ensembles of chains containing these reads at various percentages (10,20,50 etc)
- Re-run pipeline

  Results revealed this was true for a handful of tested reads

Side note 2: Second testing framework for testing algorithm validity
- Impose multiway constraints on nodes (4-5 concatemers per sample) and generate chains
- Identify which chains contain the actual concatemers
- Run read-prioritization pipeline on those chains to see if these reads pop out

  There seemed to be very high overlap and false positive rates between samples, indicating to me that the reads were truly just randomly being generated. Someone else will have to take over and re-generate the chains to make sure that the constraints are actually being met.

### Order of operations for real PoreC data
Workflow:
- Pre-process csv depending on which lab generated pore-c data
- Process input csv and convert into hypergraph
- Evaluate interesting reads
- Compare across samples
- Isolate differential reads and generate clusters
- Visualize the clusters

Mostly all done and working smoothly on (relatively) high throughput data. What's left is potentially the hardest part - connecting to good RNA-seq or TSA-seq data to get gene expression, splicing information, and nuclear positioning. Ideally should be able to build a predictive model of one of these three based on just pore-c data, but I don't know how good the signal is when the data is not quite multimodal and collected under drastically different conditions.
