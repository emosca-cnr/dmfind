# dmfind

Network diffusion-based analysis of "omics" data for the identification of differentially enriched modules. The tool implements functions for:

- performing network diffusion (also known as network propagation) of one or more vectors of gene-weights throughout an interactome of gene-gene interactions;
- defining ranked gene lists, where genes are prioritized according to the Network Smoothing Index (NSI), which is proportional to gene-weights and the network proximity of high gene-weights;
- assessing the characteristics of the "top" networks defined by the top ranked genes, based on proximity of high scoring genes, enrichment in high scoring genes and communities that compose the top networks;
- obtaining the functional cartography of the top networks;
- comparing multiple networks by topology;
- visualizing the top networks and the relations between communities.


If you use this package please cite:

Bersanelli, M., Mosca, E., Remondini, D. et al. Network diffusion-based analysis of high-throughput data for the detection of differentially enriched modules. Sci Rep 6, 34841 (2016). https://doi.org/10.1038/srep34841

