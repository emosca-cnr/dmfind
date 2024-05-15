# dmfind - Network diffusion-based analysis of "omics" data for the identification of differentially enriched modules.

This tool implements functions for:

- performing network diffusion (also known as network propagation) of one or more vectors of gene-weights throughout an interactome of gene-gene interactions;
- defining ranked gene lists, where genes are prioritized according to the Network Smoothing Index (NSI), which is proportional to gene-weights and the network proximity of high gene-weights;
- assessing the characteristics of the (top) networks defined by the top ranked genes, based on proximity of high scoring genes, enrichment in high scoring genes, modularity and number of communities;
- obtaining the functional cartography of the top networks;
- comparing multiple networks by topology;
- visualizing the top networks and the relations between communities.


If you use this package please cite:

Bersanelli, M., Mosca, E., Remondini, D. et al. Network diffusion-based analysis of high-throughput data for the detection of differentially enriched modules. Sci Rep 6, 34841 (2016). https://doi.org/10.1038/srep34841


## Installation 

To install this package, start R and enter:

```{r, eval=FALSE}
if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages(devtools)
}
devtools::install_github("emosca-cnr/NPATools", build_vignettes=T)
devtools::install_github("emosca-cnr/dmfind", build_vignettes=T)
```
