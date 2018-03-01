MetaNeighbor: a method to rapidly assess cell type identity using both functional and random gene sets
================
MetaNeighbor is a replication framework, that allows researchers to quantify the degree to which cell types replicate across datasets, and to rapidly identify clusters with high similarity for further testing. It calculates correlations between all pairs of cells that user aims to compare across datasets based on the expression of a set of genes. It then performs cross-dataset validation by hiding all type labels for one dataset at a time. Finally it predicts the cell type labels of the test set by running a neighbor voting algorithm to predict the identity of the held-out cells based on the similarity to the training data.

Please refer to its [documentation](./Documentation.md) and cite [Crow et al (2018) Nature Communications](https://www.nature.com/articles/s41467-018-03282-0) if you find MetaNeighbor useful in your research.
