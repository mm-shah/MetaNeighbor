## ----eval = TRUE---------------------------------------------------------
library("MetaNeighbor")
data(mn_data)
data(gene_set)

## ----eval=TRUE,fig.width=4,fig.height=3----------------------------------
AUROC_scores = MetaNeighbor(data = mn_data,
                            experiment_labels = as.numeric(factor(mn_data$study_id)),
                            celltype_labels = mn_data@colData@metadata$cell_labels,
                            genesets = gene_set,
                            bplot = TRUE)

## ----eval= TRUE----------------------------------------------------------
head(AUROC_scores)

## ----eval = TRUE---------------------------------------------------------
library(MetaNeighbor)
data(mn_data)

## ----eval = TRUE---------------------------------------------------------
var_genes = variableGenes(data = mn_data, exp_labels = mn_data$study_id)
head(var_genes)

## ----eval = TRUE---------------------------------------------------------
length(var_genes)

## ----eval=TRUE-----------------------------------------------------------
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             data = mn_data, 
                             study_id = mn_data$study_id,
                             cell_type = mn_data$cell_type)

## ----eval=TRUE,fig.width=7,fig.height=6.5--------------------------------
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
gplots::heatmap.2(celltype_NV,
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  cexRow = 0.6,
                  cexCol = 0.6)

## ----eval = TRUE---------------------------------------------------------
top_hits = topHits(cell_NV = celltype_NV,
                   data = mn_data,
                   study_id = mn_data$study_id,
                   cell_type = mn_data$cell_type,
                   threshold = 0.9)
top_hits

